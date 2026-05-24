#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Loader: ingest Guardant Excel reports into a SQLite DB.

Usage:
    python load_guardant.py panels.db /path/to/GuardantDir [panel_name]

- The DB schema is initialised by executing schema.sql located next to this script.
- Each Excel file is treated as one case.
- The filename is split on "_":
    the 2nd element is stored as report_id
    the 3rd element is stored as patient_id
  in the cases table.
- Expected sheet names are "SNV", "Indels", "CNAs", "Fusion", and "MSI".
  (Column-name matching is case-insensitive.)
"""

import os
import sys
import sqlite3
import re
from pathlib import Path
from typing import Dict, Any, List, Tuple

import pandas as pd


def float_or_none(x):
    if x is None:
        return None
    try:
        s = str(x).strip()
        if not s:
            return None
        return float(s)
    except Exception:
        return None


def int_or_none(x):
    if x is None:
        return None
    try:
        s = str(x).strip()
        if not s:
            return None
        return int(float(s))
    except Exception:
        return None


def init_db(conn: sqlite3.Connection, schema_path: Path) -> None:
    """
    Initialise the DB by executing schema.sql.
    Idempotent: schema uses CREATE TABLE IF NOT EXISTS.
    """
    with schema_path.open("r", encoding="utf-8") as f:
        sql = f.read()
    conn.executescript(sql)
    conn.commit()


def parse_filename_for_case(path: str | Path) -> Tuple[str, str | None]:
    """
    Split the filename on "_" and return:
      the 2nd element as report_id
      the 3rd element as patient_id
    Falls back to using the full stem as report_id when the pattern is not met.
    """
    stem = Path(path).stem
    parts = stem.split("_")
    report_id = stem
    patient_id = None
    if len(parts) >= 2:
        report_id = parts[1]
    if len(parts) >= 3:
        patient_id = parts[2]
    return report_id, patient_id


def get_or_create_case(conn: sqlite3.Connection, case: Dict[str, Any]) -> int:
    cur = conn.cursor()
    cur.execute(
        """
        SELECT case_id FROM cases
        WHERE panel_name = ? AND report_id = ?
        """,
        (case.get("panel_name"), case.get("report_id")),
    )
    row = cur.fetchone()
    if row:
        return int(row[0])

    columns = [
        "panel_name", "panel_type", "vendor", "report_id", "patient_id",
        "sex", "age", "disease", "disease_ontology", "tissue_of_origin",
        "pathology_diagnosis", "specimen_id", "test_type",
        "percent_tumor_nuclei", "purity", "msi_status",
        "tmb_score", "tmb_status", "tmb_unit",
        "non_human_content",
        "other_info",
    ]
    values = [case.get(c) for c in columns]
    placeholders = ",".join(["?"] * len(columns))
    cur.execute(
        f"INSERT INTO cases ({','.join(columns)}) VALUES ({placeholders})",
        values,
    )
    conn.commit()
    return int(cur.lastrowid)


def insert_variants(conn: sqlite3.Connection, case_id: int, variants: List[Dict[str, Any]]) -> None:
    if not variants:
        return
    cur = conn.cursor()
    # Canonical column list matching the variants table in schema.sql
    columns = [
        "case_id", "gene", "variant_type", "variant_subtype", "chrom", "pos", "pos2",
        "ref", "alt", "cds_effect", "protein_effect",
        "strand", "transcript",
        "functional_effect",
        "effect", "status", "origin", "classification", "allele_fraction",
        "depth", "copy_number", "cnv_ratio", "cnv_type",
        "other_gene", "in_frame", "supporting_read_pairs",
        "tpm",
        "read_count", "sample_name", "raw_panel_type", "extra",
        "clinvar_id", "clinvar_url", "clinvar_sig", "clinvar_match",
        "clinvar_benign", "clinvar_likely_benign", "clinvar_uncertain",
        "maf_1kg", "maf_hgvd", "maf_tommo",
        "tpm_normal_n", "tpm_normal_mean", "tpm_normal_sd",
    ]
    placeholders = ",".join(["?"] * len(columns))
    for v in variants:
        values = [case_id] + [v.get(c) for c in columns[1:]]
        cur.execute(
            f"INSERT INTO variants ({','.join(columns)}) VALUES ({placeholders})",
            values,
        )
    conn.commit()


def infer_functional_effect_from_mut_aa(mut_aa: str | None) -> str | None:
    """
    Given a protein-level variant string such as "p.G12D", return:
      - "synonymous"  when the leading and trailing amino acids are identical
      - "nonsense"    when the trailing residue is "X" or "*"
      - "missense"    otherwise
    Returns None if the input does not match the expected pattern.
    """
    if not mut_aa:
        return None
    # Examples: p.G12D, G12D, p.W26X, etc.
    m = re.search(r'([A-Za-z\*])\d+([A-Za-z\*])', mut_aa)
    if not m:
        return None
    aa_from, aa_to = m.group(1), m.group(2)
    if aa_from == aa_to:
        return "synonymous"
    if aa_to in ("X", "*"):
        return "nonsense"
    return "missense"


def split_mut_nt(mut_nt: str | None) -> Tuple[str | None, str | None]:
    """
    Split a mut_nt string on ">" and return (ref, alt).
    """
    if not mut_nt or ">" not in mut_nt:
        return None, None
    left, right = mut_nt.split(">", 1)
    left = left.strip()
    right = right.strip()
    return (left if left else None, right if right else None)


def get_sheet(df_dict: Dict[str, pd.DataFrame], name: str) -> pd.DataFrame | None:
    """
    Look up a DataFrame in the dict returned by read_excel(sheet_name=None)
    using a case-insensitive sheet-name match. Returns None if not found.
    """
    for key, df in df_dict.items():
        if key.strip().lower() == name.lower():
            return df
    return None


def parse_guardant_xlsx(path: str | Path, panel_name: str) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    """
    Parse a Guardant Excel report and return (case_dict, variants_list).
    """
    path = Path(path)
    df_dict = pd.read_excel(path, sheet_name=None)

    report_id, patient_id = parse_filename_for_case(path)

    # Read msi_score / msi_status from the MSI sheet
    tmb_score = None
    msi_status = None
    msi_df = get_sheet(df_dict, "MSI")
    if msi_df is not None and not msi_df.empty:
        df = msi_df.rename(columns=lambda c: str(c).strip().lower())
        if "msi_score" in df.columns:
            val = df["msi_score"].dropna()
            if not val.empty:
                tmb_score = float_or_none(val.iloc[0])
        if "msi_status" in df.columns:
            val = df["msi_status"].dropna()
            if not val.empty:
                msi_status = str(val.iloc[0])

    case: Dict[str, Any] = {
        "panel_name": panel_name,
        "panel_type": "Guardant",
        "vendor": "Guardant",
        "report_id": report_id,
        "patient_id": patient_id,
        "sex": None,
        "age": None,
        "disease": None,
        "disease_ontology": None,
        "tissue_of_origin": None,
        "pathology_diagnosis": None,
        "specimen_id": None,
        "test_type": "Guardant",
        "percent_tumor_nuclei": None,
        "purity": None,
        "msi_status": msi_status,
        "tmb_score": tmb_score,
        "tmb_status": None,
        "tmb_unit": None,
        "non_human_content": None,
        "other_info": f"file={path.name}",
    }

    variants: List[Dict[str, Any]] = []

    # --- SNV sheet ---
    snv_df = get_sheet(df_dict, "SNV")
    if snv_df is not None and not snv_df.empty:
        df = snv_df.rename(columns=lambda c: str(c).strip().lower())
        for _, row in df.iterrows():
            gene = row.get("gene")
            if pd.isna(gene):
                continue
            gene = str(gene).strip()
            if not gene:
                continue

            chrom = row.get("chrom")
            chrom = None if pd.isna(chrom) else str(chrom).strip()

            pos = int_or_none(row.get("position"))
            mut_nt = row.get("mut_nt")
            mut_nt = None if pd.isna(mut_nt) else str(mut_nt)

            ref, alt = split_mut_nt(mut_nt)

            mut_aa = row.get("mut_aa")
            mut_aa = None if pd.isna(mut_aa) else str(mut_aa)

            cds = row.get("cdna")
            cds = None if pd.isna(cds) else str(cds)

            transcript_id = row.get("transcript_id")
            transcript_id = None if pd.isna(transcript_id) else str(transcript_id)

            func_effect = infer_functional_effect_from_mut_aa(mut_aa)

            variants.append(
                {
                    "gene": gene,
                    "variant_type": "short_variant",
                    "variant_subtype": "SNV",
                    "chrom": chrom,
                    "pos": pos,
                    "pos2": None,
                    "ref": ref,
                    "alt": alt,
                    "cds_effect": cds,
                    "protein_effect": mut_aa,
                    "strand": None,
                    "transcript": transcript_id,
                    "functional_effect": func_effect,
                    "effect": None,
                    "status": None,
                    "origin": None,
                    "classification": None,
                    "allele_fraction": None,
                    "depth": None,
                    "copy_number": None,
                    "cnv_ratio": None,
                    "cnv_type": None,
                    "other_gene": None,
                    "in_frame": None,
                    "supporting_read_pairs": None,
                    "tpm": None,
                    "read_count": None,
                    "sample_name": None,
                    "raw_panel_type": "SNV",
                    "extra": None,
                    "clinvar_id": None,
                    "clinvar_url": None,
                    "clinvar_sig": None,
                    "clinvar_match": None,
                    "clinvar_benign": None,
                    "clinvar_likely_benign": None,
                    "clinvar_uncertain": None,
                    "maf_1kg": None,
                    "maf_hgvd": None,
                    "maf_tommo": None,
                    "tpm_normal_n": None,
                    "tpm_normal_mean": None,
                    "tpm_normal_sd": None,
                }
            )

    # --- Indels sheet ---
    indel_df = get_sheet(df_dict, "Indels")
    if indel_df is not None and not indel_df.empty:
        df = indel_df.rename(columns=lambda c: str(c).strip().lower())
        for _, row in df.iterrows():
            gene = row.get("gene")
            if pd.isna(gene):
                continue
            gene = str(gene).strip()
            if not gene:
                continue

            chrom = row.get("chrom")
            chrom = None if pd.isna(chrom) else str(chrom).strip()

            pos = int_or_none(row.get("position"))
            mut_nt = row.get("mut_nt")
            mut_nt = None if pd.isna(mut_nt) else str(mut_nt)
            ref, alt = split_mut_nt(mut_nt)

            cds = row.get("cdna")
            cds = None if pd.isna(cds) else str(cds)

            indel_type = row.get("type")
            indel_type = None if pd.isna(indel_type) else str(indel_type)

            transcript_id = row.get("transcript_id")
            transcript_id = None if pd.isna(transcript_id) else str(transcript_id)

            status = row.get("reporting_category")
            status = None if pd.isna(status) else str(status)

            variants.append(
                {
                    "gene": gene,
                    "variant_type": "short_variant",
                    "variant_subtype": "Indel",
                    "chrom": chrom,
                    "pos": pos,
                    "pos2": None,
                    "ref": ref,
                    "alt": alt,
                    "cds_effect": cds,
                    "protein_effect": None,
                    "strand": None,
                    "transcript": transcript_id,
                    "functional_effect": indel_type,
                    "effect": None,
                    "status": status,
                    "origin": None,
                    "classification": None,
                    "allele_fraction": None,
                    "depth": None,
                    "copy_number": None,
                    "cnv_ratio": None,
                    "cnv_type": None,
                    "other_gene": None,
                    "in_frame": None,
                    "supporting_read_pairs": None,
                    "tpm": None,
                    "read_count": None,
                    "sample_name": None,
                    "raw_panel_type": "Indel",
                    "extra": None,
                    "clinvar_id": None,
                    "clinvar_url": None,
                    "clinvar_sig": None,
                    "clinvar_match": None,
                    "clinvar_benign": None,
                    "clinvar_likely_benign": None,
                    "clinvar_uncertain": None,
                    "maf_1kg": None,
                    "maf_hgvd": None,
                    "maf_tommo": None,
                    "tpm_normal_n": None,
                    "tpm_normal_mean": None,
                    "tpm_normal_sd": None,
                }
            )

    # --- CNAs sheet ---
    cna_df = get_sheet(df_dict, "CNAs")
    if cna_df is not None and not cna_df.empty:
        df = cna_df.rename(columns=lambda c: str(c).strip().lower())
        for _, row in df.iterrows():
            gene = row.get("gene")
            if pd.isna(gene):
                continue
            gene = str(gene).strip()
            if not gene:
                continue

            chrom = row.get("chrom")
            chrom = None if pd.isna(chrom) else str(chrom).strip()

            copy_number = float_or_none(row.get("copy_number"))

            variants.append(
                {
                    "gene": gene,
                    "variant_type": "cnv",
                    "variant_subtype": "CNA",
                    "chrom": chrom,
                    "pos": None,
                    "pos2": None,
                    "ref": None,
                    "alt": None,
                    "cds_effect": None,
                    "protein_effect": None,
                    "strand": None,
                    "transcript": None,
                    "functional_effect": None,
                    "effect": None,
                    "status": None,
                    "origin": None,
                    "classification": None,
                    "allele_fraction": None,
                    "depth": None,
                    "copy_number": copy_number,
                    "cnv_ratio": None,
                    "cnv_type": None,
                    "other_gene": None,
                    "in_frame": None,
                    "supporting_read_pairs": None,
                    "tpm": None,
                    "read_count": None,
                    "sample_name": None,
                    "raw_panel_type": "CNA",
                    "extra": None,
                    "clinvar_id": None,
                    "clinvar_url": None,
                    "clinvar_sig": None,
                    "clinvar_match": None,
                    "clinvar_benign": None,
                    "clinvar_likely_benign": None,
                    "clinvar_uncertain": None,
                    "maf_1kg": None,
                    "maf_hgvd": None,
                    "maf_tommo": None,
                    "tpm_normal_n": None,
                    "tpm_normal_mean": None,
                    "tpm_normal_sd": None,
                }
            )

    # --- Fusion sheet ---
    fusion_df = get_sheet(df_dict, "Fusion")
    if fusion_df is not None and not fusion_df.empty:
        df = fusion_df.rename(columns=lambda c: str(c).strip().lower())
        for _, row in df.iterrows():
            # Only when the "percentage" column is non-zero
            if "percentage" not in df.columns:
                break
            pct = row.get("percentage")
            if pd.isna(pct):
                continue
            pct_val = float_or_none(pct)
            if not pct_val or pct_val == 0.0:
                continue

            gene_a = row.get("gene_a")
            gene_b = row.get("gene_b")
            if pd.isna(gene_a):
                continue
            gene_a = str(gene_a).strip()
            gene_b = None if pd.isna(gene_b) else str(gene_b).strip()

            allele_fraction = pct_val / 100.0

            variants.append(
                {
                    "gene": gene_a,
                    "variant_type": "rearrangement",
                    "variant_subtype": "fusion",
                    "chrom": None,
                    "pos": None,
                    "pos2": None,
                    "ref": None,
                    "alt": None,
                    "cds_effect": None,
                    "protein_effect": None,
                    "strand": None,
                    "transcript": None,
                    "functional_effect": None,
                    "effect": None,
                    "status": None,
                    "origin": None,
                    "classification": None,
                    "allele_fraction": allele_fraction,
                    "depth": None,
                    "copy_number": None,
                    "cnv_ratio": None,
                    "cnv_type": None,
                    "other_gene": gene_b,
                    "in_frame": None,
                    "supporting_read_pairs": None,
                    "tpm": None,
                    "read_count": None,
                    "sample_name": None,
                    "raw_panel_type": "Fusion",
                    "extra": None,
                    "clinvar_id": None,
                    "clinvar_url": None,
                    "clinvar_sig": None,
                    "clinvar_match": None,
                    "clinvar_benign": None,
                    "clinvar_likely_benign": None,
                    "clinvar_uncertain": None,
                    "maf_1kg": None,
                    "maf_hgvd": None,
                    "maf_tommo": None,
                    "tpm_normal_n": None,
                    "tpm_normal_mean": None,
                    "tpm_normal_sd": None,
                }
            )

    return case, variants


def iter_excel_files(directory: str | Path):
    for root, dirs, files in os.walk(directory):
        for name in files:
            lower = name.lower()
            if lower.endswith(".xlsx") or lower.endswith(".xls"):
                yield os.path.join(root, name)


def main() -> None:
    if len(sys.argv) < 3:
        print("Usage: python load_guardant.py <db_path> <directory> [panel_name]", file=sys.stderr)
        sys.exit(1)

    db_path = os.path.abspath(sys.argv[1])
    directory = os.path.abspath(sys.argv[2])
    if len(sys.argv) >= 4:
        panel_name = sys.argv[3]
    else:
        panel_name = "Guardant"

    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON;")
    # Keep journal_mode=DELETE (do not enable WAL)

    schema_path = Path(__file__).with_name("schema.sql")
    if not schema_path.exists():
        raise SystemExit(f"schema.sql not found next to {__file__}")
    init_db(conn, schema_path)

    for path in iter_excel_files(directory):
        print(f"[INFO] Processing {path} (panel_name={panel_name})", file=sys.stderr)
        try:
            case, variants = parse_guardant_xlsx(path, panel_name)
            case_id = get_or_create_case(conn, case)
            insert_variants(conn, case_id, variants)
        except Exception as e:
            print(f"[ERROR] Failed to parse {path}: {e}", file=sys.stderr)

    conn.close()


if __name__ == "__main__":
    main()

