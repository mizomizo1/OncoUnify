#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Loader: ingest GenMineTOP (Todai OncoPanel) XML reports into a SQLite DB (full version).

Features:
- Initialise DB by executing schema.sql
- Populate cases.date with the date extracted from <accepted> (YYYY-MM-DD)
- Populate variants.transcript with <transcript> (strand is left as None)
- Strip the leading 'p.' prefix (if any) before storing variants.protein_effect
- Infer variants.functional_effect from the normalised protein_effect
  (missense / nonsense / frameshift / splice / synonymous)
- Populate ClinVar, MAF, and normal-tissue TPM (reference/normal-expression/tpm)

Usage:
    python load_genminetop.py panels.db /path/to/GenMineTOP [panel_name]
"""

import os
import re
import sys
import sqlite3
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Optional, Tuple, Dict, Any, List


# -------------------------
# helpers: number parsing
# -------------------------
def float_or_none(x) -> Optional[float]:
    if x is None:
        return None
    try:
        s = str(x).strip()
        if not s:
            return None
        return float(s)
    except Exception:
        return None


def int_or_none(x) -> Optional[int]:
    if x is None:
        return None
    try:
        s = str(x).strip()
        if not s:
            return None
        return int(float(s))
    except Exception:
        return None


# -------------------------
# helpers: variant type normalize
# -------------------------
def normalize_variant_type(raw_type: Optional[str]) -> Optional[str]:
    if not raw_type:
        return None
    t = raw_type.lower()
    if t in (
        "snv",
        "substitution",
        "insertion",
        "deletion",
        "indel",
        "frameshift",
        "splice",
        "splicing-variant",
    ):
        return "short_variant"
    if t.startswith("cnv") or "copy-number" in t or "amplification" in t:
        return "cnv"
    if t in ("fusion", "rearrangement"):
        return "rearrangement"
    if t == "expression":
        return "expression"
    if t in ("tmb", "msi"):
        return "biomarker"
    return raw_type


# -------------------------
# helpers: schema init
# -------------------------
def init_db(conn: sqlite3.Connection, schema_path: Path) -> None:
    with schema_path.open("r", encoding="utf-8") as f:
        sql = f.read()
    conn.executescript(sql)
    conn.commit()


# -------------------------
# helpers: accepted -> date
# -------------------------
def extract_accepted_date_iso(report_node: Optional[ET.Element]) -> Optional[str]:
    """
    Extract a date from <accepted>...</accepted> and normalise to YYYY-MM-DD.
    Examples:
      2025-12-11
      2025/12/11
      2025-12-11T16:20:00
    """
    if report_node is None:
        return None
    txt = report_node.findtext(".//accepted")
    if not txt:
        return None
    txt = txt.strip()
    m = re.search(r"(\d{4})[-/](\d{1,2})[-/](\d{1,2})", txt)
    if not m:
        return None
    y = int(m.group(1))
    mo = int(m.group(2))
    d = int(m.group(3))
    return f"{y:04d}-{mo:02d}-{d:02d}"


# -------------------------
# helpers: protein normalization + functional_effect inference
# -------------------------
def normalize_protein_effect(pe: Optional[str]) -> Optional[str]:
    """
    Strip leading 'p.' or 'p' prefix from protein_effect.
    Examples:
      p.G12D        -> G12D
      p.T887Rfs*19 -> T887Rfs*19
      G12D         -> G12D
    """
    if not pe:
        return None
    pe = pe.strip()
    pe = re.sub(r"^p\.?\s*", "", pe, flags=re.IGNORECASE)
    pe = pe.strip()
    return pe or None


def infer_functional_effect_from_protein(pe: Optional[str]) -> Optional[str]:
    """
    Infer functional_effect from a normalised protein_effect.
    Target categories: missense / nonsense / frameshift / splice / synonymous
    """
    if not pe:
        return None

    s = pe.strip()

    # frameshift: contains "fs" (e.g., T887Rfs*19)
    if re.search(r"fs", s, re.IGNORECASE):
        return "frameshift"

    # nonsense: ends with "*" or "X", or contains "Ter" (e.g., W26*, W26X, Trp26Ter)
    if re.search(r"(\*|X)$", s) or re.search(r"Ter", s, re.IGNORECASE):
        return "nonsense"

    # splice: any occurrence of the substring 'splice' (safety net)
    if "splice" in s.lower():
        return "splice"

    # synonymous / missense: typical A12B form (A == B means synonymous)
    m = re.match(r"^([A-Z])(\d+)([A-Z])$", s)
    if m:
        aa1, _, aa2 = m.group(1), m.group(2), m.group(3)
        if aa1 == aa2:
            return "synonymous"
        return "missense"

    # Three-letter HGVS form (e.g., "Gly12Asp"), if present
    m3 = re.match(r"^([A-Za-z]{3})(\d+)([A-Za-z]{3})$", s)
    if m3:
        a1, _, a2 = m3.group(1).lower(), m3.group(2), m3.group(3).lower()
        if a1 == a2:
            return "synonymous"
        return "missense"

    return None


# -------------------------
# DB: cases
# -------------------------
def get_or_create_case(conn: sqlite3.Connection, case: Dict[str, Any]) -> int:
    cur = conn.cursor()
    cur.execute(
        "SELECT case_id FROM cases WHERE panel_name = ? AND report_id = ?",
        (case.get("panel_name"), case.get("report_id")),
    )
    row = cur.fetchone()
    if row:
        return int(row[0])

    # Match the cases table defined in schema.sql (including date)
    columns = [
        "panel_name", "panel_type", "vendor", "report_id", "patient_id",
        "sex", "age", "disease", "disease_ontology", "tissue_of_origin",
        "pathology_diagnosis", "specimen_id", "test_type",
        "date",
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


# -------------------------
# DB: variants
# -------------------------
def insert_variants(conn: sqlite3.Connection, case_id: int, variants: List[Dict[str, Any]]) -> None:
    if not variants:
        return
    cur = conn.cursor()

    # Canonical column order matching the variants table in schema.sql
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


# -------------------------
# XML helpers
# -------------------------
def extract_gene_from_item(item: ET.Element) -> Optional[str]:
    """
    <gene>ATRX</gene> or <gene><item>PTPRZ1</item><item>MET</item></gene>
    """
    g = item.find("gene")
    if g is None:
        return None
    children = list(g)
    if children:
        names = [c.text for c in children if c.text]
        return "-".join(names) if names else None
    return g.text


# -------------------------
# XML parsing
# -------------------------
def parse_todai_xml(path: str | Path, panel_name: str):
    tree = ET.parse(path)
    root = tree.getroot()

    report = root.find("report")
    patient = report.find("patient") if report is not None else None
    specimen = report.find("specimen") if report is not None else None
    result = report.find("result") if report is not None else None
    marker = result.find("marker") if result is not None else None

    report_id = report.findtext("id") if report is not None else None
    pathology = specimen.findtext("pathology") if specimen is not None else None

    accepted_date = extract_accepted_date_iso(report)

    case: Dict[str, Any] = {
        "panel_name": panel_name,
        "panel_type": report.findtext("test/type") if report is not None else None,
        "vendor": "TodaiOncopanel",
        "report_id": report_id,
        "patient_id": patient.findtext("id") if patient is not None else None,
        "sex": patient.findtext("sex") if patient is not None else None,
        "age": int_or_none(patient.findtext("age") if patient is not None else None),
        "disease": pathology,
        "disease_ontology": None,
        "tissue_of_origin": None,
        "pathology_diagnosis": pathology,
        "specimen_id": specimen.findtext("id") if specimen is not None else None,
        "test_type": report.findtext("test/type") if report is not None else None,
        "date": accepted_date,
        "percent_tumor_nuclei": float_or_none(
            report.findtext("qc/tumor-content/nuclei/value") if report is not None else None
        ),
        "purity": float_or_none(
            report.findtext("qc/tumor-content/estimated/value") if report is not None else None
        ),
        "msi_status": None,
        "tmb_score": float_or_none(
            marker.findtext("tmb/exon/frequency-non-synonymous-alterations") if marker is not None else None
        ),
        "tmb_status": None,
        "tmb_unit": "mutations-per-megabase" if marker is not None else None,
        "non_human_content": None,
        "other_info": None,
    }

    variants: List[Dict[str, Any]] = []
    if result is not None:
        alts_root = result.find("alterations")
        if alts_root is not None:
            for it in alts_root.findall("item"):
                gene = extract_gene_from_item(it)
                raw_type = it.findtext("type")
                vtype = normalize_variant_type(raw_type or "")

                locus = it.findtext("locus")
                if locus and ":" in locus:
                    chrom, pos = locus.split(":", 1)
                else:
                    chrom, pos = None, None

                transcript = it.findtext("transcript")

                af_str = it.findtext("allele-frequency")
                allele_fraction = None
                if af_str and "/" in af_str:
                    num, den = af_str.split("/", 1)
                    n = float_or_none(num)
                    d = float_or_none(den)
                    if d:
                        allele_fraction = n / d

                # --- ClinVar URL (dbs/clinvar/item) ---
                clinvar_url = None
                dbs_node = it.find("dbs")
                if dbs_node is not None:
                    cv = dbs_node.find("clinvar")
                    if cv is not None:
                        clinvar_url = cv.findtext("item")

                # --- clinical-relevance ---
                clinvar_id = clinvar_sig = clinvar_match = None
                clinvar_benign = clinvar_likely_benign = clinvar_uncertain = None
                maf_1kg = maf_hgvd = maf_tommo = None

                cr = it.find("clinical-relevance")
                if cr is not None:
                    cv2 = cr.find("dbs/clinvar")
                    if cv2 is not None:
                        clinvar_id = cv2.findtext("id/item")
                        clinvar_sig = cv2.findtext("clinical-significance/item")
                        clinvar_match = cv2.findtext("match/item")

                    content = cr.find("content")
                    if content is not None:
                        clinvar_benign = int_or_none(content.findtext("benign"))
                        clinvar_likely_benign = int_or_none(content.findtext("likely-benign"))
                        clinvar_uncertain = int_or_none(content.findtext("uncertain-significance"))

                    maf = cr.find("minor-allele-frequency")
                    if maf is not None:
                        for key in maf.findall("key"):
                            if key.attrib.get("name") == "1000genomes":
                                maf_1kg = float_or_none(key.text)
                        maf_hgvd = float_or_none(maf.findtext("hgvd"))
                        maf_tommo = float_or_none(maf.findtext("tommo-8p3kjpn"))

                # --- normal tpm ---
                tpm_normal_n = tpm_normal_mean = tpm_normal_sd = None
                ref_node = it.find("reference")
                if ref_node is not None:
                    tpm_node = ref_node.find("normal-expression/tpm")
                    if tpm_node is not None:
                        tpm_normal_n = int_or_none(tpm_node.findtext("n"))
                        tpm_normal_mean = float_or_none(tpm_node.findtext("mean"))
                        tpm_normal_sd = float_or_none(tpm_node.findtext("sd"))

                # case tpm
                tpm_case = float_or_none(it.findtext("tpm"))

                # Normalise protein_effect and infer functional_effect
                raw_pe = it.findtext("protein-alteration")
                norm_pe = normalize_protein_effect(raw_pe)
                func_eff = infer_functional_effect_from_protein(norm_pe)

                variants.append(
                    {
                        "gene": gene,
                        "variant_type": vtype,
                        "variant_subtype": raw_type,
                        "chrom": chrom,
                        "pos": int_or_none(pos),
                        "pos2": None,
                        "ref": it.findtext("ref"),
                        "alt": it.findtext("alt"),
                        "cds_effect": it.findtext("coding-dna-alteration"),
                        "protein_effect": norm_pe,           # "p." prefix removed
                        "strand": None,
                        "transcript": transcript,
                        "functional_effect": func_eff,       # inferred
                        "effect": None,
                        "status": it.findtext("status"),
                        "origin": it.findtext("origin"),
                        "classification": it.findtext("ag-class"),
                        "allele_fraction": allele_fraction,
                        "depth": None,
                        "copy_number": float_or_none(it.findtext("num-copy")),
                        "cnv_ratio": float_or_none(it.findtext("ratio")),
                        "cnv_type": raw_type if raw_type and "cnv" in raw_type.lower() else None,
                        "other_gene": None,
                        "in_frame": None,
                        "supporting_read_pairs": None,
                        "tpm": tpm_case,
                        "read_count": int_or_none(it.findtext("num-reads")),
                        "sample_name": None,
                        "raw_panel_type": raw_type,
                        "extra": None,
                        "clinvar_id": clinvar_id,
                        "clinvar_url": clinvar_url,
                        "clinvar_sig": clinvar_sig,
                        "clinvar_match": clinvar_match,
                        "clinvar_benign": clinvar_benign,
                        "clinvar_likely_benign": clinvar_likely_benign,
                        "clinvar_uncertain": clinvar_uncertain,
                        "maf_1kg": maf_1kg,
                        "maf_hgvd": maf_hgvd,
                        "maf_tommo": maf_tommo,
                        "tpm_normal_n": tpm_normal_n,
                        "tpm_normal_mean": tpm_normal_mean,
                        "tpm_normal_sd": tpm_normal_sd,
                    }
                )

    return case, variants


def iter_xml_files(directory: str | Path):
    for root, _, files in os.walk(directory):
        for name in files:
            if name.lower().endswith(".xml"):
                yield os.path.join(root, name)


def main() -> None:
    if len(sys.argv) < 3:
        print("Usage: python load_genminetop.py <db_path> <directory> [panel_name]", file=sys.stderr)
        sys.exit(1)

    db_path = os.path.abspath(sys.argv[1])
    directory = os.path.abspath(sys.argv[2])
    panel_name = sys.argv[3] if len(sys.argv) >= 4 else "GenMineTOP"

    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON;")
    # Keep journal_mode=DELETE (do not enable WAL)

    schema_path = Path(__file__).with_name("schema.sql")
    if not schema_path.exists():
        raise SystemExit(f"schema.sql not found next to {__file__}")
    init_db(conn, schema_path)

    for path in iter_xml_files(directory):
        print(f"[INFO] Processing {path} (panel_name={panel_name})", file=sys.stderr)
        try:
            case, variants = parse_todai_xml(path, panel_name)
            case_id = get_or_create_case(conn, case)
            insert_variants(conn, case_id, variants)
        except Exception as e:
            print(f"[ERROR] Failed to parse {path}: {e}", file=sys.stderr)

    conn.close()


if __name__ == "__main__":
    main()

