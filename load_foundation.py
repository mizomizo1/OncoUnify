#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Loader: ingest FoundationOne / FoundationOneLiquid XML reports into a SQLite DB.

Usage:
    python load_foundation.py panels.db /path/to/FoundationDir

- Even if a directory mixes FoundationOne and FoundationOneLiquid reports,
  the panel_name is auto-detected from the test-type attribute in each XML.
- The DB schema is initialised by executing schema.sql located next to this script.
"""

import os
import sys
import sqlite3
import xml.etree.ElementTree as ET
from pathlib import Path


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


def parse_pos(pos_str):
    """Split a string such as 'chr19:41727824' into (chrom, pos)."""
    if not pos_str:
        return None, None
    if ":" in pos_str:
        chrom, pos = pos_str.split(":", 1)
        return chrom, int_or_none(pos)
    else:
        return None, int_or_none(pos_str)


def infer_panel_name_from_test_type(test_type: str) -> str:
    """
    Infer the canonical panel_name from the XML test-type attribute.
    Simple rule: if the attribute contains 'liquid' (case-insensitive) return
    'FoundationOneLiquid', otherwise return 'FoundationOne'.
    """
    if not test_type:
        return "FoundationOne"
    t = test_type.lower()
    if "liquid" in t:
        return "FoundationOneLiquid"
    else:
        return "FoundationOne"


def init_db(conn: sqlite3.Connection, schema_path: Path) -> None:
    """
    Initialise the DB by executing schema.sql.
    The schema relies on CREATE TABLE IF NOT EXISTS, so this is idempotent and
    can safely be invoked multiple times.
    """
    with schema_path.open("r", encoding="utf-8") as f:
        sql = f.read()
    conn.executescript(sql)
    conn.commit()


def get_or_create_case(conn: sqlite3.Connection, case: dict) -> int:
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


def insert_variants(conn: sqlite3.Connection, case_id: int, variants: list[dict]) -> None:
    if not variants:
        return
    cur = conn.cursor()
    # Canonical column list matching the variants table defined in schema.sql
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


def insert_non_humans(conn: sqlite3.Connection, case_id: int, non_humans: list[dict]) -> None:
    if not non_humans:
        return
    cur = conn.cursor()
    columns = ["case_id", "organism", "reads_per_million", "status", "sample"]
    placeholders = ",".join(["?"] * len(columns))
    for nh in non_humans:
        values = [case_id] + [nh.get(c) for c in columns[1:]]
        cur.execute(
            f"INSERT INTO non_human_contents ({','.join(columns)}) VALUES ({placeholders})",
            values,
        )
    conn.commit()


def parse_foundation_xml(path: str | Path):
    """
    Parse a FoundationOne / FoundationOneLiquid XML report and return:
      - case dict
      - variants list
      - non-human-content list
    """

    tree = ET.parse(path)
    root = tree.getroot()

    ns_rr = {"rr": "http://integration.foundationmedicine.com/reporting"}
    ns_vr = {"vr": "http://foundationmedicine.com/compbio/variant-report-external"}

    cust = root.find("rr:CustomerInformation", ns_rr)
    ref_id = cust.findtext("rr:ReferenceID", default="", namespaces=ns_rr) if cust is not None else None
    mrn = cust.findtext("rr:MRN", default="", namespaces=ns_rr) if cust is not None else None

    payload = root.find("rr:ResultsPayload", ns_rr)
    vr = None
    if payload is not None:
        vr = payload.find("vr:variant-report", ns_vr)

    attrs = vr.attrib if vr is not None else {}

    panel_type = attrs.get("test-type", "")
    panel_name = infer_panel_name_from_test_type(panel_type)
    vendor = "FoundationMedicine"

    # The non-human-content attribute on variant-report is a case-level scalar
    non_human_attr = float_or_none(attrs.get("non-human-content"))

    case: dict = {
        "panel_name": panel_name,
        "panel_type": panel_type,
        "vendor": vendor,
        "report_id": ref_id,
        "patient_id": mrn or None,
        "sex": attrs.get("gender"),
        "age": None,
        "disease": attrs.get("disease"),
        "disease_ontology": attrs.get("disease-ontology"),
        "tissue_of_origin": attrs.get("tissue-of-origin"),
        "pathology_diagnosis": attrs.get("pathology-diagnosis"),
        "specimen_id": attrs.get("specimen"),
        "test_type": panel_type,
        "percent_tumor_nuclei": float_or_none(attrs.get("percent-tumor-nuclei")),
        "purity": float_or_none(attrs.get("purity-assessment")),
        "msi_status": None,
        "tmb_score": None,
        "tmb_status": None,
        "tmb_unit": None,
        "non_human_content": non_human_attr,
        "other_info": f"flowcell={attrs.get('flowcell-analysis','')}; pipeline={attrs.get('pipeline-version','')}",
    }

    variants: list[dict] = []
    non_humans: list[dict] = []

    if vr is not None:
        # biomarkers -> MSI/TMB
        biomarkers = vr.find("vr:biomarkers", ns_vr)
        if biomarkers is not None:
            msi = biomarkers.find("vr:microsatellite-instability", ns_vr)
            tmb = biomarkers.find("vr:tumor-mutation-burden", ns_vr)
            if msi is not None:
                case["msi_status"] = msi.attrib.get("status")
            if tmb is not None:
                case["tmb_score"] = float_or_none(tmb.attrib.get("score"))
                case["tmb_status"] = tmb.attrib.get("status")
                case["tmb_unit"] = tmb.attrib.get("unit")

        # short-variants
        svs = vr.find("vr:short-variants", ns_vr)
        if svs is not None:
            for sv in svs.findall("vr:short-variant", ns_vr):
                gene = sv.attrib.get("gene")
                pos_str = sv.attrib.get("position")
                chrom, pos = parse_pos(pos_str)
                variants.append(
                    {
                        "gene": gene,
                        "variant_type": "short_variant",
                        "variant_subtype": sv.attrib.get("functional-effect"),
                        "chrom": chrom,
                        "pos": pos,
                        "pos2": None,
                        "ref": None,
                        "alt": None,
                        "cds_effect": sv.attrib.get("cds-effect"),
                        "protein_effect": sv.attrib.get("protein-effect"),
                        "strand": sv.attrib.get("strand"),
                        "transcript": sv.attrib.get("transcript"),
                        "functional_effect": sv.attrib.get("functional-effect"),
                        "effect": None,
                        "status": sv.attrib.get("status"),
                        "origin": None,
                        "classification": None,
                        "allele_fraction": float_or_none(sv.attrib.get("allele-fraction")),
                        "depth": int_or_none(sv.attrib.get("depth")),
                        "copy_number": None,
                        "cnv_ratio": None,
                        "cnv_type": None,
                        "other_gene": None,
                        "in_frame": None,
                        "supporting_read_pairs": None,
                        "tpm": None,
                        "read_count": None,
                        "sample_name": None,
                        "raw_panel_type": "short-variant",
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

        # copy-number-alterations
        cnv_root = vr.find("vr:copy-number-alterations", ns_vr)
        if cnv_root is not None:
            for cnv in cnv_root.findall("vr:copy-number-alteration", ns_vr):
                gene = cnv.attrib.get("gene")
                pos_str = cnv.attrib.get("position")
                if pos_str and "-" in pos_str:
                    p1 = pos_str.split("-")[0]
                else:
                    p1 = pos_str
                chrom, pos = parse_pos(p1)
                variants.append(
                    {
                        "gene": gene,
                        "variant_type": "cnv",
                        "variant_subtype": cnv.attrib.get("type"),
                        "chrom": chrom,
                        "pos": pos,
                        "pos2": None,
                        "ref": None,
                        "alt": None,
                        "cds_effect": None,
                        "protein_effect": None,
                        "strand": None,
                        "transcript": None,
                        "functional_effect": None,
                        "effect": None,
                        "status": cnv.attrib.get("status"),
                        "origin": None,
                        "classification": None,
                        "allele_fraction": None,
                        "depth": None,
                        "copy_number": float_or_none(cnv.attrib.get("copy-number")),
                        "cnv_ratio": float_or_none(cnv.attrib.get("ratio")),
                        "cnv_type": cnv.attrib.get("type"),
                        "other_gene": None,
                        "in_frame": None,
                        "supporting_read_pairs": None,
                        "tpm": None,
                        "read_count": None,
                        "sample_name": None,
                        "raw_panel_type": "copy-number-alteration",
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

        # rearrangements
        rearr_root = vr.find("vr:rearrangements", ns_vr)
        if rearr_root is not None:
            for relem in rearr_root.findall("vr:rearrangement", ns_vr):
                gene = relem.attrib.get("targeted-gene") or relem.attrib.get("other-gene")

                pos1 = relem.attrib.get("pos1")
                chrom, pos = parse_pos(pos1)

                pos2_str = relem.attrib.get("pos2")
                _, pos2 = parse_pos(pos2_str) if pos2_str else (None, None)

                other_gene = relem.attrib.get("other-gene")
                in_frame   = relem.attrib.get("in-frame")
                supp_pairs = int_or_none(relem.attrib.get("supporting-read-pairs"))

                sample_name = None
                dna = relem.find("vr:dna-evidence", ns_vr)
                if dna is not None:
                    sample_name = dna.attrib.get("sample")

                variants.append(
                    {
                        "gene": gene,
                        "variant_type": "rearrangement",
                        "variant_subtype": relem.attrib.get("type"),
                        "chrom": chrom,
                        "pos": pos,
                        "pos2": pos2,
                        "ref": None,
                        "alt": None,
                        "cds_effect": None,
                        "protein_effect": None,
                        "strand": None,
                        "transcript": None,
                        "functional_effect": None,
                        "effect": relem.attrib.get("description"),
                        "status": relem.attrib.get("status"),
                        "origin": None,
                        "classification": None,
                        "allele_fraction": float_or_none(relem.attrib.get("allele-fraction")),
                        "depth": None,
                        "copy_number": None,
                        "cnv_ratio": None,
                        "cnv_type": None,
                        "other_gene": other_gene,
                        "in_frame": in_frame,
                        "supporting_read_pairs": supp_pairs,
                        "tpm": None,
                        "read_count": None,
                        "sample_name": sample_name,
                        "raw_panel_type": "rearrangement",
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

        # Per-organism rows from the non-human-content element
        nh_root = vr.find("vr:non-human-content", ns_vr)
        if nh_root is not None:
            for nh in nh_root.findall("vr:non-human", ns_vr):
                organism = nh.attrib.get("organism")
                rpm      = float_or_none(nh.attrib.get("reads-per-million"))
                status   = nh.attrib.get("status")
                sample   = None
                dna = nh.find("vr:dna-evidence", ns_vr)
                if dna is not None:
                    sample = dna.attrib.get("sample")
                non_humans.append(
                    {
                        "organism": organism,
                        "reads_per_million": rpm,
                        "status": status,
                        "sample": sample,
                    }
                )

    return case, variants, non_humans


def iter_xml_files(directory: str | Path):
    for root, dirs, files in os.walk(directory):
        for name in files:
            if name.lower().endswith(".xml"):
                yield os.path.join(root, name)


def main() -> None:
    if len(sys.argv) < 3:
        print("Usage: python load_foundation.py <db_path> <directory>", file=sys.stderr)
        sys.exit(1)

    db_path = os.path.abspath(sys.argv[1])
    directory = os.path.abspath(sys.argv[2])

    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON;")
    # Keep journal_mode=DELETE (do not enable WAL); CGI readers may otherwise

    schema_path = Path(__file__).with_name("schema.sql")
    if not schema_path.exists():
        raise SystemExit(f"schema.sql not found next to {__file__}")
    init_db(conn, schema_path)

    for path in iter_xml_files(directory):
        print(f"[INFO] Processing {path}", file=sys.stderr)
        try:
            case, variants, non_humans = parse_foundation_xml(path)
            case_id = get_or_create_case(conn, case)
            insert_variants(conn, case_id, variants)
            insert_non_humans(conn, case_id, non_humans)
        except Exception as e:
            print(f"[ERROR] Failed to parse {path}: {e}", file=sys.stderr)

    conn.close()


if __name__ == "__main__":
    main()

