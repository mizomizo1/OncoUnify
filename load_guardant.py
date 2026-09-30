#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Loader: Guardant360 CDx Excel (.xlsx) deliverables -> OncoUnify.

Usage:
    python3 load_guardant.py panels.db <file-or-directory> [...]
        [--panel-name Guardant] [--panel-version "Guardant360 CDx"]
        [--filename-pattern REGEX] [--include-uncalled]
        [--case-metadata curation.csv] [--genome-build GRCh37|GRCh38]
        [--fail-fast] [--quiet]

* Each workbook is one report.  Guardant workbooks carry no case-level
  identifiers, so report_id (and optionally patient_id) are taken from the
  file name with a configurable regular expression with named groups
  `report_id` and `patient_id`.  The default,
      ^[^_]+_(?P<patient_id>[^_]+)_(?P<report_id>[^_]+)$
  reads 'Interim_<patient_id>_<report_id>.xlsx'.  A file name that does not
  match is reported as a failure; it is never loaded with guessed identifiers.
* Sheets SNV, Indels and CNAs are required; Fusion (or Fusions) and MSI are
  optional.  Sheet and column names are matched case-insensitively.
* The workbooks list every candidate the pipeline evaluated.  Only rows with
  call = 1 (reported alterations) are loaded unless --include-uncalled is
  given, in which case the other rows are loaded with status 'not_called'.
* percentage (% cfDNA) is stored as allele_fraction = percentage / 100.
* The MSI sheet populates the `biomarkers` table (MSI score and call); it is
  never written to a TMB field.
* Coordinates are GRCh37, recorded as a loader default.
* Disease labels are not part of the deliverable; supply them through the
  case-metadata sidecar (--case-metadata).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

import oncounify_core as oc

LOADER = "load_guardant.py"
DEFAULT_PATTERN = r"^[^_]+_(?P<patient_id>[^_]+)_(?P<report_id>[^_]+)$"
DEFAULT_BUILD = "GRCh37"

REQUIRED_SHEETS = {
    "snv": ("gene", "chrom", "position", "mut_nt"),
    "indels": ("gene", "chrom", "position", "mut_nt"),
    "cnas": ("gene", "copy_number"),
}
OPTIONAL_SHEETS = {
    "fusion": ("gene_a", "gene_b"),
    "msi": ("msi_status",),
}


def split_mut_nt(mut_nt: Optional[str]) -> Tuple[Optional[str], Optional[str]]:
    """'C>T' -> ('C', 'T');  'AGT>A' -> ('AGT', 'A')."""
    s = oc.clean_str(mut_nt)
    if not s or ">" not in s:
        return None, None
    left, right = s.split(">", 1)
    return oc.clean_str(left), oc.clean_str(right)


def identifiers_from_filename(path: Path, pattern: str) -> Tuple[str, Optional[str]]:
    m = re.match(pattern, path.stem)
    if not m or not m.groupdict().get("report_id"):
        raise oc.LoaderFormatError(
            f"file name {path.name!r} does not match --filename-pattern {pattern!r}; "
            "no identifiers were guessed")
    return m.group("report_id"), m.groupdict().get("patient_id")


def _sheets(path: Path) -> Dict[str, pd.DataFrame]:
    try:
        book = pd.read_excel(path, sheet_name=None)
    except Exception as exc:
        raise oc.LoaderFormatError(f"cannot read workbook: {exc}")
    out: Dict[str, pd.DataFrame] = {}
    for name, df in book.items():
        key = str(name).strip().lower()
        if key == "fusions":
            key = "fusion"
        out[key] = df.rename(columns=lambda c: str(c).strip().lower())
    return out


def _rows(df: Optional[pd.DataFrame]):
    if df is None or df.empty:
        return []
    return [r for _, r in df.iterrows()]


def _call_filter(row, has_call: bool, include_uncalled: bool) -> Optional[str]:
    """Return the status to store, or None to skip the row."""
    if not has_call:
        return "reported"
    called = oc.int_or_none(row.get("call")) == 1
    if called:
        return "called"
    return "not_called" if include_uncalled else None


def parse_guardant_xlsx(path: Path, args: Any = None) -> oc.ParsedReport:
    pattern = getattr(args, "filename_pattern", None) or DEFAULT_PATTERN
    include_uncalled = bool(getattr(args, "include_uncalled", False))
    panel_name = getattr(args, "panel_name", None) or "Guardant"
    panel_version = getattr(args, "panel_version", None) or "Guardant360 CDx"

    report_id, patient_id = identifiers_from_filename(path, pattern)
    sheets = _sheets(path)
    warnings: List[str] = []

    for name, cols in REQUIRED_SHEETS.items():
        if name not in sheets:
            raise oc.LoaderFormatError(f"required sheet {name!r} is missing (found: {sorted(sheets)})")
        missing = [c for c in cols if c not in sheets[name].columns]
        if missing:
            raise oc.LoaderFormatError(f"sheet {name!r} lacks column(s) {missing}")
    for name, cols in OPTIONAL_SHEETS.items():
        if name not in sheets:
            warnings.append(f"optional sheet {name!r} is absent")
            continue
        missing = [c for c in cols if c not in sheets[name].columns]
        if missing:
            raise oc.LoaderFormatError(f"sheet {name!r} lacks column(s) {missing}")
    for name in ("snv", "indels", "cnas", "fusion"):
        if name in sheets and "call" not in sheets[name].columns:
            warnings.append(f"sheet {name!r} has no 'call' column; all rows were loaded as reported")

    variants: List[Dict[str, Any]] = []
    n_skipped = 0

    def status_for(row, sheet):
        nonlocal n_skipped
        st = _call_filter(row, "call" in sheets[sheet].columns, include_uncalled)
        if st is None:
            n_skipped += 1
        return st

    def pct(row) -> Optional[float]:
        v = oc.float_or_none(row.get("percentage"))
        return None if v is None else v / 100.0

    # --- SNV ----------------------------------------------------------------
    for row in _rows(sheets["snv"]):
        gene = oc.clean_str(row.get("gene"))
        if not gene:
            continue
        status = status_for(row, "snv")
        if status is None:
            continue
        ref, alt = split_mut_nt(row.get("mut_nt"))
        variants.append(oc.short_variant(
            gene=gene,
            protein=row.get("mut_aa"),
            cds=row.get("cdna"),
            hint=row.get("reporting_category"),
            variant_subtype="SNV",
            chrom=oc.normalize_chrom(row.get("chrom")), pos=oc.int_or_none(row.get("position")),
            ref=ref, alt=alt,
            transcript=oc.clean_str(row.get("transcript_id")),
            allele_fraction=pct(row),
            status=status,
            raw_panel_type="SNV",
            extra=oc.to_json({
                "exon": oc.clean_str(row.get("exon")),
                "reporting_category": oc.clean_str(row.get("reporting_category")),
                "rm_reportable": oc.clean_str(row.get("rm_reportable")),
                "percentage": oc.clean_str(row.get("percentage")),
            }),
        ))

    # --- Indels ---------------------------------------------------------------
    for row in _rows(sheets["indels"]):
        gene = oc.clean_str(row.get("gene"))
        if not gene:
            continue
        status = status_for(row, "indels")
        if status is None:
            continue
        ref, alt = split_mut_nt(row.get("mut_nt"))
        protein = oc.clean_str(row.get("mut_aa")) or oc.clean_str(row.get("mut_aa_short"))
        variants.append(oc.short_variant(
            gene=gene,
            protein=protein,
            cds=row.get("cdna"),
            hint=row.get("reporting_category"),
            variant_subtype=oc.clean_str(row.get("type")) or "Indel",
            chrom=oc.normalize_chrom(row.get("chrom")), pos=oc.int_or_none(row.get("position")),
            ref=ref, alt=alt,
            transcript=oc.clean_str(row.get("transcript_id")),
            allele_fraction=pct(row),
            status=status,
            raw_panel_type="Indels",
            extra=oc.to_json({
                "exon": oc.clean_str(row.get("exon")),
                "length": oc.clean_str(row.get("length")),
                "reporting_category": oc.clean_str(row.get("reporting_category")),
                "rm_reportable": oc.clean_str(row.get("rm_reportable")),
                "percentage": oc.clean_str(row.get("percentage")),
            }),
        ))

    # --- CNAs (Guardant360 CDx reports copy-number amplifications) -------------
    for row in _rows(sheets["cnas"]):
        gene = oc.clean_str(row.get("gene"))
        if not gene:
            continue
        status = status_for(row, "cnas")
        if status is None:
            continue
        variants.append({
            "gene": gene, "variant_type": "cnv", "variant_subtype": "CNA",
            "chrom": oc.normalize_chrom(row.get("chrom")),
            "copy_number": oc.float_or_none(row.get("copy_number")),
            "cnv_type": "amplification" if status in ("called", "reported") else None,
            "status": status,
            "raw_panel_type": "CNAs",
        })

    # --- Fusions ----------------------------------------------------------------
    fusion = sheets.get("fusion")
    for row in _rows(fusion):
        gene_a = oc.clean_str(row.get("gene_a"))
        if not gene_a:
            continue
        if "call" in fusion.columns:
            status = status_for(row, "fusion")
        else:  # older deliverables: a non-zero percentage marks a detected fusion
            status = "reported" if (oc.float_or_none(row.get("percentage")) or 0) > 0 else None
            if status is None:
                n_skipped += 1
        if status is None:
            continue
        variants.append({
            "gene": gene_a, "variant_type": "rearrangement", "variant_subtype": "fusion",
            "other_gene": oc.clean_str(row.get("gene_b")),
            "chrom": oc.normalize_chrom(row.get("chrom_a")), "pos": oc.int_or_none(row.get("pos_a")) or None,
            "chrom2": oc.normalize_chrom(row.get("chrom_b")), "pos2": oc.int_or_none(row.get("pos_b")) or None,
            "allele_fraction": pct(row),
            "status": status,
            "raw_panel_type": "Fusion",
            "extra": oc.to_json({"downstream_gene": oc.clean_str(row.get("downstream_gene"))}),
        })

    # --- MSI -> biomarkers ------------------------------------------------------
    biomarkers: List[Dict[str, Any]] = []
    msi_info: Dict[str, Any] = {}
    msi_rows = _rows(sheets.get("msi"))
    if msi_rows:
        m = msi_rows[0]
        biomarkers.append({
            "name": "MSI",
            "value": oc.float_or_none(m.get("msi_score")),
            "unit": "Guardant MSI score",
            "call": oc.normalize_msi_call(m.get("msi_status")),
            "call_raw": oc.clean_str(m.get("msi_status")),
            "assay": "Guardant360 CDx MSI (plasma cfDNA)",
            "source_field": "MSI sheet: msi_score, msi_status",
        })
        msi_info = {k: oc.clean_str(m.get(k)) for k in ("runid", "run_sample_id", "max_maf", "mean_maf")}

    qc = sheets.get("qc")
    qc_status = qc["status"].astype(str).value_counts().to_dict() if qc is not None and "status" in qc.columns else None

    case: Dict[str, Any] = {
        "panel_name": panel_name,
        "panel_version": panel_version,
        "panel_type": panel_version,
        "test_type": panel_version,
        "vendor": "Guardant Health",
        "report_id": report_id,
        "patient_id": patient_id,
        "genome_build": DEFAULT_BUILD,
        "genome_build_source": "loader-default",
        "other_info": oc.to_json({
            **msi_info,
            "qc_status_counts": qc_status,
            "uncalled_rows_skipped": n_skipped or None,
        }),
    }
    signature = ",".join(sorted(sheets))
    case.update(oc.provenance(path, LOADER, f"Guardant XLSX sheets [{signature}]"))
    return oc.ParsedReport(case, variants, [], biomarkers, warnings)


def main(argv=None) -> int:
    parser = oc.base_argument_parser(LOADER, "Load Guardant360 CDx Excel deliverables.")
    parser.add_argument("--panel-name", default="Guardant", help="panel_name to record (default: Guardant)")
    parser.add_argument("--panel-version", default="Guardant360 CDx",
                        help="panel_version to record; must match panels/*.tsv for statistics")
    parser.add_argument("--filename-pattern", default=DEFAULT_PATTERN,
                        help="regular expression applied to the file stem, with named groups "
                             "report_id (required) and patient_id (optional)")
    parser.add_argument("--include-uncalled", action="store_true",
                        help="also load rows with call = 0, stored with status 'not_called'")
    return oc.run_loader(prog=LOADER, parser=parser, suffixes=(".xlsx",),
                         parse_file=parse_guardant_xlsx, argv=argv)


if __name__ == "__main__":
    sys.exit(main())
