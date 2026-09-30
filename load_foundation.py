#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Loader: FoundationOne CDx / FoundationOne Liquid CDx XML reports -> OncoUnify.

Usage:
    python3 load_foundation.py panels.db <file-or-directory> [...]
        [--case-metadata curation.csv] [--genome-build GRCh37|GRCh38]
        [--fail-fast] [--quiet]

* FoundationOne and FoundationOne Liquid reports may be mixed; the panel is
  recognised from variant-report/@test-type.
* Foundation Medicine XML does not state the reference assembly; coordinates
  are GRCh37 (hg19), which is recorded as a loader default
  (cases.genome_build_source = 'loader-default').  Use --genome-build to
  override it for a whole run.
* Re-running the loader replaces each report and its child rows; it never
  duplicates variants.  The process exits with status 1 if any file failed.
"""

from __future__ import annotations

import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, List, Optional

import oncounify_core as oc

LOADER = "load_foundation.py"
NS_RR = "http://integration.foundationmedicine.com/reporting"
NS_VR = "http://foundationmedicine.com/compbio/variant-report-external"
NS = {"rr": NS_RR, "vr": NS_VR}
KNOWN_FORMATS = {"ResultsReport.2.1", "variant-report-external-2.2"}
DEFAULT_BUILD = "GRCh37"


def infer_panel_name(test_type: Optional[str]) -> str:
    """'FoundationOneLiquidDx' -> 'FoundationOneLiquid'; anything else -> 'FoundationOne'."""
    return "FoundationOneLiquid" if test_type and "liquid" in test_type.lower() else "FoundationOne"


def _xsd_names(elem: Optional[ET.Element]) -> List[str]:
    if elem is None:
        return []
    loc = elem.attrib.get("{http://www.w3.org/2001/XMLSchema-instance}schemaLocation", "")
    return [re.sub(r"\.xsd$", "", Path(tok).name) for tok in loc.split() if tok.endswith(".xsd")]


def _evidence_sample(elem: ET.Element) -> Optional[str]:
    dna = elem.find("vr:dna-evidence", NS)
    return dna.attrib.get("sample") if dna is not None else None


def parse_foundation_xml(path: Path, args: Any = None) -> oc.ParsedReport:
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as exc:
        raise oc.LoaderFormatError(f"not well-formed XML: {exc}")
    if root.tag != f"{{{NS_RR}}}ResultsReport":
        raise oc.LoaderFormatError(f"root element is {root.tag!r}, expected rr:ResultsReport")
    vr = root.find("rr:ResultsPayload/vr:variant-report", NS)
    if vr is None:
        raise oc.LoaderFormatError("rr:ResultsPayload/variant-report not found")

    warnings: List[str] = []
    formats = _xsd_names(root) + _xsd_names(vr)
    unknown = [f for f in formats if f not in KNOWN_FORMATS]
    if unknown:
        warnings.append(f"untested Foundation Medicine schema version(s): {', '.join(unknown)}")

    a = vr.attrib
    cust = root.find("rr:CustomerInformation", NS)
    report_id = oc.clean_str(cust.findtext("rr:ReferenceID", namespaces=NS)) if cust is not None else None
    mrn = oc.clean_str(cust.findtext("rr:MRN", namespaces=NS)) if cust is not None else None

    test_type = oc.clean_str(a.get("test-type"))
    panel_name = infer_panel_name(test_type)
    liquid = panel_name == "FoundationOneLiquid"
    assay_label = "FoundationOne Liquid CDx" if liquid else "FoundationOne CDx"
    pipeline = oc.clean_str(a.get("pipeline-version"))

    qc = vr.find("vr:quality-control", NS)
    samples = [
        {k: s.attrib.get(k) for k in ("name", "nucleic-acid-type", "bait-set", "mean-exon-depth") if s.attrib.get(k)}
        for s in vr.findall("vr:samples/vr:sample", NS)
    ]

    case: Dict[str, Any] = {
        "panel_name": panel_name,
        "panel_version": test_type,
        "panel_type": test_type,
        "test_type": test_type,
        "vendor": "Foundation Medicine",
        "report_id": report_id,
        "patient_id": mrn,
        "sex": oc.clean_str(a.get("gender")),
        "disease": oc.clean_str(a.get("disease")),
        "disease_ontology": oc.clean_str(a.get("disease-ontology")),
        "tissue_of_origin": oc.clean_str(a.get("tissue-of-origin")),
        "pathology_diagnosis": oc.clean_str(a.get("pathology-diagnosis")),
        "specimen_id": oc.clean_str(a.get("specimen")),
        "percent_tumor_nuclei": oc.float_or_none(a.get("percent-tumor-nuclei")),
        "purity": oc.float_or_none(a.get("purity-assessment")),
        "non_human_content": oc.float_or_none(a.get("non-human-content")),
        "genome_build": DEFAULT_BUILD,
        "genome_build_source": "loader-default",
        "other_info": oc.to_json({
            "flowcell_analysis": a.get("flowcell-analysis"),
            "pipeline_version": pipeline,
            "study": a.get("study"),
            "test_request": a.get("test-request"),
            "quality_control": qc.attrib.get("status") if qc is not None else None,
            "samples": samples,
        }),
    }
    case.update(oc.provenance(path, LOADER, "; ".join(formats + ([f"pipeline {pipeline}"] if pipeline else [])) or None))

    variants: List[Dict[str, Any]] = []

    # --- short variants ------------------------------------------------------
    for sv in vr.findall("vr:short-variants/vr:short-variant", NS):
        s = sv.attrib
        chrom, pos, _ = oc.parse_locus(s.get("position"))
        variants.append(oc.short_variant(
            gene=s.get("gene"),
            protein=s.get("protein-effect"),
            cds=s.get("cds-effect"),
            vendor_term=s.get("functional-effect"),
            variant_subtype=oc.clean_str(s.get("functional-effect")),
            chrom=chrom, pos=pos,
            strand=oc.clean_str(s.get("strand")),
            transcript=oc.clean_str(s.get("transcript")),
            status=oc.clean_str(s.get("status")),
            allele_fraction=oc.float_or_none(s.get("allele-fraction")),
            depth=oc.int_or_none(s.get("depth")),
            raw_panel_type="short-variant",
            extra=oc.to_json({
                "equivocal": s.get("equivocal"),
                "subclonal": s.get("subclonal"),
                "percent_reads": s.get("percent-reads"),
                "dna_evidence": _evidence_sample(sv),
            }),
        ))

    # --- copy-number alterations --------------------------------------------
    for cnv in vr.findall("vr:copy-number-alterations/vr:copy-number-alteration", NS):
        s = cnv.attrib
        chrom, pos, end = oc.parse_locus(s.get("position"))
        variants.append({
            "gene": s.get("gene"),
            "variant_type": "cnv",
            "variant_subtype": oc.clean_str(s.get("type")),
            "chrom": chrom, "pos": pos, "pos2": end,
            "cnv_type": oc.normalize_cnv_type(s.get("type")),
            "copy_number": oc.float_or_none(s.get("copy-number")),
            "cnv_ratio": oc.float_or_none(s.get("ratio")),
            "status": oc.clean_str(s.get("status")),
            "raw_panel_type": "copy-number-alteration",
            "extra": oc.to_json({
                "equivocal": s.get("equivocal"),
                "number_of_exons": s.get("number-of-exons"),
                "dna_evidence": _evidence_sample(cnv),
            }),
        })

    # --- rearrangements -------------------------------------------------------
    for rr in vr.findall("vr:rearrangements/vr:rearrangement", NS):
        s = rr.attrib
        chrom, pos, _ = oc.parse_locus(s.get("pos1"))
        chrom2, pos2, _ = oc.parse_locus(s.get("pos2"))
        gene = oc.clean_str(s.get("targeted-gene")) or oc.clean_str(s.get("other-gene"))
        other = oc.clean_str(s.get("other-gene"))
        variants.append({
            "gene": gene,
            "variant_type": "rearrangement",
            "variant_subtype": oc.clean_str(s.get("type")),
            "chrom": chrom, "pos": pos, "chrom2": chrom2, "pos2": pos2,
            "other_gene": other if other != gene else None,
            "in_frame": oc.normalize_in_frame(s.get("in-frame")),
            "supporting_read_pairs": oc.int_or_none(s.get("supporting-read-pairs")),
            "allele_fraction": oc.float_or_none(s.get("allele-fraction")),
            "status": oc.clean_str(s.get("status")),
            "effect": oc.clean_str(s.get("description")),
            "sample_name": _evidence_sample(rr),
            "raw_panel_type": "rearrangement",
            "extra": oc.to_json({"equivocal": s.get("equivocal"), "percent_reads": s.get("percent-reads")}),
        })

    # --- biomarkers -----------------------------------------------------------
    biomarkers: List[Dict[str, Any]] = []
    msi = vr.find("vr:biomarkers/vr:microsatellite-instability", NS)
    if msi is not None and oc.clean_str(msi.attrib.get("status")):
        biomarkers.append({
            "name": "MSI",
            "call": oc.normalize_msi_call(msi.attrib.get("status")),
            "call_raw": oc.clean_str(msi.attrib.get("status")),
            "assay": f"{assay_label} MSI",
            "source_field": "biomarkers/microsatellite-instability/@status",
        })
    tmb = vr.find("vr:biomarkers/vr:tumor-mutation-burden", NS)
    if tmb is not None and (oc.clean_str(tmb.attrib.get("score")) or oc.clean_str(tmb.attrib.get("status"))):
        biomarkers.append({
            "name": "TMB",
            "value": oc.float_or_none(tmb.attrib.get("score")),
            "unit": oc.normalize_tmb_unit(tmb.attrib.get("unit")),
            "call": oc.normalize_tmb_call(tmb.attrib.get("status")),
            "call_raw": oc.clean_str(tmb.attrib.get("status")),
            "assay": f"{assay_label} {'blood TMB (bTMB)' if liquid else 'tissue TMB'}",
            "source_field": "biomarkers/tumor-mutation-burden/@score",
        })

    # --- non-human content ----------------------------------------------------
    non_humans = [{
        "organism": oc.clean_str(nh.attrib.get("organism")),
        "reads_per_million": oc.float_or_none(nh.attrib.get("reads-per-million")),
        "status": oc.clean_str(nh.attrib.get("status")),
        "sample": _evidence_sample(nh),
    } for nh in vr.findall("vr:non-human-content/vr:non-human", NS)]

    return oc.ParsedReport(case, variants, non_humans, biomarkers, warnings)


def main(argv=None) -> int:
    parser = oc.base_argument_parser(LOADER, "Load FoundationOne CDx / FoundationOne Liquid CDx XML reports.")
    return oc.run_loader(prog=LOADER, parser=parser, suffixes=(".xml",),
                         parse_file=parse_foundation_xml, argv=argv)


if __name__ == "__main__":
    sys.exit(main())
