#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Loader: GenMineTOP Cancer Genome Profiling System XML reports -> OncoUnify.

GenMineTOP (GenMine Labs, Inc., formerly Konica Minolta REALM, Inc.;
co-developed with the University of Tokyo and the National Cancer Center
Research Institute) is the productized successor of the Todai OncoPanel.  Its
XML declares vendor-id 'todai-oncopanel'.

Usage:
    python3 load_genminetop.py panels.db <file-or-directory> [...]
        [--panel-name GenMineTOP] [--case-metadata curation.csv]
        [--genome-build GRCh37|GRCh38] [--fail-fast] [--quiet]

* The reference assembly is read from report/result/reference-genome
  (fallback: input/options/genome-version).
* Tumor-normal paired: each alteration carries origin somatic/germline.
* Fusions and RNA exon-skipping events become `rearrangement` rows with both
  breakpoints in pos/pos2 (chrom/chrom2); exon skipping is classified as the
  Sequence Ontology term exon_loss_variant.
* specimen/pathology is stored as pathology_diagnosis only; disease and
  tissue_of_origin are left for the case-metadata sidecar.
"""

from __future__ import annotations

import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import oncounify_core as oc

LOADER = "load_genminetop.py"
KNOWN_VERSIONS = {"2"}

SHORT_TYPES = {"snv", "substitution", "mnv", "insertion", "deletion", "indel", "delins", "duplication"}
CNV_PREFIXES = ("cnv", "copy-number")


def _text(e: Optional[ET.Element], path: str) -> Optional[str]:
    return oc.clean_str(e.findtext(path)) if e is not None else None


def _items(e: Optional[ET.Element], tag: str) -> List[str]:
    """<tag>X</tag> -> [X];  <tag><item>X</item><item>Y</item></tag> -> [X, Y]."""
    if e is None:
        return []
    node = e.find(tag)
    if node is None:
        return []
    children = list(node)
    if children:
        return [oc.clean_str(c.text) for c in children if oc.clean_str(c.text)]
    t = oc.clean_str(node.text)
    return [t] if t else []


def accepted_date(report: ET.Element) -> Optional[str]:
    txt = _text(report, "date/accepted")
    if not txt:
        return None
    m = re.search(r"(\d{4})[-/](\d{1,2})[-/](\d{1,2})", txt)
    return f"{int(m.group(1)):04d}-{int(m.group(2)):02d}-{int(m.group(3)):02d}" if m else None


def allele_fraction(raw: Optional[str]) -> Tuple[Optional[float], Optional[int]]:
    """'123/456' (alt reads / depth) -> (0.2697, 456)."""
    if not raw or "/" not in raw:
        return oc.float_or_none(raw), None
    num, den = raw.split("/", 1)
    n, d = oc.float_or_none(num), oc.float_or_none(den)
    if not d:
        return None, None
    return (n / d if n is not None else None), int(d)


def clinical_relevance(it: ET.Element) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    out["clinvar_url"] = _text(it, "dbs/clinvar/item")
    cr = it.find("clinical-relevance")
    if cr is None:
        return out
    cv = cr.find("dbs/clinvar")
    if cv is not None:
        out["clinvar_id"] = _text(cv, "id/item")
        out["clinvar_sig"] = _text(cv, "clinical-significance/item")
        out["clinvar_match"] = _text(cv, "match/item")
    content = cr.find("content")
    if content is not None:
        out["clinvar_benign"] = oc.int_or_none(content.findtext("benign"))
        out["clinvar_likely_benign"] = oc.int_or_none(content.findtext("likely-benign"))
        out["clinvar_uncertain"] = oc.int_or_none(content.findtext("uncertain-significance"))
    maf = cr.find("minor-allele-frequency")
    if maf is not None:
        for key in maf.findall("key"):
            if key.attrib.get("name") == "1000genomes":
                out["maf_1kg"] = oc.float_or_none(key.text)
        out["maf_hgvd"] = oc.float_or_none(maf.findtext("hgvd"))
        out["maf_tommo"] = oc.float_or_none(maf.findtext("tommo-8p3kjpn"))
    return out


def breakpoints(it: ET.Element) -> List[Dict[str, Optional[str]]]:
    return [{
        "gene": _text(b, "gene"), "transcript": _text(b, "transcript"),
        "chr": _text(b, "chr"), "pos": _text(b, "pos"),
        "region": _text(b, "region"), "index": _text(b, "index"),
    } for b in it.findall("breakpoint/item")]


def parse_genminetop_xml(path: Path, args: Any = None) -> oc.ParsedReport:
    panel_name = getattr(args, "panel_name", None) or "GenMineTOP"
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as exc:
        raise oc.LoaderFormatError(f"not well-formed XML: {exc}")
    report = root.find("report")
    if root.tag != "root" or report is None:
        raise oc.LoaderFormatError("expected <root><report>...; this is not a GenMineTOP report")

    warnings: List[str] = []
    version = _text(root, "version")
    vendor_id = _text(root, "vendor-id")
    if version not in KNOWN_VERSIONS:
        warnings.append(f"untested GenMineTOP XML version {version!r}")
    if vendor_id and vendor_id != "todai-oncopanel":
        warnings.append(f"unexpected vendor-id {vendor_id!r}")

    patient = report.find("patient")
    specimen = report.find("specimen")
    result = report.find("result")
    marker = result.find("marker") if result is not None else None

    build_raw = _text(result, "reference-genome") or _text(root, "input/options/genome-version")
    build, build_source = None, None
    if build_raw:
        try:
            build, build_source = oc.normalize_genome_build(build_raw), "vendor-file"
        except ValueError as exc:
            raise oc.LoaderFormatError(str(exc))

    panel_items = [oc.clean_str(i.text) for i in report.findall("test/panel/item") if oc.clean_str(i.text)]
    test_type = _text(report, "test/type")
    tmb_count = _text(marker, "tmb/exon/num-non-synonymous-alterations")

    case: Dict[str, Any] = {
        "panel_name": panel_name,
        "panel_version": "+".join(panel_items) or None,
        "panel_type": test_type,
        "test_type": test_type,
        "vendor": "GenMine Labs",
        "report_id": _text(report, "id"),
        "patient_id": _text(patient, "id"),
        "sex": _text(patient, "sex"),
        "age": oc.int_or_none(_text(patient, "age")),
        "date": accepted_date(report),
        "genome_build": build,
        "genome_build_source": build_source,
        "pathology_diagnosis": _text(specimen, "pathology"),
        "specimen_id": _text(specimen, "id"),
        "percent_tumor_nuclei": oc.float_or_none(_text(report, "qc/tumor-content/nuclei/value")),
        "purity": oc.float_or_none(_text(report, "qc/tumor-content/estimated/value")),
        "other_info": oc.to_json({
            "c_cat_id": _text(patient, "c-cat-id"),
            "germline_disclosure": _text(report, "preference/germline-disclosure"),
            "tmb_num_non_synonymous_alterations": tmb_count,
            "annotation_dbs": {i.findtext("name"): i.findtext("version") for i in report.findall("reference/db/item")},
            "population_dbs": {i.findtext("name"): i.findtext("version") for i in report.findall("reference/snp-db/item")},
        }),
    }
    case.update(oc.provenance(path, LOADER, f"{vendor_id or 'unknown'} XML v{version}"))

    biomarkers: List[Dict[str, Any]] = []
    tmb_val = oc.float_or_none(_text(marker, "tmb/exon/frequency-non-synonymous-alterations"))
    if tmb_val is not None:
        biomarkers.append({
            "name": "TMB",
            "value": tmb_val,
            "unit": "mutations/Mb",
            "assay": "GenMineTOP exonic non-synonymous alteration frequency (tumor-normal paired)",
            "source_field": "result/marker/tmb/exon/frequency-non-synonymous-alterations",
        })

    variants: List[Dict[str, Any]] = []
    for it in (result.findall("alterations/item") if result is not None else []):
        raw_type = _text(it, "type") or ""
        t = raw_type.lower()
        genes = _items(it, "gene")
        gene = genes[0] if genes else None
        transcripts = _items(it, "transcript")
        loci = _items(it, "locus")
        chrom, pos, end = oc.parse_locus(loci[0]) if loci else (None, None, None)
        common = {
            "variant_subtype": raw_type or None,
            "transcript": transcripts[0] if transcripts else None,
            "status": _text(it, "status"),
            "origin": oc.normalize_origin(_text(it, "origin")),
            "classification": _text(it, "ag-class"),
            "raw_panel_type": raw_type or None,
        }
        cr = clinical_relevance(it)
        extra: Dict[str, Any] = {
            "cytoband": _items(it, "cytoband"),
            "vendor_alteration_id": _text(it, "id/vendor"),
            "cosmic": _items(it.find("dbs"), "cosmic"),
        }

        if t in SHORT_TYPES:
            af, depth = allele_fraction(_text(it, "allele-frequency"))
            extra["allele_frequency_raw"] = _text(it, "allele-frequency")
            extra["allele_frequency_tumor"] = _text(it, "allele-frequency-tumor")
            variants.append(oc.short_variant(
                gene=gene,
                protein=_text(it, "protein-alteration"),
                cds=_text(it, "coding-dna-alteration"),
                hint=raw_type,
                chrom=chrom, pos=pos,
                ref=_text(it, "ref"), alt=_text(it, "alt"),
                allele_fraction=af, depth=depth,
                extra=oc.to_json(extra),
                **common, **cr,
            ))
        elif t.startswith(CNV_PREFIXES):
            variants.append({
                "gene": gene, "variant_type": "cnv",
                "chrom": chrom, "pos": pos, "pos2": end,
                "cnv_type": oc.normalize_cnv_type(raw_type),
                "copy_number": oc.float_or_none(_text(it, "num-copy")),
                "cnv_ratio": oc.float_or_none(_text(it, "ratio")),
                "extra": oc.to_json(extra),
                **common, **cr,
            })
        elif t in ("fusion", "rearrangement", "splicing-variant"):
            bps = breakpoints(it)
            c1, p1 = (oc.normalize_chrom(bps[0]["chr"]), oc.int_or_none(bps[0]["pos"])) if bps else (chrom, pos)
            if len(bps) > 1:
                c2, p2 = oc.normalize_chrom(bps[1]["chr"]), oc.int_or_none(bps[1]["pos"])
            elif len(loci) > 1:
                c2, p2, _ = oc.parse_locus(loci[1])
            else:
                c2, p2 = None, None
            extra.update({
                "breakpoints": bps,
                "num_wt_reads": _items(it, "num-wt-reads"),
                "sv_type": _text(it, "sv-type"),
                "transcripts": transcripts,
            })
            row = {
                "gene": gene, "variant_type": "rearrangement",
                "chrom": c1, "pos": p1, "chrom2": c2, "pos2": p2,
                "other_gene": genes[1] if len(genes) > 1 and genes[1] != gene else None,
                "in_frame": oc.normalize_in_frame(_text(it, "frame")),
                "read_count": oc.int_or_none(_text(it, "num-reads")),
                "extra": oc.to_json(extra),
                **common, **cr,
            }
            if t == "splicing-variant":
                # RNA-level exon skipping (e.g. MET exon 14 skipping)
                idx = [b["index"] for b in bps if b.get("index")]
                row["effect"] = (f"{_text(it, 'sv-type') or 'splicing'}: exon {idx[0]} -> exon {idx[-1]} junction"
                                 if len(idx) >= 2 else _text(it, "sv-type"))
                row.update({
                    "functional_effect": "exon_loss_variant",
                    "functional_effect_so": oc.SO_TERMS["exon_loss_variant"],
                    "functional_effect_raw": _text(it, "sv-type") or raw_type,
                    "functional_effect_source": "vendor",
                })
            variants.append(row)
        elif t == "expression":
            ref = it.find("reference/normal-expression/tpm")
            variants.append({
                "gene": gene, "variant_type": "expression",
                "tpm": oc.float_or_none(_text(it, "tpm")),
                "read_count": oc.int_or_none(_text(it, "num-reads")),
                "tpm_normal_n": oc.int_or_none(ref.findtext("n")) if ref is not None else None,
                "tpm_normal_mean": oc.float_or_none(ref.findtext("mean")) if ref is not None else None,
                "tpm_normal_sd": oc.float_or_none(ref.findtext("sd")) if ref is not None else None,
                "extra": oc.to_json(extra),
                **common, **cr,
            })
        else:
            warnings.append(f"alteration type {raw_type!r} is not supported; item for {gene} was not loaded")

    return oc.ParsedReport(case, variants, [], biomarkers, warnings)


def main(argv=None) -> int:
    parser = oc.base_argument_parser(LOADER, "Load GenMineTOP XML reports.")
    parser.add_argument("--panel-name", default="GenMineTOP", help="panel_name to record (default: GenMineTOP)")
    return oc.run_loader(prog=LOADER, parser=parser, suffixes=(".xml",),
                         parse_file=parse_genminetop_xml, argv=argv)


if __name__ == "__main__":
    sys.exit(main())
