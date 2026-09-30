#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
oncounify_core — shared ingestion library for the OncoUnify loaders.

Every vendor loader (load_foundation.py, load_genminetop.py,
load_guardant.py, and any third-party load_<vendor>.py) delegates the
vendor-independent work to this module, so that the normalization rules are
applied identically regardless of the source format:

* canonical protein keys (``hgvs_p``) and ``c.``-prefixed coding changes;
* functional consequence as a Sequence Ontology (SO) term, from the vendor's
  own category when one is supplied and from a single shared rule set
  otherwise (``classify_consequence``);
* controlled vocabularies for chromosome names, genome builds, CNV types,
  fusion frame, variant origin, and TMB/MSI calls;
* transactional, idempotent replacement of a report and all its child rows
  (re-running a loader never duplicates variants);
* an optional case-metadata sidecar (CSV) that lets curators supply or
  override disease labels, OncoTree codes, and patient identifiers;
* a common command-line driver that validates input, prints a per-run
  summary, and exits non-zero when any input file fails.

The module depends only on the Python standard library.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import hashlib
import json
import os
import re
import sqlite3
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple

__version__ = "2.0.0"
SCHEMA_VERSION = 2

# ---------------------------------------------------------------------------
# Canonical column lists (single source of truth for all loaders)
# ---------------------------------------------------------------------------

CASE_COLUMNS: Tuple[str, ...] = (
    "panel_name", "panel_version", "panel_type", "vendor", "report_id",
    "patient_id", "sex", "age", "date",
    "genome_build", "genome_build_source",
    "disease", "disease_ontology", "oncotree_code",
    "tissue_of_origin", "pathology_diagnosis",
    "specimen_id", "test_type", "percent_tumor_nuclei", "purity",
    "non_human_content", "other_info", "curated_fields",
    "source_file", "source_sha256", "format_version", "loader", "loaded_at",
)

VARIANT_COLUMNS: Tuple[str, ...] = (
    "gene", "variant_type", "variant_subtype",
    "chrom", "pos", "pos2", "chrom2", "ref", "alt",
    "transcript", "strand", "cds_effect", "protein_effect",
    "hgvs_c", "hgvs_p",
    "functional_effect", "functional_effect_so",
    "functional_effect_raw", "functional_effect_source",
    "status", "origin", "classification",
    "allele_fraction", "depth",
    "copy_number", "cnv_ratio", "cnv_type",
    "other_gene", "in_frame", "supporting_read_pairs",
    "tpm", "read_count", "sample_name",
    "effect", "raw_panel_type", "extra",
    "clinvar_id", "clinvar_url", "clinvar_sig", "clinvar_match",
    "clinvar_benign", "clinvar_likely_benign", "clinvar_uncertain",
    "maf_1kg", "maf_hgvd", "maf_tommo",
    "tpm_normal_n", "tpm_normal_mean", "tpm_normal_sd",
)

NON_HUMAN_COLUMNS: Tuple[str, ...] = ("organism", "reads_per_million", "status", "sample")

BIOMARKER_COLUMNS: Tuple[str, ...] = (
    "name", "value", "unit", "call", "call_raw", "assay", "source_field",
)

VARIANT_TYPES = ("short_variant", "cnv", "rearrangement", "expression")

# Fields a curator may supply or override through the case-metadata sidecar.
CURATABLE_CASE_FIELDS: Tuple[str, ...] = (
    "patient_id", "sex", "age", "date",
    "disease", "oncotree_code", "tissue_of_origin", "pathology_diagnosis",
)


# ---------------------------------------------------------------------------
# Errors and parsed-report container
# ---------------------------------------------------------------------------

class LoaderFormatError(Exception):
    """The input file does not match the vendor format this loader supports."""


@dataclass
class ParsedReport:
    case: Dict[str, Any]
    variants: List[Dict[str, Any]] = field(default_factory=list)
    non_humans: List[Dict[str, Any]] = field(default_factory=list)
    biomarkers: List[Dict[str, Any]] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Scalar helpers
# ---------------------------------------------------------------------------

_MISSING = {"", "na", "n/a", "nan", "<na>", "nat", "none", "null", "-", "."}


def clean_str(x: Any) -> Optional[str]:
    """Strip a value and map empty / NA-style sentinels to None."""
    if x is None:
        return None
    try:
        # pandas NaN
        if isinstance(x, float) and x != x:
            return None
    except Exception:
        pass
    s = str(x).strip()
    if s.lower() in _MISSING:
        return None
    return s


def float_or_none(x: Any) -> Optional[float]:
    s = clean_str(x)
    if s is None:
        return None
    try:
        v = float(s)
    except ValueError:
        return None
    if v != v:  # NaN
        return None
    return v


def int_or_none(x: Any) -> Optional[int]:
    v = float_or_none(x)
    return None if v is None else int(v)


def to_json(d: Optional[Dict[str, Any]]) -> Optional[str]:
    """Serialize an overflow dict, dropping empty values; None if nothing left."""
    if not d:
        return None
    kept = {k: v for k, v in d.items() if v is not None and v != "" and v != [] and v != {}}
    if not kept:
        return None
    return json.dumps(kept, ensure_ascii=False, sort_keys=True)


def _read_symbol_map(path: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line and not line.startswith("#") and not line.startswith("previous_symbol\t"):
            prev, cur = line.split("\t")
            out[prev] = cur
    return out


# Previous HGNC symbols of the genes in the supported assays -> current
# approved symbol, from gene_symbol_map.tsv (written by tools/make_symbol_map.py
# from the HGNC complete set).  init_db copies it into the gene_symbol_map table.
GENE_SYMBOL_MAP: Dict[str, str] = _read_symbol_map(Path(__file__).with_name("gene_symbol_map.tsv"))


def current_symbol(symbol: Optional[str]) -> Optional[str]:
    s = clean_str(symbol)
    return GENE_SYMBOL_MAP.get(s, s) if s else None


def harmonize_gene_symbols(variants: List[Dict[str, Any]]) -> int:
    """
    Replace previous HGNC symbols in gene / other_gene by the current symbol,
    keeping the vendor's symbol in the JSON column `extra`
    (vendor_gene / vendor_other_gene).  Returns the number of changed fields.
    """
    changed = 0
    for v in variants:
        notes = {}
        for col, key in (("gene", "vendor_gene"), ("other_gene", "vendor_other_gene")):
            old = clean_str(v.get(col))
            new = current_symbol(old)
            if old and new != old:
                v[col] = new
                notes[key] = old
                changed += 1
        if notes:
            extra = json.loads(v["extra"]) if v.get("extra") else {}
            extra.update(notes)
            v["extra"] = to_json(extra)
    return changed


def normalize_chrom(x: Any) -> Optional[str]:
    """Return a UCSC-style chromosome name ('chr1' .. 'chr22', 'chrX', 'chrY', 'chrM')."""
    s = clean_str(x)
    if s is None or s == "0":
        return None
    s = re.sub(r"^chr", "", s, flags=re.IGNORECASE)
    if s.upper() in ("M", "MT"):
        return "chrM"
    if s.upper() in ("X", "Y"):
        return "chr" + s.upper()
    return "chr" + s


def normalize_genome_build(x: Any) -> Optional[str]:
    s = clean_str(x)
    if s is None:
        return None
    t = s.lower().replace("_", "").replace("-", "").replace(" ", "")
    t = re.sub(r"\.?p\d+$", "", t)          # patch level, e.g. GRCh38.p14
    if t in ("hg19", "grch37", "b37", "hs37d5", "37"):
        return "GRCh37"
    if t in ("hg38", "grch38", "b38", "38"):
        return "GRCh38"
    raise ValueError(f"unrecognized genome build: {s!r}")


def parse_locus(x: Any) -> Tuple[Optional[str], Optional[int], Optional[int]]:
    """Parse 'chr12:25398284' or 'chr17:37844167-37886679' into (chrom, start, end)."""
    s = clean_str(x)
    if s is None:
        return None, None, None
    m = re.match(r"^(?:chr)?([0-9XYMT]+)[:\s]+([0-9,]+)(?:\s*-\s*([0-9,]+))?$", s, re.IGNORECASE)
    if not m:
        return None, None, None
    chrom = normalize_chrom(m.group(1))
    start = int(m.group(2).replace(",", ""))
    end = int(m.group(3).replace(",", "")) if m.group(3) else None
    return chrom, start, end


def normalize_in_frame(x: Any) -> Optional[str]:
    s = clean_str(x)
    if s is None:
        return None
    t = s.lower()
    if t in ("yes", "true", "in-frame", "inframe", "in frame", "1"):
        return "yes"
    if t in ("no", "false", "out-of-frame", "out of frame", "frameshift", "0"):
        return "no"
    return "unknown"


def normalize_origin(x: Any) -> Optional[str]:
    s = clean_str(x)
    if s is None:
        return None
    t = s.lower()
    if t.startswith("somatic"):
        return "somatic"
    if t.startswith("germline"):
        return "germline"
    return None


def normalize_cnv_type(x: Any) -> Optional[str]:
    s = clean_str(x)
    if s is None:
        return None
    t = s.lower()
    if "amplif" in t or t in ("gain", "cnv-gain"):
        return "amplification"
    if "loss" in t or "delet" in t or "homozygous" in t:
        return "deletion"
    return "other"


# ---------------------------------------------------------------------------
# Protein / coding notation
# ---------------------------------------------------------------------------

_AA3 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "Sec": "U", "Pyl": "O", "Xaa": "X",
}
_AA3_RE = re.compile("|".join(sorted(_AA3, key=len, reverse=True)))


def canonical_protein(pe: Any) -> Optional[str]:
    """
    Return a canonical protein-level search key: 'p.' + one-letter HGVS body.

      p.G12D, G12D, pG12D, p.(G12D), p.Gly12Asp   -> p.G12D
      T887Rfs*19, p.Thr887ArgfsTer19              -> p.T887Rfs*19
      W26X, p.Trp26Ter                            -> p.W26*
      A999A, p.Ala999=                            -> p.A999=
      splice site 2034+1G>A, empty, NaN           -> None

    Only the notation is normalized; the variant is not re-validated against a
    transcript model and indels are not re-justified (3' rule).  The vendor's
    verbatim string is kept separately in variants.protein_effect.
    """
    s = clean_str(pe)
    if s is None:
        return None
    if "splice" in s.lower():
        return None
    s = s.replace(" ", "")
    # 'p.' prefix; a bare lower-case 'p' before an amino-acid letter is also a
    # prefix.  An upper-case 'P' is proline and is never stripped.
    s = re.sub(r"^\(?p\.", "", s)
    s = re.sub(r"^p(?=[A-Z*(])", "", s)
    s = s.strip("()")
    s = _AA3_RE.sub(lambda m: _AA3[m.group(0)], s)
    # legacy 'X' for a stop codon
    s = re.sub(r"(?<=\d)X$", "*", s)
    s = s.replace("fsX", "fs*")
    # synonymous written as A999A -> A999=
    m = re.match(r"^([A-Z*])(\d+)([A-Z*])$", s)
    if m and m.group(1) == m.group(3):
        s = f"{m.group(1)}{m.group(2)}="
    # must look like a protein change (e.g. not 'promoter-124C>T')
    if not re.match(r"^(?:[A-Z*]\d+|[=?0]$|0\?$)", s):
        return None
    return "p." + s


def canonical_cds(cds: Any) -> Optional[str]:
    """Return the coding change with a single 'c.' prefix ('35G>A' -> 'c.35G>A')."""
    s = clean_str(cds)
    if s is None:
        return None
    s = s.replace(" ", "")
    s = re.sub(r"^c\.", "", s)
    return "c." + s


# ---------------------------------------------------------------------------
# Functional consequence (Sequence Ontology)
# ---------------------------------------------------------------------------

# SO term name -> SO accession.  Term names are the ones used by Ensembl VEP.
SO_TERMS: Dict[str, str] = {
    "missense_variant": "SO:0001583",
    "stop_gained": "SO:0001587",
    "frameshift_variant": "SO:0001589",
    "synonymous_variant": "SO:0001819",
    "stop_retained_variant": "SO:0001567",
    "start_lost": "SO:0002012",
    "stop_lost": "SO:0001578",
    "inframe_deletion": "SO:0001822",
    "inframe_insertion": "SO:0001821",
    "inframe_indel": "SO:0001820",
    "splice_donor_variant": "SO:0001575",
    "splice_acceptor_variant": "SO:0001574",
    "splice_region_variant": "SO:0001630",
    "exon_loss_variant": "SO:0001572",
    "protein_altering_variant": "SO:0001818",
    "coding_sequence_variant": "SO:0001580",
    "5_prime_UTR_variant": "SO:0001623",
    "3_prime_UTR_variant": "SO:0001624",
    "UTR_variant": "SO:0001622",
    "intron_variant": "SO:0001627",
    "upstream_gene_variant": "SO:0001631",
    "non_coding_transcript_variant": "SO:0001619",
    "sequence_variant": "SO:0001060",
}

# Display / statistics groups (used by panel_stats.cgi).
SO_GROUP: Dict[str, str] = {
    "missense_variant": "missense",
    "stop_gained": "nonsense",
    "frameshift_variant": "frameshift",
    "splice_donor_variant": "splice",
    "splice_acceptor_variant": "splice",
    "splice_region_variant": "splice",
    "exon_loss_variant": "splice",
    "inframe_deletion": "inframe",
    "inframe_insertion": "inframe",
    "inframe_indel": "inframe",
    "start_lost": "other",
    "stop_lost": "other",
    "protein_altering_variant": "other",
    "coding_sequence_variant": "other",
    "upstream_gene_variant": "other",
    "sequence_variant": "other",
    "synonymous_variant": "silent",
    "stop_retained_variant": "silent",
    "5_prime_UTR_variant": "noncoding",
    "3_prime_UTR_variant": "noncoding",
    "UTR_variant": "noncoding",
    "intron_variant": "noncoding",
    "non_coding_transcript_variant": "noncoding",
}

# Foundation Medicine functional-effect keywords -> SO (checked in this order;
# 'nonframeshift' must precede 'frameshift', which it contains).
_VENDOR_TERM_RULES: Tuple[Tuple[str, str], ...] = (
    ("nonframeshift", "@inframe"),
    ("non-frameshift", "@inframe"),
    ("frameshift", "frameshift_variant"),
    ("nonsense", "stop_gained"),
    ("stop_gained", "stop_gained"),
    ("stopgain", "stop_gained"),
    ("splice", "@splice"),
    ("inframe", "@inframe"),
    ("in-frame", "@inframe"),
    ("start", "start_lost"),
    ("stoploss", "stop_lost"),
    ("stop_lost", "stop_lost"),
    ("missense", "missense_variant"),
    ("promoter", "upstream_gene_variant"),
    ("synonymous", "synonymous_variant"),
    ("silent", "synonymous_variant"),
)

_POS_TOKEN = re.compile(r"^([-*]?)(\d+)(?:([+-])(\d+))?$")


def _parse_cds_positions(cds_body: str) -> List[Tuple[str, int, int]]:
    """
    Return [(region, base, offset), ...] for the one or two positions at the
    start of a coding HGVS body, where region is 'cds', '5utr' or '3utr' and
    offset is the signed intronic offset (0 if exonic).
    """
    m = re.match(r"^([-*]?\d+(?:[+-]\d+)?)(?:_([-*]?\d+(?:[+-]\d+)?))?", cds_body)
    if not m:
        return []
    out = []
    for tok in (m.group(1), m.group(2)):
        if not tok:
            continue
        pm = _POS_TOKEN.match(tok)
        if not pm:
            continue
        region = {"": "cds", "-": "5utr", "*": "3utr"}[pm.group(1)]
        base = int(pm.group(2))
        off = int(pm.group(4)) if pm.group(4) else 0
        if pm.group(3) == "-":
            off = -off
        out.append((region, base, off))
    return out


def _splice_class(positions: List[Tuple[str, int, int]]) -> Optional[str]:
    """Classify intronic involvement: donor / acceptor / region / intron / None."""
    if not positions:
        return None
    offs = [o for (_, _, o) in positions]
    if all(o == 0 for o in offs):
        return None
    if len(positions) == 2:
        (_, b1, o1), (_, b2, o2) = positions
        # A range that crosses an exon/intron boundary covers offsets +/-1,2.
        if (o1 == 0) != (o2 == 0):
            o = o1 if o1 != 0 else o2
            return "donor" if o > 0 else "acceptor"
        # intron -> whole exon(s) -> intron
        if o1 < 0 < o2:
            return "exon_loss"
        # c.99+5_200-5: two different introns with at least one exon between
        if o1 > 0 > o2 and b2 > b1 + 1:
            return "exon_loss"
        # otherwise both endpoints lie in the same intron: nearest offset wins
    nearest = min((o for o in offs if o != 0), key=abs)
    a = abs(nearest)
    if a <= 2:
        return "donor" if nearest > 0 else "acceptor"
    if a <= 8:
        return "region"
    return "intron"


def _indel_len_from_cds(cds_body: str, positions) -> Optional[int]:
    """Net length change (inserted - deleted) of an exonic coding indel, if computable."""
    if not positions or any(o != 0 for (_, _, o) in positions):
        return None
    span = positions[-1][1] - positions[0][1] + 1 if len(positions) == 2 else 1
    # Foundation Medicine writes deletion-insertions as substitutions:
    # 2034G>CA (G replaced by CA), 2237_2255>T (19 bases replaced by T)
    gt = re.match(r"^[-*]?\d+(?:_[-*]?\d+)?([ACGTN]*)>([ACGTN]+)$", cds_body)
    if gt:
        ref_len = len(gt.group(1)) if gt.group(1) else span
        return len(gt.group(2)) - ref_len
    if "delins" in cds_body:
        ins = re.search(r"delins([ACGTN]+)$", cds_body)
        if not ins:
            return None
        return len(ins.group(1)) - span
    if "del" in cds_body:
        return -span
    if "dup" in cds_body:
        return span
    if "ins" in cds_body:
        ins = re.search(r"ins([ACGTN]+)$", cds_body)
        if ins:
            return len(ins.group(1))
        n = re.search(r"ins(\d+)$", cds_body)
        return int(n.group(1)) if n else None
    return None


# p.E746_S752delinsV, and the Foundation Medicine form E746_S752>V
_PROT_DELINS = re.compile(r"^([A-Z*])(\d+)(?:_([A-Z*])(\d+))?(?:delins|>)([A-Z*]+)$")

_INFRAME_OR_MORE_SPECIFIC = ("inframe_deletion", "inframe_insertion", "inframe_indel",
                             "stop_gained", "stop_lost", "start_lost", "missense_variant")


def _protein_indel_class(p: str) -> Optional[str]:
    """Class of a non-frameshift protein deletion / insertion / deletion-insertion."""
    m = _PROT_DELINS.match(p)
    if m:
        first, start, end, inserted = m.group(1), int(m.group(2)), m.group(4), m.group(5)
        deleted = (int(end) - start + 1) if end else 1
        if "*" in inserted:
            return "stop_gained"
        if first == "*":
            return "stop_lost"
        if first == "M" and start == 1:
            return "start_lost"
        if len(inserted) < deleted:
            return "inframe_deletion"
        if len(inserted) > deleted:
            return "inframe_insertion"
        return "missense_variant"          # multi-residue substitution of equal length
    if "delins" in p or ">" in p:
        return "inframe_indel"
    if "del" in p:
        return "inframe_deletion"
    if "ins" in p or "dup" in p:
        return "inframe_insertion"
    return None


def _classify_from_protein(p: str) -> Optional[str]:
    """p: canonical protein body without the 'p.' prefix."""
    if "fs" in p:
        return "frameshift_variant"
    if "ext" in p:
        return "stop_lost"
    if re.match(r"^\*\d+(=|\*)$", p):
        return "stop_retained_variant"
    # stop_gained before start_lost: p.M1* is reported as nonsense (and ranks
    # above start_lost in the Ensembl VEP severity order)
    if re.match(r"^[A-Z]\d+\*$", p):
        return "stop_gained"
    if re.match(r"^\*\d+[A-Z]", p):
        return "stop_lost"
    if re.match(r"^M1(?!\d)(?!=)", p):
        return "start_lost"
    if re.match(r"^[A-Z*]\d+=$", p):
        return "synonymous_variant"
    if re.match(r"^[A-Z]\d+[A-Z]$", p):
        return "missense_variant"
    return _protein_indel_class(p)


def classify_consequence(
    hgvs_p: Optional[str],
    cds: Optional[str],
    *,
    vendor_term: Optional[str] = None,
    hint: Optional[str] = None,
    ref: Optional[str] = None,
    alt: Optional[str] = None,
) -> Tuple[str, str, str]:
    """
    Assign one SO term to a short variant.

    Returns (so_term, so_accession, source) where source is 'vendor' when the
    vendor supplied the category, 'inferred' when it was derived from the
    protein / coding notation by the shared rule set, and 'unclassified' when
    neither was informative (so_term is then 'sequence_variant').

    Rule order: vendor category -> protein change -> coding change (splice
    offsets, UTR positions, exonic indel length) -> vendor hint (e.g. Guardant
    reporting_category) -> VCF-style ref/alt length -> 'sequence_variant'.
    """
    p = hgvs_p[2:] if hgvs_p and hgvs_p.startswith("p.") else (hgvs_p or "")
    c = canonical_cds(cds)
    c_body = c[2:] if c else ""
    positions = _parse_cds_positions(c_body) if c_body else []
    splice = _splice_class(positions)

    def done(term: str, source: str) -> Tuple[str, str, str]:
        return term, SO_TERMS[term], source

    # 1. vendor-supplied category
    vt = (clean_str(vendor_term) or "").lower()
    if vt:
        for key, target in _VENDOR_TERM_RULES:
            if key in vt:
                if target == "@splice":
                    if splice in ("donor", "acceptor", "region", "exon_loss"):
                        return done({"donor": "splice_donor_variant",
                                     "acceptor": "splice_acceptor_variant",
                                     "region": "splice_region_variant",
                                     "exon_loss": "exon_loss_variant"}[splice], "vendor")
                    return done("splice_region_variant", "vendor")
                if target == "@inframe":
                    # the vendor says "not a frameshift"; refine from the notation
                    t = _classify_from_protein(p) if p else None
                    if t in _INFRAME_OR_MORE_SPECIFIC:
                        return done(t, "vendor")
                    net = _indel_len_from_cds(c_body, positions)
                    if net:
                        return done("inframe_deletion" if net < 0 else "inframe_insertion", "vendor")
                    return done("inframe_indel", "vendor")
                return done(target, "vendor")

    # 2. protein change
    if p:
        term = _classify_from_protein(p)
        if term:
            return done(term, "inferred")

    # 3. coding change
    if positions:
        if splice == "donor":
            return done("splice_donor_variant", "inferred")
        if splice == "acceptor":
            return done("splice_acceptor_variant", "inferred")
        if splice == "exon_loss":
            return done("exon_loss_variant", "inferred")
        regions = {r for (r, _, _) in positions}
        if splice is None and regions == {"5utr"}:
            return done("5_prime_UTR_variant", "inferred")
        if splice is None and regions == {"3utr"}:
            return done("3_prime_UTR_variant", "inferred")
        if splice == "region":
            return done("splice_region_variant", "inferred")
        if splice == "intron":
            return done("intron_variant", "inferred")
        net = _indel_len_from_cds(c_body, positions)
        if net is not None and net != 0:
            if net % 3:
                return done("frameshift_variant", "inferred")
            return done("inframe_deletion" if net < 0 else "inframe_insertion", "inferred")
        if regions == {"cds"}:
            return done("coding_sequence_variant", "inferred")

    # 4. vendor hint (e.g. Guardant reporting_category, GenMineTOP type)
    h = (clean_str(hint) or "").lower()
    if h:
        if "promoter" in h:
            return done("upstream_gene_variant", "inferred")
        if "utr" in h:
            return done("UTR_variant", "inferred")
        if "non_coding" in h or "noncoding" in h or "intron" in h:
            return done("non_coding_transcript_variant", "inferred")
        if "skipping" in h or "splicing-variant" in h:
            return done("exon_loss_variant", "inferred")

    # 5. VCF-style alleles of a coding variant with no other information
    r, a = clean_str(ref), clean_str(alt)
    if r and a and len(r) != len(a) and h in ("protein_coding", "coding"):
        net = len(a) - len(r)
        if net % 3:
            return done("frameshift_variant", "inferred")
        return done("inframe_deletion" if net < 0 else "inframe_insertion", "inferred")

    return done("sequence_variant", "unclassified")


# ---------------------------------------------------------------------------
# Biomarkers (TMB / MSI)
# ---------------------------------------------------------------------------

def normalize_msi_call(x: Any) -> Optional[str]:
    """Controlled MSI call: 'MSI-H', 'MSS' (MSS and MSI-L), or 'indeterminate'."""
    s = clean_str(x)
    if s is None:
        return None
    t = s.lower().replace("_", "-")
    if "not detected" in t or "negative" in t:
        return "MSS"
    if any(k in t for k in ("ambig", "indeterm", "cannot", "unknown", "equivocal", "qns", "fail")):
        return "indeterminate"
    if "msi-h" in t or "msi high" in t or "msi-high" in t or t == "high":
        return "MSI-H"
    if "mss" in t or "stable" in t or "msi-l" in t or "msi low" in t or "msi-low" in t:
        return "MSS"
    return "indeterminate"


def normalize_tmb_call(x: Any) -> Optional[str]:
    """Controlled TMB call: 'high', 'intermediate', 'low', or 'indeterminate'."""
    s = clean_str(x)
    if s is None:
        return None
    t = s.lower()
    for k in ("high", "intermediate", "low"):
        if k in t:
            return k
    return "indeterminate"


def normalize_tmb_unit(x: Any) -> Optional[str]:
    s = clean_str(x)
    if s is None:
        return None
    t = s.lower().replace(" ", "")
    if t in ("mutations-per-megabase", "muts/mb", "mut/mb", "mutations/mb", "mutationspermegabase"):
        return "mutations/Mb"
    return s


# ---------------------------------------------------------------------------
# Case-metadata sidecar (curation path)
# ---------------------------------------------------------------------------

def load_case_metadata(path: Optional[str]) -> Dict[Tuple[str, str], Dict[str, str]]:
    """
    Read a curator-maintained CSV keyed by (panel_name, report_id).

    Required columns: panel_name, report_id.  Optional columns: any of
    CURATABLE_CASE_FIELDS.  Non-empty cells override the vendor-derived value
    at load time; the overridden field names are recorded in
    cases.curated_fields for auditability.
    """
    if not path:
        return {}
    p = Path(path)
    if not p.is_file():
        raise SystemExit(f"[FATAL] case-metadata file not found: {p}")
    table: Dict[Tuple[str, str], Dict[str, str]] = {}
    with p.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        cols = set(reader.fieldnames or [])
        missing = {"panel_name", "report_id"} - cols
        if missing:
            raise SystemExit(f"[FATAL] case-metadata file {p} lacks column(s): {', '.join(sorted(missing))}")
        unknown = cols - {"panel_name", "report_id", "note", *CURATABLE_CASE_FIELDS}
        if unknown:
            print(f"[WARN] case-metadata: ignoring unknown column(s): {', '.join(sorted(unknown))}", file=sys.stderr)
        for row in reader:
            key = (clean_str(row.get("panel_name")) or "", clean_str(row.get("report_id")) or "")
            if not key[0] or not key[1]:
                continue
            table[key] = {k: row[k] for k in CURATABLE_CASE_FIELDS if k in row and clean_str(row[k]) is not None}
    return table


def apply_case_metadata(case: Dict[str, Any], table: Dict[Tuple[str, str], Dict[str, str]]) -> List[str]:
    rec = table.get((case.get("panel_name") or "", case.get("report_id") or ""))
    if not rec:
        return []
    applied = []
    for k, v in rec.items():
        val: Any = clean_str(v)
        if k == "age":
            val = int_or_none(val)
        if val is None:
            continue
        case[k] = val
        applied.append(k)
    if applied:
        case["curated_fields"] = ",".join(sorted(applied))
    return applied


# ---------------------------------------------------------------------------
# Database
# ---------------------------------------------------------------------------

def schema_path() -> Path:
    return Path(__file__).with_name("schema.sql")


def connect(db_path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON;")
    # journal_mode is left at the default (DELETE): read-only CGI readers
    # cannot create the -wal/-shm files that WAL mode would require.
    return conn


def init_db(conn: sqlite3.Connection) -> None:
    """Create the schema on an empty database; refuse to write to an old one."""
    has_cases = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='cases'"
    ).fetchone()
    version = conn.execute("PRAGMA user_version").fetchone()[0]
    if has_cases and version < SCHEMA_VERSION:
        raise SystemExit(
            f"[FATAL] database schema version {version} is older than {SCHEMA_VERSION}; "
            "run `python3 migrate_db.py <db>` first (or rebuild the database from the source reports)."
        )
    if has_cases and version > SCHEMA_VERSION:
        raise SystemExit(f"[FATAL] database schema version {version} is newer than this loader ({SCHEMA_VERSION}).")
    conn.executescript(schema_path().read_text(encoding="utf-8"))
    with conn:
        conn.execute("DELETE FROM gene_symbol_map")
        conn.executemany("INSERT INTO gene_symbol_map (previous_symbol, symbol) VALUES (?, ?)",
                         sorted(GENE_SYMBOL_MAP.items()))


def _check_keys(rows: Iterable[Dict[str, Any]], allowed: Sequence[str], what: str) -> None:
    allowed_set = set(allowed)
    for r in rows:
        bad = set(r) - allowed_set
        if bad:
            raise ValueError(f"{what}: unknown column(s) {sorted(bad)}")


def replace_report(conn: sqlite3.Connection, rep: ParsedReport) -> Tuple[int, str]:
    """
    Insert a report, or replace an existing report with the same
    (panel_name, report_id), together with all of its child rows, inside a
    single transaction.  Returns (case_id, 'inserted' | 'replaced').

    Replacing (rather than skipping) keeps the database consistent with the
    source files: re-running a loader never duplicates variants, and an
    amended vendor report or an edited sidecar entry is picked up on the next
    run.
    """
    case = rep.case
    _check_keys([case], CASE_COLUMNS, "case")
    _check_keys(rep.variants, VARIANT_COLUMNS, "variant")
    _check_keys(rep.non_humans, NON_HUMAN_COLUMNS, "non_human")
    _check_keys(rep.biomarkers, BIOMARKER_COLUMNS, "biomarker")
    for v in rep.variants:
        if v.get("variant_type") not in VARIANT_TYPES:
            raise ValueError(f"variant_type must be one of {VARIANT_TYPES}, got {v.get('variant_type')!r}")

    with conn:  # BEGIN ... COMMIT (ROLLBACK on exception)
        row = conn.execute(
            "SELECT case_id FROM cases WHERE panel_name = ? AND report_id = ?",
            (case["panel_name"], case["report_id"]),
        ).fetchone()
        cols = list(CASE_COLUMNS)
        vals = [case.get(c) for c in cols]
        if row:
            case_id = int(row[0])
            conn.execute(
                f"UPDATE cases SET {', '.join(c + ' = ?' for c in cols)} WHERE case_id = ?",
                vals + [case_id],
            )
            for t in ("variants", "non_human_contents", "biomarkers"):
                conn.execute(f"DELETE FROM {t} WHERE case_id = ?", (case_id,))
            action = "replaced"
        else:
            cur = conn.execute(
                f"INSERT INTO cases ({', '.join(cols)}) VALUES ({', '.join('?' * len(cols))})",
                vals,
            )
            case_id = int(cur.lastrowid)
            action = "inserted"

        for table, columns, rows in (
            ("variants", VARIANT_COLUMNS, rep.variants),
            ("non_human_contents", NON_HUMAN_COLUMNS, rep.non_humans),
            ("biomarkers", BIOMARKER_COLUMNS, rep.biomarkers),
        ):
            if not rows:
                continue
            sql = (f"INSERT INTO {table} (case_id, {', '.join(columns)}) "
                   f"VALUES (?, {', '.join('?' * len(columns))})")
            conn.executemany(sql, [[case_id] + [r.get(c) for c in columns] for r in rows])
    return case_id, action


# ---------------------------------------------------------------------------
# Row builders
# ---------------------------------------------------------------------------

def short_variant(
    *,
    gene: str,
    protein: Optional[str] = None,
    cds: Optional[str] = None,
    vendor_term: Optional[str] = None,
    hint: Optional[str] = None,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Build a canonical short-variant row, deriving hgvs_p/hgvs_c and the SO term."""
    hgvs_p = canonical_protein(protein)
    term, so, source = classify_consequence(
        hgvs_p, cds, vendor_term=vendor_term, hint=hint,
        ref=kwargs.get("ref"), alt=kwargs.get("alt"),
    )
    row: Dict[str, Any] = {
        "gene": gene,
        "variant_type": "short_variant",
        "cds_effect": clean_str(cds),
        "protein_effect": clean_str(protein),
        "hgvs_c": canonical_cds(cds),
        "hgvs_p": hgvs_p,
        "functional_effect": term,
        "functional_effect_so": so,
        "functional_effect_raw": clean_str(vendor_term) or clean_str(hint),
        "functional_effect_source": source,
    }
    row.update(kwargs)
    return row


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def provenance(path: Path, loader_name: str, format_version: Optional[str]) -> Dict[str, Any]:
    return {
        "source_file": path.name,
        "source_sha256": file_sha256(path),
        "format_version": format_version,
        "loader": f"{loader_name} {__version__}",
        "loaded_at": _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat(),
    }


# ---------------------------------------------------------------------------
# Command-line driver
# ---------------------------------------------------------------------------

def base_argument_parser(prog: str, description: str) -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog=prog, description=description)
    ap.add_argument("db", help="SQLite database file (created from schema.sql if absent)")
    ap.add_argument("inputs", nargs="+", help="input files and/or directories (searched recursively)")
    ap.add_argument("--case-metadata", metavar="CSV",
                    help="curator sidecar CSV keyed by panel_name,report_id (see docs/CURATION.md)")
    ap.add_argument("--genome-build", choices=("GRCh37", "GRCh38"),
                    help="override the reference assembly recorded for every report in this run")
    ap.add_argument("--fail-fast", action="store_true", help="stop at the first file that fails")
    ap.add_argument("--quiet", action="store_true", help="only print warnings, errors and the summary")
    return ap


def iter_input_files(inputs: Sequence[str], suffixes: Tuple[str, ...]) -> List[Path]:
    files: List[Path] = []
    for inp in inputs:
        p = Path(inp)
        if p.is_dir():
            for root, _dirs, names in os.walk(p):
                for n in sorted(names):
                    if n.startswith(".") or n.startswith("~$"):
                        continue
                    if n.lower().endswith(suffixes):
                        files.append(Path(root) / n)
        elif p.is_file():
            files.append(p)
        else:
            raise SystemExit(f"[FATAL] input not found: {p}")
    return sorted(dict.fromkeys(files))


def run_loader(
    *,
    prog: str,
    parser: argparse.ArgumentParser,
    suffixes: Tuple[str, ...],
    parse_file: Callable[[Path, argparse.Namespace], ParsedReport],
    argv: Optional[Sequence[str]] = None,
) -> int:
    """
    Shared main(): parse arguments, load every input file inside its own
    transaction, print a per-run summary to stderr, and return the process
    exit code (0 = all files loaded, 1 = at least one file failed,
    2 = nothing to do / fatal).
    """
    args = parser.parse_args(argv)
    files = iter_input_files(args.inputs, suffixes)
    if not files:
        print(f"[FATAL] {prog}: no input files with suffix {suffixes} found", file=sys.stderr)
        return 2

    metadata = load_case_metadata(args.case_metadata)
    conn = connect(args.db)
    init_db(conn)

    stats = {"inserted": 0, "replaced": 0, "failed": 0, "variants": 0, "biomarkers": 0, "curated": 0}
    failures: List[Tuple[Path, str]] = []
    warnings: List[Tuple[Path, str]] = []

    for path in files:
        try:
            rep = parse_file(path, args)
            if args.genome_build:
                rep.case["genome_build"] = args.genome_build
                rep.case["genome_build_source"] = "cli"
            if not rep.case.get("report_id"):
                raise LoaderFormatError("no report identifier could be determined")
            if not rep.case.get("genome_build"):
                raise LoaderFormatError("no genome build could be determined; pass --genome-build")
            if apply_case_metadata(rep.case, metadata):
                stats["curated"] += 1
            harmonize_gene_symbols(rep.variants)
            case_id, action = replace_report(conn, rep)
            stats[action] += 1
            stats["variants"] += len(rep.variants)
            stats["biomarkers"] += len(rep.biomarkers)
            for w in rep.warnings:
                warnings.append((path, w))
            if not args.quiet:
                print(f"[INFO] {action:8s} case_id={case_id} report_id={rep.case['report_id']} "
                      f"variants={len(rep.variants)} <- {path}", file=sys.stderr)
        except Exception as exc:  # recorded, reported, and reflected in the exit code
            stats["failed"] += 1
            failures.append((path, f"{type(exc).__name__}: {exc}"))
            print(f"[ERROR] {path}: {type(exc).__name__}: {exc}", file=sys.stderr)
            if args.fail_fast:
                break

    conn.close()

    for path, w in warnings:
        print(f"[WARN] {path}: {w}", file=sys.stderr)
    loaded = stats["inserted"] + stats["replaced"]
    print(
        f"[SUMMARY] {prog}: files={len(files)} loaded={loaded} "
        f"(inserted={stats['inserted']} replaced={stats['replaced']}) failed={stats['failed']} "
        f"warnings={len(warnings)} variants={stats['variants']} biomarkers={stats['biomarkers']} "
        f"curated_cases={stats['curated']}",
        file=sys.stderr,
    )
    for path, msg in failures:
        print(f"[FAILED] {path}: {msg}", file=sys.stderr)
    if stats["failed"]:
        return 1
    return 0
