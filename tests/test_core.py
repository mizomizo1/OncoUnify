#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for the shared normalization rules in oncounify_core."""

from __future__ import annotations

import json
import sqlite3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import oncounify_core as oc  # noqa: E402


class CanonicalProtein(unittest.TestCase):
    CASES = {
        "p.G12D": "p.G12D",
        "G12D": "p.G12D",
        "pG12D": "p.G12D",
        "p.(G12D)": "p.G12D",
        "p.Gly12Asp": "p.G12D",
        "T887Rfs*19": "p.T887Rfs*19",
        "p.Thr887ArgfsTer19": "p.T887Rfs*19",
        "W26X": "p.W26*",
        "p.Trp26Ter": "p.W26*",
        "A999A": "p.A999=",
        "p.Ala999=": "p.A999=",
        "P312S": "p.P312S",          # proline must not be taken for a 'p' prefix
        "E746_A750del": "p.E746_A750del",
        "I679fs*21": "p.I679fs*21",
        "splice site 904+1G>A": None,
        "promoter -124C>T": None,
        "": None,
        None: None,
        "NaN": None,
    }

    def test_cases(self):
        for raw, expected in self.CASES.items():
            with self.subTest(raw=raw):
                self.assertEqual(oc.canonical_protein(raw), expected)

    def test_canonical_cds(self):
        self.assertEqual(oc.canonical_cds("35G>A"), "c.35G>A")
        self.assertEqual(oc.canonical_cds("c.35G>A"), "c.35G>A")
        self.assertIsNone(oc.canonical_cds(" "))


class Consequence(unittest.TestCase):
    def term(self, protein=None, cds=None, **kw):
        return oc.classify_consequence(oc.canonical_protein(protein), cds, **kw)[0]

    def test_vendor_terms(self):
        self.assertEqual(self.term("G12D", "35G>A", vendor_term="missense"), "missense_variant")
        self.assertEqual(self.term("R1450*", "4348C>T", vendor_term="nonsense"), "stop_gained")
        self.assertEqual(self.term("S1982fs*22", "5946delT", vendor_term="frameshift"), "frameshift_variant")
        # 'nonframeshift' contains 'frameshift' and must not be taken for it
        self.assertEqual(self.term("E746_A750del", "2235_2249del", vendor_term="nonframeshift"), "inframe_deletion")
        self.assertEqual(self.term("splice site 904+1G>A", "904+1G>A", vendor_term="splice"), "splice_donor_variant")
        self.assertEqual(self.term("splice site 905-2A>G", "905-2A>G", vendor_term="splice"),
                         "splice_acceptor_variant")
        self.assertEqual(self.term("promoter -124C>T", "-124C>T", vendor_term="promoter"), "upstream_gene_variant")
        src = oc.classify_consequence("p.G12D", "c.35G>A", vendor_term="missense")[2]
        self.assertEqual(src, "vendor")

    def test_inferred_from_protein(self):
        # reviewer 1, T3: a Guardant frameshift was stored as missense in version 1
        self.assertEqual(self.term("G302fs*4"), "frameshift_variant")
        self.assertEqual(self.term("G12D"), "missense_variant")
        self.assertEqual(self.term("R1450*"), "stop_gained")
        self.assertEqual(self.term("A999A"), "synonymous_variant")
        self.assertEqual(self.term("M1?"), "start_lost")
        self.assertEqual(self.term("M1I"), "start_lost")
        self.assertEqual(self.term("*110Qext*17"), "stop_lost")
        self.assertEqual(self.term("E746_A750del"), "inframe_deletion")
        self.assertEqual(self.term("E746_T751delinsA"), "inframe_deletion")
        self.assertEqual(self.term("A767_V769dup"), "inframe_insertion")
        self.assertEqual(self.term("D770_N771insG"), "inframe_insertion")
        src = oc.classify_consequence("p.G12D", None)[2]
        self.assertEqual(src, "inferred")

    def test_inferred_from_cds(self):
        # reviewer 1, M4: in-frame deletions, splice-region changes without a
        # protein change, and start loss must not be left unassigned
        self.assertEqual(self.term(None, "c.905-9_905del"), "splice_acceptor_variant")
        self.assertEqual(self.term(None, "c.904_904+5del"), "splice_donor_variant")
        self.assertEqual(self.term(None, "c.2921+5G>A"), "splice_region_variant")
        self.assertEqual(self.term(None, "c.2921+14T>C"), "intron_variant")
        self.assertEqual(self.term(None, "c.99+99_100-99del"), "intron_variant")
        self.assertEqual(self.term(None, "c.99-10_200+10del"), "exon_loss_variant")
        self.assertEqual(self.term(None, "c.-124C>T"), "5_prime_UTR_variant")
        self.assertEqual(self.term(None, "c.*999G>C"), "3_prime_UTR_variant")
        self.assertEqual(self.term(None, "c.3466del"), "frameshift_variant")
        self.assertEqual(self.term(None, "c.2235_2249del"), "inframe_deletion")
        self.assertEqual(self.term(None, "c.123_124insGA"), "frameshift_variant")
        self.assertEqual(self.term(None, "c.123_125dup"), "inframe_insertion")
        self.assertEqual(self.term(None, "c.35G>A"), "coding_sequence_variant")

    def test_foundation_substitution_style_delins(self):
        # discordances found by tools/qc_report.py on institutional data
        self.assertEqual(self.term("E746_S752>V"), "inframe_deletion")
        self.assertEqual(self.term("E746_S752>V", "2236_2256>GTT", vendor_term="nonframeshift"), "inframe_deletion")
        self.assertEqual(self.term("D770_N771>DSVDN"), "inframe_insertion")
        self.assertEqual(self.term("G12_G13>VC"), "missense_variant")
        self.assertEqual(self.term("Q12_K13>*"), "stop_gained")
        self.assertEqual(self.term("Q12_K13delins*"), "stop_gained")
        self.assertEqual(self.term(None, "c.2237_2255>T"), "inframe_deletion")
        self.assertEqual(self.term(None, "c.2034G>CA"), "frameshift_variant")
        self.assertEqual(self.term(None, "c.2034G>A"), "coding_sequence_variant")

    def test_nonsense_at_start_codon_and_nonframeshift_refinement(self):
        self.assertEqual(self.term("M1*"), "stop_gained")
        self.assertEqual(self.term("M1*", "1A>T", vendor_term="nonsense"), "stop_gained")
        self.assertEqual(self.term("M1?"), "start_lost")
        self.assertEqual(self.term("*1165Yext*9", vendor_term="nonframeshift"), "stop_lost")
        self.assertEqual(self.term(None, "2310_2311insGGT", vendor_term="nonframeshift"), "inframe_insertion")
        self.assertEqual(self.term("V600_K601>E", vendor_term="nonframeshift"), "inframe_deletion")

    def test_hints_and_fallback(self):
        self.assertEqual(self.term(None, None, hint="promoter"), "upstream_gene_variant")
        self.assertEqual(self.term(None, None, hint="splicing-variant"), "exon_loss_variant")
        self.assertEqual(self.term(None, None, hint="protein_coding", ref="CA", alt="C"), "frameshift_variant")
        term, so, source = oc.classify_consequence(None, None)
        self.assertEqual((term, so, source), ("sequence_variant", "SO:0001060", "unclassified"))

    def test_gene_symbol_map_matches_schema(self):
        conn = sqlite3.connect(":memory:")
        conn.executescript((ROOT / "schema.sql").read_text(encoding="utf-8"))
        self.assertEqual(dict(conn.execute("SELECT previous_symbol, symbol FROM gene_symbol_map")), oc.GENE_SYMBOL_MAP)
        rows = [{"gene": "WHSC1L1", "other_gene": "MRE11A", "extra": '{"x": 1}'}, {"gene": "TP53"}]
        self.assertEqual(oc.harmonize_gene_symbols(rows), 2)
        self.assertEqual((rows[0]["gene"], rows[0]["other_gene"]), ("NSD3", "MRE11"))
        self.assertEqual(json.loads(rows[0]["extra"]), {"x": 1, "vendor_gene": "WHSC1L1", "vendor_other_gene": "MRE11A"})
        self.assertEqual(rows[1], {"gene": "TP53"})

    def test_every_term_has_a_group(self):
        self.assertEqual(set(oc.SO_TERMS), set(oc.SO_GROUP))

    def test_schema_vocabulary_matches_code(self):
        conn = sqlite3.connect(":memory:")
        conn.executescript((ROOT / "schema.sql").read_text(encoding="utf-8"))
        db = {t: (a, g) for t, a, g in conn.execute("SELECT term, accession, display_group FROM so_terms")}
        code = {t: (oc.SO_TERMS[t], oc.SO_GROUP[t]) for t in oc.SO_TERMS}
        self.assertEqual(db, code)


class Documentation(unittest.TestCase):
    def test_field_mapping_covers_every_canonical_column(self):
        rows = (ROOT / "docs" / "field_mapping.tsv").read_text(encoding="utf-8").splitlines()[1:]
        listed = {tuple(r.split("\t")[:2]) for r in rows if r}
        expected = ({("cases", c) for c in oc.CASE_COLUMNS} | {("variants", c) for c in oc.VARIANT_COLUMNS}
                    | {("biomarkers", c) for c in oc.BIOMARKER_COLUMNS}
                    | {("non_human_contents", c) for c in oc.NON_HUMAN_COLUMNS})
        self.assertEqual(expected - listed, set(), "columns missing from docs/field_mapping.tsv")
        self.assertEqual(listed - expected, set(), "docs/field_mapping.tsv lists unknown columns")


class ControlledVocabularies(unittest.TestCase):
    def test_msi(self):
        self.assertEqual(oc.normalize_msi_call("MSI-H"), "MSI-H")
        self.assertEqual(oc.normalize_msi_call("MSI-High"), "MSI-H")
        self.assertEqual(oc.normalize_msi_call("MSS/MSI-L"), "MSS")
        self.assertEqual(oc.normalize_msi_call("MSS"), "MSS")
        self.assertEqual(oc.normalize_msi_call("MSI-High NOT DETECTED"), "MSS")
        self.assertEqual(oc.normalize_msi_call("Cannot Be Determined"), "indeterminate")
        self.assertIsNone(oc.normalize_msi_call(""))

    def test_tmb(self):
        self.assertEqual(oc.normalize_tmb_call("High"), "high")
        self.assertEqual(oc.normalize_tmb_call("cannot be determined"), "indeterminate")
        self.assertEqual(oc.normalize_tmb_unit("mutations-per-megabase"), "mutations/Mb")

    def test_genome_build(self):
        self.assertEqual(oc.normalize_genome_build("hg19"), "GRCh37")
        self.assertEqual(oc.normalize_genome_build("hg38"), "GRCh38")
        self.assertEqual(oc.normalize_genome_build("GRCh38.p14"), "GRCh38")
        with self.assertRaises(ValueError):
            oc.normalize_genome_build("hg17")

    def test_misc(self):
        self.assertEqual(oc.normalize_chrom("12"), "chr12")
        self.assertEqual(oc.normalize_chrom("chrX"), "chrX")
        self.assertEqual(oc.normalize_chrom("MT"), "chrM")
        self.assertIsNone(oc.normalize_chrom("0"))
        self.assertEqual(oc.parse_locus("chr17:37844167-37886679"), ("chr17", 37844167, 37886679))
        self.assertEqual(oc.normalize_in_frame("Yes"), "yes")
        self.assertEqual(oc.normalize_in_frame("in-frame"), "yes")
        self.assertEqual(oc.normalize_cnv_type("loss"), "deletion")
        self.assertEqual(oc.normalize_cnv_type("cnv-amplification"), "amplification")
        self.assertEqual(oc.normalize_origin("Germline"), "germline")
        self.assertIsNone(oc.float_or_none("N/A"))
        self.assertEqual(oc.to_json({"a": None, "b": ""}), None)


if __name__ == "__main__":
    unittest.main()
