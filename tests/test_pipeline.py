#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
End-to-end tests: load the synthetic fixtures with every loader, then check
the stored rows, idempotency, error handling, the schema migration, and
(when Perl with CGI, DBI and DBD::SQLite is available) the CGI endpoints.

Run from the repository root:
    python3 -m unittest discover -s tests -v
"""

from __future__ import annotations

import json
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data"
INVALID = ROOT / "tests" / "data_invalid"
SIDECAR = DATA / "case_metadata.csv"
PY = sys.executable

# Expected content of the synthetic fixtures (report_id -> number of variant rows)
EXPECTED_VARIANTS = {
    "SYN-F1-0001": 7, "SYN-F1-0002": 6, "SYN-F1L-0003": 4,
    "SYN-GM-0001": 9, "SYN-GM-0002": 4,
    "SYN-G360-0001": 5, "SYN-G360-0002": 5,
}
EXPECTED_BUILD = {"FoundationOne": "GRCh37", "FoundationOneLiquid": "GRCh37",
                  "GenMineTOP": "GRCh38", "Guardant": "GRCh37"}


def run(*args, check=True):
    p = subprocess.run([PY, *map(str, args)], cwd=ROOT, capture_output=True, text=True)
    if check and p.returncode != 0:
        raise AssertionError(f"{args} exited {p.returncode}\nSTDERR:\n{p.stderr}")
    return p


def load_all(db: Path, *extra):
    run(ROOT / "load_foundation.py", db, DATA / "foundation", "--case-metadata", SIDECAR, "--quiet", *extra)
    run(ROOT / "load_genminetop.py", db, DATA / "genminetop", "--case-metadata", SIDECAR, "--quiet", *extra)
    run(ROOT / "load_guardant.py", db, DATA / "guardant", "--case-metadata", SIDECAR, "--quiet", *extra)
    run(ROOT / "load_panel_genes.py", db, DATA / "panels", "--allow-synthetic")


class Pipeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="oncounify-test-"))
        cls.db = cls.tmp / "panels.db"
        load_all(cls.db)
        cls.conn = sqlite3.connect(cls.db)

    @classmethod
    def tearDownClass(cls):
        cls.conn.close()
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def q(self, sql, *args):
        return self.conn.execute(sql, args).fetchall()

    # -- structure ---------------------------------------------------------
    def test_schema_version_and_integrity(self):
        self.assertEqual(self.q("PRAGMA user_version")[0][0], 2)
        self.assertEqual(self.q("PRAGMA integrity_check")[0][0], "ok")
        self.assertEqual(self.q("PRAGMA foreign_key_check"), [])

    def test_variant_counts_per_report(self):
        got = dict(self.q("SELECT c.report_id, COUNT(v.variant_id) FROM cases c "
                          "LEFT JOIN variants v USING(case_id) GROUP BY c.report_id"))
        self.assertEqual(got, EXPECTED_VARIANTS)

    # -- reviewer 1, M2 / reviewer 3: genome build ------------------------------
    def test_every_report_has_a_genome_build(self):
        for panel, build, source in self.q("SELECT panel_name, genome_build, genome_build_source FROM cases"):
            self.assertEqual(build, EXPECTED_BUILD[panel], panel)
            self.assertIn(source, ("vendor-file", "loader-default"))
        self.assertEqual(self.q("SELECT genome_build_source FROM cases WHERE panel_name='GenMineTOP' LIMIT 1")[0][0],
                         "vendor-file")

    # -- reviewer 1, M4 / T3: one controlled consequence vocabulary ------------
    def test_every_short_variant_has_an_so_term(self):
        rows = self.q("SELECT v.functional_effect, v.functional_effect_so, v.functional_effect_source, t.term "
                      "FROM variants v LEFT JOIN so_terms t ON t.term = v.functional_effect "
                      "WHERE v.variant_type = 'short_variant'")
        self.assertTrue(rows)
        for fe, so, source, known in rows:
            self.assertIsNotNone(known, f"{fe} not in so_terms")
            self.assertTrue(so.startswith("SO:"))
            self.assertIn(source, ("vendor", "inferred", "unclassified"))

    def test_consequences_of_interest(self):
        got = dict(((r, g, p), fe) for r, g, p, fe in self.q(
            "SELECT c.report_id, v.gene, COALESCE(v.hgvs_p, v.hgvs_c), v.functional_effect "
            "FROM variants v JOIN cases c USING(case_id) WHERE v.functional_effect IS NOT NULL"))
        self.assertEqual(got[("SYN-F1-0002", "EGFR", "p.E746_A750del")], "inframe_deletion")
        self.assertEqual(got[("SYN-F1-0001", "SMAD4", "c.904+1G>A")], "splice_donor_variant")
        self.assertEqual(got[("SYN-GM-0001", "SMAD4", "c.905-9_905del")], "splice_acceptor_variant")
        self.assertEqual(got[("SYN-GM-0001", "BRCA1", "p.Q1756Pfs*74")], "frameshift_variant")
        self.assertEqual(got[("SYN-G360-0001", "APC", "c.3466del")], "frameshift_variant")
        self.assertEqual(got[("SYN-G360-0002", "EGFR", "p.E746_A750del")], "inframe_deletion")
        skip = self.q("SELECT variant_type, functional_effect, pos, pos2 FROM variants "
                      "WHERE variant_subtype = 'splicing-variant'")
        self.assertEqual(skip, [("rearrangement", "exon_loss_variant", 116771976, 116774880)])

    def test_protein_keys_are_canonical_across_vendors(self):
        rows = self.q("SELECT c.panel_name FROM variants v JOIN cases c USING(case_id) "
                      "WHERE v.gene = 'KRAS' AND v.hgvs_p = 'p.G12D' ORDER BY 1")
        self.assertEqual([r[0] for r in rows], ["FoundationOne", "GenMineTOP", "Guardant"])
        verbatim = dict(self.q("SELECT c.panel_name, v.protein_effect FROM variants v JOIN cases c USING(case_id) "
                               "WHERE v.gene = 'KRAS'"))
        self.assertEqual(verbatim["GenMineTOP"], "p.G12D")   # vendor string kept verbatim
        self.assertEqual(verbatim["FoundationOne"], "G12D")

    # -- reviewer 1, M7 / T1: biomarkers with assay, never pooled ------------
    def test_biomarkers(self):
        rows = self.q("SELECT c.panel_name, b.name, b.value, b.unit, b.call, b.assay FROM biomarkers b "
                      "JOIN cases c USING(case_id)")
        for panel, name, value, unit, call, assay in rows:
            self.assertTrue(assay)
            if panel == "Guardant":
                self.assertEqual(name, "MSI", "Guardant MSI score must never be stored as TMB")
        g = self.q("SELECT c.report_id, b.value, b.call FROM biomarkers b JOIN cases c USING(case_id) "
                   "WHERE c.panel_name = 'Guardant' ORDER BY c.report_id")
        self.assertEqual(g, [("SYN-G360-0001", 1.0, "MSS"), ("SYN-G360-0002", 24.0, "MSI-H")])
        assays = {a for (a,) in self.q("SELECT DISTINCT assay FROM biomarkers WHERE name = 'TMB'")}
        self.assertEqual(len(assays), 3)  # tissue TMB, bTMB and GenMineTOP are kept apart

    # -- Guardant specifics ------------------------------------------------
    def test_guardant_uncalled_rows_are_not_loaded(self):
        self.assertEqual(self.q("SELECT COUNT(*) FROM variants WHERE status = 'not_called'")[0][0], 0)
        genes = {g for (g,) in self.q("SELECT v.gene FROM variants v JOIN cases c USING(case_id) "
                                      "WHERE c.report_id = 'SYN-G360-0001'")}
        self.assertEqual(genes, {"KRAS", "APC", "TP53", "ERBB2"})

    def test_guardant_include_uncalled(self):
        db = self.tmp / "uncalled.db"
        run(ROOT / "load_guardant.py", db, DATA / "guardant", "--include-uncalled", "--quiet")
        n = sqlite3.connect(db).execute("SELECT COUNT(*) FROM variants WHERE status = 'not_called'").fetchone()[0]
        self.assertEqual(n, 6)   # 3 SNV, 1 indel, 1 CNA and 1 fusion candidate with call = 0

    def test_guardant_fusion_both_sheet_spellings_and_af(self):
        rows = self.q("SELECT v.gene, v.other_gene, v.chrom, v.pos, v.pos2, v.allele_fraction "
                      "FROM variants v JOIN cases c USING(case_id) "
                      "WHERE c.panel_name='Guardant' AND v.variant_type='rearrangement'")
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0][:5], ("EML4", "ALK", "chr2", 42522656, 29446394))
        self.assertAlmostEqual(rows[0][5], 0.0061)
        af = self.q("SELECT allele_fraction FROM variants WHERE gene='KRAS' AND chrom='chr12' AND pos=25398284 "
                    "AND case_id=(SELECT case_id FROM cases WHERE report_id='SYN-G360-0001')")[0][0]
        self.assertAlmostEqual(af, 0.0421)

    # -- curation sidecar (reviewer 1, M6) ----------------------------------
    def test_sidecar_curation(self):
        row = self.q("SELECT patient_id, oncotree_code, curated_fields FROM cases WHERE report_id='SYN-F1-0001'")[0]
        self.assertEqual(row[:2], ("SYNPT-0001", "COAD"))
        self.assertIn("patient_id", row[2])
        g = self.q("SELECT disease, tissue_of_origin FROM cases WHERE report_id='SYN-G360-0001'")[0]
        self.assertEqual(g, ("Colon adenocarcinoma", "Colon"))

    # -- reviewer 1, T2: re-running never duplicates -----------------------
    def test_rerun_is_idempotent(self):
        db = self.tmp / "rerun.db"
        load_all(db)
        before = sqlite3.connect(db).execute(
            "SELECT (SELECT COUNT(*) FROM cases), (SELECT COUNT(*) FROM variants), "
            "(SELECT COUNT(*) FROM biomarkers), (SELECT COUNT(*) FROM non_human_contents), "
            "(SELECT COUNT(*) FROM panel_genes)").fetchone()
        load_all(db)
        after = sqlite3.connect(db).execute(
            "SELECT (SELECT COUNT(*) FROM cases), (SELECT COUNT(*) FROM variants), "
            "(SELECT COUNT(*) FROM biomarkers), (SELECT COUNT(*) FROM non_human_contents), "
            "(SELECT COUNT(*) FROM panel_genes)").fetchone()
        self.assertEqual(before, after)
        p = run(ROOT / "load_foundation.py", db, DATA / "foundation", "--quiet")
        self.assertIn("replaced=3", p.stderr)

    # -- reviewer 1, minor 5 / reviewer 3: failures are visible ---------------
    def test_invalid_inputs_fail_loudly(self):
        db = self.tmp / "invalid.db"
        p = run(ROOT / "load_foundation.py", db, INVALID / "foundation", "--quiet", check=False)
        self.assertEqual(p.returncode, 1)
        self.assertIn("failed=2", p.stderr)
        p = run(ROOT / "load_guardant.py", db, INVALID / "guardant", "--quiet", check=False)
        self.assertEqual(p.returncode, 1)
        self.assertIn("does not match --filename-pattern", p.stderr)
        self.assertIn("required sheet 'indels' is missing", p.stderr)
        self.assertEqual(sqlite3.connect(db).execute("SELECT COUNT(*) FROM cases").fetchone()[0], 0)

    def test_filename_pattern_is_configurable(self):
        db = self.tmp / "pattern.db"
        run(ROOT / "load_guardant.py", db, DATA / "guardant" / "Interim_SYNPT-0001_SYN-G360-0001.xlsx",
            "--filename-pattern", r"^Interim_(?P<report_id>[^_]+)_(?P<patient_id>[^_]+)$", "--quiet")
        row = sqlite3.connect(db).execute("SELECT report_id, patient_id FROM cases").fetchone()
        self.assertEqual(row, ("SYNPT-0001", "SYN-G360-0001"))

    def test_default_filename_pattern_is_patient_then_report(self):
        row = self.q("SELECT report_id, patient_id FROM cases WHERE panel_name = 'Guardant' ORDER BY report_id")
        self.assertEqual(row, [("SYN-G360-0001", "SYNPT-0001"), ("SYN-G360-0002", "SYNPT-0004")])

    # -- gene symbols: previous HGNC symbols are harmonized -------------------------
    def test_gene_symbols_are_harmonized(self):
        row = self.q("SELECT gene, extra FROM variants WHERE gene = 'NSD3'")
        self.assertEqual(len(row), 1)
        self.assertEqual(json.loads(row[0][1])["vendor_gene"], "WHSC1L1")
        self.assertEqual(self.q("SELECT COUNT(*) FROM variants WHERE gene = 'WHSC1L1'")[0][0], 0)
        listed = {g for (g,) in self.q("SELECT gene FROM panel_genes")}
        self.assertIn("NSD3", listed)
        self.assertNotIn("WHSC1L1", listed)
        tested = dict(self.q("SELECT gene, COUNT(DISTINCT case_id) FROM v_case_genes_tested GROUP BY gene"))
        self.assertEqual(tested["NSD3"], 2)          # both FoundationOne CDx reports

    # -- reviewer 3, major 8: reliability of rule-based consequences ------------
    def test_qc_report(self):
        p = run(ROOT / "tools" / "qc_report.py", self.db, "--json")
        r = json.loads(p.stdout)
        i = r["inference_vs_vendor"]
        self.assertEqual(i["variants_with_vendor_category"], 13)
        self.assertEqual(i["same_term"], 12)   # TERT c.-124C>T: vendor 'promoter', rules '5_prime_UTR_variant'
        g = r["gene_lists"]
        self.assertEqual((g["reports_with_gene_list"], g["reports"]), (7, 7))
        self.assertEqual(g["synthetic_lists_loaded"], 4)
        self.assertEqual(r["genes_outside_gene_list"], [])

    # -- reviewer 1, M5: denominators -----------------------------------------
    def test_panel_denominators(self):
        tested = dict(self.q("SELECT gene, COUNT(DISTINCT case_id) FROM v_case_genes_tested "
                             "WHERE short_variants = 1 GROUP BY gene"))
        self.assertEqual(tested["TP53"], 7)
        self.assertEqual(tested["SMAD4"], 4)   # not in the synthetic Guardant / FoundationOne Liquid lists


class GeneLists(unittest.TestCase):
    """Official lists load; synthetic test lists are refused unless explicitly allowed."""

    def test_official_and_synthetic_lists(self):
        tmp = Path(tempfile.mkdtemp(prefix="oncounify-panels-"))
        self.addCleanup(shutil.rmtree, tmp, True)
        db = tmp / "p.db"
        p = run(ROOT / "load_panel_genes.py", db, DATA / "panels", check=False)
        self.assertEqual(p.returncode, 1)
        self.assertIn("SYNTHETIC", p.stderr)
        conn = sqlite3.connect(db)
        self.assertEqual(conn.execute("SELECT COUNT(*) FROM panel_versions").fetchone()[0], 0)
        run(ROOT / "load_panel_genes.py", db, ROOT / "panels")
        counts = dict(conn.execute("SELECT pv.panel_name || ' / ' || pv.panel_version, SUM(pg.short_variants) "
                                   "FROM panel_versions pv JOIN panel_genes pg USING(panel_version_id) GROUP BY 1"))
        self.assertEqual(counts, {"FoundationOne / FoundationOneDx": 311, "FoundationOneLiquid / FoundationOneLiquidDx": 311,
                                  "GenMineTOP / TDv1.1.0+TRv6.4.3": 737, "GenMineTOP / TDv6+TRv6": 737,
                                  "Guardant / Guardant360 CDx": 74})
        conn.close()


class Migration(unittest.TestCase):
    """A version-1 database (as written by the original loaders) upgrades cleanly."""

    def test_v1_to_v2(self):
        tmp = Path(tempfile.mkdtemp(prefix="oncounify-mig-"))
        try:
            db = tmp / "v1.db"
            c = sqlite3.connect(db)
            c.executescript((ROOT / "tests" / "legacy" / "schema_v1.sql").read_text(encoding="utf-8"))
            c.execute("INSERT INTO cases (panel_name, panel_type, vendor, report_id, msi_status, tmb_score, test_type) "
                      "VALUES ('Guardant','Guardant','Guardant','G1','MSS/MSI-L',1.0,'Guardant')")
            c.execute("INSERT INTO cases (panel_name, panel_type, vendor, report_id, msi_status, tmb_score, "
                      "tmb_status, tmb_unit, test_type) VALUES ('FoundationOne','FoundationOneDx','FoundationMedicine',"
                      "'F1','MSS',5.0,'low','mutations-per-megabase','FoundationOneDx')")
            for _ in range(2):  # version-1 re-runs appended duplicates
                c.execute("INSERT INTO variants (case_id, gene, variant_type, variant_subtype, chrom, pos, "
                          "protein_effect, functional_effect) VALUES (1,'KRAS','short_variant','SNV','12',25398284,"
                          "'G302fs*4','missense')")
                c.execute("INSERT INTO variants (case_id, gene, variant_type, protein_effect, cds_effect, "
                          "functional_effect) VALUES (2,'EGFR','short_variant','E746_A750del',"
                          "'2235_2249del','nonframeshift')")
            c.commit()
            c.close()
            p = run(ROOT / "migrate_db.py", db)
            self.assertIn("duplicates removed=2", p.stderr)
            m = sqlite3.connect(db)
            self.assertEqual(m.execute("PRAGMA user_version").fetchone()[0], 2)
            self.assertEqual(m.execute("SELECT COUNT(*) FROM variants").fetchone()[0], 2)
            self.assertEqual(m.execute("SELECT chrom, functional_effect FROM variants WHERE gene='KRAS'").fetchone(),
                             ("chr12", "frameshift_variant"))
            self.assertEqual(m.execute("SELECT functional_effect FROM variants WHERE gene='EGFR'").fetchone()[0],
                             "inframe_deletion")
            self.assertEqual(m.execute("SELECT name, value, call FROM biomarkers WHERE case_id=1").fetchall(),
                             [("MSI", 1.0, "MSS")])
            self.assertEqual(m.execute("SELECT genome_build FROM cases ORDER BY case_id").fetchall(),
                             [("GRCh37",), ("GRCh37",)])
            self.assertTrue(list(tmp.glob("v1.db.v1-backup-*")))
            # a version-2 loader now writes to the migrated file
            run(ROOT / "load_foundation.py", db, DATA / "foundation", "--quiet")
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_loader_refuses_v1_database(self):
        tmp = Path(tempfile.mkdtemp(prefix="oncounify-mig-"))
        try:
            db = tmp / "v1.db"
            sqlite3.connect(db).executescript((ROOT / "tests" / "legacy" / "schema_v1.sql").read_text())
            p = run(ROOT / "load_foundation.py", db, DATA / "foundation", check=False)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("migrate_db.py", p.stderr)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


def perl_ok() -> bool:
    if not shutil.which("perl"):
        return False
    p = subprocess.run(["perl", "-MCGI", "-MDBI", "-MDBD::SQLite", "-MJSON::PP", "-e", "1"], capture_output=True)
    return p.returncode == 0


@unittest.skipUnless(perl_ok(), "Perl with CGI, DBI, DBD::SQLite and JSON::PP is required for CGI tests")
class CGI(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="oncounify-cgi-"))
        cls.db = cls.tmp / "panels.db"
        load_all(cls.db)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def cgi(self, script, query):
        env = dict(os.environ, ONCOUNIFY_DB=str(self.db), REQUEST_METHOD="GET", QUERY_STRING=query)
        p = subprocess.run(["perl", "-w", str(ROOT / script)], env=env, capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(p.stderr.strip(), "", f"Perl warnings: {p.stderr}")
        # text=True applies universal newlines, so the CRLF header terminator arrives as "\n\n"
        head, _, body = p.stdout.partition("\n\n")
        return head, body

    def tsv(self, query):
        head, body = self.cgi("panel_search.cgi", query + "&download=tsv")
        self.assertIn("text/tab-separated-values", head)
        lines = [l for l in body.splitlines() if l]
        cols = lines[0].split("\t")
        return [dict(zip(cols, l.split("\t"))) for l in lines[1:]]

    def test_flagship_query(self):
        # "KRAS p.G12D, any vendor, diagnosis containing 'colon'" (Introduction)
        rows = self.tsv("gene=KRAS&protein_effect=p.Gly12Asp&disease=colon&view_mode=case")
        self.assertEqual(sorted(r["panel_name"] for r in rows), ["FoundationOne", "GenMineTOP", "Guardant"])
        rows = self.tsv("gene=KRAS&protein_effect=G12D&disease=colon&view_mode=patient")
        self.assertEqual(sorted((r["patient_group"], r["case_count"]) for r in rows),
                         [("SYNPT-0001", "2"), ("SYNPT-0002", "1")])

    def test_gene_is_exact_by_default(self):
        self.assertEqual({r["gene"] for r in self.tsv("gene=AR")}, {"AR"})
        self.assertIn("ARID1A", {r["gene"] for r in self.tsv("gene=AR&gene_match=substring")})
        partners = {(r["gene"], r["other_gene"]) for r in self.tsv("gene=ALK")}
        self.assertEqual(partners, {("EML4", "ALK")})
        # a previous HGNC symbol finds the harmonized row
        self.assertEqual({r["gene"] for r in self.tsv("gene=WHSC1L1")}, {"NSD3"})

    def test_protein_query_forms(self):
        self.assertEqual(len(self.tsv("protein_effect=G12")), 3)             # residue query
        self.assertEqual(len(self.tsv("protein_effect=S1982Rfs*22")), 1)     # any frameshift at S1982
        self.assertEqual(len(self.tsv("protein_effect=G12Dfs")), 0)

    def test_like_metacharacters_are_literal(self):
        self.assertEqual(self.tsv("patient_id=%25"), [])
        self.assertEqual(self.tsv("disease=_"), [])

    def test_tsv_has_build_and_so(self):
        r = self.tsv("gene=KRAS&panel_name=GenMineTOP")[0]
        self.assertEqual((r["genome_build"], r["functional_effect"], r["functional_effect_so"]),
                         ("GRCh38", "missense_variant", "SO:0001583"))

    def test_html_pages(self):
        for script, query in (("panel_search.cgi", "gene=TP53&view_mode=variant"),
                              ("panel_search.cgi", "gene=TP53&view_mode=patient"),
                              ("panel_stats.cgi", ""), ("panel_stats.cgi", "src=oncotree&organ=COAD&genes=common"),
                              ("case_detail.cgi", "case_id=1")):
            head, body = self.cgi(script, query)
            self.assertIn("text/html", head)
            self.assertIn('<html lang="en">', body)
        head, body = self.cgi("panel_stats.cgi", "")
        self.assertIn("2 / 4", body)       # SMAD4 counted against the four reports whose assay covers it

    def test_case_detail_escapes_output(self):
        db = sqlite3.connect(self.db)
        db.execute("UPDATE cases SET pathology_diagnosis = '<script>alert(1)</script>' WHERE case_id = 2")
        db.commit()
        db.close()
        _, body = self.cgi("case_detail.cgi", "case_id=2")
        self.assertNotIn("<script>alert(1)</script>", body)
        self.assertIn("&lt;script&gt;", body)

    def test_suggest(self):
        _, body = self.cgi("suggest.cgi", "type=gene&q=AL")
        self.assertIn("ALK", json.loads(body)["items"])
        _, body = self.cgi("suggest.cgi", "type=protein&q=Gly12")
        self.assertEqual(json.loads(body)["items"], ["p.G12D"])
        _, body = self.cgi("suggest.cgi", "type=panel")
        self.assertEqual(len(json.loads(body)["items"]), 4)


if __name__ == "__main__":
    unittest.main()
