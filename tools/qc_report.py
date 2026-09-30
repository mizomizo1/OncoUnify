#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Registry quality report for an OncoUnify database.

    python3 tools/qc_report.py panels.db [--json]

    python3 tools/qc_report.py panels.db --examples 3     # show discordant notations
    python3 tools/qc_report.py panels.db --benchmark      # time the use-case queries

Reports, per assay: number of reports and patients, genome builds and their
source, variant rows by class, how short-variant consequences were assigned
(vendor / inferred / unclassified), biomarker assays, curation coverage,
gene-list coverage, Guardant candidate rows that were not loaded — and a
reliability estimate for the rule-based consequence inference: for every
short variant whose consequence came from the vendor (Foundation Medicine
supplies one), the term derived from the vendor category and the term
inferred from the protein and coding notation alone (as for a vendor that
supplies none) are both recomputed with the current rules and compared.
No identifier is printed.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sqlite3
import statistics
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import oncounify_core as oc  # noqa: E402


def q(conn, sql, *args):
    return conn.execute(sql, args).fetchall()


def report(db: str, examples: int = 0) -> dict:
    conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    out: dict = {}
    out["reports"] = [dict(zip(("panel_name", "panel_version", "genome_build", "genome_build_source", "n"), r))
                      for r in q(conn, "SELECT panel_name, panel_version, genome_build, genome_build_source, COUNT(*) "
                                       "FROM cases GROUP BY 1, 2, 3, 4 ORDER BY 1, 2")]
    out["variant_rows"] = [dict(zip(("panel_name", "variant_type", "n"), r))
                           for r in q(conn, "SELECT c.panel_name, v.variant_type, COUNT(*) FROM variants v "
                                            "JOIN cases c USING(case_id) GROUP BY 1, 2 ORDER BY 1, 2")]
    out["consequence_source"] = [dict(zip(("panel_name", "source", "n"), r))
                                 for r in q(conn, "SELECT c.panel_name, v.functional_effect_source, COUNT(*) "
                                                  "FROM variants v JOIN cases c USING(case_id) "
                                                  "WHERE v.variant_type = 'short_variant' GROUP BY 1, 2 ORDER BY 1, 2")]
    out["biomarkers"] = [dict(zip(("assay", "name", "n", "with_value", "with_call"), r))
                         for r in q(conn, "SELECT assay, name, COUNT(*), COUNT(value), COUNT(call) FROM biomarkers "
                                          "GROUP BY 1, 2 ORDER BY 1, 2")]
    n_cases = q(conn, "SELECT COUNT(*) FROM cases")[0][0]
    out["patients"] = {
        "distinct_patient_ids": q(conn, "SELECT COUNT(DISTINCT patient_id) FROM cases WHERE patient_id IS NOT NULL")[0][0],
        "reports_without_patient_id": q(conn, "SELECT COUNT(*) FROM cases WHERE patient_id IS NULL")[0][0],
        "patients_with_reports_from_2_or_more_assays": q(
            conn, "SELECT COUNT(*) FROM (SELECT patient_id FROM cases WHERE patient_id IS NOT NULL "
                  "GROUP BY patient_id HAVING COUNT(DISTINCT panel_name) >= 2)")[0][0],
    }
    skipped = 0
    for (info,) in q(conn, "SELECT other_info FROM cases WHERE panel_name LIKE 'Guardant%' AND other_info IS NOT NULL"):
        try:
            skipped += int(json.loads(info).get("uncalled_rows_skipped") or 0)
        except (ValueError, TypeError):
            pass
    out["guardant"] = {
        "reports": q(conn, "SELECT COUNT(*) FROM cases WHERE panel_name LIKE 'Guardant%'")[0][0],
        "variant_rows_loaded": q(conn, "SELECT COUNT(*) FROM variants v JOIN cases c USING(case_id) "
                                       "WHERE c.panel_name LIKE 'Guardant%'")[0][0],
        "uncalled_candidate_rows_not_loaded": skipped,
    }
    out["curation"] = {
        "reports": n_cases,
        "with_disease_label": q(conn, "SELECT COUNT(*) FROM cases WHERE COALESCE(disease, pathology_diagnosis, "
                                      "tissue_of_origin) IS NOT NULL")[0][0],
        "with_oncotree_code": q(conn, "SELECT COUNT(*) FROM cases WHERE oncotree_code IS NOT NULL")[0][0],
        "with_patient_id": q(conn, "SELECT COUNT(*) FROM cases WHERE patient_id IS NOT NULL")[0][0],
        "with_curated_fields": q(conn, "SELECT COUNT(*) FROM cases WHERE curated_fields IS NOT NULL")[0][0],
    }
    out["gene_lists"] = {
        "reports_with_gene_list": q(conn, "SELECT COUNT(*) FROM cases c WHERE EXISTS (SELECT 1 FROM panel_versions pv "
                                          "WHERE pv.panel_name = c.panel_name AND pv.panel_version = c.panel_version)")[0][0],
        "reports": n_cases,
    }

    # short variants in genes that are not in the assay's gene list (outdated list or symbol mismatch)
    outside = q(conn, """
        SELECT c.panel_name, c.panel_version, v.gene, COUNT(*) FROM variants v JOIN cases c USING(case_id)
        JOIN panel_versions pv ON pv.panel_name = c.panel_name AND pv.panel_version = c.panel_version
        WHERE v.variant_type = 'short_variant' AND (v.status IS NULL OR v.status != 'not_called')
          AND NOT EXISTS (SELECT 1 FROM panel_genes pg WHERE pg.panel_version_id = pv.panel_version_id
                          AND pg.gene = v.gene AND pg.short_variants = 1)
        GROUP BY 1, 2, 3 ORDER BY 4 DESC, 3""")
    out["genes_outside_gene_list"] = [dict(zip(("panel_name", "panel_version", "gene", "n"), r)) for r in outside]
    out["gene_symbols_harmonized"] = q(conn, "SELECT COUNT(*) FROM variants WHERE extra LIKE '%\"vendor_gene\"%' "
                                             "OR extra LIKE '%\"vendor_other_gene\"%'")[0][0]

    # reliability of rule-based inference against vendor-supplied categories
    group = {t: g for t, g in q(conn, "SELECT term, display_group FROM so_terms")}
    rows = q(conn, "SELECT v.hgvs_p, v.cds_effect, v.functional_effect_raw, v.protein_effect FROM variants v "
                   "WHERE v.variant_type = 'short_variant' AND v.functional_effect_source = 'vendor'")
    pairs = Counter()
    shown = {}
    for hgvs_p, cds, raw, protein in rows:
        vendor_term = oc.classify_consequence(hgvs_p, cds, vendor_term=raw)[0]
        inferred = oc.classify_consequence(hgvs_p, cds)[0]
        pairs[(vendor_term, inferred)] += 1
        if vendor_term != inferred and examples:
            ex = shown.setdefault((vendor_term, inferred), [])
            if len(ex) < examples and (raw, protein, cds) not in ex:
                ex.append((raw, protein, cds))
    n = sum(pairs.values())
    term_agree = sum(k for (a, b), k in pairs.items() if a == b)
    group_agree = sum(k for (a, b), k in pairs.items() if group.get(a) == group.get(b))
    out["inference_vs_vendor"] = {
        "variants_with_vendor_category": n,
        "same_term": term_agree,
        "same_display_group": group_agree,
        "same_term_pct": round(100 * term_agree / n, 1) if n else None,
        "same_group_pct": round(100 * group_agree / n, 1) if n else None,
        "discordant": [{"vendor": a, "inferred": b, "n": k,
                        "examples": [{"vendor_category": e[0], "protein_effect": e[1], "cds_effect": e[2]}
                                     for e in shown.get((a, b), [])]}
                       for (a, b), k in sorted(pairs.items(), key=lambda x: -x[1]) if a != b],
    }
    conn.close()
    return out


def _median_ms(fn, repeat=5):
    times = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn()
        times.append((time.perf_counter() - t0) * 1000)
    return round(statistics.median(times), 1)


def benchmark(db: str) -> dict:
    """Time the three use-case queries: SQL only, and end to end through the CGI scripts."""
    conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    pat = q(conn, "SELECT patient_id FROM cases WHERE patient_id IS NOT NULL GROUP BY patient_id "
                  "ORDER BY COUNT(*) DESC LIMIT 1")
    patient = pat[0][0] if pat else ""
    sql = {
        "recruitment (KRAS p.G12D, disease ~ colon, per patient)": (
            "SELECT COALESCE(c.patient_id, c.case_id), COUNT(DISTINCT c.case_id) FROM variants v "
            "JOIN cases c ON c.case_id = v.case_id WHERE (v.gene IN ('KRAS') OR v.other_gene IN ('KRAS')) "
            "AND v.hgvs_p = 'p.G12D' AND (c.disease LIKE '%colon%' OR c.disease_ontology LIKE '%colon%' "
            "OR c.oncotree_code LIKE '%colon%' OR c.tissue_of_origin LIKE '%colon%' "
            "OR c.pathology_diagnosis LIKE '%colon%') GROUP BY 1", ()),
        "longitudinal (all variants of the most-tested patient)": (
            "SELECT c.case_id, v.gene, v.hgvs_p FROM variants v JOIN cases c ON c.case_id = v.case_id "
            "WHERE c.patient_id = ?", (patient,)),
        "institutional (panel-aware gene frequencies, all reports)": (
            "WITH tested AS (SELECT case_id, gene FROM v_case_genes_tested WHERE short_variants = 1), "
            "hits AS (SELECT DISTINCT v.case_id, v.gene FROM variants v JOIN tested t ON t.case_id = v.case_id "
            "AND t.gene = v.gene WHERE v.functional_effect IN (SELECT term FROM so_terms WHERE display_group "
            "NOT IN ('silent','noncoding'))) SELECT t.gene, COUNT(DISTINCT t.case_id), COUNT(DISTINCT h.case_id) "
            "FROM tested t LEFT JOIN hits h ON h.case_id = t.case_id AND h.gene = t.gene GROUP BY t.gene "
            "ORDER BY 3 DESC LIMIT 20", ()),
    }
    out = {"database_bytes": os.path.getsize(db),
           "rows": {t: q(conn, f"SELECT COUNT(*) FROM {t}")[0][0] for t in ("cases", "variants", "biomarkers")},
           "sql_median_ms": {k: _median_ms(lambda s=s, a=a: conn.execute(s, a).fetchall()) for k, (s, a) in sql.items()}}
    conn.close()

    root = Path(__file__).resolve().parents[1]
    if shutil.which("perl") and subprocess.run(["perl", "-MCGI", "-MDBI", "-MDBD::SQLite", "-e", "1"],
                                               capture_output=True).returncode == 0:
        cgi = {
            "recruitment": ("panel_search.cgi", "gene=KRAS&protein_effect=G12D&disease=colon&view_mode=patient"),
            "longitudinal": ("panel_search.cgi", "patient_id=" + patient + "&view_mode=case"),
            "institutional": ("panel_stats.cgi", ""),
        }
        res = {}
        for name, (script, query) in cgi.items():
            env = dict(os.environ, ONCOUNIFY_DB=os.path.abspath(db), REQUEST_METHOD="GET", QUERY_STRING=query)
            res[name] = _median_ms(lambda script=script, env=env: subprocess.run(
                ["perl", str(root / script)], env=env, capture_output=True, check=True), repeat=3)
        out["cgi_end_to_end_median_ms"] = res
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Quality report for an OncoUnify database.")
    ap.add_argument("db")
    ap.add_argument("--json", action="store_true", help="machine-readable output")
    ap.add_argument("--examples", type=int, default=0, metavar="N",
                    help="show up to N vendor notations for each discordant consequence pair")
    ap.add_argument("--benchmark", action="store_true", help="time the use-case queries (SQL and CGI)")
    args = ap.parse_args(argv)
    r = report(args.db, examples=args.examples)
    if args.benchmark:
        r["benchmark"] = benchmark(args.db)
    if args.json:
        print(json.dumps(r, indent=2, ensure_ascii=False))
        return 0
    print("== Reports")
    for x in r["reports"]:
        print(f"  {x['panel_name']:<20} {str(x['panel_version']):<24} {x['genome_build']:<8} "
              f"({x['genome_build_source']})  n={x['n']}")
    pt = r["patients"]
    print(f"== Patients: {pt['distinct_patient_ids']} distinct patient IDs; {pt['reports_without_patient_id']} reports "
          f"without one; {pt['patients_with_reports_from_2_or_more_assays']} patients tested with 2 or more assays")
    print("== Variant rows")
    for x in r["variant_rows"]:
        print(f"  {x['panel_name']:<20} {x['variant_type']:<14} {x['n']}")
    print("== Short-variant consequence source")
    for x in r["consequence_source"]:
        print(f"  {x['panel_name']:<20} {x['source']:<13} {x['n']}")
    print("== Biomarkers")
    for x in r["biomarkers"]:
        print(f"  {x['name']:<4} {x['assay']:<75} n={x['n']} value={x['with_value']} call={x['with_call']}")
    c = r["curation"]
    print(f"== Curation: {c['reports']} reports; disease label {c['with_disease_label']}, OncoTree "
          f"{c['with_oncotree_code']}, patient id {c['with_patient_id']}, sidecar-curated {c['with_curated_fields']}")
    g = r["gene_lists"]
    print(f"== Gene lists: {g['reports_with_gene_list']} of {g['reports']} reports have an assay gene list")
    og = r["genes_outside_gene_list"]
    print(f"== Short-variant rows in genes outside the assay gene list: {sum(x['n'] for x in og)} "
          f"({len(og)} panel/gene pairs); gene symbols harmonized to HGNC: {r['gene_symbols_harmonized']} rows")
    for x in og[:25]:
        print(f"     {x['panel_name']:<20} {str(x['panel_version']):<22} {x['gene']:<12} n={x['n']}")
    i = r["inference_vs_vendor"]
    print(f"== Rule-based consequence vs vendor category: {i['variants_with_vendor_category']} variants; "
          f"same SO term {i['same_term']} ({i['same_term_pct']}%), same display group "
          f"{i['same_display_group']} ({i['same_group_pct']}%)")
    for d in i["discordant"]:
        print(f"     vendor {d['vendor']:<26} inferred {d['inferred']:<26} n={d['n']}")
        for e in d["examples"]:
            print(f"         e.g. category={e['vendor_category']!r} protein={e['protein_effect']!r} cds={e['cds_effect']!r}")
    gd = r["guardant"]
    print(f"== Guardant: {gd['reports']} reports; {gd['variant_rows_loaded']} variant rows loaded; "
          f"{gd['uncalled_candidate_rows_not_loaded']} uncalled candidate rows not loaded")
    if "benchmark" in r:
        b = r["benchmark"]
        print(f"== Benchmark: database {b['database_bytes'] / 1e6:.1f} MB; rows {b['rows']}")
        for k, v in b["sql_median_ms"].items():
            print(f"     SQL  {k:<60} {v} ms")
        for k, v in b.get("cgi_end_to_end_median_ms", {}).items():
            print(f"     CGI  {k:<60} {v} ms (includes Perl start-up)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
