#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Registry quality report for an OncoUnify database.

    python3 tools/qc_report.py panels.db [--json]

Reports, per assay: number of reports, genome builds and their source,
variant rows by class, how short-variant consequences were assigned (vendor /
inferred / unclassified), biomarker assays, curation coverage, gene-list
coverage — and a reliability estimate for the rule-based consequence
inference: for every short variant whose consequence came from the vendor
(Foundation Medicine supplies one), the term is recomputed from the protein
and coding notation alone, as it would be for a vendor that supplies none,
and the two are compared.
"""

from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import oncounify_core as oc  # noqa: E402


def q(conn, sql, *args):
    return conn.execute(sql, args).fetchall()


def report(db: str) -> dict:
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

    # reliability of rule-based inference against vendor-supplied categories
    group = {t: g for t, g in q(conn, "SELECT term, display_group FROM so_terms")}
    rows = q(conn, "SELECT v.hgvs_p, v.cds_effect, v.functional_effect FROM variants v "
                   "WHERE v.variant_type = 'short_variant' AND v.functional_effect_source = 'vendor'")
    pairs = Counter()
    for hgvs_p, cds, vendor_term in rows:
        inferred = oc.classify_consequence(hgvs_p, cds)[0]
        pairs[(vendor_term, inferred)] += 1
    n = sum(pairs.values())
    term_agree = sum(k for (a, b), k in pairs.items() if a == b)
    group_agree = sum(k for (a, b), k in pairs.items() if group.get(a) == group.get(b))
    out["inference_vs_vendor"] = {
        "variants_with_vendor_category": n,
        "same_term": term_agree,
        "same_display_group": group_agree,
        "same_term_pct": round(100 * term_agree / n, 1) if n else None,
        "same_group_pct": round(100 * group_agree / n, 1) if n else None,
        "discordant": [{"vendor": a, "inferred": b, "n": k} for (a, b), k in sorted(pairs.items(), key=lambda x: -x[1])
                       if a != b],
    }
    conn.close()
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Quality report for an OncoUnify database.")
    ap.add_argument("db")
    ap.add_argument("--json", action="store_true", help="machine-readable output")
    args = ap.parse_args(argv)
    r = report(args.db)
    if args.json:
        print(json.dumps(r, indent=2, ensure_ascii=False))
        return 0
    print("== Reports")
    for x in r["reports"]:
        print(f"  {x['panel_name']:<20} {str(x['panel_version']):<24} {x['genome_build']:<8} "
              f"({x['genome_build_source']})  n={x['n']}")
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
    i = r["inference_vs_vendor"]
    print(f"== Rule-based consequence vs vendor category: {i['variants_with_vendor_category']} variants; "
          f"same SO term {i['same_term']} ({i['same_term_pct']}%), same display group "
          f"{i['same_display_group']} ({i['same_group_pct']}%)")
    for d in i["discordant"]:
        print(f"     vendor {d['vendor']:<26} inferred {d['inferred']:<26} n={d['n']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
