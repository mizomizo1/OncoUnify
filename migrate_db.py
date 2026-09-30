#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Upgrade an OncoUnify database from schema version 1 to version 2.

Usage:
    python3 migrate_db.py panels.db [--no-backup]

The preferred upgrade path is to rebuild the database by re-running the
loaders on the original vendor files, because several version-2
improvements need the source files (Guardant call filtering, GenMineTOP
fusion breakpoints and panel version, vendor fields kept in `extra`).  This
tool exists for sites whose source files are no longer at hand.

What it does:
  * builds a new version-2 database next to the old one and swaps it in only
    after every row has been copied (the old file is kept as a backup);
  * preserves case_id values;
  * removes the duplicate variant and non-human rows that version-1 loaders
    appended when they were re-run;
  * records genome_build from the vendor (Foundation Medicine and Guardant:
    GRCh37; GenMineTOP: GRCh38) with genome_build_source = 'migration';
  * moves TMB/MSI into the `biomarkers` table with unit and assay, and moves
    the Guardant MSI score that version 1 stored in cases.tmb_score to an MSI
    record;
  * recomputes hgvs_p / hgvs_c and the Sequence Ontology functional_effect
    with the shared version-2 rules;
  * reclassifies GenMineTOP RNA exon-skipping rows as rearrangements;
  * replaces previous HGNC gene symbols by current ones (vendor symbol kept in extra).
"""

from __future__ import annotations

import argparse
import datetime as _dt
import os
import shutil
import sqlite3
import sys
from typing import Any, Dict, List

import oncounify_core as oc

PROG = "migrate_db.py"


def vendor_profile(panel_name: str, vendor: str) -> Dict[str, Any]:
    p = f"{panel_name or ''} {vendor or ''}".lower()
    if "foundation" in p:
        return {"vendor": "Foundation Medicine", "build": "GRCh37", "family": "foundation"}
    if "genmine" in p or "todai" in p:
        return {"vendor": "GenMine Labs", "build": "GRCh38", "family": "genminetop"}
    if "guardant" in p:
        return {"vendor": "Guardant Health", "build": "GRCh37", "family": "guardant"}
    return {"vendor": vendor or None, "build": "unknown", "family": "other"}


def migrate(db_path: str, backup: bool = True) -> Dict[str, int]:
    old = sqlite3.connect(db_path)
    old.row_factory = sqlite3.Row
    version = old.execute("PRAGMA user_version").fetchone()[0]
    tables = {r[0] for r in old.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    if version >= oc.SCHEMA_VERSION:
        print(f"[INFO] {db_path} is already at schema version {version}; nothing to do.", file=sys.stderr)
        return {}
    if not {"cases", "variants"} <= tables:
        raise SystemExit(f"[FATAL] {db_path} does not look like an OncoUnify database")

    stamp = _dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    if backup:
        bak = f"{db_path}.v1-backup-{stamp}"
        shutil.copy2(db_path, bak)
        print(f"[INFO] backup written to {bak}", file=sys.stderr)

    tmp = f"{db_path}.migrating"
    if os.path.exists(tmp):
        os.remove(tmp)
    new = oc.connect(tmp)
    oc.init_db(new)

    now = _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat()
    stats = {"cases": 0, "variants_in": 0, "variants_out": 0, "duplicates_removed": 0, "skipped_rows": 0,
             "biomarkers": 0, "guardant_msi_moved": 0, "non_human_in": 0, "non_human_out": 0}
    families = set()
    case_cols = {r[1] for r in old.execute("PRAGMA table_info(cases)")}

    def g(row, key):
        return row[key] if key in row.keys() else None

    with new:
        for c in old.execute("SELECT * FROM cases ORDER BY case_id"):
            prof = vendor_profile(c["panel_name"], g(c, "vendor"))
            families.add(prof["family"])
            legacy = {"legacy_other_info": g(c, "other_info")}
            for k in ("ccat_date", "ccat_cancer"):
                if k in case_cols:
                    legacy[k] = g(c, k)
            panel_version = None
            if prof["family"] == "foundation":
                panel_version = oc.clean_str(g(c, "test_type")) or oc.clean_str(g(c, "panel_type"))
            elif prof["family"] == "guardant":
                panel_version = "Guardant360 CDx"
            case = {
                "panel_name": c["panel_name"],
                "panel_version": panel_version,
                "panel_type": g(c, "panel_type"),
                "test_type": g(c, "test_type"),
                "vendor": prof["vendor"],
                "report_id": oc.clean_str(g(c, "report_id")) or f"legacy-case-{c['case_id']}",
                "patient_id": g(c, "patient_id"),
                "sex": g(c, "sex"),
                "age": g(c, "age"),
                "date": g(c, "date"),
                "genome_build": prof["build"],
                "genome_build_source": "migration",
                "disease": g(c, "disease"),
                "disease_ontology": g(c, "disease_ontology"),
                "tissue_of_origin": g(c, "tissue_of_origin"),
                "pathology_diagnosis": g(c, "pathology_diagnosis"),
                "specimen_id": g(c, "specimen_id"),
                "percent_tumor_nuclei": g(c, "percent_tumor_nuclei"),
                "purity": g(c, "purity"),
                "non_human_content": g(c, "non_human_content"),
                "other_info": oc.to_json(legacy),
                "format_version": "migrated from schema v1",
                "loader": f"{PROG} {oc.__version__}",
                "loaded_at": now,
            }
            cols = ["case_id"] + list(case)
            new.execute(f"INSERT INTO cases ({', '.join(cols)}) VALUES ({', '.join('?' * len(cols))})",
                        [c["case_id"]] + list(case.values()))
            stats["cases"] += 1

            # --- biomarkers (T1: Guardant's tmb_score held the MSI score) ---
            bms: List[Dict[str, Any]] = []
            msi_raw, tmb = g(c, "msi_status"), g(c, "tmb_score")
            if prof["family"] == "guardant":
                if msi_raw is not None or tmb is not None:
                    bms.append({"name": "MSI", "value": tmb, "unit": "Guardant MSI score",
                                "call": oc.normalize_msi_call(msi_raw), "call_raw": msi_raw,
                                "assay": "Guardant360 CDx MSI (plasma cfDNA)",
                                "source_field": "migrated: v1 cases.msi_status, cases.tmb_score (MSI score)"})
                    if tmb is not None:
                        stats["guardant_msi_moved"] += 1
            else:
                label = {"foundation": "FoundationOne Liquid CDx" if "liquid" in (c["panel_name"] or "").lower()
                         else "FoundationOne CDx", "genminetop": "GenMineTOP"}.get(prof["family"], c["panel_name"])
                if msi_raw is not None:
                    bms.append({"name": "MSI", "call": oc.normalize_msi_call(msi_raw), "call_raw": msi_raw,
                                "assay": f"{label} MSI", "source_field": "migrated: v1 cases.msi_status"})
                if tmb is not None or g(c, "tmb_status") is not None:
                    assay = ("GenMineTOP exonic non-synonymous alteration frequency (tumor-normal paired)"
                             if prof["family"] == "genminetop" else
                             f"{label} {'blood TMB (bTMB)' if 'Liquid' in label else 'tissue TMB'}")
                    bms.append({"name": "TMB", "value": tmb, "unit": oc.normalize_tmb_unit(g(c, "tmb_unit")),
                                "call": oc.normalize_tmb_call(g(c, "tmb_status")),
                                "call_raw": g(c, "tmb_status"), "assay": assay,
                                "source_field": "migrated: v1 cases.tmb_score"})
            for b in bms:
                new.execute(f"INSERT INTO biomarkers (case_id, {', '.join(oc.BIOMARKER_COLUMNS)}) "
                            f"VALUES (?, {', '.join('?' * len(oc.BIOMARKER_COLUMNS))})",
                            [c["case_id"]] + [b.get(k) for k in oc.BIOMARKER_COLUMNS])
                stats["biomarkers"] += 1

            # --- variants (T2: drop exact duplicates) ---
            seen = set()
            for v in old.execute("SELECT * FROM variants WHERE case_id = ? ORDER BY variant_id", (c["case_id"],)):
                stats["variants_in"] += 1
                key = tuple(v[k] for k in v.keys() if k != "variant_id")
                if key in seen:
                    stats["duplicates_removed"] += 1
                    continue
                seen.add(key)
                vt = g(v, "variant_type")
                sub = g(v, "variant_subtype")
                row: Dict[str, Any] = {k: g(v, k) for k in oc.VARIANT_COLUMNS if k in v.keys()}
                row["chrom"] = oc.normalize_chrom(g(v, "chrom"))
                row["origin"] = oc.normalize_origin(g(v, "origin"))
                row["cnv_type"] = oc.normalize_cnv_type(g(v, "cnv_type") or (sub if vt == "cnv" else None))
                row["in_frame"] = oc.normalize_in_frame(g(v, "in_frame"))
                for k in ("functional_effect", "functional_effect_so", "functional_effect_raw",
                          "functional_effect_source", "hgvs_p", "hgvs_c"):
                    row[k] = None
                if (sub or "").lower() == "splicing-variant":
                    vt = "rearrangement"
                    row.update({"functional_effect": "exon_loss_variant",
                                "functional_effect_so": oc.SO_TERMS["exon_loss_variant"],
                                "functional_effect_raw": sub, "functional_effect_source": "vendor"})
                if vt not in oc.VARIANT_TYPES:
                    stats["skipped_rows"] += 1
                    continue
                row["variant_type"] = vt
                if vt == "short_variant":
                    old_fe = g(v, "functional_effect")
                    vendor_term = old_fe if prof["family"] == "foundation" else None
                    sv = oc.short_variant(gene=g(v, "gene"), protein=g(v, "protein_effect"), cds=g(v, "cds_effect"),
                                          vendor_term=vendor_term, hint=sub, ref=g(v, "ref"), alt=g(v, "alt"))
                    for k in ("hgvs_p", "hgvs_c", "functional_effect", "functional_effect_so",
                              "functional_effect_source"):
                        row[k] = sv[k]
                    row["functional_effect_raw"] = oc.clean_str(old_fe) or oc.clean_str(sub)
                oc.harmonize_gene_symbols([row])
                cols = ["case_id"] + list(row)
                new.execute(f"INSERT INTO variants ({', '.join(cols)}) VALUES ({', '.join('?' * len(cols))})",
                            [c["case_id"]] + list(row.values()))
                stats["variants_out"] += 1

            if "non_human_contents" in tables:
                seen_nh = set()
                for n in old.execute("SELECT organism, reads_per_million, status, sample FROM non_human_contents "
                                     "WHERE case_id = ? ORDER BY id", (c["case_id"],)):
                    stats["non_human_in"] += 1
                    k = tuple(n)
                    if k in seen_nh:
                        continue
                    seen_nh.add(k)
                    new.execute("INSERT INTO non_human_contents (case_id, organism, reads_per_million, status, sample) "
                                "VALUES (?, ?, ?, ?, ?)", (c["case_id"], *k))
                    stats["non_human_out"] += 1
        # keep AUTOINCREMENT counters beyond the preserved ids
        new.execute("UPDATE sqlite_sequence SET seq = (SELECT MAX(case_id) FROM cases) WHERE name = 'cases'")

    new.close()
    old.close()
    os.replace(tmp, db_path)

    print(f"[SUMMARY] {PROG}: cases={stats['cases']} variants {stats['variants_in']} -> {stats['variants_out']} "
          f"(duplicates removed={stats['duplicates_removed']}, unsupported rows skipped={stats['skipped_rows']}) "
          f"biomarkers={stats['biomarkers']} (Guardant MSI scores moved out of tmb_score={stats['guardant_msi_moved']}) "
          f"non_human {stats['non_human_in']} -> {stats['non_human_out']}", file=sys.stderr)
    print("[WARN] genome builds were assigned from vendor defaults (genome_build_source='migration').", file=sys.stderr)
    if "guardant" in families:
        print("[WARN] Guardant reports were loaded by version 1 without call filtering; re-load the Guardant "
              "workbooks with load_guardant.py to drop candidates that were not reported (call = 0).",
              file=sys.stderr)
    if "genminetop" in families:
        print("[WARN] GenMineTOP panel_version, fusion breakpoints and vendor extras need the source XML; "
              "re-load with load_genminetop.py to populate them.", file=sys.stderr)
    return stats


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog=PROG, description="Upgrade an OncoUnify database to schema version 2.")
    ap.add_argument("db")
    ap.add_argument("--no-backup", action="store_true", help="do not keep a copy of the version-1 file")
    args = ap.parse_args(argv)
    if not os.path.isfile(args.db):
        print(f"[FATAL] {args.db} not found", file=sys.stderr)
        return 2
    migrate(args.db, backup=not args.no_backup)
    return 0


if __name__ == "__main__":
    sys.exit(main())
