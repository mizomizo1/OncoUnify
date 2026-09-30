#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Load panel gene lists (assay content) into OncoUnify.

Gene lists are the denominators for cross-panel statistics: a gene counts as
"tested" in a case only if the case's assay (cases.panel_name,
cases.panel_version) interrogates it.  Without them, mutation frequencies
pooled across panels are confounded by panel composition.

Usage:
    python3 load_panel_genes.py panels.db panels/  [more files or directories]

File format (tab-separated, UTF-8), one file per assay version:

    # panel_name: FoundationOne
    # panel_version: FoundationOneDx
    # description: FoundationOne CDx (324 genes)
    # source: FoundationOne CDx Technical Information, RAL-0003 v33.0
    gene	short_variants	copy_number	rearrangements
    ABL1	1	1	0
    ...

`panel_name` and `panel_version` must match the values the vendor loader
records in `cases`.  Re-loading a file replaces that version's gene list.
Exit status is 1 if any file is malformed.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import oncounify_core as oc

LOADER = "load_panel_genes.py"


def parse_panel_file(path: Path) -> Tuple[Dict[str, str], List[Tuple[str, int, int, int]]]:
    meta: Dict[str, str] = {}
    rows: List[Tuple[str, int, int, int]] = []
    header = None
    for n, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        if not line.strip():
            continue
        if line.startswith("#"):
            if ":" in line:
                k, v = line[1:].split(":", 1)
                meta[k.strip().lower()] = v.strip()
            continue
        cells = [c.strip() for c in line.split("\t")]
        if header is None:
            header = [c.lower() for c in cells]
            if not header or header[0] != "gene":
                raise ValueError(f"line {n}: the first column must be 'gene'")
            continue
        rec = dict(zip(header, cells))
        gene = oc.clean_str(rec.get("gene"))
        if not gene:
            raise ValueError(f"line {n}: empty gene symbol")

        def flag(col: str, default: int) -> int:
            v = rec.get(col)
            if v is None or v == "":
                return default
            if v not in ("0", "1"):
                raise ValueError(f"line {n}: column {col!r} must be 0 or 1, got {v!r}")
            return int(v)

        rows.append((gene, flag("short_variants", 1), flag("copy_number", 0), flag("rearrangements", 0)))
    for k in ("panel_name", "panel_version"):
        if not meta.get(k):
            raise ValueError(f"missing '# {k}: ...' header line")
    genes = [r[0] for r in rows]
    dups = sorted({g for g in genes if genes.count(g) > 1})
    if dups:
        raise ValueError(f"duplicated gene symbol(s): {', '.join(dups)}")
    if not rows:
        raise ValueError("no genes listed")
    # previous HGNC symbols -> current symbol, as for variants; merge flags if both appear
    merged: Dict[str, List[int]] = {}
    for gene, sv, cn, rr in rows:
        cur = oc.current_symbol(gene)
        f = merged.setdefault(cur, [0, 0, 0])
        merged[cur] = [max(f[0], sv), max(f[1], cn), max(f[2], rr)]
    return meta, [(g, *flags) for g, flags in merged.items()]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog=LOADER, description="Load assay gene lists (panel content).")
    ap.add_argument("db")
    ap.add_argument("inputs", nargs="+", help="TSV files or directories containing *.tsv")
    ap.add_argument("--allow-synthetic", action="store_true",
                    help="also load lists marked SYNTHETIC (tests/data/panels; tests and demos only)")
    args = ap.parse_args(argv)

    files = oc.iter_input_files(args.inputs, (".tsv",))
    if not files:
        print(f"[FATAL] {LOADER}: no .tsv files found", file=sys.stderr)
        return 2
    conn = oc.connect(args.db)
    oc.init_db(conn)
    failed = 0
    for path in files:
        try:
            meta, rows = parse_panel_file(path)
            synthetic = "SYNTHETIC" in (meta.get("description") or "").upper()
            if synthetic and not args.allow_synthetic:
                raise ValueError("SYNTHETIC test list, not the vendor's assay content: not loaded. The official "
                                 "lists are in panels/ (use --allow-synthetic only for tests and demos)")
            with conn:
                conn.execute(
                    "INSERT INTO panel_versions (panel_name, panel_version, description, source) VALUES (?, ?, ?, ?) "
                    "ON CONFLICT(panel_name, panel_version) DO UPDATE SET description = excluded.description, "
                    "source = excluded.source",
                    (meta["panel_name"], meta["panel_version"], meta.get("description"), meta.get("source")),
                )
                pv_id = conn.execute(
                    "SELECT panel_version_id FROM panel_versions WHERE panel_name = ? AND panel_version = ?",
                    (meta["panel_name"], meta["panel_version"]),
                ).fetchone()[0]
                conn.execute("DELETE FROM panel_genes WHERE panel_version_id = ?", (pv_id,))
                conn.executemany(
                    "INSERT INTO panel_genes (panel_version_id, gene, short_variants, copy_number, rearrangements) "
                    "VALUES (?, ?, ?, ?, ?)",
                    [(pv_id, *r) for r in rows],
                )
            print(f"[INFO] {meta['panel_name']} / {meta['panel_version']}: {len(rows)} genes <- {path}",
                  file=sys.stderr)
            if synthetic:
                print(f"[WARNING] {path.name} is a SYNTHETIC test list, not the vendor's assay content",
                      file=sys.stderr)
        except Exception as exc:
            failed += 1
            print(f"[ERROR] {path}: {exc}", file=sys.stderr)
    conn.close()
    print(f"[SUMMARY] {LOADER}: files={len(files)} failed={failed}", file=sys.stderr)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
