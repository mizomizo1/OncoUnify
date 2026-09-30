#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Write gene_symbol_map.tsv — previous HGNC symbols of the genes in the
supported assays and their current approved symbol.  The loaders replace
these symbols in reports and gene lists (keeping the vendor's symbol), and
the search interface resolves them.

    python3 tools/make_symbol_map.py hgnc_complete_set.txt [--date YYYY-MM-DD]

hgnc_complete_set.txt is the HGNC complete set (tab-separated), from
https://www.genenames.org/download/archive/ .  A previous symbol is used
only if it is unambiguous: it belongs to exactly one approved gene, is not
the approved symbol of another gene and is not an alias of another gene
(so that, for example, a search for ERK is not answered with EPHB2).  The
exception is a synonym that the vendors' own documentation uses for the
gene (VENDOR_SYNONYMS).  Aliases are not used as a source of mappings.
Genes are those of every panels/*.tsv file and panels/*_expression_genes.txt.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "gene_symbol_map.tsv"

# previous symbols that are also an alias of another gene, but that the vendor
# gene lists give for this gene: "KMT2D (MLL2)" in the FoundationOne CDx and
# FoundationOne Liquid CDx technical information and the GenMineTOP brochure
VENDOR_SYNONYMS = {"MLL2": "KMT2D"}


def panel_symbols() -> set:
    genes = set()
    for f in sorted((ROOT / "panels").glob("*.tsv")) + sorted((ROOT / "panels").glob("*_expression_genes.txt")):
        for line in f.read_text(encoding="utf-8").splitlines():
            g = line.split("\t")[0].strip()
            if g and not g.startswith("#") and g != "gene":
                genes.add(g)
    return genes


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("hgnc_complete_set")
    ap.add_argument("--date", default=datetime.date.today().isoformat(), help="download date of the HGNC file")
    args = ap.parse_args(argv)

    approved, prev_of, alias_of = set(), defaultdict(set), defaultdict(set)
    with open(args.hgnc_complete_set, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("status", "Approved") != "Approved":
                continue
            approved.add(r["symbol"])
            for p in filter(None, (r.get("prev_symbol") or "").split("|")):
                prev_of[p].add(r["symbol"])
            for a in filter(None, (r.get("alias_symbol") or "").split("|")):
                alias_of[a].add(r["symbol"])
    unique_prev = {p: next(iter(s)) for p, s in prev_of.items() if len(s) == 1 and p not in approved}

    listed = panel_symbols()
    unknown = sorted(g for g in listed if g not in approved and g not in unique_prev)
    if unknown:
        print(f"[ERROR] not an approved or unambiguous previous HGNC symbol: {', '.join(unknown)}", file=sys.stderr)
        return 1
    current = {unique_prev.get(g, g) for g in listed}
    pairs, excluded = [], []
    for p, s in sorted(unique_prev.items()):
        if s not in current:
            continue
        if alias_of[p] - {s} and VENDOR_SYNONYMS.get(p) != s:
            excluded.append(p)
            continue
        pairs.append((p, s))
    missing = sorted(g for g in listed if g not in approved and g not in dict(pairs))
    if missing:
        print(f"[ERROR] gene-list symbols that would not be harmonized: {', '.join(missing)}; "
              "add them to VENDOR_SYNONYMS if the vendor documents the synonym", file=sys.stderr)
        return 1

    lines = [
        "# Previous HGNC symbols of the genes in the supported assays -> current approved symbol.",
        f"# source: HGNC complete set, downloaded {args.date} (https://www.genenames.org/download/archive/)",
        "# rule: previous symbols of exactly one gene that are neither an approved symbol nor an alias of another "
        "gene (exception: synonyms in the vendor gene lists, " + ", ".join(f"{p}->{s}" for p, s in VENDOR_SYNONYMS.items())
        + "); aliases are not used",
        f"# genes: {len(current)} current symbols from panels/; regenerate with tools/make_symbol_map.py",
        "previous_symbol\tsymbol",
    ] + [f"{p}\t{s}" for p, s in pairs]
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {OUT}: {len(pairs)} previous symbols for {len({s for _, s in pairs})} genes "
          f"({len(excluded)} excluded as aliases of another gene)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
