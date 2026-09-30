# Assay gene lists (panel content)

Cross-panel statistics need to know which genes each assay interrogates:
a gene is "tested" in a report only if the report's assay covers it.
Place one tab-separated file per assay version in this directory and load
them with:

```bash
python3 load_panel_genes.py panels.db panels/
```

## File format

```
# panel_name: FoundationOne
# panel_version: FoundationOneDx
# description: FoundationOne CDx
# source: <document, version and date the list was taken from>
gene	short_variants	copy_number	rearrangements
ABL1	1	0	0
...
```

* `panel_name` and `panel_version` must equal the values the loader writes
  to `cases` (check with
  `SELECT DISTINCT panel_name, panel_version FROM cases;`).
  Defaults: `FoundationOne` / `FoundationOneDx`,
  `FoundationOneLiquid` / `FoundationOneLiquidDx`, `GenMineTOP` / `TDv6+TRv6`
  (from the report), `Guardant` / `Guardant360 CDx`.
* The flags state whether the assay reports short variants, copy-number
  changes and rearrangements in that gene.
* Record the provenance in `# source:`.  Take the list from the regulatory
  labeling of the assay version used at your institution (e.g. the Japanese
  package insert or the FDA Technical Information), because gene content
  changes between versions.

## Lists shipped with OncoUnify

| file | assay | genes (short variants / CN / rearrangements) | source |
|---|---|---|---|
| FoundationOne_FoundationOneDx.tsv | FoundationOne CDx | 324 (311 / 309 / 36) | FMI Technical Information RAL-0003 v34.0, Tables 2–3 |
| FoundationOneLiquid_FoundationOneLiquidDx.tsv | FoundationOne Liquid CDx | 324 (311 / 3 / 6) | FMI Technical Information RAL-0035-19, Tables 2–3 |
| GenMineTOP_TDv6+TRv6.tsv, GenMineTOP_TDv1.1.0+TRv6.4.3.tsv | GenMineTOP | 1110 (737 / 737 / 457) | GenMine Labs brochure Ver.4.0 (DNA 737, fusion 455, exon skipping 5) |
| Guardant_Guardant360CDx.tsv | Guardant360 CDx (Japan) | 74 (74 / 18 / 6) | Guardant Health Japan product information |

They were extracted on 2026-09-30 from public documents (URLs in each file's
`# source:` line); every symbol was checked against the HGNC complete set.
Copy-number and rearrangement flags follow what the document states is
reported; only the short-variant flag is used for frequencies.  Both
GenMineTOP version labels are assumed to have the approved content of the
brochure.  `GenMineTOP_expression_genes.txt` lists the 27 expression genes
for reference.  Verify the lists against the labeling of the assay versions
used at your institution; `tools/qc_report.py` lists reported short variants
that fall outside the declared content.

`template.tsv.example` is an empty template.  Synthetic lists used by the
tests live in `tests/data/panels/` and must not be used for real data.
