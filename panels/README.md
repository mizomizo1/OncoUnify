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

`template.tsv.example` is an empty template.  Synthetic lists used by the
tests live in `tests/data/panels/` and must not be used for real data.
