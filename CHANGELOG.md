# Changelog

## 2.0.0 (2026-09)

Schema version 2; databases created by 1.x must be rebuilt from the source
reports or upgraded with `migrate_db.py`.

### Data correctness
- Re-running a loader replaces a report and its child rows in one
  transaction; 1.x appended a duplicate copy of every variant on each run.
- Guardant: the MSI score is stored as an MSI biomarker; 1.x wrote it to
  `cases.tmb_score`.
- Guardant: only reported alterations (`call = 1`) are loaded by default;
  1.x loaded every candidate row of the workbook, including germline
  polymorphisms.  `--include-uncalled` keeps them as `not_called`.
- Guardant: fusions are read from both `Fusion` and `Fusions` sheets;
  `percentage` populates `allele_fraction`; chromosomes are `chr`-prefixed.
- One Sequence Ontology consequence classifier for all vendors; 1.x used two
  inconsistent regular expressions (Guardant frameshifts were stored as
  missense) plus raw vendor values.
- Foundation Medicine `nonframeshift` is an in-frame change.
- GenMineTOP: RNA exon skipping is a rearrangement row (`exon_loss_variant`)
  rather than a short variant; fusion breakpoints are stored in pos/pos2.

### Schema
- `cases.genome_build` (NOT NULL) and `genome_build_source`; `panel_version`;
  `oncotree_code`; `curated_fields`; provenance columns.
- `variants.hgvs_p` / `hgvs_c` canonical keys next to the verbatim vendor
  strings; SO accession, raw value and source of `functional_effect`;
  `chrom2`; controlled vocabularies for variant type, origin, CNV type and
  fusion frame; JSON overflow columns.
- New tables `biomarkers`, `so_terms`, `panel_versions`, `panel_genes`;
  typed views; `ON DELETE CASCADE`; index on `variants(case_id)`.
- C-CAT placeholder columns removed.

### Loaders and tools
- Shared library `oncounify_core.py`; input format validation; per-run
  summary; exit status 1 when any file fails; `--fail-fast`,
  `--genome-build`, `--case-metadata`; Guardant `--filename-pattern`
  (no silent fallback), `--panel-version`, `--include-uncalled`.
- `load_panel_genes.py` (assay content), `migrate_db.py` (1.x upgrade).
- `tools/qc_report.py` (registry quality report and benchmark of the use-case
  queries) and `tools/scale_benchmark.py` (synthetic registries of 10^2-10^4
  reports loaded and queried; results in `docs/BENCHMARK.md`).

### Web interface
- Configuration through environment variables (`ONCOUNIFY_DB`, …).
- Exact gene matching including fusion partners; canonical protein matching;
  consequence filter; per-report view; LIKE and TSV escaping; HTML escaping
  on every page; genome build shown and exported.
- Panel-aware statistics (`n_mutated / n_tested`).
- `logout.cgi` removed.

### Gene content and symbols
- Official gene lists for FoundationOne CDx, FoundationOne Liquid CDx,
  GenMineTOP and Guardant360 CDx (Japan) extracted from public vendor
  documents and checked against HGNC.
- Previous HGNC symbols are harmonized to current symbols
  (`gene_symbol_map`, 599 unambiguous previous symbols of the genes in the
  supported assays, generated from the HGNC complete set by
  `tools/make_symbol_map.py`), in variants and gene lists; searches resolve
  them.

### Packaging
- MIT LICENSE file, Dockerfile and docker-compose demo, GitHub Actions CI,
  synthetic fixtures and automated tests, complete field mapping
  (`docs/field_mapping.tsv`), worked example, curation guide.
- Requests for new assays: README section and GitHub issue form
  (`.github/ISSUE_TEMPLATE/new_assay.yml`); the maintainers write the loader
  from a de-identified or synthetic example or the vendor's format
  specification, and add synthetic fixtures to the test suite.
