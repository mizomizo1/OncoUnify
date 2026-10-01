# OncoUnify

[![tests](https://github.com/mizomizo1/OncoUnify/actions/workflows/tests.yml/badge.svg)](https://github.com/mizomizo1/OncoUnify/actions/workflows/tests.yml)

OncoUnify turns the heterogeneous deliverables of commercial cancer genomic
panels into one normalized, queryable SQLite database that stays inside the
institution.  It is bioinformatics infrastructure: installed once (by IT or a
bioinformatician), run on a schedule, and queried by anyone with SQL, R,
Python — or through the bundled reference web interface.

| Vendor deliverable | Loader | Genome build |
|---|---|---|
| FoundationOne CDx / FoundationOne Liquid CDx (Foundation Medicine, XML) | `load_foundation.py` | GRCh37 (loader default) |
| GenMineTOP Cancer Genome Profiling System (GenMine Labs, XML) | `load_genminetop.py` | read from the report (GRCh38) |
| Guardant360 CDx (Guardant Health, Excel `.xlsx`) | `load_guardant.py` | GRCh37 (loader default) |

## What is normalized

* **One report = one `cases` row** keyed by `(panel_name, report_id)`, with
  the reference assembly (`genome_build`), assay version, provenance (source
  file, SHA-256, loader version) and any vendor field without a canonical
  home kept as JSON.
* **Variants** of every class (SNV/indel, CNV, fusion/rearrangement/exon
  skipping, expression) in one table with typed views.  Vendor strings are
  stored verbatim; canonical keys sit next to them: `hgvs_p` (one-letter HGVS
  protein key, e.g. `p.G12D` for `G12D`, `p.G12D`, `p.Gly12Asp`) and `hgvs_c`.
* **Functional consequence** as a Sequence Ontology term (`missense_variant`,
  `frameshift_variant`, `splice_acceptor_variant`, `inframe_deletion`, …),
  assigned by one shared rule set in every loader, with the SO accession and
  whether it came from the vendor or was inferred.
* **Biomarkers** (TMB, MSI) in their own table with value, unit, controlled
  call, verbatim call and assay — values of different assays are never pooled.
* **Assay content** (gene lists) as data, so that cross-panel mutation
  frequencies use the right denominator (`n_mutated / n_tested`).
* **Curation**: a sidecar CSV supplies what vendor files lack (disease label,
  OncoTree code, patient identifier linking reports of one patient).

See [docs/SCHEMA.md](docs/SCHEMA.md) for the data dictionary and
[docs/field_mapping.tsv](docs/field_mapping.tsv) for the complete
vendor-to-canonical field mapping.

## Quick start (Docker demo with synthetic data)

```bash
git clone https://github.com/mizomizo1/OncoUnify.git
cd OncoUnify
docker compose up --build
# open http://localhost:8080/  (demo database built from tests/data)
```

## Quick start (command line)

Python ≥ 3.9 with `pandas` and `openpyxl` (Guardant only):

```bash
python3 load_foundation.py  panels.db /reports/foundation  --case-metadata curation.csv
python3 load_genminetop.py  panels.db /reports/genminetop  --case-metadata curation.csv
python3 load_guardant.py    panels.db /reports/guardant    --case-metadata curation.csv
python3 load_panel_genes.py panels.db panels/
sqlite3 panels.db "SELECT panel_name, report_id, gene, hgvs_p, functional_effect
                   FROM v_short_variants WHERE gene = 'KRAS' AND hgvs_p = 'p.G12D';"
```

Every loader replaces a report and its child rows in one transaction (re-runs
never duplicate data), prints a per-run summary, and exits with status 1 if
any input file failed — suitable for a nightly `cron` job.  To deploy the web
interface on Apache see [docs/INSTALL.md](docs/INSTALL.md).  To upgrade a
database created by OncoUnify 1.x run `python3 migrate_db.py panels.db`.

## Repository layout

```
OncoUnify/
├── schema.sql              canonical schema (version 2)
├── oncounify_core.py       shared normalization library used by every loader
├── gene_symbol_map.tsv     previous HGNC symbols -> current symbol (tools/make_symbol_map.py)
├── load_foundation.py      FoundationOne CDx / FoundationOne Liquid CDx XML
├── load_genminetop.py      GenMineTOP XML
├── load_guardant.py        Guardant360 CDx XLSX
├── load_panel_genes.py     assay gene lists (denominators)
├── migrate_db.py           schema v1 -> v2 upgrade
├── panel_search.cgi        reference web interface: cross-vendor search, TSV export
├── case_detail.cgi         per-report drill-down
├── panel_stats.cgi         registry overview with panel-aware frequencies
├── suggest.cgi             autocomplete endpoint (JSON)
├── search.html             search form
├── panels/                 official assay gene lists (see panels/README.md)
├── tools/qc_report.py      registry quality report (incl. consequence-inference check)
├── tools/make_symbol_map.py  regenerates gene_symbol_map.tsv from the HGNC complete set
├── tools/scale_benchmark.py  loads and queries synthetic registries of 10^2-10^4 reports
├── docs/                   INSTALL, SCHEMA, LOADER_DEV, CURATION, field_mapping.tsv, WORKED_EXAMPLE, BENCHMARK
├── docker/, Dockerfile, docker-compose.yml
└── tests/                  synthetic fixtures, fixture generator, unit and end-to-end tests
```

## Tests

```bash
python3 -m unittest discover -s tests -v
```

The fixtures under `tests/data/` are fully synthetic (invented identifiers,
public hotspot variants).  **Never place real vendor reports inside this
repository.**

## Requesting support for a new assay

You do not need to write code to have another assay supported.  Open an
issue with the
[**Request support for a new assay**](https://github.com/mizomizo1/OncoUnify/issues/new?template=new_assay.yml)
template and provide one of the following:

* a de-identified example report (every identifier, name, date and
  free-text clinical field replaced with invented values);
* a synthetic report in the vendor's format; or
* the vendor's format specification (for example an XSD or a data
  dictionary).

We write the loader, derive synthetic test fixtures from the example and add
them to the test suite, so that the new assay is tested like the existing
ones.  **Never attach real patient reports: issues are public.**  If a
de-identified file cannot be posted publicly, say so in the issue and we will
arrange another route.

## Extending OncoUnify

Contributed loaders are equally welcome.  A new assay needs one
`load_<vendor>.py` that parses the vendor file into the canonical
dictionaries and hands them to `oncounify_core`; the schema, the
normalization rules and the web layer are shared.  A contributed loader
comes with synthetic fixtures and tests.  See
[docs/LOADER_DEV.md](docs/LOADER_DEV.md).

## License and citation

MIT License (see [LICENSE](LICENSE)).  If you use OncoUnify, please cite:
Mizoue H, Higasa K. OncoUnify: an open-source ingestion and normalization
layer for multi-vendor cancer genomic panel reports (manuscript in
preparation).
