# OncoUnify — Canonical schema reference

OncoUnify stores all ingested data in three relational tables backed by a
single SQLite file (`panels.db`). The schema is intentionally narrow: no
vendor-specific table exists, and every loader maps its native fields onto
the same set of canonical columns.

The authoritative definition is [`schema.sql`](../schema.sql); this document
is its annotated companion.

---

## 1. `cases` — one row per genomic report

Each case represents a single completed panel test. The `(panel_name,
report_id)` pair is unique.

| Column | Type | Notes |
|---|---|---|
| `case_id` | INTEGER PK | autoincrement |
| `panel_name` | TEXT NOT NULL | stable vendor + assay identifier, e.g. `FoundationOne`, `FoundationOneLiquid`, `GenMineTOP`, `Guardant360` |
| `panel_type` | TEXT | raw `test-type` attribute as recorded by the vendor |
| `vendor` | TEXT | manufacturer name |
| `report_id` | TEXT | unique within `panel_name` |
| `patient_id` | TEXT | institutional MRN or pseudonym |
| `sex` | TEXT | as reported by the vendor |
| `age` | INTEGER | at the time of testing |
| `disease` | TEXT | free-text disease label from the vendor report |
| `disease_ontology` | TEXT | ontology identifier (e.g., NCIt term) when supplied |
| `tissue_of_origin` | TEXT | as reported by the vendor |
| `pathology_diagnosis` | TEXT | as reported by the vendor |
| `specimen_id` | TEXT | vendor-side specimen identifier |
| `test_type` | TEXT | duplicates `panel_type`; preserved for backwards compatibility |
| `percent_tumor_nuclei` | REAL | per-specimen QC |
| `purity` | REAL | tumour purity (0–1 or %, vendor-dependent) |
| `msi_status` | TEXT | microsatellite-instability status |
| `tmb_score` | REAL | tumour mutational burden score |
| `tmb_status` | TEXT | TMB call (low / intermediate / high / ...) |
| `tmb_unit` | TEXT | unit of `tmb_score` (e.g., mutations/Mb) |
| `date` | TEXT | report date in `YYYY-MM-DD` |
| `ccat_date` | TEXT | C-CAT specimen collection date |
| `ccat_cancer` | TEXT | C-CAT cancer-type label |
| `non_human_content` | REAL | case-level scalar from FoundationOne(Liquid) |
| `other_info` | TEXT | free-form `key=value;...` for vendor fields without a canonical home |

`UNIQUE(panel_name, report_id)` enforces loader idempotency.

---

## 2. `variants` — one row per reported molecular finding

A variant is anything the vendor reports as a finding: a short variant
(SNV/indel), a copy-number alteration, a rearrangement/fusion, an
expression call, or a biomarker. The `variant_type` column distinguishes
them.

| Column | Type | Notes |
|---|---|---|
| `variant_id` | INTEGER PK | autoincrement |
| `case_id` | INTEGER NOT NULL | FK → `cases(case_id)` |
| `gene` | TEXT | HGNC symbol when possible |
| `variant_type` | TEXT | `short_variant` / `cnv` / `rearrangement` / `expression` / `biomarker` |
| `variant_subtype` | TEXT | finer label (e.g., `amplification`, `deletion`, `fusion`) |
| `chrom` | TEXT | UCSC-style (`chr1`...) |
| `pos` | INTEGER | 1-based; first breakpoint for `rearrangement` |
| `pos2` | INTEGER | second breakpoint for `rearrangement` |
| `ref` / `alt` | TEXT | reference / alternate alleles |
| `cds_effect` | TEXT | HGVS coding-level expression (`c.35G>A`) |
| `protein_effect` | TEXT | HGVS protein-level expression with leading `p.` stripped (`G12D`) |
| `strand` | TEXT | `+` or `-` |
| `transcript` | TEXT | NM_ / ENST / RefSeq accession |
| `functional_effect` | TEXT | one of `missense` / `nonsense` / `frameshift` / `splice` / `synonymous` (loader-inferred when not vendor-supplied) |
| `effect` | TEXT | vendor-supplied effect description (free text) |
| `status` | TEXT | vendor call (known / likely / unknown) |
| `origin` | TEXT | germline / somatic when reported |
| `classification` | TEXT | vendor / curator interpretation |
| `allele_fraction` | REAL | 0–1 |
| `depth` | INTEGER | read depth |
| `copy_number` | REAL | for `cnv` |
| `cnv_ratio` | REAL | log2 or linear ratio (vendor-dependent) |
| `cnv_type` | TEXT | `amplification` / `loss` / ... |
| `other_gene` | TEXT | partner gene for `rearrangement` |
| `in_frame` | TEXT | `yes` / `no` for fusions |
| `supporting_read_pairs` | INTEGER | fusion evidence |
| `tpm` | REAL | for `expression` |
| `read_count` | INTEGER | for `expression` |
| `sample_name` | TEXT | vendor sample identifier when needed |
| `raw_panel_type` | TEXT | the vendor element name from which this row originated (audit trail) |
| `extra` | TEXT | free-form `key=value;...` overflow for vendor fields |
| `clinvar_id` | TEXT | GenMineTOP-only |
| `clinvar_url` | TEXT | GenMineTOP-only |
| `clinvar_sig` | TEXT | clinical significance |
| `clinvar_match` | TEXT | match level |
| `clinvar_benign` | INTEGER | count of benign assertions |
| `clinvar_likely_benign` | INTEGER | count |
| `clinvar_uncertain` | INTEGER | count |
| `maf_1kg` | REAL | 1000 Genomes minor allele frequency |
| `maf_hgvd` | REAL | HGVD MAF (Japanese population) |
| `maf_tommo` | REAL | ToMMo 8.3KJPN MAF |
| `tpm_normal_n` | INTEGER | reference normal-tissue cohort size |
| `tpm_normal_mean` | REAL | normal-tissue TPM mean |
| `tpm_normal_sd` | REAL | normal-tissue TPM standard deviation |

---

## 3. `non_human_contents` — one row per detected organism

Populated by `load_foundation.py` when the FoundationOne(Liquid) XML
contains a `<non-human-content>` element with per-organism children.

| Column | Type | Notes |
|---|---|---|
| `id` | INTEGER PK | |
| `case_id` | INTEGER NOT NULL | FK → `cases(case_id)` |
| `organism` | TEXT | e.g., HHV-4, HHV-8, HPV-16 |
| `reads_per_million` | REAL | |
| `status` | TEXT | `present` / `unknown` / ... |
| `sample` | TEXT | sample identifier from `<dna-evidence>` |

---

## 4. Indexes

The schema declares secondary indexes on the columns most often queried by
the CGIs:

- `variants(gene)`, `variants(variant_type)`, `variants(protein_effect)`
- `cases(panel_name)`, `cases(report_id)`, `cases(patient_id)`,
  `cases(date)`, `cases(ccat_date)`, `cases(ccat_cancer)`,
  `cases(disease)`, `cases(tissue_of_origin)`, `cases(pathology_diagnosis)`
- `non_human_contents(case_id)`

No full-text-search index is created by default; SQLite `LIKE` is fast
enough at the registry sizes (≤10⁵ cases) for which OncoUnify is intended.

---

## 5. Migration

`CREATE TABLE IF NOT EXISTS` does **not** add columns to an existing table.
When the canonical schema grows, ship an incremental migration as a separate
`migrations/YYYY-MM-DD-<description>.sql` file containing only
`ALTER TABLE` statements, and document it in [INSTALL.md §6](INSTALL.md).
