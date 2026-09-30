# OncoUnify — schema reference and data dictionary (schema version 2)

The authoritative definition is [`schema.sql`](../schema.sql)
(`PRAGMA user_version = 2`).  This document explains every column.  The
vendor-by-vendor origin of each column is tabulated in
[`field_mapping.tsv`](field_mapping.tsv).

## Conventions

* **One row in `cases` per vendor report.** `(panel_name, report_id)` is
  unique; loaders replace an existing report and all its child rows in one
  transaction.
* **Coordinates** (`chrom`, `pos`, `pos2`, `chrom2`) are 1-based and
  vendor-native on the assembly in `cases.genome_build`.  Chromosomes are
  UCSC-style (`chr1` … `chrX`, `chrY`, `chrM`).  Indels are **not**
  re-normalized (left/right justification and HGVS 3′ shifting are left as
  reported); compare indels across vendors by protein key or after
  normalizing with an external tool.
* **Verbatim and canonical side by side.** `protein_effect` and `cds_effect`
  hold exactly what the vendor wrote; `hgvs_p` and `hgvs_c` hold the
  canonical keys derived from them.  Nothing the vendor emits is dropped:
  fields without a canonical column go to the JSON objects `cases.other_info`
  and `variants.extra`.
* **Controlled vocabularies** are enforced with `CHECK` constraints or
  lookup tables, as noted below.  `NULL` means "not reported / not
  assessed", never zero.

## `so_terms` — consequence vocabulary

Sequence Ontology term names (as used by Ensembl VEP) allowed in
`variants.functional_effect`, with accession and the display group used for
colouring and filtering.

| term | accession | group |
|---|---|---|
| missense_variant | SO:0001583 | missense |
| stop_gained | SO:0001587 | nonsense |
| frameshift_variant | SO:0001589 | frameshift |
| splice_donor_variant | SO:0001575 | splice |
| splice_acceptor_variant | SO:0001574 | splice |
| splice_region_variant | SO:0001630 | splice |
| exon_loss_variant | SO:0001572 | splice |
| inframe_deletion | SO:0001822 | inframe |
| inframe_insertion | SO:0001821 | inframe |
| inframe_indel | SO:0001820 | inframe |
| start_lost | SO:0002012 | other |
| stop_lost | SO:0001578 | other |
| protein_altering_variant | SO:0001818 | other |
| coding_sequence_variant | SO:0001580 | other |
| upstream_gene_variant | SO:0001631 | other |
| sequence_variant | SO:0001060 | other |
| synonymous_variant | SO:0001819 | silent |
| stop_retained_variant | SO:0001567 | silent |
| 5_prime_UTR_variant | SO:0001623 | noncoding |
| 3_prime_UTR_variant | SO:0001624 | noncoding |
| UTR_variant | SO:0001622 | noncoding |
| intron_variant | SO:0001627 | noncoding |
| non_coding_transcript_variant | SO:0001619 | noncoding |

### How a term is assigned (`oncounify_core.classify_consequence`)

The same function is used by every loader:

1. **Vendor category**, when the vendor supplies one (Foundation Medicine
   `functional-effect`): `missense` → missense_variant, `nonsense` →
   stop_gained, `frameshift` → frameshift_variant, `nonframeshift` →
   inframe_deletion / inframe_insertion from the protein change, `splice` →
   splice_donor / splice_acceptor / splice_region / exon_loss from the
   intronic offsets in the coding change, `promoter` → upstream_gene_variant.
   `functional_effect_source = 'vendor'`.
2. **Protein change** (`hgvs_p`): `fs` → frameshift; `ext` or `*n<aa>` →
   stop_lost; `<aa>n*` → stop_gained (also for `p.M1*`, which ranks above
   start_lost in the Ensembl VEP severity order); other `M1…` → start_lost;
   `<aa>n=` → synonymous; `<aa>n<aa>` → missense; `del`/`ins`/`dup`/`delins`,
   including Foundation Medicine's substitution-style form `E746_S752>V` →
   in-frame deletion/insertion by net length, stop_gained if the replacement
   contains a stop codon, missense if the lengths are equal.
3. **Coding change** (`hgvs_c`): offsets ±1–2 → donor/acceptor (also for
   ranges that cross the exon boundary, e.g. `c.905-9_905del`), ±3–8 →
   splice_region, larger → intron; a range from one intron across whole exons
   to another → exon_loss; `c.-n` → 5′UTR; `c.*n` → 3′UTR; exonic indels
   (including substitution-style `c.2237_2255>T` and `c.2034G>CA`) →
   frameshift or in-frame from their net length; an exonic substitution with
   no protein change → coding_sequence_variant.

When the vendor category is "nonframeshift", the protein change refines it
(in-frame deletion or insertion, stop gained/lost, start lost).

`tools/qc_report.py` measures how often the rules in steps 2–6 reproduce
the vendor category for Foundation Medicine variants and lists the
discordant notations (`--examples`).
4. **Vendor hint** (Guardant `reporting_category`, GenMineTOP `type`):
   `promoter` → upstream_gene_variant, `utr` → UTR_variant, `non_coding` →
   non_coding_transcript_variant, RNA exon skipping → exon_loss_variant.
5. VCF-style ref/alt length of a coding variant → frameshift / in-frame.
6. Otherwise `sequence_variant` with `functional_effect_source =
   'unclassified'` — never `NULL` for a short variant.

Steps 2–5 give `functional_effect_source = 'inferred'`.  The rules operate on
notation only (no transcript model); see Limitations in the manuscript.

## `gene_symbol_map` — gene symbol harmonization

Vendors do not all use current HGNC symbols; Foundation Medicine, for
example, reports `WHSC1L1` where GenMineTOP reports `NSD3`, and older
Foundation Medicine reports use `MLL2` for `KMT2D`.  The table holds the
previous HGNC symbols of every gene in the supported assays' gene lists with
their current approved symbol — 599 symbols of 448 genes (HGNC complete set,
downloaded 2026-09-30).  A previous symbol is included only if it belongs to
exactly one gene and is neither the approved symbol nor an alias of another
gene (so that a search for `ERK` is not answered with `EPHB2`), except where
the vendors' gene lists give it as the synonym (`KMT2D (MLL2)`).  Aliases are
not resolved.

The table is filled by the loaders from `gene_symbol_map.tsv`, which
`tools/make_symbol_map.py` regenerates from a newer HGNC file or after gene
lists are added.  Loaders store the current symbol in `variants.gene` /
`other_gene` and `panel_genes.gene` (keeping the vendor symbol in
`variants.extra`), and the search interface resolves previous symbols, so a
gene is counted once across vendors.  After the map changes, reload the
reports (the loaders are idempotent) so that stored rows use it.

## `cases` — one row per report

| column | type | meaning |
|---|---|---|
| case_id | INTEGER PK | surrogate key |
| panel_name | TEXT NOT NULL | assay family used throughout the interface: `FoundationOne`, `FoundationOneLiquid`, `GenMineTOP`, `Guardant` |
| panel_version | TEXT | assay/gene-list version; with panel_name, joins `panel_versions` (e.g. `FoundationOneDx`, `TDv6+TRv6`, `Guardant360 CDx`) |
| panel_type | TEXT | vendor's own test-type label, verbatim |
| vendor | TEXT NOT NULL | `Foundation Medicine`, `GenMine Labs`, `Guardant Health` |
| report_id | TEXT NOT NULL | vendor report identifier (Guardant: from the file name) |
| patient_id | TEXT | institutional patient identifier or pseudonym; may be supplied by the sidecar |
| sex, age | TEXT, INTEGER | as reported |
| date | TEXT | report date `YYYY-MM-DD` (GenMineTOP `<accepted>`; others via sidecar) |
| genome_build | TEXT NOT NULL | `GRCh37` \| `GRCh38` \| `unknown` (the last only for migrated legacy rows) |
| genome_build_source | TEXT | `vendor-file` (read from the report), `loader-default` (vendor's documented assembly), `cli` (`--genome-build`), `migration` |
| disease | TEXT | vendor free-text disease label, or curated via sidecar |
| disease_ontology | TEXT | the vendor's own disease classification label (Foundation Medicine `disease-ontology`); not an ontology identifier |
| oncotree_code | TEXT | OncoTree code supplied by a curator through the sidecar |
| tissue_of_origin | TEXT | as reported, or curated |
| pathology_diagnosis | TEXT | as reported (GenMineTOP `specimen/pathology`; Foundation Medicine `pathology-diagnosis`) |
| specimen_id, test_type | TEXT | vendor specimen identifier; `test_type` duplicates `panel_type` for compatibility |
| percent_tumor_nuclei, purity | REAL | percent (0–100) |
| non_human_content | REAL | Foundation Medicine `variant-report/@non-human-content` (case-level scalar); per-organism detail is in `non_human_contents` |
| other_info | TEXT (JSON) | vendor case-level fields without a canonical column |
| curated_fields | TEXT | comma-separated names of fields supplied or overridden by the sidecar |
| source_file, source_sha256 | TEXT | provenance of the input file |
| format_version | TEXT | detected vendor format / schema / pipeline version |
| loader, loaded_at | TEXT | loader name and version; UTC timestamp of the last (re)load |

## `variants` — one row per reported finding

| column | type | meaning |
|---|---|---|
| variant_id | INTEGER PK | |
| case_id | INTEGER NOT NULL | FK → cases (ON DELETE CASCADE) |
| gene | TEXT | current HGNC symbol (5′ partner for fusions); previous symbols used by vendors are replaced through `gene_symbol_map` and the vendor symbol is kept in `extra.vendor_gene` |
| variant_type | TEXT NOT NULL | `short_variant` \| `cnv` \| `rearrangement` (fusions, rearrangements, RNA exon skipping) \| `expression` |
| variant_subtype | TEXT | vendor-native class label, verbatim (e.g. `amplification`, `fusion`, `splicing-variant`, `SNV`, `Deletion`) |
| chrom, pos | TEXT, INTEGER | position; start of a CNV segment; first breakpoint |
| pos2, chrom2 | INTEGER, TEXT | end of a CNV segment; second breakpoint |
| ref, alt | TEXT | alleles as reported (VCF-style with anchor base for GenMineTOP and Guardant) |
| transcript, strand | TEXT | as reported |
| cds_effect, protein_effect | TEXT | vendor strings, verbatim |
| hgvs_c | TEXT | `c.`-prefixed coding change (not validated against a transcript) |
| hgvs_p | TEXT | canonical protein key: `p.` + one-letter HGVS body (`p.G12D`, `p.T887Rfs*19`, `p.A999=`); NULL when no protein change is reported |
| functional_effect | TEXT | Sequence Ontology term (FK → so_terms); set for short variants and exon skipping |
| functional_effect_so | TEXT | SO accession |
| functional_effect_raw | TEXT | vendor category or hint the term was derived from |
| functional_effect_source | TEXT | `vendor` \| `inferred` \| `unclassified` |
| status | TEXT | vendor call status: Foundation Medicine `known`/`likely`/`unknown`; GenMineTOP `finding`/`notice`; Guardant `called` (or `not_called` with `--include-uncalled`) |
| origin | TEXT | `somatic` \| `germline`; NULL = not assessed (tumor-only assays) |
| classification | TEXT | vendor interpretation class (GenMineTOP `ag-class`) |
| allele_fraction | REAL | 0–1 (Guardant `percentage`/100; GenMineTOP alt/depth) |
| depth | INTEGER | read depth at the locus |
| copy_number, cnv_ratio | REAL | vendor-native values; **not harmonized** across vendors |
| cnv_type | TEXT | `amplification` \| `deletion` \| `other` |
| other_gene | TEXT | fusion / rearrangement partner (3′ partner); harmonized like `gene` |
| in_frame | TEXT | `yes` \| `no` \| `unknown` |
| supporting_read_pairs | INTEGER | DNA evidence for rearrangements (Foundation Medicine) |
| read_count | INTEGER | RNA reads supporting a fusion/exon skipping event, or expression read count |
| tpm | REAL | expression, transcripts per million (tumor) |
| sample_name | TEXT | vendor sample identifier for the evidence |
| effect | TEXT | vendor free-text description (Foundation Medicine rearrangement `description`; GenMineTOP exon-skipping junction) |
| raw_panel_type | TEXT | vendor element or sheet the row came from (audit trail) |
| extra | TEXT (JSON) | vendor variant fields without a canonical column (e.g. equivocal, subclonal, exon, reporting_category, cytoband, breakpoints) |
| clinvar_* | | GenMineTOP ClinVar annotation (id, URL, significance, match level, assertion counts) |
| maf_1kg, maf_hgvd, maf_tommo | REAL | GenMineTOP population allele frequencies |
| tpm_normal_n/mean/sd | | GenMineTOP normal-tissue expression reference |

## `biomarkers` — TMB and MSI

| column | meaning |
|---|---|
| name | `TMB` \| `MSI` |
| value | numeric score if reported (Foundation Medicine TMB; GenMineTOP exonic non-synonymous alteration frequency; Guardant MSI score) |
| unit | e.g. `mutations/Mb`, `Guardant MSI score` |
| call | TMB: `high` / `intermediate` / `low` / `indeterminate`; MSI: `MSI-H` / `MSS` (MSS and MSI-L) / `indeterminate` |
| call_raw | vendor verbatim call |
| assay | method identifier; values are comparable only within one assay (e.g. `FoundationOne CDx tissue TMB`, `FoundationOne Liquid CDx blood TMB (bTMB)`, `GenMineTOP exonic non-synonymous alteration frequency (tumor-normal paired)`, `Guardant360 CDx MSI (plasma cfDNA)`) |
| source_field | vendor field the value came from |

## `non_human_contents`

Per-organism rows from Foundation Medicine `non-human-content/non-human`
(organism, reads per million, status, sample).

## `panel_versions`, `panel_genes` — assay content

One `panel_versions` row per `(panel_name, panel_version)`; `panel_genes`
lists the genes it interrogates, with flags for short variants, copy number
and rearrangements.  Loaded from `panels/*.tsv` by `load_panel_genes.py`.
These are the denominators of the frequencies in `panel_stats.cgi`.

## Views

| view | content |
|---|---|
| v_short_variants | short variants joined to report identifiers and genome build |
| v_copy_number | CNV rows |
| v_rearrangements | fusions, rearrangements, exon skipping |
| v_expression | expression rows |
| v_case_genes_tested | one row per (case, gene) interrogated by the case's assay |

## Indexes

`variants(case_id)`, `variants(gene)`, `variants(other_gene)`,
`variants(hgvs_p)`, `variants(variant_type)`, `variants(functional_effect)`,
`cases(panel_name, panel_version)`, `cases(report_id)`, `cases(patient_id)`,
`cases(date)`, `cases(disease)`, `cases(oncotree_code)`,
`cases(tissue_of_origin)`, `cases(pathology_diagnosis)`,
`biomarkers(case_id)`, `non_human_contents(case_id)`, `panel_genes(gene)`.
Gene and protein filters are exact-match and use these indexes; free-text
disease filters are substring matches and scan `cases`.

## Upgrading

A version-1 database is refused by version-2 loaders.  Rebuild it from the
source reports (preferred), or upgrade it in place with
`python3 migrate_db.py panels.db` (a backup is kept).
