-- schema.sql
-- OncoUnify canonical schema, version 2 (PRAGMA user_version = 2).
--
-- The loaders execute this file on every run; every statement is idempotent
-- (CREATE ... IF NOT EXISTS).  CREATE TABLE IF NOT EXISTS cannot add columns
-- to an existing table, so a version-1 database must be upgraded with
-- `python3 migrate_db.py <db>` (or rebuilt from the source reports) before a
-- version-2 loader will write to it.
--
-- Conventions (see docs/SCHEMA.md for the full data dictionary):
--   * one row in `cases` per vendor report; (panel_name, report_id) is unique;
--   * coordinates (chrom/pos/pos2/chrom2) are 1-based, vendor-native, on the
--     assembly recorded in cases.genome_build; indels are NOT re-justified;
--   * `protein_effect`/`cds_effect` hold the vendor's verbatim strings, while
--     `hgvs_p`/`hgvs_c` hold canonical search keys derived from them;
--   * `functional_effect` is a Sequence Ontology term (SO accession in
--     `functional_effect_so`), assigned identically by every loader;
--   * `other_info` / `extra` are JSON objects holding vendor fields that have
--     no canonical column, so that nothing the vendor emits is dropped.

PRAGMA foreign_keys = ON;

-- ---------------------------------------------------------------------------
-- so_terms: controlled vocabulary for variants.functional_effect
-- (Sequence Ontology term names as used by Ensembl VEP).  display_group is
-- the coarse class used for colouring and filtering in the web interface.
-- Kept identical to oncounify_core.SO_TERMS / SO_GROUP (checked by the tests).
-- ---------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS so_terms (
    term           TEXT PRIMARY KEY,
    accession      TEXT NOT NULL,
    display_group  TEXT NOT NULL
                   CHECK (display_group IN ('missense', 'nonsense', 'frameshift', 'splice',
                                            'inframe', 'other', 'silent', 'noncoding'))
);
INSERT OR REPLACE INTO so_terms (term, accession, display_group) VALUES
    ('missense_variant',              'SO:0001583', 'missense'),
    ('stop_gained',                   'SO:0001587', 'nonsense'),
    ('frameshift_variant',            'SO:0001589', 'frameshift'),
    ('splice_donor_variant',          'SO:0001575', 'splice'),
    ('splice_acceptor_variant',       'SO:0001574', 'splice'),
    ('splice_region_variant',         'SO:0001630', 'splice'),
    ('exon_loss_variant',             'SO:0001572', 'splice'),
    ('inframe_deletion',              'SO:0001822', 'inframe'),
    ('inframe_insertion',             'SO:0001821', 'inframe'),
    ('inframe_indel',                 'SO:0001820', 'inframe'),
    ('start_lost',                    'SO:0002012', 'other'),
    ('stop_lost',                     'SO:0001578', 'other'),
    ('protein_altering_variant',      'SO:0001818', 'other'),
    ('coding_sequence_variant',       'SO:0001580', 'other'),
    ('upstream_gene_variant',         'SO:0001631', 'other'),
    ('sequence_variant',              'SO:0001060', 'other'),
    ('synonymous_variant',            'SO:0001819', 'silent'),
    ('stop_retained_variant',         'SO:0001567', 'silent'),
    ('5_prime_UTR_variant',           'SO:0001623', 'noncoding'),
    ('3_prime_UTR_variant',           'SO:0001624', 'noncoding'),
    ('UTR_variant',                   'SO:0001622', 'noncoding'),
    ('intron_variant',                'SO:0001627', 'noncoding'),
    ('non_coding_transcript_variant', 'SO:0001619', 'noncoding');

-- ---------------------------------------------------------------------------
-- cases: one row per vendor report
-- ---------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS cases (
    case_id               INTEGER PRIMARY KEY AUTOINCREMENT,
    panel_name            TEXT NOT NULL,     -- FoundationOne | FoundationOneLiquid | GenMineTOP | Guardant | ...
    panel_version         TEXT,              -- assay / gene-list version; joins panel_versions
    panel_type            TEXT,              -- vendor's own test-type label (verbatim)
    vendor                TEXT NOT NULL,     -- manufacturer / laboratory
    report_id             TEXT NOT NULL,
    patient_id            TEXT,
    sex                   TEXT,
    age                   INTEGER,
    date                  TEXT,              -- report date, YYYY-MM-DD

    genome_build          TEXT NOT NULL
                          CHECK (genome_build IN ('GRCh37', 'GRCh38', 'unknown')),
    genome_build_source   TEXT,              -- vendor-file | loader-default | cli | migration

    disease               TEXT,              -- vendor free-text disease label
    disease_ontology      TEXT,              -- vendor's own disease classification label (e.g. FMI disease-ontology)
    oncotree_code         TEXT,              -- OncoTree code supplied by a curator (case-metadata sidecar)
    tissue_of_origin      TEXT,
    pathology_diagnosis   TEXT,

    specimen_id           TEXT,
    test_type             TEXT,              -- duplicate of panel_type kept for compatibility
    percent_tumor_nuclei  REAL,              -- percent (0-100)
    purity                REAL,              -- percent (0-100)
    non_human_content     REAL,              -- FMI variant-report/@non-human-content (case-level scalar);
                                             -- per-organism detail is in non_human_contents

    other_info            TEXT,              -- JSON object: vendor case-level fields without a canonical column
    curated_fields        TEXT,              -- comma-separated fields supplied/overridden by the sidecar

    source_file           TEXT,              -- provenance: input file name
    source_sha256         TEXT,              -- provenance: SHA-256 of the input file
    format_version        TEXT,              -- provenance: vendor format / pipeline version detected
    loader                TEXT,              -- provenance: loader name and version
    loaded_at             TEXT,              -- provenance: UTC timestamp of the last (re)load

    UNIQUE (panel_name, report_id)
);

-- ---------------------------------------------------------------------------
-- variants: one row per reported molecular finding
-- ---------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS variants (
    variant_id               INTEGER PRIMARY KEY AUTOINCREMENT,
    case_id                  INTEGER NOT NULL REFERENCES cases(case_id) ON DELETE CASCADE,

    gene                     TEXT,
    variant_type             TEXT NOT NULL
                             CHECK (variant_type IN ('short_variant', 'cnv', 'rearrangement', 'expression')),
    variant_subtype          TEXT,           -- vendor-native class label (verbatim)

    chrom                    TEXT,           -- UCSC style (chr1 .. chrX, chrY, chrM)
    pos                      INTEGER,        -- 1-based; start of a CNV segment; first breakpoint
    pos2                     INTEGER,        -- end of a CNV segment; second breakpoint
    chrom2                   TEXT,           -- chromosome of the second breakpoint
    ref                      TEXT,
    alt                      TEXT,

    transcript               TEXT,
    strand                   TEXT,
    cds_effect               TEXT,           -- verbatim coding change from the vendor
    protein_effect           TEXT,           -- verbatim protein change from the vendor
    hgvs_c                   TEXT,           -- 'c.'-prefixed coding change (not transcript-validated)
    hgvs_p                   TEXT,           -- canonical protein key: 'p.' + one-letter HGVS body

    functional_effect        TEXT REFERENCES so_terms(term),  -- Sequence Ontology term name
    functional_effect_so     TEXT,           -- Sequence Ontology accession
    functional_effect_raw    TEXT,           -- vendor verbatim category, if any
    functional_effect_source TEXT
                             CHECK (functional_effect_source IN ('vendor', 'inferred', 'unclassified')),

    status                   TEXT,           -- vendor call status (e.g. known/likely/unknown; finding/notice)
    origin                   TEXT CHECK (origin IN ('somatic', 'germline')),  -- NULL = not assessed
    classification           TEXT,           -- vendor/curator interpretation class
    allele_fraction          REAL,           -- 0-1
    depth                    INTEGER,

    copy_number              REAL,           -- vendor-native estimate (not harmonized across vendors)
    cnv_ratio                REAL,
    cnv_type                 TEXT CHECK (cnv_type IN ('amplification', 'deletion', 'other')),

    other_gene               TEXT,           -- fusion / rearrangement partner
    in_frame                 TEXT CHECK (in_frame IN ('yes', 'no', 'unknown')),
    supporting_read_pairs    INTEGER,

    tpm                      REAL,
    read_count               INTEGER,
    sample_name              TEXT,

    effect                   TEXT,           -- vendor free-text description (e.g. FMI rearrangement/@description)
    raw_panel_type           TEXT,           -- vendor element / sheet the row came from (audit trail)
    extra                    TEXT,           -- JSON object: vendor variant fields without a canonical column

    clinvar_id               TEXT,
    clinvar_url              TEXT,
    clinvar_sig              TEXT,
    clinvar_match            TEXT,
    clinvar_benign           INTEGER,
    clinvar_likely_benign    INTEGER,
    clinvar_uncertain        INTEGER,

    maf_1kg                  REAL,
    maf_hgvd                 REAL,
    maf_tommo                REAL,

    tpm_normal_n             INTEGER,
    tpm_normal_mean          REAL,
    tpm_normal_sd            REAL
);

-- ---------------------------------------------------------------------------
-- biomarkers: TMB / MSI with their unit and assay.  Values are comparable
-- only within the same `assay`; they are never pooled into one column.
-- ---------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS biomarkers (
    biomarker_id   INTEGER PRIMARY KEY AUTOINCREMENT,
    case_id        INTEGER NOT NULL REFERENCES cases(case_id) ON DELETE CASCADE,
    name           TEXT NOT NULL CHECK (name IN ('TMB', 'MSI')),
    value          REAL,             -- numeric score, if reported
    unit           TEXT,             -- e.g. 'mutations/Mb', 'Guardant MSI score'
    call           TEXT,             -- TMB: high/intermediate/low/indeterminate; MSI: MSI-H/MSS/indeterminate
    call_raw       TEXT,             -- vendor verbatim call
    assay          TEXT NOT NULL,    -- method identifier, e.g. 'FoundationOne CDx tissue TMB'
    source_field   TEXT              -- vendor field the value came from (audit trail)
);

-- ---------------------------------------------------------------------------
-- non_human_contents: per-organism rows (FoundationOne Liquid CDx)
-- ---------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS non_human_contents (
    id                    INTEGER PRIMARY KEY AUTOINCREMENT,
    case_id               INTEGER NOT NULL REFERENCES cases(case_id) ON DELETE CASCADE,
    organism              TEXT,              -- e.g. HHV-4, HHV-8, HPV-16
    reads_per_million     REAL,
    status                TEXT,
    sample                TEXT
);

-- ---------------------------------------------------------------------------
-- Panel content (denominators for cross-panel statistics).  Loaded from
-- panels/*.tsv by load_panel_genes.py.  A case is joined to its gene list
-- through (cases.panel_name, cases.panel_version).
-- ---------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS panel_versions (
    panel_version_id  INTEGER PRIMARY KEY AUTOINCREMENT,
    panel_name        TEXT NOT NULL,
    panel_version     TEXT NOT NULL,
    description       TEXT,
    source            TEXT,             -- provenance of the gene list (document, version, date)
    UNIQUE (panel_name, panel_version)
);

CREATE TABLE IF NOT EXISTS panel_genes (
    panel_version_id  INTEGER NOT NULL REFERENCES panel_versions(panel_version_id) ON DELETE CASCADE,
    gene              TEXT NOT NULL,
    short_variants    INTEGER NOT NULL DEFAULT 1,   -- SNV/indel interrogated
    copy_number       INTEGER NOT NULL DEFAULT 0,   -- CNA interrogated
    rearrangements    INTEGER NOT NULL DEFAULT 0,   -- fusion/rearrangement interrogated
    PRIMARY KEY (panel_version_id, gene)
);

-- ---------------------------------------------------------------------------
-- Indexes
-- ---------------------------------------------------------------------------
CREATE INDEX IF NOT EXISTS idx_variants_case         ON variants(case_id);
CREATE INDEX IF NOT EXISTS idx_variants_gene         ON variants(gene);
CREATE INDEX IF NOT EXISTS idx_variants_other_gene   ON variants(other_gene);
CREATE INDEX IF NOT EXISTS idx_variants_hgvs_p       ON variants(hgvs_p);
CREATE INDEX IF NOT EXISTS idx_variants_vtype        ON variants(variant_type);
CREATE INDEX IF NOT EXISTS idx_variants_feffect      ON variants(functional_effect);

CREATE INDEX IF NOT EXISTS idx_cases_panel           ON cases(panel_name, panel_version);
CREATE INDEX IF NOT EXISTS idx_cases_report          ON cases(report_id);
CREATE INDEX IF NOT EXISTS idx_cases_patient         ON cases(patient_id);
CREATE INDEX IF NOT EXISTS idx_cases_date            ON cases(date);
CREATE INDEX IF NOT EXISTS idx_cases_disease         ON cases(disease);
CREATE INDEX IF NOT EXISTS idx_cases_oncotree        ON cases(oncotree_code);
CREATE INDEX IF NOT EXISTS idx_cases_tissue          ON cases(tissue_of_origin);
CREATE INDEX IF NOT EXISTS idx_cases_pathology       ON cases(pathology_diagnosis);

CREATE INDEX IF NOT EXISTS idx_biomarkers_case       ON biomarkers(case_id);
CREATE INDEX IF NOT EXISTS idx_nonhuman_case         ON non_human_contents(case_id);
CREATE INDEX IF NOT EXISTS idx_panel_genes_gene      ON panel_genes(gene);

-- ---------------------------------------------------------------------------
-- Typed read-only views over the sparse `variants` table
-- ---------------------------------------------------------------------------
CREATE VIEW IF NOT EXISTS v_short_variants AS
SELECT c.panel_name, c.report_id, c.patient_id, c.genome_build,
       v.variant_id, v.case_id, v.gene, v.variant_subtype, v.chrom, v.pos, v.ref, v.alt,
       v.transcript, v.cds_effect, v.protein_effect, v.hgvs_c, v.hgvs_p,
       v.functional_effect, v.functional_effect_so, v.functional_effect_source,
       v.status, v.origin, v.allele_fraction, v.depth,
       v.clinvar_id, v.clinvar_sig, v.maf_1kg, v.maf_hgvd, v.maf_tommo
FROM variants v JOIN cases c ON c.case_id = v.case_id
WHERE v.variant_type = 'short_variant';

CREATE VIEW IF NOT EXISTS v_copy_number AS
SELECT c.panel_name, c.report_id, c.patient_id, c.genome_build,
       v.variant_id, v.case_id, v.gene, v.variant_subtype, v.chrom, v.pos, v.pos2,
       v.cnv_type, v.copy_number, v.cnv_ratio, v.status
FROM variants v JOIN cases c ON c.case_id = v.case_id
WHERE v.variant_type = 'cnv';

CREATE VIEW IF NOT EXISTS v_rearrangements AS
SELECT c.panel_name, c.report_id, c.patient_id, c.genome_build,
       v.variant_id, v.case_id, v.gene, v.other_gene, v.variant_subtype,
       v.chrom, v.pos, v.chrom2, v.pos2, v.in_frame, v.supporting_read_pairs,
       v.allele_fraction, v.read_count, v.functional_effect, v.effect, v.status
FROM variants v JOIN cases c ON c.case_id = v.case_id
WHERE v.variant_type = 'rearrangement';

CREATE VIEW IF NOT EXISTS v_expression AS
SELECT c.panel_name, c.report_id, c.patient_id,
       v.variant_id, v.case_id, v.gene, v.tpm, v.read_count,
       v.tpm_normal_n, v.tpm_normal_mean, v.tpm_normal_sd, v.status
FROM variants v JOIN cases c ON c.case_id = v.case_id
WHERE v.variant_type = 'expression';

-- One row per (case, gene) that the case's assay interrogates for short variants.
CREATE VIEW IF NOT EXISTS v_case_genes_tested AS
SELECT c.case_id, pg.gene, pg.short_variants, pg.copy_number, pg.rearrangements
FROM cases c
JOIN panel_versions pv ON pv.panel_name = c.panel_name AND pv.panel_version = c.panel_version
JOIN panel_genes pg    ON pg.panel_version_id = pv.panel_version_id;

PRAGMA user_version = 2;
