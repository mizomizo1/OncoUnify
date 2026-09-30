-- schema.sql
-- OncoUnify: unified cancer genomic panel database schema
--
-- NOTE:
--   CREATE TABLE IF NOT EXISTS does not add new columns to an existing DB.
--   If you need to add columns to an existing database, run the bundled
--   migrate_db.py to apply the necessary ALTER TABLE statements.

PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS cases (
    case_id               INTEGER PRIMARY KEY AUTOINCREMENT,
    panel_name            TEXT NOT NULL,
    panel_type            TEXT,
    vendor                TEXT,
    report_id             TEXT,
    patient_id            TEXT,
    sex                   TEXT,
    age                   INTEGER,

    -- Disease / tissue information supplied by the institution (e.g., parsed from each vendor XML)
    disease               TEXT,
    disease_ontology      TEXT,
    tissue_of_origin      TEXT,
    pathology_diagnosis   TEXT,

    specimen_id           TEXT,
    test_type             TEXT,

    percent_tumor_nuclei  REAL,
    purity                REAL,
    msi_status            TEXT,
    tmb_score             REAL,
    tmb_status            TEXT,
    tmb_unit              TEXT,

    -- Test date (e.g., the report-creation or <accepted> date from each vendor)
    date                  TEXT,

    -- C-CAT-derived information (e.g., extracted from c-cat-f PDF)
    ccat_date             TEXT,   -- specimen collection date
    ccat_cancer           TEXT,   -- cancer type

    -- FoundationOne(Liquid) non-human-content attribute (case-level scalar)
    non_human_content     REAL,

    other_info            TEXT,

    UNIQUE(panel_name, report_id)
);

CREATE TABLE IF NOT EXISTS variants (
    variant_id            INTEGER PRIMARY KEY AUTOINCREMENT,
    case_id               INTEGER NOT NULL,

    -- Generic variant attributes
    gene                  TEXT,
    variant_type          TEXT,   -- short_variant / cnv / rearrangement / expression / biomarker, etc.
    variant_subtype       TEXT,
    chrom                 TEXT,
    pos                   INTEGER,
    pos2                  INTEGER,    -- second breakpoint for rearrangement

    ref                   TEXT,
    alt                   TEXT,
    cds_effect            TEXT,
    protein_effect        TEXT,

    -- FoundationOne / GenMineTOP: strand and transcript
    strand                TEXT,
    transcript            TEXT,

    -- Predicted protein-level consequence (missense / nonsense / frameshift / splice / synonymous, etc.)
    functional_effect     TEXT,
    effect                TEXT,
    status                TEXT,
    origin                TEXT,
    classification        TEXT,
    allele_fraction       REAL,
    depth                 INTEGER,
    copy_number           REAL,
    cnv_ratio             REAL,
    cnv_type              TEXT,

    -- Rearrangement-specific (FoundationOne)
    other_gene            TEXT,
    in_frame              TEXT,
    supporting_read_pairs INTEGER,

    -- TPM (case-side expression)
    tpm                   REAL,
    read_count            INTEGER,
    sample_name           TEXT,

    raw_panel_type        TEXT,
    extra                 TEXT,

    -- GenMineTOP: ClinVar-related
    clinvar_id            TEXT,
    clinvar_url           TEXT,
    clinvar_sig           TEXT,
    clinvar_match         TEXT,
    clinvar_benign        INTEGER,
    clinvar_likely_benign INTEGER,
    clinvar_uncertain     INTEGER,

    -- GenMineTOP: MAF (population allele frequencies)
    maf_1kg               REAL,
    maf_hgvd              REAL,
    maf_tommo             REAL,

    -- GenMineTOP: normal-tissue TPM reference distribution
    tpm_normal_n          INTEGER,
    tpm_normal_mean       REAL,
    tpm_normal_sd         REAL,

    FOREIGN KEY(case_id) REFERENCES cases(case_id)
);

-- Per-organism rows derived from the FoundationOne(Liquid) non-human-content element
CREATE TABLE IF NOT EXISTS non_human_contents (
    id                    INTEGER PRIMARY KEY AUTOINCREMENT,
    case_id               INTEGER NOT NULL,
    organism              TEXT,   -- e.g., HHV-4, HHV-8, HPV-16
    reads_per_million     REAL,   -- reads-per-million
    status                TEXT,   -- unknown, present, etc.
    sample                TEXT,   -- sample ID from dna-evidence

    FOREIGN KEY(case_id) REFERENCES cases(case_id)
);

-- indexes
CREATE INDEX IF NOT EXISTS idx_variants_gene         ON variants(gene);
CREATE INDEX IF NOT EXISTS idx_variants_vtype        ON variants(variant_type);
CREATE INDEX IF NOT EXISTS idx_variants_protein      ON variants(protein_effect);

CREATE INDEX IF NOT EXISTS idx_cases_panel           ON cases(panel_name);
CREATE INDEX IF NOT EXISTS idx_cases_report          ON cases(report_id);
CREATE INDEX IF NOT EXISTS idx_cases_patient         ON cases(patient_id);
CREATE INDEX IF NOT EXISTS idx_cases_date            ON cases(date);
CREATE INDEX IF NOT EXISTS idx_cases_ccat_date       ON cases(ccat_date);
CREATE INDEX IF NOT EXISTS idx_cases_ccat_cancer     ON cases(ccat_cancer);

CREATE INDEX IF NOT EXISTS idx_cases_disease         ON cases(disease);
CREATE INDEX IF NOT EXISTS idx_cases_tissue          ON cases(tissue_of_origin);
CREATE INDEX IF NOT EXISTS idx_cases_pathology       ON cases(pathology_diagnosis);

CREATE INDEX IF NOT EXISTS idx_nonhuman_case         ON non_human_contents(case_id);
