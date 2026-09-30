#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Write docs/field_mapping.tsv — the complete vendor-to-canonical field mapping
(Supplementary Table S1 of the manuscript).  tests/test_core.py checks that
every canonical column of oncounify_core is listed here, so that the
documentation cannot silently drift from the code.

    python3 docs/make_field_mapping.py
"""

from pathlib import Path

NA = "—"
FM = "Foundation Medicine XML"
GM = "GenMineTOP XML"
GD = "Guardant360 CDx XLSX"

ROWS = [
    # table, column, Foundation Medicine, GenMineTOP, Guardant, notes
    ("cases", "panel_name", "FoundationOneLiquid if variant-report/@test-type contains 'liquid', else FoundationOne",
     "--panel-name (default GenMineTOP)", "--panel-name (default Guardant)", ""),
    ("cases", "panel_version", "variant-report/@test-type", "report/test/panel/item (joined with '+')",
     "--panel-version (default 'Guardant360 CDx')", "joins panel_versions for denominators"),
    ("cases", "panel_type", "variant-report/@test-type", "report/test/type", "= panel_version", "verbatim"),
    ("cases", "test_type", "= panel_type", "= panel_type", "= panel_type", "kept for compatibility"),
    ("cases", "vendor", "'Foundation Medicine'", "'GenMine Labs'", "'Guardant Health'", "constant"),
    ("cases", "report_id", "rr:CustomerInformation/rr:ReferenceID", "report/id",
     "file name, named group report_id of --filename-pattern", "required; a report without it fails"),
    ("cases", "patient_id", "rr:CustomerInformation/rr:MRN", "report/patient/id",
     "file name, named group patient_id of --filename-pattern", "sidecar may supply/override"),
    ("cases", "sex", "variant-report/@gender", "report/patient/sex", NA, "sidecar may supply"),
    ("cases", "age", NA, "report/patient/age", NA, "sidecar may supply"),
    ("cases", "date", NA, "report/date/accepted (normalized to YYYY-MM-DD)", NA, "sidecar may supply"),
    ("cases", "genome_build", "GRCh37 (loader default; not stated in the file)",
     "report/result/reference-genome (fallback input/options/genome-version)", "GRCh37 (loader default)",
     "--genome-build overrides; NOT NULL"),
    ("cases", "genome_build_source", "'loader-default'", "'vendor-file'", "'loader-default'",
     "'cli' when --genome-build is used"),
    ("cases", "disease", "variant-report/@disease", NA, NA, "sidecar may supply/override"),
    ("cases", "disease_ontology", "variant-report/@disease-ontology", NA, NA,
     "vendor's own classification label, not an ontology identifier"),
    ("cases", "oncotree_code", NA, NA, NA, "curator-supplied through the sidecar"),
    ("cases", "tissue_of_origin", "variant-report/@tissue-of-origin", NA, NA, "sidecar may supply/override"),
    ("cases", "pathology_diagnosis", "variant-report/@pathology-diagnosis", "report/specimen/pathology", NA,
     "sidecar may supply/override"),
    ("cases", "specimen_id", "variant-report/@specimen", "report/specimen/id", NA, ""),
    ("cases", "percent_tumor_nuclei", "variant-report/@percent-tumor-nuclei", "report/qc/tumor-content/nuclei/value",
     NA, "percent"),
    ("cases", "purity", "variant-report/@purity-assessment", "report/qc/tumor-content/estimated/value", NA, "percent"),
    ("cases", "non_human_content", "variant-report/@non-human-content", NA, NA, "case-level scalar"),
    ("cases", "other_info",
     "JSON: @flowcell-analysis, @pipeline-version, @study, @test-request, quality-control/@status, samples/sample",
     "JSON: patient/c-cat-id, preference/germline-disclosure, marker/tmb/exon/num-non-synonymous-alterations, "
     "reference/db and reference/snp-db versions",
     "JSON: MSI sheet runid, run_sample_id, max_maf, mean_maf; QC status counts; number of uncalled rows skipped", ""),
    ("cases", "curated_fields", "sidecar", "sidecar", "sidecar", "names of fields taken from the sidecar"),
    ("cases", "source_file", "loader", "loader", "loader", "provenance"),
    ("cases", "source_sha256", "loader", "loader", "loader", "provenance"),
    ("cases", "format_version", "xsd names of ResultsReport and variant-report; @pipeline-version",
     "vendor-id and root/version", "sheet names present", "provenance; unknown versions raise a warning"),
    ("cases", "loader", "loader", "loader", "loader", "provenance"),
    ("cases", "loaded_at", "loader", "loader", "loader", "provenance"),

    ("variants", "gene", "short-variant/@gene; copy-number-alteration/@gene; rearrangement/@targeted-gene",
     "alterations/item/gene (first item for fusions)", "SNV, Indels, CNAs: gene; Fusion(s): gene_a",
     "previous HGNC symbols replaced by the current symbol (gene_symbol_map); vendor symbol kept in extra"),
    ("variants", "variant_type", "element: short-variant, copy-number-alteration, rearrangement",
     "item/type: snv, insertion, deletion, indel, delins, mnv -> short_variant; cnv-* -> cnv; fusion, "
     "splicing-variant -> rearrangement; expression -> expression",
     "sheet: SNV, Indels -> short_variant; CNAs -> cnv; Fusion(s) -> rearrangement", "controlled"),
    ("variants", "variant_subtype", "short-variant/@functional-effect; copy-number-alteration/@type; rearrangement/@type",
     "item/type", "'SNV'; Indels: type; 'CNA'; 'fusion'", "verbatim"),
    ("variants", "chrom", "@position / @pos1 (chr part)", "item/locus; fusion: breakpoint/item[1]/chr",
     "chrom (prefixed with chr); Fusion(s): chrom_a", "UCSC style"),
    ("variants", "pos", "short-variant/@position; CNA @position start; rearrangement/@pos1",
     "item/locus; CNV locus start; breakpoint/item[1]/pos", "position; Fusion(s): pos_a (0 -> NULL)", "1-based"),
    ("variants", "pos2", "CNA @position end; rearrangement/@pos2", "CNV locus end; breakpoint/item[2]/pos",
     "Fusion(s): pos_b (0 -> NULL)", ""),
    ("variants", "chrom2", "rearrangement/@pos2 (chr part)", "breakpoint/item[2]/chr", "Fusion(s): chrom_b", ""),
    ("variants", "ref", NA, "item/ref", "mut_nt (left of '>')", ""),
    ("variants", "alt", NA, "item/alt", "mut_nt (right of '>')", ""),
    ("variants", "transcript", "short-variant/@transcript", "item/transcript (first item)", "transcript_id", ""),
    ("variants", "strand", "short-variant/@strand", NA, NA, ""),
    ("variants", "cds_effect", "short-variant/@cds-effect", "item/coding-dna-alteration", "cdna", "verbatim"),
    ("variants", "protein_effect", "short-variant/@protein-effect", "item/protein-alteration",
     "SNV: mut_aa; Indels: mut_aa or mut_aa_short", "verbatim"),
    ("variants", "hgvs_c", "derived from cds_effect", "derived from cds_effect", "derived from cds_effect",
     "'c.' prefix; not transcript-validated"),
    ("variants", "hgvs_p", "derived from protein_effect", "derived from protein_effect", "derived from protein_effect",
     "canonical one-letter key; search key"),
    ("variants", "functional_effect", "@functional-effect mapped to SO (splice refined from cds offsets)",
     "inferred from protein / coding change / item type; splicing-variant -> exon_loss_variant",
     "inferred from mut_aa / cdna / reporting_category", "Sequence Ontology term (so_terms)"),
    ("variants", "functional_effect_so", "derived", "derived", "derived", "SO accession"),
    ("variants", "functional_effect_raw", "@functional-effect", "item/type or sv-type", "reporting_category",
     "what the term was derived from"),
    ("variants", "functional_effect_source", "'vendor'", "'inferred' ('vendor' for exon skipping)", "'inferred'",
     "'unclassified' when nothing was informative"),
    ("variants", "status", "@status (known / likely / unknown)", "item/status (finding / notice)",
     "'called' (call = 1); 'not_called' only with --include-uncalled", ""),
    ("variants", "origin", NA, "item/origin (somatic / germline)", NA, "NULL = not assessed (tumor-only / plasma)"),
    ("variants", "classification", NA, "item/ag-class", NA, ""),
    ("variants", "allele_fraction", "@allele-fraction", "item/allele-frequency 'alt/depth'", "percentage / 100", "0-1"),
    ("variants", "depth", "short-variant/@depth", "denominator of item/allele-frequency", NA, ""),
    ("variants", "copy_number", "copy-number-alteration/@copy-number", "item/num-copy", "CNAs: copy_number",
     "vendor-native; not harmonized"),
    ("variants", "cnv_ratio", "copy-number-alteration/@ratio", "item/ratio", NA, "vendor-native"),
    ("variants", "cnv_type", "copy-number-alteration/@type (amplification / loss -> deletion)",
     "item/type (cnv-amplification -> amplification)", "'amplification' for called CNAs", "controlled"),
    ("variants", "other_gene", "rearrangement/@other-gene", "alterations/item/gene (second item)", "gene_b",
     "fusion partner; symbol harmonized as for gene"),
    ("variants", "in_frame", "rearrangement/@in-frame", "item/frame", NA, "yes / no / unknown"),
    ("variants", "supporting_read_pairs", "rearrangement/@supporting-read-pairs", NA, NA, "DNA evidence"),
    ("variants", "tpm", NA, "item/tpm (expression)", NA, ""),
    ("variants", "read_count", NA, "item/num-reads (fusion, exon skipping, expression)", NA, "RNA evidence"),
    ("variants", "sample_name", "rearrangement/dna-evidence/@sample", NA, NA, ""),
    ("variants", "effect", "rearrangement/@description", "exon-skipping junction (breakpoint exon indices)", NA,
     "free text"),
    ("variants", "raw_panel_type", "element name", "item/type", "sheet name", "audit trail"),
    ("variants", "extra", "JSON: @equivocal, @subclonal, @percent-reads, @number-of-exons, dna-evidence/@sample",
     "JSON: cytoband, id/vendor, dbs/cosmic, raw allele-frequency(-tumor), breakpoints, num-wt-reads, sv-type, "
     "transcripts", "JSON: exon, length, reporting_category, rm_reportable, percentage, downstream_gene", ""),
    ("variants", "clinvar_id", NA, "item/clinical-relevance/dbs/clinvar/id/item", NA, ""),
    ("variants", "clinvar_url", NA, "item/dbs/clinvar/item", NA, ""),
    ("variants", "clinvar_sig", NA, "item/clinical-relevance/dbs/clinvar/clinical-significance/item", NA, ""),
    ("variants", "clinvar_match", NA, "item/clinical-relevance/dbs/clinvar/match/item", NA, ""),
    ("variants", "clinvar_benign", NA, "item/clinical-relevance/content/benign", NA, "count"),
    ("variants", "clinvar_likely_benign", NA, "item/clinical-relevance/content/likely-benign", NA, "count"),
    ("variants", "clinvar_uncertain", NA, "item/clinical-relevance/content/uncertain-significance", NA, "count"),
    ("variants", "maf_1kg", NA, "item/clinical-relevance/minor-allele-frequency/key[@name='1000genomes']", NA, ""),
    ("variants", "maf_hgvd", NA, "item/clinical-relevance/minor-allele-frequency/hgvd", NA, ""),
    ("variants", "maf_tommo", NA, "item/clinical-relevance/minor-allele-frequency/tommo-8p3kjpn", NA, ""),
    ("variants", "tpm_normal_n", NA, "item/reference/normal-expression/tpm/n", NA, ""),
    ("variants", "tpm_normal_mean", NA, "item/reference/normal-expression/tpm/mean", NA, ""),
    ("variants", "tpm_normal_sd", NA, "item/reference/normal-expression/tpm/sd", NA, ""),

    ("biomarkers", "name", "'MSI', 'TMB'", "'TMB'", "'MSI'", "controlled"),
    ("biomarkers", "value", "TMB: biomarkers/tumor-mutation-burden/@score",
     "TMB: result/marker/tmb/exon/frequency-non-synonymous-alterations", "MSI: MSI sheet msi_score", ""),
    ("biomarkers", "unit", "TMB: tumor-mutation-burden/@unit (-> mutations/Mb)", "TMB: 'mutations/Mb'",
     "MSI: 'Guardant MSI score'", ""),
    ("biomarkers", "call", "MSI: microsatellite-instability/@status; TMB: tumor-mutation-burden/@status", NA,
     "MSI: MSI sheet msi_status", "MSI-H / MSS / indeterminate; high / intermediate / low / indeterminate"),
    ("biomarkers", "call_raw", "verbatim @status", NA, "verbatim msi_status", ""),
    ("biomarkers", "assay", "'FoundationOne CDx tissue TMB', 'FoundationOne Liquid CDx blood TMB (bTMB)', "
     "'FoundationOne (Liquid) CDx MSI'", "'GenMineTOP exonic non-synonymous alteration frequency (tumor-normal paired)'",
     "'Guardant360 CDx MSI (plasma cfDNA)'", "values comparable only within one assay"),
    ("biomarkers", "source_field", "loader", "loader", "loader", "audit trail"),

    ("non_human_contents", "organism", "non-human-content/non-human/@organism", NA, NA, ""),
    ("non_human_contents", "reads_per_million", "non-human-content/non-human/@reads-per-million", NA, NA, ""),
    ("non_human_contents", "status", "non-human-content/non-human/@status", NA, NA, ""),
    ("non_human_contents", "sample", "non-human-content/non-human/dna-evidence/@sample", NA, NA, ""),
]


def main():
    out = Path(__file__).with_name("field_mapping.tsv")
    lines = ["\t".join(("table", "column", FM, GM, GD, "notes"))]
    for r in ROWS:
        assert len(r) == 6, r
        lines.append("\t".join(c.replace("\t", " ") for c in r))
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {out} ({len(ROWS)} rows)")


if __name__ == "__main__":
    main()
