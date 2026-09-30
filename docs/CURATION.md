# Curating case metadata (the sidecar file)

Vendor deliverables do not always carry what cross-vendor questions need:
Guardant workbooks contain no disease label, Foundation Medicine reports may
leave the medical record number empty, and free-text diagnoses differ between
vendors and sites.  OncoUnify therefore accepts a curator-maintained CSV that
the loaders read at load time:

```bash
python3 load_guardant.py panels.db /reports/guardant --case-metadata curation.csv
```

## Format

UTF-8 CSV with a header.  `panel_name` and `report_id` identify the report
(the same values the loader records).  All other columns are optional;
non-empty cells override the vendor value, empty cells leave it unchanged.

| column | use |
|---|---|
| panel_name, report_id | key (required) |
| patient_id | link reports of the same patient (e.g. tissue and liquid biopsy) |
| disease | harmonized disease label |
| oncotree_code | OncoTree code, e.g. `COAD`, `LUAD` |
| tissue_of_origin, pathology_diagnosis | harmonized values |
| sex, age, date | when missing from the vendor file |
| note | free text for curators; ignored by the loaders |

Example (`tests/data/case_metadata.csv`):

```csv
panel_name,report_id,patient_id,disease,oncotree_code,tissue_of_origin,note
FoundationOne,SYN-F1-0001,SYNPT-0001,,COAD,,MRN absent in the vendor file
Guardant,SYN-G360-0001,,Colon adenocarcinoma,COAD,Colon,disease is not part of the Guardant deliverable
```

Fields taken from the sidecar are listed in `cases.curated_fields` and are
marked "(curated)" on the case-detail page.  Because loading a report
replaces it, editing the sidecar and re-running the loader updates the
database; removing a row restores the vendor values.

## Scope

The sidecar is a manual curation path.  OncoUnify does not yet map free-text
diagnoses to OncoTree or Disease Ontology automatically; see the project
roadmap.
