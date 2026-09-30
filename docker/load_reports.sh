#!/bin/sh
# Load every vendor directory found under $ONCOUNIFY_REPORTS into $ONCOUNIFY_DB.
#   $ONCOUNIFY_REPORTS/foundation/   FoundationOne CDx / Liquid CDx XML
#   $ONCOUNIFY_REPORTS/genminetop/   GenMineTOP XML
#   $ONCOUNIFY_REPORTS/guardant/     Guardant360 CDx XLSX
#   $ONCOUNIFY_REPORTS/panels/       assay gene lists (else ./panels/*.tsv)
#   $ONCOUNIFY_REPORTS/case_metadata.csv   curation sidecar (optional)
# Exits 1 if any loader reported a failed file.
set -u
cd "$(dirname "$0")/.."
DB="${ONCOUNIFY_DB:-/data/panels.db}"
R="${ONCOUNIFY_REPORTS:-/reports}"
set --
[ -f "$R/case_metadata.csv" ] && set -- --case-metadata "$R/case_metadata.csv"
status=0
if [ -d "$R/foundation" ]; then python3 load_foundation.py "$DB" "$R/foundation" "$@" --quiet || status=1; fi
if [ -d "$R/genminetop" ]; then python3 load_genminetop.py "$DB" "$R/genminetop" "$@" --quiet || status=1; fi
if [ -d "$R/guardant" ];   then python3 load_guardant.py   "$DB" "$R/guardant"   "$@" --quiet || status=1; fi
if [ -d "$R/panels" ]; then
    python3 load_panel_genes.py "$DB" "$R/panels" || status=1
elif ls panels/*.tsv >/dev/null 2>&1; then
    python3 load_panel_genes.py "$DB" panels || status=1
fi
exit $status
