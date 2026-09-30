#!/bin/sh
set -eu
DB="${ONCOUNIFY_DB:-/data/panels.db}"
export ONCOUNIFY_DB="$DB"
mkdir -p "$(dirname "$DB")"
cd /opt/oncounify

if [ "${ONCOUNIFY_DEMO:-0}" = "1" ] && [ ! -s "$DB" ]; then
    echo "[oncounify] building the demo database from the synthetic fixtures"
    ONCOUNIFY_REPORTS=/opt/oncounify/tests/data ./docker/load_reports.sh \
        || echo "[oncounify] WARNING: the demo load reported failures"
fi

if [ "${ONCOUNIFY_LOAD:-0}" = "1" ]; then
    echo "[oncounify] loading reports from ${ONCOUNIFY_REPORTS:-/reports}"
    ./docker/load_reports.sh || echo "[oncounify] WARNING: some reports failed to load; see the summary above"
fi

if [ -f /etc/oncounify/htpasswd ]; then
    cat > /etc/apache2/conf-enabled/oncounify-auth.conf <<'CONF'
<Location />
    AuthType Basic
    AuthName "OncoUnify"
    AuthUserFile /etc/oncounify/htpasswd
    Require valid-user
</Location>
CONF
fi

if [ -f "$DB" ]; then
    chown www-data:www-data "$DB" || true
else
    echo "[oncounify] WARNING: $DB does not exist yet; set ONCOUNIFY_DEMO=1 or ONCOUNIFY_LOAD=1"
fi
exec "$@"
