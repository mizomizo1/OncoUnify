# OncoUnify — installation and deployment

OncoUnify has two parts that share one SQLite file:

1. **Ingestion** (the core of the package): Python loaders that normalize
   vendor reports into `panels.db`.  Needs Python ≥ 3.9 (`pandas` and
   `openpyxl` for Guardant workbooks).  Runs anywhere; no server required.
2. **Reference web interface** (optional): four read-only Perl CGI scripts
   and a static page.  Needs Perl 5 with `CGI`, `DBI`, `DBD::SQLite`
   (`JSON::PP` is part of core Perl) and any CGI-capable web server.

The intended operator is the institution's IT staff or a bioinformatician,
who installs OncoUnify once and schedules the loaders.  End users query the
database through the web interface, SQL, R or Python.

---

## Option A — Docker (demo or single-host deployment)

```bash
docker compose up --build            # demo: builds a database from tests/data
# open http://localhost:8080/
```

For institutional data, mount the reports read-only and let the container
load them at start-up:

```bash
docker build -t oncounify .
docker run -d --name oncounify -p 8080:80 \
  -e ONCOUNIFY_LOAD=1 \
  -v /secure/reports:/reports:ro \        # foundation/ genminetop/ guardant/ panels/ case_metadata.csv
  -v /secure/oncounify-db:/data \          # panels.db lives here
  -v /secure/htpasswd:/etc/oncounify/htpasswd:ro \   # enables HTTP Basic authentication
  oncounify
```

Re-run the loaders later with
`docker exec oncounify /opt/oncounify/docker/load_reports.sh`.

The container serves plain HTTP.  Put it behind the institution's TLS
reverse proxy / single sign-on, and keep it inside the hospital network.

## Option B — manual installation on Debian/Ubuntu with Apache

### 1. Packages

```bash
sudo apt-get install -y python3 python3-pandas python3-openpyxl sqlite3 \
    apache2 libcgi-pm-perl libdbi-perl libdbd-sqlite3-perl
sudo a2enmod cgid
```

### 2. Files

```bash
sudo git clone https://github.com/mizomizo1/OncoUnify.git /opt/oncounify
sudo install -d /var/www/oncounify/cgi-bin /var/www/oncounify/panel /var/lib/oncounify
sudo install -m 0755 /opt/oncounify/*.cgi /var/www/oncounify/cgi-bin/
sudo install -m 0644 /opt/oncounify/search.html /var/www/oncounify/panel/
```

If the CGIs are not served under `/cgi-bin`, edit the single constant
`CGI_BASE` near the end of `search.html`.

### 3. Load the reports

```bash
cd /opt/oncounify
DB=/var/lib/oncounify/panels.db
python3 load_foundation.py  $DB /data/foundation  --case-metadata /data/case_metadata.csv
python3 load_genminetop.py  $DB /data/genminetop  --case-metadata /data/case_metadata.csv
python3 load_guardant.py    $DB /data/guardant    --case-metadata /data/case_metadata.csv
python3 load_panel_genes.py $DB panels/            # assay gene lists (see panels/README.md)
sudo chown -R www-data:www-data /var/lib/oncounify && sudo chmod 640 $DB
```

Loader behaviour you can rely on in scripts:

* a report is identified by `(panel_name, report_id)`; loading it again
  **replaces** it and all its child rows in one transaction (no duplicates;
  amended reports and edited sidecar entries are picked up);
* every run ends with a `[SUMMARY]` line on stderr;
* the exit status is `0` when every file loaded, `1` when any file failed
  (the failed files are listed), `2` when there was nothing to do;
* `--fail-fast` stops at the first failure; `--genome-build` overrides the
  recorded assembly; Guardant accepts `--filename-pattern`,
  `--panel-version` and `--include-uncalled` (see `--help`).

A nightly job, for example:

```cron
30 2 * * *  cd /opt/oncounify && ONCOUNIFY_DB=/var/lib/oncounify/panels.db ONCOUNIFY_REPORTS=/data ./docker/load_reports.sh >> /var/log/oncounify-load.log 2>&1 || mail -s "OncoUnify load failed" it@example.org < /dev/null
```

### 4. Apache

`/etc/apache2/sites-available/oncounify.conf` (adapt the paths):

```apache
<VirtualHost *:443>
    ServerName oncounify.example.org
    DocumentRoot /var/www/oncounify
    Alias       /panel/   /var/www/oncounify/panel/
    ScriptAlias /cgi-bin/ /var/www/oncounify/cgi-bin/
    RedirectMatch ^/$ /panel/search.html

    SetEnv ONCOUNIFY_DB       /var/lib/oncounify/panels.db
    SetEnv ONCOUNIFY_CGI_URL  /cgi-bin
    SetEnv ONCOUNIFY_HTML_URL /panel
    # SetEnv ONCOUNIFY_MAX_ROWS 100000

    <Location />
        AuthType Basic
        AuthName "OncoUnify"
        AuthUserFile /etc/apache2/oncounify.htpasswd
        Require valid-user
    </Location>
    <Directory /var/www/oncounify/cgi-bin>
        Options +ExecCGI
        AllowOverride None
    </Directory>

    SSLEngine on
    SSLCertificateFile    /etc/ssl/certs/oncounify.crt
    SSLCertificateKeyFile /etc/ssl/private/oncounify.key
</VirtualHost>
```

```bash
sudo htpasswd -c /etc/apache2/oncounify.htpasswd alice
sudo a2enmod ssl && sudo a2ensite oncounify && sudo systemctl reload apache2
```

### 5. Smoke test

```bash
curl -u alice https://oncounify.example.org/panel/search.html
curl -u alice "https://oncounify.example.org/cgi-bin/panel_search.cgi?gene=TP53&view_mode=patient&download=tsv"
```

## Security and governance notes

* Every CGI opens the database with `SQLITE_OPEN_READONLY`; all queries are
  parameterized and all output is HTML-escaped.
* Read-only access prevents modification but is **not** an audit trail.
  Authentication and per-user query logging must come from the web server or
  the institutional SSO in front of it (e.g. `mod_auth_openidc` with
  `LogFormat` including `%u` and the query string).
* A search with no filter returns the whole registry up to
  `ONCOUNIFY_MAX_ROWS` rows, identifiers included, and can be exported as TSV.
  Restrict access accordingly.
* Encrypt the volume holding `panels.db` at rest with operating-system
  facilities if required by local policy.

## Known limitations of the reference interface

* Autocompletion suggests values for the last comma-separated term of a
  field only.
* Gene search resolves previous HGNC symbols (`gene_symbol_map`) but not
  aliases.
* The disease filter is a substring match across the disease fields and
  cannot use an index; gene and protein filters are exact and indexed.

## Upgrading from OncoUnify 1.x

Version-2 loaders refuse a version-1 database.  Either rebuild the database
from the source reports (recommended: this also applies Guardant call
filtering and recovers GenMineTOP breakpoints), or upgrade in place:

```bash
python3 migrate_db.py /var/lib/oncounify/panels.db     # keeps panels.db.v1-backup-<timestamp>
```
