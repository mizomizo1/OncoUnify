# OncoUnify — Installation and Deployment

This document describes how to deploy OncoUnify on a single Linux host.
The reference stack is:

- **Python 3.10+** (loaders)
- **Perl 5.20+** (CGI scripts)
- **SQLite 3.30+**
- **Apache 2.4** with `mod_cgi` (or `mod_cgid` + `suEXEC`)

OncoUnify has no external service dependencies; everything runs inside the
host filesystem. A single `panels.db` file is the only persistent store.

---

## 1. System prerequisites

### 1.1 Operating system packages (Ubuntu/Debian)

```bash
sudo apt-get update
sudo apt-get install -y \
    python3 python3-pip python3-venv \
    perl libcgi-pm-perl libdbi-perl libdbd-sqlite3-perl libjson-perl \
    sqlite3 \
    apache2 apache2-utils
sudo a2enmod cgi rewrite
```

### 1.2 Python dependencies

The loaders use only the Python standard library plus:

```bash
pip3 install --user openpyxl pandas
```

These are required by `load_guardant.py` to read `.xlsx` files. The
FoundationOne and GenMineTOP loaders depend only on the standard library
(`xml.etree.ElementTree`, `sqlite3`).

---

## 2. Filesystem layout (reference)

```
/var/www/oncounify/
├── panel/               (DocumentRoot)
│   └── search.html
└── cgi-bin/
    ├── panel_search.cgi
    ├── case_detail.cgi
    ├── panel_stats.cgi
    ├── suggest.cgi
    └── logout.cgi

/var/lib/oncounify/
├── panels.db            (shared SQLite database)
└── schema.sql           (canonical schema; bundled in the repo)
```

The database path is read by every CGI from a single environment variable
`ONCOUNIFY_DB`. Set it in your Apache configuration (see §4).

---

## 3. Database initialisation and ingestion

```bash
# 3.1 Initialise an empty database from the canonical schema
sqlite3 /var/lib/oncounify/panels.db < /var/lib/oncounify/schema.sql

# 3.2 Ingest vendor reports (any order; loaders are idempotent on report_id)
python3 load_foundation.py  /var/lib/oncounify/panels.db /data/foundation/
python3 load_genminetop.py  /var/lib/oncounify/panels.db /data/genminetop/
python3 load_guardant.py    /var/lib/oncounify/panels.db /data/guardant/

# 3.3 Optional: tighten file permissions
sudo chown -R www-data:www-data /var/lib/oncounify
sudo chmod 640 /var/lib/oncounify/panels.db
```

A loader can be re-run safely. Each `(panel_name, report_id)` pair is unique
in the `cases` table; existing cases are skipped rather than duplicated.

---

## 4. Apache configuration (example)

`/etc/apache2/sites-available/oncounify.conf`:

```apache
<VirtualHost *:443>
    ServerName  oncounify.example.org

    DocumentRoot /var/www/oncounify/panel
    Alias /panel/ /var/www/oncounify/panel/
    ScriptAlias /cgi-bin/ /var/www/oncounify/cgi-bin/

    # Database location is read by every CGI from this env var
    SetEnv ONCOUNIFY_DB /var/lib/oncounify/panels.db

    <Directory "/var/www/oncounify/cgi-bin">
        Options +ExecCGI
        AddHandler cgi-script .cgi
        AllowOverride None
        Require valid-user
        AuthType Basic
        AuthName "OncoUnify"
        AuthUserFile /etc/apache2/oncounify.htpasswd
    </Directory>

    <Directory "/var/www/oncounify/panel">
        AllowOverride None
        Require valid-user
        AuthType Basic
        AuthName "OncoUnify"
        AuthUserFile /etc/apache2/oncounify.htpasswd
    </Directory>

    SSLEngine on
    SSLCertificateFile      /etc/ssl/certs/oncounify.crt
    SSLCertificateKeyFile   /etc/ssl/private/oncounify.key
</VirtualHost>
```

Then:

```bash
sudo htpasswd -c /etc/apache2/oncounify.htpasswd alice
sudo a2ensite oncounify
sudo systemctl reload apache2
```

---

## 5. Smoke test

```bash
curl -u alice https://oncounify.example.org/panel/search.html
curl -u alice "https://oncounify.example.org/cgi-bin/panel_stats.cgi"
curl -u alice "https://oncounify.example.org/cgi-bin/panel_search.cgi?gene=TP53&view=patient"
```

Each request must return HTTP 200 with valid HTML. If the page reports
`DB error`, verify that `ONCOUNIFY_DB` is exported to the CGI environment
and that the file is readable by the Apache user.

---

## 6. Upgrading the schema

OncoUnify uses `CREATE TABLE IF NOT EXISTS`, so re-running the schema file
on an existing DB does **not** add new columns. When the schema gains new
columns:

```bash
sqlite3 /var/lib/oncounify/panels.db < migrations/2026-01-add-column.sql
```

The `migrations/` directory under the repository contains incremental
`ALTER TABLE` statements; apply them in chronological order.
