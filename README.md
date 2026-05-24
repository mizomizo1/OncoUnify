# OncoUnify

A deployable, vendor-agnostic database and web interface for the integrated
re-use of multi-vendor cancer genomic panel reports.

OncoUnify ingests reports produced by **FoundationOne / FoundationOneLiquid**
(Foundation Medicine, XML), **GenMineTOP** (XML), and
**Guardant360** (Excel `.xlsx`), normalises them into a single relational
schema, and exposes them through a lightweight CGI-based web interface for
cross-vendor search, per-case drill-down, and registry-level statistics.

The package is intentionally minimal — Python 3, Perl 5, SQLite, and Apache
(or any other CGI-capable web server) are the only runtime dependencies — so
that each institution can self-host an instance with the data kept on
premises.

## Quick start

```bash
# 1. Initialise the database and ingest reports
python3 load_foundation.py  panels.db /path/to/foundation_xml_dir
python3 load_genminetop.py  panels.db /path/to/genminetop_xml_dir
python3 load_guardant.py    panels.db /path/to/guardant_xlsx_dir

# 2. Deploy the CGIs and search.html behind a web server
#    (see docs/INSTALL.md for an Apache + suEXEC example)
```

## Repository layout

```
OncoUnify/
├── schema.sql                  # canonical relational schema
├── load_foundation.py          # FoundationOne / FoundationOneLiquid XML loader
├── load_genminetop.py          # GenMineTOP XML loader
├── load_guardant.py            # Guardant360 Excel loader
├── panel_search.cgi            # cross-vendor search
├── case_detail.cgi             # per-case drill-down
├── panel_stats.cgi             # registry-level statistics
├── suggest.cgi                 # autocomplete endpoint
├── logout.cgi                  # Basic-auth sign-out helper
├── search.html                 # static search form
├── docs/
│   ├── INSTALL.md              # deployment guide
│   ├── LOADER_DEV.md           # how to add a new loader
│   └── SCHEMA.md               # canonical schema reference
└── tests/
    └── data/                   # synthetic vendor sample reports
```

## Extending OncoUnify

Adding support for a new panel vendor is intentionally local: implement a
new `load_<vendor>.py` script that parses the vendor's native deliverable
and inserts rows into `cases`, `variants`, and (optionally)
`non_human_contents`. No CGI or schema change is required as long as the
vendor's fields can be mapped onto the canonical columns documented in
[docs/SCHEMA.md](docs/SCHEMA.md). See [docs/LOADER_DEV.md](docs/LOADER_DEV.md)
for a step-by-step recipe.

## License

OncoUnify is released under the MIT License.

## Citation

If you use OncoUnify in your research, please cite the accompanying paper:

> OncoUnify: a deployable, vendor-agnostic database and web interface for
> the integrated re-use of multi-vendor cancer genomic panel reports.
> *Database (Oxford)*, submitted.
