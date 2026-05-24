# OncoUnify — Writing a new loader

OncoUnify is designed so that supporting a new cancer-panel vendor never
requires touching the schema, the CGIs, or the search interface. The entire
vendor-specific logic lives in a single Python script, `load_<vendor>.py`.

This document explains the contract a new loader must satisfy.

---

## 1. The contract

A loader is any executable Python script that, given:

- a path to a SQLite database (created from `schema.sql`), and
- a path to a directory containing native vendor deliverables,

inserts one row into `cases` for every report and zero or more rows into
`variants` (and optionally `non_human_contents`) for every reported finding.

Loaders are independent processes. They share **only** the schema and the
unique key `(panel_name, report_id)`. Two loaders never block each other.

---

## 2. Required helper functions

Every loader in this repository exposes four helpers, in this order:

| Function | Responsibility |
|---|---|
| `init_db(conn, schema_path)` | Idempotently execute `schema.sql` on the connection. |
| `get_or_create_case(conn, case)` | Return an existing `case_id` for `(panel_name, report_id)`, or insert and return a new one. |
| `insert_variants(conn, case_id, variants)` | Bulk-insert variant rows. |
| `parse_<vendor>(path)` | Parse a single native deliverable and return `(case_dict, variants_list[, non_humans_list])`. |

You can copy the helpers verbatim from `load_foundation.py`; only
`parse_<vendor>()` is vendor-specific.

---

## 3. Canonical column mapping

The full canonical schema is documented in [SCHEMA.md](SCHEMA.md). A loader
must map vendor fields onto canonical columns wherever a sensible mapping
exists. Fields that have no counterpart in the canonical schema should be
stored verbatim in `cases.other_info` or `variants.extra` as a serialised
`key=value;...` string, **not** dropped silently.

The minimum required population per row:

### `cases`

| Column | Required | Notes |
|---|---|---|
| `panel_name` | yes | Stable vendor + assay identifier (e.g., `FoundationOne`, `GenMineTOP`, `Guardant360`) |
| `vendor` | yes | Manufacturer name |
| `report_id` | yes | Unique within `panel_name` |
| `patient_id` | recommended | Institutional MRN or pseudonym |
| `date` | recommended | YYYY-MM-DD |

### `variants`

| Column | Required | Notes |
|---|---|---|
| `case_id` | yes | From `get_or_create_case()` |
| `gene` | yes | Official HGNC symbol when possible |
| `variant_type` | yes | One of `short_variant`, `cnv`, `rearrangement`, `expression`, `biomarker` |
| `protein_effect` | recommended | Normalised: leading `p.` stripped |
| `functional_effect` | recommended | One of `missense`, `nonsense`, `frameshift`, `splice`, `synonymous` |

---

## 4. Step-by-step recipe

The fastest way to start a new loader is to copy one of the existing files
and rewire the parser:

```bash
cp load_guardant.py load_yourvendor.py
```

Then:

1. Rewrite `parse_<vendor>(path)` to return the canonical `(case, variants)`
   pair. Keep everything else identical.
2. Replace the file-iteration helper if the vendor's deliverable is not a
   flat directory of files (e.g., a single ZIP containing one PDF per case).
3. Add a test fixture under `tests/data/<vendor>/` and verify the loader
   produces non-empty `cases` and `variants` rows when run against it.
4. Add a row to the vendor-to-canonical mapping table in the manuscript /
   `SCHEMA.md` so reviewers and downstream users can see exactly which
   vendor field populates which canonical column.

---

## 5. Idempotency and re-ingestion

Re-running a loader on the same input directory must not duplicate cases.
The existing loaders rely on the `UNIQUE(panel_name, report_id)` constraint
in the `cases` table; `get_or_create_case()` short-circuits when the pair
is already present.

If you need to **replace** an existing case (for example, because a vendor
has issued an amended report), delete the case first:

```sql
DELETE FROM variants            WHERE case_id IN (SELECT case_id FROM cases WHERE panel_name='Foo' AND report_id='12345');
DELETE FROM non_human_contents  WHERE case_id IN (SELECT case_id FROM cases WHERE panel_name='Foo' AND report_id='12345');
DELETE FROM cases               WHERE panel_name='Foo' AND report_id='12345';
```

Then re-run the loader.

---

## 6. Conventions

- **Stdout is silent**, stderr carries human-readable progress messages.
- Loaders should be runnable as `python3 load_<vendor>.py <db> <dir>` with
  no other arguments. Extra options (e.g., a panel-name override) are fine,
  but the two positional arguments must stay.
- Loaders must **not** modify `journal_mode`. The CGI side assumes
  `journal_mode=DELETE`.
- Numeric coercion goes through `int_or_none()` / `float_or_none()`. Empty
  strings, whitespace, and `"N/A"`-style sentinels must become SQL `NULL`,
  not `0`.

---

## 7. Submitting a new loader upstream

When contributing a loader back to the project:

1. Open a pull request that touches **only** `load_<vendor>.py`,
   `tests/data/<vendor>/`, and `docs/SCHEMA.md` (the mapping table).
2. Include at least one synthetic sample report under `tests/data/<vendor>/`
   with all identifying information removed.
3. Update the README's vendor list.

Contributions are reviewed for: canonical mapping correctness, idempotency,
and absence of vendor-proprietary information in the test fixtures.
