# Synthetic vendor sample reports

This directory holds **fully synthetic** sample reports used to smoke-test
the three loaders bundled with OncoUnify.

No real patient data is present. Identifiers, dates, and variants have been
generated from scratch and bear no relationship to any real patient,
laboratory, or institution.

## Layout

```
tests/data/
├── foundation/      synthetic FoundationOne / FoundationOneLiquid XML
├── genminetop/      synthetic GenMineTOP XML
└── guardant/        synthetic Guardant360 Excel reports
```

## Regenerating

The fixtures are static. If you add a new field to the canonical schema or
to a loader's parser, append a new synthetic case to the relevant directory
rather than editing existing fixtures, so that the existing tests remain
reproducible.

## Smoke test

```bash
cd ..
python3 ../load_foundation.py  /tmp/test.db tests/data/foundation/
python3 ../load_genminetop.py  /tmp/test.db tests/data/genminetop/
python3 ../load_guardant.py    /tmp/test.db tests/data/guardant/
sqlite3 /tmp/test.db "SELECT panel_name, COUNT(*) FROM cases GROUP BY panel_name;"
```

Each loader must add at least one row to `cases` and ≥1 row to `variants`.
