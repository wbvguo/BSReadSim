# Test fixtures

This directory is intentionally flat. `fixtures.json` records the format,
purpose, validity, and expected use of every fixture. Add metadata there in the
same change whenever a fixture is added, renamed, or removed.

The hexadecimal files contain lowercase hexadecimal representations of frozen
binary frames. Whitespace is not part of their decoded bytes.

`tests/unit/test_fixture_registry.py` enforces one metadata entry per file and
requires every invalid fixture to declare its expected rejection result.
