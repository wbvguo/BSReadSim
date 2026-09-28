# Python and integration tests

Everything here is a maintained correctness contract for the current product.

- `semantic/` defines the scientific behavior of simulated data. It drives
  only the public CLI and reads BAM truth tags and exported truth files, so it
  must keep passing, unchanged, across any reimplementation of the internals.
  Its statistical checks use tolerances of at least four standard deviations;
  `BSREADSIM_SEMANTIC_SEED` reruns every simulation with another seed.
- `unit/` contains Python unit tests.
- `integration/` contains executable and cross-language checks.
- `helpers/` contains support code shared by tests.
- `fixtures/` contains small, frozen inputs described by `fixtures.json`.

The native core has no separate unit tests: the semantic suite checks its
behavior through the public interface.

Exploratory tests and large source data belong in the ignored `dev/`
directory, not here.

## Running the tests

The root CMake build is the single test entry point. It builds an isolated
Python package and the htsim core and runs the Python unit, semantic, and
cross-language integration suites:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Individual test modules assume the package under test is already importable.
Some cross-language integration scripts also add the repository root to
`sys.path` so they can import shared helpers from `tests/`; the package itself
still comes from the isolated `PYTHONPATH` configured by CMake/CTest.
