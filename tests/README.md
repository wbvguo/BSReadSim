# Python and integration tests

Everything here is a maintained correctness contract for the current product.

- `unit/` contains Python unit tests.
- `integration/` contains executable and cross-language checks.
- `helpers/` contains support code shared by tests.
- `fixtures/` contains small, frozen inputs described by `fixtures.json`.

C++ unit tests owned by the native core live under `htsim/tests/` and use the
same root CMake/CTest entry point as these suites.

Human-maintained exploratory tests and large source data belong under `dev/`.
Agent experiments and their generated data belong in a dated Git worktree
under `workspace/worktree/`, not in this directory.

## Running the tests

The root CMake build is the single test entry point. It builds an isolated
Python package with the C extension and runs the C++, Python unit, and
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
