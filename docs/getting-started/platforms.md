# Supported platforms

BSReadSim's production support contract is intentionally narrower than the
set of systems on which the source might compile.

| Environment | Status | Notes |
| --- | --- | --- |
| Linux x86-64 | Supported | Tested in CI with Python 3.10, 3.11, 3.12, 3.13, and 3.14 |
| Windows through WSL2 | Supported | Install and run entirely inside the Linux distribution |
| Native Windows | Not supported | The pinned HTSlib build requires a POSIX Autotools and Make environment |
| macOS | Not currently supported | No release CI or compatibility guarantee |
| Linux ARM64 | Not currently supported | No release CI or compatibility guarantee |

## Build environment

A supported source build needs:

- Python 3.10 or newer;
- CMake 3.20 or newer;
- a C++17 compiler and a Python C-extension toolchain;
- zlib development headers;
- Autoconf and Make; and
- Git with recursive submodule support.

The Python package builds and bundles its pinned `htsim-core` executable,
HTSlib, and private C extension. It does not use a separately installed
system HTSlib.

## WSL2

Clone the repository into the WSL filesystem when practical, install the Linux
build dependencies in that distribution, and run Python and BSReadSim there.
Do not mix a native Windows Python interpreter with Linux build products.

Python installation from source requires the complete toolchain listed above.
The planned Bioconda package will supply its compiled components through Conda
and will not depend on a PyPI wheel.
