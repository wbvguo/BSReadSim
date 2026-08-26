# Installation

BSReadSim supports Linux x86-64 with Python 3.10 or newer, including Windows through WSL2.
See [supported platforms](platforms.md) for details.

## Install from source

### Requirements

- CMake 3.20 or newer
- A C++17 compiler and Python development headers
- Git, Autoconf, Make, and zlib development files

On Ubuntu or WSL2, install the required build tools with:

```bash
sudo apt-get update
sudo apt-get install --yes \
  autoconf build-essential cmake git python3-dev zlib1g-dev
```

### Build and install

Clone the repository with its pinned submodules and install the package:

```bash
git clone --recurse-submodules --shallow-submodules \
  https://github.com/wbvguo/BSReadSim.git
cd BSReadSim
python3 -m pip install .
```

HTSlib and htscodecs are built from pinned submodules; a separate system HTSlib installation is not required.

??? tip "Optional installation in a virtual environment"
    On Ubuntu or WSL2, create and activate a virtual environment before installing:

    ```bash
    sudo apt-get install --yes python3-venv
    python3 -m venv .venv
    source .venv/bin/activate
    python -m pip install .
    ```

    Activate it again with `source .venv/bin/activate` when opening a new shell.

??? tip "Installing from an existing checkout"
    Initialize any missing submodules before installing:

    ```bash
    git submodule update --init --recursive --depth 1
    python3 -m pip install .
    ```

## Bioconda package (coming soon)

Installation instructions via Bioconda will be added here when it becomes available.

## Verify the installation

Both commands should complete without errors:

```bash
bsreadsim --version
bsreadsim export --help
```

<div class="next-step" markdown>

**Next:** [Run the quick start](quickstart.md) to create and verify a small WGBS
simulation.

</div>
