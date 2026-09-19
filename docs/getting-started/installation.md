# Installation

## Supported platforms

| Platform | Environment |
| --- | --- |
| Linux x86-64 | Python 3.10 or later |
| Windows (WSL2) | Linux Python 3.10 or later inside WSL2 |

BSReadSim supports Linux x86-64 and Windows through WSL2. Native Windows,
macOS, and Linux ARM64 are not currently supported.

## Install from Bioconda

!!! info "Coming soon"
    Bioconda installation is coming soon. It's the recommended installation 
    method, especially on shared systems and HPC clusters.

## Install from source

### Requirements

The following dependencies are needed to build BSReadSim from source.

- Python 3.10 or later
- GCC or Clang with C++17 support
- CMake 3.20 or later
- GNU Make, Autoconf, and zlib

Choose one of the following ways to prepare an isolated build environment.

??? tip highlight "Using Conda or Mamba"
    Conda can provide the compilers, libraries, and Python environment in user
    space without administrator access:

    ```bash
    conda create -n bsenv -c conda-forge \
      python pip cmake c-compiler cxx-compiler make autoconf zlib git
    conda activate bsenv
    ```

    `mamba` can be used in place of `conda`.

??? tip highlight "Using venv with system tools"
    `venv` isolates the Python packages; the build tools must be installed on
    the system or made available separately.

    On Ubuntu or WSL2 with sudo access, install the required tools with:

    ```bash
    sudo apt-get update
    sudo apt-get install --yes \
      autoconf build-essential cmake git python3-dev python3-venv zlib1g-dev
    ```

    On HPC or shared server without sudo access, load or install the dependencies in your user 
    environment first, then create and activate the virtual environment:

    ```bash
    python3 -m venv bsenv
    source bsenv/bin/activate
    ```

    If the virtual environment was created with HPC modules, load the same 
    modules before activating it in future sessions or batch jobs.

### Build and install

With the selected environment active, clone the repository and install BSReadSim:

```bash
# For Conda: conda activate bsenv
# For venv: source bsenv/bin/activate

git clone --recurse-submodules https://github.com/wbvguo/BSReadSim.git
cd BSReadSim
python -m pip install .
```

## Verify the installation

Check the installed version with:

```bash
bsreadsim --version
```

<div class="next-step" markdown>

**Next:** [Run the quick start](quickstart.md) to create a demo WGBS simulation.

</div>
