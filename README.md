BSReadSim: a versatile and efficient simulator to generate realistic bisulfite sequencing reads

Utility under construction...

## Build (0.1.0)

htsim needs a C/C++ compiler, GSL, zlib, and the libbz2, liblzma, and libcurl
development files used by the bundled HTSLIB 1.10.2:

    make -C HTSLIB lib-static lib-shared
    make -C HTSIM

`HTSIM/htsim` simulates WGS reads, or with `-O 1` the reads of WGBS, RRBS
(`-T 1`, with a BED of restriction fragments from `HTSIM/rrcut`), and targeted
bisulfite sequencing (`-T 2`) for the Python modules. Run either tool without
arguments to list its options.

Known limitations of 0.1.0:

- The Python modules in `Simulate/` still use the htsim options of March 2023;
  they are not connected to this htsim.
- Fragment positions are drawn from generators seeded by `std::random_device`,
  so runs are not reproducible.
- `-V` (an input snpmeth file) is not read yet.
