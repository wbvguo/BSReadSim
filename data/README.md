# Product data

This directory is the canonical source for immutable data released with
BSReadSim: bundled models, supported profiles, and small user-facing examples.

Use `models/` for released pretrained models, `profiles/` for supported
runtime parameter sets, and `example/` for inputs used by user
documentation. Register every bundled resource in `registry.json` with its
checksum, provenance, and license. The containing BSReadSim release and exact
SHA-256 identify the resource bytes.

External research datasets and local benchmark corpora do not belong here;
they live in the ignored `workspace/datasets/` directory.

The source distribution includes this directory directly. Wheel builds copy
the registry and each registered resource into `bsreadsim/data/`; the top-level
file remains the only source copy. Runtime code accesses installed resources
through `importlib.resources` rather than assuming a site-packages path.

Keep large optional models outside the wheel. They may still use registry
metadata in a future download/cache workflow, but must not be added as bundled
resources unless package-size and redistribution requirements have been
reviewed.
