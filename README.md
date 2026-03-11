# FlamMap

`FlamMap` is a flammability workflow package to assess the flammability of chemical species using calculated adiabatic flame temperature (CAFT).

It is designed around the public interface:

- Input:
  - one floating-point `tae_hartree`
  - one `orca.out` file for the frequency calculation
- Output:
  - one flammability phase-diagram PDF
  - one `.dat` file containing the CAFT grid
  - one YAML file containing the species thermochemistry and metadata
  - one JSON summary containing `tae_hartree`, `orca_out`, `Hf`, `LFL`, and `UFL`

This repository is organized as an independent package. It provides:

- a packaged Python CLI
- a config file for a single species run
- one documented pipeline entry point
- a fixed workflow that generates the YAML, `.dat`, diagram PDF, and 1600 K flammability limits in one run

## Layout

- `pyproject.toml`: Python package metadata
- `configs/example.yaml`: example configuration
- `data/reference/elem_enthalpies.json`: vendored elemental enthalpy corrections
- `examples/orca/orca.out`: bundled example ORCA frequency output
- `src/flammability/`: package source

## Quick Start

1. Create a config file by copying `configs/example.yaml`.
2. Update `tae_hartree`, `orca_out`, and `output_dir`.
   The repository already includes one example ORCA output at `examples/orca/orca.out`.
   Paths in the config can be written relative to the `FlamMap` repository root.
3. Run the pipeline:

```bash
python FlamMap/run.py FlamMap/configs/example.yaml
```

To validate the config only:

```bash
python FlamMap/run.py --validate-only FlamMap/configs/example.yaml
```

## Current Scope

The package contains:

- internal Python modules for thermochemistry-to-YAML generation
- internal Python modules for sequential CAFT grid computation
- vendored reference data required for thermo generation

## Outputs

By default the pipeline writes:

- `<species>.yaml`
- `<species>.dat`
- `<species>.pdf`
- `<species>.json`

The workflow uses fixed values:

- `c_bond = 718.1`
- CAFT flammability threshold: `1600 K`
