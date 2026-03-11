# FlamMap

`FlamMap` is a flammability workflow package to assess the flammability of chemical species using calculated adiabatic flame temperature (CAFT).

It is designed around the public interface:

- Input:
  - one floating-point `tae_hartree`
  - one geometry source: `orca.out` or `xyz_geom`
  - either vibrational frequencies parsed from `orca.out` or a direct `freqs_cm1` array in `cm^-1`
- Output:
  - one flammability phase-diagram PDF
  - one `.dat` file containing the CAFT grid
  - one YAML file containing the species thermochemistry and metadata
  - one JSON summary containing `tae_hartree`, `orca_out`, `xyz_geom`, `freqs_cm1`, `Hf`, `LFL`, and `UFL`

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
- `examples/`: example inputs
- `src/flammability/`: package source

## Installation

This project uses `uv` for environment and dependency management.
Recommended Python version: `3.11`.

1. Create a virtual environment and install dependencies:

```bash
uv python install 3.11
uv venv --python 3.11
source .venv/bin/activate
uv sync
```

2. Run the CLI:

```bash
uv run python run.py configs/example.yaml
```

`pymsym` is installed as part of `uv sync` and is used to determine the rotational symmetry number directly from the ORCA geometry.
Python `3.11` is recommended for compatibility with the scientific Python stack used here.

## Quick Start

1. Create a config file by copying `configs/example.yaml`.
2. Update `tae_hartree`, geometry input, and `output_dir`.
   The repository already includes one example ORCA output at `examples/orca/orca.out`.
   If you want to override the vibrational frequencies, also set `freqs_cm1` to a list of frequencies in `cm^-1`.
   If `orca_out` is omitted, you must provide `xyz_geom` for geometry, formula inference, and symmetry-number detection.
   Paths in the config can be written relative to the `FlamMap` repository root.
3. Run the pipeline:

```bash
uv run python run.py configs/example.yaml
```

To validate the config only:

```bash
uv run python run.py --validate-only configs/example.yaml
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

## Frequency Input Options

- Default: omit `freqs_cm1` and the code will parse vibrational frequencies from `orca.out`.
- Override: set `inputs.freqs_cm1` to a YAML list of positive frequencies in `cm^-1`.
- The number of positive vibrational frequencies must match the geometry: `3N-5` for linear molecules or `3N-6` for nonlinear molecules.

## Geometry Input Rules

- `inputs.orca_out` can provide both geometry and vibrational frequencies.
- `inputs.xyz_geom` can provide geometry for formula inference, inertia, and symmetry-number detection.
- If `inputs.freqs_cm1` is provided and `inputs.orca_out` is omitted, `inputs.xyz_geom` is required.
- If `inputs.freqs_cm1` is omitted, `inputs.orca_out` is required because frequencies must be parsed from ORCA output.
