# FlamMap

`FlamMap` is a flammability workflow package that computes flammability limits and phase diagrams for chemical species using calculated adiabatic flame temperature (CAFT).

The default workflow takes a single XYZ geometry, predicts the total atomization energy and vibrational frequencies via ML models (skala + hip), and produces:

- a flammability phase-diagram PDF
- a `.dat` file containing the CAFT grid
- a YAML file containing the species thermochemistry
- a JSON summary including TAE, frequencies, Hf, LFL, and UFL

## Installation

This project uses `uv` for environment and dependency management. Python 3.11 is recommended.

```bash
uv python install 3.11
uv venv --python 3.11
source .venv/bin/activate
uv sync
```

`hip` is installed as a path dependency from `/home/jiaqi/git/hip`.
`skala` is not packaged; its source dir is loaded via `sys.path` at runtime.

## Usage

### Default: ML mode

Pass an XYZ file. Everything else is auto-resolved from defaults.

```bash
uv run python run.py examples/molecule.xyz
```

ML mode requires:

- `HIP_CKPT` env var or `--hip-checkpoint PATH` (default: `/home/jiaqi/git/hip/ckpt/hip_v2.ckpt`)
- `SKALA_TAE_DIR` env var or `--skala-dir PATH` (default: `/home/jiaqi/git/Skala_TAE`)

### AB mode

Provide TAE and frequencies via a small YAML alongside the XYZ:

```bash
uv run python run.py examples/molecule.xyz --mode ab --yaml examples/molecule_data.yaml
```

The case-data YAML contains exactly two fields:

```yaml
tae:   0.669583                                    # total atomization energy in Hartree
freqs: [1340.40, 1340.57, 1340.65, 1557.35, 1557.45, 3029.19, 3131.50, 3131.75, 3131.78]  # vibrational frequencies in cm^-1
```

All frequencies must be positive. The count must match the geometry: `3N-5` for linear molecules, `3N-6` for nonlinear.

### Other flags

- `--validate-only` — check inputs and exit
- `--output-dir DIR` — override `outputs/`
- `--no-plot` — skip the phase-diagram PDF
- `--hip-device cpu|cuda` — torch device for hip (default `cpu`)

## Fixed workflow constants

- `c_bond = 718.1`
- CAFT flammability threshold: `1600 K`

## Layout

- `src/flammability/` — package source
- `data/reference/elem_enthalpies.json` — vendored elemental enthalpy corrections
- `data/cantera/` — vendored Cantera reference and product mechanisms
- `examples/molecule.xyz` — sample CH4 geometry
- `examples/molecule_data.yaml` — sample ab-mode case data for CH4
