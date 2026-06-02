# FlamMap

`FlamMap` is a flammability workflow package that computes flammability limits and phase diagrams for chemical species using calculated adiabatic flame temperature (CAFT).

The default workflow takes a single XYZ geometry, predicts the total atomization energy and vibrational frequencies via ML models (skala + hip), and produces:

- a flammability phase-diagram PDF
- a `.dat` file containing the CAFT grid
- a YAML file containing the species thermochemistry
- a JSON summary including TAE, frequencies, Hf, the CAFT threshold, LFL, and UFL

## Installation

This project uses `uv` for environment and dependency management. Python 3.12 is required (pinned by `hip`).

`hip` and `skala` are pinned to specific upstream commits in [pyproject.toml](pyproject.toml) and pulled by `uv sync` — no sibling clones required. Bump the `rev` SHAs under `[tool.uv.sources]` to track newer upstream code.

GPU install (CUDA 12.1 wheels, suitable for NVIDIA drivers ≥ 530):

```bash
git clone <flammap-url>
cd FlamMap
uv python install 3.12
uv sync
```

`uv sync` pulls the cu121 PyTorch wheels (`torch==2.5.1+cu121`) plus the matching `pyg-lib`, `torch-scatter`, `torch-sparse`, `torch-cluster` GPU wheels from the PyG wheel index, `cupy-cuda12x` and `gpu4pyscf-cuda12x` for the skala GPU SCF, and clones+builds `hip` and `skala` from GitHub at the pinned commits.

For local development against a live checkout, edit `[tool.uv.sources]` in [pyproject.toml](pyproject.toml) to swap `{ git = ..., rev = ... }` for `{ path = "../hip", editable = true }` (and likewise for skala).

Verify the GPU stack:

```bash
uv run python -c "import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))"
```

> **Driver note.** The pin is cu121 because CUDA 12.6 wheels require driver ≥ 560. If you have a newer driver and want cu126, change the `pytorch-cu121` index URL to `https://download.pytorch.org/whl/cu126`, bump the find-links URL to `torch-2.11.0+cu126`, and use `hip[cuda126]` instead of `hip[cuda121]` in `pyproject.toml`.

`compute_tae` calls the installed `skala` package directly: `skala.gpu4pyscf.SkalaKS` when `--skala-device cuda` (end-to-end gpu4pyscf), `skala.pyscf.SkalaKS` when `--skala-device cpu`. Per-element reference SCFs (H/C/N/O) are cached to `~/.cache/flammap/atomic_refs_<xc>_<basis>.json` (override with `FLAMMAP_CACHE`) so they run once per (xc, basis) pair.

To run the full pipeline on GPU:

```bash
uv run python run.py examples/molecule.xyz --skala-device cuda --hip-device cuda
```

## Usage

### Default: ML mode

Pass an XYZ file. Everything else is auto-resolved from defaults.

```bash
uv run python run.py examples/molecule.xyz
```

ML mode needs the hip checkpoint and the installed `skala` package. Both are resolved automatically:

- hip checkpoint — if `--hip-checkpoint PATH` and `HIP_CKPT` are both unset, `hip_v3_cf.ckpt` (the conservative-forces variant, trained with `model.direct_forces=False` — appropriate for Hessian-based frequency analysis) is downloaded from [`huggingface.co/andreasburger/hip`](https://huggingface.co/andreasburger/hip) on first use and cached under `$HF_HOME` (default `~/.cache/huggingface/hub/`). Subsequent runs hit the cache. Set `--hip-checkpoint` / `HIP_CKPT` to point at a local file (e.g. `ckpt/hip_v3.ckpt`) to skip the download or use a different variant.
- `skala` — provided by `uv sync`. Functional weights (`skala-1.0`, `skala-1.1`) are themselves pulled from HuggingFace on first use by the `skala` package.

### AB mode

"AB" stands for *ab initio*: this mode is the escape hatch from the ML predictions. Compute TAE and harmonic frequencies yourself with a high-level ab initio method (e.g. CCSD(T)/CBS, MP2, or a hybrid DFT), drop them into a small YAML, and the rest of the pipeline (thermochemistry → CAFT grid → LFL/UFL) runs unchanged. Use it when you want reference numbers instead of skala+hip, or to compare ML against ab initio on the same downstream code path.

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
- `--threshold K` — CAFT cutoff temperature (K) for LFL/UFL extraction (default `1600`). Recorded per run in the JSON summary and the `threshold_K` column of `SUMMARY.csv`. Applies to single, `--input-list` batch, and HyperQueue modes alike.
- `--no-plot` — skip the phase-diagram PDF
- `--skip-existing` — exit early if `<output_dir>/json/<xyz-stem>.json` already exists; works in both single-molecule and batch mode (HyperQueue retry guard)
- `--hip-device auto|cpu|cuda` — torch device for hip (default `auto`)
- `--skala-device auto|cpu|cuda` — torch device for the skala SCF (default `auto`)

`auto` reads the XYZ and picks `cuda` if the molecule has more than 10 heavy (non-H) atoms and CUDA is actually available; otherwise `cpu`. Small molecules go faster on CPU because of the ~10-15 s cupy/gpu4pyscf init cost. Override with `--skala-device cuda` etc. if you want a specific device.

Per-molecule outputs use the same layout in single and batch modes (`yaml/<stem>.yaml`, `dat/<stem>.dat`, `pdf/<stem>.pdf`, `json/<stem>.json`) — see the tree in [Batch mode](#batch-mode) below. `SUMMARY.csv` / `FAILED.csv` are batch-only.

### Batch mode

Two paths, pick by scale: `--input-list` (below) for single-node runs up to a few hundred molecules — one Python process per worker amortizes the ~5–15 s torch/hip/skala import once across many molecules. For multi-node runs or thousands of molecules, jump to [HyperQueue](#multi-node-batches-via-hyperqueue) — dynamic load balancing across nodes outweighs the per-task interpreter startup cost.

```bash
ls molecules/*.xyz > inputs.txt
uv run python run.py --input-list inputs.txt --jobs 12 --no-plot --skip-existing
```

- `--input-list FILE` — text file with one XYZ path per line. Blank lines and lines starting with `#` are ignored.
- `--jobs N` — parallel worker processes. Workers use `multiprocessing` with `spawn`, so CUDA contexts are safe. Each worker is automatically capped at `nproc / jobs` OMP threads to avoid oversubscription. For a 96-core box, `--jobs 12` gives each worker 8 threads, which is the sweet spot for small CHON molecules.
- `--skip-existing` — skip XYZ files whose output JSON already exists in `--output-dir`. Use this to resume after a crash.
- `--no-plot` — recommended for big batches (~5 s/molecule saved).

Template SLURM script: [job.sh.example](job.sh.example). Copy next to your data, edit `FLAMMAP_DIR` / `INPUT_LIST` / `OUTPUT_DIR` / `JOBS`, `sbatch`.

Per-molecule outputs land in `--output-dir`, sorted into per-type subdirs and named by the **XYZ filename stem** (not the chemical formula — avoids collisions between isomers):

```
<output_dir>/
├── yaml/<name>.yaml      # Cantera mechanism with the species added
├── dat/<name>.dat        # CAFT temperature grid
├── pdf/<name>.pdf        # phase-diagram plot (only when --plot)
├── json/<name>.json      # summary (TAE, freqs, Hf, threshold_K, LFL, UFL, ...)
├── SUMMARY.csv           # batch only: one row per success
└── FAILED.csv            # batch only: one row per failure
```

- `SUMMARY.csv` — columns: `stem, formula, tae_Ha, Hf_298K_kJ, threshold_K, LFL_percent, UFL_percent, n_freqs, elapsed_s`. Appended across runs; `threshold_K` records the CAFT cutoff each row's LFL/UFL was extracted at.
- `FAILED.csv` — columns: `stem, formula, error_type, error, elapsed_s`. Appended across runs. Re-run a single offending molecule with `run.py <xyz>` to get a full traceback for debugging.

`stem` is the trailing `_`-separated piece of the XYZ filename (e.g. `C2H2_ca2cc2cc.xyz` → `ca2cc2cc`) — the unique-hash convention used by most isomer datasets. Per-molecule output files still use the full filename stem.

Failures (parse errors, SCF non-convergence, imaginary modes …) are logged and the loop continues; one bad molecule never kills the batch.

### Multi-node batches via HyperQueue

For batches that need more than one node (e.g. 12k molecules across 12×128-core nodes), use [HyperQueue](https://it4innovations.github.io/hyperqueue/stable/) instead of `--jobs`. HQ schedules one task per molecule across all allocated nodes dynamically, so load balance is naturally better than static SLURM array sharding.

Workflow:

1. Allocate nodes via SLURM, start HQ server on the head node, launch one HQ worker per node with `srun --overlap`.
2. Submit per-molecule tasks: `hq submit --each-line xyzlist.txt -- python run.py {entry} --skip-existing --no-plot`.
3. After HQ finishes, aggregate the per-molecule JSONs into the standard summary CSV: `python run.py --collect-summary xyzlist.txt --output-dir <dir>`.

Template script: [job_hq.sh.example](job_hq.sh.example). Copy next to your data, edit `FLAMMAP_DIR` / `INPUT_DIR` / `OUTPUT_DIR` / `CPUS_PER_TASK`, then run it on the login node (it submits tasks to a HQ server you've already started in a separate SLURM worker job).

The `--skip-existing` flag works in single-molecule mode too — that's the HQ retry guard. Re-submitting a failed HQ job re-runs only the molecules whose JSONs don't exist yet. `--collect-summary` fills `elapsed_s` and `threshold_K` from each per-molecule JSON (recorded by `run_pipeline`); failed/missing molecules have no JSON, so their `elapsed_s` stays blank — check `hq job info` for those. JSONs written before `--threshold` existed have no `threshold_K` field, so that column is left blank for them (those runs all used `1600 K`).

## Fixed workflow constants

- Elemental atomization enthalpies dH_f(X,g,298.15K), kJ/mol per atom, from CODATA Key Values for Thermodynamics (Cox, Wagman & Medvedev, 1989; redistributed by NIST WebBook): H = 217.998, C = 716.68, N = 472.68, O = 249.18. Stored in `data/reference/elem_enthalpies.json` with per-element citations.
- CAFT flammability threshold: defaults to `1600 K`; override per run with `--threshold K`.

## Layout

- `src/flammability/` — package source
- `data/reference/elem_enthalpies.json` — vendored elemental enthalpy corrections
- `data/cantera/` — vendored Cantera reference and product mechanisms
- `examples/molecule.xyz` — sample CH4 geometry
- `examples/molecule_data.yaml` — sample ab-mode case data for CH4
