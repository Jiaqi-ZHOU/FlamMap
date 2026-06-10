from __future__ import annotations

from pathlib import Path
import os
import warnings
from contextlib import redirect_stdout, redirect_stderr

import cantera as ct
import numpy as np


def compute_ternary_phase_diagram(yaml_file: str | Path, output_dir: str | Path, n_points: int = 101):
    yaml_file = Path(yaml_file)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fuel = yaml_file.stem

    try:
        gas = ct.Solution(str(yaml_file))
    except Exception as exc:
        return fuel, False, f"Failed to load YAML: {exc}"

    # Condensed-phase carbon (soot) sink: a SEPARATE solid phase, not a gas-phase
    # product. Equilibrium below is a multiphase Gibbs minimisation (gas + graphite),
    # the standard treatment for fuel-rich equilibria; gas-only omits soot and
    # overestimates rich-side T_ad. Fuel YAMLs are unchanged. Uses Cantera's bundled
    # graphite.yaml (NASA/McBride C(gr), a reference element with Delta_fH298 = 0).
    # To revert to gas-only: `git show d9c55ed^:src/flammability/caft.py`.
    try:
        carbon = ct.Solution("graphite.yaml")
    except Exception as exc:
        return fuel, False, f"Failed to load graphite phase: {exc}"

    o2_range = np.linspace(0.0, 1.0, n_points)
    n2_range = np.linspace(0.0, 1.0, n_points)
    o2_grid, n2_grid = np.meshgrid(o2_range, n2_range)
    fuel_grid = 1.0 - o2_grid - n2_grid
    valid_mask = fuel_grid >= -1e-3
    o2_valid = o2_grid[valid_mask]
    n2_valid = n2_grid[valid_mask]
    fuel_valid = fuel_grid[valid_mask]

    # Multiphase HP equilibrium with a vcs -> gibbs fallback chain. vcs is the more
    # robust multiphase solver (handles condensed-phase appearance/disappearance and
    # converges on exotic/energetic fuels where gibbs fails); on the rare point where
    # vcs fails, retry with gibbs; only then record NaN. Both solvers give identical T
    # where they converge, so this is a robustness measure, not a physics change.
    # Cantera's multiphase solvers emit non-fatal "FAILURE its=..." diagnostics on the
    # output stream even when they converge; silence them around the grid sweep so a
    # large batch does not flood logs. n_fail below is the real convergence signal.
    temperatures = np.full_like(o2_valid, np.nan)
    n_fail = 0
    with open(os.devnull, "w") as _dn, redirect_stdout(_dn), redirect_stderr(_dn):
        for i, (n_o2, n_n2, n_fuel) in enumerate(zip(o2_valid, n2_valid, fuel_valid)):
            comp = {"O2": float(n_o2), "N2": float(n_n2), fuel: float(n_fuel)}
            result = np.nan
            for solver in ("vcs", "gibbs"):
                try:
                    gas.TPX = 300.0, ct.one_atm, comp  # reset to clean initial state each attempt
                    mix = ct.Mixture([(gas, 1.0), (carbon, 0.0)])
                    mix.T = 300.0
                    mix.P = ct.one_atm
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        mix.equilibrate("HP", solver=solver, max_steps=5000)
                    result = mix.T
                    break
                except Exception:
                    continue
            if not np.isfinite(result):
                n_fail += 1
            temperatures[i] = result

    data = np.zeros((len(o2_valid), 4), dtype=float)
    data[:, 0] = o2_valid
    data[:, 1] = n2_valid
    data[:, 2] = fuel_valid
    data[:, 3] = temperatures

    out_path = output_dir / f"{fuel}.dat"
    np.savetxt(out_path, data, fmt="%.3f\t%.3f\t%.3f\t%.1f", header="O2\tN2\tFuel\tTemperature (K)")
    return fuel, True, f"{out_path}  (fail_pts={n_fail}/{len(o2_valid)})"


def generate_all_caft(yaml_dir: str | Path, output_dir: str | Path, *, n_points: int):
    yaml_dir = Path(yaml_dir)
    output_dir = Path(output_dir)
    yaml_files = sorted(p for p in yaml_dir.iterdir() if p.suffix in {".yaml", ".yml"})
    total = len(yaml_files)
    if total == 0:
        print(f"No YAML files found in {yaml_dir}")
        return

    failed = 0
    for index, yaml_file in enumerate(yaml_files, start=1):
        fuel, ok, msg = compute_ternary_phase_diagram(yaml_file, output_dir, n_points)
        status = "OK" if ok else "FAILED"
        if not ok:
            failed += 1
        print(f"[{index:>4}/{total}] {status:6} {fuel} -> {msg}", flush=True)

    if failed:
        raise RuntimeError(f"{failed} CAFT jobs failed.")
