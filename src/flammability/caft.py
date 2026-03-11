from __future__ import annotations

from pathlib import Path
import warnings

import cantera as ct
import numpy as np


def compute_ternary_phase_diagram(yaml_file: str | Path, output_dir: str | Path, n_points: int = 101):
    yaml_file = Path(yaml_file)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fuel = yaml_file.stem

    try:
        gas = ct.ThermoPhase(str(yaml_file))
    except Exception as exc:
        return fuel, False, f"Failed to load YAML: {exc}"

    o2_range = np.linspace(0.0, 1.0, n_points)
    n2_range = np.linspace(0.0, 1.0, n_points)
    o2_grid, n2_grid = np.meshgrid(o2_range, n2_range)
    fuel_grid = 1.0 - o2_grid - n2_grid
    valid_mask = fuel_grid >= -1e-3
    o2_valid = o2_grid[valid_mask]
    n2_valid = n2_grid[valid_mask]
    fuel_valid = fuel_grid[valid_mask]

    temperatures = np.full_like(o2_valid, np.nan)
    for i, (n_o2, n_n2, n_fuel) in enumerate(zip(o2_valid, n2_valid, fuel_valid)):
        try:
            gas.X = {"O2": float(n_o2), "N2": float(n_n2), fuel: float(n_fuel)}
            gas.TP = 300.0, ct.one_atm
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                gas.equilibrate("HP", rtol=1e-6, max_steps=5000, solver="auto")
            temperatures[i] = gas.T
        except Exception:
            temperatures[i] = np.nan

    data = np.zeros((len(o2_valid), 4), dtype=float)
    data[:, 0] = o2_valid
    data[:, 1] = n2_valid
    data[:, 2] = fuel_valid
    data[:, 3] = temperatures

    out_path = output_dir / f"{fuel}.dat"
    np.savetxt(out_path, data, fmt="%.3f\t%.3f\t%.3f\t%.1f", header="O2\tN2\tFuel\tTemperature (K)")
    return fuel, True, str(out_path)


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
