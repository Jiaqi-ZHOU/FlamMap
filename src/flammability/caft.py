from __future__ import annotations

from pathlib import Path
import multiprocessing as mp
import os
import warnings
from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout, redirect_stderr

import cantera as ct
import numpy as np

# Products are NOT baked into the fuel YAML. The per-fuel YAML carries only the one
# species we computed (its NASA7 thermo); the combustion products are brought in here,
# at equilibrium time, from Cantera's bundled mechanisms. This keeps the generated YAML
# a clean record of "the thermochemistry I produced" and puts every product (gas + soot)
# on the same footing instead of splitting them between the YAML and this file.
PRODUCTS_YAML = "gri30.yaml"     # Cantera-bundled GRI-3.0 gas-phase product species
GRAPHITE_YAML = "graphite.yaml"  # Cantera-bundled condensed carbon (soot), NASA/McBride C(gr)

# Toggle for the condensed-carbon (soot) sink. False (default) = gas-only HP equilibrium.
# True = multiphase gas+graphite equilibrium, which adds the condensed-carbon sink that
# lowers rich-side T_ad but is much slower to converge on some fuels.
INCLUDE_SOOT = False


def _build_gas_and_carbon(yaml_file, include_soot):
    """Build the runtime gas phase and (optionally) the condensed-carbon phase.

    The gas phase is GRI-3.0 product species + the single fuel species from this YAML,
    merged into one homogeneous ideal-gas phase (they must share one phase so the
    multiphase solver gets the gas-gas mixing entropy right; the fuel is NOT a separate
    phase). On the rare name clash with a GRI-3.0 species (e.g. a bare "CH4" fuel), the
    fuel is suffixed "_fuel" and coexists with GRI-3.0's own species rather than replacing
    it. Normal `formula_stem` names never clash, so no dedup.

    Returns ``(gas, carbon, fuel_key)`` where ``carbon`` is None when ``include_soot`` is
    False. Raises RuntimeError with a human-readable message on any failure.
    """
    try:
        products = ct.Species.list_from_file(PRODUCTS_YAML)
        fuel_species = ct.Species.list_from_file(str(yaml_file))
    except Exception as exc:
        raise RuntimeError(f"Failed to load species: {exc}")
    if len(fuel_species) != 1:
        raise RuntimeError(
            f"Expected exactly one fuel species in {Path(yaml_file).name}, "
            f"found {len(fuel_species)}. This YAML predates the runtime-merge format "
            "(it bundled the GRI-3.0 products); regenerate it."
        )
    fuel_sp = fuel_species[0]
    product_names = {sp.name for sp in products}
    fuel_key = fuel_sp.name
    if fuel_key in product_names:
        fuel_key = f"{fuel_key}_fuel"
        renamed = ct.Species(fuel_key, fuel_sp.composition)
        renamed.thermo = fuel_sp.thermo
        fuel_sp = renamed
    try:
        gas = ct.Solution(
            thermo="ideal-gas", species=list(products) + [fuel_sp], transport_model="none"
        )
    except Exception as exc:
        raise RuntimeError(f"Failed to build gas phase: {exc}")

    # Condensed-phase carbon (soot) sink: a SEPARATE solid phase, not a gas-phase product.
    # Equilibrium is a multiphase Gibbs minimisation (gas + graphite), the standard
    # treatment for fuel-rich equilibria; gas-only omits soot and overestimates rich-side
    # T_ad. Disabled when INCLUDE_SOOT is False (gas-only run).
    carbon = None
    if include_soot:
        try:
            carbon = ct.Solution(GRAPHITE_YAML)
        except Exception as exc:
            raise RuntimeError(f"Failed to load graphite phase: {exc}")
    return gas, carbon, fuel_key


def _equilibrate_point(gas, carbon, fuel_key, include_soot, n_o2, n_n2, n_fuel):
    """Multiphase HP equilibrium temperature for one composition, NaN if it fails.

    vcs -> gibbs fallback: vcs is the more robust multiphase solver (handles
    condensed-phase appearance/disappearance and converges on exotic/energetic fuels
    where gibbs fails); on the rare point where vcs fails, retry with gibbs; only then
    record NaN. Both solvers give identical T where they converge.
    """
    comp = {"O2": float(n_o2), "N2": float(n_n2), fuel_key: float(n_fuel)}
    for solver in ("vcs", "gibbs"):
        try:
            gas.TPX = 300.0, ct.one_atm, comp  # reset to clean initial state each attempt
            phases = [(gas, 1.0)] + ([(carbon, 0.0)] if include_soot else [])
            mix = ct.Mixture(phases)
            mix.T = 300.0
            mix.P = ct.one_atm
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                mix.equilibrate("HP", solver=solver, max_steps=5000)
            return mix.T
        except Exception:
            continue
    return float("nan")


def _equilibrate_grid(gas, carbon, fuel_key, include_soot, points):
    """Run a list of (index, n_o2, n_n2, n_fuel) points; return (results, n_fail).

    Cantera's multiphase solvers emit non-fatal "FAILURE its=..." diagnostics on the
    output stream even when they converge; silence them around the sweep so a large batch
    does not flood logs. The returned n_fail is the real convergence signal.
    """
    results = []
    n_fail = 0
    with open(os.devnull, "w") as _dn, redirect_stdout(_dn), redirect_stderr(_dn):
        for idx, n_o2, n_n2, n_fuel in points:
            t = _equilibrate_point(gas, carbon, fuel_key, include_soot, n_o2, n_n2, n_fuel)
            if not np.isfinite(t):
                n_fail += 1
            results.append((idx, t))
    return results, n_fail


# --- parallel worker plumbing (spawn-based; each worker builds its own cantera phases,
# never imports torch, so it stays ~tens of MB regardless of molecule) ---------------
_WORKER = {}


def _worker_init(yaml_file, include_soot):
    try:
        gas, carbon, fuel_key = _build_gas_and_carbon(yaml_file, include_soot)
        _WORKER.update(gas=gas, carbon=carbon, fuel_key=fuel_key, include_soot=include_soot, err=None)
    except Exception as exc:  # surface the build error via every chunk this worker runs
        _WORKER.update(err=str(exc))


def _worker_chunk(points):
    if _WORKER.get("err"):
        return [(idx, float("nan")) for idx, *_ in points], len(points)
    return _equilibrate_grid(
        _WORKER["gas"], _WORKER["carbon"], _WORKER["fuel_key"], _WORKER["include_soot"], points
    )


def compute_ternary_phase_diagram(
    yaml_file: str | Path, output_dir: str | Path, n_points: int = 101, n_workers: int = 1
):
    """Compute the CAFT ternary grid for one fuel and write ``<fuel>.dat``.

    The ~n_points^2/2 grid points are independent multiphase HP equilibria. With
    ``n_workers > 1`` they are spread across that many spawn processes; each worker
    rebuilds the cantera phases once (it never touches torch, so it costs only tens of MB)
    and the per-molecule memory peak — which is set elsewhere by the ML/SCF stage, not by
    CAFT — is unaffected. ``n_workers == 1`` keeps the original serial path exactly.
    """
    yaml_file = Path(yaml_file)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fuel = yaml_file.stem
    include_soot = INCLUDE_SOOT

    # Validate the phase build up front in this process so a bad YAML (e.g. the legacy
    # fat format) returns the precise error message rather than a grid of silent NaNs.
    try:
        gas, carbon, fuel_key = _build_gas_and_carbon(yaml_file, include_soot)
    except Exception as exc:
        return fuel, False, str(exc)

    o2_range = np.linspace(0.0, 1.0, n_points)
    n2_range = np.linspace(0.0, 1.0, n_points)
    o2_grid, n2_grid = np.meshgrid(o2_range, n2_range)
    fuel_grid = 1.0 - o2_grid - n2_grid
    valid_mask = fuel_grid >= -1e-3
    o2_valid = o2_grid[valid_mask]
    n2_valid = n2_grid[valid_mask]
    fuel_valid = fuel_grid[valid_mask]
    npts = len(o2_valid)

    temperatures = np.full(npts, np.nan)
    n_fail = 0
    points = list(zip(range(npts), o2_valid.tolist(), n2_valid.tolist(), fuel_valid.tolist()))

    if n_workers and n_workers > 1 and npts > 0:
        # Round-robin chunks keep the per-worker load balanced. Spawn (not fork) gives
        # fresh torch-free processes that match the ~tens-of-MB footprint measured for a
        # cantera-only worker; fork would COW-inherit the parent's loaded ML models.
        chunks = [points[j::n_workers] for j in range(n_workers)]
        chunks = [c for c in chunks if c]
        try:
            ctx = mp.get_context("spawn")
            with ProcessPoolExecutor(
                max_workers=len(chunks),
                mp_context=ctx,
                initializer=_worker_init,
                initargs=(str(yaml_file), include_soot),
            ) as pool:
                for results, nf in pool.map(_worker_chunk, chunks):
                    n_fail += nf
                    for idx, t in results:
                        temperatures[idx] = t
        except Exception as exc:
            return fuel, False, f"Parallel CAFT failed: {exc}"
    else:
        results, n_fail = _equilibrate_grid(gas, carbon, fuel_key, include_soot, points)
        for idx, t in results:
            temperatures[idx] = t

    data = np.zeros((npts, 4), dtype=float)
    data[:, 0] = o2_valid
    data[:, 1] = n2_valid
    data[:, 2] = fuel_valid
    data[:, 3] = temperatures

    out_path = output_dir / f"{fuel}.dat"
    np.savetxt(out_path, data, fmt="%.3f\t%.3f\t%.3f\t%.1f", header="O2\tN2\tFuel\tTemperature (K)")
    return fuel, True, f"{out_path}  (fail_pts={n_fail}/{npts}, workers={max(1, n_workers)})"


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
