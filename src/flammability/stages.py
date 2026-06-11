from __future__ import annotations

import json
import time
from pathlib import Path

from .caft import compute_ternary_phase_diagram
from .config import DEFAULT_THRESHOLD_TEMPERATURE, ProjectConfig
from .pd_analysis import compute_split_flammability_limits, plot_phase_diagram
from .species import infer_case_name
from .thermo_extract import extract_atoms, formula_from_atoms
from .thermo_process import get_thermo
from .yaml_generation import gen_custom_yaml

# Backward-compatible alias; live values come from cfg.lfl_threshold / cfg.ufl_threshold per run.
THRESHOLD_TEMPERATURE = DEFAULT_THRESHOLD_TEMPERATURE


def _output_paths(cfg: ProjectConfig, case_name: str) -> dict[str, Path]:
    """Per-molecule outputs go into per-type subdirs under ``cfg.output_dir``.

    Layout:
        <output_dir>/yaml/<case>.yaml
        <output_dir>/dat/<case>.dat
        <output_dir>/pdf/<case>.pdf   (only when --plot)
        <output_dir>/json/<case>.json

    Batch-level files (SUMMARY.csv, FAILED.csv) stay at the root so they're
    easy to find without descending into a subdir.
    """
    output_dir = cfg.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    subdirs = {kind: output_dir / kind for kind in ("yaml", "dat", "pdf", "json")}
    for d in subdirs.values():
        d.mkdir(parents=True, exist_ok=True)
    return {
        "yaml": subdirs["yaml"] / f"{case_name}.yaml",
        "dat": subdirs["dat"] / f"{case_name}.dat",
        "pdf": subdirs["pdf"] / f"{case_name}.pdf",
        "json": subdirs["json"] / f"{case_name}.json",
    }


def _resolve_tae_and_freqs(
    cfg: ProjectConfig, *, quiet: bool = False
) -> tuple[float, list[float], str]:
    if cfg.mode == "ml":
        from .ml import compute_freqs, compute_tae, resolve_device_for_xyz

        skala_device = resolve_device_for_xyz(cfg.skala_device, cfg.xyz_geom)
        hip_device = resolve_device_for_xyz(cfg.hip_device, cfg.xyz_geom)

        # Order matters: hip Hessian (~10–30 s) runs first so molecules with
        # imaginary modes (non-minimum geometries) fail fast without paying for
        # the much more expensive skala SCFs (~30–180 s/molecule).
        if not quiet:
            print(
                f"ML mode: computing vibrational frequencies via hip (device={hip_device})..."
            )
        freqs = compute_freqs(
            cfg.xyz_geom,
            checkpoint_path=cfg.hip_checkpoint,
            device=hip_device,
            quiet=quiet,
        )
        if not quiet:
            print(f"   {len(freqs)} vibrational frequencies (cm^-1):")
            for i in range(0, len(freqs), 6):
                chunk = "  ".join(f"{f:8.2f}" for f in freqs[i : i + 6])
                print(f"     {chunk}")

        if not quiet:
            print(f"ML mode: computing TAE via skala (device={skala_device})...")
        tae = compute_tae(cfg.xyz_geom, device=skala_device, quiet=quiet)
        if not quiet:
            print(f"   TAE = {tae:.6f} Ha")
        return tae, freqs, "ml (hip)"

    return cfg.tae, cfg.freqs, "case YAML"


def run_pipeline(
    cfg: ProjectConfig,
    *,
    quiet: bool = False,
    case_name_override: str | None = None,
) -> dict:
    t0 = time.perf_counter()
    geometry_file = cfg.xyz_geom
    formula = formula_from_atoms(extract_atoms(geometry_file))
    case_name = case_name_override or infer_case_name(formula)
    species_name = case_name
    paths = _output_paths(cfg, case_name)

    if not quiet:
        print(f"Running flammability pipeline for {case_name} (mode={cfg.mode})")
        print(f"Using geometry source: {geometry_file}")
        print(f"Inferred molecular formula from geometry: {formula}")

    tae, freqs, freq_source = _resolve_tae_and_freqs(cfg, quiet=quiet)
    if not quiet:
        print(f"Using TAE = {tae:.6f} Ha and frequency source: {freq_source}")

    if not quiet:
        print("1. Fitting thermochemistry and generating a Cantera YAML mechanism...")
    yaml_file = gen_custom_yaml(
        species_name,
        formula=formula,
        tae=tae,
        geometry_file=geometry_file,
        output_dir=paths["yaml"].parent,
        polys_temps=cfg.polys_temps,
        bond_enthalpy_json=cfg.bond_enthalpy_json,
        freqs=freqs,
    )
    if yaml_file != paths["yaml"]:
        raise RuntimeError(f"Unexpected YAML output path: {yaml_file}")

    if not quiet:
        print("2. Calculating the standard enthalpy of formation at 298.15 K...")
    _cp, hf_kj, _s = get_thermo(
        geometry_file,
        298.15,
        formula=formula,
        tae=tae,
        bond_enthalpy_json=cfg.bond_enthalpy_json,
        freqs=freqs,
    )
    hf_kj = float(hf_kj)
    if not quiet:
        print(f"   Hf(298.15 K) = {hf_kj:.3f} kJ/mol")
        print(
            "3. Computing the CAFT ternary phase diagram grid and writing the .dat file..."
        )
    _fuel, ok, msg = compute_ternary_phase_diagram(
        paths["yaml"],
        paths["dat"].parent,
        n_points=cfg.npoints,
    )
    if not ok:
        raise RuntimeError(msg)

    if not quiet:
        if cfg.lfl_threshold == cfg.ufl_threshold:
            print(
                f"4. Extracting lower and upper flammability limits at {cfg.lfl_threshold:g} K..."
            )
        else:
            print(
                "4. Extracting flammability limits "
                f"(LFL at {cfg.lfl_threshold:g} K, UFL at {cfg.ufl_threshold:g} K)..."
            )
    lfl, ufl, _lfl_segs, _ufl_segs = compute_split_flammability_limits(
        paths["dat"], cfg.lfl_threshold, cfg.ufl_threshold
    )

    diagram_pdf: str | None = str(paths["pdf"])
    if cfg.plot_map:
        if not quiet:
            print("5. Plotting the flammability phase diagram PDF...")
        plot_phase_diagram(
            paths["dat"],
            paths["pdf"],
            formula=formula,
            lfl_threshold_temperature=cfg.lfl_threshold,
            ufl_threshold_temperature=cfg.ufl_threshold,
            lfl_percent=lfl,
            ufl_percent=ufl,
        )
    else:
        diagram_pdf = None
        if not quiet:
            print("5. plot_map is false; skipping phase diagram PDF generation.")

    summary = {
        "mode": cfg.mode,
        "tae": tae,
        "xyz_geom": str(cfg.xyz_geom),
        "freqs": freqs,
        "case_name": case_name,
        "formula": formula,
        "Hf_298K_kJ": hf_kj,
        "lfl_threshold_K": cfg.lfl_threshold,
        "ufl_threshold_K": cfg.ufl_threshold,
        "LFL_percent": lfl,
        "UFL_percent": ufl,
        "diagram_pdf": diagram_pdf,
        "yaml_file": str(paths["yaml"]),
        "dat_file": str(paths["dat"]),
        "elapsed_s": round(time.perf_counter() - t0, 3),
    }
    with open(paths["json"], "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    if quiet:
        return summary

    print("Finished. Summary:")
    print(f"LFL (CAFT threshold = {cfg.lfl_threshold:g} K): {lfl:.3f}%")
    print(f"UFL (CAFT threshold = {cfg.ufl_threshold:g} K): {ufl:.3f}%")
    print(f"Generated YAML: {paths['yaml']}")
    print(f"CAFT DAT: {paths['dat']}")
    if cfg.plot_map:
        print(f"Diagram PDF: {paths['pdf']}")
    else:
        print("Diagram PDF: skipped")
    print(f"Summary JSON: {paths['json']}")
    return summary
