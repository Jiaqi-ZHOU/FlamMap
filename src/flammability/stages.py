from __future__ import annotations

import json
from pathlib import Path

from .caft import compute_ternary_phase_diagram
from .config import ProjectConfig
from .pd_analysis import compute_flammability_limits, plot_phase_diagram
from .species import infer_case_name
from .thermo_extract import extract_atoms, formula_from_atoms
from .thermo_process import get_dft_therms
from .yaml_generation import gen_custom_yaml

C_BOND = 718.1
THRESHOLD_TEMPERATURE = 1600


def _output_paths(cfg: ProjectConfig, case_name: str) -> dict[str, Path]:
    output_dir = cfg.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    return {
        "yaml": output_dir / f"{case_name}.yaml",
        "dat": output_dir / f"{case_name}.dat",
        "pdf": output_dir / f"{case_name}.pdf",
        "json": output_dir / f"{case_name}.json",
    }


def _compute_hf(cfg: ProjectConfig, formula: str) -> float:
    _cp, hf_kj, _s = get_dft_therms(
        cfg.orca_out,
        298.15,
        formula=formula,
        tae_hartree=cfg.tae_hartree,
        symmetry_number_json=cfg.symmetry_number_json,
        bond_enthalpy_json=cfg.bond_enthalpy_json,
        c_bond=C_BOND,
        freqs_cm1=cfg.freqs_cm1,
    )
    return float(hf_kj)


def run_pipeline(cfg: ProjectConfig) -> None:
    formula = formula_from_atoms(extract_atoms(cfg.orca_out))
    case_name = infer_case_name(formula)
    species_name = case_name
    paths = _output_paths(cfg, case_name)
    freq_source = "inputs.freqs_cm1" if cfg.freqs_cm1 is not None else str(cfg.orca_out)
    print(f"Running flammability pipeline for {case_name}")
    print(f"Using TAE = {cfg.tae_hartree:.6f} Ha and frequency source: {freq_source}")
    print(f"Inferred molecular formula from ORCA structure: {formula}")

    print("1. Fitting thermochemistry and generating a Cantera YAML mechanism...")
    yaml_file = gen_custom_yaml(
        species_name,
        formula=formula,
        tae_hartree=cfg.tae_hartree,
        orca_out=cfg.orca_out,
        ref_yaml=cfg.ref_yaml,
        prod_yaml=cfg.prod_yaml,
        output_dir=cfg.output_dir,
        polys_temps=cfg.polys_temps,
        symmetry_number_json=cfg.symmetry_number_json,
        bond_enthalpy_json=cfg.bond_enthalpy_json,
        c_bond=C_BOND,
        freqs_cm1=cfg.freqs_cm1,
    )
    if yaml_file != paths["yaml"]:
        raise RuntimeError(f"Unexpected YAML output path: {yaml_file}")

    print("2. Calculating the standard enthalpy of formation at 298.15 K...")
    hf_kj = _compute_hf(cfg, formula)
    print(f"   Hf(298.15 K) = {hf_kj:.3f} kJ/mol")

    print(
        "3. Computing the CAFT ternary phase diagram grid and writing the .dat file..."
    )
    _fuel, ok, msg = compute_ternary_phase_diagram(
        paths["yaml"],
        cfg.output_dir,
        n_points=cfg.npoints,
    )
    if not ok:
        raise RuntimeError(msg)

    print(
        f"4. Extracting lower and upper flammability limits at {THRESHOLD_TEMPERATURE} K..."
    )
    lfl, ufl, _segments = compute_flammability_limits(
        paths["dat"], THRESHOLD_TEMPERATURE
    )

    print("5. Plotting the flammability phase diagram PDF...")
    plot_phase_diagram(
        paths["dat"],
        paths["pdf"],
        formula=formula,
        threshold_temperature=THRESHOLD_TEMPERATURE,
        lfl_percent=lfl,
        ufl_percent=ufl,
    )

    summary = {
        "tae_hartree": cfg.tae_hartree,
        "orca_out": str(cfg.orca_out),
        "freqs_cm1": cfg.freqs_cm1,
        "case_name": case_name,
        "formula": formula,
        "Hf_298K_kJ": hf_kj,
        "LFL_percent": lfl,
        "UFL_percent": ufl,
        "diagram_pdf": str(paths["pdf"]),
        "yaml_file": str(paths["yaml"]),
        "dat_file": str(paths["dat"]),
    }
    with open(paths["json"], "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print("Finished. Summary:")
    print(f"LFL (CAFT threshold = {THRESHOLD_TEMPERATURE} K): {lfl:.3f}%")
    print(f"UFL (CAFT threshold = {THRESHOLD_TEMPERATURE} K): {ufl:.3f}%")
    print(f"Generated YAML: {paths['yaml']}")
    print(f"CAFT DAT: {paths['dat']}")
    print(f"Diagram PDF: {paths['pdf']}")
    print(f"Summary JSON: {paths['json']}")
