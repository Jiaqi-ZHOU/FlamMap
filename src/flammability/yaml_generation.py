from __future__ import annotations

from datetime import datetime
from pathlib import Path

import numpy as np
import yaml

from .nasa import fit_nasa7_poly
from .species import canonical_species_name
from .thermo_extract import count_atoms
from .thermo_process import get_thermo_temps


def fit_temp_range(
    *,
    formula,
    tae,
    geometry_file,
    polys_temps,
    bond_enthalpy_json,
    freqs,
):
    temps_low = np.arange(polys_temps[0], polys_temps[1] + 1, 100)
    temps_high = np.arange(polys_temps[1], polys_temps[2] + 1, 100)

    def get_fit_coeffs(temps):
        cp_kj, hf_kj, s_kj = get_thermo_temps(
            geometry_file,
            temps,
            formula=formula,
            tae=tae,
            bond_enthalpy_json=bond_enthalpy_json,
            freqs=freqs,
        )
        cp = cp_kj * 1000
        hf = hf_kj * 1000
        entropy = s_kj * 1000
        return [float(x) for x in fit_nasa7_poly(temps, cp, hf, entropy)]

    return {"fit_low": get_fit_coeffs(temps_low), "fit_high": get_fit_coeffs(temps_high)}


def gen_custom_yaml(
    species_id,
    *,
    formula,
    tae,
    geometry_file,
    output_dir,
    polys_temps,
    bond_enthalpy_json,
    freqs,
):
    formula = canonical_species_name(formula)
    atom_counts = count_atoms(formula)
    fits = fit_temp_range(
        formula=formula,
        tae=tae,
        geometry_file=geometry_file,
        polys_temps=polys_temps,
        bond_enthalpy_json=bond_enthalpy_json,
        freqs=freqs,
    )

    temps_low = np.arange(polys_temps[0], polys_temps[1] + 1, 100)
    temps_high = np.arange(polys_temps[1], polys_temps[2] + 1, 100)
    entry = {
        "name": species_id,
        "composition": atom_counts,
        "thermo": {
            "model": "NASA7",
            "temperature-ranges": [int(temps_low[0]), int(temps_low[-1]), int(temps_high[-1])],
            "data": [fits["fit_low"], fits["fit_high"]],
        },
        "note": "Custom W1-F12/B3LYP",
    }

    # Emit a "thin" YAML: only the species we computed, no combustion products. The
    # products (GRI-3.0 gas + graphite soot) are merged in at equilibrium time by
    # caft.compute_ternary_phase_diagram from Cantera's bundled mechanisms, so this file
    # stays a clean record of the produced thermochemistry. A species-only document (no
    # `phases` block) is valid input for ct.Species.list_from_file.
    document = {
        "date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "species": [entry],
    }

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    outfile = output_dir / f"{species_id}.yaml"
    with open(outfile, "w", encoding="utf-8") as handle:
        yaml.dump(document, handle, sort_keys=False, width=1000, default_flow_style=None)
    return outfile
