from __future__ import annotations

from datetime import datetime
from pathlib import Path

import numpy as np
import yaml

from .nasa import fit_nasa7_poly
from .species import canonical_species_name
from .thermo_extract import count_atoms
from .thermo_process import get_dft_therms_temps


def fit_temp_range(
    *,
    formula,
    tae_hartree,
    geometry_file,
    polys_temps,
    bond_enthalpy_json,
    c_bond,
    orca_outfile=None,
    freqs_cm1=None,
):
    temps_low = np.arange(polys_temps[0], polys_temps[1] + 1, 100)
    temps_high = np.arange(polys_temps[1], polys_temps[2] + 1, 100)

    def get_fit_coeffs(temps):
        cp_kj, hf_kj, s_kj = get_dft_therms_temps(
            geometry_file,
            temps,
            formula=formula,
            tae_hartree=tae_hartree,
            bond_enthalpy_json=bond_enthalpy_json,
            c_bond=c_bond,
            orca_outfile=orca_outfile,
            freqs_cm1=freqs_cm1,
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
    tae_hartree,
    geometry_file,
    ref_yaml,
    prod_yaml,
    output_dir,
    polys_temps,
    bond_enthalpy_json,
    c_bond,
    orca_out=None,
    freqs_cm1=None,
):
    formula = canonical_species_name(formula)
    atom_counts = count_atoms(formula)
    fits = fit_temp_range(
        formula=formula,
        tae_hartree=tae_hartree,
        geometry_file=geometry_file,
        polys_temps=polys_temps,
        bond_enthalpy_json=bond_enthalpy_json,
        c_bond=c_bond,
        orca_outfile=orca_out,
        freqs_cm1=freqs_cm1,
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

    with open(prod_yaml, "r", encoding="utf-8") as handle:
        base_yaml = yaml.safe_load(handle)

    filtered_species = []
    for sp in base_yaml.get("species", []):
        sp.pop("transport", None)
        if sp.get("name") != species_id:
            filtered_species.append(sp)
    base_yaml["species"] = filtered_species

    phase_species = ["NO" if s is False else s for s in base_yaml["phases"][0]["species"]]
    phase_species = [name for name in phase_species if name != species_id]
    base_yaml["phases"][0]["species"] = phase_species
    base_yaml["phases"][0].pop("kinetics", None)
    base_yaml["phases"][0].pop("transport", None)

    for key in ["description", "generator", "input-files", "cantera-version", "units"]:
        base_yaml.pop(key, None)
    base_yaml["date"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    base_yaml["phases"][0]["name"] = f"GRI3.0Products+{species_id}"
    base_yaml["phases"][0]["species"].append(species_id)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    outfile = output_dir / f"{species_id}.yaml"
    with open(outfile, "w", encoding="utf-8") as handle:
        yaml.dump(base_yaml, handle, sort_keys=False, width=1000, default_flow_style=None)
        yaml.dump([entry], handle, sort_keys=False, width=1000, default_flow_style=None)
    return outfile
