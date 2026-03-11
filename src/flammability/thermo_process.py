from __future__ import annotations

import json

import numpy as np
import pandas as pd

from .thermo_constants import Na, R, c_cm, cal2J, h, ha2kJ, kB
from .thermo_extract import count_atoms, parse_mol_struct, parse_orca_output


def calc_enthalpy_trans_rot(is_linear, temperature):
    h_trans = 5 / 2 * R * temperature
    h_rot = R * temperature if is_linear else (3 / 2) * R * temperature
    return h_trans, h_rot


def calc_entropy_trans_rot(molecular_weight, sigma, inertia, temperature, pressure=1.0):
    amu2kg = 1.66053906660e-27
    angstrom2m = 1e-10
    atm2pa = 101325

    mass = molecular_weight * amu2kg
    pressure_pa = pressure * atm2pa

    lambda_factor = (2 * np.pi * mass * kB * temperature / h**2) ** 1.5
    v_per_molecule = kB * temperature / pressure_pa
    q_trans = lambda_factor * v_per_molecule
    s_trans = R * (np.log(q_trans) + 5 / 2)

    if np.isscalar(inertia) or np.allclose(inertia[0], 0, atol=1e-6):
        inertia_value = inertia * amu2kg * angstrom2m**2
        q_rot = 8 * np.pi**2 * inertia_value * kB * temperature / (sigma * h**2)
        s_rot = R * (np.log(q_rot) + 1)
    else:
        inertia_value = np.array(inertia) * amu2kg * angstrom2m**2
        q_rot = ((np.pi * np.prod(inertia_value)) ** 0.5 * (8 * np.pi**2 * kB * temperature) ** 1.5) / (
            sigma * h**3
        )
        s_rot = R * (np.log(q_rot) + 3 / 2)

    return s_trans, s_rot


def calc_heat_capacity_trans_rot(is_linear):
    cv_trans = (3 / 2) * R
    cp_trans = cv_trans + R
    cv_rot = R if is_linear else (3 / 2) * R
    cp_rot = cv_rot
    return cp_trans, cp_rot


def calc_vibrational_thermo(freqs_cm1, temperature):
    freqs_hz = np.array(freqs_cm1) * c_cm
    theta_vib = h * freqs_hz / kB
    x = theta_vib / temperature
    cp_vib = R * np.sum((x**2 * np.exp(x)) / (np.exp(x) - 1) ** 2)
    s_vib = R * np.sum(x / (np.exp(x) - 1) - np.log(1 - np.exp(-x)))
    h_vib = R * temperature * np.sum(x / (np.exp(x) - 1))
    return cp_vib, h_vib, s_vib


def calc_zpe_from_freqs(freqs_cm1):
    zpe_j_per_mol = 0.5 * h * c_cm * Na * np.sum(freqs_cm1)
    return zpe_j_per_mol / (1000 * cal2J)


def expected_vibrational_mode_count(natoms, is_linear):
    if natoms < 2:
        return 0
    return 3 * natoms - (5 if is_linear else 6)


def _resolve_freqs_cm1(orca_outfile, freqs_cm1):
    if freqs_cm1 is not None:
        parsed_freqs_cm = [float(freq) for freq in freqs_cm1]
        return parsed_freqs_cm, calc_zpe_from_freqs(parsed_freqs_cm), None, None

    if orca_outfile is None:
        raise ValueError("inputs.orca_out is required when vibrational frequencies are not provided directly.")

    parsed_freqs_cm, zpe_kcal, total_mass, e0 = parse_orca_output(orca_outfile)
    assert parsed_freqs_cm, "No vibrational frequencies found in ORCA output."
    return parsed_freqs_cm, zpe_kcal, total_mass, e0


def get_dft_quantities(geometry_file, tae_hartree, *, orca_outfile=None, freqs_cm1=None):
    is_linear, molecular_weight, sigma, inertia, natoms = parse_mol_struct(geometry_file)
    freqs_cm, zpe_kcal, total_mass, _e0 = _resolve_freqs_cm1(orca_outfile, freqs_cm1)
    if total_mass is not None:
        assert abs(total_mass - molecular_weight) < 1e-2, "Total mass mismatch in ORCA output."
    expected_nfreq = expected_vibrational_mode_count(natoms, is_linear)
    assert len(freqs_cm) == expected_nfreq, (
        f"Vibrational frequency count mismatch: got {len(freqs_cm)} non-zero modes, "
        f"expected {expected_nfreq} for {'linear' if is_linear else 'nonlinear'} {natoms}-atom geometry."
    )
    return freqs_cm, zpe_kcal, tae_hartree, is_linear, molecular_weight, sigma, inertia


def get_dft_therms(
    geometry_file,
    temperature,
    *,
    formula,
    tae_hartree,
    bond_enthalpy_json,
    c_bond,
    orca_outfile=None,
    freqs_cm1=None,
):
    freqs_cm, zpe_kcal, tae, is_linear, molecular_weight, sigma, inertia = get_dft_quantities(
        geometry_file, tae_hartree, orca_outfile=orca_outfile, freqs_cm1=freqs_cm1
    )
    zpe_kj = 0.0 if zpe_kcal is None else zpe_kcal * cal2J
    cp_vib, h_vib, s_vib = calc_vibrational_thermo(freqs_cm, temperature)
    cp_trans, cp_rot = calc_heat_capacity_trans_rot(is_linear)
    s_trans, s_rot = calc_entropy_trans_rot(molecular_weight, sigma, inertia, temperature, pressure=1.0)
    h_trans, h_rot = calc_enthalpy_trans_rot(is_linear, temperature)
    h_mol_corr_kj = (h_trans + h_rot + h_vib) / 1e3

    atom_counts = count_atoms(formula)

    with open(bond_enthalpy_json, "r", encoding="utf-8") as handle:
        elem_dict = json.load(handle)
    elem_corr_kj = 0.0
    for elem, count in atom_counts.items():
        if elem == "C":
            atom_corr = c_bond
        elif elem == "H":
            atom_corr = elem_dict["H-H"]
        elif elem == "O":
            atom_corr = elem_dict["O-O"]
        elif elem == "N":
            atom_corr = elem_dict["N-N"]
        else:
            raise KeyError(f"Unsupported element in formula: {elem}")
        elem_corr_kj += atom_corr * count

    num_atom = sum(atom_counts.values())
    h_atom_corr_kj = (5 / 2 * R * 298.15 * num_atom) / 1e3
    tae_kj = tae * ha2kJ
    hf_kj = -tae_kj + zpe_kj + elem_corr_kj + h_mol_corr_kj - h_atom_corr_kj
    s_kj = (s_trans + s_rot + s_vib) / 1e3
    cp_kj = (cp_trans + cp_rot + cp_vib) / 1e3
    return cp_kj, hf_kj, s_kj


def get_dft_therms_temps(
    geometry_file,
    temps,
    *,
    formula,
    tae_hartree,
    bond_enthalpy_json,
    c_bond,
    orca_outfile=None,
    freqs_cm1=None,
):
    cp_kj_array = np.zeros(len(temps))
    hf_kj_array = np.zeros(len(temps))
    s_kj_array = np.zeros(len(temps))
    for i, temperature in enumerate(temps):
        cp, hf, entropy = get_dft_therms(
            geometry_file,
            temperature,
            formula=formula,
            tae_hartree=tae_hartree,
            bond_enthalpy_json=bond_enthalpy_json,
            c_bond=c_bond,
            orca_outfile=orca_outfile,
            freqs_cm1=freqs_cm1,
        )
        cp_kj_array[i], hf_kj_array[i], s_kj_array[i] = cp, hf, entropy
    return cp_kj_array, hf_kj_array, s_kj_array
