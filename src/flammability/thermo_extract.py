from __future__ import annotations

import collections
import collections.abc
import json
import os
import re
import tempfile
from collections import defaultdict

# Compatibility shim for legacy molsym imports on modern Python.
if not hasattr(collections, "Iterable"):
    collections.Iterable = collections.abc.Iterable

import molsym
import numpy as np
from ase.units import Bohr
from molsym.molecule import Molecule


def extract_atoms(orca_outfile):
    with open(orca_outfile, "r", encoding="utf-8") as handle:
        lines = handle.readlines()

    last_idx = None
    for idx, line in enumerate(lines):
        if "CARTESIAN COORDINATES (A.U.)" in line:
            last_idx = idx

    if last_idx is None:
        raise ValueError("No 'CARTESIAN COORDINATES (A.U.)' section found.")

    atoms = []
    found_section = False
    for i in range(last_idx, len(lines)):
        line = lines[i]
        if not found_section:
            if "CARTESIAN COORDINATES (A.U.)" in line:
                found_section = True
            continue
        if "MASS" in line and "X" in line and "Y" in line and "Z" in line:
            continue
        parts = line.strip().split()
        if len(parts) == 8 and parts[0].isdigit():
            symbol = parts[1]
            mass = float(parts[4])
            x, y, z = map(float, parts[5:8])
            atoms.append({"symbol": symbol, "mass": mass, "xyz": np.array([x, y, z]) * Bohr})
        elif found_section and atoms and (len(parts) < 8 or not parts[0].isdigit()):
            break

    if not atoms:
        raise ValueError("No atoms found in the file.")
    return atoms


def compute_center_of_mass(atoms):
    total_mass = sum(atom["mass"] for atom in atoms)
    com = sum(atom["mass"] * atom["xyz"] for atom in atoms) / total_mass
    return com, total_mass


def compute_inertia_matrix(atoms, com):
    inertia_matrix = np.zeros((3, 3))
    for atom in atoms:
        r = atom["xyz"] - com
        x, y, z = r
        inertia_matrix[0, 0] += atom["mass"] * (y**2 + z**2)
        inertia_matrix[1, 1] += atom["mass"] * (x**2 + z**2)
        inertia_matrix[2, 2] += atom["mass"] * (x**2 + y**2)
        inertia_matrix[0, 1] -= atom["mass"] * x * y
        inertia_matrix[0, 2] -= atom["mass"] * x * z
        inertia_matrix[1, 2] -= atom["mass"] * y * z
    inertia_matrix[1, 0] = inertia_matrix[0, 1]
    inertia_matrix[2, 0] = inertia_matrix[0, 2]
    inertia_matrix[2, 1] = inertia_matrix[1, 2]
    return inertia_matrix


def write_temp_xyz(atoms):
    with tempfile.NamedTemporaryFile("w", suffix=".xyz", delete=False) as tmp:
        tmp.write(f"{len(atoms)}\n\n")
        for atom in atoms:
            x, y, z = atom["xyz"]
            tmp.write(f"{atom['symbol']} {x:.8f} {y:.8f} {z:.8f}\n")
        return tmp.name


def get_inertia_and_linearity(inertia_matrix):
    inertia_values = np.sort(np.linalg.eigvalsh(inertia_matrix))
    is_linear = np.isclose(inertia_values[0], 0, atol=1e-6)
    inertia = float(inertia_values[2]) if is_linear else inertia_values.tolist()
    return inertia, is_linear


def get_sigma(xyz_file, symmetry_number_json):
    mol = Molecule.from_file(xyz_file)
    mol.tol = 1e-2
    pg, _axes = molsym.find_point_group(mol)
    with open(symmetry_number_json, "r", encoding="utf-8") as handle:
        symmetry_numbers = json.load(handle)
    return symmetry_numbers.get(pg, 1)


def parse_mol_struct(orca_outfile, symmetry_number_json):
    atoms = extract_atoms(orca_outfile)
    com, total_mass = compute_center_of_mass(atoms)
    inertia_matrix = compute_inertia_matrix(atoms, com)
    inertia, is_linear = get_inertia_and_linearity(inertia_matrix)
    tmp_path = write_temp_xyz(atoms)
    sigma = get_sigma(tmp_path, symmetry_number_json)
    os.remove(tmp_path)
    return is_linear, total_mass, sigma, inertia


def parse_orca_output(orca_outfile):
    with open(orca_outfile, "r", encoding="utf-8") as handle:
        lines = handle.readlines()
    freqs_cm, zpe, total_mass, e0 = [], 0, 0, 0
    for i, line in enumerate(lines):
        if "VIBRATIONAL FREQUENCIES" in line:
            j = i + 5
            while True:
                if not re.match(r"\s*\d+:\s+([0-9.]+)", lines[j]):
                    break
                match = re.search(r"\d+:\s+([0-9.]+)", lines[j])
                val = float(match.group(1))
                if val > 10:
                    freqs_cm.append(val)
                j += 1
        if "Zero point energy" in line:
            match = re.search(r"(-?[0-9.]+) kcal/mol", line)
            if match:
                zpe = float(match.group(1))
        if "Total Mass" in line:
            match = re.search(r"(-?[0-9.]+) AMU", line)
            if match:
                total_mass = float(match.group(1))
        if "Electronic energy" in line:
            match = re.search(r"(-?[0-9.]+) Eh", line)
            if match:
                e0 = float(match.group(1))
    return freqs_cm, zpe, total_mass, e0


def count_atoms(formula):
    pattern = r"([A-Z][a-z]*)(\d*)"
    counts = defaultdict(int)
    for element, count in re.findall(pattern, formula):
        counts[element] += int(count) if count else 1
    return dict(counts)


def formula_from_atoms(atoms) -> str:
    counts = defaultdict(int)
    for atom in atoms:
        counts[atom["symbol"]] += 1

    ordered = []
    if "C" in counts:
        ordered.append("C")
    if "H" in counts:
        ordered.append("H")
    ordered.extend(sorted(elem for elem in counts if elem not in {"C", "H"}))

    formula = ""
    for elem in ordered:
        count = counts[elem]
        formula += elem if count == 1 else f"{elem}{count}"
    return formula
