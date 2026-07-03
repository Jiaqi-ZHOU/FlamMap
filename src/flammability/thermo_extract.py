from __future__ import annotations

import re
from collections import defaultdict

import numpy as np
import pymsym
from ase.data import atomic_numbers
from ase.data import atomic_masses


def extract_atoms_from_xyz(xyz_file):
    with open(xyz_file, "r", encoding="utf-8") as handle:
        lines = [line.rstrip() for line in handle]

    if len(lines) < 3:
        raise ValueError("XYZ file is too short.")

    try:
        natoms = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError("First line of XYZ file must be the atom count.") from exc

    atom_lines = lines[2 : 2 + natoms]
    if len(atom_lines) != natoms:
        raise ValueError("XYZ file does not contain the declared number of atoms.")

    atoms = []
    for line in atom_lines:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Malformed XYZ atom line: {line}")
        symbol = parts[0]
        x, y, z = map(float, parts[1:4])
        atoms.append(
            {
                "symbol": symbol,
                "mass": float(atomic_masses[atomic_numbers[symbol]]),
                "xyz": np.array([x, y, z], dtype=float),
            }
        )

    if not atoms:
        raise ValueError("No atoms found in the XYZ file.")
    return atoms


def extract_atoms(geometry_file):
    return extract_atoms_from_xyz(geometry_file)


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


def get_inertia_and_linearity(inertia_matrix):
    inertia_values = np.sort(np.linalg.eigvalsh(inertia_matrix))
    # Linear molecules have a vanishing smallest principal moment of inertia (all atoms
    # on one axis), while the other two are large and equal. Test the smallest RELATIVE
    # to the largest, not against an absolute atol: coordinate noise (~1e-4 A) inflates
    # the smallest moment to ~1e-6-1e-5 amu*A^2, which an absolute 1e-6 cutoff misreads
    # as non-linear -- then the 3N-5 vibration count disagrees with this linearity and
    # the pipeline raises FrequencyCountMismatchError. The relative test mirrors the
    # rot_thresh fix in hip's frequency analysis so the two determinations agree.
    largest = inertia_values[2]
    is_linear = bool(largest > 0 and inertia_values[0] < 1e-3 * largest)
    inertia = float(inertia_values[2]) if is_linear else inertia_values.tolist()
    return inertia, is_linear


def get_sigma(atoms):
    atom_numbers = [atomic_numbers[atom["symbol"]] for atom in atoms]
    positions = [atom["xyz"].tolist() for atom in atoms]
    return int(pymsym.get_symmetry_number(atom_numbers, positions))


def parse_mol_struct(geometry_file):
    atoms = extract_atoms(geometry_file)
    com, total_mass = compute_center_of_mass(atoms)
    inertia_matrix = compute_inertia_matrix(atoms, com)
    inertia, is_linear = get_inertia_and_linearity(inertia_matrix)
    sigma = get_sigma(atoms)
    return is_linear, total_mass, sigma, inertia, len(atoms)


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
