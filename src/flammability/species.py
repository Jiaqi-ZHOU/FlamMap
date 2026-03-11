from __future__ import annotations


def canonical_species_name(formula: str) -> str:
    formula = str(formula).strip()
    return "CH3OH" if formula == "CH4O" else formula


def infer_case_name(formula: str) -> str:
    """Use the inferred chemical formula as the output stem."""
    return canonical_species_name(formula)
