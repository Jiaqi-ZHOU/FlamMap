from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass(slots=True)
class ProjectConfig:
    tae_hartree: float
    orca_out: Path | None
    xyz_geom: Path | None
    freqs_cm1: list[float] | None
    bond_enthalpy_json: Path
    ref_yaml: Path
    prod_yaml: Path
    output_dir: Path
    polys_temps: list[int]
    npoints: int


def _pathify(value: str | Path) -> Path:
    return Path(value).expanduser().resolve()


def _pathify_from_repo(value: str | Path, repo_root: Path) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path.resolve()
    return (repo_root / path).resolve()


def load_config(config_path: str | Path) -> ProjectConfig:
    config_path = _pathify(config_path)
    repo_root = config_path.parent.parent
    with open(config_path, "r", encoding="utf-8") as handle:
        raw: dict[str, Any] = yaml.safe_load(handle)

    inputs = raw["inputs"]
    outputs = raw["outputs"]
    parameters = raw["parameters"]

    return ProjectConfig(
        tae_hartree=float(inputs["tae_hartree"]),
        orca_out=None if inputs.get("orca_out") is None else _pathify_from_repo(inputs["orca_out"], repo_root),
        xyz_geom=None if inputs.get("xyz_geom") is None else _pathify_from_repo(inputs["xyz_geom"], repo_root),
        freqs_cm1=None if inputs.get("freqs_cm1") is None else [float(x) for x in inputs["freqs_cm1"]],
        bond_enthalpy_json=_pathify_from_repo(inputs["bond_enthalpy_json"], repo_root),
        ref_yaml=_pathify_from_repo(inputs["ref_yaml"], repo_root),
        prod_yaml=_pathify_from_repo(inputs["prod_yaml"], repo_root),
        output_dir=_pathify_from_repo(outputs["output_dir"], repo_root),
        polys_temps=[int(x) for x in parameters["polys_temps"]],
        npoints=int(parameters["npoints"]),
    )


def validate_config(cfg: ProjectConfig) -> list[str]:
    errors: list[str] = []
    for path in [cfg.bond_enthalpy_json, cfg.ref_yaml, cfg.prod_yaml]:
        if not path.exists():
            errors.append(f"Missing required path: {path}")
    if cfg.orca_out is not None and not cfg.orca_out.exists():
        errors.append(f"Missing required path: {cfg.orca_out}")
    if cfg.xyz_geom is not None and not cfg.xyz_geom.exists():
        errors.append(f"Missing required path: {cfg.xyz_geom}")

    if len(cfg.polys_temps) != 3:
        errors.append("parameters.polys_temps must have exactly 3 entries.")
    if cfg.npoints <= 1:
        errors.append("parameters.npoints must be greater than 1.")
    if not isinstance(cfg.tae_hartree, float):
        errors.append("inputs.tae_hartree must be a float.")
    if cfg.freqs_cm1 is not None:
        if not cfg.freqs_cm1:
            errors.append("inputs.freqs_cm1 must not be empty when provided.")
        elif any(freq <= 0 for freq in cfg.freqs_cm1):
            errors.append("inputs.freqs_cm1 must contain only positive frequencies in cm^-1.")
    if cfg.orca_out is None and cfg.xyz_geom is None:
        errors.append("Provide at least one geometry source: inputs.orca_out or inputs.xyz_geom.")
    if cfg.orca_out is None and cfg.freqs_cm1 is None:
        errors.append("inputs.orca_out is required when inputs.freqs_cm1 is not provided.")
    if cfg.orca_out is None and cfg.freqs_cm1 is not None and cfg.xyz_geom is None:
        errors.append("inputs.xyz_geom is required when using inputs.freqs_cm1 without inputs.orca_out.")

    return errors
