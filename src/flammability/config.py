from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

DEFAULT_BOND_ENTHALPY_JSON = REPO_ROOT / "data" / "reference" / "elem_enthalpies.json"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "outputs"
DEFAULT_POLYS_TEMPS = [200, 1000, 3600]
DEFAULT_NPOINTS = 101
DEFAULT_THRESHOLD_TEMPERATURE = 1600  # CAFT cutoff (K) for LFL/UFL extraction


@dataclass(slots=True)
class ProjectConfig:
    xyz_geom: Path
    mode: str
    tae: float | None
    freqs: list[float] | None
    bond_enthalpy_json: Path
    output_dir: Path
    polys_temps: list[int]
    npoints: int
    lfl_threshold: float
    ufl_threshold: float
    plot_map: bool
    hip_checkpoint: Path | None
    hip_device: str
    skala_device: str


def _resolve_path(value: str | Path | None, fallback: Path) -> Path:
    if value is None:
        return fallback
    return Path(value).expanduser().resolve()


def _load_case_yaml(path: Path) -> tuple[float, list[float]]:
    with open(path, "r", encoding="utf-8") as handle:
        raw: dict[str, Any] = yaml.safe_load(handle) or {}
    if "tae" not in raw:
        raise ValueError(f"{path} is missing required field 'tae' (total atomization energy in Hartree).")
    if "freqs" not in raw:
        raise ValueError(f"{path} is missing required field 'freqs' (vibrational frequencies in cm^-1).")
    return float(raw["tae"]), [float(x) for x in raw["freqs"]]


def build_config(
    *,
    xyz_geom: str | Path,
    mode: str = "ml",
    case_yaml: str | Path | None = None,
    output_dir: str | Path | None = None,
    threshold: float = DEFAULT_THRESHOLD_TEMPERATURE,
    lfl_threshold: float | None = None,
    ufl_threshold: float | None = None,
    plot_map: bool = True,
    hip_checkpoint: str | Path | None = None,
    hip_device: str = "auto",
    skala_device: str = "auto",
) -> ProjectConfig:
    tae: float | None = None
    freqs: list[float] | None = None
    if case_yaml is not None:
        tae, freqs = _load_case_yaml(Path(case_yaml).expanduser().resolve())

    # `threshold` is the shared shorthand; per-limit overrides fall back to it.
    lfl_t = lfl_threshold if lfl_threshold is not None else threshold
    ufl_t = ufl_threshold if ufl_threshold is not None else threshold

    return ProjectConfig(
        xyz_geom=Path(xyz_geom).expanduser().resolve(),
        mode=mode,
        tae=tae,
        freqs=freqs,
        bond_enthalpy_json=DEFAULT_BOND_ENTHALPY_JSON,
        output_dir=_resolve_path(output_dir, DEFAULT_OUTPUT_DIR),
        polys_temps=list(DEFAULT_POLYS_TEMPS),
        npoints=DEFAULT_NPOINTS,
        lfl_threshold=lfl_t,
        ufl_threshold=ufl_t,
        plot_map=plot_map,
        hip_checkpoint=(
            Path(hip_checkpoint).expanduser().resolve()
            if hip_checkpoint is not None
            else None
        ),
        hip_device=hip_device,
        skala_device=skala_device,
    )


def validate_config(cfg: ProjectConfig) -> list[str]:
    errors: list[str] = []
    if not cfg.xyz_geom.exists():
        errors.append(f"Missing required path: {cfg.xyz_geom}")
    if not cfg.bond_enthalpy_json.exists():
        errors.append(f"Missing required path: {cfg.bond_enthalpy_json}")

    if cfg.mode not in {"ab", "ml"}:
        errors.append(f"--mode must be 'ab' or 'ml', got {cfg.mode!r}.")

    for label, value in (("lfl", cfg.lfl_threshold), ("ufl", cfg.ufl_threshold)):
        if value <= 0:
            errors.append(f"--{label}-threshold must be a positive temperature in K, got {value!r}.")

    for label, value in (("skala-device", cfg.skala_device), ("hip-device", cfg.hip_device)):
        if value not in {"auto", "cpu", "cuda"}:
            errors.append(f"--{label} must be 'auto', 'cpu', or 'cuda', got {value!r}.")

    if cfg.mode == "ab":
        if cfg.tae is None or cfg.freqs is None:
            errors.append(
                "--mode ab requires --yaml DATA.yaml with 'tae' (Hartree) and 'freqs' (cm^-1)."
            )
        else:
            if not cfg.freqs:
                errors.append("freqs must not be empty.")
            elif any(freq <= 0 for freq in cfg.freqs):
                errors.append("freqs must contain only positive frequencies in cm^-1.")
    else:
        if cfg.hip_checkpoint is not None and not cfg.hip_checkpoint.exists():
            errors.append(
                f"Missing hip checkpoint: {cfg.hip_checkpoint} "
                "(unset --hip-checkpoint / HIP_CKPT to auto-download from HuggingFace)."
            )

    return errors
