from __future__ import annotations

import argparse

from .config import build_config, validate_config
from .stages import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the flammability pipeline. Default mode 'ml' computes TAE and "
        "frequencies from an XYZ file using skala+hip. Mode 'ab' reads them from a small YAML."
    )
    parser.add_argument("xyz", help="Path to the XYZ geometry file.")
    parser.add_argument(
        "--mode",
        choices=["ml", "ab"],
        default="ml",
        help="ml: compute TAE/freqs from XYZ via skala+hip (default). "
        "ab: read TAE/freqs from --yaml.",
    )
    parser.add_argument(
        "--yaml",
        dest="case_yaml",
        default=None,
        help="Case-data YAML with 'tae' (Hartree) and 'freqs' (cm^-1); required for --mode ab.",
    )
    parser.add_argument("--output-dir", default=None, help="Override output directory.")
    parser.add_argument("--no-plot", action="store_true", help="Skip writing the phase-diagram PDF.")
    parser.add_argument("--hip-checkpoint", default=None, help="Path to hip checkpoint (or set HIP_CKPT).")
    parser.add_argument("--hip-device", default="cpu", help="hip torch device (default: cpu).")
    parser.add_argument("--skala-dir", default=None, help="Path to Skala_TAE source dir (or set SKALA_TAE_DIR).")
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate inputs and exit without running the pipeline.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    cfg = build_config(
        xyz_geom=args.xyz,
        mode=args.mode,
        case_yaml=args.case_yaml,
        output_dir=args.output_dir,
        plot_map=not args.no_plot,
        hip_checkpoint=args.hip_checkpoint,
        hip_device=args.hip_device,
        skala_dir=args.skala_dir,
    )

    errors = validate_config(cfg)
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        raise SystemExit(1)

    if args.validate_only:
        print("Inputs are valid.")
        return

    run_pipeline(cfg)


if __name__ == "__main__":
    main()
