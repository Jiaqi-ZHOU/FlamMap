from __future__ import annotations

import sys
from pathlib import Path

from .thermo_extract import parse_mol_struct


def compute_tae(
    xyz_path: Path,
    *,
    skala_dir: Path,
    basis: str = "def2-tzvp",
    xc: str = "skala-1.1",
) -> float:
    sys.path.insert(0, str(skala_dir))
    try:
        from src.energy_calculator import EnergyCalculator
        from src.utils import read_xyz
    finally:
        if str(skala_dir) in sys.path:
            sys.path.remove(str(skala_dir))

    EnergyCalculator.XC = xc
    atoms = read_xyz(str(xyz_path))
    return float(EnergyCalculator(basis=basis).compute_atomization_energy(atoms))


def _apply_hip_torch_patches() -> None:
    import types

    import torch

    if not hasattr(torch.ops, "torch_scatter"):
        torch.ops.torch_scatter = types.SimpleNamespace()

    if not hasattr(torch.ops.torch_scatter, "segment_sum_coo"):
        def segment_sum_coo(src, index, out=None, dim_size=None):
            if dim_size is None:
                dim_size = int(index.max().item()) + 1 if index.numel() > 0 else 0
            if out is None:
                shape = (dim_size,) if src.dim() == 1 else (dim_size, src.size(-1))
                out = torch.zeros(shape, device=src.device, dtype=src.dtype)
            else:
                out.zero_()
            out.index_add_(0, index, src)
            return out

        torch.ops.torch_scatter.segment_sum_coo = segment_sum_coo

    if not hasattr(torch.ops.torch_scatter, "segment_sum_csr"):
        def segment_sum_csr(src, indptr, out=None):
            n = int(indptr.numel() - 1)
            if out is None:
                shape = (n,) if src.dim() == 1 else (n, src.size(-1))
                out = torch.zeros(shape, device=src.device, dtype=src.dtype)
            else:
                out.zero_()
            for i in range(n):
                start = int(indptr[i].item())
                end = int(indptr[i + 1].item())
                if start < end:
                    out[i] = src[start:end].sum(dim=0)
            return out

        torch.ops.torch_scatter.segment_sum_csr = segment_sum_csr

    import torch_geometric.nn.pool as pyg_pool

    def radius_graph_naive(x, r, batch=None, loop=False, max_num_neighbors=32, **_):
        if x.dim() != 2:
            raise ValueError(f"x must be [N, F], got {x.shape}")
        n = x.size(0)
        device = x.device
        if batch is None:
            batch = torch.zeros(n, dtype=torch.long, device=device)
        dist = torch.cdist(x, x)
        same = batch.view(-1, 1).eq(batch.view(1, -1))
        mask = same & (dist <= r)
        if not loop:
            mask.fill_diagonal_(False)
        rows, cols = [], []
        for i in range(n):
            js = torch.nonzero(mask[i], as_tuple=False).view(-1)
            if js.numel() == 0:
                continue
            if max_num_neighbors is not None and js.numel() > max_num_neighbors:
                d = dist[i, js]
                _, idx = torch.topk(d, k=max_num_neighbors, largest=False)
                js = js[idx]
            rows.append(js)
            cols.append(torch.full_like(js, i))
        if not rows:
            return torch.empty((2, 0), dtype=torch.long, device=device)
        return torch.stack([torch.cat(rows), torch.cat(cols)], dim=0)

    pyg_pool.radius_graph = radius_graph_naive


def _hessian_to_3n(hessian, natoms: int):
    import numpy as np
    import torch

    target = 3 * natoms
    if isinstance(hessian, torch.Tensor):
        hessian = hessian.detach().cpu().numpy()
    if hessian.ndim == 2 and hessian.shape == (target, target):
        return hessian
    if hessian.ndim == 4 and hessian.shape == (natoms, 3, natoms, 3):
        return hessian.reshape(target, target)
    if hessian.ndim == 1 and hessian.size == target * target:
        return hessian.reshape(target, target)
    if hessian.ndim == 3 and hessian.shape[1:] == (target, target):
        return hessian[0]
    if hessian.ndim == 5 and hessian.shape[1:] == (natoms, 3, natoms, 3):
        return hessian[0].reshape(target, target)
    raise ValueError(f"Unexpected Hessian shape: {hessian.shape} (natoms={natoms})")


def _eigvals_eV_A2_to_cm1(eigvals):
    import numpy as np
    import scipy.constants as const

    scale = (const.e / (1e-10 * 1e-10)) / const.physical_constants["atomic mass constant"][0]
    omega = np.sqrt(np.abs(eigvals) * scale)
    wn = omega / (2.0 * np.pi * const.c * 100.0)
    return np.sign(eigvals) * wn


def compute_freqs(xyz_path: Path, *, checkpoint_path: Path, device: str = "cpu") -> list[float]:
    _apply_hip_torch_patches()

    import sys

    import numpy as np
    from ase.io import read
    from hip.equiformer_torch_calculator import EquiformerTorchCalculator
    from hip.frequency_analysis import analyze_frequencies_np

    is_linear, *_ = parse_mol_struct(xyz_path)
    atoms = read(str(xyz_path))
    symbols = atoms.get_chemical_symbols()
    natoms = len(symbols)
    pos_ang = atoms.get_positions()

    calc = EquiformerTorchCalculator(
        checkpoint_path=str(checkpoint_path),
        hessian_method="predict",
        device=device,
    )
    results = calc.predict(coords=pos_ang, atomic_nums=atoms.get_atomic_numbers())
    hessian = _hessian_to_3n(results["hessian"], natoms)

    fa = analyze_frequencies_np(hessian, pos_ang, symbols)
    wn_cm1 = _eigvals_eV_A2_to_cm1(fa["eigvals"])

    n_drop = 5 if is_linear else 6
    order = np.argsort(np.abs(wn_cm1))
    vibrational = np.sort(wn_cm1[order[n_drop:]])

    negative = vibrational[vibrational < 0]
    if negative.size > 0:
        formatted = ", ".join(f"{f:.2f}" for f in negative)
        print(
            f"WARNING: hip predicted {negative.size} imaginary mode(s) at {formatted} cm^-1 "
            f"for {xyz_path}; geometry is not a minimum.",
            file=sys.stderr,
            flush=True,
        )
        raise ValueError(
            f"hip predicted negative vibrational frequencies for {xyz_path}; aborting."
        )
    return [float(x) for x in vibrational]
