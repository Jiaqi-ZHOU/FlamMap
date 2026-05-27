from __future__ import annotations

import json
import os
import warnings
from pathlib import Path


def _silence_benign_warnings() -> None:
    # gpu4pyscf prints this whenever cutensor isn't installed alongside cupy;
    # it just means contractions fall back to cupy. Harmless, but fires even
    # on CPU runs because skala's backend.py imports gpu4pyscf in a try/except.
    warnings.filterwarnings(
        "ignore",
        message=r"using \S+ as the tensor contraction engine\.",
        category=UserWarning,
    )
    # torch-sparse 0.6.18+pt25cu121 has a Python wrapper whose Optional[Tensor]
    # annotation disagrees with the C++ op's Tensor signature, so PyG's
    # JIT-compile probe fails and torch-sparse gets disabled at import time.
    # Hip's Hessian path doesn't go through torch-sparse, so the warning is
    # noise — silence it. (Fix has to come from upstream torch-scatter.)
    warnings.filterwarnings(
        "ignore",
        message=r"An issue occurred while importing 'torch-sparse'.*",
        category=UserWarning,
    )
    # hf_xet (the Rust-based Xet transfer backend bundled with huggingface_hub
    # 1.x) writes "Warning: You are sending unauthenticated requests..." plus
    # its own tracing events directly to fd 2 from native code — Python's
    # warnings filter and sys.stderr override can't catch them. Disabling the
    # Xet backend (HF_HUB_DISABLE_XET=1) also doesn't silence the warning in
    # huggingface_hub 1.16; the only reliable knob is `_silenced_fd2` below,
    # used around `hf_hub_download`. Leaving Xet off is still desirable as a
    # belt-and-braces measure — the hip checkpoint is a one-shot 190 MB
    # download where plain HTTP is fine.
    os.environ.setdefault("HF_HUB_DISABLE_XET", "1")


_silence_benign_warnings()


import contextlib


@contextlib.contextmanager
def _silenced_fd2():
    """Redirect file descriptor 2 (OS-level stderr) to /dev/null.

    Needed because hf_xet's Rust code writes "Warning: unauthenticated…" and
    tracing events directly to fd 2, bypassing `sys.stderr`. Python-level
    `warnings.filterwarnings` and `contextlib.redirect_stderr` can't catch
    those — only an fd-level redirect works.

    Real failures from `hf_hub_download` surface as Python exceptions, so
    suppressing fd 2 chatter during the call doesn't hide anything important.
    """
    devnull = os.open(os.devnull, os.O_WRONLY)
    saved = os.dup(2)
    try:
        os.dup2(devnull, 2)
        yield
    finally:
        os.dup2(saved, 2)
        os.close(devnull)
        os.close(saved)


_ATOMIC_SPINS = {"H": 1, "C": 2, "N": 3, "O": 2}
_ALLOWED_ELEMENTS = set(_ATOMIC_SPINS)

# Heuristic: below this many heavy atoms (non-H), CPU is faster wall-time
# than GPU because of the ~10–15 s cupy / gpu4pyscf init overhead.
AUTO_DEVICE_HEAVY_ATOM_THRESHOLD = 10


def _count_heavy_atoms(xyz_path: Path) -> int:
    # Use a lightweight parser instead of ase.io.read to avoid pulling heavy
    # deps just for atom counting in the auto-device decision.
    with xyz_path.open(encoding="utf-8") as fh:
        lines = [ln.strip() for ln in fh if ln.strip()]
    if len(lines) < 3:
        return 0
    return sum(1 for ln in lines[2:] if ln.split() and ln.split()[0] != "H")


def resolve_device_for_xyz(device_arg: str, xyz_path: Path) -> str:
    """Resolve 'auto'/'cpu'/'cuda' against the molecule size.

    'auto' picks 'cuda' when the molecule has more heavy atoms than the
    threshold *and* CUDA is actually available; otherwise 'cpu'.
    """
    if device_arg in ("cpu", "cuda"):
        return device_arg
    if device_arg != "auto":
        raise ValueError(f"Unsupported device: {device_arg!r}")

    try:
        import torch

        cuda_ok = torch.cuda.is_available()
    except Exception:
        cuda_ok = False
    if not cuda_ok:
        return "cpu"

    n_heavy = _count_heavy_atoms(xyz_path)
    return "cuda" if n_heavy > AUTO_DEVICE_HEAVY_ATOM_THRESHOLD else "cpu"

_CUDA_LIBS_PRELOADED = False


def _preload_cuda_libs() -> None:
    # cupy / gpu4pyscf dlopen libcusolver.so.11 etc. The nvidia-*-cu12 wheels
    # ship them under site-packages/nvidia/<name>/lib/ but that directory is
    # not on the loader search path, so we preload them explicitly with
    # RTLD_GLOBAL before the first cupy import. We skip libnvblas (a BLAS
    # interceptor that breaks numpy/torch BLAS when present in the process).
    global _CUDA_LIBS_PRELOADED
    if _CUDA_LIBS_PRELOADED:
        return
    import ctypes

    try:
        import nvidia  # type: ignore[import-not-found]
    except ImportError:
        _CUDA_LIBS_PRELOADED = True
        return

    nvidia_root = Path(nvidia.__path__[0])
    for name in (
        "cuda_runtime",
        "nvjitlink",
        "cuda_nvrtc",
        "cublas",
        "cusparse",
        "cusolver",
        "cufft",
        "curand",
    ):
        lib_dir = nvidia_root / name / "lib"
        if not lib_dir.is_dir():
            continue
        for so in sorted(lib_dir.glob("lib*.so.*")):
            if "nvblas" in so.name.lower():
                continue
            try:
                ctypes.CDLL(str(so), mode=ctypes.RTLD_GLOBAL)
            except OSError:
                pass
    _CUDA_LIBS_PRELOADED = True


def _atomic_refs_cache_path(xc: str, basis: str) -> Path:
    cache_dir = Path(
        os.environ.get("FLAMMAP_CACHE", Path.home() / ".cache" / "flammap")
    )
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir / f"atomic_refs_{xc}_{basis.replace('/', '_')}.json"


def _skala_energy(
    atom: str, *, basis: str, xc: str, device, spin: int = 0, charge: int = 0
) -> float:
    from pyscf import gto

    mol = gto.M(
        atom=atom,
        basis=basis,
        charge=charge,
        spin=spin,
        unit="Angstrom",
        verbose=0,
    )
    if device.type == "cuda":
        # skala.pyscf.SkalaKS with device=cuda mixes pyscf grids into
        # gpu4pyscf numint and breaks. The supported GPU path is
        # skala.gpu4pyscf.SkalaKS, which extends gpu4pyscf classes end-to-end.
        from skala.gpu4pyscf import SkalaKS as SkalaKS_gpu

        ks = SkalaKS_gpu(mol, xc=xc)
    else:
        from skala.pyscf import SkalaKS

        ks = SkalaKS(mol, xc=xc, device=device)
    return float(ks.kernel())


def _atomic_reference_energies(
    *, basis: str, xc: str, device
) -> dict[str, float]:
    cache = _atomic_refs_cache_path(xc, basis)
    if cache.exists():
        with cache.open(encoding="utf-8") as fh:
            cached = json.load(fh)
        if set(cached) >= _ALLOWED_ELEMENTS:
            return {el: float(cached[el]) for el in _ATOMIC_SPINS}

    energies: dict[str, float] = {}
    for element, spin in _ATOMIC_SPINS.items():
        energies[element] = _skala_energy(
            f"{element} 0 0 0", basis=basis, xc=xc, device=device, spin=spin
        )

    with cache.open("w", encoding="utf-8") as fh:
        json.dump(energies, fh, indent=2)
    return energies


def _read_xyz_atoms(path: Path) -> list[tuple[str, float, float, float]]:
    from ase.io import read

    atoms_obj = read(str(path))
    coords: list[tuple[str, float, float, float]] = []
    for sym, (x, y, z) in zip(
        atoms_obj.get_chemical_symbols(), atoms_obj.get_positions()
    ):
        if sym not in _ALLOWED_ELEMENTS:
            raise ValueError(
                f"Unsupported element {sym!r}. Only H, C, N, O are allowed."
            )
        coords.append((sym, float(x), float(y), float(z)))
    return coords


def compute_tae(
    xyz_path: Path,
    *,
    basis: str = "def2-tzvp",
    xc: str = "skala-1.0",
    device: str = "cpu",
) -> float:
    from collections import Counter

    import torch

    torch_device = torch.device(device)
    if torch_device.type == "cuda":
        _preload_cuda_libs()

    atoms = _read_xyz_atoms(xyz_path)
    atom_block = "\n".join(f"{sym} {x} {y} {z}" for sym, x, y, z in atoms)

    atomic_energies = _atomic_reference_energies(
        basis=basis, xc=xc, device=torch_device
    )
    molecular_energy = _skala_energy(
        atom_block, basis=basis, xc=xc, device=torch_device
    )

    counts = Counter(symbol for symbol, *_ in atoms)
    separated = sum(
        counts.get(el, 0) * atomic_energies[el] for el in _ALLOWED_ELEMENTS
    )
    tae = separated - molecular_energy

    print(f"   Element reference energies (xc={xc}, basis={basis}, device={device}):")
    for el in ("H", "C", "N", "O"):
        if el in atomic_energies:
            print(
                f"     {el}: {atomic_energies[el]:.6f} Ha   (count in molecule: {counts.get(el, 0)})"
            )
    print(f"   Molecular energy: {molecular_energy:.6f} Ha")
    print(
        f"   TAE = sum(n_i * E_atom_i) - E_mol "
        f"= {separated:.6f} - ({molecular_energy:.6f}) = {tae:.6f} Ha"
    )

    return float(tae)


def _apply_hip_torch_patches() -> None:
    import os
    import types

    import torch

    os.environ.setdefault("HIP_STRICT_PATHS", "0")

    if not getattr(torch.load, "_flammap_unsafe_patched", False):
        _orig_load = torch.load

        def _trusted_load(*args, **kwargs):
            kwargs["weights_only"] = False
            return _orig_load(*args, **kwargs)

        _trusted_load._flammap_unsafe_patched = True
        torch.load = _trusted_load

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

    import torch_geometric.nn as pyg_nn

    pyg_nn.radius_graph = radius_graph_naive


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

    scale = (const.e / (1e-10 * 1e-10)) / const.physical_constants[
        "atomic mass constant"
    ][0]
    omega = np.sqrt(np.abs(eigvals) * scale)
    wn = omega / (2.0 * np.pi * const.c * 100.0)
    return np.sign(eigvals) * wn


HIP_HF_REPO = "andreasburger/hip"
HIP_HF_FILENAME = "ckpt/hip_v3_cf.ckpt"

# hip's pyproject.toml lacks `package-data` / `include-package-data`, so
# setuptools silently drops these non-Python files when uv builds hip from
# a git source. Keep this rev in sync with the `hip` entry in pyproject.toml.
HIP_GIT_REPO = "BurgerAndreas/hip"
HIP_GIT_REV = "06fedad91a558ed1589c79a13bc2b5dca7f733de"
_HIP_DATA_FILES = (
    "nets/equiformer_v2/Jd.pt",
    "ocpmodels/models/scn/Jd.pt",
)


def _patch_hip_config_path() -> None:
    """Redirect hip's `find_project_root` to a vendored config dir.

    hip's `training_module.PotentialModule.__init__` calls
    `find_project_root()` and then reads `configs/equiformer_v2.yaml`
    relative to it. The function walks up from `hip/path_config.py`
    looking for `pyproject.toml` or `.git`, which works for editable
    installs but lands at unpredictable directories for wheel installs.
    We ship the required yaml inside the flammability package and patch
    the lookup to return our data dir.
    """
    import hip.training_module as hip_tm

    hip_data_dir = Path(__file__).parent / "_hip_data"

    def _flammap_find_project_root(
        start_path=None, markers=("pyproject.toml", ".git")
    ) -> str:
        return str(hip_data_dir)

    hip_tm.find_project_root = _flammap_find_project_root


def _ensure_hip_packaged_data() -> None:
    """Backfill hip's missing data files in the installed wheel.

    Editable / path installs already have these files next to the source.
    Wheel installs from git lack them due to an upstream packaging gap; we
    fetch them once from raw.githubusercontent.com at the pinned rev.

    We can't use `importlib.util.find_spec` to locate `nets` because
    `nets/__init__.py` eagerly imports submodules that load these files,
    so triggering the import is exactly what we're trying to avoid. Instead
    we look directly in the site-packages dir; if `nets` lives elsewhere
    (e.g. an editable install pointing at a source checkout), the file is
    already next to the source and no backfill is needed.
    """
    import sysconfig
    import urllib.request

    purelib = Path(sysconfig.get_paths()["purelib"])
    for rel in _HIP_DATA_FILES:
        target = purelib / rel
        # Skip if the surrounding package isn't installed under purelib (e.g.
        # editable install pointing at a source checkout that already has the
        # data file). We detect "installed here" by the presence of any .py
        # file in the target's parent dir — `nets/equiformer_v2/` is a
        # namespace package so we can't rely on __init__.py.
        if not any(target.parent.glob("*.py")):
            continue
        if target.is_file():
            continue
        url = f"https://raw.githubusercontent.com/{HIP_GIT_REPO}/{HIP_GIT_REV}/{rel}"
        print(f"   Fetching missing {rel} from hip @ {HIP_GIT_REV[:7]} ...")
        target.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(url, timeout=30) as r, target.open("wb") as f:
            f.write(r.read())


def resolve_hip_checkpoint(checkpoint_path: Path | None) -> Path:
    """Return a usable hip checkpoint path.

    Resolution order:
      1. ``checkpoint_path`` argument (typically from --hip-checkpoint).
      2. ``$HIP_CKPT`` environment variable.
      3. Download ``hip_v3.ckpt`` from HuggingFace (cached under
         ``$HF_HOME``/``~/.cache/huggingface``).
    """
    if checkpoint_path is not None:
        p = Path(checkpoint_path).expanduser()
        if not p.is_file():
            raise FileNotFoundError(f"hip checkpoint not found: {p}")
        return p

    env = os.environ.get("HIP_CKPT")
    if env:
        p = Path(env).expanduser()
        if not p.is_file():
            raise FileNotFoundError(f"HIP_CKPT points at a missing file: {p}")
        return p

    from huggingface_hub import hf_hub_download

    print(
        f"   No --hip-checkpoint / HIP_CKPT set; fetching {HIP_HF_FILENAME} "
        f"from huggingface.co/{HIP_HF_REPO} (cached in HF hub cache)..."
    )
    with _silenced_fd2():
        cached = hf_hub_download(repo_id=HIP_HF_REPO, filename=HIP_HF_FILENAME)
    return Path(cached)


def compute_freqs(
    xyz_path: Path, *, checkpoint_path: Path | None = None, device: str = "cpu"
) -> list[float]:
    _ensure_hip_packaged_data()
    _apply_hip_torch_patches()

    import sys

    import numpy as np
    from ase.io import read
    from hip.equiformer_torch_calculator import EquiformerTorchCalculator
    from hip.frequency_analysis import analyze_frequencies_np

    _patch_hip_config_path()
    resolved_ckpt = resolve_hip_checkpoint(checkpoint_path)

    atoms = read(str(xyz_path))
    symbols = atoms.get_chemical_symbols()
    natoms = len(symbols)
    pos_ang = atoms.get_positions()

    calc = EquiformerTorchCalculator(
        checkpoint_path=str(resolved_ckpt),
        hessian_method="predict",
        device=device,
    )
    results = calc.predict(coords=pos_ang, atomic_nums=atoms.get_atomic_numbers())
    hessian = _hessian_to_3n(results["hessian"], natoms)

    fa = analyze_frequencies_np(hessian, pos_ang, symbols)
    wn_cm1 = _eigvals_eV_A2_to_cm1(fa["eigvals"])
    vibrational = np.sort(wn_cm1)

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
