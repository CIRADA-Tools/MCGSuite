from __future__ import annotations

from typing import Literal, Tuple


def resolve_dtype(dtype: str):
    import numpy as _np

    dtype = dtype.lower()
    if dtype in ("float32", "fp32"):  # pragma: no cover - mapping
        return _np.float32
    if dtype in ("float64", "fp64", "double"):
        return _np.float64
    raise ValueError(f"Unsupported dtype: {dtype}")


def get_array_module(device: Literal["cpu", "cuda"]) -> Tuple[object, bool]:
    """Return the array module (NumPy or CuPy) and a boolean indicating GPU use."""
    if device == "cuda":
        try:
            import cupy as cp  # type: ignore

            return cp, True
        except Exception as exc:  # pragma: no cover - runtime env dependent
            raise RuntimeError(
                "Device 'cuda' requested but CuPy is not available.\n"
                "Install a matching CuPy wheel, e.g. 'pip install cupy-cuda12x'."
            ) from exc
    else:
        import numpy as np

        return np, False


