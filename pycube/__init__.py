"""Python implementation of the Mock Cube Generator (MCG) pipeline.

This package will progressively replace the Fortran executable with a
pure-Python implementation supporting both CPU (NumPy) and GPU (CuPy).

Public API (initial):
 - make_cube: entry point mirroring the current MakeCubes flow.
"""

from .make_cube import make_cube  # noqa: F401

__all__ = [
    "make_cube",
]


