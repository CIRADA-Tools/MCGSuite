from __future__ import annotations

import numpy as np
from typing import Tuple


def gaussian_beam_kernel_pixels(
    sigma_x_pix: float,
    sigma_y_pix: float,
    angle_deg: float,
    n_sigma: float = 5.0,
) -> np.ndarray:
    """Create a normalized 2D Gaussian beam kernel in pixel units.

    This mirrors Calculate2DBeamKernel: rotation uses negative angle to align
    with the beam axes, and normalization ensures sum(kernel) == 1.
    """
    if sigma_x_pix <= 0 or sigma_y_pix <= 0:
        raise ValueError("Beam sigmas must be positive")

    radius = int(np.ceil(n_sigma * max(sigma_x_pix, sigma_y_pix)))
    size = 2 * radius + 1

    y, x = np.mgrid[-radius: radius + 1, -radius: radius + 1]
    theta = -np.deg2rad(angle_deg)
    c, s = np.cos(theta), np.sin(theta)

    xp = x * c - y * s
    yp = x * s + y * c

    r2 = (xp / sigma_x_pix) ** 2 + (yp / sigma_y_pix) ** 2
    kernel = np.exp(-0.5 * r2) / (2.0 * np.pi * sigma_x_pix * sigma_y_pix)
    kernel /= kernel.sum()
    return kernel.astype(np.float32)


def pad_and_wrap_kernel(
    kernel: np.ndarray,
    pad_shape: Tuple[int, int],
) -> np.ndarray:
    """Pad kernel to pad_shape and wrap such that center moves to [0,0].

    This matches MakeWrappedArray: central value moves to index (1,1) in Fortran,
    which corresponds to (0,0) in 0-based Python indexing.
    """
    ny_pad, nx_pad = map(int, pad_shape)
    wrapped = np.zeros((ny_pad, nx_pad), dtype=np.float64)

    ky, kx = kernel.shape
    cy, cx = ky // 2, kx // 2

    for i in range(ky):
        for j in range(kx):
            k = (i - cy) % ny_pad
            l = (j - cx) % nx_pad
            wrapped[k, l] = float(kernel[i, j])

    return wrapped


def complex_beam_fft(
    kernel: np.ndarray,
    pad_shape: Tuple[int, int],
) -> np.ndarray:
    """Compute rfft2 of the wrapped, padded kernel.

    Returns the complex spectrum with shape (ny_pad, nx_pad//2 + 1).
    """
    wrapped = pad_and_wrap_kernel(kernel, pad_shape)
    ker_fft = np.fft.rfft2(wrapped, s=pad_shape, axes=(-2, -1))
    return ker_fft


