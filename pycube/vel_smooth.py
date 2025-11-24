from __future__ import annotations

from typing import Optional

import numpy as np
from scipy.ndimage import gaussian_filter1d


def apply_velocity_smoothing(
    cube: np.ndarray,
    sigma_kms: Optional[float],
    channel_kms: float,
) -> np.ndarray:
    if sigma_kms is None or sigma_kms <= 0:
        return cube
    sigma_ch = float(sigma_kms) / float(channel_kms)
    if sigma_ch <= 0:
        return cube
    # Apply along spectral axis (axis=0 for (nz, ny, nx))
    return gaussian_filter1d(cube, sigma=sigma_ch, axis=0, mode='nearest')


def smooth_fits_cube_inplace(
    filename: str,
    *,
    sigma_kms: Optional[float],
    channel_kms: float,
    chunk_size: int,
) -> None:
    """Apply velocity smoothing along spectral axis to a FITS cube on disk.

    Processes data in chunks with overlap so memory stays bounded.
    """
    if sigma_kms is None or sigma_kms <= 0:
        return
    sigma_ch = float(sigma_kms) / float(channel_kms)
    if sigma_ch <= 0:
        return

    from astropy.io import fits

    margin = int(np.ceil(3.0 * sigma_ch))  # 3-sigma overlap
    with fits.open(filename, mode='update', memmap=True) as hdul:
        data = hdul[0].data  # shape (nz, ny, nx)
        nz = data.shape[0]
        for z0 in range(0, nz, chunk_size):
            z1 = min(nz, z0 + chunk_size)
            a = max(0, z0 - margin)
            b = min(nz, z1 + margin)
            window = np.array(data[a:b, :, :], copy=False)
            # Smooth the window fully, then extract the center region
            smoothed_window = gaussian_filter1d(window, sigma=sigma_ch, axis=0, mode='nearest')
            out_start = z0 - a
            out_end = out_start + (z1 - z0)
            data[z0:z1, :, :] = smoothed_window[out_start:out_end, :, :]
        hdul.flush()


