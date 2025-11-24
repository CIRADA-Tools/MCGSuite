from __future__ import annotations

from typing import Optional

import numpy as np

from . import cube_convolution as cc

def add_gaussian_noise(
    cube: np.ndarray,
    noise_mjy_per_beam: float,
    ker_fft: np.ndarray,
    pad_shape: tuple,
    *,
    seed: Optional[int] = None,
) -> np.ndarray:
    # FITS cubes and code use Jy units; input noise is in mJy/beam
    std_jy = float(noise_mjy_per_beam) / 1000.0
    if std_jy <= 0:
        return cube
    rng = np.random.default_rng(None if seed is None or seed <= 0 else int(seed))
    
    noise=rng.normal(0.0, std_jy, size=cube.shape).astype(cube.dtype, copy=False)
    print("Noise",np.nansum(noise))
    
    noise=cc.standard_conv(noise,ker_fft,pad_shape)
    return cube + noise#rng.normal(0.0, std_jy, size=cube.shape).astype(cube.dtype, copy=False)


