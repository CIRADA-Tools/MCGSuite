from __future__ import annotations

from typing import Any, Dict, Optional

from .config import options_from_basicio
from .backend import get_array_module, resolve_dtype
from .beam import gaussian_beam_kernel_pixels, complex_beam_fft, pad_and_wrap_kernel
from .fft_conv import conv2d_fft_real_batched
from .fitsio import write_cube_fits, create_empty_cube_fits, write_cube_chunk

import numpy as np

def get_kernel(
    DataCube: Any):

    import numpy as _np
    
    print("calculating Kernel")
    
    fwhm_major_as = float(DataCube.beam_dimensions[0])
    fwhm_minor_as = float(DataCube.beam_dimensions[1])
    angle_deg = float(DataCube.beam_dimensions[2])
    pix_as = float(DataCube.pixel_size)
    to_sigma = 1.0 / (2.0 * _np.sqrt(2.0 * _np.log(2.0)))
    sigma_x_pix = (fwhm_major_as / pix_as) * to_sigma
    sigma_y_pix = (fwhm_minor_as / pix_as) * to_sigma
    BeamArea=2.*np.pi*sigma_x_pix*sigma_y_pix
    kernel = gaussian_beam_kernel_pixels(sigma_x_pix, sigma_y_pix, angle_deg, n_sigma=float(DataCube.nSigma))
    print("Ini kernel total", np.nansum(kernel))
    # Pad shape and kernel FFT (CPU: compute spectrum; GPU path computes on-device)
    ny, nx = int(DataCube.cube_shape[0]), int(DataCube.cube_shape[1])
    pad_shape = (int(_np.ceil(ny + kernel.shape[0])), int(_np.ceil(nx + kernel.shape[1])))
    ker_fft = complex_beam_fft(kernel, pad_shape)
    return ker_fft,pad_shape,BeamArea


def chunked_convolution(GalaxyIO,DataCube,chunk_size,ker_fft,pad_shape,source,seed):
    
    nz = int(DataCube.cube_shape[2])
    noiseless_path = create_empty_cube_fits(GalaxyIO, DataCube, suffix='ConvolvedSourceCube')
    final_path = create_empty_cube_fits(GalaxyIO, DataCube, suffix='')

    for z0 in range(0, nz, chunk_size):
        z1 = min(nz, z0 + chunk_size)
        src_chunk = source[z0:z1, :, :]
        conv_chunk = conv2d_fft_real_batched(src_chunk, ker_fft, pad_shape)
        write_cube_chunk(noiseless_path, z0, z1, conv_chunk)
        # Add noise on the fly for final
        noisy_chunk = add_gaussian_noise(conv_chunk, float(DataCube.noise), seed=seed)
        write_cube_chunk(final_path, z0, z1, noisy_chunk)


    return noiseless_path,final_path


def standard_conv(source,ker_fft,pad_shape):

    convolved = conv2d_fft_real_batched(source, ker_fft, pad_shape)
    
    return convolved
