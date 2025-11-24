from __future__ import annotations

from typing import Any, Dict, Optional

from .config import options_from_basicio
from .backend import get_array_module, resolve_dtype
from .beam import gaussian_beam_kernel_pixels, complex_beam_fft, pad_and_wrap_kernel
from .fft_conv import conv2d_fft_real_batched
from .source_cube import fill_tilted_ring_source_cube
from .vel_smooth import apply_velocity_smoothing, smooth_fits_cube_inplace
from .noise import add_gaussian_noise
from .fitsio import write_cube_fits, create_empty_cube_fits, write_cube_chunk

from . import cube_convolution as cc

import numpy as np

def make_cube(
    GalaxyIO: Any,
    DataCube: Any,
    TiltedRing: Any,
    *,
    device: str = "cpu",
    dtype: str = "float32",
    chunk_size: Optional[int] = None,
    threads: Optional[int] = None,
    seed: Optional[int] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Create a mock cube using the pure-Python backend.

    End-to-end CPU pipeline:
      1) Build source cube (particles -> voxel binning)
      2) Create beam kernel and its FFT
      3) Batched 2D FFT-based spatial convolution across channels
      4) Optional velocity smoothing
      5) Write noiseless convolved cube
      6) Add noise and write final cube
    """
    # Resolve runtime options (prefers explicit args; falls back to BasicIO)
    if options is None:
        opts = options_from_basicio(GalaxyIO)
        device = device or opts.device
        dtype = dtype or opts.dtype
        chunk_size = chunk_size if chunk_size is not None else opts.chunk_size
        threads = threads if threads is not None else opts.threads
        seed = seed if seed is not None else opts.seed

    # Resolve backend (M2 enables 'cuda')
    xp, is_gpu = get_array_module(device if device in ('cpu','cuda') else 'cpu')
    np_dtype = resolve_dtype(dtype)

    # 1) Build source cube in NumPy, shape (nz, ny, nx)
    source = fill_tilted_ring_source_cube(DataCube, TiltedRing, seed=seed)
    source = source.astype(np_dtype, copy=False)
    print("UnConvolved Source Cube Flux", np.nansum(source))

    # 2) Beam kernel in pixel units
    # Convert FWHM to sigma in pixels: sigma = FWHM / (2*sqrt(2*ln 2))
    
    ker_fft,pad_shape,BeamArea=cc.get_kernel(DataCube)
    NoiseUse=float(DataCube.noise)*np.sqrt(2*BeamArea)

    import numpy as _np
    """
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
    """
    # 3) Batched spatial convolution (CPU or GPU), with optional chunking along spectral axis
    nz = int(DataCube.cube_shape[2])
    if chunk_size and chunk_size > 0 and not is_gpu:
    
        noiseless_path,final_path=cc.chunked_convolution(GalaxyIO,DataCube,chunk_size,ker_fft,pad_shape,source,seed)
        
        """
        # CPU chunked path
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
        """
        # Apply velocity smoothing in-place on both cubes if requested
        if int(getattr(DataCube, 'velocity_smoothing_switch', 0)) == 1:
            smooth_fits_cube_inplace(
                noiseless_path,
                sigma_kms=getattr(DataCube, 'velocity_smooth_sigma', None),
                channel_kms=float(DataCube.channel_size),
                chunk_size=chunk_size,
            )
            smooth_fits_cube_inplace(
                final_path,
                sigma_kms=getattr(DataCube, 'velocity_smooth_sigma', None),
                channel_kms=float(DataCube.channel_size),
                chunk_size=chunk_size,
            )

        return {
            'noiseless_convolved_cube': noiseless_path,
            'final_cube': final_path,
            'dtype': str(np_dtype),
            'device': 'cpu',
        }
    else:
        # Full in-memory path (CPU or GPU)
        if is_gpu:
            import cupy as cp  # type: ignore
            from .fft_conv_gpu import conv2d_fft_real_batched_gpu
            src_gpu = cp.asarray(source)
            ker_wrapped_host = pad_and_wrap_kernel(kernel, pad_shape)
            ker_wrapped = cp.asarray(ker_wrapped_host, dtype=src_gpu.dtype)
            conv_gpu = conv2d_fft_real_batched_gpu(src_gpu, ker_wrapped, pad_shape)
            convolved = cp.asnumpy(conv_gpu)
        else:
            print("Base Convert")
            convolved = cc.standard_conv(source,ker_fft,pad_shape)
    print("Post Convolution",np.nansum(convolved),np.nansum(source))
    #   Now do a final conversion to go from Jy/pixel to Jy/beam
    convolved=convolved*BeamArea
    print("After unit conversion", np.nansum(convolved),BeamArea)
    # 4) Optional velocity smoothing
    if int(getattr(DataCube, 'velocity_smoothing_switch', 0)) == 1:
        convolved = apply_velocity_smoothing(
            convolved,
            sigma_kms=getattr(DataCube, 'velocity_smooth_sigma', None),
            channel_kms=float(DataCube.channel_size),
        ).astype(convolved.dtype, copy=False)

    # 5) Write noiseless, convolved cube
    noiseless_path = write_cube_fits(convolved, GalaxyIO, DataCube, suffix='ConvolvedSourceCube')

    print(type(ker_fft),type(pad_shape))
    # 6) Add Gaussian noise and write final cube
    final = add_gaussian_noise(convolved, NoiseUse,ker_fft,pad_shape, seed=seed)
    final_path = write_cube_fits(final, GalaxyIO, DataCube, suffix='')

    return {
        'noiseless_convolved_cube': noiseless_path,
        'final_cube': final_path,
        'dtype': str(np_dtype),
        'device': 'cuda' if is_gpu else 'cpu',
    }


