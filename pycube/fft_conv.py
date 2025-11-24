from __future__ import annotations

from typing import Tuple, Union

import numpy as np


def conv2d_fft_real(
    real_image: np.ndarray,
    kernel_fft: np.ndarray,
    pad_shape: Tuple[int, int],
) -> np.ndarray:
    """Convolve 2D real_image with a kernel represented by its rfft2.

    Matches the Fortran flow:
    - zero-pad real input to pad_shape
    - forward rfft2
    - multiply by kernel_fft
    - inverse irfft2
    - normalize by pad_size product
    - crop back to input shape
    """
    ny, nx = real_image.shape
    ny_pad, nx_pad = map(int, pad_shape)

    # pad input
    padded = np.zeros((ny_pad, nx_pad), dtype=np.float64)
    padded[:ny, :nx] = real_image.astype(np.float64, copy=False)

    img_fft = np.fft.rfft2(padded, s=pad_shape, axes=(-2, -1))
    out_fft = img_fft * kernel_fft
    real_conv = np.fft.irfft2(out_fft, s=pad_shape, axes=(-2, -1))
    #real_conv /= (ny_pad * nx_pad)

    return real_conv[:ny, :nx].astype(real_image.dtype, copy=False)


def conv2d_fft_real_batched(
    real_stack: np.ndarray,
    kernel_fft: Union[np.ndarray, np.ndarray],
    pad_shape: Tuple[int, int],
) -> np.ndarray:
    """Batched 2D convolution for shape (nz, ny, nx). Kernel FFT can broadcast.

    Performs a single batched rfft2/irfft2 across the last two dimensions.
    """
    nz, ny, nx = real_stack.shape
    ny_pad, nx_pad = map(int, pad_shape)

    padded = np.zeros((nz, ny_pad, nx_pad), dtype=np.float64)
    padded[:, :ny, :nx] = real_stack.astype(np.float64, copy=False)

    img_fft = np.fft.rfft2(padded, s=pad_shape, axes=(-2, -1))
    out_fft = img_fft * kernel_fft  # broadcast over batch
    real_conv = np.fft.irfft2(out_fft, s=pad_shape, axes=(-2, -1))
    print("After numpy fft", np.nansum(real_conv))
    #real_conv /= (ny_pad * nx_pad)


    return real_conv[:, :ny, :nx].astype(real_stack.dtype, copy=False)


