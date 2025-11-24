from __future__ import annotations

from typing import Tuple


def conv2d_fft_real_batched_gpu(
    real_stack,  # cupy.ndarray, shape (nz, ny, nx)
    kernel_wrapped,  # cupy.ndarray, shape (ny_pad, nx_pad)
    pad_shape: Tuple[int, int],
):
    """Batched 2D convolution on GPU using cuFFT via CuPy.

    Steps:
      - zero-pad each channel to pad_shape
      - rfft2 batched across last two axes
      - multiply with kernel spectrum (broadcasted)
      - irfft2 back and normalize by pad_size product
      - crop to original ny, nx
    Returns cupy.ndarray (nz, ny, nx)
    """
    import cupy as cp  # type: ignore

    nz, ny, nx = real_stack.shape
    ny_pad, nx_pad = map(int, pad_shape)

    padded = cp.zeros((nz, ny_pad, nx_pad), dtype=real_stack.dtype)
    padded[:, :ny, :nx] = real_stack

    img_fft = cp.fft.rfft2(padded, s=pad_shape, axes=(-2, -1))
    ker_fft = cp.fft.rfft2(kernel_wrapped, s=pad_shape, axes=(-2, -1))
    out_fft = img_fft * ker_fft  # broadcast over batch
    real_conv = cp.fft.irfft2(out_fft, s=pad_shape, axes=(-2, -1))
    #real_conv /= (ny_pad * nx_pad)

    return real_conv[:, :ny, :nx]


