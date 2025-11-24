from __future__ import annotations

import os
from typing import Any, Dict

import numpy as np
from astropy.io import fits


def ensure_output_folder(folder: str) -> None:
    os.makedirs(folder, exist_ok=True)


def write_cube_fits(
    data_cube: np.ndarray,  # shape (nz, ny, nx)
    GalaxyIO: Any,
    DataCube: Any,
    *,
    suffix: str = "",
) -> str:
    ensure_output_folder(GalaxyIO.GalaxyName)
    name = GalaxyIO.GalaxyName
    if suffix:
        filename = f"{name}/{name}_{suffix}.fits"
    else:
        filename = f"{name}/{name}.fits"

    hdr = fits.Header()
    # Required
    hdr['BITPIX'] = -32
    hdr['NAXIS'] = 3
    hdr['NAXIS1'] = int(DataCube.cube_shape[0])
    hdr['NAXIS2'] = int(DataCube.cube_shape[1])
    hdr['NAXIS3'] = int(DataCube.cube_shape[2])

    # Axis definitions
    # Match Fortran WCS projection and spectral axis naming
    hdr['CTYPE1'] = 'RA---SIN'
    hdr['CTYPE2'] = 'DEC--SIN'
    hdr['CTYPE3'] = 'VELO'
    hdr['CRPIX1'] = int(DataCube.reference_locations[0]) + 1
    hdr['CRPIX2'] = int(DataCube.reference_locations[1]) + 1
    hdr['CRPIX3'] = int(DataCube.reference_locations[2]) + 1
    hdr['CDELT1'] = float(DataCube.cube_dimensions[0]) / 3600.0
    hdr['CDELT2'] = float(DataCube.cube_dimensions[1]) / 3600.0
    hdr['CDELT3'] = float(DataCube.cube_dimensions[2]) * 1000.0
    # reference_values[0,1] are in degrees in the Python pipeline
    hdr['CRVAL1'] = float(DataCube.reference_values[0])
    hdr['CRVAL2'] = float(DataCube.reference_values[1])
    hdr['CRVAL3'] = float(DataCube.reference_values[2]) * 1000.0
    hdr['CUNIT1'] = 'deg'
    hdr['CUNIT2'] = 'deg'
    hdr['CUNIT3'] = 'm/s'

    # Units and meta
    hdr['BUNIT'] = 'Jy/Beam'
    hdr['BTYPE'] = 'intensity'

    # Beam keywords (degrees)
    bmaj = float(DataCube.beam_dimensions[0]) / 3600.0
    bmin = float(DataCube.beam_dimensions[1]) / 3600.0
    bpa = float(DataCube.beam_dimensions[2])
    hdr['BMAJ'] = bmaj
    hdr['BMIN'] = bmin
    hdr['BPA'] = bpa

    # Extra keys to mirror Fortran output
    hdr['RestFreq'] = 1.42040575179e9
    hdr['SPECSYS'] = 'BARYCENT'
    hdr['RADESYS'] = 'FK5'
    hdr['PC01_01'] = 1.0
    hdr['PC02_01'] = 0.0
    hdr['PC03_01'] = 0.0
    hdr['PC01_02'] = 0.0
    hdr['PC02_02'] = 1.0
    hdr['PC03_02'] = 0.0
    hdr['PC01_03'] = 0.0
    hdr['PC02_03'] = 0.0
    hdr['PC03_03'] = 1.0
    hdr['LONPOLE'] = 180.0
    hdr['LATPOLE'] = 0.0

    # Write primary HDU
    hdu = fits.PrimaryHDU(data=data_cube.astype(np.float32, copy=False), header=hdr)
    hdu.writeto(filename, overwrite=True)
    return filename


def create_empty_cube_fits(GalaxyIO: Any, DataCube: Any, *, suffix: str = "") -> str:
    """Create an empty FITS cube with correct header and data shape, filled with zeros.

    Returns the filename. Use write_cube_chunk to fill slices incrementally.
    """
    ensure_output_folder(GalaxyIO.GalaxyName)
    name = GalaxyIO.GalaxyName
    if suffix:
        filename = f"{name}/{name}_{suffix}.fits"
    else:
        filename = f"{name}/{name}.fits"

    hdr = fits.Header()
    hdr['BITPIX'] = -32
    hdr['NAXIS'] = 3
    hdr['NAXIS1'] = int(DataCube.cube_shape[0])
    hdr['NAXIS2'] = int(DataCube.cube_shape[1])
    hdr['NAXIS3'] = int(DataCube.cube_shape[2])
    hdr['CTYPE1'] = 'RA---SIN'
    hdr['CTYPE2'] = 'DEC--SIN'
    hdr['CTYPE3'] = 'VELO'
    hdr['CRPIX1'] = int(DataCube.reference_locations[0]) + 1
    hdr['CRPIX2'] = int(DataCube.reference_locations[1]) + 1
    hdr['CRPIX3'] = int(DataCube.reference_locations[2]) + 1
    hdr['CDELT1'] = float(DataCube.cube_dimensions[0]) / 3600.0
    hdr['CDELT2'] = float(DataCube.cube_dimensions[1]) / 3600.0
    hdr['CDELT3'] = float(DataCube.cube_dimensions[2]) * 1000.0
    hdr['CRVAL1'] = float(DataCube.reference_values[0])
    hdr['CRVAL2'] = float(DataCube.reference_values[1])
    hdr['CRVAL3'] = float(DataCube.reference_values[2]) * 1000.0
    hdr['CUNIT1'] = 'deg'
    hdr['CUNIT2'] = 'deg'
    hdr['CUNIT3'] = 'm/s'
    hdr['BUNIT'] = 'Jy/Beam'
    hdr['BTYPE'] = 'intensity'
    hdr['RestFreq'] = 1.42040575179e9
    hdr['SPECSYS'] = 'BARYCENT'
    hdr['RADESYS'] = 'FK5'
    hdr['PC01_01'] = 1.0
    hdr['PC02_01'] = 0.0
    hdr['PC03_01'] = 0.0
    hdr['PC01_02'] = 0.0
    hdr['PC02_02'] = 1.0
    hdr['PC03_02'] = 0.0
    hdr['PC01_03'] = 0.0
    hdr['PC02_03'] = 0.0
    hdr['PC03_03'] = 1.0
    hdr['LONPOLE'] = 180.0
    hdr['LATPOLE'] = 0.0
    hdr['BMAJ'] = float(DataCube.beam_dimensions[0]) / 3600.0
    hdr['BMIN'] = float(DataCube.beam_dimensions[1]) / 3600.0
    hdr['BPA'] = float(DataCube.beam_dimensions[2])

    nz = int(DataCube.cube_shape[2])
    ny = int(DataCube.cube_shape[0])
    nx = int(DataCube.cube_shape[1])
    data = np.zeros((nz, ny, nx), dtype=np.float32)
    fits.PrimaryHDU(data=data, header=hdr).writeto(filename, overwrite=True)
    return filename


def write_cube_chunk(filename: str, z_start: int, z_end: int, chunk_data: np.ndarray) -> None:
    """Write a chunk [z_start:z_end] into an existing FITS cube file in-place."""
    with fits.open(filename, mode='update', memmap=True) as hdul:
        hdul[0].data[z_start:z_end, :, :] = chunk_data.astype(np.float32, copy=False)
        hdul.flush()


