from __future__ import annotations

import math
from typing import Any, Tuple

import numpy as np


def _rng(seed: int | None) -> np.random.Generator:
    if seed is None or seed <= 0:
        return np.random.default_rng()
    return np.random.default_rng(int(seed))


def fill_tilted_ring_source_cube(
    DataCube: Any,
    TiltedRing: Any,
    *,
    seed: int | None = None,
) -> np.ndarray:
    """Generate a source cube (pre-convolution, no noise) by sampling particles per ring.

    The implementation mirrors the logic in the Fortran SingleRingGeneration/FillDataCubeByTiltedRing:
    - Determine number of particles per ring proportional to ring area and surface density
    - Sample particle positions in the ring (sqrt radius distribution), azimuth, and vertical z from atanh distribution
    - Project to image plane using inclination and position angle
    - Assign a velocity per particle from rotation, radial/vertical components, dispersion noise
    - Bin into integer pixel/channel indices with nearest rounding
    - Accumulate flux equally split across particles in the ring
    """
    ny, nx, nz = int(DataCube.cube_shape[0]), int(DataCube.cube_shape[1]), int(DataCube.cube_shape[2])

    # Output array shape matches FITS convention used elsewhere: (nz, ny, nx)
    source = np.zeros((nz, ny, nx), dtype=np.float32)

    # Pixel scaling (arcsec per pixel)
    pix_as = float(DataCube.pixel_size)

    # Reference pixel indices (0-based in Python; provided as ints)
    ref_x = int(DataCube.reference_locations[0])
    ref_y = int(DataCube.reference_locations[1])
    ref_ch = int(DataCube.reference_locations[2])

    # Reference physical values
    ref_v = float(DataCube.reference_values[2])  # km/s
    dv = float(DataCube.channel_size)  # km/s per channel

    # Angles to radians
    inc_rad = math.radians(float(TiltedRing.inclination))
    pa_rad = math.radians(float(TiltedRing.position_angle))+np.pi/2.

    cpa, spa = math.cos(pa_rad), math.sin(pa_rad)
    cinc = math.cos(inc_rad)
    sinc = math.sin(inc_rad)

    rng = _rng(seed)
    
    dV=DataCube.channel_size*1000.
    BeamArea=2*np.pi*(DataCube.beam_fwhm*DataCube.beam_fwhm*DataCube.beam_flattening)/2.355**2.

    for i in range(int(TiltedRing.nRings)):
        # Ring parameters in pixel units
        rmid_pix = float(TiltedRing.R_array[i]) / pix_as
        rwidth_pix = float(TiltedRing.Rwidth) / pix_as
        rmin = rmid_pix - 0.5 * rwidth_pix
        rmax = rmid_pix + 0.5 * rwidth_pix
        if rmin < 0:
            rmin = 0.0

        # Surface brightness (Jy km/s arcsec^-2), scale height (arcsec)
        sigma_sb = float(TiltedRing.sigma_array[i])
        #   Convert to Jy/Beam
        #print("INi SB", sigma_sb)
        sigma_sb=sigma_sb/(dV/1000.)*pix_as**2.
        #print("Rad - sigma_b", sigma_sb)
        #   Now get the scale height
        z0_pix = float(TiltedRing.z_scale[i]) / pix_as if getattr(TiltedRing, 'z_scale', None) is not None else 0.0
        # Velocities (km/s)
        vrot = float(TiltedRing.v_tangential_array[i])
        vrad = float(getattr(TiltedRing, 'v_radial', 0.0))
        vvert = float(getattr(TiltedRing, 'v_vertical', 0.0))
        vdisp = float(getattr(TiltedRing, 'v_dispersion', 0.0))
        vsys = float(TiltedRing.vsys)
        dvdz = float(getattr(TiltedRing, 'dvdz', 0.0))

        # Number of particles in ring
        pixel_ring_area = math.pi * (rmax * rmax - rmin * rmin)  # in pixel^2
        cmode = int(getattr(TiltedRing, 'cmode', 0))
        cloud_surf_dens = float(getattr(TiltedRing, 'CloudSurfDens', 100.0))
        dens_mult = cloud_surf_dens * (sigma_sb ** cmode)
        n_particles = int(dens_mult * pixel_ring_area) + 1
        if n_particles <= 0:
            continue

        # Sample cylindrical coordinates
        u_r = rng.random(n_particles)
        rr = np.sqrt(u_r * (rmax * rmax - rmin * rmin) + rmin * rmin)
        theta = rng.random(n_particles) * (2.0 * math.pi)

        if z0_pix > 0:
            # Avoid infs at |u|=1
            u = rng.random(n_particles)
            u = np.clip(u, 1e-6, 1.0 - 1e-6)
            z = np.arctanh(2.0 * u - 1.0) * z0_pix
        else:
            z = np.zeros(n_particles, dtype=np.float64)

        # Cartesian 3D (pixels)
        x = rr * np.cos(theta)
        y = rr * np.sin(theta)

        # Incline about x-axis into major-axis-aligned frame
        y_temp = y * cinc + z * sinc
        x_temp = x

        # Rotate by position angle into observed frame
        x_proj = x_temp * cpa - y_temp * spa
        y_proj = x_temp * spa + y_temp * cpa

        # Shift to cube center (pixel indices)
        x_pix = x_proj + ref_x
        y_pix = y_proj + ref_y

        # Velocities per particle
        cth = np.cos(theta)
        sth = np.sin(theta)
        # Apply vertical gradient beyond threshold height
        vrot_eff = vrot - np.maximum(0.0, np.abs(z) - (5.0 * z0_pix)) * dvdz
        v_from_rot = vrot_eff * cth * sinc
        v_from_rad = vrad * sth * sinc
        v_from_vert = vvert * cinc
        v = vsys + v_from_rot + v_from_rad + v_from_vert
        if vdisp > 0:
            v = v + rng.normal(0.0, vdisp, size=n_particles)

        # Map to discrete indices with nearest rounding
        ix = np.floor(x_pix + 0.5).astype(int)
        iy = np.floor(y_pix + 0.5).astype(int)
        ich = np.floor((v - ref_v) / dv + ref_ch + 0.5).astype(int)

        # Bounds check
        mask = (
            (ix >= 0) & (ix < nx) &
            (iy >= 0) & (iy < ny) &
            (ich >= 0) & (ich < nz)
        )
        if not mask.any():
            continue

        ix = ix[mask]
        iy = iy[mask]
        ich = ich[mask]

        # Flux per particle: Sigma * area / n_particles
        area_pix2 = pixel_ring_area
        flux_pp = (sigma_sb * area_pix2) / float(n_particles)

        # Accumulate
        np.add.at(source, (ich, iy, ix), flux_pp)
    return source


