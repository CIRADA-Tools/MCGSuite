#!/usr/bin/env python3
import argparse
import os
import shutil
import subprocess
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
from astropy.io import fits

# Ensure project root is on sys.path so we can import make_galaxy_code
HERE = os.path.abspath(os.path.dirname(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

# Import to discover OutputFolder from config
from make_galaxy_code import Inputs as IN


def list_dirs(path: str) -> List[str]:
    if not os.path.isdir(path):
        return []
    return [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]


def run_and_capture(cmd: List[str], env_overrides: Optional[Dict[str, str]] = None) -> int:
    env = os.environ.copy()
    if env_overrides:
        env.update(env_overrides)
    print("Running:", " ".join(cmd))
    try:
        res = subprocess.run(cmd, check=True, env=env)
        return res.returncode
    except subprocess.CalledProcessError as e:
        print("Command failed:", e)
        return e.returncode if e.returncode is not None else 1


def detect_new_dir(base: str, before: List[str]) -> Optional[str]:
    after = list_dirs(base)
    new = [d for d in after if d not in before]
    if len(new) == 1:
        return os.path.join(base, new[0])
    if len(new) > 1:
        # Pick the most recently modified among new ones
        new_sorted = sorted(new, key=lambda d: os.path.getmtime(os.path.join(base, d)), reverse=True)
        return os.path.join(base, new_sorted[0])
    # Fallback: pick the most recently modified directory in base
    if after:
        after_sorted = sorted(after, key=lambda d: os.path.getmtime(os.path.join(base, d)), reverse=True)
        return os.path.join(base, after_sorted[0])
    return None


def rename_output_folder(found_path: str, suffix: str) -> str:
    base_dir = os.path.dirname(found_path)
    new_path = found_path + f"_{suffix}"
    if os.path.exists(new_path):
        shutil.rmtree(new_path, ignore_errors=True)
    shutil.move(found_path, new_path)
    return new_path


def find_cube_files(run_dir: str) -> Tuple[str, str]:
    files = [f for f in os.listdir(run_dir) if f.lower().endswith('.fits')]
    noiseless = None
    final = None
    for f in files:
        if f.endswith('_ConvolvedSourceCube.fits'):
            noiseless = os.path.join(run_dir, f)
        elif f.endswith('.fits'):
            final = os.path.join(run_dir, f)
    if noiseless is None or final is None:
        # Fallback: try to detect by substrings
        for f in files:
            if 'Convolved' in f:
                noiseless = noiseless or os.path.join(run_dir, f)
            else:
                final = final or os.path.join(run_dir, f)
    if noiseless is None or final is None:
        raise FileNotFoundError(f"Could not locate cube files in {run_dir}")
    return noiseless, final


def header_compare(h1: fits.Header, h2: fits.Header) -> List[str]:
    keys = [
        'NAXIS','NAXIS1','NAXIS2','NAXIS3',
        'CTYPE1','CTYPE2','CTYPE3',
        'CRPIX1','CRPIX2','CRPIX3',
        'CDELT1','CDELT2','CDELT3',
        'CRVAL1','CRVAL2','CRVAL3',
        'CUNIT1','CUNIT2','CUNIT3',
        'BUNIT','BTYPE',
        'BMAJ','BMIN','BPA'
    ]
    diffs = []
    for k in keys:
        v1 = h1.get(k)
        v2 = h2.get(k)
        if isinstance(v1, float) or isinstance(v2, float):
            if v1 is None or v2 is None:
                diffs.append(f"{k}: {v1} vs {v2}")
            else:
                if not np.isfinite(v1) or not np.isfinite(v2) or abs(v1 - v2) > 1e-6:
                    diffs.append(f"{k}: {v1} vs {v2}")
        else:
            if v1 != v2:
                diffs.append(f"{k}: {v1} vs {v2}")
    return diffs


def data_report(a_path: str, b_path: str, label: str) -> None:
    a = fits.getdata(a_path).astype(np.float64)
    b = fits.getdata(b_path).astype(np.float64)
    ha = fits.getheader(a_path)
    hb = fits.getheader(b_path)
    print(f"[{label}] shapes: {a.shape} vs {b.shape}")
    diffs = header_compare(ha, hb)
    if diffs:
        print(f"[{label}] header diffs:")
        for d in diffs:
            print("  ", d)
    d = a - b
    with np.errstate(invalid='ignore'):
        max_abs = np.nanmax(np.abs(d))
        rms = np.sqrt(np.nanmean(d**2))
        sra = np.nansum(a)
        srb = np.nansum(b)
        ratio = (sra / srb) if srb != 0 else np.nan
    print(f"[{label}] max_abs={max_abs:.6g} rms={rms:.6g} sum_ratio={ratio:.6g}")


def main():
    parser = argparse.ArgumentParser(description="Run and compare Fortran vs Python (CPU/GPU) outputs")
    parser.add_argument("--compare", choices=["cpu","gpu","both"], default="both")
    parser.add_argument("--dtype", choices=["float32","float64"], default="float32")
    parser.add_argument("--chunk-size", type=int, default=0)
    parser.add_argument("--keep-temp", action="store_true", help="Keep intermediate output folders")
    args = parser.parse_args()

    # Discover base output folder from config
    Galaxy, DataCube, TiltedRing, Profiles, GalaxyIO = IN.GetMakeGalaxyInputs()
    out_base = GalaxyIO.OutputFolder
    os.makedirs(out_base, exist_ok=True)

    # Prepare paths before runs
    before_dirs = set(list_dirs(out_base))

    # FORTRAN baseline
    print("=== Running Fortran baseline ===")
    env_no_py = os.environ.copy()
    env_no_py.pop('MCG_BACKEND', None)
    rc = run_and_capture([sys.executable, 'make_galaxy_MCG.py'], env_overrides=env_no_py)
    if rc != 0:
        print("Fortran run failed; aborting compare")
        sys.exit(1)
    # Attempt to detect new folder by comparing before/after and by inferring name
    f_dir = detect_new_dir(out_base, list(before_dirs))
    if not f_dir:
        # Try to infer from GalaxyIO.GalaxyName inside the Fortran path
        # The Fortran run moves the galaxy folder into OutputFolder
        # We'll search for a folder starting with 'ba_' or matching patterns
        candidates = [d for d in list_dirs(out_base) if d.startswith('ba_')]
        if candidates:
            candidates_sorted = sorted(candidates, key=lambda d: os.path.getmtime(os.path.join(out_base, d)), reverse=True)
            f_dir = os.path.join(out_base, candidates_sorted[0])
    if not f_dir:
        print("Could not detect Fortran output folder; ensure the Fortran run produced output inside", out_base)
        sys.exit(1)
    f_dir = rename_output_folder(f_dir, 'fortran')
    before_dirs.add(os.path.basename(f_dir))

    # PY CPU
    cpu_dir = None
    if args.compare in ("cpu","both"):
        print("=== Running Python CPU ===")
        rc = run_and_capture([sys.executable, 'make_galaxy_MCG.py', '--backend','python','--device','cpu','--dtype', args.dtype] + ([] if args.chunk_size<=0 else ['--chunk-size', str(args.chunk_size)]))
        if rc == 0:
            cpu_dir = detect_new_dir(out_base, list(before_dirs))
            if not cpu_dir:
                print("Could not detect Python CPU output folder")
            else:
                cpu_dir = rename_output_folder(cpu_dir, 'python_cpu')
                before_dirs.add(os.path.basename(cpu_dir))
        else:
            print("Python CPU run failed")

    # PY GPU
    gpu_dir = None
    if args.compare in ("gpu","both"):
        print("=== Running Python GPU ===")
        rc = run_and_capture([sys.executable, 'make_galaxy_MCG.py', '--backend','python','--device','cuda','--dtype', args.dtype])
        if rc == 0:
            gpu_dir = detect_new_dir(out_base, list(before_dirs))
            if not gpu_dir:
                print("Could not detect Python GPU output folder")
            else:
                gpu_dir = rename_output_folder(gpu_dir, 'python_gpu')
                before_dirs.add(os.path.basename(gpu_dir))
        else:
            print("Python GPU run failed (ensure CUDA/CuPy are installed and configured)")

    # Compare
    print("=== Comparisons ===")
    f_noiseless, f_final = find_cube_files(f_dir)

    if cpu_dir:
        c_noiseless, c_final = find_cube_files(cpu_dir)
        print("-- Fortran vs Python CPU: noiseless --")
        data_report(c_noiseless, f_noiseless, "noiseless")
        print("-- Fortran vs Python CPU: final --")
        data_report(c_final, f_final, "final")

    if gpu_dir:
        g_noiseless, g_final = find_cube_files(gpu_dir)
        print("-- Fortran vs Python GPU: noiseless --")
        data_report(g_noiseless, f_noiseless, "noiseless")
        print("-- Fortran vs Python GPU: final --")
        data_report(g_final, f_final, "final")

    if not args.keep_temp:
        # keep the Fortran dir; remove python dirs by default
        for d in [cpu_dir, gpu_dir]:
            if d and os.path.isdir(d):
                shutil.rmtree(d, ignore_errors=True)
        print("Removed temporary Python output folders (use --keep-temp to retain)")


if __name__ == '__main__':
    main()


