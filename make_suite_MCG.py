#!/usr/bin/env python3
import argparse
from make_galaxy_code import *

#from joblib import Parallel, delayed
import multiprocessing as mp
import numpy as np
from multiprocessing import freeze_support


def Main():
    print("Making a suite of galaxies using MCG")

    # CLI overrides for backend/device/runtime knobs
    parser = argparse.ArgumentParser(description="Generate a suite of mock galaxy cubes")
    parser.add_argument("--backend", "-b", choices=["fortran", "python"], help="Select backend (default: fortran unless MCG_BACKEND env set)")
    parser.add_argument("--device", "-d", choices=["cpu", "cuda"], help="Device for Python backend")
    parser.add_argument("--dtype", choices=["float32", "float64"], help="Numeric dtype for Python backend")
    parser.add_argument("--chunk-size", type=int, help="Chunk size along spectral axis for Python backend")
    parser.add_argument("--threads", type=int, help="Thread count for CPU FFT (Python backend)")
    parser.add_argument("--seed", type=int, help="Random seed forwarded to Python backend")
    parser.add_argument("--n-proc", type=int, help="Override number of worker processes for the suite")
    args = parser.parse_args()

    #   Get the inputs and set up the objects
    Suite=IN.Get_MakeSuite_Inputs()    #From config files

    # Apply CLI overrides to Suite IO and pool size
    if args.n_proc is not None:
        Suite.nProcessors = int(args.n_proc)
    if args.backend:
        Suite.SuiteIO.backend = args.backend
    if args.device:
        Suite.SuiteIO.device = args.device
    if args.dtype:
        Suite.SuiteIO.dtype = args.dtype
    if args.chunk_size is not None:
        Suite.SuiteIO.chunk_size = args.chunk_size
    if args.threads is not None:
        Suite.SuiteIO.threads = args.threads
    if args.seed is not None:
        Suite.SuiteIO.seed = args.seed

    #   Do the safety checks to make sure all parameters are within acceptable limits
    SC.SuiteChecks(Suite,Suite.Templates[1])
    #   Create a combined array of all combinations of the suite parameters
    OC.SuiteConfig(Suite)

    #   Set up parallel processing
    print("n Processors", Suite.nProcessors)

    #   Run the main loop for the processing
    with mp.Pool(processes=Suite.nProcessors) as pool:
        Suite.DBTable=pool.starmap(SO.SuiteMainLoopFn, [(i,Suite) for i in range(Suite.n_galaxies)])

    SO.CatalogueOutput(Suite)


if __name__=="__main__":
    freeze_support()
    Main()


