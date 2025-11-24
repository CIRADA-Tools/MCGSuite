#!/usr/bin/env python3
import argparse
from make_galaxy_code import *

#   CLI: allow selecting backend/device and runtime knobs for Python backend
parser = argparse.ArgumentParser(description="Generate a mock galaxy cube")
parser.add_argument("--backend", "-b", choices=["fortran", "python"], help="Select backend (default: fortran unless MCG_BACKEND env set)")
parser.add_argument("--device", "-d", choices=["cpu", "cuda"], help="Device for Python backend")
parser.add_argument("--dtype", choices=["float32", "float64"], help="Numeric dtype for Python backend")
parser.add_argument("--chunk-size", type=int, help="Chunk size along spectral axis for Python backend")
parser.add_argument("--threads", type=int, help="Thread count for CPU FFT (Python backend)")
parser.add_argument("--seed", type=int, help="Random seed forwarded to Python backend")
args = parser.parse_args()

#   Get the inputs and set up the objects
Galaxy,DataCube,TiltedRing,Profiles,GalaxyIO=IN.GetMakeGalaxyInputs()    #From config files
#   Check the various limits to see if the code will work.
SC.FirstChecks(Galaxy,DataCube)
#   Calculate all the various galaxy parameters from the HI Mass
Galaxy=MG.MakeGalaxy(Galaxy,DataCube)
#   Configure the various output objects
DataCube,TiltedRing,Profiles,GalaxyIO=OC.ConfigObjects(Galaxy,DataCube,TiltedRing,Profiles,GalaxyIO)

#   Calculate the profiles based on the galaxy parameters
Profiles=MP.MakeProfiles(Profiles,Galaxy,GalaxyIO)
#   Calculate the full tilted ring model.
TiltedRing=MTR.MakeTiltedRing(Galaxy,DataCube,TiltedRing)
#   Apply CLI overrides before generating cubes
if args.backend:
    GalaxyIO.backend = args.backend
if args.device:
    GalaxyIO.device = args.device
if args.dtype:
    GalaxyIO.dtype = args.dtype
if args.chunk_size is not None:
    GalaxyIO.chunk_size = args.chunk_size
if args.threads is not None:
    GalaxyIO.threads = args.threads
if args.seed is not None:
    GalaxyIO.seed = args.seed

#   Make the cubes (currently using MCG)
MC.MakeCubes(GalaxyIO,DataCube,TiltedRing)
#   Make all Plots
DP.MakeAllPlots(GalaxyIO,Galaxy,DataCube,Profiles,TiltedRing)
#   Clean up outputs
OutClean.CleanOutput(GalaxyIO)
