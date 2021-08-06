#!/usr/bin/env python3
from src.make_galaxy_code import *

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
#   Make the cubes (currently using MCG)
MC.MakeCubes(GalaxyIO,DataCube,TiltedRing)
#   Make all Plots
DP.MakeAllPlots(GalaxyIO,Galaxy,DataCube,Profiles,TiltedRing)
#   Clean up outputs
OutClean.CleanOutput(GalaxyIO)
