#!/usr/bin/env python3
import numpy as np
import os

def CleanOutput(GalaxyIO):
# Cleans up the outputs from make_galaxy and MCG by either removing files or placing
# them into a user specified folder.
#-->INPUT: GalaxyIO = a GalaxyIO object
#--------------------
    #   Make the target folder to hold all the outputs
    if os.path.isdir(GalaxyIO.OutputFolder) == False:
        os.system("mkdir " +GalaxyIO.OutputFolder)

#

    os.system("mv "+ GalaxyIO.ProfileFile+ " " + GalaxyIO.GalaxyName)

    if GalaxyIO.plot_verbose:
        os.system("mv "+ GalaxyIO.ProfilePlotName+ " " + GalaxyIO.GalaxyName)
        os.system("mv "+ GalaxyIO.MapPlotName+ " " + GalaxyIO.GalaxyName)

    if GalaxyIO.file_verbose:
        os.system("mv "+ GalaxyIO.MCG_Main_InputFile+ " " + GalaxyIO.GalaxyName)
        os.system("mv "+ GalaxyIO.MCG_DataCube_InputFile+ " " + GalaxyIO.GalaxyName)
        os.system("mv "+ GalaxyIO.MCG_TiltedRing_InputFile+ " " + GalaxyIO.GalaxyName)
    else:
        os.system("rm "+ GalaxyIO.MCG_Main_InputFile)
        os.system("rm "+ GalaxyIO.MCG_DataCube_InputFile)
        os.system("rm "+ GalaxyIO.MCG_TiltedRing_InputFile)

    if os.path.isdir(GalaxyIO.OutputFolder+"/"+GalaxyIO.GalaxyName) == False:
        os.system("mv "+ GalaxyIO.GalaxyName+ " " + GalaxyIO.OutputFolder)
    else:
        os.system("rm -r "+ GalaxyIO.OutputFolder+"/"+GalaxyIO.GalaxyName)
        os.system("mv "+ GalaxyIO.GalaxyName + " " + GalaxyIO.OutputFolder)
