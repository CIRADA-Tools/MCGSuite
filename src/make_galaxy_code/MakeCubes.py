#!/usr/bin/env python3
import numpy as np
import os
from . import ObjectDefinitions as OD


def MakeCubes(GalaxyIO,DataCube,TiltedRing):
    print("Making Cubes")

    #   Write an MCG main input file
    WriteMCG_MainInput(GalaxyIO,DataCube)

    #   Write an MCG data cube input file
    WriteMCG_DataCubeInput(GalaxyIO,DataCube)

    #   Write an MCG tilted ring input file
    WriteMCG_TiltedRingInput(GalaxyIO,TiltedRing)

    #   Run MCG using the input files
    os.system("./Programs/MockCubeGenerator "+GalaxyIO.MCG_Main_InputFile)

def WriteMCG_MainInput(GalaxyIO,DataCube):
    file=open(GalaxyIO.MCG_Main_InputFile,"w")
    file.write("#    Base name for the output folder and all file names \n")
    file.write(GalaxyIO.GalaxyName+"\n")
    
    file.write("#    Minimal, Moderate, of Full Outputs (0,1, and 2 respectively)\n")
    file.write(str(GalaxyIO.output_volume_switch)+"\n")
    
    file.write("#    Name of the file containing the cube and beam definitions\n")
    file.write(GalaxyIO.MCG_DataCube_InputFile+"\n")

    file.write("#    Name of the file containing the underlying ring model\n")
    file.write(GalaxyIO.MCG_TiltedRing_InputFile+"\n")

    file.write("#    Noise Unit Switch (0=mJy/beam)\n")
    file.write("0 \n")
    
    file.write("#    Noise Value\n")
    file.write(str(DataCube.noise)+"\n")
    
    file.write("#    Random Seed\n")
    file.write("1 \n")

    file.close()

def WriteMCG_DataCubeInput(GalaxyIO,DataCube):
    file=open(GalaxyIO.MCG_DataCube_InputFile,"w")

    file.write("#    number of pixels and channel\n")
    file.write(str(DataCube.cube_shape[0]) +"\t" +str(DataCube.cube_shape[1]) +"\t"+str(DataCube.cube_shape[2]) +"\t"+" \n")
    
    file.write("#    Pixel Size Units (0=degrees, 1=arcsec) - Reference Pixel Units (0=degrees) - Channel Size Units (0=m/s, 1=km/s) - Reference Channel Units (0=km/s)  — Beam axis Units (0=arcsec) — Beam Rotation Angle Units (0=degrees)\n")
    file.write("1\t 0 \t 1 \t 1 \t 1 \t 0\n")
    
    file.write("#    Pixel dimensions and channel size\n")
    file.write(str(DataCube.cube_dimensions[0]) +"\t" +str(DataCube.cube_dimensions[1]) +"\t"+str(DataCube.cube_dimensions[2]) +"\t"+" \n")
    
    file.write("#    Reference Pixel(channel) in each dimension (CRPIX values)\n")
    file.write(str(DataCube.reference_locations[0]) +"\t" +str(DataCube.reference_locations[1]) +"\t"+str(DataCube.reference_locations[2]) +"\t"+" \n")
    
    file.write("#    Reference value in each dimension (CRVAL values)\n")
    file.write(str(DataCube.reference_values[0]) +"\t" +str(DataCube.reference_values[1]) +"\t"+str(DataCube.reference_values[2]) +"\t"+" \n")
    print("ref value check", DataCube.reference_values)
    
    file.write("#    Beam dimensions (major, minor, position angle)\n")
    file.write(str(DataCube.beam_dimensions[0]) +"\t" +str(DataCube.beam_dimensions[1]) +"\t"+str(DataCube.beam_dimensions[2]) +"\t"+" \n")
    
    file.write("#    number of sigma lengths to reach with the beam\n")
    file.write(str(DataCube.nSigma)+" \n")
    
    file.write("#    Type of velocity smoothing to use (0=none, 1=Gaussian)\n")
    file.write(str(DataCube.velocity_smoothing_switch)+" \n")
    
    file.write("#    The velocity smoothing sigma if using Gaussian smoothing (km/s)\n")
    file.write(str(DataCube.velocity_smooth_sigma)+" \n")

    file.close()

def WriteMCG_TiltedRingInput(GalaxyIO,TiltedRing):
    file=open(GalaxyIO.MCG_TiltedRing_InputFile,"w")
    
    #file.write("#        The format of the input file — 1=row, 2= column\n")
    #file.write("2 \n")
    
    file.write("#   Number of Rings in the model\n")
    file.write(str(TiltedRing.nRings)+"\n")


    file.write("#    The cloud mode you are using\n")
    file.write(str(TiltedRing.cmode)+"\n")
    
    file.write("#    The base cloud surface density\n")
    file.write(str(TiltedRing.CloudSurfDens)+"\n")
    
    file.write("#    Central Position Switch (0=degrees, 1=arcsec), Inc/PA Unit switch (0=degrees, 1=arcsec), Velocity Units (0=m/s, 1=km/s), Brightness Units (0=Jy km/s arcsec^-2)\n")
    file.write("0 \t 0 \t 1 \t  0 \n")
    
    file.write("#    The parameters in each radial bin\n")
    file.write("#    Rmid    Rwidth    Xcent    Ycent    Inc    PA    VSys    VRot    VRad    Vvert    VDisp    dvdz    Sigma        z0    zGradStart\n")


    for i in range(TiltedRing.nRings):
        TiltedRing.z_gradiantstart[i]=5.*TiltedRing.z_scale[i]
        file.write(str(TiltedRing.R_array[i])+"\t"+str(TiltedRing.Rwidth) +"\t" + str(TiltedRing.ra_center) + "\t"+ str(TiltedRing.dec_center) +"\t" + str(TiltedRing.inclination) + "\t" +str(TiltedRing.position_angle) +"\t" + str(TiltedRing.vsys)+ "\t" + str(TiltedRing.v_tangential_array[i]) +"\t" + str(TiltedRing.v_radial) + "\t" + str(TiltedRing.v_vertical) +"\t" + str(TiltedRing.v_dispersion) + "\t"+ str(TiltedRing.dvdz) + "\t" + str(TiltedRing.sigma_array[i])+"\t" + str(TiltedRing.z_scale[i])+ "\t" + str(TiltedRing.z_gradiantstart[i]) +   "\n")

    file.close()

