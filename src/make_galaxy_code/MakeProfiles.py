#!/usr/bin/env python3
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from . import ObjectDefinitions as OD
from . import StandardRelations as SR

def MakeProfiles(Profiles,Galaxy,GalaxyIO):
     print("Generating Profiles")
     #
     #################################################
     #
     # Set the profile sampling and generate profiles:
     #
     #################################################
     #

     #   Radial profile in kpc:
     Profiles.R = Profiles.R*Galaxy.RHI
     #   Calculate Vrot for that profile:
     Profiles.Vrot=SR.CalcV(Galaxy.RC,Profiles.R)
     #   Calculate the surface brightness for that profile:
     Profiles.SB=SR.CalcSB(Galaxy.SB,Galaxy.VHI,Profiles.R,Galaxy.UDG_switch,Galaxy.logMHI,Galaxy.RHI)
     #
     #
     ####################################           
     #
     # Generate the output profile file:
     #
     ####################################
     #
     file=open(GalaxyIO.ProfileFile,"w")
     #
     #  Output user inputs:
     #
     file.write("# \n")
     file.write("# Galaxy Model: MCG v1.0, copyright 2020 \n")
     file.write("# \n")
     file.write("# User Inputs: \n")
     file.write("#------------- \n")
     file.write("# log10(HI mass/Mo)= {:.2f} \n".format(Galaxy.logMHI))
     file.write("# Inclination (deg)= {:.2f} \n".format(Galaxy.inclination))
     file.write("# Position angle (deg)= {:.2f} \n".format(Galaxy.pa))
     file.write("# Velocity dispersion (km/s)= {:.2f} \n".format(Galaxy.veldisp))
     file.write("# Beams across major axis= {:.2f} \n".format(Galaxy.nBeams))
     file.write("# \n")
     #
     #  Output calculated model parameters:
     #
     file.write("# Calculated Model Parameters: \n")
     file.write("#----------------------- \n")
     file.write("# HI radius at NHI = 1 Mo/pc^2 (kpc)= {:.2f} \n".format(Galaxy.RHI))
     file.write("# Rotation velocity at NHI = 1 Mo/pc^2 (km/s)= {:.2f} \n".format(Galaxy.VHI))
     file.write("# Distance (Mpc)= {:.2f} \n ".format(Galaxy.distance))
     file.write("# Rotation curve parameters, where V = vPE * (1-exp(-r/rPE) * (1 + aPE*r/rPE) \n")
     file.write("#          vPE (km/s)= {:.2f} \n".format(Galaxy.RC.vPE))
     file.write("#          rPE (kpc)= {:.2f} \n".format(Galaxy.RC.rPE))
     file.write("#          aPE (unitless)= {:.2f} \n".format(Galaxy.RC.aPE))
     file.write("# Surface brightness profile parameters, where SB = XXX \n")
     file.write("#          sbmax (Mo/pc^2)= {:.2f} \n".format(Galaxy.SB.sbmax))
     file.write("#          Grmax (kpc)= {:.2f} \n".format(Galaxy.SB.Grmax))
     file.write("#          Gsig (kpc)= {:.2f} \n".format(Galaxy.SB.Gsig))
     file.write("#          Er (kpc)= {:.2f} \n".format(Galaxy.SB.Er))
     file.write("#          Vsig (km/s)= {:.2f} \n".format(Galaxy.SB.Vsig))
     file.write("# \n")
     file.write("# \n")
     #
     # Output analytical model profiles:
     #
     file.write("# Analytical Rotation Curve, Surface Brightness profiles from r = [0,RHI]: \n")
     file.write("#--------------------------------------------------------- \n")
     file.write("#   r        Vrot(r)    SB(r) \n")
     file.write("# (kpc)      (km/s)   (Mo/pc^2) \n")
     file.write("#--------------------------------------------------------- \n")
     for i in range(len(Profiles.R)):
         file.write("  {:.2f}     {:.2f}     {:.2f} \n".format(Profiles.R[i],Profiles.Vrot[i],Profiles.SB[i]))
     #
     file.close()

     return Profiles


