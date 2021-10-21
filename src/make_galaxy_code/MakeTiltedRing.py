#!/usr/bin/env python3
import numpy as np

from . import ObjectDefinitions as OD
from . import StandardRelations as SR

def MakeTiltedRing(Galaxy,DataCube,TiltedRing):
    print("Making the tilted ring object")
    #   First set the radii in arcsec.
    TiltedRing.R_array=CalcRarray(TiltedRing.R_array,TiltedRing.Rwidth)
    #   Next convert the radii from arcsec to kpc
    R_phys=SR.ArcsecToKpc(TiltedRing.R_array,Galaxy.distance)
    
    RTest=SR.KpcToArcSec(Galaxy.RHI,Galaxy.distance)
    #print("R Checks",RTest,Galaxy.RHI, Galaxy.distance)
    RTest2=SR.ArcsecToKpc(RTest,Galaxy.distance)
    #print(RTest2)
    #print(TiltedRing.R_array)
    
    #   Now calculate the velocity at each radial point
    TiltedRing.v_tangential_array=SR.CalcV(Galaxy.RC,R_phys)
    
    #   Next get the surface brightness at each radial point
    TiltedRing.sigma_array=SR.CalcSB(Galaxy.SB,Galaxy.VHI,R_phys)
    TiltedRing.sigma_array=SR.SBUnitConversion(TiltedRing.sigma_array)
    #TiltedRing.sigma_array=R_phys/R_phys
    #   Now get the scale height
    TiltedRing.z_scale=SR.CalcScaleHeight(Galaxy.veldisp,Galaxy.RC,R_phys)
    TiltedRing.z_scale=SR.KpcToArcSec(TiltedRing.z_scale,Galaxy.distance)
    #   Use the input inclination and position angle for the galaxy
    TiltedRing.inclination=Galaxy.inclination
    TiltedRing.position_angle=Galaxy.pa
    #   Use the input velocity dispersions for the galaxy
    TiltedRing.v_dispersion=Galaxy.veldisp
    
    return TiltedRing


def CalcRarray(R,Rwidth):
    for i in range(len(R)):
        R[i]=(float(i+0.5))*Rwidth
    return R

