#!/usr/bin/env python3
import numpy as np
from . import ObjectDefinitions as OD

def CalcV(RCObj,R):
# Returns the value of the polyex rotation velocity at a galactocentric
# Radius R given the parameters in RCObj of the Galaxy.RotationCurve subclass.
#-->INPUT: RCObj = rotation curve object of subclass Galaxy.RotationCurve
#          R = galactocentric radius R in kpc. 
#-->OUTPUT: V = rotation velocity at R, in km/s.
#--------------------
    V = RCObj.vPE * ( 1. - np.exp(-R/RCObj.rPE) ) * ( 1. + RCObj.aPE*R/RCObj.rPE )
    return V


def CalcSB(SBObj,VHI,R):
# Returns the value of the Gaussian+Exponential surface brightness profile
# at a galactocentric radius R given the parameters in SBObj of the
# Galaxy.SurfaceBrightness subclass.
#-->INPUT: SBObj = surface brightness object of subclass Galaxy.SurfaceBrightness
#          VHI = V(RHI), km/s.
#          R = galactocentric radius R in kpc. 
#-->OUTPUT: SB = surface brightness at R, in Mo/pc^2.
#--------------------
    Gauss = np.exp( -1.0 * (R - SBObj.Grmax)**2 / (2*SBObj.Gsig**2) )
    Exp = ( 1.0 - np.sqrt(VHI / SBObj.Vsig) ) * np.exp( -1.0*R/SBObj.Er )
    SBval = SBObj.sbmax * ( Gauss + Exp )
    Zarray = np.zeros_like(SBval)
    SB = np.maximum(SBval,Zarray)
    return SB



def CalcScaleHeight(VDisp,RCObj,R):
# TO BE CHECKED LATER.
    dR=0.1
    if dR < R[1]-R[0]:
        dR=R[1]-R[0]
    
    RNeg=R-dR
    if RNeg[0] < 0:
        RNeg[0]=0.
    RPos=R+dR
    #   Get the mass of the inner and outer sphere in units of 2.325e9 M_sol
    MNeg=MSphereFromVel(RCObj,RNeg)
    MPos=MSphereFromVel(RCObj,RPos)
    #   Get the mass of the shell around the radius
    MShell=MPos-MNeg
    #   Calculate the volume in kpc^3
    Vol=4./3.*np.pi*(RPos**3.-RNeg**3.)
    #   Get the spherical approximation for the denisty (units of 2.325e9 M_sol/kpc^3)
    rhoApproximation=MShell/Vol
    #   Get the Gaussian scale height in kpc using Equation 7 from Puche et al. 1992
    ZTest=(VDisp/100.)/(4.*np.pi*rhoApproximation)**0.5
    #   Convert to a sech^2 scale height.
    ScaleHeight=ZTest*0.9333
    #print(R)
    #print(ScaleHeight)
    #ScaleHeight=ScaleHeight/ScaleHeight*0.

    return ScaleHeight

def MSphereFromVel(RCObj,R):
    #   Get the velocity, but convert to 100 km/s units to set G=1
    V=CalcV(RCObj,R)/100.
    #   In kpc and 100 km/s units, the mass will be in units of 2.325e9 M_sol
    M=V*V*R
    return M

def ArcsecToKpc(size_arcsec,dist):
# Returns the conversion factor from arcsec to kpc.
#-->INPUT: size_arcsec = a size in arcsec
#          dist = distance in Mpc.
#-->OUTPUT: a_to_k, where size_kpc = size_arcsec * a_to_k.
#--------------------
    #a_to_k=size_arcsec/3600.*np.pi/180.*(dist*1000.)
    a_to_k=( size_arcsec/206265. ) * ( dist*1000. )
    return a_to_k

def KpcToArcSec(size_kpc,dist):
    # Returns the conversion factor from kpc to arcsec.
    #-->INPUT: size_kpc = a size in kpc
    #          dist = distance in Mpc.
    #-->OUTPUT: k_to_a, where size_arc = size_arcsec * a_to_k.
    #--------------------
    k_to_a=( size_kpc) / ( dist*1000. )*206265.
    return k_to_a

def SBUnitConversion(SB):
#   Converts a surface brightness from M_sol/pc^2 to Jy km/s arcsec^-2
#-->INPUT: SB = a surface brightness in M_sol/pc^2
#-->OUTPUT: SB = a surface brightness in Jy km/s arcsec^-2
    SB=SB*1.24756e+20/(6.0574E5*1.823E18*(2.*np.pi/np.log(256.)))
    return SB
