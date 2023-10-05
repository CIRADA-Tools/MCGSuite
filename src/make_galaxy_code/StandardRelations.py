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


def SD_UDG_exp(R,A,b,RHI):
    #Functional form in Eqn. 1 in arxiv: 2204.05981
    #Both A and b are free parameters that are scaled based off MHI and RHI
    n=4
    m=4
    SD=A*(1+(b*R/RHI)**m)**(-n/m)
    return SD

def scale_exp(R_integrate,RHI,VHI,log_MHI):
    #This function finds the optimal A and b parameter values that satisfies the following relations
    # i) Surface density at R=RHI is equal to 1 Mo/pc^2
    # ii) Integral of surface density over all space is equal to MHI
    
    N=50  #number of steps used in algorithm finding optimal A and b values for SD functional form
    A=4.1
    b=np.linspace(1,2,N)
    #Looping over potential b-values
    for i in range(N):
        #Checking that surface density at RHI is equal to 1 mo/pc^2
        if SD_UDG_exp(RHI,A,b[i],RHI)>0.99 and SD_UDG_exp(RHI,A,b[i],RHI)<1.01:
            #Copying b-value that satisfies above criteria
            b1=b[i]
            break
        
    thresh=0.05  #threshold that determines how much b value can vary from best-fit value
    low_b1=b1-thresh*b1
    high_b1=b1+thresh*b1
    #Varying A-parameter to ensure integral of SD is equal to MHI
    A_array=np.linspace(3.9,5.5,N)
    new_b_array=np.linspace(low_b1,high_b1,N)
    for j in range(N):
        for k in range(N):
            #integrating surface density to check MHI condition
            integrand=np.trapz(SD_UDG_exp(R_integrate,A_array[j],new_b_array[k],RHI)*R_integrate*1000,R_integrate*1000)
            ratio=abs(integrand*2*np.pi)/(10**log_MHI)
            if ratio>0.95 and ratio<1.05:
                return A_array[j],new_b_array[k]
 

def MP22_SD(r):
    #This is the functional form published in "No need for DM in AGC 114905"
    #https://arxiv.org/abs/2112.00017
    SD0=5.2
    R1=1.1
    R2=16.5
    alpha=18
    SD=SD0*(np.exp(-r/R1))*(1+r/R2)**alpha
    print('Using the surface density functional form in MP22 publication')
    return SD/1.33333
    
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
    #If UDG switch is on, recalculate surface density based on functional form in Eqn. 1 in arxiv: 2204.05981
    if UDG_switch==True:
        R_integrate=np.linspace(0.01,RHI*10.5,len(R))
        A,b=scale_exp(R_integrate,RHI,VHI,log_MHI)
        if A==None:
          print('Invalid Input: MHI is outside of UDG Module capabilities')
        SB=SD_UDG_exp(R,A,b,RHI)
    #If the user wants to use the functional form in "No need for DM in AGC 114905"
    if UDG_switch=='114905':
        SB=MP22_SD(R)
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
