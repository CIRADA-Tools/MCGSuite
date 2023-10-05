#!/usr/bin/env python3
import numpy as np
from scipy.optimize import curve_fit
from scipy import integrate
import copy as copy

from . import ObjectDefinitions as OD
from . import StandardRelations as SR


def MakeGalaxy(gal,dc):
#Define the galaxy model.
#--> INPUT: galaxy object and datacube object as defined in OD
#--> OUTPUT: galaxy object, with RHI, distance,
#            VHI,veldisp,RotationCurve, SurfaceBrightness
#            attributes populated
#--------------------
#
    print("Making the Galaxy Model")
    #
    #Set the HI radius RHI in kpc from MHI:
    gal.RHI = setRHI( gal.logMHI )
    #
    #Set the distance in Mpc from RHI and nBeams: 
    gal.distance = ( 2.* gal.RHI * 206265. )/( gal.nBeams * dc.beam_fwhm * 1000. )
    #
    #Set the HI RC velocity at RHI, in kms, from MHI:
    gal.VHI = setVHI(gal.logMHI,gal.UDG_switch,gal.v_HI)
    #
    #Set the Polyex rotation curve parameters:
        Polyex_V, Polyex_r, Polyex_a = setRC(gal.RHI,gal.VHI,gal.UDG_switch)
    gal.RC.vPE = copy.deepcopy(Polyex_V)
    gal.RC.rPE = copy.deepcopy(Polyex_r)
    gal.RC.aPE = copy.deepcopy(Polyex_a)
    #
    #Set the constant velocity dispersion:
    #gal.veldisp = 10.0
    #
    #Set the surface brightness profile:
    SB_Grmax, SB_Er, SB_Gsig, SB_Vsig, SB_SBmax = setSB(gal.logMHI,gal.RHI,gal.VHI)
    gal.SB.Grmax = SB_Grmax
    gal.SB.Er = SB_Er
    gal.SB.Gsig = SB_Gsig
    gal.SB.Vsig = SB_Vsig
    gal.SB.sbmax = SB_SBmax


    ########Sanity Checks to go here#########
    
    return gal



def setRHI(logMHI):    
#Use the HI mass - diameter relation from
#Wang+2016 to estimate RHI from logMHI.
#-->INPUT: log(MHI/Mo)
#-->OUTPUT: RHI in kpc
#--------------------
    slope = 0.506
    intercept = -3.293
    logDHI = slope*logMHI + intercept
    RHI = (10**logDHI)/2.
    return RHI


def setVHI(logMHI,UDG_switch,v_HI):
#Use relationship derived from abundance-matching
#the a.40 HIMF and VF to find VHI from logMHI.
#--> INPUT: log(MHI/Mo)
#--> OUTPUT: VHI in km/s
#--------------------
     z = np.array([0.0345, -0.955, 9.134, -27.99])
     z = np.array([0.0012992, -0.0620705, 1.210101, -11.98959, 60.46594, -122.66267])
     y = np.poly1d(z)
     V = 10**y(logMHI)
     # Set max at 300km/s:
     VHI = np.amin([V,300.0])
     #If UDG mode is activated, the user inputs the v_HI they want
     if UDG_switch==True:
          VHI=v_HI
     return VHI

 
def setRC(RHI,VHI,UDG_switch):
#Define the RC through the parameters of the
#Polyex function:
# Vrot = VPE * (1-exp(-r/rPE))*(1+aPE*r/rPE)
#--> INPUT: RHI in kpc, VHI in km/s
#--> OUTPUT: VPE in km/s, rPE in kpc, aPE (unitless)
#--------------------
    #
    # Functional forms:
    def linear(x,m,b):
         return m*x + b
    # 
    def exp_plus_lin(x,a,b,c,d):
         return a*np.exp(-x*b-c) +d*x
    # 
    # Set Vopt (= Vrot at Ropt as defined by Catinella+06) from VHI:
    Vopt = 0.95*VHI
    #
    # Fit Catinella+06 table 2 values to get relations between polyex params and magnitude:
    CVPE  = np.array([275.,255.,225.,200.,170.,148.,141.,122.,103.,85.])  #kms
    CdVPE = np.array([6.,2.,1.,1.,1.,2.,2.,2.,2.,5.])  #kms
    CrPE  = np.array([0.126,0.132,0.149,0.164,0.178,0.201,0.244,0.261,0.260,0.301])  #NOTE: /Ropt. CONVERSION TO KPC AT END.
    CdrPE = np.array([0.007,0.003,0.003,0.002,0.003,0.004,0.005,0.008,0.008,0.002])  #NOTE: /Ropt. CONVERSION TO KPC AT END.
    CaPE  = np.array([0.008,0.002,0.003,0.002,0.011,0.022,0.010,0.020,0.029,0.019])
    CdaPE = np.array([0.003,0.001,0.001,0.001,0.001,0.002,0.003,0.005,0.005,0.015])
    CMag  = np.array([-23.76,-23.37,-22.98,-22.60,-22.19,-21.80,-21.41,-21.02,-20.48,-19.38])  #Iband
    #
    CVPE_best, foo = curve_fit(exp_plus_lin,CMag,CVPE,sigma=CdVPE)
    CrPE_best, foo = curve_fit(linear, CMag, CrPE,sigma=CdrPE)
    CaPE_best, foo = curve_fit(linear, CMag, CaPE,sigma=CdaPE)
    #
    Magvec = np.arange(-24.,0.,0.001)
    CVPEvec = exp_plus_lin(Magvec,*CVPE_best)
    CrPEvec = linear(Magvec,*CrPE_best)
    CaPEvec = linear(Magvec,*CaPE_best)
    #
    Vrot_at_Ropt = CVPEvec * (1. - np.exp(-1/CrPEvec))*(1 + CaPEvec/CrPEvec)
    ind = np.argmin(abs(Vrot_at_Ropt - Vopt))
    VPE_temp = CVPEvec[ind]
    rPE = CrPEvec[ind]
    aPE_temp = CaPEvec[ind]
    
    #If UDG switch is on, just set the parameters to the lowest mass option in$
    if UDG_switch==True:
        VPE_temp=275 #units?
        rPE_temp=0.126  #ratio of optical radii
        aPE_temp=0.008
    
    # Now find the outer slope corresponding to RHI from Dutton+19:
    logDHI = np.log10(2*RHI)
    Vslope = -0.385*logDHI + 0.666
    # Find the value of aPE that produces this slope, assuming that RHI = 1.7Ropt as in Broeils+97:
    aPErange = np.arange(0.,0.8,0.001)
    logRin=np.log10((1-np.exp(-0.85/rPE))*(1+0.85*aPErange/rPE)) 
    logRout=np.log10((1-np.exp(-1.7/rPE))*(1+1.7*aPErange/rPE)) 
    log_slope_range = (logRout - logRin)/np.log10(2)
    ind = np.argmin(abs(log_slope_range - Vslope))
    aPE = aPErange[ind]
    # Finally, renormalize to find VPE where VHI is at RHI, with RHI = 1.7Ropt as in Broeils+97:
    VPE = VHI/((1-np.exp(-1.7/rPE))*(1 + 1.7*aPE/rPE))
    #
    return VPE, rPE*(RHI/1.7), aPE


def setSB(logMHI,RHI,VHI):
#Define the surface density profile
#using a Gaussian and exponential:
#SB/SBmax = exp(-(r - Rsigm)^2/2sigsig^2) +
#           (1 - sqrt(VHI/Vsig))exp(-r/Rs)
#--> INPUT: log(MHI/Mo), RHI in kpc, VHI in km/s
#--> OUTPUT: Grmax/kpc, Gsig/kpc, Vsig/kms, Er/kpc, SBmax/Mopc-2
#--------------------
#   
    # Vsig is fixed:
    Vsig = 120.0
    # Rsigm is fixed:
    Grmax = 0.4*RHI
    Er = 0.2*RHI
    Gsig = 0.33*RHI  #obtained by trial and error
    # Now find SBmax such that SB at RHI = 1 Mo/pc^2:
    SBmax = 1/(np.exp(-1.0*(RHI - Grmax)**2/(2*Gsig**2)) + (1.0 - np.sqrt(VHI/Vsig)) * np.exp(-1.0*RHI/Er))
    #
    # All of this is superfluous now that I have Gsig:
    #
    # Define the profile Gsig:
    #radvals = np.arange(0.,RHI,RHI/1000.)
    #SBr = np.zeros_like(radvals)
    #for i in range(len(radvals)):
    #    SBval = SBmax*(np.exp(-1.0*(radvals[i] - Grmax)**2/(2*Gsig**2)) + (1.0 - np.sqrt(VHI/Vsig)) * np.exp(-1.0*radvals[i]/Er))
    #    #Make sure that surface density isn't negative in the centre:
    #    SBr[i] = np.amax([SBval,0.])
    #    
    #Use this line to integrate the profile to check the HI mass. Note that the integral out to RHI should be 0.85*MHI:
    #Mass = (integrate.trapz(SBr*2*np.pi*radvals,radvals)*1000.**2.)
    #print(np.log10(Mass)-logMHI)
    return Grmax, Er, Gsig, Vsig, SBmax






