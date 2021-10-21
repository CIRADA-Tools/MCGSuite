#!/usr/bin/env python3
import numpy as np


class BasicIO:
    #   Set the full dictionary of BasicIO attributes
    __slots__ = ['Beta_InputsFile', 'OutputFolder','GalaxyName','ProfileFile','OutputCubeNames','MCG_Main_InputFile','MCG_DataCube_InputFile','MCG_TiltedRing_InputFile','plot_verbose','file_verbose','output_volume_switch','MapPlotName','ProfilePlotName']
    #   Basic IO initializing function
    def __init__(self):
        #   By Default no Beta file is specified
        self.Beta_InputsFile=None
        #   Set the default output folder name
        self.OutputFolder='DefaultFolder'
        #   Set the default name for the galaxy -- will be determined by various parameters
        self.GalaxyName=None
  
        #   Set the default name for the MCG input file
        self.MCG_Main_InputFile='MCGMain.in'
        #   Set the default name for the MCG datacube file
        self.MCG_DataCube_InputFile='MCG_DataCube.in'
        #   Set the default name for the MCG tilted ring file
        self.MCG_TiltedRing_InputFile='MCG_TiltedRing.in'
        #   Set the default level for the number of output plots
        self.plot_verbose=False
        #   Set the default level for the number of output files
        self.file_verbose=False

        #   Set the default output MCG file volume switch.  This will be determined by the plot and file verbose flags
        self.output_volume_switch=2


class Suite:
    #   Set the full dictionary of Suite attributes
    __slots__=['mass_array','beams_array','inclination_array','pa_array','veldisp_array','num_realizations','num_array','CatalogueArray','n_galaxies','Templates','SuiteIO','nProcessors','DBTable']
    def __init__(self):
        self.mass_array=None

        self.beams_array=None

        self.inclination_array=None

        self.pa_array=None
        
        self.veldisp_array=None

        self.num_realizations=0
        
        self.num_array=None

        self.n_galaxies=None

        self.CatalogueArray=None

        #      Create a BasicIO instance for the suite
        self.SuiteIO=BasicIO()
        #   Create a profiles, tilted ring, and data cube template for the suite that will be used for initialzing each instance of a galaxy.
        ProfilesTemplate=Profiles()
        TiltedRingTemplate=TiltedRing()
        DataCubeTemplate=DataCube()
        self.Templates=[ProfilesTemplate,DataCubeTemplate,TiltedRingTemplate]
            #By default use 1 processor
        self.nProcessors=1

        self.DBTable=None


class Galaxy:
    #   Set the full dictionary of Galaxy attributes
    __slots__=['RC','SB','logMHI','nBeams','inclination','pa','RHI','distance','VHI','SizeLimits','MassLimits','veldisp','dispersionLimits','version_number','ID']
    #   Galaxy Initializing function
    def __init__(self):
        #   Initialize the Rotation Curve sub class
        self.RC=self.RotationCurve()
        #   Initialize the Surface Brightness sub class
        self.SB=self.SurfaceBrightness()
        #   Set the default HI mass to 10^10
        self.logMHI=10.
        #   Set the default number of beams across the galaxy
        self.nBeams=5
        #   Set the default inclination
        self.inclination=0.
        #   Set the default position angle
        self.pa=0.
        #   Set the default RHI size (kpc) -- will be recalculated
        self.RHI=5.
        #   Set the default distance (Mpc) -- will be recalculated
        self.distance=10.
        #   Set the default velocity disperions (km/s) -- will be calculated.
        self.veldisp=8.
        #   Set the default velocity at RHI (km/s) -- will be set by the user
        self.VHI=100.
        #   Set the minimum and maximum number of beams
        self.SizeLimits=[2,50]
          #   Set the minimum and maximum logarithmic HI mass
        self.MassLimits=[7.5,10.5]
        #   Set the minimum and maximum velocity dispersion in km/s
        self.dispersionLimits=[-1.,20.]
        #   Set the version number of the galaxy object to None (only used for making a full suite)
        self.version_number=None
        #   Set the ID number of the galaxy object to None (only used for making a full suite)
        self.ID=None

    class RotationCurve:
        #   Rotation Curve Initialization
        __slots__=['rPE','vPE','aPE']
        
        def __init__(self):
            #   Set the default scale length (kpc).
            self.rPE=1.
                #   Default scale velocity (km/s)
            self.vPE=100.
            #   Default shape parameter (unitless)
            self.aPE=1.


    class SurfaceBrightness:
        #   Surface Brightness Initialzation
        __slots__=['sbmax','Grmax','Gsig','Er','Vsig']
        
        def __init__(self):
            #   Set the default peak scale brightness to 1 (M_sol/pc^2).
            self.sbmax=1.
            #   Set the default radius of the gaussian peak in kpc.
            self.Grmax=1.
            #   Set the default standard deviation of the gaussian peak in kpc.
            self.Gsig=1.
            #   Set the default inner exponential scale length in kpc.
            self.Er=1.
            #   Set the default threshold velocity in km/s.
            self.Vsig=1.


class Profiles:
    #   Profiles Initialization
    __slots__=['nBinsPerRHI','limR_RHI','step','R','Vrot','SB']
    def __init__(self):
        #   Set the number of profile points per HI radius
        self.nBinsPerRHI=100.
        #   Set the maximum profile radius in units of HI radius
        self.limR_RHI=1.2


class DataCube:
    #   Data cube object initialization
    __slots__=['cube_shape','pixel_size','channel_size','cube_dimensions','beam_fwhm','beam_flattening','beam_angle','noise','cube_margin','cubesize_minimum','velocity_margin','velocity_smoothing_switch','nSigma','velocity_smooth_sigma','ChannelSizeLimits','BeamPixelRatioLimits','beam_dimensions','reference_locations','reference_values']
    def __init__(self):
        #   Default cube size
        self.cube_shape=None
        #   Default pixel size in arcsec
        self.pixel_size=4.
        #   Default channel size in km/s
        self.channel_size=4.
        #   Default beam FWHM radius (arcsec)
        self.beam_fwhm=30.
        #   Default beam flattening (0<flattening<1)
        self.beam_flattening=1.
        #   Default beam position angle (degrees)
        self.beam_angle=0.
        #   The default noise level in mJy/Beam
        self.noise=1.6
        #   The default size of the margins in beam size
        self.cube_margin=10
        #   The minimum size of the cube in terms of beams
        self.cubesize_minimum=10.
        #   The default size of the velocity margin in km/s
        self.velocity_margin=50.
        #   Default switch for velocity smoothing (0=none, 1=Gaussian)
        self.velocity_smoothing_switch=0


        #   Default value for the calculated extent of the beam in sigma units -- can be reset using BetaConfig
        self.nSigma=5.
        #   Default value for the sigma when velocity smoothing is turned on -- can be reset using BetaConfig
        self.velocity_smooth_sigma=None
    
        #   Limits on the possible channel size and beam/pixel ratio
        self.ChannelSizeLimits=[1.,20.]
        self.BeamPixelRatioLimits=[2.,10.]


class TiltedRing:
    
    __slots__=['cmode','CloudSurfDens','nRingsPerBeam','RHI_Limit','ra_center','dec_center','inclination','position_angle','vsys','v_radial','v_vertical','v_dispersion','dvdz','ranseed','nRings','Rwidth','R_array','v_tangential_array','sigma_array','z_scale','z_gradiantstart']

    def __init__(self):
        #   Default cloud density mode (allowable values are 0, 1, 2)
        self.cmode=0
        #   Default cloud # surface density (#/pixel)
        self.CloudSurfDens=100
        
        #   Default number of rings per beam
        self.nRingsPerBeam=5.
        
        #   Set the limit to calculate the tilted ring model to in terms of RHI
        self.RHI_Limit=1.5


        #   Default RA center (degrees)
        self.ra_center=0.
        #   Default declination center (degrees)
        self.dec_center=0.
        #   Default inclination (degrees)
        self.inclination=45.
        #   Default position angle (degrees)
        self.position_angle=0.
        #   Default systemic velocity (km/s)
        self.vsys=0.
        #   Default radial velocity (km/s)
        self.v_radial=0.
        #   Default vertical velocity (km/s)
        self.v_vertical=0.
        #   Default velocity dispersion (km/s)
        self.v_dispersion=0.
        #   Default tangential velocity gradient as a function of vertical height (km/s arcsec^-1)
        self.dvdz=0.

        #   Default random seed for MCG code -- positive value will use the time.
        self.ranseed=1
