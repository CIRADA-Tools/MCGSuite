#!/usr/bin/env python3
import numpy as np
import importlib

from Inputs.observatory_config_MCG import *
from . import ObjectDefinitions as OD


def Get_MakeSuite_Inputs():
    #   Import the suite config
    SuiteIn=importlib.import_module('Inputs.suite_config_MCG')
    #   Set up the Suite object
    Suite=OD.Suite()
    #   Get the suite inputs
    Suite=GetSuiteInputs(Suite,SuiteIn)
    

    #   And get their inputs from suite_config_MCG.py
    Suite.SuiteIO=GetBasicIO(Suite.SuiteIO,SuiteIn)
    print("SuiteIO", Suite.SuiteIO.OutputFolder, Suite.SuiteIO.Beta_InputsFile)
    #   ...and get the inputs from observatory_config_MCG.py
    Suite.Templates[1]=GetDataCubeInputs(Suite.Templates[1])
    

    #   Finally Check on beta configuation options
    BetaConfigInputs(Suite.SuiteIO,Suite.Templates[1],Suite.Templates[2],Suite.Templates[0])
    
    return Suite

def GetSuiteInputs(Suite,ModName):
    #   Define the error message to go with the exit signal
    ErrorMsg="= in "+ModName.__file__
    #   Try to get the array of masses
    Suite.mass_array=ModuleVarCheck_Exit(Suite.mass_array,'Masses',ModName,ErrorMsg)
    #   Try to get the array of beams
    Suite.beams_array=ModuleVarCheck_Exit(Suite.beams_array,'Beams',ModName,ErrorMsg)
    #   Try to get the array of inclinations
    Suite.inclination_array=ModuleVarCheck_Exit(Suite.inclination_array,'Inclinations',ModName,ErrorMsg)
    #   Try to get the array of position angles
    Suite.pa_array=ModuleVarCheck_Exit(Suite.pa_array,'PositionAngles',ModName,ErrorMsg)
    #   Try to get the array of velocity dispersions
    Suite.veldisp_array=ModuleVarCheck_Exit(Suite.veldisp_array,'veldisps',ModName,ErrorMsg)
    #   Try to get the number of realizations for each combination of parameters
    Suite.num_realizations=ModuleVarCheck_Exit(Suite.num_realizations,'NumRealizations',ModName,ErrorMsg)
    #   Try to get the number of processors to use
    Suite.nProcessors=ModuleVarCheck_Exit(Suite.nProcessors,'n_Processors',ModName,ErrorMsg)
    return Suite

def GetMakeGalaxyInputs():
    #   Import the galaxy config file
    GalaxyIn=importlib.import_module('Inputs.galaxyconfig_MCG')
    #   Set up the galaxy object
    Galaxy=OD.Galaxy()
    #   Now set the specific galaxy parameters based on the config file
    Galaxy=GetGalaxyInputs(Galaxy,GalaxyIn)

    #   Now set up the basic cube object
    DataCube=OD.DataCube()
    #   And read in the data cube values
    DataCube=GetDataCubeInputs(DataCube)
    
    Profiles=OD.Profiles()
       
    #   Next initialize the various I/O options
    GalaxyIO=OD.BasicIO()
    #   And read in the options
    GalaxyIO=GetBasicIO(GalaxyIO,GalaxyIn)
    
    #   Also Initialize the Titled Ring
    TiltedRing=OD.TiltedRing()
    
    #   Finally Check on beta configuation options
    BetaConfigInputs(GalaxyIO,DataCube,TiltedRing,Profiles)
    
    return Galaxy,DataCube,TiltedRing,Profiles,GalaxyIO


def GetGalaxyInputs(Galaxy,ModName):
    #   Set the generic error message
    ErrorMsg="= in "+ModName.__file__
    #   Try to set the HI Mass
    Galaxy.logMHI=ModuleVarCheck_Exit(Galaxy.logMHI,'Mass',ModName,ErrorMsg)
    #   Try to set the number of beams across the major axis
    Galaxy.nBeams=ModuleVarCheck_Exit(Galaxy.nBeams,'Beams',ModName,ErrorMsg)
    #   Try to set the Inclination
    Galaxy.inclination=ModuleVarCheck_Exit(Galaxy.inclination,'Inclination',ModName,ErrorMsg)
    #   Try to set the Position Angle
    Galaxy.pa=ModuleVarCheck_Exit(Galaxy.pa,'PositionAngle',ModName,ErrorMsg)
    #   Try to set the velocity dispersion
    Galaxy.veldisp=ModuleVarCheck_Exit(Galaxy.veldisp,'veldisp',ModName,ErrorMsg)
    #  Try to set the UDG_switch
    Galaxy.UDG_switch=ModuleVarCheck_Exit(Galaxy.UDG_switch,'UDG_switch',ModName,ErrorMsg)
    #  Try to set the v_HI value if UDG_switch=True
    Galaxy.v_HI=ModuleVarCheck_Exit(Galaxy.v_HI,'v_HI',ModName,ErrorMsg)
    
    return Galaxy


def GetDataCubeInputs(DataCube):
    #   Try to set the cell dimensions
    ErrorMsg="Pixel size is not set.  Set PixelSize= in observatory_config_MCG.py"
    DataCube.pixel_size=GlobalVarCheck(DataCube.pixel_size,'PixelSize',ErrorMsg,True)
    
    ErrorMsg="Channel size is not set.  Set ChannelSize= in observatory_config_MCG.py"
    DataCube.channel_size=GlobalVarCheck(DataCube.channel_size,'ChannelSize',ErrorMsg,True)
    
    #   Try to set the beam FWHM
    ErrorMsg="Pixel/Channel scales are not set.  Set cell_dimensions= in observatory_config_MCG.py"
    DataCube.beam_fwhm=GlobalVarCheck(DataCube.beam_fwhm,'beam_fwhm',ErrorMsg,True)
    
    #   Try to set the noise value
    ErrorMsg="Noise value not set.  Set noise_value= in observatory_config_MCG.py"
    DataCube.noise=GlobalVarCheck(DataCube.noise,'noise_value',ErrorMsg,True)

    return DataCube


def GetBasicIO(GalaxyIO,ModName):
    #   Set the generic error message
    ErrorMsg="= in "+ModName.__file__
    #   Try to set output folder
    GalaxyIO.OutputFolder=ModuleVarCheck_Exit(GalaxyIO.OutputFolder,'OutFolder',ModName,ErrorMsg)
    #   Figure out if there will be few or many plots
    GalaxyIO.plot_verbose=ModuleVarCheck_Exit(GalaxyIO.plot_verbose,'PlotVerbose',ModName,ErrorMsg)
    #   Figure out if there will be few or many files
    GalaxyIO.file_verbose=ModuleVarCheck_Exit(GalaxyIO.file_verbose,'FileVerbose',ModName,ErrorMsg)
    #   Check if there is a beta configuration file
    Msg=".py  Will be using beta configuration options."
    GalaxyIO.Beta_InputsFile=ModuleVarCheck(GalaxyIO.Beta_InputsFile,'BetaConfigFile',ModName,Msg)
    return GalaxyIO

def BetaConfigInputs(GalaxyIO,DataCube,TiltedRing,Profiles):
    #   Only use these configs if the beta config file is present
    if GalaxyIO.Beta_InputsFile != None :
        print("Getting beta configuration options from ", GalaxyIO.Beta_InputsFile)
        BetaIn=importlib.import_module(GalaxyIO.Beta_InputsFile)
        
        #   Write the general message
        Message=" based on "+GalaxyIO.Beta_InputsFile+".py."
        
        #   Check on the cloud mode for the tilted ring model
        TiltedRing.cmode=ModuleVarCheck(TiltedRing.cmode,'cmode',BetaIn,Message)
        
        #   Check on the cloud surface density for the tilted ring model.
        TiltedRing.CloudSurfDens=ModuleVarCheck(TiltedRing.CloudSurfDens,'CloudSurfDens',BetaIn,Message)

        #   Check on the random seed
        TiltedRing.ranseed=ModuleVarCheck(TiltedRing.ranseed,'ranseed',BetaIn,Message)
        
        #   Check on the central positions
        TiltedRing.ra_center=ModuleVarCheck(TiltedRing.ra_center,'RA',BetaIn,Message)
        TiltedRing.dec_center=ModuleVarCheck(TiltedRing.dec_center,'DEC',BetaIn,Message)
        print("beta config Ra", TiltedRing.ra_center, TiltedRing.dec_center)
        
        #   Check on the systemic, radial and vertical velocites
        TiltedRing.vsys=ModuleVarCheck(TiltedRing.vsys,'vsys',BetaIn,Message)
        TiltedRing.v_radial=ModuleVarCheck(TiltedRing.v_radial,'vrad',BetaIn,Message)
        TiltedRing.v_vertical=ModuleVarCheck(TiltedRing.v_vertical,'vvert',BetaIn,Message)
        
        #   Check on the vertical gradiant of the tangential velocity
        TiltedRing.dvdz=ModuleVarCheck(TiltedRing.dvdz,'dvdz',BetaIn,Message)
        
        #   Check on the number of bins in RHI and the extent over which to calculate the profiles
        Profiles.nBinsPerRHI=ModuleVarCheck(Profiles.nBinsPerRHI,'nBinsPerRHI',BetaIn,Message)
        Profiles.limR_RHI=ModuleVarCheck(Profiles.limR_RHI,'limR_RHI',BetaIn,Message)
        
        #   Check on the cube shape
        DataCube.cube_shape=ModuleVarCheck(DataCube.cube_shape,'cube_shape',BetaIn,Message)
        
        #   Check on the beam flattening and position angle
        DataCube.beam_flattening=ModuleVarCheck(DataCube.beam_flattening,'Beam_Flattening',BetaIn,Message)
        DataCube.beam_angle=ModuleVarCheck(DataCube.beam_angle,'Beam_PositionAngle',BetaIn,Message)
        
        #   Check on whether the size of the sigma smoothing should change
        DataCube.nSigma=ModuleVarCheck(DataCube.nSigma,'nSigma',BetaIn,Message)
        #   Check on whether to apply velocity smoothing and the value of sigma
        DataCube.velocity_smoothing_switch=ModuleVarCheck(DataCube.velocity_smoothing_switch,'velocity_smooth_switch',BetaIn,Message)
        DataCube.velocity_smooth_sigma=ModuleVarCheck(DataCube.velocity_smooth_sigma,'velocity_smooth_sigma',BetaIn,Message)

        return
    else:
        return

def GlobalVarCheck(Var1,Var2Str,ErrorMsg,ExitCheck):
    try:
        Var1=globals()[Var2Str]
    except:
        print(ErrorMsg)
        if ExitCheck:
            exit()
    return Var1

def ModuleVarCheck(Var1,Var2Str,Module,Msg):
    try:
        Var1=vars(Module)[Var2Str]
        print("Setting " + Var2Str +" = " + str(Var1) +Msg)
    except:
        pass
    return Var1

def ModuleVarCheck_Exit(Var1,Var2Str,Module,Msg):
    try:
        Var1=vars(Module)[Var2Str]
    except:
        print("No value for " + Var2Str +". Set " +Var2Str +Msg)
        exit()
    return Var1





