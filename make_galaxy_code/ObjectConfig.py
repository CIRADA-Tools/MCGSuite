#!/usr/bin/env python3
import numpy as np
import importlib

from . import ObjectDefinitions as OD


def SuiteConfig(Suite):
    #   Create an array of numbers to keep track of the different realizations of each object
    Suite.num_array=np.linspace(0,Suite.num_realizations-1,Suite.num_realizations)
    #   Create an array that consists of the full catalogue of galaxy parameter combinations
    if Suite.UDG_switch==False:
        Suite.CatalogueArray=np.array(np.meshgrid(Suite.mass_array, Suite.beams_array,Suite.inclination_array,Suite.pa_array,Suite.veldisp_array,Suite.num_array)).T.reshape(-1,6)
    elif Suite.UDG_switch:
        Suite.CatalogueArray=np.array(np.meshgrid(Suite.mass_array, Suite.beams_array,Suite.inclination_array,Suite.pa_array,Suite.veldisp_array,Suite.num_array,Suite.VHI_array)).T.reshape(-1,7)
    
    #   Get the total number of galaxies that will be made
    Suite.n_galaxies=np.shape(Suite.CatalogueArray)[0]

    Suite.DBTable=[None]*Suite.n_galaxies

    MaxCubeMemoryEstimate=200*200*200*4/1.e9        #The maximum cube size estimate assuming 200x200x200 cube
    SuiteMaxMemoryEstimate=MaxCubeMemoryEstimate*3*Suite.n_galaxies
    print("The estimated maximum memory usage (assuming 200x200x200 cube): ",\
          SuiteMaxMemoryEstimate, " Gb")

    if SuiteMaxMemoryEstimate > 10:
        Check=input("Warning: This suite may store more than 10 Gb of memory. Enter 'Y' to proceed ")
        if Check != 'Y' and Check !='y':
            print ("Exiting the program")
            exit()


            
def RandomSuiteConfig(Suite):

    Suite.n_galaxies=Suite.SuiteDict['nTot']
    #   Create an array that consists of the full catalogue of galaxy parameter combinations
    DictList=['MassDict','BeamDict','IncDict','PADict','VelDispDict']
    if Suite.UDG_switch:
        DictList.append('VHIDict')
        VHIs=[None]*Suite.n_galaxies
    Suite.CatalogueArray=[None]*Suite.n_galaxies
    Mass=[None]*Suite.n_galaxies
    Beams=[None]*Suite.n_galaxies
    Inc=[None]*Suite.n_galaxies
    PA=[None]*Suite.n_galaxies
    veldisp=[None]*Suite.n_galaxies
    for i in range(Suite.n_galaxies):
        for x in DictList:
            if x=='MassDict':
                Mass[i]=Set_Specific_Param(Suite.SuiteDict[x])
            elif x=='BeamDict':
                Beams[i]=Set_Specific_Param(Suite.SuiteDict[x])
            elif x=='IncDict':
                Inc[i]=Set_Specific_Param(Suite.SuiteDict[x])
            elif x=='PADict':
                PA[i]=Set_Specific_Param(Suite.SuiteDict[x])
            elif x=='VelDispDict':
                veldisp[i]=Set_Specific_Param(Suite.SuiteDict[x])
            elif x=='VHIDict':
                VHIs[i]=Set_Specific_Param(Suite.SuiteDict[x])
        if Suite.UDG_switch:
            Suite.CatalogueArray[i]=[Mass[i],Beams[i],Inc[i],PA[i],veldisp[i],VHIs[i]]
        else:
            Suite.CatalogueArray[i]=[Mass[i],Beams[i],Inc[i],PA[i],veldisp[i]]
    if Suite.UDG_switch==False:
        Suite.CatArrDict={'Mass':Mass,'Beam':Beams,'Inc':Inc,'PA':PA}
    else:
        Suite.CatArrDict={'Mass':Mass,'Beam':Beams,'Inc':Inc,'PA':PA,'VHI':VHIs}

 
    print("Catalogue Array shape", np.shape(Suite.CatalogueArray))
    Suite.DBTable=[None]*Suite.n_galaxies

    MaxCubeMemoryEstimate=200*200*200*4/1.e9        #The maximum cube size estimate assuming 200x200x200 cube
    SuiteMaxMemoryEstimate=MaxCubeMemoryEstimate*3*Suite.n_galaxies
    print("The estimated maximum memory usage (assuming 200x200x200 cube): ",\
          SuiteMaxMemoryEstimate, " Gb")

    if SuiteMaxMemoryEstimate > 10:
        Check=input("Warning: This suite may store more than 10 Gb of memory. Enter 'Y' to proceed ")
        if Check != 'Y' and Check !='y':
            print ("Exiting the program")
            exit()

def Set_Specific_Param(ParamDict):
    if ParamDict['Flag']:
        Val=GetRanVal(ParamDict['Range'],ParamDict['RanType'])
    else:
        Val=ParamDict['FixedVal']
    return Val

def GetRanVal(Limits,RanType):
    if RanType==0:
        Del=Limits[1]-Limits[0]
        Val=np.random.rand()*Del + Limits[0]
    elif RanType==1:
        LogLow=np.log10(Limits[0])
        LogHigh=np.log10(Limits[1])
        Del=LogHigh-LogLow
        Val=np.random.rand()*Del + LogLow
        Val=10.**Val
    else:
        print("No valid random selection type picked")
        exit()
    return Val

def ConfigObjects(Galaxy,DataCube,TiltedRing,Profiles,GalaxyIO):

    DataCube=DataCubeConfig(DataCube,Galaxy.nBeams,Galaxy.VHI,Galaxy.veldisp,TiltedRing)
    
    Profiles=ProfilesConfig(Profiles)
    
    TiltedRing=TiltedRingConfig(TiltedRing,Galaxy.nBeams,DataCube.beam_fwhm)
    
    GalaxyIO=BasicIOConfig(GalaxyIO,Galaxy.nBeams,Galaxy.logMHI,Galaxy.inclination,Galaxy.pa,Galaxy.veldisp,Galaxy.v_HI,Galaxy.version_number,Galaxy.ID,Galaxy.UDG_switch,DataCube.noise)
    
    return DataCube,TiltedRing,Profiles,GalaxyIO


def BasicIOConfig(GalaxyIO,beams, mass, inc, pa,veldisp,v_HI,Version=None,ID=None,UDG_switch=False,noise=None):
    #   Name the Output galaxy
    GalaxyIO.GalaxyName="ba_"+str(beams)+".mass_"+str(round(mass,5))+".inc_"+str(inc)+".pa_"+ str(pa)+".veldisp_"+str(veldisp)+".noise_"+str(round(noise,3))
    #       If there's a version number, include that in the name
    if Version !=None:
        GalaxyIO.GalaxyName="ba_"+str(beams)+".mass_"+str(round(mass,5))+".inc_"+str(inc)+".pa_"+ str(pa)+".veldisp_"+str(veldisp)+".noise_"+str(round(noise,3))+ ".version_"+str(Version)
    if UDG_switch:
        GalaxyIO.GalaxyName+=".UDG_True.v_HI_"+str(v_HI)
    #   Name the diagnostic moment maps plot
    GalaxyIO.MapPlotName=GalaxyIO.GalaxyName+"_MomentMaps.png"
        #   Name the diagnostic profiles plot
    GalaxyIO.ProfilePlotName=GalaxyIO.GalaxyName+"_Profiles.png"
        #   Name the profile plot
    GalaxyIO.ProfileFile=GalaxyIO.GalaxyName+"_Profile.txt"

    #   Also, if there's a version number, name the MCG input files differently
    if ID !=None:
        GalaxyIO.MCG_Main_InputFile="MCGMain_"+str(ID)+".in"
        #   Set the default name for the MCG datacube file
        GalaxyIO.MCG_DataCube_InputFile="MCG_DataCube_"+str(ID)+".in"
        #   Set the default name for the MCG tilted ring file
        GalaxyIO.MCG_TiltedRing_InputFile="MCG_TiltedRing_"+str(ID)+".in"
    
    return GalaxyIO


def ProfilesConfig(Profiles):
    Profiles.step = 1 / ( Profiles.nBinsPerRHI*Profiles.limR_RHI )
    #   Allocate an array galactocentric radii (kpc)
    Profiles.R=np.arange(0.,(Profiles.limR_RHI + Profiles.step),Profiles.step)
        #   Allocate an array of zeros for the velocity (km/s)
    Profiles.Vrot=np.zeros_like(Profiles.R)
        #   Allocate an array of zeros for the surface brightness
    Profiles.SB=np.zeros_like(Profiles.R)
    return Profiles

def DataCubeConfig(DataCube,nBeams,VHI,veldisp,TiltedRing):
     #   Set the cube dimensions (in arcsec, arcsec, and km/s)
    DataCube.cube_dimensions=[-DataCube.pixel_size,DataCube.pixel_size,DataCube.channel_size]
    
    #   Now check if the data cube shape has been defined.  If it has use that shape.
    if DataCube.cube_shape == None:
        #   If it hasn't set the cube size to be at least 10 beams across of the size+margin across
        if nBeams+DataCube.cube_margin <= DataCube.cubesize_minimum:
            SpatialLim=DataCube.cubesize_minimum
        else:
            SpatialLim=nBeams+DataCube.cube_margin
        #   Convert the spatial lim from units of beams to units of pixels
        SpatialLim=int(SpatialLim*DataCube.beam_fwhm/DataCube.pixel_size)+1
        #   Set the spectral limit to be 2* (VHI+1.5*v_disp)+margin in channels
        SpectralLim=int((2*(VHI+1.5*veldisp)+DataCube.velocity_margin)/DataCube.channel_size)+1
        #   Set the cubes shape
        DataCube.cube_shape=[SpatialLim,SpatialLim,SpectralLim]
        #   Calculate the beam dimensions
    DataCube.beam_dimensions=[DataCube.beam_fwhm,DataCube.beam_fwhm*DataCube.beam_flattening,DataCube.beam_angle]
        #   reference pixel location (CRPIX in normal FITS files)
    DataCube.reference_locations=[int(DataCube.cube_shape[0]/2),int(DataCube.cube_shape[1]/2),int(DataCube.cube_shape[2]/2)]
        #   Default reference pixel values (CRVAL in normal FITS files) -- will be recalculated when the centers of the cube are known.
    DataCube.reference_values=[TiltedRing.ra_center,TiltedRing.dec_center,TiltedRing.vsys]
    print("DataCube Config", TiltedRing.ra_center)
    
    #   Set the velocity smoothing to be equal to the channel width
    if DataCube.velocity_smooth_sigma == None:
        DataCube.velocity_smooth_sigma=DataCube.cube_dimensions[2]
    return DataCube

def TiltedRingConfig(TiltedRing,nBeams,beam_fwhm):
    #   Calculate the number of rings based on the number of rings/beam and how many RHI to extend the calculation to.
    TiltedRing.nRings=int((TiltedRing.RHI_Limit*(nBeams))*TiltedRing.nRingsPerBeam)+1
    #   Default ring width (arcsec) (based on the number of rings/beam
    TiltedRing.Rwidth=beam_fwhm/TiltedRing.nRingsPerBeam
    #   Default array of zeros for the radii -- values calculated in profiles (arsec).
    TiltedRing.R_array=np.zeros(TiltedRing.nRings)
    #   Default array of zeros for the tangential velocity --values calculated in profiles (km/s).
    TiltedRing.v_tangential_array=np.zeros(TiltedRing.nRings)
    #   Default array of zeros for the surface brightness -- valuess calculated in profiles (Jy km/s arcsec^-2).
    TiltedRing.sigma_array=np.zeros(TiltedRing.nRings)
    #   Default array of zeros for the scale height in observed coordinates(arcsec)
    TiltedRing.z_scale=np.zeros(TiltedRing.nRings)
    #   Default velocity gradient starting height in observed coordinates (arcsec)
    TiltedRing.z_gradiantstart=5.*TiltedRing.z_scale
    return TiltedRing

