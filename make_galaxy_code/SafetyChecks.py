#!/usr/bin/env python3
import numpy as np
from . import ObjectDefinitions as OD

ErrorMsg="Ending program"

def SuiteChecks(Suite,DataCube):
    #   Check on the data cube first as this will be constant for each galaxy
    DataCubeChecks(DataCube)
    #   Next do the suite parameter checks
    SuiteLimitChecks(Suite)

def FirstChecks(Galaxy,DataCube):
    #   Check the Galaxy
    GalaxyChecks(Galaxy)
    #   Next check the data cube
    DataCubeChecks(DataCube)


def SuiteLimitChecks(Suite):
    #   Check the mass array
    MassMsgStrings=['mass', 'Masses', 'suite_config_MCG.py']
    ArrayCheck(Suite.mass_array,OD.Galaxy().MassLimits,MassMsgStrings)
    
    #   Check the beam size array
    SizeMsgStrings=['beam size', 'Beams', 'suite_config_MCG.py']
    ArrayCheck(Suite.beams_array,OD.Galaxy().SizeLimits,SizeMsgStrings)
    
    #   Check the velocity dispersion array
    DispersionMsgStrings=['velocity dispersion', 'veldisps', 'suite_config_MCG.py']
    ArrayCheck(Suite.veldisp_array,OD.Galaxy().dispersionLimits,DispersionMsgStrings)

    #   Check the inclinations array
    InclinationMsgStrings=['inclination', 'Inclinations', 'suite_config_MCG.py']
    ArrayCheck(Suite.inclination_array,[0.,90.],InclinationMsgStrings)

    #   Check the position angles array
    PositionAngleMsgStrings=['position angle', 'PositionAngles', 'suite_config_MCG.py']
    ArrayCheck(Suite.pa_array,[0.,360.],PositionAngleMsgStrings)



def ArrayCheck(Arr,Lims,MsgStrings):
    #   Check that all elements of the array are between the acceptable limits
    for i in range(len(Arr)):
        if Arr[i] < Lims[0] or Arr[i] > Lims[1]:
            print("The "+str(i+1)+" element of the " +MsgStrings[0]+" array is out of the acceptable range:", Arr[i])
            print("Set all elements of "+ MsgStrings[1]+"in "+MsgStrings[2] +"are between:", Lims )
            print(ErrorMsg)
            exit()
    #   Check that all elements of the array are unique.
    if len(np.unique(Arr)) != len(Arr):
        print("Some elements of the "+ MsgStrings[0] +" array are repeated.  Make sure all elements of "+ MsgStrings[1] +" are unique in " + MsgStrings[2])
        print(ErrorMsg)
        exit()


def GalaxyChecks(Galaxy):
    #   Check the size of the galaxy in beams
    if Galaxy.nBeams < Galaxy.SizeLimits[0] or Galaxy.nBeams > Galaxy.SizeLimits[1]:
        print("Galaxy size in beams is out of acceptable range:", Galaxy.nBeams)
        print("Set the nBeams in galaxyconfig_MCG between:", Galaxy.SizeLimits)
        print(ErrorMsg)
        exit()

    if Galaxy.logMHI < Galaxy.MassLimits[0] or Galaxy.logMHI > Galaxy.MassLimits[1]:
        print("Galaxy HI mass is out of acceptable range:", Galaxy.logMHI)
        print("Set logMHI in galaxyconfig_MCG between:", Galaxy.MassLimits)
        print(ErrorMsg)
        exit()

    if Galaxy.veldisp < Galaxy.dispersionLimits[0] or Galaxy.veldisp > Galaxy.dispersionLimits[1]:
        print("Galaxy velocity dispersion is out of acceptable range:", Galaxy.veldisp)
        print("Set veldisp in galaxyconfig_MCG between:", Galaxy.dispersionLimits)
        print(ErrorMsg)
        exit()


def DataCubeChecks(DataCube):
    #   Check the beam/pixel size ratio
    BeamPixelRatio=DataCube.beam_fwhm/DataCube.pixel_size
    if BeamPixelRatio < DataCube.BeamPixelRatioLimits[0] or BeamPixelRatio > DataCube.BeamPixelRatioLimits[1]:
        print("The beam/pixel size ratio is out of the acceptable range: ",BeamPixelRatio)
        print("Adjust the values of PixelSize and beam_fwhm in observatory_config_MCG to have a ratio between",  DataCube.BeamPixelRatioLimits)
        exit()
    #   Check the channel size
    if DataCube.channel_size < DataCube.ChannelSizeLimits[0] or DataCube.channel_size > DataCube.ChannelSizeLimits[1]:
            print("The channel size is out of the acceptable range: ",DataCube.channel_size)
            print("Set the ChannelSize in observatory_config_MCG between:", DataCube.ChannelSizeLimits)
            exit()


