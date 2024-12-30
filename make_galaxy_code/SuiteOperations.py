#!/usr/bin/env python3
import numpy as np
import multiprocessing as mp
import copy as copy
import os

from . import ObjectDefinitions as OD
from . import Inputs as IN
from . import ObjectConfig as OC
from . import SafetyChecks as SC
from . import MakeGalaxy as MG
from . import MakeProfiles as MP
from . import MakeTiltedRing as MTR
from . import MakeCubes as MC
from . import DiagnosticPlots as DP
from . import OutputClean as OutClean


import resource

import gc

def SuiteMainLoopFn(i,Suite):
#   This function runs through the loop needed to create a specific galaxy from the
#   suite of input parameters
#---->  INPUT:  i == the step in the suite.
#               Suite == Suite object containing the array of parameters.
#---->  OUTPUT: DBEntry == An entry for the final database.

    #   Keep track of the current iteration, the processor, and galaxy parameters to be used.
    current = mp.current_process()
    print("Step and Processor ID",i,Suite.CatalogueArray[i],current._identity)

    #   Create the Galaxy object from the suite
    Galaxy=CreateGalaxyInstance(i,Suite)
    #   Copy the suite templates for the profiles, data cube, tilted ring, and SuiteIO objects
    Profiles=copy.deepcopy(Suite.Templates[0])
    DataCube=copy.deepcopy(Suite.Templates[1])
    TiltedRing=copy.deepcopy(Suite.Templates[2])
    GalaxyIO=copy.deepcopy(Suite.SuiteIO)
  
    #   Calculate all the various galaxy parameters from the HI Mass
    Galaxy=MG.MakeGalaxy(Galaxy,DataCube)
    
    #   Configure the various objects
    DataCube,TiltedRing,Profiles,GalaxyIO=OC.ConfigObjects(Galaxy,DataCube,TiltedRing,Profiles,GalaxyIO)
    #   Calculate the profiles based on the galaxy parameters
    Profiles=MP.MakeProfiles(Profiles,Galaxy,GalaxyIO)
    #   Calculate the full tilted ring model.
    TiltedRing=MTR.MakeTiltedRing(Galaxy,DataCube,TiltedRing)
    #   Make the cubes (currently using MCG)
    MC.MakeCubes(GalaxyIO,DataCube,TiltedRing)
    #   Make all Plots
    DP.MakeAllPlots(GalaxyIO,Galaxy,DataCube,Profiles,TiltedRing)
    #   Clean up outputs
    OutClean.CleanOutput(GalaxyIO)

    #   Write out the database entry for this galaxy
    DBEntry=DatabaseEntries(i,Galaxy,GalaxyIO)

    #   Remove from memory the Galaxy, DataCube, TiltedRing, Profiles, and GalaxyIO instances
    del Galaxy, Profiles,DataCube,TiltedRing,GalaxyIO
    #   Make sure the memory is cleared
    gc.collect()

    #   Track the RAM memory usage for potential issues.
    MemCheck=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print("Process memory check", i,current._identity, MemCheck/1e9)

    return DBEntry

def CreateGalaxyInstance(step,Suite):
#   Generates the intial galaxy object for an element from the suite - catalogue array

    #   Initialize the galaxy object
    Galaxy=OD.Galaxy()
    #   Set the galaxy parameters to the correct element from the catalogue array
    Galaxy.logMHI=Suite.CatalogueArray[step,0]
    Galaxy.nBeams=Suite.CatalogueArray[step,1]
    Galaxy.inclination=Suite.CatalogueArray[step,2]
    Galaxy.pa=Suite.CatalogueArray[step,3]
    Galaxy.veldisp=Suite.CatalogueArray[step,4]
    Galaxy.version_number=int(Suite.CatalogueArray[step,5])
    Galaxy.ID=int(step)
    Galaxy.UDG_switch=Suite.UDG_switch
    if Suite.UDG_switch:
        Galaxy.v_HI=Suite.CatalogueArray[step,6]
        #print(Suite.VHI_array)
        #print("UDG COnfig Calc Checks", Galaxy.v_HI)
        #print(Suite.CatalogueArray[step,:])

    return Galaxy
    

def SuiteMainLoopFn_random(i,Suite):
#   This function runs through the loop needed to create a specific galaxy from the
#   suite of input parameters
#---->  INPUT:  i == the step in the suite.
#               Suite == Suite object containing the array of parameters.
#---->  OUTPUT: DBEntry == An entry for the final database.

    #   Keep track of the current iteration, the processor, and galaxy parameters to be used.
    current = mp.current_process()
    print("Step and Processor ID",i,Suite.CatalogueArray[i],current._identity)
    print("Suite check", Suite.UDG_switch)

    #   Copy the suite templates for the profiles, data cube, tilted ring, and SuiteIO objects
    Profiles=copy.deepcopy(Suite.Templates[0])
    DataCube=copy.deepcopy(Suite.Templates[1])
    TiltedRing=copy.deepcopy(Suite.Templates[2])
    GalaxyIO=copy.deepcopy(Suite.SuiteIO)
    
    #   Create the Galaxy object from the suite
    Galaxy,TiltedRing,GalaxyIO=CreateGalaxyInstance_random(i,Suite,TiltedRing,GalaxyIO)
  
    #   Calculate all the various galaxy parameters from the HI Mass
    Galaxy=MG.MakeGalaxy(Galaxy,DataCube)
    
    #   Configure the various objects
    DataCube,TiltedRing,Profiles,GalaxyIO=OC.ConfigObjects(Galaxy,DataCube,TiltedRing,Profiles,GalaxyIO)
    #   Calculate the profiles based on the galaxy parameters
    Profiles=MP.MakeProfiles(Profiles,Galaxy,GalaxyIO)
    #   Calculate the full tilted ring model.
    TiltedRing=MTR.MakeTiltedRing(Galaxy,DataCube,TiltedRing)
    #   Make the cubes (currently using MCG)
    MC.MakeCubes(GalaxyIO,DataCube,TiltedRing)
    #   Make all Plots
    DP.MakeAllPlots(GalaxyIO,Galaxy,DataCube,Profiles,TiltedRing)
    #   Clean up outputs
    OutClean.CleanOutput(GalaxyIO)

    #   Write out the database entry for this galaxy
    DBEntry=DatabaseEntries(i,Galaxy,GalaxyIO)

    #   Remove from memory the Galaxy, DataCube, TiltedRing, Profiles, and GalaxyIO instances
    del Galaxy, Profiles,DataCube,TiltedRing,GalaxyIO
    #   Make sure the memory is cleared
    gc.collect()

    #   Track the RAM memory usage for potential issues.
    MemCheck=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print("Process memory check", i,current._identity, MemCheck/1e9)

    return DBEntry
    
def CreateGalaxyInstance_random(step,Suite,TiltedRing,GalaxyIO):
#   Generates the intial galaxy object for an element from the suite - catalogue array
    print("Creating random galaxy",Suite.CatalogueArray[step])
    #   Initialize the galaxy object
    Galaxy=OD.Galaxy()
    #   Set the galaxy parameters to the correct element from the catalogue array
    Galaxy.logMHI=Suite.CatalogueArray[step][0]
    Galaxy.nBeams=Suite.CatalogueArray[step][1]
    Galaxy.inclination=Suite.CatalogueArray[step][2]
    Galaxy.pa=Suite.CatalogueArray[step][3]
    Galaxy.veldisp=Suite.CatalogueArray[step][4]
    Galaxy.version_number=0
    Galaxy.ID=int(step)
    Galaxy.UDG_switch=Suite.UDG_switch
    print("Making a random instance")
    print("UDG Switch", Galaxy.UDG_switch)
    if Suite.UDG_switch:
        Galaxy.v_HI=Suite.CatalogueArray[step][5]

    return Galaxy,TiltedRing,GalaxyIO

def DatabaseEntries(step,Galaxy,GalaxyIO):
    DBEntry=[step,GalaxyIO.GalaxyName,Galaxy.logMHI,Galaxy.nBeams,Galaxy.veldisp,
             Galaxy.inclination,Galaxy.pa,\
             Galaxy.RHI,Galaxy.VHI,Galaxy.distance,\
             Galaxy.RC.rPE,Galaxy.RC.vPE,Galaxy.RC.aPE,\
             Galaxy.SB.sbmax, Galaxy.SB.Grmax,Galaxy.SB.Gsig,Galaxy.SB.Er, Galaxy.SB.Vsig]

    return DBEntry

def CatalogueOutput(Suite):
    AsciiName=Suite.SuiteIO.OutputFolder+".txt"
    SQLName=Suite.SuiteIO.OutputFolder+".sql"

    WriteTextCatalogue(Suite,AsciiName)
    WriteSQLCatalogue_Entry(Suite,SQLName)

    os.system("mv "+ AsciiName+ " " + Suite.SuiteIO.OutputFolder)
    os.system("mv "+ SQLName+ " " + Suite.SuiteIO.OutputFolder)


def WriteTextCatalogue(Suite,fName):
    Header=["ID","Name", "Mass", "Beams", "VelocityDispersion", "Inclination", "PositionAngle",\
            "RHI","VHI","Distance","rPE","vPE","aPE","SBmax", "Grmax", "Gsig", "Er", "Vsig"]
    HeaderStr='\t'.join(map(str, Header))+"\n"
    file=open(fName,"w")
    file.write(HeaderStr)

    for i in range(Suite.n_galaxies):
        EntryStr='\t'.join(map(str, Suite.DBTable[i]))+"\n"
        file.write(EntryStr)

    file.close()

def WriteSQLCatalogue_Entry(Suite,fName):
    OriginString="--MCG suite maker (version 0.9.1)\n \n"
    file=open(fName,"w")
    file.write(OriginString)
    #SQLMode='SET SQL_MODE = "NO_AUTO_VALUE_ON_ZERO";'+'\n\n'
    #file.write(SQLMode)
    
    TbleName="MCG-Catalogue"
    TbleCreateStr=" CREATE TABLE IF NOT EXISTS '"+TbleName+"'(\n"
    
    HeaderNames=["'ID'","'Name'", "'Mass'", "'Beams'", "'VelocityDispersion'",\
                 "'Inclination'", "'PositionAngle'",\
                 "'RHI'","'VHI'","'Distance'","'rPE'","'vPE'","'aPE'",\
                 "'SBmax'", "'Grmax'", "'Gsig'", "'Er'", "'Vsig'"]
    HeaderFmt=["int","text","double","double","double","double","double",\
               "double","double","double","double","double","double",\
               "double","double","double","double","double"]
               
               
    HeaderStr=TbleCreateStr
    for i in range(len(HeaderNames)):
        if i ==0:
            HeaderStr+=" "+HeaderNames[i]+" "+HeaderFmt[i]+" NOT NULL PRIMARY KEY,\n"
        elif i ==len(HeaderNames)-1:
            HeaderStr+=" "+HeaderNames[i]+" "+HeaderFmt[i]+" NOT NULL);\n\n"
        else:
            HeaderStr+=" "+HeaderNames[i]+" "+HeaderFmt[i]+" NOT NULL,\n"

#HeaderStr+=" PRIMARY KEY ('ID'),\n  KEY ('ID')\n) DEFAULT CHARSET=utf8 COMMENT=\'MCG-Catalogue\';\n\n"
    file.write(HeaderStr)

    HeaderList=', '.join(map(str, HeaderNames))
    for i in range(Suite.n_galaxies):
        Suite.DBTable[i][1]="'"+Suite.DBTable[i][1]+"'"
        ValueStr=', '.join(map(str, Suite.DBTable[i]))
        EntryStr="INSERT INTO '"+TbleName+"' ("+HeaderList+") VALUES\n("+ValueStr+");\n"
        file.write(EntryStr)

    file.close()


def CatalogueOutput_csv(Suite):
    import pandas as pd
    
    CSVName=Suite.SuiteIO.OutputFolder+".csv"
    
    #print(Suite.CatArrDict)
    Name=[None]*Suite.SuiteDict['nTot']
    ID=[None]*Suite.SuiteDict['nTot']
    RHI=[None]*Suite.SuiteDict['nTot']
    VHI=[None]*Suite.SuiteDict['nTot']
    Distance=[None]*Suite.SuiteDict['nTot']
    rPE=[None]*Suite.SuiteDict['nTot']
    vPE=[None]*Suite.SuiteDict['nTot']
    aPE=[None]*Suite.SuiteDict['nTot']
    SBmax=[None]*Suite.SuiteDict['nTot']
    Grmax=[None]*Suite.SuiteDict['nTot']
    Gsig=[None]*Suite.SuiteDict['nTot']
    Er=[None]*Suite.SuiteDict['nTot']
    Vsig=[None]*Suite.SuiteDict['nTot']
    for i in range(len(Suite.DBTable)):
        #print(Suite.DBTable[i][6],Suite.CatArrDict['PA'][i])
        ID[i]=Suite.DBTable[i][0]
        Name[i]=str(Suite.DBTable[i][1])
        RHI[i]=Suite.DBTable[i][7]
        VHI[i]=Suite.DBTable[i][8]
        Distance[i]=Suite.DBTable[i][9]
        rPE[i]=Suite.DBTable[i][10]
        vPE[i]=Suite.DBTable[i][11]
        aPE[i]=Suite.DBTable[i][12]
        SBmax[i]=Suite.DBTable[i][13]
        Grmax[i]=Suite.DBTable[i][14]
        Gsig[i]=Suite.DBTable[i][15]
        Er[i]=Suite.DBTable[i][16]
        Vsig[i]=Suite.DBTable[i][17]
        
    
    Suite.CatArrDict['ID']=ID
    Suite.CatArrDict['Name']=Name
    Suite.CatArrDict['RHI']=RHI
    Suite.CatArrDict['VHI']=VHI
    Suite.CatArrDict['Distance']=Distance
    Suite.CatArrDict['rPE']=rPE
    Suite.CatArrDict['vPE']=vPE
    Suite.CatArrDict['aPE']=aPE
    Suite.CatArrDict['SBmax']=SBmax
    Suite.CatArrDict['Grmax']=Grmax
    Suite.CatArrDict['Gsig']=Gsig
    Suite.CatArrDict['Er']=Er
    Suite.CatArrDict['Vsig']=Vsig

    
    
    DF=pd.DataFrame.from_dict(Suite.CatArrDict)
    DF.to_csv(CSVName,index=False)
    os.system("mv "+ CSVName+ " " + Suite.SuiteIO.OutputFolder)




