#!/usr/bin/env python3
import numpy as np
import importlib
import astropy
from astropy.io import fits
import gc


def GetProfilesAndMaps(GalaxyIO,PositionAngle,BeamFWHM,noise,dist):
    #   Get the profiles and maps from a data cube
    #INPUT-->   GalaxyIO -- the Galaxy IO object defined in ObjectDefinitions
    #OUTPUT-->  Velocities -- a 1D array of velocities from the cube
    #           VelocityProfile -- a 1D array of fluxes from the profiles
    #           Moment0 --The moment 0 map
    #           Moment1 -- The moment 1 map
    #           Moment2 -- The moment 2 map
    
    #   Get the Cube
    CubeFile=GalaxyIO.GalaxyName+"/"+GalaxyIO.GalaxyName+".fits"
    FullCube = fits.open(CubeFile)
    FullData=FullCube[0].data
    header=FullCube[0].header
    
    #   Get the noiseless convolved cube to act as a mask
    NoiselessFile=GalaxyIO.GalaxyName+"/"+GalaxyIO.GalaxyName+"_ConvolvedSourceCube.fits"
    #NoiselessFile=GalaxyIO.GalaxyName+"/"+GalaxyIO.GalaxyName+"_SourceCube.fits"
    NoiselessCube = fits.open(NoiselessFile)
    NoiselessData=NoiselessCube[0].data
    #   Use the noiseless cube as a mask
    MaskedFullData=SimpleMask(FullData,NoiselessData,noise)
    print("convolved data flux check", np.nansum(NoiselessData),np.nansum(FullData))
    
    Velocities=GetVels(header)
    PX,PY=GetPixels(header)
    BS=1./2.*(BeamFWHM/((PX[1]-PX[0])*3600.))
    
    FullResult=GeneralProfileAndMaps(Velocities,MaskedFullData)
    NoiselessResult=GeneralProfileAndMaps(Velocities,NoiselessData)
    PVMajorData,PVMinorData=PVCalculatorV2(MaskedFullData,Velocities,PositionAngle,BS)
    PVMajorSource,PVMinorSource=PVCalculatorV2(NoiselessData,Velocities,PositionAngle,BS)
    
    
    PixSize=[PX[1]-PX[0],PY[1]-PY[0]]
    ChannelSize=Velocities[1]-Velocities[0]
    MTot=MassCalc(NoiselessData,PixSize,BeamFWHM,dist,ChannelSize)
    print("Total Mass Check:", MTot,np.log10(MTot))
    print((10.**9.6-MTot)/10**9.6)
    
    #FRot=SpatialRotate(NoiselessData,PositionAngle)
    #NoiselessResult=GeneralProfileAndMaps(Velocities,FRot)
    
    
    FullCube.close()
    NoiselessCube.close()
    
    del FullCube, FullData, NoiselessCube,NoiselessData
    gc.collect()
    
    return PX,PY,Velocities,FullResult[0],FullResult[1],FullResult[2],FullResult[3],NoiselessResult[0],NoiselessResult[1],NoiselessResult[2],NoiselessResult[3],PVMajorData,PVMinorData,PVMajorSource,PVMinorSource

def MassCalc(cube,PixSize,beam_size,dist,ChannelSize):
    #print(PixSize,beam_size,dist,ChannelSize)
    beamarea=(np.pi*beam_size**2.)/(4.*np.log(2.))
    pixperbeam=beamarea/(abs(PixSize[0]*3600.)*abs(PixSize[0])*3600.)
    totalsignal = np.sum(cube)/pixperbeam
    Mtest1 = 0.236*(dist*1000.)**2*totalsignal*ChannelSize
    #print(Mtest1,np.log10(Mtest1),totalsignal,pixperbeam,beamarea)
    return Mtest1


def GeneralProfileAndMaps(Velocities,Data):
    VelocityProfile=Calculate1DProfile(Data)
    Moment0,Moment1,Moment2=RadioMomentMapCalc(Data,Velocities)
    
    return VelocityProfile,Moment0,Moment1,Moment2


def SimpleMask(Arr1,Arr2,noise):
    
    MaskLim=0.5*noise/1000.
    #MaskLim=1.e-5
    TestArr=np.copy(Arr1)
    channels=np.shape(Arr1)[0]
    rows=np.shape(Arr1)[1]
    cols=np.shape(Arr1)[2]
    for i in range(channels):
        for j in range(rows):
            for k in range(cols):
                if Arr2[i,j,k] <= MaskLim:
                    TestArr[i,j,k]=0.

    return TestArr


def GetVels(hdr):
    #   Get a 1D array of velocities for use in the moment calculations
    #INPUT--> hdr -- a fits file header
    #OUTPUT--> Vels -- a 1D array of velocities corresponding to each channel in the fits file
    Vels=np.zeros(hdr['NAXIS3'])
    RefChannel=hdr['CRPIX3']
    RefVel=hdr['CRVAL3']/1000.
    dV=hdr['CDELT3']/1000.
    for i in range(hdr['NAXIS3']):
        Vels[i]=(i-RefChannel+1)*dV+RefVel
    return Vels

def GetPixels(hdr):
    #   Get a 1D array of velocities for use in the moment calculations
    #INPUT--> hdr -- a fits file header
    #OUTPUT--> PixelsX -- a 1D array of velocities corresponding to each channel in the fits file
    PixelsX=np.zeros(hdr['NAXIS1'])
    PixelsY=np.zeros(hdr['NAXIS2'])
    RefX=hdr['CRPIX1']
    RefY=hdr['CRPIX2']
    RefValX=hdr['CRVAL1']
    RefValY=hdr['CRVAL2']
    dX=hdr['CDELT1']
    dY=hdr['CDELT2']
    for i in range(hdr['NAXIS1']):
        PixelsX[i]=(i-RefX+1)*dX+RefValX
    for i in range(hdr['NAXIS2']):
        PixelsY[i]=(i-RefY+1)*dY+RefValY
        #PixelsX[i]=(i-RefX+1)*dX
        #PixelsY[i]=(i-RefY+1)*dY
#    print(PixelsX)
#    print(PixelsY)
#    print(dX,dY,RefValX,RefValY)
    return PixelsX,PixelsY

def Calculate1DProfile(Fluxes):
    #   Get a 1D spectrum for the cube
    #INPUT--> Fluxes -- a 3D array with channels as the first dimension
    #OUTPUT-->Profile -- a 1D array with the summed up flux of the profile
    Profile=np.zeros(np.shape(Fluxes)[0])
    for i in range(np.shape(Fluxes)[0]):
        Profile[i]=np.sum(Fluxes[i,:,:])
    return Profile


def RadioMomentMapCalc(Fluxes,Velocities):
    #   Initialize the maps to zero
    Map0=np.zeros([np.shape(Fluxes)[1],np.shape(Fluxes)[2]])
    Map1=np.zeros([np.shape(Fluxes)[1],np.shape(Fluxes)[2]])
    Map2=np.zeros([np.shape(Fluxes)[1],np.shape(Fluxes)[2]])
    #   First get the Moment 0 and Moment 1 maps
    for i in range(np.shape(Fluxes)[1]):
        for j in range(np.shape(Fluxes)[2]):
            Map1[i,j]=np.sum(Fluxes[:,i,j]*(Velocities[:]))
            Map0[i,j]=np.sum(Fluxes[:,i,j])
            Map1[i,j]=Map1[i,j]/Map0[i,j]
            #   Make sure the flux isn't zero
            if Map0[i,j] <= 1.e-5:
                Map1[i,j]=np.nan
    #   Now get the Moment 2 map using the Moment 1 map
    for i in range(np.shape(Fluxes)[1]):
        for j in range(np.shape(Fluxes)[2]):
            Map2[i,j]=np.sum(Fluxes[:,i,j]*(Velocities[:]-Map1[i,j])**2.)
            #   Make sure the flux isn't zero
            if Map0[i,j] <= 1.e-5:
                Map2[i,j]=np.nan
    Map2=np.sqrt(Map2)
    return Map0,Map1,Map2


def PVCalculatorV2(Fluxes,Velocities,PositionAngle,SliceThickness):
    #   Adjust the position angle so that 0 points in the y direction
    paRad_Major=(90-PositionAngle)*np.pi/180.
    paRad_Minor=paRad_Major+np.pi/2.
    #   Set the dimensions of the PV digram
    PVDimensions=[np.shape(Fluxes)[0],np.shape(Fluxes)[1]]
    #   Get the value of the central pixel
    CentPix=[int(np.shape(Fluxes)[1]/2)-1,int(np.shape(Fluxes)[2]/2)-1]
    #   Initialize the PV diagrams to zero.
    PVMajor=np.zeros(PVDimensions)
    PVMinor=np.zeros(PVDimensions)
    #   Set the thickness of the slice
    SliceThicknessInt=int(SliceThickness)+1
    
    for i in range(PVDimensions[1]):
        #   Get x for the major diagram and yy for the minor
        x=float(i-CentPix[0])
        yy=float(i-CentPix[1])
        for j in range(-SliceThicknessInt,SliceThicknessInt+1):
            #   Get y for the major diagram and xx for the minor
            y=float(j)
            xx=float(j)
            #   Rotate the points from PV space to the observed orientation.
            xP=x*np.cos(paRad_Major)-y*np.sin(paRad_Major)
            yP=x*np.sin(paRad_Major)+y*np.cos(paRad_Major)
            #   Get the correct indices for the 'minor' axis point in the observed orientation.
            k=int(round(xP+CentPix[0]))
            l=int(round(yP+CentPix[1]))
            #   Make sure the points are inside the observed image
            if k >= 0 and k < np.shape(Fluxes)[1]:
                if l >= 0 and l < np.shape(Fluxes)[1]:
                    #   Get the flux for each channel
                    for m in range(PVDimensions[0]):
                        PVMinor[m,i]+=Fluxes[m,k,l]
            #   Repeat the process for the major axis
            xP=xx*np.cos(paRad_Major)-yy*np.sin(paRad_Major)
            yP=xx*np.sin(paRad_Major)+yy*np.cos(paRad_Major)
            k=int(round(xP+CentPix[0]))
            l=int(round(yP+CentPix[1]))
            if k >= 0 and k < np.shape(Fluxes)[1]:
                if l >= 0 and l < np.shape(Fluxes)[1]:
                    for m in range(PVDimensions[0]):
                        PVMajor[m,i]+=Fluxes[m,k,l]

    return PVMajor,PVMinor


def SpatialRotate(Fluxes,PositionAngle):
    
    paRad_Major=(-90+PositionAngle)*np.pi/180.
    NewFlux=np.zeros(np.shape(Fluxes))
    CentPix=[int(np.shape(Fluxes)[1]/2),int(np.shape(Fluxes)[2]/2)]
    
    for i in range(np.shape(Fluxes)[1]):
        for j in range(np.shape(Fluxes)[2]):
            x=float(i-CentPix[0])
            y=float(j-CentPix[1])
            xP=x*np.cos(paRad_Major)-y*np.sin(paRad_Major)
            yP=x*np.sin(paRad_Major)+y*np.cos(paRad_Major)
            kF=xP+CentPix[0]
            lF=yP+CentPix[1]
            k=int(round(kF))
            l=int(round(lF))
            if k >=0 and k < np.shape(Fluxes)[1]:
                if l >=0 and l < np.shape(Fluxes)[1]:
                    for m in range(np.shape(Fluxes)[0]):
                        NewFlux[m,k,l]=Fluxes[m,i,j]
    return NewFlux

