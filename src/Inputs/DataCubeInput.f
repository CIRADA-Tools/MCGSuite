cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines for reading in a
c       datacube.  It also contains routines for getting
c       a text version of a data cube object header.
c       Since the beam parameters are usually contained in the
c       datacube headers, this file also gets those in combined
c       input files
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module DataCubeInputMod

      use DataCubeMod
      use BeamMod
      use CommonConsts
      use CalcBeamKernelMod
      use InputUnitConversionsMod

      implicit none

      contains

ccccccc
c           This routine reads in a text version of a data cube object header
      subroutine ReadTextDataCubeHeader(DC,Beam,HeaderFileName)
      implicit none
      Type(DataCube) DC
      Type(Beam2D) Beam
      character(*) HeaderFileName
      logical FileCheck

      print*, "Getting the data cube input file"
      INQUIRE(file=HeaderFileName,EXIST=FileCheck)
      if(FileCheck .eqv. .False.) then
        print*, "Data cube header file does not exist "
     &          , trim(HeaderFileName)
        stop
      endif

      open(10,file=trim(HeaderFileName),status='old')
c       Get the number of pixels and channels
      read(10,*)
      read(10,*) DC%DH%nPixels(0:1), DC%DH%nChannels
c       Get the units for the datacube sizes, reference positions, beam axis, and beam angle
      read(10,*)
      read(10,*) DC%DH%DimensionUnitSwitch(0:3)
     &                  ,Beam%BeamUnitsSwitch(0:1)
c       Get the dimensions of the pixels and channels
      read(10,*)
      read(10,*) DC%DH%PixelSize(0:1), DC%DH%ChannelSize
c       Get the reference locations of each dimension
      read(10,*)
      read(10,*) DC%DH%RefLocation(0:2)
c       Get the reference values for each dimension
      read(10,*)
      read(10,*) DC%DH%RefVal(0:2)
c      DC%DH%RefVal(0:1)=DC%DH%RefVal(0:1)*3600
c       Get the beam size/shape
      read(10,*)
      read(10,*)  Beam%BeamMajorAxis, Beam%BeamMinorAxis
     &              ,Beam%BeamPositionAngle

c      Beam%BeamPositionAngle=Beam%BeamPositionAngle*Pi/180.


c       Get the number of sigma lengths to calculate the beam to
      read(10,*)
      read(10,*) Beam%SigmaLengths
c       Get the type of velocity smoothing being used
      read(10,*)
      read(10,*) Beam%VelocitySmoothSwitch
c       If using Gaussian smoothing read in the PSF
      if(Beam%VelocitySmoothSwitch .eq. 1) then
        read(10,*)
        read(10,*) Beam%VelocitySmoothSigma
      endif
      close(10)


c       Do the datacube and beam unit conversions
      call TextDCHeader_AngularConversions(DC,Beam)
c      print*, "Cube Vel Dimensions", DC%DH%ChannelSize
c     &          , DC%DH%RefVal(2)


c       Finally allocate the data cube and the beam
      call AllocateDataCube(DC)     !/src/ObjectDefinitions/DataCube.f
c               It is necessary to set the beam pixel size to datacube pixel size
      Beam%PixelSize=DC%DH%PixelSize
      call Allocate_Beam2D(Beam,DC%DH%nPixels)   !/src/ObjectDefinitions/Beam.f
c       Also calculate the real kernel for the beam
c      call Calculate2DBeamKernel(Beam,Beam%PixelSize)

      return
      end subroutine
cccccccc


      end module
