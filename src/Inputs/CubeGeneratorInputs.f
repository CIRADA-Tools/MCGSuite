cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines for getting the
c       inputs to the cube generator code
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module CubeGeneratorInputMod

      use CubeGeneratorGlobals
      use DataCubeInputMod
      use TiltedRingInputMod
      use InputUnitConversionsMod

      implicit none

      contains

ccccccc
c           This routine is gets the top level inputs for generating a mock cube
      subroutine MockCubeGeneratorIn()
      implicit none
      character(500) TopLevelInputFile

      integer i
      real BeamPixels,BeamArea
      real arcSecPixRatio
      integer timearray(3)

      print*, "Getting the cube inputs"
      call getarg(1,TopLevelInputFile)
      if(TopLevelInputFile .eq. " ") then
        print*, "Cube Generator Infile is necessary"
        stop
      endif

      open(10, file=TopLevelInputFile,status='old')
c           Get the base name that will be used for the output folder and all file names
      read(10,*)
      read(10,'(A)') OutputBaseNames
c           Get the switch determining how much output to generate
      read(10,*)
      read(10,*) OutputFileVolumeSwitch
c           Get the cube header file
      read(10,*)
      read(10,'(A)') CubeHeaderFile
c           Get the tilted ring file
      read(10,*)
      read(10,'(A)') TiltedRingInputFile
c           Get the units for the uncertainty
      read(10,*)
      read(10,*) ModelDC%DH%UncertaintyUnitSwitch
c           Get the noise per pixel
      read(10,*)
      read(10,*) ModelDC%DH%Uncertainty
c           Get the random seed
      read(10,*)
      read(10,*) idum


      close(10)
c           If idum > 0 use the time to get the random seed
      if(idum .ge. 0) then
        print*, "using time to generate random seed"
        call itime(timeArray)
        idum=abs(timeArray(1)*idum)
        idum=-int(idum*timeArray(2))-timeArray(3)
      endif

c       Read in the data cube and header info
      call ReadTextDataCubeHeader(ModelDC,ModelBeam
     &              ,CubeHeaderFile)
c       Now get the tilted ring model
      call ReadTiltedRingModel(ModelTiltedRing,TiltedRingInputFile)
c       Include a unit conversion for surface densities

      BeamArea=Pi*ModelDC%DH%PixelSize(0)*ModelDC%DH%PixelSize(1)/4.
      arcSecPixRatio=BeamArea/1.

      BeamArea=2.*Pi*ModelBeam%BeamSigmaVector(0)**2.
      BeamPixels=BeamArea/(ModelDC%DH%PixelSize(0)
     &              *ModelDC%DH%PixelSize(1))
      ModelDC%DH%Uncertainty=ModelDC%DH%Uncertainty/1000!.*BeamPixels

      ModelDC%DH%Uncertainty=ModelDC%DH%Uncertainty*2.*sqrt(Pi
     &                      *ModelBeam%BeamSigmaVector(0)
     &                      *ModelBeam%BeamSigmaVector(1))
c     &                      /ModelDC%DH%PixelSize(0)

      call TR_UnitConversions(ModelTiltedRing
     &              ,ModelDC%DH%ChannelSize,BeamArea
     &              ,ModelDC)

      print*, "Uncertainty", BeamArea,BeamPixels,ModelDC%DH%Uncertainty
c     &                      ,ModelDC%DH%ChannelSize

      return
      end subroutine
cccccccc



      end module
