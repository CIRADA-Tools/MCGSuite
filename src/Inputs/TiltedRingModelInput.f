cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines for reading in a
c       some tilted ring model parameters
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module TiltedRingInputMod

      use TiltedRingMod
      use CommonConsts

      implicit none

      contains



c           This routine reads in some tilted ring parameters in row format
      subroutine ReadTiltedRingModel(TR,HeaderFileName)
      implicit none
      Type(TiltedRingModel) TR
      character(*) HeaderFileName
      logical FileCheck

      integer unit, FormatSwitch

      print*, "Reading the tilted ring input file"
      INQUIRE(file=HeaderFileName,EXIST=FileCheck)
      if(FileCheck .eqv. .False.) then
        print*, "Tilted Ring Model file does not exist "
     &          , trim(HeaderFileName)
        stop
      endif

      unit=10
      open(unit,file=trim(HeaderFileName),status='old')

      call ReadTiltedRingModel_Column(TR,unit)

      close(unit)

      return
      end subroutine
cccccccc

cccccccc
c           This routine reads in some tilted ring parameters in column format
      subroutine ReadTiltedRingModel_Column(TR,unit)
      implicit none
      Type(TiltedRingModel),INTENT(INOUT) :: TR
      integer,INTENT(IN) :: unit
      integer i

c           The file will already have been opend in the main routine so just
c               start reading
c       First get the number of rings
      read(unit,*)
      read(unit,*) TR%nRings
c       Then the cloud mode for getting the points
      read(unit,*)
      read(unit,*) TR%cmode
c       And next get the base cloud surface density
      read(unit,*)
      read(unit,*) TR%CloudBaseSurfDens
c       Get a switch for all the different input units
      read(unit,*)
      read(unit,*) TR%UnitSwitchs(0:3)


c       Now we need to allocate the rings
      call TiltRing_Allocate(TR)

      read(unit,*)
      read(unit,*)
      do i=0, TR%nRings-1
        read(unit,*) TR%R(i)%Rmid,TR%R(i)%Rwidth
     &          ,TR%R(i)%CentPos(0),TR%R(i)%CentPos(1)
     &          ,TR%R(i)%Inclination
     &          ,TR%R(i)%PositionAngle, TR%R(i)%VSys
     &          ,TR%R(i)%VRot
     &          ,TR%R(i)%VRad, TR%R(i)%Vvert,TR%R(i)%VDisp
     &          ,TR%R(i)%dvdz, TR%R(i)%Sigma, TR%R(i)%z0
     &          ,TR%R(i)%zGradiantStart
c        print*, "Input SD", TR%R(i)%Sigma
c        print*, "Input PA",TR%R(i)%PositionAngle
c        print*, "Central position input", TR%R(i)%CentPos
      enddo


      return
      end subroutine
cccccccc



      end module
