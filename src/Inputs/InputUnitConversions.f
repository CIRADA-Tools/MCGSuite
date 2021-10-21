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

      module InputUnitConversionsMod

      use DataCubeMod
      use BeamMod
      use TiltedRingMod
      use CalcBeamKernelMod

      use CommonConsts

      implicit none

      contains


ccccc
c     This routine converts angles to radians
      subroutine GeneralAngularConversion(Switch,A)
      implicit none
      integer, INTENT(IN) :: Switch
      real, INTENT(INOUT) :: A

      if(Switch .eq. 0) then        !0=degrees to radians
        A=A*Pi/180.
      endif

      return
      end subroutine
ccccccccc

cccccccc
c       This routine converts distances (angular) to arcseconds
      subroutine GeneralDistanceConversion(Switch,A)
      implicit none
      integer, INTENT(IN) :: Switch
      real, INTENT(INOUT) :: A

      if(Switch .eq. 0) then        !0=degrees to arcseconds
        A=A*3600.
      endif

      return
      end subroutine
cccccccc

ccccc
c
c       This routine converts the tilted ring
c           center coordinates to a pixel value
c
c       For this routine to work, StartVal must be the
c       pixel location (0) and have units of arcseconds
      subroutine PixelPositionConversion(Switch,A
     &                  ,PixelSize,StartVal)
      implicit none
      integer, INTENT(IN) :: Switch
      real, INTENT(IN) :: StartVal,PixelSize
      real, INTENT(INOUT) :: A
c
c      print*, "Before position Conversion", A
c     &          ,Switch, StartVal,PixelSize
      if(Switch .eq. 2) then        !2==pixels already
        A=A
      elseif(Switch .eq. 1) then    !1==arcseconds to pixels
        A=(A-StartVal)/PixelSize
      elseif(Switch .eq. 0) then    !0==degrees to pixels
        A=A*3600.
        A=(A-StartVal)/PixelSize
      endif
c      print*, "After Conversion", A

      return
      end subroutine
ccccccc


ccccc
c
c       This routine converts the tilted ring
c           radial values (radius and widths) to pixel values
c       It assumes that the radius is in arcseconds
      subroutine Radii_PixelConversion(R
     &                  ,PixelSize)
      implicit none
      real, INTENT(IN) :: PixelSize
      real, INTENT(INOUT) :: R
c
      R=R/abs(PixelSize)
      return
      end subroutine
ccccccc




cccccc
c       This routine does velocity conversions
      subroutine GeneralVelocityConversion(Switch,A)
      implicit none
      integer, INTENT(IN) :: Switch
      real, INTENT(INOUT) :: A

      if(Switch .eq. 0) then        !0=m/s to km/s
        A=A/1000.
      endif

      return
      end subroutine
ccccccc

ccccccc
c     This routine does brightness conversion
      subroutine GeneralBrightnessConversion(Switch,A,ChannelSize
     &              ,BeamArea,PixelSize)
      implicit none
      integer, INTENT(IN) :: Switch
      real, INTENT(INOUT) :: A
      real, INTENT(IN) :: ChannelSize,BeamArea,PixelSize

c      print*, "SD Convert", Switch,A,ChannelSize,BeamArea
c     &                  ,PixelSize
c       The surface brightness needs to be converted to Jy pixel^2
      if(Switch .eq. 0) then        !0=Jy km/s arcsec^-2 to Jy arcsec^-2
c        print*, "SD Convert Factor",abs(PixelSize)**2./abs(ChannelSize)
        A=A/abs(ChannelSize)
        A=A*abs(PixelSize)**2.           !Then Jy arcsec^2 to Jy pixel
      elseif(Switch .eq. 1) then    !1=mJy/beam to Jy/Beam
        A=A/1000.
      elseif(Switch .eq. 2) then    !2=Jy arcsec^-2 to Jy/Beam
        A=A*BeamArea
      endif
      return
      end subroutine
ccccccc

cccccc
c       This routine does the angular conversions for tilted ring parameters
      subroutine TR_UnitConversions(TR,ChannelSize,BeamArea
     &                  ,DC)
      implicit none
      Type(TiltedRingModel),INTENT(INOUT) :: TR
      Type(DataCube),INTENT(IN) :: DC
      real, INTENT(IN) :: ChannelSize,BeamArea
      integer i
c       Loop through all rings
      do i=0, TR%nRings-1
c       Do the central position conversions first  --> The final positions should be in pixels
        call PixelPositionConversion(TR%UnitSwitchs(0)
     &          ,TR%R(i)%CentPos(0)
     &          ,DC%DH%PixelSize(0),DC%DH%Start(0))
        call PixelPositionConversion(TR%UnitSwitchs(0)
     &          ,TR%R(i)%CentPos(1)
     &          ,DC%DH%PixelSize(1),DC%DH%Start(1))
c        print*, "Cent Pos",TR%R(i)%CentPos(0:1)
c       Now do the angular conversions
c       Do the geometric conversions first
        call GeneralAngularConversion(TR%UnitSwitchs(1)
     &          ,TR%R(i)%Inclination)
        call GeneralAngularConversion(TR%UnitSwitchs(1)
     &          ,TR%R(i)%PositionAngle)
c           Once the position angle is in radians, add a 90^degree rotation to get the standard astronomy
c               definition for position angle (North-South with approaching velocity in the North)
c                   Note that in outputs we'll need to undo this rotation to record the position angle
        TR%R(i)%PositionAngle=TR%R(i)%PositionAngle+Pi/2.
c
c       Next do all the velocities
        call GeneralVelocityConversion(TR%UnitSwitchs(2)
     &          ,TR%R(i)%VSys)
        call GeneralVelocityConversion(TR%UnitSwitchs(2)
     &          ,TR%R(i)%VRot)
        call GeneralVelocityConversion(TR%UnitSwitchs(2)
     &          ,TR%R(i)%VRad)
        call GeneralVelocityConversion(TR%UnitSwitchs(2)
     &          ,TR%R(i)%Vvert)
        call GeneralVelocityConversion(TR%UnitSwitchs(2)
     &          ,TR%R(i)%VDisp)
        call GeneralVelocityConversion(TR%UnitSwitchs(2)
     &          ,TR%R(i)%dvdz)
c       Do the surface brightness conversion
        call GeneralBrightnessConversion(TR%UnitSwitchs(3)
     &          ,TR%R(i)%Sigma,ChannelSize,BeamArea
     &          ,DC%DH%PixelSize(0))
c       Finally convert the radii to pixels
        call Radii_PixelConversion(TR%R(i)%Rmid
     &              ,DC%DH%PixelSize(0))
        call Radii_PixelConversion(TR%R(i)%Rwidth
     &              ,DC%DH%PixelSize(0))
c       Also make sure to convert z0 and z0Grad to pixels
        call Radii_PixelConversion(TR%R(i)%z0
     &              ,DC%DH%PixelSize(0))
        call Radii_PixelConversion(TR%R(i)%zGradiantStart
     &              ,DC%DH%PixelSize(0))
      enddo


      return
      end subroutine
ccccccccc


ccccc
c       This routine does the angular conversions for the data cube parameters
      subroutine TextDCHeader_AngularConversions(DC,Beam)
      implicit none
      Type(DataCube),INTENT(INOUT) :: DC
      Type(Beam2D), INTENT(INOUT) :: Beam

c       Do the data cube unit conversions
c           Start with the pixel sizes
      call GeneralDistanceConversion(DC%DH%DimensionUnitSwitch(0)
     &          ,DC%DH%PixelSize(0))    !This File
      call GeneralDistanceConversion(DC%DH%DimensionUnitSwitch(0)
     &          ,DC%DH%PixelSize(1))    !This File
c           Next the channel size
      call GeneralVelocityConversion(DC%DH%DimensionUnitSwitch(2)
     &          ,DC%DH%ChannelSize)     !This File
c       Now do up the reference values
      call GeneralDistanceConversion(DC%DH%DimensionUnitSwitch(1)
     &          ,DC%DH%RefVal(0))       !This Fil
      call GeneralDistanceConversion(DC%DH%DimensionUnitSwitch(1)
     &          ,DC%DH%RefVal(1))       !This File
      call GeneralVelocityConversion(DC%DH%DimensionUnitSwitch(3)
     &          ,DC%DH%RefVal(2))       !This File
c      print*, "New Ref Values", DC%DH%RefVal(0:2)


c      Next do the Beam unit conversions
c           Want to have the beam in pixel dimensions

      call Radii_PixelConversion(Beam%BeamMajorAxis
     &              ,DC%DH%PixelSize(0))
      call Radii_PixelConversion(Beam%BeamMinorAxis
     &              ,DC%DH%PixelSize(0))

c      call GeneralAngularConversion(Beam%BeamUnitsSwitch(0)
c     &          ,Beam%BeamMajorAxis)     !This File
c      call GeneralDistanceConversion(Beam%BeamUnitsSwitch(0)
c     &          ,Beam%BeamMinorAxis)    !This File
      call GeneralAngularConversion(Beam%BeamUnitsSwitch(1)
     &          ,Beam%BeamPositionAngle)    !This File

      return
      end subroutine
ccccccc

      end module
