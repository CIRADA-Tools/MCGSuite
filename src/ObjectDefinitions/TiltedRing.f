cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the definition of a tilted
c       ring model.  It requires the Particle module found
c       in /src/ObjectDefinitions/Particle.f
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module TiltedRingMod
      use ParticleMod

      implicit none

      Type Ring
        integer nParticles
        real Rmid, Rwidth
        real CentPos(0:1), Inclination, PositionAngle
        real VSys, VRot,VRad, VDisp, Vvert, dvdz
        real dVRot_dR(0:1), dVRad_dR(0:1), dVDisp_dR(0:1)
        real Sigma,dSigma_dR(0:1)

        real z0,zGradiantStart
        Type(Particle), ALLOCATABLE :: P(:)
      end Type

      Type TiltedRingModel
        integer nRings,cmode
        real CloudBaseSurfDens
        Type(Ring),ALLOCATABLE :: R(:)
        integer UnitSwitchs(0:3)  !0->CentPos, 1->Inclination;PositionAngle, 2->Velocity Quantities, 3->Sigma
      end Type


      Type TiltedRingFittingOptions
        integer nRings,nFittedParamsTotal
        integer nFittedRadialParams,nFixedRadialParams
        integer nFittedConstantParams,nFixedConstantParams
        logical ConstParams(0:12)
        logical FixedParams(0:12)
        real ParamLowerLims(0:12), ParamUpperLims(0:12)
        Type(Ring),ALLOCATABLE,dimension(:) :: RadialProfiles

c       0=  Xcent
c       1=  Ycent
c       2=  Inclination
c       3=  Position Angle
c       4=  VSys
c       5=  VRot
c       6=  VRad
c       7=  VDisp
c       8=  Vvert
c       9=  dvdz
c       10= Sigma
c       11= z0
c       12= zGradiantStart

      end Type




      contains


cccccc
      subroutine Ring_ParticleAllocation(R)
      implicit none
      Type(Ring), INTENT(INOUT) :: R
c
      ALLOCATE(R%P(0:R%nParticles-1))       !Keep indexing consistent with C
      return
      end subroutine
cccccccc

cccccc
      subroutine Ring_ParticleDeAllocation(R)
      implicit none
      Type(Ring) R
c
      DEALLOCATE(R%P)
      return
      end subroutine
cccccccc


ccccccc
      subroutine TiltRing_Allocate(TR)
      implicit none
      Type(TiltedRingModel), INTENT(INOUT) :: TR
      ALLOCATE(TR%R(0:TR%nRings-1))
      return
      end subroutine
cccccccc


ccccccc
      subroutine TiltRing_DeAllocate(TR)
      implicit none
      Type(TiltedRingModel), INTENT(INOUT) :: TR
      integer i

      do i=0, TR%nRings-1
        call Ring_ParticleDeAllocation(TR%R(i))
      enddo
      DEALLOCATE(TR%R)
      return
      end subroutine
cccccccc

ccccccc
      subroutine TiltRing_DeAllocateStruct(TR)
      implicit none
      Type(TiltedRingModel), INTENT(INOUT) :: TR
      integer i

      DEALLOCATE(TR%R)
      return
      end subroutine
cccccccc


ccccccc
      subroutine TiltRingFittingOptions_Allocate(TRFO)
      implicit none
      Type(TiltedRingFittingOptions), INTENT(INOUT) :: TRFO
      integer i
c           Allocate the rings
      ALLOCATE(TRFO%RadialProfiles(0:TRFO%nRings-1))
cc           Set all logicals to False as a default (indicating that all should be fit)
c      TRFO%ConstParams=.False.
c      TRFO%FixedParams=.False.
c           Initialize all parameters to 0
      do i=0,TRFO%nRings-1
        TRFO%RadialProfiles(i)%CentPos=0.
        TRFO%RadialProfiles(i)%Inclination=0.
        TRFO%RadialProfiles(i)%PositionAngle=0.
        TRFO%RadialProfiles(i)%VSys=0.
        TRFO%RadialProfiles(i)%VRot=0.
        TRFO%RadialProfiles(i)%VRad=0.
        TRFO%RadialProfiles(i)%VDisp=0.
        TRFO%RadialProfiles(i)%Vvert=0.
        TRFO%RadialProfiles(i)%dvdz=0.
        TRFO%RadialProfiles(i)%Sigma=0.
        TRFO%RadialProfiles(i)%z0=0.
        TRFO%RadialProfiles(i)%zGradiantStart=0.


      enddo

      return
      end subroutine
cccccccc


ccccccc
      subroutine TiltRingFittingOptions_DeAllocate(TRFO)
      implicit none
      Type(TiltedRingFittingOptions), INTENT(INOUT) :: TRFO
      DEALLOCATE(TRFO%RadialProfiles)
      return
      end subroutine
cccccccc


cccccccccc
c           This routine figures out the indexing of parameters and
c           gets the number of free parameters for some set of
c           logical flags and rings
      subroutine LogicalTiltedRingIndexing(TRFO)
      implicit none
      Type(TiltedRingFittingOptions), INTENT(INOUT) :: TRFO
      integer i

      TRFO%nFittedConstantParams=0
      TRFO%nFittedRadialParams=0
      do i=0,12
        if(TRFO%ConstParams(i)) then
            if(TRFO%FixedParams(i)) then
                TRFO%nFixedConstantParams=TRFO%nFixedConstantParams+1
            else
                TRFO%nFittedConstantParams=TRFO%nFittedConstantParams+1
            endif
        else
            if(TRFO%FixedParams(i)) then
                TRFO%nFixedRadialParams=TRFO%nFixedRadialParams+1
            else
                TRFO%nFittedRadialParams=TRFO%nFittedRadialParams+1
            endif
        endif
      enddo

c      print*, "Truth Check", TRFO%ConstParams
c      print*, "number check",TRFO%nFittedConstantParams
c     &              ,TRFO%nFittedRadialParams,TRFO%nRings

      TRFO%nFittedParamsTotal=
     &                  TRFO%nFittedConstantParams
     &                  +TRFO%nFittedRadialParams
     &                  *TRFO%nRings

      return
      end subroutine
cccccccc



      end module
