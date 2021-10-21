cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines needed to
c       fill a single ring with the particles necessary for
c       a tilted ring model
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module TiltedRingOutputsMod
      use TiltedRingMod
      use CommonConsts


      implicit none

      contains

ccccccc
      subroutine TiltedRing_ParticleOutputs(TR,fname)
      implicit none
      Type(TiltedRingModel), INTENT(IN):: TR
      character(*), INTENT(IN) :: fname

      integer i, j, nTot
      print*, "Outputting Tilted Ring Particles"

      nTot=sum(TR%R(0:TR%nRings-1)%nParticles)

      open(10,file=fname,status='replace')
      write(10,*) TR%nRings, nTot

c      i=TR%nRings-1
      do i=0,TR%nRings-1
        do j=0, TR%R(i)%nParticles-1
            write(10,*) TR%R(i)%P(j)%ProjectedPos(0:1)
     &                   ,TR%R(i)%P(j)%ProjectedVel(2)

        enddo
      enddo

      close(10)

      return
      end subroutine
cccccccc

cccccccc
c           This routine outputs the parameters from a tilted ring model
      subroutine TiltedRing_ParamOutputs(TR,fname)
      implicit none
      Type(TiltedRingModel), INTENT(IN):: TR
      character(*), INTENT(IN) :: fname
      integer i

      open(10,file=fname,status='replace')
      write(10,*) "Number of Rings"
      write(10,*) TR%nRings

      write(10,*) "X   Y   I  PA  VSys   VRot  VRad  VDisp"
     &              //"   Vvert   dvdz    Sigma    z0 "
     &              //"zGradientStart"

      do i=0, TR%nRings-1
        write(10,*) TR%R(i)%CentPos, TR%R(i)%Inclination*180./Pi
     &          ,TR%R(i)%PositionAngle*180./Pi, TR%R(i)%VSys
     &          ,TR%R(i)%VRot
     &          ,TR%R(i)%VRad, TR%R(i)%VDisp, TR%R(i)%Vvert
     &          ,TR%R(i)%dvdz, TR%R(i)%Sigma, TR%R(i)%z0
     &          ,TR%R(i)%zGradiantStart
      enddo




      close(10)


      return
      end subroutine


      end module
