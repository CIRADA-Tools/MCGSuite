cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines needed to
c       build a tilted ring model full of particles.
c          The routines assume the tilted ring has been allocated
c           and the individual ring parameters have been set.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module TiltedRingGenerationMod
      use TiltedRingMod
      use SingleRingGenerationMod


      implicit none
      contains

ccccc
      subroutine BuildTiltedRingModel(TR,pixelarea,idum)
      implicit none
      Type(TiltedRingModel), INTENT(INOUT) :: TR
      real,INTENT(IN) :: pixelarea      !The pixel area is needed to figure out how many particles are needed
      integer,INTENT(INOUT) :: idum   !A seed for the random number generation

      integer i

c      print*, "Building Tilted Ring"
c      print*, "Surf Dens Sanity", TR%R(1)%Sigma

      do i=0,TR%nRings-1
        call Ring_CalcNumParticles(TR%R(i),pixelarea        !Get the # of particles in the ring
     &              , TR%cmode,TR%CloudBaseSurfDens)
        call Ring_ParticleAllocation(TR%R(i))               !Allocate the particle array
        call Ring_ParticleGeneration(TR%R(i),idum)          !Generate the particles in the ring
      enddo

      return
      end subroutine
cccccc



      end module
