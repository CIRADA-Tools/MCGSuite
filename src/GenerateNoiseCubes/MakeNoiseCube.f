cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains a set of functions for generating
c       cubes that are full of noise.  These are for both generating
c       realistic mock data and bootstrapping samples
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module NoiseCubeMod
      use DataCubeMod
      use BasicRanNumGen
      use CubeKernelConvolutionMod

      contains
cccccccc
c       This routine makes a source cube full of point source noise
      subroutine MakePointSourceNoiseCube(DC,idum)
      implicit none
      Type(DataCube),INTENT(INOUT) :: DC
      integer, INTENT(INOUT) :: idum
      integer i,j,k

      do i=0, DC%DH%nPixels(0)-1
        do j=0, DC%DH%nPixels(1)-1
            do k=0,DC%DH%nChannels-1
c                   Draw random values from a Gaussian for the noise
                DC%Flux(i,j,k)=DC%DH%Uncertainty*gasdev(idum)
            enddo
        enddo
      enddo
      return
      end subroutine
cccccc


ccccccc
c       This routine makes a point source cube and convolves it with a beam
      subroutine MakeConvolvedNoiseCube(DC,Beam,idum)
      implicit none
      Type(DataCube),INTENT(INOUT) :: DC
      Type(Beam2D),INTENT(INOUT) :: Beam
      integer, INTENT(INOUT) :: idum
c
c       Begin by making the point source cube
      call MakePointSourceNoiseCube(DC,idum)    !In this file.

c       Now convolve that with the beam
      call CubeBeamConvolution(DC,Beam)     !/src/ConvolveCube/CubeKernelConvolution.f


      return
      end subroutine
ccccccc

      end module
