cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines needed to
c       convolve a datacube with a 2D beam
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module CubeKernelConvolutionMod
      use BeamMod
      use DataCubeMod
      use CalcBeamKernelMod
      use TwoDConvolutionMod



      implicit none

      contains

cccccccc
c           This is the main routine for convolving the cube with the beam
      subroutine CubeBeamConvolution(DC,B)
      implicit none
      Type(DataCube),INTENT(INOUT) :: DC
      Type(Beam2D), INTENT(IN) :: B

      integer i

c       First check if the complex beam kernel has been generated yet
      if(B%ComplexKernelCreated .eqv. .False.) then
c           This routine assumes that the 2D real kernel has been calculated
c           Get the complex kernel that matches the real kernel
        call CalculateComplex2DKernel(B)        !/src/ConvolveCube/CalculateBeamKernel.f
      endif
c
c           Loop through all channels
      do i=0, DC%DH%nChannels-1
c               Use the 2D convolution routine found in /src/ConvolveCube/TwoDConvolution.f
        call Convolve2D(DC%Flux(:,:,i),DC%DH%nPixels
     &          ,B%ComplexKernel
     &          ,B%ComplexSize, B%PaddedSize
     &          ,DC%Flux(:,:,i))

      enddo


      return
      end subroutine
cccccc


      end module
