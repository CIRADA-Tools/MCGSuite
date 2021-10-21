cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the definitions of the Beam object
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module BeamMod
      implicit none

      Type Beam2D
        real BeamMajorAxis,BeamMinorAxis,BeamPositionAngle
        real BeamFWHM,SigmaLengths
        real BeamSigmaVector(0:2)
        real PixelSize(0:1)
        integer nRadialCells
        real,dimension(:,:),ALLOCATABLE :: Kernel

        integer VelocitySmoothSwitch
        real VelocitySmoothSigma

        integer BeamUnitsSwitch(0:1)


        integer PaddedSize(2),ComplexSize(2)
        double complex,dimension(:,:),ALLOCATABLE :: ComplexKernel
        logical ComplexKernelCreated
      end Type




      contains
ccccccc
      subroutine Allocate_Beam2D(B,nCubePixels)
      implicit none
      Type(Beam2D),INTENT(INOUT) :: B
      integer nCubePixels(2)
c       If we are using the FWHM, set the appropriate values (this is indicated by the BeamFWHM >0)
c           The 2.355 factor is to convert from FWHM values to sigma values
      if(B%BeamFWHM .gt. 0.) then
        B%BeamSigmaVector(0)=B%BeamFWHM/2.355
        B%BeamSigmaVector(1)=B%BeamFWHM/2.355
        B%BeamSigmaVector(2)=0.
      else
c           If BeamFWHM <0 then the axis values should be set prior to the allocation
        B%BeamSigmaVector(0)=B%BeamMajorAxis/2.355
        B%BeamSigmaVector(1)=B%BeamMinorAxis/2.355
        B%BeamSigmaVector(2)=B%BeamPositionAngle
      endif
c      print*,"Beam Allocation", B%BeamSigmaVector, B%PixelSize
c     &          ,B%SigmaLengths

c           Calculate the number of cells in the largest direction
c      B%nRadialCells=((B%BeamSigmaVector(0)/minval(abs(B%PixelSize)))
c     &                  *B%SigmaLengths) +1
      B%nRadialCells=int(B%BeamSigmaVector(0)*B%SigmaLengths)+1

c      print*, "Cells in one radial direction", B%nRadialCells
c           Allocate a square kernel to account for all the cells
      ALLOCATE(B%Kernel(-B%nRadialCells:B%nRadialCells
     &                  ,-B%nRadialCells:B%nRadialCells))

      B%ComplexKernelCreated=.False.
      B%PaddedSize=2*B%nRadialCells+1+nCubePixels
      B%ComplexSize(1)=B%PaddedSize(1)/2+1
      B%ComplexSize(2)=B%PaddedSize(2)
      ALLOCATE(B%ComplexKernel(B%ComplexSize(1),B%ComplexSize(2)))  !Note the array is centered on 1 to make
c                                                                       sure it works with the fftw library
c      print*, "Beam Size", B%nRadialCells,B%PixelSize
c     &          , B%PaddedSize, B%ComplexSize

      return
      end subroutine
ccccccc

cccccccc
      subroutine DeAllocate_Beam2D(B)
      implicit none
      Type(Beam2D) B
c
      DEALLOCATE(B%Kernel)
      DEALLOCATE(B%ComplexKernel)
      B%ComplexKernelCreated=.False.
      return
      end subroutine
ccccccccc

      end module
