cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines needed to
c       fill a single ring with the particles necessary for
c       a tilted ring model
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module CalcBeamKernelMod
      use BeamMod
      use CommonConsts


      implicit none

      contains

ccccccc
      subroutine Calculate2DBeamKernel(B,PixelSizes)
      implicit none
      Type(Beam2d), INTENT(INOUT):: B
      real,INTENT(IN):: PixelSizes(0:1)

      integer i, j, nTot
      real x,y,R2
      real xp,yp,cpa,spa
      real cellArea
c      print*, "About to calculate kernel"
c      print*, "Calculating Kernel", B%BeamSigmaVector,B%nRadialCells

      cpa=cos(-B%BeamSigmaVector(2))     !Will need to use a negative rotation for the transform
      spa=sin(-B%BeamSigmaVector(2))     ! to the beam's axes
      do i=-B%nRadialCells,B%nRadialCells
        do j=-B%nRadialCells,B%nRadialCells
c               Get the positions in x and y for the particle pixel
c            x=real(i)*PixelSizes(0)
c            y=real(j)*PixelSizes(1)
            x=real(i)
            y=real(j)
c               Now rotate these to the major axis
            xp=x*cpa-y*spa
            yp=x*spa+y*cpa
c               Get the normalized squared radius
            R2=(xp/B%BeamSigmaVector(0))**2.
     &              +(yp/B%BeamSigmaVector(1))**2.
c               Calculate the kernel value for this cell
            B%Kernel(i,j)=1./sqrt(2.*Pi*B%BeamSigmaVector(0)
     &                  *B%BeamSigmaVector(1))
     &                  *exp(-R2/2)
        enddo
      enddo
c           Renormalize so that the sum of the kernel is 1.
c      cellArea=PixelSizes(0)*PixelSizes(1)
c      B%Kernel=B%Kernel/(sum(B%Kernel)*cellArea)
c      print*, sum(B%kernel)
      B%Kernel=B%Kernel/(sum(B%Kernel))
c      print*, sum(B%kernel)

      return
      end subroutine
cccccccc


cccccccc
      subroutine CalculateComplex2DKernel(B)
      use, intrinsic :: iso_c_binding
      implicit none
      include 'fftw3.f'

      Type(Beam2d) B
      double precision,ALLOCATABLE :: PaddedKernel(:,:)
      double precision,ALLOCATABLE::WrappedPaddedKernel(:,:)
      integer i,j,k,l,CentKernel(2)
      integer*8 ArrPlan_r2c

c       Allocate  padded and padded & wrapped arrays
      ALLOCATE(PaddedKernel(B%PaddedSize(1),B%PaddedSize(2)))
      ALLOCATE(WrappedPaddedKernel(B%PaddedSize(1),B%PaddedSize(2)))

c       Set up padded kernel
      PaddedKernel=0.
      do i=1, 2*B%nRadialCells+1
        do j=1, 2*B%nRadialCells+1
            k=i-B%nRadialCells-1
            l=j-B%nRadialCells-1
            PaddedKernel(i,j)=B%Kernel(k,l)
c            print*, i,j,k,l
        enddo
      enddo

c           Wrap the padded kernel
      CentKernel=(2*B%nRadialCells+1)/2
      call MakeWrappedArray(B%PaddedSize,CentKernel
     &              ,PaddedKernel,WrappedPaddedKernel)  !In this File

c           Make the fftw plan to get the complex kernel
      call dfftw_plan_dft_r2c_2d(ArrPlan_r2c,B%PaddedSize(1)
     &                      ,B%PaddedSize(2)
     &                      ,PaddedKernel
     &                      ,B%ComplexKernel,FFTW_ESTIMATE
     &                      ,FFTW_PRESERVE_INPUT)       !FFTW3 Routine
c           Get the complex kernel
      call dfftw_execute_dft_r2c(ArrPlan_r2c, WrappedPaddedKernel
     &                          , B%ComplexKernel)      !FFTW3 Routine

c       Get rid of the fftw plan
      call dfftw_destroy_plan(ArrPlan_r2c)          !FFTW3 Routine

c       Deallocate the unneeded arrays
      DEALLOCATE(PaddedKernel)
      DEALLOCATE(WrappedPaddedKernel)

c       Note that the complex kernel has been created
      B%ComplexKernelCreated=.True.

      return
      end subroutine

cccccccccc




ccccccc
c           Wrap a double precision array such that the central value is at (1,1) in the new array
      subroutine MakeWrappedArray(SA,CV,Arr,WrappedArr)
      implicit none

      integer,INTENT(IN) :: SA(2),CV(2)     !The size of the array and the 'central values'
      double precision,INTENT(IN) :: Arr(SA(1),SA(2))
      double precision,INTENT(INOUT) :: WrappedArr(SA(1),SA(2))

      integer i,j,k,l

c           Loop over all cells
      do i=1, SA(1)
        do j=1, SA(2)
c               Get the indices of the wrapped array
            k=i-CV(1)
            l=j-CV(2)
c               Wrap in x if k<=0
            if(k .le. 0) then
                k=SA(1)+k
            endif
c               Wrap in y if l<=0
            if(l .le. 0) then
                l=SA(2)+l
            endif
c            print*, i,j,k,l
            WrappedArr(k,l)=Arr(i,j)
        enddo
      enddo


      return
      end subroutine
ccccccccccc


      end module
