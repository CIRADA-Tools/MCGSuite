cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines needed to
c       convolve some 2D array with a kernel
c
c           The routines assume that the complex kernel has
c           already been calculated
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module TwoDConvolutionMod




      implicit none

      contains

ccccccc
      subroutine Convolve2D(Arr,SizeArray,ComplexKernel,SC,SizePad
     &          ,ConvolvedArray)

      use, intrinsic :: iso_c_binding
      implicit none
      include 'fftw3.f'

      integer,INTENT(IN) :: SizeArray(2)
      integer,INTENT(IN) :: SizePad(2)

      real,INTENT(IN) :: Arr(SizeArray(1),SizeArray(2))
      real,INTENT(OUT) :: ConvolvedArray(SizeArray(1),SizeArray(2))

      double precision,ALLOCATABLE :: PaddedArray(:,:),PaddedKernel(:,:)
      double precision,ALLOCATABLE :: RealConvolve(:,:)
      double precision, ALLOCATABLE :: WrappedPaddedKernel(:,:)

      integer i,j,k,l

      integer,INTENT(IN) :: SC(2)       !The size of the complex arrays
c      integer SizeComplex(2)
      double complex,ALLOCATABLE :: ComplexArr(:,:)
      double complex, INTENT(IN) :: ComplexKernel(SC(1),SC(2))
      double complex,ALLOCATABLE :: ComplexConvolve(:,:)

      integer*8 ArrPlan_r2c,ArrPlan_c2r
ccccccccc
c           Allocate Arrays
      call AllocateConvolutionArrays(SizeArray,Arr
     &              ,SizePad,SC,PaddedArray
     &              ,RealConvolve
     &              ,ComplexArr
     &              ,ComplexConvolve)           !In this file

      call SetupPaddedArrays(SizeArray,SizePad
     &                      ,Arr
     &                      ,PaddedArray)       !In this file

c       Make the real-to-complex fftw plan
      call dfftw_plan_dft_r2c_2d(ArrPlan_r2c,SizePad(1),SizePad(2)
     &                      ,PaddedArray
     &                      ,ComplexArr,FFTW_ESTIMATE
     &                      ,FFTW_PRESERVE_INPUT)                     !FFTW3 Routine
c       Make the complex-to-real fftw plan
      call dfftw_plan_dft_c2r_2d(ArrPlan_c2r,SizePad(1),SizePad(2)
     &              ,ComplexArr,PaddedArray
     &              ,FFTW_ESTIMATE,FFTW_PRESERVE_INPUT)                !FFTW3 Routine

c           Do the fft transform of the padded image
      call dfftw_execute_dft_r2c(ArrPlan_r2c, PaddedArray, ComplexArr)  !FFTW3 Routine

c           Convolve the transformed arrays
      do i=1, SizePad(1)/2+1
        do j=1, SizePad(2)
            ComplexConvolve(i,j)=ComplexArr(i,j)*ComplexKernel(i,j)
        enddo
      enddo

c           Do the backwards fft transform fo the convolved array

      call dfftw_execute_dft_c2r(ArrPlan_c2r, ComplexConvolve
     &                  ,RealConvolve)               !FFTW3 Routine
c           Normalize the real array
      RealConvolve=RealConvolve/SizePad(1)/SizePad(2)

c           Set the convolved array (of the correct size) to the real-convolve array
      do i=1, SizeArray(1)
        do j=1, SizeArray(2)
            ConvolvedArray(i,j)=RealConvolve(i,j)
        enddo
      enddo
c           Destroy the fftw plans
      call dfftw_destroy_plan(ArrPlan_r2c)      !FFTW3 Routine
      call dfftw_destroy_plan(ArrPlan_c2r)      !FFTW3 Routine

c       Free up the padded array space
      DEALLOCATE(PaddedArray)
      DEALLOCATE(ComplexArr)
      DEALLOCATE(ComplexConvolve)
      DEALLOCATE(RealConvolve)

      return
      end subroutine
cccccccc


cccccc
c           This function allocates all the different arrays needed
      subroutine AllocateConvolutionArrays(SizeArray
     &              ,Arr
     &              ,SizePad,SizeComplex,PaddedArray
     &              ,RealConvolve
     &              ,ComplexArr
     &              ,ComplexConvolve)
      implicit none
      integer,INTENT(IN) :: SizeArray(2), SizeComplex(2)
      integer,INTENT(IN) :: SizePad(2)

      real,INTENT(IN) :: Arr(SizeArray(1),SizeArray(2))

      double precision,ALLOCATABLE,INTENT(OUT) :: PaddedArray(:,:)
      double precision,ALLOCATABLE,INTENT(OUT) :: RealConvolve(:,:)

      double complex,ALLOCATABLE,INTENT(OUT) :: ComplexArr(:,:)
      double complex,ALLOCATABLE,INTENT(OUT) :: ComplexConvolve(:,:)


c           Allocate the real/double precision arrays
      ALLOCATE(PaddedArray(SizePad(1),SizePad(2)))
      ALLOCATE(RealConvolve(SizePad(1),SizePad(2)))
c           Allocate the complex arrays
      ALLOCATE(ComplexArr(SizeComplex(1),SizeComplex(2)))
      ALLOCATE(ComplexConvolve(SizeComplex(1),SizeComplex(2)))

      return
      end subroutine
cccccccc

ccccccccc
c           This routine initializes the padded arrays
      subroutine SetupPaddedArrays(SizeArray,SP
     &                      ,Arr
     &                      ,PaddedArray)
      implicit none

      integer, INTENT(IN) :: SizeArray(2), SP(2)
      real,INTENT(IN) :: Arr(SizeArray(1),SizeArray(2))
      double precision,INTENT(INOUT)::PaddedArray(SP(1),SP(2))

      integer i, j
c          Make the padded array
      PaddedArray=0.
      do i=1, SizeArray(1)
        do j=1, SizeArray(2)
            PaddedArray(i,j)=Arr(i,j)
        enddo
      enddo

      return
      end subroutine
cccccccc



      end module
