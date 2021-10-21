cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the general procedures for modeling
c       a galaxy once the pre-analysis is completed.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module GalaxyFitMod
      use ParameterVectorMod
      use DataCubeMod
      use SoFiACatalogueMod
      use BeamMod
      use CalcBeamKernelMod
      use DownhillSimplexMod
      use FullModelComparisonMod
      implicit none

      PROCEDURE(GeneralFitInterface),POINTER :: GalaxyFit =>null()

      ABSTRACT INTERFACE
        subroutine GeneralFitInterface(CatItem)
            import :: CatalogueItem
            Type(CatalogueItem),INTENT(IN) :: CatItem
        END subroutine GeneralFitInterface
      END INTERFACE

      contains
cccccccc
c
      subroutine GalaxyFit_Simple(CatItem)
      use PipelineGlobals
      implicit none
      Type(CatalogueItem),INTENT(IN) :: CatItem

      real,ALLOCATABLE :: paramGuesses(:,:),chiArray(:)
      real chi2
      integer i,iter

c      print*, "Fitting Galaxy",PID
c
c       Set up the beam
c
c      call Allocate_Beam2D(ObservedBeam,ObservedDC%DH%nPixels)
      call Calculate2DBeamKernel(ObservedBeam,ObservedDC%DH%PixelSize)

c       Allocate a model datacube with the same dimensions as the observed cube
      ModelDC%DH=ObservedDC%DH
      call AllocateDataCube(ModelDC)
c       Copy the initial parameter guess into the model parameter vector
      PVModel%nParams=PVIni%nParams
      call AllocateParamVector(PVModel)
      PVModel%Param=PVIni%Param
      PVModel%ParamLowerLims=PVIni%ParamLowerLims
      PVModel%ParamUpperLims=PVIni%ParamUpperLims
c           Set up the array of parameter guesses and chi^2 values needed
      ALLOCATE(paramGuesses(PVModel%nParams+1,
     &          PVModel%nParams))
      ALLOCATE(chiArray(PVModel%nParams+1))
c           Create an array of parameter guesses
      IniGuessWidth=0.1
      call MakeParamGuessArray(PVModel,ParamGuesses
     &          ,PVIni%nParams,idum,IniGuessWidth)
c       Get the initial gueses chi^2 values
      chiArray=0.
c      do i=1,1
      do i=1,PVModel%nParams+1
        PVModel%Param(0:PVModel%nParams-1)=
     &                  ParamGuesses(i,1:PVModel%nParams)
        call TiltedRingModelComparison(PVModel%Param,chi2)
        chiArray(i)=chi2
c        print*, PID,i,paramGuesses(i,:),chi2
c        print*, "chitest", chiArray
      enddo
c      print*, "hmmm", PID,chiArray

      ftol=2.5e-2

c       Run the main fitting routine
      call amoeba(paramGuesses,chiArray
     &                  ,PVModel%nParams+1
     &                  ,PVModel%nParams
     &                  ,PVModel%nParams,ftol
     &                  ,TiltedRingModelComparison,iter,PID)

      print*, "Best Fit", PID
     &          ,paramGuesses(1,1:PVModel%nParams),chiArray(1)

      PVModel%Param(0:PVModel%nParams-1)=
     &              ParamGuesses(1,1:PVModel%nParams)
c      print*, PVModel%Param


c       Deallocate the model vector at the end
c      call DeAllocateParamVector(PVModel)
      DEALLOCATE(chiArray)
      DEALLOCATE(paramGuesses)
      return
      end subroutine
ccccccccc



ccccccc
c


ccccccc
      subroutine MakeParamGuessArray(PredictedPV,ParamGuesses
     &                      ,ndim,idum,lambda)

      use BasicRanNumGen
      implicit none

      integer idum
      integer ndim
      Type(ParameterVector) PredictedPV
      real ParamGuesses(ndim+1,ndim)
      real lambda,lambdaPar

      integer i,j,k


c      print*, "Param Guess Array Creation"
      do k=1,ndim+1
        if(k .eq. 1) then
            ParamGuesses(k,1:ndim)=PredictedPV%Param(0:ndim-1)
        else
            do i=1,ndim
                j=i-1
100             lambdaPar=lambda*(PredictedPV%ParamUpperLims(j)
     &                      -PredictedPV%ParamLowerLims(j))
                lambdaPar=(2*ran2(idum)-1.)*lambdaPar
                ParamGuesses(k,i)=PredictedPV%Param(j)+lambdaPar
c            print*, i,ParamGuesses(k,i),PredictedPV%ParamLowerLims(j)
c     &                  ,PredictedPV%ParamUpperLims(j)
                if(ParamGuesses(k,i) .lt.
     &                   PredictedPV%ParamLowerLims(j)) goto 100
                if(ParamGuesses(k,i) .gt.
     &                   PredictedPV%ParamUpperLims(j)) goto 100
            enddo
        endif
c        print*, k,ParamGuesses(k,:)
      enddo

      return
      end subroutine
cccccc

      end module
