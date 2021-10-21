cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the definition of the
c       general model parameter vector.  It is used in
c       the fitting routines and interfaces with the models
c       using pointer functions for generality.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module ParameterVectorMod


      implicit none

      Type ParameterVector
        integer nParams
        real,dimension(:),ALLOCATABLE :: Param
        real,dimension(:),ALLOCATABLE :: BestParam, ParamErr
        real,dimension(:),ALLOCATABLE :: ParamLowerLims
        real,dimension(:),ALLOCATABLE :: ParamUpperLims
        real CurrLike,BestLike
      end Type

      contains
cccccccc
c           This routine allocates the parameter vector
      subroutine AllocateParamVector(PV)
      implicit none
      Type(ParameterVector),INTENT(INOUT) :: PV

      ALLOCATE(PV%Param(0:PV%nParams-1))
      ALLOCATE(PV%BestParam(0:PV%nParams-1))
      ALLOCATE(PV%ParamErr(0:PV%nParams-1))
      ALLOCATE(PV%ParamLowerLims(0:PV%nParams-1))
      ALLOCATE(PV%ParamUpperLims(0:PV%nParams-1))
      PV%Param=-1.

      return
      end subroutine
c
ccccccccc

cccccccc
c           This routine deallocates the parameter vector
      subroutine DeAllocateParamVector(PV)
      implicit none
      Type(ParameterVector),INTENT(INOUT) :: PV

      DEALLOCATE(PV%Param)
      DEALLOCATE(PV%BestParam)
      DEALLOCATE(PV%ParamErr)
      DEALLOCATE(PV%ParamLowerLims)
      DEALLOCATE(PV%ParamUpperLims)

      return
      end subroutine
c
ccccccccc


      end module
