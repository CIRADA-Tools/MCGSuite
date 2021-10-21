cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines needed to
c       fill a single ring with the particles necessary for
c       a tilted ring model
c
c       All the routines assume that some basic parameters
c       of the tilted ring model have been set:
c           nRings, cmode, CloudBaseSurfDens, Rmid, Rwidth
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module ParameterVectorToTiltedRingMod
      use TiltedRingMod
      use ParameterVectorMod

      implicit none

      Procedure(P_To_TR_Interface),POINTER:: ParamToTiltedRing=>null()

      ABSTRACT INTERFACE
        subroutine P_To_TR_Interface(PV,TR,TRFO)
            import :: ParameterVector
            import :: TiltedRingModel
            import :: TiltedRingFittingOptions
            implicit none
            Type(ParameterVector), INTENT(IN) :: PV
            Type(TiltedRingModel), target, INTENT(INOUT) :: TR
            Type(TiltedRingFittingOptions), target,INTENT(IN) :: TRFO
        end subroutine P_To_TR_Interface
      END INTERFACE


      contains

ccccccccccc
c
      subroutine SimpleFlatDiskParamToTiltedRing(PV,TR,TRFO)
      implicit none
      Type(ParameterVector), INTENT(IN) :: PV
      Type(TiltedRingModel), target,INTENT(INOUT) :: TR
      Type(TiltedRingFittingOptions), target,INTENT(IN) :: TRFO
      real Inclination, PositionAngle
      real ConstVSys,ConstVDisp
      real ConstXCent,ConstYCent
      integer i,j

      real derivLow,derivHigh

      print*, "Converting Param vector to flat tilted ring"
c           Set the constant geometry parameters
      Inclination=PV%Param(0)
      PositionAngle=PV%Param(1)
      ConstXCent=PV%Param(2)
      ConstYCent=PV%Param(3)
c           Set the constant velocity parameters
      ConstVSys=PV%Param(4)
      ConstVDisp=PV%Param(5)


      do i=0, TR%nRings-1
        j=6+i*2
c               Set the free ring parameters
        TR%R(i)%VRot=PV%Param(j)
        TR%R(i)%Sigma=PV%Param(j+1)
c               Zero out the simple parameters
        TR%R(i)%VRad=0.
        TR%R(i)%Vvert=0.
        TR%R(i)%dvdz=0.
c               Set the sech^2 scale height to a constant
        TR%R(i)%z0=0.1
        TR%R(i)%zGradiantStart=5.*TR%R(j)%z0
c               Set the parameters that are constant across all rings
        TR%R(i)%VSys=ConstVSys
        TR%R(i)%VDisp=ConstVDisp
        TR%R(i)%CentPos(0)=ConstXCent
        TR%R(i)%CentPos(1)=ConstYCent
        TR%R(i)%Inclination=Inclination
        TR%R(i)%PositionAngle=PositionAngle
      enddo
c       Calculate the various important derivatives on both sides of each point
      call TiltedRingDerivativeCalc(TR)


      return
      end subroutine
ccccccccc

cccccccc
c       This routine does a generalized conversion of a parameter vector to a tilted ring model
      subroutine GeneralizedParamVectorToTiltedRing(PV,TR,TRFO)
      implicit none
      Type(ParameterVector), INTENT(IN) :: PV
      Type(TiltedRingModel), target,INTENT(INOUT) :: TR
      Type(TiltedRingFittingOptions), target, INTENT(IN) :: TRFO
      integer i, CurrParam
      real,pointer,dimension(:) :: TempVector,VariableVector

c      print*, "in generalize param vector"

      CurrParam=0
      ALLOCATE(TempVector(0:TR%nRings-1))
      ALLOCATE(VariableVector(0:TR%nRings-1))

c      print*, "Param Vec", PV%Param(0:PV%nParams-1), PV%nParams

c           Loop through all the specifc parameters
      do i=0,12
c           Associate the temporary pointers with the correct array
        if(i .eq.0) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%CentPos(0)
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%CentPos(0)
        elseif(i .eq.1) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%CentPos(1)
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%CentPos(1)
        elseif(i .eq.2) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%Inclination
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%Inclination
        elseif(i .eq.3) then
            VariableVector=>
     &            TRFO%RadialProfiles(0:TRFO%nRings-1)%PositionAngle
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%PositionAngle
        elseif(i .eq.4) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%VSys
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%VSys
        elseif(i .eq.5) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%VRot
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%VRot
        elseif(i .eq.6) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%VRad
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%VRad
        elseif(i .eq.7) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%VDisp
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%VDisp
        elseif(i .eq.8) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%Vvert
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%Vvert
        elseif(i .eq.9) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%dvdz
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%dvdz
        elseif(i .eq.10) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%Sigma
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%Sigma
        elseif(i .eq.11) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%z0
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%z0
        elseif(i .eq.12) then
            VariableVector=>
     &            TRFO%RadialProfiles(0:TRFO%nRings-1)%zGradiantStart
            TempVector=>
     &                  TR%R(0:TRFO%nRings-1)%zGradiantStart
c
        endif
c               Pass all the vectors to the vector assignment
        call SetSpecificVector(PV%nParams,TR%nRings,TRFO%ConstParams(i)
     &          ,TRFO%FixedParams(i),CurrParam
     &          ,PV%Param,VariableVector,TempVector)

c        print*, "Fit Check",TRFO%ConstParams(i),TRFO%FixedParams(i)
c        print*, "vari", VariableVector
c        print*, i, TempVector

c        DEALLOCATE(VariableVector)
c        DEALLOCATE(TempVector)
      enddo
c       Calculate the various important derivatives on both sides of each point
      call TiltedRingDerivativeCalc(TR)



      return
      end subroutine
ccccccc

cccccc
c           This routine sets a specifc vector based on the fitting options
      subroutine SetSpecificVector(nParam,nRings,Const,Fixed,CurrParam
     &                  ,Param,AlternateSourceVector
     &                  ,TargVec)
      implicit none
      integer,INTENT(IN) :: nRings,nParam
      logical,INTENT(IN) :: Const,Fixed
      real,INTENT(IN) :: Param(0:nParam-1)
      real, INTENT(IN) :: AlternateSourceVector(0:nRings-1)
      integer,INTENT(INOUT) :: CurrParam
      real,INTENT(INOUT) :: TargVec(0:nRings-1)

      integer i

c      print*, "Setting Vec", CurrParam, Param(Currparam),Const,Fixed
      if(Const) then
        if(Fixed) then
            TargVec(0:nRings-1)=AlternateSourceVector(0)
        else
            TargVec(0:nRings-1)=Param(CurrParam)
            CurrParam=CurrParam+1
        endif
      else
        if(Fixed) then
            TargVec(0:nRings-1)=AlternateSourceVector(0:nRings-1)
        else
            TargVec(0:nRings-1)=Param(CurrParam:CurrParam+nRings-1)
            CurrParam=CurrParam+nRings
        endif
      endif
c      print*,TargVec(0:nRings-1)


      return
      end subroutine
cccccccc

cccccccc
c       This subroutine calculate all the derivatives used in a tilted ring model
      subroutine TiltedRingDerivativeCalc(TR)
      implicit none
      Type(TiltedRingModel), INTENT(INOUT) :: TR
      integer i

      do i=0, TR%nRings-1
        if(i .eq. 0) then
            TR%R(i)%dVRot_dR(0)=(TR%R(i)%VRot-0.)/(TR%R(i)%Rwidth/2)
            TR%R(i)%dVRad_dR(0)=(TR%R(i)%VRad-0.)/(TR%R(i)%Rwidth/2)
c            TR%R(i)%dVDisp_dR(0)=(TR%R(i)%VDisp-0.)/(TR%R(i)%Rwidth/2)
            TR%R(i)%dVDisp_dR(0)=0.
            TR%R(i)%dSigma_dR(0)=0.
            TR%R(i)%dVRot_dR(1)=(TR%R(i+1)%VRot-TR%R(i)%VRot)
     &                      /TR%R(i)%Rwidth
            TR%R(i)%dVRad_dR(1)=(TR%R(i+1)%VRad-TR%R(i)%VRad)
     &                      /TR%R(i)%Rwidth
            TR%R(i)%dVDisp_dR(1)=(TR%R(i+1)%VDisp-TR%R(i)%VDisp)
     &                      /TR%R(i)%Rwidth
            TR%R(i)%dSigma_dR(1)=(TR%R(i+1)%Sigma-TR%R(i)%Sigma)
     &                      /TR%R(i)%Rwidth
        elseif(i .eq. TR%nRings-1) then
            TR%R(i)%dVRot_dR(0)=TR%R(i-1)%dVRot_dR(1)
            TR%R(i)%dVRad_dR(0)=TR%R(i-1)%dVRad_dR(1)
            TR%R(i)%dVDisp_dR(0)=TR%R(i-1)%dVDisp_dR(1)
            TR%R(i)%dSigma_dR(0)=TR%R(i-1)%dSigma_dR(1)

            TR%R(i)%dVRot_dR(1)=TR%R(i-1)%dVRot_dR(1)
            TR%R(i)%dVRad_dR(1)=TR%R(i-1)%dVRad_dR(1)
            TR%R(i)%dVDisp_dR(1)=TR%R(i-1)%dVDisp_dR(1)
            TR%R(i)%dSigma_dR(1)=TR%R(i-1)%dSigma_dR(1)
        else
            TR%R(i)%dVRot_dR(0)=TR%R(i-1)%dVRot_dR(1)
            TR%R(i)%dVRad_dR(0)=TR%R(i-1)%dVRad_dR(1)
            TR%R(i)%dVDisp_dR(0)=TR%R(i-1)%dVDisp_dR(1)
            TR%R(i)%dSigma_dR(0)=TR%R(i-1)%dSigma_dR(1)

            TR%R(i)%dVRot_dR(1)=(TR%R(i+1)%VRot-TR%R(i)%VRot)
     &                      /TR%R(i)%Rwidth
            TR%R(i)%dVRad_dR(1)=(TR%R(i+1)%VRad-TR%R(i)%VRad)
     &                      /TR%R(i)%Rwidth
            TR%R(i)%dVDisp_dR(1)=(TR%R(i+1)%VDisp-TR%R(i)%VDisp)
     &                      /TR%R(i)%Rwidth
            TR%R(i)%dSigma_dR(1)=(TR%R(i+1)%Sigma-TR%R(i)%Sigma)
     &                      /TR%R(i)%Rwidth

        endif
      enddo



      return
      end subroutine





ccccccc
c       This routine does a generalized conversion of a parameter vector to a tilted ring model
      subroutine TiltedRingOptionsToPV(PV,TRFO)
      implicit none
      Type(ParameterVector), INTENT(INOUT) :: PV
      Type(TiltedRingFittingOptions), target, INTENT(IN) :: TRFO
      integer i, CurrParam
      real,pointer,dimension(:) :: VariableVector

c      print*, "in generalize param vector"

      CurrParam=0
      ALLOCATE(VariableVector(0:TRFO%nRings-1))
c           Loop through all the specifc parameters
      do i=0,12
c           Associate the temporary pointers with the correct array
        if(i .eq.0) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%CentPos(0)
        elseif(i .eq.1) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%CentPos(1)
        elseif(i .eq.2) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%Inclination
        elseif(i .eq.3) then
            VariableVector=>
     &            TRFO%RadialProfiles(0:TRFO%nRings-1)%PositionAngle
        elseif(i .eq.4) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%VSys
        elseif(i .eq.5) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%VRot
        elseif(i .eq.6) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%VRad
        elseif(i .eq.7) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%VDisp
        elseif(i .eq.8) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%Vvert
        elseif(i .eq.9) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%dvdz
        elseif(i .eq.10) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%Sigma
        elseif(i .eq.11) then
            VariableVector=>
     &                  TRFO%RadialProfiles(0:TRFO%nRings-1)%z0
        elseif(i .eq.12) then
            VariableVector=>
     &            TRFO%RadialProfiles(0:TRFO%nRings-1)%zGradiantStart
        endif
c               Pass all the vectors to the vector assignment
        call SetParamLimsFromFittingOptions(i,PV%nParams,TRFO%nRings
     &          ,TRFO%ConstParams(i)
     &          ,TRFO%FixedParams(i),CurrParam
     &          ,PV%ParamLowerLims
     &          ,TRFO%ParamLowerLims)

        call SetParamLimsFromFittingOptions(i,PV%nParams,TRFO%nRings
     &          ,TRFO%ConstParams(i)
     &          ,TRFO%FixedParams(i),CurrParam
     &          ,PV%ParamUpperLims
     &          ,TRFO%ParamUpperLims)

        call SetParamFromFittingOptions(PV%nParams,TRFO%nRings
     &          ,TRFO%ConstParams(i)
     &          ,TRFO%FixedParams(i),CurrParam
     &          ,PV%Param,VariableVector)

c        print*, "Fit Check Ini",i,TRFO%ConstParams(i)
c     &                  ,TRFO%FixedParams(i)
c        print*, VariableVector

      enddo
c      print*, "Final Parameter vector", PV%Param(0:PV%nParams-1)
c     &              ,PV%nParams


      return
      end subroutine
ccccccc


cccccc
c           This routine steps through a parameter vector and sets the
c           values from a target vector based on the logical switchs
      subroutine SetParamLimsFromFittingOptions(ParamNum,nParam,nRings
     &                  ,Const,Fixed,CurrParam
     &                  ,ParamLims
     &                  ,TargLims)
      implicit none
      integer,INTENT(IN) :: ParamNum,nRings,nParam
      logical,INTENT(IN) :: Const,Fixed
      real,INTENT(INOUT) :: ParamLims(0:nParam-1)
      integer,INTENT(INOUT) :: CurrParam
      real,INTENT(IN) :: TargLims(0:12)

      integer i

      if(Fixed .eqv. .False.) then
        if(Const) then
            ParamLims(CurrParam)=TargLims(ParamNum)
        else
            ParamLims(CurrParam:CurrParam+nRings-1)=TargLims(ParamNum)
        endif
      endif

      return
      end subroutine
cccccccc




cccccc
c           This routine steps through a parameter vector and sets the
c           values from a target vector based on the logical switchs
      subroutine SetParamFromFittingOptions(nParam,nRings
     &                  ,Const,Fixed,CurrParam
     &                  ,Param
     &                  ,TargVec)
      implicit none
      integer,INTENT(IN) :: nRings,nParam
      logical,INTENT(IN) :: Const,Fixed
      real,INTENT(INOUT) :: Param(0:nParam-1)
      integer,INTENT(INOUT) :: CurrParam
      real,INTENT(INOUT) :: TargVec(0:nRings-1)

      integer i

      if(Fixed .eqv. .False.) then
        if(Const) then
            Param(CurrParam)=TargVec(0)
            CurrParam=CurrParam+1
        else
            Param(CurrParam:CurrParam+nRings-1)=TargVec(0:nRings-1)
            CurrParam=CurrParam+nRings
        endif
      endif

      return
      end subroutine
cccccccc



      end module
