cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the procedures for modeling
c       a galaxy using BBarolo once the pre-analysis is
c       complete.  It is designed to be an option
c       for the GeneralFitInterface (see GalaxyFit.f) to point
c       towards.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module BBaroloFitMod
      use ParameterVectorMod
      use DataCubeMod
      use SoFiACatalogueMod
      use BeamMod
      use CalcBeamKernelMod
      use DownhillSimplexMod
      use FullModelComparisonMod
      use GalaxyFitMod
      implicit none


      contains
cccccccc
c
      subroutine BBaroloFit(CatItem)
      use PipelineGlobals
      implicit none
      Type(CatalogueItem),INTENT(IN) :: CatItem

      real,ALLOCATABLE :: paramGuesses(:,:),chiArray(:)
      real chi2
      integer i,iter
      character(200) BBaroloName, script
      character(5) IDStr

      write(IDStr,"(I3.3)") CatItem%ObID
      BBaroloName=trim(ObjectBaseName)
     &      //trim(IDStr)//"_BBarolo.par"

      print*, "Fitting Galaxy with BBarolo"
      call MakeBBaroloInput(CatItem,BBaroloName)

      script="Programs/BBarolo -p "//trim(BBaroloName)
      call system(script)
c       Clean up the results
c      script='rm '//trim(BBaroloName)
c      call system(script)

      return
      end subroutine
ccccccccc


cccccccc
c
      subroutine MakeBBaroloInput(CatItem,BBaroloName)
      use PipelineGlobals
      implicit none
      Type(CatalogueItem),INTENT(IN) :: CatItem
      character(200) BBaroloName,LineStr
      integer nRings
      character(5) IDStr


      print*, "Making BBarolo Input File"

      open(10,file=BBaroloName,status='replace')
c       The name of the data cube file to analyze
      LineStr="FITSFILE    "//trim(DataCubeFile)
      write(10,*) trim(LineStr)
c       The namke of the mask file to use
      LineStr="MASK   "//trim(MaskName)
      write(10,*) trim(LineStr)
c       The number of threads to use
      LineStr="THREADS     1"
      write(10,*) trim(LineStr)
c       The output folder for the BBarolo analysis
      write(IDStr,"(I3.3)") CatItem%ObID
      LineStr="OUTFOLDER   "//trim(ObjectBaseName)//"_Outputs_"
     &      //trim(IDStr)
      write(10,*) trim(LineStr)
c       Make sure we are using the 3D Fit
      LineStr="3DFIT       true"
      write(10,*) trim(LineStr)
c       Turn off the errors
      LineStr="flagErrors  false"
      write(10,*) trim(LineStr)
c       The Weight function
      LineStr="WFUNC       0"
      write(10,*) trim(LineStr)
c       The L-type
      LineStr="LTYPE       2"
      write(10,*) trim(LineStr)
c       The F type
      LineStr="FTYPE       2"
      write(10,*) trim(LineStr)
c       Use two-stage fitting
      LineStr="TWOSTAGE    true"
      write(10,*) trim(LineStr)
c       The number of radial bins to use
      nRings=TR_FittingOptions%nRings
      LineStr="NRADII"
      write(10,*) trim(LineStr), TR_FittingOptions%nRings
c       The size of the radial bins
      LineStr="RADSEP"
      write(10,*) trim(LineStr), ObservedBeam%BeamMajorAxis
c       Initial guess for Vsys
      LineStr="VSYS"
      write(10,*) trim(LineStr)
     &          , TR_FittingOptions%RadialProfiles(0)%VSys
c       Initial guess for inclination
      LineStr="INC"
      write(10,*) trim(LineStr)
     &          , TR_FittingOptions%RadialProfiles(0)%Inclination
     &          *180./Pi
c       Initial guess for Position Angle
      LineStr="PA"
      write(10,*) trim(LineStr)
     &          , TR_FittingOptions%RadialProfiles(0)%PositionAngle
     &          *180./Pi-90.
c       Initial Guess for central x-position in pixels
      LineStr="XPOS"
      write(10,*) trim(LineStr)
     &              ,ObservedDC%DH%PixelCenterIndx(0)
c     &          , (TR_FittingOptions%RadialProfiles(0)%CentPos(0)
c     &          -ObservedDC%DH%Start(0))
c     &          /abs(ObservedDC%DH%PixelSize(0))
c      print*, "Cent check"
c     &      ,TR_FittingOptions%RadialProfiles(0)%CentPos(0)
c     &      ,ObservedDC%DH%Start(0),ObservedDC%DH%PixelSize(0)
c     &      ,ObservedDC%DH%nPixels(0)
c       Initial Guess for central y-position in pixels
      LineStr="YPOS"
      write(10,*) trim(LineStr)
     &              ,ObservedDC%DH%PixelCenterIndx(1)
c     &          , TR_FittingOptions%RadialProfiles(0)%CentPos(1)
c&          *180./Pi-90.
c       Initial Guess for VRot
      LineStr="VROT"
      write(10,*) trim(LineStr)
     &          , TR_FittingOptions%RadialProfiles(nRings-1)%VRot

c       Set up the free parameters
      LineStr="FREE"
      if(TR_FittingOptions%FixedParams(0) .eqv. .False.) then
        LineStr=trim(LineStr)//" XPOS"
      endif
      if(TR_FittingOptions%FixedParams(1) .eqv. .False.) then
        LineStr=trim(LineStr)//" YPOS"
      endif
      if(TR_FittingOptions%FixedParams(2) .eqv. .False.) then
        LineStr=trim(LineStr)//" INC"
      endif
      if(TR_FittingOptions%FixedParams(3) .eqv. .False.) then
        LineStr=trim(LineStr)//" PA"
      endif
      if(TR_FittingOptions%FixedParams(4) .eqv. .False.) then
        LineStr=trim(LineStr)//" VSYS"
      endif
      if(TR_FittingOptions%FixedParams(5) .eqv. .False.) then
        LineStr=trim(LineStr)//" VROT"
      endif
      if(TR_FittingOptions%FixedParams(6) .eqv. .False.) then
        LineStr=trim(LineStr)//" VRAD"
      endif
      if(TR_FittingOptions%FixedParams(7) .eqv. .False.) then
        LineStr=trim(LineStr)//" VDISP"
      endif


      write(10,*) trim(LineStr)
      close(10)


      return
      end subroutine
ccccccccc




      end module
