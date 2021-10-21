ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   This is the main routine for the the WALLABY Resolved Kinematic Pipeline.
c   It is open source, but requires the openmpi and fftw3 libraries.
c
c   The code is meant to take in a catalogue of resolved detections and
c   generate a rotating disk model for each object.
c
c   The code is focused on low-resolution images, but it is meant to be
c   flexible and easy for the user to modify.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program MakeMockCube
c       Modules neeeded at this level
      use DataCubeOutputsMod
      use CalcBeamKernelMod
      use CubeGeneratorGlobals
      use CubeGeneratorInputMod

      implicit none
      character(500) script

      real BeamArea,BeamPixels

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      print*, "Mock Cube Generator"
c       Get all the inputs
      call MockCubeGeneratorIn()    !/src/Inputs/CubeGeneratorInputs.f
c       Initialize the various output files
      call IniOutputFiles()         !This File
c           Write out the names of the cube output units
      call NameCubeUnits(ModelDC)   !This File
c       Make the cube itself
      call BuildSourceCube()        !This File
c       Make the Beam
      call Calculate2DBeamKernel(ModelBeam,ModelDC%DH%PixelSize)    !/src/ConvolveCube/CalculateBeamKernel.f
c       Do the cube convolutions (velocity and spatial)
      call MakeConvolvedCube()      !This File
c       Make the Noise Cube
      call CreateNoiseCube()        !This File
c       Sum the noise with the flux
      ModelDC%Flux=ModelDC%Flux+NoiseCube%Flux
c       Output the new signal+noise cube
      call WriteDataCubeToFITS(ModelDC,ModelBeam,OutputFileNames(3))    !/src/Outputs/DataCubeOutputs/f
      script="mv "//trim(OutputFileNames(3))//" "//trim(OutputBaseNames)
c     &              //"/"
      call system(trim(script))


      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccc
c           This routine builds the source cube
      subroutine BuildSourceCube()
      use CubeGeneratorGlobals
      use TiltedRingGenerationMod
      use FillDataCubeWithTiltedRingMod
      use DataCubeOutputsMod

c      use TiltedRingOutputsMod

      implicit none
      real pixelarea
      character(500) script

c       Calculate the area of a pixel in arcsec^2
      pixelarea=abs(ModelDC%DH%PixelSize(0)*ModelDC%DH%PixelSize(1))
c           First buld the model
      call BuildTiltedRingModel(ModelTiltedRing,pixelarea,idum) !/src/TiltedRingModelGeneration/TiltedRingModelGeneration.f

c      call TiltedRing_ParticleOutputs(ModelTiltedRing
c     &          ,"TiltedRingParticles.txt")
c       Create the point-source data cube
      call FillDataCubeWithTiltedRing(ModelDC,ModelTiltedRing) !/src/TiltedRingToDataCube/FillDataCubeByTiltedRing.f
c       Save the source cube to a file
      if(OutputFileVolumeSwitch .eq. 0) then
        return
      else
        call WriteDataCubeToFITS(ModelDC,ModelBeam,OutputFileNames(0)) !/src/Outputs/DataCubeOutputs/f
c       Move the source cube into the target directory
        script="mv "//trim(OutputFileNames(0))
     &                 //" "//trim(OutputBaseNames)//"/"
        call system(trim(script))
      endif
      return
      end subroutine
cccccc

cccccc
c       This routine runs the convolutions and outputs the results (depending on the switch)
      subroutine MakeConvolvedCube()
      use CubeGeneratorGlobals
      use TwoDConvolutionMod
      use CubeKernelConvolutionMod
      use VelSmoothMod
      use DataCubeOutputsMod
      implicit none
      character(500) script
      real BeamArea, BeamPixels
c       Do the velocity smoothing
      call SmoothVelChannels(ModelDC,ModelBeam) !/src/ConvolveCube/VelocitySmoothing.f
c       If being very verbose output the velocity smoothed source cube
      if(OutputFileVolumeSwitch .eq. 3) then
        call WriteDataCubeToFITS(ModelDC,ModelBeam,OutputFileNames(1)) !/src/Outputs/DataCubeOutputs/f
c       Move the source cube into the target directory
        script="mv "//trim(OutputFileNames(1))//
     &                  " "//trim(OutputBaseNames)//"/"
        call system(trim(script))
      endif
c       Now convolve the cube with the beam
      call CubeBeamConvolution(ModelDC,ModelBeam) !/src/ConvolveCube/CubeKernelConvolution.f
c       Re-scale the source cube flux
      BeamArea=2.*Pi*(ModelBeam%BeamSigmaVector(0)
     &                      *ModelBeam%BeamSigmaVector(1))
c      BeamPixels=BeamArea/(ModelDC%DH%PixelSize(0)
c     &              *ModelDC%DH%PixelSize(1))
      BeamPixels=BeamArea
      ModelDC%Flux=ModelDC%Flux*BeamPixels

c       If being verbose output the convolved source cube
      if(OutputFileVolumeSwitch .ge. 2 ) then
        call WriteDataCubeToFITS(ModelDC,ModelBeam,OutputFileNames(2))
c       Move the source cube into the target directory
        script="mv "//trim(OutputFileNames(2))
     &                  //" "//trim(OutputBaseNames)//"/"
        call system(trim(script))
      endif

      return
      end subroutine
cccccc

cccccc
c           This routine makes the convolved noise cube
      subroutine CreateNoiseCube()
      use CubeGeneratorGlobals
      use NoiseCubeMod
      implicit none
      NoiseCube%DH=ModelDC%DH
      call AllocateDataCube(NoiseCube) !/src/ObjectDefinitions/DataCube.f
      call MakeConvolvedNoiseCube(NoiseCube,MOdelBeam,idum) !/src/GenerateNoiseCubes/MakeNoiseCube.f
      return
      end subroutine
ccccccc


ccccccc
c           This routine makes the output directory and names the various output cubes
      subroutine IniOutputFiles()
      use CubeGeneratorGlobals
      implicit none
      character(500) script
c       Make the output directory
      script="mkdir "//trim(OutputBaseNames)
      call system(trim(script))
c       Name each of the output cubes
      OutputFileNames(0)=trim(OutputBaseNames)
     &                  //"_SourceCube.fits"        !The source cube
      OutputFileNames(1)=trim(OutputBaseNames)
     &                  //"_VelocitySmoothedSourceCube.fits"    !Convolved with the velocity smooth
      OutputFileNames(2)=trim(OutputBaseNames)
     &                  //"_ConvolvedSourceCube.fits"       !Convolved with the PSF
      OutputFileNames(3)=trim(OutputBaseNames)          !With the convolved noise added to the source
     &                  //".fits"
      return
      end subroutine
ccccccccccc


ccccccccc
c           This routine sets the units for the output cube header
      subroutine NameCubeUnits(DC)
      use DataCubeMod
      implicit none
      Type(DataCube), INTENT(INOUT) :: DC

      DC%DH%Units(0)='DEGREES'
      DC%DH%Units(1)='DEGREES'
      DC%DH%Units(2)='m/s    '

      DC%DH%AxisType(0)='RA---SIN'
      DC%DH%AxisType(1)='DEC--SIN'
      DC%DH%AxisType(2)='VELO-LSR'
      DC%DH%FUnit='Jy/Beam '
      DC%DH%FType='intensity'
      DC%DH%Epoch=2000.
      return
      end subroutine
cccccc
