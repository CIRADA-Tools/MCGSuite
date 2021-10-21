

      module CubeGeneratorGlobals
c               Use the various object definition routines
      use TiltedRingMod
      use DataCubeMod
      use BeamMod
      use ParameterVectorMod

      implicit none

      Type(DataCube) ModelDC, DataCubeMask, NoiseCube

      Type(TiltedRingModel) ModelTiltedRing

      Type(Beam2D) ModelBeam

      integer idum

      integer OutputFileVolumeSwitch

      character(500) CubeHeaderFile, TiltedRingInputFile,BeamInputFile
      character(200) OutputBaseNames
      character(300) OutputFileNames(0:3)

      real MaskingValue

      end module CubeGeneratorGlobals
