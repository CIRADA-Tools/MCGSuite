mFlags=makeflags 
include $(mFlags)


all: LocalLib Programs

First:
	cd thirdparty….

LocalLib: 
	cd StandardMath; make;
	cd ObjectDefinitions; make;
	cd TiltedRingModelGeneration; make;
	cd TiltedRingToDataCube; make;
	cd ConvolveCube; make;
	cd GenerateNoiseCubes; make;
	cd ParameterToTiltedRingInterface; make;
	cd ParticleGeometryTransformations; make;
	cd Globals; make;
	cd Inputs; make;
	cd Outputs; make;
	cd ProgramMains; make;


Programs:
	cd $(LocBinLoc); make; make install



clean:
	cd StandardMath; make clean;
	cd ObjectDefinitions; make clean;
	cd TiltedRingModelGeneration; make clean;
	cd TiltedRingToDataCube; make clean;
	cd ConvolveCube; make clean;
	cd GenerateNoiseCubes; make clean;
	cd ParameterToTiltedRingInterface; make clean;
	cd ParticleGeometryTransformations; make clean;
	cd Inputs; make clean;
	cd Globals; make clean;
	cd GalaxyAnalysis; make clean;
	cd Outputs; make clean;
	cd ProgramMains; make clean;


	cd $(LocBinLoc); rm *.o;
	cd $(LocModLoc); rm *.mod
	cd ../Programs; rm $(AllProgs)

	

