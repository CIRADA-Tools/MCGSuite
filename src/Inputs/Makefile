mFlags=../makeflags  
include $(mFlags)

all: $(InputObj) install

InputUnitConversions.o:
	$(MakeFObj) InputUnitConversions.f

DataCubeInput.o:
	$(MakeFObj) DataCubeInput.f

TiltedRingModelInput.o:
	$(MakeFObj) TiltedRingModelInput.f

SofiaInputs.o:
	$(MakeFObj) SofiaInputs.f

FittingOptionsInputs.o:
	$(MakeFObj) FittingOptionsInputs.f


CubeGeneratorInputs.o: DataCubeInput.o TiltedRingModelInput.o InputUnitConversions.o
	$(MakeFObj) CubeGeneratorInputs.f


FullGalaxyInput.o: DataCubeInput.o
	$(MakeFObj) FullGalaxyInput.f

PipelineRuntimeInputs.o: SofiaInputs.o FittingOptionsInputs.o
	$(MakeFObj) PipelineRuntimeInputs.f


clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
