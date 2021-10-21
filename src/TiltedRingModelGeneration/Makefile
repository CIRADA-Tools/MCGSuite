mFlags=../makeflags  
include $(mFlags)

all: $(TiltedRingGenerationObj) install

SingleRingGeneration.o:
	$(MakeFObj) SingleRingGeneration.f


TiltedRingModelGeneration.o: SingleRingGeneration.o
	$(MakeFObj) TiltedRingModelGeneration.f

clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
