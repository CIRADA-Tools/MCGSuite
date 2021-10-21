mFlags=../makeflags  
include $(mFlags)

all: $(TiltedRingInterfaceObj) install

ParameterToTiltedRingVector.o:
	$(MakeFObj) ParameterToTiltedRingVector.f



clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
