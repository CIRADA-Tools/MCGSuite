mFlags=../makeflags  
include $(mFlags)

all: $(OutputObj) install

SingleRingGeneration.o:
	$(MakeFObj) SingleRingGeneration.f



clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
