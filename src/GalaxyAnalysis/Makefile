mFlags=../makeflags  
include $(mFlags)

all: $(GalaxyAnalysisObj) install


GalaxyFit.o:
	$(MakeFObj) GalaxyFit.f

BBaroloFit.o: GalaxyFit.o
	$(MakeFObj) BBaroloFit.f


clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
