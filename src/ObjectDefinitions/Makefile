mFlags=../makeflags  
include $(mFlags)

all: $(ObjectDefObj) install

Particle.o:
	$(MakeFObj) Particle.f

TiltedRing.o: Particle.o
	$(MakeFObj) TiltedRing.f

clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
