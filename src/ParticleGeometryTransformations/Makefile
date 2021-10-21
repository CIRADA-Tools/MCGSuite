mFlags=../makeflags  
include $(mFlags)

all: $(ParticleGeometryObj) install

GalObjDef.o:
	$(MakeFObj) GalObjDef.f

clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
