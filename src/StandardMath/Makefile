mFlags=../makeflags  
include $(mFlags)

all: $(StandardMathObj) install


clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
