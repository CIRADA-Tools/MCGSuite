mFlags=../makeflags  
include $(mFlags)

all: $(ConvolutionObj) install


CalculateBeamKernel.o:
	$(MakeFObj) CalculateBeamKernel.f

TwoDConvolution.o:
	$(MakeFObj) TwoDConvolution.f

CubeKernelConvolution.o: CalculateBeamKernel.o TwoDConvolution.o
	$(MakeFObj) CubeKernelConvolution.f

clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
