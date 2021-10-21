mFlags=../makeflags 
include $(mFlags)

VPATH :=$(BinLoc)

ProgramList= $(AllProgs)

all: $(ProgramList)


MockCubeGenerator : GenerateMockCube.o
	$(F77) $(FLAGS) -o MockCubeGenerator  GenerateMockCube.o $(AllObj) $(BSObj) $(FitsLibLoc) $(FFTW_LIB)


install:
	mv $(ProgramList) ../../Programs


clean:
	rm -f $(ProgramList);
	rm -f *.o;


