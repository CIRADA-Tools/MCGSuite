MCGSuite

version 1.0

---
authors: N. Deg, K. Spekkens
Developed by CIRADA (Canadian Institute for Radio Astronomy Data Analysis)
References: Lewis 2019, Spekkens et al. (in prep)

----
Mock Cube Generator Suite (MCGSuite) is a set of three codes designed to generate mock observations of tilted ring models.  These mock observations mimic observations of axisymmetric HI gas disks.  The codes are:
    
    MockCubeGenerator  (herafter MCG) -- a Fortran code that generates cubes from tilted-ring models
    
    make_galaxy_MCG -- a python code that constructs a flat-disk tilted ring model from scaling relations and uses MCG to produce a corresponding mock cube.
    
    make_suite_MCG -- a python code that constructs a suite of flat-disk tilted ring models from scaling relations and uses MCG to produce corresponding mock cubes. 

Full documentation can be found in Documentation/
In that folder, the file InstallationGuide.txt explains how to install the programs, while QuickstartGuide.txt explains the basics on how to run each code.  The MCG_Guide.pdf contains an outline of tilted ring models and explains how the MCG code works.  The MakeGalaxyGuide.pdf file explains how the make_galaxy_MCG and make_suite_MCG codes combine scaling relations with MCG to produce realistic mock observations.

-----
DEPENDENCIES:
MockCubeGenerator -- gfortran, cfitsio*, fftw3*
make_galaxy_MCG -- python, numpy, scipy, matplotlib, astropy
make_suite_MCG -- python, numpy, scipy, matplotlib, astropy

*For convenience, versions of the cfitsio and fftw3 libraries are contained in PSOFT14/MCGSuite/third_party

-----
Installation

0) Configure and install cfitsio and fftw3 according to their instructions.

1) Edit /src/makeflags so that 
    a) FitsLibLoc points to libcfitsio.a
    b) FFTW_* points to the fftw3 library.
    c) F77 is the Fortran compiler to be used.
    
2) In terminal run:
    cd src/
    make clean
    make 
    
    
----
Sample Outputs

A set of sample outputs can be accessed at:
https://www.dropbox.com/s/g0l1iu5jdzrgxr0/SampleOutputs.zip?dl=0

----
Basic Operations

For the purposes of this document, the user is assumed to be in main directory where README.md file is there.

--MCG--
MCG requires three text input files; a main input file, a data cube file, and a tilted ring model file.  Examples of these files can be found in /Inputs.

Assuming these files exist, the code is run in terminal via:
    ./Programs/MockCubeGenerator “MainInFile”

Using the sample files, this would be:
    ./Programs/MockCubeGenerator Inputs/MockCubeGeneratorInputs.in

As a test, the output from the sample input files should agree with the contents of
/SampleOutputs/MCG_SampleOut

--make_galaxy_MCG--
make_galaxy_MCG requires a number of things to run properly:

1) The MCG program is located in the Programs subfolder
2) That galaxy_config_MCG.py and observatory_config_MCG.py are in Inputs/ 
3) If using a beta configuration file, that it is specified in galaxy_config_MCG.py 

Assuming these conditions are met, the code is run in terminal via:

    ./make_galaxy_MCG

As a test, the given input files should produce outputs agreeing with 
/SampleOutputs/MakeGalaxy_SampleOutput/


--make_suite_MCG--

make_suite_MCG has similar requirements to make_galaxy_MCG:

1) The MCG program is located in the Programs subfolder
2) That suite_config_MCG.py and observatory_config_MCG.py are in Inputs/
3) If using a beta configuration file, that it is specified in suite_config_MCG.py 

Assuming these conditions are met, the code is run in terminal via:

    ./make_suite_MCG

As a test, the given input files should produce outputs agreeing with 
/SampleOutputs/MakeSuite_SampleOutput/




