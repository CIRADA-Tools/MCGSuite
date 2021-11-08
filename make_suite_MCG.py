#!/usr/bin/env python3
from src.make_galaxy_code import *

#from joblib import Parallel, delayed
import multiprocessing as mp
import numpy as np
from multiprocessing import freeze_support

def Main():
    print("Making a suite of galaxies using MCG")
    #   Get the inputs and set up the objects
    Suite=IN.Get_MakeSuite_Inputs()    #From config files

    #   Do the safety checks to make sure all parameters are within acceptable limits
    SC.SuiteChecks(Suite,Suite.Templates[1])
    #   Create a combined array of all combinations of the suite parameters
    OC.SuiteConfig(Suite)

    #   Set up parallel processing
    print("n Processors", Suite.nProcessors)

        #   Run the main loop for the processing
    pool=mp.Pool(processes=Suite.nProcessors)
    Suite.DBTable=pool.starmap(SO.SuiteMainLoopFn, [(i,Suite) for i in range(Suite.n_galaxies)])


    SO.CatalogueOutput(Suite)


if __name__=="__main__":
    freeze_support()
    Main()


