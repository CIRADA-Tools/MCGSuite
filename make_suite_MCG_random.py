#!/usr/bin/env python3
from make_galaxy_code import *

#from joblib import Parallel, delayed
import multiprocessing as mp
import numpy as np
from multiprocessing import freeze_support

def Main():
    print("Making a suite of galaxies using MCG")
    #   Get the inputs and set up the objects
    Suite=IN.Get_MakeSuite_Inputs_random()    #From config files
    #   Create a combined array of all combinations of the suite parameters
    OC.RandomSuiteConfig(Suite)

    """
    #   Set up parallel processing
    print("n Processors", Suite.nProcessors)

        #   Run the main loop for the processing
    pool=mp.Pool(processes=Suite.nProcessors)
    Suite.DBTable=pool.starmap(SO.SuiteMainLoopFn_random, [(i,Suite) for i in range(Suite.n_galaxies)])
    #i=0
    #SO.SuiteMainLoopFn_random(i,Suite)

    SO.CatalogueOutput(Suite)
    SO.CatalogueOutput_csv(Suite)
    """

if __name__=="__main__":
    freeze_support()
    Main()

