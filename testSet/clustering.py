import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir+"/Clustering") 

import runClustering

for i in range(3):
    filename = f"sim{i+1}.root"
    outputFilename = f"sim{i+1}.csv"
    fEMinDamage =  5
    fEMaxDamage = 37.5
    probIndirect = 0.405
    sugarPosFilename = "sugarPos_4x4_300nm.csv"

    runClustering.runClustering(filename, outputFilename, fEMinDamage, fEMaxDamage, probIndirect, sugarPosFilename)