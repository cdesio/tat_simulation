import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import runClustering

filename = "testFile.root"
outputFilename = "results.csv"
fEMinDamage =  5
fEMaxDamage = 37.5
probIndirect = 0.405
sugarPosFilename = "sugarPosTest.csv"

runClustering.runClustering(filename, outputFilename, fEMinDamage, fEMaxDamage, probIndirect, sugarPosFilename)