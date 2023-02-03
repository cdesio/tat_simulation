"""
Calculate expected range for testSet results

10 x 1,000 simulations

Error bars are the maximum deviation from the mean
"""

import numpy as np
import os
import sys
import inspect

clusteringDir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))) + "/Clustering/"
sys.path.insert(0, clusteringDir) 

from readResults import *

refValues = {
# sim number: [mean value, maximum difference from mean]
"sim1": {"LET": [], "direct": [], "indirect": [], "total": []}, # 9 MeV
"sim2": {"LET": [], "direct": [], "indirect": [], "total": []}, # 8 MeV
"sim3": {"LET": [], "indirect": [], "total": []} # 5 MeV
}

energy = {"sim1": "9", "sim2": "8", "sim3": "5"}

path = "."

for key in energy:
    dataset = [] 
    datasetDSB = [] 
    files = [a for a in os.listdir(path) if (f"alpha_{energy[key]}MeV" in a) and ('DSB' not in a) and ('csv' in a)]
    filesDSB = [a for a in os.listdir(path) if (f"alpha_{energy[key]}MeV" in a) and ('DSB' in a) and ('csv' in a)]
    files.sort()
    filesDSB.sort()

    for j in files:
        dataset.append(readIN('{}/{}'.format(path, j)))
    for j in filesDSB:
        datasetDSB.append(readIN('{}/{}'.format(path, j)))


    refValues[key]["LET"]=np.sum([getLET(dataset[a]) for a in range(len(dataset))])/len(dataset)

    totalSB_ = []
    directSB_ = []
    indirectSB_ = []


    for m in range(len(dataset)):
            totalSB_.append(calcBreaksperDose(f"TotalSBtotal",dataset[m]))
            directSB_.append(calcBreaksperDose(f"TotalSBdirect",dataset[m]))
            indirectSB_.append(calcBreaksperDose(f"TotalSBindirect",dataset[m]))

    refValues[key]["total"] = [np.mean(totalSB_), max([np.mean(totalSB_) - min(totalSB_),max(totalSB_)-np.mean(totalSB_)]) ] # mean, max variation from mean
    refValues[key]["direct"] = [np.mean(directSB_), max([np.mean(directSB_) - min(directSB_),max(directSB_)-np.mean(totalSB_)]) ] # mean, max variation from mean
    refValues[key]["indirect"] = [np.mean(indirectSB_), max([np.mean(indirectSB_) - min(indirectSB_),max(indirectSB_)-np.mean(indirectSB_)]) ] # mean, max variation from mean


print(refValues)


