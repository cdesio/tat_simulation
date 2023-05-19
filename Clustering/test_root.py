from ROOT import TFile
import numpy as np

fFile = TFile.Open('../output/test_At100k_160k_boxes_10.5um/out_AtDNA_100k_spacing_2um_5115.root')
fDirectory = fFile.Get("ntuple")
tEdep = fDirectory.Get("EventEdep")
tDirect = fDirectory.Get("Direct")
num_entries = tDirect.GetEntries()
values = np.empty(num_entries, dtype=float)

for i in range(num_entries):
    tDirect.GetEntry(i)
    values[i] = tDirect.x
print(len(values))
print("done.")