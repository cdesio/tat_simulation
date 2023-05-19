import os
import sys
import inspect
from typing import Union
import numpy as np
from scipy.spatial import cKDTree
from ROOT import TFile

# To import clustering from build folder
currentdir = os.path.dirname(os.path.abspath(
    inspect.getfile(inspect.currentframe())))
builddir = currentdir+r"/build"
sys.path.insert(0, builddir)


from clustering import clustering

def IsEdepSufficient(pEdep: float, fEMinDamage: float, fEMaxDamage: float) -> bool:
    if (pEdep < fEMinDamage):
        return False
    if (pEdep > fEMaxDamage):
        return True
    else:
        proba = (pEdep - fEMinDamage)/(fEMaxDamage-fEMinDamage)
        return (proba > np.random.rand())


def checkPoint(point: list, T0: cKDTree, T1: cKDTree) -> list:
    result = [-1, -1]
    Rdirect = 0.35  # nm
    d0, idx0 = T0.query(point, k=1)
    d1, idx1 = T1.query(point, k=1)

    if d0 < d1:
        currentDist = d0
        currentStrand = 0
        currentCopy = idx0
    else:
        currentDist = d1
        currentStrand = 1
        currentCopy = idx1

    if currentDist < Rdirect:
        result[0] = currentStrand
        result[1] = currentCopy

    return result


def getIndex(point: list, T0: cKDTree, T1: cKDTree) -> list:

    d0, idx0 = T0.query(point, k=1)
    d1, idx1 = T1.query(point, k=1)

    d = min(d0, d1)
    idx = idx0 if d0 < d1 else idx1
    strand = 0 if d0 < d1 else 1

    if d > 1e-10:
        raise ValueError

    return [idx, strand]


def readClusteringGitHash(builddir):
    with open(builddir+"/git-state.txt", "r") as f:
        data = f.readlines()
    return data[0]


def runClustering(filename: str, outputFilename: str, fEMinDamage: float, fEMaxDamage: float, probIndirect: float, sugarPosFilename: str, filenamePhoton: Union[str, bool] = False, separate_r=False, n_boxes = 3200, boxes_per_R = 800):
    print("runClustering OLD VERSION")
    with open(sugarPosFilename, "r") as f:
        data = f.readlines()

    data = [[float(b)*1E6 for b in a.split("\t")] for a in data]

    sugar0 = []
    sugar1 = []

    for line in data:
        sugar0.append(line[0:3])
        sugar1.append(line[3:6])

    sugar0 = np.array(sugar0)
    sugar1 = np.array(sugar1)

    T0 = cKDTree(sugar0)
    T1 = cKDTree(sugar1)

    # Load data to analyse
    fFile = TFile.Open(filename)
    fDirectory = fFile.Get("ntuple")

    tInfo = fDirectory.Get("Info")
    tInfo.GetEntry(0)
    chromatinVolume = tInfo.ChromatinVolume_m3  # in m3
    numBP = tInfo.NumBasepairs
    if filenamePhoton:
        LET = "N/A"
    else:
        if hasattr(tInfo,"MeanLET"):
            LET = tInfo.MeanLET
        else:
            LET = "N/A"

    gitHash = tInfo.GitHash

    assert (len(data) == numBP)

    tEdep = fDirectory.Get("EventEdep")

    pathLength = {}
    dosePerEvent = {}
    meanKEperEvent = {}

    entryEdepNumber = tEdep.GetEntries()
    
    # build mapping to decay event number
    eventInfo = fDirectory.Get("Events")
    nentriesMapping = eventInfo.GetEntries() #get number of events for mapping
    mapping = {}  # photon event ID to simulation event ID dict
    for irow in range(0, nentriesMapping): #loop on nentries
        eventInfo.GetEntry(irow)            
        mapping[eventInfo.EventNo] = eventInfo.PhotonEventID #save photon event number into key simulation event number
    # end mapping 

    numEvt = 0
    for i in range(entryEdepNumber):
        tEdep.GetEntry(i)
        if tEdep.EventNo > numEvt:
            simEvt = tEdep.EventNo
            if simEvt in mapping.keys():
                numEvt = mapping[simEvt]
                if hasattr(tEdep, "PathLengthChromatin"):  # older results do not have path length
                    pathLength[numEvt] = tEdep.PathLengthChromatin
                else:
                    pathLength[numEvt] = "N/A"

                if hasattr(tEdep, "PrimaryKEEntrance"):  # older results do not have path length
                    meanKEperEvent[numEvt] = (
                        tEdep.PrimaryKEEntrance+tEdep.PrimaryKEExit)/2

                    #meanKEperEvent[tEdep.EventNo] = "N/A"

                if tEdep.Edep_J > 0:
                    if numEvt in dosePerEvent.keys():
                        dosePerEvent[numEvt] += tEdep.Edep_J / \
                        (1000 * chromatinVolume)
                        meanKEperEvent[numEvt] += tEdep.Edep_J
                    else:
                        dosePerEvent[numEvt] = tEdep.Edep_J / \
                        (1000 * chromatinVolume)
                        meanKEperEvent[numEvt] = tEdep.Edep_J
            else:
                if tEdep.Edep_J>0:
                    print(f"evt {simEvt} not in mapping but Edep>0")

    if np.any(list(meanKEperEvent.values()) != "N/A"):
        energy = sum(meanKEperEvent.values())/len(meanKEperEvent.values())
    else:
        energy = "N/A"

    numEvt += 1  # zero indexed
    
    #print(meanKEperEvent)
    if filenamePhoton:
        dosePerEventPhoton = {}
        NumIntersecting = 0
        fFilePhoton = TFile.Open(filenamePhoton)
        fDirectoryPhoton = fFilePhoton.Get("ntuple")
        tEdepPhoton = fDirectoryPhoton.Get("EventEdep")
        entryEdepNumberPhoton = tEdepPhoton.GetEntries()
        for i in range(0, entryEdepNumberPhoton):
            tEdepPhoton.GetEntry(i)
            if tEdepPhoton.Edep_J > 0:
                dosePerEventPhoton[tEdepPhoton.EventNo] = tEdepPhoton.Edep_J / \
                    (1000 * chromatinVolume)
        fInfoPhoton = fDirectoryPhoton.Get("Info")
        entryEdepNumberPhoton = fInfoPhoton.GetEntries()
        for i in range(0, entryEdepNumberPhoton):
            fInfoPhoton.GetEntry(i)
            NumIntersecting += fInfoPhoton.NumIntersecting

    # // Read out direct damage
    input_tree = fDirectory.Get("Direct")

    nentries = input_tree.GetEntries()

    cumulatedEnergyDep = {}

    for irow in range(nentries):
        input_tree.GetEntry(irow)
        result = checkPoint([input_tree.x, input_tree.y, input_tree.z], T0, T1)
        if (result[0] != -1):
            if input_tree.EventNo in mapping.keys():
                key = (mapping[input_tree.EventNo], result[0], result[1])
                if key in cumulatedEnergyDep:
                    cumulatedEnergyDep[key] += input_tree.eDep_eV
                else:
                    cumulatedEnergyDep[key] = input_tree.eDep_eV
            else:

                print(f"DIRECT: evt:{input_tree.EventNo} not in mapping but results[0]: {result[0]}")
        #  Add energy deposition from photon simulation
    if filenamePhoton:
        input_tree_photon = fDirectoryPhoton.Get("Direct")
        nentriesPhoton = input_tree_photon.GetEntries()

        cumulatedEnergyDepPhoton = {}

        for irow in range(nentriesPhoton):
            input_tree_photon.GetEntry(irow)

            result = checkPoint(
                [input_tree_photon.x, input_tree_photon.y, input_tree_photon.z], T0, T1)
            if (result[0] != -1):
                key = (input_tree_photon.EventNo, result[0], result[1])
                if key in cumulatedEnergyDepPhoton:
                    cumulatedEnergyDepPhoton[key] += input_tree_photon.eDep_eV
                else:
                    cumulatedEnergyDepPhoton[key] = input_tree_photon.eDep_eV

        eventInfo = fDirectory.Get("Events")
        nentriesMapping = eventInfo.GetEntries()
        mapping = {}  # photon event ID to simulation event ID
        for irow in range(0, nentriesMapping):
            eventInfo.GetEntry(irow)
            mapping[eventInfo.EventNo] = eventInfo.PhotonEventID
        #print(f"mapping photon: {mapping}")
        for key in cumulatedEnergyDep:
            if ((mapping[key[0]], key[1], key[2]) in cumulatedEnergyDepPhoton):
                cumulatedEnergyDep[key] += cumulatedEnergyDepPhoton[mapping[key[0]], key[1], key[2]]
        #print(f"cumulatedEnergyDep: {cumulatedEnergyDep}")
        # Add photon dose
        for key in dosePerEvent:
            if mapping[key] in dosePerEventPhoton:
                dosePerEvent[key] += dosePerEventPhoton[mapping[key]]

    eventsListDirect = []
    copyListDirect = []
    strandListDirect = []

    
    for key in cumulatedEnergyDep:
        if IsEdepSufficient(cumulatedEnergyDep[key], fEMinDamage, fEMaxDamage):

            eventsListDirect.append(key[0])
            
            copyListDirect.append(key[2])
            strandListDirect.append(key[1])
    
    # // Read out indirect damage
    input_tree_indirect = fDirectory.Get("Indirect")

    eventsListIndirect = []
    copyListIndirect = []
    strandListIndirect = []

    nentries_indirect = input_tree_indirect.GetEntries()

    for irow in range(nentries_indirect):
        input_tree_indirect.GetEntry(irow)
        if input_tree_indirect.EventNo in mapping.keys():
            if hasattr(input_tree_indirect, "copyNum"):
                if np.random.rand() <= probIndirect:
                    eventsListIndirect.append(mapping[input_tree_indirect.EventNo])
                    copyListIndirect.append(input_tree_indirect.copyNum)
            elif not hasattr(input_tree_indirect, "DNAmolecule"):
                if np.random.rand() <= probIndirect:
                    eventsListIndirect.append(mapping[input_tree_indirect.EventNo])
                    point = [input_tree_indirect.x,
                            input_tree_indirect.y, input_tree_indirect.z]
                    c, s = getIndex(point, T0, T1)
                    copyListIndirect.append(c)
                    strandListIndirect.append(s)
            else:
                if ((input_tree_indirect.DNAmolecule[0] == "D") and (input_tree_indirect.radical == "OH^0")):
                    if np.random.rand() <= probIndirect:
                        eventsListIndirect.append(mapping[input_tree_indirect.EventNo])
                        point = [input_tree_indirect.x,
                                input_tree_indirect.y, input_tree_indirect.z]
                        c, s = getIndex(point, T0, T1)
                        copyListIndirect.append(c)
                        strandListIndirect.append(s)
        else:
            print("INDIRECT check evt")
    # clustering
    evts = list(set(eventsListDirect+eventsListIndirect))
    

    tempResults = clustering(evts, eventsListDirect, copyListDirect,
                             strandListDirect, eventsListIndirect, copyListIndirect, strandListIndirect)

    clusteringResults = tempResults[0]  # strand break number results
    clusterSize = tempResults[1]  # DSB cluster size, cluster of size 1-10

    events = list(dosePerEvent.keys())
    eventsWithClusteringResults = [a[0] for a in clusteringResults]
    
    #print(eventsWithClusteringResults)
    
    clusteringGitHash = readClusteringGitHash(builddir)

    print("Finished: {}".format(filename))

    # selecting copynumbers to save results in separate files
    eventInfo = fDirectory.Get("Events")
    nentriesMapping = eventInfo.GetEntries()
    mapping_copyNo = {}  # copyNo to simulation event ID
    for irow in range(0, nentriesMapping):
        eventInfo.GetEntry(irow)
        mapping_copyNo[eventInfo.PhotonEventID] = eventInfo.copyNo
    ranges_radii = {}
    cnumbers = np.arange(0, n_boxes+boxes_per_R, boxes_per_R)
    n_files = n_boxes/boxes_per_R
    for r, (cmin, cmax) in enumerate(zip(cnumbers, cnumbers[1:])):
        ranges_radii[r] = range(cmin, cmax)
    assert(len(ranges_radii)==n_files)
    #print(f"dose: {dosePerEvent}\n pathlength: {pathLength}\n meanKE: {meanKEperEvent}\n")
    #print(f"mapping_copyno: {mapping_copyNo}, ranges_radii: {ranges_radii}")
    # end
    if separate_r:
        print("separating output files per radius")
        for r in ranges_radii:
            with open(os.path.splitext(outputFilename)[0]+f"_{r}.csv", "w") as f:
                f.write(
                    "Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,GitHash,ClusteringGitHash\n")
                if filenamePhoton:
                    f.write("{},{},{},{},{},{},{},{},{}\n".format(filename, energy, LET, NumIntersecting,
                            chromatinVolume, numBP, sugarPosFilename, gitHash, clusteringGitHash))
                else:
                    f.write("{},{},{},{},{},{},{},{},{}\n".format(filename, energy, LET, len(
                        pathLength.keys()), chromatinVolume, numBP, sugarPosFilename, gitHash, clusteringGitHash))
                f.write("EventNo,DoseGy,PathLength_nm,MeanKEMeV,TotalSBdirect,SSBdirect,cSSBdirect,DSBdirect,TotalSBindirect,SSBindirect,cSSBindirect,DSBindirect,TotalSBtotal,SSBtotal,cSSBtotal,DSBtotal\n")
                for event in events:
                    if mapping_copyNo[event] in ranges_radii[r]:
                        if event in eventsWithClusteringResults:
                            text = str([clusteringResults[eventsWithClusteringResults.index(
                                event)][1:]]).strip("[").strip("]")
                        else:
                            text = "0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0"
                        f.write("{},{},{},{},{}".format(
                            event, 
                            dosePerEvent[event], 
                            pathLength[event], 
                            meanKEperEvent[event], text))
                        f.write('\n')

            with open(os.path.splitext(outputFilename)[0]+f"_DSB_{r}.csv", "w") as f:
                f.write(
                    "Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,SimGitHash,ClusteringGitHash\n")
                if filenamePhoton:
                    f.write("{},{},{},{},{},{},{},{},{}\n".format(filename, energy, LET, NumIntersecting,
                            chromatinVolume, numBP, sugarPosFilename, gitHash, clusteringGitHash))
                else:
                    f.write("{},{},{},{},{},{},{},{},{}\n".format(filename, energy, LET, len(
                        events), chromatinVolume, numBP, sugarPosFilename, gitHash, clusteringGitHash))
                f.write("EventNo,DoseGy,PathLength_nm,Direct_1_SBperDSBcluster,Direct_2_SBperDSBcluster,Direct_3_SBperDSBcluster,Direct_4_SBperDSBcluster,Direct_5_SBperDSBcluster,Direct_6_SBperDSBcluster,Direct_7_SBperDSBcluster,Direct_8_SBperDSBcluster,Direct_9_SBperDSBcluster,Direct_10_SBperDSBcluster,Indirect_1_SBperDSBcluster,Indirect_2_SBperDSBcluster,Indirect_3_SBperDSBcluster,Indirect_4_SBperDSBcluster,Indirect_5_SBperDSBcluster,Indirect_6_SBperDSBcluster,Indirect_7_SBperDSBcluster,Indirect_8_SBperDSBcluster,Indirect_9_SBperDSBcluster,Indirect_10_SBperDSBcluster,Hybrid_1_SBperDSBcluster,Hybrid_2_SBperDSBcluster,Hybrid_3_SBperDSBcluster,Hybrid_4_SBperDSBcluster,Hybrid_5_SBperDSBcluster,Hybrid_6_SBperDSBcluster,Hybrid_7_SBperDSBcluster,Hybrid_8_SBperDSBcluster,Hybrid_9_SBperDSBcluster,Hybrid_10_SBperDSBcluster,Mixed_1_SBperDSBcluster,Mixed_2_SBperDSBcluster,Mixed_3_SBperDSBcluster,Mixed_4_SBperDSBcluster,Mixed_5_SBperDSBcluster,Mixed_6_SBperDSBcluster,Mixed_7_SBperDSBcluster,Mixed_8_SBperDSBcluster,Mixed_9_SBperDSBcluster,Mixed_10_SBperDSBcluster,Total_1_SBperDSBcluster,Total_2_SBperDSBcluster,Total_3_SBperDSBcluster,Total_4_SBperDSBcluster,Total_5_SBperDSBcluster,Total_6_SBperDSBcluster,Total_7_SBperDSBcluster,Total_8_SBperDSBcluster,Total_9_SBperDSBcluster,Total_10_SBperDSBcluster\n")
                for event in events:
                    if mapping_copyNo[event] in ranges_radii[r]:
                        if event in eventsWithClusteringResults:
                            text = str([clusterSize[eventsWithClusteringResults.index(event)][1:]]).strip(
                                "[").strip("]")
                        else:
                            text = "0, "*49+"0"

                        f.write("{},{},{},{}".format(
                            event, dosePerEvent[event], pathLength[event], text))
                        f.write('\n')
    else:
        with open(outputFilename, "w") as f:
            f.write(
                "Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,GitHash,ClusteringGitHash\n")
            if filenamePhoton:
                f.write("{},{},{},{},{},{},{},{},{}\n".format(filename, energy, LET, NumIntersecting,
                        chromatinVolume, numBP, sugarPosFilename, gitHash, clusteringGitHash))
            else:
                f.write("{},{},{},{},{},{},{},{},{}\n".format(filename, energy, LET, len(
                    pathLength.keys()), chromatinVolume, numBP, sugarPosFilename, gitHash, clusteringGitHash))
            f.write("EventNo,DoseGy,PathLength_nm,MeanKEMeV,TotalSBdirect,SSBdirect,cSSBdirect,DSBdirect,TotalSBindirect,SSBindirect,cSSBindirect,DSBindirect,TotalSBtotal,SSBtotal,cSSBtotal,DSBtotal\n")
            for event in events:
                if event in eventsWithClusteringResults:
                    text = str([clusteringResults[eventsWithClusteringResults.index(
                        event)][1:]]).strip("[").strip("]")
                else:
                    text = "0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0"
                f.write("{},{},{},{},{}".format(
                    event, dosePerEvent[event], pathLength[event], meanKEperEvent[event], text))
                f.write('\n')

        with open(os.path.splitext(outputFilename)[0]+"_DSB.csv", "w") as f:
            f.write(
                "Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,SimGitHash,ClusteringGitHash\n")
            if filenamePhoton:
                f.write("{},{},{},{},{},{},{},{},{}\n".format(filename, energy, LET, NumIntersecting,
                        chromatinVolume, numBP, sugarPosFilename, gitHash, clusteringGitHash))
            else:
                f.write("{},{},{},{},{},{},{},{},{}\n".format(filename, energy, LET, len(
                    events), chromatinVolume, numBP, sugarPosFilename, gitHash, clusteringGitHash))
            f.write("EventNo,DoseGy,PathLength_nm,Direct_1_SBperDSBcluster,Direct_2_SBperDSBcluster,Direct_3_SBperDSBcluster,Direct_4_SBperDSBcluster,Direct_5_SBperDSBcluster,Direct_6_SBperDSBcluster,Direct_7_SBperDSBcluster,Direct_8_SBperDSBcluster,Direct_9_SBperDSBcluster,Direct_10_SBperDSBcluster,Indirect_1_SBperDSBcluster,Indirect_2_SBperDSBcluster,Indirect_3_SBperDSBcluster,Indirect_4_SBperDSBcluster,Indirect_5_SBperDSBcluster,Indirect_6_SBperDSBcluster,Indirect_7_SBperDSBcluster,Indirect_8_SBperDSBcluster,Indirect_9_SBperDSBcluster,Indirect_10_SBperDSBcluster,Hybrid_1_SBperDSBcluster,Hybrid_2_SBperDSBcluster,Hybrid_3_SBperDSBcluster,Hybrid_4_SBperDSBcluster,Hybrid_5_SBperDSBcluster,Hybrid_6_SBperDSBcluster,Hybrid_7_SBperDSBcluster,Hybrid_8_SBperDSBcluster,Hybrid_9_SBperDSBcluster,Hybrid_10_SBperDSBcluster,Mixed_1_SBperDSBcluster,Mixed_2_SBperDSBcluster,Mixed_3_SBperDSBcluster,Mixed_4_SBperDSBcluster,Mixed_5_SBperDSBcluster,Mixed_6_SBperDSBcluster,Mixed_7_SBperDSBcluster,Mixed_8_SBperDSBcluster,Mixed_9_SBperDSBcluster,Mixed_10_SBperDSBcluster,Total_1_SBperDSBcluster,Total_2_SBperDSBcluster,Total_3_SBperDSBcluster,Total_4_SBperDSBcluster,Total_5_SBperDSBcluster,Total_6_SBperDSBcluster,Total_7_SBperDSBcluster,Total_8_SBperDSBcluster,Total_9_SBperDSBcluster,Total_10_SBperDSBcluster\n")
            for event in events:
                if event in eventsWithClusteringResults:
                    text = str([clusterSize[eventsWithClusteringResults.index(event)][1:]]).strip(
                        "[").strip("]")
                else:
                    text = "0, "*49+"0"

                f.write("{},{},{},{}".format(
                    event, dosePerEvent[event], pathLength[event], text))
                f.write('\n')


if __name__ == "__main__":
    # Test clustering
    eventsListDirect = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2]
    copyListDirect = [6, 7, 208, 6000, 820, 822, 600,
                      620, 650, 670, 6, 7, 211, 295, 250, 8546]
    strandListDirect = [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0]
    eventsListIndirect = [0, 0, 0, 0, 0, 0, 1, 2, 2, 2]
    copyListIndirect = [6, 400, 410, 217, 823, 826, 200, 200, 300, 302]
    strandListIndirect = [1, 0, 1, 1, 1, 0, 1, 1, 0, 0]

    test = {
        # DoseGy,PathLength_nm,MeanKEMeV, TotalSBdirect,SSBdirect,cSSBdirect,DSBdirect,TotalSBindirect,SSBindirect,cSSBindirect,DSBindirect,TotalSBtotal,SSBtotal,cSSBtotal,DSBtotal
        0: ["N/A", "N/A", "N/A", 6, 2, 1, 1, 6, 2, 0, 2, 12, 1, 1, 3],
        1: ["N/A", "N/A", "N/A", 7, 5, 0, 1, 1, 1, 0, 0, 8, 6, 0, 1],
        2: ["N/A", "N/A", "N/A", 3, 3, 0, 0, 3, 1, 1, 0, 6, 3, 0, 1]
    }

    clusterSizetest = {
        0: ["N/A", "N/A", "N/A",
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  # direct
            0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  # indirect
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  # hybrid
            0, 0, 1, 1, 0, 0, 0, 0, 0, 0,  # mixed
            0, 1, 1, 1, 0, 0, 0, 0, 0, 0],  # total

        1: ["N/A", "N/A", "N/A",
            0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  # direct
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  # indirect
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  # hybrid
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  # mixed
            0, 1, 0, 0, 0, 0, 0, 0, 0, 0],  # total

        2: ["N/A", "N/A", "N/A",
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  # direct
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  # indirect
            0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  # hybrid
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  # mixed
            0, 0, 1, 0, 0, 0, 0, 0, 0, 0],  # total
    }

    tempResults = clustering(3, eventsListDirect, copyListDirect, strandListDirect,
                             eventsListIndirect, copyListIndirect, strandListIndirect)

    clusteringResults = tempResults[0]  # strand break number results
    clusterSize = tempResults[1]  # DSB cluster size, cluster of size 1-10
    for i in range(3):
        assert (clusteringResults[i][1:] == test[i][3:])
        assert (clusterSize[i][1:] == clusterSizetest[i][3:])

    print("tests passed")
