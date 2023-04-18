import os
import sys
import inspect
from typing import Union
import numpy as np
from scipy.spatial import cKDTree
from ROOT import TFile, TTree
from typing import Union, Tuple
from collections import defaultdict

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


def readSugarFile(sugarFname: str) -> Tuple[cKDTree, cKDTree]:
    """ Function to read sugar geometry files
    Args: 
        sugarFname (str): Path to sugar geometry file

    Returns: 
        Tuple[cKDTree, cKDTree]: KDTrees containing the positions of the sugars on stand 0 and strand 1

    """
    # Read sugar csv file
    with open(sugarFname, "r") as f:
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

    return T0, T1

def map_DNA_to_decay(eventInfo: TTree):
     # build map to decay event number
    nentriesMapping = eventInfo.GetEntries() #get number of events for mapping
    mapping = {}  # photon event ID to simulation event ID dict
    for irow in range(0, nentriesMapping): #loop on nentries
        eventInfo.GetEntry(irow)            
        mapping[eventInfo.EventNo] = eventInfo.PhotonEventID #save photon event number into key simulation event number
    # end mapping 
    return mapping

def map_decay_to_copyno(eventInfo: TTree):
    # build map to decay event number
    nentriesMapping = eventInfo.GetEntries() #get number of events for mapping
    mapping_copyNo = {}  # photon event ID to simulation event ID dict
    for irow in range(0, nentriesMapping): #loop on nentries
        eventInfo.GetEntry(irow)            
        mapping_copyNo[eventInfo.PhotonEventID] = eventInfo.copyNo #save photon event number into key simulation event number
    # end mapping 
    return mapping_copyNo

def map_PID_to_DNA(eventInfo: TTree):
    # build map to decay event number
    nentriesMapping = eventInfo.GetEntries() #get number of events for mapping
    mapping_PID = defaultdict(list)  # photon event ID to simulation event ID dict
    for irow in range(0, nentriesMapping): #loop on nentries
        eventInfo.GetEntry(irow)            
        mapping_PID[eventInfo.particleID].append(eventInfo.EventNo) #save photon event number into key simulation event number
    # end mapping 
    return mapping_PID


def map_radius_copyno(n_boxes, boxes_per_R):
    ranges_radii = {}
    cnumbers = np.arange(0, n_boxes+boxes_per_R, boxes_per_R)
    n_files = n_boxes/boxes_per_R
    for r, (cmin, cmax) in enumerate(zip(cnumbers, cnumbers[1:])):
        ranges_radii[r] = range(cmin, cmax)
    assert(len(ranges_radii)==n_files)
    return ranges_radii

def calculateDose(tEdep: TTree, pathLength: dict, dosePerEvent: dict, meanKEperEvent: dict, chromatinVolume: float, mapping_evts: dict, mapping_PID: dict, primaryParticle: int):
    print("start dose calculation")
    sim_evts = mapping_PID[primaryParticle]
    entryEdepNumber = tEdep.GetEntries()
    print(entryEdepNumber)

    #numEvt = 0
    for i in range(entryEdepNumber):
        
        tEdep.GetEntry(i)
        simEvt = tEdep.EventNo
        if simEvt not in sim_evts:
            continue
        decayEvt = mapping_evts[simEvt]
        #print(simEvt, decayEvt)
        #print(decayEvt, simEvt, primaryParticle)
        if tEdep.Edep_J > 0:
            dosePerEvent.setdefault(decayEvt, 0)
            dosePerEvent[decayEvt] += tEdep.Edep_J / (1000 * chromatinVolume)
            #print(f"{decayEvt}: {dosePerEvent[decayEvt]}")
            if hasattr(tEdep, "PathLengthChromatin"):  # older results do not have path length
                pathLength[decayEvt] = tEdep.PathLengthChromatin
            else:
                pathLength[decayEvt] = "N/A"

            if hasattr(tEdep, "PrimaryKEEntrance"):  # older results do not have path length
                meanKEperEvent[decayEvt] = (
                    tEdep.PrimaryKEEntrance+tEdep.PrimaryKEExit)/2
            else:
                meanKEperEvent[decayEvt] = "N/A"
    
        else:
            if tEdep.Edep_J>0:
                print(f"evt {simEvt} not in mapping but Edep>0")

    if meanKEperEvent and list(meanKEperEvent.values())[0] != "N/A":
        energy = sum(meanKEperEvent.values())/len(meanKEperEvent.values())
    else:
        energy = "N/A"

    return energy, dosePerEvent, meanKEperEvent, pathLength

def AccumulateEdep(input_tree: TTree, cumulatedEnergyDep: dict, T0: cKDTree, T1: cKDTree, mapping_evts: dict, mapping_PID: dict, primaryParticle: int):
    sim_evts = mapping_PID[primaryParticle]
    print("start AccumulateEdep")
    nentries = input_tree.GetEntries()
    
    for irow in range(nentries):
        input_tree.GetEntry(irow)
        simEvt = input_tree.EventNo
        if simEvt not in sim_evts:
            continue
            #print([input_tree.x, input_tree.y, input_tree.z])
        decayEvt = mapping_evts[simEvt]
        result = checkPoint([input_tree.x, input_tree.y, input_tree.z], T0, T1)
        if (result[0] != -1):
            key = (decayEvt, result[0], result[1])
            cumulatedEnergyDep.setdefault(key, 0)
            cumulatedEnergyDep[key] += input_tree.eDep_eV

    return cumulatedEnergyDep

def calcDirectDamage(cumulatedEnergyDep: dict, fEMinDamage: float, fEMaxDamage: float):
    eventsListDirect = []
    copyListDirect = []
    strandListDirect = []
    print("start DirectDamage")
    for key in cumulatedEnergyDep:
        if IsEdepSufficient(cumulatedEnergyDep[key], fEMinDamage, fEMaxDamage):
            eventsListDirect.append(key[0])
            copyListDirect.append(key[2])
            strandListDirect.append(key[1])

    return eventsListDirect, copyListDirect, strandListDirect

def calcIndirectDamage(input_tree_indirect: TTree, probIndirect: float, T0: cKDTree, T1: cKDTree, mapping_evts: dict, mapping_PID: dict, primaryParticle: int):

    print("start IndirectDamage")
    eventsListIndirect = []
    copyListIndirect = []
    strandListIndirect = []

    nentries_indirect = input_tree_indirect.GetEntries()
    sim_evts = mapping_PID[primaryParticle]
    for irow in range(nentries_indirect):
        input_tree_indirect.GetEntry(irow)
        simEvt = input_tree_indirect.EventNo
        if simEvt not in sim_evts:
            continue
        decayEvt = mapping_evts[simEvt]

        if hasattr(input_tree_indirect, "copyNum"):
            if np.random.rand() <= probIndirect:
                eventsListIndirect.append(decayEvt)
                copyListIndirect.append(input_tree_indirect.copyNum)
        elif not hasattr(input_tree_indirect, "DNAmolecule"):
            if np.random.rand() <= probIndirect:
                eventsListIndirect.append(decayEvt)
                point = [input_tree_indirect.x,
                        input_tree_indirect.y, input_tree_indirect.z]
                c, s = getIndex(point, T0, T1)
                copyListIndirect.append(c)
                strandListIndirect.append(s)
        else:
            if ((input_tree_indirect.DNAmolecule[0] == "D") and (input_tree_indirect.radical == "OH^0")):
                if np.random.rand() <= probIndirect:
                    eventsListIndirect.append(decayEvt)
                    point = [input_tree_indirect.x,
                            input_tree_indirect.y, input_tree_indirect.z]
                    c, s = getIndex(point, T0, T1)
                    copyListIndirect.append(c)
                    strandListIndirect.append(s)
    
    return eventsListIndirect, copyListIndirect, strandListIndirect


def runClustering(filename_DNA: str, outputFilename: str, fEMinDamage: float, fEMaxDamage: float, probIndirect: float, sugarFname: str, primaryParticle: str = None, filenamePhoton: Union[str, bool] = False, separate_r=False, n_boxes = 3200, boxes_per_R = 800):
    """ Run clustering on the DNA root file and output DNA damage
    
    Args:
        filename_DNA (str): DNA root file
        outputFilename (str): csv output file
        fEMinDamage (float): minimimum energy deposit for a strand break
        fEMaxDamage (float): maximum energy above which a strand break always occurs
        probIndirect (float): probability that OH+sugar reaction leads to strand damage
        sugarFname (str): text file containing the sugar molecule positions, must match the one used for simulation
        #simulationType (str, optional): Type of simulation standalone, decay, photon. Defaults to "standalone".
        filenamePhoton (Union[str,bool], optional): If photon, path to root file from photon simulation. Defaults to False.
        #continuous (bool, optional): Continuous DNA structure from fractal DNA - True, original DNA structure with discontinuities - False. Defaults to True.
        #part1_CopyNum (Union[int,bool], optional): If decay, which copy number to use for clustering. Defaults to False.
        primaryParticle (Union[str,bool], optional): If decay, particle to use for clustering alpha or e-, if False all. Defaults to False.
        separate_r (bool, optional): save DNA damage all together or separating boxes rings
        n_boxes (int): total number of boxes
        boxes_per_R (int): total number of boxes in a ring
    
    """
    particleMap = { "alpha": 1 ,
                    "gamma": 2,
                    "e-": 3, 
                    "nu_e": 4,
                      "At211": 5,
                      "Po211": 6, 
                      "Bi207": 7, 
                      "Pb207": 8
    }
    print("runClustering REFACTOR")
    if primaryParticle:
        if primaryParticle in particleMap.keys():
            primaryParticleID = particleMap[primaryParticle]
            print("primary: ", primaryParticle, primaryParticleID)
        else:
            raise ValueError(f"particle not in particleMap: {primaryParticle}")

    T0, T1 = readSugarFile(sugarFname)

    # Load data to analyse
    fFile = TFile.Open(filename_DNA)
    fDirectory = fFile.Get("ntuple")
    tInfo = fDirectory.Get("Info")
    tInfo.GetEntry(0)
    tEdep = fDirectory.Get("EventEdep")
    
    tDirect = fDirectory.Get("Direct")
    tIndirect = fDirectory.Get("Indirect")
    eventInfo = fDirectory.Get("Events")

    chromatinVolume = tInfo.ChromatinVolume_m3  # in m3
    numBP = tInfo.NumBasepairs
    gitHash = tInfo.GitHash
    assert T0.n==numBP, "T0 data does not match numGBP"

    if filenamePhoton:
        LET = "N/A"
    else:
        if hasattr(tInfo,"MeanLET"):
            LET = tInfo.MeanLET
        else:
            LET = "N/A"

    pathLength = {}
    dosePerEvent = {}
    meanKEperEvent = {}
    cumulatedEnergyDep = {}

    mapping_evts = map_DNA_to_decay(eventInfo)
    mapping_PID = map_PID_to_DNA(eventInfo)

    # if primaryParticleID:
    energy, dosePerEvent, _, _ = calculateDose(tEdep, pathLength, dosePerEvent, meanKEperEvent, chromatinVolume, mapping_evts, mapping_PID, primaryParticleID)   
    print("done")
    cumulatedEnergyDep = AccumulateEdep(tDirect, cumulatedEnergyDep, T0, T1, mapping_evts, mapping_PID, primaryParticleID)
    print("done")
    #print(cumulatedEnergyDep)
    # else:
    #     energy, dosePerEvent, _, _ = calculateDose(tEdep, pathLength, dosePerEvent, meanKEperEvent, chromatinVolume, mapping_evts)   
    #     cumulatedEnergyDep = AccumulateEdep(input_tree, cumulatedEnergyDep, T0, T1, mapping_evts)
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

    #  Add energy deposition from photon simulation
    # if filenamePhoton:
    #     input_tree_photon = fDirectoryPhoton.Get("Direct")
    #     nentriesPhoton = input_tree_photon.GetEntries()

    #     cumulatedEnergyDepPhoton = {}

    #     for irow in range(nentriesPhoton):
    #         input_tree_photon.GetEntry(irow)

    #         result = checkPoint(
    #             [input_tree_photon.x, input_tree_photon.y, input_tree_photon.z], T0, T1)
    #         if (result[0] != -1):
    #             key = (input_tree_photon.EventNo, result[0], result[1])
    #             if key in cumulatedEnergyDepPhoton:
    #                 cumulatedEnergyDepPhoton[key] += input_tree_photon.eDep_eV
    #             else:
    #                 cumulatedEnergyDepPhoton[key] = input_tree_photon.eDep_eV

    #     for key in cumulatedEnergyDep:
    #         if ((mapping[key[0]], key[1], key[2]) in cumulatedEnergyDepPhoton):
    #             cumulatedEnergyDep[key] += cumulatedEnergyDepPhoton[mapping[key[0]], key[1], key[2]]
    #     #print(f"cumulatedEnergyDep: {cumulatedEnergyDep}")
    #     # Add photon dose
    #     for key in dosePerEvent:
    #         if mapping[key] in dosePerEventPhoton:
    #             dosePerEvent[key] += dosePerEventPhoton[mapping[key]]

    
    # // Read out indirect damage
    eventsListDirect, copyListDirect, strandListDirect = calcDirectDamage(cumulatedEnergyDep, fEMinDamage, fEMaxDamage)
    print("done")
    eventsListIndirect, copyListIndirect, strandListIndirect = calcIndirectDamage(tIndirect, probIndirect, T0, T1, mapping_evts, mapping_PID, primaryParticleID)
    print("done")
    numEvt = list(set(eventsListDirect+eventsListIndirect))
    
    # clustering
    tempResults = clustering(numEvt, eventsListDirect, copyListDirect,
                             strandListDirect, eventsListIndirect, copyListIndirect, strandListIndirect)

    clusteringResults = tempResults[0]  # strand break number results
    clusterSize = tempResults[1]  # DSB cluster size, cluster of size 1-10

    events = list(dosePerEvent.keys())
    eventsWithClusteringResults = [a[0] for a in clusteringResults]
    
    #print(eventsWithClusteringResults)
    
    clusteringGitHash = readClusteringGitHash(builddir)

    print("Finished: {}".format(filename_DNA))

    # selecting copynumbers to save results in separate files
    mapping_copyNo = map_decay_to_copyno(eventInfo)
    ranges_radii = map_radius_copyno(n_boxes, boxes_per_R)
    
    if separate_r:
        print("separating output files per radius")
        for r in ranges_radii:
            outfile = os.path.splitext(outputFilename)[0]+f"_{r}.csv" if not primaryParticle else os.path.splitext(outputFilename)[0]+f"_{primaryParticle}_{r}.csv"
            with open(outfile, "w") as f:
                f.write(
                    "Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,GitHash,ClusteringGitHash\n")
                if filenamePhoton:
                    f.write("{},{},{},{},{},{},{},{},{}\n".format(filename_DNA, energy, LET, NumIntersecting,
                            chromatinVolume, numBP, sugarFname, gitHash, clusteringGitHash))
                else:
                    f.write("{},{},{},{},{},{},{},{},{}\n".format(filename_DNA, energy, LET, len(
                        pathLength.keys()), chromatinVolume, numBP, sugarFname, gitHash, clusteringGitHash))
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
            outfile_DSB = os.path.splitext(outputFilename)[0]+f"_DSB_{r}.csv" if not primaryParticle else os.path.splitext(outputFilename)[0]+f"_{primaryParticle}_DSB_{r}.csv"
            with open(outfile_DSB, "w") as f:
                f.write(
                    "Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,SimGitHash,ClusteringGitHash\n")
                if filenamePhoton:
                    f.write("{},{},{},{},{},{},{},{},{}\n".format(filename_DNA, energy, LET, NumIntersecting,
                            chromatinVolume, numBP, sugarFname, gitHash, clusteringGitHash))
                else:
                    f.write("{},{},{},{},{},{},{},{},{}\n".format(filename_DNA, energy, LET, len(
                        events), chromatinVolume, numBP, sugarFname, gitHash, clusteringGitHash))
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
        outfile_tot = outputFilename if not primaryParticle else os.path.splitext(outputFilename)[0]+f"_{primaryParticle}.csv"
        with open(outfile_tot, "w") as f:
            f.write(
                "Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,GitHash,ClusteringGitHash\n")
            if filenamePhoton:
                f.write("{},{},{},{},{},{},{},{},{}\n".format(filename_DNA, energy, LET, NumIntersecting,
                        chromatinVolume, numBP, sugarFname, gitHash, clusteringGitHash))
            else:
                f.write("{},{},{},{},{},{},{},{},{}\n".format(filename_DNA, energy, LET, len(
                    pathLength.keys()), chromatinVolume, numBP, sugarFname, gitHash, clusteringGitHash))
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
        outfile_DSB_tot = os.path.splitext(outputFilename)[0]+"_DSB.csv" if not primaryParticle else os.path.splitext(outputFilename)[0]+f"_{primaryParticle}_DSB.csv"
        with open(outfile_DSB_tot, "w") as f:
            f.write(
                "Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,SimGitHash,ClusteringGitHash\n")
            if filenamePhoton:
                f.write("{},{},{},{},{},{},{},{},{}\n".format(filename_DNA, energy, LET, NumIntersecting,
                        chromatinVolume, numBP, sugarFname, gitHash, clusteringGitHash))
            else:
                f.write("{},{},{},{},{},{},{},{},{}\n".format(filename_DNA, energy, LET, len(
                    events), chromatinVolume, numBP, sugarFname, gitHash, clusteringGitHash))
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
    pass
