import os
import sys
import inspect
import numpy as np
from scipy.spatial import cKDTree
from typing import Union, Tuple
from collections import defaultdict
import uproot as up
from functools import partial

# To import clustering from build folder
currentdir = os.path.dirname(os.path.abspath(
    inspect.getfile(inspect.currentframe())))
builddir = currentdir+r"/build"
sys.path.insert(0, builddir)

up.default_library='np'

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

    return np.asarray(result)


def getResults(d0, idx0, d1, idx1):
    result = [-1, -1]

    Rdirect = 0.35  # nm

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

    return np.asarray(result)


def getIndex(point: list, T0: cKDTree, T1: cKDTree) -> list:

    d0, idx0 = T0.query(point, k=1)
    d1, idx1 = T1.query(point, k=1)

    d = min(d0, d1)
    idx = idx0 if d0 < d1 else idx1
    strand = 0 if d0 < d1 else 1

    if d > 1e-10:
        raise ValueError

    return [idx, strand]

def getIndex_v(d0, idx0, d1, idx1) -> list:

    d = min(d0, d1)
    idx = idx0 if d0 < d1 else idx1
    strand = 0 if d0 < d1 else 1

    if d > 1e-10:
        raise ValueError

    return np.asarray([idx, strand])

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

# def map_DNA_to_decay(eventInfo: dict):
#      # build map to decay event number
#     mapping = {}  # photon event ID to simulation event ID dict
#     step1_evts = eventInfo['PhotonEventID'].array(library='np')
#     step2_evts = eventInfo['EventNo'].array(library='np')
#     for evt_step1, evt_step2 in zip(step1_evts, step2_evts ): #loop on nentries                 
#         mapping[evt_step2] = evt_step1 #save photon event number into key simulation event number
#     return mapping

def map_decay_to_copyno(eventInfo):
    # build map to decay event number
    mapping_copyNo = {}  # photon event ID to simulation event ID dict
    step_1_evts = eventInfo['PhotonEventID'].array(library='np')
    copynos = eventInfo['copyNo'].array(library='np')
    for evt_step1, copyNo in zip(step_1_evts, copynos): #loop on nentries            
        mapping_copyNo[evt_step1] = copyNo #save photon event number into key simulation event number
    # end mapping 
    return mapping_copyNo

# def map_PID_to_DNA(eventInfo: dict, pid: int):
#     # # build map to decay event number
#     # mapping_PID = defaultdict(list)  # photon event ID to simulation event ID dict
#     # for PID, evt_step2 in zip(eventInfo['particleID'], eventInfo['EventNo']): #loop on nentries           
#     #     mapping_PID[PID].append(evt_step2) #save photon event number into key simulation event number
#     # # end mapping 
#     # return mapping_PID
#     return lambda pid: eventInfo['EventNo'][eventInfo['particleID']==pid]

def map_pid_idx(eventInfo_pid, pid: int):
    if pid==-1:
        return np.where(eventInfo_pid!=-1)[0]
    else:
        return np.where(eventInfo_pid==pid)[0]

def map_radius_copyno(n_boxes, boxes_per_R):
    ranges_radii = {}
    cnumbers = np.arange(0, n_boxes+boxes_per_R, boxes_per_R)
    n_files = n_boxes/boxes_per_R
    for r, (cmin, cmax) in enumerate(zip(cnumbers, cnumbers[1:])):
        ranges_radii[r] = range(cmin, cmax)
    assert(len(ranges_radii)==n_files)
    return ranges_radii

def calculateDose(eventEdep, chromatinVolume: float, evts_step1, evts_step2):

    print("start calculate dose")
    evts_step1_uniq, idx, _ = np.unique(evts_step1, return_counts=True,return_inverse=True)
    edep_J = eventEdep['Edep_J'].array(library='numpy')
    dosePerEvent = np.bincount(idx, edep_J[evts_step2]/(1000*chromatinVolume))
    meanKEperEvent = np.bincount(idx, edep_J[evts_step2])
    energy = np.sum(meanKEperEvent)/len(meanKEperEvent)
    del edep_J
    return energy, evts_step1_uniq, dosePerEvent, meanKEperEvent

def AccumulateEdep(direct, T0: cKDTree, T1: cKDTree, evts_step1, evts_step2, out_path):
    print("start accumulate Edep")
    
    eventNo = direct['EventNo'].array(library='numpy')
    meventNo = np.memmap(out_path+'_temp_mevt.dat', mode='w+', shape = eventNo.shape, dtype=eventNo.dtype)
    meventNo[:] = eventNo[:]
    del eventNo

    selected_direct = np.isin(meventNo, evts_step2)

    x = direct['x'].array(library='np')
    mx = np.memmap(out_path+'_temp_x.dat', mode='w+', shape = x.shape, dtype=x.dtype)
    mx[:] = x[:]
    del x
    y = direct['y'].array(library='np')
    my = np.memmap(out_path+'_temp_y.dat', mode='w+', shape = y.shape, dtype=y.dtype)
    my[:] = y[:]
    del y
    z = direct['z'].array(library='np')
    mz = np.memmap(out_path+'_temp_z.dat', mode='w+', shape = z.shape, dtype=z.dtype)
    mz[:] = z[:]
    del z
    
    xsel = mx[selected_direct]
    ysel = my[selected_direct]
    zsel = mz[selected_direct]

    xyz = np.vstack((xsel, ysel, zsel)).swapaxes(0, 1)

    eDep_eV = direct['eDep_eV'].array(library='numpy')
    medep = np.memmap(out_path+'_temp_edep.dat', mode='w+', shape = eDep_eV.shape, dtype=eDep_eV.dtype)
    medep[:] = eDep_eV[:]
    del eDep_eV

    edep = medep[selected_direct]
    if len(edep)==0:
        return [], []
    d0, idx0 = T0.query(xyz)
    d1, idx1 = T1.query(xyz)
    get_results_v = np.vectorize(getResults, otypes=[int], signature='(),(),(),()->(2)')
    results = get_results_v(d0, idx0, d1, idx1)
    filter_index = np.where(results[:,0]!=-1)[0]
    results_filtered = results[filter_index]
    evt2_filtered = meventNo[selected_direct][filter_index]
    edep_filtered = medep[selected_direct][filter_index]

    #change to evt_step1
    p2_to_p1 = dict(zip(evts_step2, evts_step1))
    map_p2_to_p1 = lambda eid: p2_to_p1[eid]
    mapv = np.vectorize(map_p2_to_p1)
    evt1_filtered = mapv(evt2_filtered)

    keys = np.vstack([evt1_filtered, results_filtered[:,0], results_filtered[:,1]]).swapaxes(0,1)
    
    # use bincount to accumulate edep
    _, idx, inv, _= np.unique(keys,axis=0, return_counts=True,return_inverse=True, return_index=True)
    cumEdep = np.bincount(inv, edep_filtered)
    keys_inv = keys[idx]

    # for decayEvt, simEvt in zip(evts_step1, evts_step2):
    #     x = direct['x'][direct['EventNo']==simEvt]
    #     y = direct['y'][direct['EventNo']==simEvt]
    #     z = direct['z'][direct['EventNo']==simEvt]
    #     for idx, (i, j ,k) in enumerate(zip(x, y, z)):
    #         result = checkPoint([i, j, k], T0, T1)
    #         if (result[0] != -1):
    #             key = (decayEvt, result[0], result[1])
    #             cumulatedEnergyDep.setdefault(key, 0)
    #             cumulatedEnergyDep[key] += direct['eDep_eV'][direct['EventNo']==simEvt][idx]
    # return cumulatedEnergyDep  
    # 
    os.remove(out_path+'_temp_mevt.dat')
    os.remove(out_path+'_temp_x.dat')
    os.remove(out_path+'_temp_y.dat')
    os.remove(out_path+'_temp_z.dat')
    os.remove(out_path+'_temp_edep.dat')
    return keys_inv, cumEdep

def calcDirectDamage(keys_edep, cumulatedEnergyDep, fEMinDamage: float, fEMaxDamage: float):
    print("start direct damage")
    
    if len(keys_edep)==0 or len(cumulatedEnergyDep)==0:
        return [], [], []
    check_edep = partial(IsEdepSufficient, fEMinDamage=fEMinDamage, fEMaxDamage=fEMaxDamage)
    check_edep_v = np.vectorize(check_edep)
    filtered_index = check_edep_v(cumulatedEnergyDep)
    eventsListDirect = keys_edep[filtered_index][:,0]
    copyListDirect = keys_edep[filtered_index][:,2]
    strandListDirect = keys_edep[filtered_index][:,1]
    # for key, edep in zip(keys_edep, cumulatedEnergyDep):
    #     if IsEdepSufficient(edep, fEMinDamage, fEMaxDamage):
    #         eventsListDirect.append(key[0])
    #         copyListDirect.append(key[2])
    #         strandListDirect.append(key[1])
    
    return eventsListDirect, copyListDirect, strandListDirect

def calcIndirectDamage(indirect: dict, probIndirect: float, T0: cKDTree, T1: cKDTree, evts_step1, evts_step2, out_path):


    eventsListIndirect = []
    copyListIndirect = []
    strandListIndirect = []
    
    # for decayEvt, simEvt in zip(evts_step1, evts_step2):
        
    #     if "copyNum" in indirect:
    #         pass
    #         # if np.random.rand() <= probIndirect:
    #         #     eventsListIndirect.append(decayEvt)
    #         #     copyListIndirect.append(c for c in input_tree_indirect['copyNum'][input_tree_indirect['EventNo']==simEvt])
    #     elif "DNAmolecule" not in indirect:
    #         pass
    #         # #this is wrong
    #         # if np.random.rand() <= probIndirect:
    #         #     eventsListIndirect.append(decayEvt)
    #         #     x = input_tree_indirect['x'][input_tree_indirect['EventNo']==simEvt]
    #         #     y = input_tree_indirect['y'][input_tree_indirect['EventNo']==simEvt] 
    #         #     z = input_tree_indirect['z'][input_tree_indirect['EventNo']==simEvt]
    #         #     for i, j, k in zip(x, y, z):
    #         #         c, s = getIndex([i, j, k], T0, T1)
    #         #         copyListIndirect.append(c)
    #         #         strandListIndirect.append(s) #this is wrong - lists do not have the same length - check
    #     else:
    eventNo = indirect['EventNo'].array(library='np')
    meventNo = np.memmap(out_path+'_temp_eventno_ind.dat', shape = eventNo.shape, mode='w+', dtype = eventNo.dtype)
    meventNo[:] = eventNo[:]
    del eventNo

    selected_indirect = np.isin(meventNo, evts_step2)
    evt2 = meventNo[selected_indirect]
    
    dnamolecule = indirect['DNAmolecule'].array(library='np')
    mdnamol = np.memmap(out_path+'_temp_mdnamol.dat', shape=dnamolecule.shape, dtype = dnamolecule.dtype, mode='w+')
    mdnamol[:] = dnamolecule[:]
    del dnamolecule

    selected_dnamol = mdnamol[selected_indirect]

    radicals = indirect['radical'].array(library='np')
    mradicals = np.memmap(out_path+'_temp_radicals.dat', shape = radicals.shape, dtype = radicals.dtype, mode='w+')
    mradicals[:] = radicals[:]
    del radicals
    selected_radicals = mradicals[selected_indirect]

    if len(selected_dnamol)==0 or len(selected_radicals)==0:
        return [], [], []
    selected_radical_idx = np.where((selected_radicals=='OH^0') & (selected_dnamol=='Deoxyribose^0'))[0]
    x = indirect['x'].array(library='np')[selected_indirect][selected_radical_idx]
    y = indirect['y'].array(library='np')[selected_indirect][selected_radical_idx]
    z = indirect['z'].array(library='np')[selected_indirect][selected_radical_idx]

    rand_arr = np.random.rand(selected_radical_idx.shape[0])
    rand_indx = rand_arr<=probIndirect

    points = np.vstack([x, y, z]).swapaxes(0,1)[rand_indx]
    d0, idx0 = T0.query(points)
    d1, idx1 = T1.query(points)
    get_index = np.vectorize(getIndex_v, otypes=[int], signature='(),(),(),()->(2)')
    results = get_index(d0, idx0, d1, idx1)            
    copyListIndirect = results[:,0]
    strandListIndirect = results[:,1]

    evt2_filtered = evt2[selected_radical_idx][rand_indx]
    p2_to_p1 = dict(zip(evts_step2, evts_step1))
    map_p2_to_p1 = lambda eid: p2_to_p1[eid]
    mapv = np.vectorize(map_p2_to_p1)
    eventsListIndirect = mapv(evt2_filtered)
        
        # for i, _ in enumerate(indirect['EventNo'][indirect['EventNo']==simEvt]):
        #     if ((indirect['DNAmolecule'][[indirect['EventNo']==simEvt]][i].startswith("D")) and (indirect['radical'][[indirect['EventNo']==simEvt]][i] == "OH^0")):
        #         if np.random.rand() <= probIndirect:
        #             eventsListIndirect.append(decayEvt)
        #             point = [indirect['x'][[indirect['EventNo']==simEvt]][i],
        #                     indirect['y'][[indirect['EventNo']==simEvt]][i], 
        #                     indirect['z'][[indirect['EventNo']==simEvt]][i]]
        #             c, s = getIndex(point, T0, T1)
        #             copyListIndirect.append(c)
        #             strandListIndirect.append(s)
    os.remove(out_path+'_temp_eventno_ind.dat')
    os.remove(out_path+'_temp_mdnamol.dat')
    os.remove(out_path+'_temp_radicals.dat')
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
    print("runClustering UP")
    if primaryParticle:
        if primaryParticle in particleMap.keys():
            primaryParticleID = particleMap[primaryParticle]
            print("primary: ", primaryParticle, primaryParticleID)
        elif primaryParticle == "all":
            primaryParticleID = -1
            print("primary: ", primaryParticle, primaryParticleID)
        else:
            raise ValueError(f"particle not in particleMap: {primaryParticle}")
    out_path = os.path.splitext(outputFilename)[0] if not primaryParticle else os.path.splitext(outputFilename)[0]+f"_{primaryParticle}"

    T0, T1 = readSugarFile(sugarFname)

    
    # Load data to analyse
    ufile = up.open(filename_DNA)
    
    if up.default_library!='np':
        up.default_library='np'

    info = ufile["output/Info"]
    eventInfo = ufile["ntuple/Events"]
    eventInfo_pid = eventInfo['particleID'].array(library='np')

    evt_id_pid = map_pid_idx(eventInfo_pid, primaryParticleID)
    
    evts_step1 = eventInfo['PhotonEventID'].array(library='np')[evt_id_pid]
    evts_step2 = eventInfo['EventNo'].array(library='np')[evt_id_pid]
    copyNo_pid = eventInfo['copyNo'].array(library='np')[evt_id_pid]

    pathLength={}
    chromatinVolume =  info['ChromatinVolume_m3'].array(library='np')[0]  # in m3
    numBP =  info['NumBasepairs'].array(library='np')[0]
    gitHash = info['GitHash'].array(library='np')[0]

    assert T0.n==numBP, "T0 data does not match numGBP"


    if filenamePhoton:
        LET = "N/A"
    else:
        if hasattr(info,"MeanLET"):
            LET = float(info['MeanLET'])
        else:
            LET = "N/A"

  
    #cumulatedEnergyDep = {}
    
    eventEdep = ufile["ntuple/EventEdep"]
    
    energy, evt_step1_unique_dose, dosePerEvent, meanKEperEvent = calculateDose(eventEdep, chromatinVolume, evts_step1, evts_step2)   
    #print(evt_step1_unique, dosePerEvent)
    print("calculateDose Done")

    print("loading Direct tree")
    direct = ufile["ntuple/Direct"]
    
    keys_edep, cumulatedEnergyDep = AccumulateEdep(direct, T0, T1, evts_step1, evts_step2, out_path)
    

    print("AccumulateEdep done")
 
    
    
    # // Read out indirect damage
    eventsListDirect, copyListDirect, strandListDirect = calcDirectDamage(keys_edep, cumulatedEnergyDep, fEMinDamage, fEMaxDamage)
    
    indirect = ufile["ntuple/Indirect"]
    eventsListIndirect, copyListIndirect, strandListIndirect = calcIndirectDamage(indirect, probIndirect, T0, T1, evts_step1, evts_step2, out_path)
  

    
    
    numEvt = list(set(np.concatenate([eventsListDirect, eventsListIndirect])))
    #(len(eventsListDirect), len(eventsListIndirect), len(copyListDirect), len(strandListDirect), len(numEvt))
    # clustering
    tempResults = clustering(numEvt, eventsListDirect, copyListDirect,
                             strandListDirect, eventsListIndirect, copyListIndirect, strandListIndirect)

    clusteringResults = np.asarray(tempResults[0], dtype=int)  # strand break number results
    clusterSize = np.asarray(tempResults[1])  # DSB cluster size, cluster of size 1-10
    events = evt_step1_unique_dose

    #events = list(dosePerEvent.keys())
    if len(clusteringResults):
        eventsWithClusteringResults = clusteringResults[:,0]
    else:
        raise ValueError("no clustering results. Done")
        
    _, idx1, idx2 = np.intersect1d(events, eventsWithClusteringResults, return_indices=True)
    #complement = events[np.isin(events, eventsWithClusteringResults, invert=True)]

    text_SB = np.zeros(shape = (events.shape[0], 12), dtype=int)
    text_SB[idx1] = clusteringResults[:, 1:][idx2]
    
    text_DSB = np.zeros(shape=(events.shape[0], 50), dtype=int)
    text_DSB[idx1] = clusterSize[:,1:][idx2]
    clusteringGitHash = readClusteringGitHash(builddir)

    print("Finished: {}".format(filename_DNA))

    # selecting copynumbers to save results in separate files
    mapping_copyNo = map_decay_to_copyno(eventInfo)
    ranges_radii = map_radius_copyno(n_boxes, boxes_per_R)
    
    if separate_r:
        print("separating output files per radius")
        for r in ranges_radii:
            outfile = out_path+f"_{r}.csv" 
            with open(outfile, "w") as f:
                
                f.write(
                    "Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,GitHash,ClusteringGitHash\n")
                # if filenamePhoton:
                #     f.write("{},{},{},{},{},{},{},{},{}\n".format(filename_DNA, energy, LET, NumIntersecting,
                #             chromatinVolume, numBP, sugarFname, gitHash, clusteringGitHash))
                # else:
                
                f.write("{},{},{},{},{},{},{},{},{}\n".format(filename_DNA, energy, LET, len(
                    pathLength.keys()), chromatinVolume, numBP, sugarFname, gitHash, clusteringGitHash))
                f.write("EventNo,DoseGy,PathLength_nm,MeanKEMeV,TotalSBdirect,SSBdirect,cSSBdirect,DSBdirect,TotalSBindirect,SSBindirect,cSSBindirect,DSBindirect,TotalSBtotal,SSBtotal,cSSBtotal,DSBtotal\n")
                for i, event in enumerate(events):
                    if mapping_copyNo[event] in ranges_radii[r]:
                        #if event in eventsWithClusteringResults:
                        #    text = str([clusteringResults[list(eventsWithClusteringResults).index(
                        #        event)][1:]]).strip("[").strip("]")
                        #else:
                        #    text = "0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0"
                        f.write("{},{},{},{},{}".format(
                            event, 
                            dosePerEvent[i], 
                            "N/A", 
                            meanKEperEvent[i], 
                            str(list(text_SB[i])).strip("[").strip("]")))
                        f.write('\n')
            
                
            outfile_DSB = out_path+f"_DSB_{r}.csv"
            with open(outfile_DSB, "w") as f:
                f.write(
                    "Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,SimGitHash,ClusteringGitHash\n")
                # if filenamePhoton:
                #     f.write("{},{},{},{},{},{},{},{},{}\n".format(filename_DNA, energy, LET, NumIntersecting,
                #             chromatinVolume, numBP, sugarFname, gitHash, clusteringGitHash))
                # else:
                f.write("{},{},{},{},{},{},{},{},{}\n".format(filename_DNA, energy, LET, len(
                        events), chromatinVolume, numBP, sugarFname, gitHash, clusteringGitHash))
                f.write("EventNo,DoseGy,PathLength_nm,Direct_1_SBperDSBcluster,Direct_2_SBperDSBcluster,Direct_3_SBperDSBcluster,Direct_4_SBperDSBcluster,Direct_5_SBperDSBcluster,Direct_6_SBperDSBcluster,Direct_7_SBperDSBcluster,Direct_8_SBperDSBcluster,Direct_9_SBperDSBcluster,Direct_10_SBperDSBcluster,Indirect_1_SBperDSBcluster,Indirect_2_SBperDSBcluster,Indirect_3_SBperDSBcluster,Indirect_4_SBperDSBcluster,Indirect_5_SBperDSBcluster,Indirect_6_SBperDSBcluster,Indirect_7_SBperDSBcluster,Indirect_8_SBperDSBcluster,Indirect_9_SBperDSBcluster,Indirect_10_SBperDSBcluster,Hybrid_1_SBperDSBcluster,Hybrid_2_SBperDSBcluster,Hybrid_3_SBperDSBcluster,Hybrid_4_SBperDSBcluster,Hybrid_5_SBperDSBcluster,Hybrid_6_SBperDSBcluster,Hybrid_7_SBperDSBcluster,Hybrid_8_SBperDSBcluster,Hybrid_9_SBperDSBcluster,Hybrid_10_SBperDSBcluster,Mixed_1_SBperDSBcluster,Mixed_2_SBperDSBcluster,Mixed_3_SBperDSBcluster,Mixed_4_SBperDSBcluster,Mixed_5_SBperDSBcluster,Mixed_6_SBperDSBcluster,Mixed_7_SBperDSBcluster,Mixed_8_SBperDSBcluster,Mixed_9_SBperDSBcluster,Mixed_10_SBperDSBcluster,Total_1_SBperDSBcluster,Total_2_SBperDSBcluster,Total_3_SBperDSBcluster,Total_4_SBperDSBcluster,Total_5_SBperDSBcluster,Total_6_SBperDSBcluster,Total_7_SBperDSBcluster,Total_8_SBperDSBcluster,Total_9_SBperDSBcluster,Total_10_SBperDSBcluster\n")
                for i, event in enumerate(events):
                    if mapping_copyNo[event] in ranges_radii[r]:
                        #if event in eventsWithClusteringResults:
                        #    text = str([clusterSize[eventsWithClusteringResults.index(event)][1:]]).strip(
                        #        "[").strip("]")
                        #else:
                        #    text = "0, "*49+"0"

                        f.write("{},{},{},{}".format(
                            event, dosePerEvent[i],"N/A", str(list(text_DSB[i])).strip("[").strip("]")))
                        f.write('\n')
    else:
        outfile_tot = outputFilename
        with open(outfile_tot, "w") as f:
            f.write(
                "Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,GitHash,ClusteringGitHash\n")
            # if filenamePhoton:
            #     f.write("{},{},{},{},{},{},{},{},{}\n".format(filename_DNA, energy, LET, NumIntersecting,
            #             chromatinVolume, numBP, sugarFname, gitHash, clusteringGitHash))
            # else:
            f.write("{},{},{},{},{},{},{},{},{}\n".format(filename_DNA, energy, LET, len(
                pathLength.keys()), chromatinVolume, numBP, sugarFname, gitHash, clusteringGitHash))
            f.write("EventNo,DoseGy,PathLength_nm,MeanKEMeV,TotalSBdirect,SSBdirect,cSSBdirect,DSBdirect,TotalSBindirect,SSBindirect,cSSBindirect,DSBindirect,TotalSBtotal,SSBtotal,cSSBtotal,DSBtotal\n")
            for i, event in enumerate(events):
                # if event in eventsWithClusteringResults:
                #     text = str([clusteringResults[eventsWithClusteringResults.index(
                #         event)][1:]]).strip("[").strip("]")
                # else:
                #     text = "0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0"
                f.write("{},{},{},{},{}".format(
                    event, dosePerEvent[i], "N/A", meanKEperEvent[i], str(list(text_SB[i])).strip("[").strip("]")))
                f.write('\n')
        outfile_DSB_tot = out_path+"_DSB.csv"
        with open(outfile_DSB_tot, "w") as f:
            f.write(
                "Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,SimGitHash,ClusteringGitHash\n")
            # if filenamePhoton:
            #     f.write("{},{},{},{},{},{},{},{},{}\n".format(filename_DNA, energy, LET, NumIntersecting,
            #             chromatinVolume, numBP, sugarFname, gitHash, clusteringGitHash))
            # else:
            f.write("{},{},{},{},{},{},{},{},{}\n".format(filename_DNA, energy, LET, len(
                    events), chromatinVolume, numBP, sugarFname, gitHash, clusteringGitHash))
            f.write("EventNo,DoseGy,PathLength_nm,Direct_1_SBperDSBcluster,Direct_2_SBperDSBcluster,Direct_3_SBperDSBcluster,Direct_4_SBperDSBcluster,Direct_5_SBperDSBcluster,Direct_6_SBperDSBcluster,Direct_7_SBperDSBcluster,Direct_8_SBperDSBcluster,Direct_9_SBperDSBcluster,Direct_10_SBperDSBcluster,Indirect_1_SBperDSBcluster,Indirect_2_SBperDSBcluster,Indirect_3_SBperDSBcluster,Indirect_4_SBperDSBcluster,Indirect_5_SBperDSBcluster,Indirect_6_SBperDSBcluster,Indirect_7_SBperDSBcluster,Indirect_8_SBperDSBcluster,Indirect_9_SBperDSBcluster,Indirect_10_SBperDSBcluster,Hybrid_1_SBperDSBcluster,Hybrid_2_SBperDSBcluster,Hybrid_3_SBperDSBcluster,Hybrid_4_SBperDSBcluster,Hybrid_5_SBperDSBcluster,Hybrid_6_SBperDSBcluster,Hybrid_7_SBperDSBcluster,Hybrid_8_SBperDSBcluster,Hybrid_9_SBperDSBcluster,Hybrid_10_SBperDSBcluster,Mixed_1_SBperDSBcluster,Mixed_2_SBperDSBcluster,Mixed_3_SBperDSBcluster,Mixed_4_SBperDSBcluster,Mixed_5_SBperDSBcluster,Mixed_6_SBperDSBcluster,Mixed_7_SBperDSBcluster,Mixed_8_SBperDSBcluster,Mixed_9_SBperDSBcluster,Mixed_10_SBperDSBcluster,Total_1_SBperDSBcluster,Total_2_SBperDSBcluster,Total_3_SBperDSBcluster,Total_4_SBperDSBcluster,Total_5_SBperDSBcluster,Total_6_SBperDSBcluster,Total_7_SBperDSBcluster,Total_8_SBperDSBcluster,Total_9_SBperDSBcluster,Total_10_SBperDSBcluster\n")
            for i, event in enumerate(events):
                # if event in eventsWithClusteringResults:
                #     text = str([clusterSize[eventsWithClusteringResults.index(event)][1:]]).strip(
                #         "[").strip("]")
                # else:
                #     text = "0, "*49+"0"

                f.write("{},{},{},{}".format(
                    event, dosePerEvent[i], "N/A", str(list(text_DSB[i])).strip("[").strip("]")))
                f.write('\n')


if __name__ == "__main__":
    pass
