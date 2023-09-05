import os
import sys
import inspect
import numpy as np
from scipy.spatial import cKDTree
from typing import Union, Tuple
from collections import defaultdict
import uproot as up
from functools import partial
from scipy import constants
import pandas as pd

# To import clustering from build folder
currentdir = os.path.dirname(os.path.abspath(
    inspect.getfile(inspect.currentframe())))
builddir = currentdir+r"/build"
sys.path.insert(0, builddir)

up.default_library = 'np'

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



def calculateDose(eventEdep, chromatinVolume: float):

    print("start calculate dose")

    edep = eventEdep.arrays(library='np')
    edep_df = pd.DataFrame.from_dict(edep)

    ke_dose = edep_df

    ke_dose = ke_dose.groupby(['step1_eventID', 'step1_copyNo', 'step1_PID', 'step1_primaryID'], as_index=False)[
        ['edep_J', 'edep_MeV']].sum()
    
    ke_dose['dose'] = ke_dose['edep_J']/(1000*chromatinVolume)
    # ke_dose.loc[(ke_dose['step1_PID'] == 3) & (
    #     ke_dose['step1_primaryID'] != 1), 'step1_PID'] = 1

    mean_energy = ke_dose.groupby('step1_eventID')['edep_MeV'].mean().mean()

    # return mean energy per event, unique pairs, dose, meanke
    return mean_energy, ke_dose


def AccumulateEdep(direct, T0: cKDTree, T1: cKDTree, keys_df, out_path):
    print("start accumulate Edep")
    # load event number (step2) from direct tree
    eventNo = direct['step2_eventID'].array(library='numpy')
    meventNo = np.memmap(out_path+'_temp_mevt.dat', mode='w+',
                         shape=eventNo.shape, dtype='int64')
    meventNo[:] = eventNo[:]
    del eventNo
   
    #get x, y, z
    x = direct['x'].array(library='np')
    mx = np.memmap(out_path+'_temp_x.dat', mode='w+',
                   shape=x.shape, dtype=x.dtype)
    mx[:] = x[:]
    del x
    y = direct['y'].array(library='np')
    my = np.memmap(out_path+'_temp_y.dat', mode='w+',
                   shape=y.shape, dtype=y.dtype)
    my[:] = y[:]
    del y
    z = direct['z'].array(library='np')
    mz = np.memmap(out_path+'_temp_z.dat', mode='w+',
                   shape=z.shape, dtype=z.dtype)
    mz[:] = z[:]
    del z

    xyz = np.vstack((mx, my, mz)).swapaxes(0, 1)
    #get energy dep.
    eDep_eV = direct['eDep_eV'].array(library='numpy')
    medep = np.memmap(out_path+'_temp_edep.dat', mode='w+',
                      shape=eDep_eV.shape, dtype=eDep_eV.dtype)
    medep[:] = eDep_eV[:]
    del eDep_eV
    
    
    if len(medep) == 0:
        return [], []
    # get strand no. and copy
    d0, idx0 = T0.query(xyz)
    d1, idx1 = T1.query(xyz)

    get_results_v = np.vectorize(
        getResults, otypes=[int], signature='(),(),(),()->(2)')
    results = get_results_v(d0, idx0, d1, idx1)
    #get points with result != -1
    filter_index = np.where(results[:, 0] != -1)[0]
    results_filtered = results[filter_index]
    #filter events (step2) and edep
    evt2_filtered = meventNo[filter_index]
    edep_filtered = medep[filter_index]
    #change to evt_step1
    keys_map = keys_df.loc[keys_df['step2_eventID'].isin(
        evt2_filtered)].set_index('step2_eventID')

    #keys_map_ind = keys_map.set_index('step2')

    evt1_cp_filtered = keys_map.loc[np.array(evt2_filtered)]
    #evt1_filtered = mapv(evt2_filtered) #evt step1 corresponding to filtered step2 evts (with result != -1)

    evt1_cp_filtered['strand'] = results_filtered[:, 0]
    # evt ID (1), copyNo, strand, copy
    evt1_cp_filtered['copy'] = results_filtered[:, 1]
    evt1_cp_filtered['edep'] = edep_filtered

    cumEdep = evt1_cp_filtered.groupby(
        ['step1_eventID', 'step1_copyNo', 'strand', 'copy', 'step1_particleID', 'step1_primaryID'], as_index=False)['edep'].sum()
    # print(cumEdep[cumEdep['step1_eventID'] == 4].groupby(['step1_particleID', 'step1_copyNo'])['edep'].sum())
    # cumEdep.loc[(cumEdep['step1_particleID'] == 3) & (
    #     cumEdep['step1_primaryID'] != 1), 'step1_particleID'] = 1
    os.remove(out_path+'_temp_mevt.dat')
    os.remove(out_path+'_temp_x.dat')
    os.remove(out_path+'_temp_y.dat')
    os.remove(out_path+'_temp_z.dat')
    os.remove(out_path+'_temp_edep.dat')
    return cumEdep


def calcDirectDamage(cumulatedEnergyDep, fEMinDamage: float, fEMaxDamage: float):
    print("start direct damage")
    # if len(cumulatedEnergyDep[0])==0:
    #     return pd.DataFrame(columns = ['None'])

    check_edep = partial(
        IsEdepSufficient, fEMinDamage=fEMinDamage, fEMaxDamage=fEMaxDamage)
    cumulatedEnergyDep['doesdamage'] = cumulatedEnergyDep['edep'].apply(
        check_edep)
    filtered_cumul_edep = cumulatedEnergyDep.loc[cumulatedEnergyDep['doesdamage']]
    return filtered_cumul_edep


def calcIndirectDamage(indirect: dict, probIndirect: float, T0: cKDTree, T1: cKDTree, keys_df, out_path):

    eventNo = indirect['step2_eventID'].array(library='np')
    meventNo = np.memmap(out_path+'_temp_eventno_ind.dat',
                         shape=eventNo.shape, mode='w+', dtype=eventNo.dtype)
    meventNo[:] = eventNo[:]
    del eventNo
    # evts_step2 = keys_df['step2'].to_numpy()

    # selected_indirect = np.isin(meventNo, evts_step2)
    # evt2 = meventNo[selected_indirect]

    dnamolecule = indirect['DNAmolecule'].array(library='np')
    mdnamol = np.memmap(out_path+'_temp_mdnamol.dat',
                        shape=dnamolecule.shape, dtype=dnamolecule.dtype, mode='w+')
    mdnamol[:] = dnamolecule[:]
    del dnamolecule

    # selected_dnamol = mdnamol[selected_indirect]

    radicals = indirect['radical'].array(library='np')
    mradicals = np.memmap(out_path+'_temp_radicals.dat',
                          shape=radicals.shape, dtype=radicals.dtype, mode='w+')
    mradicals[:] = radicals[:]
    del radicals
    # selected_radicals = mradicals[mradicals]

    if len(mdnamol) == 0 or len(mradicals) == 0:
        return [], [], []
    selected_radical_idx = np.where((mradicals == 'OH^0') & (
        mdnamol == 'Deoxyribose^0'))[0]
    x = indirect['x'].array(library='np')[selected_radical_idx]
    y = indirect['y'].array(library='np')[selected_radical_idx]
    z = indirect['z'].array(library='np')[selected_radical_idx]

    rand_arr = np.random.rand(selected_radical_idx.shape[0])
    rand_indx = rand_arr <= probIndirect

    points = np.vstack([x, y, z]).swapaxes(0, 1)[rand_indx]
    d0, idx0 = T0.query(points)
    d1, idx1 = T1.query(points)
    get_index = np.vectorize(
        getIndex_v, otypes=[int], signature='(),(),(),()->(2)')
    results = get_index(d0, idx0, d1, idx1)
    copyListIndirect = results[:, 0]
    strandListIndirect = results[:, 1]

    evt2_filtered = meventNo[selected_radical_idx][rand_indx]
    keys_map = keys_df.loc[keys_df['step2_eventID'].isin(
        evt2_filtered)].set_index('step2_eventID')

    events_indirect = keys_map.loc[np.array(evt2_filtered)]
    # events_indirect.loc[(events_indirect['step1_particleID'] == 3) & (
    #     events_indirect['step1_primaryID'] != 1), 'step1_particleID'] = 1
    events_indirect['strand'] = strandListIndirect
    events_indirect['copy'] = copyListIndirect
    
    os.remove(out_path+'_temp_eventno_ind.dat')
    os.remove(out_path+'_temp_mdnamol.dat')
    os.remove(out_path+'_temp_radicals.dat')

    return events_indirect


def runClustering(filename_DNA: str, outputFilename: str, fEMinDamage: float, fEMaxDamage: float, probIndirect: float, sugarFname: str, ndiv_R: int, ndiv_Z: int, spacing: float, primaryParticle: str = None, filenamePhoton: Union[str, bool] = False,  separate_r=False, start_R=10.5):
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
    particleMap = {"alpha": 1,
                   "gamma": 2,
                   "e-": 3,
                   "nu_e": 4,
                   "At211": 5,
                   "Po211": 6,
                   "Bi207": 7,
                   "Pb207": 8,
                   "all": -1
                   }
    print("runClustering UP")
    #get particle ID to use to select data
    if primaryParticle:
        if primaryParticle in particleMap.keys():
            primaryParticleID = particleMap[primaryParticle]
            print("primary: ", primaryParticle, primaryParticleID)
        elif primaryParticle == "all":
            primaryParticleID = -1
            print("primary: ", primaryParticle, primaryParticleID)
        else:
            raise ValueError(f"particle not in particleMap: {primaryParticle}")

    # if not primaryParticle else os.path.splitext(outputFilename)[0]+f"_{primaryParticle}"
    out_path = os.path.splitext(outputFilename)[0]
    #read positions of chromatin segment components
    T0, T1 = readSugarFile(sugarFname)

    # Load data to analyse
    ufile = up.open(filename_DNA)

    info = ufile["output/Info"]
    ps_data = ufile["output/PS_data"].arrays(library='np')
    
    keys_df = pd.DataFrame.from_dict(ps_data)

    pathLength = {}
    chromatinVolume = info['ChromatinVolume_m3'].array(library='np')[
        0]  # in m3
    numBP = info['NumBasepairs'].array(library='np')[0]
    gitHash = info['GitHash'].array(library='np')[0]
    clusteringGitHash = readClusteringGitHash(builddir)

    assert T0.n == numBP, "T0 data does not match numGBP"

    if filenamePhoton:
        LET = "N/A"
    else:
        if hasattr(info, "MeanLET"):
            LET = float(info['MeanLET'])
        else:
            LET = "N/A"

    eventEdep = ufile["output/EventEdep"]
    #calculate dose: output: df with various ID, mean energy and dose
    energy, evt_copyNo_Ke_dose_df = calculateDose(eventEdep, chromatinVolume)

    print("calculateDose Done")

    # this can take long, so just checking when it's done
    print("loading Direct tree")
    direct = ufile["output/Direct"]

    #calculate cumulatedEnergyDep: output is a DF
    cumulatedEnergyDep_df = AccumulateEdep(direct, T0, T1, keys_df, out_path)

    print("AccumulateEdep done")

    # Direct and indirect damage: output DF
    events_strand_copy_direct_df = calcDirectDamage(
        cumulatedEnergyDep_df, fEMinDamage, fEMaxDamage)
    if not len(events_strand_copy_direct_df):
        raise ValueError("no events in direct list.")

    indirect = ufile["output/Indirect"]
    events_strand_copy_indirect_df = calcIndirectDamage(
        indirect, probIndirect, T0, T1, keys_df, out_path)
    if not len(events_strand_copy_indirect_df):
        raise ValueError("no events in indirect list.")
    print("Direct and indirect damage calculated. Running clustering.\n")

    if separate_r:
        from ranges_radii_dense import get_ranges_radii
        ranges_radii =get_ranges_radii(ndiv_R = ndiv_R, ndiv_Z = ndiv_Z, spacing = spacing, start_R = start_R)[0]

        print("separating results per radius")
        for r in ranges_radii:

            check_copynop = partial(
                lambda n, range_r: n in range_r, range_r=ranges_radii[r])
            # select events per copyNo

            events_direct = events_strand_copy_direct_df.loc[events_strand_copy_direct_df['step1_copyNo'].apply(
                check_copynop)]
            events_indirect = events_strand_copy_indirect_df.loc[events_strand_copy_indirect_df['step1_copyNo'].apply(
                check_copynop)]
            events_tot = evt_copyNo_Ke_dose_df.loc[evt_copyNo_Ke_dose_df['step1_copyNo'].apply(
                check_copynop)]

            # get all events (direct+indirect)
            # events_tot_ids = events_tot['step1'].to_numpy()

            unique_pids = list(set(np.concatenate(
                [np.unique(events_direct['step1_primaryID']), np.unique(events_indirect['step1_primaryID'])])))
            unique_pids.append(-1)
            #run clustering on selected events
            for particle in unique_pids:
                if particle == -1:
                    events_direct_part = events_direct[events_direct['step1_primaryID'] != -1]
                    events_indirect_part = events_indirect[events_indirect['step1_primaryID'] != -1]
                else:
                    events_direct_part = events_direct[events_direct['step1_primaryID'] == particle]
                    events_indirect_part = events_indirect[events_indirect['step1_primaryID'] == particle]

                numEvts = list(set(np.concatenate(
                    [events_direct_part['step1_eventID'], events_indirect_part['step1_eventID']])))
                tempResults_part = clustering(numEvts, np.array(events_direct_part['step1_eventID']), np.array(events_direct_part['copy']), np.array(events_direct_part['strand']),
                                              np.array(events_indirect_part['step1_eventID']), np.array(events_indirect_part['copy']), np.array(events_indirect_part['strand']))
                #print(tempResults)
                # strand break number results
                clusteringResults = np.asarray(tempResults_part[0], dtype=int)
                # DSB cluster size, cluster of size 1-10
                clusterSize = np.asarray(tempResults_part[1])
                #prepare text to write to file
                if len(clusteringResults):
                    eventsWithClusteringResults = clusteringResults[:, 0]
                    # _, idx1, idx2 = np.intersect1d(events_tot_ids, eventsWithClusteringResults, return_indices=True)

                    # text_SB = np.zeros(shape = (events_tot_ids.shape[0], 12), dtype=int)
                    # text_SB[idx1] = clusteringResults[:, 1:][idx2]

                    # text_DSB = np.zeros(shape=(events_tot_ids.shape[0], 50), dtype=int)
                    # text_DSB[idx1] = clusterSize[:,1:][idx2]
                    # _, idx1, idx2 = np.intersect1d(events_tot_ids, eventsWithClusteringResults, return_indices=True)

                    text_SB = np.zeros(
                        shape=(eventsWithClusteringResults.shape[0], 12), dtype=int)
                    text_SB[:] = clusteringResults[:, 1:]

                    text_DSB = np.zeros(
                        shape=(eventsWithClusteringResults.shape[0], 50), dtype=int)
                    text_DSB[:] = clusterSize[:, 1:]

                    particle_name_map = {v: k for k, v in particleMap.items()}
                    # print(particle, unique_pids)
                    outfile_SB = out_path + \
                        f"_{particle_name_map[particle]}"+f"_{r}.csv"
                    #write SB and DSB to file
                    with open(outfile_SB, "w") as f:
                        f.write(
                            "Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,GitHash,ClusteringGitHash\n")
                        f.write("{},{},{},{},{},{},{},{},{}\n".format(filename_DNA, np.mean(energy), LET, len(
                            numEvts), chromatinVolume, numBP, sugarFname, gitHash, clusteringGitHash))
                        f.write("EventNo,DoseGy,ParticleID,MeanKEMeV,TotalSBdirect,SSBdirect,cSSBdirect,DSBdirect,TotalSBindirect,SSBindirect,cSSBindirect,DSBindirect,TotalSBtotal,SSBtotal,cSSBtotal,DSBtotal\n")
                        for i, event in enumerate(sorted(eventsWithClusteringResults)):
                            if particle == -1:

                                # print(event, evt_copyNo_Ke_dose_df['dose'][(evt_copyNo_Ke_dose_df['step1']==event)&(evt_copyNo_Ke_dose_df['pid']!=particle)].sum())
                                f.write("{},{},{},{},{}".format(
                                    event,
                                    events_tot['dose'][(events_tot['step1_eventID'] == event) & (
                                        events_tot['step1_primaryID'] != particle)].sum(),#/events_tot['dose'][(events_tot['step1_eventID'] == event) & (events_tot['step1_primaryID'] != particle)].count(),
                                    # evt_copyNo_Ke_dose_df['pid'][evt_copyNo_Ke_dose_df['step1']==event],
                                    particle,
                                    events_tot['edep_MeV'][(events_tot['step1_eventID'] == event) & (
                                        events_tot['step1_primaryID'] != particle)].sum(),
                                    str(list(text_SB[i])).strip("[").strip("]")))
                                f.write('\n')
                            else:
                                
                                # print(event, particle, events_tot['dose'][(events_tot['step1']==event)&(events_tot['pid']==particle)])
                                f.write("{},{},{},{},{}".format(
                                    event,
                                    events_tot['dose'][(events_tot['step1_eventID'] == event) & (
                                        events_tot['step1_primaryID'] == particle)].sum(),#/events_tot['dose'][(events_tot['step1_eventID'] == event) & (events_tot['step1_primaryID'] == particle)].count(),
                                    # evt_copyNo_Ke_dose_df['pid'][evt_copyNo_Ke_dose_df['step1']==event],
                                    particle,
                                    events_tot['edep_MeV'][(events_tot['step1_eventID'] == event) & (
                                        events_tot['step1_primaryID'] == particle)].sum(),
                                    str(list(text_SB[i])).strip("[").strip("]")))
                                f.write('\n')

                    outfile_DSB = out_path + \
                        f"_{particle_name_map[particle]}"+f"_DSB_{r}.csv"
                    with open(outfile_DSB, "w") as f:
                        f.write(
                            "Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,SimGitHash,ClusteringGitHash\n")
                        f.write("{},{},{},{},{},{},{},{},{}\n".format(filename_DNA, np.mean(energy), LET, len(
                            numEvts), chromatinVolume, numBP, sugarFname, gitHash, clusteringGitHash))
                        f.write("EventNo,DoseGy,ParticleID,Direct_1_SBperDSBcluster,Direct_2_SBperDSBcluster,Direct_3_SBperDSBcluster,Direct_4_SBperDSBcluster,Direct_5_SBperDSBcluster,Direct_6_SBperDSBcluster,Direct_7_SBperDSBcluster,Direct_8_SBperDSBcluster,Direct_9_SBperDSBcluster,Direct_10_SBperDSBcluster,Indirect_1_SBperDSBcluster,Indirect_2_SBperDSBcluster,Indirect_3_SBperDSBcluster,Indirect_4_SBperDSBcluster,Indirect_5_SBperDSBcluster,Indirect_6_SBperDSBcluster,Indirect_7_SBperDSBcluster,Indirect_8_SBperDSBcluster,Indirect_9_SBperDSBcluster,Indirect_10_SBperDSBcluster,Hybrid_1_SBperDSBcluster,Hybrid_2_SBperDSBcluster,Hybrid_3_SBperDSBcluster,Hybrid_4_SBperDSBcluster,Hybrid_5_SBperDSBcluster,Hybrid_6_SBperDSBcluster,Hybrid_7_SBperDSBcluster,Hybrid_8_SBperDSBcluster,Hybrid_9_SBperDSBcluster,Hybrid_10_SBperDSBcluster,Mixed_1_SBperDSBcluster,Mixed_2_SBperDSBcluster,Mixed_3_SBperDSBcluster,Mixed_4_SBperDSBcluster,Mixed_5_SBperDSBcluster,Mixed_6_SBperDSBcluster,Mixed_7_SBperDSBcluster,Mixed_8_SBperDSBcluster,Mixed_9_SBperDSBcluster,Mixed_10_SBperDSBcluster,Total_1_SBperDSBcluster,Total_2_SBperDSBcluster,Total_3_SBperDSBcluster,Total_4_SBperDSBcluster,Total_5_SBperDSBcluster,Total_6_SBperDSBcluster,Total_7_SBperDSBcluster,Total_8_SBperDSBcluster,Total_9_SBperDSBcluster,Total_10_SBperDSBcluster\n")
                        for i, event in enumerate(sorted(eventsWithClusteringResults)):
                            if particle == -1:
                                f.write("{},{},{},{}".format(
                                    event, events_tot['dose'][(events_tot['step1_eventID'] == event) & (
                                        events_tot['step1_primaryID'] != particle)].sum(),#/events_tot['dose'][(events_tot['step1_eventID'] == event) & (events_tot['step1_primaryID'] != particle)].count(),
                                    particle,
                                    str(list(text_DSB[i])).strip("[").strip("]")))
                                f.write('\n')
                            else:
                                f.write("{},{},{},{}".format(
                                    event, events_tot['dose'][(events_tot['step1_eventID'] == event) & (
                                        events_tot['step1_primaryID'] == particle)].sum(),#/events_tot['dose'][(events_tot['step1_eventID'] == event) & (events_tot['step1_primaryID'] == particle)].count(),
                                    particle,
                                    str(list(text_DSB[i])).strip("[").strip("]")))
                                f.write('\n')

                else:
                    print(f"no clustering results for radius {r}")
                    continue

    else:
        pass
        # print("Calculate clustering without sepration on radius.")
        # events_direct = events_strand_copy_direct_df
        # events_indirect = events_strand_copy_indirect_df
        # events_tot = evt_copyNo_Ke_dose_df

        # numEvts = list(set(np.concatenate([events_direct['step1'], events_indirect['step1']])))
        # events_tot_ids = events_tot['step1'].to_numpy()
        # tempResults = clustering(numEvts, np.array(events_direct['step1']), np.array(events_direct['copy']), np.array(events_direct['strand']),
        #                              np.array(events_indirect['step1']), np.array(events_indirect['copy']), np.array(events_indirect['strand']))
        # clusteringResults = np.asarray(tempResults[0], dtype=int)  # strand break number results
        # clusterSize = np.asarray(tempResults[1])  # DSB cluster size, cluster of size 1-10

        # if len(clusteringResults):
        #     eventsWithClusteringResults = clusteringResults[:,0]
        #     _, idx1, idx2 = np.intersect1d(events_tot_ids, eventsWithClusteringResults, return_indices=True)

        #     text_SB = np.zeros(shape = (events_tot_ids.shape[0], 12), dtype=int)
        #     text_SB[idx1] = clusteringResults[:, 1:][idx2]

        #     text_DSB = np.zeros(shape=(events_tot_ids.shape[0], 50), dtype=int)
        #     text_DSB[idx1] = clusterSize[:,1:][idx2]

        #     outfile_tot = outputFilename
        #     print(f"output_file: {outfile_tot}")

        #     with open(outfile_tot, "w") as f:
        #         f.write(
        #             "Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,GitHash,ClusteringGitHash\n")
        #         f.write("{},{},{},{},{},{},{},{},{}\n".format(filename_DNA, energy, LET, len(
        #             events_tot), chromatinVolume, numBP, sugarFname, gitHash, clusteringGitHash))
        #         f.write("EventNo,DoseGy,ParticleID,MeanKEMeV,TotalSBdirect,SSBdirect,cSSBdirect,DSBdirect,TotalSBindirect,SSBindirect,cSSBindirect,DSBindirect,TotalSBtotal,SSBtotal,cSSBtotal,DSBtotal\n")
        #         for i, event in enumerate(events_tot_ids):
        #             f.write("{},{},{},{},{}".format(
        #                 event,
        #                 evt_copyNo_Ke_dose_df['dose'][i], evt_copyNo_Ke_dose_df['pid'][i],
        #                 evt_copyNo_Ke_dose_df['edep_J'][i],
        #                 str(list(text_SB[i])).strip("[").strip("]")))
        #             f.write('\n')

        #     outfile_DSB_tot = out_path+"_DSB.csv"
        #     with open(outfile_DSB_tot, "w") as f:
        #         f.write(
        #             "Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,SimGitHash,ClusteringGitHash\n")
        #         f.write("{},{},{},{},{},{},{},{},{}\n".format(filename_DNA, energy, LET, len(
        #                 events_tot), chromatinVolume, numBP, sugarFname, gitHash, clusteringGitHash))
        #         f.write("EventNo,DoseGy,ParticleID,Direct_1_SBperDSBcluster,Direct_2_SBperDSBcluster,Direct_3_SBperDSBcluster,Direct_4_SBperDSBcluster,Direct_5_SBperDSBcluster,Direct_6_SBperDSBcluster,Direct_7_SBperDSBcluster,Direct_8_SBperDSBcluster,Direct_9_SBperDSBcluster,Direct_10_SBperDSBcluster,Indirect_1_SBperDSBcluster,Indirect_2_SBperDSBcluster,Indirect_3_SBperDSBcluster,Indirect_4_SBperDSBcluster,Indirect_5_SBperDSBcluster,Indirect_6_SBperDSBcluster,Indirect_7_SBperDSBcluster,Indirect_8_SBperDSBcluster,Indirect_9_SBperDSBcluster,Indirect_10_SBperDSBcluster,Hybrid_1_SBperDSBcluster,Hybrid_2_SBperDSBcluster,Hybrid_3_SBperDSBcluster,Hybrid_4_SBperDSBcluster,Hybrid_5_SBperDSBcluster,Hybrid_6_SBperDSBcluster,Hybrid_7_SBperDSBcluster,Hybrid_8_SBperDSBcluster,Hybrid_9_SBperDSBcluster,Hybrid_10_SBperDSBcluster,Mixed_1_SBperDSBcluster,Mixed_2_SBperDSBcluster,Mixed_3_SBperDSBcluster,Mixed_4_SBperDSBcluster,Mixed_5_SBperDSBcluster,Mixed_6_SBperDSBcluster,Mixed_7_SBperDSBcluster,Mixed_8_SBperDSBcluster,Mixed_9_SBperDSBcluster,Mixed_10_SBperDSBcluster,Total_1_SBperDSBcluster,Total_2_SBperDSBcluster,Total_3_SBperDSBcluster,Total_4_SBperDSBcluster,Total_5_SBperDSBcluster,Total_6_SBperDSBcluster,Total_7_SBperDSBcluster,Total_8_SBperDSBcluster,Total_9_SBperDSBcluster,Total_10_SBperDSBcluster\n")
        #         for i, event in enumerate(events_tot_ids):
        #             f.write("{},{},{},{}".format(
        #                 event, evt_copyNo_Ke_dose_df['dose'][i], evt_copyNo_Ke_dose_df['pid'][i], str(list(text_DSB[i])).strip("[").strip("]")))
        #             f.write('\n')
        # else:
        #     raise ValueError(f"no clustering results.")

    print("Finished: {}".format(filename_DNA))


if __name__ == "__main__":
    pass
