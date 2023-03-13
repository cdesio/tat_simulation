import matplotlib.pyplot as plt
from readResults import *
import argparse
import os
import numpy as np

def plot_results(folder, fname_prefix, spacing, nevents, savefig=True, n_div_r=20, seed=None, out_folder = "./"):
    radii = range(n_div_r)
    if spacing:
        if seed:
            fnames_SSB = [os.path.join(
            folder, f"{fname_prefix}_{spacing}um_"+f"{seed}"+f"_{r}.csv") for r in radii]
            fnames_DSB = [os.path.join(
            folder, f"{fname_prefix}_{spacing}um_"+f"{seed}"+f"_DSB_{r}.csv") for r in radii]
        else:    
            fnames_SSB = [os.path.join(
                folder, f"{fname_prefix}_{spacing}um_{r}.csv") for r in radii]
            fnames_DSB = [os.path.join(
                folder, f"{fname_prefix}_{spacing}um_DSB_{r}.csv") for r in radii]
    else:
        fnames_SSB = [os.path.join(
            folder, f"{fname_prefix}_{r}.csv") for r in radii]
        fnames_DSB = [os.path.join(
            folder, f"{fname_prefix}_DSB_{r}.csv") for r in radii]
    
    dataset_SSB = {}
    dataset_DSB = {}


    for i, fSSB, fDSB in zip(radii, fnames_SSB, fnames_DSB):
        to_combine_SSB = []
        to_combine_DSB = []
        to_combine_SSB.append(readIN(fSSB))
        to_combine_DSB.append(readIN(fDSB))
        if len(to_combine_SSB) != 0:
            dataset_SSB[i] = combine(to_combine_SSB)
        if len(to_combine_DSB) != 0:
            dataset_DSB[i] = combine(to_combine_DSB)
    # plot strand breaks vs LET & energy per dose
    
    total_SSB = []
    direct_SSB = []
    indirect_SSB = []

    dose_SSB_tot = []
    dose_SSB_direct = []
    dose_SSB_indirect = []

    total_DSB = []
    direct_DSB = []
    indirect_DSB = []

    dose_DSB_tot = []
    dose_DSB_direct = []
    dose_DSB_indirect = []

    totalErr = []
    directErr = []
    indirectErr = []

    radii_out = []

    markers = ['s', '^', "o", "*", "x"]


    for i in radii:

        #en.append(getEnergy(dataset[i]))
        nSSB, dose, nGBP = calcBreaksperDose("TotalSBtotal", (dataset_SSB[i]))
        total_SSB.append(nSSB/nGBP)
        dose_SSB_tot.append(dose if dose else 0)
        nSSB, dose, nGBP = calcBreaksperDose("TotalSBdirect", (dataset_SSB[i]))
        direct_SSB.append(nSSB/nGBP)
        dose_SSB_direct.append(dose if dose else 0)
        nSSB, dose, nGBP = calcBreaksperDose("TotalSBindirect", (dataset_SSB[i]))
        indirect_SSB.append(nSSB/nGBP)
        dose_SSB_indirect.append(dose if dose else 0)
        
        nDSB = calcSimpleDSBperDose("Total", (dataset_DSB[i]))
        total_DSB.append(nDSB)
        #dose_DSB_tot.append(dose if dose else 0)
        nDSB = calcSimpleDSBperDose("Direct", (dataset_DSB[i]))
        direct_DSB.append(nDSB)
        #dose_DSB_direct.append(dose if dose else 0)
        nDSB = calcSimpleDSBperDose("Indirect", (dataset_DSB[i]))
        indirect_DSB.append(nDSB)
        #dose_DSB_indirect.append(dose if dose else 0)


        radii_out.append(10.5+i*spacing)  


    total_SSB = np.asarray(total_SSB)
    direct_SSB = np.asarray(direct_SSB)
    indirect_SSB = np.asarray(indirect_SSB)

    dose_SSB_tot = np.asarray(dose_SSB_tot)
    dose_SSB_direct = np.asarray(dose_SSB_direct)
    dose_SSB_indirect = np.asarray(dose_SSB_indirect)

    total_DSB = np.asarray(total_DSB)
    direct_DSB = np.asarray(direct_DSB)
    indirect_DSB = np.asarray(indirect_DSB)

    dose_DSB_tot = np.asarray(dose_DSB_tot)
    dose_DSB_direct = np.asarray(dose_DSB_direct)
    dose_DSB_indirect = np.asarray(dose_DSB_indirect)

    # not dividing by dose
    plt.figure(figsize=(8, 5))
    plt.plot(radii_out, total_SSB, color='k',
                label=f'Total SB', marker='x', linewidth=.6, markersize=3)
    plt.plot(radii_out, indirect_SSB,
                color='mediumvioletred', label=f'Indirect SB',linewidth=.6, marker='o', markersize=3)
    plt.plot(radii_out, direct_SSB, color='cornflowerblue',
            label=f'Direct SB', linewidth=.6,  marker='s', markersize=3)
    plt.legend()
    spacing_str = spacing if spacing else 10
    plt.title(f"DNA damage vs distance, spacing: {spacing_str} um")
    plt.xlabel("Radial distance (um)")
    #plt.ylabel('Number of strand breaks $Gbp^{-1} per decay$)', fontsize=11)
    plt.ylabel('Number of strand breaks ($Gbp^{-1}$)', fontsize=11)
    if savefig:
        plt.savefig(os.path.join(out_folder, f"SSB_{fname_prefix}_{spacing}um_{seed}.png"))
    else:
        plt.show()

    # Dividing by dose
    plt.figure(figsize=(8, 5))
    plt.plot(radii_out, total_SSB/dose_SSB_tot, color='k',
                label=f'Total SB', marker='x', linewidth=.6, markersize=3)
    plt.plot(radii_out, indirect_SSB/dose_SSB_indirect,
                color='mediumvioletred', label=f'Indirect SB',linewidth=.6, marker='o', markersize=3)
    plt.plot(radii_out, direct_SSB/dose_SSB_direct, color='cornflowerblue',
            label=f'Direct SB', linewidth=.6,  marker='s', markersize=3)
    plt.legend()
    spacing_str = spacing if spacing else 10
    plt.title(f"DNA damage vs distance, spacing: {spacing_str} um")
    plt.xlabel("Radial distance (um)")
    #plt.ylabel('Number of strand breaks $Gbp^{-1} per decay$)', fontsize=11)
    plt.ylabel('Number of strand breaks ($Gy^{-1} Gbp^{-1}$)', fontsize=11)
    if savefig:
        plt.savefig(os.path.join(out_folder, f"SSB_byDose_{fname_prefix}_{spacing}um_{seed}.png"))
    else:
        plt.show()


 # Dividing by dose
    plt.figure(figsize=(8, 5))
    plt.plot(radii_out, total_DSB, color='k',
                label=f'Total DSB', marker='x', linewidth=.6, markersize=3)
    plt.plot(radii_out, indirect_DSB,
                color='mediumvioletred', label=f'Indirect DSB',linewidth=.6, marker='o', markersize=3)
    plt.plot(radii_out, direct_DSB, color='cornflowerblue',
            label=f'Direct DSB', linewidth=.6,  marker='s', markersize=3)
    plt.legend()
    spacing_str = spacing if spacing else 10
    plt.title(f"DNA damage vs distance, spacing: {spacing_str} um")
    plt.xlabel("Radial distance (um)")
    #plt.ylabel('Number of strand breaks $Gbp^{-1} per decay$)', fontsize=11)
    plt.ylabel('Number of Double Strand Breaks ($Gy^{-1} Gbp^{-1}$)', fontsize=11)
    if savefig:
        plt.savefig(os.path.join(out_folder, f"DSB_byDose_{fname_prefix}_{spacing}um_{seed}.png"))
    else:
        plt.show()


    #DOSE

    fig2 = plt.figure(figsize=(8, 5))
    plt.plot(radii_out, dose_SSB_tot, color='darkblue', marker='s', markersize=3, linestyle='-', linewidth=.9,
                label=f'Dose')
    plt.legend()
    # plt.xlabel('energy (MeV/$\mu$m)',fontsize=11)
    # plt.xlabel('LET (keV/$\mu$m)',fontsize=11)
    plt.title(f"Dose vs distance, spacing: {spacing_str} um")
    plt.xlabel("Radial distance (um)")
    # plt.ylabel('Number of strand breaks $Gbp^{-1} per decay$)', fontsize=11)
    plt.ylabel('Dose ($Gy$)', fontsize=11)
    if savefig:
        plt.savefig(os.path.join(out_folder, f"dose_{fname_prefix}_{spacing}um.png"))
    else:
        plt.show()

    return



if __name__ == "__main__":


    parser = argparse.ArgumentParser()

    parser.add_argument("--folder",  type=str ,default = '/', 
                        help="input folder for plots: clustering folder")
    parser.add_argument("--fname_prefix",  type=str, 
                        help="input filename prefix")
    parser.add_argument("--out_folder",  type=str, default = "./",
                        help="output folder path: def: ./ ")
    parser.add_argument("--n",  type=int, help="number of events")
    parser.add_argument("--n_div_R",  type=int, default = 20, help="number of events")
    parser.add_argument("--spacing", nargs="+",  type=float,default = [2.0, 5.0, 10.0], 
                    help="spacing in micrometer for the dna voxels")
    parser.add_argument("--savefig",  type=bool, default = True)  
    parser.add_argument("--seed", type=int, default = None, help='seed of simulation')
    
    args = parser.parse_args()
    
    if args.folder:
        folder = args.folder
    if args.fname_prefix:
        fname_prefix = args.fname_prefix
    if args.spacing:
        spacing = args.spacing
    if args.n:
        nevents = args.n
    if args.savefig:
        savefig = args.savefig
    if args.n_div_R:
        n_div_r = args.n_div_R
    if args.seed:
        seed = args.seed
    if args.out_folder:
        out_folder = args.out_folder

    
    if len(spacing)>1:
        for s in spacing:
            plot_results(folder=folder, fname_prefix=fname_prefix, 
                         nevents=nevents, spacing=int(s), seed=seed, out_folder=out_folder, n_div_r=n_div_r)
    else:
        plot_results(folder=folder, fname_prefix=fname_prefix, 
                         nevents=nevents, spacing=int(spacing[0]), seed=seed, out_folder=out_folder, n_div_r=n_div_r)



