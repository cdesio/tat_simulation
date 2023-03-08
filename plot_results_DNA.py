import matplotlib.pyplot as plt
from readResults import *
import argparse
import os
import numpy as np

def plot_results(folder, fname_prefix, spacing, nevents, savefig=True, n_div_r=20, seed=None, out_folder = "./"):
    radii = range(n_div_r)
    if spacing:
        if seed:
            fnames = [os.path.join(
            folder, f"{fname_prefix}_{spacing}um_"+f"{seed}"+f"_{r}.csv") for r in radii]
        else:    
            fnames = [os.path.join(
                folder, f"{fname_prefix}_{spacing}um_{r}.csv") for r in radii]
    else:
        fnames = [os.path.join(
            folder, f"{fname_prefix}_{r}.csv") for r in radii]
    dataset = {}

    for i, f in zip(radii, fnames):
        to_combine = []
        to_combine.append(readIN(f))
        if len(to_combine) != 0:
            dataset[i] = combine(to_combine)
    # plot strand breaks vs LET & energy per dose
    en = []
    total = []
    direct = []
    indirect = []

    totalErr = []
    directErr = []
    indirectErr = []

    dose_tot = []
    dose_direct = []
    dose_indirect = []

    radii_out = []

    markers = ['s', '^', "o", "*", "x"]


    for i in radii:

        #en.append(getEnergy(dataset[i]))
        a, d, n = calcBreaksperDose("TotalSBtotal", (dataset[i]))
        total.append(a/d/n if d else 0)
        dose_tot.append(d if d else 0)
        a,d, n = calcBreaksperDose("TotalSBdirect", (dataset[i]))
        direct.append(a/d/n if d else 0)
        dose_direct.append(d if d else 0)
        a,d, n = calcBreaksperDose("TotalSBindirect", (dataset[i]))
        indirect.append(a/d/n if d else 0)
        dose_indirect.append(d if d else 0)
        
        radii_out.append(10.5+i*spacing)   

    dose_tot = np.asarray(dose_tot)
    dose_direct = np.asarray(dose_direct)
    dose_indirect = np.asarray(dose_indirect)

    fig1 = plt.figure(figsize=(8, 5))
    plt.plot(radii_out, np.array(total), color='k',
                label=f'Total SB', marker='x', linewidth=.6, markersize=3)
    plt.plot(radii_out, np.array(indirect),
                color='mediumvioletred', label=f'Indirect SB',linewidth=.6, marker='o', markersize=3)
    plt.plot(radii_out, np.array(direct), color='cornflowerblue',
            label=f'Direct SB', linewidth=.6,  marker='s', markersize=3)
    plt.legend()
    spacing_str = spacing if spacing else 10
    plt.title(f"DNA damage vs distance, spacing: {spacing_str} um")
    plt.xlabel("Radial distance (um)")
    #plt.ylabel('Number of strand breaks $Gbp^{-1} per decay$)', fontsize=11)
    plt.ylabel('Number of strand breaks ($Gy^{-1} Gbp^{-1}$)', fontsize=11)
    if savefig:
        plt.savefig(os.path.join(out_folder, f"SSB_{fname_prefix}_{spacing}um.png"))
    plt.show()

    #DOSE

    fig2 = plt.figure(figsize=(8, 5))
    plt.plot(radii_out, dose_tot, color='darkblue', marker='s', markersize=3, linestyle='-', linewidth=.9,
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



