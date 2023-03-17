import matplotlib.pyplot as plt
from readResults import *
import argparse
import os
import numpy as np

def plot_results(folder, fname_prefix, spacing, nevents, savefig=True, n_div_r=20, seed=None, out_folder = "./", damage_type='SSB', particle=None):
    radii = range(n_div_r)
    
    fname_base = f"{fname_prefix}_{spacing}um_"+f"{seed}"
    if particle: 
        fname_base+=f"_{particle}"

    if damage_type == "SSB":
        fname = fname_base
    elif damage_type == "DSB":
        fname = fname_base+f"_DSB"
    fnames = [os.path.join(folder, fname+f"_{r}.csv") for r in radii]
    
    dataset = {}

    for i, fname in zip(radii, fnames):
        to_combine = []

        to_combine.append(readIN(fname))

        if len(to_combine) != 0:
            dataset[i] = combine(to_combine)

    # plot strand breaks vs LET & energy per dose
    
    total = np.zeros(n_div_r)
    direct = np.zeros(n_div_r)
    indirect = np.zeros(n_div_r)

    dose_out = np.zeros(n_div_r)
    radii_out = np.zeros(n_div_r)


    for i in radii:
        if damage_type=='SSB':
        #en.append(getEnergy(dataset[i]))
            nSSB, dose, nGBP = calcBreaksperDose("TotalSBtotal", (dataset[i]))
            total[i]=nSSB/nGBP
            nSSB, dose, nGBP = calcBreaksperDose("TotalSBdirect", (dataset[i]))
            direct[i] = nSSB/nGBP
            nSSB, dose, nGBP = calcBreaksperDose("TotalSBindirect", (dataset[i]))
            indirect[i]= nSSB/nGBP
        
        elif damage_type=='DSB':
            nDSB, dose, nGBP = calcSimpleDSBperDose("Total", (dataset[i]))
            total[i] = nDSB/nGBP
            nDSB, dose, nGBP = calcSimpleDSBperDose("Direct", (dataset[i]))
            direct[i] = nDSB/nGBP
            nDSB, dose, nGBP = calcSimpleDSBperDose("Indirect", (dataset[i]))
            indirect[i] = nDSB/nGBP
        
        radii_out[i] = 10.5+i*spacing  
        dose_out[i] = dose if dose else 0


    # not dividing by dose
    plt.figure(figsize=(8, 5))
    plt.plot(radii_out, total, color='k',
                label=f'Total SB', marker='x', linewidth=.6, markersize=3)
    plt.plot(radii_out, indirect,
                color='mediumvioletred', label=f'Indirect SB',linewidth=.6, marker='o', markersize=3)
    plt.plot(radii_out, direct, color='cornflowerblue',
            label=f'Direct SB', linewidth=.6,  marker='s', markersize=3)
    plt.legend()
    spacing_str = spacing if spacing else 10
    plt.title(f"DNA damage vs distance, spacing: {spacing_str} um")
    plt.xlabel("Radial distance (um)")
    #plt.ylabel('Number of strand breaks $Gbp^{-1} per decay$)', fontsize=11)
    plt.ylabel(f'Number of {damage_type} ($Gbp^{-1}$)', fontsize=11)
    if savefig:
        plt.savefig(os.path.join(out_folder, f"{damage_type}"+(f"_{particle}" if particle else None)+f"_abs_{fname_prefix}_{spacing}um_{seed}.png"))
    else:
        plt.show()

    # Dividing by dose
    plt.figure(figsize=(8, 5))
    plt.plot(radii_out, np.divide(total, dose_out, out=np.zeros_like(total),where=dose_out!=0), color='k',
                label=f'Total SB', marker='x', linewidth=.6, markersize=3)
    plt.plot(radii_out, np.divide(indirect, dose_out, out=np.zeros_like(indirect),where=dose_out!=0),
                color='mediumvioletred', label=f'Indirect SB',linewidth=.6, marker='o', markersize=3)
    plt.plot(radii_out, np.divide(direct, dose_out, out=np.zeros_like(direct),where=dose_out!=0), color='cornflowerblue',
            label=f'Direct SB', linewidth=.6,  marker='s', markersize=3)
    plt.legend()
    spacing_str = spacing if spacing else 10
    plt.title(f"DNA damage vs distance, spacing: {spacing_str} um")
    plt.xlabel("Radial distance (um)")
    #plt.ylabel('Number of strand breaks $Gbp^{-1} per decay$)', fontsize=11)
    plt.ylabel('Number of {damage_type} ($Gy^{-1} Gbp^{-1}$)', fontsize=11)
    if savefig:
        plt.savefig(os.path.join(out_folder, f"{damage_type}"+(f"_{particle}" if particle else None)+f"_byDose_{fname_prefix}_{spacing}um_{seed}.png"))
    else:
        plt.show()



    #DOSE

    fig2 = plt.figure(figsize=(8, 5))
    plt.plot(radii_out, dose_out, color='darkblue', marker='s', markersize=3, linestyle='-', linewidth=.9,
                label=f'Dose')
    plt.legend()
    plt.title(f"Dose vs distance, spacing: {spacing_str} um")
    plt.xlabel("Radial distance (um)")
    plt.ylabel('Dose ($Gy$)', fontsize=11)
    if savefig:
        plt.savefig(os.path.join(out_folder, f"dose_{damage_type}"+(f"_{particle}" if particle else None)+f"_{fname_prefix}_{spacing}um_{seed}.png"))
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
    parser.add_argument('--type', type=str, help="SSB or DSB")
    parser.add_argument("--particle", type=str,default=None)
    
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
    if args.type:
        damage_type= args.type
    if args.particle:
        particle = args.particle

    
    if len(spacing)>1:
        for s in spacing:
            plot_results(folder=folder, fname_prefix=fname_prefix, 
                         nevents=nevents, spacing=int(s), seed=seed, out_folder=out_folder, n_div_r=n_div_r, damage_type=damage_type, particle=particle)
    else:
        plot_results(folder=folder, fname_prefix=fname_prefix, 
                         nevents=nevents, spacing=int(spacing[0]), seed=seed, out_folder=out_folder, n_div_r=n_div_r,damage_type=damage_type, particle=particle)



