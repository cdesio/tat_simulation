"""
Create all files required to run alpha simulation(s) with set seeds and energies. simulation & clustering run separately so that multithreading can be used on simulation
"""

import numpy as np
import os
import random
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--spacing", nargs="+",  type=float ,default = [2.0, 5.0, 10.0], 
                    help="spacing in micrometer for the dna voxels")

parser.add_argument("--n_div_R", type=int,
                    help="n of divisions in R", default = 20)
parser.add_argument("--n_div_Z", type=int,
                    help="n of divisions in R", default = 10)
parser.add_argument("--n_div_theta", type=int,
                    help="n of divisions in R", default = 80)

parser.add_argument("--startR", type=float, default = 10.5)

parser.add_argument("--gpsHalfZ", type = float, default = 3.0)
parser.add_argument("--n", type=int, default = 100)
parser.add_argument("--slurm", type=bool, required=False, default = False)

args = parser.parse_args()



# parameters to set
numThread = 4

if args.slurm:
    slurm = args.slurm
else:
    slurm = False

seed = random.randint(1, 10000)

if slurm:
    simulation_parent = os.path.join(
        "/", "user", "work", "yw18581", "DaRT", "TAT", "tat_simu_2steps"
    )
else:
    simulation_parent = os.path.join(
        "/", "home", "cdesio", "TAT", "tat_simu_2steps"
    )

histoneFile = os.path.join(
    simulation_parent, "geometryFiles", "histonePositions_4x4_300nm.csv"
)
sugarFile = os.path.join(
    simulation_parent, "geometryFiles", "sugarPos_4x4_300nm.csv"
)

time = "0-10:00"
mem = 10
targetSize = 300  # nanometers
gpsRadius = 10 # micrometers

if args.gpsHalfZ:
    gpsHalfZ = args.gpsHalfZ  # micrometers

if args.n_div_R:
    n_div_R = args.n_div_R

if args.n_div_theta:
    n_div_theta = args.n_div_theta
if args.n_div_Z:
    n_div_Z = args.n_div_Z

nboxes = n_div_R*n_div_theta*n_div_Z
nboxes_string = f"{int(nboxes/1000)}k" if nboxes/1000 >=1 else f"{nboxes}"

if args.spacing:
    spacing = args.spacing

if args.startR:
    startR = args.startR

vessellength = 10

if args.n:
    numIons = args.n

printProgress = int(numIons/10)

n_string = f"{int(numIons/1000)}k" if numIons/1000 >=1 else f"{numIons}"

folder = f"test_At{n_string}_{nboxes_string}_boxes_{startR}um"

print(f"Creating scripts in folder: {folder} with seed {seed}.")
####

decay_sim_folder = os.path.join(simulation_parent, "decay_simulation")
dna_sim_folder = os.path.join(simulation_parent, "simulation")
clustering_folder = os.path.join(simulation_parent, "Clustering")

makerundir = lambda d: os.path.join(d, "build")

#### create folder for current test
test_dir = os.path.join(simulation_parent, f"output/{folder}")
if not os.path.exists(test_dir):
    os.makedirs(test_dir)

#### create one file per spacing
for s in spacing:
    s_string = f"{int(s)}"

    #  write mac file

    filename_mac = os.path.join(
        test_dir, f"input_Atdecay_{n_string}_spacing_{s_string}um_{seed}.in"
    )  # to be run from folder created
    with open(filename_mac, "w") as f:
        f.write("/run/verbose 2\n")
        f.write("/control/verbose 2\n")
        f.write("\n")
        f.write(f"/det/set_ndiv_R {n_div_R}\n")
        f.write(f"/det/set_ndiv_theta {n_div_theta}\n")
        f.write(f"/det/set_ndiv_Z {n_div_Z}\n")
        f.write("\n")
        f.write(f"/det/set_spacing {s} um\n")
        f.write(f"/det/set_startR {startR} um\n")
        f.write(f"/det/set_vessellength {vessellength} um\n")
        f.write("\n")
        f.write("/run/initialize\n")
        f.write("\n")
    
        f.write("/gps/pos/type Surface\n")
        f.write("/gps/pos/shape Cylinder\n")
        f.write("/gps/pos/centre 0 0 0 um\n")
        f.write(f"/gps/pos/radius {gpsRadius} um\n")
        f.write(f"/gps/pos/halfz {gpsHalfZ} um\n")
        f.write("\n")

        f.write("/gps/ang/type cos\n")
        f.write("/gps/particle ion\n")
        f.write("/gps/ion 85 211\n")
        f.write("/gps/ene/type Mono\n")
        f.write(f"/gps/energy 0 MeV\n")
        f.write("\n")
        f.write(f"/run/printProgress {printProgress}\n")
        f.write(f"/run/beamOn {numIons}\n")

    filename_decay = os.path.join(
    test_dir, f"run_Atdecay_{n_string}_spacing_{s_string}um_{seed}.sh")  # to be run from folder created


    if slurm:
        print(f"sbatch run_Atdecay_{n_string}_spacing_{s_string}um_{seed}.sh")
    with open(filename_decay, "w") as f:
        f.write("#!/bin/bash --login\n")
        if slurm:
            f.write(f"#SBATCH --job-name=Atdecay_spacing_{s_string}um_{seed}\n")
            f.write(f"#SBATCH --output=Atdecay_spacing_{s_string}um_{seed}.out.%J\n")
            f.write(f"#SBATCH --error=Atdecay_spacing_{s_string}um_{seed}.err.%J\n")
            # maximum job time in D-HH:MM
            f.write(f"#SBATCH --time={time}\n")
            f.write("#SBATCH --nodes=1\n")
            f.write("#SBATCH -p short\n")
            f.write(f"#SBATCH --ntasks-per-node={numThread}\n")
            f.write(f"#SBATCH --mem {mem}GB\n")
            f.write("\n")
            f.write("module load apps/geant/4.11.1.0\n")
            f.write("module load apps/root/6.26.00\n")
            f.write("source /user/home/yw18581/.bash_profile\n")
            f.write("source activate dart\n")

        else:
            f.write("conda activate rootpy\n")
            f.write("source /opt/geant4-v11.1.0-install/bin/geant4.sh\n")
        # run decay simulation
        f.write(
            f"time {makerundir(decay_sim_folder)}/decaySim -mac {filename_mac} -out {test_dir}/out_Atdecay_{n_string}_spacing_{s_string}um -seed {seed}\n"
        )
        

        filename_DNA = os.path.join(
        test_dir, f"run_AtDNA_{n_string}_spacing_{s_string}um_{seed}.sh"
    )  # to be run from folder created
    if slurm:
        print(f"sbatch run_AtDNA_{n_string}_spacing_{s_string}um_{seed}.sh")
    with open(filename_DNA, "w") as f:
        f.write("#!/bin/bash --login\n")
        if slurm:
            f.write(f"#SBATCH --job-name=AtDNA_spacing_{s_string}um_{seed}\n")
            f.write(f"#SBATCH --output=AtDNA_spacing_{s_string}um_{seed}.out.%J\n")
            f.write(f"#SBATCH --error=AtDNA_spacing_{s_string}um_{seed}.err.%J\n")
            # maximum job time in D-HH:MM
            f.write(f"#SBATCH --time={time}\n")
            f.write("#SBATCH --nodes=1\n")
            f.write("#SBATCH -p short\n")
            f.write(f"#SBATCH --ntasks-per-node={numThread}\n")
            f.write(f"#SBATCH --mem {mem}GB\n")
            f.write("\n")
            f.write("module load apps/geant/4.11.1.0\n")
            f.write("module load apps/root/6.26.00\n")
            f.write("source /user/home/yw18581/.bash_profile\n")
            f.write("source activate dart\n")

        else:
            f.write("conda activate rootpy\n")
            f.write("source /opt/geant4-v11.1.0-install/bin/geant4.sh\n")
        # run DNA simulation
        f.write(
            f"time {makerundir(dna_sim_folder)}/tat -PS {test_dir}/out_Atdecay_{n_string}_spacing_{s_string}um.bin -out {test_dir}/out_AtDNA_{n_string}_spacing_{s_string}um.root -sugar {sugarFile} -histone {histoneFile} -seed {seed} \n"
        )

        

        # combinedFilename = os.path.join(test_dir, f"alpha_{angle}deg_{seed}_18.5.root")

        # output = "hadd " + combinedFilename

        # for i in range(0, numThread):
        #     toCombine = f" alpha_{angle}deg_{seed}_18.5_t{i}.root"
        #     output += toCombine

        # f.write(output + "\n")

        # f.write(f"if test -f {combinedFilename}; \n")
        # f.write("\t then \n")
        # fileToDelete = os.path.join(test_dir, f"alpha_{angle}deg_{seed}_18.5_t*.root")
        # f.write("rm {}; fi\n".format(fileToDelete))

    filename_clustering = os.path.join(
        test_dir, f"clustering_At_{n_string}_spacing_{s_string}um_{seed}.sh"
    )  # to be run from folder created

    with open(filename_clustering, "w") as f:
        f.write("#!/bin/bash --login\n")

        if slurm:
            f.write(f"#SBATCH --job-name=clustering_At_{n_string}_spacing_{s_string}um\n")
            f.write(
                f"#SBATCH --output=clustering_At_{n_string}_spacing_{s_string}um.out.%J\n"
            )
            f.write(f"#SBATCH --error=clustering_At_{n_string}_spacing_{s_string}um.err.%J\n")
            # maximum job time in D-HH:MM
            f.write(f"#SBATCH --time=0-02:00\n")
            f.write("#SBATCH --nodes=1\n")
            f.write("#SBATCH -p short\n")
            f.write("#SBATCH --ntasks-per-node=1\n")
            f.write(f"#SBATCH --mem 2GB\n")
            f.write("\n")
            f.write("module load apps/root/6.26.00\n")

        f.write("\n")

        # clustering
        clustering_dir = os.path.join(simulation_parent, "Clustering")
        clustering_outdir = os.path.join(test_dir, "clustering_out")
        if not os.path.exists(clustering_outdir):
            os.makedirs(clustering_outdir)

        if slurm:
            f.write("source /user/home/yw18581/.bash_profile\n")
            f.write("source activate clustering\n")
            f.write(
                f"/user/home/yw18581/.conda/envs/clustering/bin/python3.6 {clustering_dir}/run.py --filename {test_dir}/out_AtDNA_{n_string}_spacing_{s_string}um.root --output {clustering_outdir}/out_AtDNA_{n_string}_spacing_{s_string}um_{seed}.csv --sugar {sugarFile}  --sepR True --n_boxes {nboxes} --boxes_per_R {n_div_theta*n_div_Z}\n"
            )
        else:
            f.write("conda activate clustering\n")
            f.write(
                f"python {clustering_dir}/run.py --filename {test_dir}/out_AtDNA_{n_string}_spacing_{s_string}um.root --output {clustering_outdir}/out_AtDNA_{n_string}_spacing_{s_string}um_{seed}.csv --sugar {sugarFile}  --sepR True --n_boxes {nboxes} --boxes_per_R {n_div_theta*n_div_Z}\n"
            )
    filename_runscript = os.path.join(test_dir, f"run_script_At_{n_string}_spacing_{s_string}um_{seed}.sh")
    with open(filename_runscript, "w") as f:
        f.write("#!/bin/bash --login\n")

        if slurm:
            f.write(f"#SBATCH --job-name=run_At{n_string}_sp_{s_string}um\n")
            f.write(f"#SBATCH --output=run_At_{n_string}_spacing_{s_string}um.out.%J\n")
            f.write(f"#SBATCH --error=run_At_{n_string}_spacing_{s_string}um.err.%J\n")
            f.write(f"#SBATCH --time=0-10:00\n")
            f.write("#SBATCH --nodes=1\n")
            f.write("#SBATCH -p short\n")
            f.write("#SBATCH --ntasks-per-node=1\n")
            f.write(f"#SBATCH --mem 2GB\n")
            f.write("\n")
            f.write("source /user/home/yw18581/.bash_profile\n")
            f.write("module load apps/root/6.26.00\n")

        f.write("\n")

        
        if not slurm:
            f.write(f"source {filename_decay}\n")
            f.write(f"source {filename_DNA}\n")
            f.write(f"source {filename_clustering}\n")
        if  slurm:
            f.write(f"sbatch {filename_decay}\n")
            f.write(f"sbatch {filename_DNA}\n")
            f.write(f"sbatch {filename_clustering}\n")
        

print("Done.")