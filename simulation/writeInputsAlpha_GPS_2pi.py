"""
Create all files required to run alpha simulation(s) with set seeds and energies. simulation & clustering run separately so that multithreading can be used on simulation
"""

import numpy as np
import os
import random

# parameters to set
folder = "test_2GPS_2pi"
numThread = 4
slurm = True
seed = random.randint(1, 10000)
angles = [0, 22.5, 45, 67.5, 90]
shifts_x = [10.50, 9.70, 7.42, 4.02, 0.00]
shifts_y = [0.00, 4.02, 7.42, 9.70, 10.50]
shifts_nm_x = np.array(shifts_x) * 1000
shifts_nm_y = np.array(shifts_y) * 1000

if slurm:
    simulation_parent = os.path.join(
        "/", "user", "work", "yw18581", "DaRT", "TAT", "tat_simulation_dna_vertical"
    )
else:
    simulation_parent = os.path.join(
        "/", "Users", "yw18581", "UoB", "DaRT", "TAT", "tat_simulation_DNA"
    )

histoneFile = os.path.join(
    simulation_parent, "geometryFiles", "histonePositions_4x4_300nm.csv"
)
sugarFile = os.path.join(simulation_parent, "geometryFiles", "sugarPos_4x4_300nm.csv")

time = "0-10:00"
mem = 10
targetSize = 300  # nanometers
gpsRadius = 10
gpsHalfZ = 20  # micrometers
printProgress = 100000
numAlpha = 10000000

####
simulation_folder = os.path.join(simulation_parent, "simulation")
run_dir = os.path.join(simulation_folder, "build")

#### create folder for current test
test_dir = os.path.join(simulation_folder, f"output/{folder}")
if not os.path.exists(test_dir):
    os.makedirs(test_dir)


#### create one file per direction/shift
for angle, sh_x, sh_y, sh_x_nm, sh_y_nm in zip(
    angles, shifts_x, shifts_y, shifts_nm_x, shifts_nm_y
):

    filename = os.path.join(
        test_dir, f"run_alpha_{angle}deg_{seed}.sh"
    )  # to be run from folder created
    if slurm:
        print(f"sbatch run_alpha_{angle}deg_{seed}.sh")
    with open(filename, "w") as f:
        f.write("#!/bin/bash --login\n")
        if slurm:
            f.write(f"#SBATCH --job-name=alpha_{angle}deg_{seed}\n")
            f.write(f"#SBATCH --output=alpha_{angle}deg_{seed}.out.%J\n")
            f.write(f"#SBATCH --error=alpha_{angle}deg_{seed}.err.%J\n")
            # maximum job time in D-HH:MM
            f.write(f"#SBATCH --time={time}\n")
            f.write("#SBATCH --nodes=1\n")
            f.write("#SBATCH -p short\n")
            f.write(f"#SBATCH --ntasks-per-node={numThread}\n")
            f.write(f"#SBATCH --mem {mem}GB\n")
            f.write("\n")
            f.write("module load apps/geant/4.11.0.0\n")
            f.write("module load apps/root/6.26.00\n")
            f.write("source /user/home/yw18581/.bash_profile\n")
            f.write("source activate dart\n")

        else:
            f.write("conda activate dart\n")
        # run simulation
        f.write(
            f"time {run_dir}/rbe -mac {test_dir}/alpha_{angle}deg_{seed}.in -out {test_dir}/alpha_{angle}deg_{seed}.root -sugar {sugarFile} -histone {histoneFile} -seed {seed}\n"
        )

        combinedFilename = os.path.join(test_dir, f"alpha_{angle}deg_{seed}.root")

        output = "hadd " + combinedFilename

        for i in range(0, numThread):
            toCombine = f" alpha_{angle}deg_{seed}_t{i}.root"
            output += toCombine

        f.write(output + "\n")

        f.write(f"if test -f {combinedFilename}; \n")
        f.write("\t then \n")
        fileToDelete = os.path.join(test_dir, f"alpha_{angle}deg_{seed}_t*.root")
        f.write("rm {}; fi\n".format(fileToDelete))

        if not slurm:
            f.write(f"chmod 755 clustering_alpha_{angle}deg_{seed}.sh\n")
            f.write(f"./clustering_alpha_{angle}deg_{seed}.sh\n")
        if slurm:
            f.write(f"sbatch clustering_alpha_{angle}deg_{seed}.sh\n")

    filename = os.path.join(
        test_dir, f"clustering_alpha_{angle}deg_{seed}.sh"
    )  # to be run from folder created

    with open(filename, "w") as f:
        f.write("#!/bin/bash --login\n")

        if slurm:
            f.write(f"#SBATCH --job-name=clustering_alpha_{angle}deg_{seed}\n")
            f.write(f"#SBATCH --output=clustering_alpha_{angle}deg_{seed}.out.%J\n")
            f.write(f"#SBATCH --error=clustering_alpha_{angle}deg_{seed}.err.%J\n")
            # maximum job time in D-HH:MM
            f.write(f"#SBATCH --time=0-01:00\n")
            f.write("#SBATCH --nodes=1\n")
            f.write("#SBATCH -p short\n")
            f.write("#SBATCH --ntasks-per-node=1\n")
            f.write(f"#SBATCH --mem 2GB\n")
            f.write("\n")
            f.write("module load apps/root/6.26.00\n")

        f.write("\n")

        # clustering
        clustering_dir = os.path.join(simulation_parent, "Clustering")
        if slurm:
            f.write("source /user/home/yw18581/.bash_profile\n")
            f.write("source activate clustering\n")
            f.write(
                f"/user/home/yw18581/.conda/envs/clustering/bin/python3.6 {clustering_dir}/run.py --filename {test_dir}/alpha_{angle}deg_{seed}.root --output {test_dir}/alpha_{angle}deg_{seed}.csv --sugar {sugarFile} --shift {sh_x_nm} {sh_y_nm} 0.00\n"
            )
        else:
            f.write("conda activate clustering\n")
            f.write(
                f"python {clustering_dir}/run.py --filename alpha_{angle}deg_{seed}.root --output alpha_{angle}deg_{seed}.csv --sugar {sugarFile} --shift {sh_x_nm} {sh_y_nm} 0.00\n"
            )

    #  write mac file

    filename = os.path.join(
        test_dir, f"alpha_{angle}deg_{seed}.in"
    )  # to be run from folder created
    with open(filename, "w") as f:
        f.write("/run/verbose 2\n")
        f.write("/control/verbose 2\n")
        f.write(f"/det/displ_X {sh_x} um\n")
        f.write(f"/det/displ_Y {sh_y} um\n")
        f.write(f"/run/numberOfThreads {numThread}\n")
        f.write("/run/initialize\n")
        f.write("\n")
        for en in [5.9, 7.5]:
            f.write("/gps/source/add 1.\n")
            f.write("/gps/pos/type Surface\n")
            f.write("/gps/pos/shape Cylinder\n")
            f.write("/gps/pos/centre 0 0 0 um\n")
            f.write(f"/gps/pos/radius {gpsRadius} um\n")
            f.write(f"/gps/pos/halfz {gpsHalfZ} um\n")
            f.write("\n")
            f.write("/gps/ang/type cos\n")
            f.write("/gps/particle alpha\n")
            f.write("/gps/ene/type Mono\n")
            f.write(f"/gps/ene/mono {en} MeV\n")
            # f.write(f"/gps/ang/mintheta {0} deg\n")
            # f.write(f"/gps/ang/maxtheta {360} deg\n")
            f.write("\n")
        f.write(f"/run/printProgress {int(numAlpha/10)}\n")
        f.write(f"/run/beamOn {numAlpha}\n")
