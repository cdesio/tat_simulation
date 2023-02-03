'''
Create all files required to run alpha simulation(s) with set seeds and energies. simulation & clustering run separately so that multithreading can be used on simulation
'''

import numpy as np
import os
import random

# parameters to set
folder = "testing6thSept"
numThread = 4
slurm = True
seeds = [random.randint(1, 10000) for a in range(1)]
energies = ["5","0p5"]
histoneFile = "../../../../geometryFiles/histonePositions_4x4_300nm.csv"
sugarFile = "../../../../geometryFiles/sugarPos_4x4_300nm.csv"
time = "0-10:00"
mem = 2
targetSize = 300 #nanometers
gpsRadius = 1 #micrometers
printProgress = 1
numAlpha = 4

####

os.makedirs(f"./results/{folder}")
for energy in energies:
    for seed in seeds:
        filename = f"results/{folder}/alpha_{energy}MeV_{seed}.sh" #to be run from folder created
        if slurm:
            print(f"sbatch alpha_{energy}MeV_{seed}.sh")
        with open(filename , 'w') as f:
            f.write("#!/bin/bash --login\n")
            if slurm:
                f.write(f"#SBATCH --job-name=alpha_{energy}MeV_{seed}\n")
                f.write(f"#SBATCH --output=alpha_{energy}MeV_{seed}.out.%J\n")
                f.write(f"#SBATCH --error=alpha_{energy}MeV_{seed}.err.%J\n")
                # maximum job time in D-HH:MM
                f.write(f"#SBATCH --time={time}\n")
                f.write("#SBATCH --nodes=1\n")
                f.write("#SBATCH -p short\n")
                f.write(f"#SBATCH --ntasks-per-node={numThread}\n") 
                f.write(f"#SBATCH --mem {mem}GB\n")
                f.write("\n")
                f.write("module load apps/geant/4.11.0.0\n")

            # secondary particle simulation
            f.write(f"time ../../rbe -mac alpha_{energy}MeV.in -out alpha_{energy}MeV_{seed} -sugar {sugarFile} -histone {histoneFile} -seed {seed}\n")
            
            if not slurm:
                f.write(f"chmod 755 clustering_alpha_{energy}MeV_{seed}.sh\n")
                f.write(f"./clustering_alpha_{energy}MeV_{seed}.sh\n")
            if  slurm:
                f.write(f"sbatch clustering_alpha_{energy}MeV_{seed}.sh\n")

        filename = f"./results/{folder}/clustering_alpha_{energy}MeV_{seed}.sh" #to be run from folder created

        with open(filename , 'w') as f:
            f.write("#!/bin/bash --login\n")

            if slurm:
                f.write(f"#SBATCH --job-name=clustering_alpha_{energy}MeV_{seed}\n")
                f.write(f"#SBATCH --output=clustering_alpha_{energy}MeV_{seed}.out.%J\n")
                f.write(f"#SBATCH --error=clustering_alpha_{energy}MeV_{seed}.err.%J\n")
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
            f.write("conda activate clustering\n")
            f.write(f"python ../../../../Clustering/run.py --filename alpha_{energy}MeV_{seed}.root --output alpha_{energy}MeV_{seed}.csv --sugar {sugarFile} \n")

    #  write mac file

    filename = f"./results/{folder}/alpha_{energy}MeV.in" #to be run from folder created
    with open(filename , 'w') as f:
        f.write("/run/verbose 2\n")
        f.write("/control/verbose 2\n")
        f.write(f"/det/setSize {targetSize} nm\n")
        f.write(f"/run/numberOfThreads {numThread}\n")
        f.write("/run/initialize\n")

        f.write("/gps/particle alpha\n")
        f.write("/gps/pos/type Surface\n")
        f.write("/gps/pos/shape Sphere\n")
        f.write("/gps/pos/centre 0 0 0\n")
        f.write(f"/gps/pos/radius {gpsRadius} um\n")
        f.write("/gps/ene/type Mono\n")
        en = energy.replace("p",".")
        f.write(f"/gps/ene/mono {en} MeV\n")
        f.write("/gps/ang/type cos\n")
        maxTheta = np.arcsin((3**0.5*targetSize/2)/(gpsRadius*1e3))
        f.write(f"/gps/ang/maxtheta {maxTheta} rad\n")
        f.write(f"/scheduler/endTime 5 nanosecond\n")
        f.write(f"/run/printProgress {printProgress}\n")
        f.write(f"/run/beamOn {numAlpha}\n")

