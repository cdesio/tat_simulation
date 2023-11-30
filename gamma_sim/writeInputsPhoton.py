'''
Create all files required to run photon simulation(s) with set seeds. Photon and secondaries simulation run separately so that multithreading can be used on secondary simulation
'''

import numpy as np
import os
import random

# parameters to set
folder = "testing"
numThread = 1
slurm = True
seeds = [random.randint(1, 10000) for a in range(1)]
histoneFile = "../../../../geometryFiles/histonePositions_4x4_300nm.csv"
sugarFile = "../../../../geometryFiles/sugarPos_4x4_300nm.csv"
timePhoton = "0-01:00"
memPhoton = 1
timeSec = "0-01:00"
memSec = 2
targetSize = 300 #nanometers
gpsRadius = 3 #millimeters
printProgress = 100000
numPhotons = 10000000

####

os.makedirs(f"./results/{folder}")
for seed in seeds:
    filename = f"./results/{folder}/runPhotonScript_{seed}.sh" #to be run from folder created
    if slurm:
            print(f"sbatch runPhotonScript_{seed}.sh")
    with open(filename , 'w') as f:
        f.write("#!/bin/bash --login\n")

        if slurm:
            f.write(f"#SBATCH --job-name=photon_{seed}\n")
            f.write(f"#SBATCH --output=photon_{seed}.out.%J\n")
            f.write(f"#SBATCH --error=photon_{seed}.err.%J\n")
            # maximum job time in D-HH:MM
            f.write(f"#SBATCH --time={timePhoton}\n")
            f.write("#SBATCH --nodes=1\n")
            f.write("#SBATCH -p short\n")
            f.write("#SBATCH --ntasks-per-node=1\n") #photon sim is single threaded
            f.write(f"#SBATCH --mem {memPhoton}GB\n")
            f.write("\n")
            f.write("module load apps/geant/4.11.0.0\n")
        
        f.write("\n")

        # photon simulation
        f.write(f"time ../../photonSim -mac Co60.in -out photon_{seed} -seed {seed}\n")

        f.write("\n")

        # secondary particle simulation
        f.write(f"time ../../../../simulation/build/rbe -mac secondary.in -PS photon_{seed}.bin -out secondary_{seed} -sugar {sugarFile} -histone {histoneFile} -seed {seed}\n")
        


        if slurm:
            f.write("module load apps/root/6.26.00\n")
        
        f.write("\n")

        # clustering
        f.write("conda activate clustering\n")
        f.write(f"python ../../../../Clustering/run.py --filename secondary_{seed}.root --output photon_{seed}.csv --sugar {sugarFile} --filenamePhoton photon_{seed}.root \n")


#  write mac files
filename = f"./results/{folder}/Co60.in" 
with open(filename , 'w') as f:
    f.write("/run/verbose 2\n")
    f.write("/control/verbose 2\n")
    f.write(f"/det/setSize {targetSize} nm\n")
    f.write("/run/initialize\n")

    f.write("/gps/particle gamma\n")
    f.write("/gps/pos/type Surface\n")
    f.write("/gps/pos/shape Sphere\n")
    f.write("/gps/pos/centre 0 0 0\n")
    f.write(f"/gps/pos/radius {gpsRadius} mm\n")
    f.write("/gps/ang/type cos\n")
    maxTheta = np.arcsin((3**0.5*targetSize/2)/(gpsRadius*1e6))
    f.write(f"/gps/ang/maxtheta {maxTheta} rad\n")
    f.write("/gps/ene/type Mono\n")
    f.write(f"/run/printProgress {printProgress}\n")
    f.write(f"/run/beamOn {numPhotons}\n")

filename = f"./results/{folder}/secondary.in"
with open(filename , 'w') as f:
    f.write("/run/verbose 2\n")
    f.write("/control/verbose 2\n")
    f.write(f"/det/setSize {targetSize} nm\n")
    f.write(f"/run/numberOfThreads {numThread}\n")




