#!/bin/bash --login
#SBATCH --job-name=photon_9124
#SBATCH --output=photon_9124.out.%J
#SBATCH --error=photon_9124.err.%J
#SBATCH --time=0-01:00
#SBATCH --nodes=1
#SBATCH -p short
#SBATCH --ntasks-per-node=1
#SBATCH --mem 1GB

module load apps/geant/4.11.0.0

time ../../photonSim -mac Co60.in -out photon_9124 -seed 9124
sbatch runSecondaryScript_9124.sh
