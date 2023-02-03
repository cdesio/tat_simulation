#!/bin/bash --login
#SBATCH --job-name=secondary_9124
#SBATCH --output=secondary_9124.out.%J
#SBATCH --error=secondary_9124.err.%J
#SBATCH --time=0-01:00
#SBATCH --nodes=1
#SBATCH -p short
#SBATCH --ntasks-per-node=1
#SBATCH --mem 2GB

module load apps/geant/4.11.0.0
time ../../../../simulation/build/rbe -mac secondary.in -PS photon_9124.bin -out secondary_9124 -sugar ../../../../geometryFiles/sugarPos_4x4_300nm.csv -histone ../../../../geometryFiles/histonePositions_4x4_300nm.csv -seed 9124
module load apps/root/6.26.00
hadd secondary_9124.root secondary_9124_t0.root
rm secondary_9124_t0.root
sbatch clustering_9124.sh
