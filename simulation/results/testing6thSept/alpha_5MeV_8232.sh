#!/bin/bash --login
#SBATCH --job-name=alpha_5MeV_8232
#SBATCH --output=alpha_5MeV_8232.out.%J
#SBATCH --error=alpha_5MeV_8232.err.%J
#SBATCH --time=0-10:00
#SBATCH --nodes=1
#SBATCH -p short
#SBATCH --ntasks-per-node=4
#SBATCH --mem 2GB

module load apps/geant/4.11.0.0
time ../../rbe -mac alpha_5MeV.in -out alpha_5MeV_8232 -sugar ../../../../geometryFiles/sugarPos_4x4_300nm.csv -histone ../../../../geometryFiles/histonePositions_4x4_300nm.csv -seed 8232
module load apps/root/6.26.00
hadd alpha_5MeV_8232.root alpha_5MeV_8232_t0.root alpha_5MeV_8232_t1.root alpha_5MeV_8232_t2.root alpha_5MeV_8232_t3.root
rm alpha_5MeV_8232_t0.root
rm alpha_5MeV_8232_t1.root
rm alpha_5MeV_8232_t2.root
rm alpha_5MeV_8232_t3.root
sbatch clustering_alpha_5MeV_8232.sh
