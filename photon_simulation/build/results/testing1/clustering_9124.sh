#!/bin/bash --login
#SBATCH --job-name=clustering_9124
#SBATCH --output=clustering_9124.out.%J
#SBATCH --error=clustering_9124.err.%J
#SBATCH --time=0-01:00
#SBATCH --nodes=1
#SBATCH -p short
#SBATCH --ntasks-per-node=1
#SBATCH --mem 2GB

module load apps/root/6.26.00

conda activate clustering
python ../../../../Clustering/run.py --filename secondary_9124.root --output photon_9124.csv --sugar ../../../../geometryFiles/sugarPos_4x4_300nm.csv --filenamePhoton photon_9124.root 
