#!/bin/bash --login
#SBATCH --job-name=clustering_alpha_0p5MeV_8232
#SBATCH --output=clustering_alpha_0p5MeV_8232.out.%J
#SBATCH --error=clustering_alpha_0p5MeV_8232.err.%J
#SBATCH --time=0-01:00
#SBATCH --nodes=1
#SBATCH -p short
#SBATCH --ntasks-per-node=1
#SBATCH --mem 2GB

module load apps/root/6.26.00

conda activate clustering
python ../../../../Clustering/run.py --filename alpha_0p5MeV_8232.root --output alpha_0p5MeV_8232.csv --sugar ../../../../geometryFiles/sugarPos_4x4_300nm.csv 
