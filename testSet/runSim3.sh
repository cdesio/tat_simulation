#!/bin/bash --login
#SBATCH --job-name=testSet
#SBATCH --output=testSet.out.%J
#SBATCH --error=testSet.err.%J
#SBATCH --time=0-02:00
#SBATCH --nodes=1
#SBATCH -p short
#SBATCH --ntasks-per-node=8
#SBATCH --mem 2GB

module add lang/gcc/9.1.0
module add apps/geant/4.11.0.0

time ../../simulation/build/rbe -mac ../sim3_bp.in -out sim3 -sugar sugarPos_4x4_300nm.csv -histone histonePositions_4x4_300nm.csv -seed 1

module load apps/root/6.26.00

conda activate clustering
python ../../Clustering/run.py --filename sim3.root --output sim3.csv --sugar sugarPos_4x4_300nm.csv


