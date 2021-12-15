#!/bin/bash -l
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 21
#SBATCH --time=1-00:00:00
#SBATCH -p batch

parallel -j 21 'zip -m {}.zip {}' ::: img_data/*/Sequence/*/reads_*.fq
