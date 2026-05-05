#!/bin/bash -l
#$ -l h_rt=48:00:00
#$ -l mem=80G
#$ -N scvelo_v2
#$ -pe smp 4
#$ -wd /home/regmddy/Scratch/scRNASeq
#$ -o scvelo_v2.out
#$ -e scvelo_v2.err


module -f unload compilers mpi gcc-libs

source activate /home/regmddy/ACFS/Programmes/miniconda/envs/scvelo
export PATH="/home/regmddy/ACFS/Programmes/miniconda/envs/monocle3_env/bin:$PATH"
export OMP_NUM_THREADS=4

echo "starting Python Script"
python scvelo_v2.py

