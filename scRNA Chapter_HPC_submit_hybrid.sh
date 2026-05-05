#!/bin/bash -l
#$ -l h_rt=2:00:00
#$ -l mem=64G
#$ -N Hybrid_Cell
#$ -pe smp 4
#$ -wd /home/regmddy/Scratch/scRNASeq
#$ -o hybrid.out
#$ -e hybrid.err

module -f unload compilers mpi gcc-libs

source activate /home/regmddy/ACFS/Programmes/miniconda/envs/sc_env
export PATH="/home/regmddy/ACFS/Programmes/miniconda/envs/sc_env/bin:$PATH"

# Set number of threads
export OMP_NUM_THREADS=4


# Run the Hybrid Python logic
python run_hybrid_model.py

# Run the Individual R plotting
Rscript visualize_individual.R

