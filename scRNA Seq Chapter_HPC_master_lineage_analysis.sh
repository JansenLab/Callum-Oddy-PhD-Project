#!/bin/bash -l
#$ -l h_rt=4:00:00
#$ -l mem=80G
#$ -N Pub_Viz
#$ -pe smp 4
#$ -wd /home/regmddy/Scratch/scRNASeq
#$ -o master_lineage.out
#$ -e master_lineage.err


module -f unload compilers mpi gcc-libs

source activate /home/regmddy/ACFS/Programmes/miniconda/envs/sc_env
export PATH="/home/regmddy/ACFS/Programmes/miniconda/envs/sc_env/bin:$PATH"
export OMP_NUM_THREADS=4

R CMD BATCH --no-save master_lineage_analysis.R

