#!/bin/bash -l
#$ -l h_rt=48:00:00
#$ -l mem=80G
#$ -N cytotrace_anndata
#$ -pe smp 4
#$ -wd /home/regmddy/Scratch/scRNASeq
#$ -o cytotrace_export_phate.out
#$ -e cytotrace_export_phate.err


module -f unload compilers mpi gcc-libs

source activate /home/regmddy/ACFS/Programmes/miniconda/envs/monocle3_env
export PATH="/home/regmddy/ACFS/Programmes/miniconda/envs/monocle3_env/bin:$PATH"
export OMP_NUM_THREADS=4


echo "Starting R Script"
Rscript cytotrace_to_anndata_for_phate.R

echo "starting Python Script"
python cytotrace_to_anndata_for_phate.py
