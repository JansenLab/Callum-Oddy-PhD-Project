#!/bin/bash -l

#$ -l h_rt=6:00:00              # Request 6 hours runtime
#$ -l mem=64G                   # Request 64GB RAM
#$ -N seurat_processing         # Job name
#$ -wd /home/regmddy/Scratch/scRNASeq
#$ -pe smp 4                    # Request 4 cores
#$ -o seurat_processing.out     # Standard output
#$ -e seurat_processing.err     # Standard error


# Load required modules
module -f unload compilers mpi gcc-libs


source activate /home/regmddy/ACFS/Programmes/miniconda/envs/sc_env
export PATH="/home/regmddy/ACFS/Programmes/miniconda/envs/sc_env/bin:$PATH"


# Set number of threads
export OMP_NUM_THREADS=4

# Print job information
echo "=========================================="
echo "Job started on: $(date)"
echo "Job ID: $JOB_ID"
echo "Running on node: $(hostname)"
echo "Working directory: $(pwd)"
echo "=========================================="
echo ""

# Debugging: Prove we are using the right R
echo "=========================================="
echo "We are running this R executable:"
which R
echo "=========================================="

# Run R script
Rscript seurat_processing.R

echo ""
echo "=========================================="
echo "Job finished on: $(date)"
echo "=========================================="

