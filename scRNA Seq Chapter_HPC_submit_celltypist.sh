#!/bin/bash -l

# 1. Scheduler Resources
#$ -l h_rt=02:00:00            
#$ -l mem=64G                 
#$ -N CellTypist_Annotation    
#$ -pe smp 4                  
#$ -wd /home/regmddy/Scratch/scRNASeq
#$ -o celltypist.out           
#$ -e celltypist.err     

# 2. Load Modules
module purge
module load python/miniconda3

# 3. Activate Environment
source activate /home/regmddy/ACFS/Programmes/miniconda/envs/sc_env

# 4. Critical: Ensure the correct Python and local Models are found
export PATH="/home/regmddy/ACFS/Programmes/miniconda/envs/sc_env/bin:$PATH"
export PYTHONPATH=$PYTHONPATH:$(pwd)

# 5. Run the script
echo "Starting CellTypist annotation at $(date)"
python run_celltypist.py
echo "Finished at $(date)"
