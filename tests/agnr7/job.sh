#!/bin/bash -l
#SBATCH --job-name="bse"
#SBATCH --account="hck"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jiacao@ethz.ch
#SBATCH --time=02:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=256
####SBATCH --constraint=gpu
###SBATCH --partition=debug
#SBATCH --uenv=prgenv-gnu/24.7:v3
#SBATCH --view=modules
####SBATCH --reservation=eurohack24
set -e -u

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=1
export MPICH_GPU_SUPPORT_ENABLED=1
#export OMP_PROC_BIND=true
#export OMP_PLACES=cores
#export NSYS=1
#export NSYS_FILE=bse_dist_${SLURM_JOBID}_${OMP_NUM_THREADS}_numRanks${SLURM_NPROCS}.qdrep


source ~/load_modules.sh
conda activate magi

srun python bse-sparse.py
