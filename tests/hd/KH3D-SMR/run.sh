#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --gpus=1
#SBATCH --partition=gpu_h100
#SBATCH --mem=0
#SBATCH --job-name=NVHPC
#SBATCH --mail-type=ALL
#SBATCH --time=00-00:30:00
#SBATCH -o slurms/single_out-%j # STDOUT 
#SBATCH -e slurms/single_err-%j # STDERR

source $AMRVAC_DIR/.venv/bin/activate

#ulimit -s
#ulimit -s 131072

source load_modules.sh
# my load_modules.sh contains:
# module purge
# module load 2023 OpenMPI/4.1.5-NVHPC-24.5-CUDA-12.1.1

#these do not work presently
#module purge
#module load 2025
#module load OpenMPI/5.0.7-NVHPC-25.3-CUDA-12.8.0
#module load 2023
#module load OpenMPI/4.1.5-NVHPC-23.7-CUDA-12.1.1

#export CUDA_LAUNCH_BLOCKING=1
#export CUDA_MEMCHECK=1  # or use cuda-memcheck tool

###make
make clean-all
make clean
make -j16 arch=nvidia OPENACC=1 
#NOGPUDIRECT=1

#export PGI_ACC_NOTIFY=16

srun ./amrvac

#cuda-memcheck ./amrvac

