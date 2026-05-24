#!/bin/bash
#SBATCH -p gpu_a100
#SBATCH -N 1
#SBATCH --gpus-per-node 2
#SBATCH --ntasks-per-node 2
#SBATCH --cpus-per-task 1
#SBATCH --gpu-bind=map_gpu:0,1,2,3
#SBATCH -t 00:12:00
#SBATCH -o out


module purge
module load 2023
module load OpenMPI/4.1.5-NVHPC-24.5-CUDA-12.1.1

source $AGILE_DIR/.venv/bin/activate

make clean-all
make -j8 arch=nvidia OPENACC=1 NOGPUDIRECT=1

srun ./agile -i agile.par
