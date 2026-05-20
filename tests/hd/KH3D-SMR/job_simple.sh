#!/bin/bash
#SBATCH -p standard-g
#SBATCH -N 1
#SBATCH --account=project_465002953 
#SBATCH --gpus-per-node 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH -t 00:30:00
#SBATCH -o job_mpi_lumi.out
#
# Note on partitions: standard-g should be used for full nodes, small-g for smaller allocations

module load LUMI/25.03 partition/G cray-python/3.11.7 PrgEnv-cray lumi-CrayPath

## See https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/lumig-job/
#
## Enable GPU-aware MPI support
#export MPICH_GPU_SUPPORT_ENABLED=1
## Enable GPU-based NIC selection
#export MPICH_OFI_NIC_POLICY=GPU
#
## Ensure proper GPU to CPU binding, only works when reserving full nodes.
## If not, remove this part and run "srun ./amrvac"
#cat << EOF > select_gpu
##!/bin/bash
#
#export ROCR_VISIBLE_DEVICES=\$SLURM_LOCALID
#exec \$*
#EOF
#
#chmod +x ./select_gpu
#
## GPU 0 is bound to CPU 49, GPU 1 to CPU 57, etc.
#CPU_BIND="map_cpu:49,57,17,25,1,9,33,41"

srun ./amrvac
