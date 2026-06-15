#!/bin/bash -l
#SBATCH --job-name=OT3Dbig
#SBATCH --output=OT3Dbig.o%j
#SBATCH --error=OT3Dbig.e%j
#SBATCH --account=project_465002541
#SBATCH --partition=standard-g
#SBATCH --nodes=64
#SBATCH --gpus-per-node=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --mem=480G

module --force purge
module load LUMI/25.09 partition/G cray-python/3.11.7 PrgEnv-cray lumi-CrayPath

AGILE_DIR="$HOME/FORT/AGILE"
export AGILE_DIR
echo $AGILE_DIR

cd /flash/project_465002541/keppensr/OT3D

# See https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/lumig-job/
# Enable GPU-aware MPI support
export MPICH_GPU_SUPPORT_ENABLED=1
# Enable GPU-based NIC selection
export MPICH_OFI_NIC_POLICY=GPU

# Ensure proper GPU to CPU binding, only works when reserving full nodes.
# If not, remove this part and run "srun ./amrvac"
cat << EOF > select_gpu
#!/bin/bash

export ROCR_VISIBLE_DEVICES=\$SLURM_LOCALID
exec \$*
EOF

chmod +x ./select_gpu

# GPU 0 is bound to CPU 49, GPU 1 to CPU 57, etc.
CPU_BIND="map_cpu:49,57,17,25,1,9,33,41"

srun --cpu-bind=${CPU_BIND} ./select_gpu ./agile_ot3d -i agile_ot3d_4l.par 
rm -f ./select_gpu
