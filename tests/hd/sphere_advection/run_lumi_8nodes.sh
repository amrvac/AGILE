#!/bin/bash
#SBATCH -p dev-g
#SBATCH -N 1
#SBATCH --account project_465002953
#SBATCH --gpus-per-node 8
#SBATCH --ntasks-per-node 8
#SBATCH --cpus-per-task 1
#SBATCH -t 00:20:00
#SBATCH -o slurms/job_mpi_lumi-%j.out

set -euo pipefail

#module load LUMI/25.03 partition/G cray-python/3.11.7 PrgEnv-cray lumi-CrayPath
module load LUMI/25.09 partition/G cray-python/3.11.7 PrgEnv-cray lumi-CrayPath #gdb4hpc rocm

export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_OFI_NIC_POLICY=GPU

# --- Where you submitted from (contains amrvac + parfile) ---
#SUBMIT_DIR="/users/kellyadr/AGILE/AGILE-experimental/tests/hd/sphere_advection"
SUBMIT_DIR="/users/jessevos/agile_lumi/AGILE-experimental/tests/hd/sphere_advection"

# --- Make a per-job run dir in scratch ---
SCRATCH_BASE="/scratch/project_465002541"
RUN_DIR="${SCRATCH_BASE}/jvos/${SLURM_JOB_ID}"
mkdir -p "${RUN_DIR}"

# --- Stage inputs ---
# Put your actual parfile name here
PARFILE="amrvac.par"

cp -a "${SUBMIT_DIR}/amrvac" "${RUN_DIR}/"
cp -a "${SUBMIT_DIR}/${PARFILE}" "${RUN_DIR}/"

# Optional: if your parfile expects extra files (tables, cooling curves, etc.)
# cp -a "${SUBMIT_DIR}/additional_data.txt" "${RUN_DIR}/"

cd "${RUN_DIR}"
mkdir -p output

## --- GPU selection wrapper ---
#cat <<'EOF' > select_gpu
##!/bin/bash
#export ROCR_VISIBLE_DEVICES=${SLURM_LOCALID}
#exec "$@"
#EOF
#chmod +x ./select_gpu
#
## GPU 0 is bound to CPU 49, GPU 1 to CPU 57, etc.
#CPU_BIND="map_cpu:49,57,17,25,1,9,33,41"

# --- Run ---
#srun -l --cpu-bind=${CPU_BIND} ./select_gpu ./amrvac
srun -l ./amrvac

# --- Copy results back to submit dir (or wherever you want) ---
# rsync -av --info=stats2 "${RUN_DIR}/output/" "${SUBMIT_DIR}/output_${SLURM_JOB_ID}/"

echo "Done. Run dir: ${RUN_DIR}"
