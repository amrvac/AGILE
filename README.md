[![tests (CPU)](https://github.com/amrvac/AGILE/actions/workflows/tests.yml/badge.svg)](https://github.com/amrvac/AGILE/actions/workflows/tests.yml)
[![tests (GPU)](https://github.com/amrvac/AGILE/actions/workflows/gpu.yml/badge.svg)](https://github.com/amrvac/AGILE/actions/workflows/gpu.yml)

# AGILE

This is the public development repository of AGILE, a GPU enabled fork of MPI-AMRVAC. 

To install, follow these steps:
- install the uv package manager `pip install uv` or `curl -LsSf https://astral.sh/uv/install.sh | sh`
- make sure `$AGILE_DIR` points to the repository root folder
- install the required python packages: `cd $AGILE_DIR` and run `uv sync` and activate them `source $AGILE_DIR/.venv/bin/activate`
- go into a test, e.g. `cd $AGILE_DIR/tests/hd/KH3D`
- to compile, load the appropriate modules, e.g. on snellius:
```
module purge
module load 2023
module load OpenMPI/4.1.5-NVHPC-24.5-CUDA-12.1.1
```
- compile with nvfortran and activated OPENACC via `make arch=nvidia OPENACC=1`

## Currently supported features on master
- Cartesian grids
- Physics modules: hd, mhd [glm], ffhd
- Source terms (gravity, radiative cooling, hyperbolic thermal conduction, user defined) and boundary conditions (`symm, asymm, cont` etc. but also `special`)
- Multi-GPU (MPI)
- Uniform grid, static mesh refinement (SMR) and adaptive mesh refinement (AMR)
