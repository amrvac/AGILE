# MP5 and high-order FV integration tests

Requires Python with NumPy and f90nml, Fypp on PATH (or in `.venv/bin`), GNU
Fortran and MPI. Run from the repository root:

```sh
export AGILE_DIR=$PWD
python tests/mp5/test_kernel.py
make -C tests/weno_js5/hd -j4
python tests/mp5/test_faces.py
python tests/mp5/test_layouts.py --output /tmp/agile-layout-check
python tests/mp5/test_layout_accuracy.py --output /tmp/agile-spatial-check
python tests/mp5/run_validation.py --layout per_cell --output /tmp/agile-per-cell
python tests/mp5/run_validation.py --layout shared --output /tmp/agile-shared
```

On macOS, some MPI installations also require `HWLOC_COMPONENTS=-opencl`.
Use fresh output directories. `--quick` omits the longer dynamic-AMR and shock
checks. A debug build can be made in the same case with `DEBUG=1`; GNU flags
include array bounds checks and floating-point exception traps.

The scalar reference derives the fifth-order candidate from cell-integral
constraints and implements limiting through interval clipping. It covers both
zero and default tolerance. Face tests call the actual single-face interface
in a temporary public build, including donor replacement and shared faces.
`test_layouts.py` compares four methods, three axes and three fluxes bitwise;
`--previous /path/to/old/agile` additionally compares a prior executable.
The spatial test extracts the RHS using one Euler step and checks fifth order
for WENO/MP5, with bitwise agreement between layouts.

The solver suite covers halo widths 3/4, MPI, reflection, mixed limiters,
static/dynamic AMR and invalid configuration. Shock smoke checks use a short-run
constant-boundary estimate, not a general numerical-flux conservation audit.
The optional `--convergence` retains a strict finest-step halving criterion;
MP5 can fail that criterion near roundoff even when spatial order passes.

Production diagnostic counters have been removed. Unmeasured fallback counts
are null, not zero. Donor replacement itself is tested through returned states.
Full GPU runtime acceptance and performance comparison are separate tasks.
