# WENO-JS5 mathematical tests

Requires Python with NumPy, Fypp on PATH (or in `.venv/bin`), GNU Fortran,
and MPI for the face/solver tests. Run from the repository root:

```sh
python tests/weno_js5/test_kernel.py
export AGILE_DIR=$PWD
make -C tests/weno_js5/hd -j4
python tests/weno_js5/test_faces.py
```

The kernel reference solves cell-average polynomial constraints and integrates
squared derivatives independently of the Fortran formulas. Face tests compile
a temporary public copy of the production module, call the actual six-point
single-face implementation, and check both interfaces, nonlinear state
conversion and donor replacement. Production visibility is unchanged.

For both layouts, MPI/AMR, configuration and spatial-order tests, see
[MP5 integration tests](../mp5/README.md). No generated output is a baseline
that should be refreshed automatically.
