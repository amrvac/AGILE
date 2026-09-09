title: Getting started

# Installing AGILE {: #install }
(from the repository [README](https://github.com/amrvac/AGILE))

* Install the [uv](https://docs.astral.sh/uv/) package manager, which manages
  the Python dependencies (`fypp`, `fortdepend`, `f90nml`) needed by the
  build system:

        pip install uv
        # or
        curl -LsSf https://astral.sh/uv/install.sh | sh

* Make sure `$AGILE_DIR` points to the repository root folder:

        export AGILE_DIR=/path/to/AGILE

* Install the required Python packages and activate them:

        cd $AGILE_DIR
        uv sync
        source $AGILE_DIR/.venv/bin/activate

* Go into a test, e.g.

        cd $AGILE_DIR/tests/hd/KH3D

* To compile, load the appropriate modules, e.g. on Snellius:

        module purge
        module load 2023
        module load OpenMPI/4.1.5-NVHPC-24.5-CUDA-12.1.1

* Compile with `nvfortran` and OpenACC activated via

        make arch=nvidia OPENACC=1

  For a plain CPU build (`gfortran`/OpenMPI, no GPU toolchain required),
  use `make arch=gnu` instead.

# Running the compiled test {: #first-run }
Once compiled, each test directory has a `./agile` executable and a
parameter file (`agile.par` or `test.par`). Run it with MPI, for example:

    mpirun -np 4 ./agile -i test.par

This is the same pattern the automated test suite uses (see
`tests/test_rules.make`) -- it produces `.dat`/`.log`/`.vtu` output in the
test directory, which can be compared against the reference output already
present in that test's `correct_output/` folder.
