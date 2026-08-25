# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

AGILE is a GPU-enabled fork of MPI-AMRVAC, a Fortran finite-volume code for
solving hyperbolic PDEs (HD, MHD, FFHD, SRHD) with adaptive mesh refinement.
Master supports 3D Cartesian and 3D spherical `(r, theta, phi)` grids; see
"Coordinate systems" below for the limits of the spherical support. Domain
decomposition and GPU offload use MPI + OpenACC (`!$acc` directives throughout
the hot loops).

Source is written as `.fpp` files (Fortran + `fypp` preprocessor directives,
e.g. `#:if PHYS == 'hd'`) which get preprocessed into `.f90` before
compilation. Unlike upstream MPI-AMRVAC, this fork does **not** use the old
VACPP/LASY dimension-independent `.t` preprocessor (`doc/vacpp.md` and much of
`doc/*.md` describe upstream MPI-AMRVAC and are only partially applicable
here — e.g. `setup.pl` and `arch/*.arch` files referenced there don't exist in
this repo). Loop bounds are written out explicitly per-dimension (e.g.
`ixGmin1,ixGmin2,ixGmin3`) directly in the `.fpp` source.

## Environment setup

```bash
pip install uv                       # or curl -LsSf https://astral.sh/uv/install.sh | sh
export AGILE_DIR=/path/to/this/repo  # must point at the repo root
cd $AGILE_DIR && uv sync
source $AGILE_DIR/.venv/bin/activate
```

`AGILE_DIR` and an activated venv (providing `fypp`, `fortdepend`, `f90nml`)
are required for every build. Requires Python >=3.13.

## Building

Each test/run case is its own directory with a `Makefile` and an `agile.par`
parameter file. Build from inside such a directory:

```bash
cd tests/hd/KH3D
make arch=gnu              # compiles and links ./agile executable
make arch=gnu OPENACC=1    # enable GPU/OpenACC build
make arch=gnu DEBUG=1      # debug flags (-g -O0 -fcheck=all, etc.)
make clean                 # remove this case's build products
```

- `arch` selects a file from `arch/*.mk` (`gnu`, `ifx`, `nvidia`, `cray`,
  `llvm`) which sets the compiler/linker and flags. Default is `gnu`
  (`mpif90`/gfortran).
- Physics-relevant compile-time options (which physics module, the coordinate
  system, tracers, gravity, cooling, thermal conduction, etc.) are declared in
  `agile.par` and turned into `config.mk`/fypp defines by
  `make/config_reader.py`, validated
  against `make/config_schema.toml`. Editing `agile.par` and rerunning `make`
  regenerates `config.mk` automatically — don't hand-edit `config.mk`.
- Builds are cached under `build/<arch>-<hash>/`, keyed by a hash of the
  enabled compile-time flags; `build/latest` symlinks to the most recent one.
  `make clean-all` wipes the whole `build/` cache.
- The actual make logic lives in `make/*.mk`, included in numeric order by
  `make/include-all.mk` (00-prelude, 05-find-agile, 10-select-arch,
  15-read-config, 20-build-dir, 25-store-flags, 30-fypp, 35-dependencies,
  40-compile, 50-link, 90-clean).

## Running tests

Tests live under `tests/{hd,mhd,ffhd,srhd}/<case>/` and each case has a
`test.make` including `tests/test_rules.make`, which builds the case fresh
and runs it (via `mpirun`/`srun`), then numerically diffs the log output
against `tests/<physics>/<case>/correct_output/` using
`tools/fortran/compare_logs` (tolerances `1.0e-5` relative / `1.0e-8`
absolute).

```bash
cd tests
make hd            # run all hd/* test cases
make mhd           # run all mhd/* test cases
make ffhd          # run all ffhd/* test cases
make -j4 hd         # run in parallel
make -s hd          # suppress "Entering/leaving directory" noise
make clean          # clean all test build artifacts
```

To run a single case's test directly:

```bash
cd tests/hd/KH3D
make -f test.make arch=gnu
```

This is also what CI does, via `tests/test_runner.sh <target>` (used as
`arch=gnu ./test_runner.sh hd`), which fails if the build fails or if
`**FAILED**` appears in the output log.

Add a new test case by copying an existing case directory (`mod_usr.fpp`,
`agile.par`/`test.par`, `test.make`, `correct_output/`) and registering it in
`tests/Makefile` under the appropriate `*_DIRS` variable.

## Coordinate systems

The coordinate system is a **compile-time** choice, because the finite-volume
kernels are fypp-specialized for it and the Cartesian kernel must stay free of
curvilinear bookkeeping. Set `geometry` in `&meshlist` of `agile.par`:

```fortran
 &meshlist
   geometry = 'spherical'   ! or 'Cartesian' (the default)
 /
```

That becomes the `GEOM` fypp define via `make/config_schema.toml`, and it is
the only place the coordinate system needs stating. `main()` in
`src/agile.fpp` calls `set_coordinate_system_from_config()` before `usr_init`,
so a case does **not** have to call `set_coordinate_system` itself. A case may
still do so — the call is idempotent — but the name then has to agree with
`GEOM`, and `set_coordinate_system` calls `mpistop` on a mismatch rather than
silently running the wrong kernel.

Deriving it centrally matters because `set_coordinate_system` is the sole
writer of `coordinate`, `ndir`, `r_`, `phi_` and `z_`. When it was left to
each case, forgetting the call left `coordinate` at its `-1` default, and the
failure was silent rather than loud: `read_par_files` skips the conversion of
the angular domain bounds, and the `select case (coordinate)` in
`get_surface_area` and in the cell-volume fill match nothing, so `surfaceC`
and `dvolume` are never written and the run completes with garbage.
`tests/hd/spherical_uniform_flow` omits the call and
`tests/hd/spherical_blast` keeps it, so both paths stay covered.

Two things to know when writing a spherical case:

- **Angular coordinates in `agile.par` are given in units of 2*pi** (an
  MPI-AMRVAC convention, applied in `src/io/mod_input_output.fpp`). So
  `xprobmin2 = 0.125d0`, `xprobmax2 = 0.375d0` means `theta` from `pi/4` to
  `3*pi/4`.
- `ndim` is a compile-time `parameter` of 3, so spherical means the full
  3D `(r, theta, phi)` system with `ndir = 3`.

How it works: in a curvilinear build the finite-volume update in
`src/mod_finite_volume.fpp` divides face fluxes by the cell volume weighted by
the face areas (`ps(n)%surfaceC`, `ps(n)%dvolume`, filled per block by
`fillgeo`/`get_surface_area` and copied to the device in
`src/amr/mod_amr_solution_node.fpp`), instead of by the uniform `rnode` cell
spacing. The leftover curvature terms of the momentum equations are added by
`addsource_geometry`, a fypp macro each physics template defines (see
`src/hd/mod_hd_templates.fpp`); `src/physics/mod_physics_dummies.fpp` provides
a version that `mpistop`s, so a physics module without one fails loudly.
`setdt` likewise switches from `rnode` spacing to the physical cell sizes
`ps(igrid)%ds`.

Current limits of the spherical support:

- Only `phys = 'hd'` and `phys = 'mhd'` implement `addsource_geometry`; ffhd
  and srhd will stop at startup if built with `geometry = 'spherical'`. MHD's
  curvature terms (`src/mhd/mod_mhd_templates.fpp`) were ported from upstream
  MPI-AMRVAC's `mhd_add_source_geom` (cell-centred GLM-MHD branch; this fork
  has no `stagger_grid`/constrained-transport support), with the isotropic
  (`ptot`, `psi`) terms rewritten to use the discrete `dAdV` well-balancing
  factor instead of the continuous `2/r`, `cot(theta)/r` prefactors upstream
  uses directly, for consistency with this fork's HD implementation.
- No polar-axis handling. `set_pole`/`poleB` and the pi-periodic root-neighbor
  lookup in `src/amr/mod_amr_neighbors.fpp` survive from upstream, but AGILE's
  rewritten `getbc` in `src/mod_ghostcells_update.fpp` never performs the pole
  copies (`pole_buf` is allocated and otherwise unused). Keep the domain away
  from `theta = 0` and `theta = pi`.
- AMR across curvilinear levels is untested here; prolongation in
  `src/amr/mod_refine.fpp` has the `slab_uniform` branch that uses `dvolume`,
  but `fix_conserve` is commented out in `src/mod_advance.fpp` for all
  geometries.

Validated by `tests/hd/spherical_uniform_flow` (a uniform Cartesian flow
written in spherical components, which converges at second order),
`tests/hd/spherical_blast`, and `tests/mhd/spherical_uniform_flow` (the same
uniform-flow idea extended with a uniform Cartesian magnetic field, to
exercise the induction equation's and GLM psi's curvature terms too).

## Writing a new simulation case (`mod_usr.fpp`)

Every case directory has a `mod_usr.fpp` defining `module mod_usr` with a
`usr_init()` subroutine that:
1. optionally calls `set_coordinate_system(...)` — the coordinate system is
   otherwise taken from `geometry` in `agile.par`, see "Coordinate systems",
2. registers callbacks by assigning the procedure pointers declared in
   `src/mod_usr_methods.fpp` (e.g. `usr_init_one_grid`, `usr_special_bc`,
   `usr_source`, `usr_process_grid`, `usr_refine_threshold`, ...),
3. calls `phys_activate()` last.

The physics module itself (hd/mhd/ffhd/srhd) is selected via `phys` in
`agile.par`'s `&methodlist`, which the config system turns into a `PHYS`
fypp define consumed by `src/physics/mod_physics.fpp` and the per-physics
`mod_*_templates.fpp` files.

## Code architecture

- `src/agile.fpp` — program entry point (`program agile`): MPI init, OpenACC
  device selection, then `main()` which reads parameters, calls
  `usr_init()`, initializes the AMR tree, and runs `timeintegration()` (the
  main time-stepping loop: `setdt` → optional user process hooks → I/O →
  `advance` → AMR regrid → loop).
- `src/mod_global_parameters.fpp` — the large module of shared global state
  (grid geometry, ghost cells, timers, I/O settings) used throughout.
- `src/mod_advance.fpp`, `src/mod_finite_volume.fpp` — the finite-volume
  update step; these contain the performance-critical OpenACC-annotated
  loops (`!$acc parallel loop`, `!$acc routine seq`, etc.).
- `src/amr/` — the block-based octree AMR machinery: forest/tree bookkeeping
  (`mod_forest.fpp`), refinement/coarsening, load balancing, space-filling
  curve ordering, flux correction at refinement boundaries
  (`mod_amr_fct.fpp`, `mod_fix_conserve.fpp`).
- `src/physics/` — physics-agnostic scaffolding: `mod_physics.fpp` defines
  the procedure-pointer interface (`phys_get_flux`, `phys_to_conserved`,
  etc.) that each physics module fills in via `phys_activate()`; also
  radiative cooling, thermal emission, B0-splitting helpers shared across
  physics modules.
- `src/hd/`, `src/mhd/`, `src/ffhd/`, `src/srhd/` — one `mod_*_templates.fpp`
  each, guarded by a top-level `#:if PHYS == '...'` fypp block. These are
  **not** standalone Fortran modules: they define fypp macros (physics-specific
  flux functions, primitive/conservative conversions) that get textually
  `#:include`d — via `src/physics/mod_physics_templates.fpp` — into the real
  module files that use them (`mod_physics.fpp`, `mod_physics_vars.fpp`,
  `mod_finite_volume.fpp`, `mod_dt.fpp`, `mod_source.fpp`,
  `mod_radiative_cooling.fpp`). Only one `PHYS` value's branch is active per
  build, so e.g. `mod_physics.fpp` compiles to a different `mod_physics`
  module depending on which physics module was selected.
- `src/io/` — parameter/config file reading (`mod_config.fpp`), snapshot and
  log I/O, `.dat`/VTU conversion (`mod_convert*.fpp`), slices, collapsed
  views.
- `src/particle/` — optional Lagrangian particle module (advection, Lorentz,
  guiding-center, sampling).
- `src/limiter/`, `src/modules/` — flux limiters and small shared utilities
  (interpolation, ODE integration, RNG, timing, LU solve).
- `src/m_octree_mg_*d.fpp` — vendored multigrid solver (from the
  `octree-mg` submodule, `external_libs/octree-mg`) used for e.g. Poisson
  solves in source terms.
- `make/config_reader.py` reads `agile.par` and, together with
  `make/config_schema.toml`, produces the fypp defines / build hash inputs
  in `config.mk`. `make/fypp-deps.py` and `fortdepend` generate `.fpp`/`.f90`
  dependency rules for incremental rebuilds.

## API documentation (FORD)

Source-level API docs are generated with [FORD](https://forddocs.readthedocs.io/)
from the `!>`/`!<` doc comments (see `doc/documentation.md`), replacing the
old Doxygen setup:

```bash
uv sync --extra docs
source .venv/bin/activate
doc/ford/build_docs.sh   # preprocesses src/*.fpp with fypp, then runs ford
firefox doc/ford/html/index.html
```

`doc/ford/docgen.par` picks one representative build (`PHYS='mhd'`, most
optional features on) to preprocess with, since the `#:if PHYS == '...'`
splicing described above means only one physics module's variant of the
shared files can be documented per pass. `doc/ford/src/` and
`doc/ford/html/` are generated/gitignored, not committed. The narrative
`doc/*.md` pages are not yet wired into the FORD site (`page_dir`) — they
still use Doxygen-specific `{#label}`/`@ref`/`[TOC]` syntax that needs
converting first.

## Code style (from `doc/code_style_guide.md`)

- Fortran 90/95/2003-era style; no implicit typing, no common blocks.
- snake_case names, self-explanatory beyond trivial loop indices.
- 2-space indent for `program`/`module`/`subroutine`/`function`/`associate`;
  3-space for `type`/`interface`/`do`/`if`/`select case`/`where`/`forall`.
- Use `dp` kind and `_dp` suffixes for double precision, not
  `double precision`/bare literals.
- Array loops ordered so the leftmost index varies fastest (Fortran
  column-major order).
- All procedure arguments must declare `intent`.
