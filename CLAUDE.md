# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

AGILE is a GPU-enabled fork of MPI-AMRVAC, a Fortran finite-volume code for
solving hyperbolic PDEs (scalar advection RHO, HD, FFHD, SRHD, MHD, SRMHD) with
adaptive mesh refinement.
Master supports 3D Cartesian, 3D spherical `(r, theta, phi)` and 3D
cylindrical `(r, z, phi)` grids, each of the latter two also in a
logarithmically stretched radial variant (`logSpherical`, `logCylindrical`);
see "Coordinate systems" below for the limits of the curvilinear support.
Domain decomposition uses MPI; GPU offload supports two selectable backends,
OpenACC and OpenMP target offload. Every offload site in the hot loops is
written once as a call to a macro from `src/mod_gpu_directives.fpp` (e.g.
`${GPU_PARALLEL_LOOP(...)}$`), which fypp expands to the matching
`!$acc ...` or `!$omp target ...` directive depending on which backend is
selected at build time (see `OPENACC=1`/`OPENMP=1` below) — don't write raw
`!$acc`/`!$omp` directives directly in hot-loop code, add a macro call
instead.

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
are required for every build. Requires Python >=3.13. `tools/vscode/setup.sh` sets
up VS Code for the `.fpp` sources — fypp syntax highlighting, `fortls` for
navigation (the `ide` extra), and a linter that runs fypp before gfortran; see
`tools/vscode/README.md`.

## Building

Each test/run case is its own directory with a `Makefile` and an `agile.par`
parameter file. Build from inside such a directory:

```bash
cd tests/hd/KH3D
make arch=gnu              # compiles and links ./agile executable
make arch=gnu OPENACC=1    # enable GPU build with the OpenACC backend
make arch=gnu OPENMP=1     # enable GPU build with the OpenMP target-offload backend
make arch=gnu DEBUG=1      # debug flags (-g -O0 -fcheck=all, etc.)
make clean                 # remove this case's build products
```

- `arch` selects a file from `arch/*.mk` (`gnu`, `ifx`, `nvidia`, `cray`,
  `llvm`) which sets the compiler/linker and flags. Default is `gnu`
  (`mpif90`/gfortran).
- `OPENACC` and `OPENMP` are mutually exclusive GPU offload backend switches
  (`make` errors out if both are given); neither flag gives a plain CPU
  build, where every `GPU_*` directive macro expands to nothing.
- `NOGPUDIRECT=1` (only meaningful with `OPENACC` or `OPENMP`) stages MPI
  buffers through the host instead of handing device pointers to MPI. **Use
  it when all ranks share a single GPU.** Without it, a workstation with one
  GTX 1660 Ti (nvhpc 26.9, HPC-X) runs `tests/hd/KH3D` correctly on one
  rank but loses mass deterministically on two or more (`mean(rho)` 5.42 at
  2 ranks, 5.37 at 4, against 5.5), whichever MPI build or UCX transport is
  used; with `NOGPUDIRECT=1` it passes at 1, 2 and 4 ranks. The cause is not
  understood, so treat this as a workaround, not a known limit. Multi-GPU
  machines (e.g. the GPU CI) keep GPU-direct.
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

Tests live under `tests/{rho,hd,ffhd,srhd,mhd,srmhd}/<case>/` and each case has a
`test.make` including `tests/test_rules.make`, which builds the case fresh
and runs it (via `mpirun`/`srun`), then numerically diffs the log output
against `tests/<physics>/<case>/correct_output/` using
`tools/fortran/compare_logs` (tolerances `1.0e-5` relative / `1.0e-8`
absolute).

```bash
cd tests
make rho           # run the basic rho/* test cases
make hd            # run the basic hd/* test cases
make ffhd          # run the basic ffhd/* test cases
make srhd          # run the basic srhd/* test cases
make mhd           # run the basic mhd/* test cases
make srmhd         # run the basic srmhd/* test cases
make SUITE=all hd  # basic and extended hd/* cases
make SUITE=all list # print the directories a SUITE selects
make -j4 hd         # run in parallel
make -s hd          # suppress "Entering/leaving directory" noise
make clean          # clean all test build artifacts
```

The test directories are split into two suites in `tests/Makefile`: each
physics module has a `*_BASIC` and a `*_EXTENDED` list. `SUITE` (a make
variable or environment variable, like `arch`) selects `basic` (the default),
`extended` or `all`. The split is per directory, because a directory is one
build and builds dominate the cost. The basic suite keeps one Cartesian case
per physics plus the curvilinear cases whose checks nothing else provides
(`hd/spherical` for refluxing, `hd/log_cylindrical_pole` for `log_r0`, ...).
Naming a single directory (`make hd/spherical`) runs it regardless of `SUITE`.

In CI (`.github/workflows/tests.yml` and `gpu.yml`) pushes and PRs run the
basic suite. The full suite (`all`) runs weekly on master (Sunday 02:00
UTC), on a PR carrying the `full-tests` label (adding the label triggers the
run), and from a manual run, whose `suite` input chooses between
`basic`, `extended` and `all` (`gh workflow run tests.yml -f suite=all --ref <branch>`).

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
`tests/Makefile` under the appropriate `*_BASIC` or `*_EXTENDED` variable (extended unless
it covers something no basic case does).

`tests/hd/snapshot_zfp` (ZFP-compressed snapshots through the parallel and the
default writer, byte for byte) is in no suite: ZFP is a compile-time switch
that needs a ZFP build, which CI lacks. Run it by hand with
`make -f test.make arch=... ZFP_DIR=/path/to/zfp` (see
`doc/par.md#par_snapshots` for the ZFP build); without `ZFP_DIR` it prints
SKIPPED.

## Coordinate systems

The coordinate system is a **compile-time** choice (the finite-volume kernels
are fypp-specialized on it). Set `geometry` in `&meshlist` of `agile.par` to
`'Cartesian'` (default), `'spherical'`, `'cylindrical'`, `'logSpherical'` or
`'logCylindrical'`. It becomes the `GEOM` fypp define, and `main()` calls
`set_coordinate_system_from_config()` before `usr_init`, so cases need not call
`set_coordinate_system` (if they do, the name must agree with `GEOM` or it
`mpistop`s).

Things to know when writing a curvilinear case:

- **Angular coordinates in `agile.par` are in units of 2*pi**
  (`xprobmin2 = 0.125d0` is `theta = pi/4`). Cylindrical is ordered
  `(r, z, phi)`, so only dimension 3 gets that conversion there.
- `ndim = 3` always, so curvilinear means full 3D with `ndir = 3`.
- **`ps(igrid)%x` is the cell's volume barycentre, not its face midpoint**, so
  `x +/- dx/2` is not a face; rebuild faces from `rnode` (block corner + index)
  as `mod_finite_volume` and `calc_x` do. `bgeo%ds` (the CFL length used by
  `setdt`) deliberately stays on the midpoint. `divvector`/`curlvector` in
  `src/mod_geometry.fpp` still assume midpoints; they are not on any exercised
  path.

### Geometry storage

Metrics live in `bgeo` (and `bgeoc` for the coarse representatives), a `geo_t`
in `src/mod_physicaldata.fpp` holding `x`, `ds`, `dvolume`, `surfaceC`, `dx`
for all `1:max_blocks` blocks with the grid index last. `ps(igrid)%x` etc. are
bounds-remapped pointer views into it (`point_at_geometry`). The metrics are
computed analytically on the device by `fill_geometry_device`
(`src/amr/mod_amr_solution_node.fpp`), called from `alloc_node`; all of `geo_t`
is device-resident and sized for `max_blocks` up front, which can exhaust GPU
memory for large `max_blocks` — don't add members nothing reads (`dsC` and
`surface` were removed for this; `s%dsC`/`s%surface` are never associated and
their former readers `mpistop`).

General grid stretching (`stretch_dim` etc.) is removed; only the log radius
below is supported.

**Host sync.** `alloc_node` copies nothing back. `sync_positions_host([igrid])`
and `sync_geometry_host()` in `src/mod_geometry.fpp` fetch the geometry at the
sites that read it (output in `saveamrfile`, initial conditions, user process
hooks, `prolongprimitive`, ...). **Adding a host reader of `ps(igrid)%x` or of
any metric means adding the matching sync call** — a missed one silently reads
stale data from whichever block previously held that `igrid` slot.

### Logarithmic radius

`logSpherical`/`logCylindrical` use a mesh uniform in the logical coordinate
`s = ln(1 + r/r0)` (`r0 = log_r0` in `&meshlist`, runtime; `log_r0 = 0` means
plain `s = ln(r)`). They emit the same `GEOM` plus a `LOG_RADIUS` flag, so
every `GEOM == 'spherical'` branch applies unchanged and the AMR machinery works
in `s` untouched.

- **`xprobmin1`/`xprobmax1` are physical radii in the par file but logical
  after `read_par_files`.** A `mod_usr.fpp` that needs a physical radius must
  map back with `r_of_s` (inverse `s_of_r`) from `mod_geometry`. Older
  `log_*` tests write `exp(xprobmin1)` inline, which is only right for
  `log_r0 = 0` — don't copy that.
- `slicecoord` is *not* converted; it is physical throughout the slice path.
- `log_r0 > 0` lets the grid reach `r = 0` (cylindrical axis). The map is
  odd (`r = sign(s) r0 (exp|s| - 1)`) so ghost cells across the axis exactly
  mirror the interior; `logSpherical` still has no usable origin.
- Reconstruction assumes uniform cells, so keep the stretch modest
  (tests use `q ~ 1.04-1.09`).

### Physics in curvilinear coordinates

- All six physics modules implement `addsource_geometry` per `GEOM` in their
  `mod_*_templates.fpp`. The geometric prefactors come from the discrete
  `dAdV = (A_upper - A_lower)/dV`, which are exactly the volume averages of the
  continuous `2/r`, `cot(theta)/r`, `1/r` — this makes the scheme
  well-balanced and independent of the stored position `x` (unused there).
  rho and ffhd have empty source terms (scalars). ffhd's optional `PDIVB` term
  is **not** correct in curvilinear coordinates.
- SRHD/SRMHD primitives hold `u^i = lfac*v^i`; the SRMHD stress uses the
  lab-frame `B^2` in `W = rho h lfac^2 + B^2` (correct, not a slip).
- **Refluxing** (`src/mod_fix_conserve.fpp`) stores `qdt*f/dx` in Cartesian and
  the extensive `qdt*f*A` in curvilinear builds (divided by the coarse
  `dvolume`), selected by fypp. It breaks exact well-balancing at coarse-fine
  faces (second order). Messages are aggregated per destination rank. Level
  jumps across a polar axis are not refluxed. **Check reflux changes against
  `tests/hd/spherical/blast_amr_fixgrid.par`** (closed domain, frozen mesh),
  not `blast_amr.par`, whose log is insensitive to refluxing.

### The polar axis

A spherical domain may span `theta` 0 to `pi`, a cylindrical one may start at
`r = 0`. Ghost cells across the axis are copied from the block half a turn away
in `phi` (in `getbc`, `src/mod_ghostcells_update.fpp`, via `neighbor_pole`).
`check_pole_setup` enforces:

```fortran
 &boundlist
   typeboundary_min2 = 5*'pole'       ! axis face (min1 for cylindrical)
   typeboundary_min3 = 5*'periodic'   ! phi periodic over a full turn
   typeboundary_max3 = 5*'periodic'
 /
```

plus an even number of level-1 blocks in `phi`. Sign rule: the component along
the mirrored direction and the `phi` component flip, everything else is
symmetric (for cylindrical that includes `m_r`/`B_r` — upstream gets this
wrong). ffhd needs no flips. The timestep collapses near the axis
(`ds(3) -> 0`); keep resolution modest.

**Volume averages cannot validate the pole copy** (axis cells have vanishing
volume), so pole tests assert on the ghost cells directly with
`check_pole_ghosts` at `it = 0`. **A new pole test needs that check, not just a
reference log.** Across level jumps, point-vs-analytic only works where the
state is smooth on the coarse-cell scale.

### The frozen field (`phys = 'ffhd'`)

ffhd's field direction `b-hat` (slots `b1..b3`, the `nwextra` variables) is set
analytically by the by-name device hook `usr_set_nwextra(x, wextra)` (must
return a unit vector), called by `fill_nwextra_device` from `alloc_node` for
every cell including all ghosts (`FILL_NWEXTRA_ANALYTIC` flag). It is not
exchanged by `getbc` (`nwgc = nwflux`), and `initonegrid_usr`/`usr_special_bc`
must **not** set it. Across the axis the mirrored `x` makes the hook produce
the correct signs on its own.

### Curvilinear tests

Directories are per (physics, geometry) build: `tests/<phys>/{spherical,
cylindrical}`, `.../log_{spherical,cylindrical}`, `.../{spherical,
cylindrical}_pole`, plus `tests/hd/log_{spherical,cylindrical}_pole`. Runs in
one directory share a build and select the setup via `setup` in `&usr_list`
(`uflow.par` uniform flow — the well-balancing test; `blast.par`; AMR
variants). `agile.par` is the build reference and must declare the union of
compile-time options every par file in the directory needs. **Adding a par file
is free; adding a directory costs a CI rebuild.**

Log-geometry tests register `check_log_grid`, which checks the mesh against its
closed form cell by cell; `tests/hd/log_cylindrical_pole` (`log_r0 = 0.2`)
adds `check_pole_mesh`. **A new log-geometry test needs that check, not just a
reference log** — a subtly wrong mesh barely moves the volume averages.

## Data-driven boundaries (`bc_data`, `bc_icarus`)

`src/mod_bc_data.fpp` reads one legacy VTK file (axes: the face's two
tangential coordinates, then time) into a single table, `bc_values`, with
one column per (variable, face), numbered by `bc_data_ix`. It is copied to
the device once, and `bc_data_set` is one per-cell `GPU_ROUTINE_VECTOR`
routine called from `bc_phys` on host and device alike, reading the table
through `bc_locate`/`bc_interp` (`GPU_ROUTINE_SEQ`). Keep it one code path:
there used to be a device-only copy of the data and interpolation. The table
is plain module arrays on purpose: an `LT3_t` (derived type with an
allocatable component, put on the device like `bgeo` and passed whole into a
`routine seq`) worked with nvhpc 26.9 but faulted with
`CUDA_ERROR_ILLEGAL_ADDRESS` under the CI's 25.11, while compute-sanitizer
was clean locally. `bc_icarus` is `bc_data`
followed by `specialbound_usr`. Options are in `doc/par.md`
(`#par_bc_data`); tests are `tests/rho/vtk_boundary` (basic) and
`tests/hd/vtk_boundary` (extended, mixed faces and primitive data), both
asserting the ghost cells against the table.

## Time integration

Set `time_stepper` and `time_integrator` in `&methodlist`. `advect` in
`src/mod_advance.fpp` implements these; the default `time_integrator` of each
`time_stepper` is in bold, and the overall default is `twostep`, so
`Predictor_Corrector`:

| `time_stepper` | `time_integrator` | order | stages | time levels |
|---|---|---|---|---|
| `onestep` | **`Forward_Euler`** | 1 | 1 | 2 |
| `twostep` | **`Predictor_Corrector`** | 2 | 2 | 2 |
| `twostep` | `ssprk2` | 2 | 2 | 2 |
| `twostep` | `RK2_alfa` (`rk2_alfa` in (0,1]) | 2 | 2 | 2 |
| `threestep` | **`ssprk3`** (`ssprk_order` 3 or 2) | 3 / 2 | 3 | 3 |
| `threestep` | `RK3_BT` (`rk3_switch` 1-4: Ralston, Wray, Heun, Nyström) | 3 | 3 | 4 |
| `fourstep` | **`ssprk4`** (`ssprk_order` 3: Spiteri-Ruuth, 2: SSP(4,2)) | 3 / 2 | 4 | 3 |
| `fourstep` | `rk4` (classical) | 4 | 4 | 4 |
| `fivestep` | **`ssprk5`** (Spiteri-Ruuth (5,4)) | 4 | 5 | 4 |

That is every explicit scheme upstream MPI-AMRVAC has. `read_par_files` also
parses the IMEX family and fills its coefficients, but stops at startup for
it: there is nothing implicit to integrate here, since every
`phys_implicit_update` is the dummy in `src/physics/mod_physics.fpp`.

The SSP schemes are stable to different Courant numbers relative to forward
Euler (ssprk4: 2 at order 3, 3 at order 2), but AGILE does not rescale
`courantpar` per scheme, so a run has to raise it itself to profit.

**Stages and time levels are different numbers, and only the second costs
memory.** `nstep` counts sub-steps (flux evaluations). `nlevel` counts the
solution copies the scheme holds at once, its registers, and `initialize_vars`
allocates `bg(1:nlevel)`. Each level is `nw` times the block size times
`max_blocks` doubles, all device-resident, so it often decides how large
`max_blocks` can be on a GPU. The two counts disagree in both directions:
ssprk4 is 4 stages in 3 levels, while RK3_BT is 3 stages needing 4 levels.
`alloc_node` used to size the levels from the integrator in one place while
`bg` was sized from `nstep` in another, which left RK3_BT and the IMEX schemes
pointing past the end of `bg`.

Both numbers now come from one table, `time_levels` in
`src/mod_global_parameters.fpp`. **A new scheme adds its row there**, a `case`
branch in `advect`, and its name to the accepted list in `read_par_files`.
Everything else follows from `nlevel`: `bg` allocation, the `ps1..ps4` views
that `alloc_node` points at `bg(2..nlevel)` through `ps_level`
(`src/mod_physicaldata.fpp`), and the guard in `alloc_state` that stops on a
view beyond `size(bg)`.

Within `advect` a scheme is written as upstream writes it, with `ps`, `ps1`,
`ps2`, `ps3` becoming time levels 1, 2, 3, 4:

- `advect1(dtfactor, ..., qtC, ka, qt, kb)` adds `dtfactor*dt` times the
  update evaluated on level `ka` to level `kb`. The two must differ, because
  the flux stencil reads neighbouring cells.
- `bg_assign(bgd, a1, bg1[, a2, bg2[, a3, bg3]])` (`bgd = sum`) and
  `bg_accum(bgd, s, a1, bg1[, ...])` (`bgd = s*bgd + sum`) are the pointwise
  combinations between stages. They are generic names over device kernels
  that fypp generates for one to three sources, so each reads only the
  arrays it uses. There are two shapes so that a destination is never also
  passed as a source argument, and is never read where its coefficient would
  be zero.

The schemes are tested on the AMR case of `tests/hd/KH3D-SMR`, not on a
uniform grid, because the static refinement is what puts refluxing at every
stage under test. `tests/hd/KH3D-SMR` itself (basic suite) runs ssprk3,
ssprk4 and the default, Predictor_Corrector: its `pc.par` sets no
`time_stepper`, so the same test also checks the default.
`tests/hd/KH3D-SMR-integrators` (extended suite) is the same case with one
par file per scheme and coefficient variant: all twelve, including
`ssprk_order = 2` for ssprk3 and ssprk4, and two `rk3_switch` tableaus
(Wray, Nyström). It is a directory of its own, rather than par files added
to `KH3D-SMR`, because the suite split is per directory.

Refluxing runs after every stage except the predictor of
`Predictor_Corrector` (as upstream), and only when the mesh has more than one
level (`fix_conserve_at_step`).

## MHD units

- `unit_magneticfield = sqrt(4 pi unit_pressure)` is set in `phys_units`
  (it was left at 1 before, which scaled the Braginskii kappa_perp through
  `htc_Cchi` by the field unit). `H_ion_fr`, `He_ion_fr`, `He_ion_fr2` in
  `&mhd_list` enter b = 1 + H_ion_fr + He_abundance (...), the particles per
  hydrogen nucleus in p = b n kB T; the defaults give the old
  b = 2 + 3 He_abundance, and `He_abundance = 0, H_ion_fr = 0` gives
  Legolas' p = n k T.

## Scalar advection (`phys = 'rho'`)

`src/rho/mod_rho_templates.fpp` solves d(rho)/dt + div(rho v) = 0 for a
constant Cartesian velocity `rho_v(1:3)` (`&rho_list`; `rho_v = 0` is
refused, since the CFL step would be infinite). It is a port of upstream
MPI-AMRVAC's `src/rho`, minus `rho_particles`: one variable
(`nw_phys = nw_flux = 1`, no auxiliaries), primitive equal to conservative.
`tvdlf` and `hll` both reduce to exact upwinding, because `get_cmax` and
`estimate_speeds_minmax` use the exact speed |v_n|; `hllc` is refused in
`read_par_files`. Curvilinear meshes and the polar axis are supported, see
the `addsource_geometry` bullet under "Coordinate systems".

Tests: `tests/rho/ball_3D` (basic) advects upstream's discontinuous ball on a
periodic cube with Lohner AMR (`test.par`), the same with a frozen mesh so the
ball crosses level jumps (`test-fixgrid.par`, the refluxing test) and on a
uniform grid with HLL/koren/twostep (`hll.par`). Each logs at t = 0 too
(`itsave(1,1) = 0`; `dtsave_log` alone first fires one interval in), so a
mass drift above the `compare_logs` tolerance fails the test.
`tests/rho/cylindrical_pole` (basic), `rho/spherical` and
`rho/spherical_pole` (extended) each run `uflow.par` (uniform rho) and
`blob.par` (a Gaussian carried across the domain, through the axis in the pole
cases, where `check_pole_ghosts` asserts the ghost cells at it = 0). Note that
uflow's O(h^2) free-stream deviation (<= 5e-5) is near or below the
`compare_logs` tolerance of about 1e-5 for values near one, so the blob runs
and the pole checks carry the sensitivity there.

`tests/rho/convergence/` is a convergence study, not part of any suite:
`run_convergence.sh` builds once and sweeps limiter x time stepper x flux x
resolution for a smooth wave whose initial condition and reference are exact
cell averages, with all output in `tests/rho/output/convergence/`;
`convergence.py` prints the tables and orders and draws `convergence.png`
(matplotlib: `uv sync --extra analysis`). Measured L1 orders at N = 200-400
with ssprk3: koren 2.4 (third order away from the extrema, clipped at them),
mcbeta and vanleer 2.0, minmod 1.9; L-infinity orders are 1.3-1.6, set by the
extrema. cada3 (default `cada3_radius`) is third order in L1, L2 and
L-infinity alike once resolved (3.00 at N = 800-1600, pre-asymptotic orders
of 4-7 between 200 and 400), because it does not clip smooth extrema.
`onestep` caps every limiter at order 1, and `twostep` caps koren and cada3
near 2, so only ssprk3 lets their accuracy through. hll equals tvdlf to
round-off except under onestep, which amplifies round-off. A single plane
wave, `setup = 'diag'` included, is blind to the O(h^2) face-quadrature error
of the 3D update: its errors match the quasi-1D wave's to 0.1%.

## Writing a new simulation case (`mod_usr.fpp`)

Every case directory has a `mod_usr.fpp` defining `module mod_usr` with a
`usr_init()` subroutine that:
1. registers callbacks by assigning the procedure pointers declared in
   `src/mod_usr_methods.fpp` (e.g. `usr_init_one_grid`, `usr_special_bc`,
   `usr_process_grid`, `usr_refine_threshold`, ...),
2. calls `phys_activate()` last.

`use mod_agile` brings in the modules a case usually needs (global
parameters, user methods, geometry, `mpistop`). `mod_amrvac` is kept only as
an alias for it, so that cases ported from upstream compile unchanged; new
cases use `mod_agile`.

Hooks that run on the GPU are not procedure pointers but routines the core
calls by name from `mod_usr`, compiled in by a switch in `agile.par`:
`gravity_field` (`*_gravity`), `addsource_usr` (`*_source_usr`; the
`usr_source` pointer is a leftover and never called), `usr_refine_grid`
(`refine_usr`), `usr_internal_bc` (`internalboundary` in `&boundlist`) and
`usr_set_nwextra` (ffhd). `usr_internal_bc(level, qt, x, w)` is a per-cell
internal boundary: `getintbc` applies it to every interior cell of the time
level `getbc` is exchanging, before the exchange, writing back only that
call's variable range (see `doc/par.md`, `tests/hd/internal_bc`). They are marked with
`${GPU_ROUTINE_SEQ()}$`, so `mod_usr.fpp` starts with

```fortran
#:include "mod_gpu_directives.fpp"
```

by bare file name: `make/30-fypp.mk` passes `-I $(agile)/src` to fypp (and
`make/fypp-deps.py` searches the same directories for its dependency
rules), so the include works wherever the case directory is. The tests used
to write `"../../../src/mod_gpu_directives.fpp"`, which only worked at that
depth inside the repository.

The coordinate system comes from `geometry` in `agile.par` (see "Coordinate
systems" above); `usr_init()` does not need to call `set_coordinate_system`
itself. Older cases may still call it explicitly — the call is idempotent as
long as the name agrees with `geometry` — but new cases shouldn't.

The physics module itself (rho/hd/ffhd/srhd/mhd/srmhd) is selected via `phys` in
`agile.par`'s `&methodlist`, which the config system turns into a `PHYS`
fypp define consumed by `src/physics/mod_physics.fpp` and the per-physics
`mod_*_templates.fpp` files.

## Code architecture

- `src/agile.fpp` — program entry point (`program agile`): MPI init, GPU
  device selection (`set_openacc_device`/`set_openmp_device`, whichever
  backend is active), then `main()` which reads parameters, calls
  `usr_init()`, initializes the AMR tree, and runs `timeintegration()` (the
  main time-stepping loop: `setdt` → optional user process hooks → I/O →
  `advance` → AMR regrid → loop).
- `src/mod_global_parameters.fpp` — the large module of shared global state
  (grid geometry, ghost cells, timers, I/O settings) used throughout.
- `src/mod_gpu_directives.fpp` — fypp macro definitions (`GPU_PARALLEL_LOOP`,
  `GPU_ROUTINE_SEQ`, `GPU_ENTER_DATA_COPYIN`, `GPU_HOST_DATA_USE_DEVICE`,
  etc.) that expand to OpenACC or OpenMP-target directives depending on the
  active backend; `src/mod_gpu_utils.fpp` (module `gpu_utils`) provides the
  `copy_or_update`/`copy_or_update_pointer`/`copy_or_update_alloc` helpers
  used when re-allocating AMR grid data on the device.
- `src/mod_advance.fpp`, `src/mod_finite_volume.fpp` — the finite-volume
  update step; these contain the performance-critical GPU-offloaded loops,
  written via the `GPU_*` directive macros (`GPU_PARALLEL_LOOP`,
  `GPU_ROUTINE_SEQ`, etc.) rather than raw `!$acc`/`!$omp` directives.
  `FV_KERNEL` is instantiated per flux scheme *and* per layout, `per_cell`
  (the default) and `shared`, picked at run time by `fv_face_flux` in
  `&methodlist` (see `doc/par.md`). The two differ only in where a face
  flux is computed; both build every face through the `FACE_FLUX` macro and
  update every cell through `FV_CELL_UPDATE`, so a change to the
  reconstruction, the face position or the update/flux-fixing stores goes
  in those macros once and the layouts stay bit-identical. Check that they
  do: the logs of the two must match exactly, not just within tolerance.
- `src/amr/` — the block-based octree AMR machinery: forest/tree bookkeeping
  (`mod_forest.fpp`), refinement/coarsening, load balancing, space-filling
  curve ordering, flux correction at refinement boundaries
  (`mod_amr_fct.fpp`, `mod_fix_conserve.fpp`).
- `src/physics/` — physics-agnostic scaffolding: `mod_physics.fpp` defines
  the procedure-pointer interface (`phys_get_flux`, `phys_to_conserved`,
  etc.) that each physics module fills in via `phys_activate()`; also
  radiative cooling, thermal emission, B0-splitting helpers shared across
  physics modules.
- `src/rho/`, `src/hd/`, `src/ffhd/`, `src/srhd/`, `src/mhd/`, `src/srmhd/` — one `mod_*_templates.fpp`
  each, guarded by a top-level `#:if PHYS == '...'` fypp block. These are
  **not** standalone Fortran modules: they define fypp macros (physics-specific
  flux functions, primitive/conservative conversions) that get textually
  `#:include`d — via `src/physics/mod_physics_templates.fpp` — into the real
  module files that use them (`mod_physics.fpp`, `mod_physics_vars.fpp`,
  `mod_finite_volume.fpp`, `mod_dt.fpp`, `mod_source.fpp`,
  `mod_radiative_cooling.fpp`). Only one `PHYS` value's branch is active per
  build, so e.g. `mod_physics.fpp` compiles to a different `mod_physics`
  module depending on which physics module was selected.
- `src/io/` — parameter/config file reading (`mod_config.fpp`), restart-file,
  snapshot and log I/O, `.dat`/VTU conversion (`mod_convert*.fpp`), slices,
  collapsed views.
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

## Bug-hunting notes

- **A reduction variable must not be passed by reference into a routine
  called from the reduced loop.** `forcedrefine_grid` combined the per-block
  refine/coarsen flags with `reduction(.or.:refineflag) ...` on the vector
  loop over a block's cells and passed those same flags as `intent(inout)`
  arguments to `usr_refine_grid`. The callee then wrote to the shared
  variable, not the lane's reduction copy, and only the first cell's value
  survived: a `usr_refine_grid` rule was applied to one cell per block, the
  mesh came out lopsided, and nothing complained. CPU builds were unaffected,
  which is what made it look like a physics asymmetry. The call now takes
  per-cell `private` flags that the loop body combines into the reduction
  variables itself; the loop stays vectorised.

- **When CPU is correct and one GPU compiler is wrong and the failure is
  deterministic, it is a miscompile — bisect the kernel's directives
  mechanically, do not reason about which construct "should" be at fault.**
  Toggle one thing at a time in the offending kernel (collapse level;
  array-section assignment vs explicit scalar loop; `vector` placement; loop
  bounds; `-O` level) and watch a bucketed error metric — e.g. split the
  ghost-region error into face / edge / corner, or per variable. This is far
  faster than theorising: the polar-axis GPU bug was chased through several
  wrong hypotheses (the MPI path, a reversed array section) before mechanical
  bisection pinned it to `collapse` in one `!$acc loop` in three toggles.
  Cross-check with `arch=gnu` (gfortran) and `arch=nvidia` *without*
  `OPENACC=1` (nvfortran, host) to confirm it is the OpenACC offload and not
  the front end.

- **In a `collapse`d `!$acc loop` inside an `!$acc routine`, a `call` in the
  loop body needs every one of the body's `private` automatic arrays to have
  a compile-time size.** nvfortran's OpenACC miscompiles a collapsed `vector`
  loop — reached from a gang-level `!$acc parallel loop` through a `!$acc
  routine vector` — when the body both `call`s another `!$acc routine` *and*
  declares a `private` automatic array whose size is a run-time value: it
  gets that array's per-lane stride wrong and scatters the results into the
  wrong cells (the ghost layer comes back index-rotated; results drift rather
  than crash). Both conditions are needed. Remove the `call` (inline the
  callee) *or* give every such `private` array a `parameter` bound and the
  same `collapse(3)` is generated correctly. The reduced reproducer is on
  branch `reproducer/issue-154` (`tests/regression/issue154/repro154.f90`);
  its variant 4 keeps `collapse(3)` and the `call` and only changes
  `wpt(1:nw)` to `wpt(1:nw_phys)`.

  This bit every hd curvilinear `*_pole` case's `specialbound_usr`, whose
  `wpt(1:nw)` private buffer is sized by the run-time `nw` — logs stayed put
  but `check_pole_ghosts` failed by ~0.4 on nvfortran while gfortran and
  nvfortran-without-OpenACC were exact. **The fix now in the tree** is
  `wpt(1:nw_phys)` (the compile-time conserved-variable count from the
  physics templates, `= nw` for these cases), `collapse(3)` kept. The
  mhd/srhd pole cases call a helper too but their `private` arrays
  (`v`, `b`, `x_loc`) are already compile-time sized, and the ffhd pole
  cases have no `call` in the loop at all, so all six also keep
  `collapse(3)`. Dropping the `collapse` (`!$acc loop vector` on the
  innermost loop, outer loops sequential) also stops the corruption and was
  the first fix shipped, but the reproducer shows it still silently reorders
  a short vector loop that has a `call` and a run-time-sized `private` — it
  only looked correct because the pole cases' boundary state does not vary
  along that axis. Collapse is fine at a top-level `!$acc parallel loop`.
  This is very likely the same defect as issue #154
  (`to_primitive`/`to_conservative` called from a `collapse`d loop in
  `specialbound_usr` returning nonsense).

- **Under NVHPC the OpenMP gang/vector macros expand to `target teams loop`
  and `loop bind(parallel)`, not `teams distribute` and `parallel do`.**
  nvfortran compiles `teams distribute` with a nested `parallel do` in
  generic mode: the inner region becomes an outlined device function (the
  `NVFORTRAN-W-0155` "Passing external procedure as argument" warning), and
  every kernel takes ~255 registers per thread and ~9 KB of shared memory,
  however simple. On a GTX 1660 Ti, KH3D at 192^3 took 58.0 s against
  OpenACC's 45.3 s, with the ghost-cell kernels 4x slower; with `loop` it
  took 46.4 s. `bind(parallel)` is not optional: without it nvfortran ran its
  own dependence analysis on the finite-volume kernel, gave up on its
  private temporaries, and launched it with one thread per team, 25x slower,
  without any error. A returning `W-0155` warning means the generic-mode
  pattern is back. Other compilers keep `distribute`/`parallel do`, and
  should until benchmarked: LLVM flang (and so AMD's `amdflang`) lowers a
  nested `loop bind(parallel)` to a worksharing `do` with no enclosing
  `parallel` (`GenericLoopConversion.cpp`), which leaves one thread per
  team, and Cray CCE documents the same `bind(parallel)` -> `do` mapping
  with no `loop`-specific optimisation. MFC's `omp_macros.fpp` likewise uses
  `loop` for NVHPC only.

- **A declared derived type with an allocatable component, read inside an
  `acc routine`, crashed the CI's nvfortran 25.11 (OpenACC); keep such data
  in plain module arrays.** `bc_data` once kept its table in an `LT3_t`
  (`GPU_DECLARE_CREATE`, `rows_cols` attached with enter data) and read it
  in `bc_data_set`, an `acc routine` below the `getbc` kernel. nvfortran
  26.9 ran it correctly and compute-sanitizer was clean, but 25.11 died with
  `CUDA_ERROR_ILLEGAL_ADDRESS` in the first `getbc` in four variants: the
  struct updated whole, updated member by member, only its components
  (`t%rows_cols`, `t%x_min`, ...) passed on to array-level `routine seq`
  lookups, and that last one with the table declared `target`. The same
  code with plain module arrays passed (CI run 45819295535). Why `bgeo`
  works there (declared `target`, components attached, and
  `bgeo%x(:,:,:,:,igrid)` read in the routine `fill_boundary_before_gc`) is
  not established; it is a different module's variable, and its routine is
  in that other module too.

- **Only five slope limiters exist in the kernel.** `muscl_reconstruct_prim`
  in `src/mod_finite_volume.fpp` implements minmod, vanleer, mcbeta, koren and
  cada3; `read_par_files` stops on any other limiter.

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

`build_docs.sh` preprocesses every file once per physics module
(`doc/ford/docgen_{rho,hd,ffhd,srhd,mhd,srmhd}.par`), since the `#:if PHYS == '...'`
splicing described above means one pass only shows one physics module's
variant of the shared files; `merge_physics_variants.py` keeps files that
come out identical once and suffixes the rest (`mod_physics_hd`, ...).
`doc/ford/src/` and `doc/ford/html/` are generated/gitignored, not committed.
The narrative `doc/*.md` pages are part of the site (`page_dir: ..` in
`doc/ford/project.md`), as FORD pages with a `title:` header; the time
integrators are described in `doc/time_integration.md`. Page navigation on the site is
a "Contents" drawer rendered from `doc/contents.md` itself
(`doc/ford/templates/info_page.html`), so a new page must be listed there to
be reachable.

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
