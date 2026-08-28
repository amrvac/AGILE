# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

AGILE is a GPU-enabled fork of MPI-AMRVAC, a Fortran finite-volume code for
solving hyperbolic PDEs (HD, MHD, FFHD, SRHD) with adaptive mesh refinement.
Master supports 3D Cartesian, 3D spherical `(r, theta, phi)` and 3D
cylindrical `(r, z, phi)` grids; see "Coordinate systems" below for the
limits of the curvilinear support. Domain
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
make srhd          # run all srhd/* test cases
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
   geometry = 'spherical'   ! or 'cylindrical', or 'Cartesian' (the default)
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

Two things to know when writing a spherical or cylindrical case:

- **Angular coordinates in `agile.par` are given in units of 2*pi** (an
  MPI-AMRVAC convention, applied in `src/io/mod_input_output.fpp`). So for
  spherical, `xprobmin2 = 0.125d0`, `xprobmax2 = 0.375d0` means `theta` from
  `pi/4` to `3*pi/4`. For cylindrical, `set_coordinate_system` assigns
  `r_ = 1`, `z_ = 2`, `phi_ = 3` for the 3D case (`"cylindrical"` /
  `"cylindrical_3D"`), so only `xprobmin3`/`xprobmax3` get the 2*pi
  conversion; `xprobmin2`/`xprobmax2` (`z`) is a plain length like `r`.
- `ndim` is a compile-time `parameter` of 3, so spherical or cylindrical
  means the full 3D `(r, theta, phi)` or `(r, z, phi)` system with
  `ndir = 3`.

How it works: in a curvilinear build the finite-volume update in
`src/mod_finite_volume.fpp` divides face fluxes by the cell volume weighted by
the face areas (`bgeo%surfaceC`, `bgeo%dvolume`, filled per block by
`fillgeo`/`get_surface_area` in `src/amr/mod_amr_solution_node.fpp`), instead
of by the uniform `rnode` cell spacing. The leftover curvature terms of the
momentum equations are added by `addsource_geometry`, a fypp macro each physics
template defines (see `src/hd/mod_hd_templates.fpp`);
`src/physics/mod_physics_dummies.fpp` provides a version that `mpistop`s, so a
physics module without one fails loudly. `setdt` likewise switches from `rnode`
spacing to the physical cell sizes `bgeo%ds`.

A block's geometry lives in `bgeo` (and `bgeoc` for the one-level-coarser
representatives), a `geo_t` declared in `src/mod_physicaldata.fpp` that holds
`x`, `ds`, `dvolume`, `surfaceC` and `dx` for *all* blocks with the grid index
last, exactly as `block_grid_t` holds the solution in
`bg%w`. `initialize_vars` allocates them once for `1:max_blocks`, so GPU
kernels index them directly (`bgeo%x(ix1,ix2,ix3,1:ndim,igrid)`) rather than
chasing a per-block pointer, and `alloc_node` only has to push the one block's
slice it just filled. The matching `state` components (`ps(igrid)%x`, `%dx`,
...) are bounds-remapped pointer views into those arrays (see
`point_at_geometry` in `src/amr/mod_amr_solution_node.fpp`) — remapping rather
than plain pointer assignment, because `surfaceC` starts at index 0 and
`dvolume`/`ds` do too when `nghostcells` is odd. Keeping the views is what lets
the large body of host code that uses `ps(igrid)%x` and friends, including
every case's `mod_usr.fpp`, work unchanged.

A block's metrics are built on the device, by `fill_geometry_device` in
`src/amr/mod_amr_solution_node.fpp`, which `alloc_node` calls once per block.
There is a single path for all three coordinate systems: a uniform block is an
analytic function of nothing but its corner and spacing in `rnode`, so every
member of `geo_t` is derived on the GPU rather than computed on the CPU and
shipped across the bus on each (re)allocation of a node — which, with AMR, is
most regrids. The per-cell expressions are fypp-specialized on `GEOM`, so a
build only carries the branch it needs.

Because the device is the producer, *all* of `geo_t` is device-resident in
every build, including members no kernel reads (`dx`, and in a Cartesian build
the metrics whose values its kernels take from `rnode` instead). They have to
live where they are written. Only the positions are returned to the host by
`alloc_node` itself, since `ps(igrid)%x` is read immediately for the initial
condition, the user hooks and the output.

That residency is not free: every member is a block-sized array for all
`max_blocks` blocks a rank could ever own, allocated up front whether or not
the blocks exist. On a 12 GB card, `max_blocks = 4096` with `block_nx = 16`
puts `bgeo` plus `bgeoc` at over 4 GB before a single block is allocated, and
`initialize_vars` fails with `CUDA_ERROR_OUT_OF_MEMORY` at its `enter data`.
That is a limit of the card, not a bug — the same case runs at
`max_blocks = 2048` — but it is why a member that nothing reads does not stay
in `geo_t` on the argument that it is cheap.

Two members were dropped for exactly that reason, and neither the fine nor the
coarse representative carries them any more:

- `dsC`, the cell-face lengths. Its only reader was
  `b_from_vector_potentialA` in `src/mod_constrained_transport.fpp`, which
  needs `stagger_grid` (unsupported here) and in fact has no callers at all.
- `surface`, the cell-centre face areas. Its only reader was the
  `Stokesbased` branch of `curlvector` in `src/mod_geometry.fpp`, reachable
  only via `typecurl = 'Stokesbased'` in `&methodlist`.

The matching `state` components `s%dsC` and `s%surface` are still declared but
are never associated, so that the two routines above still compile. Both now
call `mpistop` before they would dereference them, naming what was removed —
re-enabling either means giving it back a source for its metric, not just
deleting the guard.

Two consequences of building the metrics analytically:

- **Grid stretching is rejected.** `stretch_dim` in `&meshlist` makes the cell
  spacing vary within a block, which the analytic derivation cannot express,
  so `initialize_vars` calls `mpistop` if any dimension is stretched. The
  ~1300-line host path that used to handle it (`fill_geometry_host`), and
  `get_surface_area` in `src/mod_geometry.fpp`, are gone. Upstream
  MPI-AMRVAC's stretched-grid machinery (`qstretch`, `dxfirst`, ...) still
  sits in `mod_global_parameters` and `read_par_files` but now has no
  consumer.
- `Cartesian_expansion` was already unreachable here — `set_coordinate_system`
  `mpistop`s on it unless `ndim == 1`, and `ndim` is a compile-time 3 — so the
  `usr_set_surface` hook in `src/mod_usr_methods.fpp` is likewise now unused.

### Getting the geometry back to the host

Nothing is copied back by `alloc_node`. That routine runs on every regrid,
mostly for blocks no host routine will look at, so the fetch is placed at the
sites that actually read the geometry instead. Two routines in
`src/mod_geometry.fpp` do it:

- `sync_positions_host` — just the cell-centre positions, which are what host
  code asks for most often. Takes an optional `igrid` for one block, and
  defaults to every grid in `igrids`.
- `sync_geometry_host` — the whole of `geo_t`, positions and metrics alike.
  Takes no arguments: it always covers every grid in `igrids`.

Neither is guarded by a "host is already up to date" flag, deliberately: the
geometry is per block, so any single flag would have to be cleared after
refreshing whatever set of blocks happened to be listed at the time — after
which the next regrid or `selectgrids` can put a block that was not in that
set straight back in front of host code. The bulk form is not cheaper per
block either (the arrays are strided by `igrid`, so it issues one `update` per
block regardless); pick whichever form matches the caller's shape, and only
avoid calling the bulk form from inside a per-block loop.

`sync_geometry_host` is called where `bg%w` is pulled back for output —
`saveamrfile`, next to its `!$acc update host(ps(igrid)%w)` loop — which
covers the log, the snapshot, the slices, the collapsed views and
`autoconvert`; standalone `convert` mode bypasses `saveamrfile`, so
`src/agile.fpp` syncs before `generate_plotfile` too. Individual output
readers (`get_volume_average`, `calc_grid`, ...) therefore do not sync
themselves. `fix_conserve` is the one non-output metric reader, and calls it
itself.

`sync_positions_host` is called wherever host code reads `ps(igrid)%x`:
`initial_condition` (per block, in `initlevelone`'s loop, in
`coarsen_grid_siblings`, and in `refine_grids` for each newly created child),
`modify_IC`, `getintbc` (when `usr_internal_bc` is
set), `process`/`process_advanced` (when the `usr_process_*_grid` hooks are
set), `usr_modify_output` in `src/agile.fpp`, `prolong_grid`/`prolong_2nd`
(when `prolongprimitive` or `fix_small_values` is on), and `alloc_node` itself
just before `set_B0_grid`/`phys_set_equi_vars`, which read positions on the
host.

**Adding a host reader of `ps(igrid)%x` or of any metric means adding the
matching sync call.** A missed one reads whatever the host copy held for the
block that previously occupied that `igrid` slot — plausible-looking numbers,
no crash. Readers in code paths this fork does not currently exercise —
`mod_thermal_emission`, `mod_particle_*`, `mod_point_searching`,
`mod_interpolation`, `mod_functions_bfield`, `mod_constrained_transport` —
will need the same treatment when they are brought back into use.

Current limits of the curvilinear (spherical and cylindrical) support:

- `phys = 'hd'`, `phys = 'mhd'`, `phys = 'srhd'` and `phys = 'ffhd'` all
  implement `addsource_geometry` for both `geometry = 'spherical'` and
  `geometry = 'cylindrical'`, each guarded by its own `#:if GEOM == '...'`
  branch in the physics module's `mod_*_templates.fpp` (a build only compiles
  the branch matching its `GEOM` define). MHD's and SRHD's curvature terms
  (`src/mhd/mod_mhd_templates.fpp`, `src/srhd/mod_srhd_templates.fpp`) were
  ported from upstream MPI-AMRVAC's `mhd_add_source_geom` (cell-centred
  GLM-MHD branch; this fork has no `stagger_grid`/constrained-transport
  support) and `srhd_add_source_geom`, with the isotropic pressure-like terms
  (`ptot`/`psi` for MHD, `p` for SRHD) rewritten to use the discrete `dAdV`
  well-balancing factor instead of the continuous `2/r`, `cot(theta)/r`
  (spherical) or `1/r` (cylindrical) prefactors upstream uses directly, for
  consistency with this fork's HD implementation. For cylindrical, `x(1)*
  dAdV(1)` works out to exactly `1/r` rather than merely converging to it,
  because a cylindrical radial face area is linear in `r` (unlike spherical's
  `r^2`), so the well-balancing rewrite is exact there, not just consistent
  in the continuum limit. SRHD's primitive `mom(:)` slot holds the spatial
  four-velocity `u^i = lfac*v^i` rather than `v^i` itself (see
  `to_primitive`/`to_conservative` and `src/srhd/mod_con2prim.fpp`'s
  `xi = tau + D + p`, `v^2 = S^2/xi^2`), so the curvature terms use the
  primitive `xi_`/`lfac_` auxiliaries directly rather than reconstructing the
  conserved momentum density. Cylindrical's curvature terms only touch the
  `m_r`/`m_phi` (and, for MHD, `B_r`/`B_phi`) components — `m_z`/`B_z` pick up
  none, since a cylindrical volume element's `z`-extent does not depend on
  `r`. FFHD's `addsource_geometry` (`src/ffhd/mod_ffhd_templates.fpp`) is
  empty for both spherical and cylindrical — matching upstream's
  `ffhd_add_source_geom`, which is likewise empty for every coordinate system
  — because its conserved quantities (`rho`, the field-aligned momentum
  `m_par`, the field-aligned energy) are genuine scalars advected along a
  user-supplied field direction `b-hat`, and Gauss's theorem means a scalar's
  flux divergence needs no curvature correction in any coordinate system,
  unlike a vector quantity's components. What *is* geometry-sensitive for
  `ffhd` is unrelated to `addsource_geometry`: the optional `PDIVB` source
  term (`p * div(b-hat)`, compile-time `ffhd_pdivb=T`, default off) computes
  `div(b-hat)` in `addsource_nonlocal` as a plain finite difference over
  `dx(idir)`, which is only correct on a slab-uniform mesh — it is missing
  the curvilinear metric factors a true divergence needs in
  `geometry = 'spherical'` or `geometry = 'cylindrical'`, and this fork has
  not fixed that.
- The polar axis is supported for `phys = 'hd'`, `phys = 'mhd'`,
  `phys = 'srhd'` and `phys = 'ffhd'`; see "The polar axis" below. `ffhd`
  reaches the axis by a different route than the others: its conserved
  quantities are all scalars, so they take no sign flip and the ordinary pole
  copy in `getbc` carries them unchanged, while its frozen field is not
  exchanged through `getbc` at all — `fill_nwextra_device` re-derives it
  analytically in every ghost cell, the axis ghosts included (see "The frozen
  field" below).
- AMR across curvilinear levels is untested here; prolongation in
  `src/amr/mod_refine.fpp` has the `slab_uniform` branch that uses `dvolume`,
  but `fix_conserve` is commented out in `src/mod_advance.fpp` for all
  geometries.

Validated by eight test directories, one per (physics, geometry) pair, which
is as few as the compile-time parameters allow: `phys` and `geometry` are both
fypp defines, so hd cannot share a build with mhd nor spherical with
cylindrical. Within a pair the cases differ only at run time, so they share one
build and pick their initial condition through the `setup` entry of the
`&usr_list` namelist (read in `usr_init` via `params_read_user(par_files)`, the
mechanism `tests/hd/cloud_crushing` uses):

- `uflow.par`, `setup = 'uniform'` — a uniform Cartesian state, which is an
  exact steady solution. Written in curvilinear components it is a non-trivial
  function of position, so keeping it uniform exercises every curvilinear flux
  and every geometric source term at once, and is the well-balancing test.
  What "uniform" means per physics: a Cartesian velocity for hd; that plus a
  uniform Cartesian magnetic field for mhd, which exercises the induction
  equation and the GLM `psi` too; a constant field-aligned speed along a
  uniform frozen field for ffhd; and a sub-luminal Cartesian velocity for
  srhd, whose Lorentz factor stays a single global constant because it depends
  only on `|v|`, which a rotation leaves invariant.
- `blast.par`, `setup = 'blast'` — an over-pressured sphere in gas at rest,
  deliberately off-centre in all three coordinates so every momentum component
  is exercised. It is the same setup with the velocity zeroed and a hot spot
  added, which is what lets one initial-condition routine serve both.
- `blast_amr.par` — `tests/hd/spherical` only: the same blast, Lohner-refined.
  The shell keeps moving outward in `r`, so the grid is re-made throughout the
  run rather than settling after the first regrid, which is what puts the
  curvilinear prolongation, coarsening and ghost-cell prolongation paths under
  test.

`agile.par` in each directory is the build reference: `make/config_reader.py`
takes the compile-time parameters from *that file alone*, so it has to declare
the union of what every par file in the directory needs. `specialboundary` is
a runtime logical as well as a fypp define, so `blast.par` leaves it out while
`specialbound_usr` stays compiled in. **Adding a par file to such a directory
is free; adding a directory costs a full rebuild in CI.**

These domains all stay away from the polar axis and from `r = 0`; the axis
itself is covered by the separate `*_pole` directories below. Pole and non-pole
cases are kept as separate builds deliberately, even though their compile-time
parameters would have allowed merging them too.

### The polar axis

A spherical domain may span the full `theta` from 0 to `pi`, and a cylindrical
one may start at `r = 0`. The axis is not a boundary: a ghost cell just past
`theta = 0` at azimuth `phi` *is* the interior cell at `theta` and `phi + pi`,
so it is filled from the block half a revolution away rather than by a
boundary condition. Setting it up takes three things in `agile.par`, and
`check_pole_setup` in `src/mod_geometry.fpp` stops the run with a specific
message if any of them is missing:

```fortran
 &boundlist
   typeboundary_min2 = 5*'pole'       ! the axis face(s)
   typeboundary_min3 = 5*'periodic'   ! phi, over a full turn
   typeboundary_max3 = 5*'periodic'
 /
 &meshlist
   xprobmin2 = 0.0d0, xprobmax2 = 0.5d0   ! theta from 0 to pi
   xprobmin3 = 0.0d0, xprobmax3 = 1.0d0   ! phi over 2*pi (2*pi units)
   domain_nx3 = 16, block_nx3 = 8         ! ng3(1) = 2, even
 /
```

- **`'pole'` on the axis face.** `read_par_files` expands it into ordinary
  per-variable `symm`/`asymm` entries in `typeboundary`, which is where the
  ghost-cell exchange later reads the sign from. It is all-or-nothing per
  face. Only `typeboundary_min2`/`max2` may be `'pole'` in spherical and only
  `typeboundary_min1` in cylindrical.
- **`phi` periodic over exactly one full turn.** `set_pole` refuses to
  recognise an axis at all unless `periodB(phi_)`, and the partner block is
  found by shifting the `phi` block index by `ng(1)/2` — half a revolution
  only if `phi` really spans `2*pi`. Upstream MPI-AMRVAC has no check for
  this and silently shifts by the wrong angle on a partial-turn domain;
  `check_pole_setup` does check.
- **An even number of level-1 blocks across `phi`**, for the same `ng(1)/2`
  reason. This one `set_pole` already enforced.

The sign table follows one rule in both coordinate systems: **the component
along the mirrored direction flips, the `phi` component flips, everything else
is symmetric.** For spherical that is `m_theta`, `m_phi` (and `B_theta`,
`B_phi`), with `m_r`/`B_r` symmetric — the same table upstream uses. For
cylindrical it is `m_r`, `m_phi` (and `B_r`, `B_phi`). Note the radial entry:
**upstream marks only `m_phi`/`B_phi` antisymmetric at the cylindrical axis,
which is wrong** — the axis maps `(r, phi)` to `(r, phi + pi)`, under which
`e_r(phi+pi) = -e_r(phi)` exactly as `e_phi(phi+pi) = -e_phi(phi)`, so the
radial component is odd too. With upstream's table the first ghost cell of
`tests/hd/cylindrical_pole_uniform_flow` holds `-0.9713` where the exact value
is `+0.9713`; with this one it is exact to round-off. The indices come from
`iw_mom`/`iw_mag` rather than upstream's positional `rho-mom-[e]-B`
arithmetic, because a `HYPERTC` build registers `q_` between `e_` and `mag`
and shifts every field index.

Everything outside `getbc` was already pole-aware, inherited from upstream:
`set_pole` sets `poleB`, `find_root_neighbor` does the pi-periodic root lookup,
`find_neighbor` exports a `pole(ndim)` flag and suppresses the sibling-index
flip across it, `build_connectivity` records
`neighbor_pole(i1,i2,i3,igrid)` — 0, or the dimension the pole lies across —
and `alloc_node` clears `is_physical_boundary` on a pole face so `bc_phys` is
correctly skipped there. What had to be added is the copy itself, at the six
places in `getbc` (`src/mod_ghostcells_update.fpp`) that move data between
blocks: the same-rank and MPI variants of the same-level, restricting and
prolonging paths. Each now reads `neighbor_pole` and, when it is non-zero:

- **walks the source backwards along that dimension** (`sbase`/`sstep` in the
  code), because the send range runs outward from the axis while the ghost
  range runs inward toward it;
- **multiplies by the sign** from `typeboundary(iw, iB_pole)`, where
  `iB_pole = 2*(ipole-1) + iside`;
- **does not mirror the destination offset in that dimension**, because the
  partner faces you from the same side. In the 0..3 child-offset index this
  reads as: the restricting path swaps 1 with 2 and leaves 0 and 3 alone
  (instead of swapping 0 with 3), and the prolonging path does not swap at all.

For the MPI paths the mirror and the sign are applied on the **send** side, as
upstream does, so the unpack stays a plain strided copy. The same-level send
therefore also has to ship the *destination* offset in its info record rather
than its own, since the receiver cannot see the mirror; the coarse and fine
paths already shipped explicit `inc` indices and only needed the pole variant
of the swap. The prolongation interpolation itself needs no pole branch: by
then `bgc(1)%w` already holds mirrored, signed values.

Two consequences worth knowing:

- **The timestep collapses near the axis.** `ds(3) = r*sin(theta)*dphi` goes
  to zero there, so `setdt` gives a severe CFL restriction. No axis filtering
  or ring averaging is implemented; keep pole cases at modest resolution.
- **Volume averages cannot validate the pole copy.** A cell at the axis has
  vanishing volume and the face lying on the axis has zero area, so a wrong
  ghost value is multiplied by zero on its way into the interior and the jump
  it leaves behind is clipped by the TVD limiter. Reversing *every* sign at
  the pole moves the `mean(w)` the log reports only in the fifth digit. The
  pole test cases therefore assert on the ghost cells directly rather than on
  the log: `check_pole_ghosts` in each case's `mod_usr.fpp` compares them
  against the analytic state at `it = 0`, where the interior is still exactly
  analytic and the copy of it therefore has to be exact to round-off. **A new
  pole test needs that check, not just a reference log.**

Enabled for `phys = 'hd'`, `phys = 'mhd'`, `phys = 'srhd'` and
`phys = 'ffhd'`. srhd needs no change to the copy itself: its primitive
`mom(:)` slot holds the spatial four-velocity `u^i = lfac*v^i` rather than
`v^i`, but `lfac` is a scalar and therefore invariant under the pole's
pi-rotation, so `u^i` transforms exactly like an ordinary vector and takes
the same `iw_mom`-driven sign table hd's momentum does. ffhd needs no sign
table at all: `rho`, the field-aligned momentum `m_par` and the energy are
scalars — `m_par = rho*(v·b-hat)` is a contraction of two vectors, invariant
under the proper (det +1) rotation the pole map is — so `read_par_files`
skips the `iw_mom` lines for `physics_type == 'ffhd'` and every
getbc-exchanged variable stays `symm`. The frozen field, the one genuine
vector, is handled entirely outside `getbc`; see "The frozen field" below.

Validated by eight test directories, laid out exactly like the off-axis ones
above and for the same reason — the suite's cost is dominated by compilation,
so cases that agree on the fypp defines share one build and differ only in
their par file. They are kept separate from the off-axis directories even
though their compile-time parameters would have allowed merging: a pole case
and an ordinary curvilinear case are different test families, and mixing them
would make each directory harder to reason about than the extra build is
worth.

- `tests/hd/spherical_pole` — three runs from one build, selected by the
  `setup` entry of the `&usr_list` namelist (read in `usr_init` via
  `params_read_user(par_files)`, the same mechanism
  `tests/hd/cloud_crushing` uses): `uflow.par` is a uniform Cartesian flow on
  a single level, `uflow_amr.par` the same with half the domain in `phi`
  refined so blocks meet across the axis at different levels, and
  `blast_amr.par` a Lohner-refined blast that starts off the axis and expands
  across it. `movie.par` is a longer, visualisable variant of the blast.
- `tests/hd/cylindrical_pole` — the same three runs onto the cylindrical axis
  at `r = 0`, with the hot spot placed off the axis in `r` and carried across
  it by a background velocity pointing along `-y`.
- `tests/mhd/spherical_pole` and `tests/mhd/cylindrical_pole` — the uniform
  flow with a uniform Cartesian magnetic field, which puts `B_theta`/`B_r`,
  `B_phi` and the GLM `psi` through the same treatment. The cylindrical cases
  are what pin down the radial sign discussed above.
- `tests/srhd/spherical_pole` and `tests/srhd/cylindrical_pole` — a uniform,
  sub-luminal Cartesian flow, single level, one `uflow.par` each. `mean(rho)`
  in the log is `rho0*lfac0`, not `rho0`, since srhd's conserved density is
  `D = rho*lfac`.
- `tests/ffhd/spherical_pole` and `tests/ffhd/cylindrical_pole` — a uniform
  field-aligned flow along a uniform Cartesian frozen field, single level, one
  `uflow.par` each. `check_pole_ghosts` here checks two things at a same-level
  axis neighbour: the fluid ghost cells against the analytic constant (the
  ordinary `getbc` pole copy, which for ffhd only has to deliver a constant),
  and the frozen-field ghost cells rebuilt in Cartesian against the uniform
  `b0` they must reduce to (`fill_nwextra_device` at the axis, which is
  where the sign flips a vector picks up across the pole actually come from).

`agile.par` in each directory is the build reference: `make/config_reader.py`
takes the compile-time parameters from *that file alone*, so it has to declare
the union of what every par file in the directory needs. `specialboundary` and
`refine_usr` are both runtime logicals as well as fypp defines, so an
individual par file can switch them off again while the routines they guard
stay compiled in. **Adding a par file to such a directory is free; adding a
directory costs a full rebuild in CI.**

`check_pole_ghosts` compares the ghost cells against the analytic state at
`it = 0`, where the interior is still exactly analytic and the pole copy of it
therefore has to be exact. How far that can be pushed across a *level jump*
depends on the setup, because there the ghost has been restricted or prolonged
on the way and comparing it against the analytic value at a point is only
meaningful where the field is smooth on the scale of a coarse cell:

- the `uniform` setups are smooth everywhere, so both branches are checked:
  round-off (1e-10) for a same-level pole neighbour and loose (1e-1) across a
  level jump, where the measured error is 3.4e-2 — far below the 0.18 a wrong
  sign produces.
- the `blast` setup has a discontinuous spot surface, and a coarse cell
  straddling it holds the 2:1 average of hot and cold gas, which differs from
  the analytic value at its centre by a large fraction of the jump however
  correct the copy is. Widening the radial domain for `movie.par` put a pole
  face on that surface and the check fired at exactly `13.5/4`, a quarter of
  the initial jump in `e`, on a `neighbor_fine` face. So that setup checks
  same-level neighbours only and covers the restricting and prolonging paths
  through its log instead — which it can afford to do because its log is the
  one that is actually sensitive to them: reversing the `theta` sign makes
  `compare_logs` fail on `mean(m_r)`, `mean(m_theta)` and `mean(m_phi)`, where
  the same break leaves the uniform-flow logs comfortably inside tolerance.

The lesson for a new pole test: point-versus-analytic only works where the
analytic state is smooth on the scale of a coarse cell. Where it is not, check
same-level neighbours and lean on a log that has teeth.

### The frozen field (`phys = 'ffhd'`)

ffhd advects its scalar conserved quantities along a static, user-supplied
unit vector `b-hat`, stored in the `w` slots `b1,b2,b3` (registered by
`var_set_extravar` in `src/ffhd/mod_ffhd_templates.fpp`, so they are the
`nwextra` variables). It is a pure function of position and never a boundary
condition.

This is handled by a generic mechanism, not by ffhd-specific code in the AMR
core. `config_schema.toml`'s `implies` gives `phys = 'ffhd'` the
`FILL_NWEXTRA_ANALYTIC` compile flag, which switches on `fill_nwextra_device`
(`src/amr/mod_amr_solution_node.fpp`) and its call sites. That routine, on the
device, evaluates a by-name user hook

```fortran
pure subroutine usr_set_nwextra(x, wextra)   ! x(1:ndim) in, wextra(:) out
  !$acc routine seq                          ! for ffhd: b-hat, already unit
```

for every cell of every block — interior, inter-block ghosts,
physical-boundary ghosts and polar-axis ghosts alike — after every change to
the grid, and writes the result verbatim into the extra slots (so a case must
return `b-hat` already normalised; `to_spherical_unit`/`to_cylindrical_unit`
do). The hook is called by name, like `usr_refine_grid` and `gravity_field`,
so it is a compile-time dependency of any `FILL_NWEXTRA_ANALYTIC` build. The
call sites are `initlevelone`, `modify_IC` and the end of `amr_coarsen_refine`
(each looping over `igrids`); `initonegrid_usr` and `usr_special_bc` must
**not** set `b1,b2,b3`.

Consequences:

- `nwgc = nwflux` for ffhd (not `nwflux + nwextra`): the frozen field is not
  in the ghost exchange, so `getbc`, prolongation and coarsening no longer
  need to carry it and `typeboundary` needs no rows for it. RK substeps keep
  it because they copy `1:nw` wholesale and the flux update never writes past
  `nwflux`.
- **The polar axis needs no special handling for the frozen field.** `bgeo%x`
  in the ghost layer beyond `theta = 0` (or `r = 0`) carries the mirrored
  coordinate — a negative `theta`, or a negative `r` — so evaluating the
  user's analytic field there reproduces on its own the sign flips a vector
  picks up across the axis (`b_r` symmetric, `b_theta`/`b_z` and `b_phi`
  antisymmetric). This is the same identity the `analytic_state` reference in
  the hd/srhd pole tests relies on.
- The old design exchanged `b` through `getbc` (`nwgc = nwflux + nwextra`) but
  filled the *physical-boundary* ghost cells only if the case happened to fill
  the whole `ixI` in `initonegrid_usr` or rewrite `b` in `usr_special_bc`. A
  case with `cont` boundaries and an interior-only IC read uninitialised
  memory there; `tests/ffhd/spherical` and `tests/ffhd/cylindrical` had
  exactly that latent bug in their `blast.par` runs, which is why their
  `correct_output/blast.log` changed when the frozen field moved to
  `usr_set_nwextra`.

## Writing a new simulation case (`mod_usr.fpp`)

Every case directory has a `mod_usr.fpp` defining `module mod_usr` with a
`usr_init()` subroutine that:
1. registers callbacks by assigning the procedure pointers declared in
   `src/mod_usr_methods.fpp` (e.g. `usr_init_one_grid`, `usr_special_bc`,
   `usr_source`, `usr_process_grid`, `usr_refine_threshold`, ...),
2. calls `phys_activate()` last.

The coordinate system comes from `geometry` in `agile.par` (see "Coordinate
systems" above); `usr_init()` does not need to call `set_coordinate_system`
itself. Older cases may still call it explicitly — the call is idempotent as
long as the name agrees with `geometry` — but new cases shouldn't.

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
