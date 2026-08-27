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
- The polar axis is supported for `phys = 'hd'` and `phys = 'mhd'` only; see
  "The polar axis" below. `srhd` and `ffhd` must keep the domain away from
  `theta = 0` and `theta = pi` (spherical) and from `r = 0` (cylindrical), and
  `check_pole_setup` stops them loudly if they do not.
- AMR across curvilinear levels is untested here; prolongation in
  `src/amr/mod_refine.fpp` has the `slab_uniform` branch that uses `dvolume`,
  but `fix_conserve` is commented out in `src/mod_advance.fpp` for all
  geometries.

Validated by `tests/hd/spherical_uniform_flow` (a uniform Cartesian flow
written in spherical components, which converges at second order),
`tests/hd/spherical_blast`, `tests/mhd/spherical_uniform_flow` (the same
uniform-flow idea extended with a uniform Cartesian magnetic field, to
exercise the induction equation's and GLM psi's curvature terms too),
`tests/mhd/spherical_blast`, `tests/srhd/spherical_uniform_flow` (a uniform,
sub-luminal Cartesian flow, which also converges at second order — the
Lorentz factor only depends on `|v|`, which a spherical rotation leaves
invariant, so it stays a single global constant), `tests/srhd/spherical_blast`,
`tests/ffhd/spherical_uniform_flow` (a uniform state along a uniform frozen
field, which also converges at second order — this test is what caught a
pre-existing, geometry-independent bug: `phys_init` set `nwgc=nwflux`, but
the frozen field `iw_b1`/`iw_b2`/`iw_b3` is registered via
`var_set_extravar`, outside `nwflux`, so it was silently excluded from
inter-block ghost-cell exchange whenever the field is not spatially uniform;
fixed to `nwgc=nwflux+nwextra`, matching how `srhd` already includes its own
auxiliaries via `nwgc=nwflux+nwaux`), and `tests/ffhd/spherical_blast`.

The cylindrical counterparts follow the same pattern, one `cylindrical_*`
test directory per `spherical_*` one: `tests/hd/cylindrical_uniform_flow`
(the analogous uniform Cartesian flow written in cylindrical `(r, z, phi)`
components — only the radial and azimuthal velocity are non-trivial
functions of `phi`, since the axial component stays literally constant;
confirmed at second-order convergence by doubling the resolution),
`tests/hd/cylindrical_blast`, `tests/mhd/cylindrical_uniform_flow` (adding a
uniform Cartesian magnetic field, also confirmed at second order),
`tests/mhd/cylindrical_blast`, `tests/srhd/cylindrical_uniform_flow` (a
uniform, sub-luminal Cartesian flow, also confirmed at roughly second-order
convergence), `tests/srhd/cylindrical_blast`,
`tests/ffhd/cylindrical_uniform_flow` (a uniform state along a uniform
frozen field), and `tests/ffhd/cylindrical_blast`.

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
  pole test needs that check, not just a `correct_output/test.log`.**

Validated by `tests/hd/spherical_pole_uniform_flow` and
`tests/mhd/spherical_pole_uniform_flow` (a uniform Cartesian flow, and field,
on a domain running onto both poles), `tests/hd/cylindrical_pole_uniform_flow`
and `tests/mhd/cylindrical_pole_uniform_flow` (the same onto the cylindrical
axis — these are what pin down the radial sign discussed above), and
`tests/hd/spherical_pole_amr`, which refines half the domain in `phi` so that
blocks meet across the axis at different levels and the restricting and
prolonging pole paths are exercised too. The AMR case checks its pole ghosts
at a loose tolerance, since across a level jump they carry the scheme's own
second-order error (measured at 3.4e-2) rather than being exact; that is still
far below the O(1) error a wrong sign or offset produces.

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
