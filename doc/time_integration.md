title: Time integration

# Choosing a scheme {: #ti_choose }

The time integration scheme is a **run-time** choice, made with two entries
of `&methodlist` in the parameter file:

```fortran
&methodlist
   time_stepper    = 'fourstep'   ! 'onestep' | 'twostep' | 'threestep' |
                                  ! 'fourstep' | 'fivestep'
   time_integrator = 'ssprk4'     ! depends on time_stepper, see below
/
```

`time_stepper` fixes the number of stages (flux evaluations) per time step,
and `time_integrator` picks one of the schemes with that many stages. Leaving
`time_integrator` out, or setting it to `'default'`, selects the default for
that `time_stepper`, marked in bold below. Leaving out `time_stepper` as well
gives `'twostep'`, so the overall default is `Predictor_Corrector`, as in
MPI-AMRVAC.

AGILE implements every explicit scheme of MPI-AMRVAC:

| `time_stepper` | `time_integrator` | variant | method | order | stages | time levels | SSP | CFL limit |
|---|---|---|---|---|---|---|---|---|
| `onestep` | **`Forward_Euler`** (default) | | forward Euler | 1 | 1 | 2 | yes | 1 |
| `twostep` | **`Predictor_Corrector`** (default) | | explicit midpoint | 2 | 2 | 2 | no | |
| `twostep` | `ssprk2` | | Heun, SSPRK(2,2) | 2 | 2 | 2 | yes | 1 |
| `twostep` | `RK2_alfa` | `rk2_alfa` | two-stage RK family | 2 | 2 | 2 | no | |
| `threestep` | **`ssprk3`** (default) | `ssprk_order = 3` (default) | SSPRK(3,3), Shu-Osher | 3 | 3 | 3 | yes | 1 |
| `threestep` | `ssprk3` | `ssprk_order = 2` | SSP(3,2) | 2 | 3 | 3 | yes | 2 |
| `threestep` | `RK3_BT` | `rk3_switch` | third-order Butcher tableaus | 3 | 3 | 4 | no | |
| `fourstep` | **`ssprk4`** (default) | `ssprk_order = 3` (default) | SSPRK(4,3), Spiteri-Ruuth | 3 | 4 | 3 | yes | 2 |
| `fourstep` | `ssprk4` | `ssprk_order = 2` | SSP(4,2) | 2 | 4 | 3 | yes | 3 |
| `fourstep` | `rk4` | | classical RK4 | 4 | 4 | 4 | no | |
| `fivestep` | **`ssprk5`** (default) | `ssprk_order = 3` (default) | SSPRK(5,4), Gottlieb's coefficients | 4 | 5 | 4 | yes | about 1.5 |
| `fivestep` | `ssprk5` | `ssprk_order = 2` | SSPRK(5,4), Spiteri-Ruuth's coefficients | 4 | 5 | 4 | yes | about 1.5 |

The *variant* column names the parameter that selects among a scheme's
versions. Each [`ssprk_order`](#ti_ssprk_order) value of the SSP integrators
is a row of its own; `rk2_alfa` and `rk3_switch` are described
[below](#ti_variants).

*SSP* marks the strong-stability-preserving schemes: each is a convex
combination of forward-Euler steps, so it keeps the total-variation and
positivity properties of the spatial discretization up to the CFL limit
given, which is relative to that of a single forward-Euler step. The others
carry no such guarantee, so no such limit is listed for them. AGILE does
**not** rescale `courantpar` per scheme, so to profit from a larger limit,
raise `courantpar` yourself.

## ssprk_order {: #ti_ssprk_order }

```fortran
&methodlist
   time_stepper    = 'fourstep'
   time_integrator = 'ssprk4'
   ssprk_order     = 2           ! 3 (default) or 2
/
```

`ssprk_order` in `&methodlist` selects the coefficient set of the SSP
integrators `ssprk3`, `ssprk4` and `ssprk5`. It takes the values 3, the
default, and 2; any other value stops the run at startup. It has no effect
on the other integrators.

- For `ssprk3` and `ssprk4`, 2 selects a *different method*: SSP(3,2) or
  SSP(4,2) in place of SSPRK(3,3) or SSPRK(4,3). These have the same number
  of stages and time levels but are only second order. What they buy in
  exchange is a larger CFL limit: 2 instead of 1 for three stages, and 3
  instead of 2 for four.
- For `ssprk5` both values give the same fourth-order method, SSPRK(5,4).
  They select between two published sets of its coefficients, Gottlieb's (3)
  and Spiteri and Ruuth's (2), which agree to about nine significant digits.
  The value 3 is a label here, not an order.

The switch is applied in `read_par_files`, which fills the coefficients
`rk_beta*`, `rk_alfa*` and `rk_c*` from it. `advect` runs one stage sequence
per integrator for both values, because each pair of variants has the same
structure and only its coefficients differ. `ssprk_order = 2` for `ssprk3`
and `ssprk4` is tested by `ssp32.par` and `ssp42.par` in
`tests/hd/KH3D-SMR-integrators`.

## Other variant parameters {: #ti_variants }

- `rk2_alfa` (default 0.5) is the parameter of `RK2_alfa`, in (0, 1]. 0.5 is
  the explicit midpoint rule, the same method as `Predictor_Corrector`, and 1
  is Heun's method, the same method as `ssprk2`.
- `rk3_switch` (default 3) selects the Butcher tableau of `RK3_BT`:
  1 Ralston, 2 Wray, 3 Heun, 4 Nyström. All are third order.

The IMEX schemes of MPI-AMRVAC (`IMEX_Euler`, `IMEX_SP`, `IMEX_Midpoint`,
`IMEX_Trapezoidal`, `IMEX_222`, `IMEX_ARS3`, `IMEX_232`, `IMEX_CB3a`) are
still recognised by the parameter reader, but AGILE stops at startup if one
is requested. There is nothing for them to integrate implicitly yet: no
physics module provides an implicit update. The same message is given for
any scheme not in the table above.

The run prints the scheme it uses at startup, for example

    time_integrator:  ssprk4 (4 stages, 3 time levels)

# Stages and time levels {: #ti_levels }

A scheme has two different sizes, and only the second one costs memory:

- **stages** (`nstep`): flux evaluations per time step, which sets the cost
  of a step.
- **time levels** (`nlevel`): solution copies the scheme has to hold at the
  same time, its registers.

The solution of all blocks lives in `bg(1:nlevel)` (a `block_grid_t` from
`mod_physicaldata`), allocated once at startup and resident on the GPU.
Each time level is `nw` times the block size (ghost cells included) times
`max_blocks` double-precision values, so on a GPU it often decides how large
`max_blocks` can be. In the `hd/KH3D-SMR` test, for instance, one level is
819 MB per MPI rank, and `ssprk2` peaks 780 MiB lower on the device than
`ssprk3`.

The two numbers differ in both directions: `ssprk4` runs 4 stages in 3
levels, while `RK3_BT` needs 4 levels for its 3 stages. Both come from one
table, `time_levels` in `mod_global_parameters`, which `read_par_files`
consults after reading `&methodlist`. Everything that depends on the number
of levels follows from `nlevel`: the allocation of `bg` in `initialize_vars`,
and the state views `ps1`..`ps4` that `alloc_node` points at
`bg(2..nlevel)`.

# How a scheme is written {: #ti_code }

Each scheme is a `case` of `advect` in `mod_advance`, written in terms of two
operations on the time levels. The first is the stage itself,

    advect1(dtfactor, ..., qtC, ka, qt, kb):
        level kb += dtfactor*dt * R(level ka)

where `R` is the full spatial operator: the flux divergence and the unsplit
sources. It is followed by refluxing at coarse-fine interfaces and by the
ghost-cell exchange. The two levels must differ, because the flux stencil
reads neighbouring cells. The second operation is a pointwise linear
combination between stages:

    bg_assign(bgd, a1, bg1 [, a2, bg2 [, a3, bg3]])    bgd = a1*bg1 + ...
    bg_accum (bgd, s, a1, bg1 [, ...])                 bgd = s*bgd + a1*bg1 + ...

These are generic names over GPU kernels generated for one to three sources.
They come in two shapes so that a destination is never also passed as a
source, and is never read when its coefficient is zero.

Level 1 holds the solution, and `advect` starts by copying it into level 2.
For instance, `ssprk3` (Shu and Osher's form) is

    advect1(1,   1 -> 2)              ! level 2 = u + dt R(u)
    bg_assign(bg(3), 3/4, bg(1), 1/4, bg(2))
    advect1(1/4, 2 -> 3)              ! level 3 = 3/4 u + 1/4 (u1 + dt R(u1))
    bg_accum(bg(1), 1/3, 2/3, bg(3))
    advect1(2/3, 3 -> 1)              ! u = 1/3 u + 2/3 (u2 + dt R(u2))

and needs three levels because `u` itself is still read in the last
combination.

Refluxing runs after every stage, so that the coarse side of every
coarse-fine interface is corrected with the fluxes of that stage, except for
the predictor of `Predictor_Corrector`, as in MPI-AMRVAC. It is skipped
altogether when the mesh has a single level. MPI-AMRVAC's `typepred1`, a
separate flux scheme for that predictor, is not supported, so the predictor
uses `flux_scheme` like every other stage.

## Adding a scheme {: #ti_adding }

1. Give it a row in `time_levels` (`mod_global_parameters`), with its number
   of stages and the number of distinct levels its sequence touches.
2. Write its sequence as a `case` of `advect`, from `advect1`, `bg_assign`
   and `bg_accum`.
3. Add it to the list of accepted integrators in `read_par_files`, after
   the `time_stepper` parsing, along with any coefficients it needs.
4. Add a par file for it to `tests/hd/KH3D-SMR-integrators`, with a
   reference log.

# Tests {: #ti_tests }

All schemes are tested on the statically refined case `hd/KH3D-SMR`. The
refinement is what puts refluxing at every stage under test, and the
periodic domain makes mass and momentum conserved quantities, so a stage
that is refluxed wrongly shows up as drift.

- The basic suite runs `ssprk3`, `ssprk4` and the default
  (`Predictor_Corrector`) in `tests/hd/KH3D-SMR`.
- The extended suite runs every scheme and variant in
  `tests/hd/KH3D-SMR-integrators`: the nine schemes, `ssprk_order = 2` for
  `ssprk3` and `ssprk4`, and the Wray and Nyström tableaus of `RK3_BT`.

Each scheme's sequence was also checked for its order of accuracy, by
replaying it on a nonlinear ordinary differential equation and on the
Kelvin-Helmholtz case with a fixed time step.
