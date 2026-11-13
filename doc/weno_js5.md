# WENO-JS5 finite-volume reconstruction

`limiter = 20*'weno5'` selects component-wise Jiang--Shu reconstruction of
conserved cell averages. Each reconstructed face state is converted to
primitive variables and passed to the existing LLF (`tvdlf`), HLL or HLLC
solver. The supported initial scope is Cartesian HD/MHD with the energy
equation. Existing solver/physics restrictions still apply. Geometry other
than Cartesian, other physics and `flux_adaptive_diffusion` are rejected.

```fortran
&methodlist
  limiter = 20*'weno5'
  weno_js_epsilon = 1.d-6
  flux_scheme = 20*'hll'
  time_stepper = 'threestep'
/
```

Mixed MUSCL/WENO levels are allowed. The halo width is raised globally to at
least three; setting `nghostcells=4` is also supported. Unsupported limiter
names now fail explicitly rather than silently returning first-order states.
The default limiter is unchanged. Rebuild after updating the source.

## Mathematical contract and sources

The reconstruction is JS5, not WENO-Z, mapped WENO, WENO-NM or point-value
interpolation. For five conserved cell averages `v(1:5)`, with `v(3)` the donor
cell, the three right-face candidates are

\[
q_0=(2v_1-7v_2+11v_3)/6,\quad
q_1=(-v_2+5v_3+2v_4)/6,\quad
q_2=(2v_3+5v_4-v_5)/6.
\]

\[
\begin{aligned}
\beta_0&=\tfrac{13}{12}(v_1-2v_2+v_3)^2+
          \tfrac14(v_1-4v_2+3v_3)^2,\\
\beta_1&=\tfrac{13}{12}(v_2-2v_3+v_4)^2+
          \tfrac14(v_2-v_4)^2,\\
\beta_2&=\tfrac{13}{12}(v_3-2v_4+v_5)^2+
          \tfrac14(3v_3-4v_4+v_5)^2.
\end{aligned}
\]

With `d=(0.1,0.6,0.3)`, `alpha=d/(epsilon+beta)^2`, the normalized weights are
`omega=alpha/sum(alpha)` and the result is `sum(omega*q)`. The implementation
multiplies every alpha by the same `min(epsilon+beta)^2` before normalization;
this cancels algebraically and avoids overflowing inverse squared weights.
It does not alter epsilon or the JS definition.

`weno_js_epsilon` is finite, positive, defaults to `1e-6`, and acts on conserved
variables in code units. There is no automatic amplitude or grid-spacing
normalization. JS behavior at critical points and with very small amplitudes
is retained; universal fifth-order convergence at every critical point is not
an acceptance condition.

Sources consulted:

- Jiang & Shu (1996), *Efficient Implementation of Weighted ENO Schemes*,
  JCP 126, 202–228, DOI [10.1006/jcph.1996.0130](https://doi.org/10.1006/jcph.1996.0130).
  Tables I–II specify candidates and optimal weights; equations (2.15)–(2.16)
  specify normalized weights; section 3, starting with (3.1), specifies the
  derivative-integral smoothness indicators. The numerical-test epsilon is
  `1e-6`. Zotero item `QMU2G4TR`, attachment `LDGTZQ7S`.
- Shu (2009), *High Order Weighted Essentially Nonoscillatory Schemes for
  Convection Dominated Problems*, SIAM Review 51, 82–126, DOI
  [10.1137/070679065](https://doi.org/10.1137/070679065).
  Sections 2.2 and 3.1–3.4 distinguish reconstruction, FV/FD input-output
  meanings and multidimensional face integration. Zotero item `TW7H5PTC`,
  attachment `L45UVT6K`.
- Local AMRVAC reference: `2607_FIP/simulations/amrvac/src/limiter/mod_weno.t`,
  `WENO5limiter`, `var=1`, inspected 2026-11-13. SHA256
  `d26cee760f9ed6abde0643d44e9580b7a3bd184b6a0db72de4195f4b7043db23`.
  This is a source snapshot without Git metadata, not a verified upstream
  release. Its candidate coefficients, JS indicators and mirrored right
  stencil agree with the formulas above. Its epsilon is `1e-18`, deliberately
  not adopted. Its `fix_limiter` postprocessing is not transplanted.

The 1996 paper's FD flux-splitting solver is not copied wholesale. The shared
JS reconstruction is used in the FV sense explained in the 2009 review.

## Kernel and face indexing

Fypp generates separate reconstruction/flux/layout specializations. MUSCL
uses four-point primitive single-face input, WENO/MP5 six-point conserved
single-face input. Selection is outside the cell loop and filters blocks by
their level's reconstruction. Both `per_cell` and `shared` use the same face
adapter. For face `i+1/2`, input covers `i-2:i+3`:

| Side | Input to the left-biased scalar function | Donor |
|---|---|---|
| L | `u(1:5)` | `i` (`u(3)`) |
| R | `u(6:2:-1)` | `i+1` (`u(4)`) |

The lower face uses the same adapter shifted by one cell. The union of both
faces still spans `i-3:i+3`, requiring three halos.
Conserved means are read directly from the input stage. The original
primitive block buffer is retained only for MUSCL or enabled source terms.
Nonlocal sources retain their five-point buffer, compact sources their
three-point directional buffers. Neither is enlarged to seven points.

Nonfinite or nonpositive reconstructed density is rejected before conversion.
After conversion the existing density/pressure floors and auxiliary-state
repair are applied; a still-invalid face state falls back to its donor cell.
Production builds do not accumulate or report diagnostic fallback counts.
This is not a positivity-preserving update of
cell averages; donor cells themselves must be physically admissible, and the
usual CFL and physical-model constraints still apply.

## Halo and AMR audit

No new halo allocator or communication protocol is introduced. The existing
parameterized infrastructure was checked at these points:

- `read_par_files`: radius selection, physics extra layers, staggered even-width
  adjustment, host/device update and block-size validation.
- `initialize_vars`, `set_geometry_index_ranges`: full and coarse extents,
  including the extra geometry layer for an odd halo.
- `init_comm_types`, `init_bc`: MPI datatypes, sibling send/receive extents,
  restricted and prolonged ranges, and coarse-halo ceiling division.
- `getbc`, physical boundary filling and refinement/coarsening: three/four
  layers, actual coarse/fine interfaces and diagonal/full-buffer paths.
- `advect`/`advect1`: stage input and output ghost refresh around FV updates;
  existing diagonal requirements remain controlled by the physics/source path.
- FV boundary flux storage and reflux: physical interior indices stay derived
  from `nghostcells`, not fixed offsets.

A block must contain at least the halo width in each direction; AMR requires
at least twice that width. Existing even/minimum-four constraints still apply.
User special-boundary routines must fill every index in their supplied range.
The framework cannot repair a user callback that hard-codes two layers.

## Accuracy and validation limits

The verified fifth-order claim concerns the uniform Cartesian one-dimensional
spatial operator and smooth face reconstruction. Integration tests represent
one-dimensional profiles on thin 3D boxes, as in AGILE's existing shock tests.
SSPRK3 is third order in time; convergence runs decrease dt sufficiently and
repeat with dt halved.

This change does not implement characteristic reconstruction, high-order AMR
prolongation, multidimensional face quadrature or high-order source quadrature.
AMR and multidimensional configurations are supported operationally and tested
for conservation and seams; their entire discretization is not claimed to be
fifth order. Curve-coordinate WENO coefficients are not approximated by
Cartesian coefficients: such configurations are rejected.

See [WENO tests](../tests/weno_js5/README.md) and
[solver integration tests](../tests/mp5/README.md). GPU runtime validation is
still outstanding; CPU evidence does not establish device correctness.
