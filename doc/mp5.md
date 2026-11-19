# MP5 finite-volume reconstruction

Select `limiter='mp5'` in `&methodlist` (set every enabled AMR level as needed).
The default limiter is unchanged. `mp5_tolerance=1.d-10` is finite and nonnegative;
zero selects the strict early-acceptance test used by the inspected AMRVAC code.
It is a tolerance on a product of reconstructed-variable differences (variable
units squared), not the WENO smoothness epsilon. Alpha is fixed at four.

The scalar algorithm follows Suresh & Huynh (1997), *Accurate
Monotonicity-Preserving Schemes with Runge–Kutta Time Stepping*, JCP 136, 83–99:
Eq. (2.1) gives the fifth-order candidate; (2.27) the four-argument minmod
curvature; (2.24), (2.28), (2.29) the interval; (2.26) its median projection;
(2.30) permits early acceptance with the above tolerance. AMRVAC source cross-check:
`src/limiter/mod_mp5.t`, commit `6f5e6f10e5da5719dfc160da4591ae369450b603`.
The independent Python reference solves cell-integral polynomial constraints;
the AMRVAC scalar routine is separately preprocessed and compiled, not used as
an oracle for correctness.

MP5 reconstructs conserved cell averages component by component, converts the
face states to primitives and invokes the existing LLF/HLL/HLLC solver. Existing
face-state repairs and donor fallback remain in place. Production code does
not accumulate or report diagnostic fallback counts.
This is not the paper's complete characteristic-reconstruction/Roe Euler scheme.

Supported: Cartesian HD/MHD with an energy equation and the existing flux
solver's physical restrictions. Other physics, curvilinear geometry, no-energy
models and FLUX_AD are rejected. At least three ghost layers are selected
automatically; wider physics stencils and AMR stagger-grid even-width rules
still apply. MUSCL uses four-point single-face input; WENO/MP5 use distinct
six-point conserved single-face kernels selected outside the cell loop.
Both `fv_face_flux='per_cell'` (default) and `'shared'` are supported.
The former evaluates each cell's two faces separately; the latter caches
each face once per block. Halo width remains three for either layout. Mixed per-level limiters
are supported. AMR prolongation is unchanged and is not globally fifth order.

For alpha=4, the paper's sufficient scalar monotonicity CFL bound is 0.2.
CFL=0.3 is a separately tested practical configuration, not a general guarantee
for nonlinear systems. Smooth extrema can retain high accuracy under MP5;
positivity and nonoscillatory behavior for every under-resolved system are not
guaranteed by the scalar theorem.

Validation commands: [MP5 tests](../tests/mp5/README.md). CPU results do not
establish GPU correctness. Scalar kernels compiled with GNU OpenACC/OpenMP;
full-program builds on the development environment were blocked by GPU utility
and derived-type mapping errors. GPU execution remains unverified.

The smooth spatial operator achieves approximately fifth order. One finest
short-time MP5 test reached errors near 1e-12 and failed a 10% time-step-halving
criterion; additional step-count experiments supported roundoff accumulation.
This limitation does not change the default tolerance or acceptance criterion.
