title: Coordinate systems

# Available geometries {: #geom_overview }

AGILE supports three coordinate systems: Cartesian (x, y, z), spherical
(r, θ, φ) and cylindrical (r, z, φ). All three are always 3D -- `ndim` is a
compile-time constant of 3 -- so a spherical or cylindrical run is the full
3D system, never an axisymmetric or 1D/2D reduction.

The coordinate system is a **compile-time** choice: the finite-volume
kernels are specialized for it, so that the (much cheaper) Cartesian kernel
stays free of curvilinear bookkeeping. It is selected via `geometry` in
`&meshlist` of the parameter file:

```fortran
&meshlist
   geometry = 'spherical'   ! or 'cylindrical', or 'Cartesian' (the default)
/
```

Two conventions to know when setting up a spherical or cylindrical case:

* **Angular domain bounds are given in units of 2π.** So for spherical,
  `xprobmin2 = 0.125d0`, `xprobmax2 = 0.375d0` places θ between π/4 and
  3π/4. For cylindrical, only the φ bound (`xprobmin3`/`xprobmax3`) gets
  this treatment; the z bound (`xprobmin2`/`xprobmax2`) is a plain length,
  exactly like r.
* **Axis order** is (r, θ, φ) for spherical and (r, z, φ) for cylindrical --
  in both cases φ is the third coordinate direction.

## Why curvilinear coordinates need extra source terms {: #geom_why_sources }

AGILE discretizes the generic conservation law

\begin{equation}
    \partial_t \mathbf{U} + \nabla \cdot \mathbf{F}(\mathbf{U}) = \mathbf{S}(\mathbf{U},\nabla \mathbf{U},\mathbf{x},t)
\end{equation}

(see [Physics modules and equations](equations.html)) with a finite-volume
update that divides the net flux through a cell's faces by the cell volume.
In spherical or cylindrical coordinates the face areas and cell volume
depend on position (e.g. a spherical shell's outer face is larger than its
inner face), so this ratio is what actually implements the curvilinear
divergence for scalar quantities -- no extra term is needed for the mass,
energy, or (for `ffhd`) the field-aligned quantities, since these are
advected as genuine scalars and Gauss's theorem then guarantees the finite
volume update alone reproduces \( \nabla \cdot \mathbf{F} \) correctly in
any coordinate system.

Momentum (and, for `mhd`, the magnetic field) is different: it is a
*vector* quantity, and its flux is a tensor built from a position-dependent
(r, θ, φ) or (r, z, φ) basis. Writing the vector Laplacian/divergence out in
curvilinear coordinates leaves over extra terms beyond the flux-divergence
that the finite-volume update captures -- physically, these account for how
the basis vectors themselves change from one face of a cell to the other.
AGILE adds these leftover terms explicitly as a per-cell geometric source,
in the `addsource_geometry` routine of each physics module
(`src/{hd,mhd,srhd,ffhd}/mod_*_templates.fpp`); a physics module that does
not implement one aborts at startup rather than silently omitting these
terms. The explicit expressions below are exactly what these routines
evaluate (rewritten into physical notation matching
[Physics modules and equations](equations.html); see the source files for
the exact discretization, in particular how the isotropic-pressure-like
prefactor is evaluated from the actual discrete face areas rather than
plugging in the continuous coordinate factor directly, so that a uniform
pressure produces no net force to the same precision the code otherwise
resolves -- for cylindrical this recovers the continuous prefactor exactly,
since the radial face area of a cylindrical cell is linear in r).

# Geometric source terms {: #geom_sources }

The equations below give the additional source term \( \mathbf{S}_{\mathrm{geom}} \)
that a curvilinear run adds to the corresponding equation from
[Physics modules and equations](equations.html); everywhere else the
equations are unchanged. `ffhd` needs none, for any geometry (see above), so
it is omitted below beyond noting \( \mathbf{S}_{\mathrm{geom}} = 0 \).

## Spherical (r, θ, φ) {: #geom_spherical_sources }

### Hydrodynamics: hd {: #geom_spherical_hd }

\begin{eqnarray}
    S_{\mathrm{geom}}[m_r] &=& \frac{2p+\rho\left(v_\theta^2+v_\varphi^2\right)}{r} \,,\\
    S_{\mathrm{geom}}[m_\theta] &=& \frac{p\cot\theta+\rho\left(v_\varphi^2\cot\theta-v_r v_\theta\right)}{r} \,,\\
    S_{\mathrm{geom}}[m_\varphi] &=& -\frac{\rho v_\varphi\left(v_r+v_\theta\cot\theta\right)}{r} \,.
\end{eqnarray}

### Magnetohydrodynamics: mhd {: #geom_spherical_mhd }

With the total (thermal + magnetic) pressure \( p_{\mathrm{tot}} = p + \frac{1}{2}B^2 \):

\begin{eqnarray}
    S_{\mathrm{geom}}[m_r] &=& \frac{2p_{\mathrm{tot}}+\rho\left(v_\theta^2+v_\varphi^2\right)-\left(B_\theta^2+B_\varphi^2\right)}{r} \,,\\
    S_{\mathrm{geom}}[m_\theta] &=& \frac{p_{\mathrm{tot}}\cot\theta+\cot\theta\left(\rho v_\varphi^2-B_\varphi^2\right)-\rho v_r v_\theta+B_r B_\theta}{r} \,,\\
    S_{\mathrm{geom}}[m_\varphi] &=& \frac{-\rho v_\varphi\left(v_r+v_\theta\cot\theta\right)+B_\varphi\left(B_r+B_\theta\cot\theta\right)}{r} \,,\\
    S_{\mathrm{geom}}[B_r] &=& \frac{2\psi}{r} \,,\\
    S_{\mathrm{geom}}[B_\theta] &=& \frac{v_r B_\theta-v_\theta B_r+\psi\cot\theta}{r} \,,\\
    S_{\mathrm{geom}}[B_\varphi] &=& \frac{v_r B_\varphi-v_\varphi B_r+\cot\theta\left(v_\theta B_\varphi-v_\varphi B_\theta\right)}{r} \,,
\end{eqnarray}

where ψ is the GLM field (see [Divergence B source
treatments](equations.html#eq_divB_fix)); ψ itself has no geometric source.

### Special relativistic hydrodynamics: srhd {: #geom_spherical_srhd }

Structurally identical to the `hd` case, with the momentum density
\( \mathbf{m} = \xi\mathbf{v} \) (rather than \( \rho\mathbf{v} \)) and
using the spatial four-velocity components \( \mathbf{u}\equiv\Gamma\mathbf{v} \)
(the primitive momentum-slot variables; \( \mathbf{v}=\mathbf{u}/\Gamma \)),
so that \( \xi v_iv_j = (\xi/\Gamma^2)\,u_iu_j \):

\begin{eqnarray}
    S_{\mathrm{geom}}[m_r] &=& \frac{2p+(\xi/\Gamma^2)\left(u_\theta^2+u_\varphi^2\right)}{r} \,,\\
    S_{\mathrm{geom}}[m_\theta] &=& \frac{p\cot\theta+(\xi/\Gamma^2)\left(u_\varphi^2\cot\theta-u_r u_\theta\right)}{r} \,,\\
    S_{\mathrm{geom}}[m_\varphi] &=& -\frac{(\xi/\Gamma^2)\,u_\varphi\left(u_r+u_\theta\cot\theta\right)}{r} \,,
\end{eqnarray}

with \( \xi=\Gamma^2\rho h \) as in [Physics modules and
equations](equations.html#eq_srhd).

## Cylindrical (r, z, φ) {: #geom_cylindrical_sources }

Only the \(m_r\)/\(m_\varphi\) (and, for `mhd`, \(B_r\)/\(B_\varphi\))
components pick up a geometric source; \(m_z\)/\(B_z\) do not, since a
cylindrical volume element's z-extent does not depend on r.

### Hydrodynamics: hd {: #geom_cylindrical_hd }

\begin{eqnarray}
    S_{\mathrm{geom}}[m_r] &=& \frac{p+\rho v_\varphi^2}{r} \,,\\
    S_{\mathrm{geom}}[m_\varphi] &=& -\frac{\rho v_\varphi v_r}{r} \,.
\end{eqnarray}

### Magnetohydrodynamics: mhd {: #geom_cylindrical_mhd }

With \( p_{\mathrm{tot}} = p + \frac{1}{2}B^2 \):

\begin{eqnarray}
    S_{\mathrm{geom}}[m_r] &=& \frac{p_{\mathrm{tot}}-B_\varphi^2+\rho v_\varphi^2}{r} \,,\\
    S_{\mathrm{geom}}[m_\varphi] &=& \frac{-\rho v_\varphi v_r+B_\varphi B_r}{r} \,,\\
    S_{\mathrm{geom}}[B_r] &=& \frac{\psi}{r} \,,\\
    S_{\mathrm{geom}}[B_\varphi] &=& \frac{B_\varphi v_r-B_r v_\varphi}{r} \,.
\end{eqnarray}

### Special relativistic hydrodynamics: srhd {: #geom_cylindrical_srhd }

As for spherical, with \( \mathbf{u}\equiv\Gamma\mathbf{v} \) and
\( \xi=\Gamma^2\rho h \):

\begin{eqnarray}
    S_{\mathrm{geom}}[m_r] &=& \frac{p+(\xi/\Gamma^2)\,u_\varphi^2}{r} \,,\\
    S_{\mathrm{geom}}[m_\varphi] &=& -\frac{(\xi/\Gamma^2)\,u_\varphi u_r}{r} \,.
\end{eqnarray}

# Current limitations {: #geom_limitations }

* **No polar-axis/cylindrical-axis handling.** Keep the domain away from
  θ = 0 and θ = π (spherical), or from r = 0 (cylindrical).
* **AMR across curvilinear refinement levels is untested.**

See `CLAUDE.md` in the repository root for the implementation-level detail
behind both of these (which files are involved, and which tests exercise
what).
