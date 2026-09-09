title: Physics modules and equations

# List of physics modules {: #eq_list_physics }
This document describes the equations implemented.
Information about user defined source terms are in the user module. In
principle, the code handles anything of generic form

\begin{equation}
    \partial_t \mathbf{U} + \nabla \cdot \mathbf{F}(\mathbf{U}) = \mathbf{S}(\mathbf{U},\nabla \mathbf{U},\mathbf{x},t)
\end{equation}

The code is configured to use the specified set of equations by activating it via 
the code's parameter file in the methodlist

```fortran
&methodlist
   phys            = 'mhd'
/
```

where `phys` is one of the implemented physics modules (`hd`, `mhd`, `ffhd`, `srhd`),
see below. The equations and notation below follow
[Porth et al., "Astrophysics on GPUs: introducing AGILE 1.0"
(submitted to RASTI, in review; arXiv:2607.19277)](https://arxiv.org/abs/2607.19277),
Section 2.3, which is also the citable reference for AGILE 1.0.

| Module | Available solvers |
| --- | --- |
| `hd`   | `LF`, `HLL`, `HLLC` |
| `ffhd` | `LF` |
| `mhd`  | `LF`, `HLL` |
| `srhd` | `LF` |

## Hydrodynamics: hd {: #eq_hd }

```fortran
&methodlist
   phys            = 'hd'
/
```

The hydrodynamic module solves the partial differential equations that
govern compressible gas dynamics. The conservative variables are density,
vectorial momentum density, hydrodynamic energy density and a
user-selected set of `n_tr` tracers, collected in
**U** = [ρ, **m**≡ρ**v**, e_hd, D_tr^i], obeying

\begin{eqnarray}
 \partial_t \rho + \nabla \cdot (\rho\mathbf{v}) & =  & 0  \,,\\
  \partial_t \mathbf{m} + \nabla \cdot \left(\mathbf{v}\mathbf{m}+p\mathbf{I}\right) & = & \rho\mathbf{g} \,,\\
  \partial_t e_{\mathrm{hd}}+\nabla\cdot\left[\mathbf{v}(e_{\mathrm{hd}}+p)\right]
     & = &
     \rho\mathbf{v}\cdot\mathbf{g} -n_e  n_H\Lambda(T) \,,\\
  \partial_t D_{tr}^i+\nabla \cdot(\mathbf{v} D_{tr}^i) & = & 0 \,.
 \end{eqnarray}

Here *p* is the gas pressure, **v** the velocity vector and **I** the
identity tensor. Tracers are handled using their conservative variable
D_tr=ρϑ_tr, such that the actual tracer value ϑ_tr obeys a pure advection
equation. The pre-implemented sources on the right-hand-side allow the
user to set the external gravitational acceleration **g** (see
[Gravity](#eq_gravity) below) and select the treatment of local and
instantaneous optically thin radiative losses (see
[Radiative losses](#eq_radloss) below).

The pressure is a derived quantity, obtained from the energy closure

\begin{equation}
    e_{\mathrm{hd}}=\frac{p}{\hat{\gamma}-1}+\frac{1}{2}\rho v^2 \,,
\end{equation}

which introduces the equation parameter **hd_gamma** representing the
constant ratio of specific heats γ̂ (typical value 5/3 for a mono-atomic
gas). An ideal gas law relates *p* = 𝓡ρT, with a fixed gas constant
𝓡 = k_B/(μ̃m_p), where k_B is the Boltzmann constant, m_p the proton mass,
and the constant mean molecular weight μ̃ = (1+4He_a)/(2+3He_a) is set by
the fixed Helium abundance **He_abundance** (He_a ≡ n_He/n_H).

Parameters of hydrodynamics are read in the **hd_list** of the parameter
file, including the equation parameter **hd_gamma**.

This equation module can be combined with physical sources for
(local) optically thin radiative losses by setting **hd_radiative_cooling=.true.**,
see the [radiative cooling](radiative_cooling.html) page. It can also be
combined with the external gravity module by setting **hd_gravity=.true.**,
see [Gravity](#eq_gravity) below.

### Radiative losses {: #eq_radloss }
When **hd_radiative_cooling**/**mhd_radiative_cooling**/**ffhd_radiative_cooling**
is enabled, the energy equation of the corresponding module is augmented
with the optically thin radiative loss term

\begin{equation}
    -n_e n_H \Lambda(T) \,,
\end{equation}

where Λ(T) is a tabulated radiative loss function depending only on the
local temperature T (see [Radiative cooling](radiative_cooling.html) for
the available curves and the exact integration scheme of Townsend 2009
used to advance this term). These losses scale with the squared density:
ρ is converted to electron and Hydrogen number densities under the
assumption of constant, fully ionized conditions with a fixed Helium
abundance He_a ≡ n_He/n_H, so that charge neutrality sets

\begin{equation}
    n_e=n_H(1+2\,\mathrm{He}_a) \,,\qquad
    \rho=n_H m_p(1+4\,\mathrm{He}_a) \,.
\end{equation}

### Gravity {: #eq_gravity }
The external gravitational acceleration **g** enters the momentum and
energy sources (ρ**g** and ρ**v**·**g** respectively) of the `hd`, `mhd`
and `ffhd` modules alike, activated by setting **hd_gravity**/**mhd_gravity**/**ffhd_gravity**
`= .true.` in the corresponding `_list` of the parameter file. AGILE does
not ship a built-in gravity profile (e.g. no separate "uniform" vs "point
mass" switches); instead, the user supplies **g** directly for the
current cell state by implementing a `gravity_field` function in
`mod_usr.fpp`, for example a constant vertical gravity as used in the
`hd/RT3D` test case:

```fortran
pure real(dp) function gravity_field(wCT, x, idim) result(field)
  !$acc routine seq
  real(dp), intent(in)       :: wCT(nw_phys)
  real(dp), intent(in)       :: x(1:ndim)
  integer, value, intent(in) :: idim

  if (idim == 1) field =  0.0_dp
  if (idim == 2) field =  0.0_dp
  if (idim == 3) field = -1.0_dp
end function gravity_field
```

`idim` selects the spatial direction (1, 2 or 3) for which the
acceleration component is requested; a spatially- or time-varying field
(e.g. a point-mass potential) can be implemented by using `x` (and any
module-level state the user maintains) inside this function.

## Magnetohydrodynamics: mhd {: #eq_mhd }

```fortran
&methodlist
   phys            = 'mhd'
/
```

The MHD module combines mass conservation and the tracer equation above
with the following generalizations of momentum and energy, and adds the
induction equation for the magnetic field **B**:

\begin{equation}
\partial_t \mathbf{m} + \nabla \cdot \left(\mathbf{v}\mathbf{m}+(p+{\textstyle{\frac{1}{2}}}B^2)\mathbf{I}-\mathbf{B}\mathbf{B}\right)  = \rho\mathbf{g}
\end{equation}

\begin{equation}
  \partial_t e_{\mathrm{mhd}}+\nabla\cdot\left[\mathbf{v}(e_{\mathrm{mhd}}+p+{\textstyle{\frac{1}{2}}}B^2)-\mathbf{B}(\mathbf{v}\cdot\mathbf{B})\right]
  = \rho\mathbf{v}\cdot\mathbf{g} -n_e  n_H\Lambda(T) + \nabla \cdot \left[ \mathbf{B} \times \eta \mathbf{J}\right]\,,
\end{equation}

\begin{equation}
\partial_t\mathbf{B}+\nabla\cdot\left[\mathbf{v}\mathbf{B}-\mathbf{B}\mathbf{v} +\psi\mathbf{I}\right] = -\nabla\times\eta \mathbf{J} \,,
\end{equation}

where the total energy density closure is

\begin{equation}
    e_{\mathrm{mhd}}=\frac{p}{\hat{\gamma}-1}+\frac{1}{2}\rho v^2+\frac{1}{2}B^2 \,,
\end{equation}

**J** = ∇×**B** is the current density, and the magnetic field is
measured in units for which the magnetic permeability is 1. For
non-zero resistivity η, the resistive terms are added as "non-local"
source terms using the compact stencil detailed in Porth & Xia et al.
(2014). The GLM field ψ used to control the divergence of **B**
numerically is described in [Divergence B source
treatments](#eq_divB_fix) below.

Parameters of magnetohydrodynamics are read in the **mhd_list** of parameter file.
The source terms on the right hand side with **eta** in them are the resistive
terms.

There are two equation parameters: the polytropic index **mhd_gamma**
(which must be larger or equal to 1), and the resistivity **mhd_eta**. 

This equation module can be combined with physical sources for
(local) optically thin [radiative losses](radiative_cooling.html) by set **mhd_radiative_cooling=.true.**. 
It can also be combined with the external gravity modules by set **mhd_gravity=.true.**.

### Divergence B source treatments {: #eq_divB_fix }

To control magnetic monopole buildup numerically, AGILE uses the
generalized Lagrangian multiplier (GLM) approach of Dedner et al. (2002),
where ψ is an added evolved variable obeying

\begin{equation}
    \partial_t \psi+c_{h}^{2} \nabla \cdot \mathbf{B} = -\frac{c_h^2}{c_p^2}\psi \,.
\end{equation}

This equation has a hyperbolic advection speed c_h and a parabolic
factor c_p. In practice, the left-hand side is dealt with using the flux
prescription c²_max**B**, while the right-hand side (the damping term) is
solved exactly over the time step Δt_n through

\begin{equation}
    \psi^{n+1}=\psi^n\exp\left(-\Delta t_n c_{\mathrm{max}}\alpha_{\mathrm{GLM}}/\min_{i=1,3}(\Delta x_i)\right)
\end{equation}

using second-order operator splitting. This implies availability of the
maximal, global physical propagation speed c_max that also sets the time
step, while

\begin{equation}
    \alpha_{\mathrm{GLM}}=\min_{i=1,3}(\Delta x_i)\,c_h/c_p^2
\end{equation}

becomes a fixed parameter set to 0.5 by default, but generally taken
within [0,1] (Mignone, Tzeferacos & Bodo 2010).

## Frozen-field hydrodynamics: ffhd {: #eq_ffhd }

```fortran
&methodlist
   phys            = 'ffhd'
/
```

As an intermediate towards full 3D MHD applications, the frozen-field
hydrodynamic (FFHD) module introduced by Zhou, Li & Keppens (2024, 2025)
solves hydrodynamics along a given *frozen* 3D magnetic topology
**B**(**x**) that does not evolve in time. Introducing the local
magnetic unit vector b̂ = **B**/B, this module evolves the set
**U** = [ρ, m_∥=ρv_∥, e_hd∥, q_∥], the magnetic-field-line-projected
variant of the hydrodynamic mass, momentum and energy equations above,
augmented with a field-aligned heat flux q_∥:

\begin{eqnarray}
    \partial_t \rho+\nabla \cdot  \left (  \rho v_{\parallel} \mathbf{\hat{b}}\right ) &=& 0\,,\\
    \partial_t m_{\parallel}+\nabla \cdot \left [  \left (  \rho v^2_{\parallel} +p\right ){\mathbf{\hat{b}}}\right ] &=& \rho \mathbf{g}\cdot\mathbf{\hat{b}} + p (\nabla\cdot \mathbf{\hat{b}})\,,\\
    \partial_t e_{hd\parallel} + \nabla \cdot \left[\left((e_{hd\parallel} + p) v_{\parallel}-q_\parallel\right) \mathbf{\hat{b}}\right] &=& \rho v_{\parallel} \mathbf{g}\cdot\mathbf{\hat{b}}  -n_e  n_H\Lambda(T) + H_{\mathrm{user}}(\mathbf{x})\,,\\
    \partial_t q_\parallel & = & -\frac{\left(q_\parallel-\kappa_\parallel\mathbf{\hat{b}}\cdot \nabla T\right)}{\tau} \,,
\end{eqnarray}

where now e_hd∥ = p/(γ̂-1) + ρv_∥²/2. This uses a hyperbolic approach to
handle anisotropic (purely field-aligned) thermal conduction, where
κ_∥ = κ_∥,0 T^(5/2), with κ_∥,0 the equation parameter **hypertc_kappa**.
This prescription forces the parallel heat flux q_∥ to approach its
target expression κ_∥ b̂·∇T within at most four CFL-limited timesteps
(hence setting τ = max(4Δt_n, τ_reduced), where various prescriptions for
τ_reduced are possible; Rempel 2017, Warnecke & Bingert 2020, Navarro et
al. 2022, Zhou et al. 2025). A purely local user-defined heating
prescription can be coded up by the user in H_user(**x**) (see
`usr_source`/`usr_source_usr` in `mod_usr.fpp`).

Note that FFHD re-uses the external gravity module from the hydro case
(see [Gravity](#eq_gravity) above), but here only the effective gravity
projected along the fixed field, **g**·b̂, enters. The radiative loss
treatment (see [Radiative losses](#eq_radloss) above) is shared with all
other physics modules, activated with **ffhd_radiative_cooling=.true.**.
The equation parameter **ffhd_gamma** sets γ̂, analogous to **hd_gamma**.

## Special relativistic hydrodynamics: srhd {: #eq_srhd }

```fortran
&methodlist
   phys            = 'srhd'
/
```

The special relativistic hydrodynamics module solves the following set
of conservation laws in flat (Minkowski) spacetime

\begin{eqnarray}
    \partial_t D + \nabla\cdot (D \mathbf{v}) &=& 0 \,,\\
    \partial_t \mathbf{m} + \nabla \cdot (\mathbf{v}\mathbf{m} + p\mathbf{I}) &=& 0 \,,\\
    \partial_t \tau + \nabla \cdot\left(\mathbf{m} - D \mathbf{v}\right) &=& 0 \,,
\end{eqnarray}

where the conserved quantities are the density D, the momentum-density
**m** and the rest-mass-subtracted total energy density τ:=e-D, all
evaluated in the lab-frame. The conserved variables are

\begin{equation}
    \mathbf{U} = \left(\begin{array}{c}
        D  \\
        \mathbf{m} \\
        \tau
    \end{array}
    \right)
    = \left(\begin{array}{c}
        \Gamma \rho \\
        \Gamma^2 \rho h \mathbf{v} \\
        \Gamma^2 \rho h - p - D
    \end{array}
    \right) \,,
\end{equation}

with h being the specific enthalpy of the gas and Γ the Lorentz factor
Γ = 1/√(1-v²). The primitive variables are

\begin{equation}
    \mathbf{P} = \left(\begin{array}{c}
        \rho  \\
        \Gamma \mathbf{v} \\
        p
    \end{array}
    \right) \,,
\end{equation}

which correspond to the fluid-frame density ρ, the spatial components of
the four-velocity Γ**v** and the fluid-frame pressure p. Since primitive
variables are required to evaluate the fluxes, a non-linear inversion
**P**(**U**) is needed; AGILE uses the one introduced by Bergmans,
Keppens, van Odyck & Achterberg (2005), which expresses the energy
equation as a 1D root-finding problem in ξ = Γ²ρh. AGILE also saves the
set of auxiliary variables

\begin{equation}
    \mathbf{A} = \left(\begin{array}{c}
        \xi  \\
        \Gamma
    \end{array}
    \right) \,,
\end{equation}

useful for quick conversion between **P** and **U** and to initialize an
initial guess for the root-finder.

As with all Euler-type equations, this set needs to be closed with an
equation of state (EOS) of the form h=h(ρ,p). The two options currently
implemented (selected via the **srhd_eos** logical parameter) are the
*ideal gas* (**srhd_eos=.false.**, the default)

\begin{equation}
    h(\rho,p) = 1 + \frac{\hat{\gamma}}{\hat{\gamma}-1}\frac{p}{\rho} \,,
\end{equation}

and the *Synge gas* (**srhd_eos=.true.**) of the single-species
trans-relativistic gas (Synge 1957), defined as

\begin{equation}
    h(\rho,p) = \frac{K_3(\theta^{-1})}{K_2(\theta^{-1})} \,,
\end{equation}

with the relativistic temperature θ := p/ρ and K_n the modified Bessel
function of the second kind. As in previous work (Meliani & Sauty 2004;
Mignone & Plewa 2005; Ryu, Chattopadhyay & Choi 2006; Porth & Olivares et
al. 2017), AGILE uses an approximation to the previous expression that
avoids the expensive evaluation of the Bessel functions (Mathews 1971;
Goedbloed, Keppens & Poedts 2010). The equation parameter **srhd_gamma**
sets γ̂ for the ideal-gas EOS (and the adiabatic index used in the Synge
gas approximation).


# Positivity fixes {: #eq_positivity_fixes }

TODO
