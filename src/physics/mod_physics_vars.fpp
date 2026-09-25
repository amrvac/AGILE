#:mute
#:include "../mod_gpu_directives.fpp"
#:include 'mod_physics_templates.fpp'
#:endmute

module mod_physics_vars
  use mod_global_parameters, only: name_len, ndim

  implicit none
  public
  
  double precision :: phys_gamma=5.d0/3.d0
  ${GPU_DECLARE_COPYIN('phys_gamma')}$

  !> String describing the physics type of the simulation
  character(len=name_len) :: physics_type = "${PHYS}$"
  ${GPU_DECLARE_COPYIN('physics_type')}$

  !> To use wider stencils in flux calculations. A value of 1 will extend it by
  !> one cell in both directions, in any dimension
  integer :: phys_wider_stencil = 0
  ${GPU_DECLARE_COPYIN('phys_wider_stencil')}$

  !> Array per direction per variable, which can be used to specify that certain
  !> fluxes have to be treated differently
  integer, allocatable :: flux_type(:, :)
  ${GPU_DECLARE_CREATE('flux_type')}$

  !> Indicates a normal flux
  integer, parameter   :: flux_default        = 0
  !> Indicates the flux should be treated with tvdlf
  integer, parameter   :: flux_tvdlf          = 1
  
  !> Whether the physics routines require diagonal ghost cells, for example for
  !> computing a curl.
  logical :: phys_req_diagonal = .true.
  ${GPU_DECLARE_COPYIN('phys_req_diagonal')}$

  !> Solve energy equation or not
  logical :: phys_energy=.false.
  ${GPU_DECLARE_COPYIN('phys_energy')}$
  
  !> Solve total energy equation or not
  logical :: phys_total_energy=.false.
  ${GPU_DECLARE_COPYIN('phys_total_energy')}$

  !> Solve internal energy instead of total energy
  logical :: phys_internal_e=.false.
  ${GPU_DECLARE_COPYIN('phys_internal_e')}$

  !> Solve partially ionized one-fluid plasma
  logical :: phys_partial_ionization=.false.
  ${GPU_DECLARE_COPYIN('phys_partial_ionization')}$

  !> if equilibrium pressure is splitted
  logical :: phys_equi_pe=.false.
  ${GPU_DECLARE_COPYIN('phys_equi_pe')}$

  @:phys_vars()

end module mod_physics_vars
