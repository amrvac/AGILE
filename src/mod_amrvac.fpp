!> Compatibility alias for mod_agile, so that a mod_usr.fpp ported from
!> upstream MPI-AMRVAC compiles with its `use mod_amrvac` unchanged. New cases
!> should `use mod_agile`.
module mod_amrvac
  use mod_agile

  implicit none
  public

end module mod_amrvac
