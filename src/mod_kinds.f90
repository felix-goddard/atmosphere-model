module mod_kinds

    use iso_fortran_env, only: int32, real64
  
    implicit none
  
    private
    public :: ik, rk
  
    integer, parameter :: rk = real64
    integer, parameter :: ik = int32
  
  end module mod_kinds