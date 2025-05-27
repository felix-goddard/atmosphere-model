module mod_constants

    use mod_kinds, only: ik, rk

    implicit none

    ! mathematical constants
    real(rk), parameter :: pi = 3.1415926535897932384626433

    ! physical constants
    real(rk), parameter :: R_dry = 287.052874 ! specific gas constant of dry air, J/kg/K
    real(rk), parameter :: cp_dry = 1005 ! constant pressure heat capacity of dry air, J/kg/K
    real(rk), parameter :: cv_dry = cp_dry - R_dry ! constant volume heat capacity of dry air, J/kg/K

    ! model specific constants
    real(rk) :: top_pressure
    real(rk) :: gravity
    real(rk) :: coriolis_parameter

    ! derived constants
    real(rk) :: kappa
    
contains

    subroutine set_constants(ptop, g, f)
        real(rk), intent(in) :: ptop, g, f

        top_pressure = ptop
        gravity = g
        coriolis_parameter = f
    
        kappa = R_dry / cp_dry
        
    end subroutine set_constants
    
end module mod_constants