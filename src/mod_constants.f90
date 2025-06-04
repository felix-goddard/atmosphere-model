module mod_constants

    use mod_kinds, only: ik, rk

    implicit none

    ! mathematical constants
    real(rk), parameter :: pi = 3.1415926535897932384626433

    ! physical constants
    real(rk), parameter :: R_dry = 287.052874 ! specific gas constant of dry air, J/kg/K
    real(rk), parameter :: cp_dry = 1005 ! constant pressure heat capacity of dry air, J/kg/K
    real(rk), parameter :: cv_dry = cp_dry - R_dry ! constant volume heat capacity of dry air, J/kg/K
    real(rk), parameter :: kappa = R_dry / cp_dry

    ! model specific constants
    real(rk) :: gravity
    real(rk) :: coriolis_parameter
    real(rk) :: divergence_damping_coeff
    real(rk) :: vorticity_damping_coeff
    
contains

    subroutine set_constants(g, f, divergence_damping, vorticity_damping)
        real(rk), intent(in) :: g, f, divergence_damping, vorticity_damping

        gravity = g
        coriolis_parameter = f
        divergence_damping_coeff = divergence_damping
        vorticity_damping_coeff = vorticity_damping
        
    end subroutine set_constants
    
end module mod_constants