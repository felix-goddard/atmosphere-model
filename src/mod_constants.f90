module mod_constants

    use mod_kinds, only: ik, rk

    implicit none

    real(rk) :: reference_pressure
    real(rk) :: top_pressure
    real(rk) :: gravity
    real(rk) :: coriolis_parameter
    real(rk) :: dry_gas_constant
    real(rk) :: dry_heat_capacity

    ! derived constants
    real(rk) :: kappa
    real(rk) :: hydrostatic_constant
    
contains

    subroutine set_constants(pref, ptop, g, f, R, cP)
        real(rk), intent(in) :: pref, ptop, g, f, R, cP

        reference_pressure = pref
        top_pressure = ptop
        gravity = g
        coriolis_parameter = f
        dry_gas_constant = R
        dry_heat_capacity = cP
    
        kappa = dry_gas_constant / dry_heat_capacity
        hydrostatic_constant = dry_heat_capacity / reference_pressure**kappa
        
    end subroutine set_constants
    
end module mod_constants