module mod_solar

   use mod_kinds, only: ik, rk
   use mod_fields, only: radius
   use mod_gas_optics, only: n_g_points

   implicit none
   private

   public :: solar_irradiance, irradiance_fraction

   ! This array gives the fraction of the total solar irradiance represented
   ! by each point in the integration of the k-distribution; this depends on
   ! the solar spectrum and can be calculated using the data in the CKD tables.
   ! Currently we assume the solar spectrum is a blackbody spectrum at 5777 K.
   real(rk), parameter :: irradiance_fraction(n_g_points) = &
                          [5.432097281e-02, 2.913444484e-02, 2.150502170e-01, &
                           1.194550602e-01, 4.343948153e-02, 2.444170963e-02, &
                           9.025298172e-03, 1.325478654e-02, 2.346586433e-03, &
                           2.243590817e-03, 1.116597039e-01, 1.287852342e-01, &
                           1.259809419e-01, 8.776475907e-02, 1.398279548e-02, &
                           1.911441813e-02]

contains

   function solar_irradiance(time, i, j) result(irradiance)
      real(rk), intent(in) :: time
      integer(ik), intent(in) :: i, j
      real(rk) :: irradiance

      !   irradiance = 1370.

      irradiance = 600*(1.+.2/(1.+(radius(i, j)/100000e3)**2))

   end function solar_irradiance

end module mod_solar
