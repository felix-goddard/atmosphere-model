module mod_radiation

   use mod_kinds, only: ik, rk
   use mod_log, only: logger => main_logger, log_str
   use mod_config, only: config => main_config
   use mod_constants, only: pi, kappa, gravity
   use mod_tiles, only: is, ie, js, je, isd, ied, jsd, jed
   use mod_sync, only: halo_exchange
   use mod_fields, only: dp, pt, ts, play, playkap, plev, pkap, net_flux
   use mod_gas_optics, only: n_g_points, &
                             shortwave_absorption_coefficient, &
                             shortwave_scattering_coefficient, &
                             longwave_absorption_coefficient, &
                             planck_function
   use mod_solar, only: solar_irradiance, irradiance_fraction

   implicit none
   private

   public :: allocate_radiation_arrays, calculate_radiative_heating, &
             radiation_halo_exchange

   integer(ik) :: isl, iel, jsl, jel

   real(rk), allocatable :: upward_longwave_flux(:, :, :)
   real(rk), allocatable :: downward_longwave_flux(:, :, :)
   real(rk), allocatable :: upward_shortwave_flux(:, :, :)
   real(rk), allocatable :: downward_shortwave_flux(:, :, :)

   real(rk), allocatable :: layer_temperature(:, :, :)
   real(rk), allocatable :: level_temperature(:, :, :)

   real(rk), allocatable :: dm(:)
   real(rk), allocatable :: flux_up(:), flux_dn(:)
   real(rk), allocatable :: transmittance(:), reflectance(:)
   real(rk), allocatable :: albedo(:), emission_up(:), beta(:)
   real(rk), allocatable :: source_up(:), source_dn(:), solar_beam(:)

contains

   subroutine calculate_radiative_heating(time)
      real(rk), intent(in) :: time
      integer(ik) :: i, j, k, g
      real(rk), allocatable :: p1(:, :), p2(:, :)
      real(rk) :: irradiance, solar_zenith_angle, &
                  surface_albedo, surface_emissivity, plog(isl:iel, jsl:jel)

      ! ========================================================================
      ! Calculate layer properties of the entire atmosphere
      ! This populates the arrays plev, play, plog, layer_temperature, and
      ! level_temperature

      do k = config%nlay, 1, -1
         plev(isl:iel, jsl:jel, k) = &
            plev(isl:iel, jsl:jel, k + 1) + dp(isl:iel, jsl:jel, k)

         pkap(isl:iel, jsl:jel, k) = plev(isl:iel, jsl:jel, k)**kappa

         plog(isl:iel, jsl:jel) = &
            log(plev(isl:iel, jsl:jel, k)/plev(isl:iel, jsl:jel, k + 1))

         play(isl:iel, jsl:jel, k) = dp(isl:iel, jsl:jel, k)/plog

         playkap(isl:iel, jsl:jel, k) = &
            (pkap(isl:iel, jsl:jel, k + 1) - pkap(isl:iel, jsl:jel, k)) &
            /(-plog*kappa)

         layer_temperature(isl:iel, jsl:jel, k) = &
            pt(isl:iel, jsl:jel, k)*playkap(isl:iel, jsl:jel, k)

         ! linear interpolation of T in log p
         if (k < config%nlay) then
            p1 = log(play(isl:iel, jsl:jel, k))
            p2 = log(play(isl:iel, jsl:jel, k + 1))
            plog = log(plev(isl:iel, jsl:jel, k + 1))

            level_temperature(isl:iel, jsl:jel, k + 1) = &
               (layer_temperature(isl:iel, jsl:jel, k)*(p2 - plog) &
                + layer_temperature(isl:iel, jsl:jel, k + 1) &
                *(plog - p1))/(p2 - p1)
         end if
      end do

      ! linear extrapolation of T in log p for the model top
      p1 = log(play(isl:iel, jsl:jel, config%nlay - 1))
      p2 = log(play(isl:iel, jsl:jel, config%nlay))
      plog = log(plev(isl:iel, jsl:jel, config%nlay + 1))

      level_temperature(isl:iel, jsl:jel, config%nlay + 1) = &
         (layer_temperature(isl:iel, jsl:jel, config%nlay - 1)*(p2 - plog) &
          + layer_temperature(isl:iel, jsl:jel, config%nlay) &
          *(plog - p1))/(p2 - p1)

      ! linear extrapolation of T in log p for the model bottom
      p1 = log(play(isl:iel, jsl:jel, 1))
      p2 = log(play(isl:iel, jsl:jel, 2))
      plog = log(plev(isl:iel, jsl:jel, 1))

      level_temperature(isl:iel, jsl:jel, 1) = &
         (layer_temperature(isl:iel, jsl:jel, 1)*(p2 - plog) &
          + layer_temperature(isl:iel, jsl:jel, 2)*(plog - p1))/(p2 - p1)

      ! ========================================================================
      ! Now solve the two-stream equations; this is done one column at a time
      ! via a loop over the whole horizontal domain

      solar_zenith_angle = 0. ! radians
      surface_albedo = .3
      surface_emissivity = .99

      upward_shortwave_flux(:, :, :) = 0.
      downward_shortwave_flux(:, :, :) = 0.
      upward_longwave_flux(:, :, :) = 0.
      downward_longwave_flux(:, :, :) = 0.
      net_flux(:, :, :) = 0.

      do i = isl, iel
         do j = jsl, jel

            dm(:) = dp(i, j, :)/gravity

            irradiance = solar_irradiance(time, i, j)

            do g = 1, n_g_points

               call solve_shortwave_column( &
                  g, play(i, j, :), dm(:), layer_temperature(i, j, :), &
                  irradiance_fraction(g)*irradiance, solar_zenith_angle, &
                  surface_albedo)

               upward_shortwave_flux(i, j, :) = &
                  upward_shortwave_flux(i, j, :) + flux_up(:)

               downward_shortwave_flux(i, j, :) = &
                  downward_shortwave_flux(i, j, :) + flux_dn(:)

               call solve_longwave_column( &
                  g, play(i, j, :), dm(:), layer_temperature(i, j, :), &
                  level_temperature(i, j, :), ts(i, j), surface_emissivity)

               upward_longwave_flux(i, j, :) = &
                  upward_longwave_flux(i, j, :) + flux_up(:)

               downward_longwave_flux(i, j, :) = &
                  downward_longwave_flux(i, j, :) + flux_dn(:)

            end do

            net_flux(i, j, :) = (downward_shortwave_flux(i, j, :) &
                                 + downward_longwave_flux(i, j, :)) &
                                - (upward_shortwave_flux(i, j, :) &
                                   + upward_longwave_flux(i, j, :))
         end do
      end do

   end subroutine calculate_radiative_heating

   subroutine solve_shortwave_column( &
      g_point, pressure, layer_mass, temperature, &
      irradiance, sza, surface_albedo)

      integer(ik), intent(in) :: g_point
      real(rk), intent(in) :: pressure(:), layer_mass(:), &
                              temperature(:), irradiance, &
                              sza, surface_albedo
      real(rk) :: optical_thickness, single_scattering_albedo, asymmetry, &
                  k_ext, k_sca, k_abs, &
                  backscatter, gamma1, gamma2, gamma3, gamma4, alpha1, alpha2, &
                  tmp1, tmp2, lambda, e_minus, e_2minus, mu0, asy2, lam_mu, &
                  direct_transmittance, diffuse_reflectance, &
                  diffuse_transmittance
      integer(ik) :: k

      mu0 = cos(sza)
      solar_beam(config%nlay + 1) = mu0*irradiance

      do k = config%nlay, 1, -1

         k_abs = shortwave_absorption_coefficient( &
                 g_point, pressure(k), temperature(k))
         k_sca = shortwave_scattering_coefficient(g_point)
         k_ext = k_abs + k_sca

         optical_thickness = k_ext*layer_mass(k)
         single_scattering_albedo = k_sca/k_ext
         asymmetry = 0.

         call delta_eddington_rescaling( &
            optical_thickness, single_scattering_albedo, asymmetry)

         ! Two-stream coefficients for direct to diffuse according to the
         ! hybrid modified Eddington-delta function of Meador & Weaver (1980)
         ! with the Eddington approximation for the backscatter fraction

         asy2 = asymmetry**2

         backscatter = (2.-3.*asymmetry*mu0)/4.
         tmp1 = 4.*(1.-asy2*(1.-mu0))

         gamma3 = single_scattering_albedo*asy2
         gamma4 = 4.*backscatter + 3.*asymmetry
         tmp2 = 3.*asymmetry

         gamma1 = (7.-3.*asy2 - single_scattering_albedo*(4.+tmp2) &
                   + gamma3*gamma4)/tmp1

         gamma2 = (-1.+asy2 + single_scattering_albedo*(4.-tmp2) &
                   + gamma3*(gamma4 - 4.))/tmp1

         gamma3 = backscatter
         gamma4 = 1.-gamma3

         alpha1 = gamma1*gamma4 + gamma2*gamma3
         alpha2 = gamma2*gamma3 + gamma2*gamma4
         lambda = sqrt(gamma1**2 - gamma2**2)

         lam_mu = lambda*mu0
         e_minus = exp(-lambda*optical_thickness)
         e_2minus = e_minus**2
         direct_transmittance = exp(-optical_thickness/mu0)

         tmp1 = single_scattering_albedo &
                /((1.-lam_mu**2) &
                  *(lambda + gamma1 + (lambda - gamma1)*e_2minus))

         ! Direct to diffuse upward
         tmp2 = lambda*gamma3
         diffuse_reflectance = &
            tmp1*( &
            (1.-lam_mu)*(alpha2 + tmp2) - (1.+lam_mu)*(alpha2 - tmp2)*e_2minus &
            - 2.*lambda*(gamma3 - alpha2*mu0)*e_minus*direct_transmittance)

         ! Direct to diffuse downward
         tmp2 = lambda*gamma4
         diffuse_transmittance = &
            tmp1*(2.*lambda*(gamma4 + alpha1*mu0)*e_minus &
                  - direct_transmittance*( &
                  (1.+lam_mu)*(alpha1 + tmp2) &
                  - (1.-lam_mu)*(alpha1 - tmp2)*e_2minus))

         solar_beam(k) = direct_transmittance*solar_beam(k + 1)
         source_up(k) = diffuse_reflectance*solar_beam(k + 1)
         source_dn(k) = diffuse_transmittance*solar_beam(k + 1)

         ! Two-stream coefficients for diffuse to diffuse according to the
         ! hemispheric constant approximation of Meador & Weaver (1980), with
         ! the backscatter coefficient for diffuse radiation of Ritter
         ! & Geleyn (1992)

         backscatter = (4.+asymmetry)/(8.*(1.+asymmetry))
         gamma1 = 2.*(1.-single_scattering_albedo*(1.-backscatter))
         gamma2 = 2.*single_scattering_albedo*backscatter

         lambda = sqrt(gamma1**2 - gamma2**2)
         e_2minus = e_minus**2
         tmp1 = 1./(lambda + gamma1 + (lambda - gamma1)*e_2minus)

         ! Diffuse to diffuse reflectance and transmittance
         reflectance(k) = tmp1*gamma2*(1.-e_2minus)
         transmittance(k) = 2.*tmp1*lambda*e_minus

      end do

      call calculate_fluxes(solar_beam(config%nlay + 1), 0._rk, surface_albedo)

      flux_dn(:) = flux_dn(:) + solar_beam(:)

   end subroutine solve_shortwave_column

   subroutine solve_longwave_column( &
      g_point, pressure, layer_mass, temperature, &
      level_temperature, surface_temperature, surface_emissivity)

      integer(ik), intent(in) :: g_point
      real(rk), intent(in) :: pressure(:), layer_mass(:), &
                              temperature(:), level_temperature(:), &
                              surface_temperature, surface_emissivity
      real(rk) :: optical_thickness, single_scattering_albedo, asymmetry, &
                  k_ext, k_sca, k_abs, &
                  backscatter, gamma1, gamma2, lambda, e_minus, e_2minus, &
                  tmp1, surface_planckian, planckian_difference
      real(rk) :: level_planckian(config%nlay + 1)
      integer(ik) :: k

      surface_planckian = planck_function(g_point, surface_temperature)
      level_planckian(1) = planck_function(g_point, level_temperature(1))

      do k = 1, config%nlay

         k_abs = longwave_absorption_coefficient( &
                 g_point, pressure(k), temperature(k))
         k_sca = 0.
         k_ext = k_abs + k_sca

         optical_thickness = k_ext*layer_mass(k)
         single_scattering_albedo = k_sca/k_ext
         asymmetry = 0.

         call delta_eddington_rescaling( &
            optical_thickness, single_scattering_albedo, asymmetry)

         ! Two-stream coefficients according to the hemispheric constant
         ! approximation of Meador & Weaver (1980), with the backscatter
         ! coefficient for diffuse radiation of Ritter & Geleyn (1992)
         ! (corrected for use with the delta-Eddington rescaling)

         backscatter = (4.-3.*asymmetry)/8.
         gamma1 = 2.*(1.-single_scattering_albedo*(1.-backscatter))
         gamma2 = 2.*single_scattering_albedo*backscatter

         lambda = sqrt(gamma1**2 - gamma2**2)
         e_minus = exp(-lambda*optical_thickness)
         e_2minus = e_minus**2
         tmp1 = 1./(lambda + gamma1 + (lambda - gamma1)*e_2minus)

         reflectance(k) = tmp1*gamma2*(1.-e_2minus)
         transmittance(k) = 2.*tmp1*lambda*e_minus

         ! Thermal source function (Planck function linear in optical depth)

         level_planckian(k + 1) = &
            planck_function(g_point, level_temperature(k + 1))

         planckian_difference = &
            (1.+reflectance(k) - transmittance(k)) &
            *(level_planckian(k) - level_planckian(k + 1)) &
            /(optical_thickness*(gamma1 + gamma2))

         tmp1 = 1.-reflectance(k)

         source_up(k) = &
            pi*(tmp1*level_planckian(k + 1) &
                - transmittance(k)*level_planckian(k) &
                + planckian_difference)

         source_dn(k) = &
            pi*(tmp1*level_planckian(k) &
                - transmittance(k)*level_planckian(k + 1) &
                - planckian_difference)

      end do

      call calculate_fluxes(0._rk, &
                            surface_emissivity*pi*surface_planckian, &
                            1.-surface_emissivity)

   end subroutine solve_longwave_column

   subroutine calculate_fluxes(toa_downward_flux, sfc_upward_flux, sfc_albedo)

      ! Calculate fluxes given the transmittance, reflectance, and source
      ! functions of the atmosphere, following the method of Hogan & Bozzo (2018)

      real(rk) :: toa_downward_flux, sfc_albedo, sfc_upward_flux
      integer(ik) :: k, nlay

      nlay = config%nlay

      albedo(1) = sfc_albedo
      emission_up(1) = sfc_upward_flux

      do k = 1, nlay
         beta(k) = 1./(1.-reflectance(k)*albedo(k))
         albedo(k + 1) = reflectance(k) &
                         + (transmittance(k)**2)*beta(k)*albedo(k)
         emission_up(k + 1) = source_up(k) + transmittance(k)*beta(k)*( &
                              emission_up(k) + source_dn(k)*albedo(k))
      end do

      flux_dn(nlay + 1) = toa_downward_flux
      flux_up(nlay + 1) = albedo(nlay + 1)*toa_downward_flux &
                          + emission_up(nlay + 1)

      do k = nlay + 1, 2, -1
         flux_dn(k - 1) = beta(k - 1)*(source_dn(k - 1) &
                                       + transmittance(k - 1)*flux_dn(k) &
                                       + reflectance(k - 1)*emission_up(k - 1))
         flux_up(k - 1) = albedo(k - 1)*flux_dn(k - 1) + emission_up(k - 1)
      end do

   end subroutine calculate_fluxes

   subroutine delta_eddington_rescaling( &
      optical_thickness, single_scattering_albedo, asymmetry)

      real(rk), intent(inout) :: optical_thickness, single_scattering_albedo, &
                                 asymmetry
      real(rk) :: forward_scatter, coefficient

      forward_scatter = asymmetry**2
      coefficient = 1.-single_scattering_albedo*forward_scatter

      optical_thickness = optical_thickness*coefficient
      single_scattering_albedo = &
         single_scattering_albedo*(1.-forward_scatter)/coefficient
      asymmetry = asymmetry/(1.+asymmetry)

   end subroutine delta_eddington_rescaling

   subroutine radiation_halo_exchange()

      call halo_exchange(net_flux)

   end subroutine radiation_halo_exchange

   subroutine allocate_radiation_arrays()
      integer(ik) :: nlay

      nlay = config%nlay

      ! We want to restrict the size of the domain we operate on as much as
      ! possible since the radiation calculations are slow; this is the smallest
      ! domain we have to work on, defined by the domain needed for the
      ! potential temperature heating rate in the D grid shallow water update
      isl = is + 4
      iel = ie - 4
      jsl = js + 4
      jel = je - 4

      if (.not. allocated(upward_longwave_flux)) &
         allocate (upward_longwave_flux(isl:iel, jsl:jel, nlay + 1))
      if (.not. allocated(downward_longwave_flux)) &
         allocate (downward_longwave_flux(isl:iel, jsl:jel, nlay + 1))
      if (.not. allocated(upward_shortwave_flux)) &
         allocate (upward_shortwave_flux(isl:iel, jsl:jel, nlay + 1))
      if (.not. allocated(downward_shortwave_flux)) &
         allocate (downward_shortwave_flux(isl:iel, jsl:jel, nlay + 1))

      if (.not. allocated(layer_temperature)) &
         allocate (layer_temperature(isl:iel, jsl:jel, nlay))
      if (.not. allocated(level_temperature)) &
         allocate (level_temperature(isl:iel, jsl:jel, nlay + 1))

      if (.not. allocated(dm)) allocate (dm(nlay))
      if (.not. allocated(flux_up)) allocate (flux_up(nlay + 1))
      if (.not. allocated(flux_dn)) allocate (flux_dn(nlay + 1))
      if (.not. allocated(transmittance)) allocate (transmittance(nlay))
      if (.not. allocated(reflectance)) allocate (reflectance(nlay))
      if (.not. allocated(albedo)) allocate (albedo(nlay + 1))
      if (.not. allocated(emission_up)) allocate (emission_up(nlay + 1))
      if (.not. allocated(beta)) allocate (beta(nlay))
      if (.not. allocated(source_up)) allocate (source_up(nlay + 1))
      if (.not. allocated(source_dn)) allocate (source_dn(nlay + 1))
      if (.not. allocated(solar_beam)) allocate (solar_beam(nlay + 1))

   end subroutine allocate_radiation_arrays

end module mod_radiation
