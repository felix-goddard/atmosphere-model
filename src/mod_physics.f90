module mod_physics

   use mod_kinds, only: ik, rk
   use mod_config, only: config => main_config
   use mod_constants, only: gravity, kappa, cp_dry, cv_dry
   use mod_tiles, only: is, ie, js, je, isd, ied, jsd, jed
   use mod_sync, only: halo_exchange
   use mod_fields, only: dp, pt, ud, vd, ts, gz, plev, pkap, playkap, net_flux

   implicit none
   private

   public :: allocate_physics_arrays, physics_halo_exchange, &
             calculate_physics_tendencies, apply_physics_tendencies

   real(rk), allocatable :: dp_tend(:, :, :) ! pressure thickness tendency
   real(rk), allocatable :: pt_tend(:, :, :) ! potential temperature tendency
   real(rk), allocatable :: u_tend(:, :, :) ! u wind tendency (on A grid)
   real(rk), allocatable :: v_tend(:, :, :) ! u wind tendency (on A grid)
   real(rk), allocatable :: ts_tend(:, :) ! surface temperature tendency

   real(rk), allocatable :: ua(:, :, :) ! u wind on A grid
   real(rk), allocatable :: va(:, :, :) ! v wind on A grid

contains

   subroutine calculate_physics_tendencies()
      integer(ik) :: i, j, k

      dp_tend(:, :, :) = 0.
      pt_tend(:, :, :) = 0.
      u_tend(:, :, :) = 0.
      v_tend(:, :, :) = 0.

      ! Physics tendencies are calculated on the A-grid; for the winds, the
      ! tendencies themselves are interpolated back to the D-grid.

      do concurrent(i=is + 1:ie - 1, j=js + 1:je - 1)
         ua(i, j, :) = (.25*(ud(i, j, :) + ud(i, j + 1, :)) &
                        + .125*(ud(i - 1, j, :) + ud(i - 1, j + 1, :) &
                                + ud(i + 1, j, :) + ud(i + 1, j + 1, :)))

         va(i, j, :) = (.25*(vd(i, j, :) + vd(i + 1, j, :)) &
                        + .125*(vd(i, j - 1, :) + vd(i + 1, j - 1, :) &
                                + vd(i, j + 1, :) + vd(i + 1, j + 1, :)))
      end do

      do k = config%nlay, 1, -1
         plev(is + 1:ie - 1, js + 1:je - 1, k) = &
            plev(is + 1:ie - 1, js + 1:je - 1, k + 1) &
            + dp(is + 1:ie - 1, js + 1:je - 1, k)

         pkap(is + 1:ie - 1, js + 1:je - 1, k) = &
            plev(is + 1:ie - 1, js + 1:je - 1, k)**kappa
      end do

      do k = 1, config%nlay
         gz(is + 1:ie - 1, js + 1:je - 1, k + 1) = &
            gz(is + 1:ie - 1, js + 1:je - 1, k) &
            + cp_dry*pt(is + 1:ie - 1, js + 1:je - 1, k)*( &
            pkap(is + 1:ie - 1, js + 1:je - 1, k) &
            - pkap(is + 1:ie - 1, js + 1:je - 1, k + 1))

         playkap(is + 1:ie - 1, js + 1:je - 1, k) = &
            (pkap(is + 1:ie - 1, js + 1:je - 1, k + 1) &
             - pkap(is + 1:ie - 1, js + 1:je - 1, k)) &
            /log(plev(is + 1:ie - 1, js + 1:je - 1, k + 1) &
                 /plev(is + 1:ie - 1, js + 1:je - 1, k))/kappa
      end do

      call dry_convective_adjustment()

      call surface_heating()

   end subroutine calculate_physics_tendencies

   subroutine dry_convective_adjustment()
      integer(ik) :: i, j, k
      real(rk) :: ri, exchange, mixing_fraction, &
                  energy(config%nlay), u_adj(config%nlay), v_adj(config%nlay), &
                  gz_adj(config%nlay + 1), pt_adj(config%nlay)
      real(rk), parameter :: min_v2 = 1e-4
      real(rk), parameter :: critical_richardson_number = .25
      real(rk), parameter :: mixing_timescale = 600 ! 10 minutes

      do concurrent(i=isd:ied, j=jsd:jed)

         ! we mix the energy instead of the potential temperature to ensure
         ! we conserve energy

         energy(:) = cv_dry*pt(i, j, :)*playkap(i, j, :) &
                     + .5*(gz(i, j, 2:) + gz(i, j, :config%nlay) &
                           + ua(i, j, :)**2 + va(i, j, :)**2)
         u_adj(:) = ua(i, j, :)
         v_adj(:) = va(i, j, :)

         do k = 1, config%nlay - 1

            ! local richardson number on the interface between layers k and k+1
            ri = (pt(i, j, k + 1) - pt(i, j, k)) &
                 *(dp(i, j, k)*(gz(i, j, k + 1) - gz(i, j, k)) &
                   + dp(i, j, k + 1)*(gz(i, j, k + 2) - gz(i, j, k + 1))) &
                 /((dp(i, j, k)*pt(i, j, k) + dp(i, j, k + 1)*pt(i, j, k + 1)) &
                   *((ua(i, j, k + 1) - ua(i, j, k))**2 &
                     + (va(i, j, k + 1) - va(i, j, k))**2 + min_v2))

            if (ri < critical_richardson_number) then

               mixing_fraction = &
                  (1.-max(0.0, ri/critical_richardson_number))**2 &
                  *dp(i, j, k)*dp(i, j, k + 1)/(dp(i, j, k) + dp(i, j, k + 1))

               exchange = mixing_fraction*(u_adj(k) - u_adj(k + 1))
               u_adj(k) = u_adj(k) - exchange/dp(i, j, k)
               u_adj(k + 1) = u_adj(k + 1) + exchange/dp(i, j, k + 1)

               exchange = mixing_fraction*(v_adj(k) - v_adj(k + 1))
               v_adj(k) = v_adj(k) - exchange/dp(i, j, k)
               v_adj(k + 1) = v_adj(k + 1) + exchange/dp(i, j, k + 1)

               exchange = mixing_fraction*(energy(k) - energy(k + 1))
               energy(k) = energy(k) - exchange/dp(i, j, k)
               energy(k + 1) = energy(k + 1) + exchange/dp(i, j, k + 1)

            end if

         end do

         u_tend(i, j, :) = &
            u_tend(i, j, :) + (u_adj(:) - ua(i, j, k))/mixing_timescale
         v_tend(i, j, :) = &
            v_tend(i, j, :) + (v_adj(:) - va(i, j, k))/mixing_timescale

         ! reconstruct the potential temperature from the new energy

         gz_adj(1) = gz(i, j, 1)
         do k = 1, config%nlay

            pt_adj(k) = (energy(k) - gz_adj(k) &
                         - .5*(u_adj(k)**2 + v_adj(k)**2)) &
                        /(cv_dry*playkap(i, j, k) &
                          + .5*cp_dry*(pkap(i, j, k) - pkap(i, j, k + 1)))

            gz_adj(k + 1) = gz_adj(k) + cp_dry*pt_adj(k) &
                            *(pkap(i, j, k) - pkap(i, j, k + 1))
         end do

         pt_tend(i, j, :) = pt_tend(i, j, :) + &
                            (pt_adj(:) - pt(i, j, :))/mixing_timescale
      end do

   end subroutine dry_convective_adjustment

   subroutine surface_heating()
      real(rk), parameter :: surface_heat_capacity = 41680*5 ! heat capacity of 5 meters of water, J/m^2/K

      ts_tend(isd:ied, jsd:jed) = &
         net_flux(isd:ied, jsd:jed, 1)/surface_heat_capacity

   end subroutine surface_heating

   subroutine apply_physics_tendencies(dt)
      real(rk), intent(in) :: dt

      dp(isd:ied, jsd:jed, :) = dp(isd:ied, jsd:jed, :) + dt*dp_tend(:, :, :)
      pt(isd:ied, jsd:jed, :) = pt(isd:ied, jsd:jed, :) + dt*pt_tend(:, :, :)

      ts(isd:ied, jsd:jed) = ts(isd:ied, jsd:jed) + dt*ts_tend(:, :)

   end subroutine apply_physics_tendencies

   subroutine physics_halo_exchange()

      call halo_exchange(ts)

   end subroutine physics_halo_exchange

   subroutine allocate_physics_arrays()

      if (.not. allocated(dp_tend)) &
         allocate (dp_tend(isd:ied, jsd:jed, config%nlay))
      if (.not. allocated(pt_tend)) &
         allocate (pt_tend(isd:ied, jsd:jed, config%nlay))
      if (.not. allocated(u_tend)) &
         allocate (u_tend(isd:ied, jsd:jed, config%nlay))
      if (.not. allocated(v_tend)) &
         allocate (v_tend(isd:ied, jsd:jed, config%nlay))
      if (.not. allocated(ts_tend)) allocate (ts_tend(isd:ied, jsd:jed))

      if (.not. allocated(ua)) allocate (ua(is:ie, js:je, config%nlay))
      if (.not. allocated(va)) allocate (va(is:ie, js:je, config%nlay))

   end subroutine allocate_physics_arrays

end module mod_physics
