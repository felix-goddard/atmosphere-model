module mod_sw_dyn

   use mod_kinds, only: ik, rk
   use mod_log, only: logger => main_logger, log_str
   use mod_config, only: config => main_config
   use mod_constants, only: gravity, coriolis_parameter, kappa, cp_dry, &
                            divergence_damping_coeff, vorticity_damping_coeff
   use mod_tiles, only: is, ie, js, je, isd, ied, jsd, jed
   use mod_fields, only: dp, pt, ud, vd, plev, pkap, playkap, gz, net_flux, &
                         coord_Ak, coord_Bk
   use mod_sync, only: halo_exchange
   use mod_util, only: abort_now

   implicit none
   private

   public :: init_sw_dyn, &
             cgrid_dynamics_step, cgrid_halo_exchange, &
             dgrid_dynamics_step, dgrid_halo_exchange

   real(rk) :: rdx, rdx2, rdy, rdy2, rdx2rdy2, rA

   real(rk), allocatable :: &
      dpc(:, :, :), & ! pressure thickness for the C grid half-step
      ptc(:, :, :), & ! potential temperature for the C grid half-step
      uc(:, :, :), & ! u wind on C grid
      vc(:, :, :) ! v wind on C grid

   real(rk), allocatable :: vorticity(:, :), kinetic_energy(:, :)

   real(rk) :: divdamp_coeff, vortdamp_xcoeff, vortdamp_ycoeff
   real(rk), allocatable :: divergence(:, :), laplacian(:, :)

   real(rk), allocatable :: &
      pgfx(:, :), pgfy(:, :), & ! pressure gradient forces
      pkapb(:, :, :), gzb(:, :, :), & ! pkap and gz interpolated to B grid
      heating_rate(:, :, :) ! potential temperature heating rate (K/s)

   real(rk), allocatable :: &
      ua(:, :), & ! u wind on A grid
      va(:, :), & ! v wind on A grid
      ub(:, :), & ! u wind on B grid
      vb(:, :) ! v wind on B grid

   ! These variables are correlated as follows:
   !
   !        │                          │
   ! ─── ub(i,j) ──── ud(i,j) ──── ub(i+1,j) ──── ud(i+1,j) ────
   !     vb(i,j)      vc(i,j)      vb(i+1,j)      vc(i+1,j)
   !        │                          │
   !     vd(i,j)      ua(i,j)      vd(i+1,j)      ua(i+1,j)
   !     uc(i,j)      va(i,j)      uc(i+1,j)      va(i+1,j)
   !        │                          │
   ! ─── ub(i,j+1) ── ud(i,j+1) ── ub(i+1,j+1) ── ud(i+1,j+1) ──
   !     vb(i,j+1)    vc(i,j+1)    vb(i+1,j+1)    vc(i+1,j+1)
   !        │                          │
   !     vd(i,j+1)    ua(i,j+1)    vd(i+1,j+1)    ua(i+1,j+1)
   !     uc(i,j+1)    va(i,j+1)    uc(i+1,j+1)    va(i+1,j+1)
   !        │                          │
   !
   ! i.e. the A grid winds (ua and va, also dp) are located at mass
   ! points (cell centers); the B grid winds are located at cell
   ! corners; the C grid winds are perpendicular to cell edges,
   ! and the prognostic D grid winds are parallel to cell edges.

   real(rk), allocatable :: tmp(:, :)

   real(rk), allocatable :: fx(:), fy(:), cx(:, :), cy(:, :)
   real(rk), allocatable :: dpfx(:, :), dpfy(:, :) ! mass fluxes
   real(rk), allocatable :: denom_x(:, :), denom_y(:, :)
   real(rk), allocatable :: edge_L(:), edge_R(:) ! values calculated by the PPM

   integer(ik), parameter :: PIECEWISE_CONSTANT = 0
   integer(ik), parameter :: PPM_UNCONSTRAINED = 1
   integer(ik), parameter :: PPM_CONSTRAINED = 2
   integer(ik), parameter :: PPM_NONNEGATIVE = 3

contains

   subroutine cgrid_dynamics_step(dt)
      real(rk), intent(in) :: dt
      integer(ik) :: i, j, k
      real(rk) :: dtdx, dtdy

      dtdx = dt*rdx
      dtdy = dt*rdy

      ! calculate pressure on the layer interfaces
      do k = config%nlay, 1, -1
         plev(is:ie, js:je, k) = plev(is:ie, js:je, k + 1) + dp(is:ie, js:je, k)
         pkap(is:ie, js:je, k) = plev(is:ie, js:je, k)**kappa
      end do

      do k = 1, config%nlay

         ! Calculate the rate of change of potential temperature for this layer.
         ! We only need this over is+4:ie-4 (and same for j) because that's the
         ! range that the D grid potential temperature gets updated on
         playkap(is:ie, js:je, k) = &
            (pkap(is:ie, js:je, k + 1) - pkap(is:ie, js:je, k)) &
            /log(pkap(is:ie, js:je, k + 1)/pkap(is:ie, js:je, k))

         heating_rate(is:ie, js:je, k) = &
            (net_flux(is:ie, js:je, k + 1) - net_flux(is:ie, js:je, k)) &
            /(cp_dry*dp(is:ie, js:je, k)/gravity) &
            /playkap(is:ie, js:je, k)

         ! Calculate the winds
         do concurrent(i=is + 1:ie - 1, j=js + 1:je - 1)
            ub(i, j) = .5*(ud(i, j, k) + ud(i - 1, j, k))
            vb(i, j) = .5*(vd(i, j, k) + vd(i, j - 1, k))

            uc(i, j, k) = .5*(ub(i, j) + .5*( &
                              ud(i, j + 1, k) + ud(i - 1, j + 1, k)))

            vc(i, j, k) = .5*(vb(i, j) + .5*( &
                              vd(i + 1, j, k) + vd(i + 1, j - 1, k)))

            ua(i, j) = .5*(uc(i, j, k) + .25*( &
                           ud(i + 1, j, k) + ud(i, j, k) &
                           + ud(i + 1, j + 1, k) + ud(i, j + 1, k)))

            va(i, j) = .5*(vc(i, j, k) + .25*( &
                           vd(i, j + 1, k) + vd(i, j, k) &
                           + vd(i + 1, j + 1, k) + vd(i + 1, j, k)))
         end do

         ! ===========================================================================
         ! Calculate pressure thickness field update from the C grid winds (`dpc`)

         do concurrent(i=is + 1:ie - 1, j=js + 1:je - 1)
            cx(i, j) = dtdx*uc(i, j, k)
            cy(i, j) = dtdx*vc(i, j, k)
         end do

         do concurrent(i=is + 2:ie - 2, j=js + 2:je - 2)
            denom_x(i, j) = 1./(1.-dtdx*(uc(i + 1, j, k) - uc(i, j, k)))
            denom_y(i, j) = 1./(1.-dtdy*(vc(i, j + 1, k) - vc(i, j, k)))
         end do

         ! calculate the inner step (in the y-direction) for the outer x-direction step
         do i = is + 2, ie - 2
            call flux(fy(js + 1:je - 1), dp(i, js + 1:je - 1, k), &
                      cy(i, js + 1:je - 1), js + 1, je - 1, &
                      variant=PIECEWISE_CONSTANT)

            do j = js + 2, je - 2
               tmp(i, j) = .5*(dp(i, j, k) + denom_y(i, j)*( &
                               dp(i, j, k) - rA*(fy(j + 1) - fy(j))))

               laplacian(i, j) = -2.*dp(i, j, k)*rdx2rdy2 &
                                 + (dp(i + 1, j, k) + dp(i - 1, j, k))*rdx2 &
                                 + (dp(i, j + 1, k) + dp(i, j - 1, k))*rdy2
            end do
         end do

         ! calculate the outer step in the x-direction
         do j = js + 5, je - 5
            call flux(fx(is + 2:ie - 2), tmp(is + 2:ie - 2, j), &
                      cx(is + 2:ie - 2, j), is + 2, ie - 2, &
                      variant=PPM_UNCONSTRAINED)

            dpfx(is + 5, j) = fx(is + 5)*cx(is + 5, j)
            do i = is + 5, ie - 5
               dpfx(i + 1, j) = fx(i + 1)*cx(i + 1, j)

               dpc(i, j, k) = dp(i, j, k) - ( &
                              dpfx(i + 1, j) - dpfx(i, j) &
                              + dtdx*vortdamp_xcoeff*(laplacian(i + 1, j) &
                                                      - 2.*laplacian(i, j) &
                                                      + laplacian(i - 1, j)))
            end do
         end do

         ! calculate the inner step (in the x-direction) for the outer y-direction step
         do j = js + 2, je - 2
            call flux(fx(is + 1:ie - 1), dp(is + 1:ie - 1, j, k), &
                      cx(is + 1:ie - 1, j), is + 1, ie - 1, &
                      variant=PIECEWISE_CONSTANT)

            do i = is + 2, ie - 2
               tmp(i, j) = .5*(dp(i, j, k) + denom_x(i, j)*( &
                               dp(i, j, k) - rA*(fx(i + 1) - fx(i))))
            end do
         end do

         ! calculate the outer step in the y-direction
         do i = is + 5, ie - 5
            call flux(fy(js + 2:je - 2), tmp(i, js + 2:je - 2), &
                      cy(i, js + 2:je - 2), js + 2, je - 2, &
                      variant=PPM_UNCONSTRAINED)

            dpfy(i, js + 5) = fy(js + 5)*cy(i, js + 5)
            do j = js + 5, je - 5
               dpfy(i, j + 1) = fy(j + 1)*cy(i, j + 1)

               dpc(i, j, k) = dpc(i, j, k) - ( &
                              dpfy(i, j + 1) - dpfy(i, j) &
                              + dtdy*vortdamp_ycoeff*(laplacian(i, j + 1) &
                                                      - 2.*laplacian(i, j) &
                                                      + laplacian(i, j - 1)))
            end do
         end do

         ! ===========================================================================
         ! Calculate potential temperature update from the C grid winds (`ptc`)

         ! calculate the inner step (in the y-direction) for the outer x-direction step
         do i = is + 2, ie - 2
            call flux(fy(js + 1:je - 1), pt(i, js + 1:je - 1, k), &
                      cy(i, js + 1:je - 1), js + 1, je - 1, &
                      variant=PIECEWISE_CONSTANT)

            do j = js + 2, je - 2
               tmp(i, j) = .5*(pt(i, j, k) + denom_y(i, j) &
                               *(pt(i, j, k) - rA*(fy(j + 1) - fy(j))))

               laplacian(i, j) = -2.*pt(i, j, k)*rdx2rdy2 &
                                 + (pt(i + 1, j, k) + pt(i - 1, j, k))*rdx2 &
                                 + (pt(i, j + 1, k) + pt(i, j - 1, k))*rdy2
            end do
         end do

         ! calculate the outer step in the x-direction
         do j = js + 5, je - 5
            call flux(fx(is + 2:ie - 2), tmp(is + 2:ie - 2, j), &
                      cx(is + 2:ie - 2, j), is + 2, ie - 2, &
                      variant=PPM_UNCONSTRAINED)

            do i = is + 5, ie - 5
               ptc(i, j, k) = &
                  dp(i, j, k)*pt(i, j, k) &
                  - (fx(i + 1)*dpfx(i + 1, j) - fx(i)*dpfx(i, j) &
                     + vortdamp_xcoeff*( &
                     (laplacian(i + 1, j) - laplacian(i, j))*dpfx(i + 1, j) &
                     - (laplacian(i, j) - laplacian(i - 1, j))*dpfx(i, j)))
            end do
         end do

         ! calculate the inner step (in the x-direction) for the outer y-direction step
         do j = js + 2, je - 2
            call flux(fx(is + 1:ie - 1), pt(is + 1:ie - 1, j, k), &
                      cx(is + 1:ie - 1, j), is + 1, ie - 1, &
                      variant=PIECEWISE_CONSTANT)

            do i = is + 2, ie - 2
               tmp(i, j) = .5*(pt(i, j, k) + denom_x(i, j) &
                               *(pt(i, j, k) - rA*(fx(i + 1) - fx(i))))
            end do
         end do

         ! calculate the outer step in the y-direction
         do i = is + 5, ie - 5
            call flux(fy(js + 2:je - 2), tmp(i, js + 2:je - 2), &
                      cy(i, js + 2:je - 2), js + 2, je - 2, &
                      variant=PPM_UNCONSTRAINED)

            do j = js + 5, je - 5
               ptc(i, j, k) = &
                  (ptc(i, j, k) &
                   - (fy(j + 1)*dpfy(i, j + 1) - fy(j)*dpfy(i, j) &
                      + vortdamp_ycoeff*( &
                      (laplacian(i, j + 1) - laplacian(i, j))*dpfy(i, j + 1) &
                      - (laplacian(i, j) - laplacian(i, j - 1))*dpfy(i, j))) &
                   + dt*dp(i, j, k)*heating_rate(i, j, k) &
                   )/dpc(i, j, k)
            end do
         end do

      end do ! levels loop

      ! ===========================================================================
      ! Pressure and geopotential

      ! calculate pressure on the layer interfaces
      do k = config%nlay, 1, -1
         plev(is + 5:ie - 5, js + 5:je - 5, k) = &
            plev(is + 5:ie - 5, js + 5:je - 5, k + 1) &
            + dpc(is + 5:ie - 5, js + 5:je - 5, k)

         pkap(is + 5:ie - 5, js + 5:je - 5, k) = &
            plev(is + 5:ie - 5, js + 5:je - 5, k)**kappa
      end do

      ! calculate the geopotential height of the interfaces
      do k = 1, config%nlay
         gz(is + 5:ie - 5, js + 5:je - 5, k + 1) = &
            gz(is + 5:ie - 5, js + 5:je - 5, k) &
            + cp_dry*ptc(is + 5:ie - 5, js + 5:je - 5, k)*( &
            pkap(is + 5:ie - 5, js + 5:je - 5, k) &
            - pkap(is + 5:ie - 5, js + 5:je - 5, k + 1))
      end do

      do k = 1, config%nlay

         ! ===========================================================================
         ! Kinetic energy, vorticity, and pressure gradient force

         do i = is + 5, ie - 5
            call flux(fy(js + 2:je - 2), vc(i, js + 3:je - 1, k), &
                      dtdy*va(i, js + 2:je - 2), js + 2, je - 2, &
                      variant=PPM_UNCONSTRAINED)

            do j = js + 5, je - 5
               kinetic_energy(i, j) = va(i, j)*fy(j)
            end do
         end do

         do j = js + 5, je - 5
            call flux(fx(is + 2:ie - 2), uc(is + 3:ie - 1, j, k), &
                      dtdx*ua(is + 2:ie - 2, j), is + 2, ie - 2, &
                      variant=PPM_UNCONSTRAINED)

            do i = is + 5, ie - 5
               kinetic_energy(i, j) = .5*(kinetic_energy(i, j) + ua(i, j)*fx(i))
            end do
         end do

         do concurrent(i=is + 2:ie - 2, j=js + 2:je - 2)
            vorticity(i, j) = coriolis_parameter &
                              - (uc(i, j, k) - uc(i, j - 1, k))*rdy &
                              + (vc(i, j, k) - vc(i - 1, j, k))*rdx
         end do

         do concurrent(i=is + 7:ie - 7, j=js + 7:je - 7)
            pgfx(i, j) = ((pkap(i, j, k) - pkap(i - 1, j, k + 1)) &
                          *(gz(i - 1, j, k) - gz(i, j, k + 1)) &
                          + (pkap(i - 1, j, k) - pkap(i, j, k + 1)) &
                          *(gz(i - 1, j, k + 1) - gz(i, j, k))) &
                         /(pkap(i, j, k) - pkap(i - 1, j, k + 1) &
                           + pkap(i - 1, j, k) - pkap(i, j, k + 1))

            pgfy(i, j) = ((pkap(i, j, k) - pkap(i, j - 1, k + 1)) &
                          *(gz(i, j - 1, k) - gz(i, j, k + 1)) &
                          + (pkap(i, j - 1, k) - pkap(i, j, k + 1)) &
                          *(gz(i, j - 1, k + 1) - gz(i, j, k))) &
                         /(pkap(i, j, k) - pkap(i, j - 1, k + 1) &
                           + pkap(i, j - 1, k) - pkap(i, j, k + 1))
         end do

         ! ===========================================================================
         ! Calculate C grid u (`uc`)

         do concurrent(i=is + 1:ie - 1, j=js + 1:je - 1)
            denom_x(i, j) = 1./(1.-dtdx*(ud(i, j, k) - ud(i - 1, j, k)))
            denom_y(i, j) = 1./(1.-dtdy*(vd(i, j, k) - vd(i, j - 1, k)))

            cx(i, j) = dtdx*ud(i, j, k)
            cy(i, j) = dtdx*vd(i, j, k)
         end do

         ! calculate the inner step (in the x-direction) for the outer y-direction step
         do j = js + 3, je - 3
            call flux(fx(is + 2:ie - 2), vorticity(is + 2:ie - 2, j), &
                      cx(is + 1:ie - 3, j), is + 2, ie - 2, &
                      variant=PIECEWISE_CONSTANT)

            do i = is + 3, ie - 3
               tmp(i, j) = .5*(vorticity(i, j) + denom_x(i, j) &
                               *(vorticity(i, j) - rA*(fx(i + 1) - fx(i))))

               laplacian(i, j) = &
                  -2.*vorticity(i, j)*rdx2rdy2 &
                  + (vorticity(i + 1, j) + vorticity(i - 1, j))*rdx2 &
                  + (vorticity(i, j + 1) + vorticity(i, j - 1))*rdy2
            end do
         end do

         do i = is + 7, ie - 7
            call flux(fy(js + 4:je - 4), tmp(i, js + 5:je - 3), &
                      cy(i, js + 4:je - 4), js + 4, je - 4, &
                      variant=PPM_UNCONSTRAINED)

            do j = js + 7, je - 7
               uc(i, j, k) = &
                  uc(i, j, k) &
                  + dt*(vd(i, j, k)*fy(j) + vortdamp_ycoeff*( &
                        laplacian(i, j + 1) - laplacian(i, j))) &
                  - dtdx*(kinetic_energy(i, j) - kinetic_energy(i - 1, j) &
                          - pgfx(i, j))
            end do
         end do

         ! ===========================================================================
         ! Calculate C grid v (`vc`)

         ! calculate the inner step (in the y-direction) for the outer x-direction step
         do i = is + 3, ie - 3
            call flux(fy(js + 2:je - 2), vorticity(i, js + 2:je - 2), &
                      cy(i, js + 1:je - 3), js + 2, je - 2, &
                      variant=PIECEWISE_CONSTANT)

            do j = js + 3, je - 3
               tmp(i, j) = .5*(vorticity(i, j) + denom_y(i, j) &
                               *(vorticity(i, j) - rA*(fy(j + 1) - fy(j))))
            end do
         end do

         do j = js + 7, je - 7
            call flux(fx(is + 4:ie - 4), tmp(is + 5:ie - 3, j), &
                      cx(is + 4:ie - 4, j), is + 4, ie - 4, &
                      variant=PPM_UNCONSTRAINED)

            do i = is + 7, ie - 7
               vc(i, j, k) = &
                  vc(i, j, k) &
                  - dt*(ud(i, j, k)*fx(i) + vortdamp_xcoeff*( &
                        laplacian(i + 1, j) - laplacian(i, j))) &
                  - dtdy*(kinetic_energy(i, j) - kinetic_energy(i, j - 1) &
                          - pgfy(i, j))
            end do
         end do

      end do ! levels loop

   end subroutine cgrid_dynamics_step

   subroutine cgrid_halo_exchange()

      call halo_exchange(dpc, ptc, uc, vc)

   end subroutine cgrid_halo_exchange

   subroutine dgrid_dynamics_step(dt)
      real(rk), intent(in) :: dt
      integer(ik) :: i, j, k
      real(rk) :: dtdx, dtdy

      dtdx = dt*rdx
      dtdy = dt*rdy

      ! calculate pressure on the layer interfaces
      do k = config%nlay, 1, -1
         plev(is:ie, js:je, k) = plev(is:ie, js:je, k + 1) &
                                 + dpc(is:ie, js:je, k)
         pkap(is:ie, js:je, k) = plev(is:ie, js:je, k)**kappa
      end do

      do k = 1, config%nlay

         do concurrent(i=is + 1:ie, j=js + 1:je)
            ub(i, j) = .5*(uc(i, j, k) + uc(i, j - 1, k))
            vb(i, j) = .5*(vc(i, j, k) + vc(i - 1, j, k))
         end do

         do concurrent(i=is:ie, j=js:je)
            cx(i, j) = dtdx*uc(i, j, k)
            cy(i, j) = dtdx*vc(i, j, k)
         end do

         do concurrent(i=is + 1:ie - 1, j=js + 1:je - 1)
            denom_x(i, j) = 1./(1.-dtdx*(uc(i + 1, j, k) - uc(i, j, k)))
            denom_y(i, j) = 1./(1.-dtdy*(vc(i, j + 1, k) - vc(i, j, k)))
         end do

         ! Calculate the rate of change of potential temperature for this layer.
         playkap(is:ie, js:je, k) = &
            (pkap(is:ie, js:je, k + 1) - pkap(is:ie, js:je, k)) &
            /log(pkap(is:ie, js:je, k + 1)/pkap(is:ie, js:je, k))

         heating_rate(is:ie, js:je, k) = &
            (net_flux(is:ie, js:je, k + 1) - net_flux(is:ie, js:je, k)) &
            /(cp_dry*dpc(is:ie, js:je, k)/gravity) &
            /playkap(is:ie, js:je, k)

         ! ===========================================================================
         ! Update the prognostic dp field

         ! calculate the inner step (in the y-direction) for the outer x-direction step
         do i = is + 1, ie - 1
            call flux(fy(js:je), dp(i, js:je, k), cy(i, js:je), js, je, &
                      variant=PIECEWISE_CONSTANT)

            do j = js + 1, je - 1
               tmp(i, j) = .5*(dp(i, j, k) + denom_y(i, j) &
                               *(dp(i, j, k) - rA*(fy(j + 1) - fy(j))))

               laplacian(i, j) = -2.*dp(i, j, k)*rdx2rdy2 &
                                 + (dp(i + 1, j, k) + dp(i - 1, j, k))*rdx2 &
                                 + (dp(i, j + 1, k) + dp(i, j - 1, k))*rdy2
            end do
         end do

         ! calculate the outer step in the x-direction
         do j = js + 4, je - 4
            call flux(fx(is + 1:ie - 1), tmp(is + 1:ie - 1, j), &
                      cx(is + 1:ie - 1, j), is + 1, ie - 1, &
                      variant=PPM_NONNEGATIVE)

            dpfx(is + 4, j) = fx(is + 4)*cx(is + 4, j)
            do i = is + 4, ie - 4
               dpfx(i + 1, j) = fx(i + 1)*cx(i + 1, j)

               dpc(i, j, k) = dp(i, j, k) - ( &
                              dpfx(i + 1, j) - dpfx(i, j) &
                              + dtdx*vortdamp_xcoeff*(laplacian(i + 1, j) &
                                                      - 2.*laplacian(i, j) &
                                                      + laplacian(i - 1, j)))
            end do
         end do

         ! calculate the inner step (in the x-direction) for the outer y-direction step
         do j = js + 1, je - 1
            call flux(fx(is:ie), dp(is:ie, j, k), cx(is:ie, j), is, ie, &
                      variant=PIECEWISE_CONSTANT)

            do i = is + 1, ie - 1
               tmp(i, j) = .5*(dp(i, j, k) + denom_x(i, j) &
                               *(dp(i, j, k) - rA*(fx(i + 1) - fx(i))))
            end do
         end do

         ! calculate the outer step in the y-direction
         do i = is + 4, ie - 4
            call flux(fy(js + 1:je - 1), tmp(i, js + 1:je - 1), &
                      cy(i, js + 1:je - 1), js + 1, je - 1, &
                      variant=PPM_NONNEGATIVE)

            dpfy(i, js + 4) = fy(js + 4)*cy(i, js + 4)
            do j = js + 4, je - 4
               dpfy(i, j + 1) = fy(j + 1)*cy(i, j + 1)

               dpc(i, j, k) = dpc(i, j, k) - ( &
                              dpfy(i, j + 1) - dpfy(i, j) &
                              + dtdy*vortdamp_ycoeff*(laplacian(i, j + 1) &
                                                      - 2.*laplacian(i, j) &
                                                      + laplacian(i, j - 1)))
            end do
         end do

         ! ===========================================================================
         ! Update the prognostic potential temperature field

         ! calculate the inner step (in the y-direction) for the outer x-direction step
         do i = is + 1, ie - 1
            call flux(fy(js:je), pt(i, js:je, k), cy(i, js:je), js, je, &
                      variant=PIECEWISE_CONSTANT)

            do j = js + 1, je - 1
               tmp(i, j) = .5*(pt(i, j, k) + denom_y(i, j) &
                               *(pt(i, j, k) - rA*(fy(j + 1) - fy(j))))

               laplacian(i, j) = -2.*pt(i, j, k)*rdx2rdy2 &
                                 + (pt(i + 1, j, k) + pt(i - 1, j, k))*rdx2 &
                                 + (pt(i, j + 1, k) + pt(i, j - 1, k))*rdy2
            end do
         end do

         ! calculate the outer step in the x-direction
         do j = js + 4, je - 4
            call flux(fx(is + 1:ie - 1), tmp(is + 1:ie - 1, j), &
                      cx(is + 1:ie - 1, j), is + 1, ie - 1, &
                      variant=PPM_NONNEGATIVE)

            do i = is + 4, ie - 4
               ptc(i, j, k) = &
                  dp(i, j, k)*pt(i, j, k) &
                  - (fx(i + 1)*dpfx(i + 1, j) - fx(i)*dpfx(i, j) &
                     + vortdamp_xcoeff*( &
                     (laplacian(i + 1, j) - laplacian(i, j))*dpfx(i + 1, j) &
                     - (laplacian(i, j) - laplacian(i - 1, j))*dpfx(i, j)))
            end do
         end do

         ! calculate the inner step (in the x-direction) for the outer y-direction step
         do j = js + 1, je - 1
            call flux(fx(is:ie), pt(is:ie, j, k), cx(is:ie, j), is, ie, &
                      variant=PIECEWISE_CONSTANT)

            do i = is + 1, ie - 1
               tmp(i, j) = .5*(pt(i, j, k) + denom_x(i, j) &
                               *(pt(i, j, k) - rA*(fx(i + 1) - fx(i))))
            end do
         end do

         ! calculate the outer step in the y-direction
         do i = is + 4, ie - 4
            call flux(fy(js + 1:je - 1), tmp(i, js + 1:je - 1), &
                      cy(i, js + 1:je - 1), js + 1, je - 1, &
                      variant=PPM_NONNEGATIVE)

            do j = js + 4, je - 4
               ptc(i, j, k) = &
                  (ptc(i, j, k) &
                   - (fy(j + 1)*dpfy(i, j + 1) - fy(j)*dpfy(i, j) &
                      + vortdamp_ycoeff*( &
                      (laplacian(i, j + 1) - laplacian(i, j))*dpfy(i, j + 1) &
                      - (laplacian(i, j) - laplacian(i, j - 1))*dpfy(i, j))) &
                   + dt*dp(i, j, k)*heating_rate(i, j, k) &
                   )/dpc(i, j, k)
            end do
         end do

         dp(is + 4:ie - 4, js + 4:je - 4, k) = &
            dpc(is + 4:ie - 4, js + 4:je - 4, k)

         pt(is + 4:ie - 4, js + 4:je - 4, k) = &
            ptc(is + 4:ie - 4, js + 4:je - 4, k)

      end do ! levels loop

      ! ===========================================================================
      ! Pressure and geopotential

      ! calculate pressure on the layer interfaces
      do k = config%nlay, 1, -1
         plev(is + 4:ie - 4, js + 4:je - 4, k) = &
            plev(is + 4:ie - 4, js + 4:je - 4, k + 1) &
            + dp(is + 4:ie - 4, js + 4:je - 4, k)

         pkap(is + 4:ie - 4, js + 4:je - 4, k) = &
            plev(is + 4:ie - 4, js + 4:je - 4, k)**kappa

         call interpolate_to_corners( &
            pkapb(is + 4:ie - 4, js + 4:je - 4, k), &
            pkap(is + 4:ie - 4, js + 4:je - 4, k), &
            is + 4, ie - 4, js + 4, je - 4)
      end do

      ! calculate the geopotential height of the interfaces
      do k = 1, config%nlay
         gz(is + 4:ie - 4, js + 4:je - 4, k + 1) = &
            gz(is + 4:ie - 4, js + 4:je - 4, k) &
            + cp_dry*pt(is + 4:ie - 4, js + 4:je - 4, k)*( &
            pkap(is + 4:ie - 4, js + 4:je - 4, k) &
            - pkap(is + 4:ie - 4, js + 4:je - 4, k + 1))

         call interpolate_to_corners( &
            gzb(is + 4:ie - 4, js + 4:je - 4, k + 1), &
            gz(is + 4:ie - 4, js + 4:je - 4, k + 1), &
            is + 4, ie - 4, js + 4, je - 4)
      end do

      do k = 1, config%nlay

         ! ===========================================================================
         ! Kinetic energy, vorticity, and pressure gradient force

         do i = is + 4, ie - 4
            call flux(fy(js + 1:je - 1), vd(i, js + 1:je - 1, k), &
                      dtdy*vb(i, js + 1:je - 1), js + 1, je - 1, &
                      variant=PPM_CONSTRAINED)

            do j = js + 4, je - 4
               kinetic_energy(i, j) = vb(i, j)*fy(j)
            end do
         end do

         do j = js + 4, je - 4
            call flux(fx(is + 1:ie - 1), ud(is + 1:ie - 1, j, k), &
                      dtdx*ub(is + 1:ie - 1, j), is + 1, ie - 1, &
                      variant=PPM_CONSTRAINED)

            do i = is + 4, ie - 4
               kinetic_energy(i, j) = .5*(kinetic_energy(i, j) + ub(i, j)*fx(i))
            end do
         end do

         do concurrent(i=is + 1:ie - 1, j=js + 1:je - 1)
            vorticity(i, j) = coriolis_parameter &
                              - (ud(i, j + 1, k) - ud(i, j, k))*rdy &
                              + (vd(i + 1, j, k) - vd(i, j, k))*rdx
         end do

         do concurrent(i=is + 7:ie - 7, j=js + 7:je - 7)
            pgfx(i, j) = ((pkapb(i + 1, j, k) - pkapb(i, j, k + 1)) &
                          *(gzb(i, j, k) - gzb(i + 1, j, k + 1)) &
                          + (pkapb(i, j, k) - pkapb(i + 1, j, k + 1)) &
                          *(gzb(i, j, k + 1) - gzb(i + 1, j, k))) &
                         /(pkapb(i + 1, j, k) - pkapb(i, j, k + 1) &
                           + pkapb(i, j, k) - pkapb(i + 1, j, k + 1))

            pgfy(i, j) = ((pkapb(i, j + 1, k) - pkapb(i, j, k + 1)) &
                          *(gzb(i, j, k) - gzb(i, j + 1, k + 1)) &
                          + (pkapb(i, j, k) - pkapb(i, j + 1, k + 1)) &
                          *(gzb(i, j, k + 1) - gzb(i, j + 1, k))) &
                         /(pkapb(i, j + 1, k) - pkapb(i, j, k + 1) &
                           + pkapb(i, j, k) - pkapb(i, j + 1, k + 1))
         end do

         ! compute the divergence of the prognostic winds
         do concurrent(i=is + 4:ie - 4, j=js + 4:je - 4)
            divergence(i, j) = rdx*(ud(i, j, k) - ud(i - 1, j, k)) &
                               + rdy*(vd(i, j, k) - vd(i, j - 1, k))
         end do

         ! calculate ∇²(divergence)
         do concurrent(i=is + 5:ie - 5, j=js + 5:je - 5)
            tmp(i, j) = -2.*divergence(i, j)*rdx2rdy2 &
                        + (divergence(i + 1, j) + divergence(i - 1, j))*rdx2 &
                        + (divergence(i, j + 1) + divergence(i, j - 1))*rdy2
         end do

         ! calculate ∇⁴(divergence)
         do concurrent(i=is + 6:ie - 6, j=js + 6:je - 6)
            divergence(i, j) = -2.*tmp(i, j)*rdx2rdy2 &
                               + (tmp(i + 1, j) + tmp(i - 1, j))*rdx2 &
                               + (tmp(i, j + 1) + tmp(i, j - 1))*rdy2
         end do

         divergence(is + 6:ie - 6, js + 6:je - 6) = &
            divdamp_coeff*divergence(is + 6:ie - 6, js + 6:je - 6)

         ! ===========================================================================
         ! Update the prognostic D grid u (`ud`)

         ! calculate the inner step (in the x-direction) for the outer y-direction step
         do j = js + 4, je - 4
            call flux(fx(is + 1:ie - 1), vorticity(is + 1:ie - 1, j), &
                      cx(is + 1:ie - 1, j), is + 1, ie - 1, &
                      variant=PIECEWISE_CONSTANT)

            do i = is + 4, ie - 4
               tmp(i, j) = .5*(vorticity(i, j) + denom_x(i, j)*( &
                               vorticity(i, j) - rA*(fx(i + 1) - fx(i))))

               laplacian(i, j) = &
                  -2.*vorticity(i, j)*rdx2rdy2 &
                  + (vorticity(i + 1, j) + vorticity(i - 1, j))*rdx2 &
                  + (vorticity(i, j + 1) + vorticity(i, j - 1))*rdy2
            end do
         end do

         do i = is + 7, ie - 7
            call flux(fy(js + 4:je - 4), tmp(i, js + 4:je - 4), &
                      cy(i, js + 4:je - 4), js + 4, je - 4, &
                      variant=PPM_CONSTRAINED)

            do j = js + 7, je - 7
               ud(i, j, k) = &
                  ud(i, j, k) &
                  + dt*(vc(i, j, k)*fy(j) + vortdamp_ycoeff*( &
                        laplacian(i, j) - laplacian(i, j - 1))) &
                  - dtdx*(kinetic_energy(i + 1, j) - kinetic_energy(i, j) &
                          - divergence(i + 1, j) + divergence(i, j) &
                          - pgfx(i, j))
            end do
         end do

         ! ===========================================================================
         ! Update the prognostic D grid v (`vd`)

         ! calculate the inner step (in the y-direction) for the outer x-direction step
         do i = is + 4, ie - 4
            call flux(fy(js + 1:je - 1), vorticity(i, js + 1:je - 1), &
                      cy(i, js + 1:je - 1), js + 1, je - 1, &
                      variant=PIECEWISE_CONSTANT)

            do j = js + 4, je - 4
               tmp(i, j) = .5*(vorticity(i, j) + denom_y(i, j)*( &
                               vorticity(i, j) - rA*(fy(j + 1) - fy(j))))
            end do
         end do

         do j = js + 7, je - 7
            call flux(fx(is + 4:ie - 4), tmp(is + 4:ie - 4, j), &
                      cx(is + 4:ie - 4, j), is + 4, ie - 4, &
                      variant=PPM_CONSTRAINED)

            do i = is + 7, ie - 7
               vd(i, j, k) = &
                  vd(i, j, k) &
                  - dt*(uc(i, j, k)*fx(i) + vortdamp_xcoeff*( &
                        laplacian(i, j) - laplacian(i - 1, j))) &
                  - dtdy*(kinetic_energy(i, j + 1) - kinetic_energy(i, j) &
                          - divergence(i, j + 1) + divergence(i, j) &
                          - pgfy(i, j))
            end do
         end do

      end do ! levels loop

   end subroutine dgrid_dynamics_step

   subroutine dgrid_halo_exchange()

      call halo_exchange(dp, pt, ud, vd)

   end subroutine dgrid_halo_exchange

   subroutine flux(f, q, c, isl, iel, variant)
      ! Calculate numerical fluxes of q given an array of Courant numbers c, assuming
      ! c is defined on points staggered relative to q. This is done by calculating
      ! a piecewise parabolic interpolation of q, and then analytically integrating
      ! this subgrid distribution.

      integer(ik), intent(in) :: isl, iel, variant
      real(rk), intent(in)    :: q(isl:iel) ! quantity to calculate fluxes of
      real(rk), intent(in)    :: c(isl:iel) ! Courant numbers
      real(rk), intent(inout) :: f(isl:iel)
      integer(ik) :: i
      real(rk) :: s, k1, k2

      if (variant == PPM_UNCONSTRAINED &
          .or. variant == PPM_CONSTRAINED &
          .or. variant == PPM_NONNEGATIVE) then

         call ppm(q, isl, iel, variant)

         ! ppm defines the edges values over isl+2:iel-2 so we can calculate
         ! fluxes over isl+3:iel-2
         do i = isl + 3, iel - 2
            s = abs(c(i))
            k1 = (1.-2.*s)*(1.-s)
            k2 = s*(s - 1.)

            ! apply monotonicity constraint: if
            !   3 | left + right | > | left - right |
            ! holds for both cells on this interface, we fall back to a
            ! linear upwind scheme; otherwise we calculate the upwind PPM
            if (variant == PPM_CONSTRAINED &
                .and. (3.*abs(edge_L(i) + edge_R(i)) &
                       > abs(edge_L(i) - edge_R(i))) &
                .and. (3.*abs(edge_L(i - 1) + edge_R(i - 1)) &
                       > abs(edge_L(i - 1) - edge_R(i - 1)))) then

               if (c(i) > 0.) then
                  f(i) = q(i - 1)
               else
                  f(i) = q(i)
               end if
            else
               if (c(i) > 0.) then
                  f(i) = q(i - 1) + k1*edge_R(i - 1) + k2*edge_L(i - 1)
               else
                  f(i) = q(i) + k1*edge_L(i) + k2*edge_R(i)
               end if
            end if
         end do

      else if (variant == PIECEWISE_CONSTANT) then
         do i = isl + 1, iel
            if (c(i) > 0.) then
               f(i) = q(i - 1)
            else
               f(i) = q(i)
            end if
         end do

      else
         write (log_str, '(a,i4)') 'Unknown numerical flux variant:', variant
         call logger%fatal('ppm', log_str)
         call abort_now()
      end if

   end subroutine flux

   subroutine ppm(q, isl, iel, variant)
      ! Calculate, by a piecewise parabolic reconstruction, the left and right cell
      ! edge anomaly values given the cell averages q(isl:iel).
      ! The `variant` argument selects between different reconstruction algorithms,
      ! currently the options are:
      !  - PPM_UNCONSTRAINED: implements then simplified method of Lin (2004) without
      !                       any constraints on the edge values.
      !  - PPM_CONSTRAINED: also implements the simplified method of Lin (2004) with
      !                     the addition of a monotonicity constraint.

      integer(ik), intent(in) :: isl, iel
      real(rk), intent(in) :: q(isl:iel)
      integer(ik), intent(in) :: variant
      real(rk) :: dqi_min, dqi_max, a, b
      real(rk) :: dqi(isl + 1:iel - 1)
      real(rk) :: dqi_mono(isl + 1:iel - 1)
      integer(ik) :: i

      if (variant == PPM_UNCONSTRAINED &
          .or. variant == PPM_CONSTRAINED &
          .or. variant == PPM_NONNEGATIVE) then

         do i = isl + 1, iel - 1
            dqi(i) = (q(i + 1) - q(i - 1))/4.
            dqi_max = max(q(i + 1), q(i), q(i - 1)) - q(i)
            dqi_min = q(i) - min(q(i + 1), q(i), q(i - 1))
            dqi_mono(i) = sign(min(abs(dqi(i)), dqi_min, dqi_max), dqi(i))
         end do

         do i = isl + 2, iel - 2
            edge_L(i) = (q(i - 1) - q(i))/2.+(dqi_mono(i - 1) - dqi_mono(i))/3.
            edge_R(i) = (q(i + 1) - q(i))/2.+(dqi_mono(i) - dqi_mono(i + 1))/3.

            if (variant == PPM_CONSTRAINED &
                .or. variant == PPM_NONNEGATIVE) then

               edge_L(i) = -sign( &
                           min(2.*abs(dqi(i)), abs(edge_L(i))), dqi(i))

               edge_R(i) = +sign( &
                           min(2.*abs(dqi(i)), abs(edge_R(i))), dqi(i))

               if (variant == PPM_NONNEGATIVE) then
                  ! filter negative edge values
                  edge_L(i) = max(edge_L(i), -q(i))
                  edge_R(i) = max(edge_R(i), -q(i))

                  ! filter negative minima
                  a = edge_R(i) - edge_L(i)
                  b = -3.*(edge_R(i) + edge_L(i))
                  if ((abs(a) < -b) &
                      .and. (q(i) + a**2/(4*b) + b/12. < 0.)) then
                     if (edge_L(i)*edge_R(i) > 0.) then
                        edge_L(i) = 0.
                        edge_R(i) = 0.
                     else if (edge_L(i) > edge_R(i)) then
                        edge_L(i) = 2.*abs(edge_R(i))
                     else
                        edge_R(i) = 2.*abs(edge_L(i))
                     end if
                  end if
               end if
            end if
         end do

      else
         write (log_str, '(a,i4)') 'Unknown PPM variant:', variant
         call logger%fatal('ppm', log_str)
         call abort_now()
      end if

   end subroutine ppm

   subroutine interpolate_to_corners(out, q, isl, iel, jsl, jel)
      ! Interpolate a quantity defined at cell centers (the A grid)
      ! to cell corners (the B grid) using an average of Lagrange and PPM
      ! interpolations along alternate directions

      integer(ik), intent(in) :: isl, iel, jsl, jel
      real(rk), intent(inout) :: out(isl:iel, jsl:jel)
      real(rk), intent(in) :: q(isl:iel, jsl:jel)
      real(rk) :: tmpl(isl:iel, jsl:jel)

      ! PPM interpolation along x-direction
      tmpl(isl + 2:iel - 1, jsl:jel) = &
         (7./12.)*(q(isl + 2:iel - 1, jsl:jel) + q(isl + 1:iel - 2, jsl:jel)) &
         - (1./12.)*(q(isl + 3:iel, jsl:jel) + q(isl:iel - 3, jsl:jel))

      ! Lagrange interpolation along y-direction
      out(isl + 2:iel - 1, jsl + 2:jel - 1) = &
         (9./32.)*(tmpl(isl + 2:iel - 1, jsl + 2:jel - 1) &
                   + tmpl(isl + 2:iel - 1, jsl + 1:jel - 2)) &
         - (1./32.)*(tmpl(isl + 2:iel - 1, jsl + 3:jel) &
                     + tmpl(isl + 2:iel - 1, jsl:jel - 3))

      ! PPM interpolation along y-direction
      tmpl(isl:iel, jsl + 2:jel - 1) = &
         (7./12.)*(q(isl:iel, jsl + 2:jel - 1) + q(isl:iel, jsl + 1:jel - 2)) &
         - (1./12.)*(q(isl:iel, jsl + 3:jel) + q(isl:iel, jsl:jel - 3))

      ! Lagrange interpolation along x-direction
      out(isl + 2:iel - 1, jsl + 2:jel - 1) = &
         out(isl + 2:iel - 1, jsl + 2:jel - 1) + &
         (9./32.)*(tmpl(isl + 2:iel - 1, jsl + 2:jel - 1) &
                   + tmpl(isl + 1:iel - 2, jsl + 2:jel - 1)) &
         - (1./32.)*(tmpl(isl + 3:iel, jsl + 2:jel - 1) &
                     + tmpl(isl:iel - 3, jsl + 2:jel - 1))

   end subroutine interpolate_to_corners

   subroutine init_sw_dyn()
      real(rk) :: vortdamp_coeff

      rdx = 1./config%dx
      rdx2 = 1./(config%dx*config%dx)
      rdy = 1./config%dy
      rdy2 = 1./(config%dy*config%dy)
      rdx2rdy2 = rdx2 + rdy2
      rA = 1./(config%dx*config%dy)

      if (.not. allocated(dpc)) allocate (dpc(is:ie, js:je, 1:config%nlay))
      if (.not. allocated(ptc)) allocate (ptc(is:ie, js:je, 1:config%nlay))
      if (.not. allocated(uc)) allocate (uc(is:ie, js:je, 1:config%nlay))
      if (.not. allocated(vc)) allocate (vc(is:ie, js:je, 1:config%nlay))

      if (.not. allocated(vorticity)) allocate (vorticity(is:ie, js:je))
      if (.not. allocated(kinetic_energy)) &
         allocate (kinetic_energy(is:ie, js:je))

      if (.not. allocated(divergence)) allocate (divergence(is:ie, js:je))
      if (.not. allocated(laplacian)) allocate (laplacian(is:ie, js:je))

      if (.not. allocated(heating_rate)) &
         allocate (heating_rate(is:ie, js:je, 1:config%nlay))

      if (.not. allocated(pgfx)) allocate (pgfx(is:ie, js:je))
      if (.not. allocated(pgfy)) allocate (pgfy(is:ie, js:je))
      if (.not. allocated(pkapb)) &
         allocate (pkapb(is:ie, js:je, 1:config%nlay + 1))
      if (.not. allocated(gzb)) &
         allocate (gzb(is:ie, js:je, 1:config%nlay + 1))

      ! the top pressure is constant so no need to interpolate
      pkapb(is:ie, js:je, config%nlay + 1) = coord_Ak(1)**kappa

      ! the surface geopotential is constant so we can interpolate it once
      ! at the start and be done
      call interpolate_to_corners(gzb(is:ie, js:je, 1), gz(is:ie, js:je, 1), &
                                  is, ie, js, je)

      if (.not. allocated(ua)) allocate (ua(is:ie, js:je))
      if (.not. allocated(va)) allocate (va(is:ie, js:je))
      if (.not. allocated(ub)) allocate (ub(is:ie, js:je))
      if (.not. allocated(vb)) allocate (vb(is:ie, js:je))

      if (.not. allocated(tmp)) allocate (tmp(is:ie, js:je))

      if (.not. allocated(fx)) allocate (fx(min(is, js):max(ie, je)))
      if (.not. allocated(fy)) allocate (fy(min(is, js):max(ie, je)))
      if (.not. allocated(cx)) allocate (cx(is:ie, js:je))
      if (.not. allocated(cy)) allocate (cy(is:ie, js:je))
      if (.not. allocated(dpfx)) allocate (dpfx(is:ie, js:je))
      if (.not. allocated(dpfy)) allocate (dpfy(is:ie, js:je))
      if (.not. allocated(denom_x)) allocate (denom_x(is:ie, js:je))
      if (.not. allocated(denom_y)) allocate (denom_y(is:ie, js:je))

      if (.not. allocated(edge_L)) allocate (edge_L(min(is, js):max(ie, je)))
      if (.not. allocated(edge_R)) allocate (edge_R(min(is, js):max(ie, je)))

      divdamp_coeff = (divergence_damping_coeff*config%dx*config%dy)**2

      vortdamp_coeff = vorticity_damping_coeff*config%dx*config%dy
      vortdamp_xcoeff = config%dy*rdx*vortdamp_coeff
      vortdamp_ycoeff = config%dx*rdy*vortdamp_coeff

   end subroutine init_sw_dyn

end module mod_sw_dyn
