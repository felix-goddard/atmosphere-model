module mod_remap

   use mod_kinds, only: ik, rk
   use mod_config, only: config => main_config
   use mod_tiles, only: is, ie, js, je, isd, ied, jsd, jed
   use mod_constants, only: dry_heat_capacity, kappa, top_pressure
   use mod_fields, only: dp, pt, ud, vd, gz, plev, pkap, playkap

   implicit none
   private

   public :: allocate_remapping_arrays, perform_vertical_remapping

   real(rk), allocatable :: u_remap(:, :, :)
   real(rk), allocatable :: v_remap(:, :, :)
   real(rk), allocatable :: pt_remap(:, :, :)

   real(rk), allocatable :: target_plev(:, :, :), target_pkap(:, :, :), &
                            target_playkap(:, :, :)

   real(rk), allocatable :: btm_edge(:), top_edge(:), curvature(:)

contains

   subroutine perform_vertical_remapping()
      integer(ik) :: i, j, k, nlay, isl, iel, jsl, jel
      real(rk) :: src_plev(config%nlay + 1), dst_plev(config%nlay + 1), &
                  plog(isd - 3:ied + 3, jsd - 3:jed + 3), &
                  energy(config%nlay), energy_remap(config%nlay), &
                  gz_remap(config%nlay + 1), ua(config%nlay), va(config%nlay)

      nlay = config%nlay

      isl = isd - 3
      iel = ied + 3
      jsl = jsd - 3
      jel = jed + 3

      ! calculate the current vertical grid (pressure levels) and the target
      ! grid; we must ensure that the current and target grids coincide
      ! at the top and bottom
      do k = nlay, 1, -1
         plev(isl:iel, jsl:jel, k) = &
            plev(isl:iel, jsl:jel, k + 1) + dp(isl:iel, jsl:jel, k)

         pkap(isl:iel, jsl:jel, k) = plev(isl:iel, jsl:jel, k)**kappa

         plog(isl:iel, jsl:jel) = &
            log(pkap(isl:iel, jsl:jel, k + 1)/pkap(isl:iel, jsl:jel, k))

         playkap(isl:iel, jsl:jel, k) = &
            (pkap(isl:iel, jsl:jel, k + 1) - pkap(isl:iel, jsl:jel, k)) &
            /plog(isl:iel, jsl:jel)

         ! currently: remap so that every dp is equal
         target_plev(isl:iel, jsl:jel, k) = &
            target_plev(isl:iel, jsl:jel, k + 1) &
            + (plev(isl:iel, jsl:jel, 1) &
               - plev(isl:iel, jsl:jel, nlay + 1))/nlay

         target_pkap(isl:iel, jsl:jel, k) = &
            target_plev(isl:iel, jsl:jel, k)**kappa

         plog(isl:iel, jsl:jel) = &
            log(target_pkap(isl:iel, jsl:jel, k + 1) &
                /target_pkap(isl:iel, jsl:jel, k))

         target_playkap(isl:iel, jsl:jel, k) = &
            (target_pkap(isl:iel, jsl:jel, k + 1) &
             - target_pkap(isl:iel, jsl:jel, k))/plog(isl:iel, jsl:jel)
      end do

      do k = 1, nlay
         gz(isl:iel, jsl:jel, k + 1) = &
            gz(isl:iel, jsl:jel, k) &
            + dry_heat_capacity*pt(isl:iel, jsl:jel, k)*( &
            pkap(isl:iel, jsl:jel, k) - pkap(isl:iel, jsl:jel, k + 1))
      end do

      do i = isd - 1, ied + 1
         do j = jsd - 1, jed + 1

            ! to remap the winds we interpolate the current and target grids
            ! to cell faces

            src_plev(:) = (7./12.)*(plev(i, j, :) + plev(i, j - 1, :)) &
                          - (1./12.)*(plev(i, j + 1, :) + plev(i, j - 2, :))

            dst_plev(:) = &
               (7./12.)*(target_plev(i, j, :) + target_plev(i, j - 1, :)) &
               - (1./12.)*(target_plev(i, j + 1, :) + target_plev(i, j - 2, :))

            call remap_column(ud(i, j, :), u_remap(i, j, :), src_plev, dst_plev)

            src_plev(:) = (7./12.)*(plev(i, j, :) + plev(i - 1, j, :)) &
                          - (1./12.)*(plev(i + 1, j, :) + plev(i - 2, j, :))

            dst_plev(:) = &
               (7./12.)*(target_plev(i, j, :) + target_plev(i - 1, j, :)) &
               - (1./12.)*(target_plev(i + 1, j, :) + target_plev(i - 2, j, :))

            call remap_column(vd(i, j, :), v_remap(i, j, :), src_plev, dst_plev)

         end do
      end do

      ! we have to finish remapping the winds before we can remap the potential
      ! temperature, since the latter requires the kinetic energy, for which
      ! we must interpolate the remapped winds to the A grid

      do i = isd, ied
         do j = jsd, jed

            ! to remap the potential temperature we instead remap the total
            ! energy density, and then re-calculate the potential temperature
            ! from that, to ensure energy conservation

            ua(:) = (.25*(ud(i, j, :) + ud(i, j + 1, :)) &
                     + .125*(ud(i - 1, j, :) + ud(i - 1, j + 1, :) &
                             + ud(i + 1, j, :) + ud(i + 1, j + 1, :)))

            va(:) = (.25*(vd(i, j, :) + vd(i + 1, j, :)) &
                     + .125*(vd(i, j - 1, :) + vd(i + 1, j - 1, :) &
                             + vd(i, j + 1, :) + vd(i + 1, j + 1, :)))

            energy(:) = dry_heat_capacity*pt(i, j, :)*playkap(i, j, :) &
                        + .5*(gz(i, j, 2:) + gz(i, j, :nlay) &
                              + ua(:)**2 + va(:)**2)

            call remap_column(energy, energy_remap, &
                              plev(i, j, :), target_plev(i, j, :))

            ua(:) = (.25*(u_remap(i, j, :) + u_remap(i, j + 1, :)) &
                     + .125*(u_remap(i - 1, j, :) + u_remap(i - 1, j + 1, :) &
                             + u_remap(i + 1, j, :) + u_remap(i + 1, j + 1, :)))

            va(:) = (.25*(v_remap(i, j, :) + v_remap(i + 1, j, :)) &
                     + .125*(v_remap(i, j - 1, :) + v_remap(i + 1, j - 1, :) &
                             + v_remap(i, j + 1, :) + v_remap(i + 1, j + 1, :)))

            gz_remap(1) = gz(i, j, 1)
            do k = 1, nlay
               pt_remap(i, j, k) = &
                  (energy_remap(k) - gz_remap(k) - .5*(ua(k)**2 + va(k)**2)) &
                  /(dry_heat_capacity*target_playkap(i, j, k) &
                    + .5*dry_heat_capacity*(target_pkap(i, j, k) &
                                            - target_pkap(i, j, k + 1)))

               gz_remap(k + 1) = &
                  gz_remap(k) + dry_heat_capacity*pt_remap(i, j, k) &
                  *(target_pkap(i, j, k) - target_pkap(i, j, k + 1))
            end do

         end do
      end do

      dp(isd:ied, jsd:jed, :) = target_plev(isd:ied, jsd:jed, :nlay) &
                                - target_plev(isd:ied, jsd:jed, 2:)
      ud(isd:ied, jsd:jed, :) = u_remap(isd:ied, jsd:jed, :)
      vd(isd:ied, jsd:jed, :) = v_remap(isd:ied, jsd:jed, :)
      pt(isd:ied, jsd:jed, :) = pt_remap(isd:ied, jsd:jed, :)

   end subroutine perform_vertical_remapping

   subroutine remap_column(src_var, dst_var, src_coord, dst_coord)
      real(rk), intent(in) :: src_var(config%nlay), src_coord(config%nlay + 1)
      real(rk), intent(inout) :: dst_var(config%nlay), &
                                 dst_coord(config%nlay + 1)
      real(rk) :: s1, s2, integral_factor, tmp1, tmp2, tmp3, dsrc_coord
      integer(ik) :: k, m

      dst_var(:) = 0.

      call vertical_ppm(src_var, src_coord)

      do k = 1, config%nlay
         do m = 1, config%nlay
            ! calculate the contribution from source layer k to destination
            ! layer m

            ! Each layer is given a local coordinate s such that s=0 at
            ! the layer bottom and s=1 at the layer top; we perform the
            ! integrations in the s coordinate for convenience but really
            ! we want to be taking averages in our real coordinate space,
            ! i.e. (𝛿p)⁻¹ ∫ (...) dp = (𝛿p)⁻¹ ∫ (...) (dp/ds) ds
            !                        = (𝛿p source)/(𝛿p destination) ∫ (...) ds

            dsrc_coord = src_coord(k + 1) - src_coord(k)

            s1 = (dst_coord(m) - src_coord(k))/dsrc_coord
            s2 = (dst_coord(m + 1) - src_coord(k))/dsrc_coord
            integral_factor = dsrc_coord/(dst_coord(m + 1) - dst_coord(m))

            ! There are six possible configurations
            ! (note that because we work in pressure coordinates "above"
            !  and "below" in these descriptions seem flipped vs. the actual
            !  comparisons carried out, e.g. "above" means lower pressure)

            ! 1. The destination layer is entirely above the source layer
            ! 2. The source layer is entirely above the destination layer
            ! These cases contribute nothing

            if (src_coord(k + 1) >= dst_coord(m) &
                .or. src_coord(k) <= dst_coord(m + 1)) cycle

            ! 3. The destination layer overlaps the source layer from above:
            !      ┌------ top of destination layer
            !     ┌------- top of source layer (s=1)
            !     │└------ bottom of destination layer (s=s1)
            !     └------- bottom of source layer (s=0)
            ! The source layer contributes the integral from s=s1 to s=1

            tmp3 = 2.*s1**2

            if (dst_coord(m + 1) <= src_coord(k + 1) &
                .and. src_coord(k + 1) < dst_coord(m) &
                .and. dst_coord(m) < src_coord(k)) then

               tmp1 = s1 - 1.
               tmp2 = s1 + 1.

               dst_var(m) = dst_var(m) + tmp1/6.*( &
                            3.*(btm_edge(k)*tmp1 &
                                - top_edge(k)*tmp2) &
                            + curvature(k)*(tmp3 - tmp2) &
                            )*integral_factor

            end if

            ! 4. The destination layer overlaps the source layer from below:
            !     ┌------- top of source layer (s=1)
            !     │┌------ top of destination layer (s=s1)
            !     └------ bottom of source layer (s=0)
            !      └------ bottom of destination layer
            ! The source layer contributes the integral from s=0 to s=s1

            if (src_coord(k + 1) < dst_coord(m + 1) &
                .and. dst_coord(m + 1) < src_coord(k) &
                .and. src_coord(k) <= dst_coord(m)) then

               dst_var(m) = dst_var(m) + (s2/2.)*( &
                            btm_edge(k)*(2.-s2) &
                            + s2*(top_edge(k) + curvature(k)*(1.-(2./3.)*s2)) &
                            )*integral_factor

            end if

            ! 5. The destination layer entirely contains the source layer
            !     ┌------- top of destination layer
            !     │┌------ top of source layer (s=1)
            !     │└------ bottom of source layer (s=0)
            !     └------- bottom of destination layer
            ! The source layer contributes its full value (from s=0 to s=1)

            if (dst_coord(m + 1) <= src_coord(k + 1) &
                .and. src_coord(k) <= dst_coord(m)) then

               dst_var(m) = dst_var(m) + src_var(k)*integral_factor

            end if

            ! 6. The source layer entirely contains the destination layer
            !     ┌------- top of source layer (s=1)
            !     │┌------ top of destination layer (s=s2)
            !     │└------ bottom of destination layer (s=s1)
            !     └------- bottom of source layer (s=0)
            ! The source layer contributes the integral from s=s1 to s=s2

            if (src_coord(k + 1) < dst_coord(m + 1) &
                .and. dst_coord(m) < src_coord(k)) then

               tmp1 = s1 + s2

               dst_var(m) = dst_var(m) + (s1 - s2)/6.*( &
                            3.*(btm_edge(k)*(tmp1 - 2.) &
                                - top_edge(k)*tmp1) &
                            + curvature(k)*(tmp3 + tmp1*(2.*s2 - 3.)) &
                            )*integral_factor

            end if

         end do
      end do

   end subroutine remap_column

   subroutine vertical_ppm(var, coord)
      real(rk), intent(in) :: var(config%nlay), coord(config%nlay + 1)
      integer(ik) :: k, nlay
      real(rk) :: alpha, beta, gam(config%nlay), level_var(config%nlay + 1), &
                  pmp1, lac1

      nlay = config%nlay

      ! top layer
      alpha = (coord(nlay) - coord(nlay - 1))/(coord(nlay + 1) - coord(nlay))
      beta = alpha*(.5 + alpha)
      level_var(nlay + 1) = (2.*alpha*(1.+alpha)*var(nlay) + var(nlay - 1))/beta
      gam(nlay) = (1.+alpha*(1.5 + alpha))/beta

      ! interior layers
      do k = nlay, 2, -1
         alpha = (coord(k + 1) - coord(k))/(coord(k) - coord(k - 1))
         beta = 2.+2.*alpha - gam(k)
         level_var(k) = (3.*(var(k) + alpha*var(k - 1)) - level_var(k + 1))/beta
         gam(k - 1) = alpha/beta
      end do

      ! bottom layer
      beta = 1.+alpha*(1.5 + alpha)
      level_var(1) = (2*alpha*(1 + alpha)*var(1) + var(2) - beta*level_var(2)) &
                     /(alpha*(.5 + alpha) - beta*gam(1))

      do k = 2, nlay + 1
         level_var(k) = level_var(k) - gam(k - 1)*level_var(k - 1)
      end do

      ! clamp boundary values
      level_var(nlay) = min(level_var(nlay), max(var(nlay), var(nlay - 1)))
      level_var(nlay) = max(level_var(nlay), min(var(nlay), var(nlay - 1)))
      level_var(2) = min(level_var(2), max(var(2), var(1)))
      level_var(2) = max(level_var(2), min(var(2), var(1)))

      ! monotonicity constraint
      do k = 2, nlay - 1
         if (gam(k + 2)*gam(k) > 0) then
            level_var(k + 1) = min(level_var(k + 1), max(var(k), var(k + 1)))
            level_var(k + 1) = max(level_var(k + 1), min(var(k), var(k + 1)))
         else if (gam(k + 2) > 0) then
            level_var(k + 1) = max(level_var(k + 1), min(var(k), var(k + 1)))
         else
            level_var(k + 1) = min(level_var(k + 1), max(var(k), var(k + 1)))
         end if
      end do

      top_edge(:) = level_var(2:)
      btm_edge(:) = level_var(:nlay)

      ! apply the Huynh constraint
      do k = 2, nlay - 1
         pmp1 = var(k + 1) - 2.*gam(k)
         lac1 = pmp1 + 1.5*gam(k - 1)
         top_edge(k + 1) = &
            min(max(top_edge(k + 1), min(var(k + 1), pmp1, lac1)), &
                max(var(k + 1), pmp1, lac1))

         pmp1 = var(k + 1) + 2.*gam(k + 1)
         lac1 = pmp1 - 1.5*gam(k)
         btm_edge(k + 1) = &
            min(max(btm_edge(k + 1), min(var(k + 1), pmp1, lac1)), &
                max(var(k + 1), pmp1, lac1))
      end do

      curvature(:) = 6.*var - 3.*(top_edge + btm_edge)

   end subroutine vertical_ppm

   subroutine allocate_remapping_arrays()
      integer(ik) :: nlay

      nlay = config%nlay

      ! we have to remap u and v over a slightly larger domain so that we
      ! can calculate the kinetic energy of the remapped winds over the actual
      ! remapping domain
      if (.not. allocated(u_remap)) &
         allocate (u_remap(isd - 1:ied + 1, jsd - 1:jed + 1, nlay))
      if (.not. allocated(v_remap)) &
         allocate (v_remap(isd - 1:ied + 1, jsd - 1:jed + 1, nlay))

      if (.not. allocated(pt_remap)) allocate (pt_remap(isd:ied, jsd:jed, nlay))

      if (.not. allocated(target_plev)) &
         allocate (target_plev(isd - 3:ied + 3, jsd - 3:jed + 3, nlay + 1))
      if (.not. allocated(target_pkap)) &
         allocate (target_pkap(isd - 3:ied + 3, jsd - 3:jed + 3, nlay + 1))
      if (.not. allocated(target_playkap)) &
         allocate (target_playkap(isd - 3:ied + 3, jsd - 3:jed + 3, nlay))

      target_plev(:, :, nlay + 1) = top_pressure
      target_pkap(:, :, nlay + 1) = top_pressure**kappa

      if (.not. allocated(btm_edge)) allocate (btm_edge(nlay))
      if (.not. allocated(top_edge)) allocate (top_edge(nlay))
      if (.not. allocated(curvature)) allocate (curvature(nlay))

   end subroutine allocate_remapping_arrays

end module mod_remap
