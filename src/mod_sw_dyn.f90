module mod_sw_dyn

    use ieee_arithmetic, only: ieee_is_finite
    use mod_kinds, only: ik, rk
    use mod_log, only: logger => main_logger, log_str
    use mod_config, only: config => main_config
    use mod_tiles, only: is, ie, js, je, isd, ied, jsd, jed
    use mod_fields, only: h, ud, vd
    use mod_util, only: abort_now

    implicit none
    private

    public :: sw_dynamics_step, allocate_sw_dyn_arrays, is_stable

    real(rk), allocatable :: hc(:,:) ! height for the C grid half-step

    real(rk), allocatable :: vorticity(:,:)
    real(rk), allocatable :: kinetic_energy(:,:)
    real(rk), allocatable :: energy(:,:)

    real(rk), allocatable :: ua(:,:) ! u wind on A grid
    real(rk), allocatable :: va(:,:) ! v wind on A grid
    real(rk), allocatable :: ub(:,:) ! u wind on B grid
    real(rk), allocatable :: vb(:,:) ! v wind on B grid
    real(rk), allocatable :: uc(:,:) ! u wind on C grid
    real(rk), allocatable :: vc(:,:) ! v wind on C grid

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
    ! i.e. the A grid winds (ua and va, also h) are located at mass
    ! points (cell centers); the B grid winds are located at cell
    ! corners; the C grid winds are perpendicular to cell edges,
    ! and the prognostic D grid winds are parallel to cell edges.

    real(rk), allocatable :: tmp(:,:)

    real(rk), allocatable :: fx(:), fy(:), courant(:)
    real(rk), allocatable :: edge_L(:), edge_R(:) ! values calculated by the PPM

    integer(ik), parameter :: PPM_UNCONSTRAINED = 0
    integer(ik), parameter :: PPM_CONSTRAINED = 1

contains

    function is_stable()
        logical :: is_stable

        is_stable = (all(ieee_is_finite(h(isd:ied, jsd:jed))) &
            .and. all(ieee_is_finite(ud(isd:ied, jsd:jed)))   &
            .and. all(ieee_is_finite(vd(isd:ied, jsd:jed))))

    end function is_stable

    subroutine sw_dynamics_step(dt)
        ! Solve the shallow water equations following the method of Lin & Rood (1997)
        ! (i.e. on the CD grid using the finite-volume multidimensional transport scheme
        ! of Lin & Rood (1996)).

        real(rk), intent(in) :: dt
        integer(ik) :: i, j
        real(rk) :: dx, dy, dtdx, dtdy, dt2, dt2dx, dt2dy, rA
        real(rk) :: denom_x(is:ie, js:je), denom_y(is:ie, js:je)

        dx = config % dx
        dy = config % dy
        dtdx = dt / dx
        dtdy = dt / dy
        dt2 = dt / 2.
        dt2dx = dt2 / dx
        dt2dy = dt2 / dy

        rA = 1. / (dx * dy)

        ! Calculate the winds
        do concurrent (i=is+1:ie-1, j=js+1:je-1)
            ub(i,j) = .5 * (ud(i,j) + ud(i-1,j))
            vb(i,j) = .5 * (vd(i,j) + vd(i,j-1))

            uc(i,j) = .5 * (ub(i,j) + .5 * (ud(i,j+1) + ud(i-1,j+1)))
            vc(i,j) = .5 * (vb(i,j) + .5 * (vd(i+1,j) + vd(i+1,j-1)))

            ua(i,j) = .5 * (uc(i,j) + .25 * (ud(i+1,j) + ud(i,j) + ud(i+1,j+1) + ud(i,j+1)))
            va(i,j) = .5 * (vc(i,j) + .25 * (vd(i,j+1) + vd(i,j) + vd(i+1,j+1) + vd(i+1,j)))
        end do

        ! ===========================================================================
        ! Half-step for calculating C grid h

        do concurrent (i=is+2:ie-2, j=js+2:je-2)
            denom_x(i,j) = 1. / (1. - dt2dx * (uc(i+1,j) - uc(i,j)))
            denom_y(i,j) = 1. / (1. - dt2dy * (vc(i,j+1) - vc(i,j)))
        end do

        ! calculate the inner step (in the y-direction) for the outer x-direction step
        do i = is+2, ie-2
            courant(js+1:je-1) = dt2dy * vc(i,js+1:je-1)
            call ppm_flux(fy(js+1:je-1), h(i,js+1:je-1), courant(js+1:je-1), &
                js+1, je-1, variant=PPM_UNCONSTRAINED)

            do j = js+4, je-4
                tmp(i,j) = .5 * (h(i,j) + denom_y(i,j) * (h(i,j) - rA * (fy(j+1) - fy(j))))
            end do
        end do

        ! calculate the outer step in the x-direction
        do j = js+5, je-5
            courant(is+2:ie-2) = dt2dx * uc(is+2:ie-2,j)
            call ppm_flux(fx(is+2:ie-2), tmp(is+2:ie-2,j), courant(is+2:ie-2), &
                is+2, ie-2, variant=PPM_UNCONSTRAINED)

            do i = is+5, ie-5
                hc(i,j) = h(i,j) - (fx(i+1) * courant(i+1) - fx(i) * courant(i))
            end do
        end do

        ! calculate the inner step (in the x-direction) for the outer y-direction step
        do j = js+2, je-2
            courant(is+2:ie-2) = dt2dx * uc(is+2:ie-2,j)
            call ppm_flux(fx(is+2:ie-2), h(is+2:ie-2,j), courant(is+2:ie-2), &
                is+2, ie-2, variant=PPM_UNCONSTRAINED)

            do i = is+4, ie-4
                tmp(i,j) = .5 * (h(i,j) + denom_x(i,j) * (h(i,j) - rA * (fx(i+1) - fx(i))))
            end do
        end do

        ! calculate the outer step in the y-direction
        do i = is+5, ie-5
            courant(js+2:je-2) = dt2dy * vc(i,js+2:je-2)

            call ppm_flux(fy(js+2:je-2), tmp(i,js+2:je-2), courant(js+2:je-2), &
                js+2, je-2, variant=PPM_UNCONSTRAINED)

            do j = js+5, je-5
                hc(i,j) = hc(i,j) - (fy(j+1) * courant(j+1) - fy(j) * courant(j))
            end do
        end do

        ! ===========================================================================
        ! Energy and vorticity

        do i = is+2, ie-2
            courant(js+1:je-1) = dt2dy * va(i,js+1:je-1)
            call ppm_flux(fy(js+1:je-1), vc(i,js+1:je-1), courant(js+1:je-1), &
                js+1, je-1, variant=PPM_UNCONSTRAINED)
        end do

        do j = js+2, je-2
            courant(is+1:ie-1) = dt2dx * ua(is+1:ie-1,j)
            call ppm_flux(fx(is+1:ie-1), uc(is+1:ie-1,j), courant(is+1:ie-1), &
                is+1, ie-1, variant=PPM_UNCONSTRAINED)
        end do

        do concurrent (i=is+3:ie-3, j=js+3:je-3)
            vorticity(i,j) = config % coriolis &
                - (uc(i,j) - uc(i,j-1)) / dy   &
                + (vc(i,j) - vc(i-1,j)) / dx

            kinetic_energy(i,j) = .5 * (ua(i,j) * fx(i+1) + va(i,j) * fy(j+1))

            energy(i,j) = kinetic_energy(i,j) + config % gravity * hc(i,j)
        end do

        ! ===========================================================================
        ! Half-step for calculating C grid u

        ! calculate the inner step (in the x-direction) for the outer y-direction step
        do j = js+2, je-2
            courant(is+2:ie-2) = dt2dx * ub(is+2:ie-2,j)
            call ppm_flux(fx(is+2:ie-2), vorticity(is+2:ie-2,j), courant(is+2:ie-2), &
                is+2, ie-2, variant=PPM_UNCONSTRAINED)

            do i = is+4, ie-4
                tmp(i,j) = .5 * (vorticity(i,j) + &
                    denom_x(i,j) * (vorticity(i,j) - rA * (fx(i+1) - fx(i))))
            end do
        end do

        do i = is+6, ie-6
            courant(js+2:je-2) = dt2dy * vd(i,js+2:je-2)
            call ppm_flux(fy(js+2:je-2), tmp(i,js+2:je-2), courant(js+2:je-2), &
                js+2, je-2, variant=PPM_UNCONSTRAINED)

            do j = js+6, je-6
                uc(i,j) = uc(i,j)                         &
                    + dt2 * vd(i,j) * fy(j+1)             &
                    - dt2dx * (energy(i,j) - energy(i-1,j))
            end do
        end do

        ! ===========================================================================
        ! Half-step for calculating C grid v

        ! calculate the inner step (in the y-direction) for the outer x-direction step
        do i = is+2, ie-2
            courant(js+2:je-2) = dt2dy * vb(i,js+2:je-2)
            call ppm_flux(fy(js+2:je-2), vorticity(i,js+2:je-2), courant(js+2:je-2), &
                js+2, je-2, variant=PPM_UNCONSTRAINED)

            do j = js+4, je-4
                tmp(i,j) = .5 * (vorticity(i,j) + &
                    denom_y(i,j) * (vorticity(i,j) - rA * (fy(j+1) - fy(j))))
            end do
        end do

        do j = js+6, je-6
            courant(is+2:ie-2) = dt2dx * ud(is+2:ie-2,j)
            call ppm_flux(fx(is+2:ie-2), tmp(is+2:ie-2,j), courant(is+2:ie-2), &
                is+2, ie-2, variant=PPM_UNCONSTRAINED)

            do i = is+6, ie-6
                vc(i,j) = vc(i,j)                         &
                    - dt2 * ud(i,j) * fx(i+1)             &
                    - dt2dy * (energy(i,j) - energy(i,j-1))
            end do
        end do

        ! Update the other winds based on the new C grid winds; since uc and vc
        ! get updated on is/js+6 to ie/je-6 we can only calculate these on is+7 to ie-7
        do concurrent (i=is+7:ie-7, j=js+7:je-7)
            ua(i,j) = .5 * (uc(i,j) + uc(i+1,j))
            va(i,j) = .5 * (vc(i,j) + vc(i,j+1))

            ub(i,j) = .5 * (uc(i,j) + uc(i,j-1))
            vb(i,j) = .5 * (vc(i,j) + vc(i-1,j))
        end do

        do concurrent (i=is+2:ie-2, j=js+2:je-2)
            denom_x(i,j) = 1. / (1. - dtdx * (uc(i+1,j) - uc(i,j)))
            denom_y(i,j) = 1. / (1. - dtdy * (vc(i,j+1) - vc(i,j)))
        end do

        ! ===========================================================================
        ! Prognostic step for calculating h

        ! calculate the inner step (in the y-direction) for the outer x-direction step
        do i = is+2, ie-2
            courant(js+1:je-1) = dtdy * va(i,js+1:je-1)
            call ppm_flux(fy(js+1:je-1), h(i,js+1:je-1), courant(js+1:je-1), &
                js+1, je-1, variant=PPM_CONSTRAINED)

            do j = js+4, je-4
                tmp(i,j) = .5 * (h(i,j) + &
                    denom_y(i,j) * (h(i,j) - rA * (fy(j+1) - fy(j))))
            end do
        end do

        ! calculate the outer step in the x-direction
        do j = js+5, je-5
            courant(is+2:ie-2) = dtdx * uc(is+2:ie-2,j)
            call ppm_flux(fx(is+2:ie-2), tmp(is+2:ie-2,j), courant(is+2:ie-2), &
                is+2, ie-2, variant=PPM_CONSTRAINED)

            do i = is+5, ie-5
                hc(i,j) = h(i,j) - (fx(i+1) * courant(i+1) - fx(i) * courant(i))
            end do
        end do

        ! calculate the inner step (in the x-direction) for the outer y-direction step
        do j = js+2, je-2
            courant(is+2:ie-2) = dtdx * ua(is+2:ie-2,j)
            call ppm_flux(fx(is+2:ie-2), h(is+2:ie-2,j), courant(is+2:ie-2), &
                is+2, ie-2, variant=PPM_CONSTRAINED)

            do i = is+4, ie-4
                tmp(i,j) = .5 * (h(i,j) + &
                    denom_x(i,j) * (h(i,j) - rA * (fx(i+1) - fx(i))))
            end do
        end do

        ! calculate the outer step in the y-direction
        do i = is+5, ie-5
            courant(js+2:je-2) = dtdy * vc(i,js+2:je-2)
            call ppm_flux(fy(js+2:je-2), tmp(i,js+2:je-2), courant(js+2:je-2), &
                js+2, je-2, variant=PPM_CONSTRAINED)

            do j = js+5, je-5
                hc(i,j) = hc(i,j) - (fy(j+1) * courant(j+1) - fy(j) * courant(j))
            end do
        end do

        h(is+9:ie-9, js+9:je-9) = hc(is+9:ie-9, js+9:je-9)

        ! ===========================================================================
        ! Energy and vorticity

        do i = is+2, ie-2
            courant(js+1:je-1) = dtdy * vb(i,js+1:je-1)
            call ppm_flux(kinetic_energy(i,js+1:je-1), vd(i,js+1:je-1), courant(js+1:je-1), &
                js+1, je-1, variant=PPM_CONSTRAINED)
        end do

        do j = js+2, je-2
            courant(is+1:ie-1) = dtdx * ub(is+1:ie-1,j)
            call ppm_flux(fx(is+1:ie-1), ud(is+1:ie-1,j), courant(is+1:ie-1), &
                is+1, ie-1, variant=PPM_CONSTRAINED)

            kinetic_energy(is+1:ie-1,j) = .5 * ( &
                ub(is+1:ie-1,j) * fx(is+1:ie-1) &
                + vb(is+1:ie-1,j) * kinetic_energy(is+1:ie-1,j))
        end do

        do concurrent (i=is+3:ie-3, j=js+3:je-3)
            vorticity(i,j) = config % coriolis &
                - (ud(i,j+1) - ud(i,j)) / dy   &
                + (vd(i+1,j) - vd(i,j)) / dx

            energy(i,j) = kinetic_energy(i,j) + config % gravity * ( &
                h(i,j) + h(i-1,j) + h(i,j-1) + h(i-1,j-1) ) / 4.
        end do

        ! ===========================================================================
        ! Prognostic for calculating u

        ! calculate the inner step (in the x-direction) for the outer y-direction step
        do j = js+2, je-2
            courant(is+2:ie-2) = dtdx * ua(is+2:ie-2,j)
            call ppm_flux(fx(is+2:ie-2), vorticity(is+2:ie-2,j), courant(is+2:ie-2), &
                is+2, ie-2, variant=PPM_CONSTRAINED)

            do i = is+4, ie-4
                tmp(i,j) = .5 * (vorticity(i,j) + &
                    denom_x(i,j) * (vorticity(i,j) - rA * (fx(i+1) - fx(i))))
            end do
        end do

        do i = is+6, ie-6
            courant(js+2:je-2) = dtdy * vc(i,js+2:je-2)
            call ppm_flux(fy(js+2:je-2), tmp(i,js+2:je-2), courant(js+2:je-2), &
                is+2, ie-2, variant=PPM_CONSTRAINED)

            do j = js+6, je-6
                ud(i,j) = ud(i,j)                        &
                    + dt * vc(i,j) * fy(j)               &
                    - dtdx * (energy(i+1,j) - energy(i,j))
            end do
        end do

        ! ===========================================================================
        ! Prognostic for calculating v

        ! calculate the inner step (in the y-direction) for the outer x-direction step
        do i = is+2, ie-2
            courant(js+2:je-2) = dtdy * va(i,js+2:je-2)
            call ppm_flux(fy(js+2:je-2), vorticity(i,js+2:je-2), courant(js+2:je-2), &
                js+2, je-2, variant=PPM_CONSTRAINED)

            do j = js+4, je-4
                tmp(i,j) = .5 * (vorticity(i,j) + &
                    denom_y(i,j) * (vorticity(i,j) - rA * (fy(j+1) - fy(j))))
            end do
        end do

        do j = js+6, je-6
            courant(is+2:ie-2) = dtdx * uc(is+2:ie-2,j)
            call ppm_flux(fx(is+2:ie-2), tmp(is+2:ie-2,j), courant(is+2:ie-2), &
                is+2, ie-2, variant=PPM_CONSTRAINED)

            do i = is+6, ie-6
                vd(i,j) = vd(i,j)                        &
                    - dt * uc(i,j) * fx(i)               &
                    - dtdy * (energy(i,j+1) - energy(i,j))
            end do
        end do

    end subroutine sw_dynamics_step

    subroutine ppm_flux(f, q, c, isl, iel, variant)
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

        call ppm(q, isl, iel, variant)

        ! ppm defines the edges values over isl+2:iel-2 so we can calculate
        ! fluxes over isl+3:iel-2
        do i = isl+3, iel-2
            s = abs(c(i))
            k1 = (1. - 2. * s) * (1. - s)
            k2 = s * (s - 1.)

            if (c(i) > 0.) then
                f(i) = q(i-1) + k1 * edge_R(i-1) + k2 * edge_L(i-1)
            else
                f(i) = q(i) + k1 * edge_L(i) + k2 * edge_R(i)
            end if
        end do

    end subroutine ppm_flux

    subroutine ppm(q, isl, iel, variant)
        ! Calculate, by a piecewise parabolic reconstruction, the left and right cell
        ! edge anomaly values given the cell averages q(isl:iel).
        ! The `variant` argument selects between different reconstruction algorithms,
        ! currently the options are:
        !  - PPM_UNCONSTRAINED: implements then simplified method of Lin (2004) without
        !                       any constraints on the edge values.
        !  - PPM_CONSTRAINED: also implements the simplified method of Lin (2004) with
        !                     the addition of a monotonicity constraint.

        real(rk), intent(in) :: q(isl:iel)
        integer(ik), intent(in) :: isl, iel
        integer(ik), intent(in) :: variant
        real(rk) :: dqi_min, dqi_max
        real(rk) :: dqi(isl+1:iel-1)
        real(rk) :: dqi_mono(isl+1:iel-1)
        integer(ik) :: i

        if (variant == PPM_UNCONSTRAINED .or. variant == PPM_CONSTRAINED) then

            do i = isl+1, iel-1
                dqi(i) = (q(i+1) - q(i-1)) / 4.
                dqi_max = max(q(i+1), q(i), q(i-1)) - q(i)
                dqi_min = q(i) - min(q(i+1), q(i), q(i-1))
                dqi_mono(i) = sign(min(abs(dqi(i)), dqi_min, dqi_max), dqi(i))
            end do

            do i = isl+2, iel-2
                edge_L(i) = (q(i-1) - q(i)) / 2. + (dqi_mono(i-1) - dqi_mono(i)) / 3.
                edge_R(i) = (q(i+1) - q(i)) / 2. + (dqi_mono(i) - dqi_mono(i+1)) / 3.

                if (variant == PPM_CONSTRAINED) then
                    edge_L(i) = -sign( &
                        min(2.*abs(dqi(i)), abs(edge_L(i))), dqi(i))

                    edge_R(i) = +sign( &
                        min(2.*abs(dqi(i)), abs(edge_R(i))), dqi(i))
                end if
            end do

        else
            write (log_str, '(a,i4)') 'Unknown PPM variant:', variant
            call logger % fatal('ppm', log_str)
            call abort_now()
        end if

    end subroutine ppm

    subroutine allocate_sw_dyn_arrays()

        if (.not. allocated(vorticity))      allocate(vorticity     (is:ie, js:je))
        if (.not. allocated(kinetic_energy)) allocate(kinetic_energy(is:ie, js:je))
        if (.not. allocated(energy))         allocate(energy        (is:ie, js:je))

        if (.not. allocated(ua)) allocate(ua(is:ie, js:je))
        if (.not. allocated(va)) allocate(va(is:ie, js:je))
        if (.not. allocated(ub)) allocate(ub(is:ie, js:je))
        if (.not. allocated(vb)) allocate(vb(is:ie, js:je))
        if (.not. allocated(uc)) allocate(uc(is:ie, js:je))
        if (.not. allocated(vc)) allocate(vc(is:ie, js:je))

        if (.not. allocated(hc)) allocate(hc(is:ie, js:je))

        if (.not. allocated(tmp)) allocate(tmp(is:ie, js:je))

        if (.not. allocated(fx))      allocate(fx(min(is,js):max(ie,je)))
        if (.not. allocated(fy))      allocate(fy(min(is,js):max(ie,je)))
        if (.not. allocated(courant)) allocate(courant(min(is,js):max(ie,je)))

        if (.not. allocated(edge_L)) allocate(edge_L(min(is,js):max(ie,je)))
        if (.not. allocated(edge_R)) allocate(edge_R(min(is,js):max(ie,je)))

    end subroutine allocate_sw_dyn_arrays

end module mod_sw_dyn
