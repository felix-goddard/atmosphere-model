module mod_sw_dyn

    use mod_kinds, only: ik, rk
    use mod_log, only: logger => main_logger, log_str
    use mod_config, only: config => main_config
    use mod_tiles, only: is, ie, js, je, isd, ied, jsd, jed
    use mod_fields, only: h, ud, vd

    implicit none
    private

    public :: sw_dynamics_step, allocate_sw_dyn_arrays

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

    real(rk), allocatable :: intermediate(:)

    real(rk), allocatable :: left_edge(:), right_edge(:) ! edge values as calculated by the PPM

    integer, parameter :: PPM_UNCONSTRAINED = 0
    integer, parameter :: PPM_CONSTRAINED = 1
        
contains

    subroutine sw_dynamics_step(dt)
        real(rk), intent(in) :: dt
        integer(ik) :: i, j
        real(rk) :: dx, dy, dtdx, dtdy, dt2, dt2dx, dt2dy

        dx = config % dx
        dy = config % dy
        dtdx = dt / dx
        dtdy = dt / dy
        dt2 = dt / 2.
        dt2dx = dt2 / dx
        dt2dy = dt2 / dy

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
        
        do j=js+3,je-3
            ! Calculate the inner advective step for the x-direction flux
            ! of the h half-step; this uses a simple first-order upwind estimate.
            do i=is+1,ie-1
                if (va(i,j) > 0) then
                    intermediate(i) = h(i,j) - .5 * dt2dy * va(i,j) * (h(i,j) - h(i,j-1))
                else
                    intermediate(i) = h(i,j) - .5 * dt2dy * va(i,j) * (h(i,j+1) - h(i,j))
                end if
            end do

            ! Calculate the outer conservative step for the x-direction flux; this
            ! call to `ppm` sets the arrays `left_edge` and `right_edge`;
            call ppm(intermediate(is+1:ie-1), is+1, ie-1, variant=PPM_UNCONSTRAINED)
            
            ! We restrict the domain by two cells on each side each time we do a
            ! PPM interpolation so this is is+3 to ie-3
            do i=is+3,ie-3
                hc(i,j) = h(i,j) - dt2dx * (uc(i+1,j) * right_edge(i) - uc(i,j) * left_edge(i))
            end do
        end do

        ! Now repeat for the y-direction flux of the h half-step
        do i=is+3,ie-3
            do j=js+1,je-1
                if (ua(i,j) > 0) then
                    intermediate(j) = h(i,j) - .5 * dt2dx * ua(i,j) * (h(i,j) - h(i-1,j))
                else
                    intermediate(j) = h(i,j) - .5 * dt2dx * ua(i,j) * (h(i+1,j) - h(i,j))
                end if
            end do

            call ppm(intermediate(js+1:je-1), js+1, je-1, variant=PPM_UNCONSTRAINED)
            
            do j=js+3,je-3
                hc(i,j) = hc(i,j) - dt2dy * (vc(i,j+1) * right_edge(i) - vc(i,j) * left_edge(i))
            end do
        end do

        ! At this point, ch is valid on is+3 to ie-3 (and same for js/je)
        !                uc/vc are valid on is+1 to ie-1 (and same for js/je)

        ! ===========================================================================
        ! Half-step for calculating C grid u and v

        do concurrent (i=is+3:ie-3, j=js+3:je-3)
            vorticity(i,j) = config % coriolis &
                - (uc(i,j) - uc(i,j-1)) / dy   &
                + (vc(i,j) - vc(i-1,j)) / dx

            kinetic_energy(i,j) = .25 * ( ua(i,j) * (uc(i,j) + uc(i+1,j)) &
                                        + va(i,j) * (vc(i,j) + vc(i,j+1)) )

            energy(i,j) = kinetic_energy(i,j) + config % gravity * hc(i,j)
        end do

        do i=is+6,ie-6
            do j=js+4,je-4
                if (ub(i,j) > 0) then
                    intermediate(j) = vorticity(i,j) &
                        - .5 * dt2dx * ub(i,j) * (vorticity(i,j) - vorticity(i-1,j))
                else
                    intermediate(j) = vorticity(i,j) &
                        - .5 * dt2dx * ub(i,j) * (vorticity(i+1,j) - vorticity(i,j))
                end if
            end do

            call ppm(intermediate(js+4:je-4), js+4, je-4, variant=PPM_UNCONSTRAINED)
            
            do j=js+6,je-6
                uc(i,j) = uc(i,j)                         &
                    + dt2 * vd(i,j) * right_edge(j)       &
                    - dt2dx * (energy(i,j) - energy(i-1,j))
            end do
        end do

        do j=js+6,je-6
            do i=is+4,ie-4
                if (vb(i,j) > 0) then
                    intermediate(i) = vorticity(i,j) &
                        - .5 * dt2dy * vb(i,j) * (vorticity(i,j) - vorticity(i,j-1))
                else
                    intermediate(i) = vorticity(i,j) &
                        - .5 * dt2dy * vb(i,j) * (vorticity(i,j+1) - vorticity(i,j))
                end if
            end do

            call ppm(intermediate(is+4:ie-4), is+4, ie-4, variant=PPM_UNCONSTRAINED)
            
            do i=is+6,ie-6
                vc(i,j) = vc(i,j)                         &
                    - dt2 * ud(i,j) * right_edge(i)       &
                    - dt2dy * (energy(i,j) - energy(i,j-1))
            end do
        end do

        ! Update the other grid winds based on the new C grid winds; since uc and vc
        ! get updated on is/js+6 to ie/je-6 we can only calculate these on is+7 to ie-7
        do concurrent (i=is+7:ie-7, j=js+7:je-7)
            ua(i,j) = .5 * (uc(i,j) + uc(i+1,j))
            va(i,j) = .5 * (vc(i,j) + vc(i,j+1))

            ub(i,j) = .5 * (uc(i,j) + uc(i,j-1))
            vb(i,j) = .5 * (vc(i,j) + vc(i-1,j))
        end do

        ! ===========================================================================
        ! Prognostic step for h
        
        do j=js+9,je-9
            do i=is+7,ie-7
                if (va(i,j) > 0) then
                    intermediate(i) = h(i,j) - dt2dy * va(i,j) * (h(i,j) - h(i,j-1))
                else
                    intermediate(i) = h(i,j) - dt2dy * va(i,j) * (h(i,j+1) - h(i,j))
                end if
            end do

            call ppm(intermediate(is+7:ie-7), is+7, ie-7, variant=PPM_CONSTRAINED)
            
            do i=is+9,ie-9
                hc(i,j) = h(i,j) - dt2dx * (                     &
                    uc(i+1,j) * (left_edge(i+1) + right_edge(i)) &
                    - uc(i,j) * (left_edge(i) + right_edge(i-1)) )
            end do
        end do

        do i=is+9,ie-9
            do j=js+7,je-7
                if (ua(i,j) > 0) then
                    intermediate(j) = h(i,j) - dt2dx * ua(i,j) * (h(i,j) - h(i-1,j))
                else
                    intermediate(j) = h(i,j) - dt2dx * ua(i,j) * (h(i+1,j) - h(i,j))
                end if
            end do

            call ppm(intermediate(js+7:je-7), js+7, je-7, variant=PPM_CONSTRAINED)
            
            do j=js+9,je-9
                hc(i,j) = hc(i,j) - dt2dy * (                    &
                    vc(i,j+1) * (left_edge(j+1) + right_edge(j)) &
                    - vc(i,j) * (left_edge(j) + right_edge(j-1)) )
            end do
        end do

        h(is+9:ie-9, js+9:je-9) = hc(is+9:ie-9, js+9:je-9)

        ! ===========================================================================
        ! Prognostic step for u and v

        do concurrent (i=is+7:ie-7, j=js+7:je-7)
            vorticity(i,j) = config % coriolis &
                - (ud(i,j+1) - ud(i,j)) / dy   &
                + (vd(i+1,j) - vd(i,j)) / dx

            kinetic_energy(i,j) = .25 * ( ub(i,j) * (ud(i,j) + ud(i-1,j))   &
                                        + vb(i,j) * (vd(i,j) + vd(i,j-1)) )
        end do

        do concurrent (i=is+10:ie-10, j=js+10:je-10)
            energy(i,j) = kinetic_energy(i,j) + config % gravity * ( &
                h(i,j) + h(i-1,j) + h(i,j-1) + h(i-1,j-1) ) / 4.
        end do

        do i=is+11,ie-11
            do j=js+8,je-8
                if (ua(i,j) > 0) then
                    intermediate(j) = vorticity(i,j) &
                        - dt2dx * ua(i,j) * (vorticity(i,j) - vorticity(i-1,j))
                else
                    intermediate(j) = vorticity(i,j) &
                        - dt2dx * ua(i,j) * (vorticity(i+1,j) - vorticity(i,j))
                end if
            end do

            call ppm(intermediate(js+8:je-8), js+8, je-8, variant=PPM_CONSTRAINED)
            
            do j=js+11,je-11
                ud(i,j) = ud(i,j)                                      &
                    + dt2 * vc(i,j) * (left_edge(j) + right_edge(j-1)) &
                    - dtdx * (energy(i+1,j) - energy(i,j))
            end do
        end do

        do j=js+11,je-11
            do i=is+8,ie-8
                if (va(i,j) > 0) then
                    intermediate(i) = vorticity(i,j) &
                        - dt2dy * va(i,j) * (vorticity(i,j) - vorticity(i,j-1))
                else
                    intermediate(i) = vorticity(i,j) &
                        - dt2dy * va(i,j) * (vorticity(i,j+1) - vorticity(i,j))
                end if
            end do

            call ppm(intermediate(is+8:ie-8), is+8, ie-8, variant=PPM_CONSTRAINED)
            
            do i=is+11,ie-11
                vd(i,j) = vd(i,j)                                      &
                    - dt2 * uc(i,j) * (left_edge(i) + right_edge(i-1)) &
                    - dtdy * (energy(i,j+1) - energy(i,j))
            end do
        end do
        
    end subroutine sw_dynamics_step

    subroutine ppm(q, isl, iel, variant)
        real(rk), intent(in) :: q(isl:iel)
        integer(ik), intent(in) :: isl, iel
        integer :: variant
        real(rk) :: dqi, dqi_min, dqi_max
        real(rk) :: dqi_mono(min(is,js):max(ie,je))
        integer(ik) :: i

        if (variant == PPM_UNCONSTRAINED .or. variant == PPM_CONSTRAINED) then

            do i = isl+1, iel-1
                dqi = (q(i+1) - q(i-1)) / 4
                dqi_max = max(q(i+1), q(i), q(i-1)) - q(i)
                dqi_min = q(i) - min(q(i+1), q(i), q(i-1))
                dqi_mono(i) = sign(min(abs(dqi_min), dqi_min, dqi_max), dqi)
            end do

            do i = isl+1, iel-1
                if (i > isl+1) &
                    left_edge(i) = (q(i) + q(i-1)) / 2. + (dqi_mono(i-1) - dqi_mono(i)) / 3.
                
                if (i < iel-1) &
                    right_edge(i) = (q(i+1) + q(i)) / 2. + (dqi_mono(i) - dqi_mono(i+1)) / 3.

                if (variant == PPM_CONSTRAINED) then
                    if (i > isl+1) &
                        left_edge(i) = q(i) - sign( &
                            min(2*abs(dqi_mono(i)), abs(left_edge(i) - q(i))), dqi_mono(i))

                    if (i < iel-1) &
                        right_edge(i) = q(i) + sign( &
                            min(2*abs(dqi_mono(i)), abs(right_edge(i) - q(i))), dqi_mono(i))
                end if
            end do

        else
            write (log_str, '(a,i4)') 'Unknown PPM variant:', variant
            call logger % fatal('ppm', log_str)
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

        if (.not. allocated(intermediate)) allocate(intermediate(min(is,js):max(ie,je)))

        if (.not. allocated(left_edge)) allocate(left_edge(min(is,js):max(ie,je)))
        if (.not. allocated(right_edge)) allocate(right_edge(min(is,js):max(ie,je)))
        
    end subroutine allocate_sw_dyn_arrays
    
end module mod_sw_dyn