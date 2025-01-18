module mod_fields

    use mod_kinds, only: ik, rk
    use mod_config, only: config => main_config
    use mod_tiles, only: tile_indices

    implicit none
    private

    integer(ik), parameter :: halo_width = 12
    integer(ik) :: is, ie, js, je ! data bounds (including halo)
    integer(ik) :: isd, ied, jsd, jed ! domain bounds (not including halo)

    public :: init_prognostic_fields, &
              is, ie, js, je, isd, ied, jsd, jed, halo_width

    real(rk), allocatable :: h(:,:)  ! prognostic height
    real(rk), allocatable :: ud(:,:) ! prognostic u wind (on D grid)
    real(rk), allocatable :: vd(:,:) ! prognostic v wind (on D grid)

    public :: h, ud, vd
    
contains

    subroutine init_prognostic_fields()
        integer(ik) :: indices(4), upper(2), lower(2)
        integer(ik) :: i, j

        real(rk), parameter :: height = 10
        real(rk), parameter :: ic = .5, jc = .5
        real(rk), parameter :: decay = 500
        
        indices = tile_indices([config % nx, config % ny])
        lower = indices([1, 3])
        upper = indices([2, 4])

        isd = lower(1)
        ied = upper(1)
        jsd = lower(2)
        jed = upper(2)

        is = isd - halo_width
        ie = ied + halo_width
        js = jsd - halo_width
        je = jed + halo_width

        if (.not. allocated(h))  allocate(h (is:ie, js:je))
        if (.not. allocated(ud)) allocate(ud(is:ie, js:je))
        if (.not. allocated(vd)) allocate(vd(is:ie, js:je))

        ! initialize a gaussian blob in the center
        h(:,:) = 0.
        do concurrent (i=isd:ied, j=jsd:jed)
            h(i,j) = 10000. + height * exp( &
                -decay * ( (real(i-1)/(config%nx-1) - ic)**2 &
                         + (real(j-1)/(config%ny-1) - jc)**2))
        end do

        ud(:,:) = 0.
        vd(:,:) = 0.

    end subroutine init_prognostic_fields
    
end module mod_fields