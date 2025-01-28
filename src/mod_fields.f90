module mod_fields

    use mod_kinds, only: ik, rk
    use mod_log, only: logger => main_logger, log_str
    use mod_config, only: config => main_config
    use mod_netcdf, only: netcdf_file, open_netcdf
    use mod_tiles, only: is, ie, js, je, isd, ied, jsd, jed
    use mod_util, only: abort_now

    implicit none
    private

    public :: init_prognostic_fields

    real(rk), allocatable :: h(:,:)  ! prognostic height
    real(rk), allocatable :: ud(:,:) ! prognostic u wind (on D grid)
    real(rk), allocatable :: vd(:,:) ! prognostic v wind (on D grid)

    public :: h, ud, vd

contains

    subroutine init_prognostic_fields()
        integer(ik) :: i, j
        character(len=:), allocatable :: names(:)
        real(rk), allocatable :: values(:,:)
        type(netcdf_file) :: initial_nc

        if (.not. allocated(h))  allocate(h (is:ie, js:je))
        if (.not. allocated(ud)) allocate(ud(is:ie, js:je))
        if (.not. allocated(vd)) allocate(vd(is:ie, js:je))

        initial_nc = open_netcdf(config % initial_filename)

        ! Check we have all the expected coordinates
        names = ['xc', 'yc', 'xf', 'yf']
        axis_loop: do i = 1, size(names)
            do j = 1, size(initial_nc % axes)
                if (trim(initial_nc % axes(j) % name) == trim(names(i))) &
                    cycle axis_loop
            end do
            write (log_str, '(a)') 'Could not find axis `' // trim(names(i)) &
                // '` in initial condition netCDF.'
            call logger % fatal('init_prognostic_fields', log_str)
            call abort_now()
        end do axis_loop

        ! Check we have all the expected variables
        names = ['h', 'u', 'v']
        var_loop: do i = 1, size(names)
            do j = 1, size(initial_nc % vars)
                if (trim(initial_nc % vars(j) % name) == trim(names(i))) &
                    cycle var_loop
            end do
            write (log_str, '(a)') 'Could not find variable `' // trim(names(i)) &
                // '` in initial condition netCDF.'
            call logger % fatal('init_prognostic_fields', log_str)
            call abort_now()
        end do var_loop

        allocate(values(1:config % nx, 1:config % ny))

        call initial_nc % read_var('h', values(1:config % nx, 1:config % ny))
        h(isd:ied, jsd:jed) = values(isd:ied, jsd:jed)

        call initial_nc % read_var('u', values(1:config % nx, 1:config % ny))
        ud(isd:ied, jsd:jed) = values(isd:ied, jsd:jed)

        call initial_nc % read_var('v', values(1:config % nx, 1:config % ny))
        vd(isd:ied, jsd:jed) = values(isd:ied, jsd:jed)

        deallocate(values)
        call initial_nc % close()

    end subroutine init_prognostic_fields

end module mod_fields
