module mod_netcdf

    use iso_fortran_env, only: stderr => error_unit
    use mod_kinds, only: ik, rk
    use mod_log, only: logger => main_logger, log_str
    use mod_config, only: config => main_config
    use mod_util, only: abort_now
    use netcdf, only: &
        nf90_create, nf90_def_dim, nf90_def_var, nf90_float, nf90_put_var, &
        nf90_enddef, nf90_redef, nf90_close, nf90_unlimited, nf90_clobber, &
        nf90_strerror, nf90_noerr, nf90_sync

    implicit none

    private
    public :: netcdf_file, create_netcdf

    integer(ik), parameter :: max_var_name_length = 10

    type :: netcdf_file
        private
        integer(ik) :: ncid
        integer(ik) :: tdim = -1, tvar = -1, tidx = -1 ! time axis
        type(netcdf_axis), allocatable :: axes(:)
        type(netcdf_var), allocatable :: vars(:)
    contains
        procedure :: close => close_netcdf
        procedure :: create_time_axis
        procedure :: advance_time
        procedure :: create_axis
        procedure :: create_var
        procedure :: write_var
    end type

    type :: netcdf_axis
        character(len=:), allocatable :: name
        integer(ik) :: dimid, varid
    end type

    type :: netcdf_var
        character(len=:), allocatable :: name
        integer(ik) :: varid
        integer(ik), allocatable :: dimids(:)
    end type

contains

    function create_netcdf(filename) result(out)
        character(len=*), intent(in) :: filename
        type(netcdf_file) :: out
        integer(ik) :: status, ncid

        status = nf90_create(path=filename, cmode=nf90_clobber, ncid=ncid)
        call handle_netcdf_error(status, 'create_netcdf')

        status = nf90_enddef(ncid)
        call handle_netcdf_error(status, 'create_netcdf')

        out = netcdf_file(ncid)

        out % axes = [netcdf_axis ::]
        out % vars = [netcdf_var ::]

    end function

    subroutine close_netcdf(self)
        class(netcdf_file), intent(in) :: self
        integer(ik) :: status

        status = nf90_close(self % ncid)
        call handle_netcdf_error(status, 'close_netcdf')

    end subroutine close_netcdf

    subroutine create_time_axis(self)
        class(netcdf_file), intent(inout) :: self
        integer(ik) :: status

        status = nf90_redef(self % ncid)
        call handle_netcdf_error(status, 'create_time_axis')

        status = nf90_def_dim(self % ncid, 't', nf90_unlimited, self % tdim)
        call handle_netcdf_error(status, 'create_time_axis')

        status = nf90_def_var(self % ncid, 't', nf90_float, self % tdim, self % tvar)
        call handle_netcdf_error(status, 'create_time_axis')

        self % tidx = 0

        status = nf90_enddef(self % ncid)
        call handle_netcdf_error(status, 'create_time_axis')

    end subroutine create_time_axis

    subroutine advance_time(self, time)
        class(netcdf_file), intent(inout) :: self
        real(rk), intent(in) :: time
        integer(ik) :: status

        self % tidx = self % tidx + 1

        status = nf90_put_var(self % ncid, self % tvar, time, start=[self % tidx] )
        call handle_netcdf_error(status, 'advance_time')

    end subroutine advance_time

    subroutine create_axis(self, name,  points)
        class(netcdf_file), intent(inout) :: self
        character(len=*), intent(in) :: name
        real(rk), intent(in) :: points(:)
        integer(ik) :: status, dimid, varid, i
        type(netcdf_axis), allocatable :: axes(:)
        type(netcdf_axis) :: axis
        ! character(len=:), allocatable :: axis_name(:)
        ! integer(ik), allocatable :: axis_dimid(:), axis_varid(:)

        status = nf90_redef(self % ncid)
        call handle_netcdf_error(status, 'create_axis')

        do i = 1, size(self % axes)
            if (trim(self % axes(i) % name) == trim(name)) then
                write (log_str, '(a)') 'Axis by name `' // trim(name) // '` already exists.'
                call logger % fatal('create_axis', log_str)
                call abort_now()
            end if
        end do

        status = nf90_def_dim(self % ncid, trim(name), size(points), dimid)
        call handle_netcdf_error(status, 'create_axis')

        status = nf90_def_var(self % ncid, trim(name), nf90_float, dimid, varid)
        call handle_netcdf_error(status, 'create_axis')

        axis = netcdf_axis(name, dimid, varid)

        status = nf90_enddef(self % ncid)
        call handle_netcdf_error(status, 'create_axis')

        status = nf90_put_var(self % ncid, varid, points)
        call handle_netcdf_error(status, 'create_axis')

        axes = [self % axes, axis]
        self % axes = axes

    end subroutine create_axis

    subroutine create_var(self, name, dims)
        class(netcdf_file), intent(inout) :: self
        character(len=*), intent(in) :: name
        character(len=*), intent(in) :: dims(:)
        type(netcdf_var) :: var
        type(netcdf_var), allocatable :: vars(:)
        integer(ik) :: dimids(size(dims))
        integer(ik) :: status, varid, i, j, n

        do i = 1, size(self % vars)
            if (self % vars(i) % name == trim(name)) then
                write (log_str, '(a)') 'Variable by name `' // trim(name) // '` already exists.'
                call logger % fatal('create_var', log_str)
                call abort_now()
            end if
        end do

        ! Translate dimension names into dimension IDs; do time first
        ! since we have to put it as the last dimension if it's present

        dimids(:) = -1

        do i = 1, size(dims)
            if (trim(dims(i)) == 't') then
                dimids(size(dims)) = self % tdim
            end if
        end do

        n = 1
        do i = 1, size(dims)
            if (trim(dims(i)) == 't') cycle
            do j = 1, size(self % axes)
                if (trim(dims(i)) == trim(self % axes(j) % name)) then
                    dimids(n) = self % axes(j) % dimid
                end if
            end do
            if (dimids(n) == -1) then
                write (log_str, '(a)') 'Unknown dimension `' // trim(dims(i)) // '`.'
                call logger % fatal('create_var', log_str)
                call abort_now()
            end if
            n = n + 1
        end do

        status = nf90_redef(self % ncid)
        call handle_netcdf_error(status, 'create_var')

        status = nf90_def_var(self % ncid, name, nf90_float, dimids, varid)
        call handle_netcdf_error(status, 'create_var')

        status = nf90_enddef(self % ncid)
        call handle_netcdf_error(status, 'create_var')

        var = netcdf_var(name, varid, dimids)

        vars = [self % vars, var]
        self % vars = vars

    end subroutine create_var

    subroutine write_var(self, name, data)
        class(netcdf_file), intent(inout) :: self
        character(len=*), intent(in) :: name
        real(rk), intent(in) :: data(:,:)
        type(netcdf_var) :: var
        logical :: found_var
        integer(ik) :: status, i

        found_var = .false.
        do i = 1, size(self % vars)
            if (trim(name) == trim(self % vars(i) % name)) then
                var = self % vars(i)
                found_var = .true.
            end if
        end do

        if (.not. found_var) then
            write (log_str, '(a)') 'Unknown variable `' // trim(name) // '`.'
            call logger % fatal('write_var', log_str)
            call abort_now()
        end if

        if (any(var % dimids == self % tdim)) then
            status = nf90_put_var(                          &
                self % ncid, var % varid, data,             &
                start = [ 1, 1, self % tidx ],              &
                count = [ size(data, 1), size(data, 2), 1 ] )
        else
            status = nf90_put_var(                       &
                self % ncid, var % varid, data,          &
                start = [ 1, 1 ],                        &
                count = [ size(data, 1), size(data, 2) ] )
        end if

        call handle_netcdf_error(status, 'write_var')

    end subroutine write_var

    subroutine handle_netcdf_error(errcode, source)
        integer(ik), intent(in) :: errcode
        character(len=*), intent(in) :: source

        if (errcode /= nf90_noerr) then
            write (log_str, *) trim(nf90_strerror(errcode))
            call logger % fatal(source, log_str)
            call abort_now()
        end if

    end subroutine handle_netcdf_error
    
end module mod_netcdf