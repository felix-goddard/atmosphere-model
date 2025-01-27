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
        integer(ik) :: tdim = -1, tvar = -1, tidx = -1
        character(len=max_var_name_length), allocatable :: axis_name(:), var_name(:)
        integer(ik), allocatable :: axis_dimid(:), axis_varid(:), var_id(:)
    contains
        procedure :: close => close_netcdf
        procedure :: create_time_axis
        procedure :: advance_time
        procedure :: create_axis
        procedure :: create_var
        procedure :: put_timeslice
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

        out % axis_name = [character ::]
        out % var_name = [character ::]
        out % axis_dimid = [integer(ik) ::]
        out % axis_varid = [integer(ik) ::]
        out % var_id = [integer(ik) ::]

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
        character(len=:), allocatable :: axis_name(:)
        integer(ik), allocatable :: axis_dimid(:), axis_varid(:)

        status = nf90_redef(self % ncid)
        call handle_netcdf_error(status, 'create_axis')

        do i = 1, size(self % axis_varid)
            if (trim(self % axis_name(i)) == trim(name)) then
                write (log_str, '(a)') 'Axis by name `' // trim(name) // '` already exists.'
                call logger % fatal('create_axis', log_str)
                call abort_now()
            end if
        end do

        status = nf90_def_dim(self % ncid, trim(name), size(points), dimid)
        call handle_netcdf_error(status, 'create_axis')

        axis_name = [self % axis_name, trim(name)]
        axis_dimid = [self % axis_dimid, dimid]
        axis_varid = [self % axis_varid, -1]

        status = nf90_def_var(self % ncid, trim(name), nf90_float, dimid, varid)
        call handle_netcdf_error(status, 'create_axis')

        axis_varid(i+1) = varid

        status = nf90_enddef(self % ncid)
        call handle_netcdf_error(status, 'create_axis')

        status = nf90_put_var(self % ncid, varid, points)
        call handle_netcdf_error(status, 'create_axis')

        self % axis_name = axis_name
        self % axis_dimid = axis_dimid
        self % axis_varid = axis_varid

    end subroutine create_axis

    subroutine create_var(self, name, dims)
        class(netcdf_file), intent(inout) :: self
        character(len=*), intent(in) :: name
        character(len=*), intent(in) :: dims(:)
        character(len=:), allocatable :: var_name(:)
        integer(ik), allocatable :: var_id(:)
        integer(ik) :: dimids(size(dims))
        integer(ik) :: status, varid, i, j, n

        do i = 1, size(self % var_name)
            if (self % var_name(i) == trim(name)) then
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
            do j = 1, size(self % axis_name)
                if (trim(dims(i)) == trim(self % axis_name(j))) then
                    dimids(n) = self % axis_dimid(j)
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

        var_name = [self % var_name, name]
        var_id = [self % var_id, varid]

        self % var_name = var_name
        self % var_id = var_id

    end subroutine create_var

    subroutine put_timeslice(self, name, data)
        class(netcdf_file), intent(inout) :: self
        character(len=*), intent(in) :: name
        real(rk), intent(in) :: data(:,:)
        integer(ik) :: status, i, varid

        varid = -1
        do i = 1, size(self % var_name)
            if (trim(name) == trim(self % var_name(i))) &
                varid = self % var_id(i)
        end do

        if (varid == -1) then
            write (log_str, '(a)') 'Unknown variable `' // trim(name) // '`.'
            call logger % fatal('put_timeslice', log_str)
            call abort_now()
        end if

        status = nf90_put_var(                          &
            self % ncid, varid, data,                   &
            start = [ 1, 1, self % tidx ],              &
            count = [ size(data, 1), size(data, 2), 1 ] )

        call handle_netcdf_error(status, 'put_timeslice')

    end subroutine put_timeslice

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