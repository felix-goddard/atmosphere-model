module mod_io

    use iso_fortran_env, only: stderr => error_unit
    use mod_kinds, only: ik, rk
    use mod_log, only: logger => main_logger, log_str
    use mod_config, only: config => main_config
    use netcdf, only: &
      nf90_create, nf90_def_dim, nf90_def_var, nf90_float, nf90_put_var,    &
      nf90_enddef, nf90_close, nf90_unlimited, nf90_clobber, nf90_strerror, &
      nf90_noerr, nf90_sync
  
    implicit none
  
    private
    public :: init_io, finalise_io, write_time_slice, advance_time

    integer(ik) :: ncid
    integer(ik) :: dimid_t, dimid_x, dimid_y
    integer(ik) :: varid_t, varid_x, varid_y, varid_h, varid_u, varid_v
    integer(ik) :: time_index
  
  contains

    subroutine init_io()
      integer(ik) :: i, status
      integer(ik) :: xsize, ysize
      real(rk) :: dx, dy, Lx, Ly

      xsize = config % nx
      ysize = config % ny
      dx = config % dx
      Lx = config % Lx
      dy = config % dy
      Ly = config % Ly
      time_index = 0 ! since we call advance_time before ever writing, start it at 0

      status = nf90_create(path='output/output.nc', cmode=nf90_clobber, ncid=ncid)
      call handle_netcdf_error(status, 'init_io')

      status = nf90_def_dim(ncid, 't', nf90_unlimited, dimid_t)
      call handle_netcdf_error(status, 'init_io')

      status = nf90_def_var(ncid, 't', nf90_float, dimid_t, varid_t)
      call handle_netcdf_error(status, 'init_io')

      status = nf90_def_dim(ncid, 'x', xsize, dimid_x)
      call handle_netcdf_error(status, 'init_io')

      status = nf90_def_var(ncid, 'x', nf90_float, dimid_x, varid_x)
      call handle_netcdf_error(status, 'init_io')

      status = nf90_def_dim(ncid, 'y', ysize, dimid_y)
      call handle_netcdf_error(status, 'init_io')

      status = nf90_def_var(ncid, 'y', nf90_float, dimid_y, varid_y)
      call handle_netcdf_error(status, 'init_io')

      status = nf90_def_var(ncid, 'h', nf90_float, [dimid_x, dimid_y, dimid_t], varid_h)
      call handle_netcdf_error(status, 'init_io')

      status = nf90_def_var(ncid, 'u', nf90_float, [dimid_x, dimid_y, dimid_t], varid_u)
      call handle_netcdf_error(status, 'init_io')

      status = nf90_def_var(ncid, 'v', nf90_float, [dimid_x, dimid_y, dimid_t], varid_v)
      call handle_netcdf_error(status, 'init_io')

      status = nf90_enddef(ncid)
      call handle_netcdf_error(status, 'init_io')

      status = nf90_put_var(ncid, varid_x, [((dx-Lx)/2. + i*dx, i = 0, xsize-1)])
      call handle_netcdf_error(status, 'init_io')

      status = nf90_put_var(ncid, varid_y, [((dy-Ly)/2. + i*dy, i = 0, ysize-1)])
      call handle_netcdf_error(status, 'init_io')

    end subroutine init_io

    subroutine finalise_io()
      integer(ik) :: status
      status = nf90_close(ncid)
      call handle_netcdf_error(status, 'finalise_io')
    end subroutine finalise_io

    subroutine handle_netcdf_error(errcode, source)
      integer(ik), intent(in) :: errcode
      character(len=*), intent(in) :: source

      if (errcode /= nf90_noerr) then
        write (log_str, *) trim(nf90_strerror(errcode))
        call logger % fatal(source, log_str)
      end if
    end subroutine handle_netcdf_error

    subroutine advance_time(time)
      real(rk), intent(in) :: time
      integer(ik) :: status

      time_index = time_index + 1

      status = nf90_put_var( ncid, varid_t, time, start=[time_index] )
      call handle_netcdf_error(status, 'advance_time')
      
    end subroutine advance_time
  
    subroutine write_time_slice(field, fieldname)
      ! Writes a field into a binary file.
      real(rk), intent(in) :: field(:,:)
      character(len=*), intent(in) :: fieldname
      integer(ik) :: status, varid

      select case (fieldname)
      case ('h')
        varid = varid_h
      case ('u')
        varid = varid_u
      case ('v')
        varid = varid_v
      case default
        write (stderr, *) ('Unknown output variable ' // fieldname)
        stop "Stopped"
      end select

      status = nf90_put_var(                          &
        ncid, varid, field,                           &
        start = [ 1, 1, time_index ],                 &
        count = [ size(field, 1), size(field, 2), 1 ] )

      call handle_netcdf_error(status, 'write_time_slice')

    end subroutine write_time_slice
  
  end module mod_io