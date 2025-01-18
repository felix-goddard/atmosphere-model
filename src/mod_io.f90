module mod_io

    use iso_fortran_env, only: stderr => error_unit
    use mod_kinds, only: ik, rk
    use mod_config, only: config => main_config
    use netcdf, only: &
      nf90_create, nf90_def_dim, nf90_def_var, nf90_float, nf90_put_var,    &
      nf90_enddef, nf90_close, nf90_unlimited, nf90_clobber, nf90_strerror, &
      nf90_noerr, nf90_sync
  
    implicit none
  
    private
    public :: init_io, finalise_io, write_field

    integer(ik) :: ncid
    integer(ik) :: dimid_t, dimid_x, dimid_y
    integer(ik) :: varid_t, varid_x, varid_y, varid_h, varid_u, varid_v
  
  contains

    subroutine init_io()
      integer(ik) :: i, status
      integer(ik) :: xsize, ysize
      real(rk) :: dx, dy

      xsize = config % nx
      ysize = config % ny
      dx = config % Lx / (config % nx - 1)
      dy = config % Ly / (config % ny - 1)

      status = nf90_create(path='output/tsunami.nc', cmode=nf90_clobber, ncid=ncid)
      call handle_netcdf_error(status)

      status = nf90_def_dim(ncid, 't', nf90_unlimited, dimid_t)
      call handle_netcdf_error(status)

      status = nf90_def_var(ncid, 't', nf90_float, dimid_t, varid_t)
      call handle_netcdf_error(status)

      status = nf90_def_dim(ncid, 'x', xsize, dimid_x)
      call handle_netcdf_error(status)

      status = nf90_def_var(ncid, 'x', nf90_float, dimid_x, varid_x)
      call handle_netcdf_error(status)

      status = nf90_def_dim(ncid, 'y', ysize, dimid_y)
      call handle_netcdf_error(status)

      status = nf90_def_var(ncid, 'y', nf90_float, dimid_y, varid_y)
      call handle_netcdf_error(status)

      status = nf90_def_var(ncid, 'h', nf90_float, [dimid_x, dimid_y, dimid_t], varid_h)
      call handle_netcdf_error(status)

      status = nf90_def_var(ncid, 'u', nf90_float, [dimid_x, dimid_y, dimid_t], varid_u)
      call handle_netcdf_error(status)

      status = nf90_def_var(ncid, 'v', nf90_float, [dimid_x, dimid_y, dimid_t], varid_v)
      call handle_netcdf_error(status)

      status = nf90_enddef(ncid)
      call handle_netcdf_error(status)

      status = nf90_put_var(ncid, varid_x, [(i*dx, i = 0, xsize-1)])
      call handle_netcdf_error(status)

      status = nf90_put_var(ncid, varid_y, [(i*dy, i = 0, ysize-1)])
      call handle_netcdf_error(status)

    end subroutine init_io

    subroutine finalise_io()
      integer(ik) :: status
      status = nf90_close(ncid)
      call handle_netcdf_error(status)
    end subroutine finalise_io

    subroutine handle_netcdf_error(errcode)
      integer(ik), intent(in) :: errcode

      if (errcode /= nf90_noerr) then
        write (stderr, *) trim(nf90_strerror(errcode))
        stop "Stopped"
      end if
    end subroutine handle_netcdf_error
  
    subroutine write_field(field, nx, ny, fieldname, time)
      ! Writes a field into a binary file.
      real(rk),     intent(in) :: field(:,:)
      integer(ik),   intent(in) :: nx, ny
      character(len=*), intent(in) :: fieldname
      integer(ik),   intent(in) :: time
      integer(ik) :: status, varid

      status = nf90_put_var( ncid, varid_t, time * config % dt, start=[time] )
      call handle_netcdf_error(status)

      if (fieldname == 'h') then
        varid = varid_h
      else if (fieldname == 'u') then
        varid = varid_u
      else if (fieldname == 'v') then
        varid = varid_v
      else
        write (stderr, *) ('Unknown output variable ' // fieldname)
        stop "Stopped"
      end if

      status = nf90_put_var(       &
        ncid, varid, field,        &
        start = [ 1, 1, time ],    &
        count = [ nx, ny, 1 ]      )

      call handle_netcdf_error(status)

    end subroutine write_field
  
  end module mod_io