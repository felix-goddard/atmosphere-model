module mod_input

   use mod_kinds, only: ik, rk
   use mod_log, only: logger => main_logger, log_str
   use mod_config, only: config => main_config
   use mod_netcdf, only: netcdf_file, open_netcdf
   use mod_tiles, only: is, ie, js, je, isd, ied, jsd, jed
   use mod_sync, only: halo_exchange
   use mod_constants, only: kappa
   use mod_fields, only: allocate_fields, dp, pt, ud, vd, gz, ts, plev, pkap, &
                         radius, coord_Ak, coord_Bk
   use mod_util, only: abort_now

   implicit none
   private

   public :: read_initial_file, init_prognostic_fields, initial_halo_exchange

contains

   function read_initial_file() result(initial_nc)
      ! Load the initial conditions file.
      ! The initial conditions file should have dimensions for the
      ! cell centers and faces in the x and y directions, and variables
      ! for each of the dimensions (xc and yc for the cell centers, xf
      ! and yf for the cell faces) and the prognostic fields (h, pt, u, v, ts).

      integer(ik) :: i, idx, xsize, ysize, zsize
      real(rk) :: dx, dy
      character(len=:), allocatable :: names(:)
      real(rk), allocatable :: points(:)
      type(netcdf_file) :: initial_nc
      real(rk), parameter :: tolerance = 1e-3

      initial_nc = open_netcdf(config%initial_filename)

      xsize = -1
      ysize = -1
      zsize = -1

      ! Check we have all the expected coordinates

      names = ['xc ', 'yc ', 'xf ', 'yf ', 'lay', 'lev']
      do i = 1, size(names)

         idx = initial_nc%get_axis_index(names(i))

         if (idx == -1) then
            write (log_str, '(a)') 'Could not find axis `'//trim(names(i)) &
               //'` in initial condition netCDF.'
            call logger%fatal('read_initial_file', log_str)
            call abort_now()
         end if

         if (names(i) == 'xc' .or. names(i) == 'xf') then

            if (xsize == -1) then
               xsize = initial_nc%axes(idx)%size
            else if (initial_nc%axes(idx)%size /= xsize) then
               write (log_str, '(a)') &
                  'Inconsistent x axis sizes in initial condition netCDF.'
               call logger%fatal('read_initial_file', log_str)
               call abort_now()
            end if

         else if (names(i) == 'yc ' .or. names(i) == 'yf ') then

            if (ysize == -1) then
               ysize = initial_nc%axes(idx)%size
            else if (initial_nc%axes(idx)%size /= ysize) then
               write (log_str, '(a)') &
                  'Inconsistent y axis sizes in initial condition netCDF.'
               call logger%fatal('read_initial_file', log_str)
               call abort_now()
            end if

         else if (names(i) == 'lay') then

            if (zsize == -1) then
               zsize = initial_nc%axes(idx)%size
            else if (initial_nc%axes(idx)%size /= zsize) then
               write (log_str, '(a)') &
                  'Inconsistent vertical axis sizes ' &
                  //'in initial condition netCDF.'
               call logger%fatal('read_initial_file', log_str)
               call abort_now()
            end if

         else if (names(i) == 'lev') then

            if (zsize == -1) then
               zsize = initial_nc%axes(idx)%size - 1
            else if (initial_nc%axes(idx)%size /= zsize + 1) then
               write (log_str, '(a)') &
                  'Inconsistent vertical axis sizes ' &
                  //'in initial condition netCDF.'
               call logger%fatal('read_initial_file', log_str)
               call abort_now()
            end if

         end if
      end do

      ! Extract grid and time properties and put them in the config

      config%nx = xsize
      config%ny = ysize
      config%nlay = zsize

      call initial_nc%read_axis('xc', points)
      dx = points(2) - points(1)
      if (.not. all(abs(points(2:xsize) - points(1:xsize - 1)) - abs(dx) &
                    < tolerance)) then
         write (log_str, '(a)') &
            'Spacing of x-coordinate points in initial condition netCDF ' &
            //'is not uniform.'
         call logger%fatal('read_initial_file', log_str)
         call abort_now()
      end if

      config%Lx = maxval(points) - minval(points) + dx
      config%dx = dx

      call initial_nc%read_axis('yc', points)
      dy = points(2) - points(1)
      if (.not. all(abs(points(2:ysize) - points(1:ysize - 1)) - abs(dy) &
                    < tolerance)) then
         write (log_str, '(a)') &
            'Spacing of y-coordinate points in initial condition netCDF ' &
            //'is not uniform.'
         call logger%fatal('read_initial_file', log_str)
         call abort_now()
      end if

      config%Ly = maxval(points) - minval(points) + dy
      config%dy = dy

      if (initial_nc%t_axis%dimid == -1) then
         write (log_str, '(a)') &
            'Could not find axis `t` in initial condition netCDF.'
         call logger%fatal('read_initial_file', log_str)
         call abort_now()
      end if

      call initial_nc%read_axis('t', points)
      config%t_initial = points(1)
      config%t_final = config%t_final + config%t_initial

   end function read_initial_file

   subroutine init_prognostic_fields(initial_nc)
      type(netcdf_file), intent(in) :: initial_nc
      real(rk), allocatable :: values(:, :, :), points(:)
      character(len=:), allocatable :: names(:)
      real(rk) :: top_pressure
      integer(ik) :: i, idx, j, x, y

      call allocate_fields()

      ! Load the coordinate coefficients

      allocate(points(1:config%nlay+1))

      call initial_nc%read_variable('ak', points)
      coord_Ak(:) = points(:)

      call initial_nc%read_variable('bk', points)
      coord_Bk(:) = points(:)

      deallocate(points)

      ! Load 3D fields

      allocate (values(1:config%nx, 1:config%ny, 1:config%nlay))

      names = ['dp', 'pt', 'u ', 'v ']
      do i = 1, size(names)

         idx = initial_nc%get_variable_index(names(i))

         if (idx == -1) then
            write (log_str, '(a)') 'Could not find variable `' &
               //trim(names(i))//'` in initial condition netCDF.'
            call logger%fatal('init_prognostic_fields', log_str)
            call abort_now()
         end if

         call initial_nc%read_variable( &
            names(i), values(1:config%nx, 1:config%ny, 1:config%nlay))

         select case (trim(names(i)))
         case ('dp')
            dp(isd:ied, jsd:jed, :) = values(isd:ied, jsd:jed, :)
         case ('pt')
            pt(isd:ied, jsd:jed, :) = values(isd:ied, jsd:jed, :)
         case ('u')
            ud(isd:ied, jsd:jed, :) = values(isd:ied, jsd:jed, :)
         case ('v')
            vd(isd:ied, jsd:jed, :) = values(isd:ied, jsd:jed, :)
         end select
      end do

      ! Load 2D fields
      names = ['gzs', 'ts ']
      do i = 1, size(names)

         idx = initial_nc%get_variable_index(names(i))

         if (idx == -1) then
            write (log_str, '(a)') 'Could not find variable `' &
               //trim(names(i))//'` in initial condition netCDF.'
            call logger%fatal('init_prognostic_fields', log_str)
            call abort_now()
         end if

         call initial_nc%read_variable( &
            names(i), values(1:config%nx, 1:config%ny, 1))

         select case (trim(names(i)))
         case ('gzs')
            gz(isd:ied, jsd:jed, 1) = values(isd:ied, jsd:jed, 1)
         case ('ts')
            ts(isd:ied, jsd:jed) = values(isd:ied, jsd:jed, 1)
         end select
      end do

      deallocate (values)

      top_pressure = coord_Ak(1)
      plev(is:ie, js:je, config%nlay + 1) = top_pressure
      pkap(is:ie, js:je, config%nlay + 1) = top_pressure**kappa

      call initial_nc%read_axis('xc', points)
      do j = jsd, jed
         radius(isd:ied, j) = points(isd:ied)**2
      end do

      call initial_nc%read_axis('yc', points)
      do i = isd, ied
         radius(i, jsd:jed) = radius(i, jsd:jed) + points(jsd:jed)**2
      end do

      radius(:, :) = sqrt(radius(:, :))

   end subroutine init_prognostic_fields

   subroutine initial_halo_exchange()

      call halo_exchange(dp, pt, gz, ud, vd)
      call halo_exchange(ts)

   end subroutine initial_halo_exchange

end module mod_input
