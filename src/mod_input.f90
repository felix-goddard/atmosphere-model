module mod_input

   use mod_kinds, only: ik, rk
   use mod_log, only: logger => main_logger, log_str
   use mod_config, only: config => main_config
   use mod_netcdf, only: netcdf_file, open_netcdf
   use mod_util, only: abort_now

   implicit none

contains

   function read_initial_file() result(initial_nc)
      ! Load the initial conditions file.
      ! The initial conditions file should have dimensions for the
      ! cell centers and faces in the x and y directions, and variables
      ! for each of the dimensions (xc and yc for the cell centers, xf
      ! and yf for the cell faces) and the prognostic fields (h, u, v).

      integer(ik) :: i, idx
      real(rk) :: dx, dy
      integer(ik) :: xsize, ysize, zsize
      character(len=:), allocatable :: names(:)
      real(rk), allocatable :: points(:)
      type(netcdf_file) :: initial_nc

      initial_nc = open_netcdf(config%initial_filename)

      xsize = -1
      ysize = -1
      zsize = -1

      ! Check we have all the expected coordinates

      names = ['xc ', 'yc ', 'xf ', 'yf ', 'lev']
      do i = 1, size(names)

         idx = initial_nc%get_axis_index(names(i))

         if (idx == -1) then
            write (log_str, '(a)') 'Could not find axis `'//trim(names(i)) &
               //'` in initial condition netCDF.'
            call logger%fatal('init_prognostic_fields', log_str)
            call abort_now()
         end if

         if (names(i) == 'xc' .or. names(i) == 'xf') then

            if (xsize == -1) then
               xsize = initial_nc%axes(idx)%size
            else if (initial_nc%axes(idx)%size /= xsize) then
               write (log_str, '(a)') &
                  'Inconsistent x axis sizes in initial condition netCDF.'
               call logger%fatal('init_prognostic_fields', log_str)
               call abort_now()
            end if

         else if (names(i) == 'yc ' .or. names(i) == 'yf ') then

            if (ysize == -1) then
               ysize = initial_nc%axes(idx)%size
            else if (initial_nc%axes(idx)%size /= ysize) then
               write (log_str, '(a)') &
                  'Inconsistent y axis sizes in initial condition netCDF.'
               call logger%fatal('init_prognostic_fields', log_str)
               call abort_now()
            end if

         else if (names(i) == 'lev') then

            if (zsize == -1) then
               zsize = initial_nc%axes(idx)%size
            else if (initial_nc%axes(idx)%size /= zsize) then
               write (log_str, '(a)') &
                  'Inconsistent vertical axis sizes ' &
                  //'in initial condition netCDF.'
               call logger%fatal('init_prognostic_fields', log_str)
               call abort_now()
            end if

         end if
      end do

      ! Extract grid and time properties and put them in the config

      config%nx = xsize
      config%ny = ysize
      config%nlev = zsize

      points = initial_nc%read_axis('xc')
      dx = points(2) - points(1)
      if (.not. all(points(2:xsize) - points(1:xsize - 1) == dx)) then
         write (log_str, '(a)') &
            'Spacing of x-coordinate points in initial condition netCDF ' &
            //'is not uniform.'
         call logger%fatal('init_prognostic_fields', log_str)
         call abort_now()
      end if

      config%Lx = maxval(points) - minval(points) + dx
      config%dx = dx

      points = initial_nc%read_axis('yc')
      dy = points(2) - points(1)
      if (.not. all(points(2:ysize) - points(1:ysize - 1) == dy)) then
         write (log_str, '(a)') &
            'Spacing of y-coordinate points in initial condition netCDF ' &
            //'is not uniform.'
         call logger%fatal('init_prognostic_fields', log_str)
         call abort_now()
      end if

      config%Ly = maxval(points) - minval(points) + dy
      config%dy = dy

      if (initial_nc%t_axis%dimid == -1) then
         write (log_str, '(a)') &
            'Could not find axis `t` in initial condition netCDF.'
         call logger%fatal('init_prognostic_fields', log_str)
         call abort_now()
      end if

      points = initial_nc%read_axis('t')
      config%t_initial = points(1)
      config%t_final = config%t_final + points(1)

   end function read_initial_file

end module mod_input
