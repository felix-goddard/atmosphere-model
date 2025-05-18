module mod_output

   use mod_kinds, only: ik, rk
   use mod_config, only: config => main_config
   use mod_netcdf, only: netcdf_file, create_netcdf
   use mod_tiles, only: isd, ied, jsd, jed
   use mod_fields, only: dp, pt, ud, vd, gz, ts

   implicit none

   private
   public :: init_output, finalise_output, &
             accumulate_output, write_output, write_restart_file

   type(netcdf_file) :: output_nc
   type(netcdf_file) :: restart_nc

   integer(ik), parameter :: n_output_fields = 5 ! number of outputs

   real(rk), allocatable :: output_field(:, :, :, :)
   real(rk), allocatable :: compensation(:, :, :, :)
   real(rk)              :: accumulation_time
   real(rk), allocatable :: gather_coarray(:, :, :) [:]
   real(rk), allocatable :: gather(:, :, :)

   ! these indices into `output_field` determine which variable we are
   ! accumulating or outputting; make sure these are unique
   integer(ik), parameter :: DP_IDX = 1
   integer(ik), parameter :: PT_IDX = 2
   integer(ik), parameter :: U_IDX = 3
   integer(ik), parameter :: V_IDX = 4
   integer(ik), parameter :: TS_IDX = 5

contains

   subroutine init_output()
      integer(ik) :: i

      ! Allocate output arrays

      if (.not. allocated(output_field)) &
         allocate (output_field(isd:ied, jsd:jed, config%nlev, n_output_fields))

      if (.not. allocated(compensation)) &
         allocate (compensation(isd:ied, jsd:jed, config%nlev, n_output_fields))

      output_field(:, :, :, :) = 0.
      compensation(:, :, :, :) = 0.
      accumulation_time = 0.

      if (.not. allocated(gather_coarray)) &
         allocate (gather_coarray(config%nx, config%ny, config%nlev) [*])

      if (this_image() == 1) then

         if (.not. allocated(gather)) &
            allocate (gather(config%nx, config%ny, config%nlev))

         ! Create output netCDF

         output_nc = create_netcdf('output/output.nc')

         call output_nc%create_time_axis()

         call output_nc%create_axis( &
            'x', [(-.5*config%Lx + (i + .5)*config%dx, i=0, config%nx - 1)])

         call output_nc%create_axis( &
            'y', [(-.5*config%Ly + (i + .5)*config%dy, i=0, config%ny - 1)])

         call output_nc%create_axis('lev', [(i, i=1, config%nlev)])

         call output_nc%create_variable('dp', ['t  ', 'x  ', 'y  ', 'lev'])
         call output_nc%create_variable('pt', ['t  ', 'x  ', 'y  ', 'lev'])
         call output_nc%create_variable('u', ['t  ', 'x  ', 'y  ', 'lev'])
         call output_nc%create_variable('v', ['t  ', 'x  ', 'y  ', 'lev'])
         call output_nc%create_variable('ts', ['t  ', 'x  ', 'y  '])

      end if

   end subroutine init_output

   subroutine finalise_output()

      call output_nc%close()

   end subroutine finalise_output

   subroutine accumulate_output(dt)
      real(rk), intent(in) :: dt
      integer(ik) :: i, j, k

      do i = isd, ied
         do j = jsd, jed
            ! surface temperature
            call compensated_sum(dt*ts(i, j), i, j, 1, TS_IDX)

            do k = 1, config%nlev
               ! pressure thickness
               call compensated_sum(dt*dp(i, j, k), i, j, k, DP_IDX)

               ! potential temperature
               call compensated_sum(dt*pt(i, j, k), i, j, k, PT_IDX)

               ! u wind
               call compensated_sum(dt*.5*(ud(i, j, k) + ud(i, j + 1, k)), &
                                    i, j, k, U_IDX)

               ! v wind
               call compensated_sum(dt*.5*(vd(i, j, k) + vd(i + 1, j, k)), &
                                    i, j, k, V_IDX)
            end do
         end do
      end do

      accumulation_time = accumulation_time + dt

   end subroutine accumulate_output

   subroutine compensated_sum(val, i, j, k, idx)
      integer(ik), intent(in) :: i, j, k, idx
      real(rk), intent(in) :: val
      real(rk) :: y, sum, c

      y = val - compensation(i, j, k, idx)
      sum = output_field(i, j, k, idx) + y
      c = (sum - output_field(i, j, k, idx)) - y
      output_field(i, j, k, idx) = sum

   end subroutine compensated_sum

   subroutine write_output(time)
      real(rk), intent(in) :: time

      if (this_image() == 1) call output_nc%advance_time(time)

      output_field(:, :, :, :) = output_field(:, :, :, :)/accumulation_time

      ! Output pressure thickness field
      call write (output_nc, output_field(:, :, :, DP_IDX), 'dp')

      ! Output potential temperature field
      call write (output_nc, output_field(:, :, :, PT_IDX), 'pt')

      ! Output u wind
      call write (output_nc, output_field(:, :, :, U_IDX), 'u')

      ! Output v wind
      call write (output_nc, output_field(:, :, :, V_IDX), 'v')

      ! Output surface temperature
      call write (output_nc, output_field(:, :, 1, TS_IDX), 'ts')

      output_field(:, :, :, :) = 0.
      compensation(:, :, :, :) = 0.
      accumulation_time = 0.

   end subroutine write_output

   subroutine write_restart_file(filename, time)
      character(len=*), intent(in) :: filename
      real(rk), intent(in) :: time
      integer(ik) :: i

      if (this_image() == 1) then

         restart_nc = create_netcdf(filename)

         call restart_nc%create_time_axis()
         call restart_nc%advance_time(time)

         ! x centers
         call restart_nc%create_axis( &
            'xc', [(-.5*config%Lx + (i + .5)*config%dx, i=0, config%nx - 1)])

         ! x faces
         call restart_nc%create_axis( &
            'xf', [(-.5*config%Lx + i*config%dx, i=0, config%nx - 1)])

         ! y centers
         call restart_nc%create_axis( &
            'yc', [(-.5*config%Ly + (i + .5)*config%dy, i=0, config%ny - 1)])

         ! y faces
         call restart_nc%create_axis( &
            'yf', [(-.5*config%Ly + i*config%dy, i=0, config%ny - 1)])

         ! levels
         call restart_nc%create_axis('lev', [(i, i=1, config%nlev)])

         call restart_nc%create_variable('dp', ['xc ', 'yc ', 'lev'])
         call restart_nc%create_variable('pt', ['xc ', 'yc ', 'lev'])
         call restart_nc%create_variable('u', ['xc ', 'yf ', 'lev'])
         call restart_nc%create_variable('v', ['xf ', 'yc ', 'lev'])
         call restart_nc%create_variable('gzs', ['xc ', 'yc '])
         call restart_nc%create_variable('ts', ['xc ', 'yc '])

      end if

      call write (restart_nc, dp(isd:ied, jsd:jed, :), 'dp')
      call write (restart_nc, pt(isd:ied, jsd:jed, :), 'pt')
      call write (restart_nc, ud(isd:ied, jsd:jed, :), 'u')
      call write (restart_nc, vd(isd:ied, jsd:jed, :), 'v')
      call write (restart_nc, gz(isd:ied, jsd:jed, 1), 'gzs')
      call write (restart_nc, ts(isd:ied, jsd:jed), 'ts')

      if (this_image() == 1) call restart_nc%close()

   end subroutine write_restart_file

   subroutine write (netcdf, field, name)
      type(netcdf_file), intent(inout) :: netcdf
      real(rk), intent(in) :: field(isd:ied, jsd:jed, config%nlev)
      character(len=*), intent(in) :: name

      sync all
      gather_coarray(isd:ied, jsd:jed, :) [1] = field(:, :, :)
      sync all

      if (this_image() == 1) then
         gather(:, :, :) = gather_coarray(:, :, :)
         call netcdf%write_variable(name, gather)
      end if
   end subroutine write

end module mod_output
