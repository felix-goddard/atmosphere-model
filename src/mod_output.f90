module mod_output

   use mod_kinds, only: ik, rk
   use mod_config, only: config => main_config
   use mod_netcdf, only: netcdf_file, create_netcdf
   use mod_tiles, only: isd, ied, jsd, jed
   use mod_fields, only: h, ud, vd

   implicit none

   private
   public :: init_output, finalise_output, &
             accumulate_output, write_output, write_restart_file

   type(netcdf_file) :: output_nc

   integer(ik), parameter :: n_output_fields = 3 ! number of outputs

   real(rk), allocatable :: output_field(:, :, :)
   real(rk), allocatable :: compensation(:, :, :)
   real(rk)              :: accumulation_time
   real(rk), allocatable :: gather_coarray(:, :) [:]
   real(rk), allocatable :: gather(:, :)

   ! these indices into `output_field` determine which variable we are
   ! accumulating or outputting; make sure these are unique
   integer(ik), parameter :: H_IDX = 1
   integer(ik), parameter :: U_IDX = 2
   integer(ik), parameter :: V_IDX = 3

contains

   subroutine init_output()
      integer(ik) :: i

      ! Allocate output arrays

      if (.not. allocated(output_field)) &
         allocate (output_field(isd:ied, jsd:jed, n_output_fields))

      if (.not. allocated(compensation)) &
         allocate (compensation(isd:ied, jsd:jed, n_output_fields))

      accumulation_time = 0.

      if (.not. allocated(gather_coarray)) &
         allocate (gather_coarray(1:config%nx, 1:config%ny) [*])

      if (this_image() == 1) then

         if (.not. allocated(gather)) &
            allocate (gather(1:config%nx, 1:config%ny))

         ! Create output netCDF

         output_nc = create_netcdf('output/output.nc')

         call output_nc%create_time_axis()

         call output_nc%create_axis( &
            'x', [(-.5*config%Lx + (i + .5)*config%dx, i=0, config%nx - 1)])

         call output_nc%create_axis( &
            'y', [(-.5*config%Ly + (i + .5)*config%dy, i=0, config%ny - 1)])

         call output_nc%create_variable('h', ['t', 'x', 'y'])
         call output_nc%create_variable('u', ['t', 'x', 'y'])
         call output_nc%create_variable('v', ['t', 'x', 'y'])

      end if

   end subroutine init_output

   subroutine finalise_output()

      call output_nc%close()

   end subroutine finalise_output

   subroutine accumulate_output(dt)
      real(rk), intent(in) :: dt
      integer(ik) :: i, j

      do i = isd, ied
         do j = jsd, jed
            ! height
            call compensated_sum(dt*h(i, j), i, j, H_IDX)

            ! u wind
            call compensated_sum(dt*.5*(ud(i, j) + ud(i, j + 1)), i, j, U_IDX)

            ! v wind
            call compensated_sum(dt*.5*(vd(i, j) + vd(i + 1, j)), i, j, V_IDX)
         end do
      end do

      accumulation_time = accumulation_time + dt

   end subroutine accumulate_output

   subroutine compensated_sum(val, i, j, idx)
      integer(ik), intent(in) :: i, j, idx
      real(rk), intent(in) :: val
      real(rk) :: y, sum, c

      y = val - compensation(i, j, idx)
      sum = output_field(i, j, idx) + y
      c = (sum - output_field(i, j, idx)) - y
      output_field(i, j, idx) = sum

   end subroutine compensated_sum

   subroutine write_output(time)
      real(rk), intent(in) :: time

      if (this_image() == 1) call output_nc%advance_time(time)

      output_field(:, :, :) = output_field(:, :, :)/accumulation_time

      ! Output height field
      call write (output_nc, output_field(:, :, H_IDX), 'h')

      ! Output u wind
      call write (output_nc, output_field(:, :, U_IDX), 'u')

      ! Output v wind
      call write (output_nc, output_field(:, :, V_IDX), 'v')

      output_field(:, :, :) = 0.
      compensation(:, :, :) = 0.
      accumulation_time = 0.

   end subroutine write_output

   subroutine write_restart_file(filename, time)
      character(len=*), intent(in) :: filename
      real(rk), intent(in) :: time
      type(netcdf_file) :: restart_nc
      integer(ik) :: i

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

      call restart_nc%create_variable('h', ['xc', 'yc'])
      call restart_nc%create_variable('u', ['xc', 'yf'])
      call restart_nc%create_variable('v', ['xf', 'yc'])

      call write (restart_nc, h(isd:ied, jsd:jed), 'h')
      call write (restart_nc, ud(isd:ied, jsd:jed), 'u')
      call write (restart_nc, vd(isd:ied, jsd:jed), 'v')

      call restart_nc%close()

   end subroutine write_restart_file

   subroutine write (netcdf, field, name)
      type(netcdf_file), intent(inout) :: netcdf
      real(rk), intent(in) :: field(isd:ied, jsd:jed)
      character(len=*), intent(in) :: name

      sync all
      gather_coarray(isd:ied, jsd:jed) [1] = field(:, :)
      sync all

      if (this_image() == 1) then
         gather(:, :) = gather_coarray(:, :)
         call netcdf%write_variable(name, gather)
      end if
   end subroutine write

end module mod_output
