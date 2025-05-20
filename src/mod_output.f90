module mod_output

   use mod_kinds, only: ik, rk
   use mod_config, only: config => main_config
   use mod_netcdf, only: netcdf_file, create_netcdf
   use mod_tiles, only: isd, ied, jsd, jed
   use mod_fields, only: dp, pt, ud, vd, gz, ts, net_flux

   implicit none

   private
   public :: init_output, finalise_output, &
             accumulate_output, write_output, write_restart_file

   type(netcdf_file) :: output_nc
   type(netcdf_file) :: restart_nc

   integer(ik), parameter :: n_output_fields = 6 ! number of outputs

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
   integer(ik), parameter :: NFLX_IDX = 6

contains

   subroutine init_output()
      integer(ik) :: i, nlay

      nlay = config%nlay

      ! Allocate output arrays

      if (.not. allocated(output_field)) &
         allocate (output_field(isd:ied, jsd:jed, nlay + 1, n_output_fields))

      if (.not. allocated(compensation)) &
         allocate (compensation(isd:ied, jsd:jed, nlay + 1, n_output_fields))

      output_field(:, :, :, :) = 0.
      compensation(:, :, :, :) = 0.
      accumulation_time = 0.

      if (.not. allocated(gather_coarray)) &
         allocate (gather_coarray(config%nx, config%ny, nlay + 1) [*])

      if (this_image() == 1) then

         if (.not. allocated(gather)) &
            allocate (gather(config%nx, config%ny, nlay + 1))

         ! Create output netCDF

         output_nc = create_netcdf('output/output.nc')

         call output_nc%create_time_axis()

         call output_nc%create_axis( &
            'x', [(-.5*config%Lx + (i + .5)*config%dx, i=0, config%nx - 1)])

         call output_nc%create_axis( &
            'y', [(-.5*config%Ly + (i + .5)*config%dy, i=0, config%ny - 1)])

         call output_nc%create_axis('lay', [(i, i=1, nlay)])
         call output_nc%create_axis('lev', [(i, i=1, nlay + 1)])

         call output_nc%create_variable('dp', ['t  ', 'x  ', 'y  ', 'lay'])
         call output_nc%create_variable('pt', ['t  ', 'x  ', 'y  ', 'lay'])
         call output_nc%create_variable('u', ['t  ', 'x  ', 'y  ', 'lay'])
         call output_nc%create_variable('v', ['t  ', 'x  ', 'y  ', 'lay'])
         call output_nc%create_variable('ts', ['t  ', 'x  ', 'y  '])
         call output_nc%create_variable('nflx', ['t  ', 'x  ', 'y  ', 'lev'])

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

            do k = 1, config%nlay
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

               ! net radiative flux
               call compensated_sum(dt*net_flux(i, j, k), i, j, k, NFLX_IDX)
            end do

            ! net radiative flux
            call compensated_sum(dt*net_flux(i, j, config%nlay + 1), &
                                 i, j, config%nlay + 1, NFLX_IDX)
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
      integer(ik) :: nlay

      nlay = config%nlay

      if (this_image() == 1) call output_nc%advance_time(time)

      output_field(:, :, :, :) = output_field(:, :, :, :)/accumulation_time

      ! Output pressure thickness field
      call write_3d(output_nc, output_field(:, :, 1:nlay, DP_IDX), nlay, 'dp')

      ! Output potential temperature field
      call write_3d(output_nc, output_field(:, :, 1:nlay, PT_IDX), nlay, 'pt')

      ! Output u wind
      call write_3d(output_nc, output_field(:, :, 1:nlay, U_IDX), nlay, 'u')

      ! Output v wind
      call write_3d(output_nc, output_field(:, :, 1:nlay, V_IDX), nlay, 'v')

      ! Output surface temperature
      call write_2d(output_nc, output_field(:, :, 1, TS_IDX), 'ts')

      ! Output net radiative flux
      call write_3d(output_nc, output_field(:, :, 1:nlay + 1, NFLX_IDX), &
                    nlay + 1, 'nflx')

      output_field(:, :, :, :) = 0.
      compensation(:, :, :, :) = 0.
      accumulation_time = 0.

   end subroutine write_output

   subroutine write_restart_file(filename, time)
      character(len=*), intent(in) :: filename
      real(rk), intent(in) :: time
      integer(ik) :: i, nlay

      nlay = config%nlay

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
         call restart_nc%create_axis('lay', [(i, i=1, nlay)])

         call restart_nc%create_variable('dp', ['xc ', 'yc ', 'lay'])
         call restart_nc%create_variable('pt', ['xc ', 'yc ', 'lay'])
         call restart_nc%create_variable('u', ['xc ', 'yf ', 'lay'])
         call restart_nc%create_variable('v', ['xf ', 'yc ', 'lay'])
         call restart_nc%create_variable('gzs', ['xc ', 'yc '])
         call restart_nc%create_variable('ts', ['xc ', 'yc '])

      end if

      call write_3d(restart_nc, dp(isd:ied, jsd:jed, :), nlay, 'dp')
      call write_3d(restart_nc, pt(isd:ied, jsd:jed, :), nlay, 'pt')
      call write_3d(restart_nc, ud(isd:ied, jsd:jed, :), nlay, 'u')
      call write_3d(restart_nc, vd(isd:ied, jsd:jed, :), nlay, 'v')
      call write_2d(restart_nc, gz(isd:ied, jsd:jed, 1), 'gzs')
      call write_2d(restart_nc, ts(isd:ied, jsd:jed), 'ts')

      if (this_image() == 1) call restart_nc%close()

   end subroutine write_restart_file

   subroutine write_3d(netcdf, field, k_max, name)
      type(netcdf_file), intent(inout) :: netcdf
      real(rk), intent(in) :: field(isd:ied, jsd:jed, 1:k_max)
      integer(ik) :: k_max
      character(len=*), intent(in) :: name

      sync all
      gather_coarray(isd:ied, jsd:jed, 1:k_max) [1] = field(:, :, 1:k_max)
      sync all

      if (this_image() == 1) then
         gather(:, :, 1:k_max) = gather_coarray(:, :, 1:k_max)
         call netcdf%write_variable(name, gather(:, :, 1:k_max))
      end if
   end subroutine write_3d

   subroutine write_2d(netcdf, field, name)
      type(netcdf_file), intent(inout) :: netcdf
      real(rk), intent(in) :: field(isd:ied, jsd:jed)
      character(len=*), intent(in) :: name

      sync all
      gather_coarray(isd:ied, jsd:jed, 1) [1] = field(:, :)
      sync all

      if (this_image() == 1) then
         gather(:, :, 1) = gather_coarray(:, :, 1)
         call netcdf%write_variable(name, gather(:, :, 1))
      end if
   end subroutine write_2d

end module mod_output
