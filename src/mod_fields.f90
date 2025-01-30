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

   real(rk), allocatable :: h(:, :)  ! prognostic height
   real(rk), allocatable :: ud(:, :) ! prognostic u wind (on D grid)
   real(rk), allocatable :: vd(:, :) ! prognostic v wind (on D grid)

   public :: h, ud, vd

contains

   subroutine init_prognostic_fields(initial_nc)
      ! Allocate the arrays for the prognostic fields and fill them
      ! from the initial condition file.

      type(netcdf_file), intent(in) :: initial_nc
      real(rk), allocatable :: values(:, :)
      character(len=:), allocatable :: names(:)
      integer(ik) :: i, idx

      if (.not. allocated(h)) allocate (h(is:ie, js:je))
      if (.not. allocated(ud)) allocate (ud(is:ie, js:je))
      if (.not. allocated(vd)) allocate (vd(is:ie, js:je))

      allocate (values(1:config%nx, 1:config%ny))

      names = ['h', 'u', 'v']
      do i = 1, size(names)

         idx = initial_nc%get_variable_index(names(i))

         if (idx == -1) then
            write (log_str, '(a)') 'Could not find variable `' &
               //trim(names(i))//'` in initial condition netCDF.'
            call logger%fatal('init_prognostic_fields', log_str)
            call abort_now()
         end if

         call initial_nc%read_variable( &
            names(i), values(1:config%nx, 1:config%ny))

         select case (names(i))
         case ('h')
            h(isd:ied, jsd:jed) = values(isd:ied, jsd:jed)
         case ('u')
            ud(isd:ied, jsd:jed) = values(isd:ied, jsd:jed)
         case ('v')
            vd(isd:ied, jsd:jed) = values(isd:ied, jsd:jed)
         end select
      end do

      deallocate (values)

   end subroutine init_prognostic_fields

end module mod_fields
