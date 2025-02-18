module mod_fields

   use mod_kinds, only: ik, rk
   use mod_log, only: logger => main_logger, log_str
   use mod_config, only: config => main_config
   use mod_constants, only: top_pressure, kappa
   use mod_netcdf, only: netcdf_file, open_netcdf
   use mod_tiles, only: is, ie, js, je, isd, ied, jsd, jed
   use mod_sync, only: halo_exchange
   use mod_util, only: abort_now

   implicit none
   private

   public :: init_prognostic_fields, initial_halo_exchange

   ! prognostic fields
   real(rk), public, allocatable :: dp(:, :, :) ! prognostic pressure thickness
   real(rk), public, allocatable :: pt(:, :, :) ! prognostic potential temperature
   real(rk), public, allocatable :: ud(:, :, :) ! prognostic u wind (on D grid)
   real(rk), public, allocatable :: vd(:, :, :) ! prognostic v wind (on D grid)

   ! diagnostic and auxiliary fields
   real(rk), public, allocatable :: plev(:, :, :) ! pressure on interfaces
   real(rk), public, allocatable :: pkap(:, :, :) ! pressure^kappa on interfaces
   real(rk), public, allocatable :: gz(:, :, :) ! geopotential height on interfaces
   real(rk), public, allocatable :: heat_source(:, :, :)

contains

   subroutine allocate_fields()
      integer(ik) :: nlev

      nlev = config%nlev

      ! prognostic fields
      if (.not. allocated(dp)) allocate (dp(is:ie, js:je, nlev))
      if (.not. allocated(pt)) allocate (pt(is:ie, js:je, nlev))
      if (.not. allocated(ud)) allocate (ud(is:ie, js:je, nlev))
      if (.not. allocated(vd)) allocate (vd(is:ie, js:je, nlev))

      ! diagnostic and auxiliary fields
      if (.not. allocated(plev)) allocate (plev(is:ie, js:je, nlev + 1))
      if (.not. allocated(pkap)) allocate (pkap(is:ie, js:je, nlev + 1))
      if (.not. allocated(gz)) allocate (gz(is:ie, js:je, nlev + 1))
      if (.not. allocated(heat_source)) &
         allocate (heat_source(is:ie, js:je, nlev))

   end subroutine allocate_fields

   subroutine init_prognostic_fields(initial_nc)
      type(netcdf_file), intent(in) :: initial_nc
      real(rk), allocatable :: values(:, :, :)
      character(len=:), allocatable :: names(:)
      integer(ik) :: i, idx, j, x, y

      call allocate_fields()

      heat_source(:, :, :) = 0.

      allocate (values(1:config%nx, 1:config%ny, 1:config%nlev))

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
            names(i), values(1:config%nx, 1:config%ny, 1:config%nlev))

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

      plev(is:ie, js:je, config%nlev + 1) = top_pressure
      pkap(is:ie, js:je, config%nlev + 1) = top_pressure**kappa

      idx = initial_nc%get_variable_index('gzs')

      if (idx == -1) then
         write (log_str, '(a)') 'Could not find variable `gzs` ' &
            //'in initial condition netCDF.'
         call logger%fatal('init_prognostic_fields', log_str)
         call abort_now()
      end if

      call initial_nc%read_variable( &
         'gzs', values(1:config%nx, 1:config%ny, 1))

      gz(isd:ied, jsd:jed, 1) = values(isd:ied, jsd:jed, 1)

      deallocate (values)

   end subroutine init_prognostic_fields

   subroutine initial_halo_exchange()

      call halo_exchange(dp, pt, gz, ud, vd)
      
   end subroutine initial_halo_exchange

end module mod_fields
