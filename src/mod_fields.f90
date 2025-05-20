module mod_fields

   use ieee_arithmetic, only: ieee_is_finite
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

   public :: init_prognostic_fields, initial_halo_exchange, check_stability

   ! prognostic fields
   real(rk), public, allocatable :: dp(:, :, :) ! prognostic pressure thickness
   real(rk), public, allocatable :: pt(:, :, :) ! prognostic potential temperature
   real(rk), public, allocatable :: ud(:, :, :) ! prognostic u wind (on D grid)
   real(rk), public, allocatable :: vd(:, :, :) ! prognostic v wind (on D grid)
   real(rk), public, allocatable :: ts(:, :)    ! surface temperature

   ! diagnostic and auxiliary fields
   real(rk), public, allocatable :: play(:, :, :) ! average pressure of layers
   real(rk), public, allocatable :: playkap(:, :, :) ! average pressure^kappa of layers
   real(rk), public, allocatable :: plev(:, :, :) ! pressure on interfaces
   real(rk), public, allocatable :: pkap(:, :, :) ! pressure^kappa on pr
   real(rk), public, allocatable :: gz(:, :, :) ! geopotential height on interfaces

   real(rk), public, allocatable :: net_flux(:, :, :) ! net radiative flux (W/m2)

   real(rk), public, allocatable :: radius(:, :)

contains

   subroutine allocate_fields()
      integer(ik) :: nlay

      nlay = config%nlay

      ! prognostic fields
      if (.not. allocated(dp)) allocate (dp(is:ie, js:je, nlay))
      if (.not. allocated(pt)) allocate (pt(is:ie, js:je, nlay))
      if (.not. allocated(ud)) allocate (ud(is:ie, js:je, nlay))
      if (.not. allocated(vd)) allocate (vd(is:ie, js:je, nlay))
      if (.not. allocated(ts)) allocate (ts(is:ie, js:je))

      ! diagnostic and auxiliary fields
      if (.not. allocated(play)) allocate (play(is:ie, js:je, nlay))
      if (.not. allocated(playkap)) allocate (playkap(is:ie, js:je, nlay))
      if (.not. allocated(plev)) allocate (plev(is:ie, js:je, nlay + 1))
      if (.not. allocated(pkap)) allocate (pkap(is:ie, js:je, nlay + 1))
      if (.not. allocated(gz)) allocate (gz(is:ie, js:je, nlay + 1))

      if (.not. allocated(net_flux)) allocate (net_flux(is:ie, js:je, nlay + 1))

      if (.not. allocated(radius)) allocate (radius(is:ie, js:je))

   end subroutine allocate_fields

   subroutine init_prognostic_fields(initial_nc)
      type(netcdf_file), intent(in) :: initial_nc
      real(rk), allocatable :: values(:, :, :), coord_points(:)
      character(len=:), allocatable :: names(:)
      integer(ik) :: i, idx, j, x, y

      call allocate_fields()

      allocate (values(1:config%nx, 1:config%ny, 1:config%nlay))

      ! Load 3D fields
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

      plev(is:ie, js:je, config%nlay + 1) = top_pressure
      pkap(is:ie, js:je, config%nlay + 1) = top_pressure**kappa

      call initial_nc%read_axis('xc', coord_points)
      do j = js, je
         radius(is:ie, j) = coord_points(is:ie)**2
      end do

      call initial_nc%read_axis('yc', coord_points)
      do i = is, ie
         radius(i, js:je) = radius(i, js:je) + coord_points(js:je)**2
      end do

      radius(:, :) = sqrt(radius(:, :))

   end subroutine init_prognostic_fields

   subroutine initial_halo_exchange()

      call halo_exchange(dp, pt, gz, ud, vd)
      call halo_exchange(ts)

   end subroutine initial_halo_exchange

   function check_stability() result(is_stable)
      logical :: is_stable

      is_stable = .true.

      if (.not. all(ieee_is_finite(dp(isd:ied, jsd:jed, :)))) then
         is_stable = .false.
         write (log_str, '(2(a,f8.4))') &
            'Instability detected in pressure thickness; minval = ', &
            minval(dp), '; maxval = ', maxval(dp)
         call logger%fatal('run_model', log_str)
      end if

      if (.not. all(ieee_is_finite(pt(isd:ied, jsd:jed, :)))) then
         is_stable = .false.
         write (log_str, '(2(a,f8.4))') &
            'Instability detected in potential temperature; minval = ', &
            minval(pt), '; maxval = ', maxval(pt)
         call logger%fatal('run_model', log_str)
      end if

      if (.not. all(ieee_is_finite(ud(isd:ied, jsd:jed, :)))) then
         is_stable = .false.
         write (log_str, '(2(a,f8.4))') &
            'Instability detected in u wind; minval = ', &
            minval(ud), '; maxval = ', maxval(ud)
         call logger%fatal('run_model', log_str)
      end if

      if (.not. all(ieee_is_finite(vd(isd:ied, jsd:jed, :)))) then
         is_stable = .false.
         write (log_str, '(2(a,f8.4))') &
            'Instability detected in v wind; minval = ', &
            minval(vd), '; maxval = ', maxval(vd)
         call logger%fatal('run_model', log_str)
      end if

   end function check_stability

end module mod_fields
