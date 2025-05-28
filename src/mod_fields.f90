module mod_fields

   use ieee_arithmetic, only: ieee_is_finite
   use mod_kinds, only: ik, rk
   use mod_log, only: logger => main_logger, log_str
   use mod_config, only: config => main_config
   use mod_netcdf, only: netcdf_file, open_netcdf
   use mod_tiles, only: is, ie, js, je, isd, ied, jsd, jed

   implicit none
   private

   public :: allocate_fields, check_stability

   ! prognostic fields
   real(rk), public, allocatable :: &
      dp(:, :, :), & ! prognostic pressure thickness
      pt(:, :, :), & ! prognostic potential temperature
      ud(:, :, :), & ! prognostic u wind (on D grid)
      vd(:, :, :), & ! prognostic v wind (on D grid)
      ts(:, :)       ! surface temperature

   ! diagnostic and auxiliary fields
   real(rk), public, allocatable :: &
      play(:, :, :), &    ! average pressure of layers
      playkap(:, :, :), & ! average pressure^kappa of layers
      plev(:, :, :), &    ! pressure on interfaces
      pkap(:, :, :), &    ! pressure^kappa on interfaces
      gz(:, :, :)         ! geopotential height on interfaces

   real(rk), public, allocatable :: &
      net_flux(:, :, :) ! net radiative flux (W/m2)

   real(rk), public, allocatable :: &
      coord_Ak(:), coord_Bk(:) ! coefficients for the hybrid σ-p coordinate

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

      if (.not. allocated(coord_Ak)) allocate (coord_Ak(nlay + 1))
      if (.not. allocated(coord_Bk)) allocate (coord_Bk(nlay + 1))

      if (.not. allocated(radius)) allocate (radius(is:ie, js:je))

   end subroutine allocate_fields

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
