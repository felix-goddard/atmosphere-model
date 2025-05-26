module mod_model

   use mod_kinds, only: ik, rk
   use mod_log, only: logger => main_logger, log_str
   use mod_timing, only: timing_on, timing_off
   use mod_config, only: config => main_config, read_config_file
   use mod_netcdf, only: netcdf_file
   use mod_input, only: read_initial_file
   use mod_tiles, only: init_tiles
   use mod_fields, only: init_prognostic_fields, initial_halo_exchange, &
                         check_stability
   use mod_sw_dyn, only: allocate_sw_dyn_arrays, &
                         cgrid_dynamics_step, cgrid_halo_exchange, &
                         dgrid_dynamics_step, dgrid_halo_exchange
   use mod_radiation, only: allocate_radiation_arrays, &
                            calculate_radiative_heating, radiation_halo_exchange
   use mod_physics, only: allocate_physics_arrays, physics_halo_exchange, &
                          calculate_physics_tendencies, apply_physics_tendencies
   use mod_remap, only: allocate_remapping_arrays, perform_vertical_remapping
   use mod_sync, only: halo_exchange, allocate_sync_buffers
   use mod_output, only: accumulate_output, write_output, write_restart_file

   implicit none

   private
   public :: init_model, run_model

   type(netcdf_file) :: initial_netcdf
   real(rk) :: previous_write_time, previous_radiation_time, &
               previous_physics_time, previous_remap_time
   logical :: write_this_step, calculate_radiation_this_step, &
              calculate_physics_this_step, remap_this_step, stable

contains

   subroutine init_model()

      ! Load the config file; this leaves computational domain info
      ! like nx and dx unset, which is then filled when we load the
      ! initial condition netCDF.
      call read_config_file('config.nml')
      if (this_image() == 1) &
         call logger%info('main', 'Loaded config')

      initial_netcdf = read_initial_file()
      call init_tiles()

      ! After this point, the main config is fully initialised

      call init_prognostic_fields(initial_netcdf)
      call initial_netcdf%close()

      if (this_image() == 1) &
         call logger%info('main', 'Loaded initial condition from ' &
                          //config%initial_filename)

      ! Allocate the sync buffers and perform the initial sync; this
      ! synchronises the prognostic variables and gz across all processors.
      ! We must do this before calling `allocate_sw_dyn_arrays` as that
      ! uses the values in gz
      call allocate_sync_buffers()

      call timing_on('HALO EXCHANGE')
      call initial_halo_exchange()
      call timing_off('HALO EXCHANGE')

      call allocate_sw_dyn_arrays()
      call allocate_radiation_arrays()
      call allocate_physics_arrays()
      call allocate_remapping_arrays()

      previous_radiation_time = config%t_initial
      calculate_radiation_this_step = .true.

      previous_physics_time = config%t_initial
      calculate_physics_this_step = .true.

      previous_remap_time = config%t_initial
      remap_this_step = .false.

      previous_write_time = config%t_initial
      write_this_step = .false.

      stable = .true.

      if (this_image() == 1) &
         call logger%info('main', 'Initialised model fields')

   end subroutine init_model

   subroutine run_model()
      integer(ik) :: n
      real(rk) :: time, dt

      if (this_image() == 1) then
         write (log_str, '(2(a,f6.1),a,f6.3,a)') &
            'Δx = ', (config%dx/1e3), &
            ' km; Δy = ', (config%dy/1e3), &
            ' km; Δt max. = ', (config%dt_max/60.), ' min.'
         call logger%info('run_model', log_str)

         write (log_str, '(a,f10.2,a)') &
            'Max. wavespeed is ', &
            min(config%dx, config%dy)/config%dt_max, ' m/s'
         call logger%info('run_model', log_str)

         write (log_str, '(a,f10.2,a)') &
            'Running to t = ', (config%t_final/3600.), ' h'
         call logger%info('run_model', log_str)
      end if

      n = 0
      time = config%t_initial

      main_loop: do while (time < config%t_final)

         if (calculate_radiation_this_step) then
            call timing_on('RADIATION TOTAL')
            call calculate_radiative_heating(time)
            call timing_off('RADIATION TOTAL')

            call timing_on('HALO EXCHANGE')
            call radiation_halo_exchange()
            call timing_off('HALO EXCHANGE')

            calculate_radiation_this_step = .false.
         end if

         n = n + 1
         dt = calculate_dt(time)

         if (this_image() == 1) then
            write (log_str, '(a,i6,a,f6.2)') &
               'Computing time step ', n, ' with dt=', dt
            call logger%info('run_model', log_str)
         end if

         if (calculate_physics_this_step) then
            call timing_on('PHYSICS')
            call calculate_physics_tendencies()
            call timing_off('PHYSICS')

            calculate_physics_this_step = .false.
         end if

         call timing_on('CGRID DYNAMICS')
         call cgrid_dynamics_step(dt/2.)
         call timing_off('CGRID DYNAMICS')

         call timing_on('HALO EXCHANGE')
         call cgrid_halo_exchange()
         call timing_off('HALO EXCHANGE')

         call timing_on('DGRID DYNAMICS')
         call dgrid_dynamics_step(dt)
         call timing_off('DGRID DYNAMICS')

         call timing_on('PHYSICS')
         call apply_physics_tendencies(dt)
         call timing_off('PHYSICS')

         if (remap_this_step) then
            call timing_on('REMAP')
            call perform_vertical_remapping()
            call timing_off('REMAP')

            remap_this_step = .false.
         end if

         call timing_on('HALO EXCHANGE')
         call dgrid_halo_exchange()
         call physics_halo_exchange()
         call timing_off('HALO EXCHANGE')

         stable = check_stability()
         call co_reduce(stable, and_func)
         if (.not. stable) exit main_loop

         time = time + dt

         call timing_on('ACCUMULATE OUTPUT')
         call accumulate_output(dt)
         call timing_off('ACCUMULATE OUTPUT')

         if (write_this_step) then
            call timing_on('WRITE OUTPUT')
            call write_output(.5*(previous_write_time + time))
            call timing_off('WRITE OUTPUT')

            write_this_step = .false.
         end if

      end do main_loop

      if (config%save_restart_file .and. stable) then
         call write_restart_file(config%restart_filename, time)
      end if

   end subroutine run_model

   function calculate_dt(time) result(dt)
      real(rk), intent(in) :: time
      real(rk) :: dt, candidate_dt(6)
      integer(ik) :: i
      real(rk), parameter :: time_eps = 1e-2

      candidate_dt = [config%dt_max, &
                      config%t_final - time, &
                      previous_write_time + config%dt_output - time, &
                      previous_radiation_time + config%dt_radiation - time, &
                      previous_physics_time + config%dt_physics - time, &
                      previous_remap_time + config%dt_remap - time]

      i = minloc(candidate_dt, 1, mask=candidate_dt > 0)
      dt = candidate_dt(i)

      if (abs(dt - candidate_dt(3)) < time_eps) then
         write_this_step = .true.
         previous_write_time = time + dt
      end if

      if (abs(dt - candidate_dt(4)) < time_eps) then
         calculate_radiation_this_step = .true.
         previous_radiation_time = time + dt
      end if

      if (abs(dt - candidate_dt(5)) < time_eps) then
         calculate_physics_this_step = .true.
         previous_physics_time = time + dt
      end if

      if (abs(dt - candidate_dt(6)) < time_eps) then
         remap_this_step = .true.
         previous_remap_time = time + dt
      end if

   end function calculate_dt

   pure function and_func(a, b) result(res)
      logical, intent(in) :: a, b
      logical :: res

      res = a .and. b

   end function and_func

end module mod_model
