module mod_model

   use mod_kinds, only: ik, rk
   use mod_log, only: logger => main_logger, log_str
   use mod_timing, only: timing_on, timing_off
   use mod_config, only: config => main_config, read_config_file
   use mod_netcdf, only: netcdf_file
   use mod_input, only: read_initial_file
   use mod_tiles, only: init_tiles
   use mod_fields, only: init_prognostic_fields
   use mod_sw_dyn, only: allocate_sw_dyn_arrays, cgrid_dynamics_step, &
                         dgrid_dynamics_step, is_stable
   use mod_sync, only: halo_exchange, allocate_sync_buffers
   use mod_output, only: accumulate_output, write_output, write_restart_file

   implicit none

   private
   public :: init_model, run_model

   type(netcdf_file) :: initial_netcdf
   real(rk) :: previous_write_time
   logical :: write_this_step, stable

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

      call allocate_sw_dyn_arrays()
      call allocate_sync_buffers()

      call timing_on('HALO EXCHANGE')
      call halo_exchange()
      call timing_off('HALO EXCHANGE')

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

         n = n + 1
         dt = calculate_dt(time)

         if (this_image() == 1) then
            write (log_str, '(a,i6,a,f6.2)') &
               'Computing time step ', n, ' with dt=', dt
            call logger%info('run_model', log_str)
         end if

         call timing_on('CGRID DYNAMICS')
         call cgrid_dynamics_step(dt/2.)
         call timing_off('CGRID DYNAMICS')

         call timing_on('DGRID DYNAMICS')
         call dgrid_dynamics_step(dt)
         call timing_off('DGRID DYNAMICS')

         stable = is_stable()

         call timing_on('HALO EXCHANGE')
         call halo_exchange()
         call timing_off('HALO EXCHANGE')

         call co_reduce(stable, and_func)
         if (.not. stable) then
            if (this_image() == 1) then
               write (log_str, '(a)') 'Instability detecting, aborting.'
               call logger%fatal('run_model', log_str)
            end if
            exit main_loop
         end if

         time = time + dt

         call timing_on('ACCUMULATE OUTPUT')
         call accumulate_output(dt)
         call timing_off('ACCUMULATE OUTPUT')

         if (write_this_step) then
            call timing_on('WRITE OUTPUT')
            call write_output(.5*(previous_write_time + time))
            call timing_off('WRITE OUTPUT')

            previous_write_time = time
            write_this_step = .false.
         end if

      end do main_loop

      if (config%save_restart_file .and. stable) then
         call write_restart_file(config%restart_filename, time)
      end if

   end subroutine run_model

   function calculate_dt(time) result(dt)
      real(rk), intent(in) :: time
      real(rk) :: dt

      dt = min(config%dt_max, config%t_final - time)

      if (previous_write_time + config%dt_output <= config%t_final &
          .and. time + dt >= previous_write_time + config%dt_output) then

         dt = previous_write_time + config%dt_output - time
         write_this_step = .true.
      end if

   end function calculate_dt

   pure function and_func(a, b) result(res)
      logical, intent(in) :: a, b
      logical :: res

      res = a .and. b

   end function and_func

end module mod_model
