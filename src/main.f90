program main

   use mod_log, only: init_logging, logger => main_logger, &
                      log_error, log_warning
   use mod_output, only: init_output, finalise_output
   use mod_timing, only: init_timing, print_timing
   use mod_model, only: init_model, run_model

   implicit none

   if (this_image() == 1) then
      ! Initialise logging first.
      call init_logging('log/log.out', 'log/log.err', &
                        output_level=log_error, error_level=log_warning)
      call logger%info('main', 'Initialised logger')
   end if

   ! Initialise timing
   call init_timing()
   if (this_image() == 1) &
      call logger%info('main', 'Initialised timer')

   ! Initialise the model itself; this loads the config file, allocates
   ! arrays, sets up initial conditions, etc.
   call init_model()
   sync all

   ! Initialise model output; this allocates the output accumulation
   ! arrays and on the first image creates the output netCDF file.
   ! This needs to know model parameters hence why we do it after
   ! model initialisation.
   call init_output()
   if (this_image() == 1) &
      call logger%info('main', 'Initialised output')

   ! Run the model according to the config.
   call run_model()

   ! Close the output netCDF file cleanly and print timing information.
   if (this_image() == 1) then
      call finalise_output()
      call print_timing()
   end if

end program main
