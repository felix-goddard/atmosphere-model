program main
  
    use mod_io, only: init_io, finalise_io
    use mod_log, only: init_logging, logger => main_logger, log_str, log_error
    use mod_config, only: init_config
    use mod_model, only: init_model, run_model
  
    implicit none
  
    if (this_image() == 1) then
      ! Initialise logging first.
      call init_logging('log/log.out', 'log/log.err', output_level=log_error)
      call logger % info('main', 'Initialised logger')
    end if

    ! Load the configuration file; this is needed for model setup.
    call init_config('config.nml')
    if (this_image() == 1) &
      call logger % info('main', 'Loaded config file')

    ! Initialise the model itself; this allocates arrays, sets up initial
    ! conditions, etc.
    call init_model()
    sync all
    if (this_image() == 1) then
      call logger % info('main', 'Initialised model fields')

      ! Initialise the netCDF I/O; this creates the output file and its
      ! structure (dimensions, coordinate values, etc.); this needs to
      ! know model parameters hence why we do it after model initialisation.
      call init_io()
      call logger % info('main', 'Initialised file I/O')
    end if

    ! Run the model according to the config.
    call run_model()

    ! Close the output netCDF file cleanly.
    if (this_image() == 1) &
      call finalise_io()
  
  end program main