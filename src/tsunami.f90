program main

    ! Tsunami simulator.
    !
    ! Solves the non-linear 2-d shallow water equation system:
    !
    !     du/dt + u du/dx + v du/dy + g dh/dx = 0
    !     dv/dt + u dv/dx + v dv/dy + g dh/dy = 0
    !     dh/dt + d(hu)/dx + d(hv)/dy = 0
    !
    ! This version is parallelized, and uses derived types.
  
    use mod_kinds, only: ik, rk
    use mod_io, only: init_io, finalise_io
    use mod_log, only: init_logging, logger => main_logger, log_str
    use mod_config, only: init_config, config => main_config
    use mod_model, only: init_model, run_model
  
    implicit none
  
    if (this_image() == 1) then
      call init_logging('log/log.out', 'log/log.err')
      call logger % info('main', 'Initialised logger')

      call init_config('config.nml')
      call logger % info('main', 'Loaded config file')

      call init_io()
      call logger % info('main', 'Initialised file I/O')
    end if

    call init_model()
    sync all
    if (this_image() == 1) &
      call logger % info('main', 'Initialised model fields')

    call run_model()

    if (this_image() == 1) &
      call finalise_io()
  
  end program main