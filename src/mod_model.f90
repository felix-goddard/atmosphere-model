module mod_model

    use mod_kinds, only: ik, rk
    use mod_log, only: logger => main_logger, log_str
    use mod_timing, only: timing_on, timing_off
    use mod_config, only: config => main_config
    use mod_tiles, only: init_tiles
    use mod_fields, only: init_prognostic_fields
    use mod_sw_dyn, only: allocate_sw_dyn_arrays, sw_dynamics_step, is_stable
    use mod_sync, only: halo_exchange, allocate_sync_buffers
    use mod_writer, only: allocate_writer, accumulate_output, clear_accumulator, write_output

    implicit none

    private
    public :: init_model, run_model

    real(rk) :: previous_write_time
    logical :: write_this_step, stable

contains
    
    subroutine init_model()

        call init_tiles()
        call init_prognostic_fields()
        call allocate_sw_dyn_arrays()
        call allocate_sync_buffers()
        call allocate_writer()

        call timing_on('HALO EXCHANGE')
        call halo_exchange()
        call timing_off('HALO EXCHANGE')

        previous_write_time = 0.
        write_this_step = .false.
        stable = .true.

    end subroutine init_model

    subroutine run_model()
        integer(ik) :: n
        real(rk) :: time, dt

        if (this_image() == 1) then
            write (log_str, '(2(a,f6.1),a,f6.3,a)') &
                'Δx = ', (config % dx / 1e3),       &
                ' km; Δy = ', (config % dy / 1e3),  &
                ' km; Δt max. = ', (config % dt_max / 60.), ' min.'
            call logger % info('run_model', log_str)

            write (log_str, '(a,f10.2,a)') &
                'Max. wavespeed is ', min(config % dx, config % dy) / config % dt_max, ' m/s'
            call logger % info('run_model', log_str)

            write (log_str, '(a,f10.2,a)') 'Running to t = ', (config % t_final / 3600.), ' h'
            call logger % info('run_model', log_str)
        end if

        n = 0
        time = 0.

        ! Write the initial state
        call accumulate_output(1._rk)
        call write_output(0._rk, 0._rk)

        main_loop: do while (time < config % t_final)

            n = n + 1
            dt = calculate_dt(time)
  
            if (this_image() == 1) then
                write (log_str, '(a,i6,a,f6.2)') 'Computing time step ', n, ' with dt=', dt
                call logger % info('run_model', log_str)
            end if

            call timing_on('SW DYNAMICS')
            call sw_dynamics_step(dt)
            call timing_off('SW DYNAMICS')

            stable = is_stable()

            call timing_on('HALO EXCHANGE')
            call halo_exchange()
            call timing_off('HALO EXCHANGE')

            call co_reduce(stable, and_func)
            if (.not. stable) then
                if (this_image() == 1) then
                    write (log_str, '(a)') 'Instability detecting, aborting.'
                    call logger % fatal('run_model', log_str)
                end if
                exit main_loop
            end if

            time = time + dt

            call timing_on('ACCUMULATE OUTPUT')
            call accumulate_output(dt)
            call timing_off('ACCUMULATE OUTPUT')

            if (write_this_step) then
                call timing_on('WRITE OUTPUT')
                call write_output(previous_write_time, time)
                call timing_off('WRITE OUTPUT')

                previous_write_time = time
                write_this_step = .false.
            end if
        
        end do main_loop

    end subroutine run_model

    function calculate_dt(time) result(dt)
        real(rk), intent(in) :: time
        real(rk) :: dt

        dt = min(config % dt_max, config % t_final - time)

        if (previous_write_time + config % dt_output <= config % t_final &
            .and. time + dt >= previous_write_time + config % dt_output) then

            dt = previous_write_time + config % dt_output - time
            write_this_step = .true.
        end if
        
    end function calculate_dt

    pure function and_func(a, b) result(res)
        logical, intent(in) :: a, b
        logical :: res
        res = a .and. b
    end function and_func

end module mod_model