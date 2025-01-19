module mod_model

    use mod_kinds, only: ik, rk
    use mod_log, only: logger => main_logger, log_str
    use mod_timing, only: timing_on, timing_off
    use mod_config, only: config => main_config
    use mod_tiles, only: init_tiles
    use mod_fields, only: init_prognostic_fields
    use mod_sw_dyn, only: allocate_sw_dyn_arrays, sw_dynamics_step
    use mod_sync, only: halo_exchange, allocate_sync_buffers
    use mod_writer, only: allocate_writer, accumulate_output, write_output

    implicit none

    private
    public :: init_model, run_model

    real(rk) :: previous_write_time
    logical :: write_this_step

contains
    
    subroutine init_model()

        call init_tiles()
        call init_prognostic_fields()
        call allocate_sw_dyn_arrays()
        call allocate_sync_buffers()
        call allocate_writer()

        call timing_on('HALO')
        call halo_exchange()
        call timing_off('HALO')

        previous_write_time = 0.
        write_this_step = .false.

    end subroutine init_model

    subroutine run_model()
        integer(ik) :: n = 0
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

        ! Write the initial state
        call accumulate_output(1._rk)
        call write_output(0._rk, 0._rk)

        do while (time < config % t_final)

            n = n + 1
            dt = calculate_dt(time)
  
            if (this_image() == 1) then
                write (log_str, '(a,i6)') 'Computing time step', n
                call logger % info('run_model', log_str)
            end if

            call timing_on('SW_DYN')
            call sw_dynamics_step(dt)
            call timing_off('SW_DYN')

            call timing_on('HALO')
            call halo_exchange()
            call timing_off('HALO')

            time = time + dt

            call timing_on('ACCUMULATE_OUTPUT')
            call accumulate_output(dt)
            call timing_off('ACCUMULATE_OUTPUT')

            if (write_this_step) then
                call timing_on('WRITE_OUTPUT')
                call write_output(previous_write_time, time)
                call timing_off('WRITE_OUTPUT')

                previous_write_time = time
                write_this_step = .false.
            end if
        
        end do

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

end module mod_model