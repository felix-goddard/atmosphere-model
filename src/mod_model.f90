module mod_model

    use mod_kinds, only: ik, rk
    use mod_log, only: logger => main_logger, log_str
    use mod_config, only: config => main_config
    use mod_fields, only: init_prognostic_fields
    use mod_sw_dyn, only: allocate_sw_dyn_arrays, sw_dynamics_step
    use mod_sync, only: halo_exchange, init_halo_sync
    use mod_writer, only: write_output, allocate_writer

    implicit none

    private
    public :: init_model, run_model

contains
    
    subroutine init_model()

        call init_prognostic_fields()
        call allocate_sw_dyn_arrays()
        call init_halo_sync()
        call halo_exchange()
        call allocate_writer()

    end subroutine init_model

    subroutine run_model()
        integer(ik) :: n

        if (this_image() == 1) then
            write (log_str, '(2(a,f6.1),a,f6.3,a)') &
                'Δx = ', (config % dx / 1e3),       &
                ' km; Δy = ', (config % dy / 1e3),  &
                ' km; Δt = ', (config % dt / 60.), ' min.'
            call logger % info('run_model', log_str)

            write (log_str, '(a,f10.2,a)') &
                'Max. wavespeed is ', min(config % dx, config % dy) / config % dt, ' m/s'
            call logger % info('run_model', log_str)

            write (log_str, '(a,f6.2,a,i6,a)') &
                'Running to t = ', (config % nt * config % dt / 3600.), &
                ' h; will perform ', config % nt, ' steps'
            call logger % info('run_model', log_str)
        end if

        time_loop: do n = 1, config % nt
  
            if (this_image() == 1) then
                write (log_str, '(2(a,i6))') 'Computing time step', n, ' /', config % nt
                call logger % info('run_model', log_str)
            end if

            call sw_dynamics_step(config % dt)

            call halo_exchange()

            call write_output(n * config % dt)
        
        end do time_loop

    end subroutine run_model

end module mod_model