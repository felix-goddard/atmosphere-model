module mod_config

    use mod_kinds, only: ik, rk
    use mod_log, only: logger => main_logger, log_str
    use mod_util, only: abort_now, parse_duration

    implicit none

    private
    public :: main_config, init_config

    type :: Config
        ! io control
        logical :: save_restart_file
        character(len=:), allocatable :: initial_filename, restart_filename

        ! time control
        real(rk) :: dt_max, t_final, dt_output

        ! domain parameters
        integer(ik) :: nx, ny
        real(rk) :: Lx, Ly, dx, dy

        ! physics parameters
        real(rk) :: gravity, coriolis
    end type Config

    type(Config) :: main_config

contains

    subroutine init_config(filename)
        character(len=*), intent(in) :: filename
        integer(ik) :: fileunit, iostat
        real(rk) :: dx, dy, t_final, dt_max, dt_output

        logical :: save_restart_file
        character(len=99) :: initial_filename, restart_filename
        namelist /io_control/ initial_filename, save_restart_file, restart_filename

        character(len=99) :: max_timestep, run_duration, output_interval
        namelist /time_control/ max_timestep, run_duration, output_interval

        integer(ik) :: nx, ny
        real(rk) :: Lx, Ly
        namelist /domain_parameters/ nx, ny, Lx, Ly

        real(rk) :: g, f
        namelist /physics_parameters/ g, f

        open(newunit=fileunit, file=filename, status='old', action='read')

        read(fileunit, iostat=iostat, nml=io_control)
        if (iostat /= 0) then
            write (log_str, '(a)') 'Could not read `io_control` namelist'
            call logger % fatal('init_config', log_str)
            call abort_now()
        end if

        read(fileunit, iostat=iostat, nml=time_control)
        if (iostat /= 0) then
            write (log_str, '(a)') 'Could not read `time_control` namelist'
            call logger % fatal('init_config', log_str)
            call abort_now()
        end if

        read(fileunit, iostat=iostat, nml=domain_parameters)
        if (iostat /= 0) then
            write (log_str, '(a)') 'Could not read `domain_parameters` namelist'
            call logger % fatal('init_config', log_str)
            call abort_now()
        end if

        read(fileunit, iostat=iostat, nml=physics_parameters)
        if (iostat /= 0) then
            write (log_str, '(a)') 'Could not read `physics_parameters` namelist'
            call logger % fatal('init_config', log_str)
            call abort_now()
        end if

        close(fileunit)

        dx = Lx / nx
        dy = Ly / ny

        t_final = parse_duration(trim(run_duration))
        dt_max = parse_duration(trim(max_timestep))
        dt_output = parse_duration(trim(output_interval))

        main_config = Config(                                      &
            save_restart_file, initial_filename, restart_filename, &
            dt_max, t_final, dt_output,                            &
            nx, ny, Lx, Ly, dx, dy,                                &
            g, f                                                   )

    end subroutine init_config

end module mod_config
