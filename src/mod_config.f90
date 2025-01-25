module mod_config

    use mod_kinds, only: ik, rk
    use mod_log, only: logger => main_logger, log_str
    use mod_util, only: abort_now, parse_duration
  
    implicit none
  
    private
    public :: main_config, init_config
  
    type :: Config
      integer(ik) :: nx, ny
      real(rk) :: Lx, Ly
      real(rk) :: dx, dy
      real(rk) :: t_final, dt_max
      real(rk) :: dt_output
      real(rk) :: gravity, coriolis
    end type Config

    type(Config) :: main_config
  
  contains
  
    subroutine init_config(filename)
      character(len=*), intent(in) :: filename
      integer(ik) :: fileunit, iostat

      integer(ik) :: nx, ny
      character(len=99) :: max_timestep, run_duration, output_interval
      real(rk) :: dx, dy, t_final, dt_max, dt_output
      real(rk) :: Lx, Ly
      real(rk) :: g, f

      namelist /domain/ nx, ny, Lx, Ly
      namelist /time_control/ run_duration, output_interval, max_timestep
      namelist /phys_param/ g, f

      open(newunit=fileunit, file=filename, status='old', action='read')

      read(fileunit, iostat=iostat, nml=domain)
      if (iostat /= 0) then
        write (log_str, '(a)') 'Could not read `domain` namelist'
        call logger % fatal('init_config', log_str)
        call abort_now()
      end if

      read(fileunit, iostat=iostat, nml=time_control)
      if (iostat /= 0) then
        write (log_str, '(a)') 'Could not read `time_control` namelist'
        call logger % fatal('init_config', log_str)
        call abort_now()
      end if

      read(fileunit, iostat=iostat, nml=phys_param)
      if (iostat /= 0) then
        write (log_str, '(a)') 'Could not read `phys_param` namelist'
        call logger % fatal('init_config', log_str)
        call abort_now()
      end if

      close(fileunit)

      dx = Lx / nx
      dy = Ly / ny

      t_final = parse_duration(trim(run_duration))
      dt_max = parse_duration(trim(max_timestep))
      dt_output = parse_duration(trim(output_interval))

      main_config = Config(nx, ny, Lx, Ly, dx, dy, t_final, dt_max, dt_output, g, f)
    end subroutine init_config
  
  end module mod_config