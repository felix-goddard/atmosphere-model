module mod_config

    use mod_kinds, only: ik, rk
    use mod_log, only: logger => main_logger, log_str
    use mod_util, only: parse_duration
  
    implicit none
  
    private
    public :: main_config, init_config
  
    type :: Config
      integer(ik) :: nx, ny
      real(rk) :: t_final, dt_max
      real(rk) :: dx, dy
      real(rk) :: Lx, Ly
      real(rk) :: gravity, coriolis
    end type Config

    type(Config) :: main_config
  
  contains
  
    subroutine init_config(filename)
      character(len=*), intent(in) :: filename
      integer(ik) :: fileunit, iostat

      type(Config) :: conf
      integer(ik) :: nx, ny
      character(len=99) :: run_duration
      real(rk) :: dx, dy, dt_max, t_final
      real(rk) :: Lx, Ly
      real(rk) :: g, f

      namelist /domain/ nx, ny, Lx, Ly
      namelist /time_control/ run_duration, dt_max
      namelist /phys_param/ g, f

      open(newunit=fileunit, file=filename, status='old', action='read')

      read(fileunit, iostat=iostat, nml=domain)
      if (iostat /= 0) then
        write (log_str, '(a)') 'Could not read `domain` namelist'
        call logger % fatal('init_config', log_str)
      end if

      read(fileunit, iostat=iostat, nml=time_control)
      if (iostat /= 0) then
        write (log_str, '(a)') 'Could not read `time_control` namelist'
        call logger % fatal('init_config', log_str)
      end if

      read(fileunit, iostat=iostat, nml=phys_param)
      if (iostat /= 0) then
        write (log_str, '(a)') 'Could not read `phys_param` namelist'
        call logger % fatal('init_config', log_str)
      end if

      close(fileunit)

      dx = Lx / nx
      dy = Ly / ny
      t_final = parse_duration(trim(run_duration))

      main_config = Config(nx, ny, t_final, dt_max, dx, dy, Lx, Ly, g, f)
    end subroutine init_config
  
  end module mod_config