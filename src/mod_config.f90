module mod_config

    use mod_kinds, only: ik, rk
    use mod_log, only: logger => main_logger, log_str
  
    implicit none
  
    private
    public :: main_config, init_config
  
    type :: Config
      integer(ik) :: nx, ny, nt
      real(rk)    :: dx, dy, dt
      real(rk) :: Lx, Ly
      real(rk) :: gravity, coriolis
    end type Config

    type(Config) :: main_config
  
  contains
  
    subroutine init_config(filename)
      character(len=*), intent(in) :: filename
      type(Config) :: conf
      integer(ik) :: fileunit, iostat
      integer(ik) :: nx, ny, nt
      real(rk)    :: dx, dy, dt
      real(rk) :: Lx, Ly
      real(rk) :: g, f

      namelist /domain/ nx, ny, nt, dt, Lx, Ly
      namelist /phys_param/ g, f

      open(newunit=fileunit, file=filename, status='old', action='read')

      read(fileunit, iostat=iostat, nml=domain)
      if (iostat /= 0) then
        write (log_str, '(a)') 'Could not read `domain` namelist'
        call logger % fatal('init_config', log_str)
      end if

      read(fileunit, iostat=iostat, nml=phys_param)
      if (iostat /= 0) then
        write (log_str, '(a)') 'Could not read `phys_param` namelist'
        call logger % fatal('init_config', log_str)
      end if

      close(fileunit)

      dx = Lx / (nx - 1)
      dy = Ly / (ny - 1)

      main_config = Config(nx, ny, nt, dx, dy, dt, Lx, Ly, g, f)
    end subroutine init_config
  
  end module mod_config