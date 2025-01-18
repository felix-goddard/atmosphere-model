module mod_config

    use mod_kinds, only: ik, rk
  
    implicit none
  
    private
    public :: main_config, init_config
  
    type :: Config
      integer(ik) :: nx, ny, nt
      real(rk) :: dt, Lx, Ly
      real(rk) :: gravity, mean_depth
    end type Config

    type(Config) :: main_config
  
  contains
  
    subroutine init_config(filename)
      character(len=*), intent(in) :: filename
      type(Config) :: conf
      integer(ik) :: fileunit
      integer(ik) :: nx, ny, nt
      real(rk) :: dt, Lx, Ly
      real(rk) :: g, hm

      namelist /domain/ nx, ny, nt, dt, Lx, Ly
      namelist /phys_param/ g, hm

      open(newunit=fileunit, file=filename, status='old', action='read')
      read(fileunit, nml=domain)
      read(fileunit, nml=phys_param)
      close(fileunit)

      main_config = Config(nx, ny, nt, dt, Lx, Ly, g, hm)
    end subroutine init_config
  
  end module mod_config