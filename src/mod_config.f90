module mod_config

   use mod_kinds, only: ik, rk
   use mod_log, only: logger => main_logger, log_str
   use mod_constants, only: set_constants
   use mod_util, only: abort_now, parse_duration

   implicit none

   private
   public :: main_config, read_config_file

   type :: Config
      ! io control
      logical :: save_restart_file
      character(len=:), allocatable :: initial_filename, restart_filename

      ! time control
      real(rk) :: dt_max, dt_output, dt_radiation, dt_physics, dt_remap, &
                  t_initial, t_final

      ! domain parameters
      integer(ik) :: nx, ny, nlay
      real(rk) :: Lx, Ly, dx, dy
   end type Config

   type(Config) :: main_config

contains

   subroutine read_config_file(filename)
      character(len=*), intent(in) :: filename
      integer(ik) :: fileunit, iostat
      real(rk) :: t_final, dt_max, dt_output, dt_radiation, dt_physics, dt_remap

      logical :: save_restart_file
      character(len=99) :: initial_filename, restart_filename
      namelist /io_control/ initial_filename, save_restart_file, &
         restart_filename

      character(len=99) :: max_timestep, run_duration, output_interval, &
                           radiation_interval, physics_interval, remap_interval
      namelist /time_control/ max_timestep, run_duration, output_interval, &
         radiation_interval, physics_interval, remap_interval

      real(rk) :: top_pressure, gravity, coriolis_parameter, &
                  divergence_damping_coefficient
      namelist /physics_parameters/ &
         top_pressure, gravity, coriolis_parameter, &
         divergence_damping_coefficient

      open (newunit=fileunit, file=filename, status='old', action='read')

      read (fileunit, iostat=iostat, nml=io_control)
      if (iostat /= 0) then
         write (log_str, '(a)') 'Could not read `io_control` namelist'
         call logger%fatal('init_config', log_str)
         call abort_now()
      end if

      read (fileunit, iostat=iostat, nml=time_control)
      if (iostat /= 0) then
         write (log_str, '(a)') 'Could not read `time_control` namelist'
         call logger%fatal('init_config', log_str)
         call abort_now()
      end if

      read (fileunit, iostat=iostat, nml=physics_parameters)
      if (iostat /= 0) then
         write (log_str, '(a)') 'Could not read `physics_parameters` namelist'
         call logger%fatal('init_config', log_str)
         call abort_now()
      end if

      close (fileunit)

      t_final = parse_duration(trim(run_duration))
      dt_max = parse_duration(trim(max_timestep))
      dt_output = parse_duration(trim(output_interval))
      dt_radiation = parse_duration(trim(radiation_interval))
      dt_physics = parse_duration(trim(physics_interval))
      dt_remap = parse_duration(trim(remap_interval))

      main_config = Config( &
                    save_restart_file, initial_filename, restart_filename, &
                    dt_max, dt_output, dt_radiation, dt_physics, dt_remap, &
                    0, t_final, -1, -1, -1, -1, -1, -1, -1)

      call set_constants(top_pressure, gravity, coriolis_parameter, &
                         divergence_damping_coefficient)

   end subroutine read_config_file

end module mod_config
