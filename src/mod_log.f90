module mod_log

   use iso_fortran_env, only: stdout => output_unit, stderr => error_unit
   use mod_kinds, only: ik

   implicit none

   private
   public :: init_logging, main_logger, log_str, &
             log_fatal, log_error, log_warning, log_info, log_debug

   type :: Logger

      integer(ik) :: output_level, error_level
      integer(ik) :: output_unit, error_unit

   contains

      procedure, pass(self)         :: handle_message
      procedure, public, pass(self) :: fatal
      procedure, public, pass(self) :: error
      procedure, public, pass(self) :: warning
      procedure, public, pass(self) :: info
      procedure, public, pass(self) :: debug

   end type

   type(Logger) :: main_logger
   character(len=999) :: log_str

   integer(ik), parameter :: log_fatal = 40
   integer(ik), parameter :: log_error = 30
   integer(ik), parameter :: log_warning = 20
   integer(ik), parameter :: log_info = 10
   integer(ik), parameter :: log_debug = 0

contains

   subroutine init_logging(output_file, error_file, output_level, error_level)
      character(*), intent(in) :: output_file, error_file
      integer(ik), intent(in), optional :: output_level
      integer(ik), intent(in), optional :: error_level
      integer(ik) :: level_1, level_2
      integer(ik) :: output_unit, error_unit

      if (present(output_level)) then
         level_1 = output_level
      else
         level_1 = log_info
      end if

      if (present(error_level)) then
         level_2 = error_level
      else
         level_2 = log_error
      end if

      open (newunit=output_unit, file=output_file, status='replace')
      open (newunit=error_unit, file=error_file, status='replace')

      main_logger = Logger(level_1, level_2, &
                           output_unit, error_unit)

   end subroutine init_logging

   subroutine handle_message(self, level, message, preamble)
      class(Logger), intent(in out) :: self
      character(*), intent(in) :: message
      integer(ik), intent(in) :: level
      logical, intent(in) :: preamble
      integer(ik) :: datetime(8)

      call date_and_time(values=datetime)

      if (preamble) then
         write (log_str, '(a,i4,2(a,i2.2),3(a,i2.2),a,i3.3,2a)') &
            '[', datetime(1), '-', datetime(2), '-', datetime(3), &
            ' ', datetime(5), ':', datetime(6), ':', datetime(7), &
            '.', datetime(8), '] ', trim(message)
      else
         write (log_str, '(a)') trim(message)
      end if

      if (level >= self%error_level) then
         write (self%error_unit, '(a)') trim(log_str)
         flush (self%error_unit)
      else
         write (self%output_unit, '(a)') trim(log_str)
         flush (self%output_unit)
      end if

      if (level >= self%output_level) then
         if (level >= self%error_level) then
            write (stderr, '(a)') trim(log_str)
            flush (stderr)
         else
            write (stdout, '(a)') trim(log_str)
            flush (stdout)
         end if
      end if

   end subroutine handle_message

   subroutine fatal(self, source, message, preamble)
      class(Logger), intent(in out) :: self
      logical, intent(in), optional :: preamble
      character(*), intent(in) :: source, message
      logical :: print_preamble

      if (present(preamble) .and. .not. preamble) then
         call self%handle_message(log_fatal, message, .false.)
      else
         call self%handle_message( &
            log_fatal, '['//source//'] <fatal> '//trim(message), .true.)
      end if

   end subroutine fatal

   subroutine error(self, source, message, preamble)
      class(Logger), intent(in out) :: self
      logical, intent(in), optional :: preamble
      character(*), intent(in) :: source, message
      logical :: print_preamble

      if (present(preamble) .and. .not. preamble) then
         call self%handle_message(log_error, message, .false.)
      else
         call self%handle_message( &
            log_error, '['//source//'] <error> '//trim(message), .true.)
      end if

   end subroutine error

   subroutine warning(self, source, message, preamble)
      class(Logger), intent(in out) :: self
      logical, intent(in), optional :: preamble
      character(*), intent(in) :: source, message

      if (present(preamble) .and. .not. preamble) then
         call self%handle_message(log_warning, message, .false.)
      else
         call self%handle_message( &
            log_warning, '['//source//'] <warning> '//trim(message), .true.)
      end if

   end subroutine warning

   subroutine info(self, source, message, preamble)
      class(Logger), intent(in out) :: self
      logical, intent(in), optional :: preamble
      character(*), intent(in) :: source, message

      if (present(preamble) .and. .not. preamble) then
         call self%handle_message(log_info, message, .false.)
      else
         call self%handle_message( &
            log_info, '['//source//'] <info> '//trim(message), .true.)
      end if

   end subroutine info

   subroutine debug(self, source, message, preamble)
      class(Logger), intent(in out) :: self
      logical, intent(in), optional :: preamble
      character(*), intent(in) :: source, message

      if (present(preamble) .and. .not. preamble) then
         call self%handle_message(log_debug, message, .false.)
      else
         call self%handle_message( &
            log_debug, '['//source//'] <debug> '//trim(message), .true.)
      end if

   end subroutine debug

end module mod_log
