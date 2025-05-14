module mod_timing

   use mod_kinds, only: ik, rk
   use mod_log, only: logger => main_logger, log_str

   implicit none

   private
   public :: timing_on, timing_off, init_timing, print_timing

   integer(ik), parameter :: max_blocks = 100
   integer(ik), parameter :: max_label_length = 20
   integer(ik) :: n_blocks = 0

   integer(ik) :: n_times(max_blocks)

   real(rk) :: usr_total_time(max_blocks)
   real(rk) :: usr_start_time(max_blocks)

   real(rk) :: sys_total_time(max_blocks)
   real(rk) :: sys_start_time(max_blocks)

   character(len=max_label_length) :: labels(max_blocks)

contains

   subroutine init_timing()

      n_times(:) = 0
      usr_total_time(:) = 0
      sys_total_time(:) = 0
      usr_start_time(:) = -1
      sys_start_time(:) = -1

   end subroutine init_timing

   function get_index(label) result(idx)
      character(len=*), intent(in) :: label
      integer(ik) :: idx, i

      idx = -1
      do i = 1, n_blocks
         if (labels(i) == trim(label)) idx = i
      end do

      if (idx < 1) then
         ! We didn't find this label already, so add it
         n_blocks = n_blocks + 1
         labels(n_blocks) = trim(label)
         idx = n_blocks
      end if

   end function

   subroutine timing_on(label)
      character(len=*), intent(in) :: label
      integer(ik) :: idx
      real :: times(2), total_time

      call etime(times, total_time)

      idx = get_index(label)
      usr_start_time(idx) = times(1)
      sys_start_time(idx) = times(2)

   end subroutine timing_on

   subroutine timing_off(label)
      character(len=*), intent(in) :: label
      integer(ik) :: idx
      real :: times(2), total_time

      call etime(times, total_time)

      idx = get_index(label)

      if (usr_start_time(idx) < 0 .or. sys_start_time(idx) < 0) then
         write (log_str, '(a,a)') &
            'timing_off called before timing_on for ', trim(label)
         call logger%warning('timing_off', log_str)
         return
      end if

      n_times(idx) = n_times(idx) + 1
      usr_total_time(idx) = &
         usr_total_time(idx) + (times(1) - usr_start_time(idx))
      sys_total_time(idx) = &
         sys_total_time(idx) + (times(2) - sys_start_time(idx))

   end subroutine timing_off

   subroutine print_timing()
      integer(ik) :: i
      character(len=1) :: nl
      character(len=8) :: format_str

      nl = new_line('a')

      write (log_str, '(a)') 'Timing information.'//nl &
         //'Block                  n      User         Sys          Total'//nl &
         //'------------------------------------------------------------------'

      call logger%info('print_timing', log_str)

      do i = 1, n_blocks
         write (log_str, '(2a,i4)') adjustl(labels(i)), '   ', n_times(i)

         write (format_str, '(f8.3)') usr_total_time(i)
         write (log_str, '(a)') trim(log_str)//'   '//adjustl(format_str)//' s'

         write (format_str, '(f8.3)') sys_total_time(i)
         write (log_str, '(a)') trim(log_str)//'   '//adjustl(format_str)//' s'

         write (format_str, '(f8.3)') usr_total_time(i) + sys_total_time(i)
         write (log_str, '(a)') trim(log_str)//'   '//adjustl(format_str)//' s'

         call logger%info('print_timing', log_str, preamble=.false.)
      end do

   end subroutine print_timing

end module mod_timing
