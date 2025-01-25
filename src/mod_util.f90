module mod_util

    use mod_kinds, only: ik, rk
    use mod_log, only: logger => main_logger, log_str

    implicit none

    private
    public :: abort, parse_duration
    
contains

    subroutine abort()

        close (unit=logger % output_unit)
        close (unit=logger % error_unit)
        stop 'Aborted'
        
    end subroutine abort

    function parse_duration(duration_str) result(duration)
        ! Read an ISO8601 formatted duration string and output an equivalent
        ! number of seconds. We assume
        !   1 year is exactly 365 days
        !   1 month is exactly 30 days
        !   1 week is exactly 7 days
        !   1 day is exactly 24 hours

        character(len=*), intent(in) :: duration_str
        real(rk) :: duration

        integer(ik) :: prev

        integer(ik) :: i_period, i_year, i_month, i_week, i_day
        integer(ik) :: i_time, i_hour, i_minute, i_second

        real(rk) :: n_year, n_month, n_week, n_day
        real(rk) :: n_hour, n_minute, n_second

        i_period = 0
        i_year   = 0
        i_month  = 0
        i_week   = 0
        i_day    = 0
        i_time   = 0
        i_hour   = 0
        i_minute = 0
        i_second = 0

        n_year   = 0.
        n_month  = 0.
        n_week   = 0.
        n_day    = 0.
        n_hour   = 0.
        n_minute = 0.
        n_second = 0.

        i_period = scan(duration_str, 'P')
        if (i_period /= 1) goto 10

        i_time = scan(duration_str, 'T')
        if (i_time > 0) then
            i_year  = scan(duration_str(:i_time), 'Y')
            i_month = scan(duration_str(:i_time), 'M')
            i_week  = scan(duration_str(:i_time), 'W')
            i_day   = scan(duration_str(:i_time), 'D')

            i_hour   = scan(duration_str(i_time+1:), 'H')
            i_minute = scan(duration_str(i_time+1:), 'M')
            i_second = scan(duration_str(i_time+1:), 'S')

            if (i_hour   /= 0) i_hour   = i_time + i_hour
            if (i_minute /= 0) i_minute = i_time + i_minute
            if (i_second /= 0) i_second = i_time + i_second
        else
            i_year  = scan(duration_str, 'Y')
            i_month = scan(duration_str, 'M')
            i_week  = scan(duration_str, 'W')
            i_day   = scan(duration_str, 'D')
        end if

        prev = i_period

        if (i_year /= 0) then
            n_year = read_real(duration_str(prev+1:i_year-1))
            if (n_year == -1) goto 10
            prev = i_year
        end if

        if (i_month /= 0) then
            n_month = read_real(duration_str(prev+1:i_month-1))
            if (n_month == -1) goto 10
            prev = i_month
        end if

        if (i_week /= 0) then
            n_week = read_real(duration_str(prev+1:i_week-1))
            if (n_week == -1) goto 10
            prev = i_week
        end if

        if (i_day /= 0) then
            n_day = read_real(duration_str(prev+1:i_day-1))
            if (n_day == -1) goto 10
            prev = i_day
        end if

        if (i_time /= 0) then
            prev = i_time

            if (i_hour /= 0) then
                n_hour = read_real(duration_str(prev+1:i_hour-1))
                if (n_hour == -1) goto 10
                prev = i_hour
            end if

            if (i_minute /= 0) then
                n_minute = read_real(duration_str(prev+1:i_minute-1))
                if (n_minute == -1) goto 10
                prev = i_minute
            end if

            if (i_second /= 0) then
                n_second = read_real(duration_str(prev+1:i_second-1))
                if (n_second == -1) goto 10
                prev = i_second
            end if
        end if

        duration = (               &
              31536000. * n_year   &
            +  2592000. * n_month  &
            +   604800. * n_week   &
            +    86400. * n_day    &
            +     3600. * n_hour   &
            +       60. * n_minute &
            +        1. * n_second )
        return

        10 write (log_str, '(a)') 'Could not read ISO8601 duration `' // duration_str // '`'
        call logger % fatal('parse_duration', log_str)
        call abort()

    end function parse_duration

    function read_real(s) result(res)
        character(len=*), intent(in) :: s
        character(len=10) :: format_str
        real(rk) :: res

        write (format_str, '(a,i2.2,a)') '(f', len(trim(s)), '.0)'
        read (s, trim(format_str), err=10) res
        return

        10 res = -1
        
    end function read_real
    
end module mod_util