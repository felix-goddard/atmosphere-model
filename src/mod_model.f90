module mod_model

    use mod_kinds, only: ik, rk
    use mod_log, only: logger => main_logger, log_str
    use mod_config, only: config => main_config
    use mod_field, only: Field, diffx, diffy

    implicit none

    private
    public :: init_model, run_model

    type(Field) :: h, u, v
contains
    
    subroutine init_model()
        integer(ik) :: i, j
        integer(ik), parameter :: ic = 51, jc = 51
        real(rk), parameter :: decay = 0.02

        u = Field('u', [config % nx, config % ny])
        v = Field('v', [config % nx, config % ny])
        h = Field('h', [config % nx, config % ny])

        ! initialize a gaussian blob in the center
        do concurrent(i = h % lb(1)-1:h % ub(1)+1,&
            j = h % lb(2)-1:h % ub(2)+1)
            h % data(i, j) = exp(-decay * ((i - ic)**2 + (j - jc)**2))
        end do
        call h % sync_edges()

    end subroutine init_model

    subroutine run_model()
        integer(ik) :: n
        real(rk) :: dt, dx, dy
        real(rk) :: g, hm

        dt = config % dt
        dx = config % Lx / (config % nx - 1)
        dy = config % Ly / (config % ny - 1)
        g = config % gravity
        hm = config % mean_depth

        time_loop: do n = 1, config % nt
  
            if (this_image() == 1) then
                write (log_str, '(2(a,i6))') 'Computing time step', n, ' /', config % nt
                call logger % info('main', log_str)
            end if
        
            u = u - (u * diffx(u) / dx + v * diffy(u) / dy &
                + g * diffx(h) / dx) * dt
            call u % sync_edges()
        
            v = v - (u * diffx(v) / dx + v * diffy(v) / dy &
                + g * diffy(h) / dy) * dt
            call v % sync_edges()
        
            h = h - (diffx(u * (hm + h)) / dx + diffy(v * (hm + h)) / dy) * dt
            call h % sync_edges()
        
            call h % write(n)
        
        end do time_loop

    end subroutine run_model

end module mod_model