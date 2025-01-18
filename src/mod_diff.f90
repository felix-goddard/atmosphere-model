module mod_diff
    use mod_kinds, only: ik, rk

    implicit none

    private
    public :: diff_centered
    public :: diffx
    public :: diffy
    
contains

pure function diff_centered(x) result(dx)
    real(rk), intent(in) :: x(:)
    real(rk) :: dx(size(x))
    integer(ik) :: im

    im = size(x)
    dx(1)      = 0.5 * (x(2) - x(im))
    dx(im)     = 0.5 * (x(1) - x(im-1))
    dx(2:im-1) = 0.5 * (x(3:im) - x(1:im-2))
end function diff_centered

pure function diffx(x) result(dx)
    real(rk), intent(in) :: x(:,:)
    real(rk) :: dx(size(x, dim=1), size(x, dim=2))
    integer(ik) :: i

    i = size(x, dim=1)
    dx = 0
    dx(2:i-1,:) = 0.5 * (x(3:i,:) - x(1:i-2,:))
end function diffx

pure function diffy(x) result(dx)
    real(rk), intent(in) :: x(:,:)
    real(rk) :: dx(size(x, dim=1), size(x, dim=2))
    integer(ik) :: i

    i = size(x, dim=2)
    dx = 0
    dx(:,2:i-1) = 0.5 * (x(:,3:i) - x(:,1:i-2))
end function diffy
    
end module mod_diff