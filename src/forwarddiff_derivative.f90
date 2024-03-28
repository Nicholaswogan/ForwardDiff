module forwarddiff_derivative
  use forwarddiff_const, only: wp
  use forwarddiff_dual, only: dual
  implicit none
  private

  public :: derivative, derivative_sig
  public :: grad, grad_sig

  abstract interface
    function derivative_sig(x) result(res)
      import :: dual
      type(dual), intent(in) :: x
      type(dual) :: res
    end function

    function grad_sig(x) result(res)
      import :: dual
      type(dual), intent(in) :: x(:)
      type(dual) :: res
    end function
  end interface

contains

  subroutine derivative(fcn, x, fcn_x, dfcn_x)
    procedure(derivative_sig) :: fcn
    real(wp), intent(in) :: x
    real(wp), intent(out) :: fcn_x
    real(wp), intent(out) :: dfcn_x
    type(dual) :: f
    f = fcn(dual(x, [1.0_wp]))
    fcn_x = f%val
    dfcn_x = f%der(1)
  end subroutine

  subroutine grad(fcn, x, fcn_x, dfcn_x)
    procedure(grad_sig) :: fcn
    real(wp), intent(in) :: x(:)
    real(wp), intent(out) :: fcn_x
    real(wp), intent(out) :: dfcn_x(:)

    type(dual) :: xx(size(x))
    type(dual) :: f
    integer :: i

    xx%val = x
    do i = 1,size(x)
      allocate(xx(i)%der(size(x)))
      xx(i)%der = 0.0_wp
      xx(i)%der(i) = 1.0_wp
    enddo

    f = fcn(xx)
    fcn_x = f%val
    dfcn_x = f%der(:)
  end subroutine




end module