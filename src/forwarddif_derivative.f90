module forwarddif_derivative
  use forwarddif_const, only: wp
  use forwarddif_dual, only: dual
  implicit none
  private

  public :: derivative, derivative_sig

  abstract interface
    function derivative_sig(x) result(res)
      import :: dual
      type(dual), intent(in) :: x
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

end module