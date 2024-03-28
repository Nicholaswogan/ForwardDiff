program test_forwarddif
  use dnadmod, only: dual, dp => wp
  implicit none
  abstract interface
    function fcn_normal_sig(x) result(res)
      import :: dp
      real(dp), intent(in) :: x
      real(dp) :: res
    end function
  end interface

  call test()

contains

  subroutine test()
    use dnadmod, only: derivative
    implicit none
    real(dp) :: x, fcn_x, dfcn_x
    integer :: i, n
    real(dp) :: t(2), tmp, res(2)

    x = 2.0_dp
    n = 10000000

    ! Autodiff
    tmp = 0.0_dp
    call cpu_time(t(1))
    do i = 1,n
      call derivative(func_dual, x, fcn_x, dfcn_x)
      tmp = tmp + fcn_x
    enddo
    call cpu_time(t(2))
    print*,fcn_x, dfcn_x, (t(2) - t(1))/real(n,dp), tmp

    n = 100000000

    tmp = 0.0_dp
    call cpu_time(t(1))
    do i = 1,n
      fcn_x = func_normal(x)
      dfcn_x = dfunc_normal(x)
      tmp = tmp + fcn_x
    enddo
    call cpu_time(t(2))
    print*,fcn_x, dfcn_x, (t(2) - t(1))/real(n,dp), tmp

    tmp = 0.0_dp
    call cpu_time(t(1))
    do i = 1,n
      call finite_difference(func_normal, x, 1.0e-4_dp, fcn_x, dfcn_x)
      tmp = tmp + fcn_x
    enddo
    call cpu_time(t(2))
    print*,fcn_x, dfcn_x, (t(2) - t(1))/real(n,dp), tmp

  end subroutine

  subroutine finite_difference(fcn, x, eps, fcn_x, dfcn_x)
    use dnadmod, only: fcn_sig
    procedure(fcn_normal_sig) :: fcn
    real(dp), intent(in) :: x
    real(dp), intent(in) :: eps
    real(dp), intent(out) :: fcn_x
    real(dp), intent(out) :: dfcn_x
    real(dp) :: h
    fcn_x = fcn(x)
    h = x*eps
    dfcn_x = (fcn(x+h) - fcn_x)/h
  end subroutine

  function func_normal(x) result(res)
    real(dp), intent(in) :: x
    real(dp) :: res
    res = sin(x + x**2.0_dp)
  end function 

  function func_dual(x) result(res)
    use dnadmod
    type(dual), intent(in) :: x
    type(dual) :: res
    res = sin(x + x**2.0_dp)
  end function

  function dfunc_normal(x) result(res)
    real(dp), intent(in) :: x
    real(dp) :: res
    res = (2.0_dp*x + 1.0_dp)*cos(x + x**2.0_dp)
  end function
  
end program