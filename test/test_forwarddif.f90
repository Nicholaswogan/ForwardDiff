program test_forwarddif
  use forwarddif, only: wp, derivative, grad
  implicit none
  call test()

contains

  subroutine test()
    real(wp) :: x, f, dfdx
    real(wp) :: xx(2), dfdx1(2)

    x = 10.0_wp
    call derivative(func_operators, x, f, dfdx)
    print*,f, dfdx

    x = 10.0_wp
    call derivative(func_intrinsics1, x, f, dfdx)
    print*,f, dfdx

    x = 0.1_wp
    call derivative(func_intrinsics2, x, f, dfdx)
    print*,f, dfdx

    xx = [1.0_wp, 2.0_wp]
    call grad(func_grad1, xx, f, dfdx1)
    print*,f, dfdx1

    xx = [1.0_wp, 2.0_wp]
    call grad(func_grad2, xx, f, dfdx1)
    print*,f, dfdx1

  end subroutine

  function func_operators(x) result(res)
    use forwarddif_dual
    type(dual), intent(in) :: x
    type(dual) :: res

    res = x + x + 3.0_wp + x
    res = res - x - 3.0_wp - x
    res = res*x + 2.0_wp*x + x*(-5.0_wp)
    res = res/x + 2.0_wp/x + x/5.0_wp
    res = res**x + res**1.5_wp

  end function

  function func_intrinsics1(x) result(res)
    use forwarddif_dual
    type(dual), intent(in) :: x
    type(dual) :: res

    res = abs(x)
    res = res + cos(x)
    res = res + exp(x)
    res = res + log(x)
    res = res + log10(x)
    res = res + sin(x)
    res = res + tan(x)
    res = res + sqrt(x)

  end function

  function func_intrinsics2(x) result(res)
    use forwarddif_dual
    type(dual), intent(in) :: x
    type(dual) :: res

    res = acos(x)
    res = res + asin(x)
    res = res + atan(x)
    res = max(res, x)
    res = max(res, 1.0_wp)
    res = max(1.0_wp, res)
    res = min(res, res)
    res = min(res, 2.0_wp)
    res = min(2.0_wp, res)

  end function

  function func_grad1(x) result(res)
    use forwarddif_dual
    type(dual), intent(in) :: x(:)
    type(dual) :: res
    res = x(1)*x(1)*x(2) + x(1) + x(2)
  end function

  function func_grad2(x) result(res)
    use forwarddif_dual
    type(dual), intent(in) :: x(:)
    type(dual) :: res
    res = sum(x*3.14_wp)
  end function

end program