program test_forwarddif
  use forwarddif, only: wp, derivative
  implicit none
  call test()

contains

  subroutine test()
    real(wp) :: fcn_x, dfcn_x
    real(wp) :: t(2)

    call cpu_time(t(1))
    call derivative(func_operators, 10.0_wp, fcn_x, dfcn_x)
    call cpu_time(t(2))
    print*,fcn_x, dfcn_x, (t(2) - t(1))

    call cpu_time(t(1))
    call derivative(func_intrinsics, 10.0_wp, fcn_x, dfcn_x)
    call cpu_time(t(2))
    print*,fcn_x, dfcn_x, (t(2) - t(1))
    
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

  function func_intrinsics(x) result(res)
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

  

end program