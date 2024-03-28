module mymod
  use iso_fortran_env, only: dp => real64
  implicit none
  public

  type :: dual
    real(dp) :: val
    real(dp), allocatable :: der(:)
  end type
  ! interface dual
  !   module procedure :: create_dual_1
  ! end interface

  interface operator (*)
    module procedure :: mul_dd
  end interface
  
  interface operator (+)
    module procedure :: add_dd
  end interface

  interface sin
    module procedure :: sin_d
  end interface

contains

  pure function create_dual_1(val, der) result(b)
    real(dp), intent(in) :: val, der(:)
    type(dual) :: b
    b%val = val
    b%der = der
  end function

  pure elemental function mul_dd(a1, a2) result(b)
    type(dual), intent(in) :: a1
    type(dual), intent(in) :: a2
    type(dual) :: b
    b%val = a1%val*a2%val
    b%der = a1%val*a2%der + a1%der*a2%val
  end function

  pure elemental function add_dd(a1, a2) result(b)
    type(dual), intent(in) :: a1
    type(dual), intent(in) :: a2
    type(dual) :: b
    b%val = a1%val + a2%val
    b%der = a2%der + a1%der
  end function

  pure elemental function sin_d(a) result(b)
    type(dual), intent(in) :: a
    type(dual) :: b
    b%val = sin(a%val)
    b%der = a%der*cos(a%val)
  end function

end module


program test1
  use mymod, only: dp
  implicit none
  call test()

contains
  subroutine test()
    use mymod
    real(dp) :: a, b
    type(dual) :: xx, yy, f(2), x(2)

    a = 1.0_dp
    b = 2.0_dp

    xx = dual(a, [1.0_dp, 0.0_dp])
    yy = dual(b, [0.0_dp, 1.0_dp])
    x(1) = xx
    x(2) = yy

 
    call fcn1(x,f)

    print*,f(1)%val
    print*,f(2)%val

    print*,''
    print*,f(1)%der
    print*,f(2)%der

  end subroutine

  ! function fcn(x, y) result(res)
  !   use mymod
  !   type(dual), intent(in) :: x, y
  !   type(dual) :: res
  !   res = x*x*y + x + y
  ! end function

  subroutine fcn1(xx, res)
    use mymod
    type(dual), intent(in) :: xx(:)
    type(dual) :: res(:)
    type(dual) :: x, y
    x = xx(1)
    y = xx(2)
    res(1) = x*x + y*y*sin(y)
    res(2) = x + y*sin(x)
  end subroutine


end program