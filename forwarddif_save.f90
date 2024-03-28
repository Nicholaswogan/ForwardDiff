module forwarddif
  use iso_fortran_env, only: dp => real64
  implicit none
  private

  public :: derivative, fcn_sig
  
  public :: dp
  public :: dual
  public :: sin, cos
  public :: operator (**)
  public :: operator (*), operator (/)
  public :: operator (+)
  
  type :: dual
    real(dp) :: val
    real(dp) :: der
  end type
  
  interface dual
    module procedure :: create_dual_1
  end interface
  
  interface sin
    module procedure :: sin_d
  end interface
  
  interface cos
    module procedure :: cos_d
  end interface
  
  interface operator (**)
    module procedure :: pow_dr
  end interface
  
  interface operator (*)
    module procedure :: mul_dd
    module procedure :: mul_dr
    module procedure :: mul_rd
  end interface
  
  interface operator (/)
    module procedure :: div_dd
    module procedure :: div_dr
    module procedure :: div_rd
  end interface
  
  interface operator (+)
    module procedure :: add_dd
    module procedure :: add_dr
    module procedure :: add_rd
  end interface
  
  ! interface operator (-)
  !   module procedure :: sub_dd
  !   module procedure :: sub_dr
  !   module procedure :: sub_rd
  ! end interface

  abstract interface
    function fcn_sig(x) result(res)
      import :: dual
      type(dual), intent(in) :: x
      type(dual) :: res
    end function
  end interface
  
contains

  subroutine derivative(fcn, x, fcn_x, dfcn_x)
    procedure(fcn_sig) :: fcn
    real(dp), intent(in) :: x
    real(dp), intent(out) :: fcn_x
    real(dp), intent(out) :: dfcn_x
    type(dual) :: f
    f = fcn(dual(x))
    fcn_x = f%val
    dfcn_x = f%der
  end subroutine

  pure elemental function create_dual_1(val) result(b)
    real(dp), intent(in) :: val
    type(dual) :: b
    b%val = val
    b%der = 1.0_dp
  end function
  
  pure elemental function sin_d(a) result(b)
    type(dual), intent(in) :: a
    type(dual) :: b
    b%val = sin(a%val)
    b%der = a%der*cos(a%val)
  end function
  
  pure elemental function cos_d(a) result(b)
    type(dual), intent(in) :: a
    type(dual) :: b
    b%val = cos(a%val)
    b%der = -a%der*sin(a%val)
  end function
  
  pure elemental function pow_dr(a1, a2) result(b)
    type(dual), intent(in) :: a1
    real(dp), intent(in) :: a2
    type(dual) :: b
    b%val = a1%val**a2
    b%der = a2*(a1%val**(a2-1))*a1%der
  end function
  
  pure elemental function mul_dd(a1, a2) result(b)
    type(dual), intent(in) :: a1
    type(dual), intent(in) :: a2
    type(dual) :: b
    b%val = a1%val*a2%val
    b%der = a1%val*a2%der + a1%der*a2%val
  end function
  
  pure elemental function mul_dr(a1, a2) result(b)
    type(dual), intent(in) :: a1
    real(dp), intent(in) :: a2
    type(dual) :: b
    b%val = a1%val*a2
    b%der = a1%der*a2
  end function
  
  pure elemental function mul_rd(a1, a2) result(b)
    real(dp), intent(in) :: a1
    type(dual), intent(in) :: a2
    type(dual) :: b
    b%val = a1*a2%val
    b%der = a1*a2%der
  end function
  
  pure elemental function div_dd(a1, a2) result(b)
    type(dual), intent(in) :: a1
    type(dual), intent(in) :: a2
    type(dual) :: b
    b%val = a1%val/a2%val
    b%der = (a1%val*a2%der - a1%der*a2%val)/(a2%val**2)
  end function
  
  pure elemental function div_dr(a1, a2) result(b)
    type(dual), intent(in) :: a1
    real(dp), intent(in) :: a2
    type(dual) :: b
    b%val = a1%val/a2
    b%der = a1%der/a2
  end function
  
  pure elemental function div_rd(a1, a2) result(b)
    real(dp), intent(in) :: a1
    type(dual), intent(in) :: a2
    type(dual) :: b
    b%val = a1/a2%val
    b%der = a1/a2%der
  end function
  
  pure elemental function add_dd(a1, a2) result(b)
    type(dual), intent(in) :: a1
    type(dual), intent(in) :: a2
    type(dual) :: b
    b%val = a1%val + a2%val
    b%der = a2%der + a1%der
  end function
  
  pure elemental function add_dr(a1, a2) result(b)
    type(dual), intent(in) :: a1
    real(dp), intent(in) :: a2
    type(dual) :: b
    b%val = a1%val + a2
    b%der = a1%der
  end function
  
  pure elemental function add_rd(a1, a2) result(b)
    real(dp), intent(in) :: a1
    type(dual), intent(in) :: a2
    type(dual) :: b
    b = a2 + a1
  end function
  
end module