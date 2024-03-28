module forwarddiff_derivative
  use forwarddiff_const, only: wp
  use forwarddiff_dual, only: dual
  implicit none
  private

  public :: derivative, derivative_sig
  public :: grad, grad_sig
  public :: jacobian, jacobian_sig

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

    subroutine jacobian_sig(x, f)
      import :: dual
      type(dual), intent(in) :: x(:)
      type(dual), intent(out) :: f(:)
    end subroutine
  end interface

contains

  subroutine derivative(fcn, x, f, dfdx)
    procedure(derivative_sig) :: fcn
    real(wp), intent(in) :: x
    real(wp), intent(out) :: f
    real(wp), intent(out) :: dfdx
    type(dual) :: ff
    ff = fcn(dual(x, [1.0_wp]))
    f = ff%val
    dfdx = ff%der(1)
  end subroutine

  subroutine grad(fcn, x, f, dfdx, err)
    procedure(grad_sig) :: fcn
    real(wp), intent(in) :: x(:)
    real(wp), intent(out) :: f
    real(wp), intent(out) :: dfdx(:)
    character(:), allocatable :: err

    type(dual) :: xx(size(x))
    type(dual) :: ff
    integer :: i

    if (size(x) /= size(dfdx)) then
      err = 'Output dfdx array is not the right size.'
      return
    endif

    ! Set x
    xx%val = x

    ! Seed the dual number
    do i = 1,size(x)
      allocate(xx(i)%der(size(x)))
      xx(i)%der = 0.0_wp
      xx(i)%der(i) = 1.0_wp
    enddo

    ! Do differentiation
    ff = fcn(xx)

    ! Unpack f(x)
    f = ff%val

    ! Unpack gradient
    dfdx(:) = ff%der(:)

  end subroutine

  subroutine jacobian(fcn, x, f, dfdx, err)
    procedure(jacobian_sig) :: fcn
    real(wp), intent(in) :: x(:)
    real(wp), intent(out) :: f(:)
    real(wp), intent(out) :: dfdx(:,:)
    character(:), allocatable :: err

    type(dual) :: xx(size(x))
    type(dual) :: ff(size(x))
    integer :: i

    if (size(x) /= size(f)) then
      err = 'Output f array is not the right size.'
      return
    endif
    if (size(x) /= size(dfdx,1) .or. size(x) /= size(dfdx,2)) then
      err = 'Output dfdx array is not the right size.'
      return
    endif

    xx%val = x
    do i = 1,size(x)
      allocate(xx(i)%der(size(x)))
      xx(i)%der = 0.0_wp
      xx(i)%der(i) = 1.0_wp
    enddo

    ! Do differentiation
    call fcn(xx, ff)

    ! Unpack f(x)
    f = ff%val
    
    ! Unpack jacobian
    do i = 1,size(x)
      dfdx(i,:) = ff(i)%der(:)
    enddo

  end subroutine




end module