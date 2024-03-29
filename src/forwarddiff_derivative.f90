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
      err = 'Output `dfdx` array is not the right size.'
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

  subroutine jacobian(fcn, x, f, dfdx, bandwidth, err)
    procedure(jacobian_sig) :: fcn
    real(wp), intent(in) :: x(:)
    real(wp), intent(out) :: f(:)
    real(wp), intent(out) :: dfdx(:,:)
    integer, optional, intent(in) :: bandwidth
    character(:), allocatable :: err

    ! real(wp), allocatable :: JacT(:,:) 
    integer :: bandwidth_, hbw
    type(dual) :: xx(size(x))
    type(dual) :: ff(size(x))
    integer :: i, j, ii, jj, kk

    bandwidth_ = size(x)
    if (present(bandwidth)) then
      bandwidth_ = bandwidth
      if (bandwidth_ > size(x)) then
        err = '`bandwidth` can not be > size(x).'
        return
      endif
      if (bandwidth_ < 1) then
        err = '`bandwidth` can not be < 1.'
        return
      endif
      if (mod(bandwidth_,2) == 0) then
        err = '`bandwidth` must be odd.'
        return
      endif
    endif

    if (size(x) /= size(f)) then
      err = 'Output `f` array is not the right size.'
      return
    endif
    if (bandwidth_ /= size(dfdx,1) .or. size(x) /= size(dfdx,2)) then
      err = 'Output `dfdx` array is not the right size.'
      return
    endif

    ! allocate(JacT(size(x),bandwidth_))

    ! Set x
    xx%val = x

    ! Allocate dual number
    do i = 1,size(x)
      allocate(xx(i)%der(bandwidth_))
      xx(i)%der = 0.0_wp
    enddo

    ! Seed the dual number
    j = 1
    outer : do
      do i = 1,bandwidth_
        if (j > size(x)) exit outer 
        xx(j)%der(i) = 1.0_wp
        j = j + 1
      enddo
    enddo outer

    do i = 1,size(x)
      print*,xx(i)%der(:)
    enddo
    print*,"hi"

    ! Do differentiation
    call fcn(xx, ff)

    ! Unpack f(x)
    f = ff%val

    ! Unpack jacobian
    ! do j = 1,bandwidth_
    !   do i = 1,size(x)
    !     JacT(i,j) = ff(i)%der(j)
    !     ! dfdx(j,i) = ff(i)%der(j)
    !     ! dfdx(1,:) is the first column of the jacobian
    !   enddo
    ! enddo
    ! JacT(:,1) is the first column of the jacobian

    ! We now want to disentangle the jacobian
    ! dfdx(1,:) is a row. It should contain diagonals. Row 1 is the highest diagonal

    ! First row of dfdx
    ! dfdx(1,1:hbw) = 0.0_dp, unused
    ! JacT(1,hbw+1) should be at dfdx(1,hbw+1)
    ! JacT(2,hbw+2) should be at dfdx(1,hbw+2)
    ! JacT(3,hbw) should be at dfdx(1,hbw+3)
    ! JacT(4,hbw+1) should be at dfdx(1,hbw+4)
    ! JacT(5,hbw+2) should be at dfdx(1,hbw+5)
    ! JacT(6,hbw-1) should be at dfdx(1,hbw+6)

    ! Second row

    ! JacT(1,bw) is first diagonal entry

    if (present(bandwidth)) then

      hbw = (bandwidth_ - 1)/2 ! halfbandwidth 
      j = 1
      outer1 : do
        do i = 1,bandwidth_
          if (j > size(x)) exit outer1
          do jj = -hbw,hbw
            kk = jj + hbw + 1
            ii = j + jj
            if (ii < 1) then
              dfdx(kk,j) = 0.0_wp
            elseif (ii > size(x)) then
              dfdx(kk,j) = 0.0_wp
            else
              ! dfdx(kk,j) = JacT(ii,i)
              dfdx(kk,j) = ff(ii)%der(i)
            endif
          enddo
          j = j + 1
        enddo
      enddo outer1

    else
      
      ! Dense matrix
      do j = 1,bandwidth_
        do i = 1,size(x)
          dfdx(i,j) = ff(i)%der(j)
        enddo
      enddo
    
    endif

  end subroutine

end module