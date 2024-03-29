
program main
  use forwarddiff, only: wp, jacobian
  implicit none
  integer, parameter :: nz = 6

  ! call test_banded()
  call test_blockdiagonal1()
  call test_blockdiagonal2()

contains

  subroutine test_banded()
    use forwarddiff, only: BandedJacobian
    integer, parameter :: bandwidth = 3
    real(wp) :: u(nz), f(nz), dfdu(bandwidth,nz)
    real(wp) :: f1(nz), dfdu1(nz,nz)
    character(:), allocatable :: err
    integer :: i
    
    do i = 1,nz
      u(i) = i
    enddo
    call jacobian(rhs_banded_dual, u, f, dfdu, jt=BandedJacobian, bandwidth=bandwidth, err=err)
    if (allocated(err)) then
      print*,err
      stop
    endif

    call rhs_banded(u, f1)
    call jac_banded(u, dfdu1)

    do i = 1,bandwidth
      print*,dfdu(i,:)
    enddo
    print*,''

    do i = 1,nz
      print*,dfdu1(i,:)
    enddo

  end subroutine

  subroutine rhs_banded_dual(u, du)
    use forwarddiff
    type(dual), intent(in) :: u(:)
    type(dual), intent(out) :: du(:)
    integer :: i

    du(1) = 3*u(2) - u(1)
    do i = 2,nz-1
      du(i) = 3*u(i+1) - 2.0_wp*u(i) + u(i-1)
    enddo
    du(nz) = - u(nz) + u(nz-1)

  end subroutine

  subroutine rhs_banded(u, du)
    real(wp), intent(in) :: u(:)
    real(wp), intent(out) :: du(:)
    integer :: i

    du(1) = 3*u(2) - u(1)
    do i = 2,nz-1
      du(i) = 3*u(i+1) - 2.0_wp*u(i) + u(i-1)
    enddo
    du(nz) = - u(nz) + u(nz-1)

  end subroutine
  
  subroutine jac_banded(u, pd)
    real(wp), intent(in) :: u(:)
    real(wp), intent(out) :: pd(:,:)
    integer :: i

    pd = 0.0_wp
    pd(1,1) = -1.0_wp
    pd(2,1) = 1.0_wp
    do i = 2,nz-1
      pd(i,i) = -2.0_wp
      pd(i+1,i) = 1.0_wp
      pd(i-1,i) = 3.0_wp
    enddo
    pd(nz,nz) = -1.0_wp
    pd(nz-1,nz) = 3.0_wp
    
  end subroutine

  subroutine test_blockdiagonal1()
    use forwarddiff, only: BlockDiagonalJacobian
    integer, parameter :: blocksize = 2
    real(wp) :: u(nz), f(nz), dfdu(blocksize,nz)
    real(wp) :: dfdu1(nz,nz)
    character(:), allocatable :: err
    integer :: i

    do i = 1,nz
      u(i) = i
    enddo
    call jacobian(rhs_blocked1_dual, u, f, dfdu, jt=BlockDiagonalJacobian, blocksize=blocksize, err=err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    print*,''
    do i = 1,blocksize
      print*,dfdu(i,:) 
    enddo

    call jacobian(rhs_blocked1_dual, u, f, dfdu1, err=err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    print*,''
    do i = 1,nz
      print*,dfdu1(i,:) 
    enddo

  end subroutine

  subroutine test_blockdiagonal2()
    use forwarddiff, only: BlockDiagonalJacobian
    integer, parameter :: blocksize = 3
    real(wp) :: u(nz), f(nz), dfdu(blocksize,nz)
    real(wp) :: dfdu1(nz,nz)
    character(:), allocatable :: err
    integer :: i

    do i = 1,nz
      u(i) = i
    enddo
    call jacobian(rhs_blocked2_dual, u, f, dfdu, jt=BlockDiagonalJacobian, blocksize=blocksize, err=err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    print*,''
    do i = 1,blocksize
      print*,dfdu(i,:) 
    enddo

    call jacobian(rhs_blocked2_dual, u, f, dfdu1, err=err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    print*,''
    do i = 1,nz
      print*,dfdu1(i,:) 
    enddo

  end subroutine

  subroutine rhs_blocked1_dual(u, du)
    use forwarddiff
    type(dual), intent(in) :: u(:)
    type(dual), intent(out) :: du(:)
    integer :: i

    du(1) = u(1) + u(2)*u(1)
    du(2) = u(2) + u(2)*u(1)

    du(3) = u(3)*u(4) + u(3)
    du(4) = u(4) + u(3)*u(3)

    du(5) = 2.0_wp*u(5) + 0.5_wp*u(6)
    du(6) = u(5) + u(6)*u(6)

  end subroutine

  subroutine rhs_blocked2_dual(u, du)
    use forwarddiff
    type(dual), intent(in) :: u(:)
    type(dual), intent(out) :: du(:)
    integer :: i

    du(1) = u(1) + u(2)*u(1)
    du(2) = u(2) + u(2)*u(1)
    du(3) = u(3) + u(3)*u(1)

    du(4) = u(4) + u(5)*u(6)
    du(5) = 2.0_wp*u(5) + 0.5_wp*u(6)
    du(6) = u(5) + u(6)*u(6)

  end subroutine

end program