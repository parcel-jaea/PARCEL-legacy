!
!     Copyright (C) 1998-1999 Japan Atomic Energy Research Institute.
!
! **********************************************************
subroutine set_matrix_reduced_Frank( NNZR, ndim, nzero, matA, RowIndex, ColIndex, &
     &                               vecB, X, X_exact, nzero_out, &
     &                               istart, iend, matA_is_symmetric )
! **********************************************************
!                                  !!   outputs are in the coordinate system 
  implicit none
!--------------
! inout
!--------------
  integer, intent(in) :: ndim, nzero, NNZR
  integer, intent(in) :: istart, iend
  logical, intent(in) :: matA_is_symmetric
!
  integer, intent(out) :: nzero_out
!
  integer, intent(out) :: RowIndex( nzero ), ColIndex( nzero )
  real*8, intent(out) :: matA( nzero )
  real*8, intent(out) :: X(ndim), X_exact(ndim), vecB(ndim)

!--------------
! local
!--------------
  integer i, j, ix, iy, ip, ierr
  integer n1, i_min, i_max
  real*8 h, z, ctmp1, ctmp2, value
  real*8 alpha, beta, gamma, delta
!--------------
! begin
!--------------
  ip = 1

  n1 = ( NNZR -1 )/2
!  write(*,*) 'n1 = ', n1

  do j=istart, iend
     i_min = j -n1;  i_max = j +n1

     if ( i_min <= 0 ) i_min = 1
     if ( i_max >= ndim ) i_max = ndim

     do i=i_min, i_max
        if ( matA_is_symmetric .and. i<j ) cycle

        value = ndim + 1 - max(i,j)
        call pushco( matA, RowIndex, ColIndex, ip, value, j, i, nzero )
     end do
  end do

  nzero_out = ip - 1

! --------------
! set (global) rhs vector
! --------------
  alpha = 0.7870; beta = 0.923; gamma = -0.347; delta = 0.613

  Do i=1, ndim
     value = ndim +1 -max(i,1)
     if ( abs(i-1) > n1 ) value = 0.0
     VecB(i) = VecB(i) + alpha *value

     value = ndim +1 -max(i,2)
     if ( abs(i-2) > n1 ) value = 0.0
     VecB(i) = VecB(i) + beta *value

     value = ndim +1 -max(i,3)
     if ( abs(i-3) > n1 ) value = 0.0
     VecB(i) = VecB(i) + gamma *value

     value = ndim +1 -max(i,4)
     if ( abs(i-4) > n1 ) value = 0.0
     VecB(i) = VecB(i) + delta *value
  End do

! --------------
! boundary conditions 
! --------------

! --------------
! set exact solution
! --------------
  X_exact(1) = alpha;  X_exact(2) = beta;  X_exact(3) = gamma;  X_exact(4) = delta

! --------------
! set itital guess
! --------------
  X(:) = 0.4;  X(1) = 0.7;  X(2) = -0.3
!  do i=1, ndim
!     x(i) = i-1
!  end do

contains

  subroutine pushco( val_co, row_co, col_co, ip, q, ir, ic, nnz0 )
    integer, intent(in) ::  nnz0, ir, ic
    integer, intent(inout) ::  ip
    integer, intent(out) :: row_co(nnz0), col_co(nnz0)
    real*8, intent(in) :: q
    real*8, intent(out) :: val_co(nnz0)

    val_co(ip) = q
    row_co(ip) = ir
    col_co(ip) = ic
    ip = ip + 1
  end subroutine pushco

end subroutine set_matrix_reduced_Frank
