!
!     Copyright (C) 1998-1999 Japan Atomic Energy Research Institute.
!
! **********************************************************
subroutine set_matrix_glued_wilkinson( blocksize, nrepeat, glue,  &
     &                                 ndim, nzero, matA, RowIndex, ColIndex, &
     &                                 vecB, X, X_exact, nzero_out, &
     &                                 istart, iend, matA_is_symmetric )
! **********************************************************
!                                  !!   outputs are in the coordinate system 
  implicit none
!--------------
! inout
!--------------
  integer, intent(in) :: blocksize, nrepeat, ndim, nzero
  integer, intent(in) :: istart, iend
  logical, intent(in) :: matA_is_symmetric
  real*8, intent(in) :: glue
!
  integer, intent(out) :: nzero_out
!
  integer, intent(out) :: RowIndex( nzero ), ColIndex( nzero )
  real*8, intent(out) :: matA( nzero )
  real*8, intent(out) :: X(ndim), X_exact(ndim), vecB(ndim)
!--------------
! local
!--------------
  integer ib, ii, jb, jj, ip, i, j
  real*8 ctmp1, alpha, beta, gamma, delta
  real*8 ycenter
  real*8, parameter :: eta = 1.0E-8

!--------------
! begin
!--------------
  ip = 1

  ycenter = ( blocksize -1 )/2

  do i=istart, iend
! --------
! address of the block
! ------    
     ib = ( i -1 ) /blocksize +1

! --------
! position inside of the block
! --------
     ii = mod( i-1, blocksize ) +1

     if ( .not. matA_is_symmetric ) then
        j = i -1;  jb = ( j-1 )/blocksize +1;  jj = mod( j-1,blocksize ) +1
        if ( j >= 1 ) then
           if ( jb == ib ) then
              call pushco( matA, RowIndex, ColIndex, ip, 1.0d0, i, j, nzero )
           else if ( jb == ib -1 ) then
              call pushco( matA, RowIndex, ColIndex, ip, glue, i, j, nzero )
           endif
        endif
     endif

     ctmp1 = abs( ycenter -( ii -1 ) )
     if ( abs(ctmp1) < eta ) ctmp1 = eta

     call pushco( matA, RowIndex, ColIndex, ip, ctmp1, i, i, nzero )

     j = i +1;  jb = ( j-1 )/blocksize +1;  jj = mod( j-1,blocksize ) +1
     if ( j <= ndim ) then
        if ( jb == ib ) then
           call pushco( matA, RowIndex, ColIndex, ip, 1.0d0, i, j, nzero )
        else if ( jb == ib +1 ) then
           call pushco( matA, RowIndex, ColIndex, ip, glue, i, j, nzero )
        endif
     endif
  end do

  nzero_out = ip - 1

! --------------
! set (global) rhs vector
! --------------
  alpha = 0.7870; beta = 0.923; gamma = -0.347; delta = 0.613
!  write(*,*) 'ycenter ', ycenter

  vecB(1) = alpha *ycenter + beta*1.0
  vecB(2) = alpha *1.0  + beta*( ycenter-1.0 ) +gamma *1.0
  vecB(3) = beta *1.0 + gamma *( ycenter-2.0 ) +delta *1.0
  vecB(4) = gamma *1.0 + delta *( ycenter-3.0 ) 
  vecB(5) = delta *1.0

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
  
end subroutine set_matrix_glued_wilkinson
