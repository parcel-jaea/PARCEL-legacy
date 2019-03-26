!
!     Copyright (C) 1998-1999 Japan Atomic Energy Research Institute.
!
! ==========================================================
!>
!>@file
!> matrix and vector elements for the random matrix
!>
! DESCRIPTION:
!>  <B>  set matrix A and vector B, X ( A*X = b ) 
!>   in the case of the random matrix </B>
!>
!>@brief     matrix A : <B> double </B> precision, 
!>           vector b : <B> double </B> precision,
!>           vector X : <B> double </B> precision 
!>
!>@param[in]  ndim
!>@param[out] matA, vecB, X (initial value), X_exact (exact solution)
!
! ==========================================================
! **********************************************************
subroutine set_matrix_random( NNZR, ndim, nzero, matA, RowIndex, ColIndex, &
     &                        vecB, X, X_exact, nzero_out, &
     &                        istart, iend, matA_is_symmetric )
! **********************************************************
!                                  !!   outputs are in the coordinate system 
  implicit none
!--------------
! inout
!--------------
  integer, intent(in) :: NNZR, ndim, nzero
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
  integer, allocatable :: indrc(:), nzrow(:), iwork(:)
  real*8, allocatable :: work(:)
!
  integer :: seed = 12357

!--------------
! begin
!--------------
  allocate( indrc( NNZR ) )
  allocate( work( nzero ) )
  allocate( iwork( nzero ) )
  allocate( nzrow( ndim ) )

  call genmxran2( nzero, ndim, NNZR, nzero_out, &
       &          matA, RowIndex, ColIndex, X_exact, vecB, &
       &          nzrow, indrc, work, iwork, istart, iend, &
       &          matA_is_symmetric, seed )

  deallocate( indrc, work, iwork, nzrow )

  X = 1.0

end subroutine set_matrix_random

