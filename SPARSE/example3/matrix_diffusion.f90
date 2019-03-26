!
!     Copyright (C) 1998-1999 Japan Atomic Energy Research Institute.
!
! **********************************************************
subroutine set_matrix_diffusion_2D( mm, ndim, nzero, matA, RowIndex, ColIndex, &
     &                              vecB, X, X_exact, nzero_out, &
     &                              istart, iend, matA_is_symmetric )
! **********************************************************
!                                  !!   outputs are in the coordinate system 
  implicit none
!--------------
! inout
!--------------
  integer, intent(in) :: mm, ndim, nzero
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
  integer i, ix, iy, ip, ierr
  real*8 h, z, ctmp1, ctmp2, dtmp1
!
  real*8, parameter ::  PI = 3.1415926535D0

!--------------
! init
!--------------
  h = 1.0D0 / (mm+1);   z = h**2

!--------------
! begin
!--------------
  ip = 1

  do i=istart, iend
     ix = mod( i-1, mm ) +1
     iy = ( i -1 ) /mm +1

     if ( .not. matA_is_symmetric ) then
        if ( iy > 1 ) then
           call pushco( matA, RowIndex, ColIndex, ip, -1.0D0, i, i-mm, nzero )
        endif
        if ( ix > 1 ) then
           call pushco( matA, RowIndex, ColIndex, ip, -1.0D0, i, i-1, nzero )

        endif
     endif

     call pushco( matA, RowIndex, ColIndex, ip, 4.0D0+z, i, i, nzero )
     if ( ix < mm ) then
        call pushco( matA, RowIndex, ColIndex, ip, -1.0D0, i, i+1, nzero )
     endif
     if ( iy < mm ) then
        call pushco( matA, RowIndex, ColIndex, ip, -1.0D0, i, i+mm, nzero )
     endif
  end do
  nzero_out = ip - 1

! --------------
! set (global) rhs vector
! --------------
  h = 1.0D0 / (mm+1);   z = h**2
  do iy=1, mm
     do ix=1, mm
        i = mm*(iy-1) + ix
        ctmp1 = sin( PI*ix*h );  ctmp2 = sin( PI*iy*h )
        vecB(i) = z * ctmp1 *ctmp2
     end do
  end do

! --------------
! boundary conditions 
! --------------

! --------------
! set exact solution
! --------------
  z = h**2

  do iy=1, mm
     do ix=1, mm
        i = mm *(iy-1) +ix
        ctmp1 = sin( PI*ix*h );  ctmp2 = sin( PI*iy*h )
        dtmp1 = 1.0d0 + 8.0d0 /z *sin( PI*h/2.0d0 )**2
        x_exact(i) = ctmp1 *ctmp2 /dtmp1
     end do
  end do

! --------------
! set itital guess
! --------------
  x(1:ndim) = 1.0;   x(1) = 0.0

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

end subroutine set_matrix_diffusion_2D
