!
!     Copyright (C) 1998-1999 Japan Atomic Energy Research Institute.
!
! **********************************************************
subroutine set_matrix_poisson_3D_lm( mm, ndim, nzero, matA, RowIndex, ColIndex, &
     &                            vecB, X, X_exact, nzero_out, &
     &                            istart, iend, matA_is_symmetric )
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
!  real*8, intent(out) :: X(ndim), X_exact(ndim), vecB(ndim)
  real*8, intent(out) :: X(istart:iend)
  real*8, intent(out) :: X_exact(istart:iend)
  real*8, intent(out) :: vecB(istart:iend)

!--------------
! local
!--------------
  integer i, ix, iy, iz, ip, ierr
  real*8 h, z, ctmp1, ctmp2, ctmp3

!--------------
! begin
!--------------
  ip = 1

  do i=istart, iend
     ix = mod( i-1, mm ) +1
     iy = ( i -1 ) /mm +1
     iy = mod( iy -1, mm ) +1
     iz = ( i-1 )/mm**2 +1

!     write(*,*) i, ix, iy, iz
     if ( .not. matA_is_symmetric ) then
        if ( iz > 1 ) then
           call pushco( matA, RowIndex, ColIndex, ip, 1.0D0, i, i-mm**2, nzero )
        endif
        if ( iy > 1 ) then
           call pushco( matA, RowIndex, ColIndex, ip, 1.0D0, i, i-mm, nzero )
        endif
        if ( ix > 1 ) then
           call pushco( matA, RowIndex, ColIndex, ip, 1.0D0, i, i-1, nzero )
        endif
     endif

     call pushco( matA, RowIndex, ColIndex, ip, -6.0D0, i, i, nzero )
     if ( ix < mm ) then
        call pushco( matA, RowIndex, ColIndex, ip, 1.0D0, i, i+1, nzero )
     endif
     if ( iy < mm ) then
        call pushco( matA, RowIndex, ColIndex, ip, 1.0D0, i, i+mm, nzero )
     endif
     if ( iz < mm ) then
        call pushco( matA, RowIndex, ColIndex, ip, 1.0D0, i, i+mm**2, nzero )
     endif
  end do
  nzero_out = ip - 1
!  stop

! --------------
! set (global) rhs vector
! --------------
  h = 1.0D0 / (mm+1);   z = 2.0d0*h**6
  do iz=1, mm
     do iy=1, mm
        do ix=1, mm
           i = mm**2*(iz-1) +mm*(iy-1) + ix
           ctmp1 = real(ix) **2;    ctmp2 = real(iy) **2;  ctmp3 = real(iz)**2
           if (i < istart .or. iend < i) cycle ! lm
           vecB(i) = z *( ctmp1*ctmp2 +ctmp2*ctmp3 +ctmp3*ctmp1 )
        end do
     end do
  end do

! --------------
! boundary conditions 
! --------------
  ix = mm
  do iz=1, mm
     do iy=1, mm
        i = mm**2 *(iz-1) +mm*(iy-1) +ix
        if (i < istart .or. iend < i) cycle ! lm
        vecB(i) = vecB(i) -(h*iy)**2 *(h*iz)**2
     end do
  end do

  iy = mm
  do iz=1, mm
     do ix=1, mm
        i = mm**2 *(iz-1) +mm*(iy-1) +ix
        if (i < istart .or. iend < i) cycle ! lm
        vecB(i) = vecB(i) - (h*ix)**2 * (h*iz)**2
     end do
  end do

  iz = mm
  do iy=1, mm
     do ix=1, mm
        i = mm**2 *(iz-1) +mm*(iy-1) +ix
        if (i < istart .or. iend < i) cycle ! lm
        vecB(i) = vecB(i) - (h*ix)**2 * (h*iy)**2
     end do
  end do

! --------------
! set exact solution
! --------------
  z = h**6
  do iz=1, mm
     do iy=1, mm
        do ix=1, mm
           i = mm**2 *(iz-1) +mm *(iy-1) +ix
           if (i < istart .or. iend < i) cycle ! lm
           ctmp1 = real(ix) **2;    ctmp2 = real(iy) **2;  ctmp3 = real(iz)**2
           x_exact(i) = z * ( ctmp1 *ctmp2 *ctmp3 )
        end do
     end do
  end do

! --------------
! set itital guess
! --------------
  x = 0.0d0
  if (1 >= istart .and. iend >= 1) then
     X(1) = 1.0
  end if

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

end subroutine set_matrix_poisson_3D_lm
