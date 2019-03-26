! directly gives DIA formmat
! **********************************************************

subroutine set_matrix_poisson_2D_ddm2( mm, ncomp_grids, lnum_grid_axis,  &
     &                                 lsta_grid_axis, &
     &                                 margin, nnd, mdist, &
     &                                 ioff_grids, ioff_mat, &
     &                                 matA, vecB, vecX, vecX_exact,&
     &                                 rank_x, rank_y, &
     &                                 npe_x, npe_y )
! **********************************************************
!                                  !!   outputs are in the DIA
  implicit none
!--------------
! inout
!--------------
  integer, intent(in) :: mm, ncomp_grids, margin, nnd, mdist
  integer, intent(in) :: rank_x, rank_y, npe_x, npe_y
  integer, intent(in) :: lnum_grid_axis( ncomp_grids )
  integer, intent(in) :: lsta_grid_axis( ncomp_grids )
  integer, intent(out) :: ioff_grids( nnd, ncomp_grids )
  integer, intent(out) :: ioff_mat( nnd )
  real*8, intent(out) :: matA( nnd*mdist )
  real*8, intent(out) :: vecX(mdist), vecX_exact(mdist), vecB(mdist)

!--------------
! local
!--------------
  integer nx, ny, nz, mx, my, mz, m, bw
  integer jx, jy, jz, ix, iy, iz, j, num
  real*8 h, z1, ctmp1, ctmp2, ctmp3


!--------------
! init
!--------------
  if ( ncomp_grids /= 2 ) stop "ncomp_grids is not 2"
  bw = margin

  nx = lnum_grid_axis( 1 )   ! local
  ny = lnum_grid_axis( 2 )

  mx = nx +bw *2;     my = ny +bw *2; 
  m = mx *my

  matA = 0.0d0;  vecB = 0.0d0;  vecX = 0.0d0;  vecX_exact = 0.0d0

  ioff_grids = 0;
  ioff_mat = 0
!--------------
! begin
!--------------
  Do jy=0, ny-1
     Do jx=0, nx-1
        j = ( jx +bw ) +mx *( jy +bw ) +1;
           
        num = 0;

        if ( ( rank_y == 0 .and. jy > 0 ) .or. ( rank_y > 0 ) ) then
           matA( j +num*m ) = 1.0d0
        endif
        num = num +1;
        
        if ( ( rank_x == 0 .and. jx > 0 ) .or. ( rank_x > 0 ) ) then
           matA( j +num*m ) = 1.0d0
        endif
        num = num +1;

        matA( j +num*m ) = -4.0d0        ! diagonal
        num = num +1;

        if ( ( rank_x == npe_x -1 .and. jx < nx-1 ) .or. ( rank_x < npe_x -1 ) ) then
           matA( j +num*m ) = 1.0d0
!           matA( j +num*m ) = 0.9d0
        endif
        num = num +1;
        
        if ( ( rank_y == npe_y -1 .and. jy < ny-1 ) .or. ( rank_y < npe_y -1 ) ) then
           matA( j +num*m ) = 1.0d0
!           matA( j +num*m ) = 1.1d0
        endif
        num = num +1;

     End do
  End Do

! --------------
! set offset
! --------------
#if 1
  ioff_grids(1,1) = 0;    ioff_grids(1,2) =-1;
  ioff_grids(2,1) =-1;    ioff_grids(2,2) = 0;
  ioff_grids(3,1) = 0;    ioff_grids(3,2) = 0;
  ioff_grids(4,1) = 1;    ioff_grids(4,2) = 0;
  ioff_grids(5,1) = 0;    ioff_grids(5,2) = 1;
#else
  ioff_grids(1,1) = 0;    ioff_grids(1,2) =-1;
  ioff_grids(2,1) =-1;    ioff_grids(2,2) = 0;

  ioff_grids(3,1) = 0;    ioff_grids(3,2) = 1;
  ioff_grids(4,1) = 1;    ioff_grids(4,2) = 0;
  ioff_grids(5,1) = 0;    ioff_grids(5,2) = 0;
#endif
  Do j=1, nnd
     ioff_mat(j) = mx *ioff_grids(j,2) &
          &      + ioff_grids(j,1)
  End Do

! --------------
! set (local) rhs vector
! --------------
  h = 1.0D0 / (mm+1);   z1 = 2.0d0 * h**4

  Do jy=0, ny-1
     Do jx=0, nx-1
        ix = jx +lsta_grid_axis(1) -1
        iy = jy +lsta_grid_axis(2) -1
        j = ( jx +bw ) +mx *( jy +bw ) +1;

        ctmp1 = ( ix+1 )**2;   ctmp2 = (iy+1)**2;
!        vecB(j) = z1 *( ctmp1 +ctmp2 )
        vecB(j) = -1.0d0
     End do
  End Do

! --------------
! boundary conditions
! --------------
  if ( rank_x == npe_x -1 ) then
     jx = nx -1
     Do jy=0, ny-1
        ix = jx +lsta_grid_axis( 1 ) -1
        iy = jy +lsta_grid_axis( 2 ) -1

        j = ( jx +bw ) +mx *( jy +bw ) +1;
!        vecB(j) = vecB(j) -( h*(iy+1) )**2
        vecB(j) = -1.0d0
     End Do
  endif

  if ( rank_y == npe_y -1 ) then
     jy = ny -1
     Do jx=0, nx -1
        ix = jx +lsta_grid_axis( 1 ) -1
        iy = jy +lsta_grid_axis( 2 ) -1

        j = ( jx +bw ) +mx *( jy +bw )  +1;
!        vecB(j) = vecB(j) -( h*(ix+1) )**2
        vecB(j) = -1.0d0
     End Do
  endif

! --------------
! set exact solution
! --------------
  z1 = h**4

  Do jy=0, ny -1
     Do jx=0, nx -1
        ix = jx +lsta_grid_axis( 1 ) -1
        iy = jy +lsta_grid_axis( 2 ) -1

        j = ( jx +bw ) +mx *( jy +bw ) +1;

        ctmp1 = ( ix+1 )**2;   ctmp2 = (iy+1)**2;
        vecx_exact(j) = z1 *( ctmp1*ctmp2 );
     End Do
  End Do

! --------------
! set itital guess
! --------------
  vecx = 0.0d0
  if ( rank_x ==0 .and. rank_y == 0 ) then
     jx = 0;    jy = 0
     j = ( jx +bw ) +mx *( jy +bw ) +1;
!     vecx(j) = 1.0d0
     vecx(j) = 0.0d0
  endif

end subroutine set_matrix_poisson_2D_ddm2

subroutine set_matrix_poisson_2D_9stencil_ddm2( mm, ncomp_grids, lnum_grid_axis,  &
     &                                 lsta_grid_axis, &
     &                                 margin, nnd, mdist, &
     &                                 ioff_grids, ioff_mat, &
     &                                 matA, vecB, vecX, vecX_exact,&
     &                                 rank_x, rank_y, &
     &                                 npe_x, npe_y )
! **********************************************************
!                                  !!   outputs are in the DIA
  implicit none
!--------------
! inout
!--------------
  integer, intent(in) :: mm, ncomp_grids, margin, nnd, mdist
  integer, intent(in) :: rank_x, rank_y, npe_x, npe_y
  integer, intent(in) :: lnum_grid_axis( ncomp_grids )
  integer, intent(in) :: lsta_grid_axis( ncomp_grids )
  integer, intent(out) :: ioff_grids( nnd, ncomp_grids )
  integer, intent(out) :: ioff_mat( nnd )
  real*8, intent(out) :: matA( nnd*mdist )
  real*8, intent(out) :: vecX(mdist), vecX_exact(mdist), vecB(mdist)

!--------------
! local
!--------------
  integer nx, ny, nz, mx, my, mz, m, bw
  integer jx, jy, jz, ix, iy, iz, j, num
  real*8 h, z1, ctmp1, ctmp2, ctmp3


!--------------
! init
!--------------
  if ( ncomp_grids /= 2 ) stop "ncomp_grids is not 2"
  bw = margin

  nx = lnum_grid_axis( 1 )   ! local
  ny = lnum_grid_axis( 2 )

  mx = nx +bw *2;     my = ny +bw *2; 
  m = mx *my

  matA = 0.0d0;  vecB = 0.0d0;  vecX = 0.0d0;  vecX_exact = 0.0d0

  ioff_grids = 0;
  ioff_mat = 0
!--------------
! begin
!--------------
  Do jy=0, ny-1
     Do jx=0, nx-1
        j = ( jx +bw ) +mx *( jy +bw ) +1;
           
        num = 0;

        if ( ( rank_y == 0 .and. jy > 1 ) .or. ( rank_y > 0 ) ) then
           matA( j +num*m ) = -1.0d0/12.0d0
        endif
        num = num +1;
        if ( ( rank_y == 0 .and. jy > 0 ) .or. ( rank_y > 0 ) ) then
           matA( j +num*m ) = 4.0d0/3.0d0
        endif
        num = num +1;        
        if ( ( rank_x == 0 .and. jx > 1 ) .or. ( rank_x > 0 ) ) then
           matA( j +num*m ) = -1.0d0/12.0d0
        endif
        num = num +1;
        if ( ( rank_x == 0 .and. jx > 0 ) .or. ( rank_x > 0 ) ) then
           matA( j +num*m ) = 4.0d0/3.0d0
        endif
        num = num +1;

        matA( j +num*m ) = -5.0d0        ! diagonal
        num = num +1;

        if ( ( rank_x == npe_x -1 .and. jx < nx-1 ) .or. ( rank_x < npe_x -1 ) ) then
           matA( j +num*m ) = 4.0d0/3.0d0
        endif
        num = num +1;
        if ( ( rank_x == npe_x -1 .and. jx < nx-2 ) .or. ( rank_x < npe_x -1 ) ) then
           matA( j +num*m ) = -1.0d0/12.0d0
        endif
        num = num +1;
        
        if ( ( rank_y == npe_y -1 .and. jy < ny-1 ) .or. ( rank_y < npe_y -1 ) ) then
           matA( j +num*m ) = 4.0d0/3.0d0
        endif
        num = num +1;
        if ( ( rank_y == npe_y -1 .and. jy < ny-2 ) .or. ( rank_y < npe_y -1 ) ) then
           matA( j +num*m ) = -1.0d0/12.0d0
        endif
        num = num +1;

     End do
  End Do

! --------------
! set offset
! --------------
#if 0
  ioff_grids(1,1) = 0;    ioff_grids(1,2) =-1;
  ioff_grids(2,1) =-1;    ioff_grids(2,2) = 0;
  ioff_grids(3,1) = 0;    ioff_grids(3,2) = 0;
  ioff_grids(4,1) = 1;    ioff_grids(4,2) = 0;
  ioff_grids(5,1) = 0;    ioff_grids(5,2) = 1;
#else
  ioff_grids(1,1) = 0;    ioff_grids(1,2) =-2;
  ioff_grids(2,1) = 0;    ioff_grids(2,2) =-1;
  ioff_grids(3,1) =-2;    ioff_grids(3,2) = 0;
  ioff_grids(4,1) =-1;    ioff_grids(4,2) = 0;
  ioff_grids(5,1) = 0;    ioff_grids(5,2) = 0;
  ioff_grids(6,1) = 1;    ioff_grids(6,2) = 0;
  ioff_grids(7,1) = 2;    ioff_grids(7,2) = 0;
  ioff_grids(8,1) = 0;    ioff_grids(8,2) = 1;
  ioff_grids(9,1) = 0;    ioff_grids(9,2) = 2;
#endif
  Do j=1, nnd
     ioff_mat(j) = mx *ioff_grids(j,2) &
          &      + ioff_grids(j,1)
  End Do

! --------------
! set (local) rhs vector
! --------------
  h = 1.0D0 / (mm+1);   z1 = 2.0d0 * h**4

  Do jy=0, ny-1
     Do jx=0, nx-1
        ix = jx +lsta_grid_axis(1) -1
        iy = jy +lsta_grid_axis(2) -1
        j = ( jx +bw ) +mx *( jy +bw ) +1;

        ctmp1 = ( ix+1 )**2;   ctmp2 = (iy+1)**2;
!        vecB(j) = z1 *( ctmp1 +ctmp2 )
        vecB(j) = 1.0d0
     End do
  End Do

! --------------
! boundary conditions
! --------------
  if ( rank_x == npe_x -1 ) then
     jx = nx -1
     Do jy=0, ny-1
        ix = jx +lsta_grid_axis( 1 ) -1
        iy = jy +lsta_grid_axis( 2 ) -1

        j = ( jx +bw ) +mx *( jy +bw ) +1;
!        vecB(j) = vecB(j) -( h*(iy+1) )**2
        vecB(j) =1.0d0
     End Do
  endif

  if ( rank_y == npe_y -1 ) then
     jy = ny -1
     Do jx=0, nx -1
        ix = jx +lsta_grid_axis( 1 ) -1
        iy = jy +lsta_grid_axis( 2 ) -1

        j = ( jx +bw ) +mx *( jy +bw )  +1;
!        vecB(j) = vecB(j) -( h*(ix+1) )**2
        vecB(j) = 1.0d0
     End Do
  endif

! --------------
! set exact solution
! --------------
  z1 = h**4

  Do jy=0, ny -1
     Do jx=0, nx -1
        ix = jx +lsta_grid_axis( 1 ) -1
        iy = jy +lsta_grid_axis( 2 ) -1

        j = ( jx +bw ) +mx *( jy +bw ) +1;

        ctmp1 = ( ix+1 )**2;   ctmp2 = (iy+1)**2;
        vecx_exact(j) = z1 *( ctmp1*ctmp2 );
     End Do
  End Do

! --------------
! set itital guess
! --------------
  vecx = 0.0d0
  if ( rank_x ==0 .and. rank_y == 0 ) then
     jx = 0;    jy = 0
     j = ( jx +bw ) +mx *( jy +bw ) +1;
!     vecx(j) = 1.0d0
     vecx(j) = 0.0d0
  endif

end subroutine set_matrix_poisson_2D_9stencil_ddm2
