module ddm_setup
  implicit none
contains
  subroutine set_neighbor_ranks( npe_x, npe_y, npe_z, rank_x, rank_y, rank_z, &
                                num_neighbor_ranks, &
                                neighbor_ranks, max_num_neighbor_ranks )
    implicit none
    ! inout
    integer, intent(in) :: npe_x, npe_y, npe_z, max_num_neighbor_ranks
    integer, intent(in) :: rank_x, rank_y, rank_z
    integer, intent(out) :: num_neighbor_ranks
    integer, intent(out) :: neighbor_ranks( max_num_neighbor_ranks )
    ! local
    integer ixmax, iymax, izmax
    integer ix, iy, iz, num
    integer n1x, n1y, n1z
    integer iirank

    ! begin
    ixmax = 0;  iymax = 0; izmax = 0;
    num = 0;

    neighbor_ranks = -99

    if ( npe_x > 1 ) ixmax=1
    if ( npe_y > 1 ) iymax=1
    if ( npe_z > 1 ) izmax=1

    Do iz=-izmax, izmax
       Do iy=-iymax, iymax
          Do ix=-ixmax, ixmax
             n1x = rank_x +ix;
             n1y = rank_y +iy;
             n1z = rank_z +iz;

             if ( ix*ix +iy*iy + iz*iz == 0 ) cycle
             num = num +1

             if ( n1z >=0 .and. n1z < npe_z ) then
                if ( n1y >=0 .and. n1y < npe_y ) then
                   if ( n1x >=0 .and. n1x < npe_x ) then
                      iirank = npe_x *npe_y *n1z +npe_x *n1y +n1x;
                      neighbor_ranks(num) = iirank
                   endif
                endif
             endif
          End do
       End Do
    End Do
    num_neighbor_ranks = num

  end subroutine set_neighbor_ranks

  subroutine set_rank_xyz( myproc, npe_x, npe_y, npe_z, rank_x, rank_y, rank_z )
    implicit none
    !
    integer, intent(in) :: myproc
    integer, intent(in) :: npe_x, npe_y, npe_z
    integer, intent(out) :: rank_x, rank_y, rank_z
    !
    rank_z = myproc /( npe_x *npe_y )
    rank_x = mod( mod( myproc, npe_x *npe_y ), npe_x )
    rank_y = mod( myproc, npe_x *npe_y ) /npe_x
  end subroutine set_rank_xyz

  subroutine set_grids_decompostion( &
       npe_x, npe_y, npe_z, rank_x, rank_y, rank_z, &
       ncomp_grids, gnum_grids, ndiv_grids, lnum_grids, &
       ista_grids, margin, mdist )
    use PARCEL_integer
    implicit none
    !
    integer npe_x, npe_y, npe_z, rank_x, rank_y, rank_z
    integer ncomp_grids
    integer ndiv_grids( ncomp_grids )
    integer gnum_grids( ncomp_grids )
    integer lnum_grids( ncomp_grids )
    integer ista_grids( ncomp_grids )
    integer margin
    integer(kind=intType) mdist
    !
    integer i, n1, n2, mx, my, mz
    !  
    ndiv_grids(1) = npe_x
    n1 = gnum_grids(1) /npe_x;  n2 = gnum_grids(1) -n1 *npe_x;
    Do i=0, npe_x -1
       if ( rank_x /= npe_x -1 ) then
          lnum_grids(1) = n1
       else
          lnum_grids(1) = n1 +n2
       endif
    ENd Do
    ista_grids(1) = n1 *rank_x +1

    if ( ncomp_grids >= 2 ) then
       ndiv_grids(2) = npe_y
       n1 = gnum_grids(2) /npe_y;  n2 = gnum_grids(2) -n1 *npe_y;
       Do i=0, npe_y -1
          if ( rank_y /= npe_y -1 ) then
             lnum_grids(2) = n1
          else
             lnum_grids(2) = n1 +n2
          endif
       ENd Do
       ista_grids(2) = n1 *rank_y +1
    endif

    if ( ncomp_grids >= 3 ) then
       ndiv_grids(3) = npe_z
       n1 = gnum_grids(3) /npe_z;  n2 = gnum_grids(3) -n1 *npe_z;
       Do i=0, npe_z -1
          if ( rank_z /= npe_z -1 ) then
             lnum_grids(3) = n1
          else
             lnum_grids(3) = n1 +n2
          endif
       ENd Do
       ista_grids(3) = n1 *rank_z +1
    endif

    mx = lnum_grids(1) +2*margin
    mdist = mx
    if ( ncomp_grids >= 2 ) then
       my = lnum_grids(2) +2*margin
       mdist = mdist *my
    endif
    if ( ncomp_grids >= 3 ) then
       mz = lnum_grids(3) +2*margin
       mdist = mdist *mz
    endif
  end subroutine set_grids_decompostion
  
  subroutine make_poisson_3D_DDM( &
       gnx,gny,gnz, &
       ndim, n_grid,  &
       istart_grid, &
       margin, nnd, m, &
       offset_dim, offset, &
       matA, vecB, vecX, vecX_exact,&
       rank_x, rank_y, rank_z, &
       npe_x, npe_y, npe_z )
    ! **********************************************************
    !                                  !!   outputs are in the DIA
    use PARCEL_integer
    implicit none
    !--------------
    ! inout
    !--------------
    integer(kind=intType), intent(in) :: gnx,gny,gnz
    integer, intent(in) :: ndim, margin
    integer(kind=intType), intent(in) :: nnd, m
    integer, intent(in) :: rank_x, rank_y, rank_z, npe_x, npe_y, npe_z
    integer, intent(in) :: n_grid( ndim )
    integer, intent(in) :: istart_grid( ndim )
    integer offset_dim( nnd, ndim )
    integer offset( nnd )
    real*8 matA( nnd*m )
    real*8 vecX(m), vecX_exact(m), vecB(m)

    !--------------
    ! local
    !--------------
    integer nx, ny, nz, mx, my, mz, bw
    integer jx, jy, jz, ix, iy, iz, j, num
    real*8 hx,hy,hz, z1, ctmp1, ctmp2, ctmp3


    !--------------
    ! init
    !--------------
    if ( ndim /= 3 ) stop "ndim is not 3"
    bw = margin

    nx = n_grid( 1 )   ! local
    ny = n_grid( 2 )
    nz = n_grid( 3 )

    mx = nx +bw *2;     my = ny +bw *2;        mz = nz +bw *2;
!    m = mx *my *mz

    matA = 0.0d0;  vecB = 0.0d0;  vecX = 0.0d0;  vecX_exact = 0.0d0
    offset_dim = 0;
    offset = 0
    !--------------
    ! begin
    !--------------
    Do jz=0, nz-1
       Do jy=0, ny-1
          Do jx=0, nx-1
             j = ( jx +bw ) +mx *( jy +bw )  +mx *my *( jz +bw ) +1;

             num = 0;

             if ( ( rank_z == 0 .and. jz > 0 ) .or. ( rank_z > 0 ) ) then
                matA( j +num*m ) = 1.0d0
             endif
             num = num +1;

             if ( ( rank_y == 0 .and. jy > 0 ) .or. ( rank_y > 0 ) ) then
                matA( j +num*m ) = 1.0d0
             endif
             num = num +1;

             if ( ( rank_x == 0 .and. jx > 0 ) .or. ( rank_x > 0 ) ) then
                matA( j +num*m ) = 1.0d0
             endif
             num = num +1;

             matA( j +num*m ) = -6.0d0        ! diagonal
             num = num +1;

             if ( ( rank_x == npe_x -1 .and. jx < nx-1 ) .or. ( rank_x < npe_x -1 ) ) then
                matA( j +num*m ) = 1.0d0
             endif
             num = num +1;

             if ( ( rank_y == npe_y -1 .and. jy < ny-1 ) .or. ( rank_y < npe_y -1 ) ) then
                matA( j +num*m ) = 1.0d0
             endif
             num = num +1;

             if ( ( rank_z == npe_z -1 .and. jz < nz-1 ) .or. ( rank_z < npe_z -1 ) ) then
                matA( j +num*m ) = 1.0d0
             endif
             num = num +1;

          End do
       End do
    End Do

    ! --------------
    ! set offset
    ! --------------
    offset_dim(1,1) = 0;    offset_dim(1,2) = 0;    offset_dim(1,3) = -1;  
    offset_dim(2,1) = 0;    offset_dim(2,2) =-1;    offset_dim(2,3) =  0;  
    offset_dim(3,1) =-1;    offset_dim(3,2) = 0;    offset_dim(3,3) =  0;  
    offset_dim(4,1) = 0;    offset_dim(4,2) = 0;    offset_dim(4,3) =  0;  
    offset_dim(5,1) = 1;    offset_dim(5,2) = 0;    offset_dim(5,3) =  0;  
    offset_dim(6,1) = 0;    offset_dim(6,2) = 1;    offset_dim(6,3) =  0;  
    offset_dim(7,1) = 0;    offset_dim(7,2) = 0;    offset_dim(7,3) =  1;  

    Do j=1, nnd
       offset(j) = &
            ( mx *my )*offset_dim(j,3) + mx *offset_dim(j,2) &
            + offset_dim(j,1)
    End Do

    ! --------------
    ! set (local) rhs vector
    ! --------------
    hx = 1.0d0/(gnx+1)
    hy = 1.0d0/(gny+1)
    hz = 1.0d0/(gnz+1)
!    z1 = 2.0d0 * h**6
    z1 = 2.0d0*(hx**2)*(hy**2)*(hz**2)
    Do jz=0, nz-1
       Do jy=0, ny-1
          Do jx=0, nx-1
             ix = jx +istart_grid(1) -1
             iy = jy +istart_grid(2) -1
             iz = jz +istart_grid(3) -1
             j = ( jx +bw ) +mx *( jy +bw )  +mx *my *( jz +bw ) +1


             ctmp1 = (ix+1)**2
             ctmp2 = (iy+1)**2
             ctmp3 = (iz+1)**2
             vecB(j) = z1 *( ctmp1*ctmp2 +ctmp2*ctmp3 +ctmp3*ctmp1 )
          End do
       End do
    End Do

    ! --------------
    ! boundary conditions
    ! --------------
    if ( rank_x == npe_x -1 ) then
       jx = nx -1
       Do jz=0, nz-1
          Do jy=0, ny-1
             ix = jx +istart_grid( 1 ) -1
             iy = jy +istart_grid( 2 ) -1
             iz = jz +istart_grid( 3 ) -1

             j = ( jx +bw ) +mx *( jy +bw )  +mx *my *( jz +bw ) +1

             vecB(j) = vecB(j) -( hy*(iy+1) )**2 *( hz*(iz+1) )**2
          End Do
       End Do
    endif

    if ( rank_y == npe_y -1 ) then
       jy = ny -1
       Do jz=0, nz -1
          Do jx=0, nx -1
             ix = jx +istart_grid( 1 ) -1
             iy = jy +istart_grid( 2 ) -1
             iz = jz +istart_grid( 3 ) -1

             j = ( jx +bw ) +mx *( jy +bw )  +mx *my *( jz +bw ) +1

             vecB(j) = vecB(j) -( hx*(ix+1) )**2 *( hz*(iz+1) )**2
          End Do
       End Do
    endif

    if ( rank_z == npe_z -1 ) then
       jz = nz -1
       Do jy=0, ny -1
          Do jx=0, nx -1
             ix = jx +istart_grid( 1 ) -1
             iy = jy +istart_grid( 2 ) -1
             iz = jz +istart_grid( 3 ) -1

             j = ( jx +bw ) +mx *( jy +bw )  +mx *my *( jz +bw ) +1

             vecB(j) = vecB(j) -( hx*(ix+1) )**2 *( hy*(iy+1) )**2
          End Do
       End Do
    endif

    ! --------------
    ! set exact solution
    ! --------------
!    z1 = h**6
    z1 = (hx**2)*(hy**2)*(hz**2) 
    Do jz=0, nz -1
       Do jy=0, ny -1
          Do jx=0, nx -1
             ix = jx +istart_grid( 1 ) -1
             iy = jy +istart_grid( 2 ) -1
             iz = jz +istart_grid( 3 ) -1
             j = ( jx +bw ) +mx *( jy +bw )  +mx *my *( jz +bw ) +1


             ctmp1 = (ix+1)**2
             ctmp2 = (iy+1)**2
             ctmp3 = (iz+1)**2
             vecx_exact(j) = z1 *( ctmp1*ctmp2*ctmp3 )

          End Do
       End Do
    End Do

    ! --------------
    ! set itital guess
    ! --------------
    vecx = 0.0d0

  end subroutine make_poisson_3D_DDM

end module ddm_setup
