module underlap
  use ddm_matv
  use ddmMat_wrap
  implicit none

#if 0
  integer,allocatable :: ioff_x(:), ioff_y(:), ioff_z(:)  
  integer,allocatable :: region_interior(:)
#endif

contains
  
  subroutine finish_underlap(underlap)
    
    implicit none
    type(under_lap) :: underlap 
    !  deallocate(Diag)
    deallocate(underlap%ioff_x)
    deallocate(underlap%ioff_y)
    deallocate(underlap%ioff_z)
    deallocate(underlap%region_interior)

  end subroutine finish_underlap

  subroutine init_underlap(a,prm,s_unl)
    use ddmMat_wrap
    
    implicit none

    type(ddmmatrix_wrap) a
    type(ddm_param) prm
    integer,intent(in) :: s_unl

    integer :: nx
    integer :: ny, nz
    integer :: mx, my, mz
    integer :: mxy
    integer :: stmx, stmy, stmz
    integer :: nxy, nxyz

    integer :: npe_x, npe_y, npe_z

    integer jz, jy,jx, j
    integer nxs, nxe, nys, nye, nzs, nze

    integer :: nnd
    integer i

    nx = prm%nx
    ny = prm%ny
    nz = prm%nz
    mx = prm%mx
    my = prm%my
    mz = prm%mz
    mxy= prm%mxy
    stmx = prm%stmx
    stmy = prm%stmy
    stmz = prm%stmz
    nxy = prm%nxy
    nxyz = prm%nxyz

    npe_x = prm%npe_x
    npe_y = prm%npe_y
    npe_z = prm%npe_z
    nnd = a%nnd_dia
    allocate( a%underlap%region_interior(a%nrow) )
    a%underlap%region_interior = 0
    nxs =0;   nxe = nx;    nys = 0;   nye = ny;    nzs = 0;   nze = nz
    if ( nx > 1 ) then
       nxs = s_unl
       nxe = nx -s_unl
    endif

    if ( ny > 1 ) then
       nys = s_unl
       nye = ny -s_unl
    endif

    if ( nz > 1 ) then
       nzs = s_unl
       nze = nz -s_unl
    endif

    !$omp parallel do collapse(2)
    do jz = nzs,nze-1
       do jy = nys,nye-1
          do jx = nxs,nxe-1
             j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz )+1
             a%underlap%region_interior(j) = 1
          enddo
       enddo
    enddo
    !$omp end parallel do


    allocate( a%underlap%ioff_x(nnd) )
    allocate( a%underlap%ioff_y(nnd) )
    allocate( a%underlap%ioff_z(nnd) )

    if ( prm%ncomp_grids == 1 ) then
       do i=1,nnd
          a%underlap%ioff_x(i) = a%ioff_grids(i,1)
          a%underlap%ioff_y(i) = 0
          a%underlap%ioff_z(i) = 0
       enddo
    endif
    if ( prm%ncomp_grids == 2 ) then
       do i=1,nnd
          a%underlap%ioff_x(i) = a%ioff_grids(i,1)
          a%underlap%ioff_y(i) = a%ioff_grids(i,2)
          a%underlap%ioff_z(i) = 0
       enddo
    endif
    if ( prm%ncomp_grids == 3 ) then
       do i=1,nnd
          a%underlap%ioff_x(i) = a%ioff_grids(i,1)
          a%underlap%ioff_y(i) = a%ioff_grids(i,2)
          a%underlap%ioff_z(i) = a%ioff_grids(i,3)
       enddo
    endif
  end subroutine init_underlap

  subroutine set_diag_underlap(a,prm,s_unl, Diag, iflag_DDM_independ)
    !$ use omp_lib
    use ddmMat_wrap
    use PARCEL_define,only : KSP_PRECONDITION_DDM_ILUindepend
    use PARCEL_integer    
    use common_kernel,only : get_maxthreads,set_precon_maxthreads    
    implicit none

    type(ddmmatrix_wrap),intent(in) :: a
    type(ddm_param) prm
    integer,intent(in) :: s_unl
    double precision,intent(inout) :: Diag(a%nrow)
    integer,intent(in) :: iflag_DDM_independ

    integer :: nx, ny, nz
    integer :: mx, my, mz
    integer :: mxy
    integer :: stmx, stmy, stmz
    integer :: nxy, nxyz
    integer :: div_direc_th

    integer :: npe_x, npe_y,  npe_z

    integer :: nnd
    integer, allocatable :: ioff_a(:)
    double precision, allocatable :: val(:)

    integer jz, jy,jx, j
    integer i, idiag, ip, k1, offset

    integer num_threads, my_thread
    integer nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th
    integer nxs, nxe, nys, nye, nzs, nze
    integer N
    integer(kind=intType) n_tmp
    integer max_numth_0

    nx = prm%nx
    ny = prm%ny
    nz = prm%nz
    mx = prm%mx
    my = prm%my
    mz = prm%mz
    mxy= prm%mxy
    stmx = prm%stmx
    stmy = prm%stmy
    stmz = prm%stmz
    nxy = prm%nxy
    nxyz = prm%nxyz
    div_direc_th = prm%div_direc_th

    npe_x = prm%npe_x
    npe_y = prm%npe_y
    npe_z = prm%npe_z
    nnd = a%nnd_dia
    N = a%nrow
    allocate(ioff_a(nnd))
    allocate(val(N*nnd))
    val = a%val_dia
    ioff_a = a%ioff_dia

    !$omp parallel do private(j) 
    do j=1,N
       Diag(j) = 1.0
    enddo
    !$omp end parallel do


    do i=1,nnd
       if ( ioff_a(i) == 0 ) then
          idiag = i
       endif
    enddo
    nxs = 0;   nxe = nx;    nys = 0;   nye = ny;    nzs = 0;   nze = nz
    if ( nx > 1 ) then
       nxs = s_unl
       nxe = nx -s_unl
    endif

    if ( ny > 1 ) then
       nys = s_unl
       nye = ny -s_unl
    endif

    if ( nz > 1 ) then
       nzs = s_unl
       nze = nz -s_unl
    endif

    !$omp parallel do collapse(2) private(jx,jy,jz,j) 
    do jz=nzs,nze-1
       do jy=nys,nye-1
          do jx=nxs,nxe-1
             j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
             Diag(j) = val( N *(idiag-1) +j )
          enddo
       enddo
    enddo
    !$omp end parallel do

    call get_maxthreads(max_numth_0) 
    if(div_direc_th.eq.1)then
       n_tmp = nx
    endif
    if(div_direc_th.eq.2)then
       n_tmp = ny
    endif
    if(div_direc_th.eq.3)then
       n_tmp = nz
    endif
    call set_precon_maxthreads(n_tmp,prm%precon_thblock)

    !$omp parallel private ( num_threads, my_thread ) &
    !$omp private( jz, jy, jx, j, offset, ip ) &
    !$omp private( nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th )

    num_threads = omp_get_num_threads()
    my_thread = omp_get_thread_num()

    call set_region_for_threads2(&
!    call set_region_for_threads(&
         my_thread, num_threads, &
         nxs, nxe, nys, nye, nzs, nze, &
         nxs_th, nxe_th, &
         nys_th, nye_th, &
         nzs_th, nze_th, &
         div_direc_th )
!    write(*,*)"init underlap --- ",my_thread,nzs_th,nze_th-1
    do jz=nzs_th,nze_th-1
       do jy=nys_th,nye_th-1
          do jx=nxs_th,nxe_th-1
             j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz )+1
             Diag(j) = 1.0 / Diag(j)

             do ip=1,nnd
                if ( ioff_a(ip) <= 0 ) cycle
                if ( jx +a%underlap%ioff_x(ip) < nxs_th .OR. jx +a%underlap%ioff_x(ip) > nxe_th -1 ) cycle
                if ( jy +a%underlap%ioff_y(ip) < nys_th .OR. jy +a%underlap%ioff_y(ip) > nye_th -1 ) cycle
                if ( jz +a%underlap%ioff_z(ip) < nzs_th .OR. jz +a%underlap%ioff_z(ip) > nze_th -1 ) cycle

                if(iflag_DDM_independ.eq.KSP_PRECONDITION_DDM_ILUindepend)then
                   if( ( a%underlap%ioff_y(ip).eq.0 ).and.( a%underlap%ioff_z(ip).eq.0 ) )then
                      cycle
                   endif
                endif

                offset = ioff_a(ip)
                if ( a%underlap%region_interior(j+offset) == 0 ) cycle

                Diag(j+offset) = Diag(j+offset) -val( N*(ip-1) +j )*val( N*(ip-1) +j ) *Diag(j)
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel

    ! do jz=nzs,nze-1
    !    do jy=nys,nye-1
    !       do jx=nxs,nxe-1
    !          j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
    !          write(*,*)"Diag",jx,jy,jz,Diag(j)
    !       enddo
    !    enddo
    ! enddo

    call omp_set_num_threads(max_numth_0)

  end subroutine set_diag_underlap

  ! アンダーラップ前処理
  ! ブロックヤコビ前処理で中心領域を計算後に周囲をポイントヤコビで計算する
  subroutine underlap_ddm( &
       n, nnd, underlap,&
       val_dia, ioff_dia, diag, vec1, vec2,&
       prm, ncomp_grids, num_grids, ioff_grids )

    use ddmMat_wrap
    use PARCEL_integer    
    use common_kernel,only : get_maxthreads,set_precon_maxthreads    
    implicit none

    integer,intent(in) :: n, nnd
    type(under_lap) :: underlap
    integer,intent(in) :: ioff_dia(nnd), ncomp_grids, num_grids(ncomp_grids), ioff_grids(nnd,ncomp_grids)
    double precision,intent(in) :: val_dia(n*nnd), diag(n)
    type(ddm_param),intent(in) :: prm
    double precision,intent(inout) :: vec1(n)
    double precision,intent(inout) :: vec2(n)
!    integer,intent(in) :: s_unl
!    integer,intent(in) :: s_ovp
    integer s_unl,s_ovp
    integer div_direc_th
    integer(kind=intType) n_tmp
    integer max_numth_0
    s_unl = prm%underlap_step
    s_ovp = 0
    div_direc_th = prm%div_direc_th

    call get_maxthreads(max_numth_0) 
    if(div_direc_th.eq.1)then
       n_tmp = prm%nx
    endif
    if(div_direc_th.eq.2)then
       n_tmp = prm%ny
    endif
    if(div_direc_th.eq.3)then
       n_tmp = prm%nz
    endif
    call set_precon_maxthreads(n_tmp,prm%precon_thblock)
    if ( div_direc_th == 0 ) then
       call ilu_ddm_case_nodiv_underlap( &
            nnd, n, underlap,&
            val_dia, ioff_dia, vec1, vec2,&
            Diag, ncomp_grids, prm, &
            s_unl, s_ovp, div_direc_th )
    else if ( div_direc_th == 1 ) then
       call ilu_ddm_case_xdiv_underlap( &
            nnd, n, underlap,&
            val_dia, ioff_dia, vec1, vec2,&
            Diag, ncomp_grids, prm, &
            s_unl, s_ovp, div_direc_th )
    else if ( div_direc_th == 2 ) then
       call ilu_ddm_case_ydiv_underlap( &
            nnd, n, underlap,&
            val_dia, ioff_dia, vec1, vec2,&
            Diag, ncomp_grids, prm, &
            s_unl, s_ovp, div_direc_th )
    else if ( div_direc_th == 3 ) then
       call ilu_ddm_case_zdiv_underlap( &
            nnd, n, underlap,&
            val_dia, ioff_dia, vec1, vec2,&
            Diag, ncomp_grids, prm, &
            s_unl, s_ovp, div_direc_th )
    endif
    call omp_set_num_threads(max_numth_0)

  end subroutine underlap_ddm

  subroutine ilu_ddm_case_nodiv_underlap( nnd, N, underlap,&
       a, ioff_a, vec1, vec2,&
       Diag, ncomp_grids, prm, &
       s_unl, s_ovp, div_direc_th )
    use ddmMat_wrap
    
    implicit none

    integer,intent(in) :: nnd
    integer,intent(in) :: N
    type(under_lap),intent(in) :: underlap
    double precision,intent(in) :: a(:)
    integer,intent(in) :: ioff_a(:)
    double precision,intent(in) :: vec1(:)
    double precision,intent(inout) :: vec2(:)
    double precision,intent(in) :: Diag(:)
    integer,intent(in) :: ncomp_grids
    type(ddm_param),intent(in) :: prm
    integer,intent(in) :: s_unl, s_ovp, div_direc_th

    integer :: nx, ny, nz
    integer :: mx, my, mz
    integer :: mxy
    integer :: stmx, stmy, stmz
    integer :: nxy, nxyz

    integer :: npe_x, npe_y,  npe_z

    integer jx, jy, jz, j, ip, k1, offset
    integer num_threads, my_thread, i
    integer nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th
    integer nxs, nxe, nys, nye, nzs, nze

    double precision sum
    integer ip_diag


    nx = prm%nx
    ny = prm%ny
    nz = prm%nz
    mx = prm%mx
    my = prm%my
    mz = prm%mz
    mxy= prm%mxy
    stmx = prm%stmx
    stmy = prm%stmy
    stmz = prm%stmz
    nxy = prm%nxy
    nxyz = prm%nxyz

    npe_x = prm%npe_x
    npe_y = prm%npe_y
    npe_z = prm%npe_z

    !$omp parallel do
    do j=1,N
       !vec2(j) = 0.0;
       vec2(j) = vec1(j)
    enddo
    !$omp end parallel do

    nxs = 0;   nxe = nx;    nys = 0;   nye = ny;    nzs = 0;   nze = nz
    if ( nx > 1 ) then
       nxs = s_unl;     nxe = nx -s_unl
    endif
    if ( ny > 1 ) then 
       nys = s_unl;     nye = ny -s_unl
    endif
    if ( nz > 1 ) then 
       nzs = s_unl;     nze = nz -s_unl
    endif

    nxs_th = nxs;   nxe_th = nxe
    nys_th = nys;   nye_th = nye
    nzs_th = nzs;   nze_th = nze


    !write(*,'"ppp ",2I10') nxs_th, nxe_th
    !write(*,'"ppp ",2I10') nys_th, nye_th
    !write(*,'"ppp ",2I10') nzs_th, nze_th


    do jz=nzs_th,nze_th-1
       do jy=nys_th,nye_th-1
          do jx=nxs_th,nxe_th-1
             j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
             sum = 0.0
             vec2(j) = 0.0

             do ip=1,nnd
                offset = ioff_a(ip)
                if ( offset > 0 ) cycle
                if ( underlap%region_interior( j+offset ) == 0 ) cycle

                sum = sum +a( N *(ip-1) +j ) *vec2(j+offset)
             enddo
             vec2(j) = Diag(j) *( vec1(j) -sum )
          enddo
       enddo
    enddo

    do jz=nze_th-1,nzs_th,-1
       do jy=nye_th-1,nys_th,-1
          do jx=nxe_th-1,nxs_th,-1
             j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz )+1

             sum = 0.0;
             do ip=1,nnd
                offset = ioff_a(ip)
                if ( offset <= 0 ) cycle
                if ( underlap%region_interior( j+offset ) == 0 ) cycle
                sum = sum +a( N *(ip-1) +j ) *vec2(j+offset)
             enddo
             vec2(j) = vec2(j) -Diag(j) *sum
          enddo
       enddo
    enddo


    do ip=1,nnd
       offset = ioff_a(ip)
       if ( offset == 0 ) then
          ip_diag = ip
       endif
    enddo

    nxs = 0;   nxe = nx;    nys = 0;   nye = ny;    nzs = 0;   nze = nz
    if ( nx > 1 ) then
       nxs = -s_ovp;     nxe = nx +s_ovp
    endif
    if ( ny > 1 ) then 
       nys = -s_ovp;     nye = ny +s_ovp
    endif
    if ( nz > 1 ) then   
       nzs = -s_ovp;     nze = nz +s_ovp
    endif

    nxs_th = nxs;   nxe_th = nxe
    nys_th = nys;   nye_th = nye
    nzs_th = nzs;   nze_th = nze

    do jz=nzs_th,nze_th-1
       do jy=nys_th,nye_th-1
          do jx=nxs_th,nxe_th-1
             j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz )+1
             if ( underlap%region_interior(j) == 1 ) cycle

             if ( abs( a(N*(ip_diag-1) +j) ) > 1.0e-6 ) then
                vec2(j) = vec1(j) / a( N*(ip_diag-1) +j)
             else
                vec2(j) = vec1(j) 
             endif
          enddo
       enddo
    enddo

  end subroutine ilu_ddm_case_nodiv_underlap

  subroutine ilu_ddm_case_xdiv_underlap( nnd, N, underlap,&
       a, ioff_a, vec1, vec2, &
       Diag, ncomp_grids, prm, &
       s_unl, s_ovp, div_direc_th )
    !$ use omp_lib
    use ddmMat_wrap
    
    implicit none

    integer,intent(in) ::N, nnd
    type(under_lap),intent(in) :: underlap
    double precision,intent(in) :: a(:)
    integer,intent(in) :: ioff_a(:)
    double precision,intent(in) :: vec1(:)
    double precision,intent(inout) :: vec2(:)
    integer,intent(in) :: ncomp_grids
    type(ddm_param),intent(in) :: prm
    integer,intent(in) :: s_unl, s_ovp, div_direc_th

    double precision,intent(in) :: Diag(:)
    integer :: nx, ny, nz
    integer :: mx, my, mz
    integer :: mxy
    integer :: stmx, stmy, stmz
    integer :: nxy, nxyz

    integer :: npe_x, npe_y,  npe_z

    integer jx, jy, jz, j, ip, k1, offset
    integer num_threads, my_thread, i
    integer nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th
    integer nxs, nxe, nys, nye, nzs, nze

    double precision sum
    integer ip_diag

    nx = prm%nx
    ny = prm%ny
    nz = prm%nz
    mx = prm%mx
    my = prm%my
    mz = prm%mz
    mxy= prm%mxy
    stmx = prm%stmx
    stmy = prm%stmy
    stmz = prm%stmz
    nxy = prm%nxy
    nxyz = prm%nxyz

    npe_x = prm%npe_x
    npe_y = prm%npe_y
    npe_z = prm%npe_z

    !$omp parallel do
    do j=1,N
       !vec2(j) = 0.0;
       vec2(j) = vec1(j)
    enddo
    !$omp end parallel do

    nxs = 0;   nxe = nx;    nys = 0;   nye = ny;    nzs = 0;   nze = nz
    if ( nx > 1 ) then
       nxs = s_unl;     nxe = nx -s_unl
    endif
    if ( ny > 1 ) then
       nys = s_unl;     nye = ny -s_unl
    endif
    if ( nz > 1 ) then
       nzs = s_unl;     nze = nz -s_unl
    endif

    !$omp parallel private ( num_threads, my_thread ) &
    !$omp private( i, jz, jy, jx, j, offset, ip )     &
    !$omp private( nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th ,sum)

    num_threads = omp_get_num_threads()
    my_thread = omp_get_thread_num()

    call set_region_for_threads2( my_thread, num_threads, nxs, nxe, nys, nye, nzs, nze, &
         nxs_th, nxe_th, nys_th, nye_th, &
         nzs_th, nze_th, div_direc_th )

    do jz=nzs_th,nze_th-1
       do jy=nys_th,nye_th-1
          do jx=nxs_th,nxe_th-1
             j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz )+1

             sum = 0.0
             vec2(j) = 0.0

             do ip=1,nnd
                offset = ioff_a(ip)
                if ( offset > 0 ) cycle
                if ( jx +underlap%ioff_x(ip) < nxs_th .OR. jx +underlap%ioff_x(ip) > nxe_th -1 ) cycle
                !if ( jy +ioff_y(ip) < nys_th .OR. jy +ioff_y(ip) > nye_th -1 ) cycle
                !if ( jz +ioff_z(ip) < nzs_th .OR. jz +ioff_z(ip) > nze_th -1 ) cycle
                if ( underlap%region_interior( j+offset ) == 0 ) cycle

                sum = sum +a( N *(ip-1) +j ) *vec2(j+offset)
             enddo
             vec2(j) = Diag(j) *( vec1(j) -sum )
          enddo
       enddo
    enddo
    !$omp end parallel

    !$omp parallel private ( num_threads, my_thread ) &
    !$omp private( i, jz, jy, jx, j, offset, ip )     &
    !$omp private( nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th , sum)

    num_threads = omp_get_num_threads()
    my_thread = omp_get_thread_num()

    call set_region_for_threads2( my_thread, num_threads, nxs, nxe, nys, nye, nzs, nze, &
         nxs_th, nxe_th, nys_th, nye_th, &
         nzs_th, nze_th, div_direc_th )

    do jz=nze_th-1,nzs_th,-1
       do jy=nye_th-1,nys_th,-1
          do jx=nxe_th-1,nxs_th,-1
             j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz )+1

             sum = 0.0
             do ip=1,nnd
                offset = ioff_a(ip)
                if ( offset <= 0 ) cycle
                if ( jx +underlap%ioff_x(ip) < nxs_th .OR. jx +underlap%ioff_x(ip) > nxe_th -1 ) cycle
                !if ( jy +ioff_y(ip) < nys_th .OR. jy +ioff_y(ip) > nye_th -1 ) cycle
                !if ( jz +ioff_z(ip) < nzs_th .OR. jz +ioff_z(ip) > nze_th -1 ) cycle
                if ( underlap%region_interior( j+offset ) == 0 ) cycle

                sum = sum +a( N *(ip-1) +j ) *vec2(j+offset)
             enddo
             vec2(j) = vec2(j) -Diag(j) *sum
          enddo
       enddo
    enddo
    !$omp end parallel

    do ip=1,nnd
       offset = ioff_a(ip)
       if ( offset == 0 ) then
          ip_diag = ip
       endif
    enddo

    nxs = 0;   nxe = nx;    nys = 0;   nye = ny;    nzs = 0;   nze = nz
    if ( nx > 1 ) then
       nxs = -s_ovp;     nxe = nx +s_ovp
    endif
    if ( ny > 1 ) then
       nys = -s_ovp;     nye = ny +s_ovp
    endif
    if ( nz > 1 ) then
       nzs = -s_ovp;     nze = nz +s_ovp
    endif

    !$omp parallel private ( num_threads, my_thread ) &
    !$omp private( i, jz, jy, jx, j, offset, ip )     &
    !$omp private( nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th )

    num_threads = omp_get_num_threads();
    my_thread = omp_get_thread_num();

    call set_region_for_threads2( my_thread, num_threads, nxs, nxe, nys, nye, nzs, nze, &
         nxs_th, nxe_th, nys_th, nye_th, &
         nzs_th, nze_th, div_direc_th )

    do jz=nzs_th,nze_th-1
       do jy=nys_th,nye_th-1
          do jx=nxs_th,nxe_th-1
             j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz )+1
             if ( underlap%region_interior(j) == 1 ) cycle

             if ( abs( a(N*(ip_diag-1) +j) ) > 1.0e-6 ) then
                vec2(j) = vec1(j) / a( N*(ip_diag-1) +j)
             else
                vec2(j) = vec1(j) 
             endif
          enddo
       enddo
    enddo
    !$omp end parallel

  end subroutine ilu_ddm_case_xdiv_underlap

  subroutine ilu_ddm_case_ydiv_underlap( nnd, N, underlap, &
       a, ioff_a, vec1, vec2, &
       Diag, ncomp_grids, prm, &
       s_unl, s_ovp, div_direc_th )
    !$ use omp_lib
    
    use ddmMat_wrap
    use PARCEL_define,only : KSP_PRECONDITION_DDM_ILUindepend
    implicit none

    integer,intent(in) :: N,nnd
    type(under_lap),intent(in) :: underlap
    double precision,intent(in) :: a(N*nnd)
    integer,intent(in) :: ioff_a(nnd)
    double precision,intent(in) :: vec1(N)
    double precision,intent(inout) :: vec2(N)
    integer,intent(in) :: ncomp_grids
    type(ddm_param),intent(in) :: prm
    integer,intent(in) :: s_unl, s_ovp, div_direc_th

    double precision,intent(in) :: Diag(n)
    integer :: nx, ny, nz
    integer :: mx, my, mz
    integer :: mxy
    integer :: stmx, stmy, stmz
    integer :: nxy, nxyz

    integer :: npe_x, npe_y,  npe_z

    integer jx, jy, jz, j, ip, k1, offset
    integer num_threads, my_thread, i
    integer nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th
    integer nxs, nxe, nys, nye, nzs, nze

    double precision sum
    integer ip_diag

    nx = prm%nx
    ny = prm%ny
    nz = prm%nz
    mx = prm%mx
    my = prm%my
    mz = prm%mz
    mxy= prm%mxy
    stmx = prm%stmx
    stmy = prm%stmy
    stmz = prm%stmz
    nxy = prm%nxy
    nxyz = prm%nxyz

    npe_x = prm%npe_x
    npe_y = prm%npe_y
    npe_z = prm%npe_z

    !$omp parallel do
    do j=1,N
       !vec2(j) = 0.0;
       vec2(j) = vec1(j)
    enddo
    !$omp end parallel do
    nxs = 0;   nxe = nx;    nys = 0;   nye = ny;    nzs = 0;   nze = nz
    if ( nx > 1 ) then
       nxs = s_unl;     nxe = nx -s_unl
    endif
    if ( ny > 1 ) then
       nys = s_unl;     nye = ny -s_unl
    endif
    if ( nz > 1 ) then
       nzs = s_unl;     nze = nz -s_unl
    endif

    !$omp parallel private ( num_threads, my_thread ) &
    !$omp private( i, jz, jy, jx, j, offset, ip )     &
    !$omp private( nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th ,sum )

    num_threads = omp_get_num_threads()
    my_thread = omp_get_thread_num()

    call set_region_for_threads2( my_thread, num_threads, nxs, nxe, nys, nye, nzs, nze, &
         nxs_th, nxe_th, nys_th, nye_th, &
         nzs_th, nze_th, div_direc_th )
    do jz=nzs_th,nze_th-1
       do jy=nys_th,nye_th-1
          do jx=nxs_th,nxe_th-1
             j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz )+1

             sum = 0.0;
             vec2(j) = 0.0;

             do ip=1,nnd
                offset = ioff_a(ip)
                if(prm%iflag_DDM_independ.eq.KSP_PRECONDITION_DDM_ILUindepend)then
                   if( (underlap%ioff_y(ip).eq.0).and.(underlap%ioff_z(ip).eq.0) )then
                      cycle
                   endif
                endif
 
                if ( offset > 0 ) cycle
                !if ( jx +ioff_x(ip) < nxs_th .OR. jx +ioff_x(ip) > nxe_th -1 ) cycle
                if ( jy +underlap%ioff_y(ip) < nys_th .OR. jy +underlap%ioff_y(ip) > nye_th -1 ) cycle
                !if ( jz +ioff_z(ip) < nzs_th .OR. jz +ioff_z(ip) > nze_th -1 ) cycle
                if ( underlap%region_interior( j+offset ) == 0 ) cycle

                sum = sum +a( N *(ip-1) +j ) *vec2(j+offset)
             enddo
             vec2(j) = Diag(j) *( vec1(j) -sum )
          enddo
       enddo
    enddo
    !$omp end parallel


    !$omp parallel private ( num_threads, my_thread ) &
    !$omp private( i, jz, jy, jx, j, offset, ip )     &
    !$omp private( nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th ,sum)

    num_threads = omp_get_num_threads();
    my_thread = omp_get_thread_num();

    call set_region_for_threads2( my_thread, num_threads, nxs, nxe, nys, nye, nzs, nze, &
         nxs_th, nxe_th, nys_th, nye_th, &
         nzs_th, nze_th, div_direc_th )

    do jz=nze_th-1,nzs_th,-1
       do jy=nye_th-1,nys_th,-1
          do jx=nxe_th-1,nxs_th,-1
             j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz )+1

             sum = 0.0
             do ip=1,nnd
                offset = ioff_a(ip)
                if(prm%iflag_DDM_independ.eq.KSP_PRECONDITION_DDM_ILUindepend)then
                   if( (underlap%ioff_y(ip).eq.0).and.(underlap%ioff_z(ip).eq.0) )then
                      cycle
                   endif
                endif

                if ( offset <= 0 ) cycle
                !if ( jx +ioff_x(ip) < nxs_th .OR. jx +ioff_x(ip) > nxe_th -1 ) cycle
                if ( jy +underlap%ioff_y(ip) < nys_th .OR. jy +underlap%ioff_y(ip) > nye_th -1 ) cycle
                !if ( jz +ioff_z(ip) < nzs_th .OR. jz +ioff_z(ip) > nze_th -1 ) cycle
                if ( underlap%region_interior( j+offset ) == 0 ) cycle

                sum = sum +a( N *(ip-1) +j ) *vec2(j+offset)
             enddo
             vec2(j) = vec2(j) -Diag(j) *sum
          enddo
       enddo
    enddo
    !$omp end parallel

    do ip=1,nnd
       offset = ioff_a(ip)
       if ( offset == 0 ) then
          ip_diag = ip
       endif
    enddo

    nxs = 0;   nxe = nx;    nys = 0;   nye = ny;    nzs = 0;   nze = nz
    if ( nx > 1 ) then
       nxs = -s_ovp;     nxe = nx +s_ovp
    endif
    if ( ny > 1 ) then
       nys = -s_ovp;     nye = ny +s_ovp
    endif
    if ( nz > 1 ) then
       nzs = -s_ovp;     nze = nz +s_ovp
    endif

    !$omp parallel private ( num_threads, my_thread ) &
    !$omp private( i, jz, jy, jx, j, offset, ip )     &
    !$omp private( nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th )

    num_threads = omp_get_num_threads()
    my_thread = omp_get_thread_num()

    call set_region_for_threads2( my_thread, num_threads, nxs, nxe, nys, nye, nzs, nze, &
         nxs_th, nxe_th, nys_th, nye_th, &
         nzs_th, nze_th, div_direc_th )
    do jz=nzs_th,nze_th-1
       do jy=nys_th,nye_th-1
          do jx=nxs_th,nxe_th-1
             j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz )+1
             if ( underlap%region_interior(j) == 1 ) cycle

             if ( abs( a(N*(ip_diag-1) +j) ) > 1.0e-6 ) then
                vec2(j) = vec1(j) / a( N*(ip_diag-1) +j)
             else
                vec2(j) = vec1(j) 
             endif
          enddo
       enddo
    enddo
    !$omp end parallel
  end subroutine ilu_ddm_case_ydiv_underlap

  subroutine ilu_ddm_case_zdiv_underlap( nnd, N, underlap,&
       a, ioff_a, vec1, vec2, &
       Diag, ncomp_grids, prm, &
       s_unl, s_ovp, div_direc_th )
    !$ use omp_lib

    use ddmMat_wrap
    use PARCEL_define,only : KSP_PRECONDITION_DDM_ILUindepend
    implicit none

    integer,intent(in) :: N,nnd
    type(under_lap),intent(in) :: underlap
    double precision,intent(in) :: a(:)
    integer,intent(in) :: ioff_a(:)
    double precision,intent(in) :: vec1(:)
    double precision,intent(inout) :: vec2(:)
    integer,intent(in) :: ncomp_grids
    type(ddm_param),intent(in) :: prm
    integer,intent(in) :: s_unl, s_ovp, div_direc_th

    double precision,intent(in) :: Diag(:)
    integer :: nx, ny, nz
    integer :: mx, my, mz
    integer :: mxy
    integer :: stmx, stmy, stmz
    integer :: nxy, nxyz

    integer :: npe_x, npe_y,  npe_z

    integer jx, jy, jz, j, ip, k1, offset
    integer num_threads, my_thread, i
    integer nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th
    integer nxs, nxe, nys, nye, nzs, nze

    double precision sum
    integer ip_diag
!    write(*,*)"underlap zdiv threads -- ",prm%iflag_DDM_independ,KSP_PRECONDITION_DDM_ILUindepend

    nx = prm%nx
    ny = prm%ny
    nz = prm%nz
    mx = prm%mx
    my = prm%my
    mz = prm%mz
    mxy= prm%mxy
    stmx = prm%stmx
    stmy = prm%stmy
    stmz = prm%stmz
    nxy = prm%nxy
    nxyz = prm%nxyz

    npe_x = prm%npe_x
    npe_y = prm%npe_y
    npe_z = prm%npe_z

    !$omp parallel do
    do j=1,N
       !vec2(j) = 0.0;
       vec2(j) = vec1(j)
    enddo
    !$omp end parallel do

    nxs = 0;   nxe = nx;    nys = 0;   nye = ny;    nzs = 0;   nze = nz
    if ( nx > 1 ) then
       nxs = s_unl;     nxe = nx -s_unl
    endif
    if ( ny > 1 ) then
       nys = s_unl;     nye = ny -s_unl
    endif
    if ( nz > 1 ) then
       nzs = s_unl;     nze = nz -s_unl
    endif

       !$omp parallel private ( num_threads, my_thread ) &
       !$omp private( i, jz, jy, jx, j, offset, ip )     &
       !$omp private( nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th ,sum)

       num_threads = omp_get_num_threads()
       my_thread = omp_get_thread_num()

       call set_region_for_threads2( my_thread, num_threads, nxs, nxe, nys, nye, nzs, nze, &
            nxs_th, nxe_th, nys_th, nye_th, &
            nzs_th, nze_th, div_direc_th )

       do jz=nzs_th,nze_th-1
          do jy=nys_th,nye_th-1
             do jx=nxs_th,nxe_th-1
                j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz )+1

                sum = 0.0
                vec2(j) = 0.0

                do ip=1,nnd
                   offset = ioff_a(ip)
                   if(prm%iflag_DDM_independ.eq.KSP_PRECONDITION_DDM_ILUindepend)then
                      if( (underlap%ioff_y(ip).eq.0).and.(underlap%ioff_z(ip).eq.0) )then
                         cycle
                      endif
                   endif
                   if ( offset > 0 ) cycle
                   !if ( jx +ioff_x(ip) < nxs_th .OR. jx +ioff_x(ip) > nxe_th -1 ) cycle
                   !if ( jy +ioff_y(ip) < nys_th .OR. jy +ioff_y(ip) > nye_th -1 ) cycle
                   if ( jz +underlap%ioff_z(ip) < nzs_th .OR. jz +underlap%ioff_z(ip) > nze_th -1 ) cycle
                   if ( underlap%region_interior( j+offset ) == 0 ) cycle

                   sum = sum +a( N *(ip-1) +j ) *vec2(j+offset)
                enddo
                vec2(j) = Diag(j) *( vec1(j) -sum );
             enddo
          enddo
       enddo
       !$omp end parallel

       !$omp parallel private ( num_threads, my_thread ) &
       !$omp private( i, jz, jy, jx, j, offset, ip )     &
       !$omp private( nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th , sum)

       num_threads = omp_get_num_threads();
       my_thread = omp_get_thread_num();

       call set_region_for_threads2( my_thread, num_threads, nxs, nxe, nys, nye, nzs, nze, &
            nxs_th, nxe_th, nys_th, nye_th, &
            nzs_th, nze_th, div_direc_th )

       do jz=nze_th-1,nzs_th,-1
          do jy=nye_th-1,nys_th,-1
             do jx=nxe_th-1,nxs_th,-1
                j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz )+1

                sum = 0.0;
                do ip=1,nnd
                   offset = ioff_a(ip)
                   if(prm%iflag_DDM_independ.eq.KSP_PRECONDITION_DDM_ILUindepend)then
                      if( (underlap%ioff_y(ip).eq.0).and.(underlap%ioff_z(ip).eq.0) )then
                         cycle
                      endif
                   endif
                   if ( offset <= 0 ) cycle
                   !if ( jx +ioff_x(ip) < nxs_th .OR. jx +ioff_x(ip) > nxe_th -1 ) cycle
                   !if ( jy +ioff_y(ip) < nys_th .OR. jy +ioff_y(ip) > nye_th -1 ) cycle
                   if ( jz +underlap%ioff_z(ip) < nzs_th .OR. jz +underlap%ioff_z(ip) > nze_th -1 ) cycle
                   if ( underlap%region_interior( j+offset ) == 0 ) cycle

                   sum = sum +a( N *(ip-1) +j ) *vec2(j+offset)
                enddo
                vec2(j) = vec2(j) -Diag(j) *sum
             enddo
          enddo
       enddo
       !$omp end parallel

    do ip=1,nnd
       offset = ioff_a(ip)
       if ( offset == 0 ) then
          ip_diag = ip
       endif
    enddo

    nxs = 0;   nxe = nx;    nys = 0;   nye = ny;    nzs = 0;   nze = nz
    if ( nx > 1 ) then
       nxs = -s_ovp;     nxe = nx +s_ovp
    endif
    if ( ny > 1 ) then
       nys = -s_ovp;     nye = ny +s_ovp
    endif
    if ( nz > 1 ) then
       nzs = -s_ovp;     nze = nz +s_ovp
    endif

    !$omp parallel private ( num_threads, my_thread ) &
    !$omp private( i, jz, jy, jx, j, offset, ip )     &
    !$omp private( nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th )

    num_threads = omp_get_num_threads()
    my_thread = omp_get_thread_num()

    call set_region_for_threads2( my_thread, num_threads, nxs, nxe, nys, nye, nzs, nze, &
         nxs_th, nxe_th, nys_th, nye_th, &
         nzs_th, nze_th, div_direc_th )

    do jz=nzs_th,nze_th-1
       do jy=nys_th,nye_th-1
          do jx=nxs_th,nxe_th-1
             j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz )+1
             if ( underlap%region_interior(j) == 1 ) cycle

             if ( abs( a(N*(ip_diag-1) +j) ) > 1.0e-6 ) then
                vec2(j) = vec1(j) / a( N*(ip_diag-1) +j)
             else
                vec2(j) = vec1(j) 
             endif
          enddo
       enddo
    enddo
    !$omp end parallel

  end subroutine ilu_ddm_case_zdiv_underlap

end module underlap
