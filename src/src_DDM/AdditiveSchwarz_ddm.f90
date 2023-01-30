
#include "ddfunc.h"

module param_AS_ddm
  use mpi
  use PARCEL_define,only : KSP_DOUBLE_DOUBLE_PRECISION,KSP_DOUBLE_PRECISION
contains
  subroutine get_xyz_neighbor_ranks(&
       num_neighbor_ranks, &
       neighbor_ranks, &
       npe_x,npe_y,npe_z, &
       neighbor_ranks_xyz)
    implicit none
    integer num_neighbor_ranks 
    integer neighbor_ranks(num_neighbor_ranks)
    integer npe_x,npe_y,npe_z
    integer neighbor_ranks_xyz(6)
    integer xpos(6), ypos(6), zpos(6)
    integer i
    integer k,kx,ky,kz
    integer kxmax,kymax,kzmax
    !
    ! 通信相手の番号 rank_* の取得。 相手がいない場合は -1
    !
    xpos(1) =-1;  ypos(1) = 0;  zpos(1) = 0  ! -x
    xpos(2) = 1;  ypos(2) = 0;  zpos(2) = 0  ! +x
    xpos(3) = 0;  ypos(3) =-1;  zpos(3) = 0  ! -y
    xpos(4) = 0;  ypos(4) = 1;  zpos(4) = 0  ! +y
    xpos(5) = 0;  ypos(5) = 0;  zpos(5) =-1  ! -z
    xpos(6) = 0;  ypos(6) = 0;  zpos(6) = 1  ! +z
    neighbor_ranks_xyz(:) = -1

    do i=1, 6
       k = 1;
       kxmax = 0;  kymax = 0;  kzmax = 0;
       if ( npe_z > 1 )  kzmax = 1;
       if ( npe_y > 1 )  kymax = 1; 
       if ( npe_x > 1 )  kxmax = 1; 
       Do kz=-kzmax, kzmax
          Do ky=-kymax, kymax
             Do kx=-kxmax, kxmax
                if ( kx*kx +ky*ky +kz*kz == 0 ) cycle
                if (neighbor_ranks(k) > -1) then
                   if ( kx==xpos(i) .and. ky==ypos(i) .and. kz==zpos(i) ) then
                      neighbor_ranks_xyz(i) = neighbor_ranks(k)
                   end if
                end if
                k = k +1
             End Do
          End Do
       End Do
    end do
  end subroutine get_xyz_neighbor_ranks

  subroutine set_AS_METHOD_TYPE_DDM(precon,astype)
    use param_AS_type_ddm,only : precon_AS_wrap
    implicit none
    integer astype
    type(precon_AS_wrap) precon
    precon%AS_METHOD_TYPE=astype
  end subroutine set_AS_METHOD_TYPE_DDM

  subroutine get_AS_region_DDM(&
       nx,ny,nz,&
       div_direc_th,&
       nrk_th_direc_p,nrk_th_direc_m, &
       addL,&
       precon,precon_thblock,&
       PRECISION_CALC, &
       PRECISION_PRECONDITION)

    use param_AS_type_ddm,only : precon_AS_wrap, &
         KSP_DDM_DIVIDE_THREADS_DIRECT_X, &
         KSP_DDM_DIVIDE_THREADS_DIRECT_Y, &
         KSP_DDM_DIVIDE_THREADS_DIRECT_Z
    use omp_lib
    use common_kernel,only : divide_region,get_maxthreads,set_precon_maxthreads
    use PARCEL_integer
    implicit none
    integer nx,ny,nz
    integer div_direc_th
    integer precon_thblock
    integer nrk_th_direc_p,nrk_th_direc_m
    integer addL
    integer nb
    integer mx,my,mz,m
    integer i,j
    integer max_threads
    type(precon_AS_wrap) precon

    integer(kind=intType) n
    integer my_thread
    integer num_threads
    integer max_numth_0
    integer(kind=intType) nxs_th,nxe_th
    integer(kind=intType) nys_th,nye_th
    integer(kind=intType) nzs_th,nze_th
    integer(kind=intType) nx0s,ny0s,nz0s
    integer(kind=intType) nx0e,ny0e,nz0e
    
    integer,optional :: PRECISION_CALC
    integer,optional :: PRECISION_PRECONDITION
    integer PRECISION

    PRECISION = KSP_DOUBLE_PRECISION
    if(  present(PRECISION_CALC)  )then
       PRECISION = PRECISION_CALC
    endif
  
    nx0s=1
    ny0s=1
    nz0s=1
    nx0e=nx
    ny0e=ny
    nz0e=nz
    if(nrk_th_direc_m > -1)then
       select case(div_direc_th)
       case(KSP_DDM_DIVIDE_THREADS_DIRECT_X)
          nx0s=nx0s+addL
       case(KSP_DDM_DIVIDE_THREADS_DIRECT_Y)
          ny0s=ny0s+addL
       case(KSP_DDM_DIVIDE_THREADS_DIRECT_Z)
          nz0s=nz0s+addL
       end select
    endif
    if(nrk_th_direc_p > -1)then
       select case(div_direc_th)
       case(KSP_DDM_DIVIDE_THREADS_DIRECT_X)
          nx0e=nx0e-addL
       case(KSP_DDM_DIVIDE_THREADS_DIRECT_Y)
          ny0e=ny0e-addL
       case(KSP_DDM_DIVIDE_THREADS_DIRECT_Z)
          nz0e=nz0e-addL
       end select
    endif

    call get_maxthreads(max_numth_0) 
    n = nx * ny * nz 
    call set_precon_maxthreads(n,precon_thblock)
!    write(*,*)"max_numth_0,precon_thblock",max_numth_0,precon_thblock
    call get_maxthreads(max_threads) 
    precon%max_threads = max_threads

    !    allocate(precon%AS_DDM(max_threads)
    allocate(precon%AS_DDM(max_threads))

    !$OMP Parallel private( num_threads, my_thread ) &
    !$OMP private( nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th )
    num_threads = omp_get_num_threads()
    my_thread   = omp_get_thread_num()

    nxs_th=nx0s-1
    nxe_th=nx0e
    nys_th=ny0s-1
    nye_th=ny0e
    nzs_th=nz0s-1
    nze_th=nz0e
    select case(div_direc_th)
    case(KSP_DDM_DIVIDE_THREADS_DIRECT_X)
       call divide_region(&
            my_thread,num_threads,nx0s,nx0e,nxs_th,nxe_th)
       nxs_th=nxs_th -1
       nxe_th=nxe_th 
    case(KSP_DDM_DIVIDE_THREADS_DIRECT_Y)
       call divide_region(&
            my_thread,num_threads,ny0s,ny0e,nys_th,nye_th)
       nys_th=nys_th -1
       nye_th=nye_th 
    case(KSP_DDM_DIVIDE_THREADS_DIRECT_Z)
       call divide_region(&
            my_thread,num_threads,nz0s,nz0e,nzs_th,nze_th)
       nzs_th=nzs_th -1
       nze_th=nze_th 
!    case(KSP_DDM_DIVIDE_THREADS_DIRECT_NONE)
!
    end select

    precon%AS_DDM(my_thread+1)%nxs_local_th=nxs_th
    precon%AS_DDM(my_thread+1)%nxe_local_th=nxe_th
    precon%AS_DDM(my_thread+1)%nys_local_th=nys_th
    precon%AS_DDM(my_thread+1)%nye_local_th=nye_th
    precon%AS_DDM(my_thread+1)%nzs_local_th=nzs_th
    precon%AS_DDM(my_thread+1)%nze_local_th=nze_th

    select case(div_direc_th)
    case(KSP_DDM_DIVIDE_THREADS_DIRECT_X)
       nxs_th=nxs_th-addL
       nxe_th=nxe_th+addL
       if(nxs_th<0)then
          nxs_th=0
       endif
       if(nx<nxe_th)then
          nxe_th=nx
       endif
    case(KSP_DDM_DIVIDE_THREADS_DIRECT_Y)
       nys_th=nys_th-addL
       nye_th=nye_th+addL
       if(nys_th<0)then
          nys_th=0
       endif
       if(ny<nye_th)then
          nye_th=ny
       endif
    case(KSP_DDM_DIVIDE_THREADS_DIRECT_Z)
       nzs_th=nzs_th-addL
       nze_th=nze_th+addL
       if(nzs_th<0)then
          nzs_th=0
       endif
       if(nz<nze_th)then
          nze_th=nz
       endif

    end select

    precon%AS_DDM(my_thread+1)%nxs_th = nxs_th
    precon%AS_DDM(my_thread+1)%nxe_th = nxe_th
    precon%AS_DDM(my_thread+1)%nys_th = nys_th
    precon%AS_DDM(my_thread+1)%nye_th = nye_th
    precon%AS_DDM(my_thread+1)%nzs_th = nzs_th
    precon%AS_DDM(my_thread+1)%nze_th = nze_th
    precon%AS_DDM(my_thread+1)%nx_th  = nxe_th-nxs_th
    precon%AS_DDM(my_thread+1)%ny_th  = nye_th-nys_th
    precon%AS_DDM(my_thread+1)%nz_th  = nze_th-nzs_th

    ! 加法シュワルツは袖なしの配列を利用するため
    precon%AS_DDM(my_thread+1)%mx_th  = &
         precon%AS_DDM(my_thread+1)%nx_th
    precon%AS_DDM(my_thread+1)%my_th  = &
         precon%AS_DDM(my_thread+1)%ny_th
    precon%AS_DDM(my_thread+1)%mz_th  = &
         precon%AS_DDM(my_thread+1)%nz_th
    precon%AS_DDM(my_thread+1)%mxy_th = &
         precon%AS_DDM(my_thread+1)%mx_th * &
         precon%AS_DDM(my_thread+1)%my_th 
    precon%AS_DDM(my_thread+1)%m_th = &
         precon%AS_DDM(my_thread+1)%mx_th * &
         precon%AS_DDM(my_thread+1)%my_th * &
         precon%AS_DDM(my_thread+1)%mz_th 

    !$OMP END Parallel

    do i=1,max_threads
       mx = precon%AS_DDM(i)%mx_th
       my = precon%AS_DDM(i)%my_th
       mz = precon%AS_DDM(i)%mz_th
       m  = mx*my*mz

       allocate(precon%AS_DDM(i)%Dinv(m))
       allocate(precon%AS_DDM(i)%rvec(m))
       allocate(precon%AS_DDM(i)%svec(m))

       !$OMP PARALLEL DO PRIVATE(j)
       do j=1,m
          precon%AS_DDM(i)%Dinv(j)=0.0d0
          precon%AS_DDM(i)%rvec(j)=0.0d0
          precon%AS_DDM(i)%svec(j)=0.0d0
       enddo
       !$OMP END PARALLEL DO

       ! call system("sleep 1")
       ! write(*,*)"PRECISION,KSP_DOUBLE_DOUBLE_PRECISION : ",PRECISION,KSP_DOUBLE_DOUBLE_PRECISION
       ! call system("sleep 2")
       if(PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
          if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
             allocate(precon%AS_DDM(i)%Dinv_lo(m))
             !$OMP PARALLEL DO PRIVATE(j)
             do j=1,m
                precon%AS_DDM(i)%Dinv_lo(j)=0.0d0
             enddo
          endif

          allocate(precon%AS_DDM(i)%rvec_lo(m))
          allocate(precon%AS_DDM(i)%svec_lo(m))
          !$OMP PARALLEL DO PRIVATE(j)
          do j=1,m
             precon%AS_DDM(i)%rvec_lo(j)=0.0d0
             precon%AS_DDM(i)%svec_lo(j)=0.0d0
          enddo
          !$OMP END PARALLEL DO
       endif
    enddo

    call omp_set_num_threads(max_numth_0)

  end subroutine get_AS_region_DDM

  subroutine free_AS_region_DDM(precon)
    use param_AS_type_ddm,only : precon_AS_wrap
    implicit none
    integer i
    integer max_threads
    type(precon_AS_wrap) precon

    max_threads=precon%max_threads
    do i=1,max_threads
       deallocate(precon%AS_DDM(i)%Dinv)
       deallocate(precon%AS_DDM(i)%rvec)
       deallocate(precon%AS_DDM(i)%svec)
       if(allocated(precon%AS_DDM(i)%Dinv_lo))then
          deallocate(precon%AS_DDM(i)%Dinv_lo)
       endif
       if(allocated(precon%AS_DDM(i)%rvec_lo))then
          deallocate(precon%AS_DDM(i)%rvec_lo)
       endif
       if(allocated(precon%AS_DDM(i)%svec_lo))then
          deallocate(precon%AS_DDM(i)%svec_lo)
       endif

    enddo
    deallocate(precon%AS_DDM)
    
  end subroutine free_AS_region_DDM

  ! DILU分解 (最大3次元までしか対応しない)
  subroutine init_AS_DILU_DDM(&
       m,&
       nnd,&
       ddmOffset,&
       ddmA,&
       ioff_stencil, &
       prm,&
       ddmDim, &
       precon, &
       PRECISION_CALC , &
       PRECISION_A, &
       PRECISION_PRECONDITION, &
       ddmA_lo &
       )
    use PARCEL_define
    use PARCEL_integer
    use param_AS_type_ddm,only : precon_AS_wrap
    use common_kernel,only : get_maxthreads,set_precon_maxthreads
    use ddmMat_wrap
    use omp_lib
!    use mpi
    implicit none
    !   include "ddm.h"

    integer(kind=intType) m
    integer(kind=intType) n
    integer nnd
    integer ddmDim
    integer ddmOffset(nnd)
    real*8 ddmA(m*nnd)
    integer ioff_stencil(nnd,ddmDim)
    type(ddm_param) :: prm

    integer(kind=intType) i
    integer(kind=intType) jz,jy,jx
    integer(kind=intType) ip,ipp
    integer nx,ny,nz
    integer mx,my,mxy
    integer stmx,stmy,stmz

    integer(kind=intType) j_th,j
    integer num_threads,my_thread
    integer max_numth_0
    integer nxs_th,nxe_th
    integer nys_th,nye_th
    integer nzs_th,nze_th
    integer mx_th,my_th,mz_th
    integer mxy_th,m_th

    integer ioff_x(nnd), ioff_y(nnd), ioff_z(nnd)
    integer diag_offset
    integer offset
    integer offset_m
    type(precon_AS_wrap) precon
    integer,optional :: PRECISION_CALC
    integer,optional :: PRECISION_A
    integer,optional :: PRECISION_PRECONDITION
    real*8,optional  :: ddmA_lo(m*nnd)
    integer myrank,ierr

    integer PRECISION
    real*8 dddiv_var1
    real*8 ddmul_var1
    real*8 ddmul_var2
    real*8 ddsub_var1
    real*8 c_hi,c_lo

!    write(*,*)"AS m = ",m

    PRECISION = KSP_DOUBLE_PRECISION
    if(  present(PRECISION_CALC)  )then
       PRECISION = PRECISION_CALC
    endif

    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr )

    do i=1,nnd
       ioff_x(i)=0
       ioff_y(i)=0
       ioff_z(i)=0
    enddo

    select case(ddmDim)
    case(1)
       do i=1, nnd
          ioff_x(i) = ioff_stencil(i,1)
       end do
    case(2)
       do i=1, nnd
          ioff_x(i) = ioff_stencil(i,1)
          ioff_y(i) = ioff_stencil(i,2)
       end do
    case(3)
       do i=1, nnd
          ioff_x(i) = ioff_stencil(i,1)
          ioff_y(i) = ioff_stencil(i,2)
          ioff_z(i) = ioff_stencil(i,3)
       end do
    end select

    nx   = prm%nx
    ny   = prm%ny
    nz   = prm%nz
    mx   = prm%mx
    my   = prm%my
    mxy  = prm%mxy
    stmx = prm%stmx
    stmy = prm%stmy
    stmz = prm%stmz
    diag_offset=-1
    do i=1,nnd
       if(ddmoffset(i).eq.0)then
          diag_offset=i
          exit
       endif
    enddo

    call get_maxthreads(max_numth_0) 
    n = nx * ny * nz 
    call set_precon_maxthreads(n,prm%precon_thblock)
! 加法シュバルツは袖なしのデータ構造なので
!    mx = nx, mxy=nx*ny : nx,ny,nzは拡張領域の長さ
!    stmx,stmy,stmz = 0
! である.
!    write(*,*)" AS :  mx,mxy ",mx,mxy
!    write(*,*)" AS : stmx,stmy,stmz : ",stmx,stmy,stmz
    !$OMP PARALLEL &
    !$OMP PRIVATE( &
    !$OMP num_threads,my_thread,&
    !$OMP nxs_th,nxe_th,nys_th,nye_th,nzs_th,nze_th, &
    !$OMP mx_th,my_th,mz_th,m_th,mxy_th, &
    !$OMP offset,offset_m, &
    !$OMP j_th,&
    !$OMP jz,jy,jx,j,ip,ipp, &
    !$OMP c_hi,c_lo, &
    !$OMP dddiv_var1, &
    !$OMP ddmul_var1, &
    !$OMP ddmul_var2, &
    !$OMP ddsub_var1  &
    !$OMP ) 

    num_threads = omp_get_num_threads()
    my_thread   = omp_get_thread_num()

    nxs_th = precon%AS_DDM(my_thread+1)%nxs_th
    nxe_th = precon%AS_DDM(my_thread+1)%nxe_th
    nys_th = precon%AS_DDM(my_thread+1)%nys_th
    nye_th = precon%AS_DDM(my_thread+1)%nye_th
    nzs_th = precon%AS_DDM(my_thread+1)%nzs_th
    nze_th = precon%AS_DDM(my_thread+1)%nze_th
    mx_th  = precon%AS_DDM(my_thread+1)%mx_th
    my_th  = precon%AS_DDM(my_thread+1)%my_th
    mz_th  = precon%AS_DDM(my_thread+1)%mz_th
    mxy_th = mx_th*my_th
    m_th   = mx_th * my_th * mz_th 

    ! write(*,*) "nxs_th, nxe_th = ", nxs_th, nxe_th
    ! write(*,*) "nys_th, nye_th = ", nys_th, nye_th
    ! write(*,*) "nzs_th, nze_th = ",  nzs_th, nze_th
    ! write(*,*) "mx_th , my_th , mz_th = ", mx_th , my_th , mz_th 
    ! write(*,*) "mxy_th, m_th = ", mxy_th, m_th  
    do jz=1,m_th
       precon%AS_DDM(my_thread+1)%Dinv(jz) = 1.0d0
    enddo

    if(PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then 
       if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
          do jz=1,m_th
             precon%AS_DDM(my_thread+1)%Dinv_lo(jz) = 0.0d0
          enddo
       endif
    endif

    do jz=nzs_th,nze_th-1
       do jy=nys_th,nye_th-1
          do jx=nxs_th,nxe_th-1
             j     = &
                  ( jx +stmx ) +  &
                  mx  * ( jy +stmy ) +  &
                  mxy * ( jz +stmz ) + 1
             j_th  = &
                  ( jx-nxs_th +stmx ) +  &
                  mx_th *( jy-nys_th +stmy ) +  &
                  mxy_th*( jz-nzs_th +stmz ) + 1
             precon%AS_DDM(my_thread+1)%Dinv(j_th) = ddmA( m*(diag_offset-1) +j )                     
             if(PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then 
                if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
!                   if(  present(ddmA_lo)  )then
                   if(PRECISION_A.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                      precon%AS_DDM(my_thread+1)%Dinv_lo(j_th) = ddmA_lo( m*(diag_offset-1) +j )
                   else
                      precon%AS_DDM(my_thread+1)%Dinv_lo(j_th) = 0.0d0
                   endif
                endif
             endif
          enddo
       enddo
    enddo

    do jz=nzs_th,nze_th-1
       do jy=nys_th,nye_th-1
          do jx=nxs_th,nxe_th-1
             j     = &
                  ( jx +stmx ) +  &
                  mx  * ( jy +stmy ) +  &
                  mxy * ( jz +stmz ) + 1
             j_th  = &
                  ( jx-nxs_th +stmx ) +  &
                  mx_th *( jy-nys_th +stmy ) +  &
                  mxy_th*( jz-nzs_th +stmz ) + 1
             if(PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then                
                   dddiv_calc1(1.0d0,0.0d0,precon%AS_DDM(my_thread+1)%Dinv(j_th),precon%AS_DDM(my_thread+1)%Dinv_lo(j_th),precon%AS_DDM(my_thread+1)%Dinv(j_th),precon%AS_DDM(my_thread+1)%Dinv_lo(j_th))
                else
                   dddiv_calc1(1.0d0,0.0d0,precon%AS_DDM(my_thread+1)%Dinv(j_th),0.0d0,c_hi,c_lo)
                   precon%AS_DDM(my_thread+1)%Dinv(j_th) = c_hi
                endif
             else
                precon%AS_DDM(my_thread+1)%Dinv(j_th) = 1.0d0/precon%AS_DDM(my_thread+1)%Dinv(j_th)                  
             endif
             do ip=1, nnd
                if ( ddmoffset(ip) <= 0 ) cycle
                offset = &
                     ioff_x(ip)         + &
                     ioff_y(ip) * mx_th + &
                     ioff_z(ip) * mxy_th 
                if ( jx +ioff_x(ip) < nxs_th .or. jx +ioff_x(ip) > nxe_th -1 ) cycle
                if ( jy +ioff_y(ip) < nys_th .or. jy +ioff_y(ip) > nye_th -1 ) cycle
                if ( jz +ioff_z(ip) < nzs_th .or. jz +ioff_z(ip) > nze_th -1 ) cycle
                if ( j +offset > m ) cycle


                if(prm%iflag_DDM_independ.eq.KSP_PRECONDITION_DDM_ILUindepend)then
                   if( ( ioff_y(ip).eq.0 ).and.( ioff_z(ip).eq.0 ) )then
                      cycle
                   endif
                endif

                do ipp=1,nnd
                   offset_m = &
                        ioff_x(ipp)         + &
                        ioff_y(ipp) * mx_th + &
                        ioff_z(ipp) * mxy_th 
                   offset_m =-offset_m 
                   if( offset.eq.offset_m )then
                      if(PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
!                         if( present(ddmA_lo)  )then
                         if(PRECISION_A.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                            ddmul_calc1(ddmA(m*(ip-1)+j),ddmA_lo(m*(ip-1)+j),ddmA(m*(ipp-1)+j+offset),ddmA_lo(m*(ipp-1)+j+offset),c_hi,c_lo)
                         else
                            ddmul_calc1(ddmA(m*(ip-1)+j),0.0d0,ddmA(m*(ipp-1)+j+offset),0.0d0,c_hi,c_lo)
                         endif
                         
                         if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then                
                            ddmul_calc1(c_hi,c_lo,precon%AS_DDM(my_thread+1)%Dinv(j_th),precon%AS_DDM(my_thread+1)%Dinv_lo(j_th),c_hi,c_lo)
                            ddsub_calc1(precon%AS_DDM(my_thread+1)%Dinv(j_th+offset),precon%AS_DDM(my_thread+1)%Dinv_lo(j_th+offset),c_hi,c_lo,precon%AS_DDM(my_thread+1)%Dinv(j_th+offset),precon%AS_DDM(my_thread+1)%Dinv_lo(j_th+offset))
                         else
                            ddmul_calc1(c_hi,c_lo,precon%AS_DDM(my_thread+1)%Dinv(j_th),0.0d0,c_hi,c_lo)
                            ddsub_calc1(precon%AS_DDM(my_thread+1)%Dinv(j_th+offset),0.0d0,c_hi,c_lo,precon%AS_DDM(my_thread+1)%Dinv(j_th+offset),c_lo)
                         endif
                      else
                         precon%AS_DDM(my_thread+1)%Dinv(j_th+offset) = &
                              precon%AS_DDM(my_thread+1)%Dinv(j_th+offset) - &
                              ddmA( m*(ip-1)  +j        ) * &
                              ddmA( m*(ipp-1) +j+offset ) * &
                              precon%AS_DDM(my_thread+1)%Dinv(j_th)
                      endif
                      exit
                   endif
                enddo
             enddo

          enddo
       enddo
    enddo
 !    if(myrank.eq.0)then
 !    do jz=1,m_th
 !       write(*,*)my_thread+1," Dinv - ",precon%AS_DDM(my_thread+1)%Dinv(jz)
 !    enddo
 !    do jz=nzs_th,nze_th-1
 !       do jy=nys_th,nye_th-1
 !          do jx=nxs_th,nxe_th-1
 !             j_th  = &
 !                  ( jx-nxs_th +stmx ) +  &
 !                  mx_th *( jy-nys_th +stmy ) +  &
 !                  mxy_th*( jz-nzs_th +stmz ) + 1
 !             write(*,*)my_thread+1,jx,jy,jz,&
 !                  precon%AS_DDM(my_thread+1)%Dinv(j_th)
 !          enddo
 !       enddo
 !    enddo
 ! endif
    !$OMP END PARALLEL 

    ! call MPI_Barrier(MPI_COMM_WORLD,ierr)
    ! stop

    call omp_set_num_threads(max_numth_0)

  end subroutine init_AS_DILU_DDM

  ! AS前処理

  subroutine solve_ILU_DDM(&
!       r, &
       m,nnd,&
       mx,mxy,&
       stmx,stmy,stmz,&
       m_th,&
       mx_th,mxy_th, &
       nxs_th,nxe_th,&
       nys_th,nye_th,&
       nzs_th,nze_th,&
       ioff_x,ioff_y,ioff_z, &
       ddmOffset,&
       ddmA, &
       Dinv,rvec,svec,&
       iflag_DDM_independ, &
       PRECISION_CALC, &
       PRECISION_PRECONDITION, &
       ddmA_lo, &
       Dinv_lo, &
       rvec_lo, &
       svec_lo  &
       )
    use PARCEL_define
    use PARCEL_integer
    implicit none

    integer(kind=intType) m
    integer nnd
!    real*8 r(m)
    integer mx,mxy
    integer stmx,stmy,stmz

    integer m_th
    integer mx_th,mxy_th
    integer nxs_th,nxe_th
    integer nys_th,nye_th
    integer nzs_th,nze_th

    integer ioff_x(nnd),ioff_y(nnd),ioff_z(nnd)
    integer ddmOffset(nnd)
    real*8 ddmA(m*nnd)

    real*8 Dinv(m_th)
    real*8 rvec(m_th)
    real*8 svec(m_th)

    integer,optional :: PRECISION_CALC
    integer,optional :: PRECISION_PRECONDITION
    real*8,optional :: ddmA_lo(m*nnd)
    real*8,optional :: Dinv_lo(m_th)
    real*8,optional :: rvec_lo(m_th)
    real*8,optional :: svec_lo(m_th)

    integer(kind=intType) j,j_th,jx,jy,jz,ip
    integer(kind=intType) offset
!    integer jx0(nnd)
    real*8 sum_lx0,sum_lx(nxs_th:nxe_th -1)
    real*8 sum_lx0_lo,sum_lx_lo(nxs_th:nxe_th -1)
    integer iflag_DDM_independ
    integer ipL,ipU,ipLn,ipUn
    integer  ip_diaU(nnd),ip_diaL(nnd)
    integer(kind=intType)  j0,j0_th
    real*8 ddmul_var1
    real*8 ddsub_var1
    real*8 ddadd_var1
    real*8 c_hi,c_lo
    integer PRECISION

    PRECISION = KSP_DOUBLE_PRECISION
    if(  present(PRECISION_CALC)  )then
       PRECISION = PRECISION_CALC
    endif
!    PRECISION = KSP_DOUBLE_PRECISION

    do ip=1, nnd
       ip_diaU(ip) = 0
       ip_diaL(ip) = 0
    enddo
    ipL = 1
    ipU = 1
    ipLn = 0
    ipUn = 0
    do ip=1, nnd
       offset = ddmOffset(ip)
       if(offset.eq.0) cycle
       if ( offset < 0 )then
          ip_diaL(ipL) = ip
          ipL  = ipL  + 1
          ipLn = ipLn + 1
       else
          ip_diaU(ipU) = ip
          ipU  = ipU  + 1
          ipUn = ipUn + 1
       endif
    enddo
    
    if(iflag_DDM_independ.eq.KSP_PRECONDITION_DDM_ILUindepend)then
       ! jx0 = -1
       ! do ipL=1,ipLn
       !    ip = ip_diaL(ipL)
       !    do jx=nxs_th, nxe_th -1
       !       if ( jx +ioff_x(ip) < nxs_th )then
       !          cycle
       !       endif
       !       jx0(ipL) = jx
       !       exit
       !    enddo
       ! enddo

       do jz=nzs_th, nze_th -1
          do jy=nys_th, nye_th -1
             
             do jx=nxs_th, nxe_th -1
                j_th  = &
                     ( jx-nxs_th +stmx ) +  &
                     mx_th *( jy-nys_th +stmy ) +  &
                     mxy_th*( jz-nzs_th +stmz ) + 1
                sum_lx(jx) = 0.0d0
                if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                   sum_lx_lo(jx) = 0.0d0
                endif
             enddo
             
             do ipL=1,ipLn
                ip = ip_diaL(ipL)
                offset = &
                     ioff_x(ip)         + &
                     ioff_y(ip) * mx_th + &
                     ioff_z(ip) * mxy_th 
                if( (ioff_y(ip).eq.0).and.(ioff_z(ip).eq.0) )then
                   cycle
                endif
                if ( jy +ioff_y(ip) < nys_th ) then
                   cycle
                endif
                if ( jz +ioff_z(ip) < nzs_th ) then
                   cycle
                endif
                if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                   if( PRECISION_PRECONDITION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                      do jx=nxe_th -1, nxs_th, -1
                         j =( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                         j_th  = &
                              ( jx-nxs_th +stmx ) +  &
                              mx_th *( jy-nys_th +stmy ) +  &
                              mxy_th*( jz-nzs_th +stmz ) + 1
                         ! sum_lx(jx) = sum_lx(jx) +ddmA( m*(ip-1) +j ) *svec(j_th+offset)
                         ddmul_calc1(ddmA(m*(ip-1)+j),ddmA_lo(m*(ip-1)+j),svec(j_th+offset),svec_lo(j_th+offset),c_hi,c_lo)
                         ddadd_calc1(c_hi,c_lo,sum_lx(jx),sum_lx_lo(jx),sum_lx(jx),sum_lx_lo(jx))
                      end do
                   else
                      do jx=nxe_th -1, nxs_th, -1
                         j =( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                         j_th  = &
                              ( jx-nxs_th +stmx ) +  &
                              mx_th *( jy-nys_th +stmy ) +  &
                              mxy_th*( jz-nzs_th +stmz ) + 1
                         ! sum_lx(jx) = sum_lx(jx) +ddmA( m*(ip-1) +j ) *svec(j_th+offset)
                         ddmul_calc1(ddmA(m*(ip-1)+j),0.0d0,svec(j_th+offset),svec_lo(j_th+offset),c_hi,c_lo)
                         ddadd_calc1(c_hi,c_lo,sum_lx(jx),sum_lx_lo(jx),sum_lx(jx),sum_lx_lo(jx))
                      end do
                   endif
                else
                   do jx=nxe_th -1, nxs_th, -1
                      j =( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                      j_th  = &
                           ( jx-nxs_th +stmx ) +  &
                           mx_th *( jy-nys_th +stmy ) +  &
                           mxy_th*( jz-nzs_th +stmz ) + 1
                      sum_lx(jx) = sum_lx(jx) +ddmA( m*(ip-1) +j ) *svec(j_th+offset)
                   end do
                endif
             end do
             if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                if( PRECISION_PRECONDITION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                   do jx=nxs_th, nxe_th -1
                      j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                      j_th  = &
                           ( jx-nxs_th +stmx ) +  &
                           mx_th *( jy-nys_th +stmy ) +  &
                           mxy_th*( jz-nzs_th +stmz ) + 1
                      ! svec(j_th) =Dinv(j_th) *( rvec(j_th) -sum_lx(jx) )
                      ddsub_calc1(rvec(j_th),rvec_lo(j_th) ,sum_lx(jx),sum_lx_lo(jx),c_hi,c_lo)
                      ddmul_calc1(Dinv(j_th),Dinv_lo(j_th),c_hi,c_lo,svec(j_th),svec_lo(j_th))
                   end do
                else
                   do jx=nxs_th, nxe_th -1
                      j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                      j_th  = &
                           ( jx-nxs_th +stmx ) +  &
                           mx_th *( jy-nys_th +stmy ) +  &
                           mxy_th*( jz-nzs_th +stmz ) + 1
                      ! svec(j_th) =Dinv(j_th) *( rvec(j_th) -sum_lx(jx) )
                      ddsub_calc1(rvec(j_th),rvec_lo(j_th) ,sum_lx(jx),sum_lx_lo(jx),c_hi,c_lo)
                      ddmul_calc1(Dinv(j_th),0.0d0,c_hi,c_lo,svec(j_th),svec_lo(j_th))
                   end do
                endif
             else
                do jx=nxs_th, nxe_th -1
                   j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                   j_th  = &
                        ( jx-nxs_th +stmx ) +  &
                        mx_th *( jy-nys_th +stmy ) +  &
                        mxy_th*( jz-nzs_th +stmz ) + 1
                   svec(j_th) =Dinv(j_th) *( rvec(j_th) -sum_lx(jx) )
                end do
             endif
          end do
       end do

       do jz=nze_th -1, nzs_th, -1
          do jy=nye_th -1, nys_th, -1

             do jx=nxe_th -1, nxs_th, -1
                sum_lx(jx) = 0.0d0
                if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                   sum_lx_lo(jx) = 0.0d0
                endif
             enddo
             
             j0 = &
                  ( -nxs_th +stmx ) +  &
                  mx_th *( jy-nys_th +stmy ) +  &
                  mxy_th*( jz-nzs_th +stmz ) + 1                     
             j0_th = ( stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
             do ipU=1,ipUn
                ip = ip_diaU(ipU)
                offset = &
                     ioff_x(ip)         + &
                     ioff_y(ip) * mx_th + &
                     ioff_z(ip) * mxy_th 
                if( (ioff_y(ip).eq.0).and.(ioff_z(ip).eq.0) )then
                   cycle
                endif
                if(  jy +ioff_y(ip) > nye_th -1 ) then
                   cycle
                endif
                if(  jz +ioff_z(ip) > nze_th -1 ) then
                   cycle
                endif

                if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                   if( PRECISION_PRECONDITION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                      do jx=nxe_th -1, nxs_th, -1
                         j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                         j_th  = &
                              ( jx-nxs_th +stmx ) +  &
                              mx_th *( jy-nys_th +stmy ) +  &
                              mxy_th*( jz-nzs_th +stmz ) + 1                   
                         ! sum_lx(jx) = sum_lx(jx) +ddmA( m*(ip-1) +j ) *svec(j_th+offset)
                         ddmul_calc1(ddmA( m*(ip-1) +j ),ddmA_lo( m*(ip-1) +j ) ,svec(j_th+offset),svec_lo(j_th+offset),c_hi,c_lo)
                         ddadd_calc1(c_hi,c_lo,sum_lx(jx),sum_lx_lo(jx),sum_lx(jx),sum_lx_lo(jx))
                      end do
                   else
                      do jx=nxe_th -1, nxs_th, -1
                         j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                         j_th  = &
                              ( jx-nxs_th +stmx ) +  &
                              mx_th *( jy-nys_th +stmy ) +  &
                              mxy_th*( jz-nzs_th +stmz ) + 1                   
                         ! sum_lx(jx) = sum_lx(jx) +ddmA( m*(ip-1) +j ) *svec(j_th+offset)
                         ddmul_calc1(ddmA( m*(ip-1) +j ),0.0d0 ,svec(j_th+offset),svec_lo(j_th+offset),c_hi,c_lo)
                         ddadd_calc1(c_hi,c_lo,sum_lx(jx),sum_lx_lo(jx),sum_lx(jx),sum_lx_lo(jx))

                      end do
                   endif
                else
                   do jx=nxe_th -1, nxs_th, -1
                      j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                      j_th  = &
                           ( jx-nxs_th +stmx ) +  &
                           mx_th *( jy-nys_th +stmy ) +  &
                           mxy_th*( jz-nzs_th +stmz ) + 1                   
                      sum_lx(jx) = sum_lx(jx) +ddmA( m*(ip-1) +j ) *svec(j_th+offset)
                   end do
                endif
             end do
             if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                if( PRECISION_PRECONDITION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                   do jx=nxe_th -1, nxs_th, -1
                      j_th  = &
                           ( jx-nxs_th +stmx ) +  &
                           mx_th *( jy-nys_th +stmy ) +  &
                           mxy_th*( jz-nzs_th +stmz ) + 1
                      ! svec(j_th) = svec(j_th) - Dinv(j_th)*sum_lx(jx)
                      ddmul_calc1(Dinv(j_th),Dinv_lo(j_th),sum_lx(jx),sum_lx_lo(jx),c_hi,c_lo)
                      ddsub_calc1(svec(j_th),svec_lo(j_th),c_hi,c_lo,svec(j_th),svec_lo(j_th))
                   end do                  
                else
                   do jx=nxe_th -1, nxs_th, -1
                      j_th  = &
                           ( jx-nxs_th +stmx ) +  &
                           mx_th *( jy-nys_th +stmy ) +  &
                           mxy_th*( jz-nzs_th +stmz ) + 1
                      ! svec(j_th) = svec(j_th) - Dinv(j_th)*sum_lx(jx)
                      ddmul_calc1(Dinv(j_th),0.0d0,sum_lx(jx),sum_lx_lo(jx),c_hi,c_lo)
                      ddsub_calc1(svec(j_th),svec_lo(j_th),c_hi,c_lo,svec(j_th),svec_lo(j_th))
                   end do
                endif
             else
                do jx=nxe_th -1, nxs_th, -1
                   j_th  = &
                        ( jx-nxs_th +stmx ) +  &
                        mx_th *( jy-nys_th +stmy ) +  &
                        mxy_th*( jz-nzs_th +stmz ) + 1
                   ! j_th = j0_th + jx
                   svec(j_th) = svec(j_th) - Dinv(j_th)*sum_lx(jx)
                end do
             endif
          end do
       end do
    else
       do jz=nzs_th, nze_th -1
          do jy=nys_th, nye_th -1
             do jx=nxs_th, nxe_th -1
                j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                j_th  = &
                     ( jx-nxs_th +stmx ) +  &
                     mx_th *( jy-nys_th +stmy ) +  &
                     mxy_th*( jz-nzs_th +stmz ) + 1

                sum_lx0 = 0.0d0
                svec(j_th) = 0.0d0
                if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                   sum_lx0_lo = 0.0d0
                   svec_lo(j_th) = 0.0d0
                endif
                if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then               
                   do ipL=1,ipLn
                      ip = ip_diaL(ipL)
                      if ( jx +ioff_x(ip) < nxs_th ) cycle
                      if ( jy +ioff_y(ip) < nys_th ) cycle
                      if ( jz +ioff_z(ip) < nzs_th ) cycle
                      offset = &
                           ioff_x(ip)         + &
                           ioff_y(ip) * mx_th + &
                           ioff_z(ip) * mxy_th 
                      if( PRECISION_PRECONDITION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                         ! sum_lx0 = sum_lx0 +ddmA( m*(ip-1) +j ) *svec(j_th+offset)
                         ddmul_calc1(ddmA(m*(ip-1)+j),ddmA_lo(m*(ip-1)+j),svec(j_th+offset),svec_lo(j_th+offset),c_hi,c_lo)
                         ddadd_calc1(sum_lx0,sum_lx0_lo,c_hi,c_lo,sum_lx0,sum_lx0_lo)
                      else
                         ! sum_lx0 = sum_lx0 +ddmA( m*(ip-1) +j ) *svec(j_th+offset)
                         ddmul_calc1(ddmA(m*(ip-1)+j),0.0d0,svec(j_th+offset),svec_lo(j_th+offset),c_hi,c_lo)
                         ddadd_calc1(sum_lx0,sum_lx0_lo,c_hi,c_lo,sum_lx0,sum_lx0_lo)
                      endif
                   end do
                   ! svec(j_th) =Dinv(j_th) *( rvec(j_th) -sum_lx0 )
                   ddsub_calc1(rvec(j_th),rvec_lo(j_th),sum_lx0,sum_lx0_lo,c_hi,c_lo)
                   if( PRECISION_PRECONDITION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                      ddmul_calc1(c_hi,c_lo,Dinv(j_th),Dinv_lo(j_th),svec(j_th),svec_lo(j_th))
                   else
                      ddmul_calc1(c_hi,c_lo,Dinv(j_th),0.0d0,svec(j_th),svec_lo(j_th))
                   endif
                else
                   do ipL=1,ipLn
                      ip = ip_diaL(ipL)
                      if ( jx +ioff_x(ip) < nxs_th ) cycle
                      if ( jy +ioff_y(ip) < nys_th ) cycle
                      if ( jz +ioff_z(ip) < nzs_th ) cycle
                      offset = &
                           ioff_x(ip)         + &
                           ioff_y(ip) * mx_th + &
                           ioff_z(ip) * mxy_th 
                      sum_lx0 = sum_lx0 +ddmA( m*(ip-1) +j ) *svec(j_th+offset)
                   end do
                   svec(j_th) =Dinv(j_th) *( rvec(j_th) -sum_lx0 )
                endif
                !                svec(j_th) =Dinv(j_th) *( r(j) -sum_lx0 )
             end do
          end do
       end do
       do jz=nze_th -1, nzs_th, -1
          do jy=nye_th -1, nys_th, -1
             do jx=nxe_th -1, nxs_th, -1
                j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                j_th  = &
                     ( jx-nxs_th +stmx ) +  &
                     mx_th *( jy-nys_th +stmy ) +  &
                     mxy_th*( jz-nzs_th +stmz ) + 1

                sum_lx0 = 0.0d0
                if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                   sum_lx0_lo = 0.0d0
                endif
                if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then               
                   do ipU=1,ipUn
                      ip = ip_diaU(ipU)
                      if (  jx +ioff_x(ip) > nxe_th -1 ) cycle
                      if (  jy +ioff_y(ip) > nye_th -1 ) cycle
                      if (  jz +ioff_z(ip) > nze_th -1 ) cycle
                      offset = &
                           ioff_x(ip)         + &
                           ioff_y(ip) * mx_th + &
                           ioff_z(ip) * mxy_th 
                      if( PRECISION_PRECONDITION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                         ! sum_lx0 = sum_lx0 +ddmA( m*(ip-1) +j ) *svec(j_th+offset)
                         ddmul_calc1(ddmA( m*(ip-1) +j ) ,ddmA_lo( m*(ip-1) +j ) ,svec(j_th+offset),svec_lo(j_th+offset),c_hi,c_lo)
                         ddadd_calc1(sum_lx0,sum_lx0_lo,c_hi,c_lo, sum_lx0,sum_lx0_lo)
                      else
                         ! sum_lx0 = sum_lx0 +ddmA( m*(ip-1) +j ) *svec(j_th+offset)
                         ddmul_calc1(ddmA( m*(ip-1) +j ) ,0.0d0 ,svec(j_th+offset),svec_lo(j_th+offset),c_hi,c_lo)
                         ddadd_calc1(sum_lx0,sum_lx0_lo,c_hi,c_lo, sum_lx0,sum_lx0_lo)
                      endif
                   end do
                   if( PRECISION_PRECONDITION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                      ! svec(j_th) = svec(j_th) - Dinv(j_th)*sum_lx0
                      ddmul_calc1(Dinv(j_th),Dinv_lo(j_th),sum_lx0,sum_lx0_lo,c_hi,c_lo)
                      ddsub_calc1(svec(j_th),svec_lo(j_th),c_hi,c_lo,svec(j_th),svec_lo(j_th))
                   else
                      ddmul_calc1(Dinv(j_th),0.0d0,sum_lx0,sum_lx0_lo,c_hi,c_lo)
                      ddsub_calc1(svec(j_th),svec_lo(j_th),c_hi,c_lo,svec(j_th),svec_lo(j_th))
                   endif
                else
                   do ipU=1,ipUn
                      ip = ip_diaU(ipU)
                      if (  jx +ioff_x(ip) > nxe_th -1 ) cycle
                      if (  jy +ioff_y(ip) > nye_th -1 ) cycle
                      if (  jz +ioff_z(ip) > nze_th -1 ) cycle
                      offset = &
                           ioff_x(ip)         + &
                           ioff_y(ip) * mx_th + &
                           ioff_z(ip) * mxy_th 
                      sum_lx0 = sum_lx0 +ddmA( m*(ip-1) +j ) *svec(j_th+offset)
                   end do
                   svec(j_th) = svec(j_th) - Dinv(j_th)*sum_lx0
                endif
             end do
          end do
       end do
    endif

  end subroutine solve_ILU_DDM

  subroutine AdditiveSchwarz_ILU_DDM(&
       m,nnd,&
!       r, &
       prm,precon,&
       ddmDim, &
       ddmOffset, &
       ddmA, &
       ioff_stencil, &
       PRECISION_CALC, &
       PRECISION_PRECONDITION, &
       ddmA_lo &
       )
    use omp_lib
    use ddmMat_wrap
    use common_kernel,only : get_maxthreads,set_precon_maxthreads
    use param_AS_type_ddm,only : precon_AS_wrap, &
         KSP_DDM_DIVIDE_THREADS_DIRECT_NONE, &
         KSP_DDM_DIVIDE_THREADS_DIRECT_X, &
         KSP_DDM_DIVIDE_THREADS_DIRECT_Y, &
         KSP_DDM_DIVIDE_THREADS_DIRECT_Z
    use PARCEL_integer

    implicit none
    !      include "ddm.h"
    integer(kind=intType) m
    integer nnd
    real*8 r(m)
    type(ddm_param) prm
    type(precon_AS_wrap) precon
    integer ddmDim
    integer ddmOffset(nnd)
    real*8 ddmA(m*nnd)
    integer ioff_stencil(nnd,ddmDim)

    integer(kind=intType) n_tmp
    integer max_numth_0
    integer mx,mxy
    integer stmx,stmy,stmz

    integer m_th
    integer mx_th,mxy_th
    integer nxs_th,nxe_th
    integer nys_th,nye_th
    integer nzs_th,nze_th
    integer myth
    integer ioff_x(nnd), ioff_y(nnd), ioff_z(nnd)
    integer(kind=intType) i
    integer div_direc_th
    integer PRECISION
    integer,optional :: PRECISION_CALC
    integer,optional :: PRECISION_PRECONDITION
    real*8,optional ::  ddmA_lo(m*nnd)

    PRECISION = KSP_DOUBLE_PRECISION
    if(  present(PRECISION_CALC)  )then
       PRECISION = PRECISION_CALC
    endif
  
    div_direc_th = prm%div_direc_th;

    do i=1,nnd
       ioff_x(i)=0
       ioff_y(i)=0
       ioff_z(i)=0
    enddo
    select case(ddmDim)
    case(1)
       do i=1, nnd
          ioff_x(i) = ioff_stencil(i,1)
       end do
    case(2)
       do i=1, nnd
          ioff_x(i) = ioff_stencil(i,1)
          ioff_y(i) = ioff_stencil(i,2)
       end do
    case(3)
       do i=1, nnd
          ioff_x(i) = ioff_stencil(i,1)
          ioff_y(i) = ioff_stencil(i,2)
          ioff_z(i) = ioff_stencil(i,3)
       end do
    end select

    mx   = prm%mx
    mxy  = prm%mxy
    stmx = prm%stmx
    stmy = prm%stmy
    stmz = prm%stmz

    call get_maxthreads(max_numth_0) 
    n_tmp = m
    call set_precon_maxthreads(n_tmp,prm%precon_thblock)

    !$OMP PARALLEL &
    !$OMP PRIVATE( &
    !$OMP myth,&
    !$OMP m_th,&
    !$OMP mx_th,mxy_th, &
    !$OMP nxs_th,nxe_th,&
    !$OMP nys_th,nye_th,&
    !$OMP nzs_th,nze_th &
    !$OMP )

    myth   = omp_get_thread_num()
    call set_AS_region_ddm(&
         precon, &
         myth,   &
         m_th,   &
         mx_th,  &
         mxy_th, &
         nxs_th,nxe_th, &
         nys_th,nye_th, &
         nzs_th,nze_th  &
         )
!    if ( div_direc_th == 0 ) then

    if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
       call solve_ILU_DDM(&
!            r, &
            m,nnd,&
            mx,mxy,&
            stmx,stmy,stmz,&
            m_th,&
            mx_th,mxy_th, &
            nxs_th,nxe_th,&
            nys_th,nye_th,&
            nzs_th,nze_th,&
            ioff_x,ioff_y,ioff_z, &
            ddmOffset,&
            ddmA, &
            precon%AS_DDM(myth+1)%Dinv,&
            precon%AS_DDM(myth+1)%rvec,&
            precon%AS_DDM(myth+1)%svec,&
            prm%iflag_DDM_independ, &
            PRECISION,&
            PRECISION_PRECONDITION, &
            ddmA_lo, &
            precon%AS_DDM(myth+1)%Dinv_lo,&
            precon%AS_DDM(myth+1)%rvec_lo,&
            precon%AS_DDM(myth+1)%svec_lo &
            )
    else
       call solve_ILU_DDM(&
!            r, &
            m,nnd,&
            mx,mxy,&
            stmx,stmy,stmz,&
            m_th,&
            mx_th,mxy_th, &
            nxs_th,nxe_th,&
            nys_th,nye_th,&
            nzs_th,nze_th,&
            ioff_x,ioff_y,ioff_z, &
            ddmOffset,&
            ddmA, &
            precon%AS_DDM(myth+1)%Dinv,&
            precon%AS_DDM(myth+1)%rvec,&
            precon%AS_DDM(myth+1)%svec,&
            prm%iflag_DDM_independ &
            )
    endif

    !$OMP END PARALLEL
    call omp_set_num_threads(max_numth_0)


  end subroutine AdditiveSchwarz_ILU_DDM

#if 0
  subroutine init_AS_Mat( &
       m,m_AS, &
       nnd, &
       prm,prm_AS,prm_g_AS,&
       ioff_stencil,&
       comm,A,A_AS,addL)
    use mpi
    use ddmMat_wrap,only : ddm_param,AS_geom

    implicit none
    integer m,m_AS
    integer nnd
    integer n
    type(ddm_param) :: prm
    type(ddm_param) :: prm_AS
    type(AS_geom) :: prm_g_AS ! 計算領域のオーバーラップに関する幾何情報
    integer ioff_stencil(nnd,prm%ncomp_grids)
    integer comm ! MPIコミュニケータ
    real*8 A(m*nnd)
    real*8 A_AS(m_AS*nnd)
    integer addL

    integer ip
    integer i,j,j_AS
    integer jx,jy,jz
    integer num_neighbor_ranks

    integer stmx, stmy, stmz
    integer nx, ny, nz
    integer mx, my, mz, mxy
    integer nxy,nxz,nyz

    integer stmx_AS,stmy_AS,stmz_AS
    integer nx_AS,ny_AS,nz_AS
    integer mx_AS,my_AS,mz_AS
    integer nxy_AS,nxz_AS,nyz_AS
    integer mxy_AS,mxz_AS,myz_AS

    integer addL_mx,addL_px
    integer addL_my,addL_py
    integer addL_mz,addL_pz

    integer ioff_x(nnd),ioff_y(nnd),ioff_z(nnd)
    integer ndata

    type comm_buff
       integer rank
       integer ndata
       real*8,allocatable :: buff(:)
    end type comm_buff
    type(comm_buff) comm_sbuff(prm%num_neighbor_ranks,2)
    type(comm_buff) comm_rbuff(prm%num_neighbor_ranks,2)


    type mpi_icomm_stat
       integer stat(MPI_STATUS_SIZE)
    end type mpi_icomm_stat
    type(mpi_icomm_stat) :: send_stat(prm%num_neighbor_ranks,2)
    type(mpi_icomm_stat) :: recv_stat(prm%num_neighbor_ranks,2)
    integer send_req(prm%num_neighbor_ranks,2)
    integer recv_req(prm%num_neighbor_ranks,2)

    integer ierr
    integer myrank
    integer mpi_req
    integer ncomp_grids
    integer np,nm
    integer neighbor_ranks_xyz(6)
    integer comm_ndata(3,3)
    call MPI_Comm_rank(comm,myrank,ierr)


    stmx = prm%stmx
    stmy = prm%stmy
    stmz = prm%stmz

    nx = prm%nx
    ny = prm%ny
    nz = prm%nz
    mx = prm%mx
    my = prm%my
    mz = prm%mz
    mxy = mx*my
    nxy = nx * ny
    nxz = nx * nz
    nyz = ny * nz
    ! ----------
    stmx_AS = prm_AS%stmx
    stmy_AS = prm_AS%stmy
    stmz_AS = prm_AS%stmz
    nx_AS = prm_AS%nx
    ny_AS = prm_AS%ny
    nz_AS = prm_AS%nz
    mx_AS = prm_AS%mx
    my_AS = prm_AS%my
    mz_AS = prm_AS%mz
    mxy_AS = mx_AS*my_AS
    nxy_AS = nx_AS * ny_AS
    nxz_AS = nx_AS * nz_AS
    nyz_AS = ny_AS * nz_AS
    mxy_AS = mx_AS * my_AS
    mxz_AS = mx_AS * mz_AS
    myz_AS = my_AS * mz_AS
    ! -----
    ! AS_g%sleeve_mx,AS_g%rank_s --> これらはここで作れる?
    ! sleeve_mxは拡張領域であれがaddLにするだけなのでできるのでは？
    ! rank_sは通信相手記憶するだけでは?
    ! 問題があるとすればAS前処理でベクトルの通信をする処理.自作すれば問題は解消する?

    ! 隣接ランク数(斜めを含む)
    num_neighbor_ranks = prm%num_neighbor_ranks

    ! 問題の次元数
    ncomp_grids = prm%ncomp_grids

    ! 3次元格子まで対応 
    call get_xyz_neighbor_ranks(&
         num_neighbor_ranks, &
         prm%neighbor_ranks, &
         prm%npe_x,prm%npe_y,prm%npe_z, &
         neighbor_ranks_xyz)

    ! ---------
    ! 送受信データの個数
    comm_ndata = 0
    ! 1次元格子
    comm_ndata(1,1) = nnd * addL
    ! 2次元格子
    comm_ndata(2,1) = my_AS * nnd * addL
    comm_ndata(2,2) = mx_AS * nnd * addL
    ! 3次元格子
    comm_ndata(3,1) = myz_AS * nnd * addL
    comm_ndata(3,2) = mxz_AS * nnd * addL
    comm_ndata(3,3) = mxy_AS * nnd * addL
    do i = 1, ncomp_grids
       ! 1 : -方向, 2 : +方向
       do j = 1,2 
          comm_sbuff(i,j)%rank = neighbor_ranks_xyz(j + (i-1)*2)
          comm_rbuff(i,j)%rank = neighbor_ranks_xyz(j + (i-1)*2)
          ndata = comm_ndata(ncomp_grids,i)
          comm_sbuff(i,j)%ndata = ndata
          comm_rbuff(i,j)%ndata = ndata
          allocate(comm_sbuff(i,j)%buff(ndata))
          allocate(comm_rbuff(i,j)%buff(ndata))
          comm_sbuff(i,j)%buff = 0.0d0
          comm_rbuff(i,j)%buff = 0.0d0
       enddo
    enddo

    ! ---------
    ! ローカル領域のコピー
    addL_mx = prm_g_AS%sleeve_mx 
    addL_px = prm_g_AS%sleeve_px 
    addL_my = prm_g_AS%sleeve_my 
    addL_py = prm_g_AS%sleeve_py 
    addL_mz = prm_g_AS%sleeve_mz 
    addL_pz = prm_g_AS%sleeve_pz 

    do jz = 0, nz-1
       do jy = 0, ny-1
          do jx = 0, nx-1
             j    = ( jx + stmx              ) + mx    * ( jy + stmy             ) + mxy    * ( jz + stmz             ) + 1
             j_AS = ( jx + stmx_AS + addL_mx ) + mx_AS * ( jy + stmy_AS + addL_my) + mxy_AS * ( jz + stmz_AS + addL_mz) + 1
             do ip = 1,nnd                
                A_AS(m_AS*(ip-1)+j_AS)=A(m*(ip-1)+j)             
             enddo
          enddo
       enddo
    enddo
    ! ---------

    ! -------
    ! irecvの発行
    do j = 1, 2 
       do i = 1, ncomp_grids 
          if(comm_rbuff(i,j)%rank>-1)then
             call MPI_Irecv( &
                  comm_rbuff(i,j)%buff , &
                  comm_rbuff(i,j)%ndata , &
                  MPI_REAL8, &
                  comm_rbuff(i,j)%rank , &
                  comm_rbuff(i,j)%rank , &
                  comm,&
                  recv_req(i,j), &
                  ierr)
          endif
       enddo
    enddo


    ! -----
    ! x方向のデータパック,アンパック
    if (ncomp_grids>0)then
       np = 0
       nm = 0
       do jz = 0, mz_AS-1
          do jy = 0, my_AS-1
             if(comm_sbuff(1,1)%rank>-1)then
                do jx = 1,addL
                   j_AS = ( jx-1 + stmx_AS + addL_mx  ) + mx_AS * ( jy ) + mxy_AS * ( jz ) + 1
                   do ip = 1,nnd                
                      np = np + 1
                      comm_sbuff(1,1)%buff(np) = A_AS(m_AS*(ip-1)+j_AS)
                   enddo
                enddo
             endif
             if(comm_sbuff(1,2)%rank>-1)then
                do jx = nx - addL + 1,nx
                   j_AS = ( jx-1 + stmx_AS + addL_mx  ) + mx_AS * ( jy ) + mxy_AS * ( jz ) + 1
                   do ip = 1,nnd                
                      nm = nm + 1
                      comm_sbuff(1,2)%buff(nm) = A_AS(m_AS*(ip-1)+j_AS)
                   enddo
                enddo
             endif
          enddo
       enddo
       i = 1
       do j = 1, 2 
          if(comm_sbuff(i,j)%rank>-1)then
             call MPI_Isend( &
                  comm_sbuff(i,j)%buff ,&
                  comm_sbuff(i,j)%ndata ,&
                  MPI_REAL8,&
                  comm_sbuff(i,j)%rank ,&
                  myrank,&
                  comm,&
                  send_req(i,j),&
                  ierr)
          endif
       enddo
       i = 1
       do j = 1, 2 
          if(comm_rbuff(i,j)%rank>-1)then
             call MPI_Wait(recv_req(i,j), recv_stat(i,j)%stat,ierr)
          endif
       enddo
       i = 1
       do j = 1, 2 
          if(comm_sbuff(i,j)%rank>-1)then
             call MPI_Wait(send_req(i,j),send_stat(i,j)%stat,ierr)
          endif
       enddo
       np = 0
       nm = 0
       do jz = 0, mz_AS-1
          do jy = 0, my_AS-1
             if(comm_rbuff(1,1)%rank>-1)then
                do jx = 1,addL
                   j_AS = ( jx-1 + stmx_AS ) + mx_AS * ( jy ) + mxy_AS * ( jz ) + 1
                   do ip = 1,nnd                
                      np = np + 1
                      A_AS(m_AS*(ip-1)+j_AS) = comm_rbuff(1,1)%buff(np) 
                   enddo
                enddo
             endif
             if(comm_rbuff(1,2)%rank>-1)then
                do jx = nx_AS - addL + 1,nx_AS
                   j_AS = ( jx-1 + stmx_AS ) + mx_AS * ( jy ) + mxy_AS * ( jz ) + 1
                   do ip = 1,nnd                
                      nm = nm + 1
                      A_AS(m_AS*(ip-1)+j_AS) = comm_rbuff(1,2)%buff(nm) 
                   enddo
                enddo
             endif
          enddo
       enddo

    endif

    ! y方向のデータパック,アンパック
    if (ncomp_grids>1)then
       np = 0
       nm = 0
       do jz = 0, mz_AS-1
          if(comm_sbuff(2,1)%rank>-1)then
             do jy = 1,addL
                do jx = 0, mx_AS-1
                   j_AS = ( jx ) + mx_AS * ( jy-1 + stmy_AS + addL_my) + mxy_AS * ( jz ) + 1
                   do ip = 1,nnd                
                      np = np + 1
                      comm_sbuff(2,1)%buff(np) = A_AS(m_AS*(ip-1)+j_AS)
                   enddo
                enddo
             enddo
          endif
          if(comm_sbuff(2,2)%rank>-1)then
             do jy = ny - addL + 1,ny
                do jx = 0, mx_AS-1
                   j_AS = ( jx ) + mx_AS * ( jy-1 + stmy_AS + addL_my) + mxy_AS * ( jz ) + 1
                   do ip = 1,nnd                
                      nm = nm + 1
                      comm_sbuff(2,2)%buff(nm) = A_AS(m_AS*(ip-1)+j_AS)
                   enddo
                enddo
             enddo
          endif
       enddo
       i = 2
       do j = 1, 2 
          if(comm_sbuff(i,j)%rank>-1)then
             call MPI_Isend( &
                  comm_sbuff(i,j)%buff ,&
                  comm_sbuff(i,j)%ndata ,&
                  MPI_REAL8,&
                  comm_sbuff(i,j)%rank ,&
                  myrank,&
                  comm,&
                  send_req(i,j),&
                  ierr)
          endif
       enddo
       i = 2
       do j = 1, 2 
          if(comm_rbuff(i,j)%rank>-1)then
             call MPI_Wait(recv_req(i,j), recv_stat(i,j)%stat,ierr)
          endif
       enddo
       i = 2
       do j = 1, 2 
          if(comm_sbuff(i,j)%rank>-1)then
             call MPI_Wait(send_req(i,j),send_stat(i,j)%stat,ierr)
          endif
       enddo
       np = 0
       nm = 0
       do jz = 0, mz_AS-1
          if(comm_rbuff(2,1)%rank>-1)then
             do jy = 1,addL
                do jx = 0, mx_AS-1
                   j_AS = ( jx ) + mx_AS * ( jy-1 + stmy_AS ) + mxy_AS * ( jz ) + 1
                   do ip = 1,nnd                
                      np = np + 1
                      A_AS(m_AS*(ip-1)+j_AS) = comm_rbuff(2,1)%buff(np)
                   enddo
                enddo
             enddo
          endif
          if(comm_rbuff(2,2)%rank>-1)then
             do jy = ny_AS - addL + 1,ny_AS
                do jx = 0, mx_AS-1
                   j_AS = ( jx ) + mx_AS * ( jy-1 + stmy_AS ) + mxy_AS * ( jz ) + 1
                   do ip = 1,nnd                
                      nm = nm + 1
                      A_AS(m_AS*(ip-1)+j_AS) = comm_rbuff(2,2)%buff(nm)
                   enddo
                enddo
             enddo
          endif
       enddo
    endif

    ! z方向のデータパック,アンパック
    if (ncomp_grids>2)then
       np = 0
       nm = 0
       if(comm_sbuff(3,1)%rank>-1)then
          do jz = 1,addL
             do jy = 0, my_AS-1
                do jx = 0, mx_AS-1
                   j_AS = ( jx ) + mx_AS * ( jy ) + mxy_AS * ( jz-1 + stmz_AS + addL_mz) + 1
                   do ip = 1,nnd                
                      np = np + 1
                      comm_sbuff(3,1)%buff(np) = A_AS(m_AS*(ip-1)+j_AS)
                   enddo
                enddo
             enddo
          enddo
       endif
       if(comm_sbuff(3,2)%rank>-1)then
          do jz = nz - addL + 1,nz
             do jy = 0, my_AS-1
                do jx = 0, mx_AS-1
                   j_AS = ( jx ) + mx_AS * ( jy ) + mxy_AS * ( jz-1 + stmz_AS + addL_mz) + 1
                   do ip = 1,nnd                
                      nm = nm + 1
                      comm_sbuff(3,2)%buff(nm) = A_AS(m_AS*(ip-1)+j_AS)
                   enddo
                enddo
             enddo
          enddo
       endif
       i = 3
       do j = 1, 2 
          if(comm_sbuff(i,j)%rank>-1)then
             call MPI_Isend( &
                  comm_sbuff(i,j)%buff ,&
                  comm_sbuff(i,j)%ndata ,&
                  MPI_REAL8,&
                  comm_sbuff(i,j)%rank ,&
                  myrank,&
                  comm,&
                  send_req(i,j),&
                  ierr)
          endif
       enddo
       i = 3
       do j = 1, 2 
          if(comm_rbuff(i,j)%rank>-1)then
             call MPI_Wait(recv_req(i,j), recv_stat(i,j)%stat,ierr)
          endif
       enddo
       i = 3
       do j = 1, 2 
          if(comm_sbuff(i,j)%rank>-1)then
             call MPI_Wait(send_req(i,j),send_stat(i,j)%stat,ierr)
          endif
       enddo

       np = 0
       nm = 0
       if(comm_rbuff(3,1)%rank>-1)then
          do jz = 1,addL
             do jy = 0, my_AS-1
                do jx = 0, mx_AS-1
                   j_AS = ( jx ) + mx_AS * ( jy ) + mxy_AS * ( jz-1 + stmz_AS ) + 1
                   do ip = 1,nnd                
                      np = np + 1
                      A_AS(m_AS*(ip-1)+j_AS) = comm_rbuff(3,1)%buff(np)
                   enddo
                enddo
             enddo
          enddo
       endif
       if(comm_rbuff(3,2)%rank>-1)then
          do jz = nz_AS - addL + 1,nz_AS
             do jy = 0, my_AS-1
                do jx = 0, mx_AS-1
                   j_AS = ( jx ) + mx_AS * ( jy ) + mxy_AS * ( jz-1 + stmz_AS ) + 1

                   do ip = 1,nnd                
                      nm = nm + 1
                       A_AS(m_AS*(ip-1)+j_AS) = comm_rbuff(3,2)%buff(nm) 
                   enddo
                enddo
             enddo
          enddo
       endif


    endif

    do i = 1, ncomp_grids
       do j = 1,2
          deallocate(comm_sbuff(i,j)%buff)
          deallocate(comm_rbuff(i,j)%buff)
       enddo
    enddo

  end subroutine init_AS_Mat
#endif

  subroutine set_AS_vec_DDM(&
       prm,precon,m,r,  &
       PRECISION_CALC,r_lo &
       )
    use omp_lib
    use ddmMat_wrap
    use PARCEL_define
    use param_AS_type_ddm,only : precon_AS_wrap
    use PARCEL_integer
    use common_kernel,only : get_maxthreads,set_precon_maxthreads
    implicit none
    !      include "ddm.h"
    type(ddm_param) :: prm
    type(precon_AS_wrap) precon
    integer(kind=intType) m
    real*8 r(m)
    real*8,optional :: r_lo(m)
    integer,optional :: PRECISION_CALC

    integer mx,mxy
    integer stmx,stmy,stmz

    integer myth
    integer m_th
    integer mx_th
    integer mxy_th
    integer nxs_th,nxe_th
    integer nys_th,nye_th
    integer nzs_th,nze_th
    integer(kind=intType) jz,jy,jx,j,j_th
    integer nxs_local,nxe_local
    integer nys_local,nye_local
    integer nzs_local,nze_local
    integer iflagAS 
    integer precon_thblock
    integer(kind=intType) n_tmp
    integer max_numth_0
    integer PRECISION

    PRECISION = KSP_DOUBLE_PRECISION
    if(  present(PRECISION_CALC)  )then
       PRECISION = PRECISION_CALC
    endif

    mx   = prm%mx
    mxy  = prm%mxy
    stmx = prm%stmx
    stmy = prm%stmy
    stmz = prm%stmz
    iflagAS=precon%AS_METHOD_TYPE
    precon_thblock=prm%precon_thblock

    call get_maxthreads(max_numth_0) 
    n_tmp =m
    call set_precon_maxthreads(n_tmp,precon_thblock)
    if( (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_BASIC)&
         .or. &
         (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_RESTRICT) )then
       !$OMP PARALLEL &
       !$OMP PRIVATE( &
       !$OMP jz,jy,jx,j,j_th,&
       !$OMP myth,&
       !$OMP m_th,&
       !$OMP mx_th,&
       !$OMP mxy_th, &
       !$OMP nxs_th,nxe_th,&
       !$OMP nys_th,nye_th,&
       !$OMP nzs_th,nze_th,&
       !$OMP nxs_local,nxe_local,&
       !$OMP nys_local,nye_local,&
       !$OMP nzs_local,nze_local &
       !$OMP )
       myth   = omp_get_thread_num()
       call set_AS_region_ddm(&
            precon, &
            myth,   &
            m_th,   &
            mx_th,  &
            mxy_th, &
            nxs_th,nxe_th, &
            nys_th,nye_th, &
            nzs_th,nze_th  &
            )

       do j_th=1,m_th
          precon%AS_DDM(myth+1)%rvec(j_th)=0.0d0
       enddo
       if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
          do j_th=1,m_th
             precon%AS_DDM(myth+1)%rvec_lo(j_th)=0.0d0
          enddo
       endif
       do jz=nzs_th,nze_th-1
          do jy=nys_th,nye_th-1
             do jx=nxs_th,nxe_th-1
                j     = &
                     ( jx +stmx ) +  &
                     mx  * ( jy +stmy ) +  &
                     mxy * ( jz +stmz ) + 1
                j_th  = &
                     ( jx-nxs_th +stmx ) +  &
                     mx_th *( jy-nys_th +stmy ) +  &
                     mxy_th*( jz-nzs_th +stmz ) + 1
                precon%AS_DDM(myth+1)%rvec(j_th)=r(j)
             enddo
             if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                do jx=nxs_th,nxe_th-1
                   j     = &
                        ( jx +stmx ) +  &
                        mx  * ( jy +stmy ) +  &
                        mxy * ( jz +stmz ) + 1
                   j_th  = &
                        ( jx-nxs_th +stmx ) +  &
                        mx_th *( jy-nys_th +stmy ) +  &
                        mxy_th*( jz-nzs_th +stmz ) + 1
                   precon%AS_DDM(myth+1)%rvec_lo(j_th)=r_lo(j)
                enddo
             endif
          enddo
       enddo
       !$OMP END PARALLEL 

    else
       !$OMP PARALLEL &
       !$OMP PRIVATE( &
       !$OMP jz,jy,jx,j,j_th,&
       !$OMP myth,&
       !$OMP m_th,&
       !$OMP mx_th,mxy_th, &
       !$OMP nxs_th,nxe_th,&
       !$OMP nys_th,nye_th,&
       !$OMP nzs_th,nze_th, &
       !$OMP nxs_local,nxe_local,&
       !$OMP nys_local,nye_local,&
       !$OMP nzs_local,nze_local &
       !$OMP ) 
       myth   = omp_get_thread_num()
       call set_AS_region_ddm(&
            precon, &
            myth,   &
            m_th,   &
            mx_th,  &
            mxy_th, &
            nxs_th,nxe_th, &
            nys_th,nye_th, &
            nzs_th,nze_th  &
            )
       do j_th=1,m_th
          precon%AS_DDM(myth+1)%rvec(j_th)=0.0d0
       enddo
       if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
          do j_th=1,m_th
             precon%AS_DDM(myth+1)%rvec_lo(j_th)=0.0d0
          enddo
       endif
       nxs_local=precon%AS_DDM(myth+1)%nxs_local_th
       nxe_local=precon%AS_DDM(myth+1)%nxe_local_th
       nys_local=precon%AS_DDM(myth+1)%nys_local_th
       nye_local=precon%AS_DDM(myth+1)%nye_local_th
       nzs_local=precon%AS_DDM(myth+1)%nzs_local_th
       nze_local=precon%AS_DDM(myth+1)%nze_local_th
       do jz=nzs_local,nze_local-1
          do jy=nys_local,nye_local-1
             do jx=nxs_local,nxe_local-1
                j     = &
                     ( jx +stmx ) +  &
                     mx  * ( jy +stmy ) +  &
                     mxy * ( jz +stmz ) + 1
                j_th  = &
                     ( jx-nxs_th +stmx ) +  &
                     mx_th *( jy-nys_th +stmy ) +  &
                     mxy_th*( jz-nzs_th +stmz ) + 1         
                precon%AS_DDM(myth+1)%rvec(j_th)=r(j) 
             enddo
             if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                do jx=nxs_local,nxe_local-1
                   j     = &
                        ( jx +stmx ) +  &
                        mx  * ( jy +stmy ) +  &
                        mxy * ( jz +stmz ) + 1
                   j_th  = &
                        ( jx-nxs_th +stmx ) +  &
                        mx_th *( jy-nys_th +stmy ) +  &
                        mxy_th*( jz-nzs_th +stmz ) + 1         
                   precon%AS_DDM(myth+1)%rvec_lo(j_th)=r_lo(j) 
                enddo                
             endif
          enddo
       enddo
       !$OMP END PARALLEL 
    endif

    call omp_set_num_threads(max_numth_0)

  end subroutine set_AS_vec_DDM

  subroutine get_AS_vec_DDM(&
       prm,precon,m,s,addL, &
       PRECISION_CALC,s_lo &
       )
    use omp_lib
    use ddmMat_wrap
    use PARCEL_define
!    use param_AS_type_ddm,only : precon_AS_wrap
    use PARCEL_integer
    use common_kernel,only : get_maxthreads,set_precon_maxthreads
    use param_AS_type_ddm
    implicit none
    !      include "ddm.h"
    type(ddm_param) :: prm
    type(precon_AS_wrap) precon
    integer(kind=intType) m
    real*8 s(m)
    integer,optional :: PRECISION_CALC
    real*8,optional ::  s_lo(m)

    integer addL
    integer mx,mxy
    integer stmx,stmy,stmz

    integer myth
    integer m_th
    integer mx_th
    integer mxy_th
    integer nxs_th,nxe_th
    integer nys_th,nye_th
    integer nzs_th,nze_th

    integer nxs_local,nxe_local
    integer nys_local,nye_local
    integer nzs_local,nze_local

    integer(kind=intType) jz,jy,jx,j,j_th
    integer iflagAS
    integer precon_thblock
    integer div_direc_th
    integer(kind=intType) n_tmp
    integer max_numth_0

    integer nxs_th0,nxe_th0
    integer nys_th0,nye_th0
    integer nzs_th0,nze_th0
    real*8 ddadd_var1

    integer PRECISION

    PRECISION = KSP_DOUBLE_PRECISION
    if(  present(PRECISION_CALC)  )then
       PRECISION = PRECISION_CALC
    endif
!    PRECISION = KSP_DOUBLE_PRECISION

    mx   = prm%mx
    mxy  = prm%mxy
    stmx = prm%stmx
    stmy = prm%stmy
    stmz = prm%stmz
    iflagAS=precon%AS_METHOD_TYPE
    precon_thblock=prm%precon_thblock
    div_direc_th=prm%div_direc_th

    call get_maxthreads(max_numth_0) 
    n_tmp =m
    call set_precon_maxthreads(n_tmp,precon_thblock)
    if( (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_BASIC)&
         .or. &
         (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_INTERPOLATE) )then
       !$OMP PARALLEL DO PRIVATE(j)
       do j=1,m
          s(j)=0.0d0
       end do
       if(PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
          !$OMP PARALLEL DO PRIVATE(j)
          do j=1,m
             s_lo(j)=0.0d0
          end do
       endif
#if 1
       !$OMP PARALLEL &
       !$OMP PRIVATE( &
       !$OMP jz,jy,jx,j,j_th,&
       !$OMP myth,&
       !$OMP m_th,&
       !$OMP mx_th,mxy_th, &
       !$OMP nxs_th,nxe_th,&
       !$OMP nys_th,nye_th,&
       !$OMP nzs_th,nze_th, &
       !$OMP ddadd_var1 &
       !$OMP ) 
       myth   = omp_get_thread_num()
       call set_AS_region_ddm(&
            precon, &
            myth,   &
            m_th,   &
            mx_th,  &
            mxy_th, &
            nxs_th,nxe_th, &
            nys_th,nye_th, &
            nzs_th,nze_th  &
            )
       nxs_th0 = nxs_th
       nxe_th0 = nxe_th
       nys_th0 = nys_th
       nye_th0 = nye_th
       nzs_th0 = nzs_th
       nze_th0 = nze_th
       select case(div_direc_th)
       case(KSP_DDM_DIVIDE_THREADS_DIRECT_X)
          do jz=nzs_th,nze_th-1
             do jy=nys_th,nye_th-1
                if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                   do jx=nxs_th,nxs_th+addL-1
                      j     = &
                           ( jx +stmx ) +  &
                           mx  * ( jy +stmy ) +  &
                           mxy * ( jz +stmz ) + 1
                      j_th  = &
                           ( jx-nxs_th +stmx ) +  &
                           mx_th *( jy-nys_th +stmy ) +  &
                           mxy_th*( jz-nzs_th +stmz ) + 1         
                      ! s(j)=s(j)+ precon%AS_DDM(myth+1)%svec(j_th)
                      ddadd_calc1(s(j),s_lo(j),precon%AS_DDM(myth+1)%svec(j_th),precon%AS_DDM(myth+1)%svec_lo(j_th),s(j),s_lo(j))
                   enddo
                else
                   do jx=nxs_th,nxs_th+addL-1
                      j     = &
                           ( jx +stmx ) +  &
                           mx  * ( jy +stmy ) +  &
                           mxy * ( jz +stmz ) + 1
                      j_th  = &
                           ( jx-nxs_th +stmx ) +  &
                           mx_th *( jy-nys_th +stmy ) +  &
                           mxy_th*( jz-nzs_th +stmz ) + 1         
                      s(j)=s(j)+ precon%AS_DDM(myth+1)%svec(j_th)
                   enddo
                endif
             enddo
          enddo
          !$OMP BARRIER          
          do jz=nzs_th,nze_th-1
             do jy=nys_th,nye_th-1
!                do jx=nxs_th,nxe_th-1
                if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                   do jx=nxs_th+addL,nxe_th-1-addL
                      j     = &
                           ( jx +stmx ) +  &
                           mx  * ( jy +stmy ) +  &
                           mxy * ( jz +stmz ) + 1
                      j_th  = &
                           ( jx-nxs_th +stmx ) +  &
                           mx_th *( jy-nys_th +stmy ) +  &
                           mxy_th*( jz-nzs_th +stmz ) + 1         
                      ! s(j)=s(j)+ precon%AS_DDM(myth+1)%svec(j_th)
                      ddadd_calc1(s(j),s_lo(j),precon%AS_DDM(myth+1)%svec(j_th),precon%AS_DDM(myth+1)%svec_lo(j_th),s(j),s_lo(j))
                   enddo
                else
                   do jx=nxs_th+addL,nxe_th-1-addL
                      j     = &
                           ( jx +stmx ) +  &
                           mx  * ( jy +stmy ) +  &
                           mxy * ( jz +stmz ) + 1
                      j_th  = &
                           ( jx-nxs_th +stmx ) +  &
                           mx_th *( jy-nys_th +stmy ) +  &
                           mxy_th*( jz-nzs_th +stmz ) + 1         
                      s(j)=s(j)+ precon%AS_DDM(myth+1)%svec(j_th)
                   enddo
                endif
             enddo
          enddo
          !$OMP BARRIER          
          do jz=nzs_th,nze_th-1
             do jy=nys_th,nye_th-1
!                do jx=nxs_th,nxe_th-1
                if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                   do jx=nxe_th-1-addL+1,nxe_th-1
                      j     = &
                           ( jx +stmx ) +  &
                           mx  * ( jy +stmy ) +  &
                           mxy * ( jz +stmz ) + 1
                      j_th  = &
                           ( jx-nxs_th +stmx ) +  &
                           mx_th *( jy-nys_th +stmy ) +  &
                           mxy_th*( jz-nzs_th +stmz ) + 1         
                      ! s(j)=s(j)+ precon%AS_DDM(myth+1)%svec(j_th)
                      ddadd_calc1(s(j),s_lo(j),precon%AS_DDM(myth+1)%svec(j_th),precon%AS_DDM(myth+1)%svec_lo(j_th),s(j),s_lo(j))
                   enddo
                else
                   do jx=nxe_th-1-addL+1,nxe_th-1
                      j     = &
                           ( jx +stmx ) +  &
                           mx  * ( jy +stmy ) +  &
                           mxy * ( jz +stmz ) + 1
                      j_th  = &
                           ( jx-nxs_th +stmx ) +  &
                           mx_th *( jy-nys_th +stmy ) +  &
                           mxy_th*( jz-nzs_th +stmz ) + 1         
                      s(j)=s(j)+ precon%AS_DDM(myth+1)%svec(j_th)
                   enddo
                endif
             enddo
          enddo
       case(KSP_DDM_DIVIDE_THREADS_DIRECT_Y)
          do jz=nzs_th,nze_th-1
!             do jy=nys_th,nye_th-1
             do jy=nys_th,nys_th+addL-1
                if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                   do jx=nxs_th,nxe_th-1
                      j     = &
                           ( jx +stmx ) +  &
                           mx  * ( jy +stmy ) +  &
                           mxy * ( jz +stmz ) + 1
                      j_th  = &
                        ( jx-nxs_th +stmx ) +  &
                        mx_th *( jy-nys_th +stmy ) +  &
                        mxy_th*( jz-nzs_th +stmz ) + 1         
                      ! s(j)=s(j)+ precon%AS_DDM(myth+1)%svec(j_th)
                      ddadd_calc1(s(j),s_lo(j),precon%AS_DDM(myth+1)%svec(j_th),precon%AS_DDM(myth+1)%svec_lo(j_th),s(j),s_lo(j))
                   enddo
                else
                   do jx=nxs_th,nxe_th-1
                      j     = &
                           ( jx +stmx ) +  &
                           mx  * ( jy +stmy ) +  &
                           mxy * ( jz +stmz ) + 1
                      j_th  = &
                        ( jx-nxs_th +stmx ) +  &
                        mx_th *( jy-nys_th +stmy ) +  &
                        mxy_th*( jz-nzs_th +stmz ) + 1         
                      s(j)=s(j)+ precon%AS_DDM(myth+1)%svec(j_th)
                   enddo
                endif
             enddo
          enddo
          !$OMP BARRIER          
          do jz=nzs_th,nze_th-1
!             do jy=nys_th,nye_th-1
             do jy=nys_th+addL,nye_th-1-addL
                if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                   do jx=nxs_th,nxe_th-1
                      j     = &
                           ( jx +stmx ) +  &
                           mx  * ( jy +stmy ) +  &
                           mxy * ( jz +stmz ) + 1
                      j_th  = &
                           ( jx-nxs_th +stmx ) +  &
                           mx_th *( jy-nys_th +stmy ) +  &
                           mxy_th*( jz-nzs_th +stmz ) + 1         
                     !  s(j)=s(j)+ precon%AS_DDM(myth+1)%svec(j_th)
                      ddadd_calc1(s(j),s_lo(j),precon%AS_DDM(myth+1)%svec(j_th),precon%AS_DDM(myth+1)%svec_lo(j_th),s(j),s_lo(j))
                   enddo
                else
                   do jx=nxs_th,nxe_th-1
                      j     = &
                           ( jx +stmx ) +  &
                           mx  * ( jy +stmy ) +  &
                           mxy * ( jz +stmz ) + 1
                      j_th  = &
                           ( jx-nxs_th +stmx ) +  &
                           mx_th *( jy-nys_th +stmy ) +  &
                           mxy_th*( jz-nzs_th +stmz ) + 1         
                      s(j)=s(j)+ precon%AS_DDM(myth+1)%svec(j_th)
                   enddo
                endif
             enddo
          enddo
          !$OMP BARRIER          
          do jz=nzs_th,nze_th-1
!             do jy=nys_th,nye_th-1
             do jy=nye_th-1-addL+1,nye_th-1
                if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                   do jx=nxs_th,nxe_th-1
                      j     = &
                           ( jx +stmx ) +  &
                           mx  * ( jy +stmy ) +  &
                           mxy * ( jz +stmz ) + 1
                      j_th  = &
                           ( jx-nxs_th +stmx ) +  &
                           mx_th *( jy-nys_th +stmy ) +  &
                           mxy_th*( jz-nzs_th +stmz ) + 1         
                      ! s(j)=s(j)+ precon%AS_DDM(myth+1)%svec(j_th)
                      ddadd_calc1(s(j),s_lo(j),precon%AS_DDM(myth+1)%svec(j_th),precon%AS_DDM(myth+1)%svec_lo(j_th),s(j),s_lo(j))
                   enddo
                else
                   do jx=nxs_th,nxe_th-1
                      j     = &
                           ( jx +stmx ) +  &
                           mx  * ( jy +stmy ) +  &
                           mxy * ( jz +stmz ) + 1
                      j_th  = &
                           ( jx-nxs_th +stmx ) +  &
                           mx_th *( jy-nys_th +stmy ) +  &
                           mxy_th*( jz-nzs_th +stmz ) + 1         
                      s(j)=s(j)+ precon%AS_DDM(myth+1)%svec(j_th)
                   enddo
                endif
             enddo
          enddo
       case(KSP_DDM_DIVIDE_THREADS_DIRECT_Z)
!          do jz=nzs_th,nze_th-1
          do jz=nzs_th,nzs_th+addL-1
             do jy=nys_th,nye_th-1
                if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                   do jx=nxs_th,nxe_th-1
                      j     = &
                           ( jx +stmx ) +  &
                           mx  * ( jy +stmy ) +  &
                           mxy * ( jz +stmz ) + 1
                      j_th  = &
                           ( jx-nxs_th +stmx ) +  &
                           mx_th *( jy-nys_th +stmy ) +  &
                           mxy_th*( jz-nzs_th +stmz ) + 1         
                      ! s(j)=s(j)+ precon%AS_DDM(myth+1)%svec(j_th)
                      ddadd_calc1(s(j),s_lo(j),precon%AS_DDM(myth+1)%svec(j_th),precon%AS_DDM(myth+1)%svec_lo(j_th),s(j),s_lo(j))
                   enddo
                else
                   do jx=nxs_th,nxe_th-1
                      j     = &
                           ( jx +stmx ) +  &
                           mx  * ( jy +stmy ) +  &
                           mxy * ( jz +stmz ) + 1
                      j_th  = &
                           ( jx-nxs_th +stmx ) +  &
                           mx_th *( jy-nys_th +stmy ) +  &
                           mxy_th*( jz-nzs_th +stmz ) + 1         
                      s(j)=s(j)+ precon%AS_DDM(myth+1)%svec(j_th)
                   enddo
                endif
             enddo
          enddo
          !$OMP BARRIER          
!          do jz=nzs_th,nze_th-1
          do jz=nzs_th+addL,nze_th-1-addL
             do jy=nys_th,nye_th-1
                if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                   do jx=nxs_th,nxe_th-1
                      j     = &
                           ( jx +stmx ) +  &
                           mx  * ( jy +stmy ) +  &
                           mxy * ( jz +stmz ) + 1
                      j_th  = &
                           ( jx-nxs_th +stmx ) +  &
                           mx_th *( jy-nys_th +stmy ) +  &
                           mxy_th*( jz-nzs_th +stmz ) + 1         
                      ! s(j)=s(j)+ precon%AS_DDM(myth+1)%svec(j_th)
                      ddadd_calc1(s(j),s_lo(j),precon%AS_DDM(myth+1)%svec(j_th),precon%AS_DDM(myth+1)%svec_lo(j_th),s(j),s_lo(j))
                   enddo
                else
                   do jx=nxs_th,nxe_th-1
                      j     = &
                           ( jx +stmx ) +  &
                           mx  * ( jy +stmy ) +  &
                           mxy * ( jz +stmz ) + 1
                      j_th  = &
                           ( jx-nxs_th +stmx ) +  &
                           mx_th *( jy-nys_th +stmy ) +  &
                           mxy_th*( jz-nzs_th +stmz ) + 1         
                      s(j)=s(j)+ precon%AS_DDM(myth+1)%svec(j_th)
                   enddo
                endif
             enddo
          enddo
          !$OMP BARRIER          
!          do jz=nzs_th,nze_th-1
          do jz=nze_th-1-addL+1,nze_th-1
             do jy=nys_th,nye_th-1
                if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                   do jx=nxs_th,nxe_th-1
                      j     = &
                           ( jx +stmx ) +  &
                           mx  * ( jy +stmy ) +  &
                           mxy * ( jz +stmz ) + 1
                      j_th  = &
                           ( jx-nxs_th +stmx ) +  &
                           mx_th *( jy-nys_th +stmy ) +  &
                           mxy_th*( jz-nzs_th +stmz ) + 1         
                      ! s(j)=s(j)+ precon%AS_DDM(myth+1)%svec(j_th)
                      ddadd_calc1(s(j),s_lo(j),precon%AS_DDM(myth+1)%svec(j_th),precon%AS_DDM(myth+1)%svec_lo(j_th),s(j),s_lo(j))
                   enddo
                else
                   do jx=nxs_th,nxe_th-1
                      j     = &
                           ( jx +stmx ) +  &
                           mx  * ( jy +stmy ) +  &
                           mxy * ( jz +stmz ) + 1
                      j_th  = &
                           ( jx-nxs_th +stmx ) +  &
                           mx_th *( jy-nys_th +stmy ) +  &
                           mxy_th*( jz-nzs_th +stmz ) + 1         
                      s(j)=s(j)+ precon%AS_DDM(myth+1)%svec(j_th)
                   enddo
                endif
             enddo
          enddo
       end select
       
       !$OMP END PARALLEL 
#else
       !$OMP PARALLEL &
       !$OMP PRIVATE( &
       !$OMP jz,jy,jx,j,j_th,&
       !$OMP myth,&
       !$OMP m_th,&
       !$OMP mx_th,mxy_th, &
       !$OMP nxs_th,nxe_th,&
       !$OMP nys_th,nye_th,&
       !$OMP nzs_th,nze_th &
       !$OMP ) &
       !$OMP REDUCTION(+:s)

       myth   = omp_get_thread_num()
       call set_AS_region_ddm(&
            precon, &
            myth,   &
            m_th,   &
            mx_th,  &
            mxy_th, &
            nxs_th,nxe_th, &
            nys_th,nye_th, &
            nzs_th,nze_th  &
            )
       !     write(*,*)nys_th,nye_th
       do jz=nzs_th,nze_th-1
          do jy=nys_th,nye_th-1
             do jx=nxs_th,nxe_th-1
                j     = &
                     ( jx +stmx ) +  &
                     mx  * ( jy +stmy ) +  &
                     mxy * ( jz +stmz ) + 1
                j_th  = &
                     ( jx-nxs_th +stmx ) +  &
                     mx_th *( jy-nys_th +stmy ) +  &
                     mxy_th*( jz-nzs_th +stmz ) + 1         
                s(j)=s(j)+ precon%AS_DDM(myth+1)%svec(j_th)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL 
#endif
       
    else
       !$OMP PARALLEL DO PRIVATE(j)
       do j=1,m
          s(j)=0.0d0
       end do

       if(PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
          !$OMP PARALLEL DO PRIVATE(j)
          do j=1,m
             s_lo(j)=0.0d0
          end do
       endif

       !$OMP PARALLEL &
       !$OMP PRIVATE( &
       !$OMP jz,jy,jx,j,j_th,&
       !$OMP myth,&
       !$OMP m_th,&
       !$OMP mx_th,mxy_th, &
       !$OMP nxs_th,nxe_th,&
       !$OMP nys_th,nye_th,&
       !$OMP nzs_th,nze_th, &
       !$OMP nxs_local,nxe_local,&
       !$OMP nys_local,nye_local,&
       !$OMP nzs_local,nze_local &
       !$OMP ) 
       myth   = omp_get_thread_num()
       call set_AS_region_ddm(&
            precon, &
            myth,   &
            m_th,   &
            mx_th,  &
            mxy_th, &
            nxs_th,nxe_th, &
            nys_th,nye_th, &
            nzs_th,nze_th  &
            )

       nxs_local=precon%AS_DDM(myth+1)%nxs_local_th
       nxe_local=precon%AS_DDM(myth+1)%nxe_local_th
       nys_local=precon%AS_DDM(myth+1)%nys_local_th
       nye_local=precon%AS_DDM(myth+1)%nye_local_th
       nzs_local=precon%AS_DDM(myth+1)%nzs_local_th
       nze_local=precon%AS_DDM(myth+1)%nze_local_th
       do jz=nzs_local,nze_local-1
          do jy=nys_local,nye_local-1
             do jx=nxs_local,nxe_local-1
                j     = &
                     ( jx +stmx ) +  &
                     mx  * ( jy +stmy ) +  &
                     mxy * ( jz +stmz ) + 1
#if 0
                j_th  = &
                     ( jx-nxs_th +stmx ) +  &
                     mx_th *( jy-nys_th +stmy ) +  &
                     mxy_th*( jz-nzs_th +stmz ) + 1 
#else
                j_th  = &
                     
                     ( jx - nxs_local + (nxs_local - nxs_th ) + stmx ) +  &
                     mx_th *( jy - nys_local + (nys_local - nys_th ) + stmy ) +  &
                     mxy_th*( jz - nzs_local + (nzs_local - nzs_th ) + stmz ) + 1 
#endif
                s(j) = precon%AS_DDM(myth+1)%svec(j_th)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL 

    endif

    call omp_set_num_threads(max_numth_0)

  end subroutine get_AS_vec_DDM

  subroutine set_AS_region_ddm(&
       precon, &
       myth,   &
       m_th,   &
       mx_th,  &
       mxy_th, &
       nxs_th,nxe_th, &
       nys_th,nye_th, &
       nzs_th,nze_th  &
       )
    use param_AS_type_ddm,only : precon_AS_wrap
    implicit none
    type(precon_AS_wrap) precon
    integer myth
    integer m_th
    integer mx_th
    integer mxy_th
    integer nxs_th,nxe_th
    integer nys_th,nye_th
    integer nzs_th,nze_th
    m_th   = precon%AS_DDM(myth+1)%m_th
    mx_th  = precon%AS_DDM(myth+1)%mx_th
    mxy_th = precon%AS_DDM(myth+1)%mxy_th
    nxs_th = precon%AS_DDM(myth+1)%nxs_th
    nxe_th = precon%AS_DDM(myth+1)%nxe_th
    nys_th = precon%AS_DDM(myth+1)%nys_th
    nye_th = precon%AS_DDM(myth+1)%nye_th
    nzs_th = precon%AS_DDM(myth+1)%nzs_th
    nze_th = precon%AS_DDM(myth+1)%nze_th

  end subroutine set_AS_region_ddm

  subroutine  AdditiveSchwarz_DDM(&
       m,nnd, &
       r, &
       s, &
       prm, &
       precon, &
       ddmDim, &
       ddmOffset, &
       ddmA, &
       ioff_stencil, &
       addL, &
       PRECISION_CALC, &
       PRECISION_PRECONDITION, &
       ddmA_lo, &
       r_lo, &
       s_lo  &
       )
    use param_AS_type_ddm,only : precon_AS_wrap
    use omp_lib
    use vec_kernel
    use ddmMat_wrap
    use PARCEL_time,only : precon_ddm_time
    implicit none
    !      include "ddm.h"
    integer(kind=intType) m
    integer nnd
    real*8 r(m)
    real*8 s(m)
    type(ddm_param) prm
    type(precon_AS_wrap) precon
    integer ddmDim
    integer ddmOffset(nnd)
    real*8 ddmA(m*nnd)
    integer ioff_stencil(nnd,ddmDim)
    integer(kind=intType) i
    integer addL
    real*8 stime,etime

    integer,optional :: PRECISION_CALC
    integer,optional :: PRECISION_PRECONDITION
    real*8,optional ::  ddmA_lo(m*nnd)
    real*8,optional ::  r_lo(m)
    real*8,optional ::  s_lo(m)
    integer PRECISION
    PRECISION = KSP_DOUBLE_PRECISION
    if(  present(PRECISION_CALC)  )then
       PRECISION = PRECISION_CALC
    endif

#if 1
       stime=MPI_Wtime()
       if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
          call set_AS_vec_DDM(&
               prm,precon,m,r, &
               PRECISION_CALC, &
               r_lo &
               )
       else
          call set_AS_vec_DDM(&
               prm,precon,m,r &
               )
       endif
       etime=MPI_Wtime()
       precon_ddm_time(11) = precon_ddm_time(11) + etime - stime
#endif

#if 1
       stime=MPI_Wtime()
       if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then

          call AdditiveSchwarz_ILU_DDM(&
               m,nnd,&
!               r, &
               prm,precon,&
               ddmDim, &
               ddmOffset, &
               ddmA, &
               ioff_stencil, &
               PRECISION, &
               PRECISION_PRECONDITION, &
               ddmA_lo &
               )
       else
          call AdditiveSchwarz_ILU_DDM(&
               m,nnd,&
!               r, &
               prm,precon,&
               ddmDim, &
               ddmOffset, &
               ddmA, &
               ioff_stencil &
               )
       endif
       etime=MPI_Wtime()
       precon_ddm_time(12) = precon_ddm_time(12) + etime - stime
#endif
    !  write(*,*)" @ 3 @ "
    ! call set_AS_METHOD_TYPE_DDM(precon,KSP_ADDITIVESCHWARZ_TYPE_RESTRICT)
    !  call set_AS_METHOD_TYPE_DDM(precon,KSP_ADDITIVESCHWARZ_TYPE_BASIC)
    !  write(*,*)" @ 4 @ "
#if 1
       stime=MPI_Wtime()
       if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
          call get_AS_vec_DDM(&
               prm,precon,m,s,addL, &
               PRECISION,s_lo &
               )
       else
          call get_AS_vec_DDM(&
               prm,precon,m,s,addL &
               )
       endif
       etime=MPI_Wtime()
       precon_ddm_time(13) = precon_ddm_time(13) + etime - stime
#endif    

    !   write(*,*)" @ 5 @ "
    !   do i=1,m
    !      write(*,*)i,s(i)
    !   enddo
  end subroutine AdditiveSchwarz_DDM

  subroutine set_AS_DDM(&
       addL, &
       iflagAS, &
       prm,prm_g, &
       precon, &
       m, &
       nnd, &
       ddmOffset, &
       ddmA, &
       ioff_stencil, &
       ddmDim, &
       PRECISION_CALC, &
       PRECISION_A, &
       PRECISION_PRECONDITION, &
       ddmA_lo &
       )
    use param_AS_type_ddm,only : precon_AS_wrap, &
         KSP_DDM_DIVIDE_THREADS_DIRECT_X, &
         KSP_DDM_DIVIDE_THREADS_DIRECT_Y, &
         KSP_DDM_DIVIDE_THREADS_DIRECT_Z

    use ddmMat_wrap
    implicit none
    !     include "ddm.h"

    integer addL
    integer iflagAS
    type(ddm_param) :: prm
    type(AS_geom) :: prm_g ! 計算領域のオーバーラップに関する幾何情報 --> これは不要?
    type(precon_AS_wrap) precon
    integer precon_thblock
    integer(kind=intType) m
    integer nnd
    integer ddmDim
    integer ddmOffset(nnd)
    real*8 ddmA(m*nnd)
    integer ioff_stencil(nnd,ddmDim)

    integer nx,ny,nz
    integer div_direc_th
    integer nrk_xp,nrk_xm
    integer nrk_yp,nrk_ym
    integer nrk_zp,nrk_zm
    integer nrk_th_direc_p,nrk_th_direc_m

    integer,optional :: PRECISION_CALC
    integer,optional :: PRECISION_A
    integer,optional :: PRECISION_PRECONDITION
    real*8,optional  :: ddmA_lo(m*nnd)
    integer PRECISION

    PRECISION = KSP_DOUBLE_PRECISION
    if(  present(PRECISION_CALC)  )then
       PRECISION = PRECISION_CALC
    endif

    nx=prm%nx
    ny=prm%ny
    nz=prm%nz
    div_direc_th=prm%div_direc_th
    precon_thblock=prm%precon_thblock
    nrk_xm=prm_g%rank_s(1)
    nrk_xp=prm_g%rank_s(2)
    nrk_ym=prm_g%rank_s(3)
    nrk_yp=prm_g%rank_s(4)
    nrk_zm=prm_g%rank_s(5)
    nrk_zp=prm_g%rank_s(6)

    nrk_th_direc_m=-1
    nrk_th_direc_p=-1
    select case(div_direc_th)
    case(KSP_DDM_DIVIDE_THREADS_DIRECT_X)
       nrk_th_direc_m=prm_g%rank_s(1)
       nrk_th_direc_p=prm_g%rank_s(2)
    case(KSP_DDM_DIVIDE_THREADS_DIRECT_Y)
       nrk_th_direc_m=prm_g%rank_s(3)
       nrk_th_direc_p=prm_g%rank_s(4)
    case(KSP_DDM_DIVIDE_THREADS_DIRECT_Z)
       nrk_th_direc_m=prm_g%rank_s(5) 
       nrk_th_direc_p=prm_g%rank_s(6) 
    end select

    ! write(*,*)" m = ",m
    ! write(*,*)"prm_g%rank_s(1)",prm_g%rank_s(1)
    ! write(*,*)"prm_g%rank_s(2)",prm_g%rank_s(2)
    ! write(*,*)"prm_g%rank_s(3)",prm_g%rank_s(3)
    ! write(*,*)"prm_g%rank_s(4)",prm_g%rank_s(4)
    ! write(*,*)"prm_g%rank_s(5)",prm_g%rank_s(5)
    ! write(*,*)"prm_g%rank_s(6)",prm_g%rank_s(6)
    ! write(*,*)"ddmOffset" ,ddmOffset

    !     write(*,*)" nx, ny, nz", nx, ny, nz
    if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       call get_AS_region_DDM(&
            nx,ny,nz,&
            div_direc_th,&
            nrk_th_direc_p,nrk_th_direc_m, &
            addL,&
            precon,precon_thblock, &
            PRECISION_CALC,&
            PRECISION_PRECONDITION)
    else
       call get_AS_region_DDM(&
            nx,ny,nz,&
            div_direc_th,&
            nrk_th_direc_p,nrk_th_direc_m, &
            addL,&
            precon,precon_thblock)
    endif
    call set_AS_METHOD_TYPE_DDM(precon,iflagAS)

    if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       if( PRECISION_A.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
          call init_AS_DILU_DDM(&
               m,&
               nnd,&
               ddmOffset,&
               ddmA,&
               ioff_stencil, &
               prm,&
               ddmDim, &
               precon, &
               PRECISION_CALC , &
               PRECISION_A, &
               PRECISION_PRECONDITION, &
               ddmA_lo )
       else
          call init_AS_DILU_DDM(&
               m,&
               nnd,&
               ddmOffset,&
               ddmA,&
               ioff_stencil, &
               prm,&
               ddmDim, &
               precon, &
               PRECISION_PRECONDITION)
       endif
    else
       call init_AS_DILU_DDM(&
            m,&
            nnd,&
            ddmOffset,&
            ddmA,&
            ioff_stencil, &
            prm,&
            ddmDim, &
            precon)
    endif
  end subroutine set_AS_DDM

end module param_AS_ddm
