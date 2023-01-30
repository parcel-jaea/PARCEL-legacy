#include "ddfunc.h"
! #define CUDA_AWARE 1

module ddm_matv
  use PARCEL_define,only : KSP_DOUBLE_DOUBLE_PRECISION,KSP_DOUBLE_PRECISION
contains

  subroutine matv_ddm( &
       NEWCOMM, NND, mdist, A, IOFF_A, P, Q, &
       myrank, NPES, prm, commHide, &
       PRECISION_CALC, &
       PRECISION_A, &
       A_lo,&
       P_lo, Q_lo)
    ! use ddm_tmp
    use ddmMat_wrap
    implicit none
    !
    integer(kind=intType), intent(in) :: mdist
    integer, intent(in) :: newcomm, nnd
    integer, intent(in) :: npes, myrank
    integer, intent(in) :: ioff_a(nnd)
    real*8, intent(in) ::  A(nnd*mdist)
    real*8 P(mdist)
    real*8, intent(inout) :: Q(mdist)
    type(ddm_param), intent(in) :: prm
    integer, intent(in) :: commHide
    integer,optional :: PRECISION_CALC
    integer,optional :: PRECISION_A
    real*8,optional  :: A_lo(nnd*mdist)
    real*8,optional :: P_lo(mdist)
    real*8,optional :: Q_lo(mdist)

    !
    integer(kind=intType) i
    real*8, allocatable :: vec(:)
    real*8, allocatable :: vec_lo(:)
    logical wideSleeve
    integer PRECISION

    PRECISION = KSP_DOUBLE_PRECISION
    if(  present(PRECISION_CALC)  )then
       PRECISION = PRECISION_CALC
    endif

   ! !$OMP PARALLEL do
   ! Do i=1, mdist
   !    Q(i) = 0.0d0
   ! End Do
   ! !$OMP end PARALLEL do

    !
    ! プロセス数が1の時の処理
    !
    if ( npes == 1 .or. prm%leng_halo == 0 ) then
       if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
!          write(*,*)"dd spaxv_ddm_commhide_0"
          call spaxv_ddm_commhide_0( &
               NND, mdist, A, IOFF_A, P, Q, prm, &
               PRECISION_CALC, &
               PRECISION_A,    &
               A_lo,&
               P_lo,Q_lo)
       else
          call spaxv_ddm_commhide_0( NND, mdist, A, IOFF_A, P, Q, prm )    
       endif
       ! ↑通信隠蔽を行わずに Q = A * P を計算する

       return

    endif


    allocate( vec(mdist) )
    !$OMP PARALLEL do
    Do i=1, mdist
       vec(i) = P(i)
    End do
    !$OMP end PARALLEL do
    if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
       allocate( vec_lo(mdist) )
       !$OMP PARALLEL do
       Do i=1, mdist
          vec_lo(i) = P_lo(i)
       End do
       !$OMP end PARALLEL do
    endif

    if ( commHide==0 ) then
       !
       ! 通信隠蔽なし
       !
       call sleev_comm_novec2( mdist, vec, prm, myrank, newcomm )
       if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
          call sleev_comm_novec2( mdist, vec_lo, prm, myrank, newcomm )
       endif
       ! ↑袖領域のデータを（袖領域と重なる領域を担当するプロセスから）取得
       !
!       call spaxv_ddm_commhide_0( NND, mdist, A, IOFF_A, vec, Q, prm )    
       if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
          call spaxv_ddm_commhide_0( &
               NND, mdist, A, IOFF_A, vec, Q, prm, &
               PRECISION_CALC, &
               PRECISION_A,    &
               A_lo,&
               vec_lo,Q_lo)
       else
          call spaxv_ddm_commhide_0( NND, mdist, A, IOFF_A, vec, Q, prm )    
       endif

       ! ↑ Q = A * vec を計算
       !
    else if ( commHide == 1 ) then
       !
       ! 通信隠蔽（一括通信）
       !
       write(*,*)"not support : matv_ddm, ovl = ",commHide

       ! if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
       !    call spaxv_ddm_commhide_1( NND, mdist, A, IOFF_A, vec, Q, prm, myrank, newcomm ,&
       !         PRECISION_CALC, &
       !         PRECISION_A,    &
       !         A_lo,&
       !         vec_lo,Q_lo)
       ! else
       !    call spaxv_ddm_commhide_1( NND, mdist, A, IOFF_A, vec, Q, prm, myrank, newcomm )
       ! endif
       !
    else if ( commHide == 2 ) then
       !
       ! 通信隠蔽（小分け通信）
       !
       write(*,*)"not support : matv_ddm, ovl = ",commHide
       if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
          ! call spaxv_ddm_commhide_2( NND, mdist, A, IOFF_A, vec, Q, prm, myrank, newcomm, &
          !      PRECISION_CALC, &
          !      PRECISION_A,    &
          !      A_lo,&
          !      vec_lo,Q_lo)
       else
         !  call spaxv_ddm_commhide_2( NND, mdist, A, IOFF_A, vec, Q, prm, myrank, newcomm )    
       endif
       !
    endif

    deallocate( vec )
    if(allocated(vec_lo))deallocate( vec_lo )

  end subroutine matv_ddm


  subroutine matv_ddm_GPU( &
       method,&
       x, offset_x,&
       y, offset_y)
    ! use ddm_tmp
    use krylov
    use ddmMat_wrap
    use cu_Matrix,only : Matvec_ddm_FP64_GPU
    use PARCEL_VEC_kernel,only : PARCEL_VEC
    implicit none
    type(KSP) :: method
!    type(PARCEL_VEC) :: x
!    type(PARCEL_VEC) :: y
    type(devicemem) :: x
    type(devicemem) :: y
    integer(kind=intType) offset_x
    integer(kind=intType) offset_y

    integer(kind=intType) m
    integer(kind=intType) nx,ny,nz
    integer(kind=intType) mx,my,mxy
    integer(kind=intType) stmx,stmy,stmz
    integer(kind=intType) nnd
    integer(kind=intType) npe

    m    = method%ddmMat%n
    nx   = method%ddmMat%prm%nx
    ny   = method%ddmMat%prm%ny
    nz   = method%ddmMat%prm%nz
  
    mx   = method%ddmMat%prm%mx
    my   = method%ddmMat%prm%my
    mxy  = method%ddmMat%prm%mxy

    stmx = method%ddmMat%prm%stmx
    stmy = method%ddmMat%prm%stmy
    stmz = method%ddmMat%prm%stmz

    nnd = method%ddmMat%mat%nnd_dia
    npe = method%ddmMat%ncpu

    ! write(*,*) m
    ! write(*,*) nx
    ! write(*,*) ny
    ! write(*,*) nz
    ! write(*,*) mx
    ! write(*,*) my
    ! write(*,*) mxy
    ! write(*,*) stmx
    ! write(*,*) stmy
    ! write(*,*) stmz
    ! write(*,*) nnd
    ! write(*,*) npe
    ! write(*,*) " -------------- "
    ! write(*,*) 

    if ( npe .ne. 1 ) then
#if 0
       call cudaMemcpy_d2d_fp64( &
            x, &
            method%ddmMat%C_work,m)
       call sleev_comm_GPU( method, method%ddmMat%C_work )
       call Matvec_ddm_FP64_GPU(&
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            nx, &
            ny, &
            nz, &
            stmx, &
            stmy, &
            stmz, &
            mx, &
            mxy, &
            m, &
            nnd, &
            method%ddmMat%mat%C_val_dia, &
            method%ddmMat%C_work,&
            y, &		
            method%ddmMat%mat%C_ioff_dia &
            )
#else
       call sleev_comm_GPU(method,x,offset_x)
       call Matvec_ddm_FP64_GPU(&
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            nx, &
            ny, &
            nz, &
            stmx, &
            stmy, &
            stmz, &
            mx, &
            mxy, &
            m, &
            nnd, &
            method%ddmMat%mat%C_val_dia, &
            x, offset_x, &
            y, offset_y, &
            method%ddmMat%mat%C_ioff_dia &
            )
       call sleev_zero_GPU(method,x,offset_x)
#endif
    else
       call Matvec_ddm_FP64_GPU(&
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            nx, &
            ny, &
            nz, &
            stmx, &
            stmy, &
            stmz, &
            mx, &
            mxy, &
            m, &
            nnd, &
            method%ddmMat%mat%C_val_dia, &
            x, offset_x, &
            y, offset_y, &
            method%ddmMat%mat%C_ioff_dia &
            )
    endif

    return    
  end subroutine matv_ddm_GPU

  subroutine matv_Neumann_ddm_GPU( &
       method,&
       x, offset_x,&
       y, offset_y,&
       nocomm)
    ! use ddm_tmp
    use krylov
    use ddmMat_wrap
    use cu_Matrix,only : Matvec_ddm_FP64_GPU
    use PARCEL_VEC_kernel,only : PARCEL_VEC
    implicit none
    type(KSP) :: method
!    type(PARCEL_VEC) :: x
!    type(PARCEL_VEC) :: y
    type(devicemem) :: x
    integer(kind=intType) offset_x
    type(devicemem) :: y
    integer(kind=intType) offset_y

    integer(kind=intType) m
    integer(kind=intType) nx,ny,nz
    integer(kind=intType) mx,my,mxy
    integer(kind=intType) stmx,stmy,stmz
    integer(kind=intType) nnd
    integer(kind=intType) npe
    integer,optional :: nocomm

    m    = method%ddmMat%n
    nx   = method%ddmMat%C_DDMMat_Neumann%nx
    ny   = method%ddmMat%C_DDMMat_Neumann%ny
    nz   = method%ddmMat%C_DDMMat_Neumann%nz
  
    mx   = method%ddmMat%C_DDMMat_Neumann%mx
    my   = method%ddmMat%C_DDMMat_Neumann%my
    mxy  = method%ddmMat%C_DDMMat_Neumann%mxy

    stmx = method%ddmMat%C_DDMMat_Neumann%stmx
    stmy = method%ddmMat%C_DDMMat_Neumann%stmy
    stmz = method%ddmMat%C_DDMMat_Neumann%stmz

    nnd = method%ddmMat%C_DDMMat_Neumann%nnd
    npe = method%ddmMat%ncpu

    if ( npe .ne. 1 ) then
       if(  present(nocomm)  )then
          call sleev_zero_GPU(method,x,offset_x)
       else
          call sleev_comm_GPU(method,x,offset_x)
       endif
       call Matvec_ddm_FP64_GPU(&
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            nx, &
            ny, &
            nz, &
            stmx, &
            stmy, &
            stmz, &
            mx, &
            mxy, &
            m, &
            nnd, &
            method%ddmMat%C_DDMMat_Neumann%ddmA, &
            x, offset_x,&
            y, offset_y,&
            method%ddmMat%C_DDMMat_Neumann%offset &
            )
       call sleev_zero_GPU(method,x,offset_x)
    else
       call Matvec_ddm_FP64_GPU(&
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            nx, &
            ny, &
            nz, &
            stmx, &
            stmy, &
            stmz, &
            mx, &
            mxy, &
            m, &
            nnd, &
            method%ddmMat%C_DDMMat_Neumann%ddmA, &
            x, offset_x,&
            y, offset_y,&		
            method%ddmMat%C_DDMMat_Neumann%offset &
            )
    endif

    return    
  end subroutine matv_Neumann_ddm_GPU

  subroutine jacobi_iter1_ddm_GPU( &
       method,&
       x, offset_x,&
       y, offset_y, &
       nocomm)
    ! use ddm_tmp
    use krylov
    use ddmMat_wrap
    use cu_Matrix,only : jacobi_iter1_ddm_FP64_GPU
    use PARCEL_VEC_kernel,only : PARCEL_VEC
    implicit none
    type(KSP) :: method
!    type(PARCEL_VEC) :: x
!    type(PARCEL_VEC) :: y
    type(devicemem) :: x
    integer(kind=intType) offset_x
    type(devicemem) :: y
    integer(kind=intType) offset_y
    integer,optional :: nocomm

    integer(kind=intType) m
    integer(kind=intType) nx,ny,nz
    integer(kind=intType) mx,my,mxy
    integer(kind=intType) stmx,stmy,stmz
    integer(kind=intType) nnd
    integer(kind=intType) npe
    integer(kind=intType) nnd_

    m    = method%ddmMat%n
    nx   = method%ddmMat%prm%nx
    ny   = method%ddmMat%prm%ny
    nz   = method%ddmMat%prm%nz
  
    mx   = method%ddmMat%prm%mx
    my   = method%ddmMat%prm%my
    mxy  = method%ddmMat%prm%mxy

    stmx = method%ddmMat%prm%stmx
    stmy = method%ddmMat%prm%stmy
    stmz = method%ddmMat%prm%stmz

    nnd = method%ddmMat%mat%nnd_dia
    npe = method%ddmMat%ncpu

    if ( npe .ne. 1 ) then
       if(  present(nocomm)  )then
          call sleev_zero_GPU(method,x,offset_x)
       else
          call sleev_comm_GPU(method,x,offset_x)
       endif
       nnd_ = method%ddmMat%C_DDMMat_Neumann%nnd
       call jacobi_iter1_ddm_FP64_GPU(&
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            nx, &
            ny, &
            nz, &
            stmx, &
            stmy, &
            stmz, &
            mx, &
            mxy, &
            m, &
            nnd_, &
            method%ddmMat%C_DDMMat_Neumann%ddmA, &
            x, offset_x,&
            y, offset_y,&		
            method%ddmMat%C_DDMMat_Neumann%offset &
            )
       call sleev_zero_GPU(method,x,offset_x)
    else
       nnd_ = method%ddmMat%C_DDMMat_Neumann%nnd
       call jacobi_iter1_ddm_FP64_GPU(&
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            nx, &
            ny, &
            nz, &
            stmx, &
            stmy, &
            stmz, &
            mx, &
            mxy, &
            m, &
            nnd_, &
            method%ddmMat%C_DDMMat_Neumann%ddmA, &
            x, offset_x,&
            y, offset_y,&		
            method%ddmMat%C_DDMMat_Neumann%offset &
            )
    endif

    return    
  end subroutine jacobi_iter1_ddm_GPU

  subroutine spaxv_ddm_commhide_0( &
       NND, mdist, A, IOFF_A, P, Q, prm , &
       PRECISION_CALC, &
       PRECISION_A,&
       A_lo,&
       P_lo,Q_lo)

    use omp_lib
    ! use ddm_tmp
    use ddmMat_wrap
    use common_kernel,only : set_region_for_threads
    implicit none
    !
    integer(kind=intType), intent(in) ::  mdist
    integer, intent(in) :: nnd
    real*8, intent(in) ::  A(nnd*mdist)
    integer, intent(in) :: ioff_a(nnd)
    real*8  P(mdist)
    real*8, intent(inout) :: Q(mdist)
    type(ddm_param) :: prm
    integer,optional :: PRECISION_CALC
    integer,optional :: PRECISION_A
    real*8,optional,intent(in)  :: A_lo(nnd*mdist)
    real*8,optional :: P_lo(mdist)
    real*8,optional :: Q_lo(mdist)

    !
    integer(kind=intType) j, jx, jy, jz, ip
    integer nx, ny, nz, nxy, mx, my, mxy
    integer stmx, stmy, stmz
    integer PRECISION
    real*8 ddmul_var1
    real*8 ddadd_var1
    real*8 c_hi,c_lo

    PRECISION = KSP_DOUBLE_PRECISION
    if(  present(PRECISION_CALC)  )then
       PRECISION = PRECISION_CALC
    endif

    nx=prm%nx;  ny=prm%ny; nz=prm%nz;  nxy=prm%nxy;
    mx=prm%mx;  my=prm%my; mxy=prm%mxy;
    stmx=prm%stmx;  stmy=prm%stmy; stmz=prm%stmz;
    
    if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then

       if( present(PRECISION_A) )then
          if( PRECISION_A.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
             !$OMP parallel do private(jx,jy,jz,ip,j,&
             !$OMP ddmul_var1, &
             !$OMP ddadd_var1, &
             !$OMP c_hi,c_lo &
             !$OMP ) schedule(dynamic, 1) 
             do jz=0,nz-1
                do jy=0,ny-1
                   do jx=0,nx-1
                      j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1
                      Q(j)    = 0.0d0
                      Q_lo(j) = 0.0d0
                   enddo
                   do ip=1,nnd
                      do jx=0,nx-1
                         j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1
                         !                   Q(j) = Q(j) + A( mdist*(ip-1) +j ) *P( j +ioff_A(ip) )
                       !  ddmul_calc1(A(mdist*(ip-1)+j),A_lo(mdist*(ip-1)+j),P(j+ioff_A(ip)),P_lo(j+ioff_A(ip)),c_hi,c_lo)
                       !  ddadd_calc1(c_hi,c_lo,Q(j),Q_lo(j),Q(j),Q_lo(j))

                         ddmul_calc1(A(mdist*(ip-1)+j),A_lo(mdist*(ip-1)+j),P(j+ioff_A(ip)),P_lo(j+ioff_A(ip)),c_hi,c_lo)
                         ddadd_calc1(c_hi,c_lo,Q(j),Q_lo(j),Q(j),Q_lo(j))

                      end do
                   end do
                end do
             end do

             return             
          endif
       endif

       !$OMP parallel do private(jx,jy,jz,ip,j,&
       !$OMP ddmul_var1, &
       !$OMP ddadd_var1, &
       !$OMP c_hi,c_lo &
       !$OMP ) schedule(dynamic, 1) 
       do jz=0,nz-1
          do jy=0,ny-1
             do jx=0,nx-1
                j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1
                Q(j)    = 0.0d0
                Q_lo(j) = 0.0d0
             enddo
             do ip=1,nnd
                do jx=0,nx-1
                   j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1
                   !                   Q(j) = Q(j) + A( mdist*(ip-1) +j ) *P( j +ioff_A(ip) )
                   ddmul_calc1(A(mdist*(ip-1)+j),0.0d0,P(j+ioff_A(ip)),P_lo(j+ioff_A(ip)),c_hi,c_lo)
                   ddadd_calc1(c_hi,c_lo,Q(j),Q_lo(j),Q(j),Q_lo(j))
                end do
             end do
          end do
       end do
    else
       !$OMP parallel do private(jx,jy,jz,ip,j) schedule(dynamic, 1) 
       do jz=0,nz-1
          do jy=0,ny-1
             do jx=0,nx-1
                j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1
                Q(j) = 0.0d0
             enddo
             do ip=1,nnd
                do jx=0,nx-1
                   j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1
                   Q(j) = Q(j) + A( mdist*(ip-1) +j ) *P( j +ioff_A(ip) )
                end do
             end do
          end do
       end do
    endif

  end subroutine spaxv_ddm_commhide_0

  !!
  !! commhide_1
  !!
  subroutine spaxv_ddm_commhide_1( NND, mdist, A, IOFF_A, P, Q, prm, myrank, newcomm,&
       PRECISION_CALC, &
       PRECISION_A,&
       A_lo,&
       P_lo,Q_lo)
    ! use ddm_tmp
    use ddmMat_wrap
    implicit none
    !
    integer(kind=intType), intent(in) ::  mdist
    integer, intent(in) :: nnd
    real*8, intent(in) ::  A(nnd*mdist)
    integer, intent(in) :: ioff_a(nnd)
    real*8 P(mdist)
    real*8, intent(inout) :: Q(mdist)
    type(ddm_param), intent(in) :: prm
    integer, intent(in) :: newcomm, myrank
    integer,optional :: PRECISION_CALC
    integer,optional :: PRECISION_A
    real*8,optional,intent(in)  :: A_lo(nnd*mdist)
    real*8,optional :: P_lo(mdist)
    real*8,optional :: Q_lo(mdist)

    !
    integer ip
    integer j, jx, jy, jz
    integer nx, ny, nz, nxy, mx, my, mxy
    integer stmx, stmy, stmz
    integer div_direc_th
    integer halox, haloy, haloz
    integer ncomp_grids
    !
    integer i, offset
    integer nx_center, ny_center, nz_center, nxy_center ! 通信なしで計算できる部分のグリッド数
    integer size_center, size_margin ! 中心部のグリッド数、残りのグリッド数
    integer offset_bottom, offset_top, size_mid
    integer offset_bottom_xy, offset_top_xy, size_mid_xy
    logical is_center
    real*8 ddmul_var1
    real*8 ddadd_var1
    real*8 c_hi,c_lo
    integer PRECISION
    PRECISION = KSP_DOUBLE_PRECISION
    if(  present(PRECISION_CALC)  )then
       PRECISION = PRECISION_CALC
    endif

    nx=prm%nx;  ny=prm%ny; nz=prm%nz;  nxy=prm%nxy;
    mx=prm%mx;  my=prm%my; mxy=prm%mxy;
    stmx=prm%stmx;  stmy=prm%stmy; stmz=prm%stmz;
    div_direc_th = prm%div_direc_th;
    ncomp_grids = prm%ncomp_grids
    !
    halox = prm%leng_halo;  haloy = 0;  haloz = 0
    if ( prm%ncomp_grids >= 2 ) haloy = prm%leng_halo
    if ( prm%ncomp_grids == 3 ) haloz = prm%leng_halo
    nx_center = max(0, nx-2*halox)
    ny_center = max(0, ny-2*haloy)
    nz_center = max(0, nz-2*haloz)
    nxy_center = nx_center * ny_center
    size_center = nx_center * ny_center * nz_center
    size_margin = nx*ny*nz - size_center
    offset_bottom = nx*ny*haloz
    offset_top = size_margin - offset_bottom
    size_mid = nxy - nxy_center
    offset_bottom_xy = nx*haloy
    offset_top_xy = nx * ny - nx_center * ny_center - offset_bottom_xy
    size_mid_xy = nx - nx_center


    !$OMP parallel private( jz, jy, jx, j, ip, i, offset, is_center, &
    !$OMP ddmul_var1, &
    !$OMP ddadd_var1, &
    !$OMP c_hi,c_lo &
    !$OMP )

    !$OMP master
    !
    ! 袖領域を取得する
    !
    call sleev_comm_novec2( mdist, P, prm, myrank, newcomm )
    
    if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
       call sleev_comm_novec2( mdist, P_lo, prm, myrank, newcomm )
    endif
    !
    !$OMP end master

    !
    ! 行列AとベクトルPの積を計算（通信なしで計算できる部分のみ）
    ! 後でマスタースレッドが合流する
    !
    if (prm%ncomp_grids == 3) then
       !$OMP do schedule(dynamic, 1) 
       !
       do jz=haloz, haloz+nz_center-1
          if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
             do jy=haloy, haloy+ny_center-1
                do jx=halox, halox+nx_center-1
                   j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1
                   Q(j)    = 0.0d0
                   Q_lo(j) = 0.0d0
                enddo
             enddo
             if( present(PRECISION_A) )then
                if( PRECISION_A.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                   do ip=1, nnd
                      do jy=haloy, haloy+ny_center-1
                         do jx=halox, halox+nx_center-1
                            j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1
                            !  Q(j) = Q(j) +A( mdist*(ip-1) +j ) *P( j +ioff_A(ip) )
                            ddmul_calc1(A(mdist*(ip-1)+j),0.0d0,P(j+ioff_A(ip)),P_lo(j+ioff_A(ip)),c_hi,c_lo)
                            ddadd_calc1(c_hi,c_lo,Q(j),Q_lo(j),Q(j),Q_lo(j))
                         end do
                      end do
                   end do                   
                   cycle
                endif
             endif
             do ip=1, nnd
                do jy=haloy, haloy+ny_center-1
                   do jx=halox, halox+nx_center-1
                      j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1
                      !  Q(j) = Q(j) +A( mdist*(ip-1) +j ) *P( j +ioff_A(ip) )
                      ddmul_calc1(A(mdist*(ip-1)+j),A_lo( mdist*(ip-1) +j ),P(j+ioff_A(ip)),P_lo(j+ioff_A(ip)),c_hi,c_lo)
                      ddadd_calc1(c_hi,c_lo,Q(j),Q_lo(j),Q(j),Q_lo(j))
                   end do
                end do
             end do
          else
             do jy=haloy, haloy+ny_center-1
                do jx=halox, halox+nx_center-1
                   j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1
                   Q(j)    = 0.0d0
                enddo
             enddo
             do ip=1, nnd
                do jy=haloy, haloy+ny_center-1
                   do jx=halox, halox+nx_center-1
                      j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1
                      Q(j) = Q(j) +A( mdist*(ip-1) +j ) *P( j +ioff_A(ip) )
                   end do
                end do
             end do
          endif
       end do
       !
       !$OMP end do
    else
       !$OMP do schedule(dynamic, 1) 
       !
       do jz=haloz, haloz+nz_center-1
          if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
             do jy=haloy, haloy+ny_center-1
                do jx=halox, halox+nx_center-1
                   j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1
                   Q(j) = 0.0d0
                   Q_lo(j) = 0.0d0
                end do
             end do
             if( present(PRECISION_A) )then
                if( PRECISION_A.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                   do jy=haloy, haloy+ny_center-1
                      do ip=1, nnd
                         do jx=halox, halox+nx_center-1
                            j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1
                            ! Q(j) = Q(j) +A( mdist*(ip-1) +j ) *P( j +ioff_A(ip) )
                            ddmul_calc1(A(mdist*(ip-1)+j),A_lo(mdist*(ip-1)+j),P(j+ioff_A(ip)),P_lo(j+ioff_A(ip)),c_hi,c_lo)
                            ddadd_calc1(c_hi,c_lo,Q(j),Q_lo(j),Q(j),Q_lo(j))
                         end do
                      end do
                   end do
                   cycle
                endif
             endif
             do jy=haloy, haloy+ny_center-1
                do ip=1, nnd
                   do jx=halox, halox+nx_center-1
                      j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1
                      ! Q(j) = Q(j) +A( mdist*(ip-1) +j ) *P( j +ioff_A(ip) )
                      ddmul_calc1(A(mdist*(ip-1)+j),0.0d0,P(j+ioff_A(ip)),P_lo(j+ioff_A(ip)),c_hi,c_lo)
                      ddadd_calc1(c_hi,c_lo,Q(j),Q_lo(j),Q(j),Q_lo(j))
                   end do
                end do
             end do
          else
             do jy=haloy, haloy+ny_center-1
                do jx=halox, halox+nx_center-1
                   j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1
                   Q(j) = 0.0d0
                end do
             end do
             do jy=haloy, haloy+ny_center-1
                do ip=1, nnd
                   do jx=halox, halox+nx_center-1
                      j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1
                      Q(j) = Q(j) +A( mdist*(ip-1) +j ) *P( j +ioff_A(ip) )
                   end do
                end do
             end do
          endif
       end do
       !
       !$OMP end do
    end if

    ! OMP doには暗黙のバリアがあるので、この段階で全てのスレッドが並んでいる

    !
    ! 行列AとベクトルPの積を計算（通信が必要な部分）
    !
    !$OMP do schedule(dynamic, 1) collapse(2)
    !
    do jz=0, nz-1
       do jy=0, ny-1
          is_center = .false.
          if (haloz<=jz .and. jz<nz-haloz &
               .and. haloy<=jy .and. jy<ny-haloy) then
             is_center = .true.
          end if
          jx = 0
          do while (jx < nx)
             j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1;
             
             if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                if( present(PRECISION_A) )then
                   if( PRECISION_A.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                      do ip=1, nnd
                         ! Q(j) = Q(j) +A( mdist*(ip-1) +j ) *P( j +ioff_A(ip) )
                         ddmul_calc1(A( mdist*(ip-1) +j ) ,A_lo( mdist*(ip-1) +j ) ,P( j +ioff_A(ip) ),P( j +ioff_A(ip) ),c_hi,c_lo)
                         ddadd_calc1(c_hi,c_lo,Q(j),Q_lo(j),Q(j),Q_lo(j))
                      end do

                      jx = jx + 1
                      if (jx==halox .and. is_center) then
                         jx = jx + nx_center
                      end if
                      cycle
                   endif
                endif
                do ip=1, nnd
                   ! Q(j) = Q(j) +A( mdist*(ip-1) +j ) *P( j +ioff_A(ip) )
                   ddmul_calc1(A( mdist*(ip-1) +j ) , 0.0d0 ,P( j +ioff_A(ip) ),P( j +ioff_A(ip) ),c_hi,c_lo)
                   ddadd_calc1(c_hi,c_lo,Q(j),Q_lo(j),Q(j),Q_lo(j))
                end do
                jx = jx + 1
                if (jx==halox .and. is_center) then
                   jx = jx + nx_center
                end if
             else
                do ip=1, nnd
                   Q(j) = Q(j) +A( mdist*(ip-1) +j ) *P( j +ioff_A(ip) )
                end do
                jx = jx + 1
                if (jx==halox .and. is_center) then
                   jx = jx + nx_center
                end if
             endif
          end do
       end do
    end do
    !
    !$OMP end do

    !$OMP end PARALLEL

  end subroutine spaxv_ddm_commhide_1


  !!
  !! commhide_2
  !!
  subroutine spaxv_ddm_commhide_2( NND, mdist, A, IOFF_A, P, Q, prm, myrank, newcomm )
    ! use ddm_tmp
    use ddmMat_wrap
    use mpi
    implicit none
!    include "mpif.h"
    !
    integer, intent(in) :: nnd, mdist
    real*8, intent(in) ::  A(nnd*mdist)
    integer, intent(in) :: ioff_a(nnd)
    real*8 P(mdist)
    real*8, intent(inout) :: Q(mdist)
    type(ddm_param), intent(in) :: prm
    integer, intent(in) :: newcomm, myrank
    !
    integer j, jx, jy, jz, i, k, ip, ix, iy, iz
    integer kx, ky, kz, kxmax, kymax, kzmax
    integer halox, haloy, haloz, info
    !
    type array
       real*8, allocatable :: v(:)
    end type array
    type(array) sbuf(6), rbuf(6) ! 6 = ±x, ±y, ±z   
    integer :: ireqs(6), ireqr(6)
    integer istatus (mpi_status_size)
    integer xyz, pm
    !
    ! A*Pの計算領域
    !
    integer x0_AP(8), y0_AP(8), z0_AP(8)  ! 8 = ダミー, 中心, -x, +x, -y, +y, -z, +z
    integer xlen_AP(8), ylen_AP(8), zlen_AP(8)
    integer size_AP(8)
    !
    ! 通信用
    !
    integer :: grids(6)  ! 通信領域の全グリッド数
    integer :: rank_s (6), rank_r (6)          ! 通信相手のプロセスのランク
    integer :: x0s    (6), y0s    (6), z0s (6) ! 通信領域の原点（s = send, r = receive）
    integer :: x0r    (6), y0r    (6), z0r (6) !
    integer :: xlen   (6), ylen   (6), zlen(6) ! 通信領域の幅
    integer :: xpos   (6), ypos   (6), zpos(6) !

    halox = prm%leng_halo;  haloy = 0;  haloz = 0
    if ( prm%ncomp_grids >= 2 ) haloy = prm%leng_halo
    if ( prm%ncomp_grids == 3 ) haloz = prm%leng_halo
    !
    x0_AP(1) = 0           ;  y0_AP(1) = 0           ;  z0_AP(1) = 0       
    x0_AP(2) = halox       ;  y0_AP(2) = haloy       ;  z0_AP(2) = haloz  
    x0_AP(3) = 0           ;  y0_AP(3) = haloy       ;  z0_AP(3) = haloz  
    x0_AP(4) = prm%nx-halox;  y0_AP(4) = haloy       ;  z0_AP(4) = haloz  
    x0_AP(5) = 0           ;  y0_AP(5) = 0           ;  z0_AP(5) = haloz  
    x0_AP(6) = 0           ;  y0_AP(6) = prm%ny-haloy;  z0_AP(6) = haloz  
    x0_AP(7) = 0           ;  y0_AP(7) = 0           ;  z0_AP(7) = 0      
    x0_AP(8) = 0           ;  y0_AP(8) = 0           ;  z0_AP(8) = prm%nz-haloz  
    !
    xlen_AP(1) = 0             ;  ylen_AP(1) = 0             ;  zlen_AP(1) = 0     
    xlen_AP(2) = prm%nx-2*halox;  ylen_AP(2) = prm%ny-2*haloy;  zlen_AP(2) = prm%nz-2*haloz         
    xlen_AP(3) = halox         ;  ylen_AP(3) = prm%ny-2*haloy;  zlen_AP(3) = prm%nz-2*haloz         
    xlen_AP(4) = halox         ;  ylen_AP(4) = prm%ny-2*haloy;  zlen_AP(4) = prm%nz-2*haloz         
    xlen_AP(5) = prm%nx        ;  ylen_AP(5) = haloy         ;  zlen_AP(5) = prm%nz-2*haloz         
    xlen_AP(6) = prm%nx        ;  ylen_AP(6) = haloy         ;  zlen_AP(6) = prm%nz-2*haloz         
    xlen_AP(7) = prm%nx        ;  ylen_AP(7) = prm%ny        ;  zlen_AP(7) = haloz         
    xlen_AP(8) = prm%nx        ;  ylen_AP(8) = prm%ny        ;  zlen_AP(8) = haloz         
    !
    do i=1,8
       size_AP(i) = xlen_AP(i) * ylen_AP(i) * zlen_AP(i)
    end do


    xpos(1) =-1;  ypos(1) = 0;  zpos(1) = 0  ! -x
    xpos(2) = 1;  ypos(2) = 0;  zpos(2) = 0  ! +x
    xpos(3) = 0;  ypos(3) =-1;  zpos(3) = 0  ! -y
    xpos(4) = 0;  ypos(4) = 1;  zpos(4) = 0  ! +y
    xpos(5) = 0;  ypos(5) = 0;  zpos(5) =-1  ! -z
    xpos(6) = 0;  ypos(6) = 0;  zpos(6) = 1  ! +z
    rank_s(:) = -1
    rank_r(:) = -1
    do i=1, 6
       k = 1;
       kxmax = 0;  kymax = 0;  kzmax = 0;
       if ( prm%npe_z > 1 )  kzmax = 1;
       if ( prm%npe_y > 1 )  kymax = 1; 
       if ( prm%npe_x > 1 )  kxmax = 1; 
       Do kz=-kzmax, kzmax
          Do ky=-kymax, kymax
             Do kx=-kxmax, kxmax
                if ( kx*kx +ky*ky +kz*kz == 0 ) cycle
                if (prm%neighbor_ranks(k) > -1 .and. prm%comm_flag(k) == 1 ) then
                   if ( kx==xpos(i) .and. ky==ypos(i) .and. kz==zpos(i) ) then
                      rank_s(i) = prm%neighbor_ranks(k)
                   end if
                   if ( kx==-xpos(i) .and. ky==-ypos(i) .and. kz==-zpos(i) ) then
                      rank_r(i) = prm%neighbor_ranks(k)
                   end if
                end if
                k = k +1
             End Do
          End Do
       End Do
    end do


    !
    ! 通信領域の原点と形状 *0*, *len, grids の取得
    !
    x0s(1)    = 0
    x0s(2)    = prm%nx - prm%leng_halo
    x0s(3:6)  = -prm%leng_halo
    x0r(1)    = prm%nx
    x0r(2:6)  = -prm%leng_halo
    xlen(1:2) = prm%leng_halo
    xlen(3:6) = prm%nx + 2 * prm%leng_halo
    !
    y0s(1:3)  = 0
    y0s(4)    = prm%ny - prm%leng_halo
    y0s(5:6)  = -prm%leng_halo
    y0r(1:2)  = 0
    y0r(3)    = prm%ny
    y0r(4:6)  = -prm%leng_halo
    ylen(1:2) = prm%ny
    ylen(3:4) = prm%leng_halo
    ylen(5:6) = prm%ny + 2 * prm%leng_halo
    !
    z0s(1:5)  = 0
    z0s(6)    = prm%nz - prm%leng_halo
    z0r(1:4)  = 0
    z0r(5)    = prm%nz
    z0r(6)    = -prm%leng_halo
    zlen(1:4) = prm%nz
    zlen(5:6) = prm%leng_halo
    !
    do i=1, 6
       if ( prm%ncomp_grids == 1 ) then
          grids(i) = xlen(i)
       else if ( prm%ncomp_grids == 2 ) then
          grids(i) = xlen(i) * ylen(i)
       else
          grids(i) = xlen(i) * ylen(i) * zlen(i)
       end if
    end do


    !
    ! 通信処理
    !
    do i=1, 6
       if ( rank_s(i) /= -1 ) allocate( sbuf(i)%v(grids(i)) )
       if ( rank_r(i) /= -1 ) allocate( rbuf(i)%v(grids(i)) )
    end do
    !
    do xyz=0, 3  ! x, y, z, (3のときは通信せず計算のみ)
       !
       !$OMP PARALLEL private(pm, i, ix, iy, iz, k, j, ip)
       !
       do pm=0, 1  ! +xyz, -xyz
          if(xyz == 3) exit
          i = 2 * xyz + pm + 1
          !
          ! 通信データの設定
          !
          if ( rank_s(i) /= -1 ) then
             !
             !$OMP do
             !
             do j=0, grids(i)-1
                ix = x0s(i) + mod(j, xlen(i))
                iy = y0s(i) + mod(j/xlen(i),  ylen(i))
                iz = z0s(i) + j / (xlen(i) * ylen(i))
                k = (ix+prm%stmx) + prm%mx*(iy+prm%stmy) + prm%mxy*(iz+prm%stmz) +1;
                sbuf(i)%v(j + 1) = P(k)
             end do
             !
             !$OMP end do
             !
          end if
       end do
       !
       ! 通信
       !
       !$OMP master
       !
       do pm=0, 1  ! +xyz, -xyz
          if(xyz == 3) exit
          i = 2 * xyz + pm + 1
          if ( rank_s(i) /= -1 ) then
             call MPI_Isend( sbuf(i)%v, size(sbuf(i)%v), mpi_double_precision, &
                  & rank_s(i), rank_s(i), newcomm, ireqs(i), info )
             call MPI_Wait( ireqs(i), istatus, info )
          end if
          if ( rank_r(i) /= -1 ) then
             call MPI_Irecv( rbuf(i)%v, size(rbuf(i)%v), mpi_double_precision, &
                  & rank_r(i), myrank, newcomm, ireqr(i), info )
             call MPI_Wait( ireqr(i), istatus, info )
          end if
       end do
       !
       !$OMP end master

       !
       ! A*P
       !
       if (prm%ncomp_grids == 3) then
          i = 2 * xyz + 1
          !$OMP do schedule(dynamic, 1) 
          !
          do iz=z0_AP(i), z0_AP(i)+zlen_AP(i)-1
             Do ip=1, nnd
                do iy=y0_AP(i), y0_AP(i)+ylen_AP(i)-1
                   do ix=x0_AP(i), x0_AP(i)+xlen_AP(i)-1
                      j = ( ix +prm%stmx ) +prm%mx *( iy +prm%stmy )  +prm%mxy *( iz +prm%stmz ) +1;
                      Q(j) = Q(j) +A( mdist*(ip-1) +j ) *P( j +ioff_A(ip) )
                   end do
                end do
             end do
          end do
          !
          !$OMP end do
          !
          i = 2 * xyz + 2
          !$OMP do schedule(dynamic, 1) 
          !
          do iz=z0_AP(i), z0_AP(i)+zlen_AP(i)-1
             Do ip=1, nnd
                do iy=y0_AP(i), y0_AP(i)+ylen_AP(i)-1
                   do ix=x0_AP(i), x0_AP(i)+xlen_AP(i)-1
                      j = ( ix +prm%stmx ) +prm%mx *( iy +prm%stmy )  +prm%mxy *( iz +prm%stmz ) +1;
                      Q(j) = Q(j) +A( mdist*(ip-1) +j ) *P( j +ioff_A(ip) )
                   end do
                end do
             end do
          end do
          !
          !$OMP end do
       else
          i = 2 * xyz + 1
          !$OMP do schedule(dynamic, 1) collapse(2)
          !
          do iz=z0_AP(i), z0_AP(i)+zlen_AP(i)-1
             do iy=y0_AP(i), y0_AP(i)+ylen_AP(i)-1
                Do ip=1, nnd
                   do ix=x0_AP(i), x0_AP(i)+xlen_AP(i)-1
                      j = ( ix +prm%stmx ) +prm%mx *( iy +prm%stmy )  +prm%mxy *( iz +prm%stmz ) +1;
                      Q(j) = Q(j) +A( mdist*(ip-1) +j ) *P( j +ioff_A(ip) )
                   end do
                end do
             end do
          end do
          !
          !$OMP end do
          !
          i = 2 * xyz + 2
          !$OMP do schedule(dynamic, 1) collapse(2)
          !
          do iz=z0_AP(i), z0_AP(i)+zlen_AP(i)-1
             do iy=y0_AP(i), y0_AP(i)+ylen_AP(i)-1
                Do ip=1, nnd
                   do ix=x0_AP(i), x0_AP(i)+xlen_AP(i)-1
                      j = ( ix +prm%stmx ) +prm%mx *( iy +prm%stmy )  +prm%mxy *( iz +prm%stmz ) +1;
                      Q(j) = Q(j) +A( mdist*(ip-1) +j ) *P( j +ioff_A(ip) )
                   end do
                end do
             end do
          end do
          !
          !$OMP end do
       end if
       !
       ! 受信したデータをコピー
       !
       do pm=0, 1
          if(xyz == 3) exit
          i = 2 * xyz + pm + 1
          !
          ! 受信したデータをコピー
          !
          if ( rank_r(i) /= -1 ) then
             !
             !$OMP do
             !
             do j=0, grids(i)-1
                ix = x0r(i) + mod(j, xlen(i))
                iy = y0r(i) + mod(j/xlen(i),  ylen(i))
                iz = z0r(i) + j / (xlen(i) * ylen(i))
                k = (ix+prm%stmx) + prm%mx*(iy+prm%stmy) + prm%mxy*(iz+prm%stmz) +1;
                P(k) = rbuf(i)%v(j + 1) 
             end do
             !
             !$OMP end do
             !
          end if
       end do
       !
       !$OMP end PARALLEL

    end do

    do i=1, 6
       if (allocated(sbuf(i)%v)) deallocate( sbuf(i)%v )
       if (allocated(rbuf(i)%v)) deallocate( rbuf(i)%v )
    end do

  end subroutine spaxv_ddm_commhide_2


  subroutine init_sleev_comm( &
       ncomp_grids, num_grids, ndiv_grids, &
       num_neighbor_ranks, neighbor_ranks, &
       leng_halo, margin, prm )
    ! use ddm_tmp
    use ddmMat_wrap
    implicit none

    ! inout
    integer, intent(in) :: num_neighbor_ranks
    integer, intent(in) :: neighbor_ranks( num_neighbor_ranks )
    integer, intent(in) :: ncomp_grids, leng_halo
    integer, intent(in) :: ndiv_grids( ncomp_grids )
    integer, intent(in) :: num_grids( ncomp_grids )
    integer, intent(in) :: margin
    type(ddm_param), intent(inout) :: prm

    ! local
    integer i

    ! begin
    prm%stmx = 0;   prm%stmy = 0;   prm%stmz = 0;
    prm%nx = 1;     prm%ny = 1;     prm%nz = 1;
    prm%mx = 1;     prm%my = 1;     prm%mz = 1;
    prm%npe_x = 1;  prm%npe_y = 1;  prm%npe_z = 1;

    if ( ncomp_grids == 1 ) then
       prm%stmx = margin;
       prm%nx   = num_grids(1);
       prm%mx   = prm%nx + 2 *prm%stmx;
       prm%npe_x = ndiv_grids(1);
    endif
    if ( ncomp_grids == 2 ) then
       prm%stmx = margin;       prm%stmy = margin;
       prm%nx   = num_grids(1); prm%ny   = num_grids(2);
       prm%mx   = prm%nx + 2 *prm%stmx;
       prm%my   = prm%ny + 2 *prm%stmy;
       prm%npe_x = ndiv_grids(1);
       prm%npe_y = ndiv_grids(2);
    endif
    if ( ncomp_grids == 3 ) then
       prm%stmx = margin;        prm%stmy = margin;      prm%stmz = margin;
       prm%nx   = num_grids(1);  prm%ny = num_grids(2);  prm%nz = num_grids(3);
       prm%mx   = prm%nx + 2 *prm%stmx;
       prm%my   = prm%ny + 2 *prm%stmy;
       prm%mz   = prm%nz + 2 *prm%stmz;
       prm%npe_x = ndiv_grids(1);
       prm%npe_y = ndiv_grids(2);
       prm%npe_z = ndiv_grids(3);
    endif

    prm%mxy = prm%mx *prm%my;
    prm%nxy = prm%nx *prm%ny;   prm%nxyz = prm%nxy *prm%nz;

    prm%num_neighbor_ranks = num_neighbor_ranks;

    allocate( prm%neighbor_ranks( num_neighbor_ranks ) )
    Do i=1, num_neighbor_ranks
       prm%neighbor_ranks(i) = neighbor_ranks(i);
    End Do
    prm%ncomp_grids = ncomp_grids;
    prm%leng_halo = leng_halo;

  end subroutine init_sleev_comm

  subroutine sleev_comm_novec2( mdist, vec, prm, myrank, comm )
    ! use ddm_tmp
    use ddmMat_wrap
    use mpi
    implicit none
!    include "mpif.h"

    ! inout
    integer(kind=intType), intent(in) :: mdist
    real*8,intent(inout) :: vec( mdist )
    type(ddm_param), intent(in) :: prm
    integer, intent(in) :: comm, myrank

    ! local
    integer(kind=intType) kx, ky, kz, kxmax, kymax, kzmax
    integer(kind=intType) ix, iy, iz, i, j, k
    integer info
    type array
       real*8, allocatable :: v(:)
    end type array
    type(array) sbuf(6), rbuf(6) ! 6 = ±x, ±y, ±z   
    integer :: ireqs(6), ireqr(6)
    integer istatus (mpi_status_size)
    integer xyz, pm
    !
    integer :: grids(6)  ! 通信領域の全グリッド数
    integer :: rank_s (6), rank_r (6)          ! 通信相手のプロセスのランク
    integer :: x0s    (6), y0s    (6), z0s (6) ! 通信領域の原点（s = send, r = receive）
    integer :: x0r    (6), y0r    (6), z0r (6) !
    integer :: xlen   (6), ylen   (6), zlen(6) ! 通信領域の幅
    integer :: xpos   (6), ypos   (6), zpos(6) !

    if ( prm%npe_x * prm%npe_y * prm%npe_z == 1 ) return;


    xpos(1) =-1;  ypos(1) = 0;  zpos(1) = 0  ! -x
    xpos(2) = 1;  ypos(2) = 0;  zpos(2) = 0  ! +x
    xpos(3) = 0;  ypos(3) =-1;  zpos(3) = 0  ! -y
    xpos(4) = 0;  ypos(4) = 1;  zpos(4) = 0  ! +y
    xpos(5) = 0;  ypos(5) = 0;  zpos(5) =-1  ! -z
    xpos(6) = 0;  ypos(6) = 0;  zpos(6) = 1  ! +z
    rank_s(:) = -1
    rank_r(:) = -1
    do i=1, 6
       k = 1;
       kxmax = 0;  kymax = 0;  kzmax = 0;
       if ( prm%npe_z > 1 )  kzmax = 1;
       if ( prm%npe_y > 1 )  kymax = 1; 
       if ( prm%npe_x > 1 )  kxmax = 1; 
       Do kz=-kzmax, kzmax
          Do ky=-kymax, kymax
             Do kx=-kxmax, kxmax
                if ( kx*kx +ky*ky +kz*kz == 0 ) cycle
                if (prm%neighbor_ranks(k) > -1 .and. prm%comm_flag(k) == 1 ) then
                   if ( kx==xpos(i) .and. ky==ypos(i) .and. kz==zpos(i) ) then
                      rank_s(i) = prm%neighbor_ranks(k)
                   end if
                   if ( kx==-xpos(i) .and. ky==-ypos(i) .and. kz==-zpos(i) ) then
                      rank_r(i) = prm%neighbor_ranks(k)
                   end if
                end if
                k = k +1
             End Do
          End Do
       End Do
    end do

    !
    ! 通信領域の原点と形状 *0*, *len, grids の取得
    !
    x0s(1)    = 0
    x0s(2)    = prm%nx - prm%leng_halo
    x0s(3:6)  = -prm%leng_halo
    x0r(1)    = prm%nx
    x0r(2:6)  = -prm%leng_halo
    xlen(1:2) = prm%leng_halo
    xlen(3:6) = prm%nx + 2 * prm%leng_halo
    !
    y0s(1:3)  = 0
    y0s(4)    = prm%ny - prm%leng_halo
    y0s(5:6)  = -prm%leng_halo
    y0r(1:2)  = 0
    y0r(3)    = prm%ny
    y0r(4:6)  = -prm%leng_halo
    ylen(1:2) = prm%ny
    ylen(3:4) = prm%leng_halo
    ylen(5:6) = prm%ny + 2 * prm%leng_halo
    !
    z0s(1:5)  = 0
    z0s(6)    = prm%nz - prm%leng_halo
    z0r(1:4)  = 0
    z0r(5)    = prm%nz
    z0r(6)    = -prm%leng_halo
    zlen(1:4) = prm%nz
    zlen(5:6) = prm%leng_halo
    !
    do i=1, 6
       if ( prm%ncomp_grids == 1 ) then
          grids(i) = xlen(i)
       else if ( prm%ncomp_grids == 2 ) then
          grids(i) = xlen(i) * ylen(i)
       else
          grids(i) = xlen(i) * ylen(i) * zlen(i)
       end if
    end do

    ! 通信処理
    do i=1, 6
       if ( rank_s(i) /= -1 ) allocate( sbuf(i)%v(grids(i)) )
       if ( rank_r(i) /= -1 ) allocate( rbuf(i)%v(grids(i)) )
    end do

    do xyz=0, 2  ! x, y, z
       ! 同じ軸方向の領域のみ同時に通信できる
       do pm=0, 1  ! +xyz, -xyz
          i = 2 * xyz + pm + 1
          ! 通信データの設定
          if ( rank_s(i) /= -1 ) then
             !$OMP PARALLEL do private(ix, iy, iz, k)
             do j=0, grids(i)-1
                ix = x0s(i) + mod(j, xlen(i))
                iy = y0s(i) + mod(j/xlen(i),  ylen(i))
                iz = z0s(i) + j / (xlen(i) * ylen(i))
                k = (ix+prm%stmx) + prm%mx*(iy+prm%stmy) + prm%mxy*(iz+prm%stmz) +1;
                sbuf(i)%v(j + 1) = vec(k)
             end do
             !$OMP end PARALLEL do
          end if
          !
          ! 通信開始
          !
          if ( rank_s(i) /= -1 ) then
             call MPI_Isend( sbuf(i)%v, size(sbuf(i)%v), mpi_double_precision, &
                  rank_s(i), rank_s(i), comm, ireqs(i), info )
          end if
          if ( rank_r(i) /= -1 ) then
             call MPI_Irecv( rbuf(i)%v, size(rbuf(i)%v), mpi_double_precision, &
                  rank_r(i), myrank, comm, ireqr(i), info )
          end if
       end do
       !
       ! 通信を完了させて、受信したデータをコピー
       !
       do pm=0, 1
          i = 2 * xyz + pm + 1
          if ( rank_s(i) /= -1 ) call MPI_Wait( ireqs(i), istatus, info )
          if ( rank_r(i) /= -1 ) call MPI_Wait( ireqr(i), istatus, info )
          !
          ! 受信したデータをコピー
          !
          if ( rank_r(i) /= -1 ) then
             !$OMP PARALLEL do private(ix, iy, iz, k)
             do j=0, grids(i)-1
                ix = x0r(i) + mod(j, xlen(i))
                iy = y0r(i) + mod(j/xlen(i),  ylen(i))
                iz = z0r(i) + j / (xlen(i) * ylen(i))
                k = (ix+prm%stmx) + prm%mx*(iy+prm%stmy) + prm%mxy*(iz+prm%stmz) +1;
                vec(k) = rbuf(i)%v(j + 1) 
             end do
             !$OMP end PARALLEL do
          end if
       end do
    end do

    do i=1, 6
       if (allocated(sbuf(i)%v)) deallocate( sbuf(i)%v )
       if (allocated(rbuf(i)%v)) deallocate( rbuf(i)%v )
    end do

  end subroutine sleev_comm_novec2

  subroutine sleev_comm_GPU( &
       method, &
       vec,offset_vec)
    ! use ddm_tmp
    use krylov,only : KSP
    use ddmMat_wrap
    use mpi
    use PARCEL_VEC_kernel,only : PARCEL_VEC
    use cu_comm,only : &
         C_MPI_Wait,C_MPI_Irecv_FP64,C_MPI_Isend_FP64,&
         set_sbuf_DDM_FP64_GPU,set_rbuf_DDM_FP64_GPU, &
         C_MPI_Isend_FP64buf_FP64,C_MPI_Irecv_FP64buf_FP64, &
         comm_ddm_GPU,comm_ddm_GPU_host
    use cu_common_wrapper,only : &
         cudaDeviceSynchronize, &
         cudaMemcpy_d2h_fp64,cudaMemcpy_h2d_fp64

    implicit none

    ! inout
    type(KSP) :: method
    type(devicemem) :: vec
    integer(kind=intType) offset_vec
    integer(kind=intType),parameter :: offset_zero=0

    ! local
    integer myrank
    integer(kind=intType) kx, ky, kz, kxmax, kymax, kzmax
    integer(kind=intType) ix, iy, iz, i, j, k
    integer info

    integer :: ireqs(6), ireqr(6)
    integer istatus (mpi_status_size)
    integer(kind=intType) xyz, pm
    !
    integer :: grids(6)  ! 通信領域の全グリッド数
    integer :: rank_s (6), rank_r (6)          ! 通信相手のプロセスのランク
    integer :: x0s    (6), y0s    (6), z0s (6) ! 通信領域の原点（s = send, r = receive）
    integer :: x0r    (6), y0r    (6), z0r (6) !
    integer :: xlen   (6), ylen   (6), zlen(6) ! 通信領域の幅
    integer :: xpos   (6), ypos   (6), zpos(6) !

!    real*8 vec_cpu(method%ddmMat%n)
!    real*8 sbuf(method%ddmMat%prm%nx*method%ddmMat%prm%ny*method%ddmMat%prm%nz,6)
!    real*8 rbuf(method%ddmMat%prm%nx*method%ddmMat%prm%ny*method%ddmMat%prm%nz,6)
    integer offset 

    if ( method%ddmMat%prm%npe_x * method%ddmMat%prm%npe_y * method%ddmMat%prm%npe_z == 1 ) return
    myrank = method%commVec%myrank
    xpos(1) =-1;  ypos(1) = 0;  zpos(1) = 0  ! -x
    xpos(2) = 1;  ypos(2) = 0;  zpos(2) = 0  ! +x
    xpos(3) = 0;  ypos(3) =-1;  zpos(3) = 0  ! -y
    xpos(4) = 0;  ypos(4) = 1;  zpos(4) = 0  ! +y
    xpos(5) = 0;  ypos(5) = 0;  zpos(5) =-1  ! -z
    xpos(6) = 0;  ypos(6) = 0;  zpos(6) = 1  ! +z
    rank_s(:) = -1
    rank_r(:) = -1
    do i=1, 6
       k = 1;
       kxmax = 0;  kymax = 0;  kzmax = 0;
       if ( method%ddmMat%prm%npe_z > 1 )  kzmax = 1;
       if ( method%ddmMat%prm%npe_y > 1 )  kymax = 1; 
       if ( method%ddmMat%prm%npe_x > 1 )  kxmax = 1; 
       Do kz=-kzmax, kzmax
          Do ky=-kymax, kymax
             Do kx=-kxmax, kxmax
                if ( kx*kx +ky*ky +kz*kz == 0 ) cycle
                if (method%ddmMat%prm%neighbor_ranks(k) > -1 .and. method%ddmMat%prm%comm_flag(k) == 1 ) then
                   if ( kx==xpos(i) .and. ky==ypos(i) .and. kz==zpos(i) ) then
                      rank_s(i) = method%ddmMat%prm%neighbor_ranks(k)
                   end if
                   if ( kx==-xpos(i) .and. ky==-ypos(i) .and. kz==-zpos(i) ) then
                      rank_r(i) = method%ddmMat%prm%neighbor_ranks(k)
                   end if
                end if
                k = k +1
             End Do
          End Do
       End Do
    end do

    !
    ! 通信領域の原点と形状 *0*, *len, grids の取得
    !
    x0s(1)    = 0
    x0s(2)    = method%ddmMat%prm%nx - method%ddmMat%prm%leng_halo
    x0s(3:6)  = -method%ddmMat%prm%leng_halo
    x0r(1)    = method%ddmMat%prm%nx
    x0r(2:6)  = -method%ddmMat%prm%leng_halo
    xlen(1:2) = method%ddmMat%prm%leng_halo
    xlen(3:6) = method%ddmMat%prm%nx + 2 * method%ddmMat%prm%leng_halo
    !
    y0s(1:3)  = 0
    y0s(4)    = method%ddmMat%prm%ny - method%ddmMat%prm%leng_halo
    y0s(5:6)  = -method%ddmMat%prm%leng_halo
    y0r(1:2)  = 0
    y0r(3)    = method%ddmMat%prm%ny
    y0r(4:6)  = -method%ddmMat%prm%leng_halo
    ylen(1:2) = method%ddmMat%prm%ny
    ylen(3:4) = method%ddmMat%prm%leng_halo
    ylen(5:6) = method%ddmMat%prm%ny + 2 * method%ddmMat%prm%leng_halo
    !
    z0s(1:5)  = 0
    z0s(6)    = method%ddmMat%prm%nz - method%ddmMat%prm%leng_halo
    z0r(1:4)  = 0
    z0r(5)    = method%ddmMat%prm%nz
    z0r(6)    = -method%ddmMat%prm%leng_halo
    zlen(1:4) = method%ddmMat%prm%nz
    zlen(5:6) = method%ddmMat%prm%leng_halo
    !
    do i=1, 6
       if ( method%ddmMat%prm%ncomp_grids == 1 ) then
          grids(i) = xlen(i)
       else if ( method%ddmMat%prm%ncomp_grids == 2 ) then
          grids(i) = xlen(i) * ylen(i)
       else
          grids(i) = xlen(i) * ylen(i) * zlen(i)
       end if
    end do

    call cudaDeviceSynchronize()     

    ! 通信処理
    offset = 1
    do xyz=0, 2  ! x, y, z
       do pm=0, 1  ! +xyz, -xyz
          i = 2 * xyz + pm + 1
          if ( rank_s(i) /= -1 ) then
             call set_sbuf_DDM_FP64_GPU( &
                  CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
                  vec,offset_vec, &
                  method%ddmMat%mat%C_sbuf, offset_zero, &
                  offset, &
                  grids(i), &
                  method%ddmMat%prm%mx,method%ddmMat%prm%mxy, &
                  x0s(i),y0s(i),z0s(i), &
                  xlen(i),ylen(i), &
                  method%ddmMat%prm%stmx, &
                  method%ddmMat%prm%stmy, &
                  method%ddmMat%prm%stmz &
                  )
          end if
          offset = offset + grids(i)
       enddo
    enddo
          !
          ! 通信開始
          !
#if CUDA_AWARE
    call cudaDeviceSynchronize()     
    call comm_ddm_GPU( &
         method%ddmMat%prm%nx*method%ddmMat%prm%ny*method%ddmMat%prm%nz, &
         method%ddmMat%mat%C_sbuf, offset_zero,&
         method%ddmMat%mat%C_rbuf, offset_zero,&
         grids, &
         myrank, &
         rank_s,rank_r, &
         method%commVec%C_comm, &
         method%commVec%C_s_sta,method%commVec%C_r_sta, &
         method%commVec%C_s_req,method%commVec%C_r_req &
         )
    call cudaDeviceSynchronize()     
#else
    call cudaDeviceSynchronize()     
    call comm_ddm_GPU_host( &
         method%ddmMat%prm%nx*method%ddmMat%prm%ny*method%ddmMat%prm%nz, &
         method%ddmMat%mat%C_sbuf, offset_zero,&
         method%ddmMat%mat%C_rbuf, offset_zero,&
         grids, &
         myrank, &
         rank_s,rank_r, &
         method%commVec%C_comm, &
         method%commVec%C_s_sta,method%commVec%C_r_sta, &
         method%commVec%C_s_req,method%commVec%C_r_req &
         )
    call cudaDeviceSynchronize()     
#endif

    offset = 1
    do xyz=0, 2  ! x, y, z
       do pm=0, 1
          i = 2 * xyz + pm + 1
          !
          ! 受信したデータをコピー
          !
          if ( rank_r(i) /= -1 ) then
             call set_rbuf_DDM_FP64_GPU( &
                  CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
                  vec, offset_vec,&
                  method%ddmMat%mat%C_rbuf, offset_zero,&
                  offset, &
                  grids(i), &
                  method%ddmMat%prm%mx,method%ddmMat%prm%mxy, &
                  x0r(i),y0r(i),z0r(i), &
                  xlen(i),ylen(i), &
                  method%ddmMat%prm%stmx, &
                  method%ddmMat%prm%stmy, &
                  method%ddmMat%prm%stmz &
                  )
          end if
          offset = offset + grids(i)
       end do
    end do

    call cudaDeviceSynchronize()     

  end subroutine sleev_comm_GPU

  subroutine sleev_zero_GPU( &
       method,&
       vec,offset_vec)
    ! use ddm_tmp
    use krylov,only : KSP
    use ddmMat_wrap
    use mpi
    use PARCEL_VEC_kernel,only : PARCEL_VEC
    use cu_comm,only : sleev_zero_DDM_FP64_GPU
    use cu_common_wrapper,only : &
         cudaDeviceSynchronize

    implicit none

    ! inout
    type(KSP) :: method
    type(devicemem) :: vec
    integer(kind=intType) offset_vec
    integer(kind=intType),parameter :: offset_zero=0

    ! local
    integer myrank
    integer(kind=intType) kx, ky, kz, kxmax, kymax, kzmax
    integer(kind=intType) ix, iy, iz, i, j, k
    integer info
    integer(kind=intType) xyz, pm
    !
    integer :: grids(6)  ! 通信領域の全グリッド数

    integer :: rank_r (6)          ! 通信相手のプロセスのランク
    integer :: x0r    (6), y0r    (6), z0r (6) !
    integer :: xlen   (6), ylen   (6), zlen(6) ! 通信領域の幅
    integer :: xpos   (6), ypos   (6), zpos(6) !

    integer offset 
    if ( method%ddmMat%prm%npe_x * method%ddmMat%prm%npe_y * method%ddmMat%prm%npe_z == 1 ) return
    myrank = method%commVec%myrank
    xpos(1) =-1;  ypos(1) = 0;  zpos(1) = 0  ! -x
    xpos(2) = 1;  ypos(2) = 0;  zpos(2) = 0  ! +x
    xpos(3) = 0;  ypos(3) =-1;  zpos(3) = 0  ! -y
    xpos(4) = 0;  ypos(4) = 1;  zpos(4) = 0  ! +y
    xpos(5) = 0;  ypos(5) = 0;  zpos(5) =-1  ! -z
    xpos(6) = 0;  ypos(6) = 0;  zpos(6) = 1  ! +z

    rank_r(:) = -1
    do i=1, 6
       k = 1;
       kxmax = 0;  kymax = 0;  kzmax = 0;
       if ( method%ddmMat%prm%npe_z > 1 )  kzmax = 1;
       if ( method%ddmMat%prm%npe_y > 1 )  kymax = 1; 
       if ( method%ddmMat%prm%npe_x > 1 )  kxmax = 1; 
       do kz=-kzmax, kzmax
          do ky=-kymax, kymax
             do kx=-kxmax, kxmax
                if ( kx*kx +ky*ky +kz*kz == 0 ) cycle
                if (method%ddmMat%prm%neighbor_ranks(k) > -1 .and. method%ddmMat%prm%comm_flag(k) == 1 ) then
                   if ( kx==-xpos(i) .and. ky==-ypos(i) .and. kz==-zpos(i) ) then
                      rank_r(i) = method%ddmMat%prm%neighbor_ranks(k)
                   end if
                end if
                k = k +1
             enddo
          enddo
       enddo
    enddo

    !
    ! 通信領域の原点と形状 *0*, *len, grids の取得
    !
    x0r(1)    = method%ddmMat%prm%nx
    x0r(2:6)  = -method%ddmMat%prm%leng_halo
    xlen(1:2) = method%ddmMat%prm%leng_halo
    xlen(3:6) = method%ddmMat%prm%nx + 2 * method%ddmMat%prm%leng_halo
    !
    y0r(1:2)  = 0
    y0r(3)    = method%ddmMat%prm%ny
    y0r(4:6)  = -method%ddmMat%prm%leng_halo
    ylen(1:2) = method%ddmMat%prm%ny
    ylen(3:4) = method%ddmMat%prm%leng_halo
    ylen(5:6) = method%ddmMat%prm%ny + 2 * method%ddmMat%prm%leng_halo
    !
    z0r(1:4)  = 0
    z0r(5)    = method%ddmMat%prm%nz
    z0r(6)    = -method%ddmMat%prm%leng_halo
    zlen(1:4) = method%ddmMat%prm%nz
    zlen(5:6) = method%ddmMat%prm%leng_halo
    !
    do i=1, 6
       if ( method%ddmMat%prm%ncomp_grids == 1 ) then
          grids(i) = xlen(i)
       else if ( method%ddmMat%prm%ncomp_grids == 2 ) then
          grids(i) = xlen(i) * ylen(i)
       else
          grids(i) = xlen(i) * ylen(i) * zlen(i)
       end if
    end do

    call cudaDeviceSynchronize()     

    offset = 1
    do xyz=0, 2  ! x, y, z
       do pm=0, 1
          i = 2 * xyz + pm + 1
          if ( rank_r(i) /= -1 ) then
             ! 受信データを詰めた場所のみゼロにする
             call sleev_zero_DDM_FP64_GPU( &
                  CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
                  vec, offset_vec, &
                  offset, &
                  grids(i), &
                  method%ddmMat%prm%mx,method%ddmMat%prm%mxy, &
                  x0r(i),y0r(i),z0r(i), &
                  xlen(i),ylen(i), &
                  method%ddmMat%prm%stmx, &
                  method%ddmMat%prm%stmy, &
                  method%ddmMat%prm%stmz &
                  )
          endif
          offset = offset + grids(i)
       end do
    end do

    call cudaDeviceSynchronize()     

  end subroutine sleev_zero_GPU

  subroutine sleev_comm_init_GPU(method)
    ! use ddm_tmp
    use krylov,only : KSP
    use ddmMat_wrap
    use cu_common_wrapper,only : cudaMalloc_fp64,cudaMemcpy_h2d_fp64

    use mpi
    implicit none

    ! inout
    type(KSP) :: method

    ! local
    integer kx, ky, kz, kxmax, kymax, kzmax
    integer ix, iy, iz, i, j, info, k

    integer :: ireqs(6), ireqr(6)
    integer istatus (mpi_status_size)
    integer xyz, pm
    !
    integer :: grids(6)  ! 通信領域の全グリッド数
    integer :: rank_s (6), rank_r (6)          ! 通信相手のプロセスのランク
    integer :: x0s    (6), y0s    (6), z0s (6) ! 通信領域の原点（s = send, r = receive）
    integer :: x0r    (6), y0r    (6), z0r (6) !
    integer :: xlen   (6), ylen   (6), zlen(6) ! 通信領域の幅
    integer :: xpos   (6), ypos   (6), zpos(6) !
    integer grids_max
    real*8, allocatable :: vec(:)

    integer(kind=intType) parcel_int

    if ( method%ddmMat%prm%npe_x * method%ddmMat%prm%npe_y * method%ddmMat%prm%npe_z == 1 ) return;


    xpos(1) =-1;  ypos(1) = 0;  zpos(1) = 0  ! -x
    xpos(2) = 1;  ypos(2) = 0;  zpos(2) = 0  ! +x
    xpos(3) = 0;  ypos(3) =-1;  zpos(3) = 0  ! -y
    xpos(4) = 0;  ypos(4) = 1;  zpos(4) = 0  ! +y
    xpos(5) = 0;  ypos(5) = 0;  zpos(5) =-1  ! -z
    xpos(6) = 0;  ypos(6) = 0;  zpos(6) = 1  ! +z
    rank_s(:) = 1
    rank_r(:) = 1
    do i=1, 6
       k = 1;
       kxmax = 0;  kymax = 0;  kzmax = 0;
       if ( method%ddmMat%prm%npe_z > 1 )  kzmax = 1;
       if ( method%ddmMat%prm%npe_y > 1 )  kymax = 1; 
       if ( method%ddmMat%prm%npe_x > 1 )  kxmax = 1; 
       Do kz=-kzmax, kzmax
          Do ky=-kymax, kymax
             Do kx=-kxmax, kxmax
                if ( kx*kx +ky*ky +kz*kz == 0 ) cycle
                if (method%ddmMat%prm%neighbor_ranks(k) > -1 .and. method%ddmMat%prm%comm_flag(k) == 1 ) then
                   if ( kx==xpos(i) .and. ky==ypos(i) .and. kz==zpos(i) ) then
                      rank_s(i) = method%ddmMat%prm%neighbor_ranks(k)
                   end if
                   if ( kx==-xpos(i) .and. ky==-ypos(i) .and. kz==-zpos(i) ) then
                      rank_r(i) = method%ddmMat%prm%neighbor_ranks(k)
                   end if
                end if
                k = k +1
             End Do
          End Do
       End Do
    end do

    ! 通信領域の原点と形状 *0*, *len, grids の取得
    x0s(1)    = 0
    x0s(2)    = method%ddmMat%prm%nx - method%ddmMat%prm%leng_halo
    x0s(3:6)  = -method%ddmMat%prm%leng_halo
    x0r(1)    = method%ddmMat%prm%nx
    x0r(2:6)  = -method%ddmMat%prm%leng_halo
    xlen(1:2) = method%ddmMat%prm%leng_halo
    xlen(3:6) = method%ddmMat%prm%nx + 2 * method%ddmMat%prm%leng_halo

    y0s(1:3)  = 0
    y0s(4)    = method%ddmMat%prm%ny - method%ddmMat%prm%leng_halo
    y0s(5:6)  = -method%ddmMat%prm%leng_halo
    y0r(1:2)  = 0
    y0r(3)    = method%ddmMat%prm%ny
    y0r(4:6)  = -method%ddmMat%prm%leng_halo
    ylen(1:2) = method%ddmMat%prm%ny
    ylen(3:4) = method%ddmMat%prm%leng_halo
    ylen(5:6) = method%ddmMat%prm%ny + 2 * method%ddmMat%prm%leng_halo

    z0s(1:5)  = 0
    z0s(6)    = method%ddmMat%prm%nz - method%ddmMat%prm%leng_halo
    z0r(1:4)  = 0
    z0r(5)    = method%ddmMat%prm%nz
    z0r(6)    = -method%ddmMat%prm%leng_halo
    zlen(1:4) = method%ddmMat%prm%nz
    zlen(5:6) = method%ddmMat%prm%leng_halo

    do i=1, 6
       if ( method%ddmMat%prm%ncomp_grids == 1 ) then
          grids(i) = xlen(i)
       else if ( method%ddmMat%prm%ncomp_grids == 2 ) then
          grids(i) = xlen(i) * ylen(i)
       else
          grids(i) = xlen(i) * ylen(i) * zlen(i)
       end if
    end do

    grids_max = 1
    do i=1, 6
       grids_max = grids_max + grids(i)      
    enddo
    allocate(vec(grids_max))
    vec = 0.0d0

    parcel_int = grids_max
    call cudaMalloc_fp64(method%ddmMat%mat%C_sbuf,parcel_int)
    call cudaMalloc_fp64(method%ddmMat%mat%C_rbuf,parcel_int)
#if 0
    call cudaMemcpy_h2d_fp64( &
         vec, &
         method%ddmMat%mat%C_sbuf,parcel_int,0,0)
    call cudaMemcpy_h2d_fp64( &
         vec, &
         method%ddmMat%mat%C_rbuf,parcel_int,0,0)
#endif

    deallocate(vec)
    
  end subroutine sleev_comm_init_GPU
  
  subroutine set_sleev_comm_flag( nnd, ncomp_grids, ioff_grids, prm )
    ! use ddm_tmp
    use ddmMat_wrap
    implicit none
    ! inout
    integer, intent(in) :: nnd 
    integer, intent(in) :: ncomp_grids
    integer, intent(in) :: ioff_grids(nnd,ncomp_grids)
    type(ddm_param), intent(inout) :: prm

    ! local
    integer :: num_neighbor_ranks
    integer :: kx, ky, kz, kxmax, kymax, kzmax
    integer :: i, num, npe_x, npe_y, npe_z, found
    integer, allocatable :: ioff_x(:), ioff_y(:), ioff_z(:)

    ! -------------------------------------------
    ! init

    npe_x = prm%npe_x;  npe_y = prm%npe_y;  npe_z = prm%npe_z;
    num_neighbor_ranks = prm%num_neighbor_ranks
    allocate( prm%comm_flag( num_neighbor_ranks ) )

    allocate( ioff_x(nnd) );   allocate( ioff_y(nnd) );   allocate( ioff_z(nnd) ); 
    ioff_x = 0;  ioff_y = 0;  ioff_z = 0

    if ( ncomp_grids == 1 ) then
       Do i=1, nnd
          ioff_x(i) = ioff_grids(i,1)
       End Do
    endif
    if ( ncomp_grids == 2 ) then
       Do i=1, nnd
          ioff_x(i) = ioff_grids(i,1)
          ioff_y(i) = ioff_grids(i,2)
       End Do
    endif
    if ( ncomp_grids == 3 ) then
       Do i=1, nnd
          ioff_x(i) = ioff_grids(i,1)
          ioff_y(i) = ioff_grids(i,2)
          ioff_z(i) = ioff_grids(i,3)
       End Do
    endif

    kxmax = 0;  kymax = 0;  kzmax = 0;
    if ( npe_z > 1 )  kzmax = 1;
    if ( npe_y > 1 )  kymax = 1; 
    if ( npe_x > 1 )  kxmax = 1; 

    prm%comm_flag = 0

    ! ---------------
    num = 1;
    Do kz=-kzmax, kzmax
       Do ky=-kymax, kymax
          Do kx=-kxmax, kxmax
             if ( kx*kx +ky*ky +kz*kz == 0 ) cycle

             found = 0
             Loop1: Do i=1, nnd
                if ( kx == ioff_x(i) .and. ky == ioff_y(i) .and. kz == ioff_z(i) ) then
                   found = 1
                   exit Loop1
                endif
             End Do Loop1
             if ( found == 1 ) prm%comm_flag(num) = 1
             num = num +1
          End do
       End do
    End Do
    !  comm_flag =1
    deallocate( ioff_x );  deallocate( ioff_y );  deallocate( ioff_z )

  end subroutine set_sleev_comm_flag

  subroutine set_div_dir_threads( prm, div_direc_in )
    ! use ddm_tmp
    use ddmMat_wrap
    use omp_lib
    implicit none
    ! include 'omp_lib.h'
    !
    type(ddm_param), intent(inout) :: prm
    integer, intent(in) :: div_direc_in
    !
    integer :: num_threads
    logical :: xdiv, ydiv, zdiv
    integer :: div_direc

    div_direc = div_direc_in

    !$OMP Parallel
    num_threads = omp_get_num_threads();
    !$OMP end Parallel

    ! if ( num_threads == 1 ) div_direc = 0

    if ( prm%ncomp_grids == 1 ) then
       if ( div_direc >=2 ) stop "div_direc should be lower than 2"
    endif
    if ( prm%ncomp_grids == 2 ) then
       if ( div_direc >=3 ) stop "div_direc should be lower than 3"
    endif
    
    if ( (div_direc >= 1) .and. (div_direc <=3) ) then
       prm%div_direc_th = div_direc
       return
    endif

    ! --- automatic setting  ---
    xdiv = .false.;   ydiv = .false.;   zdiv = .false.
    if ( prm%ncomp_grids == 1 ) then
       xdiv = .true.
    endif
    if ( prm%ncomp_grids == 2 ) then
       if ( prm%ny >= prm%nx ) then
          ydiv = .true.
       else
          xdiv = .true.
       endif
    endif
    if ( prm%ncomp_grids == 3 ) then
       if ( prm%nz >= prm%nx .and. prm%nz >= prm%ny ) then
          zdiv = .true.
       else if ( prm%ny >= prm%nx .and. prm%ny >= prm%nz ) then
          ydiv = .true.
       else if ( prm%nx >= prm%ny .and. prm%nx >= prm%nz ) then
          xdiv = .true.
       endif
    endif

    if ( xdiv ) div_direc = 1
    if ( ydiv ) div_direc = 2
    if ( zdiv ) div_direc = 3

    prm%div_direc_th = div_direc
    return

  end subroutine set_div_dir_threads

  subroutine set_region_for_threads2( &
       my_thread, num_threads, nxmin, nxmax, nymin, nymax, nzmin, nzmax,&
       nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th, div_dir )
    ! use ddm_tmp
    use ddmMat_wrap
    implicit none
    ! inout
    integer, intent(in) :: my_thread, num_threads
    integer, intent(in) :: nxmin, nxmax, nymin, nymax, nzmin, nzmax, div_dir
    !
    integer, intent(inout) :: nxs_th, nys_th, nzs_th
    integer, intent(inout) :: nxe_th, nye_th, nze_th
    ! local
    integer n1, nxw, nyw, nzw
    !
    nxs_th = nxmin;  nys_th = nymin;  nzs_th = nzmin;
    nxe_th = nxmax;  nye_th = nymax;  nze_th = nzmax;

    nxw = nxmax - nxmin
    nyw = nymax - nymin
    nzw = nzmax - nzmin

    if ( div_dir == 1 ) then
       n1 = nxw /num_threads;
       nxs_th = nxmin + n1 *my_thread;     nxe_th = nxs_th +n1;
       if ( my_thread == num_threads -1 ) nxe_th = nxmax
    else if ( div_dir == 2 ) then
       n1 = nyw /num_threads;
       nys_th = nymin + n1 *my_thread;     nye_th = nys_th +n1;
       if ( my_thread == num_threads -1 ) nye_th = nymax
    else if ( div_dir == 3 ) then
       n1 = nzw /num_threads;
       nzs_th = nzmin + n1 *my_thread;     nze_th = nzs_th +n1;
       if ( my_thread == num_threads -1 ) nze_th = nzmax
    endif
    if ( div_dir == 0 ) then
       if ( my_thread > 0 ) then
          nxs_th = 0;    nys_th = 0;  nzs_th = 0
          nxe_th = 0;    nye_th = 0;  nze_th = 0
       endif
    endif

    return
  end subroutine set_region_for_threads2

end module ddm_matv
