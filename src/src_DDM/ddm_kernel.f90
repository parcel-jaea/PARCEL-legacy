#define PRINT_TIME 0

#include "ddfunc.h"
module ddm_kernel
  integer,parameter :: DDM_CUDA_BJ = 21
  integer,parameter :: DDM_CUDA_AS = 22
contains

  subroutine wrap_free_ddm(DDM_MAT_)
    use DDM_TYPE
    implicit none
    type (DDM_ALL_TYPE) DDM_MAT_

    if(allocated(DDM_MAT_%AS_v1))then
       deallocate(DDM_MAT_%AS_v1)
    endif

    if(allocated(DDM_MAT_%AS_v2))then
       deallocate(DDM_MAT_%AS_v2)
    endif

    if(allocated(DDM_MAT_%AS_v1_lo))then
       deallocate(DDM_MAT_%AS_v1_lo)
    endif

    if(allocated(DDM_MAT_%AS_v2_lo))then
       deallocate(DDM_MAT_%AS_v2_lo)
    endif

    if(allocated(DDM_MAT_%DDMMat_Neumann%ddmA))then
       deallocate(DDM_MAT_%DDMMat_Neumann%ddmA)
    endif
    if(allocated(DDM_MAT_%DDMMat_Neumann%offset))then
       deallocate(DDM_MAT_%DDMMat_Neumann%offset)
    endif

  end subroutine wrap_free_ddm

  subroutine wrap_init_ddm(&
       DDM_MAT_, &
       precon, &
       ndim, &
       n, &
       n0, &
       nnd, &
       mat, &
       ! b, &
       ! x, &
       ncpu, &
       myrank, &
       icomm, &
       diag, &
       oloop, &
       tol, &
       ipc, &
       nres, &
       res, &
       iret, &
       commHide, &
       AS_len, &
       AS_mean, &
       AS_n, &
       AS_g, &
       AS_val_dia, &
       AS_ioff_dia, &
       AS_prm ,&
       diag_lo, &
       AS_val_dia_lo &
       )  
    use DDM_TYPE
    use ddmMat_wrap
    use param_AS_type_ddm
    use PARCEL_integer
    implicit none

    type(precon_AS_wrap),target :: precon
    integer ndim !dimension of global matrix
    integer(kind=intType) n !number of rows of local matrix  
    integer n0
    integer nnd
    type(ddmmatrix_wrap),target :: mat
    ! double precision,target :: b(n)
    ! double precision,target :: x(n) 
    integer ncpu
    integer myrank
    integer icomm !MPI communicator
    double precision,target :: diag(n)
    double precision,optional,target :: diag_lo(n)
    integer,target :: oloop !outer loop
    double precision tol !tolerance
    integer ipc !precondition flag
    integer nres !number of residual history
    double precision,target :: res(nres) !residual history
    integer iret !result status
    integer commHide         ! 0: no hiding, 1: one time, 2: each time
    integer AS_len           ! additive Schwarz overlap length
    integer AS_mean          ! 1 => 前処理の最後にプロセス間平均をとる
    integer(kind=intType) AS_n
    type(AS_geom),target :: AS_g
    double precision,target :: AS_val_dia(AS_n*nnd)
    double precision,optional,target :: AS_val_dia_lo(AS_n*nnd)
    integer,target :: AS_ioff_dia(nnd)
    type(ddm_param),target :: AS_prm

    ! ---

    type (DDM_ALL_TYPE) DDM_MAT_

    DDM_MAT_%precon=>precon
    DDM_MAT_%ndim=ndim
    DDM_MAT_%n=n
    DDM_MAT_%n0=n0
    DDM_MAT_%nnd=nnd
    DDM_MAT_%istart=-1

    DDM_MAT_%ncpu=ncpu
    DDM_MAT_%myrank=myrank
    DDM_MAT_%icomm=icomm
    DDM_MAT_%diag=>diag
    DDM_MAT_%oloop=>oloop
    DDM_MAT_%tol=tol
    DDM_MAT_%ipc=ipc
    DDM_MAT_%nres=nres
    DDM_MAT_%res=>res
    DDM_MAT_%iret=iret
    DDM_MAT_%commHide=commHide
    DDM_MAT_%AS_len=AS_len
    DDM_MAT_%AS_mean=AS_mean
    DDM_MAT_%AS_n=AS_n
    DDM_MAT_%AS_val_dia=>AS_val_dia
    DDM_MAT_%AS_ioff_dia=>AS_ioff_dia

    if(  present(diag_lo)  )then
       DDM_MAT_%diag_lo=>diag_lo
    endif
    if(  present(AS_val_dia_lo)  )then
       DDM_MAT_%AS_val_dia_lo=>AS_val_dia_lo
    endif

    allocate(DDM_MAT_%AS_v1(AS_n))
    allocate(DDM_MAT_%AS_v2(AS_n))
    allocate(DDM_MAT_%AS_v1_lo(AS_n))
    allocate(DDM_MAT_%AS_v2_lo(AS_n))

  end subroutine wrap_init_ddm

  subroutine wrap_init_ddm_GPU(method)
    use cu_common_wrapper
    use krylov  
    use ddm_matv,only : sleev_comm_init_GPU
    use PARCEL_VEC_kernel,only : parcel_sync
    use PARCEL_integer
    implicit none
    type(KSP) :: method
    integer myrank
    integer (kind=intType) n
    integer (kind=intType) parcel_int

    myrank=method%commVec%myrank
    ! call set_GPU(myrank)
    call set_cublas(method%cublas_solver)

#if 0
    if(allocated(method%ddmMat%mat%ioff_dia))then
       n = size(method%ddmMat%mat%ioff_dia)
       call cudaMalloc_int(method%ddmMat%mat%C_ioff_dia,n)
       call cudaMemcpy_h2d_int( &
            method%ddmMat%mat%ioff_dia, &
            method%ddmMat%mat%C_ioff_dia,n,0,0)
    endif

    if(allocated(method%ddmMat%mat%val_dia))then
       n = size(method%ddmMat%mat%val_dia)
       call cudaMalloc_fp64(method%ddmMat%mat%C_val_dia,n)
       call cudaMemcpy_h2d_fp64( &
            method%ddmMat%mat%val_dia, &
            method%ddmMat%mat%C_val_dia,n,0,0)
    endif

#endif

    n = method%prm%n
    call cudaMalloc_fp64(method%ddmMat%C_work,n)

    ! call parcel_sync(method);write(*,*)"@ 5 @"
    ! MPI communication
    call sleev_comm_init_GPU(method)
    call malloc_MPI_Comm(method%commVec%C_comm,1) 
    call MPI_Comm_f2c(method%commVec%C_comm,method%commVec%comm) 

    call malloc_MPI_Status( method%commVec%C_s_sta,6)
    call malloc_MPI_Request(method%commVec%C_s_req,6)

    call malloc_MPI_Status( method%commVec%C_r_sta,6)
    call malloc_MPI_Request(method%commVec%C_r_req,6)

  end subroutine wrap_init_ddm_GPU

  subroutine wrap_free_ddm_GPU(method)
    use cu_common_wrapper
    use krylov  
    use PARCEL_integer
    implicit none
    type(KSP) :: method
    integer (kind=intType) i
    integer (kind=intType) n

    if(allocated(method%ddmMat%mat%val_dia))then
       call cudaFree(method%ddmMat%mat%C_val_dia)
    endif
    if(allocated(method%ddmMat%mat%ioff_dia))then
       call cudaFree(method%ddmMat%mat%C_ioff_dia)
    endif

    !    if(allocated(method%ddmMat%diag))then
    !    call cudaFree(method%ddmMat%C_diag)
    call cudaFree(method%ddmMat%C_work)
    !    endif

    ! MPI communication
    call cudaFree(method%ddmMat%mat%C_sbuf)
    call cudaFree(method%ddmMat%mat%C_rbuf)
 
    call free(method%commVec%C_comm);
    call free(method%commVec%C_s_sta);
    call free(method%commVec%C_s_req);
    call free(method%commVec%C_r_sta);
    call free(method%commVec%C_r_req);

    call free_cublas(method%cublas_solver)

  end subroutine wrap_free_ddm_GPU

  subroutine precon_main_ddm(&
       ipc, n, ndim, nnd, diag, val_dia, ioff_dia,&
       v1, v2, prm, ncomp_grids, ioff_grids, s, &
       NEWCOMM, my_rank, AS_n, AS_prm,  &
       AS_mean, AS_g, AS_val_dia, AS_ioff_dia, AS_v1, AS_v2,precon, &
       PRECISION_CALC, &
       PRECISION_PRECONDITION, &
       diag_lo, &
       val_dia_lo, &
       v1_lo, &
       v2_lo, &
       AS_val_dia_lo, &
       AS_v1_lo, &
       AS_v2_lo &
       )
    use param_AS_type_ddm
    use mpi
    use ddmMat_wrap
    use underlap,only : underlap_ddm
    use AS_util_ddm,only : AS_util_copy_to_AS,AS_util_copy_from_AS
    use as_ddm,only : aslu_ddm,aslu_ddm_th
    use PARCEL_integer

    implicit none
    ! inout
    integer(kind=intType) n,ndim,AS_n

    integer ipc, nnd,  myrank, icomm, s
    integer ncomp_grids, num_grids(ncomp_grids)
    integer ioff_dia(nnd), ioff_grids(nnd,ncomp_grids)
    integer margin
    double precision diag(n), val_dia(nnd*n), v1(n), v2(n)
    type(ddm_param) :: prm, AS_prm
!    integer :: div_direc_th
    !
    integer NEWCOMM, my_rank, AS_mean
    type(AS_geom) :: AS_g
    double precision AS_val_dia(AS_n*nnd), AS_v1(AS_n), AS_v2(AS_n)
    integer AS_ioff_dia(nnd)
    type(precon_AS_wrap) precon
    integer iflag_AS

    integer,optional :: PRECISION_CALC
    integer,optional :: PRECISION_PRECONDITION
    double precision,optional :: diag_lo(n)
    double precision,optional :: val_dia_lo(nnd*n)
    double precision,optional :: v1_lo(n)
    double precision,optional :: v2_lo(n)
    double precision,optional :: AS_val_dia_lo(AS_n*nnd)
    double precision,optional :: AS_v1_lo(AS_n)
    double precision,optional :: AS_v2_lo(AS_n)

    ! local
    integer i, iend
#if 1
    integer nx, ny, nz, mx, my, mz, mxy, jx, jy, jz
    integer nxy, nxyz, j
    integer stmx, stmy, stmz      
#endif
    real*8 ddmul_var1
    integer PRECISION
    PRECISION = KSP_DOUBLE_PRECISION
    if(  present(PRECISION_CALC)  )then
       PRECISION = PRECISION_CALC
    endif

!    write(*,*)"@ AS_n,n = ",AS_n,n
    if (ipc .eq. 0) then            !! none
       if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
          !$omp parallel do
          do i = 1, n
             v2(i)    = v1(i)
             v2_lo(i) = v1_lo(i)
          end do
          !$omp end parallel do
       else
          !$omp parallel do
          do i = 1, n
             v2(i) = v1(i)
          end do
          !$omp end parallel do
       endif
    else if (ipc .eq. 1) then       !! point jacobi
       if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
          if( PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
             !$omp parallel do private( &
             !$omp i, &
             !$omp ddmul_var1  &
             !$omp )
             do i = 1, n
                !v2(i) = v1(i) * diag(i)
                ddmul_calc1(v1(i),v1_lo(i),diag(i),diag_lo(i),v2(i),v2_lo(i))
             end do
             !$omp end parallel do
          else
             !$omp parallel do private( &
             !$omp i, &
             !$omp ddmul_var1  &
             !$omp )
             do i = 1, n
                !v2(i) = v1(i) * diag(i)
                ddmul_calc1(v1(i),v1_lo(i),diag(i),0.0d0,v2(i),v2_lo(i))
             end do
             !$omp end parallel do
          endif
       else
          !$omp parallel do
          do i = 1, n
             v2(i) = v1(i) * diag(i)
          end do
          !$omp end parallel do
       endif
    else if (ipc .eq. 2) then !! block jacobi
       if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
          if( PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
             call aslu_ddm( n, nnd, val_dia, ioff_dia,&
                  diag, v1, v2, prm, ncomp_grids, num_grids, ioff_grids, precon, &
                  PRECISION_CALC, &
                  PRECISION_PRECONDITION , &
                  val_dia_lo, &
                  v1_lo, &
                  v2_lo, &
                  diag_lo &
                  )
          else
             call aslu_ddm( n, nnd, val_dia, ioff_dia,&
                  diag, v1, v2, prm, ncomp_grids, num_grids, ioff_grids, precon, &
                  PRECISION_CALC, &
                  PRECISION_PRECONDITION , &
                  val_dia_lo, &
                  vec1_lo=v1_lo, &
                  vec2_lo=v2_lo  &
                  )
          endif
       else
          call aslu_ddm( n, nnd, val_dia, ioff_dia,&
               diag, v1, v2, prm, ncomp_grids, num_grids, ioff_grids, precon )
       endif
    else if ( (ipc.eq.3) .or. (ipc.eq.DDM_CUDA_AS) ) then  !! additive schwartz
       !     write(*,*)"AS_n, nnd=",AS_n, nnd
       if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
!       if( PRECISION.eq.-1 )then

             iflag_AS=precon%AS_METHOD_TYPE
             call AS_util_copy_to_AS ( NEWCOMM, my_rank, n, prm, AS_n, AS_prm, &
                  1, AS_g, v1, AS_v1,iflag_AS )
             call AS_util_copy_to_AS ( NEWCOMM, my_rank, n, prm, AS_n, AS_prm, &
                  1, AS_g, v1_lo, AS_v1_lo,iflag_AS )
#if 1
             call aslu_ddm_th( AS_n, nnd, AS_val_dia, AS_ioff_dia, diag, AS_v1, AS_v2, &
                  AS_prm, ncomp_grids, num_grids, ioff_grids, AS_g%addL ,precon , &
                  PRECISION_CALC, &
                  PRECISION_PRECONDITION, &
                  AS_val_dia_lo, &
                  AS_v1_lo, &
                  AS_v2_lo  &
             )
#else
             call aslu_ddm_th( AS_n, nnd, AS_val_dia, AS_ioff_dia, diag, AS_v1, AS_v2, &
                  AS_prm, ncomp_grids, num_grids, ioff_grids, AS_g%addL ,precon  &
             )
             
#endif             

             if(&
                  (iflag_AS.eq.KSP_ADDITIVESCHWARZ_TYPE_INTERPOLATE)&
                  .or. &
                  (iflag_AS.eq.KSP_ADDITIVESCHWARZ_TYPE_BASIC)&
                  )then
                AS_mean=1
             else
                AS_mean=0
             endif

             !$omp parallel do
             do i = 1, n
                v2(i) = v1(i)
                v2_lo(i) = v1_lo(i)
             end do
             !$omp end parallel do
             
#if 1
             call AS_util_copy_from_AS ( NEWCOMM, my_rank, n, prm, AS_n, AS_prm, &
                  1, AS_mean, AS_g, v2, AS_v2 , &
                  PRECISION, & 
                  v2_lo, &
                  AS_v2_lo)
#else
             call AS_util_copy_from_AS ( NEWCOMM, my_rank, n, prm, AS_n, AS_prm, &
                  1, AS_mean, AS_g, v2, AS_v2 )
#endif
       else
          iflag_AS=precon%AS_METHOD_TYPE
          call AS_util_copy_to_AS ( NEWCOMM, my_rank, n, prm, AS_n, AS_prm, &
               1, AS_g, v1, AS_v1,iflag_AS )

          call aslu_ddm_th( AS_n, nnd, AS_val_dia, AS_ioff_dia, diag, AS_v1, AS_v2, &
               AS_prm, ncomp_grids, num_grids, ioff_grids, AS_g%addL ,precon )

          if(&
               (iflag_AS.eq.KSP_ADDITIVESCHWARZ_TYPE_INTERPOLATE)&
               .or. &
               (iflag_AS.eq.KSP_ADDITIVESCHWARZ_TYPE_BASIC)&
               )then
             AS_mean=1
          else
             AS_mean=0
          endif

          !$omp parallel do
          do i = 1, n
             v2(i) = v1(i)
          end do
          !$omp end parallel do

          call AS_util_copy_from_AS ( NEWCOMM, my_rank, n, prm, AS_n, AS_prm, &
               !            & AS_diag_neighbor, AS_mean, AS_g, v2, AS_v2 )
               1, AS_mean, AS_g, v2, AS_v2 )
       endif

    else if (ipc .eq. 4 ) then      !! underlap
       write(*,*)"not support : precon_main_ddm "
       ! call underlap_ddm( n, nnd, val_dia, ioff_dia, diag, v1, v2,&
       !      prm, ncomp_grids, num_grids, ioff_grids) 

    end if

  end subroutine precon_main_ddm

  subroutine parcel_ddm_main( &
       icomm,ndim, nrow, &
       nnd,&
#if CUDA_SOLVER
       ioff_dia_ptr, val_dia_ptr, &
#else
       ioff_dia, val_dia, &
#endif
       num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
       ndiv_grids, num_grids, ioff_grids, div_direc_th, precon_thblock, &
       iovlflag, nas, iflag_AS,&
#if CUDA_SOLVER
       bptr, xptr, &
#else
       b, x, &
#endif
       s,t, &
       tol, abstol, bflag, qrflag, pc, &
       nres, res, &
       Neuman_itr,&
       iflag_DDM_independ, &
       nBlock, &
       iret, &
       GMRES_m, &
       GMRES_GSflag, &
       CBCG_kstep, &
       CBCG_Eigenflag, &
       power_method_itrmax, &
       CAARNOLDI_sstep, &
       CAARNOLDI_tstep, &
       CAARNOLDI_basis, &
       CAARNOLDI_QRflag, &
       CAGMRES_sstep, &
       CAGMRES_tstep, &
       CAGMRES_basis, &
       CAGMRES_QRflag, &
       solver, &
! -----
! optional
       PRECISION_KRYLOV, & 
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION, &
       val_dia_lo,b_lo, x_lo, & 
       EigenType &
       )
    use krylov  
    use DDM_TYPE
    use init_ddm_param, only: &
         parcel_ddm_init, &
         parcel_ddm_free

    use PARCEL_define
    use vec_kernel
    use vec_ddm_kernel
!    use krylov_kernel,only :get_EigenPair_count_KSP,get_EigenPair_KSP
    use EigenSolver
    use mpi
    use precon_DILU
    use PARCEL_integer
    use cu_Vec,only : &
         cuda_xcpy_ddm_FP64_GPU,zero_init_FP64_GPU, &
         cuda_get_diamatrix_neumanndiag_FP64_GPU, &
         cuda_get_diamatrix_diag_FP64_GPU

    use PARCEL_VEC_kernel,only : parcel_sync

    implicit none

    !arguments
    integer, intent(in) :: icomm !MPI communicator
    integer, intent(in) :: ndim !dimension of global matrix
    integer(kind=intType), intent(in) :: nrow !number of local matrix rows

#ifdef CUDA_SOLVER
    ! type(c_ptr) b !b vector
    ! type(c_ptr) x !x vector
    type(c_ptr) bptr !b vector
    type(c_ptr) xptr !x vector
    double precision b(nrow) !b vector
    double precision x(nrow) !x vector
#else
    double precision, intent(in) :: b(nrow) !b vector
    double precision, intent(inout) :: x(nrow) !x vector
#endif
    integer, intent(in) :: s ,t
!    integer, intent(inout) :: oloop
    double precision, intent(inout) :: tol
    double precision, intent(inout) :: abstol
    integer, intent(in) :: bflag 
    integer, intent(in) :: qrflag 
    integer  pc !pre-conditioner flag
    integer, intent(inout) :: nres !number of residual history
    double precision,intent(inout) :: res(nres) !residual history

    integer, intent(in) :: nnd 
#ifdef CUDA_SOLVER
    type(c_ptr) ioff_dia_ptr ! offset
    type(c_ptr) val_dia_ptr  ! values
    integer ioff_dia(nnd) !offset
    integer ioff_dia_neuman(nnd) !offset
    double precision val_dia(nrow*nnd) !values
#else
    integer, intent(in) :: ioff_dia(nnd) !offset
    double precision, intent(in) :: val_dia(nrow*nnd) !values
#endif
    integer, intent(in) :: num_neighbor_ranks
    integer, intent(in) :: neighbor_ranks( num_neighbor_ranks )
    integer, intent(in) :: margin
    integer, intent(in) :: ncomp_grids
    integer, intent(in) :: ndiv_grids(ncomp_grids)
    integer, intent(in) :: num_grids(ncomp_grids)
    integer, intent(in) :: ioff_grids(nnd,ncomp_grids)
    integer, intent(in) :: div_direc_th
    integer, intent(in) :: precon_thblock
!    integer, intent(in) :: commhide ! 0: no hiding, 1: one time, 2: each time
    integer, intent(in) :: iovlflag ! 0: no hiding, 1: one time, 2: each time
    integer, intent(inout) :: nas
    integer, intent(in) :: iflag_AS         
    integer, intent(in) ::  iflag_DDM_independ
    integer, intent(in) :: nBlock
    integer, intent(inout) :: iret
    integer, intent(in) :: Neuman_itr
    integer leng_halo
    integer AS_ioff_dia(nnd)

    integer  nprocs, myrank, ncpu, newcomm, ierr,nmax
    integer(kind=intType)  iend,  i, j
    integer  as_n(4), nn_l, nn_u, &
         nsend1, nsend2, nrecv1, nrecv2, n_blk, nbuf
    double precision dummy,t_all,t_precon_init,t_cagmres_main

    integer GMRES_m         
    integer GMRES_GSflag    
    integer solver          
    integer CBCG_kstep
    integer CBCG_Eigenflag
    integer power_method_itrmax
    integer CAARNOLDI_sstep
    integer CAARNOLDI_tstep
    integer CAARNOLDI_basis
    integer CAARNOLDI_QRflag
    integer CAGMRES_sstep
    integer CAGMRES_tstep
    integer CAGMRES_basis
    integer CAGMRES_QRflag

    integer,optional :: PRECISION_KRYLOV
    integer,optional :: PRECISION_A 
    integer,optional :: PRECISION_b 
    integer,optional :: PRECISION_x 
    integer,optional :: PRECISION_PRECONDITION 
    double precision,optional :: val_dia_lo( nrow*nnd) ! Matrix values low bit
    double precision,optional :: b_lo(nrow) ! b vector low bit
    double precision,optional :: x_lo(nrow) ! x vector low bit
    type(Eigen_Type),optional :: EigenType    

    type(KSP) :: method
    integer(kind=intType) n_tmp
    
    integer PRECISION_A_DDM
    integer PRECISION_b_DDM
    integer PRECISION_x_DDM
    integer PRECISION_KRYLOV_DDM
    integer PRECISION_PRECONDITION_DDM ! 前処理行列の精度 (前処理の精度ではない)
    integer ncount
    integer commhide
    integer(kind=intType) parcel_int
    integer(kind=intType),parameter :: offset_zero=0

    double precision stime0,etime0
    double precision stime,etime

#ifdef CUDA_SOLVER
    real*8,pointer :: vecx_Fptr(:)
    real*8,pointer :: vecb_Fptr(:)
    real*8,pointer :: diaA_Fptr(:)
!    integer(kind=intType),pointer :: offset_Fptr(:)
    integer,pointer :: offset_Fptr(:)
    type(devicemem) cu_work
#endif
    integer(kind=intType) idiag 
    integer(kind=intType) n_parcel_int
    integer(kind=intType) nnd_parcel_int

#if CUDA_SOLVER

#else
    call MPI_Comm_rank(icomm, myrank, ierr)

    if( (pc.eq.KSP_PRECONDITION_NEUMANN) .or. &
         (pc.eq.KSP_PRECONDITION_NEUMANN_LOCAL) )then 
       if(myrank.eq.0)then
          write(*,*)"not support Neuman Preconditione"
       endif
       nres = 1
       return
    endif
#endif

#if CUDA_SOLVER
!#if 1

    call MPI_Comm_rank(icomm, myrank, ierr)

    stime0 = MPI_Wtime()
    stime = MPI_Wtime()

#ifdef CUDA_SOLVER

    ! DDM Matrix
    parcel_int = nrow*nnd
    call c_f_pointer(val_dia_ptr,diaA_Fptr,[parcel_int])
    call cudaMalloc_fp64(method%ddmMat%mat%C_val_dia,parcel_int)
    call cudaMemcpy_dptr2d_double(diaA_Fptr,method%ddmMat%mat%C_val_dia,parcel_int)
    call cudaMemcpy_d2h_fp64(method%ddmMat%mat%C_val_dia,val_dia,parcel_int,offset_zero,offset_zero)

    parcel_int = nnd
    call c_f_pointer(ioff_dia_ptr,offset_Fptr,[parcel_int])
    call cudaMalloc_int(method%ddmMat%mat%C_ioff_dia,parcel_int)
    call cudaMemcpy_dptr2d_int(offset_Fptr,method%ddmMat%mat%C_ioff_dia,parcel_int)
    call cudaMemcpy_d2h_int(method%ddmMat%mat%C_ioff_dia,ioff_dia,parcel_int,offset_zero,offset_zero)

!    call cudafree(method%ddmMat%mat%C_val_dia)
!    call cudafree(method%C_diaMat%offset)

#endif

    etime = MPI_Wtime()
#if PRINT_TIME
    if(myrank.eq.0)then
       write(*,*)"     @ 1-1 @",etime - stime
    endif
#endif
    stime = MPI_Wtime()

    ! --------------
    ! flag
    if(PARCEL_MACHINE_FLAG_DEFAULT.eq.PARCEL_GPU_FLAG_ON)then
       if(  present(PRECISION_KRYLOV)  )then
       else
          if(pc.eq.KSP_PRECONDITION_BLOCKJACOBI)then
             pc = DDM_CUDA_BJ
          endif
          if(pc.eq.KSP_PRECONDITION_ADDITIVESCHWARZ)then
             pc = DDM_CUDA_AS
          endif
       endif
    else
       if( (pc.eq.KSP_PRECONDITION_NEUMANN) .or. &
           (pc.eq.KSP_PRECONDITION_NEUMANN_LOCAL) )then 
          pc = KSP_PRECONDITION_NONE
       endif
    endif

    ! --------------
    ! init
    commhide = 0
    if(  present(PRECISION_KRYLOV)  )then
       ! 4倍精度
       call parcel_ddm_init( &
            method%ddmMat, &
            icomm,&
            ndim,  &
            nrow, &
            nnd, &
            ioff_dia, &
            val_dia, &
            num_neighbor_ranks, &
            neighbor_ranks, &
            margin, &
            ncomp_grids, &
            ndiv_grids, num_grids, ioff_grids, &
            div_direc_th, precon_thblock, &
            commhide, nas, iflag_AS,&
            b, x, s, t, nres, tol, bflag, qrflag, pc, &
            iflag_DDM_independ, &
            nres, res,&
            PRECISION_KRYLOV, & 
            PRECISION_A, &
            PRECISION_b, &
            PRECISION_x, &
            PRECISION_PRECONDITION, &
            val_dia_lo,b_lo, x_lo & 
            )
    else
       ! 倍精度
       call parcel_ddm_init( &
            method%ddmMat, &
            icomm,&
            ndim,  &
            nrow, &
            nnd, &
            ioff_dia, &
            val_dia, &
            num_neighbor_ranks, &
            neighbor_ranks, &
            margin, &
            ncomp_grids, &
            ndiv_grids, num_grids, ioff_grids, &
            div_direc_th, precon_thblock, &
            commhide, nas, iflag_AS,&
            b, x, s, t, nres, tol, bflag, qrflag, pc, &
            iflag_DDM_independ, &
            nres, res &
            )
    endif

    ! --------------
    ! solverの精度決定
    PRECISION_A_DDM            = KSP_DOUBLE_PRECISION
    PRECISION_b_DDM            = KSP_DOUBLE_PRECISION
    PRECISION_x_DDM            = KSP_DOUBLE_PRECISION
    PRECISION_KRYLOV_DDM       = KSP_DOUBLE_PRECISION
    PRECISION_PRECONDITION_DDM = KSP_DOUBLE_PRECISION

    if(  present(PRECISION_A)  )then
       PRECISION_A_DDM = PRECISION_A
    endif
    if(  present(PRECISION_b)  )then
       PRECISION_b_DDM = PRECISION_b
    endif
    if(  present(PRECISION_x)  )then
       PRECISION_x_DDM = PRECISION_x
    endif
    if(  present(PRECISION_PRECONDITION)  )then
       PRECISION_PRECONDITION_DDM = PRECISION_PRECONDITION
    endif
    if(  present(PRECISION_KRYLOV)  )then
       PRECISION_KRYLOV_DDM = PRECISION_KRYLOV
    endif
    if(PRECISION_KRYLOV_DDM.eq.KSP_DOUBLE_PRECISION)then
       PRECISION_A_DDM            = KSP_DOUBLE_PRECISION
       PRECISION_b_DDM            = KSP_DOUBLE_PRECISION
       PRECISION_x_DDM            = KSP_DOUBLE_PRECISION
       PRECISION_PRECONDITION_DDM = KSP_DOUBLE_PRECISION
    endif

    etime = MPI_Wtime()
#if PRINT_TIME
    if(myrank.eq.0)then
       write(*,*)"     @ 1-2 @",etime - stime
    endif
#endif
    stime = MPI_Wtime()

    ! ------------
    call PARCEL_KSP_Default_Setting(method)

    call MPI_Comm_size(icomm, nprocs, ierr)
    call MPI_Comm_rank(icomm, myrank, ierr)
    ncpu = nprocs

    if(PRECISION_PRECONDITION_DDM.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       call wrap_init_ddm(&
            method%ddmMat, &
            method%ddmMat%precon_ddm, &
            ndim, nrow, nmax, nnd, &
            method%ddmMat%mat, &
            ncpu, myrank, icomm, &           
            method%ddmMat%ASval%diag_init, &
            -1, tol, pc, nres, res, iret, &
            commHide, nas, iflag_AS, &
            method%ddmMat%ASval%AS_num, &
            method%ddmMat%AS_g, &
            method%ddmMat%ASval%AS_val_dia_init, &
            method%ddmMat%ASval%AS_ioff_dia_init, &
            method%ddmMat%AS_prm , &
            method%ddmMat%ASval%diag_init_lo, &
            method%ddmMat%ASval%AS_val_dia_init_lo &
            )
    else
       call wrap_init_ddm(&
            method%ddmMat, &
            method%ddmMat%precon_ddm, &
            ndim, nrow, nmax, nnd,  &
            method%ddmMat%mat, &
            ncpu, myrank, icomm, &
            method%ddmMat%ASval%diag_init, &
            -1, tol, pc, nres, res, iret, &
            commHide, nas, iflag_AS, &
            method%ddmMat%ASval%AS_num, &
            method%ddmMat%AS_g, &
            method%ddmMat%ASval%AS_val_dia_init, &
            method%ddmMat%ASval%AS_ioff_dia_init, &
            method%ddmMat%AS_prm  &
            )
    endif

    etime = MPI_Wtime()
#if PRINT_TIME
    if(myrank.eq.0)then
       write(*,*)"     @ 1-3 @",etime - stime
    endif
#endif
    stime = MPI_Wtime()

#ifdef CUDA_SOLVER

    ! vector input
    ! xcpy_ddm --> ddm 形式の入力(x,b)は袖領域をゼロにする必要がある.
    parcel_int = nrow
    call c_f_pointer(xptr,vecx_Fptr,[parcel_int])
    call cudaMalloc_fp64(cu_work,parcel_int)

    call cudaMemcpy_dptr2d_double(vecx_Fptr,cu_work,parcel_int)
    call cudaMalloc_fp64(method%C_x,parcel_int)
    call parcel_sync(method)

    call zero_init_FP64_GPU(parcel_int,method%C_x,offset_zero)
    call cuda_xcpy_ddm_FP64_GPU( &
         method%ddmMat%prm%nx  , &
         method%ddmMat%prm%ny  , &
         method%ddmMat%prm%nz  , &
         method%ddmMat%prm%stmx, &
         method%ddmMat%prm%stmy, &
         method%ddmMat%prm%stmz, &
         method%ddmMat%prm%mx  , &
         method%ddmMat%prm%mxy , &
         cu_work,offset_zero,    &
         method%C_x,offset_zero  &
         )
    ! call parcel_sync(method)
    ! parcel_int = nrow
    ! call cudaMemcpy_d2h_fp64(method%C_x,x,parcel_int,offset_zero,offset_zero)

    parcel_int = nrow
    call c_f_pointer(bptr,vecb_Fptr,[parcel_int])
    call cudaMemcpy_dptr2d_double(vecb_Fptr,cu_work,parcel_int)
    call cudaMalloc_fp64(method%C_b,parcel_int)
    call parcel_sync(method)

    call zero_init_FP64_GPU(parcel_int,method%C_b,offset_zero)
    call cuda_xcpy_ddm_FP64_GPU( &
         method%ddmMat%prm%nx  , &
         method%ddmMat%prm%ny  , &
         method%ddmMat%prm%nz  , &
         method%ddmMat%prm%stmx, &
         method%ddmMat%prm%stmy, &
         method%ddmMat%prm%stmz, &
         method%ddmMat%prm%mx  , &
         method%ddmMat%prm%mxy , &
         cu_work   ,offset_zero, &
         method%C_b,offset_zero  &
         )

    ! call parcel_sync(method)
    ! call cudaMemcpy_d2h_fp64(method%C_b,b,parcel_int,offset_zero,offset_zero)
    ! call parcel_sync(method)
#if 0
    do i = 1,nrow
       write(*,*)"vecb : ",i,b(i) ,x(i) 
    enddo
    do j = 0,nnd-1
       do i = 1,nrow
          if(val_dia(i+nrow*j) .ne. 0.0d0 )then
             write(*,*)"diaA : ",j,i,val_dia(i+nrow*j)
          endif
       enddo
    enddo
#endif

    call cudafree(cu_work)

#endif

    etime = MPI_Wtime()
#if PRINT_TIME
    if(myrank.eq.0)then
       write(*,*)"     @ 1-4 @",etime - stime
    endif
#endif
    stime = MPI_Wtime()

    method%ddmMat%k = s
    method%ddmMat%sstep  = s
    method%ddmMat%tstep  = t
    method%ddmMat%qrflag = qrflag
    method%ddmMat%bflag  = bflag
    method%ddmMat%CBCG_Eigenflag = CBCG_Eigenflag
    method%ddmMat%solver         = solver

    ! -----
    method%PRECISION_A            = PRECISION_A_DDM
    method%PRECISION_b            = PRECISION_b_DDM
    method%PRECISION_x            = PRECISION_x_DDM
    method%PRECISION_KRYLOV       = PRECISION_KRYLOV_DDM
    method%PRECISION_PRECONDITION = PRECISION_PRECONDITION_DDM
    method%CAGMRES_sstep  = s
    method%CAGMRES_tstep  = t
    method%CAGMRES_QRflag = qrflag
    method%CAGMRES_basis  = bflag
    method%CBCG_Eigenflag = CBCG_Eigenflag
    method%solver         = solver
    method%GMRES_m        = GMRES_m
    method%GMRES_GSflag   = GMRES_GSflag

    method%commVec%myrank=myrank
    method%commVec%comm=icomm
    method%iovlflag = commhide
    method%itrmax=nres

#if 1
    if(pc.eq.KSP_PRECONDITION_NEUMANN_LOCAL)then
       pc = KSP_PRECONDITION_NEUMANN
       method%C_preconNeuman%local_flag = -1
    endif
    method%C_preconNeuman%itr = Neuman_itr
    method%ipreflag = pc
#else
    method%ipreflag = pc
#endif

    allocate(method%reshistory(max(nres+100,100)))

#ifdef CUDA_SOLVER

    select case (method%ipreflag)
    case(KSP_PRECONDITION_NONE)

    case(KSP_PRECONDITION_POINTJACOBI)

       call cudaMalloc_fp64(method%C_D,nrow)
       call cudaMalloc_fp64(method%C_Dinv,nrow)
       idiag = 1
       do i = 1,nnd
          if(ioff_dia(i).eq.0)then
             idiag = i
             exit
          endif
       enddo

       n_parcel_int   = nrow
       nnd_parcel_int = nnd
       call cuda_get_diamatrix_diag_FP64_GPU( &
          n_parcel_int, &
          nnd_parcel_int, &
          idiag, &
          method%ddmMat%mat%C_val_dia, &
          method%C_D, &
          offset_zero, &
          method%C_Dinv, &
          offset_zero &
          )

    case(DDM_CUDA_BJ)
       ! 問題サイズの除算値が1ブロックサイズの目安
       call init_DILU_DDM( &
            method,&
            method%ddmMat%prm%nx/4 ,&
            method%ddmMat%prm%ny/4 ,&
            method%ddmMat%prm%nz/4  &
            )
       call get_DILU_DDM(method)
    case(DDM_CUDA_AS)

       ! do i =1,num_neighbor_ranks
       !    write(*,*)i,neighbor_ranks(i)
       ! enddo

       call init_AS_DDM( &
            method,&
#if 0
            method%ddmMat%prm%nx ,&
            method%ddmMat%prm%ny ,&
            method%ddmMat%prm%nz  &
#else
            method%ddmMat%prm%nx/4 ,&
            method%ddmMat%prm%ny/4 ,&
            method%ddmMat%prm%nz/4  &
#endif
            )
       call get_DILU_AS_DDM(method)
    case(KSP_PRECONDITION_NEUMANN)

       parcel_int = nrow*nnd
       call cudaMalloc_fp64(method%ddmMat%C_DDMMat_Neumann%ddmA   , parcel_int)
       call cudaMemcpy_d2d_FP64(&
            method%ddmMat%mat%C_val_dia, &
            method%ddmMat%C_DDMMat_Neumann%ddmA , parcel_int,offset_zero,offset_zero)

       parcel_int = nnd
       call cudaMalloc_int(method%ddmMat%C_DDMMat_Neumann%offset , parcel_int)

       parcel_int = nrow
       call cudaMalloc_fp64(method%C_preconNeuman%Dinv , parcel_int)
       call cudaMalloc_fp64(method%C_preconNeuman%vecx , parcel_int)
       call cudaMalloc_fp64(method%C_preconNeuman%vecy , parcel_int) 

       parcel_int = method%ddmMat%C_DDMMat_Neumann%m
       parcel_int = parcel_int * method%ddmMat%C_DDMMat_Neumann%nnd

       parcel_int = method%ddmMat%C_DDMMat_Neumann%nnd

       method%ddmMat%C_DDMMat_Neumann%nx   = method%ddmMat%prm%nx  
       method%ddmMat%C_DDMMat_Neumann%ny   = method%ddmMat%prm%ny  
       method%ddmMat%C_DDMMat_Neumann%nz   = method%ddmMat%prm%nz  
       method%ddmMat%C_DDMMat_Neumann%stmx = method%ddmMat%prm%stmx
       method%ddmMat%C_DDMMat_Neumann%stmy = method%ddmMat%prm%stmy
       method%ddmMat%C_DDMMat_Neumann%stmz = method%ddmMat%prm%stmz
       
       method%ddmMat%C_DDMMat_Neumann%mx   = method%ddmMat%prm%mx  
       method%ddmMat%C_DDMMat_Neumann%my   = method%ddmMat%prm%my  
       method%ddmMat%C_DDMMat_Neumann%mz   = method%ddmMat%prm%mz  
       method%ddmMat%C_DDMMat_Neumann%mxy  = method%ddmMat%prm%mxy 
       method%ddmMat%C_DDMMat_Neumann%m    = nrow

       idiag = -1
       do i = 1, nnd
          if(ioff_dia(i).eq.0)then
             idiag = i
          endif
       enddo

       call cuda_get_ddmmatrix_neumanndiag_FP64_GPU( &
            method%ddmMat%C_DDMMat_Neumann%nx, &
            method%ddmMat%C_DDMMat_Neumann%ny, &
            method%ddmMat%C_DDMMat_Neumann%nz, &
            method%ddmMat%C_DDMMat_Neumann%stmx, &
            method%ddmMat%C_DDMMat_Neumann%stmy, &
            method%ddmMat%C_DDMMat_Neumann%stmz, &
            method%ddmMat%C_DDMMat_Neumann%mx, &
            method%ddmMat%C_DDMMat_Neumann%mxy, &
            method%ddmMat%C_DDMMat_Neumann%m, &
            nnd  , &
            idiag, &
            method%ddmMat%C_DDMMat_Neumann%ddmA, &
            method%C_preconNeuman%Dinv,offset_zero, &
            method%ddmMat%C_DDMMat_Neumann%offset  &
            )
       ioff_dia_neuman        = ioff_dia
       ioff_dia_neuman(idiag) = ioff_dia_neuman(nnd) 

       parcel_int = nnd
       call cudaMemcpy_h2d_int( &
            ioff_dia_neuman, &
            method%ddmMat%C_DDMMat_Neumann%offset, &
            parcel_int, &
            offset_zero, &
            offset_zero)

       method%ddmMat%C_DDMMat_Neumann%nnd  = nnd -1

    end select

#endif

    n_tmp=nres
    call zero_initialize(n_tmp,method%reshistory)

    method%rtolmax   = tol
    method%abstolmax = abstol
    method%prm%n     = nrow

    method%CBCG_kstep     = CBCG_kstep
    method%CBCG_Eigenflag = CBCG_Eigenflag
    method%power_method_itrmax = power_method_itrmax
    method%CAARNOLDI_sstep  = CAARNOLDI_sstep
    method%CAARNOLDI_tstep  = CAARNOLDI_tstep
    method%CAARNOLDI_basis  = CAARNOLDI_basis
    method%CAARNOLDI_QRflag = CAARNOLDI_QRflag
    method%CAGMRES_sstep    = CAGMRES_sstep
    method%CAGMRES_tstep    = CAGMRES_tstep
    method%CAGMRES_basis    = CAGMRES_basis
    method%CAGMRES_QRflag   = CAGMRES_QRflag
    method%nBlock = nBlock

    etime = MPI_Wtime()
#if PRINT_TIME
    if(myrank.eq.0)then
       write(*,*)"     @ 1-5 @",etime - stime
    endif
#endif
    stime = MPI_Wtime()

    method%prm%MatFormat    = MATRIX_FORMAT_DDM
    method%prm_AS%MatFormat = MATRIX_FORMAT_DDM

#ifdef CUDA_SOLVER

#else

    if(method%PRECISION_KRYLOV.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       allocate(method%x_hi(method%prm%n))
       allocate(method%b_hi(method%prm%n))
       call zero_initialize(method%prm%n,method%x_hi)
       call zero_initialize(method%prm%n,method%b_hi)
    else
       allocate(method%x(method%prm%n))
       allocate(method%b(method%prm%n))
       call zero_initialize(method%prm%n,method%x)
       call zero_initialize(method%prm%n,method%b)       
    endif

    if(method%PRECISION_KRYLOV.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       call xcpy_ddm(method%prm%n,x,method%x_hi,method%ddmMat%prm)
       call xcpy_ddm(method%prm%n,b,method%b_hi,method%ddmMat%prm)
    else
       call xcpy_ddm(method%prm%n,x,method%x,method%ddmMat%prm)
       call xcpy_ddm(method%prm%n,b,method%b,method%ddmMat%prm)
    endif

    if(method%PRECISION_x.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       allocate(method%x_lo(method%prm%n))
       call zero_initialize(method%prm%n,method%x_lo)
       if( present(x_lo) )then
          call xcpy_ddm(method%prm%n,x_lo,method%x_lo,method%ddmMat%prm)
       endif
    endif
    if(method%PRECISION_b.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       allocate(method%b_lo(method%prm%n))
       call zero_initialize(method%prm%n,method%b_lo)
       if( present(b_lo) )then
          call xcpy_ddm(method%prm%n,b_lo,method%b_lo,method%ddmMat%prm)
       endif
    endif
 
#endif

    etime = MPI_Wtime()
#if PRINT_TIME
    if(myrank.eq.0)then
       write(*,*)"     @ 1-6 @",etime - stime
    endif
#endif
    stime = MPI_Wtime()

#ifdef CUDA_SOLVER
    if(method%GPUflag.eq.PARCEL_GPU_FLAG_ON)then
       call wrap_init_ddm_GPU(method)
       if(method%ipreflag.eq.DDM_CUDA_BJ)then
          call init_DILU_DDM_GPU(method)
       endif
       if(method%ipreflag.eq.DDM_CUDA_AS)then
          call init_DILU_AS_DDM_GPU(method)
       endif
    endif
#endif

    etime = MPI_Wtime()
#if PRINT_TIME
    if(myrank.eq.0)then
       write(*,*)"     @ 1-7 @",etime - stime
    endif
#endif
    stime = MPI_Wtime()

    if(method%PRECISION_KRYLOV.eq.KSP_DOUBLE_DOUBLE_PRECISION)then

       select case (solver)
       case(KSP_SOLVER_CG)
          call ddKSP_pCG(method)
       case(KSP_SOLVER_BiCG)
          call ddKSP_pBiCG(method)
       case(KSP_SOLVER_GMRES)
          call ddKSP_pGMRES(method)
       case default
          if(myrank.eq.0)then
             write(*,*)"error :: not support solver = ",solver
          end if
       end select

    else

       select case (solver)
       case(KSP_SOLVER_CG)
          call KSP_pCG(method)
       case(KSP_SOLVER_BiCG)
          call KSP_pBiCG(method)
       case(KSP_SOLVER_GMRES)
          call KSP_pGMRES(method)
       case(KSP_SOLVER_CAGMRES)
          call KSP_pCAGMRES(method)
       case(KSP_SOLVER_CBCG)
          call KSP_pCBCG(method)
       case(KSP_SOLVER_POWER_METHOD)
          call power_method_symMatrix(method)
       case(KSP_SOLVER_CAARNOLDI)
          call KSP_pCAARNOLDI(method)
       case default
          if(myrank.eq.0)then
             write(*,*)"error :: not support solver = ",solver
          end if
       end select

#ifdef CUDA_SOLVER
       parcel_int = nrow
       call cudaMemcpy_d2dptr_double(method%C_x,vecx_Fptr,parcel_int)
#endif

    endif

    etime = MPI_Wtime()
#if PRINT_TIME
    if(myrank.eq.0)then
       write(*,*)"     @ 1-8 @",etime - stime
    endif
#endif
    stime = MPI_Wtime()

!    write(*,*)"addL     = ", method%addL
!    write(*,*)"ipreflag = ", method%ipreflag

    if( (solver.eq.KSP_SOLVER_POWER_METHOD).or.(solver.eq.KSP_SOLVER_CAARNOLDI) )then
       call get_EigenPair_count_KSP(ncount,method)
       if(allocated(method%EigenValueReal)) deallocate(method%EigenValueReal)
       if(allocated(method%EigenValueImag)) deallocate(method%EigenValueImag)
       if(allocated(method%EigenVector)) deallocate(method%EigenVector)
       if(allocated(method%EigenErrorNorm)) deallocate(method%EigenErrorNorm)

       allocate( method%EigenValueReal(ncount) )
       allocate( method%EigenValueImag(ncount) )
       allocate( method%EigenVector(ncount*nrow) )
       allocate( method%EigenErrorNorm(ncount) )

       call get_EigenPair_KSP( &           
            ncount,&                             
            method%EigenValueReal,method%EigenValueImag, &     
            method%EigenVector, &                       
            method%EigenErrorNorm, &                    
            method%MaxEigenID,method%MinEigenID,method)


       allocate( EigenType%EigenValueReal(ncount) )
       allocate( EigenType%EigenValueImag(ncount) )
       allocate( EigenType%EigenVector(ncount*nrow) )
       allocate( EigenType%EigenErrorNorm(ncount) )
       EigenType%EigenValueReal = method%EigenValueReal
       EigenType%EigenValueImag = method%EigenValueImag
       EigenType%EigenVector = method%EigenVector
       EigenType%EigenErrorNorm = method%EigenErrorNorm
       EigenType%MaxEigenID = method%MaxEigenID
       EigenType%MinEigenID = method%MinEigenID

    endif

    iret=method%iret
    n_tmp=nres

    !  call xcpy(nres,method%reshistory,res)
    call xcpy(n_tmp,method%reshistory,res)

#ifdef CUDA_SOLVER
    call cudafree(method%C_x)
    call cudafree(method%C_b)
#else
    if(method%PRECISION_KRYLOV.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       call xcpy_ddm(method%prm%n,method%x_hi,x,method%ddmMat%prm)
    else
       call xcpy_ddm(method%prm%n,method%x,x,method%ddmMat%prm)
    endif
    if( present(x_lo) )then
       if( method%PRECISION_x .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
          call xcpy_ddm(method%prm%n,method%x_lo,x_lo,method%ddmMat%prm)
       endif
    endif
#endif

!    oloop=method%niter
    nres = method%niter

    if(method%GPUflag.eq.PARCEL_GPU_FLAG_ON)then
       select case (method%ipreflag)
       case(KSP_PRECONDITION_NONE)

       case(KSP_PRECONDITION_POINTJACOBI)
          call cudafree(method%C_D)
          call cudafree(method%C_Dinv)
       case(DDM_CUDA_BJ)
          call free_DILU_DDM_GPU(method)
       case(DDM_CUDA_AS)
          call free_DILU_AS_DDM_GPU(method)
       case(KSP_PRECONDITION_NEUMANN)
       call cudafree(method%ddmMat%C_DDMMat_Neumann%ddmA  )
       call cudafree(method%ddmMat%C_DDMMat_Neumann%offset)
       call cudafree(method%C_preconNeuman%Dinv)
       call cudafree(method%C_preconNeuman%vecx)
       call cudafree(method%C_preconNeuman%vecy)
    end select
    
    call wrap_free_ddm_GPU(method)
    endif

#if 0
    if(method%ipreflag.eq.DDM_CUDA_BJ)then
       call free_DILU_DDM(method)
    endif

    if ( (method%ipreflag.eq.3) .or. &
         (method%ipreflag.eq.DDM_CUDA_AS) ) then  !! additive schwartz
       call free_DILU_DDM(method)
    endif
#endif

    call wrap_free_ddm(method%ddmMat)
!    if(allocated(@)) deallocate( @ )
    if(allocated(method%reshistory)) deallocate( method%reshistory )
    if(allocated(method%x)) deallocate( method%x )
    if(allocated(method%b)) deallocate( method%b )

    if(allocated(method%x_hi)) deallocate( method%x_hi )
    if(allocated(method%b_hi)) deallocate( method%b_hi )

    if(allocated(method%x_lo)) deallocate( method%x_lo )
    if(allocated(method%b_lo)) deallocate( method%b_lo )

    call parcel_ddm_free(pc,method%ddmMat)

    etime = MPI_Wtime()
#if PRINT_TIME
    if(myrank.eq.0)then
       write(*,*)"     @ 1-9 @",etime - stime
    endif
#endif
    stime = MPI_Wtime()

    etime0 = MPI_Wtime()
#if PRINT_TIME
    if(myrank.eq.0)then
       write(*,*)"@ 1 @",etime0 - stime0
    endif
#endif

#else

    ! --------------
    ! init
    commhide = 0
    if(  present(PRECISION_KRYLOV)  )then
       ! 4倍精度
       call parcel_ddm_init( &
            method%ddmMat, &
            icomm,&
            ndim,  &
            nrow, &
            nnd, &
            ioff_dia, &
            val_dia, &
            num_neighbor_ranks, &
            neighbor_ranks, &
            margin, &
            ncomp_grids, &
            ndiv_grids, num_grids, ioff_grids, &
            div_direc_th, precon_thblock, &
            commhide, nas, iflag_AS,&
            b, x, s, t, nres, tol, bflag, qrflag, pc, &
            iflag_DDM_independ, &
            nres, res,&
            PRECISION_KRYLOV, & 
            PRECISION_A, &
            PRECISION_b, &
            PRECISION_x, &
            PRECISION_PRECONDITION, &
            val_dia_lo,b_lo, x_lo & 
            )
    else
       ! 倍精度
       call parcel_ddm_init( &
            method%ddmMat, &
            icomm,&
            ndim,  &
            nrow, &
            nnd, &
            ioff_dia, &
            val_dia, &
            num_neighbor_ranks, &
            neighbor_ranks, &
            margin, &
            ncomp_grids, &
            ndiv_grids, num_grids, ioff_grids, &
            div_direc_th, precon_thblock, &
            commhide, nas, iflag_AS,&
            b, x, s, t, nres, tol, bflag, qrflag, pc, &
            iflag_DDM_independ, &
            nres, res &
            )
    endif

    ! --------------
    ! solverの精度決定
    PRECISION_A_DDM            = KSP_DOUBLE_PRECISION
    PRECISION_b_DDM            = KSP_DOUBLE_PRECISION
    PRECISION_x_DDM            = KSP_DOUBLE_PRECISION
    PRECISION_KRYLOV_DDM       = KSP_DOUBLE_PRECISION
    PRECISION_PRECONDITION_DDM = KSP_DOUBLE_PRECISION

    if(  present(PRECISION_A)  )then
       PRECISION_A_DDM = PRECISION_A
    endif
    if(  present(PRECISION_b)  )then
       PRECISION_b_DDM = PRECISION_b
    endif
    if(  present(PRECISION_x)  )then
       PRECISION_x_DDM = PRECISION_x
    endif
    if(  present(PRECISION_PRECONDITION)  )then
       PRECISION_PRECONDITION_DDM = PRECISION_PRECONDITION
    endif
    if(  present(PRECISION_KRYLOV)  )then
       PRECISION_KRYLOV_DDM = PRECISION_KRYLOV
    endif
    if(PRECISION_KRYLOV_DDM.eq.KSP_DOUBLE_PRECISION)then
       PRECISION_A_DDM            = KSP_DOUBLE_PRECISION
       PRECISION_b_DDM            = KSP_DOUBLE_PRECISION
       PRECISION_x_DDM            = KSP_DOUBLE_PRECISION
       PRECISION_PRECONDITION_DDM = KSP_DOUBLE_PRECISION
    endif

    call PARCEL_KSP_Default_Setting(method)

    call MPI_Comm_size(icomm, nprocs, ierr)
    call MPI_Comm_rank(icomm, myrank, ierr)
    ncpu = nprocs

#if 1
    if(PRECISION_PRECONDITION_DDM.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       call wrap_init_ddm(&
            method%ddmMat, &
            method%ddmMat%precon_ddm, &
            ndim, nrow, nmax, nnd, &
            method%ddmMat%mat, &
            ncpu, myrank, icomm, &           
            method%ddmMat%ASval%diag_init, &
            -1, tol, pc, nres, res, iret, &
            commHide, nas, iflag_AS, &
            method%ddmMat%ASval%AS_num, &
            method%ddmMat%AS_g, &
            method%ddmMat%ASval%AS_val_dia_init, &
            method%ddmMat%ASval%AS_ioff_dia_init, &
            method%ddmMat%AS_prm , &
            method%ddmMat%ASval%diag_init_lo, &
            method%ddmMat%ASval%AS_val_dia_init_lo &
            )
    else
       call wrap_init_ddm(&
            method%ddmMat, &
            method%ddmMat%precon_ddm, &
            ndim, nrow, nmax, nnd,  &
            method%ddmMat%mat, &
            ncpu, myrank, icomm, &
            method%ddmMat%ASval%diag_init, &
            -1, tol, pc, nres, res, iret, &
            commHide, nas, iflag_AS, &
            method%ddmMat%ASval%AS_num, &
            method%ddmMat%AS_g, &
            method%ddmMat%ASval%AS_val_dia_init, &
            method%ddmMat%ASval%AS_ioff_dia_init, &
            method%ddmMat%AS_prm  &
            )
    endif

#else

    if(PRECISION_PRECONDITION_DDM.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       call wrap_init_ddm(&
            method%ddmMat, &
            precon_ddm, &
            ndim, nrow, nmax, nnd, &
            mat, &
            ncpu, myrank, icomm, &
            diag_init, &
            -1, tol, pc, nres, res, iret, &
            commHide, nas, iflag_AS, &
            AS_num, AS_g, AS_val_dia_init, AS_ioff_dia_init, AS_prm ,&
            diag_init_lo, AS_val_dia_init_lo &
            )
    else
       call wrap_init_ddm(&
            method%ddmMat, &
!            precon_ddm, &
            method%ddmMat%precon_ddm, &
            ndim, nrow, nmax, nnd,  &
            mat, &
            ncpu, myrank, icomm, &
            diag_init, &
            -1, tol, pc, nres, res, iret, &
            commHide, nas, iflag_AS, &
            AS_num, AS_g, AS_val_dia_init, AS_ioff_dia_init, AS_prm &
            )
    endif
#endif

    method%ddmMat%k=s
    method%ddmMat%sstep=s
    method%ddmMat%tstep=t
    method%ddmMat%qrflag=qrflag
    method%ddmMat%bflag=bflag
    method%ddmMat%CBCG_Eigenflag=CBCG_Eigenflag
    method%ddmMat%solver=solver

    ! -----
    method%PRECISION_A            = PRECISION_A_DDM
    method%PRECISION_b            = PRECISION_b_DDM
    method%PRECISION_x            = PRECISION_x_DDM
    method%PRECISION_KRYLOV       = PRECISION_KRYLOV_DDM
    method%PRECISION_PRECONDITION = PRECISION_PRECONDITION_DDM

    ! -----
    method%CAGMRES_sstep  = s
    method%CAGMRES_tstep  = t
    method%CAGMRES_QRflag = qrflag
    method%CAGMRES_basis  = bflag
    method%CBCG_Eigenflag = CBCG_Eigenflag
    method%solver         = solver
    method%GMRES_m        = GMRES_m
    method%GMRES_GSflag   = GMRES_GSflag
    ! -----

    method%commVec%myrank=myrank
    method%commVec%comm=icomm
    method%iovlflag = commhide
    method%itrmax=nres
    allocate(method%reshistory(nres))

    n_tmp=nres
    call zero_initialize(n_tmp,method%reshistory)

    method%rtolmax=tol
    method%abstolmax=abstol

    method%prm%n=nrow

    method%CBCG_kstep = CBCG_kstep
    method%CBCG_Eigenflag = CBCG_Eigenflag
    method%power_method_itrmax = power_method_itrmax
    method%CAARNOLDI_sstep = CAARNOLDI_sstep
    method%CAARNOLDI_tstep = CAARNOLDI_tstep
    method%CAARNOLDI_basis = CAARNOLDI_basis
    method%CAARNOLDI_QRflag = CAARNOLDI_QRflag
    method%CAGMRES_sstep = CAGMRES_sstep
    method%CAGMRES_tstep = CAGMRES_tstep
    method%CAGMRES_basis = CAGMRES_basis
    method%CAGMRES_QRflag = CAGMRES_QRflag
    method%nBlock = nBlock

    if(method%PRECISION_KRYLOV.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       allocate(method%x_hi(method%prm%n))
       allocate(method%b_hi(method%prm%n))
       call zero_initialize(method%prm%n,method%x_hi)
       call zero_initialize(method%prm%n,method%b_hi)
    else
       allocate(method%x(method%prm%n))
       allocate(method%b(method%prm%n))
       call zero_initialize(method%prm%n,method%x)
       call zero_initialize(method%prm%n,method%b)       
    endif

    method%prm%MatFormat=MATRIX_FORMAT_DDM
    method%prm_AS%MatFormat=MATRIX_FORMAT_DDM

    if(method%PRECISION_KRYLOV.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       call xcpy_ddm(method%prm%n,x,method%x_hi,method%ddmMat%prm)
       call xcpy_ddm(method%prm%n,b,method%b_hi,method%ddmMat%prm)
    else
       call xcpy_ddm(method%prm%n,x,method%x,method%ddmMat%prm)
       call xcpy_ddm(method%prm%n,b,method%b,method%ddmMat%prm)
    endif

    if(method%PRECISION_x.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       allocate(method%x_lo(method%prm%n))
       call zero_initialize(method%prm%n,method%x_lo)
       if( present(x_lo) )then
          call xcpy_ddm(method%prm%n,x_lo,method%x_lo,method%ddmMat%prm)
       endif
    endif
    if(method%PRECISION_b.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       allocate(method%b_lo(method%prm%n))
       call zero_initialize(method%prm%n,method%b_lo)
       if( present(b_lo) )then
          call xcpy_ddm(method%prm%n,b_lo,method%b_lo,method%ddmMat%prm)
       endif
    endif

!    write(*,*)"solver = ",solver
!    write(*,*)" PRECISION_KRYLOV - KSP_DOUBLE_DOUBLE_PRECISION ",method%PRECISION_KRYLOV,KSP_DOUBLE_DOUBLE_PRECISION

!    mat%val_dia_lo = val_dia_lo
    if(method%PRECISION_A.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       ! do i = 1,nrow*nnd
       !    if(mat%val_dia_lo(i).ne.val_dia_lo(i))then
       !       write(*,*)i,mat%val_dia_lo(i) , val_dia_lo(i)
       !    endif
       ! enddo
    endif
    if(method%PRECISION_KRYLOV.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       select case (solver)
       case(KSP_SOLVER_CG)
          call ddKSP_pCG(method)
       case(KSP_SOLVER_BiCG)
          call ddKSP_pBiCG(method)
       case(KSP_SOLVER_GMRES)
          call ddKSP_pGMRES(method)
       case default
          if(myrank.eq.0)then
             write(*,*)"error :: not support solver = ",solver
          end if
       end select
    else
       select case (solver)
       case(KSP_SOLVER_CG)
          call KSP_pCG(method)
       case(KSP_SOLVER_BiCG)
          call KSP_pBiCG(method)
       case(KSP_SOLVER_GMRES)
          call KSP_pGMRES(method)
       case(KSP_SOLVER_CAGMRES)
          call KSP_pCAGMRES(method)
       case(KSP_SOLVER_CBCG)
          call KSP_pCBCG(method)
       case(KSP_SOLVER_POWER_METHOD)
          call power_method_symMatrix(method)
       case(KSP_SOLVER_CAARNOLDI)
          call KSP_pCAARNOLDI(method)
       case default
          if(myrank.eq.0)then
             write(*,*)"error :: not support solver = ",solver
          end if
       end select
    endif

    if( (solver.eq.KSP_SOLVER_POWER_METHOD).or.(solver.eq.KSP_SOLVER_CAARNOLDI) )then
       call get_EigenPair_count_KSP(ncount,method)
       if(allocated(method%EigenValueReal)) deallocate(method%EigenValueReal)
       if(allocated(method%EigenValueImag)) deallocate(method%EigenValueImag)
       if(allocated(method%EigenVector)) deallocate(method%EigenVector)
       if(allocated(method%EigenErrorNorm)) deallocate(method%EigenErrorNorm)

       allocate( method%EigenValueReal(ncount) )
       allocate( method%EigenValueImag(ncount) )
       allocate( method%EigenVector(ncount*nrow) )
       allocate( method%EigenErrorNorm(ncount) )

       call get_EigenPair_KSP( &           
            ncount,&                             
            method%EigenValueReal,method%EigenValueImag, &     
            method%EigenVector, &                       
            method%EigenErrorNorm, &                    
            method%MaxEigenID,method%MinEigenID,method)


       allocate( EigenType%EigenValueReal(ncount) )
       allocate( EigenType%EigenValueImag(ncount) )
       allocate( EigenType%EigenVector(ncount*nrow) )
       allocate( EigenType%EigenErrorNorm(ncount) )
       EigenType%EigenValueReal = method%EigenValueReal
       EigenType%EigenValueImag = method%EigenValueImag
       EigenType%EigenVector = method%EigenVector
       EigenType%EigenErrorNorm = method%EigenErrorNorm
       EigenType%MaxEigenID = method%MaxEigenID
       EigenType%MinEigenID = method%MinEigenID

    endif


    iret=method%iret
    n_tmp=nres

    !  call xcpy(nres,method%reshistory,res)
    call xcpy(n_tmp,method%reshistory,res)

    if(method%PRECISION_KRYLOV.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       call xcpy_ddm(method%prm%n,method%x_hi,x,method%ddmMat%prm)
    else
       call xcpy_ddm(method%prm%n,method%x,x,method%ddmMat%prm)
    endif
    if( present(x_lo) )then
       if( method%PRECISION_x .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
          call xcpy_ddm(method%prm%n,method%x_lo,x_lo,method%ddmMat%prm)
       endif
    endif
!    oloop=method%niter
    nres = method%niter
    call wrap_free_ddm(method%ddmMat)
!    if(allocated(@)) deallocate( @ )
    if(allocated(method%reshistory)) deallocate( method%reshistory )
    if(allocated(method%x)) deallocate( method%x )
    if(allocated(method%b)) deallocate( method%b )

    if(allocated(method%x_hi)) deallocate( method%x_hi )
    if(allocated(method%b_hi)) deallocate( method%b_hi )

    if(allocated(method%x_lo)) deallocate( method%x_lo )
    if(allocated(method%b_lo)) deallocate( method%b_lo )

    if( (solver.eq.KSP_SOLVER_POWER_METHOD).or.(solver.eq.KSP_SOLVER_CAARNOLDI) )then 
    else
       call parcel_ddm_free(pc,method%ddmMat)
    endif


#endif

  end subroutine parcel_ddm_main

  subroutine parcel_ddm_main_dd( &
       icomm,ndim, nrow,  &
       nnd, ioff_dia, val_dia, &
       num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
       ndiv_grids, num_grids, ioff_grids, div_direc_th, precon_thblock, &
       commhide, nas, iflag_AS,&
       b, x, &
       ! s,t,oloop, tol, bflag, qrflag, pc, &
       tol, abstol, pc, &
       nres, res, &
       Neuman_itr,&
       iflag_DDM_independ, &
       nBlock, &
       iret, &
       GMRES_m, &
       GMRES_GSflag, &
       solver, &
       PRECISION_KRYLOV, & 
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION, &
       val_dia_lo,b_lo, x_lo & 
       )
    use krylov  
    use DDM_TYPE

    use PARCEL_define
    use vec_kernel
    use vec_ddm_kernel
    use mpi
    use PARCEL_integer
    implicit none

    !arguments
    integer, intent(in) :: icomm !MPI communicator
    integer, intent(in) :: ndim !dimension of global matrix
    integer(kind=intType), intent(in) :: nrow !number of local matrix rows

    double precision, intent(in) :: b(nrow) !b vector
    double precision, intent(inout) :: x(nrow) !x vector
!    integer, intent(in) :: s ,t
!    integer, intent(inout) :: oloop
    double precision, intent(inout) :: tol
    double precision, intent(inout) :: abstol
!    integer, intent(in) :: bflag 
!    integer, intent(in) :: qrflag 
    integer, intent(in) :: pc !pre-conditioner flag
    integer, intent(inout) :: nres !number of residual history
    double precision,intent(inout) :: res(nres) !residual history

    integer, intent(in) :: nnd 
    integer, intent(in) :: ioff_dia(nnd) !offset
    double precision, intent(in) :: val_dia(nrow*nnd) !values
    integer, intent(in) :: num_neighbor_ranks
    integer, intent(in) :: neighbor_ranks( num_neighbor_ranks )
    integer, intent(in) :: margin
    integer, intent(in) :: ncomp_grids
    integer, intent(in) :: ndiv_grids(ncomp_grids)
    integer, intent(in) :: num_grids(ncomp_grids)
    integer, intent(in) :: ioff_grids(nnd,ncomp_grids)
    integer, intent(in) :: div_direc_th
    integer, intent(in) :: precon_thblock
    integer, intent(in) :: commhide ! 0: no hiding, 1: one time, 2: each time
    integer, intent(inout) :: nas
    integer, intent(in) :: iflag_AS         
    integer, intent(in) :: Neuman_itr
    integer, intent(in) :: iflag_DDM_independ
    integer, intent(in) :: nBlock
    integer, intent(inout) :: iret

    integer leng_halo
    integer AS_ioff_dia(nnd)

    integer  nprocs, myrank, ncpu, newcomm, ierr, nmax
    integer(kind=intType)  iend, i, j
    integer  as_n(4), nn_l, nn_u, &
         nsend1, nsend2, nrecv1, nrecv2, n_blk, nbuf
    double precision dummy,t_all,t_precon_init,t_cagmres_main

    integer GMRES_m         
    integer GMRES_GSflag    
    integer solver          
    integer :: PRECISION_KRYLOV
    integer :: PRECISION_A 
    integer :: PRECISION_b 
    integer :: PRECISION_x 
    integer :: PRECISION_PRECONDITION 
    double precision :: val_dia_lo( nrow*nnd) ! Matrix values low bit
    double precision :: b_lo(nrow) ! b vector low bit
    double precision :: x_lo(nrow) ! x vector low bit
    
    type(KSP) :: method
    integer(kind=intType) n_tmp
    
    integer PRECISION_A_DDM
    integer PRECISION_b_DDM
    integer PRECISION_x_DDM
    integer PRECISION_KRYLOV_DDM
    integer PRECISION_PRECONDITION_DDM ! 前処理行列の精度 (前処理の精度ではない)

#ifdef CUDA_SOLVER

#else
    call parcel_ddm_main( &
       icomm,ndim, nrow,  &
       nnd, ioff_dia, val_dia, &
       num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
       ndiv_grids, num_grids, ioff_grids, div_direc_th, precon_thblock, &
       commhide, nas, iflag_AS,&
       b, x, -1,-1, &
!       s,t,oloop, tol, bflag, qrflag, pc, 
       tol, abstol, -1,-1, pc ,&
       nres, res, &       
       Neuman_itr,&
       iflag_DDM_independ, &
       nBlock, &
       iret, &
       GMRES_m, &
       GMRES_GSflag, &
       -1, &
       -1, &
       -1, &
       -1, &
       -1, &
       -1, &
       -1, &
       -1, &
       -1, &
       -1, &
       -1, &
       solver, &
       PRECISION_KRYLOV, & 
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION, &
       val_dia_lo,b_lo, x_lo & 
       )
#endif

  end subroutine parcel_ddm_main_dd

end module ddm_kernel

module ddm_wrapper

contains
  subroutine parcel_ddcg_ddm(&
       icomm, &
       x, x_lo,&
       b, b_lo, &
       nrow, &
       nnd, ioff_dia, &
       val_dia,val_dia_lo, &
       num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
       ndiv_grids, num_grids, ioff_grids,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS,&      
       itrmax, rtolmax, reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr, &
       independ_nvec, &
       nBlock, &
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION, &
       iret &
       )
    use PARCEL_define
    use ddm_kernel
    use PARCEL_integer
    implicit none

    !arguments
    integer, intent(in) :: icomm !MPI communicator
    integer(kind=intType), intent(in) :: nrow !number of local matrix rows

    double precision, intent(inout) :: x(nrow) !x vector
    double precision, intent(inout) :: x_lo(nrow) !x vector
    double precision, intent(in) :: b(nrow) !b vector
    double precision, intent(in) :: b_lo(nrow) !b vector
!    integer, intent(inout) :: oloop
    double precision, intent(inout) :: rtolmax
    !    integer, intent(in) :: s ,t
    !    integer, intent(in) :: bflag 
    !    integer, intent(in) :: qrflag 
    integer  ipreflag !pre-conditioner flag
    integer, intent(inout) :: itrmax !number of reshistoryidual history
    double precision,intent(inout) :: reshistory(itrmax) !reshistoryidual history

    integer, intent(in) :: nnd 
    integer, intent(in) :: ioff_dia(nnd) !offset
    double precision, intent(in) :: val_dia(nrow*nnd) !values
    double precision, intent(in) :: val_dia_lo(nrow*nnd) !values
    integer, intent(in) :: num_neighbor_ranks
    integer, intent(in) :: neighbor_ranks( num_neighbor_ranks )
    integer, intent(in) :: margin
    integer, intent(in) :: ncomp_grids
    integer, intent(in) :: ndiv_grids(ncomp_grids)
    integer, intent(in) :: num_grids(ncomp_grids)
    integer, intent(in) :: ioff_grids(nnd,ncomp_grids)
    integer, intent(in) :: div_direc_th
    integer, intent(in) :: precon_thblock
    integer, intent(in) :: iovlflag ! 0: no hiding, 1: one time, 2: each time
    integer, intent(inout) :: addL
    integer, intent(in) :: iflagAS         
    integer, intent(in) :: Neuman_itr
    integer, intent(in) :: independ_nvec
    integer, intent(in) :: nBlock
    integer PRECISION_A
    integer PRECISION_b
    integer PRECISION_x
    integer PRECISION_PRECONDITION

    integer, intent(inout) :: iret

    integer solver
    double precision absrtolmax

    solver = KSP_SOLVER_CG    
    absrtolmax = 0.0d0

#ifdef CUDA_SOLVER

#else
    call parcel_ddm_main( &
         icomm,-1 , nrow, &
         nnd, ioff_dia, val_dia, &
         num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
         ndiv_grids, num_grids, ioff_grids, div_direc_th, precon_thblock, &
         iovlflag, addL, iflagAS,&
         b, x, -1,-1, &
         rtolmax, absrtolmax, -1, -1, ipreflag, &
         itrmax, reshistory, &
         Neuman_itr, &
         independ_nvec, &
         nBlock, &
         iret, &
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1, &
         -1,&
         -1,&
         -1,&
         -1, &
         solver, &
         KSP_DOUBLE_DOUBLE_PRECISION, &
         PRECISION_A, &
         PRECISION_b, &
         PRECISION_x, &
         PRECISION_PRECONDITION, &
         val_dia_lo,b_lo,x_lo &
         )
#endif

  end subroutine parcel_ddcg_ddm

  subroutine parcel_ddbicgstab_ddm( &
       icomm, &
       x, x_lo,&
       b, b_lo, &
       nrow, &
       nnd, ioff_dia, &
       val_dia,val_dia_lo, &
       num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
       ndiv_grids, num_grids, ioff_grids,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS,&      
       itrmax, rtolmax, reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr, &
       independ_nvec, &
       nBlock, &
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION, &
       iret &
       )
    use PARCEL_define
    use ddm_kernel
    use PARCEL_integer
    implicit none

    !arguments
    integer, intent(in) :: icomm !MPI communicator
    integer(kind=intType), intent(in) :: nrow !number of local matrix rows

    double precision, intent(inout) :: x(nrow) !x vector
    double precision, intent(inout) :: x_lo(nrow) !x vector
    double precision, intent(in) :: b(nrow) !b vector
    double precision, intent(in) :: b_lo(nrow) !b vector
!    integer, intent(inout) :: oloop
    double precision, intent(inout) :: rtolmax
    !    integer, intent(in) :: s ,t
    !    integer, intent(in) :: bflag 
    !    integer, intent(in) :: qrflag 
    integer  ipreflag !pre-conditioner flag
    integer, intent(inout) :: itrmax !number of reshistoryidual history
    double precision,intent(inout) :: reshistory(itrmax) !reshistoryidual history

    integer, intent(in) :: nnd 
    integer, intent(in) :: ioff_dia(nnd) !offset
    double precision, intent(in) :: val_dia(nrow*nnd) !values
    double precision, intent(in) :: val_dia_lo(nrow*nnd) !values
    integer, intent(in) :: num_neighbor_ranks
    integer, intent(in) :: neighbor_ranks( num_neighbor_ranks )
    integer, intent(in) :: margin
    integer, intent(in) :: ncomp_grids
    integer, intent(in) :: ndiv_grids(ncomp_grids)
    integer, intent(in) :: num_grids(ncomp_grids)
    integer, intent(in) :: ioff_grids(nnd,ncomp_grids)
    integer, intent(in) :: div_direc_th
    integer, intent(in) :: precon_thblock
    integer, intent(in) :: iovlflag ! 0: no hiding, 1: one time, 2: each time
    integer, intent(inout) :: addL
    integer, intent(in) :: iflagAS         
    integer, intent(in) :: Neuman_itr
    integer, intent(in) :: independ_nvec
    integer, intent(in) :: nBlock
    integer, intent(inout) :: iret
    integer PRECISION_A
    integer PRECISION_b
    integer PRECISION_x
    integer PRECISION_PRECONDITION

    integer solver
    double precision absrtolmax

    solver = KSP_SOLVER_BiCG 
    absrtolmax = 0.0d0

#ifdef CUDA_SOLVER

#else
    call parcel_ddm_main( &
         icomm,-1 , nrow, &
         nnd, ioff_dia, val_dia, &
         num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
         ndiv_grids, num_grids, ioff_grids, div_direc_th, precon_thblock, &
         iovlflag, addL, iflagAS,&
         b, x, -1,-1, &
         rtolmax, absrtolmax, -1, -1, ipreflag, &
         itrmax, reshistory, &
       Neuman_itr, &
         independ_nvec, &
         nBlock, &
         iret, &
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1, &
         -1,&
         -1,&
         -1,&
         -1, &
         solver, &
         KSP_DOUBLE_DOUBLE_PRECISION, &
         PRECISION_A, &
         PRECISION_b, &
         PRECISION_x, &
         PRECISION_PRECONDITION, &
         val_dia_lo,b_lo,x_lo &
         )
#endif

  end subroutine parcel_ddbicgstab_ddm

  subroutine parcel_ddgmres_ddm( &
       icomm, &
       x, x_lo,&
       b, b_lo, &
       nrow, &
       nnd, ioff_dia, &
       val_dia,val_dia_lo, &
       num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
       ndiv_grids, num_grids, ioff_grids,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS,&      
       itrmax, rtolmax, reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr, &
       independ_nvec, &
       nBlock, &
       GMRES_m, &
       GMRES_GSflag, &
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION, &
       iret &
       )
    use PARCEL_define
    use ddm_kernel
    use PARCEL_integer
    implicit none

    !arguments
    integer, intent(in) :: icomm !MPI communicator
    integer(kind=intType), intent(in) :: nrow !number of local matrix rows

    double precision, intent(inout) :: x(nrow) !x vector
    double precision, intent(inout) :: x_lo(nrow) !x vector
    double precision, intent(in) :: b(nrow) !b vector
    double precision, intent(in) :: b_lo(nrow) !b vector
!    integer, intent(inout) :: oloop
    double precision, intent(inout) :: rtolmax
    !    integer, intent(in) :: s ,t
    !    integer, intent(in) :: bflag 
    !    integer, intent(in) :: qrflag 
    integer  ipreflag !pre-conditioner flag
    integer, intent(inout) :: itrmax !number of reshistoryidual history
    double precision,intent(inout) :: reshistory(itrmax) !reshistoryidual history

    integer, intent(in) :: nnd 
    integer, intent(in) :: ioff_dia(nnd) !offset
    double precision, intent(in) :: val_dia(nrow*nnd) !values
    double precision, intent(in) :: val_dia_lo(nrow*nnd) !values
    integer, intent(in) :: num_neighbor_ranks
    integer, intent(in) :: neighbor_ranks( num_neighbor_ranks )
    integer, intent(in) :: margin
    integer, intent(in) :: ncomp_grids
    integer, intent(in) :: ndiv_grids(ncomp_grids)
    integer, intent(in) :: num_grids(ncomp_grids)
    integer, intent(in) :: ioff_grids(nnd,ncomp_grids)
    integer, intent(in) :: div_direc_th
    integer, intent(in) :: precon_thblock
    integer, intent(in) :: iovlflag ! 0: no hiding, 1: one time, 2: each time
    integer, intent(inout) :: addL
    integer, intent(in) :: iflagAS         
    integer, intent(in) :: Neuman_itr
    integer, intent(in) :: independ_nvec
    integer, intent(in) :: nBlock
    integer PRECISION_A
    integer PRECISION_b
    integer PRECISION_x
    integer PRECISION_PRECONDITION
    integer, intent(inout) :: iret
    integer GMRES_m         
    integer GMRES_GSflag    

    integer solver
    double precision absrtolmax

    absrtolmax = 0.0d0
    solver = KSP_SOLVER_GMRES

#ifdef CUDA_SOLVER

#else
    call parcel_ddm_main( &
         icomm,-1 , nrow, &
         nnd, ioff_dia, val_dia, &
         num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
         ndiv_grids, num_grids, ioff_grids, div_direc_th, precon_thblock, &
         iovlflag, addL, iflagAS,&
         b, x, -1,-1, &
         rtolmax, absrtolmax, -1, -1, ipreflag, &
         itrmax, reshistory, &
         Neuman_itr, &
         independ_nvec, &
         nBlock, &
         iret, &
         GMRES_m, &
         GMRES_GSflag, &
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1, &
         -1,&
         -1,&
         -1,&
         -1, &
         solver, &
         KSP_DOUBLE_DOUBLE_PRECISION, &
         PRECISION_A, &
         PRECISION_b, &
         PRECISION_x, &
         PRECISION_PRECONDITION, &
         val_dia_lo,b_lo,x_lo &
         )
#endif

  end subroutine parcel_ddgmres_ddm

#ifdef CUDA_SOLVER
  subroutine parcel_dcg_ddm(&
       icomm, &
       x,b, &
       nrow, &
       nnd, ioff_dia, val_dia, &
       num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
       ndiv_grids, num_grids, ioff_grids,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS,&      
       itrmax, rtolmax, reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr,&
       independ_nvec, &
       nBlock, &
       iret &
       )
    use PARCEL_define
    use ddm_kernel
    use PARCEL_integer
    use iso_c_binding
    implicit none

    !arguments
    integer, intent(in) :: icomm !MPI communicator
    integer(kind=intType), intent(in) :: nrow !number of local matrix rows

    type(c_ptr) x !x vector
    type(c_ptr) b !b vector
    double precision, intent(inout) :: rtolmax
    !    integer, intent(in) :: s ,t
    !    integer, intent(in) :: bflag 
    !    integer, intent(in) :: qrflag 
    integer  ipreflag !pre-conditioner flag
    integer, intent(inout) :: itrmax !number of reshistoryidual history
    double precision,intent(inout) :: reshistory(itrmax) !reshistoryidual history

    integer, intent(in) :: nnd 
    type(c_ptr) ioff_dia !offset
    type(c_ptr) val_dia  !values

    integer, intent(in) :: num_neighbor_ranks
    integer, intent(in) :: neighbor_ranks( num_neighbor_ranks )
    integer, intent(in) :: margin
    integer, intent(in) :: ncomp_grids
    integer, intent(in) :: ndiv_grids(ncomp_grids)
    integer, intent(in) :: num_grids(ncomp_grids)
    integer, intent(in) :: ioff_grids(nnd,ncomp_grids)
    integer, intent(in) :: div_direc_th
    integer, intent(in) :: precon_thblock
    integer, intent(in) :: iovlflag ! 0: no hiding, 1: one time, 2: each time
    integer, intent(inout) :: addL
    integer, intent(in) :: iflagAS         
    integer, intent(in) :: Neuman_itr
    integer, intent(in) :: independ_nvec
    integer, intent(in) :: nBlock
    integer, intent(inout) :: iret

    integer solver
    double precision absrtolmax

    solver = KSP_SOLVER_CG    
    absrtolmax = 0.0d0

    call parcel_ddm_main( &
         icomm,-1 , nrow, &
         nnd, ioff_dia, val_dia, &
         num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
         ndiv_grids, num_grids, ioff_grids, div_direc_th, precon_thblock, &
         iovlflag, addL, iflagAS,&
         b, x, -1,-1, &
         rtolmax, absrtolmax, -1, -1, ipreflag, &
         itrmax, reshistory, &
         Neuman_itr,&
         independ_nvec, &
         nBlock, &
         iret, &
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1, &
         -1,&
         -1,&
         -1,&
         -1, &
         solver &
         )

  end subroutine parcel_dcg_ddm

  subroutine parcel_dbicgstab_ddm( &
       icomm, &
       x,b, &
       nrow, &
       nnd, ioff_dia, val_dia, &
       num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
       ndiv_grids, num_grids, ioff_grids,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS,&      
       itrmax, rtolmax, reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr,&
       independ_nvec, &
       nBlock, &
       iret &
       )
    use PARCEL_define
    use ddm_kernel
    use PARCEL_integer
    use iso_c_binding
    implicit none

    !arguments
    integer, intent(in) :: icomm !MPI communicator
    integer(kind=intType), intent(in) :: nrow !number of local matrix rows

    type(c_ptr) x !x vector
    type(c_ptr) b !b vector
    double precision, intent(inout) :: rtolmax
    !    integer, intent(in) :: s ,t
    !    integer, intent(in) :: bflag 
    !    integer, intent(in) :: qrflag 
    integer  ipreflag !pre-conditioner flag
    integer, intent(inout) :: itrmax !number of reshistoryidual history
    double precision,intent(inout) :: reshistory(itrmax) !reshistoryidual history

    integer, intent(in) :: nnd 
    type(c_ptr) ioff_dia !offset
    type(c_ptr) val_dia  !values
    integer, intent(in) :: num_neighbor_ranks
    integer, intent(in) :: neighbor_ranks( num_neighbor_ranks )
    integer, intent(in) :: margin
    integer, intent(in) :: ncomp_grids
    integer, intent(in) :: ndiv_grids(ncomp_grids)
    integer, intent(in) :: num_grids(ncomp_grids)
    integer, intent(in) :: ioff_grids(nnd,ncomp_grids)
    integer, intent(in) :: div_direc_th
    integer, intent(in) :: precon_thblock
    integer, intent(in) :: iovlflag ! 0: no hiding, 1: one time, 2: each time
    integer, intent(inout) :: addL
    integer, intent(in) :: iflagAS         
    integer, intent(in) :: Neuman_itr
    integer, intent(in) :: independ_nvec
    integer, intent(in) :: nBlock
    integer, intent(inout) :: iret

    integer solver
    double precision absrtolmax

    solver = KSP_SOLVER_BiCG 
    absrtolmax = 0.0d0

    call parcel_ddm_main( &
         icomm,-1 , nrow, &
         nnd, ioff_dia, val_dia, &
         num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
         ndiv_grids, num_grids, ioff_grids, div_direc_th, precon_thblock, &
         iovlflag, addL, iflagAS,&
         b, x, -1,-1, &
         rtolmax, absrtolmax, -1, -1, ipreflag, &
         itrmax, reshistory, &
         Neuman_itr,&
         independ_nvec, &
         nBlock, &
         iret, &
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1, &
         -1,&
         -1,&
         -1,&
         -1, &
         solver &
         )
  end subroutine parcel_dbicgstab_ddm

  subroutine parcel_dgmres_ddm( &
       icomm, &
       x,b, &
       nrow, &
       nnd, ioff_dia, val_dia, &
       num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
       ndiv_grids, num_grids, ioff_grids,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS,&      
       itrmax, rtolmax, reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr,&
       independ_nvec, &
       nBlock, &
       GMRES_m, &
       GMRES_GSflag, &
       iret &
       )
    use PARCEL_define
    use ddm_kernel
    use PARCEL_integer
    use iso_c_binding
    implicit none

    !arguments
    integer, intent(in) :: icomm !MPI communicator
    integer(kind=intType), intent(in) :: nrow !number of local matrix rows

    type(c_ptr) x !x vector
    type(c_ptr) b !b vector
    double precision, intent(inout) :: rtolmax
    !    integer, intent(in) :: s ,t
    !    integer, intent(in) :: bflag 
    !    integer, intent(in) :: qrflag 
    integer  ipreflag !pre-conditioner flag
    integer, intent(inout) :: itrmax !number of reshistoryidual history
    double precision,intent(inout) :: reshistory(itrmax) !reshistoryidual history

    integer, intent(in) :: nnd 
    type(c_ptr) ioff_dia !offset
    type(c_ptr) val_dia  !values
    integer, intent(in) :: num_neighbor_ranks
    integer, intent(in) :: neighbor_ranks( num_neighbor_ranks )
    integer, intent(in) :: margin
    integer, intent(in) :: ncomp_grids
    integer, intent(in) :: ndiv_grids(ncomp_grids)
    integer, intent(in) :: num_grids(ncomp_grids)
    integer, intent(in) :: ioff_grids(nnd,ncomp_grids)
    integer, intent(in) :: div_direc_th
    integer, intent(in) :: precon_thblock
    integer, intent(in) :: iovlflag ! 0: no hiding, 1: one time, 2: each time
    integer, intent(inout) :: addL
    integer, intent(in) :: iflagAS         
    integer, intent(in) :: Neuman_itr
    integer, intent(in) :: independ_nvec
    integer, intent(in) :: nBlock
    integer, intent(inout) :: iret
    integer GMRES_m         
    integer GMRES_GSflag    

    integer solver
    double precision absrtolmax

    absrtolmax = 0.0d0
    solver = KSP_SOLVER_GMRES

    call parcel_ddm_main( &
         icomm,-1 , nrow, &
         nnd, ioff_dia, val_dia, &
         num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
         ndiv_grids, num_grids, ioff_grids, div_direc_th, precon_thblock, &
         iovlflag, addL, iflagAS,&
         b, x, -1,-1, &
         rtolmax, absrtolmax, -1, -1, ipreflag, &
         itrmax, reshistory, &
         Neuman_itr,&
         independ_nvec, &
         nBlock, &
         iret, &
         GMRES_m, &
         GMRES_GSflag, &
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1, &
         -1,&
         -1,&
         -1,&
         -1, &
         solver &
         )
  end subroutine parcel_dgmres_ddm

   subroutine parcel_dcagmres_ddm(&
       icomm, &
       x,b, &
       nrow, &
       nnd, ioff_dia, val_dia, &
       num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
       ndiv_grids, num_grids, ioff_grids,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS,&      
       itrmax, rtolmax, reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr,&
       independ_nvec, &
       nBlock, &
       CAGMRES_sstep, &
       CAGMRES_tstep, &
       CAGMRES_basis, &
       CAGMRES_QRflag, &
       iret &
       )
    use PARCEL_define
    use ddm_kernel
    use PARCEL_integer
    use iso_c_binding
    implicit none

    !arguments
    integer, intent(in) :: icomm !MPI communicator
    integer(kind=intType), intent(in) :: nrow !number of local matrix rows

    type(c_ptr) x !x vector
    type(c_ptr) b !b vector
    double precision, intent(inout) :: rtolmax
    !    integer, intent(in) :: s ,t
    !    integer, intent(in) :: bflag 
    !    integer, intent(in) :: qrflag 
    integer  ipreflag !pre-conditioner flag
    integer, intent(inout) :: itrmax !number of reshistoryidual history
    double precision,intent(inout) :: reshistory(itrmax) !reshistoryidual history

    integer, intent(in) :: nnd 
    type(c_ptr) ioff_dia !offset
    type(c_ptr) val_dia  !values
    integer, intent(in) :: num_neighbor_ranks
    integer, intent(in) :: neighbor_ranks( num_neighbor_ranks )
    integer, intent(in) :: margin
    integer, intent(in) :: ncomp_grids
    integer, intent(in) :: ndiv_grids(ncomp_grids)
    integer, intent(in) :: num_grids(ncomp_grids)
    integer, intent(in) :: ioff_grids(nnd,ncomp_grids)
    integer, intent(in) :: div_direc_th
    integer, intent(in) :: precon_thblock
    integer, intent(in) :: iovlflag ! 0: no hiding, 1: one time, 2: each time
    integer, intent(inout) :: addL
    integer, intent(in) :: iflagAS         
    integer, intent(in) :: Neuman_itr
    integer, intent(in) :: independ_nvec
    integer, intent(in) :: nBlock
    integer, intent(inout) :: iret
    integer CAGMRES_sstep
    integer CAGMRES_tstep
    integer CAGMRES_basis
    integer CAGMRES_QRflag
    integer solver
    double precision absrtolmax

    absrtolmax = 0.0d0
    solver = KSP_SOLVER_CAGMRES

    call parcel_ddm_main( &
         icomm,-1 , nrow, &
         nnd, ioff_dia, val_dia, &
         num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
         ndiv_grids, num_grids, ioff_grids, div_direc_th, precon_thblock, &
         iovlflag, addL, iflagAS,&
         b, x, -1,-1, &
         rtolmax, absrtolmax, -1, -1, ipreflag, &
         itrmax, reshistory, &
         Neuman_itr,&
         independ_nvec, &
         nBlock, &
         iret, &
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         CAGMRES_sstep, &
         CAGMRES_tstep, &
         CAGMRES_basis, &
         CAGMRES_QRflag, &
         solver &
         )
  end subroutine parcel_dcagmres_ddm

  subroutine parcel_dcaarnoldi_ddm(&
       icomm, &
       x,b, &
       nrow, &
       nnd, ioff_dia, val_dia, &
       num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
       ndiv_grids, num_grids, ioff_grids,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS,&      
       iovlflag, &
       precon_thblock, &
       Neuman_itr,&
       independ_nvec, &
       nBlock, &
       CAARNOLDI_sstep, &
       CAARNOLDI_tstep, &
       CAARNOLDI_basis, &
       CAARNOLDI_QRflag, &
       EigenValueReal_DDM, &
       EigenValueImag_DDM, &
       EigenVector_DDM, &
       EigenErrorNorm_DDM, &
       EmaxID_DDM, &
       EminID_DDM, &
       iret &
       )
    use PARCEL_define
    use init_ddm_param, only:&
         parcel_ddm_free
    use EigenSolver
    use ddm_kernel
    use PARCEL_integer
    use iso_c_binding
    implicit none

    !arguments
    integer, intent(in) :: icomm !MPI communicator
    integer(kind=intType), intent(in) :: nrow !number of local matrix rows

    type(c_ptr) x !x vector
    type(c_ptr) b !b vector
    !    integer, intent(in) :: s ,t
    !    integer, intent(in) :: bflag 
    !    integer, intent(in) :: qrflag 
    integer  ipreflag !pre-conditioner flag
    integer, intent(in) :: nnd 
    type(c_ptr) ioff_dia !offset
    type(c_ptr) val_dia  !values
    integer, intent(in) :: num_neighbor_ranks
    integer, intent(in) :: neighbor_ranks( num_neighbor_ranks )
    integer, intent(in) :: margin
    integer, intent(in) :: ncomp_grids
    integer, intent(in) :: ndiv_grids(ncomp_grids)
    integer, intent(in) :: num_grids(ncomp_grids)
    integer, intent(in) :: ioff_grids(nnd,ncomp_grids)
    integer, intent(in) :: div_direc_th
    integer, intent(in) :: precon_thblock
    integer, intent(in) :: iovlflag ! 0: no hiding, 1: one time, 2: each time
    integer, intent(inout) :: addL
    integer, intent(in) :: iflagAS         
    integer, intent(in) :: Neuman_itr
    integer, intent(in) :: independ_nvec
    integer, intent(in) :: nBlock
    integer, intent(inout) :: iret

    integer CAARNOLDI_sstep
    integer CAARNOLDI_tstep
    integer CAARNOLDI_basis
    integer CAARNOLDI_QRflag
    double precision EigenValueReal_DDM(CAARNOLDI_sstep*CAARNOLDI_tstep)
    double precision EigenValueImag_DDM(CAARNOLDI_sstep*CAARNOLDI_tstep)
    double precision EigenVector_DDM(CAARNOLDI_sstep*CAARNOLDI_tstep*nrow)
    double precision EigenErrorNorm_DDM(CAARNOLDI_sstep*CAARNOLDI_tstep)
    integer EmaxID_DDM
    integer EminID_DDM

    integer solver
    double precision absrtolmax
    double precision  rtolmax
    double precision reshistory(CAARNOLDI_sstep*CAARNOLDI_tstep+1) 
    integer itrmax
    type(Eigen_Type) :: EigenType_

    absrtolmax = 0.0d0
    rtolmax = 0.0d0
    reshistory = 0.0d0
    itrmax = CAARNOLDI_sstep*CAARNOLDI_tstep
    solver = KSP_SOLVER_CAARNOLDI

    call parcel_ddm_main( &
         icomm,-1 , nrow, &
         nnd, ioff_dia, val_dia, &
         num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
         ndiv_grids, num_grids, ioff_grids, div_direc_th, precon_thblock, &
         iovlflag, addL, iflagAS,&
         b, x, -1,-1, &
         rtolmax, absrtolmax, -1, -1, ipreflag, &
         itrmax, reshistory, &
         Neuman_itr,&
         independ_nvec, &
         nBlock, &
         iret, &
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         CAARNOLDI_sstep, &
         CAARNOLDI_tstep, &
         CAARNOLDI_basis, &
         CAARNOLDI_QRflag, &
         -1, &
         -1, &
         -1, &
         -1, &
         solver, &
         EigenType=EigenType_ &
         )

    EigenValueReal_DDM(1:CAARNOLDI_sstep*CAARNOLDI_tstep)   = EigenType_%EigenValueReal(1:CAARNOLDI_sstep*CAARNOLDI_tstep)
    EigenValueImag_DDM(1:CAARNOLDI_sstep*CAARNOLDI_tstep)   = EigenType_%EigenValueImag(1:CAARNOLDI_sstep*CAARNOLDI_tstep)
    EigenErrorNorm_DDM(1:CAARNOLDI_sstep*CAARNOLDI_tstep)   = EigenType_%EigenErrorNorm(1:CAARNOLDI_sstep*CAARNOLDI_tstep)
    EigenVector_DDM(1:CAARNOLDI_sstep*CAARNOLDI_tstep*nrow) = EigenType_%EigenVector(1:CAARNOLDI_sstep*CAARNOLDI_tstep*nrow)
    EmaxID_DDM = EigenType_%MaxEigenID 
    EminID_DDM = EigenType_%MinEigenID

  end subroutine parcel_dcaarnoldi_ddm

  subroutine parcel_dcbcg_ddm(&
       icomm, &
       x,b, &
       nrow, &
       nnd, ioff_dia, val_dia, &
       num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
       ndiv_grids, num_grids, ioff_grids,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS,&      
       itrmax, rtolmax, reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr,&
       independ_nvec, &
       nBlock, &
       CBCG_kstep, &
       CBCG_Eigenflag, &
       power_method_itrmax, &
       CAARNOLDI_sstep, &
       CAARNOLDI_tstep, &
       CAARNOLDI_basis, &
       CAARNOLDI_QRflag, &
       iret &
       )
    use PARCEL_define
    use ddm_kernel
    use PARCEL_integer
    use iso_c_binding
    implicit none

    !arguments
    integer, intent(in) :: icomm !MPI communicator
    integer(kind=intType), intent(in) :: nrow !number of local matrix rows

    type(c_ptr) x !x vector
    type(c_ptr) b !b vector
!    integer, intent(inout) :: oloop
    double precision, intent(inout) :: rtolmax
    !    integer, intent(in) :: s ,t
    !    integer, intent(in) :: bflag 
    !    integer, intent(in) :: qrflag 
    integer  ipreflag !pre-conditioner flag
    integer, intent(inout) :: itrmax !number of reshistoryidual history
    double precision,intent(inout) :: reshistory(itrmax) !reshistoryidual history

    integer, intent(in) :: nnd 
    type(c_ptr) ioff_dia !offset
    type(c_ptr) val_dia  !values

    integer, intent(in) :: num_neighbor_ranks
    integer, intent(in) :: neighbor_ranks( num_neighbor_ranks )
    integer, intent(in) :: margin
    integer, intent(in) :: ncomp_grids
    integer, intent(in) :: ndiv_grids(ncomp_grids)
    integer, intent(in) :: num_grids(ncomp_grids)
    integer, intent(in) :: ioff_grids(nnd,ncomp_grids)
    integer, intent(in) :: div_direc_th
    integer, intent(in) :: precon_thblock
    integer, intent(in) :: iovlflag ! 0: no hiding, 1: one time, 2: each time
    integer, intent(inout) :: addL
    integer, intent(in) :: iflagAS         
    integer, intent(in) :: Neuman_itr
    integer, intent(in) :: independ_nvec
    integer, intent(in) :: nBlock
    integer, intent(inout) :: iret
    integer CBCG_kstep
    integer CBCG_Eigenflag
    integer power_method_itrmax
    integer CAARNOLDI_sstep
    integer CAARNOLDI_tstep
    integer CAARNOLDI_basis
    integer CAARNOLDI_QRflag
    integer solver
    double precision absrtolmax

    absrtolmax = 0.0d0
    solver = KSP_SOLVER_CBCG

    call parcel_ddm_main( &
         icomm,-1 , nrow, &
         nnd, ioff_dia, val_dia, &
         num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
         ndiv_grids, num_grids, ioff_grids, div_direc_th, precon_thblock, &
         iovlflag, addL, iflagAS,&
         b, x, -1,-1, &
         rtolmax, absrtolmax, -1, -1, ipreflag, &
         itrmax, reshistory, &
       Neuman_itr,&
         independ_nvec, &
         nBlock, &
         iret, &
         -1,&
         -1,&
         CBCG_kstep, &
         CBCG_Eigenflag, &
         power_method_itrmax, &
         CAARNOLDI_sstep, &
         CAARNOLDI_tstep, &
         CAARNOLDI_basis, &
         CAARNOLDI_QRflag, &
         -1, &
         -1, &
         -1, &
         -1, &
         solver &
         )
  end subroutine parcel_dcbcg_ddm
#else
  subroutine parcel_dcg_ddm(&
       icomm, &
       x,b, &
       nrow, &
       nnd, ioff_dia, val_dia, &
       num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
       ndiv_grids, num_grids, ioff_grids,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS,&      
       itrmax, rtolmax, reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr,&
       independ_nvec, &
       nBlock, &
       iret &
       )
    use PARCEL_define
    use ddm_kernel
    use PARCEL_integer
    implicit none

    !arguments
    integer, intent(in) :: icomm !MPI communicator
    integer(kind=intType), intent(in) :: nrow !number of local matrix rows

    double precision, intent(in) :: b(nrow) !b vector
    double precision, intent(inout) :: x(nrow) !x vector
!    integer, intent(inout) :: oloop
    double precision, intent(inout) :: rtolmax
    !    integer, intent(in) :: s ,t
    !    integer, intent(in) :: bflag 
    !    integer, intent(in) :: qrflag 
    integer  ipreflag !pre-conditioner flag
    integer, intent(inout) :: itrmax !number of reshistoryidual history
    double precision,intent(inout) :: reshistory(itrmax) !reshistoryidual history

    integer, intent(in) :: nnd 
    integer, intent(in) :: ioff_dia(nnd) !offset
    double precision, intent(in) :: val_dia(nrow*nnd) !values
    integer, intent(in) :: num_neighbor_ranks
    integer, intent(in) :: neighbor_ranks( num_neighbor_ranks )
    integer, intent(in) :: margin
    integer, intent(in) :: ncomp_grids
    integer, intent(in) :: ndiv_grids(ncomp_grids)
    integer, intent(in) :: num_grids(ncomp_grids)
    integer, intent(in) :: ioff_grids(nnd,ncomp_grids)
    integer, intent(in) :: div_direc_th
    integer, intent(in) :: precon_thblock
    integer, intent(in) :: iovlflag ! 0: no hiding, 1: one time, 2: each time
    integer, intent(inout) :: addL
    integer, intent(in) :: iflagAS         
    integer, intent(in) :: Neuman_itr
    integer, intent(in) :: independ_nvec
    integer, intent(in) :: nBlock
    integer, intent(inout) :: iret

    integer solver
    double precision absrtolmax

    solver = KSP_SOLVER_CG    
    absrtolmax = 0.0d0

    call parcel_ddm_main( &
         icomm,-1 , nrow, &
         nnd, ioff_dia, val_dia, &
         num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
         ndiv_grids, num_grids, ioff_grids, div_direc_th, precon_thblock, &
         iovlflag, addL, iflagAS,&
         b, x, -1,-1, &
         rtolmax, absrtolmax, -1, -1, ipreflag, &
         itrmax, reshistory, &
         Neuman_itr,&
         independ_nvec, &
         nBlock, &
         iret, &
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1, &
         -1,&
         -1,&
         -1,&
         -1, &
         solver &
         )

  end subroutine parcel_dcg_ddm

  subroutine parcel_dbicgstab_ddm( &
       icomm, &
       x,b, &
       nrow, &
       nnd, ioff_dia, val_dia, &
       num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
       ndiv_grids, num_grids, ioff_grids,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS,&      
       itrmax, rtolmax, reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr,&
       independ_nvec, &
       nBlock, &
       iret &
       )
    use PARCEL_define
    use ddm_kernel
    use PARCEL_integer
    implicit none

    !arguments
    integer, intent(in) :: icomm !MPI communicator
    integer(kind=intType), intent(in) :: nrow !number of local matrix rows

    double precision, intent(in) :: b(nrow) !b vector
    double precision, intent(inout) :: x(nrow) !x vector
!    integer, intent(inout) :: oloop
    double precision, intent(inout) :: rtolmax
    !    integer, intent(in) :: s ,t
    !    integer, intent(in) :: bflag 
    !    integer, intent(in) :: qrflag 
    integer  ipreflag !pre-conditioner flag
    integer, intent(inout) :: itrmax !number of reshistoryidual history
    double precision,intent(inout) :: reshistory(itrmax) !reshistoryidual history

    integer, intent(in) :: nnd 
    integer, intent(in) :: ioff_dia(nnd) !offset
    double precision, intent(in) :: val_dia(nrow*nnd) !values
    integer, intent(in) :: num_neighbor_ranks
    integer, intent(in) :: neighbor_ranks( num_neighbor_ranks )
    integer, intent(in) :: margin
    integer, intent(in) :: ncomp_grids
    integer, intent(in) :: ndiv_grids(ncomp_grids)
    integer, intent(in) :: num_grids(ncomp_grids)
    integer, intent(in) :: ioff_grids(nnd,ncomp_grids)
    integer, intent(in) :: div_direc_th
    integer, intent(in) :: precon_thblock
    integer, intent(in) :: iovlflag ! 0: no hiding, 1: one time, 2: each time
    integer, intent(inout) :: addL
    integer, intent(in) :: iflagAS         
    integer, intent(in) :: Neuman_itr
    integer, intent(in) :: independ_nvec
    integer, intent(in) :: nBlock
    integer, intent(inout) :: iret

    integer solver
    double precision absrtolmax

    solver = KSP_SOLVER_BiCG 
    absrtolmax = 0.0d0

    call parcel_ddm_main( &
         icomm,-1 , nrow, &
         nnd, ioff_dia, val_dia, &
         num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
         ndiv_grids, num_grids, ioff_grids, div_direc_th, precon_thblock, &
         iovlflag, addL, iflagAS,&
         b, x, -1,-1, &
         rtolmax, absrtolmax, -1, -1, ipreflag, &
         itrmax, reshistory, &
         Neuman_itr,&
         independ_nvec, &
         nBlock, &
         iret, &
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1, &
         -1,&
         -1,&
         -1,&
         -1, &
         solver &
         )
  end subroutine parcel_dbicgstab_ddm

  subroutine parcel_dgmres_ddm( &
       icomm, &
       x,b, &
       nrow, &
       nnd, ioff_dia, val_dia, &
       num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
       ndiv_grids, num_grids, ioff_grids,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS,&      
       itrmax, rtolmax, reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr,&
       independ_nvec, &
       nBlock, &
       GMRES_m, &
       GMRES_GSflag, &
       iret &
       )
    use PARCEL_define
    use ddm_kernel
    use PARCEL_integer
    implicit none

    !arguments
    integer, intent(in) :: icomm !MPI communicator
    integer(kind=intType), intent(in) :: nrow !number of local matrix rows

    double precision, intent(in) :: b(nrow) !b vector
    double precision, intent(inout) :: x(nrow) !x vector
!    integer, intent(inout) :: oloop
    double precision, intent(inout) :: rtolmax
    !    integer, intent(in) :: s ,t
    !    integer, intent(in) :: bflag 
    !    integer, intent(in) :: qrflag 
    integer  ipreflag !pre-conditioner flag
    integer, intent(inout) :: itrmax !number of reshistoryidual history
    double precision,intent(inout) :: reshistory(itrmax) !reshistoryidual history

    integer, intent(in) :: nnd 
    integer, intent(in) :: ioff_dia(nnd) !offset
    double precision, intent(in) :: val_dia(nrow*nnd) !values
    integer, intent(in) :: num_neighbor_ranks
    integer, intent(in) :: neighbor_ranks( num_neighbor_ranks )
    integer, intent(in) :: margin
    integer, intent(in) :: ncomp_grids
    integer, intent(in) :: ndiv_grids(ncomp_grids)
    integer, intent(in) :: num_grids(ncomp_grids)
    integer, intent(in) :: ioff_grids(nnd,ncomp_grids)
    integer, intent(in) :: div_direc_th
    integer, intent(in) :: precon_thblock
    integer, intent(in) :: iovlflag ! 0: no hiding, 1: one time, 2: each time
    integer, intent(inout) :: addL
    integer, intent(in) :: iflagAS         
    integer, intent(in) :: Neuman_itr
    integer, intent(in) :: independ_nvec
    integer, intent(in) :: nBlock
    integer, intent(inout) :: iret
    integer GMRES_m         
    integer GMRES_GSflag    

    integer solver
    double precision absrtolmax

    absrtolmax = 0.0d0
    solver = KSP_SOLVER_GMRES

    call parcel_ddm_main( &
         icomm,-1 , nrow, &
         nnd, ioff_dia, val_dia, &
         num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
         ndiv_grids, num_grids, ioff_grids, div_direc_th, precon_thblock, &
         iovlflag, addL, iflagAS,&
         b, x, -1,-1, &
         rtolmax, absrtolmax, -1, -1, ipreflag, &
         itrmax, reshistory, &
         Neuman_itr,&
         independ_nvec, &
         nBlock, &
         iret, &
         GMRES_m, &
         GMRES_GSflag, &
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1, &
         -1,&
         -1,&
         -1,&
         -1, &
         solver &
         )
  end subroutine parcel_dgmres_ddm

   subroutine parcel_dcagmres_ddm(&
       icomm, &
       x,b, &
       nrow, &
       nnd, ioff_dia, val_dia, &
       num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
       ndiv_grids, num_grids, ioff_grids,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS,&      
       itrmax, rtolmax, reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr,&
       independ_nvec, &
       nBlock, &
       CAGMRES_sstep, &
       CAGMRES_tstep, &
       CAGMRES_basis, &
       CAGMRES_QRflag, &
       iret &
       )
    use PARCEL_define
    use ddm_kernel
    use PARCEL_integer
    implicit none

    !arguments
    integer, intent(in) :: icomm !MPI communicator
    integer(kind=intType), intent(in) :: nrow !number of local matrix rows

    double precision, intent(in) :: b(nrow) !b vector
    double precision, intent(inout) :: x(nrow) !x vector
!    integer, intent(inout) :: oloop
    double precision, intent(inout) :: rtolmax
    !    integer, intent(in) :: s ,t
    !    integer, intent(in) :: bflag 
    !    integer, intent(in) :: qrflag 
    integer  ipreflag !pre-conditioner flag
    integer, intent(inout) :: itrmax !number of reshistoryidual history
    double precision,intent(inout) :: reshistory(itrmax) !reshistoryidual history

    integer, intent(in) :: nnd 
    integer, intent(in) :: ioff_dia(nnd) !offset
    double precision, intent(in) :: val_dia(nrow*nnd) !values
    integer, intent(in) :: num_neighbor_ranks
    integer, intent(in) :: neighbor_ranks( num_neighbor_ranks )
    integer, intent(in) :: margin
    integer, intent(in) :: ncomp_grids
    integer, intent(in) :: ndiv_grids(ncomp_grids)
    integer, intent(in) :: num_grids(ncomp_grids)
    integer, intent(in) :: ioff_grids(nnd,ncomp_grids)
    integer, intent(in) :: div_direc_th
    integer, intent(in) :: precon_thblock
    integer, intent(in) :: iovlflag ! 0: no hiding, 1: one time, 2: each time
    integer, intent(inout) :: addL
    integer, intent(in) :: iflagAS         
    integer, intent(in) :: Neuman_itr
    integer, intent(in) :: independ_nvec
    integer, intent(in) :: nBlock
    integer, intent(inout) :: iret
    integer CAGMRES_sstep
    integer CAGMRES_tstep
    integer CAGMRES_basis
    integer CAGMRES_QRflag
    integer solver
    double precision absrtolmax

    absrtolmax = 0.0d0
    solver = KSP_SOLVER_CAGMRES

    call parcel_ddm_main( &
         icomm,-1 , nrow, &
         nnd, ioff_dia, val_dia, &
         num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
         ndiv_grids, num_grids, ioff_grids, div_direc_th, precon_thblock, &
         iovlflag, addL, iflagAS,&
         b, x, -1,-1, &
         rtolmax, absrtolmax, -1, -1, ipreflag, &
         itrmax, reshistory, &
         Neuman_itr,&
         independ_nvec, &
         nBlock, &
         iret, &
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         CAGMRES_sstep, &
         CAGMRES_tstep, &
         CAGMRES_basis, &
         CAGMRES_QRflag, &
         solver &
         )
  end subroutine parcel_dcagmres_ddm

  subroutine parcel_dcaarnoldi_ddm(&
       icomm, &
       x,b, &
       nrow, &
       nnd, ioff_dia, val_dia, &
       num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
       ndiv_grids, num_grids, ioff_grids,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS,&      
       iovlflag, &
       precon_thblock, &
       Neuman_itr,&
       independ_nvec, &
       nBlock, &
       CAARNOLDI_sstep, &
       CAARNOLDI_tstep, &
       CAARNOLDI_basis, &
       CAARNOLDI_QRflag, &
       EigenValueReal_DDM, &
       EigenValueImag_DDM, &
       EigenVector_DDM, &
       EigenErrorNorm_DDM, &
       EmaxID_DDM, &
       EminID_DDM, &
       iret &
       )
    use PARCEL_define
    use init_ddm_param, only:&
         parcel_ddm_free
    use EigenSolver
    use ddm_kernel
    use PARCEL_integer
    implicit none

    !arguments
    integer, intent(in) :: icomm !MPI communicator
    integer(kind=intType), intent(in) :: nrow !number of local matrix rows

    double precision, intent(in) :: b(nrow) !b vector
    double precision, intent(inout) :: x(nrow) !x vector
!    integer, intent(inout) :: oloop
    !    integer, intent(in) :: s ,t
    !    integer, intent(in) :: bflag 
    !    integer, intent(in) :: qrflag 
    integer  ipreflag !pre-conditioner flag
    integer, intent(in) :: nnd 
    integer, intent(in) :: ioff_dia(nnd) !offset
    double precision, intent(in) :: val_dia(nrow*nnd) !values
    integer, intent(in) :: num_neighbor_ranks
    integer, intent(in) :: neighbor_ranks( num_neighbor_ranks )
    integer, intent(in) :: margin
    integer, intent(in) :: ncomp_grids
    integer, intent(in) :: ndiv_grids(ncomp_grids)
    integer, intent(in) :: num_grids(ncomp_grids)
    integer, intent(in) :: ioff_grids(nnd,ncomp_grids)
    integer, intent(in) :: div_direc_th
    integer, intent(in) :: precon_thblock
    integer, intent(in) :: iovlflag ! 0: no hiding, 1: one time, 2: each time
    integer, intent(inout) :: addL
    integer, intent(in) :: iflagAS         
    integer, intent(in) :: Neuman_itr
    integer, intent(in) :: independ_nvec
    integer, intent(in) :: nBlock
    integer, intent(inout) :: iret

    integer CAARNOLDI_sstep
    integer CAARNOLDI_tstep
    integer CAARNOLDI_basis
    integer CAARNOLDI_QRflag
    double precision EigenValueReal_DDM(CAARNOLDI_sstep*CAARNOLDI_tstep)
    double precision EigenValueImag_DDM(CAARNOLDI_sstep*CAARNOLDI_tstep)
    double precision EigenVector_DDM(CAARNOLDI_sstep*CAARNOLDI_tstep*nrow)
    double precision EigenErrorNorm_DDM(CAARNOLDI_sstep*CAARNOLDI_tstep)
    integer EmaxID_DDM
    integer EminID_DDM

    integer solver
    double precision absrtolmax
    double precision  rtolmax
    double precision reshistory(CAARNOLDI_sstep*CAARNOLDI_tstep+1) 
    integer itrmax
    type(Eigen_Type) :: EigenType_

    absrtolmax = 0.0d0
    rtolmax = 0.0d0
    reshistory = 0.0d0
    itrmax = CAARNOLDI_sstep*CAARNOLDI_tstep
    solver = KSP_SOLVER_CAARNOLDI

    call parcel_ddm_main( &
         icomm,-1 , nrow, &
         nnd, ioff_dia, val_dia, &
         num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
         ndiv_grids, num_grids, ioff_grids, div_direc_th, precon_thblock, &
         iovlflag, addL, iflagAS,&
         b, x, -1,-1, &
         rtolmax, absrtolmax, -1, -1, ipreflag, &
         itrmax, reshistory, &
         Neuman_itr,&
         independ_nvec, &
         nBlock, &
         iret, &
         -1,&
         -1,&
         -1,&
         -1,&
         -1,&
         CAARNOLDI_sstep, &
         CAARNOLDI_tstep, &
         CAARNOLDI_basis, &
         CAARNOLDI_QRflag, &
         -1, &
         -1, &
         -1, &
         -1, &
         solver, &
         EigenType=EigenType_ &
         )

    EigenValueReal_DDM(1:CAARNOLDI_sstep*CAARNOLDI_tstep)   = EigenType_%EigenValueReal(1:CAARNOLDI_sstep*CAARNOLDI_tstep)
    EigenValueImag_DDM(1:CAARNOLDI_sstep*CAARNOLDI_tstep)   = EigenType_%EigenValueImag(1:CAARNOLDI_sstep*CAARNOLDI_tstep)
    EigenErrorNorm_DDM(1:CAARNOLDI_sstep*CAARNOLDI_tstep)   = EigenType_%EigenErrorNorm(1:CAARNOLDI_sstep*CAARNOLDI_tstep)
    EigenVector_DDM(1:CAARNOLDI_sstep*CAARNOLDI_tstep*nrow) = EigenType_%EigenVector(1:CAARNOLDI_sstep*CAARNOLDI_tstep*nrow)
    EmaxID_DDM = EigenType_%MaxEigenID 
    EminID_DDM = EigenType_%MinEigenID

  end subroutine parcel_dcaarnoldi_ddm

  subroutine parcel_dcbcg_ddm(&
       icomm, &
       x,b, &
       nrow, &
       nnd, ioff_dia, val_dia, &
       num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
       ndiv_grids, num_grids, ioff_grids,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS,&      
       itrmax, rtolmax, reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr,&
       independ_nvec, &
       nBlock, &
       CBCG_kstep, &
       CBCG_Eigenflag, &
       power_method_itrmax, &
       CAARNOLDI_sstep, &
       CAARNOLDI_tstep, &
       CAARNOLDI_basis, &
       CAARNOLDI_QRflag, &
       iret &
       )
    use PARCEL_define
    use ddm_kernel
    use PARCEL_integer
    implicit none

    !arguments
    integer, intent(in) :: icomm !MPI communicator
    integer(kind=intType), intent(in) :: nrow !number of local matrix rows

    double precision, intent(in) :: b(nrow) !b vector
    double precision, intent(inout) :: x(nrow) !x vector
!    integer, intent(inout) :: oloop
    double precision, intent(inout) :: rtolmax
    !    integer, intent(in) :: s ,t
    !    integer, intent(in) :: bflag 
    !    integer, intent(in) :: qrflag 
    integer  ipreflag !pre-conditioner flag
    integer, intent(inout) :: itrmax !number of reshistoryidual history
    double precision,intent(inout) :: reshistory(itrmax) !reshistoryidual history

    integer, intent(in) :: nnd 
    integer, intent(in) :: ioff_dia(nnd) !offset
    double precision, intent(in) :: val_dia(nrow*nnd) !values
    integer, intent(in) :: num_neighbor_ranks
    integer, intent(in) :: neighbor_ranks( num_neighbor_ranks )
    integer, intent(in) :: margin
    integer, intent(in) :: ncomp_grids
    integer, intent(in) :: ndiv_grids(ncomp_grids)
    integer, intent(in) :: num_grids(ncomp_grids)
    integer, intent(in) :: ioff_grids(nnd,ncomp_grids)
    integer, intent(in) :: div_direc_th
    integer, intent(in) :: precon_thblock
    integer, intent(in) :: iovlflag ! 0: no hiding, 1: one time, 2: each time
    integer, intent(inout) :: addL
    integer, intent(in) :: iflagAS         
    integer, intent(in) :: Neuman_itr
    integer, intent(in) :: independ_nvec
    integer, intent(in) :: nBlock
    integer, intent(inout) :: iret
    integer CBCG_kstep
    integer CBCG_Eigenflag
    integer power_method_itrmax
    integer CAARNOLDI_sstep
    integer CAARNOLDI_tstep
    integer CAARNOLDI_basis
    integer CAARNOLDI_QRflag
    integer solver
    double precision absrtolmax

    absrtolmax = 0.0d0
    solver = KSP_SOLVER_CBCG

    call parcel_ddm_main( &
         icomm,-1 , nrow, &
         nnd, ioff_dia, val_dia, &
         num_neighbor_ranks, neighbor_ranks, margin, ncomp_grids,&
         ndiv_grids, num_grids, ioff_grids, div_direc_th, precon_thblock, &
         iovlflag, addL, iflagAS,&
         b, x, -1,-1, &
         rtolmax, absrtolmax, -1, -1, ipreflag, &
         itrmax, reshistory, &
         Neuman_itr,&
         independ_nvec, &
         nBlock, &
         iret, &
         -1,&
         -1,&
         CBCG_kstep, &
         CBCG_Eigenflag, &
         power_method_itrmax, &
         CAARNOLDI_sstep, &
         CAARNOLDI_tstep, &
         CAARNOLDI_basis, &
         CAARNOLDI_QRflag, &
         -1, &
         -1, &
         -1, &
         -1, &
         solver &
         )
  end subroutine parcel_dcbcg_ddm
#endif

end module ddm_wrapper
