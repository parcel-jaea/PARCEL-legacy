! #define CUDA_AWARE 0
! #define PRINT_TIME 0
! #define PRINT_TIME 1
! #define USE_C_MPI 0  ! // C言語MPIのテスト用

subroutine PARCEL_Init(myrank)
  use cu_common_wrapper,only : set_GPU,init_cudahandle
  implicit none
  integer myrank
#ifdef CUDA_SOLVER
  call set_GPU(myrank)
  call init_cudahandle()
#endif
end subroutine PARCEL_Init

subroutine PARCEL_Free()
  use cu_common_wrapper,only : free_cudahandle
  implicit none
#ifdef CUDA_SOLVER
  call free_cudahandle()
#endif
end subroutine PARCEL_Free

! フラグを読み込み設定されていない場合は入力されたCPU配列を利用する
subroutine init_csrMatrix_GPU()
  implicit none

end subroutine init_csrMatrix_GPU

subroutine set_csrMatrix_GPU()
  implicit none

end subroutine set_csrMatrix_GPU

module cu_Matrix
  use PARCEL_integer
  use c_wrapper
  implicit none
  interface

     subroutine Matvec_ddm_FP64_GPU(&
		       numGPUblock, &
		       numGPUthread, &
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
		       A, &
		       x, offset_x, &
		       y, offset_y, &		
		       ioff_A &
          ) &
          bind(c,name = "cuda_MatVec_ddm")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value :: nx
       integer(kind=intType),value :: ny
       integer(kind=intType),value :: nz
       integer(kind=intType),value :: stmx
       integer(kind=intType),value :: stmy
       integer(kind=intType),value :: stmz
       integer(kind=intType),value :: mx
       integer(kind=intType),value :: mxy
       integer(kind=intType),value :: m
       integer(kind=intType),value :: nnd
       type(devicemem)  A
       type(devicemem)  x
       integer(kind=intType),value :: offset_x
       type(devicemem)  y		
       integer(kind=intType),value :: offset_y
       type(devicemem)  ioff_A 
       
     end subroutine Matvec_ddm_FP64_GPU

     subroutine jacobi_iter1_ddm_FP64_GPU(&
		       numGPUblock, &
		       numGPUthread, &
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
                        A, &
                        x, offset_x, &
                        y, offset_y, &		
                        ioff_A &
          ) &
          bind(c,name = "cuda_jacobi_iter1_ddm")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value :: nx
       integer(kind=intType),value :: ny
       integer(kind=intType),value :: nz
       integer(kind=intType),value :: stmx
       integer(kind=intType),value :: stmy
       integer(kind=intType),value :: stmz
       integer(kind=intType),value :: mx
       integer(kind=intType),value :: mxy
       integer(kind=intType),value :: m
       integer(kind=intType),value :: nnd
       type(devicemem)  A
       type(devicemem)  x
       integer(kind=intType),value :: offset_x
       type(devicemem)  y		
       integer(kind=intType),value :: offset_y
       type(devicemem)  ioff_A 
       
     end subroutine jacobi_iter1_ddm_FP64_GPU

     subroutine cuda_MatVec_dia_init_GPU( &
          nnd,offset &			    
          )&
          bind(c,name = "cuda_MatVec_dia_init")
       import :: intType
       implicit none
       integer(kind=intType),value ::  nnd
       integer(kind=intType) offset(nnd)

     end subroutine cuda_MatVec_dia_init_GPU

     subroutine MatVec_crs_FP64_GPU( &
          numGPUblock,&
          numGPUthread,&
          n,&
          istart,&
          iend,&
          crsA, &
          crsRow_ptr, &
          crsCol, &
          iLs, &
          xL,offset_xL,&
          xD,offset_xD,&
          xU,offset_xU,&
          Ax,offset_Ax &
          )&
          bind(c,name = "cuda_MatVec_crs")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value ::  n
       integer(kind=intType),value ::  istart
       integer(kind=intType),value ::  iend
       type(devicemem) crsA 
       type(devicemem) crsRow_ptr 
       type(devicemem) crsCol 
       type(devicemem) iLs 
       type(devicemem) xL
       integer(kind=intType),value :: offset_xL
       type(devicemem) xD
       integer(kind=intType),value :: offset_xD
       type(devicemem) xU
       integer(kind=intType),value :: offset_xU
       type(devicemem) Ax 
       integer(kind=intType),value :: offset_Ax
         
     end subroutine MatVec_crs_FP64_GPU

     subroutine Diag_Transfer_FP64toBF16_crs_FP64_GPU( &
          numGPUblock,&
          numGPUthread,&
          n,&
          crsA_BF16, &
          crsA_FP64, &
          crsRow_ptr, &
          crsCol, &
          iLs, &
          mode &
          )&
          bind(c,name = "cuda_Diag_Transfer_FP64toBF16_crs")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value ::  n
       type(devicemem) crsA_BF16
       type(devicemem) crsA_FP64
       type(devicemem) crsRow_ptr 
       type(devicemem) crsCol 
       type(devicemem) iLs 
       integer,value :: mode
         
     end subroutine Diag_Transfer_FP64toBF16_crs_FP64_GPU

     subroutine Diag_Transfer_FP64toBF16_localcrs_FP64_GPU( &
          numGPUblock,&
          numGPUthread,&
          n,&
          crsA_FP64, &
          crsRow_ptr, &
          crsCol, &
          mode, &
          mode_diag &
          )&
          bind(c,name = "cuda_Diag_Transfer_FP64toBF16_localcrs")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value ::  n
       type(devicemem) crsA_FP64
       type(devicemem) crsRow_ptr 
       type(devicemem) crsCol 
       integer,value :: mode
       integer,value :: mode_diag
         
     end subroutine Diag_Transfer_FP64toBF16_localcrs_FP64_GPU

     subroutine JACOBI_MatBF16_crs_FP64_GPU( &
          numGPUblock,&
          numGPUthread,&
          n,&
          istart,&
          iend,&
          crsA, &
          crsRow_ptr, &
          crsCol, &
          iLs, &
          b ,offset_b , &
          xL,offset_xL,&
          xD,offset_xD,&
          xU,offset_xU,&
          Ax,offset_Ax &
          )&
          bind(c,name = "cuda_JACOBI_MatBF16_crs")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value ::  n
       integer(kind=intType),value ::  istart
       integer(kind=intType),value ::  iend
       type(devicemem) crsA 
       type(devicemem) crsRow_ptr 
       type(devicemem) crsCol 
       type(devicemem) iLs 
       type(devicemem) b
       integer(kind=intType),value :: offset_b
       type(devicemem) xL
       integer(kind=intType),value :: offset_xL
       type(devicemem) xD
       integer(kind=intType),value :: offset_xD
       type(devicemem) xU
       integer(kind=intType),value :: offset_xU
       type(devicemem) Ax 
       integer(kind=intType),value :: offset_Ax
         
     end subroutine JACOBI_MatBF16_crs_FP64_GPU

     subroutine JACOBI_MatFP64_crs_FP64_GPU( &
          numGPUblock,&
          numGPUthread,&
          n,&
          istart,&
          iend,&
          crsA, &
          crsRow_ptr, &
          crsCol, &
          iLs, &
          b ,offset_b , &
          xL,offset_xL,&
          xD,offset_xD,&
          xU,offset_xU,&
          Ax,offset_Ax &
          )&
          bind(c,name = "cuda_JACOBI_MatFP64_crs")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value ::  n
       integer(kind=intType),value ::  istart
       integer(kind=intType),value ::  iend
       type(devicemem) crsA 
       type(devicemem) crsRow_ptr 
       type(devicemem) crsCol 
       type(devicemem) iLs 
       type(devicemem) b
       integer(kind=intType),value :: offset_b
       type(devicemem) xL
       integer(kind=intType),value :: offset_xL
       type(devicemem) xD
       integer(kind=intType),value :: offset_xD
       type(devicemem) xU
       integer(kind=intType),value :: offset_xU
       type(devicemem) Ax 
       integer(kind=intType),value :: offset_Ax
         
     end subroutine JACOBI_MatFP64_crs_FP64_GPU

     ! CRS形式はcusparceを利用する実装にしたので不要
     subroutine MatVec_crs_neumann_FP64_GPU( &
          numGPUblock,&
          numGPUthread,&
          n,&
          istart,&
          iend,&
          crsA, &
          crsRow_ptr, &
          crsCol, &
          iLs, &
          xL,offset_xL,&
          xD,offset_xD,&
          xU,offset_xU,&
          Ax,offset_Ax &
          )&
          bind(c,name = "cuda_MatVec_neumann_crs")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value ::  n
       integer(kind=intType),value ::  istart
       integer(kind=intType),value ::  iend
       type(devicemem) crsA 
       type(devicemem) crsRow_ptr 
       type(devicemem) crsCol 
       type(devicemem) iLs 
       type(devicemem) xL
       integer(kind=intType),value :: offset_xL
       type(devicemem) xD
       integer(kind=intType),value :: offset_xD
       type(devicemem) xU
       integer(kind=intType),value :: offset_xU
       type(devicemem) Ax 
       integer(kind=intType),value :: offset_Ax
         
     end subroutine MatVec_crs_neumann_FP64_GPU

     subroutine MatVec_dia_FP64_GPU( &
          numGPUblock,&
          numGPUthread,&
          n,&
          nnd,&
          istart,&
          nBlock,&
          x,offset_x, &
          y,offset_y, &
          diaA, &
          offset , &
          iLs, &
          rbuf &
          )&
          bind(c,name = "cuda_MatVec_dia")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value ::  n
       integer(kind=intType),value ::  nnd
       integer(kind=intType),value ::  istart
       integer(kind=intType),value ::  nBlock
       type(devicemem) x
       integer(kind=intType),value :: offset_x
       type(devicemem) y
       integer(kind=intType),value :: offset_y

       type(devicemem) diaA
       type(devicemem) offset
       type(devicemem) iLs
       type(devicemem) rbuf
         
     end subroutine MatVec_dia_FP64_GPU

     subroutine MatVec_cbcg1_dia_FP64_GPU( &
          numGPUblock,&
          numGPUthread,&
          n,&
          nnd,&
          istart,&
          nBlock,&
          x ,offset_x, &
          Ax,offset_Ax, &
          y ,offset_y, &
          z ,offset_z, &
          alpha, &
          beta , &
          diaA, &
          offset , &
          iLs, &
          rbuf &
          )&
          bind(c,name = "cuda_MatVec_cbcg1_dia")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value ::  n
       integer(kind=intType),value ::  nnd
       integer(kind=intType),value ::  istart
       integer(kind=intType),value ::  nBlock
       type(devicemem) x
       integer(kind=intType),value :: offset_x
       type(devicemem) Ax
       integer(kind=intType),value :: offset_Ax
       type(devicemem) y
       integer(kind=intType),value :: offset_y
       type(devicemem) z
       integer(kind=intType),value :: offset_z

       real*8,value :: alpha
       real*8,value :: beta

       type(devicemem) diaA
       type(devicemem) offset
       type(devicemem) iLs
       type(devicemem) rbuf
         
     end subroutine MatVec_cbcg1_dia_FP64_GPU

     subroutine MatVec_cbcg2_dia_FP64_GPU( &
          numGPUblock,&
          numGPUthread,&
          n,&
          nnd,&
          istart,&
          nBlock,&
          x ,offset_x, &
          Ax,offset_Ax, &
          y ,offset_y, &
          z ,offset_z, &
          w ,offset_w, &
          alpha, &
          beta , &
          gamma , &
          diaA, &
          offset , &
          iLs, &
          rbuf &
          )&
          bind(c,name = "cuda_MatVec_cbcg2_dia")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value ::  n
       integer(kind=intType),value ::  nnd
       integer(kind=intType),value ::  istart
       integer(kind=intType),value ::  nBlock
       type(devicemem) x
       integer(kind=intType),value :: offset_x
       type(devicemem) Ax
       integer(kind=intType),value :: offset_Ax
       type(devicemem) y
       integer(kind=intType),value :: offset_y
       type(devicemem) z
       integer(kind=intType),value :: offset_z
       type(devicemem) w
       integer(kind=intType),value :: offset_w

       real*8,value :: alpha
       real*8,value :: beta
       real*8,value :: gamma

       type(devicemem) diaA
       type(devicemem) offset
       type(devicemem) iLs
       type(devicemem) rbuf
         
     end subroutine MatVec_cbcg2_dia_FP64_GPU

     subroutine MatVec_dia_neumann_FP64_GPU( &
          numGPUblock,&
          numGPUthread,&
          n,&
          nnd,&
          istart,&
          nBlock,&
          x,offset_x, &
          y,offset_y, &
          diaA, &
          offset , &
          iLs, &
          rbuf &
          )&
          bind(c,name = "cuda_MatVec_neumann_dia")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value ::  n
       integer(kind=intType),value ::  nnd
       integer(kind=intType),value ::  istart
       integer(kind=intType),value ::  nBlock
       type(devicemem) x
       integer(kind=intType),value :: offset_x
       type(devicemem) y
       integer(kind=intType),value :: offset_y
       type(devicemem) diaA
       type(devicemem) offset
       type(devicemem) iLs
       type(devicemem) rbuf
         
     end subroutine MatVec_dia_neumann_FP64_GPU

  end interface

!contains

end module cu_Matrix


module cu_krylov
  use cu_common_wrapper
  use cu_Matrix,only:cuda_MatVec_dia_init_GPU
  implicit none

contains

  subroutine set_CUSPARCE(method)
    !    use krylov
    !    use matrix_kernel
    use krylov,only : KSP
    use matrix,only : CRS_Mat
    use matrix_kernel,only : &
         Localcrs_to_Clang, &
         Localcrs_to_Fortran, &
         crs2Localcrs, &
         dia2crs, &
         free_crsMat, &
         get_crs_Row, &
         init_MatVec_crs, &
         init_crsMat, &
         make_DiagBlock_Localcrs, &
         get_crscol_ncol
    use PARCEL_define,only : &
         MATRIX_FORMAT_DIA, &
         KSP_PRECONDITION_NONE, &
         KSP_PRECONDITION_BLOCKJACOBI, &
         KSP_PRECONDITION_NEUMANN, &
         KSP_PRECONDITION_ADDITIVESCHWARZ
    use cu_Matrix,only : Diag_Transfer_FP64toBF16_localcrs_FP64_GPU

    implicit none

    type(KSP) :: method
    integer myrank
    integer(kind=intType) n,n_precon,nn
    integer(kind=intType) istart,iend
    integer(kind=intType) istart_AS,iend_AS
    integer(kind=intType) nnz,nnz_precon,nnd
    integer(kind=intType) nCol,nColL,nColU,nColD
    integer(kind=intType) k,j
    integer(kind=intType) parcel_int

    type(CRS_Mat) :: crsMat ! DIA形式からのCRS形式の変換で利用 
    !    type(Local_CRS_Mat) :: LocalcrsMat
    ! type(CRS_Mat) :: LocalcrsMat
    real*8,allocatable  :: LocalcrsMat_crsA(:)
    integer,allocatable :: LocalcrsMat_crsRow_ptr(:)
    integer,allocatable :: LocalcrsMat_crsCol(:)
    integer MatFormat
    integer(kind=intType),allocatable :: Colsort(:)
    type(devicemem) :: crsA_FP64
    type(devicemem) :: crsA_iLs
    integer(kind=intType),parameter :: offset_zero=0
    integer(kind=intType) work_int_size
    integer,allocatable :: work_int(:)
    real*8 stime,etime

    myrank = method%commVec%myrank
    MatFormat = method%prm%MatFormat 

    ! CUSPARCE initialize
    n      = method%prm%n
    istart = method%prm%istart
    iend   = method%prm%iend
    nnz    = method%prm%nnz
    nCol   = n + method%crsMat%nColL + method%crsMat%nColU

    ! ! write(*,*)"nnz  = ",nnz
    ! init

    !  ! write(*,*)" set_CUSPARCE : @ 1 @"
    !    call init_LocalcrsMat(LocalcrsMat,n,nnz)
    !@@@@@@    call init_crsMat(LocalcrsMat,n,nnz)

    stime=MPI_Wtime()
    allocate( LocalcrsMat_crsA(nnz) )
    allocate( LocalcrsMat_crsRow_ptr(n+1) )
    allocate( LocalcrsMat_crsCol(nnz) )

    if(MatFormat.eq.MATRIX_FORMAT_DIA)then
       n      = method%prm%n
       istart = method%prm%istart
       nnz    = method%prm%nnz
       nnd    =  method%prm%nnd

       !  ! ! write(*,*)"set_CUSPARCE : @ 2-0 @",n,istart,nnz,nnd

       ! ---------------------
       call init_crsMat(crsMat,n,nnz)      
       call dia2crs(method%diaMat,crsMat,n,istart,nnd,nnz)
       ! write(*,*)" set_CUSPARCE : @ 2 @"

       allocate( Colsort(nnz) )
       ! write(*,*)" set_CUSPARCE : @ 2-1 @"
       call get_crsCol_nCol(&
            method%prm, crsMat%crsCol, Colsort,&
            nCol,nColL,nColU,nColD)
       ! write(*,*)" set_CUSPARCE : @ 2-2 @"

       crsMat%nColL = nColL
       crsMat%nColU = nColU
       call get_crs_Row(method%prm,crsMat,nCol,Colsort)
       ! write(*,*)" set_CUSPARCE : @ 2-3 @"

       allocate(crsMat%iLs(nnz))
       ! write(*,*)" set_CUSPARCE : @ 2-4 @"
       call init_MatVec_crs(method%prm,crsMat,method%nBlock)       
       ! write(*,*)" set_CUSPARCE : @ 2-5 @"
       nCol   = n + crsMat%nColL + crsMat%nColU
       ! ---------------------
       ! write(*,*)" set_CUSPARCE : @ 3 @"

       call crs2Localcrs( &
            crsMat, &
            LocalcrsMat_crsA, &
            LocalcrsMat_crsRow_ptr, &
            LocalcrsMat_crsCol, &
            n,istart,nnz, &
            n, nnz)

       deallocate(crsMat%iLs)
       deallocate(Colsort)
       call free_crsMat(crsMat)

       ! write(*,*)" set_CUSPARCE :         @ 4 @"
    else
       call crs2Localcrs( &
            method%crsMat, &
            LocalcrsMat_crsA, &
            LocalcrsMat_crsRow_ptr, &
            LocalcrsMat_crsCol, &
            n,istart,nnz, &
            n, nnz)
    endif
    etime=MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"        @ 4-6-1 @",etime-stime
    endif
#endif
    stime=MPI_Wtime()

    !    call view_localcrs(n,LocalcrsMat)
    call Localcrs_to_Clang( &
         LocalcrsMat_crsA, &
         LocalcrsMat_crsRow_ptr, &
         LocalcrsMat_crsCol, &
         n,istart,nnz)

    etime=MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"        @ 4-6-2 @",etime-stime
    endif
#endif
    stime=MPI_Wtime()

    ! write(*,*)" set_CUSPARCE : @ 5 @", size(LocalcrsMat_crsA)
!    nn = size(LocalcrsMat_crsA)
    nn = method%prm%nnz
    call cudaMalloc_fp64(method%cusparce_crsMatrix%crsA,nn)
    call cudaMemcpy_h2d_fp64(LocalcrsMat_crsA,method%cusparce_crsMatrix%crsA,nn,offset_zero,offset_zero)

!    work_int_size = max( n+1 , nnz )
    work_int_size = nnz + 1
    if(allocated(work_int))then
       deallocate( work_int )  
    endif
    allocate(work_int(work_int_size))

    nn = n+1
    !$OMP PARALLEL DO PRIVATE(k) 
    do k = 1,nn
       work_int(k) = LocalcrsMat_crsRow_ptr(k)
    enddo
    ! write(*,*)" set_CUSPARCE : @ 5-1 @"
    call cudaMalloc_int(method%cusparce_crsMatrix%crsRow_ptr,nn)
    call cudaMemcpy_h2d_int(work_int,method%cusparce_crsMatrix%crsRow_ptr,nn,offset_zero,offset_zero)

!    nn = size(LocalcrsMat_crsCol)
    nn = method%prm%nnz
    !$OMP PARALLEL DO PRIVATE(k) 
    do k = 1,nn
       work_int(k) = LocalcrsMat_crsCol(k)
    enddo
    ! write(*,*)" set_CUSPARCE : @ 5-2 @"
    call cudaMalloc_int(method%cusparce_crsMatrix%crsCol,nn)
    call cudaMemcpy_h2d_int(LocalcrsMat_crsCol,method%cusparce_crsMatrix%crsCol,nn,offset_zero,offset_zero)
    if(allocated(work_int))then
       deallocate( work_int )  
    endif

    etime=MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"        @ 4-6-3 @",etime-stime
    endif
#endif
    stime=MPI_Wtime()

    ! write(*,*)" set_CUSPARCE : @ 5-3 @"
    ! CUDA用の前処理行列作成
    nnz_precon = method%prm%nnz
    if(method%ipreflag.eq.KSP_PRECONDITION_NEUMANN)then
       ! Neumann用の行列をコピーする
       n_precon   = n
       nnz_precon = method%prm%nnz

       if(MatFormat.eq.MATRIX_FORMAT_DIA)then

       else
          call crs2Localcrs( &
               method%preconNeuman%crsMat, &
               LocalcrsMat_crsA, &
               LocalcrsMat_crsRow_ptr, &
               LocalcrsMat_crsCol, &
               n,istart,nnz, &
               n, nnz)
       endif

    else
       ! cusparseのILU領域を加法シュヴァルツ,ブロックヤコビで切り替える
       if(method%ipreflag.eq.KSP_PRECONDITION_ADDITIVESCHWARZ)then

          if(method%commVec%npes.eq.1)then
             ! 拡張領域がなければBlock Jacobiで計算 : CPU側で1プロセスの場合前処理フラグをブロックヤコビに切り替えているため動作未確認.
             n_precon = n
             nnz_precon = method%prm%nnz
             call Localcrs_to_Fortran( &
                  LocalcrsMat_crsA, &
                  LocalcrsMat_crsRow_ptr, &
                  LocalcrsMat_crsCol, &
                  n_precon,istart,nnz_precon)
             call make_DiagBlock_Localcrs( &
                  n,nnz_precon,&
                  LocalcrsMat_crsA, &
                  LocalcrsMat_crsRow_ptr, &
                  LocalcrsMat_crsCol &
                  )
          else
             ! このmethod%prm_AS%nはプロセスごと変えなければcusparceで正しく計算できないのでは?
             n_precon   = method%prm_AS%n
             nnz_precon = method%prm_AS%nnz

             ! 加法シュヴァルツの場合LocalcrsMatではメモリが足りないためメモリ確保を再度行う
             if(allocated(LocalcrsMat_crsA))then
                deallocate( LocalcrsMat_crsA )  
             endif
             if(allocated(LocalcrsMat_crsRow_ptr))then
                deallocate( LocalcrsMat_crsRow_ptr )  
             endif
             if(allocated(LocalcrsMat_crsCol))then
                deallocate( LocalcrsMat_crsCol )  
             endif
             allocate( LocalcrsMat_crsA(nnz_precon) )
             allocate( LocalcrsMat_crsRow_ptr(n_precon+1) )
             allocate( LocalcrsMat_crsCol(nnz_precon) )

             !  istart がaddL領域分小さくなる
             istart_AS = istart - method%addL 
             if(istart_AS<=0)then
                istart_AS = 1
             endif
             iend_AS = iend + method%addL
             if(iend_AS>n)then
                iend_AS = n
             endif

             !  加法シュワルツ法のローカルインデックスを作成
             !  拡張領域のデータしか存在しないため拡張領域の開始位置に相当する行番号でインデックスをずらせば良い
             allocate(method%precon_crsMatAS%iLs(nnz_precon))
             do k = 1,n_precon
                do j = method%precon_crsMatAS%crsRow_ptr(k),method%precon_crsMatAS%crsRow_ptr(k+1)-1
                   method%precon_crsMatAS%iLs(j) = method%precon_crsMatAS%crsCol(j)-1+istart_AS
                enddo
             enddo
#if 0
             myrank = method%commVec%myrank
             if( myrank.eq. 1 )then
!                do k = 1,n_precon
                do k = 1,10
                   do j = method%precon_crsMatAS%crsRow_ptr(k),method%precon_crsMatAS%crsRow_ptr(k+1)-1
                      ! ! write(*,*)k,method%precon_crsMatAS%crsCol(j) - istart_AS+1 ,method%precon_crsMatAS%crsA(j)
                   enddo
                enddo
             endif
#endif
             call crs2Localcrs( &
                  method%precon_crsMatAS, &
                  LocalcrsMat_crsA, &
                  LocalcrsMat_crsRow_ptr, &
                  LocalcrsMat_crsCol, &
                  n_precon,&
                  istart_AS, &
                  nnz_precon, &
                  n_precon, &
                  nnz_precon)

#if 0
             myrank = method%commVec%myrank
             if( myrank.eq. 1 )then
                do k = 1,n_precon
                   do j = LocalcrsMat_crsRow_ptr(k),LocalcrsMat_crsRow_ptr(k+1)-1
                      ! ! write(*,*)k,LocalcrsMat_crsCol(j),LocalcrsMat_crsA(j)
                   enddo
                enddo
             endif
#endif
             call make_DiagBlock_Localcrs( &
                  n_precon,nnz_precon,&
                  LocalcrsMat_crsA, &
                  LocalcrsMat_crsRow_ptr, &
                  LocalcrsMat_crsCol &
                  )
             if(allocated(method%precon_crsMatAS%iLs))then
                deallocate( method%precon_crsMatAS%iLs )  
             endif
          endif

          ! AS用ベクトル確保
          call cudaMalloc_fp64(method%C_ivec_AS,n_precon)
          call cusparseCreateDnVec_FP64(method%C_ivec_AS,n_precon,n_precon,offset_zero)
          call cudaMalloc_fp64(method%C_ovec_AS,n_precon)
          call cusparseCreateDnVec_FP64(method%C_ovec_AS,n_precon,n_precon,offset_zero)

          ! 通信バッファ
          parcel_int = max(1,method%addL)
          call cudaMalloc_FP64(method%C_s_buffL,parcel_int)
          call cudaMalloc_FP64(method%C_r_buffL,parcel_int)
          call cudaMalloc_FP64(method%C_s_buffU,parcel_int)
          call cudaMalloc_FP64(method%C_r_buffU,parcel_int)

          ! MPIリクエスト
          call malloc_MPI_Request(method%C_s_reqL,1) 
          call malloc_MPI_Request(method%C_r_reqL,1) 
          call malloc_MPI_Request(method%C_s_reqU,1) 
          call malloc_MPI_Request(method%C_r_reqU,1) 

          ! MPIステータス               
          call malloc_MPI_Status(method%C_s_staL,1)
          call malloc_MPI_Status(method%C_r_staL,1)
          call malloc_MPI_Status(method%C_s_staU,1)
          call malloc_MPI_Status(method%C_r_staU,1)

       else

!          if(method%ipreflag.eq.KSP_PRECONDITION_BLOCKJACOBI)then
             ! Block Jacobi
             n_precon = n
             nnz_precon = method%prm%nnz
             call Localcrs_to_Fortran( &
                  LocalcrsMat_crsA, &
                  LocalcrsMat_crsRow_ptr, &
                  LocalcrsMat_crsCol, &
                  n_precon,istart,nnz_precon)
             call make_DiagBlock_Localcrs( &
                  n,nnz_precon,&
                  LocalcrsMat_crsA, &
                  LocalcrsMat_crsRow_ptr, &
                  LocalcrsMat_crsCol &
                  )    
!          endif

       endif
    endif
    etime=MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"        @ 4-6-4 @",etime-stime
    endif
#endif
    stime=MPI_Wtime()
!     stop
    ! write(*,*)" set_CUSPARCE : @ 5-4 @"
    if(  &
         (method%ipreflag.eq.KSP_PRECONDITION_NEUMANN) .and.  &
         (MatFormat.eq.MATRIX_FORMAT_DIA)  &
         )then
       method%cusparce_solver%ipreflag = KSP_PRECONDITION_NONE
    else
       ! write(*,*)" set_CUSPARCE : @ 5-4-1 @"
       call Localcrs_to_Clang( &
            LocalcrsMat_crsA, &
            LocalcrsMat_crsRow_ptr, &
            LocalcrsMat_crsCol, &
            n_precon,istart,nnz_precon)

       ! write(*,*)" set_CUSPARCE : @ 5-4-2 @"
       nn = nnz_precon
       call cudaMalloc_fp64(method%cusparce_preconcrsMatrix%crsA,nn)
       call cudaMemcpy_h2d_fp64(LocalcrsMat_crsA,method%cusparce_preconcrsMatrix%crsA,nn,offset_zero,offset_zero)

       ! write(*,*)" set_CUSPARCE : @ 5-4-3 @"
!       work_int_size = max( (n_precon + 1), nnz_precon )
       work_int_size = nnz_precon+1
       if(allocated(work_int))then
          deallocate( work_int )  
       endif
       allocate(work_int(work_int_size))

       nn = n_precon + 1
       !$OMP PARALLEL DO PRIVATE(k) 
       do k = 1,nn
          work_int(k) = LocalcrsMat_crsRow_ptr(k)
       enddo
       call cudaMalloc_int(method%cusparce_preconcrsMatrix%crsRow_ptr,nn)
       call cudaMemcpy_h2d_int(work_int,method%cusparce_preconcrsMatrix%crsRow_ptr,nn,offset_zero,offset_zero)

       ! write(*,*)" set_CUSPARCE : @ 5-4-4 @"
       nn = nnz_precon
       !$OMP PARALLEL DO PRIVATE(k) 
       do k = 1,nn
          work_int(k) = LocalcrsMat_crsCol(k)
       enddo
       ! write(*,*)" set_CUSPARCE : @ 5-4-5 @"
       call cudaMalloc_int(method%cusparce_preconcrsMatrix%crsCol,nn)
       call cudaDeviceSynchronize() 
       ! write(*,*)" set_CUSPARCE : @ 5-4-6 @"
       call cudaMemcpy_h2d_int(work_int,method%cusparce_preconcrsMatrix%crsCol,nn,offset_zero,offset_zero)
       call cudaDeviceSynchronize() 
       ! write(*,*)" set_CUSPARCE : @ 5-4-7 @"

       if(allocated(work_int))then
          deallocate( work_int )  
       endif
         ! write(*,*)" set_CUSPARCE : @ 9 @"
       if(method%ipreflag.eq.KSP_PRECONDITION_NEUMANN)then
          method%cusparce_preconcrsMatrix%n    = n;
          method%cusparce_preconcrsMatrix%nCol = nCol;
          method%cusparce_preconcrsMatrix%nnz  = nnz_precon;
       else
          method%cusparce_preconcrsMatrix%n    = n_precon;
          method%cusparce_preconcrsMatrix%nCol = n_precon;
          method%cusparce_preconcrsMatrix%nnz  = nnz_precon;
       endif

       !    endif
       ! CUSPARCEのSpMV,ILUの初期化
       method%cusparce_solver%ipreflag = method%ipreflag
    endif
    etime=MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"        @ 4-6-5 @",etime-stime
    endif
#endif
    stime=MPI_Wtime()

    ! write(*,*)" set_CUSPARCE : @ 10 @"
    !    ! ! write(*,*)"precon n    = ",method%cusparce_preconcrsMatrix%n
    !    ! ! write(*,*)"precon nCol = ",method%cusparce_preconcrsMatrix%nCol

#if 0
    call cudaDeviceSynchronize() 
    ! call cuda_view_FP64(10,method%C_crsMat%crsA)
    ! call cuda_view_PARCEL_INT(10,method%C_crsMat%crsCol)
    ! call cuda_view_PARCEL_INT(10,method%C_crsMat%crsRow_ptr)
    if(myrank.eq.1)then
       call cuda_view_FP64(nnz,method%cusparce_crsMatrix%crsA)
       call cuda_view_PARCEL_INT(nnz,method%cusparce_crsMatrix%crsCol)
       call cuda_view_PARCEL_INT(n+1,method%cusparce_crsMatrix%crsRow_ptr)
    endif
    call cudaDeviceSynchronize() 
#endif
    
    call cusparse_init( &
         method%cusparce_solver, &
         ! デバイス側のCRS 行列 (ローカル座標)
         method%cusparce_crsMatrix%crsA, &
         method%cusparce_crsMatrix%crsCol, &
         method%cusparce_crsMatrix%crsRow_ptr, &
         n, &
         nCol, &
         nnz, &
         method%cusparce_preconcrsMatrix%crsA, &
         method%cusparce_preconcrsMatrix%crsCol, &
         method%cusparce_preconcrsMatrix%crsRow_ptr, &
         method%cusparce_preconcrsMatrix%n, &
         method%cusparce_preconcrsMatrix%nCol, &
         method%cusparce_preconcrsMatrix%nnz &
         )

    deallocate( LocalcrsMat_crsA )
    deallocate( LocalcrsMat_crsRow_ptr )
    deallocate( LocalcrsMat_crsCol )

    etime=MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"        @ 4-6-6 @",etime-stime
    endif
#endif
    stime=MPI_Wtime()

  end subroutine set_CUSPARCE

  subroutine set_VEC_GPU(method)
    use krylov
!    use vec_kernel
    implicit none
    type(KSP) :: method
    integer(kind=intType) n
    integer(kind=intType),parameter :: offset_zero=0
    ! --------
    ! Vector
    ! type(devicemem) :: C_r
    ! type(devicemem) :: C_r_hi
    ! type(devicemem) :: C_r_lo
    ! type(devicemem) :: C_x
    ! type(devicemem) :: C_x_hi
    ! type(devicemem) :: C_x_lo
    ! type(devicemem) :: C_b
    ! type(devicemem) :: C_b_hi
    ! type(devicemem) :: C_b_lo
    ! --------
    ! if(allocated(method%@))then
    !    n = size(method%@)
    !    call cudaMalloc_fp64(method%C_@,n)
    !    call cudaMemcpy_h2d_fp64(method%@,method%C_@,n)
    ! endif

    if(allocated(method%r))then
!       n = size(method%r)
       n = method%prm%n
       call cudaMalloc_fp64(method%C_r,n)
       call cudaMemcpy_h2d_fp64(method%r,method%C_r,n,offset_zero,offset_zero)
    endif
    if(allocated(method%r_hi))then
!       n = size(method%r_hi)
       n = method%prm%n
       call cudaMalloc_fp64(method%C_r_hi,n)
       call cudaMemcpy_h2d_fp64(method%r_hi,method%C_r_hi,n,offset_zero,offset_zero)
    endif
    if(allocated(method%r_lo))then
!       n = size(method%r_lo)
       n = method%prm%n
       call cudaMalloc_fp64(method%C_r_lo,n)
       call cudaMemcpy_h2d_fp64(method%r_lo,method%C_r_lo,n,offset_zero,offset_zero)
    endif

    if(allocated(method%x))then
!       n = size(method%x)
       n = method%prm%n
       call cudaMalloc_fp64(method%C_x,n)
       call cudaMemcpy_h2d_fp64(method%x,method%C_x,n,offset_zero,offset_zero)
    endif
    if(allocated(method%x_hi))then
!       n = size(method%x_hi)
       n = method%prm%n
       call cudaMalloc_fp64(method%C_x_hi,n)
       call cudaMemcpy_h2d_fp64(method%x_hi,method%C_x_hi,n,offset_zero,offset_zero)
    endif
    if(allocated(method%x_lo))then
!       n = size(method%x_lo)
       n = method%prm%n
       call cudaMalloc_fp64(method%C_x_lo,n)
       call cudaMemcpy_h2d_fp64(method%x_lo,method%C_x_lo,n,offset_zero,offset_zero)
    endif

    if(allocated(method%b))then
!       n = size(method%b)
       n = method%prm%n
       call cudaMalloc_fp64(method%C_b,n)
       call cudaMemcpy_h2d_fp64(method%b,method%C_b,n,offset_zero,offset_zero)
    endif
    if(allocated(method%b_hi))then
!       n = size(method%b_hi)
       n = method%prm%n
       call cudaMalloc_fp64(method%C_b_hi,n)
       call cudaMemcpy_h2d_fp64(method%b_hi,method%C_b_hi,n,offset_zero,offset_zero)
    endif
    if(allocated(method%b_lo))then
!       n = size(method%b_lo)
       n = method%prm%n
       call cudaMalloc_fp64(method%C_b_lo,n)
       call cudaMemcpy_h2d_fp64(method%b_lo,method%C_b_lo,n,offset_zero,offset_zero)
    endif

    if(allocated(method%Dinv))then
!       n = size(method%Dinv)
       n = method%prm%n
!       call view_vec(n,method%Dinv)
       call cudaMalloc_fp64(method%C_Dinv,n)
       call cudaMemcpy_h2d_fp64(method%Dinv,method%C_Dinv,n,offset_zero,offset_zero)
    endif

    if(allocated(method%Dinv_hi))then
!       n = size(method%Dinv_hi)
       n = method%prm%n
       call cudaMalloc_fp64(method%C_Dinv_hi,n)
       call cudaMemcpy_h2d_fp64(method%Dinv_hi,method%C_Dinv_hi,n,offset_zero,offset_zero)
    endif

    if(allocated(method%Dinv_lo))then
!       n = size(method%Dinv_lo)
       n = method%prm%n
       call cudaMalloc_fp64(method%C_Dinv_lo,n)
       call cudaMemcpy_h2d_fp64(method%Dinv_lo,method%C_Dinv_lo,n,offset_zero,offset_zero)
    endif

  end subroutine set_VEC_GPU

  subroutine set_COMM_GPU(method)
    use krylov
    implicit none
    type(KSP) :: method
    integer i
    integer(kind=intType) n,nn
    integer MatFormat
    integer(kind=intType),parameter :: offset_zero=0
    MatFormat = method%prm%MatFormat 

    ! --------
    ! communication
    !     integer(kind=intType),allocatable :: istart(:)
    if(allocated(method%commVec%istart))then
!       n = size(method%commVec%istart)
       n = method%commVec%npes+1
       call cudaMalloc_PARCEL_INT(method%commVec%C_istart,n)
       call cudaMemcpy_h2d_PARCEL_INT(method%commVec%istart,method%commVec%C_istart,n,offset_zero,offset_zero)
    endif

    ! type(nrk),allocatable :: nrk_prm(:) 
    if(allocated(method%commVec%nrk_prm))then
!       nn = size(method%commVec%nrk_prm)
       nn = method%commVec%nrk_rank
       do i = 1, nn
          if(allocated(method%commVec%nrk_prm(i)%offset))then
!             n = size(method%commVec%nrk_prm(i)%offset)
             n = method%commVec%nrk_prm(i)%ndata
             call cudaMalloc_PARCEL_INT(method%commVec%nrk_prm(i)%C_offset,n)
             call cudaMemcpy_h2d_PARCEL_INT(method%commVec%nrk_prm(i)%offset,method%commVec%nrk_prm(i)%C_offset,n,offset_zero,offset_zero)
          endif
          if(allocated(method%commVec%nrk_prm(i)%nrk_ptr))then
!             n = size(method%commVec%nrk_prm(i)%nrk_ptr)
             n = method%prm%n + 1
             call cudaMalloc_PARCEL_INT(method%commVec%nrk_prm(i)%C_nrk_ptr,n)
             call cudaMemcpy_h2d_PARCEL_INT(method%commVec%nrk_prm(i)%nrk_ptr,method%commVec%nrk_prm(i)%C_nrk_ptr,n,offset_zero,offset_zero)
          endif
          if(allocated(method%commVec%nrk_prm(i)%nrkRow))then
!             n = size(method%commVec%nrk_prm(i)%nrkRow)
             n = method%commVec%nrk_prm(i)%nnz
             call cudaMalloc_PARCEL_INT(method%commVec%nrk_prm(i)%C_nrkRow,n)
             call cudaMemcpy_h2d_PARCEL_INT(method%commVec%nrk_prm(i)%nrkRow,method%commVec%nrk_prm(i)%C_nrkRow,n,offset_zero,offset_zero)
          endif
          if(allocated(method%commVec%nrk_prm(i)%nrkCol))then
!             n = size(method%commVec%nrk_prm(i)%nrkCol)
             n = method%commVec%nrk_prm(i)%nnz
             call cudaMalloc_PARCEL_INT(method%commVec%nrk_prm(i)%C_nrkCol,n)
             call cudaMemcpy_h2d_PARCEL_INT(method%commVec%nrk_prm(i)%nrkCol,method%commVec%nrk_prm(i)%C_nrkCol,n,offset_zero,offset_zero)
          endif

          ! if(allocated(method%commVec%nrk_prm(i)%@))then
          !    n = size(method%commVec%nrk_prm(i)%@)
          !    call cudaMalloc_int(method%commVec%nrk_prm(i)%C_@,n)
          !    call cudaMemcpy_h2d_int(method%commVec%nrk_prm(i)%@,method%commVec%nrk_prm(i)%C_@,n)
          ! endif

       enddo
    endif

    !     type(nrk),allocatable :: nrk_prm_send(:)
    if(allocated(method%commVec%nrk_prm_send))then
!       nn = size(method%commVec%nrk_prm_send)
       nn = method%commVec%nrk_rank_send
       do i = 1, nn
          if(allocated(method%commVec%nrk_prm_send(i)%offset))then
!             n = size(method%commVec%nrk_prm_send(i)%offset)
             n = method%commVec%nrk_prm_send(i)%ndata
             call cudaMalloc_PARCEL_INT(method%commVec%nrk_prm_send(i)%C_offset,n)
             call cudaMemcpy_h2d_PARCEL_INT(method%commVec%nrk_prm_send(i)%offset,method%commVec%nrk_prm_send(i)%C_offset,n,offset_zero,offset_zero)
          endif
          if(allocated(method%commVec%nrk_prm_send(i)%nrk_ptr))then
!             n = size(method%commVec%nrk_prm_send(i)%nrk_ptr)
             n =  method%prm%n + 1
             call cudaMalloc_PARCEL_INT(method%commVec%nrk_prm_send(i)%C_nrk_ptr,n)
             call cudaMemcpy_h2d_PARCEL_INT(method%commVec%nrk_prm_send(i)%nrk_ptr,method%commVec%nrk_prm_send(i)%C_nrk_ptr,n,offset_zero,offset_zero)
          endif
          if(allocated(method%commVec%nrk_prm_send(i)%nrkRow))then
!             n = size(method%commVec%nrk_prm_send(i)%nrkRow)
             n = method%commVec%nrk_prm(i)%nnz
            call cudaMalloc_PARCEL_INT(method%commVec%nrk_prm_send(i)%C_nrkRow,n)
             call cudaMemcpy_h2d_PARCEL_INT(method%commVec%nrk_prm_send(i)%nrkRow,method%commVec%nrk_prm_send(i)%C_nrkRow,n,offset_zero,offset_zero)
          endif
          if(allocated(method%commVec%nrk_prm_send(i)%nrkCol))then
!             n = size(method%commVec%nrk_prm_send(i)%nrkCol)
             n = method%commVec%nrk_prm(i)%nnz
             call cudaMalloc_PARCEL_INT(method%commVec%nrk_prm_send(i)%C_nrkCol,n)
             call cudaMemcpy_h2d_PARCEL_INT(method%commVec%nrk_prm_send(i)%nrkCol,method%commVec%nrk_prm_send(i)%C_nrkCol,n,offset_zero,offset_zero)
          endif

          ! if(allocated(method%commVec%nrk_prm_send(i)%@))then
          !    n = size(method%commVec%nrk_prm_send(i)%@)
          !    call cudaMalloc_int(method%commVec%nrk_prm_send(i)%C_@,n)
          !    call cudaMemcpy_h2d_int(method%commVec%nrk_prm_send(i)%@,method%commVec%nrk_prm_send(i)%C_@,n)
          ! endif

       enddo
    endif

    call malloc_MPI_Comm(method%commVec%C_comm,1) 
    call MPI_Comm_f2c(method%commVec%C_comm,method%commVec%comm) 

    if(method%commVec%nrk_rank_send>0)then
       call malloc_MPI_Status( &
            method%commVec%C_s_sta,method%commVec%nrk_rank_send) 
       call malloc_MPI_Request( &
            method%commVec%C_s_req,method%commVec%nrk_rank_send)
    endif
    if(method%commVec%nrk_rank>0)then
       call malloc_MPI_Status( &
            method%commVec%C_r_sta,method%commVec%nrk_rank) 
       call malloc_MPI_Request( &
            method%commVec%C_r_req,method%commVec%nrk_rank)
    endif

    ! --------
    ! MPI buffer
    ! type(devicemem) :: C_buff
    ! type(devicemem) :: C_buff_hi
    ! type(devicemem) :: C_buff_lo
    ! type(devicemem) :: C_rbuffL
    ! type(devicemem) :: C_rbuffL_hi
    ! type(devicemem) :: C_rbuffL_lo
    ! type(devicemem) :: C_rbuffU
    ! type(devicemem) :: C_rbuffU_hi
    ! type(devicemem) :: C_rbuffU_lo
    ! if(allocated(method%@))then
    !    n = size(method%@)
    !    call cudaMalloc_fp64(method%C_@,n)
    !    call cudaMemcpy_h2d_fp64(method%@,method%C_@,n)
    ! endif

    n = max(method%commVec%nrk_data_send,1)
!    if(n>0)then
       call cudaMalloc_fp64(method%C_sbuff,n)
!    endif
    if(allocated(method%buff))then
!       n = size(method%buff)
       n = method%commVec%nrk_data+1
!       write(*,*)"n, n_ = ",size(method%buff),method%commVec%nrk_data+1 
       call cudaMalloc_fp64(method%C_buff,n)
       call cudaMemcpy_h2d_fp64(method%buff,method%C_buff,n,offset_zero,offset_zero)
    endif
    if(allocated(method%buff_hi))then
       n = size(method%buff_hi)
!       n = method%commVec%nrk_data+1
       call cudaMalloc_fp64(method%C_buff_hi,n)
       call cudaMemcpy_h2d_fp64(method%buff_hi,method%C_buff_hi,n,offset_zero,offset_zero)
    endif
    if(allocated(method%buff_lo))then
       n = size(method%buff_lo)
!       n = method%commVec%nrk_data+1
       call cudaMalloc_fp64(method%C_buff_lo,n)
       call cudaMemcpy_h2d_fp64(method%buff_lo,method%C_buff_lo,n,offset_zero,offset_zero)
    endif

  end subroutine set_COMM_GPU

  subroutine set_DIA_Matrix_GPU(method,diaMat,C_diaMat)
    use krylov
    implicit none
    type(KSP) :: method
    type(DIA_Mat)    :: diaMat 
    type(C_DIA_Mat)  :: C_diaMat 
    integer(kind=intType) n
    integer(kind=intType) nnd
    integer(kind=intType),parameter :: offset_zero=0

    nnd = diaMat%nnd
    C_diaMat%nnd   = diaMat%nnd
    C_diaMat%nCol  = diaMat%nCol
    C_diaMat%nColL = diaMat%nColL
    C_diaMat%nColU = diaMat%nColU
    C_diaMat%nColD = diaMat%nColD

    call cuda_MatVec_dia_init_GPU(nnd,diaMat%offset)

    if(allocated(diaMat%diaA))then
!       n = size(diaMat%diaA)
       n = diaMat%nnd * method%prm%n
       call cudaMalloc_fp64(C_diaMat%diaA,n)
       call cudaMemcpy_h2d_fp64(diaMat%diaA,C_diaMat%diaA,n,offset_zero,offset_zero)
    endif

    if(allocated(diaMat%diaA_hi))then
!       n = size(diaMat%diaA_hi)
       n = diaMat%nnd * method%prm%n
       call cudaMalloc_fp64(C_diaMat%diaA_hi,n)
       call cudaMemcpy_h2d_fp64(diaMat%diaA_hi,C_diaMat%diaA_hi,n,offset_zero,offset_zero)
    endif

    if(allocated(diaMat%diaA_lo))then
!       n = size(diaMat%diaA_lo)
       n = diaMat%nnd * method%prm%n
       call cudaMalloc_fp64(C_diaMat%diaA_lo,n)
       call cudaMemcpy_h2d_fp64(diaMat%diaA_lo,C_diaMat%diaA_lo,n,offset_zero,offset_zero)
    endif

    if(allocated(diaMat%offset))then
!       n = size(diaMat%offset)
       n = diaMat%nnd
       call cudaMalloc_PARCEL_INT(C_diaMat%offset,n)
       call cudaMemcpy_h2d_PARCEL_INT(diaMat%offset,C_diaMat%offset,n,offset_zero,offset_zero)
    endif

    if(allocated(diaMat%iLs))then
!       n = size(diaMat%iLs)
       n = diaMat%nnd * method%prm%n
       call cudaMalloc_PARCEL_INT(C_diaMat%iLs,n)
       call cudaMemcpy_h2d_PARCEL_INT(diaMat%iLs,C_diaMat%iLs,n,offset_zero,offset_zero)
    endif

    if(allocated(diaMat%offsetVec))then
!       n = size(diaMat%offsetVec)
       n =  diaMat%nnd 
       call cudaMalloc_PARCEL_INT(C_diaMat%offsetVec,n)
       call cudaMemcpy_h2d_PARCEL_INT(diaMat%offsetVec,C_diaMat%offsetVec,n,offset_zero,offset_zero)
    endif

    if(allocated(diaMat%Vec))then
!       n = size(diaMat%Vec)
       n = method%diaMat%nCol
       call cudaMalloc_fp64(C_diaMat%Vec,n)
       call cudaMemcpy_h2d_fp64(diaMat%Vec,C_diaMat%Vec,n,offset_zero,offset_zero)
    endif

    if(allocated(diaMat%Vec_hi))then
!       n = size(diaMat%Vec_hi)
       n = method%diaMat%nCol
       call cudaMalloc_fp64(C_diaMat%Vec_hi,n)
       call cudaMemcpy_h2d_fp64(diaMat%Vec_hi,C_diaMat%Vec_hi,n,offset_zero,offset_zero)
    endif

    if(allocated(diaMat%Vec_lo))then
!       n = size(diaMat%Vec_lo)
       n = method%diaMat%nCol
       call cudaMalloc_fp64(C_diaMat%Vec_lo,n)
       call cudaMemcpy_h2d_fp64(diaMat%Vec_lo,C_diaMat%Vec_lo,n,offset_zero,offset_zero)
    endif

    if(allocated(diaMat%RowL))then
!       n = size(diaMat%RowL)
       n = diaMat%nColL
       call cudaMalloc_PARCEL_INT(C_diaMat%RowL,n)
       call cudaMemcpy_h2d_PARCEL_INT(diaMat%RowL,C_diaMat%RowL,n,offset_zero,offset_zero)
    endif

    if(allocated(diaMat%RowD))then
!       n = size(diaMat%RowD)
       n = diaMat%nColD
       call cudaMalloc_PARCEL_INT(C_diaMat%RowD,n)
       call cudaMemcpy_h2d_PARCEL_INT(diaMat%RowD,C_diaMat%RowD,n,offset_zero,offset_zero)
    endif

    if(allocated(diaMat%RowU))then
!       n = size(diaMat%RowU)
       n = diaMat%nColU
       call cudaMalloc_PARCEL_INT(C_diaMat%RowU,n)
       call cudaMemcpy_h2d_PARCEL_INT(diaMat%RowU,C_diaMat%RowU,n,offset_zero,offset_zero)
    endif

    ! if(allocated(diaMat%@))then
    !    n = size(diaMat%@)
    !    call cudaMalloc_PARCEL_INT(C_diaMat%@,n)
    !    call cudaMemcpy_h2d_PARCEL_INT(diaMat%@,C_diaMat%@,n)
    ! endif
  end subroutine set_DIA_Matrix_GPU

#if 0
  subroutine set_precon_DILU(method)
    use krylov
    implicit none
    type(KSP) :: method
    integer MatFormat
    integer(kind=intType) n
    MatFormat = method%prm%MatFormat 
    if(MatFormat.eq.MATRIX_FORMAT_CRS)then
       if(allocated(method%preconILU))then
@          n = size(method%preconILU)
!          call cudaMalloc_PARCEL_INT(C_diaMat%@,n)
!          call cudaMemcpy_h2d_PARCEL_INT(diaMat%@,C_diaMat%@,n)
       endif
    endif

  end subroutine set_precon_DILU
#endif

  subroutine set_precon_Neumann_GPU(method)
    use krylov
    implicit none
    type(KSP) :: method
    integer MatFormat
    integer(kind=intType) n
    integer(kind=intType) vecLen
    integer(kind=intType) nColL
    integer(kind=intType) nColU
    integer(kind=intType),parameter :: offset_zero=0

    ! method%C_preconNeuman%itr = 1
    ! method%C_preconNeuman%local_flag = KSP_PRECONDITION_NEUMANN_TRUE_COMM 

    nColL = 0
    nColU = 0

    MatFormat = method%prm%MatFormat 
    ! dia形式ではcuSPARCEを利用せずに通信用バッファでSpMVを計算するためベクトルに通信用のバッファが不要
    if(MatFormat.eq.MATRIX_FORMAT_CRS)then
       nColL  = method%crsMat%nColL
       nColU  = method%crsMat%nColU
    endif
    vecLen = method%prm%n+nColL+nColU

    if(MatFormat.eq.MATRIX_FORMAT_CRS)then
       call set_CRS_Matrix_GPU( &
            method,&
            method%preconNeuman%crsMat,&
            method%C_preconNeuman%crsMat)
    endif
    if(MatFormat.eq.MATRIX_FORMAT_DIA)then
       call set_DIA_Matrix_GPU( &
            method, &
            method%preconNeuman%diaMat,&
            method%C_preconNeuman%diaMat)
    endif

    n = vecLen
    call cudaMalloc_fp64(method%C_preconNeuman%vecx,n)
    call cusparseCreateDnVec_fp64(method%C_preconNeuman%vecx,method%prm%n,vecLen,offset_zero)

    n = vecLen
    call cudaMalloc_fp64(method%C_preconNeuman%vecx_hi,n)
    call cusparseCreateDnVec_fp64(method%C_preconNeuman%vecx_hi,method%prm%n,vecLen,offset_zero)

    n = vecLen
    call cudaMalloc_fp64(method%C_preconNeuman%vecx_lo,n)
    call cusparseCreateDnVec_fp64(method%C_preconNeuman%vecx_lo,method%prm%n,vecLen,offset_zero)

    n = vecLen
    call cudaMalloc_fp64(method%C_preconNeuman%vecy,n)
    call cusparseCreateDnVec_fp64(method%C_preconNeuman%vecy,method%prm%n,vecLen,offset_zero)

    n = vecLen
    call cudaMalloc_fp64(method%C_preconNeuman%vecy_hi,n)
    call cusparseCreateDnVec_fp64(method%C_preconNeuman%vecy_hi,method%prm%n,vecLen,offset_zero)

    n = vecLen
    call cudaMalloc_fp64(method%C_preconNeuman%vecy_lo,n)
    call cusparseCreateDnVec_fp64(method%C_preconNeuman%vecy_lo,method%prm%n,vecLen,offset_zero)

    n = vecLen
    call cudaMalloc_fp64(method%C_preconNeuman%D,n)
    if(allocated(method%preconNeuman%D))then
       call cudaMemcpy_h2d_fp64(method%preconNeuman%D,method%C_preconNeuman%D,method%prm%n,offset_zero,offset_zero)
    endif
    call cusparseCreateDnVec_fp64(method%C_preconNeuman%D,method%prm%n,vecLen,offset_zero)

    n = vecLen
    call cudaMalloc_fp64(method%C_preconNeuman%D_hi,n)
    if(allocated(method%preconNeuman%D_hi))then
       call cudaMemcpy_h2d_fp64(method%preconNeuman%D_hi,method%C_preconNeuman%D_hi,method%prm%n,offset_zero,offset_zero)
    endif
    call cusparseCreateDnVec_fp64(method%C_preconNeuman%D_hi,method%prm%n,vecLen,offset_zero)

    n = vecLen
    call cudaMalloc_fp64(method%C_preconNeuman%D_lo,n)
    if(allocated(method%preconNeuman%D_lo))then
       call cudaMemcpy_h2d_fp64(method%preconNeuman%D_lo,method%C_preconNeuman%D_lo,method%prm%n,offset_zero,offset_zero)
    endif
    call cusparseCreateDnVec_fp64(method%C_preconNeuman%D_lo,method%prm%n,vecLen,offset_zero)

    n = vecLen
    call cudaMalloc_fp64(method%C_preconNeuman%Dinv,n)
    if(allocated(method%preconNeuman%Dinv))then
       call cudaMemcpy_h2d_fp64(method%preconNeuman%Dinv,method%C_preconNeuman%Dinv,method%prm%n,offset_zero,offset_zero)
    endif
    call cusparseCreateDnVec_fp64(method%C_preconNeuman%Dinv,method%prm%n,vecLen,offset_zero)

    n = vecLen
    call cudaMalloc_fp64(method%C_preconNeuman%Dinv_hi,n)
    if(allocated(method%preconNeuman%Dinv_hi))then
       call cudaMemcpy_h2d_fp64(method%preconNeuman%Dinv_hi,method%C_preconNeuman%Dinv_hi,method%prm%n,offset_zero,offset_zero)
    endif
    call cusparseCreateDnVec_fp64(method%C_preconNeuman%Dinv_hi,method%prm%n,vecLen,offset_zero)

    n = vecLen
    call cudaMalloc_fp64(method%C_preconNeuman%Dinv_lo,n)
    if(allocated(method%preconNeuman%Dinv_lo))then
       call cudaMemcpy_h2d_fp64(method%preconNeuman%Dinv_lo,method%C_preconNeuman%Dinv_lo,method%prm%n,offset_zero,offset_zero)
    endif
    call cusparseCreateDnVec_fp64(method%C_preconNeuman%Dinv_lo,method%prm%n,vecLen,offset_zero)

  end subroutine set_precon_Neumann_GPU

  subroutine set_CRS_Matrix_GPU(method,crsMat,C_crsMat)
    use krylov
    implicit none
    type(KSP) :: method
    type(CRS_Mat)    :: crsMat 
    type(C_CRS_Mat)  :: C_crsMat 
    integer(kind=intType) k
    integer(kind=intType) n
    integer(kind=intType) work_int_size
    integer(kind=intType),parameter :: parcel_int_zero=0
    integer,allocatable :: work_int(:)

    C_crsMat%nnz   = crsMat%nnz
    C_crsMat%nColL = crsMat%nColL 
    C_crsMat%nColU = crsMat%nColU 

    ! write(*,*)" @ 1 @ "
    if(allocated(crsMat%crsA))then
!       n = size(crsMat%crsA)
       n = crsMat%nnz
       call cudaMalloc_fp64(C_crsMat%crsA,n)
       call cudaMemcpy_h2d_fp64(crsMat%crsA,C_crsMat%crsA,n,parcel_int_zero,parcel_int_zero)
    endif

    ! write(*,*)" @ 2 @ "
    if(allocated(crsMat%crsA_hi))then
!       n = size(crsMat%crsA_hi)
       n = crsMat%nnz
       call cudaMalloc_fp64(C_crsMat%crsA_hi,n)
       call cudaMemcpy_h2d_fp64(crsMat%crsA_hi,C_crsMat%crsA_hi,n,parcel_int_zero,parcel_int_zero)
    endif
    if(allocated(crsMat%crsA_lo))then
!       n = size(crsMat%crsA_lo)
       n = crsMat%nnz
       call cudaMalloc_fp64(C_crsMat%crsA_lo,n)
       call cudaMemcpy_h2d_fp64(crsMat%crsA_lo,C_crsMat%crsA_lo,n,parcel_int_zero,parcel_int_zero)
    endif
! 
    ! write(*,*)" @ 3 @ "
    if(allocated(crsMat%crsRow_ptr))then
!       n = size(crsMat%crsRow_ptr)
       n = method%prm%n + 1
       call cudaMalloc_PARCEL_int(C_crsMat%crsRow_ptr,n)
       call cudaMemcpy_h2d_PARCEL_INT(crsMat%crsRow_ptr,C_crsMat%crsRow_ptr,n,parcel_int_zero,parcel_int_zero)
    endif

    ! write(*,*)"      @ 4 @ "
    if(allocated(crsMat%crsRow_ptr_Diag_s))then
!       n = size(crsMat%crsRow_ptr_Diag_s)
       n = method%prm%n 
       call cudaMalloc_PARCEL_INT(C_crsMat%crsRow_ptr_Diag_s,n)
       call cudaMemcpy_h2d_PARCEL_INT(crsMat%crsRow_ptr_Diag_s,C_crsMat%crsRow_ptr_Diag_s,n,parcel_int_zero,parcel_int_zero)
    endif
    if(allocated(crsMat%crsRow_ptr_Diag_e))then
!       n = size(crsMat%crsRow_ptr_Diag_e)
       n = method%prm%n 
       call cudaMalloc_PARCEL_INT(C_crsMat%crsRow_ptr_Diag_e,n)
       call cudaMemcpy_h2d_PARCEL_INT(crsMat%crsRow_ptr_Diag_e,C_crsMat%crsRow_ptr_Diag_e,n,parcel_int_zero,parcel_int_zero)
    endif

    ! write(*,*)" @ 5 @ "
    if(allocated(crsMat%crsCol))then
       n = crsMat%nnz
       call cudaMalloc_PARCEL_INT(C_crsMat%crsCol,n)
       call cudaMemcpy_h2d_PARCEL_INT(crsMat%crsCol,C_crsMat%crsCol,n,parcel_int_zero,parcel_int_zero)
    endif

    ! write(*,*)" @ 6 @ "
    if(allocated(crsMat%RowL))then
!       n = size(crsMat%RowL)
       n = crsMat%nColL
       call cudaMalloc_PARCEL_INT(C_crsMat%RowL,n)
       call cudaMemcpy_h2d_PARCEL_INT(crsMat%RowL,C_crsMat%RowL,n,parcel_int_zero,parcel_int_zero)
    endif
    if(allocated(crsMat%RowD))then
!       n = size(crsMat%RowD)
       n = method%prm%n 
       call cudaMalloc_PARCEL_INT(C_crsMat%RowD,n)
       call cudaMemcpy_h2d_PARCEL_INT(crsMat%RowD,C_crsMat%RowD,n,parcel_int_zero,parcel_int_zero)
    endif
    if(allocated(crsMat%RowU))then
       !       n = size(crsMat%RowU)
       n = crsMat%nColU
       call cudaMalloc_PARCEL_INT(C_crsMat%RowU,n)
       call cudaMemcpy_h2d_PARCEL_INT(crsMat%RowU,C_crsMat%RowU,n,parcel_int_zero,parcel_int_zero)
    endif

    ! write(*,*)" @ 7 @ "
    if(allocated(crsMat%iLs))then
!       n = size(crsMat%iLs)
       n = crsMat%nnz
       call cudaMalloc_PARCEL_INT(C_crsMat%iLs,n)
       call cudaMemcpy_h2d_PARCEL_INT(crsMat%iLs,C_crsMat%iLs,n,parcel_int_zero,parcel_int_zero)
    endif
  end subroutine set_CRS_Matrix_GPU

  subroutine set_CRS_Matrix_FP64toBF16_GPU(&
       method, &
       n, &
       crsMat,&
       C_crsMat_FP64,&
       C_crsMat_BF16,&
       mode,mode_diag)
    use krylov
    use cu_Matrix,only : Diag_Transfer_FP64toBF16_crs_FP64_GPU
    implicit none
    type(KSP) :: method
    type(CRS_Mat)    :: crsMat 
    type(C_CRS_Mat)  :: C_crsMat_FP64
    type(C_CRS_Mat)  :: C_crsMat_BF16 
    integer mode,mode_diag
    integer(kind=intType) :: n
    integer(kind=intType) nn
    integer(kind=intType),parameter :: offset_zero=0

    call set_CRS_Matrix_GPU( &
         method, &
         crsMat,C_crsMat_BF16)
!    nn = size(crsMat%crsA)
    nn = crsMat%nnz
    call cudaFree(C_crsMat_BF16%crsA)
    call cudaMalloc_bf16(C_crsMat_BF16%crsA,nn)
    call cudaTransfer_h2d_fp64toBF16( &
         crsMat%crsA,C_crsMat_BF16%crsA,nn,mode)

    call Diag_Transfer_FP64toBF16_crs_FP64_GPU( &
         CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
          n,&
          C_crsMat_BF16%crsA, &
          C_crsMat_FP64%crsA, &
          C_crsMat_FP64%crsRow_ptr, &
          C_crsMat_FP64%crsCol, &
          C_crsMat_FP64%iLs, &
          mode_diag &
          )

  end subroutine set_CRS_Matrix_FP64toBF16_GPU

  subroutine free_CRS_Matrix_GPU(crsMat,C_crsMat)
    use krylov
    implicit none
    type(CRS_Mat)    :: crsMat 
    type(C_CRS_Mat)  :: C_crsMat 

    ! --------
    ! CRS Matrix
    ! if(allocated(crsMat%crsA))then
       call cudaFree(C_crsMat%crsA)
    ! endif
    ! if(allocated(crsMat%crsA_hi))then
       call cudaFree(C_crsMat%crsA_hi)
    ! endif
    ! if(allocated(crsMat%crsA_lo))then
       call cudaFree(C_crsMat%crsA_lo)
    ! endif
    ! if(allocated(crsMat%crsRow_ptr))then
       call cudaFree(C_crsMat%crsRow_ptr)
    ! endif
    ! if(allocated(crsMat%crsRow_ptr_Diag_s))then
       call cudaFree(C_crsMat%crsRow_ptr_Diag_s)
    ! endif
    ! if(allocated(crsMat%crsRow_ptr_Diag_e))then
       call cudaFree(C_crsMat%crsRow_ptr_Diag_e)
    ! endif
    ! if(allocated(crsMat%crsCol))then
       call cudaFree(C_crsMat%crsCol)
    ! endif
    ! if(allocated(crsMat%RowL))then
       call cudaFree(C_crsMat%RowL)
    ! endif
    ! if(allocated(crsMat%RowD))then
       call cudaFree(C_crsMat%RowD)
    ! endif
    ! if(allocated(crsMat%RowU))then
       call cudaFree(C_crsMat%RowU)
    ! endif
    ! if(allocated(crsMat%iLs))then
       call cudaFree(C_crsMat%iLs)
    ! endif

  end subroutine free_CRS_Matrix_GPU

  subroutine free_DIA_Matrix_GPU(diaMat,C_diaMat)
    use krylov
    implicit none
    type(DIA_Mat)    :: diaMat 
    type(C_DIA_Mat)  :: C_diaMat 

    ! --------
    ! DIA Matrix
    ! if(allocated(diaMat%diaA))then
       call cudaFree(C_diaMat%diaA)
    ! endif
    ! if(allocated(diaMat%diaA_hi))then
       call cudaFree(C_diaMat%diaA_hi)
    ! endif
    ! if(allocated(diaMat%diaA_lo))then
       call cudaFree(C_diaMat%diaA_lo)
    ! endif
    ! if(allocated(diaMat%offset))then
       call cudaFree(C_diaMat%offset)
    ! endif
    ! if(allocated(diaMat%iLs))then
       call cudaFree(C_diaMat%iLs)
    ! endif
    ! if(allocated(diaMat%offsetVec))then
       call cudaFree(C_diaMat%offsetVec)
    ! endif
    ! if(allocated(diaMat%Vec))then
       call cudaFree(C_diaMat%Vec)
    ! endif
    ! if(allocated(diaMat%Vec_hi))then
       call cudaFree(C_diaMat%Vec_hi)
    ! endif
    ! if(allocated(diaMat%Vec_lo))then
       call cudaFree(C_diaMat%Vec_lo)
    ! endif
    ! if(allocated(diaMat%RowL))then
       call cudaFree(C_diaMat%RowL)
    ! endif
    ! if(allocated(diaMat%RowD))then
       call cudaFree(C_diaMat%RowD)
    ! endif
    ! if(allocated(diaMat%RowU))then
       call cudaFree(C_diaMat%RowU)
    ! endif

  end subroutine free_DIA_Matrix_GPU

  subroutine free_VEC_GPU(method)
    use krylov
    implicit none
    type(KSP) :: method
    ! --------
    ! Vector
    if(allocated(method%r))then
       call cudaFree(method%C_r)
    endif
    if(allocated(method%r_hi))then
       call cudaFree(method%C_r_hi)
    endif
    if(allocated(method%r_lo))then
       call cudaFree(method%C_r_lo)
    endif

    if(allocated(method%x))then
       call cudaFree(method%C_x)
    endif
    if(allocated(method%x_hi))then
       call cudaFree(method%C_x_hi)
    endif
    if(allocated(method%x_lo))then
       call cudaFree(method%C_x_lo)
    endif

    if(allocated(method%b))then
       call cudaFree(method%C_b)
    endif
    if(allocated(method%b_hi))then
       call cudaFree(method%C_b_hi)
    endif
    if(allocated(method%b_lo))then
       call cudaFree(method%C_b_lo)
    endif

    if(allocated(method%Dinv))then
       call cudaFree(method%C_Dinv)
    endif

    if(allocated(method%Dinv_hi))then
       call cudaFree(method%C_Dinv_hi)
    endif

    if(allocated(method%Dinv_lo))then
       call cudaFree(method%C_Dinv_lo)
    endif

  end subroutine free_VEC_GPU

  subroutine free_precon_Neumann_GPU(method)
    use krylov
    implicit none
    type(KSP) :: method
    integer MatFormat
    MatFormat = method%prm%MatFormat 

    if(MatFormat.eq.MATRIX_FORMAT_CRS)then
       call free_CRS_Matrix_GPU(&
            method%preconNeuman%crsMat, &
            method%C_preconNeuman%crsMat)
    endif
    if(MatFormat.eq.MATRIX_FORMAT_DIA)then
       call free_DIA_Matrix_GPU( &
            method%preconNeuman%diaMat,&
            method%C_preconNeuman%diaMat)
    endif

    ! if(allocated(method%preconNeuman%vecx))then
    call cudaFree(method%C_preconNeuman%vecx)
    call cusparseDestroyDnVec(method%C_preconNeuman%vecx)
    ! endif
    ! if(allocated(method%preconNeuman%vecx_hi))then
    call cudaFree(method%C_preconNeuman%vecx_hi)
    call cusparseDestroyDnVec(method%C_preconNeuman%vecx_hi)
    ! endif
    ! if(allocated(method%preconNeuman%vecx_lo))then
    call cudaFree(method%C_preconNeuman%vecx_lo)
    call cusparseDestroyDnVec(method%C_preconNeuman%vecx_lo)
    ! endif

    ! if(allocated(method%preconNeuman%vecy))then
    call cudaFree(method%C_preconNeuman%vecy)
    call cusparseDestroyDnVec(method%C_preconNeuman%vecy)
    ! endif
    ! if(allocated(method%preconNeuman%vecy_hi))then
    call cudaFree(method%C_preconNeuman%vecy_hi)
    call cusparseDestroyDnVec(method%C_preconNeuman%vecy_hi)
    ! endif
    ! if(allocated(method%preconNeuman%vecy_lo))then
    call cudaFree(method%C_preconNeuman%vecy_lo)
    call cusparseDestroyDnVec(method%C_preconNeuman%vecy_lo)
    ! endif

    ! if(allocated(method%preconNeuman%D))then
    call cudaFree(method%C_preconNeuman%D)
    call cusparseDestroyDnVec(method%C_preconNeuman%D)
    ! endif
    ! if(allocated(method%preconNeuman%D_hi))then
    call cudaFree(method%C_preconNeuman%D_hi)
    call cusparseDestroyDnVec(method%C_preconNeuman%D_hi)
    ! endif
    ! if(allocated(method%preconNeuman%D_lo))then
    call cudaFree(method%C_preconNeuman%D_lo)
    call cusparseDestroyDnVec(method%C_preconNeuman%D_lo)
    ! endif

    ! if(allocated(method%preconNeuman%Dinv))then
    call cudaFree(method%C_preconNeuman%Dinv)
    call cusparseDestroyDnVec(method%C_preconNeuman%Dinv)
    ! endif
    ! if(allocated(method%preconNeuman%Dinv_hi))then
    call cudaFree(method%C_preconNeuman%Dinv_hi)
    call cusparseDestroyDnVec(method%C_preconNeuman%Dinv_hi)
    ! endif
    ! if(allocated(method%preconNeuman%Dinv_lo))then
    call cudaFree(method%C_preconNeuman%Dinv_lo)
    call cusparseDestroyDnVec(method%C_preconNeuman%Dinv_lo)
    ! endif

  end subroutine free_precon_Neumann_GPU
 
  ! init , free ともにベクトル側のメモリの指定を行う
  subroutine init_GPU(method)
    use krylov
    use matrix_kernel
    implicit none
    type(KSP) :: method
    integer(kind=intType) i
    integer(kind=intType) n
    integer(kind=intType) nn
    integer(kind=intType) nnd
    integer myrank
    integer (kind=intType) istart
    integer (kind=intType) nnz,nnz_
    integer (kind=intType) nCol,nColL,nColU,nColD
    type(Local_CRS_Mat) :: LocalcrsMat
    type(CRS_Mat) :: crsMat ! DIA形式からのCRS形式の変換で利用 
    integer(kind=intType),allocatable :: Colsort(:)
    integer MatFormat
    real*8 stime,etime

    MatFormat = method%prm%MatFormat 

    ! --------
    ! GPU Setting
    myrank=method%commVec%myrank

    stime=MPI_Wtime()
!    call set_GPU(myrank)
    etime=MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"     @ 4-1 @",etime-stime
    endif
#endif

    stime=MPI_Wtime()
    call set_cublas(method%cublas_solver)
    etime=MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"     @ 4-2 @",etime-stime
    endif
#endif

    ! --------
    ! CRS Matrix
    stime=MPI_Wtime()
    if(MatFormat.eq.MATRIX_FORMAT_CRS)then
       call set_CRS_Matrix_GPU( &
            method, &
            method%crsMat,method%C_crsMat)
    endif

    etime=MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"     @ 4-3 @",etime-stime
    endif
#endif

    ! --------
    ! DIA Matrix
    stime=MPI_Wtime()
    if(MatFormat.eq.MATRIX_FORMAT_DIA)then
       call set_DIA_Matrix_GPU(&
                        method, &
                        method%diaMat,method%C_diaMat)
    endif
    etime=MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"     @ 4-4 @",etime-stime
    endif
#endif

    stime=MPI_Wtime()
    if(method%ipreflag.eq.KSP_PRECONDITION_NEUMANN)then
       call set_precon_Neumann_GPU(method)
    endif
    call set_VEC_GPU(method)

    ! C MPI    
    call set_COMM_GPU(method)
    etime=MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"     @ 4-5 @",etime-stime
    endif
#endif

    stime=MPI_Wtime()
    ! CUSPARCE initialize (SpMV.ILU)
    !    if(method%ipreflag.eq.KSP_PRECONDITION_BLOCKJACOBI)then
#if 1
    call set_CUSPARCE(method)
#else
    call set_CUSPARCE_GPUmem(method)
#endif
!   
    !    endif

    etime=MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"     @ 4-6 @",etime-stime
    endif
#endif
 
  end subroutine init_GPU

 subroutine set_CRS_Matrix_val_GPU( &
       crsMat, &
       crsMat_val_GPU, &
       n,istart,nnz &
       )
    use matrix
    use PARCEL_integer
    implicit none
    type(CRS_Mat) :: crsMat
    type(devicemem) :: crsMat_val_GPU
    integer(kind=intType) n
    integer(kind=intType) nnz
    real*8 crsA(nnz)
 
    integer(kind=intType) istart
    integer(kind=intType) DiagBlock_Row_ptr_start
    integer(kind=intType) DiagBlock_Row_ptr_end
    integer(kind=intType) k,j,j_,j0
    integer(kind=intType) iend
    integer(kind=intType),parameter :: offset_zero=0    
    real*8 stime,etime

    iend   = istart + n -1

    stime=MPI_Wtime()
    !$OMP PARALLEL DO PRIVATE(k,j,j0,j_,DiagBlock_Row_ptr_start,DiagBlock_Row_ptr_end )    
    do k=1,n
       DiagBlock_Row_ptr_start = crsMat%crsRow_ptr(k)
       DiagBlock_Row_ptr_end   = crsMat%crsRow_ptr(k+1)-1

       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          if(crsMat%crsCol(j)>=istart)then
             DiagBlock_Row_ptr_start = j
             exit
          endif
       enddo
       do j=crsMat%crsRow_ptr(k+1)-1,crsMat%crsRow_ptr(k),-1
          if(crsMat%crsCol(j)<=iend)then
             DiagBlock_Row_ptr_end = j
             exit
          endif
       enddo
       
       ! ------------
       j0 = 0
       ! 対角ブロック
       do j=DiagBlock_Row_ptr_start,DiagBlock_Row_ptr_end
          j_ = crsMat%crsRow_ptr(k) + j0
          j0 = j0 + 1
          crsA(j_)   = crsMat%crsA(j)
       enddo
       
       ! L ブロック
       do j=crsMat%crsRow_ptr(k),DiagBlock_Row_ptr_start-1
          j_ = crsMat%crsRow_ptr(k) + j0
          j0 = j0 + 1
          crsA(j_)   =  crsMat%crsA(j) 
       enddo

       ! U ブロック
       do j=DiagBlock_Row_ptr_end+1,crsMat%crsRow_ptr(k+1)-1
          j_ = crsMat%crsRow_ptr(k) + j0
          j0 = j0 + 1
          crsA(j_)   = crsMat%crsA(j) 
       enddo
       
    enddo
    !$OMP END PARALLEL DO    

    etime=MPI_Wtime()
#if PRINT_TIME_
   if(method%commVec%myrank.eq.0)then
       write(*,*)"@ reset - set_CRS_Matrix_val_GPU - CPU time @",etime-stime
   endif
#endif

    stime=MPI_Wtime()
    call cudaMemcpy_h2d_fp64(crsA,crsMat_val_GPU,nnz,offset_zero,offset_zero)
    etime=MPI_Wtime()
#if PRINT_TIME_
    if(method%commVec%myrank.eq.0)then
      write(*,*)"@ reset - set_CRS_Matrix_val_GPU - cudaMemcpy_h2d time @",etime-stime
    endif
#endif

  end subroutine set_CRS_Matrix_val_GPU

  subroutine reset_GPU(method,ipreflag0)
    use krylov
    use matrix_kernel
    implicit none
    type(KSP) :: method
    integer ipreflag
    integer ipreflag0
    integer(kind=intType) n
    integer (kind=intType) nnz
    integer (kind=intType) nCol
    integer MatFormat
    real*8 stime,etime
    ! GPUメモリを入力とする方式に未対応
#if 0
    MatFormat = method%prm%MatFormat 
    ipreflag = method%ipreflag
    ! call cusparse_free(method%cusparce_solver)
    call cudaDeviceSynchronize() 
    call cusparse_free(method%cusparce_solver)

    ! if( (ipreflag0.eq.KSP_PRECONDITION_BLOCKJACOBI) .or. &
    !      (ipreflag0.eq.KSP_PRECONDITION_ADDITIVESCHWARZ) &
    !      )then
    !    call cusparse_free(method%cusparce_solver)          
    ! endif

    ! write(*,*)"@ reset_GPU : debug 1 @"
#if 1
    call cudaFree(method%cusparce_crsMatrix%crsA)
    call cudaFree(method%cusparce_crsMatrix%crsRow_ptr)
    call cudaFree(method%cusparce_crsMatrix%crsCol)
    
    if(ipreflag0.eq.KSP_PRECONDITION_NEUMANN)then
       !       call cusparseDestroyDnVec(@); call cudaFree(@)
       if(ipreflag.ne.KSP_PRECONDITION_NEUMANN)then
          call free_precon_Neumann_GPU(method)
       endif

    else

       if(ipreflag0.eq.KSP_PRECONDITION_ADDITIVESCHWARZ)then
          call cusparseDestroyDnVec(method%C_ivec_AS)
          call cusparseDestroyDnVec(method%C_ovec_AS)
          call cudaFree(method%C_ivec_AS)
          call cudaFree(method%C_ovec_AS)
          call cudaFree(method%C_s_buffL)
          call cudaFree(method%C_r_buffL)
          call cudaFree(method%C_s_buffU)
          call cudaFree(method%C_r_buffU)
          
          call free(method%C_s_reqL)  
          call free(method%C_r_reqL)  
          call free(method%C_s_reqU)  
          call free(method%C_r_reqU)  
             
          ! MPIステータス
          call free(method%C_s_staL)   
          call free(method%C_r_staL)   
          call free(method%C_s_staU)   
          call free(method%C_r_staU) 
       endif  

    endif

    if(  &
         (ipreflag0.eq.KSP_PRECONDITION_NEUMANN) .and.  &
         (MatFormat.eq.MATRIX_FORMAT_DIA)  &
         )then

    else
       call cudaFree(method%cusparce_preconcrsMatrix%crsA)
       call cudaFree(method%cusparce_preconcrsMatrix%crsRow_ptr)
       call cudaFree(method%cusparce_preconcrsMatrix%crsCol)
    endif

#endif

    ! write(*,*)"@ reset_GPU : debug 2 @"
    if(MatFormat.eq.MATRIX_FORMAT_CRS)then
       stime=MPI_Wtime()
       if(allocated(method%crsMat%crsA))then
          n = method%crsMat%nnz
          call cudaMemcpy_h2d_fp64(method%crsMat%crsA,method%C_crsMat%crsA,n,offset_zero,offset_zero)
       endif

       if(allocated(method%Dinv))then
          n = method%prm%n
          call cudaMemcpy_h2d_fp64(&
               method%Dinv,&
               method%C_Dinv,n,offset_zero,offset_zero)
       endif

       etime=MPI_Wtime()
#if PRINT_TIME
       if(method%commVec%myrank.eq.0)then
          write(*,*)"@ reset - CRS - Matrix @",etime-stime
       endif
#endif

       stime=MPI_Wtime()
       ! call set_CRS_Matrix_val_GPU( &
       !      method%crsMat,&
       !      method%cusparce_crsMatrix%crsA, &
       !      method%prm%n,method%prm%istart,method%prm%nnz &
       !      )
       etime=MPI_Wtime()

#if PRINT_TIME
       if(method%commVec%myrank.eq.0)then
          write(*,*)"@ reset - CRS - set_CRS_Matrix_val_GPU @",etime-stime
       endif
#endif
       
!        return 
       stime=MPI_Wtime()
       call set_CUSPARCE(method)
       etime=MPI_Wtime()

#if PRINT_TIME
       if(method%commVec%myrank.eq.0)then
          write(*,*)"@ reset - CRS - CUSPARCE @",etime-stime
       endif
#endif

!        return 
       stime=MPI_Wtime()
       if(method%ipreflag.eq.KSP_PRECONDITION_NEUMANN)then
          if(ipreflag0.ne.KSP_PRECONDITION_NEUMANN)then
             call set_precon_Neumann_GPU(method)
          endif
          if(allocated(method%preconNeuman%Dinv))then
!             n = size(method%preconNeuman%Dinv)
             n = method%prm%n
             call cudaMemcpy_h2d_fp64(&
                  method%preconNeuman%Dinv,&
                  method%C_preconNeuman%Dinv,n,offset_zero,offset_zero)
          endif
          if(allocated(method%preconNeuman%D))then
!             n = size(method%preconNeuman%D)
             n = method%prm%n
             call cudaMemcpy_h2d_fp64(&
                  method%preconNeuman%D,&
                  method%C_preconNeuman%D,n,offset_zero,offset_zero)
          endif

          call set_CRS_Matrix_val_GPU( &
               method%preconNeuman%crsMat, &
               method%cusparce_preconcrsMatrix%crsA, &
               method%prm%n,method%prm%istart,method%prm%nnz &
               )
       endif
       etime=MPI_Wtime()
!       return 

#if PRINT_TIME
       if(method%commVec%myrank.eq.0)then
          write(*,*)"@ reset - CRS - Neumann @",etime-stime
       endif
#endif

    endif

    ! write(*,*)"@ reset_GPU : debug 3 @"
    call cudaDeviceSynchronize() 
    if(MatFormat.eq.MATRIX_FORMAT_DIA)then
       call cudaDeviceSynchronize() 
       if(allocated(method%diaMat%diaA))then
          n = method%diaMat%nnd * method%prm%n
          call cudaMemcpy_h2d_fp64(method%diaMat%diaA,method%C_diaMat%diaA,n,offset_zero,offset_zero)
       endif

       call cudaDeviceSynchronize() 
       if(allocated(method%Dinv))then
!          n = size(method%Dinv)
          n = method%prm%n
          call cudaMemcpy_h2d_fp64(&
               method%Dinv,&
               method%C_Dinv,n,offset_zero,offset_zero)
       endif

       ! write(*,*)"@ reset_GPU : debug 4 @"
       call set_CUSPARCE(method)
       ! write(*,*)"@ reset_GPU : debug 5 @"

       call cudaDeviceSynchronize() 
       if(method%ipreflag.eq.KSP_PRECONDITION_NEUMANN)then
          if(ipreflag0.ne.KSP_PRECONDITION_NEUMANN)then
             call set_precon_Neumann_GPU(method)
          endif
          if(allocated(method%preconNeuman%Dinv))then
!             n = size(method%preconNeuman%Dinv)
             n = method%prm%n
             call cudaMemcpy_h2d_fp64(&
                  method%preconNeuman%Dinv,&
                  method%C_preconNeuman%Dinv,n,offset_zero,offset_zero)
          endif
          if(allocated(method%preconNeuman%D))then
!             n = size(method%preconNeuman%D)
             n = method%prm%n
             call cudaMemcpy_h2d_fp64(&
                  method%preconNeuman%D,&
                  method%C_preconNeuman%D,n,offset_zero,offset_zero)
          endif
          if(allocated(method%preconNeuman%diaMat%diaA))then
!             n = size(method%preconNeuman%diaMat%diaA)
             n = method%diaMat%nnd * method%prm%n
             call cudaMemcpy_h2d_fp64(&
                  method%preconNeuman%diaMat%diaA,&
                  method%C_preconNeuman%diaMat%diaA,n,offset_zero,offset_zero)
          endif         

       endif
       ! write(*,*)"@ reset_GPU : debug 6 @"

    endif
#endif
  end subroutine reset_GPU

  subroutine free_GPU(method)
    use krylov
    implicit none
    type(KSP) :: method
    integer(kind=intType) i
    integer(kind=intType) nn
    ! initがcpu側のメモリが存在する場合にGPU側のメモリを確保するため解放する場合もCPU側で確認している.
    ! CPU側の解放の前にGPU側のメモリを解放する必要がある
    call free_CRS_Matrix_GPU(method%crsMat,method%C_crsMat)
    call free_CRS_Matrix_GPU(method%crsMat,method%C_precon_CRS_lo_precision)
    call free_DIA_Matrix_GPU(method%diaMat,method%C_diaMat)
    call free_VEC_GPU(method)
    call free_COMM_GPU(method)
    if(method%ipreflag.eq.KSP_PRECONDITION_NEUMANN)then
       call free_precon_Neumann_GPU(method)
    endif
    ! cusparce free
    !    if(method%ipreflag.eq.KSP_PRECONDITION_BLOCKJACOBI)then
    call cusparse_free(method%cusparce_solver)
    call cudaFree(method%cusparce_crsMatrix%crsA)
    call cudaFree(method%cusparce_crsMatrix%crsRow_ptr)
    call cudaFree(method%cusparce_crsMatrix%crsCol)
    call cudaFree(method%cusparce_preconcrsMatrix%crsA)
    call cudaFree(method%cusparce_preconcrsMatrix%crsRow_ptr)
    call cudaFree(method%cusparce_preconcrsMatrix%crsCol)

    if(method%ipreflag.eq.KSP_PRECONDITION_ADDITIVESCHWARZ)then
       call cusparseDestroyDnVec(method%C_ivec_AS)
       call cudaFree(method%C_ivec_AS)
       call cusparseDestroyDnVec(method%C_ovec_AS)
       call cudaFree(method%C_ovec_AS)

       call cudaFree(method%C_s_buffL)
       call cudaFree(method%C_s_buffU)
       call cudaFree(method%C_r_buffL)
       call cudaFree(method%C_r_buffL)

       call free(method%C_s_reqL) 
       call free(method%C_r_reqL) 
       call free(method%C_s_reqU) 
       call free(method%C_r_reqU) 
       call free(method%C_s_staL)
       call free(method%C_r_staL)
       call free(method%C_s_staU)
       call free(method%C_r_staU)
    endif

    call free_cublas(method%cublas_solver)

  end subroutine free_GPU

  subroutine free_COMM_GPU(method)
    use krylov
    implicit none
    type(KSP) :: method
    integer(kind=intType) i
    integer(kind=intType) nn
    if(allocated(method%commVec%istart))then
       call cudaFree(method%commVec%C_istart)
    endif
    if(allocated(method%commVec%nrk_prm))then
!       nn = size(method%commVec%nrk_prm)
       nn =  method%commVec%nrk_rank
       do i = 1, nn
          if(allocated(method%commVec%nrk_prm(i)%offset))then
             call cudaFree(method%commVec%nrk_prm(i)%C_offset)
          endif
          if(allocated(method%commVec%nrk_prm(i)%nrk_ptr))then
             call cudaFree(method%commVec%nrk_prm(i)%C_nrk_ptr)
          endif
          if(allocated(method%commVec%nrk_prm(i)%nrkRow))then
             call cudaFree(method%commVec%nrk_prm(i)%C_nrkRow)
          endif
          if(allocated(method%commVec%nrk_prm(i)%nrkCol))then
             call cudaFree(method%commVec%nrk_prm(i)%C_nrkCol)
          endif
       enddo
    end if

    if(allocated(method%commVec%nrk_prm_send))then
!       nn = size(method%commVec%nrk_prm_send)
       nn = method%commVec%nrk_rank_send
       do i = 1, nn
          if(allocated(method%commVec%nrk_prm_send(i)%offset))then
             call cudaFree(method%commVec%nrk_prm_send(i)%C_offset)
          endif
          if(allocated(method%commVec%nrk_prm_send(i)%nrk_ptr))then
             call cudaFree(method%commVec%nrk_prm_send(i)%C_nrk_ptr)
          endif
          if(allocated(method%commVec%nrk_prm_send(i)%nrkRow))then
             call cudaFree(method%commVec%nrk_prm_send(i)%C_nrkRow)
          endif
          if(allocated(method%commVec%nrk_prm_send(i)%nrkCol))then
             call cudaFree(method%commVec%nrk_prm_send(i)%C_nrkCol)
          endif
       enddo
    end if

    ! C MPI    
    call free(method%commVec%C_comm);
    if(method%commVec%nrk_rank_send>0)then
       call free(method%commVec%C_s_sta);
       call free(method%commVec%C_s_req);
    endif
    if(method%commVec%nrk_rank>0)then
       call free(method%commVec%C_r_sta);
       call free(method%commVec%C_r_req);
    endif

    ! if(allocated(method%@))then
    !    call cudaFree(method%C_@)
    ! endif

    ! --------
    ! MPI buffer
    nn = method%commVec%nrk_data_send
!    if(nn>0)then
       call cudaFree(method%C_sbuff)
!    endif
    if(allocated(method%buff))then
       call cudaFree(method%C_buff)
    endif
    if(allocated(method%buff_hi))then
       call cudaFree(method%C_buff_hi)
    endif
    if(allocated(method%buff_lo))then
       call cudaFree(method%C_buff_lo)
    endif

  end subroutine free_COMM_GPU

end module cu_krylov

module cu_comm
  use cu_common_wrapper
  implicit none

  interface

     subroutine set_sbuf_DDM_FP64_GPU( &
          numGPUblock, &
          numGPUthread, &
          vec,  offset_vec , &
          sbuf, offset_sbuf, &
          offset, &
          grids, &
          mx,mxy, &
          x0s,y0s,z0s, &
          xlen,ylen, &
          stmx,stmy,stmz &
          ) bind(c,name = "cuda_set_sbuf_DDM_FP64_GPU")
       import :: intType
       import :: devicemem
       implicit none

       integer,value :: numGPUblock
       integer,value :: numGPUthread
       type(devicemem) vec
       integer(kind=intType),value :: offset_vec
       type(devicemem) sbuf
       integer(kind=intType),value :: offset_sbuf
       integer,value :: offset
       integer,value :: grids
       integer,value :: mx
       integer,value :: mxy
       integer,value :: x0s
       integer,value :: y0s
       integer,value :: z0s
       integer,value :: xlen
       integer,value :: ylen
       integer,value :: stmx
       integer,value :: stmy
       integer,value :: stmz

     end subroutine set_sbuf_DDM_FP64_GPU

     subroutine set_rbuf_DDM_FP64_GPU( &
          numGPUblock, &
          numGPUthread, &
          vec, offset_vec, &
          rbuf, offset_rbuf, &
          offset, &
          grids, &
          mx,mxy, &
          x0r,y0r,z0r, &
          xlen,ylen, &
          stmx,stmy,stmz &
          ) bind(c,name = "cuda_set_rbuf_DDM_FP64_GPU")
       import :: intType
       import :: devicemem
       implicit none

       integer,value :: numGPUblock
       integer,value :: numGPUthread
       type(devicemem) vec
       integer(kind=intType),value ::  offset_vec
       type(devicemem) rbuf
       integer(kind=intType),value :: offset_rbuf
       integer,value :: offset
       integer,value :: grids
       integer,value :: mx
       integer,value :: mxy
       integer,value :: x0r
       integer,value :: y0r
       integer,value :: z0r
       integer,value :: xlen
       integer,value :: ylen
       integer,value :: stmx
       integer,value :: stmy
       integer,value :: stmz

     end subroutine set_rbuf_DDM_FP64_GPU

     subroutine sleev_zero_DDM_FP64_GPU( &
          numGPUblock, &
          numGPUthread, &
          vec, offset_vec, &
          offset, &
          grids, &
          mx,mxy, &
          x0r,y0r,z0r, &
          xlen,ylen, &
          stmx,stmy,stmz &
          ) bind(c,name = "cuda_sleev_zero_DDM_FP64_GPU")
       import :: intType
       import :: devicemem
       implicit none

       integer,value :: numGPUblock
       integer,value :: numGPUthread
       type(devicemem) vec
       integer(kind=intType),value :: offset_vec
       integer,value :: offset
       integer,value :: grids
       integer,value :: mx
       integer,value :: mxy
       integer,value :: x0r
       integer,value :: y0r
       integer,value :: z0r
       integer,value :: xlen
       integer,value :: ylen
       integer,value :: stmx
       integer,value :: stmy
       integer,value :: stmz

     end subroutine sleev_zero_DDM_FP64_GPU

#if 0
     subroutine comm_ddm( &
          n, &
          sbuf, rbuf, &
          grids, &
          myrank, &
          rank_s,rank_r, &
          comm_, &
          s_stat_,r_stat_, &
          s_req_,r_req_ &
          ) bind(c,name = "comm_ddm")
       import :: devicemem
       implicit none
       integer,value :: n
       real*8 sbuf(*)
       real*8 rbuf(*)
       integer grids(*)
       integer,value :: myrank
       integer rank_s(*)
       integer rank_r(*)
       type(devicemem) comm_
       type(devicemem) s_stat_
       type(devicemem) r_stat_
       type(devicemem) s_req_
       type(devicemem) r_req_
     end subroutine comm_ddm
#endif

     subroutine comm_ddm_GPU( &
          n, &
          sbuf, offset_sbuf,&
          rbuf, offset_rbuf,&
          grids, &
          myrank, &
          rank_s,rank_r, &
          comm_, &
          s_stat_,r_stat_, &
          s_req_,r_req_ &
          ) bind(c,name = "comm_ddm_FP64_GPU")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: n
       type(devicemem) sbuf
       integer(kind=intType),value :: offset_sbuf
       type(devicemem) rbuf
       integer(kind=intType),value :: offset_rbuf
       integer grids(*)
       integer,value :: myrank
       integer rank_s(*)
       integer rank_r(*)
       type(devicemem) comm_
       type(devicemem) s_stat_
       type(devicemem) r_stat_
       type(devicemem) s_req_
       type(devicemem) r_req_
     end subroutine comm_ddm_GPU

     subroutine comm_ddm_GPU_host( &
          n, &
          sbuf, offset_sbuf,&
          rbuf, offset_rbuf,&
          grids, &
          myrank, &
          rank_s,rank_r, &
          comm_, &
          s_stat_,r_stat_, &
          s_req_,r_req_ &
          ) bind(c,name = "comm_ddm_FP64_GPU_host")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: n
       type(devicemem) sbuf
       integer(kind=intType),value :: offset_sbuf
       type(devicemem) rbuf
       integer(kind=intType),value :: offset_rbuf
       integer grids(*)
       integer,value :: myrank
       integer rank_s(*)
       integer rank_r(*)
       type(devicemem) comm_
       type(devicemem) s_stat_
       type(devicemem) r_stat_
       type(devicemem) s_req_
       type(devicemem) r_req_
     end subroutine comm_ddm_GPU_host

     subroutine get_nrk_vecLU_add_FP64_GPU( &
               rbuff, offset_rbuff,&
               n, &
               nColL, &
               nColU, &
               vec, offset_vec &
               ) bind(c,name = "get_nrk_vecLU_add")
       import :: intType
       import :: devicemem
       implicit none
       type(devicemem) rbuff
       integer(kind=intType),value :: offset_rbuff
       integer(kind=intType),value :: n
       integer(kind=intType),value :: nColL
       integer(kind=intType),value :: nColU
       type(devicemem) vec
       integer(kind=intType),value :: offset_vec
       
     end subroutine get_nrk_vecLU_add_FP64_GPU
     
     subroutine get_nrk_vecLU_FP64_GPU( &
          numGPUblock,numGPUthread, &
          rbuff,  offset_rbuff, &
          nColL,nColU, &
          rbuffL, offset_rbuffL, &
          rbuffU, offset_rbuffU  &
          ) bind(c,name = "cuda_get_nrk_vecLU")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       type(devicemem) rbuff
       integer(kind=intType),value :: offset_rbuff
       type(devicemem) rbuffL
       integer(kind=intType),value :: offset_rbuffL
       type(devicemem) rbuffU
       integer(kind=intType),value :: offset_rbuffU
       integer(kind=intType),value :: nColL,nColU

     end subroutine get_nrk_vecLU_FP64_GPU

     subroutine nrk_vec_pack_FP64_GPU( &
          numGPUblock,numGPUthread, &
          vec,  offset_vec,&
          sbuf, offset_sbuf,&
          offset,&
          ndata,i0) &
          bind(c,name = "cuda_nrk_vec_pack")
       import :: intType
       import :: devicemem
       implicit none

       integer,value :: numGPUblock
       integer,value :: numGPUthread
       type(devicemem) vec
       type(devicemem) sbuf
       integer(kind=intType),value :: offset_vec
       type(devicemem) offset 
       integer(kind=intType),value :: offset_sbuf
       integer(kind=intType),value :: ndata
       integer(kind=intType),value :: i0
       
     end subroutine nrk_vec_pack_FP64_GPU
       
       ! subroutine nrk_vec_comm(sbuf,rbuf,) &
       !      bind(c,name = "cuda_nrk_vec_comm")
       !   import :: intType
       !   import :: devicemem
       !   implicit none
       ! end subroutine nrk_vec_comm
     
     subroutine nrk_vec_comm_FP64_C( &
          f_comm, &
          sbuf, &
          rbuf, &
          nrk_rank_send, &
          nrk_prm_send_offset, &
          nrk_prm_send_ndata, &
          nrk_prm_send_rank, &
          nrk_rank_recv, &
          nrk_prm_recv_offset, &
          nrk_prm_recv_ndata, &
          nrk_prm_recv_rank, &
          myrank &
          ) bind(c,name="cuda_nrk_vec_comm_FP64")
       import :: devicemem
       implicit none
       integer,value :: f_comm
       type(devicemem) sbuf
       type(devicemem) rbuf
       integer,value ::  nrk_rank_send
       integer nrk_prm_send_offset(*)
       integer nrk_prm_send_ndata(*)
       integer nrk_prm_send_rank(*)
       integer,value ::  nrk_rank_recv
       integer nrk_prm_recv_offset(*)
       integer nrk_prm_recv_ndata(*)
       integer nrk_prm_recv_rank(*)
       integer,value :: myrank

     end subroutine nrk_vec_comm_FP64_C

    end interface

contains

  subroutine nrk_vec_comm_FP64_GPU(commVec,sbuf_,rbuf_)
    use mpi
    use comm
    use PARCEL_integer
    implicit none
    type(comm_Vec) :: commVec
    type(devicemem) sbuf_
    type(devicemem) rbuf_

#if CUDA_AWARE

#else
    real*8 sbuf(commVec%nrk_data_send)
    real*8 rbuf(commVec%nrk_data)

    integer s_req(commVec%nrk_rank_send)
    integer r_req(commVec%nrk_rank)

    type send_stat
       integer stat(MPI_STATUS_SIZE)
    end type send_stat
    type recv_stat
       integer stat(MPI_STATUS_SIZE)
    end type recv_stat
    type(send_stat) :: s_sta(commVec%nrk_rank_send)
    type(recv_stat) :: r_sta(commVec%nrk_rank)
#endif

    integer(kind=intType) i,j,k
    integer(kind=intType) i0
    integer(kind=intType) parcel_int
    integer ierr
    integer int1,int2

    integer nrk_prm_send_offset( commVec%nrk_rank_send )
    integer nrk_prm_send_ndata( commVec%nrk_rank_send )
    integer nrk_prm_send_rank( commVec%nrk_rank_send )

    integer nrk_prm_recv_offset( commVec%nrk_rank )
    integer nrk_prm_recv_ndata( commVec%nrk_rank )
    integer nrk_prm_recv_rank( commVec%nrk_rank )
    integer(kind=intType),parameter :: offset_zero=0

    call cudaDeviceSynchronize() 

#if CUDA_AWARE
    do i=1,commVec%nrk_rank_send
       i0=0
       do k=1,i-1
          i0=i0+commVec%nrk_prm_send(k)%ndata
       enddo
       nrk_prm_send_offset(i) = i0
       nrk_prm_send_ndata(i)  = commVec%nrk_prm_send(i)%ndata
       nrk_prm_send_rank(i)   = commVec%nrk_prm_send(i)%rank
    enddo
    do i=commVec%nrk_rank,1,-1
       i0=0
       do k=1,i-1
          i0=i0+commVec%nrk_prm(k)%ndata
       enddo
       nrk_prm_recv_offset(i) = i0
       nrk_prm_recv_ndata(i)  = commVec%nrk_prm(i)%ndata
       nrk_prm_recv_rank(i)   = commVec%nrk_prm(i)%rank
    enddo

    call nrk_vec_comm_FP64_C( &
          commVec%comm, &
          sbuf_, &
          rbuf_, &
          commVec%nrk_rank_send, &
          nrk_prm_send_offset, &
          nrk_prm_send_ndata, &
          nrk_prm_send_rank, &
          commVec%nrk_rank, &
          nrk_prm_recv_offset, &
          nrk_prm_recv_ndata, &
          nrk_prm_recv_rank, &
          commVec%myrank &
          )
    return
#else
    sbuf = 0.0d0
    rbuf = 0.0d0
    parcel_int = commVec%nrk_data_send
    call cudaMemcpy_d2h_fp64(sbuf_,sbuf,parcel_int,offset_zero,offset_zero)
    call cudaDeviceSynchronize() 
    do i=1,commVec%nrk_rank_send
       i0=0
       do k=1,i-1
          i0=i0+commVec%nrk_prm_send(k)%ndata
       enddo

#ifdef USE_C_MPI
       call C_MPI_Isend_FP64buf_FP64( &
            sbuf, &
            commVec%nrk_prm_send(i)%ndata ,&
            commVec%nrk_prm_send(i)%rank,&
            commVec%myrank,&
            commVec%C_comm, &
            commVec%C_s_req, &
            i,(1+i0) &
            )
#else
       call MPI_Isend( &
            sbuf(1+i0),&
            commVec%nrk_prm_send(i)%ndata ,&
            MPI_REAL8,&
            commVec%nrk_prm_send(i)%rank,&
            commVec%myrank,&
            commVec%comm,&
            s_req(i),&
            ierr)
#endif

    enddo
    !    do i=1,commVec%nrk_rank
    do i=commVec%nrk_rank,1,-1
       i0=0
       do k=1,i-1
          i0=i0+commVec%nrk_prm(k)%ndata
       enddo

#ifdef USE_C_MPI
       call C_MPI_Irecv_FP64buf_FP64( &
            rbuf,&
            commVec%nrk_prm(i)%ndata ,&
            commVec%nrk_prm(i)%rank,&
            commVec%nrk_prm(i)%rank,&
            commVec%C_comm, &
            commVec%C_r_req, &
            i,(1+i0) &
            )
#else

       call MPI_Irecv( &
            rbuf(1+i0),&
            commVec%nrk_prm(i)%ndata ,&
            MPI_REAL8, &
            commVec%nrk_prm(i)%rank, &
            commVec%nrk_prm(i)%rank, &
            commVec%comm, &
            r_req(i), &
            ierr)
#endif

    enddo

    do i=1,commVec%nrk_rank_send
#ifdef USE_C_MPI
       call C_MPI_Wait(commVec%C_s_req,commVec%C_s_sta,i)
#else
       call MPI_Wait(s_req(i),s_sta(i)%stat,ierr)
#endif
    enddo

    do i=commVec%nrk_rank,1,-1
#ifdef USE_C_MPI
       call C_MPI_Wait( commVec%C_r_req,commVec%C_r_sta,i)
#else
       call MPI_Wait(r_req(i),r_sta(i)%stat,ierr)
#endif
    enddo

    parcel_int = commVec%nrk_data
    call cudaMemcpy_h2d_fp64(rbuf,rbuf_,parcel_int,offset_zero,offset_zero) 
    call cudaDeviceSynchronize() 
#endif

  end subroutine nrk_vec_comm_FP64_GPU

  subroutine get_nrk_vec_FP64_GPU( &
       commVec, &
       vec,  offset_vec,  &
       sbuf, offset_sbuf, &
       rbuf, offset_rbuf)

    use comm
    use PARCEL_integer
    implicit none
    type(comm_Vec) :: commVec
!    real*8 vec(commVec%n)
!    real*8 sbuf(commVec%nrk_data_send)
!    real*8 rbuf(commVec%nrk_data)
    type(devicemem) vec
    integer(kind=intType) offset_vec
    type(devicemem) sbuf
    integer(kind=intType) offset_sbuf
    type(devicemem) rbuf
    integer(kind=intType) offset_rbuf
    integer(kind=intType) j,k,i0
    integer(kind=intType) parcel_int

    ! 送信用データをバッファに詰める ( nrk_vec_pack )
    do j=1,commVec%nrk_rank_send
       i0=0
       do k=1,j-1
          i0=i0+commVec%nrk_prm_send(k)%ndata
       enddo

       parcel_int = commVec%nrk_prm_send(j)%ndata
       call nrk_vec_pack_FP64_GPU( &
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            vec, offset_vec, &
            sbuf, offset_sbuf, &
            commVec%nrk_prm_send(j)%C_offset, &
            parcel_int, &
            i0)
    enddo
    call nrk_vec_comm_FP64_GPU(commVec,sbuf,rbuf)

    ! nrk_vec_comm

  end subroutine get_nrk_vec_FP64_GPU

end module cu_comm
