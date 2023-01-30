#define PRINT_TIME 0

module krylov_kernel_gpu
  use PARCEL_define
  use omp_lib
  use cu_krylov
  implicit none
contains

  subroutine set_KSP_DIA_cuda(&
       icomm, &
       vecx,vecb, &
       n,gn,istart, &
       diaA,offset,nnd,&
       method)
    use krylov
    use matrix_kernel
    use comm_kernel,only : set_comm_dia_GPUmem
    use cu_Vec,only : &
         cuda_get_diamatrix_neumanndiag_FP64_GPU, &
         cuda_get_diamatrix_diag_FP64_GPU
    use iso_c_binding

    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) istart ! ベクトルの開始位置
    integer(kind=intType) iend ! ベクトルの終了位置
    type(c_ptr) vecx ! 解ベクトル
    type(c_ptr) vecb ! 右辺ベクトル
    integer(kind=intType) nnd
    type(c_ptr) diaA
    type(c_ptr) offset
    type(KSP) :: method

    real*8,pointer :: vecx_Fptr(:)
    real*8,pointer :: vecb_Fptr(:)
    real*8,pointer :: diaA_Fptr(:)
    integer(kind=intType),pointer :: offset_Fptr(:)
    integer(kind=intType),parameter :: offset_zero=0

    integer(kind=intType) parcel_int
    integer(kind=intType) nColL,nColU
    integer(kind=intType) i
    integer(kind=intType) offset_diag,idiag

    integer ierr
    real*8 stime,etime
    real*8 stime0,etime0
    integer myrank
    integer(kind=intType) offset_neuman(nnd)

    stime0 = MPI_Wtime()
    stime  = MPI_Wtime()

    call set_cublas(method%cublas_solver)

    ! ------
    ! c_ptrをFortranポインタに変換
    call c_f_pointer(vecx,vecx_Fptr,[n])
    call c_f_pointer(vecb,vecb_Fptr,[n])
    call c_f_pointer(diaA,diaA_Fptr,[nnd*n])
    call c_f_pointer(offset,offset_Fptr,[nnd])

    ! ------
    ! DIA Format
    iend = istart + n -1
    method%prm%MatFormat = MATRIX_FORMAT_DIA 

    method%commVec%comm = icomm
    method%prm%n        = n
    method%prm%gn       = gn
    method%prm%istart   = istart
    method%prm%iend     = iend
    method%prm%nnd      = nnd
    call MPI_Comm_rank(icomm,myrank,ierr)
    method%commVec%myrank = myrank

    method%C_diaMat%nnd = nnd
    if( allocated(method%reshistory) )then
       deallocate(method%reshistory)
    endif
    allocate(method%reshistory(max(method%itrmax+100,100)))
    !$OMP PARALLEL DO PRIVATE(i)
    do i=1,max(method%itrmax+100,100)
       method%reshistory(i)=0.0d0
    enddo

    etime = MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"   @ 1-1 @",etime - stime
    endif
#endif
    stime  = MPI_Wtime()

    ! Vector
    call cudaMalloc_fp64(method%C_r,n)
    call cudaMalloc_fp64(method%C_x,n)
    call cudaMalloc_fp64(method%C_b,n)

    ! DIA Matrix
    parcel_int = nnd*n
    call cudaMalloc_fp64(method%C_diaMat%diaA        ,parcel_int)
    call cudaMalloc_PARCEL_int(method%C_diaMat%offset,nnd)

    etime = MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"      @ 1-2-1 @",etime - stime
    endif
#endif
    stime  = MPI_Wtime()

    ! 入力ベクトルをコピー
    call cudaMemcpy_dptr2d_double(vecx_Fptr,method%C_x,n)
    call cudaMemcpy_dptr2d_double(vecb_Fptr,method%C_b,n)

    ! DIA Matrixをコピー
    parcel_int = nnd*n
    call cudaMemcpy_dptr2d_double(diaA_Fptr      ,method%C_diaMat%diaA  ,parcel_int)
    call cudaMemcpy_dptr2d_parcel_int(offset_Fptr,method%C_diaMat%offset,nnd)

    etime = MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"      @ 1-2-2 @",etime - stime
    endif
#endif
    stime  = MPI_Wtime()

    if( allocated(method%diaMat%diaA) )then
       deallocate(method%diaMat%diaA)
    endif
    parcel_int = nnd * n
    allocate(method%diaMat%diaA(parcel_int))

    if( allocated(method%diaMat%iLs) )then
       deallocate(method%diaMat%iLs)
    endif
    parcel_int = nnd * n
    allocate(method%diaMat%iLs(parcel_int))

    ! !$OMP PARALLEL DO PRIVATE(i)
    ! do i = 1, parcel_int
    !    method%diaMat%iLs(i) = -1
    ! enddo

    if( allocated(method%diaMat%offset) )then
       deallocate(method%diaMat%offset)
    endif
    allocate(method%diaMat%offset(nnd))

    etime = MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"      @ 1-2-3 @",etime - stime
    endif
#endif
    stime  = MPI_Wtime()

#if 0
    parcel_int = nnd * n
    call cudaMemcpy_d2h_fp64(method%C_diaMat%diaA        , method%diaMat%diaA  ,parcel_int,offset_zero,offset_zero)
#endif
    call cudaMemcpy_d2h_PARCEL_INT(method%C_diaMat%offset, method%diaMat%offset,nnd       ,offset_zero,offset_zero)

    etime = MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"      @ 1-2-4 @",etime - stime
    endif
#endif
    stime  = MPI_Wtime()

    ! setting communication
    ! 各プロセスのistart<-->iendを超える場合は除外される
    nColL = -1
    nColU = -1
    call set_comm_dia_GPUmem( &
       method%prm, &
       method%commVec, &
       istart, &
       n, nnd, &
       nColL, nColU, &
       method%diaMat%diaA, method%diaMat%offset&
       )

    etime = MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"   @ 1-3 @",etime - stime
    endif
#endif
    stime  = MPI_Wtime()

    ! iLs計算
    call init_MatVec_dia( &
         method%prm, &
         method%diaMat,  &
         method%commVec, &
         method%nBlock)

    etime = MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"   @ 1-4 @",etime - stime
    endif
#endif
    stime  = MPI_Wtime()

    parcel_int = nnd * n
    call cudaMalloc_PARCEL_INT(  method%C_diaMat%iLs,parcel_int)
    call cudaMemcpy_h2d_PARCEL_INT(method%diaMat%iLs,method%C_diaMat%iLs,parcel_int,offset_zero,offset_zero)

    call set_KSP_MPI_malloc_cuda(method)

    etime = MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"   @ 1-5 @",etime - stime
    endif
#endif
    stime  = MPI_Wtime()

    ! precon
    method%prm_AS%MatFormat = MATRIX_FORMAT_DIA
    if(method%ipreflag.eq.KSP_PRECONDITION_NEUMANN_LOCAL)then
       method%ipreflag = KSP_PRECONDITION_NEUMANN
       method%C_preconNeuman%local_flag = -1
    endif
    method%cusparce_solver%ipreflag = method%ipreflag

    select case (method%ipreflag)
    case(KSP_PRECONDITION_NONE)

    case(KSP_PRECONDITION_POINTJACOBI)
       ! cusparseを利用しないのでベクトル記述子は不要
       ! write(*,*)" n = ",n
       call cudaMalloc_fp64(method%C_D,n)
       call cudaMalloc_fp64(method%C_Dinv,n)
       idiag = 1
       do i = 1,nnd
          if(method%diaMat%offset(i).eq.0)then
             idiag = i
             exit
          endif
       enddo
       call cuda_get_diamatrix_diag_FP64_GPU( &
          n, &
          nnd, &
          idiag, &
          method%C_diaMat%diaA, &
          method%C_D, &
          offset_zero, &
          method%C_Dinv, &
          offset_zero &
          )

    case(KSP_PRECONDITION_BLOCKJACOBI)
    case(KSP_PRECONDITION_ADDITIVESCHWARZ)

    case(KSP_PRECONDITION_NEUMANN)
       ! 対角位置さえわかればGPU側で計算可能 + 対角ゼロにしているけど,元対角位置のデータに値を詰めてnndを-1したほうが効率的. 問題は動作するかどうか?

       parcel_int = nnd * n
       call cudaMalloc_fp64(method%C_preconNeuman%diaMat%diaA  ,parcel_int)
       call cudaMalloc_PARCEL_int(method%C_preconNeuman%diaMat%offset,nnd)
       parcel_int = nnd * n
       call  cudaMalloc_PARCEL_int(method%C_preconNeuman%diaMat%iLs   ,parcel_int)

       parcel_int = nnd * n
       call cudaMemcpy_d2d_FP64( &
            method%C_diaMat%diaA,method%C_preconNeuman%diaMat%diaA,parcel_int,offset_zero,offset_zero)
       call cudaMemcpy_d2d_PARCEL_INT(&
            method%C_diaMat%offset,method%C_preconNeuman%diaMat%offset,nnd,offset_zero,offset_zero)

       parcel_int = nnd * n
       call cudaMemcpy_d2d_PARCEL_INT(&
            method%C_diaMat%iLs,method%C_preconNeuman%diaMat%iLs,parcel_int,offset_zero,offset_zero)

       call cudaMalloc_fp64(method%C_preconNeuman%vecx,n)
       call cudaMalloc_fp64(method%C_preconNeuman%vecy,n)
       call cudaMalloc_fp64(method%C_preconNeuman%D,n)
       call cudaMalloc_fp64(method%C_preconNeuman%Dinv,n)

       idiag = 1
       do i = 1,nnd
          if(method%diaMat%offset(i).eq.0)then
             idiag = i
             exit
          endif
       enddo

       call cuda_get_diamatrix_neumanndiag_FP64_GPU( &
            n, &
          nnd, &
          idiag, &
          method%C_preconNeuman%diaMat%diaA, &
          method%C_preconNeuman%diaMat%offset, &
          method%C_preconNeuman%diaMat%iLs, &
          method%C_preconNeuman%D, &
          offset_zero, & 
          method%C_preconNeuman%Dinv, &
          offset_zero & 
          )
#if 1
       do i = 1,nnd
          offset_neuman(i) = method%diaMat%offset(i)
       enddo
       offset_neuman(idiag) = method%diaMat%offset(nnd)
       offset_neuman(nnd)   = 0 
       method%C_preconNeuman%diaMat%nnd = nnd - 1
       call cudaMemcpy_h2d_PARCEL_INT(offset_neuman,method%C_preconNeuman%diaMat%offset,nnd,offset_zero,offset_zero)
#else
       method%C_preconNeuman%diaMat%nnd = nnd 
#endif

    end select

    deallocate(method%diaMat%diaA)
    deallocate(method%diaMat%iLs)
    deallocate(method%diaMat%offset)

    etime = MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"   @ 1-6 @",etime - stime
    endif
#endif
    stime  = MPI_Wtime()

    etime0 = MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"@ 1 @",etime0 - stime0
    endif
#endif

  end subroutine set_KSP_DIA_cuda

  ! method%C_crsMat%crsCol,method%C_crsMat%crsRow_ptr : メモリ領域はPARCEL_INTだがcusparceに利用するのでデータ型はint
  subroutine set_KSP_CRS_cuda(&
       icomm, &
       vecx,vecb, &
       n,gn,nnz,istart, &
       crsA,crsRow_ptr,crsCol, &
       method)

    use krylov
    use cu_krylov
    use matrix_kernel
    use common_kernel
    use cu_Vec,only : &
         cuda_get_crsmatrix_neumanndiag_FP64_GPU, &
         cuda_get_crsmatrix_diag_FP64_GPU
    use comm_kernel,only : set_comm_crs_GPUmem
    use iso_c_binding

    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    integer(kind=intType) iend ! ベクトルの終了位置
    type(c_ptr) vecx ! 解ベクトル
    type(c_ptr) vecb ! 右辺ベクトル
    type(c_ptr) crsA ! 行列値
    type(c_ptr) crsRow_ptr ! ポインタ 
    type(c_ptr) crsCol ! 列番号

    real*8,pointer :: vecx_Fptr(:)
    real*8,pointer :: vecb_Fptr(:)
    real*8,pointer :: crsA_Fptr(:)
    integer(kind=intType),pointer :: crsRow_ptr_Fptr(:)
    integer(kind=intType),pointer :: crsCol_Fptr(:)

    integer itrmax ! 最大反復回数
    type(KSP) :: method

    integer ierr

    ! ------
    real*8  crsA_cpu(nnz) ! 行列値
    integer(kind=intType) crsRow_ptr_cpu(n+1) ! ポインタ 
    integer(kind=intType) crsCol_cpu(nnz)     ! 列番号

    integer crsRow_ptr_cpu_int4(n+1) ! ポインタ 
    integer crsCol_cpu_int4(nnz)     ! 列番号

    integer(kind=intType) nCol
    integer(kind=intType) nColL
    integer(kind=intType) nColU
    integer(kind=intType) i,j
    integer(kind=intType) parcel_int
    integer(kind=intType),parameter :: offset_zero=0

    integer(kind=intType) n_precon 
    integer(kind=intType) nCol_precon 
    integer(kind=intType) nnz_precon

    real*8 stime,etime
    real*8 stime0,etime0
    integer myrank

    stime0 = MPI_Wtime()
    stime  = MPI_Wtime()

    call set_cublas(method%cublas_solver)

    ! ------
    ! c_ptrをFortranポインタに変換
    call c_f_pointer(vecx,vecx_Fptr,[n])
    call c_f_pointer(vecb,vecb_Fptr,[n])
    call c_f_pointer(crsA,crsA_Fptr,[nnz])
    call c_f_pointer(crsRow_ptr,crsRow_ptr_Fptr,[n+1])
    call c_f_pointer(crsCol,crsCol_Fptr,[nnz])

    ! ------
    ! CRS Format
    iend  = istart + n -1
    method%prm%MatFormat  = MATRIX_FORMAT_CRS 

    method%commVec%comm = icomm
    method%prm%n        = n
    method%prm%gn       = gn
    method%prm%nnz      = nnz
    method%prm%istart   = istart
    method%prm%iend     = iend
    call MPI_Comm_rank(icomm,myrank,ierr)
    method%commVec%myrank = myrank

    if( allocated(method%reshistory) )then
       deallocate(method%reshistory)
    endif
    allocate(method%reshistory(max(method%itrmax+100,100)))
    !$OMP PARALLEL DO PRIVATE(i)
    do i=1,max(method%itrmax+100,100)
       method%reshistory(i)=0.0d0
    enddo

    call MPI_Barrier(icomm,ierr)

    etime = MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"   @ 1-1 @",etime - stime
    endif
#endif
    stime  = MPI_Wtime()

    ! Vector
    call cudaMalloc_fp64(method%C_r,n)
    call cudaMalloc_fp64(method%C_x,n)
    call cudaMalloc_fp64(method%C_b,n)

    ! 入力ベクトルをコピー
    call cudaDeviceSynchronize() 
    call cudaMemcpy_dptr2d_double(vecx_Fptr,method%C_x,n)
    call cudaMemcpy_dptr2d_double(vecb_Fptr,method%C_b,n)
    call cudaDeviceSynchronize() 

    ! CRS Matrix
    call cudaMalloc_fp64(method%C_crsMat%crsA,nnz)
    parcel_int = n + 1
    call cudaMalloc_PARCEL_int(method%C_crsMat%crsRow_ptr,parcel_int)
    call cudaMalloc_PARCEL_int(method%C_crsMat%crsCol,nnz)

    ! CRS Matrixをコピー
    call cudaDeviceSynchronize() 
    call cudaMemcpy_dptr2d_double(&
         crsA_Fptr, method%C_crsMat%crsA, &
         nnz &
         )
    parcel_int = n + 1
    call cudaMemcpy_dptr2d_parcel_int(&
         crsRow_ptr_Fptr, method%C_crsMat%crsRow_ptr, &
         parcel_int &
         )
    call cudaMemcpy_dptr2d_parcel_int(&
         crsCol_Fptr, method%C_crsMat%crsCol, &
         nnz &
         )
    call cudaDeviceSynchronize() 

    etime = MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"   @ 1-2 @",etime - stime
    endif
#endif
    stime  = MPI_Wtime()

    ! CPU側にコピー 
    call cudaMemcpy_d2h_fp64      (method%C_crsMat%crsA      ,crsA_cpu      ,nnz       ,offset_zero,offset_zero)
    parcel_int = n + 1
    call cudaMemcpy_d2h_PARCEL_INT(method%C_crsMat%crsRow_ptr,crsRow_ptr_cpu,parcel_int,offset_zero,offset_zero)
    call cudaMemcpy_d2h_PARCEL_INT(method%C_crsMat%crsCol    ,crsCol_cpu    ,nnz       ,offset_zero,offset_zero)
    call cudaDeviceSynchronize() 

#if 0
    if(myrank.eq.0)then
       write(*,*)" Original CRS"
!       do i = 1,5
       do i = n-5,n
          write(*,*)i," : ",crsRow_ptr_cpu(i)," --> ",crsRow_ptr_cpu(i+1)-1
          do j = crsRow_ptr_cpu(i),crsRow_ptr_cpu(i+1)-1
             write(*,*)"     ",i,crsCol_cpu(j),crsA_cpu(j)
          enddo
       enddo
    endif
#endif

    etime = MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"   @ 1-3 @",etime - stime
    endif
#endif
    stime  = MPI_Wtime()

    ! setting communication
    nColL = -1
    nColU = -1
    call set_comm_crs_GPUmem( &
       method%prm, &
       method%commVec, &
       istart, &
       n, nnz, &
       nColL, nColU, &
       crsA_cpu,crsCol_cpu,crsRow_ptr_cpu &
       )   

    method%crsMat%nnz = nnz
    method%crsMat%nColL = nColL
    method%crsMat%nColU = nColU

    method%C_crsMat%nnz   = nnz
    method%C_crsMat%nColL = nColL
    method%C_crsMat%nColU = nColU

    etime = MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"   @ 1-4 @",etime - stime
    endif
#endif
    stime  = MPI_Wtime()

#ifdef PARCEL_INT8
    !$OMP PARALLEL DO PRIVATE(i) 
    do i = 1,nnz
       crsCol_cpu_int4(i) = crsCol_cpu(i) 
    enddo
    !$OMP PARALLEL DO PRIVATE(i) 
    do i = 1,n+1
       crsRow_ptr_cpu_int4(i) = crsRow_ptr_cpu(i) 
    enddo

    ! GPU側にcusparse行列をコピー
    call cudaDeviceSynchronize() 
    call cudaMemcpy_h2d_fp64(crsA_cpu      ,method%C_crsMat%crsA      ,nnz,offset_zero,offset_zero)
    call cudaMemcpy_h2d_INT(crsRow_ptr_cpu_int4,method%C_crsMat%crsRow_ptr,n+1,offset_zero,offset_zero)
    call cudaMemcpy_h2d_INT(crsCol_cpu_int4    ,method%C_crsMat%crsCol    ,nnz,offset_zero,offset_zero)
    call cudaDeviceSynchronize() 

#else
    ! GPU側にcusparse行列をコピー
    call cudaDeviceSynchronize() 
    call cudaMemcpy_h2d_fp64(crsA_cpu      ,method%C_crsMat%crsA      ,nnz,offset_zero,offset_zero)
    call cudaMemcpy_h2d_INT(crsRow_ptr_cpu,method%C_crsMat%crsRow_ptr,n+1,offset_zero,offset_zero)
    call cudaMemcpy_h2d_INT(crsCol_cpu    ,method%C_crsMat%crsCol    ,nnz,offset_zero,offset_zero)
    call cudaDeviceSynchronize() 
#endif

    etime = MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"   @ 1-5 @",etime - stime
    endif
#endif
    stime  = MPI_Wtime()

    call set_KSP_MPI_malloc_cuda(method)

    etime = MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"   @ 1-6 @",etime - stime
    endif
#endif
    stime  = MPI_Wtime()

    ! precon
    method%prm_AS%MatFormat = MATRIX_FORMAT_CRS
    if(method%ipreflag.eq.KSP_PRECONDITION_NEUMANN_LOCAL)then
       method%ipreflag = KSP_PRECONDITION_NEUMANN
       method%C_preconNeuman%local_flag = -1
    endif
    ! write(*,*)" local_flag = ",method%C_preconNeuman%local_flag
    method%cusparce_solver%ipreflag = method%ipreflag
    nCol = n + nColL + nColU

    select case (method%ipreflag)
    case(KSP_PRECONDITION_NONE)

    case(KSP_PRECONDITION_POINTJACOBI)

       call cudaMalloc_fp64(method%C_D,nCol)
       call cusparseCreateDnVec_fp64(method%C_D,method%prm%n,nCol,offset_zero)
       call cudaMalloc_fp64(method%C_Dinv,nCol)
       call cusparseCreateDnVec_fp64(method%C_Dinv,method%prm%n,nCol,offset_zero)
       call cuda_get_crsmatrix_diag_FP64_GPU( &
            n, &
            method%C_crsMat%crsA, &
            method%C_crsMat%crsRow_ptr, &
            method%C_crsMat%crsCol, &
            method%C_D   , offset_zero, &
            method%C_Dinv, offset_zero  &          
            )
       method%cusparce_preconcrsMatrix%n    = n
       method%cusparce_preconcrsMatrix%nCol = nCol
       method%cusparce_preconcrsMatrix%nnz  = nnz

    case(KSP_PRECONDITION_BLOCKJACOBI)

       n_precon = n
       nnz_precon = nnz

       !$OMP PARALLEL DO PRIVATE(i) 
       do i = 1,nnz
          crsCol_cpu_int4(i) = crsCol_cpu(i) + 1
       enddo
       !$OMP PARALLEL DO PRIVATE(i) 
       do i = 1,n+1
          crsRow_ptr_cpu_int4(i) = crsRow_ptr_cpu(i) + 1
       enddo 

       call make_DiagBlock_Localcrs( &
            n_precon,nnz_precon,&
            crsA_cpu, &
            crsRow_ptr_cpu_int4, &
            crsCol_cpu_int4 &
            )
       !$OMP PARALLEL DO PRIVATE(i) 
       do i = 1,nnz
          crsCol_cpu_int4(i) = crsCol_cpu_int4(i) - 1
       enddo
       !$OMP PARALLEL DO PRIVATE(i) 
       do i = 1,n+1
          crsRow_ptr_cpu_int4(i) = crsRow_ptr_cpu_int4(i) - 1
       enddo 
    
       method%cusparce_preconcrsMatrix%n    = n
       method%cusparce_preconcrsMatrix%nCol = n
       method%cusparce_preconcrsMatrix%nnz  = nnz_precon

       call cudaMalloc_fp64(method%cusparce_preconcrsMatrix%crsA,nnz_precon)
       call cudaMalloc_PARCEL_int(method%cusparce_preconcrsMatrix%crsCol,nnz_precon)
       call cudaMalloc_PARCEL_int(method%cusparce_preconcrsMatrix%crsRow_ptr,n+1)

       ! call cudaDeviceSynchronize() 
       call cudaMemcpy_d2d_FP64(&
            method%C_crsMat%crsA,method%cusparce_preconcrsMatrix%crsA,nnz_precon,offset_zero,offset_zero)
       call cudaMemcpy_h2d_INT(&
            crsCol_cpu_int4    ,method%cusparce_preconcrsMatrix%crsCol,nnz_precon,offset_zero,offset_zero)       
       call cudaMemcpy_h2d_INT(&
            crsRow_ptr_cpu_int4,method%cusparce_preconcrsMatrix%crsRow_ptr,n+1,offset_zero,offset_zero)

   case(KSP_PRECONDITION_ADDITIVESCHWARZ)

   case(KSP_PRECONDITION_NEUMANN)
       call cudaMalloc_fp64(method%cusparce_preconcrsMatrix%crsA,nnz)
       call cudaMalloc_int(method%cusparce_preconcrsMatrix%crsCol,nnz)
       call cudaMalloc_int(method%cusparce_preconcrsMatrix%crsRow_ptr,n+1)

       call cudaMemcpy_d2d_FP64(&
            method%C_crsMat%crsA,method%cusparce_preconcrsMatrix%crsA,nnz,offset_zero,offset_zero)
       call cudaMemcpy_d2d_INT(&
            method%C_crsMat%crsCol,method%cusparce_preconcrsMatrix%crsCol,nnz,offset_zero,offset_zero)       
       call cudaMemcpy_d2d_INT(&
            method%C_crsMat%crsRow_ptr,method%cusparce_preconcrsMatrix%crsRow_ptr,n+1,offset_zero,offset_zero)
       method%cusparce_preconcrsMatrix%n    = n
       method%cusparce_preconcrsMatrix%nCol = nCol
       method%cusparce_preconcrsMatrix%nnz  = nnz

       call cudaMalloc_fp64(method%C_preconNeuman%vecx,nCol)
       call cusparseCreateDnVec_fp64(method%C_preconNeuman%vecx,method%prm%n,nCol,offset_zero)

       call cudaMalloc_fp64(method%C_preconNeuman%vecy,nCol)
       call cusparseCreateDnVec_fp64(method%C_preconNeuman%vecy,method%prm%n,nCol,offset_zero)
       
       call cudaMalloc_fp64(method%C_preconNeuman%D,nCol)
       call cusparseCreateDnVec_fp64(method%C_preconNeuman%D,method%prm%n,nCol,offset_zero)

       call cudaMalloc_fp64(method%C_preconNeuman%Dinv,nCol)
       call cusparseCreateDnVec_fp64(method%C_preconNeuman%Dinv,method%prm%n,nCol,offset_zero)
       
       call cuda_get_crsmatrix_neumanndiag_FP64_GPU( &
          n, &
          method%cusparce_preconcrsMatrix%crsA, &
          method%cusparce_preconcrsMatrix%crsRow_ptr, &
          method%cusparce_preconcrsMatrix%crsCol, &
          method%C_preconNeuman%D, offset_zero, &
          method%C_preconNeuman%Dinv, offset_zero &
          )
       n_precon    = n
       nnz_precon  = nnz
       nCol_precon = nCol
    end select

    etime = MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"   @ 1-7 @",etime - stime
    endif
#endif
    stime  = MPI_Wtime()

    ! cusparse

    call cudaDeviceSynchronize() 
#if 0
    if(myrank.eq.1)then
       call cuda_view_FP64(nnz,method%C_crsMat%crsA)
       call cuda_view_PARCEL_INT(nnz,method%C_crsMat%crsCol)
       call cuda_view_PARCEL_INT(n+1,method%C_crsMat%crsRow_ptr)
    endif
    call cudaDeviceSynchronize() 
#endif

    call cusparse_init( &
         method%cusparce_solver, &
         ! デバイス側のCRS 行列 (ローカル座標)
         method%C_crsMat%crsA, &
         method%C_crsMat%crsCol, &
         method%C_crsMat%crsRow_ptr, &
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
    ! 
    ! Neumann Precon - preconで作成
    ! method%C_preconNeuman%Dinv
    ! method%C_preconNeuman%vecx
    ! method%C_preconNeuman%vecy

    ! cusparce Matrix 
    ! call cudaMalloc_PARCEL_int(method%C_crsMat%crsRow_ptr,n+1)

    etime = MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"   @ 1-8 @",etime - stime
    endif
#endif
    stime  = MPI_Wtime()

#if 0
!    write(*,*)crsA_cpu(1:3)
    call cudaDeviceSynchronize() 
    call cuda_view_FP64(10,method%C_x)
    call cuda_view_FP64(10,method%C_b)
    call cuda_view_FP64(10,method%C_crsMat%crsA)
    call cuda_view_PARCEL_INT(10,method%C_crsMat%crsCol)
    call cuda_view_PARCEL_INT(10,method%C_crsMat%crsRow_ptr)
    call cudaDeviceSynchronize() 
#endif

    etime0 = MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"@ 1 @",etime0 - stime0
    endif
#endif

  end subroutine set_KSP_CRS_cuda

  subroutine set_KSP_MPI_malloc_cuda(&
       method &
       )
    use krylov,only : KSP
    use cu_Vec,only : zero_init_FP64_GPU

    type(KSP) :: method
    integer(kind=intType) i,j
    integer(kind=intType) parcel_int
    integer(kind=intType),parameter :: offset_zero=0

    ! MPI buffer

    ! 受信バッファ
    parcel_int = method%commVec%nrk_data + 1
    call cudaMalloc_fp64(method%C_buff  ,parcel_int)
    call zero_init_FP64_GPU(parcel_int,method%C_buff,offset_zero)

    ! 送信バッファ
    parcel_int = method%commVec%nrk_data_send + 1
    call cudaMalloc_fp64(method%C_sbuff ,parcel_int)
    call zero_init_FP64_GPU(parcel_int,method%C_sbuff,offset_zero)

    do i = 1,method%commVec%nrk_rank_send
       parcel_int = method%commVec%nrk_prm_send(i)%ndata
       call cudaMalloc_PARCEL_INT( &
            method%commVec%nrk_prm_send(i)%C_offset, &
            parcel_int)
       parcel_int = method%commVec%nrk_prm_send(i)%ndata
       call cudaMemcpy_h2d_PARCEL_INT(&
            method%commVec%nrk_prm_send(i)%offset, &
            method%commVec%nrk_prm_send(i)%C_offset, &
            parcel_int, &
            offset_zero,offset_zero)
    enddo

#if 0

    do i = 1,method%commVec%nrk_rank_send
       do j = 1, method%commVec%nrk_prm(i)%ndata
          write(*,*) &
               "@ offset @ ", &
               method%commVec%myrank,&
               method%commVec%nrk_prm_send(i)%offset(j)
       enddo
    enddo
    
#endif

#if 0
    write(*,*)method%commVec%myrank,"buffer size = ", &
         method%commVec%nrk_data,&
         method%commVec%nrk_data_send
#endif

    ! MPI コミュニケータ
    call malloc_MPI_Comm(method%commVec%C_comm,1) 
    call MPI_Comm_f2c(method%commVec%C_comm,method%commVec%comm) 

    ! MPI リクエスト,ステータス
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

  end subroutine set_KSP_MPI_malloc_cuda

  subroutine free_KSP_MPI_malloc_cuda(&
       method &
       )
    use krylov,only : KSP
    type(KSP) :: method
    integer(kind=intType) i
    integer(kind=intType) parcel_int
    integer(kind=intType),parameter :: offset_zero=0

    ! MPI
    call cudafree(method%C_buff  )
    call cudafree(method%C_sbuff )

    do i = 1,method%commVec%nrk_rank_send
       call cudafree(method%commVec%nrk_prm_send(i)%C_offset)
    enddo

    ! MPI コミュニケータ
    call free(method%commVec%C_comm);

    ! MPI リクエスト,ステータス
    if(method%commVec%nrk_rank_send>0)then
       call free(method%commVec%C_s_sta);
       call free(method%commVec%C_s_req);
    endif
    if(method%commVec%nrk_rank>0)then
       call free(method%commVec%C_r_sta);
       call free(method%commVec%C_r_req);
    endif

  end subroutine free_KSP_MPI_malloc_cuda

  subroutine free_KSP_CRS_cuda(&
       method)

    use krylov
    use cu_krylov
    use comm_kernel,only : free_comm
    implicit none
    type(KSP) :: method
    integer(kind=intType) i

    ! ----------------
    ! free

    ! CRS Matrix
    call cudaFree(method%C_crsMat%crsA)
    call cudaFree(method%C_crsMat%crsRow_ptr)
    call cudaFree(method%C_crsMat%crsCol)

    ! Vector
    call cudafree(method%C_r)
    call cudafree(method%C_x)
    call cudafree(method%C_b)

    call free_KSP_MPI_malloc_cuda(method)

    ! free precon
    select case (method%ipreflag)
    case(KSP_PRECONDITION_NONE)
    case(KSP_PRECONDITION_POINTJACOBI)
       call cudafree(method%C_D   ) ; call cusparseDestroyDnVec( method%C_D    ) 
       call cudafree(method%C_Dinv) ; call cusparseDestroyDnVec( method%C_Dinv ) 

    case(KSP_PRECONDITION_BLOCKJACOBI)
       call cudafree(method%cusparce_preconcrsMatrix%crsA)
       call cudafree(method%cusparce_preconcrsMatrix%crsCol)
       call cudafree(method%cusparce_preconcrsMatrix%crsRow_ptr)

    case(KSP_PRECONDITION_ADDITIVESCHWARZ)
    case(KSP_PRECONDITION_NEUMANN)

       call cudafree(method%cusparce_preconcrsMatrix%crsA)
       call cudafree(method%cusparce_preconcrsMatrix%crsCol)
       call cudafree(method%cusparce_preconcrsMatrix%crsRow_ptr)

       call cudafree(method%C_preconNeuman%vecx) ; call cusparseDestroyDnVec( method%C_preconNeuman%vecx ) 
       call cudafree(method%C_preconNeuman%vecy) ; call cusparseDestroyDnVec( method%C_preconNeuman%vecy ) 
       call cudafree(method%C_preconNeuman%D   ) ; call cusparseDestroyDnVec( method%C_preconNeuman%D    ) 
       call cudafree(method%C_preconNeuman%Dinv) ; call cusparseDestroyDnVec( method%C_preconNeuman%Dinv ) 
       
    end select

    ! free cusparse
    call cusparse_free(method%cusparce_solver)

    ! free communication
    call free_comm(method%commVec)

    ! free cublas
    call free_cublas(method%cublas_solver)

  end subroutine free_KSP_CRS_cuda

  subroutine free_KSP_DIA_cuda(&
       method)

    use krylov
    use cu_krylov
    use comm_kernel,only : free_comm
    implicit none
    type(KSP) :: method
    integer(kind=intType) i

    ! ----------------
    ! free

    ! CRS Matrix
    call cudaFree(method%C_diaMat%diaA)
    call cudaFree(method%C_diaMat%offset)

    ! Vector
    call cudafree(method%C_r)
    call cudafree(method%C_x)
    call cudafree(method%C_b)

    call free_KSP_MPI_malloc_cuda(method)

    ! free precon
    select case (method%ipreflag)
    case(KSP_PRECONDITION_NONE)
    case(KSP_PRECONDITION_POINTJACOBI)
       call cudafree(method%C_D)
       call cudafree(method%C_Dinv)
    case(KSP_PRECONDITION_BLOCKJACOBI)
    case(KSP_PRECONDITION_ADDITIVESCHWARZ)
    case(KSP_PRECONDITION_NEUMANN)

      call cudafree(method%C_preconNeuman%diaMat%diaA)
      call cudafree(method%C_preconNeuman%diaMat%offset)
      call cudafree(method%C_preconNeuman%diaMat%iLs)

       call cudafree(method%C_preconNeuman%vecx) 
       call cudafree(method%C_preconNeuman%vecy) 
       call cudafree(method%C_preconNeuman%D   ) 
       call cudafree(method%C_preconNeuman%Dinv) 
       
    end select

    ! free cusparse
    ! call cusparse_free(method%cusparce_solver)

    ! free communication
    call free_comm(method%commVec)

    ! free cublas
    call free_cublas(method%cublas_solver)

  end subroutine free_KSP_DIA_cuda

end module krylov_kernel_gpu
