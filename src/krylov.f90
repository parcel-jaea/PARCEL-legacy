#define JACOBI_LIKE 0
! #define PRINT_TIME 0
#define PRINT_TIME 0

module krylov_kernel
  use PARCEL_define
  use omp_lib
  use cu_krylov
  implicit none
contains
  subroutine SET_KSP_PARAM(method,input_param)
    use krylov
    implicit none

    type(KSP) :: method
    type(KSP_INPUT_PARAM) :: input_param

    call set_KSP_ipreflag(method,input_param%ipreflag)
    call set_KSP_ILU_method(method,input_param%ILU_method)
    call set_KSP_addL(method,input_param%addL)
    call set_KSP_iflagAS(method,input_param%iflagAS)
    call set_KSP_solver(method,input_param%solver)
    call set_KSP_iovlflag(method,input_param%iovlflag)
    call set_KSP_GMRES_m(method,input_param%GMRES_m)
    call set_KSP_GMRES_GSflag(method,input_param%GMRES_GSflag)
    call set_KSP_CAGMRES_sstep(method,input_param%CAGMRES_sstep)
    call set_KSP_CAGMRES_tstep(method,input_param%CAGMRES_tstep)
    call set_KSP_CAGMRES_basis(method,input_param%CAGMRES_basis)
    call set_KSP_CAGMRES_QRflag(method,input_param%CAGMRES_QRflag)
    call set_KSP_CBCG_kstep(method,input_param%CBCG_kstep)
    call set_KSP_CBCG_Eigenflag(method,input_param%CBCG_Eigenflag)
    call set_KSP_power_method_itrmax(method,input_param%power_method_itrmax)
    call set_KSP_CAARNOLDI_sstep(method,input_param%CAARNOLDI_sstep)
    call set_KSP_CAARNOLDI_tstep(method,input_param%CAARNOLDI_tstep)
    call set_KSP_CAARNOLDI_basis(method,input_param%CAARNOLDI_basis)
    call set_KSP_CAARNOLDI_QRflag(method,input_param%CAARNOLDI_QRflag)

    call set_KSP_Neumanitr(method,input_param%Neuman_itr)

    ! call set_KSP_dataflag(method,input_param%dataflag)\



    ! call set_KSP_calcflag(method,input_param%calcflag)

    call set_KSP_dataflag(method,KSP_DOUBLE_PRECISION)
    call set_KSP_calcflag(method,KSP_DOUBLE_PRECISION)

    call set_KSP_itrmax(method,input_param%itrmax)
    call set_KSP_rtolmax(method,input_param%rtolmax)
    call set_KSP_abstolmax(method,input_param%abstolmax)
    call set_KSP_independ_nvec(method,input_param%independ_nvec)
    call set_KSP_precon_thblock(method,input_param%precon_thblock)
    call set_KSP_nBlock(method,input_param%nBlock)

    call set_PRECISION_A(method,input_param%PRECISION_A)
    call set_PRECISION_b(method,input_param%PRECISION_b)
    call set_PRECISION_x(method,input_param%PRECISION_x)
    call set_PRECISION_PRECONDITION(method,input_param%PRECISION_PRECONDITION)

  end subroutine SET_KSP_PARAM

  subroutine view_KSP_PARAM(method)
    use krylov
    implicit none
    type(KSP) :: method
    write(*,*)
    write(*,*)"----------------------"
    write(*,*)"ipreflag              ",method%ipreflag              
    write(*,*)"ILU_method            ",method%ILU_method            
    write(*,*)"addL                  ",method%addL                  
    write(*,*)"iflagAS               ",method%iflagAS               
    write(*,*)"solver                ",method%solver                
    write(*,*)"iovlflag              ",method%iovlflag              
    write(*,*)"GMRES_m               ",method%GMRES_m               
    write(*,*)"GMRES_GSflag          ",method%GMRES_GSflag          
    write(*,*)"CAGMRES_sstep         ",method%CAGMRES_sstep         
    write(*,*)"CAGMRES_tstep         ",method%CAGMRES_tstep         
    write(*,*)"CAGMRES_basis         ",method%CAGMRES_basis         
    write(*,*)"CAGMRES_QRflag        ",method%CAGMRES_QRflag        
    write(*,*)"CBCG_kstep            ",method%CBCG_kstep            
    write(*,*)"CBCG_Eigenflag        ",method%CBCG_Eigenflag        
    write(*,*)"power_method_itrmax   ",method%power_method_itrmax   
    write(*,*)"CAARNOLDI_sstep       ",method%CAARNOLDI_sstep       
    write(*,*)"CAARNOLDI_tstep       ",method%CAARNOLDI_tstep       
    write(*,*)"CAARNOLDI_basis       ",method%CAARNOLDI_basis       
    write(*,*)"CAARNOLDI_QRflag      ",method%CAARNOLDI_QRflag      

    write(*,*)"itrmax                ",method%itrmax      
    write(*,*)"rtolmax               ",method%rtolmax      
    write(*,*)"abstolmax             ",method%abstolmax      
    write(*,*)"independ_nvec         ",method%independ_nvec      
    write(*,*)"precon_thblock        ",method%precon_thblock    
  
    write(*,*)"PRECISION_A           ",method%PRECISION_A           
    write(*,*)"PRECISION_b           ",method%PRECISION_b           
    write(*,*)"PRECISION_x           ",method%PRECISION_x           
    write(*,*)"PRECISION_PRECONDITION",method%PRECISION_PRECONDITION
    write(*,*)"----------------------"
    write(*,*)

  end subroutine view_KSP_PARAM

  subroutine START_KSP_SOLVER(method)
    use krylov
    use PARCEL_time
    implicit none
    type(KSP) :: method
    integer solver
    integer myrank
    solver=method%solver
    myrank=method%commVec%myrank
    precon_time = 0.0d0
!     if(myrank.eq.0)then
!        write(*,*)"MatFormat = " ,method%prm%MatFormat
!     endif

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
    case(-99)
       call KSP_ARNOLDI(method)
    case(KSP_CHECK_QR)
       call MPK_QR_TEST(method)
    case default
       if(myrank.eq.0)then
          write(*,*)"error :: not support solver = ",solver
       end if
    end select

! #ifdef timer
!        if(myrank.eq.0)then
!           write(*,*)
!           write(*,*)" ----------------------- "
!           write(*,*)" ----- Precon Time ----- "
!           write(*,*)" set_AS_vec = ",precon_time(1)
!           write(*,*)" ILU solve  = ",precon_time(2)
!           write(*,*)" AS_Local   = ",precon_time(3)
!           write(*,*)" xcpy       = ",precon_time(4)
!           write(*,*)" get_AS_vec = ",precon_time(5)
!           write(*,*)" --------------------- "
!        endif
! #endif
    
  end subroutine START_KSP_SOLVER

  subroutine KSP_DIA( &
       icomm, &
       vecx,vecb, &
       n,gn,&
       istart, &
       diaA,offset,nnd,&
       input_param,&
       niter,reshistory,iret,&
       EigenValueReal, &
       EigenValueImag, &
       EigenVector, &
       EigenErrorNorm, &
       MaxEigenID, &
       MinEigenID &
       )

    use krylov
    use mpi
    use krylov_kernel_gpu,only : set_KSP_DIA_cuda,free_KSP_DIA_cuda
    use iso_c_binding

    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    integer(kind=intType) iend ! ベクトルの終了位置
    integer(kind=intType) nnd ! 対角数

#ifdef CUDA_SOLVER
    type(c_ptr) vecx   ! 解ベクトル
    type(c_ptr) vecb   ! 右辺ベクトル
    type(c_ptr) diaA   ! 行列値
    type(c_ptr) offset ! 対角数

    real*8,allocatable :: vecx_cpu(:)
    real*8,allocatable :: vecb_cpu(:)
    real*8,allocatable :: diaA_cpu(:)
    integer(kind=intType),allocatable :: offset_cpu(:)

    real*8,pointer :: vecx_Fptr(:)
    real*8,pointer :: vecb_Fptr(:)
    real*8,pointer :: diaA_Fptr(:)
    integer(kind=intType),pointer :: offset_Fptr(:)

    type(devicemem) :: C_diaA
    type(devicemem) :: C_offset

#else
    real*8 vecx(n) ! 解ベクトル
    real*8 vecb(n) ! 右辺ベクトル
    real*8 diaA(nnd*n)  ! 行列値
    integer(kind=intType) offset(nnd) ! 対角オフセット
#endif

    type(KSP_INPUT_PARAM) :: input_param ! KSPの設定パラメータ
    integer niter
    integer iret 

    real*8,optional :: EigenValueReal(:)
    real*8,optional :: EigenValueImag(:)
    real*8,optional :: EigenVector(:)
    real*8,optional :: EigenErrorNorm(:)
    integer,optional :: MaxEigenID
    integer,optional :: MinEigenID

    real*8 reshistory(input_param%itrmax) ! 残差履歴
    integer precon_thblock
    type(KSP) :: method

    integer(kind=intType) parcel_int
    integer myrank
    integer ierr
    integer i
    integer ipreflag,ILU_method,addL
    real*8 stime,etime
    real*8 initialize_time
    real*8 solver_time
    real*8 finalize_time 
    integer solver
    integer ncount
    integer(kind=intType),parameter :: offset_zero=0

    solver=method%solver

    call MPI_Barrier(icomm,ierr)
    stime=MPI_Wtime()
    call PARCEL_KSP_Default_Setting(method)
    call SET_KSP_PARAM(method,input_param)

    call MPI_Comm_rank( icomm, myrank, ierr )
    if( (method%ipreflag.eq.KSP_PRECONDITION_NEUMANN) .or. &
         (method%ipreflag.eq.KSP_PRECONDITION_NEUMANN_LOCAL) )then
       if( method%flag%machineflag.eq.PARCEL_CPU_FLAG_ON )then
          if(myrank.eq.0)then
             write(*,*)"not support Neuman Preconditione"
          endif
          return 
       endif
    endif

#if 1
    if( &
         (method%ipreflag.eq.KSP_PRECONDITION_ADDITIVESCHWARZ) .or. &
         (method%ipreflag.eq.KSP_PRECONDITION_BLOCKJACOBI)     .or. &
         (method%flag%machineflag.eq.PARCEL_CPU_FLAG_ON) &
         )then
#else
       if(1.eq.1)then
#endif

#ifdef CUDA_SOLVER

       allocate(vecx_cpu(n))
       allocate(vecb_cpu(n))
       parcel_int = nnd*n
       allocate(diaA_cpu(parcel_int))
       allocate(offset_cpu(nnd))

       call cudaMalloc_fp64(method%C_x,n)
       call cudaMalloc_fp64(method%C_b,n)
       parcel_int = nnd*n
       call cudaMalloc_fp64(C_diaA,parcel_int)
       call cudaMalloc_PARCEL_int(C_offset,nnd)

       ! CポインタをFortranポインタに割り当て
       call c_f_pointer(vecx,vecx_Fptr,[n])
       call c_f_pointer(vecb,vecb_Fptr,[n])
       parcel_int = nnd*n
       call c_f_pointer(diaA  ,diaA_Fptr,[parcel_int])
       call c_f_pointer(offset,offset_Fptr,[nnd])

       ! CポインタのデータをGPUメモリにコピー
       call cudaMemcpy_dptr2d_double(vecx_Fptr,method%C_x,n)
       call cudaMemcpy_dptr2d_double(vecb_Fptr,method%C_b,n)
       parcel_int = nnd*n
       call cudaMemcpy_dptr2d_double(diaA_Fptr,C_diaA,parcel_int)
       call cudaMemcpy_dptr2d_parcel_int(offset_Fptr,C_offset,nnd)

       call cudaMemcpy_d2h_fp64(method%C_x,vecx_cpu,n,offset_zero,offset_zero)
       call cudaMemcpy_d2h_fp64(method%C_b,vecb_cpu,n,offset_zero,offset_zero)

       parcel_int = nnd*n
       call cudaMemcpy_d2h_fp64(C_diaA,diaA_cpu,parcel_int,offset_zero,offset_zero)    
       call cudaMemcpy_d2h_PARCEL_INT(C_offset,offset_cpu,nnd,offset_zero,offset_zero)

       call cudaFree(C_diaA)
       call cudaFree(C_offset)

       call set_KSP_DIA(&
            icomm, &
            vecx_cpu,vecb_cpu,&
            n,gn,istart, &
            diaA_cpu,offset_cpu,nnd,&
            method)           

       deallocate(vecb_cpu)
       deallocate(diaA_cpu)
       deallocate(offset_cpu)

#else

       call set_KSP_DIA(&
            icomm, &
            vecx,vecb,&
            n,gn,istart, &
            diaA,offset,nnd,&
            method)           

#endif
       call MPI_Barrier(icomm,ierr)
       etime=MPI_Wtime()
       initialize_time = etime-stime

       call MPI_Barrier(icomm,ierr)
       stime=MPI_Wtime()

       myrank=method%commVec%myrank
       call START_KSP_SOLVER(method)

       call MPI_Barrier(icomm,ierr)
       etime=MPI_Wtime()
       solver_time=etime-stime

       call MPI_Barrier(icomm,ierr)
       stime=MPI_Wtime()
       ! -------
       !$OMP PARALLEL DO PRIVATE(i)
       do i=1,method%niter
          reshistory(i)=method%reshistory(i)
       enddo
       iret=method%iret
       niter=method%niter

#ifdef CUDA_SOLVER
       call cudaMemcpy_d2dptr_double(method%C_x,vecx_Fptr,n)
#else
       call get_vecx_KSP(n,vecx,method)
#endif

       solver=method%solver
       if(20<solver)then
          ncount = size( method%EigenValueReal )
          do i=1,ncount
             EigenValueReal(i) = method%EigenValueReal(i)
             EigenValueImag(i) = method%EigenValueImag(i)
             EigenVector(i) = method%EigenVector(i)
             EigenErrorNorm(i) = method%EigenErrorNorm(i)
          enddo
          MaxEigenID = method%MaxEigenID
          MinEigenID = method%MinEigenID
       endif

       ! -------

#ifdef CUDA_SOLVER
       deallocate(vecx_cpu)
       call cudafree(method%C_x)
       call cudafree(method%C_b)
#endif

       call free_KSP(method)
       call MPI_Barrier(icomm,ierr)
       etime=MPI_Wtime()
       finalize_time=etime-stime

    else

#ifdef CUDA_SOLVER
       call set_KSP_DIA_cuda(&
            icomm, &
            vecx,vecb, &
            n,gn,istart, &
            diaA, &
            offset,nnd, &
            method)

       etime=MPI_Wtime()
       initialize_time = etime-stime

       stime=MPI_Wtime()

       call START_KSP_SOLVER(method)
       etime=MPI_Wtime()
       solver_time=etime-stime

       stime=MPI_Wtime()
       !$OMP PARALLEL DO PRIVATE(i)
       do i=1,method%niter
          reshistory(i)=method%reshistory(i)
       enddo
       iret=method%iret
       niter=method%niter

#if CUDA_SOLVER
       call c_f_pointer(vecx,vecx_Fptr,[n])
       call cudaMemcpy_d2dptr_double(method%C_x,vecx_Fptr,n)
#else
       call get_vecx_KSP(n,vecx,method)
#endif
       ! ------------------------
       call free_KSP_DIA_cuda(method)

       etime=MPI_Wtime()
       finalize_time=etime-stime

       if(method%commVec%myrank.eq.0)then
          write(*,*)"initialize time = ",initialize_time,"[s]"
          write(*,*)"solver time     = ",solver_time,"[s]"
          write(*,*)"finalize time   = ",finalize_time,"[s]"
       endif

       return 
#endif

    endif

#if 0
    if(myrank.eq.0)then
       write(*,*)"initialize time = ",initialize_time,"[s]"
       write(*,*)"solver time     = ",solver_time,"[s]"
       write(*,*)"finalize time   = ",finalize_time,"[s]"
    endif
#endif

  end subroutine KSP_DIA
  
  subroutine KSP_CRS( &
       icomm, &
       vecx,vecb, &
       n,gn,nnz,&
       istart, &
       crsA,crsRow_ptr, crsCol,&
       input_param,&
       niter,reshistory,iret, &
       EigenValueReal, &
       EigenValueImag, &
       EigenVector, &
       EigenErrorNorm, &
       MaxEigenID, &
       MinEigenID &
       )

    use mpi
    use krylov
!    use comm_kernel,only : set_comm_GPUmem,free_comm
    use krylov_kernel_gpu,only : set_KSP_CRS_cuda,free_KSP_CRS_cuda
    use iso_c_binding

    implicit none

    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    integer(kind=intType) iend ! ベクトルの終了位置

#ifdef CUDA_SOLVER
    type(c_ptr) vecx ! 解ベクトル
    type(c_ptr) vecb ! 右辺ベクトル
    type(c_ptr) crsA ! 行列値
    type(c_ptr) crsRow_ptr ! ポインタ 
    type(c_ptr) crsCol     ! 列番号

    real*8,allocatable :: vecx_cpu(:)
    real*8,allocatable :: vecb_cpu(:)
    real*8,allocatable :: crsA_cpu(:)
    integer(kind=intType),allocatable :: crsRow_ptr_cpu(:)
    integer(kind=intType),allocatable :: crsCol_cpu(:)

    real*8,pointer :: vecx_Fptr(:)
    real*8,pointer :: vecb_Fptr(:)
    real*8,pointer :: crsA_Fptr(:)
    integer(kind=intType),pointer :: crsRow_ptr_Fptr(:)
    integer(kind=intType),pointer :: crsCol_Fptr(:)

     type(devicemem) :: C_crsA
     type(devicemem) :: C_crsRow_ptr
     type(devicemem) :: C_crsCol

#else

    real*8 vecx(n) ! 解ベクトル
    real*8 vecb(n) ! 右辺ベクトル
    real*8 crsA(nnz) ! 行列値
    integer(kind=intType) crsRow_ptr(n+1) ! ポインタ 
    integer(kind=intType) crsCol(nnz) ! 列番号

#endif

    !    integer itrmax ! 最大反復回数
    !    real*8 rtolmax ! 相対残差 (収束条件)
    !    real*8 abstolmax  !  絶対残差 (収束条件)
    type(KSP_INPUT_PARAM) :: input_param ! KSPの設定パラメータ

    integer(kind=intType) parcel_int
    integer niter
    integer iret 

    real*8,optional :: EigenValueReal(:)
    real*8,optional :: EigenValueImag(:)
    real*8,optional :: EigenVector(:)
    real*8,optional :: EigenErrorNorm(:)
    integer,optional :: MaxEigenID
    integer,optional :: MinEigenID

    real*8 reshistory(input_param%itrmax) ! 残差履歴

    type(KSP) :: method

    integer myrank
    integer ierr
    integer(kind=intType) i,j
    integer ipreflag,ILU_method,addL

    real*8 stime,etime

    real*8 initialize_time
    real*8 solver_time
    real*8 finalize_time 
    integer precon_thblock
    integer solver
    integer ncount
    integer(kind=intType),parameter :: offset_zero=0

    call MPI_Barrier(icomm,ierr)
    stime=MPI_Wtime()
    call PARCEL_KSP_Default_Setting(method)
    call SET_KSP_PARAM(method,input_param)

    call MPI_Comm_rank( icomm, myrank, ierr )
    if( (method%ipreflag.eq.KSP_PRECONDITION_NEUMANN) .or. &
         (method%ipreflag.eq.KSP_PRECONDITION_NEUMANN_LOCAL) )then
       if( method%flag%machineflag.eq.PARCEL_CPU_FLAG_ON )then
          if(myrank.eq.0)then
             write(*,*)"not support Neuman Preconditione"
          endif
          return 
       endif
    endif
    ! write(*,*)"ipreflag = ",method%ipreflag
    ! call view_KSP_PARAM(method)
#if 1
    if( &
         (method%ipreflag.eq.KSP_PRECONDITION_ADDITIVESCHWARZ) .or. &
         (method%ipreflag.eq.KSP_PRECONDITION_BLOCKJACOBI) .or. &
         (method%flag%machineflag.eq.PARCEL_CPU_FLAG_ON) &
         )then
#else
       if(1.eq.1)then
#endif

#ifdef CUDA_SOLVER

       allocate(vecx_cpu(n))
       allocate(vecb_cpu(n))
       allocate(crsA_cpu(nnz))
       parcel_int = n + 1
       allocate(crsRow_ptr_cpu(parcel_int))
       allocate(crsCol_cpu(nnz))

       call cudaMalloc_fp64(method%C_x,n)
       call cudaMalloc_fp64(method%C_b,n)
       call cudaMalloc_fp64(C_crsA,nnz)
       parcel_int = n + 1
       call cudaMalloc_PARCEL_int(C_crsRow_ptr,parcel_int)
       call cudaMalloc_PARCEL_int(C_crsCol,nnz)

       ! CポインタをFortranポインタに割り当て
       call c_f_pointer(vecx,vecx_Fptr,[n])
       call c_f_pointer(vecb,vecb_Fptr,[n])
       call c_f_pointer(crsA,crsA_Fptr,[nnz])
       parcel_int = n + 1
       call c_f_pointer(crsRow_ptr,crsRow_ptr_Fptr,[parcel_int])
       call c_f_pointer(crsCol,crsCol_Fptr,[nnz])

       ! CポインタのデータをGPUメモリにコピー
       call cudaMemcpy_dptr2d_double(vecx_Fptr,method%C_x,n)
       call cudaMemcpy_dptr2d_double(vecb_Fptr,method%C_b,n)
       
       call cudaMemcpy_dptr2d_double(crsA_Fptr,C_crsA,nnz)
       call cudaMemcpy_dptr2d_parcel_int(crsCol_Fptr,C_crsCol,nnz)
       parcel_int = n + 1
       call cudaMemcpy_dptr2d_parcel_int(crsRow_ptr_Fptr,C_crsRow_ptr,parcel_int)

       ! GPU側のデータをCPU側にコピー
       call cudaMemcpy_d2h_fp64(method%C_x,vecx_cpu,n,offset_zero,offset_zero)
       call cudaMemcpy_d2h_fp64(method%C_b,vecb_cpu,n,offset_zero,offset_zero)

       call cudaMemcpy_d2h_fp64(C_crsA,crsA_cpu,nnz,offset_zero,offset_zero)    
       call cudaMemcpy_d2h_PARCEL_INT(C_crsRow_ptr,crsRow_ptr_cpu,n+1,offset_zero,offset_zero)
       call cudaMemcpy_d2h_PARCEL_INT(C_crsCol,crsCol_cpu,nnz,offset_zero,offset_zero)
       call cudaFree(C_crsA)
       call cudaFree(C_crsRow_ptr)
       call cudaFree(C_crsCol)

       call set_KSP_CRS(&
            icomm, &
            vecx_cpu,vecb_cpu, &
            n,gn,nnz,istart, &
            crsA_cpu,crsRow_ptr_cpu,crsCol_cpu, &
            method)

       deallocate(vecb_cpu)
       deallocate(crsA_cpu)
       deallocate(crsRow_ptr_cpu)
       deallocate(crsCol_cpu)

#else
       call set_KSP_CRS(&
            icomm, &
            vecx,vecb, &
            n,gn,nnz,istart, &
            crsA,crsRow_ptr,crsCol, &
            method)
#endif

       call MPI_Barrier(icomm,ierr)
       etime=MPI_Wtime()
       initialize_time = etime-stime

       call MPI_Barrier(icomm,ierr)
       stime=MPI_Wtime()
       myrank=method%commVec%myrank

       call START_KSP_SOLVER(method)
       call MPI_Barrier(icomm,ierr)
       etime=MPI_Wtime()
       solver_time=etime-stime

       call MPI_Barrier(icomm,ierr)
       stime=MPI_Wtime()

       !$OMP PARALLEL DO PRIVATE(i)
       do i=1,method%niter
          reshistory(i)=method%reshistory(i)
       enddo
       iret=method%iret
       niter=method%niter
#ifdef CUDA_SOLVER
       call cudaMemcpy_d2dptr_double(method%C_x,vecx_Fptr,n)
#else
       call get_vecx_KSP(n,vecx,method)
#endif
       solver=method%solver
       if(20<solver)then
          call get_EigenPair_count_KSP(ncount,method)
          call get_EigenPair_KSP( &
               ncount,&
               EigenValueReal,EigenValueImag, &
               EigenVector, &
               EigenErrorNorm, &
               MaxEigenID,MinEigenID,method)
       endif

       ! -------

#ifdef CUDA_SOLVER
       deallocate(vecx_cpu)
       call cudafree(method%C_x)
       call cudafree(method%C_b)
#endif

       call free_KSP(method)
       call MPI_Barrier(icomm,ierr)
       etime=MPI_Wtime()

       finalize_time=etime-stime

#if 0
       if(myrank.eq.0)then
          write(*,*)"initialize time = ",initialize_time,"[s]"
          write(*,*)"solver time     = ",solver_time,"[s]"
          write(*,*)"finalize time   = ",finalize_time,"[s]"
       endif
#endif
       return 

    else

#ifdef CUDA_SOLVER
       ! GPU側の入力
       call set_KSP_CRS_cuda(&
            icomm, &
            vecx,vecb, &
            n,gn,nnz,istart, &
            crsA, &
            crsRow_ptr, &
            crsCol, &
            method)

       call START_KSP_SOLVER(method)

       !$OMP PARALLEL DO PRIVATE(i)
       do i=1,method%niter
          reshistory(i)=method%reshistory(i)
       enddo
       iret=method%iret
       niter=method%niter

#if CUDA_SOLVER
       call c_f_pointer(vecx,vecx_Fptr,[n])
       call cudaMemcpy_d2dptr_double(method%C_x,vecx_Fptr,n)
#else
       call get_vecx_KSP(n,vecx,method)
#endif

       ! ------------------------
       call free_KSP_CRS_cuda(method)
       return 
#endif

    endif

  end subroutine KSP_CRS

  subroutine set_KSP_CRS(&
       icomm, &
       vecx,vecb, &
       n,gn,nnz,istart, &
       crsA,crsRow_ptr,crsCol, &
!       itrmax,&
!       itrmax,rtolmax,abstolmax,&
!       ipreflag,ILU_method,addL,precon_thblock, &
       method)

    use krylov
    use cu_krylov
    use matrix_kernel
    use common_kernel

    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    integer(kind=intType) iend ! ベクトルの終了位置
    real*8 vecx(n) ! 解ベクトル
    real*8 vecb(n) ! 右辺ベクトル
    real*8 crsA(nnz) ! 行列値
    integer(kind=intType) crsRow_ptr(n+1) ! ポインタ 
    integer(kind=intType) crsCol(nnz) ! 列番号
    integer itrmax ! 最大反復回数
!    integer ipreflag ! 前処理方法　
!    integer ILU_method ! 不完全LU分解(D-ILU分解)の対角行列の決定方法
!    integer addL ! 加法シュワルツのオーバラップ領域の大きさ
!    integer iflagAS ! 加法シュワルツのオーバラップ領域補間方法
!    integer precon_thblock !  BJ,AS前処理のプロセス内ブロック数
    type(KSP) :: method
    integer i
    integer ierr
    real*8 stime,etime

    ! ------
    ! CRS Format
    iend  =istart + n -1
    method%prm%MatFormat  = MATRIX_FORMAT_CRS 

    allocate(method%reshistory(max(method%itrmax+100,100)))
    !$OMP PARALLEL DO PRIVATE(i)
    do i=1,max(method%itrmax+100,100)
       method%reshistory(i)=0.0d0
    enddo

    call MPI_Barrier(icomm,ierr)
    stime=MPI_Wtime()

    call init_KSP_CRS(&
         icomm, &
         n,gn,nnz,istart,iend, &
         crsA, crsRow_ptr,crsCol,&
         method)

    call MPI_Barrier(icomm,ierr)
    etime=MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"@ 1 @",etime-stime
    endif
#endif
    stime=MPI_Wtime()
    call set_vecx_KSP(n,vecx,method)
    call set_vecb_KSP(n,vecb,method)
    call MPI_Barrier(icomm,ierr)
    etime=MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"@ 2 @",etime-stime
    endif
#endif

    call get_maxthreads(method%max_numth)

    call MPI_Barrier(icomm,ierr)
    stime=MPI_Wtime()
    call set_precon_CRS(method%ipreflag,method)

    call MPI_Barrier(icomm,ierr)
    etime=MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"@ 3 @",etime-stime
    endif
#endif

    ! --------    
    ! GPU計算が有効ならGPU側のメモリを確保する
    call MPI_Barrier(icomm,ierr)
    stime=MPI_Wtime()
    if(method%GPUflag .eq. PARCEL_GPU_FLAG_ON)then
       call init_GPU(method)
    endif

    call MPI_Barrier(icomm,ierr)
    etime=MPI_Wtime()
#if PRINT_TIME
    if(method%commVec%myrank.eq.0)then
       write(*,*)"@ 4 @",etime-stime
    endif
#endif

  end subroutine set_KSP_CRS
    
  subroutine set_KSP_DIA(&
       icomm, &
       vecx,vecb, &
       n,gn,istart, &
       diaA,offset,nnd,&
!       itrmax,rtolmax,abstolmax, &
!       itrmax,ipreflag,ILU_method,addL,precon_thblock, &
       method)
    use krylov
    use matrix_kernel
    use comm_kernel

    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    integer(kind=intType) iend ! ベクトルの終了位置
    real*8 vecx(n) ! 解ベクトル
    real*8 vecb(n) ! 右辺ベクトル
    integer(kind=intType) nnd
    real*8 diaA(nnd*n)
    integer(kind=intType) offset(nnd)
!    integer itrmax ! 最大反復回数
!    integer ipreflag ! 前処理方法　
!    integer ILU_method ! 不完全LU分解(D-ILU分解)の対角行列の決定方法
!    integer addL ! 加法シュワルツのオーバラップ領域の大きさ
!    integer iflagAS ! 加法シュワルツのオーバラップ領域補間方法
!    integer precon_thblock !  BJ,AS前処理のプロセス内ブロック数
    type(KSP) :: method
    integer i

    ! ------
    ! DIA Format
    iend  =istart + n -1
    method%prm%MatFormat  = MATRIX_FORMAT_DIA 
    allocate(method%reshistory(max(method%itrmax+100,100)))

    !$OMP PARALLEL DO PRIVATE(i)
    do i=1,max(method%itrmax+100,100)
       method%reshistory(i)=0.0d0
    enddo

    call get_nnz_dia(diaA,offset,n,istart,nnd,nnz)
    call init_KSP_DIA(&
         icomm, &
         n,gn,nnz,istart,iend, &
         diaA,offset,nnd,&
         method)
    call set_vecx_KSP(n,vecx,method)
    call set_vecb_KSP(n,vecb,method)

    call get_maxthreads(method%max_numth)
    call set_precon_DIA(method%ipreflag,method)

    if(method%GPUflag .eq. PARCEL_GPU_FLAG_ON)then
       call init_GPU(method)
    endif

  end subroutine set_KSP_DIA

#if 0
  subroutine reset_CRS_Mat(&
       itrmax, &
       ipreflag, &
       ILU_method, &
       addL,&
       crsA,&
       method)

    use krylov
    use matrix_kernel
    use vec_kernel
    implicit none
    type(KSP) :: method
    integer itrmax
    integer ipreflag
    integer ipreflag0
    integer ILU_method
    integer addL
    real*8  crsA(method%prm%nnz) ! 行列値
    integer(kind=intType) nnz
    integer(kind=intType) i

    call free_GPU(method)

    ipreflag0 = method%ipreflag
    nnz=method%prm%nnz
    if(method%ipreflag.eq.KSP_PRECONDITION_ADDITIVESCHWARZ)then
       deallocate(method%ivec_AS)
       deallocate(method%ovec_AS)
       call free_crsMat(method%precon_crsMatAS)
       call free_precon_AS(method)
    endif
    call free_ILUregion(method)
    call free_precon_ILU_crs(method)

    if(ipreflag.ne.PARCEL_RESET_INHERIT)then
       method%ipreflag=ipreflag
    endif
    if(itrmax.ne.PARCEL_RESET_INHERIT)then
       method%itrmax=itrmax
       deallocate(method%reshistory)
       allocate(method%reshistory(max(method%itrmax+100,100)))
       !$OMP PARALLEL DO PRIVATE(i)
       do i=1,max(method%itrmax+100,100)
          method%reshistory(i)=0.0d0
       enddo
    endif
    if(ILU_method.ne.PARCEL_RESET_INHERIT)then
       method%ILU_method=ILU_method
    endif
    if(addL.ne.PARCEL_RESET_INHERIT)then
       method%addL=addL
    endif

    call xcpy(nnz,crsA,method%crsMat%crsA)
    call set_precon_CRS(method%ipreflag,method)

    call init_GPU(method)

  end subroutine reset_CRS_Mat

  subroutine reset_DIA_Mat(&
       itrmax, &
       ipreflag, &
       ILU_method, &
       addL,&
       diaMat,&
       method)
    use krylov
    use matrix_kernel
    use vec_kernel
    implicit none
    type(KSP) :: method
    integer itrmax
    integer ipreflag
    integer ipreflag0
    integer ILU_method
    integer addL
    integer(kind=intType) n,nnd
    integer(kind=intType) i
    real*8 diaMat(method%prm%nnd * method%prm%n)

    call free_GPU(method)
    call free_crsMat(method%precon_crsMatAS)

    ipreflag0 = method%ipreflag
    n=method%prm%n
    nnd=method%prm%nnd
    if(method%ipreflag.eq.KSP_PRECONDITION_ADDITIVESCHWARZ)then
       deallocate(method%ivec_AS)
       deallocate(method%ovec_AS)
       call free_diaMat(method%precon_diaMatAS)
       call free_crsMat(method%precon_crsMatAS)
       call free_precon_AS(method)
    endif

    call free_ILUregion(method)
    call free_precon_ILU_dia(method)

    if(ipreflag.ne.PARCEL_RESET_INHERIT)then
       method%ipreflag=ipreflag
    endif
    if(itrmax.ne.PARCEL_RESET_INHERIT)then
       method%itrmax=itrmax
       deallocate(method%reshistory)
       allocate(method%reshistory(itrmax+100))
       !$OMP PARALLEL DO PRIVATE(i)
       do i=1,itrmax
          method%reshistory(i)=0.0d0
       enddo
    endif
    if(ILU_method.ne.PARCEL_RESET_INHERIT)then
       method%ILU_method=ILU_method
    endif
    if(addL.ne.PARCEL_RESET_INHERIT)then
       method%addL=addL
    endif

    call xcpy(n*nnd,diaMat,method%diaMat%diaA)
    call set_precon_DIA(method%ipreflag,method)

    call init_GPU(method)

  end subroutine reset_DIA_Mat

#else
  subroutine reset_CRS_Mat(&
       itrmax, &
       ipreflag, &
       ILU_method, &
       addL,&
       crsA,&
       method)

    use krylov
    use matrix_kernel
    use vec_kernel
    implicit none
    type(KSP) :: method
    integer itrmax
    integer ipreflag
    integer ipreflag0
    integer ILU_method
    integer addL
    real*8  crsA(method%prm%nnz) ! 行列値
    integer(kind=intType) nnz
    integer(kind=intType) i

    ipreflag0 = method%ipreflag
    nnz=method%prm%nnz
    if(method%ipreflag.eq.KSP_PRECONDITION_ADDITIVESCHWARZ)then
       deallocate(method%ivec_AS)
       deallocate(method%ovec_AS)
       call free_crsMat(method%precon_crsMatAS)
       call free_precon_AS(method)
    endif
    call free_ILUregion(method)
    call free_precon_ILU_crs(method)

    if(ipreflag.ne.PARCEL_RESET_INHERIT)then
       method%ipreflag=ipreflag
    endif
    if(itrmax.ne.PARCEL_RESET_INHERIT)then
       method%itrmax=itrmax
       deallocate(method%reshistory)
       allocate(method%reshistory(max(method%itrmax+100,100)))
       !$OMP PARALLEL DO PRIVATE(i)
       do i=1,max(method%itrmax+100,100)
          method%reshistory(i)=0.0d0
       enddo
    endif
    if(ILU_method.ne.PARCEL_RESET_INHERIT)then
       method%ILU_method=ILU_method
    endif
    if(addL.ne.PARCEL_RESET_INHERIT)then
       method%addL=addL
    endif

    call xcpy(nnz,crsA,method%crsMat%crsA)
    call set_precon_CRS(method%ipreflag,method)

    if(method%GPUflag .eq. PARCEL_GPU_FLAG_ON)then
       call reset_GPU(method,ipreflag0)
    endif

  end subroutine reset_CRS_Mat

  subroutine reset_DIA_Mat(&
       itrmax, &
       ipreflag, &
       ILU_method, &
       addL,&
       diaMat,&
       method)
    use krylov
    use matrix_kernel
    use vec_kernel
    implicit none
    type(KSP) :: method
    integer itrmax
    integer ipreflag
    integer ipreflag0
    integer ILU_method
    integer addL
    integer(kind=intType) n,nnd
    integer(kind=intType) i
    real*8 diaMat(method%prm%nnd * method%prm%n)

    if(method%GPUflag .eq. PARCEL_GPU_FLAG_ON)then
       call free_crsMat(method%precon_crsMatAS)
    endif
    call cudaDeviceSynchronize() 

    ipreflag0 = method%ipreflag
    n=method%prm%n
    nnd=method%prm%nnd
    if(method%ipreflag.eq.KSP_PRECONDITION_ADDITIVESCHWARZ)then
       deallocate(method%ivec_AS)
       deallocate(method%ovec_AS)
       call free_diaMat(method%precon_diaMatAS)
       call free_crsMat(method%precon_crsMatAS)
       call free_precon_AS(method)
    endif

    call free_ILUregion(method)
    call free_precon_ILU_dia(method)

    if(ipreflag.ne.PARCEL_RESET_INHERIT)then
       method%ipreflag=ipreflag
    endif
    if(itrmax.ne.PARCEL_RESET_INHERIT)then
       method%itrmax=itrmax
       deallocate(method%reshistory)
       allocate(method%reshistory(itrmax+100))
       !$OMP PARALLEL DO PRIVATE(i)
       do i=1,itrmax
          method%reshistory(i)=0.0d0
       enddo
    endif
    if(ILU_method.ne.PARCEL_RESET_INHERIT)then
       method%ILU_method=ILU_method
    endif
    if(addL.ne.PARCEL_RESET_INHERIT)then
       method%addL=addL
    endif

    call xcpy(n*nnd,diaMat,method%diaMat%diaA)
    call set_precon_DIA(method%ipreflag,method)

    call cudaDeviceSynchronize() 
    if(method%GPUflag .eq. PARCEL_GPU_FLAG_ON)then
       call reset_GPU(method,ipreflag0)
    endif

  end subroutine reset_DIA_Mat
#endif
  subroutine set_precon_DIA(ipreflag,method)
    use krylov
    use matrix_kernel
    use vec_kernel
    implicit none
    integer ipreflag
    type(KSP) :: method
    integer(kind=intType) n,nnz,istart,nnd
    integer(kind=intType) n_nnd

    n      = method%prm%n
    istart = method%prm%istart
    nnd    = method%prm%nnd
    nnz    = method%prm%nnz

    method%prm_AS%MatFormat=MATRIX_FORMAT_DIA
!     write(*,*)ipreflag,method%addL,nnz
!    method%independ_nvec=CONST_nvec
!    method%precon_thblock=CONST_thblock

    method%prm_AS%MatFormat=MATRIX_FORMAT_CRS

    select case (ipreflag)
    case(KSP_PRECONDITION_NONE)

    case(KSP_PRECONDITION_POINTJACOBI)
       call get_DInv_dia(method%prm,method%diaMat,method%Dinv)
    case(KSP_PRECONDITION_BLOCKJACOBI)
       if(method%ILU_method.eq.KSP_PRECONDITION_ILU)then
          call get_ILU_LU_dia(method,method%diaMat)
       else
          call get_DILU_DInv_dia( &
               method%prm, &
               method%diaMat, &
               method%Dinv, &
               method%ILU_method, &
               method%independ_nvec, &
               method%precon_thblock)
       endif
    case(KSP_PRECONDITION_ADDITIVESCHWARZ)
       call init_crsMat(method%crsMat,n,nnz)
       call dia2crs(method%diaMat,method%crsMat,n,istart,nnd,nnz)
       call get_DILU_DInv_AS_crs(method,method%crsMat)

       if(method%addL.ne.0)then
          call get_nnd_crs_without_sort(&
               method%precon_diaMatAS, &
               method%prm_AS%n,&
               method%prm_AS%nnz,&
               method%precon_crsMatAS%crsCol,&
               method%precon_crsMatAS%crsRow_ptr,&
               method%prm_AS%istart,&
               method%prm_AS%iend,&
               method%prm_AS%nnd &
               )
          call crs2dia_without_sort( &
               method%precon_crsMatAS, &
               method%precon_diaMatAS, &
               method%prm_AS%nnd, &
               method%prm_AS%n, &
               method%prm_AS%nnz,&
               method%prm_AS%istart)
          ! addLが0の場合はブロックヤコビ前処理なので加法シュヴァルツ前処理用の行列は作成しない
          ! CUDAの場合この結果が必要
          if(method%GPUflag .eq. PARCEL_GPU_FLAG_ON)then

          else
             call free_crsMat(method%precon_crsMatAS)
          endif
       endif
       call free_crsMat(method%crsMat)

    case(KSP_PRECONDITION_NEUMANN)
       call xcpy_diaMat(method,method%diaMat,method%preconNeuman%diaMat)

       if(allocated(method%preconNeuman%D))deallocate(method%preconNeuman%D)
       if(allocated(method%preconNeuman%Dinv))deallocate(method%preconNeuman%Dinv)
       if(allocated(method%preconNeuman%vecx))deallocate(method%preconNeuman%vecx)
       if(allocated(method%preconNeuman%vecy))deallocate(method%preconNeuman%vecy)

       allocate(method%preconNeuman%D(n));
       allocate(method%preconNeuman%Dinv(n));
       allocate(method%preconNeuman%vecx(n)); 
       allocate(method%preconNeuman%vecy(n)); 

       call zero_initialize(n,method%preconNeuman%D)
       call zero_initialize(n,method%preconNeuman%Dinv)
       call zero_initialize(n,method%preconNeuman%vecx)
       call zero_initialize(n,method%preconNeuman%vecy)

       ! A = I - F 
       ! (I - F)^-1 = I + F + F^2 + F^3 + ....
       ! Fを計算する 
       call Neumann_scaling_dia( &
            method%prm, &
            method%preconNeuman%diaMat, &
            method%preconNeuman%D , &
            method%preconNeuman%Dinv)

    case default
       if(method%commVec%myrank.eq.0)then
          write(*,*)"error :: not support preconditioning = ",ipreflag
       end if
    end select

    select case (method%ipreflag)
    case(KSP_PRECONDITION_BLOCKJACOBI)
       call init_blockJacobi_dia(method)
    case(KSP_PRECONDITION_ADDITIVESCHWARZ)
       call init_AdditiveSchwarz_dia(method)
    end select

  end subroutine set_precon_DIA

  subroutine set_precon_CRS(ipreflag,method)
    use krylov
    use vec_kernel,only : zero_initialize
    use matrix_kernel
    implicit none
    integer ipreflag
    type(KSP) :: method
    integer(kind=intType) n
    integer(kind=intType) nnz

    method%prm_AS%MatFormat=MATRIX_FORMAT_CRS

!    method%independ_nvec=CONST_nvec
!    method%precon_thblock=CONST_thblock

    select case (ipreflag)
    case(KSP_PRECONDITION_NONE)

    case(KSP_PRECONDITION_POINTJACOBI)
       call get_DInv_crs(method%prm,method%crsMat,method%Dinv)
    case(KSP_PRECONDITION_BLOCKJACOBI)
       if(method%ILU_method.eq.KSP_PRECONDITION_ILU)then
          call get_ILU_LU_crs(method,method%crsMat)
       else
          call get_DILU_DInv_crs(&
               method%prm, &
               method%crsMat, &
               method%Dinv, &
               method%ILU_method, &
               method%independ_nvec, &
               method%precon_thblock)
       endif
    case(KSP_PRECONDITION_ADDITIVESCHWARZ)
       call get_DILU_DInv_AS_crs(method,method%crsMat)

    case(KSP_PRECONDITION_NEUMANN)

       n = method%prm%n
       call xcpy_crsMat(method,method%crsMat,method%preconNeuman%crsMat)

       if(allocated(method%preconNeuman%D))deallocate(method%preconNeuman%D)
       if(allocated(method%preconNeuman%Dinv))deallocate(method%preconNeuman%Dinv)
       if(allocated(method%preconNeuman%vecx))deallocate(method%preconNeuman%vecx)
       if(allocated(method%preconNeuman%vecy))deallocate(method%preconNeuman%vecy)

       allocate(method%preconNeuman%D(n));
       allocate(method%preconNeuman%Dinv(n));
       allocate(method%preconNeuman%vecx(n)); 
       allocate(method%preconNeuman%vecy(n)); 

       call zero_initialize(n,method%preconNeuman%D)
       call zero_initialize(n,method%preconNeuman%Dinv)
       call zero_initialize(n,method%preconNeuman%vecx)
       call zero_initialize(n,method%preconNeuman%vecy)

#if JACOBI_LIKE 
       call init_jacobilike_crs( &
            method%prm, &
            method%preconNeuman%crsMat, &
            method%preconNeuman%D , &
            method%preconNeuman%Dinv)
#else
       ! A = I - F 
       ! (I - F)^-1 = I + F + F^2 + F^3 + ....
       ! Fを計算する 
       call Neumann_scaling_crs( &
            method%prm, &
            method%preconNeuman%crsMat, &
            method%preconNeuman%D , &
            method%preconNeuman%Dinv)
#endif

    case default
       if(method%commVec%myrank.eq.0)then
          write(*,*)"error :: not support preconditioning = ",ipreflag
       end if
    end select

    select case (method%ipreflag)
    case(KSP_PRECONDITION_BLOCKJACOBI)
       if(method%ILU_method.eq.KSP_PRECONDITION_ILU)then
          
       else
          call init_blockJacobi_crs(method)
       endif
    case(KSP_PRECONDITION_ADDITIVESCHWARZ)
       call init_AdditiveSchwarz_crs(method)
    end select

  end subroutine set_precon_CRS

  subroutine set_CRS_Mat(n,nnz, crsA,crsRow_ptr,crsCol, crsMat)
    use krylov
    implicit none
    integer(kind=intType) n,nnz
    real*8  crsA(nnz) ! 行列値
    integer(kind=intType) crsRow_ptr(n+1) ! ポインタ 
    integer(kind=intType) crsCol(nnz) ! 列番号
    type(CRS_Mat) :: crsMat
    integer(kind=intType) i

    !$OMP PARALLEL DO PRIVATE(i)  
    do i=1,n+1
       crsMat%crsRow_ptr(i)=crsRow_ptr(i)
    enddo
    !$OMP END PARALLEL DO 

!    !$OMP PARALLEL DO PRIVATE(i)  
    do i=1,nnz
       crsMat%crsA(i)   = crsA(i)
       crsMat%crsCol(i) = crsCol(i)
    enddo
!    !$OMP END PARALLEL DO 

  end subroutine set_CRS_Mat

  subroutine init_KSP_Mat_CRS(crsA,crsRow_ptr,crsCol, method)
    use krylov
    use comm_kernel
    use matrix_kernel
    implicit none
    type(KSP) :: method
    real*8 crsA(method%prm%nnz) ! 行列値
    integer(kind=intType) crsRow_ptr(method%prm%n+1) ! ポインタ 
    integer(kind=intType) crsCol(method%prm%nnz) ! 列番号

    integer(kind=intType) i
    integer(kind=intType) n,gn,nnz,istart,iend
    integer icomm

    integer(kind=intType) nCol,nColL,nColU,nColD
    integer myrank,ierr
    integer(kind=intType),allocatable :: Colsort(:)

    real*8 ts,te
    n     =method%prm%n
    gn    =method%prm%gn
    nnz   =method%prm%nnz
    istart=method%prm%istart
    iend  =method%prm%iend
    icomm  =method%commVec%comm

    call MPI_Comm_rank(icomm,myrank,ierr)

    call MPI_Barrier(icomm,ierr)
    ts=MPI_Wtime()
    allocate( Colsort(nnz) )
    call init_crsMat(method%crsMat,n,nnz)
    call MPI_Barrier(icomm,ierr)
    te=MPI_Wtime()
#if PRINT_TIME
     if(myrank.eq.0)then
        write(*,*)"@ 1-1 @",te-ts
     endif
#endif
!
    call MPI_Barrier(icomm,ierr)
    ts=MPI_Wtime()
    call set_CRS_Mat(n,nnz, crsA,crsRow_ptr,crsCol, method%crsMat)
    call MPI_Barrier(icomm,ierr)
    te=MPI_Wtime()


#if PRINT_TIME
     if(myrank.eq.0)then
        write(*,*)"@ 1-2 @",te-ts
     endif
#endif

    call MPI_Barrier(icomm,ierr)
    ts=MPI_Wtime()
    call get_crsCol_nCol(method%prm, method%crsMat%crsCol, Colsort,&
         nCol,nColL,nColU,nColD)
    call MPI_Barrier(icomm,ierr)
    te=MPI_Wtime()
#if PRINT_TIME
     if(myrank.eq.0)then
        write(*,*)"@ 1-3 @",te-ts
     endif
#endif
    method%crsMat%nColL = nColL
    method%crsMat%nColU = nColU

    call MPI_Barrier(icomm,ierr)
    ts=MPI_Wtime()
    call get_crs_Row(method%prm,method%crsMat,nCol,Colsort)
    call MPI_Barrier(icomm,ierr)
    te=MPI_Wtime()

#if PRINT_TIME
    if(myrank.eq.0)then
       write(*,*)"@ 1-4 @",te-ts
    endif
#endif

    call MPI_Barrier(icomm,ierr)
    ts=MPI_Wtime()
    call set_comm(method%commVec,icomm,method%prm,Colsort,nCol)
    call MPI_Barrier(icomm,ierr)
    te=MPI_Wtime()
#if PRINT_TIME
    if(myrank.eq.0)then
       write(*,*)"@ 1-5 @",te-ts
    endif
#endif
    ! ----
   !  call MPI_Barrier(icomm,ierr)
    ts=MPI_Wtime()
!    write(*,*)" @ nnz @", nnz
    allocate(method%crsMat%iLs(nnz))
    call init_MatVec_crs(method%prm,method%crsMat,method%nBlock)
    if(method%iovlflag.eq.KSP_SPMV_OVL_RANK)then
       call init_MatVec_crs_ovl2(method%prm,method%crsMat,method%commVec,method%nBlock)
    endif

!     do i=1,method%commVec%nrk_rank
!        allocate(method%commVec%nrk_prm(i)%iLs_Ovl2(nnz))
!        allocate(method%commVec%nrk_prm(i)%iLe_Ovl2(nnz))
!     enddo
    te=MPI_Wtime()
#if PRINT_TIME
     if(myrank.eq.0)then
        write(*,*)"@ 1-6 @",te-ts
     endif
#endif
   !  call MPI_Barrier(icomm,ierr)
    ts=MPI_Wtime()
!    call init_MatVec_crs_ovl2(method%prm,method%crsMat,method%commVec)
   !  call MPI_Barrier(icomm,ierr)
    te=MPI_Wtime()
!     if(myrank.eq.0)then
!        write(*,*)"@ 7 @",te-ts
!     endif
    ! ----

    deallocate( Colsort )

  end subroutine init_KSP_Mat_CRS

  subroutine set_DIA_Mat(n,nnd, diaMat,offset, method)
    use krylov
    use vec_kernel

    implicit none
    integer(kind=intType) n,nnd
    real*8 diaMat(nnd*n)
    integer(kind=intType) offset(nnd)
    type(KSP) :: method
    integer(kind=intType) i
    integer(kind=intType) n_nnd
    n_nnd=n*nnd
    call xcpy(n_nnd,diaMat,method%diaMat%diaA)
    !$OMP PARALLEL DO PRIVATE(i)  
    do i=1,nnd
       method%diaMat%offset(i)=offset(i)
    enddo
    !$OMP END PARALLEL DO 

  end subroutine set_DIA_Mat

  subroutine init_KSP_Mat_DIA(diaA,offset,method)

    use krylov
    use comm_kernel
    use matrix_kernel
    implicit none

    type(KSP) :: method
    real*8 diaA(method%prm%nnd * method%prm%n)
    integer(kind=intType) offset(method%prm%nnd)
    integer(kind=intType),allocatable :: Colsort(:)
    integer icomm
    integer(kind=intType) n,gn,nnz,istart,iend,nnd
    integer(kind=intType) i
    integer(kind=intType) n_nnd
    real*8 stime0,etime0
    real*8 stime,etime
    real*8 rtime(40)

    rtime=0.0d0
    stime0=MPI_Wtime()

    n      = method%prm%n
    gn     = method%prm%gn
    nnz    = method%prm%nnz
    istart = method%prm%istart
    iend   = method%prm%iend
    icomm  = method%commVec%comm
    nnd    = method%prm%nnd

    stime=MPI_Wtime()
    call init_diaMat(method%diaMat,n,nnd)
    etime=MPI_Wtime()
    rtime(1)=etime-stime

    stime=MPI_Wtime()
    call set_DIA_Mat(n,nnd, diaA,offset,method)
    etime=MPI_Wtime()
    rtime(2)=etime-stime

    stime=MPI_Wtime()
    call set_dia_Row(method%prm,method%diaMat)
    etime=MPI_Wtime()
    rtime(3)=etime-stime

    stime=MPI_Wtime()
    deallocate(method%diaMat%Vec)
    allocate(method%diaMat%Vec(method%diaMat%nCol))
    allocate( Colsort(nnz) )
    call get_nnz_Colsort_dia( &
         method%prm, &
         method%diaMat, &
         Colsort, &
         method%diaMat%nCol,  &
         method%diaMat%nColL, &
         method%diaMat%nColU, &
         method%diaMat%nColD  &
         )
    etime=MPI_Wtime()
    rtime(5)=etime-stime

    stime=MPI_Wtime()
    call set_comm(method%commVec, &
         icomm, &
         method%prm, &
         Colsort, &
         method%diaMat%nCol &
         )
    etime=MPI_Wtime()
    rtime(6)=etime-stime

    stime=MPI_Wtime()
    n_nnd=n*nnd
    allocate(method%diaMat%iLs( n_nnd ))
    call init_MatVec_dia(method%prm,method%diaMat,method%commVec,method%nBlock)
    if(method%iovlflag.eq.KSP_SPMV_OVL_RANK)then
       call init_MatVec_dia_ovl2(method%prm,method%diaMat,method%commVec,method%nBlock)
    endif
!    do i=1,method%commVec%nrk_rank
!       allocate(method%commVec%nrk_prm(i)%iLs_Ovl2( n*nnd ) )
!       allocate(method%commVec%nrk_prm(i)%iLe_Ovl2( n*nnd ) )
!    enddo
!    call init_MatVec_dia2(method%prm,method%diaMat,method%commVec)

    deallocate( Colsort )
    etime=MPI_Wtime()
    rtime(7)=etime-stime

    etime0=MPI_Wtime()

! #ifdef timer_init
#if 0
    if(method%commVec%myrank.eq.0)then
       write(*,*)"------------------------------------------"
       write(*,*)"init_KSP_Mat_DIA         ",etime0 - stime0
       write(*,*)"     init_diaMat         ",rtime(1)
       write(*,*)"     set_DIA_Mat         ",rtime(2)
       write(*,*)"     set_dia_Row         ",rtime(3)
       write(*,*)"     get_nnz_Colsort_dia ",rtime(5)
       write(*,*)"     set_comm            ",rtime(6)
       write(*,*)"     init_MatVec         ",rtime(7)
       write(*,*)"------------------------------------------"
    endif
#endif

  end subroutine init_KSP_Mat_DIA

  subroutine init_KSP_CRS(&
       icomm, &
       n,gn,nnz,&
       istart,iend, &
       crsA,crsRow_ptr,&
       crsCol, method)
    use krylov
    use comm_kernel
    use matrix_kernel
    use vec_kernel
    implicit none
    integer icomm
    integer(kind=intType) n,gn,nnz,istart,iend
    real*8  crsA(nnz) ! 行列値
    integer(kind=intType) crsRow_ptr(n+1) ! ポインタ 
    integer(kind=intType) crsCol(nnz) ! 列番号

    type(KSP) :: method
    integer(kind=intType) i

    method%commVec%comm=icomm
    method%prm%n=n
    method%prm%gn=gn
    method%prm%nnz=nnz
    method%prm%istart=istart
    method%prm%iend=iend

    allocate(method%r(n))
    allocate(method%x(n))
    allocate(method%b(n))
    allocate(method%Dinv(n))

    ! do i=1,n
    !    method%r(i)=0.0d0
    !    method%x(i)=0.0d0
    !    method%b(i)=0.0d0     
    ! enddo
    call zero_initialize(n,method%r)
    call zero_initialize(n,method%x)
    call zero_initialize(n,method%b)
    call zero_initialize(n,method%Dinv)

    call init_KSP_Mat_CRS(crsA,crsRow_ptr,crsCol, method)

!    allocate(method%buff(method%commVec%nrk_data))

    allocate(method%buff(method%commVec%nrk_data+1))
    allocate(method%rbuffL(method%crsMat%nColL))
    allocate(method%rbuffU(method%crsMat%nColU))

    !$OMP PARALLEL DO PRIVATE(i)  
!    do i=1,method%commVec%nrk_data
    do i=1,method%commVec%nrk_data+1
       method%buff(i)=0.0d0
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i)  
    do i=1,method%crsMat%nColL
       method%rbuffL(i)=0.0d0
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i)  
    do i=1,method%crsMat%nColU
       method%rbuffU(i)=0.0d0
    enddo
    !$OMP END PARALLEL DO

  end subroutine init_KSP_CRS

  subroutine init_KSP_DIA(&
       icomm, &
       n,gn,&
       nnz,istart,iend, &
       diaA,offset,nnd,&
       method)
    use krylov
    use comm_kernel
    use matrix_kernel
    use vec_kernel
    implicit none
    integer icomm
    integer(kind=intType) n,gn,nnz,istart,iend
    integer(kind=intType) nnd
    real*8  diaA(nnd*n) ! 行列値
    integer(kind=intType) offset(nnd)
    type(KSP) :: method
    integer(kind=intType) i

    method%commVec%comm=icomm
    method%prm%n=n
    method%prm%gn=gn
    method%prm%nnz=nnz
    method%prm%istart=istart
    method%prm%iend=iend
    method%prm%nnd=nnd

    allocate(method%r(n))
    allocate(method%x(n))
    allocate(method%b(n))
    allocate(method%Dinv(n))

    call zero_initialize(n,method%r)
    call zero_initialize(n,method%x)
    call zero_initialize(n,method%b)
    call zero_initialize(n,method%Dinv)
    call init_KSP_Mat_DIA(diaA,offset,method)

     !     method%commVec%nrk_data, &
     !     method%diaMat%nColL, &
     !     method%diaMat%nColU

!    allocate(method%buff(method%commVec%nrk_data))
    allocate(method%buff(method%commVec%nrk_data+1))
    allocate(method%rbuffL(method%diaMat%nColL))
    allocate(method%rbuffU(method%diaMat%nColU))
       
    !$OMP PARALLEL DO PRIVATE(i)  
    do i=1,method%commVec%nrk_data+1
!     do i=1,method%commVec%nrk_data
       method%buff(i)=0.0d0
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i)  
    do i=1,method%diaMat%nColL
       method%rbuffL(i)=0.0d0
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i)  
    do i=1,method%diaMat%nColU
       method%rbuffU(i)=0.0d0
    enddo
    !$OMP END PARALLEL DO
  end subroutine init_KSP_DIA

  subroutine free_KSP(method)
    use krylov
    use cu_krylov
    use matrix_kernel,only : free_crsMat
    implicit none
    type(KSP) :: method
    integer MatFormat 
    MatFormat = method%prm%MatFormat
    if(method%GPUflag .eq. PARCEL_GPU_FLAG_ON)then
       call free_GPU(method)
       select case (MatFormat)
       case(MATRIX_FORMAT_CRS)

       case(MATRIX_FORMAT_DIA)
          call free_crsMat(method%precon_crsMatAS)
       end select
    endif
    
    select case (MatFormat)
    case(MATRIX_FORMAT_CRS)
       call free_KSP_CRS(method)
    case(MATRIX_FORMAT_DIA)
       call free_KSP_DIA(method)
    end select

    if(allocated(method%buff))deallocate(method%buff)
    if(allocated(method%buff_hi))deallocate(method%buff_hi)
    if(allocated(method%buff_lo))deallocate(method%buff_lo)
    if(allocated(method%rbuffL))deallocate(method%rbuffL)
    if(allocated(method%rbuffL_hi))deallocate(method%rbuffL_hi)
    if(allocated(method%rbuffL_lo))deallocate(method%rbuffL_lo)
    if(allocated(method%rbuffU))deallocate(method%rbuffU)
    if(allocated(method%rbuffU_hi))deallocate(method%rbuffU_hi)
    if(allocated(method%rbuffU_lo))deallocate(method%rbuffU_lo)
    if(allocated(method%EigenValueReal))deallocate(method%EigenValueReal)
    if(allocated(method%EigenValueImag))deallocate(method%EigenValueImag)
    if(allocated(method%EigenVector))deallocate(method%EigenVector)
    if(allocated(method%EigenErrorNorm))deallocate(method%EigenErrorNorm)
    if(allocated(method%r))deallocate(method%r)
    if(allocated(method%r_hi))deallocate(method%r_hi)
    if(allocated(method%r_lo))deallocate(method%r_lo)
    if(allocated(method%x))deallocate(method%x)
    if(allocated(method%x_hi))deallocate(method%x_hi)
    if(allocated(method%x_lo))deallocate(method%x_lo)
    if(allocated(method%b))deallocate(method%b)
    if(allocated(method%b_hi))deallocate(method%b_hi)
    if(allocated(method%b_lo))deallocate(method%b_lo)
    if(allocated(method%reshistory))deallocate(method%reshistory)
    if(allocated(method%Dinv))deallocate(method%Dinv)
    if(allocated(method%Dinv_hi))deallocate(method%Dinv_hi)
    if(allocated(method%Dinv_lo))deallocate(method%Dinv_lo)
    if(allocated(method%ILUregionMat))deallocate(method%ILUregionMat)
    if(allocated(method%ivec_AS))deallocate(method%ivec_AS)
    if(allocated(method%ivec_AS_hi))deallocate(method%ivec_AS_hi)
    if(allocated(method%ivec_AS_lo))deallocate(method%ivec_AS_lo)
    if(allocated(method%ovec_AS))deallocate(method%ovec_AS)
    if(allocated(method%ovec_AS_hi))deallocate(method%ovec_AS_hi)
    if(allocated(method%ovec_AS_lo))deallocate(method%ovec_AS_lo)
    if(allocated(method%preconILU))deallocate(method%preconILU)
    if(allocated(method%preconAS))deallocate(method%preconAS)

  end subroutine free_KSP

  subroutine free_ILUregion(method)
    use krylov
    use common_kernel
    implicit none
    type(KSP) :: method
    integer i,max_numth
    integer numth
    ! -----
    call get_maxthreads(max_numth)
    call set_precon_maxthreads(method%prm%n,method%precon_thblock)
    call get_maxthreads(numth)

    if(allocated(method%ILUregionMat))then
!       do i=1,max_numth
       do i=1,numth
          ! 配列が割り当てられていれば.true.
          if(allocated(method%ILUregionMat(i)%ILU_Lstart))then
             deallocate(method%ILUregionMat(i)%ILU_Lstart)
          endif
          if(allocated(method%ILUregionMat(i)%ILU_Lend))then
             deallocate(method%ILUregionMat(i)%ILU_Lend)
          endif
          if(allocated(method%ILUregionMat(i)%ILU_Ustart))then
             deallocate(method%ILUregionMat(i)%ILU_Ustart)
          endif
          if(allocated(method%ILUregionMat(i)%ILU_Uend))then
             deallocate(method%ILUregionMat(i)%ILU_Uend)
          endif
       enddo
       deallocate(method%ILUregionMat)
    endif
    call omp_set_num_threads(max_numth)

  end subroutine free_ILUregion
    ! -----

  subroutine free_KSP_CRS(method)
    use krylov
    use matrix_kernel
    use comm_kernel
    implicit none
    type(KSP) :: method

    ! C側ポインタの後片付けを行う    
    ! CPU側のメモリで解放するメモリを指定するためCPUメモリ解放前に解放する
    
    deallocate(method%r)
    deallocate(method%x)
    deallocate(method%b)
    deallocate(method%Dinv)

    deallocate(method%buff)
    deallocate(method%rbuffL)
    deallocate(method%rbuffU)
    deallocate(method%crsMat%iLs)

    call free_comm(method%commVec)
    call free_crsMat(method%crsMat)

    call free_ILUregion(method)
    call free_precon_ILU_crs(method)
    ! -----
    if(method%ipreflag.eq.KSP_PRECONDITION_ADDITIVESCHWARZ)then
       deallocate(method%ivec_AS)
       deallocate(method%ovec_AS)
       call free_crsMat(method%precon_crsMatAS)
       call free_precon_AS(method)
    endif

    if(method%ipreflag.eq.KSP_PRECONDITION_NEUMANN)then
       call free_crsMat(method%preconNeuman%crsMat)
       call free_diaMat(method%preconNeuman%diaMat)

       if(allocated(method%preconNeuman%vecx))then
          deallocate(method%preconNeuman%vecx)
       endif
       if(allocated(method%preconNeuman%vecx_hi))then
          deallocate(method%preconNeuman%vecx_hi)
       endif
       if(allocated(method%preconNeuman%vecx_lo))then
          deallocate(method%preconNeuman%vecx_lo)
       endif

       if(allocated(method%preconNeuman%vecy))then
          deallocate(method%preconNeuman%vecy)
       endif
       if(allocated(method%preconNeuman%vecy_hi))then
          deallocate(method%preconNeuman%vecy_hi)
       endif
       if(allocated(method%preconNeuman%vecy_lo))then
          deallocate(method%preconNeuman%vecy_lo)
       endif

       if(allocated(method%preconNeuman%D))then
          deallocate(method%preconNeuman%D)
       endif
       if(allocated(method%preconNeuman%D_hi))then
          deallocate(method%preconNeuman%D_hi)
       endif
       if(allocated(method%preconNeuman%D_lo))then
          deallocate(method%preconNeuman%D_lo)
       endif

       if(allocated(method%preconNeuman%Dinv))then
          deallocate(method%preconNeuman%Dinv)
       endif
       if(allocated(method%preconNeuman%Dinv_hi))then
          deallocate(method%preconNeuman%Dinv_hi)
       endif
       if(allocated(method%preconNeuman%Dinv_lo))then
          deallocate(method%preconNeuman%Dinv_lo)
       endif

    end if

    if(allocated(method%reshistory))deallocate(method%reshistory)

  end subroutine free_KSP_CRS

  subroutine free_KSP_DIA(method)
    use krylov
    use comm_kernel
    use matrix_kernel
    implicit none
    type(KSP) :: method
    integer i,max_numth

    deallocate(method%r)
    deallocate(method%x)
    deallocate(method%b)
    deallocate(method%Dinv)
    deallocate(method%buff)
    deallocate(method%rbuffL)
    deallocate(method%rbuffU)
    deallocate(method%diaMat%iLs)

    call free_comm(method%commVec)
    call free_diaMat(method%diaMat)

    call free_ILUregion(method)
    call free_precon_ILU_dia(method)
    if(method%ipreflag.eq.KSP_PRECONDITION_ADDITIVESCHWARZ)then
       deallocate(method%ivec_AS)
       deallocate(method%ovec_AS)
       call free_diaMat(method%precon_diaMatAS)
       call free_precon_AS(method)
    endif

    if(method%ipreflag.eq.KSP_PRECONDITION_NEUMANN)then
       call free_crsMat(method%preconNeuman%crsMat)
       call free_diaMat(method%preconNeuman%diaMat)

       if(allocated(method%preconNeuman%vecx))then
          deallocate(method%preconNeuman%vecx)
       endif
       if(allocated(method%preconNeuman%vecx_hi))then
          deallocate(method%preconNeuman%vecx_hi)
       endif
       if(allocated(method%preconNeuman%vecx_lo))then
          deallocate(method%preconNeuman%vecx_lo)
       endif

       if(allocated(method%preconNeuman%vecy))then
          deallocate(method%preconNeuman%vecy)
       endif
       if(allocated(method%preconNeuman%vecy_hi))then
          deallocate(method%preconNeuman%vecy_hi)
       endif
       if(allocated(method%preconNeuman%vecy_lo))then
          deallocate(method%preconNeuman%vecy_lo)
       endif

       if(allocated(method%preconNeuman%D))then
          deallocate(method%preconNeuman%D)
       endif
       if(allocated(method%preconNeuman%D_hi))then
          deallocate(method%preconNeuman%D_hi)
       endif
       if(allocated(method%preconNeuman%D_lo))then
          deallocate(method%preconNeuman%D_lo)
       endif

       if(allocated(method%preconNeuman%Dinv))then
          deallocate(method%preconNeuman%Dinv)
       endif
       if(allocated(method%preconNeuman%Dinv_hi))then
          deallocate(method%preconNeuman%Dinv_hi)
       endif
       if(allocated(method%preconNeuman%Dinv_lo))then
          deallocate(method%preconNeuman%Dinv_lo)
       endif

    end if

    deallocate(method%reshistory)
  end subroutine free_KSP_DIA

  subroutine set_vecb_KSP(n,b, method)
    use krylov
    use vec_kernel
    implicit none
    integer(kind=intType) n
    real*8 b(n)
    type(KSP) :: method

    call xcpy(n,b,method%b)

  end subroutine set_vecb_KSP

  subroutine set_vecx_KSP(n,x, method)
    use krylov
    use vec_kernel
    implicit none
    integer(kind=intType) n
    real*8 x(n)
    type(KSP) :: method

    call xcpy(n,x,method%x)

  end subroutine set_vecx_KSP

  subroutine get_reshistory_KSP(niter,reshistory, method)
    use krylov
    use vec_kernel
    implicit none
    type(KSP) :: method
    integer niter
    real*8 reshistory(method%niter)
    integer(kind=intType) niter_intType

    niter         = method%niter
    niter_intType = method%niter
    call xcpy(niter_intType,method%reshistory,reshistory)

  end subroutine get_reshistory_KSP

  subroutine get_vecx_KSP(n,x, method)
    use krylov
    use vec_kernel
    implicit none
    integer(kind=intType) n
    real*8 x(n)
    type(KSP) :: method
    call xcpy(n,method%x,x)

  end subroutine get_vecx_KSP

  ! r = b-Ax
  subroutine calc_res(method,x,b,r)
    use krylov
    implicit none
    type(KSP) :: method
    real*8 x(method%prm%n)
    real*8 b(method%prm%n)
    real*8 r(method%prm%n)
    integer MatFormat 
    MatFormat = method%prm%MatFormat 
    select case (MatFormat)
    case(MATRIX_FORMAT_CRS)
       call calc_res_CRS(method,x,b,r)
    case(MATRIX_FORMAT_DIA)
       call calc_res_DIA(method,x,b,r)
   case(MATRIX_FORMAT_DDM)
       call calc_res_DDM(method%ddmMat,x,b,r)
    end select
  end subroutine calc_res

  subroutine calc_Matvec(method,x,y)
    use krylov
    implicit none
    type(KSP) :: method
    real*8 x(method%prm%n)
    real*8 y(method%prm%n)

    integer MatFormat 
    MatFormat = method%prm%MatFormat 
    select case (MatFormat)
    case(MATRIX_FORMAT_CRS)
       call calc_Matvec_CRS(method,x,y)
    case(MATRIX_FORMAT_DIA)
       call calc_Matvec_DIA(method,x,y)
    case(MATRIX_FORMAT_DDM)
       call calc_Matvec_DDM(method%ddmMat,x,y)
    end select

  end subroutine calc_Matvec

  subroutine calc_res_CRS(method,x,b,r)
    use krylov
    use comm_kernel
    use matrix_kernel
    use vec_kernel

    implicit none
    type(KSP) :: method
    real*8 x(method%prm%n)
    real*8 b(method%prm%n)
    real*8 r(method%prm%n)
    real*8 alpha,beta
    call calc_Matvec_CRS(method,x,r)
    alpha= 1.0d0
    beta =-1.0d0
    call axpy(method%prm%n,alpha,beta,b,r)

  end subroutine calc_res_CRS

  ! r = b-Ax
  subroutine calc_res_DIA(method,x,b,r)
    use krylov
    use comm_kernel
    use matrix_kernel
    use vec_kernel

    implicit none
    type(KSP) :: method
    real*8 x(method%prm%n)
    real*8 b(method%prm%n)
    real*8 r(method%prm%n)
    real*8 alpha,beta

    call calc_Matvec_DIA(method,x,r)
    alpha= 1.0d0
    beta =-1.0d0
    call axpy(method%prm%n,alpha,beta,b,r)

  end subroutine calc_res_DIA

  ! y = Ax
  subroutine calc_Matvec_CRS(method,x,y)
    use krylov
    use comm_kernel
    use matrix_kernel
    use vec_kernel

    implicit none
    type(KSP) :: method
    real*8 x(method%prm%n)
    real*8 y(method%prm%n)
    integer iovlflag
    iovlflag=method%iovlflag
    select case (iovlflag)
    case(KSP_SPMV_OVL_NONE)
       call get_nrk_vec(method%commVec,x,method%buff)
       call get_nrk_vecLU(&
            method%commVec, method%buff, method%crsMat%nColL, &
            method%crsMat%nColU, method%rbuffL, method%rbuffU&
            )
       call MatVec_crs(&
            method%prm, method%crsMat, method%rbuffL, x,&
            method%rbuffU, y,&
            method%nBlock)
    case(KSP_SPMV_OVL_CORE)
       call MatVec_crs_ovl1(method%prm, method%crsMat, method%rbuffL,&
            x, method%rbuffU, y, method%commVec,&
            method%nBlock)
    case(KSP_SPMV_OVL_RANK)
       call MatVec_crs_ovl2(&
            method%prm, method%crsMat, x, y, method%commVec,&
            method%nBlock &
            )
    case default
       if(method%commVec%myrank.eq.0)then
          write(*,*)"error :: not support  ovl = ",iovlflag
       endif
    end select

  end subroutine calc_Matvec_CRS

  subroutine calc_Matvec_DIA(method,x,y)
    use krylov
    use comm_kernel
    use matrix_kernel
    use vec_kernel

    implicit none
    type(KSP) :: method
    real*8 x(method%prm%n)
    real*8 y(method%prm%n)
    integer(kind=intType) n
    integer iovlflag

    iovlflag=method%iovlflag
    n=method%prm%n

    select case (iovlflag)
    case(KSP_SPMV_OVL_NONE)
#if 0
       call get_nrk_vec(method%commVec,x,method%buff)
       call get_nrk_diaVec(method%commVec,method%buff,method%diaMat)
       call xcpy(n,x,&
            method%diaMat%Vec(1+method%diaMat%nColL))
       call MatVec_dia0(method%prm,method%diaMat,y)
#else
       call get_nrk_vec(method%commVec,x,method%buff)
       call MatVec_dia(method%prm,method%diaMat,x,y,method%commVec,method%buff,method%nBlock)
#endif
    case(KSP_SPMV_OVL_CORE)
       call MatVec_dia_ovl1(method%prm,method%diaMat,x,y,method%commVec,method%buff,method%nBlock)
    case(KSP_SPMV_OVL_RANK)
       call MatVec_dia_ovl2(method%prm,method%diaMat,x,y,method%commVec,method%nBlock)
    case default
       if(method%commVec%myrank.eq.0)then
          write(*,*)"error :: not support  ovl = ",iovlflag
       endif
    end select

  end subroutine calc_Matvec_DIA

  subroutine calc_precon(method,r,s)
    use krylov
    use precon_DILU,only : solve_DILU_DDM
    implicit none
    type(KSP) :: method
    real*8 s(method%prm%n)
    real*8 r(method%prm%n)
    integer MatFormat 

    MatFormat = method%prm_AS%MatFormat 
    select case (MatFormat)
    case(MATRIX_FORMAT_CRS)
       call precon_crs(method,r,s)
    case(MATRIX_FORMAT_DIA)
       call precon_dia(method,r,s)
    case(MATRIX_FORMAT_DDM)
       if(method%ipreflag.eq.21)then
          call solve_DILU_DDM(method,r,s)
       else
          call calc_precon_DDM(method%ddmMat,r,s)
       endif
    end select

  end subroutine calc_precon

  !  subroutine precon_crs(prm,crsMat,r,s,iflag,addL)
  subroutine precon_crs(method,r,s)
    use krylov
    use vec_kernel
    use matrix_kernel
    implicit none
    type(KSP) :: method
    real*8 s(method%prm%n)
    real*8 r(method%prm%n)
    integer iflag
    integer myrank

    myrank=method%commVec%myrank
    iflag=method%ipreflag

    ! ------
    select case (iflag)
    case(KSP_PRECONDITION_NONE)
       ! iflag = 0 :: none preconditioning
       call xcpy(method%prm%n,r,s)
       return
    case(KSP_PRECONDITION_POINTJACOBI)
       ! iflag = 1 :: point Jacobi
       call pointJacobi(method%prm%n,r,method%Dinv,s)
       return
    case(KSP_PRECONDITION_BLOCKJACOBI)
       ! iflag = 2 :: Block Jacobi, Additive Schwarz(addL=0)
!       call blockJacobi_crs(method%prm,method%crsMat,method%Dinv,r,s)
       if(method%ILU_method.eq.KSP_PRECONDITION_ILU)then
          ! ILU(0)
          call blockJacobi_ILU_crs(method,r,s)
       else
          ! D-ILU(0),M-ILU(0)
          call blockJacobi_crs(method,r,s)
       endif
       return
    case(KSP_PRECONDITION_ADDITIVESCHWARZ)
       ! iflag = 3 :: Additive Schwarz (addL>0)
       call AdditiveSchwarz_crs(method,r,s)
       return
    case(KSP_PRECONDITION_NEUMANN)
       ! iflag = 4 :: Neumann 
#if JACOBI_LIKE
       call Jacobilike_crs(method,r,s)
#else
       call Neumann_crs(method,r,s)
#endif
       return

#if 0
    case(4) ! 反復形式の加法シュワルツ (未確認)
       call zero_initialize(n,r1)
       call zero_initialize(n,As)
       call zero_initialize(n,sn)
       call blockJacobi_crs(method%prm,method%crsMat,method%Dinv,r,s)
       ! 反復ループ
       do i=1,4
          call calc_Matvec_CRS(method,s,As)
          call axpy3(n,-1.0d0,1.0d0,As,r,r1)
          call blockJacobi_crs(method%prm,method%crsMat,method%Dinv,r1,sn)
          call axpy(n,1.0d0,1.0d0,sn,s)
       enddo
       return
#endif
    end select
  end subroutine precon_crs

  subroutine precon_dia(method,r,s)
    use krylov
    use vec_kernel
    use matrix_kernel
    implicit none
    type(KSP) :: method
    real*8 s(method%prm%n)
    real*8 r(method%prm%n)
    integer iflag
    integer myrank

    myrank=method%commVec%myrank
    iflag=method%ipreflag

    ! ------
    select case (iflag)
    case(KSP_PRECONDITION_NONE)
       ! iflag = 0 :: none preconditioning
       call xcpy(method%prm%n,r,s)
       return
    case(KSP_PRECONDITION_POINTJACOBI)
       ! iflag = 1 :: point Jacobi
       call pointJacobi(method%prm%n,r,method%Dinv,s)
       return
    case(KSP_PRECONDITION_BLOCKJACOBI)
       ! iflag = 2 :: Block Jacobi, Additive Schwarz(addL=0)
       if(method%ILU_method.eq.KSP_PRECONDITION_ILU)then
          call blockJacobi_ILU_dia(method,r,s)
       else
          call blockJacobi_dia( & 
               method,r,s)
       endif
       return
    case(KSP_PRECONDITION_ADDITIVESCHWARZ)
       ! iflag = 3 :: Additive Schwarz (addL>0)
       call AdditiveSchwarz_dia(method,r,s)
       return
    case(KSP_PRECONDITION_NEUMANN)
       ! iflag = 4 :: Neumann 
       call Neumann_dia(method,r,s)
       return
    end select
  end subroutine precon_dia

! ----
  subroutine calc_Matvec_DDM(DDM_MAT_,x,y,PRECISION_CALC,PRECISION_A,x_lo,y_lo)
    use DDM_TYPE
    use ddm_matv,only : matv_ddm

    implicit none

    type (DDM_ALL_TYPE) DDM_MAT_
    double precision x(DDM_MAT_%n)
    double precision y(DDM_MAT_%n)
    integer ,optional :: PRECISION_CALC
    integer ,optional :: PRECISION_A
    double precision ,optional ::  x_lo(DDM_MAT_%n)
    double precision ,optional ::  y_lo(DDM_MAT_%n)

    if(  present(PRECISION_CALC)  )then
       call matv_ddm( &
            DDM_MAT_%icomm, &
            DDM_MAT_%mat%nnd_dia, &
            DDM_MAT_%n, &
            DDM_MAT_%mat%val_dia, &
            DDM_MAT_%mat%ioff_dia, &
            x, &
            y,  &
            DDM_MAT_%myrank, &
            DDM_MAT_%ncpu,   &
            DDM_MAT_%prm,&
            DDM_MAT_%commHide,&
            PRECISION_CALC,&
            PRECISION_A,&
            DDM_MAT_%mat%val_dia_lo, &
            x_lo,y_lo)
    else
       call matv_ddm( &
            DDM_MAT_%icomm, &
            DDM_MAT_%mat%nnd_dia, &
            DDM_MAT_%n, &
            DDM_MAT_%mat%val_dia, &
            DDM_MAT_%mat%ioff_dia, &
            x, &
            y,  &
            DDM_MAT_%myrank, &
            DDM_MAT_%ncpu,   &
            DDM_MAT_%prm,&
            DDM_MAT_%commHide )
    endif
  end subroutine calc_Matvec_DDM

  subroutine calc_res_DDM(DDM_MAT_,x,b,r)

    use DDM_TYPE
    use vec_kernel
    use PARCEL_integer
    implicit none

    type (DDM_ALL_TYPE) DDM_MAT_
    real*8 x(DDM_MAT_%n)
    real*8 b(DDM_MAT_%n)
    real*8 r(DDM_MAT_%n)
    real*8 alpha,beta
    integer(kind=intType) n_intType
    
    call calc_Matvec_DDM(DDM_MAT_,x,r)
    alpha= 1.0d0
    beta =-1.0d0
    n_intType=DDM_MAT_%n
    call axpy(n_intType,alpha,beta,b,r)

  end subroutine calc_res_DDM

  subroutine calc_precon_DDM(&
       DDM_MAT_,r,s,&
       PRECISION_CALC, &
       PRECISION_PRECONDITION, &
       r_lo, &
       s_lo  &
       )
    use DDM_TYPE
    use ddm_kernel,only : precon_main_ddm
    implicit none    
    type (DDM_ALL_TYPE) DDM_MAT_
    double precision r(DDM_MAT_%n)
    double precision s(DDM_MAT_%n)
    integer,optional :: PRECISION_CALC
    integer,optional :: PRECISION_PRECONDITION
    double precision,optional :: r_lo(DDM_MAT_%n)
    double precision,optional :: s_lo(DDM_MAT_%n)
    integer(kind=intType) parcel_int
    integer PRECISION
    PRECISION = KSP_DOUBLE_PRECISION
    if(  present(PRECISION_CALC)  )then
       PRECISION = PRECISION_CALC
    endif
    if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
       if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
          parcel_int = -1
          call precon_main_ddm(&
               DDM_MAT_%ipc, &
               DDM_MAT_%n, &
               parcel_int, &
               DDM_MAT_%nnd, &
               DDM_MAT_%diag, &
               DDM_MAT_%mat%val_dia,&
               DDM_MAT_%mat%ioff_dia, &
               r, &
               s, &
               DDM_MAT_%prm, &
               DDM_MAT_%prm%ncomp_grids, &
               DDM_MAT_%mat%ioff_grids, & 
               DDM_MAT_%k, &
               DDM_MAT_%icomm, &
               DDM_MAT_%myrank, &
               DDM_MAT_%AS_n, &
               DDM_MAT_%AS_prm, &
               DDM_MAT_%AS_mean, &
               DDM_MAT_%AS_g, &
               DDM_MAT_%AS_val_dia, &
               DDM_MAT_%AS_ioff_dia, &
               DDM_MAT_%AS_v1, &
               DDM_MAT_%AS_v2 , &
               DDM_MAT_%precon,&
               PRECISION_CALC, &
               PRECISION_PRECONDITION, &
               DDM_MAT_%diag_lo, &
               DDM_MAT_%mat%val_dia_lo,&
               r_lo, &
               s_lo, &
               DDM_MAT_%AS_val_dia_lo, &
               DDM_MAT_%AS_v1_lo, &
               DDM_MAT_%AS_v2_lo  &
               )
       else
          parcel_int = -1
          call precon_main_ddm(&
               DDM_MAT_%ipc, &
               DDM_MAT_%n, &
               parcel_int, &
               DDM_MAT_%nnd, &
               DDM_MAT_%diag, &
               DDM_MAT_%mat%val_dia,&
               DDM_MAT_%mat%ioff_dia, &
               r, &
               s, &
               DDM_MAT_%prm, &
               DDM_MAT_%prm%ncomp_grids, &
               DDM_MAT_%mat%ioff_grids, & 
               DDM_MAT_%k, &
               DDM_MAT_%icomm, &
               DDM_MAT_%myrank, &
               DDM_MAT_%AS_n, &
               DDM_MAT_%AS_prm, &
               DDM_MAT_%AS_mean, &
               DDM_MAT_%AS_g, &
               DDM_MAT_%AS_val_dia, &
               DDM_MAT_%AS_ioff_dia, &
               DDM_MAT_%AS_v1, &
               DDM_MAT_%AS_v2 , &
               DDM_MAT_%precon,&
               PRECISION_CALC, &
               PRECISION_PRECONDITION, &
!               DDM_MAT_%diag_lo, &
!               DDM_MAT_%mat%val_dia_lo,&
               v1_lo=r_lo, &
               v2_lo=s_lo, &
!               DDM_MAT_%AS_val_dia_lo, &
               AS_v1_lo=DDM_MAT_%AS_v1_lo,&
               AS_v2_lo=DDM_MAT_%AS_v2_lo & 
               )
       endif
    else
       parcel_int = -1
       call precon_main_ddm(&
         DDM_MAT_%ipc, &
         DDM_MAT_%n, &
         parcel_int, &
         DDM_MAT_%nnd, &
         DDM_MAT_%diag, &
         DDM_MAT_%mat%val_dia,&
         DDM_MAT_%mat%ioff_dia, &
         r, &
         s, &
         DDM_MAT_%prm, &
         DDM_MAT_%prm%ncomp_grids, &
         DDM_MAT_%mat%ioff_grids, & 
         DDM_MAT_%k, &
         DDM_MAT_%icomm, &
         DDM_MAT_%myrank, &
         DDM_MAT_%AS_n, &
         DDM_MAT_%AS_prm, &
         DDM_MAT_%AS_mean, &
         DDM_MAT_%AS_g, &
         DDM_MAT_%AS_val_dia, &
         DDM_MAT_%AS_ioff_dia, &
         DDM_MAT_%AS_v1, &
         DDM_MAT_%AS_v2 , &
         DDM_MAT_%precon )
    endif
    
   end subroutine calc_precon_DDM

end module krylov_kernel



