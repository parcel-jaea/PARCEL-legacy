module krylov_kernel
  use PARCEL_define
  implicit none
contains
  subroutine SET_KSP_PARAM(method,input_param)
    use krylov
    implicit none
    type(KSP) :: method
    type(KSP_INPUT_PARAM) :: input_param
    method%ipreflag               = input_param%ipreflag
    method%ILU_method             = input_param%ILU_method
    method%addL                   = input_param%addL
    method%iflagAS                = input_param%iflagAS
    method%solver                 = input_param%solver
    method%iovlflag               = input_param%iovlflag
    method%GMRES_m                = input_param%GMRES_m
    method%GMRES_GSflag           = input_param%GMRES_GSflag
    method%CAGMRES_sstep          = input_param%CAGMRES_sstep
    method%CAGMRES_tstep          = input_param%CAGMRES_tstep
    method%CAGMRES_basis          = input_param%CAGMRES_basis
    method%CAGMRES_QRflag         = input_param%CAGMRES_QRflag
    method%CBCG_kstep             = input_param%CBCG_kstep
    method%CBCG_Eigenflag         = input_param%CBCG_Eigenflag
    method%power_method_itrmax    = input_param%power_method_itrmax
    method%CAARNOLDI_sstep        = input_param%CAARNOLDI_sstep
    method%CAARNOLDI_tstep        = input_param%CAARNOLDI_tstep
    method%CAARNOLDI_basis        = input_param%CAARNOLDI_basis
    method%CAARNOLDI_QRflag       = input_param%CAARNOLDI_QRflag
    method%PRECISION_A            = input_param%PRECISION_A
    method%PRECISION_b            = input_param%PRECISION_b
    method%PRECISION_x            = input_param%PRECISION_x
    method%PRECISION_PRECONDITION = input_param%PRECISION_PRECONDITION
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
    write(*,*)"PRECISION_A           ",method%PRECISION_A           
    write(*,*)"PRECISION_b           ",method%PRECISION_b           
    write(*,*)"PRECISION_x           ",method%PRECISION_x           
    write(*,*)"PRECISION_PRECONDITION",method%PRECISION_PRECONDITION
    write(*,*)"----------------------"
    write(*,*)
  end subroutine view_KSP_PARAM

  subroutine START_KSP_SOLVER(method)
    use krylov
    implicit none
    type(KSP) :: method
    integer solver
    integer myrank
    solver=method%solver
    myrank=method%commVec%myrank

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
    case(KSP_CHECK_QR)
       call MPK_QR_TEST(method)
    case default
       if(myrank.eq.0)then
          write(*,*)"error :: not support solver = ",solver
       end if
    end select
  end subroutine START_KSP_SOLVER

  subroutine KSP_DIA( &
       icomm, &
       vecx,vecb, &
       n,gn,&
       istart, &
       diaA,offset,nnd,&
       itrmax,rtolmax,abstolmax, &
       input_param,&
       niter,reshistory,iret)
    use krylov
    use mpi
    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    integer(kind=intType) iend ! ベクトルの終了位置
    integer(kind=intType) nnd ! 対角数
    real*8 vecx(n) ! 解ベクトル
    real*8 vecb(n) ! 右辺ベクトル
    real*8 diaA(nnd*n)  ! 行列値
    integer(kind=intType) offset(nnd) ! 対角オフセット
    integer itrmax ! 最大反復回数
    real*8 rtolmax ! 相対残差 (収束条件)
    real*8 abstolmax  !  絶対残差 (収束条件)
    type(KSP_INPUT_PARAM) :: input_param ! KSPの設定パラメータ
    integer niter
    integer iret 
    real*8 reshistory(itrmax) ! 残差履歴

    type(KSP) :: method

    integer myrank
    integer ierr
    integer i
    integer ipreflag,ILU_method,addL
    real*8 stime,etime
    real*8 initialize_time
    real*8 solver_time
    real*8 finalize_time 

    call MPI_Barrier(icomm,ierr)
    stime=MPI_Wtime()

    call PARCEL_KSP_Default_Setting(method)
    call SET_KSP_PARAM(method,input_param)

    method%rtolmax       = rtolmax
    method%abstolmax     = abstolmax
    ipreflag             = method%ipreflag
    ILU_method           = method%ILU_method
    addL                 = method%addL

    call set_KSP_DIA(&
         icomm, &
         vecx,vecb,&
         n,gn,istart, &
         diaA,offset,nnd,&
!         itrmax,rtolmax,abstolmax,&
         itrmax,ipreflag,ILU_method,addL,&
         method)           
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
    call get_vecx_KSP(n,vecx,method)

    ! -------
    call free_KSP(method)
    call MPI_Barrier(icomm,ierr)
    etime=MPI_Wtime()
    finalize_time=etime-stime

    if(myrank.eq.0)then
       write(*,*)"initialize time = ",initialize_time,"[s]"
       write(*,*)"solver time     = ",solver_time,"[s]"
       write(*,*)"finalize time   = ",finalize_time,"[s]"
    endif

  end subroutine KSP_DIA
  
  subroutine KSP_CRS( &
       icomm, &
       vecx,vecb, &
       n,gn,nnz,&
       istart, &
       crsA,crsRow_ptr, crsCol,&
       itrmax,rtolmax,abstolmax,&
       input_param,&
       niter,reshistory,iret)
    use mpi
    use krylov
    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    integer(kind=intType) iend ! ベクトルの終了位置
    real*8 vecx(n) ! 解ベクトル
    real*8 vecb(n) ! 右辺ベクトル
    real*8  crsA(nnz) ! 行列値
    integer(kind=intType) crsRow_ptr(n+1) ! ポインタ 
    integer(kind=intType) crsCol(nnz) ! 列番号
    integer itrmax ! 最大反復回数
    real*8 rtolmax ! 相対残差 (収束条件)
    real*8 abstolmax  !  絶対残差 (収束条件)
    type(KSP_INPUT_PARAM) :: input_param ! KSPの設定パラメータ

    integer niter
    integer iret 
    real*8 reshistory(itrmax) ! 残差履歴

    type(KSP) :: method

    integer myrank
    integer ierr
    integer i
    integer ipreflag,ILU_method,addL

    real*8 stime,etime

    real*8 initialize_time
    real*8 solver_time
    real*8 finalize_time 

    call MPI_Barrier(icomm,ierr)
    stime=MPI_Wtime()
    call PARCEL_KSP_Default_Setting(method)
    call SET_KSP_PARAM(method,input_param)

    method%rtolmax       = rtolmax
    method%abstolmax     = abstolmax
    ipreflag             = method%ipreflag
    ILU_method           = method%ILU_method
    addL                 = method%addL
    call set_KSP_CRS(&
         icomm, &
         vecx,vecb, &
         n,gn,nnz,istart, &
         crsA,crsRow_ptr,crsCol, &
         itrmax,ipreflag,ILU_method,addL,&
         method)

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
    call get_vecx_KSP(n,vecx,method)
    ! -------

    call free_KSP(method)
    call MPI_Barrier(icomm,ierr)
    etime=MPI_Wtime()

    finalize_time=etime-stime
    if(myrank.eq.0)then
       write(*,*)"initialize time = ",initialize_time,"[s]"
       write(*,*)"solver time     = ",solver_time,"[s]"
       write(*,*)"finalize time   = ",finalize_time,"[s]"
    endif

  end subroutine KSP_CRS

  subroutine set_KSP_CRS(&
       icomm, &
       vecx,vecb, &
       n,gn,nnz,istart, &
       crsA,crsRow_ptr,crsCol, &
!       itrmax,rtolmax,abstolmax,&
       itrmax,ipreflag,ILU_method,addL, &
       method)

    use krylov
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
    real*8  crsA(nnz) ! 行列値
    integer(kind=intType) crsRow_ptr(n+1) ! ポインタ 
    integer(kind=intType) crsCol(nnz) ! 列番号
    integer itrmax ! 最大反復回数
    integer ipreflag ! 前処理方法　
    integer ILU_method ! 不完全LU分解(D-ILU分解)の対角行列の決定方法
    integer addL ! 加法シュワルツのオーバラップ領域の大きさ
    integer iflagAS ! 加法シュワルツのオーバラップ領域補間方法
    type(KSP) :: method
    integer i
    integer ierr
    real*8 stime,etime

    ! ------
    ! CRS Format
    iend  =istart + n -1
    method%prm%MatFormat = MATRIX_FORMAT_CRS 
    method%itrmax        = itrmax
    method%ipreflag      = ipreflag
    method%ILU_method    = ILU_method
    method%addL          = addL

!    method%rtolmax       = rtolmax
!    method%abstolmax     = abstolmax

    allocate(method%reshistory(itrmax))
    !$OMP PARALLEL DO PRIVATE(i)
    do i=1,itrmax
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
!     if(method%commVec%myrank.eq.0)then
!        write(*,*)"@ 1 @",etime-stime
!     endif

    stime=MPI_Wtime()
    call set_vecx_KSP(n,vecx,method)
    call set_vecb_KSP(n,vecb,method)
    call MPI_Barrier(icomm,ierr)
    etime=MPI_Wtime()
!     if(method%commVec%myrank.eq.0)then
!        write(*,*)"@ 2 @",etime-stime
!     endif

    method%ipreflag = ipreflag 
    method%addL     = addL
    iflagAS  = method%iflagAS
    call get_maxthreads(method%max_numth)

    call MPI_Barrier(icomm,ierr)
    stime=MPI_Wtime()
    call set_precon_CRS(ipreflag,method)
    call MPI_Barrier(icomm,ierr)
    etime=MPI_Wtime()
!     if(method%commVec%myrank.eq.0)then
!        write(*,*)"@ 3 @",etime-stime
!     endif

  end subroutine set_KSP_CRS
  
  subroutine set_KSP_DIA(&
       icomm, &
       vecx,vecb, &
       n,gn,istart, &
       diaA,offset,nnd,&
!       itrmax,rtolmax,abstolmax, &
       itrmax,ipreflag,ILU_method,addL,&
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
    integer itrmax ! 最大反復回数
    integer ipreflag ! 前処理方法　
    integer ILU_method ! 不完全LU分解(D-ILU分解)の対角行列の決定方法
    integer addL ! 加法シュワルツのオーバラップ領域の大きさ
    integer iflagAS ! 加法シュワルツのオーバラップ領域補間方法
    type(KSP) :: method
    integer i

    ! ------
    ! DIA Format
    iend  =istart + n -1
    method%prm%MatFormat = MATRIX_FORMAT_DIA 
    method%itrmax        = itrmax
    method%ipreflag      = ipreflag
    method%ILU_method    = ILU_method
    method%addL          = addL

!    method%rtolmax       = rtolmax
!    method%abstolmax     = abstolmax
    allocate(method%reshistory(itrmax))
    !$OMP PARALLEL DO PRIVATE(i)
    do i=1,itrmax
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

    method%ipreflag = ipreflag
    method%addL     = addL
    iflagAS  = method%iflagAS
    call get_maxthreads(method%max_numth)
    ! -----
#if 0
    method%prm_AS%MatFormat=MATRIX_FORMAT_DIA

    select case (ipreflag)
    case(KSP_PRECONDITION_NONE)

    case(KSP_PRECONDITION_POINTJACOBI)
       call get_DInv_dia(method%prm,method%diaMat,method%Dinv)
    case(KSP_PRECONDITION_BLOCKJACOBI)
       if(method%ILU_method.eq.KSP_PRECONDITION_ILU)then
          call init_crsMat(method%crsMat,n,nnz)
          call dia2crs(method%diaMat,method%crsMat,n,istart,nnd,nnz)
          call get_ILU_LU_crs(method,method%crsMat)
          call free_crsMat(method%crsMat)
       else
          call get_DILU_DInv_dia(method%prm,method%diaMat,method%Dinv,method%ILU_method)
       endif
    case(KSP_PRECONDITION_ADDITIVESCHWARZ)
       call init_crsMat(method%crsMat,n,nnz)
       call dia2crs(method%diaMat,method%crsMat,n,istart,nnd,nnz)
       call get_DILU_DInv_AS_crs(method,method%crsMat)
       call get_nnd_crs( &
            method%precon_crsMatAS, &
            method%prm_AS%n, &
            method%prm_AS%istart, &
            method%prm_AS%nnz, &
            method%prm_AS%nnd)
       call init_diaMat( &
            method%precon_diaMatAS,&
            method%prm_AS%n,method%prm_AS%nnd)
       call crs2dia( &
            method%precon_crsMatAS, &
            method%precon_diaMatAS, &
            method%prm_AS%n, &
            method%prm_AS%istart, &
            method%prm_AS%nnd, &
            method%prm_AS%nnz)
       if(method%addL.ne.0)then
          ! addLが0の場合はブロックヤコビ前処理なので加法シュヴァルツ前処理用の行列は作成しない
          call free_crsMat(method%precon_crsMatAS)
       endif
       call free_crsMat(method%crsMat)
    case default
       if(method%commVec%myrank.eq.0)then
          write(*,*)"error :: not support preconditioning = ",ipreflag
       end if
    end select
#else
    call set_precon_DIA(ipreflag,method)
#endif

  end subroutine set_KSP_DIA

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
    integer ILU_method
    integer addL
    real*8  crsA(method%prm%nnz) ! 行列値
    integer(kind=intType) nnz
    integer(kind=intType) i
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
       allocate(method%reshistory(itrmax))
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

    call xcpy(nnz,crsA,method%crsMat%crsA)
    call set_precon_CRS(method%ipreflag,method)

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
    integer ILU_method
    integer addL
    integer(kind=intType) n,nnd
    integer(kind=intType) i
    real*8 diaMat(method%prm%nnd * method%prm%n)

    n=method%prm%n
    nnd=method%prm%nnd
    if(method%ipreflag.eq.KSP_PRECONDITION_ADDITIVESCHWARZ)then
       deallocate(method%ivec_AS)
       deallocate(method%ovec_AS)
       call free_diaMat(method%precon_diaMatAS)
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
       allocate(method%reshistory(itrmax))
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

  end subroutine reset_DIA_Mat

  subroutine set_precon_DIA(ipreflag,method)
    use krylov
    use matrix_kernel
    implicit none
    integer ipreflag
    type(KSP) :: method
    integer(kind=intType) n,nnz,istart,nnd

    n      = method%prm%n
    istart = method%prm%istart
    nnd    = method%prm%nnd
    nnz    = method%prm%nnz

    method%prm_AS%MatFormat=MATRIX_FORMAT_DIA
!     write(*,*)ipreflag,method%addL,nnz
    select case (ipreflag)
    case(KSP_PRECONDITION_NONE)

    case(KSP_PRECONDITION_POINTJACOBI)
       call get_DInv_dia(method%prm,method%diaMat,method%Dinv)
    case(KSP_PRECONDITION_BLOCKJACOBI)
       if(method%ILU_method.eq.KSP_PRECONDITION_ILU)then
#if 0
          call init_crsMat(method%crsMat,n,nnz)
          call dia2crs(method%diaMat,method%crsMat,n,istart,nnd,nnz)
          call get_ILU_LU_crs(method,method%crsMat)
          call free_crsMat(method%crsMat)
!           call view_crs(n,istart,method%preconILU(1)%crsMatL)
!           call view_crs(n,istart,method%preconILU(1)%crsMatU)
!           stop
#else
          call get_ILU_LU_dia(method,method%diaMat)
!           call view_dia(&
!                n,&
!                method%preconILU(1)%nndL_th,&
!                istart,&
!                method%preconILU(1)%diaMatL)
!           call view_dia(&
!                n,&
!                method%preconILU(1)%nndU_th,&
!                istart,&
!                method%preconILU(1)%diaMatU)
!           stop
#endif

       else
          call get_DILU_DInv_dia(method%prm,method%diaMat,method%Dinv,method%ILU_method)
       endif
    case(KSP_PRECONDITION_ADDITIVESCHWARZ)
       call init_crsMat(method%crsMat,n,nnz)
       call dia2crs(method%diaMat,method%crsMat,n,istart,nnd,nnz)
       call get_DILU_DInv_AS_crs(method,method%crsMat)

       if(method%addL.ne.0)then
#if 0
          call get_nnd_crs( &
               method%precon_crsMatAS, &
               method%prm_AS%n, &
               method%prm_AS%istart, &
               method%prm_AS%nnz, &
               method%prm_AS%nnd)
          call init_diaMat( &
               method%precon_diaMatAS,&
               method%prm_AS%n,method%prm_AS%nnd)
          call crs2dia( &
               method%precon_crsMatAS, &
               method%precon_diaMatAS, &
               method%prm_AS%n, &
               method%prm_AS%istart, &
               method%prm_AS%nnd, &
               method%prm_AS%nnz)
#else
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
#endif
          ! addLが0の場合はブロックヤコビ前処理なので加法シュヴァルツ前処理用の行列は作成しない
          call free_crsMat(method%precon_crsMatAS)
       endif
       call free_crsMat(method%crsMat)

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
    use matrix_kernel
    implicit none
    integer ipreflag
    type(KSP) :: method

    method%prm_AS%MatFormat=MATRIX_FORMAT_CRS
    select case (ipreflag)
    case(KSP_PRECONDITION_NONE)

    case(KSP_PRECONDITION_POINTJACOBI)
       call get_DInv_crs(method%prm,method%crsMat,method%Dinv)
    case(KSP_PRECONDITION_BLOCKJACOBI)
       if(method%ILU_method.eq.KSP_PRECONDITION_ILU)then
          call get_ILU_LU_crs(method,method%crsMat)
       else
          call get_DILU_DInv_crs(method%prm,method%crsMat,method%Dinv,method%ILU_method)
       endif
    case(KSP_PRECONDITION_ADDITIVESCHWARZ)
       call get_DILU_DInv_AS_crs(method,method%crsMat)
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

  subroutine set_CRS_Mat(n,nnz, crsA,crsRow_ptr,crsCol, method)
    use krylov
    implicit none
    integer(kind=intType) n,nnz
    real*8  crsA(nnz) ! 行列値
    integer(kind=intType) crsRow_ptr(n+1) ! ポインタ 
    integer(kind=intType) crsCol(nnz) ! 列番号
    type(KSP) :: method
    integer(kind=intType) i

    !$OMP PARALLEL DO PRIVATE(i)  
    do i=1,n+1
       method%crsMat%crsRow_ptr(i)=crsRow_ptr(i)
    enddo
    !$OMP END PARALLEL DO 

    !$OMP PARALLEL DO PRIVATE(i)  
    do i=1,nnz
       method%crsMat%crsA(i)=crsA(i)
       method%crsMat%crsCol(i)=crsCol(i)
    enddo
    !$OMP END PARALLEL DO 

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
!      if(myrank.eq.0)then
!         write(*,*)"@ 1-1 @",te-ts
!      endif

    call MPI_Barrier(icomm,ierr)
    ts=MPI_Wtime()
    call set_CRS_Mat(n,nnz, crsA,crsRow_ptr,crsCol, method)
    call MPI_Barrier(icomm,ierr)
    te=MPI_Wtime()
!      if(myrank.eq.0)then
!         write(*,*)"@ 1-2 @",te-ts
!      endif

    call MPI_Barrier(icomm,ierr)
    ts=MPI_Wtime()
    call get_crsCol_nCol(method%prm, method%crsMat%crsCol, Colsort,&
         & nCol,nColL,nColU,nColD)
    call MPI_Barrier(icomm,ierr)
    te=MPI_Wtime()
!      if(myrank.eq.0)then
!         write(*,*)"@ 1-3 @",te-ts
!      endif

    method%crsMat%nColL = nColL
    method%crsMat%nColU = nColU

    call MPI_Barrier(icomm,ierr)
    ts=MPI_Wtime()
    call get_crs_Row(method%prm,method%crsMat,nCol,Colsort)
    call MPI_Barrier(icomm,ierr)
    te=MPI_Wtime()
!     if(myrank.eq.0)then
!        write(*,*)"@ 1-4 @",te-ts
!     endif

    call MPI_Barrier(icomm,ierr)
    ts=MPI_Wtime()
    call set_comm(method%commVec,icomm,method%prm,Colsort,nCol)
    call MPI_Barrier(icomm,ierr)
    te=MPI_Wtime()
!     if(myrank.eq.0)then
!        write(*,*)"@ 1-5 @",te-ts
!     endif

    ! ----
   !  call MPI_Barrier(icomm,ierr)
    ts=MPI_Wtime()
!    write(*,*)" @ nnz @", nnz
    allocate(method%crsMat%iLs(nnz))
    call init_MatVec_crs(method%prm,method%crsMat)
    if(method%iovlflag.eq.KSP_SPMV_OVL_RANK)then
       call init_MatVec_crs_ovl2(method%prm,method%crsMat,method%commVec)
    endif

!     do i=1,method%commVec%nrk_rank
!        allocate(method%commVec%nrk_prm(i)%iLs_Ovl2(nnz))
!        allocate(method%commVec%nrk_prm(i)%iLe_Ovl2(nnz))
!     enddo
    te=MPI_Wtime()
!      if(myrank.eq.0)then
!         write(*,*)"@ 1-6 @",te-ts
!      endif

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

    n      = method%prm%n
    gn     = method%prm%gn
    nnz    = method%prm%nnz
    istart = method%prm%istart
    iend   = method%prm%iend
    icomm  = method%commVec%comm
    nnd    = method%prm%nnd

    call init_diaMat(method%diaMat,n,nnd)
    call set_DIA_Mat(n,nnd, diaA,offset,method)
    call set_dia_Row(method%prm,method%diaMat)
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
    call set_comm(method%commVec, &
         icomm, &
         method%prm, &
         Colsort, &
         method%diaMat%nCol &
         )

    n_nnd=n*nnd
    allocate(method%diaMat%iLs( n_nnd ))
    call init_MatVec_dia(method%prm,method%diaMat,method%commVec)
    if(method%iovlflag.eq.KSP_SPMV_OVL_RANK)then
       call init_MatVec_dia_ovl2(method%prm,method%diaMat,method%commVec)
    endif
!    do i=1,method%commVec%nrk_rank
!       allocate(method%commVec%nrk_prm(i)%iLs_Ovl2( n*nnd ) )
!       allocate(method%commVec%nrk_prm(i)%iLe_Ovl2( n*nnd ) )
!    enddo
!    call init_MatVec_dia2(method%prm,method%diaMat,method%commVec)

    deallocate( Colsort )

  end subroutine init_KSP_Mat_DIA

  subroutine init_KSP_CRS(icomm, n,gn,nnz,istart,iend, crsA,crsRow_ptr,&
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
    implicit none
    type(KSP) :: method
    integer MatFormat 
    MatFormat = method%prm%MatFormat 
    select case (MatFormat)
    case(MATRIX_FORMAT_CRS)
       call free_KSP_CRS(method)
    case(MATRIX_FORMAT_DIA)
       call free_KSP_DIA(method)
    end select

  end subroutine free_KSP

  subroutine free_ILUregion(method)
    use krylov
    use common_kernel
    implicit none
    type(KSP) :: method
    integer i,max_numth
    ! -----
    call get_maxthreads(max_numth)
    if(allocated(method%ILUregionMat))then
       do i=1,max_numth
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
  end subroutine free_ILUregion
    ! -----

  subroutine free_KSP_CRS(method)
    use krylov
    use matrix_kernel
    use comm_kernel
    implicit none
    type(KSP) :: method

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

#if 0
    ! -----
    call get_maxthreads(max_numth)
! ここはResetで計算可能にする
    do i=1,max_numth
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
#else
    call  free_ILUregion(method)
#endif
    call free_precon_ILU_crs(method)
    ! -----
    if(method%ipreflag.eq.KSP_PRECONDITION_ADDITIVESCHWARZ)then
       deallocate(method%ivec_AS)
       deallocate(method%ovec_AS)
       call free_crsMat(method%precon_crsMatAS)
       call free_precon_AS(method)
    endif
    deallocate(method%reshistory)

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

#if 0
    ! -----
    call get_maxthreads(max_numth)
! ここはResetで計算可能にする
    do i=1,max_numth
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
    ! -----
#else
    call  free_ILUregion(method)
#endif
    call free_precon_ILU_dia(method)
    if(method%ipreflag.eq.KSP_PRECONDITION_ADDITIVESCHWARZ)then
       deallocate(method%ivec_AS)
       deallocate(method%ovec_AS)
       call free_diaMat(method%precon_diaMatAS)
       call free_precon_AS(method)
    endif
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
    use ddm_kernel
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
    use ddm_kernel
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
            method%rbuffU, y)
    case(KSP_SPMV_OVL_CORE)
       call MatVec_crs_ovl1(method%prm, method%crsMat, method%rbuffL,&
            x, method%rbuffU, y, method%commVec)
    case(KSP_SPMV_OVL_RANK)
       call MatVec_crs_ovl2(method%prm, method%crsMat, x, y, method%commVec &
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
       call MatVec_dia(method%prm,method%diaMat,x,y,method%commVec,method%buff)
#endif
    case(KSP_SPMV_OVL_CORE)
       call MatVec_dia_ovl1(method%prm,method%diaMat,x,y,method%commVec,method%buff)
    case(KSP_SPMV_OVL_RANK)
       call MatVec_dia_ovl2(method%prm,method%diaMat,x,y,method%commVec)
    case default
       if(method%commVec%myrank.eq.0)then
          write(*,*)"error :: not support  ovl = ",iovlflag
       endif
    end select

  end subroutine calc_Matvec_DIA

  subroutine calc_precon(method,r,s)
    use ddm_kernel
    use krylov
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
       call calc_precon_DDM(method%ddmMat,r,s)
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
    end select
  end subroutine precon_dia

! ----
  subroutine calc_Matvec_DDM(DDM_MAT,x,y)
    use DDM_TYPE
    implicit none

    type (DDM_ALL_TYPE) DDM_MAT
    double precision x(DDM_MAT%n)
    double precision y(DDM_MAT%n)

    write(*,*)"not support"
!     call matv_ddm( &
!          DDM_MAT%icomm, &
!          DDM_MAT%mat%nnd_dia, &
!          DDM_MAT%n, &
!          DDM_MAT%mat%val_dia, &
!          DDM_MAT%mat%ioff_dia, &
!          x, &
!          y,  &
!          DDM_MAT%myrank, &
!          DDM_MAT%ncpu,   &
!          DDM_MAT%mat%prm,&
!          DDM_MAT%commHide )

  end subroutine calc_Matvec_DDM

  subroutine calc_res_DDM(DDM_MAT,x,b,r)

    use DDM_TYPE
    use vec_kernel
    use PARCEL_integer
    implicit none

    type (DDM_ALL_TYPE) DDM_MAT
    real*8 x(DDM_MAT%n)
    real*8 b(DDM_MAT%n)
    real*8 r(DDM_MAT%n)
    real*8 alpha,beta
    integer(kind=intType) n_intType
    
    call calc_Matvec_DDM(DDM_MAT,x,r)
    alpha= 1.0d0
    beta =-1.0d0
    n_intType=DDM_MAT%n
    call axpy(n_intType,alpha,beta,b,r)

  end subroutine calc_res_DDM

  subroutine calc_precon_DDM(DDM_MAT,r,s)
    use DDM_TYPE
    implicit none
    type (DDM_ALL_TYPE) DDM_MAT
    double precision r(DDM_MAT%n)
    double precision s(DDM_MAT%n)

    write(*,*)"not support"

!     call precon_ddm(&
!          DDM_MAT%ipc, &
!          DDM_MAT%n, &
!          -1, &
!          DDM_MAT%nnd, &
!          DDM_MAT%diag, &
!          DDM_MAT%mat%val_dia,&
!          DDM_MAT%mat%ioff_dia, &
!          r, &
!          s, &
!          DDM_MAT%mat%prm, &
!          DDM_MAT%mat%prm%ncomp_grids, &
!          DDM_MAT%mat%ioff_grids, & 
!          DDM_MAT%k, &
!          DDM_MAT%icomm, &
!          DDM_MAT%myrank, &
!          DDM_MAT%AS_n, &
!          DDM_MAT%AS_prm, &
!          DDM_MAT%AS_diag_neighbor, &
!          DDM_MAT%AS_mean, &
!          DDM_MAT%AS_g, &
!          DDM_MAT%AS_val_dia, &
!          DDM_MAT%AS_ioff_dia, &
!          DDM_MAT%AS_v1, &
!          DDM_MAT%AS_v2 , &
!          DDM_MAT%precon )
   end subroutine calc_precon_DDM


end module krylov_kernel

