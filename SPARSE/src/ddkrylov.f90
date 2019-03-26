! Todo 
! 名前の付け方が適切でない
! dd の位置を統一したほうが良い
module ddkrylov_kernel
  use PARCEL_define
  implicit none
contains
  subroutine START_DDKSP_SOLVER(method)
    use krylov
    implicit none
    type(KSP) :: method
    integer solver
    integer myrank
    solver=method%solver
    myrank=method%commVec%myrank
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
  end subroutine START_DDKSP_SOLVER

  subroutine ddKSP_DIA( &
       icomm,&
       vecx_hi,vecx_lo,&
       vecb_hi,vecb_lo,&
       n,gn,istart,&
       diaA_hi,diaA_lo,offset,nnd,&
       itrmax,rtolmax,abstolmax,&
       input_param,&
       niter,reshistory,iret)
    use krylov
    use krylov_kernel
    use mpi
    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    integer(kind=intType) iend ! ベクトルの終了位置
    integer(kind=intType) nnd ! 対角数
    real*8 vecx_hi(n) ! 解ベクトル
    real*8 vecx_lo(n) ! 解ベクトル
    real*8 vecb_hi(n) ! 右辺ベクトル
    real*8 vecb_lo(n) ! 右辺ベクトル
    real*8 diaA_hi(nnd*n)  ! 行列値
    real*8 diaA_lo(nnd*n)  ! 行列値
    integer(kind=intType) offset(nnd) ! 対角オフセット

    type(KSP_INPUT_PARAM) :: input_param ! KSPの設定パラメータ
    integer itrmax ! 最大反復回数
    real*8 rtolmax ! 相対残差 (収束条件)
    real*8 abstolmax  !  絶対残差 (収束条件)
    integer niter
    integer iret 
    real*8 reshistory(itrmax) ! 残差履歴
    type(KSP) :: method

    integer myrank
    integer ierr
    integer i
    real*8 stime,etime
    real*8 initialize_time
    real*8 solver_time
    real*8 finalize_time 

    call PARCEL_KSP_Default_Setting(method) 
    call SET_KSP_PARAM(method,input_param)

    call MPI_Barrier(icomm,ierr)
    stime=MPI_Wtime()

    call set_ddKSP_DIA( &
         icomm, &
         vecx_hi,vecx_lo, &
         vecb_hi,vecb_lo, &
         n,gn,istart, &
         diaA_hi,diaA_lo,offset,nnd, &
         itrmax,rtolmax,abstolmax, &
         method)           

    call MPI_Barrier(icomm,ierr)
    etime=MPI_Wtime()
    initialize_time = etime-stime
    call MPI_Barrier(icomm,ierr)
    stime=MPI_Wtime()

    myrank=method%commVec%myrank
    call START_DDKSP_SOLVER(method)
    call MPI_Barrier(icomm,ierr)
    etime=MPI_Wtime()
    solver_time=etime-stime

    call MPI_Barrier(icomm,ierr)
    stime=MPI_Wtime()
    ! -------
    do i=1,method%niter
       reshistory(i)=method%reshistory(i)
    enddo
    iret=method%iret
    niter=method%niter
    call ddget_vecx_KSP(n,vecx_hi,vecx_lo,method)
    ! -------
    call free_ddKSP(method)
    call MPI_Barrier(icomm,ierr)
    etime=MPI_Wtime()
    finalize_time=etime-stime

    if(myrank.eq.0)then
       write(*,*)"initialize time = ",initialize_time,"[s]"
       write(*,*)"solver time     = ",solver_time,"[s]"
       write(*,*)"finalize time   = ",finalize_time,"[s]"
    endif

  end subroutine ddKSP_DIA
  
  subroutine ddKSP_CRS( &
       icomm, &
       vecx_hi,vecx_lo, &
       vecb_hi,vecb_lo, &
       n, &
       gn, &
       nnz, &
       istart, &
       crsA_hi, &
       crsA_lo, &
       crsRow_ptr, &
       crsCol, &
       itrmax,rtolmax,abstolmax,&
       input_param,&
       niter, &
       reshistory, &
       iret)
    use mpi
    use krylov
    use krylov_kernel
    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    integer(kind=intType) iend ! ベクトルの終了位置
    real*8 vecx_hi(n) ! 解ベクトル
    real*8 vecx_lo(n) ! 解ベクトル
    real*8 vecb_hi(n) ! 右辺ベクトル
    real*8 vecb_lo(n) ! 右辺ベクトル
    real*8  crsA_hi(nnz) ! 行列値
    real*8  crsA_lo(nnz) ! 行列値
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

    real*8 stime,etime

    real*8 initialize_time
    real*8 solver_time
    real*8 finalize_time 

    call PARCEL_KSP_Default_Setting(method)
    call SET_KSP_PARAM(method,input_param)

    call MPI_Barrier(icomm,ierr)
    stime=MPI_Wtime()
    call set_ddKSP_CRS( &
         icomm, &
         vecx_hi,vecx_lo, &
         vecb_hi,vecb_lo, &
         n,gn,nnz, &
         istart, &
         crsA_hi,crsA_lo,crsRow_ptr,crsCol, &
         itrmax,rtolmax,abstolmax, &
         method)
    call MPI_Barrier(icomm,ierr)
    etime=MPI_Wtime()
    initialize_time = etime-stime

    call MPI_Barrier(icomm,ierr)
    stime=MPI_Wtime()
    myrank=method%commVec%myrank
    call START_DDKSP_SOLVER(method)

    call MPI_Barrier(icomm,ierr)
    etime=MPI_Wtime()
    solver_time=etime-stime

    call MPI_Barrier(icomm,ierr)
    stime=MPI_Wtime()
    ! -------
    do i=1,method%niter
       reshistory(i)=method%reshistory(i)
    enddo
    iret=method%iret
    niter=method%niter
    call ddget_vecx_KSP(n,vecx_hi,vecx_lo,method)
    ! -------

    call free_ddKSP(method)
    call MPI_Barrier(icomm,ierr)
    etime=MPI_Wtime()

    finalize_time=etime-stime
    if(myrank.eq.0)then
       write(*,*)"initialize time = ",initialize_time,"[s]"
       write(*,*)"solver time     = ",solver_time,"[s]"
       write(*,*)"finalize time   = ",finalize_time,"[s]"
    endif
    
  end subroutine ddKSP_CRS

  subroutine set_ddKSP_CRS( &
       icomm, &
       vecx_hi, vecx_lo, &
       vecb_hi, vecb_lo, &
       n,gn,nnz,istart, &
       crsA_hi,crsA_lo,crsRow_ptr,crsCol, &
       itrmax,rtolmax,abstolmax, &
       method)

    use krylov
    use ddmatrix_kernel
    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    integer(kind=intType) iend ! ベクトルの終了位置
    real*8 vecx_hi(n) ! 解ベクトル
    real*8 vecx_lo(n) ! 解ベクトル
    real*8 vecb_hi(n) ! 右辺ベクトル
    real*8 vecb_lo(n) ! 右辺ベクトル
    real*8  crsA_hi(nnz) ! 行列値
    real*8  crsA_lo(nnz) ! 行列値
    integer(kind=intType) crsRow_ptr(n+1) ! ポインタ 
    integer(kind=intType) crsCol(nnz) ! 列番号
    integer ipreflag ! 前処理方法　
    integer addL ! 加法シュワルツのオーバラップ領域の大きさ
    integer iflagAS ! 加法シュワルツのオーバラップ領域補間方法
    integer itrmax ! 最大反復回数
    real*8 rtolmax ! 相対残差 (収束条件)
    real*8 abstolmax  !  絶対残差 (収束条件)
    type(KSP) :: method
    integer MatFormat
    integer i
    ! ------
    ! CRS Format
    iend  =istart + n -1
    method%prm%MatFormat=MATRIX_FORMAT_CRS ! --> CRS Format
    method%itrmax        = itrmax
    method%rtolmax       = rtolmax
    method%abstolmax     = abstolmax
    allocate(method%reshistory(itrmax))
    !$OMP PARALLEL DO PRIVATE(i)
    do i=1,itrmax
       method%reshistory(i)=0.0d0
    enddo

    call ddinit_KSP_CRS( &
         icomm, &
         n, &
         gn, &
         nnz, &
         istart, &
         iend, &
         crsA_hi,crsA_lo, &
         crsRow_ptr, &
         crsCol, &
         method)
    call set_ddvecx_KSP(n,vecx_hi,vecx_lo,method)
    call set_ddvecb_KSP(n,vecb_hi,vecb_lo,method)

    ipreflag = method%ipreflag
    addL     = method%addL
    iflagAS  = method%iflagAS
    call get_maxthreads(method%max_numth)

#if 0
    method%prm_AS%MatFormat=MATRIX_FORMAT_CRS
    select case (ipreflag)
    case(KSP_PRECONDITION_NONE)

    case(KSP_PRECONDITION_POINTJACOBI)
       call ddget_DInv_crs(method%prm,method%crsMat,method%Dinv_hi,method%Dinv_lo)
    case(KSP_PRECONDITION_BLOCKJACOBI)
       call ddget_DILU_DInv_crs(method%prm,method%crsMat,&
            method%Dinv_hi,method%Dinv_lo,&
            method%ILU_method)
    case(KSP_PRECONDITION_ADDITIVESCHWARZ)
       call ddget_DILU_DInv_AS_crs(method,method%crsMat)
    case default
       if(method%commVec%myrank.eq.0)then
          write(*,*)"error :: not support preconditioning = ",ipreflag
       end if
    end select
#else
    call set_ddprecon_CRS(ipreflag,method)
#endif

  end subroutine set_ddKSP_CRS

  subroutine set_ddprecon_CRS(ipreflag,method)
    use krylov
    use ddmatrix_kernel
    implicit none
    integer ipreflag
    type(KSP) :: method

    method%prm_AS%MatFormat=MATRIX_FORMAT_CRS
    select case (ipreflag)
    case(KSP_PRECONDITION_NONE)

    case(KSP_PRECONDITION_POINTJACOBI)
       call ddget_DInv_crs(method%prm,method%crsMat,method%Dinv_hi,method%Dinv_lo)
    case(KSP_PRECONDITION_BLOCKJACOBI)
       if(method%ILU_method.eq.KSP_PRECONDITION_ILU)then
          call ddget_ILU_LU_crs(method,method%crsMat)
       else
          call ddget_DILU_DInv_crs(method%prm,method%crsMat,&
               method%Dinv_hi,method%Dinv_lo,&
               method%ILU_method)
       endif
    case(KSP_PRECONDITION_ADDITIVESCHWARZ)
          call ddget_DILU_DInv_AS_crs(method,method%crsMat)
    case default
       if(method%commVec%myrank.eq.0)then
          write(*,*)"error :: not support preconditioning = ",ipreflag
       end if
    end select

  end subroutine set_ddprecon_CRS

  subroutine set_ddprecon_DIA(ipreflag,method)
    use krylov
    use matrix_kernel,only:get_nnd_crs
    use ddmatrix_kernel
    implicit none
    integer ipreflag
    type(KSP) :: method
    integer(kind=intType) n,nnz,istart,nnd

    n      = method%prm%n
    nnz    = method%prm%nnz
    istart = method%prm%istart
    nnd    = method%prm%nnd
    method%prm_AS%MatFormat=MATRIX_FORMAT_DIA

    select case (ipreflag)
    case(KSP_PRECONDITION_NONE)

    case(KSP_PRECONDITION_POINTJACOBI)
       call ddget_DInv_dia(&
            method%prm,method%diaMat,method%Dinv_hi,method%Dinv_lo)
    case(KSP_PRECONDITION_BLOCKJACOBI)
       if(method%ILU_method.eq.KSP_PRECONDITION_ILU)then
          call ddget_ILU_LU_dia(method,method%diaMat)
       else
          call ddget_DILU_DInv_dia(method%prm,method%diaMat,&
               method%Dinv_hi,method%Dinv_lo,method%ILU_method)
       endif
    case(KSP_PRECONDITION_ADDITIVESCHWARZ)
       call init_ddcrsMat(method%crsMat,n,nnz)
       call dddia2crs(method%diaMat,method%crsMat,n,istart,nnd,nnz)
       call ddget_DILU_DInv_AS_crs(method,method%crsMat)
       ! nndなので列ポインタのみで計算可能
       call get_nnd_crs( &
            method%precon_crsMatAS, &
            method%prm_AS%n, &
            method%prm_AS%istart, &
            method%prm_AS%nnz, &
            method%prm_AS%nnd)
       call init_dddiaMat( &
            method%precon_diaMatAS,&
            method%prm_AS%n,method%prm_AS%nnd)
       call ddcrs2dia( &
            method%precon_crsMatAS, &
            method%precon_diaMatAS, &
            method%prm_AS%n, &
            method%prm_AS%istart, &
            method%prm_AS%nnd, &
            method%prm_AS%nnz)
       if(method%addL.ne.0)then
          ! addLが0の場合はブロックヤコビ前処理なので加法シュヴァルツ前処理用の行列は作成しない
          call ddfree_crsMat(method%precon_crsMatAS)
       endif
       call ddfree_crsMat(method%crsMat)
    case default
       if(method%commVec%myrank.eq.0)then
          write(*,*)"error :: not support preconditioning = ",ipreflag
       end if
    end select

  end subroutine set_ddprecon_DIA

  subroutine set_ddKSP_DIA( &
       icomm, &
       vecx_hi,vecx_lo, &
       vecb_hi,vecb_lo, &
       n,gn,istart, &
       diaA_hi,diaA_lo,offset,nnd, &
       itrmax,rtolmax,abstolmax, &
       method)
    use krylov
!    use matrix_kernel,only:get_nnd_crs
    use ddmatrix_kernel

    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    integer(kind=intType) iend ! ベクトルの終了位置
    integer(kind=intType) nnd

    real*8 vecx_hi(n) ! 解ベクトル
    real*8 vecx_lo(n) ! 解ベクトル
    real*8 vecb_hi(n) ! 右辺ベクトル
    real*8 vecb_lo(n) ! 右辺ベクトル
    real*8 diaA_hi(nnd*n)
    real*8 diaA_lo(nnd*n)
    integer(kind=intType) offset(nnd)

    integer ipreflag ! 前処理方法　
    integer addL ! 加法シュワルツのオーバラップ領域の大きさ
    integer iflagAS ! 加法シュワルツのオーバラップ領域補間方法
    integer itrmax ! 最大反復回数
    real*8 rtolmax ! 相対残差 (収束条件)
    real*8 abstolmax  !  絶対残差 (収束条件)
    type(KSP) :: method
    integer MatFormat
    integer i
    ! ------
    ! DIA Format
    iend  =istart + n -1
    method%prm%MatFormat = MATRIX_FORMAT_DIA ! --> DIA Format
    method%itrmax        = itrmax
    method%rtolmax       = rtolmax
    method%abstolmax     = abstolmax
    allocate(method%reshistory(itrmax))
    !$OMP PARALLEL DO PRIVATE(i)
    do i=1,itrmax
       method%reshistory(i)=0.0d0
    enddo
    call ddget_nnz_dia(&
         diaA_hi,diaA_lo,&
         offset,n,istart,nnd,nnz)
    call ddinit_KSP_DIA(&
         icomm,&
         n,gn,nnz,istart,iend,&
         diaA_hi,diaA_lo,offset,nnd,&
         method)
    call set_ddvecx_KSP(n,vecx_hi,vecx_lo,method)
    call set_ddvecb_KSP(n,vecb_hi,vecb_lo,method)

    ipreflag = method%ipreflag
    addL     = method%addL
    iflagAS  = method%iflagAS
    call get_maxthreads(method%max_numth)
    ! -----
#if 0
    method%prm_AS%MatFormat=MATRIX_FORMAT_DIA

    select case (ipreflag)
    case(KSP_PRECONDITION_NONE)

    case(KSP_PRECONDITION_POINTJACOBI)
       call ddget_DInv_dia(&
            method%prm,method%diaMat,method%Dinv_hi,method%Dinv_lo)
    case(KSP_PRECONDITION_BLOCKJACOBI)
       call ddget_DILU_DInv_dia(method%prm,method%diaMat,&
            method%Dinv_hi,method%Dinv_lo,method%ILU_method)
    case(KSP_PRECONDITION_ADDITIVESCHWARZ)
       call init_ddcrsMat(method%crsMat,n,nnz)
       call dddia2crs(method%diaMat,method%crsMat,n,istart,nnd,nnz)
       call ddget_DILU_DInv_AS_crs(method,method%crsMat)
       ! nndなので列ポインタのみで計算可能
       call get_nnd_crs( &
            method%precon_crsMatAS, &
            method%prm_AS%n, &
            method%prm_AS%istart, &
            method%prm_AS%nnz, &
            method%prm_AS%nnd)
       call init_dddiaMat( &
            method%precon_diaMatAS,&
            method%prm_AS%n,method%prm_AS%nnd)
       call ddcrs2dia( &
            method%precon_crsMatAS, &
            method%precon_diaMatAS, &
            method%prm_AS%n, &
            method%prm_AS%istart, &
            method%prm_AS%nnd, &
            method%prm_AS%nnz)
       if(method%addL.ne.0)then
          ! addLが0の場合はブロックヤコビ前処理なので加法シュヴァルツ前処理用の行列は作成しない
          call ddfree_crsMat(method%precon_crsMatAS)
       endif
       call ddfree_crsMat(method%crsMat)
    case default
       if(method%commVec%myrank.eq.0)then
          write(*,*)"error :: not support preconditioning = ",ipreflag
       end if
    end select
#else
    call set_ddprecon_DIA(ipreflag,method)
#endif

  end subroutine set_ddKSP_DIA

! resetは整備しているがサンプルを整備していないので動作は未確認
  subroutine reset_ddCRS_Mat(&
       ipreflag,crsA_hi,crsA_lo,method)
    use krylov
!    use matrix_kernel
    use ddmatrix_kernel
    use vec_kernel
    implicit none
    type(KSP) :: method
    integer ipreflag
    real*8  crsA_hi(method%prm%nnz) ! 行列値上位ビット
    real*8  crsA_lo(method%prm%nnz) ! 行列値下位ビット
    integer(kind=intType) nnz
    nnz=method%prm%nnz
    if(method%ipreflag.eq.KSP_PRECONDITION_ADDITIVESCHWARZ)then
       deallocate(method%ivec_AS_hi)
       deallocate(method%ivec_AS_lo)
       deallocate(method%ovec_AS_hi)
       deallocate(method%ovec_AS_lo)
       call ddfree_crsMat(method%precon_crsMatAS)
       call ddfree_precon_AS(method)
    endif
    call ddfree_precon_ILU_crs(method)
    method%ipreflag=ipreflag
    call xcpy(nnz,crsA_hi,method%crsMat%crsA_hi)
    call xcpy(nnz,crsA_lo,method%crsMat%crsA_lo)
    call set_ddprecon_CRS(method%ipreflag,method)

  end subroutine reset_ddCRS_Mat

! resetは整備しているがサンプルを整備していないので動作は未確認
  subroutine reset_DIA_Mat(&
       ipreflag,diaMat_hi,diaMat_lo,method)
    use krylov
    use ddmatrix_kernel
    use vec_kernel
    implicit none
    type(KSP) :: method
    integer ipreflag
    real*8 diaMat_hi(method%prm%nnd * method%prm%n)
    real*8 diaMat_lo(method%prm%nnd * method%prm%n)
    integer(kind=intType) offset(method%prm%nnd)
    integer(kind=intType) n,nnd
    integer(kind=intType) n_nnd

    n_nnd=n*nnd
    n=method%prm%n
    nnd=method%prm%nnd
    if(method%ipreflag.eq.KSP_PRECONDITION_ADDITIVESCHWARZ)then
       deallocate(method%ivec_AS_hi)
       deallocate(method%ivec_AS_lo)
       deallocate(method%ovec_AS_hi)
       deallocate(method%ovec_AS_lo)
       call ddfree_diaMat(method%precon_diaMatAS)
       call ddfree_precon_AS(method)
    endif
    call ddfree_precon_ILU_dia(method)
    method%ipreflag=ipreflag
    call xcpy(n_nnd,diaMat_hi,method%diaMat%diaA_hi)
    call xcpy(n_nnd,diaMat_lo,method%diaMat%diaA_lo)
    call set_ddprecon_DIA(method%ipreflag,method)

  end subroutine reset_DIA_Mat

  subroutine set_ddCRS_Mat(n,nnz,crsA_hi,crsA_lo,crsRow_ptr,crsCol, method)
    use krylov
    implicit none
    integer(kind=intType) n,nnz
    real*8  crsA_hi(nnz) ! 行列値
    real*8  crsA_lo(nnz) ! 行列値
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
       method%crsMat%crsA_hi(i)=crsA_hi(i)
       method%crsMat%crsA_lo(i)=crsA_lo(i)
       method%crsMat%crsCol(i)=crsCol(i)
    enddo
    !$OMP END PARALLEL DO 

  end subroutine set_ddCRS_Mat

  subroutine ddinit_KSP_Mat_CRS(crsA_hi,crsA_lo,crsRow_ptr,crsCol, method)
    use krylov
    use comm_kernel
    use matrix_kernel,only:&
         get_crsCol_nCol,get_crs_Row,init_MatVec_crs
    use ddmatrix_kernel
    implicit none
    type(KSP) :: method
    real*8  crsA_hi(method%prm%nnz) ! 行列値
    real*8  crsA_lo(method%prm%nnz) ! 行列値
    integer(kind=intType) crsRow_ptr(method%prm%n+1) ! ポインタ 
    integer(kind=intType) crsCol(method%prm%nnz) ! 列番号

    integer(kind=intType) i
    integer(kind=intType) n,gn,nnz,istart,iend
    integer icomm

    integer(kind=intType) nCol,nColL,nColU,nColD

    integer(kind=intType),allocatable :: Colsort(:)

!    real*8 ts,te
    n     =method%prm%n
    gn    =method%prm%gn
    nnz   =method%prm%nnz
    istart=method%prm%istart
    iend  =method%prm%iend
    icomm  =method%commVec%comm

    allocate( Colsort(nnz) )
    call init_ddcrsMat(method%crsMat,n,nnz)
    call set_ddCRS_Mat(n,nnz, crsA_hi,crsA_lo,crsRow_ptr,crsCol, method)
    call get_crsCol_nCol(method%prm, method%crsMat%crsCol, Colsort,&
         & nCol,nColL,nColU,nColD)
    method%crsMat%nColL = nColL
    method%crsMat%nColU = nColU

    call get_crs_Row(method%prm,method%crsMat,nCol,Colsort)
    call set_comm(method%commVec,icomm,method%prm,Colsort,nCol)

    allocate(method%crsMat%iLs(nnz))
    call init_MatVec_crs(method%prm,method%crsMat)
!    do i=1,method%commVec%nrk_rank
!       allocate(method%commVec%nrk_prm(i)%iLs_Ovl2(n))
! ----
!        allocate(method%commVec%nrk_prm(i)%iLs_Ovl2(nnz))
!        allocate(method%commVec%nrk_prm(i)%iLe_Ovl2(nnz))

!    enddo
!    call init_MatVec_crs_ovl2(method%prm,method%crsMat,method%commVec)

    deallocate( Colsort )

  end subroutine ddinit_KSP_Mat_CRS

  subroutine set_ddDIA_Mat(n,nnd,diaMat_hi,diaMat_lo,offset, method)
    use krylov
    use ddvec_kernel

    implicit none
    integer(kind=intType) n,nnd
    real*8 diaMat_hi(nnd*n)
    real*8 diaMat_lo(nnd*n)
    integer(kind=intType) offset(nnd)
    type(KSP) :: method
    integer(kind=intType) i

    call ddxcpy(n*nnd,diaMat_hi,diaMat_lo,method%diaMat%diaA_hi,method%diaMat%diaA_lo)
    !$OMP PARALLEL DO PRIVATE(i)  
    do i=1,nnd
       method%diaMat%offset(i)=offset(i)
    enddo
    !$OMP END PARALLEL DO 

  end subroutine set_ddDIA_Mat

  subroutine ddinit_KSP_Mat_DIA(diaA_hi,diaA_lo,offset,method)

    use krylov
    use comm_kernel
    use matrix_kernel,only:get_dia_Row,init_MatVec_dia
    use ddmatrix_kernel
    implicit none

    type(KSP) :: method
    integer(kind=intType),allocatable :: Colsort(:)
    integer(kind=intType) i
    integer(kind=intType) n,gn,nnz,istart,iend,nnd
    real*8 diaA_hi(method%prm%nnd * method%prm%n)
    real*8 diaA_lo(method%prm%nnd * method%prm%n)
    integer(kind=intType) offset(method%prm%nnd)
    integer icomm

    n     =method%prm%n
    gn    =method%prm%gn
    nnz   =method%prm%nnz
    istart=method%prm%istart
    iend  =method%prm%iend
    icomm  =method%commVec%comm
    nnd   =method%prm%nnd

    call init_dddiaMat(method%diaMat,n,nnd)
    call set_ddDIA_Mat(n,nnd,diaA_hi,diaA_lo,offset,method)
    call set_dddia_Row(method%prm,method%diaMat)
    deallocate(method%diaMat%Vec_hi)
    deallocate(method%diaMat%Vec_lo)
    allocate(method%diaMat%Vec_hi(method%diaMat%nCol))
    allocate(method%diaMat%Vec_lo(method%diaMat%nCol))
    allocate( Colsort(nnz) )
    call ddget_nnz_Colsort_dia( &
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

    allocate(method%diaMat%iLs(n*nnd))
    call init_MatVec_dia(method%prm,method%diaMat,method%commVec)
!    do i=1,method%commVec%nrk_rank
!        allocate(method%commVec%nrk_prm(i)%iLs_Ovl2( n*nnd ) )
!        allocate(method%commVec%nrk_prm(i)%iLe_Ovl2( n*nnd ) )
!    enddo
!    call init_MatVec_dia(method%prm,method%diaMat,method%commVec)
    deallocate( Colsort )

  end subroutine ddinit_KSP_Mat_DIA

  subroutine ddinit_KSP_CRS( &
       icomm, &
       n, &
       gn, &
       nnz, &
       istart, &
       iend, &
       crsA_hi,crsA_lo, &
       crsRow_ptr, &
       crsCol, &
       method)
    use krylov
    use comm_kernel
    use ddmatrix_kernel
    use ddvec_kernel
    use vec_kernel,only:zero_initialize


    implicit none
    integer icomm
    integer(kind=intType) n,gn,nnz,istart,iend
    real*8  crsA_hi(nnz) ! 行列値
    real*8  crsA_lo(nnz) ! 行列値
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

    allocate(method%r_hi(n))
    allocate(method%r_lo(n))
    allocate(method%x_hi(n))
    allocate(method%x_lo(n))
    allocate(method%b_hi(n))
    allocate(method%b_lo(n))
    allocate(method%Dinv_hi(n))
    allocate(method%Dinv_lo(n))

    ! do i=1,n
    !    method%r(i)=0.0d0
    !    method%x(i)=0.0d0
    !    method%b(i)=0.0d0     
    ! enddo
    call zero_initialize(n,method%r_hi)
    call zero_initialize(n,method%r_lo)
    call zero_initialize(n,method%x_hi)
    call zero_initialize(n,method%x_lo)
    call zero_initialize(n,method%b_hi)
    call zero_initialize(n,method%b_lo)
    call zero_initialize(n,method%Dinv_hi)
    call zero_initialize(n,method%Dinv_lo)

    call ddinit_KSP_Mat_CRS(crsA_hi,crsA_lo,crsRow_ptr,crsCol, method)

    allocate(method%buff_hi(method%commVec%nrk_data))
    allocate(method%buff_lo(method%commVec%nrk_data))
    allocate(method%rbuffL_hi(method%crsMat%nColL))
    allocate(method%rbuffL_lo(method%crsMat%nColL))
    allocate(method%rbuffU_hi(method%crsMat%nColU))
    allocate(method%rbuffU_lo(method%crsMat%nColU))

    !$OMP PARALLEL DO PRIVATE(i)  
    do i=1,method%commVec%nrk_data
       method%buff_hi(i)=0.0d0
       method%buff_lo(i)=0.0d0
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i)  
    do i=1,method%crsMat%nColL
       method%rbuffL_hi(i)=0.0d0
       method%rbuffL_lo(i)=0.0d0
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i)  
    do i=1,method%crsMat%nColU
       method%rbuffU_hi(i)=0.0d0
       method%rbuffU_lo(i)=0.0d0
    enddo
    !$OMP END PARALLEL DO

  end subroutine ddinit_KSP_CRS

  subroutine ddinit_KSP_DIA( &
       icomm, &
       n,gn,nnz,istart,iend, &
       diaA_hi, diaA_lo,offset,nnd,&
        method)
    use krylov
    use comm_kernel
    use matrix_kernel
    use vec_kernel
    implicit none
    integer icomm
    integer(kind=intType) n,gn,nnz,istart,iend
    integer(kind=intType) nnd
    real*8 diaA_hi(nnd*n) ! 行列値
    real*8 diaA_lo(nnd*n) ! 行列値
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

    allocate(method%r_hi(n))
    allocate(method%r_lo(n))
    allocate(method%x_hi(n))
    allocate(method%x_lo(n))
    allocate(method%b_hi(n))
    allocate(method%b_lo(n))
    allocate(method%Dinv_hi(n))
    allocate(method%Dinv_lo(n))

    call zero_initialize(n,method%r_hi)
    call zero_initialize(n,method%r_lo)
    call zero_initialize(n,method%x_hi)
    call zero_initialize(n,method%x_lo)
    call zero_initialize(n,method%b_hi)
    call zero_initialize(n,method%b_lo)
    call zero_initialize(n,method%Dinv_hi)
    call zero_initialize(n,method%Dinv_lo)

    call ddinit_KSP_Mat_DIA(diaA_hi,diaA_lo,offset,method)

     !     method%commVec%nrk_data, &
     !     method%diaMat%nColL, &
     !     method%diaMat%nColU

    allocate(method%buff_hi(method%commVec%nrk_data))
    allocate(method%buff_lo(method%commVec%nrk_data))
    allocate(method%rbuffL_hi(method%diaMat%nColL))
    allocate(method%rbuffL_lo(method%diaMat%nColL))
    allocate(method%rbuffU_hi(method%diaMat%nColU))
    allocate(method%rbuffU_lo(method%diaMat%nColU))
       
    !$OMP PARALLEL DO PRIVATE(i)  
    do i=1,method%commVec%nrk_data
       method%buff_hi(i)=0.0d0
       method%buff_lo(i)=0.0d0
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i)  
    do i=1,method%diaMat%nColL
       method%rbuffL_hi(i)=0.0d0
       method%rbuffL_lo(i)=0.0d0
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i)  
    do i=1,method%diaMat%nColU
       method%rbuffU_hi(i)=0.0d0
       method%rbuffU_lo(i)=0.0d0
    enddo
    !$OMP END PARALLEL DO
  end subroutine ddinit_KSP_DIA

  subroutine free_ddKSP(method)
    use krylov
    implicit none
    type(KSP) :: method
    integer MatFormat 
    MatFormat = method%prm%MatFormat 
    select case (MatFormat)
    case(MATRIX_FORMAT_CRS)
       call free_ddKSP_CRS(method)
    case(MATRIX_FORMAT_DIA)
       call free_ddKSP_DIA(method)
    end select

  end subroutine free_ddKSP

  subroutine free_ddKSP_CRS(method)
    use krylov
    use comm_kernel
    use ddmatrix_kernel
    implicit none
    type(KSP) :: method
    deallocate(method%r_hi)
    deallocate(method%r_lo)
    deallocate(method%x_hi)
    deallocate(method%x_lo)
    deallocate(method%b_hi)
    deallocate(method%b_lo)
    deallocate(method%Dinv_hi)
    deallocate(method%Dinv_lo)

    deallocate(method%buff_hi)
    deallocate(method%buff_lo)
    deallocate(method%rbuffL_hi)
    deallocate(method%rbuffL_lo)
    deallocate(method%rbuffU_hi)
    deallocate(method%rbuffU_lo)

    call ddfree_precon_ILU_crs(method)
    !    call free_window(method%commVec)
    call free_comm(method%commVec)
    call ddfree_crsMat(method%crsMat)
    if(method%ipreflag.eq.KSP_PRECONDITION_ADDITIVESCHWARZ)then
       deallocate(method%ivec_AS_hi)
       deallocate(method%ivec_AS_lo)
       deallocate(method%ovec_AS_hi)
       deallocate(method%ovec_AS_lo)
       call ddfree_crsMat(method%precon_crsMatAS)
       call ddfree_precon_AS(method)
    endif
    deallocate(method%reshistory)

  end subroutine free_ddKSP_CRS

  subroutine free_ddKSP_DIA(method)
    use krylov
    use comm_kernel
    use ddmatrix_kernel
    implicit none
    type(KSP) :: method
    deallocate(method%r_hi)
    deallocate(method%r_lo)
    deallocate(method%x_hi)
    deallocate(method%x_lo)
    deallocate(method%b_hi)
    deallocate(method%b_lo)
    deallocate(method%Dinv_hi)
    deallocate(method%Dinv_lo)
    deallocate(method%buff_hi)
    deallocate(method%buff_lo)
    deallocate(method%rbuffL_hi)
    deallocate(method%rbuffL_lo)
    deallocate(method%rbuffU_hi)
    deallocate(method%rbuffU_lo)

    call ddfree_precon_ILU_dia(method)
!    call free_window(method%commVec)
    call free_comm(method%commVec)
    call ddfree_diaMat(method%diaMat)

    ! -----
    if(method%ipreflag.eq.KSP_PRECONDITION_ADDITIVESCHWARZ)then
       deallocate(method%ivec_AS_hi)
       deallocate(method%ivec_AS_lo)
       deallocate(method%ovec_AS_hi)
       deallocate(method%ovec_AS_lo)
       call ddfree_diaMat(method%precon_diaMatAS)
       call ddfree_precon_AS(method)
    endif
    deallocate(method%reshistory)

  end subroutine free_ddKSP_DIA
  
  subroutine set_ddvecb_KSP(n,b_hi,b_lo,method)
    use krylov
    use vec_kernel
    implicit none
    integer(kind=intType) n
    real*8 b_hi(n)
    real*8 b_lo(n)
    type(KSP) :: method

    integer(kind=intType) i
    call xcpy(n,b_hi,method%b_hi)
    call xcpy(n,b_lo,method%b_lo)

  end subroutine set_ddvecb_KSP

  subroutine set_ddvecx_KSP(n,x_hi,x_lo,method)
    use krylov
    use vec_kernel
    implicit none
    integer(kind=intType) n
    real*8 x_hi(n)
    real*8 x_lo(n)
    type(KSP) :: method

    integer(kind=intType) i
    ! do i=1,n
    !    method%x(i)=x(i)
    ! enddo
    call xcpy(n,x_hi,method%x_hi)
    call xcpy(n,x_lo,method%x_lo)

  end subroutine set_ddvecx_KSP

  subroutine ddget_vecx_KSP(n,x_hi,x_lo, method)
    use krylov
    use vec_kernel
    implicit none
    integer(kind=intType) n
    real*8 x_hi(n)
    real*8 x_lo(n)
    type(KSP) :: method

    integer(kind=intType) i
    ! do i=1,n
    !    x(i)=method%x(i)
    ! enddo
!    call ddxcpy(n,method%x_hi,method%x_lo,x_hi,x_lo)
    call xcpy(n,method%x_hi,x_hi)
    call xcpy(n,method%x_lo,x_lo)


  end subroutine ddget_vecx_KSP

  ! r = b-Ax
  subroutine ddcalc_res(method,x_hi,x_lo,b_hi,b_lo,r_hi,r_lo)
    use krylov
    implicit none
    type(KSP) :: method
    real*8 x_hi(method%prm%n)
    real*8 x_lo(method%prm%n)
    real*8 b_hi(method%prm%n)
    real*8 b_lo(method%prm%n)
    real*8 r_hi(method%prm%n)
    real*8 r_lo(method%prm%n)
    integer MatFormat 

    ! 入力データが倍精度でも計算は4倍精度で行う
    MatFormat = method%prm%MatFormat 
    select case (MatFormat)
    case(MATRIX_FORMAT_CRS)
       call ddcalc_res_CRS(method,&
            x_hi,x_lo,b_hi,b_lo,r_hi,r_lo)
    case(MATRIX_FORMAT_DIA)
       call ddcalc_res_DIA(method,&
            x_hi,x_lo,b_hi,b_lo,r_hi,r_lo)
    end select
  end subroutine ddcalc_res

  subroutine ddcalc_Matvec(method,x_hi,x_lo,y_hi,y_lo)
    use krylov
    implicit none
    type(KSP) :: method
    real*8 x_hi(method%prm%n)
    real*8 x_lo(method%prm%n)
    real*8 y_hi(method%prm%n)
    real*8 y_lo(method%prm%n)

    integer MatFormat 
    MatFormat = method%prm%MatFormat 
    select case (MatFormat)
    case(MATRIX_FORMAT_CRS)
       call ddcalc_Matvec_CRS(method,x_hi,x_lo,y_hi,y_lo)
    case(MATRIX_FORMAT_DIA)
       call ddcalc_Matvec_DIA(method,x_hi,x_lo,y_hi,y_lo)
    end select

  end subroutine ddcalc_Matvec

  subroutine ddcalc_res_CRS(method,x_hi,x_lo,b_hi,b_lo,r_hi,r_lo)
    use krylov
    use comm_kernel
    use matrix_kernel
    use ddvec_kernel

    implicit none
    type(KSP) :: method
    real*8 x_hi(method%prm%n)
    real*8 x_lo(method%prm%n)
    real*8 b_hi(method%prm%n)
    real*8 b_lo(method%prm%n)
    real*8 r_hi(method%prm%n)
    real*8 r_lo(method%prm%n)

    real*8 alpha_hi,alpha_lo
    real*8 beta_hi,beta_lo
    call ddcalc_Matvec_CRS(method,&
         x_hi,x_lo,r_hi,r_lo)
    alpha_hi= 1.0d0
    alpha_lo= 0.0d0
    beta_hi =-1.0d0
    beta_lo = 0.0d0
    call ddaxpy(method%prm%n,&
         alpha_hi,alpha_lo, &
         beta_hi,beta_lo, &
         b_hi,b_lo, &
         r_hi,r_lo)

  end subroutine ddcalc_res_CRS

  ! r = b-Ax
  subroutine ddcalc_res_DIA(method,&
       x_hi,x_lo,b_hi,b_lo,r_hi,r_lo)
    use krylov
    use ddcomm_kernel
    use ddmatrix_kernel
    use ddvec_kernel

    implicit none
    type(KSP) :: method
    real*8 x_hi(method%prm%n)
    real*8 x_lo(method%prm%n)
    real*8 b_hi(method%prm%n)
    real*8 b_lo(method%prm%n)
    real*8 r_hi(method%prm%n)
    real*8 r_lo(method%prm%n)
    real*8 alpha_hi,alpha_lo
    real*8 beta_hi,beta_lo
    
    call ddcalc_Matvec_DIA(method,x_hi,x_lo,r_hi,r_lo)
    alpha_hi= 1.0d0
    alpha_lo= 0.0d0
    beta_hi =-1.0d0
    beta_lo = 0.0d0
    call ddaxpy(method%prm%n,&
         alpha_hi,alpha_lo,&
         beta_hi,beta_lo,&
         b_hi,b_lo,&
         r_hi,r_lo)

  end subroutine ddcalc_res_DIA

  ! y = Ax
  subroutine ddcalc_Matvec_CRS(method,x_hi,x_lo,y_hi,y_lo)
    use krylov
    use ddcomm_kernel
    use ddmatrix_kernel
    use vec_kernel

    implicit none
    type(KSP) :: method
    real*8 x_hi(method%prm%n)
    real*8 x_lo(method%prm%n)
    real*8 y_hi(method%prm%n)
    real*8 y_lo(method%prm%n)
    integer iovlflag

    iovlflag=method%iovlflag
    select case (iovlflag)
    case(KSP_SPMV_OVL_NONE)
       call ddget_nrk_vec(&
            method%commVec,x_hi,x_lo,method%buff_hi,method%buff_lo)
       call ddget_nrk_vecLU(&
            method%commVec, &
            method%buff_hi, method%buff_lo, &
            method%crsMat%nColL, &
            method%crsMat%nColU, &
            method%rbuffL_hi,method%rbuffL_lo,&
            method%rbuffU_hi,method%rbuffU_lo &
       )
       call ddMatVec_crs(&
            method%prm, method%crsMat, &
            method%rbuffL_hi,method%rbuffL_lo, &
            x_hi            ,x_lo, &
            method%rbuffU_hi,method%rbuffU_lo, &
            y_hi            ,y_lo,&
            method%PRECISION_A)

    case(KSP_SPMV_OVL_CORE)
       call ddMatVec_crs_ovl1( &
            method%prm, method%crsMat, &
            method%rbuffL_hi,method%rbuffL_lo,&
            x_hi            ,x_lo, &
            method%rbuffU_hi,method%rbuffU_lo, &
            y_hi            ,y_lo, &
            method%commVec, &
            method%PRECISION_A)
    case(KSP_SPMV_OVL_RANK)
#if 0
       call ddMatVec_crs_ovl1( &
            method%prm, method%crsMat, &
            method%rbuffL_hi,method%rbuffL_lo,&
            x_hi            ,x_lo, &
            method%rbuffU_hi,method%rbuffU_lo, &
            y_hi            ,y_lo, &
            method%commVec, &
            method%PRECISION_A)
#else
        call ddMatVec_crs_ovl2(&
             method%prm, method%crsMat,&
             x_hi,x_lo, &
             y_hi,y_lo, &
             method%commVec,&
             method%PRECISION_A &
             )
#endif
    case default
       if(method%commVec%myrank.eq.0)then
          write(*,*)"error :: not support  ovl = ",iovlflag
       endif
    end select

  end subroutine ddcalc_Matvec_CRS

  subroutine ddcalc_Matvec_DIA(method,x_hi,x_lo,y_hi,y_lo)
    use krylov
    use ddcomm_kernel
    use ddmatrix_kernel
    use vec_kernel

    implicit none
    type(KSP) :: method
    real*8 x_hi(method%prm%n)
    real*8 x_lo(method%prm%n)
    real*8 y_hi(method%prm%n)
    real*8 y_lo(method%prm%n)
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
    call MatVec_dia(method%prm,method%diaMat,y)
#else 
    call ddget_nrk_vec( &
         method%commVec,x_hi,x_lo,method%buff_hi,method%buff_lo)
    call ddMatVec_dia(method%prm,method%diaMat,&
         x_hi,x_lo,&
         y_hi,y_lo,&
         method%commVec,&
         method%buff_hi,method%buff_lo,&
         method%PRECISION_A &
    )
#endif
    case(KSP_SPMV_OVL_CORE)
       call ddMatVec_dia_ovl1(method%prm, &
            method%diaMat, &
            x_hi,x_lo, &
            y_hi,y_lo, &
            method%commVec, &
            method%PRECISION_A)
    case(KSP_SPMV_OVL_RANK)
       call ddMatVec_dia_ovl2( &
            method%prm, method%diaMat, &
            x_hi,x_lo, &
            y_hi,y_lo, &
            method%commVec, &
            method%PRECISION_A)
    case default
       if(method%commVec%myrank.eq.0)then
          write(*,*)"error :: not support  ovl = ",iovlflag
       endif
    end select

  end subroutine ddcalc_Matvec_DIA

  subroutine ddcalc_precon(method,r_hi,r_lo,s_hi,s_lo)
    use krylov
    implicit none
    type(KSP) :: method
    real*8 s_hi(method%prm%n)
    real*8 s_lo(method%prm%n)
    real*8 r_hi(method%prm%n)
    real*8 r_lo(method%prm%n)
    integer MatFormat 

    MatFormat = method%prm_AS%MatFormat 
    select case (MatFormat)
    case(MATRIX_FORMAT_CRS)
       call ddprecon_crs(method,r_hi,r_lo,s_hi,s_lo)
    case(MATRIX_FORMAT_DIA)
       call ddprecon_dia(method,r_hi,r_lo,s_hi,s_lo)
    end select

  end subroutine ddcalc_precon

  !  subroutine precon_crs(prm,crsMat,r,s,iflag,addL)
  subroutine ddprecon_crs(method,r_hi,r_lo,s_hi,s_lo)
    use krylov
    use vec_kernel
    use ddmatrix_kernel
    implicit none
    type(KSP) :: method
    real*8 s_hi(method%prm%n)
    real*8 s_lo(method%prm%n)
    real*8 r_hi(method%prm%n)
    real*8 r_lo(method%prm%n)
    integer iflag
    integer myrank

    myrank=method%commVec%myrank
    iflag=method%ipreflag

    ! ------
    select case (iflag)
    case(KSP_PRECONDITION_NONE)
       ! iflag = 0 :: none preconditioning
       call xcpy(method%prm%n,r_hi,s_hi)
       call xcpy(method%prm%n,r_lo,s_lo)
       return
    case(KSP_PRECONDITION_POINTJACOBI)
       ! iflag = 1 :: point Jacobi
       call ddpointJacobi( &
            method%prm%n, &
            r_hi,r_lo, &
            method%Dinv_hi,method%Dinv_lo, &
            s_hi,s_lo,&
            method%PRECISION_PRECONDITION)
       return
    case(KSP_PRECONDITION_BLOCKJACOBI)
       ! iflag = 2 :: Block Jacobi, Additive Schwarz(addL=0)
       if(method%ILU_method.eq.KSP_PRECONDITION_ILU)then
          ! ILU(0)
          call ddblockJacobi_ILU_crs(method,r_hi,r_lo,s_hi,s_lo)
       else
       call ddblockJacobi_crs(&
            method%prm,method%crsMat, &
            method%Dinv_hi,method%Dinv_lo,&
            r_hi,r_lo,&
            s_hi,s_lo,&
            method%PRECISION_PRECONDITION)
    endif
       return
    case(KSP_PRECONDITION_ADDITIVESCHWARZ)
       ! iflag = 3 :: Additive Schwarz (addL>0)
       ! method%iflagAS = 0 :: BASIC 
       !   rの要素をオーバーラップさせて局所解の加算を解として採用する
       ! method%iflagAS = 1 :: RESTRICT
       !   rの要素をオーバーラップさせて局所解を解として採用する
       ! method%iflagAS = 2 :: INTERPOLATE
       !   rの要素をオーバーラップさせずに局所解の加算を解として採用する
       ! method%iflagAS = 3 :: NONE
       !   rの要素をオーバーラップさせずに局所解を採用する
       ! ---
       call ddAdditiveSchwarz_crs(method,r_hi,r_lo,s_hi,s_lo)
       return
    end select
  end subroutine ddprecon_crs

  subroutine ddprecon_dia(method,r_hi,r_lo,s_hi,s_lo)
    use krylov
    use vec_kernel
    use ddmatrix_kernel
    implicit none
    type(KSP) :: method
    real*8 s_hi(method%prm%n)
    real*8 s_lo(method%prm%n)
    real*8 r_hi(method%prm%n)
    real*8 r_lo(method%prm%n)
    integer iflag
    integer myrank

    myrank=method%commVec%myrank
    iflag=method%ipreflag

    ! ------
    select case (iflag)
    case(KSP_PRECONDITION_NONE)
       ! iflag = 0 :: none preconditioning
       call xcpy(method%prm%n,r_hi,s_hi)
       call xcpy(method%prm%n,r_lo,s_lo)
       return
    case(KSP_PRECONDITION_POINTJACOBI)
       ! iflag = 1 :: point Jacobi
       call ddpointJacobi( &
            method%prm%n, &
            r_hi,r_lo, &
            method%Dinv_hi,method%Dinv_lo, &
            s_hi,s_lo,&
            method%PRECISION_PRECONDITION)
       return
    case(KSP_PRECONDITION_BLOCKJACOBI)
       ! iflag = 2 :: Block Jacobi, Additive Schwarz(addL=0)
       if(method%ILU_method.eq.KSP_PRECONDITION_ILU)then
          ! ILU(0)
          call ddblockJacobi_ILU_dia(method,r_hi,r_lo,s_hi,s_lo)
       else
          call ddblockJacobi_dia(&
               method%prm,method%diaMat, &
               method%Dinv_hi,method%Dinv_lo,&
               r_hi,r_lo,&
               s_hi,s_lo,&
               method%PRECISION_PRECONDITION)
       endif

       return
    case(KSP_PRECONDITION_ADDITIVESCHWARZ)
       ! iflag = 3 :: Additive Schwarz (addL>0)
       ! method%iflagAS = 0 :: BASIC 
       !   rの要素をオーバーラップさせて局所解の加算を解として採用する
       ! method%iflagAS = 1 :: RESTRICT
       !   rの要素をオーバーラップさせて局所解を解として採用する
       ! method%iflagAS = 2 :: INTERPOLATE
       !   rの要素をオーバーラップさせずに局所解の加算を解として採用する
       ! method%iflagAS = 3 :: NONE
       !   rの要素をオーバーラップさせずに局所解を採用する
       ! ---
       !       call AdditiveSchwarz_crs(method,r,s)
       call ddAdditiveSchwarz_dia(method,r_hi,r_lo,s_hi,s_lo)
       return
    end select
  end subroutine ddprecon_dia

end module ddkrylov_kernel

