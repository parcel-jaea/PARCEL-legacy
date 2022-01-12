! module PARCEL_ddwrapper
! contains
  ! BiCG法(CRS)
  subroutine parcel_ddbicgstab(&
       icomm, &
       vecx_hi,vecx_lo,&
       vecb_hi,vecb_lo,&
       n,gn,nnz,istart,&
       crsA_hi,crsA_lo,crsRow_ptr,crsCol,&
       ipreflag,ILU_method,addL,iflagAS,&
       itrmax,rtolmax,&
       reshistory,&
       iovlflag,&
       precon_thblock, &
       independ_nvec, &
       nBlock, &
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION, &
       iret &
       )
    use PARCEL_define
    use krylov
!    use PARCEL_wrapper,only : INIT_PARCEL_INPUT
    implicit none
    ! -----
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    real*8 vecx_hi(n) ! 解ベクトル
    real*8 vecx_lo(n) ! 解ベクトル
    real*8 vecb_hi(n) ! 右辺ベクトル
    real*8 vecb_lo(n) ! 右辺ベクトル
    real*8  crsA_hi(nnz) ! 行列値
    real*8  crsA_lo(nnz) ! 行列値
    integer(kind=intType) crsRow_ptr(n+1) ! ポインタ 
    integer(kind=intType) crsCol(nnz) ! 列番号
    integer ipreflag ! 前処理方法　
    integer ILU_method ! D-ILU分解の方法
    integer addL ! 加法シュワルツのオーバラップ領域の大きさ
    integer iflagAS ! 加法シュワルツのオーバラップ領域補間方法
    integer itrmax ! 最大反復回数
    real*8 rtolmax ! 相対残差 (収束条件)
    real*8 abstolmax  !  絶対残差 (収束条件)
    integer solver ! 線形連立方程式の解法
    integer iovlflag ! 疎行列ベクトル積の通信・計算オーバーラップ方法
    integer iret
    real*8 reshistory(itrmax) ! 残差履歴
    integer niter
    integer precon_thblock
    integer independ_nvec
    integer nBlock

    integer PRECISION_A
    integer PRECISION_b
    integer PRECISION_x
    integer PRECISION_PRECONDITION
    type(KSP_INPUT_PARAM) :: PARCEL_INPUT
    abstolmax           = 0.0d0
    solver              = KSP_SOLVER_BiCG
    call INIT_PARCEL_INPUT(PARCEL_INPUT)
    PARCEL_INPUT%ipreflag=ipreflag
    PARCEL_INPUT%ILU_method=ILU_method
    PARCEL_INPUT%addL=addL
    PARCEL_INPUT%iflagAS=iflagAS
    PARCEL_INPUT%solver=solver
    PARCEL_INPUT%iovlflag=iovlflag
    PARCEL_INPUT%itrmax = itrmax
    PARCEL_INPUT%rtolmax = rtolmax
    PARCEL_INPUT%abstolmax = abstolmax
    PARCEL_INPUT%precon_thblock=precon_thblock
    PARCEL_INPUT%independ_nvec=independ_nvec
    PARCEL_INPUT%nBlock=nBlock

    PARCEL_INPUT%PRECISION_A=PRECISION_A
    PARCEL_INPUT%PRECISION_b=PRECISION_b
    PARCEL_INPUT%PRECISION_x=PRECISION_x
    PARCEL_INPUT%PRECISION_PRECONDITION=PRECISION_PRECONDITION
    call PARCEL_ddKSP_CRS( &
         icomm, &
         vecx_hi,vecx_lo,vecb_hi,vecb_lo, &
         n,gn,nnz,istart, &
         crsA_hi,crsA_lo,crsRow_ptr, crsCol,&
!         itrmax,rtolmax,abstolmax,&
         PARCEL_INPUT,&
         niter,reshistory,iret)
    itrmax=niter
  end subroutine parcel_ddbicgstab
  ! BiCG法(DIA)
  subroutine parcel_ddbicgstab_dia(&
       icomm, &
       vecx_hi,vecx_lo,&
       vecb_hi,vecb_lo,&
       n,gn,istart,&
       diaA_hi,diaA_lo,offset,nnd,&
       ipreflag,ILU_method,addL,iflagAS,&
       itrmax,rtolmax,&
       reshistory,&
       iovlflag,&
       precon_thblock, &
       independ_nvec, &
       nBlock, &
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION, &
       iret &
       )
    use PARCEL_define
    use krylov
!    use PARCEL_wrapper,only : INIT_PARCEL_INPUT
    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    integer(kind=intType) nnd ! 対角数
    real*8 vecx_hi(n) ! 解ベクトル
    real*8 vecx_lo(n) ! 解ベクトル
    real*8 vecb_hi(n) ! 右辺ベクトル
    real*8 vecb_lo(n) ! 右辺ベクトル
    real*8 diaA_hi(nnd*n)  ! 行列値
    real*8 diaA_lo(nnd*n)  ! 行列値
    integer(kind=intType) offset(nnd) ! 対角オフセット
    integer ipreflag ! 前処理方法　
    integer ILU_method ! D-ILU分解の方法
    integer addL ! 加法シュワルツのオーバラップ領域の大きさ
    integer iflagAS ! 加法シュワルツのオーバラップ領域補間方法
    integer itrmax ! 最大反復回数
    real*8 rtolmax ! 相対残差 (収束条件)
    real*8 abstolmax  !  絶対残差 (収束条件)
    integer solver ! 線形連立方程式の解法
    integer iovlflag ! 疎行列ベクトル積の通信・計算オーバーラップ方法
    real*8 reshistory(itrmax) ! 残差履歴
    integer niter
    integer iret
    integer precon_thblock
    integer independ_nvec
    integer nBlock

    integer PRECISION_A
    integer PRECISION_b
    integer PRECISION_x
    integer PRECISION_PRECONDITION
    type(KSP_INPUT_PARAM) :: PARCEL_INPUT
    abstolmax           = 0.0d0
    solver              = KSP_SOLVER_BiCG
    call INIT_PARCEL_INPUT(PARCEL_INPUT)
    PARCEL_INPUT%ipreflag=ipreflag
    PARCEL_INPUT%ILU_method=ILU_method
    PARCEL_INPUT%addL=addL
    PARCEL_INPUT%iflagAS=iflagAS
    PARCEL_INPUT%solver=solver
    PARCEL_INPUT%iovlflag=iovlflag
    PARCEL_INPUT%itrmax = itrmax
    PARCEL_INPUT%rtolmax = rtolmax
    PARCEL_INPUT%abstolmax = abstolmax
    PARCEL_INPUT%precon_thblock=precon_thblock
    PARCEL_INPUT%independ_nvec=independ_nvec
    PARCEL_INPUT%nBlock=nBlock

    PARCEL_INPUT%PRECISION_A=PRECISION_A
    PARCEL_INPUT%PRECISION_b=PRECISION_b
    PARCEL_INPUT%PRECISION_x=PRECISION_x
    PARCEL_INPUT%PRECISION_PRECONDITION=PRECISION_PRECONDITION
    call PARCEL_ddKSP_DIA( &
         icomm, &
         vecx_hi,vecx_lo,vecb_hi,vecb_lo, &
         n,gn,istart, &
         diaA_hi,diaA_lo,offset,nnd,&
!         itrmax,rtolmax,abstolmax, &
         PARCEL_INPUT,&
         niter,reshistory,iret)
    itrmax=niter
  end subroutine parcel_ddbicgstab_dia
  ! CG法(CRS)
  subroutine parcel_ddcg(&
       icomm, &
       vecx_hi,vecx_lo,&
       vecb_hi,vecb_lo,&
       n,gn,nnz,istart,&
       crsA_hi,crsA_lo,crsRow_ptr,crsCol,&
       ipreflag,ILU_method,addL,iflagAS,&
       itrmax,rtolmax,&
       reshistory,&
       iovlflag,&
       precon_thblock, &
       independ_nvec, &
       nBlock, &
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION, &
       iret &
       )
    use PARCEL_define
    use krylov
!    use PARCEL_wrapper,only : INIT_PARCEL_INPUT
    implicit none
    ! -----
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    real*8 vecx_hi(n) ! 解ベクトル
    real*8 vecx_lo(n) ! 解ベクトル
    real*8 vecb_hi(n) ! 右辺ベクトル
    real*8 vecb_lo(n) ! 右辺ベクトル
    real*8  crsA_hi(nnz) ! 行列値
    real*8  crsA_lo(nnz) ! 行列値
    integer(kind=intType) crsRow_ptr(n+1) ! ポインタ 
    integer(kind=intType) crsCol(nnz) ! 列番号
    integer ipreflag ! 前処理方法　
    integer ILU_method ! D-ILU分解の方法
    integer addL ! 加法シュワルツのオーバラップ領域の大きさ
    integer iflagAS ! 加法シュワルツのオーバラップ領域補間方法
    integer itrmax ! 最大反復回数
    real*8 rtolmax ! 相対残差 (収束条件)
    real*8 abstolmax  !  絶対残差 (収束条件)
    integer solver ! 線形連立方程式の解法
    integer iovlflag ! 疎行列ベクトル積の通信・計算オーバーラップ方法
    integer iret
    real*8 reshistory(itrmax) ! 残差履歴
    integer niter
    integer precon_thblock
    integer independ_nvec
    integer nBlock

    integer PRECISION_A
    integer PRECISION_b
    integer PRECISION_x
    integer PRECISION_PRECONDITION
    type(KSP_INPUT_PARAM) :: PARCEL_INPUT
!    write(*,*)"rtolmax = ",rtolmax
    abstolmax           = 0.0d0
    solver              = KSP_SOLVER_CG
    call INIT_PARCEL_INPUT(PARCEL_INPUT)
    PARCEL_INPUT%ipreflag=ipreflag
    PARCEL_INPUT%ILU_method=ILU_method
    PARCEL_INPUT%addL=addL
    PARCEL_INPUT%iflagAS=iflagAS
    PARCEL_INPUT%solver=solver
    PARCEL_INPUT%iovlflag=iovlflag
    PARCEL_INPUT%itrmax = itrmax
    PARCEL_INPUT%rtolmax = rtolmax
    PARCEL_INPUT%abstolmax = abstolmax
    PARCEL_INPUT%precon_thblock=precon_thblock
    PARCEL_INPUT%independ_nvec=independ_nvec
    PARCEL_INPUT%nBlock=nBlock

    PARCEL_INPUT%PRECISION_A=PRECISION_A
    PARCEL_INPUT%PRECISION_b=PRECISION_b
    PARCEL_INPUT%PRECISION_x=PRECISION_x
    PARCEL_INPUT%PRECISION_PRECONDITION=PRECISION_PRECONDITION
!    write(*,*)"PARCEL_INPUT%rtolmax -- ",PARCEL_INPUT%rtolmax
    call PARCEL_ddKSP_CRS( &
         icomm, &
         vecx_hi,vecx_lo,&
         vecb_hi,vecb_lo, &
         n,gn,nnz,istart, &
         crsA_hi,crsA_lo,crsRow_ptr, crsCol,&
!         itrmax,rtolmax,abstolmax,&
         PARCEL_INPUT,&
         niter,reshistory,iret)
    itrmax=niter
  end subroutine parcel_ddcg
  ! CG法(DIA)
  subroutine parcel_ddcg_dia(&
       icomm, &
       vecx_hi,vecx_lo,&
       vecb_hi,vecb_lo,&
       n,gn,istart,&
       diaA_hi,diaA_lo,offset,nnd,&
       ipreflag,ILU_method,addL,iflagAS,&
       itrmax,rtolmax,&
       reshistory,&
       iovlflag,&
       precon_thblock, &
       independ_nvec, &
       nBlock, &
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION , &
       iret &
       )
    use PARCEL_define
    use krylov
!    use PARCEL_wrapper,only : INIT_PARCEL_INPUT
    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    integer(kind=intType) nnd ! 対角数
    real*8 vecx_hi(n) ! 解ベクトル
    real*8 vecx_lo(n) ! 解ベクトル
    real*8 vecb_hi(n) ! 右辺ベクトル
    real*8 vecb_lo(n) ! 右辺ベクトル
    real*8 diaA_hi(nnd*n)  ! 行列値
    real*8 diaA_lo(nnd*n)  ! 行列値
    integer(kind=intType) offset(nnd) ! 対角オフセット
    integer ipreflag ! 前処理方法　
    integer ILU_method ! D-ILU分解の方法
    integer addL ! 加法シュワルツのオーバラップ領域の大きさ
    integer iflagAS ! 加法シュワルツのオーバラップ領域補間方法
    integer itrmax ! 最大反復回数
    real*8 rtolmax ! 相対残差 (収束条件)
    real*8 abstolmax  !  絶対残差 (収束条件)
    integer solver ! 線形連立方程式の解法
    integer iovlflag ! 疎行列ベクトル積の通信・計算オーバーラップ方法
    real*8 reshistory(itrmax) ! 残差履歴
    integer niter
    integer iret
    integer precon_thblock
    integer independ_nvec
    integer nBlock

    integer PRECISION_A
    integer PRECISION_b
    integer PRECISION_x
    integer PRECISION_PRECONDITION
    type(KSP_INPUT_PARAM) :: PARCEL_INPUT
    abstolmax           = 0.0d0
    solver              = KSP_SOLVER_CG
    call INIT_PARCEL_INPUT(PARCEL_INPUT)
    PARCEL_INPUT%ipreflag=ipreflag
    PARCEL_INPUT%ILU_method=ILU_method
    PARCEL_INPUT%addL=addL
    PARCEL_INPUT%iflagAS=iflagAS
    PARCEL_INPUT%solver=solver
    PARCEL_INPUT%iovlflag=iovlflag
    PARCEL_INPUT%itrmax = itrmax
    PARCEL_INPUT%rtolmax = rtolmax
    PARCEL_INPUT%abstolmax = abstolmax
    PARCEL_INPUT%precon_thblock=precon_thblock
    PARCEL_INPUT%independ_nvec=independ_nvec
    PARCEL_INPUT%nBlock=nBlock

    PARCEL_INPUT%PRECISION_A=PRECISION_A
    PARCEL_INPUT%PRECISION_b=PRECISION_b
    PARCEL_INPUT%PRECISION_x=PRECISION_x
    PARCEL_INPUT%PRECISION_PRECONDITION=PRECISION_PRECONDITION
    call PARCEL_ddKSP_DIA( &
         icomm, &
         vecx_hi,vecx_lo,vecb_hi,vecb_lo, &
         n,gn,istart, &
         diaA_hi,diaA_lo,offset,nnd,&
!         itrmax,rtolmax,abstolmax, &
         PARCEL_INPUT,&
         niter,reshistory,iret)
    itrmax=niter
  end subroutine parcel_ddcg_dia
  ! GMRES法(CRS)
  subroutine parcel_ddgmres( &
       icomm, &
       vecx_hi,vecx_lo,&
       vecb_hi,vecb_lo,&
       n,gn,nnz,istart,&
       crsA_hi,crsA_lo,crsRow_ptr,crsCol,&
       ipreflag,ILU_method,addL,iflagAS,&
       itrmax,rtolmax,&
       reshistory,&
       iovlflag, &
       precon_thblock, &
       independ_nvec, &
       nBlock, &
       GMRES_m,GMRES_GSflag, &
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION, &
       iret &
       )
    use PARCEL_define
    use krylov
!    use PARCEL_wrapper,only : INIT_PARCEL_INPUT
    implicit none
    ! -----
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    real*8 vecx_hi(n) ! 解ベクトル
    real*8 vecx_lo(n) ! 解ベクトル
    real*8 vecb_hi(n) ! 右辺ベクトル
    real*8 vecb_lo(n) ! 右辺ベクトル
    real*8  crsA_hi(nnz) ! 行列値
    real*8  crsA_lo(nnz) ! 行列値
    integer(kind=intType) crsRow_ptr(n+1) ! ポインタ 
    integer(kind=intType) crsCol(nnz) ! 列番号
    integer ipreflag ! 前処理方法　
    integer ILU_method ! D-ILU分解の方法
    integer addL ! 加法シュワルツのオーバラップ領域の大きさ
    integer iflagAS ! 加法シュワルツのオーバラップ領域補間方法
    integer itrmax ! 最大反復回数
    real*8 rtolmax ! 相対残差 (収束条件)
    real*8 abstolmax  !  絶対残差 (収束条件)
    integer solver ! 線形連立方程式の解法
    integer iovlflag ! 疎行列ベクトル積の通信・計算オーバーラップ方法
    integer GMRES_m ! GMRES(m)のリスタート数  
    integer GMRES_GSflag ! GMRES(m)のグラムシュミット直交化方法
    integer iret
    real*8 reshistory(itrmax) ! 残差履歴
    integer niter
    integer precon_thblock
    integer independ_nvec
    integer nBlock

    integer PRECISION_A
    integer PRECISION_b
    integer PRECISION_x
    integer PRECISION_PRECONDITION
    type(KSP_INPUT_PARAM) :: PARCEL_INPUT
    abstolmax           = 0.0d0
    solver              = KSP_SOLVER_GMRES
    call INIT_PARCEL_INPUT(PARCEL_INPUT)
    PARCEL_INPUT%ipreflag=ipreflag
    PARCEL_INPUT%ILU_method=ILU_method
    PARCEL_INPUT%addL=addL
    PARCEL_INPUT%iflagAS=iflagAS
    PARCEL_INPUT%solver=solver
    PARCEL_INPUT%iovlflag=iovlflag
    PARCEL_INPUT%itrmax = itrmax
    PARCEL_INPUT%rtolmax = rtolmax
    PARCEL_INPUT%abstolmax = abstolmax
    PARCEL_INPUT%GMRES_m=GMRES_m
    !  PARCEL_INPUT%GMRES_m=10
    PARCEL_INPUT%GMRES_GSflag=GMRES_GSflag
    PARCEL_INPUT%precon_thblock=precon_thblock
    PARCEL_INPUT%independ_nvec=independ_nvec
    PARCEL_INPUT%nBlock=nBlock

    PARCEL_INPUT%PRECISION_A=PRECISION_A
    PARCEL_INPUT%PRECISION_b=PRECISION_b
    PARCEL_INPUT%PRECISION_x=PRECISION_x
    PARCEL_INPUT%PRECISION_PRECONDITION=PRECISION_PRECONDITION
    !  write(*,*)""
    call PARCEL_ddKSP_CRS( &
         icomm, &
         vecx_hi,vecx_lo,vecb_hi,vecb_lo, &
         n,gn,nnz,istart, &
         crsA_hi,crsA_lo,crsRow_ptr,crsCol,&
!         itrmax,rtolmax,abstolmax,&
         PARCEL_INPUT,&
         niter,reshistory,iret)
    itrmax=niter
  end subroutine parcel_ddgmres
  ! GMRES法(DIA)
  subroutine parcel_ddgmres_dia(&
       icomm, &
       vecx_hi,vecx_lo,&
       vecb_hi,vecb_lo,&
       n,gn,istart,&
       diaA_hi,diaA_lo,offset,nnd,&
       ipreflag,ILU_method,addL,iflagAS,&
       itrmax,rtolmax,&
       reshistory,&
       iovlflag, &
       precon_thblock, &
       independ_nvec, &
       nBlock, &
       GMRES_m,GMRES_GSflag, &
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION, &
       iret &
       )
    use PARCEL_define
    use krylov
!    use PARCEL_wrapper,only : INIT_PARCEL_INPUT
    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    integer(kind=intType) nnd ! 対角数
    real*8 vecx_hi(n) ! 解ベクトル
    real*8 vecx_lo(n) ! 解ベクトル
    real*8 vecb_hi(n) ! 右辺ベクトル
    real*8 vecb_lo(n) ! 右辺ベクトル
    real*8 diaA_hi(nnd*n)  ! 行列値
    real*8 diaA_lo(nnd*n)  ! 行列値
    integer(kind=intType) offset(nnd) ! 対角オフセット
    integer ipreflag ! 前処理方法　
    integer ILU_method ! D-ILU分解の方法
    integer addL ! 加法シュワルツのオーバラップ領域の大きさ
    integer iflagAS ! 加法シュワルツのオーバラップ領域補間方法
    integer itrmax ! 最大反復回数
    real*8 rtolmax ! 相対残差 (収束条件)
    real*8 abstolmax  !  絶対残差 (収束条件)
    integer solver ! 線形連立方程式の解法
    integer iovlflag ! 疎行列ベクトル積の通信・計算オーバーラップ方法
    integer GMRES_m ! GMRES(m)のリスタート数  
    integer GMRES_GSflag ! GMRES(m)のグラムシュミット直交化方法
    real*8 reshistory(itrmax) ! 残差履歴
    integer niter
    integer iret
    integer precon_thblock
    integer independ_nvec
    integer nBlock

    integer PRECISION_A
    integer PRECISION_b
    integer PRECISION_x
    integer PRECISION_PRECONDITION
    type(KSP_INPUT_PARAM) :: PARCEL_INPUT
    abstolmax           = 0.0d0
    solver              = KSP_SOLVER_GMRES
    call INIT_PARCEL_INPUT(PARCEL_INPUT)
    PARCEL_INPUT%ipreflag=ipreflag
    PARCEL_INPUT%ILU_method=ILU_method
    PARCEL_INPUT%addL=addL
    PARCEL_INPUT%iflagAS=iflagAS
    PARCEL_INPUT%solver=solver
    PARCEL_INPUT%iovlflag=iovlflag
    PARCEL_INPUT%itrmax = itrmax
    PARCEL_INPUT%rtolmax = rtolmax
    PARCEL_INPUT%abstolmax = abstolmax
    PARCEL_INPUT%GMRES_m=GMRES_m
    PARCEL_INPUT%GMRES_GSflag=GMRES_GSflag
    PARCEL_INPUT%precon_thblock=precon_thblock
    PARCEL_INPUT%independ_nvec=independ_nvec
    PARCEL_INPUT%nBlock=nBlock

    PARCEL_INPUT%PRECISION_A=PRECISION_A
    PARCEL_INPUT%PRECISION_b=PRECISION_b
    PARCEL_INPUT%PRECISION_x=PRECISION_x
    PARCEL_INPUT%PRECISION_PRECONDITION=PRECISION_PRECONDITION
    call PARCEL_ddKSP_DIA( &
         icomm, &
         vecx_hi,vecx_lo,vecb_hi,vecb_lo, &
         n,gn,istart, &
         diaA_hi,diaA_lo,offset,nnd,&
!         itrmax,rtolmax,abstolmax, &
         PARCEL_INPUT,&
         niter,reshistory,iret)
    itrmax=niter
  end subroutine parcel_ddgmres_dia
  ! -----------
  ! @@@@@@@@@@@
  ! -----------
  subroutine PARCEL_ddKSP_CRS( &
       icomm, &
       vecx_hi,vecx_lo,&
       vecb_hi,vecb_lo,&
       n,gn,nnz,istart, &
       crsA_hi,crsA_lo,crsRow_ptr, crsCol,&
       !     itrmax,rtolmax,abstolmax,&
       PARCEL_INPUT,&
       niter,reshistory,iret)
    use krylov
    use krylov_kernel
    use ddkrylov_kernel
    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    real*8 vecx_hi(n) ! 解ベクトル 上位ビット
    real*8 vecx_lo(n) ! 解ベクトル 下位ビット
    real*8 vecb_hi(n) ! 右辺ベクトル 上位ビット
    real*8 vecb_lo(n) ! 右辺ベクトル 下位ビット
    real*8  crsA_hi(nnz) ! 行列値 上位ビット
    real*8  crsA_lo(nnz) ! 行列値 下位ビット
    integer(kind=intType) crsRow_ptr(n+1) ! ポインタ 
    integer(kind=intType) crsCol(nnz) ! 列番号
    integer itrmax ! 最大反復回数
    real*8 rtolmax ! 相対残差 (収束条件)
    real*8 abstolmax  !  絶対残差 (収束条件)
    type(KSP_INPUT_PARAM) :: PARCEL_INPUT
    integer niter
    integer iret ! 残差履歴
    real*8 reshistory(PARCEL_INPUT%itrmax)
!    write(*,*)"PARCEL_INPUT%rtolmax = ",PARCEL_INPUT%rtolmax
    call ddKSP_CRS(&
         icomm, &
         vecx_hi,vecx_lo,&
         vecb_hi,vecb_lo,&
         n,gn,nnz,&
         istart, &
         crsA_hi,crsA_lo,crsRow_ptr, crsCol,&
         !       itrmax,rtolmax,abstolmax,&
         PARCEL_INPUT,&
         niter,reshistory,iret)
  end subroutine PARCEL_ddKSP_CRS
  subroutine PARCEL_ddKSP_DIA( &
       icomm, &
       vecx_hi,vecx_lo,vecb_hi,vecb_lo, &
       n,gn,istart, &
       diaA_hi,diaA_lo,offset,nnd,&
       !     itrmax,rtolmax,abstolmax, &
       PARCEL_INPUT,&
       niter,reshistory,iret)
    use krylov
    use krylov_kernel
    use ddkrylov_kernel
    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    integer(kind=intType) nnd ! 対角数
    real*8 vecx_hi(n) ! 解ベクトル
    real*8 vecx_lo(n) ! 解ベクトル
    real*8 vecb_hi(n) ! 右辺ベクトル
    real*8 vecb_lo(n) ! 右辺ベクトル
    real*8 diaA_hi(nnd*n)  ! 行列値
    real*8 diaA_lo(nnd*n)  ! 行列値
    integer(kind=intType) offset(nnd) ! 対角オフセット
    integer itrmax ! 最大反復回数
    real*8 rtolmax ! 相対残差 (収束条件)
    real*8 abstolmax  !  絶対残差 (収束条件)
    type(KSP_INPUT_PARAM) :: PARCEL_INPUT
    real*8 reshistory(PARCEL_INPUT%itrmax) ! 残差履歴
    integer niter
    integer iret
    ! 未対応
    call ddKSP_DIA( &
         icomm, &
         vecx_hi,vecx_lo,&
         vecb_hi,vecb_lo, &
         n,gn,istart, &
         diaA_hi,diaA_lo,offset,nnd,&
!         itrmax,rtolmax,abstolmax,&
         PARCEL_INPUT,&
         niter,reshistory,iret)
  end subroutine PARCEL_ddKSP_DIA
! end module PARCEL_ddwrapper
