! module PARCEL_wrapper
! contains
  subroutine INIT_PARCEL_INPUT(PARCEL_INPUT)
    use krylov
    implicit none
    type(KSP_INPUT_PARAM) :: PARCEL_INPUT

    PARCEL_INPUT%ipreflag               = -1 
    PARCEL_INPUT%ILU_method             = -1 
    PARCEL_INPUT%addL                   = -1 
    PARCEL_INPUT%iflagAS                = -1 
    PARCEL_INPUT%solver                 = -1 
    PARCEL_INPUT%iovlflag               = -1 
    PARCEL_INPUT%GMRES_m                = -1 
    PARCEL_INPUT%GMRES_GSflag           = -1 
    PARCEL_INPUT%CAGMRES_sstep          = -1 
    PARCEL_INPUT%CAGMRES_tstep          = -1 
    PARCEL_INPUT%CAGMRES_basis          = -1 
    PARCEL_INPUT%CAGMRES_QRflag         = -1 
    PARCEL_INPUT%CBCG_kstep             = -1 
    PARCEL_INPUT%CBCG_Eigenflag         = -1 
    PARCEL_INPUT%power_method_itrmax    = -1 
    PARCEL_INPUT%CAARNOLDI_sstep        = -1 
    PARCEL_INPUT%CAARNOLDI_tstep        = -1 
    PARCEL_INPUT%CAARNOLDI_basis        = -1 
    PARCEL_INPUT%precon_thblock         = -1
    PARCEL_INPUT%independ_nvec          = -1
    PARCEL_INPUT%nBlock                 = -1
    PARCEL_INPUT%PRECISION_A            = -1
    PARCEL_INPUT%PRECISION_b            = -1
    PARCEL_INPUT%PRECISION_x            = -1
    PARCEL_INPUT%PRECISION_PRECONDITION = -1

  end subroutine INIT_PARCEL_INPUT
! end module PARCEL_wrapper

  subroutine  PARCEL_KSP_CRS( &
       icomm, &
       vecx,vecb, &
       n,gn,nnz,istart, &
       crsA,crsRow_ptr, crsCol,&
       !     itrmax,rtolmax,abstolmax,&
       PARCEL_INPUT,&
       niter,reshistory,iret, &
       EigenValueReal, &
       EigenValueImag, &
       EigenVector, &
       EigenErrorNorm, &
       MaxEigenID, &
       MinEigenID &
       )
    use krylov
    use krylov_kernel
    use PARCEL_integer
    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    real*8 vecx(n) ! 解ベクトル
    real*8 vecb(n) ! 右辺ベクトル
    real*8  crsA(nnz) ! 行列値
    integer(kind=intType) crsRow_ptr(n+1) ! ポインタ 
    integer(kind=intType) crsCol(nnz) ! 列番号
    integer itrmax ! 最大反復回数
    real*8 rtolmax ! 相対残差 (収束条件)
    real*8 abstolmax  !  絶対残差 (収束条件)
    type(KSP_INPUT_PARAM) :: PARCEL_INPUT
    integer niter
    integer iret ! 残差履歴
    real*8 reshistory(PARCEL_INPUT%itrmax)

    real*8,optional :: EigenValueReal(PARCEL_INPUT%CAARNOLDI_sstep*PARCEL_INPUT%CAARNOLDI_tstep)
    real*8,optional :: EigenValueImag(PARCEL_INPUT%CAARNOLDI_sstep*PARCEL_INPUT%CAARNOLDI_tstep)
    real*8,optional :: EigenVector(PARCEL_INPUT%CAARNOLDI_sstep*PARCEL_INPUT%CAARNOLDI_tstep)
    real*8,optional :: EigenErrorNorm(PARCEL_INPUT%CAARNOLDI_sstep*PARCEL_INPUT%CAARNOLDI_tstep)
    integer,optional :: MaxEigenID
    integer,optional :: MinEigenID
       
    if(  present(MinEigenID)  )then
       call KSP_CRS(&
            icomm, &
            vecx,vecb, &
            n,gn,nnz,&
            istart, &
            crsA,crsRow_ptr, crsCol,&
            PARCEL_INPUT,&
            niter,reshistory,iret,&
            EigenValueReal, &
            EigenValueImag, &
            EigenVector, &
            EigenErrorNorm, &
            MaxEigenID, &
            MinEigenID &
            )
    else
       call KSP_CRS(&
            icomm, &
            vecx,vecb, &
            n,gn,nnz,&
            istart, &
            crsA,crsRow_ptr, crsCol,&
            PARCEL_INPUT,&
            niter,reshistory,iret &
            )
    endif
  end subroutine PARCEL_KSP_CRS

  subroutine PARCEL_KSP_DIA( &
       icomm, &
       vecx,vecb, &
       n,gn,istart, &
       diaA,offset,nnd,&
       !     itrmax,rtolmax,abstolmax, &
       PARCEL_INPUT,&
       niter,reshistory,iret,&
       EigenValueReal, &
       EigenValueImag, &
       EigenVector, &
       EigenErrorNorm, &
       MaxEigenID, &
       MinEigenID &
       )
    use krylov
    use krylov_kernel
    use PARCEL_integer
    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    real*8 vecx(n) ! 解ベクトル
    real*8 vecb(n) ! 右辺ベクトル
    integer(kind=intType) nnd ! 対角数
    real*8 diaA(nnd*n)  ! 行列値
    integer(kind=intType) offset(nnd) ! 対角オフセット
    integer itrmax ! 最大反復回数
    real*8 rtolmax ! 相対残差 (収束条件)
    real*8 abstolmax  !  絶対残差 (収束条件)
    type(KSP_INPUT_PARAM) :: PARCEL_INPUT
    real*8 reshistory(PARCEL_INPUT%itrmax) ! 残差履歴
    integer niter
    integer iret

    real*8,optional :: EigenValueReal(PARCEL_INPUT%CAARNOLDI_sstep*PARCEL_INPUT%CAARNOLDI_tstep)
    real*8,optional :: EigenValueImag(PARCEL_INPUT%CAARNOLDI_sstep*PARCEL_INPUT%CAARNOLDI_tstep)
    real*8,optional :: EigenVector(PARCEL_INPUT%CAARNOLDI_sstep*PARCEL_INPUT%CAARNOLDI_tstep)
    real*8,optional :: EigenErrorNorm(PARCEL_INPUT%CAARNOLDI_sstep*PARCEL_INPUT%CAARNOLDI_tstep)
    integer,optional :: MaxEigenID
    integer,optional :: MinEigenID

    if(  present(MinEigenID)  )then
       call KSP_DIA( &
            icomm, &
            vecx,vecb, &
            n,gn,istart, &
            diaA,offset,nnd,&
            PARCEL_INPUT,&
            niter,reshistory,iret, &
            EigenValueReal, &
            EigenValueImag, &
            EigenVector, &
            EigenErrorNorm, &
            MaxEigenID, &
            MinEigenID &
            )
    else
       call KSP_DIA( &
            icomm, &
            vecx,vecb, &
            n,gn,istart, &
            diaA,offset,nnd,&
            PARCEL_INPUT,&
            niter,reshistory,iret &
            )
    endif

  end subroutine PARCEL_KSP_DIA


  ! BiCG法(CRS)
  subroutine parcel_dbicgstab(&
       icomm, &
       vecx,vecb,&
       n,gn,nnz,istart,&
       crsA,crsRow_ptr,crsCol,&
       ipreflag,ILU_method,addL,iflagAS,&
       itrmax,rtolmax,&
       reshistory,&
       iovlflag, &
       precon_thblock, &
       independ_nvec,  &
       nBlock, &
       iret &
       )
    use PARCEL_define
    use krylov
    use PARCEL_integer
    implicit none
    ! -----
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    real*8 vecx(n) ! 解ベクトル
    real*8 vecb(n) ! 右辺ベクトル
    real*8  crsA(nnz) ! 行列値
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
    integer precon_thblock
    integer independ_nvec
    integer nBlock
    integer iret
    real*8 reshistory(itrmax) ! 残差履歴
    integer niter
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
    PARCEL_INPUT%precon_thblock = precon_thblock
    PARCEL_INPUT%independ_nvec = independ_nvec
    PARCEL_INPUT%nBlock = nBlock 

    call PARCEL_KSP_CRS( &
         icomm, &
         vecx,vecb, &
         n,gn,nnz,istart, &
         crsA,crsRow_ptr, crsCol,&
!         itrmax,rtolmax,abstolmax,&
         PARCEL_INPUT,&
         niter,reshistory,iret)
    itrmax=niter
  end subroutine parcel_dbicgstab

  ! BiCG法(DIA)
  subroutine parcel_dbicgstab_dia(&
       icomm, &
       vecx,vecb,&
       n,gn,istart,&
       diaA,offset,nnd,&
       ipreflag,ILU_method,addL,iflagAS,&
       itrmax,rtolmax,&
       reshistory,&
       iovlflag,&
       precon_thblock, &
       independ_nvec,  &
       nBlock, &
       iret&
       )
    use PARCEL_define
    use krylov
    use PARCEL_integer
    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    real*8 vecx(n) ! 解ベクトル
    real*8 vecb(n) ! 右辺ベクトル
    integer(kind=intType) nnd ! 対角数
    real*8 diaA(nnd*n)  ! 行列値
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
    integer precon_thblock
    integer independ_nvec
    integer nBlock
    integer iret

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
    PARCEL_INPUT%precon_thblock = precon_thblock
    PARCEL_INPUT%independ_nvec = independ_nvec
    PARCEL_INPUT%nBlock = nBlock
    
    call PARCEL_KSP_DIA( &
         icomm, &
         vecx,vecb, &
         n,gn,istart, &
         diaA,offset,nnd,&
!         itrmax,rtolmax,abstolmax, &
         PARCEL_INPUT,&
         niter,reshistory,iret)
    itrmax=niter

  end subroutine parcel_dbicgstab_dia

  ! CG法(CRS)
  subroutine parcel_dcg(&
       icomm, &
       vecx,vecb,&
       n,gn,nnz,istart,&
       crsA,crsRow_ptr,crsCol,&
       ipreflag,ILU_method,addL,iflagAS,&
       itrmax,rtolmax,&
       reshistory,&
       iovlflag,&
       precon_thblock, &
       independ_nvec,  &
       nBlock, &
       iret&
       )
    use PARCEL_define
    use krylov
    use PARCEL_integer
    implicit none
    ! -----
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    real*8 vecx(n) ! 解ベクトル
    real*8 vecb(n) ! 右辺ベクトル
    real*8 crsA(nnz) ! 行列値
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
    integer precon_thblock
    integer independ_nvec
    integer nBlock
    integer iret
    real*8 reshistory(itrmax) ! 残差履歴
    integer niter

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
    PARCEL_INPUT%precon_thblock = precon_thblock
    PARCEL_INPUT%independ_nvec = independ_nvec
    PARCEL_INPUT%nBlock = nBlock 

    call PARCEL_KSP_CRS( &
         icomm, &
         vecx,vecb, &
         n,gn,nnz,istart, &
         crsA,crsRow_ptr, crsCol,&
!         itrmax,rtolmax,abstolmax,&
         PARCEL_INPUT,&
         niter,reshistory,iret)
    itrmax=niter
  end subroutine parcel_dcg

  ! CG法(DIA)
  subroutine parcel_dcg_dia(&
       icomm, &
       vecx,vecb,&
       n,gn,istart,&
       diaA,offset,nnd,&
       ipreflag,ILU_method,addL,iflagAS,&
       itrmax,rtolmax,&
       reshistory,&
       iovlflag,&
       precon_thblock, &
       independ_nvec,  &
       nBlock, &
       iret&
       )
    use PARCEL_define
    use krylov
    use PARCEL_integer
    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    real*8 vecx(n) ! 解ベクトル
    real*8 vecb(n) ! 右辺ベクトル
    integer(kind=intType) nnd ! 対角数
    real*8 diaA(nnd*n)  ! 行列値
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
    integer precon_thblock
    integer independ_nvec
    integer nBlock
    integer iret

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
    PARCEL_INPUT%precon_thblock = precon_thblock
    PARCEL_INPUT%independ_nvec = independ_nvec
    PARCEL_INPUT%nBlock = nBlock 

    call PARCEL_KSP_DIA( &
         icomm, &
         vecx,vecb, &
         n,gn,istart, &
         diaA,offset,nnd,&
!         itrmax,rtolmax,abstolmax, &
         PARCEL_INPUT,&
         niter,reshistory,iret)
    itrmax=niter

  end subroutine parcel_dcg_dia

  ! GMRES法(CRS)
  subroutine parcel_dgmres( &
       icomm, &
       vecx,vecb,&
       n,gn,nnz,istart,&
       crsA,crsRow_ptr,crsCol,&
       ipreflag,ILU_method,addL,iflagAS,&
       itrmax,rtolmax,&
       reshistory,&
       iovlflag, &
       precon_thblock, &
       independ_nvec,  &
       nBlock, &
       GMRES_m,GMRES_GSflag, &
       iret &
       )
    use PARCEL_define
    use krylov
    use PARCEL_integer
    implicit none
    ! -----
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    real*8 vecx(n) ! 解ベクトル
    real*8 vecb(n) ! 右辺ベクトル
    real*8 crsA(nnz) ! 行列値
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
    integer precon_thblock
    integer independ_nvec
    integer nBlock
    integer iret
    real*8 reshistory(itrmax) ! 残差履歴
    integer niter

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
    PARCEL_INPUT%precon_thblock = precon_thblock
    PARCEL_INPUT%independ_nvec = independ_nvec
    PARCEL_INPUT%nBlock = nBlock 
    PARCEL_INPUT%GMRES_m=GMRES_m
    PARCEL_INPUT%GMRES_GSflag=GMRES_GSflag

    call PARCEL_KSP_CRS( &
         icomm, &
         vecx,vecb, &
         n,gn,nnz,istart, &
         crsA,crsRow_ptr, crsCol,&
!         itrmax,rtolmax,abstolmax,&
         PARCEL_INPUT,&
         niter,reshistory,iret)
    itrmax=niter
  end subroutine parcel_dgmres

  ! GMRES法(DIA)
  subroutine parcel_dgmres_dia(&
       icomm, &
       vecx,vecb,&
       n,gn,istart,&
       diaA,offset,nnd,&
       ipreflag,ILU_method,addL,iflagAS,&
       itrmax,rtolmax,&
       reshistory,&
       iovlflag,&
       precon_thblock, &
       independ_nvec,  &
       nBlock, &
       GMRES_m,GMRES_GSflag, &
       iret &
       )
    use PARCEL_define
    use krylov
    use PARCEL_integer
    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    integer(kind=intType) nnd ! 対角数
    real*8 vecx(n) ! 解ベクトル
    real*8 vecb(n) ! 右辺ベクトル
    real*8 diaA(nnd*n)  ! 行列値
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
    integer precon_thblock
    integer independ_nvec
    integer nBlock
    integer iret

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
    PARCEL_INPUT%precon_thblock = precon_thblock
    PARCEL_INPUT%independ_nvec = independ_nvec
    PARCEL_INPUT%nBlock = nBlock 
    PARCEL_INPUT%GMRES_m=GMRES_m
    PARCEL_INPUT%GMRES_GSflag=GMRES_GSflag

    call PARCEL_KSP_DIA( &
         icomm, &
         vecx,vecb, &
         n,gn,istart, &
         diaA,offset,nnd,&
!         itrmax,rtolmax,abstolmax, &
         PARCEL_INPUT,&
         niter,reshistory,iret)
    itrmax=niter
  end subroutine parcel_dgmres_dia

  ! POWER法(CRS)
  subroutine parcel_dpower(&
       icomm, &
       vecx,vecb,&
       n,gn,nnz,istart,&
       crsA,crsRow_ptr,crsCol,&
       ipreflag,ILU_method,addL,iflagAS,&
       itrmax,&
       iovlflag, &
       precon_thblock, &
       independ_nvec,  &
       nBlock, &
       power_method_itrmax, &
       EigenValueReal, &
       EigenValueImag, &
       EigenVector, &
       EigenErrorNorm, &
       MaxEigenID, &
       MinEigenID, &
       iret &
       )
    use PARCEL_define
    use krylov
    use PARCEL_integer
    implicit none
    ! -----
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    real*8  crsA(nnz) ! 行列値
    integer(kind=intType) crsRow_ptr(n+1) ! ポインタ 
    integer(kind=intType) crsCol(nnz) ! 列番号
    real*8 vecx(n) ! 解ベクトル
    real*8 vecb(n) ! 右辺ベクトル
    integer ipreflag ! 前処理方法　
    integer ILU_method ! D-ILU分解の方法
    integer addL ! 加法シュワルツのオーバラップ領域の大きさ
    integer iflagAS ! 加法シュワルツのオーバラップ領域補間方法
    integer itrmax ! 最大反復回数
    real*8 reshistory(itrmax) ! 残差履歴
    real*8 rtolmax ! 相対残差 (収束条件)
    real*8 abstolmax  !  絶対残差 (収束条件)
    integer solver ! 線形連立方程式の解法
    integer iovlflag ! 疎行列ベクトル積の通信・計算オーバーラップ方法
    integer power_method_itrmax ! べき乗法の反復回数
    integer precon_thblock
    integer independ_nvec
    integer nBlock
    real*8,optional :: EigenValueReal(1)
    real*8,optional :: EigenValueImag(1)
    real*8,optional :: EigenVector(1)
    real*8,optional :: EigenErrorNorm(1)
    integer,optional :: MaxEigenID
    integer,optional :: MinEigenID
    integer iret

    integer niter

    type(KSP_INPUT_PARAM) :: PARCEL_INPUT
    abstolmax           = 0.0d0
    rtolmax             = 0.0d0
    solver              = KSP_SOLVER_POWER_METHOD
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
    PARCEL_INPUT%precon_thblock = precon_thblock
    PARCEL_INPUT%independ_nvec = independ_nvec
    PARCEL_INPUT%nBlock = nBlock 
    PARCEL_INPUT%power_method_itrmax=power_method_itrmax

    call PARCEL_KSP_CRS( &
         icomm, &
         vecx,vecb, &
         n,gn,nnz,istart, &
         crsA,crsRow_ptr, crsCol,&
!         itrmax,rtolmax,abstolmax,&
         PARCEL_INPUT,&
         niter,reshistory,iret,&
       EigenValueReal, &
       EigenValueImag, &
       EigenVector, &
       EigenErrorNorm, &
       MaxEigenID, &
       MinEigenID)

    itrmax=niter
  end subroutine parcel_dpower

  ! POWER法(DIA)
  subroutine parcel_dpower_dia(&
       icomm, &
       vecx,vecb,&
       n,gn,istart,&
       diaA,offset,nnd,&
       ipreflag,ILU_method,addL,iflagAS,&
       itrmax,&
       iovlflag, &
       precon_thblock, &
       independ_nvec,  &
       nBlock, &
       power_method_itrmax, &
       EigenValueReal, &
       EigenValueImag, &
       EigenVector, &
       EigenErrorNorm, &
       MaxEigenID, &
       MinEigenID, &
       iret &
       )
    use PARCEL_define
    use krylov
    use PARCEL_integer
    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    real*8 vecx(n) ! 解ベクトル
    real*8 vecb(n) ! 右辺ベクトル
    integer(kind=intType) nnd ! 対角数
    real*8 diaA(nnd*n)  ! 行列値
    integer(kind=intType) offset(nnd) ! 対角オフセット
    integer ipreflag ! 前処理方法　
    integer ILU_method ! D-ILU分解の方法
    integer addL ! 加法シュワルツのオーバラップ領域の大きさ
    integer iflagAS ! 加法シュワルツのオーバラップ領域補間方法
    integer itrmax ! 最大反復回数
    real*8 reshistory(itrmax) ! 残差履歴
    real*8 rtolmax ! 相対残差 (収束条件)
    real*8 abstolmax  !  絶対残差 (収束条件)
    integer solver ! 線形連立方程式の解法
    integer iovlflag ! 疎行列ベクトル積の通信・計算オーバーラップ方法
    integer power_method_itrmax ! べき乗法の反復回数

    integer niter
    integer precon_thblock
    integer independ_nvec
    integer nBlock
    real*8,optional :: EigenValueReal(1)
    real*8,optional :: EigenValueImag(1)
    real*8,optional :: EigenVector(1)
    real*8,optional :: EigenErrorNorm(1)
    integer,optional :: MaxEigenID
    integer,optional :: MinEigenID
    integer iret

    type(KSP_INPUT_PARAM) :: PARCEL_INPUT
    abstolmax           = 0.0d0
    rtolmax             = 0.0d0
    solver              = KSP_SOLVER_POWER_METHOD
    call INIT_PARCEL_INPUT(&
         PARCEL_INPUT&
         )
    PARCEL_INPUT%ipreflag=ipreflag
    PARCEL_INPUT%ILU_method=ILU_method
    PARCEL_INPUT%addL=addL
    PARCEL_INPUT%iflagAS=iflagAS
    PARCEL_INPUT%solver=solver
    PARCEL_INPUT%iovlflag=iovlflag
    PARCEL_INPUT%itrmax = itrmax
    PARCEL_INPUT%rtolmax = rtolmax
    PARCEL_INPUT%abstolmax = abstolmax
    PARCEL_INPUT%precon_thblock = precon_thblock
    PARCEL_INPUT%independ_nvec = independ_nvec
    PARCEL_INPUT%nBlock = nBlock 
    PARCEL_INPUT%power_method_itrmax=power_method_itrmax

    call PARCEL_KSP_DIA( &
         icomm, &
         vecx,vecb, &
         n,gn,istart, &
         diaA,offset,nnd,&
!         itrmax,rtolmax,abstolmax, &
         PARCEL_INPUT,&
         niter,reshistory,iret, &
       EigenValueReal, &
       EigenValueImag, &
       EigenVector, &
       EigenErrorNorm, &
       MaxEigenID, &
       MinEigenID )
    itrmax=niter

  end subroutine parcel_dpower_dia


  ! CBCG法(CRS)
  subroutine parcel_dcbcg(&
       icomm, &
       vecx,vecb,&
       n,gn,nnz,istart,&
       crsA,crsRow_ptr,crsCol,&
       ipreflag,ILU_method,addL,iflagAS,&
       itrmax,rtolmax,&
       reshistory,&
       iovlflag, &
       precon_thblock, &
       independ_nvec,  &
       nBlock, &
       CBCG_kstep,CBCG_Eigenflag,power_method_itrmax, &
       CAARNOLDI_sstep,CAARNOLDI_tstep, &
       CAARNOLDI_basis,CAARNOLDI_QRflag, &
       iret &
       )
    use PARCEL_define
    use krylov
    use PARCEL_integer
    implicit none
    ! -----
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    real*8  crsA(nnz) ! 行列値
    integer(kind=intType) crsRow_ptr(n+1) ! ポインタ 
    integer(kind=intType) crsCol(nnz) ! 列番号
    real*8 vecx(n) ! 解ベクトル
    real*8 vecb(n) ! 右辺ベクトル
    integer ipreflag ! 前処理方法　
    integer ILU_method ! D-ILU分解の方法
    integer addL ! 加法シュワルツのオーバラップ領域の大きさ
    integer iflagAS ! 加法シュワルツのオーバラップ領域補間方法
    integer itrmax ! 最大反復回数
    real*8 rtolmax ! 相対残差 (収束条件)
    real*8 abstolmax  !  絶対残差 (収束条件)
    integer solver ! 線形連立方程式の解法
    integer iovlflag ! 疎行列ベクトル積の通信・計算オーバーラップ方法
    integer CBCG_kstep ! CBCGのステップ数
    integer CBCG_Eigenflag ! CBCG法の固有値計算方法
    integer power_method_itrmax ! べき乗法の反復回数
    integer CAARNOLDI_sstep  ! CAARNOLDIのs step数 
    integer CAARNOLDI_tstep  ! CAARNOLDIのt step数 
    integer CAARNOLDI_basis  ! CAARNOLDIの基底生成方法
    integer CAARNOLDI_QRflag ! CAARNOLDIのQR分解方法
    integer precon_thblock
    integer independ_nvec
    integer nBlock
    integer iret
    real*8 reshistory(itrmax) ! 残差履歴
    integer niter

    type(KSP_INPUT_PARAM) :: PARCEL_INPUT
    abstolmax           = 0.0d0
    solver              = KSP_SOLVER_CBCG
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
    PARCEL_INPUT%precon_thblock = precon_thblock
    PARCEL_INPUT%independ_nvec = independ_nvec
    PARCEL_INPUT%nBlock = nBlock 
    PARCEL_INPUT%CBCG_kstep=CBCG_kstep
    PARCEL_INPUT%CBCG_Eigenflag=CBCG_Eigenflag
    PARCEL_INPUT%power_method_itrmax=power_method_itrmax
    PARCEL_INPUT%CAARNOLDI_sstep=CAARNOLDI_sstep
    PARCEL_INPUT%CAARNOLDI_tstep=CAARNOLDI_tstep
    PARCEL_INPUT%CAARNOLDI_basis=CAARNOLDI_basis
    PARCEL_INPUT%CAARNOLDI_QRflag=CAARNOLDI_QRflag

    call PARCEL_KSP_CRS( &
         icomm, &
         vecx,vecb, &
         n,gn,nnz,istart, &
         crsA,crsRow_ptr, crsCol,&
!         itrmax,rtolmax,abstolmax,&
         PARCEL_INPUT,&
         niter,reshistory,iret)
    itrmax=niter
  end subroutine parcel_dcbcg

  ! CBCG法(DIA)
  subroutine parcel_dcbcg_dia(&
       icomm, &
       vecx,vecb,&
       n,gn,istart,&
       diaA,offset,nnd,&
       ipreflag,ILU_method,addL,iflagAS,&
       itrmax,rtolmax,&
       reshistory,&
       iovlflag, &
       precon_thblock, &
       independ_nvec,  &
       nBlock, &
       CBCG_kstep,CBCG_Eigenflag,power_method_itrmax, &
       CAARNOLDI_sstep,CAARNOLDI_tstep, &
       CAARNOLDI_basis,CAARNOLDI_QRflag, &
       iret &
       )
    use PARCEL_define
    use krylov
    use PARCEL_integer
    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    real*8 vecx(n) ! 解ベクトル
    real*8 vecb(n) ! 右辺ベクトル
    integer(kind=intType) nnd ! 対角数
    real*8 diaA(nnd*n)  ! 行列値
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
    integer CBCG_kstep ! CBCGのステップ数
    integer CBCG_Eigenflag ! CBCG法の固有値計算方法
    integer power_method_itrmax ! べき乗法の反復回数
    integer CAARNOLDI_sstep  ! CAARNOLDIのs step数 
    integer CAARNOLDI_tstep  ! CAARNOLDIのt step数 
    integer CAARNOLDI_basis  ! CAARNOLDIの基底生成方法
    integer CAARNOLDI_QRflag ! CAARNOLDIのQR分解方法

    real*8 reshistory(itrmax) ! 残差履歴
    integer niter
    integer precon_thblock
    integer independ_nvec
    integer nBlock
    integer iret

    type(KSP_INPUT_PARAM) :: PARCEL_INPUT
    abstolmax           = 0.0d0
    solver              = KSP_SOLVER_CBCG
    call INIT_PARCEL_INPUT(&
         PARCEL_INPUT&
         )
    PARCEL_INPUT%ipreflag=ipreflag
    PARCEL_INPUT%ILU_method=ILU_method
    PARCEL_INPUT%addL=addL
    PARCEL_INPUT%iflagAS=iflagAS
    PARCEL_INPUT%solver=solver
    PARCEL_INPUT%iovlflag=iovlflag
    PARCEL_INPUT%itrmax = itrmax
    PARCEL_INPUT%rtolmax = rtolmax
    PARCEL_INPUT%abstolmax = abstolmax
    PARCEL_INPUT%precon_thblock = precon_thblock
    PARCEL_INPUT%independ_nvec = independ_nvec
    PARCEL_INPUT%nBlock = nBlock 
    PARCEL_INPUT%CBCG_kstep=CBCG_kstep
    PARCEL_INPUT%CBCG_Eigenflag=CBCG_Eigenflag
    PARCEL_INPUT%power_method_itrmax=power_method_itrmax
    PARCEL_INPUT%CAARNOLDI_sstep=CAARNOLDI_sstep
    PARCEL_INPUT%CAARNOLDI_tstep=CAARNOLDI_tstep
    PARCEL_INPUT%CAARNOLDI_basis=CAARNOLDI_basis
    PARCEL_INPUT%CAARNOLDI_QRflag=CAARNOLDI_QRflag

    call PARCEL_KSP_DIA( &
         icomm, &
         vecx,vecb, &
         n,gn,istart, &
         diaA,offset,nnd,&
!         itrmax,rtolmax,abstolmax, &
         PARCEL_INPUT,&
         niter,reshistory,iret)
    itrmax=niter

  end subroutine parcel_dcbcg_dia

  ! CAGMRES法(CRS)
  subroutine parcel_dcagmres(&
       icomm, &
       vecx,vecb, &
       n,gn,nnz,istart, &
       crsA,crsRow_ptr,crsCol, &
       ipreflag,ILU_method,addL,iflagAS, &
       itrmax,rtolmax, &
       reshistory, &
       iovlflag, &
       precon_thblock, &
       independ_nvec,  &
       nBlock, &
       CAGMRES_sstep, &
       CAGMRES_tstep, &
       CAGMRES_basis, &
       CAGMRES_QRflag, &
       iret &
       )
    use PARCEL_define
    use krylov
    use PARCEL_integer
    implicit none
    ! -----
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    real*8 vecx(n) ! 解ベクトル
    real*8 vecb(n) ! 右辺ベクトル
    real*8  crsA(nnz) ! 行列値
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
    integer CAGMRES_sstep  ! CAGMRES のs step数
    integer CAGMRES_tstep  ! CAGMRES のt step数
    integer CAGMRES_basis  ! CAGMRES の基底生成方法
    integer CAGMRES_QRflag ! CAGMRES のQR分解方法
    integer precon_thblock
    integer independ_nvec
    integer iret
    integer nBlock
    real*8 reshistory(itrmax) ! 残差履歴
    integer niter
    type(KSP_INPUT_PARAM) :: PARCEL_INPUT

    abstolmax           = 0.0d0
    solver              = KSP_SOLVER_CAGMRES
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
    PARCEL_INPUT%precon_thblock = precon_thblock
    PARCEL_INPUT%independ_nvec = independ_nvec
    PARCEL_INPUT%nBlock = nBlock 
    PARCEL_INPUT%CAGMRES_sstep=CAGMRES_sstep
    PARCEL_INPUT%CAGMRES_tstep=CAGMRES_tstep
    PARCEL_INPUT%CAGMRES_basis=CAGMRES_basis
    PARCEL_INPUT%CAGMRES_QRflag=CAGMRES_QRflag

    call PARCEL_KSP_CRS( &
         icomm, &
         vecx,vecb, &
         n,gn,nnz,istart, &
         crsA,crsRow_ptr, crsCol,&
!         itrmax,rtolmax,abstolmax,&
         PARCEL_INPUT,&
         niter,reshistory,iret)
    itrmax=niter
  end subroutine parcel_dcagmres

  ! CAGMRES法(DIA)
  subroutine parcel_dcagmres_dia(&
       icomm, &
       vecx,vecb, &
       n,gn,istart, &
       diaA,offset,nnd, &
       ipreflag,ILU_method,addL,iflagAS, &
       itrmax,rtolmax, &
       reshistory, &
       iovlflag, &
       precon_thblock, &
       independ_nvec,  &
       nBlock, &
       CAGMRES_sstep, &
       CAGMRES_tstep, &
       CAGMRES_basis, &
       CAGMRES_QRflag, &
       iret &
       )
    use PARCEL_define
    use krylov
    use PARCEL_integer
    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    real*8 vecx(n) ! 解ベクトル
    real*8 vecb(n) ! 右辺ベクトル
    integer(kind=intType) nnd ! 対角数
    real*8 diaA(nnd*n)  ! 行列値
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
    integer CAGMRES_sstep  ! CAGMRES のs step数
    integer CAGMRES_tstep  ! CAGMRES のt step数
    integer CAGMRES_basis  ! CAGMRES の基底生成方法
    integer CAGMRES_QRflag ! CAGMRES のQR分解方法
    real*8 reshistory(itrmax) ! 残差履歴
    integer niter
    integer precon_thblock
    integer independ_nvec
    integer nBlock
    integer iret

    type(KSP_INPUT_PARAM) :: PARCEL_INPUT

    abstolmax           = 0.0d0
    solver              = KSP_SOLVER_CAGMRES

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
    PARCEL_INPUT%precon_thblock = precon_thblock
    PARCEL_INPUT%independ_nvec = independ_nvec
    PARCEL_INPUT%nBlock = nBlock 
    PARCEL_INPUT%CAGMRES_sstep=CAGMRES_sstep
    PARCEL_INPUT%CAGMRES_tstep=CAGMRES_tstep
    PARCEL_INPUT%CAGMRES_basis=CAGMRES_basis
    PARCEL_INPUT%CAGMRES_QRflag=CAGMRES_QRflag

    call PARCEL_KSP_DIA( &
         icomm, &
         vecx,vecb, &
         n,gn,istart, &
         diaA,offset,nnd,&
!         itrmax,rtolmax,abstolmax, &
         PARCEL_INPUT,&
         niter,reshistory,iret)
    itrmax=niter

  end subroutine parcel_dcagmres_dia

  ! CAARNOLDI法(CRS)
  subroutine parcel_dcaarnoldi(&
       icomm, &
       vecx,vecb, &
       n,gn,nnz,istart, &
       crsA,crsRow_ptr,crsCol, &
       ipreflag,ILU_method,addL,iflagAS, &
!       itrmax, &
       iovlflag, &
       precon_thblock, &
       independ_nvec,  &
       nBlock, &
       CAARNOLDI_sstep, &
       CAARNOLDI_tstep, &
       CAARNOLDI_basis, &
       CAARNOLDI_QRflag, &
       EigenValueReal, &
       EigenValueImag, &
       EigenVector, &
       EigenErrorNorm, &
       MaxEigenID, &
       MinEigenID, &
       iret &
       )
    use PARCEL_define
    use krylov
    use PARCEL_integer
    implicit none
    ! -----
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    real*8 vecx(n) ! 解ベクトル
    real*8 vecb(n) ! 右辺ベクトル
    real*8  crsA(nnz) ! 行列値
    integer(kind=intType) crsRow_ptr(n+1) ! ポインタ 
    integer(kind=intType) crsCol(nnz) ! 列番号
    integer ipreflag ! 前処理方法　
    integer ILU_method ! D-ILU分解の方法
    integer addL ! 加法シュワルツのオーバラップ領域の大きさ
    integer iflagAS ! 加法シュワルツのオーバラップ領域補間方法
    real*8 rtolmax  !  絶対残差 (収束条件)
    real*8 abstolmax  !  絶対残差 (収束条件)
    integer solver ! 線形連立方程式の解法
    integer iovlflag ! 疎行列ベクトル積の通信・計算オーバーラップ方法
    integer CAARNOLDI_sstep  ! CAARNOLDI のs step数
    integer CAARNOLDI_tstep  ! CAARNOLDI のt step数
    integer CAARNOLDI_basis  ! CAARNOLDI の基底生成方法
    integer CAARNOLDI_QRflag ! CAARNOLDI のQR分解方法
    integer itrmax
    real*8 reshistory(CAARNOLDI_sstep*CAARNOLDI_tstep+1) ! 残差履歴
    integer precon_thblock
    integer independ_nvec
    integer iret
    integer nBlock

    real*8 EigenValueReal(CAARNOLDI_sstep*CAARNOLDI_tstep)
    real*8 EigenValueImag(CAARNOLDI_sstep*CAARNOLDI_tstep)
    real*8 EigenVector(CAARNOLDI_sstep*CAARNOLDI_tstep)
    real*8 EigenErrorNorm(CAARNOLDI_sstep*CAARNOLDI_tstep)
    integer MaxEigenID
    integer MinEigenID

    integer niter
    type(KSP_INPUT_PARAM) :: PARCEL_INPUT

    abstolmax           = 0.0d0
    rtolmax             = 0.0d0
    itrmax              = CAARNOLDI_sstep * CAARNOLDI_tstep
    solver              = KSP_SOLVER_CAARNOLDI
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
    PARCEL_INPUT%precon_thblock = precon_thblock
    PARCEL_INPUT%independ_nvec = independ_nvec
    PARCEL_INPUT%nBlock = nBlock 
    PARCEL_INPUT%CAARNOLDI_sstep=CAARNOLDI_sstep
    PARCEL_INPUT%CAARNOLDI_tstep=CAARNOLDI_tstep
    PARCEL_INPUT%CAARNOLDI_basis=CAARNOLDI_basis
    PARCEL_INPUT%CAARNOLDI_QRflag=CAARNOLDI_QRflag

    call PARCEL_KSP_CRS( &
         icomm, &
         vecx,vecb, &
         n,gn,nnz,istart, &
         crsA,crsRow_ptr, crsCol,&
         PARCEL_INPUT,&
         niter,reshistory,iret,&
         EigenValueReal, &
         EigenValueImag, &
         EigenVector, &
         EigenErrorNorm, &
         MaxEigenID, &
         MinEigenID  )

!    itrmax=niter
  end subroutine parcel_dcaarnoldi

  ! CAARNOLDI法(DIA)
  subroutine parcel_dcaarnoldi_dia(&
       icomm, &
       vecx,vecb, &
       n,gn,istart, &
       diaA,offset,nnd, &
       ipreflag,ILU_method,addL,iflagAS, &
!       itrmax,&
       iovlflag, &
       precon_thblock, &
       independ_nvec,  &
       nBlock, &
       CAARNOLDI_sstep, &
       CAARNOLDI_tstep, &
       CAARNOLDI_basis, &
       CAARNOLDI_QRflag, &
       EigenValueReal, &
       EigenValueImag, &
       EigenVector, &
       EigenErrorNorm, &
       MaxEigenID, &
       MinEigenID, &
       iret &
       )
    use PARCEL_define
    use krylov
    use PARCEL_integer
    implicit none
    integer icomm ! MPI コミュニケータ
    integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
    integer(kind=intType) gn ! ベクトル全体のサイズ
    integer(kind=intType) nnz ! 非零要素数
    integer(kind=intType) istart ! ベクトルの開始位置
    real*8 vecx(n) ! 解ベクトル
    real*8 vecb(n) ! 右辺ベクトル
    integer(kind=intType) nnd ! 対角数
    real*8 diaA(nnd*n)  ! 行列値
    integer(kind=intType) offset(nnd) ! 対角オフセット
    integer ipreflag ! 前処理方法　
    integer ILU_method ! D-ILU分解の方法
    integer addL ! 加法シュワルツのオーバラップ領域の大きさ
    integer iflagAS ! 加法シュワルツのオーバラップ領域補間方法
    real*8 rtolmax  !  絶対残差 (収束条件)
    real*8 abstolmax  !  絶対残差 (収束条件)
    integer solver ! 線形連立方程式の解法
    integer iovlflag ! 疎行列ベクトル積の通信・計算オーバーラップ方法
    integer CAARNOLDI_sstep  ! CAARNOLDI のs step数
    integer CAARNOLDI_tstep  ! CAARNOLDI のt step数
    integer CAARNOLDI_basis  ! CAARNOLDI の基底生成方法
    integer CAARNOLDI_QRflag ! CAARNOLDI のQR分解方法

    integer niter
    integer precon_thblock
    integer independ_nvec
    integer nBlock

    real*8 EigenValueReal(CAARNOLDI_sstep*CAARNOLDI_tstep)
    real*8 EigenValueImag(CAARNOLDI_sstep*CAARNOLDI_tstep)
    real*8 EigenVector(CAARNOLDI_sstep*CAARNOLDI_tstep)
    real*8 EigenErrorNorm(CAARNOLDI_sstep*CAARNOLDI_tstep)
    integer MaxEigenID
    integer MinEigenID

    integer iret
    integer itrmax ! 最大反復回数
    real*8 reshistory(CAARNOLDI_sstep*CAARNOLDI_tstep) ! 残差履歴

    type(KSP_INPUT_PARAM) :: PARCEL_INPUT

    abstolmax           = 0.0d0
    rtolmax             = 0.0d0
    reshistory          = 0.0d0
    itrmax              = CAARNOLDI_sstep * CAARNOLDI_tstep
    solver              = KSP_SOLVER_CAARNOLDI

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
    PARCEL_INPUT%precon_thblock = precon_thblock
    PARCEL_INPUT%independ_nvec = independ_nvec
    PARCEL_INPUT%nBlock = nBlock 
    PARCEL_INPUT%CAARNOLDI_sstep=CAARNOLDI_sstep
    PARCEL_INPUT%CAARNOLDI_tstep=CAARNOLDI_tstep
    PARCEL_INPUT%CAARNOLDI_basis=CAARNOLDI_basis
    PARCEL_INPUT%CAARNOLDI_QRflag=CAARNOLDI_QRflag

    call PARCEL_KSP_DIA( &
         icomm, &
         vecx,vecb, &
         n,gn,istart, &
         diaA,offset,nnd,&
!         itrmax,rtolmax,abstolmax, &
         PARCEL_INPUT,&
         niter,reshistory,iret,&
         EigenValueReal, &
         EigenValueImag, &
         EigenVector, &
         EigenErrorNorm, &
         MaxEigenID, &
         MinEigenID  )
!    itrmax=niter

  end subroutine parcel_dcaarnoldi_dia


