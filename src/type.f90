! プログラムのリファクタリング
! matrix,krylovはDIA,CRS分ける
! todo :: DDMを8バイト整数型に対応させる

! -----
module ddmMat_wrap
  implicit none
  type ddm_param
     integer stmx, stmy, stmz
     integer nx, ny, nz, nxy, nxyz
     integer mx, my, mz, mxy
     integer num_neighbor_ranks
     integer, allocatable :: neighbor_ranks(:)
     integer, allocatable :: comm_flag(:)
     integer npe_x, npe_y, npe_z
     integer leng_halo
     integer ncomp_grids
     integer div_direc_th
     integer precon_thblock
     integer iflag_DDM_independ
     integer underlap_step
     
  end type ddm_param

  type ddmmatrix_wrap
     integer mat_type ! 1:CRS 2:DIA 3:DDM 
     integer ndim
     integer nrow

     !parameters for DDM
     integer nnd_dia
     integer, allocatable :: ioff_dia(:)
     double precision, allocatable :: val_dia(:)
     double precision, allocatable :: val_dia_lo(:)
     integer, allocatable :: ioff_grids(:,:)
     type(ddm_param) :: prm
  end type ddmmatrix_wrap

  ! 加法シュワルツ前処理に特有のパラメータ
  !
  ! 配列の6成分は それぞれ -x, +x, -y, +y, -z, +z 方向の隣接領域に対応する
  !
  type comm_buff_id
     integer, allocatable :: id(:)
  end type comm_buff_id

  type AS_geom
     integer ndist
     integer addL
     integer sleeve_mx , sleeve_px  ! -x, +x 方向の重なり幅（隣接する領域がない場合は0になる）
     integer sleeve_my , sleeve_py  !
     integer sleeve_mz , sleeve_pz  !
     integer grids_s(6), grids_r(6) ! 通信領域の全グリッド数
     integer rank_s (6), rank_r (6) ! 通信相手のプロセスのランク
     integer x0s    (6), y0s    (6), z0s   (6) ! 通信領域の原点（s = send, r = receive）
     integer x0r    (6), y0r    (6), z0r   (6) !
     integer xlen_s (6), ylen_s (6), zlen_s(6) ! 通信領域の幅
     integer xlen_r (6), ylen_r (6), zlen_r(6) !
     type(comm_buff_id) send_from_id(6),recv_from_id(6)
     type(comm_buff_id) send_to_id(6),recv_to_id(6)

  end type AS_geom

end module ddmMat_wrap

module param_AS_type_ddm
  use PARCEL_define
  implicit none

  ! DDM前処理方向
  ! (x->y->z :: 1->2->3)
  integer,parameter :: KSP_DDM_DIVIDE_THREADS_DIRECT_X=1 ! X方向 
  integer,parameter :: KSP_DDM_DIVIDE_THREADS_DIRECT_Y=2 ! Y方向
  integer,parameter :: KSP_DDM_DIVIDE_THREADS_DIRECT_Z=3 ! Z方向
  integer,parameter :: KSP_DDM_DIVIDE_THREADS_DIRECT_NONE=0 ! スレッド分割なし

  type precon_AS_DDM
     integer nxs_local_th,nxe_local_th
     integer nys_local_th,nye_local_th
     integer nzs_local_th,nze_local_th

     integer nxs_th,nxe_th
     integer nys_th,nye_th
     integer nzs_th,nze_th

     integer nz_th
     integer ny_th
     integer nx_th

     integer mz_th
     integer my_th
     integer mx_th

     integer mxy_th
     integer m_th

     real*8,allocatable  :: Dinv(:)
     real*8,allocatable  :: Dinv_lo(:)
     ! 前処理前のベクトル
     real*8,allocatable  :: rvec(:)
     real*8,allocatable  :: rvec_lo(:)
     ! 前処理後のベクトル
     real*8,allocatable  :: svec(:)
     real*8,allocatable  :: svec_lo(:)
  end type precon_AS_DDM

  type precon_AS_wrap
     integer max_stencil
     integer max_threads
!     integer precon_thblock
     integer AS_METHOD_TYPE
     type(precon_AS_DDM),allocatable :: AS_DDM(:)
  end type precon_AS_wrap

end module param_AS_type_ddm

module DDM_TYPE
  use ddmMat_wrap
  use param_AS_type_ddm
  implicit none

  type DDM_ALL_TYPE
     type(precon_AS_wrap), pointer ::  precon
     integer ndim !dimension of global matrix
     integer n !number of rows of local matrix  
     integer n0
     integer nnd
     integer istart !start row index in global matrix
     type(ddmmatrix_wrap), pointer :: mat
     double precision, pointer :: b(:)
     double precision, pointer :: x(:) 
     double precision, pointer :: b_lo(:)
     double precision, pointer :: x_lo(:) 
     integer ncpu
     integer myrank
     integer icomm !MPI communicator
     double precision, pointer :: diag(:)
     double precision, pointer :: diag_lo(:)
     integer, pointer :: oloop !outer loop
     double precision tol !tolerance
     integer ipc !precondition flag
     integer nres !number of residual history
     double precision,pointer :: res(:) !residual history
     integer iret !result status
     integer commHide         ! 0: no hiding, 1: one time, 2: each time
     integer AS_len           ! additive Schwarz overlap length
     integer AS_mean
     integer AS_n
     type(AS_geom), pointer :: AS_g
     double precision, pointer :: AS_val_dia(:)
     double precision, pointer :: AS_val_dia_lo(:)
     integer,pointer :: AS_ioff_dia(:)
     type(ddm_param), pointer :: AS_prm
     double precision,pointer ::  AS_v1(:), AS_v2(:)
     double precision,pointer ::  AS_v1_lo(:), AS_v2_lo(:)

     integer solver
     ! -----
     ! CBCG
     integer k
     integer bflag 
     double precision EigenMax_Real,EigenMax_Imag
     double precision EigenMin_Real,EigenMin_Imag
     integer CBCG_Eigenflag

     ! -----
     ! CAGMRES
     integer sstep
     integer tstep
     integer qrflag
     integer basis

  end type DDM_ALL_TYPE

end module DDM_TYPE

! -----

module comm
  ! MPI,OpenMPに付随するものは4バイト整数型とする
  ! myrank,
  ! Urank_AS,Lrank_AS 
  ! rank
  ! comm
  ! 
  ! numth,myth,max_numth
  ! 
  ! -----

  use mpi
  use PARCEL_integer
  implicit none

  type nrk
     integer rank
     !     integer(kind=intType) ndata 
     integer ndata 
     integer(kind=intType),allocatable :: offset(:)

     integer(kind=intType) nnz
     integer(kind=intType),allocatable :: nrk_ptr(:)
     integer(kind=intType),allocatable :: nrkRow(:)
     integer(kind=intType),allocatable :: nrkCol(:)

  end type nrk
  
  type comm_Vec
     integer comm
     integer myrank
     integer npes
     integer(kind=intType) n
     !     integer(kind=intType) idwin
     !     integer(kind=MPI_ADDRESS_KIND) imem_size ! バイト
     !     real*8,allocatable  :: buff(:) ! mem window 
     integer(kind=intType),allocatable :: istart(:)

     ! 受信データ情報
     integer nrk_rank
     integer nrk_data
     integer nrk_data_max
     type(nrk),allocatable :: nrk_prm(:)

     ! 送信データ情報
     integer nrk_rank_send
     integer nrk_data_send
     integer nrk_data_max_send
     type(nrk),allocatable :: nrk_prm_send(:)
     !     integer(kind=intType) win_info
  end type comm_Vec

end module comm

module matrix
  use PARCEL_integer
  implicit none
  ! ブロックサイズ --> 未調整
  !integer,parameter :: nBlock
  !  integer,save :: nBlock=200
  ! nBlock0は一時的なもの4倍精度の修正時に廃止する
  !  integer,save :: nBlock0=2000
  real*8,parameter :: minReal8=0.0d0

  type prmMat
     ! ローカル行サイズ
     integer(kind=intType) n
     ! グローバル行サイズ
     integer(kind=intType) gn
     ! 受け持つグローバル行の開始,終了番号    
     integer(kind=intType) istart,iend
     ! 非零要素数
     integer(kind=intType) nnz
     ! 対角行の本数
     integer(kind=intType) nnd
     ! -----
     ! MatFormat = 1 :: CRS
     ! MatFormat = 2 :: DIA
     ! MatFormat = 3 :: DDA
     integer MatFormat 

  end type prmMat

  type DEN_Mat
     ! 行指向 
     integer(kind=intType) RowN,ColN
     real*8,allocatable  :: denA(:) 

     real*8,allocatable  :: denA_hi(:) 
     real*8,allocatable  :: denA_lo(:) 
  end type DEN_Mat

  type COO_Mat
     real*8,allocatable  :: cooA(:)   ! 行列値
     real*8,allocatable  :: cooA_hi(:)   ! 行列値
     real*8,allocatable  :: cooA_lo(:)   ! 行列値
     integer(kind=intType),allocatable :: cooRow(:) ! 行番号 
     integer(kind=intType),allocatable :: cooCol(:) ! 列番号
  end type COO_Mat

  type CRS_Mat
     real*8,allocatable  :: crsA(:)       ! 行列値
     real*8,allocatable  :: crsA_hi(:)    ! 行列値
     real*8,allocatable  :: crsA_lo(:)    ! 行列値
     integer(kind=intType),allocatable :: crsRow_ptr(:) ! ポインタ 
     integer(kind=intType),allocatable :: crsRow_ptr_Diag_s(:) ! ポインタ 
     integer(kind=intType),allocatable :: crsRow_ptr_Diag_e(:) ! ポインタ 
     integer(kind=intType),allocatable :: crsCol(:)     ! 列番号

     integer(kind=intType),allocatable :: RowL(:)
     integer(kind=intType),allocatable :: RowD(:)
     integer(kind=intType),allocatable :: RowU(:)

     integer(kind=intType),allocatable :: iLs(:)
     integer(kind=intType) nColL,nColU ! L,U側の列幅
  end type CRS_Mat

  type HybCRS_Mat
     real*8,allocatable  :: crsAL(:),crsAD(:),crsAU(:) ! 行列値
     real*8,allocatable  :: crsAL_hi(:),crsAD_hi(:),crsAU_hi(:) ! 行列値
     real*8,allocatable  :: crsAL_lo(:),crsAD_lo(:),crsAU_lo(:) ! 行列値
     integer(kind=intType),allocatable :: crsRow_ptr_L(:),crsRow_ptr_D(:),crsRow_ptr_U(:) ! ポインタ
     integer(kind=intType),allocatable :: crsColL(:),crsColD(:),crsColU(:)  ! 列番号

     integer(kind=intType),allocatable :: RowL(:),RowD(:),RowU(:) ! 行番号
     integer(kind=intType) nnzL,nnzD,nnzU
     integer(kind=intType) nColL,nColD,nColU

  end type HybCRS_Mat

  type DIA_Mat
     real*8,allocatable  :: diaA(:)    ! 行列値
     real*8,allocatable  :: diaA_hi(:) ! 行列値
     real*8,allocatable  :: diaA_lo(:) ! 行列値
     integer(kind=intType),allocatable :: offset(:)  ! 対角位置からのオフセット
     integer(kind=intType),allocatable :: iLs(:)

     integer(kind=intType),allocatable :: offsetVec(:) ! ベクトルのオフセット
     real*8,allocatable :: Vec(:)    ! DIA形式の行列ベクトル積のベクトル
     real*8,allocatable :: Vec_hi(:) ! DIA形式の行列ベクトル積のベクトル
     real*8,allocatable :: Vec_lo(:) ! DIA形式の行列ベクトル積のベクトル
     integer(kind=intType),allocatable :: RowL(:),RowD(:),RowU(:) ! 行番号
     integer(kind=intType) nCol
     integer(kind=intType) nColL,nColD,nColU

  end type DIA_Mat

end module matrix

module precon
  use matrix
  use PARCEL_integer
  implicit none

  type ILU_regionMat
     integer(kind=intType),allocatable :: ILU_Lstart(:),ILU_Lend(:)
     integer(kind=intType),allocatable :: ILU_Ustart(:),ILU_Uend(:)
  end type ILU_regionMat

  type precon_AS
     integer(kind=intType) istart_th 
     integer(kind=intType) iend_th   
     integer(kind=intType) n_startth 
     integer(kind=intType) n_endth   

     integer(kind=intType) istart0_th 
     integer(kind=intType) iend0_th   
     integer(kind=intType) n0_startth 
     integer(kind=intType) n0_endth   

     real*8,allocatable  :: Dinv(:)
     real*8,allocatable  :: Dinv_hi(:)
     real*8,allocatable  :: Dinv_lo(:)
     ! 前処理前のベクトル
     real*8,allocatable  :: rvec(:)
     real*8,allocatable  :: rvec_hi(:)
     real*8,allocatable  :: rvec_lo(:)
     ! 前処理後のベクトル
     real*8,allocatable  :: svec(:)
     real*8,allocatable  :: svec_hi(:)
     real*8,allocatable  :: svec_lo(:)
  end type precon_AS

  type precon_ILU

     integer(kind=intType) n_th 
     integer(kind=intType) istart_th 
     integer(kind=intType) iend_th 

     type(CRS_Mat) crsMatLU
     integer(kind=intType) nnz_th     
     type(CRS_Mat) crsMatL
     integer(kind=intType) nnzL_th
     type(CRS_Mat) crsMatU
     integer(kind=intType) nnzU_th

     type(DIA_Mat) diaMatLU
     integer(kind=intType) nnd_th     
     type(DIA_Mat) diaMatL
     integer(kind=intType) nndL_th
     type(DIA_Mat) diaMatU
     integer(kind=intType) nndU_th

  end type precon_ILU

end module precon

module krylov
  use comm
  use matrix
  use precon
  use DDM_TYPE,only : DDM_ALL_TYPE
  implicit none

  type KSP
     type (prmMat) prm
     type(comm_Vec) :: commVec
     real*8,allocatable :: buff(:)
     real*8,allocatable :: buff_hi(:)
     real*8,allocatable :: buff_lo(:)
     real*8,allocatable :: rbuffL(:)
     real*8,allocatable :: rbuffL_hi(:)
     real*8,allocatable :: rbuffL_lo(:)
     real*8,allocatable :: rbuffU(:)
     real*8,allocatable :: rbuffU_hi(:)
     real*8,allocatable :: rbuffU_lo(:)
     integer max_numth
     real*8 EigenMax_Real,EigenMax_Imag
     real*8 EigenMin_Real,EigenMin_Imag
     real*8,allocatable :: EigenValueReal(:)
     real*8,allocatable :: EigenValueImag(:)
     real*8,allocatable :: EigenVector(:)
     real*8,allocatable :: EigenErrorNorm(:)
     integer MaxEigenID
     integer MinEigenID

     ! -----
     ! 精度選択
     integer PRECISION_A
     integer PRECISION_b
     integer PRECISION_x
     integer PRECISION_KRYLOV 
     integer PRECISION_PRECONDITION ! 前処理行列の精度 (前処理の精度ではない)

     ! クリロフ部分空間法の選択
     integer solver
     ! 行列ベクトル積の通信,計算オーバーラップ有無
     integer iovlflag

     ! 対応する疎行列形式
     type(CRS_Mat)    :: crsMat 
     type(DIA_Mat)    :: diaMat 
     type(HybCRS_Mat) :: HybcrsMat
     type(DDM_ALL_TYPE) :: ddmMat

     real*8,allocatable :: r(:)
     real*8,allocatable :: r_hi(:)
     real*8,allocatable :: r_lo(:)
     real*8,allocatable :: x(:)
     real*8,allocatable :: x_hi(:)
     real*8,allocatable :: x_lo(:)
     real*8,allocatable :: b(:)
     real*8,allocatable :: b_hi(:)
     real*8,allocatable :: b_lo(:)

     ! 最大反復回数
     integer itrmax
     ! 収束判定条件 :: 相対残差ノルム ||b-Ax||/||b||
     real*8 rtolmax
     ! 収束判定条件 :: 絶対残差ノルム ||b-Ax||
     real*8 abstolmax

     ! 残差履歴
     real*8,allocatable :: reshistory(:)

     ! 反復回数
     integer niter

     ! 前処理フラグ
     integer ipreflag
     real*8,allocatable :: Dinv(:)
     real*8,allocatable :: Dinv_hi(:)
     real*8,allocatable :: Dinv_lo(:)

     ! 不完全LU分解の対角行列の決定方法
     integer ILU_method

     ! ILU分解の行列の参照範囲
     type(ILU_regionMat),allocatable :: ILUregionMat(:)

     ! 加法シュワルツの重なり幅
     integer addL
     real*8,allocatable :: ivec_AS(:)
     real*8,allocatable :: ivec_AS_hi(:)
     real*8,allocatable :: ivec_AS_lo(:)
     real*8,allocatable :: ovec_AS(:)
     real*8,allocatable :: ovec_AS_hi(:)
     real*8,allocatable :: ovec_AS_lo(:)

     ! SpMVの際に1スレッドが担当するタスクの最小行数.
     !     integer,save :: nBlock=2000
     integer nBlock
     integer nvecBlock

     !  加法シュワルツの通信相手
     integer Urank_AS,Lrank_AS ! 通信相手がいない場合は-1

     ! ILU分解用構造体
     type(precon_ILU),allocatable ::  preconILU(:) 

     ! 細分化前処理で独立に計算可能な長さ
     integer independ_nvec

     ! 前処理でプロセス内で分割するブロック数(1ブロック1スレッドが担当)
     integer precon_thblock

     ! スレッド並列加法シュワルツ計算用構造体
     type(precon_AS),allocatable ::  preconAS(:) 
     type(prmMat)  :: prm_AS
     type(dia_Mat) :: precon_diaMatAS
     type(CRS_Mat) :: precon_crsMatAS
     integer iflagAS
     ! -----
     integer GMRES_m ! GMRES(m) リスタート数
     integer GMRES_GSflag ! GMRES(m)のグラムシュミット直交化方法

     integer CAGMRES_sstep  ! CAGMRESのs step数 
     integer CAGMRES_tstep  ! CAGMRESのt step数 
     integer CAGMRES_basis  ! CAGMRESの基底生成方法
     integer CAGMRES_QRflag ! CAGMRESのQR分解方法

     integer CAARNOLDI_sstep  ! CAARNOLDIのs step数 
     integer CAARNOLDI_tstep  ! CAARNOLDIのt step数 
     integer CAARNOLDI_basis  ! CAARNOLDIの基底生成方法
     integer CAARNOLDI_QRflag ! CAARNOLDIのQR分解方法

     integer CALANCZOS_sstep  ! CALANCZOSのs step数 
     integer CALANCZOS_tstep  ! CALANCZOSのt step数 
     integer CALANCZOS_basis  ! CALANCZOSの基底生成方法
     integer CALANCZOS_QRflag ! CALANCZOSのQR分解方法
     integer EIGflag ! これは仮の変数
     integer EIG_calcflag ! これは仮の数変
     integer CBCG_kstep ! CBCGのステップ数
     integer power_method_itrmax ! べき乗法の反復回数
     integer CBCG_Eigenflag ! CBCG法の固有値計算方法

     integer iret ! リターンコード
  end type KSP

  type KSP_INPUT_PARAM
     integer ipreflag
     integer ILU_method
     integer addL
     integer iflagAS
     integer solver
     integer iovlflag
     integer GMRES_m
     integer GMRES_GSflag
     integer CAGMRES_sstep
     integer CAGMRES_tstep
     integer CAGMRES_basis
     integer CAGMRES_QRflag
     integer CBCG_kstep
     integer CBCG_Eigenflag
     integer power_method_itrmax
     integer CAARNOLDI_sstep
     integer CAARNOLDI_tstep
     integer CAARNOLDI_basis
     integer CAARNOLDI_QRflag

     integer itrmax
     real*8  rtolmax
     real*8  abstolmax
     integer independ_nvec
     integer precon_thblock
     integer nBlock
     integer PRECISION_A
     integer PRECISION_b
     integer PRECISION_x
     integer PRECISION_PRECONDITION

  end type KSP_INPUT_PARAM
contains

  subroutine PARCEL_KSP_Default_Setting(method)
    use PARCEL_define
    implicit none
    type(KSP) :: method
    integer solver_default
    integer addL_default
    integer iflagAS_default
    integer iovlflag_default
    integer GMRES_m_default
    integer GMRES_GSflag_default
    integer CAGMRES_sstep_default 
    integer CAGMRES_tstep_default 
    integer CAGMRES_basis_default 
    integer CAGMRES_QRflag_default

    integer CAARNOLDI_sstep_default 
    integer CAARNOLDI_tstep_default 
    integer CAARNOLDI_basis_default 
    integer CAARNOLDI_QRflag_default

    integer CALANCZOS_sstep_default 
    integer CALANCZOS_tstep_default 
    integer CALANCZOS_basis_default 
    integer CALANCZOS_QRflag_default

    integer nBlock_default
    integer nvecBlock_default
    integer CBCG_kstep_default
    integer power_method_itrmax_default
    integer CBCG_Eigenflag_default
    integer ipreflag_default
    integer ILU_method_default
    integer independ_nvec_default
    integer precon_thblock_default
    integer itrmax_default
    real*8  rtolmax_default
    real*8  abstolmax_default
    integer PRECISION_A_default
    integer PRECISION_b_default
    integer PRECISION_x_default
    integer PRECISION_PRECONDITION_default
    integer PRECISION_KRYLOV_default

    solver_default       = KSP_SOLVER_GMRES
    addL_default         = 0
    iflagAS_default      = KSP_ADDITIVESCHWARZ_TYPE_NONE
    iovlflag_default     = KSP_SPMV_OVL_NONE
    GMRES_m_default      = 4
    GMRES_GSflag_default = KSP_GMRES_GS_TYPE_MGS

    CAGMRES_sstep_default  = 4
    CAGMRES_tstep_default  = 1
    CAGMRES_basis_default  = KSP_MPK_MONOMIAL_BASIS
    CAGMRES_QRflag_default = KSP_CAGMRES_QR_TYPE_CHOLESKYQR

    CAARNOLDI_sstep_default  = 4
    CAARNOLDI_tstep_default  = 10
    !    CAARNOLDI_tstep_default  = 1
    CAARNOLDI_basis_default  = KSP_MPK_NEWTON_BASIS
    !    CAARNOLDI_basis_default  = KSP_MPK_MONOMIAL_BASIS
    CAARNOLDI_QRflag_default = KSP_CAARNOLDI_QR_TYPE_TSQR
    !    CAARNOLDI_QRflag_default = KSP_QRDECOMPOSITION_CHOLESKYQR

    CALANCZOS_sstep_default  = 4
    CALANCZOS_tstep_default  = 10
    CALANCZOS_basis_default  = KSP_MPK_NEWTON_BASIS
    CALANCZOS_QRflag_default = KSP_CALANCZOS_QR_TYPE_TSQR

    CBCG_kstep_default = 4
    power_method_itrmax_default = 100
    CBCG_Eigenflag_default = KSP_EIGENVALUE_TYPE_CAARNOLDI

    ipreflag_default     = KSP_PRECONDITION_NONE
    ILU_method_default   = KSP_PRECONDITION_ILUD 
    !    ILU_method_default   = KSP_PRECONDITION_MLUD 
    independ_nvec_default = 1
    precon_thblock_default = -1

    itrmax_default       = 100
    rtolmax_default      = 1.0d0-8
    abstolmax_default    = 1.0d0-50

    nBlock_default = 2000
    nvecBlock_default = 100
#if 0
    PRECISION_A_default            = KSP_DOUBLE_PRECISION
    PRECISION_b_default            = KSP_DOUBLE_PRECISION
    PRECISION_x_default            = KSP_DOUBLE_PRECISION
    PRECISION_PRECONDITION_default = KSP_DOUBLE_PRECISION
#else
    PRECISION_A_default            = KSP_DOUBLE_DOUBLE_PRECISION
    PRECISION_b_default            = KSP_DOUBLE_DOUBLE_PRECISION
    PRECISION_x_default            = KSP_DOUBLE_DOUBLE_PRECISION
    PRECISION_PRECONDITION_default = KSP_DOUBLE_DOUBLE_PRECISION
!    PRECISION_KRYLOV_default       = KSP_DOUBLE_DOUBLE_PRECISION
    PRECISION_KRYLOV_default       = KSP_DOUBLE_PRECISION
#endif

    ! -----

    call set_KSP_solver(method,solver_default)
    call set_PRECISION_PRECONDITION(method,PRECISION_PRECONDITION_default)
    call set_PRECISION_A(method,PRECISION_A_default)
    call set_PRECISION_x(method,PRECISION_x_default)
    call set_PRECISION_b(method,PRECISION_b_default)
    call set_PRECISION_KRYLOV(method,PRECISION_KRYLOV_default)

    call set_KSP_addL(method,addL_default)
    call set_KSP_nBlock(method,nBlock_default)
    call set_KSP_nvecBlock(method,nvecBlock_default)
    call set_KSP_iflagAS(method,iflagAS_default)
    call set_KSP_iovlflag(method,iovlflag_default)
    call set_KSP_GMRES_m(method,GMRES_m_default)
    call set_KSP_GMRES_GSflag(method,GMRES_GSflag_default)

    call set_KSP_CAGMRES_sstep(method,CAGMRES_sstep_default)
    call set_KSP_CAGMRES_tstep(method,CAGMRES_tstep_default)
    call set_KSP_CAGMRES_basis(method,CAGMRES_basis_default)
    call set_KSP_CAGMRES_QRflag(method,CAGMRES_QRflag_default)

    call set_KSP_CAARNOLDI_sstep(method,CAARNOLDI_sstep_default)
    call set_KSP_CAARNOLDI_tstep(method,CAARNOLDI_tstep_default)
    call set_KSP_CAARNOLDI_basis(method,CAARNOLDI_basis_default)
    call set_KSP_CAARNOLDI_QRflag(method,CAARNOLDI_QRflag_default)

    call set_KSP_CALANCZOS_sstep(method,CALANCZOS_sstep_default)
    call set_KSP_CALANCZOS_tstep(method,CALANCZOS_tstep_default)
    call set_KSP_CALANCZOS_basis(method,CALANCZOS_basis_default)
    call set_KSP_CALANCZOS_QRflag(method,CALANCZOS_QRflag_default)

    call set_KSP_CBCG_kstep(method,CBCG_kstep_default)
    call set_KSP_power_method_itrmax(method,power_method_itrmax_default)
    call set_KSP_CBCG_Eigenflag(method,CBCG_Eigenflag_default)

    call set_KSP_ipreflag(method,ipreflag_default)
    call set_KSP_ILU_method(method,ILU_method_default)

    call set_KSP_independ_nvec(method,independ_nvec_default)
    call set_KSP_precon_thblock(method,precon_thblock_default)

    call set_KSP_itrmax(method,itrmax_default)
    call set_KSP_rtolmax(method,rtolmax_default)
    call set_KSP_abstolmax(method,abstolmax_default)

  end subroutine PARCEL_KSP_Default_Setting

  subroutine set_KSP_solver(method,solver)
    implicit none
    type(KSP) :: method
    integer solver
    method%solver=solver
  end subroutine set_KSP_solver

  subroutine set_PRECISION_A(method,PRECISION_A)
    implicit none
    type(KSP) :: method
    integer PRECISION_A
    method%PRECISION_A=PRECISION_A
  end subroutine set_PRECISION_A

  subroutine set_PRECISION_x(method,PRECISION_x)
    implicit none
    type(KSP) :: method
    integer PRECISION_x
    method%PRECISION_x=PRECISION_x
  end subroutine set_PRECISION_x

  subroutine set_PRECISION_b(method,PRECISION_b)
    implicit none
    type(KSP) :: method
    integer PRECISION_b
    method%PRECISION_b=PRECISION_b
  end subroutine set_PRECISION_b

  subroutine set_PRECISION_KRYLOV(method,PRECISION_KRYLOV)
    implicit none
    type(KSP) :: method
    integer PRECISION_KRYLOV
    method%PRECISION_KRYLOV=PRECISION_KRYLOV
  end subroutine set_PRECISION_KRYLOV

  subroutine set_PRECISION_PRECONDITION(method,PRECISION_PRECONDITION)
    implicit none
    type(KSP) :: method
    integer PRECISION_PRECONDITION
    method%PRECISION_PRECONDITION=PRECISION_PRECONDITION
  end subroutine set_PRECISION_PRECONDITION

  subroutine set_KSP_nBlock(method,nBlock)
    implicit none
    type(KSP) :: method
    integer nBlock
    if(0<nBlock)then
       method%nBlock=nBlock
    else
       method%nBlock=2000
    endif
  end subroutine set_KSP_nBlock

  subroutine set_KSP_nvecBlock(method,nvecBlock)
    implicit none
    type(KSP) :: method
    integer nvecBlock
    method%nvecBlock=nvecBlock
  end subroutine set_KSP_nvecBlock

  subroutine set_KSP_addL(method,addL)
    implicit none
    type(KSP) :: method
    integer addL
    method%addL=addL
  end subroutine set_KSP_addL

  subroutine set_KSP_iflagAS(method,iflagAS)
    implicit none
    type(KSP) :: method
    integer iflagAS
    method%iflagAS=iflagAS
  end subroutine set_KSP_iflagAS

  subroutine set_KSP_iovlflag(method,iovlflag)
    implicit none
    type(KSP) :: method
    integer iovlflag
    method%iovlflag=iovlflag
  end subroutine set_KSP_iovlflag

  subroutine set_KSP_CAGMRES_sstep(method,CAGMRES_sstep)
    implicit none
    type(KSP) :: method
    integer CAGMRES_sstep
    method%CAGMRES_sstep=CAGMRES_sstep
  end subroutine set_KSP_CAGMRES_sstep

  subroutine set_KSP_CAGMRES_tstep(method,CAGMRES_tstep)
    implicit none
    type(KSP) :: method
    integer CAGMRES_tstep
    method%CAGMRES_tstep=CAGMRES_tstep
  end subroutine set_KSP_CAGMRES_tstep

  subroutine set_KSP_CAGMRES_basis(method,CAGMRES_basis)
    implicit none
    type(KSP) :: method
    integer CAGMRES_basis
    method%CAGMRES_basis=CAGMRES_basis
  end subroutine set_KSP_CAGMRES_basis

  subroutine set_KSP_CAGMRES_QRflag(method,CAGMRES_QRflag)
    implicit none
    type(KSP) :: method
    integer CAGMRES_QRflag
    method%CAGMRES_QRflag=CAGMRES_QRflag
  end subroutine set_KSP_CAGMRES_QRflag

  subroutine set_KSP_CAARNOLDI_sstep(method,CAARNOLDI_sstep)
    implicit none
    type(KSP) :: method
    integer CAARNOLDI_sstep
    method%CAARNOLDI_sstep=CAARNOLDI_sstep
  end subroutine set_KSP_CAARNOLDI_sstep

  subroutine set_KSP_CAARNOLDI_tstep(method,CAARNOLDI_tstep)
    implicit none
    type(KSP) :: method
    integer CAARNOLDI_tstep
    method%CAARNOLDI_tstep=CAARNOLDI_tstep
  end subroutine set_KSP_CAARNOLDI_tstep

  subroutine set_KSP_CAARNOLDI_basis(method,CAARNOLDI_basis)
    implicit none
    type(KSP) :: method
    integer CAARNOLDI_basis
    method%CAARNOLDI_basis=CAARNOLDI_basis
  end subroutine set_KSP_CAARNOLDI_basis

  subroutine set_KSP_CAARNOLDI_QRflag(method,CAARNOLDI_QRflag)
    implicit none
    type(KSP) :: method
    integer CAARNOLDI_QRflag
    method%CAARNOLDI_QRflag=CAARNOLDI_QRflag
  end subroutine set_KSP_CAARNOLDI_QRflag

  subroutine set_KSP_CALANCZOS_sstep(method,CALANCZOS_sstep)
    implicit none
    type(KSP) :: method
    integer CALANCZOS_sstep
    method%CALANCZOS_sstep=CALANCZOS_sstep
  end subroutine set_KSP_CALANCZOS_sstep

  subroutine set_KSP_CALANCZOS_tstep(method,CALANCZOS_tstep)
    implicit none
    type(KSP) :: method
    integer CALANCZOS_tstep
    method%CALANCZOS_tstep=CALANCZOS_tstep
  end subroutine set_KSP_CALANCZOS_tstep

  subroutine set_KSP_CALANCZOS_basis(method,CALANCZOS_basis)
    implicit none
    type(KSP) :: method
    integer CALANCZOS_basis
    method%CALANCZOS_basis=CALANCZOS_basis
  end subroutine set_KSP_CALANCZOS_basis

  subroutine set_KSP_CALANCZOS_QRflag(method,CALANCZOS_QRflag)
    implicit none
    type(KSP) :: method
    integer CALANCZOS_QRflag
    method%CALANCZOS_QRflag=CALANCZOS_QRflag
  end subroutine set_KSP_CALANCZOS_QRflag

  subroutine set_KSP_CBCG_kstep(method,CBCG_kstep)
    implicit none
    type(KSP) :: method
    integer CBCG_kstep
    method%CBCG_kstep=CBCG_kstep
  end subroutine set_KSP_CBCG_kstep

  subroutine set_KSP_power_method_itrmax(method,power_method_itrmax)
    implicit none
    type(KSP) :: method
    integer power_method_itrmax
    method%power_method_itrmax=power_method_itrmax
  end subroutine set_KSP_power_method_itrmax

  subroutine set_KSP_CBCG_Eigenflag(method,CBCG_Eigenflag)
    implicit none
    type(KSP) :: method
    integer CBCG_Eigenflag
    method%CBCG_Eigenflag=CBCG_Eigenflag
  end subroutine set_KSP_CBCG_Eigenflag

  subroutine set_KSP_GMRES_m (method,GMRES_m)
    implicit none
    type(KSP) :: method
    integer GMRES_m 
    method%GMRES_m=GMRES_m
  end subroutine set_KSP_GMRES_m

  subroutine set_KSP_GMRES_GSflag (method,GMRES_GSflag)
    implicit none
    type(KSP) :: method
    integer GMRES_GSflag 
    method%GMRES_GSflag=GMRES_GSflag
  end subroutine set_KSP_GMRES_GSflag

  subroutine set_KSP_ipreflag(method,ipreflag)
    implicit none
    type(KSP) :: method
    integer ipreflag
    method%ipreflag=ipreflag
  end subroutine set_KSP_ipreflag

  subroutine set_KSP_ILU_method(method,ILU_method)
    implicit none
    type(KSP) :: method
    integer ILU_method
    method%ILU_method=ILU_method
  end subroutine set_KSP_ILU_method

  subroutine set_KSP_independ_nvec(method,independ_nvec)
    implicit none
    type(KSP) :: method
    integer independ_nvec
    if(independ_nvec<1)then
       method%independ_nvec=1
    else
       method%independ_nvec=independ_nvec
    endif
  end subroutine set_KSP_independ_nvec

  subroutine set_KSP_precon_thblock(method,precon_thblock)
    implicit none
    type(KSP) :: method
    integer precon_thblock
    method%precon_thblock=precon_thblock
  end subroutine set_KSP_precon_thblock

  subroutine set_KSP_itrmax(method,itrmax)
    implicit none
    type(KSP) :: method
    integer itrmax
    method%itrmax=itrmax
  end subroutine set_KSP_itrmax

  subroutine set_KSP_rtolmax(method,rtolmax)
    implicit none
    type(KSP) :: method
    real*8  rtolmax
    method%rtolmax=rtolmax
  end subroutine set_KSP_rtolmax

  subroutine set_KSP_abstolmax(method,abstolmax)
    implicit none
    type(KSP) :: method
    real*8  abstolmax
    method%abstolmax=abstolmax
  end subroutine set_KSP_abstolmax

end module krylov
