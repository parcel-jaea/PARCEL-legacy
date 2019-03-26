module PARCEL_define
  implicit none
  ! 疎行列格納形式
  integer,parameter :: MATRIX_FORMAT_CRS=1  ! 圧縮行  格納形式
  integer,parameter :: MATRIX_FORMAT_DIA=2  ! 圧縮対角格納形式
  integer,parameter :: MATRIX_FORMAT_DDM=3  ! DDM形式

  ! 倍精度,4倍精度フラグ
  integer,parameter :: KSP_DOUBLE_PRECISION=1        !  倍精度
  integer,parameter :: KSP_DOUBLE_DOUBLE_PRECISION=2 ! 4倍精度

  ! クリロフ部分空間法
  integer,parameter :: KSP_SOLVER_CG=1      ! CG法
  integer,parameter :: KSP_SOLVER_BiCG=2    ! BiCG法
  integer,parameter :: KSP_SOLVER_GMRES=3   ! GMRES法
  integer,parameter :: KSP_SOLVER_CAGMRES=4 ! CAGMRES法
  integer,parameter :: KSP_SOLVER_CBCG=5    ! CBCG法
  integer,parameter :: KSP_CHECK_QR=-1      ! QR分解の直交性確認

  ! 疎行列ベクトル積の通信・計算オーバーラップ
  integer,parameter :: KSP_SPMV_OVL_NONE=0 ! 通信・計算オーバーラップ無し
  integer,parameter :: KSP_SPMV_OVL_CORE=1 ! ローカルな計算と通信をオーバーラップ
  integer,parameter :: KSP_SPMV_OVL_RANK=2 ! ランクごとの計算と通信をオーバーラップ

  ! クリロフ部分空間法の前処理方法
  integer,parameter :: KSP_PRECONDITION_NONE=0            ! 前処理無し
  integer,parameter :: KSP_PRECONDITION_POINTJACOBI=1     ! ポイントヤコビ前処理
  integer,parameter :: KSP_PRECONDITION_BLOCKJACOBI=2     ! ブロックヤコビ前処理
  integer,parameter :: KSP_PRECONDITION_ADDITIVESCHWARZ=3 ! 加法シュヴァルツ前処理

  ! 不完全LU分解(D-ILU分解)の対角行列の決定方法
  integer,parameter :: KSP_PRECONDITION_ILU  =0  ! ILU(0)
  integer,parameter :: KSP_PRECONDITION_ILUD =1  ! D-ILU 対角成分が一致
  integer,parameter :: KSP_PRECONDITION_MLUD =2  ! D-ILU 行内の要素和が一致

  ! 加法シュヴァルツ前処理方法
  integer,parameter :: KSP_ADDITIVESCHWARZ_TYPE_BASIC=1       ! BASIC
  integer,parameter :: KSP_ADDITIVESCHWARZ_TYPE_RESTRICT=2    ! RESTRICT
  integer,parameter :: KSP_ADDITIVESCHWARZ_TYPE_INTERPOLATE=3 ! INTERPOLATE
  integer,parameter :: KSP_ADDITIVESCHWARZ_TYPE_NONE=4        ! NONE

  ! QR分解方法
  integer,parameter :: KSP_QRDECOMPOSITION_MGS=1         ! 修正グラムシュミット
  integer,parameter :: KSP_QRDECOMPOSITION_CGS=2         ! 古典グラムシュミット
  integer,parameter :: KSP_QRDECOMPOSITION_TSQR=3        ! TSQR
  integer,parameter :: KSP_QRDECOMPOSITION_CHOLESKYQR=4  ! コレスキーQR
  integer,parameter :: KSP_QRDECOMPOSITION_CHOLESKYQR2=5 ! コレスキーQR2

  ! QR分解の直交行列計算有無
  integer,parameter :: KSP_QRDECOMPOSITION_FALSE_ORTHOGONAL_MATRIX=0 ! 直交行列の計算無
  integer,parameter :: KSP_QRDECOMPOSITION_TRUE_ORTHOGONAL_MATRIX=1  ! 直交行列の計算有

  ! GMRES(m)のグラムシュミット直交化方法
  integer,parameter :: KSP_GMRES_GS_TYPE_MGS=KSP_QRDECOMPOSITION_MGS ! 修正グラムシュミット
  integer,parameter :: KSP_GMRES_GS_TYPE_CGS=KSP_QRDECOMPOSITION_CGS ! 古典グラムシュミット

  ! CAGMRES(s,t)のQR分解方法
  integer,parameter :: KSP_CAGMRES_QR_TYPE_MGS         = KSP_QRDECOMPOSITION_MGS         ! 修正グラムシュミット
  integer,parameter :: KSP_CAGMRES_QR_TYPE_CGS         = KSP_QRDECOMPOSITION_CGS         ! 古典グラムシュミット
  integer,parameter :: KSP_CAGMRES_QR_TYPE_TSQR        = KSP_QRDECOMPOSITION_TSQR        ! TSQR
  integer,parameter :: KSP_CAGMRES_QR_TYPE_CHOLESKYQR  = KSP_QRDECOMPOSITION_CHOLESKYQR  ! コレスキーQR 
  integer,parameter :: KSP_CAGMRES_QR_TYPE_CHOLESKYQR2 = KSP_QRDECOMPOSITION_CHOLESKYQR2 ! コレスキーQR2

  ! CAARNOLDI(s,t)のQR分解方法
  integer,parameter :: KSP_CAARNOLDI_QR_TYPE_MGS         = KSP_QRDECOMPOSITION_MGS         ! 修正グラムシュミット
  integer,parameter :: KSP_CAARNOLDI_QR_TYPE_CGS         = KSP_QRDECOMPOSITION_CGS         ! 古典グラムシュミット
  integer,parameter :: KSP_CAARNOLDI_QR_TYPE_TSQR        = KSP_QRDECOMPOSITION_TSQR        ! TSQR
  integer,parameter :: KSP_CAARNOLDI_QR_TYPE_CHOLESKYQR  = KSP_QRDECOMPOSITION_CHOLESKYQR  ! コレスキーQR 
  integer,parameter :: KSP_CAARNOLDI_QR_TYPE_CHOLESKYQR2 = KSP_QRDECOMPOSITION_CHOLESKYQR2 ! コレスキーQR2

  ! 固有値計算方法
  integer,parameter :: KSP_EIGENVALUE_TYPE_POWER=1
  integer,parameter :: KSP_EIGENVALUE_TYPE_CAARNOLDI=2

  ! MPKの基底変換行列
  integer,parameter :: KSP_MPK_MONOMIAL_BASIS=0 ! 単基底
  integer,parameter :: KSP_MPK_NEWTON_BASIS  =1 ! ニュートン基底

  ! PARCEL リターンコード
  integer,parameter :: PARCEL_RETURNCODE_NORMAL=0       ! 正常終了
  integer,parameter :: PARCEL_RETURNCODE_MAXITERATION=1 ! 未収束
  ! 
  integer,parameter :: PARCEL_RESET_INHERIT=-1 ! リセット時の設定を引き継ぐ

! todo :: CG法の場合に残差ベクトルは収束したが真の残差が収束していない場合  integer,parameter :: PARCEL_RETURNCODE_MAXITERATION=2 ! 残差ベクトルは収束はしたが|| b -Ax ||/||b||の残差は収束条件よりも大きい
end module PARCEL_define

module PARCEL_integer
  use mpi
  implicit none
  ! PARCEL内部で利用する整数型のバイト数を定義
  ! MPI,LAPACKのサブルーチンについては4バイト整数も利用する
  ! -----
  ! * Fortran 2003コンパイラに対応している環境では
  !     " use iso_fortran_env + integer(int64) "
  !   を利用したほうが適切だが、JAEA ICE Xのgfortranのバージョンでは
  !   利用できないためinteger,parameter :: intTypeを定義している.

#ifdef PARCEL_INT8
! #if 1
! 8 バイト整数モード
  integer, parameter :: intType    = 8
  integer, parameter :: intMPIType = MPI_INTEGER8
#else
! 4 バイト整数モード
  integer, parameter :: intType    = 4
  integer, parameter :: intMPIType = MPI_INTEGER4
#endif

end module PARCEL_integer

