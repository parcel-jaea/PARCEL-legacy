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
  end type ddm_param

  type ddmmatrix_wrap
     integer mat_type ! 1:CRS 2:DIA 3:DDM 
     integer ndim
     integer nrow

     ! !parameters for CRS
     ! integer, allocatable :: row_ptr(:)
     ! integer, allocatable :: col_ind(:)
     ! double precision, allocatable :: val(:)

     !parameters for DDM
     integer nnd_dia
     integer, allocatable :: ioff_dia(:)
     double precision, allocatable :: val_dia(:)
     integer, allocatable :: ioff_grids(:,:)
     type(ddm_param) :: prm
  end type ddmmatrix_wrap

  ! 加法シュワルツ前処理に特有のパラメータ
  !
  ! 配列の6成分は それぞれ -x, +x, -y, +y, -z, +z 方向の隣接領域に対応する
  !
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
  end type AS_geom

end module ddmMat_wrap

module param_AS_type_ddm
  implicit none

  ! 加法シュヴァルツ前処理方法
  integer,parameter :: KSP_ADDITIVESCHWARZ_TYPE_BASIC=1       ! BASIC
  integer,parameter :: KSP_ADDITIVESCHWARZ_TYPE_RESTRICT=2    ! RESTRICT
  integer,parameter :: KSP_ADDITIVESCHWARZ_TYPE_INTERPOLATE=3 ! INTERPOLATE
  integer,parameter :: KSP_ADDITIVESCHWARZ_TYPE_NONE=4        ! NONE

  ! DDM前処理方向
  ! (x->y->z :: 1->2->3)
  integer,parameter :: KSP_DDM_DIVIDE_THREADS_DIRECT_X=1 ! X方向 
  integer,parameter :: KSP_DDM_DIVIDE_THREADS_DIRECT_Y=2 ! Y方向
  integer,parameter :: KSP_DDM_DIVIDE_THREADS_DIRECT_Z=3 ! Z方向

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
  end type precon_AS_DDM
  type precon_AS_wrap
     integer max_stencil
     integer max_threads
     integer AS_METHOD_TYPE
     type(precon_AS_DDM),allocatable :: AS_DDM(:)
  end type precon_AS_wrap

  ! contains

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
     integer ncpu
     integer myrank
     integer icomm !MPI communicator
     double precision, pointer :: diag(:)
     integer, pointer :: oloop !outer loop
     double precision tol !tolerance
     integer ipc !precondition flag
     integer nres !number of residual history
     double precision,pointer :: res(:) !residual history
     integer iret !result status
     integer commHide         ! 0: no hiding, 1: one time, 2: each time
     integer AS_len           ! additive Schwarz overlap length
     integer AS_diag_neighbor ! 1 => 斜め隣りの領域との通信を行う
     integer AS_mean
     integer AS_n
     type(AS_geom), pointer :: AS_g
     double precision, pointer :: AS_val_dia(:)
     integer,pointer :: AS_ioff_dia(:)
     type(ddm_param), pointer :: AS_prm
     double precision,pointer ::  AS_v1(:), AS_v2(:)

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
