#include "ddfunc.h"
!!|
!!| AS_util_set_size
!!|
module AS_util_ddm

contains

  subroutine AS_util_set_size( myrank, ndist, prm, AS_len, &
        AS_diag_neighbor, AS_ndist, AS_prm, AS_g)
    !  use param_AS_type_ddm
    use ddmMat_wrap
    use param_AS_type_ddm
    implicit none
    !  include "ddm.h"
    integer     , intent(in ) :: myrank   ! 自身のプロセス番号
    integer(kind=intType)        , intent(in ) :: ndist    ! 
    type(ddm_param), intent(in ) :: prm      ! 
    integer        , intent(in ) :: AS_len   ! 加法シュワルツ法の重なり領域の幅
    integer        , intent(in ) :: AS_diag_neighbor ! 斜め隣の領域の値を使うなら1、使わないなら0
    integer(kind=intType), intent(inout) :: AS_ndist ! 加法シュワルツ前処理で使用する領域のグリッド数
    type(ddm_param), intent(inout) :: AS_prm   ! 加法シュワルツ前処理用に領域を変形したもの
    type(AS_geom)  , intent(inout) :: AS_g     ! 加法シュワルツ前処理に必要な変数を予め計算しておくためのもの
    !
    integer(kind=intType) i, j, k, kx,ky,kz, kxmax,kymax,kzmax
    integer xpos(6), ypos(6), zpos(6)

    !!
    !! AS_len /= 0 の時は加法シュワルツ前処理
    !!

    !
    ! 通信相手の番号 AS_g%rank_* の取得。 相手がいない場合は -1
    !
    xpos(1) =-1;  ypos(1) = 0;  zpos(1) = 0  ! -x
    xpos(2) = 1;  ypos(2) = 0;  zpos(2) = 0  ! +x
    xpos(3) = 0;  ypos(3) =-1;  zpos(3) = 0  ! -y
    xpos(4) = 0;  ypos(4) = 1;  zpos(4) = 0  ! +y
    xpos(5) = 0;  ypos(5) = 0;  zpos(5) =-1  ! -z
    xpos(6) = 0;  ypos(6) = 0;  zpos(6) = 1  ! +z
    AS_g%rank_s(:) = -1
    AS_g%rank_r(:) = -1
    do i=1, 6
       k = 1;
       kxmax = 0;  kymax = 0;  kzmax = 0;
       if ( prm%npe_z > 1 )  kzmax = 1;
       if ( prm%npe_y > 1 )  kymax = 1; 
       if ( prm%npe_x > 1 )  kxmax = 1; 
       Do kz=-kzmax, kzmax
          Do ky=-kymax, kymax
             Do kx=-kxmax, kxmax
                if ( kx*kx +ky*ky +kz*kz == 0 ) cycle
                if (prm%neighbor_ranks(k) > -1) then
                   if ( kx==xpos(i) .and. ky==ypos(i) .and. kz==zpos(i) ) then
                      AS_g%rank_s(i) = prm%neighbor_ranks(k)
                   end if
                   if ( kx==-xpos(i) .and. ky==-ypos(i) .and. kz==-zpos(i) ) then
                      AS_g%rank_r(i) = prm%neighbor_ranks(k)
                   end if
                end if
                k = k +1
             End Do
          End Do
       End Do
    end do

    !
    ! 重なり幅 sleeve_** の取得
    !
    AS_g%sleeve_mx = 0;  AS_g%sleeve_my = 0;  AS_g%sleeve_mz = 0
    AS_g%sleeve_px = 0;  AS_g%sleeve_py = 0;  AS_g%sleeve_pz = 0

    if( AS_g%rank_s(1) /= -1 )  AS_g%sleeve_mx = AS_len
    if( AS_g%rank_s(2) /= -1 )  AS_g%sleeve_px = AS_len
    if( AS_g%rank_s(3) /= -1 )  AS_g%sleeve_my = AS_len
    if( AS_g%rank_s(4) /= -1 )  AS_g%sleeve_py = AS_len
    if( AS_g%rank_s(5) /= -1 )  AS_g%sleeve_mz = AS_len
    if( AS_g%rank_s(6) /= -1 )  AS_g%sleeve_pz = AS_len
    !
    ! AS_ndist, AS_prm の取得
    !
    AS_prm = prm
    !
    AS_prm%nx = prm%nx + AS_g%sleeve_mx + AS_g%sleeve_px
    AS_prm%ny = prm%ny + AS_g%sleeve_my + AS_g%sleeve_py
    AS_prm%nz = prm%nz + AS_g%sleeve_mz + AS_g%sleeve_pz

    !
    AS_prm%nxy = AS_prm%nx * AS_prm%ny
    AS_prm%nxyz = AS_prm%nxy* AS_prm%nz
    AS_prm%stmx=0
    AS_prm%stmy=0
    AS_prm%stmz=0
    AS_prm%mx = AS_prm%nx
    AS_prm%my = AS_prm%ny
    AS_prm%mz = AS_prm%nz
    AS_prm%mxy = AS_prm%nxy
    !
    AS_ndist = AS_prm%nx * AS_prm%ny * AS_prm%nz

    !
    ! 通信領域の原点と形状 AS_g%x0s, AS_g%xlen_s, AS_g%grids_* とその仲間の取得
    !
    AS_g%x0s = AS_g%sleeve_mx;  AS_g%y0s = AS_g%sleeve_my;  AS_g%z0s = AS_g%sleeve_mz
    AS_g%x0s(2) = AS_prm%nx - 2 * AS_g%sleeve_px
    AS_g%y0s(4) = AS_prm%ny - 2 * AS_g%sleeve_py
    AS_g%z0s(6) = AS_prm%nz - 2 * AS_g%sleeve_pz
    !
    AS_g%x0r = AS_g%sleeve_mx;  AS_g%y0r = AS_g%sleeve_my;  AS_g%z0r = AS_g%sleeve_mz
    AS_g%x0r(1) = AS_prm%nx - AS_g%sleeve_px;  AS_g%x0r(2) = 0
    AS_g%y0r(3) = AS_prm%ny - AS_g%sleeve_py;  AS_g%y0r(4) = 0
    AS_g%z0r(5) = AS_prm%nz - AS_g%sleeve_pz;  AS_g%z0r(6) = 0
    !
    AS_g%xlen_s = prm%nx;  AS_g%ylen_s = prm%ny;  AS_g%zlen_s = prm%nz
    AS_g%xlen_s(1) = AS_g%sleeve_mx;  AS_g%xlen_s(2) = AS_g%sleeve_px
    AS_g%ylen_s(3) = AS_g%sleeve_my;  AS_g%ylen_s(4) = AS_g%sleeve_py
    AS_g%zlen_s(5) = AS_g%sleeve_mz;  AS_g%zlen_s(6) = AS_g%sleeve_pz
    !
    if( AS_diag_neighbor /= 0 ) then
       do k=3, 6
          AS_g%x0s(k) = 0;  AS_g%x0r(k) = 0;  AS_g%xlen_s(k) = AS_prm%nx
       end do
       do k=5, 6
          AS_g%y0s(k) = 0;  AS_g%y0r(k) = 0;  AS_g%ylen_s(k) = AS_prm%ny
       end do
    end if
    !
    do i=1, 6
       if ( AS_prm%ncomp_grids == 1 ) then
          AS_g%grids_s(i) = AS_g%xlen_s(i)
       else if ( AS_prm%ncomp_grids == 2 ) then
          AS_g%grids_s(i) = AS_g%xlen_s(i) * AS_g%ylen_s(i)
       else
          AS_g%grids_s(i) = AS_g%xlen_s(i) * AS_g%ylen_s(i) * AS_g%zlen_s(i)
       end if
    end do
    do i=1, 6
       if ( mod(i,2) == 0 ) j = i - 1
       if ( mod(i,2) == 1 ) j = i + 1
       AS_g%grids_r(i) = AS_g%grids_s(j)
       AS_g% xlen_r(i) = AS_g% xlen_s(j)
       AS_g% ylen_r(i) = AS_g% ylen_s(j)
       AS_g% zlen_r(i) = AS_g% zlen_s(j)
    end do

  end subroutine AS_util_set_size

  !!|
  !!| AS_util_set_ioff_A
  !!|
  !!| AS_ioff_A を設定する
  !!|
  subroutine AS_util_set_ioff_A ( nnd, ncomp_grids, AS_ioff_grids, AS_prm, AS_ioff_A )
    use ddmMat_wrap
    implicit none
    !  include "ddm.h"
    integer        , intent(in ) :: nnd            ! 係数行列の非ゼロ対角線の数   
    integer        , intent(in ) :: ncomp_grids    ! 領域の次元     
    integer        , intent(in ) :: AS_ioff_grids(nnd,ncomp_grids) ! 
    type(ddm_param), intent(in ) :: AS_prm         !  
    integer        , intent(inout) :: AS_ioff_A(nnd) ! [i] = i番目の副対角部が対角部からどれだけずれているか
    !
    integer i

    do i=1, nnd

       if( ncomp_grids == 1) then
          AS_ioff_A(i) = AS_ioff_grids(i,1)
       else if ( ncomp_grids == 2) then
          AS_ioff_A(i) = AS_ioff_grids(i,1) + AS_prm%mx * AS_ioff_grids(i,2)
       else if ( ncomp_grids == 3) then
          AS_ioff_A(i) = AS_ioff_grids(i,1) + AS_prm%mx * AS_ioff_grids(i,2) + AS_prm%mxy * AS_ioff_grids(i,3)
       end if

    end do
  end subroutine AS_util_set_ioff_A



  !!|
  !!| AS_util_set_A
  !!|
  !!| A を AS_A にコピーする
  !!|
  subroutine AS_util_set_A ( &
       comm, myrank, nnd, ndist, prm, A, AS_ioff_A, &
       AS_ndist, AS_len, AS_A )
    use ddmMat_wrap
    use mpi
    use PARCEL_define,only : KSP_DOUBLE_PRECISION,KSP_DOUBLE_DOUBLE_PRECISION
    implicit none
    !   include "ddm.h"
    !   include "mpif.h"
    integer        , intent(in ) :: comm           !
    integer        , intent(in ) :: myrank         !
    integer        , intent(in ) :: nnd            ! 係数行列の非ゼロ対角線の数   
    integer(kind=intType) , intent(in ) :: ndist          !      
    type(ddm_param), intent(in ) :: prm            !    
    real*8         , intent(in ) :: A(ndist*nnd)   ! 係数行列　    
    integer        , intent(in ) :: AS_ioff_A(nnd) ! 対角要素のオフセット
    integer(kind=intType), intent(in ) :: AS_ndist       ! （AS = 加法シュワルツ）
    integer        , intent(in ) :: AS_len         !
    real*8         , intent(inout) :: AS_A(AS_ndist*nnd) ! 
    integer PRECISION
    !
    integer(kind=intType) AS_i, i, ix, iy, iz, j, k, n
    integer xyz, pm, lenxy0, lenxyz0, info
    integer idiag, val
    type array
       real*8, allocatable :: v(:)
    end type array
    type(array) sbuf(6), rbuf(6) ! 6 = ±x, ±y, ±z   
    integer :: ireqs(6), ireqr(6)
    integer istatus (mpi_status_size)
    integer(kind=intType)         AS_ndist_dummy
    type(ddm_param) AS_prm
    type(AS_geom)   AS_g

    !!|
    !!| パラメータ設定
    !!|   斜め隣りの係数行列Aも取得するようにAS_gを作り直す
    !!|
    !    call AS_util_set_size( myrank, ndist, prm, AS_len, &
    !          .true., AS_ndist_dummy, AS_prm, AS_g )
    ! call AS_util_set_size( myrank, ndist, prm, AS_len, &
    !       0, AS_ndist_dummy, AS_prm, AS_g )
    call AS_util_set_size( myrank, ndist, prm, AS_len, &
          1, AS_ndist_dummy, AS_prm, AS_g )

    !!|
    !!| 通信が不要な部分のコピー
    !!|
    lenxy0  = prm%nx * prm%ny
    lenxyz0 = prm%nx * prm%ny * prm%nz
    !
    !$OMP PARALLEL 
    !
    !$OMP do private(i, AS_i,  ix, iy, iz, j)
    !
    do k=0, lenxyz0 - 1
       ix = mod(k, prm%nx)
       iy = mod(k/prm%nx,  prm%ny)
       iz = k/lenxy0
       i = ( ix +prm%stmx ) +prm%mx *( iy +prm%stmy ) +prm%mxy *( iz +prm%stmz ) +1
       AS_i = (ix + AS_g%sleeve_mx) + AS_prm%mx * (iy + AS_g%sleeve_my) + AS_prm%mxy * (iz + AS_g%sleeve_mz) +1

       do j=0, nnd-1
          AS_A(AS_i + j * AS_ndist) = A(i + j * ndist)
       end do
    end do

    !$OMP end do

    !$OMP end PARALLEL


    !!|
    !!| 通信が必要な部分のコピー
    !!|   Aのデータの代わりに、上でコピーしたAS_Aのデータを通信する
    !!|
    !
    do i=1, 6
       allocate( sbuf(i)%v(nnd * AS_g%grids_s(i)) )
       allocate( rbuf(i)%v(nnd * AS_g%grids_r(i)) )
    end do
    !
    do xyz=0, 2  ! x, y, z
       !
       ! 同じ軸方向の領域のみ同時に通信できる
       !
       do pm=0, 1  ! +xyz, -xyz
          i = 2 * xyz + pm + 1
          !
          ! 通信データの設定
          !
          if ( AS_g%rank_s(i) /= -1 ) then
             !
             !$OMP PARALLEL do private(ix, iy, iz, k, n)
             !
             do j=0, AS_g%grids_s(i)-1
                ix = AS_g%x0s(i) + mod(j, AS_g%xlen_s(i))
                iy = AS_g%y0s(i) + mod(j/AS_g%xlen_s(i),  AS_g%ylen_s(i))
                iz = AS_g%z0s(i) + j / (AS_g%xlen_s(i) * AS_g%ylen_s(i))
                k = ix +AS_prm%mx * iy + AS_prm%mxy * iz + 1
                do n=0, nnd-1
                   sbuf(i)%v(j + 1 + n * AS_g%grids_s(i)) = AS_A(k + n * AS_ndist)
                end do
             end do
             !
             !$OMP end PARALLEL do
             !
          end if
          !
          ! 通信開始
          !
          if ( AS_g%rank_s(i) /= -1 ) then
             call MPI_Isend( sbuf(i)%v, size(sbuf(i)%v), mpi_double_precision, &
                  AS_g%rank_s(i), AS_g%rank_s(i), comm, ireqs(i), info )
          end if
          if ( AS_g%rank_r(i) /= -1 ) then
             call MPI_Irecv( rbuf(i)%v, size(rbuf(i)%v), mpi_double_precision, &
                  AS_g%rank_r(i), myrank, comm, ireqr(i), info )
          end if
       end do
       !
       ! 通信を完了させて、受信したデータをコピー
       !
       do pm=0, 1
          i = 2 * xyz + pm + 1
          if ( AS_g%rank_s(i) /= -1 ) call MPI_Wait( ireqs(i), istatus, info )
          if ( AS_g%rank_r(i) /= -1 ) call MPI_Wait( ireqr(i), istatus, info )

          !
          ! 受信したデータをコピー
          !
          if ( AS_g%rank_r(i) /= -1 ) then
             !$OMP PARALLEL do private(ix, iy, iz, k, n)
             do j=0, AS_g%grids_r(i)-1
                ix = AS_g%x0r(i) + mod(j, AS_g%xlen_r(i))
                iy = AS_g%y0r(i) + mod(j/AS_g%xlen_r(i),  AS_g%ylen_r(i))
                iz = AS_g%z0r(i) + j / (AS_g%xlen_r(i) * AS_g%ylen_r(i))
                k = ix + AS_prm%mx * iy + AS_prm%mxy * iz + 1
                do n=0, nnd-1
                   AS_A(k + n * AS_ndist) = rbuf(i)%v(j + 1 + n * AS_g%grids_r(i)) 
                end do
             end do
             !$OMP end PARALLEL do
          end if
       end do
    end do

    do i=1, 6
       deallocate( sbuf(i)%v )
       deallocate( rbuf(i)%v )
    end do

  end subroutine AS_util_set_A


  !!|
  !!| AS_util_copy_to_AS
  !!|
  !!| wv1 を AS_wv1 にコピーする
  !!|
  subroutine AS_util_copy_to_AS ( comm, myrank, ndist, prm, AS_ndist, AS_prm, &
       AS_diag_neighbor, AS_g, wv1, AS_wv1,iflag_AS)
    use param_AS_type_ddm
    use ddmMat_wrap
    use mpi
    implicit none
    !   include "ddm.h"
    !   include "mpif.h"
    integer        , intent(in ) :: comm   !
    integer        , intent(in ) :: myrank !
    integer(kind=intType)        , intent(in ) :: ndist  ! 
    type(ddm_param), intent(in ) :: prm    ! 
    integer(kind=intType) , intent(in ) :: AS_ndist         ! 領域のグリッド数（AS = 加法シュワルツ）
    type(ddm_param), intent(in ) :: AS_prm           ! 領域形状
    integer        , intent(in ) :: AS_diag_neighbor ! 斜め隣の領域の値を使うなら1、使わないなら0
    type(AS_geom)  , intent(in ) :: AS_g
    real*8         , intent(in ) :: wv1(ndist)       ! 
    real*8         , intent(inout) :: AS_wv1(AS_ndist) ! 前処理行列の作用対象となるベクトル
    integer iflag_AS ! 加法シュワルツの補正方法  !
    integer(kind=intType) k, AS_i, i, ix, iy, iz, j
    integer lenx0, leny0, lenz0, lenxy0, lenxyz0
    integer xyz, pm, info
    type array
       real*8, allocatable :: v(:)
    end type array
    type(array) sbuf(6), rbuf(6) ! 6 = ±x, ±y, ±z   
    integer     ireqs(6), ireqr(6)
    integer     istatus (mpi_status_size)

    ! 通信が不要な部分のコピー
    !$OMP PARALLEL
    !$OMP do private(i)
    do i=1, AS_ndist
       AS_wv1(i) = 0.0d0
    end do
    !$OMP end do
    !$OMP do private( i, AS_i, ix, iy, iz) collapse(2)
    do iz = 0,prm%nz -1
       do iy = 0,prm%ny -1 
          do ix = 0,prm%nx -1             
             AS_i = ( ix +AS_g%sleeve_mx ) +AS_prm%mx *( iy +AS_g%sleeve_my )  +AS_prm%mxy *( iz +AS_g%sleeve_mz ) +1
             i    = ( ix +prm%stmx ) +prm%mx *( iy +prm%stmy )  +prm%mxy *( iz +prm%stmz ) +1
             AS_wv1(AS_i) = wv1(i)
          enddo
       enddo
    enddo
    !$OMP end PARALLEL

    if(&
         (iflag_AS.eq.KSP_ADDITIVESCHWARZ_TYPE_INTERPOLATE) &
         .or. &
         (iflag_AS.eq.KSP_ADDITIVESCHWARZ_TYPE_NONE) &
         )then
       return
    endif

    !!| 隣接プロセスからデータを通信により取得する
    !!| wv1のデータの代わりに、上でコピーしたAS_wv1のデータを通信する
    do i=1, 6
       allocate( sbuf(i)%v(AS_g%grids_s(i)) )
       allocate( rbuf(i)%v(AS_g%grids_r(i)) )
    end do

    ! 斜め隣りのプロセスの値を使用する場合
    do xyz=0, 2  ! x, y, z
       ! 同じ軸方向の領域のみ同時に通信できる
       do pm=0, 1  ! +xyz, -xyz
          i = 2 * xyz + pm + 1
          ! 通信データの設定
          if ( AS_g%rank_s(i) /= -1 ) then
             !$OMP PARALLEL do private(ix, iy, iz, k)
             do j=0, size(sbuf(i)%v)-1
                ix = AS_g%x0s(i) + mod(j , AS_g%xlen_s(i)                 )
                iy = AS_g%y0s(i) + mod(j / AS_g%xlen_s(i) , AS_g%ylen_s(i))
                iz = AS_g%z0s(i) + j / (AS_g%xlen_s(i) * AS_g%ylen_s(i))
                k = ix + AS_prm%mx * iy + AS_prm%mxy * iz + 1
                sbuf(i)%v(j+1) = AS_wv1(k)
             end do
             !$OMP end PARALLEL do
          end if

          ! 通信開始
          if ( AS_g%rank_s(i) /= -1 ) then
             call MPI_Isend( sbuf(i)%v, size(sbuf(i)%v), mpi_double_precision, &
                  AS_g%rank_s(i), AS_g%rank_s(i), comm, ireqs(i), info )
          end if
          if ( AS_g%rank_r(i) /= -1 ) then
             call MPI_Irecv( rbuf(i)%v, size(rbuf(i)%v), mpi_double_precision, &
                  AS_g%rank_r(i), myrank, comm, ireqr(i), info )
          end if
       end do

       ! 通信を完了させて、受信したデータをコピー
       do pm=0, 1
          i = 2 * xyz + pm + 1
          if ( AS_g%rank_s(i) /= -1 ) call MPI_Wait( ireqs(i), istatus, info )
          if ( AS_g%rank_r(i) /= -1 ) call MPI_Wait( ireqr(i), istatus, info )
          ! 受信したデータをコピー
          if ( AS_g%rank_r(i) /= -1 ) then
             !$OMP PARALLEL do private(ix, iy, iz, k)
             do j=0, size(rbuf(i)%v)-1
                ix = AS_g%x0r(i) + mod(j, AS_g%xlen_r(i))
                iy = AS_g%y0r(i) + mod(j/AS_g%xlen_r(i),  AS_g%ylen_r(i))
                iz = AS_g%z0r(i) + j / (AS_g%xlen_r(i) * AS_g%ylen_r(i))
                k = ix + AS_prm%mx * iy + AS_prm%mxy * iz + 1
                AS_wv1(k) = rbuf(i)%v(j+1)
             end do
             !$OMP end PARALLEL do
          end if
       end do
    end do

    do i=1, 6
       deallocate( sbuf(i)%v )
       deallocate( rbuf(i)%v )
    end do

  end subroutine AS_util_copy_to_AS


  !!|
  !!| AS_util_copy_from_AS
  !!|
  !!| AS_wv4 を wv4 にコピーする
  !!|
  subroutine AS_util_copy_from_AS ( comm, myrank, ndist, prm, AS_ndist, AS_prm, &
       AS_diag_neighbor, AS_mean, AS_g, wv4, AS_wv4,&
       PRECISION_CALC, &
       wv4_lo, &
       AS_wv4_lo &
       )

    use ddmMat_wrap
    use mpi
    use PARCEL_define
    implicit none
    !   include "ddm.h"
    !   include "mpif.h"
    integer        , intent(in ) :: comm      !
    integer        , intent(in ) :: myrank    !
    integer(kind=intType)        , intent(in ) :: ndist     ! 
    type(ddm_param), intent(in ) :: prm       ! 
    integer(kind=intType)        , intent(in ) :: AS_ndist  ! 領域のグリッド数（AS = 加法シュワルツ）
    type(ddm_param), intent(in ) :: AS_prm    ! 領域形状
    integer        , intent(in ) :: AS_diag_neighbor ! 斜め隣の領域の値を使うなら1、使わないなら0
    integer        , intent(in ) :: AS_mean          ! 隣の領域との平均をとるなら1、とらないなら0
    type(AS_geom)  , intent(in ) :: AS_g             !
    real*8         , intent(inout) :: wv4(ndist)       ! 
    real*8         , intent(inout) :: AS_wv4(AS_ndist) ! 前処理行列の作用対象となるベクトル
    !
    integer(kind=intType) AS_i, i, ix, iy, iz, j, k
    integer size, x0, y0, z0, xlen, ylen, zlen, info
    real*8  coeff
    integer mpx, mpy, mpz, xyz, pm, nsize
    type array
       real*8, allocatable :: v(:)
    end type array
    type(array) sbuf(6), rbuf(6) ! 6 = ±x, ±y, ±z   
    integer :: ireqs(6), ireqr(6)
    integer istatus (mpi_status_size)

    integer,optional :: PRECISION_CALC

    real*8,optional , intent(inout) :: wv4_lo(ndist)       ! 
    real*8,optional , intent(inout) :: AS_wv4_lo(AS_ndist) ! 前処理行列の作用対象とな   

    real*8 ddadd_var1
    type(array) sbuf_lo(6), rbuf_lo(6) ! 6 = ±x, ±y, ±z   
    integer :: ireqs_lo(6), ireqr_lo(6)

    integer PRECISION
    PRECISION = KSP_DOUBLE_PRECISION
    if(  present(PRECISION_CALC)  )then
       PRECISION = PRECISION_CALC
    endif

    !!| 
    !!| 平均をとらない場合は通信せずにコピーだけ
    !!|
    if ( AS_mean .eq. 0 ) goto 100


    !!| 
    !!| 加算する場合 
    !!|
    ! 隣接プロセスからデータを通信により取得する
    !   AS_util_copy_to_AS の場合に対して send / receive の関係が逆になることに注意
    !   取得したデータをAS_wv4に加算すたのち、適当な数で割ることで平均を計算する
    do i=1, 6
       allocate( sbuf(i)%v(AS_g%grids_s(i)) ) ! こっちが受信用
       allocate( rbuf(i)%v(AS_g%grids_r(i)) ) ! こっちが送信用

       allocate( sbuf_lo(i)%v(AS_g%grids_s(i)) ) ! こっちが受信用
       allocate( rbuf_lo(i)%v(AS_g%grids_r(i)) ) ! こっちが送信用
    end do

    !
    ! 斜め隣りのプロセスの値を使用する場合
    !
    do xyz=2, 0, -1  ! x, y, z
       !
       ! 同じ軸方向の領域のみ同時に通信できる
       !
       do pm=0, 1  ! +xyz, -xyz
          i = 2 * xyz + pm + 1
          !
          ! 通信データの設定
          !
          if ( AS_g%rank_r(i) /= -1 ) then
#if 1
             !$OMP PARALLEL do private(ix, iy, iz, k)
             !             do j=0, size(rbuf(i)%v)-1
             do j=0, AS_g%grids_r(i)-1
                ix = AS_g%x0r(i) + mod(j, AS_g%xlen_r(i))
                iy = AS_g%y0r(i) + mod(j/AS_g%xlen_r(i), AS_g%ylen_r(i))
                iz = AS_g%z0r(i) + j / (AS_g%xlen_r(i) * AS_g%ylen_r(i))
                k = ix + AS_prm%mx * iy + AS_prm%mxy * iz + 1
                rbuf(i)%v(j+1) = AS_wv4(k)
                if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                   rbuf_lo(i)%v(j+1) = AS_wv4_lo(k)
                endif
             end do
             !$OMP end PARALLEL do
#endif
          end if

          ! 通信開始
          if ( AS_g%rank_r(i) /= -1 ) then
                call MPI_Isend( rbuf(i)%v, size(rbuf(i)%v), mpi_double_precision, &
                     & AS_g%rank_r(i), AS_g%rank_r(i), comm, ireqr(i), info )
                if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                   call MPI_Isend( rbuf_lo(i)%v, size(rbuf_lo(i)%v), mpi_double_precision, &
                        & AS_g%rank_r(i), AS_g%rank_r(i) + 1 , comm, ireqr_lo(i), info )
                endif
             end if
             if ( AS_g%rank_s(i) /= -1 ) then
                call MPI_Irecv( sbuf(i)%v, size(sbuf(i)%v), mpi_double_precision, &
                     & AS_g%rank_s(i), myrank, comm, ireqs(i), info )
                if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                   call MPI_Irecv( sbuf_lo(i)%v, size(sbuf_lo(i)%v), mpi_double_precision, &
                        & AS_g%rank_s(i), myrank + 1, comm, ireqs_lo(i), info )
                endif
             end if
          end do

          ! 通信を完了させて、受信したデータを加算
          do pm=0, 1
             i = 2 * xyz + pm + 1
             if ( AS_g%rank_s(i) /= -1 ) call MPI_Wait( ireqs(i), istatus, info )
             if ( AS_g%rank_r(i) /= -1 ) call MPI_Wait( ireqr(i), istatus, info )
             if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                if ( AS_g%rank_s(i) /= -1 ) call MPI_Wait( ireqs_lo(i), istatus, info )
                if ( AS_g%rank_r(i) /= -1 ) call MPI_Wait( ireqr_lo(i), istatus, info )
             endif
             ! 受信したデータを加算
             if ( AS_g%rank_s(i) /= -1 ) then
#if 1
                !$OMP PARALLEL do private(ix, iy, iz, k,&
                !$OMP ddadd_var1 &
                !$OMP )
                !                do j=0, size(sbuf(i)%v)-1
                do j=0, AS_g%grids_s(i)-1
                   ix = AS_g%x0s(i) + mod(j, AS_g%xlen_s(i))
                   iy = AS_g%y0s(i) + mod(j/AS_g%xlen_s(i),  AS_g%ylen_s(i))
                   iz = AS_g%z0s(i) + j / (AS_g%xlen_s(i) * AS_g%ylen_s(i))
                   k = ix + AS_prm%mx * iy + AS_prm%mxy * iz + 1
                   if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                      ! AS_wv4(k) = AS_wv4(k) + sbuf(i)%v(j+1)
                      ddadd_calc1(AS_wv4(k),AS_wv4_lo(k),sbuf(i)%v(j+1),sbuf_lo(i)%v(j+1),AS_wv4(k),AS_wv4_lo(k))
                   else
                      AS_wv4(k) = AS_wv4(k) + sbuf(i)%v(j+1)
                   endif
                end do
                !$OMP end PARALLEL do
#endif                
             end if
          end do
       end do

       !
       do i=1, 6
          deallocate( sbuf(i)%v )
          deallocate( rbuf(i)%v )
          deallocate( sbuf_lo(i)%v )
          deallocate( rbuf_lo(i)%v )
       end do

100 continue

    !!| 
    !!| AS_wv4 を wv4 にコピー
    !!|
    !
!    PRECISION = KSP_DOUBLE_PRECISION

    !$OMP PARALLEL
#if 1
    !$OMP do private(i)
    do i=1, ndist
       wv4(i) = 0.0
    end do
    !$OMP end do
#endif
    if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
       !$OMP do private(i)
       do i=1, ndist
          wv4_lo(i) = 0.0
       end do
       !$OMP end do
    end if

#if 1
 
#if 1
    !$OMP do private(i, AS_i, ix, iy, iz) collapse(2)
    do iz = 0,prm%nz -1
       do iy = 0,prm%ny -1 
          do ix = 0,prm%nx -1             
             AS_i = ( ix +AS_g%sleeve_mx ) +AS_prm%mx *( iy +AS_g%sleeve_my )  +AS_prm%mxy *( iz +AS_g%sleeve_mz ) +1
             i = ( ix +prm%stmx ) +prm%mx *( iy +prm%stmy )  +prm%mxy *( iz +prm%stmz ) +1
             wv4(i) = AS_wv4(AS_i)
             if( PRECISION .eq. KSP_DOUBLE_DOUBLE_PRECISION)then
                wv4_lo(i) = AS_wv4_lo(AS_i)
             endif
          end do
       end do
    end do
    !$OMP end do
#endif
    
#else
             
    !$OMP do private(i, AS_i, ix, iy, iz)
    do k=0, prm%nxyz - 1
       ix = mod(k, prm%nx)
       iy = mod(k/prm%nx,  prm%ny)
       iz = k/prm%nxy
       AS_i = ( ix +AS_g%sleeve_mx ) +AS_prm%mx *( iy +AS_g%sleeve_my )  +AS_prm%mxy *( iz +AS_g%sleeve_mz ) +1
       i = ( ix +prm%stmx ) +prm%mx *( iy +prm%stmy )  +prm%mxy *( iz +prm%stmz ) +1
       wv4(i) = AS_wv4(AS_i)
    end do
    !$OMP end do
#endif    
    !$OMP end PARALLEL

  end subroutine AS_util_copy_from_AS


end module AS_util_ddm
