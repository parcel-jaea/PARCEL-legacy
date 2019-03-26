module PARCEL_wrapper
  !  use krylov, only: KSP
  use krylov
  use krylov_kernel
  implicit none
  type(KSP),allocatable :: PARCEL_KSP(:)
  ! contains
end module PARCEL_wrapper

subroutine parcel_ksp_default_setting_wrapper(ino)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer ino
  call PARCEL_KSP_Default_Setting(PARCEL_KSP(ino))
end subroutine parcel_ksp_default_setting_wrapper

! 開始番号は0始まり
subroutine set_parcel_wrapper(imax)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer imax
  integer i
  allocate(PARCEL_KSP(0:imax-1))

  do i=1,imax
     call parcel_ksp_default_setting_wrapper(i-1)
  enddo
end subroutine set_parcel_wrapper

subroutine free_parcel_wrapper()
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer imax
  integer i
  imax=size(PARCEL_KSP)
  if(allocated(PARCEL_KSP))then
     do i=1,imax
        call free_KSP(PARCEL_KSP(i-1))
     enddo
     deallocate(PARCEL_KSP)
  endif
end subroutine free_parcel_wrapper

subroutine set_ksp_crs_wrapper(&
     icomm, &
     vecx,vecb, &
     n,gn,nnz,istart, &
     crsA,crsRow_ptr,crsCol, &
     itrmax,ipreflag,ILU_method,addL, &
     ino)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
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
  integer ipreflag ! 前処理方法　
  integer ILU_method ! 不完全LU分解(D-ILU分解)の対角行列の決定方法
  integer addL ! 加法シュワルツのオーバラップ領域の大きさ
  integer ino
  call set_KSP_CRS(&
       icomm, &
       vecx,vecb, &
       n,gn,nnz,istart, &
       crsA,crsRow_ptr,crsCol, &
       itrmax,ipreflag,ILU_method,addL, &
       PARCEL_KSP(ino))
end subroutine set_ksp_crs_wrapper

subroutine set_ksp_dia_wrapper(&
     icomm, &
     vecx,vecb, &
     n,gn,istart, &
     diaA,offset,nnd,&
     itrmax,ipreflag,ILU_method,addL,&
     ino)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer icomm ! MPI コミュニケータ
  integer(kind=intType) n ! プロセスの受け持つベクトルのサイズ
  integer(kind=intType) gn ! ベクトル全体のサイズ
  integer(kind=intType) istart ! ベクトルの開始位置
  real*8 vecx(n) ! 解ベクトル
  real*8 vecb(n) ! 右辺ベクトル
  integer(kind=intType) nnd
  real*8 diaA(nnd*n)
  integer(kind=intType) offset(nnd)
  integer itrmax ! 最大反復回数
  integer ipreflag ! 前処理方法　
  integer ILU_method ! 不完全LU分解(D-ILU分解)の対角行列の決定方法
  integer addL ! 加法シュワルツのオーバラップ領域の大きさ
  integer ino
  call set_KSP_DIA(&
       icomm, &
       vecx,vecb, &
       n,gn,istart, &
       diaA,offset,nnd,&
       itrmax,ipreflag,ILU_method,addL,&
       PARCEL_KSP(ino))

end subroutine set_ksp_dia_wrapper

subroutine set_vecx_ksp_wrapper(n,x,ino)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer(kind=intType) n
  real*8 x(n)
  integer ino
  call set_vecx_KSP(n,x,PARCEL_KSP(ino))
end subroutine set_vecx_KSP_wrapper

subroutine set_vecb_ksp_wrapper(n,b,ino)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer(kind=intType) n
  real*8 b(n)
  integer ino
  call set_vecb_KSP(n,b,PARCEL_KSP(ino))
end subroutine set_vecb_ksp_wrapper

subroutine get_vecx_ksp_wrapper(n,x,ino)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer ino
  integer(kind=intType) n
  real*8 x(n)
  call get_vecx_KSP(n,x,PARCEL_KSP(ino))
end subroutine get_vecx_ksp_wrapper

subroutine get_reshistory_ksp_wrapper(niter,reshistory,ino)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer ino
  integer niter
  real*8 reshistory(PARCEL_KSP(ino)%niter)
  integer(kind=intType) niter_intType
  call get_reshistory_KSP(niter,reshistory,PARCEL_KSP(ino))

end subroutine get_reshistory_ksp_wrapper

subroutine reset_crs_mat_wrapper(&
     itrmax, &
     ipreflag, &
     ILU_method, &
     addL,&
     crsA,&
     ino)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer itrmax
  integer ipreflag
  integer ILU_method
  integer addL
  real*8  crsA(PARCEL_KSP(ino)%prm%nnz) ! 行列値
  integer ino

  call reset_CRS_Mat(&
       itrmax, &
       ipreflag, &
       ILU_method, &
       addL,&
       crsA,&
       PARCEL_KSP(ino))

end subroutine reset_crs_mat_wrapper

subroutine reset_dia_mat_wrapper(&
     itrmax, &
     ipreflag, &
     ILU_method, &
     addL,&
     diaMat,&
     ino)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer itrmax
  integer ipreflag
  integer ILU_method
  integer addL
  real*8 diaMat(PARCEL_KSP(ino)%prm%nnd * PARCEL_KSP(ino)%prm%n)
  integer ino

  call reset_DIA_Mat(&
       itrmax, &
       ipreflag, &
       ILU_method, &
       addL,&
       diaMat,&
       PARCEL_KSP(ino))

end subroutine reset_dia_mat_wrapper

subroutine set_ksp_solver_wrapper(ino,solver)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer ino,solver
  call set_KSP_solver(PARCEL_KSP(ino),solver)
end subroutine set_ksp_solver_wrapper

subroutine set_ksp_abstolmax_wrapper(ino,abstolmax)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer ino
  real*8  abstolmax
  call set_KSP_abstolmax(PARCEL_KSP(ino),abstolmax)
end subroutine set_ksp_abstolmax_wrapper

subroutine set_ksp_rtolmax_wrapper(ino,rtolmax)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer ino
  real*8  rtolmax
  call set_KSP_rtolmax(PARCEL_KSP(ino),rtolmax)
end subroutine set_KSP_rtolmax_wrapper

subroutine set_ksp_iflagas_wrapper(ino,iflagAS)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer ino
  integer iflagAS
  call set_KSP_iflagas(PARCEL_KSP(ino),iflagAS)
end subroutine set_ksp_iflagas_wrapper

subroutine set_ksp_iovlflag_wrapper(ino,iovlflag)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer ino
  integer iovlflag
  call set_KSP_iovlflag(PARCEL_KSP(ino),iovlflag)
end subroutine set_ksp_iovlflag_wrapper

subroutine set_ksp_gmres_m_wrapper(ino,GMRES_m)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer ino
  integer GMRES_m
  call set_KSP_GMRES_m(PARCEL_KSP(ino),GMRES_m)
end subroutine set_ksp_gmres_m_wrapper

subroutine set_ksp_gmres_gsflag_wrapper(ino,GMRES_GSflag)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer ino
  integer GMRES_GSflag
  call set_KSP_GMRES_GSflag(PARCEL_KSP(ino),GMRES_GSflag)
end subroutine set_ksp_gmres_gsflag_wrapper

subroutine set_ksp_cagmres_sstep_wrapper(ino,CAGMRES_sstep)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer ino
  integer CAGMRES_sstep
  call set_KSP_CAGMRES_sstep(PARCEL_KSP(ino),CAGMRES_sstep)
end subroutine set_ksp_cagmres_sstep_wrapper

subroutine set_ksp_cagmres_tstep_wrapper(ino,CAGMRES_tstep)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer ino
  integer CAGMRES_tstep
  call set_KSP_CAGMRES_tstep(PARCEL_KSP(ino),CAGMRES_tstep)
end subroutine set_ksp_cagmres_tstep_wrapper

subroutine set_ksp_cagmres_basis_wrapper(ino,CAGMRES_basis)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer ino,CAGMRES_basis
  call set_KSP_CAGMRES_basis(PARCEL_KSP(ino),CAGMRES_basis)
end subroutine set_ksp_cagmres_basis_wrapper

subroutine set_ksp_cagmres_qrflag_wrapper(ino,CAGMRES_QRflag)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer ino
  integer CAGMRES_QRflag
  call set_KSP_CAGMRES_QRflag(PARCEL_KSP(ino),CAGMRES_QRflag)
end subroutine set_ksp_cagmres_qrflag_wrapper

subroutine set_ksp_cbcg_kstep_wrapper(ino,CBCG_kstep)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer ino
  integer CBCG_kstep
  call set_KSP_CBCG_kstep(PARCEL_KSP(ino),CBCG_kstep)
end subroutine set_ksp_cbcg_kstep_wrapper

subroutine set_ksp_cbcg_eigenflag_wrapper(ino,CBCG_Eigenflag)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer ino
  integer CBCG_Eigenflag
  call set_KSP_CBCG_Eigenflag(PARCEL_KSP(ino),CBCG_Eigenflag)
end subroutine set_ksp_cbcg_eigenflag_wrapper

subroutine set_ksp_caarnoldi_sstep_wrapper(ino,CAARNOLDI_sstep)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer ino
  integer CAARNOLDI_sstep
  call set_KSP_CAARNOLDI_sstep(PARCEL_KSP(ino),CAARNOLDI_sstep)
end subroutine set_ksp_caarnoldi_sstep_wrapper

subroutine set_ksp_caarnoldi_tstep_wrapper(ino,CAARNOLDI_tstep)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer ino
  integer CAARNOLDI_tstep
  call set_KSP_CAARNOLDI_tstep(PARCEL_KSP(ino),CAARNOLDI_tstep)
end subroutine set_ksp_caarnoldi_tstep_wrapper

subroutine set_ksp_caarnoldi_basis_wrapper(ino,CAARNOLDI_basis)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer ino
  integer CAARNOLDI_basis
  call set_KSP_CAARNOLDI_basis(PARCEL_KSP(ino),CAARNOLDI_basis)
end subroutine set_ksp_caarnoldi_basis_wrapper

subroutine set_ksp_caarnoldi_qrflag_wrapper(ino,CAARNOLDI_QRflag)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer ino
  integer CAARNOLDI_QRflag
  call set_KSP_CAARNOLDI_QRflag(PARCEL_KSP(ino),CAARNOLDI_QRflag)
end subroutine set_ksp_caarnoldi_qrflag_wrapper

subroutine set_KSP_power_method_itrmax_wrapper(ino,power_method_itrmax)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer ino
  integer power_method_itrmax
  call set_KSP_power_method_itrmax(PARCEL_KSP(ino),power_method_itrmax)
end subroutine set_KSP_power_method_itrmax_wrapper

subroutine start_ksp_solver_wrapper(ino)
  use krylov
  use krylov_kernel
  use PARCEL_wrapper
  implicit none
  integer ino
  call START_KSP_SOLVER(PARCEL_KSP(ino))
end subroutine start_ksp_solver_wrapper

