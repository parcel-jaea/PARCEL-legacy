#define AS_Debug 0
! #define CUDA_AWARE 1

module precon_DILU
  use c_wrapper
  implicit none

  interface
     subroutine solve_DILU_DDM_FP64_GPU(&
          numGPUblock, &
          numGPUthread, &
          A, &
          Dinv, offset_Dinv, &
          r,    offset_r, &
          s,    offset_s, &
          svec, offset_svec, &
          ioff_dia, &
          offset_xyz, &
          nxe, &
          nye, &
          nze, &
          nx, &
          ny, &
          nz, &
          stmx, &
          stmy, &
          stmz, &
          mx, &
          my, &
          mz, &
          nxblock, &
          nyblock, &
          nzblock, &
          nx_local, &
          ny_local, &
          nz_local, &
          mx_local, &
          my_local, &
          mz_local, &
          nnd       )&
          bind(c,name = "cuda_solve_DILU_DDM")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       type(devicemem) A
       type(devicemem) Dinv
       integer(kind=intType),value :: offset_Dinv
       type(devicemem) r
       integer(kind=intType),value :: offset_r
       type(devicemem) s
       integer(kind=intType),value :: offset_s
       type(devicemem) svec
       integer(kind=intType),value :: offset_svec

       type(devicemem) ioff_dia
       type(devicemem) offset_xyz
       type(devicemem) nxe
       type(devicemem) nye
       type(devicemem) nze
       integer,value ::  nx
       integer,value ::  ny
       integer,value ::  nz
       integer,value ::  stmx
       integer,value ::  stmy
       integer,value ::  stmz
       integer,value ::  mx
       integer,value ::  my
       integer,value ::  mz
       integer,value ::  nxblock
       integer,value ::  nyblock
       integer,value ::  nzblock
       integer,value ::  nx_local
       integer,value ::  ny_local
       integer,value ::  nz_local
       integer,value ::  mx_local
       integer,value ::  my_local
       integer,value ::  mz_local
       integer,value ::  nnd
     end subroutine solve_DILU_DDM_FP64_GPU

     subroutine solve_DILU_AS_DDM_FP64_GPU(&
          AS_v1, offset_AS_v1, &
          s    , offset_s, &
          rvec, offset_rvec, &
          svec, offset_svec, &
          svec_AS, offset_svec_AS, &
          AS_v2, offset_AS_v2, &
          AS_val_dia_init, &
          Dinv, offset_Dinv, &
          ioff_dia, &
          offset_xyz, &
          nnd,  &
          nxs, nxe, &
          nys, nye, &
          nzs, nze, &
          nxs0, nxe0, &
          nys0, nye0, &
          nzs0, nze0, &
          addL, &
          nblock, &
          nxblock,nyblock,nzblock, &
          stmx,stmy,stmz, &
          mx, mxy, &
          m_local, &
          mx_local,my_local, &
          m_AS, &
          mx_AS ,mxy_AS, &
          sleeve_mx,sleeve_my,sleeve_mz &
          )&
          bind(c,name = "cuda_solve_DILU_AS_DDM")
       import :: intType
       import :: devicemem
       implicit none

       type(devicemem) AS_v1
       integer(kind=intType),value :: offset_AS_v1
       type(devicemem) s
       integer(kind=intType),value :: offset_s
       type(devicemem) rvec
       integer(kind=intType),value :: offset_rvec
       type(devicemem) svec
       integer(kind=intType),value :: offset_svec
       type(devicemem) svec_AS
       integer(kind=intType),value :: offset_svec_AS
       type(devicemem) AS_v2
       integer(kind=intType),value :: offset_AS_v2
       type(devicemem) AS_val_dia_init 
       type(devicemem) Dinv            
       integer(kind=intType),value :: offset_Dinv
       type(devicemem) ioff_dia
       type(devicemem) offset_xyz
       integer,value :: nnd 
       type(devicemem) nxs
       type(devicemem) nxe

       type(devicemem) nys
       type(devicemem) nye

       type(devicemem) nzs
       type(devicemem) nze

       type(devicemem) nxs0
       type(devicemem) nxe0

       type(devicemem) nys0
       type(devicemem) nye0

       type(devicemem) nzs0
       type(devicemem) nze0

       integer,value :: addL

       integer,value :: nblock

       integer,value :: nxblock
       integer,value :: nyblock
       integer,value :: nzblock

       integer,value :: stmx
       integer,value :: stmy
       integer,value :: stmz

       integer,value :: mx
       integer,value :: mxy

       integer,value :: m_local
       integer,value :: mx_local
       integer,value :: my_local

       integer,value :: m_AS

       integer,value :: mx_AS 
       integer,value :: mxy_AS

       integer,value :: sleeve_mx
       integer,value :: sleeve_my
       integer,value :: sleeve_mz

     end subroutine solve_DILU_AS_DDM_FP64_GPU

     subroutine AS_util_copy_only_FP64_GPU(&
			 wv1, offset_wv1, &
			 AS_wv1, offset_AS_wv1, &
			 nx,ny,nz, &
			 mx,mxy, &
			 stmx,stmy,stmz, &
			 AS_ndist, &
			 mx_AS,mxy_AS, &
			 sleeve_mx,sleeve_my,sleeve_mz &
          ) &
          bind(c,name = "cuda_AS_util_copy_only")
       import :: intType
       import :: devicemem
       implicit none

       type(devicemem) wv1
       integer(kind=intType),value :: offset_wv1
       type(devicemem) AS_wv1
       integer(kind=intType),value :: offset_AS_wv1
       integer,value :: nx
       integer,value :: ny
       integer,value :: nz
       integer,value :: mx
       integer,value :: mxy
       integer,value :: stmx
       integer,value :: stmy
       integer,value :: stmz
       integer(kind=intType),value :: AS_ndist
       integer,value :: mx_AS
       integer,value :: mxy_AS
       integer,value :: sleeve_mx
       integer,value :: sleeve_my
       integer,value :: sleeve_mz

     end subroutine AS_util_copy_only_FP64_GPU

     subroutine util_AS_copy_only_FP64_GPU(&
			 wv1, offset_wv1, &
			 AS_wv1, offset_AS_wv1, &
			 nx,ny,nz, &
			 mx,mxy, &
			 stmx,stmy,stmz, &
			 ndist, &
			 mx_AS,mxy_AS, &
			 sleeve_mx,sleeve_my,sleeve_mz &
          ) &
          bind(c,name = "cuda_util_AS_copy_only")
       import :: intType
       import :: devicemem
       implicit none

       type(devicemem) wv1
       integer(kind=intType),value :: offset_wv1
       type(devicemem) AS_wv1
       integer(kind=intType),value :: offset_AS_wv1
       integer,value :: nx
       integer,value :: ny
       integer,value :: nz
       integer,value :: mx
       integer,value :: mxy
       integer,value :: stmx
       integer,value :: stmy
       integer,value :: stmz
       integer(kind=intType),value :: ndist
       integer,value :: mx_AS
       integer,value :: mxy_AS
       integer,value :: sleeve_mx
       integer,value :: sleeve_my
       integer,value :: sleeve_mz

     end subroutine Util_AS_copy_only_FP64_GPU

     subroutine AS_util_set_sbuf_FP64_GPU(&
			v_, offset_v, &
			AS_wv1_, offset_AS_wv1, &
			nloop, &
			mx_AS,mxy_AS, &
			xlen_s,ylen_s, &
			x0s,y0s,z0s  &
          ) &
          bind(c,name = "cuda_AS_set_sbuf")
       import :: intType
       import :: devicemem
       implicit none
       type(devicemem) v_
       integer(kind=intType),value ::  offset_v
       type(devicemem) AS_wv1_
       integer(kind=intType),value ::  offset_AS_wv1
       integer,value :: nloop
       integer,value :: mx_AS
       integer,value :: mxy_AS
       integer,value :: xlen_s
       integer,value :: ylen_s
       integer,value :: x0s
       integer,value :: y0s
       integer,value :: z0s 

     end subroutine AS_util_set_sbuf_FP64_GPU

     subroutine AS_util_set_rbuf_FP64_GPU(&
			v_, offset_v, &
			AS_wv1_, offset_AS_wv1, &
			nloop, &
			mx_AS,mxy_AS, &
			xlen_r,ylen_r, &
			x0r,y0r,z0r,&
                        alpha &
          ) &
          bind(c,name = "cuda_AS_set_rbuf")
       import :: intType
       import :: devicemem
       implicit none
       type(devicemem) v_
       integer(kind=intType),value ::  offset_v
       type(devicemem) AS_wv1_
       integer(kind=intType),value ::  offset_AS_wv1
       integer,value :: nloop
       integer,value :: mx_AS
       integer,value :: mxy_AS
       integer,value :: xlen_r
       integer,value :: ylen_r
       integer,value :: x0r
       integer,value :: y0r
       integer,value :: z0r 
       real*8,value :: alpha

     end subroutine AS_util_set_rbuf_FP64_GPU

     subroutine cuda_get_ddmmatrix_neumanndiag_FP64_GPU( &
          nx, &
          ny, &
          nz, &
          stmx, &
          stmy, &
          stmz, &
          mx, &
          mxy, &
          m, &
          nnd, &
          diag_offset, &
          A, &
          Dinv_,  offset_Dinv, &
          ioff_A &
          )bind(c,name = "cuda_get_ddmmatrix_neumanndiag")

       import :: intType
       import :: devicemem
       implicit none
#if 0
       integer(kind=intType),value ::  nx
       integer(kind=intType),value ::  ny
       integer(kind=intType),value ::  nz
       integer(kind=intType),value ::  stmx
       integer(kind=intType),value ::  stmy
       integer(kind=intType),value ::  stmz
       integer(kind=intType),value ::  mx
       integer(kind=intType),value ::  mxy
       integer(kind=intType),value ::  m
       integer(kind=intType),value ::  nnd
#else
       integer,value ::  nx
       integer,value ::  ny
       integer,value ::  nz
       integer,value ::  stmx
       integer,value ::  stmy
       integer,value ::  stmz
       integer,value ::  mx
       integer,value ::  mxy
       integer(kind=intType),value ::  m
       integer,value ::  nnd
#endif
       integer(kind=intType),value ::  diag_offset
       type(devicemem) A
       type(devicemem) Dinv_ 
       integer(kind=intType),value ::  offset_Dinv
       type(devicemem) ioff_A
       
     end  subroutine cuda_get_ddmmatrix_neumanndiag_FP64_GPU


  end interface

contains

  ! -------------------------------------
  ! -------------------------------------

  subroutine init_AS_DDM( &
       method,&
       xdiv,ydiv,zdiv &
       )
    use krylov
    implicit none
    type(KSP) :: method
    integer xdiv,ydiv,zdiv 

    integer nnd
    integer nx,ny,nz
    integer stmx, stmy, stmz

    integer lattice_dimention
    integer nx_local,ny_local,nz_local
    integer m_local
    integer mx_local,my_local,mz_local
    integer nblock
    integer nxblock,nyblock,nzblock

    integer i
    integer nxs,nxe
    integer nys,nye
    integer nzs,nze

    integer j
    integer jx,jy,jz
    integer jxblock,jyblock,jzblock

    xdiv = max(xdiv,1)
    ydiv = max(ydiv,1)
    zdiv = max(zdiv,1)

    lattice_dimention = method%ddmMat%prm%ncomp_grids

    nnd = method%ddmMat%nnd
    nx = method%ddmMat%prm%nx
    ny = method%ddmMat%prm%ny
    nz = method%ddmMat%prm%nz

    stmx = method%ddmMat%prm%stmx
    stmy = method%ddmMat%prm%stmy
    stmz = method%ddmMat%prm%stmz

    ! ------------
    if(mod(nx,xdiv).ne.0)then
       nx_local = nx/xdiv + 1
    else
       nx_local = nx/xdiv
    endif
    allocate(method%ddmMat%DDMMat_DILU%nxs(xdiv))
    allocate(method%ddmMat%DDMMat_DILU%nxe(xdiv))
    allocate(method%ddmMat%DDMMat_DILU%nxs0(xdiv))
    allocate(method%ddmMat%DDMMat_DILU%nxe0(xdiv))

    do i = 1,xdiv
       nxs = 1 + nx_local*(i-1)
       nxe = nx_local*i
       method%ddmMat%DDMMat_DILU%nxe(i)  = min(nx,nxe)

       method%ddmMat%DDMMat_DILU%nxs0(i) = nxs
       method%ddmMat%DDMMat_DILU%nxe0(i) = min(nx,nxe)
    enddo
    ! mx_local = nx_local + stmx*2

    ! ------------
    if(mod(ny,ydiv).ne.0)then
       ny_local = ny/ydiv + 1
    else
       ny_local = ny/ydiv
    endif
    allocate(method%ddmMat%DDMMat_DILU%nys(ydiv))
    allocate(method%ddmMat%DDMMat_DILU%nye(ydiv))
    allocate(method%ddmMat%DDMMat_DILU%nys0(ydiv))
    allocate(method%ddmMat%DDMMat_DILU%nye0(ydiv))

    do i = 1,ydiv
       nys = 1 + ny_local*(i-1)
       nye = ny_local*i
       method%ddmMat%DDMMat_DILU%nye(i) = min(ny,nye)

       method%ddmMat%DDMMat_DILU%nys0(i) = nys
       method%ddmMat%DDMMat_DILU%nye0(i) = min(ny,nye)
    enddo
    ! my_local = ny_local + stmy*2

    ! ------------
    if(mod(nz,zdiv).ne.0)then
       nz_local = nz/zdiv + 1
    else
       nz_local = nz/zdiv
    endif
    allocate(method%ddmMat%DDMMat_DILU%nzs(zdiv))
    allocate(method%ddmMat%DDMMat_DILU%nze(zdiv))
    allocate(method%ddmMat%DDMMat_DILU%nzs0(zdiv))
    allocate(method%ddmMat%DDMMat_DILU%nze0(zdiv))
    do i = 1,zdiv
       nzs = 1 + nz_local*(i-1)
       nze = nz_local*i
       method%ddmMat%DDMMat_DILU%nze(i) = min(nz,nze)

       method%ddmMat%DDMMat_DILU%nzs0(i) = nzs
       method%ddmMat%DDMMat_DILU%nze0(i) = min(nz,nze)
    enddo
    ! mz_local = nz_local + stmz*2

    ! ------------
    ! nx_local,ny_local,nz_local
    ! 通信で持ってくる行列データの袖の長さ : 0なら計算領域の端
    ! method%ddmMat%AS_g%AS_g%sleeve_mx 
    ! method%ddmMat%AS_g%AS_g%sleeve_my 
    ! method%ddmMat%AS_g%AS_g%sleeve_mz 
    ! method%ddmMat%AS_g%AS_g%sleeve_px 
    ! method%ddmMat%AS_g%AS_g%sleeve_py 
    ! method%ddmMat%AS_g%AS_g%sleeve_pz 

#if 0
    write(*,*) "AS_g   : addL  = ",method%ddmMat%AS_g%addL
    write(*,*) "x : sleev = ", &
         method%ddmMat%AS_g%sleeve_mx,method%ddmMat%AS_g%sleeve_px
    write(*,*) "y : sleev = ", &
         method%ddmMat%AS_g%sleeve_my,method%ddmMat%AS_g%sleeve_py
    write(*,*) "z : sleev = ", &
         method%ddmMat%AS_g%sleeve_mz,method%ddmMat%AS_g%sleeve_pz
#endif

#if AS_Debug
    method%ddmMat%AS_g%addL = 0
#endif

    ! 領域の先頭は-sleeve_m[x,y,z]でずらす --> 端で領域外アクセスにならない?
    do i = 1,xdiv
       nxs = 1 + nx_local*(i-1)
       if( (method%ddmMat%AS_g%sleeve_mx.eq.0) .and. (i.eq.1) )then
          method%ddmMat%DDMMat_DILU%nxs(i) = nxs
       else
          method%ddmMat%DDMMat_DILU%nxs(i) = nxs - method%ddmMat%AS_g%addL
       endif

       if( (method%ddmMat%AS_g%sleeve_px.eq.0) .and. (i.eq.xdiv) )then
          method%ddmMat%DDMMat_DILU%nxe(i) = method%ddmMat%DDMMat_DILU%nxe(i)
       else
          method%ddmMat%DDMMat_DILU%nxe(i) = &
               method%ddmMat%DDMMat_DILU%nxe(i) + method%ddmMat%AS_g%addL
       endif
    enddo
    do i = 1,ydiv
       nys = 1 + ny_local*(i-1)
       if( (method%ddmMat%AS_g%sleeve_my.eq.0) .and. (i.eq.1) )then
          method%ddmMat%DDMMat_DILU%nys(i) = nys
       else
          method%ddmMat%DDMMat_DILU%nys(i) = nys - method%ddmMat%AS_g%addL
       endif

       if( (method%ddmMat%AS_g%sleeve_py.eq.0) .and. (i.eq.ydiv) )then
          method%ddmMat%DDMMat_DILU%nye(i) = method%ddmMat%DDMMat_DILU%nye(i)
       else
          method%ddmMat%DDMMat_DILU%nye(i) = &
               method%ddmMat%DDMMat_DILU%nye(i) + method%ddmMat%AS_g%addL
       endif
    enddo
    do i = 1,zdiv
       nzs = 1 + nz_local*(i-1)
       if( (method%ddmMat%AS_g%sleeve_mz.eq.0) .and. (i.eq.1) )then
          method%ddmMat%DDMMat_DILU%nzs(i) = nzs
       else
          method%ddmMat%DDMMat_DILU%nzs(i) = nzs - method%ddmMat%AS_g%addL
       endif

       if( (method%ddmMat%AS_g%sleeve_pz.eq.0) .and. (i.eq.zdiv) )then
          method%ddmMat%DDMMat_DILU%nze(i) = method%ddmMat%DDMMat_DILU%nze(i)
       else
          method%ddmMat%DDMMat_DILU%nze(i) = method%ddmMat%DDMMat_DILU%nze(i) + method%ddmMat%AS_g%addL
       endif
    enddo

    nx_local = nx_local + method%ddmMat%AS_g%addL*2
    ny_local = ny_local + method%ddmMat%AS_g%addL*2
    nz_local = nz_local + method%ddmMat%AS_g%addL*2
    mx_local = nx_local + stmx*2
    my_local = ny_local + stmy*2
    mz_local = nz_local + stmz*2
    ! ------------

    m_local = mx_local * my_local * mz_local 

    nxblock = xdiv
    nyblock = ydiv
    nzblock = zdiv

    nblock  = nxblock * nyblock * nzblock
    method%ddmMat%DDMMat_DILU%lattice_dimention = lattice_dimention
    method%ddmMat%DDMMat_DILU%nxblock = nxblock
    method%ddmMat%DDMMat_DILU%nyblock = nyblock
    method%ddmMat%DDMMat_DILU%nzblock = nzblock

    method%ddmMat%DDMMat_DILU%nx_local = nx_local
    method%ddmMat%DDMMat_DILU%ny_local = ny_local
    method%ddmMat%DDMMat_DILU%nz_local = nz_local

    method%ddmMat%DDMMat_DILU%mx_local = mx_local
    method%ddmMat%DDMMat_DILU%my_local = my_local
    method%ddmMat%DDMMat_DILU%mz_local = mz_local
    method%ddmMat%DDMMat_DILU%m_local  = m_local 

    allocate(method%ddmMat%DDMMat_DILU%Dinv(m_local * nblock))
    allocate(method%ddmMat%DDMMat_DILU%svec(m_local * nblock))
    allocate(method%ddmMat%DDMMat_DILU%rvec(m_local * nblock))
    allocate(method%ddmMat%DDMMat_DILU%svec_AS(m_local * nblock))

    !$OMP PARALLEL DO
    do i = 1,m_local * nblock
       method%ddmMat%DDMMat_DILU%Dinv(i) = 0.0d0
       method%ddmMat%DDMMat_DILU%svec(i) = 0.0d0
       method%ddmMat%DDMMat_DILU%rvec(i) = 0.0d0
       method%ddmMat%DDMMat_DILU%svec_AS(i) = 0.0d0
    enddo
    !$OMP END PARALLEL DO

    allocate(method%ddmMat%DDMMat_DILU%offset_xyz(nnd*3))
    method%ddmMat%DDMMat_DILU%offset_xyz = 0 

    select case (lattice_dimention)
    case(1)
       do i = 1,nnd
          method%ddmMat%DDMMat_DILU%offset_xyz(i+nnd*0) = method%ddmMat%mat%ioff_grids(i,1)
       enddo
    case(2)
       do i = 1,nnd
          method%ddmMat%DDMMat_DILU%offset_xyz(i+nnd*0) = method%ddmMat%mat%ioff_grids(i,1)
           method%ddmMat%DDMMat_DILU%offset_xyz(i+nnd*1) = method%ddmMat%mat%ioff_grids(i,2)
       enddo
    case(3)
       do i = 1,nnd
          method%ddmMat%DDMMat_DILU%offset_xyz(i+nnd*0) = method%ddmMat%mat%ioff_grids(i,1)
          method%ddmMat%DDMMat_DILU%offset_xyz(i+nnd*1) = method%ddmMat%mat%ioff_grids(i,2)
          method%ddmMat%DDMMat_DILU%offset_xyz(i+nnd*2) = method%ddmMat%mat%ioff_grids(i,3)
       enddo
    end select

    ! write(*,*)"nxblock  = ",nxblock
    ! write(*,*)"nyblock  = ",nyblock
    ! write(*,*)"nzblock  = ",nzblock

    ! write(*,*)"mx_local = ",mx_local
    ! write(*,*)"my_local = ",my_local
    ! write(*,*)"mz_local = ",mz_local

#if 0
    do jxblock = 1,size(method%ddmMat%DDMMat_DILU%nxe)
!       nxs = 1 + nx_local*(jxblock-1)
       nxs = method%ddmMat%DDMMat_DILU%nxs(jxblock)
       nxe = method%ddmMat%DDMMat_DILU%nxe(jxblock)
       write(*,*)method%commVec%myrank," jxloop : ",jxblock,nxs,nxe

    enddo
    do jyblock = 1,size(method%ddmMat%DDMMat_DILU%nye)
       nys = method%ddmMat%DDMMat_DILU%nys(jyblock)
       nye = method%ddmMat%DDMMat_DILU%nye(jyblock)
       write(*,*)method%commVec%myrank," jyloop : ",jyblock,nys,nye
    enddo
    do jzblock = 1,size(method%ddmMat%DDMMat_DILU%nze)
       nzs = method%ddmMat%DDMMat_DILU%nzs(jzblock)
       nze = method%ddmMat%DDMMat_DILU%nze(jzblock)
       write(*,*)method%commVec%myrank," jzloop : ",jzblock,nzs,nze
    enddo

    ! if(method%commVec%myrank.eq.0)then
    !    write(*,*)"DDM BlockJACOBI local size : nx,ny,nz = ", &
    !         nx_local,ny_local,nz_local
    ! endif
#endif

  end subroutine init_AS_DDM

  subroutine get_DILU_AS_DDM(method)
    use krylov
    implicit none

    type(KSP) :: method
    integer nx,ny,nz
    integer stmx, stmy, stmz
    integer m
    integer mx,my,mz
    integer mxy

    integer m_AS
    integer nx_AS,ny_AS,nz_AS
    integer mxy_AS
    integer mx_AS,my_AS,mz_AS

    integer lattice_dimention
    integer nblock
    integer nxblock,nyblock,nzblock
    integer nx_local,ny_local,nz_local
    integer m_local
    integer mx_local,my_local,mz_local

    integer nnd,offseti,offsetj,offsetii,offsetjj
    integer offsetj_AS,offseti_AS

    integer offsetx,offsety,offsetz
    
    integer j,jj
    integer jx,jy,jz
    integer jjx,jjy,jjz
    integer block_id
    integer jxblock,jyblock,jzblock
    integer ip,ipp,ip_diag
    integer j_AS
    integer nxs,nxe
    integer nys,nye
    integer nzs,nze
    real*8 Dinv 

    nx = method%ddmMat%prm%nx
    ny = method%ddmMat%prm%ny
    nz = method%ddmMat%prm%nz

    stmx = method%ddmMat%prm%stmx
    stmy = method%ddmMat%prm%stmy
    stmz = method%ddmMat%prm%stmz

    mx  = method%ddmMat%prm%mx
    my  = method%ddmMat%prm%my
    mz  = method%ddmMat%prm%mz
    mxy = method%ddmMat%prm%mxy
    m   =  mx * my * mz

    ! ----------
    ! m_AS : 拡張領域の長さが必要 --> ASデータは袖なしなのでそれを考慮
    nx_AS = method%ddmMat%AS_prm%nx
    ny_AS = method%ddmMat%AS_prm%ny
    nz_AS = method%ddmMat%AS_prm%nz

    mx_AS = nx_AS
    my_AS = ny_AS
    mz_AS = nz_AS
    mxy_AS = mx_AS * my_AS
    m_AS   = mx_AS * my_AS * mz_AS
    ! write(*,*)"m,m_AS = ",m,m_AS
    ! ----------

    nxblock = method%ddmMat%DDMMat_DILU%nxblock
    nyblock = method%ddmMat%DDMMat_DILU%nyblock
    nzblock = method%ddmMat%DDMMat_DILU%nzblock
    nblock  = nxblock*nyblock*nzblock;

    lattice_dimention = method%ddmMat%DDMMat_DILU%lattice_dimention 
    nx_local = method%ddmMat%DDMMat_DILU%nx_local
    ny_local = method%ddmMat%DDMMat_DILU%ny_local
    nz_local = method%ddmMat%DDMMat_DILU%nz_local

    mx_local = method%ddmMat%DDMMat_DILU%mx_local
    my_local = method%ddmMat%DDMMat_DILU%my_local
    mz_local = method%ddmMat%DDMMat_DILU%mz_local
    m_local  = method%ddmMat%DDMMat_DILU%m_local

    nnd = method%ddmMat%nnd

    ip_diag = -1
    do ip = 1,nnd
       if(method%ddmMat%mat%ioff_dia(ip).eq.0)then
          ip_diag = ip
       endif
    enddo

    ! do ipp = 1,nnd
    !    write(*,*)" AS : offset ",method%ddmMat%ASval%AS_ioff_dia_init(ipp)
    ! enddo

    ! do ip = 1,nnd
    !    write(*,*)"AS offset : ", &
    !         method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*0), &
    !         method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*1), &
    !         method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*2)
    ! enddo


    ! do jz = 1,1
    !    do jy = 1,1
    !       do jx = 1,nx_AS
    !          j_AS = jx + mx_AS*jy + mxy_AS*jz + 1
    !          method%ddmMat%ASval%AS_val_dia_init(j_AS + m_AS*(ip_diag-1))
    !       enddo
    !    enddo
    ! enddo

    !$OMP PARALLEL DO PRIVATE(&
    !$OMP jzblock,jyblock,jxblock, &
    !$OMP nzs,nze, &
    !$OMP nys,nye, &
    !$OMP nxs,nxe, &
    !$OMP block_id, &
    !$OMP j,jx,jy,jz, &
    !$OMP jj,jjx,jjy,jjz, &
    !$OMP ip,ipp, &
    !$OMP j_AS, &
    !$OMP offsetj,offsetjj, &
    !$OMP offseti,offsetii, &
    !$OMP offsetx,offsety,offsetz, &
    !$OMP offsetj_AS,offseti_AS   &
    !$OMP ) &
    !$OMP COLLAPSE(3)
    do jzblock=1,nzblock
       do jyblock=1,nyblock
          do jxblock=1,nxblock
             nzs = method%ddmMat%DDMMat_DILU%nzs(jzblock)
             nze = method%ddmMat%DDMMat_DILU%nze(jzblock)
             nys = method%ddmMat%DDMMat_DILU%nys(jyblock)
             nye = method%ddmMat%DDMMat_DILU%nye(jyblock)
             nxs = method%ddmMat%DDMMat_DILU%nxs(jxblock)
             nxe = method%ddmMat%DDMMat_DILU%nxe(jxblock)
             block_id = (jxblock-1) + (jyblock-1)*nxblock + (jzblock-1)*nxblock*nyblock 

             do jz = nzs-1,nze-1
                do jy = nys-1,nye-1
                   do jx = nxs-1,nxe-1
                      j = (jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1
#if 0
                      j_AS =  jx + jy * mx_AS + jz * mxy_AS + 1
#else
                      j_AS = &
                           (jx + method%ddmMat%AS_g%sleeve_mx) + &
                           (jy + method%ddmMat%AS_g%sleeve_my) * mx_AS  + &
                           (jz + method%ddmMat%AS_g%sleeve_mz) * mxy_AS + &
                           1
#endif
                      jjx = jx - (nxs-1)
                      jjy = jy - (nys-1)
                      jjz = jz - (nzs-1)
                      jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1                      
#if AS_Debug
                      method%ddmMat%DDMMat_DILU%Dinv(block_id + jj*nblock) = &
                           method%ddmMat%mat%val_dia(m*(ip_diag-1)+j)
#else
                      ! ASの行列データは袖なしで確保されている --> 座標が0以下になるのはおかしい?
                      ! jx,jy,jzに関して補正が必要?
                      method%ddmMat%DDMMat_DILU%Dinv(block_id + jj*nblock) = &
                           method%ddmMat%ASval%AS_val_dia_init(j_AS + m_AS*(ip_diag-1))
#endif
                   enddo
                enddo
             enddo

             do jz = nzs-1,nze-1
                do jy = nys-1,nye-1
                   do jx = nxs-1,nxe-1
                      j = ( jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1
#if 0
                      j_AS =  jx + jy * mx_AS + jz * mxy_AS + 1
#else
                      j_AS = &
                           (jx + method%ddmMat%AS_g%sleeve_mx) + &
                           (jy + method%ddmMat%AS_g%sleeve_my) * mx_AS  + &
                           (jz + method%ddmMat%AS_g%sleeve_mz) * mxy_AS + &
                           1
#endif

                      jjx = jx - (nxs-1)
                      jjy = jy - (nys-1)
                      jjz = jz - (nzs-1)
                      jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1

                      ! -----
                      ! D-ILU分解
                      ! A = AL + AD + AU ~= K
                      ! K = (D+AL)D^-1(D+AU)
                      ! 条件1 対角成分が一致     :: D(j,j) = A(j,j) - Sum[i<j]          A(j,i) D(i,i)^-1 A(i,j)
                      ! 条件2 行内の要素和が一致 :: D(j,j) = A(j,j) - Sum[i<j] Sum[k>i] A(j,i) D(i,i)^-1 A(i,k)
                      ! -----
                      do ip = 1,nnd
                         offsetj = method%ddmMat%mat%ioff_dia(ip)
                         offsetj_AS = &
                              method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*0) + &
                              method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*1)*mx_AS  + &
                              method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*2)*mxy_AS 

                         offsetx  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*0)
                         offsety  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*1)
                         offsetz  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*2)
                         offsetjj = offsetx + offsety*mx_local + offsetz*mx_local*my_local

                         ! A(j,i) : method%ddmMat%mat%val_dia(j+m*(ip-1))
                         do ipp = 1,nnd
                            offseti  = method%ddmMat%mat%ioff_dia(ipp)
                            offseti_AS = &
                              method%ddmMat%DDMMat_DILU%offset_xyz(ipp+nnd*0) + &
                              method%ddmMat%DDMMat_DILU%offset_xyz(ipp+nnd*1)*mx_AS  + &
                              method%ddmMat%DDMMat_DILU%offset_xyz(ipp+nnd*2)*mxy_AS 

                            offsetx  = method%ddmMat%DDMMat_DILU%offset_xyz(ipp+nnd*0)
                            offsety  = method%ddmMat%DDMMat_DILU%offset_xyz(ipp+nnd*1)
                            offsetz  = method%ddmMat%DDMMat_DILU%offset_xyz(ipp+nnd*2)
                            ! offsetii_AS ?
                            offsetii = offsetx + offsety*mx_local + offsetz*mx_local*my_local 
                            if( jj <= jj+offsetii)then
                               cycle
                            endif
                            ! D(j,j) = A(j,j) - Sum[i<j] A(j,i) D(i,i)^-1 A(i,j) 
#if AS_Debug
                            if( offseti .eq. (-offsetj) )then
                               method%ddmMat%DDMMat_DILU%Dinv(jj*nblock + block_id) = &
                                    method%ddmMat%DDMMat_DILU%Dinv( (jj           )*nblock + block_id) - &
                                    method%ddmMat%DDMMat_DILU%Dinv( (jj + offsetii)*nblock + block_id) * &
                                    method%ddmMat%mat%val_dia(j +           m*(ipp-1)  )               * &
                                    method%ddmMat%mat%val_dia(j + offseti + m*(ip-1)   )
                            endif
#else
                            if( offseti_AS .eq. (-offsetj_AS) )then
                               method%ddmMat%DDMMat_DILU%Dinv(jj*nblock + block_id) = &
                                    method%ddmMat%DDMMat_DILU%Dinv( (jj           )*nblock + block_id) - &
                                    method%ddmMat%DDMMat_DILU%Dinv( (jj + offsetii)*nblock + block_id) * &
                                    ! 拡張領域の行列データ
                                    method%ddmMat%ASval%AS_val_dia_init(j_AS +              m_AS*(ipp-1)) * &
                                    method%ddmMat%ASval%AS_val_dia_init(j_AS + offseti_AS + m_AS*(ip-1) )
                            endif
#endif
                         enddo
                      enddo
                      if(abs(method%ddmMat%DDMMat_DILU%Dinv(block_id + jj*nblock)).eq.0.0d0)then
                         Dinv = 0.0d0
                      else
                         Dinv = 1.0d0 / method%ddmMat%DDMMat_DILU%Dinv(block_id + jj*nblock)
                      endif
                      method%ddmMat%DDMMat_DILU%Dinv(block_id + jj*nblock) = Dinv

                   enddo
                enddo
             enddo

          enddo
       enddo
    end do

#if 0
    ! do jzblock = 1,1
    !    do jyblock = 1,1
    !       do jxblock = 1,1
    do jzblock=1,nzblock
       do jyblock=1,nyblock
          do jxblock=1,nxblock
             nzs = method%ddmMat%DDMMat_DILU%nzs(jzblock)
             nze = method%ddmMat%DDMMat_DILU%nze(jzblock)
             nys = method%ddmMat%DDMMat_DILU%nys(jyblock)
             nye = method%ddmMat%DDMMat_DILU%nye(jyblock)
             nxs = method%ddmMat%DDMMat_DILU%nxs(jxblock)
             nxe = method%ddmMat%DDMMat_DILU%nxe(jxblock)
             block_id = (jxblock-1) + (jyblock-1)*nxblock + (jzblock-1)*nxblock*nyblock 

             do jz = nzs-1,nze-1
                do jy = nys-1,nye-1
                   do jx = nxs-1,nxe-1
                      j  = (jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1

                      jjx = jx - (nxs-1)
                      jjy = jy - (nys-1)
                      jjz = jz - (nzs-1)
                      jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1                      
                      if(method%commVec%myrank.eq.0)then
!                      if(method%commVec%myrank.eq.1)then
                         write(*,*)" AS : ",jx,jy,jz,method%ddmMat%DDMMat_DILU%Dinv(block_id + jj*nblock)
                      endif

                   enddo
                enddo
             enddo

          enddo
       enddo
    end do

#endif

  end subroutine get_DILU_AS_DDM

  ! オーバーラップ領域の最大値は1ブロックサイズ
  subroutine solve_DILU_AS_DDM(method,r,s)
    use krylov
    use AS_util_ddm,only : AS_util_copy_to_AS,AS_util_copy_from_AS
    implicit none
    type(KSP) :: method

    real*8 r(method%ddmmat%n)
    real*8 s(method%ddmmat%n)

    real*8 zz

    integer nx,ny,nz
    integer stmx, stmy, stmz
    integer m
    integer mx,my,mz
    integer mxy

    integer m_AS
    integer nx_AS,ny_AS,nz_AS
    integer mxy_AS
    integer mx_AS,my_AS,mz_AS

    integer nblock
    integer nxblock,nyblock,nzblock
    integer nx_local,ny_local,nz_local
    integer m_local
    integer mx_local,my_local,mz_local

    integer nnd,offsetj,offsetjj
    integer offsetx,offsety,offsetz

    integer j,jj,jjj

    integer jx,jy,jz
    integer jjx,jjy,jjz
    integer jjjx,jjjy,jjjz
    integer block_id
    integer jxblock,jyblock,jzblock
    integer jjxblock,jjyblock,jjzblock
    integer jxblockm,jxblockp
    integer jyblockm,jyblockp
    integer jzblockm,jzblockp

    integer jjz0,jjy0,jjx0
    integer jjblock_id
    integer ip,ipp,ip_diag
    integer j_AS

    integer nxs,nxe
    integer nys,nye
    integer nzs,nze

    integer nxs0,nxe0
    integer nys0,nye0
    integer nzs0,nze0
    integer addzs,addys,addxs
    integer addze,addye,addxe
    integer offset
    integer iflag_AS,AS_mean

    nx = method%ddmMat%prm%nx
    ny = method%ddmMat%prm%ny
    nz = method%ddmMat%prm%nz

    stmx = method%ddmMat%prm%stmx
    stmy = method%ddmMat%prm%stmy
    stmz = method%ddmMat%prm%stmz

    mx  = method%ddmMat%prm%mx
    my  = method%ddmMat%prm%my
    mz  = method%ddmMat%prm%mz
    mxy = method%ddmMat%prm%mxy
    m   =  mx * my * mz

    ! ----------
    ! m_AS : 拡張領域の長さが必要 --> AS行列データは袖なしなのでそれを考慮
    ! ASベクトルも袖なし?
    nx_AS = method%ddmMat%AS_prm%nx
    ny_AS = method%ddmMat%AS_prm%ny
    nz_AS = method%ddmMat%AS_prm%nz

    mx_AS = nx_AS
    my_AS = ny_AS
    mz_AS = nz_AS
    mxy_AS = mx_AS * my_AS
    m_AS   = mx_AS * my_AS * mz_AS
    ! write(*,*)"m,m_AS = ",m,m_AS
    ! ----------

    nxblock = method%ddmMat%DDMMat_DILU%nxblock
    nyblock = method%ddmMat%DDMMat_DILU%nyblock
    nzblock = method%ddmMat%DDMMat_DILU%nzblock
    nblock  = nxblock*nyblock*nzblock;

    nx_local = method%ddmMat%DDMMat_DILU%nx_local
    ny_local = method%ddmMat%DDMMat_DILU%ny_local
    nz_local = method%ddmMat%DDMMat_DILU%nz_local

    mx_local =    method%ddmMat%DDMMat_DILU%mx_local
    my_local =    method%ddmMat%DDMMat_DILU%my_local
    mz_local =    method%ddmMat%DDMMat_DILU%mz_local
    m_local  =    method%ddmMat%DDMMat_DILU%m_local

    nnd = method%ddmMat%nnd

    ! -----
    ! まずはここで拡張領域の結果を受信する.
    ! とりあえずnoneは動作したことにする.
    ! rの受信を行う
    iflag_AS = method%ddmMat%precon%AS_METHOD_TYPE
    if(&
         (iflag_AS.eq.KSP_ADDITIVESCHWARZ_TYPE_INTERPOLATE)&
         .or. &
         (iflag_AS.eq.KSP_ADDITIVESCHWARZ_TYPE_BASIC)&
         )then
       AS_mean=1
    else
       AS_mean=0
    endif

    call AS_util_copy_to_AS ( &
         method%ddmMat%icomm, &
         method%ddmMat%myrank, &
         method%ddmMat%n, &
         method%ddmMat%prm, &
         method%ddmMat%AS_n, &
         method%ddmMat%AS_prm, &
         1, &
         method%ddmMat%AS_g, &
         r, &
         method%ddmMat%AS_v1, &
         iflag_AS )

    !$OMP PARALLEL DO PRIVATE(jx)
    do jx = 1,m_local * nblock
       method%ddmMat%DDMMat_DILU%svec(jx) = 0.0d0
       method%ddmMat%DDMMat_DILU%rvec(jx) = 0.0d0
       method%ddmMat%DDMMat_DILU%svec_AS(jx) = 0.0d0
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(&
    !$OMP jzblock,jyblock,jxblock, &
    !$OMP nzs,nze, &
    !$OMP nys,nye, &
    !$OMP nxs,nxe, &
    !$OMP block_id, &
    !$OMP j,jx,jy,jz, &
    !$OMP jj,jjx,jjy,jjz, &
    !$OMP j_AS, &
    !$OMP zz, &
    !$OMP ip, &
    !$OMP offset,offsetjj, &
    !$OMP offsetx,offsety,offsetz &
    !$OMP ) &
    !$OMP COLLAPSE(3)
    do jzblock=1,nzblock
       do jyblock=1,nyblock
          do jxblock=1,nxblock
             nzs = method%ddmMat%DDMMat_DILU%nzs(jzblock)
             nze = method%ddmMat%DDMMat_DILU%nze(jzblock)
             nys = method%ddmMat%DDMMat_DILU%nys(jyblock)
             nye = method%ddmMat%DDMMat_DILU%nye(jyblock)
             nxs = method%ddmMat%DDMMat_DILU%nxs(jxblock)
             nxe = method%ddmMat%DDMMat_DILU%nxe(jxblock)
             block_id = &
                  (jxblock-1) + (jyblock-1)*nxblock + (jzblock-1)*nxblock*nyblock 
             do jz = nzs-1,nze-1
                do jy = nys-1,nye-1
                   do jx = nxs-1,nxe-1
                      j  = (jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1

                      j_AS = &
                           (jx + method%ddmMat%AS_g%sleeve_mx) + &
                           (jy + method%ddmMat%AS_g%sleeve_my) * mx_AS  + &
                           (jz + method%ddmMat%AS_g%sleeve_mz) * mxy_AS + &
                           1

                      jjx = jx - (nxs-1)
                      jjy = jy - (nys-1)
                      jjz = jz - (nzs-1)
                      jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1
#if 0
                      method%ddmMat%DDMMat_DILU%rvec(block_id + jj*nblock) = r(j)
#else
                      method%ddmMat%DDMMat_DILU%rvec(block_id + jj*nblock) = &
                           method%ddmMat%AS_v1(j_AS)
#endif
                   enddo
                enddo
             enddo

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO 

    ! -----
    ! (D+L)D^-1(D+U)s=r
    ! -----
    !  1. (D+L)z=r       
    !     z = Dinv * (r - Lz)
    !  2. (I+D^-1U)s=z
    !     s = z-D^-1 Us       

    ! write(*,*)" solve_AS "

    !$OMP PARALLEL DO PRIVATE(&
    !$OMP jzblock,jyblock,jxblock, &
    !$OMP nzs,nze, &
    !$OMP nys,nye, &
    !$OMP nxs,nxe, &
    !$OMP nzs0,nze0, &
    !$OMP nys0,nye0, &
    !$OMP nxs0,nxe0, &
    !$OMP addzs,addys,addxs, &
    !$OMP addze,addye,addxe, &
    !$OMP block_id, &
    !$OMP j,jx,jy,jz, &
    !$OMP jj,jjx,jjy,jjz, &
    !$OMP j_AS, &
    !$OMP zz, &
    !$OMP ip, &
    !$OMP offset,offsetjj, &
    !$OMP offsetx,offsety,offsetz &
    !$OMP ) &
    !$OMP COLLAPSE(3)
    do jzblock=1,nzblock
       do jyblock=1,nyblock
          do jxblock=1,nxblock
             nzs = method%ddmMat%DDMMat_DILU%nzs(jzblock)
             nze = method%ddmMat%DDMMat_DILU%nze(jzblock)
             nys = method%ddmMat%DDMMat_DILU%nys(jyblock)
             nye = method%ddmMat%DDMMat_DILU%nye(jyblock)
             nxs = method%ddmMat%DDMMat_DILU%nxs(jxblock)
             nxe = method%ddmMat%DDMMat_DILU%nxe(jxblock)
             block_id = &
                  (jxblock-1) + (jyblock-1)*nxblock + (jzblock-1)*nxblock*nyblock 

             ! 前進代入
             do jz = nzs-1,nze-1
                do jy = nys-1,nye-1
                   do jx = nxs-1,nxe-1

                      j  = (jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1
                      ! j_AS = jx + mx_AS*jy + mxy_AS*jz + 1
#if 0
                      j_AS =  jx + jy * mx_AS + jz * mxy_AS + 1
#else
                      j_AS = &
                           (jx + method%ddmMat%AS_g%sleeve_mx) + &
                           (jy + method%ddmMat%AS_g%sleeve_my) * mx_AS  + &
                           (jz + method%ddmMat%AS_g%sleeve_mz) * mxy_AS + &
                           1
#endif
                      jjx = jx - (nxs-1)
                      jjy = jy - (nys-1)
                      jjz = jz - (nzs-1)
                      jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1

                      zz = method%ddmMat%DDMMat_DILU%rvec(block_id + jj*nblock)
                      do ip = 1,nnd
                         offset   = method%ddmMat%mat%ioff_dia(ip)
                         offsetx  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*0)
                         offsety  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*1)
                         offsetz  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*2)
                         offsetjj = offsetx + offsety*mx_local + offsetz*mx_local*my_local
                         if(offset>=0)then
                            cycle
                         endif
                         ! zz = zz  & 
                         !      - method%ddmMat%DDMMat_DILU%svec(block_id + (jj+offsetjj)*nblock) * &
                         !        method%ddmMat%mat%val_dia(m*(ip-1)+j)
                         zz = zz  & 
                              - method%ddmMat%DDMMat_DILU%svec(block_id + (jj+offsetjj)*nblock) * &
                                method%ddmMat%ASval%AS_val_dia_init(m_AS*(ip-1)+j_AS)

                      enddo
                      method%ddmMat%DDMMat_DILU%svec(block_id + jj*nblock) = &
                           method%ddmMat%DDMMat_DILU%Dinv(block_id + jj*nblock)*zz
                   enddo
                enddo
             enddo

             ! 後退代入
             do jz = nze-1,nzs-1,-1
                do jy = nye-1,nys-1,-1
                   do jx = nxe-1,nxs-1,-1

                      j  = (jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1
                      ! j_AS = jx + mx_AS*jy + mxy_AS*jz + 1
#if 0
                      j_AS =  jx + jy * mx_AS + jz * mxy_AS + 1
#else
                      j_AS = &
                           (jx + method%ddmMat%AS_g%sleeve_mx) + &
                           (jy + method%ddmMat%AS_g%sleeve_my) * mx_AS  + &
                           (jz + method%ddmMat%AS_g%sleeve_mz) * mxy_AS + &
                           1
#endif
                      jjx = jx - (nxs-1)
                      jjy = jy - (nys-1)
                      jjz = jz - (nzs-1)
                      jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1

                      zz = 0.0d0
                      do ip = 1,nnd
                         offset = method%ddmMat%mat%ioff_dia(ip)
                         offsetx  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*0)
                         offsety  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*1)
                         offsetz  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*2)
                         offsetjj = offsetx + offsety*mx_local + offsetz*mx_local*my_local
                         if(offset<=0)then
                            cycle
                         endif
                         ! zz = zz  & 
                         !      + method%ddmMat%DDMMat_DILU%svec(block_id + (jj+offsetjj)*nblock) * &
                         !        method%ddmMat%mat%val_dia(m*(ip-1)+j)
                         zz = zz  & 
                              + method%ddmMat%DDMMat_DILU%svec(block_id + (jj+offsetjj)*nblock) * &
                                method%ddmMat%ASval%AS_val_dia_init(m_AS*(ip-1)+j_AS)

                      enddo
                      method%ddmMat%DDMMat_DILU%svec(block_id + jj*nblock) =  &
                           method%ddmMat%DDMMat_DILU%svec(block_id + jj*nblock) - &
                           method%ddmMat%DDMMat_DILU%Dinv(block_id + jj*nblock)*zz
                      ! sの更新は拡張領域を修正した後にする?
#if AS_Debug
                      s(j) = method%ddmMat%DDMMat_DILU%svec(block_id + jj*nblock)
#endif
                      ! 前処理無し
!                      s(j) = method%ddmMat%DDMMat_DILU%rvec(block_id + jj*nblock)
                   enddo
                enddo
             enddo
#if AS_Debug

#else

#if 0
             ! ------------------------
             ! 結果の書き込み
             nzs0 = method%ddmMat%DDMMat_DILU%nzs0(jzblock)
             nys0 = method%ddmMat%DDMMat_DILU%nys0(jyblock)
             nxs0 = method%ddmMat%DDMMat_DILU%nxs0(jxblock)

             nze0 = method%ddmMat%DDMMat_DILU%nze0(jzblock)
             nye0 = method%ddmMat%DDMMat_DILU%nye0(jyblock)
             nxe0 = method%ddmMat%DDMMat_DILU%nxe0(jxblock)

             addzs = nzs0 - nzs
             addys = nys0 - nys
             addxs = nxs0 - nxs

             do jz = nzs0-1,nze0-1
                do jy = nys0-1,nye0-1
                   do jx = nxs0-1,nxe0-1
                      j  = (jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1
                      jjx = jx - (nxs0-1) + addxs
                      jjy = jy - (nys0-1) + addys
                      jjz = jz - (nzs0-1) + addzs
                      jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1
                      s(j) = method%ddmMat%DDMMat_DILU%svec(block_id + jj*nblock)
!                      s(j) = method%ddmMat%DDMMat_DILU%rvec(block_id + jj*nblock)

                   enddo
                enddo
             enddo
             ! ------------------------
#endif

#endif

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO 

    do jzblock=1,nzblock
       do jyblock=1,nyblock
          do jxblock=1,nxblock
             nzs = method%ddmMat%DDMMat_DILU%nzs(jzblock)
             nze = method%ddmMat%DDMMat_DILU%nze(jzblock)
             nys = method%ddmMat%DDMMat_DILU%nys(jyblock)
             nye = method%ddmMat%DDMMat_DILU%nye(jyblock)
             nxs = method%ddmMat%DDMMat_DILU%nxs(jxblock)
             nxe = method%ddmMat%DDMMat_DILU%nxe(jxblock)
             block_id = &
                  (jxblock-1) + (jyblock-1)*nxblock + (jzblock-1)*nxblock*nyblock 
             jxblockm = max(jxblock-1,1)
             jxblockp = min(jxblock+1,nxblock)
             jyblockm = max(jyblock-1,1)
             jyblockp = min(jyblock+1,nyblock)
             jzblockm = max(jzblock-1,1)
             jzblockp = min(jzblock+1,nzblock)

             do jjzblock = jzblockm, jzblockp
                select case (jjzblock - jzblock)
                case(-1)
                   nzs = method%ddmMat%DDMMat_DILU%nzs(jzblock)
                   nze = method%ddmMat%DDMMat_DILU%nzs(jzblock) + method%ddmMat%AS_g%addL - 1
                   jjz0 = (method%ddmMat%DDMMat_DILU%nze(jjzblock) - method%ddmMat%DDMMat_DILU%nzs(jjzblock) + 1)  - method%ddmMat%AS_g%addL 
                case(0)
                   nzs = method%ddmMat%DDMMat_DILU%nzs(jzblock)
                   nze = method%ddmMat%DDMMat_DILU%nze(jzblock)
                   jjz0 = 0
                case(1)
                   nzs = method%ddmMat%DDMMat_DILU%nze(jzblock) - method%ddmMat%AS_g%addL + 1
                   nze = method%ddmMat%DDMMat_DILU%nze(jzblock) 
                   jjz0 = 0
                end select

                do jjyblock = jyblockm, jyblockp
                   select case (jjyblock - jyblock)
                   case(-1)
                      nys = method%ddmMat%DDMMat_DILU%nys(jyblock)
                      nye = method%ddmMat%DDMMat_DILU%nys(jyblock) + method%ddmMat%AS_g%addL - 1
                      jjy0 = (method%ddmMat%DDMMat_DILU%nye(jjyblock) - method%ddmMat%DDMMat_DILU%nys(jjyblock) + 1)  - method%ddmMat%AS_g%addL 
                   case(0)
                      nys = method%ddmMat%DDMMat_DILU%nys(jyblock)
                      nye = method%ddmMat%DDMMat_DILU%nye(jyblock)
                      jjy0 = 0
                   case(1)
                      nys = method%ddmMat%DDMMat_DILU%nye(jyblock) - method%ddmMat%AS_g%addL + 1
                      nye = method%ddmMat%DDMMat_DILU%nye(jyblock) 
                      jjy0 = 0
                   end select

                   do jjxblock = jxblockm, jxblockp
                      select case (jjxblock - jxblock)
                      case(-1)
                         nxs = method%ddmMat%DDMMat_DILU%nxs(jxblock)
                         nxe = method%ddmMat%DDMMat_DILU%nxs(jxblock) + method%ddmMat%AS_g%addL - 1
                         jjx0 = (method%ddmMat%DDMMat_DILU%nxe(jjxblock) - method%ddmMat%DDMMat_DILU%nxs(jjxblock) + 1)  - method%ddmMat%AS_g%addL 
                      case(0)
                         nxs = method%ddmMat%DDMMat_DILU%nxs(jxblock)
                         nxe = method%ddmMat%DDMMat_DILU%nxe(jxblock)
                         jjx0 = 0
                      case(1)
                         nxs = method%ddmMat%DDMMat_DILU%nxe(jxblock) - method%ddmMat%AS_g%addL + 1
                         nxe = method%ddmMat%DDMMat_DILU%nxe(jxblock) 
                         jjx0 = 0
                      end select
                      jjblock_id = (jjxblock-1) + (jjyblock-1)*nxblock + (jjzblock-1)*nxblock*nyblock 

                      ! write(*,*)jxblock,jyblock,jzblock,jjxblock,jjyblock,jjzblock
                      ! write(*,*)"jz = ", nzs-1,nze-1, " : jjz0 = ",jjz0
                      ! write(*,*)"jy = ", nys-1,nye-1, " : jjy0 = ",jjy0
                      ! write(*,*)"jx = ", nxs-1,nxe-1, " : jjx0 = ",jjx0

!                      if(jjblock_id.ne.block_id)then
                         do jz = nzs-1,nze-1
                            do jy = nys-1,nye-1
                               do jx = nxs-1,nxe-1
                                  jjx = jx - (method%ddmMat%DDMMat_DILU%nxs(jxblock)-1)
                                  jjy = jy - (method%ddmMat%DDMMat_DILU%nys(jyblock)-1)
                                  jjz = jz - (method%ddmMat%DDMMat_DILU%nzs(jzblock)-1)
                                  jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1

                                  jjjx = jx - (nxs-1) + jjx0
                                  jjjy = jy - (nys-1) + jjy0
                                  jjjz = jz - (nzs-1) + jjz0
                                  jjj = (jjjx+stmx) + mx_local*(jjjy+stmy) + mx_local*my_local*(jjjz+stmz) + 1
                                  ! write(*,*)"jx,jjjx = ",jx,jjjx
                                  ! write(*,*)"jy,jjjy = ",jy,jjjy
                                  ! write(*,*)"jz,jjjz = ",jz,jjjz
                                  ! write(*,*)"j + jj = ",block_id + jj*nblock,jjblock_id + jjj*nblock
                                  method%ddmMat%DDMMat_DILU%svec_AS(block_id + jj*nblock) = &
                                       method%ddmMat%DDMMat_DILU%svec_AS(block_id   + jj*nblock) + &
                                       method%ddmMat%DDMMat_DILU%svec(jjblock_id + jjj*nblock) 
                               enddo
                            enddo
                         enddo
!                      endif

                   enddo
                enddo
             enddo

             nzs = method%ddmMat%DDMMat_DILU%nzs(jzblock)
             nze = method%ddmMat%DDMMat_DILU%nze(jzblock)
             nys = method%ddmMat%DDMMat_DILU%nys(jyblock)
             nye = method%ddmMat%DDMMat_DILU%nye(jyblock)
             nxs = method%ddmMat%DDMMat_DILU%nxs(jxblock)
             nxe = method%ddmMat%DDMMat_DILU%nxe(jxblock)

             nzs0 = method%ddmMat%DDMMat_DILU%nzs0(jzblock)
             nys0 = method%ddmMat%DDMMat_DILU%nys0(jyblock)
             nxs0 = method%ddmMat%DDMMat_DILU%nxs0(jxblock)

             nze0 = method%ddmMat%DDMMat_DILU%nze0(jzblock)
             nye0 = method%ddmMat%DDMMat_DILU%nye0(jyblock)
             nxe0 = method%ddmMat%DDMMat_DILU%nxe0(jxblock)

             addzs = nzs0 - nzs
             addys = nys0 - nys
             addxs = nxs0 - nxs
             do jz = nzs0-1,nze0-1
                do jy = nys0-1,nye0-1
                   do jx = nxs0-1,nxe0-1
                      j  = (jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1
                      jjx = jx - (nxs0-1) + addxs
                      jjy = jy - (nys0-1) + addys
                      jjz = jz - (nzs0-1) + addzs
                      jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1
                      ! s(j) = method%ddmMat%DDMMat_DILU%svec(block_id + jj*nblock)
                      s(j) = method%ddmMat%DDMMat_DILU%svec_AS(block_id + jj*nblock)
                      ! 前処理無し
                      ! s(j) = method%ddmMat%DDMMat_DILU%rvec(block_id + jj*nblock)
                      ! method%ddmMat%AS_v2(j_AS)

                   enddo
                enddo
             enddo

             if(jzblock.eq.1)then
                nzs0  = nzs
                addzs = 0
             endif
             if(jzblock.eq.nzblock)then
                nze0  = nze
             endif

             if(jyblock.eq.1)then
                nys0  = nys
                addys = 0
             endif
             if(jyblock.eq.nyblock)then
                nye0  = nye
             endif

             if(jxblock.eq.1)then
                nxs0  = nxs
                addxs = 0
             endif
             if(jxblock.eq.nxblock)then
                nxe0  = nxe
             endif
             ! -------
             do jz = nzs0-1,nze0-1
                do jy = nys0-1,nye0-1
                   do jx = nxs0-1,nxe0-1
                      j_AS = &
                           (jx + method%ddmMat%AS_g%sleeve_mx) + &
                           (jy + method%ddmMat%AS_g%sleeve_my) * mx_AS  + &
                           (jz + method%ddmMat%AS_g%sleeve_mz) * mxy_AS + &
                           1                     
                      jjx = jx - (nxs0-1) + addxs
                      jjy = jy - (nys0-1) + addys
                      jjz = jz - (nzs0-1) + addzs
                      jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1
                      ! method%ddmMat%AS_v2(j_AS) = method%ddmMat%DDMMat_DILU%svec(block_id + jj*nblock) 
                      method%ddmMat%AS_v2(j_AS) = method%ddmMat%DDMMat_DILU%svec_AS(block_id + jj*nblock) 
                   enddo
                enddo
             enddo

          enddo
       enddo
    enddo

    call AS_util_copy_from_AS ( &
         method%ddmMat%icomm, &
         method%ddmMat%myrank, &
         method%ddmMat%n, &
         method%ddmMat%prm, &
         method%ddmMat%AS_n, &
         method%ddmMat%AS_prm, &
         1, &
         AS_mean, &
         method%ddmMat%AS_g, &
         s, &
         method%ddmMat%AS_v2 )

    
#if 0
    do jzblock=1,nzblock
       do jyblock=1,nyblock
          do jxblock=1,nxblock
             nzs = method%ddmMat%DDMMat_DILU%nzs(jzblock)
             nze = method%ddmMat%DDMMat_DILU%nze(jzblock)
             nys = method%ddmMat%DDMMat_DILU%nys(jyblock)
             nye = method%ddmMat%DDMMat_DILU%nye(jyblock)
             nxs = method%ddmMat%DDMMat_DILU%nxs(jxblock)
             nxe = method%ddmMat%DDMMat_DILU%nxe(jxblock)
             block_id = &
                  (jxblock-1) + (jyblock-1)*nxblock + (jzblock-1)*nxblock*nyblock 

             nzs0 = method%ddmMat%DDMMat_DILU%nzs0(jzblock)
             nys0 = method%ddmMat%DDMMat_DILU%nys0(jyblock)
             nxs0 = method%ddmMat%DDMMat_DILU%nxs0(jxblock)

             nze0 = method%ddmMat%DDMMat_DILU%nze0(jzblock)
             nye0 = method%ddmMat%DDMMat_DILU%nye0(jyblock)
             nxe0 = method%ddmMat%DDMMat_DILU%nxe0(jxblock)

             addzs = nzs0 - nzs
             addys = nys0 - nys
             addxs = nxs0 - nxs
             do jz = nzs0-1,nze0-1
                do jy = nys0-1,nye0-1
                   do jx = nxs0-1,nxe0-1
                      j  = (jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1
                      jjx = jx - (nxs0-1) + addxs
                      jjy = jy - (nys0-1) + addys
                      jjz = jz - (nzs0-1) + addzs
                      jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1
                      if(method%commVec%myrank.eq.1)then
                         write(*,*)j,"AS : jx,jy,jz = ",jx,jy,jz
                      endif
                   enddo
                enddo
             enddo

          enddo
       enddo
    enddo

    call system("sleep 5")
    stop

#endif

  end subroutine solve_DILU_AS_DDM

  ! -------------------------------------
  subroutine solve_DILU_AS_DDM_GPU(method,r,s)
    use krylov,only : KSP
    use vec_type,only : PARCEL_VEC
    use PARCEL_VEC_kernel,only : parcel_sync
    use PARCEL_define,only : &
         KSP_ADDITIVESCHWARZ_TYPE_INTERPOLATE, &
         KSP_ADDITIVESCHWARZ_TYPE_BASIC
!    use c_wrapper,only : devicemem
    use AS_util_ddm,only : AS_util_copy_to_AS,AS_util_copy_from_AS
    use cu_common_wrapper,only : cudamemcpy_d2h_fp64,cudamemcpy_h2d_fp64
    implicit none
    type(KSP) :: method
    integer(kind=intType),parameter :: offset_zero=0

!    real*8 r(method%ddmmat%n)
!    real*8 s(method%ddmmat%n)
    type(PARCEL_VEC) :: r
    type(PARCEL_VEC) :: s
    real*8 zz

    integer nx,ny,nz
    integer stmx, stmy, stmz
    integer m
    integer mx,my,mz
    integer mxy

    integer m_AS
    integer nx_AS,ny_AS,nz_AS
    integer mxy_AS
    integer mx_AS,my_AS,mz_AS

    integer nblock
    integer nxblock,nyblock,nzblock
    integer nx_local,ny_local,nz_local
    integer m_local
    integer mx_local,my_local,mz_local

    integer nnd,offsetj,offsetjj
    integer offsetx,offsety,offsetz

    integer j,jj,jjj

    integer jx,jy,jz
    integer jjx,jjy,jjz
    integer jjjx,jjjy,jjjz
    integer block_id
    integer jxblock,jyblock,jzblock
    integer jjxblock,jjyblock,jjzblock
    integer jxblockm,jxblockp
    integer jyblockm,jyblockp
    integer jzblockm,jzblockp

    integer jjz0,jjy0,jjx0
    integer jjblock_id
    integer ip,ipp,ip_diag
    integer j_AS

    integer nxs,nxe
    integer nys,nye
    integer nzs,nze

    integer nxs0,nxe0
    integer nys0,nye0
    integer nzs0,nze0
    integer addzs,addys,addxs
    integer addze,addye,addxe
    integer offset
    integer iflag_AS,AS_mean
    integer(kind=intType) parcel_int

    nx = method%ddmMat%prm%nx
    ny = method%ddmMat%prm%ny
    nz = method%ddmMat%prm%nz

    stmx = method%ddmMat%prm%stmx
    stmy = method%ddmMat%prm%stmy
    stmz = method%ddmMat%prm%stmz

    mx  = method%ddmMat%prm%mx
    my  = method%ddmMat%prm%my
    mz  = method%ddmMat%prm%mz
    mxy = method%ddmMat%prm%mxy
    m   =  mx * my * mz

    ! ----------
    ! m_AS : 拡張領域の長さが必要 --> AS行列データは袖なしなのでそれを考慮
    ! ASベクトルも袖なし?
    nx_AS = method%ddmMat%AS_prm%nx
    ny_AS = method%ddmMat%AS_prm%ny
    nz_AS = method%ddmMat%AS_prm%nz

    mx_AS = nx_AS
    my_AS = ny_AS
    mz_AS = nz_AS
    mxy_AS = mx_AS * my_AS
    m_AS   = mx_AS * my_AS * mz_AS
    ! ----------

    nxblock = method%ddmMat%DDMMat_DILU%nxblock
    nyblock = method%ddmMat%DDMMat_DILU%nyblock
    nzblock = method%ddmMat%DDMMat_DILU%nzblock
    nblock  = nxblock*nyblock*nzblock;

    nx_local = method%ddmMat%DDMMat_DILU%nx_local
    ny_local = method%ddmMat%DDMMat_DILU%ny_local
    nz_local = method%ddmMat%DDMMat_DILU%nz_local

    mx_local = method%ddmMat%DDMMat_DILU%mx_local
    my_local = method%ddmMat%DDMMat_DILU%my_local
    mz_local = method%ddmMat%DDMMat_DILU%mz_local
    m_local  = method%ddmMat%DDMMat_DILU%m_local

    nnd = method%ddmMat%nnd

    ! -----
    iflag_AS = method%ddmMat%precon%AS_METHOD_TYPE
    if(&
         (iflag_AS.eq.KSP_ADDITIVESCHWARZ_TYPE_INTERPOLATE)&
         .or. &
         (iflag_AS.eq.KSP_ADDITIVESCHWARZ_TYPE_BASIC)&
         )then
       AS_mean=1
    else
       AS_mean=0
    endif

    ! ! ここもGPU側のメモリを初期化
    ! !$OMP PARALLEL DO PRIVATE(jx)
    ! do jx = 1,m_local * nblock
    !    method%ddmMat%DDMMat_DILU%svec(jx) = 0.0d0
    !    method%ddmMat%DDMMat_DILU%rvec(jx) = 0.0d0
    !    method%ddmMat%DDMMat_DILU%svec_AS(jx) = 0.0d0
    ! enddo
    ! !$OMP END PARALLEL DO
 
    if(allocated(r%vec_FP64))then
       deallocate(r%vec_FP64)
    endif
    if(allocated(s%vec_FP64))then
       deallocate(s%vec_FP64)
    endif
    allocate( r%vec_FP64(m) )
    allocate( s%vec_FP64(m) )
    !$OMP PARALLEL DO PRIVATE(jx)
    do jx = 1,m
       r%vec_FP64(jx) = 0.0d0
       s%vec_FP64(jx) = 0.0d0
    enddo
    !$OMP END PARALLEL DO

#if 1
    ! GPU側のrをCPU側に送る
    parcel_int = m
    call cudaMemcpy_d2h_FP64(&
         r%C_vec_FP64, &
         r%vec_FP64, &
         parcel_int,&
         r%offset, &
         offset_zero  &
         ) 
    call parcel_sync(method)
    call AS_util_copy_to_AS ( &
         method%ddmMat%icomm, &
         method%ddmMat%myrank, &
         method%ddmMat%n, &
         method%ddmMat%prm, &
         method%ddmMat%AS_n, &
         method%ddmMat%AS_prm, &
         1, &
         method%ddmMat%AS_g, &
         r%vec_FP64, &
         method%ddmMat%AS_v1, &
         iflag_AS )
    call parcel_sync(method)

    ! CPU側のAS_v1をCPU側に送る
    parcel_int = size( method%ddmMat%AS_v1 )
    call cudaMemcpy_h2d_FP64(&
         method%ddmMat%AS_v1,method%ddmMat%C_AS_v1, parcel_int,&
         offset_zero, &
         offset_zero) 
    call parcel_sync(method)
#else
    call parcel_sync(method)
    call AS_util_copy_to_AS_GPU( &
         method, &  
         method%ddmMat%icomm,  &
         method%ddmMat%myrank,  &
         method%ddmMat%n,  &
         method%ddmMat%prm,  &
         method%ddmMat%AS_n,  &
         method%ddmMat%AS_prm, &
         1, &
         method%ddmMat%AS_g, &
         r%C_vec_FP64, r%offset, &
         method%ddmMat%C_AS_v1, &
         offset_zero, &
         iflag_AS)
#endif
          call parcel_sync(method)
          call solve_DILU_AS_DDM_FP64_GPU(&
               method%ddmMat%C_AS_v1 , offset_zero, &
               s%C_vec_FP64, s%offset, &
               method%ddmMat%C_DDMMat_DILU%rvec      , offset_zero, &
               method%ddmMat%C_DDMMat_DILU%svec      , offset_zero, &
               method%ddmMat%C_DDMMat_DILU%svec_AS   , offset_zero, &
               method%ddmMat%C_AS_v2                 , offset_zero, &
               method%ddmMat%ASval%C_AS_val_dia_init ,  &
               method%ddmMat%C_DDMMat_DILU%Dinv, offset_zero, &
               method%ddmMat%mat%C_ioff_dia, &
               method%ddmMat%C_DDMMat_DILU%offset_xyz, &
               nnd,  &
               method%ddmMat%C_DDMMat_DILU%nxs, method%ddmMat%C_DDMMat_DILU%nxe, &
               method%ddmMat%C_DDMMat_DILU%nys, method%ddmMat%C_DDMMat_DILU%nye, &
               method%ddmMat%C_DDMMat_DILU%nzs, method%ddmMat%C_DDMMat_DILU%nze, &
               method%ddmMat%C_DDMMat_DILU%nxs0, method%ddmMat%C_DDMMat_DILU%nxe0, &
               method%ddmMat%C_DDMMat_DILU%nys0, method%ddmMat%C_DDMMat_DILU%nye0, &
               method%ddmMat%C_DDMMat_DILU%nzs0, method%ddmMat%C_DDMMat_DILU%nze0, &
               method%ddmMat%AS_g%addL, &
               nblock, &
               nxblock,nyblock,nzblock, &
               stmx,stmy,stmz, &
               mx, mxy, &
               m_local, &
               mx_local,my_local, &
               m_AS, &
               mx_AS ,mxy_AS, &
               method%ddmMat%AS_g%sleeve_mx, &
               method%ddmMat%AS_g%sleeve_my, &
               method%ddmMat%AS_g%sleeve_mz &
               )

#if 1
          call parcel_sync(method)
          ! GPU側のAS_v2をCPU側に送信
          parcel_int = size(method%ddmMat%AS_v2)
           call cudaMemcpy_d2h_FP64(&
               method%ddmMat%C_AS_v2, &
               method%ddmMat%AS_v2, &
               parcel_int, &
               offset_zero, &
               offset_zero  &
               ) 
          call parcel_sync(method)
          call AS_util_copy_from_AS ( &
               method%ddmMat%icomm, &
               method%ddmMat%myrank, &
               method%ddmMat%n, &
               method%ddmMat%prm, &
               method%ddmMat%AS_n, &
               method%ddmMat%AS_prm, &
               1, &
               AS_mean, &
               method%ddmMat%AS_g, &
               s%vec_FP64, &
               method%ddmMat%AS_v2 )
          call parcel_sync(method)

          parcel_int = m
          call cudaMemcpy_h2d_FP64(&
               s%vec_FP64, &
               s%C_vec_FP64, &
               parcel_int,&
               offset_zero, &
               s%offset) 
          call parcel_sync(method)

#else
    call parcel_sync(method)
    call AS_util_copy_from_AS_GPU( &
         method, &  
         method%ddmMat%icomm, &
         method%ddmMat%myrank, &
         method%ddmMat%n,  &
         method%ddmMat%prm,  &
         method%ddmMat%AS_n,  &
         method%ddmMat%AS_prm, &
         1, &
         AS_mean, &
         method%ddmMat%AS_g, &
         s%C_vec_FP64, s%offset, &
         method%ddmMat%C_AS_v2, &
         offset_zero &
         )
#endif

    ! CPU側のsをGPU側に送信
    if(allocated(r%vec_FP64))then
       deallocate(r%vec_FP64)
    endif
    if(allocated(s%vec_FP64))then
       deallocate(s%vec_FP64)
    endif

  end subroutine solve_DILU_AS_DDM_GPU

  subroutine AS_util_copy_to_AS_GPU( &
       method, &  
       comm, myrank, ndist, prm, AS_ndist, AS_prm, &
       AS_diag_neighbor, AS_g, &
       wv1,    offset_wv1,&
       AS_wv1, offset_AS_wv1,&
       iflag_AS)
!    use param_AS_type_ddm,only : AS_geom
    use krylov,only : KSP
    use PARCEL_VEC_kernel,only : parcel_sync
    use ddmMat_wrap,only : ddm_param,AS_geom
    use c_wrapper,only : devicemem
    use cu_common_wrapper,only : &
         cudaMalloc_FP64,cudaFree, &
         C_MPI_Isend_FP64,C_MPI_Irecv_FP64,C_MPI_Wait, &
         cudamemcpy_d2h_fp64,cudamemcpy_h2d_fp64

    use PARCEL_define,only : &
         KSP_ADDITIVESCHWARZ_TYPE_INTERPOLATE, &
         KSP_ADDITIVESCHWARZ_TYPE_NONE, &
         KSP_ADDITIVESCHWARZ_TYPE_BASIC

    use mpi
    implicit none
    type(KSP) :: method
    !   include "ddm.h"
    !   include "mpif.h"
    integer        , intent(in ) :: comm   !
    integer        , intent(in ) :: myrank !
    integer(kind=intType), intent(in ) :: ndist  ! 
    type(ddm_param), intent(in ) :: prm    ! 
    integer(kind=intType), intent(in ) :: AS_ndist         ! 領域のグリッド数（AS = 加法シュワルツ）
    type(ddm_param), intent(in ) :: AS_prm           ! 領域形状
    integer        , intent(in ) :: AS_diag_neighbor ! 斜め隣の領域の値を使うなら1、使わないなら0
    type(AS_geom)  , intent(in ) :: AS_g
    ! real*8         , intent(in ) :: wv1(ndist)       ! 
    ! real*8         , intent(inout) :: AS_wv1(AS_ndist) ! 前処理行列の作用対象となるベクトル
    type(devicemem) wv1
    integer(kind=intType) offset_wv1
    type(devicemem) AS_wv1
    integer(kind=intType) offset_AS_wv1

    integer iflag_as ! 加法シュワルツの補正方法  !
    integer(kind=inttype) k, as_i, i, ix, iy, iz, j
    integer lenx0, leny0, lenz0, lenxy0, lenxyz0
    integer xyz, pm, info
    type gpu_buff
       type(devicemem) v
    end type gpu_buff
    type(gpu_buff) sbuf(6), rbuf(6) ! 6 = ±x, ±y, ±z   

    type cpu_buff
       real*8, allocatable :: v(:)
    end type cpu_buff
    type(cpu_buff) sbuf_host(6), rbuf_host(6) ! 6 = ±x, ±y, ±z   

    integer     ireqs(6), ireqr(6)
    integer     istatus (mpi_status_size)
    integer(kind=inttype),parameter :: offset_zero=0
    integer(kind=inttype) parcel_int

    call as_util_copy_only_fp64_gpu(&
         wv1   , offset_wv1, &
         as_wv1, offset_as_wv1, &
         prm%nx,prm%ny,prm%nz, &
         prm%mx,prm%mxy, &
         prm%stmx,prm%stmy,prm%stmz, &
         as_ndist, &
         as_prm%mx,as_prm%mxy, &
         as_g%sleeve_mx,as_g%sleeve_my,as_g%sleeve_mz &
         )
    if(&
         (iflag_as.eq.ksp_additiveschwarz_type_interpolate) &
         .or. &
         (iflag_as.eq.ksp_additiveschwarz_type_none) &
         )then
       return
    endif

    !!| 隣接プロセスからデータを通信により取得する
    !!| wv1のデータの代わりに、上でコピーしたas_wv1のデータを通信する
    do i=1, 6
       parcel_int = as_g%grids_s(i)
       call cudamalloc_fp64(sbuf(i)%v, parcel_int ) 
       parcel_int = as_g%grids_r(i)
       call cudamalloc_fp64(rbuf(i)%v, parcel_int ) 
    end do

    do xyz=0, 2  ! x, y, z
       do pm=0, 1  ! +xyz, -xyz
          i = 2 * xyz + pm + 1
          ! 通信データの設定
          if ( as_g%rank_s(i) /= -1 ) then
             call as_util_set_sbuf_fp64_gpu(&
			sbuf(i)%v , offset_zero, &
			as_wv1, offset_as_wv1, &
			as_g%grids_s(i), &
			as_prm%mx, as_prm%mxy, &
			as_g%xlen_s(i),as_g%ylen_s(i), &
			as_g%x0s(i),as_g%y0s(i),as_g%z0s(i)  &
                        ) 
          end if
       enddo
    enddo

    call parcel_sync(method)

#if CUDA_AWARE_

    do xyz=0, 2  ! x, y, z
       do pm=0, 1  ! +xyz, -xyz
          i = 2 * xyz + pm + 1
          ! 通信開始
          if ( AS_g%rank_s(i) /= -1 ) then
             call c_mpi_isend_fp64(&
                  sbuf(i)%v, &
                  AS_g%grids_s(i), &
                  AS_g%rank_s(i), &
                  AS_g%rank_s(i), &
                  method%commvec%c_comm, &
                  method%commvec%c_s_req, &
                  i,1+offset_zero &
                  ) 
          end if
          if ( AS_g%rank_r(i) /= -1 ) then
             call c_mpi_irecv_fp64(&
                  rbuf(i)%v, &
                  AS_g%grids_r(i), &
                  AS_g%rank_r(i), &
                  myrank, &
                  method%commvec%c_comm, &
                  method%commvec%c_r_req, &
                  i,1+offset_zero &
                  ) 
          end if
       end do
    enddo

    call parcel_sync(method)
    ! 通信を完了させて、受信したデータをコピー
    do xyz=0, 2  ! x, y, z
       do pm=0, 1
          i = 2 * xyz + pm + 1
          if ( as_g%rank_s(i) /= -1 ) &
               call c_mpi_wait(method%commvec%c_s_req,method%commvec%c_s_sta,i) 
          if ( as_g%rank_r(i) /= -1 ) &
               call c_mpi_wait(method%commvec%c_r_req,method%commvec%c_r_sta,i) 
       enddo
    enddo

#else

    do i=1, 6
       allocate( sbuf_host(i)%v(AS_g%grids_s(i)) )
       allocate( rbuf_host(i)%v(AS_g%grids_r(i)) )
    end do

    do xyz=0, 2  ! x, y, z
       do pm=0, 1  ! +xyz, -xyz
          i = 2 * xyz + pm + 1

          parcel_int = AS_g%grids_s(i)
          call cudaMemcpy_d2h_FP64(&
               sbuf(i)%v ,sbuf_host(i)%v, parcel_int,&
               offset_zero, &
               offset_zero  &
               ) 
       enddo
    enddo

    call parcel_sync(method)
    do xyz=0, 2  ! x, y, z
       do pm=0, 1  ! +xyz, -xyz
          i = 2 * xyz + pm + 1
          ! 通信開始
          if ( as_g%rank_s(i) /= -1 ) then
             call MPI_Isend( sbuf_host(i)%v, &
                  size(sbuf_host(i)%v), &
                  mpi_double_precision, &
                  AS_g%rank_s(i), &
                  AS_g%rank_s(i), &
                  comm, &
                  ireqs(i), &
                  info )
          end if
          if ( as_g%rank_r(i) /= -1 ) then
             call MPI_Irecv( rbuf_host(i)%v, &
                  size(rbuf_host(i)%v), &
                  mpi_double_precision, &
                  AS_g%rank_r(i), &
                  myrank, &
                  comm, &
                  ireqr(i), &
                  info )
          end if
       end do
    enddo
    call parcel_sync(method)
    do xyz=0, 2  ! x, y, z
       do pm=0, 1  ! +xyz, -xyz
          i = 2 * xyz + pm + 1
          if ( AS_g%rank_s(i) /= -1 ) call MPI_Wait( ireqs(i), istatus, info )
          if ( AS_g%rank_r(i) /= -1 ) call MPI_Wait( ireqr(i), istatus, info )
       enddo
    enddo
    do i=1, 6
       parcel_int = AS_g%grids_r(i)
       call cudaMemcpy_h2d_FP64(&
         rbuf_host(i)%v, &
         rbuf(i)%v, parcel_int,&
         offset_zero, &
         offset_zero) 
    enddo

    call parcel_sync(method)
    do i=1, 6
       deallocate( sbuf_host(i)%v )
       deallocate( rbuf_host(i)%v )
    end do

#endif

    call parcel_sync(method)
    do xyz=0, 2  ! x, y, z
       ! 通信を完了させて、受信したデータをコピー
       do pm=0, 1
          i = 2 * xyz + pm + 1
          ! 受信したデータをコピー
          if ( as_g%rank_r(i) /= -1 ) then
             ! !$omp parallel do private(ix, iy, iz, k)
             ! do j=0, size(rbuf(i)%v)-1
             !    ix = as_g%x0r(i) + mod(j, as_g%xlen_r(i))
             !    iy = as_g%y0r(i) + mod(j/as_g%xlen_r(i),  as_g%ylen_r(i))
             !    iz = as_g%z0r(i) + j / (as_g%xlen_r(i) * as_g%ylen_r(i))
             !    k = ix + as_prm%mx * iy + as_prm%mxy * iz + 1
             !    as_wv1(k) = rbuf(i)%v(j+1)
             ! end do
             ! !$omp end parallel do
             call AS_util_set_rbuf_FP64_GPU(&
			rbuf(i)%v , offset_zero, &
			as_wv1, offset_as_wv1, &
			as_g%grids_r(i), &
			as_prm%mx, as_prm%mxy, &
			as_g%xlen_r(i),as_g%ylen_r(i), &
			as_g%x0r(i),as_g%y0r(i),as_g%z0r(i), &
                        0.0d0 &
                        ) 
          end if
       end do
    end do
    call parcel_sync(method)

    do i=1, 6
       call cudafree( sbuf(i)%v )
       call cudafree( rbuf(i)%v )
    end do

  end subroutine AS_util_copy_to_AS_GPU

  subroutine AS_util_copy_from_AS_GPU( &
       method, &  
       comm, myrank, ndist, prm, as_ndist, as_prm, &
       & as_diag_neighbor, as_mean, as_g,&
       wv4   , offset_wv4, &
       as_wv4, offset_as_wv4 &
       )
    use krylov,only : ksp
    use parcel_vec_kernel,only : parcel_sync
    use ddmmat_wrap,only : ddm_param,as_geom
    use c_wrapper,only : devicemem
    use cu_common_wrapper,only : &
         cudamalloc_fp64,cudafree, &
         c_mpi_isend_fp64,c_mpi_irecv_fp64,c_mpi_wait, &
         cudamemcpy_d2h_fp64,cudamemcpy_h2d_fp64

    use parcel_define,only : &
         ksp_additiveschwarz_type_interpolate, &
         ksp_additiveschwarz_type_none, &
         ksp_additiveschwarz_type_basic

    use mpi
    implicit none
    type(ksp) :: method
    !   include "ddm.h"
    !   include "mpif.h"
    integer        , intent(in ) :: comm      !
    integer        , intent(in ) :: myrank    !
    integer(kind=inttype) , intent(in ) :: ndist     ! 
    type(ddm_param), intent(in ) :: prm       ! 
    integer(kind=inttype) , intent(in ) :: as_ndist  ! 領域のグリッド数（as = 加法シュワルツ）
    type(ddm_param), intent(in ) :: as_prm    ! 領域形状
    integer        , intent(in ) :: as_diag_neighbor ! 斜め隣の領域の値を使うなら1、使わないなら0
    integer        , intent(in ) :: as_mean          ! 隣の領域との平均をとるなら1、とらないなら0
    type(as_geom)  , intent(in ) :: as_g             !
!    real*8         , intent(inout) :: wv4(ndist)       ! 
!    real*8         , intent(inout) :: as_wv4(as_ndist) ! 前処理行列の作用対象となるベクトル

    type(devicemem) wv4
    integer(kind=inttype) offset_wv4
    type(devicemem) as_wv4
    integer(kind=inttype) offset_as_wv4

    integer as_i, i, ix, iy, iz, j, k
    integer size, x0, y0, z0, xlen, ylen, zlen, info
    real*8  coeff
    integer mpx, mpy, mpz, xyz, pm, nsize
    type gpu_buff
       type(devicemem) v
    end type gpu_buff
    type(gpu_buff) sbuf(6), rbuf(6) ! 6 = ±x, ±y, ±z   

    type cpu_buff
       real*8, allocatable :: v(:)
    end type cpu_buff
    type(cpu_buff) sbuf_host(6), rbuf_host(6) ! 6 = ±x, ±y, ±z   

    integer :: ireqs(6), ireqr(6)
    integer istatus (mpi_status_size)
    integer(kind=intType),parameter :: offset_zero=0
    integer(kind=inttype) parcel_int

    !!| 平均をとらない場合は通信せずにコピーだけ
    if ( as_mean .eq. 0 ) goto 100

    !!| 加算する場合 
    ! 隣接プロセスからデータを通信により取得する
    !   as_util_copy_to_as の場合に対して send / receive の関係が逆になることに注意
    !   取得したデータをas_wv4に加算

    do i=1, 6
       parcel_int = as_g%grids_s(i)
       call cudamalloc_fp64(sbuf(i)%v, parcel_int ) 
       parcel_int = as_g%grids_r(i)
       call cudamalloc_fp64(rbuf(i)%v, parcel_int ) 
    end do

    call parcel_sync(method)
    ! 斜め隣りのプロセスの値を使用する場合
    ! 同じ軸方向の領域のみ同時に通信できる
    do xyz=2, 0, -1  ! x, y, z
       do pm=0, 1  ! +xyz, -xyz
          i = 2 * xyz + pm + 1
          ! 通信データの設定
          if ( as_g%rank_r(i) /= -1 ) then
             ! do j=0, as_g%grids_r(i)-1
             !    ix = as_g%x0r(i) + mod(j, as_g%xlen_r(i))
             !    iy = as_g%y0r(i) + mod(j/as_g%xlen_r(i), as_g%ylen_r(i))
             !    iz = as_g%z0r(i) + j / (as_g%xlen_r(i) * as_g%ylen_r(i))
             !    k = ix + as_prm%mx * iy + as_prm%mxy * iz + 1
             !    rbuf(i)%v(j+1) = as_wv4(k)
             ! end do
             call as_util_set_sbuf_fp64_gpu(&
			rbuf(i)%v , offset_zero, &
			as_wv4    , offset_as_wv4, &
			as_g%grids_r(i), &
			as_prm%mx, as_prm%mxy, &
			as_g%xlen_r(i),as_g%ylen_r(i), &
			as_g%x0r(i),as_g%y0r(i),as_g%z0r(i)  &
                        ) 
          end if
       enddo
    enddo
    call parcel_sync(method)

#if CUDA_AWARE_
    do xyz=2, 0, -1  ! x, y, z
       do pm=0, 1  ! +xyz, -xyz
          i = 2 * xyz + pm + 1
          ! 通信開始
          if ( as_g%rank_r(i) /= -1 ) then
             ! call mpi_isend( rbuf(i)%v, size(rbuf(i)%v), mpi_double_precision, &
             !      as_g%rank_r(i), as_g%rank_r(i), comm, ireqr(i), info )
             call c_mpi_isend_fp64(&
                  rbuf(i)%v, &
                  as_g%grids_r(i), &
                  as_g%rank_r(i), &
                  as_g%rank_r(i), &
                  method%commvec%c_comm, &
                  method%commvec%c_r_req, &
                  i,offset_zero &
                  ) 
          end if
          if ( as_g%rank_s(i) /= -1 ) then
             ! call mpi_irecv( sbuf(i)%v, size(sbuf(i)%v), mpi_double_precision, &
             !      as_g%rank_s(i), myrank, comm, ireqs(i), info )
             call c_mpi_irecv_fp64(&
                  rbuf(i)%v, &
                  as_g%grids_s(i), &
                  as_g%rank_s(i), &
                  myrank, &
                  method%commvec%c_comm, &
                  method%commvec%c_s_req, &
                  i,offset_zero &
                  ) 

          end if
       end do
    end do

    call parcel_sync(method)
    ! 通信を完了させて、受信したデータを加算
    do xyz=2, 0, -1  ! x, y, z
       do pm=0, 1  ! +xyz, -xyz
          i = 2 * xyz + pm + 1
          ! if ( as_g%rank_s(i) /= -1 ) call mpi_wait( ireqs(i), istatus, info )
          ! if ( as_g%rank_r(i) /= -1 ) call mpi_wait( ireqr(i), istatus, info )
          if ( as_g%rank_s(i) /= -1 ) call c_mpi_wait(method%commvec%c_s_req,method%commvec%c_s_sta,i) 
          if ( as_g%rank_r(i) /= -1 ) call c_mpi_wait(method%commvec%c_r_req,method%commvec%c_r_sta,i) 
       enddo
    enddo
#else

    do i=1, 6
       allocate( sbuf_host(i)%v(AS_g%grids_s(i)) )
       allocate( rbuf_host(i)%v(AS_g%grids_r(i)) )
    end do

    do xyz=0, 2  ! x, y, z
       do pm=0, 1  ! +xyz, -xyz
          i = 2 * xyz + pm + 1
          parcel_int = AS_g%grids_r(i)
          call cudaMemcpy_d2h_FP64(&
               rbuf(i)%v ,rbuf_host(i)%v, parcel_int,&
               offset_zero, &
               offset_zero  &
               ) 
       enddo
    enddo

    call parcel_sync(method)
    do xyz=2, 0, -1  ! x, y, z
       do pm=0, 1  ! +xyz, -xyz
          i = 2 * xyz + pm + 1
          ! 通信開始
          if ( as_g%rank_r(i) /= -1 ) then
             call mpi_isend( rbuf_host(i)%v, &
                  size(rbuf_host(i)%v), &
                  mpi_double_precision, &
                  as_g%rank_r(i), &
                  as_g%rank_r(i), &
                  comm, &
                  ireqr(i), &
                  info )
          end if
          if ( as_g%rank_s(i) /= -1 ) then
             call mpi_irecv( sbuf_host(i)%v, &
                  size(sbuf_host(i)%v), &
                  mpi_double_precision, &
                  as_g%rank_s(i), &
                  myrank, &
                  comm, &
                  ireqs(i), &
                  info )

          end if
       end do
    end do

    call parcel_sync(method)
    ! 通信を完了させて、受信したデータを加算
    do xyz=2, 0, -1  ! x, y, z
       do pm=0, 1  ! +xyz, -xyz
          i = 2 * xyz + pm + 1
          if ( as_g%rank_s(i) /= -1 ) call mpi_wait( ireqs(i), istatus, info )
          if ( as_g%rank_r(i) /= -1 ) call mpi_wait( ireqr(i), istatus, info )
       enddo
    enddo

    do i=1, 6
       parcel_int = AS_g%grids_s(i)
       call cudaMemcpy_h2d_FP64(&
         sbuf_host(i)%v, &
         sbuf(i)%v, parcel_int,&
         offset_zero, &
         offset_zero) 
    enddo
    call parcel_sync(method)
    do i=1, 6
       deallocate( sbuf_host(i)%v )
       deallocate( rbuf_host(i)%v )
    end do


#endif

    call parcel_sync(method)
    do xyz=2, 0, -1  ! x, y, z
       do pm=0, 1  ! +xyz, -xyz
          i = 2 * xyz + pm + 1
          ! 受信したデータを加算
          if ( as_g%rank_s(i) /= -1 ) then
             ! do j=0, as_g%grids_s(i)-1
             !    ix = as_g%x0s(i) + mod(j, as_g%xlen_s(i))
             !    iy = as_g%y0s(i) + mod(j/as_g%xlen_s(i),  as_g%ylen_s(i))
             !    iz = as_g%z0s(i) + j / (as_g%xlen_s(i) * as_g%ylen_s(i))
             !    k = ix + as_prm%mx * iy + as_prm%mxy * iz + 1
             !    as_wv4(k) = as_wv4(k) + sbuf(i)%v(j+1)
             ! end do
             call AS_util_set_rbuf_FP64_GPU(&
                  sbuf(i)%v , offset_zero, &
                  as_wv4, offset_as_wv4, &
                  as_g%grids_s(i), &
                  as_prm%mx, as_prm%mxy, &
                  as_g%xlen_s(i),as_g%ylen_s(i), &
                  as_g%x0s(i),as_g%y0s(i),as_g%z0s(i), &
                  1.0d0 &
                  ) 
          end if
       end do
    end do

    call parcel_sync(method)

    do i=1, 6
       call cudafree( sbuf(i)%v )
       call cudafree( rbuf(i)%v )
    end do

100 continue

    call util_AS_copy_only_FP64_GPU(&
         wv4, offset_wv4, &
         AS_wv4, offset_AS_wv4, &
         prm%nx,prm%ny,prm%nz, &
         prm%mx,prm%mxy, &
         prm%stmx,prm%stmy,prm%stmz, &
         ndist, &
         as_prm%mx,as_prm%mxy, &
         as_g%sleeve_mx,as_g%sleeve_my,as_g%sleeve_mz &
         )
    call parcel_sync(method)
    
  end subroutine AS_util_copy_from_AS_GPU

  ! -------------------------------------

  subroutine solve_DILU_DDM(method,r,s)
    use krylov
    implicit none
    type(KSP) :: method

    real*8 r(method%ddmmat%n)
    real*8 s(method%ddmmat%n)

    real*8 zz

    integer nx,ny,nz
    integer stmx, stmy, stmz
    integer m
    integer mx,my,mz
    integer mxy

    integer nblock
    integer nxblock,nyblock,nzblock
    integer nx_local,ny_local,nz_local
    integer m_local
    integer mx_local,my_local,mz_local

    integer nnd,offsetj,offsetjj
    integer offsetx,offsety,offsetz

    integer j,jj
    integer jx,jy,jz
    integer jjx,jjy,jjz
    integer block_id
    integer jxblock,jyblock,jzblock
    integer ip,ipp,ip_diag

    integer nxs,nxe
    integer nys,nye
    integer nzs,nze
    integer offset

    nx = method%ddmMat%prm%nx
    ny = method%ddmMat%prm%ny
    nz = method%ddmMat%prm%nz

    stmx = method%ddmMat%prm%stmx
    stmy = method%ddmMat%prm%stmy
    stmz = method%ddmMat%prm%stmz

    mx  = method%ddmMat%prm%mx
    my  = method%ddmMat%prm%my
    mz  = method%ddmMat%prm%mz
    mxy = method%ddmMat%prm%mxy
    m   =  mx * my * mz

    nxblock = method%ddmMat%DDMMat_DILU%nxblock
    nyblock = method%ddmMat%DDMMat_DILU%nyblock
    nzblock = method%ddmMat%DDMMat_DILU%nzblock
    nblock  = nxblock*nyblock*nzblock;

    nx_local = method%ddmMat%DDMMat_DILU%nx_local
    ny_local = method%ddmMat%DDMMat_DILU%ny_local
    nz_local = method%ddmMat%DDMMat_DILU%nz_local

    mx_local =    method%ddmMat%DDMMat_DILU%mx_local
    my_local =    method%ddmMat%DDMMat_DILU%my_local
    mz_local =    method%ddmMat%DDMMat_DILU%mz_local
    m_local  =    method%ddmMat%DDMMat_DILU%m_local

    nnd = method%ddmMat%nnd

    ! -----
    ! (D+L)D^-1(D+U)s=r
    ! -----
    !  1. (D+L)z=r       
    !     z = Dinv * (r - Lz)
    !  2. (I+D^-1U)s=z
    !     s = z-D^-1 Us       

    !$OMP PARALLEL DO PRIVATE(&
    !$OMP jzblock,jyblock,jxblock, &
    !$OMP nzs,nze, &
    !$OMP nys,nye, &
    !$OMP nxs,nxe, &
    !$OMP block_id, &
    !$OMP j,jx,jy,jz, &
    !$OMP jj,jjx,jjy,jjz, &
    !$OMP zz, &
    !$OMP ip, &
    !$OMP offset,offsetjj, &
    !$OMP offsetx,offsety,offsetz &
    !$OMP ) &
    !$OMP COLLAPSE(3)
    do jzblock=1,nzblock
       do jyblock=1,nyblock
          do jxblock=1,nxblock
             nzs = 1 + nz_local*(jzblock-1)
             nze = method%ddmMat%DDMMat_DILU%nze(jzblock)
             nys = 1 + ny_local*(jyblock-1)
             nye = method%ddmMat%DDMMat_DILU%nye(jyblock)
             nxs = 1 + nx_local*(jxblock-1)
             nxe = method%ddmMat%DDMMat_DILU%nxe(jxblock)
             block_id = (jxblock-1) + (jyblock-1)*nxblock + (jzblock-1)*nxblock*nyblock 

#if 1
             ! 前進代入
             do jz = nzs-1,nze-1
                do jy = nys-1,nye-1
                   do jx = nxs-1,nxe-1
                      j  = (jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1

                      jjx = jx - (nxs-1)
                      jjy = jy - (nys-1)
                      jjz = jz - (nzs-1)
                      jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1

                      zz = r(j)
                      do ip = 1,nnd
                         offset   = method%ddmMat%mat%ioff_dia(ip)
                         offsetx  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*0)
                         offsety  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*1)
                         offsetz  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*2)
                         offsetjj = offsetx + offsety*mx_local + offsetz*mx_local*my_local
                         if(offset>=0)then
                            cycle
                         endif
                         zz = zz  & 
                              - method%ddmMat%DDMMat_DILU%svec(block_id + (jj+offsetjj)*nblock) * &
                                method%ddmMat%mat%val_dia(m*(ip-1)+j)
                      enddo
                      method%ddmMat%DDMMat_DILU%svec(block_id + jj*nblock) = &
                           method%ddmMat%DDMMat_DILU%Dinv(block_id + jj*nblock)*zz
                   enddo
                enddo
             enddo

             ! 後退代入
             do jz = nze-1,nzs-1,-1
                do jy = nye-1,nys-1,-1
                   do jx = nxe-1,nxs-1,-1

                      j  = (jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1

                      jjx = jx - (nxs-1)
                      jjy = jy - (nys-1)
                      jjz = jz - (nzs-1)
                      jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1

                      zz = 0.0d0
                      do ip = 1,nnd
                         offset = method%ddmMat%mat%ioff_dia(ip)
                         offsetx  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*0)
                         offsety  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*1)
                         offsetz  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*2)
                         offsetjj = offsetx + offsety*mx_local + offsetz*mx_local*my_local
                         if(offset<=0)then
                            cycle
                         endif
                         zz = zz  & 
                              + method%ddmMat%DDMMat_DILU%svec(block_id + (jj+offsetjj)*nblock) * &
                                method%ddmMat%mat%val_dia(m*(ip-1)+j)
                      enddo
                      method%ddmMat%DDMMat_DILU%svec(block_id + jj*nblock) =  &
                           method%ddmMat%DDMMat_DILU%svec(block_id + jj*nblock) - &
                           method%ddmMat%DDMMat_DILU%Dinv(block_id + jj*nblock)*zz
                      s(j) = method%ddmMat%DDMMat_DILU%svec(block_id + jj*nblock) 
                   enddo
                enddo
             enddo


#else
             ! 前進代入
             do jz = nzs-1,nze-1
                do jy = nys-1,nye-1
                   do jx = nxs-1,nxe-1
                      j  = (jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1

                      jjx = jx - (nxs-1)
                      jjy = jy - (nys-1)
                      jjz = jz - (nzs-1)
                      jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1

                      zz = r(j)
                      do ip = 1,nnd
                         offset   = method%ddmMat%mat%ioff_dia(ip)
                         offsetx  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*0)
                         offsety  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*1)
                         offsetz  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*2)
                         offsetjj = offsetx + offsety*mx_local + offsetz*mx_local*my_local
                         if(offset>=0)then
                            cycle
                         endif
                         zz = zz  & 
                              - method%ddmMat%DDMMat_DILU%svec(jj + m_local * block_id+offsetjj) * &
                                method%ddmMat%mat%val_dia(m*(ip-1)+j)
                      enddo
                      method%ddmMat%DDMMat_DILU%svec(jj + m_local * block_id) = &
                           method%ddmMat%DDMMat_DILU%Dinv(jj + m_local * block_id)*zz                      
                   enddo
                enddo
             enddo

             ! 後退代入
             do jz = nze-1,nzs-1,-1
                do jy = nye-1,nys-1,-1
                   do jx = nxe-1,nxs-1,-1

                      j  = (jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1

                      jjx = jx - (nxs-1)
                      jjy = jy - (nys-1)
                      jjz = jz - (nzs-1)
                      jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1

                      zz = 0.0d0
                      do ip = 1,nnd
                         offset = method%ddmMat%mat%ioff_dia(ip)
                         offsetx  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*0)
                         offsety  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*1)
                         offsetz  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*2)
                         offsetjj = offsetx + offsety*mx_local + offsetz*mx_local*my_local
                         if(offset<=0)then
                            cycle
                         endif
                         zz = zz  & 
                              + method%ddmMat%DDMMat_DILU%svec(jj + m_local * block_id+offsetjj) * &
                                method%ddmMat%mat%val_dia(m*(ip-1)+j)
                      enddo
                      method%ddmMat%DDMMat_DILU%svec(jj + m_local * block_id) =  &
                           method%ddmMat%DDMMat_DILU%svec(jj + m_local * block_id) - &
                           method%ddmMat%DDMMat_DILU%Dinv(jj + m_local * block_id)*zz
                      s(j) = method%ddmMat%DDMMat_DILU%svec(jj + m_local * block_id) 
                   enddo
                enddo
             enddo
#endif

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO 

  end subroutine solve_DILU_DDM

  subroutine get_DILU_DDM(method)
    use krylov
    implicit none

    type(KSP) :: method
    integer nx,ny,nz
    integer stmx, stmy, stmz
    integer m
    integer mx,my,mz
    integer mxy

    integer lattice_dimention
    integer nblock
    integer nxblock,nyblock,nzblock
    integer nx_local,ny_local,nz_local
    integer m_local
    integer mx_local,my_local,mz_local

    integer nnd,offseti,offsetj,offsetii,offsetjj
    integer offsetx,offsety,offsetz

    integer j,jj
    integer jx,jy,jz
    integer jjx,jjy,jjz
    integer block_id
    integer jxblock,jyblock,jzblock
    integer ip,ipp,ip_diag

    integer nxs,nxe
    integer nys,nye
    integer nzs,nze
    real*8 Dinv 

    nx = method%ddmMat%prm%nx
    ny = method%ddmMat%prm%ny
    nz = method%ddmMat%prm%nz

    stmx = method%ddmMat%prm%stmx
    stmy = method%ddmMat%prm%stmy
    stmz = method%ddmMat%prm%stmz

    mx  = method%ddmMat%prm%mx
    my  = method%ddmMat%prm%my
    mz  = method%ddmMat%prm%mz
    mxy = method%ddmMat%prm%mxy
    m   =  mx * my * mz

    nxblock = method%ddmMat%DDMMat_DILU%nxblock
    nyblock = method%ddmMat%DDMMat_DILU%nyblock
    nzblock = method%ddmMat%DDMMat_DILU%nzblock
    nblock  = nxblock*nyblock*nzblock;

    lattice_dimention = method%ddmMat%DDMMat_DILU%lattice_dimention 
    nx_local = method%ddmMat%DDMMat_DILU%nx_local
    ny_local = method%ddmMat%DDMMat_DILU%ny_local
    nz_local = method%ddmMat%DDMMat_DILU%nz_local

    mx_local =    method%ddmMat%DDMMat_DILU%mx_local
    my_local =    method%ddmMat%DDMMat_DILU%my_local
    mz_local =    method%ddmMat%DDMMat_DILU%mz_local
    m_local =    method%ddmMat%DDMMat_DILU%m_local

    nnd = method%ddmMat%nnd

    ip_diag = -1
    do ip = 1,nnd
       if(method%ddmMat%mat%ioff_dia(ip).eq.0)then
          ip_diag = ip
       endif
    enddo

    !$OMP PARALLEL DO PRIVATE(&
    !$OMP jzblock,jyblock,jxblock, &
    !$OMP nzs,nze, &
    !$OMP nys,nye, &
    !$OMP nxs,nxe, &
    !$OMP block_id, &
    !$OMP j,jx,jy,jz, &
    !$OMP jj,jjx,jjy,jjz, &
    !$OMP ip,ipp, &
    !$OMP offsetj,offsetjj, &
    !$OMP offseti,offsetii, &
    !$OMP offsetx,offsety,offsetz &
    !$OMP ) &
    !$OMP COLLAPSE(3)
    do jzblock=1,nzblock
       do jyblock=1,nyblock
          do jxblock=1,nxblock
             nzs = 1 + nz_local*(jzblock-1)
             nze = method%ddmMat%DDMMat_DILU%nze(jzblock)
             nys = 1 + ny_local*(jyblock-1)
             nye = method%ddmMat%DDMMat_DILU%nye(jyblock)
             nxs = 1 + nx_local*(jxblock-1)
             nxe = method%ddmMat%DDMMat_DILU%nxe(jxblock)
             block_id = (jxblock-1) + (jyblock-1)*nxblock + (jzblock-1)*nxblock*nyblock 

#if 1
             do jz = nzs-1,nze-1
                do jy = nys-1,nye-1
                   do jx = nxs-1,nxe-1
                      j  = (jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1

                      jjx = jx - (nxs-1)
                      jjy = jy - (nys-1)
                      jjz = jz - (nzs-1)
                      jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1                      
                      method%ddmMat%DDMMat_DILU%Dinv(block_id + jj*nblock) = &
                           method%ddmMat%mat%val_dia(m*(ip_diag-1)+j) 
                   enddo
                enddo
             enddo

             do jz = nzs-1,nze-1
                do jy = nys-1,nye-1
                   do jx = nxs-1,nxe-1
                      j = ( jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1

                      jjx = jx - (nxs-1)
                      jjy = jy - (nys-1)
                      jjz = jz - (nzs-1)
                      jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1

                      ! -----
                      ! D-ILU分解
                      ! A = AL + AD + AU ~= K
                      ! K = (D+AL)D^-1(D+AU)
                      ! 条件1 対角成分が一致     :: D(j,j) = A(j,j) - Sum[i<j]          A(j,i) D(i,i)^-1 A(i,j)
                      ! 条件2 行内の要素和が一致 :: D(j,j) = A(j,j) - Sum[i<j] Sum[k>i] A(j,i) D(i,i)^-1 A(i,k)
                      ! -----
                      do ip = 1,nnd
                         offsetj  = method%ddmMat%mat%ioff_dia(ip)
                         offsetx  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*0)
                         offsety  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*1)
                         offsetz  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*2)
                         offsetjj = offsetx + offsety*mx_local + offsetz*mx_local*my_local

                         ! A(j,i) : method%ddmMat%mat%val_dia(j+m*(ip-1))
                         do ipp = 1,nnd
                            offseti  = method%ddmMat%mat%ioff_dia(ipp)
                            offsetx  = method%ddmMat%DDMMat_DILU%offset_xyz(ipp+nnd*0)
                            offsety  = method%ddmMat%DDMMat_DILU%offset_xyz(ipp+nnd*1)
                            offsetz  = method%ddmMat%DDMMat_DILU%offset_xyz(ipp+nnd*2)
                            offsetii = offsetx + offsety*mx_local + offsetz*mx_local*my_local 
                            if( jj <= jj+offsetii)then
                               cycle
                            endif
                            ! D(j,j) = A(j,j) - Sum[i<j] A(j,i) D(i,i)^-1 A(i,j) 
                            if( offseti .eq. (-offsetj) )then
                               method%ddmMat%DDMMat_DILU%Dinv(jj*nblock + block_id) = &
                                    method%ddmMat%DDMMat_DILU%Dinv( (jj           )*nblock + block_id) - &
                                    method%ddmMat%DDMMat_DILU%Dinv( (jj + offsetii)*nblock + block_id) * &
                                    method%ddmMat%mat%val_dia(j +           m*(ipp-1)  )                * &
                                    method%ddmMat%mat%val_dia(j + offseti + m*(ip-1) )
                            endif
                         enddo
                      enddo
                      if(abs(method%ddmMat%DDMMat_DILU%Dinv(block_id + jj*nblock)).eq.0.0d0)then
                         Dinv = 0.0d0
                      else
                         Dinv = 1.0d0 / method%ddmMat%DDMMat_DILU%Dinv(block_id + jj*nblock)
                      endif
                      method%ddmMat%DDMMat_DILU%Dinv(block_id + jj*nblock) = Dinv

                   enddo
                enddo
             enddo

#else
             do jz = nzs-1,nze-1
                do jy = nys-1,nye-1
                   do jx = nxs-1,nxe-1
                      j  = (jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1

                      jjx = jx - (nxs-1)
                      jjy = jy - (nys-1)
                      jjz = jz - (nzs-1)
                      jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1                      
                      method%ddmMat%DDMMat_DILU%Dinv(jj + m_local * block_id) = &
                           method%ddmMat%mat%val_dia(m*(ip_diag-1)+j) 
                   enddo
                enddo
             enddo

             do jz = nzs-1,nze-1
                do jy = nys-1,nye-1
                   do jx = nxs-1,nxe-1
                      j = ( jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1

                      jjx = jx - (nxs-1)
                      jjy = jy - (nys-1)
                      jjz = jz - (nzs-1)
                      jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1

                      ! -----
                      ! D-ILU分解
                      ! A = AL + AD + AU ~= K
                      ! K = (D+AL)D^-1(D+AU)
                      ! 条件1 対角成分が一致     :: D(j,j) = A(j,j) - Sum[i<j]          A(j,i) D(i,i)^-1 A(i,j)
                      ! 条件2 行内の要素和が一致 :: D(j,j) = A(j,j) - Sum[i<j] Sum[k>i] A(j,i) D(i,i)^-1 A(i,k)
                      ! -----
                      do ip = 1,nnd
                         offsetj  = method%ddmMat%mat%ioff_dia(ip)
                         offsetx  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*0)
                         offsety  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*1)
                         offsetz  = method%ddmMat%DDMMat_DILU%offset_xyz(ip+nnd*2)
                         offsetjj = offsetx + offsety*mx_local + offsetz*mx_local*my_local

                         ! A(j,i) : method%ddmMat%mat%val_dia(j+m*(ip-1))
                         do ipp = 1,nnd
                            offseti  = method%ddmMat%mat%ioff_dia(ipp)
                            offsetx  = method%ddmMat%DDMMat_DILU%offset_xyz(ipp+nnd*0)
                            offsety  = method%ddmMat%DDMMat_DILU%offset_xyz(ipp+nnd*1)
                            offsetz  = method%ddmMat%DDMMat_DILU%offset_xyz(ipp+nnd*2)
                            offsetii = offsetx + offsety*mx_local + offsetz*mx_local*my_local 
                            if( jj <= jj+offsetii)then
                               cycle
                            endif
                            ! D(j,j) = A(j,j) - Sum[i<j] A(j,i) D(i,i)^-1 A(i,j) 
                            if( offseti .eq. (-offsetj) )then
                               method%ddmMat%DDMMat_DILU%Dinv(jj + m_local * block_id) = &
                                    method%ddmMat%DDMMat_DILU%Dinv(jj            + m_local * block_id) - &
                                    method%ddmMat%DDMMat_DILU%Dinv(jj + offsetii + m_local * block_id) * &
                                    method%ddmMat%mat%val_dia(j +           m*(ipp-1)  )                * &
                                    method%ddmMat%mat%val_dia(j + offseti + m*(ip-1) )
                            endif
                         enddo
                      enddo
                      if(abs(method%ddmMat%DDMMat_DILU%Dinv(jj + m_local * block_id)).eq.0.0d0)then
                         Dinv = 0.0d0
                      else
                         Dinv = 1.0d0 / method%ddmMat%DDMMat_DILU%Dinv(jj + m_local * block_id)
                      endif
                      method%ddmMat%DDMMat_DILU%Dinv(jj + m_local * block_id) = Dinv

                   enddo
                enddo
             enddo
#endif
          enddo
       enddo
    end do

#if 0
    do jzblock=1,nzblock
       do jyblock=1,nyblock
          do jxblock=1,nxblock
             nzs = 1 + nz_local*(jzblock-1)
             nze = method%ddmMat%DDMMat_DILU%nze(jzblock)
             nys = 1 + ny_local*(jyblock-1)
             nye = method%ddmMat%DDMMat_DILU%nye(jyblock)
             nxs = 1 + nx_local*(jxblock-1)
             nxe = method%ddmMat%DDMMat_DILU%nxe(jxblock)

             block_id = (jxblock-1) + (jyblock-1)*nxblock + (jzblock-1)*nxblock*nyblock 

             do jz = nzs-1,nze-1
                do jy = nys-1,nye-1
                   do jx = nxs-1,nxe-1
                      j  = (jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1

                      jjx = jx - (nxs-1)
                      jjy = jy - (nys-1)
                      jjz = jz - (nzs-1)
                      jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1                      
                      write(*,*)" BJ : ",jx,jy,jz,method%ddmMat%DDMMat_DILU%Dinv(block_id + jj*nblock)
                   enddo
                enddo
             enddo

          enddo
       enddo
    end do

#endif

  end subroutine get_DILU_DDM

  subroutine init_DILU_DDM( &
       method,&
       xdiv,ydiv,zdiv &
       )
    use krylov
    implicit none
    type(KSP) :: method
    integer xdiv,ydiv,zdiv 

    integer nnd
    integer nx,ny,nz
    integer stmx, stmy, stmz

    integer lattice_dimention
    integer nx_local,ny_local,nz_local
    integer m_local
    integer mx_local,my_local,mz_local
    integer nblock
    integer nxblock,nyblock,nzblock

    integer i
    integer nxs,nxe
    integer nys,nye
    integer nzs,nze

    integer j
    integer jx,jy,jz
    integer jxblock,jyblock,jzblock

    xdiv = max(xdiv,1)
    ydiv = max(ydiv,1)
    zdiv = max(zdiv,1)

    lattice_dimention = method%ddmMat%prm%ncomp_grids

    nnd = method%ddmMat%nnd
    nx = method%ddmMat%prm%nx
    ny = method%ddmMat%prm%ny
    nz = method%ddmMat%prm%nz

    stmx = method%ddmMat%prm%stmx
    stmy = method%ddmMat%prm%stmy
    stmz = method%ddmMat%prm%stmz

    ! ------------
    if(mod(nx,xdiv).ne.0)then
       nx_local = nx/xdiv + 1
    else
       nx_local = nx/xdiv
    endif
    allocate(method%ddmMat%DDMMat_DILU%nxe(xdiv))
    do i = 1,xdiv
       nxs = 1 + nx_local*(i-1)
       nxe = nx_local*i
       method%ddmMat%DDMMat_DILU%nxe(i) = min(nx,nxe)
    enddo
    mx_local = nx_local + stmx*2

    ! ------------
    if(mod(ny,ydiv).ne.0)then
       ny_local = ny/ydiv + 1
    else
       ny_local = ny/ydiv
    endif
    allocate(method%ddmMat%DDMMat_DILU%nye(ydiv))
    do i = 1,ydiv
       nys = 1 + ny_local*(i-1)
       nye = ny_local*i
       method%ddmMat%DDMMat_DILU%nye(i) = min(ny,nye)
    enddo
    my_local = ny_local + stmy*2

    ! ------------
    if(mod(nz,zdiv).ne.0)then
       nz_local = nz/zdiv + 1
    else
       nz_local = nz/zdiv
    endif
    allocate(method%ddmMat%DDMMat_DILU%nze(zdiv))
    do i = 1,zdiv
       nzs = 1 + nz_local*(i-1)
       nze = nz_local*i
       method%ddmMat%DDMMat_DILU%nze(i) = min(nz,nze)
    enddo
    mz_local = nz_local + stmz*2

    ! ------------

    m_local = mx_local * my_local * mz_local 

    nxblock = xdiv
    nyblock = ydiv
    nzblock = zdiv

    nblock  = nxblock * nyblock * nzblock
    method%ddmMat%DDMMat_DILU%lattice_dimention = lattice_dimention
    method%ddmMat%DDMMat_DILU%nxblock = nxblock
    method%ddmMat%DDMMat_DILU%nyblock = nyblock
    method%ddmMat%DDMMat_DILU%nzblock = nzblock

    method%ddmMat%DDMMat_DILU%nx_local = nx_local
    method%ddmMat%DDMMat_DILU%ny_local = ny_local
    method%ddmMat%DDMMat_DILU%nz_local = nz_local

    method%ddmMat%DDMMat_DILU%mx_local = mx_local
    method%ddmMat%DDMMat_DILU%my_local = my_local
    method%ddmMat%DDMMat_DILU%mz_local = mz_local
    method%ddmMat%DDMMat_DILU%m_local  = m_local 

    allocate(method%ddmMat%DDMMat_DILU%Dinv(m_local * nblock))
    allocate(method%ddmMat%DDMMat_DILU%svec(m_local * nblock))
    allocate(method%ddmMat%DDMMat_DILU%rvec(m_local * nblock))

    !$OMP PARALLEL DO
    do i = 1,m_local * nblock
       method%ddmMat%DDMMat_DILU%Dinv(i) = 0.0d0
       method%ddmMat%DDMMat_DILU%svec(i) = 0.0d0
       method%ddmMat%DDMMat_DILU%rvec(i) = 0.0d0
    enddo
    !$OMP END PARALLEL DO

    allocate(method%ddmMat%DDMMat_DILU%offset_xyz(nnd*3))
    method%ddmMat%DDMMat_DILU%offset_xyz = 0 

    select case (lattice_dimention)
    case(1)
       do i = 1,nnd
          method%ddmMat%DDMMat_DILU%offset_xyz(i+nnd*0) = method%ddmMat%mat%ioff_grids(i,1)
       enddo
    case(2)
       do i = 1,nnd
          method%ddmMat%DDMMat_DILU%offset_xyz(i+nnd*0) = method%ddmMat%mat%ioff_grids(i,1)
          method%ddmMat%DDMMat_DILU%offset_xyz(i+nnd*1) = method%ddmMat%mat%ioff_grids(i,2)
       enddo
    case(3)
       do i = 1,nnd
          method%ddmMat%DDMMat_DILU%offset_xyz(i+nnd*0) = method%ddmMat%mat%ioff_grids(i,1)
          method%ddmMat%DDMMat_DILU%offset_xyz(i+nnd*1) = method%ddmMat%mat%ioff_grids(i,2)
          method%ddmMat%DDMMat_DILU%offset_xyz(i+nnd*2) = method%ddmMat%mat%ioff_grids(i,3)
       enddo
    end select

    ! write(*,*)"nxblock  = ",nxblock
    ! write(*,*)"nyblock  = ",nyblock
    ! write(*,*)"nzblock  = ",nzblock

    ! write(*,*)"mx_local = ",mx_local
    ! write(*,*)"my_local = ",my_local
    ! write(*,*)"mz_local = ",mz_local

#if 0
    do jxblock = 1,size(method%ddmMat%DDMMat_DILU%nxe)
       nxs = 1 + nx_local*(jxblock-1)
       nxe = method%ddmMat%DDMMat_DILU%nxe(jxblock)
       write(*,*)" jxloop : ",jxblock,nxs,nxe
    enddo
    do jyblock = 1,size(method%ddmMat%DDMMat_DILU%nye)
       nys = 1 + ny_local*(jyblock-1)
       nye = method%ddmMat%DDMMat_DILU%nye(jyblock)
       write(*,*)" jyloop : ",jyblock,nys,nye
    enddo
    do jzblock = 1,size(method%ddmMat%DDMMat_DILU%nze)
       nzs = 1 + nz_local*(jzblock-1)
       nze = method%ddmMat%DDMMat_DILU%nze(jzblock)
       write(*,*)" jzloop : ",jzblock,nzs,nze
    enddo

    if(method%commVec%myrank.eq.0)then
       write(*,*)"DDM BlockJACOBI local size : nx,ny,nz = ", &
            nx_local,ny_local,nz_local
    endif

#endif

  end subroutine init_DILU_DDM

  subroutine free_DILU_DDM(method)
    use krylov
    implicit none
    type(KSP) :: method

    if(allocated(method%ddmMat%DDMMat_DILU%nxs))deallocate(method%ddmMat%DDMMat_DILU%nxs)
    if(allocated(method%ddmMat%DDMMat_DILU%nys))deallocate(method%ddmMat%DDMMat_DILU%nys)
    if(allocated(method%ddmMat%DDMMat_DILU%nzs))deallocate(method%ddmMat%DDMMat_DILU%nzs)

    if(allocated(method%ddmMat%DDMMat_DILU%nxe))deallocate(method%ddmMat%DDMMat_DILU%nxe)
    if(allocated(method%ddmMat%DDMMat_DILU%nye))deallocate(method%ddmMat%DDMMat_DILU%nye)
    if(allocated(method%ddmMat%DDMMat_DILU%nze))deallocate(method%ddmMat%DDMMat_DILU%nze)

    if(allocated(method%ddmMat%DDMMat_DILU%nxs0))deallocate(method%ddmMat%DDMMat_DILU%nxs0)
    if(allocated(method%ddmMat%DDMMat_DILU%nys0))deallocate(method%ddmMat%DDMMat_DILU%nys0)
    if(allocated(method%ddmMat%DDMMat_DILU%nzs0))deallocate(method%ddmMat%DDMMat_DILU%nzs0)

    if(allocated(method%ddmMat%DDMMat_DILU%nxe0))deallocate(method%ddmMat%DDMMat_DILU%nxe0)
    if(allocated(method%ddmMat%DDMMat_DILU%nye0))deallocate(method%ddmMat%DDMMat_DILU%nye0)
    if(allocated(method%ddmMat%DDMMat_DILU%nze0))deallocate(method%ddmMat%DDMMat_DILU%nze0)

    if(allocated(method%ddmMat%DDMMat_DILU%Dinv))deallocate(method%ddmMat%DDMMat_DILU%Dinv)
    if(allocated(method%ddmMat%DDMMat_DILU%svec))deallocate(method%ddmMat%DDMMat_DILU%svec)
    if(allocated(method%ddmMat%DDMMat_DILU%rvec))deallocate(method%ddmMat%DDMMat_DILU%rvec)
    if(allocated(method%ddmMat%DDMMat_DILU%svec_AS))deallocate(method%ddmMat%DDMMat_DILU%svec_AS)
    if(allocated(method%ddmMat%DDMMat_DILU%offset_xyz))deallocate(method%ddmMat%DDMMat_DILU%offset_xyz)

  end subroutine free_DILU_DDM

  subroutine init_DILU_DDM_GPU(method)
    use krylov
    implicit none

    type(KSP) :: method
    integer(kind=intType) n
    integer(kind=intType),parameter :: offset_zero=0

    method%ddmMat%C_DDMMat_DILU%lattice_dimention = method%ddmMat%DDMMat_DILU%lattice_dimention
    method%ddmMat%C_DDMMat_DILU%nx_local = method%ddmMat%DDMMat_DILU%nx_local
    method%ddmMat%C_DDMMat_DILU%ny_local = method%ddmMat%DDMMat_DILU%ny_local
    method%ddmMat%C_DDMMat_DILU%nz_local = method%ddmMat%DDMMat_DILU%nz_local
    method%ddmMat%C_DDMMat_DILU%mx_local = method%ddmMat%DDMMat_DILU%mx_local
    method%ddmMat%C_DDMMat_DILU%my_local = method%ddmMat%DDMMat_DILU%my_local
    method%ddmMat%C_DDMMat_DILU%mz_local = method%ddmMat%DDMMat_DILU%mz_local
    method%ddmMat%C_DDMMat_DILU%m_local  = method%ddmMat%DDMMat_DILU%m_local
    method%ddmMat%C_DDMMat_DILU%nblock   = method%ddmMat%DDMMat_DILU%nblock
    method%ddmMat%C_DDMMat_DILU%nxblock  = method%ddmMat%DDMMat_DILU%nxblock
    method%ddmMat%C_DDMMat_DILU%nyblock  = method%ddmMat%DDMMat_DILU%nyblock
    method%ddmMat%C_DDMMat_DILU%nzblock  = method%ddmMat%DDMMat_DILU%nzblock

    if(allocated(method%ddmMat%DDMMat_DILU%nxe))then
       n = size(method%ddmMat%DDMMat_DILU%nxe)
       call cudaMalloc_int(method%ddmMat%C_DDMMat_DILU%nxe,n)
       call cudaMemcpy_h2d_int(method%ddmMat%DDMMat_DILU%nxe,method%ddmMat%C_DDMMat_DILU%nxe,n,offset_zero,offset_zero)
    endif

    if(allocated(method%ddmMat%DDMMat_DILU%nye))then
       n = size(method%ddmMat%DDMMat_DILU%nye)
       call cudaMalloc_int(method%ddmMat%C_DDMMat_DILU%nye,n)
       call cudaMemcpy_h2d_int(method%ddmMat%DDMMat_DILU%nye,method%ddmMat%C_DDMMat_DILU%nye,n,offset_zero,offset_zero)
    endif

    if(allocated(method%ddmMat%DDMMat_DILU%nze))then
       n = size(method%ddmMat%DDMMat_DILU%nze)
       call cudaMalloc_int(method%ddmMat%C_DDMMat_DILU%nze,n)
       call cudaMemcpy_h2d_int(method%ddmMat%DDMMat_DILU%nze,method%ddmMat%C_DDMMat_DILU%nze,n,offset_zero,offset_zero)
    endif

    if(allocated(method%ddmMat%DDMMat_DILU%Dinv))then
       n = size(method%ddmMat%DDMMat_DILU%Dinv)
       call cudaMalloc_fp64(method%ddmMat%C_DDMMat_DILU%Dinv,n)
       call cudaMemcpy_h2d_fp64(method%ddmMat%DDMMat_DILU%Dinv,method%ddmMat%C_DDMMat_DILU%Dinv,n,offset_zero,offset_zero)
    endif

    if(allocated(method%ddmMat%DDMMat_DILU%svec))then
       n = size(method%ddmMat%DDMMat_DILU%svec)
       call cudaMalloc_fp64(method%ddmMat%C_DDMMat_DILU%svec,n)
!       call cudaMemcpy_h2d_fp64(method%ddmMat%DDMMat_DILU%svec,method%ddmMat%C_DDMMat_DILU%svec,n,offset_zero,offset_zero)
    endif
    if(allocated(method%ddmMat%DDMMat_DILU%rvec))then
       n = size(method%ddmMat%DDMMat_DILU%rvec)
       call cudaMalloc_fp64(method%ddmMat%C_DDMMat_DILU%rvec,n)
!       call cudaMemcpy_h2d_fp64(method%ddmMat%DDMMat_DILU%rvec,method%ddmMat%C_DDMMat_DILU%rvec,n,offset_zero,offset_zero)
    endif
    if(allocated(method%ddmMat%DDMMat_DILU%offset_xyz))then
       n = size(method%ddmMat%DDMMat_DILU%offset_xyz)
       call cudaMalloc_int(method%ddmMat%C_DDMMat_DILU%offset_xyz,n)
       call cudaMemcpy_h2d_int(method%ddmMat%DDMMat_DILU%offset_xyz,method%ddmMat%C_DDMMat_DILU%offset_xyz,n,offset_zero,offset_zero)
    endif

  end subroutine init_DILU_DDM_GPU


  subroutine init_DILU_AS_DDM_GPU(method)
    use krylov
    implicit none

    type(KSP) :: method
    integer(kind=intType) n
    integer(kind=intType),parameter :: offset_zero=0

    method%ddmMat%C_DDMMat_DILU%lattice_dimention = method%ddmMat%DDMMat_DILU%lattice_dimention
    method%ddmMat%C_DDMMat_DILU%nx_local = method%ddmMat%DDMMat_DILU%nx_local
    method%ddmMat%C_DDMMat_DILU%ny_local = method%ddmMat%DDMMat_DILU%ny_local
    method%ddmMat%C_DDMMat_DILU%nz_local = method%ddmMat%DDMMat_DILU%nz_local
    method%ddmMat%C_DDMMat_DILU%mx_local = method%ddmMat%DDMMat_DILU%mx_local
    method%ddmMat%C_DDMMat_DILU%my_local = method%ddmMat%DDMMat_DILU%my_local
    method%ddmMat%C_DDMMat_DILU%mz_local = method%ddmMat%DDMMat_DILU%mz_local
    method%ddmMat%C_DDMMat_DILU%m_local  = method%ddmMat%DDMMat_DILU%m_local
    method%ddmMat%C_DDMMat_DILU%nblock   = method%ddmMat%DDMMat_DILU%nblock
    method%ddmMat%C_DDMMat_DILU%nxblock  = method%ddmMat%DDMMat_DILU%nxblock
    method%ddmMat%C_DDMMat_DILU%nyblock  = method%ddmMat%DDMMat_DILU%nyblock
    method%ddmMat%C_DDMMat_DILU%nzblock  = method%ddmMat%DDMMat_DILU%nzblock

    if(allocated(method%ddmMat%DDMMat_DILU%nxs))then
       n = size(method%ddmMat%DDMMat_DILU%nxs)
       call cudaMalloc_int(method%ddmMat%C_DDMMat_DILU%nxs,n)
       call cudaMemcpy_h2d_int(method%ddmMat%DDMMat_DILU%nxs,method%ddmMat%C_DDMMat_DILU%nxs,n,offset_zero,offset_zero)
    endif
    if(allocated(method%ddmMat%DDMMat_DILU%nxe))then
       n = size(method%ddmMat%DDMMat_DILU%nxe)
       call cudaMalloc_int(method%ddmMat%C_DDMMat_DILU%nxe,n)
       call cudaMemcpy_h2d_int(method%ddmMat%DDMMat_DILU%nxe,method%ddmMat%C_DDMMat_DILU%nxe,n,offset_zero,offset_zero)
    endif

    if(allocated(method%ddmMat%DDMMat_DILU%nys))then
       n = size(method%ddmMat%DDMMat_DILU%nys)
       call cudaMalloc_int(method%ddmMat%C_DDMMat_DILU%nys,n)
       call cudaMemcpy_h2d_int(method%ddmMat%DDMMat_DILU%nys,method%ddmMat%C_DDMMat_DILU%nys,n,offset_zero,offset_zero)
    endif
    if(allocated(method%ddmMat%DDMMat_DILU%nye))then
       n = size(method%ddmMat%DDMMat_DILU%nye)
       call cudaMalloc_int(method%ddmMat%C_DDMMat_DILU%nye,n)
       call cudaMemcpy_h2d_int(method%ddmMat%DDMMat_DILU%nye,method%ddmMat%C_DDMMat_DILU%nye,n,offset_zero,offset_zero)
    endif

    if(allocated(method%ddmMat%DDMMat_DILU%nzs))then
       n = size(method%ddmMat%DDMMat_DILU%nzs)
       call cudaMalloc_int(method%ddmMat%C_DDMMat_DILU%nzs,n)
       call cudaMemcpy_h2d_int(method%ddmMat%DDMMat_DILU%nzs,method%ddmMat%C_DDMMat_DILU%nzs,n,offset_zero,offset_zero)
    endif
    if(allocated(method%ddmMat%DDMMat_DILU%nze))then
       n = size(method%ddmMat%DDMMat_DILU%nze)
       call cudaMalloc_int(method%ddmMat%C_DDMMat_DILU%nze,n)
       call cudaMemcpy_h2d_int(method%ddmMat%DDMMat_DILU%nze,method%ddmMat%C_DDMMat_DILU%nze,n,offset_zero,offset_zero)
    endif

    ! --------------------------------------------------

    if(allocated(method%ddmMat%DDMMat_DILU%nxs0))then
       n = size(method%ddmMat%DDMMat_DILU%nxs0)
       call cudaMalloc_int(method%ddmMat%C_DDMMat_DILU%nxs0,n)
       call cudaMemcpy_h2d_int(method%ddmMat%DDMMat_DILU%nxs0,method%ddmMat%C_DDMMat_DILU%nxs0,n,offset_zero,offset_zero)
    endif
    if(allocated(method%ddmMat%DDMMat_DILU%nxe0))then
       n = size(method%ddmMat%DDMMat_DILU%nxe0)
       call cudaMalloc_int(method%ddmMat%C_DDMMat_DILU%nxe0,n)
       call cudaMemcpy_h2d_int(method%ddmMat%DDMMat_DILU%nxe0,method%ddmMat%C_DDMMat_DILU%nxe0,n,offset_zero,offset_zero)
    endif

    if(allocated(method%ddmMat%DDMMat_DILU%nys0))then
       n = size(method%ddmMat%DDMMat_DILU%nys0)
       call cudaMalloc_int(method%ddmMat%C_DDMMat_DILU%nys0,n)
       call cudaMemcpy_h2d_int(method%ddmMat%DDMMat_DILU%nys0,method%ddmMat%C_DDMMat_DILU%nys0,n,offset_zero,offset_zero)
    endif
    if(allocated(method%ddmMat%DDMMat_DILU%nye0))then
       n = size(method%ddmMat%DDMMat_DILU%nye0)
       call cudaMalloc_int(method%ddmMat%C_DDMMat_DILU%nye0,n)
       call cudaMemcpy_h2d_int(method%ddmMat%DDMMat_DILU%nye0,method%ddmMat%C_DDMMat_DILU%nye0,n,offset_zero,offset_zero)
    endif

    if(allocated(method%ddmMat%DDMMat_DILU%nzs0))then
       n = size(method%ddmMat%DDMMat_DILU%nzs0)
       call cudaMalloc_int(method%ddmMat%C_DDMMat_DILU%nzs0,n)
       call cudaMemcpy_h2d_int(method%ddmMat%DDMMat_DILU%nzs0,method%ddmMat%C_DDMMat_DILU%nzs0,n,offset_zero,offset_zero)
    endif
    if(allocated(method%ddmMat%DDMMat_DILU%nze0))then
       n = size(method%ddmMat%DDMMat_DILU%nze0)
       call cudaMalloc_int(method%ddmMat%C_DDMMat_DILU%nze0,n)
       call cudaMemcpy_h2d_int(method%ddmMat%DDMMat_DILU%nze0,method%ddmMat%C_DDMMat_DILU%nze0,n,offset_zero,offset_zero)
    endif

    ! --------------------------------------------------

    if(allocated(method%ddmMat%DDMMat_DILU%Dinv))then
       n = size(method%ddmMat%DDMMat_DILU%Dinv)
       call cudaMalloc_fp64(method%ddmMat%C_DDMMat_DILU%Dinv,n)
       call cudaMemcpy_h2d_fp64(method%ddmMat%DDMMat_DILU%Dinv,method%ddmMat%C_DDMMat_DILU%Dinv,n,offset_zero,offset_zero)
    endif

    if(allocated(method%ddmMat%DDMMat_DILU%svec))then
       n = size(method%ddmMat%DDMMat_DILU%svec)
!       write(*,*)"svec size = ",n
       call cudaMalloc_fp64(method%ddmMat%C_DDMMat_DILU%svec,n)
       call cudaMemcpy_h2d_fp64(method%ddmMat%DDMMat_DILU%svec,method%ddmMat%C_DDMMat_DILU%svec,n,offset_zero,offset_zero)
    endif
    if(allocated(method%ddmMat%DDMMat_DILU%rvec))then
       n = size(method%ddmMat%DDMMat_DILU%rvec)
       call cudaMalloc_fp64(method%ddmMat%C_DDMMat_DILU%rvec,n)
       call cudaMemcpy_h2d_fp64(method%ddmMat%DDMMat_DILU%rvec,method%ddmMat%C_DDMMat_DILU%rvec,n,offset_zero,offset_zero)
    endif
    if(allocated(method%ddmMat%DDMMat_DILU%svec_AS))then
       n = size(method%ddmMat%DDMMat_DILU%svec_AS)
       call cudaMalloc_fp64(method%ddmMat%C_DDMMat_DILU%svec_AS,n)
       call cudaMemcpy_h2d_fp64(method%ddmMat%DDMMat_DILU%svec_AS,method%ddmMat%C_DDMMat_DILU%svec_AS,n,offset_zero,offset_zero)
    endif
    if(allocated(method%ddmMat%DDMMat_DILU%offset_xyz))then
       n = size(method%ddmMat%DDMMat_DILU%offset_xyz)
       call cudaMalloc_int(method%ddmMat%C_DDMMat_DILU%offset_xyz,n)
       call cudaMemcpy_h2d_int(method%ddmMat%DDMMat_DILU%offset_xyz,method%ddmMat%C_DDMMat_DILU%offset_xyz,n,offset_zero,offset_zero)
    endif

    if(allocated(method%ddmMat%AS_v2))then
       n = size(method%ddmMat%AS_v2)
       call cudaMalloc_fp64(method%ddmMat%C_AS_v2,n)
       call cudaMemcpy_h2d_fp64(method%ddmMat%AS_v2,method%ddmMat%C_AS_v2,n,offset_zero,offset_zero)
    endif

    if(allocated(method%ddmMat%AS_v1))then
       n = size(method%ddmMat%AS_v1)
       call cudaMalloc_fp64(method%ddmMat%C_AS_v1,n)
       call cudaMemcpy_h2d_fp64(method%ddmMat%AS_v1,method%ddmMat%C_AS_v1,n,offset_zero,offset_zero)
    endif

    if(allocated(method%ddmMat%ASval%AS_val_dia_init))then
       n = size(method%ddmMat%ASval%AS_val_dia_init)
       call cudaMalloc_fp64(method%ddmMat%ASval%C_AS_val_dia_init,n)
       call cudaMemcpy_h2d_fp64(method%ddmMat%ASval%AS_val_dia_init,method%ddmMat%ASval%C_AS_val_dia_init,n,offset_zero,offset_zero)
    endif

#if 0
    if(allocated(@))then
       n = size(@)
       call cudaMalloc_fp64(method%ddmMat%C_DDMMat_DILU%rvec,n)
       call cudaMemcpy_h2d_fp64(@,method%ddmMat%C_DDMMat_DILU%rvec,n,offset_zero,offset_zero)
    endif
#endif

  end subroutine init_DILU_AS_DDM_GPU

  subroutine free_DILU_DDM_GPU(method)
    use krylov
    implicit none
    type(KSP) :: method

    if(allocated(method%ddmMat%DDMMat_DILU%nxe))call cudaFree(method%ddmMat%C_DDMMat_DILU%nxe)
    if(allocated(method%ddmMat%DDMMat_DILU%nye))call cudaFree(method%ddmMat%C_DDMMat_DILU%nye)
    if(allocated(method%ddmMat%DDMMat_DILU%nze))call cudaFree(method%ddmMat%C_DDMMat_DILU%nze)

    if(allocated(method%ddmMat%DDMMat_DILU%Dinv))call cudaFree(method%ddmMat%C_DDMMat_DILU%Dinv)
    if(allocated(method%ddmMat%DDMMat_DILU%svec))call cudaFree(method%ddmMat%C_DDMMat_DILU%svec)
    if(allocated(method%ddmMat%DDMMat_DILU%rvec))call cudaFree(method%ddmMat%C_DDMMat_DILU%rvec)

    if(allocated(method%ddmMat%DDMMat_DILU%offset_xyz))call cudaFree(method%ddmMat%C_DDMMat_DILU%offset_xyz)

  end subroutine free_DILU_DDM_GPU

  subroutine free_DILU_AS_DDM_GPU(method)
    use krylov
    implicit none
    type(KSP) :: method

    if(allocated(method%ddmMat%DDMMat_DILU%nxs))call cudaFree(method%ddmMat%C_DDMMat_DILU%nxs)
    if(allocated(method%ddmMat%DDMMat_DILU%nys))call cudaFree(method%ddmMat%C_DDMMat_DILU%nys)
    if(allocated(method%ddmMat%DDMMat_DILU%nzs))call cudaFree(method%ddmMat%C_DDMMat_DILU%nzs)
    if(allocated(method%ddmMat%DDMMat_DILU%nxe))call cudaFree(method%ddmMat%C_DDMMat_DILU%nxe)
    if(allocated(method%ddmMat%DDMMat_DILU%nye))call cudaFree(method%ddmMat%C_DDMMat_DILU%nye)
    if(allocated(method%ddmMat%DDMMat_DILU%nze))call cudaFree(method%ddmMat%C_DDMMat_DILU%nze)

    if(allocated(method%ddmMat%DDMMat_DILU%nxs0))call cudaFree(method%ddmMat%C_DDMMat_DILU%nxs0)
    if(allocated(method%ddmMat%DDMMat_DILU%nys0))call cudaFree(method%ddmMat%C_DDMMat_DILU%nys0)
    if(allocated(method%ddmMat%DDMMat_DILU%nzs0))call cudaFree(method%ddmMat%C_DDMMat_DILU%nzs0)
    if(allocated(method%ddmMat%DDMMat_DILU%nxe0))call cudaFree(method%ddmMat%C_DDMMat_DILU%nxe0)
    if(allocated(method%ddmMat%DDMMat_DILU%nye0))call cudaFree(method%ddmMat%C_DDMMat_DILU%nye0)
    if(allocated(method%ddmMat%DDMMat_DILU%nze0))call cudaFree(method%ddmMat%C_DDMMat_DILU%nze0)

    if(allocated(method%ddmMat%DDMMat_DILU%Dinv))call cudaFree(method%ddmMat%C_DDMMat_DILU%Dinv)
    if(allocated(method%ddmMat%DDMMat_DILU%svec))call cudaFree(method%ddmMat%C_DDMMat_DILU%svec)
    if(allocated(method%ddmMat%DDMMat_DILU%rvec))call cudaFree(method%ddmMat%C_DDMMat_DILU%rvec)
    if(allocated(method%ddmMat%DDMMat_DILU%svec_AS))call cudaFree(method%ddmMat%C_DDMMat_DILU%svec_AS)

    if(allocated(method%ddmMat%DDMMat_DILU%offset_xyz))call cudaFree(method%ddmMat%C_DDMMat_DILU%offset_xyz)

  end subroutine free_DILU_AS_DDM_GPU

end module precon_DILU
