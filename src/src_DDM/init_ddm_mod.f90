#include "ddfunc.h"
module init_ddm_param
  use param_AS_type_ddm
  use ddmMat_wrap

  implicit none

  real*8, allocatable :: AS_val_dia_init(:)
  real*8, allocatable :: diag_init(:)
  real*8, allocatable :: AS_val_dia_init_lo(:)
  real*8, allocatable :: diag_init_lo(:)

  integer AS_num
  integer,allocatable :: AS_ioff_dia_init(:)
  type(AS_geom)   :: AS_g
  type(ddm_param) :: AS_prm
  type(ddm_param) :: prm
  type(precon_AS_wrap) :: precon_ddm
  type(ddmmatrix_wrap) :: mat
  integer PRECISION_A_DDM
  integer PRECISION_b_DDM
  integer PRECISION_x_DDM
  integer PRECISION_KRYLOV_DDM
  integer PRECISION_PRECONDITION_DDM ! 前処理行列の精度 (前処理の精度ではない)
  ! ------
  real*8,allocatable :: EigenValueReal(:)
  real*8,allocatable :: EigenValueImag(:)
  real*8,allocatable :: EigenVector(:)
  real*8,allocatable :: EigenErrorNorm(:)
  integer MaxEigenID
  integer MinEigenID
  
contains
  subroutine get_EigenPair_count_DDM_KSP(ncount)
    implicit none
    integer ncount
    if( allocated(EigenValueReal) )then
       ncount = size( EigenValueReal )
    else
       ncount = -1
    endif

  end subroutine get_EigenPair_count_DDM_KSP

  subroutine get_EigenPair_DDM_KSP( &
       ncount,&
       EigenValueReal_DDM,EigenValueImag_DDM, &
       EigenVector_DDM, &
       EigenErrorNorm_DDM, &
       EmaxID_DDM,EminID_DDM)
     use PARCEL_integer
     use vec_kernel,only : xcpy
    implicit none
    integer ncount
    real*8 EigenVector_DDM(:) 
    real*8 EigenValueReal_DDM(ncount) 
    real*8 EigenValueImag_DDM(ncount) 
    real*8 EigenErrorNorm_DDM(ncount) 
    integer EmaxID_DDM,EminID_DDM
    integer(kind=intType) isize
    integer i

    do i = 1,ncount
       EigenValueReal_DDM(i) = EigenValueReal(i)
       EigenValueImag_DDM(i) = EigenValueImag(i)
       EigenErrorNorm_DDM(i) = EigenErrorNorm(i)
    enddo
    isize = size(EigenVector)
    call xcpy(isize,EigenVector,EigenVector_DDM)

    EmaxID_DDM = MaxEigenID
    EminID_DDM = MinEigenID

  end subroutine get_EigenPair_DDM_KSP

  subroutine free_AS_g()
    implicit none
    integer i
    do i = 1,6
       if(allocated(AS_g%send_from_id(i)%id))then
          deallocate(AS_g%send_from_id(i)%id)
       endif
       if(allocated(AS_g%recv_from_id(i)%id))then
          deallocate(AS_g%recv_from_id(i)%id)
       endif
       
       if(allocated(AS_g%send_to_id(i)%id))then
          deallocate(AS_g%send_to_id(i)%id)
       endif
       if(allocated(AS_g%recv_to_id(i)%id))then
          deallocate(AS_g%recv_to_id(i)%id)
       endif
    enddo
  end subroutine free_AS_g
  
  subroutine set_init_ddm_param(&
       n,nnd,&
       input_diag, &
       input_AS_val_dia, &
       input_AS_ioff_dia, &
       PRECISION_KRYLOV, &
       PRECISION_PRECONDITION, &
       input_diag_lo, &
       input_AS_val_dia_lo &
       )
    implicit none
    integer n
    integer nnd
    real*8  input_diag(n)
    real*8  input_AS_val_dia(n*nnd)
    integer input_AS_ioff_dia(nnd)
    integer i

    integer,optional :: PRECISION_PRECONDITION
    integer,optional :: PRECISION_KRYLOV
    real*8,optional ::  input_diag_lo(n)
    real*8,optional ::  input_AS_val_dia_lo(n*nnd)
    integer PRECISION
    real*8 dddiv_var1

    PRECISION = KSP_DOUBLE_PRECISION
    if(  present(PRECISION_KRYLOV)  )then
       PRECISION = PRECISION_KRYLOV
    endif

    allocate(diag_init(n))
    allocate(AS_val_dia_init(n*nnd) )
    allocate(AS_ioff_dia_init(nnd))
    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
       diag_init(i)=input_diag(i)
    enddo
    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n*nnd
       AS_val_dia_init(i)=input_AS_val_dia(i)
    enddo

    if(PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       allocate(diag_init_lo(n))
       allocate(AS_val_dia_init_lo(n*nnd) )
       !$OMP PARALLEL DO PRIVATE(i) 
       do i=1,n
          diag_init_lo(i)=input_diag_lo(i)
       enddo
       !$OMP PARALLEL DO PRIVATE(i) 
       do i=1,n*nnd
          AS_val_dia_init_lo(i)=input_AS_val_dia_lo(i)
       enddo
    endif
    AS_ioff_dia_init=input_AS_ioff_dia

  end subroutine set_init_ddm_param

  subroutine free_init_ddm_param
    implicit none
    ! deallocate(AS_val_dia_init)
    ! deallocate(diag_init)
    ! deallocate(AS_ioff_dia_init)
    if(allocated(diag_init)         ) deallocate( diag_init )
    if(allocated(diag_init_lo)      ) deallocate( diag_init_lo )
    if(allocated(AS_val_dia_init)   ) deallocate( AS_val_dia_init )
    if(allocated(AS_val_dia_init_lo)) deallocate( AS_val_dia_init_lo )
    if(allocated(AS_ioff_dia_init)  ) deallocate( AS_ioff_dia_init )

!    if(allocated(@)  ) deallocate( @ )
    if(allocated(EigenValueReal)  ) deallocate( EigenValueReal )
    if(allocated(EigenValueImag)  ) deallocate( EigenValueImag )
    if(allocated(EigenVector)  ) deallocate( EigenVector )
    if(allocated(EigenErrorNorm)  ) deallocate( EigenErrorNorm )

  end subroutine free_init_ddm_param

  subroutine init_pjacob_ddm( ndim, nnd, mdist, val_dia, ioff_dia, diag, prm ,&
       PRECISION_KRYLOV, &
       PRECISION_A, &
       PRECISION_PRECONDITION, &
       val_dia_lo,diag_lo)
    use ddmMat_wrap
    use mpi
    implicit none
    ! inout
    integer ndim, nnd, mdist, ioff_dia(nnd)
    real*8 val_dia(nnd*mdist), diag(mdist)
    type(ddm_param) :: prm
    integer,optional :: PRECISION_KRYLOV
    integer,optional :: PRECISION_A
    integer,optional :: PRECISION_PRECONDITION
    real*8,optional  :: val_dia_lo(nnd*mdist)
    real*8,optional  :: diag_lo(mdist)

    ! local
    integer nx, ny, nz, mx, my, mz, mxy
    integer stmx, stmy, stmz
    integer jx, jy, jz, j, k, i, idiag, nxy, nxyz
    integer PRECISION
    real*8 dddiv_var1
    real*8 c_hi,c_lo

    PRECISION = KSP_DOUBLE_PRECISION
    if(  present(PRECISION_KRYLOV)  )then
       PRECISION = PRECISION_KRYLOV
    endif

    ! begin
    diag    = 0.0d0
!    if(  present(diag_lo)  )then
    if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       diag_lo = 0.0d0
    endif

    nx=prm%nx;  ny=prm%ny; nz=prm%nz;
    mx=prm%mx;  my=prm%my; mxy=prm%mxy;
    stmx=prm%stmx;  stmy=prm%stmy; stmz=prm%stmz;

    Do i=1, nnd
       if ( ioff_dia(i) == 0 ) idiag = i
    End do

    if(PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       !$OMP PARALLEL do private(jx,jy,jz,j,c_hi,c_lo, &
       !$OMP dddiv_var1) &
       !$OMP collapse(2)
       Do jz=0, nz -1
          Do jy=0, ny -1
             Do jx=0, nx -1
                j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1;
                ! diag(j) = 1.0d0 /val_dia(mdist*(idiag-1) +j)
                if(PRECISION_A.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                   dddiv_calc1(1.0d0,0.0d0,val_dia(mdist*(idiag-1)+j),val_dia_lo(mdist*(idiag-1)+j),c_hi,c_lo)                   
                else
                   dddiv_calc1(1.0d0,0.0d0,val_dia(mdist*(idiag-1)+j),0.0d0,c_hi,c_lo)
                endif
                diag(j) = c_hi
                if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                   diag_lo(j) = c_lo
                endif

             End Do
          End do
       ENd Do
       !$OMP end PARALLEL do
    else
       !$OMP PARALLEL do private(jx,jy,jz,j) collapse(2)
       Do jz=0, nz -1
          Do jy=0, ny -1
             Do jx=0, nx -1
                j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1;
                diag(j) = 1.0d0 /val_dia(mdist*(idiag-1) +j)
             End Do
          End do
       ENd Do
       !$OMP end PARALLEL do
    endif


  end subroutine init_pjacob_ddm

  subroutine as_ddm_init( mdist, nnd, ioff_dia, val_dia, Dinv, &
       prm, ncomp_grids, ioff_grids , precon_thblock, &
       PRECISION_CALC, &
       PRECISION_A, &
       PRECISION_PRECONDITION,&
       val_dia_lo,Dinv_lo)
    use ddmMat_wrap
    use mpi
    use omp_lib
    use common_kernel,only : set_region_for_threads,get_maxthreads,set_precon_maxthreads
    use PARCEL_define,only : KSP_PRECONDITION_DDM_ILUindepend
    use PARCEL_integer
    implicit none
    !   include 'mpif.h'
    !   include "omp_lib.h"
    !   include "ddm.h"

    ! inout
    integer, intent(in) :: nnd, mdist
    integer, intent(in) :: ncomp_grids
    integer precon_thblock
    integer, intent(in) :: ioff_dia(nnd), ioff_grids(nnd,ncomp_grids)
    type(ddm_param), intent(in) :: prm

    real*8, intent(in) :: val_dia(mdist*nnd)
    real*8, intent(inout) :: Dinv(mdist)
    integer,optional :: PRECISION_CALC
    integer,optional :: PRECISION_A
    integer,optional :: PRECISION_PRECONDITION
    real*8,optional  :: val_dia_lo(nnd*mdist)
    real*8,optional  :: Dinv_lo(mdist)

    ! local
    integer k1, j, ip,ipp, idiag, offset, istart
    integer(kind=intType) n
    integer nx, ny, nz, mx, my, mz, mxy, jx, jy, jz, n1, i
    integer stmx, stmy, stmz
    integer nxy, nxyz
    integer num_threads, my_thread
    integer max_numth_0
    integer nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th
    integer  div_dir_th
    integer, allocatable :: ioff_x(:), ioff_y(:), ioff_z(:)
    real*8 c1
    real*8 ddmul_var1
    real*8 dddiv_var1
    real*8 ddmul_var2
    real*8 ddsub_var1
    real*8 c_hi,c_lo

    integer PRECISION
    PRECISION = KSP_DOUBLE_PRECISION
    if(  present(PRECISION_CALC)  )then
       PRECISION = PRECISION_CALC
    endif

    !  write(*,*)" in prm%iflag_DDM_independ = ",prm%iflag_DDM_independ 

    ! init
    nx=prm%nx;  ny=prm%ny; nz=prm%nz;
    mx=prm%mx;  my=prm%my; mxy=prm%mxy;
    stmx=prm%stmx;  stmy=prm%stmy; stmz=prm%stmz;
    div_dir_th = prm%div_direc_th;
    ! ---
    allocate( ioff_x(nnd) );   allocate( ioff_y(nnd) );   allocate( ioff_z(nnd) );
    ioff_x = 0;  ioff_y = 0;  ioff_z = 0

    if ( ncomp_grids == 1 ) then
       Do i=1, nnd
          ioff_x(i) = ioff_grids(i,1)
       End Do
    endif
    if ( ncomp_grids == 2 ) then
       Do i=1, nnd
          ioff_x(i) = ioff_grids(i,1)
          ioff_y(i) = ioff_grids(i,2)
       End Do
    endif
    if ( ncomp_grids == 3 ) then
       Do i=1, nnd
          ioff_x(i) = ioff_grids(i,1)
          ioff_y(i) = ioff_grids(i,2)
          ioff_z(i) = ioff_grids(i,3)
       End Do
    endif

    !$OMP PARALLEL do
    Do i=1, mdist
       Dinv(i) = 1.0d0
    End do
    !$OMP end PARALLEL do
!    if(  present(Dinv_lo)  )then

    if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
       if( PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
          !$OMP PARALLEL do
          Do i=1, mdist
             Dinv_lo(i) = 0.0d0
          End do
          !$OMP end PARALLEL do
       endif
    endif

    ! begin
    Do ip=1, nnd
       if ( ioff_dia(ip) == 0 ) idiag = ip
    End do

    !$OMP PARALLEL do private(jx,jy,jz,j) collapse(2)
    Do jz=0, nz -1
       Do jy=0, ny -1
          Do jx=0, nx -1
             j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1;
             Dinv(j) = val_dia( mdist*(idiag-1) +j )
          End Do
       End do
    ENd Do
    !$OMP end PARALLEL do

    if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
       if( PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       !$OMP PARALLEL do private(jx,jy,jz,j) collapse(2)
          Do jz=0, nz -1
             Do jy=0, ny -1
                Do jx=0, nx -1
                   j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1
                   if( PRECISION_A.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                      Dinv_lo(j) = val_dia_lo( mdist*(idiag-1) +j )
                   else
                      Dinv_lo(j) = 0.0d0
                   endif
                End Do
             End do
          ENd Do
          !$OMP end PARALLEL do       
       endif
    endif
    
    call get_maxthreads(max_numth_0) 
    n = nx * ny * nz 
    call set_precon_maxthreads(n,precon_thblock)

    !$OMP Parallel private( num_threads, my_thread ) &
    !$OMP private( i, jz, jy, jx, j, c1, offset, ip,ipp) &
    !$OMP private( nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th ) &
    !$OMP private( c_hi,c_lo ,&
    !$OMP ddmul_var1, &
    !$OMP dddiv_var1, &
    !$OMP ddmul_var2, &
    !$OMP ddsub_var1 &
    !$OMP ) 

    num_threads = omp_get_num_threads()
    my_thread = omp_get_thread_num()

    call set_region_for_threads( my_thread, num_threads, nx, ny, nz, &
         nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th, &
         div_dir_th )

    Do jz=nzs_th, nze_th -1
       Do jy=nys_th, nye_th -1
          Do jx=nxs_th, nxe_th -1
             j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1;

             if(PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                ! Dinv(j) = 1.0d0 /Dinv(j)
                if( PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                   dddiv_calc1(1.0d0,0.0d0,Dinv(j),Dinv_lo(j),Dinv(j),Dinv_lo(j))
                else
                   dddiv_calc1(1.0d0,0.0d0,Dinv(j),0.0d0,Dinv(j),c_lo)
                endif
             else
                Dinv(j) = 1.0d0 /Dinv(j)
             endif

             Do ip=1, nnd
                if ( ioff_dia(ip) <= 0 ) cycle
                offset = ioff_dia(ip)
                if ( jx +ioff_x(ip) < nxs_th .or. jx +ioff_x(ip) > nxe_th -1 ) cycle
                if ( jy +ioff_y(ip) < nys_th .or. jy +ioff_y(ip) > nye_th -1 ) cycle
                if ( jz +ioff_z(ip) < nzs_th .or. jz +ioff_z(ip) > nze_th -1 ) cycle
                if ( j +offset > mdist ) cycle
                !              Dinv(j+offset) = Dinv(j+offset) - val_dia( mdist*(ip-1) +j )**2 *Dinv(j)

                if(prm%iflag_DDM_independ.eq.KSP_PRECONDITION_DDM_ILUindepend)then
                   !                 write(*,*)"@ Precon @ Bjacobi -- "
                   if( ( ioff_y(ip).eq.0 ).and.( ioff_z(ip).eq.0 ) )then
                      cycle
                   endif
                   do ipp=1,nnd                    
                      if( (offset.eq.( -ioff_dia(ipp) )) )then
                         if(PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                            ! Dinv(j+offset) = Dinv(j+offset) - &
                            !      val_dia( mdist*(ip-1)  +j        ) * &
                            !      val_dia( mdist*(ipp-1) +j+offset ) * &
                            !      Dinv(j)
                            if( PRECISION_A.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                               ddmul_calc1(val_dia(mdist*(ip-1)+j),val_dia_lo(mdist*(ip-1)+j),val_dia(mdist*(ipp-1)+j+offset),val_dia_lo(mdist*(ipp-1)+j+offset),c_hi,c_lo)
                            else
                               ddmul_calc1(val_dia(mdist*(ip-1)+j),0.0d0,val_dia(mdist*(ipp-1)+j+offset),0.0d0,c_hi,c_lo)
                            endif

                            if( PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                               ddmul_calc2(c_hi,c_lo,Dinv(j),Dinv_lo(j),c_hi,c_lo)
                               ddsub_calc1(Dinv(j+offset),Dinv_lo(j+offset),c_hi,c_lo,Dinv(j+offset),Dinv_lo(j+offset))
                            else
                               ddmul_calc2(c_hi,c_lo,Dinv(j),0.0d0,c_hi,c_lo)
                               ddsub_calc1(Dinv(j+offset),0.0d0,c_hi,c_lo,Dinv(j+offset),c_lo)
                            endif
                         else
                            Dinv(j+offset) = Dinv(j+offset) - &
                                 val_dia( mdist*(ip-1)  +j        ) * &
                                 val_dia( mdist*(ipp-1) +j+offset ) * &
                                 Dinv(j)
                         endif
                         exit
                      endif
                   enddo
                else
                   do ipp=1,nnd                    
                      if( (offset.eq.( -ioff_dia(ipp) )) )then
                         if(PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                            ! Dinv(j+offset) = Dinv(j+offset) - &
                            !      val_dia( mdist*(ip-1)  +j        ) * &
                            !      val_dia( mdist*(ipp-1) +j+offset ) * &
                            !      Dinv(j)
                            if( PRECISION_A.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                               ddmul_calc1(val_dia(mdist*(ip-1)+j),val_dia_lo(mdist*(ip-1)+j),val_dia(mdist*(ipp-1)+j+offset),val_dia_lo(mdist*(ipp-1)+j+offset),c_hi,c_lo)
                            else
                               ddmul_calc1(val_dia(mdist*(ip-1)+j),0.0d0,val_dia(mdist*(ipp-1)+j+offset),0.0d0,c_hi,c_lo)
                            endif

                            if( PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                               ddmul_calc2(c_hi,c_lo,Dinv(j),Dinv_lo(j),c_hi,c_lo)
                               ddsub_calc1(Dinv(j+offset),Dinv_lo(j+offset),c_hi,c_lo,Dinv(j+offset),Dinv_lo(j+offset))
                            else
                               ddmul_calc2(c_hi,c_lo,Dinv(j),0.0d0,c_hi,c_lo)
                               ddsub_calc1(Dinv(j+offset),0.0d0,c_hi,c_lo,Dinv(j+offset),c_lo)
                            endif
                         else
                            Dinv(j+offset) = Dinv(j+offset) - &
                                 val_dia( mdist*(ip-1)  +j        ) * &
                                 val_dia( mdist*(ipp-1) +j+offset ) * &
                                 Dinv(j)
                         endif
                         exit
                      endif
                   enddo
                endif

             End do
          End Do
       End Do
    End Do
    !$OMP end PARALLEL

    call omp_set_num_threads(max_numth_0)

    ! Do jz=0, nz -1
    !    Do jy=0, ny -1
    !       Do jx=0, nx -1
    !          j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1;
    !          write(*,*)"Diag",jx,jy,jz,Dinv(j)
    !       End Do
    !    End do
    ! ENd Do


    deallocate( ioff_x );  deallocate( ioff_y );  deallocate( ioff_z )

  end subroutine as_ddm_init


  subroutine parcel_ddm_init(&
       comm,&
       ndim,  &
       m, &
       nnd, &
       ioff_dia, &
       val_dia, &
       num_neighbor_ranks, &
       neighbor_ranks, &
       margin, &
       ncomp_grids, &
       ndiv_grids, num_grids, ioff_grids, &
       div_direc_th, precon_thblock, &
       commhide, nas, iflag_AS,&
       b, x, s, t, oloop, tol, bflag, qrflag, pc, &
       iflag_DDM_independ, &
       nres, res,&
! -----
! optional
       PRECISION_KRYLOV, & 
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION, &
       val_dia_lo,b_lo, x_lo & 
       )

    use ddmMat_wrap
    use ddm_matv
    use param_AS_type_ddm
    use mpi
    use underlap,only : init_underlap,set_diag_underlap
    !    use AS_util_ddm,only : AS_util_set_size,AS_util_set_ioff_A,AS_util_set_A,AS_util_get_buf_id
    use AS_util_ddm,only : AS_util_set_size,AS_util_set_ioff_A,AS_util_set_A
    use PARCEL_define,only : KSP_PRECONDITION_DDM_ILUindepend
    use param_AS_ddm,only : set_AS_DDM,init_AS_Mat
    use PARCEL_integer
    use vec_kernel,only : zero_initialize
    implicit none
    !  include 'iles.h'
    !  include 'mpif.h'

    !arguments
    integer, intent(in) :: comm !MPI communicator
    integer, intent(in) :: ndim !dimension of global matrix
    integer, intent(in) :: m 
    double precision, intent(in) :: b(m) !b vector
    double precision, intent(inout) :: x(m) !x vector
    integer, intent(in) :: s, t
    integer, intent(inout) :: oloop
    double precision, intent(inout) :: tol
    integer, intent(in) :: bflag !basis flag (1:Newton, 2:monominal)
    integer, intent(in) :: qrflag !QR flag  (1:CholeskyQR, 2:TSQR)
    integer, intent(inout) :: pc !pre-conditioner flag
    integer, intent(in) :: iflag_DDM_independ
    integer, intent(in) :: nres !number of residual history
    double precision  res(nres) !residual history

    integer, intent(in) :: nnd 
    integer, intent(in) :: ioff_dia(nnd) !offset
    double precision, intent(in) :: val_dia( m*nnd) !values
    integer, intent(in) :: num_neighbor_ranks
    integer, intent(in) :: neighbor_ranks( num_neighbor_ranks )
    integer, intent(in) :: margin
    integer, intent(in) :: ncomp_grids
    integer, intent(in) :: ndiv_grids(ncomp_grids)
    integer, intent(in) :: num_grids(ncomp_grids)
    integer, intent(in) :: ioff_grids(nnd,ncomp_grids)
    integer, intent(in) :: div_direc_th
    integer, intent(in) :: precon_thblock
    integer underlap_step
    !  type(ddmmatrix_wrap) :: mat
    integer, intent(in) :: commhide ! 0: no hiding, 1: one time, 2: each time
    integer, intent(inout) :: nas
    integer, intent(in) :: iflag_AS  
    integer,optional :: PRECISION_KRYLOV
    integer,optional :: PRECISION_A 
    integer,optional :: PRECISION_b 
    integer,optional :: PRECISION_x 
    integer,optional :: PRECISION_PRECONDITION 
    double precision,optional :: val_dia_lo( m*nnd) ! Matrix values low bit
    double precision,optional :: b_lo(m) ! b vector low bit
    double precision,optional :: x_lo(m) ! x vector low bit

    ! local variables
    integer, dimension(:), allocatable:: iwork
    double precision, dimension(:), allocatable:: dwork
    integer, dimension(:,:), allocatable:: ibk
    integer, dimension(:), allocatable:: nncs, nncr
    integer, dimension(:,:), allocatable::  reqs, reqr
    double precision , dimension(:), allocatable:: diag
    integer iret

    integer leng_halo
    real*8, allocatable :: AS_val_dia(:)
    integer AS_ioff_dia(nnd)

    integer, dimension(:), allocatable:: ptr_l, ptr_u
    double precision, dimension(:), allocatable:: aswork

    integer, dimension(:,:), allocatable:: ibuf
    double precision, dimension(:,:), allocatable:: dbuf

    double precision, dimension(:), allocatable:: wkbufs, wkbufr

    integer, allocatable :: procs_comm_flag(:)

    ! --
    integer  nprocs, myrank, ncpu, newcomm, ierr
    integer  iend, nmax, i, j
    integer  as_n(4), nn_l, nn_u, &
         &     nsend1, nsend2, nrecv1, nrecv2, n_blk, nbuf
    double precision dummy,t_all,t_precon_init,t_cagmres_main
    integer(kind=intType) n_tmp

    real*8, allocatable :: diag_lo(:)
    real*8, allocatable :: AS_val_dia_lo(:)

    !  type(precon_AS_wrap) precon
    underlap_step = 1
    call MPI_Comm_size(comm, nprocs, ierr)
    call MPI_Comm_rank(comm, myrank, ierr)
    ncpu = nprocs

    ! KSP_PRECONDITION_DDM_ILUindepend
    ! prm%iflag_DDM_independ    = -1
    ! AS_prm%iflag_DDM_independ = -1
    if(iflag_DDM_independ>0)then
       prm%iflag_DDM_independ    = KSP_PRECONDITION_DDM_ILUindepend
       AS_prm%iflag_DDM_independ = KSP_PRECONDITION_DDM_ILUindepend
    else
       prm%iflag_DDM_independ    = -1
       AS_prm%iflag_DDM_independ = -1
    endif

    !
    ! input parameters check (INSUFFICIENT !!!!)
    !
    if (pc/=0 .and. pc/=1 .and. pc/=2 .and. pc/=3 .and. pc/=4) then
       write(*,*) '*** preconditioner flag must be 0-4'
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    endif

    ! set block partitioning table
    ! sore table ibk(0:ncpu-1,4) int iwk(libk)
    ! ibk(:,1) = firtst row number of local matrix.
    ! ibk(:,2) = last row number of local matrix.
    ! ibk(:,3) = number of rows (=ibk(:,2)-ibk(:3)+1).
    ! ibk(:,4) = ibk(:,1) - 1
    !
    ! nmax = max(ibk(:,3))
    !

    ! --------------
    ! solverの精度決定
    PRECISION_A_DDM            = KSP_DOUBLE_PRECISION
    PRECISION_b_DDM            = KSP_DOUBLE_PRECISION
    PRECISION_x_DDM            = KSP_DOUBLE_PRECISION
    PRECISION_KRYLOV_DDM       = KSP_DOUBLE_PRECISION
    PRECISION_PRECONDITION_DDM = KSP_DOUBLE_PRECISION

    if(  present(PRECISION_A)  )then
       PRECISION_A_DDM = PRECISION_A
    endif
    if(  present(PRECISION_b)  )then
       PRECISION_b_DDM = PRECISION_b
    endif
    if(  present(PRECISION_x)  )then
       PRECISION_x_DDM = PRECISION_x
    endif
    if(  present(PRECISION_PRECONDITION)  )then
       PRECISION_PRECONDITION_DDM = PRECISION_PRECONDITION
    endif
    if(  present(PRECISION_KRYLOV)  )then
       PRECISION_KRYLOV_DDM = PRECISION_KRYLOV
    endif
    if(PRECISION_KRYLOV_DDM.eq.KSP_DOUBLE_PRECISION)then
       PRECISION_A_DDM            = KSP_DOUBLE_PRECISION
       PRECISION_b_DDM            = KSP_DOUBLE_PRECISION
       PRECISION_x_DDM            = KSP_DOUBLE_PRECISION
       PRECISION_PRECONDITION_DDM = KSP_DOUBLE_PRECISION
    endif
    ! --------------
    ! write(*,*)"PRECISION_A_DDM : ",PRECISION_A_DDM           
    ! write(*,*)"PRECISION_PRECONDITION_DDM : ",PRECISION_PRECONDITION_DDM
    ! write(*,*)"PRECISION_KRYLOV_DDM : ",PRECISION_KRYLOV_DDM
!    write(*,*)"num_neighbor_ranks",num_neighbor_ranks
    leng_halo = maxval( abs(ioff_grids) )
    call init_sleev_comm( ncomp_grids, num_grids, ndiv_grids, &
         &                num_neighbor_ranks, neighbor_ranks, &
         &                leng_halo, margin, prm )
    call set_sleev_comm_flag( nnd, ncomp_grids, ioff_grids, prm )
    call set_div_dir_threads( prm, div_direc_th )
    prm%precon_thblock = precon_thblock
    !
    ! set partitoning table <ptabl> and <rowtbl>, and
    ! make lookup tables for communication before (matrix)x(vector);
    ! <nncs>,<nncr>,<nnxs> and <nnxr>
    !
    allocate(nncs(0:ncpu-1))
    allocate(nncr(0:ncpu-1))
    allocate(reqs(0:ncpu-1,2))
    allocate(reqr(0:ncpu-1,2))

    !allocate(iwork(ndim))
    !  mat%mat_type =  3 !3 means "matrix is DDM type"
    mat%nnd_dia = nnd
    allocate(mat%ioff_dia(size(ioff_dia)))
    allocate(mat%val_dia(size(val_dia)))
    allocate(mat%ioff_grids(nnd,ncomp_grids))
    mat%ioff_dia = ioff_dia
    mat%val_dia = val_dia

    if(  present(val_dia_lo)  )then
       allocate(mat%val_dia_lo(size(val_dia_lo)))
       mat%val_dia_lo = val_dia_lo
    endif

    mat%ioff_grids = ioff_grids
    mat%ndim = ndim
    mat%nrow = m
    mat%prm = prm
    !
    ! prepare priconditioner
    !
    if(pc>=4)then
       write(*,*)"warning :: not support preconditioning : ",pc
       write(*,*)"           change preconditioning to Block Jacobi "
       pc = 2
    endif
    if (pc .eq. 0) then      !! no preconditoner
       AS_num=m
       allocate(diag(m))    ! これは1で良い?
       allocate( AS_val_dia(AS_num*nnd)  ) ! これは1で良い?
       if(PRECISION_KRYLOV_DDM.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
          allocate( diag_lo(m)) ! これは1で良い?
          allocate( AS_val_dia_lo(AS_num*nnd) ) ! これは1で良い?
       endif
    else if (pc .eq. 1) then
       !
       ! point jacobi
       !
       AS_num=m
       allocate(diag(m))
       allocate( AS_val_dia(AS_num*nnd) ) ! これは1で良い?
       if(PRECISION_KRYLOV_DDM.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
          allocate(diag_lo(m)) 
          allocate( AS_val_dia_lo(AS_num*nnd) ) ! これは1で良い?
          call init_pjacob_ddm( ndim, nnd, m, val_dia, ioff_dia, diag, prm , &
               PRECISION_KRYLOV_DDM, &
               PRECISION_A_DDM, &
               PRECISION_PRECONDITION_DDM,&
               val_dia_lo,diag_lo)
       else
          call init_pjacob_ddm( ndim, nnd, m, val_dia, ioff_dia, diag, prm )
       endif
       n_blk = m
    else if (pc .eq. 2  .or.  pc .eq. 3) then
       if   (pc .eq. 2) then
          nas = 0
       endif
       if ( (pc .eq. 3).and.(nas.eq.0) ) then
          nas = 0
          pc=2
       endif

       if( pc == 2 ) then
          AS_num=m
          allocate(diag(m))
          allocate( AS_val_dia(AS_num*nnd) ) 
          if(PRECISION_KRYLOV_DDM.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
             allocate(diag_lo(m)) 
             allocate(AS_val_dia_lo(AS_num*nnd) ) 
             call as_ddm_init( m, nnd, ioff_dia, val_dia, diag, prm, ncomp_grids, &
                  ioff_grids , precon_thblock,&
                  PRECISION_KRYLOV_DDM, &
                  PRECISION_A_DDM, &
                  PRECISION_PRECONDITION_DDM, &
                  AS_val_dia_lo,diag_lo)
          else
             call as_ddm_init( m, nnd, ioff_dia, val_dia, diag, prm, ncomp_grids, &
                  ioff_grids , precon_thblock)
          endif
       else

#if 1
          CALL AS_util_set_size( &
               myrank, m, prm, nas, 1, &
               AS_num, AS_prm, AS_g )
          CALL AS_util_set_ioff_A( &
               nnd, ncomp_grids, ioff_grids, AS_prm, AS_ioff_dia )

          allocate( AS_val_dia(AS_num*nnd) )
          CALL AS_util_set_A( &
               comm, myrank, nnd, m, prm, val_dia, AS_ioff_dia, &
               AS_num, nas, AS_val_dia )

          if(PRECISION_KRYLOV_DDM.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
             allocate( AS_val_dia_lo(AS_num*nnd) )
             if(PRECISION_PRECONDITION_DDM.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                CALL AS_util_set_A( &
                     comm, myrank, nnd, m, prm, val_dia_lo, AS_ioff_dia, &
                     AS_num, nas, AS_val_dia_lo )
             endif
          endif

          allocate(diag(AS_num))
          AS_g%ndist=AS_num
          AS_g%addL=nas
          AS_ioff_dia = mat%ioff_dia
          if(PRECISION_KRYLOV_DDM.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
             allocate(diag_lo(AS_num))
             call set_AS_DDM(&
                  nas, &
                  iflag_AS,&
                  AS_prm,AS_g, &
                  precon_ddm, &
                  AS_num, &
                  nnd, &
                  AS_ioff_dia, & 
                  AS_val_dia, & ! ここに値を入れて後で他の配列にコピーしている
                  ioff_grids, &
                  ncomp_grids, &
                  PRECISION_KRYLOV_DDM, &
                  PRECISION_A_DDM, &
                  PRECISION_PRECONDITION_DDM, &
                  AS_val_dia_lo &
                  )

          else
             call set_AS_DDM(&
                  nas, &
                  iflag_AS,&
                  AS_prm,AS_g, &
                  precon_ddm, &
                  AS_num, &
                  nnd, &
                  AS_ioff_dia, & 
                  AS_val_dia, & ! ここに値を入れて後で他の配列にコピーしている
                  ioff_grids, &
                  ncomp_grids &
                  )

          endif

#else
@@@@@@@@
          allocate(diag(1))
          AS_prm%div_direc_th   = prm%div_direc_th
          AS_prm%precon_thblock = prm%precon_thblock
          CALL AS_util_set_size( &
               myrank, m, prm, nas, 1, &
               AS_num, AS_prm, AS_g )
          allocate( AS_val_dia(AS_num*nnd) )
          n_tmp=AS_num*nnd
          call zero_initialize( n_tmp,AS_val_dia)

          ! 加法シュバルツの拡張領域を含んだ行列値を取得する
          ! in  : val_dia@行列値 -- 袖あり
          ! in  : prm@幾何情報
          !       prm%neighbor_ranks(k) ! 斜め方向の情報を含み通信相手がいない場合は0未満 .3軸分割では26,2軸分割では8プロセスが隣接プロセスとなる.
          !       3次元では (-1,-1,-1),( 0,-1,-1),( 1,-1,-1),(-1,0,-1),........,(1,0,1),(-1,1,1),(0,1,1),(1,1,1) となる.
          ! in  : nas@加法シュバルツの重なり幅
          ! out : AS_val_dia@拡張領域の行列値 -- 袖なし
! fujitsu環境でばぐ?
          call init_AS_Mat( &
               m,AS_num, &
               nnd, &
               prm,AS_prm,AS_g,&
               ioff_grids,&
               comm,val_dia,AS_val_dia,nas)
          ! if(prm%div_direc_th.eq.0)then
          !    call as_ddm_init( &
          !         AS_num, &
          !         nnd, &
          !         AS_ioff_dia, &
          !         AS_val_dia, &
          !         diag, &
          !         AS_prm, &
          !         ncomp_grids, &
          !         ioff_grids, precon_thblock )
          ! else
             AS_g%ndist=AS_num
             AS_g%addL=nas
             AS_ioff_dia = mat%ioff_dia
             call set_AS_DDM(&
                  nas, &
                  iflag_AS,&
                  AS_prm,AS_g, &
                  precon_ddm, &
                  AS_num, &
                  nnd, &
                  AS_ioff_dia, & 
                  AS_val_dia, & ! ここに値を入れて後で他の配列にコピーしている
                  ioff_grids, &
                  ncomp_grids &
                  )
!          endif
#endif

       end if

    else if (pc .eq. 4) then
       write(*,*)"not support "
       ! allocate(diag(m))
       ! AS_num=m
       ! allocate( AS_val_dia(AS_num*nnd) )

       ! if( underlap_step > 0 )then
       !    prm%underlap_step = underlap_step
       ! else
       !    prm%underlap_step = 1
       ! endif
       ! call init_underlap(mat, prm%underlap_step)
       ! call set_diag_underlap(mat, s-1, diag, prm%iflag_DDM_independ )

    end if

    ! if(PRECISION_PRECONDITION_DDM.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
    if(PRECISION_KRYLOV_DDM.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       call set_init_ddm_param(&
            AS_num,nnd,&
            diag, &
            AS_val_dia, &
            AS_ioff_dia, &
            PRECISION_KRYLOV_DDM, &
            PRECISION_PRECONDITION_DDM, &
            diag_lo, &
            AS_val_dia_lo &
            )
    else
       call set_init_ddm_param(&
            AS_num,nnd,&
            diag, &
            AS_val_dia, &
            AS_ioff_dia &
            )
    endif

    deallocate(nncs)
    deallocate(nncr)
    deallocate(reqs)
    deallocate(reqr)

#if 0
    deallocate(ibk)
    deallocate(iwork)
#endif

    deallocate( AS_val_dia )
    deallocate(diag)
    !  if (allocated(@)) deallocate( @ )
    if (allocated(diag_lo)) deallocate( diag_lo )
    if (allocated(AS_val_dia_lo)) deallocate( AS_val_dia_lo )

  end subroutine parcel_ddm_init

  subroutine parcel_ddm_init_dd(comm,&
       ndim, &
       m,&
       nnd, &
       ioff_dia, &
       val_dia, &
       num_neighbor_ranks, &
       neighbor_ranks, &
       margin, &
       ncomp_grids, &
       ndiv_grids, num_grids, ioff_grids, &
       div_direc_th, precon_thblock, &
       commhide, nas, iflag_AS,&
       b, x, &
!      s, t,  tol, bflag, qrflag, pc, &
       oloop,tol,pc, &
       iflag_DDM_independ, &
       nres, res,&
       PRECISION_KRYLOV, & 
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION, &
       val_dia_lo,b_lo, x_lo & 
       )
    integer, intent(in) :: comm !MPI communicator
    integer, intent(in) :: ndim !dimension of global matrix
    integer, intent(in) :: m 
    double precision, intent(in) :: b(m) !b vector
    double precision, intent(inout) :: x(m) !x vector
!    integer, intent(in) :: s, t
    integer, intent(inout) :: oloop
    double precision, intent(inout) :: tol
!    integer, intent(in) :: bflag !basis flag (1:Newton, 2:monominal)
!    integer, intent(in) :: qrflag !QR flag  (1:CholeskyQR, 2:TSQR)
    integer, intent(inout) :: pc !pre-conditioner flag
    integer, intent(in) :: iflag_DDM_independ
    integer, intent(in) :: nres !number of residual history
    double precision  res(nres) !residual history

    integer, intent(in) :: nnd 
    integer, intent(in) :: ioff_dia(nnd) !offset
    double precision, intent(in) :: val_dia( m*nnd) !values
    integer, intent(in) :: num_neighbor_ranks
    integer, intent(in) :: neighbor_ranks( num_neighbor_ranks )
    integer, intent(in) :: margin
    integer, intent(in) :: ncomp_grids
    integer, intent(in) :: ndiv_grids(ncomp_grids)
    integer, intent(in) :: num_grids(ncomp_grids)
    integer, intent(in) :: ioff_grids(nnd,ncomp_grids)
    integer, intent(in) :: div_direc_th
    integer, intent(in) :: precon_thblock
    !  type(ddmmatrix_wrap) :: mat
    integer, intent(in) :: commhide ! 0: no hiding, 1: one time, 2: each time
    integer, intent(inout) :: nas
    integer, intent(in) :: iflag_AS  
    integer :: PRECISION_KRYLOV
    integer :: PRECISION_A 
    integer :: PRECISION_b 
    integer :: PRECISION_x 
    integer :: PRECISION_PRECONDITION 
    double precision :: val_dia_lo( m*nnd) ! Matrix values low bit
    double precision :: b_lo(m) ! b vector low bit
    double precision :: x_lo(m) ! x vector low bit

    call parcel_ddm_init(comm,&
       ndim, &
       m, &
       nnd, &
       ioff_dia, &
       val_dia, &
       num_neighbor_ranks, &
       neighbor_ranks, &
       margin, &
       ncomp_grids, &
       ndiv_grids, num_grids, ioff_grids, &
       div_direc_th, precon_thblock, &
       commhide, nas, iflag_AS,&
!       b, x, s, t, oloop, tol, bflag, qrflag, pc, &
        b, x, -1, -1, oloop, tol, -1, -1, pc, &
       iflag_DDM_independ, &
       nres, res,&
       PRECISION_KRYLOV, & 
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION, &
       val_dia_lo,b_lo, x_lo & 
       )


  end subroutine parcel_ddm_init_dd

  subroutine parcel_ddm_free(ipc)
    use param_AS_ddm,only : free_AS_region_DDM
    implicit none
    integer ipc
    deallocate( prm%neighbor_ranks )
    deallocate( prm%comm_flag )
    deallocate(mat%ioff_dia)
    deallocate(mat%val_dia)
    deallocate(mat%ioff_grids)
    if (ipc .eq. 3) then 
       call free_AS_region_DDM(precon_ddm)
       call free_AS_g()
    endif
    call free_init_ddm_param()

  end subroutine parcel_ddm_free
end module init_ddm_param
