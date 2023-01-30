#include "ddfunc.h"
module init_ddm_param
  use param_AS_type_ddm
  use ddmMat_wrap
  use DDM_TYPE

  implicit none

#if 0
  ! real*8, allocatable :: AS_val_dia_init(:)
  ! real*8, allocatable :: diag_init(:)
  ! real*8, allocatable :: AS_val_dia_init_lo(:)
  ! real*8, allocatable :: diag_init_lo(:)

  ! integer AS_num
  ! integer,allocatable :: AS_ioff_dia_init(:)
  ! type(AS_geom)   :: AS_g
  ! type(ddm_param) :: AS_prm
  ! type(ddm_param) :: prm
  ! type(precon_AS_wrap) :: precon_ddm
  ! type(ddmmatrix_wrap) :: mat
  ! integer PRECISION_A_DDM
  ! integer PRECISION_b_DDM
  ! integer PRECISION_x_DDM
  ! integer PRECISION_KRYLOV_DDM
  ! integer PRECISION_PRECONDITION_DDM ! 前処理行列の精度 (前処理の精度ではない)
  ! ------
  ! real*8,allocatable :: EigenValueReal(:)
  ! real*8,allocatable :: EigenValueImag(:)
  ! real*8,allocatable :: EigenVector(:)
  ! real*8,allocatable :: EigenErrorNorm(:)
  ! integer MaxEigenID
  ! integer MinEigenID
#endif
  
contains

#if 0

  subroutine get_EigenPair_count_DDM_KSP(ncount)
    implicit none
    integer ncount
#if 0
    if( allocated(EigenValueReal) )then
       ncount = size( EigenValueReal )
    else
       ncount = -1
    endif
#endif
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

#if 0
    do i = 1,ncount
       EigenValueReal_DDM(i) = EigenValueReal(i)
       EigenValueImag_DDM(i) = EigenValueImag(i)
       EigenErrorNorm_DDM(i) = EigenErrorNorm(i)
    enddo
    isize = size(EigenVector)
    call xcpy(isize,EigenVector,EigenVector_DDM)

    EmaxID_DDM = MaxEigenID
    EminID_DDM = MinEigenID
#endif

  end subroutine get_EigenPair_DDM_KSP
#endif

  subroutine free_AS_g(DDM_MAT_)
    implicit none
    type (DDM_ALL_TYPE) DDM_MAT_
    integer i

    do i = 1,6
       if(allocated(DDM_MAT_%AS_g%send_from_id(i)%id))then
          deallocate(DDM_MAT_%AS_g%send_from_id(i)%id)
       endif
       if(allocated(DDM_MAT_%AS_g%recv_from_id(i)%id))then
          deallocate(DDM_MAT_%AS_g%recv_from_id(i)%id)
       endif
       
       if(allocated(DDM_MAT_%AS_g%send_to_id(i)%id))then
          deallocate(DDM_MAT_%AS_g%send_to_id(i)%id)
       endif
       if(allocated(DDM_MAT_%AS_g%recv_to_id(i)%id))then
          deallocate(DDM_MAT_%AS_g%recv_to_id(i)%id)
       endif
    enddo

  end subroutine free_AS_g
  
  subroutine set_init_ddm_param(&
       DDM_MAT_, &
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
    type (DDM_ALL_TYPE) DDM_MAT_
    integer(kind=intType) n
    integer(kind=intType) nnd
    real*8  input_diag(n)
    real*8  input_AS_val_dia(n*nnd)
    integer input_AS_ioff_dia(nnd)
    integer(kind=intType) i

    integer,optional :: PRECISION_PRECONDITION
    integer,optional :: PRECISION_KRYLOV
    real*8,optional ::  input_diag_lo(n)
    real*8,optional ::  input_AS_val_dia_lo(n*nnd)
    integer PRECISION
    real*8 dddiv_var1
    ! integer(kind=intType) parcel_int

    PRECISION = KSP_DOUBLE_PRECISION
    if(  present(PRECISION_KRYLOV)  )then
       PRECISION = PRECISION_KRYLOV
    endif

    if(allocated(DDM_MAT_%ASval%diag_init)) deallocate(DDM_MAT_%ASval%diag_init)
    if(allocated(DDM_MAT_%ASval%AS_val_dia_init)) deallocate(DDM_MAT_%ASval%AS_val_dia_init)
    if(allocated(DDM_MAT_%ASval%AS_ioff_dia_init)) deallocate(DDM_MAT_%ASval%AS_ioff_dia_init)

    allocate(DDM_MAT_%ASval%diag_init(n))
    allocate(DDM_MAT_%ASval%AS_val_dia_init(n*nnd) )
    allocate(DDM_MAT_%ASval%AS_ioff_dia_init(nnd))

    !$OMP PARALLEL DO PRIVATE(i) 
    do i = 1,n
       DDM_MAT_%ASval%diag_init(i) = input_diag(i)
    enddo
    !$OMP PARALLEL DO PRIVATE(i) 

    do i = 1,n*nnd
       DDM_MAT_%ASval%AS_val_dia_init(i) = input_AS_val_dia(i)
    enddo

    if(PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       allocate(DDM_MAT_%ASval%diag_init_lo(n))
       allocate(DDM_MAT_%ASval%AS_val_dia_init_lo(n*nnd) )
       !$OMP PARALLEL DO PRIVATE(i) 
       do i = 1,n
          DDM_MAT_%ASval%diag_init_lo(i) = input_diag_lo(i)
       enddo
       !$OMP PARALLEL DO PRIVATE(i) 
       do i = 1,n*nnd
          DDM_MAT_%ASval%AS_val_dia_init_lo(i) = input_AS_val_dia_lo(i)
       enddo
    endif
    DDM_MAT_%ASval%AS_ioff_dia_init = input_AS_ioff_dia

  end subroutine set_init_ddm_param

  subroutine free_init_ddm_param(DDM_MAT_)
    implicit none
    type (DDM_ALL_TYPE) DDM_MAT_

    ! deallocate(AS_val_dia_init)
    ! deallocate(diag_init)
    ! deallocate(AS_ioff_dia_init)
    if(allocated(DDM_MAT_%ASval%diag_init)         )  &
         deallocate( DDM_MAT_%ASval%diag_init )
    if(allocated(DDM_MAT_%ASval%diag_init_lo)      )  &
         deallocate( DDM_MAT_%ASval%diag_init_lo )
    if(allocated(DDM_MAT_%ASval%AS_val_dia_init)   )  &
         deallocate( DDM_MAT_%ASval%AS_val_dia_init )
    if(allocated(DDM_MAT_%ASval%AS_val_dia_init_lo))  &
         deallocate( DDM_MAT_%ASval%AS_val_dia_init_lo )
    if(allocated(DDM_MAT_%ASval%AS_ioff_dia_init)  )  &
         deallocate( DDM_MAT_%ASval%AS_ioff_dia_init )

#if 0
!    if(allocated(@)  ) deallocate( @ )
    if(allocated(EigenValueReal)  ) deallocate( EigenValueReal )
    if(allocated(EigenValueImag)  ) deallocate( EigenValueImag )
    if(allocated(EigenVector)  ) deallocate( EigenVector )
    if(allocated(EigenErrorNorm)  ) deallocate( EigenErrorNorm )
#endif

  end subroutine free_init_ddm_param

  subroutine init_pjacob_ddm( &
       ndim, nnd, mdist, val_dia, ioff_dia, diag, prm ,&
       PRECISION_KRYLOV, &
       PRECISION_A, &
       PRECISION_PRECONDITION, &
       val_dia_lo,diag_lo)
    use ddmMat_wrap
    use mpi
    implicit none
    ! inout
    integer(kind=intType) mdist
    integer ndim, nnd,  ioff_dia(nnd)
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
    if(  present(diag_lo)  )then
       if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
          diag_lo = 0.0d0
       endif
    endif

    nx=prm%nx;  ny=prm%ny; nz=prm%nz;
    mx=prm%mx;  my=prm%my; mxy=prm%mxy;
    stmx=prm%stmx;  stmy=prm%stmy; stmz=prm%stmz;

    idiag = -1
    Do i=1, nnd
       if ( ioff_dia(i) == 0 ) idiag = i
    End do

    if(idiag.eq.-1)then
       diag = 1.0d0
    endif

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

  subroutine init_neumann_ddm( &
       ndim, nnd, mdist, &
       val_dia, ioff_dia, diag,&
       prm ,&
       val_dia_Neuman, ioff_dia_Neuman)

    use ddmMat_wrap
    use mpi
    implicit none

    ! inout
    integer(kind=intType) mdist
    integer ndim, nnd,  ioff_dia(nnd)
    real*8 val_dia(nnd*mdist), diag(mdist)

    integer ioff_dia_Neuman(nnd-1)
    real*8 val_dia_Neuman((nnd-1)*mdist)
    type(ddm_param) :: prm

    ! local
    integer nx, ny, nz, mx, my, mz, mxy
    integer stmx, stmy, stmz
    integer(kind=intType) jx, jy, jz, j, k, i, idiag, nxy, nxyz
    integer(kind=intType) ii

    integer PRECISION
    real*8 dddiv_var1
    real*8 c_hi,c_lo


    ! begin
    diag    = 0.0d0

    nx=prm%nx;  ny=prm%ny; nz=prm%nz;
    mx=prm%mx;  my=prm%my; mxy=prm%mxy;
    stmx=prm%stmx;  stmy=prm%stmy; stmz=prm%stmz;

    ii = 1
    idiag = -1
    do i=1, nnd
       if ( ioff_dia(i) == 0 )then
          idiag = i
       else
          ioff_dia_Neuman(ii) = ioff_dia(i)
          ii = ii + 1
       endif
    end do

    if(idiag.eq.-1)then
       diag = 1.0d0
    endif

    !$OMP PARALLEL do private(jx,jy,jz,j,ii) collapse(2)
    do jz=0, nz -1
       do jy=0, ny -1
          do jx=0, nx -1
             j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1
!             write(*,*)j,mdist*(idiag-1) +j
             diag(j) = 1.0d0 /val_dia(mdist*(idiag-1) +j)
             ii = 1
             do i=1, nnd
                if ( ioff_dia(i) .ne. 0 )then
                   val_dia_Neuman(j+mdist*(ii-1)) = val_dia(j+mdist*(i-1)) * diag(j)
                   ii = ii + 1
                endif
             enddo
          end do
       end do
    end do
    !$OMP end PARALLEL do
    
  end subroutine init_neumann_ddm


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
    integer(kind=intType), intent(in) :: mdist
    integer, intent(in) :: nnd
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
    integer(kind=intType) k1, j, ip,ipp, idiag, offset, istart,  n1, i
    integer(kind=intType) n
    integer nx, ny, nz, mx, my, mz, mxy, jx, jy, jz
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
#if 0
    Do jz=0, nz -1
       Do jy=0, ny -1
          Do jx=0, nx -1
             j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1;
             write(*,*)"Diag",jx,jy,jz,Dinv(j)
          End Do
       End do
    ENd Do
#endif

    deallocate( ioff_x );  deallocate( ioff_y );  deallocate( ioff_z )

  end subroutine as_ddm_init


  subroutine ILU_subdividing_ddm_init( &
       mdist, nnd, ioff_dia, val_dia, Dinv, &
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
    integer(kind=intType) mdist
    integer, intent(in) :: nnd
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
    integer(kind=intType) k1, j, ip,ipp, idiag, offset, istart, n1, i
    integer(kind=intType) n
    integer nx, ny, nz, mx, my, mz, mxy, jx, jy, jz
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
    
#if 0
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

#endif
    deallocate( ioff_x );  deallocate( ioff_y );  deallocate( ioff_z )

  end subroutine ILU_subdividing_ddm_init

  subroutine parcel_ddm_init(&
       DDM_MAT_,&
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
!    use param_AS_ddm,only : set_AS_DDM,init_AS_Mat
    use param_AS_ddm,only : set_AS_DDM
    use PARCEL_integer
    use vec_kernel,only : zero_initialize
    implicit none

    type (DDM_ALL_TYPE) DDM_MAT_
    !arguments
    integer, intent(in) :: comm !MPI communicator
    integer, intent(in) :: ndim !dimension of global matrix
    integer(kind=intType) , intent(in) :: m 
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

    ! ----------------
    integer  nprocs, myrank, ncpu, newcomm, ierr
    integer  iend, nmax, i, j
    integer  as_n(4), nn_l, nn_u, &
         nsend1, nsend2, nrecv1, nrecv2, n_blk, nbuf
    double precision dummy,t_all,t_precon_init,t_cagmres_main
    integer(kind=intType) n_tmp
    integer(kind=intType) parcel_int
 
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
       DDM_MAT_%prm%iflag_DDM_independ    = KSP_PRECONDITION_DDM_ILUindepend
       DDM_MAT_%AS_prm%iflag_DDM_independ = KSP_PRECONDITION_DDM_ILUindepend
    else
       DDM_MAT_%prm%iflag_DDM_independ    = -1
       DDM_MAT_%AS_prm%iflag_DDM_independ = -1
    endif

    ! input parameters check (INSUFFICIENT !!!!)
#if 0
    if (pc/=0 .and. pc/=1 .and. pc/=2 .and. pc/=3 .and. pc/=4) then
       write(*,*) '*** preconditioner flag must be 0-4'
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    endif
#endif

    ! set block partitioning table
    ! sore table ibk(0:ncpu-1,4) int iwk(libk)
    ! ibk(:,1) = firtst row number of local matrix.
    ! ibk(:,2) = last row number of local matrix.
    ! ibk(:,3) = number of rows (=ibk(:,2)-ibk(:3)+1).
    ! ibk(:,4) = ibk(:,1) - 1
    ! nmax = max(ibk(:,3))

    ! --------------
    ! solverの精度決定
    DDM_MAT_%PRECISION_A_DDM            = KSP_DOUBLE_PRECISION
    DDM_MAT_%PRECISION_b_DDM            = KSP_DOUBLE_PRECISION
    DDM_MAT_%PRECISION_x_DDM            = KSP_DOUBLE_PRECISION
    DDM_MAT_%PRECISION_KRYLOV_DDM       = KSP_DOUBLE_PRECISION
    DDM_MAT_%PRECISION_PRECONDITION_DDM = KSP_DOUBLE_PRECISION

    if(  present(PRECISION_A)  )then
       DDM_MAT_%PRECISION_A_DDM = PRECISION_A
    endif
    if(  present(PRECISION_b)  )then
       DDM_MAT_%PRECISION_b_DDM = PRECISION_b
    endif
    if(  present(PRECISION_x)  )then
       DDM_MAT_%PRECISION_x_DDM = PRECISION_x
    endif
    if(  present(PRECISION_PRECONDITION)  )then
       DDM_MAT_%PRECISION_PRECONDITION_DDM = PRECISION_PRECONDITION
    endif
    if(  present(PRECISION_KRYLOV)  )then
       DDM_MAT_%PRECISION_KRYLOV_DDM = PRECISION_KRYLOV
    endif
    if(DDM_MAT_%PRECISION_KRYLOV_DDM.eq.KSP_DOUBLE_PRECISION)then
       DDM_MAT_%PRECISION_A_DDM            = KSP_DOUBLE_PRECISION
       DDM_MAT_%PRECISION_b_DDM            = KSP_DOUBLE_PRECISION
       DDM_MAT_%PRECISION_x_DDM            = KSP_DOUBLE_PRECISION
       DDM_MAT_%PRECISION_PRECONDITION_DDM = KSP_DOUBLE_PRECISION
    endif

    ! --------------
    ! write(*,*)"PRECISION_A_DDM : ",PRECISION_A_DDM           
    ! write(*,*)"PRECISION_PRECONDITION_DDM : ",PRECISION_PRECONDITION_DDM
    ! write(*,*)"PRECISION_KRYLOV_DDM : ",PRECISION_KRYLOV_DDM
    ! write(*,*)"num_neighbor_ranks",num_neighbor_ranks
    leng_halo = maxval( abs(ioff_grids) )
    call init_sleev_comm( ncomp_grids, num_grids, ndiv_grids, &
                        num_neighbor_ranks, neighbor_ranks, &
                        leng_halo, margin, DDM_MAT_%prm )
    call set_sleev_comm_flag( nnd, ncomp_grids, ioff_grids, DDM_MAT_%prm )
    call set_div_dir_threads( DDM_MAT_%prm, div_direc_th )
    DDM_MAT_%prm%precon_thblock = precon_thblock

#ifdef CUDA_SOLVER
    allocate(DDM_MAT_%mat%ioff_dia(size(ioff_dia)))
    allocate(DDM_MAT_%mat%val_dia(size(val_dia)))
    allocate(DDM_MAT_%mat%ioff_grids(nnd,ncomp_grids))

    DDM_MAT_%mat%ioff_dia = ioff_dia
    DDM_MAT_%mat%val_dia = val_dia
    DDM_MAT_%mat%nnd_dia = nnd

    DDM_MAT_%mat%ioff_grids = ioff_grids
    DDM_MAT_%mat%ndim = ndim
    DDM_MAT_%mat%nrow = m

    DDM_MAT_%ASval%AS_num=m

    allocate(DDM_MAT_%ASval%diag_init(m))
    allocate(DDM_MAT_%ASval%AS_val_dia_init(1) )
    allocate(DDM_MAT_%ASval%AS_ioff_dia_init(1))

    if(pc.ne.22)then
       return 
    endif

#endif

    ! set partitoning table <ptabl> and <rowtbl>, and
    ! make lookup tables for communication before (matrix)x(vector);
    ! <nncs>,<nncr>,<nnxs> and <nnxr>
    allocate(nncs(0:ncpu-1))
    allocate(nncr(0:ncpu-1))
    allocate(reqs(0:ncpu-1,2))
    allocate(reqr(0:ncpu-1,2))

    !allocate(iwork(ndim))
    !  mat%mat_type =  3 !3 means "matrix is DDM type"
    DDM_MAT_%mat%nnd_dia = nnd

    if(allocated(DDM_MAT_%mat%ioff_dia)) deallocate(DDM_MAT_%mat%ioff_dia)
    if(allocated(DDM_MAT_%mat%val_dia)) deallocate(DDM_MAT_%mat%val_dia)
    if(allocated(DDM_MAT_%mat%ioff_grids)) deallocate(DDM_MAT_%mat%ioff_grids)
    allocate(DDM_MAT_%mat%ioff_dia(size(ioff_dia)))
    allocate(DDM_MAT_%mat%val_dia(size(val_dia)))
    allocate(DDM_MAT_%mat%ioff_grids(nnd,ncomp_grids))
    DDM_MAT_%mat%ioff_dia = ioff_dia
    DDM_MAT_%mat%val_dia = val_dia

    if(  present(val_dia_lo)  )then
       allocate(DDM_MAT_%mat%val_dia_lo(size(val_dia_lo)))
       DDM_MAT_%mat%val_dia_lo = val_dia_lo
    endif

    DDM_MAT_%mat%ioff_grids = ioff_grids
    DDM_MAT_%mat%ndim = ndim
    DDM_MAT_%mat%nrow = m

    ! prepare priconditioner
#if 0
    if(pc>=5)then
       write(*,*)"warning :: not support preconditioning : ",pc
       write(*,*)"           change preconditioning to Block Jacobi "
       pc = 2
    endif
#endif

    if (pc .eq. 0) then !! no preconditoner
       DDM_MAT_%ASval%AS_num=m
       allocate(diag(m)) 
       allocate( AS_val_dia(m*nnd)  ) 
       if(DDM_MAT_%PRECISION_KRYLOV_DDM.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
          allocate( diag_lo(m)) 
          allocate( AS_val_dia_lo(m*nnd) ) 
       endif
    else if (pc .eq. 1) then
       ! point jacobi
       DDM_MAT_%ASval%AS_num=m
       allocate( diag(m) )
       allocate( AS_val_dia(m*nnd) )
       if(DDM_MAT_%PRECISION_KRYLOV_DDM.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
          allocate(diag_lo(m)) 
          allocate( AS_val_dia_lo(m*nnd) )
          call init_pjacob_ddm( &
               ndim, nnd, m, val_dia, ioff_dia, diag, DDM_MAT_%prm , &
               DDM_MAT_%PRECISION_KRYLOV_DDM, &
               DDM_MAT_%PRECISION_A_DDM, &
               DDM_MAT_%PRECISION_PRECONDITION_DDM,&
               val_dia_lo,diag_lo)
       else
          call init_pjacob_ddm( &
               ndim, nnd, m, val_dia, ioff_dia, diag, DDM_MAT_%prm )
       endif
       n_blk = m
    else if ( (pc .eq. 2)  .or.  (pc .eq. 3) .or. (pc .eq. 22) ) then
       if   (pc .eq. 2) then
          nas = 0
       endif
       if (  (pc .eq. 3).and.(nas.eq.0) ) then
          nas = 0
          pc=2
       endif

       if( pc == 2 ) then
          DDM_MAT_%ASval%AS_num=m
          allocate(diag(m))
          allocate( AS_val_dia(m*nnd) ) 
          if(DDM_MAT_%PRECISION_KRYLOV_DDM.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
             allocate(diag_lo(m)) 
             allocate(AS_val_dia_lo(m*nnd) ) 
             call as_ddm_init( m, nnd, ioff_dia, val_dia, diag, DDM_MAT_%prm, ncomp_grids, &
                  ioff_grids , precon_thblock,&
                  DDM_MAT_%PRECISION_KRYLOV_DDM, &
                  DDM_MAT_%PRECISION_A_DDM, &
                  DDM_MAT_%PRECISION_PRECONDITION_DDM, &
                  AS_val_dia_lo,diag_lo)
          else
             call as_ddm_init( m, nnd, ioff_dia, val_dia, diag, DDM_MAT_%prm, ncomp_grids, &
                  ioff_grids , precon_thblock)
          endif
       else

          CALL AS_util_set_size( &
               myrank, m, DDM_MAT_%prm, nas, 1, &
               DDM_MAT_%ASval%AS_num, DDM_MAT_%AS_prm, DDM_MAT_%AS_g )
          CALL AS_util_set_ioff_A( &
               nnd, ncomp_grids, ioff_grids, DDM_MAT_%AS_prm, AS_ioff_dia )

          ! DDM_MAT_%ASval%AS_num : 袖なしのサイズ
          ! write(*,*)" AS region size : DDM_MAT_%ASval%AS_num = ",DDM_MAT_%ASval%AS_num 
          allocate( AS_val_dia(DDM_MAT_%ASval%AS_num*nnd) )
          CALL AS_util_set_A( &
               comm, myrank, nnd, m, DDM_MAT_%prm, val_dia, AS_ioff_dia, &
               DDM_MAT_%ASval%AS_num, nas, AS_val_dia )
          ! AS_val_dia : 拡張領域の行列データ
 
          if(DDM_MAT_%PRECISION_KRYLOV_DDM.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
             allocate( AS_val_dia_lo(DDM_MAT_%ASval%AS_num*nnd) )
             if(DDM_MAT_%PRECISION_PRECONDITION_DDM.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                CALL AS_util_set_A( &
                     comm, myrank, nnd, m, DDM_MAT_%prm, val_dia_lo, AS_ioff_dia, &
                     DDM_MAT_%ASval%AS_num, nas, AS_val_dia_lo )
             endif
          endif

          allocate(diag(DDM_MAT_%ASval%AS_num))
          DDM_MAT_%AS_g%ndist=DDM_MAT_%ASval%AS_num
          DDM_MAT_%AS_g%addL=nas
          AS_ioff_dia = DDM_MAT_%mat%ioff_dia
          if(DDM_MAT_%PRECISION_KRYLOV_DDM.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
             allocate(diag_lo(DDM_MAT_%ASval%AS_num))
             call set_AS_DDM(&
                  nas, &
                  iflag_AS,&
                  DDM_MAT_%AS_prm, &
                  DDM_MAT_%AS_g  , &
                  DDM_MAT_%precon_ddm, &
                  DDM_MAT_%ASval%AS_num, &
                  nnd, &
                  AS_ioff_dia, & ! offset
                  AS_val_dia, & ! ここに値を入れて後で他の配列にコピーしている
                  ioff_grids, &
                  ncomp_grids, &
                  DDM_MAT_%PRECISION_KRYLOV_DDM, &
                  DDM_MAT_%PRECISION_A_DDM, &
                  DDM_MAT_%PRECISION_PRECONDITION_DDM, &
                  AS_val_dia_lo &
                  )

          else
             call set_AS_DDM(&
                  nas, &
                  iflag_AS,&
                  DDM_MAT_%AS_prm,DDM_MAT_%AS_g, &
                  DDM_MAT_%precon_ddm, &
                  DDM_MAT_%ASval%AS_num, &
                  nnd, &
                  AS_ioff_dia, & 
                  AS_val_dia, & ! ここに値を入れて後で他の配列にコピーしている
                  ioff_grids, &
                  ncomp_grids &
                  )

          endif


       end if

    else if ( (pc.eq.4) .or. (pc.eq.5) ) then
       ! Neumann precondition
       ! 対角要素の逆数Dinvが必要
       DDM_MAT_%ASval%AS_num=m
       allocate(diag(m))
       allocate( AS_val_dia(m*nnd) )

       allocate( DDM_MAT_%DDMMat_Neumann%ddmA(m*(nnd-1)) )
       allocate( DDM_MAT_%DDMMat_Neumann%offset(nnd-1)   )

      DDM_MAT_%DDMMat_Neumann%nx   = DDM_MAT_%prm%nx
      DDM_MAT_%DDMMat_Neumann%ny   = DDM_MAT_%prm%ny
      DDM_MAT_%DDMMat_Neumann%nz   = DDM_MAT_%prm%nz
      DDM_MAT_%DDMMat_Neumann%stmx = DDM_MAT_%prm%stmx
      DDM_MAT_%DDMMat_Neumann%stmy = DDM_MAT_%prm%stmy
      DDM_MAT_%DDMMat_Neumann%stmz = DDM_MAT_%prm%stmz
       
      DDM_MAT_%DDMMat_Neumann%mx  = DDM_MAT_%prm%mx
      DDM_MAT_%DDMMat_Neumann%my  = DDM_MAT_%prm%my
      DDM_MAT_%DDMMat_Neumann%mz  = DDM_MAT_%prm%mz
      DDM_MAT_%DDMMat_Neumann%mxy = DDM_MAT_%prm%mx * DDM_MAT_%prm%my 
      DDM_MAT_%DDMMat_Neumann%m   = DDM_MAT_%prm%mx * DDM_MAT_%prm%my * DDM_MAT_%prm%mz
      DDM_MAT_%DDMMat_Neumann%nnd = DDM_MAT_%mat%nnd_dia - 1

       call init_neumann_ddm( &
            ndim, nnd, m, &
            val_dia, ioff_dia, diag,&
            DDM_MAT_%prm , &
            DDM_MAT_%DDMMat_Neumann%ddmA, DDM_MAT_%DDMMat_Neumann%offset)

    end if

    if (allocated(diag))then
    else
       allocate(diag(m))
    endif

    if (allocated(AS_val_dia))then
    else
       allocate(AS_val_dia(1))
    endif

    ! if(PRECISION_PRECONDITION_DDM.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
    if(DDM_MAT_%PRECISION_KRYLOV_DDM.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       parcel_int = nnd
       call set_init_ddm_param(&
            DDM_MAT_, &
            DDM_MAT_%ASval%AS_num,parcel_int,&
            diag, &
            AS_val_dia, &
            AS_ioff_dia, &
            DDM_MAT_%PRECISION_KRYLOV_DDM, &
            DDM_MAT_%PRECISION_PRECONDITION_DDM, &
            diag_lo, &
            AS_val_dia_lo &
            )
    else
       parcel_int = nnd
       call set_init_ddm_param(&
            DDM_MAT_, &
            DDM_MAT_%ASval%AS_num,  parcel_int,&
            diag, &
            AS_val_dia, &
            AS_ioff_dia &
            )
    endif
 
    if(allocated(nncs))deallocate(nncs)
    if(allocated(nncr))deallocate(nncr)
    if(allocated(reqs))deallocate(reqs)
    if(allocated(reqr))deallocate(reqr)

    if(allocated(AS_val_dia))deallocate( AS_val_dia )
    if(allocated(diag))deallocate(diag)

    if (allocated(diag_lo)) deallocate( diag_lo )
    if (allocated(AS_val_dia_lo)) deallocate( AS_val_dia_lo )

  end subroutine parcel_ddm_init

  subroutine parcel_ddm_init_dd(&
       DDM_MAT_,&
       comm,&
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
    implicit none
    type (DDM_ALL_TYPE) DDM_MAT_
    integer, intent(in) :: comm !MPI communicator
    integer, intent(in) :: ndim !dimension of global matrix
    integer(kind=intType) , intent(in) :: m 
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

    call parcel_ddm_init(&
       DDM_MAT_,&
       comm, &
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

  subroutine parcel_ddm_free(ipc,DDM_MAT_)
    use param_AS_ddm,only : free_AS_region_DDM
    implicit none
    type (DDM_ALL_TYPE) DDM_MAT_
    integer ipc

    deallocate(DDM_MAT_%prm%neighbor_ranks )
    deallocate(DDM_MAT_%prm%comm_flag )
    deallocate(DDM_MAT_%mat%ioff_dia)
    deallocate(DDM_MAT_%mat%val_dia)
    deallocate(DDM_MAT_%mat%ioff_grids)
    if (ipc .eq. 3) then 
       call free_AS_region_DDM(DDM_MAT_%precon_ddm)
       call free_AS_g(DDM_MAT_)
    endif
    call free_init_ddm_param(DDM_MAT_)

  end subroutine parcel_ddm_free
end module init_ddm_param
