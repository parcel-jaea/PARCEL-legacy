#include "ddfunc.h"

module as_ddm
  use PARCEL_define,only : KSP_DOUBLE_DOUBLE_PRECISION,KSP_DOUBLE_PRECISION

contains

  subroutine aslu_ddm( mdist, nnd, val_dia, ioff_dia, &
       diag, vec1, vec2, prm, ncomp_grids, num_grids, ioff_grids,precon,&
       PRECISION_CALC, &
       PRECISION_PRECONDITION , &
       val_dia_lo, &
       vec1_lo, &
       vec2_lo, &
       diag_lo &
       )
    use param_AS_type_ddm,only : precon_AS_wrap, &
         KSP_DDM_DIVIDE_THREADS_DIRECT_NONE, &
         KSP_DDM_DIVIDE_THREADS_DIRECT_X, &
         KSP_DDM_DIVIDE_THREADS_DIRECT_Y, &
         KSP_DDM_DIVIDE_THREADS_DIRECT_Z
    use ddmMat_wrap
    use omp_lib
    use param_AS_ddm,only : AdditiveSchwarz_DDM

!    use AdditiveSchwarz_DDM
!    use precon,only : precon_AS
    implicit none

    ! input
    integer, intent(in) ::  nnd, mdist
    integer, intent(in) ::  ncomp_grids, num_grids(ncomp_grids)
    integer, intent(in) ::  ioff_dia(nnd), ioff_grids(nnd,ncomp_grids)
    real*8, intent(in) ::  val_dia(mdist*nnd), vec1(mdist), diag(mdist)
    type(ddm_param), intent(in) :: prm

    ! output
    real*8, intent(inout) :: vec2(mdist)

    ! optional 
    integer,optional :: PRECISION_CALC
    integer,optional :: PRECISION_PRECONDITION
    real*8,optional :: val_dia_lo(mdist*nnd)
    real*8,optional :: vec1_lo(mdist)
    real*8,optional :: vec2_lo(mdist)
    real*8,optional :: diag_lo(mdist)

    ! local
    integer div_direc_th
!    type(precon_AS) precon
    type(precon_AS_wrap) precon
    integer precon_thblock
    integer PRECISION
    
    PRECISION = KSP_DOUBLE_PRECISION
    if(  present(PRECISION_CALC)  )then
       PRECISION = PRECISION_CALC
    endif

    div_direc_th = prm%div_direc_th
    precon_thblock = prm%precon_thblock
    !  write(*,*)"div_direc_th = ",div_direc_th
    if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
       if ( div_direc_th == KSP_DDM_DIVIDE_THREADS_DIRECT_NONE ) then
          write(*,*)"not support"
          stop
       else if ( div_direc_th == KSP_DDM_DIVIDE_THREADS_DIRECT_X ) then
          call ilu_ddm_case_xdiv( mdist, nnd, val_dia, ioff_dia, diag, vec1, vec2, &
               prm, ncomp_grids, ioff_grids, div_direc_th , precon_thblock,&
               PRECISION_CALC, &
               PRECISION_PRECONDITION, &
               val_dia_lo, &
               vec1_lo, &
               diag_lo, &
               vec2_lo &
               )
       else if ( div_direc_th == KSP_DDM_DIVIDE_THREADS_DIRECT_Y ) then
          call ilu_ddm_case_ydiv( mdist, nnd, val_dia, ioff_dia, diag, vec1, vec2, &
               prm, ncomp_grids, ioff_grids, div_direc_th , precon_thblock,&
               PRECISION_CALC, &
               PRECISION_PRECONDITION, &
               val_dia_lo, &
               vec1_lo, &
               diag_lo, &
               vec2_lo &
               )

       else if ( div_direc_th == KSP_DDM_DIVIDE_THREADS_DIRECT_Z ) then
          call ilu_ddm_case_zdiv( mdist, nnd, val_dia, ioff_dia, diag, vec1, vec2, &
               prm, ncomp_grids, ioff_grids, div_direc_th , precon_thblock,&
               PRECISION_CALC, &
               PRECISION_PRECONDITION, &
               val_dia_lo, &
               vec1_lo, &
               diag_lo, &
               vec2_lo &
               )

       else       
          write(*,*) "not support div_direc_th = ",div_direc_th
          stop
       endif

    else
       if ( div_direc_th == KSP_DDM_DIVIDE_THREADS_DIRECT_NONE ) then
          write(*,*)"not support"
          stop
          ! call ilu_ddm_case_nodiv( mdist, nnd, val_dia, ioff_dia, diag, vec1, vec2, &
          !      prm, ncomp_grids, ioff_grids, div_direc_th )
       else if ( div_direc_th == KSP_DDM_DIVIDE_THREADS_DIRECT_X ) then
          call ilu_ddm_case_xdiv( mdist, nnd, val_dia, ioff_dia, diag, vec1, vec2, &
               prm, ncomp_grids, ioff_grids, div_direc_th , precon_thblock)
       else if ( div_direc_th == KSP_DDM_DIVIDE_THREADS_DIRECT_Y ) then
          call ilu_ddm_case_ydiv( mdist, nnd, val_dia, ioff_dia, diag, vec1, vec2, &
               prm, ncomp_grids, ioff_grids, div_direc_th , precon_thblock)
       else if ( div_direc_th == KSP_DDM_DIVIDE_THREADS_DIRECT_Z ) then
          call ilu_ddm_case_zdiv( mdist, nnd, val_dia, ioff_dia, diag, vec1, vec2, &
               prm, ncomp_grids, ioff_grids, div_direc_th , precon_thblock)
       else       
          write(*,*) "not support div_direc_th = ",div_direc_th
          stop
       endif

    endif

  end subroutine aslu_ddm

  subroutine aslu_ddm_th( mdist, nnd, val_dia, ioff_dia, &
       diag, vec1, vec2, prm, ncomp_grids, num_grids, ioff_grids, addL, precon,&
       PRECISION_CALC, &
       PRECISION_PRECONDITION, &
       val_dia_lo, &
       vec1_lo, &
       vec2_lo  &
       )
!       diag, vec1, vec2, prm, ncomp_grids, num_grids, ioff_grids )
    use param_AS_type_ddm
    use param_AS_ddm,only : AdditiveSchwarz_DDM
    use ddmMat_wrap
    use omp_lib
    implicit none
    ! input
    integer, intent(in) ::  nnd, mdist
    integer, intent(in) ::  ncomp_grids, num_grids(ncomp_grids)
    integer, intent(in) ::  ioff_dia(nnd), ioff_grids(nnd,ncomp_grids)
    real*8, intent(in) ::  val_dia(mdist*nnd), vec1(mdist), diag(mdist)
    type(ddm_param), intent(in) :: prm
    ! output
    real*8, intent(inout) :: vec2(mdist)

    ! optional
    integer,optional :: PRECISION_CALC
    integer,optional :: PRECISION_PRECONDITION
    real*8,optional , intent(in) ::  val_dia_lo(mdist*nnd)
    real*8,optional , intent(in) ::  vec1_lo(mdist)
!    real*8,optional , intent(in) ::  diag_lo(mdist)
    real*8,optional , intent(inout) :: vec2_lo(mdist)

    ! local
    integer div_direc_th
    type(precon_AS_wrap) precon
    integer i
    integer addL
    integer PRECISION
    PRECISION = KSP_DOUBLE_PRECISION
    if(  present(PRECISION_CALC)  )then
       PRECISION = PRECISION_CALC
    endif

    div_direc_th = prm%div_direc_th;
    if(PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
       call AdditiveSchwarz_DDM(&
         mdist,nnd, &
         vec1, &
         vec2, &
         prm, &
         precon, &
         ncomp_grids, &
         ioff_dia, &
         val_dia, &
         ioff_grids, &
         addL, &
         PRECISION_CALC, &
         PRECISION_PRECONDITION, &
         val_dia_lo, &
         vec1_lo, &
         vec2_lo  &
         )
    else
       call AdditiveSchwarz_DDM(&
            mdist,nnd, &
            vec1, &
            vec2, &
            prm, &
            precon, &
            ncomp_grids, &
            ioff_dia, &
            val_dia, &
            ioff_grids, &
            addL &
            )
    endif
  end subroutine aslu_ddm_th

#if 0
  subroutine ilu_ddm_case_nodiv( mdist, nnd, val_dia, ioff_dia, &
       diag, vec1, vec2, prm, ncomp_grids, ioff_grids, div_direc_th )
    use ddmMat_wrap
    use omp_lib
    implicit none

    ! input
    integer, intent(in) ::  nnd, mdist, div_direc_th
    integer, intent(in) ::  ncomp_grids
    integer, intent(in) ::  ioff_dia(nnd), ioff_grids(nnd,ncomp_grids)
    real*8, intent(in) ::  val_dia(mdist*nnd), vec1(mdist), diag(mdist)
    type(ddm_param), intent(in) :: prm

    ! output
    real*8, intent(inout) :: vec2(mdist)
    ! local  
    integer k1, j, ip, idiag, offset, istart
    integer nx, ny, nz, mx, my, mz, mxy, jx, jy, jz
    integer stmx, stmy, stmz
    integer nxy, nxyz, i
    integer num_threads, my_thread
    integer nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th
    integer, allocatable :: ioff_x(:), ioff_y(:), ioff_z(:)

    real*8 sum

    !     init
    nx=prm%nx;  ny=prm%ny; nz=prm%nz;
    mx=prm%mx;  my=prm%my; mxy=prm%mxy;
    stmx=prm%stmx;  stmy=prm%stmy; stmz=prm%stmz;
    ! ---
    allocate( ioff_x(nnd) );   allocate( ioff_y(nnd) );   allocate( ioff_z(nnd) );
    ioff_x = 0;  ioff_y = 0;  ioff_z = 0

    if ( ncomp_grids == 1 ) then
       do i=1, nnd
          ioff_x(i) = ioff_grids(i,1)
       end do
    endif
    if ( ncomp_grids == 2 ) then
       do i=1, nnd
          ioff_x(i) = ioff_grids(i,1)
          ioff_y(i) = ioff_grids(i,2)
       end do
    endif
    if ( ncomp_grids == 3 ) then
       do i=1, nnd
          ioff_x(i) = ioff_grids(i,1)
          ioff_y(i) = ioff_grids(i,2)
          ioff_z(i) = ioff_grids(i,3)
       end do
    endif

    !$OMP PARALLEL do
    do i=1, mdist
       vec2(i) = 0.0d0
    end do
    !$OMP end PARALLEL do


    nxs_th = 0;  nxe_th = nx
    nys_th = 0;  nye_th = ny
    nzs_th = 0;  nze_th = nz

    do jz=nzs_th, nze_th -1
       do jy=nys_th, nye_th -1
          do jx=nxs_th, nxe_th -1
             j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;

             sum = 0.0d0
             vec2(j) = 0.0d0

             do ip=1, nnd
                offset = ioff_dia(ip)
                if ( offset > 0 ) cycle
                if ( jx +ioff_x(ip) < nxs_th .or. jx +ioff_x(ip) > nxe_th -1 ) cycle
                if ( jy +ioff_y(ip) < nys_th .or. jy +ioff_y(ip) > nye_th -1 ) cycle
                if ( jz +ioff_z(ip) < nzs_th .or. jz +ioff_z(ip) > nze_th -1 ) cycle
                if (  j +offset < 1 ) cycle

                sum = sum +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
             end do
             vec2(j) = diag(j) *( vec1(j) -sum )
          end do
       end do
    end do

    do jz=nze_th -1, nzs_th, -1
       do jy=nye_th -1, nys_th, -1
          do jx=nxe_th -1, nxs_th, -1
             j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;

             sum = 0.0d0
             do ip=1, nnd
                offset = ioff_dia(ip)
                if ( offset <= 0 ) cycle
                if ( jx +ioff_x(ip) < nxs_th .or. jx +ioff_x(ip) > nxe_th -1 ) cycle
                if ( jy +ioff_y(ip) < nys_th .or. jy +ioff_y(ip) > nye_th -1 ) cycle
                if ( jz +ioff_z(ip) < nzs_th .or. jz +ioff_z(ip) > nze_th -1 ) cycle
                if (  j +offset > mdist ) cycle

                sum = sum +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
             end do

             vec2(j) = vec2(j) -diag(j)*SUM
          end do
       end do
    end do

    deallocate( ioff_x );  deallocate( ioff_y );  deallocate( ioff_z )
    return

  end subroutine ilu_ddm_case_nodiv
#endif

  subroutine ilu_ddm_case_xdiv( mdist, nnd, val_dia, ioff_dia, &
       diag, vec1, vec2, prm, ncomp_grids, ioff_grids, div_direc_th, precon_thblock ,&
       PRECISION_CALC, &
       PRECISION_PRECONDITION, &
       val_dia_lo, &
       vec1_lo, &
       diag_lo, &
       vec2_lo &
       )

    use ddmMat_wrap
    use omp_lib
    use common_kernel,only : set_region_for_threads,get_maxthreads,set_precon_maxthreads
    use PARCEL_define,only : KSP_PRECONDITION_DDM_ILUindepend
    use PARCEL_integer
    implicit none

    ! input
    integer, intent(in) ::  nnd, mdist, div_direc_th, precon_thblock

    integer, intent(in) ::  ncomp_grids
    integer, intent(in) ::  ioff_dia(nnd), ioff_grids(nnd,ncomp_grids)
    real*8, intent(in) ::  val_dia(mdist*nnd), vec1(mdist), diag(mdist)
    type(ddm_param), intent(in) :: prm

    ! output
    real*8, intent(inout) :: vec2(mdist)

    ! optional 
    integer,optional :: PRECISION_CALC
    integer,optional :: PRECISION_PRECONDITION
    real*8,optional :: val_dia_lo(mdist*nnd)
    real*8,optional :: vec1_lo(mdist)
    real*8,optional :: diag_lo(mdist)
    real*8,optional :: vec2_lo(mdist)

    ! local  
    integer k1, j, ip, idiag, offset, istart
    integer(kind=intType) n
    integer nx, ny, nz, mx, my, mz, mxy, jx, jy, jz
    integer stmx, stmy, stmz
    integer nxy, nxyz, i
    integer num_threads, my_thread
    integer max_numth_0
    integer nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th
    integer, allocatable :: ioff_x(:), ioff_y(:), ioff_z(:)
    integer ipL,ipU,ipLn,ipUn
    integer  ip_diaU(nnd),ip_diaL(nnd)
    integer ioffL_xmin
    integer ioffU_xmax
    integer jx0(nnd)
    real*8 sum_lx(prm%nx),sum_lx0
    real*8 sum_lx_lo(prm%nx),sum_lx0_lo

    real*8 ddmul_var1
    real*8 ddadd_var1
    real*8 ddsub_var1
    real*8 c_hi,c_lo

    integer PRECISION
    
    PRECISION = KSP_DOUBLE_PRECISION
    if(  present(PRECISION_CALC)  )then
       PRECISION = PRECISION_CALC
    endif
 
    ! init
    nx=prm%nx;  ny=prm%ny; nz=prm%nz;
    mx=prm%mx;  my=prm%my; mxy=prm%mxy;
    stmx=prm%stmx;  stmy=prm%stmy; stmz=prm%stmz;

    ! ---
    allocate( ioff_x(nnd) );   allocate( ioff_y(nnd) );   allocate( ioff_z(nnd) );
    ioff_x = 0;  ioff_y = 0;  ioff_z = 0

    if ( ncomp_grids == 1 ) then
       do i=1, nnd
          ioff_x(i) = ioff_grids(i,1)
       end do
    endif
    if ( ncomp_grids == 2 ) then
       do i=1, nnd
          ioff_x(i) = ioff_grids(i,1)
          ioff_y(i) = ioff_grids(i,2)
       end do
    endif
    if ( ncomp_grids == 3 ) then
       do i=1, nnd
          ioff_x(i) = ioff_grids(i,1)
          ioff_y(i) = ioff_grids(i,2)
          ioff_z(i) = ioff_grids(i,3)
       end do
    endif

    do ip=1, nnd
       ip_diaU(ip) = 0
       ip_diaL(ip) = 0
    enddo
    ipL = 1
    ipU = 1
    ipLn = 0
    ipUn = 0
    do ip=1, nnd
       offset = ioff_dia(ip)
       if(offset.eq.0) cycle
       if ( offset < 0 )then
          ip_diaL(ipL) = ip
          ipL  = ipL  + 1
          ipLn = ipLn + 1
       else
          ip_diaU(ipU) = ip
          ipU  = ipU  + 1
          ipUn = ipUn + 1
       endif
    enddo

    ioffL_xmin = 0
    do ipL=1,ipLn
       ip = ip_diaL(ipL)
       if(ioffL_xmin > ioff_x(ip))then
          ioffL_xmin = ioff_x(ip)
       endif
    enddo

    ioffU_xmax = 0
    do ipU=1,ipUn
       ip = ip_diaU(ipU)
       if ( ioff_x(ip) > ioffU_xmax )then
          ioffU_xmax = ioff_x(ip)
       endif
    enddo

    call get_maxthreads(max_numth_0) 
    n = mdist
    call set_precon_maxthreads(n,precon_thblock)

    !$OMP parallel private( num_threads, my_thread ) &
    !$OMP private( i, jz, jy, jx, jx0, j, sum_lx,sum_lx0, offset, ip, ipL  ) &
    !$OMP private( nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th ) &
    !$OMP private( &
    !$OMP ddmul_var1, &
    !$OMP ddadd_var1, &
    !$OMP ddsub_var1, &
    !$OMP sum_lx_lo,sum_lx0_lo, c_hi,c_lo )

    num_threads = omp_get_num_threads();
    my_thread = omp_get_thread_num();

    call set_region_for_threads( my_thread, num_threads, nx, ny, nz, &
         nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th, &
         div_direc_th )

    if(prm%iflag_DDM_independ.eq.KSP_PRECONDITION_DDM_ILUindepend)then
       jx0 = -1
       do ipL=1,ipLn
          ip = ip_diaL(ipL)
          offset = ioff_dia(ip)
          do jx=nxs_th, nxe_th -1
             if ( jx +ioff_x(ip) < nxs_th ) then
                cycle
             endif
             jx0(ipL) = jx 
             exit
          enddo
       enddo

       do jz=nzs_th, nze_th -1
          do jy=nys_th, nye_th -1
             do jx=nxs_th, nxe_th -1
                sum_lx(jx+1) = 0.0d0
             enddo
             if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                do jx=nxs_th, nxe_th -1
                   sum_lx_lo(jx+1) = 0.0d0
                enddo
             endif
             do ipL=1,ipLn
                ip = ip_diaL(ipL)
                offset = ioff_dia(ip)
                if( (ioff_y(ip).eq.0).and.(ioff_z(ip).eq.0) )then
                   cycle
                endif
!                do jx=nxs_th, nxe_th -1
                if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                   if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                      do jx=jx0(ipL), nxe_th -1
                         j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                         ! sum_lx(jx+1) = sum_lx(jx+1) +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                         ddmul_calc1(val_dia(mdist*(ip-1)+j),val_dia_lo(mdist*(ip-1)+j),vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo)
                         ddadd_calc1(c_hi,c_lo,sum_lx(jx+1),sum_lx_lo(jx+1),sum_lx(jx+1),sum_lx_lo(jx+1))
                      end do
                      cycle
                   endif
                   do jx=jx0(ipL), nxe_th -1
                      j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                      ! sum_lx(jx+1) = sum_lx(jx+1) +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                      ddmul_calc1(val_dia(mdist*(ip-1)+j),0.0d0,vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo)
                      ddadd_calc1(c_hi,c_lo,sum_lx(jx+1),sum_lx_lo(jx+1),sum_lx(jx+1),sum_lx_lo(jx+1))
                   end do
                else
                   do jx=jx0(ipL), nxe_th -1
                      j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                      sum_lx(jx+1) = sum_lx(jx+1) +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                   end do
                endif
             enddo
             if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                   do jx=nxs_th, nxe_th -1
                      j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                      ! vec2(j) = diag(j) *( vec1(j) -sum_lx(jx+1) )
                      ddsub_calc1(vec1(j),vec1_lo(j),sum_lx(jx+1),sum_lx_lo(jx+1),c_hi,c_lo)
                      ddmul_calc1(diag(j),diag_lo(j),c_hi,c_lo,vec2(j),vec2_lo(j))
                   end do
                   cycle
                endif
                do jx=nxs_th, nxe_th -1
                   j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                   ! vec2(j) = diag(j) *( vec1(j) -sum_lx(jx+1) )
                   ddsub_calc1(vec1(j),vec1_lo(j),sum_lx(jx+1),sum_lx_lo(jx+1),c_hi,c_lo)
                   ddmul_calc1(diag(j),0.0d0,c_hi,c_lo,vec2(j),vec2_lo(j))
                end do
             else
                do jx=nxs_th, nxe_th -1
                   j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                   vec2(j) = diag(j) *( vec1(j) -sum_lx(jx+1) )
                end do
             endif

          end do
       end do
    else

       do jz=nzs_th, nze_th -1
          do jy=nys_th, nye_th -1
             do jx=nxs_th, nxe_th -1
                j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                sum_lx0 = 0.0d0
                if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                   sum_lx0_lo = 0.0d0
                endif
                do ipL=1,ipLn
                   ip = ip_diaL(ipL)
                   offset = ioff_dia(ip)
                   !                   if ( jx +ioff_x(ip) < nxs_th .or. jx +ioff_x(ip) > nxe_th -1 ) cycle
                   if ( jx +ioff_x(ip) < nxs_th ) then
                      cycle
                   endif
                   if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                      if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                         ! sum_lx0 = sum_lx0 +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                         ddmul_calc1(val_dia(mdist*(ip-1)+j),val_dia_lo(mdist*(ip-1)+j),vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo)
                         ddadd_calc1(sum_lx0,sum_lx0_lo,c_hi,c_lo,sum_lx0,sum_lx0_lo)
                         cycle
                      endif
                      ! sum_lx0 = sum_lx0 +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                      ddmul_calc1(val_dia(mdist*(ip-1)+j),0.0d0,vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo)
                      ddadd_calc1(sum_lx0,sum_lx0_lo,c_hi,c_lo,sum_lx0,sum_lx0_lo)
                   else
                      sum_lx0 = sum_lx0 +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                   endif
                end do
                if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                   if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                      ! vec2(j) = diag(j) *( vec1(j) -sum_lx0 )
                      ddsub_calc1(vec1(j),vec1_lo(j),sum_lx0,sum_lx0_lo,c_hi,c_lo)
                      ddmul_calc1(diag(j),diag_lo(j),c_hi,c_lo,vec2(j),vec2_lo(j))
                      cycle
                   endif
                   ! vec2(j) = diag(j) *( vec1(j) -sum_lx0 )
                   ddsub_calc1(vec1(j),vec1_lo(j),sum_lx0,sum_lx0_lo,c_hi,c_lo)
                   ddmul_calc1(diag(j),0.0d0,c_hi,c_lo,vec2(j),vec2_lo(j))
                else
                   vec2(j) = diag(j) *( vec1(j) -sum_lx0 )
                endif
             end do
          end do
       end do
    endif
    !$OMP end parallel

    !$OMP Parallel private( num_threads, my_thread ) &
    !$OMP private( i, jz, jy, jx, jx0, j, sum_lx,sum_lx0, offset, ip, ipU ) &
    !$OMP private( nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th ) &
    !$OMP private( &
    !$OMP ddmul_var1, &
    !$OMP ddadd_var1, &
    !$OMP ddsub_var1, &
    !$OMP sum_lx_lo,sum_lx0_lo, c_hi,c_lo )

    num_threads = omp_get_num_threads();
    my_thread = omp_get_thread_num();

    call set_region_for_threads( my_thread, num_threads, nx, ny, nz, &
         nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th, &
         div_direc_th )
    if(prm%iflag_DDM_independ.eq.KSP_PRECONDITION_DDM_ILUindepend)then

       jx0 = -1
       do ipU=1,ipUn
          ip = ip_diaU(ipU)
          offset = ioff_dia(ip)
          do jx=nxe_th -1, nxs_th, -1
             if ( jx +ioff_x(ip) > nxe_th -1 ) then
                cycle
             endif
             jx0(ipU) = jx 
             exit
          enddo
       enddo
       
       do jz=nze_th -1, nzs_th, -1
          do jy=nye_th -1, nys_th, -1
             do jx=nxe_th -1, nxs_th, -1
                sum_lx(jx+1) = 0.0d0
             enddo
             if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                do jx=nxs_th, nxe_th -1
                   sum_lx_lo(jx+1) = 0.0d0
                enddo
             endif
             do ipU=1,ipUn
                ip = ip_diaU(ipU)
                offset = ioff_dia(ip)
                if( (ioff_y(ip).eq.0).and.(ioff_z(ip).eq.0) )then
                   cycle
                endif
                if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                   if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                      do jx=jx0(ipU), nxs_th, -1
                         j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                        ! sum_lx(jx+1) = sum_lx(jx+1) +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                         ddmul_calc1( val_dia(mdist*(ip-1)+j),val_dia_lo(mdist*(ip-1)+j),vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo )
                         ddadd_calc1( sum_lx(jx+1),sum_lx_lo(jx+1),c_hi,c_lo,sum_lx(jx+1),sum_lx_lo(jx+1))
                      end do
                      cycle
                   endif
                   do jx=jx0(ipU), nxs_th, -1
                      j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                      ! sum_lx(jx+1) = sum_lx(jx+1) +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                      ddmul_calc1( val_dia(mdist*(ip-1)+j),0.0d0,vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo )
                      ddadd_calc1( sum_lx(jx+1),sum_lx_lo(jx+1),c_hi,c_lo,sum_lx(jx+1),sum_lx_lo(jx+1))
                   end do
                else
                   do jx=jx0(ipU), nxs_th, -1
                      j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                      sum_lx(jx+1) = sum_lx(jx+1) +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                   end do
                endif
             end do
             if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                   do jx=nxe_th -1, nxs_th, -1
                      j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                      ! vec2(j) = vec2(j) -diag(j)*sum_lx(jx+1)
                      ddmul_calc1(diag(j),diag_lo(j),sum_lx(jx+1),sum_lx_lo(jx+1),c_hi,c_lo)
                      ddsub_calc1(vec2(j),vec2_lo(j),c_hi,c_lo,vec2(j),vec2_lo(j))
                   end do
                   cycle
                endif
                do jx=nxe_th -1, nxs_th, -1
                   j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                   !  vec2(j) = vec2(j) -diag(j)*sum_lx(jx+1)
                   ddmul_calc1(diag(j),0.0d0,sum_lx(jx+1),sum_lx_lo(jx+1),c_hi,c_lo)
                   ddsub_calc1(vec2(j),vec2_lo(j),c_hi,c_lo,vec2(j),vec2_lo(j))
                end do
             else
                do jx=nxe_th -1, nxs_th, -1
                   j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                   vec2(j) = vec2(j) -diag(j)*sum_lx(jx+1)
                end do
             endif
          end do
       end do
    else
       do jz=nze_th -1, nzs_th, -1
          do jy=nye_th -1, nys_th, -1
             do jx=nxe_th -1, nxs_th, -1
                j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                sum_lx0 = 0.0d0
                if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                   sum_lx0_lo = 0.0d0
                endif
                do ipU=1,ipUn
                   ip = ip_diaU(ipU)
                   offset = ioff_dia(ip)
                   !                if ( jx +ioff_x(ip) < nxs_th .or. jx +ioff_x(ip) > nxe_th -1 ) cycle
                   if ( jx +ioff_x(ip) > nxe_th -1 ) then
                      cycle
                   endif
                   if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                      if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                         ! sum_lx0 = sum_lx0 +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                         ddmul_calc1(val_dia(mdist*(ip-1)+j),val_dia_lo(mdist*(ip-1)+j),vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo)
                         ddadd_calc1(sum_lx0,sum_lx0_lo,c_hi,c_lo,sum_lx0,sum_lx0_lo)
                         cycle
                      endif
                      ! sum_lx0 = sum_lx0 +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)                   
                      ddmul_calc1(val_dia(mdist*(ip-1)+j),0.0d0,vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo)
                      ddadd_calc1(sum_lx0,sum_lx0_lo,c_hi,c_lo,sum_lx0,sum_lx0_lo)
                   else
                     sum_lx0 = sum_lx0 +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                   endif
                end do
                   if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                      if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                         ! vec2(j) = vec2(j) -diag(j)*sum_lx0
                         ddmul_calc1(diag(j),diag_lo(j),sum_lx0,sum_lx0_lo,c_hi,c_lo)
                         ddsub_calc1(vec2(j),vec2_lo(j),c_hi,c_lo,vec2(j),vec2_lo(j))
                         cycle
                      endif
                      ! vec2(j) = vec2(j) -diag(j)*sum_lx0
                      ddmul_calc1(diag(j),0.0d00,sum_lx0,sum_lx0_lo,c_hi,c_lo)
                      ddsub_calc1(vec2(j),vec2_lo(j),c_hi,c_lo,vec2(j),vec2_lo(j))
                   else
                      vec2(j) = vec2(j) -diag(j)*sum_lx0
                   endif
             end do
          end do
       end do
    endif
    !$OMP end Parallel

    deallocate( ioff_x );  deallocate( ioff_y );  deallocate( ioff_z )

    call omp_set_num_threads(max_numth_0)
    
  end subroutine ilu_ddm_case_xdiv

  subroutine ilu_ddm_case_ydiv( mdist, nnd, val_dia, ioff_dia, &
       diag, vec1, vec2, prm, ncomp_grids, ioff_grids, div_direc_th, precon_thblock, &
       PRECISION_CALC, &
       PRECISION_PRECONDITION, &
       val_dia_lo, &
       vec1_lo, &
       diag_lo, &
       vec2_lo &
       )


    use ddmMat_wrap
    use omp_lib
    use common_kernel,only : set_region_for_threads,get_maxthreads,set_precon_maxthreads
    use PARCEL_define,only : KSP_PRECONDITION_DDM_ILUindepend
    use PARCEL_integer
    implicit none

    ! input
    integer, intent(in) ::  nnd, mdist, div_direc_th, precon_thblock
    integer, intent(in) ::  ncomp_grids
    integer, intent(in) ::  ioff_dia(nnd), ioff_grids(nnd,ncomp_grids)
    real*8, intent(in) ::  val_dia(mdist*nnd), vec1(mdist), diag(mdist)
    type(ddm_param), intent(in) :: prm

    ! output
    real*8, intent(inout) :: vec2(mdist)

    ! optional 
    integer,optional :: PRECISION_CALC
    integer,optional :: PRECISION_PRECONDITION
    real*8,optional :: val_dia_lo(mdist*nnd)
    real*8,optional :: vec1_lo(mdist)
    real*8,optional :: diag_lo(mdist)
    real*8,optional :: vec2_lo(mdist)

    ! local  
    integer k1, j, ip, idiag, offset, istart
    integer(kind=intType) n
    integer nx, ny, nz, mx, my, mz, mxy, jx, jy, jz
    integer stmx, stmy, stmz
    integer nxy, nxyz, i
    integer num_threads, my_thread
    integer max_numth_0
    integer nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th
    integer, allocatable :: ioff_x(:), ioff_y(:), ioff_z(:)
    integer ipL,ipU,ipLn,ipUn
    integer  ip_diaU(nnd),ip_diaL(nnd)
    integer ioffL_ymin
    integer ioffU_ymax

    real*8 sum_lx(prm%nx),sum_lx0
    real*8 sum_lx_lo(prm%nx),sum_lx0_lo

    real*8 ddmul_var1
    real*8 ddadd_var1
    real*8 ddsub_var1
    real*8 c_hi,c_lo

    integer PRECISION
    
    PRECISION = KSP_DOUBLE_PRECISION
    if(  present(PRECISION_CALC)  )then
       PRECISION = PRECISION_CALC
    endif
 

    !     init
    nx=prm%nx;  ny=prm%ny; nz=prm%nz;
    mx=prm%mx;  my=prm%my; mxy=prm%mxy;
    stmx=prm%stmx;  stmy=prm%stmy; stmz=prm%stmz;

    ! ---
    allocate( ioff_x(nnd) );   allocate( ioff_y(nnd) );   allocate( ioff_z(nnd) );
    ioff_x = 0;  ioff_y = 0;  ioff_z = 0

    if ( ncomp_grids == 1 ) then
       do i=1, nnd
          ioff_x(i) = ioff_grids(i,1)
       end do
    endif
    if ( ncomp_grids == 2 ) then
       do i=1, nnd
          ioff_x(i) = ioff_grids(i,1)
          ioff_y(i) = ioff_grids(i,2)
       end do
    endif
    if ( ncomp_grids == 3 ) then
       do i=1, nnd
          ioff_x(i) = ioff_grids(i,1)
          ioff_y(i) = ioff_grids(i,2)
          ioff_z(i) = ioff_grids(i,3)
       end do
    endif

    do ip=1, nnd
       ip_diaU(ip) = 0
       ip_diaL(ip) = 0
    enddo
    ipL = 1
    ipU = 1
    ipLn = 0
    ipUn = 0
    do ip=1, nnd
       offset = ioff_dia(ip)
       if(offset.eq.0) cycle
       if ( offset < 0 )then
          ip_diaL(ipL) = ip
          ipL  = ipL  + 1
          ipLn = ipLn + 1
       else
          ip_diaU(ipU) = ip
          ipU  = ipU  + 1
          ipUn = ipUn + 1
       endif
    enddo

    ioffL_ymin = 0
    do ipL=1,ipLn
       ip = ip_diaL(ipL)
       if(ioffL_ymin > ioff_y(ip))then
          ioffL_ymin = ioff_y(ip)
       endif
    enddo

    ioffU_ymax = 0
    do ipU=1,ipUn
       ip = ip_diaU(ipU)
       if ( ioff_y(ip) > ioffU_ymax )then
          ioffU_ymax = ioff_y(ip)
       endif
    enddo


    call get_maxthreads(max_numth_0) 
    n = mdist
    call set_precon_maxthreads(n,precon_thblock)
    !$OMP Parallel private( num_threads, my_thread ) &
    !$OMP private( i, jz, jy, jx, j, sum_lx, sum_lx0, offset, ip, ipL ) &
    !$OMP private( nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th ) &
    !$OMP private( &
    !$OMP ddmul_var1, &
    !$OMP ddadd_var1, &
    !$OMP ddsub_var1, &
    !$OMP sum_lx_lo,sum_lx0_lo, c_hi,c_lo )

    num_threads = omp_get_num_threads();
    my_thread = omp_get_thread_num();

    call set_region_for_threads( my_thread, num_threads, nx, ny, nz, &
                                nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th, &
                                div_direc_th )

    if(prm%iflag_DDM_independ.eq.KSP_PRECONDITION_DDM_ILUindepend)then
       do jz=nzs_th, nze_th -1
          do jy=nys_th, nye_th -1
             do jx=nxs_th, nxe_th -1
                sum_lx(jx+1) = 0.0d0
             enddo
             if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                do jx=nxs_th, nxe_th -1
                   sum_lx_lo(jx+1) = 0.0d0
                enddo
             endif
             do ipL=1,ipLn
                ip = ip_diaL(ipL)
                offset = ioff_dia(ip)
                if( (ioff_y(ip).eq.0).and.(ioff_z(ip).eq.0) )then
                   cycle
                endif
!                if ( jy +ioff_y(ip) < nys_th .or. jy +ioff_y(ip) > nye_th -1 ) then
                if ( jy +ioff_y(ip) < nys_th) then
                   cycle
                endif
                if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                   if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                      do jx=nxs_th, nxe_th -1
                         j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                         ! sum_lx(jx+1) = sum_lx(jx+1) +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                         ddmul_calc1(val_dia(mdist*(ip-1)+j),val_dia_lo(mdist*(ip-1)+j),vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo)
                         ddadd_calc1(c_hi,c_lo,sum_lx(jx+1),sum_lx_lo(jx+1),sum_lx(jx+1),sum_lx_lo(jx+1))
                      end do
                      cycle
                   endif
                   do jx=nxs_th, nxe_th -1
                      j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                      ! sum_lx(jx+1) = sum_lx(jx+1) +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                      ddmul_calc1(val_dia(mdist*(ip-1)+j),0.0d0,vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo)
                      ddadd_calc1(c_hi,c_lo,sum_lx(jx+1),sum_lx_lo(jx+1),sum_lx(jx+1),sum_lx_lo(jx+1))
                   end do
                else
                   do jx=nxs_th, nxe_th -1
                      j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                      sum_lx(jx+1) = sum_lx(jx+1) +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                   end do
                endif
             end do
             if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                   do jx=nxs_th, nxe_th -1
                      j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                      ! vec2(j) = diag(j) *( vec1(j) -sum_lx(jx+1) )
                      ddsub_calc1(vec1(j),vec1_lo(j),sum_lx(jx+1),sum_lx_lo(jx+1),c_hi,c_lo)
                      ddmul_calc1(diag(j),diag_lo(j),c_hi,c_lo,vec2(j),vec2_lo(j))
                   end do
                   cycle
                endif
                do jx=nxs_th, nxe_th -1
                   j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                   ! vec2(j) = diag(j) *( vec1(j) -sum_lx(jx+1) )
                   ddsub_calc1(vec1(j),vec1_lo(j),sum_lx(jx+1),sum_lx_lo(jx+1),c_hi,c_lo)
                   ddmul_calc1(diag(j),0.0d0,c_hi,c_lo,vec2(j),vec2_lo(j))
                end do
             else
                do jx=nxs_th, nxe_th -1
                   j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                   vec2(j) = diag(j) *( vec1(j) -sum_lx(jx+1) )
                end do
             endif
          end do
       end do
    else
       do jz=nzs_th, nze_th -1
          do jy=nys_th, nye_th -1
             do jx=nxs_th, nxe_th -1
                j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                sum_lx0 = 0.0d0
                if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                   sum_lx0_lo = 0.0d0
                endif
                do ipL=1,ipLn
                   ip = ip_diaL(ipL)
                   offset = ioff_dia(ip)
!                   if ( jy +ioff_y(ip) < nys_th .or. jy +ioff_y(ip) > nye_th -1 ) cycle
                   if ( jy +ioff_y(ip) < nys_th) then
                      cycle
                   endif
                   if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                      if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                         ddmul_calc1(val_dia(mdist*(ip-1)+j),val_dia_lo(mdist*(ip-1)+j),vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo)
                         ddadd_calc1(sum_lx0,sum_lx0_lo,c_hi,c_lo,sum_lx0,sum_lx0_lo)
                         cycle
                      endif
                      ddmul_calc1(val_dia(mdist*(ip-1)+j),0.0d0,vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo)
                      ddadd_calc1(sum_lx0,sum_lx0_lo,c_hi,c_lo,sum_lx0,sum_lx0_lo)
                   else
                      sum_lx0 = sum_lx0 +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                   endif

                end do
                if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                   if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                      ddsub_calc1(vec1(j),vec1_lo(j),sum_lx0,sum_lx0_lo,c_hi,c_lo)
                      ddmul_calc1(diag(j),diag_lo(j),c_hi,c_lo,vec2(j),vec2_lo(j))
                      cycle
                   endif
                   ddsub_calc1(vec1(j),vec1_lo(j),sum_lx0,sum_lx0_lo,c_hi,c_lo)
                   ddmul_calc1(diag(j),0.0d0,c_hi,c_lo,vec2(j),vec2_lo(j))
                else
                   vec2(j) = diag(j) *( vec1(j) -sum_lx0 )
                endif
             end do
          end do
       end do
    endif
    !$OMP end parallel

    !$OMP Parallel private( num_threads, my_thread ) &
    !$OMP private( i, jz, jy, jx, j, sum_lx, sum_lx0, offset, ip, ipU ) &
    !$OMP private( nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th ) &
    !$OMP private( &
    !$OMP ddmul_var1, &
    !$OMP ddadd_var1, &
    !$OMP ddsub_var1, &
    !$OMP sum_lx_lo,sum_lx0_lo, c_hi,c_lo )

    num_threads = omp_get_num_threads();
    my_thread = omp_get_thread_num();

    call set_region_for_threads( my_thread, num_threads, nx, ny, nz, &
                                nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th, &
                                div_direc_th )
    if(prm%iflag_DDM_independ.eq.KSP_PRECONDITION_DDM_ILUindepend)then
       do jz=nze_th -1, nzs_th, -1
          do jy=nye_th -1, nys_th, -1
             do jx=nxe_th -1, nxs_th, -1
                sum_lx(jx+1) = 0.0d0
             enddo
             if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                do jx=nxs_th, nxe_th -1
                   sum_lx_lo(jx+1) = 0.0d0
                enddo
             endif
             do ipU=1,ipUn
                ip = ip_diaU(ipU)
                offset = ioff_dia(ip)
                if( (ioff_y(ip).eq.0).and.(ioff_z(ip).eq.0) )then
                   cycle
                endif
!                if ( jy +ioff_y(ip) < nys_th .or. jy +ioff_y(ip) > nye_th -1 ) then
                if ( jy +ioff_y(ip) > nye_th -1 ) then
                   cycle
                endif
                if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                   if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                      do jx=nxe_th -1, nxs_th, -1
                         j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                         ! sum_lx(jx+1) = sum_lx(jx+1) +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                         ddmul_calc1( val_dia(mdist*(ip-1)+j),val_dia_lo(mdist*(ip-1)+j),vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo )
                         ddadd_calc1( sum_lx(jx+1),sum_lx_lo(jx+1),c_hi,c_lo,sum_lx(jx+1),sum_lx_lo(jx+1))
                      end do
                      cycle
                   endif
                   do jx=nxe_th -1, nxs_th, -1
                      j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                      ! sum_lx(jx+1) = sum_lx(jx+1) +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                      ddmul_calc1( val_dia(mdist*(ip-1)+j),0.0d0,vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo )
                      ddadd_calc1( sum_lx(jx+1),sum_lx_lo(jx+1),c_hi,c_lo,sum_lx(jx+1),sum_lx_lo(jx+1))
                   end do
                else
                   do jx=nxe_th -1, nxs_th, -1
                      j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                      sum_lx(jx+1) = sum_lx(jx+1) +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                   end do
                endif
             enddo
             if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                   do jx=nxe_th -1, nxs_th, -1
                      j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                      ! vec2(j) = vec2(j) -diag(j)*sum_lx(jx+1)
                      ddmul_calc1(diag(j),diag_lo(j),sum_lx(jx+1),sum_lx_lo(jx+1),c_hi,c_lo)
                      ddsub_calc1(vec2(j),vec2_lo(j),c_hi,c_lo,vec2(j),vec2_lo(j))
                   end do
                   cycle
                endif
                do jx=nxe_th -1, nxs_th, -1
                   j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                   ! vec2(j) = vec2(j) -diag(j)*sum_lx(jx+1)
                   ddmul_calc1(diag(j),0.0d0,sum_lx(jx+1),sum_lx_lo(jx+1),c_hi,c_lo)
                   ddsub_calc1(vec2(j),vec2_lo(j),c_hi,c_lo,vec2(j),vec2_lo(j))
                end do
             else
                do jx=nxe_th -1, nxs_th, -1
                   j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                   vec2(j) = vec2(j) -diag(j)*sum_lx(jx+1)
                end do
             endif
          end do
       end do
    else
       do jz=nze_th -1, nzs_th, -1
          do jy=nye_th -1, nys_th, -1
             do jx=nxe_th -1, nxs_th, -1
                j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                sum_lx0 = 0.0d0
                if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                   sum_lx0_lo = 0.0d0
                endif
                do ipU=1,ipUn
                   ip = ip_diaU(ipU)
                   offset = ioff_dia(ip)
                   if ( jy +ioff_y(ip) > nye_th -1 ) then
                      cycle
                   endif
                   if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                      if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                         ddmul_calc1(val_dia(mdist*(ip-1)+j),val_dia_lo(mdist*(ip-1)+j),vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo)
                         ddadd_calc1(sum_lx0,sum_lx0_lo,c_hi,c_lo,sum_lx0,sum_lx0_lo)
                         cycle 
                      endif
                      ddmul_calc1(val_dia(mdist*(ip-1)+j),0.0d0,vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo)
                      ddadd_calc1(sum_lx0,sum_lx0_lo,c_hi,c_lo,sum_lx0,sum_lx0_lo)
                   else
                      sum_lx0 = sum_lx0 +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                   endif
                end do
                if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                   if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                      ! vec2(j) = vec2(j) -diag(j)*sum_lx0
                      ddmul_calc1(diag(j),diag_lo(j),sum_lx0,sum_lx0_lo,c_hi,c_lo)
                      ddsub_calc1(vec2(j),vec2_lo(j),c_hi,c_lo,vec2(j),vec2_lo(j))
                      cycle
                   endif
                   ! vec2(j) = vec2(j) -diag(j)*sum_lx0
                   ddmul_calc1(diag(j),0.0d00,sum_lx0,sum_lx0_lo,c_hi,c_lo)
                   ddsub_calc1(vec2(j),vec2_lo(j),c_hi,c_lo,vec2(j),vec2_lo(j))
                else
                   vec2(j) = vec2(j) -diag(j)*sum_lx0
                endif
             end do
          end do
       end do
    endif
    !$OMP end Parallel

    deallocate( ioff_x );  deallocate( ioff_y );  deallocate( ioff_z )

    call omp_set_num_threads(max_numth_0)

  end subroutine ilu_ddm_case_ydiv

  subroutine ilu_ddm_case_zdiv( mdist, nnd, val_dia, ioff_dia, &
       diag, vec1, vec2, prm, ncomp_grids, ioff_grids, div_direc_th, precon_thblock,&
       PRECISION_CALC, &
       PRECISION_PRECONDITION, &
       val_dia_lo, &
       vec1_lo, &
       diag_lo, &
       vec2_lo &
       )

    use ddmMat_wrap
    use omp_lib
    use common_kernel,only : set_region_for_threads,get_maxthreads,set_precon_maxthreads
    use PARCEL_define,only : KSP_PRECONDITION_DDM_ILUindepend
    use PARCEL_integer
    implicit none

    ! input
    integer, intent(in) ::  nnd, mdist, div_direc_th, precon_thblock
    integer, intent(in) ::  ncomp_grids
    integer, intent(in) ::  ioff_dia(nnd), ioff_grids(nnd,ncomp_grids)
    real*8, intent(in) ::  val_dia(mdist*nnd), vec1(mdist), diag(mdist)
    type(ddm_param), intent(in) :: prm

    ! output
    real*8, intent(inout) :: vec2(mdist)

    ! optional 
    integer,optional :: PRECISION_CALC
    integer,optional :: PRECISION_PRECONDITION
    real*8,optional :: val_dia_lo(mdist*nnd)
    real*8,optional :: vec1_lo(mdist)
    real*8,optional :: diag_lo(mdist)
    real*8,optional :: vec2_lo(mdist)

    ! local  
    integer k1, j, ip, idiag, offset, istart
    integer(kind=intType) n
    integer nx, ny, nz, mx, my, mz, mxy, jx, jy, jz
    integer stmx, stmy, stmz
    integer nxy, nxyz, i
    integer num_threads, my_thread
    integer max_numth_0
    integer nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th
    integer, allocatable :: ioff_x(:), ioff_y(:), ioff_z(:)
    integer ipL,ipU,ipLn,ipUn
    integer  ip_diaU(nnd),ip_diaL(nnd)
    integer ioffL_zmin
    integer ioffU_zmax

    real*8 sum_lx(prm%nx),sum_lx0
    real*8 sum_lx_lo(prm%nx),sum_lx0_lo

    real*8 ddmul_var1
    real*8 ddadd_var1
    real*8 ddsub_var1
    real*8 c_hi,c_lo

    integer PRECISION
    
    PRECISION = KSP_DOUBLE_PRECISION
    if(  present(PRECISION_CALC)  )then
       PRECISION = PRECISION_CALC
    endif

    ! init
    nx=prm%nx;  ny=prm%ny; nz=prm%nz;
    mx=prm%mx;  my=prm%my; mxy=prm%mxy;
    stmx=prm%stmx;  stmy=prm%stmy; stmz=prm%stmz;

    ! ---
    allocate( ioff_x(nnd) );   allocate( ioff_y(nnd) );   allocate( ioff_z(nnd) );
    ioff_x = 0;  ioff_y = 0;  ioff_z = 0

    if ( ncomp_grids == 1 ) then
       do i=1, nnd
          ioff_x(i) = ioff_grids(i,1)
       end do
    endif
    if ( ncomp_grids == 2 ) then
       do i=1, nnd
          ioff_x(i) = ioff_grids(i,1)
          ioff_y(i) = ioff_grids(i,2)
       end do
    endif
    if ( ncomp_grids == 3 ) then
       do i=1, nnd
          ioff_x(i) = ioff_grids(i,1)
          ioff_y(i) = ioff_grids(i,2)
          ioff_z(i) = ioff_grids(i,3)
       end do
    endif

    do ip=1, nnd
       ip_diaU(ip) = 0
       ip_diaL(ip) = 0
    enddo
    ipL = 1
    ipU = 1
    ipLn = 0
    ipUn = 0
    do ip=1, nnd
       offset = ioff_dia(ip)
       if(offset.eq.0) cycle
       if ( offset < 0 )then
          ip_diaL(ipL) = ip
          ipL  = ipL  + 1
          ipLn = ipLn + 1
       else
          ip_diaU(ipU) = ip
          ipU  = ipU  + 1
          ipUn = ipUn + 1
       endif
    enddo

    ioffL_zmin = 0
    do ipL=1,ipLn
       ip = ip_diaL(ipL)
       if(ioffL_zmin > ioff_z(ip))then
          ioffL_zmin = ioff_z(ip)
       endif
    enddo

    ioffU_zmax = 0
    do ipU=1,ipUn
       ip = ip_diaU(ipU)
       if ( ioff_z(ip) > ioffU_zmax )then
          ioffU_zmax = ioff_z(ip)
       endif
    enddo

    call get_maxthreads(max_numth_0) 
    n = mdist
    call set_precon_maxthreads(n,precon_thblock)
    !$OMP Parallel private( num_threads, my_thread ) &
    !$OMP private( i, jz, jy, jx, j, sum_lx,sum_lx0, offset, ip, ipL ) &
    !$OMP private( nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th )&
    !$OMP private( &
    !$OMP ddmul_var1, &
    !$OMP ddadd_var1, &
    !$OMP ddsub_var1, &
    !$OMP sum_lx_lo,sum_lx0_lo, c_hi,c_lo )

    num_threads = omp_get_num_threads();
    my_thread = omp_get_thread_num();

    call set_region_for_threads(  &
         my_thread, num_threads, nx, ny, nz, &
         nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th, &
         div_direc_th )

    if(prm%iflag_DDM_independ.eq.KSP_PRECONDITION_DDM_ILUindepend)then
       do jz=nzs_th , nze_th -1
          do jy=nys_th, nye_th -1
             do jx=nxs_th, nxe_th -1
                sum_lx(jx+1) = 0.0d0
             enddo
             if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                do jx=nxs_th, nxe_th -1
                   sum_lx_lo(jx+1) = 0.0d0
                enddo
             endif

             do ipL=1,ipLn
                ip     = ip_diaL(ipL)
                offset = ioff_dia(ip)
                if( (ioff_y(ip).eq.0).and.(ioff_z(ip).eq.0) )then
                   cycle
                endif
                if ( jz + ioff_z(ip) < nzs_th) then
                   cycle
                endif
                if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                   if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                      do jx=nxs_th, nxe_th -1
                         j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                         ! sum_lx(jx+1) = sum_lx(jx+1) +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                         ddmul_calc1(val_dia(mdist*(ip-1)+j),val_dia_lo(mdist*(ip-1)+j),vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo)
                         ddadd_calc1(c_hi,c_lo,sum_lx(jx+1),sum_lx_lo(jx+1),sum_lx(jx+1),sum_lx_lo(jx+1))
                      end do
                      cycle
                   endif
                   do jx=nxs_th, nxe_th -1
                      j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                     ! sum_lx(jx+1) = sum_lx(jx+1) +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                      ddmul_calc1(val_dia(mdist*(ip-1)+j),0.0d0,vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo)
                      ddadd_calc1(c_hi,c_lo,sum_lx(jx+1),sum_lx_lo(jx+1),sum_lx(jx+1),sum_lx_lo(jx+1))
                   end do
                else
                   do jx=nxs_th, nxe_th -1
                      j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                      sum_lx(jx+1) = sum_lx(jx+1) +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                   end do
                endif
             enddo
             if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                   do jx=nxs_th, nxe_th -1
                      j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                      ! vec2(j) = diag(j) *( vec1(j) - sum_lx(jx+1) )
                      ddsub_calc1(vec1(j),vec1_lo(j),sum_lx(jx+1),sum_lx_lo(jx+1),c_hi,c_lo)
                      ddmul_calc1(diag(j),diag_lo(j),c_hi,c_lo,vec2(j),vec2_lo(j))
                   end do
                   cycle
                endif
                do jx=nxs_th, nxe_th -1
                   j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                   ! vec2(j) = diag(j) *( vec1(j) - sum_lx(jx+1) )
                   ddsub_calc1(vec1(j),vec1_lo(j),sum_lx(jx+1),sum_lx_lo(jx+1),c_hi,c_lo)
                   ddmul_calc1(diag(j),0.0d0,c_hi,c_lo,vec2(j),vec2_lo(j))
                end do
             else
                do jx=nxs_th, nxe_th -1
                   j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                   vec2(j) = diag(j) *( vec1(j) - sum_lx(jx+1) )
                end do
             endif
          end do
       enddo
    else
       do jz=nzs_th , nze_th -1
          do jy=nys_th, nye_th -1
             do jx=nxs_th, nxe_th -1
                j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                sum_lx0 = 0.0d0
                if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                   sum_lx0_lo = 0.0d0
                endif
                do ipL=1,ipLn
                   ip = ip_diaL(ipL)
                   offset = ioff_dia(ip)
                   if ( jz + ioff_z(ip) < nzs_th) cycle

                   if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                      if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                         ddmul_calc1(val_dia(mdist*(ip-1)+j),val_dia_lo(mdist*(ip-1)+j),vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo)
                         ddadd_calc1(sum_lx0,sum_lx0_lo,c_hi,c_lo,sum_lx0,sum_lx0_lo)
                         cycle
                      endif
                      ddmul_calc1(val_dia(mdist*(ip-1)+j),0.0d0,vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo)
                      ddadd_calc1(sum_lx0,sum_lx0_lo,c_hi,c_lo,sum_lx0,sum_lx0_lo)
                   else
                      sum_lx0 = sum_lx0 +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                   endif
                end do
                if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                   if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                      ddsub_calc1(vec1(j),vec1_lo(j),sum_lx0,sum_lx0_lo,c_hi,c_lo)
                      ddmul_calc1(diag(j),diag_lo(j),c_hi,c_lo,vec2(j),vec2_lo(j))
                      cycle
                   endif
                   ddsub_calc1(vec1(j),vec1_lo(j),sum_lx0,sum_lx0_lo,c_hi,c_lo)
                   ddmul_calc1(diag(j),0.0d0,c_hi,c_lo,vec2(j),vec2_lo(j))
                else
                   vec2(j) = diag(j) *( vec1(j) -sum_lx0 )
                endif

             end do
          end do
       enddo
    endif

    !$OMP end parallel

    !$OMP Parallel private( num_threads, my_thread ) &
    !$OMP private( i, jz, jy, jx, j, sum_lx, sum_lx0,offset, ip , ipU) &
    !$OMP private( nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th ) &
    !$OMP private( &
    !$OMP ddmul_var1, &
    !$OMP ddadd_var1, &
    !$OMP ddsub_var1, &
    !$OMP sum_lx_lo,sum_lx0_lo, c_hi,c_lo )

    num_threads = omp_get_num_threads();
    my_thread = omp_get_thread_num();

    call set_region_for_threads( &
         my_thread, num_threads, nx, ny, nz, &
         nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th, &
         div_direc_th )

    if(prm%iflag_DDM_independ.eq.KSP_PRECONDITION_DDM_ILUindepend)then
       do jz=nze_th -1, nzs_th, -1
          do jy=nye_th -1, nys_th, -1
             do jx=nxe_th -1, nxs_th, -1
                sum_lx(jx+1) = 0.0d0
             enddo

             if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                do jx=nxs_th, nxe_th -1
                   sum_lx_lo(jx+1) = 0.0d0
                enddo
             endif
             do ipU=1,ipUn
                ip = ip_diaU(ipU)
                offset = ioff_dia(ip)
                if( (ioff_y(ip).eq.0).and.(ioff_z(ip).eq.0) )then
                   cycle
                endif
                if ( jz +ioff_z(ip) > nze_th -1 ) then
                   cycle
                endif
                if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                   if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                      do jx=nxe_th -1, nxs_th, -1
                         j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                         ! sum_lx(jx+1) = sum_lx(jx+1) +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                         ddmul_calc1( val_dia(mdist*(ip-1)+j),val_dia_lo(mdist*(ip-1)+j),vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo )
                         ddadd_calc1( sum_lx(jx+1),sum_lx_lo(jx+1),c_hi,c_lo,sum_lx(jx+1),sum_lx_lo(jx+1))
                      end do
                      cycle
                   endif
                   do jx=nxe_th -1, nxs_th, -1
                      j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                      ! sum_lx(jx+1) = sum_lx(jx+1) +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                      ddmul_calc1( val_dia(mdist*(ip-1)+j),0.0d0,vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo )
                      ddadd_calc1( sum_lx(jx+1),sum_lx_lo(jx+1),c_hi,c_lo,sum_lx(jx+1),sum_lx_lo(jx+1))
                   end do
                else
                   do jx=nxe_th -1, nxs_th, -1
                      j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                      sum_lx(jx+1) = sum_lx(jx+1) +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                   end do
                endif
             enddo
             if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                   do jx=nxe_th -1, nxs_th, -1
                      j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                      ! vec2(j) = vec2(j) -diag(j)*sum_lx(jx+1)
                      ddmul_calc1(diag(j),diag_lo(j),sum_lx(jx+1),sum_lx_lo(jx+1),c_hi,c_lo)
                      ddsub_calc1(vec2(j),vec2_lo(j),c_hi,c_lo,vec2(j),vec2_lo(j))
                   end do
                   cycle
                endif
                do jx=nxe_th -1, nxs_th, -1
                   j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                   ! vec2(j) = vec2(j) -diag(j)*sum_lx(jx+1)
                   ddmul_calc1(diag(j),0.0d0,sum_lx(jx+1),sum_lx_lo(jx+1),c_hi,c_lo)
                   ddsub_calc1(vec2(j),vec2_lo(j),c_hi,c_lo,vec2(j),vec2_lo(j))
                end do
             else
                do jx=nxe_th -1, nxs_th, -1
                   j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1;
                   vec2(j) = vec2(j) -diag(j)*sum_lx(jx+1)
                end do
             endif
          end do
       end do
    else
       do jz=nze_th -1, nzs_th, -1
          do jy=nye_th -1, nys_th, -1
             do jx=nxe_th -1, nxs_th, -1
                j = ( jx +stmx ) +mx *( jy +stmy ) +mxy *( jz +stmz ) +1
                sum_lx0 = 0.0d0
                if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                   sum_lx0_lo = 0.0d0
                endif
                do ipU=1,ipUn
                   ip = ip_diaU(ipU)
                   offset = ioff_dia(ip)
                   if ( jz +ioff_z(ip) > nze_th -1 ) cycle
                   if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                      if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                         ddmul_calc1(val_dia(mdist*(ip-1)+j),val_dia_lo(mdist*(ip-1)+j),vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo)
                         ddadd_calc1(sum_lx0,sum_lx0_lo,c_hi,c_lo,sum_lx0,sum_lx0_lo)
                         cycle 
                      endif
                      ddmul_calc1(val_dia(mdist*(ip-1)+j),0.0d0,vec2(j+offset),vec2_lo(j+offset),c_hi,c_lo)
                      ddadd_calc1(sum_lx0,sum_lx0_lo,c_hi,c_lo,sum_lx0,sum_lx0_lo)
                   else
                      sum_lx0 = sum_lx0 +val_dia( mdist*(ip-1) +j ) *vec2(j+offset)
                   endif
                end do
                if( PRECISION.eq.KSP_DOUBLE_DOUBLE_PRECISION )then
                   if(PRECISION_PRECONDITION.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
                      ! vec2(j) = vec2(j) -diag(j)*sum_lx0
                      ddmul_calc1(diag(j),diag_lo(j),sum_lx0,sum_lx0_lo,c_hi,c_lo)
                      ddsub_calc1(vec2(j),vec2_lo(j),c_hi,c_lo,vec2(j),vec2_lo(j))
                      cycle
                   endif
                   ! vec2(j) = vec2(j) -diag(j)*sum_lx0
                   ddmul_calc1(diag(j),0.0d00,sum_lx0,sum_lx0_lo,c_hi,c_lo)
                   ddsub_calc1(vec2(j),vec2_lo(j),c_hi,c_lo,vec2(j),vec2_lo(j))
                else
                   vec2(j) = vec2(j) -diag(j)*sum_lx0
                endif
             end do
          end do
       end do
    endif

    !$OMP end Parallel

    deallocate( ioff_x );  deallocate( ioff_y );  deallocate( ioff_z )

    call omp_set_num_threads(max_numth_0)


  end subroutine ilu_ddm_case_zdiv

end module as_ddm
