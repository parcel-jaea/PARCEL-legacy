module PARCEL_time
  implicit none
  real*8 precon_ddm_time(30)
  real*8 precon_time(30)
end module PARCEL_time

module KSP_STATUS
  use PARCEL_define
  implicit none

  type STATUS_FLAG
     integer machineflag
     integer dataflag
     integer calcflag
  end type STATUS_FLAG

contains
  subroutine init_staus_flag(flag)
    implicit none
    type(STATUS_FLAG) :: flag
    flag%machineflag = PARCEL_CPU_FLAG_ON
    flag%dataflag    = KSP_DOUBLE_PRECISION
    flag%calcflag    = KSP_DOUBLE_PRECISION
  end subroutine init_staus_flag

  subroutine set_machineflag(flag,machineflag)
    implicit none
    type(STATUS_FLAG) :: flag
    integer machineflag
    flag%machineflag = machineflag
  end subroutine set_machineflag

  subroutine set_dataflag(flag,dataflag)
    implicit none
    type(STATUS_FLAG) :: flag
    integer dataflag
    flag%dataflag = dataflag
  end subroutine set_dataflag

  subroutine set_calcflag(flag,calcflag)
    implicit none
    type(STATUS_FLAG) :: flag
    integer calcflag
    flag%calcflag = calcflag
  end subroutine set_calcflag

end module KSP_STATUS

module vec_type
  
  use c_wrapper
  use KSP_STATUS
  use PARCEL_define
  use PARCEL_integer
  implicit none

  type PARCEL_VEC
     real*8,allocatable :: vec_FP64(:)
     ! double - double 演算
     real*8,allocatable :: vec_hi(:)
     real*8,allocatable :: vec_lo(:)

     type(devicemem)    :: C_vec_FP64
     type(magma_memory) :: C_magma_vec_FP64     
     integer(kind=intType) n
     integer CPUflag
     integer GPUflag
     integer dataflag
     integer(kind=intType) offset
  end type PARCEL_VEC
  
end module vec_type

module common_kernel
contains
  subroutine get_maxthreads(max_numth)
    use omp_lib
    implicit none
    integer max_numth

    !$OMP PARALLEL
    !$OMP MASTER
    max_numth = omp_get_num_threads()
    !$OMP END MASTER
    !$OMP END PARALLEL

  end subroutine get_maxthreads

  subroutine set_precon_maxthreads(n,precon_thblock)
    use omp_lib
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer precon_thblock
    integer n_int4
    integer max_numth
    if(precon_thblock<1)then
       call get_maxthreads(max_numth)
       n_int4=n
       if(n_int4<=max_numth)then
          call omp_set_num_threads(n_int4)
       endif
    else
       call get_maxthreads(max_numth)
       if(max_numth<precon_thblock)then
          n_int4=n
          if(n_int4<=max_numth)then
             call omp_set_num_threads(n_int4)
!             call omp_set_num_threads(max_numth)
          else
             call omp_set_num_threads(max_numth)
          endif
       else
          n_int4=n
          if(n_int4<=precon_thblock)then
             call omp_set_num_threads(n_int4)
!             call omp_set_num_threads(precon_thblock)
          else
             call omp_set_num_threads(precon_thblock)
          endif
       endif
    endif
  end subroutine set_precon_maxthreads

  subroutine  divide_region( &
       imyrank_,inprocess_,idown_,iup_,idown_my_,iup_my_)
    use PARCEL_integer
    implicit none
    ! ランク
    integer imyrank_
    ! 分割数
    integer inprocess_
    ! 領域の下限
    integer(kind=intType) idown_
    ! 領域の上限
    integer(kind=intType) iup_
    ! ランクに割り振られた領域の下限
    integer(kind=intType) idown_my_
    ! ランクに割り振られた領域の上限
    integer(kind=intType) iup_my_

    integer imyrank,inprocess
    integer(kind=intType) idown,iup,idown_my,iup_my
    integer(kind=intType) iregion,iblock,over_region
    imyrank   = imyrank_
    inprocess = inprocess_
    idown     = idown_
    iup       = iup_
    idown_my  = 99
    iup_my    = -99
  
    if(iup==idown)then
       idown_my_=idown
       iup_my_=iup
       return 
    endif

    iregion=iup-idown+1
    iblock=(iregion)/inprocess
    over_region=mod(iregion,inprocess)
    if(iblock==0)then
       write(*,*)"error process is over much \n"
       return 
    endif

    if(imyrank < over_region)then
       idown_my=idown+iblock*imyrank+imyrank
       iup_my=idown_my+iblock
    else
       idown_my=idown+iblock*imyrank+over_region
       iup_my=idown_my+iblock-1
    endif
  
    if((iup_my-idown_my+1)<1)then
       write(*,*)"error process is over much \n"
       return 
    endif

    idown_my_=idown_my
    iup_my_=iup_my
    return 
  end subroutine divide_region

  subroutine set_Head_Tail(downmax,upmax,mydown,myup,myHead,myTail)
    ! //      upmax   :: 計算領域の上限
    ! //      downmax :: 計算領域の下限
    ! //      myup    :: 分割領域の上限
    ! //      mydown  :: 分割領域の下限
    ! //      myHead  :: 計算範囲の下限(出力)
    ! //      myTail  :: 計算範囲の上限(出力)
    use PARCEL_integer
    implicit none
    integer(kind=intType) upmax,downmax,myup,mydown,myHead,myTail

    if(mydown>downmax)then
       myHead=mydown
    else
       myHead=downmax
    endif
    if(myup<upmax)then
       myTail=myup  
    else
       myTail=upmax
    endif

    !  // 領域外
    if(mydown>upmax)then
       myHead=666
       myTail=-666
    endif
    if(myup<downmax)then
       myHead=666
       myTail=-666
    endif
  end subroutine set_Head_Tail

  ! Floor of binary logarithm
  ! If argument n is less than 1 then assign dummy(-9999) to n.
  ! If n is 2**i and i is integer(kind=intType) then assign exactly i to n.
  ! Otherwise, assign floor(log_2 (n)) to n.
  subroutine floor_log2(n)
    use PARCEL_integer
    implicit none
    integer,intent(inout) :: n
    integer(kind=intType) :: power_of_2, power

    if(n.le.0)then
       ! If n<=0, no log2(n) exists. So, assign n to dummy(-9999).
       n=-9999
       return
    endif

#ifdef PARCEL_INT4
    power_of_2=2
    do power=0,29
       if(n.lt.power_of_2)then
          n=power
          return
       endif
       power_of_2=power_of_2*2
    enddo

    ! 4 Byte integer type variable is definitely less than 2**31
    n=30
#else
    power_of_2=2
    do power=0,61
       if(n.lt.power_of_2)then
          n=power
          return
       endif
       power_of_2=power_of_2*2
    enddo
    ! 8 Byte integer(kind=intType) type variable is definitely less than 2**31
    n=62
#endif

  end subroutine floor_log2

  subroutine dlog2(x)
    implicit none
    real*8 x
    x=dlog(x)/dlog(2.0d0)
  end subroutine dlog2

  subroutine set_region_for_threads( &
       my_thread, num_threads, nx, ny, nz, &
       nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th, div_dir )
    ! use ddm_tmp
!    use ddmMat_wrap
    implicit none
    ! inout
    integer, intent(in) :: my_thread, num_threads
    integer, intent(in) :: nx, ny, nz, div_dir
    !
    integer, intent(inout) :: nxs_th, nys_th, nzs_th
    integer, intent(inout) :: nxe_th, nye_th, nze_th
    ! local
    integer n1
    !
    nxs_th =  0;  nys_th =  0;  nzs_th =  0;
    nxe_th = nx;  nye_th = ny;  nze_th = nz;

    if ( div_dir == 1 ) then
       n1 = nx /num_threads;
       nxs_th = n1 *my_thread;     nxe_th = nxs_th +n1;
       if ( my_thread == num_threads -1 ) nxe_th = nx
    else if ( div_dir == 2 ) then
       n1 = ny /num_threads;
       nys_th = n1 *my_thread;     nye_th = nys_th +n1;
       if ( my_thread == num_threads -1 ) nye_th = ny
    else if ( div_dir == 3 ) then
       n1 = nz /num_threads;
       nzs_th = n1 *my_thread;     nze_th = nzs_th +n1;
       if ( my_thread == num_threads -1 ) nze_th = nz
    endif
    if ( div_dir == 0 ) then
       if ( my_thread > 0 ) then
          nxe_th = 0;    nye_th = 0;  nze_th = 0 
       endif
    endif
    return
  end subroutine set_region_for_threads

end module common_kernel

! ----
! C(m,n)=A(m,k)*B(k,n)
module wrap_DGEMM

contains
  subroutine wrap_BASIC_DGEMM(m,n,k,A,B,C)
    use PARCEL_integer
    implicit none
    integer(kind=intType) m,n,k
    real*8 A(m*k)
    real*8 B(k*n)
    real*8 C(m*n)
    real*8 ALPHA,BETA
    integer m_int4,n_int4,k_int4
    ALPHA=1.0d0
    BETA =0.0d0
    m_int4=m
    n_int4=n
    k_int4=k
    call DGEMM('N','N',m_int4,n_int4,k_int4,ALPHA,A,m_int4,B,k_int4,BETA,C,m_int4)
  end subroutine wrap_BASIC_DGEMM
end module wrap_DGEMM
