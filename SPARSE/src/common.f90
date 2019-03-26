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

  subroutine set_precon_maxthreads(n)
    use omp_lib
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer n_int4
    integer max_numth
    call get_maxthreads(max_numth)
    n_int4=n
    if(n_int4<=max_numth)then
!    if(n<=max_numth)then
!       call omp_set_num_threads(n)
       call omp_set_num_threads(n_int4)
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

end module common_kernel

