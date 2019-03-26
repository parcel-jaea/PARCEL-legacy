module sort_kernel

contains
  subroutine arrayMerge( &
       nprocess, &
       ndata, &
       myid, &
       pairid, &
       nthblock_omp, &
       nthstart_omp, &
       val,mval, &
       jstep)
    use PARCEL_integer
    implicit none

    integer nprocess
    integer(kind=intType) ndata
    integer myid,pairid
    integer(kind=intType) nthblock_omp(nprocess)
    integer(kind=intType) nthstart_omp(nprocess)
    integer(kind=intType) val(ndata)
    integer(kind=intType) mval(ndata)
    integer(kind=intType) jstep

    integer(kind=intType) i,k
    integer(kind=intType) nthblock 
    integer(kind=intType) myvalid  
    integer(kind=intType) pairvalid 
    integer(kind=intType) myidmax,pairidmax
    integer(kind=intType) sortid
    integer(kind=intType) myval,pairval

    sortid    = nthstart_omp(myid+1) 
    nthblock  = nthblock_omp(myid+1)+nthblock_omp(pairid+1)
    myvalid   = nthstart_omp(myid+1)
    pairvalid = nthstart_omp(pairid+1)

    myidmax   = nthstart_omp(myid+1) + nthblock_omp(myid+1)   -1
    pairidmax = nthstart_omp(pairid+1) + nthblock_omp(pairid+1) -1

    ! write(*,*)jstep,"myid    ,pairid    = ",myid,pairid
    ! write(*,*)jstep,"myvalid ,pairvalid = ",myvalid,pairvalid
    ! write(*,*)jstep,"myidmax ,pairidmax = ",myidmax,pairidmax
    do i=1,nthblock
       myval  =val(  myvalid)
       pairval=val(pairvalid)

       if(myval<pairval)then
          mval(sortid)=myval
          myvalid=myvalid+1
          sortid=sortid+1
       else
          mval(sortid)=pairval
          pairvalid=pairvalid+1
          sortid=sortid+1
       endif
       if(myvalid.eq.(myidmax+1))then
          do k=pairvalid,pairidmax
             pairval=val(k)
             !           write(*,*)"  k,sortid,pairval ",k,sortid,pairval
             mval(sortid)=pairval
             sortid=sortid+1
          enddo
          exit
       endif
       if(pairvalid.eq.(pairidmax+1))then
          do k=myvalid,myidmax
             myval=val(k)
             !           write(*,*)"  k,sortid,myval ",k,sortid,myval
             mval(sortid)=myval
             sortid=sortid+1
          enddo
          exit
       endif
    enddo

    do i=nthstart_omp(myid+1),nthstart_omp(myid+1)+nthblock-1
       val(i)=mval(i)
    enddo
    nthblock_omp(myid+1) = nthblock

  end subroutine ArrayMerge

  subroutine binaryTreeMerge( &
       max_nprocess,&
       ndata, &
       nthblock_omp, &
       nthstart_omp, &
       val)
    use omp_lib
    use common_kernel
    use PARCEL_integer
    implicit none

    integer max_nprocess
    integer(kind=intType) ndata
    integer(kind=intType) nthblock_omp(max_nprocess)
    integer(kind=intType) nthstart_omp(max_nprocess)
    integer(kind=intType) val(ndata)
    integer(kind=intType) mval(ndata)
    integer(kind=intType) nthBarrier(max_nprocess)

    integer nprocess,irank
    integer(kind=intType) i,j,k
    integer(kind=intType) i_,j_
    integer(kind=intType) stage
    integer nrk
    real*8  tmp
    integer(kind=intType) itmp
    integer istage
    integer irank_stage,irank0,irank0_master

    !$OMP PARALLEL PRIVATE(i)
    !$OMP MASTER
    nprocess= omp_get_num_threads()
    !$OMP END MASTER

    !$OMP DO
    do i=1,ndata
       mval(i)=-1
    end do

    !$OMP END PARALLEL
    nthBarrier=-1

    tmp=nprocess
#if 0
    call dlog2(tmp)
    istage=tmp
    itmp=1
    do i=1,istage+1
       itmp=itmp*2
    enddo
    if(nprocess.eq.itmp)then
       istage=istage+1
    endif
#else
    istage=tmp
    call floor_log2(istage)
#endif

    ! -----
    !$OMP PARALLEL PRIVATE( &
    !$OMP nprocess, &
    !$OMP irank, &
    !$OMP irank0, &
    !$OMP irank0_master, &
    !$OMP irank_stage, &
    !$OMP nrk, &
    !$OMP stage, &
    !$OMP i,j &
    !$OMP )

    stage =0
    irank0=0
    irank0_master=0
    nprocess = omp_get_num_threads();
    irank    = omp_get_thread_num();

    irank_stage=irank
    do i=istage,0,-1
       !     !$OMP BARRIER
       if(irank_stage<2**i)then
          ! ----
          ! irank0_master --> マージしたランク
          ! irank0_master - irank0 --> マージするランク
          ! -----
          j=0
          stage=1
          do while(stage<2**i)
             j=j+1
             stage=2**j
             !           if(mod(irank,stage).eq.0)then
             if(mod(irank,stage).eq.0)then
                if(irank.eq.(nprocess-1))then
                   nrk=irank0_master
                   stage=2**i
                   !! Pairとして処理される
                   nthBarrier(irank+1)=0
                   exit
                else
                   nrk=irank+stage/2
                   !! マージを行う
                   if(nrk<(nprocess))then
                      ! nrk の配列をマージ
                      !                    write(*,*)"stage myrank pair rank = ",stage,irank,nrk
                      ! call MPI_Irecv(rbuff,s*s,MPI_DOUBLE,nrk,irank,comm1,req_recv,ierr)
                      ! call MPI_Wait(req_recv, stat_recv,ierr)

                      do 
                         !$OMP FLUSH(nthBarrier) 
                         if(nthBarrier(nrk+1).eq.0)then
                            exit
                         endif
                      enddo
                      !$OMP FLUSH(nthblock_omp,val,mval) 
                      call arrayMerge( &
                           nprocess, &
                           ndata, &
                           irank, &
                           nrk, &
                           nthblock_omp, &
                           nthstart_omp, &
                           val,mval, &
                           stage)

                   endif
                endif
             else
                nrk=irank-stage/2
                !! Pairとして処理される
                nthBarrier(irank+1)=0
                exit
             endif
          end do
          ! -----
          ! ここで1ステージ上のツリーとマージする

          if((irank0_master.eq.irank).and.(irank.ne.0))then
             if((irank+2**i)<nprocess)then
                nrk=irank+2**i
                ! nrk の配列をマージ
                ! call MPI_Irecv(rbuff, s*s, MPI_DOUBLE, nrk, irank,   comm1, req_recv,ierr)
                ! call MPI_Wait(req_recv, stat_recv,ierr)            
                do 
                   !$OMP FLUSH(nthBarrier) 
                   if(nthBarrier(nrk+1).eq.0)then
                      exit
                   endif
                enddo
                !$OMP FLUSH(nthblock_omp,val,mval) 
                call arrayMerge( &
                     nprocess, &
                     ndata, &
                     irank, &
                     nrk, &
                     nthblock_omp, &
                     nthstart_omp, &
                     val,mval, &
                     stage)
                !              write(*,*)"stage myrank pair rank = ",stage,irank,nrk

             endif
             nrk=irank0_master - irank0
             nthBarrier(irank+1)=0

             !! Pairとして処理される
             !write(*,*)"Send stage = ",stage,irank,nrk
          endif
          if((irank.eq.0).and.(2**(i)<nprocess))then
             nrk=irank+2**(i)
             ! nrk の配列をマージ
             !           write(*,*)"stage myrank pair rank = ",stage,irank,nrk
             ! call MPI_Irecv(rbuff, s*s, MPI_DOUBLE, nrk, irank,   comm1, req_recv,ierr)
             ! call MPI_Wait(req_recv, stat_recv,ierr)   

             do 
                !$OMP FLUSH(nthBarrier) 
                if(nthBarrier(nrk+1).eq.0)then
                   exit
                endif
             enddo
             !$OMP FLUSH(nthblock_omp,val,mval) 
             call arrayMerge( &
                  nprocess, &
                  ndata, &
                  irank, &
                  nrk, &
                  nthblock_omp, &
                  nthstart_omp, &
                  val,mval, &
                  stage)
          endif
          exit
       else
          irank0=2**(i)
          irank0_master=irank0_master+irank0
          irank_stage=mod(irank_stage,2**i)
       endif
    enddo
    nthBarrier(irank+1)=0

    !$OMP END PARALLEL

  end subroutine binaryTreeMerge

  subroutine MergeSort(val,n)
    use common_kernel
    use omp_lib
    use PARCEL_integer
    implicit none

    integer(kind=intType) n
    integer(kind=intType) val(n)
    integer(kind=intType),allocatable :: nthblock_omp(:)
    integer(kind=intType),allocatable :: nthstart_omp(:)

    integer max_numth
    integer numth,myth 
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) nthblock
    integer(kind=intType)  ntmp
    if(n<=1)then
       ! nが1以下なら並び替える必要がない
       return
    endif

    call get_maxthreads(max_numth)
    call set_precon_maxthreads(n)
    allocate(nthblock_omp(max_numth))
    allocate(nthstart_omp(max_numth))

    !@@    write(*,*)max_numth,n
    !$OMP PARALLEL &

    !$OMP PRIVATE( &
    !$OMP numth,myth,n_startth,n_endth, &
    !$OMP nthblock,ntmp &
    !$OMP )
    numth = omp_get_num_threads();
    myth  = omp_get_thread_num();

    ntmp=1
    call divide_region( &
         myth,numth,ntmp,n,n_startth,n_endth)

    nthblock = n_endth - n_startth+1
    nthblock_omp(myth+1)=nthblock
    nthstart_omp(myth+1)=n_startth
!    call sort_c(val(n_startth),nthblock)
    call quicksort(nthblock,val(n_startth))
    !$OMP END PARALLEL 

    call binaryTreeMerge( &
         max_numth,&
         n, &
         nthblock_omp, &
         nthstart_omp, &
         val)
    deallocate(nthblock_omp)
    deallocate(nthstart_omp)

    if(n<=max_numth)then
       call omp_set_num_threads(max_numth)
    endif

  end subroutine MergeSort

  recursive subroutine quicksort_main(n,x,l,r)
    use PARCEL_integer
    implicit none

    integer(kind=intType) n,l,r
    integer(kind=intType) x(n)
    integer(kind=intType) c, ls, i
    integer(kind=intType) tmp
    if (l < r) then
       c=(l+r)/2
       tmp=x(l)
       x(l)=x(c)
       x(c)=tmp
       ls=l
       do i=l+1, r
          if (x(i) < x(l)) then
             ls=ls+1
             tmp=x(ls)
             x(ls)=x(i)
             x(i)=tmp
          end if
       end do
       tmp=x(ls)
       x(ls)=x(l)
       x(l)=tmp
       call quicksort_main(n,x,l,ls-1)
       call quicksort_main(n,x,ls+1,r)
    end if
  end subroutine quicksort_main

  subroutine quicksort(n,x)
    use PARCEL_integer
    implicit none

    integer(kind=intType) n
    integer(kind=intType) x(n)
    call quicksort_main(n,x,1,n)
  end subroutine quicksort

end module sort_kernel

