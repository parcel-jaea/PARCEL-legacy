  subroutine MatVec_crs_ovl2(prm,crsMat,x,Ax,commVec)
    use matrix
    use comm
    use comm_kernel

    implicit none
    type(prmMat)  :: prm 
    type(CRS_Mat) :: crsMat 
    type(comm_Vec) :: commVec
    real*8 x(prm%n)

    real*8 Ax(prm%n)
    integer i,j,k,kk
    integer n,istart,iend
    integer Row,Col
    integer nColL,nColU
    real*8 sbuf(commVec%nrk_data_max_send)
    real*8 rbuf(commVec%nrk_data_max)
    real*8 rbuf0(commVec%nrk_data_max)

    integer irank,irankMax
    integer srank,srankMax
    integer rrank,rrankMax
    integer istart_rrank,iend_rrank  
    integer offset
    integer rrank0

    integer recv_rank0,send_rank0

    integer s_req
    integer s_sta(MPI_STATUS_SIZE)
    integer r_req
    integer r_sta(MPI_STATUS_SIZE)
    integer ierr
    integer isflag,irflag
    LOGICAL isendflag,irecvflag

    integer rrank_count,srank_count
    integer rrank_calc

    srankMax=commVec%nrk_rank_send
    rrankMax=commVec%nrk_rank

    n = prm%n
    istart=prm%istart
    iend=prm%iend

    ! -----
    send_rank0=1
    do irank=1,commVec%nrk_rank_send
       if(commVec%myrank < commVec%nrk_prm_send(irank)%rank)then
          send_rank0=irank
          exit
       endif
    enddo

    recv_rank0=commVec%nrk_rank
    do irank=commVec%nrk_rank,1,-1
       if(commVec%nrk_prm(irank)%rank < commVec%myrank   )then
          recv_rank0=irank
          exit
       endif
    enddo
!    write(*,*)"send_rank0,recv_rank0",send_rank0,recv_rank0

    rrank_calc=-1
    if(commVec%nrk_rank<1)then
       rrank =-1
       irflag=-1
       rrank_count=rrankMax+1
    else
       irflag=0
       rrank_count=1
#if 0
       rrank=1
#else
       rrank=recv_rank0-rrank_count+1
       if(rrank<=0)then
          rrank=rrankMax+rrank
       endif
#endif
    endif

    if(commVec%nrk_rank_send<1)then
       srank =-1
       isflag=-1
       srank_count=srankMax+1
    else
       isflag=0
       srank_count=1
#if 0
       srank=1
#else
       srank=mod(send_rank0+srank_count+1,srankMax)+1
#endif

    endif
    rrank0=0 ! 0は自分の領域
    ! -----

    ! isflag ( 全送信終了:-1,  送信完了:0, 送信中:1)
    ! irflag ( 全受信終了:-1,  受信完了:0, 受信中:1)
    do 

#if 0
       write(*,*)&
            commVec%myrank, &
            commVec%nrk_prm(rrank)%rank, &
            commVec%nrk_prm_send(srank)%rank
    call system("sleep 10")
    stop
#endif
       if( (0<srank).and.(isflag.eq.0) )then
          call nrk_vec_pack_ovl2( &
               commVec%nrk_prm_send, &
               commVec%n, &
               commVec%nrk_data_send, &
               x, &
               sbuf, &
               commVec%nrk_rank_send,&
               srank)
          call nrk_vec_comm_ovl2_isend( &
               commVec,sbuf,srank,s_req,s_sta)
          isflag=1
       endif
       if( (0<rrank).and.(irflag.eq.0) )then
          call nrk_vec_comm_ovl2_irecv(&
               commVec,rbuf,rrank,r_req,r_sta)
          irflag=1
       endif

       !$OMP PARALLEL PRIVATE( &
       !$OMP istart_rrank,iend_rrank ,j,i,k,kk,offset,Col) 

       !$OMP MASTER
       if((isflag.eq.1).or.(irflag.eq.1))then
          do 
             if(isflag.eq.1)then
                call MPI_TEST(s_req,isendflag,s_sta,ierr)
                if(isendflag)then
                   isflag=0
                   srank_count=srank_count+1
#if 1
                   srank=srank+1
#else
                   srank=mod(send_rank0+srank_count+1,srankMax)+1
#endif
                   exit
                endif
             endif
             if(irflag.eq.1)then
                call MPI_TEST(r_req,irecvflag,r_sta,ierr)
                if(irecvflag)then
                   irflag=0
                   rrank_count=rrank_count+1

                   rrank_calc=rrank
#if 1
                   rrank  = rrank+1
#else
                   rrank=recv_rank0-rrank_count+1
                   if(rrank<=0)then
                      rrank=rrankMax+rrank
                   endif
#endif
                   exit
                endif

             endif
          enddo
       endif
       !$OMP END MASTER

       if(rrank0.eq.0)then
          !$OMP DO SCHEDULE(dynamic,nBlock)
          do k=1,n
             Ax(k)=0.0d0
             do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
                if(crsMat%crsCol(j)<istart)then
                   cycle
                endif
                if( iend<crsMat%crsCol(j) )then
                   exit
                endif
                Col=crsMat%crsCol(j)-istart+1
                Ax(k)=Ax(k)+crsMat%crsA(j)*x(Col)
             enddo
          enddo
       endif
       if( 0<rrank0 )then
          !$OMP DO SCHEDULE(dynamic,nBlock)
          do k=1,n
             do kk = &
                  commVec%nrk_prm(rrank0)%nrk_ptr(k),&
                  commVec%nrk_prm(rrank0)%nrk_ptr(k+1)-1
                j=commVec%nrk_prm(rrank0)%nrkCol(kk)
                i=commVec%nrk_prm(rrank0)%nrkRow(kk)
                Ax(k)=Ax(k)+crsMat%crsA(j)*rbuf0(i)
             enddo
          enddo
       endif

       !$OMP END PARALLEL 

       if((isflag.eq.-1).and.(irflag.eq.-1))then
          if(rrank0.eq.-1)then
             exit
          endif
       endif

       if(irflag.eq.0)then
#if 0
          rrank0=rrank-1
#else
          rrank0=rrank_calc
#endif
          !$OMP PARALLEL DO PRIVATE(i)
          do i=1,commVec%nrk_data_max
             rbuf0(i)=rbuf(i)
          enddo
          !$OMP END PARALLEL DO
       else
          rrank0=-1
       endif

       if(srankMax<srank_count)then
          isflag=-1
          srank=-1
       endif
       if(rrankMax<rrank_count)then
          irflag=-1
          rrank=-1
       endif

    enddo
  end subroutine MatVec_crs_ovl2
