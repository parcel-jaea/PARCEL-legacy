! Todo 通信が上位ビット、下位ビット個別に送付するようになっている 
!          --> 一度の通信で完了するように修正
!          --> Allreduce が1プロセスに集めて計算を行う実装になっているためAllreduceの実行で計算できるように修正する
module ddcomm_kernel

contains

  subroutine ddget_nrk_vec( &
       commVec,vec_hi,vec_lo,rbuf_hi,rbuf_lo&
       )
    use comm
    implicit none
    type(comm_Vec) :: commVec
    real*8 vec_hi(commVec%n)
    real*8 vec_lo(commVec%n)
    real*8 rbuf_hi(commVec%nrk_data)
    real*8 rbuf_lo(commVec%nrk_data)
    real*8 sbuf_hi(commVec%nrk_data_send)
    real*8 sbuf_lo(commVec%nrk_data_send)

    call ddnrk_vec_pack( &
         commVec%nrk_prm_send, &
         commVec%n, &
         commVec%nrk_data_send, &
         vec_hi,vec_lo, &
         sbuf_hi,sbuf_lo, &
         commVec%nrk_rank_send)
    call ddnrk_vec_comm(commVec,&
         sbuf_hi,sbuf_lo,rbuf_hi,rbuf_lo)

  end subroutine ddget_nrk_vec

  subroutine ddnrk_vec_pack( &
       nrk_prm,n,ndata,&
       vec_hi,vec_lo,&
       sbuf_hi,sbuf_lo,&
       nrk_rank)
    use comm
    use PARCEL_integer
    implicit none
    type(nrk) :: nrk_prm(nrk_rank)
    integer(kind=intType) n
    integer ndata
    real*8 vec_hi(n)
    real*8 vec_lo(n)
    real*8 sbuf_hi(ndata)
    real*8 sbuf_lo(ndata)
    integer nrk_rank
    integer(kind=intType) offset
    integer(kind=intType) i,j,k,i0

    !$OMP PARALLEL PRIVATE(j,i0,k,i,offset) 
    do j=1,nrk_rank
       i0=0
       do k=1,j-1
          i0=i0+nrk_prm(k)%ndata
       enddo
       !$OMP DO
       do i=1,nrk_prm(j)%ndata
          offset=nrk_prm(j)%offset(i)
          sbuf_hi(i+i0)=vec_hi(1+offset)
          sbuf_lo(i+i0)=vec_lo(1+offset)
       enddo
       !$OMP END DO NOWAIT

    enddo
    !$OMP END PARALLEL 

  end subroutine ddnrk_vec_pack

  subroutine ddnrk_vec_comm(commVec ,&
         sbuf_hi,sbuf_lo,rbuf_hi,rbuf_lo)
    use mpi
    use comm
    use PARCEL_integer
    implicit none
    type(comm_Vec) :: commVec

    integer(kind=intType) i,j,k
    integer(kind=intType) i0
    integer ierr

    real*8 sbuf_hi(commVec%nrk_data_send)
    real*8 sbuf_lo(commVec%nrk_data_send)
    real*8 rbuf_hi(commVec%nrk_data)
    real*8 rbuf_lo(commVec%nrk_data)

    integer s_req(commVec%nrk_rank_send)
    integer r_req(commVec%nrk_rank)

    type send_stat
       integer stat(MPI_STATUS_SIZE)
    end type send_stat
    type recv_stat
       integer stat(MPI_STATUS_SIZE)
    end type recv_stat
    type(send_stat) :: s_sta(commVec%nrk_rank_send)
    type(recv_stat) :: r_sta(commVec%nrk_rank)

!    integer s_sta(MPI_STATUS_SIZE,commVec%nrk_rank_send)
!    integer r_sta(MPI_STATUS_SIZE,commVec%nrk_rank)

    do i=1,commVec%nrk_rank_send
       i0=0
       do k=1,i-1
          i0=i0+commVec%nrk_prm_send(k)%ndata
       enddo
       call MPI_Isend( &
            sbuf_hi(1+i0),&
            commVec%nrk_prm_send(i)%ndata ,&
            MPI_REAL8,&
            commVec%nrk_prm_send(i)%rank,&
            commVec%myrank,&
            commVec%comm,&
            s_req(i),&
            ierr)
    enddo
    do i=1,commVec%nrk_rank
       i0=0
       do k=1,i-1
          i0=i0+commVec%nrk_prm(k)%ndata
       enddo
       call MPI_Irecv( &
            rbuf_hi(1+i0),&
            commVec%nrk_prm(i)%ndata ,&
            MPI_REAL8, &
            commVec%nrk_prm(i)%rank, &
            commVec%nrk_prm(i)%rank, &
            commVec%comm, &
            r_req(i), &
            ierr)
    enddo
    do i=1,commVec%nrk_rank_send
       call MPI_Wait(s_req(i),s_sta(i)%stat,ierr)
    enddo
    do i=1,commVec%nrk_rank
       call MPI_Wait(r_req(i),r_sta(i)%stat,ierr)
    enddo

! -----

    do i=1,commVec%nrk_rank_send
       i0=0
       do k=1,i-1
          i0=i0+commVec%nrk_prm_send(k)%ndata
       enddo
       call MPI_Isend( &
            sbuf_lo(1+i0),&
            commVec%nrk_prm_send(i)%ndata ,&
            MPI_REAL8,&
            commVec%nrk_prm_send(i)%rank,&
            commVec%myrank,&
            commVec%comm,&
            s_req(i),&
            ierr)
    enddo
    do i=1,commVec%nrk_rank
       i0=0
       do k=1,i-1
          i0=i0+commVec%nrk_prm(k)%ndata
       enddo
       call MPI_Irecv( &
            rbuf_lo(1+i0),&
            commVec%nrk_prm(i)%ndata ,&
            MPI_REAL8, &
            commVec%nrk_prm(i)%rank, &
            commVec%nrk_prm(i)%rank, &
            commVec%comm, &
            r_req(i), &
            ierr)
    enddo
    do i=1,commVec%nrk_rank_send
       call MPI_Wait(s_req(i),s_sta(i)%stat,ierr)
    enddo
    do i=1,commVec%nrk_rank
       call MPI_Wait(r_req(i),r_sta(i)%stat,ierr)
    enddo


  end subroutine ddnrk_vec_comm

  subroutine ddget_nrk_vecLU( &
       commVec,&
       rbuff_hi,rbuff_lo,&
       nColL,nColU,&
       rbuffL_hi,rbuffL_lo,&
       rbuffU_hi,rbuffU_lo &
       )
    use comm
    use PARCEL_integer
    implicit none
    integer(kind=intType) nColL,nColU
    type(comm_Vec) :: commVec
    real*8 rbuff_hi(commVec%nrk_data)
    real*8 rbuff_lo(commVec%nrk_data)

    real*8 rbuffL_hi(nColL)
    real*8 rbuffL_lo(nColL)
    real*8 rbuffU_hi(nColU)
    real*8 rbuffU_lo(nColU)
    integer(kind=intType) i

    !$OMP PARALLEL DO PRIVATE(i)  
    do i=1,nColL
       rbuffL_hi(i)=rbuff_hi(i)
       rbuffL_lo(i)=rbuff_lo(i)
    enddo
    !$OMP END PARALLEL DO 
    !$OMP PARALLEL DO PRIVATE(i)  
    do i=1,nColU
       rbuffU_hi(i)=rbuff_hi(i+nColL)
       rbuffU_lo(i)=rbuff_lo(i+nColL)
    enddo
    !$OMP END PARALLEL DO 
  end subroutine ddget_nrk_vecLU
  
  subroutine ddnrk_vec_pack_ovl2( &
       nrk_prm,&
       n,&
       ndata,&
       vec_hi,vec_lo,&
       sbuf_hi,sbuf_lo,&
       nrk_rank,&
       srank)
    use comm
    use PARCEL_integer
    implicit none
    type(nrk) :: nrk_prm(nrk_rank)
    integer(kind=intType) n
    integer ndata
    real*8 vec_hi(n)
    real*8 vec_lo(n)
    real*8 sbuf_hi(ndata)
    real*8 sbuf_lo(ndata)
    integer nrk_rank
    integer srank

    integer(kind=intType) offset
    integer(kind=intType) i,j,k,i0
    
    if(srank>-1)then
       j=srank
       !$OMP PARALLEL DO PRIVATE(i,offset) 
       do i=1,nrk_prm(j)%ndata
          offset=nrk_prm(j)%offset(i)
          sbuf_hi(i)=vec_hi(1+offset)
          sbuf_lo(i)=vec_lo(1+offset)
       enddo
       !$OMP END PARALLEL DO
    endif
  end subroutine ddnrk_vec_pack_ovl2
  
  subroutine ddnrk_vec_comm_ovl2(commVec,&
       sbuf_hi,sbuf_lo,&
       rbuf_hi,rbuf_lo,&
       srank,rrank)
    use comm
    use PARCEL_integer
    implicit none
    type(comm_Vec) :: commVec

    integer(kind=intType) i,j,k
    integer(kind=intType) i0
    integer ierr

    real*8 sbuf_hi(commVec%nrk_data_max_send)
    real*8 sbuf_lo(commVec%nrk_data_max_send)
    real*8 rbuf_hi(commVec%nrk_data_max)
    real*8 rbuf_lo(commVec%nrk_data_max)
    integer srank,rrank

    integer s_req_lo
    integer s_sta_lo(MPI_STATUS_SIZE)
    integer r_req_lo
    integer r_sta_lo(MPI_STATUS_SIZE)

    integer s_req_hi
    integer s_sta_hi(MPI_STATUS_SIZE)
    integer r_req_hi
    integer r_sta_hi(MPI_STATUS_SIZE)

    if(srank>-1)then
       i=srank
       call MPI_Isend( &
            sbuf_hi,&
            commVec%nrk_prm_send(i)%ndata ,&
            MPI_REAL8,&
            commVec%nrk_prm_send(i)%rank,&
            commVec%myrank,&
            commVec%comm,&
            s_req_hi,&
            ierr)
    endif
    if(rrank>-1)then
       i=rrank
       call MPI_Irecv( &
            rbuf_hi,&
            commVec%nrk_prm(i)%ndata ,&
            MPI_REAL8, &
            commVec%nrk_prm(i)%rank, &
            commVec%nrk_prm(i)%rank, &
            commVec%comm, &
            r_req_hi, &
            ierr)
    endif

    if(srank>-1)then
       call MPI_Wait(s_req_hi,s_sta_hi,ierr)
    endif
    if(rrank>-1)then
       call MPI_Wait(r_req_hi,r_sta_hi,ierr)
    endif
! -----

    if(srank>-1)then
       i=srank
       call MPI_Isend( &
            sbuf_lo,&
            commVec%nrk_prm_send(i)%ndata ,&
            MPI_REAL8,&
            commVec%nrk_prm_send(i)%rank,&
            commVec%myrank,&
            commVec%comm,&
            s_req_lo,&
            ierr)
    endif
    if(rrank>-1)then
       i=rrank
       call MPI_Irecv( &
            rbuf_lo,&
            commVec%nrk_prm(i)%ndata ,&
            MPI_REAL8, &
            commVec%nrk_prm(i)%rank, &
            commVec%nrk_prm(i)%rank, &
            commVec%comm, &
            r_req_lo, &
            ierr)
    endif

    if(srank>-1)then
       call MPI_Wait(s_req_lo,s_sta_lo,ierr)
    endif
    if(rrank>-1)then
       call MPI_Wait(r_req_lo,r_sta_lo,ierr)
    endif
! -----

  end subroutine ddnrk_vec_comm_ovl2


end module ddcomm_kernel
