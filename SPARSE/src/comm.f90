
module comm_kernel
contains
  subroutine get_nrk_comm(commVec,nCol,Colsort)
    use comm
    use PARCEL_integer
    implicit none

    type(comm_Vec) :: commVec
    integer(kind=intType) nCol
    integer(kind=intType) Colsort(nCol)

    integer i,j,k
    integer nrank,rank
    integer istart_rank

    nrank=0
    istart_rank=1
    rank=-1
    do i=1,nCol
       do k=istart_rank,commVec%npes
          if(commVec%myrank.eq.(k-1))then
             cycle
          endif
          if( &
               (commVec%istart(k)<=Colsort(i)) .and. &
               (Colsort(i)<commVec%istart(k+1)) &
               )then
             if(rank.ne.(k-1))then
                rank=k-1
                nrank=nrank+1
                istart_rank=k
             endif
             exit
          endif
       enddo
    enddo
    commVec%nrk_rank=nrank
    allocate(commVec%nrk_prm(commVec%nrk_rank))
    do i=1,commVec%nrk_rank
       commVec%nrk_prm(i)%ndata=0
    enddo

    if(commVec%nrk_rank.ne.0)then
       nrank=0
       istart_rank=1
       rank=-1
       do i=1,nCol
          do k=istart_rank,commVec%npes
             if(commVec%myrank.eq.(k-1))then
                cycle
             endif
             if( &
                  (commVec%istart(k)<=Colsort(i)) .and. &
                  (Colsort(i)<commVec%istart(k+1)) &
                  )then
                if(rank.ne.(k-1))then
                   rank=k-1
                   nrank=nrank+1
                   istart_rank=k
                   commVec%nrk_prm(nrank)%rank=rank

                endif
                exit
             endif
          enddo
       enddo
       do i=1,nCol
          do j=1,commVec%nrk_rank
             k=commVec%nrk_prm(j)%rank+1
             if( &
                  (commVec%istart(k)<=Colsort(i)) .and. &
                  (Colsort(i)<commVec%istart(k+1)) &
                  )then
                commVec%nrk_prm(j)%ndata=commVec%nrk_prm(j)%ndata+1
                exit
             endif
          enddo
       enddo
       do i=1,commVec%nrk_rank
          allocate(  &
               commVec%nrk_prm(i)%offset(commVec%nrk_prm(i)%ndata) ) 
          commVec%nrk_prm(i)%ndata=0
       enddo

       do i=1,nCol
          do j=1,commVec%nrk_rank
             k=commVec%nrk_prm(j)%rank+1
             if( &
                  (commVec%istart(k)<=Colsort(i)) .and. &
                  (Colsort(i)<commVec%istart(k+1)) &
                  )then
                commVec%nrk_prm(j)%ndata=commVec%nrk_prm(j)%ndata+1
                commVec%nrk_prm(j)%offset(commVec%nrk_prm(j)%ndata) = &
                     Colsort(i)-commVec%istart(k)
                exit
             endif
          enddo
       enddo
    endif

    commVec%nrk_data=0
    commVec%nrk_data_max=0
    do i=1,commVec%nrk_rank
       commVec%nrk_data=commVec%nrk_data+commVec%nrk_prm(i)%ndata
       if(commVec%nrk_data_max<commVec%nrk_prm(i)%ndata )then
          commVec%nrk_data_max=commVec%nrk_prm(i)%ndata
       endif
    enddo
  end subroutine get_nrk_comm

  subroutine set_comm(&
       commVec,mpi_communicator,&
       prm, &
       Colsort,nCol&
       )
    use mpi
    use comm
    use matrix
    use PARCEL_integer
    implicit none
    type(comm_Vec) :: commVec
    type(prmMat) :: prm
    integer mpi_communicator

    integer(kind=intType) n
    integer(kind=intType) gn
    integer(kind=intType) istart
    integer(kind=intType) nCol
    integer(kind=intType) Colsort(nCol)

    !    integer(kind=intType) myrank,npes
    integer ierr
    integer(kind=intType) i,j,k
    integer nrank,rank
    integer istart_rank
    integer nrk_rank_send
    integer(kind=intType) nrk_data_send
    integer(kind=intType) nrk_data_max

    integer,allocatable :: comm_pair_send(:)
    integer,allocatable :: comm_pair_recv(:)
    integer,allocatable :: s_req(:),r_req(:)

    type send_stat
       integer stat(MPI_STATUS_SIZE)
    end type send_stat
    type recv_stat
       integer stat(MPI_STATUS_SIZE)
    end type recv_stat
    type(send_stat),allocatable :: s_sta(:)
    type(recv_stat),allocatable :: r_sta(:)
    n=prm%n
    gn=prm%gn
    istart=prm%istart

    !    allocate(commVec%buff(n))
    commVec%n=n
    commVec%comm=mpi_communicator
    !    commVec%imem_size=commVec%n*sizeReal8

    call MPI_Comm_size( commVec%comm, commVec%npes,   ierr )
    call MPI_Comm_rank( commVec%comm, commVec%myrank, ierr )
    allocate(commVec%istart(commVec%npes+1))
    do i=1,commVec%npes
       commVec%istart(i)=-1
    enddo
    commVec%istart(1)=istart
    ! 各プロセスの行の開始位置を通信
    call MPI_Allgather( &
         istart, &
         1, &
         !         MPI_INTEGER4, &
         intMPIType, &
         commVec%istart, &
         1, &
         !         MPI_INTEGER4, &
         intMPIType, &
         commVec%comm, &
         ierr)
    commVec%istart(commVec%npes+1)=gn+1
    call get_nrk_comm(commVec,nCol,Colsort)
!     if(commVec%myrank.eq.0)then
!        do i=1,commVec%npes+1
!           write(*,*)i,commVec%istart(i)
!        enddo
!     endif
!    write(*,*)commVec%myrank," nrk_data = ",commVec%nrk_data

#if 1
    ! -----
    allocate(comm_pair_send(commVec%npes))
    allocate(comm_pair_recv(commVec%npes))
    do i=1,commVec%npes
       comm_pair_send(i)=-1
       comm_pair_recv(i)=-1
    enddo
    do i=1,commVec%nrk_rank
       comm_pair_send( commVec%nrk_prm(i)%rank +1 ) = commVec%myrank
    enddo

    ! 送信データを必要とするランクを通信
    ! 送信先のランクでは無い場合は -1
    call MPI_Alltoall(&
         comm_pair_send, &
         1, &
         MPI_INTEGER4, &
         comm_pair_recv, &
         1, &
         MPI_INTEGER4, &
         commVec%comm, &
         ierr)
    nrk_rank_send=0
    do i=1,commVec%npes
       if(comm_pair_recv(i)>-1)then
          nrk_rank_send=nrk_rank_send+1
       endif
    enddo
    commVec%nrk_rank_send=nrk_rank_send
    allocate(commVec%nrk_prm_send(nrk_rank_send))

    j=0
    do i=1,commVec%npes
       if(comm_pair_recv(i)>-1)then
          j=j+1
          commVec%nrk_prm_send(j)%rank=comm_pair_recv(i)
       endif
    enddo

    do i=1,commVec%nrk_rank
       comm_pair_send( commVec%nrk_prm(i)%rank +1 ) = &
            commVec%nrk_prm(i)%ndata
    enddo
    ! 送信するデータ数を通信
    call MPI_Alltoall(comm_pair_send, &
         1, &
         MPI_INTEGER4, &
         !         intMPIType, &
         comm_pair_recv, &
         1, &
         MPI_INTEGER4, &
         !         intMPIType, &
         commVec%comm, &
         ierr)
    j=0
    nrk_data_send=0
    do i=1,commVec%npes
       if(comm_pair_recv(i)>0)then
          j=j+1
          commVec%nrk_prm_send(j)%ndata=comm_pair_recv(i)
          allocate(commVec%nrk_prm_send(j)%offset( &
               commVec%nrk_prm_send(j)%ndata))
          nrk_data_send=nrk_data_send+commVec%nrk_prm_send(j)%ndata
       endif
    enddo
    commVec%nrk_data_send=nrk_data_send
    allocate(s_req(commVec%nrk_rank) )
    allocate(r_req(commVec%nrk_rank_send) )

    allocate(s_sta(commVec%nrk_rank))
    allocate(r_sta(commVec%nrk_rank_send))
    !    allocate(s_sta(MPI_STATUS_SIZE,commVec%nrk_rank) )
    !    allocate(r_sta(MPI_STATUS_SIZE,commVec%nrk_rank_send) )

    ! 送信するデータの位置を各プロセスの行開始位置からのオフセットとして送信
    do i=1,commVec%nrk_rank
       call MPI_Isend( &
            commVec%nrk_prm(i)%offset,&
            commVec%nrk_prm(i)%ndata ,&
            !             MPI_INTEGER4,&
            intMPIType, &
            commVec%nrk_prm(i)%rank,&
            commVec%myrank,&
            commVec%comm,&
            s_req(i),&
            ierr)
    enddo
    do i=1,commVec%nrk_rank_send
       call MPI_Irecv( &
            commVec%nrk_prm_send(i)%offset,&
            commVec%nrk_prm_send(i)%ndata ,&
            !             MPI_INTEGER4, &
            intMPIType, &
            commVec%nrk_prm_send(i)%rank, &
            commVec%nrk_prm_send(i)%rank, &
            commVec%comm, &
            r_req(i), &
            ierr)
    enddo

    do i=1,commVec%nrk_rank
       !        call MPI_Wait(s_req(i),s_sta(1,i),ierr)
       call MPI_Wait(s_req(i),s_sta(i)%stat,ierr)
    enddo
    do i=1,commVec%nrk_rank_send
       call MPI_Wait(r_req(i),r_sta(i)%stat,ierr)
    enddo

    nrk_data_max=0
    do j=1,commVec%nrk_rank_send
       if(nrk_data_max<commVec%nrk_prm_send(j)%ndata)then
          nrk_data_max=commVec%nrk_prm_send(j)%ndata
       endif
    enddo
    commVec%nrk_data_max_send=nrk_data_max

    nrk_data_max=0
    do j=1,commVec%nrk_rank
       if(nrk_data_max<commVec%nrk_prm(j)%ndata)then
          nrk_data_max=commVec%nrk_prm(j)%ndata
       endif
    enddo
    commVec%nrk_data_max=nrk_data_max

    ! if(commVec%myrank.eq.1)then
    !    do j=1,commVec%nrk_rank_send
    !       do i=1,commVec%nrk_prm_send(j)%ndata
    !          write(*,*)&
    !               commVec%nrk_prm_send(j)%rank,&
    !               commVec%nrk_prm_send(j)%offset(i)
    !       enddo
    !    enddo
    ! endif
    ! -----
    !     write(*,*)commVec%myrank,commVec%nrk_rank_send ,commVec%nrk_rank
#else
    ! lock,unlockを利用しない
    call MPI_INFO_CREATE(commVec%win_info,ierr)
    call MPI_INFO_SET(commVec%win_info,"no_locks","true",ierr)   
    call MPI_Win_create( &
         commVec%buff, &
         commVec%imem_size,&
         sizeReal8,&
         commVec%win_info,&
         commVec%comm,&
         commVec%idwin,&
         ierr)
    call MPI_Win_fence(0,commVec%idwin,ierr)
#endif

    deallocate(comm_pair_send)
    deallocate(comm_pair_recv)
    deallocate(s_req)
    deallocate(s_sta)
    deallocate(r_req)
    deallocate(r_sta)

  end subroutine set_comm

#if 0
  subroutine free_window(commVec)
    use comm
    use PARCEL_integer
    implicit none
    type(comm_Vec) :: commVec
    integer i
    integer ierr

    call MPI_Win_fence(0,commVec%idwin,ierr)
    call MPI_Info_free(commVec%win_info,ierr)
    call MPI_Win_free(commVec%idwin,ierr)

  end subroutine free_window
#endif 

  subroutine free_comm(commVec)
    use comm
    use PARCEL_integer
    implicit none
    type(comm_Vec) :: commVec
    integer(kind=intType) i
    integer ierr

    !    deallocate(commVec%buff)
    deallocate(commVec%istart)

    !     do i=1,commVec%nrk_rank
    !        deallocate(commVec%nrk_prm(i)%iLs_Ovl2)
    !        deallocate(commVec%nrk_prm(i)%iLe_Ovl2)
    !     enddo
    do i=1,commVec%nrk_rank
       if(allocated(commVec%nrk_prm(i)%offset))then
          deallocate(commVec%nrk_prm(i)%offset)
       endif

       if(allocated(commVec%nrk_prm(i)%nrk_ptr))then
          deallocate(commVec%nrk_prm(i)%nrk_ptr)
       endif

       if(allocated(commVec%nrk_prm(i)%nrkRow))then
          deallocate(commVec%nrk_prm(i)%nrkRow)
       endif

       if(allocated(commVec%nrk_prm(i)%nrkCol))then
          deallocate(commVec%nrk_prm(i)%nrkCol)
       endif
    enddo
    deallocate(commVec%nrk_prm)

    do i=1,commVec%nrk_rank_send
       deallocate(commVec%nrk_prm_send(i)%offset)
    enddo
    deallocate(commVec%nrk_prm_send)

  end subroutine free_comm

  subroutine nrk_vec_pack_ovl2(nrk_prm,n,ndata,vec,sbuf,nrk_rank,srank)
    use comm
    use PARCEL_integer
    implicit none
    type(nrk) :: nrk_prm(nrk_rank)
    integer(kind=intType) n
    integer ndata
    real*8 vec(n)
    real*8 sbuf(ndata)
    integer nrk_rank
    integer srank

    integer(kind=intType) offset
    integer(kind=intType) i,j,k,i0

    if(srank>-1)then
       j=srank
       !$OMP PARALLEL DO PRIVATE(i,offset) 
       do i=1,nrk_prm(j)%ndata
          offset=nrk_prm(j)%offset(i)
          sbuf(i)=vec(1+offset)
       enddo
       !$OMP END PARALLEL DO
    endif
  end subroutine nrk_vec_pack_ovl2

  subroutine nrk_vec_pack(nrk_prm,n,ndata,vec,sbuf,nrk_rank)
    use comm
    use PARCEL_integer
    implicit none
    integer nrk_rank
    type(nrk) :: nrk_prm(nrk_rank)
    integer(kind=intType) n
    integer ndata
    real*8 vec(n)
    real*8 sbuf(ndata)
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
          sbuf(i+i0)=vec(1+offset)
       enddo
       !$OMP END DO NOWAIT

    enddo
    !$OMP END PARALLEL 

  end subroutine nrk_vec_pack

  subroutine nrk_vec_comm(commVec,sbuf,rbuf)
    use mpi
    use comm
    use PARCEL_integer
    implicit none
    type(comm_Vec) :: commVec

    integer(kind=intType) i,j,k
    integer(kind=intType) i0
    integer ierr

    real*8 sbuf(commVec%nrk_data_send)
    real*8 rbuf(commVec%nrk_data)

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

    !    integer(kind=intType) s_sta(MPI_STATUS_SIZE,commVec%nrk_rank_send)
    !    integer(kind=intType) r_sta(MPI_STATUS_SIZE,commVec%nrk_rank)

    do i=1,commVec%nrk_rank_send
       i0=0
       do k=1,i-1
          i0=i0+commVec%nrk_prm_send(k)%ndata
       enddo
       call MPI_Isend( &
            sbuf(1+i0),&
            commVec%nrk_prm_send(i)%ndata ,&
            MPI_REAL8,&
            commVec%nrk_prm_send(i)%rank,&
            commVec%myrank,&
            commVec%comm,&
            s_req(i),&
            ierr)
    enddo
    !    do i=1,commVec%nrk_rank
    do i=commVec%nrk_rank,1,-1
       i0=0
       do k=1,i-1
          i0=i0+commVec%nrk_prm(k)%ndata
       enddo
       call MPI_Irecv( &
            rbuf(1+i0),&
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
    !    do i=1,commVec%nrk_rank
    do i=commVec%nrk_rank,1,-1
       call MPI_Wait(r_req(i),r_sta(i)%stat,ierr)
    enddo

  end subroutine nrk_vec_comm

  subroutine nrk_vec_comm_ovl2(commVec,sbuf,rbuf,srank,rrank)
    use comm
    use PARCEL_integer
    implicit none
    type(comm_Vec) :: commVec

    integer(kind=intType) i,j,k
    integer(kind=intType) i0
    integer ierr

    real*8 sbuf(commVec%nrk_data_max_send)
    real*8 rbuf(commVec%nrk_data_max)

    integer srank,rrank
    integer s_req
    integer s_sta(MPI_STATUS_SIZE)
    integer r_req
    integer r_sta(MPI_STATUS_SIZE)

    if(srank>-1)then
       i=srank
       ! i0=0
       ! do k=1,i-1
       !    i0=i0+commVec%nrk_prm_send(i)%ndata
       ! enddo
       call MPI_Isend( &
            !            sbuf(1+i0),&
            sbuf,&
            commVec%nrk_prm_send(i)%ndata ,&
            MPI_REAL8,&
            commVec%nrk_prm_send(i)%rank,&
            commVec%myrank,&
            commVec%comm,&
            s_req,&
            ierr)
    endif
    if(rrank>-1)then
       i=rrank
       ! i0=0
       ! do k=1,i-1
       !    i0=i0+commVec%nrk_prm(i)%ndata
       ! enddo
       call MPI_Irecv( &
            !            rbuf(1+i0),&
            rbuf,&
            commVec%nrk_prm(i)%ndata ,&
            MPI_REAL8, &
            commVec%nrk_prm(i)%rank, &
            commVec%nrk_prm(i)%rank, &
            commVec%comm, &
            r_req, &
            ierr)
    endif

    if(srank>-1)then
       call MPI_Wait(s_req,s_sta,ierr)
    endif
    if(rrank>-1)then
       call MPI_Wait(r_req,r_sta,ierr)
    endif

  end subroutine nrk_vec_comm_ovl2

  subroutine nrk_vec_comm_ovl2_isend( &
       commVec,sbuf,srank,s_req,s_sta)
    use comm
    use PARCEL_integer
    implicit none
    type(comm_Vec) :: commVec
    integer(kind=intType) i
    integer ierr
    real*8 sbuf(commVec%nrk_data_max_send)
    integer srank
    integer s_req
    integer s_sta(MPI_STATUS_SIZE)

    if(srank>-1)then
       i=srank
       call MPI_Isend( &
            sbuf,&
            commVec%nrk_prm_send(i)%ndata ,&
            MPI_REAL8,&
            commVec%nrk_prm_send(i)%rank,&
            commVec%myrank,&
            commVec%comm,&
            s_req,&
            ierr)
    endif

  end subroutine nrk_vec_comm_ovl2_isend

  subroutine nrk_vec_comm_ovl2_irecv(&
       commVec,rbuf,rrank,r_req,r_sta)
    use comm
    use PARCEL_integer
    implicit none

    type(comm_Vec) :: commVec
    integer(kind=intType) i
    integer ierr

    real*8 rbuf(commVec%nrk_data_max)
    integer rrank
    integer r_req
    integer r_sta(MPI_STATUS_SIZE)

    if(rrank>-1)then
       i=rrank
       call MPI_Irecv( &
            rbuf,&
            commVec%nrk_prm(i)%ndata ,&
            MPI_REAL8, &
            commVec%nrk_prm(i)%rank, &
            commVec%nrk_prm(i)%rank, &
            commVec%comm, &
            r_req, &
            ierr)
    endif

  end subroutine nrk_vec_comm_ovl2_irecv

  subroutine get_nrk_vec(commVec,vec,rbuf&
       )
    use comm
    use PARCEL_integer
    implicit none
    type(comm_Vec) :: commVec
    real*8 vec(commVec%n)
    real*8 rbuf(commVec%nrk_data)
    real*8 sbuf(commVec%nrk_data_send)
    integer(kind=intType) i,j,k,i0
#if 1
    call nrk_vec_pack( &
         commVec%nrk_prm_send, &
         commVec%n, &
         commVec%nrk_data_send, &
         vec, &
         sbuf, &
         commVec%nrk_rank_send)
    call nrk_vec_comm(commVec,sbuf,rbuf)

#else
    integer(kind=intType) i,j
    integer(kind=intType) icount
    integer(kind=intType) rank
    integer(kind=intType) ierr
    integer(kind=MPI_ADDRESS_KIND) count_send
    integer(kind=MPI_ADDRESS_KIND) data_offset

    do i=1,commVec%n
       commVec%buff(i)=vec(i)
    enddo
    call MPI_Win_fence(0,commVec%idwin,ierr)
    icount=0
    do i=1,commVec%nrk_rank
       rank=commVec%nrk_prm(i)%rank
       do j=1,commVec%nrk_prm(i)%ndata
          data_offset=commVec%nrk_prm(i)%offset(j)
          count_send=1
          icount=icount+count_send
          call MPI_GET(&
               rbuff(icount), &
               count_send, &
               MPI_REAL8,&
               rank, &
               data_offset, &
               count_send, &
               MPI_REAL8,&
               commVec%idwin, &
               ierr)
       enddo
    enddo
    call MPI_Win_fence(0,commVec%idwin,ierr)

#endif

  end subroutine get_nrk_vec

  subroutine get_nrk_vec_ovl2( &
       commVec,vec,rbuff &
       )
    use comm
    use PARCEL_integer
    implicit none
    type(comm_Vec) :: commVec
    real*8 vec(commVec%n)
    real*8 rbuff(commVec%nrk_data)
    integer srank,rrank
    integer(kind=intType) i
    integer(kind=intType) imax
    integer srankmax,rrankmax
    srankmax=commVec%nrk_rank_send
    rrankmax=commVec%nrk_rank
    imax=max(srankmax,rrankmax)

    do i=1,imax
       if(i>srankmax)then
          srank=-1
       else
          srank=i
       endif
       if(i>rrankmax)then
          rrank=-1
       else
          rrank=i
       endif
       ! call get_nrk_vec_ovl_comm( &
       !      commVec,vec,rbuff,&
       !      srank,rrank &
       !      )
    enddo
  end subroutine get_nrk_vec_ovl2

  subroutine get_nrk_vecLU( &
       commVec,rbuff,&
       nColL,nColU,&
       rbuffL,rbuffU&
       )
    use comm
    use PARCEL_integer
    implicit none
    integer(kind=intType) nColL,nColU
    type(comm_Vec) :: commVec
    !    real*8 vec(commVec%n)
    real*8 rbuff(commVec%nrk_data)

    real*8 rbuffL(nColL)
    real*8 rbuffU(nColU)
    integer(kind=intType) i

    !$OMP PARALLEL DO PRIVATE(i)  
    do i=1,nColL
       rbuffL(i)=rbuff(i)
    enddo
    !$OMP END PARALLEL DO 
    !$OMP PARALLEL DO PRIVATE(i)  
    do i=1,nColU
       rbuffU(i)=rbuff(i+nColL)
    enddo
    !$OMP END PARALLEL DO 
  end subroutine get_nrk_vecLU

  subroutine get_nrk_diaVec( &
       commVec,rbuff, &
       diaMat &
       )
    use comm
    use matrix
    use PARCEL_integer
    implicit none
    type(comm_Vec) :: commVec
    real*8 rbuff(commVec%nrk_data)
    type(DIA_Mat) :: diaMat

    integer(kind=intType) i,j,k
    integer(kind=intType) j0,i0
    integer(kind=intType) kL,kU

    integer(kind=intType) nColL,nColD,nColU
    integer(kind=intType) Row

    nColL=diaMat%nColL
    nColD=diaMat%nColD
    nColU=diaMat%nColU

    !     !$OMP PARALLEL DO PRIVATE( &
    !     !$OMP )
#if 0
    i0=0
    do i=1,commVec%nrk_rank
       kL=1
       kU=1
       do j=1,commVec%nrk_prm(i)%ndata
          ! write(*,*)"get rank,offset", &
          !      commVec%nrk_prm(i)%rank,&
          !      commVec%istart( commVec%nrk_prm(i)%rank+1 ) + &
          !      commVec%nrk_prm(i)%offset(j)
          ! ------
          Row = commVec%istart( commVec%nrk_prm(i)%rank+1 ) +  commVec%nrk_prm(i)%offset(j) 

          if( commVec%nrk_prm(i)%rank < commVec%myrank )then
             do k=kL,nColL
                if( diaMat%RowL(k).eq.Row )then
                   i0=i0+1
                   diaMat%Vec(k)=rbuff(i0)
                   kL=k
                   cycle
                endif
             enddo
          else
             do k=kU,nColU
                if(diaMat%RowU(k).eq.Row)then
                   i0=i0+1
                   diaMat%Vec(nColL+nColD+k)=rbuff(i0)
                   kU=k
                   cycle
                endif
             enddo
          endif

       enddo
    enddo
#else

    !$OMP PARALLEL PRIVATE( &
    !$OMP i,j0,j,kL,kU,i0,Row,k &
    !$OMP )

    do i=1,commVec%nrk_rank
       j0=0
       do j=1,i-1
          j0 = j0 + commVec%nrk_prm(i)%ndata
       enddo
       kL=1
       kU=1
       !$OMP DO SCHEDULE(dynamic,nBlock)
       do j=1,commVec%nrk_prm(i)%ndata
          i0=j0+j
          Row = commVec%istart( commVec%nrk_prm(i)%rank+1 ) +  commVec%nrk_prm(i)%offset(j) 

          if( commVec%nrk_prm(i)%rank < commVec%myrank )then
             do k=kL,nColL
                if( diaMat%RowL(k).eq.Row )then
                   diaMat%Vec(k)=rbuff(i0)
                   kL=k
                   cycle
                endif
             enddo
          else
             do k=kU,nColU
                if(diaMat%RowU(k).eq.Row)then
                   diaMat%Vec(nColL+nColD+k)=rbuff(i0)
                   kU=k
                   cycle
                endif
             enddo
          endif
       enddo
       !$OMP END DO NOWAIT
    enddo
    !$OMP END PARALLEL

#endif
    !    write(*,*)diaMat%nColL,diaMat%nColU
    ! do i=1,commVec%nrk_data
    !    write(*,*)"@ recv @",rbuff(i)
    ! enddo

  end subroutine get_nrk_diaVec

  subroutine view_nrk(commVec)
    use comm
    use PARCEL_integer
    implicit none
    type(comm_Vec) :: commVec
    integer irank,jrank
    integer myrank
    integer npes
    integer ierr

    integer recv_rank(commVec%npes)
    integer send_rank(commVec%npes)
    integer(kind=intType) imax

    myrank=commVec%myrank
    npes=commVec%npes

    call MPI_Barrier(commVec%comm,ierr)
!    call system("sleep 1")
    call MPI_Barrier(commVec%comm,ierr)
    do irank=0,npes-1
       if(irank.eq.myrank)then
          recv_rank=-1
          send_rank=-1
          do jrank=1,commVec%nrk_rank
             recv_rank(jrank)=commVec%nrk_prm(jrank)%rank
          enddo
          do jrank=1,commVec%nrk_rank_send
             send_rank(jrank)=commVec%nrk_prm_send(jrank)%rank
          enddo
          imax=max(commVec%nrk_rank,commVec%nrk_rank_send)
          do  jrank=1,imax
             write(*,*)myrank,jrank,send_rank(jrank),recv_rank(jrank)
          enddo
          write(*,*)"-----------------------"
       endif
       call MPI_Barrier(commVec%comm,ierr)
!       call system("sleep 1")
       call MPI_Barrier(commVec%comm,ierr)
    enddo

  end subroutine view_nrk

end module comm_kernel
