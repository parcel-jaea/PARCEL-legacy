module TSQR_kernel
contains

  subroutine TSQR_local_Ronly(n,s,A,R,nvecBlock)
    use omp_lib
    use vec_kernel
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    real*8 A(n*s)
    real*8 R(s*s)
    integer nvecBlock

    integer(kind=intType) i,j,jend,jBlock ,jj
    integer(kind=intType) nQRblock
    integer iRflag
    integer myth
    integer numth
    real*8 R_th(s*s)
    integer max_numth
    integer numblock
    real*8,allocatable :: R_all(:)

    jBlock = nvecBlock + s
    !$OMP PARALLEL 
    !$OMP MASTER
    numth=omp_get_num_threads()
    !$OMP END MASTER
    !$OMP END PARALLEL 

    allocate(R_all(s*s*numth))
    R=0.0d0
    R_all=0.0d0

    !$OMP PARALLEL PRIVATE( &
    !$OMP j,jend,nQRblock,iRflag,&
    !$OMP  myth,R_th )
    myth=omp_get_thread_num()
    R_th=0.0d0

    iRflag=0

    !$OMP DO
    do j=1,n,jBlock 
       jend=min(n,j+jBlock-1)
       nQRblock=jend-j+1
       !        write(*,*)myth,j,jend,nQRblock
       if( (iRflag.eq.0).and.(nQRblock<=s))then
          !        write(*,*)myth,j,jend,nQRblock
          do i=1,s
             do jj=j,jend
                R_th(1+jj-j+(i-1)*s)=A( jj + (i-1)*n )
             enddo
          enddo
       else
          call House_local_Ronly(&
               n,s,A,R_th,&
               j,nQRblock,iRflag)
       endif
       iRflag=1
    enddo
    !$OMP END DO NOWAIT
    do i=1,s
       do j=1,s
          R_all(j+myth*s + (i-1)*s*numth)=R_th(j+(i-1)*s)
       enddo
    enddo
    !$OMP END PARALLEL
    call House_local_Ronly_omp(s*numth,s,R_all,R)

  end subroutine TSQR_local_Ronly

  subroutine TSQR_GLOBAL(comm1,s,R)
    use mpi
    use common_kernel

    use PARCEL_integer
    implicit none
    integer(kind=intType) s
    integer comm1
    integer(kind=intType) i,j,k
    integer(kind=intType) i_,j_

    integer(kind=intType) s2
    integer ierr
    integer irank,nprocess
    integer stage
    integer stat_recv(mpi_status_size),stat_send(mpi_status_size)

    integer req_send,req_recv
    integer tag,nrk
    real*8 rbuff(s*s),sbuff(s*s)
    real*8 A(2*s,s)

    real*8 tmp
    integer itmp
    integer istage,irank_stage,irank0,irank0_master
    real*8 R(s*s)

    real*8 tau(2*s*s)
    real*8 work(2*s*s)
    integer lwork
    integer ss_int4
    lwork=2*s*s

    call MPI_Comm_size(comm1,nprocess, ierr)
    call MPI_Comm_rank(comm1,irank, ierr)

    ! rbuff=-99.0d0
    ! sbuff=irank
    A=0.0d0
    do k=1,s*s
       sbuff(k)=R(k)
    enddo
    do i=1,s
       do j=1,s
          A(j,i)=R(j+s*(i-1))
       enddo
    enddo

    stage=0
    irank0=0
    irank0_master=0

    istage=nprocess
    call floor_log2(istage)

    irank_stage=irank
    do i=istage,0,-1

       if(irank_stage<2**i)then
          ! ----
          ! irank0_master --> 集約された三角行列を持つランク
          ! irank0_master - irank0 --> 集約した三角行列の送信先ランク
          ! -----
          j=0
          stage=1
          do while(stage<2**i)
             j=j+1
             stage=2**j
             !           if(mod(irank,stage).eq.0)then
             if(mod(irank,stage).eq.0)then
#if 0
!! この処理は不要?
                if(irank.eq.(nprocess-1))then
                   nrk=irank0_master
                   stage=2**i
                   !! 送信先が自プロセスでない場合
                   if(nrk.ne.irank)then
                      ! write(*,'(A15,I5.4,I5.4)')"Send stage 1st = ",irank,nrk
                      ss_int4=s*s
                      call MPI_Isend(sbuff,ss_int4, MPI_DOUBLE, nrk, nrk,  comm1, req_send,ierr)
                      call MPI_Wait(req_send, stat_send,ierr)
                   endif

                   exit
                else
#endif
                   nrk=irank+stage/2
                   !! 受信先がある場合
                   if(nrk<(nprocess))then
                      ! write(*,'(A15,I5.4,I5.4)')"Recv stage = ",irank,nrk
                      ss_int4=s*s
!                      ! write(*,'(A15,I5.4,I5.4)')"Recv stage = ",irank,nrk,rbuff
                      call MPI_Irecv(rbuff,ss_int4,MPI_DOUBLE,nrk,irank,comm1,req_recv,ierr)
                      call MPI_Wait(req_recv, stat_recv,ierr)
                      do i_=1,s
                         do j_=1,s
                            A(j_+s,i_)=rbuff(j_+s*(i_-1))
                         enddo
                      enddo
#ifdef TSQR_MGS
                      s2=s*2
                      call MGS_local_Ronly(s2,s,A,R)
#elif TSQR_GEQRF
                      ! 8バイト整数型は未対応
                      call DGEQRF(2*s,s,A,2*s, tau, work, lwork, ierr) 
                      do i_=1,s
                         do j_=1,i_ 
                            R(j_+(i_-1)*s)=A(j_,i_) 
                         enddo
                      enddo
#else
                      s2=s*2
                      !                      call House_local_Ronly_omp(2*s,s,A,R)
                      call House_local_Ronly_unroll(s2,s,A,R)
#endif
                      do i_=1,s
                         do j_=1,s
                            A(j_,i_)=R(j_+s*(i_-1))
                         enddo
                      enddo
                      do k=1,s*s
                         sbuff(k)=R(k)
                      enddo

                   endif
!                endif
             else
                nrk=irank-stage/2
                ! write(*,'(A15,I5.4,I5.4)')"Send stage 2nd = ",irank,nrk   
                ss_int4=s*s
                call MPI_Isend(sbuff, ss_int4, MPI_DOUBLE, nrk, nrk,  comm1, req_send,ierr)
                call MPI_Wait(req_send, stat_send,ierr)
                exit
             endif

             !              ! write(*,'(A15,I5.4,I5.4)')stage,nrk,irank
          end do
          ! -----
          ! ここで1ステージ上のツリーに集約を行う

          if((irank0_master.eq.irank).and.(irank.ne.0))then
             if((irank+2**i)<nprocess)then
                nrk=irank+2**i
                ss_int4=s*s
                ! write(*,'(A15,I5.4,I5.4)')"Recv stage = ",irank,nrk,rbuff
                call MPI_Irecv(rbuff,ss_int4, MPI_DOUBLE, nrk, irank,   comm1, req_recv,ierr)
                call MPI_Wait(req_recv, stat_recv,ierr)             

                do i_=1,s
                   do j_=1,s
                      A(j_+s,i_)=rbuff(j_+s*(i_-1))
                   enddo
                enddo
#ifdef TSQR_MGS
                s2=s*2
                call MGS_local_Ronly(s2,s,A,R)
#elif TSQR_GEQRF
                ! 8バイト整数型は未対応
                call DGEQRF(2*s,s,A,2*s, tau, work, lwork, ierr) 
                do i_=1,s
                   do j_=1,i_ 
                      R(j_+(i_-1)*s)=A(j_,i_) 
                   enddo
                enddo
#else
                !                call House_local_Ronly_omp(2*s,s,A,R)
                s2=s*2
                call House_local_Ronly_unroll(s2,s,A,R)
#endif
                do i_=1,s
                   do j_=1,s
                      A(j_,i_)=R(j_+s*(i_-1))
                   enddo
                enddo
                do k=1,s*s
                   sbuff(k)=R(k)
                enddo

             endif

             nrk=irank0_master - irank0
             ! write(*,'(A15,I5.4,I5.4)')"Send stage 3rd = ",irank,nrk
             ss_int4=s*s
             call MPI_Isend(sbuff, ss_int4, MPI_DOUBLE, nrk, nrk,  comm1, req_send,ierr)
             call MPI_Wait(req_send, stat_send,ierr)
          endif
          if((irank.eq.0).and.(2**(i)<nprocess))then
             nrk=irank+2**(i)
             ss_int4=s*s
             ! write(*,'(A15,I5.4,I5.4)')"Recv stage = ",irank,nrk
             call MPI_Irecv(rbuff, ss_int4, MPI_DOUBLE, nrk, irank,   comm1, req_recv,ierr)
             call MPI_Wait(req_recv, stat_recv,ierr)             

             do i_=1,s
                do j_=1,s
                   A(j_+s,i_)=rbuff(j_+s*(i_-1))
                enddo
             enddo
#ifdef TSQR_MGS
             s2=s*2
             call MGS_local_Ronly(s2,s,A,R)
#elif TSQR_GEQRF
             ! 8バイト整数型は未対応
             call DGEQRF(2*s,s,A,2*s, tau, work, lwork, ierr) 
             do i_=1,s
                do j_=1,i_ 
                   R(j_+(i_-1)*s)=A(j_,i_) 
                enddo
             enddo
#else
             s2=s*2
             !             call House_local_Ronly_omp(2*s,s,A,R)
             call House_local_Ronly_unroll(s2,s,A,R)
#endif
             do i_=1,s
                do j_=1,s
                   A(j_,i_)=R(j_+s*(i_-1))
                enddo
             enddo
             do k=1,s*s
                sbuff(k)=R(k)
             enddo

             !! write(*,'(A15,I5.4,I5.4)')"Recv stage = ",irank,nrk,rbuff
          endif
          exit
       else
          irank0=2**(i)
          irank0_master=irank0_master+irank0
          irank_stage=mod(irank_stage,2**i)
       endif
    enddo
    ss_int4=s*s
    call MPI_Bcast(R,ss_int4,MPI_DOUBLE,0,comm1, ierr)

  end subroutine TSQR_GLOBAL

  subroutine MGS_local_Ronly(n,s,A,R)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n,s
    real*8 tmp1,tmp2,tmp3,tmp4
    real*8  A(n*s),R(s*s)
    real*8  Q(n*s)
    integer(kind=intType) i,j,k

    R=0.0d0
    do j = 1,s
       tmp1=0.0d0
       do k=1,n
          tmp1=tmp1+A(k+(j-1)*n)*A(k+(j-1)*n)
       enddo
       R((j-1)*s+j) = sqrt(tmp1)
       do k=1,n
          Q(k+(j-1)*n)=A(k+(j-1)*n)/R((j-1)*s+j)
       enddo
       if(mod(s-(j+1)+1,2).eq.0)then
          do i = s,j+1,-2
             tmp1=0.0d0
             tmp2=0.0d0
             do k=1,n
                tmp1=tmp1+Q(k+(j-1)*n)*A(k+(i-1)*n)
                tmp2=tmp2+Q(k+(j-1)*n)*A(k+(i-1-1)*n)
             enddo
             R(j+(i-1)*s) = tmp1
             R(j+(i-1-1)*s) = tmp2
             do k=1,n
                A(k+(i-1)*n)=A(k+(i-1)*n)-tmp1*Q(k+(j-1)*n)
                A(k+(i-1-1)*n)=A(k+(i-1-1)*n)-tmp2*Q(k+(j-1)*n)
             enddo
          enddo
       else
          i=s
          tmp1=0.0d0
          do k=1,n
             tmp1=tmp1+Q(k+(j-1)*n)*A(k+(i-1)*n)
          enddo
          R(j+(i-1)*s) = tmp1
          do k=1,n
             A(k+(i-1)*n)=A(k+(i-1)*n)-R(j+(i-1)*s)*Q(k+(j-1)*n)
          enddo
          do i = s-1,j+1,-2
             tmp1=0.0d0
             tmp2=0.0d0
             do k=1,n
                tmp1=tmp1+Q(k+(j-1)*n)*A(k+(i-1)*n)
                tmp2=tmp2+Q(k+(j-1)*n)*A(k+(i-1-1)*n)
             enddo
             R(j+(i-1)*s) = tmp1
             R(j+(i-1-1)*s) = tmp2
             do k=1,n
                A(k+(i-1)*n)=A(k+(i-1)*n)-tmp1*Q(k+(j-1)*n)
                A(k+(i-1-1)*n)=A(k+(i-1-1)*n)-tmp2*Q(k+(j-1)*n)
             enddo
          enddo
       endif
    enddo

  end subroutine MGS_local_Ronly

  subroutine House_local_Ronly(&
       n,s,A,R,&
       kk,nblock,iRflag)
    use PARCEL_integer
    implicit none
    ! kk :: 先頭アドレス
    ! nblock :: QR分解のブロックサイズ
    !--------------------------------------------------------------------------------------------------  
    integer(kind=intType) nHouse,s,n
    real*8 R(s*s)
    real*8 A(n,s)
    integer(kind=intType) kk,nblock
    real*8 A_TSQR(nblock+s,s) ! 
    integer iRflag
    integer(kind=intType) l_
    !--------------------------------------------------------------------------------------------------  
    integer(kind=intType) ii,i,ii_, j, k,k_
    real*8 y1,t1
    real*8 tmp,tmp2,tmp3,tmp4
    real*8 tmp_t1(s),tmp_y1(s)
    !-------------------------------------------------------------------------------------------------- 

    i=1
    tmp2=0.0d0
    do k_=kk,kk+nblock-1
       A_TSQR(k_-kk+1,i)=A(k_,i)
       tmp2=tmp2+A(k_,i)*A(k_,i)
    enddo

    nHouse=nblock
    if(iRflag.eq.1)then
       do l_=1,s
          A_TSQR(l_+nblock,i)=R(l_+(i-1)*s)
          tmp2=tmp2+R(l_+(i-1)*s)*R(l_+(i-1)*s)
       enddo
       nHouse=nHouse+s
    endif
    tmp2=tmp2-A_TSQR(1,1)*A_TSQR(1,1) 
    ! -----
    i=2
    do k_=kk,kk+nblock-1
       A_TSQR(k_-kk+1,i)=A(k_,i)
    enddo
    if(iRflag.eq.1)then
       do l_=1,s
          A_TSQR(l_+nblock,i)=R(l_+(i-1)*s)
       enddo
    endif

    do i=1,s

       tmp3=tmp2+A_TSQR(i,i)*A_TSQR(i,i)
       tmp3=sqrt(tmp3)
#if 1
       if(A_TSQR(i,i)>0.0d0)then
          y1=A_TSQR(i,i)+tmp3
       else
          y1=A_TSQR(i,i)-tmp3
       endif
#else
       y1=A_TSQR(i,i)-tmp3
#endif
       tmp=tmp2+y1*A_TSQR(i,i)
       tmp2=tmp2+y1*y1 

       if(tmp2.eq.0.0d0)then
          tmp2=1.0d0-300
       endif
       t1=2.0d0/tmp2
       A_TSQR(i,i)=A_TSQR(i,i)-t1*y1*tmp
       tmp_y1(i)=y1     
       tmp_t1(i)=t1

       if(i.ne.s)then
          tmp=tmp_y1(i)*A_TSQR(i,i+1)
          do j=i+1,nHouse
             tmp=tmp+A_TSQR(j,i)*A_TSQR(j,i+1)
          enddo
          tmp=tmp_t1(i)*tmp
          A_TSQR(i,i+1)=A_TSQR(i,i+1)-tmp_y1(i)*tmp
          A_TSQR(i+1,i+1)=A_TSQR(i+1,i+1)-A_TSQR(i+1,i)*tmp

          tmp4=0.0d0
          do j=i+2,nHouse
             A_TSQR(j,i+1)=A_TSQR(j,i+1)-A_TSQR(j,i)*tmp
             tmp4=tmp4+A_TSQR(j,i+1)*A_TSQR(j,i+1)
          enddo

          k=i+2
          if(k<=s)then
             do k_=kk,kk+nblock-1
                A_TSQR(k_-kk+1,k)=A(k_,k)
             enddo
             if(iRflag.eq.1)then
                do l_=1,s
                   A_TSQR(l_+nblock,k)=R(l_+(k-1)*s)
                enddo
             endif
             ! -----
             if(mod(i,2).eq.0)then
                ii_=1
             else
                ii=1
                tmp=tmp_y1(ii)*A_TSQR(ii,k)
                do j=ii+1,nHouse
                   tmp=tmp+A_TSQR(j,ii)*A_TSQR(j,k)
                enddo
                tmp=tmp_t1(ii)*tmp
                A_TSQR(ii,k)=A_TSQR(ii,k)-tmp_y1(ii)*tmp
                do j=ii+1,nHouse
                   A_TSQR(j,k)=A_TSQR(j,k)-A_TSQR(j,ii)*tmp
                enddo
                ii_=2
             endif

             ii=ii_
             tmp=tmp_y1(ii)*A_TSQR(ii,k)
             do j=ii+1,nHouse
                tmp=tmp+A_TSQR(j,ii)*A_TSQR(j,k)
             enddo
             do ii=ii_,i,2
                tmp=tmp_t1(ii)*tmp
                A_TSQR(ii,k)=A_TSQR(ii,k)-tmp_y1(ii)*tmp
                A_TSQR(ii+1,k)=A_TSQR(ii+1,k)-A_TSQR(ii+1,ii)*tmp
                tmp2=tmp_y1(ii+1)*A_TSQR(ii+1,k)
                do j=ii+2,nHouse
                   A_TSQR(j,k)=A_TSQR(j,k)-A_TSQR(j,ii)*tmp
                   tmp2=tmp2+A_TSQR(j,ii+1)*A_TSQR(j,k)
                enddo
                tmp2=tmp_t1(ii+1)*tmp2
                if(ii+2<i)then
                   A_TSQR(ii+1,k)=A_TSQR(ii+1,k)-tmp_y1(ii+1)*tmp2
                   A_TSQR(ii+2,k)=A_TSQR(ii+2,k)-A_TSQR(ii+2,ii+1)*tmp2
                   tmp=tmp_y1(ii+2)*A_TSQR(ii+2,k)
                   do j=ii+3,nHouse
                      A_TSQR(j,k)=A_TSQR(j,k)-A_TSQR(j,ii+1)*tmp2
                      tmp=tmp+A_TSQR(j,ii+2)*A_TSQR(j,k)
                   enddo
                else
                   A_TSQR(ii+1,k)=A_TSQR(ii+1,k)-tmp_y1(ii+1)*tmp2
                   do j=ii+2,nHouse
                      A_TSQR(j,k)=A_TSQR(j,k)-A_TSQR(j,ii+1)*tmp2
                   enddo
                endif
             enddo
             ! -----
          endif
          tmp2=tmp4
       endif
    enddo

    R=0.0d0
    do i=1,s
       do j=i,s
          R(i+(j-1)*s)=A_TSQR(i,j)
       enddo
    enddo

  end subroutine House_local_Ronly

  subroutine House_local_Ronly_org(n,s,A,R)
    use PARCEL_integer
    implicit none
    !--------------------------------------------------------------------------------------------------  
    integer(kind=intType)  n,s
    real*8 R(s*s)
    real*8 A(n,s)       

    !--------------------------------------------------------------------------------------------------  
    integer(kind=intType) ii,kk,i,ii_, j, k,k_,ss,ss_
    integer(kind=intType) ii_block,k_block
    real*8 y1,t1
    real*8 y2,t2
    real*8 tmp,tmp2,tmp3,tmp4
    real*8 tmp_t1(s),tmp_tmp(s,s),tmp_y1(s)
    !-------------------------------------------------------------------------------------------------- 
    do i=1,s

       tmp2=0.0d0
       do j=i+1,n
          tmp2=tmp2+A(j,i)*A(j,i) 
       enddo
       tmp3=tmp2+A(i,i)*A(i,i)
       tmp3=sqrt(tmp3)
#if 1
       if(A(i,i)>0.0d0)then
          y1=A(i,i)+tmp3
       else
          y1=A(i,i)-tmp3
       endif
#else
       y1=A(i,i)-tmp3
#endif
       tmp=tmp2+y1*A(i,i)
       tmp2=tmp2+y1*y1 

       if(tmp2.eq.0.0d0)then
          tmp2=1.0d0-300
       endif
       t1=2.0d0/tmp2
       A(i,i)=A(i,i)-t1*y1*tmp
       tmp_y1(i)=y1     
       tmp_t1(i)=t1

       if(i.ne.s)then
          tmp=tmp_y1(i)*A(i,i+1)
          do j=i+1,n
             tmp=tmp+A(j,i)*A(j,i+1)
          enddo
          tmp=tmp_t1(i)*tmp

          A(i,i+1)=A(i,i+1)-tmp_y1(i)*tmp
          A(i+1,i+1)=A(i+1,i+1)-A(i+1,i)*tmp
          do j=i+2,n
             A(j,i+1)=A(j,i+1)-A(j,i)*tmp
          enddo
          k=i+2
          if(k<=s)then
             do ii=1,i
                tmp=tmp_y1(ii)*A(ii,k)
                do j=ii+1,n
                   tmp=tmp+A(j,ii)*A(j,k)
                enddo
                tmp=tmp_t1(ii)*tmp
                do j=ii+1,n
                   A(j,k)=A(j,k)-A(j,ii)*tmp
                enddo
                A(ii,k)=A(ii,k)-tmp_y1(ii)*tmp
             enddo
          endif
       endif
    enddo

    R=0.0d0
    do i=1,s
       do j=i,s
          R(i+(j-1)*s)=A(i,j)
       enddo
    enddo

  end subroutine House_local_Ronly_org

  subroutine House_local_Ronly_omp(n,s,A,R)
    use PARCEL_integer
    implicit none
    !--------------------------------------------------------------------------------------------------  
    integer(kind=intType)  n,s
    real*8 R(s*s)
    real*8 A(n,s)
    !--------------------------------------------------------------------------------------------------  
    integer(kind=intType) ii,kk,i,ii_, j, k,k_,ierr,ss,ss_
    integer(kind=intType) ii_block,k_block
    real*8 y1,t1
    real*8 y2,t2
    real*8 tmp,tmp2,tmp3,tmp4
    real*8 tmp_t1(s),tmp_tmp(s,s),tmp_y1(s)
    !-------------------------------------------------------------------------------------------------- 
    do i=1,s

       tmp2=0.0d0
       !$OMP PARALLEL DO PRIVATE(j) &
       !$OMP REDUCTION(+:tmp2) 
       do j=i+1,n
          tmp2=tmp2+A(j,i)*A(j,i) 
       enddo
       tmp3=tmp2+A(i,i)*A(i,i)
       tmp3=sqrt(tmp3)

#if 1
       if(A(i,i)>0.0d0)then
          y1=A(i,i)+tmp3
       else
          y1=A(i,i)-tmp3
       endif
#else
       y1=A(i,i)-tmp3
#endif

       tmp=tmp2+y1*A(i,i)
       tmp2=tmp2+y1*y1 
       if(tmp2.eq.0.0d0)then
          tmp2=1.0d0-300
       endif
       t1=2.0d0/tmp2

       A(i,i)=A(i,i)-t1*y1*tmp
       tmp_y1(i)=y1     
       tmp_t1(i)=t1

       if(i.ne.s)then
          tmp=tmp_y1(i)*A(i,i+1)
          !$OMP PARALLEL DO PRIVATE(j) &
          !$OMP REDUCTION(+:tmp) 
          do j=i+1,n
             tmp=tmp+A(j,i)*A(j,i+1)
          enddo
          tmp=tmp_t1(i)*tmp

          A(i,i+1)=A(i,i+1)-tmp_y1(i)*tmp
          A(i+1,i+1)=A(i+1,i+1)-A(i+1,i)*tmp
          !$OMP PARALLEL DO PRIVATE(j) 
          do j=i+2,n
             A(j,i+1)=A(j,i+1)-A(j,i)*tmp
          enddo
          k=i+2
          if(k<=s)then
             do ii=1,i
                tmp=tmp_y1(ii)*A(ii,k)
                !$OMP PARALLEL DO PRIVATE(j) &
                !$OMP REDUCTION(+:tmp) 
                do j=ii+1,n
                   tmp=tmp+A(j,ii)*A(j,k)
                enddo
                tmp=tmp_t1(ii)*tmp
                !$OMP PARALLEL DO PRIVATE(j) 
                do j=ii+1,n
                   A(j,k)=A(j,k)-A(j,ii)*tmp
                enddo
                A(ii,k)=A(ii,k)-tmp_y1(ii)*tmp
             enddo
          endif
       endif
    enddo
    R=0.0d0
    do i=1,s
       do j=i,s
          R(i+(j-1)*s)=A(i,j)
       enddo
    enddo

  end subroutine House_local_Ronly_omp

  subroutine House_local_Ronly_unroll(n,s,A,R)
    use PARCEL_integer
    implicit none
    !--------------------------------------------------------------------------------------------------  
    integer(kind=intType) n,s
    real*8 R(s*s)
    real*8 A(n,s)       

    !--------------------------------------------------------------------------------------------------  
    integer(kind=intType) i,j,k,ss,ss_
    real*8 y1,t1
    real*8 y2,t2
    real*8 tmp,tmp2,tmp3,tmp4
    !-------------------------------------------------------------------------------------------------- 

    if(mod((s-1+1),2).eq.0)then
       ss=1
    else
       ss=2

       i=1
       tmp2=0.0d0
       do j=n,i+1,-1
          tmp2=tmp2+A(j,i)*A(j,i) 
       enddo
       j=i
       tmp3=tmp2+A(j,i)*A(j,i)
       tmp3=sqrt(tmp3)
       y1=A(j,i)-tmp3

       tmp=tmp2+y1*A(j,j)
       tmp2=tmp2+y1*y1 

       if(tmp2.eq.0.0d0)then
          tmp2=1.0d0-300
       endif
       t1=2.0d0/tmp2
       A(j,j)=A(j,j)-t1*y1*tmp

       if(mod(s-(i+1)+1,2).eq.0)then
          do k=s,i+1,-2
             j=i
             tmp=y1*A(j,k)
             tmp2=y1*A(j,k-1)
             do j=i+1,n
                tmp=tmp+A(j,i)*A(j,k)
                tmp2=tmp2+A(j,i)*A(j,k-1)
             enddo
             j=i
             A(j,k)=A(j,k)-t1*y1*tmp
             A(j,k-1)=A(j,k-1)-t1*y1*tmp2
             do j=i+1,n
                A(j,k)=A(j,k)-t1*A(j,i)*tmp
                A(j,k-1)=A(j,k-1)-t1*A(j,i)*tmp2
             enddo
          enddo
       else
          k=s
          j=i
          tmp=y1*A(j,k)
          do j=i+1,n
             tmp=tmp+A(j,i)*A(j,k)
          enddo
          j=i
          A(j,k)=A(j,k)-t1*y1*tmp
          do j=i+1,n
             A(j,k)=A(j,k)-t1*A(j,i)*tmp
          enddo
          do k=s-1,i+1,-2
             j=i
             tmp=y1*A(j,k)
             tmp2=y1*A(j,k-1)
             do j=i+1,n
                tmp=tmp+A(j,i)*A(j,k)
                tmp2=tmp2+A(j,i)*A(j,k-1)
             enddo
             j=i
             A(j,k)=A(j,k)-t1*y1*tmp
             A(j,k-1)=A(j,k-1)-t1*y1*tmp2
             do j=i+1,n
                A(j,k)=A(j,k)-t1*A(j,i)*tmp
                A(j,k-1)=A(j,k-1)-t1*A(j,i)*tmp2
             enddo
          enddo
       endif
    endif

    do i=ss,s,2

       tmp2=0.0d0
       do j=n,i+1,-1
          tmp2=tmp2+A(j,i)*A(j,i) 
       enddo
       tmp3=tmp2+A(i,i)*A(i,i)
       tmp3=sqrt(tmp3)
       y1=A(i,i)-tmp3
       tmp=tmp2+y1*A(i,i)
       tmp2=tmp2+y1*y1 
       if(tmp2.eq.0.0d0)then
          tmp2=1.0d0-300
       endif
       t1=2.0d0/tmp2
       A(i,i)=A(i,i)-t1*y1*tmp
       tmp=y1*A(i,i+1)
       do j=i+1,n
          tmp=tmp+A(j,i)*A(j,i+1)
       enddo
       A(i,i+1)=A(i,i+1)-t1*y1*tmp
       A(i+1,i+1)=A(i+1,i+1)-t1*A(i+1,i)*tmp
       tmp2=0.0d0
       do j=i+2,n
          A(j,i+1)=A(j,i+1)-t1*A(j,i)*tmp
          tmp2=tmp2+A(j,i+1)*A(j,i+1) 
       enddo
       ! ---
       !       tmp2=0.0d0 
       !       tmp=0.0d0
       !       do j=i+1,n
       !          tmp2=tmp2+A(j,i)*A(j,i) 
       !          tmp=tmp+A(j,i)*A(j,i+1)
       !       enddo
       !       tmp3=tmp2+A(i,i)*A(i,i)
       !       tmp3=sqrt(tmp3)
       !       y1=A(i,i)-tmp3

       !       tmp3=tmp2+y1*A(i,i)
       !       tmp2=tmp2+y1*y1 
       !       t1=2.0d0/tmp2
       !       A(i,i)=A(i,i)-t1*y1*tmp3

       !       tmp=tmp+y1*A(i,i+1)
       !       A(i,i+1)=A(i,i+1)-t1*y1*tmp
       !       A(i+1,i+1)=A(i+1,i+1)-t1*A(i+1,i)*tmp
       !       tmp2=0.0d0
       !       do j=i+2,n
       !          A(j,i+1)=A(j,i+1)-t1*A(j,i)*tmp
       !          tmp2=tmp2+A(j,i+1)*A(j,i+1) 
       !       enddo

       ! -----

       tmp3=tmp2+A(i+1,i+1)*A(i+1,i+1)
       tmp3=sqrt(tmp3)
       y2=A(i+1,i+1)-tmp3
       tmp=tmp2+y2*A(i+1,i+1)
       tmp2=tmp2+y2*y2 
       if(tmp2.eq.0.0d0)then
          tmp2=1.0d0-300
       endif
       t2=2.0d0/tmp2
       A(i+1,i+1)=A(i+1,i+1)-t2*y2*tmp

       if(mod(s-(i+2)+1,2).eq.0)then
          ss_=s
       else
          ss_=s+1
          k=s
          tmp=y1*A(i,k)
          do j=i+1,n
             tmp=tmp+A(j,i)*A(j,k)
          enddo
          A(i,k)=A(i,k)-t1*y1*tmp
          A(i+1,k)=A(i+1,k)-t1*A(i+1,i)*tmp

          tmp2=0.0d0
          do j=i+2,n
             A(j,k)=A(j,k)-t1*A(j,i)*tmp
             tmp2=tmp2+A(j,i+1)*A(j,k)
          enddo
          tmp2=tmp2+y2*A(i+1,k)

          A(i+1,k)=A(i+1,k)-t2*y2*tmp2
          do j=i+2,n
             A(j,k)=A(j,k)-t2*A(j,i+1)*tmp2
          enddo
       endif
       do k=ss_,i+2,-2
          tmp=y1*A(i,k)
          tmp2=y1*A(i,k-1)
          do j=i+1,n
             tmp=tmp+A(j,i)*A(j,k)
             tmp2=tmp2+A(j,i)*A(j,k-1)
          enddo
          A(i,k)=A(i,k)-t1*y1*tmp
          A(i+1,k)=A(i+1,k)-t1*A(i+1,i)*tmp

          A(i,k-1)=A(i,k-1)-t1*y1*tmp2
          A(i+1,k-1)=A(i+1,k-1)-t1*A(i+1,i)*tmp2

          tmp3=0.0d0 
          tmp4=0.0d0
          do j=i+2,n
             A(j,k)=A(j,k)-t1*A(j,i)*tmp
             tmp3=tmp3+A(j,i+1)*A(j,k)

             A(j,k-1)=A(j,k-1)-t1*A(j,i)*tmp2
             tmp4=tmp4+A(j,i+1)*A(j,k-1)
          enddo
          tmp3=tmp3+y2*A(i+1,k)
          tmp4=tmp4+y2*A(i+1,k-1)

          A(i+1,k)=A(i+1,k)-t2*y2*tmp3
          A(i+1,k-1)=A(i+1,k-1)-t2*y2*tmp4
          do j=i+2,n
             A(j,k)=A(j,k)-t2*A(j,i+1)*tmp3
             A(j,k-1)=A(j,k-1)-t2*A(j,i+1)*tmp4
          enddo
       enddo

    enddo

    R=0.0d0
    do i=1,s
       do j=i,s
          R(i+(j-1)*s)=A(i,j)
       enddo
    enddo

  end subroutine House_local_Ronly_unroll
end module TSQR_kernel

