!   subroutine dd_add(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo)
!   subroutine dd_mul(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo)
!   subroutine dd_div(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo)
!   subroutine dd_sub(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo)

module ddvec_kernel
  use dd_func
  implicit none
  integer,parameter :: nvecBlock = 200
contains
  ! y(:)=alpha*x(:)+beta*y(:)
  subroutine ddaxpy(&
       n,&
       alpha_hi,alpha_lo, &
       beta_hi,beta_lo, &
       x_hi,x_lo, &
       y_hi,y_lo)
  use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 alpha_hi,alpha_lo
    real*8 beta_hi,beta_lo
    real*8 x_hi(n)
    real*8 x_lo(n)
    real*8 y_hi(n)
    real*8 y_lo(n)
    integer(kind=intType) i 
    real*8 b_hi,b_lo
    real*8 c_hi,c_lo

    !$OMP PARALLEL DO PRIVATE(i,b_hi,b_lo,c_hi,c_lo )
    do i=1,n
!       y(i)=alpha*x(i)+beta*y(i)
       call dd_mul(beta_hi,beta_lo,y_hi(i),y_lo(i),c_hi,c_lo)
       call dd_mul(alpha_hi,alpha_lo,x_hi(i),x_lo(i),b_hi,b_lo)
       call dd_add(b_hi,b_lo,c_hi,c_lo,y_hi(i),y_lo(i))
    enddo
    !$OMP END PARALLEL DO 
  end subroutine ddaxpy

  ! y(:)=alpha*x(:)
  subroutine ddinplace_vec(&
       n,&
       alpha_hi,alpha_lo,&
       x_hi,x_lo,&
       y_hi,y_lo)
  use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 alpha_hi,alpha_lo
    real*8 x_hi(n),x_lo(n)
    real*8 y_hi(n),y_lo(n)
    integer(kind=intType) i 

    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
!       y(i)=alpha*x(i)
       call dd_mul(alpha_hi,alpha_lo,x_hi(i),x_lo(i),y_hi(i),y_lo(i))
    enddo
    !$OMP END PARALLEL DO 
  end subroutine ddinplace_vec


  ! y(:)=x(:)
  subroutine ddxcpy(&
       n,&
       x_hi,x_lo,&
       y_hi,y_lo)
  use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 x_hi(n)
    real*8 x_lo(n)
    real*8 y_hi(n)
    real*8 y_lo(n)
    real*8 alpha_hi,alpha_lo
    real*8 beta_hi,beta_lo
    alpha_hi=1.0d0
    alpha_lo=0.0d0
    beta_hi=0.0d0
    beta_lo=0.0d0
    call ddaxpy( &
       n,&
       alpha_hi,alpha_lo, &
       beta_hi,beta_lo, &
       x_hi,x_lo, &
       y_hi,y_lo)
  end subroutine ddxcpy

  ! dotp=<x,y>
  subroutine ddcalc_dot( &
       n,icomm, &
       x_hi,x_lo, &
       y_hi,y_lo, &
       dotp_hi,dotp_lo)
    use mpi
  use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer icomm
    real*8 x_hi(n)
    real*8 x_lo(n)
    real*8 y_hi(n)
    real*8 y_lo(n)
    real*8 dotp_hi
    real*8 dotp_lo
    real*8 dotp_l_hi
    real*8 dotp_l_lo
    integer(kind=intType) i 
    integer ierr
    real*8 c_hi,c_lo
    real*8 tmp_hi,tmp_lo

    dotp_hi=0.0d0
    dotp_lo=0.0d0
    dotp_l_hi=0.0d0
    dotp_l_lo=0.0d0
    !$OMP PARALLEL PRIVATE(i,c_hi,c_lo,tmp_hi,tmp_lo)
    tmp_hi=0.0d0
    tmp_lo=0.0d0
    !$OMP DO
    do i=1,n
!       dotp_l = dotp_l+x(i)*y(i)
       call dd_mul(x_hi(i),x_lo(i),y_hi(i),y_lo(i),c_hi,c_lo)
       call dd_add(tmp_hi,tmp_lo,c_hi,c_lo,tmp_hi,tmp_lo)
    enddo
    !$OMP END DO NOWAIT

    !$OMP CRITICAL
    call dd_add(dotp_l_hi,dotp_l_lo,tmp_hi,tmp_lo,dotp_l_hi,dotp_l_lo)
    !$OMP END CRITICAL

    !$OMP END PARALLEL  

    ! call MPI_Allreduce(dotp_l,dotp,1,MPI_REAL8,MPI_SUM,comm,ierr)    
    call MPI_AllReduce_dd( &
         dotp_l_hi,dotp_l_lo,dotp_hi,dotp_lo,&
         1,MPI_REAL8,MPI_SUM,icomm,ierr)

  end subroutine ddcalc_dot

  ! dotpij=<xi,yj>
  subroutine ddcalc_gemm_dot( &
       n,s,t,icomm,&
       x_hi,x_lo,&
       y_hi,y_lo,&
       dotp_hi,dotp_lo)
    use mpi
  use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer icomm
    integer(kind=intType) s,t
    real*8 x_hi(n*s)
    real*8 x_lo(n*s)
    real*8 y_hi(n*t)
    real*8 y_lo(n*t)
    real*8 dotp_hi(s*t)
    real*8 dotp_lo(s*t)
    real*8 dotp_l_hi(s*t)
    real*8 dotp_l_lo(s*t)
    real*8 c_hi,c_lo

    integer(kind=intType) ii,i,j,k
    integer ierr
    real*8 tmp_hi(s*t),tmp_lo(s*t)
    integer st

    dotp_hi=0.0d0
    dotp_lo=0.0d0
    dotp_l_hi=0.0d0
    dotp_l_lo=0.0d0

    !$OMP PARALLEL PRIVATE(ii,k,j,i,tmp_hi,tmp_lo,c_hi,c_lo)
    tmp_hi=0.0d0
    tmp_lo=0.0d0
    !$OMP DO
    do ii=1,n,nvecBlock
       do k=1,t
          do j=1,s
             do i=ii,min(n,ii+nvecBlock-1)
#if 0
                ! tmp_hi(j+s*(k-1)) = tmp_hi(j+s*(k-1)) + &
                !      x_hi(i+n*(j-1))*y_hi(i+n*(k-1))
#else
                call dd_mul(&
                     x_hi(i+n*(j-1)),x_lo(i+n*(j-1)),&
                     y_hi(i+n*(k-1)),y_lo(i+n*(k-1)),&
                     c_hi,c_lo)
                call dd_add( &
                     tmp_hi(j+s*(k-1)), &
                     tmp_lo(j+s*(k-1)), &
                     c_hi, &
                     c_lo, &
                     tmp_hi(j+s*(k-1)), &
                     tmp_lo(j+s*(k-1)))
#endif
             enddo
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP CRITICAL
    do k=1,t
       do j=1,s
          call dd_add( &
               dotp_l_hi(j+s*(k-1)), &
               dotp_l_lo(j+s*(k-1)), &
               tmp_hi(j+s*(k-1)), &
               tmp_lo(j+s*(k-1)), &
               dotp_l_hi(j+s*(k-1)), &
               dotp_l_lo(j+s*(k-1)) )
       enddo
    enddo
    !$OMP END CRITICAL

    !$OMP END PARALLEL 
    st=s*t
    call MPI_AllReduce_dd( &
         dotp_l_hi,dotp_l_lo,dotp_hi,dotp_lo,&
         st,MPI_REAL8,MPI_SUM,icomm,ierr)
  end subroutine ddcalc_gemm_dot

! Todo 通信が上位ビット、下位ビット個別に送付するようになっている 
!          --> 一度の通信で完了するように修正
  subroutine set_ddAS_vec( &
       n,comm, &
       addL, &
       myrank,Urank,Lrank, &
       x_hi,x_lo,y_hi,y_lo)
    use mpi
    use vec_kernel
  use PARCEL_integer
    implicit none
    integer comm
    integer(kind=intType) n
    integer addL
    integer myrank,Urank,Lrank
    real*8 x_hi(n)
    real*8 x_lo(n)
    real*8 y_hi(n+addL*2)
    real*8 y_lo(n+addL*2)

    integer(kind=intType) i 
    real*8 r_buffL_hi(addL),r_buffU_hi(addL)
    real*8 r_buffL_lo(addL),r_buffU_lo(addL)
    real*8 s_buffL_hi(addL),s_buffU_hi(addL)
    real*8 s_buffL_lo(addL),s_buffU_lo(addL)

    integer s_Ureq,s_Lreq
    integer r_Ureq,r_Lreq
    integer s_Usta(MPI_STATUS_SIZE),s_Lsta(MPI_STATUS_SIZE)
    integer r_Usta(MPI_STATUS_SIZE),r_Lsta(MPI_STATUS_SIZE)
    integer ierr
    integer(kind=intType) i0
    if(Urank.ne.-1)then
       do i=1,addL
          s_buffU_hi(i)=x_hi(i)
       enddo
       r_buffU_hi=0.0d0
       call MPI_Irecv(r_buffU_hi,addL,MPI_REAL8,Urank,Urank ,comm,r_Ureq,ierr)
       call MPI_Isend(s_buffU_hi,addL,MPI_REAL8,Urank,myrank,comm,s_Ureq,ierr)
    endif
    if(Lrank.ne.-1)then
       do i=1,addL
          s_buffL_hi(i)=x_hi(n-addL+i)
       enddo
       r_buffL_hi=0.0d0
       call MPI_Irecv(r_buffL_hi,addL,MPI_REAL8,Lrank,Lrank ,comm,r_Lreq,ierr)
       call MPI_Isend(s_buffL_hi,addL,MPI_REAL8,Lrank,myrank,comm,s_Lreq,ierr)
    endif

    i0=0
    if(Urank.ne.-1)then
       call MPI_Wait(r_Ureq,r_Usta,ierr)
       call MPI_Wait(s_Ureq,s_Usta,ierr)
       do i=1,addL
          y_hi(i+i0)=r_buffU_hi(i)
       enddo
       i0=i0+addL
    endif

    call xcpy(n,x_hi,y_hi(1+i0))
    call xcpy(n,x_lo,y_lo(1+i0))
    i0=i0+n

    if(Lrank.ne.-1)then
       call MPI_Wait(r_Lreq,r_Lsta,ierr)
       call MPI_Wait(s_Lreq,s_Lsta,ierr)
       do i=1,addL
          y_hi(i0+i)=r_buffL_hi(i)
       enddo
    endif
! ------
    if(Urank.ne.-1)then
       do i=1,addL
          s_buffU_lo(i)=x_lo(i)
       enddo
       r_buffU_lo=0.0d0
       call MPI_Irecv(r_buffU_lo,addL,MPI_REAL8,Urank,Urank ,comm,r_Ureq,ierr)
       call MPI_Isend(s_buffU_lo,addL,MPI_REAL8,Urank,myrank,comm,s_Ureq,ierr)
    endif
    if(Lrank.ne.-1)then
       do i=1,addL
          s_buffL_lo(i)=x_lo(n-addL+i)
       enddo
       r_buffL_lo=0.0d0
       call MPI_Irecv(r_buffL_lo,addL,MPI_REAL8,Lrank,Lrank ,comm,r_Lreq,ierr)
       call MPI_Isend(s_buffL_lo,addL,MPI_REAL8,Lrank,myrank,comm,s_Lreq,ierr)
    endif

    i0=0
    if(Urank.ne.-1)then
       call MPI_Wait(r_Ureq,r_Usta,ierr)
       call MPI_Wait(s_Ureq,s_Usta,ierr)
       do i=1,addL
          y_lo(i+i0)=r_buffU_lo(i)
       enddo
       i0=i0+addL
    endif
    call xcpy(n,x_hi,y_hi(1+i0))
    call xcpy(n,x_lo,y_lo(1+i0))
    i0=i0+n

    if(Lrank.ne.-1)then
       call MPI_Wait(r_Lreq,r_Lsta,ierr)
       call MPI_Wait(s_Lreq,s_Lsta,ierr)
       do i=1,addL
          y_lo(i0+i)=r_buffL_lo(i)
       enddo
    endif

  end subroutine set_ddAS_vec

  subroutine get_ddAS_vec( &
       n,comm, &
       addL, &
       myrank,Urank,Lrank, &
       x_hi,x_lo, &
       vecL_hi,vecL_lo,&
       vecU_hi,vecU_lo)
!       x,y)
    use mpi
  use PARCEL_integer
    implicit none
    integer comm
    integer(kind=intType) n
    integer addL
    integer myrank,Urank,Lrank
    real*8 x_hi(n+addL*2),x_lo(n+addL*2)
    real*8 vecL_hi(addL),vecL_lo(addL)
    real*8 vecU_hi(addL),vecU_lo(addL)

!    real*8 y(n)
    integer iflag ! オーバーラップ領域の補間方法

    integer(kind=intType) i 
    real*8 r_buffL_hi(addL),r_buffU_hi(addL)
    real*8 r_buffL_lo(addL),r_buffU_lo(addL)
    real*8 s_buffL_hi(addL),s_buffU_hi(addL)
    real*8 s_buffL_lo(addL),s_buffU_lo(addL)

    integer s_Ureq,s_Lreq
    integer r_Ureq,r_Lreq
    integer s_Usta(MPI_STATUS_SIZE),s_Lsta(MPI_STATUS_SIZE)
    integer r_Usta(MPI_STATUS_SIZE),r_Lsta(MPI_STATUS_SIZE)
    integer ierr
    integer(kind=intType) i0
 
    ! ---
    ! iflag = 0 :: 相加平均
    ! iflag = 1 :: 破棄
    ! ---

!    call xcpy(n,x(1+i0),y)
    iflag=0
    if(iflag.eq.1)then
       return
    endif

    if(Urank.ne.-1)then
       i0=addL
    else
       i0=0
    endif
    if(Urank.ne.-1)then
       do i=1,addL
          s_buffU_hi(i)=x_hi(i)
       enddo
       r_buffU_hi=0.0d0
       call MPI_Irecv(r_buffU_hi,addL,MPI_REAL8,Urank,Urank ,comm,r_Ureq,ierr)
       call MPI_Isend(s_buffU_hi,addL,MPI_REAL8,Urank,myrank,comm,s_Ureq,ierr)
    endif

    if(Lrank.ne.-1)then
       do i=1,addL
          s_buffL_hi(i)=x_hi(i0+n+i)
       enddo
       r_buffL_hi=0.0d0
       call MPI_Irecv(r_buffL_hi,addL,MPI_REAL8,Lrank,Lrank ,comm,r_Lreq,ierr)
       call MPI_Isend(s_buffL_hi,addL,MPI_REAL8,Lrank,myrank,comm,s_Lreq,ierr)
    endif

    if(Urank.ne.-1)then
       call MPI_Wait(r_Ureq,r_Usta,ierr)
       call MPI_Wait(s_Ureq,s_Usta,ierr)
       do i=1,addL
!          y(i)=(y(i)+r_buffU(i))*0.5d0
          vecU_hi(i)=r_buffU_hi(i)
       enddo
    endif
    if(Lrank.ne.-1)then
       call MPI_Wait(r_Lreq,r_Lsta,ierr)
       call MPI_Wait(s_Lreq,s_Lsta,ierr)
       do i=1,addL
!          y(n-addL+i)=(y(n-addL+i)+r_buffL(i))*0.5d0
          vecL_hi(i)=r_buffL_hi(i)
       enddo
    endif
! -----
    if(Urank.ne.-1)then
       i0=addL
    else
       i0=0
    endif
    if(Urank.ne.-1)then
       do i=1,addL
          s_buffU_lo(i)=x_lo(i)
       enddo
       r_buffU_lo=0.0d0
       call MPI_Irecv(r_buffU_lo,addL,MPI_REAL8,Urank,Urank ,comm,r_Ureq,ierr)
       call MPI_Isend(s_buffU_lo,addL,MPI_REAL8,Urank,myrank,comm,s_Ureq,ierr)
    endif

    if(Lrank.ne.-1)then
       do i=1,addL
          s_buffL_lo(i)=x_lo(i0+n+i)
       enddo
       r_buffL_lo=0.0d0
       call MPI_Irecv(r_buffL_lo,addL,MPI_REAL8,Lrank,Lrank ,comm,r_Lreq,ierr)
       call MPI_Isend(s_buffL_lo,addL,MPI_REAL8,Lrank,myrank,comm,s_Lreq,ierr)
    endif

    if(Urank.ne.-1)then
       call MPI_Wait(r_Ureq,r_Usta,ierr)
       call MPI_Wait(s_Ureq,s_Usta,ierr)
       do i=1,addL
!          y(i)=(y(i)+r_buffU(i))*0.5d0
          vecU_lo(i)=r_buffU_lo(i)
       enddo
    endif
    if(Lrank.ne.-1)then
       call MPI_Wait(r_Lreq,r_Lsta,ierr)
       call MPI_Wait(s_Lreq,s_Lsta,ierr)
       do i=1,addL
!          y(n-addL+i)=(y(n-addL+i)+r_buffL(i))*0.5d0
          vecL_lo(i)=r_buffL_lo(i)
       enddo
    endif

  end subroutine get_ddAS_vec

end module ddvec_kernel

