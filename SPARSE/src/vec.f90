! todo dot,dot_localの関係を整理する
module vec_kernel
  use PARCEL_integer
  implicit none
  integer(kind=intType),parameter :: nvecBlock = 100
contains
! 
  ! x(:)=0.0d0
  subroutine zero_initialize(n,x)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 x(n)
    integer(kind=intType) i 
    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
       x(i)=0.0d0
    enddo
    !$OMP END PARALLEL DO 
  end subroutine zero_initialize

  subroutine zero_initialize_int(n,x)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) x(n)
    integer(kind=intType) i 
    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
       x(i)=0
    enddo
    !$OMP END PARALLEL DO 
  end subroutine zero_initialize_int

  subroutine zero_initialize_int4(n,x)
    use PARCEL_integer
    implicit none
    integer n
    integer x(n)
    integer i 
    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
       x(i)=0
    enddo
    !$OMP END PARALLEL DO 
  end subroutine zero_initialize_int4

  ! y(:)=alpha*x(:)+beta*y(:)
  subroutine axpy(n,alpha,beta,x,y)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 alpha,beta
    real*8 x(n)
    real*8 y(n)
    integer(kind=intType) i 
    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
       y(i)=alpha*x(i)+beta*y(i)
    enddo
    !$OMP END PARALLEL DO 
  end subroutine axpy

  ! y(:)=alpha*x(:)
  subroutine inplace_vec(n,alpha,x,y)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 alpha
    real*8 x(n)
    real*8 y(n)
    integer(kind=intType) i 
    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
       y(i)=alpha*x(i)
    enddo
    !$OMP END PARALLEL DO 
  end subroutine inplace_vec


  ! y(:)=x(:)
  subroutine xcpy(n,x,y)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 x(n)
    real*8 y(n)
    real*8 alpha,beta
    alpha=1.0d0
    beta=0.0d0
    call axpy(n,alpha,beta,x,y)
  end subroutine xcpy

  ! dotp=<x,y>
  subroutine calc_dot(n,comm,x,y,dotp)
    use mpi
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer comm
    real*8 x(n)
    real*8 y(n)
    real*8 dotp
    real*8 dotp_l
    integer(kind=intType) i 
    integer ierr
    integer ncount
    dotp=0.0d0
    dotp_l=0.0d0
    call calc_dot_local(n,x,y,dotp_l)
    ncount=1
    call MPI_Allreduce(dotp_l,dotp,ncount,MPI_REAL8,MPI_SUM,comm,ierr)    

  end subroutine calc_dot

  ! dotp=<x,y>
  subroutine calc_dot_local(n,x,y,dotp_l)
    use mpi
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 x(n)
    real*8 y(n)
    real*8 dotp_l
    integer(kind=intType) i 
    integer ierr
    dotp_l=0.0d0
    !$OMP PARALLEL DO PRIVATE(i) REDUCTION(+:dotp_l) 
    do i=1,n
       dotp_l = dotp_l+x(i)*y(i)
    enddo
    !$OMP END PARALLEL DO 
  end subroutine calc_dot_local

  ! dotpij=<xi,yj>
  subroutine calc_gemm_dot(n,s,t,comm,x,y,dotp)
    use mpi
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer comm
    integer(kind=intType) s,t
    real*8 x(n*s)
    real*8 y(n*t)
    real*8 dotp(s*t)
    real*8 dotp_l(s*t)
    integer(kind=intType) ii,i,j,k
    integer ierr
    integer st
    dotp=0.0d0
    dotp_l=0.0d0
    call calc_gemm_dot_local(n,s,t,x,y,dotp_l)
    st=s*t
    call MPI_Allreduce(dotp_l,dotp,st,MPI_REAL8,MPI_SUM,comm,ierr)   

  end subroutine calc_gemm_dot

  subroutine gemv(n,m,v,v0,y)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) m

    real*8 v(n*m)
    real*8 v0(n)
    real*8 y(m)

    integer(kind=intType) i,j,jj

    !$OMP PARALLEL DO PRIVATE(jj,j,i) 
    do jj=1,n,nvecBlock 
       do i=1,m
          do j=jj,min(n,jj+nvecBlock-1)
             v0(j)=v0(j)+y(i)*v(j+n*(i-1))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO 

  end subroutine gemv

  ! w = v - w*z
  subroutine gemm(n,m,v,w,z)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) m

    real*8 v(n*(m+1))
    real*8 w(n*(m+1))
    real*8 z(m*m)

    integer(kind=intType) i,j,k,jj

    !$OMP PARALLEL DO PRIVATE(jj,k,j,i) 
    do jj=1,n,nvecBlock 

       do k=1,m
          do i=1,m
             do j=jj,min(n,jj+nvecBlock-1)
                !                v0(j)=v0(j)+y(i)*v(j+n*(i-1))
                v(j+(k-1)*n)=v(j+(k-1)*n) - w(j+(i-1)*n)*z(i+(k-1)*m)
             enddo
          enddo
       enddo
       do k=1,m
          do j=jj,min(n,jj+nvecBlock-1)
             !                v0(j)=v0(j)+y(i)*v(j+n*(i-1))
             w(j+(k-1)*n)=v(j+(k-1)*n)
          enddo
       enddo

    enddo

    !$OMP END PARALLEL DO 

  end subroutine gemm


  ! dotpij=<xi,yj>
  subroutine calc_gemm_dot_local(n,s,t,x,y,dotp_l)
    use mpi
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s,t
    real*8 x(n*s)
    real*8 y(n*t)
    real*8 dotp_l(s*t)
    integer(kind=intType) ii,i,j,k
    integer ierr
    dotp_l=0.0d0

    !$OMP PARALLEL DO PRIVATE(ii,k,j,i) REDUCTION(+:dotp_l) 
    do ii=1,n,nvecBlock
       do k=1,t
          do j=1,s
             do i=ii,min(n,ii+nvecBlock-1)
                dotp_l(j+s*(k-1)) = &
                     dotp_l(j+s*(k-1)) +x(i+n*(j-1))*y(i+n*(k-1))
             enddo
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO 

  end subroutine calc_gemm_dot_local

  ! dotpij=<xi,yj>
  subroutine calc_sym_gemm_dot_local(n,s,x,y,dotp_l)
    use mpi
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    real*8 x(n*s)
    real*8 y(n*s)
    real*8 dotp_l(s*s)
    integer(kind=intType) ii,i,j,k
    integer ierr
    dotp_l=0.0d0

    !$OMP PARALLEL DO PRIVATE(ii,k,j,i) REDUCTION(+:dotp_l) 
    do ii=1,n,nvecBlock
       do k=1,s
          do j=1,k
             do i=ii,min(n,ii+nvecBlock-1)
                dotp_l(j+s*(k-1)) = dotp_l(j+s*(k-1)) +x(i+n*(j-1))*y(i+n*(k-1))
             enddo
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO 

  end subroutine calc_sym_gemm_dot_local

  ! B(i,j)=<xi,xj>
  subroutine calc_syrk_dot_local(n,s,x,B)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    real*8 x(n*s)
    real*8 B(s*s)
    integer(kind=intType) ii,i,j,k
    real*8 Btmp1,Btmp2,Btmp3,Btmp4
    integer(kind=intType) s_odd
    integer(kind=intType) count
    s_odd=mod(s,2)
#if 1
    B=0.0d0
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(ii,k,j,i,Btmp1,Btmp2,Btmp3,Btmp4,count) & 
    !$OMP REDUCTION(+:B) 
    do ii=1,n,nvecBlock
       count=-1
       if(s_odd.eq.0)then
          do k=1,s,2
             do j=s,k+2,-2
                Btmp1=0.0d0
                Btmp2=0.0d0
                Btmp3=0.0d0
                Btmp4=0.0d0
                do i=ii,min(n,ii+nvecBlock-1)
                   Btmp1 = Btmp1 + x(i+n*(k-1))   * x(i+n*(j-1))
                   Btmp2 = Btmp2 + x(i+n*(k-1+1)) * x(i+n*(j-1))       
                   Btmp3 = Btmp3 + x(i+n*(k-1))   * x(i+n*(j-1-1))
                   Btmp4 = Btmp4 + x(i+n*(k-1+1)) * x(i+n*(j-1-1)) 
                enddo
                B(k+(j-1)*s)     = B(k+(j-1)*s)     + Btmp1
                B(k+1+(j-1)*s)   = B(k+1+(j-1)*s)   + Btmp2    
                B(k+(j-1-1)*s)   = B(k+(j-1-1)*s)   + Btmp3
                B(k+1+(j-1-1)*s) = B(k+1+(j-1-1)*s) + Btmp4 
             enddo
             j=k+1
             Btmp1=0.0d0
             Btmp2=0.0d0
             do i=ii,min(n,ii+nvecBlock-1)
                Btmp1 = Btmp1 + x(i+n*(k-1))   * x(i+n*(j-1))
                Btmp2 = Btmp2 + x(i+n*(k-1+1)) * x(i+n*(j-1)) 
             enddo
             B(k+(j-1)*s)   = B(k+(j-1)*s)   + Btmp1
             B(k+1+(j-1)*s) = B(k+1+(j-1)*s) + Btmp2 

             j=k
             Btmp1=0.0d0
             do i=ii,min(n,ii+nvecBlock-1)
                Btmp1 = Btmp1 + x(i+n*(k-1)) * x(i+n*(j-1))
             enddo
             B(k+(j-1)*s)=B(k+(j-1)*s) + Btmp1
          enddo

       else 
          do k=1,s-s_odd,2
             count=count*(-1)
             if(count.eq.(-1))then
                j=k
                Btmp1=0.0d0
                do i=ii,min(n,ii+nvecBlock-1)
                   Btmp1 = Btmp1 + x(i+n*(k-1)) * x(i+n*(j-1))
                enddo
                B(k+(j-1)*s)=B(k+(j-1)*s) + Btmp1

                do j=k+1,s,2
                   Btmp1=0.0d0
                   Btmp2=0.0d0
                   Btmp3=0.0d0
                   Btmp4=0.0d0
                   do i=ii,min(n,ii+nvecBlock-1)
                      Btmp1 = Btmp1 + x(i+n*(k-1))   * x(i+n*(j-1))
                      Btmp2 = Btmp2 + x(i+n*(k-1+1)) * x(i+n*(j-1))       
                      Btmp3 = Btmp3 + x(i+n*(k-1))   * x(i+n*(j-1+1))
                      Btmp4 = Btmp4 + x(i+n*(k-1+1)) * x(i+n*(j-1+1)) 
                   enddo
                   B(k+(j-1)*s)     = B(k+(j-1)*s)     + Btmp1
                   B(k+1+(j-1)*s)   = B(k+1+(j-1)*s)   + Btmp2    
                   B(k+(j+1-1)*s)   = B(k+(j+1-1)*s)   + Btmp3
                   B(k+1+(j+1-1)*s) = B(k+1+(j+1-1)*s) + Btmp4 
                enddo
             else
                do j=s,k+2,-2
                   Btmp1=0.0d0
                   Btmp2=0.0d0
                   Btmp3=0.0d0
                   Btmp4=0.0d0
                   do i=ii,min(n,ii+nvecBlock-1)
                      Btmp1 = Btmp1 + x(i+n*(k-1))   * x(i+n*(j-1))
                      Btmp2 = Btmp2 + x(i+n*(k-1+1)) * x(i+n*(j-1))       
                      Btmp3 = Btmp3 + x(i+n*(k-1))   * x(i+n*(j-1-1))
                      Btmp4 = Btmp4 + x(i+n*(k-1+1)) * x(i+n*(j-1-1)) 
                   enddo
                   B(k+(j-1)*s)     = B(k+(j-1)*s)     + Btmp1
                   B(k+1+(j-1)*s)   = B(k+1+(j-1)*s)   + Btmp2    
                   B(k+(j-1-1)*s)   = B(k+(j-1-1)*s)   + Btmp3
                   B(k+1+(j-1-1)*s) = B(k+1+(j-1-1)*s) + Btmp4 
                enddo
                j=k
                Btmp1=0.0d0
                do i=ii,min(n,ii+nvecBlock-1)
                   Btmp1 = Btmp1 + x(i+n*(k-1)) * x(i+n*(j-1))
                enddo
                B(k+(j-1)*s)=B(k+(j-1)*s) + Btmp1
             endif
          enddo

          k=s
          j=s
          Btmp1=0.0d0
          do i=ii,min(n,ii+nvecBlock-1)
             Btmp1 = Btmp1 + x(i+n*(k-1))*x(i+n*(j-1))
          enddo
          B(k+(j-1)*s)=B(k+(j-1)*s) + Btmp1
       endif

    enddo
    !$OMP END PARALLEL DO 
#else
    B=0.0d0
    !$OMP PARALLEL DO PRIVATE(ii,k,j,i) REDUCTION(+:B) 
    do ii=1,n,nvecBlock

       do k=1,s
          do j=k,s
             do i=ii,min(n,ii+nvecBlock-1)
                B(k+s*(j-1)) = B(k+s*(j-1)) +x(i+n*(j-1))*x(i+n*(k-1))
             enddo
          enddo
       enddo

    enddo
    !$OMP END PARALLEL DO 
#endif

  end subroutine calc_syrk_dot_local

  subroutine view_vec(n,x)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 x(n)
    integer(kind=intType) i
    do i=1,n
       write(*,*)i,x(i)
    enddo
  end subroutine view_vec

  subroutine set_AS_vec( &
       n,comm, &
       addL, &
       myrank,Urank,Lrank, &
       x,y)
    use mpi
    use PARCEL_integer
    implicit none
    integer comm
    integer(kind=intType) n
    integer addL
    integer myrank,Urank,Lrank
    real*8 x(n)
    real*8 y(n+addL*2)

    integer(kind=intType) i 
    real*8 r_buffL(addL),r_buffU(addL)
    real*8 s_buffL(addL),s_buffU(addL)

    integer s_Ureq,s_Lreq
    integer r_Ureq,r_Lreq
    integer s_Usta(MPI_STATUS_SIZE),s_Lsta(MPI_STATUS_SIZE)
    integer r_Usta(MPI_STATUS_SIZE),r_Lsta(MPI_STATUS_SIZE)

    integer ierr
    integer(kind=intType) i0

    if(Urank.ne.-1)then
       do i=1,addL
          s_buffU(i)=x(i)
       enddo
       r_buffU=0.0d0
       call MPI_Irecv(r_buffU,addL,MPI_REAL8,Urank,Urank ,comm,r_Ureq,ierr)
       call MPI_Isend(s_buffU,addL,MPI_REAL8,Urank,myrank,comm,s_Ureq,ierr)
    endif
    if(Lrank.ne.-1)then
       do i=1,addL
          !          x(n-addL+i)=x(n-addL+i)*0.5d0
          s_buffL(i)=x(n-addL+i)
       enddo
       r_buffL=0.0d0
       call MPI_Irecv(r_buffL,addL,MPI_REAL8,Lrank,Lrank ,comm,r_Lreq,ierr)
       call MPI_Isend(s_buffL,addL,MPI_REAL8,Lrank,myrank,comm,s_Lreq,ierr)
    endif

    i0=0
    if(Urank.ne.-1)then
       call MPI_Wait(r_Ureq,r_Usta,ierr)
       call MPI_Wait(s_Ureq,s_Usta,ierr)
       do i=1,addL
          y(i+i0)=r_buffU(i)
       enddo
       i0=i0+addL
    endif
    call xcpy(n,x,y(1+i0))
    i0=i0+n

    if(Lrank.ne.-1)then
       call MPI_Wait(r_Lreq,r_Lsta,ierr)
       call MPI_Wait(s_Lreq,s_Lsta,ierr)
       do i=1,addL
          y(i0+i)=r_buffL(i)
       enddo
    endif
  end subroutine set_AS_vec

  subroutine get_AS_vec( &
       n,comm, &
       addL, &
       myrank,Urank,Lrank, &
       x, &
       vecL,vecU)
    !       x,y)
    use mpi
    use PARCEL_integer
    implicit none
    integer comm
    integer(kind=intType) n
    integer addL
    integer myrank,Urank,Lrank
    real*8 x(n+addL*2)
    real*8 vecL(addL),vecU(addL)

    !    real*8 y(n)
    integer(kind=intType) iflag ! オーバーラップ領域の補間方法

    integer(kind=intType) i 
    real*8 r_buffL(addL),r_buffU(addL)
    real*8 s_buffL(addL),s_buffU(addL)

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
          s_buffU(i)=x(i)
       enddo
       r_buffU=0.0d0
       call MPI_Irecv(r_buffU,addL,MPI_REAL8,Urank,Urank ,comm,r_Ureq,ierr)
       call MPI_Isend(s_buffU,addL,MPI_REAL8,Urank,myrank,comm,s_Ureq,ierr)
    endif

    if(Lrank.ne.-1)then
       do i=1,addL
          s_buffL(i)=x(i0+n+i)
       enddo
       r_buffL=0.0d0
       call MPI_Irecv(r_buffL,addL,MPI_REAL8,Lrank,Lrank ,comm,r_Lreq,ierr)
       call MPI_Isend(s_buffL,addL,MPI_REAL8,Lrank,myrank,comm,s_Lreq,ierr)
    endif

    if(Urank.ne.-1)then
       call MPI_Wait(r_Ureq,r_Usta,ierr)
       call MPI_Wait(s_Ureq,s_Usta,ierr)
       do i=1,addL
          !          y(i)=(y(i)+r_buffU(i))*0.5d0
          vecU(i)=r_buffU(i)
       enddo
    endif
    if(Lrank.ne.-1)then
       call MPI_Wait(r_Lreq,r_Lsta,ierr)
       call MPI_Wait(s_Lreq,s_Lsta,ierr)
       do i=1,addL
          !          y(n-addL+i)=(y(n-addL+i)+r_buffL(i))*0.5d0
          vecL(i)=r_buffL(i)
       enddo
    endif
  end subroutine get_AS_vec

end module vec_kernel

