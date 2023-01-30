! todo dot,dot_localの関係を整理する


! vec_kernelは将来的に廃止してvec_kernel_cpu,vec_kernel_gpuに分離する
! 現バージョンは互換性のために残す
module vec_kernel
  use PARCEL_integer
  implicit none
contains

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

  ! z(:)=alpha*x(:)+beta*y(:)+gamma*z(:)
  subroutine axpy_vec3(n,alpha,beta,gamma,x,y,z)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 alpha,beta,gamma
    real*8 x(n)
    real*8 y(n)
    real*8 z(n)
    integer(kind=intType) i 
    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
       z(i)=alpha*x(i)+beta*y(i)+gamma*z(i)
    enddo
    !$OMP END PARALLEL DO 
  end subroutine axpy_vec3

 ! z(:)=x(:)*y(:)
  subroutine vec_scaling(n,x,y,z)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 alpha,beta
    real*8 x(n)
    real*8 y(n)
    real*8 z(n)
    integer(kind=intType) i 

    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
       z(i)=y(i)*x(i)
    enddo
    !$OMP END PARALLEL DO 
  end subroutine vec_scaling

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

    integer(kind=intType) i 

    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
       y(i) = x(i)
    enddo
    !$OMP END PARALLEL DO 

    ! alpha=1.0d0
    ! beta=0.0d0
    ! call axpy(n,alpha,beta,x,y)

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

#if 0
  subroutine calc_complex_dot(n,comm,x_re,x_im,y_re,y_im,dotp_re,dotp_im)
    use mpi
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer comm
    real*8 x_re(n)
    real*8 x_im(n)
    real*8 y_re(n)
    real*8 y_im(n)
    real*8 dotp_re
    real*8 dotp_im
    real*8 dotp_re_l
    real*8 dotp_im_l
    integer(kind=intType) i 
    integer ierr
    integer ncount
    complex*16 sbuf,rbuf
    dotp_re   = 0.0d0
    dotp_im   = 0.0d0
    dotp_re_l = 0.0d0
    dotp_im_l = 0.0d0
    call calc_complex_dot_local(n,x_re,x_im,y_re,y_im,dotp_re_l,dotp_im_l)
    ncount=1
!    call MPI_Allreduce(dotp_l,dotp,ncount,MPI_REAL8,MPI_SUM,comm,ierr)    
    sbuf = (dotp_re_l,dotp_im_l)
    rbuf = 0.0d0
    call MPI_Allreduce(sbuf,rbuf,ncount,MPI_COMPLEX16,MPI_SUM,comm,ierr)   
    dotp_re = real(rbuf)
    dotp_im = aimag(rbuf)
  end subroutine calc_complex_dot

  ! dotp=<x,y>
  subroutine calc_complex_dot_local(n,x_re,x_im,y_re,y_im,dotp_re_l,dotp_im_l)
    use mpi
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 x_re(n)
    real*8 x_im(n)
    real*8 y_re(n)
    real*8 y_im(n)
    real*8 dotp_re_l
    real*8 dotp_im_l
    integer(kind=intType) i 
    integer ierr
    dotp_re_l=0.0d0
    dotp_im_l=0.0d0
    !$OMP PARALLEL DO PRIVATE(i) REDUCTION(+:dotp_re_l,dotp_im_l) 
    do i=1,n
!       dotp_l = dotp_l+x(i)*y(i)
       dotp_re_l = dotp_re_l + x_re(i) * y_re(i) + x_im(i) * y_im(i)
       dotp_im_l = dotp_im_l - x_re(i) * y_im(i) + x_im(i) * y_re(i)
    enddo
    !$OMP END PARALLEL DO 
  end subroutine calc_complex_dot_local
#endif

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
  subroutine calc_gemm_dot_local_CPU(&
       n,s,t,&
       x, ldx,&
       y, ldy,&
       dotp_l ,lddot, &
       nvecBlock)
    use mpi
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s,t
    integer(kind=intType) ldx
    integer(kind=intType) ldy
    integer(kind=intType) lddot
    real*8 x(ldx*s) ! (n*s)
    real*8 y(ldy*t) ! (n*t)
    real*8 dotp_l(lddot*s) ! (s*t)
    integer nvecBlock
    integer(kind=intType) ii,i,j,k

#if 1
    real*8 tmp(s*t)

    do k=1,t
       do j=1,s
          tmp(j + s*(k-1)) = 0.0d0
       enddo
    enddo

    !    write(*,*)"s,t = ",s,t
    !$OMP PARALLEL DO PRIVATE(ii,k,j,i) REDUCTION(+:tmp) 
    do ii=1,n,nvecBlock
       do k=1,t
          do j=1,s
             do i=ii,min(n,ii+nvecBlock-1)
                tmp(j + s*(k-1)) = &
                     tmp(j + s*(k-1)) + x(i+ldx*(j-1))*y(i+ldy*(k-1))
             enddo
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO 
    do k=1,t
       do j=1,s
          dotp_l(j+lddot*(k-1)) = tmp(j + s*(k-1))
       enddo
    enddo


#else
    do k=1,t
       do j=1,s
          dotp_l(j+lddot*(k-1)) = 0.0d0
       enddo
    enddo

    !    write(*,*)"s,t = ",s,t
    !$OMP PARALLEL DO PRIVATE(ii,k,j,i) REDUCTION(+:dotp_l) 
    do ii=1,n,nvecBlock
       do k=1,t
          do j=1,s
             do i=ii,min(n,ii+nvecBlock-1)
                dotp_l(j+lddot*(k-1)) = &
                     dotp_l(j+lddot*(k-1)) + x(i+ldx*(j-1))*y(i+ldy*(k-1))
             enddo
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO 
#endif

  end subroutine calc_gemm_dot_local_CPU

  ! dotpij=<xi,yj>
  subroutine calc_gemm_dot_CPU(&
       n,s,t,&
       comm,&
       x, ldx,&
       y, ldy,&
       dotp, lddot,&
       nvecBlock)
    use mpi
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer comm
    integer(kind=intType) s,t
    integer(kind=intType) ldx
    real*8 x(ldx*s) ! (n*s)
    real*8 y(ldx*t) ! (n*t)
    integer(kind=intType) ldy
    real*8 dotp(s*t) ! (s*t)
!    real*8 dotp(s*t) ! (s*t)
    integer(kind=intType) lddot
    integer nvecBlock

    real*8 dotp_l(s*t)
    integer(kind=intType) ii,i,j,k
    integer ierr
    integer st

    dotp = 0.0d0
    dotp_l = 0.0d0
#if 0
    dotp = 1.0d0
    return 
#endif
    call calc_gemm_dot_local_CPU(&
         n,s,t,&
         x, ldx,&
         y, ldy,&
         dotp_l, lddot,&
         nvecBlock)
!    write(*,*)"s,t = ",s,t
    st=s*t

    call MPI_Allreduce(dotp_l,dotp,st,MPI_REAL8,MPI_SUM,comm,ierr)   
!    dotp=1.0d0
  end subroutine calc_gemm_dot_CPU

  subroutine gemv(&
       n,m,&
       v , ldv, &
       v0, ldv0,&
       y,nvecBlock)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) m

    integer(kind=intType) ldv
    integer(kind=intType) ldv0

!    real*8 v(n*m)
    real*8 v(ldv*m)
    real*8 v0(ldv0)
    real*8 y(m)
    integer nvecBlock

    integer(kind=intType) i,j,jj

    !$OMP PARALLEL DO PRIVATE(jj,j,i) 
    do jj=1,n,nvecBlock 
       do i=1,m
          do j=jj,min(n,jj+nvecBlock-1)
             v0(j)=v0(j)+y(i)*v(j+ldv*(i-1))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO 

  end subroutine gemv

  ! xのデータが破壊される --> xのデータが壊れないように修正が必要
  ! y = y + x*alpha
  subroutine gemm( &
       x,ldx, &
       y,ldy, &
       alpha,ldalpha, &
       n,m,l, &
       nvecBlock)
    use PARCEL_integer
    implicit none

    integer(kind=intType) ldx
    integer(kind=intType) ldy
    integer(kind=intType) m
    integer(kind=intType) l
    integer(kind=intType) ldalpha
    ! real*8 x(ldx*(m+1))
    ! real*8 y(ldy*(m+1))
    ! real*8 x(*)
    ! real*8 y(*)
    real*8 x(ldx*(l+1))
    real*8 y(ldy*(m+1))

    real*8 alpha(ldalpha*m)

    integer(kind=intType) n
    integer nvecBlock

    integer(kind=intType) i,j,k,jj

    !$OMP PARALLEL DO PRIVATE(jj,k,j,i) 
    do jj=1,n,nvecBlock 

       ! do k=1,m
       !    do j=jj,min(n,jj+nvecBlock-1)
       !       z(j+(k-1)*ldz) = 0.0d0
       !    enddo
       ! enddo

       do k=1,m
!          do i=1,m
          do i=1,l
             ! v0(j)=v0(j)+y(i)*v(j+n*(i-1))
             do j=jj,min(n,jj+nvecBlock-1)
                ! z(j+(k-1)*ldz) = &
                !      z(j+(k-1)*ldz) + &
                !      x(j+(k-1)*ldx) + y(j+(i-1)*ldy)*alpha(i+(k-1)*ldalpha)
                y(j+(k-1)*ldy) = &
                     y(j+(k-1)*ldy) + &
                     x(j+(i-1)*ldx)*alpha(i+(k-1)*ldalpha)

             enddo
          enddo
       enddo

    enddo

    !$OMP END PARALLEL DO 

  end subroutine gemm


  ! dotpij=<xi,yj>
  subroutine calc_sym_gemm_dot_local(n,s,x,y,dotp_l,nvecBlock)
    use mpi
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    real*8 x(n*s)
    real*8 y(n*s)
    real*8 dotp_l(s*s)
    integer nvecBlock

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
  subroutine calc_syrk_dot_local(n,s,x,B,nvecBlock)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    real*8 x(n*s)
    real*8 B(s*s)
    integer nvecBlock

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
       write(*,*)"view_vec",i,x(i)
    enddo
  end subroutine view_vec

  subroutine view_vec_diff(n,x,y)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 x(n)
    real*8 y(n)
    integer(kind=intType) i
    do i=1,n
       write(*,*)"view_vec",i,x(i),y(i)
    enddo
  end subroutine view_vec_diff

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

module vec_ddm_kernel
  use PARCEL_integer
  implicit none
contains

  ! y(:)=alpha*x(:)+beta*y(:)
  subroutine axpy_ddm(n,alpha,beta,x,y,prm)
    use omp_lib
    use PARCEL_integer
    use ddmMat_wrap
    use common_kernel,only : set_region_for_threads

    implicit none
    integer(kind=intType) n
    real*8 alpha,beta
    real*8 x(n)
    real*8 y(n)
    integer(kind=intType) i 
    type(ddm_param), intent(in) :: prm

    integer j, jx, jy, jz
    integer nx, ny, nz, nxy, mx, my, mxy
    integer stmx, stmy, stmz
    integer num_threads, my_thread
    integer nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th
    integer div_direc_th

    nx=prm%nx
    ny=prm%ny
    nz=prm%nz
    nxy=prm%nxy
    mx=prm%mx
    my=prm%my
    mxy=prm%mxy
    stmx=prm%stmx
    stmy=prm%stmy
    stmz=prm%stmz
    div_direc_th = prm%div_direc_th
    !$OMP Parallel default(none)&
    !$OMP firstprivate(  alpha,beta )&
    !$OMP private(my_thread, num_threads)&
    !$OMP private( jz, jy, jx, j)&
    !$OMP private( nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th )&
    !$OMP shared(  nx,ny,nz,nxy,mx,my,mxy,stmx,stmy,stmz,prm)&
    !$OMP shared(  div_direc_th,x,y)
    num_threads = omp_get_num_threads();
    my_thread = omp_get_thread_num();
    call set_region_for_threads( my_thread, num_threads, nx, ny, nz, &
         nxs_th, nxe_th, nys_th, nye_th, nzs_th, &
         nze_th, div_direc_th )


    Do jz=nzs_th, nze_th -1
       Do jy=nys_th, nye_th -1
          Do jx=nxs_th, nxe_th -1
             j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1;
             y(j)= alpha*x(j)+beta*y(j)! v1(j)+a*v2(j)+b*v3(j)

          End do
       End Do
    End Do

    !$OMP end PARALLEL

  end subroutine axpy_ddm

  ! y(:)=x(:)
  subroutine xcpy_ddm(n,x,y,prm)
    use ddmMat_wrap
    use PARCEL_integer
    implicit none
    type(ddm_param):: prm
    integer(kind=intType) n
    real*8 x(n)
    real*8 y(n)
    real*8 alpha,beta
    alpha=1.0d0
    beta=0.0d0
    call axpy_ddm(n,alpha,beta,x,y,prm)
  end subroutine xcpy_ddm

  subroutine set_scala_ddm_CPU(n,alpha,x,prm)
    use omp_lib
    use PARCEL_integer
    use ddmMat_wrap
    use common_kernel,only : set_region_for_threads

    implicit none
    integer(kind=intType) n
    real*8 alpha
    real*8 x(n)
    integer(kind=intType) i 
    type(ddm_param), intent(in) :: prm

    integer j, jx, jy, jz
    integer nx, ny, nz, nxy, mx, my, mxy
    integer stmx, stmy, stmz
    integer num_threads, my_thread
    integer nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th
    integer div_direc_th

    nx=prm%nx
    ny=prm%ny
    nz=prm%nz
    nxy=prm%nxy
    mx=prm%mx
    my=prm%my
    mxy=prm%mxy
    stmx=prm%stmx
    stmy=prm%stmy
    stmz=prm%stmz
    div_direc_th = prm%div_direc_th
    !$OMP Parallel default(none)&
    !$OMP private(my_thread, num_threads)&
    !$OMP private( jz, jy, jx, j)&
    !$OMP private( nxs_th, nxe_th, nys_th, nye_th, nzs_th, nze_th )&
    !$OMP shared(  nx,ny,nz,nxy,mx,my,mxy,stmx,stmy,stmz,prm)&
    !$OMP shared(  div_direc_th,x,alpha)
    num_threads = omp_get_num_threads();
    my_thread = omp_get_thread_num();
    call set_region_for_threads( my_thread, num_threads, nx, ny, nz, &
         nxs_th, nxe_th, nys_th, nye_th, nzs_th, &
         nze_th, div_direc_th )

    Do jz=nzs_th, nze_th -1
       Do jy=nys_th, nye_th -1
          Do jx=nxs_th, nxe_th -1
             j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) +1;
             x(j)= alpha
          End do
       End Do
    End Do

    !$OMP end PARALLEL

  end subroutine set_scala_ddm_CPU

end module vec_ddm_kernel

module vec_kernel_CPU
  use PARCEL_integer
  implicit none
contains
  ! y(:)=alpha*x(:)+beta*y(:)
  subroutine axpy_FP64_CPU(n,alpha,beta,x,y)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 alpha
    real*8 beta
    real*8 x(n)
    real*8 y(n)
    integer(kind=intType) i 

    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
       y(i)=alpha*x(i)+beta*y(i)
    enddo
    !$OMP END PARALLEL DO 

  end subroutine axpy_FP64_CPU

 ! y(:)=alpha*x(:)+beta*y(:)+gamma*z(:)
  subroutine axpy_vec3_FP64_CPU(n,alpha,beta,gamma,x,y,z)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 alpha,beta,gamma
    real*8 x(n)
    real*8 y(n)
    real*8 z(n)
    integer(kind=intType) i 
    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
       z(i)=alpha*x(i)+beta*y(i)+gamma*z(i)
    enddo
    !$OMP END PARALLEL DO 
  end subroutine axpy_vec3_FP64_CPU

 ! y(:)=alpha*x(:)+beta*y(:)+gamma*z(:)
  subroutine axpy_vec4_FP64_CPU(&
       n,alpha,beta,gamma,delta,&
       x,y,z,w)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 alpha,beta,gamma,delta
    real*8 x(n)
    real*8 y(n)
    real*8 z(n)
    real*8 w(n)
    integer(kind=intType) i 
    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
       w(i)=alpha*x(i)+beta*y(i)+gamma*z(i)+delta*w(i)
    enddo
    !$OMP END PARALLEL DO 
  end subroutine axpy_vec4_FP64_CPU

  ! dotp=<x,y>
  subroutine calc_dot_local_FP64_CPU(n,x,y,dotp_l)
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
  end subroutine calc_dot_local_FP64_CPU

  ! dotp=<x,y>
  subroutine calc_dot_FP64_CPU(n,comm,x,y,dotp)
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
!    call calc_dot_local(n,x,y,dotp_l)
    call calc_dot_local_FP64_CPU(n,x,y,dotp_l)
    ncount=1
    call MPI_Allreduce(dotp_l,dotp,ncount,MPI_REAL8,MPI_SUM,comm,ierr)    

  end subroutine calc_dot_FP64_CPU

  ! z(:)=x(:)*y(:)
  subroutine vec_scaling_FP64_CPU(n,x,y,z)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 alpha,beta
    real*8 x(n)
    real*8 y(n)
    real*8 z(n)
    integer(kind=intType) i 
    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
       z(i)=y(i)*x(i)
    enddo
    !$OMP END PARALLEL DO 
  end subroutine vec_scaling_FP64_CPU

#if 0
  subroutine set_AS_vec_FP64_CPU( &
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
  end subroutine set_AS_vec_FP64_CPU
#endif

end module vec_kernel_CPU

  subroutine get_EigenPair_count_KSP(ncount,method)
    use krylov
    implicit none
    type(KSP) :: method
    integer ncount
    if( allocated(method%EigenValueReal) )then
       ncount = size( method%EigenValueReal )
    else
       ncount = -1
    endif

  end subroutine get_EigenPair_count_KSP

  subroutine get_EigenPair_KSP( &
       ncount,&
       EigenValueReal,EigenValueImag, &
       EigenVector, &
       EigenErrorNorm, &
       EmaxID,EminID,method)
    use krylov
    use vec_kernel
    implicit none
    integer ncount
    type(KSP) :: method
    real*8 EigenVector(method%prm%n*ncount) 
    real*8 EigenValueReal(ncount) 
    real*8 EigenValueImag(ncount) 
    real*8 EigenErrorNorm(ncount) 
    integer EmaxID,EminID
    integer i

    do i = 1,ncount
       EigenValueReal(i) = method%EigenValueReal(i)
       EigenValueImag(i) = method%EigenValueImag(i)
       EigenErrorNorm(i) = method%EigenErrorNorm(i)
    enddo
    call xcpy(method%prm%n*ncount,method%EigenVector,EigenVector)

    EmaxID = method%MaxEigenID
    EminID = method%MinEigenID

  end subroutine get_EigenPair_KSP
