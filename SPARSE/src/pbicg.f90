module KSP_BiCG
contains
  !   e=r-alpha*q
  subroutine BiCG_axpy3(n,r,q,e,alpha)
  use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 r(n)
    real*8 q(n)
    real*8 e(n)
    real*8 alpha
    integer(kind=intType) i 
    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
       e(i)=r(i)-alpha*q(i)
    enddo
    !$OMP END PARALLEL DO 
  end subroutine BiCG_axpy3

  !   p=r+beta*(p-c3*q)
  subroutine BiCG_axpy4_(n,p,r,q,beta,c3)
  use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 p(n)
    real*8 r(n)
    real*8 q(n) 
    real*8 beta
    real*8 c3

    integer(kind=intType) i
    !$OMP PARALLEL DO PRIVATE(i)
    do i=1,n
       p(i)=r(i)+beta*( p(i) -c3*q(i) )
    enddo
    !$OMP END PARALLEL DO 

  end subroutine BiCG_axpy4_

  !   r=e-c3*v  
  !   c1  =<r0,r>   
  !   rdot=<r,r>
  subroutine BiCG_axpy3_dot(n,icomm,r,e,v,r0,c3,c1,rdot)
    use mpi
  use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer icomm
    real*8 r(n)
    real*8 e(n)
    real*8 v(n) 
    real*8 r0(n) 
    real*8 c3
    real*8 c1
    real*8 rdot
    real*8 tmp1,tmp2
    real*8 sbuf(2),rbuf(2)
    integer(kind=intType) i
    integer ierr

    tmp1=0.0d0
    tmp2=0.0d0
    !$OMP PARALLEL DO PRIVATE(i) REDUCTION(+:tmp1,tmp2) 
    do i=1,n
       r(i)=e(i)-c3*v(i)
       tmp1=tmp1+r(i)*r0(i)
       tmp2=tmp2+r(i)*r(i)
    enddo
    !$OMP END PARALLEL DO 

    sbuf(1)=tmp1
    sbuf(2)=tmp2
    rbuf(1)=0.0d0
    rbuf(2)=0.0d0
    call MPI_Allreduce(sbuf,rbuf,2,MPI_REAL8,MPI_SUM,icomm,ierr)    
    c1=rbuf(1)
    rdot=rbuf(2)


  end subroutine BiCG_axpy3_dot

  !   x=x+alpha*pp+c3*pe   
  subroutine BiCG_axpy4(n,x,pp,pe,alpha,c3)
  use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 x(n)
    real*8 pp(n)
    real*8 pe(n) 
    real*8 alpha
    real*8 c3

    integer(kind=intType) i
    !$OMP PARALLEL DO PRIVATE(i)
    do i=1,n
       x(i)=x(i)+alpha*pp(i)+c3*pe(i)
    enddo
    !$OMP END PARALLEL DO 

  end subroutine BiCG_axpy4

  ! ev_dot = <e,v>, vv_dot=<v,v> 
  subroutine BiCG_calc_dot2(n,icomm,e,v,ev_dot,vv_dot)
    use mpi
  use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer icomm
    real*8 e(n)
    real*8 v(n)
    real*8 ev_dot
    real*8 vv_dot

    integer(kind=intType) i
    real*8 tmp1,tmp2
    real*8 sbuf(2),rbuf(2)
    integer ierr

    tmp1=0.0d0
    tmp2=0.0d0
    !$OMP PARALLEL DO PRIVATE(i) REDUCTION(+:tmp1,tmp2) 
    do i=1,n
       tmp1=tmp1+e(i)*v(i)
       tmp2=tmp2+v(i)*v(i)
    enddo
    !$OMP END PARALLEL DO 

    sbuf(1)=tmp1
    sbuf(2)=tmp2
    rbuf(1)=0.0d0
    rbuf(2)=0.0d0
    call MPI_Allreduce(sbuf,rbuf,2,MPI_REAL8,MPI_SUM,icomm,ierr)    
    ev_dot=rbuf(1)
    vv_dot=rbuf(2)

  end subroutine BiCG_calc_dot2

end module KSP_BiCG
subroutine KSP_pBiCG( pbicg &
     )
  use krylov
  use krylov_kernel
  use vec_kernel
  use matrix_kernel
  use KSP_BiCG

  use PARCEL_integer
  implicit none
  type(KSP) :: pbicg
  integer(kind=intType) n
  integer icomm
  integer itrmax ! 最大反復回数
  real*8 rtolmax ! 相対残差 (収束条件)
  real*8 abstolmax  !  絶対残差 (収束条件)
  integer myrank

  real*8,allocatable :: r0(:)
  real*8,allocatable :: p(:)
  real*8,allocatable :: q(:)
  real*8,allocatable :: pe(:)
  real*8,allocatable :: v(:)
  real*8,allocatable :: pp(:)
  real*8,allocatable :: e(:)

  real*8,allocatable :: x(:)
  real*8,allocatable :: r(:)

  real*8 alpha,beta
  real*8 c1,c2,c3
  real*8 rnorm
  real*8 bnorm
  integer itr
  real*8 tmp1,tmp2

  real*8 stime0
  real*8 stime,etime
  real*8 bicg_time
  real*8 precon_time
  real*8 matvec_time
  bicg_time  =0.0d0
  precon_time=0.0d0
  matvec_time=0.0d0

  stime0=MPI_Wtime()

  itrmax=pbicg%itrmax
  rtolmax=pbicg%rtolmax
  abstolmax=pbicg%abstolmax
  myrank=pbicg%commVec%myrank
  icomm=pbicg%commVec%comm
  n=pbicg%prm%n

  allocate(x(n))
  allocate(r(n))
  allocate(r0(n))
  allocate(p(n))
  allocate(q(n))
  allocate(pe(n))
  allocate(v(n))
  allocate(pp(n))
  allocate(e(n))

  call zero_initialize(n,r)
  call zero_initialize(n,x)
  call zero_initialize(n,r0)
  call zero_initialize(n,p)
  call zero_initialize(n,q)
  call zero_initialize(n,pe)
  call zero_initialize(n,v)
  call zero_initialize(n,pp)
  call zero_initialize(n,e)

  call xcpy(n,pbicg%x,x)
!  call write_vec(n,pbicg%Dinv)
! -----
! debug
  ! call calc_Matvec(pbicg,x,p)
  ! return
! -----

  ! ------
  ! pbicg start
  ! r=b-Ax
  ! r0=r
  ! p=r
  ! c1=<r0,r0>
  ! do itr =1,....
  !   solve Mpp=p
  !   q=App
  !   c2=<r0,q>
  !   alpha=c1/c2
  !   e=r-alpha*q
  !   solve Mpe=e
  !   v=Ape
  !   c3=<e,v>/<v,v>
  !   x=x+alpha*pp+c3*pe
  !   r=e-c3*v
  !   c1=<r0,r>
  !   if(||r||<rmax) exit
  !   beta=c1/(c3*c2)
  !   p=r+beta*(p-c3*q)
  ! enddo
  ! ------
  c1=0.0d0
  c2=0.0d0
  c3=0.0d0

  ! r=b-Ax
  call calc_res(pbicg,x,pbicg%b,r)

  ! r0=r
  call xcpy(n,r,r0)

  ! p=r
  call xcpy(n,r,p)
  call calc_dot(n,icomm,pbicg%b,pbicg%b,bnorm)

  ! c1=<r0,r0>
  call calc_dot(n,icomm,r0,r0,c1)
  rnorm=c1
  bnorm=sqrt(bnorm)
  rnorm=sqrt(rnorm)

  do itr = 1,itrmax
     pbicg%reshistory(itr)=rnorm/bnorm 
     if(myrank.eq.0)then
!        write(*,*)"@ res @",itr,rnorm,rnorm/bnorm
     endif
     !   solve Mpp=p
     stime=MPI_Wtime()
     call  calc_precon(pbicg,p,pp) 
     etime=MPI_Wtime()
     precon_time=precon_time+etime-stime

     !   q=App
     stime=MPI_Wtime()
     call calc_Matvec(pbicg,pp,q)
     etime=MPI_Wtime()
     matvec_time=matvec_time+etime-stime

     !   c2=<r0,q>
     call calc_dot(n,icomm,r0,q,c2)
     !   alpha=c1/c2 
     alpha=c1/c2

     !   e=r-alpha*q
     call BiCG_axpy3(n,r,q,e,alpha)

     !   solve Mpe=e
     stime=MPI_Wtime()
     call  calc_precon(pbicg,e,pe)
     etime=MPI_Wtime()
     precon_time=precon_time+etime-stime

     !   v=Ape
     stime=MPI_Wtime()
     call calc_Matvec(pbicg,pe,v)
     etime=MPI_Wtime()
     matvec_time=matvec_time+etime-stime

     !   c3=<e,v>/<v,v>
     call BiCG_calc_dot2(n,icomm,e,v,tmp1,tmp2)
     c3=tmp1/tmp2

     !   x=x+alpha*pp+c3*pe
     call BiCG_axpy4(n,x,pp,pe,alpha,c3)

     !   r=e-c3*v
     !   c1=<r0,r>
     call BiCG_axpy3_dot(n,icomm,r,e,v,r0,c3,c1,rnorm)

     rnorm=sqrt(rnorm)
     if( (rnorm/bnorm<rtolmax) ) then
!        pbicg%reshistory(itr+1)=rnorm/bnorm 
        if(myrank.eq.0)then
!           write(*,*)"@ res @",itr+1,rnorm,rnorm/bnorm
        endif
        exit
     endif

     !   beta=c1/(c3*c2)
     beta=c1/(c3*c2)

     !   p=r+beta*(p-c3*q)
     call BiCG_axpy4_(n,p,r,q,beta,c3)

  enddo
#if 0
  ! ------
  ! 真の残差を確認
  call calc_res(pbicg,x,pbicg%b,r)
  call calc_dot(n,icomm,r,r,c1)
  rnorm=sqrt(c1)
  if(myrank.eq.0)then
     write(*,*)"@ true res @",rnorm,rnorm/bnorm
  endif
#endif
  call xcpy(n,x,pbicg%x)

  pbicg%niter=min(itrmax,itr)
  if(itrmax<=itr)then
     pbicg%iret=PARCEL_RETURNCODE_MAXITERATION
  else
     pbicg%iret=PARCEL_RETURNCODE_NORMAL
  endif
  ! ------
  deallocate(r0)
  deallocate(p)
  deallocate(q)
  deallocate(pe)
  deallocate(v)
  deallocate(pp)
  deallocate(e)
  deallocate(x)
  deallocate(r)

  etime=MPI_Wtime()
  bicg_time=etime-stime0
#if 0
  if(myrank.eq.0)then
     write(*,*)"@ bicg     [s] @",bicg_time
     write(*,*)"@   matvec [s] @",matvec_time
     write(*,*)"@   precon [s] @",precon_time
  endif
#endif

end subroutine KSP_pBiCG
