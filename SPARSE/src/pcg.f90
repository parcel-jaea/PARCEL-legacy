module KSP_CG
contains
  !   x = x + alpha*p
  !   p = s + beta * p
  subroutine CG_axpy2(n,x,p,s,alpha,beta)
  use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 x(n)
    real*8 p(n)
    real*8 s(n)
    real*8 alpha
    real*8 beta
    integer(kind=intType) i

    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
       x(i)=x(i)+alpha*p(i)
       p(i)=s(i)+beta*p(i)
    enddo
    !$OMP END PARALLEL DO 

  end subroutine CG_axpy2
 
  !   rho1 = <s,r>
  !   rdot = <r,r>
  subroutine CG_calc_dot2(n,icomm,s,r,rho1,rdot)
    use mpi
  use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer icomm
    real*8 s(n)
    real*8 r(n)
    real*8 rho1
    real*8 rdot
    real*8 tmp1,tmp2
    real*8 sbuf(2),rbuf(2)

    integer(kind=intType) i
    integer ierr

    tmp1=0.0d0
    tmp2=0.0d0
    !$OMP PARALLEL DO PRIVATE(i) REDUCTION(+:tmp1,tmp2) 
    do i=1,n
       tmp1=tmp1+s(i)*r(i)
       tmp2=tmp2+r(i)*r(i)
    enddo
    !$OMP END PARALLEL DO 

    sbuf(1)=tmp1
    sbuf(2)=tmp2
    rbuf(1)=0.0d0
    rbuf(2)=0.0d0
    call MPI_Allreduce(sbuf,rbuf,2,MPI_REAL8,MPI_SUM,icomm,ierr)    
    rho1=rbuf(1)
    rdot=rbuf(2)

  end subroutine CG_calc_dot2
end module KSP_CG

subroutine KSP_pCG( pcg&
     )
  use KSP_CG
  use krylov
  use krylov_kernel
  use vec_kernel
#if 0
  use mpi
#endif

  use PARCEL_integer
  implicit none
  type(KSP) :: pcg
  integer(kind=intType) n
  integer icomm
  integer itrmax ! 最大反復回数
  real*8 rtolmax ! 相対残差 (収束条件)
  real*8 abstolmax  !  絶対残差 (収束条件)
  integer myrank

  real*8,allocatable :: p(:)
  real*8,allocatable :: q(:)
  real*8,allocatable :: s(:)
  real*8,allocatable :: r(:)
  real*8,allocatable :: x(:)

  real*8 alpha,beta
  real*8 gamma
  real*8 rho0,rho1
  real*8 rnorm
  real*8 bnorm
  integer(kind=intType) itr
  real*8 tmp1,tmp2
  real*8 rtime(40)
  real*8 stime0,etime0
  real*8 stime,etime
  integer ierr

#if 0
  integer i
 real*8 recvbuf(pcg%commVec%npes)
#endif

! write(*,*)"@@@@@@@@@@"

  rtime=0.0d0

  itrmax=pcg%itrmax
  rtolmax=pcg%rtolmax
  abstolmax=pcg%abstolmax
  myrank=pcg%commVec%myrank
  icomm=pcg%commVec%comm
  n=pcg%prm%n

  call MPI_Barrier(icomm,ierr)
  stime0=MPI_Wtime()
  stime=MPI_Wtime()

  allocate(x(n))
  allocate(r(n))
  allocate(p(n))
  allocate(q(n))
  allocate(s(n))

  call zero_initialize(n,r)
  call zero_initialize(n,x)
  call zero_initialize(n,p)
  call zero_initialize(n,q)
  call zero_initialize(n,s)
  call xcpy(n,pcg%x,x)

  ! ------
  ! pcg start
  ! p=0.0
  ! alpha=0.0
  ! beta=0.0
  ! solve Ms=r
  ! rho0=<s,r>
  ! do itr =1,....
  !   x = x + alpha*p
  !   p = s + beta * p
  !   q = Ap
  !   gamma = <p,q>
  !   alpha=rho0/gamma
  !   r = r - alpha * q
  !   solve Ms=r
  !   rho1 = <s,r>
  !   if(||r||<rmax) exit
  !   beta = rho1/rho0
  !   rho0=rho1
  ! enddo
  ! x = x + alpha*p
  ! ------

  call calc_res(pcg,x,pcg%b,r)
  call calc_dot(n,icomm,pcg%b,pcg%b,bnorm)
  call calc_dot(n,icomm,r,r,rnorm)
  bnorm=sqrt(bnorm)
  rnorm=sqrt(rnorm)

  ! alpha=0.0
  ! beta=0.0
  alpha=0.0d0
  beta=0.0d0
  ! solve Ms=r
  call calc_precon(pcg,r,s)

  ! rho0=<s,r>
  call calc_dot(n,icomm,s,r,rho0)

   ! call MPI_Barrier(icomm,ierr)
  etime=MPI_Wtime()
  rtime(1) =  rtime(1) + etime-stime
  do itr = 1,itrmax
     pcg%reshistory(itr)=rnorm/bnorm 
     if(myrank.eq.0)then
        ! write(*,*)"@ res @",itr,rnorm,rnorm/bnorm
     endif
     !   x = x + alpha*p
     !   p = s + beta * p
      ! call MPI_Barrier(icomm,ierr)
     stime=MPI_Wtime()
     call CG_axpy2(n,x,p,s,alpha,beta)
      ! call MPI_Barrier(icomm,ierr)
     etime=MPI_Wtime()
     rtime(5) =  rtime(5) + etime-stime

     !   q = Ap
      ! call MPI_Barrier(icomm,ierr)
     stime=MPI_Wtime()
     call calc_Matvec(pcg,p,q)
      ! call MPI_Barrier(icomm,ierr)
     etime=MPI_Wtime()
     rtime(2) =  rtime(2) + etime-stime
    
     !   gamma = <p,q>
      ! call MPI_Barrier(icomm,ierr)
     stime=MPI_Wtime()
     call calc_dot(n,icomm,p,q,gamma)
      ! call MPI_Barrier(icomm,ierr)
     etime=MPI_Wtime()
     rtime(6) =  rtime(6) + etime-stime

     !   alpha=rho0/gamma
     alpha=rho0/gamma

     !   r = r - alpha * q
     tmp1=-alpha
     tmp2=1.0d0
      ! call MPI_Barrier(icomm,ierr)
     stime=MPI_Wtime()
     call axpy(n,tmp1,tmp2,q,r)
      ! call MPI_Barrier(icomm,ierr)
     etime=MPI_Wtime()
     rtime(7) =  rtime(7) + etime-stime

     !   solve Ms=r
      ! call MPI_Barrier(icomm,ierr)
     stime=MPI_Wtime()
     call  calc_precon(pcg,r,s)
      ! call MPI_Barrier(icomm,ierr)
     etime=MPI_Wtime()
     rtime(3) =  rtime(3) + etime-stime

     !   rho1 = <s,r>
      ! call MPI_Barrier(icomm,ierr)
     stime=MPI_Wtime()
     call CG_calc_dot2(n,icomm,s,r,rho1,rnorm)

      ! call MPI_Barrier(icomm,ierr)
     etime=MPI_Wtime()
     rtime(4) =  rtime(4) + etime-stime
     rnorm=sqrt(rnorm)
     if( (rnorm/bnorm<rtolmax) ) then
!        pcg%reshistory(itr+1)=rnorm/bnorm 
        if(myrank.eq.0)then
!           write(*,*)"@ res @",itr+1,rnorm,rnorm/bnorm
        endif
        exit
     endif
     !   beta = rho1/rho0
     !   rho0=rho1
     beta = rho1/rho0
     rho0=rho1


#if 0
  call calc_res(pcg,x,pcg%b,q)
  call calc_dot(n,icomm,q,q,rnorm)
  rnorm=sqrt(rnorm)
  if(myrank.eq.0)then
     write(*,*)"@ true res @",rnorm,rnorm/bnorm
  endif
#endif
  enddo

  ! x = x + alpha*p
  tmp1=1.0d0
  call axpy(n,alpha,tmp1,p,x)

  call MPI_Barrier(icomm,ierr)
  etime0=MPI_Wtime()
  ! ------
  ! 真の残差を確認
#if 0
  call calc_res(pcg,x,pcg%b,r)
  call calc_dot(n,icomm,r,r,rnorm)
  rnorm=sqrt(rnorm)
  if(myrank.eq.0)then
     write(*,*)"@ true res @",rnorm,rnorm/bnorm
  endif
#endif

 ! if(myrank.eq.0)then
 !     write(*,*)"cg      ",etime0-stime0
 !     write(*,*)"    init     ",rtime(1)
 !     write(*,*)"    MatVec   ",rtime(2)
 !     write(*,*)"    Precon   ",rtime(3)
 !     write(*,*)"    calc_dot2",rtime(4)
 !     write(*,*)"    axpy2    ",rtime(5)
 !     write(*,*)"    calc_dot ",rtime(6)
 !     write(*,*)"    axpy     ",rtime(7)
 !  endif

#if 0
  recvbuf=0.0d0
  call MPI_Barrier(icomm,ierr)
  call MPI_Gather(&
       rtime(2), 1, MPI_REAL8, &
       recvbuf, 1, MPI_REAL8, &
       0, icomm,ierr)

  call system("sleep 1")
  if(myrank.eq.0)then
     do i=1,pcg%commVec%npes
        write(*,*)i," -- MatVec -- ",recvbuf(i)
     enddo
  endif
! -----

#if 0
  call system("sleep 1")
  rtime(1)=pcg%commVec%nrk_data_max
  call MPI_Gather(&
       rtime(1), 1, MPI_REAL8, &
       recvbuf, 1, MPI_REAL8, &
       0, icomm,ierr)
  call system("sleep 1")
  if(myrank.eq.0)then
     do i=1,pcg%commVec%npes
        write(*,*)i," -- nrk_data_max -- ",recvbuf(i)
     enddo
  endif
  call system("sleep 1")
! -----
  call system("sleep 1")
  rtime(1)=pcg%commVec%nrk_data_max_send
  call MPI_Gather(&
       rtime(1), 1, MPI_REAL8, &
       recvbuf, 1, MPI_REAL8, &
       0, icomm,ierr)
  call system("sleep 1")
  if(myrank.eq.0)then
     do i=1,pcg%commVec%npes
        write(*,*)i," -- nrk_data_max_send -- ",recvbuf(i)
     enddo
  endif
  call system("sleep 1")

! -----
  call system("sleep 1")
  rtime(1)=pcg%commVec%nrk_rank
  call MPI_Gather(&
       rtime(1), 1, MPI_REAL8, &
       recvbuf, 1, MPI_REAL8, &
       0, icomm,ierr)
  call system("sleep 1")
  if(myrank.eq.0)then
     do i=1,pcg%commVec%npes
        write(*,*)i," -- nrk_rank -- ",recvbuf(i)
     enddo 
 endif
  call system("sleep 1")
! -----
  call system("sleep 1")
  rtime(1)=pcg%commVec%nrk_rank_send
  call MPI_Gather(&
       rtime(1), 1, MPI_REAL8, &
       recvbuf, 1, MPI_REAL8, &
       0, icomm,ierr)
  call system("sleep 1")
  if(myrank.eq.0)then
     do i=1,pcg%commVec%npes
        write(*,*)i," -- nrk_rank_send -- ",recvbuf(i)
     enddo
  endif
  call system("sleep 1")
#endif

#endif

  call xcpy(n,x,pcg%x)

  pcg%niter=min(itrmax,itr)
  if(itrmax<=itr)then
     pcg%iret=PARCEL_RETURNCODE_MAXITERATION
  else
     pcg%iret=PARCEL_RETURNCODE_NORMAL
  endif

  ! ------
  deallocate(x)
  deallocate(r)
  deallocate(p)
  deallocate(q)
  deallocate(s)

end subroutine KSP_pCG
