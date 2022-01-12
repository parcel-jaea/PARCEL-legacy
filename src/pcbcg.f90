module KSP_CBCG
contains
  ! SS=eta*AS - zeta * r
  subroutine CBCG_axpy2(n,SS,AS,r,eta,zeta)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 SS(n)
    real*8 AS(n)
    real*8 r(n)
    real*8 eta,zeta
    integer(kind=intType) i
    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
       SS(i)=eta*AS(i)-zeta*r(i)
    enddo
    !$OMP END PARALLEL DO 

  end subroutine CBCG_axpy2

  ! SS(i) = 2.0*eta*AS-2.0*zeta*SS(i-1)-SS(i-2)
  subroutine CBCG_axpy3(n,SS0,AS,SS1,SS2,eta2,zeta2)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 SS0(n)
    real*8 AS(n)
    real*8 SS1(n)
    real*8 SS2(n)
    real*8 eta2,zeta2
    integer(kind=intType) i

    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
       SS0(i)=eta2*AS(i)-zeta2*SS1(i)-SS2(i)
    enddo
    !$OMP END PARALLEL DO 

  end subroutine CBCG_axpy3

  subroutine make_chebyshev_basis(&
       pcbcg,&
       n,k,&
       Emax,Emin,&
       r,S,AS,SS)
    use krylov
    use krylov_kernel
    use PARCEL_integer
    implicit none
    type(KSP) :: pcbcg
    integer(kind=intType) n,k
    real*8 Emax,Emin
    real*8 r(n)
    real*8 S(n*k)
    real*8 SS(n*k)
    real*8 AS(n*k)

    real*8 eta
    real*8 zeta

    real*8 eta2
    real*8 zeta2

    integer(kind=intType) i
    integer(kind=intType) itmp

    eta   = 2.0d0/(Emax-Emin);
    zeta  = (Emax+Emin)/(Emax-Emin)
    eta2  = eta *2.0d0
    zeta2 = zeta*2.0d0

    call calc_precon(pcbcg,r,S)
    call calc_Matvec(pcbcg,S,AS)
    if(k.eq.1)then
       return 
    endif
    call CBCG_axpy2(n,SS,AS,r,eta,zeta)
    call calc_precon(pcbcg,SS,S(1+n))
    call calc_Matvec(pcbcg,S(1+n),AS(1+n))

    do i=2,k-1
       ! write(*,*)"make basis --",i
       itmp=2
       if(i.eq.itmp)then
          call CBCG_axpy3(n,&
               SS(1+(i-1)*n),&
               AS(1+(i-1)*n),&
               SS(1+(i-2)*n),&
               r,&
               eta2,zeta2)
       else
          call CBCG_axpy3(n,&
               SS(1+(i-1)*n),&
               AS(1+(i-1)*n),&
               SS(1+(i-2)*n),&
               SS(1+(i-3)*n),&
               eta2,zeta2)
       endif
       call calc_precon(pcbcg,SS(1+(i-1)*n), S(1+(i+0)*n))
       call calc_Matvec(pcbcg, S(1+(i+0)*n),AS(1+(i+0)*n))
    enddo

  end subroutine make_chebyshev_basis

  subroutine getEigenValue(method,Eigenflag,Emax,Emin)
    use krylov
    use PARCEL_define
    type(KSP) :: method
    integer Eigenflag
    real*8 Emax,Emin
    method%EigenMax_Real=0.0d0
    method%EigenMin_Real=0.0d0
    select case (Eigenflag)
    case(KSP_EIGENVALUE_TYPE_POWER)
       call power_method_symMatrix(method)
    case(KSP_EIGENVALUE_TYPE_CAARNOLDI)
       call KSP_pCAARNOLDI(method)
    end select
    Emax=method%EigenMax_Real
    Emin=method%EigenMin_Real

  end subroutine getEigenValue

  ! -------
  ! x  = x + Q * alpha
  ! Q  = S - Q * Beta
  ! AQ = AS - AQ * Beta
  ! Qt*AQ 
  ! Qtr
  subroutine CBCG_gemm_dot_local(n,k,x,Q,S,AQ,AS,r,alpha,Beta,QtAQ,Qtr,rnorm,nvecBlock)
    use vec_kernel
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) k
    real*8 x(n)
    real*8 Q(n*k)
    real*8 S(n*k)
    real*8 AQ(n*k)
    real*8 AS(n*k)
    real*8 r(n)
    real*8 alpha(k)
    real*8 Beta(k*k)
    real*8 QtAQ(k*k)
    real*8 Qtr(k)
    real*8 rnorm
    integer nvecBlock
    integer(kind=intType) i,j,l,jj

    QtAQ=0.0d0
    Qtr=0.0d0
    rnorm=0.0d0

    !$OMP PARALLEL DO PRIVATE(jj,l,j,i)  REDUCTION(+:QtAQ,Qtr,rnorm)
    do jj=1,n,nvecBlock 

       do l=1,k
          do j=jj,min(n,jj+nvecBlock-1)
             x(j) = x(j) + Q(j+(l-1)*n)*alpha(l)
          enddo
          do i=1,k
             do j=jj,min(n,jj+nvecBlock-1)
                S(j+(i-1)*n)=S(j+(i-1)*n)-Q(j+(l-1)*n)*Beta(l+(i-1)*k)
             enddo
          enddo
       enddo
       do l=1,k
          do j=jj,min(n,jj+nvecBlock-1)
             Q(j+(l-1)*n)=S(j+(l-1)*n)
          enddo
       enddo

       do l=1,k
          do i=1,k
             do j=jj,min(n,jj+nvecBlock-1)
                AS(j+(i-1)*n)=AS(j+(i-1)*n)-AQ(j+(l-1)*n)*Beta(l+(i-1)*k)
             enddo
          enddo
       enddo
       do l=1,k
          do j=jj,min(n,jj+nvecBlock-1)
             AQ(j+(l-1)*n)=AS(j+(l-1)*n)
          enddo
       enddo

       do l=1,k
          do i=1,l
             do j=jj,min(n,jj+nvecBlock-1)
                QtAQ(i+(l-1)*k)=QtAQ(i+(l-1)*k)+Q(j+(l-1)*n)*AQ(j+(i-1)*n)
             enddo
          enddo
          do j=jj,min(n,jj+nvecBlock-1)
             Qtr(l)=Qtr(l)+Q(j+(l-1)*n)*r(j)
          enddo
       enddo
       do j=jj,min(n,jj+nvecBlock-1)
          rnorm=rnorm+r(j)*r(j)
       enddo
    enddo
    !$OMP END PARALLEL DO 

  end subroutine CBCG_gemm_dot_local

  subroutine CBCG_ALLreduce(k,icomm,Qtr,QtAQ,rnorm)
    use mpi
    use PARCEL_integer
    implicit none
    integer(kind=intType) k
    integer icomm
    real*8 Qtr(k)
    real*8 QtAQ(k*k)
    real*8 rnorm

    real*8 sbuf((k+1)*k/2+k+1)
    real*8 rbuf((k+1)*k/2+k+1)
    integer(kind=intType) i,j,jj
    integer ierr
    integer ntmp

    ! All Reduce -- QtAQ,QtR
    jj=1
    do i=1,k
       do j=1,i
          sbuf(jj)=QtAQ(j+k*(i-1))
          jj=jj+1
       enddo
    enddo
    do i=1,k
       sbuf(jj)=Qtr(i)
       jj=jj+1
    enddo
    sbuf(jj)=rnorm
    jj=jj+1

    ntmp=(k+1)*k/2+k+1
    call MPI_Allreduce(sbuf,rbuf,ntmp,MPI_REAL8,MPI_SUM,icomm,ierr)    

    jj=1
    do i=1,k
       do j=1,i
          QtAQ(j+k*(i-1))=rbuf(jj)
          jj=jj+1
       enddo
    enddo
    do i=1,k
       Qtr(i)=rbuf(jj)
       jj=jj+1
    enddo
    rnorm=rbuf(jj)
    jj=jj+1
    do i=1,k
       do j=1,i
          QtAQ(i+k*(j-1))=QtAQ(j+k*(i-1))
       enddo
    enddo

  end subroutine CBCG_ALLreduce

end module KSP_CBCG

subroutine KSP_pCBCG(pcbcg)
  use krylov
  use krylov_kernel
  use PARCEL_define
  use vec_kernel
  use KSP_CBCG
  use PARCEL_integer
  use wrap_DGEMM,only : wrap_BASIC_DGEMM
  implicit none

  type(KSP) :: pcbcg
  integer(kind=intType) n
  integer icomm
  integer itrmax ! 最大反復回数
  real*8 rtolmax ! 相対残差 (収束条件)
  real*8 abstolmax  !  絶対残差 (収束条件)
  integer myrank
  real*8 Emax,Emin
  integer(kind=intType) k ! ステップ数
  integer(kind=intType) msize
  integer nvecBlock

  real*8,allocatable :: x(:)
  real*8,allocatable :: r(:)
  real*8,allocatable :: AQ(:)
  real*8,allocatable :: Q(:)
  real*8,allocatable :: S(:)
  real*8,allocatable :: AS(:)
  real*8,allocatable :: SS(:)
  real*8,allocatable :: alpha(:)
  real*8,allocatable :: alpha_(:)
  real*8,allocatable :: Beta(:)

  real*8,allocatable :: QtAS(:)
  real*8,allocatable :: QtAQ(:)
  real*8,allocatable :: QtAQ_inv(:)
  real*8,allocatable :: Qtr(:)

  real*8 rnorm
  real*8 bnorm
  integer(kind=intType) i,j,jj
  integer(kind=intType) itr 
  integer ierr
  integer(kind=intType)  ntmp
  integer k_int4
  real*8 rtime(40)
  real*8 stime0,etime0
  real*8 stime,etime

  rtime=0.0d0

  itrmax=pcbcg%itrmax
  rtolmax=pcbcg%rtolmax
  abstolmax=pcbcg%abstolmax
  myrank=pcbcg%commVec%myrank
  icomm=pcbcg%commVec%comm
  n=pcbcg%prm%n
  k=pcbcg%CBCG_kstep 
  nvecBlock=pcbcg%nvecBlock

  call MPI_Barrier(icomm,ierr)
  stime0=MPI_Wtime()
  stime=MPI_Wtime()

  allocate(x(n))
  allocate(r(n))
  allocate(AQ(n*k))
  allocate(Q(n*k))
  allocate(S(n*k))
  allocate(SS(n*k))
  allocate(AS(n*k))
  allocate(alpha(k))
  allocate(alpha_(k))
  allocate(Beta(k*k))

  allocate(QtAS(k*k))
  allocate(QtAQ(k*k))
  allocate(QtAQ_inv(k*k))

  allocate(Qtr(k))

  call zero_initialize(n,x)
  call zero_initialize(n,r)

  msize=n*k
  call zero_initialize(msize,AQ)
  call zero_initialize(msize,Q)
  call zero_initialize(msize,S)
  call zero_initialize(msize,AS)
  call zero_initialize(msize,SS)

  alpha=0.0d0
  alpha_=0.0d0
  Beta=0.0d0

  QtAS=0.0d0
  QtAQ=0.0d0
  QtAQ_inv=0.0d0
  Qtr=0.0d0

  call xcpy(n,pcbcg%x,x)
  etime=MPI_Wtime()
  rtime(1) =  rtime(1) + etime-stime
  
  stime=MPI_Wtime()
  call getEigenValue(pcbcg,pcbcg%CBCG_Eigenflag,Emax,Emin)
  etime=MPI_Wtime()
  rtime(2) =  rtime(2) + etime-stime

  if(Emax.eq.Emin)then
     if(Emax<0.0d0)then
        Emax=0.0d0
     else
        Emin=0.0d0
     endif
  endif
  if( (Emax.eq.0.0d0).and.(Emin.eq.0.0d0) )then
     Emax = 1.0d0
  endif

  if(myrank.eq.0)then
     write(*,*)" Emax,Emin = ",Emax,Emin
  endif
  ! ------
  ! pcbcg start
  ! r = b-Ax
  ! S = T0(AM^-1)r,T1(AM^-1)r,T2(AM^-1)r,........Tk-1(AM^-1)r
  ! Q = S
  ! Qt*AQ,Qt*r
  ! do itr =1,....
  !   alpha = (Qt*AQ)^-1 * (Qt*r)
  !   r = r - AQ * alpha
  !   S = T0(AM^-1)r,T1(AM^-1)r,T2(AM^-1)r,........Tk-1(AM^-1)r
  !   QtA*S
  !   Beta=(Qt*AQ)^-1 * (Qt*AS)
  !   x  = x + Q * alpha
  !   Q  = S - Q * Beta
  !   AQ = AS - AQ * Beta
  !   Qt*AQ,Qt*r
  ! enddo
  ! ------
  
  stime=MPI_Wtime()
  ! r = b-Ax
  call calc_res(pcbcg,x,pcbcg%b,r)
  call calc_dot(n,icomm,pcbcg%b,pcbcg%b,bnorm)
  call calc_dot(n,icomm,r,r,rnorm)
  bnorm=sqrt(bnorm)
  rnorm=sqrt(rnorm)
  !   if(myrank.eq.0)then
  !      write(*,*)"@ res @",1,rnorm,rnorm/bnorm
  !   endif

  ! S = T0(AM^-1)r,T1(AM^-1)r,T2(AM^-1)r,........Tk-1(AM^-1)r
  call make_chebyshev_basis(&
       pcbcg,&
       n,k,&
       Emax,Emin,&
       r,S,AS,SS)
  !   do i=1,k
  !      call write_vec(5, S(1+(i-1)*n))
  !      call write_vec(5,AS(1+(i-1)*n))
  !      write(*,*)"----------"
  !   enddo

  ! Q = S
  msize=n*k
  call xcpy(msize, S, Q)
  call xcpy(msize,AS,AQ)

  ! Qt*AQ 
  call calc_sym_gemm_dot_local(n,k,AQ,Q,QtAQ,nvecBlock) 

  ! Qtr
  j = 1
  call calc_gemm_dot_local(n,k,j,Q,r,Qtr,nvecBlock)
  call calc_dot_local(n,r,r,rnorm)
  
  etime=MPI_Wtime()
  rtime(1) =  rtime(1) + etime-stime
  do itr = 1,itrmax/k
     ! All Reduce -- QtAQ,QtR
     call CBCG_ALLreduce(k,icomm,Qtr,QtAQ,rnorm)

     !   alpha = (Qt*AQ)^-1 * (Qt*r)
     k_int4=k
     call moore_penrose(QtAQ,QtAQ_inv,k,k)

     ntmp=1
     call wrap_BASIC_DGEMM(k,ntmp,k,QtAQ_inv,Qtr,alpha)

     !   r = r - AQ * alpha
     do i=1,k
        alpha_(i)=-alpha(i)
     enddo

     stime=MPI_Wtime()
     call gemv(n,k,AQ,r,alpha_,nvecBlock)
     etime=MPI_Wtime()
     rtime(3) =  rtime(3) + etime-stime

     ! --------------
     rnorm=sqrt(rnorm)
     !      if((myrank.eq.0))then
     !         write(*,*)"@ res @",1+k*(itr-1),rnorm,rnorm/bnorm
     !      endif
     pcbcg%reshistory(itr)=rnorm/bnorm 
     if( (rnorm/bnorm<rtolmax)  ) then
        exit
     endif
     ! --------------

     !   S = T0(AM^-1)r,T1(AM^-1)r,T2(AM^-1)r,........Tk-1(AM^-1)r
     stime=MPI_Wtime()
     call make_chebyshev_basis(&
          pcbcg,&
          n,k,&
          Emax,Emin,&
          r,S,AS,SS)
     etime=MPI_Wtime()
     rtime(4) =  rtime(4) + etime-stime

     !   QtA*S
     stime=MPI_Wtime()
     call calc_gemm_dot(n,k,k,icomm,Q,AS,QtAS,nvecBlock)
     !   Beta=(Qt*AQ)^-1 * (Qt*AS)
     call wrap_BASIC_DGEMM(k,k,k,QtAQ_inv,QtAS,Beta)
     etime=MPI_Wtime()
     rtime(5) =  rtime(5) + etime-stime

#if 0
     !   x  = x + Q * alpha
     call gemv(n,k,Q,x,alpha)
     !   Q  = S - Q * Beta
     call gemm(n,k, S, Q,Beta)
     !   AQ = AS - AQ * Beta
     call gemm(n,k,AS,AQ,Beta)
     ! Qt*AQ 
     call calc_sym_gemm_dot_local(n,k,AQ,Q,QtAQ) 
     ! Qtr
     call calc_gemm_dot_local(n,k,1,Q,r,Qtr)
     ! rnorm
     call calc_dot_local(n,r,r,rnorm)
#else
     ! x  = x + Q * alpha
     ! Q  = S - Q * Beta
     ! AQ = AS - AQ * Beta
     ! Qt*AQ 
     ! Qtr
     ! <r|r>
     stime=MPI_Wtime()
     call  CBCG_gemm_dot_local(n,k,x,Q,S,AQ,AS,r,alpha,Beta,QtAQ,Qtr,rnorm,nvecBlock)
     etime=MPI_Wtime()
     rtime(6) =  rtime(6) + etime-stime

#endif
     
  enddo

  ! ------
  ! 真の残差を確認
#if 0
  call calc_res(pcbcg,x,pcbcg%b,r)
  call calc_dot(n,icomm,r,r,rnorm)
  rnorm=sqrt(rnorm)
  if(myrank.eq.0)then
     write(*,*)"@ true res @",rnorm,rnorm/bnorm
  endif
#endif
  call xcpy(n,x,pcbcg%x)

  pcbcg%niter=min(itrmax,itr-1)
  if(itrmax<=itr)then
     pcbcg%iret=PARCEL_RETURNCODE_MAXITERATION
  else
     pcbcg%iret=PARCEL_RETURNCODE_NORMAL
  endif

  call MPI_Barrier(icomm,ierr)
  etime0=MPI_Wtime()
#ifdef timer
  if(myrank.eq.0)then
     write(*,*)"pcbcg               ",etime0-stime0
     write(*,*)"    init            ",rtime(1)
     write(*,*)"    EigenValue      ",rtime(2)
     write(*,*)"    gemv            ",rtime(3)
     write(*,*)"    chebyshev_basis ",rtime(4)
     write(*,*)"    gemm_dot        ",rtime(5)
     write(*,*)"    gemm_dot_local  ",rtime(6)

  endif
#endif
  

  deallocate(x)
  deallocate(r)
  deallocate(AQ)
  deallocate(Q)
  deallocate(S)
  deallocate(SS)
  deallocate(AS)
  deallocate(alpha)
  deallocate(alpha_)
  deallocate(Beta)

  deallocate(QtAS)
  deallocate(QtAQ)
  deallocate(QtAQ_inv)

  deallocate(Qtr)

  !  stop
end subroutine KSP_pCBCG
