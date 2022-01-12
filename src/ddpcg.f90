#include "ddfunc.h"

module ddKSP_CG
  use dd_func
  implicit none
contains
  !   x = x + alpha*p
  !   p = s + beta * p
  subroutine ddCG_axpy2(&
       n,&
       x_hi,x_lo,&
       p_hi,p_lo,&
       s_hi,s_lo,&
       alpha_hi,alpha_lo,&
       beta_hi,beta_lo)

    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 x_hi(n)
    real*8 x_lo(n)
    real*8 p_hi(n)
    real*8 p_lo(n)
    real*8 s_hi(n)
    real*8 s_lo(n)
    real*8 alpha_hi,alpha_lo
    real*8 beta_hi,beta_lo
    integer(kind=intType) i
    real*8 c_hi,c_lo

    real*8 ddmul_var1
    real*8 ddmul_var2
    real*8 ddadd_var1
    real*8 ddadd_var2

    !$OMP PARALLEL DO PRIVATE(i,c_hi,c_lo,&
    !$OMP ddmul_var1, &
    !$OMP ddadd_var1, &
    !$OMP ddmul_var2, &
    !$OMP ddadd_var2)
    do i=1,n
       ! x(i)=x(i)+alpha*p(i)
       ! p(i)=s(i)+beta*p(i)
#ifdef Manual_inline
       ddmul_calc1(alpha_hi,alpha_lo,p_hi(i),p_lo(i),c_hi,c_lo)
       ddadd_calc1(x_hi(i),x_lo(i),c_hi,c_lo,x_hi(i),x_lo(i))
       ddmul_calc2(beta_hi,beta_lo,p_hi(i),p_lo(i),c_hi,c_lo)
       ddadd_calc2(s_hi(i),s_lo(i),c_hi,c_lo,p_hi(i),p_lo(i))
#else
       call dd_mul(alpha_hi,alpha_lo,p_hi(i),p_lo(i),c_hi,c_lo)
       call dd_add(x_hi(i),x_lo(i),c_hi,c_lo,x_hi(i),x_lo(i))

       call dd_mul(beta_hi,beta_lo,p_hi(i),p_lo(i),c_hi,c_lo)
       call dd_add(s_hi(i),s_lo(i),c_hi,c_lo,p_hi(i),p_lo(i))
#endif
    enddo
    !$OMP END PARALLEL DO 

  end subroutine ddCG_axpy2
 
  !   rho1 = <s,r>
  !   rdot = <r,r>
  subroutine ddCG_calc_dot2(&
       n,icomm,&
       s_hi,s_lo,&
       r_hi,r_lo,&
       rho1_hi,rho1_lo,&
       rdot_hi,rdot_lo)

    use mpi
    use omp_lib
    use common_kernel,only : get_maxthreads
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer icomm
    real*8 s_hi(n)
    real*8 s_lo(n)
    real*8 r_hi(n)
    real*8 r_lo(n)
    real*8 rho1_hi
    real*8 rho1_lo
    real*8 rdot_hi
    real*8 rdot_lo

    real*8 tmp01_hi
    real*8 tmp01_lo
    real*8 tmp02_hi
    real*8 tmp02_lo

    real*8 tmp1_hi
    real*8 tmp1_lo
    real*8 tmp2_hi
    real*8 tmp2_lo
    real*8 sbuf_hi(2),sbuf_lo(2)
    real*8 rbuf_hi(2),rbuf_lo(2)

    integer(kind=intType) i,ii
    integer ierr
    real*8 c_hi,c_lo
    real*8 ddmul_var1
    real*8 ddmul_var2
    real*8 ddadd_var1
    real*8 ddadd_var2
    integer max_numth
    integer myth

    real*8 tmp1_block_hi(tmpblocksize)
    real*8 tmp1_block_lo(tmpblocksize)
    real*8 tmp2_block_hi(tmpblocksize)
    real*8 tmp2_block_lo(tmpblocksize)
    real*8 ,allocatable:: tmp1_th_hi(:)
    real*8 ,allocatable:: tmp1_th_lo(:)
    real*8 ,allocatable:: tmp2_th_hi(:)
    real*8 ,allocatable:: tmp2_th_lo(:)

    tmp01_hi=0.0d0
    tmp01_lo=0.0d0
    tmp02_hi=0.0d0
    tmp02_lo=0.0d0

    call get_maxthreads(max_numth)
    allocate( tmp1_th_hi(max_numth) )
    allocate( tmp1_th_lo(max_numth) )
    allocate( tmp2_th_hi(max_numth) )
    allocate( tmp2_th_lo(max_numth) )

    !$OMP PARALLEL PRIVATE(i,ii,tmp1_hi,tmp1_lo,tmp2_hi,tmp2_lo,c_hi,c_lo,myth, &
    !$OMP tmp1_block_hi,tmp1_block_lo, &
    !$OMP tmp2_block_hi,tmp2_block_lo, &
    !$OMP ddmul_var1, &
    !$OMP ddmul_var2, &
    !$OMP ddadd_var1, &
    !$OMP ddadd_var2) 

    myth  = omp_get_thread_num()
    tmp1_hi=0.0d0
    tmp1_lo=0.0d0
    tmp2_hi=0.0d0
    tmp2_lo=0.0d0

#if 0
    !$OMP DO 
    do i=1,n
       ! tmp1=tmp1+s(i)*r(i)
       ! tmp2=tmp2+r(i)*r(i)
#ifdef Manual_inline
        ddmul_calc1(r_hi(i),r_lo(i),s_hi(i),s_lo(i),c_hi,c_lo)
        ddadd_calc1(tmp1_hi,tmp1_lo,c_hi,c_lo,tmp1_hi,tmp1_lo)

        ddmul_calc2(r_hi(i),r_lo(i),r_hi(i),r_lo(i),c_hi,c_lo)
        ddadd_calc2(tmp2_hi,tmp2_lo,c_hi,c_lo,tmp2_hi,tmp2_lo)
#else
       call dd_mul(r_hi(i),r_lo(i),s_hi(i),s_lo(i),c_hi,c_lo)
       call dd_add(tmp1_hi,tmp1_lo,c_hi,c_lo,tmp1_hi,tmp1_lo)

       call dd_mul(r_hi(i),r_lo(i),r_hi(i),r_lo(i),c_hi,c_lo)
       call dd_add(tmp2_hi,tmp2_lo,c_hi,c_lo,tmp2_hi,tmp2_lo)
#endif
    enddo
    !$OMP END DO NOWAIT 
#else

    tmp1_block_hi = 0.0d0
    tmp1_block_lo = 0.0d0
    tmp2_block_hi = 0.0d0
    tmp2_block_lo = 0.0d0
    !$OMP DO 
    do ii=1,n,tmpblocksize
       do i=ii,min(n,ii+tmpblocksize-1)
          ! tmp1=tmp1+s(i)*r(i)
          ! tmp2=tmp2+r(i)*r(i)
          ddmul_calc1(r_hi(i),r_lo(i),s_hi(i),s_lo(i),c_hi,c_lo)
          ddadd_calc1(tmp1_block_hi(i-ii+1),tmp1_block_lo(i-ii+1),c_hi,c_lo,tmp1_block_hi(i-ii+1),tmp1_block_lo(i-ii+1))
          
          ddmul_calc2(r_hi(i),r_lo(i),r_hi(i),r_lo(i),c_hi,c_lo)
          ddadd_calc2(tmp2_block_hi(i-ii+1),tmp2_block_lo(i-ii+1),c_hi,c_lo,tmp2_block_hi(i-ii+1),tmp2_block_lo(i-ii+1))
       enddo
    enddo
    !$OMP END DO NOWAIT 
    do i = 1,tmpblocksize 
       ddadd_calc1(tmp1_hi,tmp1_lo,tmp1_block_hi(i),tmp1_block_lo(i),tmp1_hi,tmp1_lo)
       ddadd_calc2(tmp2_hi,tmp2_lo,tmp2_block_hi(i),tmp2_block_lo(i),tmp2_hi,tmp2_lo)
    enddo

#endif

#if 0

#ifdef Manual_inline

    !$OMP CRITICAL
    ddadd_calc1(tmp01_hi,tmp01_lo,tmp1_hi,tmp1_lo,tmp01_hi,tmp01_lo)
    !$OMP END CRITICAL

    !$OMP CRITICAL
    ddadd_calc2(tmp02_hi,tmp02_lo,tmp2_hi,tmp2_lo,tmp02_hi,tmp02_lo)
    !$OMP END CRITICAL
    
#else
    !$OMP CRITICAL
    call dd_add(tmp01_hi,tmp01_lo,tmp1_hi,tmp1_lo,tmp01_hi,tmp01_lo)
    call dd_add(tmp02_hi,tmp02_lo,tmp2_hi,tmp2_lo,tmp02_hi,tmp02_lo)
    !$OMP END CRITICAL
#endif

#else

    tmp1_th_hi(myth+1) = tmp1_hi
    tmp1_th_lo(myth+1) = tmp1_lo
    tmp2_th_hi(myth+1) = tmp2_hi
    tmp2_th_lo(myth+1) = tmp2_lo

#endif

    !$OMP END PARALLEL 

    do i = 1,max_numth
       ddadd_calc1(tmp01_hi,tmp01_lo,tmp1_th_hi(i),tmp1_th_lo(i),tmp01_hi,tmp01_lo)
       ddadd_calc2(tmp02_hi,tmp02_lo,tmp2_th_hi(i),tmp2_th_lo(i),tmp02_hi,tmp02_lo)
    enddo

    sbuf_hi(1)=tmp01_hi
    sbuf_lo(1)=tmp01_lo
    sbuf_hi(2)=tmp02_hi
    sbuf_lo(2)=tmp02_lo

    rbuf_hi(1)=0.0d0
    rbuf_lo(1)=0.0d0
    rbuf_hi(2)=0.0d0
    rbuf_lo(2)=0.0d0
    call MPI_Allreduce_dd(&
         sbuf_hi,sbuf_lo,rbuf_hi,rbuf_lo,&
         2,MPI_REAL8,MPI_SUM,icomm,ierr)    
    rho1_hi=rbuf_hi(1)
    rho1_lo=rbuf_lo(1)
    rdot_hi=rbuf_hi(2)
    rdot_lo=rbuf_lo(2)
    deallocate( tmp1_th_hi )
    deallocate( tmp1_th_lo )
    deallocate( tmp2_th_hi )
    deallocate( tmp2_th_lo )

  end subroutine ddCG_calc_dot2
end module ddKSP_CG

subroutine ddKSP_pCG( pcg&
     )
  use PARCEL_define
  use ddKSP_CG
  use krylov
  use ddkrylov_kernel
  use vec_kernel
  use ddvec_kernel
  use PARCEL_integer
!  use krylov_kernel,only : view_KSP_PARAM
  implicit none
  type(KSP) :: pcg
  integer(kind=intType) n
  integer icomm
  integer itrmax ! 最大反復回数
  real*8 rtolmax ! 相対残差 (収束条件)
  real*8 abstolmax  !  絶対残差 (収束条件)
  integer myrank

  real*8,allocatable :: p_hi(:)
  real*8,allocatable :: p_lo(:)
  real*8,allocatable :: q_hi(:)
  real*8,allocatable :: q_lo(:)
  real*8,allocatable :: s_hi(:)
  real*8,allocatable :: s_lo(:)
  real*8,allocatable :: r_hi(:)
  real*8,allocatable :: r_lo(:)
  real*8,allocatable :: x_hi(:)
  real*8,allocatable :: x_lo(:)
  real*8,allocatable :: b_hi(:)
  real*8,allocatable :: b_lo(:)

  real*8 alpha_hi,alpha_lo
  real*8 beta_hi,beta_lo
  real*8 gamma_hi,gamma_lo
  real*8 rho0_hi,rho0_lo
  real*8 rho1_hi,rho1_lo
  real*8 rnorm_hi,rnorm_lo
  real*8 bnorm_hi,bnorm_lo
  integer itr
  real*8 tmp1_hi,tmp1_lo
  real*8 tmp2_hi,tmp2_lo
  real*8 rtime(40)
  real*8 stime0,etime0
  real*8 stime,etime
  integer ierr

  integer PRECISION_A
  integer PRECISION_b
  integer PRECISION_x
!  call view_KSP_PARAM(pcg)

  rtime=0.0d0

  itrmax=pcg%itrmax
  rtolmax=pcg%rtolmax
  abstolmax=pcg%abstolmax
  myrank=pcg%commVec%myrank
  icomm=pcg%commVec%comm
  n=pcg%prm%n
  PRECISION_A=pcg%PRECISION_A
  PRECISION_b=pcg%PRECISION_b
  PRECISION_x=pcg%PRECISION_x

  call MPI_Barrier(icomm,ierr)
  stime0=MPI_Wtime()
  stime=MPI_Wtime()

  allocate(p_hi(n))
  allocate(p_lo(n))
  allocate(q_hi(n))
  allocate(q_lo(n))
  allocate(s_hi(n))
  allocate(s_lo(n))
  allocate(r_hi(n))
  allocate(r_lo(n))
  allocate(x_hi(n))
  allocate(x_lo(n))
  allocate(b_hi(n))
  allocate(b_lo(n))

  call zero_initialize(n,r_hi)
  call zero_initialize(n,r_lo)
  call zero_initialize(n,x_hi)
  call zero_initialize(n,x_lo)
  call zero_initialize(n,b_hi)
  call zero_initialize(n,b_lo)
  call zero_initialize(n,p_hi)
  call zero_initialize(n,p_lo)
  call zero_initialize(n,q_hi)
  call zero_initialize(n,q_lo)
  call zero_initialize(n,s_hi)
  call zero_initialize(n,s_lo)

  call xcpy(n,pcg%x_hi,x_hi)
  if(PRECISION_x.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
     call xcpy(n,pcg%x_lo,x_lo)
  endif
  call xcpy(n,pcg%b_hi,b_hi)
  if(PRECISION_b.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
     call xcpy(n,pcg%b_lo,b_lo)
  endif

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
  call ddcalc_res(pcg,x_hi,x_lo,b_hi,b_lo,r_hi,r_lo)
  call ddcalc_dot(n,icomm,&
       b_hi,b_lo,&
       b_hi,b_lo,&
       bnorm_hi,bnorm_lo)
  call ddcalc_dot(n,icomm,&
       r_hi,r_lo,&
       r_hi,r_lo,&
       rnorm_hi,rnorm_lo)

!  bnorm=sqrt(bnorm)
  call dd_sqrt(bnorm_hi,bnorm_lo,bnorm_hi,bnorm_lo)
!  rnorm=sqrt(rnorm)
  call dd_sqrt(rnorm_hi,rnorm_lo,rnorm_hi,rnorm_lo)

  ! alpha=0.0
  ! beta=0.0
  alpha_hi=0.0d0
  alpha_lo=0.0d0
  beta_hi=0.0d0
  beta_lo=0.0d0
  ! solve Ms=r
  call ddcalc_precon(pcg,r_hi,r_lo,s_hi,s_lo)

  ! rho0=<s,r>
  call ddcalc_dot(&
       n,icomm,&
       s_hi,s_lo,&
       r_hi,r_lo,&
       rho0_hi,rho0_lo)

  etime=MPI_Wtime()
  rtime(1) =  rtime(1) + etime-stime
  do itr = 1,itrmax
     pcg%reshistory(itr)=rnorm_hi/bnorm_hi
     if(myrank.eq.0)then
!        write(*,*)"@ res @",itr,rnorm_hi,rnorm_hi/bnorm_hi
     endif
     !   x = x + alpha*p
     !   p = s + beta * p
     stime=MPI_Wtime()
     call ddCG_axpy2(&
          n,&
          x_hi,x_lo,&
          p_hi,p_lo,&
          s_hi,s_lo,&
          alpha_hi,alpha_lo,&
          beta_hi,beta_lo)
     etime=MPI_Wtime()
     rtime(5) =  rtime(5) + etime-stime

     !   q = Ap
     stime=MPI_Wtime()
     call ddcalc_Matvec(pcg,p_hi,p_lo,q_hi,q_lo)
     etime=MPI_Wtime()
     rtime(2) =  rtime(2) + etime-stime

     !   gamma = <p,q>
     stime=MPI_Wtime()
     call ddcalc_dot(n,icomm,p_hi,p_lo,q_hi,q_lo,gamma_hi,gamma_lo)
     etime=MPI_Wtime()
     rtime(6) =  rtime(6) + etime-stime

     !   alpha=rho0/gamma
     call dd_div(rho0_hi,rho0_lo,gamma_hi,gamma_lo,alpha_hi,alpha_lo)

     !   r = r - alpha * q
     call dd_mul(&
          -1.0d0,0.0d0,&
          alpha_hi,alpha_lo,&
          tmp1_hi,tmp1_lo) 
     tmp2_hi=1.0d0
     tmp2_lo=0.0d0

     stime=MPI_Wtime()
     call ddaxpy(&
          n,&
          tmp1_hi,tmp1_lo,&
          tmp2_hi,tmp2_lo,&
          q_hi,q_lo,&
          r_hi,r_lo)
     etime=MPI_Wtime()
     rtime(7) =  rtime(7) + etime-stime

     !   solve Ms=r
     stime=MPI_Wtime()
     call ddcalc_precon(pcg,r_hi,r_lo,s_hi,s_lo)
     etime=MPI_Wtime()
     rtime(3) =  rtime(3) + etime-stime

     !   rho1 = <s,r>
     stime=MPI_Wtime()
     call ddCG_calc_dot2(&
          n,icomm,&
          s_hi,s_lo,&
          r_hi,r_lo,&
          rho1_hi,rho1_lo,&
          rnorm_hi,rnorm_lo)
     etime=MPI_Wtime()
     rtime(4) =  rtime(4) + etime-stime

!     rnorm=sqrt(rnorm)
     call dd_sqrt(rnorm_hi,rnorm_lo,rnorm_hi,rnorm_lo)

     if( (rnorm_hi/bnorm_hi<rtolmax) ) then
!        pcg%reshistory(itr+1)=rnorm/bnorm 
        if(myrank.eq.0)then
!           write(*,*)"@ res @",itr+1,rnorm,rnorm/bnorm
        endif
        exit
     endif
     !   beta = rho1/rho0
     !   rho0=rho1
     call dd_div(rho1_hi,rho1_lo,rho0_hi,rho0_lo,beta_hi,beta_lo)
     rho0_hi=rho1_hi
     rho0_lo=rho1_lo
  enddo

  ! x = x + alpha*p
  tmp1_hi=1.0d0
  tmp1_lo=0.0d0
  call ddaxpy(&
       n,&
       alpha_hi,alpha_lo, &
       tmp1_hi,tmp1_lo, &
       p_hi,p_lo, &
       x_hi,x_lo)

  ! ------
  ! 真の残差を確認
#if 0
  !call calc_res(pcg,x,pcg%b,r)
  !call calc_dot(n,icomm,r,r,rnorm)
  call ddcalc_res(pcg,x_hi,x_lo,b_hi,b_lo,r_hi,r_lo)
  call ddcalc_dot(n,icomm,r_hi,r_lo,r_hi,r_lo,rnorm_hi,rnorm_lo)
  call dd_sqrt(rnorm_hi,rnorm_lo,rnorm_hi,rnorm_lo)

  if(myrank.eq.0)then
     write(*,*)"@ true res @",rnorm_hi,rnorm_hi/bnorm_hi
  endif
#endif

  call MPI_Barrier(icomm,ierr)
  etime0=MPI_Wtime()
#ifdef timer
 if(myrank.eq.0)then
     write(*,*)"pcg      ",etime0-stime0
     write(*,*)"    init     ",rtime(1)
     write(*,*)"    MatVec   ",rtime(2)
     write(*,*)"    Precon   ",rtime(3)
     write(*,*)"    calc_dot2",rtime(4)
     write(*,*)"    axpy2    ",rtime(5)
     write(*,*)"    calc_dot ",rtime(6)
     write(*,*)"    axpy     ",rtime(7)
  endif
#endif

  pcg%niter=min(itrmax,itr)
  if(itrmax<=itr)then
     pcg%iret=PARCEL_RETURNCODE_MAXITERATION
  else
     pcg%iret=PARCEL_RETURNCODE_NORMAL
  endif

  call xcpy(n,x_hi,pcg%x_hi)
  if(PRECISION_x.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
     call xcpy(n,x_lo,pcg%x_lo)
  endif

  ! ------
  deallocate(p_hi)
  deallocate(p_lo)
  deallocate(q_hi)
  deallocate(q_lo)
  deallocate(s_hi)
  deallocate(s_lo)
  deallocate(r_hi)
  deallocate(r_lo)
  deallocate(x_hi)
  deallocate(x_lo)
  deallocate(b_hi)
  deallocate(b_lo)

end subroutine ddKSP_pCG
