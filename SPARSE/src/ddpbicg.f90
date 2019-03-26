! -----
! c = a+b 
! ddadd(a_hi, a_lo, b_hi, b_lo, c_hi, c_lo)
! -----
! c =a*b
! ddmul(a_hi, a_lo, b_hi, b_lo, c_hi, c_lo)

module ddKSP_BiCG
  use dd_func
  implicit none
contains
  !   e=r-alpha*q
  subroutine ddBiCG_axpy3( &
       n,&
       r_hi    ,r_lo, &
       q_hi    ,q_lo, &
       e_hi    ,e_lo, &
       alpha_hi,alpha_lo)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 r_hi(n)
    real*8 r_lo(n)

    real*8 q_hi(n)
    real*8 q_lo(n)
    real*8 e_hi(n)
    real*8 e_lo(n)

    real*8 alpha_hi
    real*8 alpha_lo

    integer(kind=intType) i 
    real*8 c_hi,c_lo

    !$OMP PARALLEL DO PRIVATE(i,c_hi,c_lo)
    do i=1,n
!       e(i)=r(i)-alpha*q(i)
       call dd_mul(alpha_hi,alpha_lo,q_hi(i),q_lo(i),c_hi,c_lo)
       call dd_sub(r_hi(i),r_lo(i),c_hi,c_lo,e_hi(i),e_lo(i))
    enddo
    !$OMP END PARALLEL DO 
  end subroutine ddBiCG_axpy3

  !   p=r+beta*(p-c3*q)
  subroutine ddBiCG_axpy4_( &
       n, &
       p_hi,p_lo, &
       r_hi,r_lo, &
       q_hi,q_lo, &
       beta_hi,beta_lo, &
       c3_hi,c3_lo)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 p_hi(n)
    real*8 p_lo(n)
    real*8 r_hi(n)
    real*8 r_lo(n)
    real*8 q_hi(n) 
    real*8 q_lo(n) 
    real*8 beta_hi
    real*8 beta_lo
    real*8 c3_hi
    real*8 c3_lo

    real*8 c_hi,c_lo
    integer(kind=intType) i
    !$OMP PARALLEL DO PRIVATE(i,c_hi,c_lo)
    do i=1,n
!       p(i)=r(i)+beta*( p(i) -c3*q(i) )
       call dd_mul(c3_hi,c3_lo,q_hi(i),q_lo(i),c_hi,c_lo)
       call dd_sub(p_hi(i),p_lo(i),c_hi,c_lo,c_hi,c_lo)
       call dd_mul(beta_hi,beta_lo,c_hi,c_lo,c_hi,c_lo)
       call dd_add(r_hi(i),r_lo(i),c_hi,c_lo,p_hi(i),p_lo(i))
    enddo
    !$OMP END PARALLEL DO 

  end subroutine ddBiCG_axpy4_

  !   r=e-c3*v  
  !   c1  =<r0,r>   
  !   rdot=<r,r>
  subroutine ddBiCG_axpy3_dot( &
       n,icomm,&
       r_hi,r_lo, &
       e_hi,e_lo, &
       v_hi,v_lo, &
       r0_hi,r0_lo, &
       c3_hi,c3_lo, &
       c1_hi,c1_lo, &
       rdot_hi,rdot_lo)
    use mpi
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer icomm
    real*8 r_hi(n)
    real*8 r_lo(n)
    real*8 e_hi(n)
    real*8 e_lo(n)
    real*8 v_hi(n) 
    real*8 v_lo(n) 
    real*8 r0_hi(n) 
    real*8 r0_lo(n) 
    real*8 c3_hi
    real*8 c3_lo
    real*8 c1_hi
    real*8 c1_lo
    real*8 rdot_hi
    real*8 rdot_lo

    real*8 tmp1_hi,tmp1_lo
    real*8 tmp2_hi,tmp2_lo

    real*8 tmp01_hi,tmp01_lo
    real*8 tmp02_hi,tmp02_lo

    real*8 sbuf_hi(2),sbuf_lo(2)
    real*8 rbuf_hi(2),rbuf_lo(2)

    real*8 c_hi,c_lo

    integer(kind=intType) i
    integer ierr

    tmp01_hi=0.0d0
    tmp01_lo=0.0d0
    tmp02_hi=0.0d0
    tmp02_lo=0.0d0
!!@@    !$OMP PARALLEL DO PRIVATE(i) REDUCTION(+:tmp1,tmp2)
    !$OMP PARALLEL PRIVATE(i,tmp1_hi,tmp1_lo,tmp2_hi,tmp2_lo,c_hi,c_lo)
    tmp1_hi=0.0d0
    tmp1_lo=0.0d0
    tmp2_hi=0.0d0
    tmp2_lo=0.0d0

    !$OMP DO
    do i=1,n

!       r(i)=e(i)-c3*v(i)
!       tmp1=tmp1+r(i)*r0(i)
!       tmp2=tmp2+r(i)*r(i)

       call dd_mul(c3_hi,c3_lo,v_hi(i),v_lo(i),c_hi,c_lo)
       call dd_sub(e_hi(i),e_lo(i),c_hi,c_lo,r_hi(i),r_lo(i))
       call dd_mul(r_hi(i),r_lo(i),r0_hi(i),r0_lo(i),c_hi,c_lo)
       call dd_add(tmp1_hi,tmp1_lo,c_hi,c_lo,tmp1_hi,tmp1_lo)
       call dd_mul(r_hi(i),r_lo(i),r_hi(i),r_lo(i),c_hi,c_lo)
       call dd_add(tmp2_hi,tmp2_lo,c_hi,c_lo,tmp2_hi,tmp2_lo)

    enddo
    !$OMP END DO NOWAIT 

    !$OMP CRITICAL
    call dd_add(tmp01_hi,tmp01_lo,tmp1_hi,tmp1_lo,tmp01_hi,tmp01_lo)
    call dd_add(tmp02_hi,tmp02_lo,tmp2_hi,tmp2_lo,tmp02_hi,tmp02_lo)
    !$OMP END CRITICAL

    !$OMP END PARALLEL 

    sbuf_hi(1)=tmp01_hi
    sbuf_lo(1)=tmp01_lo
    sbuf_hi(2)=tmp02_hi
    sbuf_lo(2)=tmp02_lo

    rbuf_hi(1)=0.0d0
    rbuf_lo(1)=0.0d0
    rbuf_hi(2)=0.0d0
    rbuf_lo(2)=0.0d0

!    call MPI_Allreduce(sbuf,rbuf,2,MPI_REAL8,MPI_SUM,icomm,ierr)    
    call MPI_AllReduce_dd( &
         sbuf_hi,sbuf_lo,rbuf_hi,rbuf_lo,&
         2,MPI_REAL8,MPI_SUM,icomm,ierr)

    c1_hi=rbuf_hi(1)
    c1_lo=rbuf_lo(1)

    rdot_hi=rbuf_hi(2)
    rdot_lo=rbuf_lo(2)

  end subroutine ddBiCG_axpy3_dot

  !   x=x+alpha*pp+c3*pe   
  subroutine ddBiCG_axpy4( &
       n,&
       x_hi,x_lo, &
       pp_hi,pp_lo, &
       pe_hi,pe_lo, &
       alpha_hi,alpha_lo, &
       c3_hi,c3_lo)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 x_hi(n)
    real*8 x_lo(n)
    real*8 pp_hi(n)
    real*8 pp_lo(n)
    real*8 pe_hi(n) 
    real*8 pe_lo(n) 
    real*8 alpha_hi
    real*8 alpha_lo
    real*8 c3_hi
    real*8 c3_lo

    real*8 b_hi,b_lo
    real*8 c_hi,c_lo
    integer(kind=intType) i

    !$OMP PARALLEL DO PRIVATE(i,c_hi,c_lo,b_hi,b_lo)
    do i=1,n
!       x(i)=x(i)+alpha*pp(i)+c3*pe(i)
       call dd_mul(c3_hi,c3_lo,pe_hi(i),pe_lo(i),c_hi,c_lo)
       call dd_mul(alpha_hi,alpha_lo,pp_hi(i),pp_lo(i),b_hi,b_lo)
       call dd_add(b_hi,b_lo,c_hi,c_lo,c_hi,c_lo)
       call dd_add(x_hi(i),x_lo(i),c_hi,c_lo,x_hi(i),x_lo(i))
    enddo
    !$OMP END PARALLEL DO 

  end subroutine ddBiCG_axpy4

  ! ev_dot = <e,v>, vv_dot=<v,v> 
  subroutine ddBiCG_calc_dot2(&
       n,icomm,&
       e_hi,e_lo, &
       v_hi,v_lo, &
       ev_dot_hi,ev_dot_lo, &
       vv_dot_hi,vv_dot_lo)
    use mpi
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer icomm
    real*8 e_hi(n)
    real*8 e_lo(n)

    real*8 v_hi(n)
    real*8 v_lo(n)

    real*8 ev_dot_hi
    real*8 ev_dot_lo
    real*8 vv_dot_hi
    real*8 vv_dot_lo

    integer(kind=intType) i
    real*8 tmp1_hi,tmp1_lo
    real*8 tmp2_hi,tmp2_lo
    real*8 tmp01_hi,tmp01_lo
    real*8 tmp02_hi,tmp02_lo
    real*8 sbuf_hi(2),sbuf_lo(2)
    real*8 rbuf_hi(2),rbuf_lo(2)
    integer ierr
    real*8 c_hi,c_lo
    tmp01_hi=0.0d0
    tmp01_lo=0.0d0
    tmp02_hi=0.0d0
    tmp02_lo=0.0d0

    !$OMP PARALLEL PRIVATE(i,tmp1_hi,tmp1_lo,tmp2_hi,tmp2_lo,c_hi,c_lo)
    tmp1_hi=0.0d0
    tmp1_lo=0.0d0
    tmp2_hi=0.0d0
    tmp2_lo=0.0d0

    !$OMP DO
    do i=1,n
       !       tmp1=tmp1+e(i)*v(i)
       !       tmp2=tmp2+v(i)*v(i)
       call dd_mul(e_hi(i),e_lo(i),v_hi(i),v_lo(i),c_hi,c_lo)
       call dd_add(tmp1_hi,tmp1_lo,c_hi,c_lo,tmp1_hi,tmp1_lo)
       call dd_mul(v_hi(i),v_lo(i),v_hi(i),v_lo(i),c_hi,c_lo)
       call dd_add(tmp2_hi,tmp2_lo,c_hi,c_lo,tmp2_hi,tmp2_lo)
    enddo
    !$OMP END DO NOWAIT

    !$OMP CRITICAL
    call dd_add(tmp01_hi,tmp01_lo,tmp1_hi,tmp1_lo,tmp01_hi,tmp01_lo)
    call dd_add(tmp02_hi,tmp02_lo,tmp2_hi,tmp2_lo,tmp02_hi,tmp02_lo)
    !$OMP END CRITICAL

    !$OMP END PARALLEL 

    sbuf_hi(1)=tmp01_hi
    sbuf_lo(1)=tmp01_lo

    sbuf_hi(2)=tmp02_hi
    sbuf_lo(2)=tmp02_lo

    rbuf_hi(1)=0.0d0
    rbuf_lo(1)=0.0d0
    rbuf_hi(2)=0.0d0
    rbuf_lo(2)=0.0d0

    call MPI_AllReduce_dd( &
         sbuf_hi,sbuf_lo,rbuf_hi,rbuf_lo,&
         2,MPI_REAL8,MPI_SUM,icomm,ierr)
    ev_dot_hi=rbuf_hi(1)
    ev_dot_lo=rbuf_lo(1)
    vv_dot_hi=rbuf_hi(2)
    vv_dot_lo=rbuf_lo(2)

  end subroutine ddBiCG_calc_dot2
  
end module ddKSP_BiCG

subroutine ddKSP_pBiCG( pbicg &
     )
  use PARCEL_define
  use krylov
  use ddkrylov_kernel
  use ddvec_kernel
  use vec_kernel
!  use matrix_kernel
  use ddKSP_BiCG
!  use KSP_BiCG

  use PARCEL_integer
  implicit none
  type(KSP) :: pbicg
  integer(kind=intType) n
  integer icomm
  integer itrmax ! 最大反復回数
  real*8 rtolmax ! 相対残差 (収束条件)
  real*8 abstolmax  !  絶対残差 (収束条件)
  integer myrank

  real*8,allocatable :: r0_hi(:)
  real*8,allocatable :: r0_lo(:)
  real*8,allocatable :: p_hi(:)
  real*8,allocatable :: p_lo(:)
  real*8,allocatable :: q_hi(:)
  real*8,allocatable :: q_lo(:)
  real*8,allocatable :: pe_hi(:)
  real*8,allocatable :: pe_lo(:)
  real*8,allocatable :: v_hi(:)
  real*8,allocatable :: v_lo(:)
  real*8,allocatable :: pp_hi(:)
  real*8,allocatable :: pp_lo(:)
  real*8,allocatable :: e_hi(:)
  real*8,allocatable :: e_lo(:)
  real*8,allocatable :: x_hi(:)
  real*8,allocatable :: x_lo(:)
  real*8,allocatable :: b_hi(:)
  real*8,allocatable :: b_lo(:)
  real*8,allocatable :: r_hi(:)
  real*8,allocatable :: r_lo(:)

  real*8 alpha_hi,alpha_lo
  real*8 beta_hi,beta_lo
  real*8 c1_hi,c1_lo
  real*8 c2_hi,c2_lo
  real*8 c3_hi,c3_lo
  real*8 rnorm_hi,rnorm_lo
  real*8 bnorm_hi,bnorm_lo

  integer itr
  real*8 tmp_hi ,tmp_lo
  real*8 tmp1_hi,tmp1_lo
  real*8 tmp2_hi,tmp2_lo

  real*8 stime0
  real*8 stime,etime
  real*8 bicg_time
  real*8 precon_time
  real*8 matvec_time

  integer PRECISION_A
  integer PRECISION_b
  integer PRECISION_x

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

  PRECISION_A=pbicg%PRECISION_A
  PRECISION_b=pbicg%PRECISION_b
  PRECISION_x=pbicg%PRECISION_x

  allocate(x_hi(n))
  allocate(b_hi(n))
  allocate(r_hi(n))
  allocate(r0_hi(n))
  allocate(p_hi(n))
  allocate(q_hi(n))
  allocate(pe_hi(n))
  allocate(v_hi(n))
  allocate(pp_hi(n))
  allocate(e_hi(n))

  allocate(x_lo(n))
  allocate(b_lo(n))
  allocate(r_lo(n))
  allocate(r0_lo(n))
  allocate(p_lo(n))
  allocate(q_lo(n))
  allocate(pe_lo(n))
  allocate(v_lo(n))
  allocate(pp_lo(n))
  allocate(e_lo(n))

  call zero_initialize(n,r_hi)
  call zero_initialize(n,r_lo)
  call zero_initialize(n,x_hi)
  call zero_initialize(n,x_lo)
  call zero_initialize(n,b_hi)
  call zero_initialize(n,b_lo)
  call zero_initialize(n,r0_hi)
  call zero_initialize(n,r0_lo)
  call zero_initialize(n,p_hi)
  call zero_initialize(n,p_lo)
  call zero_initialize(n,q_hi)
  call zero_initialize(n,q_lo)
  call zero_initialize(n,pe_hi)
  call zero_initialize(n,pe_lo)
  call zero_initialize(n,v_hi)
  call zero_initialize(n,v_lo)
  call zero_initialize(n,pp_hi)
  call zero_initialize(n,pp_lo)
  call zero_initialize(n,e_hi)
  call zero_initialize(n,e_lo)

  call xcpy(n,pbicg%x_hi,x_hi)
  if(PRECISION_x.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
     call xcpy(n,pbicg%x_lo,x_lo)
  endif
  call xcpy(n,pbicg%b_hi,b_hi)
  if(PRECISION_b.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
     call xcpy(n,pbicg%b_lo,b_lo)
  endif

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
  c1_hi=0.0d0
  c1_lo=0.0d0
  c2_hi=0.0d0
  c2_lo=0.0d0
  c3_hi=0.0d0
  c3_lo=0.0d0

  ! r=b-Ax
  call ddcalc_res(pbicg,x_hi,x_lo,b_hi,b_lo,r_hi,r_lo)
  ! r0=r
  call xcpy(n,r_hi,r0_hi)
  call xcpy(n,r_lo,r0_lo)

  ! p=r
  call xcpy(n,r_hi,p_hi)
  call xcpy(n,r_lo,p_lo)
  call ddcalc_dot(n,icomm,&
       b_hi,b_lo,&
       b_hi,b_lo,&
       bnorm_hi,bnorm_lo)

  ! c1=<r0,r0>
  call ddcalc_dot(n,icomm,&
       r0_hi,r0_lo, &
       r0_hi,r0_lo, &
       c1_hi,c1_lo)

  rnorm_hi=c1_hi
  rnorm_lo=c1_lo
  !  bnorm=sqrt(bnorm)
  call dd_sqrt(bnorm_hi,bnorm_lo,bnorm_hi,bnorm_lo)
  !  rnorm=sqrt(rnorm)
  call dd_sqrt(rnorm_hi,rnorm_lo,rnorm_hi,rnorm_lo)

  do itr = 1,itrmax
     pbicg%reshistory(itr)=rnorm_hi/bnorm_hi
     call dd_div(&
          rnorm_hi,rnorm_lo,&
          bnorm_hi,bnorm_lo,&
          pbicg%reshistory(itr),tmp_lo)

     if(myrank.eq.0)then
!        write(*,*)"@ res @",itr,rnorm_hi,rnorm_hi/bnorm_hi
     endif

     !   solve Mpp=p
     stime=MPI_Wtime()
     call ddcalc_precon(pbicg,p_hi,p_lo,pp_hi,pp_lo) 
     etime=MPI_Wtime()
     precon_time=precon_time+etime-stime

     !   q=App
     stime=MPI_Wtime()
     call ddcalc_Matvec(pbicg,pp_hi,pp_lo,q_hi,q_lo)
     etime=MPI_Wtime()
     matvec_time=matvec_time+etime-stime

     !   c2=<r0,q>
     call ddcalc_dot(n,icomm,r0_hi,r0_lo,q_hi,q_lo,c2_hi,c2_lo)
     !   alpha=c1/c2 
     call dd_div(c1_hi,c1_lo,c2_hi,c2_lo,alpha_hi,alpha_lo)

     ! write(*,*)"c1",c1_hi
     ! write(*,*)"c2",c2_hi
     ! write(*,*)"alpha",alpha_hi

     !   e=r-alpha*q
     call ddBiCG_axpy3(n, &
          r_hi,r_lo, &
          q_hi,q_lo, &
          e_hi,e_lo, &
          alpha_hi,alpha_lo)

     !   solve Mpe=e
     stime=MPI_Wtime()
     call ddcalc_precon(pbicg,e_hi,e_lo,pe_hi,pe_lo)
     etime=MPI_Wtime()
     precon_time=precon_time+etime-stime

     !   v=Ape
     stime=MPI_Wtime()
     call ddcalc_Matvec(pbicg,pe_hi,pe_lo,v_hi,v_lo)
     etime=MPI_Wtime()
     matvec_time=matvec_time+etime-stime

     !   c3=<e,v>/<v,v>
     call ddBiCG_calc_dot2(n,icomm,&
          e_hi,e_lo,&
          v_hi,v_lo,&
          tmp1_hi,tmp1_lo,&
          tmp2_hi,tmp2_lo)
!     write(*,*)"tmp1",tmp1_hi
!     write(*,*)"tmp2",tmp2_hi

    !  c3=tmp1/tmp2
     call dd_div(tmp1_hi,tmp1_lo,tmp2_hi,tmp2_lo,c3_hi,c3_lo)
     !   x=x+alpha*pp+c3*pe
     call ddBiCG_axpy4(n,&
          x_hi,x_lo,&
          pp_hi,pp_lo,&
          pe_hi,pe_lo, &
          alpha_hi,alpha_lo,&
          c3_hi,c3_lo)

     !   r=e-c3*v
     !   c1=<r0,r>
     call ddBiCG_axpy3_dot(n,icomm,&
          r_hi,r_lo, &
          e_hi,e_lo, &
          v_hi,v_lo, &
          r0_hi,r0_lo, &
          c3_hi,c3_lo, &
          c1_hi,c1_lo, &
          rnorm_hi,rnorm_lo)

!     rnorm=sqrt(rnorm)
     call dd_sqrt(rnorm_hi,rnorm_lo,rnorm_hi,rnorm_lo)

     if( (rnorm_hi/bnorm_hi<rtolmax) ) then
!        pbicg%reshistory(itr+1)=rnorm_hi
        if(myrank.eq.0)then
!           write(*,*)"@ res @",itr+1,rnorm,rnorm/bnorm
        endif
        exit
     endif

     !   beta=c1/(c3*c2)
     call dd_mul(c3_hi,c3_lo,c2_hi,c2_lo,tmp_hi,tmp_lo)
     call dd_div(c1_hi,c1_lo,tmp_hi,tmp_lo,beta_hi,beta_lo)


     !   p=r+beta*(p-c3*q)
     call ddBiCG_axpy4_(n,&
          p_hi,p_lo, &
          r_hi,r_lo, &
          q_hi,q_lo, &
          beta_hi,beta_lo, &
          c3_hi,c3_lo)

  enddo
#if 0
  ! ------
  ! 真の残差を確認
  call calc_res(pbicg,x,b,r)
  call calc_dot(n,icomm,r,r,c1)
  rnorm=sqrt(c1)
  if(myrank.eq.0)then
     write(*,*)"@ true res @",rnorm,rnorm/bnorm
  endif
#endif
  call xcpy(n,x_hi,pbicg%x_hi)
  if(PRECISION_x.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
     call xcpy(n,x_lo,pbicg%x_lo)
  endif

  pbicg%niter=min(itrmax,itr)
  if(itrmax<=itr)then
     pbicg%iret=PARCEL_RETURNCODE_MAXITERATION
  else
     pbicg%iret=PARCEL_RETURNCODE_NORMAL
  endif
  ! ------
  deallocate(r0_hi)
  deallocate(r0_lo)
  deallocate(p_hi)
  deallocate(p_lo)
  deallocate(q_hi)
  deallocate(q_lo)
  deallocate(pe_hi)
  deallocate(pe_lo)
  deallocate(v_hi)
  deallocate(v_lo)
  deallocate(pp_hi)
  deallocate(pp_lo)
  deallocate(e_hi)
  deallocate(e_lo)
  deallocate(x_hi)
  deallocate(x_lo)
  deallocate(b_hi)
  deallocate(b_lo)
  deallocate(r_hi)
  deallocate(r_lo)

  etime=MPI_Wtime()
  bicg_time=etime-stime0
#if 0
  if(myrank.eq.0)then
     write(*,*)"@ bicg     [s] @",bicg_time
     write(*,*)"@   matvec [s] @",matvec_time
     write(*,*)"@   precon [s] @",precon_time
  endif
#endif

!  stop
end subroutine ddKSP_pBiCG
