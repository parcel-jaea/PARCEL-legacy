! #define timer 0

module KSP_CG_CPU

contains
  !   x = x + alpha*p
  !   p = s + beta * p
  subroutine CG_axpy2_FP64_CPU(n,x,p,s,alpha,beta)
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

  end subroutine CG_axpy2_FP64_CPU
 
  !   rho1 = <s,r>
  !   rdot = <r,r>
  subroutine CG_calc_dot2_FP64_CPU(n,icomm,s,r,rho1,rdot)
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

  end subroutine CG_calc_dot2_FP64_CPU

end module KSP_CG_CPU

module KSP_CG_CUDA
  use iso_c_binding
  use PARCEL_integer
  use c_wrapper
  interface
     subroutine cuda_CG_axpy2_FP64(&
          numGPUblock,numGPUthread,&
          n, &
          x, offset_x,&
          p, offset_p, &
          s, offset_s, &
          alpha, &
          beta) &
          bind(c,name="cuda_CG_axpy2_FP64")
       import :: devicemem
       import :: intType
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value :: n
       type(devicemem) x
       integer(kind=intType),value :: offset_x
       type(devicemem) p
       integer(kind=intType),value :: offset_p
       type(devicemem) s
       integer(kind=intType),value :: offset_s
       real*8,value :: alpha
       real*8,value :: beta
     end subroutine cuda_CG_axpy2_FP64

     subroutine cuda_CG_calc_dot2_FP64(&
          icomm, &
          numGPUblock,numGPUthread,&
          n, &
          s, offset_s,&
          r, offset_r,&
          rho1, &
          rdot) &
          bind(c,name="cuda_CG_calc_dot2_FP64")
       import :: devicemem
       import :: intType
       implicit none
       integer,value :: icomm
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value :: n
       type(devicemem) s
       integer(kind=intType),value :: offset_s
       type(devicemem) r
       integer(kind=intType),value :: offset_r
       type(devicemem) rho1
       type(devicemem) rdot
     end subroutine cuda_CG_calc_dot2_FP64

  end interface

end module KSP_CG_CUDA

module KSP_CG
  use KSP_STATUS
  use PARCEL_integer
  use PARCEL_VEC_kernel
  use KSP_CG_CPU
  use KSP_CG_CUDA
contains
  !   x = x + alpha*p
  !   p = s + beta * p
  subroutine CG_axpy2(pcg,x,p,s,alpha,beta)
    implicit none
    type(KSP) :: pcg
     type(STATUS_FLAG) :: flag
     integer(kind=intType) n

    type(PARCEL_VEC) :: x
    type(PARCEL_VEC) :: p
    type(PARCEL_VEC) :: s
    type(PARCEL_VEC) :: alpha
    type(PARCEL_VEC) :: beta
    real*8 alpha_,beta_
    integer(kind=intType) parcel_int
    
    n = pcg%prm%n
    flag = pcg%flag
    select case(flag%machineflag)
    case(PARCEL_CPU_FLAG_ON)
       select case(flag%calcflag)
       case(KSP_DOUBLE_PRECISION)
          call CG_axpy2_FP64_CPU( &
               n,&
               x%vec_FP64,p%vec_FP64,s%vec_FP64,&
               alpha%vec_FP64(1),beta%vec_FP64(1))
       case(KSP_DOUBLE_DOUBLE_PRECISION)

       end select
    case(PARCEL_GPU_FLAG_ON)
       parcel_int = 1
       call get_val_FP64(pcg,alpha,alpha_,parcel_int)
       call get_val_FP64(pcg,beta ,beta_ ,parcel_int)
       call cuda_CG_axpy2_FP64( &
            !          numGPUblock,numGPUthread,&
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            n,&
            x%C_vec_FP64, x%offset, &
            p%C_vec_FP64, p%offset,&
            s%C_vec_FP64, s%offset,&
            alpha_,beta_)
    end select

  end subroutine CG_axpy2
 
  !   rho1 = <s,r>
  !   rdot = <r,r>
  subroutine CG_calc_dot2( &
       pcg, &
       icomm,s,r,rho1,rdot)
    use mpi
    use PARCEL_integer
    implicit none

    type(KSP) :: pcg
    type(STATUS_FLAG) :: flag
    integer(kind=intType) n
    integer icomm

    type(PARCEL_VEC) :: s
    type(PARCEL_VEC) :: r
    type(PARCEL_VEC) :: rho1
    type(PARCEL_VEC) :: rdot

    n = pcg%prm%n
    flag = pcg%flag
    select case(flag%machineflag)
    case(PARCEL_CPU_FLAG_ON)
       select case(flag%calcflag)
       case(KSP_DOUBLE_PRECISION)
          call CG_calc_dot2_FP64_CPU( &
               n,icomm, &
               s%vec_FP64,r%vec_FP64, &
               rho1%vec_FP64(1),rdot%vec_FP64(1))
       case(KSP_DOUBLE_DOUBLE_PRECISION)

       end select
    case(PARCEL_GPU_FLAG_ON)
       call cuda_CG_calc_dot2_FP64(&
            icomm, &
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            n, &
            s%C_vec_FP64, s%offset, &
            r%C_vec_FP64, r%offset, &
            rho1%C_vec_FP64, &
            rdot%C_vec_FP64)        
    end select

  end subroutine CG_calc_dot2

end module KSP_CG

subroutine KSP_pCG( pcg &
     )
  
  use KSP_CG
  use krylov
!  use krylov_kernel
!  use vec_kernel
  use PARCEL_KRYLOV_MATRIX
  use PARCEL_KRYLOV_PRECON
  use PARCEL_VEC_kernel
  use PARCEL_integer

!  use matrix_kernel
  implicit none
  
  type(KSP) :: pcg
  integer(kind=intType) n
  integer icomm
  integer itrmax ! 最大反復回数
  real*8 rtolmax ! 相対残差 (収束条件)
  real*8 abstolmax  !  絶対残差 (収束条件)
  real*8 rnorm_FP64
  real*8 bnorm_FP64
  real*8 rho1_FP64
  real*8 rho0_FP64
  real*8 gamma_FP64
  real*8 alpha_FP64
  real*8 beta_FP64
  real*8 tmp1_FP64
  real*8 tmp2_FP64

  integer myrank
  type(STATUS_FLAG) :: flag

  ! ------------------------
  ! real*8,allocatable :: p(:)
  ! real*8,allocatable :: q(:)
  ! real*8,allocatable :: s(:)
  ! real*8,allocatable :: r(:)
  ! real*8,allocatable :: x(:)

  ! real*8 alpha,beta
  ! real*8 gamma
  ! real*8 rho0,rho1
  ! real*8 rnorm
  ! real*8 bnorm
  ! real*8 tmp1,tmp2
  ! ------------------------
  type(PARCEL_VEC) :: p
  type(PARCEL_VEC) :: q
  type(PARCEL_VEC) :: s
  type(PARCEL_VEC) :: r
  type(PARCEL_VEC) :: x
  type(PARCEL_VEC) :: b

  type(PARCEL_VEC) :: alpha,beta
  type(PARCEL_VEC) :: gamma
  type(PARCEL_VEC) :: rho0,rho1
  type(PARCEL_VEC) :: rnorm
  type(PARCEL_VEC) :: bnorm
  type(PARCEL_VEC) :: tmp1,tmp2
  ! ------------------------

  integer(kind=intType) vecLen
  integer(kind=intType) nColL
  integer(kind=intType) nColU

  integer(kind=intType) i
  integer(kind=intType) itr
  real*8 rtime(40)
  real*8 stime0,etime0
  real*8 stime,etime
  integer ierr
  integer(kind=intType) parcel_int
  integer(kind=intType),parameter :: offset_zero=0

  rtime = 0.0d0

  itrmax = pcg%itrmax
  rtolmax = pcg%rtolmax
  abstolmax = pcg%abstolmax
  myrank = pcg%commVec%myrank
  icomm = pcg%commVec%comm
  n = pcg%prm%n

  call MPI_Barrier(icomm,ierr)
  stime0=MPI_Wtime()
  stime=MPI_Wtime()

  ! --------------

  ! 作業用のメモリを確保
  ! cusparse,magma の場合は袖領域の拡張領域も含んで領域を確保する必要がある
  call parcel_sync(pcg)

  nColL = 0
  nColU = 0
  if(pcg%prm%MatFormat.eq.MATRIX_FORMAT_CRS)then
     nColL = pcg%crsMat%nColL
     nColU = pcg%crsMat%nColU
  endif

  vecLen = n+nColL+nColU

  ! write(*,*)"@ 2 @"
  
  call init_vec(pcg,x,vecLen)
  call init_vec(pcg,b,vecLen)  
  call init_vec(pcg,r,vecLen)
  call init_vec(pcg,p,vecLen)
  call init_vec(pcg,q,vecLen)       
  call init_vec(pcg,s,vecLen)

  parcel_int = 1
  call init_vec(pcg,alpha,parcel_int)
  call init_vec(pcg,beta ,parcel_int)
  call init_vec(pcg,gamma,parcel_int)
  call init_vec(pcg,rho0 ,parcel_int)
  call init_vec(pcg,rho1 ,parcel_int)
  call init_vec(pcg,rnorm,parcel_int)
  call init_vec(pcg,bnorm,parcel_int)
  call init_vec(pcg,tmp1 ,parcel_int)
  call init_vec(pcg,tmp2 ,parcel_int)

  ! write(*,*)"@ 3 @"
  call set_val_FP64(pcg,alpha,0.0d0,parcel_int)
  call set_val_FP64(pcg,beta ,0.0d0,parcel_int)
  call set_val_FP64(pcg,gamma,0.0d0,parcel_int)
  call set_val_FP64(pcg,rho0 ,0.0d0,parcel_int)
  call set_val_FP64(pcg,rho1 ,0.0d0,parcel_int)
  call set_val_FP64(pcg,rnorm,0.0d0,parcel_int)
  call set_val_FP64(pcg,bnorm,0.0d0,parcel_int)
  call set_val_FP64(pcg,tmp1 ,0.0d0,parcel_int)
  call set_val_FP64(pcg,tmp2 ,0.0d0,parcel_int)

  ! call zero_initialize(n,r)
  ! call zero_initialize(n,x)
  ! call zero_initialize(n,p)
  ! call zero_initialize(n,q)
  ! call zero_initialize(n,s)

  ! ゼロ初期化
  ! write(*,*)"@ 4 @"

  parcel_int = 1  
  call set_vec_val(pcg,r,tmp1,parcel_int,vecLen)
  call set_vec_val(pcg,x,tmp1,parcel_int,vecLen)
  call set_vec_val(pcg,p,tmp1,parcel_int,vecLen)
  call set_vec_val(pcg,q,tmp1,parcel_int,vecLen)
  call set_vec_val(pcg,s,tmp1,parcel_int,vecLen)

  !  call xcpy(n,pcg%x,x)
  ! 初期解ベクトルをセット
#ifdef CUDA_SOLVER
  call cudaMemcpy_d2d_FP64(pcg%C_x,x%C_vec_FP64,n,offset_zero,offset_zero)
  call cudaMemcpy_d2d_FP64(pcg%C_b,b%C_vec_FP64,n,offset_zero,offset_zero)
#else
  call parcel_sync(pcg)
  call set_vec_FP64(pcg,x,pcg%x,n)
  call parcel_sync(pcg)
  call set_vec_FP64(pcg,b,pcg%b,n)
  call parcel_sync(pcg)
#endif
  ! write(*,*)"input x ",pcg%x(1:3)
  ! write(*,*)"input b ",pcg%b(1:3)

#if 0
!    write(*,*)crsA_cpu(1:3)
    call cudaDeviceSynchronize() 
!    call cuda_view_FP64(n,b%C_vec_FP64)
    do i = 1,n
       write(*,*)i,pcg%x(i),pcg%b(i)
    enddo
    call cudaDeviceSynchronize() 
#endif

  ! *********************************
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
#if 0
  call parcel_sync(pcg)
  call calc_dot(pcg,icomm,b,b,bnorm,n)
  parcel_int = 1  
  call sqrt_val(pcg,bnorm,parcel_int)
  call parcel_sync(pcg)
  call get_val_FP64(pcg,bnorm,bnorm_FP64,parcel_int)
  if(myrank.eq.0)then
     write(*,*)" bnorm = ",bnorm_FP64
  endif

!  call cuda_view_FP64(10,b%C_vec_FP64)

  return
#endif
 ! call system("sleep 1")
  ! write(*,*)"@ 6 @"
 ! call system("sleep 1")
  call parcel_sync(pcg) 
  call calc_res(pcg,x,b,r)

  call parcel_sync(pcg)    
  ! write(*,*)"@ 7 @"
   ! call system("sleep 1")
  call calc_dot(pcg,icomm,b,b,bnorm,n)
  parcel_int = 1
  call sqrt_val(pcg,bnorm,parcel_int)
  call parcel_sync(pcg)

  ! write(*,*)"@ 8 @"
  call calc_dot(pcg,icomm,r,r,rnorm,n)
  call sqrt_val(pcg,rnorm,parcel_int)
  call parcel_sync(pcg)
  
  ! write(*,*)"@ 9 @"
  call get_val_FP64(pcg,bnorm,bnorm_FP64,parcel_int)
  call get_val_FP64(pcg,rnorm,rnorm_FP64,parcel_int)
  call parcel_sync(pcg)

#if 0
  if(myrank.eq.0)then
     write(*,*)" rnorm,bnorm = ",rnorm_FP64,bnorm_FP64
  endif
!  return
#endif

  ! write(*,*)"@ 10 @"
  parcel_int = 1
  call set_val_FP64(pcg,alpha,0.0d0,parcel_int)
  call set_val_FP64(pcg,beta ,0.0d0,parcel_int)
  
  call parcel_sync(pcg) 
  call calc_precon(pcg,r,s)
  call parcel_sync(pcg) 

  call calc_dot(pcg,icomm,s,r,rho0,n)

  call parcel_sync(pcg) 
  etime=MPI_Wtime()
  rtime(1) =  rtime(1) + etime-stime

  call get_val_FP64(pcg,bnorm,bnorm_FP64,parcel_int)
  call parcel_sync(pcg) 

  ! write(*,*)"@ 11 @"
  do itr = 1,itrmax
     call parcel_sync(pcg)
     
#if 0
     call parcel_sync(pcg)
     
     parcel_int = 1
     call get_val_FP64(pcg,rho0 ,rho0_FP64 ,parcel_int)
     call get_val_FP64(pcg,rho1 ,rho1_FP64 ,parcel_int)
     call get_val_FP64(pcg,gamma,gamma_FP64,parcel_int)
     call get_val_FP64(pcg,beta ,beta_FP64 ,parcel_int)
     call get_val_FP64(pcg,alpha,alpha_FP64,parcel_int)
     call parcel_sync(pcg)
!     call get_val_FP64(@,@_FP64,1)
     if(myrank.eq.0)then
        write(*,*)"     rho0_FP64  = ",rho0_FP64 
        write(*,*)"     rho1_FP64  = ",rho1_FP64 
        write(*,*)"     beta_FP64  = ",beta_FP64 
        write(*,*)"     alpha_FP64 = ",alpha_FP64
        write(*,*)"     gamma_FP64 = ",gamma_FP64
        ! write(*,*)itr,"gamma_FP64 = ",gamma_FP64,myrank
     endif
#endif

     parcel_int = 1
     call get_val_FP64(pcg,rnorm,rnorm_FP64,parcel_int)

     call set_val_FP64(pcg,tmp1,rnorm_FP64,parcel_int)
     call div_val_FP64(pcg,tmp1,bnorm_FP64,parcel_int)
     call parcel_sync(pcg)
     call get_val_FP64(pcg,tmp1,pcg%reshistory(itr),parcel_int)
     call parcel_sync(pcg)
     if(myrank.eq.0)then
!        write(*,*)"@ res @",itr,pcg%reshistory(itr)
     endif

     !   x = x + alpha*p
     !   p = s + beta * p
     ! call MPI_Barrier(icomm,ierr)
     stime=MPI_Wtime()
#if 1
     call parcel_sync(pcg)
     call CG_axpy2(pcg,x,p,s,alpha,beta)
     call parcel_sync(pcg)
#else
!     call CG_axpy2(flag,x,p,s,alpha,beta)
     parcel_int = 1
     call set_val_FP64(pcg,tmp2,1.0d0,parcel_int)
     call axpy(pcg,alpha,tmp2,p,x,n)
     call axpy(pcg,tmp2,beta,s,p,n)
#endif
     ! call MPI_Barrier(icomm,ierr)
     call parcel_sync(pcg)
     etime=MPI_Wtime()
     rtime(5) =  rtime(5) + etime-stime

     !   q = Ap
     ! call MPI_Barrier(icomm,ierr)
     stime=MPI_Wtime()
     call calc_Matvec(pcg,p,q)
! debug     call calc_Matvec(pcg,b,q)

      ! call MPI_Barrier(icomm,ierr)
     call parcel_sync(pcg)
     etime=MPI_Wtime()
     rtime(2) =  rtime(2) + etime-stime
    
     !   gamma = <p,q>
     ! call MPI_Barrier(icomm,ierr)
     stime=MPI_Wtime()
     call parcel_sync(pcg)
     call calc_dot(pcg,icomm,p,q,gamma,n)
     call parcel_sync(pcg)
      ! call MPI_Barrier(icomm,ierr)
     etime=MPI_Wtime()
     rtime(6) =  rtime(6) + etime-stime

     !   alpha=rho0/gamma
     parcel_int = 1
     call get_val_FP64(pcg,rho0,rho0_FP64,parcel_int)
     call get_val_FP64(pcg,gamma,gamma_FP64,parcel_int)
     ! if(myrank.eq.0)then
     !    write(*,*)"     rho0_FP64  = ",rho0_FP64 
     !    write(*,*)"     gamma_FP64 = ",gamma_FP64
     ! endif

     ! call set_val_FP64(alpha,rho0%vec_FP64(1),1)
     ! call div_val_FP64(alpha,gamma%vec_FP64(1),1)
     call set_val_FP64(pcg,alpha,rho0_FP64 ,parcel_int)
     call div_val_FP64(pcg,alpha,gamma_FP64,parcel_int)

     !   r = r - alpha * q
     ! tmp1=-alpha
     ! tmp2=1.0d0
     !     call set_val_FP64(tmp1,alpha%vec_FP64(1),1)
     stime=MPI_Wtime()
#if 1
     parcel_int = 1
     call parcel_sync(pcg)
     call get_val_FP64(pcg,alpha,alpha_FP64,parcel_int)
     call parcel_sync(pcg)
     call set_val_FP64(pcg,tmp1,alpha_FP64,parcel_int)
     call parcel_sync(pcg)
     call mul_val_FP64(pcg,tmp1,-1.0d0,parcel_int)
     call parcel_sync(pcg)
     call set_val_FP64(pcg,tmp2,1.0d0,parcel_int)
     call parcel_sync(pcg)
     call axpy(pcg,tmp1,tmp2,q,r,n)
     call parcel_sync(pcg)
#else
     call set_val_FP64(pcg,tmp1,0.0d0,parcel_int)
     call set_val_FP64(pcg,tmp2,1.0d0,parcel_int)
     call axpy(pcg,tmp1,tmp2,q,r,n)
#endif
     etime=MPI_Wtime()
     rtime(7) =  rtime(7) + etime-stime

     !   solve Ms=r
      ! call MPI_Barrier(icomm,ierr)
     stime=MPI_Wtime()
     call  calc_precon(pcg,r,s)
      ! call MPI_Barrier(icomm,ierr)
     call parcel_sync(pcg)
     etime=MPI_Wtime()
     rtime(3) =  rtime(3) + etime-stime

     ! call parcel_sync(pcg)
     ! return 
     !   rho1 = <s,r>
     ! call MPI_Barrier(icomm,ierr)
     stime=MPI_Wtime()

     call parcel_sync(pcg)
#if 1
     call CG_calc_dot2(pcg,icomm,s,r,rho1,rnorm)
!     write(*,*)myrank,rnorm%vec_FP64(1)
#else
     call parcel_sync(pcg)
     call calc_dot(pcg,icomm,s,r,rho1,n)
     call parcel_sync(pcg)
     call calc_dot(pcg,icomm,r,r,rnorm,n)
     call parcel_sync(pcg)
#endif
     call parcel_sync(pcg)

      ! call MPI_Barrier(icomm,ierr)
     etime=MPI_Wtime()
     rtime(4) =  rtime(4) + etime-stime
!     rnorm=sqrt(rnorm)
     parcel_int = 1
     call sqrt_val(pcg,rnorm,parcel_int)
     call parcel_sync(pcg)
     call get_val_FP64(pcg,rnorm,rnorm_FP64,parcel_int)
     call parcel_sync(pcg)

#if 0
     call parcel_sync(pcg)
     ! call get_val_FP64(pcg,rho1 ,tmp1_FP64,1)
     call get_val_FP64(pcg,rnorm,tmp2_FP64,parcel_int)
     call parcel_sync(pcg)
     ! write(*,*)myrank,"rho1  = ",tmp2_FP64,rho1_FP64

     if(myrank.eq.0)then
!        write(*,*)myrank,"rnorm = ",tmp2_FP64/bnorm_FP64,rnorm_FP64/bnorm_FP64
     endif
     call parcel_sync(pcg)
#endif
     
     call parcel_sync(pcg)
     if( (rnorm_FP64/bnorm_FP64 < rtolmax) ) then
        !        write(*,*)myrank ,"exit !!!! ",rnorm_FP64,rnorm_FP64/bnorm_FP64
        call parcel_sync(pcg)
        exit
     endif

     if( (rnorm_FP64/bnorm_FP64 > (1.0d8) ) ) then
        if(myrank.eq.0)then
           write(*,*)" not convergence : residual norm = ",rnorm_FP64/bnorm_FP64
        endif
        call parcel_sync(pcg)
        exit
     endif

     !   beta = rho1/rho0
     !   rho0=rho1
     !     beta = rho1/rho0
     !     rho0=rho1
     parcel_int = 1
     call get_val_FP64(pcg,rho1,rho1_FP64,parcel_int)
     call get_val_FP64(pcg,rho0,rho0_FP64,parcel_int)
     call set_val_FP64(pcg,beta,rho1_FP64,parcel_int)
     call div_val_FP64(pcg,beta,rho0_FP64,parcel_int)
     call set_val_FP64(pcg,rho0,rho1_FP64,parcel_int)

     ! call set_val_FP64(beta,rho1%vec_FP64(1),1)
     ! call div_val_FP64(beta,rho0%vec_FP64(1),1)
     ! call set_val_FP64(rho0,rho1%vec_FP64(1),1)

#if 0
     call parcel_sync(pcg)
     call calc_res(pcg,x,b,q)
     call parcel_sync(pcg)
     call calc_dot(pcg,icomm,q,q,rnorm,n)
     call parcel_sync(pcg)
     parcel_int = 1
     call get_val_FP64(pcg,rnorm,rnorm_FP64,parcel_int)
     call parcel_sync(pcg)
     call sqrt_val(pcg,rnorm,parcel_int)
     if(myrank.eq.0)then
        ! write(*,*)"@ true res @",rnorm_FP64,rnorm_FP64/bnorm_FP64
     endif
#endif

#if 0
     if(itr.eq.itrmax)then
        call parcel_sync(pcg) 
        call free_vec(pcg,x)
        call free_vec(pcg,b)
        call free_vec(pcg,r)
        call free_vec(pcg,p)
        call free_vec(pcg,q)
        call free_vec(pcg,s)

        call free_vec(pcg,alpha)
        call free_vec(pcg,beta)
        call free_vec(pcg,gamma)
        call free_vec(pcg,rho0)
        call free_vec(pcg,rho1)
        call free_vec(pcg,rnorm)
        call free_vec(pcg,bnorm)
        call free_vec(pcg,tmp1)
        call free_vec(pcg,tmp2)
        call parcel_sync(pcg) 
        return
     endif
#endif

  enddo

  ! x = x + alpha*p
!  tmp1=1.0d0
  call parcel_sync(pcg)
  parcel_int = 1
  call set_val_FP64(pcg,tmp1,1.0d0,parcel_int)
  call parcel_sync(pcg)
  call axpy(pcg,alpha,tmp1,p,x,n)

  call MPI_Barrier(icomm,ierr)
  etime0=MPI_Wtime()
  ! ------
  ! 真の残差を確認
#if 1
  parcel_int = 1  
  call calc_res(pcg,x,b,r)
  call calc_dot(pcg,icomm,r,r,rnorm,n)
  call sqrt_val(pcg,rnorm,parcel_int)
  call get_val_FP64(pcg,rnorm,rnorm_FP64,parcel_int)

  if(myrank.eq.0)then
     ! write(*,*)"@ true res @",rnorm_FP64,rnorm_FP64/bnorm_FP64
  endif
#endif

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

#if 0
  recvbuf=0.0d0
  call MPI_Barrier(icomm,ierr)
  call MPI_Gather(&
       rtime(2), 1, MPI_REAL8, &
       recvbuf, 1, MPI_REAL8, &
       0, icomm,ierr)

  ! call system("sleep 1")
  if(myrank.eq.0)then
     do i=1,pcg%commVec%npes
        write(*,*)i," -- MatVec -- ",recvbuf(i)
     enddo
  endif
! -----

#if 0
  ! call system("sleep 1")
  rtime(1)=pcg%commVec%nrk_data_max
  call MPI_Gather(&
       rtime(1), 1, MPI_REAL8, &
       recvbuf, 1, MPI_REAL8, &
       0, icomm,ierr)
  ! call system("sleep 1")
  if(myrank.eq.0)then
     do i=1,pcg%commVec%npes
        write(*,*)i," -- nrk_data_max -- ",recvbuf(i)
     enddo
  endif
  ! call system("sleep 1")
! -----
  ! call system("sleep 1")
  rtime(1)=pcg%commVec%nrk_data_max_send
  call MPI_Gather(&
       rtime(1), 1, MPI_REAL8, &
       recvbuf, 1, MPI_REAL8, &
       0, icomm,ierr)
  ! call system("sleep 1")
  if(myrank.eq.0)then
     do i=1,pcg%commVec%npes
        write(*,*)i," -- nrk_data_max_send -- ",recvbuf(i)
     enddo
  endif
  ! call system("sleep 1")

! -----
  ! call system("sleep 1")
  rtime(1)=pcg%commVec%nrk_rank
  call MPI_Gather(&
       rtime(1), 1, MPI_REAL8, &
       recvbuf, 1, MPI_REAL8, &
       0, icomm,ierr)
  ! call system("sleep 1")
  if(myrank.eq.0)then
     do i=1,pcg%commVec%npes
        write(*,*)i," -- nrk_rank -- ",recvbuf(i)
     enddo 
 endif
  ! call system("sleep 1")
! -----
  ! call system("sleep 1")
  rtime(1)=pcg%commVec%nrk_rank_send
  call MPI_Gather(&
       rtime(1), 1, MPI_REAL8, &
       recvbuf, 1, MPI_REAL8, &
       0, icomm,ierr)
  ! call system("sleep 1")
  if(myrank.eq.0)then
     do i=1,pcg%commVec%npes
        write(*,*)i," -- nrk_rank_send -- ",recvbuf(i)
     enddo
  endif
  ! call system("sleep 1")
#endif

#endif

!  call xcpy(flag,x,pcg%x)
  call parcel_sync(pcg)
#ifdef CUDA_SOLVER
  call cudaMemcpy_d2d_FP64(x%C_vec_FP64,pcg%C_x,n,offset_zero,offset_zero)
#else
  call get_vec_FP64(pcg,x,pcg%x,n)
#endif
  call parcel_sync(pcg)

  pcg%niter=min(itrmax,itr)
  if(itrmax<=itr)then
     pcg%iret=PARCEL_RETURNCODE_MAXITERATION
  else
     pcg%iret=PARCEL_RETURNCODE_NORMAL
  endif

  ! ------
  ! deallocate(x)
  ! deallocate(r)
  ! deallocate(p)
  ! deallocate(q)
  ! deallocate(s))
  ! call MPI_Barrier(icomm,ierr)
  ! write(*,*)" PCG FREE "
  ! ! call system("sleep 3")
  ! call MPI_Barrier(icomm,ierr)
#if 1
  call parcel_sync(pcg) 
  call free_vec(pcg,x)
  call free_vec(pcg,b)
  call free_vec(pcg,r)
  call free_vec(pcg,p)
  call free_vec(pcg,q)
  call free_vec(pcg,s)

  call free_vec(pcg,alpha)
  call free_vec(pcg,beta)
  call free_vec(pcg,gamma)
  call free_vec(pcg,rho0)
  call free_vec(pcg,rho1)
  call free_vec(pcg,rnorm)
  call free_vec(pcg,bnorm)
  call free_vec(pcg,tmp1)
  call free_vec(pcg,tmp2)
  call parcel_sync(pcg) 
#endif
  ! call MPI_Barrier(icomm,ierr)
  ! write(*,*)" PCG END"
  ! ! call system("sleep 3")
  ! call MPI_Barrier(icomm,ierr)

end subroutine KSP_pCG

#if 0
subroutine JACOBI_METHOD( jacobi &
     )
  
  use KSP_CG
   use krylov
!  use krylov_kernel
!  use vec_kernel
  use PARCEL_KRYLOV_MATRIX
  use PARCEL_KRYLOV_PRECON
  use PARCEL_VEC_kernel
  use PARCEL_integer
  implicit none
  
  type(KSP) :: jacobi
  integer(kind=intType) n
  integer icomm
  integer itrmax ! 最大反復回数
  real*8 rtolmax ! 相対残差 (収束条件)
  real*8 abstolmax  !  絶対残差 (収束条件)
  real*8 rnorm_FP64
  real*8 bnorm_FP64

  integer myrank
  type(STATUS_FLAG) :: flag

  ! ------------------------
  ! real*8,allocatable :: p(:)
  ! real*8,allocatable :: q(:)
  ! real*8,allocatable :: s(:)
  ! real*8,allocatable :: r(:)
  ! real*8,allocatable :: x(:)

  ! real*8 alpha,beta
  ! real*8 gamma
  ! real*8 rho0,rho1
  ! real*8 rnorm
  ! real*8 bnorm
  ! real*8 tmp1,tmp2
  ! ------------------------
  type(PARCEL_VEC) :: r
  type(PARCEL_VEC) :: x
  type(PARCEL_VEC) :: b
  ! ------------------------
  type(PARCEL_VEC) :: tmp1
  type(PARCEL_VEC) :: rnorm
  type(PARCEL_VEC) :: bnorm

  integer(kind=intType) itr
  real*8 rtime(40)
  real*8 stime0,etime0
  real*8 stime,etime
  integer ierr

  integer(kind=intType) vecLen
  integer(kind=intType) nColL
  integer(kind=intType) nColU


  rtime=0.0d0

  itrmax=jacobi%itrmax
  rtolmax=jacobi%rtolmax
  abstolmax=jacobi%abstolmax
  myrank=jacobi%commVec%myrank
  icomm=jacobi%commVec%comm
  n=jacobi%prm%n

  call MPI_Barrier(icomm,ierr)
  stime0=MPI_Wtime()
  stime=MPI_Wtime()

  call init_staus_flag(flag)
#if 0
 call set_machineflag(flag,PARCEL_CPU_FLAG_ON)
#else
  call set_machineflag(flag,PARCEL_GPU_FLAG_ON)
#endif
  ! CPU - GPU 両方利用する場合
  !  call set_machineflag(flag,PARCEL_ALLARCH_FLAG_ON)
  call set_dataflag(flag,KSP_DOUBLE_PRECISION)
  call set_calcflag(flag,KSP_DOUBLE_PRECISION)
  ! ここではjacobi内のflagから設定しているが将来はKSP構造体で設定したものを受け取る
  jacobi%flag%machineflag = flag%machineflag
  jacobi%flag%dataflag = flag%dataflag 
  jacobi%flag%calcflag = flag%calcflag 

  ! --------------

  call parcel_sync(jacobi)
  nColL = jacobi%crsMat%nColL
  nColU = jacobi%crsMat%nColU
  vecLen = n+nColL+nColU
  
  parcel_int = 1
  call init_vec(jacobi,tmp1 ,parcel_int)
  call init_vec(jacobi,rnorm ,parcel_int)
  call init_vec(jacobi,bnorm ,parcel_int)
  call set_val_FP64(jacobi,tmp1 ,0.0d0,parcel_int)

  call init_vec(jacobi,x,vecLen)
  call init_vec(jacobi,b,vecLen)  
  call init_vec(jacobi,r,vecLen)

  ! ゼロ初期化
  parcel_int = 1
  call set_vec_val(jacobi,x,tmp1,parcel_int,n)
  call set_vec_val(jacobi,b,tmp1,parcel_int,n)
  call set_vec_val(jacobi,r,tmp1,parcel_int,n)

  ! 初期解ベクトルをセット
  call parcel_sync(jacobi)
  call set_vec_FP64(jacobi,x,jacobi%x,n)
  call parcel_sync(jacobi)
  call set_vec_FP64(jacobi,b,jacobi%b,n)
  call parcel_sync(jacobi)

  call parcel_sync(jacobi)
  call calc_dot(jacobi,icomm,b,b,bnorm,n)
  call sqrt_val(jacobi,bnorm,parcel_int)
  call parcel_sync(jacobi)  
  call get_val_FP64(jacobi,bnorm,bnorm_FP64,parcel_int)
  call parcel_sync(jacobi)

  call parcel_sync(jacobi) 
  call calc_res(jacobi,x,b,r)
  call parcel_sync(jacobi)    
  call calc_dot(jacobi,icomm,r,r,rnorm,n)
  call parcel_sync(jacobi)
  call get_val_FP64(jacobi,rnorm,rnorm_FP64,parcel_int)
  call sqrt_val(jacobi,rnorm,parcel_int)

  if(myrank.eq.0)then
!     write(*,*)" rnorm,bnorm = ",rnorm_FP64,bnorm_FP64
  endif
  do itr = 1,itrmax
     ! -------
     ! JACOBI
     call precon_JACOBI(jacobi,b,x,1,10)
     ! -------
     call parcel_sync(jacobi) 
     call calc_res(jacobi,x,b,r)
     call parcel_sync(jacobi)    
     call calc_dot(jacobi,icomm,r,r,rnorm,n)
     call sqrt_val(jacobi,rnorm,1)
     call parcel_sync(jacobi)

     call get_val_FP64(jacobi,rnorm,rnorm_FP64,1)
     call set_val_FP64(jacobi,tmp1,rnorm_FP64,1)
     call div_val_FP64(jacobi,tmp1,bnorm_FP64,1)
     call parcel_sync(jacobi)
     call get_val_FP64(jacobi,tmp1,jacobi%reshistory(itr),1)
     call parcel_sync(jacobi)
     if(myrank.eq.0)then
!        write(*,*)"@ res @",itr,jacobi%reshistory(itr)
     endif

     if( (rnorm_FP64/bnorm_FP64 < rtolmax) ) then
        call parcel_sync(jacobi)
        exit
     endif
     if( (rnorm_FP64/bnorm_FP64 > (1.0d4) ) ) then
        if(myrank.eq.0)then
           write(*,*)" not convergence : residual norm = ",rnorm_FP64/bnorm_FP64
        endif
        call parcel_sync(jacobi)
        exit
     endif

  enddo
  call get_vec_FP64(jacobi,x,jacobi%x,n)

#if 1
  call parcel_sync(jacobi) 
  call free_vec(jacobi,x)
  call free_vec(jacobi,b)
  call free_vec(jacobi,r)

  call free_vec(jacobi,rnorm)
  call free_vec(jacobi,bnorm)
  call free_vec(jacobi,tmp1)
  call parcel_sync(jacobi) 
#endif
  ! call MPI_Barrier(icomm,ierr)
  ! write(*,*)" JACOBI END"
  ! ! call system("sleep 3")
  ! call MPI_Barrier(icomm,ierr)
  
end subroutine JACOBI_METHOD
#endif

! #ifdef CUDA_SOLVER
!   call cudaMemcpy_d2d_FP64(@%C_x,x%C_vec_FP64,n,offset_zero,offset_zero)
!   call cudaMemcpy_d2d_FP64(@%C_b,b%C_vec_FP64,n,offset_zero,offset_zero)
! #else

! #ifdef CUDA_SOLVER
!   call cudaMemcpy_d2d_FP64(x%C_vec_FP64,@%C_x,n,offset_zero,offset_zero)
! #else

