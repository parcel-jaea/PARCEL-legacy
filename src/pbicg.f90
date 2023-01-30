! #define timer 1

module KSP_BiCG_CPU

contains

  !   r=e-c3*v  
  !   c1  =<r0,r>   
  !   rdot=<r,r>
  subroutine BiCG_axpy3_dot_FP64_CPU(n,icomm,r,e,v,r0,c3,c1,rdot)
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

  end subroutine BiCG_axpy3_dot_FP64_CPU

  ! ev_dot = <e,v>, vv_dot=<v,v> 
  subroutine BiCG_calc_dot2_FP64_CPU(n,icomm,e,v,ev_dot,vv_dot)
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

  end subroutine BiCG_calc_dot2_FP64_CPU

end module KSP_BiCG_CPU

module KSP_BiCG_CUDA
  
  use iso_c_binding
  use PARCEL_integer
  use c_wrapper

  interface
     !   r=e-c3*v  
     !   c1  =<r0,r>   
     !   rdot=<r,r>
     subroutine cuda_BiCG_axpy3_dot_FP64( &
          numGPUblock,numGPUthread,&
          n,icomm,&
          r, offset_r,&
          e, offset_e,&
          v, offset_v,&
          r0,offset_r0,&
          c3,c1,rdot) &
          bind(c,name="cuda_BiCG_axpy3_dot_FP64")
       import :: devicemem
       import :: intType
       implicit none

       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value :: n
       integer,value :: icomm
       type(devicemem) r
       integer(kind=intType),value :: offset_r
       type(devicemem) e
       integer(kind=intType),value :: offset_e
       type(devicemem) v
       integer(kind=intType),value :: offset_v
       type(devicemem) r0
       integer(kind=intType),value :: offset_r0
       real*8,value :: c3
       type(devicemem) c1
       type(devicemem) rdot

     end subroutine cuda_BiCG_axpy3_dot_FP64

     ! ev_dot = <e,v>, vv_dot=<v,v> 
     subroutine cuda_BiCG_calc_dot2_FP64(&
          numGPUblock,numGPUthread,&
          n, &
          icomm, &
          e, offset_e,&
          v, offset_v,&
          ev_dot, &
          vv_dot) &
          bind(c,name="cuda_BiCG_calc_dot2_FP64")
       import :: devicemem
       import :: intType
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value :: n
       integer,value :: icomm
       type(devicemem) e
       integer(kind=intType),value :: offset_e
       type(devicemem) v
       integer(kind=intType),value :: offset_v
       type(devicemem) ev_dot
       type(devicemem) vv_dot

     end subroutine cuda_BiCG_calc_dot2_FP64
     
  end interface

contains

end module KSP_BiCG_CUDA

module KSP_BiCG
  use KSP_STATUS
  use PARCEL_integer
  use PARCEL_VEC_kernel
  use KSP_BiCG_CPU
  use KSP_BiCG_CUDA
contains

  !   r=e-c3*v  
  !   c1  =<r0,r>   
  !   rdot=<r,r>
  subroutine BiCG_axpy3_dot(pbicg,icomm,r,e,v,r0,c3,c1,rdot)
    implicit none
    type(KSP) :: pbicg
    integer icomm

    type(PARCEL_VEC) :: r
    type(PARCEL_VEC) :: e
    type(PARCEL_VEC) :: v
    type(PARCEL_VEC) :: r0
    type(PARCEL_VEC) :: c3
    type(PARCEL_VEC) :: c1
    type(PARCEL_VEC) :: rdot

    type(STATUS_FLAG) :: flag
    integer(kind=intType) n    
    integer(kind=intType) parcel_int
     
    real*8 c3_FP64
    n = pbicg%prm%n
    flag = pbicg%flag

    select case(flag%machineflag)
    case(PARCEL_CPU_FLAG_ON)
       select case(flag%calcflag)
       case(KSP_DOUBLE_PRECISION)
          call BiCG_axpy3_dot_FP64_CPU(&
               n,icomm,&
               r%vec_FP64, &
               e%vec_FP64, &
               v%vec_FP64, &
               r0%vec_FP64, &
               c3%vec_FP64(1), &
               c1%vec_FP64(1), &
               rdot%vec_FP64(1) )
       case(KSP_DOUBLE_DOUBLE_PRECISION)

       end select
    case(PARCEL_GPU_FLAG_ON)
       select case(flag%calcflag)
       case(KSP_DOUBLE_PRECISION)

          parcel_int = 1
          call get_val_FP64(pbicg,c3,c3_FP64,parcel_int)
          call cuda_BiCG_axpy3_dot_FP64( &
               CUDA_NUM_BLOCKS,CUDA_NUM_THREADS,  &
               n,icomm, &
               r%C_vec_FP64 ,  r%offset, &
               e%C_vec_FP64 ,  e%offset, &
               v%C_vec_FP64 ,  v%offset, &
               r0%C_vec_FP64, r0%offset, &
               c3_FP64, &
               c1%C_vec_FP64, &
               rdot%C_vec_FP64)
       case(KSP_DOUBLE_DOUBLE_PRECISION)

       end select
    end select

  end subroutine BiCG_axpy3_dot

  ! ev_dot = <e,v>, vv_dot=<v,v> 
  subroutine BiCG_calc_dot2(pbicg,icomm,e,v,ev_dot,vv_dot)
    implicit none
    type(KSP) :: pbicg
    integer icomm
    type(PARCEL_VEC) :: e
    type(PARCEL_VEC) :: v
    type(PARCEL_VEC) :: ev_dot
    type(PARCEL_VEC) :: vv_dot

    type(STATUS_FLAG) :: flag
    integer(kind=intType) n

    n = pbicg%prm%n
    flag = pbicg%flag

    select case(flag%machineflag)
    case(PARCEL_CPU_FLAG_ON)
       select case(flag%calcflag)
       case(KSP_DOUBLE_PRECISION)
          call BiCG_calc_dot2_FP64_CPU( &
               n, &
               icomm, &
               e%vec_FP64, &
               v%vec_FP64, &
               ev_dot%vec_FP64(1), &
               vv_dot%vec_FP64(1))
       case(KSP_DOUBLE_DOUBLE_PRECISION)
       end select
    case(PARCEL_GPU_FLAG_ON)
       select case(flag%calcflag)
       case(KSP_DOUBLE_PRECISION)
          call cuda_BiCG_calc_dot2_FP64(&
               CUDA_NUM_BLOCKS,CUDA_NUM_THREADS,  &
               n, &
               icomm, &
               e%C_vec_FP64, e%offset,&
               v%C_vec_FP64, v%offset,&
               ev_dot%C_vec_FP64, &
               vv_dot%C_vec_FP64) 
       case(KSP_DOUBLE_DOUBLE_PRECISION)
       end select
    end select


  end subroutine BiCG_calc_dot2

end module KSP_BiCG


subroutine KSP_pBiCG( pbicg &
     )
  use KSP_BiCG
  use krylov
  use PARCEL_KRYLOV_MATRIX
  use PARCEL_KRYLOV_PRECON
  use PARCEL_VEC_kernel
  use PARCEL_integer

  implicit none
  type(KSP) :: pbicg
  integer(kind=intType) n
  integer icomm
  integer itrmax ! 最大反復回数
  real*8 rtolmax ! 相対残差 (収束条件)
  real*8 abstolmax  !  絶対残差 (収束条件)
  integer myrank
  type(STATUS_FLAG) :: flag

#if 0
  real*8,allocatable :: r0(:)
  real*8,allocatable :: p(:)
  real*8,allocatable :: q(:)
  real*8,allocatable :: pe(:)
  real*8,allocatable :: v(:)
  real*8,allocatable :: pp(:)
  real*8,allocatable :: e(:)

  real*8,allocatable :: x(:)
  real*8,allocatable :: r(:)
#else

  type(PARCEL_VEC) :: r0
  type(PARCEL_VEC) :: p
  type(PARCEL_VEC) :: q
  type(PARCEL_VEC) :: pe
  type(PARCEL_VEC) :: v
  type(PARCEL_VEC) :: pp
  type(PARCEL_VEC) :: e
  type(PARCEL_VEC) :: x
  type(PARCEL_VEC) :: b
  type(PARCEL_VEC) :: r

  type(PARCEL_VEC) :: alpha
  type(PARCEL_VEC) :: beta
  type(PARCEL_VEC) :: c1
  type(PARCEL_VEC) :: c2
  type(PARCEL_VEC) :: c3
  type(PARCEL_VEC) :: rnorm
  type(PARCEL_VEC) :: bnorm
  type(PARCEL_VEC) :: tmp1
  type(PARCEL_VEC) :: tmp2
  type(PARCEL_VEC) :: tmp3

#endif

  real*8 alpha_FP64
  real*8 beta_FP64
  real*8 c1_FP64
  real*8 c2_FP64
  real*8 c3_FP64
  real*8 rnorm_FP64
  real*8 bnorm_FP64
  real*8 tmp1_FP64
  real*8 tmp2_FP64
  real*8 tmp3_FP64

  integer(kind=intType) vecLen
  integer(kind=intType) nColL
  integer(kind=intType) nColU

  integer itr
  real*8 rtime(40)
  real*8 stime0,etime0
  real*8 stime,etime
  integer ierr
  integer(kind=intType),parameter :: offset_zero=0
  integer(kind=intType) parcel_int
  
  rtime = 0.0d0

  itrmax = pbicg%itrmax
  rtolmax = pbicg%rtolmax
  abstolmax = pbicg%abstolmax
  myrank = pbicg%commVec%myrank
  icomm = pbicg%commVec%comm
  n = pbicg%prm%n

  call MPI_Barrier(icomm,ierr)
  stime0=MPI_Wtime()
  stime=MPI_Wtime()

#if 0
  call init_staus_flag(flag)
#if 0
  call set_machineflag(flag,PARCEL_CPU_FLAG_ON)
#else
  call set_machineflag(flag,PARCEL_GPU_FLAG_ON)
#endif
  ! ここではpbicg内のflagから設定しているが将来はKSP構造体で設定したものを受け取る
  call set_dataflag(flag,KSP_DOUBLE_PRECISION)
  call set_calcflag(flag,KSP_DOUBLE_PRECISION)
  pbicg%flag%machineflag = flag%machineflag
  pbicg%flag%dataflag = flag%dataflag 
  pbicg%flag%calcflag = flag%calcflag 
#endif

  call parcel_sync(pbicg)

  nColL = 0
  nColU = 0
  if(pbicg%prm%MatFormat.eq.MATRIX_FORMAT_CRS)then
     nColL = pbicg%crsMat%nColL
     nColU = pbicg%crsMat%nColU
  endif
  vecLen = n+nColL+nColU

  call init_vec(pbicg,r0,vecLen)
  call init_vec(pbicg,p ,vecLen)
  call init_vec(pbicg,q ,vecLen)
  call init_vec(pbicg,pe,vecLen)
  call init_vec(pbicg,v ,vecLen)
  call init_vec(pbicg,pp,vecLen)
  call init_vec(pbicg,e ,vecLen)
  call init_vec(pbicg,x ,vecLen)
  call init_vec(pbicg,b ,vecLen)
  call init_vec(pbicg,r ,vecLen)

  parcel_int = 1
  call init_vec(pbicg,alpha,parcel_int)
  call init_vec(pbicg,beta,parcel_int)
  call init_vec(pbicg,c1,parcel_int)
  call init_vec(pbicg,c2,parcel_int)
  call init_vec(pbicg,c3,parcel_int)
  call init_vec(pbicg,rnorm,parcel_int)
  call init_vec(pbicg,bnorm,parcel_int)
  call init_vec(pbicg,tmp1,parcel_int)
  call init_vec(pbicg,tmp2,parcel_int)
  call init_vec(pbicg,tmp3,parcel_int)

  call set_val_FP64(pbicg,alpha,0.0d0,parcel_int)
  call set_val_FP64(pbicg,beta ,0.0d0,parcel_int)
  call set_val_FP64(pbicg,c1   ,0.0d0,parcel_int)
  call set_val_FP64(pbicg,c2   ,0.0d0,parcel_int)
  call set_val_FP64(pbicg,c3   ,0.0d0,parcel_int)
  call set_val_FP64(pbicg,rnorm,0.0d0,parcel_int)
  call set_val_FP64(pbicg,bnorm,0.0d0,parcel_int)
  call set_val_FP64(pbicg,tmp1 ,0.0d0,parcel_int)
  call set_val_FP64(pbicg,tmp2 ,0.0d0,parcel_int)
  call set_val_FP64(pbicg,tmp3 ,0.0d0,parcel_int)

  parcel_int = 1
  call set_vec_val(pbicg,r0,tmp1,parcel_int,vecLen)
  call set_vec_val(pbicg,p ,tmp1,parcel_int,vecLen)
  call set_vec_val(pbicg,q ,tmp1,parcel_int,vecLen)
  call set_vec_val(pbicg,pe,tmp1,parcel_int,vecLen)
  call set_vec_val(pbicg,v ,tmp1,parcel_int,vecLen)
  call set_vec_val(pbicg,pp,tmp1,parcel_int,vecLen)
  call set_vec_val(pbicg,e ,tmp1,parcel_int,vecLen)
  call set_vec_val(pbicg,x ,tmp1,parcel_int,vecLen)
  call set_vec_val(pbicg,b ,tmp1,parcel_int,vecLen)
  call set_vec_val(pbicg,r ,tmp1,parcel_int,vecLen)

#ifdef CUDA_SOLVER
  call cudaMemcpy_d2d_FP64(pbicg%C_x,x%C_vec_FP64,n,offset_zero,offset_zero)
  call cudaMemcpy_d2d_FP64(pbicg%C_b,b%C_vec_FP64,n,offset_zero,offset_zero)
#else
  call parcel_sync(pbicg)
  call set_vec_FP64(pbicg,x,pbicg%x,n)
  call parcel_sync(pbicg)
  call set_vec_FP64(pbicg,b,pbicg%b,n)
  call parcel_sync(pbicg)
#endif
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

  call parcel_sync(pbicg) 
  call calc_res(pbicg,x,b,r)

  call parcel_sync(pbicg)    
  call calc_dot(pbicg,icomm,b,b,bnorm,n)
  parcel_int = 1
  call sqrt_val(pbicg,bnorm,parcel_int)
  call parcel_sync(pbicg)

  call calc_dot(pbicg,icomm,r,r,rnorm,n)
  call sqrt_val(pbicg,rnorm,parcel_int)
  call parcel_sync(pbicg)
  
  call get_val_FP64(pbicg,bnorm,bnorm_FP64,parcel_int)
  call get_val_FP64(pbicg,rnorm,rnorm_FP64,parcel_int)
  call parcel_sync(pbicg)

  if(myrank.eq.0)then
!     write(*,*)" rnorm,bnorm = ",rnorm_FP64,bnorm_FP64
  endif

  parcel_int = 1
  call set_val_FP64(pbicg,c1,0.0d0,parcel_int)
  call set_val_FP64(pbicg,c2,0.0d0,parcel_int)
  call set_val_FP64(pbicg,c3,0.0d0,parcel_int)

  ! r0=r
  call xcpy(pbicg,r,r0,n)

  ! p=r
  call xcpy(pbicg,r,p,n)

  ! c1=<r0,r0>
  call calc_dot(pbicg,icomm,r0,r0,c1,n)
     
  etime=MPI_Wtime()
  rtime(1) =  rtime(1) + etime-stime

  call parcel_sync(pbicg)
  do itr = 1,itrmax
     pbicg%reshistory(itr)=rnorm_FP64/bnorm_FP64
     if(myrank.eq.0)then
!        write(*,*)"@ res @",itr,rnorm_FP64,rnorm_FP64/bnorm_FP64
     endif

     !   solve Mpp=p
     stime=MPI_Wtime()
     call  calc_precon(pbicg,p,pp) 
     etime=MPI_Wtime()
     rtime(2) =  rtime(2) + etime-stime

     !   q=App
     stime=MPI_Wtime()
     call calc_Matvec(pbicg,pp,q)
     etime=MPI_Wtime()
     rtime(3) =  rtime(3) + etime-stime

     !   c2=<r0,q>
     stime=MPI_Wtime()
     call calc_dot(pbicg,icomm,r0,q,c2,n)
     etime=MPI_Wtime()
     rtime(4) =  rtime(4) + etime-stime
     
     call parcel_sync(pbicg)
     !   alpha=c1/c2 
     parcel_int = 1
     call get_val_FP64(pbicg,c1,c1_FP64,parcel_int)
     call get_val_FP64(pbicg,c2,c2_FP64,parcel_int)
     alpha_FP64=c1_FP64/c2_FP64
     call set_val_FP64(pbicg,alpha,alpha_FP64,parcel_int)
     call parcel_sync(pbicg)

     stime=MPI_Wtime()
     parcel_int = 1
     call set_val_FP64(pbicg,tmp1 , 1.0d0    ,parcel_int)
     call get_val_FP64(pbicg,alpha, tmp2_FP64,parcel_int)
     call set_val_FP64(pbicg,tmp2 ,-tmp2_FP64,parcel_int)
     call set_val_FP64(pbicg,tmp3 , 0.0d0    ,parcel_int)
     call parcel_sync(pbicg)

     !   e=r-alpha*q
     call axpy_vec3(pbicg,tmp1,tmp2,tmp3,r,q,e,n)
     etime=MPI_Wtime()
     rtime(5) =  rtime(5) + etime-stime

     !   solve Mpe=e
     stime=MPI_Wtime()
     call  calc_precon(pbicg,e,pe)
     etime=MPI_Wtime()
     rtime(6) =  rtime(6) + etime-stime

     !   v=Ape
     stime=MPI_Wtime()
     call calc_Matvec(pbicg,pe,v)
     etime=MPI_Wtime()
     rtime(7) =  rtime(7) + etime-stime

     !   c3=<e,v>/<v,v>
     ! ここは後で書き換え
     stime=MPI_Wtime()
     call BiCG_calc_dot2(pbicg,icomm,e,v,tmp1,tmp2)
     etime=MPI_Wtime()
     rtime(8) =  rtime(8) + etime-stime

     parcel_int = 1
     call get_val_FP64(pbicg, tmp1, tmp1_FP64,parcel_int)
     call get_val_FP64(pbicg, tmp2, tmp2_FP64,parcel_int)
     c3_FP64 = tmp1_FP64 / tmp2_FP64
     parcel_int = 1
     call set_val_FP64(pbicg,c3,c3_FP64,parcel_int)

     !   x=x+alpha*pp+c3*pe
     stime=MPI_Wtime()
     parcel_int = 1
     call set_val_FP64(pbicg,tmp1 , 1.0d0,parcel_int)
     call axpy_vec3(pbicg,alpha,c3,tmp1,pp,pe,x,n)
     etime=MPI_Wtime()
     rtime(9) =  rtime(9) + etime-stime

     !   r=e-c3*v
     !   c1=<r0,r>
     stime=MPI_Wtime()
     call BiCG_axpy3_dot(pbicg,icomm,r,e,v,r0,c3,c1,rnorm)
     etime=MPI_Wtime()
     rtime(10) =  rtime(10) + etime-stime

     parcel_int = 1
     call get_val_FP64(pbicg,rnorm,rnorm_FP64,parcel_int)
     rnorm_FP64 = sqrt(rnorm_FP64)
     call set_val_FP64(pbicg,rnorm,rnorm_FP64,parcel_int)

     if( (rnorm_FP64/bnorm_FP64) < rtolmax ) then
        exit
     endif
     parcel_int = 1
     call get_val_FP64(pbicg,c1,c1_FP64,parcel_int)
     call get_val_FP64(pbicg,c2,c2_FP64,parcel_int)
     call get_val_FP64(pbicg,c3,c3_FP64,parcel_int)
     beta_FP64 = c1_FP64/(c3_FP64*c2_FP64)
     call set_val_FP64(pbicg,beta,beta_FP64,parcel_int)

     stime=MPI_Wtime()

     parcel_int = 1
     call set_val_FP64(pbicg,tmp1,1.0d0,parcel_int)
     tmp2_FP64 = -beta_FP64*c3_FP64
     call set_val_FP64(pbicg,tmp2,tmp2_FP64,parcel_int)
     call axpy_vec3(pbicg,tmp1,tmp2,beta,r,q,p,n)
     etime=MPI_Wtime()
     rtime(11) =  rtime(11) + etime-stime

  enddo

#if 1
  parcel_int = 1
  call calc_res(pbicg,x,b,r)
  call calc_dot(pbicg,icomm,r,r,rnorm,n)
  call sqrt_val(pbicg,rnorm,parcel_int)
  call get_val_FP64(pbicg,rnorm,rnorm_FP64,parcel_int)

  if(myrank.eq.0)then
!     write(*,*)"@ true res @",rnorm_FP64,rnorm_FP64/bnorm_FP64
  endif
#endif

  call MPI_Barrier(icomm,ierr)
  etime0=MPI_Wtime()
#ifdef timer
 if(myrank.eq.0)then
     write(*,*)"pbicg         ",etime0-stime0
     write(*,*)"    init      ",rtime(1)
     write(*,*)"    Precon-1  ",rtime(2)
     write(*,*)"    MatVec-1  ",rtime(3)
     write(*,*)"    calc_dot  ",rtime(4)
     write(*,*)"    axpy3     ",rtime(5)
     write(*,*)"    Precon-2  ",rtime(6)
     write(*,*)"    MatVec-2  ",rtime(7)
     write(*,*)"    calc_dot2 ",rtime(8)
     write(*,*)"    axpy4     ",rtime(9)
     write(*,*)"    axpy3_dot ",rtime(10)
     write(*,*)"    axpy4_    ",rtime(11)
  endif
#endif

  call parcel_sync(pbicg)
#ifdef CUDA_SOLVER
  call cudaMemcpy_d2d_FP64(x%C_vec_FP64,pbicg%C_x,n,offset_zero,offset_zero)
#else
  call get_vec_FP64(pbicg,x,pbicg%x,n)
#endif
  call parcel_sync(pbicg)

  pbicg%niter=min(itrmax,itr)
  if(itrmax<=itr)then
     pbicg%iret=PARCEL_RETURNCODE_MAXITERATION
  else
     pbicg%iret=PARCEL_RETURNCODE_NORMAL
  endif

  call parcel_sync(pbicg)
  call free_vec(pbicg,r0)
  call free_vec(pbicg,p)
  call free_vec(pbicg,q)
  call free_vec(pbicg,pe)
  call free_vec(pbicg,v)
  call free_vec(pbicg,pp)
  call free_vec(pbicg,e)
  call free_vec(pbicg,x)
  call free_vec(pbicg,r)
  call free_vec(pbicg,alpha)
  call free_vec(pbicg,beta)
  call free_vec(pbicg,c1)
  call free_vec(pbicg,c2)
  call free_vec(pbicg,c3)
  call free_vec(pbicg,rnorm)
  call free_vec(pbicg,bnorm)
  call free_vec(pbicg,tmp1)
  call free_vec(pbicg,tmp2)
  call free_vec(pbicg,tmp3)
  call parcel_sync(pbicg)

end subroutine KSP_pBiCG
