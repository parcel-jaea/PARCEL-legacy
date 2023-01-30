module power_method_kernel_CPU
contains
  subroutine power_calc_dot2(n,icomm,x,y,xxnorm,xynorm)
    use mpi
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer icomm
    real*8 x(n)
    real*8 y(n) 
    real*8 xxnorm,xynorm
    real*8 sbuf(2),rbuf(2)
    integer(kind=intType) i
    integer ierr
    xxnorm=0.0d0
    xynorm=0.0d0
    !$OMP PARALLEL DO PRIVATE(i) REDUCTION(+:xxnorm,xynorm) 
    do i=1,n
       xxnorm = xxnorm+x(i)*x(i)
       xynorm = xynorm+x(i)*y(i)
    enddo
    !$OMP END PARALLEL DO 
    sbuf(1)=xxnorm
    sbuf(2)=xynorm
    rbuf=0.0d0
    call MPI_Allreduce(sbuf,rbuf,2,MPI_REAL8,MPI_SUM,icomm,ierr)    
    xxnorm=rbuf(1)
    xynorm=rbuf(2)

  end subroutine power_calc_dot2

subroutine power_method_symMatrix_CPU(method)
  use krylov
  use krylov_kernel
  use vec_kernel
  use vec_ddm_kernel

  use PARCEL_integer
!  use vec_ddm_kernel,only : set_scala_ddm
  implicit none
  type(KSP) :: method
  integer(kind=intType) n
  integer icomm
  integer myrank ! プロセスランク
  integer itrmax

  integer(kind=intType) i
  integer itr
  real*8,allocatable :: x(:)
  real*8,allocatable :: z(:)
  real*8,allocatable :: y(:)
  real*8 Emax
  real*8 vnorm
  real*8 tmp1,tmp2
  real*8 Ax_norm,kx_norm,errnorm
  Emax = 0.0d0
  vnorm = 0.0d0
  myrank=method%commVec%myrank
  icomm=method%commVec%comm
  n=method%prm%n
  allocate(x(n))
  allocate(z(n))
  allocate(y(n))
  call zero_initialize(n,x)
  call zero_initialize(n,z)
  call zero_initialize(n,y)

  itrmax=method%power_method_itrmax
  ! ------
  ! べき乗法
  ! x=x/||x||
  ! do itr=1,itrmax
  !   y = AM^-1x
  !   Eigen = <x|A|x>/<x|x>
  !   x = y / ||y||
  ! enddo
  ! ------
  call xcpy(n,method%x,y)
  call calc_dot(n,icomm,y,y,vnorm)

  if(vnorm.eq.0.0d0)then
     if( method%prm%MatFormat .eq. MATRIX_FORMAT_DDM )then
        call set_scala_ddm_CPU(method%prm%n,1.0d0,y,method%ddmMat%prm)
     else
        !$OMP PARALLEL DO PRIVATE(i)
        do i=1,n
           y(i)=1.0d0
        enddo
     endif
     call calc_dot(n,icomm,y,y,vnorm)        
  endif
  ! tmp1=1.0d0/sqrt(vnorm)
  ! tmp2=0.0d0
  ! call axpy(n,tmp1,tmp2,y,x)

  do itr=1,itrmax
     if(myrank.eq.0)then
        write(*,*)itr,"vnorm = ",vnorm," Emax = ",Emax
     endif

     tmp1=1.0d0/sqrt(vnorm)
     tmp2=0.0d0
     call axpy(n,tmp1,tmp2,y,x)
     call calc_precon(method,x,z)
     call calc_Matvec(method,z,y)
     call power_calc_dot2(n,icomm,y,x,vnorm,Emax)
     !     write(*,*)"Eigen,",Emax
     ! 固有値はレイリー商で評価
  enddo
  ! 固有ベクトルの規格化
  tmp1=1.0d0/sqrt(vnorm)
  tmp2=0.0d0
  call axpy(n,tmp1,tmp2,y,x)


  ! ||Ax-kx||/||kx|| を評価して固有ベクトルの精度を確認する
  call calc_Matvec(method,x,y)
  call calc_dot(n,icomm,y,y,Ax_norm)
  call axpy(n,Emax,0.0d0,x,z)
  call calc_dot(n,icomm,z,z,kx_norm)

  ! if(myrank.eq.0)then
  !    write(*,*)" Ax       , kx "
  !    call view_vec_diff(n,y,z)
  ! endif

  call axpy(n,1.0d0,-1.0d0,y,z)
  call calc_dot(n,icomm,z,z,errnorm)
  Ax_norm = sqrt(abs(Ax_norm))
  errnorm = sqrt(abs(errnorm))
  kx_norm = sqrt(abs(kx_norm))

#if 1
  if(myrank.eq.0)then

     write(*,*)
     write(*,*)"Power method iteration : ",itr
     write(*,*)"Power method eigen error "
     write(*,*)"  Emax     ||Ax-kx||/||kx||    "
     write(*,*)
     write(*,*)Emax,errnorm/kx_norm
     ! write(*,*)"     ||Ax||             = ",Ax_norm
     ! write(*,*)"     ||kx||             = ",kx_norm
     ! write(*,*)"     ||Ax-kx||          = ",errnorm
     ! write(*,*)"     ||Ax-kx|| / ||kx|| = ",errnorm/kx_norm
     ! endif
     ! if(myrank.eq.0)then
     write(*,*)" ----------------------- "
  endif
#endif

  if( allocated(method%EigenValueReal) ) deallocate(method%EigenValueReal)
  if( allocated(method%EigenValueImag) ) deallocate(method%EigenValueImag)
  if( allocated(method%EigenErrorNorm) ) deallocate(method%EigenErrorNorm)
  allocate(method%EigenValueReal(1))
  allocate(method%EigenValueImag(1))
  allocate(method%EigenErrorNorm(1))
  method%EigenValueReal(1) = Emax
  method%EigenValueImag(1) = 0.0d0
  method%EigenErrorNorm(1) = errnorm/kx_norm

  if( allocated(method%EigenVector) )then
     deallocate(method%EigenVector)
  endif
  allocate(method%EigenVector(n))
  call xcpy(n,x,method%EigenVector)

  method%EigenMax_Real = Emax
  method%EigenMax_Imag = 0.0d0
  method%EigenMin_Real = 0.0d0
  method%EigenMin_Imag = 0.0d0
  method%MaxEigenID = 1
  method%MinEigenID = 1

  deallocate(x)
  deallocate(z)
  deallocate(y)

end subroutine power_method_symMatrix_CPU

end module power_method_kernel_CPU

! module  power_method_kernel
!   implicit none
! contains

#if 1
subroutine power_method_symMatrix( method &
     )
  use krylov,only : KSP
  use PARCEL_define,only : PARCEL_CPU_FLAG_ON,PARCEL_GPU_FLAG_ON
  use KSP_STATUS,only : &
       STATUS_FLAG, &
       init_staus_flag,set_machineflag, &
       set_dataflag,set_calcflag
  use PARCEL_KRYLOV_MATRIX,only : calc_Matvec,calc_res
  use PARCEL_KRYLOV_PRECON,only : calc_precon
  use PARCEL_VEC_kernel
  use PARCEL_integer,only : intType
  implicit none

  type(KSP) :: method
  integer(kind=intType) n
  integer icomm
  integer itrmax ! 最大反復回数
  integer myrank
  type(STATUS_FLAG) :: flag
  type(PARCEL_VEC) :: x
  type(PARCEL_VEC) :: y
  type(PARCEL_VEC) :: z

  type(PARCEL_VEC) :: Emax
  type(PARCEL_VEC) :: vnorm
  type(PARCEL_VEC) :: tmp1,tmp2
  type(PARCEL_VEC) :: Ax_norm,kx_norm,errnorm

  real*8 Emax_FP64
  real*8 vnorm_FP64
  real*8 tmp1_FP64,tmp2_FP64
  real*8 Ax_norm_FP64,kx_norm_FP64,errnorm_FP64

  integer(kind=intType) vecLen
  integer(kind=intType) nColL
  integer(kind=intType) nColU
  integer itr
  integer(kind=intType) parcel_int
  integer(kind=intType),parameter :: offset_zero=0

#if 0
  method%EigenMax_Real = 1.0d0
  method%EigenMax_Imag = 0.0d0
  method%EigenMin_Real = 0.0d0
  method%EigenMin_Imag = 0.0d0
  method%MaxEigenID = 1
  method%MinEigenID = 1
  return 
#endif

  flag = method%flag
  itrmax=method%power_method_itrmax
  icomm     = method%commVec%comm
  myrank    = method%commVec%myrank
  n         = method%prm%n

  call parcel_sync(method)
  nColL = 0
  nColU = 0
  if(method%prm%MatFormat.eq.MATRIX_FORMAT_CRS)then
     nColL = method%crsMat%nColL
     nColU = method%crsMat%nColU
  endif
  vecLen = n+nColL+nColU

  call init_vec(method,x,vecLen)
  call init_vec(method,y,vecLen)
  call init_vec(method,z,vecLen)

  parcel_int = 1
  call init_vec(method,Emax,parcel_int)
  call init_vec(method,vnorm,parcel_int)
  call init_vec(method,tmp1,parcel_int)
  call init_vec(method,tmp2,parcel_int)
  call init_vec(method,Ax_norm,parcel_int)
  call init_vec(method,kx_norm,parcel_int)
  call init_vec(method,errnorm,parcel_int)

  parcel_int = 1
  call set_val_FP64(method,Emax,0.0d0,parcel_int)
  call set_val_FP64(method,vnorm,0.0d0,parcel_int)
  call set_val_FP64(method,tmp1,0.0d0,parcel_int)
  call set_val_FP64(method,tmp2,0.0d0,parcel_int)
  call set_val_FP64(method,Ax_norm,0.0d0,parcel_int)
  call set_val_FP64(method,kx_norm,0.0d0,parcel_int)
  call set_val_FP64(method,errnorm,0.0d0,parcel_int)

  parcel_int = 1
  call set_vec_val(method,x,tmp1,parcel_int,vecLen)
  call set_vec_val(method,y,tmp1,parcel_int,vecLen)
  call set_vec_val(method,z,tmp1,parcel_int,vecLen)

  ! ------
  ! べき乗法
  ! x=x/||x||
  ! do itr=1,itrmax
  !   y = AM^-1x
  !   Eigen = <x|A|x>/<x|x>
  !   x = y / ||y||
  ! enddo
  ! ------

#ifdef CUDA_SOLVER
  call cudaMemcpy_d2d_FP64(y%C_vec_FP64,method%C_x,n,offset_zero,offset_zero)
#else
  call set_vec_FP64(method,y,method%x,n)
#endif

  call calc_dot(method,icomm,y,y,vnorm,n)
  parcel_int = 1
  call get_val_FP64(method,vnorm,vnorm_FP64,parcel_int)

  if(vnorm_FP64.eq.0.0d0)then
     parcel_int = 1
     call set_val_FP64(method,tmp1,1.0d0,parcel_int)
     if( method%prm%MatFormat .eq. MATRIX_FORMAT_DDM )then
        call set_scala_ddm(method,tmp1,y,n)
     else
        call set_vec_val(method,y,tmp1,parcel_int,vecLen)
     endif
     call calc_dot(method,icomm,y,y,vnorm,n)
     ! if(myrank.eq.0)then
     !    write(*,*)"vnorm = ",vnorm_FP64
     ! endif
  endif

  do itr=1,itrmax
     parcel_int = 1
     call get_val_FP64(method,vnorm,vnorm_FP64,parcel_int)
     call get_val_FP64(method,Emax ,Emax_FP64,parcel_int)
      ! if(myrank.eq.0)then
      !    write(*,*)itr,"vnorm = ",vnorm_FP64," Emax = ",Emax_FP64
      ! endif
     tmp1_FP64 = 1.0d0/sqrt(vnorm_FP64)
     tmp2_FP64 = 0.0d0
     parcel_int = 1
     call set_val_FP64(method,tmp1,tmp1_FP64,parcel_int)
     call set_val_FP64(method,tmp2,tmp2_FP64,parcel_int)
     call axpy(method,tmp1,tmp2,y,x,n)
     call calc_precon(method,x,z)
     call calc_Matvec(method,z,y)

     !     call power_calc_dot2(n,icomm,y,x,vnorm,Emax)
     call calc_dot(method,icomm,y,y,vnorm,n)
     call calc_dot(method,icomm,x,y,Emax,n)

     !     write(*,*)"Eigen,",Emax
     ! 固有値はレイリー商で評価
  enddo
  ! 固有ベクトルの規格化
  parcel_int = 1
  call get_val_FP64(method,vnorm,vnorm_FP64,parcel_int)
  tmp1_FP64 = 1.0d0/sqrt(vnorm_FP64)
  tmp2_FP64 = 0.0d0
  call set_val_FP64(method,tmp1,tmp1_FP64,parcel_int)
  call set_val_FP64(method,tmp2,tmp2_FP64,parcel_int)

  call axpy(method,tmp1,tmp2,y,x,n)


  ! ||Ax-kx||/||kx|| を評価して固有ベクトルの精度を確認する
  call calc_Matvec(method,x,y)
  call calc_dot(method,icomm,y,y,Ax_norm,n)
  tmp1_FP64 = 0.0d0

  parcel_int = 1
  call set_val_FP64(method,tmp1,tmp1_FP64,parcel_int)
  call axpy(method,Emax,tmp1,x,z,n)
  call calc_dot(method,icomm,z,z,kx_norm,n)

  ! if(myrank.eq.0)then
  !    write(*,*)" Ax       , kx "
  !    call view_vec_diff(n,y,z)
  ! endif

  tmp1_FP64 =  1.0d0
  tmp2_FP64 = -1.0d0
  parcel_int = 1
  call set_val_FP64(method,tmp1,tmp1_FP64,parcel_int)
  call set_val_FP64(method,tmp2,tmp2_FP64,parcel_int)
  call axpy(method,tmp1,tmp2,y,z,n)

  call calc_dot(method,icomm,z,z,errnorm,n)

  parcel_int = 1
  call get_val_FP64(method,Emax,Emax_FP64,parcel_int)
  call get_val_FP64(method,Ax_norm,Ax_norm_FP64,parcel_int)
  call get_val_FP64(method,errnorm,errnorm_FP64,parcel_int)
  call get_val_FP64(method,kx_norm,kx_norm_FP64,parcel_int)
  Ax_norm_FP64 = sqrt(abs(Ax_norm_FP64))
  errnorm_FP64 = sqrt(abs(errnorm_FP64))
  kx_norm_FP64 = sqrt(abs(kx_norm_FP64))

#if 0
  if(myrank.eq.0)then

     write(*,*)
     write(*,*)"Power method iteration : ",itr
     write(*,*)"Power method eigen error "
     write(*,*)"  Emax     ||Ax-kx||/||kx||    "
     write(*,*)
     write(*,*)Emax_FP64,errnorm_FP64/kx_norm_FP64
     ! write(*,*)"     ||Ax||             = ",Ax_norm
     ! write(*,*)"     ||kx||             = ",kx_norm
     ! write(*,*)"     ||Ax-kx||          = ",errnorm
     ! write(*,*)"     ||Ax-kx|| / ||kx|| = ",errnorm/kx_norm
     ! endif
     ! if(myrank.eq.0)then
     write(*,*)" ----------------------- "
  endif
#endif

  if( allocated(method%EigenValueReal) ) deallocate(method%EigenValueReal)
  if( allocated(method%EigenValueImag) ) deallocate(method%EigenValueImag)
  if( allocated(method%EigenErrorNorm) ) deallocate(method%EigenErrorNorm)
  allocate(method%EigenValueReal(1))
  allocate(method%EigenValueImag(1))
  allocate(method%EigenErrorNorm(1))
  method%EigenValueReal(1) = Emax_FP64
  method%EigenValueImag(1) = 0.0d0
  method%EigenErrorNorm(1) = errnorm_FP64/kx_norm_FP64

  if( allocated(method%EigenVector) )then
     deallocate(method%EigenVector)
  endif
  allocate(method%EigenVector(n))
  call get_vec_FP64(method,x,method%EigenVector,n)

  method%EigenMax_Real = Emax_FP64
  method%EigenMax_Imag = 0.0d0
  method%EigenMin_Real = 0.0d0
  method%EigenMin_Imag = 0.0d0
  method%MaxEigenID = 1
  method%MinEigenID = 1

  ! ---------------------

  call free_vec(method,x)
  call free_vec(method,y)
  call free_vec(method,z)
  call free_vec(method,Emax)
  call free_vec(method,vnorm)
  call free_vec(method,tmp1)
  call free_vec(method,tmp2)
  call free_vec(method,Ax_norm)
  call free_vec(method,kx_norm)
  call free_vec(method,errnorm)

end subroutine power_method_symMatrix

#else

subroutine power_method_symMatrix(method)
  use krylov,only : KSP
  use KSP_STATUS,only : &
       STATUS_FLAG,PARCEL_CPU_FLAG_ON,PARCEL_GPU_FLAG_ON
  use power_method_kernel_CPU,only : &
       power_method_symMatrix_CPU

  implicit none
  type(KSP) :: method
  type(STATUS_FLAG) :: flag

  flag = method%flag
  select case(flag%machineflag)
  case(PARCEL_CPU_FLAG_ON)
     call power_method_symMatrix_CPU(method)
  case(PARCEL_GPU_FLAG_ON)
     call power_method_symMatrix_CPU(method)
  end select

end subroutine power_method_symMatrix
#endif
