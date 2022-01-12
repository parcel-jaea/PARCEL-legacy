module power_method_kernel
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

end module power_method_kernel

! -----
subroutine power_method_symMatrix(method)
  use power_method_kernel
  use krylov
  use krylov_kernel
  use vec_kernel
  use PARCEL_integer
  use vec_ddm_kernel,only : set_scala_ddm
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
        call set_scala_ddm(method%prm%n,1.0d0,y,method%ddmMat%mat%prm)
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

end subroutine power_method_symMatrix

