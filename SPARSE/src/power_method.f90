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
     !$OMP PARALLEL DO PRIVATE(i)
     do i=1,n
        y(i)=1.0d0
     enddo
     call calc_dot(n,icomm,y,y,vnorm)
  endif
  tmp1=1.0d0/sqrt(vnorm)
  tmp2=0.0d0

  call axpy(n,tmp1,tmp2,y,x)
  do itr=1,itrmax
     call calc_precon(method,x,z)
     call calc_Matvec(method,z,y)
     call power_calc_dot2(n,icomm,y,x,vnorm,Emax)
!     write(*,*)"Eigen,",Emax
     tmp1=1.0d0/sqrt(vnorm)
     tmp2=0.0d0
     call axpy(n,tmp1,tmp2,y,x)
     ! 固有値はレイリー商で評価
  enddo

  method%EigenMax_Real=Emax
  method%EigenMax_Imag=0.0d0
  method%EigenMin_Real=0.0d0
  method%EigenMin_Imag=0.0d0

  deallocate(x)
  deallocate(z)
  deallocate(y)

end subroutine power_method_symMatrix

