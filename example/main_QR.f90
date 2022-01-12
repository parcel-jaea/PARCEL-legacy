
program main
  use mpi
  use PARCEL_integer
  use QR_kernel
  implicit none
  integer(kind=intType),parameter :: n=10000
  integer(kind=intType),parameter :: s=10

  real*8 A(n*s)
  real*8 A0(n*s)
  real*8 R(s*s)
  real*8 Q(n*s)
  integer myrank
  integer ierr
  integer i,j
  integer nblock,jend
  integer icomm
  call MPI_Init(ierr)  
  call MPI_Comm_rank(MPI_COMM_WORLD,myrank, ierr)

  A=0.0
  A0=0.0
  Q=0.0
  R=0.0
  do i=1,myrank+1
     call random_number(A)
  enddo
  A0 = A
  icomm = MPI_COMM_WORLD

 ! KSP_QRDECOMPOSITION_MGS=1        
 ! KSP_QRDECOMPOSITION_CGS=2        
 ! KSP_QRDECOMPOSITION_TSQR=3       
 ! KSP_QRDECOMPOSITION_CHOLESKYQR=4 
 ! KSP_QRDECOMPOSITION_CHOLESKYQR2=5

  do i = 1,5
     if(myrank.eq.0)then
        write(*,*)
        select case( i )
        case(1)
           write(*,*)"QR decomposition : MGS "
        case(2)
           write(*,*)"QR decomposition : CGS "
        case(3)
           write(*,*)"QR decomposition : TSQR "
        case(4)
           write(*,*)"QR decomposition : CHOLESKYQR "
        case(5)
           write(*,*)"QR decomposition : CHOLESKYQR2 "
        end select
     endif
     A = A0 
     Q = 1.0d0
     call parcel_dqr(n,s,icomm,A,Q,R,i)
     call check_QtQ(n,s,icomm,Q)
     call check_A_QR(n,s,icomm,A0,Q,R)
  enddo

  call MPI_Finalize(ierr)

end program main
