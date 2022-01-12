! 
program main
  use mpi
  use QR_kernel
  use PARCEL_define
  use wrap_DGEMM
  implicit none
  integer,parameter :: n = 5
  real*8  HesMatrix(n,n)
  real*8  HesMatrix0(n,n)
  real*8  EigenValue(n)
  real*8  Rmat(n,n)
  real*8  Rmat_MGS(n,n)
  real*8  Rmat_CGS(n,n)
  real*8  Rmat_TSQR(n,n)
  real*8  Qmat_MGS(n,n)
  real*8  Qmat_CGS(n,n)
  real*8  Qmat_TSQR(n,n)
  real*8  QRmat_MGS(n,n)
  real*8  QRmat_CGS(n,n)
  real*8  QRmat_TSQR(n,n)

  real*8  Qmat(n,n)
  real*8  Qtmat(n,n)
  real*8  QHmat(n,n)
  real*8  QHQmat(n,n)
  real*8  RQmat(n,n)
  real*8  QRmat(n,n)
  real*8 tmp
  integer i,j
  integer ierr

  call MPI_Init( ierr )
  HesMatrix = 0.0d0
  Rmat = 0.0d0
  Qmat = 0.0d0
  QRmat =0.0d0
  RQmat =0.0d0
  QHmat = 0.0d0
  QHQmat = 0.0d0
  Qtmat = 0.0d0
  HesMatrix0 = 0.0d0
  call random_number(HesMatrix)
  call random_number(EigenValue)

  do j = 1 ,n
     do i = 1,j-2
        HesMatrix(j,i) = 0.0d0
     enddo
  enddo
  HesMatrix0 = HesMatrix

  do j = 1 ,n
     write(*,*)j,"Eigen ",EigenValue(j)
     HesMatrix(j,j) =  HesMatrix(j,j) - EigenValue(j)
  enddo
  write(*,*)
!  call QRDecomposition(n,n,MPI_COMM_SELF,HesMatrix,Qmat,Rmat,KSP_QRDECOMPOSITION_MGS,-1,1)
!  call QRDecomposition(n,n,MPI_COMM_SELF,HesMatrix,Qmat,Rmat,KSP_QRDECOMPOSITION_MGS,KSP_QRDECOMPOSITION_TRUE_ORTHOGONAL_MATRIX,1)
!  call QRDecomposition(n,n,MPI_COMM_SELF,HesMatrix,Qmat,Rmat,KSP_QRDECOMPOSITION_CGS,KSP_QRDECOMPOSITION_TRUE_ORTHOGONAL_MATRIX,1)
!  call QRDecomposition(n,n,MPI_COMM_SELF,HesMatrix,Qmat,Rmat,KSP_QRDECOMPOSITION_TSQR,KSP_QRDECOMPOSITION_TRUE_ORTHOGONAL_MATRIX,1)
!  call QRDecomposition(n,n,MPI_COMM_SELF,HesMatrix,Qmat,Rmat,KSP_QRDECOMPOSITION_CHOLESKYQR,KSP_QRDECOMPOSITION_TRUE_ORTHOGONAL_MATRIX,1)


  do j = 1 ,n
     do i = 1,n
        Qtmat(j,i) = Qmat(i,j) 
     enddo
  enddo

  call wrap_BASIC_DGEMM(n,n,n,Rmat,Qmat,RQmat)
  do j = 1 ,n
     RQmat(j,j) =  RQmat(j,j) + EigenValue(j)
  enddo

  call wrap_BASIC_DGEMM(n,n,n,Qtmat,HesMatrix,QHmat)
  call wrap_BASIC_DGEMM(n,n,n,QHmat,Qmat,QHQmat)

!  call wrap_BASIC_DGEMM(n,n,n,Qtmat,Qmat,QHQmat)
!   call wrap_BASIC_DGEMM(n,n,n,Qmat,Rmat,QRmat)
  call wrap_BASIC_DGEMM(n,n,n,Rmat,Qmat,RQmat)

  do j = 1 ,n
     do i = 1 ,n
!       write(*,*) RQmat(i,j)
!       write(*,*) Rmat(i,j)
!        write(*,*) RQmat(i,j) ,  QHQmat(i,j)
!        write(*,*) QRmat(i,j) ,  HesMatrix(i,j),  HesMatrix0(i,j)
     enddo
     write(*,*)" -------------- "
     write(*,*)
  enddo

  Rmat_MGS = 0.0d0
  Rmat_CGS = 0.0d0
  Rmat_TSQR = 0.0d0
  Qmat_MGS = 0.0d0
  Qmat_CGS = 0.0d0
  Qmat_TSQR = 0.0d0

  call random_number(HesMatrix0)
  HesMatrix=HesMatrix0
  call QRDecomposition(n,n,MPI_COMM_SELF,HesMatrix,Qmat_MGS,Rmat_MGS,KSP_QRDECOMPOSITION_MGS,KSP_QRDECOMPOSITION_TRUE_ORTHOGONAL_MATRIX,1)

  HesMatrix=HesMatrix0
  call QRDecomposition(n,n,MPI_COMM_SELF,HesMatrix,Qmat_CGS,Rmat_CGS,KSP_QRDECOMPOSITION_CGS,KSP_QRDECOMPOSITION_TRUE_ORTHOGONAL_MATRIX,1)
  HesMatrix=HesMatrix0
  call QRDecomposition(n,n,MPI_COMM_SELF,HesMatrix,Qmat_TSQR,Rmat_TSQR,KSP_QRDECOMPOSITION_TSQR,KSP_QRDECOMPOSITION_TRUE_ORTHOGONAL_MATRIX,1)
  HesMatrix=HesMatrix0

  call wrap_BASIC_DGEMM(n,n,n,Qmat_MGS,Rmat_MGS,QRmat_MGS)
  call wrap_BASIC_DGEMM(n,n,n,Qmat_CGS,Rmat_CGS,QRmat_CGS)
  call wrap_BASIC_DGEMM(n,n,n,Qmat_TSQR,Rmat_TSQR,QRmat_TSQR)

  do j = 1 ,n
     do i = 1 ,n
       write(*,*) Rmat_MGS(i,j),Rmat_CGS(i,j),Rmat_TSQR(i,j)
     enddo
     write(*,*)" -------------- "
     write(*,*)
  enddo
     write(*,*)" -------------- "
     write(*,*)" -------------- "
     write(*,*)" -------------- "
     write(*,*)" -------------- "
     write(*,*)" -------------- "
  do j = 1 ,n
     do i = 1 ,n
       write(*,*) Qmat_MGS(i,j),Qmat_CGS(i,j),Qmat_TSQR(i,j)
     enddo
     write(*,*)" -------------- "
     write(*,*)
  enddo
     write(*,*)" -------------- "
     write(*,*)" -------------- "
     write(*,*)" -------------- "
     write(*,*)" -------------- "
     write(*,*)" -------------- "
  do j = 1 ,n
     do i = 1 ,n
       write(*,*) QRmat_MGS(i,j)/HesMatrix(i,j),QRmat_CGS(i,j)/HesMatrix(i,j),QRmat_TSQR(i,j)/HesMatrix(i,j)
     enddo
     write(*,*)" -------------- "
     write(*,*)
  enddo

  call MPI_Finalize(ierr)

end program main
