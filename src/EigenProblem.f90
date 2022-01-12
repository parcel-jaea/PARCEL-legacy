module EigenProblem
  implicit none

contains
  subroutine get_EigenPair_count_KSP(ncount,method)
    use krylov
    implicit none
    type(KSP) :: method
    integer ncount
    if( allocated(method%EigenValueReal) )then
       ncount = size( method%EigenValueReal )
    else
       ncount = -1
    endif

  end subroutine get_EigenPair_count_KSP

  subroutine get_EigenPair_KSP( &
       ncount,&
       EigenValueReal,EigenValueImag, &
       EigenVector, &
       EigenErrorNorm, &
       EmaxID,EminID,method)
    use krylov
    use vec_kernel
    implicit none
    integer ncount
    type(KSP) :: method
    real*8 EigenVector(method%prm%n*ncount) 
    real*8 EigenValueReal(ncount) 
    real*8 EigenValueImag(ncount) 
    real*8 EigenErrorNorm(ncount) 
    integer EmaxID,EminID
    integer i

    do i = 1,ncount
       EigenValueReal(i) = method%EigenValueReal(i)
       EigenValueImag(i) = method%EigenValueImag(i)
       EigenErrorNorm(i) = method%EigenErrorNorm(i)
    enddo
    call xcpy(method%prm%n*ncount,method%EigenVector,EigenVector)

    EmaxID = method%MaxEigenID
    EminID = method%MinEigenID

  end subroutine get_EigenPair_KSP
end module EigenProblem
