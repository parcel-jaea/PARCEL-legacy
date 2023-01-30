module moore_penrose_kernel
contains
subroutine moore_penrose(A,Ainv,m,n)
  implicit none
  integer m,n
  real*8 A(m,n),Ainv(n,m)

  integer i,j
  real*8 sigma(m*n),sigma_inv(n*m)
  real*8 Vt_sigma(n*m)
  real*8 work(m*n*5)
  integer lwork
  integer info

  integer mn_min,mn_max
!  real*8 U(m,n),Vt(n,m)
  real*8 U(m,m),Vt(n,n)

  LWORK=n*m*5
  call DGESVD('A', 'A', &
       m, n, A,m,&
       sigma, &
       U, m, &
       Vt, n, &
       work, lwork, info)

   ! write(*,*)"-------"
   ! write(*,*)"U"
   ! do i=1,m
   !    do j=1,m
   !       write(*,*)U(j,i)
   !    enddo
   !    write(*,*)"---"
   ! enddo

   ! write(*,*)"-------"
   ! write(*,*)"sigma"
   ! do j=1,min(m,n)
   !    write(*,*)sigma(j)
   ! enddo

   ! write(*,*)"-------"
   ! write(*,*)"Vt"
   ! do i=1,n
   !    do j=1,n
   !       write(*,*)Vt(j,i)
   !    enddo
   !    write(*,*)"---"
   ! enddo

  Vt_sigma=0.0d0
  sigma_inv=0.0d0
  mn_min=min(m,n)
  mn_max=max(m,n)
  do j=1,min(m,n)
     if(sigma(j)>0.0d0)then
!        sigma_inv(j,j)=1.0d0/sigma(j)
        sigma_inv(j+(j-1)*n)=1.0d0/sigma(j)
     else
!        sigma_inv(j,j)=0.0d0
        sigma_inv(j+(j-1)*n)=0.0d0
     endif
  enddo
  
  CALL DGEMM('T','N',n,m,n,1.0d0,Vt,n,sigma_inv,n,0.0d0,Vt_sigma,n)
  CALL DGEMM('N','T',n,m,m,1.0d0,Vt_sigma,n,U,m,0.0d0,Ainv,n)
  return

! ----
  ! write(*,*)"-------"
  ! write(*,*)"V^t sigma"
  ! do j=1,3
  !    write(*,*)Ainv(j,1),Ainv(j,2)
  ! enddo

end subroutine moore_Penrose
end module moore_penrose_kernel
! program main
!   implicit none
! ! 特異値計算
! ! A = V D U^T
! ! A^+= V D^+ U^T
! !!  integer,parameter :: n=2,m=3
!   integer,parameter :: n=3,m=2
!   real*8 A(m,n),Atmp(m,n)
!   real*8 A_inv(n,m)
!   real*8 AAinv(m,m)
!   real*8 AinvA(n,n)

!   real*8 AAinvA(m,n)   ! = A
!   real*8 AinvAAinv(n,m)! = A_inv
!   integer i,j

!   A=0.0d0
!   A_inv=0.0d0

!   call random_number(A)

! ! ----
!    ! A(1,1)=2.0d0
!    ! A(2,1)=1.0d0

!    ! A(1,2)=1.0d0
!    ! A(2,2)=-1.0d0

!    ! A(1,3)=2.0d0
!    ! A(2,3)=1.0d0

!    Atmp=A
! !   do j=1,3
!    do j=1,2
!       write(*,*)A(j,1),A(j,2),A(j,3)
!    enddo
!    call Moore_Penrose(Atmp,A_inv,m,n)

!    write(*,*)"-------"
! !   do j=1,2
!    do j=1,3
!       write(*,*)A_inv(j,1),A_inv(j,2)!,A_inv(j,3)
!    enddo

!    write(*,*)"-------"
!    CALL DGEMM('N','N', m, m, n, 1.0d0, A, m,  A_inv, n, 0.0d0, AAinv, m)
!    CALL DGEMM('N','N', m, n, m, 1.0d0, AAinv, m,  A, m, 0.0d0, AAinvA, m)

!    CALL DGEMM('N','N', n, n, m, 1.0d0, A_inv, n,  A, m, 0.0d0, AinvA, n)
!    CALL DGEMM('N','N',n,m,n,1.0d0,AinvA,n,A_inv,n,0.0d0,AinvAAinv,n)


!    ! write(*,*)"------"
!    ! write(*,*)" E=A A^-1 "
!    ! do i=1,m
!    !    do j=1,m
!    !       write(*,*)AAinv(j,i)
!    !    enddo
!    !    write(*,*)"------"
!    ! enddo

!    write(*,*)"-------"
! !   do j=1,3
!    do j=1,2
!       write(*,*)AAinvA(j,1),AAinvA(j,2),AAinvA(j,3)
!    enddo
!    write(*,*)"-------"
!    do j=1,3
! !   do j=1,2
!       write(*,*)AinvAAinv(j,1),AinvAAinv(j,2)!,AinvAAinv(j,3)
!    enddo

! end program main
