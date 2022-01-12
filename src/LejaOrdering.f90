module leja_ordering_module
    implicit none
    
    type leja_complex_number
    double precision real2
    double precision img
    end type leja_complex_number
    
    double precision,allocatable::mu(:)
    double precision,allocatable::tmp(:)
    double precision,allocatable::ptmp(:)
    double precision,allocatable::ptmp_prod(:)
    type(leja_complex_number),allocatable::theta(:)
    integer,allocatable::outlist(:)
    
end module leja_ordering_module
        
subroutine init_leja_ordering(s)
    use leja_ordering_module
    implicit none
    
    integer,intent(in) :: s
    
    allocate( ptmp(0:s-1) )
    allocate( theta(0:s-1) )
    allocate( outlist(0:s-1) )
    allocate( mu(0:s-1) )
    allocate( tmp(0:s-1) )
    allocate( ptmp_prod(0:s-1) )
end subroutine init_leja_ordering

subroutine finish_leja_ordering()
    use leja_ordering_module
    
    deallocate(ptmp)
    deallocate(theta)
    deallocate(outlist)
    deallocate(mu)
    deallocate(tmp)
    deallocate(ptmp_prod)
end subroutine finish_leja_ordering
    
! Whatever sorting algorithm is used is should be stable, otherwise it must be ensured
! that complex conjugate pairs appear next to each other with the positive element before 
! the negative e.g. 10 + 2i, 10 -2i not 10 - 2i, 10 + 2i   
subroutine insertion_sort(array_length,eigen_real_part,eigen_img_part)
    implicit none
    integer,intent(in) :: array_length
    double precision :: eigen_real_part(0:array_length-1)
    double precision :: eigen_img_part(0:array_length-1)
    integer i,j
    double precision tr
    double precision ti
    
    double precision treal,timg
    
    do i=1,array_length-1
        j = i
        do while(j > 0)
            if (sqrt(eigen_real_part(j-1) * eigen_real_part(j-1) + eigen_img_part(j-1) * eigen_img_part(j-1) )&
		 < sqrt(eigen_real_part(j) * eigen_real_part(j) + eigen_img_part(j) * eigen_img_part(j) )) then
                tr = eigen_real_part(j)
                ti = eigen_img_part(j)
                eigen_real_part(j) = eigen_real_part(j - 1)
			    eigen_img_part(j) = eigen_img_part(j - 1)
			    eigen_real_part(j - 1) = tr
			    eigen_img_part(j - 1) = ti
                j = j - 1
            else
                exit
            end if
        end do
    end do
end subroutine insertion_sort
    
double precision function mleja_prod(scalk,L)
    use leja_ordering_module
    implicit none
    
    type(leja_complex_number) scalk
    
    integer L
    integer j
    
    double precision tr
    double precision ti
 
    do j=0,L-2
        tr = scalk%real2 - theta(j)%real2
        ti = scalk%img - theta(j)%img
        ptmp_prod(j) = sqrt(tr*tr+ti*ti) ** (mu(outlist(j)))
    end do
    
    mleja_prod = 1.0
    do j=0,L-2
        mleja_prod = mleja_prod * ptmp_prod(j) ! out*= ptmp_prod(j) 
    end do
end function mleja_prod

subroutine modified_leja_ordering(n,eigen_real_part,eigen_img_part)
    use leja_ordering_module
    implicit none
    
    interface
        subroutine insertion_sort(array_length,eigen_real_part,eigen_img_part)
        integer,intent(in) :: array_length
        double precision :: eigen_real_part(0:array_length-1)
        double precision :: eigen_img_part(0:array_length-1)
        end subroutine 
    end interface
    
    integer n
    double precision :: eigen_real_part(0:n-1)
    double precision :: eigen_img_part(0:n-1)
    
    double precision C
    integer L
    double precision Cold
	double precision tr
    double precision ti
    integer contained
    type(leja_complex_number) :: t
    double precision max
    
    double precision mleja_prod
    
	integer k
    integer i,j
    
    double precision t1,t2
    
    
!   1. sort the data
	call insertion_sort(n, eigen_real_part, eigen_img_part) ! stable sort the date in descending order

!   allocate
    call init_leja_ordering(n)
    
!	double *ptmp = malloc(n*sizeof(double));
	do i=0,n-1
        ptmp(i) = 1.0
    end do
    
!	leja_complex_number *theta = malloc(n*sizeof(leja_complex_number));
!	int *outlist = malloc(n*sizeof(int));
!	double *mu = malloc(n*sizeof(double));
!	double *tmp = malloc(n*sizeof(double));
    
	do i=0,n-1
        mu(i) = 1.0
    end do
	C = 1.0
	L = 0
    
!	 no need to use max because the array is already sorted in descending order
!	 so the first element must be the max
    
	if (eigen_img_part(0) > 0.0) then ! must be true if the number has a complex part because conjugate complex pairs appear next to each other
!		1. case a complex conjugate pair
		theta(0)%real2 = eigen_real_part(0)
		theta(0)%img = eigen_img_part(0)
		theta(1)%real2 = eigen_real_part(1)
		theta(1)%img = eigen_img_part(1)
		outlist(0) = 0
		outlist(1) = 1
		L = 2
	else
		theta(0)%real2 = eigen_real_part(0)
		theta(0)%img = 0.0
		outlist(0) = 0
		L = 1
	end if

	do while (L < n)
		Cold = C

		do j=0,L-2
!			double test = theta[L-1].real;
!			double test1 = theta[L-1].img;
!			double test2 = mu[outlist[j]] / L;
!			double test3 = theta[j].real;
            tr = theta(L-1)%real2 - theta(j)%real2
			ti = theta(L-1)%img - theta(j)%img
			ptmp(j + 1) = sqrt(tr*tr + ti*ti) ** (mu(outlist(j))/(L))
		end do
		
!		 prod
		C = 1.0
		do j=0,n-1
            C = C * ptmp(j)
        end do
        
!		 rescale all the shifts
		do j=0,n-1
			eigen_real_part(j) = eigen_real_part(j) / (C / Cold)
			eigen_img_part(j) = eigen_img_part(j) / (C / Cold)
			theta(j)%real2 = theta(j)%real2 / (C / Cold)
			theta(j)%img = theta(j)%img / (C / Cold)
        end do
        
!		memset(tmp, 0.0, n*sizeof(double));
		do i=0,n-1
            tmp(i) = 0.0
        end do
		do j=0,n-1
			contained = 0
			do i=0,L-1
				if (outlist(i) == j) then
                    contained = 1
                end if
			end do
			if (contained == 0) then
				
				t%real2 = eigen_real_part(j)
				t%img = eigen_img_part(j)
				tmp(j) = mleja_prod(t, L) ! tmp(j) = mleja_prod(t, theta, mu, outlist, L)
			end if
        end do
        
!	 max
		max = abs(tmp(0))
		k = 0
		do j=1,n-1
			if (abs(tmp(j)) > max) then
				max = abs(tmp(j))
				k = j
			end if
        end do
        
		if (abs(eigen_img_part(k)) > 0.0) then  ! must be true if the number has a complex part because conjugate complex pairs appear next to each other
			if (eigen_img_part(k) < 0.0) then
!				1. case a complex conjugate pair
				theta(L)%real2 = eigen_real_part(k-1)
				theta(L)%img = eigen_img_part(k-1)
				theta(L + 1)%real2 = eigen_real_part(k)
				theta(L + 1)%img = eigen_img_part(k)
				outlist(L) = k-1
				outlist(L + 1) = k
				L = L + 2
			else
				theta(L)%real2 = eigen_real_part(k)
				theta(L)%img = eigen_img_part(k)
				theta(L + 1)%real2 = eigen_real_part(k+1)
				theta(L + 1)%img = eigen_img_part(k+1)
				outlist(L) = k
				outlist(L + 1) = k+1
				L = L + 2
            end if
		else
			theta(L)%real2 = eigen_real_part(k)
			theta(L)%img = 0.0
			outlist(L) = k
			L = L + 1
		end if
	end do

!	 unscale the values
	do i=0,n-1
		theta(i)%real2 = theta(i)%real2 * C
		theta(i)%img   = theta(i)%img   * C
    end do
    
#ifdef debug    
    write(*,*) ("unscale the values")
    do i=0,n-1
        write(*,'(F12.4,F12.4)') (theta(i)%real2, theta(i)%img)
    end do 
#endif
    do i=0,n-1
        eigen_real_part(i) = theta(i)%real2
        eigen_img_part(i)  = theta(i)%img
    end do

!   free
    call finish_leja_ordering()
    
end subroutine

! program main
!     use leja_ordering_module
!     implicit none

!     integer n
!     double precision,allocatable :: eigen_real_part(:), eigen_img_part(:)
!     integer i    
!     n = 4

!     allocate( eigen_real_part(0:n-1) )
!     allocate( eigen_img_part(0:n-1)  )

! !!	eigen_real_part(0) =  0.7995
! !!	eigen_real_part(1) = -0.0994
! !!	eigen_real_part(2) = -0.0994
! !!	eigen_real_part(3) = -0.1007
! !!	eigen_img_part(0)  =  0.0000
! !!	eigen_img_part(1)  =  0.4008
! !!	eigen_img_part(2)  = -0.4008
! !!	eigen_img_part(3)  =  0.0000
    
!     eigen_real_part(0) = -0.1007
!     eigen_real_part(1) = -0.0994
!     eigen_real_part(2) = -0.0994
!     eigen_real_part(3) =  0.7995

!     eigen_img_part(0)  =  0.0000
!     eigen_img_part(1)  =  0.4008
!     eigen_img_part(2)  = -0.4008
!     eigen_img_part(3)  =  0.0000

!     do i=0,n-1
!        write(*,*)"@ org  @",i,eigen_real_part(i),eigen_img_part(i)
!     enddo
!     write(*,*)" -----" 
!     call modified_leja_ordering(n, eigen_real_part, eigen_img_part)
!     do i=0,n-1
!        write(*,*)"@ sort @",i,eigen_real_part(i),eigen_img_part(i)
!     enddo
    
!     deallocate(eigen_real_part)
!     deallocate(eigen_img_part)

! end program main
