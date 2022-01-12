module MGrid_2D
  implicit none
  type MG
     real*8 ,allocatable :: A(:,:,:)
     real*8 ,allocatable :: x(:,:)
     real*8 ,allocatable :: r(:,:)
     real*8 ,allocatable :: b(:,:)
     real*8 ,allocatable :: d(:,:)
     integer nx,ny
  end type MG
contains
  subroutine make_MGrid(&
    grid, & 
    lv,&
    nb,nx0,ny0 &
    )
    implicit none
    type(MG) grid(lv)
    integer lv
    integer nb,nx0,ny0
    integer nx,ny
    integer i

    do i=1,lv
       nx=CEILING(nx0/(2.0d0**(i-1)) )
       ny=CEILING(ny0/(2.0d0**(i-1)) )
       allocate(grid(i)%A(1-nb:nx+nb,1-nb:ny+nb,5))
       allocate(grid(i)%x(1-nb:nx+nb,1-nb:ny+nb))
       allocate(grid(i)%b(1-nb:nx+nb,1-nb:ny+nb))
       allocate(grid(i)%r(1-nb:nx+nb,1-nb:ny+nb))
       allocate(grid(i)%d(1-nb:nx+nb,1-nb:ny+nb))
    enddo

  end subroutine make_MGrid

end module MGrid_2D

module poisson2D
  implicit none
  real*8,parameter ::  Lx=20.0d0
  real*8,parameter ::  Ly=20.0d0

  integer,parameter :: nx0=20
  integer,parameter :: ny0=20
  integer,parameter :: nb=1

  real*8,parameter ::  dx0=Lx/nx0
  real*8,parameter ::  dy0=Ly/ny0

  integer nx
  integer ny

  real*8 dx
  real*8 dy

  ! real*8 A(1-nb:nx+nb,1-nb:ny+nb,5)
  ! real*8 x(1-nb:nx+nb,1-nb:ny+nb)
  ! real*8 b(1-nb:nx+nb,1-nb:ny+nb)
  ! real*8 r(1-nb:nx+nb,1-nb:ny+nb)

  real*8 ,allocatable :: A(:,:,:)
  real*8 ,allocatable :: x(:,:)
  real*8 ,allocatable :: b(:,:)
  real*8 ,allocatable :: r(:,:)

contains

  subroutine make_matrix()
    implicit none
    integer i,j
    real*8 dxi,dyi
    real*8 dxi2,dyi2
    dxi=1.0d0/dx
    dyi=1.0d0/dy
    dxi2=dxi*dxi
    dyi2=dyi*dyi
    do j=1,ny
       do i=1,nx
          A(i,j,1)=dyi2
          A(i,j,2)=dxi2
          A(i,j,3)=-2.0*(dxi2+dyi2)
          A(i,j,4)=dxi2
          A(i,j,5)=dyi2

          if(i.eq.1)then
             A(i,j,2)=0.0
          endif
          if(i.eq.nx)then
             A(i,j,4)=0.0
          endif

          if(j.eq.1)then
             A(i,j,1)=0.0
          endif
          if(j.eq.ny)then
             A(i,j,5)=0.0
          endif

       enddo
    enddo
  end subroutine make_matrix

  subroutine make_vec()
    implicit none
    integer i,j
    do j=1,ny
       do i=1,nx
          x(i,j)=0.0d0
!          b(i,j)=-1.0d0 ! 一様分布
          b(i,j)=-1.0d0 ! 一様分布
          r(i,j)=0.0d0
       enddo
    enddo
  end subroutine make_vec

  subroutine set_param(nx_i,ny_i)
    implicit none
    integer nx_i,ny_i
    nx=nx_i
    ny=ny_i
    dx=Lx/nx
    dy=Ly/ny
    allocate(A(1-nb:nx+nb,1-nb:ny+nb,5))
    allocate(x(1-nb:nx+nb,1-nb:ny+nb))
    allocate(b(1-nb:nx+nb,1-nb:ny+nb))
    allocate(r(1-nb:nx+nb,1-nb:ny+nb))
  end subroutine set_param

  subroutine free_mem()
    implicit none
    deallocate(A)
    deallocate(x)
    deallocate(b)
    deallocate(r)
  end subroutine free_mem

  subroutine initialize(nx_i,ny_i)
    implicit none
    integer nx_i,ny_i
    call set_param(nx_i,ny_i)
    A=0.0d0
    x=0.0d0
    b=0.0d0
    r=0.0d0
    call make_matrix
    call make_vec
  end subroutine initialize

  subroutine SpMV(x1,y1)
    implicit none
    real*8 x1(1-nb:nx+nb,1-nb:ny+nb)
    real*8 y1(1-nb:nx+nb,1-nb:ny+nb)
    integer i,j
    do j=1,ny
       do i=1,nx
          y1(i,j)= &
               x1(i,j-1)*A(i,j,1)+ &
               x1(i-1,j)*A(i,j,2)+ &
               x1(i  ,j)*A(i,j,3)+ &
               x1(i+1,j)*A(i,j,4)+ &
               x1(i,j+1)*A(i,j,5)
       enddo
    enddo
  end subroutine SpMV

  subroutine calc_res(x1,y1)
    real*8 x1(1-nb:nx+nb,1-nb:ny+nb)
    real*8 y1(1-nb:nx+nb,1-nb:ny+nb)
    integer i,j
    call SpMV(x1,y1)
    do j=1,ny
       do i=1,nx
          y1(i,j)=b(i,j)-y1(i,j)
       enddo
    enddo
  end subroutine calc_res

  subroutine calc_dot(x1,x2,tmp)
    implicit none
    real*8 x1(1-nb:nx+nb,1-nb:ny+nb)
    real*8 x2(1-nb:nx+nb,1-nb:ny+nb)
    real*8 tmp
    integer i,j
    tmp=0.0d0
    do j=1,ny
       do i=1,nx
         tmp=tmp+x1(i,j)*x2(i,j)
       enddo
    enddo

  end subroutine calc_dot

  subroutine calc_norm(x1,norm)
    implicit none
    real*8 x1(1-nb:nx+nb,1-nb:ny+nb)
    real*8 norm
    call calc_dot(x1,x1,norm)
    norm=sqrt(norm)
  end subroutine calc_norm

  subroutine Jacobi_method(x1,y1)
    implicit none
    real*8 x1(1-nb:nx+nb,1-nb:ny+nb)
    real*8 y1(1-nb:nx+nb,1-nb:ny+nb)

    integer i,j
    do j=1,ny
       do i=1,nx
          y1(i,j)= ( b(i,j) - (&
               x1(i,j-1)*A(i,j,1)+ &
               x1(i-1,j)*A(i,j,2)+ &
               x1(i+1,j)*A(i,j,4)+ &
               x1(i,j+1)*A(i,j,5) ) ) / A(i,j,3)
       enddo
    enddo
  end subroutine Jacobi_method

  subroutine view(x1)
    implicit none
    real*8 x1(1-nb:nx+nb,1-nb:ny+nb)
    integer i,j
    do j=1,ny
       do i=1,nx
          write(*,*)j,i,x1(i,j)
       enddo
    enddo
  end subroutine view

! y=a*x + b*y
  subroutine axpy(x1,y1,alpha,beta)
    implicit none
    real*8 x1(1-nb:nx+nb,1-nb:ny+nb)
    real*8 y1(1-nb:nx+nb,1-nb:ny+nb)
    real*8 alpha,beta
    integer i,j
    do j=1,ny
       do i=1,nx
          y1(i,j)=alpha*x1(i,j)+beta*y1(i,j)
       enddo
    enddo
  end subroutine axpy

  subroutine pJacobi(x1,y1)
    implicit none
    real*8 x1(1-nb:nx+nb,1-nb:ny+nb)
    real*8 y1(1-nb:nx+nb,1-nb:ny+nb)

    integer i,j
    do j=1,ny
       do i=1,nx
          y1(i,j)= x1(i,j)/ A(i,j,3)
       enddo
    enddo
  end subroutine PJacobi

  subroutine preconditioning(x1,y1)
    implicit none
    real*8 x1(1-nb:nx+nb,1-nb:ny+nb)
    real*8 y1(1-nb:nx+nb,1-nb:ny+nb)
!    call pJacobi(x1,y1)
    y1=x1
  end subroutine preconditioning
end module poisson2D

! f, nx, ny   -- 密
! f_,nx_,ny_  -- 粗
! 密 --> 粗
subroutine restriction_2D( &
     nb, &
     f,nx,ny,&
     f_,nx_,ny_)
  implicit none
  integer nb

  integer nx,ny
  real*8 f(1-nb:nx+nb,1-nb:ny+nb)

  integer nx_,ny_
  real*8 f_(1-nb:nx_+nb,1-nb:ny_+nb)

  integer i,j,i_,j_
  integer,parameter :: step=2

  real*8,parameter :: w0=4.0d0
  real*8,parameter :: w1=2.0d0
  real*8,parameter :: w2=1.0d0

  real*8 tmp
  real*8 weight
  real*8 ixl,ixr
  real*8 iyl,iyr

  do j=1,ny,step
     do i=1,nx,step
        tmp=f(i,j)*w0
        weight=w0
        ixl=1.0d0
        ixr=1.0d0

        iyl=1.0d0
        iyr=1.0d0
! ----
        if(i==1)then
           ixl=0.0d0
        endif
        if((i+step)>nx)then
           ixr=0.0d0
        endif
        weight=weight+ixl*w1+ixr*w1
        tmp=tmp+ixl*w1*f(i-1,j)+ixr*w1*f(i+1,j)

        if(j==1)then
           iyl=0.0d0
        endif
        if((i+step)>ny)then
           iyr=0.0d0
        endif
        weight=weight+iyl*w1+iyr*w1
        tmp=tmp+iyl*w1*f(i,j-1)+iyr*w1*f(i,j+1)
! -----
        weight=weight + &
             ixl*iyl*w2 + &
             ixl*iyr*w2 + &
             ixr*iyl*w2 + &
             ixr*iyr*w2

        tmp=tmp + &
             ixl*iyl*w2*f(i-1,j-1) + &
             ixl*iyr*w2*f(i-1,j+1) + &
             ixr*iyl*w2*f(i+1,j-1) + &
             ixr*iyr*w2*f(i+1,j+1)
        i_=CEILING(i/2.0d0)
        j_=CEILING(j/2.0d0)
        f_(i_,j_)=tmp/weight
     enddo
  enddo

end subroutine restriction_2D

! -----
! 境界条件は周期境界
! f, nx, ny   -- 密
! f_,nx_,ny_  -- 粗
! 粗 --> 密
subroutine prolongation_2D( &
     nb, &
     f,nx,ny,&
     f_,nx_,ny_)
  implicit none
  integer nb

  integer nx,ny
  real*8 f(1-nb:nx+nb,1-nb:ny+nb)

  integer nx_,ny_
  real*8 f_(1-nb:nx_+nb,1-nb:ny_+nb)

  integer i,j,i_,j_
  integer,parameter :: step=1

  real*8 tmp
  real*8 ixl,ixr
  real*8 iyl,iyr

  integer ixflag,iyflag

  do j=1,ny,step
     do i=1,nx,step
        i_=CEILING(i/2.0d0)
        j_=CEILING(j/2.0d0)
!        write(*,*)"x,y -- ",i_,j_," --> ",i,j
        ixl=1.0d0
        ixr=1.0d0

        iyl=1.0d0
        iyr=1.0d0
! ----
        if(i==1)then
           ixl=0.0d0
        endif
        if((i+step)>nx)then
           ixr=0.0d0
        endif

        if(j==1)then
           iyl=0.0d0
        endif
        if((j+step)>ny)then
           iyr=0.0d0
        endif
! -----
        if(mod(i,2)==1)then
           ixflag=0
        else
           ixflag=1
        endif

        if(mod(j,2)==1)then
           iyflag=0
        else
           iyflag=1
        endif

        if( (ixflag==0).and.(iyflag==0) )then
           tmp=f_(i_,j_)
!           write(*,*)"@ 1 @"
        endif
        if( (ixflag==0).and.(iyflag==1) )then
           tmp=(f_(i_  ,j_  )+ &
                f_(i_  ,j_+1)*iyr ) /2.0d0
!           write(*,*)"@ 2 @"
        endif
        if( (ixflag==1).and.(iyflag==0) )then
           tmp=(f_(i_  ,j_)+ &
                f_(i_+1,j_)*ixr ) /2.0d0
!           write(*,*)"@ 3 @"
        endif
        if( (ixflag==1).and.(iyflag==1) )then
           tmp=(f_(i_  ,j_  )+ &
                f_(i_+1,j_  )*ixr+ &
                f_(i_  ,j_+1)*iyr+ &
                f_(i_+1,j_+1)*ixr*iyr ) /4.0d0
!           write(*,*)"@ 4 @"
        endif
        f(i,j)=tmp
!        write(*,*)i,j,tmp
     enddo
  enddo

end subroutine prolongation_2D

subroutine call_Jacobi(itrmax,abstol,rtolmax)
  use poisson2D
  implicit none
  integer itrmax
  real*8 abstol,rtolmax
 
  integer itr
  real*8 x1(1-nb:nx+nb,1-nb:ny+nb)
  real*8 x2(1-nb:nx+nb,1-nb:ny+nb)
  real*8 r0(1-nb:nx+nb,1-nb:ny+nb)
  real*8 rnorm,bmax
  x1=0.0d0
  x2=0.0d0
  r0=0.0d0
  
  x1=x
  call calc_norm(b,bmax)
  do itr=1,itrmax
     call calc_res(x1,r0)
     call calc_norm(r0,rnorm)
!     write(*,*) itr,rnorm,rnorm/bmax
     if( (rnorm/bmax ) < rtolmax)then
        exit
     endif
     if( rnorm < abstol )then
        exit
     endif

     call Jacobi_method(x1,x2)
     x1=x2
  enddo
  x=x1
end subroutine call_Jacobi

subroutine call_Jacobi_(itrmax,abstol,rtolmax)
  use poisson2D
  implicit none
  integer itrmax
  real*8 abstol,rtolmax
 
  integer itr
  real*8 x1(1-nb:nx+nb,1-nb:ny+nb)
  real*8 x2(1-nb:nx+nb,1-nb:ny+nb)
  real*8 r0(1-nb:nx+nb,1-nb:ny+nb)
  real*8 rnorm,bmax
  x1=0.0d0
  x2=0.0d0
  r0=0.0d0
  
  x1=x
  call calc_norm(b,bmax)
  do itr=1,itrmax
     call calc_res(x1,r0)
     call calc_norm(r0,rnorm)
     write(*,*) itr,rnorm,rnorm/bmax
     if( (rnorm/bmax ) < rtolmax)then
        exit
     endif
     if( rnorm < abstol )then
        exit
     endif

     call Jacobi_method(x1,x2)
     x1=x2
  enddo
  x=x1
end subroutine call_Jacobi_

subroutine call_restriction( &
     r_1,nx_1,ny_1, &
     r_2,nx_2,ny_2)
  use poisson2D
  implicit none
  integer nx_1,ny_1,nx_2,ny_2
  real*8 r_1(1-nb:nx_1+nb,1-nb:ny_1+nb)
  real*8 r_2(1-nb:nx_2+nb,1-nb:ny_2+nb)
  r_2=0.0d0
! -----
! 制約 
  write(*,*)"restriction  ",nx_1,ny_1," --> ",nx_2,ny_2
  call restriction_2D( &
     nb, &
     r_1,nx_1,ny_1,&
     r_2,nx_2,ny_2)
  call free_mem()
  call initialize(nx_1,ny_1)
end subroutine call_restriction

subroutine call_MG(itrmax,abstol,rtolmax,mglevel,mGrid)
  use MGrid_2D
  use poisson2D
  implicit none
  integer itrmax
  real*8 abstol,rtolmax
  integer mglevel
  type(MG) mGrid(mglevel)
  integer itr
  integer i,j
  real*8 ,allocatable :: r_(:,:)

  real*8 rnorm,bmax

  integer nx_,ny_
  integer nx_h,ny_h
  integer inner_itrmax
  integer lv

  inner_itrmax=8
  call calc_norm(b,bmax)
  do lv=1,mglevel
     call free_mem()
     nx_=CEILING( nx0/(2.0d0**(lv-1) ) )
     ny_=CEILING( ny0/(2.0d0**(lv-1) ) )
     call initialize(nx_,ny_)
     mGrid(lv)%A=A
     mGrid(lv)%nx=nx
     mGrid(lv)%ny=ny
     write(*,*)"Level ",lv," -- ",nx,ny
  enddo

  mGrid(1)%x=0.0d0
  call free_mem()
  call initialize(nx0,ny0)
  do itr=1,itrmax
!  do itr=1,3
     call call_Jacobi(inner_itrmax,abstol,rtolmax)
!     call call_pcg(inner_itrmax,abstol,rtolmax)
     call calc_res(x,r)
     call calc_norm(r,rnorm)
     write(*,*) itr,rnorm,rnorm/bmax
     if( (rnorm/bmax ) < rtolmax)then
        exit
     endif
     if( rnorm < abstol )then
        exit
     endif

     mGrid(1)%x=x
     mGrid(1)%b=b
     mGrid(1)%r=r
! ------------
     do lv=2,mglevel
        call restriction_2D( &
             nb, &
             mGrid(lv-1)%r,mGrid(lv-1)%nx,mGrid(lv-1)%ny,&
             mGrid(lv)%b,mGrid(lv)%nx,mGrid(lv)%ny)
        call free_mem()
        call initialize(mGrid(lv)%nx,mGrid(lv)%ny)

        b=mGrid(lv)%b
        x=0.0d0
!        call call_Jacobi(inner_itrmax,abstol,rtolmax)
        call call_pcg(inner_itrmax,abstol,rtolmax)

        call calc_res(x,mGrid(lv)%r)
        mGrid(lv)%x=x
        mGrid(lv)%b=b
        mGrid(lv)%d=x
     enddo
!     return 

     ! ---
     ! 修正量の延長計算
     ! ---
     do lv=mglevel-1,1,-1
        call prolongation_2D( &
             nb, &
             mGrid(lv)%d,mGrid(lv)%nx,mGrid(lv)%ny,&
             mGrid(lv+1)%x,mGrid(lv+1)%nx,mGrid(lv+1)%ny)
        call free_mem()
        call initialize(mGrid(lv)%nx,mGrid(lv)%ny)
        call axpy(mGrid(lv)%d,mGrid(lv)%x,1.0d0,1.0d0)    
        x=mGrid(lv)%x
        b=mGrid(lv)%b
!        call call_Jacobi(inner_itrmax,abstol,rtolmax)
        call call_pcg(inner_itrmax,abstol,rtolmax)
        mGrid(lv)%x=x
     enddo

  enddo


  return

end subroutine call_MG

subroutine call_pcg(itrmax,abstol,rtolmax)
  use poisson2D
  implicit none
  integer itrmax
  real*8 abstol,rtolmax
  real*8 bmax,rnorm
  integer itr
  real*8 alpha,beta,gamma
  real*8 rho,rho_

  real*8  p(1-nb:nx+nb,1-nb:ny+nb)
  real*8  q(1-nb:nx+nb,1-nb:ny+nb)
  real*8 x1(1-nb:nx+nb,1-nb:ny+nb)
  real*8  s(1-nb:nx+nb,1-nb:ny+nb)
  real*8 r1(1-nb:nx+nb,1-nb:ny+nb)
  real*8 r2(1-nb:nx+nb,1-nb:ny+nb)
  real*8 const

  alpha=0.0d0
  beta=0.0d0
  gamma=0.0d0

  rho=0.0d0
  rho_=0.0d0

  p=0.0d0
  q=0.0d0
  x1=0.0d0
  s=0.0d0
  r1=0.0d0
  r2=0.0d0

  call calc_norm(b,bmax)
  x1=x
  alpha=0.0d0
  beta=0.0d0
  call calc_res(x1,r1)
  call  preconditioning(r1,s)

  call calc_dot(r1,s,rho)
  do itr=1,itrmax
! ----
     call calc_norm(r1,rnorm)
!!     call calc_res(x1,r2)
!!     call calc_norm(r2,rnorm)

     ! x = x_ + alpha*p_
     call axpy(p,x1,alpha,1.0d0)
     ! p = s_ + beta*p_
     call axpy(s,p,1.0d0,beta)

     write(*,*) itr,rnorm,rnorm/bmax
     if( (rnorm/bmax ) < rtolmax)then
        exit
     endif
     if( rnorm < abstol )then
        exit
     endif
! ----

     ! q = Ap
     call SpMV(p,q)
     ! gamma =(p,q)
     call calc_dot(p,q,gamma)
!     stop
     ! alpha= rho/gamma
!     write(*,*)"  rho,gamma  ",rho,gamma  
     alpha=rho/gamma

     ! r = r_ -alpha * q
     call axpy(q,r1,-alpha,1.0d0)

     ! solve s --> Ms=r 
     call  preconditioning(r1,s)
     ! rho=(s,r)
     rho_=rho
     call calc_dot(r1,s,rho)

     ! beta = rho/rho_
     beta = rho/rho_

  enddo
  x=x1

end subroutine call_pcg

subroutine call_mgcg(itrmax,abstol,rtolmax)
  use poisson2D
  implicit none
  integer itrmax
  real*8 abstol,rtolmax
  real*8 bmax,rnorm
  integer itr
  real*8 alpha,beta,gamma
  real*8 rho,rho_

  real*8  p(1-nb:nx+nb,1-nb:ny+nb)
  real*8  q(1-nb:nx+nb,1-nb:ny+nb)
  real*8 x1(1-nb:nx+nb,1-nb:ny+nb)
  real*8  s(1-nb:nx+nb,1-nb:ny+nb)
  real*8 r1(1-nb:nx+nb,1-nb:ny+nb)
  real*8 r2(1-nb:nx+nb,1-nb:ny+nb)
  real*8 const

  alpha=0.0d0
  beta=0.0d0
  gamma=0.0d0

  rho=0.0d0
  rho_=0.0d0

  p=0.0d0
  q=0.0d0
  x1=0.0d0
  s=0.0d0
  r1=0.0d0
  r2=0.0d0

  call calc_norm(b,bmax)
  x1=x
  alpha=0.0d0
  beta=0.0d0
  call calc_res(x1,r1)
  call  preconditioning(r1,s)

  call calc_dot(r1,s,rho)
  do itr=1,itrmax
! ----
     call calc_norm(r1,rnorm)
!!     call calc_res(x1,r2)
!!     call calc_norm(r2,rnorm)

     write(*,*) itr,rnorm,rnorm/bmax
     if( (rnorm/bmax ) < rtolmax)then
        exit
     endif
     if( rnorm < abstol )then
        exit
     endif
! ----
     ! x = x_ + alpha*p_
     call axpy(p,x1,alpha,1.0d0)
     ! p = s_ + beta*p_
     call axpy(s,p,1.0d0,beta)

     ! q = Ap
     call SpMV(p,q)

     ! gamma =(p,q)
     call calc_dot(p,q,gamma)

     ! alpha= rho/gamma
     alpha=rho/gamma

     ! r = r_ -alpha * q
     call axpy(q,r1,-alpha,1.0d0)

     ! solve s --> Ms=r 
     call  preconditioning(r1,s)
     ! rho=(s,r)
     rho_=rho
     call calc_dot(r1,s,rho)

     ! beta = rho/rho_
     beta = rho/rho_

  enddo
  x=x1

end subroutine call_mgcg


program main
  use MGrid_2D
  use poisson2D
  implicit none
  integer itr
  integer i,j
  integer itrmax
  integer nx_,ny_
  real*8 abstol,rtolmax
  integer mglevel
  real*8 ,allocatable :: x_(:,:)
  type(MG),allocatable :: mGrid(:)

! ----
  itrmax=200
  abstol=1e-20
  rtolmax=1e-20
! ----
  ! call initialize(nx0,ny0)
  ! call call_Jacobi_(itrmax,abstol,rtolmax)
  ! open(92,file='poisson2D_Jacobi.dat', status='replace')
  ! do j=1,ny
  !    do i=1,nx
  !       write(92,*)i,j,x(i,j)
  !    enddo
  !    write(92,*)
  ! enddo
  ! close(92)
  ! call free_mem()
! -----

  call initialize(nx0,ny0)
  call call_pcg(itrmax,abstol,rtolmax)
  open(92,file='poisson2D_pcg.dat', status='replace')
  do j=1,ny
     do i=1,nx
        write(92,*)i,j,x(i,j)
     enddo
     write(92,*)
  enddo
  close(92)
  call free_mem()
  stop
! -----
  write(*,*)
  write(*,*)" ------------- "
  call initialize(nx0,ny0)

! ----
  mglevel=8
  allocate(mGrid(mglevel))
  call make_MGrid(mGrid,mglevel,nb,nx0,ny0)
  do i=1,mglevel
     nx_=CEILING(nx0/(2.0d0**(i-1)) )
     ny_=CEILING(ny0/(2.0d0**(i-1)) )
     call free_mem()
     call initialize(nx_,ny_)
     mGrid(i)%A=A
     mGrid(i)%x=0.0d0
     mGrid(i)%r=0.0d0
     mGrid(i)%b=b
     mGrid(i)%d=0.0d0
  enddo
  call free_mem()
  call initialize(nx0,ny0)
! ----
  call call_MG(itrmax,abstol,rtolmax,mglevel,mGrid)

  open(92,file='poisson2D_MG.dat', status='replace')
  do j=1,ny
     do i=1,nx
        write(92,*)i,j,x(i,j)
     enddo
     write(92,*)
  enddo
  close(92)
  call free_mem()


  stop

  open(92,file='poisson2D_rest.dat', status='replace')
  nx_=CEILING(nx/2.0d0)
  ny_=CEILING(ny/2.0d0)
  allocate(x_(1-nb:nx_+nb,1-nb:ny_+nb))
  x_=0.0d0

  write(*,*)"restriction  ",nx,ny," --> ",nx_,ny_
  call restriction_2D( &
     nb, &
     x,nx,ny,&
     x_,nx_,ny_)
  do j=1,ny_
     do i=1,nx_
        write(92,*)2*i-1,2*j-1,x_(i,j)
     enddo
     write(92,*)
  enddo
  close(92)

  open(92,file='poisson2D_prol.dat', status='replace')
  r=x
  x=0.0d0
  nx_=CEILING(nx/2.0d0)
  ny_=CEILING(ny/2.0d0)

  call prolongation_2D( &
     nb, &
     x,nx,ny,&
     x_,nx_,ny_)
  write(*,*)"prolongation  ",nx_,ny_," --> ",nx,ny
  do j=1,ny
     do i=1,nx
        write(92,*)i,j,x(i,j)
!        write(92,*)i,j,( x(i,j)-r(i,j) )/ r(i,j)
     enddo
     write(92,*)
  enddo
  close(92)


  deallocate(x_)

end program main
