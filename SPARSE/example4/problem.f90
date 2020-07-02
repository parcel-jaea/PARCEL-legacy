subroutine make_poisson2D_matrix(nx,ny,nx0,ny0,gny,gnx,A,istart)
  implicit none
  integer nx,ny,nx0,ny0,gnx,gny
  integer istart
  real*8 A( (nx*ny) * (gnx*gny) ) ! nx*ny :: 行

  integer i,j,k
  real*8 dxi,dyi
  real*8 dxi2,dyi2
  real*8 dx,dy
  integer nw,ne,ns,nn,nc
  integer Row,Col
  real*8 cc,cw,ce,cs,cn
  integer n
  integer ix,iy
  n=nx*ny
  dx=1.0d0
  dy=1.0d0

  dxi=1.0d0/dx
  dyi=1.0d0/dy
  dxi2=dxi*dxi
  dyi2=dyi*dyi

  !$OMP PARALLEL DO PRIVATE(i) 
  do i=1,(nx*ny) * (gnx*gny)
     A(i)=0.0d0
  enddo
  !$OMP END PARALLEL DO 

  !$OMP PARALLEL DO & 
  !$OMP PRIVATE(i,ix,iy,Row,Col, &
  !$OMP nc,nw,ne,ns,nn, &
  !$OMP cc,cw,ce,cs,cn &
  !$OMP )
  do i=0,n-1
     ix=nx0 + mod(i,nx)
     iy=ny0 + i/nx 

     Row=i+1
     Col=Row + (istart-1)
     nc = 1+ (Row-1)*gnx*gny + (Col - 1)
     nw = nc - 1  
     ne = nc + 1  
     ns = nc - nx 
     nn = nc + nx 
        write(*,*)" nc,Row,Col --> ",nc,Row,Col
        write(*,*)" ix,iy--> ",ix,iy
#if 1
        cc=-2.0*(dxi2+dyi2)
        cw=dxi2
        ce=dxi2
        cs=dyi2
        cn=dyi2

        if(ix.eq.1)then
           nw=-1
        endif
        if(ix.eq.gnx)then
           ne=-1
        endif

        if(iy.eq.1)then
           ns=-1
        endif
        if(iy.eq.gny)then
           nn=-1
        endif

! -----
        A(nc)=cc
        if(nw > 0)then
           A(nw)=cw
        endif
        if(ne > 0)then
           A(ne)=ce
        endif
        if(ns > 0)then
           A(ns)=cs
        endif
        if(nn > 0)then
           A(nn)=cn
        endif
#endif
  enddo
  !$OMP END PARALLEL DO 

end subroutine make_poisson2D_matrix

subroutine make_poisson2D_matrix_CRS_get_nnz( &
     nx,ny,nx0,ny0,gny,gnx,&
     nnz, &
     istart)
  use matrix

  implicit none
  integer nx,ny,nx0,ny0,gnx,gny
  integer istart

  integer i,j,k
  real*8 dxi,dyi
  real*8 dxi2,dyi2
  real*8 dx,dy
  integer nw,ne,ns,nn,nc
  integer Row,Col
  real*8 cc,cw,ce,cs,cn
  integer n
  integer ix,iy

  integer nnz

  n=nx*ny
  dx=1.0d0
  dy=1.0d0

  dxi=1.0d0/dx
  dyi=1.0d0/dy
  dxi2=dxi*dxi
  dyi2=dyi*dyi

  nnz=0
  !$OMP PARALLEL DO & 
  !$OMP PRIVATE(i,ix,iy,Row,Col, &
  !$OMP nc,nw,ne,ns,nn, &
  !$OMP cc,cw,ce,cs,cn &
  !$OMP ) REDUCTION(+:nnz)
  ! do i=0,n-1
  !    ix=nx0 + mod(i,nx)
  !    iy=ny0 + i/nx 
  do iy=ny0,ny0+ny-1
  do ix=nx0,nx0+nx-1
     ! ix=nx0 + mod(i,nx)
     ! iy=ny0 + i/nx 
!     Row=i+1 + istart -1 
     Row=ix+(iy-1)*gnx

     nc = Row
     nw = nc - 1  
     ne = nc + 1  
     ns = nc - gnx 
     nn = nc + gnx 

     cc=-2.0*(dxi2+dyi2)
     cw=dxi2
     ce=dxi2
     cs=dyi2
     cn=dyi2

     if(ix.eq.1)then
        nw=-1
     endif
     if(ix.eq.gnx)then
        ne=-1
     endif
     
     if(iy.eq.1)then
        ns=-1
     endif
     if(iy.eq.gny)then
        nn=-1
     endif

! -----
     nnz=nnz+1
     !        A(nc)=cc
     if(nw > 0)then
        nnz=nnz+1
        !           A(nw)=cw
     endif
     if(ne > 0)then
        nnz=nnz+1
        !           A(ne)=ce
     endif
     if(ns > 0)then
        nnz=nnz+1
        !           A(ns)=cs
     endif
     if(nn > 0)then
        nnz=nnz+1
        !           A(nn)=cn
     endif
  enddo
  enddo
  !$OMP END PARALLEL DO 
end subroutine make_poisson2D_matrix_CRS_get_nnz

subroutine make_poisson2D_matrix_DIA_get_nnd(&
           diaMat,gn,n,nnz,crsCol,crsRow,istart,iend,nnd)
  use matrix
  use matrix_kernel, only: &
      init_diamat
  implicit none
  integer gn
  integer n
  integer nnz
  integer crsCol(nnz)
  integer crsRow(n+1)
  integer istart
  integer iend
  type(DIA_Mat) :: diaMat

  integer nnd
  
  integer i, j, icol, irow, itmp
  integer,allocatable::Flag(:)
  integer nnd_p

  integer imn,imx
  
  imn = 0
  imx = 0
  do i=1,n
    do j=crsRow(i),crsRow(i+1)-1
      irow = i+istart-1
      icol = crsCol(j)
      if(irow>=istart .and. irow<=iend) then
        itmp=icol-irow
        if (imn>itmp) imn = itmp
        if (imx<itmp) imx = itmp
      endif
    enddo
  enddo
  allocate(Flag(imn:imx))
  Flag = 0

  do i=1,n
    do j=crsRow(i),crsRow(i+1)-1
      irow = i+istart-1
      icol = crsCol(j)
      if(irow>=istart .and. irow<=iend) then
        itmp = icol-irow
        Flag(itmp) = 1
      endif
    enddo
  enddo
  nnd = 0

  do i=imn,imx
    if (Flag(i)==1) then
      nnd=nnd+1
    endif
  enddo
  call init_diaMat(diaMat,n,nnd)
  nnd_p=0
  do i=imn,imx
    if (flag(i)==1) then
      nnd_p=nnd_p+1
      diaMat%offset(nnd_p)=i
    endif
  enddo
  deallocate(flag)
  
end subroutine make_poisson2D_matrix_DIA_get_nnd




subroutine make_poisson2D_matrix_CRS( &
     nx,ny,nx0,ny0,gny,gnx,&
     crsA,crsCol,crsRow_ptr,n,nnz0, &
     istart)
  use matrix

  implicit none
  integer nx,ny,nx0,ny0,gnx,gny
  integer istart
  !  real*8 A( (nx*ny) * (gnx*gny) ) ! nx*ny :: 行

  integer i,j,k
  real*8 dxi,dyi
  real*8 dxi2,dyi2
  real*8 dx,dy
  integer nw,ne,ns,nn,nc
  integer Row,Col
  real*8 cc,cw,ce,cs,cn
  integer n
  integer ix,iy

  integer nnz0
  real*8  crsA(nnz0)
  integer crsCol(nnz0)
  integer crsRow_ptr(n+1)

  integer nnz

  n=nx*ny
  dx=1.0d0
  dy=1.0d0

  dxi=1.0d0/dx
  dyi=1.0d0/dy
  dxi2=dxi*dxi
  dyi2=dyi*dyi

  nnz=1
  !  do i=0,n-1
  i=0
  do iy=ny0,ny0+ny-1
     do ix=nx0,nx0+nx-1
        ! ix=nx0 + mod(i,nx)
        ! iy=ny0 + i/nx 
        ! Row=i+1 + istart -1 

        crsRow_ptr(i+1)=nnz
        Row=ix+(iy-1)*gnx

        nc = Row
        nw = nc - 1  
        ne = nc + 1  
        ns = nc - gnx 
        nn = nc + gnx 

        cc=-2.0*(dxi2+dyi2)
#if 1
        cw=dxi2
        ce=dxi2
        cs=dyi2
        cn=dyi2
#else
        call random_number(cw)
        call random_number(ce)
        call random_number(cs)
        call random_number(cn)

        cw=dxi2/cw
        ce=dxi2/ce
        cs=dyi2/cs
        cn=dyi2/cn
        cc=cw+ce+cs+cn
!        write(*,*)cw,ce,cc
#endif

        if(ix.eq.1)then
           nw=-1
        endif
        if(ix.eq.gnx)then
           ne=-1
        endif

        if(iy.eq.1)then
           ns=-1
        endif
        if(iy.eq.gny)then
           nn=-1
        endif

        ! -----
        if(ns > 0)then
           !           A(ns)
           crsA(nnz)=cs
           crsCol(nnz) = nc - gnx
           nnz=nnz+1
        endif
        if(nw > 0)then
           !           A(nw)
           crsA(nnz)=cw
           crsCol(nnz) = nc - 1
           nnz=nnz+1
        endif
        !        A(nc)
        crsA(nnz)=cc
        crsCol(nnz) = nc
        nnz=nnz+1

        if(ne > 0)then
           !           A(ne)
           crsA(nnz)=ce
           crsCol(nnz) = nc + 1
           nnz=nnz+1
        endif
        if(nn > 0)then
           !           A(nn)
           crsA(nnz)=cn
           crsCol(nnz) = nc + gnx
           nnz=nnz+1
        endif

        !      write(*,*)i+1-1,ns,cs
        !      write(*,*)i+1-1,nw,cw
        !      write(*,*)i+1-1,nc,cc
        !      write(*,*)i+1-1,ne,ce
        !      write(*,*)i+1-1,nn,cn
        i=i+1
     enddo
  enddo
  crsRow_ptr(n+1)=nnz

end subroutine make_poisson2D_matrix_CRS

subroutine make_poisson2D_matrix_asym_CRS( &
     nx,ny,nx0,ny0,gny,gnx,&
     crsA,crsCol,crsRow_ptr,n,nnz0, &
     istart)
  use matrix

  implicit none
  integer nx,ny,nx0,ny0,gnx,gny
  integer istart
  !  real*8 A( (nx*ny) * (gnx*gny) ) ! nx*ny :: 行

  integer i,j,k
  real*8 dxi,dyi
  real*8 dxi2,dyi2
  real*8 dx,dy
  integer nw,ne,ns,nn,nc
  integer Row,Col
  real*8 cc,cw,ce,cs,cn
  integer n
  integer ix,iy

  integer nnz0
  real*8  crsA(nnz0)
  integer crsCol(nnz0)
  integer crsRow_ptr(n+1)

  integer nnz

  n=nx*ny
  dx=1.0d0
  dy=1.0d0

  dxi=1.0d0/dx
  dyi=1.0d0/dy
  dxi2=dxi*dxi
  dyi2=dyi*dyi

  nnz=1
  !  do i=0,n-1
  i=0
  do iy=ny0,ny0+ny-1
     do ix=nx0,nx0+nx-1
        ! ix=nx0 + mod(i,nx)
        ! iy=ny0 + i/nx 
        ! Row=i+1 + istart -1 

        crsRow_ptr(i+1)=nnz
        Row=ix+(iy-1)*gnx

        nc = Row
        nw = nc - 1  
        ne = nc + 1  
        ns = nc - gnx 
        nn = nc + gnx 

        cc=-2.0*(dxi2+dyi2)
        cw=dxi2
#if 0
! symmetric
        ce=dxi2
#else
! non symmetric
         ce=dxi2*0.9d0
#endif
        cs=dyi2
#if 0
! symmetric
        cn=dyi2
#else
! non symmetric
        cn=dyi2*1.1d0
#endif
        if(ix.eq.1)then
           nw=-1
        endif
        if(ix.eq.gnx)then
           ne=-1
        endif

        if(iy.eq.1)then
           ns=-1
        endif
        if(iy.eq.gny)then
           nn=-1
        endif

        ! -----
        if(ns > 0)then
           !           A(ns)
           crsA(nnz)=cs
           crsCol(nnz) = nc - gnx
           nnz=nnz+1
        endif
        if(nw > 0)then
           !           A(nw)
           crsA(nnz)=cw
           crsCol(nnz) = nc - 1
           nnz=nnz+1
        endif
        !        A(nc)
        crsA(nnz)=cc
        crsCol(nnz) = nc
        nnz=nnz+1

        if(ne > 0)then
           !           A(ne)
           crsA(nnz)=ce
           crsCol(nnz) = nc + 1
           nnz=nnz+1
        endif
        if(nn > 0)then
           !           A(nn)
           crsA(nnz)=cn
           crsCol(nnz) = nc + gnx
           nnz=nnz+1
        endif

        !      write(*,*)i+1-1,ns,cs
        !      write(*,*)i+1-1,nw,cw
        !      write(*,*)i+1-1,nc,cc
        !      write(*,*)i+1-1,ne,ce
        !      write(*,*)i+1-1,nn,cn
        i=i+1
     enddo
  enddo
  crsRow_ptr(n+1)=nnz

end subroutine make_poisson2D_matrix_asym_CRS

subroutine make_poisson2D_matrix_DIA_from_CRS( &
     crsMat, &
     diaA,offset,nnd,n,nnz, &
     istart)
  use matrix

  implicit none
  type(CRS_Mat) :: crsMat
  type(DIA_Mat) :: diaMat
  real*8 crsA(nnz)

  integer istart
  !  real*8 A( (nx*ny) * (gnx*gny) ) ! nx*ny :: 行

  integer i,j,k
  integer n,itmp
  integer icol,irow
  integer off

  integer offset(nnd)
  
  integer nnz
  integer nnd
  real*8  diaA(n*nnd)
  
  do i=1,n
    do j=crsMat%crsRow_ptr(i),crsMat%crsRow_ptr(i+1)-1
      irow = i+istart-1
      icol = crsMat%crsCol(j)
      itmp = icol-irow
      do k = 1,nnd
         off =k 
         if (offset(k)==itmp) exit
      enddo
      diaA(n*(off-1)+i) = crsMat%crsA(j)
    enddo
  enddo
end subroutine make_poisson2D_matrix_DIA_from_CRS


subroutine make_random_matrix(n,gn,A)
  implicit none
  integer i,j
  integer gn,n
  real*8 A( gn * n ) 
  real*8 absA
  call random_number(A)

  do i=1,n
     do j=1,gn
        absA=abs(A(j+gn*(i-1)))
        if(absA>0.25d0)then
           A(j+gn*(i-1))=0.0d0
        endif
     enddo
  enddo
     
end subroutine make_random_matrix

subroutine make_poisson2D_9point_matrix_CRS_get_nnz( &
     nx,ny,nx0,ny0,gny,gnx,&
     nnz, &
     istart)
  use matrix

  implicit none
  integer nx,ny,nx0,ny0,gnx,gny
  integer istart

  integer i,j,k
  integer nw,ne,ns,nn,nc
  integer nw2,ne2,ns2,nn2
  integer Row,Col
  integer n
  integer ix,iy

  integer nnz
  n=nx*ny

  nnz=0
  !$OMP PARALLEL DO & 
  !$OMP PRIVATE(i,ix,iy,Row,Col, &
  !$OMP nc,nw,ne,ns,nn, &
  !$OMP nw2,ne2,ns2,nn2 &
  !$OMP ) REDUCTION(+:nnz)
  ! do i=0,n-1
  !    ix=nx0 + mod(i,nx)
  !    iy=ny0 + i/nx 
  do iy=ny0,ny0+ny-1
  do ix=nx0,nx0+nx-1
     ! ix=nx0 + mod(i,nx)
     ! iy=ny0 + i/nx 
!     Row=i+1 + istart -1 
     Row=ix+(iy-1)*gnx

     nc = Row
     nw = nc - 1  
     ne = nc + 1  
     ns = nc - gnx 
     nn = nc + gnx 

     nw2 = nc - 2  
     ne2 = nc + 2
     ns2 = nc - gnx*2
     nn2 = nc + gnx*2


        if(iy<=2)then
           ns2=-1
        endif
        if(iy.eq.1)then
           ns=-1
        endif

        if(ix<=2)then
           nw2=-1
        endif
        if(ix.eq.1)then
           nw=-1
        endif

        if(ix.eq.gnx)then
           ne=-1
        endif
        if(ix>=gnx-1)then
           ne2=-1
        endif
        if(iy.eq.gny)then
           nn=-1
        endif
        if(iy>=gny-1)then
           nn2=-1
        endif
! -----
        ! -----
        if(ns2 > 0)then
           nnz=nnz+1
        endif
        if(ns > 0)then
           nnz=nnz+1
        endif
        if(nw2 > 0)then
           nnz=nnz+1
        endif
        if(nw > 0)then
           nnz=nnz+1
        endif
        nnz=nnz+1

        if(ne > 0)then
           nnz=nnz+1
        endif

        if(ne2 > 0)then
           nnz=nnz+1
        endif

        if(nn > 0)then
           nnz=nnz+1
        endif

        if(nn2 > 0)then
           nnz=nnz+1
        endif
  enddo
  enddo
  !$OMP END PARALLEL DO 
end subroutine make_poisson2D_9point_matrix_CRS_get_nnz

subroutine make_poisson2D_9point_matrix_CRS( &
     nx,ny,nx0,ny0,gny,gnx,&
     crsA,crsCol,crsRow_ptr,n,nnz0, &
     istart)
  use matrix

  implicit none
  integer nx,ny,nx0,ny0,gnx,gny
  integer istart
  !  real*8 A( (nx*ny) * (gnx*gny) ) ! nx*ny :: 行

  integer i,j,k
  real*8 dxi,dyi
  real*8 dxi2,dyi2
  real*8 dx,dy
  integer nw,ne,ns,nn,nc
  integer nw2,ne2,ns2,nn2
  integer Row,Col
  real*8 cc,cw,ce,cs,cn
  real*8 cw2,ce2,cs2,cn2
  integer n
  integer ix,iy

  integer nnz0
  real*8  crsA(nnz0)
  integer crsCol(nnz0)
  integer crsRow_ptr(n+1)

  integer nnz

  n=nx*ny
  dx=1.0d0
  dy=1.0d0

  dxi=1.0d0/dx
  dyi=1.0d0/dy
  dxi2=dxi*dxi
  dyi2=dyi*dyi

  nnz=1
  !  do i=0,n-1
  i=0
  do iy=ny0,ny0+ny-1
     do ix=nx0,nx0+nx-1
        ! ix=nx0 + mod(i,nx)
        ! iy=ny0 + i/nx 
        ! Row=i+1 + istart -1 

        crsRow_ptr(i+1)=nnz
        Row=ix+(iy-1)*gnx

        nc = Row
        nw = nc - 1  
        ne = nc + 1  
        ns = nc - gnx 
        nn = nc + gnx 

        nw2 = nc - 2  
        ne2 = nc + 2
        ns2 = nc - gnx*2
        nn2 = nc + gnx*2

        cc=-5.0*(dxi2)
        cw=dxi2*4.0d0/3.0d0
        ce=dxi2*4.0d0/3.0d0
        cs=dyi2*4.0d0/3.0d0
        cn=dyi2*4.0d0/3.0d0

        cw2=-dxi2/12.0d0
        ce2=-dxi2/12.0d0
        cs2=-dyi2/12.0d0
        cn2=-dyi2/12.0d0

        if(iy<=2)then
           ns2=-1
        endif
        if(iy.eq.1)then
           ns=-1
        endif

        if(ix<=2)then
           nw2=-1
        endif
        if(ix.eq.1)then
           nw=-1
        endif

        if(ix.eq.gnx)then
           ne=-1
        endif
        if(ix>=gnx-1)then
           ne2=-1
        endif

        if(iy.eq.gny)then
           nn=-1
        endif
        if(iy>=gny-1)then
           nn2=-1
        endif

        ! -----
        if(ns2 > 0)then
           !           A(ns2)
           crsA(nnz)=cs2
           crsCol(nnz) = nc - gnx*2
           nnz=nnz+1
        endif
        if(ns > 0)then
           !           A(ns)
           crsA(nnz)=cs
           crsCol(nnz) = nc - gnx
           nnz=nnz+1
        endif
        if(nw2 > 0)then
           !           A(nw2)
           crsA(nnz)=cw2
           crsCol(nnz) = nc - 1*2
           nnz=nnz+1
        endif
        if(nw > 0)then
           !           A(nw)
           crsA(nnz)=cw
           crsCol(nnz) = nc - 1
           nnz=nnz+1
        endif
        !        A(nc)
        crsA(nnz)=cc
        crsCol(nnz) = nc
        nnz=nnz+1

        if(ne > 0)then
           !           A(ne)
           crsA(nnz)=ce
           crsCol(nnz) = nc + 1
           nnz=nnz+1
        endif

        if(ne2 > 0)then
           !           A(ne)
           crsA(nnz)=ce2
           crsCol(nnz) = nc + 1*2
           nnz=nnz+1
        endif

        if(nn > 0)then
           !           A(nn)
           crsA(nnz)=cn
           crsCol(nnz) = nc + gnx
           nnz=nnz+1
        endif

        if(nn2 > 0)then
           !           A(nn)
           crsA(nnz)=cn2
           crsCol(nnz) = nc + gnx*2
           nnz=nnz+1
        endif

        !      write(*,*)i+1-1,ns,cs
        !      write(*,*)i+1-1,nw,cw
        !      write(*,*)i+1-1,nc,cc
        !      write(*,*)i+1-1,ne,ce
        !      write(*,*)i+1-1,nn,cn
        i=i+1
     enddo
  enddo
  crsRow_ptr(n+1)=nnz

end subroutine make_poisson2D_9point_matrix_CRS

subroutine make_poisson2D_matrix_DMDA( &
     nx,ny,&
     nx0,ny0,&
     gny,gnx,&
     dmdaA &
     )
  use matrix

  implicit none
  integer nx,ny
  integer nx0,ny0,gnx,gny
  integer mx,my,mxy

  integer i,j,k
  real*8 dxi,dyi
  real*8 dxi2,dyi2
  real*8 dx,dy
  integer nw,ne,ns,nn,nc
  real*8 cc,cw,ce,cs,cn
  integer n
  integer ix,iy
  integer ixg,iyg
  real*8 dmdaA((nx+2)*(ny+2)*5)

  n=nx*ny
  dx=1.0d0
  dy=1.0d0

  dxi=1.0d0/dx
  dyi=1.0d0/dy
  dxi2=dxi*dxi
  dyi2=dyi*dyi

  dmdaA=0.0d0
  mx=nx+2
  my=ny+2
  mxy=mx*my

  do iyg=ny0,ny0+ny-1
     do ixg=nx0,nx0+nx-1
        ix=1+ixg-nx0
        iy=1+iyg-ny0
        i=(ix-1)+(iy-1)*mx+1+mx

        cc=-2.0*(dxi2+dyi2)
        cw=dxi2
#if 1
! symmetric
        ce=dxi2
#else
! non symmetric
        ce=dxi2*0.9d0
#endif
        cs=dyi2
#if 1
! symmetric
        cn=dyi2
#else
! non symmetric
        cn=dyi2*1.1d0
#endif
        if(ixg.eq.1)then
           nw=-1
        endif
        if(ixg.eq.gnx)then
           ne=-1
        endif

        if(iyg.eq.1)then
           ns=-1
        endif
        if(iyg.eq.gny)then
           nn=-1
        endif

        ! -----
        if(ns > 0)then
           dmdaA(i-mx)=cs
        endif
        if(nw > 0)then
           dmdaA(i-1)=cw
        endif
        dmdaA(i)=cc

        if(ne > 0)then
           dmdaA(i+1)=ce
        endif
        if(nn > 0)then
           dmdaA(i+mx)=cn
        endif

     enddo
  enddo

end subroutine make_poisson2D_matrix_DMDA
