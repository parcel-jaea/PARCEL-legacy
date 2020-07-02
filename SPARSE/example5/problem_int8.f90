module problem_int8

contains
subroutine make_poisson2D_CRS_get_nnz_int8( &
     nx,ny,nx0,ny0,gny,gnx,&
     nnz, &
     istart)
  use matrix

  implicit none
  integer(kind=8) nx,ny,nx0,ny0,gnx,gny
  integer(kind=8) istart

  integer(kind=8) i,j,k
  real*8 dxi,dyi
  real*8 dxi2,dyi2
  real*8 dx,dy
  integer(kind=8) nw,ne,ns,nn,nc
  integer(kind=8) Row,Col
  real*8 cc,cw,ce,cs,cn
  integer(kind=8) n
  integer(kind=8) ix,iy

  integer(kind=8) nnz

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
end subroutine make_poisson2D_CRS_get_nnz_int8

subroutine make_poisson2D_CRS_int8( &
     nx,ny,nx0,ny0,gny,gnx,&
     crsA,crsCol,crsRow_ptr,n,nnz0, &
     istart)
  use matrix

  implicit none
  integer(kind=8) nx,ny,nx0,ny0,gnx,gny
  integer(kind=8) istart
  !  real*8 A( (nx*ny) * (gnx*gny) ) ! nx*ny :: 行

  integer(kind=8) i,j,k
  real*8 dxi,dyi
  real*8 dxi2,dyi2
  real*8 dx,dy
  integer(kind=8) nw,ne,ns,nn,nc
  integer(kind=8) Row,Col
  real*8 cc,cw,ce,cs,cn
  integer(kind=8) n
  integer(kind=8) ix,iy

  integer(kind=8) nnz0
  real*8  crsA(nnz0)
  integer(kind=8) crsCol(nnz0)
  integer(kind=8) crsRow_ptr(n+1)

  integer(kind=8) nnz

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

end subroutine make_poisson2D_CRS_int8

subroutine make_poisson2D_matrix_DIA_get_nnd(&
           diaMat,gn,n,nnz,crsCol,crsRow,istart,iend,nnd)
  use matrix
  use matrix_kernel, only: &
      init_diamat
  implicit none
  integer(kind=8) gn
  integer(kind=8) n
  integer(kind=8) nnz
  integer(kind=8) crsCol(nnz)
  integer(kind=8) crsRow(n+1)
  integer(kind=8) istart
  integer(kind=8) iend
  type(DIA_Mat) :: diaMat

  integer(kind=8) nnd
  
  integer(kind=8) i, j, icol, irow, itmp
  integer(kind=8),allocatable::Flag(:)
  integer(kind=8) nnd_p

  integer(kind=8) imn,imx
  
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
  
end subroutine make_poisson2D_matrix_DIA_get_nnd

subroutine make_poisson2D_matrix_DIA_from_CRS( &
     crsMat, &
     diaA,offset,nnd,n,nnz, &
     istart)
  use matrix

  implicit none
  type(CRS_Mat) :: crsMat
  type(DIA_Mat) :: diaMat
  integer(kind=8) nnz
  real*8 crsA(nnz)

  integer(kind=8) istart
  !  real*8 A( (nx*ny) * (gnx*gny) ) ! nx*ny :: 行

  integer(kind=8) i,j,k
  integer(kind=8) n,itmp
  integer(kind=8) icol,irow
  integer(kind=8) off

  integer(kind=8) nnd
  integer(kind=8) offset(nnd)
  
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



subroutine make_poisson3D_CRS_int8( &
     nx,ny,nz,   &
     nx0,ny0,nz0,&
     gnx,gny,gnz,&
     crsA,crsCol,crsRow_ptr,n,nnz0, &
     istart)
  use matrix

  implicit none
  integer(kind=8) nx,ny,nz,nx0,ny0,nz0,gnx,gny,gnz
  integer(kind=8) istart
  !  real*8 A( (nx*ny) * (gnx*gny) ) ! nx*ny :: 行

  integer(kind=8) i,j,k
  real*8 dxi,dyi,dzi
  real*8 dxi2,dyi2,dzi2
  real*8 dx,dy,dz
  integer(kind=8) nc,nw,ne,ns,nn,nb,nt
  integer(kind=8) Row,Col
  real*8  cc,cw,ce,cs,cn,cb,ct
  integer(kind=8) n
  integer(kind=8) ix,iy,iz

  integer(kind=8) nnz0
  real*8  crsA(nnz0)
  integer(kind=8) crsCol(nnz0)
  integer(kind=8) crsRow_ptr(n+1)

  integer(kind=8) nnz

  n=nx*ny*nz
  dx=1.0d0
  dy=1.0d0
  dz=1.0d0

  dxi=1.0d0/dx
  dyi=1.0d0/dy
  dzi=1.0d0/dz

  dxi2=dxi*dxi
  dyi2=dyi*dyi
  dzi2=dzi*dzi

  nnz=1
  !  do i=0,n-1
  i=0
  do iz=nz0,nz0+nz-1
     do iy=ny0,ny0+ny-1
        do ix=nx0,nx0+nx-1
        ! ix=nx0 + mod(i,nx)
        ! iy=ny0 + i/nx 
        ! Row=i+1 + istart -1 

        crsRow_ptr(i+1)=nnz
        Row=ix + (iy-1)*gnx + (iz-1)*gnx*gny

        nc = Row
        nw = nc - 1  
        ne = nc + 1  
        ns = nc - gnx 
        nn = nc + gnx 
        nb = nc - gnx*gny
        nt = nc + gnx*gny

        cc=-2.0*(dxi2+dyi2+dzi2)
#if 1
        cw=dxi2
        ce=dxi2
        cs=dyi2
        cn=dyi2
        cb=dzi2
        ct=dzi2
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

        if(iz.eq.1)then
           nb=-1
        endif
        if(iz.eq.gnz)then
           nt=-1
        endif

        ! -----
        if(nb > 0)then
           !           A(nb)
           crsA(nnz)   = cb
           crsCol(nnz) = nb
           nnz=nnz+1
        endif

        if(ns > 0)then
           !           A(ns)
           crsA(nnz)=cs
           crsCol(nnz) = ns
           nnz=nnz+1
        endif
        if(nw > 0)then
           !           A(nw)
           crsA(nnz)=cw
           crsCol(nnz) =nw
           nnz=nnz+1
        endif
        !        A(nc)
        crsA(nnz)=cc
        crsCol(nnz) = nc
        nnz=nnz+1

        if(ne > 0)then
           !           A(ne)
           crsA(nnz)=ce
           crsCol(nnz) = ne
           nnz=nnz+1
        endif
        if(nn > 0)then
           !           A(nn)
           crsA(nnz)=cn
           crsCol(nnz) = nn
           nnz=nnz+1
        endif
        if(nt > 0)then
           !           A(nt)
           crsA(nnz)=ct
           crsCol(nnz) = nt
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
enddo
  crsRow_ptr(n+1)=nnz

end subroutine make_poisson3D_CRS_int8

subroutine make_poisson3D_CRS_get_nnz_int8( &
     nx,ny,nz,   &
     nx0,ny0,nz0,&
     gnx,gny,gnz,&
     nnz, &
     istart)
  use matrix

  implicit none
  integer(kind=8) nx,ny,nz,nx0,ny0,nz0,gnx,gny,gnz
  integer(kind=8) istart
  !  real*8 A( (nx*ny) * (gnx*gny) ) ! nx*ny :: 行

  integer(kind=8) i,j,k
  real*8 dxi,dyi,dzi
  real*8 dxi2,dyi2,dzi2
  real*8 dx,dy,dz
  integer(kind=8) nc,nw,ne,ns,nn,nb,nt
  integer(kind=8) Row,Col
  real*8  cc,cw,ce,cs,cn,cb,ct
  integer(kind=8) n
  integer(kind=8) ix,iy,iz

  integer(kind=8) nnz0
  integer(kind=8) nnz

  n=nx*ny*nz
  dx=1.0d0
  dy=1.0d0
  dz=1.0d0

  dxi=1.0d0/dx
  dyi=1.0d0/dy
  dzi=1.0d0/dz

  dxi2=dxi*dxi
  dyi2=dyi*dyi
  dzi2=dzi*dzi

  nnz=0
  !  do i=0,n-1
  i=0
  do iz=nz0,nz0+nz-1
     do iy=ny0,ny0+ny-1
        do ix=nx0,nx0+nx-1
        Row=ix + (iy-1)*gnx + (iz-1)*gnx*gny

        nc = Row
        nw = nc - 1  
        ne = nc + 1  
        ns = nc - gnx 
        nn = nc + gnx 
        nb = nc - gnx*gny
        nt = nc + gnx*gny

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

        if(iz.eq.1)then
           nb=-1
        endif
        if(iz.eq.gnz)then
           nt=-1
        endif

        ! -----
        if(nb > 0)then
           nnz=nnz+1
        endif

        if(ns > 0)then
           nnz=nnz+1
        endif
        if(nw > 0)then
           nnz=nnz+1
        endif
        nnz=nnz+1

        if(ne > 0)then
           nnz=nnz+1
        endif
        if(nn > 0)then
           nnz=nnz+1
        endif
        if(nt > 0)then
           nnz=nnz+1
        endif

     enddo
  enddo
enddo

end subroutine make_poisson3D_CRS_get_nnz_int8

subroutine make_vecb_int8(&
  mmx,mmy,mmz,vecB,x,istart,iend&
  )
integer(kind=8) ix,iy,iz
integer(kind=8) i
integer(kind=8),intent(in) :: istart,iend
integer(kind=8) mmx,mmy,mmz
real*8 h,z,ctmp1,ctmp2,ctmp3
real*8, intent(inout) :: x(istart:iend)
real*8, intent(out) :: vecB(istart:iend)


! --------------
! set (global) rhs vector
! --------------
  h = 1.0D0 / ((mmx+mmy+mmz)/3+1);   z = 2.0d0*h**6
  do iz=1, mmz
     do iy=1, mmy
        do ix=1, mmx
           i = mmz**2*(iz-1) +mmy*(iy-1) + ix
           ctmp1 = real(ix) **2;    ctmp2 = real(iy) **2;  ctmp3 = real(iz)**2
           if (i < istart .or. iend < i) cycle ! lm
           vecB(i) = z *( ctmp1*ctmp2 +ctmp2*ctmp3 +ctmp3*ctmp1 )
        end do
     end do
  end do

! --------------
! boundary conditions 
! --------------
  ix = mmx
  do iz=1, mmz
     do iy=1, mmy
        i = mmz**2 *(iz-1) +mmy*(iy-1) +ix
        if (i < istart .or. iend < i) cycle ! lm
        vecB(i) = vecB(i) -(h*iy)**2 *(h*iz)**2
     end do
  end do

  iy = mmy
  do iz=1, mmz
     do ix=1, mmx
        i = mmz**2 *(iz-1) +mmy*(iy-1) +ix
        if (i < istart .or. iend < i) cycle ! lm
        vecB(i) = vecB(i) - (h*ix)**2 * (h*iz)**2
     end do
  end do

  iz = mmz
  do iy=1, mmy
     do ix=1, mmx
        i = mmz**2 *(iz-1) +mmy*(iy-1) +ix
        if (i < istart .or. iend < i) cycle ! lm
        vecB(i) = vecB(i) - (h*ix)**2 * (h*iy)**2
     end do
  end do


end subroutine make_vecb_int8

end module problem_int8
