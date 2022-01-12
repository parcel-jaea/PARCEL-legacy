module problem
  use PARCEL_integer
contains
    
  subroutine make_ToeplizMatrix_CRS_get_nnz( &      
       n,gn, &
       nnz0, &
       istart)
    use matrix    
    implicit none
    integer(kind=intType) n,gn
    real*8 gamma
    integer(kind=intType) istart
    integer(kind=intType) nnz0
    real*8 crsA(nnz0)
    integer(kind=intType) crsCol(nnz0)
    integer(kind=intType) crsRow_ptr(n+1)

    integer(kind=intType) i
    integer(kind=intType) nnz
    
    nnz = 0
    do i = 1,n
       crsRow_ptr(i)=nnz
       if( ( (i+istart-1)-2 ) > 0 )then
          nnz = nnz + 1
       endif
       nnz = nnz + 1
       if( ( (i+istart-1) + 1 ) <= gn )then
          nnz = nnz + 1
       endif
    enddo
    nnz0 = nnz    
  end subroutine make_ToeplizMatrix_CRS_get_nnz

  subroutine make_ToeplizMatrix_CRS( &      
       n,gn,gamma, &
       crsA,crsCol,crsRow_ptr,nnz0, &
       istart)
    use matrix    
    implicit none
    integer(kind=intType) n,gn
    real*8 gamma
    integer(kind=intType) istart
    integer(kind=intType) nnz0
    real*8 crsA(nnz0)
    integer(kind=intType) crsCol(nnz0)
    integer(kind=intType) crsRow_ptr(n+1)

    integer(kind=intType) i
    integer(kind=intType) nnz
    
    nnz = 1
    do i = 1,n
       crsRow_ptr(i)=nnz
       if( ( (i+istart-1)-2 ) > 0 )then
          crsA(nnz)   = gamma
          crsCol(nnz) = i+istart-1 - 2          
          nnz = nnz + 1
       endif
       crsA(nnz)   = 2.0d0
       crsCol(nnz) = i+istart-1
       nnz = nnz + 1
       if( ( (i+istart-1) + 1 ) <= gn )then
          crsA(nnz)   = 1.0d0
          crsCol(nnz) = i+istart-1 + 1
          nnz = nnz + 1
       endif
       
    enddo
    crsRow_ptr(n+1)=nnz
    
  end subroutine make_ToeplizMatrix_CRS

  subroutine check_exact_solution(n,icomm,vecX,vecX0)
    use mpi
    implicit none
    integer(kind=intType) n
    integer icomm
    real*8 vecX(n)
    real*8 vecX0(n)
    integer(kind=intType) i
    real*8 err_max,exact_max
    real*8 tmp    
    integer myrank
    integer ierr 

    call MPI_Comm_rank(icomm,myrank,ierr)
!   call mpi_allreduce( X, Y, ndim, mpi_double_precision, mpi_sum, icomm, ierr )
!    call MPI_Reduce( zdmx, zdmx_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, icomm, ierr )
    tmp = 0.0d0
    do i = 1,n
       tmp = max(tmp, abs(vecX0(i) - vecX(i)) )
!       write(*,*)i,vecX0(i),vecX(i)
    enddo
    call MPI_Reduce(tmp,err_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,icomm,ierr)

    ! --------------------------

    tmp = 0.0d0
    do i = 1,n
       tmp = max(tmp, abs(vecX0(i) ) )
    enddo
    call MPI_Reduce(tmp,exact_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,icomm,ierr)

    if(myrank.eq.0)then
       write(*,*)"max(|numerical-exact|)                = ",err_max
       write(*,*)"max(|numerical-exact|) / max(|exact|) = ",err_max / exact_max
    endif

  end subroutine check_exact_solution

  subroutine make_poisson2D_CRS_get_nnz( &
       nx,ny,nx0,ny0,gny,gnx,&
       nnz, &
       istart)
    use matrix

    implicit none
    integer(kind=intType) nx,ny,nx0,ny0,gnx,gny
    integer(kind=intType) istart

    integer(kind=intType) i,j,k
    real*8 dxi,dyi
    real*8 dxi2,dyi2
    real*8 dx,dy
    integer(kind=intType) nw,ne,ns,nn,nc
    integer(kind=intType) Row,Col
    real*8 cc,cw,ce,cs,cn
    integer(kind=intType) n
    integer(kind=intType) ix,iy

    integer(kind=intType) nnz

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
  end subroutine make_poisson2D_CRS_get_nnz

  subroutine make_poisson2D_CRS( &
       nx,ny,nx0,ny0,gny,gnx,&
       crsA,crsCol,crsRow_ptr,n,nnz0, &
       istart)
    use matrix

    implicit none
    integer(kind=intType) nx,ny,nx0,ny0,gnx,gny
    integer(kind=intType) istart
    !  real*8 A( (nx*ny) * (gnx*gny) ) ! nx*ny :: 行

    integer(kind=intType) i,j,k
    real*8 dxi,dyi
    real*8 dxi2,dyi2
    real*8 dx,dy
    integer(kind=intType) nw,ne,ns,nn,nc
    integer(kind=intType) Row,Col
    real*8 cc,cw,ce,cs,cn
    integer(kind=intType) n
    integer(kind=intType) ix,iy

    integer(kind=intType) nnz0
    real*8  crsA(nnz0)
    integer(kind=intType) crsCol(nnz0)
    integer(kind=intType) crsRow_ptr(n+1)

    integer(kind=intType) nnz

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

  end subroutine make_poisson2D_CRS

  subroutine make_poisson2D_matrix_DIA_get_nnd(&
       diaMat,gn,n,nnz,crsCol,crsRow,istart,iend,nnd)
    use matrix
    use matrix_kernel, only: &
         init_diamat
    implicit none
    integer(kind=intType) gn
    integer(kind=intType) n
    integer(kind=intType) nnz
    integer(kind=intType) crsCol(nnz)
    integer(kind=intType) crsRow(n+1)
    integer(kind=intType) istart
    integer(kind=intType) iend
    type(DIA_Mat) :: diaMat

    integer(kind=intType) nnd

    integer(kind=intType) i, j, icol, irow, itmp
    integer(kind=intType),allocatable::Flag(:)
    integer(kind=intType) nnd_p

    integer(kind=intType) imn,imx

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
    integer(kind=intType) nnz
    real*8 crsA(nnz)

    integer(kind=intType) istart
    !  real*8 A( (nx*ny) * (gnx*gny) ) ! nx*ny :: 行

    integer(kind=intType) i,j,k
    integer(kind=intType) n,itmp
    integer(kind=intType) icol,irow
    integer(kind=intType) off

    integer(kind=intType) nnd
    integer(kind=intType) offset(nnd)

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

  subroutine make_poisson3D_CRS( &
       nx,ny,nz,   &
       nx0,ny0,nz0,&
       gnx,gny,gnz,&
       crsA,crsCol,crsRow_ptr,n,nnz0, &
       istart)
    use matrix

    implicit none
    integer(kind=intType) nx,ny,nz,nx0,ny0,nz0,gnx,gny,gnz
    integer(kind=intType) istart
    !  real*8 A( (nx*ny) * (gnx*gny) ) ! nx*ny :: 行

    integer(kind=intType) i,j,k
    real*8 dxi,dyi,dzi
    real*8 dxi2,dyi2,dzi2
    real*8 dx,dy,dz
    integer(kind=intType) nc,nw,ne,ns,nn,nb,nt
    integer(kind=intType) Row,Col
    real*8  cc,cw,ce,cs,cn,cb,ct
    integer(kind=intType) n
    integer(kind=intType) ix,iy,iz

    integer(kind=intType) nnz0
    real*8  crsA(nnz0)
    integer(kind=intType) crsCol(nnz0)
    integer(kind=intType) crsRow_ptr(n+1)

    integer(kind=intType) nnz

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
    
  end subroutine make_poisson3D_CRS

  subroutine make_poisson3D_right_vector( &
       nx,ny,nz,   &
       nx0,ny0,nz0,&
       gnx,gny,gnz,&
       vecX,vecB,&
       istart)

    implicit none
    integer(kind=intType) nx,ny,nz,nx0,ny0,nz0,gnx,gny,gnz
    real*8 vecX(nx*ny*nz)
    real*8 vecB(nx*ny*nz)
    integer(kind=intType) istart


    integer(kind=intType) i
    integer(kind=intType) ix,iy,iz
    real*8 z
    real*8 hx,hy,hz
    real*8 ctmpX,ctmpY,ctmpZ

    ! --------------
    ! set (global) rhs vector
    ! --------------
    hx = 1.0d0/(gnx+1)
    hy = 1.0d0/(gny+1)
    hz = 1.0d0/(gnz+1)

!    z = 2.0d0*h**6
    z = 2.0d0*(hx**2)*(hy**2)*(hz**2)
    do iz=1, nz
       do iy=1, ny
          do ix=1, nx
!             i = mm**2*(iz-1) +mm*(iy-1) + ix
             ctmpX = real(ix+nx0-1) **2
             ctmpY = real(iy+ny0-1) **2
             ctmpZ = real(iz+nz0-1) **2
             i = ix + (iy-1)*nx + (iz-1)*nx*ny
             vecB(i) = z *( ctmpX*ctmpY +ctmpY*ctmpZ +ctmpZ*ctmpX )
          end do
       end do
    end do

    ! --------------
    ! boundary conditions 
    ! --------------
    if( (nx0+nx-1).eq.gnx )then
       ix = nx
       do iz=1, nz
          do iy=1, ny
             i = ix + (iy-1)*nx + (iz-1)*nx*ny
             vecB(i) = vecB(i) - (hy*(iy+ny0-1))**2 * (hz*(iz+nz0-1))**2
          end do
       end do
    endif
    if( (ny0+ny-1).eq.gny )then
       iy = ny
       do iz=1, nz
          do ix=1, nx
             i = ix + (iy-1)*nx + (iz-1)*nx*ny
             vecB(i) = vecB(i) - (hx*(ix+nx0-1))**2 * (hz*(iz+nz0-1))**2
          end do
       end do
    endif

    if( (nz0+nz-1).eq.gnz )then
       iz = nz
       do iy=1, ny
          do ix=1, nx
             i = ix + (iy-1)*nx + (iz-1)*nx*ny
             vecB(i) = vecB(i) - (hx*(ix+nx0-1))**2 * (hy*(iy+ny0-1))**2
          end do
       end do
    endif

    ! --------------
    ! set exact solution
    ! --------------
    z = (hx**2)*(hy**2)*(hz**2) 
    do iz=1, nz
       do iy=1, ny
          do ix=1, nx
             i = ix + (iy-1)*nx + (iz-1)*nx*ny
             ctmpX = real(ix+nx0-1) **2
             ctmpY = real(iy+ny0-1) **2
             ctmpZ = real(iz+nz0-1) **2
             vecX(i) = z * ( ctmpX * ctmpY * ctmpZ )
!             write(*,*)ix,iy , vecX(i) 
!             write(*,*)i, vecB(i) 
          end do
       end do
!       write(*,*)" ------------ "
    end do

#if 0
    ! --------------
    ! set itital guess
    ! --------------
    x = 0.0d0
    if (1 >= istart .and. iend >= 1) then
       X(1) = 1.0
    end if
#endif

  end subroutine make_poisson3D_right_vector

  subroutine make_poisson3D_CRS_get_nnz( &
       nx,ny,nz,   &
       nx0,ny0,nz0,&
       gnx,gny,gnz,&
       nnz, &
       istart)
    use matrix

    implicit none
    integer(kind=intType) nx,ny,nz,nx0,ny0,nz0,gnx,gny,gnz
    integer(kind=intType) istart
    !  real*8 A( (nx*ny) * (gnx*gny) ) ! nx*ny :: 行

    integer(kind=intType) i,j,k
    real*8 dxi,dyi,dzi
    real*8 dxi2,dyi2,dzi2
    real*8 dx,dy,dz
    integer(kind=intType) nc,nw,ne,ns,nn,nb,nt
    integer(kind=intType) Row,Col
    real*8  cc,cw,ce,cs,cn,cb,ct
    integer(kind=intType) n
    integer(kind=intType) ix,iy,iz

    integer(kind=intType) nnz0
    integer(kind=intType) nnz

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

  end subroutine make_poisson3D_CRS_get_nnz

  subroutine make_vecb(&
       mmx,mmy,mmz,vecB,x,istart,iend&
       )
    integer(kind=intType) ix,iy,iz
    integer(kind=intType) i
    integer(kind=intType),intent(in) :: istart,iend
    integer(kind=intType) mmx,mmy,mmz
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


  end subroutine make_vecb

end module problem
