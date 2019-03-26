program main
  use krylov
  use matrix_kernel, only: &
       get_nnd_crs,init_diamat,crs2dia, &
       free_diaMat,free_crsMat
  use mpi
  implicit none

  type(CRS_Mat) :: crsMat 
  type(DIA_Mat) :: diaMat 

  real*8,allocatable  :: vecX(:)
  real*8,allocatable  :: vecb(:)
  integer ityp_eq
  integer gnx,gny
  integer gnz
  integer itrmax
  integer nx,ny
  integer nx0,ny0
  integer n,gn
  integer istart,iend

  integer nnz
  integer nnd
  integer npes,myrank
  integer ierr
  integer i
  integer ipreflag
  integer addL
  integer iflagAS
  real*8  rtolmax
  real*8  abstolmax
  integer solver
  integer iret
  real*8,allocatable  :: reshistory_CRS(:)
  real*8,allocatable  :: reshistory_DIA(:)
  ! -----
  integer ILU_method
  integer iovlflag
  integer GMRES_m
  integer GMRES_GSflag
  integer CAGMRES_sstep
  integer CAGMRES_tstep
  integer CAGMRES_basis
  integer CAGMRES_QRflag
  integer CBCG_kstep
  integer CBCG_Eigenflag
  integer power_method_itrmax
  integer CAARNOLDI_sstep
  integer CAARNOLDI_tstep
  integer CAARNOLDI_basis
  integer CAARNOLDI_QRflag
  integer PRECISION_A
  integer PRECISION_b
  integer PRECISION_x
  integer PRECISION_PRECONDITION

  namelist/input/ &
       ityp_eq, &
       gnx,&
       gny, &
       itrmax, &
       abstolmax, &
       rtolmax, &
       ipreflag, &
       ILU_method, &
       iovlflag, &
       solver, &
       addL, &
       iflagAS, &
       CAGMRES_sstep, &
       CAGMRES_tstep, &
       CAGMRES_basis, &
       CAGMRES_QRflag
  open(33,file='input_namelist')
  read(33,nml=input)
  close(33)
  solver=4

  ! -----
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,npes,ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)

  nx     = gnx
  ny     = gny/npes 
  nx0    = 1  
  ny0    = 1 + ny*myrank
  gn     = gnx*gny
  n      = nx * ny
  istart = 1+n*myrank
  iend   = istart + n -1

  if(mod(gny,npes).ne.0)then
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     if(myrank.eq.0)then
        write(*,*)"Error :: not correct process"
        write(*,*)"   process = ",npes
        write(*,*)"   gny     = ",gny
     endif
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call MPI_Finalize(ierr)
     stop
  endif
  ! ------
  call make_poisson2D_matrix_CRS_get_nnz( &
       nx,ny,nx0,ny0,gny,gnx,&       
       nnz, &
       istart)
  allocate(crsMat%crsA(nnz))
  allocate(crsMat%crsCol(nnz))
  allocate(crsMat%crsRow_ptr(n+1))
  allocate(crsMat%RowL(1))
  allocate(crsMat%RowD(1))
  allocate(crsMat%RowU(1))
  if(ityp_eq.eq.1)then
     call make_poisson2D_matrix_CRS( &
          nx,ny,nx0,ny0,gny,gnx,&       
          crsMat%crsA,crsMat%crsCol,crsMat%crsRow_ptr,n,nnz, &
          istart)
  else
     call make_poisson2D_matrix_asym_CRS( &
          nx,ny,nx0,ny0,gny,gnx,&       
          crsMat%crsA,crsMat%crsCol,crsMat%crsRow_ptr,n,nnz, &
          istart)
  endif
  ! ----
  ! CRS --> DIA
  call get_nnd_crs(crsMat,n,istart,nnz,nnd)
  call init_diaMat(diaMat,n,nnd)
  call crs2dia(crsMat,diaMat,n,istart,nnd,nnz)
  allocate(vecb(n))
  allocate(vecX(n))
  vecb    = -1.0d0
  vecX    = 0.0d0
  allocate(reshistory_CRS(itrmax))
  allocate(reshistory_DIA(itrmax))
  reshistory_CRS=0.0d0
  reshistory_DIA=0.0d0
  ! ------

  vecX    = 0.0d0
  call parcel_dcagmres( &
       MPI_COMM_WORLD, &
       vecx, &
       vecb, &
       n,gn,nnz,istart, &
       crsMat%crsA, &
       crsMat%crsRow_ptr, &
       crsMat%crsCol, &
       ipreflag, &
       ILU_method, &
       addL, &
       iflagAS, &
       itrmax, &
       rtolmax, &
       reshistory_CRS, &
       iovlflag, &
       iret, &
       CAGMRES_sstep, &
       CAGMRES_tstep, &
       CAGMRES_basis, &
       CAGMRES_QRflag )

  vecx=0.0d0
  call parcel_dcagmres_dia( &
       MPI_COMM_WORLD, &
       vecx, &
       vecb, &
       n,gn,istart, &
       diaMat%diaA, &
       diaMat%offset, &
       nnd, &
       ipreflag, &
       ILU_method, &
       addL, &
       iflagAS, &
       itrmax, &
       rtolmax, &
       reshistory_DIA, &
       iovlflag, &
       iret, &
       CAGMRES_sstep, &
       CAGMRES_tstep, &
       CAGMRES_basis, &
       CAGMRES_QRflag )
  ! ------
  if(myrank.eq.0)then
     do i=1,itrmax
        write(*,*)"@ res @",&
             i,reshistory_CRS(i),reshistory_DIA(i)
     enddo
  endif
  deallocate(reshistory_CRS)
  deallocate(reshistory_DIA)
  deallocate(vecX)
  deallocate(vecb)

  call free_crsMat(crsMat)
  call free_diaMat(diaMat)
  call MPI_Finalize(ierr)
end program main
