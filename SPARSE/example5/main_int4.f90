program main
  use problem_int8
  use krylov
  use krylov_kernel
  use matrix_kernel, only: &
       get_nnd_crs,init_diamat,crs2dia, &
       free_diaMat,free_crsMat, &
       free_cooMat,init_cooMat, &
       init_crsMat,coo2crs,view_crs
  use mpi
  implicit none
  type(CRS_Mat) :: crsMat 
  type(DIA_Mat) :: diaMat 
  type(KSP) :: method_CRS
  type(KSP) :: method_DIA

  integer ityp_eq
  integer(kind=4) gnx,gny,gnz
  integer itrmax
  integer(kind=4) nx,ny,nz
  integer(kind=4) nx0,ny0,nz0
  integer(kind=4) n,gn
  integer(kind=4) istart,iend

  integer(kind=4) nnz
  integer(kind=4) nnd
  integer npes,myrank
  integer ierr
  integer(kind=4) i
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

  real*8,allocatable  :: vecx(:)
  real*8,allocatable  :: vecb(:)

  real*8 stime,etime
  real*8 init_time_CRS,solve_time_CRS
  real*8 crs2dia_time,init_time_DIA,solve_time_DIA

  integer niter,restart

  namelist/input/ &
       ityp_eq, &
       gnx, &
       gny, &
       gnz, &
       itrmax, &
       rtolmax, &
       abstolmax, &
       ipreflag, &
       addL, &
       iflagAS, &
       ILU_method, &
       iovlflag, &
       solver, &
       restart, &
       GMRES_GSflag, &
       CBCG_kstep, &
       CBCG_Eigenflag, &
       power_method_itrmax, &
       CAARNOLDI_sstep, &
       CAARNOLDI_tstep, &
       CAARNOLDI_basis, &
       CAARNOLDI_QRflag, &
       CAGMRES_sstep, &
       CAGMRES_tstep, &
       CAGMRES_basis, &
       CAGMRES_QRflag

  open(33,file='input_namelist')
  read(33,nml=input)
  close(33)


  ! -----
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,npes,ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)

  nx     = gnx
  ny     = gny
  nz     = gnz/npes 
  nx0    = 1  
  ny0    = 1 
  nz0    = 1 + nz*myrank

  gn     = gnx*gny*gnz
  n      = nx * ny * nz
  istart = 1+n*myrank
  iend   = istart + n -1

  if(mod(gnz,npes).ne.0)then
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
  call make_poisson3D_CRS_get_nnz( &
     nx,ny,nz,   &
     nx0,ny0,nz0,&
     gnx,gny,gnz,&
     nnz, &
     istart)

  call init_crsMat(crsMat,n,nnz)
  call make_poisson3D_CRS( &
       nx,ny,nz,&
       nx0,ny0,nz0,&
       gnx,gny,gnz,&       
       crsMat%crsA,crsMat%crsCol,crsMat%crsRow_ptr,n,nnz, &
       istart)
  allocate(vecb(n))
  allocate(vecX(n))
  vecb    = -1.0d0
  vecX    = 0.0d0
  allocate(reshistory_CRS(itrmax))
  allocate(reshistory_DIA(itrmax))
  reshistory_CRS=0.0d0
  reshistory_DIA=0.0d0

  ! ------

  call PARCEL_KSP_Default_Setting(method_CRS)
  call set_KSP_rtolmax(method_CRS,rtolmax)
  call set_KSP_abstolmax(method_CRS,abstolmax)
  call set_KSP_solver(method_CRS,solver)
  call set_KSP_iovlflag(method_CRS,iovlflag)
  call set_KSP_iflagAS(method_CRS,iflagAS)
  call set_KSP_GMRES_GSflag(method_CRS,GMRES_GSflag)
  call set_KSP_GMRES_m(method_CRS,restart)
  call set_KSP_CBCG_kstep(method_CRS,CBCG_kstep)
  call set_KSP_CBCG_Eigenflag(method_CRS,CBCG_Eigenflag)
  call set_KSP_power_method_itrmax(method_CRS,power_method_itrmax)
  call set_KSP_CAARNOLDI_sstep(method_CRS,CAARNOLDI_sstep)
  call set_KSP_CAARNOLDI_tstep(method_CRS,CAARNOLDI_tstep)
  call set_KSP_CAARNOLDI_basis(method_CRS,CAARNOLDI_basis)
  call set_KSP_CAARNOLDI_QRflag(method_CRS,CAARNOLDI_QRflag)
  call set_KSP_CAGMRES_sstep(method_CRS,CAGMRES_sstep)
  call set_KSP_CAGMRES_tstep(method_CRS,CAGMRES_tstep)
  call set_KSP_CAGMRES_basis(method_CRS,CAGMRES_basis)
  call set_KSP_CAGMRES_QRflag(method_CRS,CAGMRES_QRflag)  

! ------
  call PARCEL_KSP_Default_Setting(method_DIA)
  call set_KSP_rtolmax(method_DIA,rtolmax)
  call set_KSP_abstolmax(method_DIA,abstolmax)
  call set_KSP_solver(method_DIA,solver)
  call set_KSP_iovlflag(method_DIA,iovlflag)
  call set_KSP_iflagAS(method_DIA,iflagAS)
  call set_KSP_GMRES_GSflag(method_DIA,GMRES_GSflag)
  call set_KSP_GMRES_m(method_DIA,restart)
  call set_KSP_CBCG_kstep(method_DIA,CBCG_kstep)
  call set_KSP_CBCG_Eigenflag(method_DIA,CBCG_Eigenflag)
  call set_KSP_power_method_itrmax(method_DIA,power_method_itrmax)
  call set_KSP_CAARNOLDI_sstep(method_DIA,CAARNOLDI_sstep)
  call set_KSP_CAARNOLDI_tstep(method_DIA,CAARNOLDI_tstep)
  call set_KSP_CAARNOLDI_basis(method_DIA,CAARNOLDI_basis)
  call set_KSP_CAARNOLDI_QRflag(method_DIA,CAARNOLDI_QRflag)
  call set_KSP_CAGMRES_sstep(method_DIA,CAGMRES_sstep)
  call set_KSP_CAGMRES_tstep(method_DIA,CAGMRES_tstep)
  call set_KSP_CAGMRES_basis(method_DIA,CAGMRES_basis)
  call set_KSP_CAGMRES_QRflag(method_DIA,CAGMRES_QRflag)  

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  stime=MPI_Wtime()
  call set_KSP_CRS(&
       MPI_COMM_WORLD, &
       vecx, &
       vecb, &
       n,gn,nnz,istart, &
       crsMat%crsA, &
       crsMat%crsRow_ptr, &
       crsMat%crsCol, &
       itrmax, &
       ipreflag, &
       ILU_method, &
       addL, &
       method_CRS)

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  etime=MPI_Wtime()  
  init_time_CRS=etime-stime
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  stime=MPI_Wtime()
  call START_KSP_SOLVER(method_CRS)

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  etime=MPI_Wtime()  
  solve_time_CRS=etime-stime

  call get_reshistory_KSP(niter,reshistory_CRS,method_CRS)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)


  ! ------
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  stime=MPI_Wtime()
  call get_nnd_crs(crsMat,n,istart,nnz,nnd)
  call init_diaMat(diaMat,n,nnd)
  call crs2dia(crsMat,diaMat,n,istart,nnd,nnz)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  etime=MPI_Wtime()  
  crs2dia_time=etime-stime

  write(*,*)myrank,"crs -- nnz = ",nnz, "dia -- nnd",nnd

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  stime=MPI_Wtime()
  call set_KSP_DIA(&
       MPI_COMM_WORLD, &
       vecx, &
       vecb, &
       n,gn,istart, &
       diaMat%diaA, &
       diaMat%offset,nnd, &
       itrmax, &
       ipreflag, &
       ILU_method, &
       addL, &
       method_DIA)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  etime=MPI_Wtime()  
  init_time_DIA=etime-stime

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  stime=MPI_Wtime()
  call START_KSP_SOLVER(method_DIA)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  etime=MPI_Wtime()  
  solve_time_DIA=etime-stime

  call get_reshistory_KSP(niter,reshistory_DIA,method_DIA)

  ! ------
  if(myrank.eq.0)then
     do i=1,niter
        write(*,*)"@ res @",&
             i,reshistory_CRS(i),reshistory_DIA(i)
     enddo
     write(*,*)"CRS Matrix"
     write(*,*)"   init  time ",init_time_CRS
     write(*,*)"   solve time ",solve_time_CRS

     write(*,*)"DIA Matrix"
     write(*,*)"   crs2dia  time ",crs2dia_time
     write(*,*)"   init     time ",init_time_DIA
     write(*,*)"   solve    time ",solve_time_DIA
  endif
  ! ------
 ! call view_crs(n,istart,crsMat)

  call free_crsMat(crsMat)
  call free_KSP(method_CRS)
  deallocate(reshistory_CRS)
  deallocate(vecx)
  deallocate(vecb)

  call MPI_Finalize(ierr)

end program main
