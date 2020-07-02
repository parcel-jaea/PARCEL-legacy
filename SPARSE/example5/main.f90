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
  integer(kind=8) gnx,gny,gnz
  integer MAX_NITER
  integer(kind=8) nx,ny,nz
  integer(kind=8) nx0,ny0,nz0
  integer(kind=8) n,gn
  integer(kind=8) istart,iend

  integer(kind=8) nnz
  integer(kind=8) nnd
  integer npes,myrank
  integer ierr
  integer(kind=8) i
  integer ipreflag
  integer addL
  integer iflagAS
  real*8  rtolmax
  real*8  abstolmax
  integer ityp_solver
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
  integer niter,restart
  character(len=10) ::argc
  character(len=10) :: set_file_num
  character(len=70) :: file_name
  integer imode_flg,input_num


  !call getarg(1,argc)
  !set_file_num=trim(argc)


  namelist/input/ &
       ityp_eq, &
       gnx, &
       gny, &
       gnz, &
       MAX_NITER, &
       rtolmax, &
       abstolmax, &
       ipreflag, &
       addL, &
       iflagAS, &
       ILU_method, &
       iovlflag, &
       ityp_solver, &
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
  !read(set_file_num,*,iostat=imode_flg) input_num
  !if (imode_flg==0) then
  !  if ((input_num>=1).and.(input_num<=35)) then
  !    file_name='./templete/templete_input/input_namelist'//set_file_num
  !  else
  !    print *,"invalid input file number"
  !    stop
  !  endif
  !else
    file_name='./input_namelist'
  !endif

  open(33,file=file_name)
  read(33,nml=input)
  close(33)

  if ((ipreflag==3) .and.(ILU_method==0)) then
    write(*,*)"Additive schwartz doesn't support ILU precondition. Switch ILU to DILU condition 1."
    ILU_method=1
  endif
  if ((ipreflag==3).and.(addL==0)) then
    ipreflag=2
  endif

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
  call make_poisson3D_CRS_get_nnz_int8( &
     nx,ny,nz,   &
     nx0,ny0,nz0,&
     gnx,gny,gnz,&
     nnz, &
     istart)

  call init_crsMat(crsMat,n,nnz)
  call make_poisson3D_CRS_int8( &
       nx,ny,nz,&
       nx0,ny0,nz0,&
       gnx,gny,gnz,&       
       crsMat%crsA,crsMat%crsCol,crsMat%crsRow_ptr,n,nnz, &
       istart)
  allocate(vecb(n))
  allocate(vecX(n))
  vecb    = -1.0d0
  vecX    = 0.0d0
  call make_vecb_int8(gnx,gny,gnz,vecb,vecx,istart,iend)
  allocate(reshistory_CRS(MAX_NITER))
  allocate(reshistory_DIA(MAX_NITER))
  reshistory_CRS=0.0d0
  reshistory_DIA=0.0d0

  ! ------

  call PARCEL_KSP_Default_Setting(method_CRS)
  call set_KSP_rtolmax(method_CRS,rtolmax)
  call set_KSP_abstolmax(method_CRS,abstolmax)
  call set_KSP_solver(method_CRS,ityp_solver)
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
  call set_KSP_solver(method_DIA,ityp_solver)
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
       MAX_NITER, &
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

#if 0
!  call get_nnd_crs(crsMat,n,istart,nnz,nnd)
!  call init_diaMat(diaMat,n,nnd)
!  call crs2dia(crsMat,diaMat,n,istart,nnd,nnz)
  call make_poisson2D_matrix_DIA_get_nnd( &
       diaMat,gn,n,nnz,crsMat%crsCol,crsMat%crsRow_ptr,istart,iend,nnd)
  call make_poisson2D_matrix_DIA_from_CRS( &
       crsMat,diaMat%diaA,diaMat%offset,nnd,n,nnz,istart)


  write(*,*)myrank,"crs -- nnz = ",nnz, "dia -- nnd",nnd

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call set_KSP_DIA(&
       MPI_COMM_WORLD, &
       vecx, &
       vecb, &
       n,gn,istart, &
       diaMat%diaA, &
       diaMat%offset,nnd, &
       MAX_NITER, &
       ipreflag, &
       ILU_method, &
       addL, &
       method_DIA)
  call START_KSP_SOLVER(method_DIA)
  call get_reshistory_KSP(niter,reshistory_DIA,method_DIA)
#endif

  ! ------
  if(myrank.eq.0)then
     do i=1,niter
        write(*,*)"@ res @",&
!             i,reshistory_CRS(i),reshistory_DIA(i)
             i,reshistory_CRS(i)
     enddo
     write(*,*)"CRS Matrix"
     write(*,*)"   init  time ",init_time_CRS
     write(*,*)"   solve time ",solve_time_CRS
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
