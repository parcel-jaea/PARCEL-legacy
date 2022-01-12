
#define PARCEL_DIA_Solver 1

program main
  use PARCEL_integer
  use problem
  use krylov
  use krylov_kernel,only : free_ksp
  ! use PARCEL_wrapper,only : &
  !      parcel_dcg,parcel_dbicgstab,parcel_dgmres,parcel_dcagmres,parcel_dcbcg, &
  !      parcel_dcg_dia,parcel_dbicgstab_dia,parcel_dgmres_dia,parcel_dcagmres_dia,parcel_dcbcg_dia,&
  !      parcel_dpower,parcel_dpower_dia,&
  !      parcel_dcaarnoldi,parcel_dcaarnoldi_dia
  use matrix_kernel, only: &
       get_nnd_crs,init_diamat,crs2dia,&
       get_nnd_crs_without_sort,crs2dia_without_sort, &
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
  integer(kind=intType) gnx,gny,gnz
  integer itrmax
  integer MAX_NITER
  integer(kind=intType) nx,ny,nz
  integer(kind=intType) nx0,ny0,nz0
  integer(kind=intType) n,gn
  integer(kind=intType) istart,iend

  integer(kind=intType) nnz
  integer(kind=intType) nnd
  integer npes,myrank
  integer ierr
  integer(kind=intType) i
  integer ipreflag
  integer addL
  integer iflagAS
  real*8  rtolmax
  real*8  abstolmax
  integer solver,ityp_solver
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

  real*8,allocatable  :: vecx0(:)
  real*8,allocatable  :: vecx(:)
  real*8,allocatable  :: vecb(:)

  real*8 stime,etime
  real*8 init_time_CRS,solve_time_CRS
  real*8 crs2dia_time,init_time_DIA,solve_time_DIA

  integer niter,restart
  integer myranky,myrankz
  integer npey,npez
  integer precon_thblock
  integer independ_nvec
  integer nBlock

  integer ncount
  real*8,allocatable  :: EigenValueReal(:)
  real*8,allocatable  :: EigenValueImag(:)
  real*8,allocatable  :: EigenVector(:)
  real*8,allocatable  :: EigenErrorNorm(:)
  integer EmaxID,EminID

  namelist/input/ &
       ityp_eq, &
       gnx, &
       gny, &
       gnz, &
       MAX_NITER, &
       rtolmax, &
       abstolmax, &
       ipreflag, &
       independ_nvec, &
       addL, &
       iflagAS, &
       ILU_method, &
       iovlflag, &
       ityp_solver, &
       precon_thblock, &
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
  nBlock = 2000
  
  stime = 0.0d0
  etime = 0.0d0
  init_time_CRS  = 0.0d0
  solve_time_CRS = 0.0d0
  crs2dia_time   = 0.0d0
  init_time_DIA  = 0.0d0
  solve_time_DIA = 0.0d0

  itrmax = MAX_NITER
  solver = ityp_solver
  ! -----
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,npes,ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)

  if(gnz<npes)then
     npez = gnz
     npey = npes/npez
     if(( (npez*npey).ne.(npes) ) .or. ( mod(gny,npey).ne.0) )then
        if(myrank.eq.0)then
           write(*,*)"Error :: not correct process"
           write(*,*)"   process = ",npes
           write(*,*)"   gny     = ",gny
           write(*,*)"   gnz     = ",gnz
           write(*,*)"   npey    = ",npey
           write(*,*)"   npez    = ",npez
        endif
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call MPI_Finalize(ierr)
        stop
     endif
  else
     npey = 1 
     npez = npes
     if(( (npez*npey).ne.(npes) ) .or. ( mod(gnz,npez).ne.0) )then
        if(myrank.eq.0)then
           write(*,*)"Error :: not correct process"
           write(*,*)"   process = ",npes
           write(*,*)"   gny     = ",gny
           write(*,*)"   gnz     = ",gnz
           write(*,*)"   npey    = ",npey
           write(*,*)"   npez    = ",npez
        endif
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call MPI_Finalize(ierr)
        stop
     endif
  endif

  myranky = mod(myrank,npey)
  myrankz = myrank/npey
  nx     = gnx
  ny     = gny/npey
  nz     = gnz/npez
  nx0    = 1  
  ny0    = 1 + ny*myranky
  nz0    = 1 + nz*myrankz

  gn     = gnx*gny*gnz
  n      = nx * ny * nz
  istart = 1+n*myrank
  iend   = istart + n -1

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
  allocate(vecx(n))
  allocate(vecx0(n))
  vecb  = -1.0d0
  vecx  = 0.0d0
  vecx0 = 0.0d0

  call make_poisson3D_right_vector( &
       nx,ny,nz,   &
       nx0,ny0,nz0,&
       gnx,gny,gnz,&
       vecx0,vecb,&
       istart)

  allocate(reshistory_CRS(itrmax))
  allocate(reshistory_DIA(itrmax))
  reshistory_CRS=0.0d0
  reshistory_DIA=0.0d0

  ! ------

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  stime=MPI_Wtime()
  niter = MAX_NITER
  select case( solver )
  case(1)
     call  parcel_dcg(&
          MPI_COMM_WORLD, &
          vecx, &
          vecb, &
          n,gn,nnz,istart, &
          crsMat%crsA , &
          crsMat%crsRow_ptr,&
          crsMat%crsCol,&
          ipreflag,ILU_method,addL,iflagAS, &
          niter,rtolmax,&
          reshistory_CRS,&
          iovlflag,  &
          precon_thblock, &
          independ_nvec, &
          nBlock , &
          iret &
          )
  case(2)
     call  parcel_dbicgstab(&
          MPI_COMM_WORLD, &
          vecx, &
          vecb, &
          n,gn,nnz,istart, &
          crsMat%crsA , &
          crsMat%crsRow_ptr,&
          crsMat%crsCol,&
          ipreflag,ILU_method,addL,iflagAS, &
          niter,rtolmax,&
          reshistory_CRS,&
          iovlflag,  &
          precon_thblock, &
          independ_nvec, &
          nBlock , &
          iret &
          )
  case(3)
     call  parcel_dgmres(&
          MPI_COMM_WORLD, &
          vecx, &
          vecb, &
          n,gn,nnz,istart, &
          crsMat%crsA , &
          crsMat%crsRow_ptr,&
          crsMat%crsCol,&
          ipreflag,ILU_method,addL,iflagAS, &
          niter,rtolmax,&
          reshistory_CRS,&
          iovlflag,  &
          precon_thblock, &
          independ_nvec, &
          nBlock , &
          restart,&
          GMRES_GSflag, &
          iret &
          )
  case(4)
     call  parcel_dcagmres(&
          MPI_COMM_WORLD, &
          vecx, &
          vecb, &
          n,gn,nnz,istart, &
          crsMat%crsA , &
          crsMat%crsRow_ptr,&
          crsMat%crsCol,&
          ipreflag,ILU_method,addL,iflagAS, &
          niter,rtolmax,&
          reshistory_CRS,&
          iovlflag,  &
          precon_thblock, &
          independ_nvec, &
          nBlock , &
          CAGMRES_sstep, &
          CAGMRES_tstep, &
          CAGMRES_basis, &
          CAGMRES_QRflag, &
          iret &
          )
  case(5)
     call  parcel_dcbcg(&
          MPI_COMM_WORLD, &
          vecx, &
          vecb, &
          n,gn,nnz,istart, &
          crsMat%crsA , &
          crsMat%crsRow_ptr,&
          crsMat%crsCol,&
          ipreflag,ILU_method,addL,iflagAS, &
          niter,rtolmax,&
          reshistory_CRS,&
          iovlflag,  &
          precon_thblock, &
          independ_nvec,  &
          nBlock , &
          CBCG_kstep,CBCG_Eigenflag,power_method_itrmax, &
          CAARNOLDI_sstep,CAARNOLDI_tstep, &
          CAARNOLDI_basis,CAARNOLDI_QRflag, &
          iret &
          )
  case(51)

     allocate( EigenValueReal(CAARNOLDI_sstep*CAARNOLDI_tstep) )
     allocate( EigenValueImag(CAARNOLDI_sstep*CAARNOLDI_tstep) )
     allocate( EigenVector(CAARNOLDI_sstep*CAARNOLDI_tstep*n) )
     allocate( EigenErrorNorm(CAARNOLDI_sstep*CAARNOLDI_tstep) )
     call  parcel_dcaarnoldi(&
          MPI_COMM_WORLD, &
          vecx, &
          vecb, &
          n,gn,nnz,istart, &
          crsMat%crsA , &
          crsMat%crsRow_ptr,&
          crsMat%crsCol,&
          ipreflag,ILU_method,addL,iflagAS, &
          iovlflag,  &
          precon_thblock, &
          independ_nvec, &
          nBlock , &
          CAARNOLDI_sstep, &
          CAARNOLDI_tstep, &
          CAARNOLDI_basis, &
          CAARNOLDI_QRflag, &
          EigenValueReal, &
          EigenValueImag, &
          EigenVector, &
          EigenErrorNorm, &
          EmaxID,EminID,&
          iret &
          )
  end select

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  etime=MPI_Wtime()  
  solve_time_CRS=etime-stime

  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! ------

#if PARCEL_DIA_Solver
  vecx = 0.0d0
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  stime=MPI_Wtime()
  call get_nnd_crs_without_sort(&
       diaMat,n,nnz, &
       crsMat%crsCol,&
       crsMat%crsRow_ptr,&
       istart,&
       iend,&
       nnd &
       )
  call crs2dia_without_sort(crsMat,diaMat,nnd,n,nnz,istart)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  etime=MPI_Wtime()  
  crs2dia_time=etime-stime

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  stime=MPI_Wtime()

  select case( solver )
  case(1)
     call  parcel_dcg_dia(&
          MPI_COMM_WORLD, &
          vecx, &
          vecb, &
          n,gn,istart, &
          diaMat%diaA ,  &
          diaMat%offset,nnd, &
          ipreflag,ILU_method,addL,iflagAS, &
          MAX_NITER,rtolmax,&
          reshistory_DIA,&
          iovlflag,  &
          precon_thblock, &
          independ_nvec, &
          nBlock , &
          iret &
          )
  case(2)
     call  parcel_dbicgstab_dia(&
          MPI_COMM_WORLD, &
          vecx, &
          vecb, &
          n,gn,istart, &
          diaMat%diaA ,  &
          diaMat%offset,nnd, &
          ipreflag,ILU_method,addL,iflagAS, &
          MAX_NITER,rtolmax,&
          reshistory_DIA,&
          iovlflag,  &
          precon_thblock, &
          independ_nvec, &
          nBlock , &
          iret &
          )
  case(3)
     call  parcel_dgmres_dia(&
          MPI_COMM_WORLD, &
          vecx, &
          vecb, &
          n,gn,istart, &
          diaMat%diaA ,  &
          diaMat%offset,nnd, &
          ipreflag,ILU_method,addL,iflagAS, &
          MAX_NITER,rtolmax,&
          reshistory_DIA,&
          iovlflag,  &
          precon_thblock, &
          independ_nvec, &
          nBlock , &
          restart,GMRES_GSflag, &
          iret &
          )
  case(4)
     call  parcel_dcagmres_dia(&
          MPI_COMM_WORLD, &
          vecx, &
          vecb, &
          n,gn,istart, &
          diaMat%diaA ,  &
          diaMat%offset,nnd, &
          ipreflag,ILU_method,addL,iflagAS, &
          MAX_NITER,rtolmax,&
          reshistory_DIA,&
          iovlflag,  &
          precon_thblock, &
          independ_nvec, &
          nBlock , &
          CAGMRES_sstep, &
          CAGMRES_tstep, &
          CAGMRES_basis, &
          CAGMRES_QRflag, &
          iret &
          )
  case(5)
     call  parcel_dcbcg_dia(&
          MPI_COMM_WORLD, &
          vecx, &
          vecb, &
          n,gn,istart, &
          diaMat%diaA ,  &
          diaMat%offset,nnd, &
          ipreflag,ILU_method,addL,iflagAS, &
          MAX_NITER,rtolmax,&
          reshistory_DIA,&
          iovlflag,  &
          precon_thblock, &
          independ_nvec,  &
          nBlock , &
          CBCG_kstep,CBCG_Eigenflag,power_method_itrmax, &
          CAARNOLDI_sstep,CAARNOLDI_tstep, &
          CAARNOLDI_basis,CAARNOLDI_QRflag, &
          iret &
          )
  case(51)
     EigenValueReal = -1
     EigenValueImag = -1
     EigenVector = -1 
     EigenErrorNorm = -1
     call  parcel_dcaarnoldi_dia(&
          MPI_COMM_WORLD, &
          vecx, &
          vecb, &
          n,gn,istart, &
          diaMat%diaA ,  &
          diaMat%offset,nnd, &
          ipreflag,ILU_method,addL,iflagAS, &
          iovlflag,  &
          precon_thblock, &
          independ_nvec, &
          nBlock , &
          CAARNOLDI_sstep, &
          CAARNOLDI_tstep, &
          CAARNOLDI_basis, &
          CAARNOLDI_QRflag, &
          EigenValueReal, &
          EigenValueImag, &
          EigenVector, &
          EigenErrorNorm, &
          EmaxID,&
          EminID,&
          iret &
          )
  end select

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  etime=MPI_Wtime()  
  solve_time_DIA=etime-stime

#endif

  ! ------
  if(myrank.eq.0)then
     if(ityp_solver<10)then

        write(*,*)" residual norm "
#if PARCEL_DIA_Solver
        write(*,*)"     CRS    ,      DIA"
#endif
        do i=1,niter
#if PARCEL_DIA_Solver
           write(*,*)"PARCEL : @res@",&
                i,reshistory_CRS(i),reshistory_DIA(i)
#else
           write(*,*)"PARCEL : @res@",&
                i,reshistory_CRS(i)
#endif
        enddo
     else
        write(*,*)
        write(*,*) "Max Eigen Value        = ",EigenValueReal(EmaxID),EigenValueImag(EmaxID)
        write(*,*) "Error ||Ax-kx||/||kx|| = ",EigenErrorNorm(EmaxID)
        write(*,*) " -------------------- "
        write(*,*) "Min Eigen Value        = ", EigenValueReal(EminID),EigenValueImag(EminID)
        write(*,*) "Error ||Ax-kx||/||kx|| = ",EigenErrorNorm(EminID)
        write(*,*)
     endif
     
  endif

  if(ityp_solver<10)then
     call check_exact_solution(n,MPI_COMM_WORLD,vecx,vecx0)
  endif
  if(myrank.eq.0)then
     write(*,*)"CRS Matrix"
     write(*,*)"   solve time ",solve_time_CRS
#if PARCEL_DIA_Solver
     write(*,*)"DIA Matrix"
     write(*,*)"   crs2dia  time ",crs2dia_time
     write(*,*)"   solve    time ",solve_time_DIA
#endif
  endif
  ! ------

  call free_crsMat(crsMat)
  call free_KSP(method_CRS)

#if PARCEL_DIA_Solver
  call free_diaMat(diaMat)
  call free_KSP(method_DIA)
#endif
  if(ityp_solver<10)then
  else
     deallocate(EigenValueReal)
     deallocate(EigenValueImag)
     deallocate(EigenVector)
     deallocate(EigenErrorNorm)
  endif
  deallocate(vecx)
  deallocate(vecx0)
  deallocate(vecb)

  call MPI_Finalize(ierr)

end program main
