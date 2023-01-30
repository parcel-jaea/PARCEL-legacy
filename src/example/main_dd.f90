#define PARCEL_DIA_Solver 1

! 1プロセス限定
#define ToeplizMatrix 0 

#ifdef CUDA_SOLVER

program main
  implicit none
  write(*,*)"not support dd solver"
end program main

#else

program main
  use PARCEL_integer
  use problem
  use krylov 
  use krylov_kernel,only : free_ksp
  use PARCEL_wrapper
  ! use PARCEL_ddwrapper,only :&
  !      parcel_ddcg,parcel_ddbicgstab,parcel_ddgmres, &
  !      parcel_ddcg_dia,parcel_ddbicgstab_dia,parcel_ddgmres_dia
  ! use PARCEL_wrapper,only : parcel_dcg,parcel_dbicgstab,parcel_dgmres
  use common_kernel,only : divide_region
  use matrix_kernel,only : &
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
  integer i_
  integer ipreflag
  integer Neuman_itrmax
  integer addL
  integer iflagAS
  real*8  rtolmax
  real*8  abstolmax
  integer solver,ityp_solver
  integer iret
  real*8,allocatable  :: reshistory_CRS_dd(:)
  real*8,allocatable  :: reshistory_CRS_d(:)
  real*8,allocatable  :: reshistory_DIA_dd(:)
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
  real*8,allocatable  :: vecx0(:)
  real*8,allocatable  :: vecb(:)

  real*8,allocatable  :: crsA_lo(:)
  real*8,allocatable  :: diaA_lo(:)
  real*8,allocatable  :: vecx_lo(:)
  real*8,allocatable  :: vecb_lo(:)

  real*8 stime,etime
  real*8 init_time_CRS,solve_time_CRS
  real*8 crs2dia_time,init_time_DIA,solve_time_DIA

  integer niter,restart
  integer myranky,myrankz
  integer npey,npez
  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 gamma
  real*8 rtolmax0

  integer PRECISION_A
  integer PRECISION_b
  integer PRECISION_x
  integer PRECISION_PRECONDITION

  namelist/input/ &
       ityp_eq, &
       gnx, &
       gny, &
       gnz, &
       MAX_NITER, &
       rtolmax, &
       abstolmax, &
       ipreflag, &
       Neuman_itrmax, &
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
  namelist/input_dd/ &
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION

  open(33,file='input_namelist')
  read(33,nml=input)
  read(33,nml=input_dd)
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

#if ToeplizMatrix
  gn = gnx*gny*gnz
  call divide_region(myrank,npes,1,gn,istart,iend)
  n = iend - istart + 1 

  call make_ToeplizMatrix_CRS_get_nnz( &
       n,gn, &
       nnz, &
       istart)
  gamma = 2.0d0

  call init_crsMat(crsMat,n,nnz)
  allocate(crsA_lo(nnz))

  crsA_lo = 0.0d0

  call make_ToeplizMatrix_CRS( &
       n,gn,gamma, &
       crsMat%crsA,crsMat%crsCol,crsMat%crsRow_ptr,nnz, &
       istart)
#else
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
  call make_poisson3D_CRS_get_nnz( &
     nx,ny,nz,   &
     nx0,ny0,nz0,&
     gnx,gny,gnz,&
     nnz, &
     istart)

  call init_crsMat(crsMat,n,nnz)
  allocate(crsA_lo(nnz))

  crsA_lo = 0.0d0

  call make_poisson3D_CRS( &
       nx,ny,nz,&
       nx0,ny0,nz0,&
       gnx,gny,gnz,&       
       crsMat%crsA,crsMat%crsCol,crsMat%crsRow_ptr,n,nnz, &
       istart)
#endif
  allocate(vecb(n))
  allocate(vecx(n))

  allocate(vecb_lo(n))
  allocate(vecx_lo(n))

  vecb  = -1.0d0
  vecx  = 0.0d0

  vecb_lo = 0.0d0
  vecx_lo = 0.0d0
#if ToeplizMatrix
#else
  allocate(vecx0(n))
  vecx0 = 0.0d0
  call make_poisson3D_right_vector( &
       nx,ny,nz,   &
       nx0,ny0,nz0,&
       gnx,gny,gnz,&
       vecx0,vecb,&
       istart)

#endif
  !  call make_vecb(gnx,gny,gnz,vecb,vecx,istart,iend)

  allocate(reshistory_CRS_dd(itrmax))
  allocate(reshistory_CRS_d(itrmax))
  allocate(reshistory_DIA_dd(itrmax))
  reshistory_CRS_dd=0.0d0
  reshistory_CRS_d=0.0d0
  reshistory_DIA_dd=0.0d0

  ! ------

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  stime=MPI_Wtime()

  select case( solver )
  case(1)
     call  parcel_ddcg(&
          MPI_COMM_WORLD, &
          vecx,vecx_lo,&
          vecb,vecb_lo,  &
          n,gn,nnz,istart, &
          crsMat%crsA , crsA_lo, &
          crsMat%crsRow_ptr,&
          crsMat%crsCol,&
          ipreflag,ILU_method,addL,iflagAS, &
          MAX_NITER,rtolmax,&
          reshistory_CRS_dd,&
          iovlflag,  &
          precon_thblock, &
          independ_nvec, &
          nBlock,&
          PRECISION_A, &
          PRECISION_b, &
          PRECISION_x, &
          PRECISION_PRECONDITION, &
          iret &
          )
  case(2)
     call  parcel_ddbicgstab(&
          MPI_COMM_WORLD, &
          vecx,vecx_lo,&
          vecb,vecb_lo,  &
          n,gn,nnz,istart, &
          crsMat%crsA , crsA_lo, &
          crsMat%crsRow_ptr,&
          crsMat%crsCol,&
          ipreflag,ILU_method,addL,iflagAS, &
          MAX_NITER,rtolmax,&
          reshistory_CRS_dd,&
          iovlflag,  &
          precon_thblock, &
          independ_nvec, &
          nBlock,&
          PRECISION_A, &
          PRECISION_b, &
          PRECISION_x, &
          PRECISION_PRECONDITION, &
          iret &
          )
  case(3)
     call  parcel_ddgmres(&
          MPI_COMM_WORLD, &
          vecx,vecx_lo,&
          vecb,vecb_lo,  &
          n,gn,nnz,istart, &
          crsMat%crsA , crsA_lo, &
          crsMat%crsRow_ptr,&
          crsMat%crsCol,&
          ipreflag,ILU_method,addL,iflagAS, &
          MAX_NITER,rtolmax,&
          reshistory_CRS_dd,&
          iovlflag,  &
          precon_thblock, &
          independ_nvec, &
          nBlock,&
          restart,GMRES_GSflag, &
          PRECISION_A, &
          PRECISION_b, &
          PRECISION_x, &
          PRECISION_PRECONDITION, &
          iret &
          )
  end select
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  etime=MPI_Wtime()  
  solve_time_CRS=etime-stime
  niter = MAX_NITER
  vecx    = 0.0d0
  vecx_lo = 0.0d0
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
          MAX_NITER,rtolmax,&
          reshistory_CRS_d,&
          iovlflag,  &
          precon_thblock, &
          Neuman_itrmax, &
          independ_nvec, &
          nBlock,&
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
          MAX_NITER,rtolmax,&
          reshistory_CRS_d,&
          iovlflag,  &
          precon_thblock, &
          Neuman_itrmax, &
          independ_nvec, &
          nBlock,&
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
          MAX_NITER,rtolmax,&
          reshistory_CRS_d,&
          iovlflag,  &
          precon_thblock, &
          Neuman_itrmax, &
          independ_nvec, &
          nBlock,&
          restart,GMRES_GSflag, &
          iret &
          )
  end select

#if PARCEL_DIA_Solver
  MAX_NITER = itrmax 
  vecx    = 0.0d0
  vecx_lo = 0.0d0
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
  allocate(diaA_lo(nnd*n))
  diaA_lo = 0.0d0
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  etime=MPI_Wtime()  
  crs2dia_time=etime-stime

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  stime=MPI_Wtime()

  select case( solver )
  case(1)
     call  parcel_ddcg_dia(&
          MPI_COMM_WORLD, &
          vecx,vecx_lo,&
          vecb,vecb_lo,  &
          n,gn,istart, &
          diaMat%diaA , diaA_lo, &
          diaMat%offset,nnd, &
          ipreflag,ILU_method,addL,iflagAS, &
          MAX_NITER,rtolmax,&
          reshistory_DIA_dd,&
          iovlflag,  &
          precon_thblock, &
          independ_nvec, &
          nBlock,&
          PRECISION_A, &
          PRECISION_b, &
          PRECISION_x, &
          PRECISION_PRECONDITION, &
          iret &
          )
  case(2)
     call  parcel_ddbicgstab_dia(&
          MPI_COMM_WORLD, &
          vecx,vecx_lo,&
          vecb,vecb_lo,  &
          n,gn,istart, &
          diaMat%diaA , diaA_lo, &
          diaMat%offset,nnd, &
          ipreflag,ILU_method,addL,iflagAS, &
          MAX_NITER,rtolmax,&
          reshistory_DIA_dd,&
          iovlflag,  &
          precon_thblock, &
          independ_nvec, &
          nBlock,&
          PRECISION_A, &
          PRECISION_b, &
          PRECISION_x, &
          PRECISION_PRECONDITION, &
          iret &
          )
  case(3)
     call  parcel_ddgmres_dia(&
          MPI_COMM_WORLD, &
          vecx,vecx_lo,&
          vecb,vecb_lo,  &
          n,gn,istart, &
          diaMat%diaA , diaA_lo, &
          diaMat%offset,nnd, &
          ipreflag,ILU_method,addL,iflagAS, &
          MAX_NITER,rtolmax,&
          reshistory_DIA_dd,&
          iovlflag,  &
          precon_thblock, &
          independ_nvec, &
          nBlock,&
          restart,GMRES_GSflag, &
          PRECISION_A, &
          PRECISION_b, &
          PRECISION_x, &
          PRECISION_PRECONDITION, &
          iret &
          )
  end select
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  etime=MPI_Wtime()  
  solve_time_DIA=etime-stime
#endif




  ! ------
  if(myrank.eq.0)then
     write(*,*)" residual norm CRS matrix"
     write(*,*)"     double-double    ,      double"
     do i=1,niter
        i_ = i
        write(*,*)"PARCEL : @res@",&
             i_,reshistory_CRS_dd(i),reshistory_CRS_d(i)
     enddo

#if PARCEL_DIA_Solver
     write(*,*)
     write(*,*)" -------------------------------"
     write(*,*)
     write(*,*)" residual norm CRS matrix and DIA matrix"
     write(*,*)"     CRS ddmatrix    ,       DIA ddmatrix"
     do i=1,niter
        i_ = i
        write(*,*)"PARCEL : @res@",&
             i_,reshistory_CRS_dd(i),reshistory_DIA_dd(i)
     enddo
#endif

  endif
  
#if ToeplizMatrix
#else
  call check_exact_solution(n,MPI_COMM_WORLD,vecx,vecx0)
#endif
  if(myrank.eq.0)then
     write(*,*)"CRS Matrix"
     write(*,*)"   solve time ",solve_time_CRS
#if PARCEL_DIA_Solver
     write(*,*)"DIA Matrix"
     write(*,*)"   crs2dia  time ",crs2dia_time
     write(*,*)"   solve    time ",solve_time_DIA
#endif
  endif

  call free_crsMat(crsMat)
  call free_KSP(method_CRS)

#if PARCEL_DIA_Solver
  call free_diaMat(diaMat)
  call free_KSP(method_DIA)
#endif

  deallocate(vecx)
  deallocate(vecb)
  deallocate(reshistory_CRS_dd)
  deallocate(reshistory_DIA_dd)
  deallocate(reshistory_CRS_d)

  deallocate(crsA_lo)
  deallocate(vecb_lo)
  deallocate(vecx_lo)
#if ToeplizMatrix
#else
  deallocate(vecx0)
#endif
  
  call MPI_Finalize(ierr)

end program main
#endif
