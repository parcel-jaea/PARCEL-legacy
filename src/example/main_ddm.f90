
program main
  use PARCEL_integer
  use problem
  use ddm_setup
  use krylov
  use mpi
  use ddm_wrapper
  implicit none

  integer ityp_eq
  integer(kind=intType) gnx,gny,gnz
  integer itrmax
  integer MAX_NITER
  integer(kind=intType) nx,ny,nz
  integer(kind=intType) nx0,ny0,nz0
  integer(kind=intType) n,gn

  integer(kind=intType) m
  integer(kind=intType) m_
  integer npes,myrank
  integer ierr
  integer(kind=intType) i
  integer ipreflag
  integer Neuman_itrmax
  integer addL
  integer iflagAS
  real*8  rtolmax
  real*8  abstolmax
  integer solver,ityp_solver
  integer iret
  real*8,allocatable  :: reshistory_DDM(:)
  ! -----
  integer ILU_method
  integer iovlflag
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
  real*8 init_time_DDM,solve_time_DDM

  integer niter,restart
  integer npey,npez
  integer precon_thblock
  integer independ_nvec

  integer npe_x,npe_y,npe_z
  integer div_direc_th

  integer rank_x, rank_y, rank_z
  integer,parameter :: nnd = 7
  integer(kind=intType),parameter :: nnd_ = 7
  integer,parameter :: margin = 1
  integer,parameter :: ndim=3 
  integer npe_dim(ndim) 
  integer gn_grids(ndim) 
  integer n_grids(ndim) 
  integer istart_grids(ndim) 
  integer,parameter :: num_neighbor_ranks_max = 3**ndim -1
  integer neighbor_ranks( num_neighbor_ranks_max )
  integer num_neighbor_ranks
  real*8, allocatable  :: diaA(:)
  integer, allocatable :: offset(:)
  integer, allocatable :: offset_dim(:,:)
  integer ncount
  real*8,allocatable  :: EigenValueReal(:)
  real*8,allocatable  :: EigenValueImag(:)
  real*8,allocatable  :: EigenVector(:)
  real*8,allocatable  :: EigenErrorNorm(:)
  integer nBlock
  integer EmaxID,EminID
  integer(kind=intType),parameter :: offset_zero=0
  integer(kind=intType) parcel_int

#ifdef CUDA_SOLVER
  type(devicemem) cu_vecx
  type(devicemem) cu_vecb
  type(devicemem) cu_diaA
  type(devicemem) cu_offset  
#endif

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

  namelist/input_ddm/ &
       npe_x,&
       npe_y,&
       npe_z,&
       div_direc_th

  open(33,file='input_namelist')
  read(33,nml=input)
  read(33,nml=input_ddm)
  close(33)
  nBlock = 2000

  stime = 0.0d0
  etime = 0.0d0

  itrmax = MAX_NITER
  solver = ityp_solver
  ! -----
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,npes,ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)

  call PARCEL_Init(myrank)

  if( npes.ne. ( npe_x* npe_y* npe_z) )then
     if(myrank.eq.0)then
        write(*,*)"Error :: not correct process"
        write(*,*)"   process = ",npes
        write(*,*)"   process = ",npes
        write(*,*)"   npe_x   = ",npe_x
        write(*,*)"   npe_y   = ",npe_y
        write(*,*)"   npe_z   = ",npe_z
     endif
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call MPI_Finalize(ierr)
     stop
  endif
  if(  mod(gnx,npe_x).ne.0 )then
     if(myrank.eq.0)then
        write(*,*)"Error :: not correct process"
        write(*,*)"         mod(gnx,npe_x).ne.0"
        write(*,*)"   gnx     = ",gnx
        write(*,*)"   npe_x   = ",npe_x
     endif
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call MPI_Finalize(ierr)
     stop
  endif

  if(  mod(gny,npe_y).ne.0 )then
     if(myrank.eq.0)then
        write(*,*)"Error :: not correct process"
        write(*,*)"         mod(gny,npe_y).ne.0"
        write(*,*)"   gny     = ",gny
        write(*,*)"   npe_y   = ",npe_y
     endif
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call MPI_Finalize(ierr)
     stop
  endif
  if(  mod(gnz,npe_z).ne.0 )then
     if(myrank.eq.0)then
        write(*,*)"Error :: not correct process"
        write(*,*)"         mod(gnz,npe_z).ne.0"
        write(*,*)"   gnz     = ",gnz
        write(*,*)"   npe_z   = ",npe_z
     endif
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call MPI_Finalize(ierr)
     stop
  endif

  npe_dim      = -1 
  gn_grids     = -1 
  n_grids      = -1 
  istart_grids = -1  
  gn_grids(1)  = gnx
  gn_grids(2)  = gny
  gn_grids(3)  = gnz

  call set_rank_xyz( myrank, npe_x, npe_y, npe_z, rank_x, rank_y, rank_z )
  call set_grids_decompostion( &
       npe_x, npe_y, npe_z, &
       rank_x, rank_y, rank_z,  &
       ndim, gn_grids, npe_dim, n_grids,  &
       istart_grids, margin, m )

  ! プロセス分割が2次元の場合は num_neighbor_ranks = 8
  ! プロセス分割が3次元の場合は num_neighbor_ranks = 26
  ! 問題の次元数ではなくプロセス分割数によりnum_neighbor_ranksが変化する
  ! 
  num_neighbor_ranks = -1
  call set_neighbor_ranks( &
       npe_x, npe_y, npe_z, rank_x, rank_y, rank_z, &
       num_neighbor_ranks, &
       neighbor_ranks, num_neighbor_ranks_max)


  allocate( diaA(nnd*m) )
  allocate( offset(nnd) )
  allocate( offset_dim(nnd,ndim) )
  allocate( reshistory_DDM(MAX_NITER+100) )
  allocate( vecx(m) )
  allocate( vecx0(m) )
  allocate( vecb(m) )

  m_ = m
  diaA = 0.0d0
  vecx  = 0.0d0
  vecx0 = 0.0d0
  vecb  = 0.0d0
  reshistory_DDM = 0
  offset = 0
  offset_dim = 0

!  nnd_parcelint = nnd
  call make_poisson_3D_DDM( &
       gnx,gny,gnz, &
       ndim, n_grids,  &
       istart_grids, &
       margin, nnd_, m, &
       offset_dim, offset, &
       diaA, vecb, vecx, vecx0,&
       rank_x, rank_y, rank_z, &
       npe_x, npe_y, npe_z )
#if 0
  vecb  = 1.0d0
  vecx0 = 0.0d0
  vecx  = 0.0d0
#endif

#if 0
  do i=1,m
     if( vecb(i) .ne. 0.0d0)then
        write(*,*)"vecb : ",i,vecb(i) 
     endif
  enddo
#endif
     
#ifdef CUDA_SOLVER
    call cudaMalloc_FP64(cu_vecx,m) 
    call cudaMalloc_FP64(cu_vecb,m)
    call cudaMalloc_FP64(cu_diaA,m*nnd) 
    parcel_int = nnd
    call cudaMalloc_INT(cu_offset,parcel_int) 

    call cudaMemcpy_h2d_fp64(vecx,cu_vecx,m,offset_zero,offset_zero)    
    call cudaMemcpy_h2d_fp64(vecb,cu_vecb,m,offset_zero,offset_zero)    
    call cudaMemcpy_h2d_fp64(diaA,cu_diaA,m*nnd,offset_zero,offset_zero)    

    parcel_int = nnd
    call cudaMemcpy_h2d_INT(offset,cu_offset, parcel_int,offset_zero,offset_zero)
#endif

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  stime=MPI_Wtime()
  etime=MPI_Wtime()  
  init_time_DDM=etime-stime

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  stime=MPI_Wtime()
  select case( solver )
  case(1)

     call parcel_dcg_ddm(&
          MPI_COMM_WORLD, &
#ifdef CUDA_SOLVER
          cu_vecx%mem,cu_vecb%mem,  &
          m_, &
          nnd, cu_offset%mem, cu_diaA%mem, & 
#else
          vecx,vecb,  &
          m_, &
          nnd, offset, diaA, & 
#endif
          num_neighbor_ranks, neighbor_ranks, margin, ndim,&
          npe_dim, n_grids, offset_dim,&
          div_direc_th, &
          ipreflag, &
          addL, iflagAS,&      
          MAX_NITER, rtolmax, reshistory_DDM, &
          iovlflag, &
          precon_thblock, &
          Neuman_itrmax, &
          independ_nvec, &
          nBlock, &
          iret &
          )

  case(2)
     call parcel_dbicgstab_ddm(&
          MPI_COMM_WORLD, &
#ifdef CUDA_SOLVER
          cu_vecx%mem,cu_vecb%mem,  &
          m_, &
          nnd, cu_offset%mem, cu_diaA%mem, & 
#else
          vecx,vecb,  &
          m_, &
          nnd, offset, diaA, & 
#endif
          num_neighbor_ranks, neighbor_ranks, margin, ndim,&
          npe_dim, n_grids, offset_dim,&
          div_direc_th, &
          ipreflag, &
          addL, iflagAS,&      
          MAX_NITER, rtolmax, reshistory_DDM, &
          iovlflag, &
          precon_thblock, &
          Neuman_itrmax, &
          independ_nvec, &
          nBlock, &
          iret &
          )
  case(3)
     call parcel_dgmres_ddm(&
          MPI_COMM_WORLD, &
#ifdef CUDA_SOLVER
          cu_vecx%mem,cu_vecb%mem,  &
          m_, &
          nnd, cu_offset%mem, cu_diaA%mem, & 
#else
          vecx,vecb,  &
          m_, &
          nnd, offset, diaA, & 
#endif
          num_neighbor_ranks, neighbor_ranks, margin, ndim,&
          npe_dim, n_grids, offset_dim,&
          div_direc_th, &
          ipreflag, &
          addL, iflagAS,&      
          MAX_NITER, rtolmax, reshistory_DDM, &
          iovlflag, &
          precon_thblock, &
          Neuman_itrmax, &
          independ_nvec, &
          nBlock, &
          restart, &
          GMRES_GSflag, &
          iret &
          )
  case(4)
     call parcel_dcagmres_ddm(&
          MPI_COMM_WORLD, &
#ifdef CUDA_SOLVER
          cu_vecx%mem,cu_vecb%mem,  &
          m_, &
          nnd, cu_offset%mem, cu_diaA%mem, & 
#else
          vecx,vecb,  &
          m_, &
          nnd, offset, diaA, & 
#endif
          num_neighbor_ranks, neighbor_ranks, margin, ndim,&
          npe_dim, n_grids, offset_dim,&
          div_direc_th, &
          ipreflag, &
          addL, iflagAS,&      
          MAX_NITER, rtolmax, reshistory_DDM, &
          iovlflag, &
          precon_thblock, &
          Neuman_itrmax, &
          independ_nvec, &
          nBlock, &
          CAGMRES_sstep, &
          CAGMRES_tstep, &
          CAGMRES_basis, &
          CAGMRES_QRflag, &
          iret &
          )
  case(5)
     call parcel_dcbcg_ddm(&
          MPI_COMM_WORLD, &
#ifdef CUDA_SOLVER
          cu_vecx%mem,cu_vecb%mem,  &
          m_, &
          nnd, cu_offset%mem, cu_diaA%mem, & 
#else
          vecx,vecb,  &
          m_, &
          nnd, offset, diaA, & 
#endif
          num_neighbor_ranks, neighbor_ranks, margin, ndim,&
          npe_dim, n_grids, offset_dim,&
          div_direc_th, &
          ipreflag, &
          addL, iflagAS,&      
          MAX_NITER, rtolmax, reshistory_DDM, &
          iovlflag, &
          precon_thblock, &
          Neuman_itrmax, &
          independ_nvec, &
          nBlock, &
          CBCG_kstep, &
          CBCG_Eigenflag, &
          power_method_itrmax, &
          CAARNOLDI_sstep, &
          CAARNOLDI_tstep, &
          CAARNOLDI_basis, &
          CAARNOLDI_QRflag, &
          iret &
          )
  case(51)
     allocate( EigenValueReal(CAARNOLDI_sstep*CAARNOLDI_tstep) )
     allocate( EigenValueImag(CAARNOLDI_sstep*CAARNOLDI_tstep) )
     allocate( EigenVector(CAARNOLDI_sstep*CAARNOLDI_tstep*m) )
     allocate( EigenErrorNorm(CAARNOLDI_sstep*CAARNOLDI_tstep) )

#ifdef CUDA_SOLVER
     write(*,*)"not support GPU CA-ARNOlDI"
#else
     call parcel_dcaarnoldi_ddm(&
          MPI_COMM_WORLD, &
          vecx,vecb,  &
          m_, &
          nnd, offset, diaA, & 
          num_neighbor_ranks, neighbor_ranks, margin, ndim,&
          npe_dim, n_grids, offset_dim,&
          div_direc_th, &
          ipreflag, &
          addL, iflagAS,&      
          iovlflag, &
          precon_thblock, &
          Neuman_itrmax, &
          independ_nvec, &
          nBlock, &
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
#endif
  end select

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  etime=MPI_Wtime()  
  solve_time_DDM=etime-stime

#ifdef CUDA_SOLVER
  call cudaMemcpy_d2h_fp64(cu_vecx,vecx,m,offset_zero,offset_zero)    
  call cudaFree(cu_vecx)
  call cudaFree(cu_vecb)
  call cudaFree(cu_diaA)
  call cudaFree(cu_offset)
#endif

!  call parcel_ddm_free(ipreflag)
  if(myrank.eq.0)then
     if(ityp_solver<10)then
        do i = 1,MAX_NITER
           write(*,*)i,reshistory_DDM(i)
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

  if(myrank.eq.0)then
     write(*,*)"DDM Matrix"
     write(*,*)"   init  time ",init_time_DDM
     write(*,*)"   solve time ",solve_time_DDM
  endif

  deallocate( diaA )
  deallocate( offset )
  deallocate( offset_dim )
  deallocate( reshistory_DDM )
  deallocate( vecx )
  deallocate( vecx0 )
  deallocate( vecb )

  if(ityp_solver<10)then
  else
     deallocate(EigenValueReal)
     deallocate(EigenValueImag)
     deallocate(EigenVector)
     deallocate(EigenErrorNorm)
  endif

  call PARCEL_Free()
  call MPI_Finalize(ierr)

end program main
