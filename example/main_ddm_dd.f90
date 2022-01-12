

program main
  use PARCEL_integer
  use problem
  use ddm_setup
  use mpi

  implicit none

  integer ityp_eq
  integer(kind=intType) gnx,gny,gnz
  integer itrmax
  integer MAX_NITER
  integer(kind=intType) nx,ny,nz
  integer(kind=intType) nx0,ny0,nz0
  integer(kind=intType) n,gn

  integer(kind=intType) m
  integer m_
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

  real*8,allocatable  :: vecx_lo(:)
  real*8,allocatable  :: vecb_lo(:)

  real*8 stime,etime
  real*8 init_time_CRS,solve_time_CRS

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
  real*8, allocatable  :: diaA(:), diaA_lo(:)
  integer, allocatable :: offset(:)
  integer, allocatable :: offset_dim(:,:)

  integer PRECISION_A
  integer PRECISION_b
  integer PRECISION_x
  integer PRECISION_PRECONDITION
  integer nBlock

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

  namelist/input_dd/ &
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION

  namelist/input_ddm/ &
       npe_x,&
       npe_y,&
       npe_z,&
       div_direc_th

  open(33,file='input_namelist')
  read(33,nml=input)
  read(33,nml=input_dd)
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
  num_neighbor_ranks = -1
  call set_neighbor_ranks( &
       npe_x, npe_y, npe_z, rank_x, rank_y, rank_z, &
       num_neighbor_ranks, &
       neighbor_ranks, num_neighbor_ranks_max)


  allocate( diaA(nnd*m) )
  allocate( offset(nnd) )
  allocate( offset_dim(nnd,ndim) )
  allocate( reshistory_DDM(MAX_NITER) )
  allocate( vecx(m) )
  allocate( vecx0(m) )
  allocate( vecb(m) )

  allocate( diaA_lo(nnd*m) )
  allocate( vecx_lo(m) )
  allocate( vecb_lo(m) )

  m_ = m
  diaA = 0.0d0
  vecx = 0.0d0
  vecx0 = 0.0d0
  vecb = 0.0d0
  reshistory_DDM = 0
  offset = 0
  offset_dim = 0

 diaA_lo = 0.0d0
 vecx_lo = 0.0d0
 vecb_lo = 0.0d0

  call make_poisson_3D_DDM( &
       gnx,gny,gnz, &
       ndim, n_grids,  &
       istart_grids, &
       margin, nnd_, m, &
       offset_dim, offset, &
       diaA, vecb, vecx, vecx0,&
       rank_x, rank_y, rank_z, &
       npe_x, npe_y, npe_z )
  select case( solver )
  case(1)
     call parcel_ddcg_ddm(&
          MPI_COMM_WORLD, &
          vecx,vecx_lo, &
          vecb,vecb_lo, &
          m_, &
          nnd, offset, &
          diaA, diaA_lo, &
          num_neighbor_ranks, neighbor_ranks, margin, ndim,&
          npe_dim, n_grids, offset_dim,&
          div_direc_th, &
          ipreflag, &
          addL, iflagAS,&      
          MAX_NITER, rtolmax, reshistory_DDM, &
          iovlflag, &
          precon_thblock, &
          independ_nvec, &
          nBlock, &
          PRECISION_A, &
          PRECISION_b, &
          PRECISION_x, &
          PRECISION_PRECONDITION, &
          iret &
          )
  case(2)
     call parcel_ddbicgstab_ddm(&
          MPI_COMM_WORLD, &
          vecx,vecx_lo, &
          vecb,vecb_lo, &
          m_, &
          nnd, offset, & 
          diaA, diaA_lo, &
          num_neighbor_ranks, neighbor_ranks, margin, ndim,&
          npe_dim, n_grids, offset_dim,&
          div_direc_th, &
          ipreflag, &
          addL, iflagAS,&      
          MAX_NITER, rtolmax, reshistory_DDM, &
          iovlflag, &
          precon_thblock, &
          independ_nvec, &
          nBlock, &
          PRECISION_A, &
          PRECISION_b, &
          PRECISION_x, &
          PRECISION_PRECONDITION, &
          iret &
          )
  case(3)
     call parcel_ddgmres_ddm(&
          MPI_COMM_WORLD, &
          vecx,vecx_lo, &
          vecb,vecb_lo, &
          m_, &
          nnd, offset, & 
          diaA, diaA_lo, &
          num_neighbor_ranks, neighbor_ranks, margin, ndim,&
          npe_dim, n_grids, offset_dim,&
          div_direc_th, &
          ipreflag, &
          addL, iflagAS,&      
          MAX_NITER, rtolmax, reshistory_DDM, &
          iovlflag, &
          precon_thblock, &
          independ_nvec, &
          nBlock, &
          restart, &
          GMRES_GSflag, &
          PRECISION_A, &
          PRECISION_b, &
          PRECISION_x, &
          PRECISION_PRECONDITION, &
          iret &
          )
  end select

#if 0
  call parcel_ddm_main_dd(&
       MPI_COMM_WORLD, ndim, m_, &
       nnd, offset, diaA, & 
       num_neighbor_ranks, neighbor_ranks, margin, ndim,&
       npe_dim, n_grids, offset_dim,&
       div_direc_th, precon_thblock, &
       iovlflag, addL, iflagAS,&
       vecb, vecx, &
       rtolmax, abstolmax,ipreflag,&
       MAX_NITER, reshistory_DDM, &
       restart, &
       independ_nvec, &
       GMRES_GSflag, &
       ityp_solver, &
       PRECISION_KRYLOV, & 
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION, &
       diaA_lo,vecb_lo,vecx_lo &
       )
#endif

  if(myrank.eq.0)then
     do i = 1,MAX_NITER
        write(*,*)i,reshistory_DDM(i)
     enddo
  endif
  deallocate( diaA )
  deallocate( offset )
  deallocate( offset_dim )
  deallocate( reshistory_DDM )
  deallocate( vecx )
  deallocate( vecx0 )
  deallocate( vecb )


  call MPI_Finalize(ierr)

end program main
