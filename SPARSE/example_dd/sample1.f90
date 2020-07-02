program main
! ---------------------------------------------------------
!   test program for bi-cgstab/gmres(m)/tfqmr routines
!
!   coefficient matrix is block-row-partitioned.
!   rhs and result vecotors are not partitioned.
!   numerical result is checked with exact algebraic solution.
!
! ---------------------------------------------------------
!  Internal varieble
!
!     matA_cmpr           : Matrix A of CRS format          (Max elemant:NDIM*NNZR)
!     matA_colindex_coor  : Columns index of Mat, A  (CRS)  (Max elemant:NDIM*NNZR)
!     matA_rowpointer_cmpr: Row pointer of Mat. A (CRS)     (Max elemant:NDIM+1)
!
!     matA_coor           : Matrix A of Coord form.         (Max elemant:NDIM*NNZR)
!     matA_colindex_coor  : Columns index of Mat. A (COORD) (Max elemant:NDIM*NNZR)
!     matA_rowindex_coor  : Columns index of Mat. A (COORD) (Max elemant:NDIM*NNZR)
!
!     X(NDIM)      : Resultant vector
!     vecB(NDIM)   : Initial RHS
!     X_exact(NDIM)   : Exactrical solutions of test matrix
!
!  Input data from input file
!     ityp_eq      : Sample problem type
!                  :    1: 2D Poisson's equation
!                  :    2: 2D Laplace's equation
!                       3: random matrix
!                       4: glued Wilkinson matrix
!                       5: (reduced) Frank matrix
!
!     rtolmax          : Convergence tolerance
!     ipreflag         : preconditioning option
!                  :    0: without
!                       1: point Jacobi
!                       2: block Jacobi
!                       3: additive Schwarz
!
!     ityp_solver  :    1: bi-cgstab, 2: gmres(m), 3: tfqmr, 4: gpbi-cg
!
!  Others
!     MAX_NITER    : Maximum iterration number
!     myproc       : Instance number od MPI
!     NPES         : Number of PE
!     NDIM         : dimension of matrix, vector
!     Num_nonZERO  : Non zero elements number
!
! ---------------------------------------------------------
  use mpi
  use omp_lib
  implicit none
!  include 'JAERI_iles_f90.h'
!   include 'mpif.h'
!   include "omp_lib.h"
!
  integer MAX_NDIM, NNZR, MAX_NITER
!
! -- Parameters for the problems ----
  integer ndim, ityp_eq, ityp_solver, ipreflag, addL
  integer npoints_along_x
  integer ILU_method, nres, restart
  integer blocksize, nrepeat
  real*8  glue, rtolmax

! -- info of matrix
  integer :: matA_format = 1               ! 1: CRS,   2:  DIA
  logical :: matA_is_symmetric = .false.
  integer iflagAS
  integer GMRES_GSflag
  integer CBCG_kstep
  integer CBCG_Eigenflag
  integer power_method_itrmax
  integer CAARNOLDI_sstep
  integer CAARNOLDI_tstep
  integer CAARNOLDI_basis
  integer CAARNOLDI_QRflag
  integer CAGMRES_sstep
  integer CAGMRES_tstep
  integer CAGMRES_basis
  integer CAGMRES_QRflag

!
! -- COO format --
  integer, allocatable :: matA_colindex_coor(:), matA_rowindex_coor(:)
  real*8, allocatable :: matA_coor(:)

! -- CRS format ---
  integer num_nonzero
  integer, allocatable :: matA_colindex_cmpr(:), matA_rowpointer_cmpr(:)
  integer, allocatable :: matA_coor_to_cmpr_table(:)
  real*8, allocatable :: matA_cmpr(:), matA_cmpr_lo(:)

! -- DIA format --
  integer nnd
  integer, allocatable :: matA_ioff_dia(:)
  real*8, allocatable :: matA_dia(:), matA_dia_lo(:)

! -- vectors --
  real*8, allocatable :: X(:), X_exact(:), vecB(:)
!  real*8, allocatable :: X_lo(:), vecB_lo(:)
!
! -- converence --
  real*8 eps
  real*8, allocatable :: his_resid(:)

! -- parallelization --
  integer :: istatus (mpi_status_size)
  integer :: npes, myproc, icomm, ndist, istart, iend

! -- openmp --
  integer :: num_threads

! -- timing --
  integer iloop, num_loops
  real*8 tottims, tottime, tottim, caltims, caltime, caltim
  real*8 tc0, tc1, tc2, tw0, tw1, tw2, tt(2), ttmp(2)

! -- others
  integer num_nonzero_out
  integer ierr, niter, iret
  integer I, num, j, is1, ie1
  real*8, allocatable :: Y(:), X_tmp(:)

  integer iovlflag

  real*8, allocatable :: x_lo(:)
  real*8, allocatable :: b_lo(:)
!  real*8, allocatable :: matA_cmpr_lo(:)
  integer PRECISION_A
  integer PRECISION_b
  integer PRECISION_x
  integer PRECISION_PRECONDITION
  character(len=10) ::argc
  character(len=10) :: set_file_num
  character(len=40) :: file_name
  integer imode_flg,input_num


  !call getarg(1,argc)
  !set_file_num=trim(argc)


  namelist/input/ &
       ityp_eq, &
       npoints_along_x, &
       MAX_NITER, &
       rtolmax, &
       ipreflag, &
       ILU_method, &
       iovlflag, &
       matA_format, &
       ityp_solver, &
       addL, &
       iflagAS, &
       GMRES_GSflag, &
       restart, &
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
       CAGMRES_QRflag, &
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION


! --------------
!  initialization for mpi
! --------------

  call MPI_Init( ierr )
  call MPI_Comm_size( MPI_COMM_WORLD, npes, ierr )
  call MPI_Comm_rank( MPI_COMM_WORLD, myproc, ierr )

  icomm = MPI_COMM_WORLD

!$omp parallel
  num_threads = omp_get_num_threads()
!$omp end parallel

! --------------
!  read parameters
! --------------
  call set_max_values

! ------
#if 1
! -----
! デフォルト設定
  ityp_eq                = 1
  npoints_along_x        = 10
  MAX_NITER              = 100
  rtolmax                    = 1.0-8
  ipreflag                   = 0
  ILU_method              = 1
  iovlflag                = 1
  matA_format            = 1
  ityp_solver            = 3
  addL                    = 0
  iflagAS                = 4
  GMRES_GSflag           = 1
  restart                = 8
  CBCG_kstep             = 4
  CBCG_Eigenflag         = 2
  power_method_itrmax    = 100
  CAARNOLDI_sstep        = 2
  CAARNOLDI_tstep        = 5
  CAARNOLDI_basis        = 1
  CAARNOLDI_QRflag       = 3
  CAGMRES_sstep          = 2
  CAGMRES_tstep          = 4
  CAGMRES_basis          = 1
  CAGMRES_QRflag         = 3
  PRECISION_A            = 1
  PRECISION_b            = 1
  PRECISION_x            = 1
  PRECISION_PRECONDITION = 1
! -----
  !if ( myproc == 0 ) then
  
  !read(set_file_num,*,iostat=imode_flg) input_num
  !if (imode_flg==0) then
  !  if ((input_num>=1).and.(input_num<=35)) then
  !    file_name='./templete/templete_input/input_sample_poisson_namelist'//set_file_num
  !  else
  !    print *,"invalid input file number"
  !    stop
  !  endif
  !else
    file_name='./input_sample_poisson_namelist'
  !endif
 
  open(33,file=file_name)
  read(33,nml=input)
  close(33)

  if ((ipreflag==3) .and.(ILU_method==0)) then
    write(*,*)"Additive schwartz doesn't support ILU precondition. Switch ILU to DILU condition 1."
    ILU_method =1
  endif
  if ((ipreflag==3).and.(addL==0)) then
    ipreflag=2
  endif



  select case( ityp_eq )
    case(1)
       ! 2d poisson
       ndim=npoints_along_x*npoints_along_x
    case(6)
       ! 3d poisson
       ndim=npoints_along_x*npoints_along_x*npoints_along_x
    end select
!  endif
  num_loops = 1

#else
  if ( myproc == 0 ) then

     tottims = MPI_WTIME()

     call set_equation_to_be_solved( ityp_eq )
     call set_dimnesions_of_matrix( ityp_eq, npoints_along_x, ndim, &
       &                            blocksize, nrepeat, glue, NNZR )

     call set_max_iterations( MAX_NITER )
     call set_torelance( rtolmax )
     call set_preconditioner_type( ipreflag )

     call set_ilumthod_for_precon( ipreflag, ILU_method, addL )
!
     call set_solver_type( ityp_solver )
     call set_restart_param( ityp_solver, restart )

     call set_matrix_format( matA_format )
!!!     call set_lowerbit_param( ipreflagc_matA, ipreflagc_vecB, ipreflagc_precon )
     call set_communication_hiding( iovlflag )

!!     read(*,*) num_loops
     num_loops = 1
  end if
  call broadcast_params( icomm, ityp_eq, npoints_along_x, ndim, &
       &                 blocksize, nrepeat, glue, &
       &                 MAX_NITER, rtolmax, ipreflag, ILU_method, addL, ityp_solver, restart, &
       &                 NNZR, matA_format, iovlflag)
!       &                 NNZR, matA_format )
  call mpi_bcast( num_loops, 1, mpi_integer, 0, icomm, ierr )
#endif

! -------------------------------
!  call JAERI_SET_LOCAL_ILU_METHOD( ILU_method )

!  call alloc_arrays( ndim, NNZR, npes, num_nonzero, nres )

  call set_istart_iend( ndim, npes, myproc, istart, iend )

  call alloc_arrays_lm( ndim, NNZR, npes, num_nonzero, nres, istart, iend)

! --------------
!  set coefficient matrix, rhs vector, initial guess,  and exact solution
! --------------
  select case( ityp_eq )
  case (1)
     call set_matrix_poisson_2D_lm( npoints_along_x, ndim, num_nonzero, &
          &                      matA_coor, matA_rowindex_coor, matA_colindex_coor, &
          &                      vecB, X, X_exact, num_nonzero_out, &
          &                      istart, iend, matA_is_symmetric )
  case(2)
     call set_matrix_diffusion_2D( npoints_along_x, ndim, num_nonzero, &
          &                        matA_coor, matA_rowindex_coor, matA_colindex_coor, &
     &                             vecB, X, X_exact, num_nonzero_out, &
     &                             istart, iend, matA_is_symmetric )

  case(3)             ! use with caution ( *if ndim is too samall, the program hangs )
     call set_matrix_random( NNZR, ndim, num_nonzero, &
          &                  matA_coor, matA_rowindex_coor, matA_colindex_coor, &
          &                  vecB, X, X_exact, num_nonzero_out, &
          &                  istart, iend, matA_is_symmetric )

  case(4)
     call set_matrix_glued_wilkinson( blocksize, nrepeat, glue,  &
     &                                ndim, num_nonzero, &
     &                                matA_coor, &
     &                                matA_rowindex_coor, matA_colindex_coor, &
     &                                vecB, X, X_exact, num_nonzero_out, &
     &                                istart, iend, matA_is_symmetric )

  case(5)
     call set_matrix_reduced_Frank( NNZR, ndim, num_nonzero, &
     &                              matA_coor, &
     &                              matA_rowindex_coor, matA_colindex_coor, &
     &                              vecB, X, X_exact, num_nonzero_out, &
     &                              istart, iend, matA_is_symmetric )

  case (6)
!     write(*,*) ndiv_z, ndiv_y, ndiv_z
!     stop
     call set_matrix_poisson_3D_lm( npoints_along_x, ndim, num_nonzero, &
          &                      matA_coor, matA_rowindex_coor, matA_colindex_coor, &
          &                      vecB, X, X_exact, num_nonzero_out, &
          &                      istart, iend, matA_is_symmetric )

  case default

  end select

  num_nonzero = num_nonzero_out

  ndist = iend - istart +1

  select case ( matA_format )
  case (1)        ! CRS
!     call alloc_matrixA_crs( num_nonzero )
     call alloc_matrixA_crs( num_nonzero, ndist )

     CALL JAERI_CO2CRTB( ndim, num_nonzero, istart, NDIST, &
          &          matA_rowindex_coor, matA_colindex_coor, &
          &          matA_colindex_cmpr, matA_rowpointer_cmpr, &
          &          matA_coor_to_cmpr_table )
     CALL JAERI_CO2CR( num_nonzero, matA_coor, matA_cmpr, matA_coor_to_cmpr_table )

  case (2)   ! DIA format
     call co2dia_count_nnd( ndim, num_nonzero, &
          &                    matA_rowindex_coor, &
          &                    matA_colindex_coor, istart, iend, nnd )
     call alloc_matrixA_dia( nnd, ndist )

     call co2dia( ndim, num_nonzero, &
          &          matA_rowindex_coor, matA_colindex_coor, &
          &          matA_coor, istart, iend, &
          &          ndist, nnd, matA_ioff_dia, matA_dia )
  end select

! -----------------------
  niter = MAX_NITER;   eps = rtolmax

  call MPI_Barrier( icomm, ierr)
  
  call cputime(tc1);   call wallclock(tw1)
  caltims = MPI_WTIME()

! #ifdef BX_GLOBAL
!   is1=1
!   ie1=ndim
! #else
  is1=istart
  ie1=iend
!#endif

  allocate( X_tmp(is1:ie1) );

! -------------------------------------------
!$omp parallel do
     Do j=is1, ie1
        X_tmp(j) = X(j)
     End Do
!$omp end parallel do
! -------------------------------------------

  Do iloop=1, num_loops
! -------------------------------------------
!$omp parallel do
     Do j=is1, ie1
        X(j) = X_tmp(j)
     End Do
!$omp end parallel do
     niter = max_niter;     eps = rtolmax
! -------------------------------------------

! ----------
! Solving the problem here.
! ----------
     allocate(x_lo(ndist))
     allocate(b_lo(ndist))
     allocate(matA_cmpr_lo(num_nonzero))
     allocate(matA_dia_lo(nnd*ndist))
     matA_dia_lo=0.0d0
     matA_cmpr_lo=0.0d0
     select case( ityp_solver )
     case(1)
! ---------
!c  CG
!c ---------
        if ( matA_format == 1 ) then         ! CRS
           call  parcel_ddcg(&
                icomm, &
                X(is1:ie1)   ,x_lo,&
                vecB(is1:ie1),b_lo,  &
                ndist,ndim,num_nonzero,istart, &
                matA_cmpr,matA_cmpr_lo,matA_rowpointer_cmpr,matA_colindex_cmpr, &
                ipreflag,ILU_method,addL,iflagAS, &
                niter,eps,&
                his_resid,&
                iovlflag, iret, &
                PRECISION_A, &
                PRECISION_b, &
                PRECISION_x, &
                PRECISION_PRECONDITION &
                )
        else                                 ! DIA
          call parcel_ddcg_dia( &
                icomm,&
                X(is1:ie1)   ,x_lo, &
                vecB(is1:ie1),b_lo,  &
                ndist,ndim, istart, &
                matA_dia, matA_dia_lo, matA_ioff_dia, nnd, &
                ipreflag,ILU_method,addL,iflagAS, &
                niter,eps,&
                his_resid,&
                iovlflag, iret, &
                PRECISION_A, &
                PRECISION_b, &
                PRECISION_x, &
                PRECISION_PRECONDITION &
                )
        endif
     case(2)

! ---------
!c  BI-CGSTAB
!c ---------
        if ( matA_format == 1 ) then         ! CRS
           call parcel_ddbicgstab( &
                icomm, &
                X(is1:ie1),x_lo,&
                vecB(is1:ie1),b_lo,  &
                ndist,ndim,num_nonzero,istart, &
                matA_cmpr,matA_cmpr_lo,matA_rowpointer_cmpr,matA_colindex_cmpr, &
                ipreflag,ILU_method,addL,iflagAS, &
                niter,eps,&
                his_resid,&
                iovlflag, iret,&
                PRECISION_A, &
                PRECISION_b, &
                PRECISION_x, &
                PRECISION_PRECONDITION &
                )
        else                                 ! DIA
           call parcel_ddbicgstab_dia( &
                icomm,&
                X(is1:ie1),x_lo,&
                vecB(is1:ie1),b_lo,  &
                ndist,ndim, istart, &
                matA_dia, matA_dia_lo, matA_ioff_dia, nnd, &
                ipreflag,ILU_method,addL,iflagAS, &
                niter,eps,&
                his_resid,&
                iovlflag, iret,&
                PRECISION_A, &
                PRECISION_b, &
                PRECISION_x, &
                PRECISION_PRECONDITION )
        endif
     case(3)
!c ---------
!c GMRES(m)
!c ---------
        if ( matA_format == 1 ) then         ! CRS
           call parcel_ddgmres( &
                icomm, &
                X(is1:ie1),x_lo,&
                vecB(is1:ie1),b_lo,  &
                ndist,ndim,num_nonzero,istart, &
                matA_cmpr,matA_cmpr_lo,matA_rowpointer_cmpr,matA_colindex_cmpr, &
                ipreflag,ILU_method,addL,iflagAS, &
                niter,eps,&
                his_resid,&
                iovlflag, iret,&
                restart,GMRES_GSflag,&
                PRECISION_A, &
                PRECISION_b, &
                PRECISION_x, &
                PRECISION_PRECONDITION &
                )
        else
           call parcel_ddgmres_dia( &
                icomm,&
                X(is1:ie1),x_lo,&
                vecB(is1:ie1),b_lo,  &
                ndist,ndim, istart, &
                matA_dia, matA_dia_lo, matA_ioff_dia, nnd, &
                ipreflag,ILU_method,addL,iflagAS, &
                niter,eps,&
                his_resid,&
                iovlflag, iret,&
                restart,GMRES_GSflag,&
                PRECISION_A, &
                PRECISION_b, &
                PRECISION_x, &
                PRECISION_PRECONDITION )
        endif
     end select
    End Do
  deallocate( X_tmp )

  caltime = MPI_WTIME()

  call cputime(tc2);    call wallclock(tw2)
  ttmp(1) = tc2 - tc1;  ttmp(2) = tw2 - tw1
  call MPI_Reduce( ttmp, tt, 2, MPI_DOUBLE_PRECISION, MPI_MAX, 0, icomm, ierr )

  tt = tt /dble(num_loops)           ! per loop

! ----------
! Finalization
! ----------
! #ifndef BX_GLOBAL
!   allocate( Y(ndim) )
!   Do i=1, istart -1
!      X(i) = 0.0d0
!   End do
!   Do i=iend+1, ndim
!      X(i) = 0.0d0
!   End do
!   call mpi_allreduce( X, Y, ndim, mpi_double_precision, mpi_sum, icomm, ierr )
!   X = Y
!   deallocate( Y )
! #endif

  if ( myproc==0 )  then
     call print_info;   
     call write_residuals( min(nres,niter) )

! ----------
!   Time calc.
! ----------
     tottime = MPI_WTIME()
     caltim = caltime -caltims;   tottim = tottime -tottims

     call print_time
!     call calc_diffrence_from_exact_sol( ndim, X, X_exact )

  end if

  call mpi_barrier( icomm,ierr )
  call calc_diffrence_from_exact_sol_lm( ndim, X, X_exact, is1, ie1, myproc, icomm)
  call mpi_barrier( icomm,ierr )
  
  deallocate(x_lo)
  deallocate(b_lo)
  deallocate(matA_cmpr_lo)

  call MPI_Finalize(ierr)

contains

  subroutine set_max_values
! -----------------
!   MAX_NDIM  : maximum number of varieble of LE.
!   NNZR  : maximum number of columns with CR formality : 5
! --
!   MAX_NITER : maximum number of iteration steps.
! -----------------
!
    MAX_NDIM = 225000000;     NNZR = 7;        MAX_NITER = 1000

  end subroutine set_max_values

  subroutine alloc_matrixA_crs( num_nonzero, ndist )
    integer, intent(in) :: num_nonzero
    integer, intent(in) :: ndist
!
    allocate( matA_cmpr( num_nonzero ) ); matA_cmpr = 0.0d0
    allocate( matA_colindex_cmpr(num_nonzero) );   matA_colindex_cmpr = 0
!    allocate( matA_rowpointer_cmpr(ndim+1) ); matA_rowpointer_cmpr = 0
    allocate( matA_rowpointer_cmpr(ndist+1) ); matA_rowpointer_cmpr = 0
    allocate( matA_coor_to_cmpr_table(num_nonzero) ); matA_coor_to_cmpr_table = 0

  end subroutine alloc_matrixA_crs

  subroutine alloc_matrixA_dia( nnd, ndist )
    integer, intent(in) :: nnd, ndist
!
    allocate( matA_dia( nnd*ndist ) ); matA_dia = 0.0d0
    allocate( matA_ioff_dia( nnd ) ); matA_ioff_dia = 0
!
  end subroutine alloc_matrixA_dia

  subroutine alloc_arrays( ndim, NNZR, npes, num_nonzero, nres )
    integer, intent(in) :: ndim, NNZR, npes
    integer, intent(out) :: num_nonzero, nres
! -----------------------
! local
    integer nun_nonzero, n1

! begin
    num_nonzero = ndim *NNZR           ! global
    num_nonzero = num_nonzero / npes
!
!    n1 = ndim /npes
!    if ( ndim /= n1*npes ) n1 = n1 +1
!    num_nonzero = n1 *NNZR
!
    allocate( matA_coor(num_nonzero) );      matA_coor = 0.0d0
    allocate( matA_colindex_coor(num_nonzero) );   matA_colindex_coor = 0
    allocate( matA_rowindex_coor(num_nonzero) );   matA_rowindex_coor = 0

    allocate( vecB(ndim) ); vecB = 0.0d0
    allocate( X(ndim) );   X = 0.0d0
    allocate( X_exact(ndim) ); X_exact = 0.0d0
!
    nres = MAX_NITER
    allocate(  his_resid(nres) ); his_resid = 0.0d0

  end subroutine alloc_arrays
  
  subroutine alloc_arrays_lm( ndim, NNZR, npes, num_nonzero, nres, is, ie )
    integer, intent(in) :: ndim, NNZR, npes
    integer, intent(out) :: num_nonzero, nres
    integer, intent(in) :: is, ie
! -----------------------
! local
    integer nun_nonzero, n1
! begin
    num_nonzero = ndim *NNZR           ! global
    num_nonzero = num_nonzero / npes
!
!    n1 = ndim /npes
!    if ( ndim /= n1*npes ) n1 = n1 +1
!    num_nonzero = n1 *NNZR
!
    allocate( matA_coor(num_nonzero) );      matA_coor = 0.0d0
    allocate( matA_colindex_coor(num_nonzero) );   matA_colindex_coor = 0
    allocate( matA_rowindex_coor(num_nonzero) );   matA_rowindex_coor = 0

    allocate( vecB(is:ie) ); vecB = 0.0d0
    allocate( X(is:ie) );   X = 0.0d0
    allocate( X_exact(is:ie) ); X_exact = 0.0d0
!
    nres = MAX_NITER
    allocate(  his_resid(nres) ); his_resid = 0.0d0

  end subroutine alloc_arrays_lm

  
  subroutine set_equation_to_be_solved( ityp_eq )
    integer, intent(out) :: ityp_eq
! begin
    write(*,*) 'Which matrix is solved [1-3]?  '
    write(*,*) '    ( 1: Poisson eq. 2D, ', &
         &           '6: Poisson eq. 3D)'
    read(*,*) ityp_eq
!
    select case( ityp_eq )
    case(1)
       write(*,*) 'Poisson Equation (2D) is selected'
    ! case(2)
    !    write(*,*) 'Steady equation of diffusion is selected'
    ! case(3)
    !    write(*,*) 'random matrix is selected'
    ! case(4)
    !    write(*,*) 'Glued Wilkinson matrix is selected'
    ! case(5)
    !    write(*,*) '(Reduced) Frank matrix is selected'
    case(6)
       write(*,*) 'Poisson Equation (3D) is selected'
    case default
       write(*,*) 'please enter 1 or 6'; stop

    end select

    return

  end subroutine set_equation_to_be_solved
  
 subroutine set_dimnesions_of_matrix( ityp_eq, npoints_along_x, ndim, &
       &                               blocksize, nrepeat, glue, NNZR )
    integer, intent(in) :: ityp_eq
    integer, intent(out) :: npoints_along_x, ndim
    integer, intent(out) :: blocksize, nrepeat, NNZR
    real*8, intent(out) :: glue
!
    integer ntmp
!
    npoints_along_x = 0;  ndim = 0
    blocksize = 0;  nrepeat = 0;  glue = 0.0d0
!    seed = 1
!
    select case( ityp_eq )
    case (1,2)
       write(*,*) 'Number of grid points along one direction ?'
       read(*,*) npoints_along_x
       if ( npoints_along_x <= 0 ) stop "only positive values are allowed"

       ndim = npoints_along_x**2

    case(3)
       write(*,*) 'Dimension of matrix ? '
       read(*,*) ndim
       if ( ndim <= 0 ) stop "only positive values are allowed"

       write(*,*) 'Number of nonzero elements in each row ? '
       read(*,*) ntmp
       if ( ntmp <= 0 ) stop "only positive values are allowed"
       NNZR = ntmp
       if ( NNZR > ndim ) stop "NNZR is too large"

       write(*,*) 'NNZR is set to ', NNZR

!       write(*,*) 'specify the seed for random number '
!       read(*,*) seed
!       if ( seed <= 0 ) stop "only positive values are allowed for seed"
!
!       write(*,*) 'seed is set to ', seed

    case(4)
       write(*,*) 'Number of blocksize ? '
       read(*,*) blocksize
       if ( blocksize <= 0 ) stop "only positive values are allowed for blocksize"

       write(*,*) 'Number of dimension of repeating units ? '
       read(*,*) nrepeat
       if ( nrepeat <= 0 ) stop "only positive values are allowed for nrepeat"

       write(*,*) 'specify the glue value '
       read(*,*) glue
       if ( glue <= 0 ) stop "only positive values are allowed for glue"

       ndim = blocksize *nrepeat

   case(5)
       write(*,*) 'Dimension of matrix ? '
       read(*,*) ndim
       if ( ndim <= 0 ) stop "only positive values are allowed"

       write(*,*) 'Number of nonzero elements in each row ? '
       read(*,*) ntmp
       if ( ntmp <= 0 ) stop "only positive values are allowed"
       NNZR = ntmp
       if ( NNZR > ndim ) stop "NNZR is too large"

       write(*,*) 'NNZR is set to ', NNZR

    case (6)
       write(*,*) 'Number of grid points along one direction ?'
       read(*,*) npoints_along_x
       if ( npoints_along_x <= 0 ) stop "only positive values are allowed"

       ndim = npoints_along_x**3

    end select
!
    if ( ndim > MAX_NDIM ) then
       stop "ndim is too large"
    else
       write(*,*)' ndim is set to ', ndim
    endif

  end subroutine set_dimnesions_of_matrix

  subroutine set_max_iterations( MAX_NITER_new )
    integer, intent(out) :: MAX_NITER_new
!
    integer n1
!
    write(*,*) 'Max number of iterations ?'
    read(*,*) n1
!
    if ( n1 > 0 .and. n1 < MAX_NITER ) then
       MAX_NITER_new = n1
       write(*,*) 'Max_NITER is set to ', n1
    else
       MAX_NITER_new = MAX_NITER
    endif

  end subroutine set_max_iterations

  subroutine set_torelance( rtolmax )
    real*8, intent(out) :: rtolmax
!
    real*8 c1
!
    write(*,*) 'Convergence tolerance ( 10**n ); specify n '
    read(*,*) c1
    rtolmax = 10.0**c1

    write(*,*) 'rtolmax is set to ',rtolmax

  end subroutine set_torelance

  subroutine set_preconditioner_type( ipreflag )
    integer, intent(out) :: ipreflag

    write(*,*) 'Preconditioner Type(0:None 1:Jacobi 2:Additive Schwarz) ?'
    read(*,*) ipreflag

    select case( ipreflag )
    case(0)
       write(*,*) 'Precon is not used.'
    case(1)
       write(*,*) 'Point Jacobi is used.'
    case(2)
       write(*,*) 'Block Jacobi is used.'
    case(3)
       write(*,*) 'Additive schwarz is used.'
    case default
       write(*,*) 'please enter 0-3'; stop
    end select

  end subroutine set_preconditioner_type

  subroutine  set_ilumthod_for_precon( ipreflag, ILU_method, addL )
    integer, intent(in) :: ipreflag
    integer, intent(out) :: ILU_method, addL
    
    addL = 0;     ILU_method = 1
    if ( ipreflag == 3 ) then
       write(*,*) 'ILU_method, addL ? '
       read(*,*) ILU_method, addL
       
       write(*,*) 'ILU_method = ', ILU_method
       write(*,*) 'addL  = ', addL
    else
       read(*,*)
    endif
    
  end subroutine set_ilumthod_for_precon

  subroutine set_solver_type( ityp_solver )
    integer, intent(out) :: ityp_solver
! 
    write(*,*) 'Which solver is used [1-2] ?'
    write(*,*) '       ( 1: bi-cgstab, 2: gmres(m) )'
    read(*,*)   ityp_solver
!
    select case( ityp_solver )
    case(1)
       write(*,*) 'Bi-CGSTAB    is selected.'
    case(2)
       write(*,*) 'GMRES(m)     is selected.'
    case(3)
       write(*,*) 'CG           is selected.'
    case(4)
       write(*,*) 'CBCG         is selected.'
    case(5)
       write(*,*) 'CAGMRES(s,t) is selected.'
    case default
       write(*,*) 'please enter 1-2'; stop

    end select

  end subroutine set_solver_type

  subroutine set_matrix_format( matA_format )
    integer, intent(out) :: matA_format
!
    write(*,*) 'matrix format ( 1: CRS, 2: DIA ) ?'
    read(*,*) matA_format
!
    select case( matA_format )
    case(1)
       write(*,*) 'CRS format is used.'
    case(2)
       write(*,*) 'DIA format is used'
    case default
       write(*,*) 'please enter 1 or 2'; stop
    end select

  end subroutine set_matrix_format

  subroutine set_communication_hiding( iovlflag )
    integer, intent(out) :: iovlflag
!
    write(*,*) 'Communication hiding mode ( 0:none 1: once, 2: each process ) ?'
    read(*,*) iovlflag
!
    select case( iovlflag )
    case(0)
       write(*,*) 'non Communication Hiding or single process'
    case(1)
       write(*,*) 'once'
    case(2)
       write(*,*) 'each process'
    case default
       write(*,*) 'please enter 0, 1 or 2 '; stop
    end select

  end subroutine set_communication_hiding 
  
  subroutine set_restart_param( ityp_solver, restart )
    integer, intent(in) :: ityp_solver
    integer, intent(out) :: restart
!
    restart = 0
    if ( ityp_solver == 2 ) then
       write(*,*) ' value of restart ? '
       read(*,*) restart
    else
       read(*,*)
    endif

  end subroutine set_restart_param

  subroutine broadcast_params( icomm, ityp_eq, npoints_along_x, ndim, &
       &                       blocksize, nrepeat, glue, &
       &                       MAX_NITER, rtolmax, ipreflag, ILU_method, addL, &
       &                       ityp_solver, restart, NNZR, &
       &                       matA_format, iovlflag)
    integer, intent(in) :: icomm, ityp_eq,  npoints_along_x, ndim
    integer, intent(in) :: blocksize, nrepeat, MAX_NITER, NNZR
    integer, intent(in) :: ipreflag, ILU_method, addL, ityp_solver, restart
    real*8, intent(in) :: rtolmax, glue
    integer, intent(in) :: matA_format
    integer, intent(in) :: iovlflag

    integer info
!
    CALL mpi_barrier( icomm, info )
    call mpi_bcast( ityp_eq, 1, mpi_integer, 0, icomm, info )

    call mpi_bcast( npoints_along_x, 1, mpi_integer, 0, icomm, info )
    call mpi_bcast( ndim, 1, mpi_integer, 0, icomm, info )
    call mpi_bcast( blocksize, 1, mpi_integer, 0, icomm, info )
    call mpi_bcast( nrepeat, 1, mpi_integer, 0, icomm, info )

    call mpi_bcast( MAX_NITER, 1, mpi_integer, 0, icomm, info )

    call mpi_bcast( glue, 1, mpi_double_precision, 0, icomm, info )
    call mpi_bcast( rtolmax, 1, mpi_double_precision, 0, icomm, info )

    call mpi_bcast( ipreflag, 1, mpi_integer, 0, icomm,  info )

    call mpi_bcast( ILU_method, 1, mpi_integer, 0, icomm, info )
    call mpi_bcast( addL, 1, mpi_integer, 0, icomm, info )
    call mpi_bcast( ityp_solver, 1, mpi_integer, 0, icomm, info )
    call mpi_bcast( restart, 1, mpi_integer, 0, icomm, info )

    call mpi_bcast( NNZR, 1, mpi_integer, 0, icomm, info )
!    call mpi_bcast( seed, 1, mpi_integer, 0, icomm, info )

   call mpi_bcast( matA_format, 1, mpi_integer, 0, icomm, info )
   call mpi_bcast( iovlflag, 1, mpi_integer, 0, icomm, info )

  end subroutine broadcast_params

  subroutine set_istart_iend( ndim, npes, myproc, istart, iend )
    integer, intent(in) :: ndim, npes, myproc
    integer, intent(out) :: istart, iend

    integer iw1, iw2

    iw1 = ndim / npes
    iw2 = mod( ndim, npes )
    istart = myproc *iw1 + 1 + min( myproc,iw2 )
    iend   = istart + iw1 - 1
    if (iw2 > myproc) iend = iend + 1

  end subroutine set_istart_iend

  subroutine print_info
    write(*, *)
    if (ityp_solver == 1) then
       !write(*,'(1x,a)') '== Sample Program for Bi-CGSTAB Routines =='
       write(*,'(1x,a)') '== Sample Program for CG Routines =='
    else if (ityp_solver == 2) then
       write(*,'(1x,a)') '== Sample Program for Bi-CGSTAB Routines =='
       !write(*,'(1x,a)') '== Sample Program for GMRES(m) Routines =='
    else if (ityp_solver == 3) then
       write(*,'(1x,a)') '== Sample Program for GMRES(m) Routines =='
       !write(*,'(1x,a)') '== Sample Program for TFQMR Routines =='
    else if (ityp_solver == 4) then
       write(*,'(1x,a)') '== Sample Program for GPBi-CG Routines =='
    end if
    
    write(*, *)
    write(*,'(2x,a)') '--- global infomation ---'
    write(*,'(2x,a,i3)') '# of PE:                npes=', npes
    write(*,'(2x,a,i3)') '# of threads:       nthreads=', num_threads
    write(*,'(2x,a,i6)') 'matrix dimension:       ndim=', ndim
    
    write(*, *)
    write(*,'(2x,a)') '--- conditions ---'
    write(*,'(2x,a,i1)')      'preconditioner:        ipreflag=', ipreflag
    
    if ( ipreflag == 3 ) then
       write(*,'(2x,a,i5)') 'overlap of as:         addL=', addL
    endif
    
    write(*,'(2x,a,i6)')  'max # of iterations:   MAX_NITER =', MAX_NITER
    
    if ( ityp_solver ==2 ) then
       write(*,'(2x,a,i3)') 'restart number:      restart=',restart
    endif
    
    write(*,'(2x,a,1pe15.6)') 'convergence condition: rtolmax=', rtolmax
    write(*, *)
    write(*,'(2x,a)') '--- results ---'
    write(*,'(2x,a,i3)')  'return code:           iret=', iret
    write(*,'(2x,a,i6)')  '# of iterations:       niter=', niter
    write(*,'(2x,a,1pe15.6)') 'relative residual:     eps=', eps

  end subroutine print_info

  subroutine print_time
    write(*,*)
    write(*,'(2x,a)') '--- timing  ---'

    write(*,'(2x,a,1pe15.6)') 'cputime(sec):   ', tt(1)
    write(*,'(2x,a,1pe15.6)') 'wallclock(sec): ', tt(2)

    write(*,*)
    write(*,*) '<MPI_WTIME>'
    write(6,820) caltim
820 format('  PCG(SUBROUTINE)...time :',e15.6)
    write(6,825) tottim
825 format('  PCG(PROGRAM)......time :',e15.6)


  end subroutine print_time

  subroutine calc_diffrence_from_exact_sol( ndim, X, X_exact )
    integer, intent(in) :: ndim
    real*8, intent(in) :: X(ndim), X_exact(ndim)
!
    integer i
    real*8 znmx, zdmx, zerr
!
    znmx = -1d0;    zdmx = -1d0
    do i=1,ndim
       znmx = max( znmx, abs(X(i)-X_exact(i)) )
       zdmx = max( zdmx, abs(X_exact(i)) )
    enddo
    zerr=znmx/zdmx
    
    write(*,*)
    write(*,'(2x,a)') '--- check result  ---'
    write(*,'(2x,a)') &
         &           'relative error of the numerical solution'  &
         &                   // ' with algebraic solution'
    
    write(*,*) '    max(|numerical-exact|)/max(|exact|) =',zerr
    write(*,*) '    max(|numerical-exact|)              =',znmx
    write(*,*) '                           max(|exact|) =',zdmx
    write(*,*)

  end subroutine calc_diffrence_from_exact_sol

  subroutine calc_diffrence_from_exact_sol_lm( ndim, X, X_exact, is, ie, myproc, icomm)
    integer, intent(in) :: ndim
    real*8, intent(in) :: X(is:ie), X_exact(is:ie)
    !
    integer is, ie, myproc, icomm, ierr
    !
    integer i
    real*8 znmx, zdmx, zerr
    real*8 znmx_recv, zdmx_recv
    !
    znmx = -1d0;    zdmx = -1d0
    do i=is,ie
       znmx = max( znmx, abs(X(i)-X_exact(i)) )
       zdmx = max( zdmx, abs(X_exact(i)) )
    enddo
    call MPI_Reduce( znmx, znmx_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, icomm, ierr )
    call MPI_Reduce( zdmx, zdmx_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, icomm, ierr )
    if (myproc == 0 ) then
       zerr = znmx_recv / zdmx_recv

       write(*,*)
       write(*,*) '    max(|numerical-exact|)/max(|exact|) =',zerr
       write(*,*) '    max(|numerical-exact|)              =',znmx_recv
       write(*,*) '                           max(|exact|) =',zdmx_recv
       write(*,*)
    end if
  end subroutine calc_diffrence_from_exact_sol_lm
  
  subroutine write_residuals( niter )
    integer, intent(in) :: niter
!
    integer i
!
    write(*,*)
    write(*,*) '--- Convergence history ---'
    DO I=1, niter
       write(*,'(1x,I4,1PE14.5)') I, his_resid(i)
    End do

  end subroutine write_residuals

end program main

