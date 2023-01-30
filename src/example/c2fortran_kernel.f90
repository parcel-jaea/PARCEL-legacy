subroutine parcel_init_c(myrank)
  implicit none
  integer myrank
  call PARCEL_Init(myrank)
end subroutine parcel_init_c

subroutine parcel_free_c
  implicit none
  call PARCEL_Free()
end subroutine parcel_free_c

subroutine get_namelist_input( &
     ityp_eq ,&
     gnx ,&
     gny ,&
     gnz ,&
     max_niter ,&
     rtolmax ,&
     abstolmax ,&
     ipreflag ,&
     Neuman_itrmax, &
     independ_nvec ,&
     addl ,&
     iflagas ,&
     ilu_method ,&
     iovlflag ,&
     ityp_solver ,&
     precon_thblock ,&
     restart ,&
     gmres_gsflag ,&
     cbcg_kstep ,&
     cbcg_eigenflag ,&
     power_method_itrmax ,&
     cagmres_sstep ,&
     cagmres_tstep ,&
     cagmres_basis ,&
     cagmres_qrflag ,&
     caarnoldi_sstep ,&
     caarnoldi_tstep ,&
     caarnoldi_basis ,&
     caarnoldi_qrflag &
     )

  use PARCEL_integer
  use problem
  implicit none

  integer ityp_eq
  integer(kind=inttype) gnx,gny,gnz
  integer max_niter
  integer ipreflag
  integer Neuman_itrmax
  integer addl
  integer iflagas
  real*8  rtolmax
  real*8  abstolmax
  integer ityp_solver

  integer ilu_method
  integer iovlflag
  integer gmres_m
  integer gmres_gsflag
  integer cagmres_sstep
  integer cagmres_tstep
  integer cagmres_basis
  integer cagmres_qrflag
  integer cbcg_kstep
  integer cbcg_eigenflag
  integer power_method_itrmax
  integer caarnoldi_sstep
  integer caarnoldi_tstep
  integer caarnoldi_basis
  integer caarnoldi_qrflag

  integer restart
  integer precon_thblock
  integer independ_nvec

  namelist/input/ &
       ityp_eq, &
       gnx, &
       gny, &
       gnz, &
       max_niter, &
       rtolmax, &
       abstolmax, &
       ipreflag, &
       Neuman_itrmax, &
       independ_nvec, &
       addl, &
       iflagas, &
       ilu_method, &
       iovlflag, &
       ityp_solver, &
       precon_thblock, &
       restart, &
       gmres_gsflag, &
       cbcg_kstep, &
       cbcg_eigenflag, &
       power_method_itrmax, &
       caarnoldi_sstep, &
       caarnoldi_tstep, &
       caarnoldi_basis, &
       caarnoldi_qrflag, &
       cagmres_sstep, &
       cagmres_tstep, &
       cagmres_basis, &
       cagmres_qrflag

  open(67,file='input_namelist')
  read(67,nml=input)
  close(67)

end subroutine get_namelist_input

subroutine get_namelist_input_dd( &
     PRECISION_A,  &
     PRECISION_b,  &
     PRECISION_x,  &
     PRECISION_PRECONDITION) 
  implicit none

  integer PRECISION_A
  integer PRECISION_b
  integer PRECISION_x
  integer PRECISION_PRECONDITION

  namelist/input_dd/ &
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION

  open(33,file='input_namelist')
  read(33,nml=input_dd)
  close(33)

end subroutine get_namelist_input_dd

subroutine get_namelist_input_ddm( &
     npe_x, &
     npe_y, &
     npe_z, &
     div_direc_th &
     )
  implicit none
  integer npe_x
  integer npe_y
  integer npe_z
  integer div_direc_th

  namelist/input_ddm/ &
       npe_x,&
       npe_y,&
       npe_z,&
       div_direc_th
  open(33,file='input_namelist')
  read(33,nml=input_ddm)
  close(33)

end subroutine get_namelist_input_ddm

subroutine get_init_poisson3d( &
     icomm, &
     gn, &
     gnx, &
     gny, &
     gnz, &
     istart, &
     iend, &
     nx, &
     ny, &
     nz, &
     nx0, &
     ny0, &
     nz0, &
     n, &
     npes, &
     myrank, &
     myranky, &
     myrankz, &
     npey, &
     npez &
     )
  use mpi
  use PARCEL_integer
  use problem
  implicit none
  integer icomm
  integer(kind=inttype) gn
  integer(kind=inttype) gnx,gny,gnz
  integer(kind=inttype) istart,iend

  integer(kind=inttype) nx,ny,nz
  integer(kind=inttype) nx0,ny0,nz0
  integer(kind=inttype) n
  integer npes,myrank
  integer myranky,myrankz
  integer npey,npez

  integer ierr

  if(gnz<npes)then
     npez = gnz
     npey = npes/npez
     if(( (npez*npey).ne.(npes) ) .or. ( mod(gny,npey).ne.0) )then
        if(myrank.eq.0)then
           write(*,*)"error :: not correct process"
           write(*,*)"   process = ",npes
           write(*,*)"   gny     = ",gny
           write(*,*)"   gnz     = ",gnz
           write(*,*)"   npey    = ",npey
           write(*,*)"   npez    = ",npez
        endif
        call mpi_barrier(icomm,ierr)
        call mpi_finalize(ierr)
        stop
     endif
  else
     npey = 1 
     npez = npes
     if(( (npez*npey).ne.(npes) ) .or. ( mod(gnz,npez).ne.0) )then
        if(myrank.eq.0)then
           write(*,*)"error :: not correct process"
           write(*,*)"   process = ",npes
           write(*,*)"   gny     = ",gny
           write(*,*)"   gnz     = ",gnz
           write(*,*)"   npey    = ",npey
           write(*,*)"   npez    = ",npez
        endif
        call mpi_barrier(icomm,ierr)
        call mpi_finalize(ierr)
        stop
     endif
  endif

  myranky = mod(myrank,npey)
  myrankz = myrank/npey
  nx      = gnx
  ny      = gny/npey
  nz      = gnz/npez
  nx0     = 1  
  ny0     = 1 + ny*myranky
  nz0     = 1 + nz*myrankz

  gn      = gnx*gny*gnz
  n       = nx * ny * nz
  istart  = 1+n*myrank
  iend    = istart + n -1

end subroutine get_init_poisson3d

subroutine make_poisson3d_crs_get_nnz_c( &
     nx,ny,nz, &
     nx0,ny0,nz0, &
     gnx,gny,gnz, &
     nnz, &
     istart &
     )

  use PARCEL_integer
  use problem
  implicit none
  integer(kind=inttype) nx,ny,nz
  integer(kind=inttype) nx0,ny0,nz0
  integer(kind=inttype) gnx,gny,gnz
  integer(kind=inttype) nnz
  integer(kind=inttype) istart

  call make_poisson3d_crs_get_nnz( &
       nx,ny,nz,   &
       nx0,ny0,nz0,&
       gnx,gny,gnz,&
       nnz, &
       istart)
end subroutine make_poisson3d_crs_get_nnz_c

subroutine make_poisson3d_crs_c( &
     nx,ny,nz,&
     nx0,ny0,nz0,&
     gnx,gny,gnz,&
     crsa,crscol,crsrow_ptr,n,nnz, &
     istart)
  use PARCEL_integer
  use problem
  implicit none
  integer(kind=inttype) nx,ny,nz
  integer(kind=inttype) nx0,ny0,nz0
  integer(kind=inttype) gnx,gny,gnz

  integer(kind=inttype) n
  integer(kind=inttype) nnz
  integer(kind=inttype) istart

  real*8  crsa(nnz)
  integer(kind=inttype) crscol(nnz)
  integer(kind=inttype) crsrow_ptr(n+1)

  call make_poisson3d_crs( &
       nx,ny,nz,&
       nx0,ny0,nz0,&
       gnx,gny,gnz,&       
       crsa,crscol,crsrow_ptr,n,nnz, &
       istart)

end subroutine make_poisson3d_crs_c

subroutine make_poisson3d_right_vector_c( &
     nx,ny,nz,   &
     nx0,ny0,nz0,&
     gnx,gny,gnz,&
     vecx0,vecb,&
     istart)
  use PARCEL_integer
  use problem
  implicit none
  integer(kind=inttype) nx,ny,nz
  integer(kind=inttype) nx0,ny0,nz0
  integer(kind=inttype) gnx,gny,gnz

  real*8 vecx0(nx*ny*nz)
  real*8 vecb(nx*ny*nz)
  integer(kind=inttype) istart

  call make_poisson3d_right_vector( &
       nx,ny,nz,   &
       nx0,ny0,nz0,&
       gnx,gny,gnz,&
       vecx0,vecb,&
       istart)

end subroutine make_poisson3d_right_vector_c

subroutine parcel_dcg_c(&
     icomm , &
     vecx, &
     vecb, &
     n,gn,nnz,istart, &
     crsa , &
     crsrow_ptr,&
     crscol,&
     ipreflag ,&
     ilu_method ,&
     addl ,&
     iflagas ,&
     niter ,&
     rtolmax ,&
     reshistory, &
     iovlflag ,&
     precon_thblock ,&
     Neuman_itr, &
     independ_nvec ,&
     nBlock, &
     iret &
     )
  use c_wrapper
  use PARCEL_wrapper
  use PARCEL_integer
  use iso_c_binding
  use cu_common_wrapper
  implicit none
  integer icomm
  integer(kind=inttype) n,gn,nnz,istart
#ifdef CUDA_SOLVER
  type(devicemem) vecx
  type(devicemem) vecb
  type(devicemem) crsa
  type(devicemem) crscol
  type(devicemem) crsrow_ptr
#else
  real*8 vecx(n)
  real*8 vecb(n)
  real*8  crsa(nnz)
  integer(kind=inttype) crscol(nnz)
  integer(kind=inttype) crsrow_ptr(n+1)
#endif

  integer niter
  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag
  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer iret

  ! integer(kind=inttype) crscol_host(nnz)
  ! integer(kind=inttype) crsrow_ptr_host(n+1)
  ! integer i

  ! call cudaMemcpy_d2h_parcel_int(crscol,crscol_host,nnz,0,0)    
  ! call cudaMemcpy_d2h_parcel_int(crsrow_ptr,crsrow_ptr_host,n+1,0,0)    
  ! do i = 1,10
  !    write(*,*)i,crscol_host(i),crsrow_ptr_host(i)
  ! enddo

  call  parcel_dcg(&
          icomm, &
#ifdef CUDA_SOLVER
          vecx%mem, &
          vecb%mem, &
          n,gn,nnz,istart, &
          crsa%mem , &
          crsrow_ptr%mem,&
          crscol%mem,&
#else
          vecx, &
          vecb, &
          n,gn,nnz,istart, &
          crsa , &
          crsrow_ptr,&
          crscol,&
#endif
          ipreflag,ilu_method,addl,iflagas, &
          niter,rtolmax,&
          reshistory,&
          iovlflag,  &
          precon_thblock, &
          Neuman_itr,independ_nvec, &
          nBlock, &
          iret &
          )

end subroutine parcel_dcg_c

subroutine parcel_dbicgstab_c(&
     icomm , &
     vecx, &
     vecb, &
     n,gn,nnz,istart, &
     crsa , &
     crsrow_ptr,&
     crscol,&
     ipreflag ,&
     ilu_method ,&
     addl ,&
     iflagas ,&
     niter ,&
     rtolmax ,&
     reshistory, &
     iovlflag ,&
     precon_thblock ,&
      Neuman_itr,&
      independ_nvec ,&
     nBlock, &
     iret &
     )
  use PARCEL_wrapper
  use iso_c_binding
  use PARCEL_integer

  implicit none
  integer icomm
  integer(kind=inttype) n,gn,nnz,istart
#ifdef CUDA_SOLVER
  type(c_ptr) vecx
  type(c_ptr) vecb
  type(c_ptr) crsa
  type(c_ptr) crscol
  type(c_ptr) crsrow_ptr
#else
  real*8 vecx(n)
  real*8 vecb(n)
  real*8  crsa(nnz)
  integer(kind=inttype) crscol(nnz)
  integer(kind=inttype) crsrow_ptr(n+1)
#endif

  integer niter
  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag
  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer iret

  call  parcel_dbicgstab(&
          icomm, &
#ifdef CUDA_SOLVER
          vecx, &
          vecb, &
          n,gn,nnz,istart, &
          crsa , &
          crsrow_ptr,&
          crscol,&
#else
          vecx, &
          vecb, &
          n,gn,nnz,istart, &
          crsa , &
          crsrow_ptr,&
          crscol,&
#endif
          ipreflag,ilu_method,addl,iflagas, &
          niter,rtolmax,&
          reshistory,&
          iovlflag,  &
          precon_thblock, &
          Neuman_itr,independ_nvec, &
          nBlock, &
          iret &
          )

end subroutine parcel_dbicgstab_c

subroutine parcel_dgmres_c(&
     icomm , &
     vecx, &
     vecb, &
     n,gn,nnz,istart, &
     crsa , &
     crsrow_ptr,&
     crscol,&
     ipreflag ,&
     ilu_method ,&
     addl ,&
     iflagas ,&
     niter ,&
     rtolmax ,&
     reshistory, &
     iovlflag ,&
     precon_thblock ,&
     Neuman_itr, &
     independ_nvec ,&
     nBlock, &
     restart ,&
     gmres_gsflag ,&
     iret &
     )
  use PARCEL_wrapper
  use iso_c_binding
  use PARCEL_integer

  implicit none
  integer icomm
  integer(kind=inttype) n,gn,nnz,istart
#ifdef CUDA_SOLVER
  type(c_ptr) vecx
  type(c_ptr) vecb
  type(c_ptr) crsa
  type(c_ptr) crscol
  type(c_ptr) crsrow_ptr
#else
  real*8 vecx(n)
  real*8 vecb(n)
  real*8  crsa(nnz)
  integer(kind=inttype) crscol(nnz)
  integer(kind=inttype) crsrow_ptr(n+1)
#endif

  integer niter
  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag

  integer restart
  integer gmres_gsflag

  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer iret

  call parcel_dgmres(&
          icomm, &
#ifdef CUDA_SOLVER
          vecx, &
          vecb, &
          n,gn,nnz,istart, &
          crsa , &
          crsrow_ptr,&
          crscol,&
#else
          vecx, &
          vecb, &
          n,gn,nnz,istart, &
          crsa , &
          crsrow_ptr,&
          crscol,&
#endif
          ipreflag,ilu_method,addl,iflagas, &
          niter,rtolmax,&
          reshistory,&
          iovlflag,  &
          precon_thblock, &
          Neuman_itr,independ_nvec, &
          nBlock, &
          restart ,&
          gmres_gsflag ,&
          iret &
          )

end subroutine parcel_dgmres_c

subroutine parcel_ddcg_c(&
     icomm , &
     vecx,vecx_lo, &
     vecb,vecb_lo, &
     n,gn,nnz,istart, &
     crsa ,crsa_lo, &
     crsrow_ptr,&
     crscol,&
     ipreflag ,&
     ilu_method ,&
     addl ,&
     iflagas ,&
     niter ,&
     rtolmax ,&
     reshistory, &
     iovlflag ,&
     precon_thblock ,&
     Neuman_itr,independ_nvec ,&
     nBlock, &
     PRECISION_A, &
     PRECISION_b, &
     PRECISION_x, &
     PRECISION_PRECONDITION, &
     iret &
     )
  use PARCEL_wrapper
  use iso_c_binding
  use PARCEL_integer

  implicit none
  integer icomm
  integer(kind=inttype) n,gn,nnz,istart
  real*8 vecx(n)
  real*8 vecb(n)
  real*8  crsa(nnz)

  real*8 vecx_lo(n)
  real*8 vecb_lo(n)
  real*8  crsa_lo(nnz)

  integer(kind=inttype) crscol(nnz)
  integer(kind=inttype) crsrow_ptr(n+1)

  integer niter
  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag
  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer PRECISION_A
  integer PRECISION_b
  integer PRECISION_x
  integer PRECISION_PRECONDITION

  integer iret

  call  parcel_ddcg(&
          icomm, &
          vecx,vecx_lo, &
          vecb,vecb_lo, &
          n,gn,nnz,istart, &
          crsa,crsa_lo, &
          crsrow_ptr,&
          crscol,&
          ipreflag,ilu_method,addl,iflagas, &
          niter,rtolmax,&
          reshistory,&
          iovlflag,  &
          precon_thblock, &
          Neuman_itr,&
          independ_nvec, &
          nBlock, &
          PRECISION_A, &
          PRECISION_b, &
          PRECISION_x, &
          PRECISION_PRECONDITION, &
          iret &
          )

end subroutine parcel_ddcg_c

subroutine parcel_ddbicgstab_c(&
     icomm , &
     vecx,vecx_lo, &
     vecb,vecb_lo, &
     n,gn,nnz,istart, &
     crsa ,crsa_lo, &
     crsrow_ptr,&
     crscol,&
     ipreflag ,&
     ilu_method ,&
     addl ,&
     iflagas ,&
     niter ,&
     rtolmax ,&
     reshistory, &
     iovlflag ,&
     precon_thblock ,&
     Neuman_itr,&
     independ_nvec ,&
     nBlock, &
     PRECISION_A, &
     PRECISION_b, &
     PRECISION_x, &
     PRECISION_PRECONDITION, &
     iret &
     )
  use PARCEL_wrapper
  use iso_c_binding
  use PARCEL_integer

  implicit none
  integer icomm
  integer(kind=inttype) n,gn,nnz,istart
  real*8 vecx(n)
  real*8 vecb(n)
  real*8  crsa(nnz)

  real*8 vecx_lo(n)
  real*8 vecb_lo(n)
  real*8  crsa_lo(nnz)

  integer(kind=inttype) crscol(nnz)
  integer(kind=inttype) crsrow_ptr(n+1)

  integer niter
  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag
  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer PRECISION_A
  integer PRECISION_b
  integer PRECISION_x
  integer PRECISION_PRECONDITION

  integer iret

  call  parcel_ddbicgstab(&
          icomm, &
          vecx,vecx_lo, &
          vecb,vecb_lo, &
          n,gn,nnz,istart, &
          crsa,crsa_lo, &
          crsrow_ptr,&
          crscol,&
          ipreflag,ilu_method,addl,iflagas, &
          niter,rtolmax,&
          reshistory,&
          iovlflag,  &
          precon_thblock, &
          Neuman_itr,&
          independ_nvec, &
          nBlock, &
          PRECISION_A, &
          PRECISION_b, &
          PRECISION_x, &
          PRECISION_PRECONDITION, &
          iret &
          )

end subroutine parcel_ddbicgstab_c

subroutine parcel_ddgmres_c(&
     icomm , &
     vecx,vecx_lo, &
     vecb,vecb_lo, &
     n,gn,nnz,istart, &
     crsa ,crsa_lo, &
     crsrow_ptr,&
     crscol,&
     ipreflag ,&
     ilu_method ,&
     addl ,&
     iflagas ,&
     niter ,&
     rtolmax ,&
     reshistory, &
     iovlflag ,&
     precon_thblock ,&
     Neuman_itr,&
     independ_nvec ,&
     nBlock, &
     restart ,&
     gmres_gsflag ,&
     PRECISION_A, &
     PRECISION_b, &
     PRECISION_x, &
     PRECISION_PRECONDITION, &
     iret &
     )
  use PARCEL_wrapper
  use iso_c_binding
  use PARCEL_integer

  implicit none
  integer icomm
  integer(kind=inttype) n,gn,nnz,istart
  real*8 vecx(n)
  real*8 vecb(n)
  real*8  crsa(nnz)

  real*8 vecx_lo(n)
  real*8 vecb_lo(n)
  real*8  crsa_lo(nnz)

  integer(kind=inttype) crscol(nnz)
  integer(kind=inttype) crsrow_ptr(n+1)

  integer niter
  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag

  integer restart
  integer gmres_gsflag

  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer PRECISION_A
  integer PRECISION_b
  integer PRECISION_x
  integer PRECISION_PRECONDITION

  integer iret

  call parcel_ddgmres(&
          icomm, &
          vecx,vecx_lo, &
          vecb,vecb_lo, &
          n,gn,nnz,istart, &
          crsa,crsa_lo, &
          crsrow_ptr,&
          crscol,&
          ipreflag,ilu_method,addl,iflagas, &
          niter,rtolmax,&
          reshistory,&
          iovlflag,  &
          precon_thblock, &
          Neuman_itr,&
          independ_nvec, &
          nBlock, &
          restart ,&
          gmres_gsflag ,&
          PRECISION_A, &
          PRECISION_b, &
          PRECISION_x, &
          PRECISION_PRECONDITION, &
          iret &
          )

end subroutine parcel_ddgmres_c


subroutine parcel_dcagmres_c(&
     icomm , &
     vecx, &
     vecb, &
     n,gn,nnz,istart, &
     crsa , &
     crsrow_ptr,&
     crscol,&
     ipreflag ,&
     ilu_method ,&
     addl ,&
     iflagas ,&
     niter ,&
     rtolmax ,&
     reshistory, &
     iovlflag ,&
     precon_thblock ,&
          Neuman_itr,&
     independ_nvec ,&
     nBlock, &
     cagmres_sstep ,&
     cagmres_tstep ,&
     cagmres_basis ,&
     cagmres_qrflag ,&
     iret &
     )
  use PARCEL_wrapper
  use iso_c_binding
  use PARCEL_integer

  implicit none
  integer icomm
  integer(kind=inttype) n,gn,nnz,istart
#ifdef CUDA_SOLVER
  type(c_ptr) vecx
  type(c_ptr) vecb
  type(c_ptr) crsa
  type(c_ptr) crscol
  type(c_ptr) crsrow_ptr
#else
  real*8 vecx(n)
  real*8 vecb(n)
  real*8  crsa(nnz)
  integer(kind=inttype) crscol(nnz)
  integer(kind=inttype) crsrow_ptr(n+1)
#endif

  integer niter
  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag

  integer cagmres_sstep
  integer cagmres_tstep
  integer cagmres_basis
  integer cagmres_qrflag

  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer iret

  call parcel_dcagmres(&
          icomm, &
#ifdef CUDA_SOLVER
          vecx, &
          vecb, &
          n,gn,nnz,istart, &
          crsa , &
          crsrow_ptr,&
          crscol,&

#else
          vecx, &
          vecb, &
          n,gn,nnz,istart, &
          crsa , &
          crsrow_ptr,&
          crscol,&
#endif
          ipreflag,ilu_method,addl,iflagas, &
          niter,rtolmax,&
          reshistory,&
          iovlflag,  &
          precon_thblock, &
          Neuman_itr,&
          independ_nvec, &
          nBlock, &
          cagmres_sstep ,&
          cagmres_tstep ,&
          cagmres_basis ,&
          cagmres_qrflag ,&
          iret &
          )

end subroutine parcel_dcagmres_c

subroutine parcel_dcbcg_c(&
     icomm , &
     vecx, &
     vecb, &
     n,gn,nnz,istart, &
     crsa , &
     crsrow_ptr,&
     crscol,&
     ipreflag ,&
     ilu_method ,&
     addl ,&
     iflagas ,&
     niter ,&
     rtolmax ,&
     reshistory, &
     iovlflag ,&
     precon_thblock ,&
          Neuman_itr,&
     independ_nvec ,&
     nBlock, &
     cbcg_kstep ,&
     cbcg_eigenflag ,&
     power_method_itrmax ,&
     caarnoldi_sstep ,&
     caarnoldi_tstep ,&
     caarnoldi_basis ,&
     caarnoldi_qrflag, &
     iret &
     )
  use PARCEL_wrapper
  use iso_c_binding
  use PARCEL_integer

  implicit none
  integer icomm
  integer(kind=inttype) n,gn,nnz,istart
#ifdef CUDA_SOLVER
  type(c_ptr) vecx
  type(c_ptr) vecb
  type(c_ptr) crsa
  type(c_ptr) crscol
  type(c_ptr) crsrow_ptr
#else
  real*8 vecx(n)
  real*8 vecb(n)
  real*8  crsa(nnz)
  integer(kind=inttype) crscol(nnz)
  integer(kind=inttype) crsrow_ptr(n+1)
#endif

  integer niter
  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag

  integer cbcg_kstep
  integer cbcg_eigenflag
  integer power_method_itrmax
  integer caarnoldi_sstep
  integer caarnoldi_tstep
  integer caarnoldi_basis
  integer caarnoldi_qrflag

  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer iret

  call parcel_dcbcg(&
          icomm, &
#ifdef CUDA_SOLVER
          vecx, &
          vecb, &
          n,gn,nnz,istart, &
          crsa , &
          crsrow_ptr,&
          crscol,&
#else
          vecx, &
          vecb, &
          n,gn,nnz,istart, &
          crsa , &
          crsrow_ptr,&
          crscol,&
#endif
          ipreflag,ilu_method,addl,iflagas, &
          niter,rtolmax,&
          reshistory,&
          iovlflag,  &
          precon_thblock, &
          Neuman_itr,&
          independ_nvec, &
          nBlock, &
          cbcg_kstep ,&
          cbcg_eigenflag ,&
          power_method_itrmax ,&
          caarnoldi_sstep ,&
          caarnoldi_tstep ,&
          caarnoldi_basis ,&
          caarnoldi_qrflag, &
          iret &
          )

end subroutine parcel_dcbcg_c

subroutine parcel_dcaarnoldi_c(&
     icomm , &
     vecx, &
     vecb, &
     n,gn,nnz,istart, &
     crsa , &
     crsrow_ptr,&
     crscol,&
     ipreflag ,&
     ilu_method ,&
     addl ,&
     iflagas ,&
     iovlflag ,&
     precon_thblock ,&
          Neuman_itr,&
     independ_nvec ,&
     nBlock, &
     caarnoldi_sstep ,&
     caarnoldi_tstep ,&
     caarnoldi_basis ,&
     caarnoldi_qrflag, &
     EigenValueReal,&
     EigenValueImag,&
     EigenVector   ,&
     EigenErrorNorm,&
     EmaxID,EminID, &
     iret &
     )
  use PARCEL_wrapper
  use iso_c_binding
  use PARCEL_integer

! parcel_dpower_dia,parcel_dcaarnoldi_dia

  implicit none
  integer icomm
  integer(kind=inttype) n,gn,nnz,istart
#ifdef CUDA_SOLVER
  type(c_ptr) vecx
  type(c_ptr) vecb
  type(c_ptr) crsa
  type(c_ptr) crscol
  type(c_ptr) crsrow_ptr
#else
  real*8 vecx(n)
  real*8 vecb(n)
  real*8  crsa(nnz)
  integer(kind=inttype) crscol(nnz)
  integer(kind=inttype) crsrow_ptr(n+1)
#endif

  integer niter
  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag

  integer caarnoldi_sstep
  integer caarnoldi_tstep
  integer caarnoldi_basis
  integer caarnoldi_qrflag

  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(caarnoldi_sstep*caarnoldi_tstep)

  real*8 EigenValueReal(caarnoldi_sstep * caarnoldi_tstep)
  real*8 EigenValueImag(caarnoldi_sstep * caarnoldi_tstep)
  real*8 EigenVector(caarnoldi_sstep * caarnoldi_tstep*n)
  real*8 EigenErrorNorm(caarnoldi_sstep * caarnoldi_tstep)
  integer EmaxID,EminID

  integer iret
  integer solver

  call  parcel_dcaarnoldi(&
       icomm, &
       vecx, &
       vecb, &
       n,gn,nnz,istart, &
       crsa , &
       crsrow_ptr,&
       crscol,&
       ipreflag,ilu_method,addl,iflagas, &
       iovlflag,  &
       precon_thblock, &
          Neuman_itr,&
       independ_nvec, &
       nBlock, &
       caarnoldi_sstep, &
       caarnoldi_tstep, &
       caarnoldi_basis, &
       caarnoldi_qrflag, &
       EigenValueReal,&
       EigenValueImag,&
       EigenVector   ,&
       EigenErrorNorm,&
       EmaxID,EminID, &
       iret &
       )

end subroutine parcel_dcaarnoldi_c

subroutine parcel_dcg_dia_c(&
     icomm , &
     vecx, &
     vecb, &
     n,gn,istart, &
     diaA ,  &
     offset, &
     nnd, &
     ipreflag ,&
     ilu_method ,&
     addl ,&
     iflagas ,&
     niter ,&
     rtolmax ,&
     reshistory, &
     iovlflag ,&
     precon_thblock ,&
          Neuman_itr,&
     independ_nvec ,&
     nBlock, &
     iret &
     )
  use PARCEL_wrapper
  use iso_c_binding
  use PARCEL_integer

  implicit none
  integer icomm
  integer(kind=inttype) n,gn,istart
#ifdef CUDA_SOLVER
  type(c_ptr) vecx
  type(c_ptr) vecb
  integer(kind=inttype) nnd
  type(c_ptr) diaA
  type(c_ptr) offset
#else
  real*8 vecx(n)
  real*8 vecb(n)
  integer(kind=inttype) nnd
  real*8  diaA(nnd*n)
  integer(kind=inttype)offset(nnd)
#endif
  integer niter
  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag
  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer iret

  call  parcel_dcg_dia(&
          icomm, &
#ifdef CUDA_SOLVER
          vecx, &
          vecb, &
          n,gn,istart, &
          diaA,offset,nnd,  &
#else
          vecx, &
          vecb, &
          n,gn,istart, &
          diaA,offset,nnd,  &
#endif
          ipreflag,ilu_method,addl,iflagas, &
          niter,rtolmax,&
          reshistory,&
          iovlflag,  &
          precon_thblock, &
          Neuman_itr,&
          independ_nvec, &
          nBlock, &
          iret &
          )

end subroutine parcel_dcg_dia_c

subroutine parcel_dbicgstab_dia_c(&
     icomm , &
     vecx, &
     vecb, &
     n,gn,istart, &
     diaA ,  &
     offset, &
     nnd, &
     ipreflag ,&
     ilu_method ,&
     addl ,&
     iflagas ,&
     niter ,&
     rtolmax ,&
     reshistory, &
     iovlflag ,&
     precon_thblock ,&
          Neuman_itr,&
     independ_nvec ,&
     nBlock, &
     iret &
     )
  use PARCEL_wrapper
  use iso_c_binding
  use PARCEL_integer

  implicit none
  integer icomm
  integer(kind=inttype) n,gn,istart
#ifdef CUDA_SOLVER
  type(c_ptr) vecx
  type(c_ptr) vecb
  integer(kind=inttype) nnd
  type(c_ptr) diaA
  type(c_ptr) offset
#else
  real*8 vecx(n)
  real*8 vecb(n)
  integer(kind=inttype) nnd
  real*8  diaA(nnd*n)
  integer(kind=inttype)offset(nnd)
#endif
  integer niter
  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag
  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer iret

  call  parcel_dbicgstab_dia(&
          icomm, &
#ifdef CUDA_SOLVER
          vecx, &
          vecb, &
          n,gn,istart, &
          diaA,offset,nnd,  &
#else
          vecx, &
          vecb, &
          n,gn,istart, &
          diaA,offset,nnd,  &
#endif
          ipreflag,ilu_method,addl,iflagas, &
          niter,rtolmax,&
          reshistory,&
          iovlflag,  &
          precon_thblock, &
          Neuman_itr,&
          independ_nvec, &
          nBlock, &
          iret &
          )

end subroutine parcel_dbicgstab_dia_c

subroutine parcel_dgmres_dia_c(&
     icomm , &
     vecx, &
     vecb, &
     n,gn,istart, &
     diaA ,  &
     offset, &
     nnd, &
     ipreflag ,&
     ilu_method ,&
     addl ,&
     iflagas ,&
     niter ,&
     rtolmax ,&
     reshistory, &
     iovlflag ,&
     precon_thblock ,&
          Neuman_itr,&
     independ_nvec ,&
     nBlock, &
     restart ,&
     gmres_gsflag ,&
     iret &
     )
  use PARCEL_wrapper
  use iso_c_binding
  use PARCEL_integer

  implicit none
  integer icomm
  integer(kind=inttype) n,gn,istart
#ifdef CUDA_SOLVER
  type(c_ptr) vecx
  type(c_ptr) vecb
  integer(kind=inttype) nnd
  type(c_ptr) diaA
  type(c_ptr) offset
#else
  real*8 vecx(n)
  real*8 vecb(n)
  integer(kind=inttype) nnd
  real*8  diaA(nnd*n)
  integer(kind=inttype)offset(nnd)
#endif
  integer niter
  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag

  integer restart
  integer gmres_gsflag

  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer iret

  call parcel_dgmres_dia(&
          icomm, &
#ifdef CUDA_SOLVER
          vecx, &
          vecb, &
          n,gn,istart, &
          diaA,offset,nnd,  &
#else
          vecx, &
          vecb, &
          n,gn,istart, &
          diaA,offset,nnd,  &
#endif
          ipreflag,ilu_method,addl,iflagas, &
          niter,rtolmax,&
          reshistory,&
          iovlflag,  &
          precon_thblock, &
          Neuman_itr,&
          independ_nvec, &
          nBlock, &
          restart ,&
          gmres_gsflag ,&
          iret &
          )

end subroutine parcel_dgmres_dia_c

subroutine parcel_ddcg_dia_c(&
     icomm , &
     vecx,vecx_lo, &
     vecb,vecb_lo, &
     n,gn,istart, &
     diaA ,diaA_lo ,  &
     offset, &
     nnd, &
     ipreflag ,&
     ilu_method ,&
     addl ,&
     iflagas ,&
     niter ,&
     rtolmax ,&
     reshistory, &
     iovlflag ,&
     precon_thblock ,&
          Neuman_itr,&
     independ_nvec ,&
     nBlock, &
     PRECISION_A, &
     PRECISION_b, &
     PRECISION_x, &
     PRECISION_PRECONDITION, &
     iret &
     )
  use PARCEL_wrapper
  use iso_c_binding
  use PARCEL_integer

  implicit none
  integer icomm
  integer(kind=inttype) n,gn,istart
#ifdef CUDA_SOLVER
  type(c_ptr) vecx
  type(c_ptr) vecb
  integer(kind=inttype) nnd
  type(c_ptr) diaA
  type(c_ptr) offset
#else
  real*8 vecx(n)
  real*8 vecb(n)
  integer(kind=inttype) nnd
  real*8  diaA(nnd*n)
  integer(kind=inttype)offset(nnd)
#endif

  real*8 vecx_lo(n)
  real*8 vecb_lo(n)
  real*8  diaA_lo(nnd*n)

  integer niter
  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag
  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer PRECISION_A
  integer PRECISION_b
  integer PRECISION_x
  integer PRECISION_PRECONDITION
  integer iret

  call  parcel_ddcg_dia(&
          icomm, &
          vecx,vecx_lo, &
          vecb,vecb_lo, &
          n,gn,istart, &
          diaA,diaA_lo,offset,nnd,  &
          ipreflag,ilu_method,addl,iflagas, &
          niter,rtolmax,&
          reshistory,&
          iovlflag,  &
          precon_thblock, &
	  Neuman_itr, &
          independ_nvec, &
          nBlock, &
          PRECISION_A, &
          PRECISION_b, &
          PRECISION_x, &
          PRECISION_PRECONDITION, &
          iret &
          )

end subroutine parcel_ddcg_dia_c

subroutine parcel_ddbicgstab_dia_c(&
     icomm , &
     vecx,vecx_lo, &
     vecb,vecb_lo, &
     n,gn,istart, &
     diaA ,diaA_lo ,  &
     offset, &
     nnd, &
     ipreflag ,&
     ilu_method ,&
     addl ,&
     iflagas ,&
     niter ,&
     rtolmax ,&
     reshistory, &
     iovlflag ,&
     precon_thblock ,&
          Neuman_itr,&
     independ_nvec ,&
     nBlock, &
     PRECISION_A, &
     PRECISION_b, &
     PRECISION_x, &
     PRECISION_PRECONDITION, &
     iret &
     )
  use PARCEL_wrapper
  use iso_c_binding
  use PARCEL_integer

  implicit none
  integer icomm
  integer(kind=inttype) n,gn,istart
#ifdef CUDA_SOLVER
  type(c_ptr) vecx
  type(c_ptr) vecb
  integer(kind=inttype) nnd
  type(c_ptr) diaA
  type(c_ptr) offset
#else
  real*8 vecx(n)
  real*8 vecb(n)
  integer(kind=inttype) nnd
  real*8  diaA(nnd*n)
  integer(kind=inttype)offset(nnd)
#endif

  real*8 vecx_lo(n)
  real*8 vecb_lo(n)
  real*8  diaA_lo(nnd*n)

  integer niter
  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag
  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer PRECISION_A
  integer PRECISION_b
  integer PRECISION_x
  integer PRECISION_PRECONDITION
  integer iret

  call  parcel_ddbicgstab_dia(&
          icomm, &
          vecx,vecx_lo, &
          vecb,vecb_lo, &
          n,gn,istart, &
          diaA,diaA_lo,offset,nnd,  &
          ipreflag,ilu_method,addl,iflagas, &
          niter,rtolmax,&
          reshistory,&
          iovlflag,  &
          precon_thblock, &
	  Neuman_itr, &
          independ_nvec, &
          nBlock, &
          PRECISION_A, &
          PRECISION_b, &
          PRECISION_x, &
          PRECISION_PRECONDITION, &
          iret &
          )

end subroutine parcel_ddbicgstab_dia_c

subroutine parcel_ddgmres_dia_c(&
     icomm , &
     vecx,vecx_lo, &
     vecb,vecb_lo, &
     n,gn,istart, &
     diaA ,diaA_lo ,  &
     offset, &
     nnd, &
     ipreflag ,&
     ilu_method ,&
     addl ,&
     iflagas ,&
     niter ,&
     rtolmax ,&
     reshistory, &
     iovlflag ,&
     precon_thblock ,&
     independ_nvec ,&
     nBlock, &
     restart ,&
     gmres_gsflag ,&
     PRECISION_A, &
     PRECISION_b, &
     PRECISION_x, &
     PRECISION_PRECONDITION, &
     iret &
     )
  use PARCEL_wrapper
  use iso_c_binding
  use PARCEL_integer

  implicit none
  integer icomm
  integer(kind=inttype) n,gn,istart
  real*8 vecx(n)
  real*8 vecb(n)
  integer(kind=inttype) nnd
  real*8  diaA(nnd*n)

  real*8 vecx_lo(n)
  real*8 vecb_lo(n)
  real*8  diaA_lo(nnd*n)

  integer(kind=inttype)offset(nnd)
  integer niter
  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag

  integer restart
  integer gmres_gsflag

  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer PRECISION_A
  integer PRECISION_b
  integer PRECISION_x
  integer PRECISION_PRECONDITION
  integer iret

  call parcel_ddgmres_dia(&
          icomm, &
          vecx,vecx_lo, &
          vecb,vecb_lo, &
          n,gn,istart, &
          diaA,diaA_lo,offset,nnd,  &
          ipreflag,ilu_method,addl,iflagas, &
          niter,rtolmax,&
          reshistory,&
          iovlflag,  &
          precon_thblock, &
	  Neuman_itr, &
          independ_nvec, &
          nBlock, &
          restart ,&
          gmres_gsflag ,&
          PRECISION_A, &
          PRECISION_b, &
          PRECISION_x, &
          PRECISION_PRECONDITION, &
          iret &
          )

end subroutine parcel_ddgmres_dia_c

subroutine parcel_dcagmres_dia_c(&
     icomm , &
     vecx, &
     vecb, &
     n,gn,istart, &
     diaA ,  &
     offset, &
     nnd, &
     ipreflag ,&
     ilu_method ,&
     addl ,&
     iflagas ,&
     niter ,&
     rtolmax ,&
     reshistory, &
     iovlflag ,&
     precon_thblock ,&
     Neuman_itr, &
     independ_nvec ,&
     nBlock, &
     cagmres_sstep ,&
     cagmres_tstep ,&
     cagmres_basis ,&
     cagmres_qrflag ,&
     iret &
     )
  use PARCEL_wrapper
  use iso_c_binding
  use PARCEL_integer

  implicit none
  integer icomm
  integer(kind=inttype) n,gn,istart
#ifdef CUDA_SOLVER
  type(c_ptr) vecx
  type(c_ptr) vecb
  integer(kind=inttype) nnd
  type(c_ptr) diaA
  type(c_ptr) offset
#else
  real*8 vecx(n)
  real*8 vecb(n)
  integer(kind=inttype) nnd
  real*8  diaA(nnd*n)
  integer(kind=inttype)offset(nnd)
#endif
  integer niter
  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag

  integer cagmres_sstep
  integer cagmres_tstep
  integer cagmres_basis
  integer cagmres_qrflag

  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer iret

  call parcel_dcagmres_dia(&
          icomm, &
#ifdef CUDA_SOLVER
          vecx, &
          vecb, &
          n,gn,istart, &
          diaA,offset,nnd,  &
#else
          vecx, &
          vecb, &
          n,gn,istart, &
          diaA,offset,nnd,  &
#endif
          ipreflag,ilu_method,addl,iflagas, &
          niter,rtolmax,&
          reshistory,&
          iovlflag,  &
          precon_thblock, &
	  Neuman_itr, &
          independ_nvec, &
          nBlock, &
          cagmres_sstep ,&
          cagmres_tstep ,&
          cagmres_basis ,&
          cagmres_qrflag ,&
          iret &
          )

end subroutine parcel_dcagmres_dia_c

subroutine parcel_dcbcg_dia_c(&
     icomm , &
     vecx, &
     vecb, &
     n,gn,istart, &
     diaA ,  &
     offset, &
     nnd, &
     ipreflag ,&
     ilu_method ,&
     addl ,&
     iflagas ,&
     niter ,&
     rtolmax ,&
     reshistory, &
     iovlflag ,&
     precon_thblock ,&
     Neuman_itr, &
     independ_nvec ,&
     nBlock, &
     cbcg_kstep ,&
     cbcg_eigenflag ,&
     power_method_itrmax ,&
     caarnoldi_sstep ,&
     caarnoldi_tstep ,&
     caarnoldi_basis ,&
     caarnoldi_qrflag, &
     iret &
     )
  use PARCEL_wrapper
  use iso_c_binding
  use PARCEL_integer

  implicit none
  integer icomm
  integer(kind=inttype) n,gn,istart
#ifdef CUDA_SOLVER
  type(c_ptr) vecx
  type(c_ptr) vecb
  integer(kind=inttype) nnd
  type(c_ptr) diaA
  type(c_ptr) offset
#else
  real*8 vecx(n)
  real*8 vecb(n)
  integer(kind=inttype) nnd
  real*8  diaA(nnd*n)
  integer(kind=inttype)offset(nnd)
#endif
  integer niter
  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag

  integer cbcg_kstep
  integer cbcg_eigenflag
  integer power_method_itrmax
  integer caarnoldi_sstep
  integer caarnoldi_tstep
  integer caarnoldi_basis
  integer caarnoldi_qrflag

  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer iret

  call parcel_dcbcg_dia(&
          icomm, &
#ifdef CUDA_SOLVER
          vecx, &
          vecb, &
          n,gn,istart, &
          diaA,offset,nnd,  &
#else
          vecx, &
          vecb, &
          n,gn,istart, &
          diaA,offset,nnd,  &
#endif
          ipreflag,ilu_method,addl,iflagas, &
          niter,rtolmax,&
          reshistory,&
          iovlflag,  &
          precon_thblock, &
	  Neuman_itr, &
          independ_nvec, &
          nBlock, &
          cbcg_kstep ,&
          cbcg_eigenflag ,&
          power_method_itrmax ,&
          caarnoldi_sstep ,&
          caarnoldi_tstep ,&
          caarnoldi_basis ,&
          caarnoldi_qrflag, &
          iret &
          )

end subroutine parcel_dcbcg_dia_c

subroutine parcel_dcaarnoldi_dia_c(&
     icomm , &
     vecx, &
     vecb, &
     n,gn,istart, &
     diaA ,  &
     offset, &
     nnd, &
     ipreflag ,&
     ilu_method ,&
     addl ,&
     iflagas ,&
     iovlflag ,&
     precon_thblock ,&
     independ_nvec ,&
     nBlock, &
     caarnoldi_sstep ,&
     caarnoldi_tstep ,&
     caarnoldi_basis ,&
     caarnoldi_qrflag, &
     EigenValueReal,&
     EigenValueImag,&
     EigenVector   ,&
     EigenErrorNorm,&
     EmaxID,EminID, &
     iret &
     )
  use PARCEL_wrapper
  use iso_c_binding
  use PARCEL_integer

! parcel_dpower_dia,parcel_dcaarnoldi_dia

  implicit none
  integer icomm
  integer(kind=inttype) n,gn,istart
#ifdef CUDA_SOLVER
  type(c_ptr) vecx
  type(c_ptr) vecb
  integer(kind=inttype) nnd
  type(c_ptr) diaA
  type(c_ptr) offset
#else
  real*8 vecx(n)
  real*8 vecb(n)
  integer(kind=inttype) nnd
  real*8  diaA(nnd*n)
  integer(kind=inttype)offset(nnd)
#endif
  integer niter
  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag

  integer caarnoldi_sstep
  integer caarnoldi_tstep
  integer caarnoldi_basis
  integer caarnoldi_qrflag

  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(caarnoldi_sstep*caarnoldi_tstep)

  real*8 EigenValueReal(caarnoldi_sstep * caarnoldi_tstep)
  real*8 EigenValueImag(caarnoldi_sstep * caarnoldi_tstep)
  real*8 EigenVector(caarnoldi_sstep * caarnoldi_tstep*n)
  real*8 EigenErrorNorm(caarnoldi_sstep * caarnoldi_tstep)
  integer EmaxID,EminID

  integer iret
  integer solver

#ifdef CUDA_SOLVER
     write(*,*)"not support GPU CA-ARNOlDI"
#else
  call  parcel_dcaarnoldi_dia(&
       icomm, &
       vecx, &
       vecb, &
       n,gn,istart, &
       diaA ,  &
       offset, &
       nnd, &
       ipreflag,ilu_method,addl,iflagas, &
       iovlflag,  &
       precon_thblock, &
	  Neuman_itr, &
       independ_nvec, &
       nBlock, &
       caarnoldi_sstep, &
       caarnoldi_tstep, &
       caarnoldi_basis, &
       caarnoldi_qrflag, &
       EigenValueReal,&
       EigenValueImag,&
       EigenVector   ,&
       EigenErrorNorm,&
       EmaxID,EminID, &
       iret &
       )
#endif

end subroutine parcel_dcaarnoldi_dia_c


#if 1
subroutine parcel_solve_dd_c( &
     icomm , &
     n,gn,nnz,istart, &
     reshistory, &
     vecx, &
     vecb, &
     crsa , &
     crsrow_ptr,&
     crscol,&
     niter ,&
     rtolmax ,&
     abstolmax ,&
     ipreflag ,&
     independ_nvec ,&
     nBlock, &
     addl ,&
     iflagas ,&
     ilu_method ,&
     iovlflag ,&
     solver ,&
     precon_thblock ,&
     restart ,&
     gmres_gsflag ,&
     cbcg_kstep ,&
     cbcg_eigenflag ,&
     power_method_itrmax ,&
     cagmres_sstep ,&
     cagmres_tstep ,&
     cagmres_basis ,&
     cagmres_qrflag ,&
     caarnoldi_sstep ,&
     caarnoldi_tstep ,&
     caarnoldi_basis ,&
     caarnoldi_qrflag, &
     PRECISION_A, &
     PRECISION_b, &
     PRECISION_x, &
     PRECISION_PRECONDITION, &
     vecx_lo, &
     vecb_lo, &
     crsa_lo, &
     iret &
     )

  use PARCEL_wrapper
  use iso_c_binding
  use PARCEL_integer
  use problem
  use krylov
  use krylov_kernel

  implicit none
  integer icomm
  integer(kind=inttype) n,gn,nnz,istart
  real*8 vecx(n)
  real*8 vecb(n)
  real*8 crsa(nnz)
  real*8 vecx_lo(n)
  real*8 vecb_lo(n)
  real*8 crsa_lo(nnz)
  integer(kind=inttype) crscol(nnz)
  integer(kind=inttype) crsrow_ptr(n+1)

  integer niter
  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag
  integer gmres_m
  integer gmres_gsflag
  integer cagmres_sstep
  integer cagmres_tstep
  integer cagmres_basis
  integer cagmres_qrflag
  integer cbcg_kstep
  integer cbcg_eigenflag
  integer power_method_itrmax
  integer caarnoldi_sstep
  integer caarnoldi_tstep
  integer caarnoldi_basis
  integer caarnoldi_qrflag
  integer restart
  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer PRECISION_A
  integer PRECISION_b
  integer PRECISION_x
  integer PRECISION_PRECONDITION
  integer iret
  integer solver

  select case( solver )
  case(1)
     call  parcel_ddcg(&
          icomm, &
          vecx,vecx_lo,&
          vecb,vecb_lo,  &
          n,gn,nnz,istart, &
          crsA , crsA_lo, &
          crsRow_ptr,&
          crsCol,&
          ipreflag,ILU_method,addL,iflagAS, &
          niter,rtolmax,&
          reshistory,&
          iovlflag,  &
          precon_thblock, &
	  Neuman_itr, &
          independ_nvec, &
          nBlock, &
          PRECISION_A, &
          PRECISION_b, &
          PRECISION_x, &
          PRECISION_PRECONDITION, &
          iret &
          )
  case(2)
     call  parcel_ddbicgstab(&
          icomm, &
          vecx,vecx_lo,&
          vecb,vecb_lo,  &
          n,gn,nnz,istart, &
          crsA , crsA_lo, &
          crsRow_ptr,&
          crsCol,&
          ipreflag,ILU_method,addL,iflagAS, &
          niter,rtolmax,&
          reshistory,&
          iovlflag,  &
          precon_thblock, &
	  Neuman_itr, &
          independ_nvec, &
          nBlock, &
          PRECISION_A, &
          PRECISION_b, &
          PRECISION_x, &
          PRECISION_PRECONDITION, &
          iret &
          )
  case(3)
     call  parcel_ddgmres(&
          icomm, &
          vecx,vecx_lo,&
          vecb,vecb_lo,  &
          n,gn,nnz,istart, &
          crsA , crsA_lo, &
          crsRow_ptr,&
          crsCol,&
          ipreflag,ILU_method,addL,iflagAS, &
          niter,rtolmax,&
          reshistory,&
          iovlflag,  &
          precon_thblock, &
	  Neuman_itr, &
          independ_nvec, &
          nBlock, &
          restart,GMRES_GSflag, &
          PRECISION_A, &
          PRECISION_b, &
          PRECISION_x, &
          PRECISION_PRECONDITION, &
          iret &
          )
  end select

end subroutine parcel_solve_dd_c
#endif

subroutine crs2dia_c( &
     n,nnz, &
     crsa , crsrow_ptr, crscol, &
     diaA , offset, &
     nnd, &
     istart,iend &
     )
  use PARCEL_integer
  use krylov
  use matrix_kernel, only: &
       get_nnd_crs_without_sort,crs2dia_without_sort, &
       free_diaMat,free_crsMat, &
       init_crsMat
  implicit none
  integer(kind=intType) istart,iend 
  integer(kind=inttype) n,nnz
  real*8  crsa(nnz)
  integer(kind=inttype) crscol(nnz)
  integer(kind=inttype) crsrow_ptr(n+1)
  integer(kind=intType) nnd
  real*8  diaA(nnd*n)
  integer(kind=intType) offset(nnd)
  integer(kind=intType) i

  type(CRS_Mat) :: crsMat 
  type(DIA_Mat) :: diaMat 

  call init_crsMat(crsMat,n,nnz)
  ! write(*,*)" n,nnz = ",n,nnz
  ! write(*,*)" crsA       = ",size(crsMat%crsA),size(crsa)
  ! write(*,*)" crsCol     = ",size(crsMat%crsCol),size(crscol)
  ! write(*,*)" crsRow_ptr = ",size(crsMat%crsRow_ptr),size(crsrow_ptr)
  do i = 1,n+1
     crsMat%crsRow_ptr(i) = crsrow_ptr(i)
  enddo

  do i = 1,nnz
     crsMat%crsA(i)   = crsa(i)
     crsMat%crsCol(i) = crscol(i)
  enddo

  call get_nnd_crs_without_sort(&
       diaMat,n,nnz, &
       crsMat%crsCol,&
       crsMat%crsRow_ptr,&
       istart,&
       iend,&
       nnd &
       )

  call crs2dia_without_sort(crsMat,diaMat,nnd,n,nnz,istart)
  do i = 1,nnd*n
     diaA(i) = diaMat%diaA(i)
  enddo
  do i = 1,nnd
     offset(i) = diaMat%offset(i)
  enddo
  call free_crsMat(crsMat)
  call free_diaMat(diaMat)

end subroutine crs2dia_c


subroutine parcel_solve_dia_dd_c(&
     icomm , &
     n,gn,istart, &
     reshistory, &
     vecx, &
     vecb, &
     diaA ,  &
     offset, &
     nnd, &
     niter ,&
     rtolmax ,&
     abstolmax ,&
     ipreflag ,&
     independ_nvec ,&
     nBlock, &
     addl ,&
     iflagas ,&
     ilu_method ,&
     iovlflag ,&
     solver ,&
     precon_thblock ,&
     restart ,&
     gmres_gsflag ,&
     cbcg_kstep ,&
     cbcg_eigenflag ,&
     power_method_itrmax ,&
     cagmres_sstep ,&
     cagmres_tstep ,&
     cagmres_basis ,&
     cagmres_qrflag ,&
     caarnoldi_sstep ,&
     caarnoldi_tstep ,&
     caarnoldi_basis ,&
     caarnoldi_qrflag, &
     PRECISION_A, &
     PRECISION_b, &
     PRECISION_x, &
     PRECISION_PRECONDITION, &
     vecx_lo, &
     vecb_lo, &
     diaa_lo, &
     iret &
     )
  use PARCEL_integer
  use krylov
  use krylov_kernel

  implicit none
  integer icomm
  integer(kind=inttype) n,gn,istart
  real*8 vecx(n)
  real*8 vecb(n)

  integer(kind=inttype) nnd
  real*8  diaA(nnd*n)
  integer(kind=inttype)offset(nnd)

  real*8 vecx_lo(n)
  real*8 vecb_lo(n)
  real*8  diaA_lo(nnd*n)

  integer niter
  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag
  integer gmres_m
  integer gmres_gsflag
  integer cagmres_sstep
  integer cagmres_tstep
  integer cagmres_basis
  integer cagmres_qrflag
  integer cbcg_kstep
  integer cbcg_eigenflag
  integer power_method_itrmax
  integer caarnoldi_sstep
  integer caarnoldi_tstep
  integer caarnoldi_basis
  integer caarnoldi_qrflag
  integer restart
  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer PRECISION_A
  integer PRECISION_b
  integer PRECISION_x
  integer PRECISION_PRECONDITION
  integer iret
  integer solver

  select case( solver )
  case(1)
     call  parcel_ddcg_dia(&
          icomm, &
          vecx,vecx_lo,&
          vecb,vecb_lo,  &
          n,gn,istart, &
          diaA , diaA_lo, &
          offset,nnd, &
          ipreflag,ILU_method,addL,iflagAS, &
          niter,rtolmax,&
          reshistory,&
          iovlflag,  &
          precon_thblock, &
	  Neuman_itr, &
          independ_nvec, &
          nBlock, &
          PRECISION_A, &
          PRECISION_b, &
          PRECISION_x, &
          PRECISION_PRECONDITION, &
          iret &
          )
  case(2)
     call  parcel_ddbicgstab_dia(&
          icomm, &
          vecx,vecx_lo,&
          vecb,vecb_lo,  &
          n,gn,istart, &
          diaA , diaA_lo, &
          offset,nnd, &
          ipreflag,ILU_method,addL,iflagAS, &
          niter,rtolmax,&
          reshistory,&
          iovlflag,  &
          precon_thblock, &
	  Neuman_itr, &
          independ_nvec, &
          nBlock, &
          PRECISION_A, &
          PRECISION_b, &
          PRECISION_x, &
          PRECISION_PRECONDITION, &
          iret &
          )
  case(3)
     call  parcel_ddgmres_dia(&
          icomm, &
          vecx,vecx_lo,&
          vecb,vecb_lo,  &
          n,gn,istart, &
          diaA , diaA_lo, &
          offset,nnd, &
          ipreflag,ILU_method,addL,iflagAS, &
          niter,rtolmax,&
          reshistory,&
          iovlflag,  &
          precon_thblock, &
	  Neuman_itr, &
          independ_nvec, &
          nBlock, &
          restart,GMRES_GSflag, &
          PRECISION_A, &
          PRECISION_b, &
          PRECISION_x, &
          PRECISION_PRECONDITION, &
          iret &
          )
  end select

end subroutine parcel_solve_dia_dd_c

subroutine check_exact_solution_c(n,icomm,vecx,vecx0)
  use PARCEL_integer
  use problem
  implicit none
  integer icomm
  integer(kind=inttype) n
  real*8 vecx(n)
  real*8 vecx0(n)

  call check_exact_solution(n,icomm,vecx,vecx0)

end subroutine check_exact_solution_c

subroutine parcel_dqr_c(&
     n,s,icomm,A,Q,R,iQRflag)

  use PARCEL_integer
  implicit none
  integer(kind=intType) n
  integer(kind=intType) s
  integer icomm
  real*8 A(n*s)
  real*8 Q(n*s)
  real*8 R(s*s)
  integer iQRflag

  call parcel_dqr(n,s,icomm,A,Q,R,iQRflag)

end subroutine parcel_dqr_c

subroutine check_qr_c(n,s,icomm,A,Q,R,A0)
  use QR_kernel
  use mpi
  use PARCEL_integer
  implicit none

  integer(kind=intType) n
  integer(kind=intType) s
  integer icomm
  real*8 A(n*s)
  real*8 A0(n*s)
  real*8 Q(n*s)
  real*8 R(s*s)

  call check_QtQ(n,s,icomm,Q)
  call check_A_QR(n,s,icomm,A0,Q,R)

end subroutine check_qr_c

subroutine set_neighbor_ranks_c( &
     npe_x, npe_y, npe_z, rank_x, rank_y, rank_z, &
     num_neighbor_ranks, &
     neighbor_ranks, max_num_neighbor_ranks)
  use ddm_setup
  implicit none

  integer npe_x, npe_y, npe_z
  integer max_num_neighbor_ranks
  integer rank_x, rank_y, rank_z
  integer num_neighbor_ranks
  integer neighbor_ranks( max_num_neighbor_ranks )

  call set_neighbor_ranks( &
       npe_x, npe_y, npe_z, rank_x, rank_y, rank_z, &
       num_neighbor_ranks, &
       neighbor_ranks, num_neighbor_ranks)

end subroutine set_neighbor_ranks_c

subroutine set_rank_xyz_c( &
     myrank, npe_x, npe_y, npe_z, rank_x, rank_y, rank_z)
  use ddm_setup
  implicit none
  integer myrank
  integer npe_x
  integer npe_y
  integer npe_z
  integer rank_x
  integer rank_y
  integer rank_z
  call set_rank_xyz( myrank, npe_x, npe_y, npe_z, rank_x, rank_y, rank_z )

end subroutine set_rank_xyz_c

subroutine set_grids_decompostion_c( &
     npe_x, npe_y, npe_z, rank_x, rank_y, rank_z, &
     ncomp_grids, gnum_grids, ndiv_grids, lnum_grids, &
     ista_grids, margin, mdist )
  use ddm_setup
  use PARCEL_integer
  implicit none

  integer npe_x, npe_y, npe_z, rank_x, rank_y, rank_z
  integer ncomp_grids
  integer ndiv_grids( ncomp_grids )
  integer gnum_grids( ncomp_grids )
  integer lnum_grids( ncomp_grids )
  integer ista_grids( ncomp_grids )
  integer margin
  integer(kind=intType) mdist

  call set_grids_decompostion( &
       npe_x, npe_y, npe_z, rank_x, rank_y, rank_z, &
       ncomp_grids, gnum_grids, ndiv_grids, lnum_grids, &
       ista_grids, margin, mdist )

end subroutine set_grids_decompostion_c

subroutine make_poisson_3d_ddm_c( &
     gnx,gny,gnz, &
     ndim, n_grid,  &
     istart_grid, &
     margin, nnd, m, &
     offset_dim, offset, &
     matA, vecB, vecX, vecX_exact,&
     rank_x, rank_y, rank_z, &
     npe_x, npe_y, npe_z )

  use PARCEL_integer
  use ddm_setup
  implicit none
  integer(kind=intType) gnx,gny,gnz
  integer ndim, margin
  integer(kind=intType) nnd, m
  integer rank_x, rank_y, rank_z, npe_x, npe_y, npe_z
  integer n_grid( ndim )
  integer istart_grid( ndim )
  integer offset_dim( nnd, ndim )
  integer offset( nnd )
  real*8 matA( nnd*m )
  real*8 vecX(m), vecX_exact(m), vecB(m)

  call make_poisson_3D_DDM( &
       gnx,gny,gnz, &
       ndim, n_grid,  &
       istart_grid, &
       margin, nnd, m, &
       offset_dim, offset, &
       matA, vecB, vecX, vecX_exact,&
       rank_x, rank_y, rank_z, &
       npe_x, npe_y, npe_z )

end subroutine make_poisson_3d_ddm_c


#if 0
subroutine parcel_ddm_free_c(ipreflag)
  use init_ddm_param,only : parcel_ddm_free
  implicit none
  integer ipreflag
  call parcel_ddm_free(ipreflag)
end subroutine parcel_ddm_free_c
#endif

#if 0
subroutine get_eigenpair_ddm_ksp_c(&
       ncount,&
       EigenValueReal,EigenValueImag, &
       EigenVector, &
       EigenErrorNorm, &
       EmaxID,EminID)
  use init_ddm_param,only : get_EigenPair_DDM_KSP,mat
  implicit none
  integer ncount
  real*8 EigenVector( mat%nrow ) 
  real*8 EigenValueReal(ncount) 
  real*8 EigenValueImag(ncount) 
  real*8 EigenErrorNorm(ncount) 
  integer EmaxID,EminID

  call get_EigenPair_DDM_KSP( &
          ncount,&
          EigenValueReal,EigenValueImag, &
          EigenVector, &
          EigenErrorNorm, &
          EmaxID,EminID)

end subroutine get_eigenpair_ddm_ksp_c
#endif

#if 0
subroutine get_eigenpair_count_ddm_ksp_c(ncount)
  use init_ddm_param,only : get_EigenPair_count_DDM_KSP

  implicit none
  integer ncount
  call get_EigenPair_count_DDM_KSP(ncount)

end subroutine get_eigenpair_count_ddm_ksp_c
#endif

module PARCEL_structure
  use PARCEL_integer
  use krylov
  implicit none
  type(KSP) :: method_CRS
  type(KSP) :: method_DIA
end module PARCEL_structure

subroutine get_eigenpair_count_ksp_c(ncount)
  use PARCEL_structure
  implicit none
  integer ncount
  call get_EigenPair_count_KSP(ncount,method_CRS)
end subroutine get_eigenpair_count_ksp_c

subroutine get_eigenpair_ksp_c( &
       ncount,&
       EigenValueReal,EigenValueImag,&
       EigenVector,&
       EigenErrorNorm,&
       EmaxID,EminID)

  use PARCEL_structure
  implicit none

  integer ncount
  real*8 EigenValueReal(ncount)
  real*8 EigenValueImag(ncount)
  real*8 EigenVector(ncount)
  real*8 EigenErrorNorm(ncount)
  integer EmaxID,EminID

  call get_EigenPair_KSP(&
       ncount,&
       EigenValueReal,EigenValueImag,&
       EigenVector,&
       EigenErrorNorm,&
       EmaxID,EminID,method_CRS)

end subroutine get_eigenpair_ksp_c


subroutine parcel_setting_c(& 
     ipreflag, &
     Neuman_itr, &
     addL, &
     iflagAS, &
     itrmax, &
     rtolmax, &
     abstolmax, &
     solver, &
     ILU_method, &
     iovlflag, &
     GMRES_m, &
     GMRES_GSflag, &
     CAGMRES_sstep, &
     CAGMRES_tstep, &
     CAGMRES_basis, &
     CAGMRES_QRflag, &
     CBCG_kstep, &
     CBCG_Eigenflag, &
     power_method_itrmax, &
     CAARNOLDI_sstep, &
     CAARNOLDI_tstep, &
     CAARNOLDI_basis, &
     CAARNOLDI_QRflag, &
     precon_thblock, &
     independ_nvec &
     )
  use PARCEL_integer
  use krylov
  use krylov_kernel
  use PARCEL_structure
  implicit none

  integer ipreflag
  integer Neuman_itr
  integer addL
  integer iflagAS
  integer itrmax
  real*8  rtolmax
  real*8  abstolmax
  integer solver

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

  integer precon_thblock
  integer independ_nvec

  call PARCEL_KSP_Default_Setting(method_CRS)
  call set_KSP_rtolmax(method_CRS,rtolmax)
  call set_KSP_abstolmax(method_CRS,abstolmax)
  call set_KSP_solver(method_CRS,solver)
  call set_KSP_iovlflag(method_CRS,iovlflag)
  call set_KSP_iflagAS(method_CRS,iflagAS)
  call set_KSP_Neumanitr(method_CRS,Neuman_itr)
  call set_KSP_independ_nvec(method_CRS,independ_nvec)
  call set_KSP_GMRES_GSflag(method_CRS,GMRES_GSflag)
  call set_KSP_GMRES_m(method_CRS,GMRES_m)
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
  call set_KSP_ipreflag(method_CRS,ipreflag)  
  call set_KSP_itrmax(method_CRS,itrmax)  
  call set_KSP_ILU_method(method_CRS,ILU_method)  
  call set_KSP_addL(method_CRS,addL)  
  call set_KSP_precon_thblock(method_CRS,precon_thblock)

! -----

  call PARCEL_KSP_Default_Setting(method_DIA)
  call set_KSP_rtolmax(method_DIA,rtolmax)
  call set_KSP_abstolmax(method_DIA,abstolmax)
  call set_KSP_solver(method_DIA,solver)
  call set_KSP_iovlflag(method_DIA,iovlflag)
  call set_KSP_iflagAS(method_DIA,iflagAS)
  call set_KSP_Neumanitr(method_DIA,Neuman_itr)
  call set_KSP_independ_nvec(method_DIA,independ_nvec)
  call set_KSP_GMRES_GSflag(method_DIA,GMRES_GSflag)
  call set_KSP_GMRES_m(method_DIA,GMRES_m)
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
  call set_KSP_ipreflag(method_DIA,ipreflag)  
  call set_KSP_itrmax(method_DIA,itrmax)  
  call set_KSP_ILU_method(method_DIA,ILU_method)  
  call set_KSP_addL(method_DIA,addL)  
  call set_KSP_precon_thblock(method_DIA,precon_thblock)

end subroutine parcel_setting_c


subroutine set_ksp_crs_c( &
       icomm, &
       vecx, &
       vecb, &
       n,gn,nnz,istart, &
       crsA, &
       crsRow_ptr, &
       crsCol)
  use krylov_kernel,only:set_KSP_CRS
  use PARCEL_structure
  implicit none

  integer icomm
  integer(kind=intType) n,gn
  integer(kind=intType) nnz
  integer(kind=intType) istart
  real*8 vecx(n) 
  real*8 vecb(n) 
  real*8 crsA(nnz)
  integer(kind=intType) crsRow_ptr(n+1) 
  integer(kind=intType) crsCol(nnz) 

  call set_KSP_CRS(&
       icomm, &
       vecx, &
       vecb, &
       n,gn,nnz,istart, &
       crsA, &
       crsRow_ptr, &
       crsCol, &
       method_CRS)
end subroutine set_ksp_crs_c

subroutine set_ksp_dia_c( &
       icomm, &
       vecx, &
       vecb, &
       n,gn,istart, &
       diaA, &
       offset,nnd)

  use krylov_kernel,only:set_KSP_DIA
  use PARCEL_structure
  implicit none
  integer icomm
  integer(kind=intType) n,gn
  integer(kind=intType) nnd
  integer(kind=intType) istart
  real*8 vecx(n) 
  real*8 vecb(n) 
  real*8 diaA(nnd*n)
  integer(kind=intType) offset(nnd) 

  call set_KSP_DIA(&
       icomm, &
       vecx, &
       vecb, &
       n,gn,istart, &
       diaA, &
       offset,nnd, &
       method_DIA)

end subroutine set_ksp_dia_c

subroutine get_vecx_ksp_crs_c(n,vecx)
  use krylov_kernel,only:get_vecx_KSP
  use PARCEL_structure
  implicit none
  integer(kind=intType) n
  real*8 vecx(n) 

  call get_vecx_KSP(n,vecx,method_CRS)

end subroutine get_vecx_ksp_crs_c

subroutine start_ksp_solver_crs_c()
  use krylov_kernel,only:START_KSP_SOLVER
  use PARCEL_structure
  implicit none
  call START_KSP_SOLVER(method_CRS) 
end subroutine start_ksp_solver_crs_c

subroutine start_ksp_solver_dia_c()
  use krylov_kernel,only:START_KSP_SOLVER
  use PARCEL_structure
  implicit none
  call START_KSP_SOLVER(method_DIA) 

end subroutine start_ksp_solver_dia_c

subroutine get_reshistory_ksp_dia_c(niter,reshistory_DIA)
  use PARCEL_structure
  use krylov_kernel,only:get_reshistory_KSP
  implicit none
  integer niter
  real*8 reshistory_DIA(niter)
  call get_reshistory_KSP(niter,reshistory_DIA,method_DIA)
end subroutine get_reshistory_ksp_dia_c

subroutine get_reshistory_ksp_crs_c(niter,reshistory_CRS)
  use PARCEL_structure
  use krylov_kernel,only:get_reshistory_KSP
  implicit none
  integer niter
  real*8 reshistory_CRS(niter)
  call get_reshistory_KSP(niter,reshistory_CRS,method_CRS)
end subroutine get_reshistory_ksp_crs_c


subroutine parcel_dcg_ddm_c(&
       icomm, &
       vecx, &
       vecb, &
       n, &
       nnd, offset, diaA, &
       num_neighbor_ranks, neighbor_ranks, margin, ndim,&
       npe_dim, n_grids, offset_dim,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS, &      
       niter, rtolmax,&
       reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr, &
       independ_nvec, &
       nBlock, &
       iret &
       )
  use c_wrapper
  use iso_c_binding
  use ddm_wrapper
    use PARCEL_integer
implicit none

  integer icomm
  integer n
#ifdef CUDA_SOLVER
  type(devicemem) vecx
  type(devicemem) vecb
  integer nnd
  type(devicemem) offset
  type(devicemem) diaA
#else
  real*8 vecx(n)
  real*8 vecb(n)
  integer nnd
  integer offset(nnd)
  real*8  diaA(nnd*n)
#endif
 
  integer num_neighbor_ranks
  integer neighbor_ranks( num_neighbor_ranks )
  integer margin
  integer ndim
  integer npe_dim(ndim)
  integer n_grids(ndim)
  integer offset_dim(nnd,ndim)

  integer div_direc_th

  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  integer niter 
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag
  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer iret

  integer(kind=intType) n_
  n_ = n

  call parcel_dcg_ddm(&
       icomm, &
#ifdef CUDA_SOLVER
       vecx%mem, &
       vecb%mem, &
       n_, &
       nnd, offset%mem, diaA%mem, &
#else
       vecx, &
       vecb, &
       n_, &
       nnd, offset, diaA, &
#endif
       num_neighbor_ranks, neighbor_ranks, margin, ndim,&
       npe_dim, n_grids, offset_dim,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS, &      
       niter, rtolmax,&
       reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr, &
       independ_nvec, &
       nBlock, &
       iret &
       )

end subroutine parcel_dcg_ddm_c

subroutine parcel_dbicgstab_ddm_c( &
       icomm, &
       vecx, &
       vecb, &
       n, &
       nnd, offset, diaA, &
       num_neighbor_ranks, neighbor_ranks, margin, ndim,&
       npe_dim, n_grids, offset_dim,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS, &      
       niter, rtolmax,&
       reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr, &
       independ_nvec, &
       nBlock, &
       iret &
       )
  use c_wrapper
  use iso_c_binding
  use ddm_wrapper
    use PARCEL_integer
implicit none

  integer icomm
  integer n
#ifdef CUDA_SOLVER
  type(devicemem) vecx
  type(devicemem) vecb
  integer nnd
  type(devicemem) offset
  type(devicemem) diaA
#else
  real*8 vecx(n)
  real*8 vecb(n)
  integer nnd
  integer offset(nnd)
  real*8  diaA(nnd*n)
#endif 
  integer num_neighbor_ranks
  integer neighbor_ranks( num_neighbor_ranks )
  integer margin
  integer ndim
  integer npe_dim(ndim)
  integer n_grids(ndim)
  integer offset_dim(nnd,ndim)

  integer div_direc_th

  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  integer niter 
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag
  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer iret
  integer(kind=intType) n_
  n_ = n

  call parcel_dbicgstab_ddm(&
       icomm, &
#ifdef CUDA_SOLVER
       vecx%mem, &
       vecb%mem, &
       n_, &
       nnd, offset%mem, diaA%mem, &
#else
       vecx, &
       vecb, &
       n_, &
       nnd, offset, diaA, &
#endif
       num_neighbor_ranks, neighbor_ranks, margin, ndim,&
       npe_dim, n_grids, offset_dim,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS, &      
       niter, rtolmax,&
       reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr, &
       independ_nvec, &
       nBlock, &
       iret &
       )

end subroutine parcel_dbicgstab_ddm_c

subroutine parcel_dgmres_ddm_c( &
       icomm, &
       vecx, &
       vecb, &
       n, &
       nnd, offset, diaA, &
       num_neighbor_ranks, neighbor_ranks, margin, ndim,&
       npe_dim, n_grids, offset_dim,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS, &      
       niter, rtolmax,&
       reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr, &
       independ_nvec, &
       nBlock, &
       restart ,&
       gmres_gsflag ,&
       iret &
       )
  use c_wrapper
  use iso_c_binding
  use ddm_wrapper
    use PARCEL_integer
implicit none

  integer icomm
  integer n
#ifdef CUDA_SOLVER
  type(devicemem) vecx
  type(devicemem) vecb
  integer nnd
  type(devicemem) offset
  type(devicemem) diaA
#else
  real*8 vecx(n)
  real*8 vecb(n)
  integer nnd
  integer offset(nnd)
  real*8  diaA(nnd*n)
#endif 
  integer num_neighbor_ranks
  integer neighbor_ranks( num_neighbor_ranks )
  integer margin
  integer ndim
  integer npe_dim(ndim)
  integer n_grids(ndim)
  integer offset_dim(nnd,ndim)

  integer div_direc_th

  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  integer niter 
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag
  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer restart
  integer gmres_gsflag
  integer iret
  integer(kind=intType) n_
  n_ = n

  call parcel_dgmres_ddm(&
       icomm, &
#ifdef CUDA_SOLVER
       vecx%mem, &
       vecb%mem, &
       n_, &
       nnd, offset%mem, diaA%mem, &
#else
       vecx, &
       vecb, &
       n_, &
       nnd, offset, diaA, &
#endif
       num_neighbor_ranks, neighbor_ranks, margin, ndim,&
       npe_dim, n_grids, offset_dim,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS, &      
       niter, rtolmax,&
       reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr, &
       independ_nvec, &
       nBlock, &
       restart ,&
       gmres_gsflag ,&
       iret &
       )

end subroutine parcel_dgmres_ddm_c

subroutine parcel_ddcg_ddm_c(&
       icomm, &
       vecx,vecx_lo, &
       vecb,vecb_lo, &
       n, &
       nnd, offset, diaA, diaA_lo, &
       num_neighbor_ranks, neighbor_ranks, margin, ndim,&
       npe_dim, n_grids, offset_dim,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS, &      
       niter, rtolmax,&
       reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr, &
       independ_nvec, &
       nBlock, &
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION, &
       iret &
       )
  use iso_c_binding
  use ddm_wrapper
    use PARCEL_integer
implicit none

  integer icomm
  integer n
  real*8 vecx(n)
  real*8 vecb(n)
  integer nnd
  integer offset(nnd)
  real*8  diaA(nnd*n)

  real*8 vecx_lo(n)
  real*8 vecb_lo(n)
  real*8  diaA_lo(nnd*n)
 
  integer num_neighbor_ranks
  integer neighbor_ranks( num_neighbor_ranks )
  integer margin
  integer ndim
  integer npe_dim(ndim)
  integer n_grids(ndim)
  integer offset_dim(nnd,ndim)

  integer div_direc_th

  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  integer niter 
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag
  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer PRECISION_A
  integer PRECISION_b
  integer PRECISION_x
  integer PRECISION_PRECONDITION
  integer iret
  integer(kind=intType) n_
  n_ = n

  call parcel_ddcg_ddm(&
       icomm, &
       vecx,vecx_lo, &
       vecb,vecb_lo, &
       n_, &
       nnd, offset, diaA, diaA_lo, &
       num_neighbor_ranks, neighbor_ranks, margin, ndim,&
       npe_dim, n_grids, offset_dim,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS, &      
       niter, rtolmax,&
       reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr, &
       independ_nvec, &
       nBlock, &
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION, &
       iret &
       )

end subroutine parcel_ddcg_ddm_c

subroutine parcel_ddbicgstab_ddm_c( &
       icomm, &
       vecx,vecx_lo, &
       vecb,vecb_lo, &
       n, &
       nnd, offset, diaA, diaA_lo, &
       num_neighbor_ranks, neighbor_ranks, margin, ndim,&
       npe_dim, n_grids, offset_dim,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS, &      
       niter, rtolmax,&
       reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr, &
       independ_nvec, &
       nBlock, &
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION, &
       iret &
       )
  use iso_c_binding
  use ddm_wrapper
  use PARCEL_integer
  implicit none

  integer icomm
  integer n

  real*8 vecx(n)
  real*8 vecb(n)
  integer nnd
  integer offset(nnd)
  real*8  diaA(nnd*n)
  real*8 vecx_lo(n)
  real*8 vecb_lo(n)
  real*8  diaA_lo(nnd*n)

  integer num_neighbor_ranks
  integer neighbor_ranks( num_neighbor_ranks )
  integer margin
  integer ndim
  integer npe_dim(ndim)
  integer n_grids(ndim)
  integer offset_dim(nnd,ndim)

  integer div_direc_th

  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  integer niter 
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag
  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer PRECISION_A
  integer PRECISION_b
  integer PRECISION_x
  integer PRECISION_PRECONDITION
  integer iret
  integer(kind=intType) n_
  n_ = n

  call parcel_ddbicgstab_ddm(&
       icomm, &
       vecx,vecx_lo, &
       vecb,vecb_lo, &
       n_, &
       nnd, offset, diaA, diaA_lo, &
       num_neighbor_ranks, neighbor_ranks, margin, ndim,&
       npe_dim, n_grids, offset_dim,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS, &      
       niter, rtolmax,&
       reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr, &
       independ_nvec, &
       nBlock, &
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION, &
       iret &
       )

end subroutine parcel_ddbicgstab_ddm_c

subroutine parcel_ddgmres_ddm_c( &
       icomm, &
       vecx,vecx_lo, &
       vecb,vecb_lo, &
       n, &
       nnd, offset, diaA, diaA_lo, &
       num_neighbor_ranks, neighbor_ranks, margin, ndim,&
       npe_dim, n_grids, offset_dim,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS, &      
       niter, rtolmax,&
       reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr, &
       independ_nvec, &
       nBlock, &
       restart ,&
       gmres_gsflag ,&
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION, &
       iret &
       )
  use iso_c_binding
  use ddm_wrapper
    use PARCEL_integer
implicit none

  integer icomm
  integer n
  real*8 vecx(n)
  real*8 vecb(n)
  integer nnd
  integer offset(nnd)
  real*8  diaA(nnd*n)
  real*8 vecx_lo(n)
  real*8 vecb_lo(n)
  real*8  diaA_lo(nnd*n)
 
  integer num_neighbor_ranks
  integer neighbor_ranks( num_neighbor_ranks )
  integer margin
  integer ndim
  integer npe_dim(ndim)
  integer n_grids(ndim)
  integer offset_dim(nnd,ndim)

  integer div_direc_th

  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  integer niter 
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag
  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer restart
  integer gmres_gsflag
  integer PRECISION_A
  integer PRECISION_b
  integer PRECISION_x
  integer PRECISION_PRECONDITION 
  integer iret
  integer(kind=intType) n_
  n_ = n

  call parcel_ddgmres_ddm(&
       icomm, &
       vecx,vecx_lo, &
       vecb,vecb_lo, &
       n_, &
       nnd, offset, diaA, diaA_lo, &
       num_neighbor_ranks, neighbor_ranks, margin, ndim,&
       npe_dim, n_grids, offset_dim,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS, &      
       niter, rtolmax,&
       reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr, &
       independ_nvec, &
       nBlock, &
       restart ,&
       gmres_gsflag ,&
       PRECISION_A, &
       PRECISION_b, &
       PRECISION_x, &
       PRECISION_PRECONDITION, &
       iret &
       )

end subroutine parcel_ddgmres_ddm_c

subroutine parcel_dcagmres_ddm_c( &
       icomm, &
       vecx, &
       vecb, &
       n, &
       nnd, offset, diaA, &
       num_neighbor_ranks, neighbor_ranks, margin, ndim,&
       npe_dim, n_grids, offset_dim,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS, &      
       niter, rtolmax,&
       reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr, &
       independ_nvec, &
       nBlock, &
       cagmres_sstep ,&
       cagmres_tstep ,&
       cagmres_basis ,&
       cagmres_qrflag ,&
       iret &
       )
  use c_wrapper
  use iso_c_binding
  use ddm_wrapper
    use PARCEL_integer
implicit none

  integer icomm
  integer n
#ifdef CUDA_SOLVER
  type(devicemem) vecx
  type(devicemem) vecb
  integer nnd
  type(devicemem) offset
  type(devicemem) diaA
#else
  real*8 vecx(n)
  real*8 vecb(n)
  integer nnd
  integer offset(nnd)
  real*8  diaA(nnd*n)
#endif 
  integer num_neighbor_ranks
  integer neighbor_ranks( num_neighbor_ranks )
  integer margin
  integer ndim
  integer npe_dim(ndim)
  integer n_grids(ndim)
  integer offset_dim(nnd,ndim)

  integer div_direc_th

  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  integer niter 
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag
  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer cagmres_sstep
  integer cagmres_tstep
  integer cagmres_basis
  integer cagmres_qrflag

  integer iret
  integer(kind=intType) n_
  n_ = n

  call parcel_dcagmres_ddm(&
       icomm, &
#ifdef CUDA_SOLVER
       vecx%mem, &
       vecb%mem, &
       n_, &
       nnd, offset%mem, diaA%mem, &
#else
       vecx, &
       vecb, &
       n_, &
       nnd, offset, diaA, &
#endif
       num_neighbor_ranks, neighbor_ranks, margin, ndim,&
       npe_dim, n_grids, offset_dim,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS, &      
       niter, rtolmax,&
       reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr, &
       independ_nvec, &
       nBlock, &
       cagmres_sstep ,&
       cagmres_tstep ,&
       cagmres_basis ,&
       cagmres_qrflag ,&
       iret &
       )

end subroutine parcel_dcagmres_ddm_c

subroutine parcel_dcbcg_ddm_c( &
       icomm, &
       vecx, &
       vecb, &
       n, &
       nnd, offset, diaA, &
       num_neighbor_ranks, neighbor_ranks, margin, ndim,&
       npe_dim, n_grids, offset_dim,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS, &      
       niter, rtolmax,&
       reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr, &
       independ_nvec, &
       nBlock, &
       cbcg_kstep ,&
       cbcg_eigenflag ,&
       power_method_itrmax ,&
       caarnoldi_sstep ,&
       caarnoldi_tstep ,&
       caarnoldi_basis ,&
       caarnoldi_qrflag, &
       iret &
       )
  use c_wrapper
  use iso_c_binding
  use ddm_wrapper
    use PARCEL_integer
implicit none

  integer icomm
  integer n
#ifdef CUDA_SOLVER
  type(devicemem) vecx
  type(devicemem) vecb
  integer nnd
  type(devicemem) offset
  type(devicemem) diaA
#else
  real*8 vecx(n)
  real*8 vecb(n)
  integer nnd
  integer offset(nnd)
  real*8  diaA(nnd*n)
#endif 
  integer num_neighbor_ranks
  integer neighbor_ranks( num_neighbor_ranks )
  integer margin
  integer ndim
  integer npe_dim(ndim)
  integer n_grids(ndim)
  integer offset_dim(nnd,ndim)

  integer div_direc_th

  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas
  integer niter 
  real*8  rtolmax
  real*8  abstolmax

  integer ilu_method
  integer iovlflag
  integer precon_thblock
  integer independ_nvec
  integer nBlock
  real*8 reshistory(niter)

  integer cbcg_kstep
  integer cbcg_eigenflag
  integer power_method_itrmax
  integer caarnoldi_sstep
  integer caarnoldi_tstep
  integer caarnoldi_basis
  integer caarnoldi_qrflag

  integer iret
  integer(kind=intType) n_
  n_ = n

  call parcel_dcbcg_ddm(&
       icomm, &
#ifdef CUDA_SOLVER
       vecx%mem, &
       vecb%mem, &
       n_, &
       nnd, offset%mem, diaA%mem, &
#else
       vecx, &
       vecb, &
       n_, &
       nnd, offset, diaA, &
#endif
       num_neighbor_ranks, neighbor_ranks, margin, ndim,&
       npe_dim, n_grids, offset_dim,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS, &      
       niter, rtolmax,&
       reshistory, &
       iovlflag, &
       precon_thblock, &
       Neuman_itr, &
       independ_nvec, &
       nBlock, &
       cbcg_kstep ,&
       cbcg_eigenflag ,&
       power_method_itrmax ,&
       caarnoldi_sstep ,&
       caarnoldi_tstep ,&
       caarnoldi_basis ,&
       caarnoldi_qrflag, &
       iret &
       )
  
end subroutine parcel_dcbcg_ddm_c

subroutine parcel_dcaarnoldi_ddm_c( &
       icomm, &
       vecx, &
       vecb, &
       n, &
       nnd, offset, diaA, &
       num_neighbor_ranks, neighbor_ranks, margin, ndim,&
       npe_dim, n_grids, offset_dim,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS, &      
       iovlflag, &
       precon_thblock, &
       Neuman_itr, &
       independ_nvec, &
       nBlock, &
       caarnoldi_sstep ,&
       caarnoldi_tstep ,&
       caarnoldi_basis ,&
       caarnoldi_qrflag, &
       EigenValueReal,&
       EigenValueImag,&
       EigenVector   ,&
       EigenErrorNorm,&
       EmaxID,EminID, &
       iret &
       )
  use c_wrapper
  use iso_c_binding
  use ddm_wrapper
    use PARCEL_integer
implicit none

  integer icomm
  integer n
#ifdef CUDA_SOLVER
  type(devicemem) vecx
  type(devicemem) vecb
  integer nnd
  type(devicemem) offset
  type(devicemem) diaA
#else
  real*8 vecx(n)
  real*8 vecb(n)
  integer nnd
  integer offset(nnd)
  real*8  diaA(nnd*n)
#endif 
  integer num_neighbor_ranks
  integer neighbor_ranks( num_neighbor_ranks )
  integer margin
  integer ndim
  integer npe_dim(ndim)
  integer n_grids(ndim)
  integer offset_dim(nnd,ndim)

  integer div_direc_th

  integer ipreflag
  integer Neuman_itr
  integer addl
  integer iflagas

  integer ilu_method
  integer iovlflag
  integer precon_thblock
  integer independ_nvec
  integer nBlock

  integer caarnoldi_sstep
  integer caarnoldi_tstep
  integer caarnoldi_basis
  integer caarnoldi_qrflag

  real*8 EigenValueReal(caarnoldi_sstep * caarnoldi_tstep)
  real*8 EigenValueImag(caarnoldi_sstep * caarnoldi_tstep)
  real*8 EigenVector(caarnoldi_sstep * caarnoldi_tstep*n)
  real*8 EigenErrorNorm(caarnoldi_sstep * caarnoldi_tstep)
  integer EmaxID,EminID

  integer iret

  integer(kind=intType) n_
  n_ = n
  call parcel_dcaarnoldi_ddm(&
       icomm, &
#ifdef CUDA_SOLVER
       vecx%mem, &
       vecb%mem, &
       n_, &
       nnd, offset%mem, diaA%mem, &
#else
       vecx, &
       vecb, &
       n_, &
       nnd, offset, diaA, &
#endif
       num_neighbor_ranks, neighbor_ranks, margin, ndim,&
       npe_dim, n_grids, offset_dim,&
       div_direc_th, &
       ipreflag, &
       addL, iflagAS, &      
       iovlflag, &
       precon_thblock, &
       Neuman_itr, &
       independ_nvec, &
       nBlock, &
       caarnoldi_sstep ,&
       caarnoldi_tstep ,&
       caarnoldi_basis ,&
       caarnoldi_qrflag, &
       EigenValueReal,&
       EigenValueImag,&
       EigenVector   ,&
       EigenErrorNorm,&
       EmaxID,EminID, &
       iret &
       )

end subroutine parcel_dcaarnoldi_ddm_c

