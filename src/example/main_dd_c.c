#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#ifdef CUDA_SOLVER

int main( int argc, char *argv[]){
  printf("not support dd solver \n");
  return 0;
}

#else

#define PARCEL_DIA_Solver 1

#ifdef DPARCEL_INT8
#define parcel_int long long 
#else
#define parcel_int int
#endif

int main( int argc, char *argv[]){
  int ityp_eq;
  parcel_int gnx,gny,gnz;
  int itrmax;
  int max_niter;
  parcel_int nx,ny,nz;
  parcel_int nx0,ny0,nz0;
  parcel_int n,gn;
  parcel_int istart,iend;

  parcel_int nnz;
  int npes,myrank;
  int ierr;
  parcel_int i;
  int i_;
  int ipreflag;
  int addl;
  int iflagas;
  double  rtolmax;
  double  abstolmax;
  int solver,ityp_solver;
  int iret;

  int ilu_method;
  int iovlflag;
  int gmres_gsflag;
  int cagmres_sstep;
  int cagmres_tstep;
  int cagmres_basis;
  int cagmres_qrflag;
  int cbcg_kstep;
  int cbcg_eigenflag;
  int power_method_itrmax;
  int caarnoldi_sstep;
  int caarnoldi_tstep;
  int caarnoldi_basis;
  int caarnoldi_qrflag;

  int precon_thblock;
  int Neuman_itrmax;
  int independ_nvec;
  int restart;
  int myranky,myrankz;
  int npey,npez;
  MPI_Fint icomm_fort;

  get_namelist_input_( 
		      &ityp_eq ,
		      &gnx ,
		      &gny ,
		      &gnz ,
		      &max_niter ,
		      &rtolmax ,
		      &abstolmax ,
		      &ipreflag ,
		      &Neuman_itrmax ,
		      &independ_nvec ,
		      &addl ,
		      &iflagas ,
		      &ilu_method ,
		      &iovlflag ,
		      &ityp_solver ,
		      &precon_thblock ,
		      &restart ,
		      &gmres_gsflag ,
		      &cbcg_kstep ,
		      &cbcg_eigenflag ,
		      &power_method_itrmax ,
		      &cagmres_sstep ,
		      &cagmres_tstep ,
		      &cagmres_basis ,
		      &cagmres_qrflag ,
		      &caarnoldi_sstep ,
		      &caarnoldi_tstep ,
		      &caarnoldi_basis ,
		      &caarnoldi_qrflag 
		       );

  int PRECISION_A;
  int PRECISION_b;
  int PRECISION_x;
  int PRECISION_PRECONDITION;
  get_namelist_input_dd_( 
			 &PRECISION_A, 
			 &PRECISION_b, 
			 &PRECISION_x, 
			 &PRECISION_PRECONDITION
			  );
  int nBlock = 2000;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  icomm_fort=MPI_Comm_c2f(MPI_COMM_WORLD);

  get_init_poisson3d_( 
		      &icomm_fort, 
		      &gn, 
		      &gnx, 
		      &gny, 
		      &gnz, 
		      &istart, 
		      &iend, 
		      &nx, 
		      &ny, 
		      &nz, 
		      &nx0, 
		      &ny0, 
		      &nz0, 
		      &n, 
		      &npes, 
		      &myrank, 
		      &myranky, 
		      &myrankz, 
		      &npey, 
		      &npez 
		       );

  make_poisson3d_crs_get_nnz_c_( 
				&nx,
				&ny,
				&nz,
				&nx0,
				&ny0,
				&nz0,
				&gnx,
				&gny,
				&gnz,
				&nnz, 
				&istart 
				 );
  
  double *vecx     = (double *)malloc(sizeof(double)*n);
  double *vecx_lo  = (double *)malloc(sizeof(double)*n);
  double *vecb     = (double *)malloc(sizeof(double)*n);
  double *vecb_lo  = (double *)malloc(sizeof(double)*n);
  double *vecx0 = (double *)malloc(sizeof(double)*n);
  for(i=0;i<n;i++){
    vecx[i]  = 0.0;
    vecb[i]  = 0.0;
    vecx_lo[i]  = 0.0;
    vecb_lo[i]  = 0.0;
    vecx0[i] = 0.0;
  }
  double *crsa  = (double *)malloc(sizeof(double)*nnz);
  double *crsa_lo = (double *)malloc(sizeof(double)*nnz);
  parcel_int *crscol = (parcel_int *)malloc(sizeof(parcel_int)*nnz);
  for(i=0;i<nnz;i++){
    crsa[i] = 0.0;
    crsa_lo[i] = 0.0;
    crscol[i] = 0;
  }
  parcel_int *crsrow_ptr  = (parcel_int *)malloc(sizeof(parcel_int)*(n+1));
  for(i=0;i<(n+1);i++){
    crsrow_ptr[i] = 0;
  }
  double *reshistory_crs  = (double *)malloc(sizeof(double)*max_niter);
  for(i=0;i<max_niter;i++){
    reshistory_crs[i] = 0.0;
  }
  make_poisson3d_crs_c_( 
			&nx,
			&ny,
			&nz,
			&nx0,
			&ny0,
			&nz0,
			&gnx,
			&gny,
			&gnz,
			crsa,
			crscol,
			crsrow_ptr,
			&n,
			&nnz, 
			&istart);
  make_poisson3d_right_vector_c_( 
				 &nx,
				 &ny,
				 &nz,
				 &nx0,
				 &ny0,
				 &nz0,
				 &gnx,
				 &gny,
				 &gnz,
				 vecx0,
				 vecb,
				 &istart);

  switch(ityp_solver){
  case 1:
    parcel_ddcg_(
		  &icomm_fort , 
		  vecx, vecx_lo,
		  vecb, vecb_lo,
		  &n, 
		  &gn, 
		  &nnz, 
		  &istart, 
		  crsa , crsa_lo,
		  crsrow_ptr,
		  crscol,
		  &ipreflag ,
		  &ilu_method ,
		  &addl ,		  
		  &iflagas ,
		  &max_niter ,
		  &rtolmax ,
		  reshistory_crs, 
		  &iovlflag ,
		  &precon_thblock ,
		  &independ_nvec ,
		  &nBlock,
		  &PRECISION_A,
		  &PRECISION_b,
		  &PRECISION_x,
		  &PRECISION_PRECONDITION,
		  &iret 
		  );
    break;
  case 2:
    parcel_ddbicgstab_(
		  &icomm_fort , 
		  vecx, vecx_lo,
		  vecb, vecb_lo,
		  &n, 
		  &gn, 
		  &nnz, 
		  &istart, 
		  crsa , crsa_lo,
		  crsrow_ptr,
		  crscol,
		  &ipreflag ,
		  &ilu_method ,
		  &addl ,		  
		  &iflagas ,
		  &max_niter ,
		  &rtolmax ,
		  reshistory_crs, 
		  &iovlflag ,
		  &precon_thblock ,
		  &independ_nvec ,
		  &nBlock,
		  &PRECISION_A,
		  &PRECISION_b,
		  &PRECISION_x,
		  &PRECISION_PRECONDITION,
		  &iret 
		  );
    break;
  case 3:
    parcel_ddgmres_(
		  &icomm_fort , 
		  vecx, vecx_lo,
		  vecb, vecb_lo,
		  &n, 
		  &gn, 
		  &nnz, 
		  &istart, 
		  crsa , crsa_lo,
		  crsrow_ptr,
		  crscol,
		  &ipreflag ,
		  &ilu_method ,
		  &addl ,		  
		  &iflagas ,
		  &max_niter ,
		  &rtolmax ,
		  reshistory_crs, 
		  &iovlflag ,
		  &precon_thblock ,
		  &independ_nvec ,
		  &nBlock,
		  &restart ,
		  &gmres_gsflag ,
		  &PRECISION_A,
		  &PRECISION_b,
		  &PRECISION_x,
		  &PRECISION_PRECONDITION,
		  &iret 
		  );
    break;
  default:
    printf("not support solver.\n");
  }

#if PARCEL_DIA_Solver
  parcel_int nnd=10;
  double *diaA     = (double *)malloc(sizeof(double)*n*nnd);
  double *diaA_lo  = (double *)malloc(sizeof(double)*n*nnd);
  parcel_int *offset = (parcel_int *)malloc(sizeof(parcel_int)*(nnd));
  double *reshistory_dia  = (double *)malloc(sizeof(double)*max_niter);
  for(i=0;i<n*nnd;i++){
    diaA[i] = 0.0;
    diaA_lo[i] = 0.0;
  }
  for(i=0;i<max_niter;i++){
    reshistory_dia[i] = 0.0;
  }
  crs2dia_c_( 
    &n,
    &nnz, 
    crsa , 
    crsrow_ptr, 
    crscol, 
    diaA , 
    offset, 
    &nnd,
    &istart, 
    &iend 
    );

  for(i=0;i<n;i++){
    vecx[i]  = 0.0;
    vecx_lo[i]  = 0.0;
  
  }

  switch(ityp_solver){
  case 1:
    parcel_ddcg_dia_c_(
		  &icomm_fort , 
		  vecx, vecx_lo,
		  vecb, vecb_lo,
		  &n, 
		  &gn, 
		  &istart, 
		    diaA,diaA_lo,
		    offset,
		    &nnd,
		  &ipreflag ,
		  &ilu_method ,
		  &addl ,		  
		  &iflagas ,
		  &max_niter ,
		  &rtolmax ,
		  reshistory_dia, 
		  &iovlflag ,
		  &precon_thblock ,
		  &independ_nvec ,
		  &nBlock,
		  &PRECISION_A,
		  &PRECISION_b,
		  &PRECISION_x,
		  &PRECISION_PRECONDITION,
		  &iret 
		  );
    break;
  case 2:
    parcel_ddbicgstab_dia_c_(
		  &icomm_fort , 
		  vecx, vecx_lo,
		  vecb, vecb_lo,
		  &n, 
		  &gn, 
		  &istart, 
		    diaA,diaA_lo,
		    offset,
		    &nnd,
		  &ipreflag ,
		  &ilu_method ,
		  &addl ,		  
		  &iflagas ,
		  &max_niter ,
		  &rtolmax ,
		  reshistory_dia, 
		  &iovlflag ,
		  &precon_thblock ,
		  &independ_nvec ,
		  &nBlock,
		  &PRECISION_A,
		  &PRECISION_b,
		  &PRECISION_x,
		  &PRECISION_PRECONDITION,
		  &iret 
		  );
    break;
  case 3:
    parcel_ddgmres_dia_c_(
		  &icomm_fort , 
		  vecx, vecx_lo,
		  vecb, vecb_lo,
		  &n, 
		  &gn, 
		  &istart, 
		    diaA,diaA_lo,
		    offset,
		    &nnd,
		  &ipreflag ,
		  &ilu_method ,
		  &addl ,		  
		  &iflagas ,
		  &max_niter ,
		  &rtolmax ,
		  reshistory_dia, 
		  &iovlflag ,
		  &precon_thblock ,
		  &independ_nvec ,
		  &nBlock,
		  &restart ,
		  &gmres_gsflag ,
		  &PRECISION_A,
		  &PRECISION_b,
		  &PRECISION_x,
		  &PRECISION_PRECONDITION,
		  &iret 
		  );
    break;
  default:
    printf("not support solver.\n");
  }

#endif

  if(myrank==0){

#if PARCEL_DIA_Solver
     printf("\n");
     printf(" -------------------------------\n");
     printf("\n");
     printf(" residual norm CRS matrix \n");
     printf("     CRS ddmatrix    \n");
#else
     printf("\n");
     printf(" -------------------------------\n");
     printf("\n");
     printf(" residual norm CRS matrix and DIA matrix\n");
     printf("     CRS ddmatrix    ,       DIA ddmatrix\n");
#endif

    for(i=0;i<max_niter;i++){
#if PARCEL_DIA_Solver
      printf("PARCEL : @res@ %d, %e, %e \n",i+1,reshistory_crs[i],reshistory_dia[i]);
#else
      printf("PARCEL : @res@ %d, %e \n",i+1,reshistory_crs[i]);
#endif
    }
  }
  check_exact_solution_c_(&n,&icomm_fort,vecx,vecx0);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  free(vecx);
  free(vecb);
  free(vecx_lo);
  free(vecb_lo);
  free(vecx0);
  free(crsa);
  free(crsa_lo);
  free(crscol);
  free(crsrow_ptr);
  free(reshistory_crs);
#if PARCEL_DIA_Solver

  free(diaA);
  free(offset);
  free(reshistory_dia);

#endif

  return 0;
}

#endif
