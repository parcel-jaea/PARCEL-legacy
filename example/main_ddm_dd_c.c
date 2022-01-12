#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#define PARCEL_DIA_Solver 1

#ifdef PARCEL_INT8
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
  int ityp_solver;
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
  int PRECISION_KRYLOV = 1;
  get_namelist_input_dd_( 
			 &PRECISION_A, 
			 &PRECISION_b, 
			 &PRECISION_x, 
			 &PRECISION_PRECONDITION
			  );

  int npe_x,npe_y,npe_z,div_direc_th;
  get_namelist_input_ddm_( 
			  &npe_x,
			  &npe_y,
			  &npe_z,
			  &div_direc_th
			   );
  int nBlock = 2000;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  icomm_fort=MPI_Comm_c2f(MPI_COMM_WORLD);

  parcel_int nnd = 7;

  int rank_x,rank_y,rank_z;
  set_rank_xyz_c_( 
		  &myrank,
		  &npe_x,
		  &npe_y,
		  &npe_z,
		  &rank_x,
		  &rank_y,
		  &rank_z);
  int ndim=3 ;
  int margin = 1;
  int *npe_dim = (int *)malloc(sizeof(int)*ndim);
  int *gn_grids = (int *)malloc(sizeof(int)*ndim);
  int *n_grids = (int *)malloc(sizeof(int)*ndim);
  int *istart_grids = (int *)malloc(sizeof(int)*ndim);
  for(i=0;i<ndim;i++){
    npe_dim[i] = 0;
    gn_grids[i] = 0;
    n_grids[i] = 0;
    istart_grids[i] = 0;
  }
  parcel_int m=-1;

  gn_grids[0]  = gnx;
  gn_grids[1]  = gny;
  gn_grids[2]  = gnz;
  set_grids_decompostion_c_( 
			    &npe_x, &npe_y, &npe_z, 
			    &rank_x, &rank_y, &rank_z,  
			    &ndim, gn_grids, npe_dim, n_grids,  
			    istart_grids, &margin, &m );

  int num_neighbor_ranks_max = pow(3,ndim) -1;
  int *neighbor_ranks = (int *)malloc(sizeof(int)*num_neighbor_ranks_max);

  int num_neighbor_ranks=-1;

  for(i=0;i<num_neighbor_ranks_max;i++){
    neighbor_ranks[i] = -1;
  }
  set_neighbor_ranks_c_( 
			&npe_x,
			&npe_y,
			&npe_z,
			&rank_x,
			&rank_y,
			&rank_z,
			&num_neighbor_ranks,
			neighbor_ranks,
			&num_neighbor_ranks_max);
 
  int *offset = (int *)malloc(sizeof(int)*nnd);
  for(i=0;i<nnd;i++){
    offset[i] = 0;
  }
  int *offset_dim = (int *)malloc(sizeof(int)*nnd*ndim);
  for(i=0;i<nnd*ndim;i++){
    offset_dim[i] = 0;
  }
  double *diaA    = (double *)malloc(sizeof(double)*m*nnd);
  double *diaA_lo = (double *)malloc(sizeof(double)*m*nnd);
  for(i=0;i<m*nnd;i++){
    diaA[i]    = 0.0;
    diaA_lo[i] = 0.0;
  }
  double *vecx  = (double *)malloc(sizeof(double)*m);
  double *vecx0 = (double *)malloc(sizeof(double)*m);
  double *vecb  = (double *)malloc(sizeof(double)*m);

  double *vecx_lo = (double *)malloc(sizeof(double)*m);
  double *vecb_lo = (double *)malloc(sizeof(double)*m);
  for(i=0;i<m;i++){
    vecx[i] = 0.0;
    vecx0[i] = 0.0;
    vecb[i] = 0.0;

    vecx_lo[i] = 0.0;
    vecb_lo[i] = 0.0;
  }
  double *reshistory_ddm  = (double *)malloc(sizeof(double)*max_niter);
  for(i=0;i<max_niter;i++){
    reshistory_ddm[i] = 0.0;
  }

  parcel_int nnd_ = nnd;
  make_poisson_3d_ddm_c_( 
			 &gnx,
			 &gny,
			 &gnz,
			 &ndim,
			 n_grids,
			 istart_grids,
			 &margin,
			 &nnd_,
			 &m,
			 offset_dim,
			 offset,
			 diaA,
			 vecb,
			 vecx,
			 vecx0,
			 &rank_x,
			 &rank_y,
			 &rank_z,
			 &npe_x,
			 &npe_y,
			 &npe_z);

  switch(ityp_solver){
  case 1:
    parcel_ddcg_ddm_(
		      &icomm_fort,
		      vecx,vecx_lo,
		      vecb,vecb_lo,
		      &m,
		      &nnd,
		      offset,
		      diaA,diaA_lo,
		      &num_neighbor_ranks,
		      neighbor_ranks,
		      &margin,
		      &ndim,
		      npe_dim,
		      n_grids,
		      offset_dim,
		      &div_direc_th,
		      &ipreflag,
		      &addl,
		      &iflagas,
		      &max_niter,
		      &rtolmax,
		      reshistory_ddm,
		      &iovlflag,
		      &precon_thblock,
		      &independ_nvec ,
		      &nBlock ,
		      &PRECISION_A, 
		      &PRECISION_b, 
		      &PRECISION_x, 
		      &PRECISION_PRECONDITION,
		      &iret
		      );
    break;
  case 2:
    parcel_ddbicgstab_ddm_(
		      &icomm_fort,
		      vecx,vecx_lo,
		      vecb,vecb_lo,
		      &m,
		      &nnd,
		      offset,
		      diaA,diaA_lo,
		      &num_neighbor_ranks,
		      neighbor_ranks,
		      &margin,
		      &ndim,
		      npe_dim,
		      n_grids,
		      offset_dim,
		      &div_direc_th,
		      &ipreflag,
		      &addl,
		      &iflagas,
		      &max_niter,
		      &rtolmax,
		      reshistory_ddm,
		      &iovlflag,
		      &precon_thblock,
		      &independ_nvec ,
		      &nBlock ,
		      &PRECISION_A, 
		      &PRECISION_b, 
		      &PRECISION_x, 
		      &PRECISION_PRECONDITION,
		      &iret
		      );
    break;
  case 3:
    parcel_ddgmres_ddm_(
		      &icomm_fort,
		      vecx,vecx_lo,
		      vecb,vecb_lo,
		      &m,
		      &nnd,
		      offset,
		      diaA,diaA_lo,
		      &num_neighbor_ranks,
		      neighbor_ranks,
		      &margin,
		      &ndim,
		      npe_dim,
		      n_grids,
		      offset_dim,
		      &div_direc_th,
		      &ipreflag,
		      &addl,
		      &iflagas,
		      &max_niter,
		      &rtolmax,
		      reshistory_ddm,
		      &iovlflag,
		      &precon_thblock,
		      &independ_nvec ,
		      &nBlock ,
		      &restart,
		      &gmres_gsflag,
		      &PRECISION_A, 
		      &PRECISION_b, 
		      &PRECISION_x, 
		      &PRECISION_PRECONDITION,
		      &iret
		      );

  default:
    printf("not support solver.\n");
  }

  // printf(" END PARCEL\n");
  if(myrank==0){
    printf(" residual norm \n");
    printf("     DDM    \n");
    for(i=0;i<max_niter;i++){
      printf("PARCEL : @res@ %d, %e \n",i+1,reshistory_ddm[i]);
    }
  }

  free(offset);
  free(offset_dim);
  free(diaA);
  free(vecx);
  free(vecb);
  free(vecx0);

  free(diaA_lo);
  free(vecx_lo);
  free(vecb_lo);

  free(reshistory_ddm);
  free(neighbor_ranks);
  free(npe_dim);
  free(gn_grids);
  free(n_grids);
  free(istart_grids);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}

