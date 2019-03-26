#include <stdio.h>
#include "mpi.h"
#include "PARCEL_Wrapper.h"

void make_poisson2d_matrix_crs_get_nnz_( 
				     int *nx,int *ny,int *nx0,int *ny0,int *gny,int *gnx,
				     int *nnz,
				     int *istart);
void make_poisson2d_matrix_crs_( 
			    int *nx,int *ny,int *nx0,int *ny0,int *gny,int *gnx,
			    double *crsA,int *crsCol,int *crsRow_ptr,
			    int *n,int *nnz, 
			    int *istart);

int main( int argc, char *argv[]){
  int npe;
  int rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npe);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Fint icomm_Fort;
  icomm_Fort=MPI_Comm_c2f(MPI_COMM_WORLD);

  // ----------
  int ityp_solver         = 1      ; // Algorithm :: Double precision(1:CG,2:BiCG,3:GMRES,:4:CAGMRES,:5:CBCG)
  int gnx                 = 40     ; // Problem size in x
  int gny                 = 40     ; // Problem size in y
  int itrmax              = 100    ; // Maximum number of iterations
  double rtolmax          = 1.0e-9 ; // Convergence criterion (Relative error) ||b-Ax||/||b||  
  double abstolmax        = 0.0    ; // Convergence criterion (Absolute error) ||b-Ax|| 
  int ipreflag            = 2      ; // Preconditioner :: (0:None,1:Point Jacobi,2:Block Jacobi,3:Additive schwarz)
  int addL                = 1      ; // Overlapping width of Additive schwarz
  int iflagAS             = 4      ; // Averaging method in Additive schwarz :: (1:BASIC,2:RESTRICT,3:INTERPOLATE,4:NONE)
  int ILU_method          = 1      ; // Incomplete LU factorization :: (0:ILU(0),1:対角成分が一致[DILU],2:行内の要素和が一致[MLUD])
  int iovlflag            = 0      ; // Communication overlap :: (0:None,1:Single step,2:Multiple steps)
  int GMRES_GSflag        = 2      ; // Orthogonalization method in GMRES(m) :: (1:MGS,2:CGS)
  int restart             = 10     ; // Restart length in GMRES(m)
  int CAGMRES_sstep       = 5      ; // Communication avoiding steps in CA-GMRES
  int CAGMRES_tstep       = 10     ; // Number of outer iterations in CA-GMRES
  int CAGMRES_basis       = 0      ; // Basis vector in CA-GMRES :: (0:Monomial basis,1:Newton basis)
  int CAGMRES_QRflag      = 3      ; // QR factorization in CA-GMRES :: (1:MGS,2:CGS,3:TSQR,4:Cholesky QR,5:Cholesky QR2)
  int CAARNOLDI_sstep     = 2      ; // Communication avoiding steps in CA-ARNOLDI
  int CAARNOLDI_tstep     = 5      ; // Number of outer iterations in CA-ARNOLDI
  int CAARNOLDI_basis     = 0      ; // Basis vector in CA-ARNOLDI :: (0:Monomial basis,1:Newton basis)
  int CAARNOLDI_QRflag    = 4      ; // QR factorization in CA-ARNOLDI :: (1:MGS,2:CGS,3:TSQR,4:Cholesky QR,5:Cholesky QR2)
  int CBCG_kstep          = 4      ; // Communication avoiding steps in CBCG
  int CBCG_Eigenflag      = 2      ; // Eigenvalue solver in CBCG :: (1:Power method,2:CA-ARNOLDI)
  int power_method_itrmax = 100    ; // Number of iterations in power method

  int nx     = gnx;
  int ny     = gny/npe ;
  int nx0    = 1  ;
  int ny0    = 1 + ny*rank;
  int gn     = gnx*gny;
  int n      = nx * ny;
  int istart = 1+n*rank;
  int iend   = istart + n -1;
  int nnz;

  // ----------

  make_poisson2d_matrix_crs_get_nnz_( 
				     &nx,&ny,&nx0,&ny0,&gny,&gnx,
				     &nnz,
				     &istart);
  double *crsA;
  int *crsCol;
  int *crsRow_ptr;
  crsA       = (double *)malloc(sizeof(double)*nnz);
  crsCol     = (int *)malloc(sizeof(int)*nnz);
  crsRow_ptr = (int *)malloc(sizeof(int)*(n+1));

  int i;
  for(i=0;i<nnz;i++){
    crsA[i]   = 0.0;
    crsCol[i] = 0;
  }
  for(i=0;i<n+1;i++){
    crsRow_ptr[i] = 0;
  }
  make_poisson2d_matrix_crs_( 
			    &nx,&ny,&nx0,&ny0,&gny,&gnx,
			    crsA,crsCol,crsRow_ptr,
			    &n,&nnz, 
			    &istart);

  // ----------
  int ino;
  int imax;
  double *vecx;
  double *vecb;
  vecx  = (double *)malloc(sizeof(double)*n);
  vecb  = (double *)malloc(sizeof(double)*n);

  for(i=0;i<n;i++){
    vecx[i]= 0.0;
    vecb[i]=-1.0;
  }
  imax=1;
  set_parcel_wrapper_(&imax);
  ino=0;
  set_ksp_solver_wrapper_(&ino,&ityp_solver);
  set_ksp_iovlflag_wrapper_(&ino,&iovlflag);
  set_ksp_iflagas_wrapper_(&ino,&iflagAS);
  set_ksp_gmres_gsflag_wrapper_(&ino,&GMRES_GSflag);
  set_ksp_gmres_m_wrapper_(&ino,&restart);
  set_ksp_cbcg_kstep_wrapper_(&ino,&CBCG_kstep);
  set_ksp_cbcg_eigenflag_wrapper_(&ino,&CBCG_Eigenflag);
  set_ksp_power_method_itrmax_wrapper_(&ino,&power_method_itrmax);
  set_ksp_caarnoldi_sstep_wrapper_(&ino,&CAARNOLDI_sstep);
  set_ksp_caarnoldi_tstep_wrapper_(&ino,&CAARNOLDI_tstep);
  set_ksp_caarnoldi_basis_wrapper_(&ino,&CAARNOLDI_basis);
  set_ksp_caarnoldi_qrflag_wrapper_(&ino,&CAARNOLDI_QRflag);
  set_ksp_cagmres_sstep_wrapper_(&ino,&CAGMRES_sstep);
  set_ksp_cagmres_tstep_wrapper_(&ino,&CAGMRES_tstep);
  set_ksp_cagmres_basis_wrapper_(&ino,&CAGMRES_basis);
  set_ksp_cagmres_qrflag_wrapper_(&ino,&CAGMRES_QRflag)  ;
  set_ksp_rtolmax_wrapper_(&ino,&rtolmax);
  set_ksp_abstolmax_wrapper_(&ino,&abstolmax);
  set_ksp_crs_wrapper_( 
		       &icomm_Fort,  
		       vecx,vecb,  
		       &n,&gn,&nnz,&istart,  
		       crsA,crsRow_ptr,crsCol,  
		       &itrmax,&ipreflag,&ILU_method,&addL,  
		       &ino); 
  start_ksp_solver_wrapper_(&ino);

  double *reshistory;
  reshistory  = (double *)malloc(sizeof(double)*itrmax);
  for(i=0;i<itrmax;i++){
    reshistory[i]=0.0;
  }
  int niter;
  get_reshistory_ksp_wrapper_(&niter,reshistory,&ino);
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0){
    for(i=0;i<niter;i++){
    //    for(i=0;i<itrmax;i++){
      printf("@res@ %4d , %e \n",i,reshistory[i]);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  //  test_(&icomm_Fort);
  get_vecx_ksp_wrapper_(&n,vecx,&ino);

  free_parcel_wrapper_();
  free(crsA);
  free(crsCol);
  free(crsRow_ptr);
  free(vecx);
  free(vecb);
  free(reshistory);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
