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

void make_poisson2d_matrix_crs_get_nnd_(
			int *n, int*nnz,int *crsCol,int *crsRow_ptr,int *istart,int *iend,int *nnd);

void   make_poisson2d_matrix_dia_get_offset_(
			    int *offset,int *n,int *nnz,int *crsCol,int *crsRow_ptr,int *istart,int *iend,int *nnd
			    );



void make_poisson2d_matrix_dia_from_CRS_(
			 double *crsA,int *crsCol,int *crsRow_ptr,double *diaA,int *offset,int *nnd,int *n,int *nnz,int *istart);
			
void data_input(
			    int *gnx,int *gny,int *MAX_NITER, double *rtolmax,
			    int *ipreflag,int *ILU_method,int *iovlflag,int *ityp_solver,
			    int *addL,int *iflagAS,int *matA_format,int *GMRES_GSflag,
			    int *restart,int *CBCG_kstep,int *CBCG_Eigenflag,int *power_method_itrmax,
			    int *CAARNOLDI_sstep,int *CAARNOLDI_tstep,int *CAARNOLDI_basis,
			    int *CAARNOLDI_QRflag,int *CAGMRES_sstep,int *CAGMRES_tstep,
			    int *CAGMRES_basis,int *CAGMRES_QRflag,char *fileName); 

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
  int MAX_NITER           = 100    ; // Maximum number of iterations
  double rtolmax          = 1.0e-9 ; // Convergence criterion (Relative error) ||b-Ax||/||b||  
  double abstolmax        = 0.0    ; // Convergence criterion (Absolute error) ||b-Ax|| 
  int ipreflag            = 2      ; // Preconditioner :: (0:None,1:Point Jacobi,2:Block Jacobi,3:Additive schwarz)
  int addL                = 1      ; // Overlapping width of Additive schwarz
  int iflagAS             = 4      ; // Averaging method in Additive schwarz :: (1:BASIC,2:RESTRICT,3:INTERPOLATE,4:NONE)
  int ILU_method          = 1      ; // Incomplete LU factorization :: (0:ILU(0),1:対角成分が一致[DILU],2:行内の要素和が一致[MLUD])
  int matA_format         = 1      ; //
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

  char fileName[70]="./";
  
  //if (argv[1]!=NULL){
  //    strcat(fileName,"templete/templete_input/input_list");
  //    strcat(fileName,argv[1]);
  //   }
  //else{
    strcat(fileName,"input_list");
  //}

  data_input(
    &gnx,&gny,&MAX_NITER,&rtolmax,
    &ipreflag,&ILU_method,&iovlflag,&ityp_solver,
    &addL,&iflagAS,&matA_format,&GMRES_GSflag,
    &restart,&CBCG_kstep,&CBCG_Eigenflag,&power_method_itrmax,
    &CAARNOLDI_sstep,&CAARNOLDI_tstep,&CAARNOLDI_basis,
    &CAARNOLDI_QRflag,&CAGMRES_sstep,&CAGMRES_tstep,
    &CAGMRES_basis,&CAGMRES_QRflag,fileName);

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
  //make dia start
  double *diaA;
  int *offset;
  int nnd;
 if (matA_format==2){
       make_poisson2d_matrix_dia_get_nnd_( 
        &n,&nnz,crsCol,crsRow_ptr,&istart,&iend,&nnd);
       diaA   = (double *)malloc(sizeof(double)*(nnd*n));
       offset = (int *)malloc(sizeof(int)*nnd);
       for (i=0;i<nnd*n;i++){
         diaA[i] = 0.0;
       }
       for (i=0;i<nnd;i++){
         offset[i] = 0;
       }
       make_poisson2d_matrix_dia_get_offset_(
        offset,&n,&nnz,crsCol,crsRow_ptr,&istart,&iend,&nnd);
       for (i=0;i<nnd;i++){
       //  printf("offset[%d]=%d\n",i,offset[i]);
       }
       make_poisson2d_matrix_dia_from_crs_(
        crsA,crsCol,crsRow_ptr,diaA,offset,&nnd,&n,&nnz,&istart);
 }
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
  make_vecb_(&gnx,&gny,vecb,vecx,&istart,&iend);
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

  if (matA_format==1){
    set_ksp_crs_wrapper_( 
		       &icomm_Fort,  
		       vecx,vecb,  
		       &n,&gn,&nnz,&istart,  
		       crsA,crsRow_ptr,crsCol,  
		       &MAX_NITER,&ipreflag,&ILU_method,&addL,  
		       &ino); 
 }else{
    set_ksp_dia_wrapper_(
           &icomm_Fort,
           vecx,vecb,
           &n,&gn,&istart,
           diaA,offset,&nnd,
           &MAX_NITER,&ipreflag,&ILU_method,&addL,
           &ino);
 }
  start_ksp_solver_wrapper_(&ino);
  double *reshistory;
  reshistory  = (double *)malloc(sizeof(double)*MAX_NITER);
  for(i=0;i<MAX_NITER;i++){
    reshistory[i]=0.0;
  }
  int niter;
  get_reshistory_ksp_wrapper_(&niter,reshistory,&ino);
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0){
    for(i=0;i<niter;i++){
    //    for(i=0;i<MAX_NITER;i++){
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



void data_input(
     int *gnx,int *gny,int *MAX_NITER, double *rtolmax,
     int *ipreflag,int *ILU_method,int *iovlflag,int *ityp_solver,
     int *addL,int *iflagAS,int *matA_format,int *GMRES_GSflag,
     int *restart,int *CBCG_kstep,int *CBCG_Eigenflag,int *power_method_itrmax,
     int *CAARNOLDI_sstep,int *CAARNOLDI_tstep,int *CAARNOLDI_basis,
     int *CAARNOLDI_QRflag,int *CAGMRES_sstep,int *CAGMRES_tstep,
     int *CAGMRES_basis,int *CAGMRES_QRflag,char *fileName){
        int equal   = '=';
        int sup     = '!';
        int word_len;
        int i,j;
        char *ret,*ret2;
        int range_st,range_en,range_le;
        char *st_exist_check;
        char modi_text[100];
        char cur_sent[200];
        FILE *fp;
        int int_input;
        int flg=0;
        double float_input;
        char *list[24]={
          "gnx",
          "gny",
          "MAX_NITER",
          "rtolmax",
          "ipreflag",
          "ILU_method",
          "iovlflag",
          "ityp_solver",
          "addL",
          "iflagAS",
          "matA_format",
          "GMRES_GSflag",
          "restart",
          "CBCG_kstep",
          "CBCG_Eigenflag",
          "power_method_itrmax",
          "CAARNOLDI_sstep",
          "CAARNOLDI_tsteo",
          "CAARNOLDI_basis",
          "CAARNOLDI_QRflag",
          "CAGMRES_sstep",
          "CAGMRES_tstep",
          "CAGMRES_basis",
          "CAGMRES_QRflag",
        };
        
        fp=fopen(fileName,"r");
        if (fp==NULL){
          printf("Error::invalid file. ->%s \n",fileName);
          exit(1);
        }

        while(fscanf(fp,"%[^\n]%*c",cur_sent)!=EOF){
          if (cur_sent==NULL){
            break;
          }
          for (j=0;j<24;j++){
            st_exist_check=strstr(cur_sent,list[j]);
            if (st_exist_check==NULL){
              continue;
            }else{
            }
            if (st_exist_check!=NULL){
              ret      = strchr(cur_sent,equal);
              ret2     = strchr(cur_sent,sup);
              word_len = strlen(cur_sent);
              if (ret!=NULL){
                range_st=ret-cur_sent+1;
              }else{
                printf("not found Delimiter\n");
                exit(1);
              }
              if (ret2!=NULL){
                range_en=ret2-cur_sent;
              }else{
                range_en=word_len;
              }
              range_le  = range_en-range_st;

              strncpy(modi_text,cur_sent+range_st,range_le);
              modi_text[range_le]='\0';
              if (j!=3){
                int_input=strtol(modi_text,NULL,10);
              }else{
                float_input=strtof(modi_text,NULL);
              }
              switch(j){
                case 0:
                  *gnx                  = int_input;
                  break;
                case 1:
                  *gny                  = int_input;
                  break;
                case 2:
                  *MAX_NITER            = int_input;
                  break;
                case 3:
                  *rtolmax              = float_input;
                  break;
                case 4:
                  *ipreflag             = int_input;
                  break;
                case 5:
                  *ILU_method           = int_input;
                  break;
                case 6:
                  *iovlflag             = int_input;
                  break;
                case 7:
                  *ityp_solver          = int_input;
                  break;
                case 8:
                  *addL                 = int_input;
                  break;
                case 9:
                  *iflagAS              = int_input;
                  break;
                case 10:
                  *matA_format          = int_input;
                  break;
                case 11:
                  *GMRES_GSflag         = int_input;
                  break;
                case 12:
                  *restart              = int_input;
                  break;
                case 13:
                  *CBCG_kstep           = int_input;
                  break;
                case 14:
                  *CBCG_Eigenflag       = int_input;
                  break;
                case 15:
                  *power_method_itrmax  = int_input;
                  break;
                case 16:
                  *CAARNOLDI_sstep      = int_input;
                  break;
                case 17:
                  *CAARNOLDI_tstep      = int_input;
                  break;
                case 18:
                  *CAARNOLDI_basis      = int_input;
                  break;
                case 19:
                  *CAARNOLDI_QRflag     = int_input;
                  break;
                case 20:
                  *CAGMRES_sstep        = int_input;
                  break;
                case 21:
                  *CAGMRES_tstep        = int_input;
                  break;
                case 22:
                  *CAGMRES_basis        = int_input;
                  break;
                case 23:
                  *CAGMRES_QRflag       = int_input;
                  break;
                default:
                  break;
              }
              break;
            }
          }
        }
      fclose(fp);
}

