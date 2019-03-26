void set_parcel_wrapper_(int *imax);
void set_ksp_solver_wrapper_(int *ino,int *ityp_solver);
void set_ksp_iovlflag_wrapper_(int *ino,int *iovlflag);
void set_ksp_iflagas_wrapper_(int *ino,int *iflagAS);
void set_ksp_gmres_gsflag_wrapper_(int *ino,int *GMRES_GSflag);
void set_ksp_gmres_m_wrapper_(int *ino,int *restart);
void set_ksp_cbcg_kstep_wrapper_(int *ino,int *CBCG_kstep);
void set_ksp_cbcg_eigenflag_wrapper_(int *ino,int *CBCG_Eigenflag);
void set_ksp_power_method_itrmax_wrapper_(int *ino,int *power_method_itrmax);
void set_ksp_caarnoldi_sstep_wrapper_(int *ino,int *CAARNOLDI_sstep);
void set_ksp_caarnoldi_tstep_wrapper_(int *ino,int *CAARNOLDI_tstep);
void set_ksp_caarnoldi_basis_wrapper_(int *ino,int *CAARNOLDI_basis);
void set_ksp_caarnoldi_qrflag_wrapper_(int *ino,int *CAARNOLDI_QRflag);
void set_ksp_cagmres_sstep_wrapper_(int *ino,int *CAGMRES_sstep);
void set_ksp_cagmres_tstep_wrapper_(int *ino,int *CAGMRES_tstep);
void set_ksp_cagmres_basis_wrapper_(int *ino,int *CAGMRES_basis);
void set_ksp_cagmres_qrflag_wrapper_(int *ino,int *CAGMRES_QRflag);
void set_ksp_rtolmax_wrapper_(int *ino,double *rtolmax);
void set_ksp_abstolmax_wrapper_(int *ino,double *abstolmax);
void set_ksp_crs_wrapper_(
		       MPI_Fint *icomm_Fort,  
		       double *vecx,double *vecb,  
		       int *n,int *gn,int *nnz,int *istart,  
		       double *crsA,int *crsRow_ptr,int *crsCol,  
		       int *itrmax,int *ipreflag,int *ILU_method,int *addL,  
		       int *ino); 
void start_ksp_solver_wrapper_(int *ino);
