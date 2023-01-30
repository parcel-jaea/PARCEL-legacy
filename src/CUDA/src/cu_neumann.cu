#include "cu_common.h"


template <class Tdata,class Tcalc>
__global__
void cuda_get_diamatrix_neumanndiag_(
		 PARCEL_INT n     ,
		 PARCEL_INT nnd   ,
		 PARCEL_INT nBlock,
		 PARCEL_INT diag_offset,
		 Tdata *diaA,
		 Tdata *D,
		 Tdata *Dinv,
		 // PARCEL_INT *offset_, // 対角入れ替えのために利用
		 PARCEL_INT *iLs     // 対角入れ替えのために利用
		 ){

  const int thread_id   = threadIdx.x;
  const int block_id    = blockIdx.x;
  const int block_step  = gridDim.x;
  const int thread_step = blockDim.x;

  for(PARCEL_INT kk=1+block_id*nBlock;kk<=n;kk=kk+nBlock*block_step){
    PARCEL_INT kkmax = n < kk+nBlock-1 ? n : kk+nBlock-1;
    for(PARCEL_INT kid=kk+thread_id;kid<=kkmax;kid=kid+thread_step){
      PARCEL_INT j = diag_offset;
      D[kid-1]    = 1.0;
      Dinv[kid-1] = 1.0;

      // 対角要素
      Tdata diag_val = diaA[kid + (j-1)*n -1 ];
      D[kid-1] = diag_val;
      Dinv[kid-1] = 1.0;
      if( diag_val != 0.0 ){
	Dinv[kid-1] = 1.0/diag_val;
      }
      diaA[kid + (j-1)*n -1 ] = 0.0;
    }
    
    for(PARCEL_INT j = 1;j<=nnd-1;j++){
    //    for(PARCEL_INT j = 1;j<=nnd;j++){
      if(j==diag_offset){
	PARCEL_INT jend = nnd;
	for(PARCEL_INT kid=kk+thread_id;kid<=kkmax;kid=kid+thread_step){
#if 1
	  diaA[kid + (j-1)*n -1 ] = -diaA[kid + (jend-1)*n -1 ] * Dinv[kid-1];
	   iLs[kid + (j-1)*n -1 ] =   iLs[kid + (jend-1)*n -1 ]; 
#else
	  diaA[kid + (j-1)*n -1 ] = 0.0;
#endif
	}
      }else{
	for(PARCEL_INT kid=kk+thread_id;kid<=kkmax;kid=kid+thread_step){
	  diaA[kid + (j-1)*n -1 ] = -diaA[kid + (j-1)*n -1 ] * Dinv[kid-1];
	}
      }
    }

  }
  
  return ; 
  
}


template <class Tdata,class Tcalc>
__global__
void cuda_get_crsmatrix_neumanndiag_(
			       PARCEL_INT n,
			       Tdata      *crsA,
			       int *crsRow_ptr,
			       int *crsCol,
			       Tdata    *D,
			       Tdata    *Dinv
			       ){

  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  for(PARCEL_INT i=tid;i<n;i=i+step){
    for(PARCEL_INT j=crsRow_ptr[i];j<crsRow_ptr[i+1];j++){
      if(crsCol[j] == i){
	D[i]    = crsA[j];
	if(D[i]==0.0){
	  D[i]    = 0.0;
	  Dinv[i] = 1.0;
	}else{
	  Dinv[i] = 1.0/D[i];
	}
	crsA[j] = 0.0; // 対角をゼロにする
      }
    }
    for(PARCEL_INT j=crsRow_ptr[i];j<crsRow_ptr[i+1];j++){
      crsA[j] = -crsA[j] * Dinv[i];
    }
  }
  return ;
}


template <class Tdata,class Tcalc>
__global__
void cuda_neumann_preparation_(
			       PARCEL_INT n,
			       Tdata    *x,
			       Tdata    *y,
			       Tdata    *z,
			       Tdata    *v
			       ){

  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  for(PARCEL_INT i=tid;i<n;i=i+step){
    Tdata tmp = x[i] * y[i];
    z[i] = tmp;
    v[i] = tmp;
  }
  return ;

}

template <class Tdata,class Tcalc>
__global__
void cuda_neumann_update_(
			       PARCEL_INT n,
			       Tdata    alpha,
			       Tdata    *x,
			       Tdata    *y,
			       Tdata    *z
			       ){

  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  for(PARCEL_INT i=tid;i<n;i=i+step){
    Tdata tmp = x[i];
    z[i] = z[i] + tmp * alpha; 
    y[i] = tmp;
  }
  return ;

}

template <class Tdata,class Tcalc>
__global__
void cuda_jacobi_update_(
			       PARCEL_INT n,
			       Tdata    *x,
			       Tdata    *b,
			       Tdata    *Ax,
			       Tdata    *Dinv
			       ){
  // s = s + Dinv *( r - y ) 
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  for(PARCEL_INT i=tid;i<n;i=i+step){
    x[i] = x[i] + Dinv[i] * ( b[i] - Ax[i] ); 
  }
  return ;

}

extern "C"{

  void  cuda_get_diamatrix_neumanndiag(
		 PARCEL_INT n,
		 PARCEL_INT nnd,
		 PARCEL_INT diag_offset,
		 devicemem *diaA_, 
		 devicemem *offset_, 
		 devicemem *iLs_, 
		 devicemem *D_   , PARCEL_INT offset_D,
		 devicemem *Dinv_, PARCEL_INT offset_Dinv
				){

    double* diaA = (double*)diaA_[0].a_d; 
    PARCEL_INT* offset = (PARCEL_INT*)offset_[0].a_d; 
    PARCEL_INT* iLs = (PARCEL_INT*)iLs_[0].a_d; 

    double* D    = (double*)D_[0].a_d;
    double* Dinv = (double*)Dinv_[0].a_d;

    cuda_get_diamatrix_neumanndiag_<double,double> <<<4096,512>>>(
							  n,
							  nnd,
							  512,
							  diag_offset,
							  diaA,
							  &D[offset_D],
							  &Dinv[offset_Dinv],
							  // offset,
							  iLs
							  );

    return ;

  }


  void  cuda_get_crsmatrix_neumanndiag(
		 PARCEL_INT n,
		 devicemem *crsA_, 
		 devicemem *crsRow_ptr_, 
		 devicemem *crsCol_, 
		 devicemem *D_   , PARCEL_INT offset_D,
		 devicemem *Dinv_, PARCEL_INT offset_Dinv
				){
    double* crsA           = (double*)crsA_[0].a_d;
    int* crsRow_ptr = (int*)crsRow_ptr_[0].a_d;
    int* crsCol     = (int*)crsCol_[0].a_d;

    double* D    = (double*)D_[0].a_d;
    double* Dinv = (double*)Dinv_[0].a_d;

    cuda_get_crsmatrix_neumanndiag_<double,double> <<<512,512>>>(
							  n,
							  crsA,
							  crsRow_ptr,
							  crsCol,
							  &D[offset_D],
							  &Dinv[offset_Dinv]
							  );

    return ;

  }

  void cuda_neumann_preparation(
		 int numGPUblock,
		 int numGPUthread,
		 PARCEL_INT n,
		 devicemem *x_, PARCEL_INT offset_x,
		 devicemem *y_, PARCEL_INT offset_y,
		 devicemem *z_, PARCEL_INT offset_z,
		 devicemem *v_, PARCEL_INT offset_v
				){
    double* x= (double*)x_[0].a_d;
    double* y= (double*)y_[0].a_d;
    double* z= (double*)z_[0].a_d;
    double* v= (double*)v_[0].a_d;

    cuda_neumann_preparation_<double,double> <<<numGPUblock,numGPUthread>>>(
									    n,
									    &x[offset_x],
									    &y[offset_y],
									    &z[offset_z],
									    &v[offset_v]
									    );

    return ;

  }

  void cuda_neumann_update(
		 int numGPUblock,
		 int numGPUthread,
		 PARCEL_INT n,
		 double alpha,
		 devicemem *x_, PARCEL_INT offset_x,
		 devicemem *y_, PARCEL_INT offset_y,
		 devicemem *z_, PARCEL_INT offset_z
				){

    double* x = (double*)x_[0].a_d;
    double* y = (double*)y_[0].a_d;
    double* z = (double*)z_[0].a_d;

    cuda_neumann_update_<double,double> 
      <<<numGPUblock,numGPUthread>>>(
				     n,
				     alpha,
				     &x[offset_x],
				     &y[offset_y],
				     &z[offset_z]
				     );

    return ;

  }

  void cuda_jacobi_update(
		 int numGPUblock,
		 int numGPUthread,
		 PARCEL_INT n,
		 devicemem *x_,  PARCEL_INT offset_x,
		 devicemem *b_,  PARCEL_INT offset_b,
		 devicemem *Ax_, PARCEL_INT offset_Ax,
		 devicemem *Dinv_,  PARCEL_INT offset_Dinv
				){

    double* x    = (double*)x_[0].a_d;
    double* b    = (double*)b_[0].a_d;
    double* Ax   = (double*)Ax_[0].a_d;
    double* Dinv = (double*)Dinv_[0].a_d;
    cuda_jacobi_update_<double,double> 
      <<<numGPUblock,numGPUthread>>>(
				     n,
				     &x[offset_x],
				     &b[offset_b],
				     &Ax[offset_Ax],
				     &Dinv[offset_Dinv]
				     );

    return ;

  }

}
