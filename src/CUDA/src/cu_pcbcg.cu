#include <cuda.h>
#include <cuda_runtime.h>
#include <mpi.h>
#include "cu_common.h"
#include "cu_common.hpp"

   // x  = x  + Q  * alpha
   // Q  = S  - Q  * beta --> Sをwork領域に利用可能
   // AQ = AS - AQ * beta --> ASをwork領域に利用可能
   // Qtr    = <Q|r>
   // <r|r>  
   // GPU版ではS,ASに出力してポインタをswapする.

// QtAQ, Qtr, Qtr はゼロに初期化済み
template <class Tdata,class Tcalc>
__global__
void CBCG_gemm_dot_local_(
		    PARCEL_INT n,
		    PARCEL_INT k,
		    Tdata      *x,  PARCEL_INT ldx,
		    Tdata      *Q,  PARCEL_INT ldq,
		    Tdata      *S,  PARCEL_INT lds,
		    Tdata      *AQ, PARCEL_INT ldaq,
		    Tdata      *AS, PARCEL_INT ldas,
		    Tdata      *r,  PARCEL_INT ldr,
		    Tdata      *alpha_,
		    Tdata      *beta_,
		    Tdata      *QAQ,
		    Tdata      *Qr,
		    Tdata      *rnorm
		      ){

  // k個のグループ : k*128,k*256,..... : ブロック総数
  const PARCEL_INT thread_id   = threadIdx.x;
  const PARCEL_INT thread_step = blockDim.x;
  const PARCEL_INT block_id    = blockIdx.x;
  const PARCEL_INT block_step  = gridDim.x;
  const PARCEL_INT jblock      = thread_step * block_step;
  const PARCEL_INT jj0         = block_id*thread_step;

#if 0
  extern __shared__ Tdata shared_mem[];
  Tdata *alpha = &shared_mem[0];
  Tdata *beta  = &shared_mem[k];
  for(PARCEL_INT l=thread_id;l<k;l=l+thread_step){
    alpha[l] = alpha_[l];
  }
  for(PARCEL_INT l=thread_id;l<k*k;l=l+thread_step){
    beta[l] = beta_[l];
  }
  Tdata *Qr_ = &shared_mem[k+k*k+thread_id*k];  
  for(PARCEL_INT l=0;l<k;l++){
    Qr_[l] = 0.0;
  }

  __syncthreads();
#else
  Tdata *alpha = alpha_;
  Tdata *beta  = beta_;
  Tdata Qr_[4];  
  for(PARCEL_INT l=0;l<k;l++){
    Qr_[l] = 0.0;
  }
#endif


  Tdata rnorm_ = 0.0;
  for(PARCEL_INT jj=jj0;jj<n;jj=jj+jblock){
    const PARCEL_INT jjmin = (jj+thread_step) < n ? (jj+thread_step) : n ;
    const PARCEL_INT j0 = thread_id;
    for(PARCEL_INT j=j0+jj;j<jjmin;j=j+thread_step){

      Tdata r_ = r[j];
      Tdata x_ = 0.0;

      for(PARCEL_INT l=0;l<k;l++){
	// Tdata alpha_ = alpha[l];
	x_ = x_ + Q[j+l*ldq] * alpha[l];
      }
      x[j] = x[j] + x_;

      rnorm_ = rnorm_ +r_*r_;      

      // for(int l=0;l<k;l=l+2){
      for(PARCEL_INT l=0;l<k;l++){
	Tdata S_    =  S[j+l*lds];
	Tdata AS_   = AS[j+l*ldas];
        #pragma unroll 4
	for(PARCEL_INT i=0;i<k;i++){
	  // S  = S  - Q  * beta
	  // AS = AS - AQ * beta
	  Tdata bata0 = beta[i+l*k];
	  S_  = S_  -  Q[j+i*ldq]  * bata0;
	  AS_ = AS_ - AQ[j+i*ldaq] * bata0;
	}
	S[j+l*lds]   = S_;
	AS[j+l*ldas] = AS_;
	Qr_[l] = Qr_[l] + S_*r_;
      }

    }
  }

  for(PARCEL_INT l=0;l<k;l++){
    local_sum_cuda<Tdata,Tcalc>(Qr_[l],&Qr[l],threadIdx.x);
  }
  local_sum_cuda<Tdata,Tcalc>(rnorm_,rnorm,threadIdx.x);


  return ;

}

template <class Tdata,class Tcalc>
__global__
void CBCG_gemm_dot_local_ddm_(
			      PARCEL_INT nx,PARCEL_INT ny,PARCEL_INT nz,
			      PARCEL_INT stmx,PARCEL_INT stmy,PARCEL_INT stmz,
			      PARCEL_INT mx,PARCEL_INT mxy,
			      PARCEL_INT k,
			      Tdata      *x,  PARCEL_INT ldx,
			      Tdata      *Q,  PARCEL_INT ldq,
			      Tdata      *S,  PARCEL_INT lds,
			      Tdata      *AQ, PARCEL_INT ldaq,
			      Tdata      *AS, PARCEL_INT ldas,
			      Tdata      *r,  PARCEL_INT ldr,
			      Tdata      *alpha_,
			      Tdata      *beta_,
			      Tdata      *QAQ,
			      Tdata      *Qr,
			      Tdata      *rnorm
			      ){

  extern __shared__ Tdata shared_mem[];
  Tdata *alpha = &shared_mem[0];
  Tdata *beta  = &shared_mem[k];

  const int thread_id   = threadIdx.x;
  const int thread_step = blockDim.x;
  for(PARCEL_INT l=thread_id;l<k;l=l+thread_step){
    alpha[l] = alpha_[l];
  }
#if 1
  for(PARCEL_INT l=thread_id;l<k*k;l=l+thread_step){
    beta[l] = beta_[l];
  }
#else
  for(PARCEL_INT l=0;l<k;l++){
    for(PARCEL_INT i=thread_id;i<k;i=i+thread_step){
      beta[l+i*k] = beta_[i+l*k];
    }
  }
#endif

  Tdata *Qr_ = &shared_mem[k+k*k+thread_id*k];  
  for(PARCEL_INT l=0;l<k;l++){
    Qr_[l] = 0.0;
  }
  __syncthreads();

  Tdata rnorm_ = 0.0;
  for(
	PARCEL_INT jz = threadIdx.z + blockIdx.z*blockDim.z;
	jz<nz ;
	jz = jz + blockDim.z*gridDim.z ){
    const PARCEL_INT jz0 = mxy *( jz +stmz );
    for(
	PARCEL_INT jy = threadIdx.y + blockIdx.y*blockDim.y;
	jy<ny ;
	jy = jy + blockDim.y*gridDim.y ){
      const PARCEL_INT jy0 = mx *( jy +stmy );
      const PARCEL_INT jx0 = stmx  + jy0  + jz0;
      for(
	  PARCEL_INT jx = threadIdx.x + blockIdx.x*blockDim.x;
	  jx<nx ;
	  jx = jx + blockDim.x*gridDim.x ){
	const PARCEL_INT j = jx + jx0;

	Tdata x_ = 0.0;
#pragma unroll 4
	for(PARCEL_INT l=0;l<k;l++){
	  x_ = x_ + Q[j+l*ldq] * alpha[l];
	 }
	 x[j] = x[j] + x_;
	
	Tdata r_ = r[j];
	rnorm_ = rnorm_ +r_*r_;      
	for(PARCEL_INT l=0;l<k;l++){
	  Tdata S_    =  S[j+l*lds];
	  Tdata AS_   = AS[j+l*ldas];
#pragma unroll 4
	    for(PARCEL_INT i=0;i<k;i++){
	      // S  = S  - Q  * beta
	      // AS = AS - AQ * beta
	      Tdata bata0 = beta[i+l*k];
	      //	      Tdata bata0 = beta[l+i*k];
	      S_  = S_  -  Q[j+i*ldq]  * bata0;
	      AS_ = AS_ - AQ[j+i*ldaq] * bata0;
	    }
	    S[j+l*lds]   = S_;
	    AS[j+l*ldas] = AS_;
	    Qr_[l] = Qr_[l] + S_*r_;
	}

      }
    }
  }

  for(PARCEL_INT l=0;l<k;l++){
    local_sum_cuda<Tdata,Tcalc>(Qr_[l],&Qr[l],threadIdx.x);
  }
  local_sum_cuda<Tdata,Tcalc>(rnorm_,rnorm,threadIdx.x);

  return ;

}

// z = x + y*beta
template <class Tdata,class Tcalc>
__global__
void gemm_(
		     Tdata *x  ,PARCEL_INT ldx,
		     Tdata *y  ,PARCEL_INT ldy,
		     // Tdata *z  ,PARCEL_INT ldz,
		     Tdata *beta,PARCEL_INT ldbeta,
		     PARCEL_INT n,
		     PARCEL_INT s,PARCEL_INT t
		     ){

  // 128*s,256*s,512*s,..... @ ブロック総数
  const int thread_id   = threadIdx.x;
  const int thread_step = blockDim.x;

  const int block_step  = gridDim.x/(s);  
  const int kid = blockIdx.x/(block_step);

  const int block_id = blockIdx.x%(block_step);
  const int iid      = thread_id + block_id*thread_step;

  for(PARCEL_INT i = iid;i<n;i=i+thread_step*block_step){
    Tdata y_ = y[i+kid*ldx];
    //    for(PARCEL_INT lid = 0;lid<s;lid++){
    for(PARCEL_INT lid = 0;lid<t;lid++){
      y_ = y_ + x[i+lid*ldx]*beta[lid+kid*ldbeta];
    }
    y[i+kid*ldx] = y_;
  }

  return ;
}

template <class Tdata,class Tcalc>
__global__
void CBCG_axpy2_(
	       PARCEL_INT n,
	       Tdata alpha,
	       Tdata beta,
	       Tdata    *x,
	       Tdata    *y,
	       Tdata    *z
	       ){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  for(PARCEL_INT i=tid;i<n;i=i+step){
    z[i] = alpha * x[i] + beta * y[i];
  }
  return ;
}


extern "C" {
  void cuda_gemm_dot_local(
			 devicemem *W_,PARCEL_INT offset_W,PARCEL_INT ldw,
			 devicemem *V_,PARCEL_INT offset_V,PARCEL_INT ldv,
			 devicemem *wv_,PARCEL_INT offset_wv,
			 PARCEL_INT n,
			 PARCEL_INT wrow,
			 PARCEL_INT vrow
			   );

  // void swap_ptr(double **x,double **y){
  //   double *tmp;
  //   tmp=*y;
  //   *y=*x;
  //   *x=tmp;
  // }
  void swap_ptr(void **x,void **y){
    void *tmp;
    tmp=*y;
    *y=*x;
    *x=tmp;
  }

void cuda_gemm(
		     devicemem *x_,PARCEL_INT offset_x,PARCEL_INT ldx,
		     devicemem *y_,PARCEL_INT offset_y,PARCEL_INT ldy,
		     // devicemem *z_,PARCEL_INT offset_z,PARCEL_INT ldz,
		     devicemem *beta_,PARCEL_INT offset_beta,PARCEL_INT ldbeta,
		     PARCEL_INT n,
		     PARCEL_INT s,PARCEL_INT t
		     ){

    double* x = (double*)x_[0].a_d;
    double* y = (double*)y_[0].a_d;
    // double* z = (double*)z_[0].a_d;
    double* beta = (double*)beta_[0].a_d;
    int numGPUblock  = 512*s;
    while(numGPUblock>16384){
      numGPUblock = numGPUblock /2;
    }
    const int numGPUthread = 512;

    // gemm
    gemm_<double,double> 
      <<<numGPUblock,numGPUthread>>>(
				     &x[offset_x],ldx,
				     &y[offset_y],ldy,
				     // &z[offset_z],ldz,
				     &beta[offset_beta],ldbeta,
				     n,
				     s, t
				     );

    return ;
  }

  void cuda_CBCG_gemm_dot_local(
			      PARCEL_INT n,
			      PARCEL_INT k,
			      devicemem *x_ , PARCEL_INT offset_x , PARCEL_INT ldx ,
			      devicemem *Q_ , PARCEL_INT offset_Q , PARCEL_INT ldq ,
			      devicemem *S_ , PARCEL_INT offset_S , PARCEL_INT lds ,
			      devicemem *AQ_, PARCEL_INT offset_AQ, PARCEL_INT ldaq,
			      devicemem *AS_, PARCEL_INT offset_AS, PARCEL_INT ldas,
			      devicemem *r_ , PARCEL_INT offset_r , PARCEL_INT ldr ,
			      devicemem *alpha_,
			      devicemem *beta_,
			      devicemem *QAQ_,
			      devicemem *Qr_,
			      devicemem *rnorm_
			      ){
  double* x     = (double*)x_[0].a_d;
  double* Q     = (double*)Q_[0].a_d;
  double* S     = (double*)S_[0].a_d;
  double* AQ    = (double*)AQ_[0].a_d;
  double* AS    = (double*)AS_[0].a_d;
  double* r     = (double*)r_[0].a_d;
  double* alpha = (double*)alpha_[0].a_d;
  double* beta  = (double*)beta_[0].a_d;
  double* QAQ   = (double*)QAQ_[0].a_d;
  double* Qr    = (double*)Qr_[0].a_d;
  double* rnorm = (double*)rnorm_[0].a_d;

#if 0
  int numGPUblock  = 512*k;
  while(numGPUblock>16384){
    numGPUblock = numGPUblock /2;
  }
  const int numGPUthread = 512;
#else
  const int numGPUblock  = 512;
  const int numGPUthread = 512;
#endif
  //  numGPUblock  =k;

  cuda_set_vec_val_<double,double> <<<1,32>>>(1  ,0.0,rnorm);
  cuda_set_vec_val_<double,double> <<<1,32>>>(k  ,0.0,Qr);
  cuda_set_vec_val_<double,double> <<<1,32>>>(k*k,0.0,QAQ);

  cudaDeviceSynchronize();
  // (k*threads)
  const size_t shared_mem= ( k*numGPUthread  + (k*k)+k)*sizeof(double);
  CBCG_gemm_dot_local_<double,double> 
    <<<numGPUblock,numGPUthread,shared_mem>>> (
					       n,
					       k,
					       &x[offset_x],  ldx,
					       &Q[offset_Q],  ldq,
					       &S[offset_S],  lds,
					       &AQ[offset_AQ], ldaq,
					       &AS[offset_AS], ldas,
					       &r[offset_r],  ldr,
					       alpha,
					       beta,
					       QAQ,
					       Qr,
					       rnorm
					       );
  cudaDeviceSynchronize();

   swap_ptr( &(Q_[0].a_d) ,&(S_[0].a_d)  );
   swap_ptr( &(AQ_[0].a_d),&(AS_[0].a_d) );
    
   cudaDeviceSynchronize();

  return ;
}

  void cuda_CBCG_gemm_dot_local_ddm(
		    PARCEL_INT nx,PARCEL_INT ny,PARCEL_INT nz,
		    PARCEL_INT stmx,PARCEL_INT stmy,PARCEL_INT stmz,
		    PARCEL_INT mx,PARCEL_INT mxy,
			      PARCEL_INT k,
			      devicemem *x_ , PARCEL_INT offset_x , PARCEL_INT ldx ,
			      devicemem *Q_ , PARCEL_INT offset_Q , PARCEL_INT ldq ,
			      devicemem *S_ , PARCEL_INT offset_S , PARCEL_INT lds ,
			      devicemem *AQ_, PARCEL_INT offset_AQ, PARCEL_INT ldaq,
			      devicemem *AS_, PARCEL_INT offset_AS, PARCEL_INT ldas,
			      devicemem *r_ , PARCEL_INT offset_r , PARCEL_INT ldr ,
			      devicemem *alpha_,
			      devicemem *beta_,
			      devicemem *QAQ_,
			      devicemem *Qr_,
			      devicemem *rnorm_
			      ){
  double* x     = (double*)x_[0].a_d;
  double* Q     = (double*)Q_[0].a_d;
  double* S     = (double*)S_[0].a_d;
  double* AQ    = (double*)AQ_[0].a_d;
  double* AS    = (double*)AS_[0].a_d;
  double* r     = (double*)r_[0].a_d;
  double* alpha = (double*)alpha_[0].a_d;
  double* beta  = (double*)beta_[0].a_d;
  double* QAQ   = (double*)QAQ_[0].a_d;
  double* Qr    = (double*)Qr_[0].a_d;
  double* rnorm = (double*)rnorm_[0].a_d;

  cuda_set_vec_val_<double,double> <<<1,32>>>(1  ,0.0,rnorm);
  cuda_set_vec_val_<double,double> <<<1,32>>>(k  ,0.0,Qr);
  cuda_set_vec_val_<double,double> <<<1,32>>>(k*k,0.0,QAQ);

  int nx_para = nx < 1024 ? nx : 512;
  int ny_para = ny < 1024 ? ny : 512;
  int nz_para = nz < 1024 ? nz : 512;

  cudaDeviceSynchronize();
#if 0
  dim3 block(nx,1,1); 
  dim3 grid(1,ny,nz); 
#else
  dim3 block(nx_para,1,1); 
  dim3 grid(1,ny_para,nz_para); 
#endif
  const size_t shared_mem= ( k*nx + (k*k)+k)*sizeof(double);
  CBCG_gemm_dot_local_ddm_<double,double> 
    <<<grid,block,shared_mem>>>(
				   nx,ny,nz,
				   stmx,stmy,stmz,
				   mx, mxy,
				   k,
				   &x[offset_x],  ldx,
				   &Q[offset_Q],  ldq,
				   &S[offset_S],  lds,
				   &AQ[offset_AQ], ldaq,
				   &AS[offset_AS], ldas,
				   &r[offset_r],  ldr,
				   alpha,
				   beta,
				   QAQ,
				   Qr,
				   rnorm
				   );
   cudaDeviceSynchronize();

   swap_ptr( &(Q_[0].a_d) ,&(S_[0].a_d)  );
   swap_ptr( &(AQ_[0].a_d),&(AS_[0].a_d) );
    
   cudaDeviceSynchronize();

  return ;
}

  void cuda_CBCG_axpy2(
		       int numGPUblock,
		       int numGPUthread,
		       PARCEL_INT n,
		       double alpha,
		       double beta,
		       devicemem *x_,PARCEL_INT offset_x,
		       devicemem *y_,PARCEL_INT offset_y,
		       devicemem *z_,PARCEL_INT offset_z
		       ){
    double *x = (double *)x_[0].a_d;
    double *y = (double *)y_[0].a_d;
    double *z = (double *)z_[0].a_d;
    CBCG_axpy2_<double,double> <<<numGPUblock,numGPUthread>>>(
							     n,
							     alpha,
							     beta,
							     &x[offset_x],
							     &y[offset_y],
							     &z[offset_z]
							     );

    return ;
  }


}
