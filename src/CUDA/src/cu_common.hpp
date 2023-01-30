#ifndef HEADER_CUCOMMMON
#define HEADER_CUCOMMMON

#include <cuda_runtime.h>
#include <cuda_bf16.h>
#include "reduciton_cuda.hpp"

template <class Tdata,class Tcalc>
__global__
void cuda_vec_scaling_(
	       PARCEL_INT n,
	       Tdata    *x,
	       Tdata    *y,
	       Tdata    *z
	       ){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  for(PARCEL_INT i=tid;i<n;i=i+step){
    z[i] = x[i] * y[i];
  }
  return ;
}


template <class Tdata,class Tcalc>
__global__
void cuda_axpy_(
	       PARCEL_INT n,
	       Tdata alpha,
	       Tdata beta,
	       Tdata    *x,
	       Tdata    *y
	       ){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  for(PARCEL_INT i=tid;i<n;i=i+step){
    y[i] = alpha * x[i] + beta * y[i];
  }
  return ;
}

template <class Tdata,class Tcalc>
__global__
void cuda_zero_init_(
	       PARCEL_INT n,
	       Tdata    *x
	       ){
  //  return ;
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  for(PARCEL_INT i=tid;i<n;i=i+step){
    x[i] = 0.0;
    //        printf("i =%4d \n",i);
  }
  return ;
}

template <class Tdata,class Tcalc>
__global__
void cuda_axpy_vec3_(
	       PARCEL_INT n,
	       Tdata alpha,
	       Tdata beta,
	       Tdata gamma,
	       Tdata    *x,
	       Tdata    *y,
	       Tdata    *z
	       ){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  for(PARCEL_INT i=tid;i<n;i=i+step){
    z[i] = alpha * x[i] + beta * y[i] + gamma*z[i];
  }
  return ;
}

template <class Tdata,class Tcalc>
__global__
void cuda_axpy_vec4_(
	       PARCEL_INT n,
	       Tdata alpha,
	       Tdata beta,
	       Tdata gamma,
	       Tdata delta,
	       Tdata    *x,
	       Tdata    *y,
	       Tdata    *z,
	       Tdata    *w
	       ){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  for(PARCEL_INT i=tid;i<n;i=i+step){
    w[i] = alpha*x[i] + beta*y[i] + gamma*z[i] + delta*w[i];
  }
  return ;
}

// y[1:n] = x[1:n,n+1:n*2,....] * alpha[1:kmax]
template <class Tdata,class Tcalc>
__global__
void cuda_gemv_(
		Tdata *x, PARCEL_INT ldx,
		Tdata *y, PARCEL_INT ldy,
		Tdata *alpha,
		PARCEL_INT n,
		PARCEL_INT kmax
){

  const PARCEL_INT thread_id = threadIdx.x;
  const PARCEL_INT block_id  = blockIdx.x;
  const PARCEL_INT block_step  = gridDim.x;
  const PARCEL_INT thread_step = blockDim.x;

  const PARCEL_INT nblock = blockDim.x;
  for(PARCEL_INT jj=1+block_id*nblock;jj<=n;jj=jj+block_step*nblock){
    PARCEL_INT jjmin = (jj+nblock-1) < n ? (jj+nblock-1) : n ;
    for(PARCEL_INT k=0;k<kmax;k++){
      Tdata alpha_ = alpha[k];
      for(PARCEL_INT j=jj+thread_id-1;j<=jjmin-1;j=j + thread_step){
	y[j] = y[j] + alpha_ * x[j+k*ldx];
      }
    }
  }

  return ;

}


template <class Tdata,class Tcalc>
__global__
void cuda_dot_local_(
	       PARCEL_INT n,
	       Tdata    *val,
	       Tdata    *x,
	       Tdata    *y
	       ){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  Tdata val_ = 0.0;
  //  printf(" tid = %06d, n = %06d , step = %06d \n",tid,n,step);
  for(PARCEL_INT i=tid;i<n;i=i+step){
    val_ = val_ + y[i] * x[i];
    //    printf(" x[%06d] = %f \n",i,x[i]);
  }
  local_sum_cuda<Tdata,Tcalc>(val_,val,threadIdx.x);
  //  printf("val = %f \n",val[0]);
  return ;
}

template <class Tdata,class Tcalc>
__global__
void cuda_set_vec_val_(
	       PARCEL_INT n,
	       Tdata val,
	       Tdata    *x
	       ){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  for(PARCEL_INT i=tid;i<n;i=i+step){
    x[i] = val;
  }
  return ;
}

template <class Tdata,class Tcalc>
__global__
void cuda_add_val_(
		   Tdata val,
		   Tdata *x){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  for(PARCEL_INT i=tid;i<1;i=i+step){
    x[i] = x[i] + val;
  }
  return ;
}
template <class Tdata,class Tcalc>
__global__
void cuda_sub_val_(
		   Tdata val,
		   Tdata *x){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  for(PARCEL_INT i=tid;i<1;i=i+step){
    x[i] = x[i] - val;
  }

  return ;
}

template <class Tdata,class Tcalc>
__global__
void cuda_mul_val_(
		   Tdata val,
		   Tdata *x){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  for(PARCEL_INT i=tid;i<1;i=i+step){
    x[i] = x[i] * val;
  }

  return ;
}

template <class Tdata,class Tcalc>
__global__
void cuda_div_val_(
		   Tdata val,
		   Tdata *x){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  for(PARCEL_INT i=tid;i<1;i=i+step){
    x[i] = x[i] / val;
  }

  return ;
}

template <class Tdata,class Tcalc>
__global__
void cuda_sqrt_val_(
		   Tdata *x){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  for(PARCEL_INT i=tid;i<1;i=i+step){
    x[i] = sqrt(x[i]) ;
  }

  return ;
}

#endif
