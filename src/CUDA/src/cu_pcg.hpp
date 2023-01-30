
#ifndef HEADER_CUPCG
#define HEADER_CUPCG

#include "reduciton_cuda.hpp"

template <class Tdata,class Tcalc>
__global__
void cuda_CG_axpy2_( 
		    PARCEL_INT n,  
		    Tdata      *x  ,
		    Tdata      *p  ,
		    Tdata      *s  ,
		    Tdata      alpha,
		    Tdata      beta
		     ){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  for(PARCEL_INT i=tid;i<n;i=i+step){
    x[i] = x[i] + alpha*p[i];
    p[i] = s[i] + beta*p[i];
  }
  return ;
}

template <class Tdata,class Tcalc>
__global__
void cuda_CG_calc_dot2_( 
			PARCEL_INT n,
			Tdata      *s  ,
			Tdata      *r  ,
			Tdata      *rho1,
			Tdata      *rdot
			 ){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  Tdata tmp1 = 0.0;
  Tdata tmp2 = 0.0;
  for(PARCEL_INT i=tid;i<n;i=i+step){
    tmp1 = tmp1+s[i]*r[i];
    tmp2 = tmp2+r[i]*r[i];
  }
  local_sum_cuda<Tdata,Tcalc>(tmp1,rho1,threadIdx.x);
  local_sum_cuda<Tdata,Tcalc>(tmp2,rdot,threadIdx.x);

  return ;
}

#endif
