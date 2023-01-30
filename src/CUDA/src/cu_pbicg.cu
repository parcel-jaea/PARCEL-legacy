#include <cuda.h>
#include <cuda_runtime.h>
#include <mpi.h>
#include "cu_common.h"
#include "cu_common.hpp"

template <class Tdata,class Tcalc>
__global__
void cuda_BiCG_axpy3_dot_( 
			PARCEL_INT n,
			Tdata      *r,
			Tdata      *e,
			Tdata      *v,
			Tdata      *r0,
			Tdata      c3,
			Tdata      *c1,
			Tdata      *rdot
			 ){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  Tdata tmp1 = 0.0;
  Tdata tmp2 = 0.0;
  for(PARCEL_INT i=tid;i<n;i=i+step){
    r[i] = e[i] - c3*v[i];
    tmp1 = tmp1+r0[i]*r[i];
    tmp2 = tmp2+r[i] *r[i];
  }
  local_sum_cuda<Tdata,Tcalc>(tmp1,  c1,threadIdx.x);
  local_sum_cuda<Tdata,Tcalc>(tmp2,rdot,threadIdx.x);

  return ;
}

extern "C" {
  // r=e-c3*v  
  // c1  =<r0,r>   
  // rdot=<r,r>
  void cuda_BiCG_axpy3_dot_FP64(
				int numGPUblock,int numGPUthread,
				PARCEL_INT n,MPI_Fint f_comm,
				devicemem  *r_, PARCEL_INT offset_r,
				devicemem  *e_, PARCEL_INT offset_e,
				devicemem  *v_, PARCEL_INT offset_v,
				devicemem *r0_, PARCEL_INT offset_r0,
				double c3,
				devicemem *c1,
				devicemem *rdot
				){
    cuda_set_vec_val_<double,double> <<<1,32>>>(1,0.0,(double*)c1[0].a_d);
    cuda_set_vec_val_<double,double> <<<1,32>>>(1,0.0,(double*)rdot[0].a_d);

    double *r = (double *)r_[0].a_d;
    double *e = (double *)e_[0].a_d;
    double *v = (double *)v_[0].a_d;
    double *r0 = (double *)r0_[0].a_d;
    cuda_BiCG_axpy3_dot_<double,double> <<<numGPUblock,numGPUthread>>>(
								       n,
								       &r[offset_r],
								       &e[offset_e],
								       &v[offset_v],
								       &r0[offset_r0],
								       c3,
								       (double*)c1[0].a_d,
								       (double*)rdot[0].a_d
								       );

    double sbuf[2];
    double rbuf[2];
    MPI_Comm comm = MPI_Comm_f2c(f_comm);
    const size_t commSize = sizeof(double);
    cudaMemcpy(&sbuf[0],c1[0].a_d,commSize,cudaMemcpyDeviceToHost);
    cudaMemcpy(&sbuf[1],rdot[0].a_d,commSize,cudaMemcpyDeviceToHost);
  
    cudaDeviceSynchronize();
    MPI_Allreduce(sbuf,rbuf,2,MPI_DOUBLE,MPI_SUM,comm);
  
    cudaMemcpy(c1[0].a_d,&rbuf[0],commSize,cudaMemcpyHostToDevice);
    cudaMemcpy(rdot[0].a_d,&rbuf[1],commSize,cudaMemcpyHostToDevice);

    return ;
  }

  void cuda_CG_calc_dot2_FP64(
			      MPI_Fint f_comm,
			      int numGPUblock,
			      int numGPUthread,
			      PARCEL_INT n,  
			      devicemem   *s_, PARCEL_INT offset_s,
			      devicemem   *r_, PARCEL_INT offset_r,
			      devicemem   *rho1  ,
			      devicemem   *rdot  
			      );

  // ev_dot = <e,v>, vv_dot=<v,v> 
  void cuda_BiCG_calc_dot2_FP64(
				int numGPUblock,
				int numGPUthread,
				PARCEL_INT n, 
				MPI_Fint f_comm, 
				devicemem *e_, PARCEL_INT offset_e,
				devicemem *v_, PARCEL_INT offset_v,
				devicemem *ev_dot, 
				devicemem *vv_dot
				){
    cuda_CG_calc_dot2_FP64(
			   f_comm,
			   numGPUblock,
			   numGPUthread,
			   n,  
			   e_, offset_e,
			   v_, offset_v,
			   ev_dot  ,
			   vv_dot 
			   );
    return ;
  }

}
