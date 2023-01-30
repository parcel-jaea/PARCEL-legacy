#include <cuda.h>
#include <cuda_runtime.h>
#include <mpi.h>
#include "cu_common.h"
#include "cu_common.hpp"

#include "cu_pcg.hpp"

extern "C" {
  void cuda_CG_axpy2_FP64(
			  int numGPUblock,
			  int numGPUthread,
			  PARCEL_INT n,  
			  devicemem *x_, PARCEL_INT offset_x,
			  devicemem *p_, PARCEL_INT offset_p,
			  devicemem *s_, PARCEL_INT offset_s,
			  double      alpha,
			  double      beta
			  ){
    double *x = (double *)x_[0].a_d;
    double *p = (double *)p_[0].a_d;
    double *s = (double *)s_[0].a_d;
    cuda_CG_axpy2_<double,double> <<<numGPUblock,numGPUthread>>>(
							    n,
							    &x[offset_x],
							    &p[offset_p],
							    &s[offset_s],
							    alpha,
							    beta
							    );
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
			      ){
    cuda_set_vec_val_<double,double> <<<1,32>>>(1,0.0,(double*)rho1[0].a_d);
    cuda_set_vec_val_<double,double> <<<1,32>>>(1,0.0,(double*)rdot[0].a_d);

    double *s = (double *)s_[0].a_d;
    double *r = (double *)r_[0].a_d;
    cuda_CG_calc_dot2_<double,double> <<<numGPUblock,numGPUthread>>>(
								n,
								&s[offset_s],
								&r[offset_r],
								(double *)rho1[0].a_d,
								(double *)rdot[0].a_d  
								);

    double sbuf[2];
    double rbuf[2];
    MPI_Comm comm = MPI_Comm_f2c(f_comm);
    const size_t commSize = sizeof(double);
    cudaMemcpy(&sbuf[0],rho1[0].a_d,commSize,cudaMemcpyDeviceToHost);
    cudaMemcpy(&sbuf[1],rdot[0].a_d,commSize,cudaMemcpyDeviceToHost);
  
    cudaDeviceSynchronize();
    MPI_Allreduce(sbuf,rbuf,2,MPI_DOUBLE,MPI_SUM,comm);
  
    cudaMemcpy(rho1[0].a_d,&rbuf[0],commSize,cudaMemcpyHostToDevice);
    cudaMemcpy(rdot[0].a_d,&rbuf[1],commSize,cudaMemcpyHostToDevice);
  
  }

}

