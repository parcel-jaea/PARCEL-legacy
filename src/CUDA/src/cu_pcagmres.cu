#include <cuda.h>
#include <cuda_runtime.h>
#include <mpi.h>
#include "cu_common.h"
#include "cu_common.hpp"

template <class Tdata,class Tcalc>
__global__
void cuda_sum_vec_cagmres_( 
		    PARCEL_INT n,
		    PARCEL_INT s,
		    PARCEL_INT m,
		    PARCEL_INT mm,
		    Tdata      *Q, PARCEL_INT ldq,
		    Tdata      *V, PARCEL_INT ldv,
		    Tdata      *Qy, PARCEL_INT ldqy,
		    Tdata      *y
		    ){

  const int block_id    = blockIdx.x;
  const int block_step  = gridDim.x;
  const int thread_id   = threadIdx.x;
  const int thread_step = blockDim.x;
  const int nblock = blockDim.x;

    // !$OMP PARALLEL DO PRIVATE(jj,j,i) 
    // do jj=1,n,nvecBlock 
    //    do j=jj,min(n,jj+nvecBlock-1)
    //       Qy(j)=0.0d0
    //    enddo
    //    do i=1,m-mm
    //       do j=jj,min(n,jj+nvecBlock-1)
    //          Qy(j)=Qy(j)+y(i)*Q(j+n*(i-1))
    //       enddo
    //    enddo
    //    do i=1,mm
    //       do j=jj,min(n,jj+nvecBlock-1)
    //          Qy(j)=Qy(j)+y(i+m-mm)*V(j+n*(i-1))
    //       enddo
    //    enddo
    // enddo
    // !$OMP END PARALLEL DO 

#if 1

  for(PARCEL_INT jj=1+block_id*nblock;jj<=n;jj=jj+block_step*nblock){
    PARCEL_INT jmin = n > (jj+nblock-1) ? (jj+nblock-1) : n;
    for(PARCEL_INT j=jj + thread_id ;j<=jmin;j=j + thread_step){
      Qy[j-1]=0.0;
    }
    for(PARCEL_INT i=1;i<=m-mm;i++){
      for(PARCEL_INT j=jj + thread_id ;j<=jmin;j=j + thread_step){
	Qy[j-1]=Qy[j-1]+y[i-1]*Q[j+ldq*(i-1)-1];
      }
    }
    for(PARCEL_INT i=1;i<=mm;i++){
      for(PARCEL_INT j=jj + thread_id ;j<=jmin;j=j + thread_step){
	Qy[j-1]=Qy[j-1]+y[i+m-mm-1]*V[j+ldv*(i-1)-1];
      }
    }
    
  }


#else
  for(PARCEL_INT jj=1+block_id*nblock;jj<=n;jj=jj+block_step*nblock){
    PARCEL_INT jmin = (jj+nblock-1) < n ? (jj+nblock-1) : n ;
    for(PARCEL_INT j=jj;j<=jmin;j++){
      Qy[j]=0.0;
    }
    for(PARCEL_INT i=1;i<=(m-mm);i++){
      PARCEL_INT jmin = (jj+nblock-1) < n ? (jj+nblock-1) : n ;
      for(PARCEL_INT j=jj + thread_id ;j<=jmin;j=j + thread_step){
	Qy[j]=Qy[j]+y[i]*Q[j+ldq*(i-1)];
      }
    }
    for(PARCEL_INT i=1;i<=mm;i++){
      PARCEL_INT jmin = (jj+nblock-1) < n ? (jj+nblock-1) : n ;
      for(PARCEL_INT j=jj + thread_id ;j<=jmin;j=j + thread_step){
	Qy[j]=Qy[j]+y[i+m-mm]*V[j+ldv*(i-1)];
      }
    }
  }
#endif

  return ;
}

template <class Tdata,class Tcalc>
__global__
void cuda_sum_vec_( 
		    PARCEL_INT n,  
		    PARCEL_INT m,  
		    PARCEL_INT vecLen,  
		    int nblock,
		    Tdata      *v ,
		    Tdata      *v0,
		    Tdata      *y
		    ){
  // const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  // const int step = blockDim.x*gridDim.x;
  const int block_id    = blockIdx.x;
  const int block_step  = gridDim.x;
  const int thread_id   = threadIdx.x;
  const int thread_step = blockDim.x;

  for(PARCEL_INT jj=1+block_id*nblock;jj<=n;jj=jj+block_step*nblock){
    PARCEL_INT jjmin = (jj+nblock-1) < n ? (jj+nblock-1) : n ;
    //    for(PARCEL_INT j=jj;j<=jjmin;j++){
    for(PARCEL_INT j=jj + thread_id ;j<=jjmin;j=j + thread_step){
      v0[j-1]=0.0;
    }
    for(PARCEL_INT i=1;i<=m;i++){
      for(PARCEL_INT j=jj + thread_id ;j<=jjmin;j=j + thread_step){
      //      for(PARCEL_INT j=jj;j<=jjmin;j++){
	//	v0[j-1]=v0[j-1]+y[i-1]*v[j+n*(i-1)-1];
	v0[j-1]=v0[j-1]+y[i-1]*v[j+vecLen*(i-1)-1];
      }
    }
  }

  return ;

}


extern "C" {
  void cuda_sum_vec(
		    int numGPUblock,
		    int numGPUthread,
		    PARCEL_INT n,  
		    PARCEL_INT m,  
		    PARCEL_INT vecLen,  
		    devicemem *v_, PARCEL_INT offset_v,
		    devicemem *v0_, PARCEL_INT offset_v0,
		    devicemem *y_
		    ){

    double *v  = (double *)v_[0].a_d;
    double *v0 = (double *)v0_[0].a_d;
    double *y  = (double *)y_[0].a_d;
    
    cuda_sum_vec_
      <double,double> <<<numGPUblock,numGPUthread>>>(
						     n,
						     m,
						     vecLen,  
						     numGPUthread,
						     &v[offset_v] ,
						     &v0[offset_v0],
						     y
						     );

    return ;
  }

  void cuda_sum_vec_cagmres(
		    int numGPUblock,
		    int numGPUthread,
		    PARCEL_INT n,
		    PARCEL_INT s,
		    PARCEL_INT m,
		    PARCEL_INT mm,
		    devicemem *Q_ , PARCEL_INT offset_Q , PARCEL_INT ldq,
		    devicemem *V_ , PARCEL_INT offset_V , PARCEL_INT ldv,
		    devicemem *Qy_, PARCEL_INT offset_Qy, PARCEL_INT ldqy,
		    devicemem *y_ , PARCEL_INT offset_y 
		    ){

    double *Q  = (double *)Q_[0].a_d;
    double *V  = (double *)V_[0].a_d;
    double *Qy = (double *)Qy_[0].a_d;
    double *y  = (double *)y_[0].a_d;
    
    cuda_sum_vec_cagmres_
      <double,double> <<<numGPUblock,numGPUthread>>>(
      //      <double,double> <<<1,1>>>(
						     n,
						     s,
						     m,
						     mm,
						     &Q[offset_Q], ldq,
						     &V[offset_V], ldv,
						     &Qy[offset_Qy], ldqy,
						     &y[offset_y] 
						     );

    return ;
  }

}
