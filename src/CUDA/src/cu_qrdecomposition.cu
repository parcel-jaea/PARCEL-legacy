#include <cuda.h>
#include <cuda_runtime.h>

#include <mpi.h>

#include "cu_common.h"
#include "reduciton_cuda.hpp"
#include "cu_common.hpp"

// wv = W * V
template <class Tdata,class Tcalc>
__global__
void gemm_dot_local_(
		     Tdata *W  ,PARCEL_INT ldw,
		     Tdata *V  ,PARCEL_INT ldv,
		     Tdata *wv ,
		     PARCEL_INT n,
		     PARCEL_INT wrow,
		     PARCEL_INT vrow
		     ){

  // 128*wrow*vrow,256*wrow*vrow,512*wrow*vrow,..... @ ブロック総数
  const PARCEL_INT thread_id   = threadIdx.x;
  const PARCEL_INT thread_step = blockDim.x;

  const PARCEL_INT block_step  = gridDim.x/(wrow*vrow);  
  const PARCEL_INT wvid = blockIdx.x/(block_step);
  const PARCEL_INT wid  = wvid%wrow;
  const PARCEL_INT vid  = wvid/wrow;

  const PARCEL_INT block_id = blockIdx.x%(block_step);
  const PARCEL_INT iid = thread_id + block_id*thread_step;

  Tdata wv_ = 0.0;
  for(PARCEL_INT i = iid;i<n;i=i+thread_step*block_step){
    wv_ = wv_ + W[i+wid*ldw] * V[i+vid*ldv];
  }
  local_sum_cuda<Tdata,Tcalc>(wv_,&wv[wvid],threadIdx.x);

  return ;
}


// hki(1,2,...,i) = w * v(1:n,1:n*2,...,1:n*i)
template <class Tdata,class Tcalc>
__global__
void CGS_Orth_dot_local_(
		  Tdata *W  ,PARCEL_INT ldw,
		  Tdata *V  ,PARCEL_INT ldv,
		  Tdata *hki,
		  PARCEL_INT n,
		  PARCEL_INT vrow
){
  // gridDim = (1024,512,...) * i 
  // const int block_id  = blockIdx.x;
  // const int block_step  = gridDim.x;
  const int thread_id = threadIdx.x;
  const int thread_step = blockDim.x;

  const int block_step = gridDim.x/vrow;
  const int kid = blockIdx.x / (block_step) ;
  const int block_id = blockIdx.x % block_step;
  const int iid = thread_id + block_id*thread_step;
  Tdata hki_ = 0.0;
  for(PARCEL_INT i = iid;i<n;i=i+thread_step*block_step){
    hki_ = hki_ + V[i+kid*ldv] * W[i];
  }
  local_sum_cuda<Tdata,Tcalc>(hki_,&hki[kid],threadIdx.x);

  return ;
}

// w = sum[k=1:inner] hki*v
template <class Tdata,class Tcalc>
__global__
void CGS_Orth_axpy_(
		  Tdata *W  ,PARCEL_INT ldw,
		  Tdata *V  ,PARCEL_INT ldv,
		  Tdata *hki,
		  PARCEL_INT n,
		  PARCEL_INT vrow
){

  const int thread_id = threadIdx.x;
  const int block_id  = blockIdx.x;
  const int block_step  = gridDim.x;
  const int thread_step = blockDim.x;

  const int nblock = blockDim.x;
  for(PARCEL_INT jj=1+block_id*nblock;jj<=n;jj=jj+block_step*nblock){
    PARCEL_INT jjmin = (jj+nblock-1) < n ? (jj+nblock-1) : n ;
    for(PARCEL_INT k=0;k<vrow;k++){
      Tdata hki_ = hki[k];
      //      for(int j=jj + thread_id ;j<=jjmin;j=j + thread_step){
      for(PARCEL_INT j=jj+thread_id-1;j<=jjmin-1;j=j + thread_step){
	W[j] = W[j] - hki_ * V[j+k*ldv];
      }
    }
  }

  return ;

}

template <class Tdata,class Tcalc>
__global__
void calcQ_VRinv_(
		 PARCEL_INT s,
		 PARCEL_INT n,
		 Tdata *V, PARCEL_INT ldv,
		 Tdata *Q, PARCEL_INT ldq,
		 Tdata *Rinv, PARCEL_INT ldr
		 ){
  // do jj=1,n,nvecBlock 
  //    do i=1,s
  //       nblock=min(n,jj+nvecBlock-1)
  //       do j=jj,nblock
  //          Q(j+(i-1)*n)=0.0d0
  //       enddo
  //       do k=1,i
  //          do j=jj,nblock
  //             Q(j+(i-1)*n)= &
  //                  Q(j+(i-1)*n) + V(j+(k-1)*n) * Rinv(k+(i-1)*s)
  //          enddo
  //       enddo
  //    enddo
  // enddo
  //
  // ------------------------
  //
  int nvecBlock = blockDim.x; // blockDim.x * (gridDim.y*gridDim.z);
  int jj0 = (blockIdx.y + blockIdx.z*gridDim.y)*nvecBlock;

  for(PARCEL_INT jj = 1 + jj0; jj<=n ; jj=jj+nvecBlock * (gridDim.y*gridDim.z) ){ // block : nvecBlock=num_threads*numblocks
    PARCEL_INT nblock = n > (jj+nvecBlock-1)? (jj+nvecBlock-1) : n;
    for(PARCEL_INT i = 1; i<=s ;i++){
      PARCEL_INT j0 = jj + threadIdx.x;
      for(PARCEL_INT j=j0;j<=nblock;j=j+blockDim.x){ // threads
	Q[ j+(i-1)*ldq -1] = 0.0;
      } 

      for(PARCEL_INT k=1;k<=i;k++){
	for(PARCEL_INT j=j0;j<=nblock;j=j+blockDim.x){ 
	  Q[ j+(i-1)*ldq -1] = 
	    Q[j+(i-1)*ldq-1] + V[j+(k-1)*ldv-1] * Rinv[k+(i-1)*ldr-1];
	} 
      }
    }

  }

  return ;
}


extern "C" {
  void cuda_CGS_Orth(
			 MPI_Fint f_comm,
			 devicemem *W_,PARCEL_INT offset_W,PARCEL_INT ldw,
			 devicemem *V_,PARCEL_INT offset_V,PARCEL_INT ldv,
			 devicemem *hki_,PARCEL_INT offset_h,
			 PARCEL_INT n,
			 PARCEL_INT vrow
			 ){

    double* W   = (double*)W_[0].a_d;
    double* V   = (double*)V_[0].a_d;
    double* hki = (double*)hki_[0].a_d;

    cuda_set_vec_val_<double,double> <<<1,32>>>(vrow,0.0,&hki[offset_h]);

    //    int numGPUblock  = 512*vrow;
    int64_t numGPUblock  = 512*vrow;
    while(numGPUblock>16384){
      numGPUblock = numGPUblock /2;
    }
    const int numGPUthread = 512;
    // dot
    CGS_Orth_dot_local_<double,double> 
      <<<numGPUblock,numGPUthread>>>(				  
		    &W[offset_W],ldw,
		    &V[offset_V],ldv,
		    &hki[offset_h],
		    n,
		    vrow
						  );
   double sbuf[vrow];
   double rbuf[vrow];
   for(PARCEL_INT i=0;i<vrow;i++){
     sbuf[i] = 0.0;
     rbuf[i] = 0.0;
   }
   MPI_Comm comm = MPI_Comm_f2c(f_comm);
   const size_t commSize = sizeof(double)*vrow;
   cudaMemcpy(&sbuf[0],&hki[offset_h],commSize,cudaMemcpyDeviceToHost);

   cudaDeviceSynchronize();
   //   MPI_Allreduce(&sbuf,&rbuf,vrow,MPI_DOUBLE,MPI_SUM,comm);
   MPI_Allreduce(sbuf,rbuf,vrow,MPI_DOUBLE,MPI_SUM,comm);
   cudaMemcpy(&hki[offset_h],&rbuf[0],commSize,cudaMemcpyHostToDevice);
   cudaDeviceSynchronize();

   // axpy
   CGS_Orth_axpy_<double,double> 
     <<<numGPUblock,numGPUthread>>>(				  
				    &W[offset_W],ldw,
				    &V[offset_V],ldv,
				    &hki[offset_h],
				    n,
				    vrow
								  );

  }

  void cuda_gemm_dot_local(
			 devicemem *W_,PARCEL_INT offset_W,PARCEL_INT ldw,
			 devicemem *V_,PARCEL_INT offset_V,PARCEL_INT ldv,
			 devicemem *wv_,PARCEL_INT offset_wv,
			 PARCEL_INT n,
			 PARCEL_INT wrow,
			 PARCEL_INT vrow
			 ){
#if 0
    if( (wrow==1) && (vrow==1) ){
      cuda_dot_local(
	       512,512,
	       n,
	       wv_ , offset_wv,
	       W_  , offset_W,
	       V_  , offset_V);
      return ;
    }
#endif
    double* W   = (double*)W_[0].a_d;
    double* V   = (double*)V_[0].a_d;
    double* wv = (double*)wv_[0].a_d;
    cuda_set_vec_val_<double,double> <<<1,32>>>(wrow*vrow,0.0,&wv[offset_wv]);

    //    int numGPUblock  = 512*wrow*vrow;
    int64_t numGPUblock  = 512*wrow*vrow;
    while(numGPUblock>16384){
      numGPUblock = numGPUblock /2;
    }
    const int numGPUthread = 512;
    // dot
    gemm_dot_local_<double,double> 
      <<<numGPUblock,numGPUthread>>>(
				     &W[offset_W],ldw,
				     &V[offset_V],ldv,
				     &wv[offset_wv],
				     n,
				     wrow,
				     vrow
				     );
    return ;
  }

  void cuda_gemm_dot(
			 MPI_Fint f_comm,
			 devicemem *W_,PARCEL_INT offset_W,PARCEL_INT ldw,
			 devicemem *V_,PARCEL_INT offset_V,PARCEL_INT ldv,
			 devicemem *wv_,PARCEL_INT offset_wv,
			 PARCEL_INT n,
			 PARCEL_INT wrow,
			 PARCEL_INT vrow
			 ){
#if 0
    if( (wrow==1) && (vrow==1) ){
      cuda_dot(
	       f_comm,
	       512,512,
	       n,
	       wv_ , offset_wv,
	       W_  , offset_W,
	       V_  , offset_V);
      return ;
    }
#endif

    double* W   = (double*)W_[0].a_d;
    double* V   = (double*)V_[0].a_d;
    double* wv = (double*)wv_[0].a_d;
    cuda_set_vec_val_<double,double> <<<1,32>>>(wrow*vrow,0.0,&wv[offset_wv]);

#if 1
    int64_t numGPUblock  = 512*wrow*vrow;
    while(numGPUblock>16384){
      numGPUblock = numGPUblock /2;
    }
#else
    const int numGPUblock  = 2048;
#endif
    const int numGPUthread = 512;
    // dot
    gemm_dot_local_<double,double> 
      <<<numGPUblock,numGPUthread>>>(
				     &W[offset_W]  ,ldw,
				     &V[offset_V]  ,ldv,
				     &wv[offset_wv],
				     n,
				     wrow,
				     vrow
				     );
    double sbuf[vrow*wrow];
    double rbuf[vrow*wrow];
    //    double *sbuf = (double*)malloc(sizeof(double)*vrow*wrow);
    //    double *rbuf = (double*)malloc(sizeof(double)*vrow*wrow);

   for(PARCEL_INT i=0;i<vrow*wrow;i++){
     sbuf[i] = 0.0;
     rbuf[i] = 0.0;
     //     rbuf[i] = 1.0;
   }
   MPI_Comm comm = MPI_Comm_f2c(f_comm);
   const size_t commSize = sizeof(double)*vrow*wrow;
   cudaMemcpy(sbuf,&wv[offset_wv],commSize,cudaMemcpyDeviceToHost);
   cudaDeviceSynchronize();

   //@@@@@
   MPI_Allreduce(sbuf,rbuf,vrow*wrow,MPI_DOUBLE,MPI_SUM,comm);
   //   MPI_Allreduce(sbuf,rbuf,vrow*wrow,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

   cudaMemcpy(&wv[offset_wv],rbuf,commSize,cudaMemcpyHostToDevice);
   cudaDeviceSynchronize();

   //   free(sbuf);
   //   free(rbuf);
    return ;
  }

}


#if 1
template <class Tdata,class Tcalc>
__global__
void cuda_CGS_Orthonormalization_(
				  Tdata *V,PARCEL_INT ldv,
				  Tdata *Q,PARCEL_INT ldq,
				  Tdata *R,PARCEL_INT ldr,
				  PARCEL_INT n,
				  PARCEL_INT s,
				  PARCEL_INT i
				  ){
  //  const int block_id  = blockIdx.x;
  //  const int block_step  = gridDim.x;
  //  const int thread_id = threadIdx.x;
  const int thread_step = blockDim.x;
  //  int nvecBlock = thread_step;
  //  for(int jj=1;jj<=n;jj=jj+nvecBlock){
#if 1
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  for(PARCEL_INT j=tid;j<n;j=j+step){
    for(PARCEL_INT k=1;k<=i-1;k++){
      V[j]=V[j]-R[(k-1)+ldr*(i-1)]*Q[j+ldq*(k-1)];
    }
  }
#else
  for(PARCEL_INT jj=1+block_id;jj<=n;jj=jj+nvecBlock*block_step){
    PARCEL_INT nblock = n > (jj+nvecBlock-1) ? (jj+nvecBlock-1) :n;
    for(PARCEL_INT k=1;k<=i-1;k++){
      for(PARCEL_INT j=jj+thread_id;j<=jj+nblock-1;j=j+thread_step){
	V[j-1]=V[j-1]-R[k+ldr*(i-1)-1]*Q[j+ldq*(k-1)-1];
      }
    }
  }
#endif
  return ;
}

extern "C" {
  void cuda_CGS_Orthonormalization(
				   devicemem *V_,PARCEL_INT offset_V,PARCEL_INT ldv,
				   devicemem *Q_,PARCEL_INT offset_Q,PARCEL_INT ldq,
				   devicemem *R_,PARCEL_INT offset_R,PARCEL_INT ldr,
				   PARCEL_INT n,
				   PARCEL_INT s,
				   PARCEL_INT i
				   ){

    double* V = (double*)V_[0].a_d;
    double* Q = (double*)Q_[0].a_d;
    double* R = (double*)R_[0].a_d;

    cuda_CGS_Orthonormalization_<double,double> 
      <<<512,512>>>(				  
      //      <<<1,1>>>(				  
      //<<<numGPUblock,numGPUthread>>>(				  
		    &V[offset_V],ldv,
		    &Q[offset_Q],ldq,
		    &R[offset_R],ldr,
		    n,
		    s,
		    i
						  );

    return ;
  }

  void cuda_calcQ_VRinv(
		 PARCEL_INT s,
		 PARCEL_INT n,
		 devicemem *V_   , PARCEL_INT offset_V   , PARCEL_INT ldv,
		 devicemem *Q_   , PARCEL_INT offset_Q   , PARCEL_INT ldq,
		 devicemem *Rinv_, PARCEL_INT offset_Rinv, PARCEL_INT ldrinv
			){

    double* V    = (double*)V_[0].a_d;
    double* Q    = (double*)Q_[0].a_d;
    double* Rinv = (double*)Rinv_[0].a_d;
    
    // dim3 block(512,1,1); 
    // dim3 grid(1,32,32); 

    dim3 block(512,1,1); 
    dim3 grid(1,32,32); 

    calcQ_VRinv_<double,double>
      <<<grid,block>>>(
				  s,
				  n,
				  &V[offset_V],ldv,
				  &Q[offset_Q],ldq,
				  &Rinv[offset_Rinv],ldrinv
				  );
    return ;
  }

}    

#endif
