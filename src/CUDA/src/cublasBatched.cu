#include <cuda.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>

#include <mpi.h>
#include "cu_common.h"
#include "cu_common.hpp"

#include "reduciton_cuda.hpp"

template <class Tdata,class Tcalc>
__global__
void reduceBatched_matrices_(
		    const int batchCount, 
		    const int local_nrow, const int local_ncol,
		    double *batchedMat, int ldbatbatched,
		    double *reduceMat , int ldreduce){

  // gridDim.xはlocal_nrowの倍数でなければならない
  const int nrow = batchCount * local_nrow;

  const int iid    = threadIdx.x + blockIdx.y*blockDim.x;
  const int iblock = blockDim.x * gridDim.y;

  const int reduceMat_row_id = iid%local_nrow;
  const int reduceMat_col_id = blockIdx.z;
  const int reduceMat_id     = reduceMat_row_id + reduceMat_col_id*ldreduce;
  // gridDim.x == 1でなければならない

  Tdata local_reduce = 0.0;
  for(PARCEL_INT i = iid ; i<nrow ; i = i+iblock){
    local_reduce = local_reduce + batchedMat[i + reduceMat_col_id*ldbatbatched];
  }
  // printf("iid = %06d, reduce_id(%4d,%4d) = %04d : iblock = %d, batchCount = %d  : local_reduce = %f \n",
  // 	 iid,
  // 	 reduceMat_row_id,reduceMat_col_id,
  // 	 reduceMat_id,
  // 	 iblock,batchCount,local_reduce);

  //  for(int j = 0; j<local_nrow ; j++){
  const PARCEL_INT j = reduceMat_col_id;
  for(PARCEL_INT i = 0; i<local_nrow ; i++){    
    const PARCEL_INT ij = i + j*ldreduce;
    Tdata local_reduce_ = reduceMat_id==ij ? local_reduce : 0.0;
    local_sum_cuda<Tdata,Tcalc>(local_reduce_,&reduceMat[ij],threadIdx.x);
  }
}


void reduceBatched_matrices(
		     const int batchCount, 
		     const int local_nrow, const int local_ncol,
		     double *batchedMat, int ldbatbatched,
		     double *reduceMat , int ldreduce){
  
  cuda_set_vec_val_<double,double> 
    <<<1,local_ncol*ldreduce>>>(local_ncol*ldreduce,0.0,reduceMat);

  int local_nrow_para = local_nrow < 1024 ? local_nrow : 512; 
  if(local_nrow<1024){  
    dim3 block(local_nrow_para,1,1); 
    dim3 grid(1,32,local_ncol);
    reduceBatched_matrices_<double,double>
      <<<grid,block>>>(
		       batchCount, 
		       local_nrow,   local_ncol,
		       batchedMat,  ldbatbatched,
		       reduceMat ,  ldreduce
		       );
  }else{
    dim3 block(512,1,1); 
    dim3 grid(1,32,local_nrow_para);
    reduceBatched_matrices_<double,double>
      <<<grid,block>>>(
		       batchCount, 
		       local_nrow,   local_ncol,
		       batchedMat,  ldbatbatched,
		       reduceMat ,  ldreduce
		       );
  }

  return ;
}

extern "C" {

#if 1
  void swap_ptr_ptr_(void ***x, void ***y) {
    void **tmp;
    tmp = *y;
    *y = *x;
    *x = tmp;
  }
#else
  void swap_ptr_ptr_(void **x,void **y){
    void *tmp;
    tmp=*y;
    *y=*x;
    *x=tmp;
  }
#endif
  void swap_ptr_ptr(devicemem *d_Aptr_,devicemem *d_Bptr_) {
    swap_ptr_ptr_(
		  (void ***)(d_Aptr_[0].a_d) ,
		  (void ***)(d_Bptr_[0].a_d)
		  );
  }

  void cublas_init_DgemmBatched_gemm_dot_local(
					      PARCEL_INT batchCount,
					      devicemem *d_A_, PARCEL_INT offset_A, PARCEL_INT lda, 
					      devicemem *d_Aptr_
					      ){
#if 0
    // int64_t が正しい?
    int64_t **h_Aptr = (int64_t**)malloc( batchCount*sizeof(int64_t*) );
    double *d_A = (double *)d_A_[0].a_d;
  
    for(PARCEL_INT i = 0;i < batchCount;i++){
      h_Aptr[i] = &d_A[i*lda+offset_A]; // lda = localnRow
    } 

    cudaMalloc((void**)&d_Aptr_[0].a_d, batchCount*sizeof(int64_t*));
    int64_t **d_Aptr = (int64_t **)d_Aptr_[0].a_d;
    cudaMemcpy(d_Aptr, h_Aptr, batchCount*sizeof(int64_t*), cudaMemcpyHostToDevice);
    
#else
    double **h_Aptr = (double**)malloc( batchCount*sizeof(double*) );
    double *d_A = (double *)d_A_[0].a_d;
  
    for(PARCEL_INT i = 0;i < batchCount;i++){
      h_Aptr[i] = &d_A[i*lda+offset_A]; // lda = localnRow
    } 

    cudaMalloc((void**)&d_Aptr_[0].a_d, batchCount*sizeof(double*));
    double **d_Aptr = (double **)d_Aptr_[0].a_d;
    cudaMemcpy(d_Aptr, h_Aptr, batchCount*sizeof(double*), cudaMemcpyHostToDevice);
    
#endif
    free(h_Aptr);
    return ;

  }
  
  void cublas_free_DgemmBatched_gemm_dot_local(devicemem *d_Aptr_){
    cudaFree(d_Aptr_[0].a_d);
  }

  void cublas_DgemmBatched_gemm_dot_local(
					  // devicemem *handle_,
					  f_cublas *f_cublas_,
					  PARCEL_INT batchCount,
					  PARCEL_INT nCol, PARCEL_INT localnRow,
					  devicemem *Aptr_, PARCEL_INT ldAptr,
					  devicemem *Bptr_, PARCEL_INT ldBptr,
					  devicemem *Cptr_, PARCEL_INT ldCptr,
					  devicemem *Cmat_,  
					  devicemem *reduce_mat_ ){
    double alpha = 1.0; 
    double beta  = 0.0;
    cublasHandle_t *handle = (cublasHandle_t *) f_cublas_[0].cublasHandle;

    double *Cmat  = (double *)Cmat_[0].a_d;
    double *reduce_mat  = (double *)reduce_mat_[0].a_d;

    double **Aptr = (double **)Aptr_[0].a_d;
    double **Bptr = (double **)Bptr_[0].a_d;
    double **Cptr = (double **)Cptr_[0].a_d; // Cmatのメモリポインタ 

    int batchCount_ = batchCount;
    int ldAptr_ = ldAptr;
    int ldBptr_ = ldBptr;
    int ldCptr_ = ldCptr;

    int nCol_ = nCol;
    int localnRow_ = localnRow;

    cublasDgemmBatched(
		       // handle[0],
		       PARCEL_CUBLAS_HANDLE,
		       CUBLAS_OP_T,
		       CUBLAS_OP_N,
		       nCol_, nCol_, localnRow_,
		       &alpha,
		       Aptr, ldAptr_, // nRow,
		       Bptr, ldAptr_, // nRow,
		       &beta,
		       Cptr, nCol_*batchCount_, // 
		       batchCount_);  

    reduceBatched_matrices(
			   batchCount_,
			   nCol_, nCol_,
			   Cmat, nCol_*batchCount_,
			   reduce_mat, nCol_);

    // reduce_mat --> Cmat
    //    const size_t mem_size = nCol*nCol*sizeof(double);
    //    cudaMemcpy(Cmat,reduce_mat,mem_size,cudaMemcpyDeviceToDevice);

    return ;
  }

  void cublas_DgemmBatched_gemm_dot(
				    MPI_Fint f_comm,
				    f_cublas *f_cublas_,
				    PARCEL_INT batchCount,
				    PARCEL_INT nCol, PARCEL_INT localnRow,
				    devicemem *Aptr_, PARCEL_INT ldAptr,
				    devicemem *Bptr_, PARCEL_INT ldBptr,
				    devicemem *Cptr_, PARCEL_INT ldCptr,
				    devicemem *Cmat_,  
				    devicemem *reduce_mat_ ){

    cublas_DgemmBatched_gemm_dot_local(
				       f_cublas_,
				       batchCount,
				       nCol, localnRow,
				       Aptr_, ldAptr,
				       Bptr_, ldBptr,
				       Cptr_, ldCptr,
				       Cmat_, 
				       reduce_mat_ );

   double sbuf[(nCol*nCol)];
   double rbuf[(nCol*nCol)];
   for(PARCEL_INT i=0;i<(nCol*nCol);i++){
     sbuf[i] = 0.0;
     rbuf[i] = 0.0;
   }

   MPI_Comm comm = MPI_Comm_f2c(f_comm);
   const size_t commSize = sizeof(double)*(nCol*nCol);

    double *reduce_mat  = (double *)reduce_mat_[0].a_d;
   cudaMemcpy(&sbuf[0],reduce_mat,commSize,cudaMemcpyDeviceToHost);
   cudaDeviceSynchronize();
   MPI_Allreduce(sbuf,rbuf,(nCol*nCol),MPI_DOUBLE,MPI_SUM,comm);

   cudaMemcpy(reduce_mat,&rbuf[0],commSize,cudaMemcpyHostToDevice);
   cudaDeviceSynchronize();

    return ;
  }

}
