#include <stdio.h>
//#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "cu_common.h"
#include "cu_common.hpp"
#include <cusparse.h>
#include <cublas_v2.h>
#include <mpi.h>

cublasHandle_t PARCEL_CUBLAS_HANDLE;
// cusparseHandle_t PARCEL_CUSPARSE_SPMV_HANDLE;

__global__
void cuda_FP64toBF16_(
	       PARCEL_INT n,
	       double    *x,
	       __nv_bfloat16 *y,
	       int mode
	       ){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  switch(mode){
  case FP64TOBF16_ROUND_NEAR : 
    for(PARCEL_INT i=tid;i<n;i=i+step){
      float x_fp32 = x[i];
      y[i] = __float2bfloat16_rn(x_fp32);
    }
    break;
  case FP64TOBF16_ROUND_ZERO : 
    for(PARCEL_INT i=tid;i<n;i=i+step){
      float x_fp32 = x[i];
      y[i] = __float2bfloat16_rz(x_fp32);
    }
    break;
  case FP64TOBF16_ROUND_DOWN : 
    for(PARCEL_INT i=tid;i<n;i=i+step){
      float x_fp32 = x[i];
      y[i] = __float2bfloat16_rd(x_fp32);
    }
    break;
  case FP64TOBF16_ROUND_UP : 
    for(PARCEL_INT i=tid;i<n;i=i+step){
      float x_fp32 = x[i];
      y[i] = __float2bfloat16_ru(x_fp32);
    }
    break;
  default : 
    ;
    break;
  }

  return ;

}

template <class Tdata>
__global__
void cuda_ptr_type_view_(
	       PARCEL_INT n,
	       Tdata    *x
	       ){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  for(PARCEL_INT i=tid;i<n;i=i+step){
    //    printf(" x[%d] %f \n",i,x[i]);
    printf(" x[%d] %f \n",i,x[i]);
  }
  return ;
}

__global__
void cuda_ptr_parcel_int_view_(
	       PARCEL_INT n,
	       int    *x
	       ){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  for(PARCEL_INT i=tid;i<n;i=i+step){
    printf(" x[%d] %d \n",i,x[i]);
    //    printf(" x[%d] \n",i);
  }
  return ;
}

extern "C"{
  void cudaMalloc_FP64(devicemem *a_d,PARCEL_INT n);
  void cudaFree_(devicemem *a_d);

#if 0
  void cuda_ptr_parcel_int_view_(
		     PARCEL_INT n,
		     void *x
		     ){
    cudaDeviceSynchronize();
    cuda_ptr_parcel_int_view_ <<<1,1>>>(n,(PARCEL_INT *)x);
    cudaDeviceSynchronize();
    return ;
  }
  void cuda_view_PARCEL_INT(
		 PARCEL_INT n,devicemem *x_
		 ){
    PARCEL_INT *x = (PARCEL_INT *)x_[0].a_d;
    //    cuda_ptr_parcel_int_view_<PARCEL_INT> <<<1,1>>>(n,x);
    cuda_ptr_parcel_int_view_ <<<1,1>>>(n,x);
    return ;
  }

  void cuda_ptr_double_view_(
		     PARCEL_INT n,
		     void *x
		     ){
    //    printf(" cuda_ptr_double_view : \n");
    cudaDeviceSynchronize();
    cuda_ptr_type_view_<double> <<<1,1>>>(n,(double *)x);
    cudaDeviceSynchronize();
    return ;
  }

  void cuda_view_FP64(
		 PARCEL_INT n,devicemem *x_
		 ){
    double *x = (double *)x_[0].a_d;
    cuda_ptr_type_view_<double> <<<1,1>>>(n,x);
    return ;
  }
#endif

  void cuda_sqrt_val(
		     int numGPUblock,
		     int numGPUthread,
		     devicemem *x_,PARCEL_INT offset_x){

    double *x = (double *)x_[0].a_d;
    cuda_sqrt_val_<double,double> <<<numGPUblock,numGPUthread>>>(&x[offset_x]);
    return ;
  }
  void cuda_add_val(
		     int numGPUblock,
		     int numGPUthread,
		     double val,
                     devicemem *x_,PARCEL_INT offset_x){ 
    double *x = (double *)x_[0].a_d;
    cuda_add_val_<double,double> 
      <<<numGPUblock,numGPUthread>>>(val,&x[offset_x]);
    return ;
  }
  void cuda_sub_val(
		     int numGPUblock,
		     int numGPUthread,
		     double val,
                     devicemem *x_,PARCEL_INT offset_x){ 
    double *x = (double *)x_[0].a_d;
    cuda_sub_val_<double,double> <<<numGPUblock,numGPUthread>>>(val,&x[offset_x]);
    return ;
  }
  void cuda_mul_val(
		     int numGPUblock,
		     int numGPUthread,
		     double val,
                     devicemem *x_,PARCEL_INT offset_x){ 
    double *x = (double *)x_[0].a_d;
    cuda_mul_val_<double,double> <<<numGPUblock,numGPUthread>>>(val,&x[offset_x]);
    return ;
  }
  void cuda_div_val(
		     int numGPUblock,
		     int numGPUthread,
		     double val,
                     devicemem *x_,PARCEL_INT offset_x){ 
    double *x = (double *)x_[0].a_d;
    cuda_div_val_<double,double> <<<numGPUblock,numGPUthread>>>(val,&x[offset_x]);
    return ;
  }

  void cuda_calc_dgemmNN_dot(
			     MPI_Fint f_comm,
				f_cublas *f_cublas_,
				PARCEL_INT m,
				PARCEL_INT n,
				PARCEL_INT k,
				devicemem *alpha_,PARCEL_INT offset_alpha,
				devicemem *beta_,PARCEL_INT offset_beta,
				devicemem *A_,PARCEL_INT offset_A,PARCEL_INT lda,
				devicemem *B_,PARCEL_INT offset_B,PARCEL_INT ldb,
				devicemem *C_,PARCEL_INT offset_C,PARCEL_INT ldc
				){

    cublasHandle_t *cublasHandle = 
      (cublasHandle_t *) f_cublas_[0].cublasHandle;
    double *alpha = (double *)alpha_[0].a_d;
    double *beta = (double *)beta_[0].a_d; 
    double *A = (double *)A_[0].a_d; 
    double *B = (double *)B_[0].a_d; 
    double *C = (double *)C_[0].a_d;  

    // Column Major
    // transa = CUBLAS_OP_N , transb = CUBLAS_OP_N
    // C = alpha*A[m*k]*B[k*n] + beta*C[m*n]
    cublasDgemm(
		//		cublasHandle[0],
		PARCEL_CUBLAS_HANDLE,
		//		PARCEL_CUBLAS_HANDLE,
		CUBLAS_OP_N,CUBLAS_OP_N,
		m,  n,  k,
		&alpha[offset_alpha],
		&A[offset_A],lda,
		&B[offset_B],ldb,
		&beta[offset_beta],
		&C[offset_C],ldc);
   MPI_Comm comm = MPI_Comm_f2c(f_comm);

   const size_t commSize = sizeof(double)*m*n;
   double sbuf[m*n];
   double rbuf[m*n];
   cudaMemcpy(&sbuf,&C[offset_C],commSize,cudaMemcpyDeviceToHost);
   cudaDeviceSynchronize();
   MPI_Allreduce(
		 sbuf,rbuf,
		 m*n,MPI_DOUBLE,MPI_SUM,comm);
   cudaMemcpy(&C[offset_C],&rbuf,commSize,cudaMemcpyHostToDevice);

  }

  void cuda_calc_dgemmTN_dot(
			     MPI_Fint f_comm,
			     f_cublas *f_cublas_,
			     PARCEL_INT m,
			     PARCEL_INT n,
			     PARCEL_INT k,
			     devicemem *alpha_,PARCEL_INT offset_alpha,
			     devicemem *beta_,PARCEL_INT offset_beta,
			     devicemem *A_,PARCEL_INT offset_A,PARCEL_INT lda,
			     devicemem *B_,PARCEL_INT offset_B,PARCEL_INT ldb,
			     devicemem *C_,PARCEL_INT offset_C,PARCEL_INT ldc
			     ){
    cublasHandle_t *cublasHandle = 
      (cublasHandle_t *) f_cublas_[0].cublasHandle;
    double *alpha = (double *)alpha_[0].a_d;
    double *beta = (double *)beta_[0].a_d; 
    double *A = (double *)A_[0].a_d; 
    double *B = (double *)B_[0].a_d; 
    double *C = (double *)C_[0].a_d;  

    // Column Major
    // transa = CUBLAS_OP_N , transb = CUBLAS_OP_N
    // C = alpha*A[m*k]*B[k*n] + beta*C[m*n]
    // transa = CUBLAS_OP_N , transb = CUBLAS_OP_N
    // C = alpha*A[k*m]*B[k*n] + beta*C[m*n]
    cublasDgemm(
		//		cublasHandle[0],
		       PARCEL_CUBLAS_HANDLE,
		CUBLAS_OP_T,CUBLAS_OP_N,
		m,  n,  k,
		&alpha[offset_alpha],
		&A[offset_A],lda,
		&B[offset_B],ldb,
		&beta[offset_beta],
		&C[offset_C],ldc);

   MPI_Comm comm = MPI_Comm_f2c(f_comm);
   const size_t commSize = sizeof(double)*m*n;
   double sbuf[m*n];
   double rbuf[m*n];
   cudaMemcpy(&sbuf,&C[offset_C],commSize,cudaMemcpyDeviceToHost);
   cudaDeviceSynchronize();
   MPI_Allreduce(
		 sbuf,rbuf,
		 m*n,MPI_DOUBLE,MPI_SUM,comm);
   cudaMemcpy(&C[offset_C],&rbuf,commSize,cudaMemcpyHostToDevice);
   
  }


  void cuda_axpy(
		 int numGPUblock,
		 int numGPUthread,
		 PARCEL_INT n,
		 double alpha,
		 double beta,
		 devicemem *x_,PARCEL_INT offset_x,
		 devicemem *y_,PARCEL_INT offset_y
		 ){
    double *x = (double *)x_[0].a_d;
    double *y = (double *)y_[0].a_d;
    cuda_axpy_<double,double> <<<numGPUblock,numGPUthread>>>(
							     n,
							     alpha,
							     beta,
							     &x[offset_x],
							     &y[offset_y]
							     );

    return ;
  }

  void cuda_zero_init(
		 PARCEL_INT n,
		 devicemem *x_,PARCEL_INT offset_x
		 ){
    double *x = (double *)x_[0].a_d;
    cuda_zero_init_<double,double> <<<512,512>>>(n,&x[offset_x]);
    //    cuda_zero_init_<double,double> <<<512,512>>>(1,&x[offset_x]);
    return ;
  }

  void cuda_axpy_vec3(
		 int numGPUblock,
		 int numGPUthread,
		 PARCEL_INT n,
		 double alpha,
		 double beta,
		 double gamma,
		 devicemem *x_,PARCEL_INT offset_x,
		 devicemem *y_,PARCEL_INT offset_y,
		 devicemem *z_,PARCEL_INT offset_z
		 ){
    double *x = (double *)x_[0].a_d;
    double *y = (double *)y_[0].a_d;
    double *z = (double *)z_[0].a_d;
    cuda_axpy_vec3_<double,double> <<<numGPUblock,numGPUthread>>>(
							     n,
							     alpha,
							     beta,
							     gamma,
							     &x[offset_x],
							     &y[offset_y],
							     &z[offset_z]
							     );

    return ;
  }

  void cuda_axpy_vec4(
		 int numGPUblock,
		 int numGPUthread,
		 PARCEL_INT n,
		 double alpha,
		 double beta,
		 double gamma,
		 double delta,
		 devicemem *x_,PARCEL_INT offset_x,
		 devicemem *y_,PARCEL_INT offset_y,
		 devicemem *z_,PARCEL_INT offset_z,
		 devicemem *w_,PARCEL_INT offset_w
		 ){
    double *x = (double *)x_[0].a_d;
    double *y = (double *)y_[0].a_d;
    double *z = (double *)z_[0].a_d;
    double *w = (double *)w_[0].a_d;
    cuda_axpy_vec4_<double,double> <<<numGPUblock,numGPUthread>>>(
							     n,
							     alpha,
							     beta,
							     gamma,
							     delta,
							     &x[offset_x],
							     &y[offset_y],
							     &z[offset_z],
							     &w[offset_w]
							     );

    return ;
  }

  void cuda_gemv(
		 int numGPUblock,
		 int numGPUthread,
		 devicemem *x_,PARCEL_INT offset_x, PARCEL_INT ldx,
		 devicemem *y_,PARCEL_INT offset_y, PARCEL_INT ldy,
		 devicemem *alpha_,PARCEL_INT offset_alpha, 
		 PARCEL_INT n,
		 PARCEL_INT kmax
		 ){

    double *x     = (double *)x_[0].a_d;
    double *y     = (double *)y_[0].a_d;
    double *alpha = (double *)alpha_[0].a_d;

    cuda_gemv_<double,double> 
      <<<numGPUblock,numGPUthread>>>(
				     &x[offset_x], ldx,
				     &y[offset_y], ldy,
				     &alpha[offset_alpha],
				     n,
				     kmax
				     );    
    return ;
  }

  void cuda_vec_scaling(
		 int numGPUblock,
		 int numGPUthread,
		 PARCEL_INT n,
		 devicemem *x_, PARCEL_INT offset_x,
		 devicemem *y_, PARCEL_INT offset_y,
		 devicemem *z_, PARCEL_INT offset_z
		 ){

    double* x = (double*)x_[0].a_d;
    double* y = (double*)y_[0].a_d;
    double* z = (double*)z_[0].a_d;
    cuda_vec_scaling_<double,double> <<<numGPUblock,numGPUthread>>>(
								    n,
								  &x[offset_x],
								  &y[offset_y],
							          &z[offset_z]
							     );
    return ;
  }

  void cuda_dot(
		 MPI_Fint f_comm,
		 int numGPUblock,
		 int numGPUthread,
		 PARCEL_INT n,
		 devicemem *val_,PARCEL_INT offset_val,
		 devicemem *x_,PARCEL_INT offset_x,
		 devicemem *y_,PARCEL_INT offset_y
		){
    //    printf(" cuda_dot : n = %d \n",n);

    double *val = (double *)val_[0].a_d;
    double *x = (double *)x_[0].a_d;
    double *y = (double *)y_[0].a_d;
    cuda_set_vec_val_<double,double> <<<1,32>>>(1,0.0,&val[offset_val]);
    cuda_dot_local_<double,double> <<<numGPUblock,numGPUthread>>>(
			       			        n,
							&val[offset_val],
							&x[offset_x],
							&y[offset_y]
			       				);

   double sbuf = 1.0;
   double rbuf = 1.0;
   MPI_Comm comm = MPI_Comm_f2c(f_comm);
   const size_t commSize = sizeof(double);
   cudaMemcpy(&sbuf,&val[offset_val],commSize,cudaMemcpyDeviceToHost);
   cudaDeviceSynchronize();
   //   printf("sbuf = %f \n",sbuf);
   MPI_Allreduce(&sbuf,&rbuf,1,MPI_DOUBLE,MPI_SUM,comm);
   //   MPI_Allreduce(&sbuf,&rbuf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   cudaMemcpy(&val[offset_val],&rbuf,commSize,cudaMemcpyHostToDevice);

   return ;
  }

  void cuda_dot_local(
		 int numGPUblock,
		 int numGPUthread,
		 PARCEL_INT n,
		 devicemem *val_, PARCEL_INT offset_val,
		 devicemem *x_  , PARCEL_INT offset_x,
		 devicemem *y_  , PARCEL_INT offset_y
		){
    //    printf(" cuda_dot : n = %d \n",n);

    double *val = (double *)val_[0].a_d;
    double *x   = (double *)x_[0].a_d;
    double *y   = (double *)y_[0].a_d;
    cuda_set_vec_val_<double,double> <<<1,32>>>(1,0.0,&val[offset_val]);
    cuda_dot_local_<double,double> <<<numGPUblock,numGPUthread>>>(
			       			        n,
							&val[offset_val],
							&x[offset_x],
							&y[offset_y]
			       				);

   return ;
  }

  void cuda_set_vec_val(
		 int numGPUblock,
		 int numGPUthread,
		 PARCEL_INT n,
		 double val,
		 devicemem *x_,PARCEL_INT offset_x
		 ){
    double *x = (double *)x_[0].a_d;
    cuda_set_vec_val_<double,double> <<<numGPUblock,numGPUthread>>>(
							     n,
							     val,
							     &x[offset_x]
							     );

    return ;
  }

  void init_cudahandle(){
    //  std::cout<<"init_cudahandle "<<std::endl;
    cublasCreate(&PARCEL_CUBLAS_HANDLE);
    cusparseCreate(&PARCEL_CUSPARSE_SPMV_HANDLE);
    cusparseCreate(&PARCEL_CUSPARSE_PRECON_HANDLE);
  }
  void free_cudahandle(){
    cublasDestroy(PARCEL_CUBLAS_HANDLE);
    cusparseDestroy(PARCEL_CUSPARSE_SPMV_HANDLE);
    cusparseDestroy(PARCEL_CUSPARSE_PRECON_HANDLE);
  }
 
  void init_PARCEL_Matrix_GPU( 
			  
			   ){

  }

  void set_GPU(int rank){
    int num_gpu;
    cudaGetDeviceCount(&num_gpu);
    const int selected_gpu = rank % num_gpu;
    cudaError_t cu_err = cudaSetDevice(selected_gpu);
    if (cu_err != cudaSuccess) {
      printf(" cudaSetDevice error \n");
      fprintf(stderr, "[Error] %s (error code: %d) at %s line %d\n", cudaGetErrorString(cu_err), cu_err, __FILE__, __LINE__); \
      //    return -1;
    } 
    return ;
  }
  
  //}
  //extern "C" {

  void cudaDeviceSynchronize_(){
    cudaDeviceSynchronize();
    return ;
  }

  void set_cublas(f_cublas *f_cublas_){
    f_cublas_[0].cublasHandle = (cublasHandle_t *)malloc(sizeof(cublasHandle_t));
    //cublasHandle_t * cublasHandle = (cublasHandle_t *) f_cublas_[0].cublasHandle;
    //cublasCreate(cublasHandle);

    f_cublas_[0].alpha = (devicemem *)malloc(sizeof(devicemem));
    f_cublas_[0].beta  = (devicemem *)malloc(sizeof(devicemem));

    cudaMalloc_FP64(f_cublas_[0].alpha,1);
    cudaMalloc_FP64(f_cublas_[0].beta ,1);

#if 0
    // cublasCreate(&PARCEL_CUBLAS_HANDLE);
    // cusparseCreate(&PARCEL_CUSPARSE_SPMV_HANDLE);
    init_cudahandle();
#endif

    return ;
  }

  void free_cublas(f_cublas *f_cublas_){
    cublasHandle_t * cublasHandle = 
      (cublasHandle_t *) f_cublas_[0].cublasHandle;
    //    cublasDestroy(cublasHandle[0]);

    free(f_cublas_[0].cublasHandle);
    cudaFree_(f_cublas_[0].alpha);
    cudaFree_(f_cublas_[0].beta);
    free(f_cublas_[0].alpha);
    free(f_cublas_[0].beta);

#if 0
    // cublasDestroy(PARCEL_CUBLAS_HANDLE);
    // cusparseDestroy(PARCEL_CUSPARSE_SPMV_HANDLE);
    free_cudahandle();
#endif

    return ;
  }

  void cudaFree_(devicemem *a_d){
    cudaFree(a_d[0].a_d);
    return ;
  }

  void cudaMalloc_FP64(devicemem *a_d,PARCEL_INT n){
    const size_t len = sizeof(double)*n;
    cudaError_t err_code;
    err_code = cudaMalloc(&(a_d[0].a_d), len);

#if 1
    if (err_code != cudaSuccess) {
      printf("error : get cuda Memory size = %ld : %d , %d \n",len,err_code,cudaSuccess);
      exit(-1);
    }
#else
    cuda_malloc_check_err(err_code);
#endif

    return ;
  }

  void cudaMalloc_bf16(devicemem *a_d,PARCEL_INT n){
    const size_t len= sizeof(__nv_bfloat16)*n;
    cudaError_t err_code;
    err_code = cudaMalloc(&(a_d[0].a_d), len);
#if 0
    if (err_code != cudaSuccess) {
      printf("error : get cuda Memory size = %ld : %d , %d \n",len,err_code,cudaSuccess);
      exit(-1);
    }
#else
    cuda_malloc_check_err(err_code);
#endif
    return ;
  }

  void cusparseCreateDnVec_FP64(devicemem *a_d,PARCEL_INT n,PARCEL_INT nAdd,PARCEL_INT offset){
    const cudaDataType_t cudaDataType = CUDA_R_64F;
    double *p = (double *)a_d[0].a_d;

    a_d[0].descr = 
      (cusparseDnVecDescr_t*)malloc(sizeof(cusparseDnVecDescr_t));
    cusparseDnVecDescr_t *vecp = (cusparseDnVecDescr_t *)a_d[0].descr;
    cusparseCreateDnVec(vecp, n, &p[offset], cudaDataType);

    a_d[0].descr_add = 
      (cusparseDnVecDescr_t*)malloc(sizeof(cusparseDnVecDescr_t));
    cusparseDnVecDescr_t *vecp_add = (cusparseDnVecDescr_t *)a_d[0].descr_add;
    cusparseCreateDnVec(vecp_add, nAdd, &p[offset], cudaDataType);
    return ;
  }

  void cusparseDestroyDnVec_(devicemem *a_d){

    cusparseDnVecDescr_t *vecp = (cusparseDnVecDescr_t *)a_d[0].descr;
    cusparseDestroyDnVec(vecp[0]);
    cudaDeviceSynchronize();
    free(a_d[0].descr);

    cusparseDnVecDescr_t *vecp_add = (cusparseDnVecDescr_t *)a_d[0].descr_add;
    cusparseDestroyDnVec(vecp_add[0]);
    cudaDeviceSynchronize();
    free(a_d[0].descr_add);

    return ;
  }

#if 0
  void cusparseDnVecDescr_reset_(devicemem *a_d,PARCEL_INT n,PARCEL_INT nAdd){
    cusparseDestroyDnVec_(a_d);
    cusparseCreateDnVec_FP64(a_d, n, nAdd);
  }
#endif

  void cudaMemcpy_d2h_FP64(
			   devicemem *a_d_, 
			   double *a_h, 
			   PARCEL_INT n, 
			   PARCEL_INT offset_d, 
			   PARCEL_INT offset_h ){
    double *a_d = (double *) a_d_[0].a_d;
    const size_t len = sizeof(double) * n;
    cudaMemcpy( &a_h[offset_h], &a_d[offset_d], len, cudaMemcpyDeviceToHost);
    return ;
  }

  void cudaMemcpy_d2d_FP64(
			   devicemem *a_d_,
			   devicemem *b_d_,
			   PARCEL_INT n, 
			   PARCEL_INT offset_a, 
			   PARCEL_INT offset_b ){

    double *a_d = (double *) a_d_[0].a_d;
    double *b_d = (double *) b_d_[0].a_d;

    const size_t len = sizeof(double) * n;
    cudaMemcpy(
	       &b_d[offset_b],
	       &a_d[offset_a],
	       len,
	       cudaMemcpyDeviceToDevice);
    return ;
  }

  void cudaMemcpy_d2d_int(
			   devicemem *a_d_,
			   devicemem *b_d_,
			   PARCEL_INT n, 
			   PARCEL_INT offset_a, 
			   PARCEL_INT offset_b ){

    int *a_d = (int *) a_d_[0].a_d;
    int *b_d = (int *) b_d_[0].a_d;

    const size_t len = sizeof(int) * n;
    cudaMemcpy(
	       &b_d[offset_b],
	       &a_d[offset_a],
	       len,
	       cudaMemcpyDeviceToDevice);
    return ;
  }

  void cudaMemcpy_d2d_PARCEL_INT(
			   devicemem *a_d_,
			   devicemem *b_d_,
			   PARCEL_INT n, 
			   PARCEL_INT offset_a, 
			   PARCEL_INT offset_b ){

    double *a_d = (double *) a_d_[0].a_d;
    double *b_d = (double *) b_d_[0].a_d;

    const size_t len = sizeof(PARCEL_INT) * n;
    cudaMemcpy(
	       &b_d[offset_b],
	       &a_d[offset_a],
	       len,
	       cudaMemcpyDeviceToDevice);
    return ;
  }


  void cudamemcpy_dptr2d_parcel_int_(
			  PARCEL_INT *aptr,
			  devicemem  *a_d_,
			  PARCEL_INT *n_
			 ){
    PARCEL_INT n = n_[0];
    PARCEL_INT *a_d =(PARCEL_INT *)a_d_[0].a_d;

    cudaDeviceSynchronize();
    cuda_axpy_<PARCEL_INT,PARCEL_INT> <<<512,512>>>(
	       n,
	       1,
	       0,
	       aptr,
	       a_d);
    cudaDeviceSynchronize();
    //    cuda_ptr_type_view_<PARCEL_INT> <<<1,1>>>(10,a_d);
    //    cudaDeviceSynchronize();

    return ;
  }

  void cudamemcpy_dptr2d_int_(
			  int *aptr,
			  devicemem  *a_d_,
			  PARCEL_INT *n_
			 ){
    PARCEL_INT n = n_[0];
    int *a_d =(int *)a_d_[0].a_d;

    cudaDeviceSynchronize();
    cuda_axpy_<int,int> <<<512,512>>>(
	       n,
	       1,
	       0,
	       aptr,
	       a_d);
    cudaDeviceSynchronize();

    return ;
  }

  void cudamemcpy_dptr2d_double_(
			  double    *aptr,
			  devicemem *a_d_,
			  PARCEL_INT *n_
			 ){
    PARCEL_INT n = n_[0];
    double *a_d =(double *)a_d_[0].a_d;

    cudaDeviceSynchronize();
    cuda_axpy_<double,double> <<<512,512>>>(
	       n,
	       1.0,
	       0.0,
	       aptr,
	       a_d);
    //    cudaDeviceSynchronize();
    //    cuda_ptr_type_view_<double> <<<1,1>>>(10,a_d);
    cudaDeviceSynchronize();

    return ;
  }

  void cudamemcpy_d2dptr_double_(
			  devicemem *a_d_,
			  double    *aptr,
			  PARCEL_INT *n_
			 ){
    PARCEL_INT n = n_[0];
    double *a_d =(double *)a_d_[0].a_d;

    cudaDeviceSynchronize();
    cuda_axpy_<double,double> <<<512,512>>>(
	       n,
	       1.0,
	       0.0,
	       a_d,
	       aptr
					    );
    cudaDeviceSynchronize();

    return ;
  }


#if 0
  void cudamemcpy_dptr2d_(
			  void *aptr,
			  devicemem *a_d_,
			  int typeSize,
			  PARCEL_INT n
			 ){
    const size_t len = typeSize * n;
    void *a_d =(void *)a_d_[0].a_d;

    cudaDeviceSynchronize();
    printf(" ---------- \n");
    cuda_ptr_type_view_<double> <<<1,1>>>(10,(double *)aptr);
    printf(" ---------- \n");
    cudaDeviceSynchronize();

    cudaMemcpy( a_d, aptr,  len, cudaMemcpyDeviceToDevice);

    cudaDeviceSynchronize();
    cuda_ptr_type_view_<double> <<<1,1>>>(10,(double *)a_d);
#if 1
    cudaDeviceSynchronize();
    cuda_axpy_<double,double> <<<1,1>>>(
	       10,
	       1.0,
	       0.0,
	       (double *)aptr,
	       (double *)a_d);
    cudaDeviceSynchronize();
    cuda_ptr_type_view_<double> <<<1,1>>>(10,(double *)a_d);
    cudaDeviceSynchronize();
#endif

    printf(" ---------- \n");

    return ;
  }

  void cudamemcpy_d2dptr_(
			  devicemem *a_d_,
			  void *aptr,
			  int typeSize,
			  PARCEL_INT n
			 ){
    const size_t len = typeSize * n;
    void *a_d =(void *)a_d_[0].a_d;
    cudaMemcpy( aptr, a_d,  len, cudaMemcpyDeviceToDevice);
    return ;
  }
#endif

  void cudaMemcpy_h2d_FP64(
			   double *a_h,
			   devicemem *a_d_,
			   PARCEL_INT n,
			   PARCEL_INT offset_h, 
			   PARCEL_INT offset_d ){
    //    int n = n_[0];
    double *a_d =(double *)a_d_[0].a_d;
    const size_t len = sizeof(double) * n;
    //    printf("cudaMemcpy_h2d_FP64 : %f \n",a_h[offset_h]);
    cudaMemcpy( &a_d[offset_d], &a_h[offset_h],  len, cudaMemcpyHostToDevice);
    return ;
  }

  void cudaTransfer_h2d_FP64toBF16(double *a_h_FP64,devicemem *a_d,PARCEL_INT n, int mode){
    devicemem a_d_FP64;
    cudaMalloc_FP64(&a_d_FP64,n);
    cudaMemcpy_h2d_FP64(a_h_FP64,&a_d_FP64,n,0,0);
    cuda_FP64toBF16_<<<512,512>>>(
				   n,
				   (double *)a_d_FP64.a_d,
				   (__nv_bfloat16 *)a_d[0].a_d,
				   mode);
    cudaFree_(&a_d_FP64);
    return ;
  }

  void cudaMalloc_int(devicemem *a_d,PARCEL_INT n){
    const size_t len = sizeof(int)*n;
    cudaError_t err_code;
    err_code = cudaMalloc(&(a_d[0].a_d), len);
#if 0
    if (err_code != cudaSuccess) {
      printf("error : get cuda Memory size = %ld : %d , %d \n",len,err_code,cudaSuccess);
      exit(-1);
    }
#else
    cuda_malloc_check_err(err_code);
#endif

    return ;
  }

  void cudaMemcpy_d2h_int(
			  devicemem *a_d_,int *a_h,PARCEL_INT n,
			  PARCEL_INT offset_d, 
			  PARCEL_INT offset_h ){
    int *a_d = (int *)a_d_[0].a_d;
    const size_t len = sizeof(int) * n;
    cudaMemcpy( &a_h[offset_h], &a_d[offset_d], len, cudaMemcpyDeviceToHost);
    return ;
  }

  void cudaMemcpy_h2d_int(
			  int *a_h,devicemem *a_d_,PARCEL_INT n,
			  PARCEL_INT offset_h, 
			  PARCEL_INT offset_d ){

    int* a_d = (int *)a_d_[0].a_d;
    //    int n = n_[0];
    const size_t len = sizeof(int) * n;
    cudaMemcpy(&a_d[offset_d], &a_h[offset_h],  len, cudaMemcpyHostToDevice);
    return ;
  }

  void cudaMalloc_PARCEL_INT(devicemem *a_d,PARCEL_INT n){
    const size_t len = sizeof(PARCEL_INT)*n;
    cudaError_t err_code;
    err_code = cudaMalloc(&(a_d[0].a_d), len);
#if 0
    if (err_code != cudaSuccess) {
      printf("error : get cuda Memory size = %ld : %d , %d \n",len,err_code,cudaSuccess);
      exit(-1);
    }
#else
    cuda_malloc_check_err(err_code);
#endif

    return ;
  }

  void cudaMemcpy_d2h_PARCEL_INT(
				 devicemem *a_d_,PARCEL_INT *a_h,PARCEL_INT n,
				 PARCEL_INT offset_d, 
				 PARCEL_INT offset_h ){
    PARCEL_INT *a_d = (PARCEL_INT *)a_d_[0].a_d;
    const size_t len = sizeof(PARCEL_INT) * n;
    cudaMemcpy( &a_h[offset_h], &a_d[offset_d], len, cudaMemcpyDeviceToHost);
    return ;
  }

  void cudaMemcpy_h2d_PARCEL_INT(
				 PARCEL_INT *a_h,devicemem *a_d_,PARCEL_INT n,
				 PARCEL_INT offset_h, 
				 PARCEL_INT offset_d ){

    //    PARCEL_INT n = n_[0];
    const size_t len = sizeof(PARCEL_INT) * n;
    PARCEL_INT *a_d = (PARCEL_INT *)a_d_[0].a_d;

    cudaMemcpy( &a_d[offset_d], &a_h[offset_h],  len, cudaMemcpyHostToDevice);

    // cuda_ptr_type_view_<PARCEL_INT> <<<1,1>>>(10,a_d);
    // for(PARCEL_INT i=0;i<10;i++){
    //   printf(" d2h -- %d, %d \n",i,a_h);
    // }

    return ;

  }

  // MPI 関連
  void malloc_MPI_Status(hostmem *a_d,PARCEL_INT n){
    a_d[0].a_d = (MPI_Status*)malloc(sizeof(MPI_Status)*n);
    return ;
  }

  void malloc_MPI_Request(hostmem *a_d,PARCEL_INT n){
    //  MPI_Status   stat_send[6],stat_recv[6];
    //  MPI_Request req_send[6];
    //  MPI_Comm comm;
    //  MPI_Fint comm_f;
    //     err_code = cudamalloc(&(a_d[0].a_d), len);
    //    type *A = (type*)malloc(sizeof(type)*m_*7);
    a_d[0].a_d = (MPI_Request*)malloc(sizeof(MPI_Request)*n);
    return ;
  }
  void malloc_MPI_Comm(hostmem *a_d,PARCEL_INT n){
    a_d[0].a_d = (MPI_Comm*)malloc(sizeof(MPI_Comm)*n);
    return ;
  }

  void C_free(hostmem *a_d){
    free(a_d[0].a_d);
    return ;
  }

  void C_MPI_Comm_f2c(hostmem *a_d,MPI_Fint f_comm){
#if 0
    MPI_Comm *comm = (MPI_Comm*)a_d[0].a_d;
    comm[0] = MPI_Comm_f2c(f_comm);
#else
    MPI_Comm comm = MPI_Comm_f2c(f_comm);
    memcpy(a_d[0].a_d,&comm,sizeof(MPI_Comm));
#endif
    return ;
  }

  void C_MPI_Isend_FP64(
			hostmem *a_d,
			int ndata,
			int send_rank,
			int tag,
			hostmem *comm_,
			hostmem *req_,
			int id,int offset
			){
    int ierr;
    double *sbuf     = (double *)a_d[0].a_d;
    MPI_Comm *comm   = (MPI_Comm *)comm_[0].a_d;
    MPI_Request *req = (MPI_Request *)req_[id-1].a_d;
#if 1
    ierr = MPI_Isend(
	      &sbuf[offset-1],
	      ndata,
	      MPI_DOUBLE,
	      send_rank,
	      tag,
	      comm[0],
	      //MPI_COMM_WORLD,
	      req
	      );
#else
    ierr = MPI_Isend(
	      &sbuf[offset-1],
	      ndata,
	      MPI_DOUBLE,
	      send_rank,
	      tag,
	      MPI_COMM_WORLD,
	      req
	      );
#endif
    return ;
  }

  void C_MPI_Irecv_FP64(
			hostmem *a_d,
			int ndata,
			int recv_rank,
			int tag,
			hostmem *comm_,
			hostmem *req_,
			int id,int offset
			){
    int ierr;
    double *rbuf     = (double *)a_d[0].a_d;
    MPI_Comm *comm   = (MPI_Comm *)comm_[0].a_d;
    MPI_Request *req = (MPI_Request *)req_[id-1].a_d;
#if 1
    ierr = MPI_Irecv(
	      &rbuf[offset-1],
	      ndata,
	      MPI_DOUBLE,
	      recv_rank,
	      tag,
	      comm[0],
	      //MPI_COMM_WORLD,
	      req);
#else
    ierr = MPI_Irecv(
	      &rbuf[offset-1],
	      ndata,
	      MPI_DOUBLE,
	      recv_rank,
	      tag,
	      MPI_COMM_WORLD,
	      req);
#endif
    return ;
  }

  void C_MPI_Isend_FP64buf_FP64(
			double *a_d,
			int ndata,
			int send_rank,
			int tag,
			hostmem *comm_,
			hostmem *req_,
			int id,int offset
			){
  int ierr;
  double *sbuf     = a_d;
  MPI_Comm *comm   = (MPI_Comm *)comm_[0].a_d;
  MPI_Request *req = (MPI_Request *)req_[id-1].a_d;

  ierr = MPI_Isend(
	      &sbuf[offset-1],
	      ndata,
	      MPI_DOUBLE,
	      send_rank,
	      tag,
	      comm[0],
	      //	      MPI_COMM_WORLD,
	      req
	      );

    return ;
  }

  void C_MPI_Irecv_FP64buf_FP64(
			double *a_d,
			int ndata,
			int recv_rank,
			int tag,
			hostmem *comm_,
			hostmem *req_,
			int id,int offset
			){
    int ierr;
    double *rbuf     = a_d;
    MPI_Comm *comm   = (MPI_Comm *)comm_[0].a_d;
    MPI_Request *req = (MPI_Request *)req_[id-1].a_d;
    ierr = MPI_Irecv(
	      &rbuf[offset-1],
	      ndata,
	      MPI_DOUBLE,
	      recv_rank,
	      tag,
		comm[0],
		//MPI_COMM_WORLD,
	      req);

    return ;
  }

  void C_MPI_Wait(hostmem *req_,hostmem *stat_,int id){
    int ierr;
    MPI_Request *req = (MPI_Request *)req_[id-1].a_d; 
    MPI_Status *stat = (MPI_Status *)stat_[id-1].a_d;
    ierr = MPI_Wait(req,stat);
    return ;
  }
}
