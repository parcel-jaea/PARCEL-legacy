#include <cuda.h>
#include <cublas_v2.h>
#include <cusparse.h>

#ifndef HEADER_CUCOMMMON_H
#define HEADER_CUCOMMMON_H

// const int warp_size = 32

#define WARP_SIZE 32

#include <mpi.h>

#ifdef PARCEL_INT8
#define PARCEL_INT int64_t
#else
#define PARCEL_INT int
#endif

#include <cuda_bf16.h>
#include <stdint.h>
// define.f90と一致
#define KSP_PRECONDITION_NONE 0
#define KSP_PRECONDITION_POINTJACOBI 1
#define KSP_PRECONDITION_BLOCKJACOBI 2
#define KSP_PRECONDITION_ADDITIVESCHWARZ 3 
#define KSP_PRECONDITION_NEUMANN 4 
// -----

#define FP64TOBF16_ROUND_NEAR  0
#define FP64TOBF16_ROUND_ZERO  1
#define FP64TOBF16_ROUND_DOWN  2
#define FP64TOBF16_ROUND_UP    3
#define FP64TOBF16_ROUND_ABSINF 4

// 構造体はFortran_側と整合するように定義す
typedef struct Fortran_GPUMem{
  //  double *a_d;
  void *a_d;
  void *descr;
  void *descr_add;
} devicemem;

// hostmemは定義しているがFortran側ではdevicememで定義しているため,
// devicememと同じ記述とする.
typedef struct Fortran_CPUMem{
  //  double *a_d;
  void *a_d;
  void *descr;
  void *descr_add;
} hostmem;

/* typedef struct Fortran_cublasHandle{ */
/*   //  double *a_d; */
/*   void *handle; */
/* } cublasHandle_; */

//c_common.f90 : cublas_typeと一致させる
typedef struct Fortran_cublas{
  void *cublasHandle; // cublasHandle_t
  devicemem *alpha;
  devicemem *beta;
} f_cublas;

//c_common.f90 : cusparce_typeと一致させる
typedef struct Fortran_cusparce{
  void *cublasHandle; // cublasHandle_t
  void *cusparseHandle; // cusparseHandle_t
  
  void *cusparseMatDescr; // cusparseMatDescr_t
  void *matA; // cusparseSpMatDescr_t
  
  // ------
  // BCSR形式
  int blockDim;
  
  int mb;
  int nb;
  int nnzb; 
  void *cusparseMatDescr_bsr;  // cusparseMatDescr_t
  void *bsrRowPtr; // int
  void *bsrColInd; // int
  void *bsrVal; // double

  int ipreflag;
  int precon_mb;
  int precon_nb;
  int precon_nnzb; 
  void *cusparseMatDescr_precon_bsr;  // cusparseMatDescr_t
  void *precon_bsrRowPtr; // int
  void *precon_bsrColInd; // int  
  void *precon_bsrVal; // double
  // ------
  
  void *cusparseMatDescr_precon; // cusparseMatDescr_t
  void *precon_matA; // cusparseSpMatDescr_t

  void *infoILU; // csrilu02Info_t
  void *descrL; // cusparseMatDescr_t
  void *infoL; // csrsv2Info_t
  void *descrU; // cusparseMatDescr_t 
  void *infoU; // csrsv2Info_t

  // ------
  // BCSR形式  
  void *descrILU_bcsr; // cusparseMatDescr_t @ : 廃止
  void *infoILU_bcsr;  // bsrilu02Info_t     @ : 廃止
  void *descrL_bcsr;   // cusparseMatDescr_t @ : 廃止
  void *infoL_bcsr;    // bsrsv2Info_t       @ : 廃止
  void *descrU_bcsr;   // cusparseMatDescr_t @ : 廃止
  void *infoU_bcsr;    // bsrsv2Info_t       @ : 廃止
  // ------  

  void *buffer; //  double,float
  void *precon_buffer; //  double,float

} f_cusparse;

__inline__ void cuda_malloc_check_err(cudaError_t err_code){
  if (err_code != cudaSuccess) {
    printf("%s - %d : error : CUDA_MALLOC  %d , %d \n",__FILE__,__LINE__,err_code,cudaSuccess);
    exit(-1);
  }
}

extern "C"{
  void cuda_dot(
		 MPI_Fint f_comm,
		 int numGPUblock,
		 int numGPUthread,
		 PARCEL_INT n,
		 devicemem *val_,PARCEL_INT offset_val,
		 devicemem *x_,PARCEL_INT offset_x,
		 devicemem *y_,PARCEL_INT offset_y
		);
  void cuda_dot_local(
		 int numGPUblock,
		 int numGPUthread,
		 PARCEL_INT n,
		 devicemem *val , PARCEL_INT offset_val,
		 devicemem *x_  , PARCEL_INT offset_x,
		 devicemem *y_  , PARCEL_INT offset_y
		      );
}

extern cublasHandle_t   PARCEL_CUBLAS_HANDLE;
extern cusparseHandle_t PARCEL_CUSPARSE_SPMV_HANDLE;
extern cusparseHandle_t PARCEL_CUSPARSE_PRECON_HANDLE;

#if 0
typedef struct PARCEL_INIT_MATRIX{
  int matrixflag;
  int machineflag;
  // crs
  double *crsA;
  PARCEL_INT *;
  PARCEL_INT *;

  // dia
  double *diaA;
  PARCEL_INT *offset;
  PARCEL_INT nnd;

} PARCEL_MATRIX;
extern PARCEL_MATRIX PARCEL_Mat;
#endif

// extern cusparseHandle_t PARCEL_CUSPARSE_SPMV_HANDLE;
// extern cusparseHandle_t PARCEL_CUSPARSE_PRECON_HANDLE;

#endif

