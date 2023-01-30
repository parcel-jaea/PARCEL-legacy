#define USE_BCSR_MATRIX 0
// #define USE_BCSR_MATRIX 1
#include <stdio.h>

#include <cuda_runtime.h>
// #include "cuda_errcheck.h"
#include "cu_common.h"

#include <cusparse.h>
#include <cublas_v2.h>
#include <mpi.h>
 
// ILUに関する処理は前処理フラグで変更する
cusparseHandle_t PARCEL_CUSPARSE_SPMV_HANDLE;
cusparseHandle_t PARCEL_CUSPARSE_PRECON_HANDLE;

template <class Tdata,class Tcalc>
void cusparse_init_precon_ILU0_(
			   f_cusparse *f_cusparse_,
			   cublasHandle_t *cublasHandle,
			   cusparseHandle_t *cusparseHandle,
			   PARCEL_INT N,PARCEL_INT nCol,PARCEL_INT nz,
			   // 前処理用のCRS行列
			   cusparseMatDescr_t *precon_descr,
			   cusparseSpMatDescr_t *precon_matA,
			   Tdata *d_precon_val, // ILU(0)の行列 device memory
			   int *d_precon_col, // device memory
			   int *d_precon_row, // device memory
			   PARCEL_INT precon_N,PARCEL_INT precon_nCol,PARCEL_INT precon_nz,
			   csrilu02Info_t *infoILU,
			   cusparseMatDescr_t *descrL, csrsv2Info_t *infoL,		    
			   cusparseMatDescr_t *descrU, csrsv2Info_t *infoU,
#if USE_BCSR_MATRIX
			   // ------
			   // BCSR形式  
			   cusparseMatDescr_t *descrILU_bcsr, bsrilu02Info_t *infoILU_bcsr,
			   cusparseMatDescr_t *descrL_bcsr  , bsrsv2Info_t *infoL_bcsr,
			   cusparseMatDescr_t *descrU_bcsr  , bsrsv2Info_t *infoU_bcsr,
			   // ------  		    
#endif
			   size_t *bufferSize_ // SpMV,ILU0 共通のバッファサイズ
			   ){

  const cudaDataType_t cudaDataType = CUDA_R_64F;

#ifdef PARCEL_INT8_
  const cusparseIndexType_t cusparseIndexType= CUSPARSE_INDEX_64I;
#else
  const cusparseIndexType_t cusparseIndexType= CUSPARSE_INDEX_32I;
#endif

  /* Description of the preconditioning A matrix */
  cusparseCreateMatDescr(precon_descr);
  cusparseSetMatType(precon_descr[0], CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(precon_descr[0], CUSPARSE_INDEX_BASE_ZERO);

  /* Wrap raw data into cuSPARSE generic API objects */
  cusparseCreateCsr(
		    precon_matA,
		    precon_N, precon_nCol, precon_nz,
		    d_precon_row, d_precon_col, d_precon_val,
		    cusparseIndexType,
		    cusparseIndexType,
		    CUSPARSE_INDEX_BASE_ZERO,
		    cudaDataType);

  cudaError_t err_code;

  Tdata *d_p;
  Tdata *d_omega;
  err_code=cudaMalloc((void **)&d_p, N*sizeof(Tdata));cuda_malloc_check_err(err_code);
  err_code=cudaMalloc((void **)&d_omega, N*sizeof(Tdata));cuda_malloc_check_err(err_code);

  cusparseDnVecDescr_t vecp;
  cusparseCreateDnVec(&vecp, nCol, d_p, cudaDataType);
  cusparseDnVecDescr_t vecomega;
  cusparseCreateDnVec(&vecomega, N, d_omega, cudaDataType);

  /* Create ILU(0) info object */
  cusparseCreateCsrilu02Info(infoILU);

  /* Create L factor descriptor and triangular solve info */
  cusparseCreateMatDescr(descrL);
  cusparseSetMatType(descrL[0], CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(descrL[0], CUSPARSE_INDEX_BASE_ZERO);
  cusparseSetMatFillMode(descrL[0], CUSPARSE_FILL_MODE_LOWER);

  // ILU分解のLの対角要素を1に指定する
  cusparseSetMatDiagType(descrL[0], CUSPARSE_DIAG_TYPE_UNIT);
  //  cusparseSetMatDiagType(descrL[0], CUSPARSE_DIAG_TYPE_NON_UNIT);

  cusparseCreateCsrsv2Info(infoL);

  /* Create U factor descriptor and triangular solve info */
  cusparseCreateMatDescr(descrU);
  cusparseSetMatType(descrU[0], CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(descrU[0], CUSPARSE_INDEX_BASE_ZERO);
  cusparseSetMatFillMode(descrU[0], CUSPARSE_FILL_MODE_UPPER);
  cusparseSetMatDiagType(descrU[0], CUSPARSE_DIAG_TYPE_NON_UNIT);
					 
  cusparseCreateCsrsv2Info(infoU);

  /* Allocate workspace for cuSPARSE */
  size_t bufferSize = bufferSize_[0];
  //  size_t bufferSize = 0;
  //size_t tmp = 1;
  int stmp = 1;
  void *buffer;
  //  const Tdata Tdataone  = 1.0;
  //  const Tdata Tdatazero = 0.0;

  cusparseDcsrilu02_bufferSize(
			       // cusparseHandle[0],
			       PARCEL_CUSPARSE_PRECON_HANDLE,
			       precon_N, precon_nz, precon_descr[0], 
			       d_precon_val, d_precon_row, d_precon_col, 
			       infoILU[0], &stmp);

  if (stmp > bufferSize) {
    bufferSize = stmp;
  }
  cusparseDcsrsv2_bufferSize(
			     // cusparseHandle[0],
			     PARCEL_CUSPARSE_PRECON_HANDLE,
			     CUSPARSE_OPERATION_NON_TRANSPOSE, 
			     precon_N, precon_nz, descrL[0], d_precon_val,
			     d_precon_row, d_precon_col, 
			     infoL[0], &stmp);
  if (stmp > bufferSize) {
    bufferSize = stmp;
  }
  cusparseDcsrsv2_bufferSize(
			     // cusparseHandle[0],
			     PARCEL_CUSPARSE_PRECON_HANDLE,
			     CUSPARSE_OPERATION_NON_TRANSPOSE, 
			     precon_N, precon_nz, descrU[0], 
			     d_precon_val,d_precon_row, d_precon_col, 
			     infoU[0], &stmp);

  if (stmp > bufferSize) {
    bufferSize = stmp;
  }
  if(bufferSize_[0] < bufferSize){
    bufferSize_[0] = bufferSize;
  }
  err_code=cudaMalloc(&buffer, bufferSize);cuda_malloc_check_err(err_code);

  /* Perform analysis for ILU(0) */
  cusparseDcsrilu02_analysis(
			     //			     cusparseHandle[0], 
			     PARCEL_CUSPARSE_PRECON_HANDLE,
			     precon_N, precon_nz, precon_descr[0], 
			     d_precon_val, d_precon_row, d_precon_col,
			     infoILU[0],
			     CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer);

  /* generate the ILU(0) factors */
  cusparseDcsrilu02(
		    //		    cusparseHandle[0], 
		    PARCEL_CUSPARSE_PRECON_HANDLE,
		    precon_N, precon_nz, precon_descr[0], 
		    d_precon_val, d_precon_row, d_precon_col, infoILU[0],
		    CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer);

  /* perform triangular solve analysis */
  cusparseDcsrsv2_analysis(
			   // cusparseHandle[0],
			     PARCEL_CUSPARSE_PRECON_HANDLE,
			   CUSPARSE_OPERATION_NON_TRANSPOSE, 
			   precon_N, precon_nz, descrL[0],
			   d_precon_val, d_precon_row, d_precon_col, 
			   infoL[0], CUSPARSE_SOLVE_POLICY_USE_LEVEL,
			   buffer);
  cusparseDcsrsv2_analysis(
			   //cusparseHandle[0],
			     PARCEL_CUSPARSE_PRECON_HANDLE,
			   CUSPARSE_OPERATION_NON_TRANSPOSE, 
			   precon_N, precon_nz, descrU[0],
			   d_precon_val, d_precon_row, d_precon_col, 
			   infoU[0], CUSPARSE_SOLVE_POLICY_USE_LEVEL,
			   buffer);
  
  
  cusparseDestroyDnVec(vecp);
  cusparseDestroyDnVec(vecomega);
  cudaFree(buffer);
  cudaFree(d_p);
  cudaFree(d_omega);
  if(bufferSize_[0]<bufferSize){
    bufferSize_[0] = bufferSize;
    cudaFree(f_cusparse_[0].buffer);
    //    cudaMalloc(  (void **)&(f_cusparse_[0].buffer), bufferSize_[0] );    
  }


  return ;

}

template <class Tdata,class Tcalc>
void cusparse_init_SpMV_(
		    f_cusparse *f_cusparse_,
		    cublasHandle_t *cublasHandle,
		    cusparseHandle_t *cusparseHandle,
		    // SpMV用のCRS行列
		    cusparseMatDescr_t *descr,
		    cusparseSpMatDescr_t *matA,
		    Tdata *d_val, // device memory
		    int *d_col, // device memory
		    int *d_row,  // device memory
		    PARCEL_INT N,PARCEL_INT nCol,PARCEL_INT nz,
		    size_t *bufferSize_ // SpMV,ILU0 共通のバッファサイズ
		    ){
  const cudaDataType_t cudaDataType = CUDA_R_64F;

#ifdef PARCEL_INT8_
  const cusparseIndexType_t cusparseIndexType= CUSPARSE_INDEX_64I;
#else
  const cusparseIndexType_t cusparseIndexType= CUSPARSE_INDEX_32I;
#endif

  double stime,etime;

  /* Create CUBLAS context */
  //  cublasCreate(cublasHandle);

  /* Create CUSPARSE context */
  //  cusparseCreate(cusparseHandle);

  /* Description of the A matrix */
  stime=MPI_Wtime();
  cusparseCreateMatDescr(descr);
  cusparseSetMatType(descr[0], CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(descr[0], CUSPARSE_INDEX_BASE_ZERO);
  etime=MPI_Wtime();
  // printf("     cusparseCreateMatDescr time : %f \n",etime - stime);
  stime=MPI_Wtime();

  /* Wrap raw data into cuSPARSE generic API objects */
  stime=MPI_Wtime();
  cusparseCreateCsr(
		    matA,
		    N, nCol, nz,
		    d_row, d_col, d_val,
		    cusparseIndexType,
		    cusparseIndexType,
		    CUSPARSE_INDEX_BASE_ZERO,
		    cudaDataType);
  etime=MPI_Wtime();
  // printf("     cusparseCreateCsr time : %f \n",etime - stime);

  stime=MPI_Wtime();
  cudaError_t err_code;
  Tdata *d_p;
  Tdata *d_omega;
  err_code=cudaMalloc((void **)&d_p, N*sizeof(Tdata));cuda_malloc_check_err(err_code);
  err_code=cudaMalloc((void **)&d_omega, N*sizeof(Tdata));cuda_malloc_check_err(err_code);

  cusparseDnVecDescr_t vecp;
  cusparseCreateDnVec(&vecp, nCol, d_p, cudaDataType);
  cusparseDnVecDescr_t vecomega;
  cusparseCreateDnVec(&vecomega, N, d_omega, cudaDataType);

  /* Allocate workspace for cuSPARSE */
  size_t bufferSize = 0;
  size_t tmp = 1;
  //  int stmp = 0;
  const Tdata Tdataone  = 1.0;
  const Tdata Tdatazero = 0.0;
  
  cusparseSpMV_bufferSize(
			  //			  cusparseHandle[0],
			  PARCEL_CUSPARSE_SPMV_HANDLE,
			  CUSPARSE_OPERATION_NON_TRANSPOSE,
			  &Tdataone, matA[0],
			  vecp, &Tdatazero, vecomega, 
			  cudaDataType, CUSPARSE_MV_ALG_DEFAULT,
			  &tmp);
  
  if (tmp > bufferSize) {
    bufferSize = tmp;
  }
  bufferSize_[0] = bufferSize;

  //  cudaMalloc(  (void **)&(f_cusparse_[0].buffer), bufferSize_[0] );    
  
  cusparseDestroyDnVec(vecp);
  cusparseDestroyDnVec(vecomega);
  cudaFree(d_p);
  cudaFree(d_omega);

  etime=MPI_Wtime();
  // printf("     cusparseSpMV_bufferSize time : %f \n",etime - stime);

  return ;

}

template <class Tdata,class Tcalc>
void cusparse_init_Neumann_(
		    f_cusparse *f_cusparse_,
		    cublasHandle_t *cublasHandle,
		    cusparseHandle_t *cusparseHandle,
		    // SpMV用のCRS行列
		    cusparseMatDescr_t *descr,
		    cusparseSpMatDescr_t *matA,
		    Tdata *d_val, // device memory
		    int *d_col, // device memory
		    int *d_row,  // device memory
		    PARCEL_INT N,PARCEL_INT nCol,PARCEL_INT nz,
		    size_t *bufferSize_ // SpMV,ILU0 共通のバッファサイズ
		    ){
  const cudaDataType_t cudaDataType = CUDA_R_64F;

#ifdef PARCEL_INT8_
  const cusparseIndexType_t cusparseIndexType= CUSPARSE_INDEX_64I;
#else
  const cusparseIndexType_t cusparseIndexType= CUSPARSE_INDEX_32I;
#endif

  /* Description of the A matrix */
  cusparseCreateMatDescr(descr);
  cusparseSetMatType(descr[0], CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(descr[0], CUSPARSE_INDEX_BASE_ZERO);

  /* Wrap raw data into cuSPARSE generic API objects */
  cusparseCreateCsr(
		    matA,
		    N, nCol, nz,
		    d_row, d_col, d_val,
		    cusparseIndexType,
		    cusparseIndexType,
		    CUSPARSE_INDEX_BASE_ZERO,
		    cudaDataType);

  Tdata *d_p;
  Tdata *d_omega;
  cudaError_t err_code;
  err_code=cudaMalloc((void **)&d_p, N*sizeof(Tdata));cuda_malloc_check_err(err_code);
  err_code=cudaMalloc((void **)&d_omega, N*sizeof(Tdata));cuda_malloc_check_err(err_code);

  cusparseDnVecDescr_t vecp;
  cusparseCreateDnVec(&vecp, nCol, d_p, cudaDataType);
  cusparseDnVecDescr_t vecomega;
  cusparseCreateDnVec(&vecomega, N, d_omega, cudaDataType);

  /* Allocate workspace for cuSPARSE */
  size_t bufferSize = bufferSize_[0];
  size_t tmp = 1;
  //  void *buffer;
  const Tdata Tdataone  = 1.0;
  const Tdata Tdatazero = 0.0;
  
  cusparseSpMV_bufferSize(
			  // cusparseHandle[0],
			     PARCEL_CUSPARSE_PRECON_HANDLE,
			  CUSPARSE_OPERATION_NON_TRANSPOSE,
			  &Tdataone, matA[0],
			  vecp, &Tdatazero, vecomega, 
			  cudaDataType, CUSPARSE_MV_ALG_DEFAULT,
			  &tmp);
  
  if (tmp > bufferSize) {
    bufferSize = tmp;
    bufferSize_[0] = bufferSize;
  }
  
  cusparseDestroyDnVec(vecp);
  cusparseDestroyDnVec(vecomega);
  cudaFree(d_p);
  cudaFree(d_omega);

  return ;

}


template <class Tdata,class Tcalc>
void cusparse_init_(
		    f_cusparse *f_cusparse_,
		    cublasHandle_t *cublasHandle,
		    cusparseHandle_t *cusparseHandle,
		    // SpMV用のCRS行列
		    cusparseMatDescr_t *descr,
		    cusparseSpMatDescr_t *matA,
		    Tdata *d_val, // device memory
		    int *d_col, // device memory
		    int *d_row,  // device memory
		    PARCEL_INT N,PARCEL_INT nCol,PARCEL_INT nz,
		    // 前処理用のCRS行列
		    cusparseMatDescr_t *precon_descr,
		    cusparseSpMatDescr_t *precon_matA,
		    Tdata *d_precon_val, // ILU(0)の行列 device memory
		    int *d_precon_col, // device memory
		    int *d_precon_row, // device memory
		    PARCEL_INT precon_N,PARCEL_INT precon_nCol,PARCEL_INT precon_nz,
		    csrilu02Info_t *infoILU,
		    cusparseMatDescr_t *descrL, csrsv2Info_t *infoL,		    
		    cusparseMatDescr_t *descrU, csrsv2Info_t *infoU,
		    // ------
		    // BCSR形式  
#if USE_BCSR_MATRIX
		    cusparseMatDescr_t *descrILU_bcsr, bsrilu02Info_t *infoILU_bcsr,
		    cusparseMatDescr_t *descrL_bcsr  , bsrsv2Info_t *infoL_bcsr,
		    cusparseMatDescr_t *descrU_bcsr  , bsrsv2Info_t *infoU_bcsr,
#endif
		    // ------  		    
		    size_t *bufferSize_ // SpMV,ILU0 共通のバッファサイズ
		    ){
  bufferSize_[0] = 1;
  cusparse_init_SpMV_<Tdata,Tcalc>(
		      f_cusparse_,
		      cublasHandle,
		      cusparseHandle,
		      // SpMV用のCRS行列
		      descr,
		      matA,
		      d_val, // device memory
		      d_col, // device memory
		      d_row,  // device memory
		      N, nCol, nz,
		      bufferSize_ // SpMV,ILU0 共通のバッファサイズ
		      );
  const int ipreflag = f_cusparse_[0].ipreflag;
  if( (ipreflag == KSP_PRECONDITION_BLOCKJACOBI) || 
      (ipreflag == KSP_PRECONDITION_ADDITIVESCHWARZ) ){
    cusparse_init_precon_ILU0_<Tdata,Tcalc>(
					    f_cusparse_,
					    cublasHandle,
					    cusparseHandle,
					    N, nCol, nz,
					    // 前処理用のCRS行列
					    precon_descr,
					    precon_matA,
					    d_precon_val, // ILU(0)の行列 device memory
					    d_precon_col, // device memory
					    d_precon_row, // device memory
					    precon_N, precon_nCol, precon_nz,
					    infoILU,
					    descrL, infoL,
					    descrU, infoU,
#if USE_BCSR_MATRIX
					    // ------
					    // BCSR形式  
					    descrILU_bcsr, infoILU_bcsr,
					    descrL_bcsr  , infoL_bcsr,
					    descrU_bcsr  , infoU_bcsr,
					    // ------  		    
#endif
					    bufferSize_ // SpMV,ILU0 共通のバッファサイズ
					    );
  }

  if(ipreflag == KSP_PRECONDITION_NEUMANN){
    cusparse_init_Neumann_<Tdata,Tcalc>(
					f_cusparse_,
					cublasHandle,
					cusparseHandle,
					precon_descr,
					precon_matA,
					d_precon_val, // ILU(0)の行列 device memory
					d_precon_col, // device memory
					d_precon_row, // device memory
					precon_N, precon_nCol, precon_nz,
					bufferSize_ // SpMV,ILU0 共通のバッファサイズ
					);
  }

  return ;

}


extern "C" {
  void cusparse_init(
		     f_cusparse *f_cusparse_,
		     devicemem *d_val_ ,
		     devicemem *d_col_ ,
		     devicemem *d_row_ ,
		     PARCEL_INT nRow,
		     PARCEL_INT nCol,
		     PARCEL_INT nnz,
		     devicemem *d_precon_val_ ,
		     devicemem *d_precon_col_ ,
		     devicemem *d_precon_row_ ,
		     PARCEL_INT precon_nRow,
		     PARCEL_INT precon_nCol,
		     PARCEL_INT precon_nnz,
		     PARCEL_INT blockDim_
		     ){

    f_cusparse_[0].cublasHandle       = (cublasHandle_t *)malloc(sizeof(cublasHandle_t));
    f_cusparse_[0].cusparseHandle     = (cusparseHandle_t *)malloc(sizeof(cusparseHandle_t));

    f_cusparse_[0].cusparseMatDescr   = (cusparseMatDescr_t *)malloc(sizeof(cusparseMatDescr_t));
    f_cusparse_[0].matA               = (cusparseSpMatDescr_t *)malloc(sizeof(cusparseSpMatDescr_t));

    f_cusparse_[0].cusparseMatDescr_precon   = (cusparseMatDescr_t *)malloc(sizeof(cusparseMatDescr_t));
    f_cusparse_[0].precon_matA               = (cusparseSpMatDescr_t *)malloc(sizeof(cusparseSpMatDescr_t));

    f_cusparse_[0].infoILU = (csrilu02Info_t *)malloc(sizeof(csrilu02Info_t));
    f_cusparse_[0].descrL  = (cusparseMatDescr_t *)malloc(sizeof(cusparseMatDescr_t));
    f_cusparse_[0].infoL   = (csrsv2Info_t *)malloc(sizeof(csrsv2Info_t));
    f_cusparse_[0].descrU  = (cusparseMatDescr_t *)malloc(sizeof(cusparseMatDescr_t));
    f_cusparse_[0].infoU   = (csrsv2Info_t *)malloc(sizeof(csrsv2Info_t));
    
    double *d_val     = (double *)d_val_[0].a_d;
    // PARCEL_INT *d_col = (PARCEL_INT *)d_col_[0].a_d;
    // PARCEL_INT *d_row = (PARCEL_INT *)d_row_[0].a_d;
    int *d_col = (int *)d_col_[0].a_d;
    int *d_row = (int *)d_row_[0].a_d;
    double *d_precon_val     = (double *)d_precon_val_[0].a_d;
    // PARCEL_INT *d_precon_col = (PARCEL_INT *)d_precon_col_[0].a_d;
    // PARCEL_INT *d_precon_row = (PARCEL_INT *)d_precon_row_[0].a_d;
    int *d_precon_col = (int *)d_precon_col_[0].a_d;
    int *d_precon_row = (int *)d_precon_row_[0].a_d;

    size_t bufferSize_ = 0; // SpMV,ILU0 共通のバッファサイズ 

    cusparse_init_<double,double>(
				  f_cusparse_,
				  (cublasHandle_t *) f_cusparse_[0].cublasHandle,
				  (cusparseHandle_t *) f_cusparse_[0].cusparseHandle,
				  (cusparseMatDescr_t *) f_cusparse_[0].cusparseMatDescr,
				  (cusparseSpMatDescr_t *) f_cusparse_[0].matA,
				  d_val, // device memory
				  d_col, // device memory
				  d_row, // device memory
				  nRow,nCol,nnz,
				  (cusparseMatDescr_t *) f_cusparse_[0].cusparseMatDescr_precon,
				  (cusparseSpMatDescr_t *) f_cusparse_[0].precon_matA,
				  d_precon_val, // device memory
				  d_precon_col, // device memory
				  d_precon_row, // device memory
				  precon_nRow,precon_nCol,precon_nnz,
				  (csrilu02Info_t *) f_cusparse_[0].infoILU,
				  (cusparseMatDescr_t *) f_cusparse_[0].descrL,
				  (csrsv2Info_t *) f_cusparse_[0].infoL,
				  (cusparseMatDescr_t *) f_cusparse_[0].descrU,
				  (csrsv2Info_t *) f_cusparse_[0].infoU,
				  // ----------
#if USE_BCSR_MATRIX
				  (cusparseMatDescr_t *) f_cusparse_[0].descrILU_bcsr,
				  (bsrilu02Info_t *) f_cusparse_[0].infoILU_bcsr,
				  (cusparseMatDescr_t *) f_cusparse_[0].descrL_bcsr,
				  (bsrsv2Info_t *) f_cusparse_[0].infoL_bcsr,
				  (cusparseMatDescr_t *) f_cusparse_[0].descrU_bcsr,
				  (bsrsv2Info_t *) f_cusparse_[0].infoU_bcsr,
#endif
				  // ----------
				  &bufferSize_
				  );

    //    size_t len = sizeof(double) * nRow;
    //    cudaMalloc( (void **) &(f_cusparse_[0].precon_buffer), len );
#if 1
    cudaError_t err_code;
    err_code=cudaMalloc( (void **) &(f_cusparse_[0].buffer), bufferSize_ ) ; cuda_malloc_check_err(err_code);
    size_t len = sizeof(double) * nRow;
    err_code=cudaMalloc( (void **) &(f_cusparse_[0].precon_buffer), len )  ; cuda_malloc_check_err(err_code);
#endif

  }

  void cusparse_free(
		     f_cusparse *f_cusparse_
		     ){
    cublasHandle_t *cublasHandle  = (cublasHandle_t *)f_cusparse_[0].cublasHandle ; // cublasHandle_t
    cusparseHandle_t *cusparseHandle  = (cusparseHandle_t *)f_cusparse_[0].cusparseHandle ; // cusparseHandle_t
  
    cusparseMatDescr_t *cusparseMatDescr  = (cusparseMatDescr_t *)f_cusparse_[0].cusparseMatDescr ; // cusparseMatDescr_t
    cusparseSpMatDescr_t *matA  = (cusparseSpMatDescr_t *)f_cusparse_[0].matA ; // cusparseSpMatDescr_t
  
    cusparseMatDescr_t *cusparseMatDescr_bsr  = (cusparseMatDescr_t *)f_cusparse_[0].cusparseMatDescr_bsr ;  // cusparseMatDescr_t
    cusparseMatDescr_t *cusparseMatDescr_precon_bsr  = (cusparseMatDescr_t *)f_cusparse_[0].cusparseMatDescr_precon_bsr ;  // cusparseMatDescr_t
  
    cusparseMatDescr_t *cusparseMatDescr_precon  = (cusparseMatDescr_t *)f_cusparse_[0].cusparseMatDescr_precon ; // cusparseMatDescr_t
    cusparseSpMatDescr_t *precon_matA  = (cusparseSpMatDescr_t *)f_cusparse_[0].precon_matA ; // cusparseSpMatDescr_t

    csrilu02Info_t *infoILU  = (csrilu02Info_t *)f_cusparse_[0].infoILU ; // csrilu02Info_t
    cusparseMatDescr_t *descrL  = (cusparseMatDescr_t *)f_cusparse_[0].descrL ; // cusparseMatDescr_t
    csrsv2Info_t *infoL  = (csrsv2Info_t *)f_cusparse_[0].infoL ; // csrsv2Info_t
    cusparseMatDescr_t *descrU  = (cusparseMatDescr_t *)f_cusparse_[0].descrU ; // cusparseMatDescr_t 
    csrsv2Info_t *infoU  = (csrsv2Info_t *)f_cusparse_[0].infoU ; // csrsv2Info_t

    int ipreflag = f_cusparse_[0].ipreflag;
    /* Destroy descriptors */
    /* Destroy contexts */

    // free SpMV
    //    cusparseDestroy(cusparseHandle[0]);
    cusparseDestroySpMat(matA[0]);
    cusparseDestroyMatDescr(cusparseMatDescr[0]);
    free(f_cusparse_[0].cusparseHandle);
    free(f_cusparse_[0].matA);
    free(f_cusparse_[0].cusparseMatDescr);

    // precon ILU
    if( (ipreflag == KSP_PRECONDITION_BLOCKJACOBI)  || 
	(ipreflag == KSP_PRECONDITION_ADDITIVESCHWARZ) ){
      cusparseDestroyMatDescr(cusparseMatDescr_precon[0]);
      free(f_cusparse_[0].cusparseMatDescr_precon_bsr);
      free(f_cusparse_[0].cusparseMatDescr_precon);

      cusparseDestroySpMat(precon_matA[0]);
      cusparseDestroyCsrilu02Info(infoILU[0]);
      cusparseDestroyMatDescr(descrL[0]);
      cusparseDestroyCsrsv2Info(infoL[0]);
      cusparseDestroyMatDescr(descrU[0]);
      cusparseDestroyCsrsv2Info(infoU[0]);
    }

    free(f_cusparse_[0].precon_matA);
    free(f_cusparse_[0].infoILU);
    free(f_cusparse_[0].descrL);
    free(f_cusparse_[0].infoL);
    free(f_cusparse_[0].descrU);
    free(f_cusparse_[0].infoU);
    // -----------

    cudaFree(f_cusparse_[0].buffer);
    
    //    cublasDestroy(cublasHandle[0]);
    free(f_cusparse_[0].cublasHandle);
    
    return ;
  }

  //  void swap_ptr(double **x,double **y);

  // y = Ax
  void cusparse_SpMV(
		     f_cusparse *f_cusparse_,
		     devicemem* vecx_, PARCEL_INT offset_x,
		     devicemem* vecy_, PARCEL_INT offset_y
		     ){

    double *vecx_val_ = (double*)vecx_[0].a_d;
    double *vecy_val_ = (double*)vecy_[0].a_d;
    double *vecx_val = &vecx_val_[offset_x];
    double *vecy_val = &vecy_val_[offset_y];

    cusparseHandle_t *cusparseHandle = (cusparseHandle_t *) f_cusparse_[0].cusparseHandle;
    cusparseSpMatDescr_t *matA = (cusparseSpMatDescr_t *)f_cusparse_[0].matA;
    cusparseDnVecDescr_t *vecx = (cusparseDnVecDescr_t *)vecx_[0].descr_add;
    cusparseDnVecDescr_t *vecy = (cusparseDnVecDescr_t *)vecy_[0].descr;   

    const cudaDataType_t cudaDataType = CUDA_R_64F;
    const double doubleone  = 1.0;
    const double doublezero = 0.0;

    cusparseSpMV(
		 //		 cusparseHandle[0],
		 PARCEL_CUSPARSE_SPMV_HANDLE,
		 CUSPARSE_OPERATION_NON_TRANSPOSE,
		 &doubleone,
		 matA[0], // desc
		 vecx[0], // desc
		 &doublezero,
		 vecy[0], // desc
		 cudaDataType,
		 CUSPARSE_MV_ALG_DEFAULT,
		 f_cusparse_[0].buffer);

    // cudaDeviceSynchronize();
    // vecx_[0].a_d = vecy_val_ ;
    // vecy_[0].a_d = vecx_val_ ;
    // cudaDeviceSynchronize();

    return ;
  }

  // y = Ax
  void Neumann_SpMV(
		     f_cusparse *f_cusparse_,
		     devicemem* vecx_, PARCEL_INT offset_x,
		     devicemem* vecy_, PARCEL_INT offset_y
		     ){

    double *vecx_val_ = (double*)vecx_[0].a_d;
    double *vecy_val_ = (double*)vecy_[0].a_d;
    double *vecx_val = &vecx_val_[offset_x];
    double *vecy_val = &vecy_val_[offset_y];

    cusparseHandle_t *cusparseHandle = (cusparseHandle_t *) f_cusparse_[0].cusparseHandle;
    cusparseSpMatDescr_t *matA = (cusparseSpMatDescr_t *)f_cusparse_[0].precon_matA;
    cusparseDnVecDescr_t *vecx = (cusparseDnVecDescr_t *)vecx_[0].descr_add;
    cusparseDnVecDescr_t *vecy = (cusparseDnVecDescr_t *)vecy_[0].descr;   

    const cudaDataType_t cudaDataType = CUDA_R_64F;
    const double doubleone  = 1.0;
    const double doublezero = 0.0;

    cusparseSpMV(
		 //		 cusparseHandle[0],
			     PARCEL_CUSPARSE_SPMV_HANDLE,
		 CUSPARSE_OPERATION_NON_TRANSPOSE,
		 &doubleone,
		 matA[0], // desc
		 vecx[0], // desc
		 &doublezero,
		 vecy[0], // desc
		 cudaDataType,
		 CUSPARSE_MV_ALG_DEFAULT,
		 f_cusparse_[0].buffer);
    //		 buffer_);
       
    return ;
  }

  // preconditioner application: z = U^-1 L^-1 r
  void cusparse_precon(
		       f_cusparse *f_cusparse_,
		       devicemem *vecr,PARCEL_INT offset_r,
		       devicemem *vecz,PARCEL_INT offset_z,
		       devicemem *d_precon_val_,
		       devicemem *d_precon_col_,
		       devicemem *d_precon_row_,
		       // devicemem *d_precon_bsrVal_,
		       // devicemem *d_precon_bsrColInd_,
		       // devicemem *d_precon_bsrRowPtr_,
		       PARCEL_INT precon_N,
		       PARCEL_INT precon_nnz
		       ){
    const double doubleone  = 1.0;
    const double doublezero = 0.0;
    const int ipreflag = f_cusparse_[0].ipreflag;

    cusparseHandle_t *cusparseHandle = (cusparseHandle_t *) f_cusparse_[0].cusparseHandle;
    void *buffer = f_cusparse_[0].buffer;
    double *d_work = (double *)f_cusparse_[0].precon_buffer;
    double *d_r = (double *)vecr[0].a_d;
    double *d_z = (double *)vecz[0].a_d;

    if( (ipreflag == KSP_PRECONDITION_BLOCKJACOBI)  || 
	(ipreflag == KSP_PRECONDITION_ADDITIVESCHWARZ) ){

    double *d_precon_val     = (double *)d_precon_val_[0].a_d;
    // PARCEL_INT *d_precon_col = (PARCEL_INT *)d_precon_col_[0].a_d;
    // PARCEL_INT *d_precon_row = (PARCEL_INT *)d_precon_row_[0].a_d;
    int *d_precon_col = (int *)d_precon_col_[0].a_d;
    int *d_precon_row = (int *)d_precon_row_[0].a_d;
    
    cusparseMatDescr_t *descrL = (cusparseMatDescr_t *) f_cusparse_[0].descrL;
    csrsv2Info_t *infoL = (csrsv2Info_t *) f_cusparse_[0].infoL; // csrsv2Info_t
    cusparseMatDescr_t *descrU = (cusparseMatDescr_t *) f_cusparse_[0].descrU;
    csrsv2Info_t *infoU = (csrsv2Info_t *) f_cusparse_[0].infoU; // csrsv2Info_t
    
    cusparseDcsrsv2_solve(
			  //			  cusparseHandle[0],
			     PARCEL_CUSPARSE_PRECON_HANDLE,
			  CUSPARSE_OPERATION_NON_TRANSPOSE,
			  precon_N,
			  precon_nnz,
			  &doubleone,
			  descrL[0],
			  d_precon_val,
			  d_precon_row,
			  d_precon_col,
			  infoL[0],
			  &d_r[offset_r],
			  d_work,
			  CUSPARSE_SOLVE_POLICY_USE_LEVEL,
			  buffer);
    cusparseDcsrsv2_solve(
			  //			  cusparseHandle[0],
			  PARCEL_CUSPARSE_PRECON_HANDLE,
			  CUSPARSE_OPERATION_NON_TRANSPOSE,
			  precon_N,
			  precon_nnz,
			  &doubleone,
			  descrU[0],
			  d_precon_val,
			  d_precon_row,
			  d_precon_col,
			  infoU[0],
			  d_work,
			  &d_z[offset_z],
			  CUSPARSE_SOLVE_POLICY_USE_LEVEL,
			  buffer);
  }

    return ;
  }


}
