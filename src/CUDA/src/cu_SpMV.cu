#include <stdio.h>
#include <mpi.h>
#include <cuda.h>
#include <cuda_runtime.h>

#include "head_tail.hpp"
#include "cu_common.h"


template <class Tdata,class Tcalc>
__global__
void cuda_get_crsmatrix_diag_(
			       PARCEL_INT n,
			       Tdata      *crsA,
			       int *crsRow_ptr,
			       int *crsCol,
			       Tdata    *D,
			       Tdata    *Dinv
			       ){

  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  for(PARCEL_INT i=tid;i<n;i=i+step){
    for(PARCEL_INT j=crsRow_ptr[i];j<crsRow_ptr[i+1];j++){
      if(crsCol[j] == i){
	D[i]    = crsA[j];
	if(D[i]==0.0){
	  D[i] = 0.0;
	  Dinv[i] = 1.0;
	}else{
	  Dinv[i] = 1.0/D[i];
	}

      }
    }
  }
  return ;
}

template <class Tdata,class Tcalc>
__global__
void cuda_get_diamatrix_diag_(
			       PARCEL_INT n,
			       PARCEL_INT nnd   ,
			       PARCEL_INT nBlock,
			       PARCEL_INT diag_offset,
			       Tdata *diaA,
			       Tdata    *D,
			       Tdata    *Dinv
			       ){

  const int thread_id   = threadIdx.x;
  const int block_id    = blockIdx.x;
  const int block_step  = gridDim.x;
  const int thread_step = blockDim.x;

  for(PARCEL_INT kk=1+block_id*nBlock;kk<=n;kk=kk+nBlock*block_step){
    PARCEL_INT kkmax = n < kk+nBlock-1 ? n : kk+nBlock-1;
    for(PARCEL_INT kid=kk+thread_id;kid<=kkmax;kid=kid+thread_step){
      PARCEL_INT j = diag_offset;
      // 対角要素
      Tdata diag_val = diaA[kid + (j-1)*n -1 ];
      D[kid-1]    = diag_val;
      Dinv[kid-1] = 1.0;
      if( diag_val != 0.0 ){
	Dinv[kid-1] = 1.0/diag_val;
      }
      //      printf(" %d , %f \n",kid,Dinv[kid]);
    }
  }

  return ;

}


// __constant__ PARCEL_INT offset_constmem[300];

// Tcalcは演算の精度を指定
template <class Tdata,class Tcalc>
__global__
void MatVec_crs_(
		PARCEL_INT n,
		PARCEL_INT istart,
		PARCEL_INT iend,
	        Tdata      *crsA,
	        PARCEL_INT *crsRow_ptr,
		PARCEL_INT *crsCol,
	        PARCEL_INT *iLs,
		Tdata *xL,Tdata *xD,Tdata *xU,
		Tdata *Ax){
		//		Tdata *xL,Tdata *xD,Tdata *xU,Tdata *Ax){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;

  // crsRow_ptr,crsColが1始まりのインデックスで表現されるため配列番号を-1する
  // n要素の更新
  for(PARCEL_INT k=1+tid;k<=n;k=k+step){
    Ax[k-1]   = 0.0;
    Tdata tmp = 0.0;
    //    PARCEL_INT Row=k-1+istart;
    for(PARCEL_INT j=crsRow_ptr[k-1];j<=crsRow_ptr[k+1-1]-1;j++){
      if( crsCol[j-1]<istart){
	PARCEL_INT Col=iLs[j-1];
	tmp = tmp + crsA[j-1]*xL[Col-1];
	continue;
      }
      if( iend<crsCol[j-1] ){
	PARCEL_INT Col=iLs[j-1];
	tmp = tmp + crsA[j-1]*xU[Col-1];
	continue;
      }
      {
	PARCEL_INT Col=iLs[j-1];
	tmp = tmp + crsA[j-1]*xD[Col-1];
      }
    }
    Ax[k-1] = tmp;
  }
  return ;
}

template <class Tdata,class Tcalc>
__global__
void Diag_Transfer_FP64toBF16_crs_(
		PARCEL_INT n,
	        __nv_bfloat16 *crsA_BF16,
	        double        *crsA_FP64,
	        PARCEL_INT *crsRow_ptr,
		PARCEL_INT *crsCol,
	        PARCEL_INT *iLs,
		int mode){
		//		Tdata *xL,Tdata *xD,Tdata *xU,Tdata *Ax){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;

  // crsRow_ptr,crsColが1始まりのインデックスで表現されるため配列番号を-1する
  for(PARCEL_INT k=1+tid;k<=n;k=k+step){
    //    PARCEL_INT Row=k-1+istart;
    for(PARCEL_INT j=crsRow_ptr[k-1];j<=crsRow_ptr[k+1-1]-1;j++){
      PARCEL_INT Col=iLs[j-1];
      if((Col-1)==(k-1)){
	// 対角位置
	float A_fp32 = crsA_FP64[j-1];
	switch(mode){
	case FP64TOBF16_ROUND_NEAR : 
	  crsA_BF16[j-1] = __float2bfloat16_rn(A_fp32);
	  break;
	case FP64TOBF16_ROUND_ZERO : 
	  crsA_BF16[j-1] = __float2bfloat16_rz(A_fp32);
	  break;
	case FP64TOBF16_ROUND_DOWN : 
	  crsA_BF16[j-1] = __float2bfloat16_rd(A_fp32);
	  break;
	case FP64TOBF16_ROUND_UP : 
	  crsA_BF16[j-1] = __float2bfloat16_ru(A_fp32);
	  break;
	case FP64TOBF16_ROUND_ABSINF : 
	  if(A_fp32>=0.0){
	    crsA_BF16[j-1] = __float2bfloat16_ru(A_fp32);
	  }else{
	    crsA_BF16[j-1] = __float2bfloat16_rd(A_fp32);
	  }
	  break;
	default : 
	  ;
	  break;
	}
      }
    }
  }
  return ;
}

template <class Tdata,class Tcalc>
__global__
void Diag_Transfer_FP64toBF16_localcrs_(
		PARCEL_INT n,
	        double     *crsA_FP64,
	        PARCEL_INT *crsRow_ptr,
		PARCEL_INT *crsCol,
		int mode,
		int mode_diag){
		//		Tdata *xL,Tdata *xD,Tdata *xU,Tdata *Ax){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;

  // crsRow_ptr,crsColが1始まりのインデックスで表現されるため配列番号を-1する
  for(PARCEL_INT k=1+tid;k<=n;k=k+step){
    //    PARCEL_INT Row=k-1+istart;
    for(PARCEL_INT j=crsRow_ptr[k-1];j<=crsRow_ptr[k+1-1]-1;j++){
      //      PARCEL_INT Col=iLs[j-1];
      PARCEL_INT Col=crsCol[j];
      float A_fp32 = crsA_FP64[j];
      __nv_bfloat16 A_BF16;

      if((Col)==(k-1)){
	// 対角位置
	switch(mode_diag){
	case FP64TOBF16_ROUND_NEAR : 
	  A_BF16 = __float2bfloat16_rn(A_fp32);
	  crsA_FP64[j] = (double) __bfloat162float(A_BF16);
	  break;
	case FP64TOBF16_ROUND_ZERO : 
	  A_BF16 = __float2bfloat16_rz(A_fp32);
	  crsA_FP64[j] = (double) __bfloat162float(A_BF16);
	  break;
	case FP64TOBF16_ROUND_DOWN : 
	  A_BF16 = __float2bfloat16_rd(A_fp32);
	  crsA_FP64[j] = (double) __bfloat162float(A_BF16);
	  break;
	case FP64TOBF16_ROUND_UP : 
	  A_BF16 = __float2bfloat16_ru(A_fp32);
	  crsA_FP64[j] = (double) __bfloat162float(A_BF16);
	  break;
	default : 
	  ;
	  break;
	}

      }else{
	  switch(mode){
	  case FP64TOBF16_ROUND_NEAR : 
	    A_BF16 = __float2bfloat16_rn(A_fp32);
	    crsA_FP64[j] = (double) __bfloat162float(A_BF16);
	    break;
	  case FP64TOBF16_ROUND_ZERO : 
	    A_BF16 = __float2bfloat16_rz(A_fp32);
	    crsA_FP64[j] = (double) __bfloat162float(A_BF16);
	    break;
	  case FP64TOBF16_ROUND_DOWN : 
	    A_BF16 = __float2bfloat16_rd(A_fp32);
	    crsA_FP64[j] = (double) __bfloat162float(A_BF16);
	    break;
	  case FP64TOBF16_ROUND_UP : 
	    A_BF16 = __float2bfloat16_ru(A_fp32);
	    crsA_FP64[j] = (double) __bfloat162float(A_BF16);
	    break;
	  default : 
	    ;
	    break;	    
	  }
	}
      //      printf(" %e --> %e \n",A_fp32,crsA_FP64[j]);
      }
    }

  return ;
}

// Tcalcは演算の精度を指定
template <class Tdata,class Tcalc>
__global__
void MatVec_neumann_crs_(
		PARCEL_INT n,
		PARCEL_INT istart,
		PARCEL_INT iend,
	        Tdata      *crsA,
	        PARCEL_INT *crsRow_ptr,
		PARCEL_INT *crsCol,
	        PARCEL_INT *iLs,
		Tdata *xL,Tdata *xD,Tdata *xU,
		Tdata *Ax){
		//		Tdata *xL,Tdata *xD,Tdata *xU,Tdata *Ax){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;

  // crsRow_ptr,crsColが1始まりのインデックスで表現されるため配列番号を-1する
  // n要素の更新
  for(PARCEL_INT k=1+tid;k<=n;k=k+step){
    //    Ax[k-1]   = 0.0;
    Tdata tmp = 0.0;
    //    PARCEL_INT Row=k-1+istart;
    for(PARCEL_INT j=crsRow_ptr[k-1];j<=crsRow_ptr[k+1-1]-1;j++){

#if 1
      if( crsCol[j-1]<istart){
	PARCEL_INT Col=iLs[j-1];
	tmp = tmp + crsA[j-1]*xL[Col-1];
	//	if(tid==0){
	//	printf(" L : j = %08d, %08d \n",k,Col);
	//	}
	continue;
      }
      if( iend<crsCol[j-1] ){
	PARCEL_INT Col=iLs[j-1];
	tmp = tmp + crsA[j-1]*xU[Col-1];
	//	if(tid==0){
	//	printf(" U : j = %08d, %08d \n",k,Col);
	//	}
	continue;
      }
#endif
      {
	PARCEL_INT Col=iLs[j-1];
	tmp = tmp + crsA[j-1]*xD[Col-1];
	// if(tid==0){
	//   printf(" D : j = %08d, %08d \n",k,Col);
	// }
      }
    }
    Ax[k-1] = tmp + xD[k-1];
  }
  return ;
}

template <class Tdata,class Tcalc>
__global__
void JACOBI_MatBF16_crs_(
		PARCEL_INT n,
		PARCEL_INT istart,
		PARCEL_INT iend,
		// Tdata      *crsA,
		__nv_bfloat16 *crsA,
	        PARCEL_INT *crsRow_ptr,
		PARCEL_INT *crsCol,
	        PARCEL_INT *iLs,
		Tdata *b,
		Tdata *xL,Tdata *xD,Tdata *xU,
		Tdata *y){
		//		Tdata *xL,Tdata *xD,Tdata *xU,Tdata *Ax){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;

  // crsRow_ptr,crsColが1始まりのインデックスで表現されるため配列番号を-1する
  // n要素の更新
  for(PARCEL_INT k=1+tid;k<=n;k=k+step){
    //    y[k-1]   = 0.0;
    Tdata tmp  = 0.0;
    Tdata diag = 1.0;
    //    PARCEL_INT Row=k-1+istart;
    for(PARCEL_INT j=crsRow_ptr[k-1];j<=crsRow_ptr[k+1-1]-1;j++){

      if( crsCol[j-1]<istart){
	PARCEL_INT Col=iLs[j-1];
	tmp = tmp + ( (double) __bfloat162float(crsA[j-1]) )*xL[Col-1];
	continue;
      }
      if( iend<crsCol[j-1] ){
	PARCEL_INT Col=iLs[j-1];
	tmp = tmp + ( (double) __bfloat162float(crsA[j-1]) )*xU[Col-1];
	continue;
      }

      {
	PARCEL_INT Col=iLs[j-1];
	if((Col-1)==(k-1)){
	  // 対角位置を記憶
	  diag = ( (double) __bfloat162float(crsA[j-1]) );
	}else{
	  tmp = tmp + ( (double) __bfloat162float(crsA[j-1]) )*xD[Col-1];
	}
      }

    }
    y[k-1] = (b[k-1] - tmp) / diag;
    // y[k-1] = (b[k-1] - tmp);
  }

  return ;

}

template <class Tdata,class Tcalc>
__global__
void JACOBI_MatFP64_crs_(
		PARCEL_INT n,
		PARCEL_INT istart,
		PARCEL_INT iend,
		// Tdata      *crsA,
		double *crsA,
	        PARCEL_INT *crsRow_ptr,
		PARCEL_INT *crsCol,
	        PARCEL_INT *iLs,
		Tdata *b,
		Tdata *xL,Tdata *xD,Tdata *xU,
		Tdata *y){
		//		Tdata *xL,Tdata *xD,Tdata *xU,Tdata *Ax){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;

  // crsRow_ptr,crsColが1始まりのインデックスで表現されるため配列番号を-1する
  // n要素の更新
  for(PARCEL_INT k=1+tid;k<=n;k=k+step){
    //    y[k-1]   = 0.0;
    Tdata tmp  = 0.0;
    Tdata diag = 1.0;
    //    PARCEL_INT Row=k-1+istart;
    for(PARCEL_INT j=crsRow_ptr[k-1];j<=crsRow_ptr[k+1-1]-1;j++){

      if( crsCol[j-1]<istart){
	PARCEL_INT Col=iLs[j-1];
	tmp = tmp + crsA[j-1] * xL[Col-1];
	//	printf("up 1 \n");
	continue;
      }
      if( iend<crsCol[j-1] ){
	PARCEL_INT Col=iLs[j-1];
	tmp = tmp + crsA[j-1] * xU[Col-1];
	//	printf("down 2 \n");
	continue;
      }

      {
	PARCEL_INT Col=iLs[j-1];
	if((Col-1)==(k-1)){
	//if( crsCol[j-1]==(k+istart-1)){
	  // 対角位置を記憶
	  diag = crsA[j-1];
	}else{
	  tmp = tmp + crsA[j-1] * xD[Col-1];
	}
      }
    }
    y[k-1] = (b[k-1] - tmp) / diag;
    //    y[k-1] = (b[k-1] - tmp) ;
    //    y[k-1] = b[k-1] / diag;
  }

  return ;

}


template <class Tdata,class Tcalc>
__global__
void Matvec_dia_(
		 PARCEL_INT n     ,
		 PARCEL_INT nnd   ,
		 PARCEL_INT istart,
		 PARCEL_INT nBlock,
		 Tdata *x,
		 Tdata *Ax,
		 Tdata *diaA,
		 PARCEL_INT *offset_,
		 PARCEL_INT *iLs,
		 Tdata *rbuf
		 ){

  //  PARCEL_INT nBlock = nBlock_; 

//  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
//  const int step = blockDim.x*gridDim.x;
  const int thread_id = threadIdx.x;
  const int block_id  = blockIdx.x;
  const int block_step  = gridDim.x;
  const int thread_step = blockDim.x;

#if 0
   __shared__ PARCEL_INT offset_l[30];
   if(thread_id==0){
     for(PARCEL_INT j=1;j<=nnd;j++){
       offset_l[j-1]  = offset_[j-1]; 
     }
   }
   __syncthreads();
#endif

  for(PARCEL_INT kk=1+block_id*nBlock;kk<=n;kk=kk+nBlock*block_step){
    PARCEL_INT kkmax = n < kk+nBlock-1 ? n : kk+nBlock-1;
    for(PARCEL_INT kid=kk+thread_id;kid<=kkmax;kid=kid+thread_step){
      //      Ax[kid-1]=0.0;
      Tdata y_ = 0.0;
      for(PARCEL_INT j=1;j<=nnd;j++){
	const PARCEL_INT offset = offset_[j-1];
	//	const PARCEL_INT offset = offset_constmem[j-1];

	const PARCEL_INT minCol = kk    + istart - 1 + offset;
	const PARCEL_INT maxCol = kkmax + istart - 1 + offset;
	const PARCEL_INT minRow = 1 + istart - 1;
	const PARCEL_INT maxRow = n + istart - 1;

	PARCEL_INT kmin = minCol > minRow ? minCol : minRow;
	PARCEL_INT kmax = maxCol < maxRow ? maxCol : maxRow;
	kmin = kmin-istart+1;
	kmax = kmax-istart+1;

	// 領域内
	if( ((kmin-offset)<=kid) && (kid<=(kmax-offset)) ){	  
	  y_ = y_ + diaA[kid + (j-1)*n -1] * x[kid + offset -1];	
	  continue;
	}

	// 領域外(if文は不要)
	if( (kk<=kid) && (kid<=kkmax) ){ 
	  PARCEL_INT i = iLs[ kid + (j-1)*n  -1];
	  y_ = y_ + diaA[kid + (j-1)*n -1 ] * rbuf[i-1];
	}
      }
      Ax[kid-1]=y_;
    }

  }
  
  return ; 
  
}

// SS(:) = eta*AS(:) - zeta*r(:)
// ----------
//  z(:) = alpha*A*x(:) + beta*y(:)
// Ax(:) = A*x(:)
template <class Tdata,class Tcalc>
__global__
void Matvec_cbcg1_dia_(
		 PARCEL_INT n     ,
		 PARCEL_INT nnd   ,
		 PARCEL_INT istart,
		 PARCEL_INT nBlock,
		 Tdata *x,
		 Tdata *Ax,
		 Tdata *y,
		 Tdata *z,
		 Tdata alpha,
		 Tdata beta ,
		 Tdata *diaA,
		 PARCEL_INT *offset_,
		 PARCEL_INT *iLs,
		 Tdata *rbuf
		 ){

  const int thread_id = threadIdx.x;
  const int block_id  = blockIdx.x;
  const int block_step  = gridDim.x;
  const int thread_step = blockDim.x;

  for(PARCEL_INT kk=1+block_id*nBlock;kk<=n;kk=kk+nBlock*block_step){
    PARCEL_INT kkmax = n < kk+nBlock-1 ? n : kk+nBlock-1;
    for(PARCEL_INT kid=kk+thread_id;kid<=kkmax;kid=kid+thread_step){
      Tdata Ax_ = 0.0;
      for(PARCEL_INT j=1;j<=nnd;j++){
	const PARCEL_INT offset = offset_[j-1];

	const PARCEL_INT minCol = kk    + istart - 1 + offset;
	const PARCEL_INT maxCol = kkmax + istart - 1 + offset;
	const PARCEL_INT minRow = 1 + istart - 1;
	const PARCEL_INT maxRow = n + istart - 1;

	PARCEL_INT kmin = minCol > minRow ? minCol : minRow;
	PARCEL_INT kmax = maxCol < maxRow ? maxCol : maxRow;
	kmin = kmin-istart+1;
	kmax = kmax-istart+1;

	// 領域内
	if( ((kmin-offset)<=kid) && (kid<=(kmax-offset)) ){	  
	  Ax_ = Ax_ + diaA[kid + (j-1)*n -1] * x[kid + offset -1];	
	  continue;
	}

	// 領域外(if文は不要)
	if( (kk<=kid) && (kid<=kkmax) ){ 
	  PARCEL_INT i = iLs[ kid + (j-1)*n  -1];
	  Ax_ = Ax_ + diaA[kid + (j-1)*n -1 ] * rbuf[i-1];
	}
      }
      z[kid-1]   = alpha*Ax_ + beta*y[kid-1];
      Ax[kid-1]  = Ax_;
    }

  }
  
  return ; 
  
}

// SS(i) = 2.0*eta*AS-2.0*zeta*SS(i-1)-SS(i-2)
// ----------
//  w(:) = alpha*A*x(:) + beta*y(:) + gamma*z(:)
// Ax(:) = A*x(:)
template <class Tdata,class Tcalc>
__global__
void Matvec_cbcg2_dia_(
		 PARCEL_INT n     ,
		 PARCEL_INT nnd   ,
		 PARCEL_INT istart,
		 PARCEL_INT nBlock,
		 Tdata *x,
		 Tdata *Ax,
		 Tdata *y,
		 Tdata *z,
		 Tdata *w,
		 Tdata alpha,
		 Tdata beta ,
		 Tdata gamma ,
		 Tdata *diaA,
		 PARCEL_INT *offset_,
		 PARCEL_INT *iLs,
		 Tdata *rbuf
		 ){

  const int thread_id = threadIdx.x;
  const int block_id  = blockIdx.x;
  const int block_step  = gridDim.x;
  const int thread_step = blockDim.x;

  for(PARCEL_INT kk=1+block_id*nBlock;kk<=n;kk=kk+nBlock*block_step){
    PARCEL_INT kkmax = n < kk+nBlock-1 ? n : kk+nBlock-1;
#if 0
    //    for(PARCEL_INT kid=kk+thread_id;kid<=kkmax;kid=kid+thread_step){
    PARCEL_INT kid=kk+thread_id;
      Tdata Ax_ = 0.0;
#pragma unroll 4
      for(PARCEL_INT j=1;j<=nnd;j++){
	const PARCEL_INT offset = offset_[j-1];

	const PARCEL_INT minCol = kk    + istart - 1 + offset;
	const PARCEL_INT maxCol = kkmax + istart - 1 + offset;
	const PARCEL_INT minRow = 1 + istart - 1;
	const PARCEL_INT maxRow = n + istart - 1;

	PARCEL_INT kmin = minCol > minRow ? minCol : minRow;
	PARCEL_INT kmax = maxCol < maxRow ? maxCol : maxRow;
	kmin = kmin-istart+1;
	kmax = kmax-istart+1;

	// 領域内
	if( ((kmin-offset)<=kid) && (kid<=(kmax-offset)) ){	  
	  Ax_ = Ax_ + diaA[kid + (j-1)*n -1] * x[kid + offset -1];	
	  continue;
	}

	// 領域外(if文は不要)
	if( (kk<=kid) && (kid<=kkmax) ){ 
	  PARCEL_INT i = iLs[ kid + (j-1)*n  -1];
	  Ax_ = Ax_ + diaA[kid + (j-1)*n -1 ] * rbuf[i-1];
	}
      }
      //      w[kid-1] = alpha*Ax_ + beta*y[kid-1] + gamma*z[kid-1];
      Ax[kid-1] = Ax_;
      //    }
      //    for(PARCEL_INT kid=kk+thread_id;kid<=kkmax;kid=kid+thread_step){
      w[kid-1] = alpha*Ax_ + beta*y[kid-1] + gamma*z[kid-1];
      //    }

#else
    for(PARCEL_INT kid=kk+thread_id;kid<=kkmax;kid=kid+thread_step){
      Tdata Ax_ = 0.0;
      for(PARCEL_INT j=1;j<=nnd;j++){
	const PARCEL_INT offset = offset_[j-1];

	const PARCEL_INT minCol = kk    + istart - 1 + offset;
	const PARCEL_INT maxCol = kkmax + istart - 1 + offset;
	const PARCEL_INT minRow = 1 + istart - 1;
	const PARCEL_INT maxRow = n + istart - 1;

	PARCEL_INT kmin = minCol > minRow ? minCol : minRow;
	PARCEL_INT kmax = maxCol < maxRow ? maxCol : maxRow;
	kmin = kmin-istart+1;
	kmax = kmax-istart+1;

	// 領域内
	if( ((kmin-offset)<=kid) && (kid<=(kmax-offset)) ){	  
	  Ax_ = Ax_ + diaA[kid + (j-1)*n -1] * x[kid + offset -1];	
	  continue;
	}

	// 領域外(if文は不要)
	if( (kk<=kid) && (kid<=kkmax) ){ 
	  PARCEL_INT i = iLs[ kid + (j-1)*n  -1];
	  Ax_ = Ax_ + diaA[kid + (j-1)*n -1 ] * rbuf[i-1];
	}
      }
      //      w[kid-1] = alpha*Ax_ + beta*y[kid-1] + gamma*z[kid-1];
      Ax[kid-1] = Ax_;
    }
    for(PARCEL_INT kid=kk+thread_id;kid<=kkmax;kid=kid+thread_step){
      w[kid-1] = alpha*Ax[kid-1] + beta*y[kid-1] + gamma*z[kid-1];
    }
#endif
  }
  
  return ; 
  
}



template <class Tdata,class Tcalc>
__global__
void Matvec_neumann_dia_(
		 PARCEL_INT n     ,
		 PARCEL_INT nnd   ,
		 PARCEL_INT istart,
		 PARCEL_INT nBlock,
		 Tdata *x,
		 Tdata *Ax,
		 Tdata *diaA,
		 PARCEL_INT *offset_,
		 PARCEL_INT *iLs,
		 Tdata *rbuf
		 ){

  //  PARCEL_INT nBlock = nBlock_; 

//  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
//  const int step = blockDim.x*gridDim.x;
  const int thread_id = threadIdx.x;
  const int block_id  = blockIdx.x;
  const int block_step  = gridDim.x;
  const int thread_step = blockDim.x;

#if 0
   __shared__ PARCEL_INT offset_l[30];
   if(thread_id==0){
     for(PARCEL_INT j=1;j<=nnd;j++){
       offset_l[j-1]  = offset_[j-1]; 
     }
   }
   __syncthreads();
#endif

  for(PARCEL_INT kk=1+block_id*nBlock;kk<=n;kk=kk+nBlock*block_step){
    PARCEL_INT kkmax = n < kk+nBlock-1 ? n : kk+nBlock-1;
    for(PARCEL_INT kid=kk+thread_id;kid<=kkmax;kid=kid+thread_step){
      //      Ax[kid-1]=0.0;
      Tdata y_ = 0.0;
      for(PARCEL_INT j=1;j<=nnd;j++){
	// コンスタントメモリを利用しても対して効果ない?
	const PARCEL_INT offset = offset_[j-1];
	if(offset==0)continue;
	//	const PARCEL_INT offset = offset_constmem[j-1];

	const PARCEL_INT minCol = kk    + istart - 1 + offset;
	const PARCEL_INT maxCol = kkmax + istart - 1 + offset;
	const PARCEL_INT minRow = 1 + istart - 1;
	const PARCEL_INT maxRow = n + istart - 1;

	PARCEL_INT kmin = minCol > minRow ? minCol : minRow;
	PARCEL_INT kmax = maxCol < maxRow ? maxCol : maxRow;
	kmin = kmin-istart+1;
	kmax = kmax-istart+1;

	// 領域内
	if( ((kmin-offset)<=kid) && (kid<=(kmax-offset)) ){	  
	  y_ = y_ + diaA[kid + (j-1)*n -1] * x[kid + offset -1];	
	  continue;
	}

	// 領域外(if文は不要)
	if( (kk<=kid) && (kid<=kkmax) ){ 
	  PARCEL_INT i = iLs[ kid + (j-1)*n  -1];
	  y_ = y_ + diaA[kid + (j-1)*n -1 ] * rbuf[i-1];
	}
      }
      Ax[kid-1] = y_ + x[kid-1];
    }

  }
  
  return ; 
  
}

extern "C" {

  void  cuda_get_crsmatrix_diag(
		 PARCEL_INT n,
		 devicemem *crsA_, 
		 devicemem *crsRow_ptr_, 
		 devicemem *crsCol_, 
		 devicemem *D_   , PARCEL_INT offset_D,
		 devicemem *Dinv_, PARCEL_INT offset_Dinv
				){
    double* crsA           = (double*)crsA_[0].a_d;
    int* crsRow_ptr = (int*)crsRow_ptr_[0].a_d;
    int* crsCol     = (int*)crsCol_[0].a_d;

    double* D    = (double*)D_[0].a_d;
    double* Dinv = (double*)Dinv_[0].a_d;

    cuda_get_crsmatrix_diag_<double,double> <<<512,512>>>(
							  n,
							  crsA,
							  crsRow_ptr,
							  crsCol,
							  &D[offset_D],
							  &Dinv[offset_Dinv]
							  );

    return ;

  }

  void  cuda_get_diamatrix_diag(
		 PARCEL_INT n,
		 PARCEL_INT nnd,
		 PARCEL_INT diag_offset,
		 devicemem *diaA_, 
		 devicemem *D_   , PARCEL_INT offset_D,
		 devicemem *Dinv_, PARCEL_INT offset_Dinv
				){

    double* diaA = (double*)diaA_[0].a_d; 
    double* D    = (double*)D_[0].a_d;
    double* Dinv = (double*)Dinv_[0].a_d;

    cuda_get_diamatrix_diag_<double,double> <<<4096,512>>>(
							  n,
							  nnd,
							  512,
							  diag_offset,
							  diaA,
							  &D[offset_D],
							  &Dinv[offset_Dinv]
							  );

    return ;

  }

  void cuda_MatVec_crs(
		       int numGPUblock,
		       int numGPUthread,
		       PARCEL_INT n,
		       PARCEL_INT istart,
		       PARCEL_INT iend,
		       devicemem *crsA, // FP
		       devicemem *crsRow_ptr, // PARCEL_INT
		       devicemem *crsCol, // PARCEL_INT
		       devicemem *iLs, // PARCEL_INT
		       devicemem *xL_, // FP
		       PARCEL_INT offset_xL,
		       devicemem *xD_, // FP
		       PARCEL_INT offset_xD,
		       devicemem *xU_, // FP
		       PARCEL_INT offset_xU,
		       devicemem *Ax_,  // FP
		       PARCEL_INT offset_Ax
		   ){

    double* xL = (double*)xL_[0].a_d;
    double* xD = (double*)xD_[0].a_d;
    double* xU = (double*)xU_[0].a_d;
    double* Ax = (double*)Ax_[0].a_d;
    //    printf(" istart = %08d, iend = %08d \n",istart,iend);
    MatVec_crs_<double,double> 
      //       <<<numGPUblock,numGPUthread>>>(
      <<<4096,512>>>(
							     n,
							     istart,
							     iend,
					     (double*)crsA[0].a_d,
					     (PARCEL_INT *)crsRow_ptr[0].a_d,
					     (PARCEL_INT *)crsCol[0].a_d,
					     (PARCEL_INT *)iLs[0].a_d,
					     &xL[offset_xL],
					     &xD[offset_xD],
					     &xU[offset_xU],
					     &Ax[offset_Ax]);
    
  }

  void cuda_MatVec_neumann_crs(
		       int numGPUblock,
		       int numGPUthread,
		       PARCEL_INT n,
		       PARCEL_INT istart,
		       PARCEL_INT iend,
		       devicemem *crsA, // FP
		       devicemem *crsRow_ptr, // PARCEL_INT
		       devicemem *crsCol, // PARCEL_INT
		       devicemem *iLs, // PARCEL_INT
		       devicemem *xL_, // FP
		       PARCEL_INT offset_xL,
		       devicemem *xD_, // FP
		       PARCEL_INT offset_xD,
		       devicemem *xU_, // FP
		       PARCEL_INT offset_xU,
		       devicemem *Ax_,  // FP
		       PARCEL_INT offset_Ax
		   ){
    //    printf(" istart = %08d, iend = %08d \n",istart,iend);
    double* xL = (double*)xL_[0].a_d;
    double* xD = (double*)xD_[0].a_d;
    double* xU = (double*)xU_[0].a_d;
    double* Ax = (double*)Ax_[0].a_d;
    //    printf(" istart = %08d, iend = %08d \n",istart,iend);
    MatVec_neumann_crs_<double,double>  
      //      <<<4096,1024>>>(
      //ok
      <<<64,512>>>(
      //      <<<512,128>>>(
      //      <<<numGPUblock,numGPUthread>>>(
		      n,
		      istart,
		      iend,
		      (double*)crsA[0].a_d,
		      (PARCEL_INT *)crsRow_ptr[0].a_d,
		      (PARCEL_INT *)crsCol[0].a_d,
		      (PARCEL_INT *)iLs[0].a_d,
		      &xL[offset_xL],
		      &xD[offset_xD],
		      &xU[offset_xU],
		      &Ax[offset_Ax]);
    
  }


  // コンスタントメモリにオフセット量を書き込む
  void cuda_MatVec_dia_init(
			    PARCEL_INT nnd   ,
			    PARCEL_INT *offset
			    ){
#if 0
    cudaMemcpyToSymbol(offset_constmem, offset, sizeof(PARCEL_INT)*nnd);
#endif
  }
  
  void cuda_MatVec_dia(
		       int numGPUblock,
		       int numGPUthread,
		       PARCEL_INT n     ,
		       PARCEL_INT nnd   ,
		       PARCEL_INT istart,
		       PARCEL_INT nBlock,
		       devicemem *x_ ,PARCEL_INT offset_x,
		       devicemem *Ax_,PARCEL_INT offset_Ax,
		       devicemem *diaA,
		       devicemem *offset_,
		       devicemem *iLs,
		       devicemem *rbuf
		       ){
    double *x  = (double *)x_[0].a_d;
    double *Ax = (double *)Ax_[0].a_d;
    Matvec_dia_<double,double> <<<4096,512>>>(
    //    Matvec_dia_<double,double> <<<4096,512>>>(
    //    Matvec_dia_<double,double> <<<512,512>>>(
		 n     ,
		 nnd   ,
		 istart,
		 512,
		 //		 nBlock,
		 &x[offset_x],
		 &Ax[offset_Ax],
		 (double *)diaA[0].a_d,
		 (PARCEL_INT *)offset_[0].a_d,
		 (PARCEL_INT *)iLs[0].a_d,
		 (double *)rbuf[0].a_d
					 );
    

    return ;
  }

  void cuda_MatVec_cbcg1_dia(
		       int numGPUblock_,
		       int numGPUthread_,
		       PARCEL_INT n     ,
		       PARCEL_INT nnd   ,
		       PARCEL_INT istart,
		       PARCEL_INT nBlock,
		       devicemem *x_ ,PARCEL_INT offset_x,
		       devicemem *Ax_,PARCEL_INT offset_Ax,
		       devicemem *y_ ,PARCEL_INT offset_y,
		       devicemem *z_ ,PARCEL_INT offset_z,
		       devicemem *diaA,
		       double alpha,
		       double beta ,
		       devicemem *offset_,
		       devicemem *iLs,
		       devicemem *rbuf
		       ){
    double *x  = (double *)x_[0].a_d;
    double *Ax = (double *)Ax_[0].a_d;
    double *y = (double *)y_[0].a_d;
    double *z = (double *)z_[0].a_d;

    const int numGPUthread = 512;
    Matvec_cbcg1_dia_<double,double> <<<4096,numGPUthread>>>(
    //    Matvec_dia_<double,double> <<<4096,512>>>(
    //    Matvec_dia_<double,double> <<<512,512>>>(
		 n     ,
		 nnd   ,
		 istart,
		 numGPUthread,
		 //		 nBlock,
		 &x[offset_x],
		 &Ax[offset_Ax],
		 &y[offset_y],
		 &z[offset_z],
		 alpha,
		 beta ,
		 (double *)diaA[0].a_d,
		 (PARCEL_INT *)offset_[0].a_d,
		 (PARCEL_INT *)iLs[0].a_d,
		 (double *)rbuf[0].a_d
					 );
    

    return ;
  }

  void cuda_MatVec_cbcg2_dia(
		       int numGPUblock_,
		       int numGPUthread_,
		       PARCEL_INT n     ,
		       PARCEL_INT nnd   ,
		       PARCEL_INT istart,
		       PARCEL_INT nBlock,
		       devicemem *x_ ,PARCEL_INT offset_x,
		       devicemem *Ax_,PARCEL_INT offset_Ax,
		       devicemem *y_ ,PARCEL_INT offset_y,
		       devicemem *z_ ,PARCEL_INT offset_z,
		       devicemem *w_ ,PARCEL_INT offset_w,
		       devicemem *diaA,
		       double alpha,
		       double beta ,
		       double gamma ,
		       devicemem *offset_,
		       devicemem *iLs,
		       devicemem *rbuf
		       ){
    double *x  = (double *)x_[0].a_d;
    double *Ax = (double *)Ax_[0].a_d;
    double *y = (double *)y_[0].a_d;
    double *z = (double *)z_[0].a_d;
    double *w = (double *)w_[0].a_d;

    //    Matvec_cbcg2_dia_<double,double> <<<4096,1024>>>(
    //    Matvec_cbcg2_dia_<double,double> <<<4096,1024>>>(
    const int numGPUthread = 512;
    //    Matvec_cbcg2_dia_<double,double> <<<4096,numGPUthread>>>(  
    Matvec_cbcg2_dia_<double,double> <<<2048,numGPUthread>>>(  
		 n     ,
		 nnd   ,
		 istart,
		 numGPUthread,
		 //		 nBlock,
		 &x[offset_x],
		 &Ax[offset_Ax],
		 &y[offset_y],
		 &z[offset_z],
		 &w[offset_w],
		 alpha,
		 beta ,
		 gamma ,
		 (double *)diaA[0].a_d,
		 (PARCEL_INT *)offset_[0].a_d,
		 (PARCEL_INT *)iLs[0].a_d,
		 (double *)rbuf[0].a_d
					 );
    

    return ;
  }


  void cuda_Diag_Transfer_FP64toBF16_crs(
		       int numGPUblock,
		       int numGPUthread,
		       PARCEL_INT n,
		       devicemem *crsA_BF16, // FP
		       devicemem *crsA_FP64, // FP
		       devicemem *crsRow_ptr, // PARCEL_INT
		       devicemem *crsCol, // PARCEL_INT
		       devicemem *iLs, // PARCEL_INT
		       int mode
		   ){
    //    printf(" istart = %08d, iend = %08d \n",istart,iend);
    Diag_Transfer_FP64toBF16_crs_<double,double> 
       <<<numGPUblock,numGPUthread>>>(
      //      <<<4096,1024>>>(
				      n,
				      (__nv_bfloat16*)crsA_BF16[0].a_d,
				      (double*)crsA_FP64[0].a_d,
				      (PARCEL_INT *)crsRow_ptr[0].a_d,
				      (PARCEL_INT *)crsCol[0].a_d,
				      (PARCEL_INT *)iLs[0].a_d,
				      mode);
    
  }

  void cuda_Diag_Transfer_FP64toBF16_localcrs(
		       int numGPUblock,
		       int numGPUthread,
		       PARCEL_INT n,
		       devicemem *crsA_FP64, // FP
		       devicemem *crsRow_ptr, // PARCEL_INT
		       devicemem *crsCol, // PARCEL_INT
		       int mode,
		       int mode_diag
		   ){
    //    printf(" istart = %08d, iend = %08d \n",istart,iend);
    Diag_Transfer_FP64toBF16_localcrs_<double,double> 
       <<<numGPUblock,numGPUthread>>>(
      //      <<<4096,1024>>>(
				      n,
				      (double*)crsA_FP64[0].a_d,
				      (PARCEL_INT *)crsRow_ptr[0].a_d,
				      (PARCEL_INT *)crsCol[0].a_d,
				      mode,
				      mode_diag);
    
  }
  

  void cuda_MatVec_neumann_dia(
		       int numGPUblock,
		       int numGPUthread,
		       PARCEL_INT n     ,
		       PARCEL_INT nnd   ,
		       PARCEL_INT istart,
		       PARCEL_INT nBlock,
		       devicemem *x_,
		       PARCEL_INT offset_x,
		       devicemem *Ax_,
		       PARCEL_INT offset_Ax,
		       devicemem *diaA,
		       devicemem *offset_,
		       devicemem *iLs,
		       devicemem *rbuf
		       ){
    double *x = (double *)x_[0].a_d;
    double *Ax = (double *)Ax_[0].a_d;
    //    Matvec_neumann_dia_<double,double> <<<4096,1024>>>(
    Matvec_neumann_dia_<double,double> <<<4096,512>>>(
		 n     ,
		 nnd   ,
		 istart,
		 512,
		 //		 nBlock,
		 &x[offset_x],
		 &Ax[offset_Ax],
		 (double *)diaA[0].a_d,
		 (PARCEL_INT *)offset_[0].a_d,
		 (PARCEL_INT *)iLs[0].a_d,
		 (double *)rbuf[0].a_d
					 );

    return ;
  }

  void cuda_JACOBI_MatBF16_crs(
		       int numGPUblock,
		       int numGPUthread,
		       PARCEL_INT n,
		       PARCEL_INT istart,
		       PARCEL_INT iend,
		       devicemem *crsA, // FP
		       devicemem *crsRow_ptr, // PARCEL_INT
		       devicemem *crsCol, // PARCEL_INT
		       devicemem *iLs, // PARCEL_INT
		       devicemem *b_, // FP
		       PARCEL_INT offset_b,
		       devicemem *xL_, // FP
		       PARCEL_INT offset_xL,
		       devicemem *xD_, // FP
		       PARCEL_INT offset_xD,
		       devicemem *xU_, // FP
		       PARCEL_INT offset_xU,
		       devicemem *y_,  // FP
		       PARCEL_INT offset_y
		   ){
    //    printf(" istart = %08d, iend = %08d \n",istart,iend);
    double* b  = (double*)b_[0].a_d;
    double* xL = (double*)xL_[0].a_d;
    double* xD = (double*)xD_[0].a_d;
    double* xU = (double*)xU_[0].a_d;
    double* y  = (double*)y_[0].a_d;

    JACOBI_MatBF16_crs_<double,double> 
       <<<numGPUblock,numGPUthread>>>(
      //      <<<4096,1024>>>(
							     n,
							     istart,
							     iend,
					     (__nv_bfloat16*)crsA[0].a_d,
					     (PARCEL_INT *)crsRow_ptr[0].a_d,
					     (PARCEL_INT *)crsCol[0].a_d,
					     (PARCEL_INT *)iLs[0].a_d,
					     &b[offset_b],
					     &xL[offset_xL],
					     &xD[offset_xD],
					     &xU[offset_xU],
					     &y[offset_y]);
    
  }

  void cuda_JACOBI_MatFP64_crs(
		       int numGPUblock,
		       int numGPUthread,
		       PARCEL_INT n,
		       PARCEL_INT istart,
		       PARCEL_INT iend,
		       devicemem *crsA, // FP
		       devicemem *crsRow_ptr, // PARCEL_INT
		       devicemem *crsCol, // PARCEL_INT
		       devicemem *iLs, // PARCEL_INT
		       devicemem *b_, // FP
		       PARCEL_INT offset_b,
		       devicemem *xL_, // FP
		       PARCEL_INT offset_xL,
		       devicemem *xD_, // FP
		       PARCEL_INT offset_xD,
		       devicemem *xU_, // FP
		       PARCEL_INT offset_xU,
		       devicemem *y_,  // FP
		       PARCEL_INT offset_y
		   ){
    //    printf(" istart = %08d, iend = %08d \n",istart,iend);
    double* b  = (double*)b_[0].a_d;
    double* xL = (double*)xL_[0].a_d;
    double* xD = (double*)xD_[0].a_d;
    double* xU = (double*)xU_[0].a_d;
    double* y  = (double*)y_[0].a_d;
    JACOBI_MatFP64_crs_<double,double> 
      ///       <<<numGPUblock,numGPUthread>>>(
      <<<4096,512>>>(
				      n,
				      istart,
				      iend,
				      (double *)crsA[0].a_d,
				      (PARCEL_INT *)crsRow_ptr[0].a_d,
				      (PARCEL_INT *)crsCol[0].a_d,
				      (PARCEL_INT *)iLs[0].a_d,
					     &b[offset_b],
					     &xL[offset_xL],
					     &xD[offset_xD],
					     &xU[offset_xU],
					     &y[offset_y]);
    
  }


}


