#include <cuda.h>
#include <cuda_runtime.h>

#include "cu_common.h"
#include "reduciton_cuda.hpp"
#include "cu_common.hpp"

// template <class Tdata,class Tcalc>
// __device__ 
// void House_local_Ronly_omp(){
// 
//   return ;
// }

template <class Tdata,class Tcalc>
__device__
void House_local_Ronly_(
			const int warp_id,
			const PARCEL_INT  n,
			const PARCEL_INT  s,
			Tdata *A     , PARCEL_INT lda,
			Tdata *R     , PARCEL_INT ldr,
			Tdata *A_TSQR, PARCEL_INT lda_tsqr,
			const PARCEL_INT  kk, 
			const PARCEL_INT  nblock,
			const PARCEL_INT  iRflag,
			Tdata *shared_mem
			){

  Tdata *tmp_t1 = &shared_mem[0];
  Tdata *tmp_y1 = &shared_mem[s];
  for(PARCEL_INT i = warp_id;i<s;i=i+WARP_SIZE){
    tmp_t1[i] = 0.0;
    tmp_y1[i] = 0.0;
  }
  PARCEL_INT i = 1;
  Tdata tmp2 = 0.0;
  for(PARCEL_INT k_ = kk + warp_id; k_<=kk+nblock-1; k_=k_+WARP_SIZE){
    Tdata  Ak_i=A[k_ + (i-1)*lda -1];
    A_TSQR[k_-kk+1 + (i-1)*lda_tsqr -1 ] = Ak_i;
    tmp2 = tmp2 + Ak_i*Ak_i;
    //    printf("tmp2 add Ak_i: %d %e \n",k_,Ak_i*Ak_i);
  }

  PARCEL_INT nHouse = nblock;
  if(iRflag==1){
    for(PARCEL_INT l_ = 1 + warp_id; l_<=s ;l_ = l_ + WARP_SIZE){
      Tdata Rl_i = R[l_+(i-1)*ldr -1];
      A_TSQR[l_+nblock+(i-1)*lda_tsqr -1] = Rl_i;
      tmp2 = tmp2 + Rl_i*Rl_i;
      //      printf("tmp2 add Rl_i: %e \n",Rl_i*Rl_i);
    }
    nHouse = nHouse + s;
  }

  if(warp_id==0){
    //    tmp2 = tmp2 - A_TSQR[1+(1-1)*lda_tsqr-1]*A_TSQR[1+(1-1)*lda_tsqr-1] ;
    Tdata A11 = A_TSQR[1+(1-1)*lda_tsqr-1];
    tmp2 = tmp2 - A11*A11 ;
    //      printf("tmp2 add A11: %e \n",-A11*A11);
  }

  i=2;

  //  if(warp_id==0){printf(" kk = %d \n",kk);}
  for(PARCEL_INT k_=kk+warp_id;k_<=kk+nblock-1;k_=k_+WARP_SIZE){
    A_TSQR[k_-kk+1+(i-1)*lda_tsqr -1]=A[k_ + (i-1)*lda -1];
    //    printf(" A_TSQR i = 2 : %d %e\n",k_,A_TSQR[k_-kk+1+(i-1)*lda_tsqr -1]);
  }
  if(iRflag==1){
    for(PARCEL_INT l_=1+warp_id; l_<=s;l_=l_+WARP_SIZE){
      A_TSQR[l_+nblock+(i-1)*lda_tsqr -1] = R[l_+(i-1)*ldr -1];
    }
  }

  //  Tdata tmp2_ = warpReduceSum<Tdata,Tcalc>(tmp2);
  for(PARCEL_INT i=1;i<=s;i++){
    Tdata tmp2_ = warpReduceSum<Tdata,Tcalc>(tmp2);
    tmp2 = tmp2_;
    //    if(warp_id==0){printf(" tmp2 = %e \n",tmp2);}

    // 1warpで計算するところは1warpで計算
    if(warp_id==0){
      //    if(warp_id==0){
      //      tmp3 = tmp2 + A_TSQR[i + (i-1)*lda_tsqr -1] * A_TSQR[i + (i-1)* lda_tsqr -1];
      Tdata A_TSQR_ii = A_TSQR[i + (i-1)*lda_tsqr -1];
      Tdata tmp3 = tmp2 + 
	A_TSQR_ii * 
	A_TSQR_ii;
      tmp3 = sqrt(tmp3);
      Tdata y1 = A_TSQR_ii > 0.0 ? A_TSQR_ii+tmp3 : A_TSQR_ii-tmp3;

      Tdata tmp  = tmp2+y1*A_TSQR_ii;
      tmp2 = tmp2+y1*y1;
      //      if(warp_id==0){printf(" tmp3 = %e \n",tmp3);}
      if(tmp2==0.0){
	tmp2 = 1.0e-300;
      }
      //      printf(" tmp2 = %e \n",tmp2);
      Tdata t1=2.0/tmp2;
      A_TSQR[i + (i-1)*lda_tsqr -1] = A_TSQR_ii-t1*y1*tmp;
      tmp_y1[i-1] = y1;
      tmp_t1[i-1] = t1;
      //      if(warp_id==0){printf(" y1 = %e \n",y1);}
    }

    if(i != s) {
      //      tmp = tmp_y1(i)*A_TSQR(i,i+1) --> ここは後でwarp_id==0が足し合わせる
      Tdata tmp = 0.0; 
      for(PARCEL_INT j=i+1+warp_id;j<=nHouse;j=j+WARP_SIZE){
	tmp = tmp 
	  + A_TSQR[j + (i-1)*lda_tsqr -1] * A_TSQR[j + (i+1-1)*lda_tsqr -1] ;
	//		printf(" tmp add : %d %e \n",
	//		       j,A_TSQR[j + (i-1)*lda_tsqr -1] * A_TSQR[j + (i+1-1)*lda_tsqr -1]);

      }
      Tdata tmp_ = warpReduceSum<Tdata,Tcalc>(tmp);
      tmp = tmp_;
      if(warp_id==0){
	tmp = tmp + tmp_y1[i-1]*A_TSQR[i + (i+1-1)*lda_tsqr -1 ];
	tmp = tmp_t1[i-1]*tmp;

	//		printf(" tmp add : %e \n",
	//		       tmp_y1[i-1]*A_TSQR[i + (i+1-1)*lda_tsqr -1 ]);
	//		printf(" tmp mul : %e \n",
	//		       tmp_t1[i-1]);

	A_TSQR[i+(i+1-1)*lda_tsqr -1] = 
	  A_TSQR[i+(i+1-1)*lda_tsqr -1]-tmp_y1[i-1]*tmp;
	A_TSQR[i+1 + (i+1-1) * lda_tsqr -1] = 
	  A_TSQR[i+1 + (i+1-1)*lda_tsqr -1] - A_TSQR[i+1+(i-1)*lda_tsqr-1]*tmp;
      }

      Tdata tmp4=0.0;
      bcast_warp(tmp,0);
      //      if(warp_id==0){printf(" tmp  = %e \n",tmp);}
      for(PARCEL_INT j=i+2+warp_id;j<=nHouse;j=j+WARP_SIZE){
	A_TSQR[j+(i+1-1)*lda_tsqr -1] = A_TSQR[j + (i+1-1)*lda_tsqr -1] -A_TSQR[j+(i-1)*lda_tsqr-1]*tmp;
	tmp4 = tmp4 + A_TSQR[j+(i+1-1)*lda_tsqr-1] * A_TSQR[j+(i+1-1)*lda_tsqr -1];
      }
      //  ここは65行目で足し合わせる
      //      Tdata tmp4_ = warpReduceSum<Tdata,Tcalc>(tmp4);
      //      tmp4 = tmp4_;

      PARCEL_INT k=i+2;
      if(k<=s){
	for(PARCEL_INT k_=kk+warp_id; k_<=(kk+nblock-1); k_=k_+ WARP_SIZE){
	  A_TSQR[k_-kk+1 + (k-1)*lda_tsqr -1] = A[k_ + (k-1)*lda -1];
	  //	  printf(" A_TSQR %d %e \n",k_,A_TSQR[k_-kk+1 + (k-1)*lda_tsqr -1]);
	}
	if(iRflag == 1){
	  for(PARCEL_INT l_=1+warp_id;l_<=s;l_=l_ + WARP_SIZE){
	    A_TSQR[l_+nblock + (k-1)*lda_tsqr-1] = R[l_+(k-1)*ldr-1];
	  }
	}
	// -----------------------------------
	PARCEL_INT ii_;
	if( (i%2) == 0){
	  ii_=1;
	}else{
	  PARCEL_INT  ii=1;
	
	  //	tmp=tmp_y1(ii)*A_TSQR(ii,k); // ここは後で足し合わせる
	  tmp = 0.0;
	  for(PARCEL_INT j=ii+1+warp_id;j<=nHouse;j=j+WARP_SIZE){
	    tmp = tmp+A_TSQR[j+(ii-1)*lda_tsqr-1]*A_TSQR[j + (k-1)*lda_tsqr-1];
	  }
	  Tdata tmp__ = warpReduceSum<Tdata,Tcalc>(tmp);

	  if(warp_id==0){
	    tmp = tmp__;
	    tmp = tmp + tmp_y1[ii-1]*A_TSQR[ii+(k-1)*lda_tsqr-1];
	    tmp = tmp_t1[ii-1]*tmp;
	    A_TSQR[ii+(k-1)*lda_tsqr-1]=A_TSQR[ii+(k-1)*lda_tsqr-1]-tmp_y1[ii-1]*tmp;
	  }
	  //	  if(warp_id==0){printf(" 2 tmp  = %e \n",tmp);}

	  bcast_warp(tmp,0);
	  for(PARCEL_INT j=ii+1+warp_id;j<=nHouse;j=j+WARP_SIZE){
	    A_TSQR[j + (k-1)*lda_tsqr -1] = A_TSQR[j + (k-1)*lda_tsqr -1] - A_TSQR[j + (ii-1)*lda_tsqr -1]*tmp;
	  }
	  ii_=2;
	}

	PARCEL_INT ii=ii_;
	//              tmp=tmp_y1(ii)*A_TSQR(ii,k) //後で足し合わせ
	tmp = 0.0;
	for(PARCEL_INT j=ii+1+warp_id; j<=nHouse;j=j+WARP_SIZE){
	  tmp = tmp + A_TSQR[j + (ii-1)*lda_tsqr -1] * A_TSQR[j+(k-1)*lda_tsqr -1];
	  //	  printf(" 3 tmp add %d %e %e \n",j,A_TSQR[j + (ii-1)*lda_tsqr -1],A_TSQR[j+(k-1)*lda_tsqr -1]);
	  //	  printf(" 3 tmp add %d %e \n",j,A_TSQR[j + (ii-1)*lda_tsqr -1]*A_TSQR[j+(k-1)*lda_tsqr -1]);
	}
	Tdata tmp___ = warpReduceSum<Tdata,Tcalc>(tmp);
	if(warp_id==0){
	  tmp = tmp___;
	  tmp = tmp + tmp_y1[ii-1] * A_TSQR[ii + (k-1)*lda_tsqr -1];
	}
	//	if(warp_id==0){printf(" 3 tmp  = %e %e %e\n",tmp,tmp_y1[ii-1] , A_TSQR[ii + (k-1)*lda_tsqr -1]);}
	
	for(PARCEL_INT ii=ii_;ii<=i;ii=ii+2){
	  tmp2 = 0.0;
	  if(warp_id==0){
	    tmp = tmp_t1[ii-1]*tmp;
	    A_TSQR[ii   + (k-1)*lda_tsqr-1] = A_TSQR[ii   + (k-1)*lda_tsqr-1] - tmp_y1[ii-1]*tmp;
	    A_TSQR[ii+1 + (k-1)*lda_tsqr-1] = A_TSQR[ii+1 + (k-1)*lda_tsqr-1] - A_TSQR[ii+1+(ii-1)*lda_tsqr-1]*tmp;
	    tmp2 = tmp_y1[ii+1-1]*A_TSQR[ii+1+(k-1)*lda_tsqr-1];
	  }
	  bcast_warp(tmp,0);

	  for(PARCEL_INT j=ii+2+warp_id;j<=nHouse;j=j+WARP_SIZE){
	    A_TSQR[j + (k-1)*lda_tsqr -1] = A_TSQR[j + (k-1)*lda_tsqr -1] - A_TSQR[j + (ii-1)*lda_tsqr -1]*tmp;
	    tmp2 = tmp2 + A_TSQR[j + (ii+1-1)*lda_tsqr -1] * A_TSQR[j + (k-1)*lda_tsqr -1];
	  }
	  tmp2_ = warpReduceSum<Tdata,Tcalc>(tmp2);
	  tmp2 = tmp2_;

	  if(warp_id==0){
	    tmp2 = tmp_t1[ii+1 -1]*tmp2;
	  }
	  
	  //	  if(warp_id==0){printf(" 2 tmp2 = %e \n",tmp);}
	  bcast_warp(tmp2,0);
	  if(ii+2<i){
	    tmp = 0.0;
	    if(warp_id==0){
	      A_TSQR[ii+1 + (k-1)*lda_tsqr -1] = A_TSQR[ii+1 + (k-1)*lda_tsqr -1] - tmp_y1[ii+1-1] * tmp2;
	      A_TSQR[ii+2 + (k-1)*lda_tsqr -1] = A_TSQR[ii+2 + (k-1)*lda_tsqr -1] - A_TSQR[ii+2 + (ii+1-1)*lda_tsqr -1]*tmp2;
	      tmp = tmp_y1[ii+2-1]*A_TSQR[ii+2+(k-1)*lda_tsqr-1];
	    }

	    for(PARCEL_INT j=ii+3+warp_id;j<=nHouse;j=j+WARP_SIZE){
	      A_TSQR[j + (k-1)*lda_tsqr -1] = A_TSQR[j + (k-1)*lda_tsqr -1] - A_TSQR[j + (ii+1-1) *lda_tsqr -1]*tmp2;
	      tmp = tmp + A_TSQR[j + (ii+2-1)*lda_tsqr -1] * A_TSQR[j + (k-1)*lda_tsqr -1];
	    }
	    tmp_ = warpReduceSum<Tdata,Tcalc>(tmp);
	    tmp = tmp_;
	    //	    if(warp_id==0){printf(" 4 tmp  = %e \n",tmp);}
	  }else{
	    if(warp_id==0){
	      A_TSQR[ii+1 + (k-1)*lda_tsqr -1 ] = A_TSQR[ ii+1 + (k-1)*lda_tsqr -1] - tmp_y1[ii+1-1] * tmp2;
	    }
	    for(PARCEL_INT j=ii+2+warp_id;j<=nHouse;j=j+WARP_SIZE){
	      A_TSQR[j + (k-1)*lda_tsqr -1] = A_TSQR[j + (k-1)*lda_tsqr -1] - A_TSQR[j + (ii+1-1)*lda_tsqr -1]*tmp2;
	    }
	  }

	}
	//              ! -----
      }
      tmp2=tmp4;

    }
  }

  for(PARCEL_INT i=1+warp_id;i<=s;i=i+WARP_SIZE){
    for(PARCEL_INT j=1;j<=s;j++){
      R[i+(j-1)*ldr-1]=0.0;
    }
  }

  for(PARCEL_INT i=1+warp_id;i<=s;i=i+WARP_SIZE){
    for(PARCEL_INT j=i;j<=s;j++){
      R[i+(j-1)*ldr-1] = A_TSQR[i + (j-1)*lda_tsqr -1];
      //      R[i+(j-1)*s-1] = 99;
    }
  }

  return ;

}


// warp単位にQR分解
// スレッド数はWARP_SIZEの倍数
// nBlock : ブロック数 = n / warp_size
template <class Tdata,class Tcalc>
__global__
void TSQR_local_Ronly_(
		       const PARCEL_INT n,
		       const PARCEL_INT s,
		       Tdata *A, PARCEL_INT lda,
		       Tdata *R, PARCEL_INT ldr,
		       Tdata *work,
		       const int nBlock ){
  extern __shared__ Tdata shared_mem[];

  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x * gridDim.x;

  const int warp_id       = threadIdx.x % WARP_SIZE;
  const int warp_block_id = 
    threadIdx.x / WARP_SIZE +  blockIdx.x * (blockDim.x  / WARP_SIZE);
  const int warp_nblock   = gridDim.x * (blockDim.x  / WARP_SIZE); // 全warp数
  const int ALL_STEP      = warp_nblock * WARP_SIZE; // CUDA_BLOCK,CUDA_THREADSにタスクを割り当てた場合のサイズ

  //  printf(" warp_block_id = %d, warp_nblock = %d \n",warp_block_id,warp_nblock);

  const int ldr_th = s*warp_nblock;
  //  Tdata *R_th   = &work[ (s*s) * warp_block_id ]; // All memory size : (s*s) *  ALL_STEP
  Tdata *R_th   = &work[ s * warp_block_id ]; // All memory size : (s*s) *  ALL_STEP

  Tdata *A_TSQR = &work[ ((s+WARP_SIZE)*s) * warp_block_id + (s*s) * warp_nblock ];
  const PARCEL_INT lda_tsqr = s+WARP_SIZE;

  for(PARCEL_INT j = 0; j<s ;j++){
    for(PARCEL_INT i = warp_id; i<s ;i=i+WARP_SIZE){
      R_th[i + j*ldr_th] = 0.0;
    }
  }
  PARCEL_INT iRflag = 0;

  for(
      PARCEL_INT j = 1 + warp_block_id * WARP_SIZE ;
      j<=n ; 
      j=j + ALL_STEP){
    PARCEL_INT jend = n > (j+WARP_SIZE -1) ? (j+WARP_SIZE -1): n;
    PARCEL_INT nQRblock=jend-j+1;
    if( (iRflag==0) && (nQRblock<=s)){
      for(PARCEL_INT i=1;i<=s;i++){
 	for(PARCEL_INT jj=j+warp_id;jj<=jend;jj=jj+WARP_SIZE){
	  R_th[1+jj-j+(i-1)*ldr_th-1] = A[jj + (i-1)*lda -1];
	}
      }
    }else{      
      House_local_Ronly_<Tdata,Tcalc>(
			 warp_id,
			 n,s,
			 A     , lda,
			 R_th  , ldr_th,
			 A_TSQR, lda_tsqr,
			 j,nQRblock,iRflag,
			 shared_mem );

    }
    iRflag=1;
  }


  return ;

}

extern "C" {

  void cuda_TSQR_local_Ronly(
			     PARCEL_INT n,
			     PARCEL_INT s,
			     devicemem *A_, PARCEL_INT offset_A, PARCEL_INT lda,
			     devicemem *R_, PARCEL_INT offset_R, PARCEL_INT ldr,
			     devicemem *work_, PARCEL_INT offset_work,
			     int numblock, 
			     int nBlock ){

    double *A    = (double*)A_[0].a_d;
    double *R    = (double*)R_[0].a_d;
    double *work = (double*)work_[0].a_d;
#if 0
    const int numGPUblock  = 1;
    const int numGPUthread = WARP_SIZE*16;
#else
    //const int numGPUblock  = 1;
    //    const int numGPUblock  = 2;
    const int numGPUblock  = numblock;
    const int numGPUthread = WARP_SIZE;
#endif
    // numGPUthread = WARP_SIZE ,numGPUblock = 1 であれば三角行列はs*s一つになるはず
    //    const size_t shared_mem= ( (s+s) * numGPUthread/WARP_SIZE )*sizeof(double);
    const size_t shared_mem= ( (s+s) * numGPUthread/WARP_SIZE )*sizeof(double);
    TSQR_local_Ronly_<double,double>
      <<<numGPUblock,numGPUthread,shared_mem>>>
      (
		       n,
		       s,
		       &A[offset_A], lda,
		       &R[offset_R], ldr,
		       work,
		       nBlock );
    // workに各warpの三角行列を持つ

#if 0
    cudaMemcpy(&R[offset_R],work,s*s*sizeof(double),cudaMemcpyDeviceToDevice);
#endif

    return ;

  }

}
