#include <stdio.h>
#include <mpi.h>
#include <cuda.h>
#include <cuda_runtime.h>

#include "cu_common.h"
#include "cu_comm.hpp"

template <class Tdata,class Tcalc>
__global__
void MatVec_ddm_(
		PARCEL_INT nx,
		PARCEL_INT ny,
		PARCEL_INT nz,
		PARCEL_INT stmx,
		PARCEL_INT stmy,
		PARCEL_INT stmz,
		PARCEL_INT mx,
		PARCEL_INT mxy,
		PARCEL_INT m,
		PARCEL_INT nnd,
 	        Tdata      *A,
 	        Tdata      *x,
 	        Tdata      *y,		
		int *ioff_A
		 ){

#if 0
  const int ip_off_0 = ioff_A[0];
  const int ip_off_1 = ioff_A[1];
  const int ip_off_2 = ioff_A[2];
  const int ip_off_3 = ioff_A[3];
  const int ip_off_4 = ioff_A[4];
  const int ip_off_5 = ioff_A[5];
  const int ip_off_6 = ioff_A[6];
#endif

#if 1
  for(
      PARCEL_INT jz = threadIdx.z + blockIdx.z*blockDim.z;
      jz<nz ;
      jz = jz + blockDim.z*gridDim.z ){
    const PARCEL_INT jz0 = mxy *( jz +stmz );
    for(
	//	PARCEL_INT jy = threadIdx.y + blockIdx.y*blockDim.y;
	PARCEL_INT jy = threadIdx.y + blockIdx.y*blockDim.y;
	jy<ny ;
	jy = jy + blockDim.y*gridDim.y ){      
      const PARCEL_INT jy0 = mx *( jy +stmy );
      const PARCEL_INT jx0 = stmx  + jy0  + jz0;
      for(
	  //	  int jx = threadIdx.x + blockIdx.x*blockDim.x;
	  //　この記述ではxのブロックサイズは常に1でないといけない
	  PARCEL_INT jx = threadIdx.x + blockIdx.x;
	  jx<nx ;
	  jx = jx + blockDim.x*gridDim.x ){

	const PARCEL_INT j = jx + jx0;
	Tdata y_ = 0.0;
#pragma unroll 3
	for(PARCEL_INT ip=0;ip<nnd;ip++){
	  const PARCEL_INT ip_off = ioff_A[ip];
	  const PARCEL_INT mip = m*ip;
	  y_ = y_ + A[ mip +j ] *x[ j + ip_off ];
	}
	y[j] = y_;
      }
    }
  }
  
#endif

#if 0
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  for(PARCEL_INT i=tid;i<(nx*ny*nz);i=i+step){
    const PARCEL_INT jz = i/nx/ny;
    const PARCEL_INT jy = (i%(nx*ny))/ny;
    const PARCEL_INT jx = i%nx;
    const PARCEL_INT j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) ;
#if 1
    Tdata y_ = 0.0;
    for(PARCEL_INT ip=0;ip<nnd;ip++){
      const PARCEL_INT ip_off = ioff_A[ip];
      const PARCEL_INT mip = m*ip;
      y_ = y_ + A[ mip +j ] *x[ j + ip_off ];
    }
#else
    Tdata y_ = 
      A[ m*0 +j ] *x[ j + ip_off_0 ] + 
      A[ m*1 +j ] *x[ j + ip_off_1 ] + 
      A[ m*2 +j ] *x[ j + ip_off_2 ] + 
      A[ m*3 +j ] *x[ j + ip_off_3 ] + 
      A[ m*4 +j ] *x[ j + ip_off_4 ] + 
      A[ m*5 +j ] *x[ j + ip_off_5 ] + 
      A[ m*6 +j ] *x[ j + ip_off_6 ];      
#endif
    y[j] = y_;
  }
#endif

#if 0
  const PARCEL_INT tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const PARCEL_INT step = blockDim.x*gridDim.x;
  for(PARCEL_INT i=tid;i<(nx*ny*nz);i=i+step){
    const PARCEL_INT jz = i/nx/ny;
    const PARCEL_INT jy = (i%(nx*ny))/ny;
    const PARCEL_INT jx = i%nx;
    const PARCEL_INT j = ( jx +stmx ) +mx *( jy +stmy )  +mxy *( jz +stmz ) ;
    y[j] = x[j]+1.0;
  }
#endif

  return ;
}

template <class Tdata,class Tcalc>
__global__
void cuda_get_ddmmatrix_neumanndiag_(
		int nx,
		int ny,
		int nz,
		int stmx,
		int stmy,
		int stmz,
		int mx,
		int mxy,
		PARCEL_INT m,
		int nnd,
		PARCEL_INT diag_offset,
 	        Tdata *A,
		Tdata *Dinv,
		int *ioff_A
				     ){
  for(
      PARCEL_INT jz = threadIdx.z + blockIdx.z*blockDim.z;
      jz<nz ;
      jz = jz + blockDim.z*gridDim.z ){
    const PARCEL_INT jz0 = mxy *( jz +stmz );
    for(
	//	PARCEL_INT jy = threadIdx.y + blockIdx.y*blockDim.y;
	PARCEL_INT jy = threadIdx.y + blockIdx.y*blockDim.y;
	jy<ny ;
	jy = jy + blockDim.y*gridDim.y ){      
      const PARCEL_INT jy0 = mx *( jy +stmy );
      const PARCEL_INT jx0 = stmx  + jy0  + jz0;
      for(
	  //	  int jx = threadIdx.x + blockIdx.x*blockDim.x;
	  PARCEL_INT jx = threadIdx.x + blockIdx.x;
	  jx<nx ;
	  jx = jx + blockDim.x*gridDim.x ){
	const PARCEL_INT j   = jx + jx0;
	const PARCEL_INT ip_diag  = diag_offset-1; // fortraとnの整合性のため-1
	const PARCEL_INT mip_diag = m*ip_diag;
	//	Dinv[j] = A[ mip +j ] == 0.0 ? 1.0 : 1.0/ A[ mip +j ];
	//	A[ mip +j ] = 0.0;
	//	printf("A diag = %f \n", A[ mip + j ]);

	Dinv[j] = 1.0 / A[ mip_diag +j ];
	// A[ mip_diag +j ] = 0.0;
	for(PARCEL_INT ip=0;ip<nnd;ip++){
	  const PARCEL_INT mip = m*ip;
	  A[ mip +j ] = -A[ mip +j ] * Dinv[j];
	}
      }

      for(
	  PARCEL_INT jx = threadIdx.x + blockIdx.x;
	  jx<nx ;
	  jx = jx + blockDim.x*gridDim.x ){
	const PARCEL_INT j = jx + jx0;

	const PARCEL_INT ip_  = nnd-1; 
	const PARCEL_INT mip_ = m*ip_;

	const PARCEL_INT ip_diag = diag_offset-1; 
	const PARCEL_INT mip_diag = m*ip_diag;
	// printf("A diag = %f \n", A[ mip_diag + j ]);
	A[ mip_diag + j ] = A[ mip_ + j ];
      }

    }
  }  
  
  return ;

}

template <class Tdata,class Tcalc>
__global__
void set_scala_ddm_(
		PARCEL_INT nx,
		PARCEL_INT ny,
		PARCEL_INT nz,
		PARCEL_INT stmx,
		PARCEL_INT stmy,
		PARCEL_INT stmz,
		PARCEL_INT mx,
		PARCEL_INT mxy,
 	        Tdata      *alpha,
 	        Tdata      *x
		 ){
  Tdata alpha_ = alpha[0];
  for(
      int jz = threadIdx.z + blockIdx.z*blockDim.z;
      jz<nz ;
      jz = jz + blockDim.z*gridDim.z ){
    const int jz0 = mxy *( jz +stmz );
    for(
	int jy = threadIdx.y + blockIdx.y*blockDim.y;
	jy<ny ;
	jy = jy + blockDim.y*gridDim.y ){      
      const int jy0 = mx *( jy +stmy );
      const int jx0 =stmx  + jy0  + jz0;
      for(
	  int jx = threadIdx.x + blockIdx.x*blockDim.x;
	  jx<nx ;
	  jx = jx + blockDim.x*gridDim.x ){
	const int j = jx + jx0;
	x[j] = alpha_;
      }
    }
  }

  return ;
}


template <class Tdata,class Tcalc>
__global__
void xcpy_ddm_(
		int nx,
		int ny,
		int nz,
		int stmx,
		int stmy,
		int stmz,
		int mx,
		int mxy,
 	        Tdata      *x,
 	        Tdata      *y
		 ){

  for(
      int jz = threadIdx.z + blockIdx.z*blockDim.z;
      jz<nz ;
      jz = jz + blockDim.z*gridDim.z ){
    const int jz0 = mxy *( jz +stmz );
    for(
	int jy = threadIdx.y + blockIdx.y*blockDim.y;
	jy<ny ;
	jy = jy + blockDim.y*gridDim.y ){      
      const int jy0 = mx *( jy +stmy );
      const int jx0 =stmx  + jy0  + jz0;
      for(
	  int jx = threadIdx.x + blockIdx.x*blockDim.x;
	  jx<nx ;
	  jx = jx + blockDim.x*gridDim.x ){
	const int j = jx + jx0;
	y[j] = x[j];
	//	printf(" %d, %f \n",j+1,y[j]);
      }
    }
  }

  return ;
}

template <class Tdata,class Tcalc>
__global__
void jacobi_iter1_ddm_(
		PARCEL_INT nx,
		PARCEL_INT ny,
		PARCEL_INT nz,
		PARCEL_INT stmx,
		PARCEL_INT stmy,
		PARCEL_INT stmz,
		PARCEL_INT mx,
		PARCEL_INT mxy,
		PARCEL_INT m,
		PARCEL_INT nnd,
 	        Tdata      *A,
 	        Tdata      *x,
 	        Tdata      *y,		
		int *ioff_A
		 ){

  for(
      PARCEL_INT jz = threadIdx.z + blockIdx.z*blockDim.z;
      jz<nz ;
      jz = jz + blockDim.z*gridDim.z ){
    const PARCEL_INT jz0 = mxy *( jz +stmz );
    for(
	PARCEL_INT jy = threadIdx.y + blockIdx.y*blockDim.y;
	jy<ny ;
	jy = jy + blockDim.y*gridDim.y ){      
      const PARCEL_INT jy0 = mx *( jy +stmy );
      const PARCEL_INT jx0 =stmx  + jy0  + jz0;
      for(
	  PARCEL_INT jx = threadIdx.x + blockIdx.x*blockDim.x;
	  jx<nx ;
	  jx = jx + blockDim.x*gridDim.x ){

	const PARCEL_INT j = jx + jx0;

	Tdata y_ = 0.0;
#pragma unroll 3
	for(PARCEL_INT ip=0;ip<nnd;ip++){
	  const PARCEL_INT ip_off = ioff_A[ip];
	  const PARCEL_INT mip    = m*ip;
	  y_ = y_ + A[ mip +j ] * x[ j + ip_off ];
	}
	//	y[j] = x[j] - y_ ;
	y[j] = x[j] + y_ ;
      }
    }
  }

  return ;
}

extern "C" {
  void cuda_MatVec_ddm(
		       int numGPUblock,
		       int numGPUthread,
		       PARCEL_INT nx,
		       PARCEL_INT ny,
		       PARCEL_INT nz,
		       PARCEL_INT stmx,
		       PARCEL_INT stmy,
		       PARCEL_INT stmz,
		       PARCEL_INT mx,
		       PARCEL_INT mxy,
		       PARCEL_INT m,
		       PARCEL_INT nnd,
		       devicemem *A,
		       devicemem *x_, PARCEL_INT offset_x,
		       devicemem *y_, PARCEL_INT offset_y,
		       devicemem *ioff_A
		       ){
    double *x = (double *)x_[0].a_d;
    double *y = (double *)y_[0].a_d;

    int nx_para = nx < 1024 ? nx : 512;
    int ny_para = ny < 1024 ? ny : 512;
    int nz_para = nz < 1024 ? nz : 512;

    if(nx<1024){
      // 予想より剰余計算が遅かった
      dim3 block(nx_para,1,1); 
      dim3 grid(1,ny_para,nz_para); 
      MatVec_ddm_<double,double> 
	//      <<<numGPUblock,numGPUthread>>>(
	<<<grid,block>>>(
			 //      <<<4096,512>>>(
			 nx,
			 ny,
			 nz,
			 stmx,
			 stmy,
			 stmz,
			 mx,
			 mxy,
			 m,
			 nnd,
			 (double *)A[0].a_d,
			 &x[offset_x],
			 &y[offset_y],
			 (int *) ioff_A[0].a_d
			 );
    }else{
      dim3 block(nx_para,1,1); 
      dim3 grid(1,ny_para,nz_para); 
      MatVec_ddm_<double,double> 
	//      <<<numGPUblock,numGPUthread>>>(
	<<<grid,block>>>(
			 //      <<<4096,512>>>(
			 nx,
			 ny,
			 nz,
			 stmx,
			 stmy,
			 stmz,
			 mx,
			 mxy,
			 m,
			 nnd,
			 (double *)A[0].a_d,
			 &x[offset_x],
			 &y[offset_y],
			 (int *) ioff_A[0].a_d
			 );
    }
    
    return ;
  }

  void cuda_get_ddmmatrix_neumanndiag(
		       int nx,
		       int ny,
		       int nz,
		       int stmx,
		       int stmy,
		       int stmz,
		       int mx,
		       int mxy,
		       PARCEL_INT m,
		       int nnd,
		       PARCEL_INT diag_offset,
		       devicemem *A,
		       devicemem *Dinv_, PARCEL_INT offset_Dinv,
		       devicemem *ioff_A
		       ){

    double *Dinv = (double *)Dinv_[0].a_d;

    int nx_para = nx < 1024 ? nx : 512;
    int ny_para = ny < 1024 ? ny : 512;
    int nz_para = nz < 1024 ? nz : 512;

    if(nx<1024){
      dim3 block(nx_para,1,1); 
      dim3 grid(1,ny_para,nz_para); 
      cuda_get_ddmmatrix_neumanndiag_<double,double> 
	<<<grid,block>>>(
			 nx,
			 ny,
			 nz,
			 stmx,
			 stmy,
			 stmz,
			 mx,
			 mxy,
			 m,
			 nnd,
			 diag_offset,
			 (double *)A[0].a_d,
			 &Dinv[offset_Dinv],
			 (int *) ioff_A[0].a_d
			 );
    }else{
      dim3 block(nx_para,1,1); 
      dim3 grid(1,ny_para,nz_para); 
      cuda_get_ddmmatrix_neumanndiag_<double,double> 
	//      <<<numGPUblock,numGPUthread>>>(
	<<<grid,block>>>(
			 nx,
			 ny,
			 nz,
			 stmx,
			 stmy,
			 stmz,
			 mx,
			 mxy,
			 m,
			 nnd,
			 diag_offset,
			 (double *)A[0].a_d,
			 &Dinv[offset_Dinv],
			 (int *) ioff_A[0].a_d
			 );
    }
    
    return ;
  }

  void cuda_set_scala_ddm(
		       int numGPUblock,
		       int numGPUthread,
		       PARCEL_INT nx,
		       PARCEL_INT ny,
		       PARCEL_INT nz,
		       PARCEL_INT stmx,
		       PARCEL_INT stmy,
		       PARCEL_INT stmz,
		       PARCEL_INT mx,
		       PARCEL_INT mxy,
		       devicemem *alpha_, PARCEL_INT offset_alpha,
		       devicemem *x_, PARCEL_INT offset_x
		       ){
    double *x     = (double *)x_[0].a_d;
    double *alpha = (double *)alpha_[0].a_d;

    int nx_para = nx < 1024 ? nx : 512;
    int ny_para = ny < 1024 ? ny : 512;
    int nz_para = nz < 1024 ? nz : 512;
    if(nx<1024){
      dim3 block(nx_para,1,1); 
      dim3 grid(1,ny_para,nz_para); 
      set_scala_ddm_<double,double> 
	<<<grid,block>>>(
			 nx,
			 ny,
			 nz,
			 stmx,
			 stmy,
			 stmz,
			 mx,
			 mxy,
			 &alpha[offset_x],
			 &x[offset_x]
			 );
    }else{
      dim3 block(nx_para,1,1); 
      dim3 grid(1,ny_para,nz_para); 
      set_scala_ddm_<double,double> 
	<<<grid,block>>>(
			 nx,
			 ny,
			 nz,
			 stmx,
			 stmy,
			 stmz,
			 mx,
			 mxy,
			 &alpha[offset_x],
			 &x[offset_x]
			 );
    }
    
    return ;
  }

  void cuda_xcpy_ddm(
		       int nx,
		       int ny,
		       int nz,
		       int stmx,
		       int stmy,
		       int stmz,
		       int mx,
		       int mxy,
		       devicemem *x_, PARCEL_INT offset_x,
		       devicemem *y_, PARCEL_INT offset_y
		       ){
    double *x = (double *)x_[0].a_d;
    double *y = (double *)y_[0].a_d;

    int nx_para = nx < 1024 ? nx : 512;
    int ny_para = ny < 1024 ? ny : 512;
    int nz_para = nz < 1024 ? nz : 512;
    if(nx<1024){
      dim3 block(nx_para,1,1); 
      dim3 grid(1,ny_para,nz_para); 
      xcpy_ddm_<double,double> 
	<<<grid,block>>>(
			 nx,
			 ny,
			 nz,
			 stmx,
			 stmy,
			 stmz,
			 mx,
			 mxy,
			 &x[offset_x],
			 &y[offset_y]
			 );
    }else{
      dim3 block(nx_para,1,1); 
      dim3 grid(1,ny_para,nz_para); 
      xcpy_ddm_<double,double> 
	<<<grid,block>>>(
			 nx,
			 ny,
			 nz,
			 stmx,
			 stmy,
			 stmz,
			 mx,
			 mxy,
			 &x[offset_x],
			 &y[offset_y]
			 );
    }
    
    return ;
  }

  void cuda_jacobi_iter1_ddm(
		       int numGPUblock,
		       int numGPUthread,
		       PARCEL_INT nx,
		       PARCEL_INT ny,
		       PARCEL_INT nz,
		       PARCEL_INT stmx,
		       PARCEL_INT stmy,
		       PARCEL_INT stmz,
		       PARCEL_INT mx,
		       PARCEL_INT mxy,
		       PARCEL_INT m,
		       PARCEL_INT nnd,
		       devicemem *A,
		       devicemem *x_,PARCEL_INT offset_x,
		       devicemem *y_,PARCEL_INT offset_y,
		       devicemem *ioff_A
		       ){
    double *x = (double *)x_[0].a_d;
    double *y = (double *)y_[0].a_d;

    int nx_para = nx < 1024 ? nx : 512;
    int ny_para = ny < 1024 ? ny : 512;
    int nz_para = nz < 1024 ? nz : 512;

    if(nx<1024){
      dim3 block(nx_para,1,1); 
      dim3 grid(1,ny_para,nz_para); 
      jacobi_iter1_ddm_<double,double> 
	//      <<<numGPUblock,numGPUthread>>>(
	<<<grid,block>>>(
			 //      <<<4096,512>>>(
			 nx,
			 ny,
			 nz,
			 stmx,
			 stmy,
			 stmz,
			 mx,
			 mxy,
			 m,
			 nnd,
			 (double *)A[0].a_d,
			 &x[offset_x],
			 &y[offset_y],
			 (int *) ioff_A[0].a_d
			 );
    }else{
      dim3 block(nx_para,1,1); 
      dim3 grid(1,ny_para,nz_para); 
      jacobi_iter1_ddm_<double,double> 
	//      <<<numGPUblock,numGPUthread>>>(
	<<<grid,block>>>(
			 //      <<<4096,512>>>(
			 nx,
			 ny,
			 nz,
			 stmx,
			 stmy,
			 stmz,
			 mx,
			 mxy,
			 m,
			 nnd,
			 (double *)A[0].a_d,
			 &x[offset_x],
			 &y[offset_y],
			 (int *) ioff_A[0].a_d
			 );
    }
    
    return ;
  }


}

// DDM形式はint型のみ対応
template <class Tdata,class Tcalc>
__global__ 
void set_sbuf_DDM_(
		   Tdata *vec, 
		  Tdata *sbuf, 
		  int offset,
		  int grids,
		  int mx,int mxy,
		  int x0s,int y0s,int z0s,
		  int xlen,int ylen,
		  int stmx,int stmy,int stmz
		  ){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;

  for(PARCEL_INT j=tid;j<grids;j=j+step){
       const PARCEL_INT ix = x0s + j%xlen;
       const PARCEL_INT iy = y0s + (j/xlen)%ylen;
       const PARCEL_INT iz = z0s + j / (xlen * ylen);
       const PARCEL_INT k = (ix+stmx) + mx*(iy+stmy) + mxy*(iz+stmz);
       sbuf[j+offset-1] = vec[k];
  }

}

template <class Tdata,class Tcalc>
__global__ 
void set_rbuf_DDM_(
		  Tdata *vec,
		  Tdata *rbuf,
		  int offset,
		  int grids,
		  int mx,int mxy,
		  int x0r,int y0r,int z0r,
		  int xlen,int ylen,
		  int stmx,int stmy,int stmz
		  ){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;

  for(PARCEL_INT j=tid;j<grids;j=j+step){
       const PARCEL_INT ix = x0r + j%xlen;
       const PARCEL_INT iy = y0r + (j/xlen)%ylen;
       const PARCEL_INT iz = z0r + j / (xlen * ylen);
       const PARCEL_INT k = (ix+stmx) + mx*(iy+stmy) + mxy*(iz+stmz);
       vec[k] = rbuf[j+offset-1] ;
  }

}

template <class Tdata,class Tcalc>
__global__ 
void sleev_zero_DDM_(
		  Tdata *vec,
		  int offset,
		  int grids,
		  int mx,int mxy,
		  int x0r,int y0r,int z0r,
		  int xlen,int ylen,
		  int stmx,int stmy,int stmz
		  ){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;

  for(PARCEL_INT j=tid;j<grids;j=j+step){
       const PARCEL_INT ix = x0r + j%xlen;
       const PARCEL_INT iy = y0r + (j/xlen)%ylen;
       const PARCEL_INT iz = z0r + j / (xlen * ylen);
       const PARCEL_INT k = (ix+stmx) + mx*(iy+stmy) + mxy*(iz+stmz);
       vec[k] = 0.0;
  }

}

extern "C" {
  
  void cuda_set_sbuf_DDM_FP64_GPU(
			 int numGPUblock,
			 int numGPUthread,
			 devicemem *vec_, PARCEL_INT offset_vec,
			 devicemem *sbuf_, PARCEL_INT offset_sbuf,
			 int offset,
			 int grids,
			 int mx,int mxy,
			 int x0s,int y0s,int z0s,
			 int xlen,int ylen,
			 int stmx,int stmy,int stmz
			 ){
    double *vec  = (double *)vec_[0].a_d;
    double *sbuf = (double *)sbuf_[0].a_d;

    set_sbuf_DDM_<double,double> 
      <<<numGPUblock,numGPUthread>>>(
				     &vec[offset_vec],
				     &sbuf[offset_sbuf],
				     offset,
				     grids,
				     mx, mxy,
				     x0s, y0s, z0s,
				     xlen, ylen,
				     stmx, stmy, stmz);
      }

  void cuda_set_rbuf_DDM_FP64_GPU(
			 int numGPUblock,
			 int numGPUthread,
			 devicemem *vec_, PARCEL_INT offset_vec,
			 devicemem *rbuf_, PARCEL_INT offset_rbuf,
			 int offset,
			 int grids,
			 int mx,int mxy,
			 int x0r,int y0r,int z0r,
			 int xlen,int ylen,
			 int stmx,int stmy,int stmz
			 ){
    double *vec = (double *)vec_[0].a_d;
    double *rbuf = (double *)rbuf_[0].a_d;

    set_rbuf_DDM_<double,double> 
      <<<numGPUblock,numGPUthread>>>(
				     &vec[offset_vec],
				     &rbuf[offset_rbuf],
				     offset,
				     grids,
				     mx, mxy,
				     x0r, y0r, z0r,
				     xlen, ylen,
				     stmx, stmy, stmz);
      }
  
  void cuda_sleev_zero_DDM_FP64_GPU(
			 int numGPUblock,
			 int numGPUthread,
			 devicemem *vec_,PARCEL_INT offset_vec,
			 int offset, 
			 int grids,
			 int mx,int mxy,
			 int x0r,int y0r,int z0r,
			 int xlen,int ylen,
			 int stmx,int stmy,int stmz
			 ){

    double *vec = (double *)vec_[0].a_d;
    sleev_zero_DDM_<double,double> 
      <<<numGPUblock,numGPUthread>>>(
				     &vec[offset_vec],
				     offset,
				     grids,
				     mx, mxy,
				     x0r, y0r, z0r,
				     xlen, ylen,
				     stmx, stmy, stmz);
      }

#if 0
  void comm_ddm(
		int n,
		double *sbuf,
		double *rbuf,
		int *grids,
		int myrank,
		int *rank_s,int *rank_r,
		hostmem *comm_,
		hostmem *s_stat_,hostmem *r_stat_,
		hostmem *s_req_ ,hostmem *r_req_
		){
    MPI_Status  stat_send[6],stat_recv[6];
    MPI_Request req_send[6], req_recv[6];
    MPI_Comm *comm   = (MPI_Comm *)comm_[0].a_d;

    //    int myrank_;
    //    MPI_Comm_rank(comm[0],&myrank_);
 
    for(PARCEL_INT i = 0; i<6;i++){
#if 1
      MPI_Request *s_req = &req_send[i];
      MPI_Request *r_req = &req_recv[i];
      MPI_Status *s_stat = &stat_send[i];
      MPI_Status *r_stat = &stat_recv[i];
#else
      MPI_Request *s_req = (MPI_Request *)s_req_[i].a_d; 
      MPI_Request *r_req = (MPI_Request *)r_req_[i].a_d; 
      MPI_Status *s_stat = (MPI_Status *)s_stat_[i].a_d;
      MPI_Status *r_stat = (MPI_Status *)r_stat_[i].a_d;
#endif
      if ( rank_s[i] != -1 ) {
	//	printf(" send i = %d : %d --> %d :: %d \n",i,myrank_,rank_s[i],grids[i]);
	MPI_Isend(
                  &sbuf[i*n],
		  grids[i], 
		  //1,
		  MPI_DOUBLE,
		  rank_s[i],
                  rank_s[i],
		  //MPI_COMM_WORLD,
		  comm[0],
		  //		  &req_send[i]
		  s_req
                  );
      }
      if ( rank_r[i] != -1 ) {
	//	printf(" recv i = %d : %d <-- %d :: %d \n",i,myrank_,rank_r[i],grids[i]);
	MPI_Irecv(
		  &rbuf[i*n],
		  grids[i], 
		  //1,
		  MPI_DOUBLE,
		  rank_r[i],
                  myrank,
		  //MPI_COMM_WORLD,
		  comm[0],
		  //		  &req_recv[i]
		  r_req
                  );
      }

    }

    for(PARCEL_INT i = 0; i<6;i++){
      if ( rank_s[i] != -1 ) {
	MPI_Wait(&req_send[i],&stat_send[i]);
      }
      if ( rank_r[i] != -1 ) {
	MPI_Wait(&req_recv[i],&stat_recv[i]);
      }
    }

  }
#endif

  void comm_ddm_FP64_GPU(
		int n,
		devicemem *sbuf_, PARCEL_INT offset_sbuf,
		devicemem *rbuf_, PARCEL_INT offset_rbuf,
		int *grids,
		int myrank,
		int *rank_s,int *rank_r,
		hostmem *comm_,
		hostmem *s_stat_,hostmem *r_stat_,
		hostmem *s_req_ ,hostmem *r_req_
		){
    MPI_Status  stat_send[6],stat_recv[6];
    MPI_Request req_send[6], req_recv[6];
    MPI_Comm *comm   = (MPI_Comm *)comm_[0].a_d;

    //    int myrank_;
    //    MPI_Comm_rank(comm[0],&myrank_);
 
    double *sbuf = (double *)sbuf_[0].a_d;
    double *rbuf = (double *)rbuf_[0].a_d;

    // int grids_max = 0;
    // for(int i = 0; i<6;i++){
    //   grids_max = grids_max + grids[i];
    // }
    // double *sbuf_cpu = (double*)malloc(sizeof(double)*grids_max);
    // double *rbuf_cpu = (double*)malloc(sizeof(double)*grids_max);
    // for(int i = 0; i<grids_max;i++){
    //   sbuf_cpu[i] = 0.0;
    //   rbuf_cpu[i] = 0.0;
    // }
    // const size_t len = sizeof(double) * grids_max;
    // cudaMemcpy(sbuf_cpu, sbuf, len, cudaMemcpyDeviceToHost);
    // cudaMemcpy(rbuf_cpu, rbuf, len, cudaMemcpyDeviceToHost);

    //    printf("start comm_ddm_FP64_GPU \n");
    PARCEL_INT offset = 0;
    for(PARCEL_INT i = 0; i<6;i++){

#if 1
      MPI_Request *s_req = &req_send[i];
      MPI_Request *r_req = &req_recv[i];
      MPI_Status *s_stat = &stat_send[i];
      MPI_Status *r_stat = &stat_recv[i];
#else
      MPI_Request *s_req = (MPI_Request *)s_req_[i].a_d; 
      MPI_Request *r_req = (MPI_Request *)r_req_[i].a_d; 
      MPI_Status *s_stat = (MPI_Status *)s_stat_[i].a_d;
      MPI_Status *r_stat = (MPI_Status *)r_stat_[i].a_d;
#endif

      if ( rank_s[i] != -1 ) {
	MPI_Isend(
                  &sbuf[offset+offset_sbuf],
		  grids[i], 
		  //1,
		  MPI_DOUBLE,
		  rank_s[i],
                  rank_s[i],
		  //MPI_COMM_WORLD,
		  comm[0],
		  //		  &req_send[i]
		  s_req
                  );
      }
      if ( rank_r[i] != -1 ) {
	MPI_Irecv(
                  &rbuf[offset+offset_rbuf],
		  grids[i], 
		  MPI_DOUBLE,
		  rank_r[i],
                  myrank,
		  comm[0],
		  r_req
                  );
      }
      offset = offset + grids[i];
    }

    for(PARCEL_INT i = 0; i<6;i++){
      if ( rank_s[i] != -1 ) {
	MPI_Wait(&req_send[i],&stat_send[i]);
      }
      if ( rank_r[i] != -1 ) {
	MPI_Wait(&req_recv[i],&stat_recv[i]);
      }
    }

    // cudaMemcpy(sbuf, sbuf_cpu, len, cudaMemcpyHostToDevice);
    // cudaMemcpy(rbuf, rbuf_cpu, len, cudaMemcpyHostToDevice);
    // free(sbuf_cpu);
    // free(rbuf_cpu);
  }

  void comm_ddm_FP64_GPU_host(
		int n,
		devicemem *sbuf_, PARCEL_INT offset_sbuf,
		devicemem *rbuf_, PARCEL_INT offset_rbuf,
		int *grids,
		int myrank,
		int *rank_s,int *rank_r,
		hostmem *comm_,
		hostmem *s_stat_,hostmem *r_stat_,
		hostmem *s_req_ ,hostmem *r_req_
		){
    MPI_Status  stat_send[6],stat_recv[6];
    MPI_Request req_send[6], req_recv[6];
    MPI_Comm *comm   = (MPI_Comm *)comm_[0].a_d;

    //    int myrank_;
    //    MPI_Comm_rank(comm[0],&myrank_);
 
    double *sbuf = (double *)sbuf_[0].a_d;
    double *rbuf = (double *)rbuf_[0].a_d;

    int grids_max = 0;
    for(int i = 0; i<6;i++){
      grids_max = grids_max + grids[i];
    }
    double *sbuf_cpu = (double*)malloc(sizeof(double)*grids_max);
    double *rbuf_cpu = (double*)malloc(sizeof(double)*grids_max);
    for(int i = 0; i<grids_max;i++){
      sbuf_cpu[i] = 0.0;
      rbuf_cpu[i] = 0.0;
    }
    const size_t len = sizeof(double) * grids_max;
    cudaMemcpy(sbuf_cpu, sbuf, len, cudaMemcpyDeviceToHost);

    PARCEL_INT offset = 0;
    for(PARCEL_INT i = 0; i<6;i++){

      MPI_Request *s_req = &req_send[i];
      MPI_Request *r_req = &req_recv[i];
      MPI_Status *s_stat = &stat_send[i];
      MPI_Status *r_stat = &stat_recv[i];

      if ( rank_s[i] != -1 ) {
	MPI_Isend(
                  &sbuf_cpu[offset+offset_sbuf],
		  grids[i], 
		  //1,
		  MPI_DOUBLE,
		  rank_s[i],
                  rank_s[i],
		  //MPI_COMM_WORLD,
		  comm[0],
		  //		  &req_send[i]
		  s_req
                  );
      }
      if ( rank_r[i] != -1 ) {
	MPI_Irecv(
                  &rbuf_cpu[offset+offset_rbuf],
		  grids[i], 
		  MPI_DOUBLE,
		  rank_r[i],
                  myrank,
		  comm[0],
		  r_req
                  );
      }
      offset = offset + grids[i];
    }

    for(PARCEL_INT i = 0; i<6;i++){
      if ( rank_s[i] != -1 ) {
	MPI_Wait(&req_send[i],&stat_send[i]);
      }
      if ( rank_r[i] != -1 ) {
	MPI_Wait(&req_recv[i],&stat_recv[i]);
      }
    }

    cudaMemcpy(rbuf, rbuf_cpu, len, cudaMemcpyHostToDevice);
    free(sbuf_cpu);
    free(rbuf_cpu);
  }


}
