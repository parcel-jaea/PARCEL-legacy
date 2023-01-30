#include <stdio.h>
#include <mpi.h>
#include <cuda.h>
#include <cuda_runtime.h>

#include "cu_common.h"
#include "cu_comm.hpp"
#include "cu_common.hpp"

#define loop_single 0

template <class Tdata,class Tcalc>
__global__
void solve_DILU_DDM_(
		     Tdata *A,
		     Tdata *Dinv,
		     Tdata *r,
		     Tdata *s,
		     Tdata *svec,
		     int *ioff_dia,
		     int *offset_xyz,
		     int *nxe_,
		     int *nye_,
		     int *nze_,
		     int nx,
		     int ny,
		     int nz,
		     int stmx,
		     int stmy,
		     int stmz,
		     int mx,
		     int my,
		     int mz,
		     int nxblock,
		     int nyblock,
		     int nzblock,
		     int nx_local,
		     int ny_local,
		     int nz_local,
		     int mx_local,
		     int my_local,
		     int mz_local,
		     int nnd){
  const int mxy = mx*my;
  const int m   = mx * my * mz;
  const int m_local  = mx_local * my_local * mz_local;

#if 0
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
      const int jx0 = stmx  + jy0  + jz0;
      for(
	  int jx = threadIdx.x + blockIdx.x*blockDim.x;
	  jx<nx ;
	  jx = jx + blockDim.x*gridDim.x ){
	const int j = jx + jx0;
	s[j] = r[j];
      }
    }
  }
    svec[0] = 1.0;
  return ;
#endif
  // for(int ip=1;ip<=nnd;ip++){
  //   printf("ip = %d : %d \n",ip, ioff_dia[ip-1]);
  // }
  // return ;

  for(
      int jzblock = 1+threadIdx.z + blockIdx.z*blockDim.z;
      jzblock<= nzblock ;
      jzblock = jzblock + blockDim.z*gridDim.z ){
    const int nzs = 1 + nz_local*(jzblock-1);
    const int nze = nze_[jzblock-1];
    for(
	int jyblock = 1+threadIdx.y + blockIdx.y*blockDim.y;
	jyblock<= nyblock ;
	jyblock = jyblock + blockDim.y*gridDim.y ){
      const int nys = 1 + ny_local*(jyblock-1);
      const int nye = nye_[jyblock-1];
      for(
	  int jxblock = 1+threadIdx.x + blockIdx.x*blockDim.x;
	  jxblock<= nxblock ;
	  jxblock = jxblock + blockDim.x*gridDim.x ){
	const int nxs = 1 + nx_local*(jxblock-1);
	const int nxe = nxe_[jxblock-1];
	const int block_id = 
	  (jxblock-1) + (jyblock-1)*nxblock + (jzblock-1)*nxblock*nyblock ;

	const int nblock = nxblock*nyblock*nzblock;
	//前進代入
	for(PARCEL_INT jz=nzs-1;jz<=(nze-1);jz++){
	  for(PARCEL_INT jy=nys-1;jy<=(nye-1);jy++){
	    for(PARCEL_INT jx=nxs-1;jx<=(nxe-1);jx++){
	      const PARCEL_INT j = (jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1;

	      const PARCEL_INT jjx = jx - (nxs-1);
	      const PARCEL_INT jjy = jy - (nys-1);
	      const PARCEL_INT jjz = jz - (nzs-1);
	      const PARCEL_INT jj  = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1;

	      Tdata zz = r[j-1]; 
	      for(PARCEL_INT ip=1;ip<=nnd;ip++){
		const PARCEL_INT offset   = ioff_dia[ip-1];
		const PARCEL_INT offsetx  = offset_xyz[ip+nnd*0-1];
		const PARCEL_INT offsety  = offset_xyz[ip+nnd*1-1];
		const PARCEL_INT offsetz  = offset_xyz[ip+nnd*2-1];
		const PARCEL_INT offsetjj = offsetx + offsety*mx_local + offsetz*mx_local*my_local;
	      
		if(offset>=0){
		  continue;
		} 
		//		zz = zz - svec[jj + m_local * block_id+offsetjj-1] * A[m*(ip-1)+j-1];
		zz = zz - svec[block_id+(jj + offsetjj)*nblock -1] * A[m*(ip-1)+j-1];
	      }
	      //	      svec[jj + m_local * block_id -1] = Dinv[jj + m_local * block_id -1]*zz;
	      svec[block_id + jj*nblock -1] = Dinv[block_id + jj*nblock -1]*zz;
	    }
	  }
	}


	//! 後退代入
	for(PARCEL_INT jz=(nze-1);jz>=nzs-1;jz--){
	  for(PARCEL_INT jy=(nye-1);jy>=nys-1;jy--){
	    for(PARCEL_INT jx=(nxe-1);jx>=nxs-1;jx--){
	      const PARCEL_INT j  = (jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1;
	      const PARCEL_INT jjx = jx - (nxs-1);
	      const PARCEL_INT jjy = jy - (nys-1);
	      const PARCEL_INT jjz = jz - (nzs-1);
	      const PARCEL_INT jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1;

	      Tdata zz = 0.0;
	      for(PARCEL_INT ip=1;ip<=nnd;ip++){
		const PARCEL_INT offset   = ioff_dia[ip-1];
		const PARCEL_INT offsetx  = offset_xyz[ip+nnd*0-1];
		const PARCEL_INT offsety  = offset_xyz[ip+nnd*1-1];
		const PARCEL_INT offsetz  = offset_xyz[ip+nnd*2-1];
		const PARCEL_INT offsetjj = offsetx + offsety*mx_local + offsetz*mx_local*my_local;
		if(offset<=0){
		  continue;
		}
		//		zz = zz + svec[jj+m_local*block_id+offsetjj-1] * A[m*(ip-1)+j-1];
		zz = zz + svec[ block_id + (jj+offsetjj)*nblock -1 ] * A[m*(ip-1)+j-1];
	      }
	      //	      svec[jj + m_local * block_id -1] = 
	      //		svec[jj + m_local * block_id -1] - 
	      //		Dinv[jj + m_local * block_id -1]*zz;
	      svec[block_id + jj*nblock -1] = 
		svec[block_id + jj*nblock -1] - 
		Dinv[block_id + jj*nblock -1]*zz;

	      //	      s[j-1] = svec[jj + m_local * block_id-1];
	      s[j-1] = svec[block_id + jj*nblock -1];
	    }
	  }
	}

      }
    }
  }

  return ;
}

template <class Tdata,class Tcalc>
__global__
void solve_DILU_AS_DDM_setvec_(
			Tdata *AS_v1,
			Tdata *s,
			Tdata *rvec,Tdata *svec,Tdata *svec_AS,
			Tdata *AS_v2,
			Tdata *AS_val_dia_init, // 行列データ
			Tdata *Dinv,            // 対角要素の近似逆数
			int *ioff_dia,int *offset_xyz,
			int nnd, 
			int *nxs_, int *nxe_,
			int *nys_, int *nye_,
			int *nzs_, int *nze_,
			int *nxs0_, int *nxe0_,
			int *nys0_, int *nye0_,
			int *nzs0_, int *nze0_,
			int addL,
			int nblock,
			int nxblock,int nyblock,int nzblock,
			int stmx,int stmy,int stmz,
			int mx, int mxy,
			int mx_local,int my_local,
			int m_AS,
			int mx_AS ,int mxy_AS,
			int sleeve_mx,int sleeve_my,int sleeve_mz
			){

#if loop_single
  for(PARCEL_INT jzblock=1; jzblock<=nzblock;jzblock++){
    for(PARCEL_INT jyblock=1; jyblock<=nyblock;jyblock++){
      for(PARCEL_INT jxblock=1; jxblock<=nxblock;jxblock++){
#else
  for(
      PARCEL_INT jzblock = 1+threadIdx.z + blockIdx.z*blockDim.z;
      jzblock<= nzblock ;
      jzblock = jzblock + blockDim.z*gridDim.z ){
    for(
	PARCEL_INT jyblock = 1+threadIdx.y + blockIdx.y*blockDim.y;
	jyblock<= nyblock ;
	jyblock = jyblock + blockDim.y*gridDim.y ){
      for(
	  PARCEL_INT jxblock = 1+threadIdx.x + blockIdx.x*blockDim.x;
	  jxblock<= nxblock ;
	  jxblock = jxblock + blockDim.x*gridDim.x ){
#endif
	PARCEL_INT nzs = nzs_[jzblock-1];
	PARCEL_INT nze = nze_[jzblock-1];
	PARCEL_INT nys = nys_[jyblock-1];
	PARCEL_INT nye = nye_[jyblock-1];
	PARCEL_INT nxs = nxs_[jxblock-1];
	PARCEL_INT nxe = nxe_[jxblock-1];
	PARCEL_INT block_id = (jxblock-1) + (jyblock-1)*nxblock + (jzblock-1)*nxblock*nyblock ;
	for(PARCEL_INT jz = nzs-1; jz<=nze-1;jz++){
	  for(PARCEL_INT jy = nys-1; jy<=nye-1;jy++){
	    for(PARCEL_INT jx = nxs-1; jx<=nxe-1;jx++){
	      PARCEL_INT j  = (jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1;
	      PARCEL_INT j_AS = 
		(jx + sleeve_mx) + 
		(jy + sleeve_my) * mx_AS  + 
		(jz + sleeve_mz) * mxy_AS + 
		1;
	      PARCEL_INT jjx = jx - (nxs-1);
	      PARCEL_INT jjy = jy - (nys-1);
	      PARCEL_INT jjz = jz - (nzs-1);
	      PARCEL_INT jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1;
	      rvec[block_id + jj*nblock-1] = AS_v1[j_AS-1];
	    }
	  }
	}	
      }
    }
  }
  return ;
 }

template <class Tdata,class Tcalc>
__global__
void solve_DILU_AS_DDM_solve_(
			Tdata *AS_v1,
			Tdata *s,
			Tdata *rvec,Tdata *svec,Tdata *svec_AS,
			Tdata *AS_v2,
			Tdata *AS_val_dia_init, // 行列データ
			Tdata *Dinv,            // 対角要素の近似逆数
			int *ioff_dia,int *offset_xyz,
			int nnd, 
			int *nxs_, int *nxe_,
			int *nys_, int *nye_,
			int *nzs_, int *nze_,
			int *nxs0_, int *nxe0_,
			int *nys0_, int *nye0_,
			int *nzs0_, int *nze0_,
			int addL,
			int nblock,
			int nxblock,int nyblock,int nzblock,
			int stmx,int stmy,int stmz,
			int mx, int mxy,
			int mx_local,int my_local,
			int m_AS,
			int mx_AS ,int mxy_AS,
			int sleeve_mx,int sleeve_my,int sleeve_mz
			){

#if loop_single
  for(PARCEL_INT jzblock=1; jzblock<=nzblock;jzblock++){
    for(PARCEL_INT jyblock=1; jyblock<=nyblock;jyblock++){
      for(PARCEL_INT jxblock=1; jxblock<=nxblock;jxblock++){
#else
  for(
      PARCEL_INT jzblock = 1+threadIdx.z + blockIdx.z*blockDim.z;
      jzblock<= nzblock ;
      jzblock = jzblock + blockDim.z*gridDim.z ){
    for(
	PARCEL_INT jyblock = 1+threadIdx.y + blockIdx.y*blockDim.y;
	jyblock<= nyblock ;
	jyblock = jyblock + blockDim.y*gridDim.y ){
      for(
	  PARCEL_INT jxblock = 1+threadIdx.x + blockIdx.x*blockDim.x;
	  jxblock<= nxblock ;
	  jxblock = jxblock + blockDim.x*gridDim.x ){
#endif
	PARCEL_INT nzs = nzs_[jzblock-1];
	PARCEL_INT nze = nze_[jzblock-1];
	PARCEL_INT nys = nys_[jyblock-1];
	PARCEL_INT nye = nye_[jyblock-1];
	PARCEL_INT nxs = nxs_[jxblock-1];
	PARCEL_INT nxe = nxe_[jxblock-1];
	PARCEL_INT block_id = (jxblock-1) + (jyblock-1)*nxblock + (jzblock-1)*nxblock*nyblock ;
	for(PARCEL_INT jz = nzs-1; jz<=nze-1;jz++){
	  for(PARCEL_INT jy = nys-1; jy<=nye-1;jy++){
	    for(PARCEL_INT jx = nxs-1; jx<=nxe-1;jx++){
	      PARCEL_INT j  = (jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1;
	      PARCEL_INT j_AS = 
		(jx + sleeve_mx) + 
		(jy + sleeve_my) * mx_AS  + 
		(jz + sleeve_mz) * mxy_AS + 
		1;
	      PARCEL_INT jjx = jx - (nxs-1);
	      PARCEL_INT jjy = jy - (nys-1);
	      PARCEL_INT jjz = jz - (nzs-1);
	      PARCEL_INT jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1;
	      rvec[block_id + jj*nblock-1] = AS_v1[j_AS-1];
	    }
	  }
	}	
      }
    }
  }

  // -----
  // (D+L)D^-1(D+U)s=r
  // -----
  //  1. (D+L)z=r       
  //     z = Dinv * (r - Lz)
  //  2. (I+D^-1U)s=z
  //     s = z-D^-1 Us        
#if loop_single
  for(PARCEL_INT jzblock=1;jzblock<=nzblock;jzblock++){
    for(PARCEL_INT jyblock=1;jyblock<=nyblock;jyblock++){
      for(PARCEL_INT jxblock=1;jxblock<=nxblock;jxblock++){
#else
  for(
      PARCEL_INT jzblock = 1+threadIdx.z + blockIdx.z*blockDim.z;
      jzblock<= nzblock ;
      jzblock = jzblock + blockDim.z*gridDim.z ){
    for(
	PARCEL_INT jyblock = 1+threadIdx.y + blockIdx.y*blockDim.y;
	jyblock<= nyblock ;
	jyblock = jyblock + blockDim.y*gridDim.y ){
      for(
	  PARCEL_INT jxblock = 1+threadIdx.x + blockIdx.x*blockDim.x;
	  jxblock<= nxblock ;
	  jxblock = jxblock + blockDim.x*gridDim.x ){
#endif
	PARCEL_INT nzs = nzs_[jzblock-1];
	PARCEL_INT nze = nze_[jzblock-1];
	PARCEL_INT nys = nys_[jyblock-1];
	PARCEL_INT nye = nye_[jyblock-1];
	PARCEL_INT nxs = nxs_[jxblock-1];
	PARCEL_INT nxe = nxe_[jxblock-1];
	PARCEL_INT block_id = (jxblock-1) + (jyblock-1)*nxblock + (jzblock-1)*nxblock*nyblock ;
	  
	// ! 前進代入
	for(PARCEL_INT jz = nzs-1; jz<=nze-1;jz++){
	  for(PARCEL_INT jy = nys-1; jy<=nye-1;jy++){
	    for(PARCEL_INT jx = nxs-1; jx<=nxe-1;jx++){
	      PARCEL_INT j  = (jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1;
	      PARCEL_INT j_AS = 
		(jx + sleeve_mx) + 
		(jy + sleeve_my) * mx_AS  + 
		(jz + sleeve_mz) * mxy_AS + 
		1;
	      PARCEL_INT jjx = jx - (nxs-1);
	      PARCEL_INT jjy = jy - (nys-1);
	      PARCEL_INT jjz = jz - (nzs-1);
	      PARCEL_INT jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1;
	      Tdata zz = rvec[block_id + jj*nblock -1];
	      for(PARCEL_INT ip = 1;ip<=nnd;ip++){
		PARCEL_INT offset   = ioff_dia[ip-1];
		PARCEL_INT offsetx  = offset_xyz[ip+nnd*0-1];
		PARCEL_INT offsety  = offset_xyz[ip+nnd*1-1];
		PARCEL_INT offsetz  = offset_xyz[ip+nnd*2-1];
		PARCEL_INT offsetjj = offsetx + offsety*mx_local + offsetz*mx_local*my_local;
		if(offset>=0){
		  continue;
		}
		zz = zz 
		  - svec[block_id + (jj+offsetjj)*nblock-1] * 
		  AS_val_dia_init[m_AS*(ip-1)+j_AS-1];		
	      }
	      svec[block_id + jj*nblock-1] = Dinv[block_id + jj*nblock-1]*zz;
	    }
	  }
	}

	//! 後退代入
	for(PARCEL_INT jz = nze-1; jz>=nzs-1; jz--){
	  for(PARCEL_INT jy = nye-1; jy>=nys-1; jy--){
	    for(PARCEL_INT jx = nxe-1; jx>=nxs-1; jx--){

	      PARCEL_INT j  = (jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1;
	      PARCEL_INT j_AS = 
		(jx + sleeve_mx) + 
		(jy + sleeve_my) * mx_AS  + 
		(jz + sleeve_mz) * mxy_AS + 
		1;
	      PARCEL_INT jjx = jx - (nxs-1);
	      PARCEL_INT jjy = jy - (nys-1);
	      PARCEL_INT jjz = jz - (nzs-1);
	      PARCEL_INT jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1;

	      Tdata zz = 0.0;
	      for(PARCEL_INT ip = 1;ip<=nnd;ip++){
		PARCEL_INT offset   = ioff_dia[ip-1];
		PARCEL_INT offsetx  = offset_xyz[ip+nnd*0-1];
		PARCEL_INT offsety  = offset_xyz[ip+nnd*1-1];
		PARCEL_INT offsetz  = offset_xyz[ip+nnd*2-1];
		PARCEL_INT offsetjj = offsetx + offsety*mx_local + offsetz*mx_local*my_local;
		if(offset<=0){
		  continue;
		}
		zz = zz + svec[block_id + (jj+offsetjj)*nblock-1] * AS_val_dia_init[m_AS*(ip-1)+j_AS-1];
	      }
	      svec[block_id + jj*nblock-1] =  svec[block_id + jj*nblock-1] - Dinv[block_id + jj*nblock-1]*zz;
	    }
	  }
	}

      }
    }
  }

  return ;
 }

#if 1
template <class Tdata,class Tcalc>
__global__
void solve_DILU_AS_DDM_update_(
			Tdata *AS_v1,
			Tdata *s,
			Tdata *rvec,Tdata *svec,Tdata *svec_AS,
			Tdata *AS_v2,
			Tdata *AS_val_dia_init, // 行列データ
			Tdata *Dinv,            // 対角要素の近似逆数
			int *ioff_dia,int *offset_xyz,
			int nnd, 
			int *nxs_, int *nxe_,
			int *nys_, int *nye_,
			int *nzs_, int *nze_,
			int *nxs0_, int *nxe0_,
			int *nys0_, int *nye0_,
			int *nzs0_, int *nze0_,
			int addL,
			int nblock,
			int nxblock,int nyblock,int nzblock,
			int stmx,int stmy,int stmz,
			int mx, int mxy,
			int mx_local,int my_local,
			int m_AS,
			int mx_AS ,int mxy_AS,
			int sleeve_mx,int sleeve_my,int sleeve_mz
			){

#if loop_single
  for(PARCEL_INT jzblock=1;jzblock<=nzblock;jzblock++){
    for(PARCEL_INT jyblock=1;jyblock<=nyblock;jyblock++){
      for(PARCEL_INT jxblock=1;jxblock<=nxblock;jxblock++){
#else
  for(
      PARCEL_INT jzblock = 1+threadIdx.z + blockIdx.z*blockDim.z;
      jzblock<= nzblock ;
      jzblock = jzblock + blockDim.z*gridDim.z ){
    for(
	PARCEL_INT jyblock = 1+threadIdx.y + blockIdx.y*blockDim.y;
	jyblock<= nyblock ;
	jyblock = jyblock + blockDim.y*gridDim.y ){
      for(
	  PARCEL_INT jxblock = 1+threadIdx.x + blockIdx.x*blockDim.x;
	  jxblock<= nxblock ;
	  jxblock = jxblock + blockDim.x*gridDim.x ){
#endif
	PARCEL_INT nzs = nzs_[jzblock-1];
	PARCEL_INT nze = nze_[jzblock-1];
	PARCEL_INT nys = nys_[jyblock-1];
	PARCEL_INT nye = nye_[jyblock-1];
	PARCEL_INT nxs = nxs_[jxblock-1];
	PARCEL_INT nxe = nxe_[jxblock-1];
	PARCEL_INT block_id = (jxblock-1) + (jyblock-1)*nxblock + (jzblock-1)*nxblock*nyblock ;
                  
	PARCEL_INT jxblockm = (jxblock-1)>1 ? jxblock-1 : 1;
	PARCEL_INT jxblockp = (jxblock+1)<nxblock ? jxblock+1 : nxblock;
	PARCEL_INT jyblockm = (jyblock-1)>1 ? jyblock-1 : 1;
	PARCEL_INT jyblockp = (jyblock+1)<nyblock ? jyblock+1 : nyblock;
	PARCEL_INT jzblockm = (jzblock-1)>1 ? jzblock-1 : 1;
	PARCEL_INT jzblockp = (jzblock+1)<nzblock ? jzblock+1 : nzblock;

	for(PARCEL_INT jjzblock = jzblockm; jjzblock<=jzblockp;jjzblock++){
	  PARCEL_INT jjz0;
	  switch(jjzblock - jzblock){
	  case -1:
	    nzs = nzs_[jzblock-1];
	    nze = nzs_[jzblock-1] + addL - 1;
	    jjz0 = (nze_[jjzblock-1] - nzs_[jjzblock-1] + 1)  - addL ;
	    break;
	  case  0:
	    nzs = nzs_[jzblock-1];
	    nze = nze_[jzblock-1];
	    jjz0 = 0;
	    break;
	  case 1:
	    nzs = nze_[jzblock-1] - addL + 1;
	    nze = nze_[jzblock-1] ;
	    jjz0 = 0;
	    break;
	  }

	  for(PARCEL_INT jjyblock = jyblockm; jjyblock<=jyblockp;jjyblock++){
	    PARCEL_INT jjy0;
	    switch(jjyblock - jyblock){
	    case -1:
	      nys = nys_[jyblock-1];
	      nye = nys_[jyblock-1] + addL - 1;
	      jjy0 = (nye_[jjyblock-1] - nys_[jjyblock-1] + 1)  - addL ;
	      break;
	    case  0:
	      nys = nys_[jyblock-1];
	      nye = nye_[jyblock-1];
	      jjy0 = 0;
	      break;
	    case 1:
	      nys = nye_[jyblock-1] - addL + 1;
	      nye = nye_[jyblock-1] ;
	      jjy0 = 0;
	      break;
	    }

	    for(PARCEL_INT jjxblock = jxblockm; jjxblock<=jxblockp;jjxblock++){
	      PARCEL_INT jjx0;
	      switch(jjxblock - jxblock){
	      case -1:
		nxs = nxs_[jxblock-1];
		nxe = nxs_[jxblock-1] + addL - 1;
		jjx0 = (nxe_[jjxblock-1] - nxs_[jjxblock-1] + 1)  - addL ;
		break;
	      case  0:
		nxs = nxs_[jxblock-1];
		nxe = nxe_[jxblock-1];
		jjx0 = 0;
		break;
	      case 1:
		nxs = nxe_[jxblock-1] - addL + 1;
		nxe = nxe_[jxblock-1] ;
		jjx0 = 0;
		break;
	      }


	      PARCEL_INT jjblock_id = (jjxblock-1) + (jjyblock-1)*nxblock + (jjzblock-1)*nxblock*nyblock;

	      for(PARCEL_INT jz = nzs-1;jz<=nze-1;jz++){
		for(PARCEL_INT jy = nys-1;jy<=nye-1;jy++){
		  for(PARCEL_INT jx = nxs-1;jx<=nxe-1;jx++){
		    PARCEL_INT jjx = jx - (nxs_[jxblock-1]-1);
		    PARCEL_INT jjy = jy - (nys_[jyblock-1]-1);
		    PARCEL_INT jjz = jz - (nzs_[jzblock-1]-1);
		    PARCEL_INT jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1;

		    PARCEL_INT jjjx = jx - (nxs-1) + jjx0;
		    PARCEL_INT jjjy = jy - (nys-1) + jjy0;
		    PARCEL_INT jjjz = jz - (nzs-1) + jjz0;
		    PARCEL_INT jjj = (jjjx+stmx) + mx_local*(jjjy+stmy) + mx_local*my_local*(jjjz+stmz) + 1;
		    svec_AS[block_id + jj*nblock-1] = svec_AS[block_id   + jj*nblock-1] + svec[jjblock_id + jjj*nblock-1]; 
		  }
		}
	      }
	    }
	  }
	}

	nzs = nzs_[jzblock-1];
	nze = nze_[jzblock-1];
	nys = nys_[jyblock-1];
	nye = nye_[jyblock-1];
	nxs = nxs_[jxblock-1];
	nxe = nxe_[jxblock-1];

	PARCEL_INT nzs0 = nzs0_[jzblock-1];
	PARCEL_INT nys0 = nys0_[jyblock-1];
	PARCEL_INT nxs0 = nxs0_[jxblock-1];
	
	PARCEL_INT nze0 = nze0_[jzblock-1];
	PARCEL_INT nye0 = nye0_[jyblock-1];
	PARCEL_INT nxe0 = nxe0_[jxblock-1];

	PARCEL_INT addzs = nzs0 - nzs;
	PARCEL_INT addys = nys0 - nys;
	PARCEL_INT addxs = nxs0 - nxs;
	for(PARCEL_INT jz = nzs0-1; jz<=nze0-1;jz++){
	  for(PARCEL_INT jy = nys0-1; jy<=nye0-1;jy++){
	    for(PARCEL_INT jx = nxs0-1; jx<=nxe0-1;jx++){
	      PARCEL_INT j  = (jx+stmx) + mx*(jy+stmy) + mxy*(jz+stmz) + 1;
	      PARCEL_INT jjx = jx - (nxs0-1) + addxs;
	      PARCEL_INT jjy = jy - (nys0-1) + addys;
	      PARCEL_INT jjz = jz - (nzs0-1) + addzs;
	      PARCEL_INT jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1;
	      s[j-1] = svec_AS[block_id + jj*nblock-1];
	    }
	  }
	}

	if(jzblock==1){
	  nzs0  = nzs;
	  addzs = 0;
	}
	if(jzblock==nzblock){
	  nze0  = nze;
	}

	if(jyblock==1){
	  nys0  = nys;
	  addys = 0;
	}
	if(jyblock==nyblock){
	  nye0  = nye;
	}
	
	if(jxblock==1){
	  nxs0  = nxs;
	  addxs = 0;
	}
	if(jxblock==nxblock){
	  nxe0  = nxe;
	}
        // -------
	for(PARCEL_INT jz = nzs0-1; jz<=nze0-1;jz++){
	  for(PARCEL_INT jy = nys0-1; jy<=nye0-1;jy++){
	    for(PARCEL_INT jx = nxs0-1; jx<=nxe0-1;jx++){
	      PARCEL_INT j_AS = 
		(jx + sleeve_mx) + 
		(jy + sleeve_my) * mx_AS  + 
		(jz + sleeve_mz) * mxy_AS + 
		1;
	      PARCEL_INT jjx = jx - (nxs0-1) + addxs;
	      PARCEL_INT jjy = jy - (nys0-1) + addys;
	      PARCEL_INT jjz = jz - (nzs0-1) + addzs;
	      PARCEL_INT jj = (jjx+stmx) + mx_local*(jjy+stmy) + mx_local*my_local*(jjz+stmz) + 1;
	      AS_v2[j_AS-1] = svec_AS[block_id + jj*nblock-1];
	    }
	  }
	}

      }
    }
  }

  return ;

}

#endif

template <class Tdata,class Tcalc>
__global__
void AS_copy_(
	      Tdata *wv1,Tdata *AS_wv1,
	      int nx,int ny,int nz,
	      int mx,int mxy,
	      int stmx,int stmy,int stmz,
	      int mx_AS,int mxy_AS,
	      int sleeve_mx,int sleeve_my,int sleeve_mz
		){
  // for(int i=1; i<=AS_ndist;i++){
  //   AS_wv1[i-1] = 0.0;
  // }
#if 0
  for(PARCEL_INT jz = 0;jz<=nz-1;jz++){
    for(PARCEL_INT jy = 0;jy<=ny-1;jy++){
      for(PARCEL_INT jx = 0;jx<=nx-1;jx++){
#else
  for(
      PARCEL_INT jz = threadIdx.z + blockIdx.z*blockDim.z;
      jz<nz ;
      jz = jz + blockDim.z*gridDim.z ){
    for(
	PARCEL_INT jy = threadIdx.y + blockIdx.y*blockDim.y;
	jy<ny ;
	jy = jy + blockDim.y*gridDim.y ){      
      for(
	  PARCEL_INT jx = threadIdx.x + blockIdx.x;
	  jx<nx ;
	  jx = jx + blockDim.x*gridDim.x ){
#endif
	PARCEL_INT AS_i = ( jx +sleeve_mx ) + mx_AS *( jy + sleeve_my ) + mxy_AS *( jz + sleeve_mz ) + 1;
        PARCEL_INT i    = ( jx + stmx     ) + mx    *( jy + stmy      ) + mxy    *( jz + stmz )      + 1;
	AS_wv1[AS_i-1] = wv1[i-1];
          }
       }
    }
  return ;
}


template <class Tdata,class Tcalc>
__global__
void copy_AS_(
	      Tdata *wv1,Tdata *AS_wv1,
	      int nx,int ny,int nz,
	      int mx,int mxy,
	      int stmx,int stmy,int stmz,
	      int mx_AS,int mxy_AS,
	      int sleeve_mx,int sleeve_my,int sleeve_mz
		){
  // for(int i=1; i<=AS_ndist;i++){
  //   AS_wv1[i-1] = 0.0;
  // }
#if 0
  for(PARCEL_INT jz = 0;jz<=nz-1;jz++){
    for(PARCEL_INT jy = 0;jy<=ny-1;jy++){
      for(PARCEL_INT jx = 0;jx<=nx-1;jx++){
#else
  for(
      PARCEL_INT jz = threadIdx.z + blockIdx.z*blockDim.z;
      jz<nz ;
      jz = jz + blockDim.z*gridDim.z ){
    for(
	PARCEL_INT jy = threadIdx.y + blockIdx.y*blockDim.y;
	jy<ny ;
	jy = jy + blockDim.y*gridDim.y ){      
      for(
	  PARCEL_INT jx = threadIdx.x + blockIdx.x;
	  jx<nx ;
	  jx = jx + blockDim.x*gridDim.x ){
#endif
	PARCEL_INT AS_i = ( jx +sleeve_mx ) + mx_AS *( jy + sleeve_my ) + mxy_AS *( jz + sleeve_mz ) + 1;
        PARCEL_INT i    = ( jx + stmx     ) + mx    *( jy + stmy      ) + mxy    *( jz + stmz )      + 1;
	wv1[i-1] = AS_wv1[AS_i-1];
          }
       }
    }
  return ;
}

template <class Tdata,class Tcalc>
__global__
void AS_set_sbuf_(
		  Tdata *v,Tdata *AS_wv1,
		  int nloop,
		  int mx_AS,int mxy_AS,
		  int xlen_s,int ylen_s,
		  int x0s,int y0s,int z0s		  
		  ){
  const PARCEL_INT tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const PARCEL_INT step = blockDim.x*gridDim.x;
  //  for(PARCEL_INT j=0;j<=nloop-1;j++){
  for(PARCEL_INT j=tid;j<=nloop-1;j=j+step){
    PARCEL_INT ix = x0s +  j % xlen_s;
    PARCEL_INT iy = y0s + (j / xlen_s) % ylen_s;
    PARCEL_INT iz = z0s +  j / (xlen_s * ylen_s);
    PARCEL_INT k  = ix + mx_AS * iy + mxy_AS * iz + 1;
    v[j+1-1] = AS_wv1[k-1];
  }

  return ;
}

template <class Tdata,class Tcalc>
__global__
void AS_set_rbuf_(
		  Tdata *v,Tdata *AS_wv1,
		  int nloop,
		  int mx_AS,int mxy_AS,
		  int xlen_r,int ylen_r,
		  int x0r,int y0r,int z0r,
		  Tdata alpha
		  ){

  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
  //  for(PARCEL_INT j=0;j<=nloop-1;j++){
  for(PARCEL_INT j=tid;j<=nloop-1;j=j+step){
    PARCEL_INT ix = x0r +  j%xlen_r;
    PARCEL_INT iy = y0r + (j/xlen_r)%ylen_r ;
    PARCEL_INT iz = z0r +  j / (xlen_r * ylen_r);
    PARCEL_INT k = ix + mx_AS * iy + mxy_AS * iz + 1;
    //    AS_wv1[k-1] = v[j+1-1];
    AS_wv1[k-1] = alpha*AS_wv1[k-1] + v[j+1-1] ;
  }

  return ;
}

extern "C" {
  void cuda_solve_DILU_DDM(
			    int numGPUblock,
			    int numGPUthread,
			    devicemem *A,
			    devicemem *Dinv_, PARCEL_INT offset_Dinv,
			    devicemem *r_,    PARCEL_INT offset_r,
			    devicemem *s_,    PARCEL_INT offset_s,
			    devicemem *svec_, PARCEL_INT offset_svec,
			    devicemem *ioff_dia,
			    devicemem *offset_xyz,
			    devicemem *nxe,
			    devicemem *nye,
			    devicemem *nze,
			    int nx,
			    int ny,
			    int nz,
			    int stmx,
			    int stmy,
			    int stmz,
			    int mx,
			    int my,
			    int mz,
			    int nxblock,
			    int nyblock,
			    int nzblock,
			    int nx_local,
			    int ny_local,
			    int nz_local,
			    int mx_local,
			    int my_local,
			    int mz_local,
			    int nnd){

    double *Dinv = (double *)Dinv_[0].a_d;
    double *r    = (double *)r_[0].a_d;
    double *s    = (double *)s_[0].a_d;
    double *svec = (double *)svec_[0].a_d;

    int nx_para = nx < 1024 ? nx : 512;
    int ny_para = ny < 1024 ? ny : 512;
    int nz_para = nz < 1024 ? nz : 512;

    if(nxblock<1024){
#if 1
      dim3 block(nx_para,1,1); 
      dim3 grid(1,ny_para,nz_para); 
#else
      dim3 block(1,1,1); 
      dim3 grid(1,1,1); 
#endif
      cudaDeviceSynchronize();
      solve_DILU_DDM_<double,double> 
	<<<grid,block>>>(
			    (double *)A[0].a_d,
			    &Dinv[offset_Dinv],
			    &r[offset_r],
			    &s[offset_s],
			    &svec[offset_svec],
			    (int *)ioff_dia[0].a_d,
			    (int *)offset_xyz[0].a_d,
			    (int *)nxe[0].a_d,
			    (int *)nye[0].a_d,
			    (int *)nze[0].a_d,
			    nx,
			    ny,
			    nz,
			    stmx,
			    stmy,
			    stmz,
			    mx,
			    my,
			    mz,
			    nxblock,
			    nyblock,
			    nzblock,
			    nx_local,
			    ny_local,
			    nz_local,
			    mx_local,
			    my_local,
			    mz_local,
			    nnd
			 );
    }else{
      dim3 block(nx_para,1,1); 
      dim3 grid(1,ny_para,nz_para); 
      cudaDeviceSynchronize();
      solve_DILU_DDM_<double,double> 
	<<<grid,block>>>(
			    (double *)A[0].a_d,
			    &Dinv[offset_Dinv],
			    &r[offset_r],
			    &s[offset_s],
			    &svec[offset_svec],
			    (int *)ioff_dia[0].a_d,
			    (int *)offset_xyz[0].a_d,
			    (int *)nxe[0].a_d,
			    (int *)nye[0].a_d,
			    (int *)nze[0].a_d,
			    nx,
			    ny,
			    nz,
			    stmx,
			    stmy,
			    stmz,
			    mx,
			    my,
			    mz,
			    nxblock,
			    nyblock,
			    nzblock,
			    nx_local,
			    ny_local,
			    nz_local,
			    mx_local,
			    my_local,
			    mz_local,
			    nnd
			 );
        }

    return ;
  }


  void cuda_solve_DILU_AS_DDM(
			      devicemem *AS_v1_  , PARCEL_INT offset_AS_v1,
			      devicemem *s_      , PARCEL_INT offset_s,
			      devicemem *rvec_   , PARCEL_INT offset_rvec,
			      devicemem *svec_   , PARCEL_INT offset_svec,
			      devicemem *svec_AS_, PARCEL_INT offset_svec_AS,
			      devicemem *AS_v2_  , PARCEL_INT offset_AS_v2,
			      devicemem *AS_val_dia_init_, 
			      devicemem *Dinv_    , PARCEL_INT offset_Dinv,
			      devicemem *ioff_dia_,
			      devicemem *offset_xyz_,
			      int nnd, 
			      devicemem *nxs_, devicemem *nxe_,
			      devicemem *nys_, devicemem *nye_,
			      devicemem *nzs_, devicemem *nze_,
			      devicemem *nxs0_, devicemem *nxe0_,
			      devicemem *nys0_, devicemem *nye0_,
			      devicemem *nzs0_, devicemem *nze0_,
			      int addL,
			      int nblock,
			      int nxblock,int nyblock,int nzblock,
			      int stmx,int stmy,int stmz,
			      int mx, int mxy,
			      int m_local,
			      int mx_local,int my_local,
			      int m_AS,
			      int mx_AS ,int mxy_AS,
			      int sleeve_mx,int sleeve_my,int sleeve_mz
			      ){

    double *AS_v1   = (double *)AS_v1_[offset_AS_v1].a_d;
    double *s       = (double *)s_[offset_s].a_d;
    double *rvec    = (double *)rvec_[offset_rvec].a_d;
    double *svec    = (double *)svec_[offset_svec].a_d;
    double *svec_AS = (double *)svec_AS_[offset_svec_AS].a_d;
    double *AS_v2   = (double *)AS_v2_[offset_AS_v2].a_d;
    double *AS_val_dia_init = (double *)AS_val_dia_init_[0].a_d;    
    double *Dinv = (double *)Dinv_[offset_Dinv].a_d; 

    int *ioff_dia = (int *)ioff_dia_[0].a_d;
    int *offset_xyz = (int *)offset_xyz_[0].a_d;

    int *nxs = (int *)nxs_[0].a_d;
    int *nxe = (int *)nxe_[0].a_d;
    int *nys = (int *)nys_[0].a_d;
    int *nye = (int *)nye_[0].a_d;
    int *nzs = (int *)nzs_[0].a_d;
    int *nze = (int *)nze_[0].a_d;

    int *nxs0 = (int *)nxs0_[0].a_d;
    int *nxe0 = (int *)nxe0_[0].a_d;
    int *nys0 = (int *)nys0_[0].a_d;
    int *nye0 = (int *)nye0_[0].a_d;
    int *nzs0 = (int *)nzs0_[0].a_d;
    int *nze0 = (int *)nze0_[0].a_d;

    cudaDeviceSynchronize();
    cuda_axpy_<double,double> <<<2048,512>>>(
					     m_local * nblock,
					     0.0,
					     0.0,
					     svec,
					     svec
					     );
    cuda_axpy_<double,double> <<<2048,512>>>(
					     m_local * nblock,
					     0.0,
					     0.0,
					     rvec,
					     rvec
					     );
    cuda_axpy_<double,double> <<<2048,512>>>(
					     m_local * nblock,
					     0.0,
					     0.0,
					     svec_AS,
					     svec_AS
					     );

    int nx_para = nxblock < 1024 ? nxblock : 512;
    int ny_para = nyblock < 1024 ? nyblock : 512;
    int nz_para = nzblock < 1024 ? nzblock : 512;

    if(nxblock<1024){
#if loop_single
      dim3 block(1,1,1); 
      dim3 grid(1,1,1); 
#else
      dim3 block(nx_para,1,1); 
      dim3 grid(1,ny_para,nz_para); 
#endif
      cudaDeviceSynchronize();
      solve_DILU_AS_DDM_setvec_<double,double>
	<<<grid,block>>>(
			 AS_v1,s,
			 rvec,svec,svec_AS,
			 AS_v2,
			 AS_val_dia_init, // 行列データ
			 Dinv,            // 対角要素の逆行列
			 ioff_dia,offset_xyz,
			 nnd, 
			 nxs, nxe,
			 nys, nye,
			 nzs, nze,
			 nxs0, nxe0,
			 nys0, nye0,
			 nzs0, nze0,
			 addL,
			 nblock,
			 nxblock,nyblock,nzblock,
			 stmx,stmy,stmz,
			 mx, mxy,
			 mx_local,my_local,
			 m_AS,
			 mx_AS ,mxy_AS,
			 sleeve_mx,sleeve_my,sleeve_mz
			 );
      cudaDeviceSynchronize();
      solve_DILU_AS_DDM_solve_<double,double>
	<<<grid,block>>>(
			 AS_v1,s,
			 rvec,svec,svec_AS,
			 AS_v2,
			 AS_val_dia_init, // 行列データ
			 Dinv,            // 対角要素の逆行列
			 ioff_dia,offset_xyz,
			 nnd, 
			 nxs, nxe,
			 nys, nye,
			 nzs, nze,
			 nxs0, nxe0,
			 nys0, nye0,
			 nzs0, nze0,
			 addL,
			 nblock,
			 nxblock,nyblock,nzblock,
			 stmx,stmy,stmz,
			 mx, mxy,
			 mx_local,my_local,
			 m_AS,
			 mx_AS ,mxy_AS,
			 sleeve_mx,sleeve_my,sleeve_mz
			 );
      cudaDeviceSynchronize();
      solve_DILU_AS_DDM_update_<double,double>
	<<<grid,block>>>(
			 AS_v1,s,
			 rvec,svec,svec_AS,
			 AS_v2,
			 AS_val_dia_init, // 行列データ
			 Dinv,            // 対角要素の逆行列
			 ioff_dia,offset_xyz,
			 nnd, 
			 nxs, nxe,
			 nys, nye,
			 nzs, nze,
			 nxs0, nxe0,
			 nys0, nye0,
			 nzs0, nze0,
			 addL,
			 nblock,
			 nxblock,nyblock,nzblock,
			 stmx,stmy,stmz,
			 mx, mxy,
			 mx_local,my_local,
			 m_AS,
			 mx_AS ,mxy_AS,
			 sleeve_mx,sleeve_my,sleeve_mz
			 );
      cudaDeviceSynchronize();
    }else{
#if loop_single
      dim3 block(1,1,1); 
      dim3 grid(1,1,1); 
#else
      dim3 block(nx_para,1,1); 
      dim3 grid(1,ny_para,nz_para); 
#endif
      cudaDeviceSynchronize();
      solve_DILU_AS_DDM_setvec_<double,double>
	<<<grid,block>>>(
			 AS_v1,s,
			 rvec,svec,svec_AS,
			 AS_v2,
			 AS_val_dia_init, // 行列データ
			 Dinv,            // 対角要素の逆行列
			 ioff_dia,offset_xyz,
			 nnd, 
			 nxs, nxe,
			 nys, nye,
			 nzs, nze,
			 nxs0, nxe0,
			 nys0, nye0,
			 nzs0, nze0,
			 addL,
			 nblock,
			 nxblock,nyblock,nzblock,
			 stmx,stmy,stmz,
			 mx, mxy,
			 mx_local,my_local,
			 m_AS,
			 mx_AS ,mxy_AS,
			 sleeve_mx,sleeve_my,sleeve_mz
			 );
      cudaDeviceSynchronize();
      solve_DILU_AS_DDM_solve_<double,double>
	<<<grid,block>>>(
			 AS_v1,s,
			 rvec,svec,svec_AS,
			 AS_v2,
			 AS_val_dia_init, // 行列データ
			 Dinv,            // 対角要素の逆行列
			 ioff_dia,offset_xyz,
			 nnd, 
			 nxs, nxe,
			 nys, nye,
			 nzs, nze,
			 nxs0, nxe0,
			 nys0, nye0,
			 nzs0, nze0,
			 addL,
			 nblock,
			 nxblock,nyblock,nzblock,
			 stmx,stmy,stmz,
			 mx, mxy,
			 mx_local,my_local,
			 m_AS,
			 mx_AS ,mxy_AS,
			 sleeve_mx,sleeve_my,sleeve_mz
			 );
      cudaDeviceSynchronize();
      solve_DILU_AS_DDM_update_<double,double>
	<<<grid,block>>>(
			 AS_v1,s,
			 rvec,svec,svec_AS,
			 AS_v2,
			 AS_val_dia_init, // 行列データ
			 Dinv,            // 対角要素の逆行列
			 ioff_dia,offset_xyz,
			 nnd, 
			 nxs, nxe,
			 nys, nye,
			 nzs, nze,
			 nxs0, nxe0,
			 nys0, nye0,
			 nzs0, nze0,
			 addL,
			 nblock,
			 nxblock,nyblock,nzblock,
			 stmx,stmy,stmz,
			 mx, mxy,
			 mx_local,my_local,
			 m_AS,
			 mx_AS ,mxy_AS,
			 sleeve_mx,sleeve_my,sleeve_mz
			 );
      cudaDeviceSynchronize();
    }

    return ;
  }

  void cuda_AS_util_copy_only(
			 devicemem *wv1_, PARCEL_INT offset_wv1,
			 devicemem *AS_wv1_, PARCEL_INT offset_AS_wv1,
			 int nx,int ny,int nz,
			 int mx,int mxy,
			 int stmx,int stmy,int stmz,
			 PARCEL_INT AS_ndist,
			 int mx_AS,int mxy_AS,
			 int sleeve_mx,int sleeve_my,int sleeve_mz
			 ){

    int nx_para = nx < 1024 ? nx : 512;
    int ny_para = ny < 1024 ? ny : 512;
    int nz_para = nz < 1024 ? nz : 512;

    double *AS_wv1 = (double *)AS_wv1_[0].a_d;
    double *wv1 = (double *)wv1_[0].a_d;
    cudaDeviceSynchronize();
    cuda_axpy_<double,double> <<<2048,512>>>(
					     AS_ndist,
					     0.0,
					     0.0,
					     &AS_wv1[offset_AS_wv1],
					     &AS_wv1[offset_AS_wv1]
					     );
    // if(nx<1024){
    //   dim3 block(nx,1,1); 
    //   dim3 grid(1,ny,nz); 
    // }else{
    // }

    dim3 block(nx_para,1,1); 
    dim3 grid(1,ny_para,nz_para); 
    AS_copy_<double,double> 
      <<<grid,block>>>(
      //      <<<1,1>>>(
		       &wv1[offset_wv1],
		       &AS_wv1[offset_AS_wv1],
		       nx, ny, nz,
		       mx, mxy,
		       stmx, stmy, stmz,
		       mx_AS, mxy_AS,
		       sleeve_mx, sleeve_my, sleeve_mz
		       );

    return ;
  }

  void cuda_util_AS_copy_only(
			 devicemem *wv1_, PARCEL_INT offset_wv1,
			 devicemem *AS_wv1_, PARCEL_INT offset_AS_wv1,
			 int nx,int ny,int nz,
			 int mx,int mxy,
			 int stmx,int stmy,int stmz,
			 PARCEL_INT ndist,
			 int mx_AS,int mxy_AS,
			 int sleeve_mx,int sleeve_my,int sleeve_mz
			 ){

    int nx_para = nx < 1024 ? nx : 512;
    int ny_para = ny < 1024 ? ny : 512;
    int nz_para = nz < 1024 ? nz : 512;

    double *AS_wv1 = (double *)AS_wv1_[0].a_d;
    double *wv1 = (double *)wv1_[0].a_d;
    cudaDeviceSynchronize();
    cuda_axpy_<double,double> <<<2048,512>>>(
					     ndist,
					     0.0,
					     0.0,
					     &wv1[offset_wv1],
					     &wv1[offset_wv1]
					     );
    // if(nx<1024){
    //   dim3 block(nx,1,1); 
    //   dim3 grid(1,ny,nz); 
    // }else{
    // }

    dim3 block(nx_para,1,1); 
    dim3 grid(1,ny_para,nz_para); 
    copy_AS_<double,double> 
      <<<grid,block>>>(
		       &wv1[offset_wv1],
		       &AS_wv1[offset_AS_wv1],
		       nx, ny, nz,
		       mx, mxy,
		       stmx, stmy, stmz,
		       mx_AS, mxy_AS,
		       sleeve_mx, sleeve_my, sleeve_mz
		       );

    return ;
  }


  void cuda_AS_set_sbuf(
			devicemem *v_, PARCEL_INT offset_v,
			devicemem *AS_wv1_, PARCEL_INT offset_AS_wv1,
			int nloop,
			int mx_AS,int mxy_AS,
			int xlen_s,int ylen_s,
			int x0s,int y0s,int z0s 
			){
    double *v      = (double *)v_[0].a_d;
    double *AS_wv1 = (double *)AS_wv1_[0].a_d;
    AS_set_sbuf_<double,double> 
      <<<512,512>>>(
		    &v[offset_v],
		    &AS_wv1[offset_AS_wv1],
		    nloop,
		    mx_AS,mxy_AS,
		    xlen_s,ylen_s,
		    x0s,y0s,z0s		 
		    );
    return ;
  }

  void cuda_AS_set_rbuf(
			devicemem *v_, PARCEL_INT offset_v,
			devicemem *AS_wv1_, PARCEL_INT offset_AS_wv1,
			int nloop,
			int mx_AS,int mxy_AS,
			int xlen_r,int ylen_r,
			int x0r,int y0r,int z0r ,
			double alpha
			){
    double *v      = (double *)v_[0].a_d;
    double *AS_wv1 = (double *)AS_wv1_[0].a_d;
    AS_set_rbuf_<double,double> 
      <<<512,512>>>(
		    &v[offset_v],
		    &AS_wv1[offset_AS_wv1],
		    nloop,
		    mx_AS,mxy_AS,
		    xlen_r,ylen_r,
		    x0r,y0r,z0r,
		    alpha
		    );
    return ;
  }

}

