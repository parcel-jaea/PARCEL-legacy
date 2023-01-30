#ifndef HEADER_CUCOMM
#define HEADER_CUCOMM

template <class Tdata>
__global__
void get_nrk_vecLU_( 
		    Tdata *rbuff,
		    PARCEL_INT nColL,
		    PARCEL_INT nColU,
		    Tdata *rbuffL,
		    Tdata *rbuffU
		     ){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;
#if 0
  for(PARCEL_INT i=tid;i<nColL;i=i+step){
    rbuffL[i] = 0.0;
  }
  for(PARCEL_INT i=tid;i<nColU;i=i+step){
    rbuffU[i] = 0.0;
  }
#else
  for(PARCEL_INT i=tid;i<nColL;i=i+step){
    rbuffL[i] = rbuff[i];
    //    prPARCEL_INTf(" rbuffL i =%06d %e \n",i,rbuffL[i]);
  }
  for(PARCEL_INT i=tid;i<nColU;i=i+step){
    rbuffU[i] = rbuff[i+nColL];
    //    printf(" rbuffU i =%06d %e \n",i,rbuffU[i]);
  }
#endif
  return ;
}

template <class Tdata>
__global__
void nrk_vec_pack_(
		   PARCEL_INT i0, 
		   PARCEL_INT ndata, 
		   Tdata *vec,
		   Tdata *sbuf,
		   PARCEL_INT *offset){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;

  for(PARCEL_INT i=tid; i<ndata; i=i+step){
    sbuf[i+i0]=vec[offset[i]];
  }

  return ;
}
#endif
