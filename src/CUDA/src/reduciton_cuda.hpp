
#ifndef HEADER_REDUCTIOCUDA
#define HEADER_REDUCTIOCUDA

// __device__ static type warpReduceSum(type val);
// __device__ static void local_sum_cuda(type val,type *val_ret,int thid);

template <class Tdata,class Tcalc>
__device__ 
static Tdata warpReduceSum(Tdata val){
  for(int offset = warpSize/2;offset >0; offset/=2){
    //    val+= __shfl_down_sync(0xffffffff,val,offset);
    val = val + __shfl_down_sync(0xffffffff,val,offset);
  }
  return val;
}

template <class Tdata,class Tcalc>
__device__ static void local_sum_cuda(Tdata val,Tdata *val_ret,int thid){
  Tdata warp_val = warpReduceSum<Tdata,Tcalc>(val);
  int lane = thid % warpSize;
  if(lane==0)atomicAdd(&val_ret[0], warp_val);
  return ;
}
template <class T>
__device__ __forceinline__ void
bcast_warp(T &x_, const int root){
  T x = x_;
  x = __shfl_sync( 0xffffffff, x, root, 32 );
  x_ = x;
}

#endif
