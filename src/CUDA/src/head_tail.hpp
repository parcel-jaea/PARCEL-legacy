#include "cu_common.h"

#ifndef HEADER_HEADTAIL
#define HEADER_HEADTAIL

__device__ __forceinline__
static void set_Head_Tail(
			   PARCEL_INT downmax,
			   PARCEL_INT upmax,
			   PARCEL_INT mydown,
			   PARCEL_INT myup,
			   PARCEL_INT *myHead_,
			   PARCEL_INT *myTail_
			   ){
#if 1
  PARCEL_INT myHead = mydown > downmax ? mydown: downmax;
  PARCEL_INT myTail = myup<upmax ? myup : upmax;

  if(mydown>upmax){
    myHead =  999;
    myTail = -999; 
  }
  if(myup<downmax){
    myHead =  999;
    myTail = -999;
  }

#else
  PARCEL_INT myHead = mydown>downmax ? mydown: downmax;
  PARCEL_INT myTail = myup<upmax ? myup : upmax;

  if(mydown>upmax){
    myHead =  999;
    myTail = -999; 
  }
  if(myup<downmax){
    myHead =  999;
    myTail = -999;
  }

#endif

  myHead_[0] = myHead;
  myTail_[0] = myTail;

}

#endif
