#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int compare_int(const void *a, const void *b){
  return *(int*)a - *(int*)b;
}

int sort_c_(int *A,int *n_ ){
  size_t n=n_[0];
  qsort(A,n,sizeof(int), compare_int);
  return 0;
}

#ifdef PARCEL_INT8
#define PARCEL_INT long long 
#else
#define PARCEL_INT int
#endif
PARCEL_INT binary_search_main(PARCEL_INT ary[], PARCEL_INT key, PARCEL_INT imin, PARCEL_INT imax) {
    if (imax < imin) {
      return -99; // keyが存在しない
    } else {
        PARCEL_INT imid = imin + (imax - imin) / 2;
        if (ary[imid] > key) {
            return binary_search_main(ary, key, imin, imid - 1);
        } else if (ary[imid] < key) {
            return binary_search_main(ary, key, imid + 1, imax);
        } else {
            return imid;
        }
    }
}

PARCEL_INT binary_search_(PARCEL_INT *ary, PARCEL_INT *key, PARCEL_INT *imin, PARCEL_INT *imax,PARCEL_INT *id) {
  PARCEL_INT key0  = key[0];
  PARCEL_INT imin0 = imin[0]-1;
  PARCEL_INT imax0 = imax[0]-1;
  
  id[0] = binary_search_main(ary, key0, imin0,  imax0) ;
  if(id[0] > -1){
    id[0] = id[0] +1;
  }else{
    id[0] = -1;
  }
  return 0;
}

