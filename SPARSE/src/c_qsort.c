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
