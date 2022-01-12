#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"


#ifdef DPARCEL_INT8
#define parcel_int long long 
#else
#define parcel_int int
#endif

int main( int argc, char *argv[]){
  const parcel_int n = 10000;
  const parcel_int s = 10;
  MPI_Fint icomm_fort;
  int npes,myrank;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  icomm_fort=MPI_Comm_c2f(MPI_COMM_WORLD);

  double *A  = (double *)malloc(sizeof(double)*n*s);
  double *A0 = (double *)malloc(sizeof(double)*n*s);
  double *R  = (double *)malloc(sizeof(double)*s*s);
  double *Q  = (double *)malloc(sizeof(double)*n*s);
  int i,j;

  for(j=0;j<myrank+1;j++){
    for(i=0;i<n*s;i++){
      A[i] = (double)rand()/RAND_MAX;
    }
  }

  for(i=0;i<s*s;i++){
    R[i]=0.0;
  }

  for(i=0;i<n*s;i++){
    A0[i] = A[i];
    Q[i]  = 0.0;
  }

  for(i=1;i<6;i++){
    if(myrank == 0){
      printf(" i = %d \n",i);
      switch(i){
      case 1:
	printf("QR decomposition : MGS \n");
	break;
      case 2:
	printf("QR decomposition : CGS \n");
	break;
      case 3:
	printf("QR decomposition : TSQR \n");
	break;
      case 4:
	printf("QR decomposition : CHOLESKYQR  \n");
	break;
      case 5:
	printf("QR decomposition : CHOLESKYQR2 \n");
	break;
      }
    }
    for(j=0;j<n*s;j++){
      A[j] = A0[j];
      Q[j] = 1.0;
    }
    parcel_dqr_(&n,&s,&icomm_fort,A,Q,R,&i);
    check_qr_c_(&n,&s,&icomm_fort,A,Q,R,A0);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  free(A);
  free(A0);
  free(R);
  free(Q);

  return 0;
}
