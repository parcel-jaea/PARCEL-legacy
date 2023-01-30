#ifndef HEADER_CUDAERRCHECK
#define HEADER_CUDAERRCHECK
#include <stdio.h>
#include <stdlib.h>

#define cudaErrchk(ans) { cudaAssert((ans), __FILE__, __LINE__); }
inline void cudaAssert(cudaError_t code, const char *file, int line, bool abort=true){
  if (code != cudaSuccess) {
    fprintf(stderr,"CUDA assert: %s %s %d\n", cudaGetErrorString(code), file, line);
    if (abort) exit(code);
  }
}

/* #if 1 */
/* #include <helper_functions.h>  // shared functions common to CUDA Samples */
/* #include <helper_cuda.h>       // CUDA error checking */
/* #else */
/* #define checkCudaErrors(val){ inline void dummy_func();} */
/* #endif */

#endif


