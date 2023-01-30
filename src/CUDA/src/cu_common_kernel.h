#ifndef HEADER_CUCOMMONKERNEL
#define HEADER_CUCOMMONKERNEL
extern "C"{
  void cuda_sqrt_val(
		     int numGPUblock,
		     int numGPUthread,
		     devicemem *x);

  void cuda_add_val(
		     int numGPUblock,
		     int numGPUthread,
		     double val,
		     devicemem *x);
  void cuda_sub_val(
		     int numGPUblock,
		     int numGPUthread,
		     double val,
		     devicemem *x);
  void cuda_mul_val(
		     int numGPUblock,
		     int numGPUthread,
		     double val,
		     devicemem *x);
  void cuda_div_val(
		     int numGPUblock,
		     int numGPUthread,
		     double val,
		     devicemem *x);
  
}

#endif
