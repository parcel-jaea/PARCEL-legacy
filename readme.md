# Compile
Compile command    :: make ; make install
Configuration file :: ./arch/make_config
  CC         = C compiler
  FC         = FORTRAN compiler
  CFLAGS     = C compiler options
  FFLAGS     = FORTRAN compiler options
  FP_MODE    = Option for keeping computation order in optimization of quad precision routines
  CUDAFLAG   = CUDA compile flag (ON:enable.OFF:disable)
  CUDAAWARE  = CUDA-Aware MPI compile flag (ON:enable.OFF:disable)
  INCLUDEDIR = Path of include files
  LD_MPI     = Path of MPI library
  LD_LAPACK  = Path of LAPACK library
  LD_CUDA    = CUDA library(cuBLAS,cuSPARSE) 
  FDEFINE    = 64bit integer becomes available with -DPARCEL_INT8
  INSTALLDIR = Path of PARCEL install directory

# Install
 1. Edit ./arch/make_config depending on compilation environments on your computing platforms.
 2. Execute "make; make install"
 3. The following directories are generated in INSTALLDIR, which is specified in make_config.
 -- INSTALLDIR/include contains module files and header files.
 -- INSTALLDIR/lib contains parcel_sparse.a, which is liked with your program.
 -- INSTALLDIR/example contains sample codes, which are compiled by executing "make" in this directory.

