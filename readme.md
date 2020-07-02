# Directories
MPI    :: Quad precision MPI library
SPARSE :: PARCEL

# Compile
Compile command    :: make
Configuration file :: ./arch/make_config
  CC         = C compiler
  FC         = FORTRAN compiler
  CFLAGS     = C compiler options
  FFLAGS     = FORTRAN compiler options
  FP_MODE    = Option for keeping computation order in optimization of quad precision routines
  INCLUDEDIR = Path of include files
  LD_MPI     = Path of MPI library
  LD_LAPACK  = Path of LAPACK library
