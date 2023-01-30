confdir=../../../
include ${confdir}/arch/make_config

OBJS = no_cuda.o
ifeq (${CUDAFLAG},ON)
OBJS = \
	cu_common.o cu_comm.o cu_SpMV.o \
	cu_qrdecomposition.o \
	cu_neumann.o \
	cu_DILU.o \
	cusparse.o cu_ddm_precon.o\
	cu_pcg.o cu_pbicg.o cu_pcagmres.o cu_pcbcg.o \
	cu_ddm.o cublasBatched.o cu_tsqr.o 
endif

CU_INCLUDEDIR = -I./
%.o : %.cpp
	${CXX} -c ${CXXFLAGS} ${FDEFINE} $<
%.o : %.c
	${CC} -c ${CFLAGS} ${FDEFINE} $<
%.o : %.f90
	${FC} -c ${FFLAGS} ${INCLUDEDIR} ${FDEFINE} $<
%.o : %.cu
	${NVC} -c ${NVCFLAGS} ${FDEFINE}  ${CU_INCLUDEDIR} $<

# all : ${PROG}
all :  ${OBJS}
# ${PROG} : ${OBJS}
# 	ar rc  $@ ${OBJS}
clean :
	rm  *.o *~
# rm  ${OBJS} no_cuda.o *~


