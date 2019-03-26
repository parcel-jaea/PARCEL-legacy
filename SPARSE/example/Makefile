confdir=../../
include ${confdir}/arch/make_config
LDFLAGS = ${LD_LAPACK} ${LD_MPI} ../src/parcel_sparse.a

PROG    = sample1
OBJS    = sample1.o matrix_poisson_2D.o  matrix_diffusion.o matrix_random.o matrix_wilkinson.o matrix_frank.o genmxran.o  matrix_poisson_3D.o \
	im_conv_from_coo.o cputime.o wallclock.o

%.o : %.c
	${CC} -c -DPOSIX_C ${CFLAGS} $<
%.o : %.f90
	${FC} -c ${FFLAGS} ${INCLUDEDIR} $<
%.o : %.f
	${FC} -c ${FFLAGS} ${INCLUDEDIR} $<
%.o : %.F
	${FC} -c ${FFLAGS} ${INCLUDEDIR} $<
all : ${PROG}
${PROG} : ${OBJS}
	${FC} -o $@ ${OBJS} ${LDFLAGS} ${FFLAGS} ${INCLUDEDIR}
clean :
	rm ${PROG} *o *mod *~
