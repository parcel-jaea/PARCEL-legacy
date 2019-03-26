confdir=../../
include ${confdir}/arch/make_config

INC = -I../src/ # -Dcalc_diaMat
LDFLAGS = ${LD_LAPACK} ${LD_MPI} ../src/parcel_sparse.a  ../src/dd_MPI/lib/libMpifdd1.a 

PROG1    = run
PROBLEM_FILE=problem.o 
OBJS1 = ${PROBLEM_FILE} main.o

%.o : %.c
	${CC} -c ${CFLAGS} $<
%.o : %.f90
	${FC} -c ${FFLAGS} ${INCLUDEDIR} ${INC} $<
all : ${OBJS1} 
	${FC} -o ${PROG1} ${OBJS1} ${FFLAGS} ${INCLUDEDIR} ${LDFLAGS} 
clean :
	rm ${PROG1} ${PROG2} *o *mod *~
