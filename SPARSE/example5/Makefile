confdir=../../
include ${confdir}/arch/make_config

INC = -I../src/
LDFLAGS = ${LD_LAPACK} ${LD_MPI} ../src/parcel_sparse.a  ../src/dd_MPI/lib/libMpifdd1.a \
	<path-to>/libjaeri_ilesf77.a

# PROG1    = run_int8
# PROG2    = run_int4
PROG1    = run
PROBLEM_FILE=problem.o problem_int8.o
OBJS1 = ${PROBLEM_FILE} main.o
# OBJS1 = ${PROBLEM_FILE} main_int4.o

%.o : %.c
	${CC} -c ${CFLAGS} $<
%.o : %.f90
	${FC} -c ${FFLAGS} ${INCLUDEDIR} ${INC} $<

all : ${OBJS1} 
# all : ${OBJS1} ${OBJS2}
	${FC} -o ${PROG1} ${OBJS1} ${FFLAGS} ${INCLUDEDIR} ${LDFLAGS} 
#	${FC} -o ${PROG2} ${OBJS2} ${FFLAGS} ${INCLUDEDIR} ${LDFLAGS} 

clean :
	rm ${PROG1} ${PROG2} *o *mod *~
