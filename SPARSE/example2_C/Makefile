# confdir=../../
# include ${confdir}/arch/make_config

CC     = mpiicc 
FC     = mpiifort 

CFLAGS  = -O3 
FFLAGS  = -O3 -fpp -qopenmp -xHost  -mcmodel=large 
LD_LAPACK  = -mkl 

####

INC = -I../src/
LDFLAGS = ${LD_LAPACK} ${LD_MPI} ../src/parcel_sparse.a  \
	/home/app/intel/compilers_and_libraries_2018.3.222/linux/compiler/lib/intel64/libifcore.a 
PROG1    = run
PROBLEM_FILE=problem.o 
OBJS1 = ${PROBLEM_FILE} C_wrapper.o main.o

%.o : %.c
	${CC} -c ${CFLAGS} $<
%.o : %.f90
	${FC} -c ${FFLAGS} ${INCLUDEDIR} ${INC} $<
all : ${OBJS1} 
	${CC} -o ${PROG1} ${CFLAGS} ${OBJS1} ${LDFLAGS} 
clean :
	rm ${PROG1} ${PROG2} *o *mod *~

