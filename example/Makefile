
confdir=./
include ${confdir}/make_config
PARCEL_dir = ../
INCLUDEDIR = -I${PARCEL_dir}/include/
LDFLAGS = ${PARCEL_dir}/lib/parcel_sparse.a  ${LD_LAPACK} ${LD_MPI} 

# INCLUDEDIR = -I${PARCEL_dir}/src/
# LDFLAGS = ${LD_LAPACK} ${LD_MPI} ${PARCEL_dir}/src/parcel_sparse.a 

PROG1    = run
PROG2    = run_dd
PROG3    = run_ddm
PROG4    = run_ddm_dd
PROG5    = run_qr

PROG6    = run_c
PROG7    = run_dd_c
PROG8   = run_ddm_c
PROG9   = run_ddm_dd_c
PROG10    = run_qr_c

PROBLEM_FILE=problem.o 
OBJS1 = ${PROBLEM_FILE} main_direct.o
OBJS2 = ${PROBLEM_FILE} main_dd.o
OBJS3 = ${PROBLEM_FILE} ddm_setup.o main_ddm.o
OBJS4 = ${PROBLEM_FILE} ddm_setup.o main_ddm_dd.o
OBJS5 = main_QR.o

OBJS6 = ${PROBLEM_FILE} ddm_setup.o c2fortran_kernel.o main_direct_c.o
OBJS7 = ${PROBLEM_FILE} ddm_setup.o c2fortran_kernel.o main_dd_c.o
OBJS8 = ${PROBLEM_FILE} ddm_setup.o c2fortran_kernel.o main_ddm_c.o
OBJS9 = ${PROBLEM_FILE} ddm_setup.o c2fortran_kernel.o main_ddm_dd_c.o
OBJS10 = ${PROBLEM_FILE} ddm_setup.o c2fortran_kernel.o main_qr_c.o

%.o : %.c
	${CC} -c ${CFLAGS} ${INCLUDEDIR} ${FDEFINE} $<
%.o : %.f90
	${FC} -c ${FFLAGS} ${INCLUDEDIR} ${FDEFINE} $<

all : ${OBJS1} ${OBJS2} ${OBJS3} ${OBJS4} ${OBJS5} ${OBJS6} ${OBJS7} ${OBJS8} ${OBJS9} ${OBJS10} 
	${FC} -o ${PROG1} ${OBJS1} ${FFLAGS} ${INCLUDEDIR} ${LDFLAGS}
	${FC} -o ${PROG2} ${OBJS2} ${FFLAGS} ${INCLUDEDIR} ${LDFLAGS} 
	${FC} -o ${PROG3} ${OBJS3} ${FFLAGS} ${INCLUDEDIR} ${LDFLAGS} 
	${FC} -o ${PROG4} ${OBJS4} ${FFLAGS} ${INCLUDEDIR} ${LDFLAGS} 
	${FC} -o ${PROG5} ${OBJS5} ${FFLAGS} ${INCLUDEDIR} ${LDFLAGS} 
	${CC} -o ${PROG6} ${OBJS6} ${CFLAGS} ${INCLUDEDIR} ${LDFLAGS}
	${CC} -o ${PROG7} ${OBJS7} ${CFLAGS} ${INCLUDEDIR} ${LDFLAGS} 
	${CC} -o ${PROG8} ${OBJS8} ${CFLAGS} ${INCLUDEDIR} ${LDFLAGS} 
	${CC} -o ${PROG9} ${OBJS9} ${CFLAGS} ${INCLUDEDIR} ${LDFLAGS} 
	${CC} -o ${PROG10} ${OBJS10} ${CFLAGS} ${INCLUDEDIR} ${LDFLAGS} 

clean :
	rm ${PROG1} ${PROG2} ${PROG3} ${PROG4} ${PROG5} ${PROG6} ${PROG7} ${PROG8} ${PROG9} ${PROG10} *o *mod *~
