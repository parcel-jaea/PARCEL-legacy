PROG = parcel_sparse.a 
MPI_dd_SRC = ./MPI/src/mpif_dd_sub1.o ./MPI/src/mpif_dd_func.o

all :
	cd MPI/src;make ;cd -
	make -f Makefile.double
	make -f Makefile.double_double
	ar rc ${PROG} ${MPI_dd_SRC} *o

clean :
	cd MPI/src;make clean ;cd -
	rm *mod *o *a *~
# 	make clean -f Makefile.double
# 	make clean -f Makefile.double_double

