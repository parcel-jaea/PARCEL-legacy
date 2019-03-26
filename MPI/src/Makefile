confdir=../../
include ${confdir}/arch/make_config
###
.SUFFIXES: .o .f .f90

.f.o :
	$(FC) -c -g ${INCLUDEDIR} $*.f

LIBDIR = ../lib
OBJSUB1 = mpif_dd_sub1.o
OBJSUB2 = mpif_dd_sub2.o
OBJFUNC = mpif_dd_func.o

LIB1 = $(LIBDIR)/libMpifdd1.a
LIB2 = $(LIBDIR)/libMpifdd2.a
TARGET = $(LIB1) $(LIB2) 

OBJS = $(OBJSUB1) $(OBJSUB2) $(OBJFUNC) 

all : $(TARGET)

$(LIB1) : $(OBJFUNC) $(OBJSUB1)
	ar r $@ $?
	ranlib $@

$(LIB2) : $(OBJFUNC) $(OBJSUB2)
	ar r $@ $?
	ranlib $@

clean :
	rm -f $(OBJS) $(TARGET)


