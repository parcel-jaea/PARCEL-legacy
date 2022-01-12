confdir=./
include ${confdir}/arch/make_config

all : 
	cd src ; make ;cd ../ 
clean :
	cd src ; make clean ;cd ../ 
	rm *~
install :
	cd src; . ./install.sh ;cd -  &> /dev/null
#	cd src; . ./install.sh ${INSTALLDIR} ;cd -  &> /dev/null
