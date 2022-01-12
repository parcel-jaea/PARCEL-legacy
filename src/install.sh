#!/bin/bash
WORKDIR=`pwd`
DIRNAME=`grep "INSTALLDIR" ../arch/make_config | cut -f 2 -d =`
if [ -e parcel_sparse.a ];then
    cd ../
    if [ ! -e ${DIRNAME} ];then
	mkdir  ${DIRNAME}
    fi
    cd ${DIRNAME}
    INSTALL_DIR=`pwd`
    echo "Install dir : " ${INSTALL_DIR}
    cd ${WORKDIR}
    mkdir ${INSTALL_DIR}/lib
    cp  parcel_sparse.a  ${INSTALL_DIR}/lib/
    mkdir ${INSTALL_DIR}/include
    cp  *mod  ${INSTALL_DIR}/include/
    cp  *\.h  ${INSTALL_DIR}/include/
    cp -rf ../example ${INSTALL_DIR}/
    cp ../arch/make_config ${INSTALL_DIR}/example
else
    echo "error : not found src/parcel_sparse.a  "
    echo "Please compile again "
fi

