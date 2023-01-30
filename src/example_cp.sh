#!/bin/bash

INSTALL_DIR=$1
CUDA_FLAG=`grep "CUDAFLAG" ../arch/make_config | cut -f 2 -d = | tr -d ' ' `

CPU_ONLY_FILE=(
    "main_QR.f90"
    "main_dd.f90"
    "main_dd_c.c"
    "main_ddm_dd.f90"
    "main_ddm_dd_c.c"
    "main_qr_c.c"
)
# CPU_ONLY_TEMPLATE=(
#     "STDOUT_run_qr_c"
#     "STDOUT_run_qr"
#     "STDOUT_run_ddm_dd_c"
#     "STDOUT_run_ddm_dd"
#     "STDOUT_run_dd_c"
#     "STDOUT_run_dd"
# )



# for file in "${CPU_ONLY_FILE[@]}" ; do
# echo "[ ${file} ]"
# done

# main_ddm.f90
# main_ddm_c.c
# main_direct.f90
# main_direct_c.c

cp -rf ./example ${INSTALL_DIR}/ 
if [ ${CUDA_FLAG} == "ON" ];then
    # GPU
    # echo " GPU "
    for file in "${CPU_ONLY_FILE[@]}" ; do
	rm ${INSTALL_DIR}/example/${file}
    done

    # for file in "${CPU_ONLY_TEMPLATE[@]}" ; do
    # 	rm ${INSTALL_DIR}/example/template/${file}
    # done
    rm -rf ${INSTALL_DIR}/example/template_Eigen_cpu

    rm ${INSTALL_DIR}/example/Makefile.cpu
    mv ${INSTALL_DIR}/example/Makefile.gpu ${INSTALL_DIR}/example/Makefile

    rm ${INSTALL_DIR}/example/readme_cpu.md
    mv ${INSTALL_DIR}/example/readme_gpu.md ${INSTALL_DIR}/example/readme.md

    rm -rf ${INSTALL_DIR}/example/template_cpu
    mv ${INSTALL_DIR}/example/template_gpu ${INSTALL_DIR}/example/template
else
    # CPU
    # echo " CPU "
    rm -rf ${INSTALL_DIR}/example/template_gpu

    rm ${INSTALL_DIR}/example/Makefile.gpu
    mv ${INSTALL_DIR}/example/Makefile.cpu ${INSTALL_DIR}/example/Makefile

    rm ${INSTALL_DIR}/example/readme_gpu.md
    mv ${INSTALL_DIR}/example/readme_cpu.md ${INSTALL_DIR}/example/readme.md

    rm -rf ${INSTALL_DIR}/example/template_gpu
    mv ${INSTALL_DIR}/example/template_cpu ${INSTALL_DIR}/example/template
fi
