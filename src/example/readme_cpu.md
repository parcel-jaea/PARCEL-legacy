# Example for large problem size using 64bit integer
input_namelist :: input file 

# Running tests

## CPU : CRS,DIA version 
*(Fortran)
 mpirun -n 4 run
*(C programming language)
 mpirun -n 4 run_c

## CPU : DDM version 
*(Fortran)
 mpirun -n 4 run_ddm
*(C programming language)
 mpirun -n 4 run_ddm_c

## CPU : double-double CRS,DIA version 
*(Fortran)
 mpirun -n 4 run_dd
*(C programming language)
 mpirun -n 4 run_dd_c

## CPU : double-double DDM version 
* (Fortran)
 mpirun -n 4 run_ddm_dd
*(C programming language)
 mpirun -n 4 run_ddm_dd_c

## CPU : QR decomposition 
*(Fortran)
 mpirun -n 4 run_qr
*(C programming language)
 mpirun -n 4 run_qr_c

# Templates
Template input and output files for main algorithms are stored in ./templete and ./template_Eigen, respectively. 
The output files are generated using 4 MPI processes and 20 OpenMP threads.
