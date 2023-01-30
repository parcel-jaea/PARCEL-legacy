# Example for large problem size using 64bit integer
input_namelist :: input file 

# Running tests

## GPU : CRS,DIA version 
*(Fortran)
 mpirun -n 4 run
*(C programming language)
 mpirun -n 4 run_c

## GPU : DDM version 
*(Fortran)
 mpirun -n 4 run_ddm
*(C programming language)
 mpirun -n 4 run_ddm_c

# Templates
Template input and output files for main algorithms are stored in ./templete, respectively. 
The output files are generated using 4 MPI processes and 12 OpenMP threads.
