# Example for large problem size using 64bit integer
main.f90 :: test program 
input_namelist :: input file for main.f90

In input_namelist, "-9" is substituted to unused parameters.

# Compile
This example requires 8byte integer. Remove comment out for #FDEDINE =-DPARCEL_INT8 in ../../arch/make_config before compilation.

# Running tests

* mpirun -n 4 ./sample1

# Templates
Template imput and outpu files for main algorithms are stored in ./templete/templete_input and ./templete/templete_output, respectively. The output files are generated using 8 MPI processes and 10 OpenMP threads.
