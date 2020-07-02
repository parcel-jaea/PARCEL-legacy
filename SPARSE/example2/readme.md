# Example for symmetric/non-symmetric matrices in double precision
main.f90 :: test program 
input_namelist :: input file for main.f90

In input_namelist, "-9" is substituted to unused parameters.

# Running tests
*  mpirun -n 4 ./run

# Templates
Template imput and outpu files for main algorithms are stored in ./templete/templete_input and ./templete/templete_output, respectively. The output files are generated using 8 MPI processes and 10 OpenMP threads.
