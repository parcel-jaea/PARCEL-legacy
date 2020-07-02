# Example using indirect interface for symmetric matrix in double precision
sample1.f90 :: test program
input_sample_poisson_namelist :: input file for sample1.f90

In input_namelist, "-9" is substituted to unused parameters.

# Running tests
*  mpirun -n 4 ./sample1

# Templates
Template imput and outpu files for main algorithms are stored in ./templete/templete_input and ./templete/templete_output, respectively. The output files are generated using 8 MPI processes and 10 OpenMP threads.
