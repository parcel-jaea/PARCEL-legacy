# Example in C using indirect interface for matrices in double precision
C_wrapper.f90 :: wrapper for C
main.c :: example program
input_list :: input file for main.c

In input_list, "-9" is substituted to unused parameters.

# Running tests

* mpirun -n 4 ./run

# Templates
Template imput and outpu files for main algorithms are stored in ./templete/templete_input and ./templete/templete_output, respectively. The output files are generated using 8 MPI processes and 10 OpenMP threads.
