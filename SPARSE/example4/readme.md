# Example for matrices from MatrixMerket
main.f90 :: test program 
input_namelist :: input file for main.f90

In input_namelist, "-9" is substituted to unused parameters.

# Running tests
1. Download matrix data from "https://sparse.tamu.edu"
2. Generate input data, Matrix.coo, using mm2coo.py, mm_sym2coo.py
3. running test

* mpirun -n 4 ./run

# Templates
Template imput and outpu files for main algorithms are stored in ./templete/templete_input and ./templete/templete_output, respectively. The output files are generated using 8 MPI processes and 10 OpenMP threads. A sample matrix data "fv1.mtx" is used.
