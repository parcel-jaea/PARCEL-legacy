# Example for large problem size using 64bit integer
input_namelist :: input file 

# Running tests

* mpirun -n 4 ./run

# Templates
Template input and output files for main algorithms are stored in ./templete and ./template_Eigen, respectively. 
The output files are generated using 4 MPI processes and 20 OpenMP threads.
