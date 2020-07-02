1. Error handling is added for zero diagonal elements in the D-ILU preconditioner.
2. Memory allocation in the D-ILU preconditioner is changed to from static allocation to dynamic allocation. 
3. The additive Schwaltz preconditioner with zero overlap width is unified with the block Jacobi preconditioner.
4. A configuration file for the PGI compiler is added.
5. Sample programs are modified.
   - A soting method in the data conversion routine from CRS to DIA is improved.
   - Template files for main algorithms are stored in each example.
   - The names of input parameters in each example are changed to common names.
   - The r.h.s. vector vecB in each example is changed to the same data.
   - exampleC, 2_C are modified for DIA in addition to CSR.
   - All nonsymmetric solvers become available in example2.
