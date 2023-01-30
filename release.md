PARCEL 1.2 (2023/mm/dd)
The GPU version is supported in this update.

PARCEL 1.1 (2021/12/28)
1. Domain decomposition method (DDM) data format is supported for stencil computation on structured grids.
2. The implementation is optimized for A64FX and SIMD-optimized fine-block preconditioner is added.
3. Setting of the number of blocks (threads) for the block Jacobi preconditioner and the additive Schwarz preconditioenr is supported.
4. An eigenvalue calculation routine by the communication avoiding Arnoldi method is added.
5. QR factorization routines by the classical Gram-Schmidt method, the modified Gram-Schmidt method, the tall skinny QR method, the Cholesky QR method, and the Cholesky QR2 method are added.


PARCEL 1.01 (2020/5/20)
1. Error handling is added for zero diagonal elements in the D-ILU preconditioner.
2. Memory allocation in the D-ILU preconditioner is changed to from static allocation to dynamic allocation.
3. The additive Schwarz preconditioner with zero overlap width is unified with the block Jacobi preconditioner.
4. A configuration file for the PGI compiler is added.
5. Sample programs are modified.
-- A sorting method in the data conversion routine from CRS to DIA is improved.
-- Template files for main algorithms are stored in each example.
-- The names of input parameters in each example are changed to common names.
-- The r.h.s. vector vecB in each example is changed to the same data.
-- ExampleC, 2_C are modified for DIA in addition to CSR.
-- All nonsymmetric solvers become available in example2.

PARCEL 1.00 (2019/3/29)



