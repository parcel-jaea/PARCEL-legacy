PARCEL 1.1 (2021/12/28)
The following points have been corrected in this update.
1. Supports DDM format.
2. Supports subdivision preconditioning for A64FX.
3. Supports setting of the number of blocks for block Jacobi preprocessing and additive Schwartz preprocessing.

PARCEL 1.01 (2020/5/20)
when making DIA format based on CRS, changed to execute without quick sort.
When making precondition ILU and DILU, diagonal elements don't accept zero in new version.
reduced executing time when run ILU in CRS and DIA.
reduced executing time when DILU in DIA.
Changed the allocating memory method of DILU in DIA and CRS.
exampleC and example2_C allowed to run only CRS in previous version. new version can run also DIA.
when setting overlap for additive schwartz to 0. We added proper executing processing.
initial resultant vecB in example 2,5,C,2_c differed to example and example3.We changed that to use same vecB.
Typical 5 solvers set by developer are available with command line argument.
We made parser in example_C and example2_C, for input from input list.
It was only possible to be executed CAGMRES as solver in example2 in previous version, new version provide all solver for asymmetric matrix in PARCEL.
We integrated name of variable which were different each example.

PARCEL 1.00 (2019/3/29)
Opening to the public


