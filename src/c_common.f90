module c_wrapper
  use iso_c_binding
  use PARCEL_integer
  implicit none

  ! GPU   
  type,bind(c) :: devicemem
     type(c_ptr) mem
     ! CUSPARCE
     type(c_ptr) descr ! 
     type(c_ptr) descr_add ! 袖領域を含んだ配列サイズ(SpMVの入力に使用) --> 確保した配列サイズと同じ
  end type devicemem

  ! type,bind(c) :: cuda_handle
  !    type(c_ptr) cublasHandle
  !    type(c_ptr) cusparseHandle
  ! end type cuda_handle

  type,bind(c) :: cublas_type
     type(c_ptr) cublasHandle
     type(devicemem) alpha
     type(devicemem) beta
  end type cublas_type

  ! CUSPARCE
  type,bind(c) :: cusparce_type
     type(c_ptr) cublasHandle
     type(c_ptr) cusparseHandle
     
     type(c_ptr) cusparseMatDescr
     type(c_ptr) matA
  
     !  // ------
     !  // BCSR形式
     integer blockDim

     integer mb
     integer nb
     integer nnzb 
     type(c_ptr) cusparseMatDescr_bsr 
     type(c_ptr) bsrRowPtr ! integer 
     type(c_ptr) bsrColInd ! integer 
     type(c_ptr) bsrVal    ! double

     integer ipreflag;
     integer precon_mb
     integer precon_nb
     integer precon_nnzb 
     type(c_ptr) cusparseMatDescr_precon_bsr 
     type(c_ptr) precon_bsrRowPtr ! integer 
     type(c_ptr) precon_bsrColInd ! integer 
     type(c_ptr) precon_bsrVal    ! double
     !  // ------

     type(c_ptr) cusparseMatDescr_precon
     type(c_ptr) precon_matA
     
     type(c_ptr) infoILU
     type(c_ptr) descrL
     type(c_ptr) infoL
     type(c_ptr) descrU
     type(c_ptr) infoU

     ! ------
     ! BCSR形式  
     type(c_ptr) descrILU_bcsr
     type(c_ptr) infoILU_bcsr
     type(c_ptr) descrL_bcsr
     type(c_ptr) infoL_bcsr
     type(c_ptr) descrU_bcsr
     type(c_ptr) infoU_bcsr
     ! ------  
     
     type(c_ptr) buffer
     type(c_ptr) precon_buffer
  end type cusparce_type

  ! type,bind(c) :: cusparce_vec
  !    type(c_ptr) mem 
  !    type(c_ptr) descr
  ! end type cusparce_vec

  type :: cusparce_crsMat
     type(devicemem) crsA
     type(devicemem) crsRow_ptr
     type(devicemem) crsCol
     integer(kind=intType) n
     integer(kind=intType) nCol
     integer(kind=intType) nnz
  end type cusparce_crsMat

  ! MAGMA
  type,bind(c) :: magma_handle
     type(c_ptr) opts
     type(c_ptr) queue
     integer matrix_type
  end type magma_handle
  
  type,bind(c) :: magma_memory
     type(c_ptr) host
     type(c_ptr) device
  end type magma_memory

  integer,parameter :: MAGMA_PRECON_NONE   = 0
  integer,parameter :: MAGMA_PRECON_JACOBI = 1
  integer,parameter :: MAGMA_PRECON_ILU    = 2
  integer,parameter :: MAGMA_PRECON_PARILU = 3
  
  ! --------
  ! 3D poisson 336x336x336 :: 6.613[s]
  ! integer,parameter :: CUDA_NUM_BLOCKS  = 128
  ! integer,parameter :: CUDA_NUM_THREADS = 512

  ! --------
  ! 3D poisson 336x336x336 :: 7.153[s]
  ! integer,parameter :: CUDA_NUM_BLOCKS  = 64
  ! integer,parameter :: CUDA_NUM_THREADS = 1024

  ! --------
  ! 3D poisson 336x336x336 :: 6.193[s]
  ! integer,parameter :: CUDA_NUM_BLOCKS  = 128
  ! integer,parameter :: CUDA_NUM_THREADS = 256
  ! --------
  ! 3D poisson 336x336x336 :: 6.187[s]
!  integer,parameter :: CUDA_NUM_BLOCKS  = 128
!  integer,parameter :: CUDA_NUM_BLOCKS  = 256
!  integer,parameter :: CUDA_NUM_THREADS = 1024
  ! integer,parameter :: CUDA_NUM_BLOCKS  = 512
  ! integer,parameter :: CUDA_NUM_THREADS = 512

!  integer,parameter :: CUDA_NUM_BLOCKS  = 512
!  integer,parameter :: CUDA_NUM_BLOCKS  = 2048
  integer,parameter :: CUDA_NUM_BLOCKS  = 4096
  integer,parameter :: CUDA_NUM_THREADS = 512

end module c_wrapper

module cu_common_wrapper
  use iso_c_binding
  use PARCEL_integer
  use c_wrapper
  implicit none

  ! GPU関数インターフェース
  interface
     subroutine init_cudahandle()&
          bind(c,name="init_cudahandle")
       implicit none
     end subroutine init_cudahandle

     subroutine free_cudahandle()&
          bind(c,name="free_cudahandle")
       implicit none
     end subroutine free_cudahandle

     subroutine set_GPU(rank) &
          bind(c,name="set_GPU")
       implicit none
       integer,value :: rank
     end subroutine set_GPU

     subroutine set_cublas(cublas_solver) &
          bind(c,name="set_cublas")
       import :: cublas_type
       implicit none
       type(cublas_type)   :: cublas_solver
     end subroutine set_cublas

     subroutine free_cublas(cublas_solver) &
          bind(c,name="free_cublas")
       import :: cublas_type
       implicit none
       type(cublas_type)   :: cublas_solver
     end subroutine free_cublas

     subroutine cudaDeviceSynchronize() &
          bind(c,name="cudaDeviceSynchronize_")
       implicit none
     end subroutine cudaDeviceSynchronize

     subroutine cudaFree(a_d) &
          bind(c,name="cudaFree_")
       import :: devicemem
       implicit none
       type(devicemem) a_d
     end subroutine cudaFree

     subroutine cudaMalloc_int(a_d,n) &
          bind(c,name="cudaMalloc_int")
       import :: intType
       import :: devicemem
       implicit none
       type(devicemem) a_d
       integer(kind=intType),value :: n
     end subroutine cudaMalloc_int

     subroutine cudaMalloc_PARCEL_INT(a_d,n) &
          bind(c,name="cudaMalloc_PARCEL_INT")
       import :: intType
       import :: devicemem
       implicit none
       type(devicemem) a_d
       integer(kind=intType),value :: n
     end subroutine cudaMalloc_PARCEL_INT

     subroutine cudaMalloc_FP64(a_d,n) &
          bind(c,name="cudaMalloc_FP64")
       import :: intType
       import :: devicemem
       implicit none
       type(devicemem) a_d
       integer(kind=intType),value :: n
     end subroutine cudaMalloc_FP64

     subroutine cudaMalloc_bf16(a_d,n) &
          bind(c,name="cudaMalloc_bf16")
       import :: intType
       import :: devicemem
       implicit none
       type(devicemem) a_d
       integer(kind=intType),value :: n
     end subroutine cudaMalloc_bf16

     subroutine cudaMemcpy_h2d_FP64(&
          a_h,a_d,n,&
          offset_h, &
          offset_d) &
          bind(c,name="cudaMemcpy_h2d_FP64")
       import :: intType
       import :: devicemem
       implicit none
       real*8 a_h(*)
       type(devicemem) a_d
       ! valueを指定すれば値渡しも可能
       integer(kind=intType),value :: n
       integer(kind=intType),value ::  offset_h
       integer(kind=intType),value ::  offset_d
     end subroutine CUDAmemcpy_h2d_FP64

     subroutine cudaMemcpy_d2d_FP64(&
          a_d,b_d,n,&
          offset_a, &
          offset_b  &
          ) &
          bind(c,name="cudaMemcpy_d2d_FP64")
       import :: intType
       import :: devicemem
       implicit none
       type(devicemem) a_d
       type(devicemem) b_d
       ! valueを指定すれば値渡しも可能
       integer(kind=intType),value :: n
       integer(kind=intType),value ::  offset_a
       integer(kind=intType),value ::  offset_b
     end subroutine cudaMemcpy_d2d_FP64

     subroutine cudaMemcpy_d2d_int(&
          a_d,b_d,n,&
          offset_a, &
          offset_b  &
          ) &
          bind(c,name="cudaMemcpy_d2d_int")
       import :: intType
       import :: devicemem
       implicit none
       type(devicemem) a_d
       type(devicemem) b_d
       ! valueを指定すれば値渡しも可能
       integer(kind=intType),value :: n
       integer(kind=intType),value ::  offset_a
       integer(kind=intType),value ::  offset_b
     end subroutine cudaMemcpy_d2d_int

     subroutine cudaMemcpy_d2d_PARCEL_INT(&
          a_d,b_d,n,&
          offset_a, &
          offset_b  &
          ) &
          bind(c,name="cudaMemcpy_d2d_PARCEL_INT")
       import :: intType
       import :: devicemem
       implicit none
       type(devicemem) a_d
       type(devicemem) b_d
       ! valueを指定すれば値渡しも可能
       integer(kind=intType),value :: n
       integer(kind=intType),value :: offset_a
       integer(kind=intType),value :: offset_b
     end subroutine cudaMemcpy_d2d_PARCEL_INT

     subroutine cudaTransfer_h2d_FP64toBF16(a_h,a_d,n,mode) &
          bind(c,name="cudaTransfer_h2d_FP64toBF16")
       import :: devicemem
       import :: intType
       implicit none
       real*8 a_h(*)
       type(devicemem) a_d
       integer(kind=intType),value :: n
       integer,value :: mode
     end subroutine cudaTransfer_h2d_FP64toBF16

     subroutine cuda_view_FP64(&
          n,a_d &
          ) &
          bind(c,name="cuda_view_FP64")
       import :: devicemem
       import :: intType
       implicit none
       integer(kind=intType),value :: n
       type(devicemem) a_d
     end subroutine cuda_view_FP64

     subroutine cuda_view_PARCEL_INT(&
          n,a_d &
          ) &
          bind(c,name="cuda_view_PARCEL_INT")
       import :: devicemem
       import :: intType
       implicit none
       integer(kind=intType),value :: n
       type(devicemem) a_d
     end subroutine cuda_view_PARCEL_INT

#if 0
     subroutine ptr_real_view(&
          n,aptr &
          ) &
          bind(c,name="cuda_ptr_double_view")
       import :: intType
       import :: c_ptr
       implicit none
       integer(kind=intType),value :: n
       type(c_ptr)  aptr
     end subroutine ptr_real_view

     subroutine ptr_parcel_int_view(&
          n,aptr &
          ) &
          bind(c,name="cuda_ptr_parcel_int_view")
       import :: intType
       import :: c_ptr
       implicit none
       integer(kind=intType),value :: n
       type(c_ptr)  aptr
     end subroutine ptr_parcel_int_view

     subroutine cudaMemcpy_d2dptr(&
          a_d,aptr,typeSize,n &
          ) &
          bind(c,name="cudaMemcpy_d2dptr")
       import :: devicemem
       import :: intType
       import :: c_ptr
       implicit none
       type(devicemem) a_d
       type(c_ptr)  aptr
       integer,value :: typeSize
       integer(kind=intType),value :: n
     end subroutine CudaMemcpy_d2dptr

     subroutine cudaMemcpy_dptr2d(&
          aptr,a_d,typeSize,n &
          ) &
          bind(c,name="cudaMemcpy_dptr2d")
       import :: devicemem
       import :: intType
       import :: c_ptr
       implicit none
       type(c_ptr)  aptr
       type(devicemem) a_d
       integer,value :: typeSize
       integer(kind=intType),value :: n
     end subroutine CudaMemcpy_dptr2d
#endif
     
     subroutine cudaMemcpy_h2d_int(&
          a_h,a_d,n,&
          offset_h, &
          offset_d  &
          ) &
          bind(c,name="cudaMemcpy_h2d_int")
       import :: devicemem
       import :: intType
       implicit none
       integer a_h(*)
       type(devicemem) a_d
       ! valueを指定すれば値渡しも可能
       integer(kind=intType),value :: n
       integer(kind=intType),value ::  offset_h
       integer(kind=intType),value ::  offset_d
     end subroutine cudaMemcpy_h2d_int

     subroutine cudaMemcpy_h2d_PARCEL_INT(&
          a_h,a_d,n,&
          offset_h, &
          offset_d  &
          ) &
          bind(c,name="cudaMemcpy_h2d_PARCEL_INT")
       import :: devicemem
       import :: intType
       implicit none
       integer(kind=intType) a_h(*)
       type(devicemem) a_d
       ! valueを指定すれば値渡しも可能
       integer(kind=intType),value :: n
       integer(kind=intType),value ::  offset_h
       integer(kind=intType),value ::  offset_d
     end subroutine cudaMemcpy_h2d_PARCEL_INT

     subroutine cudaMemcpy_d2h_FP64(&
          a_d,a_h,n,&
          offset_d, &
          offset_h  &
          ) &
          bind(c,name="cudaMemcpy_d2h_FP64")
       import :: devicemem
       import :: intType
       implicit none
       real*8 a_h(*)
       type(devicemem) a_d
       ! valueを指定すれば値渡しも可能
       integer(kind=intType),value :: n
       integer(kind=intType),value ::  offset_d
       integer(kind=intType),value ::  offset_h
     end subroutine CUDAmemcpy_d2h_FP64

     subroutine cudaMemcpy_d2h_int(&
          a_d,a_h,n,&
          offset_d, &
          offset_h  &
          ) &
          bind(c,name="cudaMemcpy_d2h_int")
       import :: devicemem
       import :: intType
       implicit none
       integer a_h(*)
       type(devicemem) a_d
       ! valueを指定すれば値渡しも可能
       integer(kind=intType),value :: n
       integer(kind=intType),value ::  offset_d
       integer(kind=intType),value ::  offset_h
     end subroutine CUDAmemcpy_d2h_int

     subroutine cudaMemcpy_d2h_PARCEL_INT(&
          a_d,a_h,n, &
          offset_d, &
          offset_h  &
          ) &
          bind(c,name="cudaMemcpy_d2h_PARCEL_INT")
       import :: devicemem
       import :: intType
       implicit none
       integer(kind=intType) a_h(*)
       type(devicemem) a_d
       ! valueを指定すれば値渡しも可能
       integer(kind=intType),value :: n
       integer(kind=intType),value ::  offset_d
       integer(kind=intType),value ::  offset_h
     end subroutine CUDAmemcpy_d2h_PARCEL_INT

     subroutine cublas_init_DgemmBatched_gemm_dot_local_FP64(&
          batchCount, &
          d_A_, offset_A, lda,  &
          d_Aptr_ &
          ) &
          bind(c,name="cublas_init_DgemmBatched_gemm_dot_local")
       import :: intType
       import :: devicemem
       implicit none
       integer(kind=intType),value ::  batchCount
       type(devicemem) d_A_
       integer(kind=intType),value :: offset_A
       integer(kind=intType),value :: lda
       type(devicemem) d_Aptr_

     end subroutine cublas_init_DgemmBatched_gemm_dot_local_FP64

     subroutine cublas_free_DgemmBatched_gemm_dot_local_FP64( &
          d_Aptr_) &
          bind(c,name="cublas_free_DgemmBatched_gemm_dot_local")
       import :: intType
       import :: devicemem
       implicit none
       type(devicemem) d_Aptr_
     end subroutine cublas_free_DgemmBatched_gemm_dot_local_FP64

     subroutine cublas_DgemmBatched_gemm_dot_local_FP64( &
					  cublas_solver, &
					  batchCount, &
					  nCol, localnRow, &
					  Aptr_, ldAptr, &
					  Bptr_, ldBptr, &
					  Cptr_, ldCptr, &
					  Cmat_,  &
					  reduce_mat_ ) &
          bind(c,name="cublas_DgemmBatched_gemm_dot_local")
       import :: cublas_type
       import :: intType
       import :: devicemem
       implicit none

       type(cublas_type)  :: cublas_solver
       integer(kind=intType),value ::  batchCount
       integer(kind=intType),value ::  nCol
       integer(kind=intType),value ::  localnRow
       type(devicemem)Aptr_
       integer(kind=intType),value ::  ldAptr
       type(devicemem)Bptr_
       integer(kind=intType),value ::  ldBptr
       type(devicemem)Cptr_
       integer(kind=intType),value ::  ldCptr

       type(devicemem)Cmat_  
       type(devicemem)reduce_mat_ 
    end subroutine cublas_DgemmBatched_gemm_dot_local_FP64

    subroutine cublas_DgemmBatched_gemm_dot_FP64( &
				          f_comm, &
					  cublas_solver, &
					  batchCount, &
					  nCol, localnRow, &
					  Aptr_, ldAptr, &
					  Bptr_, ldBptr, &
					  Cptr_, ldCptr, &
					  Cmat_,  &
					  reduce_mat_ ) &
          bind(c,name="cublas_DgemmBatched_gemm_dot")
       import :: cublas_type
       import :: intType
       import :: devicemem
       implicit none

       integer,value :: f_comm
       type(cublas_type)  :: cublas_solver
       integer(kind=intType),value ::  batchCount
       integer(kind=intType),value ::  nCol
       integer(kind=intType),value ::  localnRow
       type(devicemem)Aptr_
       integer(kind=intType),value ::  ldAptr
       type(devicemem)Bptr_
       integer(kind=intType),value ::  ldBptr
       type(devicemem)Cptr_
       integer(kind=intType),value ::  ldCptr

       type(devicemem)Cmat_  
       type(devicemem)reduce_mat_ 
     end subroutine cublas_DgemmBatched_gemm_dot_FP64

    subroutine swap_ptr_ptr(d_Aptr_,d_Bptr_)
      import :: devicemem
      implicit none
       type(devicemem) d_Aptr_
       type(devicemem) d_Bptr_
    end subroutine swap_ptr_ptr

     subroutine cusparseDestroyDnVec(x) &
          bind(c,name="cusparseDestroyDnVec_")
       import :: devicemem
       implicit none
       type(devicemem) x
     end subroutine cusparseDestroyDnVec

     subroutine cusparseCreateDnVec_FP64(x,n,nAdd,offset) &
          bind(c,name="cusparseCreateDnVec_FP64")
       import :: intType
       import :: devicemem
       implicit none
       type(devicemem) x
       integer(kind=intType),value :: n
       integer(kind=intType),value :: nAdd
       integer(kind=intType),value :: offset
     end subroutine cusparseCreateDnVec_FP64

     subroutine cuda_add_val_FP64(numGPUblock,numGPUthread,val,x,offset_x) &
          bind(c,name="cuda_add_val")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       real*8,value :: val
       type(devicemem) x
       integer(kind=intType),value ::  offset_x
     end subroutine cuda_add_val_FP64

     subroutine cuda_sub_val_FP64(numGPUblock,numGPUthread,val,x,offset_x) &
          bind(c,name="cuda_sub_val")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       real*8,value :: val
       type(devicemem) x
       integer(kind=intType),value ::  offset_x
     end subroutine cuda_sub_val_FP64

     subroutine cuda_mul_val_FP64(numGPUblock,numGPUthread,val,x,offset_x) &
          bind(c,name="cuda_mul_val")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       real*8,value :: val
       type(devicemem) x
       integer(kind=intType),value ::  offset_x
     end subroutine cuda_mul_val_FP64

     subroutine cuda_div_val_FP64(numGPUblock,numGPUthread,val,x,offset_x) &
          bind(c,name="cuda_div_val")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       real*8,value :: val
       type(devicemem) x
       integer(kind=intType),value ::  offset_x
     end subroutine cuda_div_val_FP64

     subroutine cuda_sqrt_val_FP64(numGPUblock,numGPUthread,x,offset_x) &
          bind(c,name="cuda_sqrt_val")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       type(devicemem) x
       integer(kind=intType),value ::  offset_x
     end subroutine cuda_sqrt_val_FP64

     subroutine set_vec_val_FP64_GPU(&
          numGPUblock,&
          numGPUthread,&
          n, &
          val, &
          x,offset_x &
          ) &
          bind(c,name = "cuda_set_vec_val")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value ::  n
       real*8,value :: val
       type(devicemem) x
       integer(kind=intType),value ::  offset_x
     end subroutine set_vec_val_FP64_GPU

  end interface

  ! CPU関数インターフェース
  interface
     subroutine malloc_MPI_Status(a_d,n) &
          bind(c,name="malloc_MPI_Status")
       import :: devicemem
       implicit none
       type(devicemem) a_d
       integer,value :: n
     end subroutine malloc_MPI_Status

     subroutine malloc_MPI_Request(a_d,n) &
          bind(c,name="malloc_MPI_Request")
       import :: devicemem
       implicit none
       type(devicemem) a_d
       integer,value :: n
     end subroutine malloc_MPI_Request

     subroutine malloc_MPI_Comm(a_d,n) &
          bind(c,name="malloc_MPI_Comm")
       import :: devicemem
       implicit none
       type(devicemem) a_d
       integer,value :: n
     end subroutine malloc_MPI_Comm

     subroutine MPI_Comm_f2c(a_d,f_comm) &
          bind(c,name="C_MPI_Comm_f2c")
       import :: devicemem
       implicit none
       type(devicemem) a_d
       integer,value :: f_comm
!       integer f_comm
     end subroutine MPI_Comm_f2c

     subroutine free(a_d) &
          bind(c,name="C_free")
       import :: devicemem
       implicit none
       type(devicemem) a_d
     end subroutine free

     subroutine C_MPI_Isend_FP64(&
          a_d, &
          ndata, &
          rank, &
          tag, &
          comm, &
          req, &
          id,offset &
          ) &
          bind(c,name="C_MPI_Isend_FP64")
       import :: devicemem
       implicit none
       type(devicemem)  a_d
       integer,value :: ndata
       integer,value :: rank
       integer,value :: tag
       type(devicemem)  comm
       type(devicemem)  req
       integer,value :: id
       integer,value :: offset
     end subroutine C_MPI_Isend_FP64

     subroutine C_MPI_Irecv_FP64(&
          a_d, &
          ndata, &
          rank, &
          tag, &
          comm, &
          req, &
          id,offset &
          ) &
          bind(c,name="C_MPI_Irecv_FP64")
       import :: devicemem
       implicit none
       type(devicemem)  a_d
       integer,value :: ndata
       integer,value :: rank
       integer,value :: tag
       type(devicemem)  comm
       type(devicemem)  req
       integer,value :: id
       integer,value :: offset
     end subroutine C_MPI_Irecv_FP64

     subroutine C_MPI_Isend_FP64buf_FP64(&
          a_d, &
          ndata, &
          rank, &
          tag, &
          comm, &
          req, &
          id,offset &
          ) &
          bind(c,name="C_MPI_Isend_FP64buf_FP64")
       import :: devicemem
       implicit none
       real*8  a_d(*)
       integer,value :: ndata
       integer,value :: rank
       integer,value :: tag
       type(devicemem)  comm
       type(devicemem)  req
       integer,value :: id
       integer,value :: offset
     end subroutine C_MPI_Isend_FP64buf_FP64

     subroutine C_MPI_Irecv_FP64buf_FP64(&
          a_d, &
          ndata, &
          rank, &
          tag, &
          comm, &
          req, &
          id,offset &
          ) &
          bind(c,name="C_MPI_Irecv_FP64buf_FP64")
       import :: devicemem
       implicit none
       real*8  a_d(*)
       integer,value :: ndata
       integer,value :: rank
       integer,value :: tag
       type(devicemem)  comm
       type(devicemem)  req
       integer,value :: id
       integer,value :: offset
     end subroutine C_MPI_Irecv_FP64buf_FP64
 
     subroutine C_MPI_Wait(&
          req,stat,id &
          ) &
          bind(c,name="C_MPI_Wait")
       import :: devicemem
       implicit none
       type(devicemem)  req
       type(devicemem)  stat
       integer,value :: id
     end subroutine C_MPI_Wait

  end interface

  ! CUSPARCEインターフェース
  interface 
     subroutine cusparse_init( &
          cusparce, &
          d_val, &
          d_col, &
          d_row, &
          nRow, &
          nCol, &
          nnz, &
          d_precon_val, &
          d_precon_col, &
          d_precon_row, &
          precon_nRow, &
          precon_nCol, &
          precon_nnz &
          ) &
          bind(c,name="cusparse_init")
       import :: cusparce_type
       import :: devicemem
       import :: intType
       implicit none
       type(cusparce_type) cusparce
       type(devicemem)   d_val
       type(devicemem)   d_col
       type(devicemem)   d_row
       integer(kind=intType), value :: nRow
       integer(kind=intType), value :: nCol
       integer(kind=intType), value :: nnz
       type(devicemem)   d_precon_val
       type(devicemem)   d_precon_col
       type(devicemem)   d_precon_row
       integer(kind=intType), value :: precon_nRow
       integer(kind=intType), value :: precon_nCol
       integer(kind=intType), value :: precon_nnz

     end subroutine cusparse_init

     subroutine cusparse_SpMV( &
          cusparce,&
          vecx, offset_x, &
          vecy, offset_y  &
          ) &
          bind(c,name="cusparse_SpMV")
       import :: cusparce_type
       import :: devicemem
       import :: intType
       implicit none
       type(cusparce_type) cusparce
       type(devicemem) vecx
       integer(kind=intType), value :: offset_x
       type(devicemem) vecy
       integer(kind=intType), value :: offset_y

     end subroutine cusparse_SpMV

     subroutine Neumann_SpMV( &
          cusparce,&
          vecx, offset_x,&
          vecy, offset_y&
          ) &
          bind(c,name="Neumann_SpMV")
       import :: cusparce_type
       import :: devicemem
       import :: intType
       implicit none
       type(cusparce_type) cusparce
       type(devicemem) vecx
       integer(kind=intType), value :: offset_x
       type(devicemem) vecy
       integer(kind=intType), value :: offset_y
     end subroutine Neumann_SpMV

     subroutine cusparse_precon( &
          cusparce, &
          vecr, offset_r, &
          vecz, offset_z, &
          d_precon_val, &
          d_precon_col, &
          d_precon_row, &
          ! d_precon_bsrVal, &
          ! d_precon_bsrColInd, &
          ! d_precon_bsrRowPtr, &
          precon_N, &
          precon_nnz &
          ) &
          bind(c,name="cusparse_precon")
       import :: cusparce_type
       import :: devicemem
       import :: intType
       implicit none
       type(cusparce_type) cusparce
       type(devicemem) vecr
       integer(kind=intType), value :: offset_r
       type(devicemem) vecz
       integer(kind=intType), value :: offset_z
       type(devicemem) d_precon_val
       type(devicemem) d_precon_col
       type(devicemem) d_precon_row
       ! type(devicemem) d_precon_bsrVal
       ! type(devicemem) d_precon_bsrColInd
       ! type(devicemem) d_precon_bsrRowPtr       
       integer(kind=intType), value :: precon_N
       integer(kind=intType), value :: precon_nnz 
     end subroutine cusparse_precon

     subroutine cusparse_free( &
          cusparce &
          ) &
          bind(c,name="cusparse_free")
       import :: cusparce_type
       implicit none
       type(cusparce_type) cusparce
     end subroutine cusparse_free

  end interface

contains

#if 0
  subroutine cusparseDnVecDescr_reset(x,n,nAdd)
    implicit none
    type(devicemem) x
    integer,value :: n
    integer,value :: nAdd

    call cusparseDestroyDnVec(x) 
    call cusparseCreateDnVec_FP64(x,n,nAdd) 

  end subroutine cusparseDnVecDescr_reset
#endif

end module cu_common_wrapper

module magma_common_wrapper
  use iso_c_binding
  use PARCEL_integer
  use c_wrapper
  implicit none

  interface
     subroutine set_magma_GPU(rank) &
          bind(c,name="set_magma_GPU_")
       implicit none
       integer,value :: rank
     end subroutine set_magma_GPU
     
     subroutine magma_init(handle,rank,matrix_type)&
          bind(c,name="magma_init_")
       import :: magma_handle
       implicit none
       type(magma_handle) handle
       integer,value :: rank
       integer,value :: matrix_type
     end subroutine magma_init
     
     subroutine magma_final(handle) &
          bind(c,name="magma_final_")
       import :: magma_handle
       implicit none
       type(magma_handle) handle
     end subroutine magma_final

     subroutine show_gpu_id() &
          bind(c,name="show_gpu_id_")
       implicit none
     end subroutine show_gpu_id

     subroutine magma_queue_sync(handle) &
          bind(c,name="magma_queue_sync_")
       import :: magma_handle
       implicit none
       type(magma_handle) handle
     end subroutine magma_queue_sync

     subroutine magmaMemcpy_vec_d2h_FP64(&
          handle, &
          mem, &
          vec, &
          n &
       ) &
          bind(c,name="magmaMemcpy_vec_d2h_FP64_")
       import :: magma_handle
       import :: magma_memory
       import :: intType
       implicit none
       type(magma_handle) handle
       type(magma_memory) mem
       real*8 vec(*)
       integer(kind=intType),value :: n
     end subroutine magmaMemcpy_vec_d2h_FP64

     subroutine magmaMemcpy_vec_h2d_FP64(&
          handle, &
          mem, &
          vec, &
          n &
          ) &
          bind(c,name="magmaMemcpy_vec_h2d_FP64_")
       import :: magma_handle
       import :: magma_memory
       import :: intType
       implicit none
       type(magma_handle) handle
       type(magma_memory) mem
       real*8 vec(*)
       integer(kind=intType),value :: n
     end subroutine magmaMemcpy_vec_h2d_FP64
     
     subroutine magma_set_vec( &
          handle, &
          mem, &
          n) &
          bind(c,name="magma_set_vec_")
       import :: magma_handle
       import :: magma_memory
       import :: intType
       implicit none
       type(magma_handle) handle
       type(magma_memory) mem
       integer(kind=intType),value :: n
     end subroutine magma_set_vec

     subroutine magma_set_precon( &
          handle, &
          A, &
          x, &
          precon ) &
          bind(c,name="magma_set_precon_")
       import :: magma_handle
       import :: magma_memory
       implicit none
       type(magma_handle) handle
       type(magma_memory) A
       type(magma_memory) x
       integer,value :: precon
     end subroutine magma_set_precon
     
     subroutine magma_precon_solver( &
          handle, &
          M, &
          r, &
          z, &
          work &
          ) &
          bind(c,name="magma_precon_solver_")
       import :: magma_handle
       import :: magma_memory
       implicit none
       type(magma_handle) handle
       type(magma_memory) M
       type(magma_memory) r
       type(magma_memory) z
       type(magma_memory) work
     end subroutine magma_precon_solver
     
     subroutine magma_spmv( &
          handle, &
          A, &
          x, &
          y,&
          nRow, &
          nCol  &
          ) &
          bind(c,name="magma_spmv_")
       import :: magma_handle
       import :: magma_memory
       import :: intType
       implicit none
       type(magma_handle) handle
       type(magma_memory) A
       type(magma_memory) x
       type(magma_memory) y
       integer(kind=intType),value :: nRow
       integer(kind=intType),value :: nCol
     end subroutine magma_spmv

#if 0
     subroutine  magma_precon_solver( &
          type(magma_handle) handle, &
          type(magma_handle) M, &
          type(magma_handle) r, &
          type(magma_handle) z, &
          type(magma_handle) work &
          ) &
          bind(c,name="magma_precon_solver_")
       import :: magma_handle
       import :: magma_memory
       import :: intType
       implicit none
       type(magma_handle) handle
       type(magma_memory) M
       type(magma_memory) r
       type(magma_memory) z
       type(magma_memory) work
     end subroutine magma_precon_solver
#endif

     subroutine magma_set_mat( &
		    handle, &
		    mem, &
		    row,col,val, &
                    nRow, &
                    nCol) &
                    bind(c,name="magma_set_mat_")
       import :: magma_handle
       import :: magma_memory
       import :: intType
       implicit none
       type(magma_handle) handle
       type(magma_memory) mem
       integer row(*)
       integer col(*)
       real*8 val(*)
       integer(kind=intType),value :: nRow
       integer(kind=intType),value :: nCol
       
     end subroutine magma_set_mat
     
     subroutine magma_dmfree(&
          handle, &
          mem) &
          bind(c,name="magma_dmfree_")          
       import :: magma_handle
       import :: magma_memory
       type(magma_handle) handle
       type(magma_memory) mem
       
     end subroutine magma_dmfree

  end interface
  
contains
  
end module magma_common_wrapper


