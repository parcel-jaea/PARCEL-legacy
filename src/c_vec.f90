
module cu_Vec

  use PARCEL_integer
  use c_wrapper
  implicit none

  interface

     subroutine calcQ_VRinv_FP64_GPU( &
		 s, &
		 n, &
		 V_   , offset_V   , ldv, &
		 Q_   , offset_Q   , ldq, &
		 Rinv_, offset_Rinv, ldrinv &
          )bind(c,name="cuda_calcQ_VRinv")
       import :: intType
       import :: devicemem
       implicit none
 
       integer(kind=intType),value ::  s
       integer(kind=intType),value ::  n

       type(devicemem) V_   
       integer(kind=intType),value ::  offset_V   
       integer(kind=intType),value ::  ldv
       
       type(devicemem) Q_   
       integer(kind=intType),value ::  offset_Q   
       integer(kind=intType),value ::  ldq

       type(devicemem) Rinv_
       integer(kind=intType),value ::  offset_Rinv
       integer(kind=intType),value ::  ldrinv
       
     end subroutine calcQ_VRinv_FP64_GPU

     subroutine CGS_Orthonormalization_FP64_GPU( &
          V, offset_V, ldv,&
          Q, offset_Q, ldq,&
          R, offset_R, ldr,&
          n,s,i &
       ) bind(c,name="cuda_CGS_Orthonormalization")
       import :: intType
       import :: devicemem
       implicit none

       type(devicemem) V
       integer(kind=intType),value :: offset_V
       integer(kind=intType),value :: ldv

       type(devicemem) Q
       integer(kind=intType),value :: offset_Q
       integer(kind=intType),value :: ldq

       type(devicemem) R
       integer(kind=intType),value :: offset_R
       integer(kind=intType),value :: ldr

       integer(kind=intType),value :: n
       integer(kind=intType),value :: s
       integer(kind=intType),value :: i

     end subroutine CGS_Orthonormalization_FP64_GPU

     subroutine calc_dgemmNN_dot_FP64_GPU( &
          icomm, &
          cublas_solver, &
          m, &
          n, &
          k, &
          alpha_, offset_alpha, &
          beta_, offset_beta, &
          A_, offset_A, lda, &
          B_, offset_B, ldb, &
          C_, offset_C, ldc &
          ) bind(c,name="cuda_calc_dgemmNN_dot")
       import :: intType
       import :: devicemem
       import :: cublas_type
       implicit none

       integer,value :: icomm
       type(cublas_type) :: cublas_solver
       integer(kind=intType),value :: m
       integer(kind=intType),value :: n
       integer(kind=intType),value :: k

       type(devicemem) alpha_
       integer(kind=intType),value :: offset_alpha

       type(devicemem) beta_
       integer(kind=intType),value :: offset_beta

       type(devicemem) A_
       integer(kind=intType),value :: offset_A
       integer(kind=intType),value :: lda

       type(devicemem) B_
       integer(kind=intType),value :: offset_B
       integer(kind=intType),value :: ldb

       type(devicemem) C_
       integer(kind=intType),value :: offset_C
       integer(kind=intType),value :: ldc

       
     end subroutine calc_dgemmNN_dot_FP64_GPU

     subroutine calc_dgemmTN_dot_FP64_GPU( &
          icomm, &
          cublas_solver, &
          m, &
          n, &
          k, &
          alpha_, offset_alpha, &
          beta_, offset_beta, &
          A_, offset_A, lda, &
          B_, offset_B, ldb, &
          C_, offset_C, ldc &
          ) bind(c,name="cuda_calc_dgemmTN_dot")
       import :: intType
       import :: devicemem
       import :: cublas_type
       implicit none

       integer,value :: icomm
       type(cublas_type) :: cublas_solver
       integer(kind=intType),value :: m
       integer(kind=intType),value :: n
       integer(kind=intType),value :: k

       type(devicemem) alpha_
       integer(kind=intType),value :: offset_alpha

       type(devicemem) beta_
       integer(kind=intType),value :: offset_beta

       type(devicemem) A_
       integer(kind=intType),value :: offset_A
       integer(kind=intType),value :: lda

       type(devicemem) B_
       integer(kind=intType),value :: offset_B
       integer(kind=intType),value :: ldb

       type(devicemem) C_
       integer(kind=intType),value :: offset_C
       integer(kind=intType),value :: ldc
       
     end subroutine calc_dgemmTN_dot_FP64_GPU

     subroutine neumann_preparation_FP64_GPU( &
          numGPUblock,&
          numGPUthread,&
          n, &
          x, offset_x, &
          y, offset_y, &
          z, offset_z, &
          v, offset_v  &
          )bind(c,name = "cuda_neumann_preparation")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value :: n
       type(devicemem) x
       integer(kind=intType),value :: offset_x
       type(devicemem) y 
       integer(kind=intType),value :: offset_y
       type(devicemem) z
       integer(kind=intType),value :: offset_z
       type(devicemem) v
       integer(kind=intType),value :: offset_v

     end subroutine neumann_preparation_FP64_GPU

     subroutine cuda_get_crsmatrix_neumanndiag_FP64_GPU( &
          n, &
          crsA, &
          crsRow_ptr, &
          crsCol, &
          D, &
          offset_D, &
          Dinv, &
          offset_Dinv &
          )bind(c,name = "cuda_get_crsmatrix_neumanndiag")
       import :: intType
       import :: devicemem
       implicit none
       integer(kind=intType),value :: n
       type(devicemem) crsA
       type(devicemem) crsRow_ptr
       type(devicemem) crsCol
       type(devicemem) D 
       integer(kind=intType),value :: offset_D
       type(devicemem) Dinv
       integer(kind=intType),value :: offset_Dinv

     end subroutine cuda_get_crsmatrix_neumanndiag_FP64_GPU

     subroutine cuda_get_diamatrix_diag_FP64_GPU( &
          n, &
          nnd, &
          diag_offset, &
          diaA, &
          D, &
          offset_D, &
          Dinv, &
          offset_Dinv &
          )bind(c,name = "cuda_get_diamatrix_diag")
       import :: intType
       import :: devicemem
       implicit none
       integer(kind=intType),value :: n
       integer(kind=intType),value :: nnd
       integer(kind=intType),value :: diag_offset
       type(devicemem) diaA
       type(devicemem) D 
       integer(kind=intType),value :: offset_D
       type(devicemem) Dinv
       integer(kind=intType),value :: offset_Dinv

     end subroutine cuda_get_diamatrix_diag_FP64_GPU

     subroutine cuda_get_diamatrix_neumanndiag_FP64_GPU( &
          n, &
          nnd, &
          diag_offset, &
          diaA, &
          offset, &
          iLs, &
          D, &
          offset_D, &
          Dinv, &
          offset_Dinv &
          )bind(c,name = "cuda_get_diamatrix_neumanndiag")
       import :: intType
       import :: devicemem
       implicit none
       integer(kind=intType),value :: n
       integer(kind=intType),value :: nnd
       integer(kind=intType),value :: diag_offset
       type(devicemem) diaA
       type(devicemem) offset
       type(devicemem) iLs
       type(devicemem) D 
       integer(kind=intType),value :: offset_D
       type(devicemem) Dinv
       integer(kind=intType),value :: offset_Dinv

     end subroutine cuda_get_diamatrix_neumanndiag_FP64_GPU

     subroutine cuda_get_crsmatrix_diag_FP64_GPU( &
          n, &
          crsA, &
          crsRow_ptr, &
          crsCol, &
          D, &
          offset_D, &
          Dinv, &
          offset_Dinv &
          )bind(c,name = "cuda_get_crsmatrix_diag")
       import :: intType
       import :: devicemem
       implicit none
       integer(kind=intType),value :: n
       type(devicemem) crsA
       type(devicemem) crsRow_ptr
       type(devicemem) crsCol
       type(devicemem) D 
       integer(kind=intType),value :: offset_D
       type(devicemem) Dinv
       integer(kind=intType),value :: offset_Dinv

     end subroutine cuda_get_crsmatrix_diag_FP64_GPU

     subroutine jacobi_update_FP64_GPU( &
          numGPUblock,&
          numGPUthread,&
          n, &
          x, &
          b, &
          Ax, &
          Dinv &
          )bind(c,name = "cuda_jacobi_update")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value :: n
       type(devicemem) x
       type(devicemem) b 
       type(devicemem) Ax
       type(devicemem) Dinv
     end subroutine jacobi_update_FP64_GPU

     subroutine neumann_update_FP64_GPU( &
          numGPUblock,&
          numGPUthread,&
          n, &
          alpha, &
          x, offset_x, &
          y, offset_y, &
          z, offset_z  &
          )bind(c,name = "cuda_neumann_update")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value :: n
       real*8,value :: alpha
       type(devicemem) x
       integer(kind=intType),value :: offset_x
       type(devicemem) y 
       integer(kind=intType),value :: offset_y
       type(devicemem) z
       integer(kind=intType),value :: offset_z

     end subroutine neumann_update_FP64_GPU

     subroutine vec_scaling_FP64_GPU( &
          numGPUblock,&
          numGPUthread,&
          n, &
          x, offset_x, &
          y, offset_y, &
          z, offset_z  &
          )bind(c,name = "cuda_vec_scaling")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value :: n
       type(devicemem) x
       integer(kind=intType),value :: offset_x
       type(devicemem) y 
       integer(kind=intType),value :: offset_y
       type(devicemem) z
       integer(kind=intType),value :: offset_z

     end subroutine vec_scaling_FP64_GPU

     subroutine zero_init_FP64_GPU(&
          n, &
          x,offset_x &
          ) &
          bind(c,name = "cuda_zero_init")
       import :: intType
       import :: devicemem
       implicit none
       integer(kind=intType),value :: n
       type(devicemem) x
       integer(kind=intType),value :: offset_x
     end subroutine zero_init_FP64_GPU

     subroutine axpy_FP64_GPU(&
          numGPUblock,&
          numGPUthread,&
          n, &
          alpha,&
          beta ,&
          x,offset_x, &
          y,offset_y  &
          ) &
          bind(c,name = "cuda_axpy")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value :: n
       real*8,value :: alpha
       real*8,value :: beta
       type(devicemem) x
       integer(kind=intType),value :: offset_x
       type(devicemem) y 
       integer(kind=intType),value :: offset_y
     end subroutine axpy_FP64_GPU

     subroutine axpy_vec3_FP64_GPU(&
          numGPUblock,&
          numGPUthread,&
          n, &
          alpha,&
          beta ,&
          gamma,&
          x,offset_x, &
          y,offset_y,  &
          z,offset_z  &
          ) &
          bind(c,name = "cuda_axpy_vec3")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value :: n
       real*8,value :: alpha
       real*8,value :: beta
       real*8,value :: gamma

       type(devicemem) x
       integer(kind=intType),value :: offset_x
       type(devicemem) y 
       integer(kind=intType),value :: offset_y
       type(devicemem) z
       integer(kind=intType),value :: offset_z

     end subroutine axpy_vec3_FP64_GPU

     subroutine axpy_vec4_FP64_GPU(&
          numGPUblock,&
          numGPUthread,&
          n, &
          alpha,&
          beta ,&
          gamma,&
          delta,&
          x,offset_x, &
          y,offset_y, &
          z,offset_z, &
          w,offset_w  &
          ) &
          bind(c,name = "cuda_axpy_vec4")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value :: n
       real*8,value :: alpha
       real*8,value :: beta
       real*8,value :: gamma
       real*8,value :: delta

       type(devicemem) x
       integer(kind=intType),value :: offset_x
       type(devicemem) y 
       integer(kind=intType),value :: offset_y
       type(devicemem) z
       integer(kind=intType),value :: offset_z
       type(devicemem) w
       integer(kind=intType),value :: offset_w

     end subroutine axpy_vec4_FP64_GPU

     subroutine gemv_FP64_GPU(&
          numGPUblock,&
          numGPUthread,&
          x,offset_x,ldx, &
          y,offset_y,ldy,  &
          alpha,offset_alpha, &
          n, &
          kmax &
          ) &
          bind(c,name = "cuda_gemv")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       type(devicemem) x
       integer(kind=intType),value :: offset_x
       integer(kind=intType),value :: ldx

       type(devicemem) y 
       integer(kind=intType),value :: offset_y
       integer(kind=intType),value :: ldy

       type(devicemem) alpha 
       integer(kind=intType),value :: offset_alpha

       integer(kind=intType),value :: n
       integer(kind=intType),value :: kmax

     end subroutine gemv_FP64_GPU     

     subroutine calc_dot_FP64_GPU( &
          icomm, &
          numGPUblock,&
          numGPUthread,&
          n, &
          val, offset_val,&
          x,offset_x, &
          y,offset_y  &
          ) &
          bind(c,name = "cuda_dot")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer,value :: icomm
       integer(kind=intType),value :: n
       type(devicemem) val
       integer(kind=intType),value :: offset_val
       type(devicemem) x
       integer(kind=intType),value :: offset_x
       type(devicemem) y 
       integer(kind=intType),value :: offset_y
     end subroutine calc_dot_FP64_GPU

     subroutine calc_dot_local_FP64_GPU( &
          numGPUblock,&
          numGPUthread,&
          n, &
          val, offset_val,&
          x,offset_x, &
          y,offset_y  &
          ) &
          bind(c,name = "cuda_dot_local")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value :: n
       type(devicemem) val
       integer(kind=intType),value :: offset_val
       type(devicemem) x
       integer(kind=intType),value :: offset_x
       type(devicemem) y 
       integer(kind=intType),value :: offset_y
     end subroutine calc_dot_local_FP64_GPU

     subroutine cuda_CGS_Orth_FP64_GPU( &
          icomm, &
          W_  , offset_W, ldw, &
          V_  , offset_V, ldv, &
          hki_, offset_h, &
          n, &
          vrow &
          ) bind(c,name = "cuda_CGS_Orth")
       import :: intType
       import :: devicemem
       implicit none

       integer,value :: icomm

       type(devicemem) W_
       integer(kind=intType),value ::  offset_W
       integer(kind=intType),value ::  ldw

       type(devicemem) V_
       integer(kind=intType),value ::  offset_V
       integer(kind=intType),value ::  ldv

       type(devicemem) hki_
       integer(kind=intType),value ::  offset_h
       integer(kind=intType),value ::  n
       integer(kind=intType),value ::  vrow

     end subroutine cuda_CGS_Orth_FP64_GPU


     subroutine cuda_gemm_dot_local_FP64_GPU( &
          W_  , offset_W, ldw, &
          V_  , offset_V, ldv, &
          wv_, offset_wv, &
          n, &
          wrow, &
          vrow &
          ) bind(c,name = "cuda_gemm_dot_local")
       import :: intType
       import :: devicemem
       implicit none

       type(devicemem) W_
       integer(kind=intType),value ::  offset_W
       integer(kind=intType),value ::  ldw

       type(devicemem) V_
       integer(kind=intType),value ::  offset_V
       integer(kind=intType),value ::  ldv

       type(devicemem) wv_
       integer(kind=intType),value ::  offset_wv
       integer(kind=intType),value ::  n
       integer(kind=intType),value ::  wrow
       integer(kind=intType),value ::  vrow

     end subroutine cuda_gemm_dot_local_FP64_GPU

     subroutine cuda_gemm_dot_FP64_GPU( &
          icomm, &
          W_  , offset_W, ldw, &
          V_  , offset_V, ldv, &
          wv_ , offset_wv, &
          n, &
          wrow, &
          vrow &
          ) bind(c,name = "cuda_gemm_dot")
       import :: intType
       import :: devicemem
       implicit none

       integer,value :: icomm
       type(devicemem) W_
       integer(kind=intType),value ::  offset_W
       integer(kind=intType),value ::  ldw

       type(devicemem) V_
       integer(kind=intType),value ::  offset_V
       integer(kind=intType),value ::  ldv

       type(devicemem) wv_
       integer(kind=intType),value ::  offset_wv
       integer(kind=intType),value ::  n
       integer(kind=intType),value ::  wrow
       integer(kind=intType),value ::  vrow

     end subroutine cuda_gemm_dot_FP64_GPU

     subroutine cuda_gemm_FP64_GPU( &
          x_  , offset_x, ldx, &
          y_  , offset_y, ldy, &
!          z_  , offset_z, ldz, &
          beta_  , offset_beta, ldbeta, &
          n, &
          s,t &
          ) bind(c,name = "cuda_gemm")
       import :: intType
       import :: devicemem
       implicit none

       type(devicemem) x_
       integer(kind=intType),value ::  offset_x
       integer(kind=intType),value ::  ldx

       type(devicemem) y_
       integer(kind=intType),value ::  offset_y
       integer(kind=intType),value ::  ldy

       ! type(devicemem) z_
       ! integer(kind=intType),value ::  offset_z
       ! integer(kind=intType),value ::  ldz

       type(devicemem) beta_
       integer(kind=intType),value ::  offset_beta
       integer(kind=intType),value ::  ldbeta

       integer(kind=intType),value ::  n
       integer(kind=intType),value ::  s
       integer(kind=intType),value ::  t

     end subroutine cuda_gemm_FP64_GPU

     subroutine cuda_sum_vec_FP64_GPU( &
          numGPUblock,&
          numGPUthread,&
          n, &
          m, &
          vecLen,&
          v,  offset_v,  &
          v0, offset_v0, &
          y &
          )bind(c,name = "cuda_sum_vec")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value :: n
       integer(kind=intType),value :: m
       integer(kind=intType),value :: vecLen
       type(devicemem) v
       integer(kind=intType),value :: offset_v
       type(devicemem) v0
       integer(kind=intType),value :: offset_v0
       type(devicemem) y 

     end subroutine cuda_sum_vec_FP64_GPU

     subroutine cuda_sum_vec_cagmres_FP64_GPU( &
          numGPUblock,&
          numGPUthread,&
          n, &
          s, &
          m, &
          mm, &
          Q ,  offset_Q , ldq,&
          V ,  offset_V , ldv,&
          Qy,  offset_Qy, ldqy,&
          y ,  offset_y &
          )bind(c,name = "cuda_sum_vec_cagmres")
       import :: intType
       import :: devicemem
       implicit none
       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value :: n
       integer(kind=intType),value :: s
       integer(kind=intType),value :: m
       integer(kind=intType),value :: mm
       type(devicemem) Q
       integer(kind=intType),value :: offset_Q
       integer(kind=intType),value :: ldq
       type(devicemem) V
       integer(kind=intType),value :: offset_V
       integer(kind=intType),value :: ldv
       type(devicemem) Qy
       integer(kind=intType),value :: offset_Qy
       integer(kind=intType),value :: ldqy
       type(devicemem) y
       integer(kind=intType),value :: offset_y

     end subroutine cuda_sum_vec_cagmres_FP64_GPU

     subroutine cuda_set_scala_ddm_FP64_GPU( &
          numGPUblock, &
          numGPUthread, &
          nx, &
          ny, &
          nz, &
          stmx, &
          stmy, &
          stmz, &
          mx, &
          mxy, &
          alpha,  offset_alpha, &
          x    ,  offset_x &
          )bind(c,name = "cuda_set_scala_ddm")
       import :: intType
       import :: devicemem
       implicit none

       integer,value :: numGPUblock
       integer,value :: numGPUthread
       integer(kind=intType),value ::  nx
       integer(kind=intType),value ::  ny
       integer(kind=intType),value ::  nz
       integer(kind=intType),value ::  stmx
       integer(kind=intType),value ::  stmy
       integer(kind=intType),value ::  stmz
       integer(kind=intType),value ::  mx
       integer(kind=intType),value ::  mxy
       type(devicemem) alpha
       integer(kind=intType),value ::  offset_alpha
       type(devicemem) x
       integer(kind=intType),value ::  offset_x
     end subroutine cuda_set_scala_ddm_FP64_GPU

     subroutine cuda_xcpy_ddm_FP64_GPU( &
          nx, &
          ny, &
          nz, &
          stmx, &
          stmy, &
          stmz, &
          mx, &
          mxy, &
          x    ,  offset_x, &
          y    ,  offset_y  &
          )bind(c,name = "cuda_xcpy_ddm")
       import :: intType
       import :: devicemem
       implicit none

#if 0
       integer(kind=intType),value ::  nx
       integer(kind=intType),value ::  ny
       integer(kind=intType),value ::  nz
       integer(kind=intType),value ::  stmx
       integer(kind=intType),value ::  stmy
       integer(kind=intType),value ::  stmz
       integer(kind=intType),value ::  mx
       integer(kind=intType),value ::  mxy
#else
       integer,value ::  nx
       integer,value ::  ny
       integer,value ::  nz
       integer,value ::  stmx
       integer,value ::  stmy
       integer,value ::  stmz
       integer,value ::  mx
       integer,value ::  mxy
#endif
       type(devicemem) x
       integer(kind=intType),value ::  offset_x
       type(devicemem) y
       integer(kind=intType),value ::  offset_y

     end subroutine cuda_xcpy_ddm_FP64_GPU

  end interface

end module cu_Vec
