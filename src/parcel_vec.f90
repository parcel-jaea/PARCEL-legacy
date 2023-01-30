! #define CUDA_HANDMADE 1

module PARCEL_VEC_kernel
  use krylov
  use vec_type
  use vec_kernel_CPU
  
  use cu_common_wrapper

contains
  subroutine parcel_sync(method)
    use krylov
    implicit none
    type(KSP) :: method
    if(method%flag%machineflag.eq.PARCEL_CPU_FLAG_ON)then
    endif
    if(method%flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
       call cudaDeviceSynchronize()
    endif

  end subroutine parcel_sync

  subroutine init_vec(method,vec,n)
    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: vec
    integer(kind=intType) n
    integer(kind=intType) nvec
    integer(kind=intType) vecLen
    type(STATUS_FLAG) :: flag
    integer(kind=intType),parameter :: offset_zero=0

    ! integer CPUflag
    ! integer GPUflag
    ! integer dataflag

    ! CPU,GPU両方有効な場合 
    ! if(flag%machineflag.eq.PARCEL_ALLARCH_FLAG_ON)then
    !    vec%CPUflag  = PARCEL_CPU_FLAG_ON
    !    vec%GPUflag  = PARCEL_GPU_FLAG_ON
    ! endif
    flag = method%flag
    vec%CPUflag  = PARCEL_CPU_FLAG_ON
    vec%GPUflag  = PARCEL_GPU_FLAG_OFF
    if(flag%machineflag.eq.PARCEL_CPU_FLAG_ON)then
       vec%CPUflag  = PARCEL_CPU_FLAG_ON
       vec%GPUflag  = PARCEL_GPU_FLAG_OFF
    endif
    if(flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
       vec%CPUflag  = PARCEL_CPU_FLAG_OFF
       vec%GPUflag  = PARCEL_GPU_FLAG_ON
    endif
    vec%dataflag = flag%dataflag

    vec%n = n
    ! CPU側の設定
    if(vec%CPUflag.eq.PARCEL_CPU_FLAG_ON)then
       if(vec%dataflag.eq.KSP_DOUBLE_PRECISION)then
          if(allocated(vec%vec_FP64))deallocate(vec%vec_FP64)
          allocate(vec%vec_FP64(n))
       endif
       if(vec%dataflag.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
          ! double - double演算用の初期化を追加
       endif
    endif

    ! GPU側の設定
    ! CUDA handmade
    if(vec%GPUflag.eq.PARCEL_GPU_FLAG_ON)then
       if(vec%dataflag.eq.KSP_DOUBLE_PRECISION)then
          call cudaMalloc_fp64(vec%C_vec_FP64,n)
          ! 記述子のサイズはprmで制限する
          ! call cusparseCreateDnVec_fp64(vec%C_vec_FP64,n)
#if CUDA_HANDMADE

#else
          ! write(*,*)method%prm%n,n
          if(method%prm%n<n)then
             nvec   = method%prm%n
             vecLen = method%prm%n 
             if(method%prm%MatFormat.eq.MATRIX_FORMAT_CRS)then
                vecLen = vecLen + method%crsMat%nColL + method%crsMat%nColU
             endif
!             call cusparseCreateDnVec_FP64(vec%C_vec_FP64,method%prm%n,n)
             call cusparseCreateDnVec_FP64(vec%C_vec_FP64,nvec,vecLen,offset_zero)
          else
             call cusparseCreateDnVec_FP64(vec%C_vec_FP64,n,n,offset_zero)
          endif
#endif          
       end if
    endif
    
    vec%offset = 0
  end subroutine init_vec

  ! ------------------------
  ! 1要素に対する操作
  subroutine set_val_FP64(method,vec,val,i)
    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: vec
    real*8 val
    real*8 val_(1)
    integer(kind=intType) i
    integer(kind=intType) parcel_int_1,parcel_int_2
    integer(kind=intType),parameter :: offset_zero=0

    if(vec%CPUflag.eq.PARCEL_CPU_FLAG_ON)then
       vec%vec_FP64(i) = val
    endif

    if(vec%GPUflag.eq.PARCEL_GPU_FLAG_ON)then
       val_(1)=val
       parcel_int_1 = 1
       parcel_int_2 = i-1
       call cudaMemcpy_h2d_FP64(val_,vec%C_vec_FP64,parcel_int_1,offset_zero,parcel_int_2)
    endif
        
  end subroutine set_val_FP64

  subroutine get_val_FP64(method,vec,val,i)
    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: vec
    real*8 val
    real*8 val_(1)
    integer(kind=intType) i
    integer(kind=intType) parcel_int_1,parcel_int_2
    integer(kind=intType),parameter :: offset_zero=0

    if(vec%CPUflag.eq.PARCEL_CPU_FLAG_ON)then
       val = vec%vec_FP64(i)
    endif

    if(vec%GPUflag.eq.PARCEL_GPU_FLAG_ON)then
       call parcel_sync(method)
       parcel_int_1 = 1
       parcel_int_2 = i-1
       call cudaMemcpy_d2h_FP64(vec%C_vec_FP64,val_, parcel_int_1,parcel_int_2,offset_zero)
       call parcel_sync(method)
       val = val_(1)
    endif

  end subroutine get_val_FP64

  subroutine add_val_FP64(method,vec,val,i)
    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: vec
    real*8 val
    integer(kind=intType) i
    integer(kind=intType),parameter :: offset_zero=0

    if(vec%CPUflag.eq.PARCEL_CPU_FLAG_ON)then
       vec%vec_FP64(i+vec%offset) = vec%vec_FP64(i+vec%offset) + val
    endif

    if(vec%GPUflag.eq.PARCEL_GPU_FLAG_ON)then
!       call cuda_add_val_FP64(1,1,val,vec%C_vec_FP64,vec%offset) 
       call cuda_add_val_FP64(1,1,val,vec%C_vec_FP64,(i-1)+vec%offset) 
    endif

  end subroutine add_val_FP64

  subroutine sub_val_FP64(method,vec,val,i)
    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: vec
    real*8 val
    integer(kind=intType) i

    if(vec%CPUflag.eq.PARCEL_CPU_FLAG_ON)then
       vec%vec_FP64(i+vec%offset) = vec%vec_FP64(i+vec%offset) - val
    endif

    if(vec%GPUflag.eq.PARCEL_GPU_FLAG_ON)then
       call cuda_sub_val_FP64(1,1,val,vec%C_vec_FP64,(i-1)+vec%offset) 
    endif

  end subroutine sub_val_FP64

  subroutine mul_val_FP64(method,vec,val,i)
    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: vec
    real*8 val
    integer(kind=intType) i

    if(vec%CPUflag.eq.PARCEL_CPU_FLAG_ON)then
       vec%vec_FP64(i+vec%offset) = vec%vec_FP64(i+vec%offset) * val
    endif

    if(vec%GPUflag.eq.PARCEL_GPU_FLAG_ON)then
       call cuda_mul_val_FP64(1,1,val,vec%C_vec_FP64,(i-1)+vec%offset) 
    endif

  end subroutine mul_val_FP64

  subroutine div_val_FP64(method,vec,val,i)
    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: vec
    real*8 val
    integer(kind=intType) i

    if(vec%CPUflag.eq.PARCEL_CPU_FLAG_ON)then
       vec%vec_FP64(i+vec%offset) = vec%vec_FP64(i+vec%offset) / val
    endif

    if(vec%GPUflag.eq.PARCEL_GPU_FLAG_ON)then
       call cuda_div_val_FP64(1,1,val,vec%C_vec_FP64,(i-1)+vec%offset) 
    endif

  end subroutine div_val_FP64

  subroutine sqrt_val(method,vec,i)
    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: vec
    integer(kind=intType) i

    if(vec%CPUflag.eq.PARCEL_CPU_FLAG_ON)then
       vec%vec_FP64(i+vec%offset) = sqrt(vec%vec_FP64(i+vec%offset))
    endif

    if(vec%GPUflag.eq.PARCEL_GPU_FLAG_ON)then
!       call cuda_sqrt_val_FP64(1,1,vec%C_vec_FP64,vec%offset) 
       call cuda_sqrt_val_FP64(1,1,vec%C_vec_FP64,(i-1)+vec%offset) 
    endif

  end subroutine sqrt_val
  ! ------------------------

  subroutine set_vec_val(method,vec,val,i,n)
    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: vec
    type(PARCEL_VEC) :: val
    integer(kind=intType),parameter :: offset_zero=0

    integer(kind=intType) i
    integer(kind=intType) j
    integer(kind=intType) n
    integer(kind=intType) parcel_int
    real*8 tmp_(1)
    if(vec%CPUflag.eq.PARCEL_CPU_FLAG_ON)then
       select case(vec%dataflag)
       case(KSP_DOUBLE_PRECISION)
          !$OMP PARALLEL DO PRIVATE(j)
          do j = 1,n
             vec%vec_FP64(j+vec%offset) = val%vec_FP64(i+val%offset)
          enddo
       case(KSP_DOUBLE_DOUBLE_PRECISION)
          !$OMP PARALLEL DO PRIVATE(j)
          do j = 1,n
             vec%vec_hi(j+vec%offset) = val%vec_hi(i+val%offset)
             vec%vec_lo(j+vec%offset) = val%vec_lo(i+val%offset)
          enddo
       end select
    endif

    ! valはCPU側のメモリ
    if(vec%GPUflag.eq.PARCEL_GPU_FLAG_ON)then
       parcel_int = 1
       call cudaMemcpy_d2h_fp64(val%C_vec_FP64,tmp_,parcel_int,offset_zero,offset_zero)
       call set_vec_val_FP64_GPU(&
            CUDA_NUM_BLOCKS, &
            CUDA_NUM_THREADS, &
            n, &
            tmp_(1), &
            vec%C_vec_FP64,vec%offset &
            )
    endif

  end subroutine set_vec_val

  subroutine set_vec_FP64(method,vec,x,n)
    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: vec
    real*8 x(vec%n)
    integer(kind=intType) n
    integer(kind=intType) j
    integer(kind=intType),parameter :: offset_zero=0

    if(vec%CPUflag.eq.PARCEL_CPU_FLAG_ON)then
       !$OMP PARALLEL DO PRIVATE(j)
       do j = 1,n
          vec%vec_FP64(j + vec%offset) = x(j)
       enddo
    endif
    if(vec%GPUflag.eq.PARCEL_GPU_FLAG_ON)then
       call cudaMemcpy_h2d_fp64(x,vec%C_vec_FP64,n,offset_zero,vec%offset)
    endif

  end subroutine set_vec_FP64

  subroutine get_vec_FP64(method,vec,x,n)
    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: vec
    real*8 x(vec%n)
    integer(kind=intType) n
    integer(kind=intType) j
    integer(kind=intType),parameter :: offset_zero=0

    if(vec%CPUflag.eq.PARCEL_CPU_FLAG_ON)then
       !$OMP PARALLEL DO PRIVATE(j)
       do j = 1,n
          x(j) = vec%vec_FP64(j + vec%offset)
       enddo
    endif

    if(vec%GPUflag.eq.PARCEL_GPU_FLAG_ON)then
!       call cudaMemcpy_d2h_fp64(vec%C_vec_FP64,x,n,0,0) 
       call cudaMemcpy_d2h_fp64(vec%C_vec_FP64,x,n,vec%offset,offset_zero) 
    endif

  end subroutine get_vec_FP64

  subroutine free_vec(method,vec)
    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: vec
    integer CPUflag
    integer GPUflag
    integer MatFormat 

    MatFormat = method%prm%MatFormat 

    CPUflag = vec%CPUflag
    GPUflag = vec%GPUflag

    if(CPUflag.eq.PARCEL_CPU_FLAG_ON)then
       if(allocated(vec%vec_FP64))deallocate(vec%vec_FP64)
       vec%CPUflag = PARCEL_CPU_FLAG_OFF
    endif

    if(GPUflag.eq.PARCEL_GPU_FLAG_ON)then
#if CUDA_HANDMADE
       
#else 
       call cusparseDestroyDnVec(vec%C_vec_FP64)
#endif       
       call cudaFree(vec%C_vec_FP64)
       vec%GPUflag = PARCEL_GPU_FLAG_OFF
    endif
    
  end subroutine free_vec
  
  subroutine calc_dot(method,comm,x,y,dotp,n)
    use cu_Vec
    use vec_kernel_CPU
    implicit none
    type(KSP) :: method
    type(STATUS_FLAG) :: flag
    integer comm
    type(PARCEL_VEC) :: x
    type(PARCEL_VEC) :: y
    type(PARCEL_VEC) :: dotp ! 1要素の配列
    integer(kind=intType) n

    flag = method%flag
    if(flag%machineflag.eq.PARCEL_CPU_FLAG_ON)then
       call calc_dot_FP64_CPU( &
            n,comm,&
            x%vec_FP64(1+x%offset),&
            y%vec_FP64(1+y%offset),&
            dotp%vec_FP64(1+dotp%offset))
    endif

    if(flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
!       write(*,*)"calc_dot_FP64_GPU : n = ",n
       call calc_dot_FP64_GPU( &
          ! numGPUblock,&
          ! numGPUthread,&
            comm, &
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            n, &
            dotp%C_vec_FP64,dotp%offset, &
            x%C_vec_FP64,x%offset, &
            y%C_vec_FP64,y%offset  &
            ) 
    endif

  end subroutine calc_dot

  subroutine calc_dot_local(method,x,y,dotp,n)
    use cu_Vec
    use vec_kernel_CPU
    implicit none
    type(KSP) :: method
    type(STATUS_FLAG) :: flag
    type(PARCEL_VEC) :: x
    type(PARCEL_VEC) :: y
    type(PARCEL_VEC) :: dotp ! 1要素の配列
    integer(kind=intType) n

    flag = method%flag
    if(flag%machineflag.eq.PARCEL_CPU_FLAG_ON)then
       call calc_dot_local_FP64_CPU( &
            n, &
            x%vec_FP64(1+x%offset),&
            y%vec_FP64(1+y%offset),&
            dotp%vec_FP64(1))
    endif

    if(flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
       call calc_dot_local_FP64_GPU( &
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            n, &
            dotp%C_vec_FP64, dotp%offset, &
            x%C_vec_FP64,x%offset, &
            y%C_vec_FP64,y%offset  &
            ) 
    endif

  end subroutine calc_dot_local

  ! axpy
  subroutine axpy( &
       method, &
       alpha,beta,&
       x,y,n ,&
       x_offset_,y_offset_)
    use cu_Vec
    implicit none
    type(KSP) :: method
    type(STATUS_FLAG) :: flag
    type(PARCEL_VEC) :: alpha
    type(PARCEL_VEC) :: beta
    type(PARCEL_VEC) :: x
    type(PARCEL_VEC) :: y
    integer(kind=intType) n
    real*8 alpha_(1)
    real*8 beta_(1)

    integer(kind=intType) x_offset
    integer(kind=intType) y_offset
    integer(kind=intType),optional :: x_offset_
    integer(kind=intType),optional :: y_offset_
    integer(kind=intType) parcel_int
    integer(kind=intType),parameter :: offset_zero=0

    x_offset = x%offset
    y_offset = y%offset

    if(  present(x_offset_)  )then
       x_offset = x_offset_
    endif
    if(  present(y_offset_)  )then
       y_offset = y_offset_
    endif

    flag = method%flag
    if(flag%machineflag.eq.PARCEL_CPU_FLAG_ON)then
       call axpy_FP64_CPU( &
            n, &
            alpha%vec_FP64(1),beta%vec_FP64(1),&
            x%vec_FP64(1+x_offset),&
            y%vec_FP64(1+y_offset) &
            )
    endif

    if(flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
       parcel_int = 1
       call cudaMemcpy_d2h_fp64(alpha%C_vec_FP64,alpha_, parcel_int,offset_zero,offset_zero) 
       call cudaMemcpy_d2h_fp64(beta%C_vec_FP64 ,beta_ , parcel_int,offset_zero,offset_zero) 
       call axpy_FP64_GPU( &
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            n, &
            alpha_(1), &
            beta_(1) ,&
            x%C_vec_FP64,x_offset, &
            y%C_vec_FP64,y_offset)
    endif

  end subroutine axpy

  ! axpy
  ! z(:)=alpha*x(:)+beta*y(:)+gamma*z(:)
  subroutine axpy_vec3(&
       method,&
       alpha,beta,gamma,&
       x,y,z,&
       n,&
       x_offset_,y_offset_,z_offset_)

    use cu_Vec
    implicit none
    type(KSP) :: method
    type(STATUS_FLAG) :: flag
    type(PARCEL_VEC) :: alpha
    type(PARCEL_VEC) :: beta
    type(PARCEL_VEC) :: gamma
    type(PARCEL_VEC) :: x
    type(PARCEL_VEC) :: y
    type(PARCEL_VEC) :: z
    integer(kind=intType) n
    real*8 alpha_(1)
    real*8 beta_(1)
    real*8 gamma_(1)

    integer(kind=intType) x_offset
    integer(kind=intType) y_offset
    integer(kind=intType) z_offset

    integer(kind=intType),optional :: x_offset_
    integer(kind=intType),optional :: y_offset_
    integer(kind=intType),optional :: z_offset_
    integer(kind=intType) parcel_int0,parcel_int1

#if 0
    real*8 x_(n)
    real*8 y_(n)
    real*8 z_(n)
#endif

    flag = method%flag

    x_offset = x%offset
    y_offset = y%offset
    z_offset = z%offset

    if(  present(x_offset_)  )then
       x_offset = x_offset_
    endif
    if(  present(y_offset_)  )then
       y_offset = y_offset_
    endif
    if(  present(z_offset_)  )then
       z_offset = z_offset_
    endif

    if(flag%machineflag.eq.PARCEL_CPU_FLAG_ON)then
       call axpy_vec3_FP64_CPU( &
            n, &
            alpha%vec_FP64(1),beta%vec_FP64(1),gamma%vec_FP64(1),&
            x%vec_FP64(1+x_offset), &
            y%vec_FP64(1+y_offset), &
            z%vec_FP64(1+z_offset)  &
            )
    endif

    if(flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
       parcel_int0 = 0
       parcel_int1 = 1
       ! alpha, &betaはvalue  --> スカラ値がCPU側にある
       call cudaMemcpy_d2h_fp64(alpha%C_vec_FP64, &
            alpha_, &
            parcel_int1, &
            parcel_int0, &
            parcel_int0) 
       call cudaMemcpy_d2h_fp64(beta%C_vec_FP64 , &
            beta_ , &
            parcel_int1, &
            parcel_int0, &
            parcel_int0) 
       call cudaMemcpy_d2h_fp64(gamma%C_vec_FP64, &
            gamma_, &
            parcel_int1, &
            parcel_int0, &
            parcel_int0) 
#if 0
       call parcel_sync(method)
       call get_vec_FP64(method,x,x_,n)
       call get_vec_FP64(method,y,y_,n)
       call get_vec_FP64(method,z,z_,n)
       call parcel_sync(method)
       call axpy_vec3_FP64_CPU( &
            n, &
            alpha_(1) ,&
            beta_(1)  ,&
            gamma_(1) ,&
            x_, &
            y_, &
            z_  &
            )
       call parcel_sync(method)
       call set_vec_FP64(method,x,x_,n)
       call set_vec_FP64(method,y,y_,n)
       call set_vec_FP64(method,z,z_,n)
       call parcel_sync(method)

#else
       call axpy_vec3_FP64_GPU( &
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            n, &
            alpha_(1), &
            beta_(1) ,&
            gamma_(1) ,&
            x%C_vec_FP64,x_offset, &
            y%C_vec_FP64,y_offset, &
            z%C_vec_FP64,z_offset  &
            )
#endif
    endif

  end subroutine axpy_vec3

  ! axpy
  ! w(:)=alpha*x(:)+beta*y(:)+gamma*z(:)+delta*w(:)
  subroutine axpy_vec4(&
       method,&
       alpha,beta,gamma,delta,&
       x,y,z,w,&
       n,&
       x_offset_,y_offset_,z_offset_,w_offset_)

    use cu_Vec
    implicit none
    type(KSP) :: method
    type(STATUS_FLAG) :: flag
    type(PARCEL_VEC) :: alpha
    type(PARCEL_VEC) :: beta
    type(PARCEL_VEC) :: gamma
    type(PARCEL_VEC) :: delta
    type(PARCEL_VEC) :: x
    type(PARCEL_VEC) :: y
    type(PARCEL_VEC) :: z
    type(PARCEL_VEC) :: w
    integer(kind=intType) n
    real*8 alpha_(1)
    real*8 beta_(1)
    real*8 gamma_(1)
    real*8 delta_(1)

    integer(kind=intType) x_offset
    integer(kind=intType) y_offset
    integer(kind=intType) z_offset
    integer(kind=intType) w_offset

    integer(kind=intType),optional :: x_offset_
    integer(kind=intType),optional :: y_offset_
    integer(kind=intType),optional :: z_offset_
    integer(kind=intType),optional :: w_offset_
    integer(kind=intType) parcel_int0,parcel_int1

#if 0
    real*8 x_(n)
    real*8 y_(n)
    real*8 z_(n)
    real*8 w_(n)
#endif
 
    flag = method%flag

    x_offset = x%offset
    y_offset = y%offset
    z_offset = z%offset
    w_offset = w%offset

    if(  present(x_offset_)  )then
       x_offset = x_offset_
    endif
    if(  present(y_offset_)  )then
       y_offset = y_offset_
    endif
    if(  present(z_offset_)  )then
       z_offset = z_offset_
    endif
    if(  present(w_offset_)  )then
       w_offset = w_offset_
    endif

    if(flag%machineflag.eq.PARCEL_CPU_FLAG_ON)then
       call axpy_vec4_FP64_CPU( &
            n, &
            alpha%vec_FP64(1), &
            beta%vec_FP64(1) , &
            gamma%vec_FP64(1), &
            delta%vec_FP64(1), &
            x%vec_FP64(1+x_offset), &
            y%vec_FP64(1+y_offset), &
            z%vec_FP64(1+z_offset), &
            w%vec_FP64(1+w_offset)  &
            )
    endif

    if(flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
       parcel_int0 = 0
       parcel_int1 = 1
       ! alpha,betaはvalue  --> スカラ値がCPU側にある
       call cudaMemcpy_d2h_fp64(alpha%C_vec_FP64, &
            alpha_, &
            parcel_int1, &
            parcel_int0, &
            parcel_int0) 
       call cudaMemcpy_d2h_fp64(beta%C_vec_FP64, &
            beta_  , &
            parcel_int1, &
            parcel_int0, &
            parcel_int0) 
       call cudaMemcpy_d2h_fp64(gamma%C_vec_FP64, &
            gamma_, &
            parcel_int1, &
            parcel_int0, &
            parcel_int0) 
       call cudaMemcpy_d2h_fp64(delta%C_vec_FP64, &
            delta_, &
            parcel_int1, &
            parcel_int0, &
            parcel_int0) 

#if 0
       call get_vec_FP64(method,x,x_,n)
       call get_vec_FP64(method,y,y_,n)
       call get_vec_FP64(method,z,z_,n)
       call get_vec_FP64(method,w,w_,n)
       call axpy_vec4_FP64_CPU( &
            n, &
            alpha_(1), &
            beta_(1) ,&
            gamma_(1) ,&
            delta_(1) ,&
            x_, &
            y_, &
            z_, &
            w_ &
            )
       call set_vec_FP64(method,x,x_,n)
       call set_vec_FP64(method,y,y_,n)
       call set_vec_FP64(method,z,z_,n)
       call set_vec_FP64(method,w,w_,n)
#else
       call axpy_vec4_FP64_GPU( &
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            n, &
            alpha_(1), &
            beta_(1) ,&
            gamma_(1) ,&
            delta_(1) ,&
            x%C_vec_FP64,x_offset, &
            y%C_vec_FP64,y_offset, &
            z%C_vec_FP64,z_offset, &
            w%C_vec_FP64,w_offset  &
            )
#endif
    endif

  end subroutine axpy_vec4

  ! y(:)=x(:)
  subroutine xcpy(method,x,y,n)
    use PARCEL_integer
    implicit none
    type(KSP) :: method
    type(STATUS_FLAG) :: flag
    type(PARCEL_VEC) :: x
    type(PARCEL_VEC) :: y
    integer(kind=intType) n
    type(PARCEL_VEC) :: alpha
    type(PARCEL_VEC) :: beta
    integer(kind=intType) parcel_int

    ! alpha=1.0d0
    ! beta=0.0d0

    flag = method%flag
    parcel_int = 1
    call init_vec(method,alpha,parcel_int)
    call init_vec(method,beta ,parcel_int)
    call set_val_FP64(method,alpha,1.0d0,parcel_int)
    call set_val_FP64(method,beta ,0.0d0,parcel_int)

    call parcel_sync(method)
    call axpy(method,alpha,beta,x,y,n)
    call parcel_sync(method)

    call free_vec(method,alpha)
    call free_vec(method,beta)

  end subroutine xcpy

#if 0
  subroutine calc_gemm_dot(&
       method,n,s,t,icomm,&
       x, ldx,&
       y, ldy,&
       dotp, lddot,&
       nvecBlock)
    use vec_kernel,only : calc_gemm_dot_CPU
    use cu_Vec,only : calc_dgemmTN_dot_FP64_GPU
    use cu_common_wrapper,only : cudaMemcpy_h2d_FP64
    use PARCEL_integer
    implicit none
    type(KSP) :: method
    integer(kind=intType) n
    integer icomm
    integer(kind=intType) s,t
    type(PARCEL_VEC) :: x ! (n*s)
    integer(kind=intType) ldx
    type(PARCEL_VEC) :: y ! (n*t)
    integer(kind=intType) ldy
    type(PARCEL_VEC) :: dotp ! (s*t)
    integer(kind=intType) lddot
    integer nvecBlock
    type(STATUS_FLAG) :: flag
    integer(kind=intType) lda,ldb,ldc
    integer(kind=intType),parameter :: offset_zero=0
    real*8  alpha(1)
    real*8  beta(1)

    flag = method%flag
    if(flag%machineflag.eq.PARCEL_CPU_FLAG_ON)then
       call calc_gemm_dot_CPU(&
            n,s,t,icomm,&
            x%vec_FP64(1+x%offset), ldx,&
            y%vec_FP64(1+y%offset), ldy,&
            dotp%vec_FP64(1+dotp%offset),lddot, &
            nvecBlock)
    endif

    if(flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
       lda = ldx
       ldb = ldy
       ldc = lddot
       ! C = alpha*A[k*m]*B[k*n] + beta*C[m*n]
       alpha = 1.0d0
       beta  = 0.0d0

       call cudaMemcpy_h2d_FP64(alpha,method%cublas_solver%alpha,1,offset_zero,offset_zero)
       call cudaMemcpy_h2d_FP64(beta ,method%cublas_solver%beta ,1,offset_zero,offset_zero)
       call calc_dgemmTN_dot_FP64_GPU(&
            icomm, &
            method%cublas_solver, &
            s, &
            t, &
            n, &
            method%cublas_solver%alpha, offset_zero, &
            method%cublas_solver%beta , offset_zero,  &
            x%C_vec_FP64, x%offset, lda, &
            y%C_vec_FP64, y%offset, ldb, &
            dotp%C_vec_FP64, dotp%offset, ldc )
    endif

  end subroutine calc_gemm_dot
#endif

  subroutine calc_gemm_dot_local(&
       method,n,s,t,&
       x    , ldx,&
       y    , ldy,&
       xydot, ldxydot,&
       nvecBlock)

    use krylov,only : KSP
    use vec_type,only : PARCEL_VEC
    use PARCEL_integer,only : intType
    use KSP_STATUS,only : STATUS_FLAG
    use PARCEL_define,only : PARCEL_CPU_FLAG_ON,PARCEL_GPU_FLAG_ON
    use cu_Vec,only : cuda_gemm_dot_local_FP64_GPU
    use vec_kernel,only : calc_gemm_dot_local_CPU
    implicit none

    type(KSP) :: method
    integer(kind=intType) n,s,t

    type(PARCEL_VEC) :: x    ! (n*s)
    integer(kind=intType) ldx

    type(PARCEL_VEC) :: y ! (n*t)
    integer(kind=intType) ldy

    type(PARCEL_VEC) :: xydot ! (n*t)
    integer(kind=intType) ldxydot

    integer nvecBlock

    type(STATUS_FLAG) :: flag
    integer(kind=intType),parameter :: offset_zero=0

    flag = method%flag
    if(flag%machineflag.eq.PARCEL_CPU_FLAG_ON)then
       call calc_gemm_dot_local_CPU(&
            n,s,t,&
            x%vec_FP64(1+x%offset), ldx,&
            y%vec_FP64(1+y%offset), ldy,&
            xydot%vec_FP64(1+xydot%offset), ldxydot, &
            nvecBlock)
    endif

    if(flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
       call cuda_gemm_dot_local_FP64_GPU( &
          x%C_vec_FP64    , x%offset, ldx, &
          y%C_vec_FP64    , y%offset, ldy, &
          xydot%C_vec_FP64, xydot%offset , &
          n, &
          s, &
          t &
          ) 

    endif

  end subroutine calc_gemm_dot_local

  subroutine calc_gemm_dot(&
       method, &
       icomm,&
       n,s,t,&
       x    , ldx,&
       y    , ldy,&
       xydot, ldxydot,&
       nvecBlock)

    use krylov,only : KSP
    use vec_type,only : PARCEL_VEC
    use PARCEL_integer,only : intType
    use KSP_STATUS,only : STATUS_FLAG
    use PARCEL_define,only : PARCEL_CPU_FLAG_ON,PARCEL_GPU_FLAG_ON
    use cu_Vec,only : cuda_gemm_dot_FP64_GPU
    use vec_kernel,only : calc_gemm_dot_CPU

    implicit none

    type(KSP) :: method
    integer icomm
    integer(kind=intType) n,s,t

    type(PARCEL_VEC) :: x    ! (n*s)
    integer(kind=intType) ldx

    type(PARCEL_VEC) :: y ! (n*t)
    integer(kind=intType) ldy

    type(PARCEL_VEC) :: xydot ! (n*t)
    integer(kind=intType) ldxydot

    integer nvecBlock

    type(STATUS_FLAG) :: flag
    integer(kind=intType),parameter :: offset_zero=0

    flag = method%flag
    if(flag%machineflag.eq.PARCEL_CPU_FLAG_ON)then
       call calc_gemm_dot_CPU(&
            n,s,t,&
            icomm, &
            x%vec_FP64(1+x%offset), ldx,&
            y%vec_FP64(1+y%offset), ldy,&
            xydot%vec_FP64(1+xydot%offset), ldxydot, &
            nvecBlock)
    endif

    if(flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then

! x : Q
! y : V
! s : i
! t : j = 1
       call cuda_gemm_dot_FP64_GPU( &
            icomm,&
            x%C_vec_FP64    , x%offset, ldx, &
            y%C_vec_FP64    , y%offset, ldy, &
            xydot%C_vec_FP64, xydot%offset , &
            n, &
            s, &
            t &
            ) 

    endif

  end subroutine calc_gemm_dot

  subroutine calc_gemv(&
       method, &
       x, ldx, &
       y, ldy, &
       alpha,  &
       n, &
       m, &
       nvecBlock)
    use krylov,only : KSP
    use vec_type,only : PARCEL_VEC
    use PARCEL_integer,only : intType
    use KSP_STATUS,only : STATUS_FLAG
    use PARCEL_define,only : PARCEL_CPU_FLAG_ON,PARCEL_GPU_FLAG_ON
    use c_wrapper,only :  CUDA_NUM_BLOCKS , CUDA_NUM_THREADS
    use vec_kernel,only : gemv
    use cu_Vec,only : gemv_FP64_GPU

    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: x 
    integer(kind=intType) ldx

    type(PARCEL_VEC) :: y
    integer(kind=intType) ldy

    type(PARCEL_VEC) :: alpha ! alpha(m)

    integer(kind=intType) n
    integer(kind=intType) m
    integer nvecBlock

    type(STATUS_FLAG) :: flag
    integer(kind=intType),parameter :: offset_zero=0
#if 1
     real*8,allocatable :: x_(:)
     real*8,allocatable :: y_(:)
     real*8 :: alpha_(m)
#endif

    flag = method%flag
    if(flag%machineflag.eq.PARCEL_CPU_FLAG_ON)then
       call gemv(&
            n,m,&
            x%vec_FP64(1+x%offset), ldx, &
            y%vec_FP64(1+y%offset), ldy, &
            alpha%vec_FP64(1+alpha%offset) ,nvecBlock)
    endif

    if(flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
#if 0
       allocate(x_(ldx*m))
       allocate(y_(ldy))
       call get_vec_FP64(method,x,x_,ldx*m)
       call get_vec_FP64(method,y,y_,ldy)
       call get_vec_FP64(method,alpha,alpha_,m)

       call gemv(&
            n,m,&
            x_, ldx, &
            y_, ldy, &
            alpha_ ,nvecBlock)

       call set_vec_FP64(method,x,x_,ldx*m)
       call set_vec_FP64(method,y,y_,ldy)

       deallocate(x_)
       deallocate(y_)
       return 
#else
       call gemv_FP64_GPU(&
            CUDA_NUM_BLOCKS, &
            CUDA_NUM_THREADS, &
            x%C_vec_FP64, x%offset, ldx, &
            y%C_vec_FP64, y%offset, ldy, &
            alpha%C_vec_FP64, alpha%offset, &
            n, &
            m  &
            )
#endif
    endif

  end subroutine calc_gemv

  subroutine calc_gemm(&
       method, &
       x, ldx, &
       y, ldy, &
       alpha, ldalpha, &
       n, &
       s, t,&
       nvecBlock)
    use krylov,only : KSP
    use vec_type,only : PARCEL_VEC
    use PARCEL_integer,only : intType
    use KSP_STATUS,only : STATUS_FLAG
    use PARCEL_define,only : PARCEL_CPU_FLAG_ON,PARCEL_GPU_FLAG_ON
    use c_wrapper,only :  CUDA_NUM_BLOCKS , CUDA_NUM_THREADS
    use vec_kernel,only : gemm
    use cu_Vec,only : cuda_gemm_FP64_GPU

    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: x 
    integer(kind=intType) ldx

    type(PARCEL_VEC) :: y
    integer(kind=intType) ldy

    type(PARCEL_VEC) :: alpha
    integer(kind=intType) ldalpha

    integer(kind=intType) n
    integer(kind=intType) s
    integer(kind=intType) t
    integer nvecBlock

    type(STATUS_FLAG) :: flag
    integer(kind=intType),parameter :: offset_zero=0

    flag = method%flag
    if(flag%machineflag.eq.PARCEL_CPU_FLAG_ON)then
       call gemm( &
            x%vec_FP64(1+x%offset),ldx, &
            y%vec_FP64(1+y%offset),ldy, &
            alpha%vec_FP64(1+alpha%offset),ldalpha, &
            n, s, t, &
            nvecBlock)
    endif

    if(flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
       call cuda_gemm_FP64_GPU( &
            x%C_vec_FP64    , x%offset   , ldx, &
            y%C_vec_FP64    , y%offset   , ldy, &
            alpha%C_vec_FP64 , alpha%offset, ldalpha, &
            n, &
            s,t &
            ) 
    endif

  end subroutine calc_gemm
  
  subroutine set_scala_ddm(method,alpha,x,n)
    use vec_ddm_kernel,only : set_scala_ddm_CPU
    use vec_type,only : PARCEL_VEC
    use PARCEL_integer,only : intType
    use krylov,only : KSP
    use KSP_STATUS,only : &
         STATUS_FLAG,PARCEL_CPU_FLAG_ON,PARCEL_GPU_FLAG_ON
    use PARCEL_integer,only : intType
    use cu_Vec,only : cuda_set_scala_ddm_FP64_GPU
    use c_wrapper,only :  CUDA_NUM_BLOCKS , CUDA_NUM_THREADS
    implicit none

    type(KSP) :: method
    type(PARCEL_VEC) :: x 
    type(PARCEL_VEC) :: alpha
    integer(kind=intType) n

    type(STATUS_FLAG) :: flag
    integer(kind=intType),parameter :: offset_zero=0
    integer(kind=intType) nx_,ny_,nz_
    integer(kind=intType) stmx_,stmy_,stmz_
    integer(kind=intType) mx_,mxy_

    flag = method%flag
    if(flag%machineflag.eq.PARCEL_CPU_FLAG_ON)then
       call set_scala_ddm_CPU( &
            n, &
            alpha%vec_FP64(1+alpha%offset), &
            x%vec_FP64(1+x%offset), &
            method%ddmMat%prm)
    endif

    if(flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
       nx_   = method%ddmMat%prm%nx
       ny_   = method%ddmMat%prm%ny
       nz_   = method%ddmMat%prm%nz
       stmx_ = method%ddmMat%prm%stmx
       stmy_ = method%ddmMat%prm%stmy
       stmz_ = method%ddmMat%prm%stmz
       mx_   = method%ddmMat%prm%mx
       mxy_  = method%ddmMat%prm%mxy

       call cuda_set_scala_ddm_FP64_GPU( &
            CUDA_NUM_BLOCKS, &
            CUDA_NUM_THREADS, &
#if 1
            nx_, &
            ny_, &
            nz_, &
            stmx_, &
            stmy_, &
            stmz_, &
            mx_, &
            mxy_, &
#else
            method%ddmMat%prm%nx, &
            method%ddmMat%prm%ny, &
            method%ddmMat%prm%nz, &
            method%ddmMat%prm%stmx, &
            method%ddmMat%prm%stmy, &
            method%ddmMat%prm%stmz, &
            method%ddmMat%prm%mx, &
            method%ddmMat%prm%mxy, &
#endif
            alpha%C_vec_FP64, alpha%offset, &
            x%C_vec_FP64    , x%offset &
            )
    endif

  end subroutine set_scala_ddm

 
end module PARCEL_VEC_kernel
