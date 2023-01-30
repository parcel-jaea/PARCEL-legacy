! #define CUDA_HANDMADE 1
#define JACOBI_MATFP64 0
! #define ITRATION_JACOBI 2
! #define CUDA_AWARE 1

module PARCEL_KRYLOV_MATRIX
  use vec_type
  use PARCEL_VEC_kernel
  implicit none
contains
  subroutine Matvec_CRS_GPU( &
       method,crsMat,&
       x,offset_x, &
       y,offset_y)  
    use krylov
    use cu_Matrix    
    use cu_comm
    use cu_Vec
    use PARCEL_VEC_kernel
    implicit none
    type(KSP) :: method
    type(C_CRS_Mat) crsMat
    type(devicemem) x,y
    integer(kind=intType) offset_x,offset_y
    integer(kind=intType),parameter :: offset_zero=0
    integer(kind=intType) n,nadd
    integer(kind=intType) nColL_,nColU_
    integer(kind=intType) parcel_int

    call get_nrk_vec_FP64_GPU(&
         method%commVec, &
         x , offset_x, &
         method%C_sbuff,  offset_zero, &
         method%C_buff ,  offset_zero)
    call parcel_sync(method)
    ! nColL_ = crsMat%nColL
    ! nColU_ = crsMat%nColU
    call get_nrk_vecLU_add_FP64_GPU( &
         method%C_buff, offset_zero, &
         method%prm%n, &
         crsMat%nColL, &
         crsMat%nColU, &
         !            nColL_,nColU_,&
         x,offset_x   &
         )

    call parcel_sync(method)
    n = method%prm%n
    nadd = &
         method%prm%n + &
         method%crsMat%nColL + &
         method%crsMat%nColU
    if(offset_x.ne.offset_zero)then
       call cusparseDestroyDnVec(x)
       call cusparseCreateDnVec_FP64(x,n,nadd,offset_x)
    endif
    if(offset_y.ne.offset_zero)then
       call cusparseDestroyDnVec(y)
       call cusparseCreateDnVec_FP64(y,n,nadd,offset_y)
    endif
    call parcel_sync(method)
    call cusparse_SpMV( &
         method%cusparce_solver, &
         x,offset_x, &
         y,offset_y  &
         ) 
    call parcel_sync(method)
    if(offset_x.ne.offset_zero)then
       call cusparseDestroyDnVec(x)
       call cusparseCreateDnVec_FP64(x,n,nadd,offset_zero)
    endif
    if(offset_y.ne.offset_zero)then
       call cusparseDestroyDnVec(y)
       call cusparseCreateDnVec_FP64(y,n,nadd,offset_zero)
    endif
    call parcel_sync(method)

  end subroutine Matvec_CRS_GPU

  subroutine Matvec_DIA_GPU( &
       method,diaMat,&
       x, offset_x,&
       y, offset_y,&
       nocomm)  
    use krylov
    use cu_Matrix    
    use cu_comm
    use PARCEL_VEC_kernel
    implicit none
    type(KSP) :: method
    type(C_DIA_Mat) diaMat
    type(devicemem) x,y
    integer(kind=intType) offset_x,offset_y
    integer(kind=intType),parameter :: offset_zero=0
    integer ,optional :: nocomm
    integer k
    integer(kind=intType) nbuff
    integer(kind=intType) parcel_int

    call parcel_sync(method)
    if(  present(nocomm)  )then
       ! 通信なし
       nbuff = 0
       do k=1,method%commVec%nrk_rank_send
          nbuff = nbuff + method%commVec%nrk_prm_send(k)%ndata
       enddo
       call set_vec_val_FP64_GPU(&
            CUDA_NUM_BLOCKS, &
            CUDA_NUM_THREADS, &
            nbuff, &
            0.0d0, &
            method%C_buff,offset_zero &
            )
    else
       call get_nrk_vec_FP64_GPU(&
            method%commVec, &
            x, offset_x, &
            method%C_sbuff, offset_zero, &
            method%C_buff , offset_zero)
    endif
    call parcel_sync(method)
    parcel_int = method%nBlock
    call MatVec_dia_FP64_GPU( &
         CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
         method%prm%n,&
!         method%prm%nnd,&
         diaMat%nnd,&
         method%prm%istart,&
         parcel_int,&
         x,offset_x, &
         y,offset_y, &
         diaMat%diaA, &
         diaMat%offset , &
         diaMat%iLs, &
         method%C_buff &
         )
    call parcel_sync(method)

  end subroutine Matvec_DIA_GPU

  ! y = Ax
  subroutine calc_Matvec_CRS(method,x,y)
    use krylov
    use comm_kernel
    use cu_Matrix
    
    use cu_comm
    use PARCEL_MATRIX_kernel_CPU
    use PARCEL_VEC_kernel

    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: x
    type(PARCEL_VEC) :: y
    integer iovlflag

!    iovlflag=method%iovlflag
    iovlflag=KSP_SPMV_OVL_NONE
    select case (iovlflag)
    case(KSP_SPMV_OVL_NONE)
       if(method%flag%machineflag.eq.PARCEL_CPU_FLAG_ON)then

          call get_nrk_vec(&
               method%commVec, &
               x%vec_FP64(1+x%offset), &
               method%buff)
          call get_nrk_vecLU(&
               method%commVec, &
               method%buff, &
               method%crsMat%nColL, &
               method%crsMat%nColU, &
               method%rbuffL, &
               method%rbuffU&
               )
          call MatVec_crs_FP64_CPU(&
               method%prm, &
               method%crsMat, &
               method%rbuffL, &
               x%vec_FP64(1+x%offset), &
               method%rbuffU, &
               y%vec_FP64(1+y%offset),&
               method%nBlock)

       endif
       
       if(method%flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
          call Matvec_CRS_GPU(method, &
               method%C_crsMat, &
               x%C_vec_FP64, x%offset, &
               y%C_vec_FP64, y%offset  &
               )
       endif
              
    case(KSP_SPMV_OVL_CORE)
       ! call MatVec_crs_ovl1( &
       !      method%prm, method%crsMat, &
       !      method%rbuffL, x%vec, method%rbuffU, &
       !      y, &
       !      method%commVec,&
       !      method%nBlock)
    case(KSP_SPMV_OVL_RANK)
       ! call MatVec_crs_ovl2(&
       !      method%prm, method%crsMat, &
       !      x, &
       !      y, &
       !      method%commVec,&
       !      method%nBlock &
       !      )
    case default
       if(method%commVec%myrank.eq.0)then
          write(*,*)"error :: not support  ovl = ",iovlflag
       endif
    end select

  end subroutine calc_Matvec_CRS

  ! y = Ax
  subroutine calc_Matvec_DDM(method,x,y)
    use krylov
    use comm_kernel
    use ddm_matv,only : matv_ddm,matv_ddm_gpu
    
    use PARCEL_MATRIX_kernel_CPU
    use PARCEL_VEC_kernel

    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: x
    type(PARCEL_VEC) :: y
    integer iovlflag

    iovlflag=method%iovlflag
    select case (iovlflag)
    case(KSP_SPMV_OVL_NONE)
       if(method%flag%machineflag.eq.PARCEL_CPU_FLAG_ON)then
          call matv_ddm( &
               method%ddmMat%icomm, &
               method%ddmMat%mat%nnd_dia, &
               method%ddmMat%n, &
               method%ddmMat%mat%val_dia, &
               method%ddmMat%mat%ioff_dia, &
               x%vec_FP64(1+x%offset),  &
               y%vec_FP64(1+y%offset),  &
               method%ddmMat%myrank, &
               method%ddmMat%ncpu,   &
               method%ddmMat%prm,&
               method%ddmMat%commHide )
       endif
       
       if(method%flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
          call matv_ddm_GPU(&
               method, &
               x%C_vec_FP64,x%offset, &
               y%C_vec_FP64,y%offset)
       endif       

    case(KSP_SPMV_OVL_CORE)

       ! call MatVec_crs_ovl1( &
       !      method%prm, method%crsMat, &
       !      method%rbuffL, x%vec, method%rbuffU, &
       !      y, &
       !      method%commVec,&
       !      method%nBlock)
    case(KSP_SPMV_OVL_RANK)
       ! call MatVec_crs_ovl2(&
       !      method%prm, method%crsMat, &
       !      x, &
       !      y, &
       !      method%commVec,&
       !      method%nBlock &
       !      )
    case default
       if(method%commVec%myrank.eq.0)then
          write(*,*)"error :: not support  ovl = ",iovlflag
       endif
    end select

  end subroutine calc_Matvec_DDM


  subroutine calc_Matvec_DIA(method,x,y)
    use krylov
    use comm_kernel
    use cu_Matrix
    
    use cu_comm
    use PARCEL_MATRIX_kernel_CPU
    use PARCEL_VEC_kernel
    use matrix_kernel,only : MatVec_dia_
    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: x
    type(PARCEL_VEC) :: y
    integer iovlflag
    integer ierr
    real*8 x_cpu(method%prm%n)
    real*8 y_cpu(method%prm%n)

    iovlflag=method%iovlflag
!    n=method%prm%n

    select case (iovlflag)
    case(KSP_SPMV_OVL_NONE)

       if(method%flag%machineflag.eq.PARCEL_CPU_FLAG_ON)then
          call get_nrk_vec( &
               method%commVec, &
               x%vec_FP64(1+x%offset), &
               method%buff)
          call MatVec_dia_FP64_CPU( &
               method%prm, &
               method%diaMat, &
               x%vec_FP64(1+x%offset), &
               y%vec_FP64(1+y%offset), &
               method%commVec, &
               method%buff, &
               method%nBlock)
       endif
       
       if(method%flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
          call Matvec_DIA_GPU(method,&
               method%C_diaMat, &
               x%C_vec_FP64, x%offset, &
               y%C_vec_FP64, y%offset)  
       endif
    
    case(KSP_SPMV_OVL_CORE)
       !       call MatVec_dia_ovl1(method%prm,method%diaMat,x,y,method%commVec,method%buff,method%nBlock)
          write(*,*)"error :: not support  ovl = ",iovlflag
    case(KSP_SPMV_OVL_RANK)
       !       call MatVec_dia_ovl2(method%prm,method%diaMat,x,y,method%commVec,method%nBlock)
       write(*,*)"error :: not support  ovl = ",iovlflag
    case default
       if(method%commVec%myrank.eq.0)then
          write(*,*)"error :: not support  ovl = ",iovlflag
       endif
    end select

  end subroutine calc_Matvec_DIA


  subroutine calc_res_CRS(method,x,b,r)
    use krylov
    use PARCEL_VEC_kernel

    implicit none

    type(KSP) :: method
    type(PARCEL_VEC) x
    type(PARCEL_VEC) b
    type(PARCEL_VEC) r
    type(PARCEL_VEC) alpha
    type(PARCEL_VEC) beta
    integer(kind=intType) n
    integer(kind=intType) parcel_int
    
    n = method%prm%n
    ! type(STATUS_FLAG) :: flag

    ! call init_staus_flag(flag)
    ! call set_machineflag(flag,PARCEL_GPU_FLAG_ON)
    ! call set_dataflag(flag,KSP_DOUBLE_PRECISION)
    ! call set_calcflag(flag,KSP_DOUBLE_PRECISION)

!    real*8 alpha,beta
!    alpha= 1.0d0
!    beta =-1.0d0

    call calc_Matvec_CRS(method,x,r)

    call parcel_sync(method)
    parcel_int = 1
    call init_vec(method,alpha,parcel_int)
    call init_vec(method,beta ,parcel_int)
    call parcel_sync(method)

    call set_val_FP64(method,alpha,1.0d0,parcel_int)
    call set_val_FP64(method,beta,-1.0d0,parcel_int)
    call parcel_sync(method)
    call axpy(method,alpha,beta,b,r,n)

    call parcel_sync(method)
    ! call system("sleep 1")
    ! write(*,*)"free_vec 1 "
    call free_vec(method,alpha)
    ! call system("sleep 1")
    ! write(*,*)"free_vec 2 "
    call free_vec(method,beta)
    ! call system("sleep 1")
    ! write(*,*)"free_vec end "
    call parcel_sync(method)
    
    return 

  end subroutine calc_res_CRS

  subroutine calc_res_DDM(method,x,b,r)
    use krylov
    use PARCEL_VEC_kernel

    implicit none

    type(KSP) :: method
    type(PARCEL_VEC) x
    type(PARCEL_VEC) b
    type(PARCEL_VEC) r
    type(PARCEL_VEC) alpha
    type(PARCEL_VEC) beta
    integer(kind=intType) n
    integer(kind=intType) parcel_int
    
    n = method%prm%n
    ! type(STATUS_FLAG) :: flag

    call calc_Matvec_DDM(method,x,r)

    call parcel_sync(method)
    parcel_int = 1
    call init_vec(method,alpha,parcel_int)
    call init_vec(method,beta ,parcel_int)
    call parcel_sync(method)

    ! CPU + FP64のみ対応
    parcel_int = 1
    call set_val_FP64(method,alpha,1.0d0,parcel_int)
    call set_val_FP64(method,beta,-1.0d0,parcel_int)
    call parcel_sync(method)
    call axpy(method,alpha,beta,b,r,n)

    call parcel_sync(method)
    ! call system("sleep 1")
    ! write(*,*)"free_vec 1 "
    call free_vec(method,alpha)
    ! call system("sleep 1")
    ! write(*,*)"free_vec 2 "
    call free_vec(method,beta)
    ! call system("sleep 1")
    ! write(*,*)"free_vec end "
    call parcel_sync(method)
    
    return 

  end subroutine calc_res_DDM

  subroutine calc_res_DIA(method,x,b,r)
    use krylov
    use PARCEL_VEC_kernel

    implicit none

    type(KSP) :: method
    type(PARCEL_VEC) x
    type(PARCEL_VEC) b
    type(PARCEL_VEC) r
    type(PARCEL_VEC) alpha
    type(PARCEL_VEC) beta
    integer(kind=intType) n
    integer(kind=intType) parcel_int
    
    n = method%prm%n

    call calc_Matvec_DIA(method,x,r)

    call parcel_sync(method)
    parcel_int = 1
    call init_vec(method,alpha,parcel_int)
    call init_vec(method,beta ,parcel_int)
    call parcel_sync(method)

    ! CPU + FP64のみ対応
    call set_val_FP64(method,alpha,1.0d0,parcel_int)
    call set_val_FP64(method,beta,-1.0d0,parcel_int)
    call parcel_sync(method)
    call axpy(method,alpha,beta,b,r,n)

    call parcel_sync(method)
    ! call system("sleep 1")
    ! write(*,*)"free_vec 1 "
    call free_vec(method,alpha)
    ! call system("sleep 1")
    ! write(*,*)"free_vec 2 "
    call free_vec(method,beta)
    ! call system("sleep 1")
    ! write(*,*)"free_vec end "
    call parcel_sync(method)
    
    return 

  end subroutine calc_res_DIA


  ! r = b-Ax
  subroutine calc_res(method,x,b,r)
    use krylov
    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: x
    type(PARCEL_VEC) :: b
    type(PARCEL_VEC) :: r
    integer MatFormat 

    ! まずはCRS形式のみ対応
    MatFormat = method%prm%MatFormat 
    select case (MatFormat)
    case(MATRIX_FORMAT_CRS)
       call calc_res_CRS(method,x,b,r)
    case(MATRIX_FORMAT_DIA)
       call calc_res_DIA(method,x,b,r)
   case(MATRIX_FORMAT_DDM)
!      call calc_res_DDM(method%ddmMat,x,b,r)
      call calc_res_DDM(method,x,b,r)
    end select
  end subroutine calc_res

  subroutine calc_Matvec(method,x,y)
    use krylov
    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: x
    type(PARCEL_VEC) :: y

    integer MatFormat 
    MatFormat = method%prm%MatFormat 
    select case (MatFormat)
    case(MATRIX_FORMAT_CRS)
       call calc_Matvec_CRS(method,x,y)
    case(MATRIX_FORMAT_DIA)
       call calc_Matvec_DIA(method,x,y)
    case(MATRIX_FORMAT_DDM)
!       call calc_Matvec_DDM(method%ddmMat,x,y)
       call calc_Matvec_DDM(method,x,y)
    end select

  end subroutine calc_Matvec

end module PARCEL_KRYLOV_MATRIX

module PARCEL_PRECON_kernel_GPU
  implicit none
contains

  ! commVec%C_comm <-- icomm
  subroutine set_AS_vec_FP64_GPU( &
       method, &
       n,icomm,addL, &
       myrank,Urank,Lrank, &
       C_s_buffL , C_s_reqL  , C_s_staL  , &
       C_r_buffL , C_r_reqL  , C_r_staL  , &
       C_s_buffU , C_s_reqU  , C_s_staU  , &
       C_r_buffU , C_r_reqU  , C_r_staU  , &
       x,  &
       C_y &
       )

    use krylov,only : KSP
    use c_wrapper,only : devicemem
    use vec_type,only : PARCEL_VEC
    use cu_Vec,only : &
         CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
         zero_init_FP64_GPU, &
         axpy_FP64_GPU
    use PARCEL_VEC_kernel,only : parcel_sync
    use PARCEL_integer,only : intType

    use cu_common_wrapper,only : &
         C_MPI_Isend_FP64,C_MPI_Irecv_FP64,C_MPI_Wait, &
         cudaMemcpy_d2h_fp64,cudaMemcpy_h2d_fp64
    use mpi

    implicit none
    type(KSP) :: method
    integer(kind=intType) n
    type(devicemem) icomm
    integer addL
    integer myrank,Urank,Lrank
    type(PARCEL_VEC) :: x
    type(devicemem)  :: C_y

    integer(kind=intType) i 

    type(devicemem) :: C_s_buffL
    type(devicemem) :: C_s_reqL
    type(devicemem) :: C_s_staL

    type(devicemem) :: C_r_buffL
    type(devicemem) :: C_r_reqL
    type(devicemem) :: C_r_staL

    type(devicemem) :: C_s_buffU
    type(devicemem) :: C_s_reqU
    type(devicemem) :: C_s_staU

    type(devicemem) :: C_r_buffU
    type(devicemem) :: C_r_reqU
    type(devicemem) :: C_r_staU

    integer ierr
    integer(kind=intType) i0
    integer(kind=intType),parameter :: offset_zero=0

    integer(kind=intType) offset
    integer(kind=intType) offset_x0
    integer(kind=intType) parcel_int
    integer int1

#if CUDA_AWARE

#else
!    integer(kind=intType) parcel_int 
    real*8 s_buffU(addL+1)
    real*8 r_buffU(addL+1)
    real*8 s_buffL(addL+1)
    real*8 r_buffL(addL+1)

    integer s_reqU
    integer r_reqU
    integer s_reqL
    integer r_reqL

    integer s_staU(MPI_STATUS_SIZE)
    integer r_staU(MPI_STATUS_SIZE)
    integer s_staL(MPI_STATUS_SIZE)
    integer r_staL(MPI_STATUS_SIZE)

#endif

    offset    = x%offset
    offset_x0 = x%offset

    parcel_int = addL
    call zero_init_FP64_GPU(parcel_int,C_r_buffU,offset_zero)
    call zero_init_FP64_GPU(parcel_int,C_r_buffL,offset_zero)

    if(Urank.ne.-1)then
       ! s_buffU(i)=x(i)
       offset = x%offset

       parcel_int = addL
       call axpy_FP64_GPU( &
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            parcel_int,  &
            1.0d0, &
            0.0d0, &
            x%C_vec_FP64,offset, &
            C_s_buffU   ,offset_zero)
       call parcel_sync(method)
#if CUDA_AWARE
       int1 = 1+offset_zero
       call C_MPI_Irecv_FP64(&
            C_r_buffU, &
            addL, &
            Urank, &
            Urank, &
            icomm, &
            C_r_reqU, &
            1,int1 &
            ) 
       call C_MPI_Isend_FP64( &
            C_s_buffU,&
            addL ,&
            Urank,&
            myrank,&
            icomm, &
            C_s_reqU, &
            1,int1 &
            )
#else
       parcel_int = addL
       call cudaMemcpy_d2h_fp64(C_s_buffU, s_buffU, parcel_int, offset_zero , offset_zero)
       call parcel_sync(method)
       call MPI_Irecv( &
            r_buffU,&
            addL ,&
            MPI_REAL8, &
            Urank,&
            Urank,&
            method%commVec%comm,&
            r_reqU, &
            ierr)

       call MPI_Isend( &
            s_buffU,&
            addL ,&
            MPI_REAL8,&
            Urank,&
            myrank,&
            method%commVec%comm,&
            s_reqU,&
            ierr)
#endif
    endif

    if(Lrank.ne.-1)then
       ! do i=1,addL
       !    s_buffL(i)=x(n-addL+i)
       ! enddo
       offset = x%offset + n-addL
       parcel_int = addL
       call axpy_FP64_GPU( &
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            parcel_int,  &
            1.0d0, &
            0.0d0, &
            x%C_vec_FP64, offset, &
            C_s_buffL   , offset_zero )

       ! r_buffL=0.0d0
       parcel_int = addL
       call zero_init_FP64_GPU(parcel_int,C_r_buffL,offset_zero)

       ! call MPI_Irecv(r_buffL,addL,MPI_REAL8,Lrank,Lrank ,comm,r_Lreq,ierr)
       ! call MPI_Isend(s_buffL,addL,MPI_REAL8,Lrank,myrank,comm,s_Lreq,ierr)
       call parcel_sync(method)
#if CUDA_AWARE
       int1 = 1+offset_zero
       call C_MPI_Irecv_FP64(&
            C_r_buffL, &
            addL, &
            Lrank, &
            Lrank, &
            icomm, &
            C_r_reqL, &
            1,int1 &
            ) 
       call C_MPI_Isend_FP64( &
            C_s_buffL,&
            addL ,&
            Lrank,&
            myrank,&
            icomm, &
            C_s_reqL, &
            1,int1 &
            )
#else

       parcel_int = addL
       call cudaMemcpy_d2h_fp64(C_s_buffL, s_buffL, parcel_int, offset_zero , offset_zero)
       call parcel_sync(method)
       call MPI_Irecv( &
            r_buffL,&
            addL ,&
            MPI_REAL8, &
            Lrank,&
            Lrank,&
            method%commVec%comm,&
            r_reqL, &
            ierr)

       call MPI_Isend( &
            s_buffL,&
            addL ,&
            MPI_REAL8,&
            Lrank,&
            myrank,&
            method%commVec%comm,&
            s_reqL,&
            ierr)

#endif

    endif


    i0 = 0
    if(Urank.ne.-1)then
       ! call MPI_Wait(r_Ureq,r_Usta,ierr)
       ! call MPI_Wait(s_Ureq,s_Usta,ierr)
#if CUDA_AWARE
       call C_MPI_Wait(C_r_reqU,C_r_staU,1)
       call C_MPI_Wait(C_s_reqU,C_s_staU,1)
#else
       parcel_int = addL
       call MPI_Wait(r_reqU,r_staU,ierr)
       call MPI_Wait(s_reqU,s_staU,ierr)
       call cudaMemcpy_h2d_fp64(r_buffL,C_r_buffL, parcel_int, offset_zero , offset_zero)
       call parcel_sync(method)
#endif
       ! do i=1,addL
       !    y(i+i0)=r_buffU(i)
       ! enddo
       offset = i0
       parcel_int = addL
       call axpy_FP64_GPU( &
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            parcel_int,  &
            1.0d0, &
            0.0d0, &
            C_r_buffU   ,offset_zero, &
            C_y       ,offset)
       i0=i0+addL
    endif
!    call xcpy(n,x,y(1+i0))
    offset = i0
    call axpy_FP64_GPU( &
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            n,  &
            1.0d0, &
            0.0d0, &
            x%C_vec_FP64, x%offset, &
            C_y       , offset)

    i0 = i0+n
    if(Lrank.ne.-1)then
       ! call MPI_Wait(r_Lreq,r_Lsta,ierr)
       ! call MPI_Wait(s_Lreq,s_Lsta,ierr)
#if CUDA_AWARE
       call C_MPI_Wait(C_r_reqL,C_r_staL,1)
       call C_MPI_Wait(C_s_reqL,C_s_staL,1)
#else
       parcel_int = addL
       call MPI_Wait(r_reqL,r_staL,ierr)
       call MPI_Wait(s_reqL,s_staL,ierr)
       call cudaMemcpy_h2d_fp64(r_buffL,C_r_buffL, parcel_int, offset_zero , offset_zero)
       call parcel_sync(method)
#endif
       ! do i=1,addL
       !    y(i0+i)=r_buffL(i)
       ! enddo
       offset = i0
       parcel_int = addL
       call axpy_FP64_GPU( &
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            parcel_int,  &
            1.0d0, &
            0.0d0, &
            C_r_buffL   ,offset_zero, &
            C_y       ,offset)
    endif
    x%offset = offset_x0 

  end subroutine set_AS_vec_FP64_GPU

  subroutine get_AS_vec_FP64_GPU( &
       method, &
       n, icomm, addL, &
       myrank, Urank, Lrank, &
       C_s_buffL , C_s_reqL  , C_s_staL  , &
       C_r_buffL , C_r_reqL  , C_r_staL  , &
       C_s_buffU , C_s_reqU  , C_s_staU  , &
       C_r_buffU , C_r_reqU  , C_r_staU  , &
       C_x, &
       y)
    use krylov,only : KSP
    use c_wrapper,only : devicemem
    use vec_type,only : PARCEL_VEC
    use cu_Vec,only : &
         CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
         zero_init_FP64_GPU, &
         axpy_FP64_GPU
    use PARCEL_VEC_kernel,only : parcel_sync
    use PARCEL_integer,only : intType
    use cu_common_wrapper,only : &
         C_MPI_Isend_FP64,C_MPI_Irecv_FP64,C_MPI_Wait, &
         cudaMemcpy_d2h_fp64,cudaMemcpy_h2d_fp64

    use mpi

  
    implicit none
    type(KSP) :: method
    integer(kind=intType) n
    type(devicemem) :: icomm
    integer addL
    integer myrank,Urank,Lrank
    type(devicemem) :: C_x
    type(PARCEL_VEC) :: y
    integer(kind=intType) iflag ! オーバーラップ領域の補間方法

    type(devicemem) :: C_s_buffL
    type(devicemem) :: C_s_reqL
    type(devicemem) :: C_s_staL
    type(devicemem) :: C_r_buffL
    type(devicemem) :: C_r_reqL
    type(devicemem) :: C_r_staL
    type(devicemem) :: C_s_buffU
    type(devicemem) :: C_s_reqU
    type(devicemem) :: C_s_staU
    type(devicemem) :: C_r_buffU
    type(devicemem) :: C_r_reqU
    type(devicemem) :: C_r_staU

    integer(kind=intType) i 
    integer ierr
    integer(kind=intType) i0
    integer(kind=intType),parameter :: offset_zero=0
    integer(kind=intType) offset
    integer(kind=intType) parcel_int
    integer int1

#if CUDA_AWARE

#else
    real*8 s_buffU(addL+1)
    real*8 r_buffU(addL+1)
    real*8 s_buffL(addL+1)
    real*8 r_buffL(addL+1)

    integer s_reqU
    integer r_reqU
    integer s_reqL
    integer r_reqL

    integer s_staU(MPI_STATUS_SIZE)
    integer r_staU(MPI_STATUS_SIZE)
    integer s_staL(MPI_STATUS_SIZE)
    integer r_staL(MPI_STATUS_SIZE)
#endif


    ! ---
    ! iflag = 0 :: 相加平均
    ! iflag = 1 :: 破棄
    ! ---
 
    iflag = 0
    ! if(iflag.eq.1)then
    !    return
    ! endif

    if(Urank.ne.-1)then
       i0=addL
    else
       i0=0
    endif

    parcel_int = addL
    call zero_init_FP64_GPU(parcel_int,C_r_buffU,offset_zero)
    call zero_init_FP64_GPU(parcel_int,C_r_buffL,offset_zero)
    if(Urank.ne.-1)then
       ! do i=1,addL
       !    s_buffU(i)=x(i)
       ! enddo
       offset = 0
       parcel_int = addL
       call axpy_FP64_GPU( &
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            parcel_int,  &
            1.0d0, &
            0.0d0, &
            C_x         ,offset, &
            C_s_buffU   ,offset_zero)
!       call MPI_Irecv(r_buffU,addL,MPI_REAL8,Urank,Urank ,icomm,r_Ureq,ierr)
!       call MPI_Isend(s_buffU,addL,MPI_REAL8,Urank,myrank,icomm,s_Ureq,ierr)
       call parcel_sync(method)

#if CUDA_AWARE
       int1 = 1+offset_zero
       call C_MPI_Irecv_FP64(&
            C_r_buffU, &
            addL, &
            Urank, &
            Urank, &
            icomm, &
            C_r_reqU, &
            1,int1 &
            ) 
       call C_MPI_Isend_FP64( &
            C_s_buffU,&
            addL ,&
            Urank,&
            myrank,&
            icomm, &
            C_s_reqU, &
            1,int1 &
            )
#else
       parcel_int = addL
       call cudaMemcpy_d2h_fp64(C_s_buffU, s_buffU, parcel_int, offset_zero , offset_zero)
       call parcel_sync(method)
       call MPI_Irecv( &
            r_buffU,&
            addL ,&
            MPI_REAL8, &
            Urank,&
            Urank,&
            method%commVec%comm,&
            r_reqU, &
            ierr)

       call MPI_Isend( &
            s_buffU,&
            addL ,&
            MPI_REAL8,&
            Urank,&
            myrank,&
            method%commVec%comm,&
            s_reqU,&
            ierr)

#endif

    endif

    if(Lrank.ne.-1)then
       ! do i=1,addL
       !    s_buffL(i)=x(i0+n+i)
       ! enddo

       offset = n+i0
       parcel_int = addL
       call axpy_FP64_GPU( &
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            parcel_int,  &
            1.0d0, &
            0.0d0, &
            C_x         ,offset, &
            C_s_buffL   ,offset_zero)
       call parcel_sync(method)

       ! call MPI_Irecv(r_buffL,addL,MPI_REAL8,Lrank,Lrank ,icomm,r_Lreq,ierr)
       ! call MPI_Isend(s_buffL,addL,MPI_REAL8,Lrank,myrank,icomm,s_Lreq,ierr)
#if CUDA_AWARE
       int1 = 1+offset_zero
       call C_MPI_Irecv_FP64(&
            C_r_buffL, &
            addL, &
            Lrank, &
            Lrank, &
            icomm, &
            C_r_reqL, &
            1,int1 &
            ) 
       call C_MPI_Isend_FP64( &
            C_s_buffL,&
            addL ,&
            Lrank,&
            myrank,&
            icomm, &
            C_s_reqL, &
            1,int1 &
            )
#else
       parcel_int = addL
       call cudaMemcpy_d2h_fp64(C_s_buffL, s_buffL, parcel_int, offset_zero , offset_zero)
       call parcel_sync(method)
       call MPI_Irecv( &
            r_buffL,&
            addL ,&
            MPI_REAL8, &
            Lrank,&
            Lrank,&
            method%commVec%comm,&
            r_reqL, &
            ierr)

       call MPI_Isend( &
            s_buffL,&
            addL ,&
            MPI_REAL8,&
            Lrank,&
            myrank,&
            method%commVec%comm,&
            s_reqL,&
            ierr)

#endif
    endif

    if(Urank.ne.-1)then
       ! call MPI_Wait(r_Ureq,r_Usta,ierr)
       ! call MPI_Wait(s_Ureq,s_Usta,ierr)
#if CUDA_AWARE
       call C_MPI_Wait(C_r_reqU,C_r_staU,1)
       call C_MPI_Wait(C_s_reqU,C_s_staU,1)
#else
       call MPI_Wait(r_reqU,r_staU,ierr)
       call MPI_Wait(s_reqU,s_staU,ierr)
       parcel_int = addL
       call cudaMemcpy_h2d_fp64(r_buffL,C_r_buffL, parcel_int, offset_zero , offset_zero)
       call parcel_sync(method)
#endif
       ! do i=1,addL
       ! s(i)= ( s(i) + r_buffU(i) ) 
       ! enddo
       parcel_int = addL
       call axpy_FP64_GPU( &
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            parcel_int,  &
            1.0d0, &
            1.0d0, &
            C_r_buffU,offset_zero,   &
            y%C_vec_FP64 ,y%offset &
            )
    endif

    if(Lrank.ne.-1)then
       ! call MPI_Wait(r_Lreq,r_Lsta,ierr)
       ! call MPI_Wait(s_Lreq,s_Lsta,ierr)
#if CUDA_AWARE
       call C_MPI_Wait(C_r_reqL,C_r_staL,1)
       call C_MPI_Wait(C_s_reqL,C_s_staL,1)
#else
       call MPI_Wait(r_reqL,r_staL,ierr)
       call MPI_Wait(s_reqL,s_staL,ierr)
       parcel_int = addL
       call cudaMemcpy_h2d_fp64(r_buffL,C_r_buffL, parcel_int , offset_zero , offset_zero)
       call parcel_sync(method)
#endif
       ! do i=1,addL
       ! s(n-addL+i)= ( s(n-addL+i) + r_buffL(i) ) 
       ! enddo
       offset = y%offset + n-addL
       parcel_int = addL
       call axpy_FP64_GPU( &
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            parcel_int,  &
            1.0d0, &
            1.0d0, &
            C_r_buffL,offset_zero,   &
            y%C_vec_FP64 , offset &
            )

    endif

  end subroutine get_AS_vec_FP64_GPU

  subroutine AdditiveSchwarz_GPU(method,r,s)
    use krylov,only : KSP
    use vec_type,only : PARCEL_VEC
    use PARCEL_integer,only : intType
    use PARCEL_define,only : &
         KSP_ADDITIVESCHWARZ_TYPE_BASIC, &
         KSP_ADDITIVESCHWARZ_TYPE_RESTRICT, &
         KSP_ADDITIVESCHWARZ_TYPE_INTERPOLATE, &
         KSP_ADDITIVESCHWARZ_TYPE_NONE
    use PARCEL_VEC_kernel,only : &
         set_val_FP64,init_vec,free_vec, &
         parcel_sync
    use cu_Vec,only : &
         CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
         zero_init_FP64_GPU, &
         axpy_FP64_GPU
    use cu_common_wrapper,only : &
         set_vec_val_FP64_GPU,cusparse_precon, &
         cudaMemcpy_d2h_FP64

    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: r
    type(PARCEL_VEC) :: s
    integer(kind=intType),parameter :: offset_zero=0
    integer(kind=intType) offset

    ! -----------------------

    integer(kind=intType) n
    integer myrank
    integer iflagAS
    integer addL
    integer(kind=intType) addL_intType
    integer(kind=intType) i,i0
    type(PARCEL_VEC) :: tmp
    real*8,allocatable :: tmpvec(:)
    integer(kind=intType) parcel_int

    parcel_int = 1
    call init_vec(method,tmp,parcel_int)
    call set_val_FP64(method,tmp,0.0d0,parcel_int)

    n =  method%prm%n
    myrank = method%commVec%myrank
    iflagAS      = method%iflagAS
    addL         = method%addL
    addL_intType = method%addL

    allocate(tmpvec(n+addL*2))
    tmpvec = 0.0d0
!    call zero_initialize(addL_intType,vecL)
!    call zero_initialize(addL_intType,vecU)
    call parcel_sync(method)
    parcel_int = addL
    call zero_init_FP64_GPU(parcel_int,method%C_s_buffL,offset_zero)
    call zero_init_FP64_GPU(parcel_int,method%C_s_buffU,offset_zero)
    call zero_init_FP64_GPU(parcel_int,method%C_r_buffL,offset_zero)
    call zero_init_FP64_GPU(parcel_int,method%C_r_buffU,offset_zero)

#if 0
    call cudaMemcpy_d2h_FP64( &
         r%C_vec_FP64,tmpvec,n, &
         offset_zero,offset_zero    )         
    do i = 1,10
       write(*,*)" r vector  : ",i,tmpvec(i)
    enddo
#endif

!    return
    if(method%Urank_AS.ne.-1)then
       i0 = addL_intType
    else
       i0 = 0
    endif

    if( (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_BASIC) .or. &
         (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_RESTRICT) )then
       
       call set_AS_vec_FP64_GPU( &
            method, &
            method%prm%n, method%commVec%C_comm, method%addL, &
            method%commVec%myrank, method%Urank_AS, method%Lrank_AS, &
            method%C_s_buffL, method%C_s_reqL, method%C_s_staL, &
            method%C_r_buffL, method%C_r_reqL, method%C_r_staL, &
            method%C_s_buffU, method%C_s_reqU, method%C_s_staU, &
            method%C_r_buffU, method%C_r_reqU, method%C_r_staU, &
            r, &
            method%C_ivec_AS &
            )

    elseif( (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_INTERPOLATE) .or. &
         (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_NONE) )then       
       offset = offset_zero
       call zero_init_FP64_GPU(addL_intType,method%C_ivec_AS,offset)
       offset = n + i0
       call zero_init_FP64_GPU(addL_intType,method%C_ivec_AS,offset)

       offset = i0
!       write(*,*)" offset = ",offset
       call axpy_FP64_GPU( &
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            method%prm%n, &
            1.0d0, &
            0.0d0, &
            r%C_vec_FP64     , r%offset, &
            method%C_ivec_AS , offset)

    endif
    ! write(*,*)"iflagAS = ",iflagAS

    call parcel_sync(method)

#if 0
    call cudaMemcpy_d2h_FP64( &
         method%C_ivec_AS,tmpvec,n, &
         offset_zero,offset_zero    )         
    do i = 1,10
       write(*,*)" input  : ",i,tmpvec(i)
    enddo
#endif

    call cusparse_precon( &
               method%cusparce_solver, &
               method%C_ivec_AS, offset_zero,&
               method%C_ovec_AS, offset_zero,&
               method%cusparce_preconcrsMatrix%crsA, &
               method%cusparce_preconcrsMatrix%crsCol, &
               method%cusparce_preconcrsMatrix%crsRow_ptr, &
               method%cusparce_preconcrsMatrix%n, &
               method%cusparce_preconcrsMatrix%nnz &
               )

    call parcel_sync(method)
    ! -------    
    ! スレッド間の加算が不要なのでいらない
    ! call AS_Local( &
    !      method%prm_AS, &
    !      method%preconAS, &
    !      method%max_numth, &
    !      iflagAS, &
    !      method%ovec_AS, &
    !      method%precon_thblock)
    ! -------
    !  call xcpy(method%prm%n,method%ovec_AS(1+i0),s) 

#if 0
    myrank = method%commVec%myrank
    call cudaMemcpy_d2h_FP64( &
         method%C_ovec_AS,tmpvec,n, &
         s%offset,offset_zero    )         
    do i = 1,10
       write(*,*)myrank," output : ",i,tmpvec(i)
    enddo
#endif

    offset = i0
    ! write(*,*)" n        = ",method%prm%n
    ! write(*,*)" precon_n = ",method%cusparce_preconcrsMatrix%n
    ! write(*,*)" offset   = ",offset
    call axpy_FP64_GPU( &
         CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
         method%prm%n, &
         1.0d0, &
         0.0d0, &
         method%C_ovec_AS       , offset, &
         s%C_vec_FP64           , s%offset)

#if 0
    call cudaMemcpy_d2h_FP64( &
         s%C_vec_FP64,tmpvec,n, &
         s%offset,offset_zero    )         
    do i = 1,10
       write(*,*)" output : ",i,tmpvec(i)
    enddo
#endif

    if((iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_RESTRICT) .or. &
         (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_NONE))then
       return 
    endif

    ! -----
    ! オーバーラップ領域の補正
    ! -----
    call get_AS_vec_FP64_GPU( &
         method, &
         method%prm%n, &
         method%commVec%C_comm, &
         method%addL, &
         method%commVec%myrank, &
         method%Urank_AS, &
         method%Lrank_AS, &
         method%C_s_buffL,method%C_s_reqL,method%C_s_staL, &
         method%C_r_buffL,method%C_r_reqL,method%C_r_staL, &
         method%C_s_buffU,method%C_s_reqU,method%C_s_staU, &
         method%C_r_buffU,method%C_r_reqU,method%C_r_staU, &
         method%C_ovec_AS, &
         s &
         )

#if 0
    if(method%Lrank_AS.ne.-1)then
       do i=1,addL
          s(n-addL+i)= ( s(n-addL+i)+vecL(i) ) 
       enddo
    endif
    if(method%Urank_AS.ne.-1)then
       do i=1,addL 
          s(i)= ( s(i)+vecU(i) ) 
       enddo
    endif

#endif

    call free_vec(method,tmp)
    deallocate(tmpvec)

  end subroutine AdditiveSchwarz_GPU

  subroutine precon_Neumann(method,r,s)
    use vec_type
    use krylov
    use cu_Vec
    use cu_Matrix
    use cu_comm
!    use PARCEL_PRECON_kernel_CPU,only : Neumann_dia_FP64_CPU
    use PARCEL_KRYLOV_MATRIX,only : Matvec_CRS_GPU,Matvec_DIA_GPU
    use ddm_matv,only : matv_ddm_gpu,matv_Neumann_ddm_GPU,jacobi_iter1_ddm_GPU
    use matrix_kernel,only : Neumann_crs,Neumann_dia

    use PARCEL_VEC_kernel,only : parcel_sync
    use PARCEL_integer
    use c_wrapper,only : CUDA_NUM_BLOCKS , CUDA_NUM_THREADS
    use cu_common_wrapper,only : &
         set_vec_val_FP64_GPU

    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: r
    type(PARCEL_VEC) :: s
    integer(kind=intType),parameter :: offset_zero=0

    integer i
    integer k
!    integer nbuff
    integer(kind=intType) nbuff
    integer itr_neumann
    integer(kind=intType) n,nadd
    integer iflag
    integer MatFormat 
    real*8 alpha
    integer(kind=intType) parcel_int
    MatFormat = method%prm%MatFormat 

    if(method%flag%machineflag.eq.PARCEL_CPU_FLAG_ON)then
       select case (MatFormat)
       case(MATRIX_FORMAT_CRS)
          call Neumann_crs(method,r%vec_FP64(1+r%offset),s%vec_FP64(1+s%offset))
       case(MATRIX_FORMAT_DIA)
          call Neumann_dia(method,r%vec_FP64(1+r%offset),s%vec_FP64(1+s%offset))
       case(MATRIX_FORMAT_DDM)
          write(*,*)"not support"
       end select

       return
    endif

    iflag = method%ipreflag
    n     = method%prm%n
    nadd = &
         method%prm%n + &
         method%crsMat%nColL + &
         method%crsMat%nColU            

!    itr_neumann = ITRATION_JACOBI
    itr_neumann = method%C_preconNeuman%itr

    !    method%C_preconNeuman%local_flag = KSP_PRECONDITION_NEUMANN_TRUE_COMM-10
    !    method%C_preconNeuman%local_flag = KSP_PRECONDITION_NEUMANN_TRUE_COMM
    !    itr_neumann = 4
    if(MatFormat.eq.MATRIX_FORMAT_DDM)then
       ! s  = Dinv * r
       if(itr_neumann.eq.1)then
          ! write(*,*)"Neumann 1 iter"
          ! x(i+1) = x(i) + Dinv*(r - Ax(i) )
          ! -----------------------------------
          ! i = 1の場合
          ! x(1) = Dinv*r
          ! x(2) = x(1) + Dinv*r - Dinv*Ax(1) 
          ! x(2) = x(1) + x(1) - Dinv*Ax(1) 
          ! x(2) =  2.0 * x(1) - Dinv*Ax(1)

          !  x(1) = Dinv*r の計算
          call vec_scaling_FP64_GPU( &
               CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
               n, &
               r%C_vec_FP64, r%offset, &
               method%C_preconNeuman%Dinv, offset_zero, &
               method%C_preconNeuman%vecx, offset_zero)
          ! x(2) = 2.0 * x(1) - Dinv*Ax(1) --> x(1) - Dinv*(L+U)x(1)
          if( method%C_preconNeuman%local_flag .eq. &
               KSP_PRECONDITION_NEUMANN_TRUE_COMM)then
             call jacobi_iter1_ddm_GPU(method, &
                  method%C_preconNeuman%vecx,offset_zero, &
                  s%C_vec_FP64, s%offset)
          else
             call jacobi_iter1_ddm_GPU(method, &
                  method%C_preconNeuman%vecx,offset_zero, &
                  s%C_vec_FP64, s%offset, &
                  1)
          endif
       else
          ! s = Dinv * r
          ! x = s
          call neumann_preparation_FP64_GPU( &
               CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
               n, &
               r%C_vec_FP64, r%offset, &
               method%C_preconNeuman%Dinv,offset_zero, &
               s%C_vec_FP64, s%offset, &
               method%C_preconNeuman%vecx,offset_zero)

          ! x(i+1) = x(i) + Dinv*(r - Ax(i) )
          alpha = 1.0
          do i = 1,itr_neumann 
             ! : r0 = Dinv * r
             ! : Dinv*A = (I+LU)
             ! s = s + r0 - (I+LU)s
             ! s = s + r0 -s -LU*s
             ! s = r0 -LU*s
             ! --------
             ! s(0) = 0
             ! s(1) = r0
             ! s(2) = r0 - LU*s(1) --> s(1) - LU*r0
             ! s(3) = r0 - LU*s(2) --> r0 - LU*(r0 - LU*r0) --> r0 - LU*r0 + LU*LU*r0 --> s(2) + LU*LU*r0
             ! s(4) = r0 - LU*s(3) --> r0 - LU*(s(2) + LU*LU*r0) -->  r0 - LU*s(2) - LU*LU*LU*r0 -->  s(3) - LU*LU*LU*r0 
             ! s(n) = s(n-1) + (LU)^(n-1) * (-1)^(n-1) * r0
             ! --------

             ! y = Ax(i)
             if(method%C_preconNeuman%local_flag.eq.KSP_PRECONDITION_NEUMANN_TRUE_COMM)then
                call matv_Neumann_ddm_GPU(method, &
                     method%C_preconNeuman%vecx,offset_zero, &
                     method%C_preconNeuman%vecy,offset_zero)
             else
                call matv_Neumann_ddm_GPU(method, &
                     method%C_preconNeuman%vecx,offset_zero, &
                     method%C_preconNeuman%vecy,offset_zero, &
                     1)
             endif
             ! x = y
             ! s = s + y
             call neumann_update_FP64_GPU( &
                  CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
                  n, &
                  alpha,&
                  method%C_preconNeuman%vecy, offset_zero,&
                  method%C_preconNeuman%vecx, offset_zero,&
                  s%C_vec_FP64, s%offset &
                  )
!             alpha = -alpha 
          enddo
       endif

    else
!       write(*,*)"itr_neumann = ",itr_neumann
       if(itr_neumann.eq.1)then
          ! x = Dinv * r
          call vec_scaling_FP64_GPU( &
               CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
               n, &
               r%C_vec_FP64, r%offset, &
               method%C_preconNeuman%Dinv,offset_zero, &
               method%C_preconNeuman%vecx,offset_zero )
          ! y = x + Ax 
          select case (MatFormat)
          case(MATRIX_FORMAT_CRS)
!             return 
             call parcel_sync(method)
             if(method%C_preconNeuman%local_flag.eq.KSP_PRECONDITION_NEUMANN_TRUE_COMM)then
                call get_nrk_vec_FP64_GPU(&
                     method%commVec, &
                     method%C_preconNeuman%vecx, offset_zero, &
                     method%C_sbuff, offset_zero, &
                     method%C_buff,  offset_zero )
                parcel_int = 0
                call get_nrk_vecLU_add_FP64_GPU( &
                     method%C_buff, parcel_int,&
                     method%prm%n, &
                     method%C_crsMat%nColL, &
                     method%C_crsMat%nColU, &
                     method%C_preconNeuman%vecx, offset_zero &
                     )
             else
                ! 通信なし
                nbuff = method%C_crsMat%nColL + method%C_crsMat%nColU
                call set_vec_val_FP64_GPU(&
                     CUDA_NUM_BLOCKS, &
                     CUDA_NUM_THREADS, &
                     nbuff, &
                     0.0d0, &
                     method%C_preconNeuman%vecx,method%prm%n &
                     )
             endif
             call parcel_sync(method)
             if(s%offset.ne.offset_zero)then
                call cusparseDestroyDnVec(s%C_vec_FP64)
                call cusparseCreateDnVec_FP64(s%C_vec_FP64,n,nadd,s%offset)
             endif
             call parcel_sync(method)
             call Neumann_SpMV( &
                  method%cusparce_solver    ,  &
                  method%C_preconNeuman%vecx, offset_zero, &
                  s%C_vec_FP64 ,s%offset &
                  )
             call parcel_sync(method)
             if(s%offset.ne.offset_zero)then
                call cusparseDestroyDnVec(s%C_vec_FP64)
                call cusparseCreateDnVec_FP64(s%C_vec_FP64,n,nadd,offset_zero)
             endif
             call parcel_sync(method)
 
             call axpy_FP64_GPU(&
                  CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
                  method%prm%n, &
                  1.0d0,&
                  1.0d0,&
                  method%C_preconNeuman%vecx,offset_zero, &
                  s%C_vec_FP64,s%offset   &
                  ) 
          case(MATRIX_FORMAT_DIA)
             if(method%C_preconNeuman%local_flag.eq.KSP_PRECONDITION_NEUMANN_TRUE_COMM)then
                call parcel_sync(method)
                call get_nrk_vec_FP64_GPU(&
                     method%commVec, &
                     method%C_preconNeuman%vecx, offset_zero, &
                     method%C_sbuff, offset_zero, &
                     method%C_buff , offset_zero )
                call parcel_sync(method)             
             else
                ! 通信なし
                nbuff = 0
                do k=1,method%commVec%nrk_rank_send
                   nbuff = nbuff + method%commVec%nrk_prm_send(k)%ndata
                enddo
                call set_vec_val_FP64_GPU(&
                     CUDA_NUM_BLOCKS, &
                     CUDA_NUM_THREADS, &
                     nbuff, &
                     0.0d0, &
                     method%C_buff,offset_zero &
                     )
             endif

             call parcel_sync(method)
!             write(*,*)"neuman : nnd = ",method%C_preconNeuman%diaMat%nnd
!             return
             parcel_int = method%nBlock
             call MatVec_dia_neumann_FP64_GPU( &
                  CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
                  method%prm%n,&
                  method%C_preconNeuman%diaMat%nnd,&
                  method%prm%istart,&
                  parcel_int, &
                  method%C_preconNeuman%vecx,offset_zero, &
                  s%C_vec_FP64,s%offset, &
                  method%C_preconNeuman%diaMat%diaA, &
                  method%C_preconNeuman%diaMat%offset , &
                  method%C_preconNeuman%diaMat%iLs, &
                  method%C_buff &
                  )
          case(MATRIX_FORMAT_DDM)
          end select

       else
          ! s = Dinv * r
          ! x = s

          call neumann_preparation_FP64_GPU( &
               CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
               n, &
               r%C_vec_FP64, r%offset,&
               method%C_preconNeuman%Dinv, offset_zero, &
               s%C_vec_FP64, s%offset,&
               method%C_preconNeuman%vecx, offset_zero)
          do i = 1,itr_neumann 
             ! y = Ax
             select case (MatFormat)
             case(MATRIX_FORMAT_CRS)

                if(method%C_preconNeuman%local_flag.eq.KSP_PRECONDITION_NEUMANN_TRUE_COMM)then
                   call parcel_sync(method)             
                   call get_nrk_vec_FP64_GPU(&
                        method%commVec, &
                        method%C_preconNeuman%vecx, offset_zero, &
                        method%C_sbuff, offset_zero, &
                        method%C_buff, offset_zero)
                   call parcel_sync(method)             
                   call get_nrk_vecLU_add_FP64_GPU( &
                        method%C_buff, offset_zero, &
                        method%prm%n, &
                        method%C_crsMat%nColL, &
                        method%C_crsMat%nColU, &
                        method%C_preconNeuman%vecx,offset_zero &
                        )
                else
                   ! 通信なし
                   ! write(*,*)"nbuff = ",nbuff
                   nbuff = method%C_crsMat%nColL + method%C_crsMat%nColU
                   call set_vec_val_FP64_GPU(&
                        CUDA_NUM_BLOCKS, &
                        CUDA_NUM_THREADS, &
                        nbuff, &
                        0.0d0, &
                        method%C_preconNeuman%vecx,method%prm%n &
                        )
                endif

                call parcel_sync(method)             
                call Neumann_SpMV( &
                     method%cusparce_solver, &
                     method%C_preconNeuman%vecx, offset_zero, &
                     method%C_preconNeuman%vecy, offset_zero  &
                     ) 
                call parcel_sync(method)             
             case(MATRIX_FORMAT_DIA)
                if(method%C_preconNeuman%local_flag.eq.KSP_PRECONDITION_NEUMANN_TRUE_COMM)then
                   call Matvec_DIA_GPU(method,&
                        method%C_preconNeuman%diaMat, &
                        method%C_preconNeuman%vecx, offset_zero, &
                        method%C_preconNeuman%vecy, offset_zero)
                else
                   ! 通信なし
                   call Matvec_DIA_GPU(method,&
                        method%C_preconNeuman%diaMat, &
                        method%C_preconNeuman%vecx, offset_zero, &
                        method%C_preconNeuman%vecy, offset_zero, &
                        -1)
                endif
             end select

             ! x = y
             ! s = s + y
             call neumann_update_FP64_GPU( &
                  CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
                  n, &
                  1.0d0, &
                  method%C_preconNeuman%vecy,offset_zero, &
                  method%C_preconNeuman%vecx,offset_zero, &
                  s%C_vec_FP64, s%offset &
                  )
          enddo
       endif

    endif
  end subroutine precon_Neumann

end module PARCEL_PRECON_kernel_GPU

module PARCEL_KRYLOV_PRECON
  use vec_type
  use PARCEL_VEC_kernel
  use PARCEL_KRYLOV_MATRIX
  implicit none
contains
  subroutine precon_crs(method,r,s)
    use krylov
    use PARCEL_PRECON_kernel_GPU,only : precon_Neumann,AdditiveSchwarz_GPU
    use cu_Vec,only : vec_scaling_FP64_GPU
    use vec_kernel,only : vec_scaling 
    use matrix_kernel,only : &
         blockJacobi_ILU_crs,blockJacobi_crs, &
         additiveschwarz_crs

    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: r
    type(PARCEL_VEC) :: s
    integer iflag
    integer myrank
    integer(kind=intType),parameter :: offset_zero=0
    integer(kind=intType) n
    
    n = method%prm%n
    myrank=method%commVec%myrank
    iflag=method%ipreflag

    ! ------
    ! write(*,*)" precon : iflag = ",iflag
    select case (iflag)
    case(KSP_PRECONDITION_NONE)
       ! iflag = 0 :: none preconditioning
!       call parcel_sync(method)
!       write(*,*)" offset = ",r%offset,s%offset
       call xcpy(method,r,s,n)
!       call parcel_sync(method)
       return
    case(KSP_PRECONDITION_POINTJACOBI)
       ! iflag = 1 :: point Jacobi
       if(method%flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
          call vec_scaling_FP64_GPU( &
               CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
               n, &
               r%C_vec_FP64 , r%offset, &
               method%C_Dinv, offset_zero,&
               s%C_vec_FP64 , s%offset)
       else
          call vec_scaling(&
               method%prm%n,&
               r%vec_FP64(1+r%offset), &
               method%Dinv(1), &
               s%vec_FP64(1+s%offset) &
               )
       endif

       return
    case(KSP_PRECONDITION_BLOCKJACOBI)
       ! iflag = 2 :: Block Jacobi, Additive Schwarz(addL=0)

       if(method%flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
          call cusparse_precon( &
               method%cusparce_solver, &
               r%C_vec_FP64, r%offset,&
               s%C_vec_FP64, s%offset,&
               method%cusparce_preconcrsMatrix%crsA, &
               method%cusparce_preconcrsMatrix%crsCol, &
               method%cusparce_preconcrsMatrix%crsRow_ptr, &
               method%cusparce_preconcrsMatrix%n, &
               method%cusparce_preconcrsMatrix%nnz &
               )
       else
          if(method%ILU_method.eq.KSP_PRECONDITION_ILU)then
             ! ILU(0)
             ! call blockJacobi_ILU_crs(method,r,s)
             ! write(*,*)" not support ILU (0) "
             ! call xcpy(method,r,s,n)
             call blockJacobi_ILU_crs(method,r%vec_FP64(1+r%offset),s%vec_FP64(1+s%offset))
          else
             ! D-ILU(0),M-ILU(0)
             ! call blockJacobi_crs(method,r,s)
             call blockJacobi_crs(method,r%vec_FP64(1+r%offset),s%vec_FP64(1+s%offset))
          endif
       endif

       return

    case(KSP_PRECONDITION_ADDITIVESCHWARZ)
       ! iflag = 3 :: Additive Schwarz (addL>0)
       ! call AdditiveSchwarz_crs(method,r,s)
#if 0
       call xcpy(method,r,s,n)
#else
       if(method%flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
          ! GPU
          call AdditiveSchwarz_GPU(method,r,s)
       else
          ! CPU 
          call AdditiveSchwarz_crs(&
               method,r%vec_FP64(1+r%offset),s%vec_FP64(1+s%offset))
       endif
#endif
       return

    case(KSP_PRECONDITION_NEUMANN)
!       write(*,*)" precon_Neumann "
       call precon_Neumann(method,r,s)
       return
    case(KSP_PRECONDITION_JACOBI_MatBF16)
!       call precon_JACOBI(method,r,s)
       return
    end select

  end subroutine precon_crs

  subroutine precon_dia(method,r,s)
    use krylov
    use cu_Vec
    use cu_Matrix
    use cu_comm
    use PARCEL_PRECON_kernel_GPU
    use PARCEL_PRECON_kernel_CPU
    use vec_kernel,only : vec_scaling 
    use matrix_kernel,only : &
         blockJacobi_ILU_dia,blockJacobi_dia,&
         additiveschwarz_dia

    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: r
    type(PARCEL_VEC) :: s
    integer(kind=intType),parameter :: offset_zero=0
    integer iflag
    integer myrank
    integer i
    integer(kind=intType) n
    integer itr_neumann

    n = method%prm%n
    myrank=method%commVec%myrank
    iflag=method%ipreflag

    ! ------
    !    write(*,*)" precon : iflag = ",iflag
    select case (iflag)
    case(KSP_PRECONDITION_NONE)
       ! iflag = 0 :: none preconditioning
       call xcpy(method,r,s,n)
       return
    case(KSP_PRECONDITION_POINTJACOBI)
       ! iflag = 1 :: point Jacobi
       ! call pointJacobi(method%prm%n,r,method%Dinv,s)
       if(method%flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
          call vec_scaling_FP64_GPU( &
               CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
               n, &
               r%C_vec_FP64,r%offset, &
               method%C_Dinv,offset_zero, &
               s%C_vec_FP64,s%offset )
       else
          call vec_scaling(&
               method%prm%n, &
               r%vec_FP64(1+r%offset), &
               method%Dinv(1),&
               s%vec_FP64(1+s%offset) &
               )
       endif
       return
    case(KSP_PRECONDITION_BLOCKJACOBI)
       if(method%flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
          call cusparse_precon( &
               method%cusparce_solver, &
               r%C_vec_FP64, r%offset,&
               s%C_vec_FP64, s%offset,&
               method%cusparce_preconcrsMatrix%crsA, &
               method%cusparce_preconcrsMatrix%crsCol, &
               method%cusparce_preconcrsMatrix%crsRow_ptr, &
               method%cusparce_preconcrsMatrix%n, &
               method%cusparce_preconcrsMatrix%nnz &
               )
          return 
       else
          if(method%ILU_method.eq.KSP_PRECONDITION_ILU)then
             call blockJacobi_ILU_dia(&
                  method,r%vec_FP64(1+r%offset),s%vec_FP64(1+s%offset))
          else
             call blockJacobi_dia( & 
                  method,r%vec_FP64(1+r%offset),s%vec_FP64(1+s%offset))
          endif
          return 
       endif

       ! iflag = 2 :: Block Jacobi, Additive Schwarz(addL=0)
    case(KSP_PRECONDITION_ADDITIVESCHWARZ)

#if 0
       call xcpy(method,r,s,n)
#else
       if(method%flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
          ! GPU
          call AdditiveSchwarz_GPU(method,r,s)
       else
          ! CPU 
          call AdditiveSchwarz_dia(&
               method,r%vec_FP64(1+r%offset),s%vec_FP64(1+s%offset))
       endif
#endif

       return

    case(KSP_PRECONDITION_NEUMANN)
       ! A = I - F ~ M
       ! (I - F)^-1 = I + F + F^2 + F^3 + ....
       ! Ms = r --> s = (M^-1)r = (I + F + F^2 + F^3 + ....)r
#if 0
       call xcpy(method,r,s,n)
#else
       call precon_Neumann(method,r,s)
#endif
       return
    end select
  end subroutine precon_dia

  subroutine precon_ddm(method,r,s)
    use krylov
    use PARCEL_PRECON_kernel_GPU,only : precon_Neumann
    use cu_Vec,only : vec_scaling_FP64_GPU
    use krylov_kernel,only : calc_precon_DDM
    use precon_DILU,only : &
         solve_DILU_DDM_FP64_GPU,solve_DILU_DDM, &
         solve_DILU_AS_DDM, &
         solve_DILU_AS_DDM_GPU

    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: r
    type(PARCEL_VEC) :: s
    integer(kind=intType),parameter :: offset_zero=0
    integer iflag
    integer myrank
    integer(kind=intType) n
    integer ierr
    
    n = method%prm%n
    myrank=method%commVec%myrank
    iflag=method%ipreflag
    
    if(method%flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
       ! ------
       ! write(*,*)" precon : iflag = ",iflag
       select case (iflag)
       case(KSP_PRECONDITION_NONE)
          ! iflag = 0 :: none preconditioning
          call xcpy(method,r,s,n)
          return
       case(KSP_PRECONDITION_POINTJACOBI)
          call vec_scaling_FP64_GPU( &
               CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
               n, &
               r%C_vec_FP64,r%offset, &
!               method%ddmMat%C_diag,offset_zero, &
               method%C_Dinv,offset_zero, &
               s%C_vec_FP64,s%offset)
          ! iflag = 1 :: point Jacobi
          ! call pointJacobi(method%prm%n,r,method%Dinv,s) 
          return
       case(KSP_PRECONDITION_BLOCKJACOBI)
          ! iflag = 2 :: Block Jacobi, Additive Schwarz(addL=0)
          write(*,*)"GPU : not support BLOCKJACOBI "
          return
       case(KSP_PRECONDITION_ADDITIVESCHWARZ)
          ! iflag = 3 :: Additive Schwarz (addL>0)
          ! call AdditiveSchwarz_crs(method,r,s)
          write(*,*)"GPU : not support ADDITIVESCHWARZ "
          call xcpy(method,r,s,n)
          return
       case(KSP_PRECONDITION_NEUMANN)
          call precon_Neumann(method,r,s)
          ! call xcpy(method,r,s,n)
          return
       case(21)
          call parcel_sync(method)
          call solve_DILU_DDM_FP64_GPU(&
               CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
               method%ddmMat%mat%C_val_dia, &
               method%ddmMat%C_DDMMat_DILU%Dinv,offset_zero, &
               r%C_vec_FP64, r%offset, &
               s%C_vec_FP64, s%offset, &
               method%ddmMat%C_DDMMat_DILU%rvec,offset_zero, &
               method%ddmMat%mat%C_ioff_dia, &
               method%ddmMat%C_DDMMat_DILU%offset_xyz, &
               method%ddmMat%C_DDMMat_DILU%nxe, &
               method%ddmMat%C_DDMMat_DILU%nye, &
               method%ddmMat%C_DDMMat_DILU%nze, &
               method%ddmMat%prm%nx, &
               method%ddmMat%prm%ny, &
               method%ddmMat%prm%nz, &
               method%ddmMat%prm%stmx, &
               method%ddmMat%prm%stmy, &
               method%ddmMat%prm%stmz, &
               method%ddmMat%prm%mx, &
               method%ddmMat%prm%my, &
               method%ddmMat%prm%mz, &
               method%ddmMat%C_DDMMat_DILU%nxblock, &
               method%ddmMat%C_DDMMat_DILU%nyblock, &
               method%ddmMat%C_DDMMat_DILU%nzblock, &
               method%ddmMat%C_DDMMat_DILU%nx_local, &
               method%ddmMat%C_DDMMat_DILU%ny_local, &
               method%ddmMat%C_DDMMat_DILU%nz_local, &
               method%ddmMat%C_DDMMat_DILU%mx_local, &
               method%ddmMat%C_DDMMat_DILU%my_local, &
               method%ddmMat%C_DDMMat_DILU%mz_local, &
               method%ddmMat%nnd)
          call parcel_sync(method)
       ! endif
       ! call MPI_Barrier(MPI_COMM_WORLD,ierr)
       ! call system("sleep 5")
       ! stop

       case(22)
          ! write(*,*)"DDM AS precon "
          call parcel_sync(method)
          call solve_DILU_AS_DDM_GPU(method,r,s)
          call parcel_sync(method)
          return
       end select

    else
       select case (iflag)

!       case(21)
       case(2)
          call solve_DILU_DDM( &
               method,r%vec_FP64(1+r%offset),s%vec_FP64(1+s%offset))
!       case(22)
       case(3)
          call solve_DILU_AS_DDM( &
               method,r%vec_FP64(1+r%offset),s%vec_FP64(1+s%offset))

       case default
          call calc_precon_DDM( &
               method%ddmMat,r%vec_FP64(1+r%offset),s%vec_FP64(1+s%offset))
       end select

       ! if(method%ipreflag.eq.21)then
       !    call solve_DILU_DDM( &
       !         method,r%vec_FP64(1+r%offset),s%vec_FP64(1+s%offset))
       ! else
       !    call calc_precon_DDM( &
       !         method%ddmMat,r%vec_FP64(1+r%offset),s%vec_FP64(1+s%offset))
       ! endif
    endif

  end subroutine precon_ddm


!  subroutine calc_precon(method,r,s,ioffet_r,ioffet_s)
  subroutine calc_precon(method,r,s)
    use krylov
    implicit none
    type(KSP) :: method
    type(PARCEL_VEC) :: r
    type(PARCEL_VEC) :: s
    integer ioffet_r
    integer ioffet_s
    integer MatFormat 

    MatFormat = method%prm%MatFormat 
    select case (MatFormat)
    case(MATRIX_FORMAT_CRS)
!       call precon_crs(method,r,s,ioffet_r,ioffet_s)
       call precon_crs(method,r,s)
    case(MATRIX_FORMAT_DIA)
!       call precon_dia(method,r,s,ioffet_r,ioffet_s)
       call precon_dia(method,r,s)
    case(MATRIX_FORMAT_DDM)
!       call calc_precon_DDM(method%ddmMat,r,s)
       call precon_DDM(method,r,s)
    end select

  end subroutine calc_precon

end module PARCEL_KRYLOV_PRECON

