module QR_kernel
  use TSQR_kernel
contains

  subroutine MGS_Orthonormalization_CPU(&
       n,s,i,icomm,&
       V, ldv,&
       Q, ldq,&
       R, ldr )
    use vec_kernel,only : calc_dot,axpy
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    integer(kind=intType) i
    integer icomm
    integer(kind=intType) ldv
    integer(kind=intType) ldq
    integer(kind=intType) ldr
    real*8 V(ldv)   ! (n)
    real*8 Q(ldq*s) ! (n*s)
    real*8 R(ldr*s)   ! (s*s)
    real*8 tmp1,tmp2
    integer(kind=intType) k

    ! -----
    ! 修正グラムシュミット
    do k=1,i-1
       ! Rki = <Qk,Vi>
       call calc_dot(&
            n,icomm ,&
            V,&
            Q(1+ldq*(k-1))   , &
            R(k+ldr*(i-1)) )

       ! Vi=Vi-Rki*Qk
       tmp1=-R(k+ldr*(i-1))
       tmp2=1.0d0
       call axpy(n,tmp1,tmp2,Q(1+ldq*(k-1)),V)
    enddo

  end subroutine MGS_Orthonormalization_CPU

#if 1
!  subroutine MGS_Orthonormalization_GPU(&
  subroutine MGS_Orthonormalization(& 
       method, &
       n,s,i,icomm,&
       V, ldv,&
       Q, ldq,&
       R, ldr)
    use krylov,only : KSP
    use vec_type,only : PARCEL_VEC
    use PARCEL_VEC_kernel,only : &
         set_val_FP64, get_val_FP64, calc_dot,axpy,init_vec,free_vec
    use PARCEL_integer
    use PARCEL_VEC_kernel,only : parcel_sync
    implicit none
    type(KSP) :: method
    integer(kind=intType) n
    integer(kind=intType) s
    integer(kind=intType) i
    integer icomm
    type(PARCEL_VEC)  :: V
    integer(kind=intType) ldv
    type(PARCEL_VEC)  :: Q
    integer(kind=intType) ldq
    type(PARCEL_VEC)  :: R
    integer(kind=intType) ldr
    integer(kind=intType) k

    type(PARCEL_VEC) :: tmp1
    type(PARCEL_VEC) :: tmp2
    real*8 tmp1_FP64
    real*8 tmp2_FP64

    integer(kind=intType) offset_Q0
    integer(kind=intType) offset_V0
    integer(kind=intType) parcel_int

    offset_Q0 = Q%offset 
    offset_V0 = V%offset 

    parcel_int = 1
    call init_vec(method,tmp1,parcel_int)
    call init_vec(method,tmp2,parcel_int)
    call set_val_FP64(method,tmp1,0.0d0,parcel_int)
    call set_val_FP64(method,tmp2,1.0d0,parcel_int)
    ! -----
    ! 修正グラムシュミット
    do k=1,i-1
       ! Rki = <Qk,Vi>
       Q%offset = offset_Q0 + ldq*(k-1)
       call parcel_sync(method)
       call calc_dot(&
            method, &
            icomm ,&
            V,&
            Q, &
            tmp1 ,&
            n)
       Q%offset = offset_Q0 
       call parcel_sync(method)
       parcel_int = 1
       call get_val_FP64(method,tmp1,tmp1_FP64,parcel_int)
       call parcel_sync(method)
       ! write(*,*)"@ 1 @ offset_V = ",offset_V0
       ! write(*,*)"@ 1 @ offset_Q = ",offset_Q0 + ldq*(k-1)
       ! write(*,*)"@ 1 @ R(",k,") = ",tmp1_FP64,1
       parcel_int = k+s*(i-1)
       call set_val_FP64(method,R   ,tmp1_FP64,parcel_int)

       ! Vi=Vi-Rki*Qk
       tmp1_FP64 = -tmp1_FP64
       parcel_int = 1
       call set_val_FP64(method,tmp1,tmp1_FP64,parcel_int)
       call set_val_FP64(method,tmp2,1.0d0    ,parcel_int)

       Q%offset = offset_Q0 + ldq*(k-1)
       call axpy(method,tmp1,tmp2,Q,V,n)
       call parcel_sync(method)
       Q%offset = offset_Q0 
    enddo

    call free_vec(method,tmp1)
    call free_vec(method,tmp2)

    Q%offset = offset_Q0
    V%offset = offset_V0

  end subroutine MGS_Orthonormalization
#else
  subroutine MGS_Orthonormalization(&
       method,&
       n,s,i,icomm,&
       V, ldv, &
       Q, ldq,&
       R, ldr)
    use krylov,only : KSP
    use vec_type,only : PARCEL_VEC
    use KSP_STATUS,only : STATUS_FLAG,PARCEL_CPU_FLAG_ON,PARCEL_GPU_FLAG_ON
    use PARCEL_integer
    implicit none
    type(KSP) :: method
    integer(kind=intType) n
    integer(kind=intType) s
    integer(kind=intType) i
    integer icomm
    type(PARCEL_VEC)  :: V
    integer(kind=intType) ldv
    type(PARCEL_VEC)  :: Q
    integer(kind=intType) ldq
    type(PARCEL_VEC)  :: R
    integer(kind=intType) ldr
    type(STATUS_FLAG) :: flag
    integer myrank

    flag = method%flag
    select case(flag%machineflag)
    case(PARCEL_CPU_FLAG_ON)
       myrank = method%commVec%myrank
       call MGS_Orthonormalization_CPU(&
            n,s,i,icomm,&
            V%vec_FP64(1+V%offset), ldv,&
            Q%vec_FP64(1+Q%offset), ldq,&
            R%vec_FP64(1+R%offset), ldr )
    case(PARCEL_GPU_FLAG_ON)
       call MGS_Orthonormalization_GPU(&
            n,s,i,icomm,&
            V, ldv,&
            Q, ldq,&
            R, ldr )
    end select

  end subroutine MGS_Orthonormalization
#endif
  subroutine CGS_Orthonormalization_CPU(&
       n,s,i,icomm,&
       V, ldv,&
       Q, ldq,&
       R, ldr,&
       nvecBlock)
    use vec_kernel,only : calc_gemm_dot_CPU
    use vec_kernel_CPU,only : calc_dot_FP64_CPU,calc_dot_local_FP64_CPU
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    integer(kind=intType) i
    integer icomm
    integer(kind=intType) ldv
    integer(kind=intType) ldq
    integer(kind=intType) ldr
    real*8 V(ldv) ! (n)
    real*8 Q(ldq*s) ! (n*s)
    real*8 R(ldr*s) ! (s*s)
    real*8 buf(i) ! (s*s)
    integer nvecBlock

    real*8 tmp1,tmp2
    integer(kind=intType) k
    integer(kind=intType) j,jj
    integer(kind=intType) nblock
    integer ierr

    ! -----
    ! 古典グラムシュミット
    ! Rki = <Qk,Vi>
    if(i>0)then
       j = 1
       call calc_gemm_dot_CPU(&
            n,i,j,&
            icomm,&
            Q, ldq,&
            V, ldv,&
            R(1+s*(i-1)), j, &
            nvecBlock)
    endif

    !$OMP PARALLEL DO PRIVATE(jj,nblock,k,j,tmp1,tmp2) 
    do jj=1,n,nvecBlock 
       nblock=min(n,jj+nvecBlock-1)-jj+1
       do k=1,i-1
          ! Vi=Vi-Rki*Qk
          !           tmp1=-R(k+s*(i-1))
          !           tmp2=1.0d0
          !           call axpy(nblock,tmp1,tmp2,Q(jj+n*(k-1)),V(jj) )
          do j=jj,jj+nblock-1
!             V(j)=V(j)-R(k+s*(i-1))*Q(j+n*(k-1))
             V(j)=V(j)-R(k+ldr*(i-1))*Q(j+ldq*(k-1))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO 
    
  end subroutine CGS_Orthonormalization_CPU

   subroutine CGS_Orthonormalization_GPU(&
       method,&
       n,s,i,icomm,&
       V, ldv,&
       Q, ldq,&
       R, ldr,&
       nvecBlock)
    use krylov,only : KSP
    use vec_type,only : PARCEL_VEC
    use PARCEL_VEC_kernel,only : calc_gemm_dot,calc_dot
    use cu_Vec,only : CGS_Orthonormalization_FP64_GPU
    use PARCEL_integer,only : intType

    use cu_common_wrapper,only : cudaMemcpy_d2h_FP64

    implicit none
    type(KSP) :: method
    integer(kind=intType) n
    integer(kind=intType) s
    integer(kind=intType) i
    integer icomm
    type(PARCEL_VEC)  :: V
    integer(kind=intType) ldv
    type(PARCEL_VEC)  :: Q
    integer(kind=intType) ldq
    type(PARCEL_VEC)  :: R
    integer(kind=intType) ldr
    integer nvecBlock

    real*8 val_(i)
    integer(kind=intType) k
    integer(kind=intType) j,jj
    integer(kind=intType) nblock
    integer(kind=intType) offset0_R
    integer(kind=intType) offset0_Q

    ! -----
    ! 古典グラムシュミット
    ! Rki = <Qk,Vi>

    offset0_R = R%offset 
    offset0_Q = Q%offset 
    if(i>0)then
       j = 1
!       R%offset = offset0_R + s*(i-1)
#if 0
       do k = 1,i
          R%offset = offset0_R + ldr*(i-1) + (k-1)
          Q%offset = offset0_Q + ldq*(k-1)
          call calc_dot( &
               method,&
               icomm,Q,V,R, n)       
          Q%offset = offset0_Q
          R%offset = offset0_R
       enddo
#else
       R%offset = offset0_R + ldr*(i-1)
       call calc_gemm_dot(&
            method,&
            icomm,&
            n,i,j,&
            Q, ldq,&
            V, ldv,&
            R, ldr,&
            nvecBlock)
#endif
       ! val_=0.0d0
       ! call cudaMemcpy_d2h_FP64(R%C_vec_FP64,val_,i,R%offset,0)       
       ! do k = 1,i
       !    write(*,*)" R = ",k,val_(k)
       ! enddo
       R%offset = offset0_R

    endif

    call CGS_Orthonormalization_FP64_GPU( &
         V%C_vec_FP64, V%offset, ldv,&
         Q%C_vec_FP64, Q%offset, ldq,&
         R%C_vec_FP64, R%offset, ldr,&
         n,s,i &
         )
    
  end subroutine CGS_Orthonormalization_GPU

  subroutine CGS_Orthonormalization( &
       method,&
       n,s,i,&
       icomm,&
       V,ldv,&
       Q,ldq,&
       R,ldr,&
       nvecBlock)
    use krylov,only : KSP
    use vec_type,only : PARCEL_VEC
    use KSP_STATUS,only : STATUS_FLAG,PARCEL_CPU_FLAG_ON,PARCEL_GPU_FLAG_ON
    use PARCEL_integer
    implicit none
    type(KSP) :: method
    integer(kind=intType) n
    integer(kind=intType) s
    integer(kind=intType) i
    integer icomm
    type(PARCEL_VEC)  :: V
    integer(kind=intType) ldv
    type(PARCEL_VEC)  :: Q
    integer(kind=intType) ldq
    type(PARCEL_VEC)  :: R
    integer(kind=intType) ldr
    integer nvecBlock
    type(STATUS_FLAG) :: flag

    flag = method%flag
    select case(flag%machineflag)
    case(PARCEL_CPU_FLAG_ON)
       call CGS_Orthonormalization_CPU(&
            n,s,i,icomm,&
            V%vec_FP64(1+V%offset), ldv,&
            Q%vec_FP64(1+Q%offset), ldq,&
            R%vec_FP64(1+R%offset), ldr,&
            nvecBlock)
    case(PARCEL_GPU_FLAG_ON)
!       write(*,*)"CGS START"
       call CGS_Orthonormalization_GPU(&
            method, &
            n,s,i,icomm,&
            V, ldv,&
            Q, ldq,&
            R, ldr,&
            nvecBlock)
!       write(*,*)"CGS END"
    end select

  end subroutine CGS_Orthonormalization

  subroutine GS_QR(n,s,icomm,V,Q,R,iGSflag,nvecBlock)
    use PARCEL_define
    use vec_kernel
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    integer icomm
    real*8 V(n*s)
    real*8 Vtmp(n)
    real*8 Q(n*s)
    real*8 R(s*s)
    integer iGSflag
    integer nvecBlock

    integer(kind=intType) i,j
    real*8 tmp1,tmp2

    do i=1,s*s
       R(i)=0.0d0
    enddo

    ! todo 古典グラムは内積よりも行列積で実装したほうが高速
    do i=1,s
       if(iGSflag.eq.KSP_QRDECOMPOSITION_MGS)then
          call MGS_Orthonormalization_CPU(&
               n,s,i,icomm,&
               V( 1+n*(i-1)),n, &
               Q,n, &
               R,s)
       else
          call CGS_Orthonormalization_CPU(&
               n,s,i,icomm,&
               V(1+n*(i-1)),n,&
               Q,n,&
               R,s,nvecBlock)
       endif
       call calc_dot(n,icomm,V(1+n*(i-1)),V(1+n*(i-1)),R(i+s*(i-1)))
       R(i+s*(i-1))=sqrt(R(i+s*(i-1)))
       tmp1=1.0d0/R(i+s*(i-1))
       tmp2=0.0d0
       call axpy(n,tmp1,tmp2,V(1+n*(i-1)),Q(1+n*(i-1)))
    enddo

  end subroutine GS_QR

  subroutine GS_QR_( &
       method, &
       n, &
       s, &
       icomm, &
       V, ldv, &
       Q, ldq, &
       R, ldr, &
       iGSflag, &
       nvecBlock)
    use krylov,only : KSP
    use PARCEL_define,only : KSP_QRDECOMPOSITION_MGS
    use PARCEL_integer,only : intType
    use PARCEL_VEC_kernel,only : &
         init_vec,free_vec,set_vec_FP64, &
         calc_dot,xcpy,axpy,sqrt_val, &
         set_val_FP64,set_vec_val,get_val_FP64

    use vec_type,only : PARCEL_VEC

    implicit none
    type(KSP) :: method
    integer(kind=intType) n
    integer(kind=intType) s
    integer icomm

    ! real*8 V(n*s)
    ! real*8 Q(n*s)
    ! real*8 R(s*s)

    type(PARCEL_VEC) :: V
    integer(kind=intType) ldv
    type(PARCEL_VEC) :: Q
    integer(kind=intType) ldq
    type(PARCEL_VEC) :: R
    integer(kind=intType) ldr

    ! real*8 V(n*s)
    ! real*8 Q(n*s)
    ! real*8 R(s*s)

    integer iGSflag
    integer nvecBlock

    integer(kind=intType) i,j
    type(PARCEL_VEC) ::  tmp1
    real*8 tmp1_FP64
    type(PARCEL_VEC) ::  tmp2

    integer(kind=intType) offset0_V
    integer(kind=intType) offset0_Q
    integer(kind=intType) offset0_R
    integer(kind=intType) parcel_int

    parcel_int = 1
    call init_vec(method,tmp1 ,parcel_int)
    call init_vec(method,tmp2 ,parcel_int)
    call set_val_FP64(method,tmp1 ,0.0d0,parcel_int)
    call set_val_FP64(method,tmp2 ,0.0d0,parcel_int)
    call set_vec_val(method,R,tmp1,parcel_int,ldr*s)

    offset0_V = V%offset
    offset0_Q = Q%offset
    offset0_R = R%offset

    do i=1,s

       if(iGSflag.eq.KSP_QRDECOMPOSITION_MGS)then
          V%offset = offset0_V + (i-1)*ldv
          call MGS_Orthonormalization(&
               method,n,s,i,icomm,&
               V, ldv, &
               Q, ldq,&
               R, ldr)
          V%offset = offset0_V
       else
          V%offset = offset0_V + (i-1)*ldv
          call CGS_Orthonormalization( &
               method,n,s,i,icomm,&
               V, ldv,&
               Q, ldq,&
               R, ldr,&
               nvecBlock)
          V%offset = offset0_V
       endif

       V%offset = offset0_V +   (i-1)*ldv
       R%offset = offset0_R + i+(i-1)*ldr - 1
       call calc_dot( &
            method,&
            icomm,V,V,R, n)
       V%offset = offset0_V
       R%offset = offset0_R
       call sqrt_val(method,R, i+(i-1)*ldr )

       call get_val_FP64(method,R,tmp1_FP64, i+(i-1)*ldr)
       tmp1_FP64 = 1.0d0/tmp1_FP64
       !       write(*,*)"@ 2 @ R(",1+i+(i-1)*ldr,") = ",tmp1_FP64
       parcel_int = 1
       call set_val_FP64(method,tmp1,tmp1_FP64,parcel_int)
       call set_val_FP64(method,tmp2,0.0d0,parcel_int)

       V%offset = offset0_V + ldv*(i-1)
       Q%offset = offset0_Q + ldq*(i-1)
       call axpy(method,tmp1,tmp2,V,Q,n)
       V%offset = offset0_V 
       Q%offset = offset0_Q 

#if 0
       Q%offset = offset0_Q + ldq*(i-1)
       call calc_dot( &
            method,&
            icomm,Q,Q,tmp1, n)
       Q%offset = offset0_Q 
       call get_val_FP64(method,tmp1,tmp1_FP64,1)
       write(*,*)"<Q|Q>",i,tmp1_FP64

       if(i.ne.s)then
          V%offset = offset0_V + ldv*i
          Q%offset = offset0_Q + ldq*(i-1)
          call calc_dot( &
               method,&
               icomm,V,Q,tmp1, n)
          Q%offset = offset0_Q 
          V%offset = offset0_V 
          call get_val_FP64(method,tmp1,tmp1_FP64,1)

          write(*,*)"offset_V",i,offset0_V + ldv*i
          write(*,*)"offset_Q",i,offset0_Q + ldq*(i-1)
          write(*,*)"<Q|V>   ",i,tmp1_FP64
       endif
#endif

    enddo

    call free_vec(method,tmp1)
    call free_vec(method,tmp2)

    V%offset = offset0_V
    Q%offset = offset0_Q
    R%offset = offset0_R

  end subroutine GS_QR_

  subroutine calcQ_VRinv(n,s,V,Q,R,nvecBlock)
    use vec_kernel ! :only,nvecBlock
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    real*8 V(n*s)
    real*8 Q(n*s)
    real*8 R(s*s)
    integer nvecBlock

    real*8 Rinv(s*s)
    integer(kind=intType) i,j,k,jj
    integer(kind=intType) nblock
    integer ierr
    integer s_int4
    Rinv=R
    ! 三角行列の逆行列を計算 (三角行列の逆行列は三角行列)
    s_int4=s
    call DTRTRI('U','N',s_int4,Rinv,s_int4,ierr)

    !$OMP PARALLEL DO PRIVATE(jj,nblock,i,j,k) collapse(2)
    do jj=1,n,nvecBlock 
       do i=1,s
          nblock=min(n,jj+nvecBlock-1)
          do j=jj,nblock
             Q(j+(i-1)*n)=0.0d0
          enddo
          do k=1,i
             do j=jj,nblock
                Q(j+(i-1)*n)= &
                     Q(j+(i-1)*n) + V(j+(k-1)*n) * Rinv(k+(i-1)*s)
             enddo
          enddo

       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine calcQ_VRinv

  subroutine calcQ_VRinv_(&
       method, &
       n, &
       s, &
       V, ldv,&
       Q, ldq,&
       R, ldr,&
       nvecBlock)
    use krylov,only : KSP
    use PARCEL_integer,only : intType
    use c_wrapper,only : devicemem
    use vec_type,only : PARCEL_VEC
    use PARCEL_VEC_kernel,only : &
         set_vec_FP64,get_vec_FP64
    use cu_Vec,only : calcQ_VRinv_FP64_GPU
    use KSP_STATUS,only : STATUS_FLAG,PARCEL_CPU_FLAG_ON,PARCEL_GPU_FLAG_ON
    implicit none

    type(KSP) :: method
    integer(kind=intType) n
    integer(kind=intType) s
    type(PARCEL_VEC) ::  V
    integer(kind=intType) ldv
    type(PARCEL_VEC) ::  Q
    integer(kind=intType) ldq
    type(PARCEL_VEC) ::  R
    integer(kind=intType) ldr
    integer nvecBlock

    real*8 Rinv(s*s)
    real*8 R_(s*s)
    integer(kind=intType) i,j,k,jj
    integer(kind=intType) nblock
    integer ierr
    integer s_int4

    type(STATUS_FLAG) :: flag
    flag = method%flag
    select case(flag%machineflag)
    case(PARCEL_CPU_FLAG_ON)
       call calcQ_VRinv( &
            n, &
            s, &
            V%vec_FP64(1+V%offset), &
            Q%vec_FP64(1+Q%offset), &
            R%vec_FP64(1+R%offset), &
            nvecBlock)
    case(PARCEL_GPU_FLAG_ON)
       
       call get_vec_FP64(method,R,Rinv,s*s)
       R_ = Rinv
       
       ! 三角行列の逆行列を計算 (三角行列の逆行列は三角行列)
       s_int4=s
       call DTRTRI('U','N',s_int4,Rinv,s_int4,ierr)
       call set_vec_FP64(method,R,Rinv,s*s)
       
       call calcQ_VRinv_FP64_GPU( &
            s, &
            n, &
            V%C_vec_FP64 , V%offset, ldv, &
            Q%C_vec_FP64 , Q%offset, ldq, &
            R%C_vec_FP64 , R%offset, ldr  &
            )
       call set_vec_FP64(method,R,R_,s*s)
       
    end select

  end subroutine calcQ_VRinv_

  subroutine CholeskyQR_Ronly(n,s,icomm,V,R,nvecBlock)
    use mpi
    use vec_kernel
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    integer icomm
    real*8 V(n*s)
    real*8 R(s*s)
    integer nvecBlock

    integer(kind=intType) i,j
    integer ierr
    real*8 B(s*s)
    integer s_int4
    integer ss_int4
    do i=1,s*s
       R(i)=0.0d0
       B(i)=0.0d0
    enddo
#if 0
    ! DSYRK
    ! 内積形式なので行列積形式に修正する 
    !   subroutine calc_gemm_dot(n,s,t,icomm,x,y,dotp) で良い? --> 対称性を考慮したほうが良い
    do i=1,s
       do j=1,i
          call calc_dot_local(n,V(1+(i-1)*n),V(1+(j-1)*n),B(j+s*(i-1)))
       enddo
    enddo
#else
    call calc_syrk_dot_local(n,s,V,B,nvecBlock)
#endif
    ! ALLREDUCE
    ! 三角行列なので(s+1)*s/2のデータ数を送受信すれば良い
    ss_int4=s*s
    call MPI_Allreduce(B,R, ss_int4,MPI_REAL8,MPI_SUM,icomm,ierr)    

    ! コレスキー分解
    s_int4=s
    call DPOTRF('U', s_int4,R, s_int4,ierr)

  end subroutine CholeskyQR_Ronly

  subroutine CholeskyQR_Ronly_( &
       method, &
       n, &
       s, &
       localnRow, &
       icomm, &
       batchCount, &
       V_batch_prt, &
       work_batch_prt, &
       V, ldv, &
       R, ldr, &
       work, &
       nvecBlock)
    use mpi
    use krylov,only : KSP
    use PARCEL_integer,only : intType
    use c_wrapper,only : devicemem
    use vec_type,only : PARCEL_VEC
    use cu_common_wrapper,only : &
         cublas_DgemmBatched_gemm_dot_local_FP64
    use PARCEL_VEC_kernel,only : &
         set_vec_FP64,get_vec_FP64
    use KSP_STATUS,only : STATUS_FLAG,PARCEL_CPU_FLAG_ON,PARCEL_GPU_FLAG_ON
    use vec_kernel,only : calc_syrk_dot_local

    implicit none
    type(KSP) :: method
    integer(kind=intType) n
    integer(kind=intType) s
    integer(kind=intType) localnRow
    integer icomm
    integer(kind=intType) batchCount
    type(devicemem) V_batch_prt
    type(devicemem) work_batch_prt
    integer(kind=intType) parcel_int

    !    real*8 V(n*s)
    !    real*8 R(s*s)
    type(PARCEL_VEC) :: V
    integer(kind=intType) ldv

    type(PARCEL_VEC) :: R
    integer(kind=intType) ldr

    type(PARCEL_VEC) :: work
    integer nvecBlock

    integer(kind=intType) i,j
    integer ierr
    real*8 B(s*s)
    real*8 R_FP64(s*s)
    integer s_int4
    integer ss_int4
    type(STATUS_FLAG) :: flag
    integer irank

    !    call init_vec(method, R, s*ldr)

    flag = method%flag
    select case(flag%machineflag)

    case(PARCEL_CPU_FLAG_ON)
       ! CPU
       do i=1,s*s
          B(i)      = 0.0d0
          R_FP64(i) = 0.0d0
       enddo

       call calc_syrk_dot_local( &
            n, &
            s, &
            V%vec_FP64(1+V%offset), &
            B, &
            nvecBlock)

       ! ALLREDUCE
       ! 三角行列なので(s+1)*s/2のデータ数を送受信すれば良い
       ss_int4=s*s
       call MPI_Allreduce(B,R_FP64, ss_int4,MPI_REAL8,MPI_SUM,icomm,ierr)    

       ! コレスキー分解
       s_int4=s
       call DPOTRF('U', s_int4,R_FP64, s_int4,ierr)

       call set_vec_FP64(method,R,R_FP64,s*s)

    case(PARCEL_GPU_FLAG_ON)
       ! GPU
       do i=1,s*s
          B(i)      = 0.0d0
          R_FP64(i) = 0.0d0
       enddo
       parcel_int = batchCount
       call cublas_DgemmBatched_gemm_dot_local_FP64( &
            method%cublas_solver, &
            parcel_int, &
            s, localnRow, &
            V_batch_prt , ldv, &
            V_batch_prt , ldv, &
            work_batch_prt, s*batchCount, &
            work%C_vec_FP64,  &
            R%C_vec_FP64 ) 
       call get_vec_FP64(method,R,B,s*s)

       ! ALLREDUCE
       ! 三角行列なので(s+1)*s/2のデータ数を送受信すれば良い
       ss_int4=s*s
       call MPI_Allreduce(B,R_FP64, ss_int4,MPI_REAL8,MPI_SUM,icomm,ierr)    

       ! コレスキー分解
       s_int4=s
       call DPOTRF('U', s_int4,R_FP64, s_int4,ierr)

#if 1
       do j = 1,s
          do i = j+1,s
             R_FP64(i+(j-1)*s)=0.0d0
          enddo
       enddo
#endif

       call set_vec_FP64(method,R,R_FP64,s*s)
    end select

#if 0
    call MPI_Comm_rank(icomm,irank, ierr)
    do j=1,s
       do i=1,s
          if(irank.eq.0)then
             write(*,*)"Rtri : ",i,j,R_FP64(j+(i-1)*(s))
          endif
       enddo
    enddo
#endif

  end subroutine CholeskyQR_Ronly_

  subroutine TSQR_Ronly(n,s,icomm,V,R,nvecBlock)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    integer icomm
    real*8 V(n*s)
    real*8 R(s*s)
    integer nvecBlock

    integer(kind=intType) i,j
    call TSQR_local_Ronly(n,s,V,R,nvecBlock)
    call TSQR_GLOBAL(icomm,s,R)

    ! 対角成分が正の数の三角行列を返す
    ! (元の三角行列とは直交行列の列ベクトルの方向が変化する)
    ! -----
    ! CAGMRESは三角行列の対角成分が正の数であることを仮定しているため三角行列の対角成分が負の数である場合は符号を反転させなければならない。
    ! QR分解V=QRの場合、対角成分が正の三角行列R'をR=ER'と定義するとQR分解の式はV=QER'となる。
    ! 行列Eを三角行列Rの対角成分が負の数の行のみ-1でそれ以外が1の対角行列と定義するとQR分解の式はV=Q'R'となる。
    ! 新しい行列Q'=QEであるため元の直行行列QとQ'の関係は列ベクトルの方向が異なるだけであり、Q'も直行行列である。
    ! 従って三角行列の対角成分が全て正の数の三角行列R'をQR分解の三角行列として採用しても良い。

    do j=1,s
       if(R(j+(j-1)*(s)) < 0.0d0)then
          do i=j,s
             R(j+(i-1)*(s))=-R(j+(i-1)*(s))
          enddo
       endif
    enddo

  end subroutine TSQR_Ronly

  subroutine TSQR_Ronly_(&
       method, &
       n,s,icomm,&
       V, ldv, &
       R, ldr, &
       work, &
       nvecBlock)
    use krylov,only : KSP
    use vec_type,only : PARCEL_VEC
    use PARCEL_integer
    use PARCEL_define,only : WARP_SIZE
    use c_wrapper,only : CUDA_NUM_BLOCKS
    use PARCEL_VEC_kernel,only : parcel_sync,set_vec_FP64, get_vec_FP64
    ! use PARCEL_VEC_kernel,only : &
    !      set_vec_FP64, get_vec_FP64, init_vec,free_vec

    implicit none
    type(KSP) :: method
    integer(kind=intType) n
    integer(kind=intType) s
    integer icomm
    ! real*8 V(n*s)
    ! real*8 R(s*s)
    type(PARCEL_VEC)  :: V
    integer(kind=intType) ldv
    type(PARCEL_VEC)  :: R
    integer(kind=intType) ldr
    type(PARCEL_VEC)  :: work
    integer nvecBlock

    real*8 R_FP64(s*s)
    integer(kind=intType) i,j
    integer numblock 
    integer irank,ierr

#if 1

    ! 必要なワークサイズ : double * (s*s+s*(s+WARP_SIZE))*(numblock)
    ! call init_vec(method,work,(s*s+s*(s+WARP_SIZE))*(numblock))
    numblock = CUDA_NUM_BLOCKS
    call TSQR_local_Ronly_( &
       method, &
       n,s,&
       V, ldv, &
       R, ldr, &
       work, &
       numblock, &
       nvecBlock)
    ! call free_vec(method,work)
    call parcel_sync(method)
    call get_vec_FP64(method,R,R_FP64,s*s)
    call parcel_sync(method)
    call TSQR_GLOBAL(icomm,s,R_FP64)

    do j=1,s
       if(R_FP64(j+(j-1)*(s)) < 0.0d0)then
          do i=j,s
             R_FP64(j+(i-1)*(s))=-R_FP64(j+(i-1)*(s))
          enddo
       endif
    enddo
    call parcel_sync(method)
    call set_vec_FP64(method,R,R_FP64,s*s)
    call parcel_sync(method)

#if 0
    call MPI_Comm_rank(icomm,irank, ierr)
    do j=1,s
       do i=1,s
          if(irank.eq.0)then
             write(*,*)"Rtri : ",i,j,R_FP64(j+(i-1)*(s))
          endif
       enddo
    enddo
#endif

#else
    ! 1プロセスあたり64次元以上ないと正常に計算できない?
    R_FP64 = 0.0d0
    do j=1,s
       V%offset = (j-1)*ldv
       call get_vec_FP64(method,V,V_FP64(1+(j-1)*n),n)
    enddo
    V%offset = 0

! #if 0
!     call TSQR_Ronly(n,s,icomm,V_FP64,R_FP64,nvecBlock)
! #endif
#if 1
    call TSQR_local_Ronly(n,s,V_FP64,R_FP64,nvecBlock)
    call TSQR_GLOBAL(icomm,s,R_FP64)

    do j=1,s
       if(R_FP64(j+(j-1)*(s)) < 0.0d0)then
          do i=j,s
             R_FP64(j+(i-1)*(s))=-R_FP64(j+(i-1)*(s))
          enddo
       endif
    enddo
#endif
    call set_vec_FP64(method,R,R_FP64,s*s)

#if 0
    call MPI_Comm_rank(icomm,irank, ierr)
    do j=1,s
       do i=1,s
          if(irank.eq.0)then
             write(*,*)"Rtri : ",i,j,R_FP64(j+(i-1)*(s))
          endif
       enddo
    enddo
#endif

#endif

  end subroutine TSQR_Ronly_

  subroutine CholeskyQR2_Ronly(n,s,icomm,V,R,nvecBlock)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    integer icomm
    real*8 V(n*s)
    real*8 R(s*s)
    integer nvecBlock

    real*8 Q1(n*s)
    real*8 R2(s*s)
    integer s_int4

    call CholeskyQR_Ronly(n,s,icomm,V,R,nvecBlock)
    call calcQ_VRinv(n,s,V,Q1,R,nvecBlock)
    call CholeskyQR_Ronly(n,s,icomm,Q1,R2,nvecBlock)

    s_int4=s
    call DTRMM('L','U','N','N', s_int4, s_int4,1.0d0,R2,s_int4,R,s_int4)

  end subroutine CholeskyQR2_Ronly

  subroutine CholeskyQR2_Ronly_(&
       method, &
       n, &
       s, &
       localnRow, &
       icomm, &
       batchCount, &
       V_batch_prt, &
       work_batch_prt, &
       V, ldv, &
       R, ldr, &
       work, &
       nvecBlock)
    use krylov,only : KSP
    use PARCEL_integer,only : intType
    use c_wrapper,only : devicemem
    use vec_type,only : PARCEL_VEC
    use PARCEL_VEC_kernel,only : &
         set_vec_FP64, get_vec_FP64, init_vec,free_vec, xcpy

    implicit none

    type(KSP) :: method
    integer(kind=intType) n
    integer(kind=intType) s
    integer(kind=intType) localnRow
    integer icomm
    integer(kind=intType) batchCount
    type(devicemem) V_batch_prt
    type(devicemem) work_batch_prt

    type(PARCEL_VEC) V ! (n*s)
    integer(kind=intType) ldv
    type(PARCEL_VEC) R ! (s*s)
    integer(kind=intType) ldr

    type(PARCEL_VEC) :: work

    integer nvecBlock
    type(PARCEL_VEC) V0 ! (n*s)
    type(PARCEL_VEC) Q1 ! (n*s)
    type(PARCEL_VEC) R2 ! (s*s)

    real*8 R_FP64(s*s)
    real*8 R2_FP64(s*s)
    integer s_int4

    call init_vec(method,Q1,s*ldv)
    call init_vec(method,V0,s*ldv)
    call init_vec(method,R2,s*ldr)
    call xcpy(method,V,V0,s*ldv)

    call CholeskyQR_Ronly_( &
         method, &
         n, &
         s, &
         localnRow, &
         icomm, &
         batchCount, &
         V_batch_prt, &
         work_batch_prt, &
         V, ldv,&
         R, ldr,&
         work, &
         nvecBlock)

    call calcQ_VRinv_(&
         method, &
         n, &
         s, &
         V, ldv,&
         Q1, ldv,&
         R, ldr,&
         nvecBlock)
    call xcpy(method,Q1,V,s*ldv)

    call CholeskyQR_Ronly_(&
         method, &
         n, &
         s, &
         localnRow, &
         icomm, &
         batchCount, &
         V_batch_prt, &
         work_batch_prt, &
         !            Q1, &
         V, ldv, &
         R2, ldr,&
         work, &
         nvecBlock)

    call get_vec_FP64(method,R,R_FP64,s*s)
    call get_vec_FP64(method,R2,R2_FP64,s*s)

    s_int4=s
    call DTRMM('L','U','N','N', s_int4, s_int4,1.0d0,R2_FP64,s_int4,R_FP64,s_int4)

    call xcpy(method,V0,V,s*ldv)
    call set_vec_FP64(method,R,R_FP64,s*s)
    call set_vec_FP64(method,R2,R2_FP64,s*s)

    call free_vec(method,Q1)
    call free_vec(method,V0)
    call free_vec(method,R2)

  end subroutine CholeskyQR2_Ronly_

  subroutine check_QtQ(n,s,icomm,Q)
    use mpi
    use vec_kernel
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    integer icomm
    real*8 Q(n*s)
    real*8 B(s*s)
    real*8 QtQ(s*s)
    integer irank,ierr
    integer(kind=intType) i,j
    real*8 Orthnorm
    integer(kind=intType) gn
    real*8 n_
    integer ss_int4 
    gn=0
    !    call MPI_Allreduce(n,gn,1,MPI_INTEGER,MPI_SUM,icomm,ierr)    
    call MPI_Allreduce(n,gn,1,intMPIType,MPI_SUM,icomm,ierr)    
    n_=gn*s
    call MPI_Comm_rank(icomm,irank, ierr)

    !     if(irank.eq.0)then
    !        write(*,*)"--------"
    !        write(*,*)"--------"
    !        write(*,*)Q
    !        write(*,*)"--------"
    !        write(*,*)"--------"
    !     endif
    B=0.0d0
    QtQ=0.0d0
    do i=1,s
       do j=1,s
          call calc_dot_local(n,Q(1+(i-1)*n),Q(1+(j-1)*n),B(j+s*(i-1)))
       enddo
    enddo
    ss_int4=s*s
    call MPI_Allreduce(B,QtQ,ss_int4,MPI_REAL8,MPI_SUM,icomm,ierr)    

    if(irank.eq.0)then
       do i=1,s
          QtQ(i+(i-1)*s)=QtQ(i+(i-1)*s)-1.0d0
       enddo

       Orthnorm=0.0d0
       do i=1,s
          do j=1,s
             Orthnorm=Orthnorm+QtQ(j+(i-1)*s)*QtQ(j+(i-1)*s)
          enddo
       enddo
       Orthnorm=sqrt(Orthnorm)
!       write(*,*)"|| QtQ - I ||           = ",Orthnorm
       write(*,*)"|| QtQ - I || / sqrt(n*s) = ", Orthnorm/sqrt(n_)
    endif

  end subroutine check_QtQ
  
  subroutine check_A_QR(n,s,icomm,A,Q,R)
    use mpi
    use vec_kernel
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    integer icomm
    real*8 A(n*s)
    real*8 Q(n*s)
    real*8 R(s*s)
    real*8 QR(n*s)

    integer irank,ierr
    integer(kind=intType) i,j,k
    real*8 Orthnorm
    integer(kind=intType) gn
    real*8 err_norm

    real*8 n_
    real*8 sbuf,rbuf
    integer ss_int4 

    gn=0
    !    call MPI_Allreduce(n,gn,1,MPI_INTEGER,MPI_SUM,icomm,ierr)    
    call MPI_Allreduce(n,gn,1,intMPIType,MPI_SUM,icomm,ierr)    
    n_=gn*s
    call MPI_Comm_rank(icomm,irank, ierr)

    QR=0.0d0
    do k=1,n
       do i=1,s
          do j=1,s
             QR(k+(i-1)*n) = QR(k+(i-1)*n) + Q(k+(j-1)*n) * R(j+(i-1)*s)
          enddo
       enddo
    enddo
    err_norm = 0.0d0
    do k=1,n
       do i=1,s
          err_norm = err_norm + ( A(k+(i-1)*n) -  QR(k+(i-1)*n) ) * ( A(k+(i-1)*n) -  QR(k+(i-1)*n) )
!          write(*,*)i,k, A(k+(i-1)*n) ,  QR(k+(i-1)*n)
!          write(*,*)i,k, A(k+(i-1)*n) ,  Q(k+(i-1)*n) ,  R(1)  
       enddo
    enddo

    ss_int4=1
    sbuf=err_norm
    call MPI_Allreduce(sbuf,rbuf,ss_int4,MPI_REAL8,MPI_SUM,icomm,ierr)    
    err_norm=rbuf

    if(irank.eq.0)then
       err_norm=sqrt(err_norm)
!       write(*,*)"|| QtQ - I ||           = ",Orthnorm
       write(*,*)"||  A - QR || / sqrt(n*s) = ", err_norm/sqrt(n_)
    endif

  end subroutine check_A_QR

  subroutine QRDecomposition(n,s,icomm,V,Q,R,iQRflag,iQflag,nvBlock)
    use PARCEL_define
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    integer icomm
    real*8 V(n*s)
    real*8 Q(n*s)
    real*8 R(s*s)
    integer iQRflag
    integer,optional :: iQflag
    integer,optional :: nvBlock
    integer iQmatflag
    integer nvecBlock
    integer s_int4

    iQmatflag = 1
    if(  present(iQflag)  )then
       iQmatflag = iQflag
    endif

    nvecBlock = 100
    if(  present(nvBlock)  )then
       nvecBlock = nvBlock
    endif

    s_int4=s
    if(s_int4.eq.1)then
       ! 1本のベクトルに対する直交化なので正規化のみ行う
       call GS_QR(n,s,icomm,V,Q,R,KSP_QRDECOMPOSITION_MGS,nvecBlock)
       return
    endif

    !           select case (iQRflag)
    !           case(KSP_QRDECOMPOSITION_MGS)
    !              write(*,*)"@ MGS @"
    !           case(KSP_QRDECOMPOSITION_CGS)
    !              write(*,*)"@ CGS @"
    !           case(KSP_QRDECOMPOSITION_TSQR)
    !              write(*,*)"@ TSQR @"
    !           case(KSP_QRDECOMPOSITION_CHOLESKYQR)
    !              write(*,*)"@ CHOLESKYQR @"
    !           case(KSP_QRDECOMPOSITION_CHOLESKYQR2)
    !              write(*,*)"@ CHOLESKYQR2 @"
    !           end select


    ! ------
    select case (iQRflag)
    case(KSP_QRDECOMPOSITION_MGS)
       !       write(*,*)"MGS " 
       call GS_QR(n,s,icomm,V,Q,R,KSP_QRDECOMPOSITION_MGS,nvecBlock)
       return
    case(KSP_QRDECOMPOSITION_CGS)
       !       write(*,*)"CGS " 
       call GS_QR(n,s,icomm,V,Q,R,KSP_QRDECOMPOSITION_CGS,nvecBlock)
       return
    case(KSP_QRDECOMPOSITION_TSQR)
       !       write(*,*)"TSQR " 
       if(n>s)then
          call TSQR_Ronly(n,s,icomm,V,R,nvecBlock)
       else
          ! n<s は縦長行列ではないためTSQRで計算することは不可
          ! n<s の場合はCholeskyQR2で計算する
          call CholeskyQR2_Ronly(n,s,icomm,V,R,nvecBlock)
       endif
    case(KSP_QRDECOMPOSITION_CholeskyQR)
       !       write(*,*)"CholeskyQR " 
       call CholeskyQR_Ronly(n,s,icomm,V,R,nvecBlock)
    case(KSP_QRDECOMPOSITION_CholeskyQR2)
       !       write(*,*)"CholeskyQR2 " 
       call CholeskyQR2_Ronly(n,s,icomm,V,R,nvecBlock)
    end select
    if(iQmatflag.eq.KSP_QRDECOMPOSITION_TRUE_ORTHOGONAL_MATRIX)then
       call calcQ_VRinv(n,s,V,Q,R,nvecBlock)
    endif

  end subroutine QRDecomposition

  subroutine QRDecomposition_( &
       method, &
       n, &
       s, &
       icomm, &
       V, ldv,&
       Q, ldq,&
       R, ldr,&
       iQRflag, &
       iQflag, &
       nvBlock,&       
       localnRow, &
       batchCount, &
       V_batch_prt, &
       work_batch_prt, &
       work &
       )
!    integer,optional :: iQflag

    use krylov,only : KSP
    use c_wrapper,only : devicemem
    use vec_type,only : PARCEL_VEC
    use PARCEL_define
    use PARCEL_integer,only : intType
    implicit none
    type(KSP) :: method
    integer(kind=intType) n
    integer(kind=intType) s
    integer icomm
    type(PARCEL_VEC) V
    integer(kind=intType) ldv
    type(PARCEL_VEC) Q
    integer(kind=intType) ldq
    type(PARCEL_VEC) R
    integer(kind=intType) ldr

    integer iQRflag
    integer,optional :: iQflag
    integer,optional :: nvBlock

    integer(kind=intType),optional ::  localnRow
    integer(kind=intType),optional ::  batchCount
    type(devicemem),optional ::  V_batch_prt
    type(devicemem),optional ::  work_batch_prt
    type(PARCEL_VEC) ::  work

    integer iQmatflag
    integer nvecBlock
    integer s_int4


    iQmatflag = 1
    if(  present(iQflag)  )then
       iQmatflag = iQflag
    endif

    nvecBlock = 100
    if(  present(nvBlock)  )then
       nvecBlock = nvBlock
    endif

    s_int4=s
    if(s_int4.eq.1)then
       ! 1本のベクトルに対する直交化なので正規化のみ行う
!       call GS_QR(n,s,icomm,V,Q,R,KSP_QRDECOMPOSITION_MGS,nvecBlock)

       call GS_QR_( &
            method, &
            n, &
            s, &
            icomm, &
            V, ldv, &
            Q, ldq, &
            R, ldr, &
            KSP_QRDECOMPOSITION_MGS, &
            nvecBlock)
       return

    endif

    ! ------
    select case (iQRflag)
    case(KSP_QRDECOMPOSITION_MGS)
!       call GS_QR(n,s,icomm,V,Q,R,KSP_QRDECOMPOSITION_MGS,nvecBlock)
       call GS_QR_( &
            method, &
            n, &
            s, &
            icomm, &
            V, ldv, &
            Q, ldq, &
            R, ldr, &
            KSP_QRDECOMPOSITION_MGS, &
            nvecBlock)

       return
    case(KSP_QRDECOMPOSITION_CGS)
       !       write(*,*)"CGS " 
!       call GS_QR(n,s,icomm,V,Q,R,KSP_QRDECOMPOSITION_CGS,nvecBlock)
       call GS_QR_( &
            method, &
            n, &
            s, &
            icomm, &
            V, ldv, &
            Q, ldq, &
            R, ldr, &
            KSP_QRDECOMPOSITION_CGS, &
            nvecBlock)

       return
    case(KSP_QRDECOMPOSITION_TSQR)
!       write(*,*)"not support : TSQR " 
       if(n>s .and. n>64)then
#if 1
          call TSQR_Ronly_(&
               method, &
               n, &
               s, &
               icomm, &
               V, ldv,&
               R, ldr,&
               work, &
               nvecBlock)
#endif
#if 0
       call CholeskyQR_Ronly_( &
            method, &
            n, &
            s, &
            localnRow, &
            icomm, &
            batchCount, &
            V_batch_prt, &
            work_batch_prt, &
            V, ldv, &
            R, ldr, &
            work, &
            nvecBlock)
#endif
       else
          ! n<s は縦長行列ではないためTSQRで計算することは不可
          ! n<s の場合はCholeskyQR2で計算する
          ! call CholeskyQR2_Ronly(n,s,icomm,V,R,nvecBlock)
          call CholeskyQR2_Ronly_(&
               method, &
               n, &
               s, &
               localnRow, &
               icomm, &
               batchCount, &
               V_batch_prt, &
               work_batch_prt, &
               V, ldv, &
               R, ldr, &
               work, &
               nvecBlock)

       endif
    case(KSP_QRDECOMPOSITION_CholeskyQR)
!       write(*,*)"CholeskyQR " 
!       call CholeskyQR_Ronly(n,s,icomm,V,R,nvecBlock)
#if 1
       call CholeskyQR_Ronly_( &
            method, &
            n, &
            s, &
            localnRow, &
            icomm, &
            batchCount, &
            V_batch_prt, &
            work_batch_prt, &
            V, ldv, &
            R, ldr, &
            work, &
            nvecBlock)
#else

       call GS_QR_( &
            method, &
            n, &
            s, &
            icomm, &
            V, ldv, &
            Q, ldq, &
            R, ldr, &
            KSP_QRDECOMPOSITION_MGS, &
            nvecBlock)

       return
#endif

    case(KSP_QRDECOMPOSITION_CholeskyQR2)
!       call CholeskyQR2_Ronly(n,s,icomm,V,R,nvecBlock)
          call CholeskyQR2_Ronly_(&
               method, &
               n, &
               s, &
               localnRow, &
               icomm, &
               batchCount, &
               V_batch_prt, &
               work_batch_prt, &
               V, ldv, &
               R, ldr, &
               work, &
               nvecBlock)
    end select

    if(iQmatflag.eq.KSP_QRDECOMPOSITION_TRUE_ORTHOGONAL_MATRIX)then
       call calcQ_VRinv_(&
            method, &
            n, &
            s, &
            V, ldv,&
            Q, ldq,&
            R, ldr,&
            nvecBlock)
    endif

  end subroutine QRDecomposition_


  subroutine QR_TEST(n,s,icomm,Vorg,nvecBlock)
    use vec_kernel
    use PARCEL_define
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    integer icomm
    real*8 Vorg(n*s)
    integer nvecBlock

    real*8 V(n*s)
    real*8 Q(n*s)
    real*8 R(s*s)
    integer(kind=intType) i,j,k
    integer irank
    integer ierr
    integer iQRflag,iQmatflag
    integer(kind=intType) ns

    call MPI_Comm_rank(icomm,irank, ierr)
    iQmatflag = KSP_QRDECOMPOSITION_TRUE_ORTHOGONAL_MATRIX
    if(irank.eq.0)then
       write(*,*)" @@@@@@@@@@@@@@@ "
       write(*,*)"Matrix (n,s) :: ",n,s
    endif
    do i=1,5
       !    do i=1,1
       R=0.0d0
       ns=n*s
       call zero_initialize(ns,Q)
       call zero_initialize(ns,V)
       call xcpy(ns,Vorg,V)

       iQRflag   = i
       if(irank.eq.0)then
          select case (iQRflag)
          case(KSP_QRDECOMPOSITION_MGS)
             write(*,*)"@ MGS @"
          case(KSP_QRDECOMPOSITION_CGS)
             write(*,*)"@ CGS @"
          case(KSP_QRDECOMPOSITION_TSQR)
             write(*,*)"@ TSQR @"
          case(KSP_QRDECOMPOSITION_CHOLESKYQR)
             write(*,*)"@ CHOLESKYQR @"
          case(KSP_QRDECOMPOSITION_CHOLESKYQR2)
             write(*,*)"@ CHOLESKYQR2 @"
          end select
       endif

       call QRDecomposition(n,s,icomm,V,Q,R,iQRflag,iQmatflag,nvecBlock)
       call check_QtQ(n,s,icomm,Q)
       call check_A_QR(n,s,icomm,Vorg,Q,R)

       ! if(irank.eq.0)then
       !    do k=1,s
       !       do j=1,s
       !          write(*,*)"R = ",k,j,R(j+(k-1)*s)
       !       enddo
       !    enddo
       ! endif
    enddo
  end subroutine QR_TEST

end module QR_kernel

subroutine parcel_dqr(n,s,icomm,V,Q,R,iQRflag)
  use PARCEL_integer
  use QR_kernel,only : QRDecomposition
  implicit none

  integer(kind=intType) n
  integer(kind=intType) s
  integer icomm
  real*8 V(n*s)
  real*8 Q(n*s)
  real*8 R(s*s)
  integer iQRflag

  call QRDecomposition(n,s,icomm,V,Q,R,iQRflag)

end subroutine parcel_dqr

! QR分解の直交性を確認するためMPKで生成する行列に対してQR分解を行う
subroutine MPK_QR_TEST(method)
  use krylov
  use krylov_kernel
  use vec_kernel
  use QR_kernel

  use PARCEL_integer
  implicit none
  type(KSP) :: method
  integer(kind=intType)n
  integer icomm
  integer(kind=intType)s
  integer myrank
  integer nvecBlock
  integer(kind=intType)k
  integer(kind=intType)i,j
  integer(kind=intType)msize
  real*8,allocatable :: x(:)
  real*8,allocatable :: r(:)
  real*8,allocatable :: pv(:)
  real*8,allocatable :: V(:)
  real*8 rnorm
  real*8 tmp1,tmp2
  icomm=method%commVec%comm
  n=method%prm%n
  s=method%CAGMRES_sstep
  myrank=method%commVec%myrank
  nvecBlock=method%nvecBlock
  allocate(x(n))
  call zero_initialize(n,x)
  allocate(r(n))
  call zero_initialize(n,r)
  allocate(pv(n))
  call zero_initialize(n,pv)

  allocate(V(n*(s+1)))
  msize=n*(s+1)
  call zero_initialize(msize,V)
  !  call xcpy(n,method%x,x)
  do i=1,myrank
     call random_number(x)
  enddo
  ! r=b-Ax
  call calc_res(method,x,method%b,r)
  call calc_dot(n,icomm,r,r,rnorm)
  rnorm=sqrt(rnorm)
  tmp1=-1.0d0/rnorm
  tmp2=0.0d0
  ! v=r/||r||
  r = 1.0d0
  call axpy(n,tmp1,tmp2,r,V)

  !  call write_vec(s,V)
  do i=1,s
     ! Mpv=v_i
     call calc_precon(method,V(1+(i-1)*n),pV )
     ! V_i+1=ApV
     call calc_Matvec(method,pv,V(1+i*n) )
     !     call write_vec(s,V(1+i*n))
  enddo
!  do i=1,s+1
  call QR_TEST(n,i,icomm,V,nvecBlock)
!  enddo

  deallocate(x)
  deallocate(r)
  deallocate(pv)
  deallocate(V)

end subroutine MPK_QR_TEST

#if 0
program main
  use TSQR_kernel
  use QR_kernel
  use mpi
  use PARCEL_integer
  implicit none
  integer,parameter :: n=100*100*100
  integer,parameter :: s=4

  real*8 A(n*s)
  real*8 A0(n*s)
  real*8 R(s*s)
  real*8 Q(n*s)
  integer irank
  integer ierr
  integer i,j
  integer nblock,jend
  integer icomm
  call MPI_Init(ierr)  
  call MPI_Comm_rank(MPI_COMM_WORLD,irank, ierr)

  A=0.0
  A0=0.0
  Q=0.0
  R=0.0
  do i=1,irank+1
     call random_number(A)
  enddo
  A0 = A
  icomm = MPI_COMM_WORLD

 ! KSP_QRDECOMPOSITION_MGS=1        
 ! KSP_QRDECOMPOSITION_CGS=2        
 ! KSP_QRDECOMPOSITION_TSQR=3       
 ! KSP_QRDECOMPOSITION_CHOLESKYQR=4 
 ! KSP_QRDECOMPOSITION_CHOLESKYQR2=5

  do i = 1,5
     A = A0 
     Q = 1.0d0
     call QRDecomposition(n,s,icomm,A,Q,R,i)
     call check_QtQ(n,s,icomm,Q)
     call check_A_QR(n,s,icomm,A0,Q,R)
  enddo

  if(irank.eq.0)then
     ! do i=1,s
     !    do j=1,s
     !       write(*,*)j,i,R(j+(i-1)*s)
     !    enddo
     ! enddo
  endif

  call MPI_Finalize(ierr)

end program main
#endif

