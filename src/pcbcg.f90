! #define timer 1

#ifdef PARCEL_INT8
#define USE_BATCHDGEMM 0
#else
#define USE_BATCHDGEMM 1
#endif

module KSP_CBCG_CPU
contains
  ! SS=eta*AS - zeta * r
  subroutine CBCG_axpy2(n,SS,AS,r,eta,zeta)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 SS(n)
    real*8 AS(n)
    real*8 r(n)
    real*8 eta,zeta
    integer(kind=intType) i
    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
       SS(i)=eta*AS(i)-zeta*r(i)
    enddo
    !$OMP END PARALLEL DO 

  end subroutine CBCG_axpy2

  ! SS(i) = 2.0*eta*AS-2.0*zeta*SS(i-1)-SS(i-2)
  subroutine CBCG_axpy3(n,SS0,AS,SS1,SS2,eta2,zeta2)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 SS0(n)
    real*8 AS(n)
    real*8 SS1(n)
    real*8 SS2(n)
    real*8 eta2,zeta2
    integer(kind=intType) i

    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
       SS0(i)=eta2*AS(i)-zeta2*SS1(i)-SS2(i)
    enddo
    !$OMP END PARALLEL DO 

  end subroutine CBCG_axpy3

  subroutine make_chebyshev_basis_CPU(&
       pcbcg,&
       n,k,&
       Emax,Emin,&
       r,S,AS,SS)
    use krylov
    use krylov_kernel
    use PARCEL_integer
    implicit none
    type(KSP) :: pcbcg
    integer(kind=intType) n,k
    real*8 Emax,Emin
    real*8 r(n)
    real*8 S(n*k)
    real*8 SS(n*k)
    real*8 AS(n*k)

    real*8 eta
    real*8 zeta

    real*8 eta2
    real*8 zeta2

    integer(kind=intType) i
    integer(kind=intType) itmp

    eta   = 2.0d0/(Emax-Emin);
    zeta  = (Emax+Emin)/(Emax-Emin)
    eta2  = eta *2.0d0
    zeta2 = zeta*2.0d0

    call calc_precon(pcbcg,r,S)
    call calc_Matvec(pcbcg,S,AS)
    if(k.eq.1)then
       return 
    endif
    call CBCG_axpy2(n,SS,AS,r,eta,zeta)
    call calc_precon(pcbcg,SS,S(1+n))
    call calc_Matvec(pcbcg,S(1+n),AS(1+n))

    do i=2,k-1
       ! write(*,*)"make basis --",i
       itmp=2
       if(i.eq.itmp)then
          call CBCG_axpy3(n,&
               SS(1+(i-1)*n),&
               AS(1+(i-1)*n),&
               SS(1+(i-2)*n),&
               r,&
               eta2,zeta2)
       else
          call CBCG_axpy3(n,&
               SS(1+(i-1)*n),&
               AS(1+(i-1)*n),&
               SS(1+(i-2)*n),&
               SS(1+(i-3)*n),&
               eta2,zeta2)
       endif
       call calc_precon(pcbcg,SS(1+(i-1)*n), S(1+(i+0)*n))
       call calc_Matvec(pcbcg, S(1+(i+0)*n),AS(1+(i+0)*n))
    enddo

  end subroutine make_chebyshev_basis_CPU

  subroutine make_chebyshev_basis_CPU_(&
       pcbcg,&
       n,k,&
       Emax,Emin,&
       r,S,AS,SS,&
       vecLen)
    use krylov
    use krylov_kernel
    use PARCEL_integer
    implicit none
    type(KSP) :: pcbcg
    integer(kind=intType) n,k
    integer(kind=intType) vecLen
    real*8 Emax,Emin
    real*8 r(vecLen)
    real*8 S(vecLen*k)
    real*8 SS(vecLen*k)
    real*8 AS(vecLen*k)

    real*8 eta
    real*8 zeta

    real*8 eta2
    real*8 zeta2

    integer(kind=intType) i
    integer(kind=intType) itmp

    eta   = 2.0d0/(Emax-Emin);
    zeta  = (Emax+Emin)/(Emax-Emin)
    eta2  = eta *2.0d0
    zeta2 = zeta*2.0d0

    call calc_precon(pcbcg,r,S)
    call calc_Matvec(pcbcg,S,AS)
    if(k.eq.1)then
       return 
    endif
    call CBCG_axpy2(n,SS,AS,r,eta,zeta)
    call calc_precon(pcbcg,SS,S(1+vecLen))
    call calc_Matvec(pcbcg,S(1+vecLen),AS(1+vecLen))

    do i=2,k-1
       ! write(*,*)"make basis --",i
       itmp=2
       if(i.eq.itmp)then
          call CBCG_axpy3(n,&
               SS(1+(i-1)*vecLen),&
               AS(1+(i-1)*vecLen),&
               SS(1+(i-2)*vecLen),&
               r,&
               eta2,zeta2)
       else
          call CBCG_axpy3(n,&
               SS(1+(i-1)*vecLen),&
               AS(1+(i-1)*vecLen),&
               SS(1+(i-2)*vecLen),&
               SS(1+(i-3)*vecLen),&
               eta2,zeta2)
       endif
       call calc_precon(pcbcg,SS(1+(i-1)*vecLen), S(1+(i+0)*vecLen))
       call calc_Matvec(pcbcg, S(1+(i+0)*vecLen),AS(1+(i+0)*vecLen))
    enddo

  end subroutine make_chebyshev_basis_CPU_

  subroutine getEigenValue(method,Eigenflag,Emax,Emin)
    use krylov
    use PARCEL_define
    type(KSP) :: method
    integer Eigenflag
    real*8 Emax,Emin
    method%EigenMax_Real=0.0d0
    method%EigenMin_Real=0.0d0
    select case (Eigenflag)
    case(KSP_EIGENVALUE_TYPE_POWER)
       call power_method_symMatrix(method)
    case(KSP_EIGENVALUE_TYPE_CAARNOLDI)
       call KSP_pCAARNOLDI(method)
    end select
    Emax=method%EigenMax_Real
    Emin=method%EigenMin_Real

  end subroutine getEigenValue

  ! -------
  ! x  = x + Q * alpha
  ! Q  = S - Q * Beta
  ! AQ = AS - AQ * Beta
  ! Qt*AQ 
  ! Qtr
  subroutine CBCG_gemm_dot_local_CPU(&
       n,k,&
       x, ldx,&
       Q, ldq,&
       S, lds,&
       AQ, ldaq,&
       AS, ldas,&
       r,  ldr,&
       alpha, &
       Beta, &
       QtAQ, &
       Qtr, &
       rnorm, &
       nvecBlock)
    use vec_kernel
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) k
    integer(kind=intType) ldx
    integer(kind=intType) ldq
    integer(kind=intType) lds
    integer(kind=intType) ldaq
    integer(kind=intType) ldas
    integer(kind=intType) ldr
    real*8 x(ldx)
    real*8 Q(ldq*k)
    real*8 S(lds*k)
    real*8 AQ(ldaq*k)
    real*8 AS(ldas*k)
    real*8 r(ldr)
    real*8 alpha(k)
    real*8 Beta(k*k)
    real*8 QtAQ(k*k)
    real*8 Qtr(k)
    real*8 rnorm
    integer nvecBlock
    integer(kind=intType) i,j,l,jj

    QtAQ=0.0d0
    Qtr=0.0d0
    rnorm=0.0d0

    !$OMP PARALLEL DO PRIVATE(jj,l,j,i)  REDUCTION(+:QtAQ,Qtr,rnorm)
    do jj=1,n,nvecBlock 

       do l=1,k
          do j=jj,min(n,jj+nvecBlock-1)
             x(j) = x(j) + Q(j+(l-1)*ldq)*alpha(l)
          enddo
          do i=1,k
             do j=jj,min(n,jj+nvecBlock-1)
                S(j+(i-1)*lds)=S(j+(i-1)*lds)-Q(j+(l-1)*ldq)*Beta(l+(i-1)*k)
             enddo
          enddo
       enddo
       do l=1,k
          do j=jj,min(n,jj+nvecBlock-1)
             Q(j+(l-1)*ldq)=S(j+(l-1)*lds)
          enddo
       enddo

       do l=1,k
          do i=1,k
             do j=jj,min(n,jj+nvecBlock-1)
                AS(j+(i-1)*ldas)=AS(j+(i-1)*ldas)-AQ(j+(l-1)*ldaq)*Beta(l+(i-1)*k)
             enddo
          enddo
       enddo
       do l=1,k
          do j=jj,min(n,jj+nvecBlock-1)
             AQ(j+(l-1)*ldaq)=AS(j+(l-1)*ldas)
          enddo
       enddo

       do l=1,k
          do i=1,l
             do j=jj,min(n,jj+nvecBlock-1)
                QtAQ(i+(l-1)*k)=QtAQ(i+(l-1)*k)+Q(j+(l-1)*ldq)*AQ(j+(i-1)*ldaq)
             enddo
          enddo
          do j=jj,min(n,jj+nvecBlock-1)
             Qtr(l)=Qtr(l)+Q(j+(l-1)*ldq)*r(j)
          enddo
       enddo
       do j=jj,min(n,jj+nvecBlock-1)
          rnorm=rnorm+r(j)*r(j)
       enddo
    enddo
    !$OMP END PARALLEL DO 

  end subroutine CBCG_gemm_dot_local_CPU

  subroutine CBCG_ALLreduce_CPU(k,icomm,Qtr,QtAQ,rnorm)
    use mpi
    use PARCEL_integer
    implicit none
    integer(kind=intType) k
    integer icomm
    real*8 Qtr(k)
    real*8 QtAQ(k*k)
    real*8 rnorm

    real*8 sbuf((k+1)*k/2+k+1)
    real*8 rbuf((k+1)*k/2+k+1)
    integer(kind=intType) i,j,jj
    integer ierr
    integer ntmp

    ! All Reduce -- QtAQ,QtR
    jj=1
    do i=1,k
       do j=1,i
          sbuf(jj)=QtAQ(j+k*(i-1))
          jj=jj+1
       enddo
    enddo
    do i=1,k
       sbuf(jj)=Qtr(i)
       jj=jj+1
    enddo
    sbuf(jj)=rnorm
    jj=jj+1

    ntmp=(k+1)*k/2+k+1
    call MPI_Allreduce(sbuf,rbuf,ntmp,MPI_REAL8,MPI_SUM,icomm,ierr)    

    jj=1
    do i=1,k
       do j=1,i
          QtAQ(j+k*(i-1))=rbuf(jj)
          jj=jj+1
       enddo
    enddo
    do i=1,k
       Qtr(i)=rbuf(jj)
       jj=jj+1
    enddo
    rnorm=rbuf(jj)
    jj=jj+1
    do i=1,k
       do j=1,i
          QtAQ(i+k*(j-1))=QtAQ(j+k*(i-1))
       enddo
    enddo

  end subroutine CBCG_ALLreduce_CPU

end module KSP_CBCG_CPU

module KSP_CBCG
  use PARCEL_integer,only : intType
  use c_wrapper,only : devicemem

  interface
     subroutine CBCG_gemm_dot_local_FP64_GPU( &
          n, &
          k, &
          x , offset_x , ldx , &
          Q , offset_Q , ldq , &
          S , offset_S , lds , &
          AQ, offset_AQ, ldaq, &
          AS, offset_AS, ldas, &
          r , offset_r , ldr , &
          alpha, &
          beta, &
          QAQ, &
          Qr, &
          rnorm &
          ) bind(c,name="cuda_CBCG_gemm_dot_local")

       import :: intType
       import :: devicemem
       implicit none

       integer(kind=intType),value ::  n
       integer(kind=intType),value ::  k 
       type(devicemem) x
       integer(kind=intType),value ::  offset_x 
       integer(kind=intType),value ::  ldx 
       type(devicemem) Q
       integer(kind=intType),value ::  offset_Q 
       integer(kind=intType),value ::  ldq 
       type(devicemem) S
       integer(kind=intType),value ::  offset_S 
       integer(kind=intType),value ::  lds 
       type(devicemem) AQ
       integer(kind=intType),value ::  offset_AQ
       integer(kind=intType),value ::  ldaq 
       type(devicemem) AS
       integer(kind=intType),value ::  offset_AS
       integer(kind=intType),value ::  ldas 
       type(devicemem) r 
       integer(kind=intType),value ::  offset_r 
       integer(kind=intType),value ::  ldr 
       type(devicemem) alpha 
       type(devicemem) beta 
       type(devicemem) QAQ
       type(devicemem) Qr 
       type(devicemem) rnorm 
        
     end subroutine CBCG_gemm_dot_local_FP64_GPU

     subroutine CBCG_gemm_dot_local_ddm_FP64_GPU( &
          nx, ny, nz, &
          stmx, stmy, stmz, &
          mx, mxy, &
          k, &
          x , offset_x , ldx , &
          Q , offset_Q , ldq , &
          S , offset_S , lds , &
          AQ, offset_AQ, ldaq, &
          AS, offset_AS, ldas, &
          r , offset_r , ldr , &
          alpha, &
          beta, &
          QAQ, &
          Qr, &
          rnorm &
          ) bind(c,name="cuda_CBCG_gemm_dot_local_ddm")

       import :: intType
       import :: devicemem
       implicit none

       integer(kind=intType),value :: nx
       integer(kind=intType),value :: ny
       integer(kind=intType),value :: nz
       integer(kind=intType),value :: stmx
       integer(kind=intType),value :: stmy
       integer(kind=intType),value :: stmz
       integer(kind=intType),value :: mx
       integer(kind=intType),value :: mxy

       integer(kind=intType),value ::  k 
       type(devicemem) x
       integer(kind=intType),value ::  offset_x 
       integer(kind=intType),value ::  ldx 
       type(devicemem) Q
       integer(kind=intType),value ::  offset_Q 
       integer(kind=intType),value ::  ldq 
       type(devicemem) S
       integer(kind=intType),value ::  offset_S 
       integer(kind=intType),value ::  lds 
       type(devicemem) AQ
       integer(kind=intType),value ::  offset_AQ
       integer(kind=intType),value ::  ldaq 
       type(devicemem) AS
       integer(kind=intType),value ::  offset_AS
       integer(kind=intType),value ::  ldas 
       type(devicemem) r 
       integer(kind=intType),value ::  offset_r 
       integer(kind=intType),value ::  ldr 
       type(devicemem) alpha 
       type(devicemem) beta 
       type(devicemem) QAQ
       type(devicemem) Qr 
       type(devicemem) rnorm 
       
     end subroutine CBCG_gemm_dot_local_ddm_FP64_GPU

     subroutine CBCG_axpy2_FP64_GPU( &
		       numGPUblock, &
		       numGPUthread, &
		       n, &
		       alpha, &
		       beta, &
		       x, offset_x, &
		       y, offset_y, &
		       z, offset_z  &
       ) bind(c,name="cuda_CBCG_axpy2")

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
       type(devicemem) z
       integer(kind=intType),value :: offset_z
       
     end subroutine CBCG_axpy2_FP64_GPU

  end interface

contains
 
  subroutine CBCG_axpy2_( &
       method,&
       alpha,beta,&
       x,y,z,&
       n &
       )
    use krylov,only : KSP
    use vec_type,only : PARCEL_VEC
    use PARCEL_integer,only : intType
    use KSP_CBCG_CPU,only : CBCG_axpy2
    use KSP_STATUS,only : &
       STATUS_FLAG,PARCEL_CPU_FLAG_ON,PARCEL_GPU_FLAG_ON

    implicit none
    type(KSP) :: method
    real*8 alpha
    real*8 beta
    type(PARCEL_VEC) :: x
    type(PARCEL_VEC) :: y
    type(PARCEL_VEC) :: z
    integer(kind=intType) n
    integer(kind=intType) x_offset
    integer(kind=intType) y_offset
    integer(kind=intType) z_offset
    type(STATUS_FLAG) :: flag

    flag = method%flag

    x_offset = x%offset
    y_offset = y%offset
    z_offset = z%offset

    if(flag%machineflag.eq.PARCEL_CPU_FLAG_ON)then
       ! SS=eta*AS - zeta * r --> z = alpha * x + beta * y
       call CBCG_axpy2(n, &
            z%vec_FP64(1+z_offset), &
            x%vec_FP64(1+x_offset), &
            y%vec_FP64(1+y_offset), &
            alpha,beta)
    endif

    if(flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
       call CBCG_axpy2_FP64_GPU( &
            512, &
            512, &
            n, &
            alpha, &
            beta , &
            x%C_vec_FP64, x_offset, &
            y%C_vec_FP64, y_offset, &
            z%C_vec_FP64, z_offset  &
            )
    endif

  end subroutine CBCG_axpy2_

  subroutine CBCG_gemm_dot_local( &
       pcbcg,  &
       n,k,  &
       x, ldx,  &
       Q, ldq,  &
       S, lds,  &
       AQ, ldaq,  &
       AS, ldas,  &
       r, ldr,  &
       alpha,  &
       beta,  &
       QAQ,  &
       Qr,  &
       rnorm, &
       nvecBlock &
       )
    use krylov,only : KSP
    use vec_type,only : PARCEL_VEC
    use PARCEL_integer,only : intType
    use KSP_CBCG_CPU,only : CBCG_gemm_dot_local_CPU
    use KSP_STATUS,only : &
       STATUS_FLAG,PARCEL_CPU_FLAG_ON,PARCEL_GPU_FLAG_ON
    use PARCEL_define,only : MATRIX_FORMAT_DDM
    use cu_common_wrapper,only : cusparseDestroyDnVec,cusparseCreateDnVec_FP64

    use PARCEL_VEC_kernel,only : &
         set_vec_FP64,get_vec_FP64, &
         parcel_sync

    implicit none
    type(KSP) :: pcbcg
    integer(kind=intType) n,k

    type(PARCEL_VEC) :: x 
    integer(kind=intType) ldx 
    type(PARCEL_VEC) :: Q 
    integer(kind=intType) ldq 
    type(PARCEL_VEC) :: S 
    integer(kind=intType) lds 
    type(PARCEL_VEC) :: AQ
    integer(kind=intType) ldaq
    type(PARCEL_VEC) :: AS
    integer(kind=intType) ldas
    type(PARCEL_VEC) :: r 
    integer(kind=intType) ldr 
    type(PARCEL_VEC) :: alpha
    type(PARCEL_VEC) :: beta
    type(PARCEL_VEC) :: QAQ
    type(PARCEL_VEC) :: Qr
    type(PARCEL_VEC) :: rnorm
    integer nvecBlock

    integer(kind=intType),parameter :: offset_zero=0
    integer(kind=intType) nadd
    type(STATUS_FLAG) :: flag
    integer MatFormat 

    integer(kind=intType) nx_,ny_,nz_
    integer(kind=intType) stmx_,stmy_,stmz_
    integer(kind=intType) mx_,mxy_

#if 0

    real*8 x_(ldx)      
    real*8 Q_(ldq*k)    
    real*8 S_(lds*k)    
    real*8 AQ_(ldaq*k) 
    real*8 AS_(ldas*k) 
    real*8 r_(ldr) 
    real*8 alpha_(k)
    real*8 beta_(k*k)
    real*8 QAQ_(k*k)
    real*8 Qr_(k)
    real*8 rnorm_(1)

#endif

    MatFormat = pcbcg%prm%MatFormat 
    flag = pcbcg%flag

    select case(flag%machineflag)
       
    case(PARCEL_CPU_FLAG_ON)
       call CBCG_gemm_dot_local_CPU(&
            n,k,&
            x%vec_FP64(1+x%offset), ldx,&
            Q%vec_FP64(1+Q%offset), ldq,&
            S%vec_FP64(1+S%offset), lds,&
            AQ%vec_FP64(1+AQ%offset), ldaq,&
            AS%vec_FP64(1+AS%offset), ldas,&
            r%vec_FP64(1+r%offset),  ldr,&
            alpha%vec_FP64 , &
            Beta%vec_FP64 , &
            QAQ%vec_FP64 , &
            Qr%vec_FP64 , &
            rnorm%vec_FP64(1) , &
            nvecBlock)


    case(PARCEL_GPU_FLAG_ON)
       if(MatFormat.eq.MATRIX_FORMAT_DDM)then

          nx_   = pcbcg%ddmMat%prm%nx
          ny_   = pcbcg%ddmMat%prm%ny
          nz_   = pcbcg%ddmMat%prm%nz
          stmx_ = pcbcg%ddmMat%prm%stmx
          stmy_ = pcbcg%ddmMat%prm%stmy
          stmz_ = pcbcg%ddmMat%prm%stmz
          mx_   = pcbcg%ddmMat%prm%mx
          mxy_  = pcbcg%ddmMat%prm%mxy

       call CBCG_gemm_dot_local_ddm_FP64_GPU( &
            nx_, &
            ny_, &
            nz_, &
            stmx_, &
            stmy_, &
            stmz_, &
            mx_, &
            mxy_, &
            k, &
            x%C_vec_FP64 , x%offset , ldx , &
            Q%C_vec_FP64 , Q%offset , ldq , &
            S%C_vec_FP64 , S%offset , lds , &
            AQ%C_vec_FP64, AQ%offset, ldaq, &
            AS%C_vec_FP64, AS%offset, ldas, &
            r%C_vec_FP64 , r%offset , ldr , &
            alpha%C_vec_FP64, &
            beta%C_vec_FP64, &
            QAQ%C_vec_FP64, &
            Qr%C_vec_FP64, &
            rnorm%C_vec_FP64 &
            )

       else
#if 0
  call parcel_sync(pcbcg)
  call get_vec_FP64(pcbcg,x,x_,ldx)
  call get_vec_FP64(pcbcg,Q,Q_,ldq*k)
  call get_vec_FP64(pcbcg,S,S_,lds*k)
  call get_vec_FP64(pcbcg,AQ,AQ_,ldaq*k)
  call get_vec_FP64(pcbcg,AS,AS_,ldas*k)
  call get_vec_FP64(pcbcg,r,r_,ldr)

  call get_vec_FP64(pcbcg,alpha,alpha_,k)
  call get_vec_FP64(pcbcg,beta,beta_,k*k)
  call get_vec_FP64(pcbcg,QAQ,QAQ_,k*k)
  call get_vec_FP64(pcbcg,Qr,Qr_,k)
  call get_vec_FP64(pcbcg,rnorm,rnorm_,1)
  call parcel_sync(pcbcg)
       call CBCG_gemm_dot_local_CPU(&
            n,k,&
            x_(1+x%offset), ldx,&
            Q_(1+Q%offset), ldq,&
            S_(1+S%offset), lds,&
            AQ_(1+AQ%offset), ldaq,&
            AS_(1+AS%offset), ldas,&
            r_(1+r%offset),  ldr,&
            alpha_ , &
            Beta_ , &
            QAQ_ , &
            Qr_ , &
            rnorm_(1) , &
            nvecBlock)
  call parcel_sync(pcbcg)
  call set_vec_FP64(pcbcg,x,x_,ldx)
  ! call set_vec_FP64(pcbcg,Q,Q_,ldq*k)
  ! call set_vec_FP64(pcbcg,S,S_,lds*k)
  ! call set_vec_FP64(pcbcg,AQ,AQ_,ldaq*k)
  ! call set_vec_FP64(pcbcg,AS,AS_,ldas*k)
  call set_vec_FP64(pcbcg,Q,S_,ldq*k)
  call set_vec_FP64(pcbcg,S,Q_,lds*k)
  call set_vec_FP64(pcbcg,AQ,AS_,ldaq*k)
  call set_vec_FP64(pcbcg,AS,AQ_,ldas*k)
  call set_vec_FP64(pcbcg,r,r_,ldr)

  call set_vec_FP64(pcbcg,alpha,alpha_,k)
  call set_vec_FP64(pcbcg,beta,beta_,k*k)
  call set_vec_FP64(pcbcg,QAQ,QAQ_,k*k)
  call set_vec_FP64(pcbcg,Qr,Qr_,k)
  call set_vec_FP64(pcbcg,rnorm,rnorm_,1)
  call parcel_sync(pcbcg)

#if 0
       swap_ptr( &(Q_[0].a_d) ,&(S_[0].a_d)  );
       swap_ptr( &(AQ_[0].a_d),&(AS_[0].a_d) );
#endif

#else
          call CBCG_gemm_dot_local_FP64_GPU( &
               n, &
               k, &
               x%C_vec_FP64 , x%offset , ldx , &
               Q%C_vec_FP64 , Q%offset , ldq , &
               S%C_vec_FP64 , S%offset , lds , &
               AQ%C_vec_FP64, AQ%offset, ldaq, &
               AS%C_vec_FP64, AS%offset, ldas, &
               r%C_vec_FP64 , r%offset , ldr , &
               alpha%C_vec_FP64, &
               beta%C_vec_FP64, &
               QAQ%C_vec_FP64, &
               Qr%C_vec_FP64, &
               rnorm%C_vec_FP64 &
               )
#endif
    endif
       nadd = n + &
            pcbcg%crsMat%nColL + &
            pcbcg%crsMat%nColU

       call cusparseDestroyDnVec(AQ%C_vec_FP64)
       call cusparseDestroyDnVec(AS%C_vec_FP64)
       call cusparseDestroyDnVec(Q%C_vec_FP64)
       call cusparseDestroyDnVec(S%C_vec_FP64)

       call cusparseCreateDnVec_FP64(AQ%C_vec_FP64 , n , nadd , offset_zero)
       call cusparseCreateDnVec_FP64(AS%C_vec_FP64 , n , nadd , offset_zero)
       call cusparseCreateDnVec_FP64(Q%C_vec_FP64  , n , nadd , offset_zero)
       call cusparseCreateDnVec_FP64(S%C_vec_FP64  , n , nadd , offset_zero)

    end select

  end subroutine CBCG_gemm_dot_local

  subroutine make_chebyshev_basis(&
       pcbcg,&
       n,vecLen,k,&
       Emax,Emin,&
       r,S,AS,SS)
    use krylov,only : KSP
    use cu_comm,only : get_nrk_vec_FP64_GPU
    use vec_type,only : PARCEL_VEC
    use cu_Matrix,only : MatVec_cbcg1_dia_FP64_GPU,MatVec_cbcg2_dia_FP64_GPU
    use PARCEL_integer,only : intType
    use PARCEL_KRYLOV_MATRIX,only : calc_Matvec
    use PARCEL_KRYLOV_PRECON,only : calc_precon
    use PARCEL_VEC_kernel,only : &
         init_vec,free_vec,set_val_FP64,axpy_vec3,axpy_vec4,&
         set_vec_FP64,get_vec_FP64, &
         xcpy, &
         parcel_sync
    use KSP_STATUS,only : STATUS_FLAG,set_machineflag
    use PARCEL_define,only : &
         PARCEL_CPU_FLAG_ON,PARCEL_GPU_FLAG_ON,&
         MATRIX_FORMAT_DIA
    use c_wrapper,only :  CUDA_NUM_BLOCKS , CUDA_NUM_THREADS

    use KSP_CBCG_CPU
    use vec_kernel_CPU
    use comm_kernel
    use PARCEL_MATRIX_kernel_CPU

    implicit none
    type(KSP) :: pcbcg
    integer(kind=intType) n,k,vecLen

    real*8 Emax,Emin
    type(PARCEL_VEC) :: r
    type(PARCEL_VEC) :: S
    type(PARCEL_VEC) :: SS
    type(PARCEL_VEC) :: AS
    integer(kind=intType),parameter :: offset_zero=0

    real*8 eta
    real*8 zeta

    real*8 eta2
    real*8 zeta2

    integer(kind=intType) i
    integer(kind=intType) itmp

    type(PARCEL_VEC) :: tmp1
    type(PARCEL_VEC) :: tmp2
    type(PARCEL_VEC) :: tmp3
    type(PARCEL_VEC) :: tmp4

    type(STATUS_FLAG) :: flag
    integer MatFormat 
    integer(kind=intType) parcel_int
    integer(kind=intType) parcel_int1,parcel_int2,parcel_int3,parcel_int4,parcel_int5

#if 0
    ! これは問題なく動作するので原因はGPU版の関数?
    real*8 r_(vecLen)
    real*8 S_(vecLen*k)
    real*8 SS_(vecLen*k)
    real*8 AS_(vecLen*k)

  ! call set_vec_FP64(pcbcg,x,pcbcg%x,n)
  ! call get_vec_FP64(pcbcg,x,pcbcg%x,n)

  call parcel_sync(pcbcg)
  call get_vec_FP64(pcbcg,r,r_,vecLen)
  call get_vec_FP64(pcbcg,S,S_,vecLen*k)
  call get_vec_FP64(pcbcg,SS,SS_,vecLen*k)
  call get_vec_FP64(pcbcg,AS,AS_,vecLen*k)
  call parcel_sync(pcbcg)
  ! call set_machineflag(flag,PARCEL_CPU_FLAG_ON)
  ! call parcel_sync(pcbcg)
  call make_chebyshev_basis_CPU_(&
         pcbcg,&
         n,k,&
         Emax,Emin,&
         r_,S_,AS_,SS_,vecLen)
  ! call parcel_sync(pcbcg)
  ! call set_machineflag(flag,PARCEL_GPU_FLAG_ON)
  call parcel_sync(pcbcg)
  call set_vec_FP64(pcbcg,r,r_,vecLen)
  call set_vec_FP64(pcbcg,S,S_,vecLen*k)
  call set_vec_FP64(pcbcg,SS,SS_,vecLen*k)
  call set_vec_FP64(pcbcg,AS,AS_,vecLen*k)
  call parcel_sync(pcbcg)

  return 
#endif

    MatFormat = pcbcg%prm%MatFormat 
    flag = pcbcg%flag

    parcel_int = 1
    call init_vec(pcbcg,tmp1,parcel_int)
    call init_vec(pcbcg,tmp2,parcel_int)
    call init_vec(pcbcg,tmp3,parcel_int)
    call init_vec(pcbcg,tmp4,parcel_int)
    call parcel_sync(pcbcg)
    call set_val_FP64(pcbcg,tmp1,0.0d0,parcel_int)
    call set_val_FP64(pcbcg,tmp2,0.0d0,parcel_int)
    call set_val_FP64(pcbcg,tmp3,0.0d0,parcel_int)
    call set_val_FP64(pcbcg,tmp4,0.0d0,parcel_int)
    call parcel_sync(pcbcg)

    eta   = 2.0d0/(Emax-Emin);
    zeta  = (Emax+Emin)/(Emax-Emin)
    eta2  = eta *2.0d0
    zeta2 = zeta*2.0d0

    call parcel_sync(pcbcg)
    call calc_precon(pcbcg,r,S)
    call parcel_sync(pcbcg)
    if(k.eq.1)then
       call parcel_sync(pcbcg)
       call calc_Matvec(pcbcg,S,AS)
       call parcel_sync(pcbcg)
       return 
    else
       if( (flag%machineflag.eq.PARCEL_GPU_FLAG_ON) .and. &
            (MatFormat.eq.MATRIX_FORMAT_DIA) )then

          call parcel_sync(pcbcg)
          call get_nrk_vec_FP64_GPU(&
               pcbcg%commVec, &
               S%C_vec_FP64, S%offset, &
               pcbcg%C_sbuff, offset_zero, &
               pcbcg%C_buff , offset_zero)
          call parcel_sync(pcbcg)

          parcel_int = pcbcg%nBlock
          call MatVec_cbcg1_dia_FP64_GPU( &
               CUDA_NUM_BLOCKS,&
               CUDA_NUM_THREADS,&
               pcbcg%prm%n,&
               pcbcg%prm%nnd,&
               pcbcg%prm%istart,&
               parcel_int,&
               S%C_vec_FP64  ,S%offset, &
               AS%C_vec_FP64 ,AS%offset, &
               r%C_vec_FP64  ,r%offset, &
               SS%C_vec_FP64 ,SS%offset, &
               eta, &
               -zeta , &
               pcbcg%C_diaMat%diaA, &
               pcbcg%C_diaMat%offset , &
               pcbcg%C_diaMat%iLs, &
               pcbcg%C_buff &
               )
          call parcel_sync(pcbcg)

          do i=2,k-1
             call parcel_sync(pcbcg)
             SS%offset = vecLen*(i-2)
             S%offset  = vecLen*(i-1)
             call calc_precon(pcbcg,SS,S)
             call parcel_sync(pcbcg)
             SS%offset = 0
             S%offset  = 0

             call parcel_sync(pcbcg)
             S%offset  = vecLen*(i-1)
!             AS%offset = vecLen*(i-1)
             ! call calc_Matvec(pcbcg,S,AS)
             call get_nrk_vec_FP64_GPU(&
                  pcbcg%commVec, &
                  S%C_vec_FP64, S%offset, &
                  pcbcg%C_sbuff, offset_zero, &
                  pcbcg%C_buff , offset_zero)
             call parcel_sync(pcbcg)
             S%offset  = 0
!             AS%offset = 0

             itmp=2
             ! SS(i) = 2.0*eta*AS-2.0*zeta*SS(i-1)-SS(i-2)
             ! z(:)  = alpha*x(:) + beta*y(:) + gamma*z(:)
             if(i.eq.itmp)then
                ! call set_val_FP64(pcbcg, tmp1,  2.0d0*eta ,1)
                ! call set_val_FP64(pcbcg, tmp2, -2.0d0*zeta,1)
                ! call set_val_FP64(pcbcg, tmp3, -1.0d0     ,1)
                ! call set_val_FP64(pcbcg, tmp4,  0.0d0     ,1)
                ! call axpy_vec4(&
                !      pcbcg,&
                !      tmp1,tmp2,tmp3,tmp4,&
                !      AS,SS,r,SS,&
                !      n,&
                !      (i-1)*vecLen, (i-2)*vecLen, offset_zero, (i-1)*vecLen)

                call parcel_sync(pcbcg)
                parcel_int = pcbcg%nBlock
                parcel_int1 = (i-1)*vecLen
                parcel_int2 = (i-1)*vecLen
                parcel_int3 = (i-2)*vecLen
                parcel_int4 = offset_zero
                parcel_int5 = (i-1)*vecLen
                call MatVec_cbcg2_dia_FP64_GPU( &
                     CUDA_NUM_BLOCKS,&
                     CUDA_NUM_THREADS,&
                     pcbcg%prm%n,&
                     pcbcg%prm%nnd,&
                     pcbcg%prm%istart,&
                     parcel_int,&
                     S%C_vec_FP64  , parcel_int1, &
                     AS%C_vec_FP64 , parcel_int2, &
                     SS%C_vec_FP64 , parcel_int3, &
                     r%C_vec_FP64  , parcel_int4, &
                     SS%C_vec_FP64 , parcel_int5, &
                     ! S%C_vec_FP64  , (i-1)*vecLen, &
                     ! AS%C_vec_FP64 , (i-1)*vecLen, &
                     ! SS%C_vec_FP64 , (i-2)*vecLen, &
                     ! r%C_vec_FP64 , offset_zero , &
                     ! SS%C_vec_FP64 , (i-1)*vecLen, &
                     2.0d0*eta, &
                     -2.0*zeta, &
                     -1.0d0, &
                     pcbcg%C_diaMat%diaA, &
                     pcbcg%C_diaMat%offset , &
                     pcbcg%C_diaMat%iLs, &
                     pcbcg%C_buff &
                     )
                call parcel_sync(pcbcg)
             else
                ! call set_val_FP64(pcbcg, tmp1,  2.0d0*eta ,1)
                ! call set_val_FP64(pcbcg, tmp2, -2.0d0*zeta,1)
                ! call set_val_FP64(pcbcg, tmp3, -1.0d0     ,1)
                ! call set_val_FP64(pcbcg, tmp4,  0.0d0     ,1)
                ! call axpy_vec4(&
                !      pcbcg,&
                !      tmp1,tmp2,tmp3,tmp4,&
                !      AS,SS,SS,SS,&
                !      n,&
                !      (i-1)*vecLen, (i-2)*vecLen, (i-3)*vecLen, (i-1)*vecLen)
                call parcel_sync(pcbcg)
                parcel_int = pcbcg%nBlock
                parcel_int1 = (i-1)*vecLen 
                parcel_int2 = (i-1)*vecLen
                parcel_int3 = (i-2)*vecLen
                parcel_int4 = (i-3)*vecLen
                parcel_int5 = (i-1)*vecLen
                call MatVec_cbcg2_dia_FP64_GPU( &
                     CUDA_NUM_BLOCKS,&
                     CUDA_NUM_THREADS,&
                     pcbcg%prm%n,&
                     pcbcg%prm%nnd,&
                     pcbcg%prm%istart,&
                     parcel_int,&
                     S%C_vec_FP64  , parcel_int1, &
                     AS%C_vec_FP64 , parcel_int2, &
                     SS%C_vec_FP64 , parcel_int3, &
                     SS%C_vec_FP64 , parcel_int4, &
                     SS%C_vec_FP64 , parcel_int5, &
                     ! S%C_vec_FP64  , (i-1)*vecLen, &
                     ! AS%C_vec_FP64 , (i-1)*vecLen, &
                     ! SS%C_vec_FP64 , (i-2)*vecLen, &
                     ! SS%C_vec_FP64 , (i-3)*vecLen, &
                     ! SS%C_vec_FP64 , (i-1)*vecLen, &
                     2.0d0*eta, &
                     -2.0*zeta, &
                     -1.0d0, &
                     pcbcg%C_diaMat%diaA, &
                     pcbcg%C_diaMat%offset , &
                     pcbcg%C_diaMat%iLs, &
                     pcbcg%C_buff &
                     )
                call parcel_sync(pcbcg)
             endif

          enddo

       else

          call parcel_sync(pcbcg)
          call calc_Matvec(pcbcg,S,AS)
          call parcel_sync(pcbcg)
          ! SS(:)=eta*AS(:)-zeta*r(:)
          parcel_int1 = 1
          call parcel_sync(pcbcg)
          call set_val_FP64(pcbcg,tmp1 , eta  , parcel_int1 )
          call set_val_FP64(pcbcg,tmp2 , -zeta, parcel_int1 )
          call set_val_FP64(pcbcg,tmp3 , 0.0d0, parcel_int1 )
          call parcel_sync(pcbcg)
#if 0
! call axpy_vec3(pcbcg,tmp1,tmp2,tmp3,AS,r,SS,n)
  call parcel_sync(pcbcg)
  call get_vec_FP64(pcbcg,r,r_,vecLen)
  call get_vec_FP64(pcbcg,SS,SS_,vecLen)
  call get_vec_FP64(pcbcg,AS,AS_,vecLen)
  call parcel_sync(pcbcg)
       call axpy_vec3_FP64_CPU( &
            n, &
            eta , -zeta, 0.0d0, &
            AS_, &
            r_, &
            SS_ &
            )
  call parcel_sync(pcbcg)
  call set_vec_FP64(pcbcg,r,r_,vecLen)
  call set_vec_FP64(pcbcg,SS,SS_,vecLen)
  call set_vec_FP64(pcbcg,AS,AS_,vecLen)
  call parcel_sync(pcbcg)
#else
  call parcel_sync(pcbcg)
  call axpy_vec3(pcbcg,tmp1,tmp2,tmp3,AS,r,SS,n)
  call parcel_sync(pcbcg)
#endif

          do i=2,k-1
             call parcel_sync(pcbcg)
             SS%offset = vecLen*(i-2)
             S%offset  = vecLen*(i-1)
             call calc_precon(pcbcg,SS,S)
             call parcel_sync(pcbcg)
             SS%offset = 0
             S%offset  = 0

             S%offset  = vecLen*(i-1)
             AS%offset = vecLen*(i-1)
             call calc_Matvec(pcbcg,S,AS)
             call parcel_sync(pcbcg)
             S%offset  = 0
             AS%offset = 0

             itmp=2
             ! SS(i) = 2.0*eta*AS-2.0*zeta*SS(i-1)-SS(i-2)
             ! z(:)  = alpha*x(:) a+ beta*y(:) + gamma*z(:)
             if(i.eq.itmp)then
                parcel_int1 = 1
                call set_val_FP64(pcbcg, tmp1,  2.0d0*eta , parcel_int1 )
                call set_val_FP64(pcbcg, tmp2, -2.0d0*zeta, parcel_int1 )
                call set_val_FP64(pcbcg, tmp3, -1.0d0     , parcel_int1 )
                call set_val_FP64(pcbcg, tmp4,  0.0d0     , parcel_int1 )
                parcel_int1 = (i-1)*vecLen
                parcel_int2 = (i-2)*vecLen
                parcel_int3 = offset_zero
                parcel_int4 = (i-1)*vecLen
                call parcel_sync(pcbcg)
                call axpy_vec4(&
                     pcbcg,&
                     tmp1,tmp2,tmp3,tmp4,&
                     AS,SS,r,SS,&
                     n,&
                     parcel_int1,parcel_int2,parcel_int3,parcel_int4)
                call parcel_sync(pcbcg)
             else
                call parcel_sync(pcbcg)
                parcel_int1 = 1
                call set_val_FP64(pcbcg, tmp1,  2.0d0*eta , parcel_int1 )
                call set_val_FP64(pcbcg, tmp2, -2.0d0*zeta, parcel_int1 )
                call set_val_FP64(pcbcg, tmp3, -1.0d0     , parcel_int1 )
                call set_val_FP64(pcbcg, tmp4,  0.0d0     , parcel_int1 )

                call parcel_sync(pcbcg)
                parcel_int1 = (i-1)*vecLen
                parcel_int2 = (i-2)*vecLen
                parcel_int3 = (i-3)*vecLen
                parcel_int4 = (i-1)*vecLen
                call axpy_vec4(&
                     pcbcg,&
                     tmp1,tmp2,tmp3,tmp4,&
                     AS,SS,SS,SS,&
                     n,&
                     parcel_int1,parcel_int2,parcel_int3,parcel_int4)
!                     (i-1)*vecLen, (i-2)*vecLen, (i-3)*vecLen, (i-1)*vecLen)
                call parcel_sync(pcbcg)
             endif
             
          enddo

       endif
    endif

    call parcel_sync(pcbcg)
    SS%offset = (k-2)*vecLen
    S%offset  = (k-1)*vecLen
    call calc_precon(pcbcg,SS,S)
    call parcel_sync(pcbcg)
    SS%offset = 0
    S%offset  = 0

    S%offset  = (k-1)*vecLen
    AS%offset = (k-1)*vecLen
    call parcel_sync(pcbcg)
    call calc_Matvec(pcbcg, S,AS)
    S%offset  = 0
    AS%offset = 0
    call parcel_sync(pcbcg)

    call free_vec(pcbcg,tmp1)
    call free_vec(pcbcg,tmp2)
    call free_vec(pcbcg,tmp3)
    call free_vec(pcbcg,tmp4)
    call parcel_sync(pcbcg)

  end subroutine make_chebyshev_basis

end module KSP_CBCG

#if 0
#else

subroutine KSP_pCBCG( pcbcg &
      )
  use krylov,only : KSP
  use KSP_CBCG
  use PARCEL_define,only : PARCEL_CPU_FLAG_ON,PARCEL_GPU_FLAG_ON

  use KSP_STATUS,only : &
       STATUS_FLAG, &
       init_staus_flag,set_machineflag, &
       set_dataflag,set_calcflag
  use KSP_CBCG_CPU,only : &
       getEigenValue,CBCG_ALLreduce_CPU,&
       make_chebyshev_basis_CPU_, &
       CBCG_gemm_dot_local_CPU
  use PARCEL_KRYLOV_MATRIX,only : calc_Matvec,calc_res
  use PARCEL_KRYLOV_PRECON,only : calc_precon
  use PARCEL_VEC_kernel
  use PARCEL_integer,only : intType
  use wrap_DGEMM,only : wrap_BASIC_DGEMM
  use cu_common_wrapper,only : &
       cublas_init_DgemmBatched_gemm_dot_local_FP64, &
       cublas_free_DgemmBatched_gemm_dot_local_FP64, &
       cublas_DgemmBatched_gemm_dot_local_FP64
  use moore_penrose_kernel

#if 1
  use vec_kernel,only : calc_gemm_dot_CPU,calc_gemm_dot_local_CPU
#endif

  implicit none

  type(KSP) :: pcbcg
  integer(kind=intType) n
  integer icomm
  integer itrmax ! 最大反復回数
  real*8 rtolmax ! 相対残差 (収束条件)
  real*8 abstolmax  !  絶対残差 (収束条件)
  integer myrank
  type(STATUS_FLAG) :: flag

  real*8 Emax,Emin
  integer(kind=intType) k ! ステップ数
  integer k_int4
  integer(kind=intType) msize
  integer nvecBlock
  integer(kind=intType),parameter :: offset_zero=0

  type(PARCEL_VEC) :: x
  type(PARCEL_VEC) :: b
  type(PARCEL_VEC) :: r
  type(PARCEL_VEC) :: r_true
  type(PARCEL_VEC) :: AQ
  type(PARCEL_VEC) :: Q
  type(PARCEL_VEC) :: S
  type(PARCEL_VEC) :: AS
  type(PARCEL_VEC) :: SS
  type(PARCEL_VEC) :: alpha
  type(PARCEL_VEC) :: alpha_
  type(PARCEL_VEC) :: beta
  type(PARCEL_VEC) :: QtAS
  type(PARCEL_VEC) :: QtAQ
  type(PARCEL_VEC) :: QtAQ_inv
  type(PARCEL_VEC) :: Qtr
  type(PARCEL_VEC) :: tmp1
  type(PARCEL_VEC) :: rnorm
  type(PARCEL_VEC) :: bnorm

  ! --------------------------------
  ! Batch GEMM
  ! Q x AS = QtAS 
  ! Q x AQ = QtAQ
  integer(kind=intType) batchCount
  integer(kind=intType) localnRow
  type(devicemem) Q_batch_prt 
  type(devicemem) S_batch_prt 
  type(devicemem) AQ_batch_prt
  type(devicemem) AS_batch_prt
  type(devicemem) work_batch_prt

  ! real*8,allocatable :: Q_batch_prt_FP64
  ! real*8,allocatable :: AQ_batch_prt_FP64
  ! real*8,allocatable :: AS_batch_prt_FP64
  ! type(PARCEL_VEC) :: Q_batch_prt
  ! type(PARCEL_VEC) :: AQ_batch_prt
  ! type(PARCEL_VEC) :: AS_batch_prt

  ! work はSSでいいのでは?
  ! 逆にSSはwork領域で良いのでは?
  type(PARCEL_VEC) :: work

  real*8,allocatable :: QtAS_FP64(:)
  real*8,allocatable :: QtAQ_FP64(:)
  real*8,allocatable :: QtAQ_inv_FP64(:)
  real*8,allocatable :: Qtr_FP64(:)
  real*8,allocatable :: alpha_FP64(:)
  real*8,allocatable :: alpha__FP64(:)
  real*8,allocatable :: beta_FP64(:)

#if 0
  real*8,allocatable :: Q_(:)
  real*8,allocatable :: AS_(:)
  real*8,allocatable :: QtAS_(:)
#endif

  real*8 tmp1_FP64
  real*8 rnorm_FP64,rnorm_FP64_
  real*8 bnorm_FP64

  integer(kind=intType) vecLen
  integer(kind=intType) nColL
  integer(kind=intType) nColU

  integer(kind=intType) batch_len
  real*8 rtime(40)
  real*8 stime0,etime0
  real*8 stime,etime
  integer ierr
  integer itr
  integer(kind=intType) i
  integer(kind=intType) parcel_int

  rtime = 0.0d0

  itrmax    = pcbcg%itrmax
  rtolmax   = pcbcg%rtolmax
  abstolmax = pcbcg%abstolmax
  myrank    = pcbcg%commVec%myrank
  icomm     = pcbcg%commVec%comm
  n         = pcbcg%prm%n
  k         = pcbcg%CBCG_kstep 
  nvecBlock = pcbcg%nvecBlock

  call MPI_Barrier(icomm,ierr)

  stime0=MPI_Wtime()
  stime=MPI_Wtime()

#if 0
  call init_staus_flag(flag)
#if 0
  call set_machineflag(flag,PARCEL_CPU_FLAG_ON)
#else
  call set_machineflag(flag,PARCEL_GPU_FLAG_ON)
#endif
  ! ここではpcbcg内のflagから設定しているが将来はKSP構造体で設定したものを受け取る
  call set_dataflag(flag,KSP_DOUBLE_PRECISION)
  call set_calcflag(flag,KSP_DOUBLE_PRECISION)
  pcbcg%flag%machineflag = flag%machineflag
  pcbcg%flag%dataflag = flag%dataflag 
  pcbcg%flag%calcflag = flag%calcflag 
#endif

  call parcel_sync(pcbcg)

  nColL = 0
  nColU = 0
  if(pcbcg%prm%MatFormat.eq.MATRIX_FORMAT_CRS)then
     nColL = pcbcg%crsMat%nColL
     nColU = pcbcg%crsMat%nColU
  endif
  vecLen = n+nColL+nColU

  ! --------------------------------
  ! CUBLAS Batch paameter set
  ! CRS,DIA -- DDM で小行列のサイズを変更する
  if(pcbcg%flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
     !     localnRow = 32
     localnRow = 64 ! s =4
     !     localnRow = 128
     batch_len = n
     if( mod(n,localnRow).ne.0 )then
        batch_len = batch_len + ( localnRow - mod(n,localnRow) )
     endif
     batchCount = batch_len / localnRow
!     vecLen = max(vecLen,batch_len)
     if(vecLen>batch_len)then
        vecLen = vecLen
     else
        vecLen = batch_len
     endif
  endif

#if 0
  ! if(myrank.eq.0)then
  !    write(*,*)"vecLen = ",vecLen     
  ! endif
  call KSP_pCG( pcbcg )
  !  call KSP_pCG( pcbcg )
#endif
!  return 
  ! --------------------------------

  call init_vec(pcbcg,x       ,vecLen)
  call init_vec(pcbcg,b       ,vecLen)
  call init_vec(pcbcg,r       ,vecLen)
  call init_vec(pcbcg,r_true  ,vecLen)
  call init_vec(pcbcg,AQ      ,vecLen*k)
  call init_vec(pcbcg,Q       ,vecLen*k)
  call init_vec(pcbcg,S       ,vecLen*k)
  call init_vec(pcbcg,SS      ,vecLen*k)
  call init_vec(pcbcg,AS      ,vecLen*k)
  call init_vec(pcbcg,work    ,vecLen*k)
  call init_vec(pcbcg,alpha   ,k)
  call init_vec(pcbcg,alpha_  ,k)
  call init_vec(pcbcg,beta    ,k*k)
  call init_vec(pcbcg,QtAS    ,k*k)
  call init_vec(pcbcg,QtAQ    ,k*k)
  call init_vec(pcbcg,QtAQ_inv,k*k)
  call init_vec(pcbcg,Qtr     ,k)

  parcel_int = 1
  call init_vec(pcbcg,tmp1,parcel_int)
  call init_vec(pcbcg,rnorm,parcel_int)
  call init_vec(pcbcg,bnorm,parcel_int)
  call set_val_FP64(pcbcg,tmp1,0.0d0,parcel_int)
  call set_val_FP64(pcbcg,rnorm,0.0d0,parcel_int)
  call set_val_FP64(pcbcg,bnorm,0.0d0,parcel_int)

  parcel_int = 1
  call set_vec_val(pcbcg,x       ,tmp1,parcel_int,vecLen)
  call set_vec_val(pcbcg,b       ,tmp1,parcel_int,vecLen)
  call set_vec_val(pcbcg,r       ,tmp1,parcel_int,vecLen)
  call set_vec_val(pcbcg,AQ      ,tmp1,parcel_int,vecLen*k)
  call set_vec_val(pcbcg,Q       ,tmp1,parcel_int,vecLen*k)
  call set_vec_val(pcbcg,S       ,tmp1,parcel_int,vecLen*k)
  call set_vec_val(pcbcg,SS      ,tmp1,parcel_int,vecLen*k)
  call set_vec_val(pcbcg,AS      ,tmp1,parcel_int,vecLen*k)
  call set_vec_val(pcbcg,work    ,tmp1,parcel_int,vecLen*k)
  call set_vec_val(pcbcg,alpha   ,tmp1,parcel_int,k)
  call set_vec_val(pcbcg,alpha_  ,tmp1,parcel_int,k)
  call set_vec_val(pcbcg,beta    ,tmp1,parcel_int,k*k)
  call set_vec_val(pcbcg,QtAS    ,tmp1,parcel_int,k*k)
  call set_vec_val(pcbcg,QtAQ    ,tmp1,parcel_int,k*k)
  call set_vec_val(pcbcg,QtAQ_inv,tmp1,parcel_int,k*k)
  call set_vec_val(pcbcg,Qtr     ,tmp1,parcel_int,k)

  ! --------------------------------
  call parcel_sync(pcbcg)
#ifdef CUDA_SOLVER
  call cudaMemcpy_d2d_FP64(pcbcg%C_x,x%C_vec_FP64,n,offset_zero,offset_zero)
  call cudaMemcpy_d2d_FP64(pcbcg%C_b,b%C_vec_FP64,n,offset_zero,offset_zero)
#else
  call set_vec_FP64(pcbcg,x,pcbcg%x,n)
  call set_vec_FP64(pcbcg,b,pcbcg%b,n)
#endif
  call parcel_sync(pcbcg)

  allocate(QtAS_FP64(k*k))
  allocate(QtAQ_FP64(k*k))
  allocate(QtAQ_inv_FP64(k*k))
  allocate(Qtr_FP64(k))
  allocate(alpha_FP64(k))
  allocate(alpha__FP64(k))
  allocate(beta_FP64(k*k))

#if 0
  allocate(Q_(vecLen*k))
  allocate(AS_(vecLen*k))
  allocate(QtAS_(k*k))
#endif

  ! ------------------------------
  ! CUBLAS Batch init
#if USE_BATCHDGEMM
  if(pcbcg%flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
     parcel_int = batchCount
     call parcel_sync(pcbcg)
     call cublas_init_DgemmBatched_gemm_dot_local_FP64(&
          parcel_int, &
          Q%C_vec_FP64, Q%offset, localnRow,&
          Q_batch_prt&
          ) 
     call parcel_sync(pcbcg)
     call cublas_init_DgemmBatched_gemm_dot_local_FP64(&
          parcel_int, &
          S%C_vec_FP64, S%offset, localnRow,&
          S_batch_prt&
          ) 
     call parcel_sync(pcbcg)
     call cublas_init_DgemmBatched_gemm_dot_local_FP64(&
          parcel_int, &
          AQ%C_vec_FP64, AQ%offset, localnRow,&
          AQ_batch_prt&
          ) 
     call parcel_sync(pcbcg)
     call cublas_init_DgemmBatched_gemm_dot_local_FP64(&
          parcel_int, &
          AS%C_vec_FP64, AS%offset, localnRow,&
          AS_batch_prt&
          ) 
     call parcel_sync(pcbcg)
     call cublas_init_DgemmBatched_gemm_dot_local_FP64(&
          parcel_int, &
          work%C_vec_FP64, work%offset, k,  &
          work_batch_prt &
          ) 
     call parcel_sync(pcbcg)
  endif
#endif
  ! ------------------------------

  QtAS_FP64 = 0.0d0
  QtAQ_FP64 = 0.0d0
  QtAQ_inv_FP64 = 0.0d0
  Qtr_FP64 = 0.0d0

  etime=MPI_Wtime()
  rtime(1) =  rtime(1) + etime-stime

  ! 固有値計算
  call parcel_sync(pcbcg)
  stime=MPI_Wtime()
  call getEigenValue(pcbcg,pcbcg%CBCG_Eigenflag,Emax,Emin)
  call parcel_sync(pcbcg)

  etime=MPI_Wtime()
  rtime(2) =  rtime(2) + etime-stime
  if(Emax.eq.Emin)then
     if(Emax<0.0d0)then
        Emax=0.0d0
     else
        Emin=0.0d0
     endif
  endif
  if( (Emax.eq.0.0d0).and.(Emin.eq.0.0d0) )then
     Emax = 1.0d0
  endif
  if(myrank.eq.0)then
!     write(*,*)" Emax,Emin = ",Emax,Emin
  endif

  ! ------
  ! pcbcg start
  ! r = b-Ax
  ! S = T0(AM^-1)r,T1(AM^-1)r,T2(AM^-1)r,........Tk-1(AM^-1)r
  ! Q = S
  ! Qt*AQ,Qt*r
  ! do itr =1,....
  !   alpha = (Qt*AQ)^-1 * (Qt*r)
  !   r = r - AQ * alpha
  !   S = T0(AM^-1)r,T1(AM^-1)r,T2(AM^-1)r,........Tk-1(AM^-1)r
  !   ASはSの基底生成時に同時に作成
  !   QtA*S 
  !   beta=(Qt*AQ)^-1 * (Qt*AS)
  !   x  = x + Q * alpha
  !   Q  = S - Q * beta
  !   AQ = AS - AQ * beta
  !   Qt*AQ,Qt*r
  ! enddo
  ! ------

  call parcel_sync(pcbcg)
  stime=MPI_Wtime()
  ! r = b-Ax
  call calc_res(pcbcg,x,b,r)
  call parcel_sync(pcbcg)
  call calc_dot(pcbcg,icomm,b,b,bnorm,n)
  call parcel_sync(pcbcg)
  call calc_dot(pcbcg,icomm,r,r,rnorm,n)
  call parcel_sync(pcbcg)
  parcel_int = 1
  call sqrt_val(pcbcg,bnorm,parcel_int)
  call sqrt_val(pcbcg,rnorm,parcel_int)
  call parcel_sync(pcbcg)
  call get_val_FP64(pcbcg,bnorm,bnorm_FP64,parcel_int)
  call get_val_FP64(pcbcg,rnorm,rnorm_FP64,parcel_int)
  call parcel_sync(pcbcg)

  ! S = T0(AM^-1)r,T1(AM^-1)r,T2(AM^-1)r,........Tk-1(AM^-1)r
  call parcel_sync(pcbcg)
  call make_chebyshev_basis(&
       pcbcg,&
       n,vecLen,k,&
       Emax,Emin,&
       r,S,AS,SS)
  call parcel_sync(pcbcg)
  call xcpy(pcbcg, S, Q,vecLen*k)
  call xcpy(pcbcg,AS,AQ,vecLen*k)
  call parcel_sync(pcbcg)

  call calc_gemm_dot_local(&
       pcbcg,&
       n,k,k,&
       Q    , vecLen,&
       AQ   , vecLen,&
       QtAQ , k,&
       nvecBlock)
  call parcel_sync(pcbcg)

  parcel_int = 1
  call calc_gemm_dot_local(&
       pcbcg,&
       n,k, parcel_int,&
       Q    , vecLen,&
       r    , vecLen,&
       Qtr  , k,&
       nvecBlock)
  call calc_dot_local(pcbcg,r,r,rnorm,n)
  
  call parcel_sync(pcbcg)
  etime=MPI_Wtime()
  rtime(1) =  rtime(1) + etime-stime
  do itr = 1,itrmax/k

     call parcel_sync(pcbcg)
     parcel_int = size(Qtr_FP64)
     call get_vec_FP64(pcbcg , Qtr  , Qtr_FP64  , parcel_int)
     parcel_int = size(QtAQ_FP64)
     call get_vec_FP64(pcbcg , QtAQ , QtAQ_FP64 , parcel_int)
     parcel_int = 1
     call get_val_FP64(pcbcg , rnorm , rnorm_FP64 , parcel_int)
     call parcel_sync(pcbcg)
     call CBCG_ALLreduce_CPU(k , icomm , Qtr_FP64 , QtAQ_FP64 , rnorm_FP64)
     call parcel_sync(pcbcg)
     parcel_int = size(QtAQ_FP64)
     call set_vec_FP64(pcbcg , QtAQ , QtAQ_FP64 , parcel_int)
     parcel_int = size(Qtr_FP64)
     call set_vec_FP64(pcbcg , Qtr  , Qtr_FP64  , parcel_int)
     call parcel_sync(pcbcg)

#if 0
     ! -------------------------------
     if(myrank.eq.0)then
        write(*,*)
        write(*,*)" -------------- "
        do i = 1,size(Qtr_FP64)
           write(*,*)"Qr   : ",i,Qtr_FP64(i)
        enddo
        do i = 1,size(QtAQ_FP64)
           write(*,*)"QtAQ : ",i,QtAQ_FP64(i)
        enddo
     endif
     ! -------------------------------
#endif

     call parcel_sync(pcbcg)
     k_int4=k
     call moore_penrose(QtAQ_FP64,QtAQ_inv_FP64,k_int4,k_int4)
     parcel_int = 1
     call wrap_BASIC_DGEMM(k,parcel_int,k,QtAQ_inv_FP64,Qtr_FP64,alpha_FP64)
     parcel_int = size(alpha_FP64)
     call set_vec_FP64(pcbcg,alpha,alpha_FP64, parcel_int)
     call parcel_sync(pcbcg)

     stime=MPI_Wtime()
     do i=1,k
        alpha__FP64(i)=-alpha_FP64(i)
     enddo
     parcel_int = size(alpha__FP64)
     call set_vec_FP64(pcbcg,alpha_,alpha__FP64,parcel_int)
     call parcel_sync(pcbcg)

#if 0
     ! -------------------------------
     if(myrank.eq.0)then
        write(*,*)
        write(*,*)" -------------- "
        do i=1,k
           write(*,*)"alpha = ",alpha_FP64(i)
        enddo
     endif
     ! -------------------------------
#endif

#if 0
     call calc_res(pcbcg,x,b,r)
#endif

     call parcel_sync(pcbcg)
     call calc_gemv(&
          pcbcg, &
          AQ, vecLen, &
          r , vecLen, &
          alpha_,  &
          n, &
          k, &
          nvecBlock)
     call parcel_sync(pcbcg)


     call parcel_sync(pcbcg)
     etime=MPI_Wtime()
     rtime(3) =  rtime(3) + etime-stime
     rnorm_FP64 = sqrt(rnorm_FP64)
     pcbcg%reshistory(itr)=rnorm_FP64/bnorm_FP64
     if( (rnorm_FP64/bnorm_FP64<rtolmax)  ) then
        exit
     endif
     ! if(myrank.eq.0)then
     !    write(*,*)"@ res @", itr, pcbcg%reshistory(itr)
     ! endif

#if 0
     parcel_int = 1
     call parcel_sync(pcbcg)
     call calc_res(pcbcg,x,b,r_true)
     call parcel_sync(pcbcg)
     call calc_dot(pcbcg,icomm,r_true,r_true,rnorm,n)
     call parcel_sync(pcbcg)
     call sqrt_val(pcbcg,rnorm,parcel_int)
     call parcel_sync(pcbcg)
     call get_val_FP64(pcbcg,rnorm,rnorm_FP64,parcel_int)
     call parcel_sync(pcbcg)

     call parcel_sync(pcbcg)
     call calc_dot(pcbcg,icomm,r,r,rnorm,n)
     call parcel_sync(pcbcg)
     call sqrt_val(pcbcg,rnorm,parcel_int)
     call parcel_sync(pcbcg)
     call get_val_FP64(pcbcg,rnorm,rnorm_FP64_,parcel_int)
     call parcel_sync(pcbcg)

     if(myrank.eq.0)then
!        write(*,*)"@ true res @",itr,rnorm_FP64/bnorm_FP64,rnorm_FP64_/bnorm_FP64
     endif
#endif

     !   S = T0(AM^-1)r,T1(AM^-1)r,T2(AM^-1)r,........Tk-1(AM^-1)r
     call parcel_sync(pcbcg)
     stime=MPI_Wtime()
     call make_chebyshev_basis(&
          pcbcg,&
          n,vecLen,k,&
          Emax,Emin,&
          r,S,AS,SS)
     call parcel_sync(pcbcg)
     etime=MPI_Wtime()
     rtime(4) =  rtime(4) + etime-stime

     !   QtA*S
     stime=MPI_Wtime()
     select case(pcbcg%flag%machineflag)
     case(PARCEL_CPU_FLAG_ON)
     call calc_gemm_dot(&
          pcbcg,&
          icomm,&
          n,k,k,&
          Q    , vecLen,&
          AS   , vecLen,&
          QtAS , k,&
          nvecBlock)
    case(PARCEL_GPU_FLAG_ON)

#if USE_BATCHDGEMM
       parcel_int = batchCount
       call cublas_DgemmBatched_gemm_dot_FP64( &
             icomm, &
             pcbcg%cublas_solver, &
             parcel_int, &
             k, localnRow, &
             Q_batch_prt   , vecLen, &
             AS_batch_prt  , vecLen, &
             work_batch_prt, k*batchCount, &
             work%C_vec_FP64,  &
             QtAS%C_vec_FP64 ) 
#else

#if 0
  call parcel_sync(pcbcg)
  call get_vec_FP64(pcbcg,Q ,Q_,vecLen*k)
  call get_vec_FP64(pcbcg,AS,AS_,vecLen*k)
  call parcel_sync(pcbcg)
       call calc_gemm_dot_CPU(&
            n,k,k,&
            icomm,&
            Q_   , vecLen,&
            AS_  , vecLen,&
            QtAS_, k,&
       nvecBlock)
  call parcel_sync(pcbcg)
  call set_vec_FP64(pcbcg,QtAS,QtAS_,k*k)
  call parcel_sync(pcbcg)
#else
     call calc_gemm_dot(&
          pcbcg,&
          icomm,&
          n,k,k,&
          Q    , vecLen,&
          AS   , vecLen,&
          QtAS , k,&
          nvecBlock)
#endif

#endif

    end select

     call parcel_sync(pcbcg)

     !   beta=(Qt*AQ)^-1 * (Qt*AS)
     parcel_int = size(QtAS_FP64)
     call get_vec_FP64(pcbcg,QtAS,QtAS_FP64,parcel_int)
     call parcel_sync(pcbcg)

     call wrap_BASIC_DGEMM(k,k,k,QtAQ_inv_FP64,QtAS_FP64,beta_FP64)
     call parcel_sync(pcbcg)
     parcel_int = size(beta_FP64)
     call set_vec_FP64(pcbcg,beta,beta_FP64,parcel_int)
     call parcel_sync(pcbcg)

     etime=MPI_Wtime()
     rtime(5) =  rtime(5) + etime-stime

#if 0

     call parcel_sync(pcbcg)
     call get_vec_FP64(pcbcg,x,x_FP64,vecLen)
     call get_vec_FP64(pcbcg,r,r_FP64,vecLen)
     call get_vec_FP64(pcbcg,Q,Q_FP64,vecLen*k)
     call get_vec_FP64(pcbcg,S,S_FP64,vecLen*k)
     call get_vec_FP64(pcbcg,AQ,AQ_FP64,vecLen*k)
     call get_vec_FP64(pcbcg,AS,AS_FP64,vecLen*k)
     call parcel_sync(pcbcg)
     call CBCG_gemm_dot_local_CPU(&
       n,k,&
       x_FP64, vecLen,&
       Q_FP64, vecLen,&
       S_FP64, vecLen,&
       AQ_FP64, vecLen,&
       AS_FP64, vecLen,&
       r_FP64,  vecLen,&
       alpha_FP64, &
       Beta_FP64 , &
       QtAQ_FP64 , &
       Qtr_FP64  , &
       rnorm_FP64, &
       nvecBlock)
     call parcel_sync(pcbcg)
     call set_vec_FP64(pcbcg,x,x_FP64,vecLen)
     call set_vec_FP64(pcbcg,r,r_FP64,vecLen)
     call set_vec_FP64(pcbcg,Q,Q_FP64,vecLen*k)
     call set_vec_FP64(pcbcg,S,S_FP64,vecLen*k)
     call set_vec_FP64(pcbcg,AQ,AQ_FP64,vecLen*k)
     call set_vec_FP64(pcbcg,AS,AS_FP64,vecLen*k)


     parcel_int = size(QtAQ_FP64)
     call set_vec_FP64(pcbcg,QtAQ,QtAQ_FP64,parcel_int)
     parcel_int = size(Qtr_FP64)
     call set_vec_FP64(pcbcg,Qtr,Qtr_FP64,parcel_int)
     parcel_int = 1
     call set_val_FP64(pcbcg,rnorm,rnorm_FP64,1)

     call parcel_sync(pcbcg)

#else
   ! x  = x  + Q  * alpha
   ! Q  = S  - Q  * beta --> Sは反復内部で利用しない  --> Sをwork領域に利用可能
   ! AQ = AS - AQ * beta --> ASは反復内部で利用しない --> ASをwork領域に利用可能
   ! Qt*AQ  = <Q|A|Q>
   ! Qtr    = <Q|r>
   ! <r|r>  
   ! GPU版ではS,ASに出力してポインタをswapする.

     stime=MPI_Wtime()
     ! -------------------------------------------------
     parcel_int = size(beta_FP64)
     call set_vec_FP64(pcbcg,beta,beta_FP64,parcel_int)
     parcel_int = size(alpha_FP64)

#if 0
     alpha_FP64 = 0.0d0
#endif

     call set_vec_FP64(pcbcg,alpha,alpha_FP64,parcel_int)
     call parcel_sync(pcbcg)
     call CBCG_gemm_dot_local( &
          pcbcg,  &
          n,k,  &
          x, vecLen,  &
          Q, vecLen,  &
          S, vecLen,  &
          AQ, vecLen,  &
          AS, vecLen,  &
          r, vecLen,  &
          alpha,  &
          beta,  &
          QtAQ,  &
          Qtr,  &
          rnorm, &
          nvecBlock &
          )
     ! call calc_dot_local(pcbcg,r,r,rnorm,n)
     ! -------------------------------------------------
     call parcel_sync(pcbcg)

     ! Qt*AQ 
     select case(pcbcg%flag%machineflag)
     case(PARCEL_CPU_FLAG_ON)
        call calc_gemm_dot_local(&
             pcbcg, &
             n, k, k, &
             Q,  vecLen,&
             AQ, vecLen,&
             QtAQ, k,&
             nvecBlock) 
    case(PARCEL_GPU_FLAG_ON)

#if USE_BATCHDGEMM
        ! Q<-->S , AQ<-->AS のポインタを入れ替える
        call parcel_sync(pcbcg)
        call swap_ptr_ptr(Q_batch_prt,S_batch_prt)
        call swap_ptr_ptr(AQ_batch_prt,AS_batch_prt)
        call parcel_sync(pcbcg)

        parcel_int = batchCount
        call cublas_DgemmBatched_gemm_dot_local_FP64( &
                pcbcg%cublas_solver, &
                parcel_int, &
                k, localnRow, &
                Q_batch_prt   , vecLen, &
                AQ_batch_prt  , vecLen, &
                work_batch_prt, k*batchCount, &
                work%C_vec_FP64,  &
                QtAQ%C_vec_FP64 ) 
#else

#if 0
  call parcel_sync(pcbcg)
  call get_vec_FP64(pcbcg,Q ,Q_,vecLen*k)
  call get_vec_FP64(pcbcg,AQ,AS_,vecLen*k)
  call parcel_sync(pcbcg)
       call calc_gemm_dot_local_CPU(&
            n,k,k,&
            Q_   , vecLen,&
            AS_  , vecLen,&
            QtAS_, k,&
            nvecBlock)
  call parcel_sync(pcbcg)
  call set_vec_FP64(pcbcg,QtAQ,QtAS_,k*k)
  call parcel_sync(pcbcg)
#else
        call calc_gemm_dot_local(&
             pcbcg, &
             n, k, k, &
             Q,  vecLen,&
             AQ, vecLen,&
             QtAQ, k,&
             nvecBlock) 
#endif

#endif
    end select

#endif

     ! -------------------------------------------------

     etime=MPI_Wtime()
     rtime(6) =  rtime(6) + etime-stime

  enddo

  ! --------------------------------

  ! 真の残差を確認
#if 1 
  parcel_int = 1
  call calc_res(pcbcg,x,b,r)
  call calc_dot(pcbcg,icomm,r,r,rnorm,n)
  call sqrt_val(pcbcg,rnorm,parcel_int)
  call get_val_FP64(pcbcg,rnorm,rnorm_FP64,parcel_int)

  if(myrank.eq.0)then
!     write(*,*)"@ true res @",rnorm_FP64,rnorm_FP64/bnorm_FP64
  endif
#endif

  call parcel_sync(pcbcg)
#ifdef CUDA_SOLVER
  call cudaMemcpy_d2d_FP64(x%C_vec_FP64,pcbcg%C_x,n,offset_zero,offset_zero)
#else
  call get_vec_FP64(pcbcg,x,pcbcg%x,n)
#endif
  call parcel_sync(pcbcg)

  pcbcg%niter=min(itrmax,itr)
  if(itrmax<=itr)then
     pcbcg%iret=PARCEL_RETURNCODE_MAXITERATION
  else
     pcbcg%iret=PARCEL_RETURNCODE_NORMAL
  endif

  call MPI_Barrier(icomm,ierr)
  etime0=MPI_Wtime()

#ifdef timer
  if(myrank.eq.0)then
     write(*,*)"pcbcg               ",etime0-stime0
     write(*,*)"    init            ",rtime(1)
     write(*,*)"    EigenValue      ",rtime(2)
     write(*,*)"    gemv            ",rtime(3)
     write(*,*)"    chebyshev_basis ",rtime(4)
     write(*,*)"    gemm_dot        ",rtime(5)
     write(*,*)"    gemm_dot_local  ",rtime(6)

  endif
#endif
  
  ! --------------------------------
#if USE_BATCHDGEMM
  if(pcbcg%flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
     call cublas_free_DgemmBatched_gemm_dot_local_FP64(Q_batch_prt)
     call cublas_free_DgemmBatched_gemm_dot_local_FP64(S_batch_prt)
     call cublas_free_DgemmBatched_gemm_dot_local_FP64(AQ_batch_prt)
     call cublas_free_DgemmBatched_gemm_dot_local_FP64(AS_batch_prt)
     call cublas_free_DgemmBatched_gemm_dot_local_FP64(work_batch_prt)
  endif
#endif
  call free_vec(pcbcg,x)
  call free_vec(pcbcg,b)
  call free_vec(pcbcg,r)
  call free_vec(pcbcg,r_true)
  call free_vec(pcbcg,AQ)
  call free_vec(pcbcg,Q)
  call free_vec(pcbcg,S)
  call free_vec(pcbcg,SS)
  call free_vec(pcbcg,AS)
  call free_vec(pcbcg,work)
  call free_vec(pcbcg,alpha)
  call free_vec(pcbcg,alpha_)
  call free_vec(pcbcg,beta)
  call free_vec(pcbcg,QtAS)
  call free_vec(pcbcg,QtAQ)
  call free_vec(pcbcg,QtAQ_inv)
  call free_vec(pcbcg,Qtr)
  call free_vec(pcbcg,tmp1)
  call free_vec(pcbcg,rnorm)
  call free_vec(pcbcg,bnorm)

  deallocate(QtAS_FP64)
  deallocate(QtAQ_FP64)
  deallocate(QtAQ_inv_FP64)
  deallocate(Qtr_FP64)
  deallocate(alpha_FP64)
  deallocate(alpha__FP64)
  deallocate(beta_FP64)
#if 0
  deallocate(Q_)
  deallocate(AS_)
  deallocate(QtAS_)
#endif

end subroutine KSP_pCBCG

#endif
