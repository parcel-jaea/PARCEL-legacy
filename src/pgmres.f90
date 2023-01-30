! #define timer 1

module KSP_GMRES
  use KSP_STATUS
  use PARCEL_integer
  use PARCEL_VEC_kernel
  use cu_Vec,only : cuda_sum_vec_FP64_GPU
  !  use KSP_GMRES
contains

  subroutine CGS_Orth( &
       pgmres,&
       icomm, &
       w,v,hki,&
       n,vecLen,vrow,&
       nvecBlock &
       )
    use KSP_GMRES_CPU,only : CGS_Orth_FP64_CPU
    use cu_Vec,only : cuda_CGS_Orth_FP64_GPU
    implicit none

    type(KSP) :: pgmres
    integer icomm
    type(PARCEL_VEC) :: w
    type(PARCEL_VEC) :: v
    type(PARCEL_VEC) :: hki

    integer(kind=intType) n
    integer(kind=intType) vecLen
    integer(kind=intType) vrow
    integer nvecBlock
    type(STATUS_FLAG) :: flag
    flag = pgmres%flag

    select case(flag%machineflag)

    case(PARCEL_CPU_FLAG_ON)
       call CGS_Orth_FP64_CPU(&
            icomm, &
            w%vec_FP64(1+w%offset),vecLen, &
            v%vec_FP64(1+v%offset),vecLen, &
            hki%vec_FP64(1+hki%offset),&
            n, &
            vrow,  &
            nvecBlock &
            )
    case(PARCEL_GPU_FLAG_ON)
       call cuda_CGS_Orth_FP64_GPU( &
            icomm,&
            w%C_vec_FP64  , w%offset, vecLen, &
            v%C_vec_FP64  , v%offset, vecLen, &
            hki%C_vec_FP64, hki%offset, &
            n, &
            vrow &
            ) 
    end select

  end subroutine CGS_Orth

end module KSP_GMRES

subroutine KSP_pGMRES(pgmres)
  use krylov
  use PARCEL_KRYLOV_MATRIX
  use PARCEL_KRYLOV_PRECON
  use PARCEL_VEC_kernel
  use PARCEL_integer
  use KSP_GMRES_CPU,only : GivensRotations,CGS_Orth_FP64_CPU
  use KSP_GMRES,only   : CGS_Orth
  use KSP_CAGMRES,only : update_xVec
!  use QR_kernel,only : MGS_Orthonormalization,CGS_Orthonormalization
  use cu_Vec,only : cuda_CGS_Orth_FP64_GPU
  implicit none

  type(KSP) :: pgmres
  integer(kind=intType) n
  integer icomm
  integer itrmax ! 最大反復回数
  real*8 rtolmax ! 相対残差 (収束条件)
  real*8 abstolmax  !  絶対残差 (収束条件)
  integer(kind=intType) m  ! GMRES(m) リスタート数
  integer iGSflag ! グラムシュミット直交化方法
  integer nvecBlock
  integer myrank
  type(STATUS_FLAG) :: flag

  integer(kind=intType) itr
  integer(kind=intType) inner

  integer(kind=intType) i,j,k
  integer ierr
  type(PARCEL_VEC) :: x
  type(PARCEL_VEC) :: b
  type(PARCEL_VEC) :: r
  type(PARCEL_VEC) :: v
  type(PARCEL_VEC) :: pv
  type(PARCEL_VEC) :: w
  type(PARCEL_VEC) :: x0
  type(PARCEL_VEC) :: v0
  type(PARCEL_VEC) :: y
  type(PARCEL_VEC) :: hki
  integer(kind=intType) msize

  real*8 bnorm_FP64
  real*8 rnorm_FP64
  real*8 tmp1_FP64
  real*8 tmp2_FP64

  type(PARCEL_VEC) :: bnorm
  type(PARCEL_VEC) :: rnorm
  type(PARCEL_VEC) :: tmp1
  type(PARCEL_VEC) :: tmp2

  real*8,allocatable :: e_FP64(:)
  real*8,allocatable :: c_FP64(:)
  real*8,allocatable :: s_FP64(:)
  real*8,allocatable :: hki_FP64(:)
  real*8,allocatable :: Rtri_FP64(:)

  integer(kind=intType) mm
  !  type(PARCEL_VEC) :: y
  real*8 rtime(40)
  real*8 stime0,etime0
  real*8 stime,etime

  integer(kind=intType) vecLen
  integer(kind=intType) nColL
  integer(kind=intType) nColU
  integer(kind=intType),parameter :: offset_zero=0
  integer(kind=intType) parcel_int

  rtime=0.0d0

  itrmax=pgmres%itrmax
  rtolmax=pgmres%rtolmax
  abstolmax=pgmres%abstolmax
  myrank=pgmres%commVec%myrank
  icomm=pgmres%commVec%comm
  n=pgmres%prm%n
  m=pgmres%GMRES_m
  iGSflag=pgmres%GMRES_GSflag
  nvecBlock=pgmres%nvecBlock

#if 0

#if 0
  call set_machineflag(flag,PARCEL_CPU_FLAG_ON)
#else
  call set_machineflag(flag,PARCEL_GPU_FLAG_ON)
#endif
  call set_dataflag(flag,KSP_DOUBLE_PRECISION)
  call set_calcflag(flag,KSP_DOUBLE_PRECISION)

  pgmres%flag%machineflag = flag%machineflag
  pgmres%flag%dataflag = flag%dataflag 
  pgmres%flag%calcflag = flag%calcflag 

#endif

  call MPI_Barrier(icomm,ierr)
  stime0=MPI_Wtime()
  stime=MPI_Wtime()

  call parcel_sync(pgmres)

  nColL = 0
  nColU = 0
  if(pgmres%prm%MatFormat.eq.MATRIX_FORMAT_CRS)then
     nColL = pgmres%crsMat%nColL
     nColU = pgmres%crsMat%nColU
  endif
  vecLen = n+nColL+nColU

  call init_vec(pgmres,x ,vecLen)
  call init_vec(pgmres,b ,vecLen)
  call init_vec(pgmres,r ,vecLen)
  call init_vec(pgmres,v ,vecLen*(m+1))

  call init_vec(pgmres,pv,vecLen)
  call init_vec(pgmres,w ,vecLen)
  call init_vec(pgmres,x0,vecLen)
  call init_vec(pgmres,v0,vecLen)
  call init_vec(pgmres,y,m+1)

  parcel_int = 1
  call init_vec(pgmres,hki,(m+1)*m)
  call init_vec(pgmres,bnorm,parcel_int)
  call init_vec(pgmres,rnorm,parcel_int)
  call init_vec(pgmres,tmp1,parcel_int)
  call init_vec(pgmres,tmp2,parcel_int)

  call set_val_FP64(pgmres,bnorm,0.0d0,parcel_int)
  call set_val_FP64(pgmres,rnorm,0.0d0,parcel_int)
  call set_val_FP64(pgmres,tmp1,0.0d0,parcel_int)
  call set_val_FP64(pgmres,tmp2,0.0d0,parcel_int)

  call set_vec_val(pgmres,x ,tmp1,parcel_int,vecLen)
  call set_vec_val(pgmres,b ,tmp1,parcel_int,vecLen)
  call set_vec_val(pgmres,r ,tmp1,parcel_int,vecLen)
  call set_vec_val(pgmres,v ,tmp1,parcel_int,vecLen*(m+1))
  call set_vec_val(pgmres,pv,tmp1,parcel_int,vecLen)
  call set_vec_val(pgmres,w ,tmp1,parcel_int,vecLen)
  call set_vec_val(pgmres,x0,tmp1,parcel_int,vecLen)
  call set_vec_val(pgmres,v0,tmp1,parcel_int,vecLen)
  call set_vec_val(pgmres,y ,tmp1,parcel_int,m+1)
  call set_vec_val(pgmres,hki ,tmp1,parcel_int,(m+1)*m)

  allocate(e_FP64(m+1))
  allocate(c_FP64(m+1))
  allocate(s_FP64(m+1))
  allocate(hki_FP64((m+1)*m ))
  allocate(Rtri_FP64((m+1)*(m+1)))

  e_FP64 = 0.0d0
  c_FP64 = 0.0d0
  s_FP64 = 0.0d0
  hki_FP64 = 0.0d0
  Rtri_FP64 = 0.0d0

  call parcel_sync(pgmres)
#ifdef CUDA_SOLVER
  call cudaMemcpy_d2d_FP64(pgmres%C_x,x%C_vec_FP64,n,offset_zero,offset_zero)
  call cudaMemcpy_d2d_FP64(pgmres%C_b,b%C_vec_FP64,n,offset_zero,offset_zero)
#else
  call set_vec_FP64(pgmres,x,pgmres%x,n)
  call set_vec_FP64(pgmres,b,pgmres%b,n)
#endif
  call parcel_sync(pgmres)

  ! ------
  ! pgmres start
  ! ------
  ! do itr=1,....
  !   r=b-Ax
  !   v_1=r/||r||
  !   do i=1,m
  !     solve Mpv=v_i
  !     w=Apv
  !     グラムシュミット直交化
  !     v_i+1 = w/||w||
  !     ギブンス回転
  !     if(||r||<rmax)exit
  !   enddo
  !   解ベクトル更新
  !   if(||r||<rmax)exit
  ! enddo

  inner=1
  parcel_int = 1
  call calc_dot(pgmres,icomm,b,b,bnorm,n)
  call sqrt_val(pgmres,bnorm,parcel_int)
  call calc_dot(pgmres,icomm,r,r,rnorm,n)
  call sqrt_val(pgmres,rnorm,parcel_int)
  call parcel_sync(pgmres)
  
  call get_val_FP64(pgmres,bnorm,bnorm_FP64,parcel_int)
  call get_val_FP64(pgmres,rnorm,rnorm_FP64,parcel_int)

  if(myrank.eq.0)then
!     write(*,*)" rnorm,bnorm = ",rnorm_FP64,bnorm_FP64
  endif

  etime=MPI_Wtime()
  rtime(1) =  rtime(1) + etime-stime
  do itr = 1,itrmax,m
     stime=MPI_Wtime()

     ! -----
     e_FP64 = 0.0d0
     c_FP64 = 0.0d0
     s_FP64 = 0.0d0
     hki_FP64 = 0.0d0
     Rtri_FP64 = 0.0d0

     parcel_int = 1
     call set_val_FP64(pgmres,tmp1,0.0d0,parcel_int)
     call set_vec_val(pgmres,hki ,tmp1,parcel_int,(m+1)*m)

     ! r = b - Ax
     call calc_res(pgmres,x,b,r)
     call calc_dot(pgmres,icomm,r,r,rnorm,n)
     parcel_int = 1
     call sqrt_val(pgmres,rnorm,parcel_int)
     call parcel_sync(pgmres)

     ! tmp1=-1.0d0/rnorm
     ! tmp2=0.0d0
     call get_val_FP64(pgmres,rnorm,rnorm_FP64,parcel_int)
     call set_val_FP64(pgmres,tmp1,-1.0d0,parcel_int)
     call div_val_FP64(pgmres,tmp1,rnorm_FP64,parcel_int)
     call set_val_FP64(pgmres,tmp2,0.0d0,parcel_int)

     ! v=r/||r||
     call axpy(pgmres,tmp1,tmp2,r,v,n)

     ! e(1)=-rnorm
     e_FP64(1) = -rnorm_FP64
     if( (rnorm_FP64<abstolmax) .or. (rnorm_FP64/bnorm_FP64<rtolmax) ) then
        exit
     endif
     if(itr.eq.1)then
        pgmres%reshistory(itr)=rnorm_FP64/bnorm_FP64 
     endif

     etime=MPI_Wtime()
     rtime(2) =  rtime(2) + etime-stime

     ! inner ループ
     do inner=1,m
        ! if(myrank.eq.0)then
        !    write(*,*)"@ inner @ = ",inner
        ! endif
        stime=MPI_Wtime()
        ! Mpv=v
        ! call calc_precon(pgmres,v(1+n*(inner-1)),pv)
        v%offset=vecLen*(inner-1)
        call calc_precon(pgmres,v,pv)
        v%offset=0
        etime=MPI_Wtime()
        rtime(6) =  rtime(6) + etime-stime

        stime=MPI_Wtime()
        ! w=Apv
        call calc_Matvec( pgmres,pv,w )
        etime=MPI_Wtime()
        rtime(3) =  rtime(3) + etime-stime

        stime=MPI_Wtime()
        ! 直交化(古典,修正グラムシュミット)

        if(iGSflag.eq.KSP_GMRES_GS_TYPE_MGS)then
           ! MGS
           do k=1,inner
              v%offset = vecLen*(k-1)
              call calc_dot(pgmres,icomm,w,v,tmp1,n)
              v%offset = 0

              parcel_int = 1
              call get_val_FP64(pgmres,tmp1,tmp1_FP64,parcel_int)
              hki_FP64(k+(inner-1)*(m+1)) = tmp1_FP64
              tmp1_FP64 = -tmp1_FP64
              tmp2_FP64 = 1.0d0
              
              parcel_int = 1
              call set_val_FP64(pgmres,tmp1,tmp1_FP64,parcel_int)
              call set_val_FP64(pgmres,tmp2,tmp2_FP64,parcel_int)

              v%offset = vecLen*(k-1)
              call axpy(pgmres,tmp1,tmp2,v,w,n)
              v%offset = 0
           enddo
        else
           ! CGS
           hki%offset = (inner-1)*(m+1)
           call CGS_Orth( &
                pgmres, &
                icomm, &
                w,v,hki,&
                n,vecLen,inner,&
                nvecBlock &
                )
           hki%offset = 0
#if 1
           do k=1,inner
              call get_val_FP64(pgmres,hki,tmp1_FP64,k+(inner-1)*(m+1))
              hki_FP64(k+(inner-1)*(m+1)) = tmp1_FP64
           enddo
#else
           hki%offset = (inner-1)*(m+1))
           call get_vec_FP64(pgmres,hki,hki_FP64(1+(inner-1)*(m+1)),inner)
           hki%offset = 0
#endif
        endif

        ! h(k+1i)=||w||
        call calc_dot(&
             pgmres,icomm ,&
             w, &
             w, &
             tmp1,n )
        parcel_int = 1
        call get_val_FP64(pgmres,tmp1,tmp1_FP64,parcel_int)
        tmp1_FP64 =sqrt(tmp1_FP64)
        !           call set_val_FP64(pgmres,hki ,tmp1_FP64,inner+1 + (inner-1)*(m+1))
        hki_FP64(inner+1 + (inner-1)*(m+1)) = tmp1_FP64

        ! v_i+1 = w/||w|| 
        tmp1_FP64 = 1.0d0/tmp1_FP64
        tmp2_FP64 = 0.0d0

        parcel_int = 1
        call set_val_FP64(pgmres,tmp1,tmp1_FP64,parcel_int)
        call set_val_FP64(pgmres,tmp2,tmp2_FP64,parcel_int)
        v%offset = vecLen*(inner-1+1)
        call axpy(pgmres,tmp1,tmp2, w,v,n)
        v%offset = 0

        etime=MPI_Wtime()
        rtime(4) =  rtime(4) + etime-stime

        ! -----
        ! ギブンス回転
        call GivensRotations(m,hki_FP64,inner,c_FP64,s_FP64,e_FP64)
        rnorm_FP64=abs(e_FP64(inner+1))
        ! ------
        mm=inner
        pgmres%reshistory(itr+inner) = rnorm_FP64/bnorm_FP64 
        if( (rnorm_FP64/bnorm_FP64<rtolmax) ) then
           exit
        endif
        if( (itr+inner).eq.itrmax)then
           exit
        endif

     enddo

     ! update
     if(m<mm)then
        mm=m
     endif
     ! do i=1,(m+1)
     !    write(*,*)"e    ",i,e_FP64(i)
     ! enddo
     ! if(myrank.eq.0)then
     !    do k=1,m
     !       do i=1,(m+1)
     !          write(*,*)"hki  ",i,k,hki_FP64(i+(k-1)*(m+1))
     !       enddo
     !    enddo
     ! endif

     stime=MPI_Wtime()
     !   解ベクトル更新
     call update_xVec(&
          pgmres, &
          n, &
          m, &
          mm, &
          vecLen, &
          x, &
          v, &
          x0, &
          v0, &
          hki_FP64, &
          e_FP64, &
          y, &
          nvecBlock)

     if( (rnorm_FP64/bnorm_FP64<rtolmax) ) then
        exit
     endif
     if( (itr+inner).eq.itrmax)then
        exit
     endif
     etime=MPI_Wtime()
     rtime(5) =  rtime(5) + etime-stime

  enddo

  call MPI_Barrier(icomm,ierr)
  etime0=MPI_Wtime()

  ! ------
  ! 真の残差を確認
#if 1
  call calc_res(pgmres,x,b,r)
  call calc_dot(pgmres,icomm,r,r,rnorm,n)
  parcel_int = 1
  call sqrt_val(pgmres,rnorm,parcel_int)
  call get_val_FP64(pgmres,rnorm,rnorm_FP64,parcel_int)

  if(myrank.eq.0)then
!     write(*,*)"@ true res @",rnorm_FP64,rnorm_FP64/bnorm_FP64
  endif
#endif
  call parcel_sync(pgmres)
#ifdef CUDA_SOLVER
  call cudaMemcpy_d2d_FP64(x%C_vec_FP64,pgmres%C_x,n,offset_zero,offset_zero)
#else
  call get_vec_FP64(pgmres,x,pgmres%x,n)
#endif
  call parcel_sync(pgmres)

  pgmres%niter=min(itrmax,itr+inner)
  if(itrmax<=itr)then
     pgmres%iret=PARCEL_RETURNCODE_MAXITERATION
  else
     pgmres%iret=PARCEL_RETURNCODE_NORMAL
  endif

#ifdef timer
  if(myrank.eq.0)then
     write(*,*)"gmres        ",etime0-stime0
     write(*,*)"    init     ",rtime(1)
     write(*,*)"    MatVec   ",rtime(3)
     write(*,*)"    Precon   ",rtime(6)
     write(*,*)"    QR       ",rtime(4)
     write(*,*)"    update   ",rtime(5)
  endif
#endif
  call parcel_sync(pgmres)
  call free_vec(pgmres,x)
  call free_vec(pgmres,b)
  call free_vec(pgmres,r)
  call free_vec(pgmres,v)
  call free_vec(pgmres,pv)
  call free_vec(pgmres,w)
  call free_vec(pgmres,x0)
  call free_vec(pgmres,v0)
  call free_vec(pgmres,y)
  call free_vec(pgmres,hki)
  call free_vec(pgmres,bnorm)
  call free_vec(pgmres,rnorm)
  call free_vec(pgmres,tmp1)
  call free_vec(pgmres,tmp2)

  deallocate(e_FP64)
  deallocate(c_FP64)
  deallocate(s_FP64)
  deallocate(hki_FP64)
  deallocate(Rtri_FP64)

  call parcel_sync(pgmres)

end subroutine KSP_pGMRES

   subroutine KSP_ARNOLDI(parnoldi)
     use krylov
     use krylov_kernel
     use vec_kernel
     use QR_kernel
     use KSP_GMRES_CPU,only : update_xvec_CPU,givensrotations
     use KSP_CAGMRES_CPU,only : getHessenbergMat_Eigen,getRealEigen
     use KSP_CAARNOLDI
     use PARCEL_integer
     implicit none

     type(KSP) :: parnoldi
     integer(kind=intType) n
     integer icomm
     integer itrmax ! 最大反復回数
     real*8 rtolmax ! 相対残差 (収束条件)
     real*8 abstolmax  !  絶対残差 (収束条件)
     integer(kind=intType) m  ! ARNOLDI(m) リスタート数
     integer iGSflag ! グラムシュミット直交化方法
     integer nvecBlock
     integer myrank

     integer(kind=intType) itr
     integer(kind=intType) inner
     ! ----
     integer(kind=intType) i,j,k
     integer ierr
     real*8,allocatable :: x(:)
     real*8,allocatable :: r(:)
     real*8,allocatable :: v(:)
     real*8,allocatable :: pv(:)
     real*8,allocatable :: w(:)
     real*8,allocatable :: x0(:)
     real*8,allocatable :: v0(:)
     integer(kind=intType) msize

     real*8 bnorm,rnorm
     real*8 tmp1,tmp2

     real*8,allocatable :: e(:)
     real*8,allocatable :: c(:)
     real*8,allocatable :: s(:)
     !  real*8,allocatable :: y(:)
     real*8,allocatable :: hki(:)
     real*8,allocatable :: Rtri(:)

     real*8,allocatable :: HessenbergMat(:)
     real*8,allocatable :: EigenReal(:)
     real*8,allocatable :: EigenImag(:)
     real*8,allocatable :: HessenbergMat_EigenVector(:)

     integer(kind=intType) mm
     real*8 rtime(40)
     real*8 stime0,etime0
     real*8 stime,etime

     rtime=0.0d0

     itrmax=parnoldi%itrmax
     rtolmax=parnoldi%rtolmax
     abstolmax=parnoldi%abstolmax
     myrank=parnoldi%commVec%myrank
     icomm=parnoldi%commVec%comm
     n=parnoldi%prm%n
     m=parnoldi%GMRES_m
     iGSflag=parnoldi%GMRES_GSflag
     nvecBlock=parnoldi%nvecBlock

     call MPI_Barrier(icomm,ierr)
     stime0=MPI_Wtime()
     stime=MPI_Wtime()

     allocate(x(n))
     allocate(r(n))
     allocate(v(n*(m+1)))
     allocate(pv(n))
     allocate(w(n))
     allocate(x0(n))
     allocate(v0(n))

     allocate(e(m+1))
     allocate(c(m+1))
     allocate(s(m+1))
     !  allocate(y(m+1))
     allocate(hki((m+1)*m ))
     allocate(Rtri((m+1)*(m+1)))
     allocate(HessenbergMat_EigenVector((m+1) * (m+1)))
     allocate(HessenbergMat((m+1)*(m+1)))
     allocate(EigenReal(m))
     allocate(EigenImag(m))

     call zero_initialize(n,x)
     call zero_initialize(n,r)

     msize=n*(m+1)
     call zero_initialize(msize,v)
     call zero_initialize(n,pv)
     call zero_initialize(n,w)
     call zero_initialize(n,x0)
     call zero_initialize(n,v0)

     call xcpy(n,parnoldi%x,x)

     ! ------
     ! parnoldi start
     ! ------
     ! do itr=1,....
     !   r=b-Ax
     !   v_1=r/||r||
     !   do i=1,m
     !     solve Mpv=v_i
     !     w=Apv
     !     グラムシュミット直交化
     !     v_i+1 = w/||w||
     !     ギブンス回転
     !     if(||r||<rmax)exit
     !   enddo
     !   解ベクトル更新
     !   if(||r||<rmax)exit
     ! enddo
     inner=1
     call calc_dot(n,icomm,parnoldi%b,parnoldi%b,bnorm)
     bnorm=sqrt(bnorm)

     etime=MPI_Wtime()
     rtime(1) =  rtime(1) + etime-stime
     !     do itr = 1,itrmax,m
     do itr = 1,1
        stime=MPI_Wtime()
        ! -----
        e=0.0d0
        c=0.0d0
        s=0.0d0
        !     y=0.0d0
        hki = 0.0d0
        Rtri = 0.0d0
        HessenbergMat = 0.0d0
        HessenbergMat_EigenVector = 0.0d0

        ! r=b-Ax
        call calc_res(parnoldi,x,parnoldi%b,r)
        call calc_dot(n,icomm,r,r,rnorm)
        if(rnorm.eq.0.0d0)then
           !$OMP PARALLEL DO PRIVATE(i)
           do i=1,n
              r(i)=1.0d0
           enddo
           call calc_dot(n,icomm,r,r,rnorm)
        endif
        rnorm=sqrt(rnorm)

        tmp1=-1.0d0/rnorm
        tmp2=0.0d0

        ! v=r/||r||
        call axpy(n,tmp1,tmp2,r,v)
        e(1)=-rnorm
        if( (rnorm<abstolmax) .or. (rnorm/bnorm<rtolmax) ) then
           exit
        endif
        if(itr.eq.1)then
           parnoldi%reshistory(itr)=rnorm/bnorm 
        endif

        if( (myrank.eq.0) )then
           !        write(*,*)"@ true res @",itr,rnorm,rnorm/bnorm
        endif

        etime=MPI_Wtime()
        rtime(2) =  rtime(2) + etime-stime

        ! inner ループ
        do inner=1,m
           stime=MPI_Wtime()
           ! Mpv=v
           call calc_precon(parnoldi,v(1+n*(inner-1)),pv )
           etime=MPI_Wtime()
           rtime(6) =  rtime(6) + etime-stime

           stime=MPI_Wtime()
           ! w=Apv
           call calc_Matvec( parnoldi,pv,w )
           etime=MPI_Wtime()
           rtime(3) =  rtime(3) + etime-stime

           stime=MPI_Wtime()
           ! 直交化(古典,修正グラムシュミット)
           ! グラムシュミット直交化
           if(iGSflag.eq.KSP_GMRES_GS_TYPE_MGS)then
              call MGS_Orthonormalization_CPU( &
                   n, &
                   m, &
                   inner, &
                   icomm, &
                   w, n,&
                   v, n, &
                   Rtri, m )
           else
              call CGS_Orthonormalization_CPU( &
                   n, &
                   m, &
                   inner, &
                   icomm, &
                   w, n,&
                   v, n,&
                   Rtri, m, &
                   nvecBlock )
           endif
           do k=1,inner-1
              hki(k + (inner-1)*(m+1))=Rtri(k+m*(inner-1))
              HessenbergMat( k + (inner-1)*(m+1) ) = hki(k + (inner-1)*(m+1))
           enddo

           call calc_dot(n,icomm ,&
                w,&
                v(1+n*(inner-1))   , &
                hki(inner+(inner-1)*(m+1) ) )
           HessenbergMat(inner+(inner-1)*(m+1)) = hki(inner+(inner-1)*(m+1))

           tmp1=-hki(inner+(inner-1)*(m+1))
           tmp2=1.0d0
           call axpy(n,tmp1,tmp2,v(1+n*(inner-1)),w)

           ! h(k+1i)=||w||
           call calc_dot(n,icomm ,&
                w, &
                w, &
                hki(inner+1 + (inner-1)*(m+1)) )
           hki(inner+1 + (inner-1)*(m+1))  = &
                sqrt(hki(inner+1 + (inner-1)*(m+1)) )
           HessenbergMat(inner+1 + (inner-1)*(m+1)) = hki(inner+1 + (inner-1)*(m+1))

           ! v_i+1 = w/||w|| 
           tmp1=1.0d0/hki(inner+1 + (inner-1)*(m+1)) 
           tmp2=0.0d0
           call axpy(n,tmp1,tmp2, w,v(1+n*(inner-1+1)))
           etime=MPI_Wtime()
           rtime(4) =  rtime(4) + etime-stime

           ! -----
           ! ギブンス回転
           call GivensRotations(m,hki,inner,c,s,e)

           !      //      printf("   %4d %4d %16.14e \n",itr,inner,fabs(e[inner+1]))
           rnorm=abs(e(inner+1))

           ! ------
           mm=inner
!           parnoldi%reshistory(itr+inner)=rnorm/bnorm 

        enddo

        EigenReal = 0.0d0
        EigenImag = 0.0d0
#if 0
        do i=1,(m+1)*m
           write(*,*)"Hes",i,HessenbergMat(i)
        enddo
#endif

        call getHessenbergMat_Eigen(&
             m,m,&
             HessenbergMat,&
             EigenReal,EigenImag, &
             HessenbergMat_EigenVector &
             )
#if 0
        call getRealEigen(&
             m,&
             EigenReal,EigenImag,&
             parnoldi%EigenMax_Real, &
             parnoldi%EigenMax_Imag, &
             parnoldi%EigenMin_Real, &
             parnoldi%EigenMin_Imag  &
             )
#endif
        !  if(iEigenVecflag.eq.1)then
        if(1.eq.1)then
           ! x = v * HessenbergMat_EigenVector : xが欲しい固有ベクトル
           call getEigenVector(n,m,v,HessenbergMat_EigenVector,nvecBlock)
           if( allocated(parnoldi%EigenVector) )then
              deallocate(parnoldi%EigenVector)
           endif
           allocate(parnoldi%EigenVector(n*(m)))
           call xcpy(n*(m),v,parnoldi%EigenVector)           
        endif
        if( allocated(parnoldi%EigenValueReal) ) deallocate(parnoldi%EigenValueReal)
        if( allocated(parnoldi%EigenValueImag) ) deallocate(parnoldi%EigenValueImag)
        if( allocated(parnoldi%EigenErrorNorm) ) deallocate(parnoldi%EigenErrorNorm)
        allocate(parnoldi%EigenValueReal(m))
        allocate(parnoldi%EigenValueImag(m))
        allocate(parnoldi%EigenErrorNorm(m))

        parnoldi%EigenValueReal = 0.0d0
        parnoldi%EigenValueImag = 0.0d0
        parnoldi%EigenErrorNorm = -1.0d0
        call checkEigenPair(parnoldi,n,m,v, &
             EigenReal,EigenImag, &
             parnoldi%MaxEigenID,parnoldi%MinEigenID,&
             parnoldi%EigenErrorNorm)
        do i=1,m
           parnoldi%EigenValueReal(i) = EigenReal(i)
           parnoldi%EigenValueImag(i) = EigenImag(i)
        enddo

        if(0<parnoldi%MaxEigenID)then
           parnoldi%EigenMax_Real = EigenReal(parnoldi%MaxEigenID) 
           parnoldi%EigenMax_Imag = EigenImag(parnoldi%MaxEigenID)
        else
           parnoldi%EigenMax_Real = 0.0d0
           parnoldi%EigenMax_Imag = 0.0d0
        endif
        if(0<parnoldi%MinEigenID)then
           parnoldi%EigenMin_Real = EigenReal(parnoldi%MinEigenID) 
           parnoldi%EigenMin_Imag = EigenImag(parnoldi%MinEigenID)
        else
           parnoldi%EigenMin_Real = 0.0d0
           parnoldi%EigenMin_Imag = 0.0d0
        endif

     enddo

     call MPI_Barrier(icomm,ierr)
     etime0=MPI_Wtime()

     ! ------
     ! 真の残差を確認
#if 0
     call calc_res(parnoldi,x,parnoldi%b,r)
     call calc_dot(n,icomm,r,r,rnorm)
     rnorm=sqrt(rnorm)
     if(myrank.eq.0)then
        write(*,*)"    @ true res @",rnorm,rnorm/bnorm
     endif
#endif
     call xcpy(n,x,parnoldi%x)

     parnoldi%niter=min(itrmax,itr+inner)
     if(itrmax<=itr)then
        parnoldi%iret=PARCEL_RETURNCODE_MAXITERATION
     else
        parnoldi%iret=PARCEL_RETURNCODE_NORMAL
     endif

#ifdef timer
     if(myrank.eq.0)then
        write(*,*)"arnoldi        ",etime0-stime0
        write(*,*)"    init     ",rtime(1)
        write(*,*)"    MatVec   ",rtime(3)
        write(*,*)"    Precon   ",rtime(6)
        write(*,*)"    QR       ",rtime(4)
     endif
#endif

     deallocate(x)
     deallocate(r)
     deallocate(v)
     deallocate(pv)
     deallocate(w)
     deallocate(x0)
     deallocate(v0)

     deallocate(e)
     deallocate(c)
     deallocate(s)
     !  deallocate(y)
     deallocate(hki)
     deallocate(Rtri)
     deallocate(HessenbergMat)
     deallocate(HessenbergMat_EigenVector)
     deallocate(EigenReal)
     deallocate(EigenImag)


   end subroutine KSP_ARNOLDI
