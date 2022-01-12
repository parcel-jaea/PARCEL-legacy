
subroutine KSP_pGMRES(pgmres)
  use krylov
  use krylov_kernel
  use vec_kernel
  use QR_kernel
  use KSP_GMRES

  use PARCEL_integer
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

  integer(kind=intType) mm

  real*8 rtime(40)
  real*8 stime0,etime0
  real*8 stime,etime

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

  call zero_initialize(n,x)
  call zero_initialize(n,r)

  msize=n*(m+1)
  call zero_initialize(msize,v)
  call zero_initialize(n,pv)
  call zero_initialize(n,w)
  call zero_initialize(n,x0)
  call zero_initialize(n,v0)

  call xcpy(n,pgmres%x,x)

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
  call calc_dot(n,icomm,pgmres%b,pgmres%b,bnorm)
  bnorm=sqrt(bnorm)

  etime=MPI_Wtime()
  rtime(1) =  rtime(1) + etime-stime
  do itr = 1,itrmax,m
     stime=MPI_Wtime()
     ! -----
     e=0.0d0
     c=0.0d0
     s=0.0d0
     !     y=0.0d0
     hki=0.0d0
     Rtri=0.0d0

     ! r=b-Ax
     call calc_res(pgmres,x,pgmres%b,r)
     call calc_dot(n,icomm,r,r,rnorm)
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
        pgmres%reshistory(itr)=rnorm/bnorm 
     endif
#if 0

#if 1
     if( (myrank.eq.0) )then
#else
        if( (myrank.eq.0).and.(itr.eq.1) )then
#endif
           write(*,*)"@ res @",itr,rnorm,rnorm/bnorm
        endif

#else
        if( (myrank.eq.0) )then
           !        write(*,*)"@ true res @",itr,rnorm,rnorm/bnorm
        endif
#endif
        etime=MPI_Wtime()
        rtime(2) =  rtime(2) + etime-stime

        ! inner ループ
        do inner=1,m
           stime=MPI_Wtime()
           ! Mpv=v
           call calc_precon(pgmres,v(1+n*(inner-1)),pv )
           etime=MPI_Wtime()
           rtime(6) =  rtime(6) + etime-stime

           stime=MPI_Wtime()
           ! w=Apv
           call calc_Matvec( pgmres,pv,w )
           etime=MPI_Wtime()
           rtime(3) =  rtime(3) + etime-stime

           stime=MPI_Wtime()
           ! 直交化(古典,修正グラムシュミット)
#if 0
           do k=1,inner
              call calc_dot(n,icomm ,&
                   w,&
                   v(1+n*(k-1))   , &
                   hki(k+(inner-1)*(m+1)) )
              tmp1=-hki(k+(inner-1)*(m+1))
              tmp2=1.0d0
              call axpy(n,tmp1,tmp2,v(1+n*(k-1)),w)
           enddo
#else

           ! グラムシュミット直交化
           if(iGSflag.eq.KSP_GMRES_GS_TYPE_MGS)then
              !           write(*,*)"MGS !!! "
              call MGS_Orthonormalization( &
                   n, &
                   m, &
                   inner, &
                   icomm, &
                   w, &
                   v, &
                   Rtri )
           else
              !           write(*,*)"CGS !!! "
              call CGS_Orthonormalization( &
                   n, &
                   m, &
                   inner, &
                   icomm, &
                   w, &
                   v, &
                   Rtri,&
                   nvecBlock )
           endif
           do k=1,inner-1
              hki(k + (inner-1)*(m+1))=Rtri(k+m*(inner-1))
           enddo

           call calc_dot(n,icomm ,&
                w,&
                v(1+n*(inner-1))   , &
                hki(inner+(inner-1)*(m+1)) )
           tmp1=-hki(inner+(inner-1)*(m+1))
           tmp2=1.0d0
           call axpy(n,tmp1,tmp2,v(1+n*(inner-1)),w)
#endif

           ! h(k+1i)=||w||
           call calc_dot(n,icomm ,&
                w, &
                w, &
                hki(inner+1 + (inner-1)*(m+1)) )
           hki(inner+1 + (inner-1)*(m+1))  = &
                sqrt(hki(inner+1 + (inner-1)*(m+1)) )

           ! v_i+1 = w/||w|| 
           tmp1=1.0d0/hki(inner+1 + (inner-1)*(m+1)) 
           tmp2=0.0d0
           call axpy(n,tmp1,tmp2, w,v(1+n*(inner-1+1)))
           etime=MPI_Wtime()
           rtime(4) =  rtime(4) + etime-stime

           ! -----
           ! ギブンス回転
#if 0
           do k=1,inner-1
              tmp1=hki(k+(inner-1)*(m+1))*c(k)-hki(k+1+(inner-1)*(m+1))*s(k)
              tmp2=hki(k+(inner-1)*(m+1))*s(k)+hki(k+1+(inner-1)*(m+1))*c(k)
              hki(k+(inner-1)*(m+1))=tmp1
              hki(k+1+(inner-1)*(m+1))=tmp2
           enddo
           tmp1=hki(inner+1+(inner-1)*(m+1))/hki(inner+(inner-1)*(m+1))
           c(inner)=1.0d0/sqrt(1.0d0+tmp1*tmp1 )
           s(inner)=-c(inner)*tmp1
           tmp1=c(inner)*e(inner)-s(inner)*e(inner+1)
           tmp2=s(inner)*e(inner)+c(inner)*e(inner+1)
           e(inner)=tmp1
           e(inner+1)=tmp2
           hki(inner+  (inner-1)*(m+1))=c(inner)*hki(inner+(inner-1)*(m+1))-s(inner)*hki(inner+1+(inner-1)*(m+1))
           hki(inner+1+(inner-1)*(m+1))=0.0
#else
           call GivensRotations(m,hki,inner,c,s,e)
#endif
           !      //      printf("   %4d %4d %16.14e \n",itr,inner,fabs(e[inner+1]))
           rnorm=abs(e(inner+1))

           ! ------
           mm=inner
           pgmres%reshistory(itr+inner)=rnorm/bnorm 
           if(myrank.eq.0)then
              !           write(*,*)"@ res @",itr+inner,rnorm,rnorm/bnorm
           endif
           if( (rnorm/bnorm<rtolmax) ) then
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
        ! if(myrank.eq.0)then
        !    write(*,*)" update mm = ",mm
        ! endif

        ! do i=1,(m+1)
        !    write(*,*)"e    ",i,e(i)
        ! enddo
        ! do k=1,m
        !    do i=1,(m+1)
        !       write(*,*)"hki  ",i,k,hki(i+(k-1)*(m+1))
        !    enddo
        ! enddo

        stime=MPI_Wtime()
        !   解ベクトル更新
        !     call update_xVec(pgmres,n,m,mm,x,v,x0,v0,y,hki,e)
        call update_xVec(pgmres,n,m,mm,x,v,x0,v0,hki,e,nvecBlock)
        if( (rnorm/bnorm<rtolmax) ) then
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
#if 0
     call calc_res(pgmres,x,pgmres%b,r)
     call calc_dot(n,icomm,r,r,rnorm)
     rnorm=sqrt(rnorm)
     if(myrank.eq.0)then
        write(*,*)"    @ true res @",rnorm,rnorm/bnorm
     endif
#endif
     call xcpy(n,x,pgmres%x)

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

   end subroutine KSP_pGMRES


   subroutine KSP_ARNOLDI(parnoldi)
     use krylov
     use krylov_kernel
     use vec_kernel
     use QR_kernel
     use KSP_GMRES
     use KSP_CAGMRES,only : getHessenbergMat_Eigen,getRealEigen
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
              call MGS_Orthonormalization( &
                   n, &
                   m, &
                   inner, &
                   icomm, &
                   w, &
                   v, &
                   Rtri )
           else
              call CGS_Orthonormalization( &
                   n, &
                   m, &
                   inner, &
                   icomm, &
                   w, &
                   v, &
                   Rtri,&
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
