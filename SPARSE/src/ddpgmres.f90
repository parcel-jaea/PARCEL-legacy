module ddKSP_GMRES
  use dd_func
  implicit none
contains
  subroutine ddsum_vec(n,m,&
       v_hi,v_lo, &
       v0_hi,v0_lo, &
       y_hi,y_lo)
    use PARCEL_integer
    use vec_kernel
    implicit none
    integer(kind=intType) n
    integer(kind=intType) m

    real*8 v_hi(n*(m+1))
    real*8 v_lo(n*(m+1))
    real*8 v0_hi(n)
    real*8 v0_lo(n)
    real*8 y_hi(m)
    real*8 y_lo(m)

    real*8 tmp1_hi,tmp1_lo
    real*8 tmp2_hi,tmp2_lo
    integer(kind=intType) i,j,jj,k

    ! v0 = sum_k y_k*v(k)
#if 0
    do i=1,m
       if(i.eq.1)then
          tmp1=y(i)
          tmp2=0.0d0
          call axpy(n,tmp1,tmp2,v(1+(i-1)*n),v0)
       else
          tmp1=y(i)
          tmp2=1.0d0
          call axpy(n,tmp1,tmp2,v(1+(i-1)*n),v0)
       endif
    enddo
#else
    !$OMP PARALLEL DO PRIVATE(jj,j,i,tmp1_hi,tmp1_lo) 
    do jj=1,n,nvecBlock 
       do j=jj,min(n,jj+nvecBlock-1)
          v0_hi(j)=0.0d0
          v0_lo(j)=0.0d0
       enddo
       do i=1,m
          do j=jj,min(n,jj+nvecBlock-1)
!             v0(j)=v0(j)+y(i)*v(j+n*(i-1))
             call dd_mul(&
                  y_hi(i),y_lo(i),&
                  v_hi(j+n*(i-1)),v_lo(j+n*(i-1)),&
                  tmp1_hi,tmp1_lo)
             call dd_add(&
                  v0_hi(j),v0_lo(j),&
                  tmp1_hi,tmp1_lo, &
                  v0_hi(j),v0_lo(j))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO 
#endif
  end subroutine ddsum_vec

  subroutine ddupdate_xVec(&
       pgmres,n,m0,m,&
       x_hi,x_lo, &
       v_hi,v_lo, &
       x0_hi,x0_lo, &
       v0_hi,v0_lo, &
       y_hi,y_lo, &
       hki_hi,hki_lo, &
       e_hi,e_lo)
    use krylov
    use ddkrylov_kernel
    use ddvec_kernel
    use PARCEL_integer
    implicit none

    type(KSP) :: pgmres
    integer(kind=intType) n
    integer(kind=intType) m
    integer(kind=intType) m0
    real*8 x_hi(n)
    real*8 x_lo(n)
    real*8 v_hi(n*(m0+1))
    real*8 v_lo(n*(m0+1))
    real*8 x0_hi(n)
    real*8 x0_lo(n)
    real*8 v0_hi(n)
    real*8 v0_lo(n)
    real*8 y_hi(m0+1)
    real*8 y_lo(m0+1)
    real*8 hki_hi( (m0+1)*m0 )
    real*8 hki_lo( (m0+1)*m0 )
    real*8 e_hi(m0+1)
    real*8 e_lo(m0+1)
    real*8 Hinv_hi( m*m )
    real*8 Hinv_lo( m*m )

    real*8 tmp1_hi,tmp1_lo
    real*8 tmp2_hi,tmp2_lo

    integer(kind=intType) i,j
    integer ierr

    ! -----
    ! y = H^-1 e 
    do i=m,1,-1
       y_hi(i)=e_hi(i)
       y_lo(i)=e_lo(i)
       do j=m,i+1,-1
!          y(i)=y(i)-hki(i+(j-1)*(m0+1))*y(j)
          call dd_mul(&
               hki_hi(i+(j-1)*(m0+1)),hki_lo(i+(j-1)*(m0+1)),&
               y_hi(j),y_lo(j),&
               tmp1_hi,tmp1_lo)
          call dd_sub(&
               y_hi(i),y_lo(i),&
               tmp1_hi,tmp1_lo,&
               y_hi(i),y_lo(i))
       enddo
!       y(i)=y(i)/hki(i+(i-1)*(m0+1))
       if( ( hki_hi(i+(i-1)*(m0+1)).eq.0.0d0 ).and.&
           ( hki_lo(i+(i-1)*(m0+1)).eq.0.0d0 ) )then
          call dd_div(&
               y_hi(i),y_lo(i),&
               hki_hi(i+(i-1)*(m0+1)),hki_lo(i+(i-1)*(m0+1)),&
               1.0d-100,0.0d0)
       else
          call dd_div(&
               y_hi(i),y_lo(i),&
               hki_hi(i+(i-1)*(m0+1)),hki_lo(i+(i-1)*(m0+1)),&
               y_hi(i),y_lo(i))
       endif
    enddo

    ! v0 = sum_k y_k*v(k)
    call ddsum_vec(n,m,&
         v_hi ,v_lo,&
         v0_hi,v0_lo,&
         y_hi ,y_lo)

    ! solve
    call ddcalc_precon(pgmres,&
         v0_hi,v0_lo,&
         x0_hi,x0_lo )

!    call write_vec(n,x0)
    tmp1_hi=1.0d0
    tmp1_lo=0.0d0
    tmp2_hi=1.0d0
    tmp2_lo=0.0d0
    ! x = x+x0
    call ddaxpy(n,&
         tmp1_hi,tmp1_lo,&
         tmp2_hi,tmp2_lo,&
         x0_hi,x0_lo,&
         x_hi,x_lo)

  end subroutine ddupdate_xVec

  subroutine ddGivensRotations(m,&
       hki_hi,hki_lo,&
       i,&
       c_hi,c_lo,&
       s_hi,s_lo,&
       e_hi,e_lo)
    use PARCEL_integer
    implicit none
    integer(kind=intType) m
    integer(kind=intType) i
    real*8 hki_hi((m+1)*m)
    real*8 hki_lo((m+1)*m)
    real*8 c_hi(m+1)
    real*8 c_lo(m+1)
    real*8 s_hi(m+1)
    real*8 s_lo(m+1)
    real*8 e_hi(m+1)
    real*8 e_lo(m+1)

    real*8 tmp1_hi,tmp1_lo
    real*8 tmp2_hi,tmp2_lo

    real*8 tmp3_hi,tmp3_lo
    real*8 tmp4_hi,tmp4_lo
    integer(kind=intType) k

    do k=1,i-1
!       tmp1=hki(k+(i-1)*(m+1))*c(k)-hki(k+1+(i-1)*(m+1))*s(k)
       call dd_mul(&
            hki_hi(k+(i-1)*(m+1)),hki_lo(k+(i-1)*(m+1)),&
            c_hi(k),c_lo(k),&
            tmp3_hi,tmp3_lo)
       call dd_mul(&
            hki_hi(k+1+(i-1)*(m+1)),hki_lo(k+1+(i-1)*(m+1)),&
            s_hi(k),s_lo(k),&
            tmp4_hi,tmp4_lo)
       call dd_sub(&
            tmp3_hi,tmp3_lo,&
            tmp4_hi,tmp4_lo,&
            tmp1_hi,tmp1_lo)

!       tmp2=hki(k+(i-1)*(m+1))*s(k)+hki(k+1+(i-1)*(m+1))*c(k)
       call dd_mul(&
            hki_hi(k+(i-1)*(m+1)),hki_lo(k+(i-1)*(m+1)),&
            s_hi(k),s_lo(k),&
            tmp3_hi,tmp3_lo)
       call dd_mul(&
            hki_hi(k+1+(i-1)*(m+1)),hki_lo(k+1+(i-1)*(m+1)),&
            c_hi(k),c_lo(k),&
            tmp4_hi,tmp4_lo)
       call dd_add(&
            tmp3_hi,tmp3_lo,&
            tmp4_hi,tmp4_lo,&
            tmp2_hi,tmp2_lo)

       hki_hi(k+(i-1)*(m+1))=tmp1_hi
       hki_lo(k+(i-1)*(m+1))=tmp1_lo

       hki_hi(k+1+(i-1)*(m+1))=tmp2_hi
       hki_lo(k+1+(i-1)*(m+1))=tmp2_lo
    enddo
!    tmp1=hki(i+1+(i-1)*(m+1))/hki(i+(i-1)*(m+1))
    if( ( hki_hi(i+(i-1)*(m+1)).eq.0.0d0 ).and.&
        ( hki_lo(i+(i-1)*(m+1)).eq.0.0d0 ) )then
       call dd_div(&
            hki_hi(i+1+(i-1)*(m+1)),hki_lo(i+1+(i-1)*(m+1)),&
            1.0d-100,0.0d0,&
            tmp1_hi,tmp1_lo)
    else
       call dd_div(&
            hki_hi(i+1+(i-1)*(m+1)),hki_lo(i+1+(i-1)*(m+1)),&
            hki_hi(i+(i-1)*(m+1))  ,hki_lo(i+(i-1)*(m+1)),&
            tmp1_hi,tmp1_lo)
    endif
!    c(i)=1.0d0/sqrt(1.0d0+tmp1*tmp1 )
    call dd_mul(&
         tmp1_hi,tmp1_lo,&
         tmp1_hi,tmp1_lo,&
         tmp2_hi,tmp2_lo )
    call dd_add(&         
         1.0d0,0.0d0,&
         tmp2_hi,tmp2_lo,&
         tmp3_hi,tmp3_lo)
    call dd_sqrt(tmp3_hi,tmp3_lo,tmp3_hi,tmp3_lo)
    call dd_div(&
         1.0d0,0.0d0,&
         tmp3_hi,tmp3_lo,&
         c_hi(i),c_lo(i))

!    s(i)=-c(i)*tmp1
    call dd_mul(&
         c_hi(i),c_lo(i),&
         tmp1_hi,tmp1_lo,&
         s_hi(i),s_lo(i))
    call dd_mul(&
         -1.0d0,0.0d0,&
         s_hi(i),s_lo(i),&
         s_hi(i),s_lo(i))

!    tmp1=c(i)*e(i)-s(i)*e(i+1)
    call dd_mul(&
         c_hi(i),c_lo(i),&
         e_hi(i),e_lo(i),&
         tmp3_hi,tmp3_lo)
    call dd_mul(&
         s_hi(i),s_lo(i),&
         e_hi(i+1),e_lo(i+1),&
         tmp4_hi,tmp4_lo)
    call dd_sub(&
         tmp3_hi,tmp3_lo,&
         tmp4_hi,tmp4_lo,&
         tmp1_hi,tmp1_lo)

!    tmp2=s(i)*e(i)+c(i)*e(i+1)
    call dd_mul(&
         s_hi(i),s_lo(i),&
         e_hi(i),e_lo(i),&
         tmp3_hi,tmp3_lo)
    call dd_mul(&
         c_hi(i),c_lo(i),&
         e_hi(i+1),e_lo(i+1),&
         tmp4_hi,tmp4_lo)
    call dd_add(&
         tmp3_hi,tmp3_lo,&
         tmp4_hi,tmp4_lo,&
         tmp2_hi,tmp2_lo)

    ! e(i)=tmp1
    ! e(i+1)=tmp2
    e_hi(i)=tmp1_hi
    e_lo(i)=tmp1_lo
    e_hi(i+1)=tmp2_hi
    e_lo(i+1)=tmp2_lo

!    hki(i+  (i-1)*(m+1))=c(i)*hki(i+(i-1)*(m+1))-s(i)*hki(i+1+(i-1)*(m+1))
    call dd_mul(&
         c_hi(i),c_lo(i),&
         hki_hi(i+(i-1)*(m+1)),hki_lo(i+(i-1)*(m+1)),&
         tmp3_hi,tmp3_lo)    
    call dd_mul(&
         s_hi(i),s_lo(i),&
         hki_hi(i+1+(i-1)*(m+1)),hki_lo(i+1+(i-1)*(m+1)),&
         tmp4_hi,tmp4_lo)
    call dd_sub(&
         tmp3_hi,tmp3_lo,&
         tmp4_hi,tmp4_lo,&
         hki_hi(i+  (i-1)*(m+1)), hki_lo(i+  (i-1)*(m+1)))
    hki_hi(i+1+(i-1)*(m+1))=0.0
    hki_lo(i+1+(i-1)*(m+1))=0.0

  end subroutine ddGivensRotations

end module ddKSP_GMRES

subroutine ddKSP_pGMRES(pgmres)
  use krylov
  use ddkrylov_kernel
  use vec_kernel
  use ddvec_kernel
  use ddQR_kernel
  use ddKSP_GMRES

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
  integer myrank

  integer(kind=intType) itr
  integer(kind=intType) inner
  ! ----
  integer(kind=intType) i,j,k
  integer ierr

  real*8,allocatable :: x_hi(:)
  real*8,allocatable :: x_lo(:)
  real*8,allocatable :: b_hi(:)
  real*8,allocatable :: b_lo(:)
  real*8,allocatable :: r_hi(:)
  real*8,allocatable :: r_lo(:)
  real*8,allocatable :: v_hi(:)
  real*8,allocatable :: v_lo(:)
  real*8,allocatable :: pv_hi(:)
  real*8,allocatable :: pv_lo(:)
  real*8,allocatable :: w_hi(:)
  real*8,allocatable :: w_lo(:)
  real*8,allocatable :: x0_hi(:)
  real*8,allocatable :: x0_lo(:)
  real*8,allocatable :: v0_hi(:)
  real*8,allocatable :: v0_lo(:)

  integer(kind=intType) msize

  real*8 bnorm_hi,bnorm_lo
  real*8 rnorm_hi,rnorm_lo
  real*8 tmp1_hi,tmp1_lo
  real*8 tmp2_hi,tmp2_lo

  real*8,allocatable :: e_hi(:)
  real*8,allocatable :: e_lo(:)
  real*8,allocatable :: c_hi(:)
  real*8,allocatable :: c_lo(:)
  real*8,allocatable :: s_hi(:)
  real*8,allocatable :: s_lo(:)
  real*8,allocatable :: y_hi(:)
  real*8,allocatable :: y_lo(:)
  real*8,allocatable :: hki_hi(:)
  real*8,allocatable :: hki_lo(:)
  real*8,allocatable :: Rtri_hi(:)
  real*8,allocatable :: Rtri_lo(:)

  integer PRECISION_A
  integer PRECISION_b
  integer PRECISION_x

  integer(kind=intType) mm

  itrmax=pgmres%itrmax
  rtolmax=pgmres%rtolmax
  abstolmax=pgmres%abstolmax
  myrank=pgmres%commVec%myrank
  icomm=pgmres%commVec%comm
  n=pgmres%prm%n
  m=pgmres%GMRES_m
  iGSflag=pgmres%GMRES_GSflag

  PRECISION_A=pgmres%PRECISION_A
  PRECISION_b=pgmres%PRECISION_b
  PRECISION_x=pgmres%PRECISION_x

  allocate(x_hi(n))
  allocate(x_lo(n))
  allocate(b_hi(n))
  allocate(b_lo(n))
  allocate(r_hi(n))
  allocate(r_lo(n))
  allocate(v_hi(n*(m+1)))
  allocate(v_lo(n*(m+1)))
  allocate(pv_hi(n))
  allocate(pv_lo(n))
  allocate(w_hi(n))
  allocate(w_lo(n))
  allocate(x0_hi(n))
  allocate(x0_lo(n))
  allocate(v0_hi(n))
  allocate(v0_lo(n))

  allocate(e_hi(m+1))
  allocate(e_lo(m+1))
  allocate(c_hi(m+1))
  allocate(c_lo(m+1))
  allocate(s_hi(m+1))
  allocate(s_lo(m+1))
  allocate(y_hi(m+1))
  allocate(y_lo(m+1))
  allocate(hki_hi((m+1)*m ))
  allocate(hki_lo((m+1)*m ))
  allocate(Rtri_hi((m+1)*(m+1)))
  allocate(Rtri_lo((m+1)*(m+1)))

  call zero_initialize(n,x_hi)
  call zero_initialize(n,x_lo)
  call zero_initialize(n,b_hi)
  call zero_initialize(n,b_lo)
  call zero_initialize(n,r_hi)
  call zero_initialize(n,r_lo)

  msize=n*(m+1)
  call zero_initialize(msize,v_hi)
  call zero_initialize(msize,v_lo)
  call zero_initialize(n,pv_hi)
  call zero_initialize(n,pv_lo)
  call zero_initialize(n,w_hi)
  call zero_initialize(n,w_lo)
  call zero_initialize(n,x0_hi)
  call zero_initialize(n,x0_lo)
  call zero_initialize(n,v0_hi)
  call zero_initialize(n,v0_lo)

  call xcpy(n,pgmres%x_hi,x_hi)
  if(PRECISION_x.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
     call xcpy(n,pgmres%x_lo,x_lo)
  endif
  call xcpy(n,pgmres%b_hi,b_hi)
  if(PRECISION_b.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
     call xcpy(n,pgmres%b_lo,b_lo)
  endif
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

  call ddcalc_dot(n,icomm,&
       b_hi,b_lo,&
       b_hi,b_lo,&
       bnorm_hi,bnorm_lo)
  call dd_sqrt(bnorm_hi,bnorm_lo,bnorm_hi,bnorm_lo)

  do itr = 1,itrmax,m
     ! -----
     e_hi=0.0d0
     e_lo=0.0d0
     c_hi=0.0d0
     c_lo=0.0d0
     s_hi=0.0d0
     s_lo=0.0d0
     y_hi=0.0d0
     y_lo=0.0d0
     hki_hi=0.0d0
     hki_lo=0.0d0
     Rtri_hi=0.0d0
     Rtri_lo=0.0d0

     ! r=b-Ax
     call ddcalc_res(& 
          pgmres, &
          x_hi,x_lo, &
          b_hi,b_lo, &
          r_hi,r_lo)
     call ddcalc_dot(n,icomm,&
          r_hi,r_lo, &
          r_hi,r_lo, &
          rnorm_hi,rnorm_lo)
!     rnorm=sqrt(rnorm)
     call dd_sqrt(rnorm_hi,rnorm_lo,rnorm_hi,rnorm_lo)

     call dd_div(&
          -1.0d0,0.0d0, &
          rnorm_hi,rnorm_lo,&
          tmp1_hi,tmp1_lo)
     tmp2_hi=0.0d0
     tmp2_lo=0.0d0

     ! v=r/||r||
     call ddaxpy(n,&
          tmp1_hi,tmp1_lo, &
          tmp2_hi,tmp2_lo, &
          r_hi,r_lo,&
          v_hi,v_lo)
     ! e(1)=-rnorm
     call dd_mul(-1.0d0,0.0d0,rnorm_hi,rnorm_lo,e_hi(1),e_lo(1))

     if( (rnorm_hi<abstolmax) ) then
        exit
     endif
     if(itr.eq.1)then
        pgmres%reshistory(itr)=rnorm_hi/bnorm_hi
     endif
#if 0
     if( (myrank.eq.0) )then
#else
     if( (myrank.eq.0).and.(itr.eq.1) )then
#endif
!        write(*,*)"@ res @",itr,rnorm_hi,rnorm_hi/bnorm_hi
     endif
     ! inner ループ
     do inner=1,m
        ! Mpv=v
        call ddcalc_precon(pgmres,&
             v_hi(1+n*(inner-1)),v_lo(1+n*(inner-1)),&
             pv_hi,pv_lo)
        ! w=Apv
        call ddcalc_Matvec( pgmres,&
             pv_hi,pv_lo,&
             w_hi,w_lo)
        ! 直交化(古典,修正グラムシュミット)
        ! グラムシュミット直交化
        if(iGSflag.eq.KSP_GMRES_GS_TYPE_MGS)then
           call ddMGS_Orthonormalization( &
                n, &
                m, &
                inner, &
                icomm, &
                w_hi,w_lo, &
                v_hi,v_lo, &
                Rtri_hi,Rtri_lo )
        else
           call ddCGS_Orthonormalization( &
                n, &
                m, &
                inner, &
                icomm, &
                w_hi,w_lo, &
                v_hi,v_lo, &
                Rtri_hi,Rtri_lo )
        endif
        do k=1,inner-1
           hki_hi(k + (inner-1)*(m+1))=Rtri_hi(k+m*(inner-1))
           hki_lo(k + (inner-1)*(m+1))=Rtri_lo(k+m*(inner-1))
        enddo

        call ddcalc_dot(n,icomm ,&
             w_hi,w_lo,&
             v_hi(1+n*(inner-1)),v_lo(1+n*(inner-1))   , &
             hki_hi(inner+(inner-1)*(m+1)),hki_lo(inner+(inner-1)*(m+1)) )

        call dd_mul(-1.0d0,0.0d0,&
             hki_hi(inner+(inner-1)*(m+1)),hki_lo(inner+(inner-1)*(m+1)),&
             tmp1_hi,tmp1_lo)

        tmp2_hi=1.0d0
        tmp2_lo=0.0d0
        call ddaxpy(n,&
             tmp1_hi,tmp1_lo, &
             tmp2_hi,tmp2_lo, &
             v_hi(1+n*(inner-1)),v_lo(1+n*(inner-1)), &
             w_hi,w_lo)

        ! h(k+1i)=||w||
        call ddcalc_dot(n,icomm ,&
             w_hi,w_lo, &
             w_hi,w_lo, &
             hki_hi(inner+1 + (inner-1)*(m+1)),hki_lo(inner+1 + (inner-1)*(m+1)) )
!        hki(inner+1 + (inner-1)*(m+1))  = &
!             sqrt(hki(inner+1 + (inner-1)*(m+1)) )
     call dd_sqrt(&
          hki_hi(inner+1 + (inner-1)*(m+1)),hki_lo(inner+1 + (inner-1)*(m+1)),&
          hki_hi(inner+1 + (inner-1)*(m+1)),hki_lo(inner+1 + (inner-1)*(m+1)))

        ! v_i+1 = w/||w|| 
!        tmp1=1.0d0/hki(inner+1 + (inner-1)*(m+1)) 
     call dd_div(&
          1.0d0,0.0d0, &
          hki_hi(inner+1 + (inner-1)*(m+1)),hki_lo(inner+1 + (inner-1)*(m+1)),  &
          tmp1_hi,tmp1_lo)
     tmp2_hi=0.0d0
     tmp2_lo=0.0d0
     call ddaxpy(n,&
          tmp1_hi,tmp1_lo,&
          tmp2_hi,tmp2_lo,&
          w_hi,w_lo,&
          v_hi(1+n*(inner-1+1)),v_lo(1+n*(inner-1+1)))

        ! -----
        ! ギブンス回転
        call ddGivensRotations(&
             m,&
             hki_hi,hki_lo,&
             inner,&
             c_hi,c_lo,&
             s_hi,s_lo,&
             e_hi,e_lo)
        !      //      printf("   %4d %4d %16.14e \n",itr,inner,fabs(e[inner+1]));
        rnorm_hi=abs(e_hi(inner+1))

        ! ------
        mm=inner
        pgmres%reshistory(itr+inner)=rnorm_hi/bnorm_hi 
        if(myrank.eq.0)then
!           write(*,*)"@ res @",itr+inner,rnorm_hi,rnorm_hi/bnorm_hi
        endif
        if( (rnorm_hi/bnorm_hi<rtolmax) ) then
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

     !   解ベクトル更新
     call ddupdate_xVec(pgmres,n,m,mm,&
          x_hi,x_lo,&
          v_hi,v_lo,&
          x0_hi,x0_lo,&
          v0_hi,v0_lo,&
          y_hi,y_lo,&
          hki_hi,hki_lo,&
          e_hi,e_lo)

     if( (rnorm_hi/bnorm_hi<rtolmax) ) then
        exit
     endif
     if( (itr+inner).eq.itrmax)then
        exit
     endif
  enddo
  ! ------
  ! 真の残差を確認
#if 0
  !call ddcalc_res(pgmres,x,b,r)
  !call ddcalc_dot(n,icomm,r,r,rnorm)
  call ddcalc_res(pgmres,x_hi,x_lo,b_hi,b_lo,r_hi,r_lo)
  call ddcalc_dot(n,icomm,r_hi,r_lo,r_hi,r_lo,rnorm_hi,rnorm_lo)
  
  call dd_sqrt(rnorm_hi,rnorm_lo,rnorm_hi,rnorm_lo)
  if(myrank.eq.0)then
     write(*,*)"    @ true res @",rnorm_hi,rnorm_hi/bnorm_hi
  endif
#endif
  call xcpy(n,x_hi,pgmres%x_hi)
  if(PRECISION_x.eq.KSP_DOUBLE_DOUBLE_PRECISION)then
     call xcpy(n,x_lo,pgmres%x_lo)
  endif

  pgmres%niter=min(itrmax,itr+inner)
  if(itrmax<=itr)then
     pgmres%iret=PARCEL_RETURNCODE_MAXITERATION
  else
     pgmres%iret=PARCEL_RETURNCODE_NORMAL
  endif

  deallocate(x_hi)
  deallocate(x_lo)
  deallocate(b_hi)
  deallocate(b_lo)
  deallocate(r_hi)
  deallocate(r_lo)
  deallocate(v_hi)
  deallocate(v_lo)
  deallocate(pv_hi)
  deallocate(pv_lo)
  deallocate(w_hi)
  deallocate(w_lo)
  deallocate(x0_hi)
  deallocate(x0_lo)
  deallocate(v0_hi)
  deallocate(v0_lo)

  deallocate(e_hi)
  deallocate(e_lo)
  deallocate(c_hi)
  deallocate(c_lo)
  deallocate(s_hi)
  deallocate(s_lo)
  deallocate(y_hi)
  deallocate(y_lo)
  deallocate(hki_hi)
  deallocate(hki_lo)
  deallocate(Rtri_hi)
  deallocate(Rtri_lo)

end subroutine ddKSP_pGMRES

