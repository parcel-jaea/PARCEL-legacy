module QR_kernel
  use TSQR_kernel
contains

  subroutine MGS_Orthonormalization(n,s,i,icomm,V,Q,R)
    use vec_kernel
            use PARCEL_integer
implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    integer(kind=intType) i
    integer icomm
    real*8 V(n)
    real*8 Q(n*s)
    real*8 R(s*s)
    real*8 tmp1,tmp2
    integer(kind=intType) k
    ! -----
    ! 修正グラムシュミット
    do k=1,i-1
       ! Rki = <Qk,Vi>
       call calc_dot(n,icomm ,&
            V,&
            Q(1+n*(k-1))   , &
            R(k+s*(i-1)) )
       ! Vi=Vi-Rki*Qk
       tmp1=-R(k+s*(i-1))
       tmp2=1.0d0
       call axpy(n,tmp1,tmp2,Q(1+n*(k-1)),V)
    enddo

  end subroutine MGS_Orthonormalization

  subroutine CGS_Orthonormalization(n,s,i,icomm,V,Q,R)
    use vec_kernel
            use PARCEL_integer
implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    integer(kind=intType) i
    integer icomm
    real*8 V(n)
    real*8 Q(n*s)
    real*8 R(s*s)
    real*8 tmp1,tmp2
    integer(kind=intType) k
    integer(kind=intType) j,jj
    integer(kind=intType) nblock
    ! -----
    ! 古典グラムシュミット
    ! Rki = <Qk,Vi>
    if(i>0)then
       call calc_gemm_dot(n,i,1,icomm,Q,V,R(1+s*(i-1)))
    endif

#if 0
    do k=1,i-1
       ! Vi=Vi-Rki*Qk
       tmp1=-R(k+s*(i-1))
       tmp2=1.0d0
       call axpy(n,tmp1,tmp2,Q(1+n*(k-1)),V)
    enddo
#else

#if 0
    do jj=1,n,nvecBlock 
       nblock=min(n,jj+nvecBlock-1)-jj+1
       do k=1,i-1
          ! Vi=Vi-Rki*Qk
          tmp1=-R(k+s*(i-1))
          tmp2=1.0d0
          call axpy(nblock,tmp1,tmp2,Q(jj+n*(k-1)),V(jj) )
       enddo
    enddo
#else

    !$OMP PARALLEL DO PRIVATE(jj,nblock,k,j,tmp1,tmp2) 
    do jj=1,n,nvecBlock 
       nblock=min(n,jj+nvecBlock-1)-jj+1
       do k=1,i-1
          ! Vi=Vi-Rki*Qk
!           tmp1=-R(k+s*(i-1))
!           tmp2=1.0d0
!           call axpy(nblock,tmp1,tmp2,Q(jj+n*(k-1)),V(jj) )
          do j=jj,jj+nblock-1
             V(j)=V(j)-R(k+s*(i-1))*Q(j+n*(k-1))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO 
#endif

#endif

  end subroutine CGS_Orthonormalization

  subroutine GS_QR(n,s,icomm,V,Q,R,iGSflag)
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

    integer(kind=intType) i,j
    real*8 tmp1,tmp2

    do i=1,s*s
       R(i)=0.0d0
    enddo

    ! todo 古典グラムは内積よりも行列積で実装したほうが高速
    do i=1,s
       if(iGSflag.eq.KSP_QRDECOMPOSITION_MGS)then
          call MGS_Orthonormalization(n,s,i,icomm,V(1+n*(i-1)),Q,R)
       else
          call CGS_Orthonormalization(n,s,i,icomm,V(1+n*(i-1)),Q,R)
       endif
       call calc_dot(n,icomm,V(1+n*(i-1)),V(1+n*(i-1)),R(i+s*(i-1)))
       R(i+s*(i-1))=sqrt(R(i+s*(i-1)))
       tmp1=1.0d0/R(i+s*(i-1))
       tmp2=0.0d0
       call axpy(n,tmp1,tmp2,V(1+n*(i-1)),Q(1+n*(i-1)))
    enddo

  end subroutine GS_QR

  subroutine calcQ_VRinv(n,s,V,Q,R)
    use vec_kernel ! :only,nvecBlock
            use PARCEL_integer
implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    real*8 V(n*s)
    real*8 Q(n*s)
    real*8 R(s*s)
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

  subroutine CholeskyQR_Ronly(n,s,icomm,V,R)
    use mpi
    use vec_kernel
            use PARCEL_integer
implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    integer icomm
    real*8 V(n*s)
    real*8 R(s*s)

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
    call calc_syrk_dot_local(n,s,V,B)
#endif
    ! ALLREDUCE
    ! 三角行列なので(s+1)*s/2のデータ数を送受信すれば良い
    ss_int4=s*s
    call MPI_Allreduce(B,R, ss_int4,MPI_REAL8,MPI_SUM,icomm,ierr)    

    ! コレスキー分解
    s_int4=s
    call DPOTRF('U', s_int4,R, s_int4,ierr)

  end subroutine CholeskyQR_Ronly

  subroutine TSQR_Ronly(n,s,icomm,V,R)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    integer icomm
    real*8 V(n*s)
    real*8 R(s*s)
    integer(kind=intType) i,j
    call TSQR_local_Ronly(n,s,V,R)
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

  subroutine CholeskyQR2_Ronly(n,s,icomm,V,R)
            use PARCEL_integer
implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    integer icomm
    real*8 V(n*s)
    real*8 R(s*s)
    real*8 Q1(n*s)
    real*8 R2(s*s)
    integer s_int4
    call CholeskyQR_Ronly(n,s,icomm,V,R)
    call calcQ_VRinv(n,s,V,Q1,R)
    call CholeskyQR_Ronly(n,s,icomm,Q1,R2)

    s_int4=s
    call DTRMM('L','U','N','N', s_int4, s_int4,1.0d0,R2,s_int4,R,s_int4)

  end subroutine CholeskyQR2_Ronly

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
       write(*,*)"|| QtQ - I ||           = ",Orthnorm
       write(*,*)"|| QtQ - I || / sqrt(n*s) = ", Orthnorm/sqrt(n_)
    endif

  end subroutine check_QtQ

  subroutine QRDecomposition(n,s,icomm,V,Q,R,iQRflag,iQmatflag)
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
    integer iQmatflag
    integer s_int4
    s_int4=s
    if(s_int4.eq.1)then
       ! 1本のベクトルに対する直交化なので正規化のみ行う
       call GS_QR(n,s,icomm,V,Q,R,KSP_QRDECOMPOSITION_MGS)
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
       call GS_QR(n,s,icomm,V,Q,R,KSP_QRDECOMPOSITION_MGS)
       return
    case(KSP_QRDECOMPOSITION_CGS)
!       write(*,*)"CGS " 
       call GS_QR(n,s,icomm,V,Q,R,KSP_QRDECOMPOSITION_CGS)
       return
    case(KSP_QRDECOMPOSITION_TSQR)
!       write(*,*)"TSQR " 
       if(n>s)then
          call TSQR_Ronly(n,s,icomm,V,R)
       else
          ! n<s は縦長行列ではないためTSQRで計算することは不可
          ! n<s の場合はCholeskyQR2で計算する
          call CholeskyQR2_Ronly(n,s,icomm,V,R)
       endif
    case(KSP_QRDECOMPOSITION_CholeskyQR)
!       write(*,*)"CholeskyQR " 
      call CholeskyQR_Ronly(n,s,icomm,V,R)
    case(KSP_QRDECOMPOSITION_CholeskyQR2)
!       write(*,*)"CholeskyQR2 " 
       call CholeskyQR2_Ronly(n,s,icomm,V,R)
    end select
    if(iQmatflag.eq.KSP_QRDECOMPOSITION_TRUE_ORTHOGONAL_MATRIX)then
       call calcQ_VRinv(n,s,V,Q,R)
    endif

  end subroutine QRDecomposition

  subroutine QR_TEST(n,s,icomm,Vorg)
    use vec_kernel
    use PARCEL_define
            use PARCEL_integer
implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    integer icomm
    real*8 Vorg(n*s)

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

       call QRDecomposition(n,s,icomm,V,Q,R,iQRflag,iQmatflag)
       call check_QtQ(n,s,icomm,Q)
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
  call axpy(n,tmp1,tmp2,r,V)

!  call write_vec(s,V)
  do i=1,s
     ! Mpv=v_i
     call calc_precon(method,V(1+(i-1)*n),pV )
     ! V_i+1=ApV
     call calc_Matvec(method,pv,V(1+i*n) )
!     call write_vec(s,V(1+i*n))
  enddo
  do i=1,s+1
     call QR_TEST(n,i,icomm,V)
  enddo

  deallocate(x)
  deallocate(r)
  deallocate(pv)
  deallocate(V)

end subroutine MPK_QR_TEST
