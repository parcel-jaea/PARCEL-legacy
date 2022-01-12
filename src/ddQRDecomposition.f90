module ddQR_kernel
  use PARCEL_define
contains
  subroutine ddMGS_Orthonormalization( &
       n,s,i,comm,&
       V_hi,V_lo,&
       Q_hi,Q_lo,&
       R_hi,R_lo)
    use PARCEL_integer
    use ddvec_kernel
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    integer(kind=intType) i
    integer comm
    real*8 V_hi(n)
    real*8 V_lo(n)
    real*8 Q_hi(n*s)
    real*8 Q_lo(n*s)
    real*8 R_hi(s*s)
    real*8 R_lo(s*s)
    real*8 tmp1_hi,tmp1_lo
    real*8 tmp2_hi,tmp2_lo
    integer(kind=intType) k
    ! -----
    ! 修正グラムシュミット
    do k=1,i-1
       ! Rki = <Qk,Vi>
       call ddcalc_dot(n,comm ,&
            V_hi,V_lo,&
            Q_hi(1+n*(k-1)) ,Q_lo(1+n*(k-1)), &
            R_hi(k+s*(i-1)) ,R_lo(k+s*(i-1)) )
       ! Vi=Vi-Rki*Qk
       call dd_mul(&
            -1.0d0,0.0d0,&
            R_hi(k+s*(i-1)),R_lo(k+s*(i-1)),&
            tmp1_hi,tmp1_lo)
       tmp2_hi=1.0d0
       tmp2_lo=0.0d0
       call ddaxpy(n,&
            tmp1_hi,tmp1_lo, &
            tmp2_hi,tmp2_lo,&
            Q_hi(1+n*(k-1)),Q_lo(1+n*(k-1)),&
            V_hi,V_lo)
    enddo

  end subroutine ddMGS_Orthonormalization

  subroutine ddCGS_Orthonormalization( &
       n,s,i,comm,&
       V_hi,V_lo, &
       Q_hi,Q_lo, &
       R_hi,R_lo,nvecBlock)
    use PARCEL_integer
    use ddvec_kernel
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    integer(kind=intType) i
    integer comm
    real*8 V_hi(n)
    real*8 V_lo(n)
    real*8 Q_hi(n*s)
    real*8 Q_lo(n*s)
    real*8 R_hi(s*s)
    real*8 R_lo(s*s)
    real*8 tmp1_hi,tmp1_lo
    real*8 tmp2_hi,tmp2_lo
    integer(kind=intType) k
    integer(kind=intType) j,jj
    integer(kind=intType) nblock
    integer nvecBlock
    ! -----
    ! 古典グラムシュミット
    ! Rki = <Qk,Vi>
    if(i>0)then
       j = 1
       call ddcalc_gemm_dot(n,i,j,comm,&
            Q_hi,Q_lo,&
            V_hi,V_lo,&
            R_hi(1+s*(i-1)),R_lo(1+s*(i-1)),nvecBlock)
    endif

#if 1
    do k=1,i-1
       ! Vi=Vi-Rki*Qk
       call dd_mul(&
            -1.0d0,0.0d0,&
            R_hi(k+s*(i-1)),R_lo(k+s*(i-1)),&
            tmp1_hi,tmp1_lo)
       tmp2_hi=1.0d0
       tmp2_lo=0.0d0
       call ddaxpy(n,&
            tmp1_hi,tmp1_lo,&
            tmp2_hi,tmp2_lo,&
            Q_hi(1+n*(k-1)),Q_lo(1+n*(k-1)),&
            V_hi,V_lo)
    enddo
#else

    do jj=1,n,nvecBlock 
       nblock=min(n,jj+nvecBlock-1)-jj+1
       do k=1,i-1
          ! Vi=Vi-Rki*Qk
          call dd_mul(&
               -1.0d0,0.0d0,&
               R_hi(k+s*(i-1)),R_lo(k+s*(i-1)),&
               tmp1_hi,tmp1_lo)
          tmp2_hi=1.0d0
          tmp2_lo=0.0d0
          call ddaxpy(nblock,&
               tmp1_hi,tmp1_lo,&
               tmp2_hi,tmp2_lo,&
               Q_hi(jj+n*(k-1)),Q_lo(jj+n*(k-1)),&
               V_hi(jj),V_lo(jj))
       enddo
    enddo

#endif

  end subroutine ddCGS_Orthonormalization

  subroutine ddQR_GS(&
       n,s,comm,&
       V_hi,V_lo,Q_hi,Q_lo,R_hi,R_lo,&
       iGSflag,nvecBlock)
    use PARCEL_integer
    use ddvec_kernel
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    integer comm
    real*8 V_hi(n*s)
    real*8 V_lo(n*s)
    real*8 Q_hi(n*s)
    real*8 Q_lo(n*s)
    real*8 R_hi(s*s)
    real*8 R_lo(s*s)
    integer iGSflag

    integer(kind=intType) i,j
    real*8 tmp1_hi,tmp1_lo
    real*8 tmp2_hi,tmp2_lo
    integer nvecBlock
    
    do i=1,s*s
       R_hi(i)=0.0d0
       R_lo(i)=0.0d0
    enddo

    do i=1,s
       if(iGSflag.eq.KSP_QRDECOMPOSITION_MGS)then
          call ddMGS_Orthonormalization(&
               n,s,i,comm,&
               V_hi(1+n*(i-1)),V_lo(1+n*(i-1)),&
               Q_hi,Q_lo,&
               R_hi,R_lo)
       else
          call ddCGS_Orthonormalization(&
               n,s,i,comm,&
               V_hi(1+n*(i-1)),V_lo(1+n*(i-1)),&
               Q_hi,Q_lo,&
               R_hi,R_lo,nvecBlock)
       endif
       call ddcalc_dot(&
            n,comm,&
            V_hi(1+n*(i-1)),V_lo(1+n*(i-1)),&
            V_hi(1+n*(i-1)),V_lo(1+n*(i-1)),&
            R_hi(i+s*(i-1)),R_lo(i+s*(i-1)))
       call dd_sqrt(&
            R_hi(i+s*(i-1)),R_lo(i+s*(i-1)),&
            R_hi(i+s*(i-1)),R_lo(i+s*(i-1)))

       call dd_div(&
            1.0d0,0.0d0, &
            R_hi(i+s*(i-1)),R_lo(i+s*(i-1)),&
            tmp1_hi,tmp1_lo)
       tmp2_hi=0.0d0
       tmp2_lo=0.0d0

       call ddaxpy(n,&
            tmp1_hi,tmp1_lo,&
            tmp2_hi,tmp2_lo,&
            V_hi(1+n*(i-1)),V_lo(1+n*(i-1)),&
            Q_hi(1+n*(i-1)),Q_lo(1+n*(i-1)))
    enddo

  end subroutine ddQR_GS

end module ddQR_kernel
