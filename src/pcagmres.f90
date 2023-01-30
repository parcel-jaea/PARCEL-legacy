! #define timer 1
module KSP_GMRES_CPU

contains
  subroutine CGS_Orth_FP64_CPU(&
       icomm, &
       w,ldw, &
       v,ldv, &
       hki,&
       n, &
       vrow,  &
       nvecBlock &
       )
    use vec_kernel,only : calc_gemm_dot_CPU,gemv
    use PARCEL_integer,only : intType
    implicit none

    integer(kind=intType) vrow
    integer(kind=intType) ldv
    integer(kind=intType) ldw
    integer icomm
    real*8 w(ldw)
    real*8 v(vrow*ldv)
    real*8 hki(vrow)
    real*8 hki_(vrow)
    integer(kind=intType) n
    integer nvecBlock
    integer i
    integer(kind=intType) parcel_int

    parcel_int = 1
    call calc_gemm_dot_CPU(&
         n,parcel_int,vrow,&
         icomm,&
         w, ldw,&
         v, ldv,&
         hki, parcel_int,&
         nvecBlock)

    do i = 1,vrow
       hki_(i) = -hki(i)
    enddo
    call gemv(&
         n, vrow,&
         v, ldv,&
         w, ldw, &
         hki_,nvecBlock)

!  subroutine gemv(n,m,v,v0,y,nvecBlock)

  end subroutine CGS_Orth_FP64_CPU

!  subroutine sum_vec_FP64_CPU(n,m,vecLen,v,v0,y,nvecBlock)
  subroutine sum_vec_FP64_CPU( &
       n,m, &
       v , ldv, &
       v0, ldv0, &
       y , ldy, &
       nvecBlock)
    use vec_kernel
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) m
    integer(kind=intType) ldv
    integer(kind=intType) ldv0
    integer(kind=intType) ldy
    real*8 v(ldv*(m+1))
    real*8 v0(ldv0)
    real*8 y(ldy)
    integer nvecBlock

    real*8 tmp1,tmp2
    integer(kind=intType) i,j,jj,k

    ! v0 = sum_k y_k*v(k)
    !$OMP PARALLEL DO PRIVATE(jj,j,i) 
    do jj=1,n,nvecBlock 
       do j=jj,min(n,jj+nvecBlock-1)
          v0(j)=0.0d0
       enddo
       do i=1,m
          do j=jj,min(n,jj+nvecBlock-1)
             v0(j)=v0(j)+y(i)*v(j+ldv*(i-1))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO 

  end subroutine sum_vec_FP64_CPU
  
  !  subroutine update_xVec(pgmres,n,m0,m,x,v,x0,v0,y,hki,e)
  subroutine update_xVec_CPU( &
       pgmres,n,m0,m,&
       x, ldx,&
       v, ldv,&
       x0, ldx0,&
       v0, ldv0,&
       hki, ldhki,&
       e, lde,&
       nvecBlock)
    use krylov
    use krylov_kernel
    use vec_kernel
    use PARCEL_integer
    implicit none

    type(KSP) :: pgmres
    integer(kind=intType) n
    integer(kind=intType) m
    integer(kind=intType) m0
    integer(kind=intType) ldx
    integer(kind=intType) ldv
    integer(kind=intType) ldx0
    integer(kind=intType) ldv0
    integer(kind=intType) ldhki
    integer(kind=intType) lde
    real*8 x(ldx)
    real*8 v(ldv*(m0+1))
    real*8 x0(ldx0) 
    real*8 v0(ldv0)
!    real*8 hki( (m0+1)*m0 )
    real*8 hki( ldhki*m0 )
!    real*8 e(m0+1)
    real*8 e(lde)

    real*8 y(m0+1)
    integer nvecBlock

    real*8 Hinv( m*m )
    real*8 tmp1,tmp2
    integer(kind=intType) i,j
    integer ierr

    ! y = H^-1 e 
#if 0
    do i=1,m
       do j=1,m
          Hinv(i+(j-1)*m)=hki(i+(j-1)*(m0+1))
       enddo
    enddo
    call DTRTRI('U', 'N', m, Hinv, m, ierr)
    do i=1,m
       do j=1,m
          hki(j+(i-1)*(m0+1))=Hinv(j+(i-1)*m)
       enddo
    enddo

    do i=1,m
       y(i)=0.0d0
       do j=1,m
          y(i)=y(i)+hki(i+(j-1)*(m0+1))*e(j)            
       enddo
       !write(*,*)"y = ",i,y(i)
    enddo
#else
    do i=m,1,-1
       y(i)=e(i)
       do j=m,i+1,-1
          y(i)=y(i)-hki(i+(j-1)*(m0+1))*y(j)
       enddo
       if(hki(i+(i-1)*(m0+1)).eq.0.0d0)then
          y(i)=y(i)/1.0d-300
       else
          y(i)=y(i)/hki(i+(i-1)*(m0+1))
       endif
    enddo
#endif

    ! v0 = sum_k y_k*v(k)
!    call sum_vec_FP64_CPU(n,m,v,v0,y,nvecBlock)
    ! ここはvecLenに変更する
    call sum_vec_FP64_CPU(&
         n, m, &
         v, ldv,&
         v0, ldv0, &
         y , m0+1 , &
         nvecBlock)

    ! solve
    call calc_precon(pgmres,v0,x0 )
    !    call write_vec(n,x0)
    tmp1=1.0
    tmp2=1.0
    ! x = x+x0
    call axpy(n,tmp1,tmp2,x0,x)

  end subroutine update_xVec_CPU

  subroutine GivensRotations(m,hki,i,c,s,e)
    use PARCEL_integer
    implicit none
    integer(kind=intType) m
    integer(kind=intType) i
    real*8 hki((m+1)*m)
    real*8 tmp1,tmp2
    real*8 c(m+1)
    real*8 s(m+1)
    real*8 e(m+1)

    integer(kind=intType) k

#if 0
    do k=1,i-1
       tmp1=hki(k+(i-1)*(m+1))*c(k)-hki(k+1+(i-1)*(m+1))*s(k)
       tmp2=hki(k+(i-1)*(m+1))*s(k)+hki(k+1+(i-1)*(m+1))*c(k)
       hki(k+(i-1)*(m+1))=tmp1
       hki(k+1+(i-1)*(m+1))=tmp2
    enddo
    tmp1=hki(i+1+(i-1)*(m+1))/hki(i+(i-1)*(m+1))
    c(i)=1.0d0/sqrt(1.0d0+tmp1*tmp1 )
    s(i)=-c(i)*tmp1

    tmp1=c(i)*e(i)-s(i)*e(i+1)
    tmp2=s(i)*e(i)+c(i)*e(i+1)
    e(i)=tmp1
    e(i+1)=tmp2
    hki(i+  (i-1)*(m+1))=c(i)*hki(i+(i-1)*(m+1))-s(i)*hki(i+1+(i-1)*(m+1))
    hki(i+1+(i-1)*(m+1))=0.0
#else
    do k=1,i-1
       tmp1=hki(k+(i-1)*(m+1))*c(k)-hki(k+1+(i-1)*(m+1))*s(k)
       tmp2=hki(k+(i-1)*(m+1))*s(k)+hki(k+1+(i-1)*(m+1))*c(k)
       hki(k+(i-1)*(m+1))=tmp1
       hki(k+1+(i-1)*(m+1))=tmp2
    enddo
    if(hki(i+(i-1)*(m+1)).eq.0.0d0)then
       tmp1=hki(i+1+(i-1)*(m+1))/1.0d0-300
    else
       tmp1=hki(i+1+(i-1)*(m+1))/hki(i+(i-1)*(m+1))
    endif
    c(i)=1.0d0/sqrt(1.0d0+tmp1*tmp1 )
    s(i)=-c(i)*tmp1

    tmp1=c(i)*e(i)-s(i)*e(i+1)
    tmp2=s(i)*e(i)+c(i)*e(i+1)
    e(i)=tmp1
    e(i+1)=tmp2
    hki(i+  (i-1)*(m+1))=c(i)*hki(i+(i-1)*(m+1))-s(i)*hki(i+1+(i-1)*(m+1))
    hki(i+1+(i-1)*(m+1))=0.0
#endif
  end subroutine GivensRotations

end module KSP_GMRES_CPU

module KSP_CAARNOLDI

contains


  subroutine checkEigenPair(pcaarnoldi,n,s,Q,ER,EI,maxEigenID,minEigenID,err_norm_all)
    use krylov
    use krylov_kernel,only : calc_Matvec,calc_precon
    use vec_kernel,only : zero_initialize,axpy,calc_dot
    use PARCEL_integer
    implicit none
    type(KSP) :: pcaarnoldi
    integer(kind=intType) n
    integer(kind=intType) s

    real*8 Q(n*s)
    real*8 ER(s)
    real*8 EI(s)
    integer maxEigenID,minEigenID
    real*8,allocatable :: px(:)
    real*8,allocatable :: px_im(:)
    real*8,allocatable :: Ax(:)
    real*8,allocatable :: Ax_im(:)
    real*8,allocatable :: kx(:)
    real*8,allocatable :: kx_im(:)
    real*8 err_norm_all(s)
    real*8 err_norm,err_im_norm
    real*8 kx_norm,kx_im_norm
    real*8 maxEigen,minEigen,tmp
    integer i
    integer icomm 
    integer myrank

    icomm  = pcaarnoldi%commVec%comm
    myrank = pcaarnoldi%commVec%myrank

    ! 固有対の精度チェックは||Ax-kx||/||kx||により行う.kは固有値.xは固有ベクトル.
    allocate(Ax(n))
    allocate(kx(n))
    allocate(px(n))
    allocate(px_im(n))
    allocate(Ax_im(n))
    allocate(kx_im(n))
    call zero_initialize(n,px)
    call zero_initialize(n,Ax)
    call zero_initialize(n,kx)
    call zero_initialize(n,px_im)
    call zero_initialize(n,Ax_im)
    call zero_initialize(n,kx_im)

    err_norm_all = 0.0d0
#if 0
    maxEigenID = 1
    minEigenID = 1
    maxEigen = EI(1) * EI(1) + ER(1) * ER(1) 
    minEigen = EI(1) * EI(1) + ER(1) * ER(1) 
    do i=1,s
       tmp = EI(i) * EI(i) + ER(i) * ER(i)
       if(maxEigen<tmp)then
!          if( ER(i)>0.0d0)then
             maxEigen = tmp
             maxEigenID = i
!          endif
       endif
       if(minEigen>tmp)then
!          if( ER(i)<0.0d0)then
             minEigen = tmp
             minEigenID = i
!          endif
       endif
    enddo
#endif

#if 0
    if(myrank.eq.0)then
       write(*,*)"  ||Ax-kx||/||kx||    EigenValue     "
       write(*,*)
    endif
#endif
     do i=1,s
       if(EI(i).eq.0.0d0)then

          ! 固有値が実数
          call calc_precon(pcaarnoldi,Q(1+n*(i-1)),px )
          call calc_Matvec(pcaarnoldi,px,Ax )
          call axpy(n,ER(i),0.0d0,Q(1+n*(i-1)),kx)
          call calc_dot(n,icomm,kx,kx,kx_norm)
          call axpy(n,1.0d0,-1.0d0,Ax,kx)
          call calc_dot(n,icomm,kx,kx,err_norm)
          err_norm = sqrt(abs(err_norm))
          kx_norm = sqrt(abs(kx_norm))
          err_norm_all(i) = err_norm/kx_norm
#if 0
          if(myrank.eq.0)then
             write(*,'(i8,e18.8,e18.8)')i,err_norm/kx_norm,ER(i)
          endif
#endif
       else
          if(0.0d0>EI(i))then
             cycle
          endif
          ! 固有値が共役複素数
          ! 行列Aが実数なのでここはcalc_Matvecが利用可能
#if 1
          call calc_precon(pcaarnoldi,Q(1+n*(i-1)),px )
          call calc_Matvec(pcaarnoldi,px,Ax )
          call calc_precon(pcaarnoldi,Q(1+n*(i-1+1)),px_im )
          call calc_Matvec(pcaarnoldi,px_im,Ax )
#else
          call calc_Matvec( pcaarnoldi,Q(1+n*(i-1)  ),Ax )
          call calc_Matvec( pcaarnoldi,Q(1+n*(i-1+1)),Ax_im )
#endif
          ! 実数部の計算
          call axpy(n, ER(i),0.0d0,Q(1+n*(i-1)  ),kx)
          call axpy(n,-EI(i),1.0d0,Q(1+n*(i-1)+1),kx)
          ! 虚数部の計算
          call axpy(n, EI(i),0.0d0,Q(1+n*(i-1)  ),kx_im)
          call axpy(n, ER(i),1.0d0,Q(1+n*(i-1)+1),kx_im)
          call calc_dot(n,icomm,kx,kx,kx_norm)
          call calc_dot(n,icomm,kx_im,kx_im,kx_im_norm)
          kx_norm = kx_norm + kx_im_norm
          kx_norm = sqrt(abs(kx_norm))

          call axpy(n,1.0d0,-1.0d0,Ax,kx)
          call axpy(n,1.0d0,-1.0d0,Ax_im,kx_im)

          call calc_dot(n,icomm,kx,kx,err_norm)
          call calc_dot(n,icomm,kx_im,kx_im,kx_im_norm)
          err_norm = err_norm + kx_im_norm
          err_norm = sqrt(abs(err_norm))

          err_norm_all(i) = err_norm/kx_norm
#if 0
          if(myrank.eq.0)then
             write(*,'(i8,e18.8,e18.8,e18.8)')i,err_norm/kx_norm,ER(i),EI(i)
          endif
#endif
          ! -------------
#if 1
          call calc_precon(pcaarnoldi,Q(1+n*(i-1)),px )
          call calc_Matvec(pcaarnoldi,px,Ax )
          call calc_precon(pcaarnoldi,Q(1+n*(i-1+1)),px_im )
          call calc_Matvec(pcaarnoldi,px_im,Ax )
#else
          call calc_Matvec( pcaarnoldi,Q(1+n*(i-1)  ),Ax )
          call calc_Matvec( pcaarnoldi,Q(1+n*(i-1+1)),Ax_im )
#endif
          call axpy(n, 0.0d0,-1.0d0,Ax_im,Ax_im)

          ! 実数部の計算
          call axpy(n, ER(i+1),0.0d0,Q(1+n*(i-1)  ),kx)
          call axpy(n, EI(i+1),1.0d0,Q(1+n*(i-1)+1),kx)

          ! 虚数部の計算
          call axpy(n, EI(i+1),0.0d0,Q(1+n*(i-1)  ),kx_im)
          call axpy(n,-ER(i+1),1.0d0,Q(1+n*(i-1)+1),kx_im)

          call calc_dot(n,icomm,kx,kx,kx_norm)
          call calc_dot(n,icomm,kx_im,kx_im,kx_im_norm)
          kx_norm = kx_norm + kx_im_norm
          kx_norm = sqrt(abs(kx_norm))

          call axpy(n,1.0d0,-1.0d0,Ax,kx)
          call axpy(n,1.0d0,-1.0d0,Ax_im,kx_im)

          call calc_dot(n,icomm,kx,kx,err_norm)
          call calc_dot(n,icomm,kx_im,kx_im,kx_im_norm)
          err_norm = err_norm + kx_im_norm
          err_norm = sqrt(abs(err_norm))

          err_norm_all(i+1) = err_norm/kx_norm
#if 0
          if(myrank.eq.0)then
             write(*,'(i8,e18.8,e18.8,e18.8)')i+1,err_norm/kx_norm,ER(i+1),EI(i+1)
          endif
#endif
       endif

    enddo

    maxEigenID = -1
    minEigenID = -1
    do i=1,s
       if(err_norm_all(i) <1.0d0)then
          maxEigenID = i
          minEigenID = i
          maxEigen = EI(i) * EI(i) + ER(i) * ER(i) 
          minEigen = EI(i) * EI(i) + ER(i) * ER(i) 
          exit
       endif
    enddo

    do i=1,s
       if(err_norm_all(i) <1.0d0)then
          tmp = EI(i) * EI(i) + ER(i) * ER(i)
          if(maxEigen<tmp)then
             !          if( ER(i)>0.0d0)then
             maxEigen = tmp
             maxEigenID = i
             !          endif
          endif
          if(minEigen>tmp)then
             !          if( ER(i)<0.0d0)then
             minEigen = tmp
             minEigenID = i
             !          endif
          endif
       endif
    enddo

#if 0
    if(myrank.eq.0)then
       write(*,*)
       write(*,*)"Max Eigen Value ,  ||Ax-kx||/||kx||  "
       if(maxEigenID.eq.-1)then
          write(*,*)"not found eigen value"
       else
          write(*,'(i8,e18.8,e18.8,e18.8)')maxEigenID,ER(maxEigenID),EI(maxEigenID),err_norm_all(maxEigenID)
       endif
       write(*,*)"Min Eigen Value ,  ||Ax-kx||/||kx||  "
       if(minEigenID.eq.-1)then
          write(*,*)"not found eigen value"
       else
          write(*,'(i8,e18.8,e18.8,e18.8)')minEigenID,ER(minEigenID),EI(minEigenID),err_norm_all(minEigenID)
       endif
       write(*,*)
    endif
#endif


    deallocate(px)
    deallocate(Ax)
    deallocate(kx)
    deallocate(px_im)
    deallocate(Ax_im)
    deallocate(kx_im)

  end subroutine checkEigenPair

  ! Q =Qy 
  ! Q : 直交行列(in),固有ベクトル(out)
  ! y : Hessenberg行列の固有ベクトル
  ! ------
  ! yの共役固有値に対する固有ベクトルとしてReal型で実数部、虚数部を独立して持つ.
  ! PARCELは実数のみに対応しているためyのまま計算しても問題ない.
  subroutine getEigenVector(n,s,Q,y,nvecBlock)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    integer nvecBlock
    real*8 Q(n*s)
    real*8 y(s*s)
    real*8 Qy(nvecBlock*s)
    integer i,j,k,jj

    !$OMP PARALLEL DO PRIVATE(jj,j,i,k,Qy) 
    do jj=1,n,nvecBlock 
       do i=1,s
          do j=1,nvecBlock
             Qy(j+(i-1)*nvecBlock)=0.0d0
          enddo
       enddo

       do k=1,s
          do i=1,s
             do j=jj,min(n,jj+nvecBlock-1)
                Qy(j-jj+1+(k-1)*nvecBlock) = &
                     Qy(j-jj+1+(k-1) * nvecBlock) + &
                     Q(j+(i-1)*n)    * y(i+(k-1)*s)
             enddo
          enddo
       enddo

       do i=1,s
          do j=jj,min(n,jj+nvecBlock-1)
             Q(j+(i-1)*n) = Qy(j-jj+1+(i-1)*nvecBlock)
          enddo
       enddo

    enddo
    !$OMP END PARALLEL DO 

  end subroutine getEigenVector

end module KSP_CAARNOLDI


module KSP_CAGMRES_CPU

contains

  subroutine sum_vec_cagmres_CPU( &
       n, &
       s, &
       m, &
       mm, &
       Q , ldq,&
       V , ldv,&
       Qy, ldqy,&
       y, &
       nvecBlock)
    use vec_kernel
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) m
    integer(kind=intType) mm
    integer(kind=intType) s
    integer(kind=intType) ldq
    integer(kind=intType) ldv
    integer(kind=intType) ldqy
    real*8 Q(ldq*(m+1))
    real*8 V(ldv*s)
    real*8 Qy(ldqy)
    real*8 y(m)
    integer nvecBlock

    integer(kind=intType) i,j,jj,k

    !$OMP PARALLEL DO PRIVATE(jj,j,i) 
    do jj=1,n,nvecBlock 
       do j=jj,min(n,jj+nvecBlock-1)
          Qy(j)=0.0d0
       enddo
       do i=1,m-mm
          do j=jj,min(n,jj+nvecBlock-1)
             Qy(j)=Qy(j)+y(i)*Q(j+ldq*(i-1))
          enddo
       enddo
       do i=1,mm
          do j=jj,min(n,jj+nvecBlock-1)
             Qy(j)=Qy(j)+y(i+m-mm)*V(j+ldv*(i-1))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO 

  end subroutine sum_vec_cagmres_CPU

  subroutine calc_NewtonBasis(s,B1,EigenReal,EigenImag)
    use PARCEL_integer
    implicit none
    integer(kind=intType) s
    real*8 B1(s*(s+1))
    real*8 EigenReal(s)
    real*8 EigenImag(s)
    integer(kind=intType) i
    i=1
    do while(i<=s)
       if(i.eq.s)then
          B1(i+(i-1)*(s+1))=EigenReal(i)
       else
          if(abs(EigenImag(i)).eq.0.0d0)then
             B1(i  +(i-1  )*(s+1))=EigenReal(i)
          else
             B1(i  +(i-1  )*(s+1))=EigenReal(i)
             B1(i+1+(i-1+1)*(s+1))=EigenReal(i)
             B1(i  +(i-1+1)*(s+1))=-EigenImag(i)*EigenImag(i)
             i=i+1
          endif
       endif
       i=i+1
    enddo

  end subroutine calc_NewtonBasis

  subroutine getHessenbergMat_Eigen(&
       st,st0,&
       HesMat0,&
       EigenReal,EigenImag,&
       EigenVector &
       )
    use PARCEL_integer
    implicit none
    integer(kind=intType) st,st0
    real*8 HesMat0((st0+1)*st0 )
    real*8,optional :: EigenVector(st*st)
    real*8 HesMat(st*st)
    real*8 EigenReal(st)
    real*8 EigenImag(st)
    real*8 WORK(st*5)
    integer LWORK
    real*8 z(1) ! 'N'を指定しているので参照されない

    integer info
    integer(kind=intType) i,j,k
    integer st_int4
    integer itmp

    real*8 Ax(st)
    real*8 kx(st)

    real*8 err_norm
    real*8 kx_norm
  
    do j=1,st
       do i=1,st
          HesMat(i+(j-1)*st) = HesMat0(i+(j-1)*(st0+1))
!           write(*,*)"Hes     @ ",i,j,HesMat(i+(j-1)*st)
       enddo
!       write(*,*)"------"
    enddo
    st_int4=st

    if(  present(EigenVector)  )then
       EigenReal = 0.0d0
       EigenImag = 0.0d0
       do i=1,st*st
          EigenVector(i) = 0.0d0
       enddo
       LWORK = size(WORK)

       call DGEEV( &
            'N', &
            'V', &
            st_int4, &
            HesMat, &
            st_int4, &
            ! 共役複素数の固有値の場合は虚数部が正の固有値、続いて負の固有値で保存される.
            EigenReal, &
            EigenImag,  &
            ! 左固有ベクトルは計算しないため参照しない
            z, &
            ! 左固有ベクトルは計算しないため参照しない
            st_int4, &
            ! 右固有ベクトル.固有値が実数であれば実数.共役複素数であればa(+-)jbのペアでそれぞれ実数部,虚数部で出力.
            !  -----> 共役複素固有ペアで実数部と虚数部の値を持つためv(:,j) + i*v(:,j+1),v(:,j) - i*v(:,j+1)がそれぞれの固有ベクトルとなる
            EigenVector, &
            ! VRの次元数
            st_int4, &
            ! LWORKの次元を持つreal*8
            WORK, &
            ! st_int4 * 4よりも大きければ良い
            LWORK, &
            info)
       ! Ritz値で誤差評価
       do j=1,st
          do i=1,st
             HesMat(i+(j-1)*st) = HesMat0(i+(j-1)*(st0+1))
          enddo
       enddo
       do k =1,st
          ! 固有値が実数の場合
          if(EigenImag(k).eq.0.0d0)then
             do j=1,st
                Ax(j) = 0.0d0
                do i=1,st
                   Ax(j) =  Ax(j) + HesMat(j+(i-1)*st) * EigenVector(i)
                enddo
             enddo
             kx_norm = 0.0d0
             do j=1,st
                kx(j) = EigenReal(k)*EigenVector(j)
                kx_norm = kx_norm + kx(j) * kx(j) 
             enddo
             kx_norm = sqrt(abs(kx_norm))

             err_norm = 0.0d0
             do j=1,st
                Ax(j) =  Ax(j) -  kx(j)
                err_norm= err_norm +  Ax(j) * Ax(j) 
             enddo
             err_norm = sqrt(abs(err_norm))
!             write(*,*)k,EigenReal(k),err_norm/kx_norm
          endif
       enddo
    else
       itmp = 1
       call dhseqr( &
            'E', 'N', &
            st_int4,  itmp, st_int4, &
            HesMat, st_int4, &
            EigenReal,EigenImag, &
            z, itmp,WORK,st_int4,info)
    endif
#if 0
    write(*,*)
    write(*,*)"----------"
    write(*,*)"DHSEQR"
    do i=1,st
       write(*,*)i,EigenReal(i),EigenImag(i)
    enddo
    write(*,*)
    write(*,*)"----------"
    write(*,*)"DGEEV"
    do i=1,st
       write(*,*)i,EigenReal_DGEEV(i),EigenImag_DGEEV(i)
    enddo
    write(*,*)
    write(*,*)
#endif

  end subroutine getHessenbergMat_Eigen

  ! --------------
  ! PARCELは実行列にのみ対応しているため固有値は実数あるいは複素共役対
  ! 最大,最小は絶対値最大,最小固有値
  subroutine getRealEigen(&
       st,&
       EigenReal,EigenImag,&
       EmaxReal,EmaxImag,&
       EminReal,EminImag &
       )
    use PARCEL_integer
    implicit none
    integer(kind=intType) st
    real*8 EigenReal(st)
    real*8 EigenImag(st)
    real*8 EmaxReal,EmaxImag
    real*8 EminReal,EminImag

    real*8 Eigen
    real*8 Emax,Emin

    real*8 absEigen
    real*8 absEmax,absEmin

    integer i
    if(st<1)then
       EmaxReal=0.0d0
       EmaxImag=0.0d0
       EminReal=0.0d0
       EminImag=0.0d0
       return
    endif
#if 1
    EmaxReal=EigenReal(1)
    EmaxImag=EigenImag(1)
    EminReal=EigenReal(1)
    EminImag=EigenImag(1)

    absEmax=EmaxReal*EmaxReal+EmaxImag*EmaxImag
    absEmin=EminReal*EminReal+EminImag*EminImag
    do i=2,st
       absEigen=EigenReal(i)*EigenReal(i) + EigenImag(i)*EigenImag(i)
       if(absEmax<absEigen)then
          EmaxReal=EigenReal(i)
          EmaxImag=EigenImag(i)
          absEmax=EmaxReal*EmaxReal+EmaxImag*EmaxImag
       endif
       if(absEigen<absEmin)then
          EminReal=EigenReal(i)
          EminImag=EigenImag(i)
          absEmin=EminReal*EminReal+EminImag*EminImag
       endif
    enddo
#else

    EmaxReal=EigenReal(1)
    EmaxImag=EigenImag(1)
    EminReal=EigenReal(1)
    EminImag=EigenImag(1)
    Emax=EigenReal(1)
    Emin=EigenReal(1)
    do i=2,st
       if(Emax<EigenReal(i))then
          EmaxReal=EigenReal(i)
          EmaxImag=EigenImag(i)
          Emax=EmaxReal
       endif
       if(EigenReal(i)<Emin)then
          EminReal=EigenReal(i)
          EminImag=EigenImag(i)
          Emin=EminReal
       endif
    enddo

#endif

  end subroutine getRealEigen

  subroutine calc_hke1eskRRinv(s,k,hk,Rkk,Rinv,Hes)
    use PARCEL_integer
    use wrap_DGEMM,only : wrap_BASIC_DGEMM
    implicit none
    integer(kind=intType) s
    integer(kind=intType) k
    real*8 hk
    real*8 Rkk( (s*(k-1))*s )
    real*8 Rinv(s*s)
    real*8 Hes(s*s)

    real*8 e1est_Rkk(s*s)
    real*8 Hes0(s*s)
    integer(kind=intType) i,j
    ! ---------
    do j=1,s*s
       e1est_Rkk(j)=0.0d0
       Hes0(j)=0.0d0
    enddo

    j=s*(k-1)
    do i=1,s
       !        write(*,*)"  @ Rkk @ ",Rkk( j + (i-1) * ( s*(k-1) ) )
       e1est_Rkk(1 + (i-1)*s )=Rkk( j + (i-1) * ( s*(k-1) ) )*hk
    enddo
    call wrap_BASIC_DGEMM(s,s,s,e1est_Rkk,Rinv,Hes0)    
    do j=1,s*s
       Hes(j)=Hes(j)-Hes0(j)
    enddo

  end subroutine calc_hke1eskRRinv

  subroutine calc_Rk_Bk_Rkinv(s,R0,B,Rinv,RBRinv)
    use PARCEL_integer
    implicit none

    integer(kind=intType) s
    real*8 R0((s+1)*(s+1))
    real*8 B((s+1)*s)
    real*8 Rinv(s*s)
    real*8 RBRinv(s*s)

    integer(kind=intType) i,j
    real*8 RB(s*s)
    real*8 B0(s*s)
    real*8 R(s*s)
    integer s_int4
    do i=1,s
       do j=1,s
          B0(j+(i-1)*s)=B(j+(i-1)*(s+1))
          R(j+(i-1)*s)=R0(j+(i-1)*(s+1))
       enddo
    enddo
    RB=0.0d0

    s_int4=s
    call DGEMM('N','N', s_int4, s_int4, s_int4, 1.0d0, R ,s_int4, B0,  s_int4, 0.0d0,RB    ,s_int4)
    call DGEMM('N','N', s_int4, s_int4, s_int4, 1.0d0, RB,s_int4, Rinv,s_int4, 0.0d0,RBRinv,s_int4)

  end subroutine calc_Rk_Bk_Rkinv

  subroutine calc_HessenbergMatrix( &
       s,t,k, &
       B1, &
       Hes_all, &
       Rtri, &
       BlockR_ubar_prime, &
       hk_pre &
       )
    use PARCEL_integer
    use wrap_DGEMM,only : wrap_BASIC_DGEMM
    implicit none
    integer(kind=intType) s,t,k

    real*8 rho_tilde,rho 
    real*8 hk
    real*8 bk
    real*8 hk_pre
    real*8 zk

    real*8 B1((s+1)*s)
    real*8 Hes_all(  (s*t+1)*s*t )
    real*8 Rtri((s+1)*(s+1))
    real*8 BlockR_ubar_prime( (s*t+1)*(s+1) ) 

    real*8 Rinv((s+1)*(s+1))
    real*8 Rtri0((s+1)*(s+1))
    real*8 BlockR_ubar( (s*t)*(s+1) ) 
    real*8 BlockR( (s*t+1)*s ) 
    real*8 Hes( (s+1)*(s+1) )
    real*8 BlockHes( (s*t+1)*s*t )
    real*8 BlockHes_R(s*(s*t))
    real*8 BlockR_ubar_B1(s*(s*t))
    real*8 R_B_( (s+1)*s )
    integer(kind=intType) i,j
    integer ierr
    integer s_int4
    integer(kind=intType) s1,s2

    Rinv=0.0d0
    Rtri0=0.0d0
    BlockR_ubar=0.0d0
    BlockR=0.0d0
    Hes=0.0d0
    BlockHes=0.0d0
    BlockHes_R=0.0d0
    BlockR_ubar_B1=0.0d0
    R_B_=0.0d0

    if(k.eq.1)then
       do j=1,s
          do i=1,s
             Rinv(i+s*(j-1))=Rtri(i+(s+1)*(j-1))
          enddo
       enddo

       s_int4=s
       call DTRTRI('U', 'N', s_int4, Rinv, s_int4, ierr)

       ! R_B_=R_*B_
       s1=s+1
       call wrap_BASIC_DGEMM(s1,s,s1,Rtri,B1,R_B_)
       call wrap_BASIC_DGEMM(s1,s,s ,R_B_,Rinv,Hes)
       do j=1,s
          do i=1,s+1
             Hes_all(i + (j-1)*(s*t+1) ) = Hes( i + (j-1)*(s+1) )
          enddo
       enddo
       hk_pre=Hes((s+1)*s )
    else
       do i=1,(s+1)*(s+1)
          Rtri0(i)=0.0d0
       enddo
       ! ----
       ! オーバーラップさせた基底の三角行列成分
       do i=1,s+1
          Rtri0(1+(i-1)*(s+1)) = &
               BlockR_ubar_prime( s*(k-1)+1+(i-1)*(s*(k-1)+1) )
       enddo
       ! -----
       ! QR分解で得られたオーバーラップさせていない三角行列成分
       do j=1,s
          do i=1,s
             Rtri0(1+i+(j-1+1)*(s+1))=Rtri(i+(j-1)*s)
          enddo
       enddo

       ! -----
       ! R_(1:s*(k-1)+1,1:s) = R_'(1:s*(k-1)+1,1:s)
       ! R :: [s*(k-1)+1]行,[s]列
       do j=1,s+1
          do i=1,s*(k-1)+1
             BlockR_ubar(i+(j-1)*(s*(k-1)))=BlockR_ubar_prime( i + (j-1)*(s*(k-1)+1) ) 
             !                 write(*,*)"BlockR_ubar ",j,i,BlockR_ubar(i+(j-1)*(s*(k-1)))
          enddo
       enddo
       ! -----
       ! R(1:s*(k-1)   ,1:s) = R_'(1:s*(k-1),1:s)
       ! R :: [s*(k-1)]行,[s]列
       do j=1,s
          do i=1,s*(k-1)
             BlockR(i+(j-1)*(s*(k-1)) )=BlockR_ubar_prime( i + (j-1)*(s*(k-1)+1) ) 
             ! write(*,*)"BlockR ",j,i,BlockR(i+(j-1)*(s*(k-1)) )
          enddo
       enddo

       ! -----
       ! いままで直交化してきた基底と直交化した基底に対する上ヘッセンベルグ行列の要素を計算する
       ! -----
       ! Rinv(1:s,1:s) = Rtri0(1:s,1:s)
       do j=1,s
          do i=1,s
             Rinv(i+s*(j-1))=Rtri0(i+(j-1)*(s+1))
          enddo
       enddo
       s_int4=s
       call DTRTRI('U', 'N', s_int4, Rinv, s_int4, ierr)
       ! -----
       ! 計算済みの上ヘッセンベルグ行列をコピー
       do j=1,s*(k-1)
          do i=1,s*(k-1)
             BlockHes(i+(j-1)*(s*(k-1))) = Hes_all(i + (j-1)*(s*t+1) ) 
          enddo
       enddo
       ! Hes[k-1] * BlockR 
       s1=s*(k-1)
       s2=s+1
       call wrap_BASIC_DGEMM(s1,s,s1,BlockHes          ,BlockR,BlockHes_R)
       ! BlockR_ * B_
       call wrap_BASIC_DGEMM(s1,s,s2,BlockR_ubar,B1    ,BlockR_ubar_B1)

       ! ------
       ! Hesk[k] = ( -Hes[k-1]*BlockR + BlockR_*B_ )* R^-1
       do i=1,s*(k-1) * s
          BlockHes_R(i)=-BlockHes_R(i)+BlockR_ubar_B1(i)
       enddo
       s1=s*(k-1)
       call wrap_BASIC_DGEMM(s1,s,s     ,BlockHes_R,Rinv,BlockR_ubar_B1)
       do i=1,s
          do j=1,s*(k-1)
             Hes_all( j + (i + s*(k-1)-1 ) * (s*t+1)  ) = BlockR_ubar_B1(j+(i-1)*s*(k-1))
          enddo
       enddo

       ! -------
       ! Hoeman 論文(3.22)
       bk=B1((s+1)*s)
       rho_tilde = Rtri0(s   + (s-1)*(s+1))
       rho       = Rtri0(s+1 +     s*(s+1))
       hk=rho/rho_tilde*bk
       Hes_all( s*k+1 +(s*k-1)*(s*t+1) ) =hk

       ! -------
       ! 生成した基底に対するヘッセンベルグ行列の要素
       ! Hoeman 論文(3.21)

       ! Hes = Rk Bk Rinv
       call calc_Rk_Bk_Rkinv(s,Rtri0,B1,Rinv,Hes)
       ! rho~^-1 * bk * zk
       do i=1,s
          zk = Rtri0(i+s*(s+1) )
          Hes(i+s*(s-1)) = Hes(i+s*(s-1))+bk/rho_tilde*zk
       enddo
       ! hk * Rk-1,k Rinv
       call calc_hke1eskRRinv(s,k,hk_pre,BlockR,Rinv,Hes)
       ! -------
       do i=1,s
          do j=1,s
             Hes_all( j+ (i-1)*(s*t+1) + s*(k-1)+ (s*(k-1))*(s*t+1) ) =  Hes(j+s*(i-1) )
          enddo
       enddo
       hk_pre=hk
    endif
  end subroutine calc_HessenbergMatrix

  ! Vk = Vk - Qk-1 R
  subroutine CAGMRES_GEMM(n,s,ms,V,Q,R,nvecBlock)
    use vec_kernel
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    integer(kind=intType) ms
    real*8 V(n*s)
    real*8 Q(n*ms)
    real*8 R(ms*s)
    integer nvecBlock

    integer(kind=intType) ii
    integer(kind=intType) i,j,k

    !$OMP PARALLEL DO PRIVATE(ii,k,j,i) 
    do ii=1,n,nvecBlock
       do k=1,s
          do j=1,ms
             do i=ii,min(ii+nvecBlock-1,n)
                V(i+n*(k-1)) = &
                     V(i+n*(k-1)) - &
                     Q(i+n*(j-1))*R(j+ms*(k-1))
             enddo
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO 
    
  end subroutine CAGMRES_GEMM

  ! V(i)=V(i)+alpha*V(i-1)+beta*V(i-2)
  subroutine CAGMRES_axpy3(n,alpha,beta,V,V1,V2)
    use PARCEL_integer
    implicit none 
    integer(kind=intType) n
    real*8 alpha,beta
    real*8 V(n)
    real*8 V1(n)
    real*8 V2(n)
    integer(kind=intType) i

    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
       V(i)=V(i)+alpha*V1(i)+beta*V2(i)
    enddo

  end subroutine CAGMRES_axpy3

  subroutine calc_NewtonBasisVec_CPU(pcagmres,n,s,i,pv,V,B1)
    use krylov
    use krylov_kernel
    use vec_kernel
    use PARCEL_integer
    implicit none
    type(KSP) :: pcagmres
    integer(kind=intType) n
    integer(kind=intType) s
    integer(kind=intType) i
    real*8 pV(n)
    real*8 V(n*(s+1))
    real*8 B1((s+1)*s)
    real*8 alpha,beta
    integer i_int4

    i_int4=1
    call calc_Matvec( pcagmres,pv,V(1+i*n) )
    if(i.ne.i_int4)then
       alpha = -B1( i + (i-1)*(s+1) )
       beta  = B1( i-1 + (i-1-1+1)*(s+1) )
       call CAGMRES_axpy3(n,alpha,beta,V(1+i*n),V(1+(i-1)*n),V(1+(i-1-1)*n))
    else
       alpha = -B1( i + (i-1)*(s+1) )
       call axpy(n,alpha,1.0d0,V(1+(i-1)*n),V(1+i*n))
    endif
  end subroutine calc_NewtonBasisVec_CPU

end module KSP_CAGMRES_CPU

module KSP_CAGMRES

contains
  subroutine update_xVec(&
       pgmres,&
       n,m0,m,vecLen,&
       x,v,x0,v0,&
       hki,e,&
       y_,nvecBlock)
    use krylov,only : KSP
    use PARCEL_integer,only : intType
    use PARCEL_KRYLOV_PRECON,only : calc_precon
    use PARCEL_VEC_kernel,only : axpy
    use KSP_GMRES_CPU,only : sum_vec_FP64_CPU
    use vec_type,only : PARCEL_VEC
    use KSP_STATUS,only : &
       STATUS_FLAG, &
       init_staus_flag,set_machineflag, &
       set_dataflag,set_calcflag
    use PARCEL_define,only : PARCEL_CPU_FLAG_ON,PARCEL_GPU_FLAG_ON
    use c_wrapper,only :  CUDA_NUM_BLOCKS , CUDA_NUM_THREADS
    use cu_common_wrapper,only : cudaMemcpy_h2d_FP64
    use cu_Vec,only : cuda_sum_vec_fp64_gpu
    use PARCEL_VEC_kernel,only : init_vec,free_vec,set_val_FP64

    implicit none

    type(KSP) :: pgmres
    integer(kind=intType) n
    integer(kind=intType) m
    integer(kind=intType) m0
    integer(kind=intType) vecLen
    type(PARCEL_VEC) :: x
    type(PARCEL_VEC) :: v
    type(PARCEL_VEC) :: x0
    type(PARCEL_VEC) :: v0
    type(PARCEL_VEC) :: y_

    real*8 y(m0+1)
    real*8 hki( (m0+1)*m0 )
    real*8 e(m0+1)
    integer nvecBlock

    type(PARCEL_VEC) :: tmp1
    type(PARCEL_VEC) :: tmp2

    real*8 Hinv( m*m )
!    real*8 tmp1,tmp2
    integer(kind=intType) i,j
    integer ierr
    type(STATUS_FLAG) :: flag
    integer(kind=intType),parameter :: offset_zero=0
    integer(kind=intType) parcel_int
    
    flag = pgmres%flag
    ! y = H^-1 e 
    do i=m,1,-1
       y(i)=e(i)
       do j=m,i+1,-1
          y(i)=y(i)-hki(i+(j-1)*(m0+1))*y(j)
       enddo
       if(hki(i+(i-1)*(m0+1)).eq.0.0d0)then
          y(i)=y(i)/1.0d-300
       else
          y(i)=y(i)/hki(i+(i-1)*(m0+1))
       endif
    enddo

    select case(flag%machineflag)
    case(PARCEL_CPU_FLAG_ON)
       call sum_vec_FP64_CPU(&
            n,m, &
            v%vec_FP64 , vecLen,&
            v0%vec_FP64, vecLen,&
            y, m0+1,&
            nvecBlock)
    case(PARCEL_GPU_FLAG_ON)
       call cudaMemcpy_h2d_FP64(&
          y,y_%C_vec_FP64,m0+1,&
          offset_zero, &
          offset_zero) 
       call cuda_sum_vec_FP64_GPU( &
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS,&
            n, &
            m, &
            vecLen,&
            v%C_vec_FP64 , v%offset, &
            v0%C_vec_FP64, v0%offset, &
            y_%C_vec_FP64 &
            )

    end select
    ! solve
    call calc_precon(pgmres,v0,x0 )
    parcel_int = 1
    call init_vec(pgmres,tmp1,parcel_int)
    call init_vec(pgmres,tmp2,parcel_int)
    call set_val_FP64(pgmres,tmp1,1.0d0,parcel_int)
    call set_val_FP64(pgmres,tmp2,1.0d0,parcel_int)
    call axpy(pgmres,tmp1,tmp2,x0,x,n)
    call free_vec(pgmres,tmp1)
    call free_vec(pgmres,tmp2)

  end subroutine update_xVec

  subroutine sum_vec_cagmres( &
       pcagmres,&
       n,s,m,mm,&
       Q, ldq,&
       V, ldv,&
       Qy, ldqy,&
       y,  ldy, &
       nvecBlock)

    use krylov,only : KSP
    use vec_type,only : PARCEL_VEC
    use PARCEL_VEC_kernel,only : get_vec_FP64,set_vec_FP64
    use PARCEL_integer,only : intType
    use wrap_DGEMM,only : wrap_BASIC_DGEMM
    use KSP_CAGMRES_CPU,only : sum_vec_cagmres_CPU
    use KSP_STATUS,only : &
       STATUS_FLAG,PARCEL_CPU_FLAG_ON,PARCEL_GPU_FLAG_ON
    use cu_Vec,only : cuda_sum_vec_cagmres_FP64_GPU
    use c_wrapper,only :  CUDA_NUM_BLOCKS , CUDA_NUM_THREADS
    implicit none

    type(KSP) :: pcagmres
    integer(kind=intType) n
    integer(kind=intType) s
    integer(kind=intType) m
    integer(kind=intType) mm

    type(PARCEL_VEC) :: Q
    integer(kind=intType) ldq
    type(PARCEL_VEC) :: V
    integer(kind=intType) ldv
    type(PARCEL_VEC) :: Qy
    integer(kind=intType) ldqy
    type(PARCEL_VEC) :: y
    integer(kind=intType) ldy

    integer nvecBlock

    integer(kind=intType) i,j,jj,k


    type(STATUS_FLAG) :: flag
    flag = pcagmres%flag

#if 0

    call get_vec_FP64(pcagmres,Q,Q_,size(Q_))
    call get_vec_FP64(pcagmres,V,V_,size(V_))
    call get_vec_FP64(pcagmres,Qy,Qy_,size(Qy_))
    call get_vec_FP64(pcagmres,y,y_,size(y_))

       call sum_vec_cagmres_CPU( &
            n, &
            s, &
            m, &
            mm, &
            Q_ ,&
            V_ ,&
            Qy_,&
            y_,&
            nvecBlock)

    call set_vec_FP64(pcagmres,Q,Q_,size(Q_))
    call set_vec_FP64(pcagmres,V,V_,size(V_))
    call set_vec_FP64(pcagmres,Qy,Qy_,size(Qy_))
    call set_vec_FP64(pcagmres,y,y_,size(y_))

    return 

#endif

    select case(flag%machineflag)
       
    case(PARCEL_CPU_FLAG_ON)
       call sum_vec_cagmres_CPU( &
            n, &
            s, &
            m, &
            mm, &
            Q%vec_FP64(1+Q%offset)   , ldq,&
            V%vec_FP64(1+V%offset)   , ldv,&
            Qy%vec_FP64(1+Qy%offset) , ldqy,&
            y%vec_FP64(1+y%offset)   ,&
            nvecBlock)
    case(PARCEL_GPU_FLAG_ON)
       call cuda_sum_vec_cagmres_FP64_GPU( &
            CUDA_NUM_BLOCKS,CUDA_NUM_THREADS, &
            n, &
            s, &
            m, &
            mm, &
            Q%C_vec_FP64 ,  Q%offset , ldq,&
            V%C_vec_FP64 ,  V%offset , ldv,&
            Qy%C_vec_FP64,  Qy%offset, ldqy,&
            y%C_vec_FP64 ,  y%offset )

    end select

  end subroutine sum_vec_cagmres

  subroutine update_xVec_cagmres(&
       pcagmres,&
       n,m0,m,s,&
       x   , ldx,&
       Q   , ldq,&
       V   , ldv,&
       x0  , ldx0,&
       v0  , ldv0,&
       hki , ldhki,&
       e   , lde,&
       Rtri, ldRtri,&
       work, &
       nvecBlock)
    !    use QR_kernel
    use krylov,only : KSP
    use vec_type,only : PARCEL_VEC
    use PARCEL_integer,only : intType
    use wrap_DGEMM,only : wrap_BASIC_DGEMM
    use PARCEL_VEC_kernel,only : &
         init_vec,free_vec, &
         axpy,set_vec_FP64,set_val_FP64
    use PARCEL_KRYLOV_PRECON,only : calc_precon

    implicit none

    type(KSP) :: pcagmres
    integer(kind=intType) n
    integer(kind=intType) m
    integer(kind=intType) m0
    integer(kind=intType) s

    ! real*8 x(n)
    ! real*8 Q(n*(m0+1))
    ! real*8 x0(n)
    ! real*8 v0(n)
    ! real*8 V(n*(s))
    ! real*8 Rtri(s*s)

    integer(kind=intType) ldx
    integer(kind=intType) ldq
    integer(kind=intType) ldx0
    integer(kind=intType) ldv0
    integer(kind=intType) ldv
    integer(kind=intType) ldrtri
    integer(kind=intType) ldhki
    integer(kind=intType) lde

    type(PARCEL_VEC) :: x
    type(PARCEL_VEC) :: Q
    type(PARCEL_VEC) :: x0
    type(PARCEL_VEC) :: v0
    type(PARCEL_VEC) :: V

    real*8  Rtri(ldrtri*s)
!    real*8 hki( (ldhki)*m0 )
    real*8 hki( (m0+1)*m0 )
    real*8 e(lde)

    type(PARCEL_VEC) :: work

    integer nvecBlock

    real*8 Rinv0(s*s)
    real*8 Rinv((s-1)*(s-1))
    real*8 y_tmp(s)

    real*8 y(m0+1)

    real*8 Hinv( m*m )
    type(PARCEL_VEC) :: tmp1
    type(PARCEL_VEC) :: tmp2
    integer(kind=intType) i,j
    integer ierr
    integer(kind=intType) Vstep
    integer s_int4
    integer(kind=intType) parcel_int

    Vstep=(s-1)-(m0-m)
    !    write(*,*)"m,Vstep",m,Vstep

!    write(*,*)"m0 = ",m0," : m = ",m
    ! y = H^-1 e 
    do i=m,1,-1
       y(i)=e(i)
       do j=m,i+1,-1
          y(i)=y(i)-hki(i+(j-1)*(m0+1))*y(j)
       enddo
       y(i)=y(i)/hki(i+(i-1)*(m0+1))
    enddo

    do i=m-Vstep+1,m
       y_tmp(1+i-(m-Vstep+1))=y(i)
    enddo
    do j=1,s
       do i=1,s
          Rinv0(i+(j-1)*(s))=Rtri(i+(j-1)*s)
       enddo
    enddo
    s_int4=s
    call DTRTRI('U', 'N', s_int4, Rinv0, s_int4, ierr)  
    Rinv=0.0d0
    do j=1,s-1
       do i=1,s-1
          Rinv(i+(j-1)*(s-1))=Rinv0(i+(j-1)*s)
       enddo
    enddo

    i=s-1
    j=1
    call wrap_BASIC_DGEMM(i,j,i,Rinv,y_tmp,y(m-Vstep+1))

!    parcel_int = size(y)
    parcel_int = m0+1
    call set_vec_FP64(pcagmres,work,y,parcel_int)

    parcel_int = 1
    call sum_vec_cagmres( &
         pcagmres,&
         n,s,m,Vstep,&
         Q , ldq , &
         V , ldv , &
         v0, ldv0, &
         work , parcel_int, &
         nvecBlock)

    ! solve Mx0=v0
    call calc_precon(pcagmres,v0,x0)

    parcel_int = 1
    call init_vec(pcagmres,tmp1,parcel_int)
    call init_vec(pcagmres,tmp2,parcel_int)
    ! tmp1=1.0d0
    ! tmp2=1.0d0
    call set_val_FP64(pcagmres, tmp1, 1.0d0,parcel_int)
    call set_val_FP64(pcagmres, tmp2, 1.0d0,parcel_int)
    ! x = x+x0
    call axpy(pcagmres,tmp1,tmp2,x0,x,n)
    call free_vec(pcagmres,tmp1)
    call free_vec(pcagmres,tmp2)

  end subroutine update_xVec_cagmres 

  subroutine calc_NewtonBasisVec( &
       pcagmres,n,s,i,vecLen,pv,V,B1)
    use krylov,only : KSP
    use PARCEL_integer,only : intType
    use PARCEL_VEC_kernel,only : axpy,axpy_vec3
    use vec_type,only : PARCEL_VEC
    use PARCEL_KRYLOV_MATRIX,only : calc_Matvec
    use PARCEL_VEC_kernel,only : set_val_FP64,init_vec,free_vec
    implicit none
    type(KSP) :: pcagmres
    integer(kind=intType) n
    integer(kind=intType) s
    integer(kind=intType) i
    integer(kind=intType) vecLen

    type(PARCEL_VEC) :: pV
    type(PARCEL_VEC) :: V

    integer(kind=intType) offset0_V
    type(PARCEL_VEC) :: alpha
    type(PARCEL_VEC) :: beta
    type(PARCEL_VEC) :: gamma

    real*8 B1((s+1)*s)
    integer i_int4
    integer(kind=intType) parcel_int

    
    parcel_int = 1
    call init_vec(pcagmres,alpha,parcel_int)
    call init_vec(pcagmres,beta ,parcel_int)
    call init_vec(pcagmres,gamma,parcel_int)
    offset0_V = V%offset

    i_int4=1
    V%offset = offset0_V + i*vecLen
!    call calc_Matvec( pcagmres,pv,V(1+i*n) )
    call calc_Matvec( pcagmres,pv,V )
    V%offset = offset0_V 
    if(i.ne.i_int4)then
       parcel_int = 1
       call set_val_FP64(pcagmres,alpha,-B1( i + (i-1)*(s+1) ),parcel_int)
       call set_val_FP64(pcagmres,beta , B1( i-1 + (i-1-1+1)*(s+1) ) ,parcel_int)
       call set_val_FP64(pcagmres,gamma, 1.0d0 ,parcel_int)
       ! alpha = -B1( i + (i-1)*(s+1) )
       ! beta  = B1( i-1 + (i-1-1+1)*(s+1) )
       ! gamma = 1.0d0

       call axpy_vec3(&
            pcagmres,&
            alpha,beta,gamma,&
            V,V,V, &
            n,&
            V%offset+(i-1)*vecLen,  V%offset+(i-1-1)*vecLen,  V%offset+i*vecLen)
    else
       parcel_int = 1
       call set_val_FP64(pcagmres,alpha, -B1( i + (i-1)*(s+1) ) ,parcel_int)
       call set_val_FP64(pcagmres,beta , 1.0d0                  ,parcel_int)
       ! call axpy(n,alpha,1.0d0,V(1+(i-1)*n),V(1+i*n))
       call axpy( &
            pcagmres,alpha,beta, &
            V,V,&
            n, &
            V%offset+(i-1)*vecLen,V%offset+i*vecLen )

    endif
    V%offset = offset0_V 

    call free_vec(pcagmres,alpha)
    call free_vec(pcagmres,beta )
    call free_vec(pcagmres,gamma)

  end subroutine calc_NewtonBasisVec

end module KSP_CAGMRES

! ------
subroutine KSP_pCAGMRES(pcagmres)

  use krylov,only : KSP
  use PARCEL_define,only : PARCEL_CPU_FLAG_ON,PARCEL_GPU_FLAG_ON
  use KSP_STATUS,only : &
       STATUS_FLAG, &
       init_staus_flag,set_machineflag, &
       set_dataflag,set_calcflag

  use PARCEL_KRYLOV_MATRIX,only : calc_Matvec,calc_res
  use PARCEL_KRYLOV_PRECON,only : calc_precon
  use PARCEL_VEC_kernel
  use vec_type,only : PARCEL_VEC
  use PARCEL_integer,only : intType
  use KSP_CAGMRES,only : update_xVec,update_xVec_cagmres,calc_NewtonBasisVec
  use PARCEL_define,only : &
       PARCEL_CPU_FLAG_ON,PARCEL_GPU_FLAG_ON, &
       KSP_QRDECOMPOSITION_MGS, &
       KSP_QRDECOMPOSITION_CGS 
  use QR_kernel,only : &
       QRDecomposition_, &
       GS_QR_,GS_QR, &
       CholeskyQR_Ronly_,calcQ_VRinv_

  use KSP_CAGMRES_CPU,only : &
       calc_HessenbergMatrix, calc_NewtonBasis,getHessenbergMat_Eigen, &
       CAGMRES_GEMM

  use KSP_GMRES_CPU,only : GivensRotations
  use cu_common_wrapper,only : &
       cublas_init_DgemmBatched_gemm_dot_local_FP64, &
       cublas_free_DgemmBatched_gemm_dot_local_FP64

  use PARCEL_define,only : WARP_SIZE
  use c_wrapper,only : CUDA_NUM_BLOCKS

  implicit none
  type(KSP) :: pcagmres
  integer(kind=intType) n
  integer icomm
  integer itrmax ! 最大反復回数
  real*8 rtolmax ! 相対残差 (収束条件)
  real*8 abstolmax  !  絶対残差 (収束条件)
  integer myrank ! プロセスランク
  integer(kind=intType) s ! 省通信ステップ 
  integer(kind=intType) t ! リスタート数 = (s*t)
  integer iQRflag ! QR分解の方法 
  integer basis ! 基底ベクトルの生成方法
  integer iQmatflag
  type(STATUS_FLAG) :: flag

  ! ----------
  integer itr,itr0
  integer(kind=intType) mm
  integer(kind=intType) k
  integer(kind=intType) i,j
  integer(kind=intType) msize
  integer nvecBlock
  ! ----------

  type(PARCEL_VEC) :: x
  type(PARCEL_VEC) :: r
  type(PARCEL_VEC) :: b
  type(PARCEL_VEC) :: x0
  type(PARCEL_VEC) :: v0
  type(PARCEL_VEC) :: pV
  type(PARCEL_VEC) :: V
  type(PARCEL_VEC) :: Q
  type(PARCEL_VEC) :: zeta
  type(PARCEL_VEC) :: y
  type(PARCEL_VEC) :: B1
  type(PARCEL_VEC) :: EigenReal
  type(PARCEL_VEC) :: EigenImag
  type(PARCEL_VEC) :: HessenbergMat
  type(PARCEL_VEC) :: HessenbergMat_Givens
  type(PARCEL_VEC) :: Rtri
  type(PARCEL_VEC) :: BlockR_ubar_prime
  type(PARCEL_VEC) :: co
  type(PARCEL_VEC) :: si
  type(PARCEL_VEC) :: bnorm,rnorm
  type(PARCEL_VEC) :: tmp1,tmp2

  ! --------------------------------
  ! Batch GEMM
  ! V x V 
  integer(kind=intType) batchCount
  integer(kind=intType) localnRow
  integer(kind=intType) batch_len
  type(devicemem) V_s_batch_prt  
  type(devicemem) V_s1_batch_prt 
  type(devicemem) work_s_batch_prt 
  type(devicemem) work_s1_batch_prt 
  ! --------------------------------
  type(PARCEL_VEC) :: work ! Batch GEMM,TSQRで使用

  real*8 bnorm_FP64,rnorm_FP64
  real*8 tmp1_FP64,tmp2_FP64
  real*8,allocatable :: zeta_FP64(:)
  real*8,allocatable :: B1_FP64(:)
  real*8,allocatable :: EigenReal_FP64(:)
  real*8,allocatable :: EigenImag_FP64(:)
  real*8,allocatable :: HessenbergMat_FP64(:)
  real*8,allocatable :: HessenbergMat_Givens_FP64(:)
  real*8,allocatable :: Rtri_FP64(:)
  real*8,allocatable :: BlockR_ubar_prime_FP64(:)
  real*8,allocatable :: co_FP64(:)
  real*8,allocatable :: si_FP64(:)

  real*8 hk_pre_FP64
  real*8 alpha_FP64,beta_FP64

  integer(kind=intType) vecLen
  integer(kind=intType) nColL
  integer(kind=intType) nColU

  integer iflag_exit
  integer iflag_basis
  integer ierr
  real*8 rtime(40)
  real*8 stime0,etime0
  real*8 stime,etime

  integer(kind=intType) s1,s2,s3
  integer s1_int4
  integer(kind=intType) parcel_int
  integer(kind=intType),parameter :: offset_zero=0

  itrmax    = pcagmres%itrmax
  rtolmax   = pcagmres%rtolmax
  abstolmax = pcagmres%abstolmax
  myrank = pcagmres%commVec%myrank
  icomm  = pcagmres%commVec%comm
  n      = pcagmres%prm%n

  s = pcagmres%CAGMRES_sstep
  t = pcagmres%CAGMRES_tstep
  iQRflag = pcagmres%CAGMRES_QRflag
  basis = pcagmres%CAGMRES_basis
  nvecBlock= pcagmres%nvecBlock


  call MPI_Barrier(icomm,ierr)

  stime0 = MPI_Wtime()
  rtime  = 0.0d0

#if 0

  call init_staus_flag(flag)
#if 1
  call set_machineflag(flag,PARCEL_CPU_FLAG_ON)
#else
  call set_machineflag(flag,PARCEL_GPU_FLAG_ON)
#endif
  ! ここではcagmres内のflagから設定しているが将来はKSP構造体で設定したものを受け取る
  call set_dataflag(flag,KSP_DOUBLE_PRECISION)
  call set_calcflag(flag,KSP_DOUBLE_PRECISION)
  pcagmres%flag%machineflag = flag%machineflag
  pcagmres%flag%dataflag = flag%dataflag 
  pcagmres%flag%calcflag = flag%calcflag 

#endif

  call parcel_sync(pcagmres)

  nColL = 0
  nColU = 0
  if(pcagmres%prm%MatFormat.eq.MATRIX_FORMAT_CRS)then
     nColL = pcagmres%crsMat%nColL
     nColU = pcagmres%crsMat%nColU
  endif
  vecLen = n+nColL+nColU

  if(pcagmres%flag%machineflag.eq.PARCEL_CPU_FLAG_ON)then
     vecLen = n
  endif

  ! --------------------------------
  ! CUBLAS Batch paameter set
  ! CRS,DIA -- DDM で小行列のサイズを変更する
  batchCount = 1
  if(pcagmres%flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
     !     localnRow = 32
     localnRow = 64 ! s =4
     !     localnRow = 128
     batch_len = n
     if( mod(n,localnRow).ne.0 )then
        batch_len = batch_len + ( localnRow - mod(n,localnRow) )
     endif
     batchCount = batch_len / localnRow
     vecLen = max(vecLen,batch_len)
  endif
  ! --------------------------------

  call init_vec(pcagmres,x,vecLen)
  call init_vec(pcagmres,r,vecLen)
  call init_vec(pcagmres,b,vecLen)
  call init_vec(pcagmres,x0,vecLen)
  call init_vec(pcagmres,v0,vecLen)
  call init_vec(pcagmres,pV,vecLen)
  call init_vec(pcagmres,V,vecLen*(s+1))
  call init_vec(pcagmres,Q,vecLen*(s*t+1))
  !   call init_vec(pcagmres,work,batchCount*(s+1)*(s+1))
  call init_vec(pcagmres,work, &
       max( &
       batchCount*(s+1)*(s+1), &
       ((s+1)*(s+1)+(s+1)*((s+1)+WARP_SIZE))*(CUDA_NUM_BLOCKS) &
       ) &
       )

  call init_vec(pcagmres,B1,(s+1)*s)
  call init_vec(pcagmres,EigenReal,s*t)
  call init_vec(pcagmres,EigenImag,s*t)
  call init_vec(pcagmres,HessenbergMat,(s*t+1)*s*t )
  call init_vec(pcagmres,HessenbergMat_Givens,(s*t+1)*s*t)
  call init_vec(pcagmres,Rtri,(s+1)*(s+1))
  call init_vec(pcagmres,BlockR_ubar_prime,(s*t+1)*(s+1))
  call init_vec(pcagmres,zeta,(s*t+1))
  call init_vec(pcagmres,y,(s*t+1))
  call init_vec(pcagmres,co,(s*t+1))
  call init_vec(pcagmres,si,(s*t+1))

  parcel_int = 1
  call init_vec(pcagmres,bnorm,parcel_int)
  call init_vec(pcagmres,rnorm,parcel_int)
  call init_vec(pcagmres,tmp1 ,parcel_int)
  call init_vec(pcagmres,tmp2 ,parcel_int)

  call set_val_FP64(pcagmres,tmp1,0.0d0,parcel_int)
  call set_val_FP64(pcagmres,bnorm,0.0d0,parcel_int)
  call set_val_FP64(pcagmres,rnorm,0.0d0,parcel_int)
  call set_val_FP64(pcagmres,tmp1 ,0.0d0,parcel_int)
  call set_val_FP64(pcagmres,tmp2 ,0.0d0,parcel_int)
 
  call set_vec_val(pcagmres, x                   , tmp1, parcel_int, vecLen)
  call set_vec_val(pcagmres, r                   , tmp1, parcel_int, vecLen)
  call set_vec_val(pcagmres, b                   , tmp1, parcel_int, vecLen)
  call set_vec_val(pcagmres, x0                  , tmp1, parcel_int, vecLen)
  call set_vec_val(pcagmres, v0                  , tmp1, parcel_int, vecLen)
  call set_vec_val(pcagmres, pV                  , tmp1, parcel_int, vecLen)
  call set_vec_val(pcagmres, V                   , tmp1, parcel_int, vecLen*(s+1))
  call set_vec_val(pcagmres, Q                   , tmp1, parcel_int, vecLen*(s*t+1))
  call set_vec_val(pcagmres, work                , tmp1, parcel_int, batchCount*(s+1)*(s+1))
  call set_vec_val(pcagmres, B1                  , tmp1, parcel_int, (s+1)*s)
  call set_vec_val(pcagmres, EigenReal           , tmp1, parcel_int, s*t)
  call set_vec_val(pcagmres, EigenImag           , tmp1, parcel_int, s*t)
  call set_vec_val(pcagmres, HessenbergMat       , tmp1, parcel_int, (s*t+1)*s*t)
  call set_vec_val(pcagmres, HessenbergMat_Givens, tmp1, parcel_int, (s*t+1)*s*t)
  call set_vec_val(pcagmres, Rtri                , tmp1, parcel_int, (s+1)*(s+1))
  call set_vec_val(pcagmres, BlockR_ubar_prime   , tmp1, parcel_int, (s*t+1)*(s+1))
  call set_vec_val(pcagmres, zeta                , tmp1, parcel_int, (s*t+1))
  call set_vec_val(pcagmres, y                   , tmp1, parcel_int, (s*t+1))
  call set_vec_val(pcagmres, co                  , tmp1, parcel_int, (s*t+1))
  call set_vec_val(pcagmres, si                  , tmp1, parcel_int, (s*t+1))

  allocate( B1_FP64((s+1)*s) )
  allocate( EigenReal_FP64(s*t) )
  allocate( EigenImag_FP64(s*t) )
  allocate( HessenbergMat_FP64((s*t+1)*s*t) )
  allocate( HessenbergMat_Givens_FP64((s*t+1)*s*t) )
  allocate( Rtri_FP64((s+1)*(s+1)) )
  allocate( BlockR_ubar_prime_FP64((s*t+1)*(s+1)) )
  allocate( zeta_FP64((s*t+1)) )
  allocate( co_FP64((s*t+1)) ) 
  allocate( si_FP64((s*t+1)) ) 

  if(pcagmres%flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
     parcel_int = batchCount
     call cublas_init_DgemmBatched_gemm_dot_local_FP64(&
          parcel_int, &
          V%C_vec_FP64, V%offset, localnRow,&
          V_s1_batch_prt&
          ) 

     V%offset = vecLen
     parcel_int = batchCount
     call cublas_init_DgemmBatched_gemm_dot_local_FP64(&
          parcel_int, &
          V%C_vec_FP64, V%offset, localnRow,&
          V_s_batch_prt&
          ) 
     V%offset = 0

     ! work size : batchCount*s*s 
     parcel_int = batchCount
     call cublas_init_DgemmBatched_gemm_dot_local_FP64(&
          parcel_int, &
          work%C_vec_FP64, work%offset, s,  &
          work_s_batch_prt &
          ) 
     ! work size : batchCount*(s+1)*(s+1) 
     parcel_int = batchCount
     call cublas_init_DgemmBatched_gemm_dot_local_FP64(&
          parcel_int, &
          work%C_vec_FP64, work%offset, s+1,  &
          work_s1_batch_prt &
          ) 
  endif

  B1_FP64 = 0.0d0
  EigenReal_FP64 = 0.0d0
  EigenImag_FP64 = 0.0d0
  HessenbergMat_FP64 = 0.0d0
  HessenbergMat_Givens_FP64 = 0.0d0
  Rtri_FP64 = 0.0d0
  BlockR_ubar_prime_FP64 = 0.0d0
  zeta_FP64 = 0.0d0
  co_FP64 = 0.0d0
  si_FP64 = 0.0d0

  ! --------------------------------
  call parcel_sync(pcagmres)
#ifdef CUDA_SOLVER
  call cudaMemcpy_d2d_FP64(pcagmres%C_x,x%C_vec_FP64,n,offset_zero,offset_zero)
  call cudaMemcpy_d2d_FP64(pcagmres%C_b,b%C_vec_FP64,n,offset_zero,offset_zero)
#else
  call set_vec_FP64(pcagmres,x,pcagmres%x,n)
  call set_vec_FP64(pcagmres,b,pcagmres%b,n)
#endif
  call parcel_sync(pcagmres)

  ! r = b-Ax
  call calc_res(pcagmres,x,b,r)
  call calc_dot(pcagmres,icomm,b,b,bnorm,n)
  call calc_dot(pcagmres,icomm,r,r,rnorm,n)
  parcel_int = 1
  call sqrt_val(pcagmres,bnorm,parcel_int)
  call sqrt_val(pcagmres,rnorm,parcel_int)
  call get_val_FP64(pcagmres,bnorm,bnorm_FP64,parcel_int)
  call get_val_FP64(pcagmres,rnorm,rnorm_FP64,parcel_int)
  call parcel_sync(pcagmres)

  ! ----------
  ! CACAGMRES のメインループ
  ! do itr = 1,itrmax,s*t
  !   r = b-Ax
  !   do k = 1,t
  !     基底ベクトルの生成
  !     do i=1,s
  !       solve Ms=V_i
  !       Newton基底,単基底
  !       V_i+1=F(As)
  !     enddo
  !     if(k.eq.1)then
  !       QR分解
  !     else
  !       ブロック古典グラムシュミット直交化
  !       QR分解
  !     endif
  !     ヘッセンベルグ行列の計算
  !     ギブンス回転
  !     if(||r||<rmax) ;exit
  !   enddo
  !   解ベクトル更新
  !   if(||r||<rmax) ;exit
  ! enddo
  ! ----------

  B1_FP64   = 0.0d0
  Rtri_FP64 = 0.0d0
  zeta_FP64 = 0.0d0

  ! ----
  ! 基底変換行列
  do i=1,s
     B1_FP64(1+i+(i-1)*(s+1)) = 1.0d0
  enddo
  parcel_int = size(B1_FP64)
  call set_vec_FP64(pcagmres,B1,B1_FP64,parcel_int)

  iflag_exit  = -1
  iflag_basis = -1

  EigenReal_FP64 = 0.0d0
  EigenImag_FP64 = 0.0d0
  alpha_FP64     = 0.0d0
  beta_FP64      = 0.0d0

  parcel_int = size(EigenReal_FP64)
  call set_vec_FP64(&
       pcagmres, EigenReal, EigenReal_FP64, parcel_int)
  parcel_int = size(EigenImag_FP64)
  call set_vec_FP64(&
       pcagmres, EigenImag, EigenImag_FP64, parcel_int)

  do itr = 1,itrmax,s*t
     itr0 = itr
     co_FP64   = 0.0d0
     si_FP64   = 0.0d0
     zeta_FP64 = 0.0d0

     HessenbergMat_FP64 = 0.0d0
     HessenbergMat_Givens_FP64 = 0.0d0

     Rtri_FP64 = 0.0d0
     BlockR_ubar_prime_FP64 = 0.0d0
     hk_pre_FP64            = 0.0d0

  parcel_int = size(co_FP64)
     call set_vec_FP64(pcagmres,co,co_FP64,parcel_int)
  parcel_int = size(si_FP64)
     call set_vec_FP64(pcagmres,si,si_FP64,parcel_int)
  parcel_int = size(zeta_FP64)
     call set_vec_FP64(pcagmres,zeta,zeta_FP64,parcel_int)

  parcel_int = size(HessenbergMat_FP64)
     call set_vec_FP64(pcagmres,HessenbergMat,HessenbergMat_FP64,parcel_int)

  parcel_int = size(HessenbergMat_Givens_FP64)
     call set_vec_FP64(pcagmres,HessenbergMat_Givens,HessenbergMat_Givens_FP64,parcel_int)

  parcel_int = size(Rtri_FP64)
     call set_vec_FP64(pcagmres,Rtri,Rtri_FP64,parcel_int)
  parcel_int = size(BlockR_ubar_prime_FP64)
     call set_vec_FP64(pcagmres,BlockR_ubar_prime,BlockR_ubar_prime_FP64,parcel_int)
!     call set_vec_FP64(pcagmres,hk_pre,hk_pre_FP64,size(hk_pre_FP64))

     ! r=b-Ax
     ! call MPI_Barrier(icomm,ierr)
     stime=MPI_Wtime()
     call calc_res(pcagmres,x,b,r)
     call calc_dot(pcagmres,icomm,r,r,rnorm,n)
     parcel_int = 1
     call get_val_FP64(pcagmres,rnorm,rnorm_FP64,parcel_int)
     rnorm_FP64 = sqrt(rnorm_FP64)

     if( (rnorm_FP64<abstolmax) .or. (rnorm_FP64/bnorm_FP64<rtolmax) ) then
        exit
     endif

     if(itr.eq.1)then
        pcagmres%reshistory(itr) = rnorm_FP64 / bnorm_FP64 
     endif

     if(myrank.eq.0)then
!        write(*,*)"@ res @",itr,rnorm_FP64,rnorm_FP64/bnorm_FP64
     endif

     zeta_FP64(1) = -rnorm_FP64
     tmp1_FP64 = -1.0d0/rnorm_FP64
     tmp2_FP64 = 0.0d0
     parcel_int = 1
     call set_val_FP64(pcagmres,tmp1,tmp1_FP64,parcel_int)
     call set_val_FP64(pcagmres,tmp2,tmp2_FP64,parcel_int)

     ! v=r/||r||
     call axpy(pcagmres,tmp1,tmp2,r,V,n)
     etime=MPI_Wtime()
     rtime(1) =  rtime(1) + etime-stime

     do k = 1,t
        ! call MPI_Barrier(icomm,ierr)
        stime=MPI_Wtime()
        ! 基底変換行列を作成
        if( (basis.eq.KSP_MPK_NEWTON_BASIS).and.(iflag_basis.eq.1) )then
#if 0
           do i =1,size(EigenReal_FP64)
              if(myrank.eq.0)then
                 write(*,*)i," : EigenReal : ",EigenReal_FP64(i)
              endif
           enddo
#endif
           call calc_NewtonBasis(s,B1_FP64,EigenReal_FP64,EigenImag_FP64)
           parcel_int = size(B1_FP64)
           call set_vec_FP64(pcagmres,B1,B1_FP64,parcel_int)
        endif

        if(k.eq.t)then
           iQmatflag=KSP_QRDECOMPOSITION_FALSE_ORTHOGONAL_MATRIX
        else
           iQmatflag=KSP_QRDECOMPOSITION_TRUE_ORTHOGONAL_MATRIX
        endif

        do i=1,s
           ! Mpv=v_i
           V%offset = (i-1)*vecLen
           call calc_precon(pcagmres,V,pV)
           V%offset = 0

#if 0
           V%offset = 0
           call calc_dot(pcagmres,icomm,V,V,tmp1,n)
           V%offset = 0
           parcel_int = 1
           call get_vec_FP64(pcagmres,tmp1,Rtri_FP64,parcel_int)
           if(myrank.eq.0)then
              write(*,*)" dot <V|V> @ 1 @ ",i,Rtri_FP64(1)
           endif
#endif
           ! V_i+1=ApV
           if(iflag_basis.eq.1)then
              select case (basis)
              case(KSP_MPK_NEWTON_BASIS)
                 call calc_NewtonBasisVec(pcagmres,n,s,i,vecLen,pv,V,B1_FP64)
              case(KSP_MPK_MONOMIAL_BASIS)
                 V%offset = i*vecLen
                 call calc_Matvec( pcagmres,pv,V )
                 V%offset = 0
              end select
           else
              ! 初回反復は固有値の情報がないため単基底
              V%offset = i*vecLen
              call calc_Matvec( pcagmres,pv,V )
              V%offset = 0
           endif
        enddo

        etime=MPI_Wtime()
        rtime(2) =  rtime(2) + etime-stime

        stime=MPI_Wtime()

        if(k.eq.1)then
           ! ------
           ! QR分解
           s1=s+1
           call QRDecomposition_( &
                pcagmres, &
                n, &
                s+1, &
                icomm, &
                V, vecLen,&
                Q, vecLen,&
                Rtri, s+1,&
                iQRflag, &
                iQmatflag, &
                nvecBlock,   &       
                localnRow, &
                batchCount, &
                V_s1_batch_prt, &
                work_s1_batch_prt, &
                work &
                )
        else

           ! ------
           !  ブロック古典グラムシュミット直交化
           ! ------
           !  R_' = Q_k-1 * Vk
           s1=s*(k-1)+1
           s2=s+1
           ! call calc_gemm_dot_CPU(&
           !      n,s1,s2,icomm,&
           !      Q, n,&
           !      V, n,&
           !      BlockR_ubar_prime, s1,&
           !      nvecBlock)
           parcel_int = size(BlockR_ubar_prime_FP64)
           call set_vec_FP64(pcagmres, &
                BlockR_ubar_prime, &
                BlockR_ubar_prime_FP64, &
                parcel_int)
           call calc_gemm_dot(&
                pcagmres, &
                icomm,&
                n,s1,s2,&
                Q    , vecLen,&
                V    , vecLen,&
                BlockR_ubar_prime, s1,&
                nvecBlock)

           !  Vk = Vk - QR_'
           s1 = s*(k-1)+1
           parcel_int = size(BlockR_ubar_prime_FP64)
           call get_vec_FP64(pcagmres,BlockR_ubar_prime,BlockR_ubar_prime_FP64,parcel_int)
           BlockR_ubar_prime_FP64 = -BlockR_ubar_prime_FP64
           parcel_int = size(BlockR_ubar_prime_FP64)
           call set_vec_FP64(pcagmres,BlockR_ubar_prime,BlockR_ubar_prime_FP64,parcel_int)

           V%offset                 = vecLen
           BlockR_ubar_prime%offset = s*(k-1)+1

              call calc_gemm(&
                   pcagmres, &
                   Q, vecLen, &
                   V, vecLen, &
                   !                V, vecLen, &
                   BlockR_ubar_prime, s1, &
                   n, &
                   s, s1, &
                   nvecBlock)
              V%offset                 = 0
              BlockR_ubar_prime%offset = 0
              
              BlockR_ubar_prime_FP64 = -BlockR_ubar_prime_FP64
              parcel_int = size(BlockR_ubar_prime_FP64)

           call set_vec_FP64(pcagmres,BlockR_ubar_prime,BlockR_ubar_prime_FP64,parcel_int)

#if 0
           do i =1,s+1
              V%offset = vecLen*(i-1)
              call calc_dot(pcagmres,icomm,V,V,tmp1,n)
              V%offset = 0
              call get_vec_FP64(pcagmres,tmp1,Rtri_FP64,1)
              if(myrank.eq.0)then
                 write(*,*)" dot <V|V> ",i,Rtri_FP64(1)
              endif
           enddo
#endif

           ! QR分解
           V%offset = vecLen
           Q%offset = vecLen*(s*(k-1)+1)
           call QRDecomposition_( &
                pcagmres, &
                n, &
                s, &
                icomm, &
                V, vecLen,&
                Q, vecLen,&
                Rtri, s,&
                iQRflag, &
                iQmatflag, &
                nvecBlock,   &       
                localnRow, &
                batchCount, &
                V_s_batch_prt, &
                work_s_batch_prt, &
                work &
                )
           V%offset = 0
           Q%offset = 0

        endif
        etime=MPI_Wtime()
        rtime(3) =  rtime(3) + etime-stime

        ! ------
        ! ヘッセンベルグ行列の計算
        parcel_int = size(Rtri_FP64)
        call get_vec_FP64(pcagmres,Rtri,Rtri_FP64,parcel_int)
        parcel_int = size(BlockR_ubar_prime_FP64)
        call get_vec_FP64(pcagmres,BlockR_ubar_prime,BlockR_ubar_prime_FP64,parcel_int)
        call calc_HessenbergMatrix( &
             s,t,k, &
             B1_FP64, &
             HessenbergMat_FP64, &
             Rtri_FP64, &
             BlockR_ubar_prime_FP64, &
             hk_pre_FP64 &
             )

     ! --------------------------------

        ! -------
        if(k.ne.t)then
           ! 次の基底生成の初期ベクトル
           Q%offset = vecLen*( s*(k-1)+1+s-1 ) 
           call xcpy(pcagmres,Q,V,n)
           Q%offset = 0
        endif

        ! -----
        ! ギブンス回転
        do i=1+s*(k-1),s*k
           do j=1,1+s*k
              HessenbergMat_Givens_FP64( j + (i-1)*(s*t+1) ) = &
                   HessenbergMat_FP64( j + (i-1)*(s*t+1) )
           enddo
        enddo
!  call parcel_sync(pcagmres)

        do i=1+s*(k-1),s*k
           mm = i
           s1=s*t
!           write(*,*)"i = ",i,s1
           call GivensRotations(s1,HessenbergMat_Givens_FP64,i,co_FP64,si_FP64,zeta_FP64)
           rnorm_FP64 = abs(zeta_FP64(i+1))

           itr0=i+(itr-1)+1
           pcagmres%reshistory(itr0)=rnorm_FP64/bnorm_FP64

#if 1
           ! コメントアウトすればb-Axで離脱する
           if( (rnorm_FP64/bnorm_FP64<rtolmax) ) then
              iflag_exit=1
              exit
           endif
           if( (i+(itr-1)+1).eq.itrmax)then
              iflag_exit=1
              exit
           endif
#endif
        enddo

 
        if(iflag_exit.eq.1)then
           exit
        endif

         if(basis.eq.KSP_MPK_NEWTON_BASIS)then
           iflag_basis=1
           ! 固有値計算
           s1=s*k
           s2=s*t
           call getHessenbergMat_Eigen(&
                s1, s2 ,&
                HessenbergMat_FP64,&
                EigenReal_FP64,EigenImag_FP64)
           s1_int4=s*k
           call modified_leja_ordering(s1_int4,EigenReal_FP64,EigenImag_FP64)
           !              if(myrank.eq.0)then
           !                do i=1,s*k
           !                   write(*,*)"Eigen ",i,EigenReal(i),EigenImag(i)
           !                enddo
           !             end if
        endif

     enddo

     parcel_int = size(Rtri_FP64)
     call set_vec_FP64(pcagmres,  Rtri, Rtri_FP64, parcel_int)
     parcel_int = size(HessenbergMat_Givens_FP64)
     call set_vec_FP64(pcagmres,  HessenbergMat_Givens, HessenbergMat_Givens_FP64, parcel_int)
     parcel_int = size(zeta_FP64)
     call set_vec_FP64(pcagmres,  zeta, zeta_FP64, parcel_int)

#if 0
     do i =1,size(HessenbergMat_Givens_FP64)
        if(myrank.eq.0)then
           write(*,*)"Hessenberg : ",i,HessenbergMat_Givens_FP64(i)
        endif
     enddo
     do i =1,size(zeta_FP64)
        if(myrank.eq.0)then
           write(*,*)"zeta : ",i,zeta_FP64(i)
        endif
     enddo
#endif


     ! 解ベクトルの更新
     stime=MPI_Wtime()
     select case (iQRflag)
     case(KSP_CAGMRES_QR_TYPE_CGS)
        ! call get_vec_FP64(pcagmres,HessenbergMat_Givens,HessenbergMat_Givens_FP64,size(HessenbergMat_Givens_FP64))
        ! call get_vec_FP64(pcagmres,zeta                ,zeta_FP64                ,size(zeta_FP64))
        s1=s*t
        call update_xVec(&
             pcagmres, &
             n, &
             s1, &
             mm, &
             vecLen,&
             x, &
             Q, &
             x0, &
             v0, &
             HessenbergMat_Givens_FP64, &
             zeta_FP64, &
             zeta , &
             nvecBlock)

     case(KSP_CAGMRES_QR_TYPE_MGS)
        s1=s*t
        ! call get_vec_FP64(pcagmres,HessenbergMat_Givens,HessenbergMat_Givens_FP64,size(HessenbergMat_Givens_FP64))
        ! call get_vec_FP64(pcagmres,zeta                ,zeta_FP64                ,size(zeta_FP64))
        call update_xVec( &
             pcagmres, &
             n, &
             s1, &
             mm, &
             vecLen,&
             x, &
             Q, &
             x0, &
             v0, &
             HessenbergMat_Givens_FP64, &
             zeta_FP64, &
             zeta, &
             nvecBlock)

     case default
        if(iQmatflag.eq.KSP_QRDECOMPOSITION_FALSE_ORTHOGONAL_MATRIX)then
           s3=1
           if(t.eq.s3)then
              ! call get_vec_FP64(pcagmres,HessenbergMat_Givens,HessenbergMat_Givens_FP64,size(HessenbergMat_Givens_FP64))
              ! call get_vec_FP64(pcagmres,zeta                ,zeta_FP64                ,size(zeta_FP64))
              ! call get_vec_FP64(pcagmres,Rtri                ,Rtri_FP64                ,size(Rtri_FP64))

              ! t=1 の場合は三角行列が(s+1) * (s+1) 
              s1=s*t
              call update_xVec_cagmres(&
                   pcagmres,&
                   n,s1,mm,s+1,&
                   x,  vecLen,&
                   Q,  vecLen,&
                   V,  vecLen,&
                   x0, vecLen,&
                   v0, vecLen,&
                   HessenbergMat_Givens_FP64, s+1, &
                   zeta_FP64, (s*t+1),&
                   Rtri_FP64, (s+1)  ,&
                   y, &
                   nvecBlock)
           else
              s1=1
              if(s.eq.s1)then
                 ! s=1の場合は解ベクトルの更新に利用するベクトルをすでに計算済み
                 s2=s*t
!                 call get_vec_FP64(pcagmres,HessenbergMat_Givens,HessenbergMat_Givens_FP64,size(HessenbergMat_Givens_FP64))
!                 call get_vec_FP64(pcagmres,zeta                ,zeta_FP64                ,size(zeta_FP64))
                 call update_xVec( &
                       pcagmres, &
                       n, &
                       s2, &
                       mm, &
                       vecLen,&
                       x, &
                       Q, &
                       x0, &
                       v0, &
                       HessenbergMat_Givens_FP64, &
                       zeta_FP64, &
                       zeta, &
                       nvecBlock)
               else
                 s2=s*t
                 V%offset = vecLen
!                 call get_vec_FP64(pcagmres,HessenbergMat_Givens,HessenbergMat_Givens_FP64,size(HessenbergMat_Givens_FP64))
!                 call get_vec_FP64(pcagmres,zeta                ,zeta_FP64                ,size(zeta_FP64))
!                 call get_vec_FP64(pcagmres,Rtri                ,Rtri_FP64                ,size(Rtri_FP64))
                 call update_xVec_cagmres(&
                      pcagmres,&
                      n,s2,mm,s,&
                      x , vecLen,&
                      Q , vecLen,&
                      V , vecLen,&
                      x0, vecLen,&
                      v0, vecLen,&
                      HessenbergMat_Givens_FP64, s, &
                      zeta_FP64, (s*t+1), &
                      Rtri_FP64, s+1, &
                      y, &
                      nvecBlock)
                 V%offset = 0
              endif
           endif
        else
           s2=s*t
           ! call get_vec_FP64(pcagmres,HessenbergMat_Givens,HessenbergMat_Givens_FP64,size(HessenbergMat_Givens_FP64))
           ! call get_vec_FP64(pcagmres,zeta                ,zeta_FP64                ,size(zeta_FP64))
           call update_xVec(&
                pcagmres, &
                n, &
                s2, &
                mm, &
                vecLen,&
                x, &
                Q, &
                x0, &
                v0, &
                HessenbergMat_Givens_FP64, &
                zeta_FP64, &
                zeta, &
                nvecBlock)
        endif
     end select

     etime=MPI_Wtime()
     rtime(4) =  rtime(4) + etime-stime

     if(iflag_exit.eq.1)then
        exit
     endif

  enddo
  etime0=MPI_Wtime()

  ! ------
  ! 真の残差を確認
  ! 真の残差を確認
#if 1
  parcel_int = 1
  call calc_res(pcagmres,x,b,r)
  call calc_dot(pcagmres,icomm,r,r,rnorm,n)
  call sqrt_val(pcagmres,rnorm,parcel_int)
  call get_val_FP64(pcagmres,rnorm,rnorm_FP64,parcel_int)

  if(myrank.eq.0)then
!     write(*,*)"@ true res @",rnorm_FP64,rnorm_FP64/bnorm_FP64
  endif
#endif

#ifdef timer
  if(myrank.eq.0)then
     write(*,*)"cagmres      ",etime0-stime0
     write(*,*)"    res      ",rtime(1)
     write(*,*)"    MatVec   ",rtime(2)
     write(*,*)"    QR       ",rtime(3)
     write(*,*)"    update x ",rtime(4)
  endif
#endif
   
  call parcel_sync(pcagmres)
#ifdef CUDA_SOLVER
  call cudaMemcpy_d2d_FP64(x%C_vec_FP64,pcagmres%C_x,n,offset_zero,offset_zero)
#else
  call get_vec_FP64(pcagmres,x,pcagmres%x,n)
#endif
  call parcel_sync(pcagmres)

  pcagmres%niter=min(itrmax,itr0)

  ! -------
  ! todo 射影した空間で残差ノルムが十分小さくなっても、
  ! 元の空間の真の残差ノルムが十分小さくなっているとは限らないので
  ! リターンコードを追加する
  if(itrmax<=itr)then
     pcagmres%iret=PARCEL_RETURNCODE_MAXITERATION
  else
     pcagmres%iret=PARCEL_RETURNCODE_NORMAL
  endif

  if(pcagmres%flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
     call cublas_free_DgemmBatched_gemm_dot_local_FP64(V_s_batch_prt)
     call cublas_free_DgemmBatched_gemm_dot_local_FP64(V_s1_batch_prt)
     call cublas_free_DgemmBatched_gemm_dot_local_FP64(work_s_batch_prt)
     call cublas_free_DgemmBatched_gemm_dot_local_FP64(work_s1_batch_prt)
  endif

  call free_vec(pcagmres,x)
  call free_vec(pcagmres,r)
  call free_vec(pcagmres,b)
  call free_vec(pcagmres,x0)
  call free_vec(pcagmres,v0)
  call free_vec(pcagmres,pV)
  call free_vec(pcagmres,V)
  call free_vec(pcagmres,Q)
  call free_vec(pcagmres,work)
  call free_vec(pcagmres,B1)
  call free_vec(pcagmres,EigenReal)
  call free_vec(pcagmres,EigenImag)
  call free_vec(pcagmres,HessenbergMat)
  call free_vec(pcagmres,HessenbergMat_Givens)
  call free_vec(pcagmres,Rtri)
  call free_vec(pcagmres,BlockR_ubar_prime)
  call free_vec(pcagmres,zeta)
  call free_vec(pcagmres,y)
  call free_vec(pcagmres,co)
  call free_vec(pcagmres,si)
  call free_vec(pcagmres,bnorm)
  call free_vec(pcagmres,rnorm)
  call free_vec(pcagmres,tmp1 )
  call free_vec(pcagmres,tmp2 )

  ! --------------------------------

  deallocate(zeta_FP64)
  deallocate(B1_FP64)
  deallocate(EigenReal_FP64)
  deallocate(EigenImag_FP64)
  deallocate(HessenbergMat_FP64)
  deallocate(HessenbergMat_Givens_FP64)
  deallocate(Rtri_FP64)
  deallocate(BlockR_ubar_prime_FP64)
  deallocate(co_FP64)
  deallocate(si_FP64)

end subroutine KSP_PCAGMRES

#ifdef CUDA_SOLVER
subroutine KSP_pCAARNOLDI(pcaarnoldi)
  ! use krylov
  ! use krylov_kernel
  ! use vec_kernel
  ! use QR_kernel
  ! use KSP_CAGMRES_CPU
  ! use KSP_CAARNOLDI

  ! use PARCEL_integer

  use krylov,only : KSP
  use PARCEL_define,only : PARCEL_CPU_FLAG_ON,PARCEL_GPU_FLAG_ON
  use KSP_STATUS,only : &
       STATUS_FLAG, &
       init_staus_flag,set_machineflag, &
       set_dataflag,set_calcflag

  use PARCEL_KRYLOV_MATRIX,only : calc_Matvec,calc_res
  use PARCEL_KRYLOV_PRECON,only : calc_precon
  use PARCEL_VEC_kernel
  use vec_type,only : PARCEL_VEC
  use PARCEL_integer,only : intType
  use KSP_CAGMRES,only : update_xVec,update_xVec_cagmres,calc_NewtonBasisVec
  use PARCEL_define,only : &
       PARCEL_CPU_FLAG_ON,PARCEL_GPU_FLAG_ON, &
       KSP_QRDECOMPOSITION_MGS, &
       KSP_QRDECOMPOSITION_CGS 
  use QR_kernel,only : &
       QRDecomposition_, &
       GS_QR_,GS_QR, &
       CholeskyQR_Ronly_,calcQ_VRinv_

  use KSP_CAGMRES_CPU,only : calc_HessenbergMatrix, calc_NewtonBasis,getHessenbergMat_Eigen

  use KSP_GMRES_CPU,only : GivensRotations
  use cu_common_wrapper,only : &
       cublas_init_DgemmBatched_gemm_dot_local_FP64, &
       cublas_free_DgemmBatched_gemm_dot_local_FP64

  use PARCEL_define,only : WARP_SIZE
  use c_wrapper,only : CUDA_NUM_BLOCKS

  implicit none

  type(KSP) :: pcaarnoldi
  integer(kind=intType) n
  integer icomm
  integer myrank ! プロセスランク
  integer(kind=intType) s ! 省通信ステップ 
  integer(kind=intType) t ! リスタート数 = (s*t)
  integer(kind=intType) ks,ke
  integer iQRflag ! QR分解の方法 
  integer basis ! 基底ベクトルの生成方法

  integer iQmatflag
  integer iEigenVecflag
  ! ----------
  integer itr
  integer(kind=intType) k
  integer(kind=intType) i,j
  integer(kind=intType) msize
  integer nvecBlock

  type(PARCEL_VEC) :: x
  type(PARCEL_VEC) :: b
  type(PARCEL_VEC) :: r

  type(PARCEL_VEC) :: pV
  type(PARCEL_VEC) :: V
  type(PARCEL_VEC) :: Q

  type(PARCEL_VEC) :: bnorm,rnorm
  type(PARCEL_VEC) :: tmp1,tmp2
  type(PARCEL_VEC) :: B1 
  type(PARCEL_VEC) :: EigenReal
  type(PARCEL_VEC) :: EigenImag

  type(PARCEL_VEC) :: HessenbergMat
  type(PARCEL_VEC) :: Rtri 
  type(PARCEL_VEC) :: BlockR_ubar_prime 
  type(PARCEL_VEC) :: HessenbergMat_EigenVector


  ! --------------------------------
  ! Batch GEMM
  ! V x V 
  integer(kind=intType) batchCount
  integer(kind=intType) localnRow
  integer(kind=intType) batch_len
  type(devicemem) V_s_batch_prt  
  type(devicemem) V_s1_batch_prt 
  type(devicemem) work_s_batch_prt 
  type(devicemem) work_s1_batch_prt 
  ! --------------------------------
  type(PARCEL_VEC) :: work ! Batch GEMM,TSQRで使用
  integer(kind=intType) parcel_int

  real*8 bnorm_FP64,rnorm_FP64
  real*8 tmp1_FP64,tmp2_FP64
  real*8,allocatable :: zeta_FP64(:)
  real*8,allocatable :: B1_FP64(:)
  real*8,allocatable :: EigenReal_FP64(:)
  real*8,allocatable :: EigenImag_FP64(:)
  real*8,allocatable :: HessenbergMat_FP64(:)
  real*8,allocatable :: HessenbergMat_Givens_FP64(:)
  real*8,allocatable :: Rtri_FP64(:)
  real*8,allocatable :: BlockR_ubar_prime_FP64(:)
  real*8,allocatable :: co_FP64(:)
  real*8,allocatable :: si_FP64(:)
  real*8,allocatable :: HessenbergMat_EigenVector_FP64(:)

  real*8 hk_pre_FP64
  real*8 alpha_FP64,beta_FP64

  integer(kind=intType) vecLen
  integer(kind=intType) nColL
  integer(kind=intType) nColU

  integer iflag_basis
  real*8 hk_pre
  integer ierr
  real*8 alpha,beta
  integer(kind=intType) s1,s2
  integer s1_int4

  integer(kind=intType),parameter :: offset_zero=0

#if 1
  ! GPU版はArnordi法に未対応
  call power_method_symMatrix(pcaarnoldi)

#else
  myrank = pcaarnoldi%commVec%myrank
  icomm = pcaarnoldi%commVec%comm
  n = pcaarnoldi%prm%n
  s = pcaarnoldi%CAARNOLDI_sstep
  t = pcaarnoldi%CAARNOLDI_tstep
  iQRflag  =  pcaarnoldi%CAARNOLDI_QRflag
  iEigenVecflag  =  1
  basis    =  pcaarnoldi%CAARNOLDI_basis
  nvecBlock =  pcaarnoldi%nvecBlock

  ! ------------------------------
  call MPI_Barrier(icomm,ierr)

  call parcel_sync(pcaarnoldi)
  nColL = 0
  nColU = 0
  if(pcaarnoldi%prm%MatFormat.eq.MATRIX_FORMAT_CRS)then
     nColL = pcaarnoldi%crsMat%nColL
     nColU = pcaarnoldi%crsMat%nColU
  endif
  vecLen = n+nColL+nColU

  if(pcaarnoldi%flag%machineflag.eq.PARCEL_CPU_FLAG_ON)then
     vecLen = n
  endif

  ! --------------------------------
  ! CUBLAS Batch paameter set
  ! CRS,DIA -- DDM で小行列のサイズを変更する
  batchCount = 1
  if(pcaarnoldi%flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
     !     localnRow = 32
     localnRow = 64 ! s =4
     !     localnRow = 128
     batch_len = n
     if( mod(n,localnRow).ne.0 )then
        batch_len = batch_len + ( localnRow - mod(n,localnRow) )
     endif
     batchCount = batch_len / localnRow
     vecLen = max(vecLen,batch_len)
  endif
  ! --------------------------------
  call init_vec(pcaarnoldi,x,vecLen)
  call init_vec(pcaarnoldi,x,vecLen)
  call init_vec(pcaarnoldi,r,vecLen)
  call init_vec(pcaarnoldi,pV,vecLen)
  call init_vec(pcaarnoldi,V,(vecLen*(s+1)) )
  call init_vec(pcaarnoldi,Q,(vecLen*(s*t+1)) )
  call init_vec(pcaarnoldi,B1,((s+1)*s))
  call init_vec(pcaarnoldi,EigenReal,(s*t))
  call init_vec(pcaarnoldi,EigenImag,(s*t))
  call init_vec(pcaarnoldi,HessenbergMat_EigenVector,(s*t * s*t))
  call init_vec(pcaarnoldi,HessenbergMat,(  (s*t+1)*s*t ))
  call init_vec(pcaarnoldi,Rtri,((s+1)*(s+1)))
  call init_vec(pcaarnoldi,BlockR_ubar_prime,( (s*t+1)*(s+1) ) )
  call init_vec(pcaarnoldi,work, &
       max( &
       batchCount*(s+1)*(s+1), &
       ((s+1)*(s+1)+(s+1)*((s+1)+WARP_SIZE))*(CUDA_NUM_BLOCKS) &
       ) &
       )

  call init_vec(pcaarnoldi,bnorm,1)
  call init_vec(pcaarnoldi,rnorm,1)
  call init_vec(pcaarnoldi,tmp1 ,1)
  call init_vec(pcaarnoldi,tmp2 ,1)

  call set_val_FP64(pcaarnoldi,tmp1,0.0d0,1)
  call set_val_FP64(pcaarnoldi,bnorm,0.0d0,1)
  call set_val_FP64(pcaarnoldi,rnorm,0.0d0,1)
  call set_val_FP64(pcaarnoldi,tmp1 ,0.0d0,1)
  call set_val_FP64(pcaarnoldi,tmp2 ,0.0d0,1)

  if( allocated(pcaarnoldi%EigenValueReal) ) deallocate(pcaarnoldi%EigenValueReal)
  if( allocated(pcaarnoldi%EigenValueImag) ) deallocate(pcaarnoldi%EigenValueImag)
  if( allocated(pcaarnoldi%EigenErrorNorm) ) deallocate(pcaarnoldi%EigenErrorNorm)
  allocate(pcaarnoldi%EigenValueReal(s*t))
  allocate(pcaarnoldi%EigenValueImag(s*t))
  allocate(pcaarnoldi%EigenErrorNorm(s*t))

  call set_vec_val(pcaarnoldi,x ,tmp1,1,vecLen)
  call set_vec_val(pcaarnoldi,r ,tmp1,1,vecLen)
  call set_vec_val(pcaarnoldi,pv,tmp1,1,vecLen)
  msize=n*(s+1)
  call set_vec_val(pcaarnoldi,V,tmp1,1,msize)
  msize=n*(s*t+1)
  call set_vec_val(pcaarnoldi,Q,tmp1,1,msize)
  msize = (s*t * s*t)
  call set_vec_val(pcaarnoldi,HessenbergMat_EigenVector,tmp1,1,msize)

  allocate( B1_FP64((s+1)*s) )
  allocate( EigenReal_FP64(s*t) )
  allocate( EigenImag_FP64(s*t) )
  allocate( HessenbergMat_FP64((s*t+1)*s*t) )
  allocate( HessenbergMat_Givens_FP64((s*t+1)*s*t) )
  allocate( Rtri_FP64((s+1)*(s+1)) )
  allocate( BlockR_ubar_prime_FP64((s*t+1)*(s+1)) )
  allocate( zeta_FP64((s*t+1)) )
  allocate( co_FP64((s*t+1)) ) 
  allocate( si_FP64((s*t+1)) ) 
  allocate( HessenbergMat_EigenVector_FP64( s*t*s*t) )

  if(pcaarnoldi%flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
     parcel_int = batchCount
     call cublas_init_DgemmBatched_gemm_dot_local_FP64(&
          parcel_int, &
          V%C_vec_FP64, V%offset, localnRow,&
          V_s1_batch_prt&
          ) 

     V%offset = vecLen
     parcel_int = batchCount
     call cublas_init_DgemmBatched_gemm_dot_local_FP64(&
          parcel_int, &
          V%C_vec_FP64, V%offset, localnRow,&
          V_s_batch_prt&
          ) 
     V%offset = 0

     ! work size : batchCount*s*s 
     parcel_int = batchCount
     call cublas_init_DgemmBatched_gemm_dot_local_FP64(&
          parcel_int, &
          work%C_vec_FP64, work%offset, s,  &
          work_s_batch_prt &
          ) 
     ! work size : batchCount*(s+1)*(s+1) 
     parcel_int = batchCount
     call cublas_init_DgemmBatched_gemm_dot_local_FP64(&
          parcel_int, &
          work%C_vec_FP64, work%offset, s+1,  &
          work_s1_batch_prt &
          ) 
  endif
  
  B1_FP64 = 0.0d0
  EigenReal_FP64 = 0.0d0
  EigenImag_FP64 = 0.0d0
  HessenbergMat_FP64 = 0.0d0
  HessenbergMat_Givens_FP64 = 0.0d0
  Rtri_FP64 = 0.0d0
  BlockR_ubar_prime_FP64 = 0.0d0
  zeta_FP64 = 0.0d0
  co_FP64 = 0.0d0
  si_FP64 = 0.0d0

  ! -------------

  call parcel_sync(pcaarnoldi)
#ifdef CUDA_SOLVER
  call cudaMemcpy_d2d_FP64(pcaarnoldi%C_x,x%C_vec_FP64,n,offset_zero,offset_zero)
  call cudaMemcpy_d2d_FP64(pcaarnoldi%C_b,b%C_vec_FP64,n,offset_zero,offset_zero)
#else
  call set_vec_FP64(pcaarnoldi,x,pcaarnoldi%x,n)
  call set_vec_FP64(pcaarnoldi,b,pcaarnoldi%b,n)
#endif
  call parcel_sync(pcaarnoldi)

  ! r = b-Ax
  call calc_res(pcaarnoldi,x,b,r)
  call calc_dot(pcaarnoldi,icomm,b,b,bnorm,n)
  call calc_dot(pcaarnoldi,icomm,r,r,rnorm,n)
  call sqrt_val(pcaarnoldi,bnorm,1)
  call sqrt_val(pcaarnoldi,rnorm,1)
  call get_val_FP64(pcaarnoldi,bnorm,bnorm_FP64,1)
  call get_val_FP64(pcaarnoldi,rnorm,rnorm_FP64,1)
  call parcel_sync(pcaarnoldi)

  ! ----------
  ! CACAARNOLDI のメインループ
  ! r = b-Ax
  ! do k = 1,t
  !   基底ベクトルの生成
  !   do i=1,s
  !     solve Ms=V_i
  !     Newton基底,単基底
  !     V_i+1=F(As)
  !   enddo
  !   if(k.eq.1)then
  !     QR分解
  !   else
  !     ブロック古典グラムシュミット直交化
  !     QR分解
  !   endif
  !   ヘッセンベルグ行列の計算
  ! enddo
  ! ----------


  B1_FP64   = 0.0d0
  Rtri_FP64 = 0.0d0
  zeta_FP64 = 0.0d0

  ! ----
  ! 基底変換行列
  do i=1,s
     B1_FP64(1+i+(i-1)*(s+1)) = 1.0d0
  enddo
  parcel_int = size(B1_FP64)
  call set_vec_FP64(pcaarnoldi,B1,B1_FP64,parcel_int)

  iflag_basis = -1

  EigenReal_FP64 = 0.0d0
  EigenImag_FP64 = 0.0d0
  alpha_FP64     = 0.0d0
  beta_FP64      = 0.0d0

  call set_vec_FP64(&
       pcaarnoldi, EigenReal, EigenReal_FP64, parcel_int)
  parcel_int = size(EigenImag_FP64)
  call set_vec_FP64(&
       pcaarnoldi, EigenImag, EigenImag_FP64, parcel_int)

  ! -----
  if(rnorm_FP64.eq.0.0d0)then
     call set_val_FP64(pcaarnoldi,tmp1,1.0d0,1)
     call set_vec_val(pcaarnoldi,r,tmp1,1,vecLen)
     ! !$OMP PARALLEL DO PRIVATE(i)
     ! do i=1,n
     !    r(i)=1.0d0
     ! enddo
     call calc_dot(pcaarnoldi,icomm,r,r,rnorm,n)
  endif
  call sqrt_val(pcaarnoldi,rnorm,1)
  call get_val_FP64(pcaarnoldi,rnorm,rnorm_FP64,1)

  tmp1_FP64 = -1.0d0/rnorm_FP64
  tmp2_FP64 = 0.0d0
  call set_val_FP64(pcaarnoldi,tmp1,tmp1_FP64,1)
  call set_val_FP64(pcaarnoldi,tmp2,tmp2_FP64,1)
  ! v=r/||r||
  call axpy(pcaarnoldi,tmp1,tmp2,r,V,n)

  do itr = 1,1
     ks = 1
     ke = t
    
     do k = ks,ke

        ! 基底変換行列を作成
        if( (basis.eq.KSP_MPK_NEWTON_BASIS).and.(iflag_basis.eq.1) )then
           call calc_NewtonBasis(s,B1_FP64,EigenReal_FP64,EigenImag_FP64)
           parcel_int = size(B1_FP64)
           call set_vec_FP64(pcaarnoldi,B1,B1_FP64,parcel_int)
        endif
        if(k.eq.t)then
           iQmatflag=KSP_QRDECOMPOSITION_FALSE_ORTHOGONAL_MATRIX
        else
           iQmatflag=KSP_QRDECOMPOSITION_TRUE_ORTHOGONAL_MATRIX
        endif
        do i=1,s
           ! Mpv=v_i
           V%offset = (i-1)*vecLen
           call calc_precon(pcaarnoldi,V,pV)
           V%offset = 0

           ! V_i+1=ApV
           if(iflag_basis.eq.1)then
              select case (basis)
              case(KSP_MPK_NEWTON_BASIS)
                 call calc_NewtonBasisVec(pcaarnoldi,n,s,i,vecLen,pv,V,B1_FP64)
              case(KSP_MPK_MONOMIAL_BASIS)
                 V%offset = i*vecLen
                 call calc_Matvec( pcaarnoldi,pv,V )
                 V%offset = 0
              end select
           else
              ! 初回反復は固有値の情報がないため単基底
              V%offset = i*vecLen
              call calc_Matvec( pcaarnoldi,pv,V )
              V%offset = 0
           endif
        enddo

        if(k.eq.1)then
           ! ------
           ! QR分解
           s1=s+1
           call QRDecomposition_( &
                pcaarnoldi, &
                n, &
                s+1, &
                icomm, &
                V, vecLen,&
                Q, vecLen,&
                Rtri, s+1,&
                iQRflag, &
                iQmatflag, &
                nvecBlock,   &       
                localnRow, &
                batchCount, &
                V_s1_batch_prt, &
                work_s1_batch_prt, &
                work &
                )
        else
           ! ------
           !  ブロック古典グラムシュミット直交化
           ! ------
           !  R_' = Q_k-1 * Vk
           s1=s*(k-1)+1
           s2=s+1
           parcel_int = size(BlockR_ubar_prime_FP64)
           call set_vec_FP64(pcaarnoldi, &
                BlockR_ubar_prime, &
                BlockR_ubar_prime_FP64, &
                parcel_int)
           call calc_gemm_dot(&
                pcaarnoldi, &
                icomm,&
                n,s1,s2,&
                Q    , vecLen,&
                V    , vecLen,&
                BlockR_ubar_prime, s1,&
                nvecBlock)

           !  Vk = Vk - QR_'
           s1=s*(k-1)+1
           s1 = s*(k-1)+1
           parcel_int = size(BlockR_ubar_prime_FP64)
           call get_vec_FP64(pcaarnoldi,BlockR_ubar_prime,BlockR_ubar_prime_FP64,parcel_int)
           BlockR_ubar_prime_FP64 = -BlockR_ubar_prime_FP64
           parcel_int = size(BlockR_ubar_prime_FP64)
           call set_vec_FP64(pcaarnoldi,BlockR_ubar_prime,BlockR_ubar_prime_FP64,parcel_int)

           V%offset                 = vecLen
           BlockR_ubar_prime%offset = s*(k-1)+1
           ! call CAGMRES_GEMM( &
           !      n,s,s1, &
           !      V,Q,BlockR_ubar_prime ,nvecBlock)

           call calc_gemm(&
                pcaarnoldi, &
                Q, vecLen, &
                V, vecLen, &
!                V, vecLen, &
                BlockR_ubar_prime, s1, &
                n, &
                s, s1, &
                nvecBlock)
           V%offset                 = 0
           BlockR_ubar_prime%offset = 0

           BlockR_ubar_prime_FP64 = -BlockR_ubar_prime_FP64
           parcel_int = size(BlockR_ubar_prime_FP64)
           call set_vec_FP64(pcaarnoldi,BlockR_ubar_prime,BlockR_ubar_prime_FP64,parcel_int)
           ! QR分解
           V%offset = vecLen
           Q%offset = vecLen*(s*(k-1)+1)
           call QRDecomposition_( &
                pcaarnoldi, &
                n, &
                s, &
                icomm, &
                V, vecLen,&
                Q, vecLen,&
                Rtri, s,&
                iQRflag, &
                iQmatflag, &
                nvecBlock,   &       
                localnRow, &
                batchCount, &
                V_s_batch_prt, &
                work_s_batch_prt, &
                work &
                )
           V%offset = 0
           Q%offset = 0

        endif

        ! ------
        ! ヘッセンベルグ行列の計算
        parcel_int = size(Rtri_FP64)
        call get_vec_FP64(pcaarnoldi,Rtri,Rtri_FP64,parcel_int)
        parcel_int = size(BlockR_ubar_prime_FP64)
        call get_vec_FP64(pcaarnoldi,BlockR_ubar_prime,BlockR_ubar_prime_FP64,parcel_int)
        call calc_HessenbergMatrix( &
             s,t,k, &
             B1_FP64, &
             HessenbergMat_FP64, &
             Rtri_FP64, &
             BlockR_ubar_prime_FP64, &
             hk_pre_FP64 &
             )

        ! -------
        if(k.ne.t)then
           ! 次の基底生成の初期ベクトル
           Q%offset = vecLen*( s*(k-1)+1+s-1 ) 
           call xcpy(pcaarnoldi,Q,V,n)
           Q%offset = 0
        endif
        if(basis.eq.KSP_MPK_NEWTON_BASIS)then
           iflag_basis=1
           ! 固有値計算
           s1=s*k
           s2=s*t
           call getHessenbergMat_Eigen(&
                s1, s2 ,&
                HessenbergMat_FP64,&
                EigenReal_FP64,EigenImag_FP64)
           s1_int4=s*k
           call modified_leja_ordering(s1_int4,EigenReal_FP64,EigenImag_FP64)
        endif

     enddo

     s1=s*t
     call getHessenbergMat_Eigen(&
          s1, s1 ,&
          HessenbergMat_FP64,&
          EigenReal_FP64,EigenImag_FP64,&
          HessenbergMat_EigenVector_FP64)
#if 0
     ! 初回は全て計算
     ! 次の反復はks-1個の固有値を元にしてリスタートする
     if(itr.eq.1)then
        ks = ks + ke/2 
        if(EigenImag(ks).ne.0.0d0)then
           if(ks+1<=ks)then
              if(EigenReal(ks).eq.EigenReal(ks+1))then
                 ks = ks - 1 
              endif
           endif
        endif
        ks = ks  + 1
     endif
#endif
     ! else
     !    call getHessenbergMat_Eigen(&
     !         s1, s1 ,&
     !         HessenbergMat,&
     !         EigenReal,EigenImag &
     !         )
     ! endif
  enddo

  call getRealEigen(&
       s1,&
       EigenReal,EigenImag,&
       pcaarnoldi%EigenMax_Real, &
       pcaarnoldi%EigenMax_Imag, &
       pcaarnoldi%EigenMin_Real, &
       pcaarnoldi%EigenMin_Imag  &
       )

  ! 固有ベクトルの計算
  !  ここの定数名はパラメータに変更
  if(iEigenVecflag.eq.1)then
     ! x = Q * HessenbergMat_EigenVector : xが欲しい固有ベクトル
     call getEigenVector(n,s1,Q,HessenbergMat_EigenVector,nvecBlock)
     if( allocated(pcaarnoldi%EigenVector) )then
        deallocate(pcaarnoldi%EigenVector)
     endif
     allocate(pcaarnoldi%EigenVector(n*(s*t)))
     call xcpy(n*(s*t),Q,pcaarnoldi%EigenVector)
  endif

  pcaarnoldi%EigenValueReal = 0.0d0
  pcaarnoldi%EigenValueImag = 0.0d0
  pcaarnoldi%EigenErrorNorm = -1.0d0
  call checkEigenPair(pcaarnoldi,n,s1,Q,&
       EigenReal,EigenImag, &
       pcaarnoldi%MaxEigenID,pcaarnoldi%MinEigenID,&
       pcaarnoldi%EigenErrorNorm)
  do i=1,s*t
     pcaarnoldi%EigenValueReal(i) = EigenReal(i)
     pcaarnoldi%EigenValueImag(i) = EigenImag(i)
  enddo

  if(0<pcaarnoldi%MaxEigenID)then
     pcaarnoldi%EigenMax_Real = EigenReal(pcaarnoldi%MaxEigenID) 
     pcaarnoldi%EigenMax_Imag = EigenImag(pcaarnoldi%MaxEigenID)
  else
     pcaarnoldi%EigenMax_Real = 0.0d0
     pcaarnoldi%EigenMax_Imag = 0.0d0
  endif
  if(0<pcaarnoldi%MinEigenID)then
     pcaarnoldi%EigenMin_Real = EigenReal(pcaarnoldi%MinEigenID) 
     pcaarnoldi%EigenMin_Imag = EigenImag(pcaarnoldi%MinEigenID)
  else
     pcaarnoldi%EigenMin_Real = 0.0d0
     pcaarnoldi%EigenMin_Imag = 0.0d0
  endif

  if(pcagmres%flag%machineflag.eq.PARCEL_GPU_FLAG_ON)then
     call cublas_free_DgemmBatched_gemm_dot_local_FP64(V_s_batch_prt)
     call cublas_free_DgemmBatched_gemm_dot_local_FP64(V_s1_batch_prt)
     call cublas_free_DgemmBatched_gemm_dot_local_FP64(work_s_batch_prt)
     call cublas_free_DgemmBatched_gemm_dot_local_FP64(work_s1_batch_prt)
  endif


  !  call view_vec(n,pcaarnoldi%EigenVector)

#if 0
  if(myrank.eq.0)then
     do i=1,s*t
        write(*,*)"Eigen ",i,EigenReal(i),EigenImag(i)
     enddo
     ! write(*,*)"EigenMax ",pcaarnoldi%EigenMax_Real,pcaarnoldi%EigenMax_Imag
     ! write(*,*)"EigenMin ",pcaarnoldi%EigenMin_Real,pcaarnoldi%EigenMin_Imag
  end if
#endif

  call free_vec(pcagmres,x)
  call free_vec(pcagmres,r)
  call free_vec(pcagmres,b)
  call free_vec(pcagmres,x0)
  call free_vec(pcagmres,v0)
  call free_vec(pcagmres,pV)
  call free_vec(pcagmres,V)
  call free_vec(pcagmres,Q)
  call free_vec(pcagmres,work)
  call free_vec(pcagmres,B1)
  call free_vec(pcagmres,EigenReal)
  call free_vec(pcagmres,EigenImag)
  call free_vec(pcagmres,HessenbergMat)
  call free_vec(pcagmres,HessenbergMat_Givens)
  call free_vec(pcagmres,Rtri)
  call free_vec(pcagmres,BlockR_ubar_prime)
  call free_vec(pcagmres,zeta)
  call free_vec(pcagmres,y)
  call free_vec(pcagmres,co)
  call free_vec(pcagmres,si)
  call free_vec(pcagmres,bnorm)
  call free_vec(pcagmres,rnorm)
  call free_vec(pcagmres,tmp1 )
  call free_vec(pcagmres,tmp2 )

  deallocate(zeta_FP64)
  deallocate(B1_FP64)
  deallocate(EigenReal_FP64)
  deallocate(EigenImag_FP64)
  deallocate(HessenbergMat_FP64)
  deallocate(HessenbergMat_Givens_FP64)
  deallocate(Rtri_FP64)
  deallocate(BlockR_ubar_prime_FP64)
  deallocate(co_FP64)
  deallocate(si_FP64)
#endif

end subroutine KSP_PCAARNOLDI

#else
subroutine KSP_pCAARNOLDI(pcaarnoldi)
  use krylov
  use krylov_kernel
  use vec_kernel
  use QR_kernel
  use KSP_CAGMRES_CPU
  use KSP_CAARNOLDI

  use PARCEL_integer
  implicit none

  type(KSP) :: pcaarnoldi
  integer(kind=intType) n
  integer icomm
  integer myrank ! プロセスランク
  integer(kind=intType) s ! 省通信ステップ 
  integer(kind=intType) t ! リスタート数 = (s*t)
  integer(kind=intType) ks,ke
  integer iQRflag ! QR分解の方法 
  integer basis ! 基底ベクトルの生成方法

  integer iQmatflag
  integer iEigenVecflag
  ! ----------
  integer itr
  integer(kind=intType) k
  integer(kind=intType) i,j
  integer(kind=intType) msize
  integer nvecBlock

  real*8,allocatable :: x(:)
  real*8,allocatable :: r(:)

  real*8,allocatable :: pV(:)
  real*8,allocatable :: V(:)
  real*8,allocatable :: Q(:)

  real*8 bnorm,rnorm
  real*8 tmp1,tmp2
  real*8,allocatable :: B1(:) 
  real*8,allocatable :: EigenReal(:)
  real*8,allocatable :: EigenImag(:)

  real*8,allocatable :: HessenbergMat(:)
  real*8,allocatable :: Rtri(:) 
  real*8,allocatable :: BlockR_ubar_prime(:) 
  real*8,allocatable :: HessenbergMat_EigenVector(:)
  integer iflag_basis
  real*8 hk_pre
  integer ierr
  real*8 alpha,beta
  integer(kind=intType) s1,s2
  integer s1_int4

  myrank=pcaarnoldi%commVec%myrank
  icomm=pcaarnoldi%commVec%comm
  n=pcaarnoldi%prm%n
  s=pcaarnoldi%CAARNOLDI_sstep
  t=pcaarnoldi%CAARNOLDI_tstep
  iQRflag = pcaarnoldi%CAARNOLDI_QRflag
  iEigenVecflag = 1
  basis   = pcaarnoldi%CAARNOLDI_basis
  nvecBlock= pcaarnoldi%nvecBlock

  allocate(x(n))
  allocate(r(n))
  allocate(pV(n))
  allocate(V(n*(s+1)))
  allocate(Q(n*(s*t+1)))
  allocate(B1((s+1)*s))
  allocate(EigenReal(s*t))
  allocate(EigenImag(s*t))
  allocate(HessenbergMat_EigenVector(s*t * s*t))
  allocate(HessenbergMat(  (s*t+1)*s*t ))
  allocate(Rtri((s+1)*(s+1)))
  allocate(BlockR_ubar_prime( (s*t+1)*(s+1) ) )

  if( allocated(pcaarnoldi%EigenValueReal) ) deallocate(pcaarnoldi%EigenValueReal)
  if( allocated(pcaarnoldi%EigenValueImag) ) deallocate(pcaarnoldi%EigenValueImag)
  if( allocated(pcaarnoldi%EigenErrorNorm) ) deallocate(pcaarnoldi%EigenErrorNorm)
  allocate(pcaarnoldi%EigenValueReal(s*t))
  allocate(pcaarnoldi%EigenValueImag(s*t))
  allocate(pcaarnoldi%EigenErrorNorm(s*t))

  call zero_initialize(n,x)
  call zero_initialize(n,r)
  call zero_initialize(n,pv)
  msize=n*(s+1)
  call zero_initialize(msize,V)
  msize=n*(s*t+1)
  call zero_initialize(msize,Q)
  HessenbergMat_EigenVector = 0.0d0

  call xcpy(n,pcaarnoldi%x,x)

  call calc_dot(n,icomm,pcaarnoldi%b,pcaarnoldi%b,bnorm)
  bnorm=sqrt(bnorm)
  B1=0.0d0
  Rtri=0.0d0

  ! ----------
  ! CACAARNOLDI のメインループ
  ! r = b-Ax
  ! do k = 1,t
  !   基底ベクトルの生成
  !   do i=1,s
  !     solve Ms=V_i
  !     Newton基底,単基底
  !     V_i+1=F(As)
  !   enddo
  !   if(k.eq.1)then
  !     QR分解
  !   else
  !     ブロック古典グラムシュミット直交化
  !     QR分解
  !   endif
  !   ヘッセンベルグ行列の計算
  ! enddo
  ! ----------

  ! ----
  ! 基底変換行列
  do i=1,s
     B1(1+i+(i-1)*(s+1))=1.0d0
  enddo

  iflag_basis=-1
  EigenReal=0.0d0
  EigenImag=0.0d0
  alpha=0.0d0
  beta =0.0d0



  HessenbergMat=0.0d0
  Rtri=0.0d0
  BlockR_ubar_prime=0.0d0
  hk_pre=0.0d0

  ! -----
  ! ここで残差を計算しているが,初期ベクトルなので変更可能
  ! r=b-Ax
  call calc_res(pcaarnoldi,x,pcaarnoldi%b,r)
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
  call axpy(n,tmp1,tmp2,r,V)
  do itr = 1,1
     ks = 1
     ke = t
     do k = ks,ke
        ! 基底変換行列を作成
        if( (basis.eq.KSP_MPK_NEWTON_BASIS).and.(iflag_basis.eq.1) )then
           call calc_NewtonBasis(s,B1,EigenReal,EigenImag)
        endif
        if(k.eq.t)then
           iQmatflag=KSP_QRDECOMPOSITION_FALSE_ORTHOGONAL_MATRIX
        else
           iQmatflag=KSP_QRDECOMPOSITION_TRUE_ORTHOGONAL_MATRIX
        endif
        do i=1,s
           ! MpV=v_i
           call calc_precon(pcaarnoldi,V(1+(i-1)*n),pV )
           ! V_i+1=ApV
           if(iflag_basis.eq.1)then
              select case (basis)
              case(KSP_MPK_NEWTON_BASIS)
                 call calc_NewtonBasisVec_CPU(pcaarnoldi,n,s,i,pv,V,B1)
              case(KSP_MPK_MONOMIAL_BASIS)
                 call calc_Matvec( pcaarnoldi,pv,V(1+i*n) )
              end select
           else
              ! 初回反復は固有値の情報がないため単基底
              call calc_Matvec( pcaarnoldi,pv,V(1+i*n) )
           endif
        enddo

        if(k.eq.1)then
           ! ------
           ! QR分解
           s1=s+1
           call QRDecomposition(n,s1,icomm,V,Q,Rtri,iQRflag,iQmatflag,nvecBlock)

        else
           ! ------
           !  ブロック古典グラムシュミット直交化
           ! ------
           !  R_' = Q_k-1 * Vk
           s1=s*(k-1)+1
           s2=s+1
           call calc_gemm_dot_CPU(&
                n,s1,s2,icomm,&
                Q, n,&
                V, n,&
                BlockR_ubar_prime, s1,&
                nvecBlock)

           !  Vk = Vk - QR_'
           s1=s*(k-1)+1
           call CAGMRES_GEMM(n,s,s1, &
                V(1+n),Q,BlockR_ubar_prime(1+s*(k-1)+1) ,nvecBlock)
           ! QR分解
           call QRDecomposition(n,s,icomm,V(1+n),Q(1+n*(s*(k-1)+1)),Rtri,iQRflag,iQmatflag,nvecBlock)
        endif

        ! ------
        ! ヘッセンベルグ行列の計算
        call calc_HessenbergMatrix( &
             s,t,k, &
             B1, &
             HessenbergMat, &
             Rtri, &
             BlockR_ubar_prime, &
             hk_pre &
             )

        ! -------
        if(k.ne.t)then
           ! 次の基底生成の初期ベクトル
           call xcpy(n,Q( 1 + n*( s*(k-1)+1+s-1 ) ),V)
        endif
        if(basis.eq.KSP_MPK_NEWTON_BASIS)then
           iflag_basis=1
           ! 固有値計算
           s1=s*k
           s2=s*t
           call getHessenbergMat_Eigen(&
                s1, s2 ,&
                HessenbergMat,&
                EigenReal,EigenImag)
           !           s1=s*k
           s1_int4=s*k
           call modified_leja_ordering(s1_int4,EigenReal,EigenImag)
        endif

     enddo

     s1=s*t
     ! !  ここの定数1はパラメータに変更
     ! if(iEigenVecflag.eq.1)then
     call getHessenbergMat_Eigen(&
          s1, s1 ,&
          HessenbergMat,&
          EigenReal,EigenImag,&
          HessenbergMat_EigenVector)
#if 0
     ! 初回は全て計算
     ! 次の反復はks-1個の固有値を元にしてリスタートする
     if(itr.eq.1)then
        ks = ks + ke/2 
        if(EigenImag(ks).ne.0.0d0)then
           if(ks+1<=ks)then
              if(EigenReal(ks).eq.EigenReal(ks+1))then
                 ks = ks - 1 
              endif
           endif
        endif
        ks = ks  + 1
     endif
#endif
     ! else
     !    call getHessenbergMat_Eigen(&
     !         s1, s1 ,&
     !         HessenbergMat,&
     !         EigenReal,EigenImag &
     !         )
     ! endif

  enddo

#if 0
  call getRealEigen(&
       s1,&
       EigenReal,EigenImag,&
       pcaarnoldi%EigenMax_Real, &
       pcaarnoldi%EigenMax_Imag, &
       pcaarnoldi%EigenMin_Real, &
       pcaarnoldi%EigenMin_Imag  &
       )
#endif

  ! 固有ベクトルの計算
  !  ここの定数名はパラメータに変更
  if(iEigenVecflag.eq.1)then
     ! x = Q * HessenbergMat_EigenVector : xが欲しい固有ベクトル
     call getEigenVector(n,s1,Q,HessenbergMat_EigenVector,nvecBlock)
     if( allocated(pcaarnoldi%EigenVector) )then
        deallocate(pcaarnoldi%EigenVector)
     endif
     allocate(pcaarnoldi%EigenVector(n*(s*t)))
     call xcpy(n*(s*t),Q,pcaarnoldi%EigenVector)

  endif

  pcaarnoldi%EigenValueReal = 0.0d0
  pcaarnoldi%EigenValueImag = 0.0d0
  pcaarnoldi%EigenErrorNorm = -1.0d0
  call checkEigenPair(pcaarnoldi,n,s1,Q,&
       EigenReal,EigenImag, &
       pcaarnoldi%MaxEigenID,pcaarnoldi%MinEigenID,&
       pcaarnoldi%EigenErrorNorm)
  do i=1,s*t
     pcaarnoldi%EigenValueReal(i) = EigenReal(i)
     pcaarnoldi%EigenValueImag(i) = EigenImag(i)
  enddo

  if(0<pcaarnoldi%MaxEigenID)then
     pcaarnoldi%EigenMax_Real = EigenReal(pcaarnoldi%MaxEigenID) 
     pcaarnoldi%EigenMax_Imag = EigenImag(pcaarnoldi%MaxEigenID)
  else
     pcaarnoldi%EigenMax_Real = 0.0d0
     pcaarnoldi%EigenMax_Imag = 0.0d0
  endif
  if(0<pcaarnoldi%MinEigenID)then
     pcaarnoldi%EigenMin_Real = EigenReal(pcaarnoldi%MinEigenID) 
     pcaarnoldi%EigenMin_Imag = EigenImag(pcaarnoldi%MinEigenID)
  else
     pcaarnoldi%EigenMin_Real = 0.0d0
     pcaarnoldi%EigenMin_Imag = 0.0d0
  endif

  !  call view_vec(n,pcaarnoldi%EigenVector)

#if 0
  if(myrank.eq.0)then
     do i=1,s*t
        write(*,*)"Eigen ",i,EigenReal(i),EigenImag(i)
     enddo
     ! write(*,*)"EigenMax ",pcaarnoldi%EigenMax_Real,pcaarnoldi%EigenMax_Imag
     ! write(*,*)"EigenMin ",pcaarnoldi%EigenMin_Real,pcaarnoldi%EigenMin_Imag
  end if
#endif

  if( allocated(x) ) deallocate(x)
  if( allocated(r) ) deallocate(r)
  if( allocated(pV) ) deallocate(pV)
  if( allocated(V) ) deallocate(V)
  if( allocated(Q) ) deallocate(Q)
  if( allocated(B1) ) deallocate(B1)
  if( allocated(EigenReal) ) deallocate(EigenReal)
  if( allocated(EigenImag) ) deallocate(EigenImag)
  if( allocated(HessenbergMat_EigenVector) ) deallocate(HessenbergMat_EigenVector)
  if( allocated(HessenbergMat) ) deallocate(HessenbergMat)
  if( allocated(Rtri) ) deallocate(Rtri)
  if( allocated(BlockR_ubar_prime) ) deallocate(BlockR_ubar_prime)

end subroutine KSP_PCAARNOLDI


#endif
