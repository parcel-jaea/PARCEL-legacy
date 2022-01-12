
module KSP_CALANCZOS
  use wrap_DGEMM,only : wrap_BASIC_DGEMM

contains
  subroutine checkEigenPair_lan(pcalanczos,n,s,x_re,x_im,e_re,e_im)
    use krylov
    use krylov_kernel,only : calc_Matvec
    use vec_kernel,only : zero_initialize,axpy,calc_dot
    use PARCEL_integer
    implicit none

    type(KSP) :: pcalanczos

    integer s
    integer n
    real*8 x_re(n*s)
    real*8 x_im(n*s)
    real*8 e_re(s)
    real*8 e_im(s)

    integer i
    integer icomm 
    integer myrank

    real*8 Ax_re(n)
    real*8 Ax_im(n)

    real*8 kx_re(n)
    real*8 kx_im(n)
    real*8 err_norm,err_im_norm
    real*8 kx_re_norm,kx_im_norm
    
    icomm  = pcalanczos%commVec%comm
    myrank = pcalanczos%commVec%myrank
    
    ! 固有対の精度チェックは||Ax-kx||/||kx||により行う.
    ! kは固有値.xは固有ベクトル.
    do i = 1,s
       call calc_Matvec( pcalanczos,x_re(1+(i-1)*n),Ax_re )
       call calc_Matvec( pcalanczos,x_im(1+(i-1)*n),Ax_im )

       ! 実数部の計算
       call axpy(n, e_re(i),0.0d0,x_re(1+(i-1)*n),kx_re)
       call axpy(n,-e_im(i),1.0d0,x_im(1+(i-1)*n),kx_re)

       ! 虚数部の計算
       call axpy(n, e_im(i),0.0d0,x_re(1+(i-1)*n),kx_im)
       call axpy(n, e_re(i),1.0d0,x_im(1+(i-1)*n),kx_im)

       call calc_dot(n,icomm,kx_re,kx_re,kx_re_norm)
       call calc_dot(n,icomm,kx_im,kx_im,kx_im_norm)
       kx_re_norm = kx_re_norm + kx_im_norm
       kx_re_norm = sqrt(abs(kx_re_norm))

       call axpy(n,1.0d0,-1.0d0,Ax_re,kx_re)
       call axpy(n,1.0d0,-1.0d0,Ax_im,kx_im)

       call calc_dot(n,icomm,kx_re,kx_re,err_norm)
       call calc_dot(n,icomm,kx_im,kx_im,kx_im_norm)
       err_norm = err_norm + kx_im_norm
       err_norm = sqrt(abs(err_norm))

       if(myrank.eq.0)then
          write(*,'(i8,e18.8,e18.8,e18.8)')i,err_norm/kx_re_norm,e_re(i),e_im(i)
       endif
       
    enddo

  end subroutine checkEigenPair_lan

  subroutine sum_vec(n, m, v, v0, y, nvecBlock)
    use vec_kernel
    use PARCEL_integer
    implicit none

    ! arguments
    integer(kind=intType) n
    integer(kind=intType) m
    integer nvecBlock
    real*8 v(n*(m+1))
    real*8 v0(n)
    real*8 y(m)

    ! local variables
    integer i
    integer j
    integer jj

    !$OMP PARALLEL DO PRIVATE(jj,j,i)
    do jj=1,n,nvecBlock
      do j=jj,min(n,jj+nvecBlock-1)
        v0(j)=0.0d0
      enddo
      do i=1,m
        do j=jj,min(n,jj+nvecBlock-1)
          v0(j)=v0(j)+y(i)*v(j+n*(i-1))
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine sum_vec

  subroutine calc_NewtonBasis_lan(s,B1,EigenReal,EigenImag)
    use PARCEL_integer
    implicit none
    integer(kind=intType) s
    real*8 B1(s*(s+1))
    real*8 EigenReal(s)
    real*8 EigenImag(s)
    integer i
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

  end subroutine calc_NewtonBasis_lan

  subroutine getTridiagonalMat_Eigen(&
       st,st0,&
       TriMat0,&
       EigenReal,EigenImag)
    use PARCEL_integer
    implicit none
    integer(kind=intType) st,st0
    real*8 TriMat0((st0+1)*st0 )
    real*8 TriMat(st*st)
    real*8 EigenReal(st)
    real*8 EigenImag(st)
#if 1
    real*8 D(st)
    real*8 E(st-1)
    real*8 work(5*st)
    integer iblock(st)
    integer isplit(st)
    integer iwork(3*st)
    integer info
    integer i,j
    integer m
    integer nsplit
    integer idum
    real*8 dum
    real*8 zero

    zero = 0.0d0

    do j=1,st
       do i=1,st
          TriMat(i+(j-1)*st)=TriMat0(i+(j-1)*(st0+1))
       enddo
       ! write(*,*)"------"
    enddo

    do j = 1, st
      D(j) = TriMat(j + (j - 1)*st)
    enddo
    do j = 2, st
      E(j-1) = TriMat(j-1 + (j - 1)*st)
      !write(*,*) "E", E(j-1)
    enddo

    call dstebz( &
         'A', 'E', &
         st, dum, dum, &
         idum, idum, &
         zero, D, E, &
         m, nsplit, EigenReal, &
         iblock, isplit, work, &
         iwork, info)
#else
    real*8 z(1) ! 'N'を指定しているので参照されない
    real*8 work(st)
    integer info
    integer i,j

    do j=1,st
       do i=1,st
          TriMat(i+(j-1)*st)=TriMat0(i+(j-1)*(st0+1))
        !  write(*,*)"Hes     @ ",i,j,HesMat(i+(j-1)*st)
       enddo
       ! write(*,*)"------"
    enddo
    call dhseqr( &
         'E', 'N', &
         st, 1, st, &
         TriMat, st, &
         EigenReal,EigenImag, &
         z,1,work,st,info)
#endif

  end subroutine getTridiagonalMat_Eigen

  subroutine getTridiagonalMat_EigenVector(&
             st,st0,&
             TriMat0,&
             EigenReal,EigenImag, &
             n,EigenVecReal,EigenVecImag)
    use PARCEL_integer
    implicit none
    integer(kind=intType) st,st0,n
    real*8 TriMat0((st0+1)*st0 )
    real*8 TriMat(st*st)
    real*8 EigenReal(st)
    real*8 EigenImag(st)
    real*8 EigenVecReal(n*st)
    real*8 EigenVecImag(n*st)
    real*8 D(st)
    real*8 E(st-1)
    real*8 work(5*st)
    integer iblock(st)
    integer isplit(st)
    integer iwork(3*st)
    integer ifail(st)
    integer info
    integer i,j

    do j=1,st
       do i=1,st
          TriMat(i+(j-1)*st)=TriMat0(i+(j-1)*(st0+1))
        !  write(*,*)"Tri     @ ",i,j,TriMat(i+(j-1)*st)
       enddo
       ! write(*,*)"------"
    enddo

    do j = 1, st
      D(j) = TriMat(j + (j - 1)*st)
    enddo
    do j = 2, st
      E(j-1) = TriMat(j-1 + (j - 1)*st)
      !write(*,*) "E", E(j-1)
    enddo

    !for eigenvector
    EigenVecReal = 0.0d0
    EigenVecImag = 0.0d0
    iblock = 1
    isplit(1) = st

    call dstein( &
         st, D, E, &
         st, EigenReal, &
         iblock, isplit, &
         EigenVecReal, n, &
         work, iwork, &
         ifail, info)

  end subroutine getTridiagonalMat_EigenVector


  ! --------------
  ! PARCELは実行列にのみ対応しているため固有値は実数あるいは複素共役対
  ! 最大,最小は絶対値最大,最小固有値
  subroutine getRealEigen_lan(&
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

    real*8 Emax,Emin

    integer i
    if(st<1)then
       EmaxReal=0.0d0
       EmaxImag=0.0d0
       EminReal=0.0d0
       EminImag=0.0d0
       return
    endif

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

  end subroutine getRealEigen_lan

  subroutine calc_betake1eskRRinv(s,k,betak,Rkk,Rinv,Tri)
    implicit none
    integer s
    integer k
    real*8 betak
    real*8 Rkk( (s*(k-1))*s )
    real*8 Rinv(s*s)
    real*8 Tri(s*s)

    real*8 e1est_Rkk(s*s)
    real*8 Tri0(s*s)
    integer i,j
    ! ---------
    do j=1,s*s
       e1est_Rkk(j)=0.0d0
       Tri0(j)=0.0d0
    enddo

    j=s*(k-1)
    do i=1,s
       e1est_Rkk(1 + (i-1)*s )=Rkk( j + (i-1) * ( s*(k-1) ) )*betak
    enddo
    call wrap_BASIC_DGEMM(s,s,s,e1est_Rkk,Rinv,Tri0)    
    do j=1,s*s
       Tri(j)=Tri(j)-Tri0(j)
    enddo

  end subroutine calc_betake1eskRRinv

  subroutine calc_Rk_Bk_Rkinv(s,R0,B,Rinv,RBRinv)
    implicit none

    integer s
    real*8 R0((s+1)*(s+1))
    real*8 B((s+1)*s)
    real*8 Rinv(s*s)
    real*8 RBRinv(s*s)

    integer i,j
    real*8 RB(s*s)
    real*8 B0(s*s)
    real*8 R(s*s)

    do i=1,s
       do j=1,s
          B0(j+(i-1)*s)=B(j+(i-1)*(s+1))
          R(j+(i-1)*s)=R0(j+(i-1)*(s+1))
       enddo
    enddo
    RB=0.0d0
    call DGEMM('N','N', s, s, s, 1.0d0, R ,s, B0,  s, 0.0d0,RB    ,s)
    call DGEMM('N','N', s, s, s, 1.0d0, RB,s, Rinv,s, 0.0d0,RBRinv,s)

  end subroutine calc_Rk_Bk_Rkinv

  subroutine calc_TridiagonalMatrix( &
       s,t,k, &
       B1, &
       Tri_all, &
       Rtri, &
       BlockR_ubar_prime, &
       betak_pre &
       )
    implicit none
    integer s,t,k

    real*8 rho_tilde,rho 
    real*8 betak
    real*8 bk
    real*8 betak_pre
    real*8 zk

    real*8 B1((s+1)*s)
    real*8 Tri_all(  (s*t+1)*s*t )
    real*8 Rtri((s+1)*(s+1))
    real*8 BlockR_ubar_prime( (s*t+1)*(s+1) ) 

    real*8 Rinv((s+1)*(s+1))
    real*8 Rtri0((s+1)*(s+1))
    real*8 BlockR_ubar( (s*t)*(s+1) ) 
    real*8 BlockR( (s*t+1)*s ) 
    real*8 Tri( (s+1)*(s+1) )
    real*8 BlockTri( (s*t+1)*s*t )
    real*8 BlockTri_R(s*(s*t))
    real*8 BlockR_ubar_B1(s*(s*t))
    real*8 R_B_( (s+1)*s )
    integer i,j
    integer ierr

    Rinv=0.0d0
    Rtri0=0.0d0
    BlockR_ubar=0.0d0
    BlockR=0.0d0
    Tri=0.0d0
    BlockTri=0.0d0
    BlockTri_R=0.0d0
    BlockR_ubar_B1=0.0d0
    R_B_=0.0d0

    if(k.eq.1)then
       do j=1,s
          do i=1,s
             Rinv(i+s*(j-1))=Rtri(i+(s+1)*(j-1))
          enddo
       enddo
       call DTRTRI('U', 'N', s, Rinv, s, ierr)
       ! R_B_=R_*B_
       call wrap_BASIC_DGEMM(s+1,s,s+1,Rtri,B1,R_B_)
       call wrap_BASIC_DGEMM(s+1,s,s,R_B_,Rinv,Tri)
       do j=1,s
          do i=1,s+1
             Tri_all(i + (j-1)*(s*t+1) ) = Tri( i + (j-1)*(s+1) )
          enddo
       enddo
       betak_pre=Tri((s+1)*s )
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
          enddo
       enddo
       ! -----
       ! R(1:s*(k-1)   ,1:s) = R_'(1:s*(k-1),1:s)
       ! R :: [s*(k-1)]行,[s]列
       do j=1,s
          do i=1,s*(k-1)
             BlockR(i+(j-1)*(s*(k-1)) )=BlockR_ubar_prime( i + (j-1)*(s*(k-1)+1) ) 
          enddo
       enddo

       ! -----
       ! Rinv(1:s,1:s) = Rtri0(1:s,1:s)
       do j=1,s
          do i=1,s
             Rinv(i+s*(j-1))=Rtri0(i+(j-1)*(s+1))
          enddo
       enddo
       call DTRTRI('U', 'N', s, Rinv, s, ierr)
#if 1
       ! -----
       ! 計算済みの上ヘッセンベルグ行列をコピー
       do j=1,s*(k-1)
          do i=1,s*(k-1)
             BlockTri(i+(j-1)*(s*(k-1))) = Tri_all(i + (j-1)*(s*t+1) )
          enddo
       enddo
       ! Hes[k-1] * BlockR
       call wrap_BASIC_DGEMM(s*(k-1),s,s*(k-1),BlockTri          ,BlockR,BlockTri_R)
       ! BlockR_ * B_
       call wrap_BASIC_DGEMM(s*(k-1),s,s+1    ,BlockR_ubar,B1    ,BlockR_ubar_B1)

       ! ------
       ! Hesk[k] = ( -Hes[k-1]*BlockR + BlockR_*B_ )* R^-1
       do i=1,s*(k-1) * s
          BlockTri_R(i)=-BlockTri_R(i)+BlockR_ubar_B1(i)
       enddo
       call wrap_BASIC_DGEMM(s*(k-1),s,s     ,BlockTri_R,Rinv,BlockR_ubar_B1)
       do i=1,s
          do j=1,s*(k-1)
             Tri_all( j + (i + s*(k-1)-1 ) * (s*t+1)  ) = BlockR_ubar_B1(j+(i-1)*s*(k-1))
          enddo
       enddo
#endif

       ! -------
       ! Hoeman 論文(4.11)
       bk=B1((s+1)*s)
       rho_tilde = Rtri0(s   + (s-1)*(s+1))
       rho       = Rtri0(s+1 +     s*(s+1))
       betak=rho/rho_tilde*bk
       Tri_all( s*k+1 +(s*k-1)*(s*t+1) ) =betak

       ! -------
       ! 生成した基底に対する三重対角行列の要素
       ! Hoeman 論文(4.10)

       ! Tri = Rk Bk Rinv
       call calc_Rk_Bk_Rkinv(s,Rtri0,B1,Rinv,Tri)
       ! rho~^-1 * bk * zk
       do i=1,s
          zk = Rtri0(i+s*(s+1) )
          Tri(i+s*(s-1)) = Tri(i+s*(s-1))+bk/rho_tilde*zk
       enddo
       ! betak * Rk-1,k Rinv
       call calc_betake1eskRRinv(s,k,betak_pre,BlockR,Rinv,Tri)

       ! -------
       do i=1,s
          do j=1,s
             Tri_all( j+ (i-1)*(s*t+1) + s*(k-1)+ (s*(k-1))*(s*t+1) ) = Tri(j+s*(i-1) )
          enddo
       enddo
       betak_pre=betak
       ! -----

#if 0
       Tri_all(s*(k-1) + (s*(k-1))*(s*t+1)) = &
               Tri_all(s*(k-1)+1 + (s*(k-1)-1)*(s*t+1))
#endif
    endif

#if 0
    do i = 1, s*t
      do j = 1, s*t+1
        write(*,*) "Tri@",j,i,Tri_all(j + (s*t+1)*(i-1))
      enddo
        write(*,*) "-----------"
    enddo
#endif

  end subroutine calc_TridiagonalMatrix

  ! Vk = Vk - Qk-1 R
  subroutine CALANCZOS_GEMM(n,s,ms,V,Q,R,nvecBlock)
    use vec_kernel
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) s
    integer(kind=intType) ms
    integer nvecBlock
    real*8 V(n*s)
    real*8 Q(n*ms)
    real*8 R(ms*s)
    integer ii
    integer i,j,k
    ! do k=1,s
    !    do j=1,ms
    !       write(*,*)"R = ",j,k,R(j+ms*(k-1))
    !    enddo
    ! enddo

    !$OMP PARALLEL DO PRIVATE(ii,k,j,i) 
    do ii=1,n,nvecBlock
       do k=1,s
          do j=1,ms
             do i=ii,min(ii+nvecBlock-1,n)
                V(i+n*(k-1)) = V(i+n*(k-1)) - Q(i+n*(j-1))*R(j+ms*(k-1))
             enddo
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO 

  end subroutine CALANCZOS_GEMM

  ! V(i)=V(i)+alpha*V(i-1)+beta*V(i-2)
  subroutine CALANCZOS_axpy3(n,alpha,beta,V,V1,V2)
    implicit none 
    integer n
    real*8 alpha,beta
    real*8 V(n)
    real*8 V1(n)
    real*8 V2(n)
    integer i

    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
       V(i)=V(i)+alpha*V1(i)+beta*V2(i)
    enddo

  end subroutine CALANCZOS_axpy3

  subroutine calc_NewtonBasis_lanVec(pcalanczos,n,s,i,pv,V,B1)
    use krylov
    use krylov_kernel
    use vec_kernel
    implicit none
    type(KSP) :: pcalanczos
    integer n
    integer s
    integer i
    real*8 pV(n)
    real*8 V(n*(s+1))
    real*8 B1((s+1)*s)
    real*8 alpha,beta
    call calc_Matvec( pcalanczos,pv,V(1+i*n) )
    if(i.ne.1)then
       alpha = -B1( i + (i-1)*(s+1) )
       beta  = B1( i-1 + (i-1-1+1)*(s+1) )
       call CALANCZOS_axpy3(n,alpha,beta,V(1+i*n),V(1+(i-1)*n),V(1+(i-1-1)*n))
    else
       alpha = -B1( i + (i-1)*(s+1) )
       call axpy(n,alpha,1.0d0,V(1+(i-1)*n),V(1+i*n))
    endif
  end subroutine calc_NewtonBasis_lanVec

  subroutine ConvergedCheck(n, st, Tri0,&
                            EigenReal, EigenImag, &
                            EigenVecReal, EigenVecImag, &
                            tol, iflag_exit)

    use vec_kernel
    implicit none

    ! arguments
    integer n
    integer st
    real*8 Tri0((st+1)*st)
    real*8 EigenReal(st)
    real*8 EigenImag(st)
    real*8 EigenVecReal(n*st)
    real*8 EigenVecImag(n*st)
    real*8 tol
    integer iflag_exit

    ! local variables
    integer i
    integer j
    real*8 ritz_norm(st)
    real*8 tmp

    call getTridiagonalMat_Eigen(&
         st, st ,&
         Tri0,&
         EigenReal,EigenImag)

    if (iflag_exit .eq. 1) return

    call getTridiagonalMat_EigenVector(&
         st, st ,&
         Tri0,&
         EigenReal,EigenImag, &
         n,EigenVecReal,EigenVecImag)

    do i = 1, st
      tmp = 0.0d0
      call calc_dot_local(n, EigenVecReal(1+n*(i-1)), &
                     EigenVecReal(1+n*(i-1)), tmp)
      tmp = sqrt(tmp)
      !$OMP PARALLEL DO
      do j = 1, n
        EigenVecReal(j+n*(i-1)) = EigenVecReal(j+n*(i-1))/tmp
      enddo
      !$OMP END PARALLEL DO
    enddo

    ritz_norm = 0.0d0
    do i = 1, st
      ritz_norm(i) = Tri0((st+1)*st)*sqrt((EigenVecReal(st+n*(i-1))**2))
      if (ritz_norm(i) < tol) iflag_exit = 1
    enddo

  end subroutine ConvergedCheck

end module KSP_CALANCZOS

! ------

subroutine KSP_pCALANCZOS(pcalanczos)
  use krylov
  use krylov_kernel
  use vec_kernel
  use QR_kernel
  use KSP_CAGMRES
  use KSP_CALANCZOS

  implicit none
  type(KSP) :: pcalanczos
  integer n
  integer icomm
  integer myrank ! プロセスランク
  integer s ! 省通信ステップ 
  integer t ! リスタート数 = (s*t)
  integer iQRflag ! QR分解の方法 
  integer basis ! 基底ベクトルの生成方法
  integer EIGflag ! 固有ベクトルを計算するかのフラグ
  integer EIG_calcflag ! 固有ベクトルの計算方法
  integer itrmax

  integer iQmatflag
  ! ----------
  integer k
  integer j
  integer i
  integer itr
  integer itr0
  integer msize
  integer iflag_exit

  real*8,allocatable :: r(:)

  real*8,allocatable :: pV(:)
  real*8,allocatable :: V(:)
  real*8,allocatable :: Q(:)

  real*8 rnorm
  real*8 tmp1,tmp2
  real*8,allocatable :: B1(:) 
  real*8,allocatable :: EigenReal(:)
  real*8,allocatable :: EigenImag(:)
  real*8,allocatable :: EigenVecReal(:)
  real*8,allocatable :: EigenVecImag(:)

  real*8,allocatable :: TridiagonalMat(:)
  real*8,allocatable :: Rtri(:) 
  real*8,allocatable :: BlockR_ubar_prime(:) 

  integer iflag_basis
  real*8 betak_pre
  real*8 alpha,beta
  real*8 tol
  integer nvecBlock

  myrank=pcalanczos%commVec%myrank
  icomm=pcalanczos%commVec%comm
  n=pcalanczos%prm%n
  s=pcalanczos%CALANCZOS_sstep
  t=pcalanczos%CALANCZOS_tstep
  iQRflag = pcalanczos%CALANCZOS_QRflag
  basis   = pcalanczos%CALANCZOS_basis
  EIGflag = pcalanczos%EIGflag
  EIG_calcflag = pcalanczos%EIG_calcflag
  itrmax = pcalanczos%itrmax
  tol = pcalanczos%rtolmax
  nvecBlock= pcalanczos%nvecBlock

  allocate(r(n))
  allocate(pV(n))
  allocate(V(n*(s+1)))
  allocate(Q(n*(s*t+1)))
  allocate(B1((s+1)*s))
  allocate(EigenReal(s*t))
  allocate(EigenImag(s*t))
  allocate(EigenVecReal(n*s*t))
  allocate(EigenVecImag(n*s*t))

  ! ------
  if(allocated(pcalanczos%EigenValueReal))then
     deallocate(pcalanczos%EigenValueReal)
  endif
  if(allocated(pcalanczos%EigenValueImag))then
     deallocate(pcalanczos%EigenValueImag)
  endif
  allocate(pcalanczos%EigenValueReal(s*t+1))
  allocate(pcalanczos%EigenValueImag(s*t+1))
  pcalanczos%EigenValueReal = 0.0d0
  pcalanczos%EigenValueImag = 0.0d0
  ! ------

  allocate(TridiagonalMat(  (s*t+1)*s*t ))
  allocate(Rtri((s+1)*(s+1)))
  allocate(BlockR_ubar_prime( (s*t+1)*(s+1) ) )

  call zero_initialize(n,r)
  call zero_initialize(n,pv)
  msize=n*(s+1)
  call zero_initialize(msize,V)
  msize=n*(s*t+1)
  call zero_initialize(msize,Q)

  B1=0.0d0
  Rtri=0.0d0

  ! ----------
  ! CALANCZOS のメインループ
  ! r = 1.0d0
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

  iflag_exit=-1
  iflag_basis=-1
  EigenReal=0.0d0
  EigenImag=0.0d0
  alpha=0.0d0
  beta =0.0d0

  ! リスタート付きCA_LANCZOS
  ! -----
  !$OMP PARALLEL DO PRIVATE(i)
  do i=1,n
    r(i)=1.0d0
  enddo

  do itr = 1, itrmax, s*t
    itr0 = itr
    V=0.0d0
    Q=0.0d0
    TridiagonalMat=0.0d0
    Rtri=0.0d0
    BlockR_ubar_prime=0.0d0
    betak_pre=0.0d0

    call calc_dot(n,icomm,r,r,rnorm)
    rnorm=sqrt(rnorm)

    tmp1=-1.0d0/rnorm
    tmp2=0.0d0
    ! v=r/||r||
    call axpy(n,tmp1,tmp2,r,V)

    do k = 1,t
      ! 基底変換行列を作成
      if( (basis.eq.KSP_MPK_NEWTON_BASIS).and.(iflag_basis.eq.1) )then
        call calc_NewtonBasis_lan(s,B1,EigenReal,EigenImag)
      endif
      if(k.eq.t)then
        iQmatflag=KSP_QRDECOMPOSITION_FALSE_ORTHOGONAL_MATRIX
      else
        iQmatflag=KSP_QRDECOMPOSITION_TRUE_ORTHOGONAL_MATRIX
      endif
      do i=1,s
        ! MpV=v_i
        call calc_precon(pcalanczos,V(1+(i-1)*n),pV )
        ! V_i+1=ApV
        if(iflag_basis.eq.1)then
          select case (basis)
          case(KSP_MPK_NEWTON_BASIS)
            call calc_NewtonBasis_lanVec(pcalanczos,n,s,i,pv,V,B1)
          case(KSP_MPK_MONOMIAL_BASIS)
            call calc_Matvec( pcalanczos,pv,V(1+i*n) )
          end select
        else
          ! 初回反復は固有値の情報がないため単基底
          call calc_Matvec( pcalanczos,pv,V(1+i*n) )
        endif
      enddo

      if(k.eq.1)then
        ! ------
        ! QR分解
        call QRDecomposition(n,s+1,icomm,V,Q,Rtri,iQRflag,iQmatflag,nvecBlock)
      else
        ! ------
        !  ブロック古典グラムシュミット直交化
        ! ------
        !  R_' = Q_k-1 * Vk
        call calc_gemm_dot(n,s*(k-1)+1,s+1,icomm,Q,V,BlockR_ubar_prime,nvecBlock)
        !  Vk = Vk - QR_'
        call CALANCZOS_GEMM(n,s,s*(k-1)+1, &
             V(1+n),Q,BlockR_ubar_prime(1+s*(k-1)+1) ,nvecBlock)
        ! QR分解
        call QRDecomposition(n,s,icomm,V(1+n),Q(1+n*(s*(k-1)+1)),Rtri,iQRflag,iQmatflag,nvecBlock)

      endif

      ! ------
      ! ヘッセンベルグ行列の計算
      call calc_TridiagonalMatrix( &
           s,t,k, &
           B1, &
           TridiagonalMat, &
           Rtri, &
           BlockR_ubar_prime, &
           betak_pre &
           )

      ! -------
      if(k.ne.t)then
        ! 次の基底生成の初期ベクトル
        call xcpy(n,Q( 1 + n*( s*(k-1)+1+s-1 ) ),V)
      endif

      if ((itr0-1)+s*k .ge. itrmax) then
        iflag_exit = 1
        itr0 = (itr0-1) + s*k
        !exit
      endif


      if(basis.eq.KSP_MPK_NEWTON_BASIS)then
        iflag_basis=1
        ! 固有値計算
        call getTridiagonalMat_Eigen(&
             s*k, s*t ,&
             TridiagonalMat,&
             EigenReal,EigenImag)
        call modified_leja_ordering(s*k,EigenReal,EigenImag)
      endif

    enddo

    call ConvergedCheck(n, s*t, TridiagonalMat, &
                        EigenReal, EigenImag, &
                        EigenVecReal, EigenVecImag, &
                        tol, iflag_exit)

    itr0 = (itr - 1)/(s*t) + s*t
    if (iflag_exit.eq.1) exit

    tmp1 = TridiagonalMat((s*t+1)*s*t)
    !$OMP PARALLEL DO PRIVATE(i)
    do i = 1, n
!      r(i) = Q(i+n*(s*t+1))/tmp1
      r(i) = Q(i+n*(s*(t-1)))/tmp1
    enddo
    !$OMP END PARALLEL DO

  enddo

  call getRealEigen_lan(&
       s*t,&
       EigenReal,EigenImag,&
       pcalanczos%EigenMax_Real, &
       pcalanczos%EigenMax_Imag, &
       pcalanczos%EigenMin_Real, &
       pcalanczos%EigenMin_Imag  &
       )

   do i = 1, s*t
     pcalanczos%EigenValueReal(i) = EigenReal(i)
     pcalanczos%EigenValueImag(i) = EigenImag(i)
   enddo

  ! calculation eigen vector
  if ( EIGflag .eq. 1 ) then
    call getTridiagonalMat_EigenVector(&
       s*t, s*t ,&
       TridiagonalMat,&
       EigenReal,EigenImag, &
       n,EigenVecReal,EigenVecImag)
#if 0
    if (EIG_calcflag .eq. 0) then
      pcalanczos%eigvecreal = 0.0d0
      pcalanczos%eigvecimg = 0.0d0
      do i = 1, s*t
        call sum_vec(n, s*t, Q, pcalanczos%eigvecreal(1, i), EigenVecReal(1+n*(i-1)))
        call sum_vec(n, s*t, Q, pcalanczos%eigvecimg(1, i), EigenVecImag(1+n*(i-1)))
      enddo
    else if (EIG_calcflag .eq. 1) then
       call KSP_pCALANCZOS_II(pcalanczos, itr0, n, s, t, &
                              EigenVecReal, EigenVecImag)
    endif
#endif
  endif

#if 1
  if(myrank.eq.0)then
     ! do j=1,s*t
     !    do i=1,s*t
     !       write(*,*)"Tri     @ ",i,j,TridiagonalMat(i+(j-1)*(s*t+1))
     !    enddo
     !    write(*,*)"------"
     ! enddo

     do i=1,s*t
        write(*,*)"Eigen ",i,EigenReal(i),EigenImag(i)
     enddo
     write(*,*)"EigenMax ",pcalanczos%EigenMax_Real,pcalanczos%EigenMax_Imag
     write(*,*)"EigenMin ",pcalanczos%EigenMin_Real,pcalanczos%EigenMin_Imag
     write(*,*)"itr0, itr, itrmax ", itr0, itr, itrmax
  end if
#endif

  deallocate(r)
  deallocate(pV)
  deallocate(V)
  deallocate(Q)
  deallocate(B1)
  deallocate(EigenReal)
  deallocate(EigenImag)
  deallocate(EigenVecReal)
  deallocate(EigenVecImag)

  deallocate(TridiagonalMat)
  deallocate(Rtri)
  deallocate(BlockR_ubar_prime)

end subroutine KSP_pCALANCZOS

