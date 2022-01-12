! subroutine dd_mul(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo)
! subroutine dd_div(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo)
! subroutine dd_sub(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo)

#include "ddfunc.h"
module ddmatrix_kernel
  use dd_func
  use PARCEL_define
  use common_kernel
  implicit none
contains
  subroutine ddget_nnz_dia(diaA_hi,diaA_lo,offset,n,istart,nnd,nnz)
    use matrix
    use PARCEL_integer
    implicit none
    type(DIA_Mat) :: diaMat 
    integer(kind=intType) n
    integer(kind=intType) istart

    integer(kind=intType) nnd
    integer(kind=intType) nnz
    integer(kind=intType) offset(nnd)
    real*8 diaA_hi(nnd*n)
    real*8 diaA_lo(nnd*n)

    integer(kind=intType) j,k
    integer(kind=intType) Row,Col
    real*8 absAmat_hi
    real*8 absAmat_lo

    nnz=0
    do k=1,n
       Row=k-1+istart
       do j=1,nnd
          Col=Row+offset(j)
          absAmat_hi=abs(diaA_hi(k+(j-1)*n))
          absAmat_lo=abs(diaA_lo(k+(j-1)*n))
          if( (absAmat_hi > minReal8).or.(absAmat_lo > minReal8) )then
             nnz=nnz+1
             !             write(*,*)"  DIA  ",Row,Col,diaA(k+(j-1)*n)
          endif
       enddo
    enddo
  end subroutine ddget_nnz_dia

  subroutine dddia2crs(diaMat,crsMat,n,istart,nnd,nnz)
    use matrix
    use vec_kernel
    use PARCEL_integer
    implicit none

    type(CRS_Mat) :: crsMat 
    type(DIA_Mat) :: diaMat 

    integer(kind=intType) nnd
    integer(kind=intType) nnz
    integer(kind=intType) istart
    integer(kind=intType) n

    integer(kind=intType) i,j,k
    integer(kind=intType) offset
    integer(kind=intType) offset_sort(nnz)
    integer(kind=intType) Row,Col
    real*8 absAmat_hi
    real*8 absAmat_lo

    i=1
    do k=1,n
       Row=k-1+istart
       crsMat%crsRow_ptr(k) = i
       do j=1,nnd
          Col=Row+diaMat%offset(j)
#if 1
          absAmat_hi=abs(diaMat%diaA_hi(k+(j-1)*n))
          absAmat_lo=abs(diaMat%diaA_lo(k+(j-1)*n))
          if((absAmat_hi > minReal8).or.(absAmat_lo > minReal8))then
             crsMat%crsA_hi(i)=diaMat%diaA_hi(k+(j-1)*n)
             crsMat%crsA_lo(i)=diaMat%diaA_lo(k+(j-1)*n)
             crsMat%crsCol(i) = Col
             i=i+1
          endif
#else
          crsMat%crsA_hi(i)=diaMat%diaA_hi(k+(j-1)*n)
          crsMat%crsA_lo(i)=diaMat%diaA_lo(k+(j-1)*n)
          crsMat%crsCol(i) = Col
          i=i+1
#endif
       enddo
    enddo
    crsMat%crsRow_ptr(n+1) = i
  end subroutine dddia2crs

  subroutine init_dddiaMat(diaMat,n,nnd)
    use matrix
    use PARCEL_integer
    implicit none
    type(DIA_Mat) :: diaMat 
    integer(kind=intType) n
    integer(kind=intType) nnd
    integer(kind=intType) i
    allocate( diaMat%diaA_hi(n*nnd) )
    allocate( diaMat%diaA_lo(n*nnd) )
    allocate( diaMat%offset(nnd) )
    allocate( diaMat%offsetVec(nnd) )
    allocate( diaMat%Vec_hi(1) )
    allocate( diaMat%Vec_lo(1) )

    allocate(diaMat%RowL(1))
    allocate(diaMat%RowD(1))
    allocate(diaMat%RowU(1))

    !$OMP PARALLEL DO PRIVATE(i)
    do i=1,n*nnd
       diaMat%diaA_hi(i)=0.0d0
       diaMat%diaA_lo(i)=0.0d0
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i)
    do i=1,nnd
       diaMat%offset(i)=0
       diaMat%offsetVec(i)=0
    enddo
    !$OMP END PARALLEL DO

  end subroutine init_dddiaMat

  ! DIAの疎行列ベクトル積のベクトル位置とベクトルサイズを計算
  subroutine set_dddia_Row(prm,diaMat)
    use matrix
    use matrix_kernel
    use vec_kernel
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(DIA_Mat) :: diaMat 
    integer(kind=intType) RowL(prm%nnd*prm%n)
    integer(kind=intType) RowU(prm%nnd*prm%n)

    integer(kind=intType) n,istart
    integer(kind=intType) nCol,nColL,nColD,nColU
    ! -----
    n = prm%n
    istart=prm%istart

    call get_dia_Row(&
         prm,diaMat,&
         RowL,RowU,&
         n,nColL,nColD,nColU,nCol)    
    ! ベクトル全体の長さ 
    ! nCol=nColL+nColD+nColU
    nCol=nColL + nColD + nColU
    diaMat%nColL=nColL
    diaMat%nColD=nColD
    diaMat%nColU=nColU
    diaMat%nCol =nCol
    call ddreplace_dia_vec(&
         diaMat,&
         istart,&
         RowL,RowU,&
         n,nColL,nColD,nColU,nCol)    
  end subroutine set_dddia_Row

  subroutine ddreplace_dia_vec(&
       diaMat,&
       istart,&
       RowL,RowU,&
       n,nColL,nColD,nColU,nCol)
    use matrix
    use vec_kernel
    use PARCEL_integer
    implicit none
    integer(kind=intType) i
    type(DIA_Mat) :: diaMat 
    integer(kind=intType) istart
    integer(kind=intType) n,nColL,nColD,nColU,nCol

    integer(kind=intType) RowL(nColL)
    integer(kind=intType) RowU(nColU)

    deallocate(diaMat%RowL)
    deallocate(diaMat%RowD)
    deallocate(diaMat%RowU)
    deallocate(diaMat%Vec_hi)
    deallocate(diaMat%Vec_lo)

    allocate(diaMat%RowL(nColL))
    allocate(diaMat%RowD(nColD))
    allocate(diaMat%RowU(nColU))
    allocate(diaMat%Vec_hi(nCol))
    allocate(diaMat%Vec_lo(nCol))

    do i=1,nColL
       diaMat%RowL(i)=RowL(i)
    enddo
    do i=1,n
       diaMat%RowD(i)=i+istart-1
    enddo
    do i=1,nColU
       diaMat%RowU(i)=RowU(i)
    enddo
    call zero_initialize(nCol,diaMat%Vec_hi)
    call zero_initialize(nCol,diaMat%Vec_lo)

  end subroutine ddreplace_dia_vec

  subroutine ddget_nnz_Colsort_dia( &
       prm,diaMat,Colsort,&
       nCol,nColL,nColU,nColD)
    use matrix
    use sort_kernel
    use PARCEL_integer
    implicit none
    type(prmMat) :: prm
    type(DIA_Mat) :: diaMat 
    integer(kind=intType) Colsort(prm%nnz)
    integer(kind=intType) nCol
    integer(kind=intType) nColL,nColU,nColD

    integer(kind=intType) n
    integer(kind=intType) istart,iend
    integer(kind=intType) nnd
    integer(kind=intType) nnz
    integer(kind=intType) i,j,k
    integer(kind=intType) Row,Col

    real*8 absAmat_hi
    real*8 absAmat_lo

    n = prm%n
    istart = prm%istart
    iend = prm%iend
    nnd = prm%nnd
    nnz = prm%nnz
    i=0
    do k=1,n
       Row=k-1+istart
       do j=1,nnd
          absAmat_hi=abs(diaMat%diaA_hi(k+(j-1)*n))
          absAmat_lo=abs(diaMat%diaA_lo(k+(j-1)*n))
          if((absAmat_hi > minReal8).or.(absAmat_lo > minReal8))then
             i=i+1
             Colsort(i)=Row+diaMat%offset(j)
          endif
       enddo
    enddo
    !    call sort_c(Colsort,nnz)
    call MergeSort(Colsort,nnz)

    nCol=1
    do i=1,nnz
       Col=Colsort(i)
       if(Colsort(nCol).ne.Col )then
          nCol=nCol+1
          Colsort(nCol)=Col
       endif
    enddo

    if(nnz.eq.0)then
       nCol=0
    endif

    nColL=0
    nColD=0
    nColU=0
    !$OMP PARALLEL DO PRIVATE(i,Col) REDUCTION(+:nColL,nColD,nColU)
    do i=1,nCol
       Col=Colsort(i)
       if(Col<istart)then
          nColL=nColL+1
       endif
       if( (istart<=Col) .and. (Col<=iend) )then
          nColD=nColD+1
       end if
       if(iend<Col)then
          nColU=nColU+1
       endif
    enddo
    !$OMP END PARALLEL DO

  end subroutine ddget_nnz_Colsort_dia

  subroutine ddget_DInv_crs(prm,crsMat,Dinv_hi,Dinv_lo)
    use matrix
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(CRS_Mat) :: crsMat 
    real*8 Dinv_hi(prm%n)
    real*8 Dinv_lo(prm%n)
    real*8 dddiv_var1

    integer(kind=intType) i,j,k
    integer(kind=intType) n,istart,iend
    integer(kind=intType) Row,Col
    integer(kind=intType) nColL,nColU

    n = prm%n
    istart=prm%istart
    iend=prm%iend
    nColL=crsMat%nColL
    nColU=crsMat%nColU

    !$OMP PARALLEL DO PRIVATE(j,k,Row,Col,&
    !$OMP dddiv_var1) 
    do k=1,n
       Row=k-1+istart
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          if( crsMat%crsCol(j)<istart ) then
             cycle
          endif
          if( iend<crsMat%crsCol(j) )then
             cycle
          endif
          Col=crsMat%crsCol(j)-istart+1
          if(k.eq.Col)then
             !             Dinv(k)=1.0d0/crsMat%crsA(j)
#ifdef Manual_inline
             if( (crsMat%crsA_hi(j).eq.0.0d0).and.&
                  (crsMat%crsA_lo(j).eq.0.0d0))then
                dddiv_calc1(1.0d0,0.0d0,1.0d-100,0.0d0,Dinv_hi(k),Dinv_lo(k))
             else
                dddiv_calc1(1.0d0,0.0d0,crsMat%crsA_hi(j),crsMat%crsA_lo(j),Dinv_hi(k),Dinv_lo(k))
             endif
#else
             if( (crsMat%crsA_hi(j).eq.0.0d0).and.&
                 (crsMat%crsA_lo(j).eq.0.0d0))then
                call dd_div( &
                     1.0d0, &
                     0.0d0, &
                     1.0d-100, &
                     0.0d0 , &
                     Dinv_hi(k), &
                     Dinv_lo(k) )
             else
                call dd_div( &
                     1.0d0, &
                     0.0d0, &
                     crsMat%crsA_hi(j) , &
                     crsMat%crsA_lo(j) , &
                     Dinv_hi(k), &
                     Dinv_lo(k) )
             endif
#endif
          endif
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine ddget_DInv_crs

  subroutine ddget_DInv_dia(prm,diaMat,Dinv_hi,Dinv_lo)
    use matrix
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(DIA_Mat) :: diaMat 
    real*8 Dinv_hi(prm%n),Dinv_lo(prm%n)
    real*8 dddiv_var1
    integer(kind=intType) i,j,k
    integer(kind=intType) n,nnd

    n = prm%n
    nnd = prm%nnd

    do i=1,nnd
       if(diaMat%offset(i).eq.0)then
          !$OMP PARALLEL DO PRIVATE(k,dddiv_var1) 
          do k=1,n
             !             Dinv(k)=1.0d0/diaMat%diaA(k+(i-1)*n)
#ifdef Manual_inline
             if( (diaMat%diaA_hi(k+(i-1)*n).eq.0.0d0).and.&
                 (diaMat%diaA_lo(k+(i-1)*n).eq.0.0d0))then
                dddiv_calc1(1.0d0,0.0d0,1.0d-100,0.0d0,Dinv_hi(k),Dinv_lo(k))
             else
                dddiv_calc1(1.0d0,0.0d0,diaMat%diaA_hi(k+(i-1)*n),diaMat%diaA_lo(k+(i-1)*n),Dinv_hi(k),Dinv_lo(k))
             endif
#else
             if( (diaMat%diaA_hi(k+(i-1)*n).eq.0.0d0).and.&
                 (diaMat%diaA_lo(k+(i-1)*n).eq.0.0d0))then
                call dd_div( &
                     1.0d0, &
                     0.0d0, &
                     1.0d-100, &
                     0.0d0 , &
                     Dinv_hi(k), &
                     Dinv_lo(k) )
             else
                call dd_div( &
                     1.0d0, &
                     0.0d0, &
                     diaMat%diaA_hi(k+(i-1)*n), &
                     diaMat%diaA_lo(k+(i-1)*n), &
                     Dinv_hi(k), &
                     Dinv_lo(k) )
             endif
#endif
          enddo
          !$OMP END PARALLEL DO
       endif
    enddo
  end subroutine ddget_DInv_dia

  subroutine ddget_DILU_DInv_dia( &
       prm, &
       diaMat, &
       Dinv_hi, &
       Dinv_lo, &
       ILU_method, &
       independ_nvec, &
       precon_thblock)
    use omp_lib
    use matrix
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(DIA_Mat) :: diaMat 
    real*8 Dinv_hi(prm%n)
    real*8 Dinv_lo(prm%n)
    integer(kind=intType) n,nnd,istart,iend
    integer numth,myth,max_numth
    integer(kind=intType) n_th
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th
    integer ILU_method
    integer independ_nvec
    integer precon_thblock

    integer(kind=intType) ntmp

    n = prm%n
    nnd = prm%nnd
    istart=prm%istart
    iend=prm%iend

    ! -----
    ! ブロック数よりもスレッド数が大きい場合はスレッド数を制限する
    call get_maxthreads(max_numth)
    call set_precon_maxthreads(n,precon_thblock)

    ! -----
    ! istart_th=istart
    ! iend_th=iend
    ! n_startth=1
    ! n_endth=n

    !$OMP PARALLEL &
    !$OMP PRIVATE( &
    !$OMP numth,myth,n_startth,n_endth,istart_th,iend_th, ntmp &
    !$OMP )
    numth = omp_get_num_threads()
    myth  = omp_get_thread_num()
    ntmp=1
    call divide_region( &
         myth,numth,ntmp,n,n_startth,n_endth)
    istart_th = istart -1 + n_startth 
    iend_th   = istart_th + (n_endth - n_startth) 
    call ddDILU_DInv_dia( &
         diaMat,n,nnd, &
         Dinv_hi(n_startth), &
         Dinv_lo(n_startth), &
         istart,iend, &
         n_startth,n_endth, &
         istart_th,iend_th, &
         ILU_method, &
         independ_nvec)
    !$OMP END PARALLEL

!    if(n<=max_numth)then
    call omp_set_num_threads(max_numth)
!    endif
  end subroutine ddget_DILU_DInv_dia

  subroutine ddDILU_DInv_dia( &
       diaMat,n,nnd, &
       Dinv_hi, &
       Dinv_lo, &
       istart,iend, &
       n_startth,n_endth, &
       istart_th,iend_th, &
       ILU_method, &
       independ_nvec)
    use omp_lib
    use matrix
    use PARCEL_integer
    implicit none
    type(DIA_Mat) :: diaMat 
    integer(kind=intType) istart,iend
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th
    real*8 Dinv_hi(n_startth:n_endth)
    real*8 Dinv_lo(n_startth:n_endth)
    integer ILU_method
    integer independ_nvec
    integer(kind=intType) n,nnd,nnd_offset0

    ! real*8 D_hi(n_startth:n_endth)
    ! real*8 D_lo(n_startth:n_endth)
    real*8,allocatable::D_hi(:)
    real*8,allocatable::D_lo(:)

    integer(kind=intType) i,j,k,kk
    integer(kind=intType) Row,Col,gCol

    integer(kind=intType) kr
    integer(kind=intType) jr,rCol,rRow
    integer(kind=intType) jstart
    integer(kind=intType) grCol
    real*8 ddmul_var1
    real*8 dddiv_var1
    real*8 ddsub_var1
    real*8 c_hi,c_lo
    integer independ_start,independ_end

    ! -----
    ! D-ILU分解
    ! A = AL + AD + AU ~= K
    ! K = (D+AL)D^-1(D+AU)
    ! 条件1 対角成分が一致     :: D(i,i) = A(i,i) - Sum[j<i]          A(i,j) D(j,j)^-1 A(j,i)
    ! 条件2 行内の要素和が一致 :: D(i,i) = A(i,i) - Sum[j<i] Sum[k>j] A(i,j) D(j,j)^-1 A(j,k)
    ! -----

    allocate(D_hi(n_startth:n_endth))
    allocate(D_lo(n_startth:n_endth))
    nnd_offset0=-1
    do j=1,nnd
       if(diaMat%offset(j)==0) then
          nnd_offset0=j
       endif
    enddo
    
    if (nnd_offset0==-1) then
       write(*,*)"Error:: DILU--There is no diagonal element."
       stop
    endif
    
    kk=-1
    do k=n_startth,n_endth
!       if(diaMat%diaA(k+(nnd_offset0-1)*n)==0.0) then
       if( (diaMat%diaA_hi(k+(nnd_offset0-1)*n)==0.0d0) .and. (diaMat%diaA_lo(k+(nnd_offset0-1)*n)==0.0d0) ) then
         write(*,*)"Error:: DILU--It's lack of diagonal element."
         stop
       endif
    enddo
    

    do k=n_startth,n_endth
       ! 独立に回すループサイズ
       independ_start = n_startth + ((k-n_startth)/independ_nvec) * independ_nvec
       independ_end   = min(independ_start+independ_nvec-1,n_endth)

       Row=k-n_startth+istart_th
       jstart=-1
       D_hi(k)=diaMat%diaA_hi(k+(nnd_offset0-1)*n)
       D_lo(k)=diaMat%diaA_lo(k+(nnd_offset0-1)*n)
       do j=1,nnd
          gCol = k-1+istart + diaMat%offset(j)
          if( gCol<istart_th )then
             cycle
          endif
          if( iend_th<gCol )then
             cycle
          endif
          Col=k + diaMat%offset(j)  
             ! 対角ブロックの開始位置
             jstart=j
             ! -----
             ! 行列の対称位置に対して値があれば総和を計算する
             ! 下側部分を探査する


             ! -----
             ! 下側列の探査
             do jr=jstart,nnd_offset0-1
                Col=diaMat%offset(jr)+k
                kr=Col
                !                write(*,*)"Lower -- > Row,Col = ",kr,Col

                if( (independ_start<=Col).and.(Col<=independ_end) )then
                   if(Col.ne.k)then                
                      cycle
                   endif
                endif

                ! 上側行(kr行)の探査
                do i=nnd_offset0+1,nnd
                   grCol = kr-1+istart+diaMat%offset(i)
                   !if( grCol<istart_th ) then
                   !   cycle
                   !endif
                   if( iend_th<grCol )then
                      !                      cycle
                      exit
                   endif
                   rCol=kr+diaMat%offset(i)

                   if(ILU_method.eq.KSP_PRECONDITION_ILUD)then
                      ! 対角成分が一致
                      ! 行列転置位置の非ゼロ要素
                      if(k.eq.rCol)then
                         ! write(*,*)k," -- Row,Col = ",k,Col," -- ",kr,rCol
                         ! 行と列の転置位置を掛ける
                         ! kr行 k列
                         !  k行kr列
                         ! D(k)=D(k) - &
                         !      diaMat%diaA( kr + (i-1)  * n ) * &
                         !      diaMat%diaA( k  + (jr-1) * n )   / D(kr)
#ifdef Manual_inline
                         ddmul_calc1(diaMat%diaA_hi(kr+(i-1)*n),diaMat%diaA_lo(kr+(i-1)*n),diaMat%diaA_hi(k+(jr-1)*n),diaMat%diaA_lo(k+(jr-1)*n),c_hi,c_lo)
                         if( (D_hi(kr).eq.0.0d0).and.&
                              (D_lo(kr).eq.0.0d0))then
                            dddiv_calc1(c_hi,c_lo,1.0d-100,0.0d0,c_hi,c_lo)
                         else                            
                            dddiv_calc1(c_hi,c_lo,D_hi(kr),D_lo(kr),c_hi,c_lo)
                         endif
                         ddsub_calc1(D_hi(k),D_lo(k),c_hi,c_lo,D_hi(k),D_lo(k))
                         exit
#else
                         call dd_mul( &
                              diaMat%diaA_hi( kr + (i-1)  * n ), &
                              diaMat%diaA_lo( kr + (i-1)  * n ), &
                              diaMat%diaA_hi( k  + (jr-1) * n ), &
                              diaMat%diaA_lo( k  + (jr-1) * n ), &
                              c_hi,c_lo)
                         if( (D_hi(kr).eq.0.0d0).and.&
                              (D_lo(kr).eq.0.0d0))then
                            call dd_div( &
                                 c_hi,c_lo,&
                                 1.0d-100,0.0d0,&
                                 c_hi,c_lo)
                         else                            
                            call dd_div( &
                                 c_hi,c_lo,&
                                 D_hi(kr),D_lo(kr),&
                                 c_hi,c_lo)
                         endif
                         call dd_sub( &
                              D_hi(k),D_lo(k),c_hi,c_lo,D_hi(k),D_lo(k))
                         exit
#endif
                      endif
                   else
                      ! 行内の要素和が一致
                      if(kr<rCol)then
                            ! D(k)=D(k) - &
                            !      diaMat%diaA( kr + (i-1)  * n ) * &
                            !      diaMat%diaA( k  + (jr-1) * n )   / D(kr)
#ifdef Manual_inline
                         ddmul_calc1(diaMat%diaA_hi(kr+(i-1)*n),diaMat%diaA_lo(kr+(i-1)*n),diaMat%diaA_hi(k+(jr-1)*n),diaMat%diaA_lo(k+(jr-1)*n),c_hi,c_lo)
                         if( (D_hi(kr).eq.0.0d0).and.&
                              (D_lo(kr).eq.0.0d0))then
                            dddiv_calc1(c_hi,c_lo,1.0d-100,0.0d0,c_hi,c_lo)
                         else
                            dddiv_calc1(c_hi,c_lo,D_hi(kr),D_lo(kr),c_hi,c_lo)
                         endif
                         ddsub_calc1(D_hi(k),D_lo(k),c_hi,c_lo,D_hi(k),D_lo(k))
#else
                         call dd_mul( &
                              diaMat%diaA_hi( kr + (i-1)  * n ), &
                              diaMat%diaA_lo( kr + (i-1)  * n ), &
                              diaMat%diaA_hi( k  + (jr-1) * n ), &
                              diaMat%diaA_lo( k  + (jr-1) * n ), &
                              c_hi,c_lo)
                         if( (D_hi(kr).eq.0.0d0).and.&
                              (D_lo(kr).eq.0.0d0))then
                            call dd_div(c_hi,c_lo,1.0d-100,0.0d0,c_hi,c_lo)
                         else
                            call dd_div(c_hi,c_lo,D_hi(kr),D_lo(kr),c_hi,c_lo)
                         endif
                         call dd_sub(D_hi(k),D_lo(k),c_hi,c_lo,D_hi(k),D_lo(k))
#endif
                      endif
                   endif

                enddo
             enddo
             exit
       enddo
    enddo

    do k=n_startth,n_endth
!      Dinv(k)=1.0d0/D(k)
#ifdef Manual_inline
       if( (D_hi(k).eq.0.0d0).and.&
            (D_lo(k).eq.0.0d0))then
          dddiv_calc1(1.0d0,0.0d0,1.0d-100,0.0d0,Dinv_hi(k),Dinv_lo(k))
       else
          dddiv_calc1(1.0d0,0.0d0,D_hi(k),D_lo(k),Dinv_hi(k),Dinv_lo(k))
       endif
#else
       if( (D_hi(k).eq.0.0d0).and.&
            (D_lo(k).eq.0.0d0))then
          call dd_div(1.0d0,0.0d0,1.0d-100,0.0d0,Dinv_hi(k),Dinv_lo(k))
       else
          call dd_div(1.0d0,0.0d0,D_hi(k),D_lo(k),Dinv_hi(k),Dinv_lo(k))
       endif
#endif
!       write(*,*)"Dinv -- ",k,Dinv_hi(k),Dinv_lo(k)
    enddo
    deallocate(D_hi)
    deallocate(D_lo)
    ! stop

  end subroutine ddDILU_DInv_dia

  subroutine ddget_DILU_DInv_AS_crs(method,crsMat)
    use mpi
    use krylov
    use vec_kernel
    use matrix_kernel,only:get_AS_region,DILU_DInv_AS_crs,check_AS_region
    use PARCEL_integer
    implicit none
    type(KSP) :: method
    type(CRS_Mat) :: crsMat
    integer icomm
    integer myrank
    integer npes
    integer(kind=intType) n
    integer(kind=intType) gn
    integer addL
    integer(kind=intType) istart,iend

    integer(kind=intType) istartL,iendL
    integer(kind=intType) istartU,iendU

    integer(kind=intType) istartAS,iendAS

    integer ierr
    integer s_Ureq,s_Lreq
    integer r_Ureq,r_Lreq
    integer s_Usta(MPI_STATUS_SIZE),s_Lsta(MPI_STATUS_SIZE)
    integer r_Usta(MPI_STATUS_SIZE),r_Lsta(MPI_STATUS_SIZE)

    integer Urank,Lrank
    integer s_nnzL,s_nnzU
    integer r_nnzL,r_nnzU
    integer(kind=intType) kk,k,j,i
    integer iflagl,iflagg
    integer(kind=intType) Col
    integer(kind=intType) i0

    real*8,allocatable  :: sU_crsA_hi(:)
    real*8,allocatable  :: sU_crsA_lo(:)
    integer(kind=intType),allocatable :: sU_crsRow_ptr(:) 
    integer(kind=intType),allocatable :: sU_crsCol(:) 

    real*8,allocatable  :: sL_crsA_hi(:)
    real*8,allocatable  :: sL_crsA_lo(:)
    integer(kind=intType),allocatable :: sL_crsRow_ptr(:) 
    integer(kind=intType),allocatable :: sL_crsCol(:) 

    real*8 s_DinvL_hi(method%addL),s_DinvU_hi(method%addL)
    real*8 s_DinvL_lo(method%addL),s_DinvU_lo(method%addL)
    real*8 r_DinvL_hi(method%addL),r_DinvU_hi(method%addL)
    real*8 r_DinvL_lo(method%addL),r_DinvU_lo(method%addL)

    type(prmMat)     :: prm_ASU
    type(prmMat)     :: prm_ASL
    type(CRS_Mat)    :: precon_crsMatASU
    type(CRS_Mat)    :: precon_crsMatASL

    integer addL1
    integer(kind=intType) nnzAS
    integer(kind=intType) n_th
    integer max_numth

    integer(kind=intType),allocatable :: istart_th_AS(:)
    integer(kind=intType),allocatable :: iend_th_AS(:)

    n      = method%prm%n
    gn     = method%prm%gn
    icomm   = method%commVec%comm
    myrank = method%commVec%myrank
    npes   = method%commVec%npes
    addL   = method%addL
    istart = method%prm%istart
    iend   = method%prm%iend

    ! ----
    ! U(上側),L(下側)はベクトルを基準としている
    Urank=myrank-1
    Lrank=myrank+1
    if(myrank.eq.(npes-1))then
       Lrank=-1
    endif
    if(myrank.eq.0)then
       Urank=-1
    endif
    method%Lrank_AS=Lrank
    method%Urank_AS=Urank

    ! ----
    ! オーバーラップ領域が零の場合または1プロセスで実行する場合は
    ! ブロックヤコビ前処理
    !    if( (addL.eq.0).or.(npes.eq.1) )then
    if( (addL.eq.0) )then
       method%ipreflag=2
       call ddget_DILU_DInv_crs(&
            method%prm,crsMat,&
            method%Dinv_hi,method%Dinv_lo,&
            method%ILU_method, &
            method%independ_nvec, &
            method%precon_thblock)
       return 
    endif

    ! ----
    ! 重なり領域がプロセス受け持ち分のベクトルの長さを超える場合の処理
    iflagl=1
    iflagg=-1
    if(n<addL)then
       iflagl=-1
    endif
    ! ----
    ! 重なり領域がスレッド受け持ち分のベクトルの長さを超える場合の処理
    call check_AS_region( &
         method%prm,addL, &
         iflagl, &
         method%precon_thblock)

    call MPI_Allreduce(iflagl,iflagg,1,MPI_INTEGER4,MPI_MIN,icomm,ierr)
    if(iflagg<1)then
       call MPI_Barrier(icomm,ierr)
       if(myrank.eq.0)then
          write(*,*)"warning :: Additive Schwarz preconditioners"
          write(*,*)"           overlap parameter too large "
          write(*,*)"           change Additive Schwarz to Block Jacobi  "
       endif
       call MPI_Barrier(icomm,ierr)
       method%addL=0
       method%ipreflag=KSP_PRECONDITION_BLOCKJACOBI
       call ddget_DILU_DInv_crs(&
            method%prm,crsMat,&
            method%Dinv_hi,method%Dinv_lo,&
            method%ILU_method,&
            method%independ_nvec, &
            method%precon_thblock)

       return
    endif

    ! -------
    ! 加法シュワルツ

    ! ----
    ! 送信する行列の非零要素数を計算
    s_nnzU=0
    s_nnzL=0
    if(myrank.ne.(npes-1))then
       istartL=method%commVec%istart(myrank+1 +1)     - addL
       iendL = method%commVec%istart(myrank+1 +2)-1   + addL
       if(iendL>gn)then
          iendL=gn
       endif
       do k=n-addL+1,n
          do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
             Col=crsMat%crsCol(j)
             if( (istartL<=Col).and.(Col<=iendL) )then
                s_nnzL=s_nnzL+1
             endif
          enddo
       enddo
    else
       istartL=  0
       iendL  = -1
    endif
    if(myrank.ne.0)then
       istartU=method%commVec%istart(myrank+1 -1)   - addL
       iendU = method%commVec%istart(myrank+1   )-1 + addL
       if(istartU<=0)then
          istartU=1
       endif
       do k=1,addL
          do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
             Col=crsMat%crsCol(j)
             if( (istartU<=Col).and.(Col<=iendU) )then
                s_nnzU=s_nnzU+1
             endif
          enddo
       enddo

    else
       istartU=  0
       iendU  = -1
    endif

    ! -----
    ! 隣接プロセスの行列要素数を受信
    ! r_nnzU :: 上追加行の非零要素　
    ! r_nnzL :: 下追加行の非零要素　
    r_nnzU=0
    r_nnzL=0
    if(myrank.ne.(npes-1))then
       call MPI_Irecv(r_nnzL,1,MPI_INTEGER4,Lrank,Lrank ,icomm,r_Lreq,ierr)
       call MPI_Isend(s_nnzL,1,MPI_INTEGER4,Lrank,myrank,icomm,s_Lreq,ierr)

       allocate(sL_crsA_hi(s_nnzL))
       allocate(sL_crsA_lo(s_nnzL))
       allocate(sL_crsRow_ptr(addL+1))
       allocate(sL_crsCol(s_nnzL))
       ! -----
       i=1
       do k=n-addL+1,n
          !          write(*,*)" @ rank @ ",myrank,k
          kk=k-(n-addL+1) + 1
          sL_crsRow_ptr(kk)=i
          do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
             Col=crsMat%crsCol(j)
             if( (istartL<=Col).and.(Col<=iendL) )then
                sL_crsA_hi(i)=crsMat%crsA_hi(j)
                sL_crsA_lo(i)=crsMat%crsA_lo(j)
                sL_crsCol(i)=Col
                i=i+1
             endif
          enddo
       enddo
       sL_crsRow_ptr(addL+1)=i

    else
       allocate(sL_crsA_hi(1))
       allocate(sL_crsA_lo(1))
       allocate(sL_crsRow_ptr(2))
       allocate(sL_crsCol(1))
    endif
    if(myrank.ne.0)then
       call MPI_Irecv(r_nnzU,1,MPI_INTEGER4,Urank,Urank ,icomm,r_Ureq,ierr)
       call MPI_Isend(s_nnzU,1,MPI_INTEGER4,Urank,myrank,icomm,s_Ureq,ierr)

       allocate(sU_crsA_hi(s_nnzU))
       allocate(sU_crsA_lo(s_nnzU))
       allocate(sU_crsRow_ptr(addL+1))
       allocate(sU_crsCol(s_nnzU))
       i=1
       do k=1,addL
          sU_crsRow_ptr(k)=i
          do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
             Col=crsMat%crsCol(j)
             if( (istartU<=Col).and.(Col<=iendU) )then
                sU_crsA_hi(i)=crsMat%crsA_hi(j)
                sU_crsA_lo(i)=crsMat%crsA_lo(j)
                sU_crsCol(i)=Col
                i=i+1
             endif
          enddo
       enddo
       sU_crsRow_ptr(addL+1)=i
    else
       allocate(sU_crsA_hi(1))
       allocate(sU_crsA_lo(1))
       allocate(sU_crsRow_ptr(2))
       allocate(sU_crsCol(1))
    endif

    ! -----
    ! メモリの確保
    if(myrank.ne.(npes-1))then
       call MPI_Wait(r_Lreq,r_Lsta,ierr)
       call MPI_Wait(s_Lreq,s_Lsta,ierr)
       allocate(precon_crsMatASL%crsA_hi(r_nnzL))
       allocate(precon_crsMatASL%crsA_lo(r_nnzL))
       allocate(precon_crsMatASL%crsRow_ptr(addL+1))
       allocate(precon_crsMatASL%crsCol(r_nnzL))
       do i=1,r_nnzL
          precon_crsMatASL%crsA_hi(i)=0.0d0
          precon_crsMatASL%crsA_lo(i)=0.0d0
          precon_crsMatASL%crsCol(i)=0
       enddo
       do i=1,addL+1
          precon_crsMatASL%crsRow_ptr(i)=0
       enddo

       prm_ASL%n      = addL
       prm_ASL%gn     = prm_ASL%n
       prm_ASL%istart = istart - addL
       prm_ASL%iend   = istart - 1
       prm_ASL%nnz    = r_nnzL
    else
       allocate(precon_crsMatASL%crsA_hi(1))
       allocate(precon_crsMatASL%crsA_lo(1))
       allocate(precon_crsMatASL%crsRow_ptr(2))
       allocate(precon_crsMatASL%crsCol(1))
       precon_crsMatASL%crsA_hi(1)      =0.0d0
       precon_crsMatASL%crsA_lo(1)      =0.0d0
       precon_crsMatASL%crsRow_ptr(1)=1
       precon_crsMatASL%crsRow_ptr(2)=1
       precon_crsMatASL%crsCol(1)    =-1
       prm_ASL%n      =  0
       prm_ASL%istart = -1
       prm_ASL%iend   = -2
       prm_ASL%nnz    =  0
    endif

    if(myrank.ne.0)then
       call MPI_Wait(r_Ureq,r_Usta,ierr)
       call MPI_Wait(s_Ureq,s_Usta,ierr)
       allocate(precon_crsMatASU%crsA_hi(r_nnzU))
       allocate(precon_crsMatASU%crsA_lo(r_nnzU))
       allocate(precon_crsMatASU%crsRow_ptr(addL+1))
       allocate(precon_crsMatASU%crsCol(r_nnzU))

       do i=1,r_nnzU
          precon_crsMatASU%crsA_hi(i)=0.0d0
          precon_crsMatASU%crsA_lo(i)=0.0d0
          precon_crsMatASU%crsCol(i)=0
       enddo
       do i=1,addL+1
          precon_crsMatASU%crsRow_ptr(i)=0
       enddo
       prm_ASU%n      = addL
       prm_ASU%gn     = prm_ASU%n 
       prm_ASU%istart = istart + addL
       prm_ASU%iend   = istart + 1
       prm_ASU%nnz    = r_nnzU
    else
       allocate(precon_crsMatASU%crsA_hi(1))
       allocate(precon_crsMatASU%crsA_lo(1))
       allocate(precon_crsMatASU%crsRow_ptr(2))
       allocate(precon_crsMatASU%crsCol(1))
       precon_crsMatASU%crsA_hi(1)      =0.0d0
       precon_crsMatASU%crsA_lo(1)      =0.0d0
       precon_crsMatASU%crsRow_ptr(1)=1
       precon_crsMatASU%crsRow_ptr(2)=1
       precon_crsMatASU%crsCol(1)    =-1
       prm_ASU%n      =  0
       prm_ASU%istart = -1
       prm_ASU%iend   = -2
       prm_ASU%nnz    =  0
    endif

    ! 隣接プロセスの行列要素を受信
    addL1=addL+1
    if(myrank.ne.0)then
       call MPI_Irecv(precon_crsMatASU%crsA_hi,r_nnzU,MPI_REAL8,Urank,Urank ,icomm,r_Ureq,ierr)
       call MPI_Isend(sU_crsA_hi     ,s_nnzU,MPI_REAL8,Urank,myrank,icomm,s_Ureq,ierr)
       call MPI_Wait(r_Ureq,r_Usta,ierr)
       call MPI_Wait(s_Ureq,s_Usta,ierr)

       call MPI_Irecv(precon_crsMatASU%crsA_lo,r_nnzU,MPI_REAL8,Urank,Urank ,icomm,r_Ureq,ierr)
       call MPI_Isend(sU_crsA_lo     ,s_nnzU,MPI_REAL8,Urank,myrank,icomm,s_Ureq,ierr)
       call MPI_Wait(r_Ureq,r_Usta,ierr)
       call MPI_Wait(s_Ureq,s_Usta,ierr)

!        call MPI_Irecv(precon_crsMatASU%crsCol,r_nnzU,MPI_INTEGER4,Urank,Urank ,icomm,r_Ureq,ierr)
!        call MPI_Isend(sU_crsCol     ,s_nnzU,MPI_INTEGER4,Urank,myrank,icomm,s_Ureq,ierr)
       call MPI_Irecv(precon_crsMatASU%crsCol,r_nnzU,intMPIType,Urank,Urank ,icomm,r_Ureq,ierr)
       call MPI_Isend(sU_crsCol              ,s_nnzU,intMPIType,Urank,myrank,icomm,s_Ureq,ierr)
       call MPI_Wait(r_Ureq,r_Usta,ierr)
       call MPI_Wait(s_Ureq,s_Usta,ierr)

!        call MPI_Irecv(precon_crsMatASU%crsRow_ptr,addL1,MPI_INTEGER4,Urank,Urank ,icomm,r_Ureq,ierr)
!        call MPI_Isend(sU_crsRow_ptr     ,addL1,MPI_INTEGER4,Urank,myrank,icomm,s_Ureq,ierr)
       call MPI_Irecv(precon_crsMatASU%crsRow_ptr,addL1,intMPIType,Urank,Urank ,icomm,r_Ureq,ierr)
       call MPI_Isend(sU_crsRow_ptr              ,addL1,intMPIType,Urank,myrank,icomm,s_Ureq,ierr)
       call MPI_Wait(r_Ureq,r_Usta,ierr)
       call MPI_Wait(s_Ureq,s_Usta,ierr)
    endif

    if(myrank.ne.(npes-1))then
       call MPI_Irecv(precon_crsMatASL%crsA_hi,r_nnzL,MPI_REAL8,Lrank,Lrank ,icomm,r_Lreq,ierr)
       call MPI_Isend(sL_crsA_hi     ,s_nnzL,MPI_REAL8,Lrank,myrank,icomm,s_Lreq,ierr)
       call MPI_Wait(r_Lreq,r_Lsta,ierr)
       call MPI_Wait(s_Lreq,s_Lsta,ierr)
       call MPI_Irecv(precon_crsMatASL%crsA_lo,r_nnzL,MPI_REAL8,Lrank,Lrank ,icomm,r_Lreq,ierr)
       call MPI_Isend(sL_crsA_lo     ,s_nnzL,MPI_REAL8,Lrank,myrank,icomm,s_Lreq,ierr)
       call MPI_Wait(r_Lreq,r_Lsta,ierr)
       call MPI_Wait(s_Lreq,s_Lsta,ierr)

!        call MPI_Irecv(precon_crsMatASL%crsCol,r_nnzL,MPI_INTEGER4,Lrank,Lrank ,icomm,r_Lreq,ierr)
!        call MPI_Isend(sL_crsCol     ,s_nnzL,MPI_INTEGER4,Lrank,myrank,icomm,s_Lreq,ierr)
       call MPI_Irecv(precon_crsMatASL%crsCol,r_nnzL,intMPIType,Lrank,Lrank ,icomm,r_Lreq,ierr)
       call MPI_Isend(sL_crsCol              ,s_nnzL,intMPIType,Lrank,myrank,icomm,s_Lreq,ierr)
       call MPI_Wait(r_Lreq,r_Lsta,ierr)
       call MPI_Wait(s_Lreq,s_Lsta,ierr)

!        call MPI_Irecv(precon_crsMatASL%crsRow_ptr,addL1,MPI_INTEGER4,Lrank,Lrank ,icomm,r_Lreq,ierr)
!        call MPI_Isend(sL_crsRow_ptr     ,addL1,MPI_INTEGER4,Lrank,myrank,icomm,s_Lreq,ierr)
       call MPI_Irecv(precon_crsMatASL%crsRow_ptr,addL1,intMPIType,Lrank,Lrank ,icomm,r_Lreq,ierr)
       call MPI_Isend(sL_crsRow_ptr              ,addL1,intMPIType,Lrank,myrank,icomm,s_Lreq,ierr)
       call MPI_Wait(r_Lreq,r_Lsta,ierr)
       call MPI_Wait(s_Lreq,s_Lsta,ierr)
    endif

    ! -----
    ! 加法シュワルツ用の前処理行列作成

    ! 対角ブロックの非零要素数
    nnzAS=0
    do k=1,n
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          if( crsMat%crsCol(j)<istart - addL ) then
             cycle
          endif
          if( iend + addL<crsMat%crsCol(j) )then
             cycle
          endif
          nnzAS=nnzAS+1
       enddo
    enddo

    if(myrank.ne.0)then
       nnzAS=nnzAS+r_nnzU
    endif
    if(myrank.ne.(npes-1))then
       nnzAS=nnzAS+r_nnzL
    endif

    method%prm_AS%istart = istart - addL
    method%prm_AS%iend   = iend   + addL    
    method%prm_AS%nnz    = nnzAS
    if(method%prm_AS%istart<=0)then
       method%prm_AS%istart=1
    endif
    if(method%prm_AS%iend>gn)then
       method%prm_AS%iend=gn
    endif
    method%prm_AS%n  = method%prm_AS%iend - method%prm_AS%istart + 1
    method%prm_AS%gn = method%prm_AS%n

    allocate(method%precon_crsMatAS%crsA_hi(method%prm_AS%nnz))
    allocate(method%precon_crsMatAS%crsA_lo(method%prm_AS%nnz))
    allocate(method%precon_crsMatAS%crsRow_ptr(method%prm_AS%n+1))
    allocate(method%precon_crsMatAS%crsCol(method%prm_AS%nnz))
    call zero_initialize(method%prm_AS%nnz      , method%precon_crsMatAS%crsA_hi)
    call zero_initialize(method%prm_AS%nnz      , method%precon_crsMatAS%crsA_lo)

    i0=method%prm_AS%n+1
    call zero_initialize_int(i0                , method%precon_crsMatAS%crsRow_ptr)
    call zero_initialize_int(method%prm_AS%nnz , method%precon_crsMatAS%crsCol)

    istartAS = method%prm_AS%istart
    iendAS   = method%prm_AS%iend  
    i=1
    kk=1
    if(myrank.ne.0)then
       do k=1,addL
          method%precon_crsMatAS%crsRow_ptr(k+kk-1)=i
          do j=precon_crsMatASU%crsRow_ptr(k),precon_crsMatASU%crsRow_ptr(k+1)-1
             Col=precon_crsMatASU%crsCol(j)
             if(Col<istartAS)then
                cycle
             endif
             if(iendAS<Col)then
                cycle
             endif
             method%precon_crsMatAS%crsA_hi(i)  = precon_crsMatASU%crsA_hi(j)
             method%precon_crsMatAS%crsA_lo(i)  = precon_crsMatASU%crsA_lo(j)
             method%precon_crsMatAS%crsCol(i)= precon_crsMatASU%crsCol(j)
             i=i+1
          enddo
       enddo
       kk=kk+addL
    endif

    do k=1,n
       method%precon_crsMatAS%crsRow_ptr(k+kk-1)=i
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          Col=crsMat%crsCol(j)
          if(Col<istartAS)then
             cycle
          endif
          if(iendAS<Col)then
             cycle
          endif
          method%precon_crsMatAS%crsA_hi(i)  = crsMat%crsA_hi(j)
          method%precon_crsMatAS%crsA_lo(i)  = crsMat%crsA_lo(j)
          method%precon_crsMatAS%crsCol(i)= crsMat%crsCol(j)
          i=i+1
       enddo
    enddo
    kk=kk+n

    if(myrank.ne.(npes-1))then
       do k=1,addL
          method%precon_crsMatAS%crsRow_ptr(k+kk-1)=i
          do j=precon_crsMatASL%crsRow_ptr(k),precon_crsMatASL%crsRow_ptr(k+1)-1
             Col=precon_crsMatASL%crsCol(j)
             if(Col<istartAS)then
                cycle
             endif
             if(iendAS<Col)then
                cycle
             endif
             method%precon_crsMatAS%crsA_hi(i)  = precon_crsMatASL%crsA_hi(j)
             method%precon_crsMatAS%crsA_lo(i)  = precon_crsMatASL%crsA_lo(j)
             method%precon_crsMatAS%crsCol(i)= precon_crsMatASL%crsCol(j)
             i=i+1
          enddo
       enddo
       kk=kk+addL+1
    endif
    method%precon_crsMatAS%crsRow_ptr(method%prm_AS%n+1)=i

    method%precon_crsMatAS%nColL=0
    method%precon_crsMatAS%nColU=0

    ! -----
    ! D-ILU分解 
    allocate(method%ivec_AS_hi(n+addL*2))
    allocate(method%ivec_AS_lo(n+addL*2))
    allocate(method%ovec_AS_hi(n+addL*2))
    allocate(method%ovec_AS_lo(n+addL*2))
!    write(*,*)"ivec_AS_hi : size = ",n+addL*2
    call zero_initialize(n+addL*2,method%ivec_AS_hi)
    call zero_initialize(n+addL*2,method%ivec_AS_lo)
    call zero_initialize(n+addL*2,method%ovec_AS_hi)
    call zero_initialize(n+addL*2,method%ovec_AS_lo)

    ! -----
    !    call get_maxthreads(max_numth)
    max_numth=method%max_numth
    allocate(method%preconAS(max_numth)) 
    allocate(istart_th_AS(max_numth)) 
    allocate(iend_th_AS(max_numth)) 
    do i=1,max_numth
       method%preconAS(i)%istart_th = -1
       method%preconAS(i)%iend_th   = -1
       istart_th_AS(i) = -1
       iend_th_AS(i)   = -1
    enddo
    call get_AS_region( &
         method%prm,addL, &
         Lrank,Urank, &
         max_numth, &
         istart_th_AS, &
         iend_th_AS, &
         method%precon_thblock)
    do i=1,max_numth
       method%preconAS(i)%istart_th = istart_th_AS(i)
       method%preconAS(i)%iend_th   = iend_th_AS(i)

       method%preconAS(i)%n_startth = istart_th_AS(i) - istartAS +1
       method%preconAS(i)%n_endth   = &
            method%preconAS(i)%n_startth + &
            ( iend_th_AS(i)-istart_th_AS(i) )

       allocate(method%preconAS(i)%Dinv_hi( &
            method%preconAS(i)%n_startth:method%preconAS(i)%n_endth))
       allocate(method%preconAS(i)%Dinv_lo( &
            method%preconAS(i)%n_startth:method%preconAS(i)%n_endth))

       allocate(method%preconAS(i)%svec_hi( &
            method%preconAS(i)%n_startth:method%preconAS(i)%n_endth)) 
       allocate(method%preconAS(i)%svec_lo( &
            method%preconAS(i)%n_startth:method%preconAS(i)%n_endth)) 
       allocate(method%preconAS(i)%rvec_hi( &
            method%preconAS(i)%n_startth:method%preconAS(i)%n_endth)) 
       allocate(method%preconAS(i)%rvec_lo( &
            method%preconAS(i)%n_startth:method%preconAS(i)%n_endth)) 

       n_th = method%preconAS(i)%n_endth-method%preconAS(i)%n_startth+1
       call zero_initialize(n_th,method%preconAS(i)%Dinv_hi)
       call zero_initialize(n_th,method%preconAS(i)%Dinv_lo)
       call zero_initialize(n_th,method%preconAS(i)%svec_hi)
       call zero_initialize(n_th,method%preconAS(i)%svec_lo)
       call zero_initialize(n_th,method%preconAS(i)%rvec_hi)
       call zero_initialize(n_th,method%preconAS(i)%rvec_lo)
    enddo

    ! ----
    ! 領域拡大前の各スレッドの計算領域を取得
    do i=1,max_numth
       method%preconAS(i)%istart0_th = -1
       method%preconAS(i)%iend0_th   = -1
       istart_th_AS(i) = -1
       iend_th_AS(i)   = -1
    enddo
    call get_AS_region( &
       method%prm,0, &
       Lrank,Urank, &
       max_numth, &
       istart_th_AS, &
       iend_th_AS, &
       method%precon_thblock)
    do i=1,max_numth
       method%preconAS(i)%istart0_th = istart_th_AS(i)
       method%preconAS(i)%iend0_th   = iend_th_AS(i)
       method%preconAS(i)%n0_startth = istart_th_AS(i) - istartAS +1
       method%preconAS(i)%n0_endth   = &
            method%preconAS(i)%n0_startth + &
            ( iend_th_AS(i)-istart_th_AS(i) )
    enddo
    ! ---

    deallocate(istart_th_AS)
    deallocate(iend_th_AS)

    call ddDILU_DInv_AS_crs( &
         method%prm_AS, &
         method%precon_crsMatAS, &
         method%preconAS,&
         max_numth, &
         method%ILU_method, &
         method%independ_nvec, &
         method%precon_thblock &
         )

    allocate(method%precon_crsMatAS%RowL(0:0))
    allocate(method%precon_crsMatAS%RowD(method%prm_AS%n))
    allocate(method%precon_crsMatAS%RowU(0:0))
    method%precon_crsMatAS%RowL(0)=-1
    !$OMP PARALLEL DO PRIVATE(i)
    do i=1,method%prm_AS%n
       method%precon_crsMatAS%RowD(i) = &
            i + method%prm_AS%istart - 1
    enddo
    !$OMP END PARALLEL DO
    method%precon_crsMatAS%RowU(0)=-1

    ! ------
    deallocate(sL_crsA_hi)
    deallocate(sL_crsA_lo)
    deallocate(sL_crsRow_ptr)
    deallocate(sL_crsCol)
    deallocate(sU_crsA_hi)
    deallocate(sU_crsA_lo)
    deallocate(sU_crsRow_ptr)
    deallocate(sU_crsCol)

    deallocate(precon_crsMatASL%crsA_hi)
    deallocate(precon_crsMatASL%crsA_lo)
    deallocate(precon_crsMatASL%crsRow_ptr)
    deallocate(precon_crsMatASL%crsCol)

    deallocate(precon_crsMatASU%crsA_hi)
    deallocate(precon_crsMatASU%crsA_lo)
    deallocate(precon_crsMatASU%crsRow_ptr)
    deallocate(precon_crsMatASU%crsCol)

  end subroutine ddget_DILU_DInv_AS_crs

  subroutine ddDILU_DInv_AS_crs( &
       prm,crsMat,preconAS,mnumth,ILU_method,independ_nvec,precon_thblock)
    use omp_lib
    use precon
    use matrix
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(CRS_Mat) :: crsMat 
    !    real*8 Dinv(prm%n)
    type(precon_AS) preconAS(mnumth)
    integer mnumth
    integer ILU_method
    integer independ_nvec
    integer precon_thblock
    integer(kind=intType) n,istart,iend
    integer numth,myth,max_numth
    integer(kind=intType) n_th
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th

    n = prm%n
    istart=prm%istart
    iend=prm%iend

    ! -----
    ! ブロック数よりもスレッド数が大きい場合はスレッド数を制限する
    call get_maxthreads(max_numth)
    call set_precon_maxthreads(n,precon_thblock)

    !$OMP PARALLEL &
    !$OMP PRIVATE( &
    !$OMP myth &
    !$OMP )
    myth  = omp_get_thread_num()
    call ddDILU_DInv_crs( &
         crsMat, &
         preconAS(myth+1)%Dinv_hi(preconAS(myth+1)%n_startth), &
         preconAS(myth+1)%Dinv_lo(preconAS(myth+1)%n_startth), &
         istart,iend, &
         preconAS(myth+1)%n_startth, &
         preconAS(myth+1)%n_endth, &
         preconAS(myth+1)%istart_th, &
         preconAS(myth+1)%iend_th, &
         ILU_method, &
         independ_nvec &
         )
    !$OMP END PARALLEL
!    if(n<=max_numth)then
    call omp_set_num_threads(max_numth)
!    endif

  end subroutine ddDILU_DInv_AS_crs

  subroutine ddget_DILU_DInv_crs(&
       prm,crsMat,Dinv_hi,Dinv_lo,ILU_method,independ_nvec,precon_thblock)
    use omp_lib
    use matrix
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(CRS_Mat) :: crsMat 
    real*8 Dinv_hi(prm%n)
    real*8 Dinv_lo(prm%n)
    integer independ_nvec
    integer precon_thblock

    integer(kind=intType) n,istart,iend
    integer numth,myth,max_numth
    integer(kind=intType) n_th
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th
    integer ILU_method
    integer(kind=intType) ntmp

    n = prm%n
    istart=prm%istart
    iend=prm%iend

    ! -----
    ! ブロック数よりもスレッド数が大きい場合はスレッド数を制限する
    call get_maxthreads(max_numth)
    call set_precon_maxthreads(n,precon_thblock)

    ! -----
    ! istart_th=istart
    ! iend_th=iend
    ! n_startth=1
    ! n_endth=n

    !$OMP PARALLEL &
    !$OMP PRIVATE( &
    !$OMP numth,myth,n_startth,n_endth,istart_th,iend_th,ntmp &
    !$OMP )
    numth = omp_get_num_threads()
    myth  = omp_get_thread_num()
    ntmp=1
    call divide_region( &
         myth,numth,ntmp,n,n_startth,n_endth)
    istart_th = istart -1 + n_startth 
    iend_th   = istart_th + (n_endth - n_startth) 
!    write(*,*)myth,istart_th,iend_th  

    call ddDILU_DInv_crs( &
         crsMat, &
         Dinv_hi(n_startth), &
         Dinv_lo(n_startth), &
         istart,iend, &
         n_startth,n_endth, &
         istart_th,iend_th, &
         ILU_method, &
         independ_nvec &
         )
    !$OMP END PARALLEL

    !    if(n<=max_numth)then
    call omp_set_num_threads(max_numth)
    !    endif

  end subroutine ddget_DILU_DInv_crs

  subroutine ddDILU_DInv_crs( &
       crsMat, &
       Dinv_hi,Dinv_lo, &
       istart,iend, &
       n_startth,n_endth, &
       istart_th,iend_th, &
       ILU_method, &
       independ_nvec &
       )
    use omp_lib
    use matrix
    use PARCEL_integer
    implicit none
    type(CRS_Mat) :: crsMat 
    integer(kind=intType) istart,iend
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th
    real*8 Dinv_hi(n_startth:n_endth)
    real*8 Dinv_lo(n_startth:n_endth)
    integer ILU_method
    integer independ_nvec

    real*8 D_hi(n_startth:n_endth)
    real*8 D_lo(n_startth:n_endth)
    real*8 c_hi,c_lo

    integer(kind=intType) i,j,k
    integer(kind=intType) Row,Col
    real*8 ddmul_var1
    real*8 dddiv_var1
    real*8 ddsub_var1

    integer(kind=intType) kr
    integer(kind=intType) jr,rCol,rRow
    integer(kind=intType) jstart
    integer independ_start,independ_end

    ! -----
    ! D-ILU分解
    ! A = AL + AD + AU ~= K
    ! K = (D+AL)D^-1(D+AU)
    ! 条件1 対角成分が一致     :: D(i,i) = A(i,i) - Sum[j<i]          A(i,j) D(j,j)^-1 A(j,i)
    ! 条件2 行内の要素和が一致 :: D(i,i) = A(i,i) - Sum[j<i] Sum[k>j] A(i,j) D(j,j)^-1 A(j,k)
    ! -----
    do k=n_startth,n_endth
       ! 独立に回すループサイズ
       independ_start = n_startth + ((k-n_startth)/independ_nvec) * independ_nvec
       independ_end   = min(independ_start+independ_nvec-1,n_endth)

       !       Row=k-1+istart_th
       Row=k-n_startth+istart_th
       D_hi(k)=0.0d0
       D_lo(k)=0.0d0
       jstart=-1
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          if( crsMat%crsCol(j)<istart_th ) then
             cycle
          endif
          if( iend_th<crsMat%crsCol(j) )then
             cycle
          endif
          ! -----
          ! 対角ブロック
          !          Col=crsMat%crsCol(j)-istart_th+1
          Col=crsMat%crsCol(j)-istart+1
          if(jstart.eq.-1)then
             ! 対角ブロックの開始位置
             jstart=j
          endif
          ! -----
          !          if(k.eq.Col)then
          !          write(*,*)"k =",k,Col
          if(k<=Col)then

             ! -----
             !@ @! 行列の対称位置に対して値があれば総和を計算する
             ! 下側部分を探査する

             ! 対角要素による初期化
             if(k.eq.Col)then
                D_hi(k)=crsMat%crsA_hi(j) 
                D_lo(k)=crsMat%crsA_lo(j) 
             else
                D_hi(k)=0.0d0
                D_lo(k)=0.0d0
             endif
             ! -----
             ! jは対角部分より大きいので -1
             ! 下側列の探査
             do jr=jstart,j-1
                !                Col=crsMat%crsCol(jr)-istart_th+1
                Col=crsMat%crsCol(jr)-istart+1
                kr=Col
                
                if( (independ_start<=Col).and.(Col<=independ_end) )then
                   if(Col.ne.k)then                
                      cycle
                   endif
                endif

                ! 上側行(kr行)の探査
                do i=crsMat%crsRow_ptr(kr),crsMat%crsRow_ptr(kr+1)-1
                   ! -----
                   ! 上三角部分の総和を計算する
                   if( crsMat%crsCol(i)<istart_th ) then
                      cycle
                   endif
                   if( iend_th<crsMat%crsCol(i) )then
                      exit
                   endif
                   ! -----
                   ! 対角ブロック
                   rCol=crsMat%crsCol(i)-istart+1

                   if(ILU_method.eq.KSP_PRECONDITION_ILUD)then
                      ! 対角成分が一致
                      if(k.eq.rCol)then
                         !                         D(k)=D(k)-crsMat%crsA(jr)*crsMat%crsA(i)/D(kr)
#ifdef Manual_inline
                         ddmul_calc1(crsMat%crsA_hi(jr),crsMat%crsA_lo(jr),crsMat%crsA_hi(i),crsMat%crsA_lo(i),c_hi,c_lo)
                         if( (D_hi(kr).eq.0.0d0).and.&
                             (D_lo(kr).eq.0.0d0))then
                            dddiv_calc1(c_hi,c_lo,1.0d-100,0.0d0,c_hi,c_lo)
                         else
                            dddiv_calc1(c_hi,c_lo,D_hi(kr),D_lo(kr),c_hi,c_lo)
                         endif
                         ddsub_calc1(D_hi(k),D_lo(k),c_hi,c_lo,D_hi(k),D_lo(k))
#else
                         call dd_mul( &
                              crsMat%crsA_hi(jr),crsMat%crsA_lo(jr), &
                              crsMat%crsA_hi(i),crsMat%crsA_lo(i), &
                              c_hi,c_lo)
                         if( (D_hi(kr).eq.0.0d0).and.&
                             (D_lo(kr).eq.0.0d0))then
                            call dd_div( &
                                 c_hi,c_lo,&
                                 1.0d-100,0.0d0, &
                                 c_hi,c_lo)
                         else
                            call dd_div( &
                                 c_hi,c_lo,&
                                 D_hi(kr),D_lo(kr), &
                                 c_hi,c_lo)
                         endif
                         call dd_sub(D_hi(k),D_lo(k),c_hi,c_lo,D_hi(k),D_lo(k))
#endif
                         exit
                      endif
                   else
                      ! 行内の要素和が一致
                      if(kr<rCol)then
                         !                         D(k)=D(k)-crsMat%crsA(jr)*crsMat%crsA(i)/D(kr)
#ifdef Manual_inline
                         ddmul_calc1(crsMat%crsA_hi(jr),crsMat%crsA_lo(jr),crsMat%crsA_hi(i),crsMat%crsA_lo(i),c_hi,c_lo)
                         if( (D_hi(kr).eq.0.0d0).and.&
                             (D_lo(kr).eq.0.0d0))then
                            dddiv_calc1(c_hi,c_lo,1.0d-100,0.0d0,c_hi,c_lo)
                         else
                            dddiv_calc1(c_hi,c_lo,D_hi(kr),D_lo(kr),c_hi,c_lo)
                         endif
                         ddsub_calc1(D_hi(k),D_lo(k),c_hi,c_lo,D_hi(k),D_lo(k))
#else
                         call dd_mul( &
                              crsMat%crsA_hi(jr),crsMat%crsA_lo(jr), &
                              crsMat%crsA_hi(i),crsMat%crsA_lo(i), &
                              c_hi,c_lo)

                         if( (D_hi(kr).eq.0.0d0).and.&
                             (D_lo(kr).eq.0.0d0))then
                            call dd_div(&
                                 c_hi,c_lo,&
                                 1.0d-100,0.0d0,&
                                 c_hi,c_lo)
                         else
                            call dd_div(&
                                 c_hi,c_lo,&
                                 D_hi(kr),D_lo(kr),&
                                 c_hi,c_lo)
                         endif
                         call dd_sub(D_hi(k),D_lo(k),c_hi,c_lo,D_hi(k),D_lo(k))
#endif
                      endif
                   endif

                enddo
             enddo
             exit
          endif

       enddo
    enddo

    ! endif
    do k=n_startth,n_endth
       !       Dinv(k)=1.0d0/D(k)
#ifdef Manual_inline
       if( (D_hi(k).eq.0.0d0).and.&
            (D_lo(k).eq.0.0d0))then
          dddiv_calc1(1.0d0,0.0d0,1.0d-100,0.0d0,Dinv_hi(k),Dinv_lo(k))
       else
          dddiv_calc1(1.0d0,0.0d0,D_hi(k),D_lo(k),Dinv_hi(k),Dinv_lo(k))
       endif
#else
       if( (D_hi(k).eq.0.0d0).and.&
            (D_lo(k).eq.0.0d0))then
          call dd_div(1.0d0,0.0d0,1.0d-100,0.0d0,Dinv_hi(k),Dinv_lo(k))
       else
          call dd_div(1.0d0,0.0d0,D_hi(k),D_lo(k),Dinv_hi(k),Dinv_lo(k))
       endif
#endif
!       write(*,*)k,"Dinv - ",D_hi(k) + D_lo(k)
    enddo

  end subroutine ddDILU_DInv_crs

  subroutine ddcrs2dia(crsMat,diaMat,n,istart,nnd,nnz)
    use sort_kernel
    use matrix
    use PARCEL_integer
    implicit none
    type(CRS_Mat) :: crsMat 
    type(DIA_Mat) :: diaMat 

    integer(kind=intType) nnd
    integer(kind=intType) nnz
    integer(kind=intType) istart
    integer(kind=intType) n

    integer(kind=intType) i,j,k
    integer(kind=intType) offset
    integer(kind=intType) offset_sort(nnz)
    integer(kind=intType) Row,Col

    ! ----
    i=1
    do k=1,n
       Row=k-1+istart
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          Col=crsMat%crsCol(j)
          offset_sort(i)=Col-Row
          i=i+1
       enddo
    enddo
    !    call sort_c(offset_sort,nnz)
    call MergeSort(offset_sort,nnz)
    nnd=1
    do i=1,nnz
       offset=offset_sort(i)
       if(offset_sort(nnd).ne.offset)then
          nnd=nnd+1
          offset_sort(nnd)=offset
       endif
    enddo

    ! -----
    do i=1,nnd
       diaMat%offset(i)=offset_sort(i)
    enddo
    do k=1,n
       Row=k-1+istart
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          Col=crsMat%crsCol(j)
          offset=Col-Row
          do i=1,nnd
             if(offset.eq.diaMat%offset(i))then
                diaMat%diaA_hi(k+(i-1)*n)=crsMat%crsA_hi(j)
                diaMat%diaA_lo(k+(i-1)*n)=crsMat%crsA_lo(j)
                exit
             endif
          enddo
       enddo
    enddo

  end subroutine ddcrs2dia

  subroutine ddfree_crsMat(crsMat)
    use matrix
    implicit none
    type(CRS_Mat)  :: crsMat
    deallocate(crsMat%crsA_hi)
    deallocate(crsMat%crsA_lo)
    deallocate(crsMat%crsRow_ptr)
    if(allocated(crsMat%crsRow_ptr_Diag_s))then
       deallocate(crsMat%crsRow_ptr_Diag_s)
    endif
    if(allocated(crsMat%crsRow_ptr_Diag_e))then
       deallocate(crsMat%crsRow_ptr_Diag_e)
    endif
    deallocate(crsMat%crsCol)
    deallocate(crsMat%RowL)
    deallocate(crsMat%RowD)
    deallocate(crsMat%RowU)
  end subroutine ddfree_crsMat

  subroutine ddfree_diaMat(diaMat)
    use matrix
    implicit none
    type(DIA_Mat)  :: diaMat
    deallocate(diaMat%diaA_hi)
    deallocate(diaMat%diaA_lo)
    deallocate(diaMat%offset)
    deallocate(diaMat%offsetVec)
    deallocate(diaMat%Vec_hi)
    deallocate(diaMat%Vec_lo)

    deallocate(diaMat%RowL)
    deallocate(diaMat%RowD)
    deallocate(diaMat%RowU) 

  end subroutine ddfree_diaMat

  subroutine ddfree_precon_AS(method)
    use krylov
    implicit none
    type(KSP) :: method
    integer max_numth
    integer i
    max_numth=method%max_numth
    do i=1,max_numth
       deallocate(method%preconAS(i)%Dinv_hi)
       deallocate(method%preconAS(i)%Dinv_lo)
       deallocate(method%preconAS(i)%rvec_hi)
       deallocate(method%preconAS(i)%rvec_lo)
       deallocate(method%preconAS(i)%svec_hi)
       deallocate(method%preconAS(i)%svec_lo)
    enddo
    deallocate(method%preconAS)

  end subroutine ddfree_precon_AS

! todo
! 倍精度のSpMVとループの回し方がちがうので修正する

!  subroutine MatVec_crs(prm,crsMat,xL,xD,xU,Ax,nBlock)
  subroutine ddMatVec_crs( &
       prm,crsMat,&
       xL_hi,xL_lo, &
       xD_hi,xD_lo, &
       xU_hi,xU_lo, &
       Ax_hi,Ax_lo, &
       nBlock, &
       PRECISION_A)
    use matrix
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(CRS_Mat) :: crsMat 
    ! real*8 xL(crsMat%nColL)
    ! real*8 xD(prm%n)
    ! real*8 xU(crsMat%nColU)

    ! real*8 tmp
    ! real*8 Ax(prm%n)
    ! real*8 Ax0(nBlock)

    real*8 xL_hi(crsMat%nColL)
    real*8 xL_lo(crsMat%nColL)
    real*8 xD_hi(prm%n)
    real*8 xD_lo(prm%n)
    real*8 xU_hi(crsMat%nColU)
    real*8 xU_lo(crsMat%nColU)

    real*8 Ax_hi(prm%n)
    real*8 Ax_lo(prm%n)
    real*8 tmp_hi,tmp_lo

    ! real*8 tmp_A_hi(tmpblocksize,tmpblocksize)
    ! real*8 tmp_A_lo(tmpblocksize,tmpblocksize)    
    ! real*8 tmp_vec_hi(tmpblocksize,tmpblocksize)
    ! real*8 tmp_vec_lo(tmpblocksize,tmpblocksize)
    
    ! real*8 tmp_block_hi(tmpblocksize,10)
    ! real*8 tmp_block_lo(tmpblocksize,10)
    integer,parameter :: BlockRowMax=50

    real*8 tmp_A_hi(tmpblocksize  ,BlockRowMax)
    real*8 tmp_A_lo(tmpblocksize  ,BlockRowMax)    
    real*8 tmp_vec_hi(tmpblocksize,BlockRowMax)
    real*8 tmp_vec_lo(tmpblocksize,BlockRowMax)
   
    real*8 Ax0_hi(nBlock)
    real*8 Ax0_lo(nBlock)

    integer(kind=intType) iptr_s,iptr_e

    integer(kind=intType) i,j,jj,k,kk,kkmin,kkstep,kk0
    integer(kind=intType) js,je,LocalRowMax,LocalRow

    integer(kind=intType) iLs,iUs,iLe,iUe
    integer(kind=intType) n,istart,iend
    integer(kind=intType) Row,Col
    integer(kind=intType) nColL,nColU
    integer nBlock
    integer PRECISION_A
    real*8 c_hi,c_lo
    real*8 ddmul_var1
    real*8 ddadd_var1

    n = prm%n
    istart=prm%istart
    iend=prm%iend
    nColL=crsMat%nColL
    nColU=crsMat%nColU

    select case (PRECISION_A)
    case(KSP_DOUBLE_PRECISION)
#if 1
       !$OMP PARALLEL DO PRIVATE(i,j,k,kk,iLs,iUs,iLe,iUe,Row,Col,tmp_hi,&
       !$OMP tmp_lo,c_hi,c_lo,js,je, &
       !$OMP kkmin,LocalRowMax,LocalRow, &
!       !$OMP tmp_block_hi,tmp_block_lo,&
       !$OMP tmp_A_hi,tmp_A_lo,tmp_vec_hi,tmp_vec_lo, &
       !$OMP ddmul_var1, &
       !$OMP ddadd_var1) &
       !$OMP SCHEDULE(dynamic)
       do kk=1,n,tmpblocksize
          kkmin = min(kk+tmpblocksize-1,n)
          LocalRowMax = -1
          do k=kk,kkmin
             LocalRow = crsMat%crsRow_ptr(k+1) - crsMat%crsRow_ptr(k) 
             if(LocalRowMax<LocalRow)then
                LocalRowMax = LocalRow
             endif
          enddo
!          write(*,*)"LocalRowMax = ",LocalRowMax                       
          !          if(LocalRowMax>tmpblocksize)then
          if(LocalRowMax>BlockRowMax)then
             do k=kk,kkmin
                Ax_hi(k)=0.0d0
                Ax_lo(k)=0.0d0
                tmp_hi=0.0d0
                tmp_lo=0.0d0
                Row=k-1+istart
                do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
                   if( crsMat%crsCol(j)<istart ) then
                      Col=crsMat%iLs(j)
                      !Ax(k)=Ax(k)+crsMat%crsA(j)*xL(Col)
                      ! tmp=tmp+crsMat%crsA(j)*xL(Col)
                      ddmul_calc1(crsMat%crsA_hi(j),0.0d0,xL_hi(Col),xL_lo(Col),c_hi,c_lo)
                      ddadd_calc1(c_hi,c_lo,tmp_hi,tmp_lo,tmp_hi,tmp_lo)
                      cycle
                   endif
                   js=j
                   exit
                enddo
                do j=crsMat%crsRow_ptr(k+1)-1,crsMat%crsRow_ptr(k),-1
                   if( iend<crsMat%crsCol(j) )then
                      Col=crsMat%iLs(j)
                      ! tmp=tmp+crsMat%crsA(j)*xU(Col)
                      ddmul_calc1(crsMat%crsA_hi(j),0.0d0,xU_hi(Col),xU_lo(Col),c_hi,c_lo)
                      ddadd_calc1(c_hi,c_lo,tmp_hi,tmp_lo,tmp_hi,tmp_lo)
                      cycle
                   endif
                   je=j
                   exit
                enddo
                
                !       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
                do j=js,je
                   !          Col=crsMat%crsCol(j)-istart+1
                   Col=crsMat%iLs(j)
                   !Ax(k)=Ax(k)+crsMat%crsA(j)*xD(Col)
                   ! tmp=tmp+crsMat%crsA(j)*xD(Col)
                   ddmul_calc1(crsMat%crsA_hi(j),0.0d0,xD_hi(Col),xD_lo(Col),c_hi,c_lo)
                   ddadd_calc1(c_hi,c_lo,tmp_hi,tmp_lo,tmp_hi,tmp_lo)
                enddo
                Ax_hi(k)=tmp_hi
                Ax_lo(k)=tmp_lo
             enddo             
          else
             do k=kk,kkmin
                Ax_hi(k)=0.0d0
                Ax_lo(k)=0.0d0                
             enddo
             tmp_A_hi = 0.0d0
             tmp_A_lo = 0.0d0    
             tmp_vec_hi = 0.0d0
             tmp_vec_lo = 0.0d0
             
             do k=kk,kkmin
                tmp_hi=0.0d0
                tmp_lo=0.0d0
                
                Row=k-1+istart
                do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
                   if( crsMat%crsCol(j)<istart ) then
                      Col=crsMat%iLs(j)
                      ! ddmul_calc1(crsMat%crsA_hi(j),0.0d0,xL_hi(Col),xL_lo(Col),c_hi,c_lo)
                      ! ddadd_calc1(c_hi,c_lo,tmp_hi,tmp_lo,tmp_hi,tmp_lo)
                      tmp_A_hi(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = crsMat%crsA_hi(j)         
                      tmp_A_lo(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = 0.0d0
                      tmp_vec_hi(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = xL_hi(Col)
                      tmp_vec_lo(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = xL_lo(Col)
                      cycle
                   endif
                   js=j
                   exit
                enddo
                do j=crsMat%crsRow_ptr(k+1)-1,crsMat%crsRow_ptr(k),-1
                   if( iend<crsMat%crsCol(j) )then
                      Col=crsMat%iLs(j)
                      ! ddmul_calc1(crsMat%crsA_hi(j),0.0d0,xU_hi(Col),xU_lo(Col),c_hi,c_lo)
                      ! ddadd_calc1(c_hi,c_lo,tmp_hi,tmp_lo,tmp_hi,tmp_lo)
                      tmp_A_hi(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = crsMat%crsA_hi(j)         
                      tmp_A_lo(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = 0.0d0
                      tmp_vec_hi(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = xU_hi(Col)
                      tmp_vec_lo(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = xU_lo(Col)                      
                      cycle
                   endif
                   je=j
                   exit
                enddo
                
                do j=js,je
                   Col=crsMat%iLs(j)
                   ! ddmul_calc1(crsMat%crsA_hi(j),0.0d0,xD_hi(Col),xD_lo(Col),c_hi,c_lo)
                   ! ddadd_calc1(c_hi,c_lo,tmp_hi,tmp_lo,tmp_hi,tmp_lo)
                   tmp_A_hi(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = crsMat%crsA_hi(j)         
                   tmp_A_lo(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = 0.0d0
                   tmp_vec_hi(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = xD_hi(Col)
                   tmp_vec_lo(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = xD_lo(Col)
                enddo
                ! Ax_hi(k)=tmp_hi
                ! Ax_lo(k)=tmp_lo
             enddo
             
             ! do k=kk,kkmin
             !    Ax_hi(k)=0.0d0
             !    Ax_lo(k)=0.0d0                
             ! enddo
             do j = 1, LocalRowMax
                do k=kk,kkmin
                   ddmul_calc1(tmp_A_hi(k-kk+1,j),tmp_A_lo(k-kk+1,j),tmp_vec_hi(k-kk+1,j),tmp_vec_lo(k-kk+1,j),c_hi,c_lo)
                   ddadd_calc1(c_hi,c_lo,Ax_hi(k),Ax_lo(k),Ax_hi(k),Ax_lo(k))
                enddo
             enddo
          endif
          
       enddo
       !$OMP END PARALLEL DO
#else

   !$OMP PARALLEL DO PRIVATE(i,j,jj,k,kk,kkstep,kk0,kkmin,Ax0_lo,Ax0_hi,iLs,iUs,iLe,iUe,Row,Col, iptr_s,iptr_e,&
   !$OMP tmp_hi,tmp_lo,c_hi,c_lo,&
   !$OMP ddmul_var1, &
   !$OMP ddadd_var1 &
   !$OMP )
    do kk = 1, n, nBlock
       kk0   = kk
       kkmin = min(kk+nBlock-1,n)       

       do
          kkstep = 0
          do k = kk0,kkmin-1 
             if ( &
                  ( crsMat%crsRow_ptr_Diag_e(k)   -  crsMat%crsRow_ptr_Diag_s(k)   + 1 ) .eq. &
                  ( crsMat%crsRow_ptr_Diag_e(k+1) -  crsMat%crsRow_ptr_Diag_s(k+1) + 1 ) )then
                kkstep=kkstep+1
             else
                exit
             endif
          enddo


          do k = kk0,kk0+kkstep
!             tmp=0.0d0
             tmp_hi=0.0d0
             tmp_lo=0.0d0
             do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
                if( crsMat%crsCol(j)<istart ) then
                   Col=crsMat%iLs(j)
!                   tmp=tmp+crsMat%crsA(j)*xL(Col)
#ifdef Manual_inline
                   ddmul_calc1(crsMat%crsA_hi(j),crsMat%crsA_lo(j),xL_hi(Col),xL_lo(Col),c_hi,c_lo)
                   ddadd_calc1(c_hi,c_lo,tmp_hi,tmp_lo,tmp_hi,tmp_lo)
#else
                   call dd_mul( &
                        crsMat%crsA_hi(j), crsMat%crsA_lo(j),&
                        xL_hi(Col),xL_lo(Col),&
                        c_hi,c_lo)
                   call dd_add( &
                        c_hi,c_lo, &
                        tmp_hi,tmp_lo,&
                        tmp_hi,tmp_lo)
#endif
                   cycle
                endif
                exit
             enddo              
             do j=crsMat%crsRow_ptr(k+1)-1,crsMat%crsRow_ptr(k),-1
                if( iend<crsMat%crsCol(j) )then
                   Col=crsMat%iLs(j)
                   ! tmp=tmp+crsMat%crsA(j)*xU(Col)
#ifdef Manual_inline
                   ddmul_calc1(crsMat%crsA_hi(j),crsMat%crsA_lo(j),xU_hi(Col),xU_lo(Col),c_hi,c_lo)
                   ddadd_calc1(c_hi,c_lo,tmp_hi,tmp_lo,tmp_hi,tmp_lo)
#else
                   call dd_mul( &
                        crsMat%crsA_hi(j),crsMat%crsA_lo(j),&
                        xU_hi(Col),xU_lo(Col),&
                        c_hi,c_lo)
                   call dd_add( &
                        c_hi,c_lo, &
                        tmp_hi,tmp_lo,&
                        tmp_hi,tmp_lo)
#endif
                   cycle
                endif
                exit
             enddo
!             Ax0(k-kk+1) = tmp
             Ax0_hi(k-kk+1) = tmp_hi
             Ax0_lo(k-kk+1) = tmp_lo
          enddo
          
          do jj=1,crsMat%crsRow_ptr_Diag_e(kk0) -  crsMat%crsRow_ptr_Diag_s(kk0) + 1
             do k = kk0,kk0+kkstep
                j   = crsMat%crsRow_ptr_Diag_s(k) + jj -1
                Col=crsMat%iLs(j)
                ! Ax0(k-kk+1)=Ax0(k-kk+1)+crsMat%crsA(j)*xD(Col)
#ifdef Manual_inline
                ddmul_calc1(crsMat%crsA_hi(j),crsMat%crsA_lo(j),xD_hi(Col),xD_lo(Col),c_hi,c_lo)
                ddadd_calc1(c_hi,c_lo,Ax0_hi(k-kk+1),Ax0_lo(k-kk+1),Ax0_hi(k-kk+1),Ax0_lo(k-kk+1))
#else
                call dd_mul( &
                     crsMat%crsA_hi(j),crsMat%crsA_lo(j),&
                     xD_hi(Col),xD_lo(Col),&
                     c_hi,c_lo)
                call dd_add( &
                     c_hi,c_lo, &
                     Ax0_hi(k-kk+1),Ax0_lo(k-kk+1),&
                     Ax0_hi(k-kk+1),Ax0_lo(k-kk+1))
#endif
             enddo
          enddo
          do k = kk0,kk0+kkstep
             Ax_hi(k) = Ax0_hi(k-kk+1)
             Ax_lo(k) = Ax0_lo(k-kk+1)
          enddo

          kk0 = kk0+kkstep+1
          if(kk0>kkmin)then
             exit
          endif
       enddo

    enddo
    !$OMP END PARALLEL DO
    
#endif

    case(KSP_DOUBLE_DOUBLE_PRECISION)
#if 1
       !$OMP PARALLEL DO PRIVATE(i,j,k,kk,iLs,iUs,iLe,iUe,Row,Col,tmp_hi,&
       !$OMP tmp_lo,c_hi,c_lo,js,je, &
       !$OMP kkmin,LocalRowMax,LocalRow, &
!       !$OMP tmp_block_hi,tmp_block_lo,&
       !$OMP tmp_A_hi,tmp_A_lo,tmp_vec_hi,tmp_vec_lo, &
       !$OMP ddmul_var1, &
       !$OMP ddadd_var1) &
       !$OMP SCHEDULE(dynamic)
       do kk=1,n,tmpblocksize
          kkmin = min(kk+tmpblocksize-1,n)
          LocalRowMax = -1
          do k=kk,kkmin
             LocalRow = crsMat%crsRow_ptr(k+1) - crsMat%crsRow_ptr(k) 
             if(LocalRowMax<LocalRow)then
                LocalRowMax = LocalRow
             endif
          enddo
!          write(*,*)"LocalRowMax = ",LocalRowMax                       
          !          if(LocalRowMax>tmpblocksize)then
          if(LocalRowMax>BlockRowMax)then
             do k=kk,kkmin
                Ax_hi(k)=0.0d0
                Ax_lo(k)=0.0d0
                tmp_hi=0.0d0
                tmp_lo=0.0d0
                Row=k-1+istart
                do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
                   if( crsMat%crsCol(j)<istart ) then
                      Col=crsMat%iLs(j)
                      !Ax(k)=Ax(k)+crsMat%crsA(j)*xL(Col)
                      ! tmp=tmp+crsMat%crsA(j)*xL(Col)
                      ddmul_calc1(crsMat%crsA_hi(j),crsMat%crsA_lo(j),xL_hi(Col),xL_lo(Col),c_hi,c_lo)
                      ddadd_calc1(c_hi,c_lo,tmp_hi,tmp_lo,tmp_hi,tmp_lo)
                      cycle
                   endif
                   js=j
                   exit
                enddo
                do j=crsMat%crsRow_ptr(k+1)-1,crsMat%crsRow_ptr(k),-1
                   if( iend<crsMat%crsCol(j) )then
                      Col=crsMat%iLs(j)
                      ! tmp=tmp+crsMat%crsA(j)*xU(Col)
                      ddmul_calc1(crsMat%crsA_hi(j),crsMat%crsA_lo(j),xU_hi(Col),xU_lo(Col),c_hi,c_lo)
                      ddadd_calc1(c_hi,c_lo,tmp_hi,tmp_lo,tmp_hi,tmp_lo)
                      cycle
                   endif
                   je=j
                   exit
                enddo
                
                !       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
                do j=js,je
                   !          Col=crsMat%crsCol(j)-istart+1
                   Col=crsMat%iLs(j)
                   !Ax(k)=Ax(k)+crsMat%crsA(j)*xD(Col)
                   ! tmp=tmp+crsMat%crsA(j)*xD(Col)
                   ddmul_calc1(crsMat%crsA_hi(j),crsMat%crsA_lo(j),xD_hi(Col),xD_lo(Col),c_hi,c_lo)
                   ddadd_calc1(c_hi,c_lo,tmp_hi,tmp_lo,tmp_hi,tmp_lo)
                enddo
                Ax_hi(k)=tmp_hi
                Ax_lo(k)=tmp_lo
             enddo             
          else
             do k=kk,kkmin
                Ax_hi(k)=0.0d0
                Ax_lo(k)=0.0d0                
             enddo
             tmp_A_hi = 0.0d0
             tmp_A_lo = 0.0d0    
             tmp_vec_hi = 0.0d0
             tmp_vec_lo = 0.0d0
             
             do k=kk,kkmin
                tmp_hi=0.0d0
                tmp_lo=0.0d0
                
                Row=k-1+istart
                do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
                   if( crsMat%crsCol(j)<istart ) then
                      Col=crsMat%iLs(j)
                      ! ddmul_calc1(crsMat%crsA_hi(j),0.0d0,xL_hi(Col),xL_lo(Col),c_hi,c_lo)
                      ! ddadd_calc1(c_hi,c_lo,tmp_hi,tmp_lo,tmp_hi,tmp_lo)
                      tmp_A_hi(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = crsMat%crsA_hi(j)         
                      tmp_A_lo(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = crsMat%crsA_lo(j)         
                      tmp_vec_hi(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = xL_hi(Col)
                      tmp_vec_lo(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = xL_lo(Col)
                      cycle
                   endif
                   js=j
                   exit
                enddo
                do j=crsMat%crsRow_ptr(k+1)-1,crsMat%crsRow_ptr(k),-1
                   if( iend<crsMat%crsCol(j) )then
                      Col=crsMat%iLs(j)
                      ! ddmul_calc1(crsMat%crsA_hi(j),0.0d0,xU_hi(Col),xU_lo(Col),c_hi,c_lo)
                      ! ddadd_calc1(c_hi,c_lo,tmp_hi,tmp_lo,tmp_hi,tmp_lo)
                      tmp_A_hi(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = crsMat%crsA_hi(j)         
                      tmp_A_lo(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = crsMat%crsA_lo(j)         
                      tmp_vec_hi(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = xU_hi(Col)
                      tmp_vec_lo(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = xU_lo(Col)                      
                      cycle
                   endif
                   je=j
                   exit
                enddo
                
                do j=js,je
                   Col=crsMat%iLs(j)
                   ! ddmul_calc1(crsMat%crsA_hi(j),0.0d0,xD_hi(Col),xD_lo(Col),c_hi,c_lo)
                   ! ddadd_calc1(c_hi,c_lo,tmp_hi,tmp_lo,tmp_hi,tmp_lo)
                   tmp_A_hi(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = crsMat%crsA_hi(j)       
                   tmp_A_lo(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = crsMat%crsA_lo(j)         
                   tmp_vec_hi(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = xD_hi(Col)
                   tmp_vec_lo(k-kk+1,j-crsMat%crsRow_ptr(k)+1) = xD_lo(Col)
                enddo
                ! Ax_hi(k)=tmp_hi
                ! Ax_lo(k)=tmp_lo
             enddo
             
             ! do k=kk,kkmin
             !    Ax_hi(k)=0.0d0
             !    Ax_lo(k)=0.0d0                
             ! enddo
             do j = 1, LocalRowMax
                do k=kk,kkmin
                   ddmul_calc1(tmp_A_hi(k-kk+1,j),tmp_A_lo(k-kk+1,j),tmp_vec_hi(k-kk+1,j),tmp_vec_lo(k-kk+1,j),c_hi,c_lo)
                   ddadd_calc1(c_hi,c_lo,Ax_hi(k),Ax_lo(k),Ax_hi(k),Ax_lo(k))
                enddo
             enddo
          endif
          
       enddo
       !$OMP END PARALLEL DO
       
#else
       
   !$OMP PARALLEL DO PRIVATE(i,j,jj,k,kk,kkstep,kk0,kkmin,Ax0_lo,Ax0_hi,iLs,iUs,iLe,iUe,Row,Col, iptr_s,iptr_e,&
   !$OMP tmp_hi,tmp_lo,c_hi,c_lo,&
   !$OMP ddmul_var1, &
   !$OMP ddadd_var1 &
   !$OMP )
    do kk = 1, n, nBlock
       kk0   = kk
       kkmin = min(kk+nBlock-1,n)       

       do
          kkstep = 0
          do k = kk0,kkmin-1 
             if ( &
                  ( crsMat%crsRow_ptr_Diag_e(k)   -  crsMat%crsRow_ptr_Diag_s(k)   + 1 ) .eq. &
                  ( crsMat%crsRow_ptr_Diag_e(k+1) -  crsMat%crsRow_ptr_Diag_s(k+1) + 1 ) )then
                kkstep=kkstep+1
             else
                exit
             endif
          enddo


          do k = kk0,kk0+kkstep
!             tmp=0.0d0
             tmp_hi=0.0d0
             tmp_lo=0.0d0
             do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
                if( crsMat%crsCol(j)<istart ) then
                   Col=crsMat%iLs(j)
!                   tmp=tmp+crsMat%crsA(j)*xL(Col)
#ifdef Manual_inline
                   ddmul_calc1(crsMat%crsA_hi(j),crsMat%crsA_lo(j),xL_hi(Col),xL_lo(Col),c_hi,c_lo)
                   ddadd_calc1(c_hi,c_lo,tmp_hi,tmp_lo,tmp_hi,tmp_lo)
#else
                   call dd_mul( &
                        crsMat%crsA_hi(j), crsMat%crsA_lo(j),&
                        xL_hi(Col),xL_lo(Col),&
                        c_hi,c_lo)
                   call dd_add( &
                        c_hi,c_lo, &
                        tmp_hi,tmp_lo,&
                        tmp_hi,tmp_lo)
#endif
                   cycle
                endif
                exit
             enddo              
             do j=crsMat%crsRow_ptr(k+1)-1,crsMat%crsRow_ptr(k),-1
                if( iend<crsMat%crsCol(j) )then
                   Col=crsMat%iLs(j)
                   ! tmp=tmp+crsMat%crsA(j)*xU(Col)
#ifdef Manual_inline
                   ddmul_calc1(crsMat%crsA_hi(j),crsMat%crsA_lo(j),xU_hi(Col),xU_lo(Col),c_hi,c_lo)
                   ddadd_calc1(c_hi,c_lo,tmp_hi,tmp_lo,tmp_hi,tmp_lo)
#else
                   call dd_mul( &
                        crsMat%crsA_hi(j),crsMat%crsA_lo(j),&
                        xU_hi(Col),xU_lo(Col),&
                        c_hi,c_lo)
                   call dd_add( &
                        c_hi,c_lo, &
                        tmp_hi,tmp_lo,&
                        tmp_hi,tmp_lo)
#endif
                   cycle
                endif
                exit
             enddo
!             Ax0(k-kk+1) = tmp
             Ax0_hi(k-kk+1) = tmp_hi
             Ax0_lo(k-kk+1) = tmp_lo
          enddo
          
          do jj=1,crsMat%crsRow_ptr_Diag_e(kk0) -  crsMat%crsRow_ptr_Diag_s(kk0) + 1
             do k = kk0,kk0+kkstep
                j   = crsMat%crsRow_ptr_Diag_s(k) + jj -1
                Col=crsMat%iLs(j)
                ! Ax0(k-kk+1)=Ax0(k-kk+1)+crsMat%crsA(j)*xD(Col)
#ifdef Manual_inline
                ddmul_calc1(crsMat%crsA_hi(j),crsMat%crsA_lo(j),xD_hi(Col),xD_lo(Col),c_hi,c_lo)
                ddadd_calc1(c_hi,c_lo,Ax0_hi(k-kk+1),Ax0_lo(k-kk+1),Ax0_hi(k-kk+1),Ax0_lo(k-kk+1))
#else
                call dd_mul( &
                     crsMat%crsA_hi(j),crsMat%crsA_lo(j),&
                     xD_hi(Col),xD_lo(Col),&
                     c_hi,c_lo)
                call dd_add( &
                     c_hi,c_lo, &
                     Ax0_hi(k-kk+1),Ax0_lo(k-kk+1),&
                     Ax0_hi(k-kk+1),Ax0_lo(k-kk+1))
#endif
             enddo
          enddo
          do k = kk0,kk0+kkstep
             Ax_hi(k) = Ax0_hi(k-kk+1)
             Ax_lo(k) = Ax0_lo(k-kk+1)
          enddo

          kk0 = kk0+kkstep+1
          if(kk0>kkmin)then
             exit
          endif
       enddo

    enddo
    !$OMP END PARALLEL DO
    
#endif

    end select

!    call fapp_start("MatVecCRS",1,0)

!    call fapp_stop("MatVecCRS",1,0)
    
  end subroutine ddMatVec_crs

! todo
! 倍精度のSpMVとループの回し方がちがうので修正する
  subroutine ddMatVec_crs_ovl1(prm,crsMat,&
       xL_hi,xL_lo, &
       xD_hi,xD_lo, &
       xU_hi,xU_lo, &
       Ax_hi,Ax_lo, &
       commVec,&
       nBlock, &
       PRECISION_A)
    use matrix
    use comm
    use ddcomm_kernel
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(CRS_Mat) :: crsMat 
    type(comm_Vec) :: commVec
    real*8 xL_hi(crsMat%nColL)
    real*8 xL_lo(crsMat%nColL)
    real*8 xD_hi(prm%n)
    real*8 xD_lo(prm%n)
    real*8 xU_hi(crsMat%nColU)
    real*8 xU_lo(crsMat%nColU)
    real*8 Ax_hi(prm%n)
    real*8 Ax_lo(prm%n)
    integer nBlock
    integer PRECISION_A

    integer(kind=intType) i,j,k,iL,iU
    integer(kind=intType) n,istart,iend
    integer(kind=intType) Row,Col
    integer(kind=intType) nColL,nColU
    real*8 sbuf_hi(commVec%nrk_data_send)
    real*8 sbuf_lo(commVec%nrk_data_send)
    real*8 rbuf_hi(commVec%nrk_data)
    real*8 rbuf_lo(commVec%nrk_data)
    real*8 c_hi,c_lo

    n = prm%n
    istart=prm%istart
    iend=prm%iend
    nColL=crsMat%nColL
    nColU=crsMat%nColU

    select case (PRECISION_A)

    case(KSP_DOUBLE_PRECISION)

       call ddnrk_vec_pack( &
            commVec%nrk_prm_send, &
            commVec%n, &
            commVec%nrk_data_send, &
            xD_hi,xD_lo, &
            sbuf_hi,sbuf_lo, &
            commVec%nrk_rank_send)

       !$OMP PARALLEL 

       !$OMP MASTER
       call ddnrk_vec_comm(commVec,&
            sbuf_hi,sbuf_lo,&
            rbuf_hi,rbuf_lo)
       !$OMP END MASTER

       !$OMP DO PRIVATE(j,k,Col,c_hi,c_lo) SCHEDULE(dynamic,nBlock)
       do k=1,n
          Ax_hi(k)=0.0d0
          Ax_lo(k)=0.0d0
          do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
             if(crsMat%crsCol(j)<istart)then
                cycle
             endif
             if( iend<crsMat%crsCol(j) )then
                exit
             endif
             Col=crsMat%crsCol(j)-istart+1
             ! Ax(k)=Ax(k)+crsMat%crsA(j)*xD(Col)
             call dd_mul( &
                  crsMat%crsA_hi(j),0.0d0,&
                  xD_hi(Col),xD_lo(Col),&
                  c_hi,c_lo)
             call dd_add( &
                  c_hi,c_lo, &
                  Ax_hi(k),Ax_lo(k),&
                  Ax_hi(k),Ax_lo(k))

          enddo
       enddo

       !$OMP END PARALLEL

       call ddget_nrk_vecLU(&
            commVec, &
            rbuf_hi,rbuf_lo, &
            nColL, &
            nColU, &
            xL_hi,xL_lo, &
            xU_hi,xU_lo &
            )

       !$OMP PARALLEL DO PRIVATE(i,j,k,iL,iU,Row,Col,c_hi,c_lo) SCHEDULE(static,nBlock)
       do k=1,n
          Row=k-1+istart
          iL=1
          iU=1
          do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
             if( crsMat%crsCol(j)<istart ) then
                do i=iL,nColL
                   if( crsMat%RowL(i).eq.crsMat%crsCol(j) )then
                      Col=i
                      ! Ax(k)=Ax(k)+crsMat%crsA(j)*xL(Col)
                      call dd_mul( &
                           crsMat%crsA_hi(j),0.0d0,&
                           xL_hi(Col),xL_lo(Col),&
                           c_hi,c_lo)
                      call dd_add( &
                           c_hi,c_lo, &
                           Ax_hi(k),Ax_lo(k),&
                           Ax_hi(k),Ax_lo(k))
                      iL=i
                      exit
                   endif
                enddo
                cycle
             endif
             if( iend<crsMat%crsCol(j) )then
                do i=iU,nColU
                   if( crsMat%RowU(i).eq.crsMat%crsCol(j) )then
                      Col=i
                      ! Ax(k)=Ax(k)+crsMat%crsA(j)*xU(Col)
                      call dd_mul( &
                           crsMat%crsA_hi(j),0.0d0,&
                           xU_hi(Col),xU_lo(Col),&
                           c_hi,c_lo)
                      call dd_add( &
                           c_hi,c_lo, &
                           Ax_hi(k),Ax_lo(k),&
                           Ax_hi(k),Ax_lo(k))

                      iU=i
                      exit
                   endif
                enddo
                cycle
             endif
          enddo
       enddo
       !$OMP END PARALLEL DO

    case(KSP_DOUBLE_DOUBLE_PRECISION)

       call ddnrk_vec_pack( &
            commVec%nrk_prm_send, &
            commVec%n, &
            commVec%nrk_data_send, &
            xD_hi,xD_lo, &
            sbuf_hi,sbuf_lo, &
            commVec%nrk_rank_send)

       !$OMP PARALLEL 

       !$OMP MASTER
       call ddnrk_vec_comm(commVec,&
            sbuf_hi,sbuf_lo,&
            rbuf_hi,rbuf_lo)
       !$OMP END MASTER

       !$OMP DO PRIVATE(j,k,Col,c_hi,c_lo) SCHEDULE(dynamic,nBlock)
       do k=1,n
          Ax_hi(k)=0.0d0
          Ax_lo(k)=0.0d0
          do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
             if(crsMat%crsCol(j)<istart)then
                cycle
             endif
             if( iend<crsMat%crsCol(j) )then
                exit
             endif
             Col=crsMat%crsCol(j)-istart+1
             ! Ax(k)=Ax(k)+crsMat%crsA(j)*xD(Col)
             call dd_mul( &
                  crsMat%crsA_hi(j),crsMat%crsA_lo(j),&
                  xD_hi(Col),xD_lo(Col),&
                  c_hi,c_lo)
             call dd_add( &
                  c_hi,c_lo, &
                  Ax_hi(k),Ax_lo(k),&
                  Ax_hi(k),Ax_lo(k))

          enddo
       enddo

       !$OMP END PARALLEL

       call ddget_nrk_vecLU(&
            commVec, &
            rbuf_hi,rbuf_lo, &
            nColL, &
            nColU, &
            xL_hi,xL_lo, &
            xU_hi,xU_lo &
            )

       !$OMP PARALLEL DO PRIVATE(i,j,k,iL,iU,Row,Col,c_hi,c_lo) SCHEDULE(static,nBlock)
       do k=1,n
          Row=k-1+istart
          iL=1
          iU=1
          do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
             if( crsMat%crsCol(j)<istart ) then
                do i=iL,nColL
                   if( crsMat%RowL(i).eq.crsMat%crsCol(j) )then
                      Col=i
                      ! Ax(k)=Ax(k)+crsMat%crsA(j)*xL(Col)
                      call dd_mul( &
                           crsMat%crsA_hi(j),crsMat%crsA_lo(j),&
                           xL_hi(Col),xL_lo(Col),&
                           c_hi,c_lo)
                      call dd_add( &
                           c_hi,c_lo, &
                           Ax_hi(k),Ax_lo(k),&
                           Ax_hi(k),Ax_lo(k))
                      iL=i
                      exit
                   endif
                enddo
                cycle
             endif
             if( iend<crsMat%crsCol(j) )then
                do i=iU,nColU
                   if( crsMat%RowU(i).eq.crsMat%crsCol(j) )then
                      Col=i
                      ! Ax(k)=Ax(k)+crsMat%crsA(j)*xU(Col)
                      call dd_mul( &
                           crsMat%crsA_hi(j),crsMat%crsA_lo(j),&
                           xU_hi(Col),xU_lo(Col),&
                           c_hi,c_lo)
                      call dd_add( &
                           c_hi,c_lo, &
                           Ax_hi(k),Ax_lo(k),&
                           Ax_hi(k),Ax_lo(k))

                      iU=i
                      exit
                   endif
                enddo
                cycle
             endif
          enddo
       enddo
       !$OMP END PARALLEL DO

    end select

  end subroutine ddMatVec_crs_ovl1

! todo
! 倍精度のSpMVとループの回し方がちがうので修正する
  subroutine ddMatVec_crs_ovl2(prm,crsMat,&
       x_hi,x_lo,&
       Ax_hi,Ax_lo,&
       commVec,&
       nBlock, &
       PRECISION_A)

    use matrix
    use comm
    use comm_kernel
    use ddcomm_kernel
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(CRS_Mat) :: crsMat 
    type(comm_Vec) :: commVec
    real*8 x_hi(prm%n)
    real*8 x_lo(prm%n)

    real*8 Ax_hi(prm%n)
    real*8 Ax_lo(prm%n)
    integer nBlock

    integer(kind=intType) i,j,k,iL,iU
    integer(kind=intType) n,istart,iend
    integer(kind=intType) Row,Col
    integer(kind=intType) nColL,nColU
    real*8 sbuf_hi(commVec%nrk_data_max_send)
    real*8 sbuf_lo(commVec%nrk_data_max_send)
    real*8 rbuf_hi(commVec%nrk_data_max)
    real*8 rbuf_lo(commVec%nrk_data_max)
    real*8 rbuf0_hi(commVec%nrk_data_max)
    real*8 rbuf0_lo(commVec%nrk_data_max)

    integer irank,irankMax
    integer srank,srankMax
    integer rrank,rrankMax
    integer istart_rrank,iend_rrank  
    integer(kind=intType) offset
    integer rrank0
    integer PRECISION_A

    srankMax=commVec%nrk_rank_send
    rrankMax=commVec%nrk_rank
    irankMax=max(srankmax,rrankmax)

    n = prm%n
    istart=prm%istart
    iend=prm%iend
    nColL=crsMat%nColL
    nColU=crsMat%nColU

    rrank0 =-1
    do irank=1,irankMax+1

       if(irank>srankMax)then
          srank=-1
       else
          srank=irank
       endif
       if(irank>rrankMax)then
          rrank=-1
       else
          rrank=irank
       endif
#if 1
       call nrk_vec_pack_ovl2( &
            commVec%nrk_prm_send, &
            commVec%n, &
            commVec%nrk_data_send, &
            x_lo, &
            sbuf_lo, &
            commVec%nrk_rank_send,&
            srank)

       call nrk_vec_pack_ovl2( &
            commVec%nrk_prm_send, &
            commVec%n, &
            commVec%nrk_data_send, &
            x_hi, &
            sbuf_hi, &
            commVec%nrk_rank_send,&
            srank)
#else
       call ddnrk_vec_pack_ovl2( &
            commVec%nrk_prm_send, &
            commVec%n, &
            commVec%nrk_data_send, &
            x_hi,x_lo, &
            sbuf_hi,sbuf_lo, &
            commVec%nrk_rank_send,&
            srank)
#endif
       !$OMP PARALLEL PRIVATE(istart_rrank,iend_rrank ,j,i,k,iL,offset,Col) 

       !$OMP MASTER
#if 0
       call nrk_vec_comm_ovl2(commVec,sbuf_hi,rbuf_hi,srank,rrank)
       call nrk_vec_comm_ovl2(commVec,sbuf_lo,rbuf_lo,srank,rrank)
#else
       call ddnrk_vec_comm_ovl2(commVec,&
            sbuf_hi,sbuf_lo,&
            rbuf_hi,rbuf_lo,&
            srank,rrank)
#endif
       !$OMP END MASTER

       if(irank.eq.1)then
          !$OMP DO SCHEDULE(dynamic,nBlock)
          do k=1,n
             Ax_hi(k)=0.0d0
             Ax_lo(k)=0.0d0
             do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
                if(crsMat%crsCol(j)<istart)then
                   cycle
                endif
                if( iend<crsMat%crsCol(j) )then
                   exit
                endif
                Col=crsMat%crsCol(j)-istart+1
                ! Ax(k)=Ax(k)+crsMat%crsA(j)*x(Col)
                Ax_hi(k)=Ax_hi(k)+crsMat%crsA_hi(j)*x_hi(Col)
             enddo
          enddo
       endif

       if(rrank0 > -1)then
          istart_rrank=commVec%istart( commVec%nrk_prm(rrank0)%rank+1    )
          iend_rrank  =commVec%istart( commVec%nrk_prm(rrank0)%rank+1 +1 )-1

          !$OMP DO SCHEDULE(dynamic,nBlock)
          do k=1,n
             iL=1
             do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
                if( (istart_rrank<=crsMat%crsCol(j)).and.(crsMat%crsCol(j)<=iend_rrank) )then
                   do i=iL,commVec%nrk_prm(rrank0)%ndata
                      offset=commVec%nrk_prm(rrank0)%offset(i)
                      Col=istart_rrank+offset
                      if( Col.eq.crsMat%crsCol(j) )then
                        ! Ax(k)=Ax(k)+crsMat%crsA(j)*rbuf0(i)
                         Ax_hi(k)=Ax_hi(k)+crsMat%crsA_hi(j)*rbuf0_hi(i)
                         iL=i
                         exit
                      endif
                   enddo
                endif
                if( iend_rrank<crsMat%crsCol(j) )then
                   exit
                endif
             enddo
          enddo

       endif
       !$OMP END PARALLEL

       !$OMP PARALLEL DO PRIVATE(i)
       do i=1,commVec%nrk_data_max
          rbuf0_hi(i)=rbuf_hi(i)
          rbuf0_lo(i)=rbuf_lo(i)
       enddo
       !$OMP END PARALLEL DO
       rrank0=rrank

    enddo

! #endif
  end subroutine ddMatVec_crs_ovl2

!  subroutine MatVec_dia(prm,diaMat,x,Ax,commVec,rbuf,nBlock)
  subroutine ddMatVec_dia(prm,diaMat,&
       x_hi,x_lo,&
       Ax_hi,Ax_lo,&
       commVec,&
       rbuf_hi,rbuf_lo,&
       nBlock, &
       PRECISION_A)
    use matrix
    use comm
!    use comm_kernel
!    use vec_kernel
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(DIA_Mat) :: diaMat 
    real*8 x_hi(prm%n)
    real*8 x_lo(prm%n)
    real*8 Ax_hi(prm%n)
    real*8 Ax_lo(prm%n)
    type(comm_Vec) :: commVec

    !    real*8 buff(commVec%nrk_data)
    integer(kind=intType) i,j,k,kk
    integer(kind=intType) n,nnd,istart
    integer(kind=intType) offset
    ! nBlock --> module matrix 
    integer(kind=intType) Row,Col
    integer(kind=intType) maxRow,minRow
    integer(kind=intType) maxCol,minCol
    integer(kind=intType) i0
    integer myrank 
    integer irank,irank_start
    integer(kind=intType) kkmax,kmax,kmin
    !    real*8 rbuf(commVec%nrk_data)
    real*8 rbuf_hi(commVec%nrk_data+1)
    real*8 rbuf_lo(commVec%nrk_data+1)

    integer(kind=intType) iLs,iLe
    integer nBlock
    integer PRECISION_A
    real*8 c_hi,c_lo
    real*8 ddmul_var1
    real*8 ddadd_var1

    n      = prm%n
    nnd    = prm%nnd
    istart = prm%istart
    myrank = commVec%myrank

    select case (PRECISION_A)

    case(KSP_DOUBLE_PRECISION)
    !$OMP PARALLEL PRIVATE( &
    !$OMP kk,kkmax,j,irank,irank_start,&
    !$OMP i0,i,k,Col,Row,offset,iLs,iLe,&
    !$OMP minCol,maxCol,minRow,maxRow,kmin,kmax, &
    !$OMP c_hi,c_lo, &
    !$OMP ddmul_var1, &
    !$OMP ddadd_var1  &
    !$OMP ) 
    
    !$OMP DO SCHEDULE(dynamic,1)
    do kk=1,n,nBlock
       kkmax=min(n,kk+nBlock-1)
       do k=kk,kkmax
          Ax_hi(k)=0.0d0
          Ax_lo(k)=0.0d0
       enddo

       do j=1,nnd
          offset = diaMat%offset(j)
          minCol = kk    + istart - 1 + offset
          maxCol = kkmax + istart - 1 + offset
          minRow = 1 + istart - 1
          maxRow = n + istart - 1
          call set_Head_Tail(minRow,maxRow,minCol,maxCol,kmin,kmax)
          kmin=kmin-istart+1
          kmax=kmax-istart+1
          do k=kmin,kmax
!             Ax(k-offset) = &
!                  Ax(k-offset) + diaMat%diaA(k-offset + (j-1)*n) * x(k)
#ifdef Manual_inline
             ddmul_calc1(diaMat%diaA_hi(k-offset+(j-1)*n),0.0d0,x_hi(k),x_lo(k),c_hi,c_lo)
             ddadd_calc1(c_hi,c_lo,Ax_hi(k-offset),Ax_lo(k-offset),Ax_hi(k-offset),Ax_lo(k-offset))
#else
             call dd_mul( &
                  diaMat%diaA_hi(k-offset + (j-1)*n),0.0d0,&
                  x_hi(k),x_lo(k),&
                  c_hi,c_lo)
             call dd_add( &
                  c_hi,c_lo, &
                  Ax_hi(k-offset),Ax_lo(k-offset),&
                  Ax_hi(k-offset),Ax_lo(k-offset))
#endif
          enddo
          do k=kk,kkmax
             i = diaMat%iLs(  k + (j-1)*n  )
             ! Ax(k)=Ax(k) + diaMat%diaA(k+(j-1)*n) * rbuf(i)
#ifdef Manual_inline
             ddmul_calc1(diaMat%diaA_hi(k+(j-1)*n),0.0d0,rbuf_hi(i),rbuf_lo(i),c_hi,c_lo)
             ddadd_calc1(c_hi,c_lo,Ax_hi(k),Ax_lo(k),Ax_hi(k),Ax_lo(k))
#else
             call dd_mul( &
                  diaMat%diaA_hi(k+(j-1)*n),0.0d0,&
                  rbuf_hi(i)  ,rbuf_lo(i) ,&
                  c_hi,c_lo)
             call dd_add( &
                  c_hi,c_lo, &
                  Ax_hi(k),Ax_lo(k),&
                  Ax_hi(k),Ax_lo(k))
#endif
          enddo
       enddo

    enddo

    !$OMP END PARALLEL

    case(KSP_DOUBLE_DOUBLE_PRECISION)
    !$OMP PARALLEL PRIVATE( &
    !$OMP kk,kkmax,j,irank,irank_start,&
    !$OMP i0,i,k,Col,Row,offset,iLs,iLe,&
    !$OMP minCol,maxCol,minRow,maxRow,kmin,kmax, &
    !$OMP c_hi,c_lo, &
    !$OMP ddmul_var1, &
    !$OMP ddadd_var1  &
    !$OMP ) 
    
    !$OMP DO SCHEDULE(dynamic,1)
    do kk=1,n,nBlock
       kkmax=min(n,kk+nBlock-1)
       do k=kk,kkmax
          Ax_hi(k)=0.0d0
          Ax_lo(k)=0.0d0
       enddo

       do j=1,nnd
          offset = diaMat%offset(j)
          minCol = kk    + istart - 1 + offset
          maxCol = kkmax + istart - 1 + offset
          minRow = 1 + istart - 1
          maxRow = n + istart - 1
          call set_Head_Tail(minRow,maxRow,minCol,maxCol,kmin,kmax)
          kmin=kmin-istart+1
          kmax=kmax-istart+1
          do k=kmin,kmax
!             Ax(k-offset) = &
!                  Ax(k-offset) + diaMat%diaA(k-offset + (j-1)*n) * x(k)
#ifdef Manual_inline
             ddmul_calc1(diaMat%diaA_hi(k-offset+(j-1)*n),diaMat%diaA_lo(k-offset+(j-1)*n),x_hi(k),x_lo(k),c_hi,c_lo)
             ddadd_calc1(c_hi,c_lo,Ax_hi(k-offset),Ax_lo(k-offset),Ax_hi(k-offset),Ax_lo(k-offset))
#else
             call dd_mul( &
                  diaMat%diaA_hi(k-offset + (j-1)*n),diaMat%diaA_lo(k-offset + (j-1)*n),&
                  x_hi(k),x_lo(k),&
                  c_hi,c_lo)
             call dd_add( &
                  c_hi,c_lo, &
                  Ax_hi(k-offset),Ax_lo(k-offset),&
                  Ax_hi(k-offset),Ax_lo(k-offset))
#endif
          enddo
          do k=kk,kkmax
             i = diaMat%iLs(  k + (j-1)*n  )
             ! Ax(k)=Ax(k) + diaMat%diaA(k+(j-1)*n) * rbuf(i)
#ifdef Manual_inline
             ddmul_calc1(diaMat%diaA_hi(k+(j-1)*n),diaMat%diaA_lo(k+(j-1)*n),rbuf_hi(i),rbuf_lo(i),c_hi,c_lo)
             ddadd_calc1(c_hi,c_lo,Ax_hi(k),Ax_lo(k),Ax_hi(k),Ax_lo(k))
#else
             call dd_mul( &
                  diaMat%diaA_hi(k+(j-1)*n),diaMat%diaA_lo(k+(j-1)*n),&
                  rbuf_hi(i)  ,rbuf_lo(i) ,&
                  c_hi,c_lo)
             call dd_add( &
                  c_hi,c_lo, &
                  Ax_hi(k),Ax_lo(k),&
                  Ax_hi(k),Ax_lo(k))
#endif
          enddo
       enddo

    enddo

    !$OMP END PARALLEL
    end select


    !    stop
  end subroutine ddMatVec_dia



! 倍精度のSpMVとループの回し方がちがうので修正する
  subroutine ddMatVec_dia_ovl1(prm,diaMat,&
       x_hi,x_lo,&
       Ax_hi,Ax_lo,&
       commVec,&
       nBlock, &
       PRECISION_A)
    use matrix
    use comm
    use ddcomm_kernel
    use vec_kernel
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(DIA_Mat) :: diaMat 
    real*8 x_hi(prm%n)
    real*8 x_lo(prm%n)
    real*8 Ax_hi(prm%n)
    real*8 Ax_lo(prm%n)
    type(comm_Vec) :: commVec
    integer nBlock

    !    real*8 buff(commVec%nrk_data)
    integer(kind=intType) i,j,k,kk
    integer(kind=intType) n,nnd,istart
    integer(kind=intType) offset
    ! nBlock --> module matrix 
    integer(kind=intType) Row,Col
    integer(kind=intType) maxRow,minRow
    integer(kind=intType) maxCol,minCol
    integer(kind=intType) i0
    integer myrank 
    integer irank,irank_start
    integer(kind=intType) kkmax,kmax,kmin

    real*8 sbuf_hi(commVec%nrk_data_send)
    real*8 sbuf_lo(commVec%nrk_data_send)
    real*8 rbuf_hi(commVec%nrk_data)
    real*8 rbuf_lo(commVec%nrk_data)
    real*8 c_hi,c_lo
    integer PRECISION_A
    n      = prm%n
    nnd    = prm%nnd
    istart = prm%istart
    myrank = commVec%myrank

    call ddnrk_vec_pack( &
         commVec%nrk_prm_send, &
         commVec%n, &
         commVec%nrk_data_send, &
         x_hi,x_lo, &
         sbuf_hi,sbuf_lo, &
         commVec%nrk_rank_send)
    select case (PRECISION_A)
    case(KSP_DOUBLE_PRECISION)

       !$OMP PARALLEL PRIVATE( &
       !$OMP kk,kkmax,j, &
       !$OMP offset,minCol,maxCol,minRow,maxRow, &
       !$OMP kmin,kmax,k, &
       !$OMP c_hi,c_lo &
       !$OMP ) 

       !$OMP MASTER
       call ddnrk_vec_comm(commVec,sbuf_hi,sbuf_lo,rbuf_hi,rbuf_lo)
       !$OMP END MASTER

       !$OMP DO SCHEDULE(dynamic,1)
       do kk=1,n,nBlock
          kkmax=min(n,kk+nBlock-1)
          do k=kk,kkmax
             Ax_hi(k)=0.0d0
             Ax_lo(k)=0.0d0
          enddo
          do j=1,nnd
             offset = diaMat%offset(j)
             minCol = kk    + istart - 1 + offset
             maxCol = kkmax + istart - 1 + offset
             minRow = 1 + istart - 1
             maxRow = n + istart - 1
             call set_Head_Tail(minRow,maxRow,minCol,maxCol,kmin,kmax)
             ! if(myrank.eq.0)then
             !    write(*,*)minCol,maxCol,minRow,maxRow,kmin,kmax
             ! endif
             kmin=kmin-istart+1
             kmax=kmax-istart+1
             do k=kmin,kmax
                ! if(myrank.eq.1)then
                !    write(*,*)k - offset +(j-1)*n,k
                ! endif
                ! Ax(k-offset)=Ax(k-offset) + diaMat%diaA(k-offset + (j-1)*n) * x(k)
                call dd_mul( &
                     diaMat%diaA_hi(k-offset + (j-1)*n), 0.0d0,&
                     x_hi(k)  ,x_lo(k) ,&
                     c_hi,c_lo)
                call dd_add( &
                     c_hi,c_lo, &
                     Ax_hi(k-offset),Ax_lo(k-offset),&
                     Ax_hi(k-offset),Ax_lo(k-offset))

             enddo
          enddo
       enddo
       !$OMP END PARALLEL

       !$OMP PARALLEL PRIVATE( &
       !$OMP kk,kkmax,j,irank,irank_start,&
       !$OMP i0,i,k,Col,Row,offset,&
       !$OMP minCol,maxCol,minRow,maxRow,kmin,kmax, &
       !$OMP c_hi,c_lo ) 

       !$OMP DO SCHEDULE(static,1)
       do kk=1,n,nBlock
          kkmax=min(n,kk+nBlock-1)
          !       write(*,*)myrank," @ @ ",kk,kkmax
          do j=1,nnd
             offset = diaMat%offset(j)
             do irank=1,commVec%nrk_rank
                i0=0
                do i=1,irank-1
                   i0=i0+commVec%nrk_prm(i)%ndata
                enddo
                irank_start=commVec%istart( commVec%nrk_prm(irank)%rank+1 )
                minCol = kk    + istart - 1 + offset
                maxCol = kkmax + istart - 1 + offset
                minRow = irank_start
                maxRow = commVec%istart( commVec%nrk_prm(irank)%rank+1 +1 ) - 1
                call set_Head_Tail(minRow,maxRow,minCol,maxCol,kmin,kmax)
                kmin=kmin-istart+1-offset
                kmax=kmax-istart+1-offset
                do k=kmin,kmax
                   Col = k + istart+diaMat%offset(j) - 1
                   do i=1,commVec%nrk_prm(irank)%ndata
                      Row = &
                           irank_start + &
                           commVec%nrk_prm(irank)%offset(i) 
                      if(Row.eq.Col)then
                         ! Ax(k)=Ax(k) + diaMat%diaA(k+(j-1)*n) * rbuf(i+i0)
                         call dd_mul( &
                              diaMat%diaA_hi(k+(j-1)*n), 0.0d0,&
                              rbuf_hi(i+i0),rbuf_lo(i+i0) ,&
                              c_hi,c_lo)
                         call dd_add( &
                              c_hi,c_lo, &
                              Ax_hi(k),Ax_lo(k),&
                              Ax_hi(k),Ax_lo(k))

                      endif
                   enddo
                enddo

             enddo
          enddo
       enddo
       !$OMP END PARALLEL

    case(KSP_DOUBLE_DOUBLE_PRECISION)

       !$OMP PARALLEL PRIVATE( &
       !$OMP kk,kkmax,j, &
       !$OMP offset,minCol,maxCol,minRow,maxRow, &
       !$OMP kmin,kmax,k, &
       !$OMP c_hi,c_lo &
       !$OMP ) 

       !$OMP MASTER
       call ddnrk_vec_comm(commVec,sbuf_hi,sbuf_lo,rbuf_hi,rbuf_lo)
       !$OMP END MASTER

       !$OMP DO SCHEDULE(dynamic,1)
       do kk=1,n,nBlock
          kkmax=min(n,kk+nBlock-1)
          do k=kk,kkmax
             Ax_hi(k)=0.0d0
             Ax_lo(k)=0.0d0
          enddo
          do j=1,nnd
             offset = diaMat%offset(j)
             minCol = kk    + istart - 1 + offset
             maxCol = kkmax + istart - 1 + offset
             minRow = 1 + istart - 1
             maxRow = n + istart - 1
             call set_Head_Tail(minRow,maxRow,minCol,maxCol,kmin,kmax)
             ! if(myrank.eq.0)then
             !    write(*,*)minCol,maxCol,minRow,maxRow,kmin,kmax
             ! endif
             kmin=kmin-istart+1
             kmax=kmax-istart+1
             do k=kmin,kmax
                ! if(myrank.eq.1)then
                !    write(*,*)k - offset +(j-1)*n,k
                ! endif
                ! Ax(k-offset)=Ax(k-offset) + diaMat%diaA(k-offset + (j-1)*n) * x(k)
                call dd_mul( &
                     diaMat%diaA_hi(k-offset + (j-1)*n), diaMat%diaA_lo(k-offset + (j-1)*n),&
                     x_hi(k)  ,x_lo(k) ,&
                     c_hi,c_lo)
                call dd_add( &
                     c_hi,c_lo, &
                     Ax_hi(k-offset),Ax_lo(k-offset),&
                     Ax_hi(k-offset),Ax_lo(k-offset))

             enddo
          enddo
       enddo
       !$OMP END PARALLEL

       !$OMP PARALLEL PRIVATE( &
       !$OMP kk,kkmax,j,irank,irank_start,&
       !$OMP i0,i,k,Col,Row,offset,&
       !$OMP minCol,maxCol,minRow,maxRow,kmin,kmax, &
       !$OMP c_hi,c_lo ) 

       !$OMP DO SCHEDULE(static,1)
       do kk=1,n,nBlock
          kkmax=min(n,kk+nBlock-1)
          !       write(*,*)myrank," @ @ ",kk,kkmax
          do j=1,nnd
             offset = diaMat%offset(j)
             do irank=1,commVec%nrk_rank
                i0=0
                do i=1,irank-1
                   i0=i0+commVec%nrk_prm(i)%ndata
                enddo
                irank_start=commVec%istart( commVec%nrk_prm(irank)%rank+1 )
                minCol = kk    + istart - 1 + offset
                maxCol = kkmax + istart - 1 + offset
                minRow = irank_start
                maxRow = commVec%istart( commVec%nrk_prm(irank)%rank+1 +1 ) - 1
                call set_Head_Tail(minRow,maxRow,minCol,maxCol,kmin,kmax)
                kmin=kmin-istart+1-offset
                kmax=kmax-istart+1-offset
                do k=kmin,kmax
                   Col = k + istart+diaMat%offset(j) - 1
                   do i=1,commVec%nrk_prm(irank)%ndata
                      Row = &
                           irank_start + &
                           commVec%nrk_prm(irank)%offset(i) 
                      if(Row.eq.Col)then
                         ! Ax(k)=Ax(k) + diaMat%diaA(k+(j-1)*n) * rbuf(i+i0)
                         call dd_mul( &
                              diaMat%diaA_hi(k+(j-1)*n), diaMat%diaA_lo(k+(j-1)*n),&
                              rbuf_hi(i+i0),rbuf_lo(i+i0) ,&
                              c_hi,c_lo)
                         call dd_add( &
                              c_hi,c_lo, &
                              Ax_hi(k),Ax_lo(k),&
                              Ax_hi(k),Ax_lo(k))

                      endif
                   enddo
                enddo

             enddo
          enddo
       enddo
       !$OMP END PARALLEL
    end select

  end subroutine ddMatVec_dia_ovl1


! 倍精度のSpMVとループの回し方がちがうので修正する
! waitで通信完了処理をしているため問題によってはデッドロックの可能性がある
  subroutine ddMatVec_dia_ovl2(prm,diaMat,&
       x_hi,x_lo,&
       Ax_hi,Ax_lo,&
       commVec,&
       nBlock, &
       PRECISION_A)

    use matrix
    use comm
    use ddcomm_kernel
    use vec_kernel
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(DIA_Mat) :: diaMat 
    real*8 x_hi(prm%n)
    real*8 x_lo(prm%n)
    real*8 Ax_hi(prm%n)
    real*8 Ax_lo(prm%n)
    type(comm_Vec) :: commVec
    integer nBlock
    real*8 sbuf_hi(commVec%nrk_data_max_send)
    real*8 sbuf_lo(commVec%nrk_data_max_send)
    real*8 rbuf_hi(commVec%nrk_data_max)
    real*8 rbuf_lo(commVec%nrk_data_max)
    real*8 rbuf0_hi(commVec%nrk_data_max)
    real*8 rbuf0_lo(commVec%nrk_data_max)
    real*8 buff_hi(commVec%nrk_data)
    real*8 buff_lo(commVec%nrk_data)
    real*8 c_hi,c_lo

    integer(kind=intType) i,j,k,kk
    integer(kind=intType) n,nnd,istart
    integer(kind=intType) offset

    integer(kind=intType) Row,Col
    integer(kind=intType) maxRow,minRow
    integer(kind=intType) maxCol,minCol
    integer(kind=intType) i0
    integer myrank 
    integer irank_start
    integer(kind=intType) kkmax,kmax,kmin

    integer irank,irankMax
    integer srank,srankMax
    integer rrank,rrankMax
    integer rrank0
    integer PRECISION_A
    srankMax=commVec%nrk_rank_send
    rrankMax=commVec%nrk_rank
    irankMax=max(srankmax,rrankmax)

    n      = prm%n
    nnd    = prm%nnd
    istart = prm%istart
    myrank = commVec%myrank

    rrank0 =-1
    do irank=1,irankMax+1
       if(irank>srankMax)then
          srank=-1
       else
          srank=irank
       endif
       if(irank>rrankMax)then
          rrank=-1
       else
          rrank=irank
       endif

       call ddnrk_vec_pack_ovl2( &
            commVec%nrk_prm_send, &
            commVec%n, &
            commVec%nrk_data_send, &
            x_hi,x_lo, &
            sbuf_hi,sbuf_lo, &
            commVec%nrk_rank_send,&
            srank)

       select case (PRECISION_A)

       case(KSP_DOUBLE_PRECISION)

          !$OMP PARALLEL PRIVATE( &
          !$OMP kk,kkmax,j,irank_start,&
          !$OMP i0,i,k,Col,Row,offset,&
          !$OMP minCol,maxCol,minRow,maxRow,kmin,kmax, &
          !$OMP c_hi,c_lo) 

          !$OMP MASTER
          call ddnrk_vec_comm_ovl2(commVec,&
               sbuf_hi,sbuf_lo,&
               rbuf_hi,rbuf_lo,&
               srank,rrank)
          !$OMP END MASTER

          if(irank.eq.1)then
             !$OMP DO SCHEDULE(dynamic,1)
             do kk=1,n,nBlock
                kkmax=min(n,kk+nBlock-1)
                do k=kk,kkmax
                   Ax_hi(k)=0.0d0
                   Ax_lo(k)=0.0d0
                enddo
                do j=1,nnd
                   offset = diaMat%offset(j)
                   minCol = kk    + istart - 1 + offset
                   maxCol = kkmax + istart - 1 + offset
                   minRow = 1 + istart - 1
                   maxRow = n + istart - 1
                   call set_Head_Tail(minRow,maxRow,minCol,maxCol,kmin,kmax)
                   kmin=kmin-istart+1
                   kmax=kmax-istart+1
                   do k=kmin,kmax
                      ! Ax(k-offset)=Ax(k-offset) + diaMat%diaA(k-offset + (j-1)*n) * x(k)
                      call dd_mul( &
                           diaMat%diaA_hi(k-offset + (j-1)*n),0.0d0,&
                           x_hi(k),x_lo(k) ,&
                           c_hi,c_lo)
                      call dd_add( &
                           c_hi,c_lo, &
                           Ax_hi(k-offset),Ax_lo(k-offset),&
                           Ax_hi(k-offset),Ax_lo(k-offset))
                   enddo
                enddo
             enddo
          endif

          if(rrank0 > -1)then
             !$OMP DO SCHEDULE(dynamic,1)
             do kk=1,n,nBlock
                kkmax=min(n,kk+nBlock-1)
                do j=1,nnd
                   offset = diaMat%offset(j)

                   i0=0
                   do i=1,rrank0-1
                      i0=i0+commVec%nrk_prm(i)%ndata
                   enddo
                   irank_start=commVec%istart( commVec%nrk_prm(rrank0)%rank+1 )
                   minCol = kk    + istart - 1 + offset
                   maxCol = kkmax + istart - 1 + offset
                   minRow = irank_start
                   maxRow = commVec%istart( commVec%nrk_prm(rrank0)%rank+1 +1 ) - 1
                   call set_Head_Tail(minRow,maxRow,minCol,maxCol,kmin,kmax)
                   kmin=kmin-istart+1-offset
                   kmax=kmax-istart+1-offset
                   do k=kmin,kmax
                      Col = k + istart+diaMat%offset(j) - 1
                      do i=1,commVec%nrk_prm(rrank0)%ndata
                         Row = &
                              irank_start + &
                              commVec%nrk_prm(rrank0)%offset(i) 
                         if(Row.eq.Col)then
                            ! Ax(k)=Ax(k) + diaMat%diaA(k+(j-1)*n) * rbuf0(i)

                            call dd_mul( &
                                 diaMat%diaA_hi(k+(j-1)*n),0.0d0,&
                                 rbuf0_hi(i),rbuf0_lo(i) ,&
                                 c_hi,c_lo)
                            call dd_add( &
                                 c_hi,c_lo, &
                                 Ax_hi(k),Ax_lo(k),&
                                 Ax_hi(k),Ax_lo(k))

                         endif
                      enddo
                   enddo
                enddo
             enddo
          endif

          !$OMP END PARALLEL

       case(KSP_DOUBLE_DOUBLE_PRECISION)

          !$OMP PARALLEL PRIVATE( &
          !$OMP kk,kkmax,j,irank_start,&
          !$OMP i0,i,k,Col,Row,offset,&
          !$OMP minCol,maxCol,minRow,maxRow,kmin,kmax, &
          !$OMP c_hi,c_lo) 

          !$OMP MASTER
          call ddnrk_vec_comm_ovl2(commVec,&
               sbuf_hi,sbuf_lo,&
               rbuf_hi,rbuf_lo,&
               srank,rrank)
          !$OMP END MASTER

          if(irank.eq.1)then
             !$OMP DO SCHEDULE(dynamic,1)
             do kk=1,n,nBlock
                kkmax=min(n,kk+nBlock-1)
                do k=kk,kkmax
                   Ax_hi(k)=0.0d0
                   Ax_lo(k)=0.0d0
                enddo
                do j=1,nnd
                   offset = diaMat%offset(j)
                   minCol = kk    + istart - 1 + offset
                   maxCol = kkmax + istart - 1 + offset
                   minRow = 1 + istart - 1
                   maxRow = n + istart - 1
                   call set_Head_Tail(minRow,maxRow,minCol,maxCol,kmin,kmax)
                   kmin=kmin-istart+1
                   kmax=kmax-istart+1
                   do k=kmin,kmax
                      ! Ax(k-offset)=Ax(k-offset) + diaMat%diaA(k-offset + (j-1)*n) * x(k)
                      call dd_mul( &
                           diaMat%diaA_hi(k-offset + (j-1)*n), diaMat%diaA_lo(k-offset + (j-1)*n),&
                           x_hi(k),x_lo(k) ,&
                           c_hi,c_lo)
                      call dd_add( &
                           c_hi,c_lo, &
                           Ax_hi(k-offset),Ax_lo(k-offset),&
                           Ax_hi(k-offset),Ax_lo(k-offset))
                   enddo
                enddo
             enddo
          endif

          if(rrank0 > -1)then
             !$OMP DO SCHEDULE(dynamic,1)
             do kk=1,n,nBlock
                kkmax=min(n,kk+nBlock-1)
                do j=1,nnd
                   offset = diaMat%offset(j)

                   i0=0
                   do i=1,rrank0-1
                      i0=i0+commVec%nrk_prm(i)%ndata
                   enddo
                   irank_start=commVec%istart( commVec%nrk_prm(rrank0)%rank+1 )
                   minCol = kk    + istart - 1 + offset
                   maxCol = kkmax + istart - 1 + offset
                   minRow = irank_start
                   maxRow = commVec%istart( commVec%nrk_prm(rrank0)%rank+1 +1 ) - 1
                   call set_Head_Tail(minRow,maxRow,minCol,maxCol,kmin,kmax)
                   kmin=kmin-istart+1-offset
                   kmax=kmax-istart+1-offset
                   do k=kmin,kmax
                      Col = k + istart+diaMat%offset(j) - 1
                      do i=1,commVec%nrk_prm(rrank0)%ndata
                         Row = &
                              irank_start + &
                              commVec%nrk_prm(rrank0)%offset(i) 
                         if(Row.eq.Col)then
                            ! Ax(k)=Ax(k) + diaMat%diaA(k+(j-1)*n) * rbuf0(i)

                            call dd_mul( &
                                 diaMat%diaA_hi(k+(j-1)*n), diaMat%diaA_lo(k+(j-1)*n),&
                                 rbuf0_hi(i),rbuf0_lo(i) ,&
                                 c_hi,c_lo)
                            call dd_add( &
                                 c_hi,c_lo, &
                                 Ax_hi(k),Ax_lo(k),&
                                 Ax_hi(k),Ax_lo(k))

                         endif
                      enddo
                   enddo
                enddo
             enddo
          endif

          !$OMP END PARALLEL
       end select

       !$OMP PARALLEL DO PRIVATE(i)
       do i=1,commVec%nrk_data_max
          rbuf0_hi(i)=rbuf_hi(i)
          rbuf0_lo(i)=rbuf_lo(i)
       enddo
       !$OMP END PARALLEL DO
       rrank0=rrank

    enddo

  end subroutine ddMatVec_dia_ovl2

  ! z(:)=y(:)*x(:)
  subroutine ddpointJacobi(n,&
       x_hi,x_lo, &
       y_hi,y_lo, &
       z_hi,z_lo, &
       PRECISION_PRECONDITION)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 x_hi(n)
    real*8 x_lo(n)
    real*8 y_hi(n)
    real*8 y_lo(n)
    real*8 z_hi(n)
    real*8 z_lo(n)
    integer PRECISION_PRECONDITION
    integer(kind=intType) i 
    real*8 ddmul_var1

    select case (PRECISION_PRECONDITION)
    case(KSP_DOUBLE_PRECISION)
       !$OMP PARALLEL DO PRIVATE(i,ddmul_var1) 
       do i=1,n
          !          z(i)=y(i)*x(i)
#ifdef Manual_inline
          ddmul_calc1(y_hi(i),0.0d0,x_hi(i),x_lo(i),z_hi(i),z_lo(i))
#else
          call dd_mul(y_hi(i),0.0d0,x_hi(i),x_lo(i),z_hi(i),z_lo(i))
#endif
       enddo
       !$OMP END PARALLEL DO

    case(KSP_DOUBLE_DOUBLE_PRECISION)
       !$OMP PARALLEL DO PRIVATE(i,ddmul_var1) 
       do i=1,n
          !          z(i)=y(i)*x(i)
#ifdef Manual_inline
          ddmul_calc1(y_hi(i),y_lo(i),x_hi(i),x_lo(i),z_hi(i),z_lo(i))
#else
          call dd_mul(y_hi(i),y_lo(i),x_hi(i),x_lo(i),z_hi(i),z_lo(i))
#endif
       enddo
       !$OMP END PARALLEL DO      
    end select

  end subroutine ddpointJacobi

  subroutine ddblockJacobi_crs( &
       method, &
       r_hi,r_lo, &
       s_hi,s_lo )
    use krylov,only : KSP
    use omp_lib
    use matrix
    use PARCEL_integer
    implicit none
    type(KSP) :: method
    real*8 r_hi(method%prm%n)
    real*8 r_lo(method%prm%n)
    real*8 s_hi(method%prm%n)
    real*8 s_lo(method%prm%n)
    integer independ_nvec
    real*8 z
    integer PRECISION_PRECONDITION
    integer(kind=intType) i,j,k
    integer(kind=intType) Row,Col
    integer(kind=intType) n
    integer(kind=intType) istart,iend

    integer  numth,myth,max_numth
    integer(kind=intType) n_th
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th
    integer(kind=intType) ntmp

    n=method%prm%n
    istart=method%prm%istart
    iend=method%prm%iend

    ! -----
    ! ブロック数よりもスレッド数が大きい場合はスレッド数を制限する
    call get_maxthreads(max_numth)
    call set_precon_maxthreads(n,method%precon_thblock)

!    call set_precon_maxthreads0(n)

    ! -----
    ! (D+L)D^-1(D+U)s=r

    ! -----
    !   z = D^-1(D+U)s
    !   (D+L)z=r
    !$OMP PARALLEL &
    !$OMP PRIVATE( &
    !$OMP numth,myth,n_startth,n_endth,istart_th,iend_th,ntmp &
    !$OMP )
    numth = omp_get_num_threads()
    myth  = omp_get_thread_num()
    ntmp=1
    call divide_region( &
         myth,numth,ntmp,n,n_startth,n_endth)
    istart_th = istart -1 + n_startth 
    iend_th   = istart_th + (n_endth - n_startth) 
    call ddsolve_ILU_crs( &
         n,method%crsMat, &
         istart,iend, &
         istart_th,iend_th, &
         n_startth,n_endth, &
         method%Dinv_hi(n_startth), &
         method%Dinv_lo(n_startth), &
         r_hi(n_startth),r_lo(n_startth), &
         s_hi(n_startth),s_lo(n_startth), &
         method%ILUregionMat(myth+1)%ILU_Lstart, &
         method%ILUregionMat(myth+1)%ILU_Lend,   &
         method%ILUregionMat(myth+1)%ILU_Ustart, &
         method%ILUregionMat(myth+1)%ILU_Uend,   &
         method%independ_nvec, &
         method%PRECISION_PRECONDITION)
    !$OMP END PARALLEL

!    if(n<=max_numth)then
    call omp_set_num_threads(max_numth)
!    endif

  end subroutine ddblockJacobi_crs

  ! subroutine ddblockJacobi_dia(&
  !      prm,diaMat,&
  !      Dinv_hi,Dinv_lo, &
  !      r_hi,r_lo, &
  !      s_hi,s_lo, &
  !      PRECISION_PRECONDITION)
  subroutine ddblockJacobi_dia(&
       method,&
       r_hi,r_lo,s_hi,s_lo)
    use krylov,only :KSP
    use omp_lib
    use matrix
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) istart,iend
    ! type(prmMat)  :: prm
    ! type(DIA_Mat) :: diaMat
    ! real*8 Dinv_hi(prm%n)
    ! real*8 Dinv_lo(prm%n)
    type(KSP) :: method
    real*8 r_hi(method%prm%n)
    real*8 r_lo(method%prm%n)
    real*8 s_hi(method%prm%n)
    real*8 s_lo(method%prm%n)
    ! real*8 z

!    integer PRECISION_PRECONDITION
    integer(kind=intType) i,j,k
    integer(kind=intType) Row,Col

    integer numth,myth,max_numth
    integer(kind=intType) n_th
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th
    integer(kind=intType) nnd
    integer(kind=intType) ntmp
    integer precon_thblock

    n=method%prm%n
    istart=method%prm%istart
    iend=method%prm%iend
    nnd=method%prm%nnd 
    precon_thblock=method%precon_thblock 

    ! -----
    ! ブロック数よりもスレッド数が大きい場合はスレッド数を制限する
    call get_maxthreads(max_numth)
    call set_precon_maxthreads(n,precon_thblock)
!    call set_precon_maxthreads0(n)

    !$OMP PARALLEL &
    !$OMP PRIVATE( &
    !$OMP numth,myth,n_startth,n_endth,istart_th,iend_th,ntmp &
    !$OMP )
    numth = omp_get_num_threads()
    myth  = omp_get_thread_num()
    ntmp=1
    call divide_region( &
         myth,numth,ntmp,n,n_startth,n_endth)
    istart_th = istart -1 + n_startth 
    iend_th   = istart_th + (n_endth - n_startth) 
    call ddsolve_ILU_dia( &
         n,nnd, &
         method%diaMat%diaA_hi,method%diaMat%diaA_lo, &
         method%diaMat%offset, &
         istart,iend, &
         istart_th,iend_th, &
         n_startth,n_endth, &
         method%Dinv_hi(n_startth), &
         method%Dinv_lo(n_startth), &
         r_hi(n_startth), &
         r_lo(n_startth), &
         s_hi(n_startth), &
         s_lo(n_startth), &
         method%ILUregionMat(myth+1)%ILU_Lstart, &
         method%ILUregionMat(myth+1)%ILU_Lend,   &
         method%ILUregionMat(myth+1)%ILU_Ustart, &
         method%ILUregionMat(myth+1)%ILU_Uend,   &
         method%independ_nvec,    &
         method%PRECISION_PRECONDITION)         
    !$OMP END PARALLEL

!    if(n<=max_numth)then
    call omp_set_num_threads(max_numth)
!    endif

  end subroutine ddblockJacobi_dia

  subroutine ddsolve_ILU_crs( &
       n,crsMat, &
       istart,iend, &
       istart_th,iend_th, &
       n_startth,n_endth, &
       Dinv_hi,Dinv_lo,&
       r_hi,r_lo,&
       s_hi,s_lo, &
       ILU_Lstart, &
       ILU_Lend,   &
       ILU_Ustart, &
       ILU_Uend,    &
       independ_nvec, &
       PRECISION_PRECONDITION)
    use matrix
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    type(CRS_Mat) :: crsMat

    integer(kind=intType) istart,iend
    integer(kind=intType) istart_th,iend_th
    integer(kind=intType) n_startth,n_endth

    real*8 r_hi(n_startth:n_endth)
    real*8 r_lo(n_startth:n_endth)
    real*8 s_hi(n_startth:n_endth)
    real*8 s_lo(n_startth:n_endth)
    real*8 Dinv_hi(n_startth:n_endth)
    real*8 Dinv_lo(n_startth:n_endth)
    integer(kind=intType) ILU_Lstart(n_startth:n_endth)
    integer(kind=intType) ILU_Lend(n_startth:n_endth)
    integer(kind=intType) ILU_Ustart(n_startth:n_endth)
    integer(kind=intType) ILU_Uend(n_startth:n_endth)

    integer independ_nvec
    integer PRECISION_PRECONDITION 
    integer(kind=intType) j,k
    integer(kind=intType) kk,kk0,kkmin,kkstep,jj
    integer(kind=intType) Row,Col

    real*8 z_hi(independ_nvec),z_lo(independ_nvec)
    real*8 zz_hi,zz_lo
    real*8 c_hi,c_lo
    real*8 ddmul_var1
    real*8 ddsub_var1
    real*8 ddadd_var1

    select case (PRECISION_PRECONDITION)

    case(KSP_DOUBLE_PRECISION)
       if(independ_nvec.eq.1)then
          ! -----
          !   z = D^-1(D+U)s
          !   (D+L)z=r
          do k=n_startth,n_endth
             zz_hi = r_hi(k)
             zz_lo = r_lo(k)
             do j=ILU_Lstart(k),ILU_Lend(k)
                Col=crsMat%crsCol(j)-istart+1 
                ! r - Lz の計算
                ! zz=zz-s(Col)*crsA(j)
#ifdef Manual_inline
                ddmul_calc1(s_hi(Col),s_lo(Col),crsMat%crsA_hi(j),0.0d0,c_hi,c_lo)
                ddsub_calc1(zz_hi,zz_lo,c_hi,c_lo,zz_hi,zz_lo)
#else
                call dd_mul( &
                     s_hi(Col),s_lo(Col),&
                     crsMat%crsA_hi(j),0.0d0, &
                     c_hi,c_lo)
                call dd_sub( &
                     zz_hi,zz_lo, &
                     c_hi,c_lo, &
                     zz_hi,zz_lo)
#endif
             enddo
             ! z = D^-1 (r - Lz ) 
             ! s(k)=Dinv(k)*zz
#ifdef Manual_inline
             ddmul_calc1(zz_hi,zz_lo,Dinv_hi(k),0.0d0,s_hi(k),s_lo(k))
#else
             call dd_mul( &
                  zz_hi   ,zz_lo,&
                  Dinv_hi(k),0.0d0, &
                  s_hi(k)   ,s_lo(k) )             
#endif
          enddo

          ! -----
          ! D^-1(D+U)s=z
          !   (I+D^-1U)s=z 
          !   s=z - D^-1 Us       
          do k=n_endth,n_startth,-1
             zz_hi=0.0d0
             zz_lo=0.0d0
             do j=ILU_Ustart(k),ILU_Uend(k),-1 
                Col=crsMat%crsCol(j)-istart+1 
                ! Us の計算
                ! zz=zz+s(Col)*crsA(j)
#ifdef Manual_inline
                ddmul_calc1(s_hi(Col),s_lo(Col),crsMat%crsA_hi(j),0.0d0,c_hi,c_lo)
                ddadd_calc1(zz_hi,zz_lo,c_hi,c_lo,zz_hi,zz_lo)

#else
                call dd_mul( &
                     s_hi(Col)        ,s_lo(Col),&
                     crsMat%crsA_hi(j),0.0d0, &
                     c_hi,c_lo)
                call dd_add( &
                     zz_hi,zz_lo, &
                     c_hi,c_lo,&
                     zz_hi,zz_lo)                
#endif
             enddo
             ! s(k)=s(k)-Dinv(k)*zz
#ifdef Manual_inline
             ddmul_calc1(Dinv_hi(k),0.0d0,zz_hi,zz_lo,c_hi,c_lo)
             ddsub_calc1(s_hi(k),s_lo(k),c_hi,c_lo,s_hi(k),s_lo(k))
#else
             call dd_mul( &
                  Dinv_hi(k),0.0d0, &
                  zz_hi,zz_lo, &
                  c_hi,c_lo)       
             call dd_sub(&
                  s_hi(k),s_lo(k),&
                  c_hi,c_lo,&
                  s_hi(k),s_lo(k))
#endif
          enddo

       else

          ! -----
          !   z = D^-1(D+U)s
          !   (D+L)z=r
          do kk=n_startth,n_endth,independ_nvec
             do k=kk,min(kk+independ_nvec-1,n_endth)
                z_hi(k-kk+1)=r_hi(k)
                z_lo(k-kk+1)=r_lo(k)
             enddo
             kk0 = kk
             kkmin = min(kk+independ_nvec-1,n_endth)
             do
                kkstep = 0
                do k = kk0,kkmin-1 
                   if ( ( ILU_Lend(k)-ILU_Lstart(k) + 1 ).eq.( ILU_Lend(k+1)-ILU_Lstart(k+1) + 1 ) )then
                      kkstep=kkstep+1
                   else
                      exit
                   endif
                enddo

                do jj=1,(ILU_Lend(kk0)-ILU_Lstart(kk0)+1)
                   do k = kk0,kk0+kkstep
                      j = ILU_Lstart(k)+jj-1                
                      Col = crsMat%crsCol(j)-istart+1 
                      ! r - Lz の計算
                      !  z(k-kk+1)=z(k-kk+1)-s(Col)*crsA(j)
#ifdef Manual_inline
                      ddmul_calc1(s_hi(Col),s_lo(Col),crsMat%crsA_hi(j),0.0d0,c_hi,c_lo)
                      ddsub_calc1(z_hi(k-kk+1),z_lo(k-kk+1),c_hi,c_lo,z_hi(k-kk+1),z_lo(k-kk+1))
#else
                      call dd_mul( &
                           s_hi(Col),s_lo(Col),&
                           crsMat%crsA_hi(j),0.0d0, &
                           c_hi,c_lo)
                      call dd_sub( &
                           z_hi(k-kk+1),z_lo(k-kk+1), &
                           c_hi,c_lo, &
                           z_hi(k-kk+1),z_lo(k-kk+1))
#endif
                   enddo
                enddo
                kk0 = kk0+kkstep+1
                if(kk0>kkmin)then
                   exit
                endif
             enddo

             do k=kk,min(kk+independ_nvec-1,n_endth)
                ! s(k)=Dinv(k)*z(k-kk+1) 
#ifdef Manual_inline
                ddmul_calc1(z_hi(k-kk+1),z_lo(k-kk+1),Dinv_hi(k),0.0d0,s_hi(k),s_lo(k))
#else
                call dd_mul( &
                     z_hi(k-kk+1)   ,z_lo(k-kk+1),&
                     Dinv_hi(k),0.0d0, &
                     s_hi(k)   ,s_lo(k) )
#endif
             enddo
          enddo
          ! -----
          ! D^-1(D+U)s=z
          !   (I+D^-1U)s=z 
          !   s=z - D^-1 Us
          do kk=n_endth,n_startth,-independ_nvec
             do k = kk,max(kk-independ_nvec+1,n_startth),-1
                z_hi(kk-k+1)=0.0d0
                z_lo(kk-k+1)=0.0d0
             enddo
             kk0 = kk
             kkmin = max(kk-independ_nvec+1,n_startth)
             do 
                kkstep = 0
                do k = kk0,kkmin+1,-1 
                   if ( ( ILU_Ustart(k)-ILU_Uend(k) + 1 ).eq.( ILU_Ustart(k-1) - ILU_Uend(k-1) + 1 ) )then
                      kkstep=kkstep+1
                   else
                      exit
                   endif
                enddo


                do jj=1,(ILU_Ustart(kk0)-ILU_Uend(kk0)+1)
                   do k = kk0,kk0-kkstep,-1
                      j = ILU_Ustart(k)-jj+1
                      Col = crsMat%crsCol(j)-istart+1 
                      !                z(kk-k+1)=z(kk-k+1)+s(Col)*crsA(j)
#ifdef Manual_inline
                      ddmul_calc1(s_hi(Col),s_lo(Col),crsMat%crsA_hi(j),0.0d0,c_hi,c_lo)
                      ddadd_calc1(z_hi(kk-k+1),z_lo(kk-k+1),c_hi,c_lo,z_hi(kk-k+1),z_lo(kk-k+1))
#else
                      call dd_mul( &
                           s_hi(Col)        ,s_lo(Col),&
                           crsMat%crsA_hi(j),0.0d0, &
                           c_hi,c_lo)
                      call dd_add( &
                           z_hi(kk-k+1),z_lo(kk-k+1), &
                           c_hi,c_lo,&
                           z_hi(kk-k+1),z_lo(kk-k+1))
#endif
                   enddo
                enddo

                kk0 = kk0-kkstep-1
                if(kk0<kkmin)then
                   exit
                endif
             enddo

             do k = kk,max(kk-independ_nvec+1,n_startth),-1
                !             s(k)=s(k)-Dinv(k)*z(kk-k+1)
#ifdef Manual_inline
                ddmul_calc1(Dinv_hi(k),0.0d0,z_hi(kk-k+1),z_lo(kk-k+1),c_hi,c_lo)
                ddsub_calc1(s_hi(k),s_lo(k),c_hi,c_lo,s_hi(k),s_lo(k))
#else
                call dd_mul( &
                     Dinv_hi(k),0.0d0, &
                     z_hi(kk-k+1),z_lo(kk-k+1), &
                     c_hi,c_lo)       
                call dd_sub(&
                     s_hi(k),s_lo(k),&
                     c_hi,c_lo,&
                     s_hi(k),s_lo(k))
#endif
             enddo
          enddo
       endif


    case(KSP_DOUBLE_DOUBLE_PRECISION)
       if(independ_nvec.eq.1)then
          ! -----
          !   z = D^-1(D+U)s
          !   (D+L)z=r
          do k=n_startth,n_endth
             zz_hi = r_hi(k)
             zz_lo = r_lo(k)
!             write(*,*)"in  : CRS -- ",k,r_hi(k),r_lo(k)
             do j=ILU_Lstart(k),ILU_Lend(k)
                Col=crsMat%crsCol(j)-istart+1 
                ! r - Lz の計算
                ! zz=zz-s(Col)*crsA(j)
#ifdef Manual_inline
                ddmul_calc1(s_hi(Col),s_lo(Col),crsMat%crsA_hi(j),crsMat%crsA_lo(j),c_hi,c_lo)
                ddsub_calc1(zz_hi,zz_lo,c_hi,c_lo,zz_hi,zz_lo)
#else
                call dd_mul( &
                     s_hi(Col),s_lo(Col),&
                     crsMat%crsA_hi(j),crsMat%crsA_lo(j), &
                     c_hi,c_lo)
                call dd_sub( &
                     zz_hi,zz_lo, &
                     c_hi,c_lo, &
                     zz_hi,zz_lo)
#endif
             enddo
             ! z = D^-1 (r - Lz ) 
             ! s(k)=Dinv(k)*zz
#ifdef Manual_inline
             ddmul_calc1(zz_hi,zz_lo,Dinv_hi(k),Dinv_lo(k),s_hi(k),s_lo(k))
#else
             call dd_mul( &
                  zz_hi   ,zz_lo,&
                  Dinv_hi(k),Dinv_lo(k), &
                  s_hi(k)   ,s_lo(k) )             
#endif
!             write(*,*)"sav : CRS -- ",k,s_hi(k),s_lo(k)
          enddo

          ! -----
          ! D^-1(D+U)s=z
          !   (I+D^-1U)s=z 
          !   s=z - D^-1 Us       
          do k=n_endth,n_startth,-1
             zz_hi=0.0d0
             zz_lo=0.0d0
             do j=ILU_Ustart(k),ILU_Uend(k),-1 
                Col=crsMat%crsCol(j)-istart+1 
                ! Us の計算
                ! zz=zz+s(Col)*crsA(j)
#ifdef Manual_inline
                ddmul_calc1(s_hi(Col),s_lo(Col),crsMat%crsA_hi(j),crsMat%crsA_lo(j),c_hi,c_lo)
                ddadd_calc1(zz_hi,zz_lo,c_hi,c_lo,zz_hi,zz_lo)
#else
                call dd_mul( &
                     s_hi(Col)        ,s_lo(Col),&
                     crsMat%crsA_hi(j),crsMat%crsA_lo(j), &
                     c_hi,c_lo)
                call dd_add( &
                     zz_hi,zz_lo, &
                     c_hi,c_lo,&
                     zz_hi,zz_lo)                
#endif
             enddo
             ! s(k)=s(k)-Dinv(k)*zz
#ifdef Manual_inline
             ddmul_calc1(Dinv_hi(k),Dinv_lo(k),zz_hi,zz_lo,c_hi,c_lo)
             ddsub_calc1(s_hi(k),s_lo(k),c_hi,c_lo,s_hi(k),s_lo(k))
#else
             call dd_mul( &
                  Dinv_hi(k),Dinv_lo(k), &
                  zz_hi,zz_lo, &
                  c_hi,c_lo)       
             call dd_sub(&
                  s_hi(k),s_lo(k),&
                  c_hi,c_lo,&
                  s_hi(k),s_lo(k))
#endif
!             write(*,*)"out : CRS -- ",k,s_hi(k),s_lo(k)
          enddo

       else

          ! -----
          !   z = D^-1(D+U)s
          !   (D+L)z=r
          do kk=n_startth,n_endth,independ_nvec
             do k=kk,min(kk+independ_nvec-1,n_endth)
                z_hi(k-kk+1)=r_hi(k)
                z_lo(k-kk+1)=r_lo(k)
             enddo
#if 0
             do k=kk,min(kk+independ_nvec-1,n_endth)
                do j=ILU_Lstart(k),ILU_Lend(k)
                   Col=crsCol(j)-istart+1 
                   ! r - Lz の計算
                   z(k-kk+1)=z(k-kk+1)-s(Col)*crsA(j)
                enddo
             enddo
#else
             kk0 = kk
             kkmin = min(kk+independ_nvec-1,n_endth)
             do
                kkstep = 0
                do k = kk0,kkmin-1 
                   if ( ( ILU_Lend(k)-ILU_Lstart(k) + 1 ).eq.( ILU_Lend(k+1)-ILU_Lstart(k+1) + 1 ) )then
                      kkstep=kkstep+1
                   else
                      exit
                   endif
                enddo

                do jj=1,(ILU_Lend(kk0)-ILU_Lstart(kk0)+1)
                   do k = kk0,kk0+kkstep
                      j = ILU_Lstart(k)+jj-1                
                      Col = crsMat%crsCol(j)-istart+1 
                      ! r - Lz の計算
                      !  z(k-kk+1)=z(k-kk+1)-s(Col)*crsA(j)
#ifdef Manual_inline
                      ddmul_calc1(s_hi(Col),s_lo(Col),crsMat%crsA_hi(j),crsMat%crsA_lo(j),c_hi,c_lo)
                      ddsub_calc1(z_hi(k-kk+1),z_lo(k-kk+1),c_hi,c_lo,z_hi(k-kk+1),z_lo(k-kk+1))
#else
                      call dd_mul( &
                           s_hi(Col),s_lo(Col),&
                           crsMat%crsA_hi(j),crsMat%crsA_lo(j), &
                           c_hi,c_lo)
                      call dd_sub( &
                           z_hi(k-kk+1),z_lo(k-kk+1), &
                           c_hi,c_lo, &
                           z_hi(k-kk+1),z_lo(k-kk+1))
#endif
                   enddo
                enddo
                kk0 = kk0+kkstep+1
                if(kk0>kkmin)then
                   exit
                endif
             enddo
#endif

             do k=kk,min(kk+independ_nvec-1,n_endth)
                ! s(k)=Dinv(k)*z(k-kk+1) 
#ifdef Manual_inline
                ddmul_calc1(z_hi(k-kk+1),z_lo(k-kk+1),Dinv_hi(k),Dinv_lo(k),s_hi(k),s_lo(k))
#else
                call dd_mul( &
                     z_hi(k-kk+1)   ,z_lo(k-kk+1),&
                     Dinv_hi(k),Dinv_lo(k), &
                     s_hi(k)   ,s_lo(k) )
#endif
             enddo
          enddo
          ! -----
          ! D^-1(D+U)s=z
          !   (I+D^-1U)s=z 
          !   s=z - D^-1 Us
          do kk=n_endth,n_startth,-independ_nvec
             do k = kk,max(kk-independ_nvec+1,n_startth),-1
                z_hi(kk-k+1)=0.0d0
                z_lo(kk-k+1)=0.0d0
             enddo
#if 0
             do k = kk,max(kk-independ_nvec+1,n_startth),-1
                do j=ILU_Ustart(k),ILU_Uend(k),-1 
                   Col=crsCol(j)-istart+1 
                   ! Us の計算
                   z(kk-k+1)=z(kk-k+1)+s(Col)*crsA(j)
                enddo
             enddo
#else
             kk0 = kk
             kkmin = max(kk-independ_nvec+1,n_startth)
             do 
                kkstep = 0
                do k = kk0,kkmin+1,-1 
                   if ( ( ILU_Ustart(k)-ILU_Uend(k) + 1 ).eq.( ILU_Ustart(k-1) - ILU_Uend(k-1) + 1 ) )then
                      kkstep=kkstep+1
                   else
                      exit
                   endif
                enddo

                do jj=1,(ILU_Ustart(kk0)-ILU_Uend(kk0)+1)
                   do k = kk0,kk0-kkstep,-1
                      j = ILU_Ustart(k)-jj+1
                      Col = crsMat%crsCol(j)-istart+1 
                      !                z(kk-k+1)=z(kk-k+1)+s(Col)*crsA(j)
#ifdef Manual_inline
                      ddmul_calc1(s_hi(Col),s_lo(Col),crsMat%crsA_hi(j),crsMat%crsA_lo(j),c_hi,c_lo)
                      ddadd_calc1(z_hi(kk-k+1),z_lo(kk-k+1),c_hi,c_lo,z_hi(kk-k+1),z_lo(kk-k+1))
#else
                      call dd_mul( &
                           s_hi(Col)        ,s_lo(Col),&
                           crsMat%crsA_hi(j),crsMat%crsA_lo(j), &
                           c_hi,c_lo)
                      call dd_add( &
                           z_hi(kk-k+1),z_lo(kk-k+1), &
                           c_hi,c_lo,&
                           z_hi(kk-k+1),z_lo(kk-k+1))
#endif
                   enddo
                enddo

                kk0 = kk0-kkstep-1
                if(kk0<kkmin)then
                   exit
                endif
             enddo

#endif
             do k = kk,max(kk-independ_nvec+1,n_startth),-1
                !             s(k)=s(k)-Dinv(k)*z(kk-k+1)
#ifdef Manual_inline
                ddmul_calc1(Dinv_hi(k),Dinv_lo(k),z_hi(kk-k+1),z_lo(kk-k+1),c_hi,c_lo)
                ddsub_calc1(s_hi(k),s_lo(k),c_hi,c_lo,s_hi(k),s_lo(k))
#else
                call dd_mul( &
                     Dinv_hi(k),Dinv_lo(k), &
                     z_hi(kk-k+1),z_lo(kk-k+1), &
                     c_hi,c_lo)       
                call dd_sub(&
                     s_hi(k),s_lo(k),&
                     c_hi,c_lo,&
                     s_hi(k),s_lo(k))
#endif
             enddo
          enddo
       endif

    end select

  end subroutine ddsolve_ILU_crs

  subroutine ddsolve_ILU_dia( &
       n,nnd, &! diaMat, &
       diaA_hi,diaA_lo,&
!       offsetA, &
       offset, &
       istart,iend, &
       istart_th,iend_th, &
       n_startth,n_endth, &
       Dinv_hi,Dinv_lo,&
       r_hi,r_lo,&
       s_hi,s_lo, &
       ILU_Lstart, &
       ILU_Lend,   &
       ILU_Ustart, &
       ILU_Uend,   &
       independ_nvec, &
       PRECISION_PRECONDITION)
    use matrix
    use PARCEL_integer
    implicit none

    integer(kind=intType) n
    integer(kind=intType) nnd
    !    type(DIA_Mat) :: diaMat
    real*8 diaA_hi(nnd*n)
    real*8 diaA_lo(nnd*n)
    integer(kind=intType) offset(nnd)

    integer(kind=intType) istart,iend
    integer(kind=intType) istart_th,iend_th
    integer(kind=intType) n_startth,n_endth

    real*8 r_hi(n_startth:n_endth)
    real*8 r_lo(n_startth:n_endth)
    real*8 s_hi(n_startth:n_endth)
    real*8 s_lo(n_startth:n_endth)
    real*8 Dinv_hi(n_startth:n_endth)
    real*8 Dinv_lo(n_startth:n_endth)

    integer(kind=intType) ILU_Lstart(n_startth:n_endth)
    integer(kind=intType) ILU_Lend(n_startth:n_endth)
    integer(kind=intType) ILU_Ustart(n_startth:n_endth)
    integer(kind=intType) ILU_Uend(n_startth:n_endth)
    integer PRECISION_PRECONDITION 

    integer(kind=intType) j,k,kk
    integer(kind=intType) kk0,kkmin,kkstep
    integer(kind=intType) Row,Col
!    integer(kind=intType) offset

    real*8 c_hi,c_lo

    real*8 z_hi(independ_nvec),z_lo(independ_nvec)
    real*8 zz_hi,zz_lo
    integer independ_nvec

    real*8 ddmul_var1
    real*8 ddsub_var1
    real*8 ddadd_var1

!    write(*,*)"PRECISION_PRECONDITION = ",PRECISION_PRECONDITION
!    write(*,*)"KSP_DOUBLE_DOUBLE_PRECISION = ",KSP_DOUBLE_DOUBLE_PRECISION

    select case (PRECISION_PRECONDITION)

    case(KSP_DOUBLE_PRECISION)
!       write(*,*)"ddsolve_ILU_dia double"
       if(independ_nvec.eq.1)then
          ! -----
          ! (D+L)D^-1(D+U)s=r

          ! -----
          !   z = D^-1(D+U)s
          !   (D+L)z=r

          do k=n_startth,n_endth
!             zz=r(k)
             zz_hi=r_hi(k)
             zz_lo=r_lo(k)
             do j=ILU_Lstart(k),ILU_Lend(k)
                ! r - Lz の計算
                ! zz=zz-s(k+offset(j))*diaA(k+(j-1)*n)

#ifdef Manual_inline
                ddmul_calc1(s_hi(k+offset(j)),s_lo(k+offset(j)),diaA_hi(k+(j-1)*n),0.0d0,c_hi,c_lo)
                ddsub_calc1(zz_hi,zz_lo,c_hi,c_lo,zz_hi,zz_lo)
#else
                call dd_mul( &
                     s_hi(k+offset(j)),s_lo(k+offset(j)) , &
                     diaA_hi(k+(j-1)*n) , 0.0d0, &
                     c_hi,c_lo)
                call dd_sub( &
                     zz_hi,zz_lo,&
                     c_hi,c_lo,&
                     zz_hi,zz_lo )
#endif
             enddo
             ! z = D^-1 (r - Lz) 
             ! s(k)=Dinv(k)*zz
#ifdef Manual_inline
             ddmul_calc1(Dinv_hi(k),0.0d0,zz_hi,zz_lo,s_hi(k),s_lo(k))
#else
             call dd_mul( &
                  Dinv_hi(k),0.0d0 , &
                  zz_hi,zz_lo, &
                  s_hi(k),s_lo(k) )
#endif
          enddo

          ! -----
          ! D^-1(D+U)s=z
          !   (I+D^-1U)s=z 
          !   s=z - D^-1 Us
          do k=n_endth,n_startth,-1
             zz_hi=0.0d0
             zz_lo=0.0d0
             do j=ILU_Ustart(k),ILU_Uend(k),-1
!                zz=zz+s(k+offset(j))*diaA(k+(j-1)*n)
#ifdef Manual_inline
                ddmul_calc1(s_hi(k+offset(j)),s_lo(k+offset(j)),diaA_hi(k+(j-1)*n),0.0d0,c_hi,c_lo)
                ddadd_calc1(zz_hi,zz_lo,c_hi,c_lo,zz_hi,zz_lo)
#else
                call dd_mul( &
                     s_hi(k+offset(j)),s_lo(k+offset(j)) , &
                     diaA_hi(k+(j-1)*n) ,0.0d0 , &
                     c_hi,c_lo)
                call dd_add( &
                     zz_hi,zz_lo,&
                     c_hi,c_lo,&
                     zz_hi,zz_lo )
#endif
             enddo
             !             s(k)=s(k)-Dinv(k)*zz
#ifdef Manual_inline
             ddmul_calc1(zz_hi,zz_lo,Dinv_hi(k),0.0d0,c_hi,c_lo)
             ddsub_calc1(s_hi(k),s_lo(k),c_hi,c_lo,s_hi(k),s_lo(k))
#else
             call dd_mul( &
                  zz_hi,zz_lo, &
                  Dinv_hi(k),0.0d0, &
                  c_hi,c_lo)
             call dd_sub( &
                  s_hi(k), s_lo(k), &
                  c_hi,c_lo, &
                  s_hi(k), s_lo(k))
#endif
!             write(*,*)"out : DIA -- ",k,s_hi(k)+s_lo(k)
          enddo

       else
          ! -----
          ! (D+L)D^-1(D+U)s=r

          ! -----
          !   z = D^-1(D+U)s
          !   (D+L)z=r
          do kk=n_startth,n_endth,independ_nvec
             do k=kk,min(kk+independ_nvec-1,n_endth)
!                z(k-kk+1)=r(k)
                z_hi(k-kk+1)=r_hi(k)
                z_lo(k-kk+1)=r_lo(k)
             enddo

             kk0 = kk
             kkmin = min(kk+independ_nvec-1,n_endth)
             do
                kkstep = 0
                do k = kk0,kkmin-1 
                   if ( ( ILU_Lstart(k).eq.ILU_Lstart(k+1) ).and.( ILU_Lend(k).eq.ILU_Lend(k+1) ) )then
                      kkstep=kkstep+1
                   else
                      exit
                   endif
                enddo

                do j=ILU_Lstart(kk0),ILU_Lend(kk0)
                   do k = kk0,kk0+kkstep
                      ! z(k-kk+1)=z(k-kk+1)-s(k+offset(j))*diaA(k+(j-1)*n)
#ifdef Manual_inline
                      ddmul_calc1(s_hi(k+offset(j)),s_lo(k+offset(j)),diaA_hi(k+(j-1)*n),0.0d0,c_hi,c_lo)
                      ddsub_calc1(z_hi(k-kk+1),z_lo(k-kk+1),c_hi,c_lo,z_hi(k-kk+1),z_lo(k-kk+1) )
#else
                      call dd_mul( &
                           s_hi(k+offset(j)),s_lo(k+offset(j)) , &
                           diaA_hi(k+(j-1)*n) , 0.0d0 , &
                           c_hi,c_lo)
                      call dd_sub( &
                           z_hi(k-kk+1),z_lo(k-kk+1),&
                           c_hi,c_lo,&
                           z_hi(k-kk+1),z_lo(k-kk+1) )
#endif
                   enddo
                enddo

                kk0 = kk0+kkstep+1
                if(kk0>kkmin)then
                   exit
                endif
             enddo

             ! z = D^-1 (r - Lz) 
             do k=kk,min(kk+independ_nvec-1,n_endth)
                ! s(k)=Dinv(k)*z(k-kk+1)
#ifdef Manual_inline
                ddmul_calc1(Dinv_hi(k),0.0d0,z_hi(k-kk+1),z_lo(k-kk+1),s_hi(k),s_lo(k))
#else
                call dd_mul( &
                     Dinv_hi(k), 0.0d0 , &
                     z_hi(k-kk+1),z_lo(k-kk+1), &
                     s_hi(k),s_lo(k) )
#endif
             enddo
          enddo

          ! -----
          ! D^-1(D+U)s=z
          !   (I+D^-1U)s=z 
          !   s=z - D^-1 Us
          do kk=n_endth,n_startth,-independ_nvec
             do k = kk,max(kk-independ_nvec+1,n_startth),-1
                z_hi(kk-k+1)=0.0d0
                z_lo(kk-k+1)=0.0d0
             enddo

             kk0 = kk
             kkmin = max(kk-independ_nvec+1,n_startth)
             do 
                kkstep = 0
                do k = kk0,kkmin+1,-1 
                   if ( ( ILU_Ustart(k).eq.ILU_Ustart(k-1) ).and.( ILU_Uend(k).eq.ILU_Uend(k-1) ) )then
                      !             if ( ( ILU_Ustart(k) - ILU_Uend(k) + 1 ).eq.( ILU_Ustart(k-1) - ILU_Uend(k-1) +1 ) )then
                      kkstep=kkstep+1
                   else
                      exit
                   endif
                enddo

                do j=ILU_Ustart(kk0),ILU_Uend(kk0),-1
                   do k = kk0,kk0-kkstep,-1
                      ! z(kk-k+1)=z(kk-k+1)+s(k+offset(j))*diaA(k+(j-1)*n)
#ifdef Manual_inline
                      ddmul_calc1(s_hi(k+offset(j)),s_lo(k+offset(j)),diaA_hi(k+(j-1)*n),0.0d0,c_hi,c_lo)
                      ddadd_calc1(z_hi(kk-k+1),z_lo(kk-k+1),c_hi,c_lo,z_hi(kk-k+1),z_lo(kk-k+1))
#else
                call dd_mul( &
                     s_hi(k+offset(j)),s_lo(k+offset(j)) , &
                     diaA_hi(k+(j-1)*n),0.0d0 , &
                     c_hi,c_lo)
                call dd_add( &
                     z_hi(kk-k+1),z_lo(kk-k+1),&
                     c_hi,c_lo,&
                     z_hi(kk-k+1),z_lo(kk-k+1) )
#endif
                   enddo
                enddo

                kk0 = kk0-kkstep-1
                if(kk0<kkmin)then
                   exit
                endif
             enddo

             do k = kk,max(kk-independ_nvec+1,n_startth),-1
                ! s(k)=s(k)-Dinv(k)*z(kk-k+1)
#ifdef Manual_inline
                ddmul_calc1(z_hi(kk-k+1),z_lo(kk-k+1),Dinv_hi(k),0.0d0,c_hi,c_lo)
                ddsub_calc1(s_hi(k),s_lo(k),c_hi,c_lo,s_hi(k),s_lo(k))
#else
             call dd_mul( &
                  z_hi(kk-k+1),z_lo(kk-k+1), &
                  Dinv_hi(k),0.0d0, &
                  c_hi,c_lo)
             call dd_sub( &
                  s_hi(k), s_lo(k), &
                  c_hi,c_lo, &
                  s_hi(k), s_lo(k))
#endif
             enddo
          enddo
       endif

!       write(*,*)" PRECON DOUBLE DOUBLE"

    case(KSP_DOUBLE_DOUBLE_PRECISION)

       if(independ_nvec.eq.1)then
          ! -----
          ! (D+L)D^-1(D+U)s=r

          ! -----
          !   z = D^-1(D+U)s
          !   (D+L)z=r

          do k=n_startth,n_endth
!             zz=r(k)
             zz_hi=r_hi(k)
             zz_lo=r_lo(k)
!             write(*,*)k,"in  : DIA -- ",k,r_hi(k)+r_lo(k)
!             write(*,*)"loop: DIA -- ",k,ILU_Lstart(k),ILU_Lend(k)
             do j=ILU_Lstart(k),ILU_Lend(k)
                ! r - Lz の計算
                ! zz=zz-s(k+offset(j))*diaA(k+(j-1)*n)
#ifdef Manual_inline
                ddmul_calc1(s_hi(k+offset(j)),s_lo(k+offset(j)),diaA_hi(k+(j-1)*n),diaA_lo(k+(j-1)*n),c_hi,c_lo)
                ddsub_calc1(zz_hi,zz_lo,c_hi,c_lo,zz_hi,zz_lo)
#else
                call dd_mul( &
                     s_hi(k+offset(j)),s_lo(k+offset(j)) , &
                     diaA_hi(k+(j-1)*n) ,diaA_lo(k+(j-1)*n) , &
                     c_hi,c_lo)
                call dd_sub( &
                     zz_hi,zz_lo,&
                     c_hi,c_lo,&
                     zz_hi,zz_lo )
#endif
                 ! write(*,*)"diaA : s -- ",&
                 !      k+offset(j), &
                 !      s_hi(k+offset(j))  + s_lo(k+offset(j))

                ! write(*,*)"diaA : DIA -- ",&
                !      k+(j-1)*n,n, &
                !       diaA_hi(k+(j-1)*n) + diaA_lo(k+(j-1)*n)  
             enddo
             ! z = D^-1 (r - Lz) 
             ! s(k)=Dinv(k)*zz
#ifdef Manual_inline
             ddmul_calc1(Dinv_hi(k),Dinv_lo(k),zz_hi,zz_lo,s_hi(k),s_lo(k))
#else
             call dd_mul( &
                  Dinv_hi(k),Dinv_lo(k) , &
                  zz_hi,zz_lo, &
                  s_hi(k),s_lo(k) )
#endif
!    write(*,*)k,"sav : DIA -- ",Dinv_hi(k)+Dinv_lo(k)  ! ,zz_hi+zz_lo
!             write(*,*)"sav : DIA -- ",k,s_hi(k)+s_lo(k)
          enddo

          ! -----
          ! D^-1(D+U)s=z
          !   (I+D^-1U)s=z 
          !   s=z - D^-1 Us
          do k=n_endth,n_startth,-1
             zz_hi=0.0d0
             zz_lo=0.0d0
             do j=ILU_Ustart(k),ILU_Uend(k),-1
!                zz=zz+s(k+offset(j))*diaA(k+(j-1)*n)
#ifdef Manual_inline
                ddmul_calc1(s_hi(k+offset(j)),s_lo(k+offset(j)),diaA_hi(k+(j-1)*n),diaA_lo(k+(j-1)*n),c_hi,c_lo)
                ddadd_calc1(zz_hi,zz_lo,c_hi,c_lo,zz_hi,zz_lo)
#else
                call dd_mul( &
                     s_hi(k+offset(j)),s_lo(k+offset(j)) , &
                     diaA_hi(k+(j-1)*n) ,diaA_lo(k+(j-1)*n) , &
                     c_hi,c_lo)
                call dd_add( &
                     zz_hi,zz_lo,&
                     c_hi,c_lo,&
                     zz_hi,zz_lo )
#endif
             enddo
             !             s(k)=s(k)-Dinv(k)*zz
#ifdef Manual_inline
             ddmul_calc1(zz_hi,zz_lo,Dinv_hi(k),Dinv_lo(k),c_hi,c_lo)
             ddsub_calc1(s_hi(k),s_lo(k),c_hi,c_lo,s_hi(k),s_lo(k))
#else
             call dd_mul( &
                  zz_hi,zz_lo, &
                  Dinv_hi(k),Dinv_lo(k), &
                  c_hi,c_lo)
             call dd_sub( &
                  s_hi(k), s_lo(k), &
                  c_hi,c_lo, &
                  s_hi(k), s_lo(k))
#endif
!             write(*,*)"out : DIA -- ",k,s_hi(k)+s_lo(k)
          enddo

       else
          ! -----
          ! (D+L)D^-1(D+U)s=r

          ! -----
          !   z = D^-1(D+U)s
          !   (D+L)z=r
          do kk=n_startth,n_endth,independ_nvec
             do k=kk,min(kk+independ_nvec-1,n_endth)
!                z(k-kk+1)=r(k)
                z_hi(k-kk+1)=r_hi(k)
                z_lo(k-kk+1)=r_lo(k)
             enddo

             kk0 = kk
             kkmin = min(kk+independ_nvec-1,n_endth)
             do
                kkstep = 0
                do k = kk0,kkmin-1 
                   if ( ( ILU_Lstart(k).eq.ILU_Lstart(k+1) ).and.( ILU_Lend(k).eq.ILU_Lend(k+1) ) )then
                      kkstep=kkstep+1
                   else
                      exit
                   endif
                enddo

                do j=ILU_Lstart(kk0),ILU_Lend(kk0)
                   do k = kk0,kk0+kkstep
                      ! z(k-kk+1)=z(k-kk+1)-s(k+offset(j))*diaA(k+(j-1)*n)
#ifdef Manual_inline
                      ddmul_calc1(s_hi(k+offset(j)),s_lo(k+offset(j)),diaA_hi(k+(j-1)*n),diaA_lo(k+(j-1)*n),c_hi,c_lo)
                      ddsub_calc1(z_hi(k-kk+1),z_lo(k-kk+1),c_hi,c_lo,z_hi(k-kk+1),z_lo(k-kk+1))
#else
                      call dd_mul( &
                           s_hi(k+offset(j)),s_lo(k+offset(j)) , &
                           diaA_hi(k+(j-1)*n) ,diaA_lo(k+(j-1)*n) , &
                           c_hi,c_lo)
                      call dd_sub( &
                           z_hi(k-kk+1),z_lo(k-kk+1),&
                           c_hi,c_lo,&
                           z_hi(k-kk+1),z_lo(k-kk+1) )
#endif
                   enddo
                enddo

                kk0 = kk0+kkstep+1
                if(kk0>kkmin)then
                   exit
                endif
             enddo

             ! z = D^-1 (r - Lz) 
             do k=kk,min(kk+independ_nvec-1,n_endth)
                ! s(k)=Dinv(k)*z(k-kk+1)
#ifdef Manual_inline
                ddmul_calc1(Dinv_hi(k),Dinv_lo(k),z_hi(k-kk+1),z_lo(k-kk+1),s_hi(k),s_lo(k))
#else
                call dd_mul( &
                     Dinv_hi(k),Dinv_lo(k) , &
                     z_hi(k-kk+1),z_lo(k-kk+1), &
                     s_hi(k),s_lo(k) )
#endif
             enddo
          enddo

          ! -----
          ! D^-1(D+U)s=z
          !   (I+D^-1U)s=z 
          !   s=z - D^-1 Us
          do kk=n_endth,n_startth,-independ_nvec
             do k = kk,max(kk-independ_nvec+1,n_startth),-1
                z_hi(kk-k+1)=0.0d0
                z_lo(kk-k+1)=0.0d0
             enddo

             kk0 = kk
             kkmin = max(kk-independ_nvec+1,n_startth)
             do 
                kkstep = 0
                do k = kk0,kkmin+1,-1 
                   if ( ( ILU_Ustart(k).eq.ILU_Ustart(k-1) ).and.( ILU_Uend(k).eq.ILU_Uend(k-1) ) )then
                      !             if ( ( ILU_Ustart(k) - ILU_Uend(k) + 1 ).eq.( ILU_Ustart(k-1) - ILU_Uend(k-1) +1 ) )then
                      kkstep=kkstep+1
                   else
                      exit
                   endif
                enddo

                do j=ILU_Ustart(kk0),ILU_Uend(kk0),-1
                   do k = kk0,kk0-kkstep,-1
                      ! z(kk-k+1)=z(kk-k+1)+s(k+offset(j))*diaA(k+(j-1)*n)
#ifdef Manual_inline
                      ddmul_calc1(s_hi(k+offset(j)),s_lo(k+offset(j)),diaA_hi(k+(j-1)*n),diaA_lo(k+(j-1)*n),c_hi,c_lo)
                      ddadd_calc1(z_hi(kk-k+1),z_lo(kk-k+1),c_hi,c_lo,z_hi(kk-k+1),z_lo(kk-k+1))
#else
                call dd_mul( &
                     s_hi(k+offset(j)),s_lo(k+offset(j)) , &
                     diaA_hi(k+(j-1)*n),diaA_lo(k+(j-1)*n) , &
                     c_hi,c_lo)
                call dd_add( &
                     z_hi(kk-k+1),z_lo(kk-k+1),&
                     c_hi,c_lo,&
                     z_hi(kk-k+1),z_lo(kk-k+1) )
#endif
                   enddo
                enddo

                kk0 = kk0-kkstep-1
                if(kk0<kkmin)then
                   exit
                endif
             enddo

             do k = kk,max(kk-independ_nvec+1,n_startth),-1
                ! s(k)=s(k)-Dinv(k)*z(kk-k+1)
#ifdef Manual_inline
                ddmul_calc1(z_hi(kk-k+1),z_lo(kk-k+1),Dinv_hi(k),Dinv_lo(k),c_hi,c_lo)
                ddsub_calc1(s_hi(k),s_lo(k),c_hi,c_lo,s_hi(k),s_lo(k))
#else
             call dd_mul( &
                  z_hi(kk-k+1),z_lo(kk-k+1), &
                  Dinv_hi(k),Dinv_lo(k), &
                  c_hi,c_lo)
             call dd_sub( &
                  s_hi(k), s_lo(k), &
                  c_hi,c_lo, &
                  s_hi(k), s_lo(k))
#endif
             enddo
          enddo
       endif

    end select

  end subroutine ddsolve_ILU_dia

  subroutine ddAdditiveSchwarz_crs(method,r_hi,r_lo,s_hi,s_lo)
    use krylov
    use ddvec_kernel
    use vec_kernel
    use PARCEL_integer
    implicit none
    type(KSP) :: method
    real*8 s_hi(method%prm%n)
    real*8 s_lo(method%prm%n)
    real*8 r_hi(method%prm%n)
    real*8 r_lo(method%prm%n)

    real*8 vecL_hi(method%addL)
    real*8 vecL_lo(method%addL)
    real*8 vecU_hi(method%addL)
    real*8 vecU_lo(method%addL)

    integer addL
    integer(kind=intType) addL_intType
    integer iflagAS
    integer(kind=intType) i,i0
    integer(kind=intType) n
    integer myrank

    myrank=method%commVec%myrank
    n = method%prm%n
    iflagAS=method%iflagAS
    addL=method%addL
    addL_intType=addL

    call zero_initialize(addL_intType,vecL_hi)
    call zero_initialize(addL_intType,vecL_lo)
    call zero_initialize(addL_intType,vecU_hi)
    call zero_initialize(addL_intType,vecU_lo)
    ! ---
    if(method%Urank_AS.ne.-1)then
       i0=addL
    else
       i0=0
    endif

    if( (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_BASIC).or.(iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_RESTRICT) )then
       call set_ddAS_vec( &
            method%prm%n, &
            method%commVec%comm, &
            method%addL, &
            method%commVec%myrank, &
            method%Urank_AS, &
            method%Lrank_AS, &
            r_hi,r_lo, &
            method%ivec_AS_hi,method%ivec_AS_lo &
            )
    elseif( (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_INTERPOLATE).or.(iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_NONE) )then
       call zero_initialize(addL_intType,method%ivec_AS_hi(1))
       call zero_initialize(addL_intType,method%ivec_AS_lo(1))
       call xcpy(method%prm%n,r_hi,method%ivec_AS_hi(1+i0))      
       call xcpy(method%prm%n,r_lo,method%ivec_AS_lo(1+i0))      
       call zero_initialize(addL_intType,method%ivec_AS_hi(1+n+i0))
       call zero_initialize(addL_intType,method%ivec_AS_lo(1+n+i0))
    endif
    ! call ddAdditiveSchwarz_ILU_crs( &
    !      method%prm_AS, &
    !      method%precon_crsMatAS, &
    !      method%preconAS, &
    !      method%max_numth, &
    !      iflagAS, &
    !      method%ivec_AS_hi,method%ivec_AS_lo,&
    !      method%independ_nvec, &
    !      method%PRECISION_PRECONDITION)
    call ddAdditiveSchwarz_ILU_crs( &
         method,iflagAS)
         ! iflagAS, &
         ! method%ivec_AS_hi,method%ivec_AS_lo)
    call ddAS_Local( &
         method%prm_AS, &
         method%preconAS, &
         method%max_numth, &
         iflagAS, &
         method%ovec_AS_hi,method%ovec_AS_lo,&
         method%precon_thblock)
    call xcpy(method%prm%n,method%ovec_AS_hi(1+i0),s_hi) 
    call xcpy(method%prm%n,method%ovec_AS_lo(1+i0),s_lo) 
    if((iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_RESTRICT).or.(iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_NONE))then
       return 
    endif

    ! -----
    ! オーバーラップ領域の補正
    ! -----
    call get_ddAS_vec( &
         method%prm%n, &
         method%commVec%comm, &
         method%addL, &
         method%commVec%myrank, &
         method%Urank_AS, &
         method%Lrank_AS, &
         method%ovec_AS_hi, &
         method%ovec_AS_lo, &
         vecL_hi,vecL_lo, &
         vecU_hi,vecU_lo &
         )

    if(method%Lrank_AS.ne.-1)then
       do i=1,addL
!          s(n-addL+i)= ( s(n-addL+i)+vecL(i) ) 
          call dd_add( &
               s_hi(n-addL+i),s_lo(n-addL+i),&
               vecL_hi(i)    ,vecL_lo(i),&
               s_hi(n-addL+i),s_lo(n-addL+i) )
       enddo
    endif
    if(method%Urank_AS.ne.-1)then
       do i=1,addL 
!          s(i)= ( s(i)+vecU(i) ) 
          call dd_add( &
               s_hi(i),s_lo(i),&
               vecU_hi(i)    ,vecU_lo(i),&
               s_hi(i),s_lo(i) )
       enddo
    endif

  end subroutine ddAdditiveSchwarz_crs

  subroutine  ddAdditiveSchwarz_dia(method,r_hi,r_lo,s_hi,s_lo)
    use krylov
    use ddvec_kernel
    use vec_kernel
    use PARCEL_integer
    implicit none
    type(KSP) :: method
    real*8 s_hi(method%prm%n)
    real*8 s_lo(method%prm%n)
    real*8 r_hi(method%prm%n)
    real*8 r_lo(method%prm%n)

    real*8 vecL_hi(method%addL)
    real*8 vecL_lo(method%addL)
    real*8 vecU_hi(method%addL)
    real*8 vecU_lo(method%addL)
    real*8 ddadd_var1
    integer addL
    integer(kind=intType) addL_intType
    integer iflagAS
    integer(kind=intType) i,i0
    integer(kind=intType) n
    integer myrank

    myrank=method%commVec%myrank
    n = method%prm%n
    iflagAS=method%iflagAS
    addL=method%addL
    addL_intType=addL

    call zero_initialize(addL_intType,vecL_hi)
    call zero_initialize(addL_intType,vecL_lo)
    call zero_initialize(addL_intType,vecU_hi)
    call zero_initialize(addL_intType,vecU_lo)
    ! ---
    if(method%Urank_AS.ne.-1)then
       i0=addL
    else
       i0=0
    endif

    if( (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_BASIC).or.(iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_RESTRICT) )then
       call set_ddAS_vec( &
            method%prm%n, &
            method%commVec%comm, &
            method%addL, &
            method%commVec%myrank, &
            method%Urank_AS, &
            method%Lrank_AS, &
            r_hi,r_lo, &
            method%ivec_AS_hi,method%ivec_AS_lo &
            )
    elseif( (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_INTERPOLATE).or.(iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_NONE) )then
       call zero_initialize(addL_intType,method%ivec_AS_hi(1))
       call zero_initialize(addL_intType,method%ivec_AS_lo(1))
       call xcpy(method%prm%n,r_hi,method%ivec_AS_hi(1+i0))
       call xcpy(method%prm%n,r_lo,method%ivec_AS_lo(1+i0))
       call zero_initialize(addL_intType,method%ivec_AS_hi(1+n+i0))
       call zero_initialize(addL_intType,method%ivec_AS_lo(1+n+i0))
    endif

    ! call ddAdditiveSchwarz_ILU_dia( &
    !      method%prm_AS, &
    !      method%precon_diaMatAS, &
    !      method%preconAS, &
    !      method%max_numth, &
    !      iflagAS, &
    !      method%ivec_AS_hi,method%ivec_AS_lo,&
    !      method%PRECISION_PRECONDITION)

    if( method%commVec%myrank .eq. 0 )then
!       call view_vec(method%prm%n,method%ivec_AS_hi)
    endif

    call ddAdditiveSchwarz_ILU_dia( &
         method,iflagAS)

    call ddAS_Local( &
         method%prm_AS, &
         method%preconAS, &
         method%max_numth, &
         iflagAS, &
         method%ovec_AS_hi,method%ovec_AS_lo,&
         method%precon_thblock)
    call xcpy(method%prm%n,method%ovec_AS_hi(1+i0),s_hi)
    call xcpy(method%prm%n,method%ovec_AS_lo(1+i0),s_lo)

    if( method%commVec%myrank .eq. 0 )then
!       call view_vec(method%prm%n,s_hi)
    endif


    if((iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_RESTRICT).or.(iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_NONE))then


       return 
    endif

    ! -----
    ! オーバーラップ領域の補正
    ! -----
    call get_ddAS_vec( &
         method%prm%n, &
         method%commVec%comm, &
         method%addL, &
         method%commVec%myrank, &
         method%Urank_AS, &
         method%Lrank_AS, &
         method%ovec_AS_hi, &
         method%ovec_AS_lo, &
         vecL_hi,vecL_lo,&
         vecU_hi,vecU_lo &
         )

    if(method%Lrank_AS.ne.-1)then
       do i=1,addL
!          s(n-addL+i)= ( s(n-addL+i)+vecL(i) ) 
#ifdef Manual_inline
          ddadd_calc1(s_hi(n-addL+i),s_lo(n-addL+i),vecL_hi(i),vecL_lo(i),s_hi(n-addL+i),s_lo(n-addL+i))
#else
          call dd_add( &
               s_hi(n-addL+i),s_lo(n-addL+i),&
               vecL_hi(i)    ,vecL_lo(i),&
               s_hi(n-addL+i),s_lo(n-addL+i) )
#endif
       enddo
    endif
    if(method%Urank_AS.ne.-1)then
       do i=1,addL 
!          s(i)= ( s(i)+vecU(i) ) 
#ifdef Manual_inline
          ddadd_calc1(s_hi(i),s_lo(i),vecU_hi(i),vecU_lo(i),s_hi(i),s_lo(i))
#else
          call dd_add( &
               s_hi(i),s_lo(i),&
               vecU_hi(i)    ,vecU_lo(i),&
               s_hi(i),s_lo(i) ) 
#endif
       enddo
    endif

  end subroutine ddAdditiveSchwarz_dia

  ! subroutine ddAdditiveSchwarz_ILU_crs( &
  !      method,&
  !      iflagAS,&
  !      r_hi,r_lo)
  subroutine ddAdditiveSchwarz_ILU_crs( &
       method,iflagAS)

    use krylov,only : KSP
    use omp_lib
    use precon
    use matrix
    use ddvec_kernel
    use PARCEL_integer
    implicit none
    type(KSP) :: method
    integer(kind=intType) n
    integer(kind=intType) istart,iend
!    type(prmMat)  :: prm
!    type(CRS_Mat) :: crsMat
!    integer mnumth
!    type(precon_AS) preconAS(mnumth)
!    real*8 r_hi(method%prm%n) ! method%ivec_AS_hi
!    real*8 r_lo(method%prm%n) ! method%ivec_AS_lo

    integer iflagAS
    integer independ_nvec
    integer PRECISION_PRECONDITION

    integer(kind=intType) k
    integer numth,myth,max_numth
    integer(kind=intType) n_th
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th

    ! n=method%prm%n
    ! istart=method%prm%istart
    ! iend=method%prm%iend

    n      = method%prm_AS%n
    istart = method%prm_AS%istart
    iend   = method%prm_AS%iend
!    iflagAS = method%iflagAS
!    precon_thblock=method%precon_thblock 

    !    call zero_initialize(method%prm%n,s)

    ! -----
    ! ブロック数よりもスレッド数が大きい場合はスレッド数を制限する
    call get_maxthreads(max_numth)
    call set_precon_maxthreads(n,method%precon_thblock)

!    call set_precon_maxthreads0(n)
    if( (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_BASIC).or.(iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_RESTRICT) )then
       !$OMP PARALLEL &
       !$OMP PRIVATE( &
       !$OMP numth,myth,n_startth,n_endth,istart_th,iend_th &
       !$OMP ) 
       numth = omp_get_num_threads()
       myth  = omp_get_thread_num()
       istart_th = method%preconAS(myth+1)%istart_th
       iend_th   = method%preconAS(myth+1)%iend_th
       n_startth = method%preconAS(myth+1)%n_startth
       n_endth   = method%preconAS(myth+1)%n_endth
       ! interpolateではオーバーラップ領域の残差をゼロにする必要がある
       call ddsolve_ILU_crs( &
            n,  & !method%crsMat, &
            method%precon_crsMatAS, &
            istart,iend, &
            istart_th,iend_th, &
            n_startth,n_endth, &
            method%preconAS(myth+1)%Dinv_hi, &
            method%preconAS(myth+1)%Dinv_lo, &
            ! r_hi(n_startth), &
            ! r_lo(n_startth), & 
            method%ivec_AS_hi(n_startth), &
            method%ivec_AS_lo(n_startth), & 
            method%preconAS(myth+1)%svec_hi, &
            method%preconAS(myth+1)%svec_lo, &
            method%ILUregionMat(myth+1)%ILU_Lstart, &
            method%ILUregionMat(myth+1)%ILU_Lend,   &
            method%ILUregionMat(myth+1)%ILU_Ustart, &
            method%ILUregionMat(myth+1)%ILU_Uend,   &
            method%independ_nvec,&
            method%PRECISION_PRECONDITION)
       !$OMP END PARALLEL
    elseif( (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_INTERPOLATE).or.(iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_NONE) )then
       !$OMP PARALLEL &
       !$OMP PRIVATE( &
       !$OMP numth,myth,n_startth,n_endth,istart_th,iend_th,k &
       !$OMP ) 
       numth = omp_get_num_threads()
       myth  = omp_get_thread_num()

       n_startth = method%preconAS(myth+1)%n0_startth
       n_endth   = method%preconAS(myth+1)%n0_endth
       do k=n_startth,n_endth
          ! method%preconAS(myth+1)%rvec_hi(k)=r_hi(k)
          ! method%preconAS(myth+1)%rvec_lo(k)=r_lo(k)
          method%preconAS(myth+1)%rvec_hi(k)=method%ivec_AS_hi(k)
          method%preconAS(myth+1)%rvec_lo(k)=method%ivec_AS_lo(k)
       enddo

       istart_th = method%preconAS(myth+1)%istart_th
       iend_th   = method%preconAS(myth+1)%iend_th
       n_startth = method%preconAS(myth+1)%n_startth
       n_endth   = method%preconAS(myth+1)%n_endth
       ! interpolateではオーバーラップ領域の残差をゼロにする必要がある
       call ddsolve_ILU_crs( &
            n, &! method%crsMat, &
            method%precon_crsMatAS, &
            istart,iend, &
            istart_th,iend_th, &
            n_startth,n_endth, &
            method%preconAS(myth+1)%Dinv_hi, &
            method%preconAS(myth+1)%Dinv_lo, &
            method%preconAS(myth+1)%rvec_hi, &
            method%preconAS(myth+1)%rvec_lo, &
            method%preconAS(myth+1)%svec_hi, &
            method%preconAS(myth+1)%svec_lo, &
            method%ILUregionMat(myth+1)%ILU_Lstart, &
            method%ILUregionMat(myth+1)%ILU_Lend,   &
            method%ILUregionMat(myth+1)%ILU_Ustart, &
            method%ILUregionMat(myth+1)%ILU_Uend,   &
            method%independ_nvec, &
            method%PRECISION_PRECONDITION)
       !$OMP END PARALLEL
    endif

!    if(n<=max_numth)then
    call omp_set_num_threads(max_numth)
!    endif

  end subroutine ddAdditiveSchwarz_ILU_crs

  ! subroutine ddAdditiveSchwarz_ILU_dia(&
  !      prm,diaMat,preconAS,mnumth,iflagAS,&
  !      r_hi,r_lo,&
  !      PRECISION_PRECONDITION)
  subroutine ddAdditiveSchwarz_ILU_dia(&
       method,iflagAS)
!       r_hi,r_lo,iflagAS)
!       PRECISION_PRECONDITION)
    use omp_lib
    use krylov,only : KSP
    use vec_kernel
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) nnd
    integer(kind=intType) istart,iend
    type(KSP) :: method

!    type(prmMat)  :: prm
!    type(DIA_Mat) :: diaMat
    integer mnumth
!    type(precon_AS) preconAS(mnumth)
    integer iflagAS
!    real*8 r_hi(method%prm%n)
!    real*8 r_lo(method%prm%n)
    integer PRECISION_PRECONDITION

    integer(kind=intType) k
    integer numth,myth,max_numth
    integer(kind=intType) n_th
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th

    n = method%prm_AS%n
    nnd = method%prm_AS%nnd
    istart = method%prm_AS%istart
    iend = method%prm_AS%iend
    PRECISION_PRECONDITION = method%PRECISION_PRECONDITION

    !    call zero_initialize(prm%n,s)
!    write(*,*)" ddAdditiveSchwarz_ILU_dia "
    ! -----
    ! ブロック数よりもスレッド数が大きい場合はスレッド数を制限する
    call get_maxthreads(max_numth)
    call set_precon_maxthreads(n,method%precon_thblock)
!    call set_precon_maxthreads0(n)
    if( (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_BASIC).or.(iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_RESTRICT) )then
       !$OMP PARALLEL &
       !$OMP PRIVATE( &
       !$OMP numth,myth,n_startth,n_endth,istart_th,iend_th &
       !$OMP ) 
       numth = omp_get_num_threads()
       myth  = omp_get_thread_num()
       istart_th = method%preconAS(myth+1)%istart_th
       iend_th   = method%preconAS(myth+1)%iend_th
       n_startth = method%preconAS(myth+1)%n_startth
       n_endth   = method%preconAS(myth+1)%n_endth

       ! interpolateではオーバーラップ領域の残差をゼロにする必要がある
       call ddsolve_ILU_dia( &
            n,nnd, &
            method%precon_diaMatAS%diaA_hi,method%precon_diaMatAS%diaA_lo, &
            method%precon_diaMatAS%offset, &
            istart,iend, &
            istart_th,iend_th, &
            n_startth,n_endth, &
            method%preconAS(myth+1)%Dinv_hi,&
            method%preconAS(myth+1)%Dinv_lo, &
!            r_hi(n_startth), &
!            r_lo(n_startth), &
            method%ivec_AS_hi(n_startth), & 
            method%ivec_AS_lo(n_startth), & 
            method%preconAS(myth+1)%svec_hi, &
            method%preconAS(myth+1)%svec_lo, &
            method%ILUregionMat(myth+1)%ILU_Lstart, &
            method%ILUregionMat(myth+1)%ILU_Lend,   &
            method%ILUregionMat(myth+1)%ILU_Ustart, &
            method%ILUregionMat(myth+1)%ILU_Uend,   &
            method%independ_nvec,    &
            PRECISION_PRECONDITION)
       !$OMP END PARALLEL
    elseif( (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_INTERPOLATE).or.(iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_NONE) )then
       !$OMP PARALLEL &
       !$OMP PRIVATE( &
       !$OMP numth,myth,n_startth,n_endth,istart_th,iend_th,k &
       !$OMP ) 

       numth = omp_get_num_threads()
       myth  = omp_get_thread_num()
       n_startth = method%preconAS(myth+1)%n0_startth
       n_endth   = method%preconAS(myth+1)%n0_endth

       do k=n_startth,n_endth
          ! method%preconAS(myth+1)%rvec_hi(k)=r_hi(k)
          ! method%preconAS(myth+1)%rvec_lo(k)=r_lo(k)
          method%preconAS(myth+1)%rvec_hi(k)=method%ivec_AS_hi(k)
          method%preconAS(myth+1)%rvec_lo(k)=method%ivec_AS_lo(k)
!          write(*,*)k,"rvec_hi",method%preconAS(myth+1)%rvec_hi(k)
       enddo

       istart_th = method%preconAS(myth+1)%istart_th
       iend_th   = method%preconAS(myth+1)%iend_th
       n_startth = method%preconAS(myth+1)%n_startth
       n_endth   = method%preconAS(myth+1)%n_endth


       ! write(*,*)"@ 1 @",method%preconAS(myth+1)%istart_th,method%preconAS(myth+1)%iend_th
       ! write(*,*)"@ 2 @",method%preconAS(myth+1)%n0_startth,method%preconAS(myth+1)%n0_endth
       ! write(*,*)"@ 3 @",method%preconAS(myth+1)%n_startth,method%preconAS(myth+1)%n_endth



       ! interpolateではオーバーラップ領域の残差をゼロにする必要がある
       call ddsolve_ILU_dia( &
            n,nnd, &
            method%precon_diaMatAS%diaA_hi,method%precon_diaMatAS%diaA_lo, &
            method%precon_diaMatAS%offset, &
            istart,iend, &
            istart_th,iend_th, &
            n_startth,n_endth, &
            method%preconAS(myth+1)%Dinv_hi, &
            method%preconAS(myth+1)%Dinv_lo, &
            method%preconAS(myth+1)%rvec_hi, &
            method%preconAS(myth+1)%rvec_lo, &
            method%preconAS(myth+1)%svec_hi, &
            method%preconAS(myth+1)%svec_lo, &
            method%ILUregionMat(myth+1)%ILU_Lstart, &
            method%ILUregionMat(myth+1)%ILU_Lend,   &
            method%ILUregionMat(myth+1)%ILU_Ustart, &
            method%ILUregionMat(myth+1)%ILU_Uend,   &
            method%independ_nvec,    &
            PRECISION_PRECONDITION)

       ! do k=n_startth,n_endth
       !    write(*,*)k,"rvec_hi",method%preconAS(myth+1)%svec_hi(k)
       ! enddo

       !$OMP END PARALLEL

    endif

!    if(n<=max_numth)then
    call omp_set_num_threads(max_numth)
!    endif

  end subroutine ddAdditiveSchwarz_ILU_dia
    
  subroutine ddAS_Local(prm,preconAS,mnumth,iflagAS,s_hi,s_lo,precon_thblock)
    use omp_lib
    use precon
    use matrix
    use vec_kernel
    use PARCEL_integer
    implicit none

    type(prmMat)  :: prm
    integer mnumth
    type(precon_AS) preconAS(mnumth)
    integer iflagAS
    real*8 s_hi(prm%n)
    real*8 s_lo(prm%n)
    real*8 ddadd_var1
    integer precon_thblock
    integer n_startth0,n_endth0  

    integer(kind=intType) k 
    integer(kind=intType) n
    integer(kind=intType) istart,iend

    integer numth,myth
    integer(kind=intType) istart_th,iend_th
    integer(kind=intType) n_startth,n_endth
    integer max_numth

    n=prm%n
    istart=prm%istart
    iend=prm%iend

    if( (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_RESTRICT).or.(iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_NONE) )then
       ! 局所解
       call get_maxthreads(max_numth)
       call set_precon_maxthreads(n,precon_thblock)
!       call set_precon_maxthreads0(n)
       !$OMP PARALLEL &
       !$OMP PRIVATE( &
       !$OMP numth,myth,n_startth,n_endth,istart_th,iend_th,k &
       !$OMP )       

       numth = omp_get_num_threads()
       myth  = omp_get_thread_num()
#if 0
       call divide_region( &
            myth,numth,1,n,n_startth,n_endth)
#else
       n_startth = preconAS(myth+1)%n0_startth
       n_endth   = preconAS(myth+1)%n0_endth
#endif
!       write(*,*)"@@n_startth,n_endth",n_startth,n_endth
       do k=n_startth,n_endth
          s_hi(k)=preconAS(myth+1)%svec_hi(k)
          s_lo(k)=preconAS(myth+1)%svec_lo(k)
!          write(*,*)k,s_hi(k)
       enddo
       !$OMP END PARALLEL 
!       if(n<=max_numth)then
       call omp_set_num_threads(max_numth)
!       endif
    else
       ! 加算補正
       call get_maxthreads(max_numth)
!       call set_precon_maxthreads0(n)
       call set_precon_maxthreads(n,precon_thblock)
       call zero_initialize(n,s_hi)
       call zero_initialize(n,s_lo)
#if 1
       !$OMP PARALLEL &
       !$OMP PRIVATE( &
       !$OMP numth,myth,n_startth,n_endth,istart_th,iend_th,k, &
       !$OMP n_startth0,n_endth0,  &
       !$OMP ddadd_var1 &
       !$omp ) 
       numth = omp_get_num_threads()
       myth  = omp_get_thread_num()
       n_startth = preconAS(myth+1)%n_startth
       n_endth   = preconAS(myth+1)%n_endth       

       n_startth0 = preconAS(myth+1)%n0_startth
       n_endth0   = preconAS(myth+1)%n0_endth

       do k=n_startth,n_startth0-1
!          s(k)=s(k)+preconAS(myth+1)%svec(k)
#ifdef Manual_inline
          ddadd_calc1(s_hi(k),s_lo(k),preconAS(myth+1)%svec_hi(k),preconAS(myth+1)%svec_lo(k),s_hi(k),s_lo(k))
#else
          ddadd_calc1(s_hi(k),s_lo(k),preconAS(myth+1)%svec_hi(k),preconAS(myth+1)%svec_lo(k),s_hi(k),s_lo(k))
#endif
       enddo
       !$OMP BARRIER
       do k=n_startth0,n_endth0
!          s(k)=s(k)+preconAS(myth+1)%svec(k)
#ifdef Manual_inline
          ddadd_calc1(s_hi(k),s_lo(k),preconAS(myth+1)%svec_hi(k),preconAS(myth+1)%svec_lo(k),s_hi(k),s_lo(k))
#else
          call dd_add( &
               s_hi(k),s_lo(k),&
               preconAS(myth+1)%svec_hi(k),preconAS(myth+1)%svec_lo(k), &
               s_hi(k),s_lo(k))
#endif
       enddo
       !$OMP BARRIER
       do k=n_endth0+1,n_endth
!          s(k)=s(k)+preconAS(myth+1)%svec(k)
#ifdef Manual_inline
          ddadd_calc1(s_hi(k),s_lo(k),preconAS(myth+1)%svec_hi(k),preconAS(myth+1)%svec_lo(k),s_hi(k),s_lo(k))
#else
          call dd_add( &
               s_hi(k),s_lo(k),&
               preconAS(myth+1)%svec_hi(k),preconAS(myth+1)%svec_lo(k), &
               s_hi(k),s_lo(k))
#endif
       enddo

       !$OMP END PARALLEL 

#else
       !$OMP PARALLEL &
       !$OMP PRIVATE( &
       !$OMP numth,myth,n_startth,n_endth,istart_th,iend_th,k, &
       !$OMP ddadd_var1 &
       !$OMP ) 
       !$OMP REDUCTION(+:s) 
       numth = omp_get_num_threads()
       myth  = omp_get_thread_num()
!        istart_th = preconAS(myth+1)%istart_th
!        iend_th   = preconAS(myth+1)%iend_th
       n_startth = preconAS(myth+1)%n_startth
       n_endth   = preconAS(myth+1)%n_endth       
       do k=n_startth,n_endth
          !          s(k)=s(k)+preconAS(myth+1)%svec(k)
          !$OMP CRITICAL
#ifdef Manual_inline
          ddadd_calc1(s_hi(k),s_lo(k),preconAS(myth+1)%svec_hi(k),preconAS(myth+1)%svec_lo(k),s_hi(k),s_lo(k))
#else
          call dd_add( &
               s_hi(k),s_lo(k),&
               preconAS(myth+1)%svec_hi(k),preconAS(myth+1)%svec_lo(k), &
               s_hi(k),s_lo(k))
#endif
          !$OMP END CRITICAL
       enddo
       !$OMP END PARALLEL 
#endif

!       if(n<=max_numth)then
       call omp_set_num_threads(max_numth)
!       endif
    endif

  end subroutine ddAS_Local
! ---

  subroutine ddinit_crsMat(crsMat,n,nnz)
    use matrix
    use PARCEL_integer
    implicit none
    type(CRS_Mat) :: crsMat 
    integer(kind=intType) n
    integer(kind=intType) nnz
    integer(kind=intType) i

    allocate( crsMat%crsA_hi(nnz) )
    allocate( crsMat%crsA_lo(nnz) )
    allocate( crsMat%crsRow_ptr(n+1) )
    allocate( crsMat%crsCol(nnz) )

    allocate( crsMat%crsRow_ptr_Diag_s(n) )
    allocate( crsMat%crsRow_ptr_Diag_e(n) )

    !$OMP PARALLEL DO PRIVATE(i)
    do i=1,n+1
       crsMat%crsRow_ptr(i)=0
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i)
    do i=1,nnz
       crsMat%crsA_hi(i)=0.0d0
       crsMat%crsA_lo(i)=0.0d0
       crsMat%crsCol(i)=0
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i)
    do i=1,n
       crsMat%crsRow_ptr_Diag_s(i)=0
       crsMat%crsRow_ptr_Diag_e(i)=0
    enddo
    !$OMP END PARALLEL DO

    allocate(crsMat%RowL(1))
    allocate(crsMat%RowD(1))
    allocate(crsMat%RowU(1))

  end subroutine ddinit_crsMat


  subroutine ddinit_diaMat(diaMat,n,nnd)
    use matrix
    use PARCEL_integer
    implicit none
    type(DIA_Mat) :: diaMat 
    integer(kind=intType) n
    integer(kind=intType) nnd
    integer(kind=intType) i
    allocate( diaMat%diaA_hi(n*nnd) )
    allocate( diaMat%diaA_lo(n*nnd) )
    allocate( diaMat%offset(nnd) )
    allocate( diaMat%offsetVec(nnd) )
    allocate( diaMat%Vec_hi(1) )
    allocate( diaMat%Vec_lo(1) )

    allocate(diaMat%RowL(1))
    allocate(diaMat%RowD(1))
    allocate(diaMat%RowU(1))

    !$OMP PARALLEL DO PRIVATE(i)
    do i=1,n*nnd
       diaMat%diaA_hi(i)=0.0d0
       diaMat%diaA_lo(i)=0.0d0
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i)
    do i=1,nnd
       diaMat%offset(i)=0
       diaMat%offsetVec(i)=0
    enddo
    !$OMP END PARALLEL DO

  end subroutine ddinit_diaMat

  subroutine dddivide_ILU_crs( &
       n,istart,&
       crsMat,  &
       crsMatL, &
       crsMatU  &
       )
    ! ILU分解した行列を上三角成分と下三角成分に分ける。
    ! 下三角行列の対角成分は1。
    
    use matrix
    use PARCEL_integer
    implicit none
    type(CRS_Mat) :: crsMat 
    type(CRS_Mat) :: crsMatL 
    type(CRS_Mat) :: crsMatU 

    integer(kind=intType) i,j,k
    integer(kind=intType) n,istart
    integer(kind=intType) nnzL,nnzU
    integer(kind=intType) Row,Col
    real*8  matA_hi,matA_lo

    nnzL=0
    nnzU=0
    do k=1,n
       Row=k-1+istart
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          if(crsMat%crsCol(j)<Row)then
             nnzL=nnzL+1
          endif
          if(Row<crsMat%crsCol(j))then
             nnzU=nnzU+1
          endif
          if(Row.eq.crsMat%crsCol(j))then
             nnzL=nnzL+1
             nnzU=nnzU+1
          endif
       enddo
    enddo
    !    write(*,*)"nnzL,nnzU",nnzL,nnzU
    call ddinit_crsMat(crsMatL,n,nnzL)
    call ddinit_crsMat(crsMatU,n,nnzU)

    nnzL=1
    nnzU=1
    do k=1,n
       crsMatL%crsRow_ptr(k)=nnzL
       crsMatU%crsRow_ptr(k)=nnzU
       Row=k-1+istart
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          Col =crsMat%crsCol(j)
          matA_hi=crsMat%crsA_hi(j)
          matA_lo=crsMat%crsA_lo(j)
          if(crsMat%crsCol(j)<Row)then
             crsMatL%crsA_hi(nnzL)   = matA_hi
             crsMatL%crsA_lo(nnzL)   = matA_lo
             crsMatL%crsCol(nnzL) = Col
             nnzL=nnzL+1
          endif
          if(Row<crsMat%crsCol(j))then
             crsMatU%crsA_hi(nnzU)   = matA_hi
             crsMatU%crsA_lo(nnzU)   = matA_lo
             crsMatU%crsCol(nnzU) = Col
             nnzU=nnzU+1
          endif
          if(Row.eq.crsMat%crsCol(j))then
             crsMatL%crsA_hi(nnzL)   = 1.0d0
             crsMatL%crsA_lo(nnzL)   = 0.0d0
             crsMatL%crsCol(nnzL) = Col
             nnzL=nnzL+1

             crsMatU%crsA_hi(nnzU)   = matA_hi 
             crsMatU%crsA_lo(nnzU)   = matA_lo 
             crsMatU%crsCol(nnzU) = Col
             nnzU=nnzU+1
          endif
       enddo
    enddo
    crsMatL%crsRow_ptr(n+1)=nnzL
    crsMatU%crsRow_ptr(n+1)=nnzU

  end subroutine dddivide_ILU_crs

  subroutine dddivide_ILU_dia( &
       n,istart,nnd,&
       diaMat,  &
       diaMatL, &
       diaMatU  &
       )
    ! ILU分解した行列を上三角成分と下三角成分に分ける。
    ! 下三角行列の対角成分は1。
    
    use matrix
    use PARCEL_integer
    implicit none
    type(DIA_Mat) :: diaMat 
    type(DIA_Mat) :: diaMatL 
    type(DIA_Mat) :: diaMatU 

    integer(kind=intType) n,istart,nnd
    integer(kind=intType) i,j,k
    integer(kind=intType) nndL,nndU
    integer(kind=intType) Row,Col
    integer(kind=intType) offset

    nndL=0
    nndU=0
    do i=1,nnd
       offset=diaMat%offset(i)
       if(offset<0)then
          nndL=nndL+1
       endif
       if(0<offset)then
          nndU=nndU+1
       endif
    enddo
    ! 主対角があるため+1
    nndL=nndL+1
    nndU=nndU+1
    call ddinit_diaMat(diaMatL,n,nndL) 
    call ddinit_diaMat(diaMatU,n,nndU) 

    nndL=1
    nndU=1
    do i=1,nnd
       offset=diaMat%offset(i)
       if(offset<0)then
          do j=1,n
             diaMatL%diaA_hi(j+(nndL-1)*n)=diaMat%diaA_hi(j+(i-1)*n)
             diaMatL%diaA_lo(j+(nndL-1)*n)=diaMat%diaA_lo(j+(i-1)*n)
!             write(*,*)" init ",j,j+(nndL-1)*n,diaMatL%diaA_hi(j+(nndL-1)*n)
          enddo
          diaMatL%offset(nndL)=offset
          nndL=nndL+1
       endif
       if(offset.eq.0)then
          do j=1,n
             diaMatU%diaA_hi(j+(nndU-1)*n)=diaMat%diaA_hi(j+(i-1)*n)
             diaMatU%diaA_lo(j+(nndU-1)*n)=diaMat%diaA_lo(j+(i-1)*n)
             diaMatL%diaA_hi(j+(nndL-1)*n)=1.0d0
             diaMatL%diaA_lo(j+(nndL-1)*n)=0.0d0
!             write(*,*)" init ",j,j+(nndL-1)*n,diaMatL%diaA_hi(j+(nndL-1)*n)
          enddo
          diaMatL%offset(nndL)=offset
          nndL=nndL+1
          diaMatU%offset(nndU)=offset
          nndU=nndU+1
       endif
       if(0<offset)then
          do j=1,n
             diaMatU%diaA_hi(j+(nndU-1)*n)=diaMat%diaA_hi(j+(i-1)*n)
             diaMatU%diaA_lo(j+(nndU-1)*n)=diaMat%diaA_lo(j+(i-1)*n)
          enddo
          diaMatU%offset(nndU)=offset
          nndU=nndU+1
       endif
    enddo

  end subroutine dddivide_ILU_dia


  subroutine ddILU_LU_crs(&
       n,istart,nnz,&
       crsA_hi,crsA_lo,&
       crsRow_ptr,crsCol,&
       crsMatL,crsMatU)
!    use krylov
    use matrix
!    use common_kernel
    use PARCEL_integer
    implicit none

  integer(kind=intType) n,nnz,istart
  real*8 crsA_hi(nnz),crsA_lo(nnz)
  integer(kind=intType) crsRow_ptr(n+1)
  integer(kind=intType) crsCol(nnz)

  type(CRS_Mat) :: crsMatL
  type(CRS_Mat) :: crsMatU
  type(CRS_Mat) :: crsMatLU

  integer(kind=intType) nnzL
  integer(kind=intType) nnzU

  integer(kind=intType) i,j,k
  integer(kind=intType) ij,jk,kj,kk,ik
  integer(kind=intType) Col
  integer(kind=intType) kCol,jCol
  real*8 ddmul_var1
  real*8 ddsub_var1
  real*8 dddiv_var1
  real*8 c_hi,c_lo

  c_hi=0.0d0
  c_lo=0.0d0

  call ddinit_crsMat(crsMatLU,n,nnz) 

  !$OMP PARALLEL DO PRIVATE(i)
  do i=1,nnz
     crsMatLU%crsA_hi(i) = crsA_hi(i)
     crsMatLU%crsA_lo(i) = crsA_lo(i)
     crsMatLU%crsCol(i)  = crsCol(i)
  enddo

  !$OMP PARALLEL DO PRIVATE(i)
  do i=1,n+1
     crsMatLU%crsRow_ptr(i)=crsRow_ptr(i)
  enddo

  do i=2,n
     do kCol=crsMatLU%crsRow_ptr(i),crsMatLU%crsRow_ptr(i+1)-1 
        k=crsMatLU%crsCol(kCol)-istart+1
        if(k>i-1)then
           exit
        endif

        kk=-1
        do j=crsMatLU%crsRow_ptr(k),crsMatLU%crsRow_ptr(k+1)-1
           Col=crsMatLU%crsCol(j) - istart+1
           if( Col.eq.k )then
              kk=j ! k 行の対角位置
              exit
           endif
        enddo
        if(kk.eq.-1)then
           write(*,*)"Error :: ILU -- 対角成分が存在しない"
           stop
        endif

        ! LU(i,k)=LU(i,k)/LU(k,k)
        do j=crsMatLU%crsRow_ptr(i),crsMatLU%crsRow_ptr(i+1)-1
           Col=crsMatLU%crsCol(j)-istart+1
           if( Col.eq.k )then
!              crsMatLU%crsA(j) = crsMatLU%crsA(j) / crsMatLU%crsA(kk)
#ifdef Manual_inline
              dddiv_calc1(crsMatLU%crsA_hi(j),crsMatLU%crsA_lo(j),crsMatLU%crsA_hi(kk),crsMatLU%crsA_lo(kk),crsMatLU%crsA_hi(j),crsMatLU%crsA_lo(j))
#else
              call dd_div( &
                   crsMatLU%crsA_hi(j) ,crsMatLU%crsA_lo(j) , &
                   crsMatLU%crsA_hi(kk),crsMatLU%crsA_lo(kk), &
                   crsMatLU%crsA_hi(j) ,crsMatLU%crsA_lo(j) )
#endif
              exit
           endif
        enddo

        ! LU(i,j)=LU(i,j)-LU(i,k)*LU(k,j)
        do jCol=crsMatLU%crsRow_ptr(i),crsMatLU%crsRow_ptr(i+1)-1
           j=crsMatLU%crsCol(jCol)-istart+1
           if(j<k+1)then
              cycle
           endif

           do ij=crsMatLU%crsRow_ptr(i),crsMatLU%crsRow_ptr(i+1)-1
              Col=crsMatLU%crsCol(ij)-istart+1
              if( Col.eq.j )then
                 ! (i,j)位置が存在する場合
                 ! (i,k)位置の探査
                 do ik=crsMatLU%crsRow_ptr(i),crsMatLU%crsRow_ptr(i+1)-1
                    Col=crsMatLU%crsCol(ik)-istart+1
                    if(Col.eq.k)then
                       ! (k,j)位置の探査
                       do kj=crsMatLU%crsRow_ptr(k),crsMatLU%crsRow_ptr(k+1)-1
                          Col=crsMatLU%crsCol(kj)-istart+1
                          if(Col.eq.j)then
!                              crsMatLU%crsA(ij) = &
!                                   crsMatLU%crsA(ij) - &
!                                   crsMatLU%crsA(ik) * crsMatLU%crsA(kj) 
! ---
#ifdef Manual_inline
                             ddmul_calc1(crsMatLU%crsA_hi(ik),crsMatLU%crsA_lo(ik),crsMatLU%crsA_hi(kj),crsMatLU%crsA_lo(kj),c_hi,c_lo)
                             ddsub_calc1(crsMatLU%crsA_hi(ij),crsMatLU%crsA_lo(ij),c_hi,c_lo,crsMatLU%crsA_hi(ij),crsMatLU%crsA_lo(ij))
#else
                             call dd_mul( &
                                  crsMatLU%crsA_hi(ik),crsMatLU%crsA_lo(ik), &
                                  crsMatLU%crsA_hi(kj),crsMatLU%crsA_lo(kj), &
                                  c_hi,c_lo)
                             call dd_sub( &
                                  crsMatLU%crsA_hi(ij),crsMatLU%crsA_lo(ij), &
                                  c_hi,c_lo, &
                                  crsMatLU%crsA_hi(ij),crsMatLU%crsA_lo(ij) )
#endif
                             exit
                          endif
                       enddo
                       exit
                    endif
                 enddo
              endif

           enddo

        enddo

     enddo
  enddo

  call dddivide_ILU_crs( &
       n,istart,&
       crsMatLU,crsMatL,crsMatU)

  call ddfree_crsMat(crsMatLU)
  
  end subroutine ddILU_LU_crs

  subroutine ddILU_LU_dia(&
       n,istart,nnd,&
       diaA_hi,diaA_lo,&
       offset,&
       diaMatL,diaMatU)
    use krylov
    use PARCEL_integer
    implicit none

    integer(kind=intType) n,istart,nnd
    real*8 diaA_hi(nnd*n)
    real*8 diaA_lo(nnd*n)
    real*8 c_hi,c_lo
    real*8 dddiv_var1
    real*8 ddmul_var1
    real*8 ddsub_var1

    integer(kind=intType) offset(nnd)
    type(DIA_Mat) :: diaMatL
    type(DIA_Mat) :: diaMatU
    type(DIA_Mat) :: diaMatLU

    integer(kind=intType) i,j,k
    integer(kind=intType) ij,jk,kj,kk,ik
    integer(kind=intType) kCol,jCol
    integer(kind=intType) Row,Col
    integer(kind=intType) nnd_ofset0

    call ddinit_diaMat(diaMatLU,n,nnd) 
    do i=1,nnd*n
       diaMatLU%diaA_hi(i)   = diaA_hi(i)
       diaMatLU%diaA_lo(i)   = diaA_lo(i)
    enddo

    nnd_ofset0=-1
    do i=1,nnd
       diaMatLU%offset(i) = offset(i)
       if( offset(i).eq.0)then
          nnd_ofset0=i ! 主対角のオフセット番号
       endif
    enddo

    if(nnd_ofset0==-1) then
       write(*,*)"error ::ILU -- There is no diagonal element."
       stop
    endif
    do i=1,n
       if ( (diaMatLU%diaA_hi(i+n*(nnd_ofset0-1)).eq.0.0d0) .and. &
            (diaMatLU%diaA_lo(i+n*(nnd_ofset0-1)).eq.0.0d0) ) then
         write(*,*)"error::ILU-- It's lack of diagonal element."
         stop
       endif
    enddo

    do i=2,n
       !        write(*,*)"---------------------"
       !        write(*,*)"  i = ",i
       !        write(*,*)"---------------------"
       do k=1,nnd
          kCol = i + diaMatLU%offset(k) ! ローカル列番号
          if(kCol>i-1)then
             cycle
          endif
          if((kCol<1).or.(n<kCol))then
             cycle
          endif

          ! 下三角成分
          ! LU(i,k)=LU(i,k)/LU(k,k)
          !           write(*,*)"i,k",i   ,i    + diaMatLU%offset(k)
          !           write(*,*)"k,k",kCol,kCol + diaMatLU%offset(nnd_ofset0)

          !if( diaMatLU%diaA(kCol+n*(nnd_ofset0-1)).ne.0.0d0 )then
          ! diaMatLU%diaA(i+n*(k-1)) = &
          !    diaMatLU%diaA(i+n*(k-1)) / &
          !    diaMatLU%diaA(kCol+n*(nnd_ofset0-1))

#ifdef Manual_inline
          dddiv_calc1(diaMatLU%diaA_hi(i+n*(k-1)),diaMatLU%diaA_lo(i+n*(k-1)),diaMatLU%diaA_hi(kCol+n*(nnd_ofset0-1)),diaMatLU%diaA_lo(kCol+n*(nnd_ofset0-1)),diaMatLU%diaA_hi(i+n*(k-1)),diaMatLU%diaA_lo(i+n*(k-1)))
#else
          call dd_div( &
               diaMatLU%diaA_hi(i+n*(k-1)), &
               diaMatLU%diaA_lo(i+n*(k-1)), &
               diaMatLU%diaA_hi(kCol+n*(nnd_ofset0-1)), &
               diaMatLU%diaA_lo(kCol+n*(nnd_ofset0-1)), &
               diaMatLU%diaA_hi(i+n*(k-1)), &
               diaMatLU%diaA_lo(i+n*(k-1)) )
#endif
          !endif
          !          write(*,*)"-----"
          ! 上三角成分
          ! LU(i,j)=LU(i,j)-LU(i,k)*LU(k,j)
          do j=1,nnd
             jCol = i + diaMatLU%offset(j) ! ローカル列番号
             if(jCol<kCol+1)then
                cycle
             endif
             if((jCol<1).or.(n<jCol))then
                cycle
             endif
             ! DIA形式ではゼロ要素の番号も含まれるため除外する
             ! 除外しない場合はfill in を少し考慮することに相当
             if( (diaMatLU%diaA_hi(i+n*(j-1)).eq.0.0d0) .and. (diaMatLU%diaA_lo(i+n*(j-1)).eq.0.0d0) )then
                cycle
             endif
             !             write(*,*)"i,j",i,i+diaMatLU%offset(j)
             do kj=nnd_ofset0+1,nnd
                if(jCol.eq.(kCol+diaMatLU%offset(kj)))then
                   !                    write(*,*)"i,k",i,i+diaMatLU%offset(k)
                   !                    write(*,*)"k,j",i+diaMatLU%offset(k),i+diaMatLU%offset(j)            
                   ! diaMatLU%diaA(i+n*(j-1)) = &
                   !      diaMatLU%diaA(i+n*(j-1)) - &
                   !      diaMatLU%diaA(i+n*(k-1)) * diaMatLU%diaA( kCol + n*(kj-1))
#ifdef Manual_inline
                   ddmul_calc1(diaMatLU%diaA_hi(i+n*(k-1)),diaMatLU%diaA_lo(i+n*(k-1)),diaMatLU%diaA_hi(kCol+n*(kj-1)),diaMatLU%diaA_lo(kCol+n*(kj-1)),c_hi,c_lo)
                   ddsub_calc1(diaMatLU%diaA_hi(i+n*(j-1)),diaMatLU%diaA_lo(i+n*(j-1)),c_hi,c_lo,diaMatLU%diaA_hi(i+n*(j-1)),diaMatLU%diaA_lo(i+n*(j-1)))
#else
                   call dd_mul( &
                        diaMatLU%diaA_hi( i+n*(k-1) ), &
                        diaMatLU%diaA_lo( i+n*(k-1) ), &
                        diaMatLU%diaA_hi( kCol + n*(kj-1) ), &
                        diaMatLU%diaA_lo( kCol + n*(kj-1) ), &
                        c_hi,c_lo)
                   call dd_sub( &
                        diaMatLU%diaA_hi(i+n*(j-1)), &
                        diaMatLU%diaA_lo(i+n*(j-1)), &
                        c_hi, &
                        c_lo, &
                        diaMatLU%diaA_hi(i+n*(j-1)),&
                        diaMatLU%diaA_lo(i+n*(j-1)) )
#endif
!                   write(*,*)i,j," LU -- " ,diaMatLU%diaA_hi(i+n*(j-1)),diaMatLU%diaA_lo(i+n*(j-1))
                   exit
                endif
             enddo
             !              diaMatLU%diaA(i+n*(j-1)) = &
             !                   diaMatLU%diaA(i+n*(j-1)) - &
             !                   diaMatLU%diaA(i+n*(k-1)) * diaMatLU%diaA( kCol + n*(j-1))
          enddo
       enddo
    enddo

    call dddivide_ILU_dia( &
         n,istart,nnd,&
         diaMatLU,  &
         diaMatL, &
         diaMatU  &
         )

    call ddfree_diaMat(diaMatLU)



  end subroutine ddILU_LU_dia

  subroutine ddget_nnd_crs_without_sort(&
           diaMat,n,nnz,crsCol,crsRow,istart,iend,nnd)
    use matrix
    implicit none
    type(DIA_Mat) :: diaMat
    integer(kind=intType) n
    integer(kind=intType) nnz
    integer(kind=intType) crsCol(nnz)
    integer(kind=intType) crsRow(n+1)
    integer(kind=intType) istart
    integer(kind=intType) iend

    integer(kind=intType) nnd
  
    integer(kind=intType) i, j, icol, irow, itmp
    integer(kind=intType),allocatable::Flag(:)
    integer(kind=intType) nnd_p

    integer imn,imx
  
    imn = 0
    imx = 0
    do i=1,n
      do j=crsRow(i),crsRow(i+1)-1
        irow = i+istart-1
        icol = crsCol(j)
        if(irow>=istart .and. irow<=iend) then
          itmp=icol-irow
          if (imn>itmp) imn = itmp
          if (imx<itmp) imx = itmp
        endif
      enddo
    enddo
    allocate(Flag(imn:imx))
    Flag = 0

    do i=1,n
      do j=crsRow(i),crsRow(i+1)-1
        irow = i+istart-1
        icol = crsCol(j)
        if(irow>=istart .and. irow<=iend) then
          itmp = icol-irow
          Flag(itmp) = 1
        endif
      enddo
    enddo
    nnd = 0

    do i=imn,imx
      if (Flag(i)==1) then
        nnd=nnd+1
      endif
    enddo

!    call init_diaMat(diaMat,n,nnd)
    call ddinit_diaMat(diaMat,n,nnd)
    nnd_p=0
    do i=imn,imx
      if(flag(i)==1) then
        nnd_p=nnd_p+1
        diaMat%offset(nnd_p)=i
      endif
    enddo
    deallocate(Flag)

  end subroutine ddget_nnd_crs_without_sort

  subroutine ddcrs2dia_without_sort(&
       crsMat,diaMat,&
       nnd,n,nnz, &
       istart)
    use matrix

    implicit none
    type(CRS_Mat) :: crsMat
    type(DIA_Mat) :: diaMat
    integer(kind=intType) n,itmp
    integer(kind=intType) nnz
    integer(kind=intType) nnd

    real*8 crsA(nnz)
    integer(kind=intType) crsCol(nnz)
    integer(kind=intType) crsRow_ptr(n+1)

    real*8 diaA(nnd*n)
    integer(kind=intType) offset(nnd)

    integer(kind=intType) istart

    integer(kind=intType) i,j,k
    integer(kind=intType) icol,irow
    integer(kind=intType) off  

    do i=1,n
       do j=crsMat%crsRow_ptr(i),crsMat%crsRow_ptr(i+1)-1
          irow = i+istart-1
          icol = crsMat%crsCol(j)
          itmp = icol-irow
          do k = 1,nnd
             off =k 
             if (diaMat%offset(k)==itmp) exit
          enddo
          diaMat%diaA_hi(n*(off-1)+i) = crsMat%crsA_hi(j)
          diaMat%diaA_lo(n*(off-1)+i) = crsMat%crsA_lo(j)
!          write(*,*)n*(off-1)+i,i,j,crsMat%crsA_hi(j) + crsMat%crsA_lo(j)
       enddo
    enddo

  end subroutine ddcrs2dia_without_sort


  subroutine ddget_ILU_LU_crs(method,crsMat)
    use krylov
    use omp_lib
    use matrix_kernel,only : get_nnz_crs

    implicit none
    type(KSP) :: method
    type(CRS_Mat) :: crsMat 

    integer max_numth
    integer numth,myth
    integer(kind=intType) n,istart
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th
    integer(kind=intType) n_tmp
    integer(kind=intType) n_th,nnz_th
    integer(kind=intType) i,j,Col,i_th,ii,iimin
    integer independ_nvec
    integer precon_thblock

#if 0
    n      = method%prm%n
    istart = method%prm%istart

    ! ----
    ! 計算領域の分割
    call get_maxthreads(max_numth)
!    call set_precon_maxthreads0(n)
    call set_precon_maxthreads(n,method%precon_thblock)

    !$OMP PARALLEL 
    !$OMP MASTER
    numth = omp_get_num_threads()
    !$OMP END MASTER
    !$OMP END PARALLEL 
    allocate(method%preconILU(numth))

    !$OMP PARALLEL &
    !$OMP PRIVATE( &
    !$OMP numth,myth,&
    !$OMP n_startth,n_endth,&
    !$OMP istart_th,iend_th,&
    !$OMP n_tmp, &
    !$OMP n_th,nnz_th, &
    !$OMP i,j,i_th,Col &
    !$OMP )
    numth = omp_get_num_threads()
    myth  = omp_get_thread_num()

    n_tmp=1
    call divide_region( &
         myth,numth,n_tmp,n,n_startth,n_endth)

    istart_th = istart -1 + n_startth 
    iend_th   = istart_th + (n_endth - n_startth) 
    n_th      = n_endth - n_startth + 1

    ! -----------
    nnz_th=0
    do i = n_startth,n_endth
       do j = crsMat%crsRow_ptr(i),crsMat%crsRow_ptr(i+1)-1
          Col =  crsMat%crsCol(j)
          if( (istart_th<=Col) .and. (Col<=iend_th) )then
             nnz_th=nnz_th+1
          endif
       enddo
    enddo
    method%preconILU(myth+1)%nnz_th=nnz_th
    call ddinit_crsMat(method%preconILU(myth+1)%crsMatLU,n_th,nnz_th)
    nnz_th=1
    do i = n_startth,n_endth
       i_th = i - n_startth +1
       method%preconILU(myth+1)%crsMatLU%crsRow_ptr(i_th)=nnz_th
       do j = crsMat%crsRow_ptr(i),crsMat%crsRow_ptr(i+1)-1
          Col =  crsMat%crsCol(j)
          if( (istart_th<=Col) .and. (Col<=iend_th) )then
             method%preconILU(myth+1)%crsMatLU%crsA_hi(nnz_th)   = crsMat%crsA_hi(j)
             method%preconILU(myth+1)%crsMatLU%crsA_lo(nnz_th)   = crsMat%crsA_lo(j)
             method%preconILU(myth+1)%crsMatLU%crsCol(nnz_th) = Col
             nnz_th=nnz_th+1
          endif
       enddo
    enddo
    method%preconILU(myth+1)%crsMatLU%crsRow_ptr(n_th+1)=nnz_th

    ! -----------

    method%preconILU(myth+1)%istart_th = istart_th
    method%preconILU(myth+1)%iend_th   = iend_th
    method%preconILU(myth+1)%n_th      = n_th

    !    write(*,*)nnz,myth,nnz_th,n_startth
    !    write(*,*)"START ",myth,n_startth,istart_th
    !    write(*,*)"END   ",myth,n_endth  ,iend_th

    call ddILU_LU_crs( & 
         method%preconILU(myth+1)%n_th,&
         method%preconILU(myth+1)%istart_th,&
         method%preconILU(myth+1)%nnz_th,&
         method%preconILU(myth+1)%crsMatLU%crsA_hi,&
         method%preconILU(myth+1)%crsMatLU%crsA_lo,&
         method%preconILU(myth+1)%crsMatLU%crsRow_ptr,&
         method%preconILU(myth+1)%crsMatLU%crsCol,&
         method%preconILU(myth+1)%crsMatL,&
         method%preconILU(myth+1)%crsMatU)   

    call get_nnz_crs( &
         method%preconILU(myth+1)%crsMatL, &
         method%preconILU(myth+1)%n_th, &
         method%preconILU(myth+1)%nnzL_th )

    call get_nnz_crs( &
         method%preconILU(myth+1)%crsMatU, &
         method%preconILU(myth+1)%n_th, &
         method%preconILU(myth+1)%nnzU_th )

    !$OMP END PARALLEL
    ! ----

    call omp_set_num_threads(max_numth)
#else

    n      = method%prm%n
    istart = method%prm%istart
    precon_thblock=method%precon_thblock 

    ! ----
    ! 計算領域の分割
    call get_maxthreads(max_numth)
    call set_precon_maxthreads(n,precon_thblock)

    call get_maxthreads(numth)
    allocate(method%preconILU(numth))

    !$OMP PARALLEL &
    !$OMP PRIVATE( &
    !$OMP numth,myth,&
    !$OMP n_startth,n_endth,&
    !$OMP istart_th,iend_th,&
    !$OMP n_tmp, &
    !$OMP n_th,nnz_th, &
    !$OMP i,j,i_th,Col,&
    !$OMP ii,iimin,independ_nvec &
    !$OMP )
    numth = omp_get_num_threads()
    myth  = omp_get_thread_num()

    n_tmp=1
    call divide_region( &
         myth,numth,n_tmp,n,n_startth,n_endth)

    istart_th = istart -1 + n_startth 
    iend_th   = istart_th + (n_endth - n_startth) 
    n_th      = n_endth - n_startth + 1

    ! -----------
    nnz_th=0
    independ_nvec = method%independ_nvec
    do ii = n_startth,n_endth , independ_nvec
       iimin = min(ii+independ_nvec-1,n_endth)
       do i = ii,min(ii+independ_nvec-1,n_endth)
          do j = crsMat%crsRow_ptr(i),crsMat%crsRow_ptr(i+1)-1
             Col =  crsMat%crsCol(j)

             if( (istart_th<=Col) .and. (Col<=iend_th) )then
                if( ( Col < (ii + istart - 1) ) .or. ( (iimin + istart - 1) < Col ) )then
                   nnz_th=nnz_th+1
                else
!                   if( Col .eq. (i + istart - 1) )then
                   if( Col .eq. (i + istart - 1) )then
                      nnz_th=nnz_th+1
                   endif
                endif
             endif

          enddo
       enddo
    enddo
    method%preconILU(myth+1)%nnz_th=nnz_th

    call ddinit_crsMat(method%preconILU(myth+1)%crsMatLU,n_th,nnz_th)
    nnz_th=1
    do ii = n_startth,n_endth , independ_nvec
       iimin = min(ii+independ_nvec-1,n_endth)
       do i = ii,min(ii+independ_nvec-1,n_endth)
          i_th = i - n_startth +1
          method%preconILU(myth+1)%crsMatLU%crsRow_ptr(i_th)=nnz_th
          do j = crsMat%crsRow_ptr(i),crsMat%crsRow_ptr(i+1)-1
             Col =  crsMat%crsCol(j)

             if( (istart_th<=Col) .and. (Col<=iend_th) )then
                if( ( Col < (ii + istart - 1) ) .or. ( (iimin + istart - 1) < Col ) )then
!                   method%preconILU(myth+1)%crsMatLU%crsA(nnz_th)   = crsMat%crsA(j) 
                   method%preconILU(myth+1)%crsMatLU%crsA_hi(nnz_th)   = crsMat%crsA_hi(j) 
                   method%preconILU(myth+1)%crsMatLU%crsA_lo(nnz_th)   = crsMat%crsA_lo(j) 
                   method%preconILU(myth+1)%crsMatLU%crsCol(nnz_th) = Col
                   nnz_th=nnz_th+1
                else

!                   if( Col .eq. (i + istart - 1) )then
                   if( Col .eq. (i + istart - 1) )then
!                      method%preconILU(myth+1)%crsMatLU%crsA(nnz_th)   = crsMat%crsA(j) 
                      method%preconILU(myth+1)%crsMatLU%crsA_hi(nnz_th)   = crsMat%crsA_hi(j) 
                      method%preconILU(myth+1)%crsMatLU%crsA_lo(nnz_th)   = crsMat%crsA_lo(j) 
                      method%preconILU(myth+1)%crsMatLU%crsCol(nnz_th) = Col
                      nnz_th=nnz_th+1
                   endif
                endif
             endif

          enddo
       enddo
    enddo
    method%preconILU(myth+1)%crsMatLU%crsRow_ptr(n_th+1)=nnz_th


    ! -----------

    method%preconILU(myth+1)%istart_th = istart_th
    method%preconILU(myth+1)%iend_th   = iend_th
    method%preconILU(myth+1)%n_th      = n_th

    call ddILU_LU_crs( & 
         method%preconILU(myth+1)%n_th,&
         method%preconILU(myth+1)%istart_th,&
         method%preconILU(myth+1)%nnz_th,&
         method%preconILU(myth+1)%crsMatLU%crsA_hi,&
         method%preconILU(myth+1)%crsMatLU%crsA_lo,&
         method%preconILU(myth+1)%crsMatLU%crsRow_ptr,&
         method%preconILU(myth+1)%crsMatLU%crsCol,&
         method%preconILU(myth+1)%crsMatL,&
         method%preconILU(myth+1)%crsMatU)   

    call get_nnz_crs( &
         method%preconILU(myth+1)%crsMatL, &
         method%preconILU(myth+1)%n_th, &
         method%preconILU(myth+1)%nnzL_th )

    call get_nnz_crs( &
         method%preconILU(myth+1)%crsMatU, &
         method%preconILU(myth+1)%n_th, &
         method%preconILU(myth+1)%nnzU_th )

    call ddfree_crsMat(method%preconILU(myth+1)%crsMatLU)

    !$OMP END PARALLEL
    call omp_set_num_threads(max_numth)

    ! ----

#endif

  end subroutine ddget_ILU_LU_crs

  subroutine ddget_ILU_LU_dia(method,diaMat)

    use krylov
    use omp_lib
    implicit none
    type(KSP) :: method
    type(DIA_Mat) :: diaMat 

    integer max_numth
    integer numth,myth
    integer(kind=intType) n,istart
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th
    integer(kind=intType) n_tmp
    integer(kind=intType) n_th,nnz_th
    integer(kind=intType) i,j,i_th
    integer(kind=intType) nnd
    integer(kind=intType) Col,Row
    integer(kind=intType) ii,iimin
    integer independ_nvec
    integer precon_thblock

#if 1

    n      = method%prm%n
    istart = method%prm%istart
    nnd    = method%prm%nnd
    precon_thblock=method%precon_thblock 

    ! ----
    ! 計算領域の分割
    call get_maxthreads(max_numth)
    call set_precon_maxthreads(n,precon_thblock)

    !$OMP PARALLEL 
    !$OMP MASTER
    numth = omp_get_num_threads()
    !$OMP END MASTER
    !$OMP END PARALLEL 
    allocate(method%preconILU(numth))

    !$OMP PARALLEL &
    !$OMP PRIVATE( &
    !$OMP numth,myth,&
    !$OMP n_startth,n_endth,&
    !$OMP istart_th,iend_th,&
    !$OMP n_tmp, &
    !$OMP n_th, &
    !$OMP i,j,i_th,Col,Row, &
    !$OMP ii,iimin,independ_nvec &
    !$OMP )
    numth = omp_get_num_threads()
    myth  = omp_get_thread_num()

    n_tmp=1
    call divide_region( &
         myth,numth,n_tmp,n,n_startth,n_endth)

    istart_th = istart -1 + n_startth 
    iend_th   = istart_th + (n_endth - n_startth) 
    n_th      = n_endth - n_startth + 1

    ! -----------
    ! 各スレッドの領域ごとに係数行列を作成
    ! 対角本数は固定で領域外の数値をゼロとして扱う
#if 0
    allocate(method%preconILU(myth+1)%diaMatLU%diaA(nnd*n_th))
    allocate(method%preconILU(myth+1)%diaMatLU%offset(nnd))
#else
    call ddinit_diaMat(&
         method%preconILU(myth+1)%diaMatLU,n_th,nnd)
#endif

    do i=1,nnd
       method%preconILU(myth+1)%diaMatLU%offset(i) = diaMat%offset(i)
    enddo
    independ_nvec = method%independ_nvec

    do j=1,nnd       
       do ii = n_startth,n_endth , independ_nvec
          iimin = min(ii+independ_nvec-1,n_endth)
          do i = ii,min(ii+independ_nvec-1,n_endth)
             Row = i-1+istart
             Col = Row+diaMat%offset(j)
             ! 特定の範囲内をゼロにする
             if( (istart_th<=Col).and.(Col<=iend_th ) )then
                if( ( Col < (ii + istart - 1) ) .or. ( (iimin + istart - 1) < Col ) )then
                   method%preconILU(myth+1)%diaMatLU%diaA_hi(i-n_startth+1+(j-1)*n_th) = diaMat%diaA_hi(i+(j-1)*n)
                   method%preconILU(myth+1)%diaMatLU%diaA_lo(i-n_startth+1+(j-1)*n_th) = diaMat%diaA_lo(i+(j-1)*n)
                else
                   if( Col .eq. (i + istart - 1) )then
                      method%preconILU(myth+1)%diaMatLU%diaA_hi(i-n_startth+1+(j-1)*n_th) = diaMat%diaA_hi(i+(j-1)*n)
                      method%preconILU(myth+1)%diaMatLU%diaA_lo(i-n_startth+1+(j-1)*n_th) = diaMat%diaA_lo(i+(j-1)*n)
                   else
                      method%preconILU(myth+1)%diaMatLU%diaA_hi(i-n_startth+1+(j-1)*n_th) = 0.0d0
                      method%preconILU(myth+1)%diaMatLU%diaA_lo(i-n_startth+1+(j-1)*n_th) = 0.0d0
                   endif
                endif
             else
                method%preconILU(myth+1)%diaMatLU%diaA_hi(i-n_startth+1+(j-1)*n_th) = 0.0d0
                method%preconILU(myth+1)%diaMatLU%diaA_lo(i-n_startth+1+(j-1)*n_th) = 0.0d0
             endif
          enddo
       enddo
    enddo

    method%preconILU(myth+1)%istart_th = istart_th
    method%preconILU(myth+1)%iend_th   = iend_th
    method%preconILU(myth+1)%n_th      = n_th
    ! -----------
    call ddILU_LU_dia(&
         method%preconILU(myth+1)%n_th,&
         method%preconILU(myth+1)%istart_th,&
         nnd,&
         method%preconILU(myth+1)%diaMatLU%diaA_hi , &
         method%preconILU(myth+1)%diaMatLU%diaA_lo , &
         method%preconILU(myth+1)%diaMatLU%offset , &
         method%preconILU(myth+1)%diaMatL , &
         method%preconILU(myth+1)%diaMatU )

    method%preconILU(myth+1)%nnd_th  = nnd
    method%preconILU(myth+1)%nndL_th =  &
         size(method%preconILU(myth+1)%diaMatL%diaA_hi) / &
         method%preconILU(myth+1)%n_th
    method%preconILU(myth+1)%nndU_th = &
         size(method%preconILU(myth+1)%diaMatU%diaA_hi) / &
         method%preconILU(myth+1)%n_th
    call ddfree_diaMat(method%preconILU(myth+1)%diaMatLU)

    !$OMP END PARALLEL
    ! ----
    call omp_set_num_threads(max_numth)

#else
    n      = method%prm%n
    istart = method%prm%istart
    nnd    = method%prm%nnd

    ! ----
    ! 計算領域の分割
    call get_maxthreads(max_numth)
!    call set_precon_maxthreads0(n)
    call set_precon_maxthreads(n,method%precon_thblock)

    !$OMP PARALLEL 
    !$OMP MASTER
    numth = omp_get_num_threads()
    !$OMP END MASTER
    !$OMP END PARALLEL 
    allocate(method%preconILU(numth))

    !$OMP PARALLEL &
    !$OMP PRIVATE( &
    !$OMP numth,myth,&
    !$OMP n_startth,n_endth,&
    !$OMP istart_th,iend_th,&
    !$OMP n_tmp, &
    !$OMP n_th, &
    !$OMP i,j,i_th,Col,Row &
    !$OMP )
    numth = omp_get_num_threads()
    myth  = omp_get_thread_num()

    n_tmp=1
    call divide_region( &
         myth,numth,n_tmp,n,n_startth,n_endth)

    istart_th = istart -1 + n_startth 
    iend_th   = istart_th + (n_endth - n_startth) 
    n_th      = n_endth - n_startth + 1

    ! -----------
    ! 各スレッドの領域ごとに係数行列を作成
    ! 対角本数は固定で領域外の数値をゼロとして扱う
#if 0
    allocate(method%preconILU(myth+1)%diaMatLU%diaA_hi(nnd*n_th))
    allocate(method%preconILU(myth+1)%diaMatLU%diaA_lo(nnd*n_th))
    allocate(method%preconILU(myth+1)%diaMatLU%offset(nnd))
#else
    call ddinit_diaMat(&
         method%preconILU(myth+1)%diaMatLU,n_th,nnd)
#endif

    do i=1,nnd
       method%preconILU(myth+1)%diaMatLU%offset(i) = diaMat%offset(i)
    enddo
    do j=1,nnd
       do i=n_startth,n_endth
          Row = i-1+istart
          Col = Row+diaMat%offset(j)
          if( (istart_th<=Col).and.(Col<=iend_th ) )then
             method%preconILU(myth+1)%diaMatLU%diaA_hi(i-n_startth+1+(j-1)*n_th) = diaMat%diaA_hi(i+(j-1)*n)
             method%preconILU(myth+1)%diaMatLU%diaA_lo(i-n_startth+1+(j-1)*n_th) = diaMat%diaA_lo(i+(j-1)*n)
          else
             method%preconILU(myth+1)%diaMatLU%diaA_hi(i-n_startth+1+(j-1)*n_th) = 0.0d0
             method%preconILU(myth+1)%diaMatLU%diaA_lo(i-n_startth+1+(j-1)*n_th) = 0.0d0
          endif
       enddo
    enddo
    method%preconILU(myth+1)%istart_th = istart_th
    method%preconILU(myth+1)%iend_th   = iend_th
    method%preconILU(myth+1)%n_th      = n_th
    ! -----------
    call ddILU_LU_dia(&
         method%preconILU(myth+1)%n_th,&
         method%preconILU(myth+1)%istart_th,&
         nnd,&
         method%preconILU(myth+1)%diaMatLU%diaA_hi , &
         method%preconILU(myth+1)%diaMatLU%diaA_lo , &
         method%preconILU(myth+1)%diaMatLU%offset , &
         method%preconILU(myth+1)%diaMatL , &
         method%preconILU(myth+1)%diaMatU )
 
    method%preconILU(myth+1)%nnd_th  = nnd
    method%preconILU(myth+1)%nndL_th =  &
         size(method%preconILU(myth+1)%diaMatL%diaA_hi) / &
         method%preconILU(myth+1)%n_th
    method%preconILU(myth+1)%nndU_th = &
         size(method%preconILU(myth+1)%diaMatU%diaA_hi) / &
         method%preconILU(myth+1)%n_th

    !$OMP END PARALLEL
    ! ----

    call omp_set_num_threads(max_numth)
#endif

  end subroutine ddget_ILU_LU_dia
  
  subroutine ddsolve_LU_crs( &
       n,istart,&
       nnzL,crsAL_hi,crsAL_lo,&
       crsRow_ptrL,crsColL, &
       nnzU,crsAU_hi,crsAU_lo,&
       crsRow_ptrU,crsColU, &
       r_hi,r_lo,&
       s_hi,s_lo,&
       independ_nvec, &
       PRECISION_PRECONDITION &
       )
    use PARCEL_integer
    implicit none
    integer(kind=intType) n,istart

    integer(kind=intType) nnzL
    real*8 crsAL_hi(nnzL),crsAL_lo(nnzL)
    integer(kind=intType) crsRow_ptrL(n+1)
    integer(kind=intType) crsColL(nnzL)

    integer(kind=intType) nnzU
    real*8 crsAU_hi(nnzU),crsAU_lo(nnzU)
    integer(kind=intType) crsRow_ptrU(n+1)
    integer(kind=intType) crsColU(nnzU)

    real*8 r_hi(n),r_lo(n)
    real*8 s_hi(n),s_lo(n)
    integer independ_nvec
    integer PRECISION_PRECONDITION

    real*8 z_hi,z_lo
    real*8 zz_hi(independ_nvec),zz_lo(independ_nvec)

    real*8 c_hi,c_lo
    real*8 ddmul_var1
    real*8 ddsub_var1
    real*8 dddiv_var1
    integer(kind=intType) i,j
    integer(kind=intType) ii0,ii,jj
    integer(kind=intType) iistep,iimin
    integer(kind=intType) Col,Row

    select case (PRECISION_PRECONDITION)

    case(KSP_DOUBLE_PRECISION)

       if(independ_nvec.eq.1)then
          ! -----
          ! z = Us 
          ! Lz = r をzについて解く
          ! 前進代入
          do i=1,n
             z_hi=r_hi(i)
             z_lo=r_lo(i)
             ! 対角成分を持つことを仮定
             ! 下三角行列では各行の最後の値が対角成分
             do j=crsRow_ptrL(i),crsRow_ptrL(i+1)-1 -1
                Col=crsColL(j)-istart+1
                !             z=z-s(Col)*crsAL(j)
#ifdef Manual_inline
                ddmul_calc1(s_hi(Col),s_lo(Col),crsAL_hi(j),0.0d0,c_hi,c_lo)
                ddsub_calc1(z_hi,z_lo,c_hi,c_lo,z_hi,z_lo)
#else
                call dd_mul( &
                     s_hi(Col),s_lo(Col),&
                     crsAL_hi(j),0.0d0,&
                     c_hi,c_lo)
                call dd_sub( &
                     z_hi,z_lo, &
                     c_hi,c_lo, &
                     z_hi,z_lo)
#endif
             enddo
             j=crsRow_ptrL(i+1)-1
             ! s(i)=z/crsAL(j)
#ifdef Manual_inline
             dddiv_calc1(z_hi,z_lo,crsAL_hi(j),0.0d0,s_hi(i),s_lo(i))
#else
             call dd_div( &
                  z_hi,z_lo, &
                  crsAL_hi(j),0.0d0,&
                  s_hi(i),s_lo(i))
#endif
             !       write(*,*)"L^-1 r = ",i,s(i)
          enddo
          !    write(*,*)"--------"
          ! -----
          ! Us=z をsについて解く
          ! 後退代入
          do i=n,1,-1
             z_hi=s_hi(i)
             z_lo=s_lo(i)
             ! 対角成分を持つことを仮定
             ! 上三角行列では最初の値が対角成分
             do j=crsRow_ptrU(i)+1,crsRow_ptrU(i+1)-1
                Col=crsColU(j)-istart+1
                !             z=z-s(Col)*crsAU(j)
#ifdef Manual_inline
                ddmul_calc1(s_hi(Col),s_lo(Col),crsAU_hi(j),0.0d0,c_hi,c_lo)
                ddsub_calc1(z_hi,z_lo,c_hi,c_lo,z_hi,z_lo)
#else
                call dd_mul( &
                     s_hi(Col),s_lo(Col),&
                     crsAU_hi(j),0.0d0,&
                     c_hi,c_lo)
                call dd_sub( &
                     z_hi,z_lo, &
                     c_hi,c_lo, &
                     z_hi,z_lo)
#endif
             enddo
             j=crsRow_ptrU(i)
             ! s(i)=z/crsAU(j)
#ifdef Manual_inline
             dddiv_calc1(z_hi,z_lo,crsAU_hi(j),0.0d0,s_hi(i),s_lo(i))
#else
             call dd_div( &
                  z_hi,z_lo, &
                  crsAU_hi(j),0.0d0,&
                  s_hi(i),s_lo(i))
#endif
             !       write(*,*)"U^-1 z = ",i,s(i)
          enddo
       else
          ! -----
          ! z = Us 
          ! Lz = r をzについて解く
          ! 前進代入
          !       do i=1,n
          do ii = 1,n,independ_nvec
             iimin = min(ii+independ_nvec-1,n)
#if 1
             ii0 = ii
             do i = ii ,iimin
                !             zz(i-ii+1)=r(i)
                zz_hi(i-ii+1)=r_hi(i)
                zz_lo(i-ii+1)=r_lo(i)
             enddo
             do 
                iistep = 0
                do i = ii0,iimin-1
                   if ( &
                        ( (crsRow_ptrL(i+1)-1 -1)   - (crsRow_ptrL(i)) + 1 ) &
                        .eq. &
                        ( (crsRow_ptrL(i+1+1)-1 -1) - (crsRow_ptrL(i+1)) + 1 ) &
                        )then
                      iistep=iistep+1
                   else
                      exit
                   endif
                enddo
                do jj = 1,(crsRow_ptrL(i+1)-1 -1) - (crsRow_ptrL(i)) + 1 
                   do i = ii0,ii0+iistep
                      j = crsRow_ptrL(i)+jj-1  
                      Col=crsColL(j)-istart+1
                      ! zz(i-ii+1)=zz(i-ii+1)-s(Col)*crsAL(j)
#ifdef Manual_inline
                      ddmul_calc1(s_hi(Col),s_lo(Col),crsAL_hi(j),0.0d0,c_hi,c_lo)
                      ddsub_calc1(zz_hi(i-ii+1),zz_lo(i-ii+1),c_hi,c_lo,zz_hi(i-ii+1),zz_lo(i-ii+1))
#else
                      call dd_mul( &
                           s_hi(Col),s_lo(Col),&
                           crsAL_hi(j),0.0d0,&
                           c_hi,c_lo)
                      call dd_sub( &
                           zz_hi(i-ii+1),zz_lo(i-ii+1), &
                           c_hi,c_lo, &
                           zz_hi(i-ii+1),zz_lo(i-ii+1))                   
#endif
                   enddo
                enddo
                ! do i = ii0,ii0+iistep
                !    write(*,*)"LU solver",i
                ! enddo
                ! write(*,*)"------------------"

                ii0=ii0+iistep+1
                if(ii0>iimin)then
                   exit
                endif
             enddo
             do i = ii ,iimin
                j=crsRow_ptrL(i+1)-1
                !             s(i)=zz(i-ii+1)/crsAL(j)
#ifdef Manual_inline
                dddiv_calc1(zz_hi(i-ii+1),zz_lo(i-ii+1),crsAL_hi(j),0.0d0,s_hi(i),s_lo(i))
#else
                call dd_div( &
                     zz_hi(i-ii+1),zz_lo(i-ii+1), &
                     crsAL_hi(j),0.0d0,&
                     s_hi(i),s_lo(i))
#endif
             enddo
#else
             do i = ii ,iimin
                zz(i-ii+1)=r(i)
                ! 対角成分を持つことを仮定
                ! 下三角行列では各行の最後の値が対角成分
                do j=crsRow_ptrL(i),crsRow_ptrL(i+1)-1 -1
                   Col=crsColL(j)-istart+1
                   zz(i-ii+1)=zz(i-ii+1)-s(Col)*crsAL(j)
                enddo
                j=crsRow_ptrL(i+1)-1
                s(i)=zz(i-ii+1)/crsAL(j)
             enddo
#endif

          enddo

          ! -----
          ! Us=z をsについて解く
          ! 後退代入
#if 1
          do ii = n,1,-independ_nvec
             ii0 = ii
             iimin = max(ii-independ_nvec+1,1) 
             do i = ii ,iimin ,-1
                !             zz(ii-i+1)=s(i)
                zz_hi(ii-i+1)=s_hi(i)
                zz_lo(ii-i+1)=s_lo(i)
             enddo
             do 
                iistep = 0
                do i = ii0,iimin+1,-1              
                   if ( &
                        ( (crsRow_ptrU(i+1)-1)   - (crsRow_ptrU(i)+1)   + 1 ) &
                        .eq. &
                        ( (crsRow_ptrU(i+1-1)-1) - (crsRow_ptrU(i-1)+1) + 1 ) &
                        )then
                      iistep=iistep+1
                   else
                      exit
                   endif
                enddo
                do jj = 1,(crsRow_ptrU(i+1)-1) - (crsRow_ptrU(i)+1) +1
                   do i = ii0,ii0-iistep,-1
                      j = crsRow_ptrU(i)+1 + jj - 1
                      Col=crsColU(j)-istart+1
                      !                   zz(ii-i+1)=zz(ii-i+1)-s(Col)*crsAU(j)
#ifdef Manual_inline
                      ddmul_calc1(s_hi(Col),s_lo(Col),crsAU_hi(j),0.0d0,c_hi,c_lo)
                      ddsub_calc1(zz_hi(ii-i+1),zz_lo(ii-i+1),c_hi,c_lo,zz_hi(ii-i+1),zz_lo(ii-i+1))
#else
                      call dd_mul( &
                           s_hi(Col),s_lo(Col),&
                           crsAU_hi(j),0.0d0,&
                           c_hi,c_lo)
                      call dd_sub( &
                           zz_hi(ii-i+1),zz_lo(ii-i+1), &
                           c_hi,c_lo, &
                           zz_hi(ii-i+1),zz_lo(ii-i+1))
#endif
                   enddo
                enddo

                ii0 = ii0-iistep-1
                if(ii0<iimin)then
                   exit
                endif

             enddo

             do i = ii ,iimin ,-1
                j=crsRow_ptrU(i)
                !             s(i)=zz(ii-i+1)/crsAU(j)
#ifdef Manual_inline
                dddiv_calc1(zz_hi(ii-i+1),zz_lo(ii-i+1),crsAU_hi(j),0.0d0,s_hi(i),s_lo(i))
#else
                call dd_div( &
                     zz_hi(ii-i+1),zz_lo(ii-i+1), &
                     crsAU_hi(j),0.0d0,&
                     s_hi(i),s_lo(i))
#endif
             enddo

          enddo

#else
          do ii = n,1,-independ_nvec
             iimin = max(ii-independ_nvec+1,1) 
             do i = ii ,iimin ,-1
                zz(ii-i+1)=s(i)
                ! 対角成分を持つことを仮定
                ! 上三角行列では最初の値が対角成分
                do j=crsRow_ptrU(i)+1,crsRow_ptrU(i+1)-1
                   Col=crsColU(j)-istart+1
                   zz(ii-i+1)=zz(ii-i+1)-s(Col)*crsAU(j)
                enddo
                j=crsRow_ptrU(i)
                s(i)=zz(ii-i+1)/crsAU(j)
                !       write(*,*)"U^-1 z = ",i,s(i)
             enddo
          enddo
#endif

       endif

    case(KSP_DOUBLE_DOUBLE_PRECISION)

       if(independ_nvec.eq.1)then
          ! -----
          ! z = Us 
          ! Lz = r をzについて解く
          ! 前進代入
          do i=1,n
             z_hi=r_hi(i)
             z_lo=r_lo(i)
             ! 対角成分を持つことを仮定
             ! 下三角行列では各行の最後の値が対角成分
             do j=crsRow_ptrL(i),crsRow_ptrL(i+1)-1 -1
                Col=crsColL(j)-istart+1
                !             z=z-s(Col)*crsAL(j)
#ifdef Manual_inline
                ddmul_calc1(s_hi(Col),s_lo(Col),crsAL_hi(j),crsAL_lo(j),c_hi,c_lo)
                ddsub_calc1(z_hi,z_lo,c_hi,c_lo,z_hi,z_lo)
#else
                call dd_mul( &
                     s_hi(Col),s_lo(Col),&
                     crsAL_hi(j), crsAL_lo(j),&
                     c_hi,c_lo)
                call dd_sub( &
                     z_hi,z_lo, &
                     c_hi,c_lo, &
                     z_hi,z_lo)
#endif
             enddo
             j=crsRow_ptrL(i+1)-1
             ! s(i)=z/crsAL(j)
#ifdef Manual_inline
             dddiv_calc1(z_hi,z_lo,crsAL_hi(j),crsAL_lo(j),s_hi(i),s_lo(i))
#else
             call dd_div( &
                  z_hi,z_lo, &
                  crsAL_hi(j),crsAL_lo(j),&
                  s_hi(i),s_lo(i))
#endif
             !       write(*,*)"L^-1 r = ",i,s(i)
          enddo
          !    write(*,*)"--------"
          ! -----
          ! Us=z をsについて解く
          ! 後退代入
          do i=n,1,-1
             z_hi=s_hi(i)
             z_lo=s_lo(i)
             ! 対角成分を持つことを仮定
             ! 上三角行列では最初の値が対角成分
             do j=crsRow_ptrU(i)+1,crsRow_ptrU(i+1)-1
                Col=crsColU(j)-istart+1
                !             z=z-s(Col)*crsAU(j)
#ifdef Manual_inline
                ddmul_calc1(s_hi(Col),s_lo(Col),crsAU_hi(j),crsAU_lo(j),c_hi,c_lo)
                ddsub_calc1(z_hi,z_lo,c_hi,c_lo,z_hi,z_lo)
#else
                call dd_mul( &
                     s_hi(Col),s_lo(Col),&
                     crsAU_hi(j), crsAU_lo(j),&
                     c_hi,c_lo)
                call dd_sub( &
                     z_hi,z_lo, &
                     c_hi,c_lo, &
                     z_hi,z_lo)
#endif
             enddo
             j=crsRow_ptrU(i)
             ! s(i)=z/crsAU(j)
#ifdef Manual_inline
             dddiv_calc1(z_hi,z_lo,crsAU_hi(j),crsAU_lo(j),s_hi(i),s_lo(i))
#else
             call dd_div( &
                  z_hi,z_lo, &
                  crsAU_hi(j),crsAU_lo(j),&
                  s_hi(i),s_lo(i))
#endif
             !       write(*,*)"U^-1 z = ",i,s(i)
          enddo
       else
          ! -----
          ! z = Us 
          ! Lz = r をzについて解く
          ! 前進代入
          !       do i=1,n
          do ii = 1,n,independ_nvec
             iimin = min(ii+independ_nvec-1,n)
#if 1
             ii0 = ii
             do i = ii ,iimin
                !             zz(i-ii+1)=r(i)
                zz_hi(i-ii+1)=r_hi(i)
                zz_lo(i-ii+1)=r_lo(i)
             enddo
             do 
                iistep = 0
                do i = ii0,iimin-1
                   if ( &
                        ( (crsRow_ptrL(i+1)-1 -1)   - (crsRow_ptrL(i)) + 1 ) &
                        .eq. &
                        ( (crsRow_ptrL(i+1+1)-1 -1) - (crsRow_ptrL(i+1)) + 1 ) &
                        )then
                      iistep=iistep+1
                   else
                      exit
                   endif
                enddo
                do jj = 1,(crsRow_ptrL(i+1)-1 -1) - (crsRow_ptrL(i)) + 1 
                   do i = ii0,ii0+iistep
                      j = crsRow_ptrL(i)+jj-1  
                      Col=crsColL(j)-istart+1
                      ! zz(i-ii+1)=zz(i-ii+1)-s(Col)*crsAL(j)
#ifdef Manual_inline
                      ddmul_calc1(s_hi(Col),s_lo(Col),crsAL_hi(j),crsAL_lo(j),c_hi,c_lo)
                      ddsub_calc1(zz_hi(i-ii+1),zz_lo(i-ii+1),c_hi,c_lo,zz_hi(i-ii+1),zz_lo(i-ii+1))
#else
                      call dd_mul( &
                           s_hi(Col),s_lo(Col),&
                           crsAL_hi(j),crsAL_lo(j),&
                           c_hi,c_lo)
                      call dd_sub( &
                           zz_hi(i-ii+1),zz_lo(i-ii+1), &
                           c_hi,c_lo, &
                           zz_hi(i-ii+1),zz_lo(i-ii+1))                   
#endif
                   enddo
                enddo
                ! do i = ii0,ii0+iistep
                !    write(*,*)"LU solver",i
                ! enddo
                ! write(*,*)"------------------"

                ii0=ii0+iistep+1
                if(ii0>iimin)then
                   exit
                endif
             enddo
             do i = ii ,iimin
                j=crsRow_ptrL(i+1)-1
                !             s(i)=zz(i-ii+1)/crsAL(j)
#ifdef Manual_inline
                dddiv_calc1(zz_hi(i-ii+1),zz_lo(i-ii+1),crsAL_hi(j),crsAL_lo(j),s_hi(i),s_lo(i))
#else
                call dd_div( &
                     zz_hi(i-ii+1),zz_lo(i-ii+1), &
                     crsAL_hi(j),crsAL_lo(j),&
                     s_hi(i),s_lo(i))
#endif
             enddo
#else
             do i = ii ,iimin
                zz(i-ii+1)=r(i)
                ! 対角成分を持つことを仮定
                ! 下三角行列では各行の最後の値が対角成分
                do j=crsRow_ptrL(i),crsRow_ptrL(i+1)-1 -1
                   Col=crsColL(j)-istart+1
                   zz(i-ii+1)=zz(i-ii+1)-s(Col)*crsAL(j)
                enddo
                j=crsRow_ptrL(i+1)-1
                s(i)=zz(i-ii+1)/crsAL(j)
             enddo
#endif

          enddo

          ! -----
          ! Us=z をsについて解く
          ! 後退代入
#if 1
          do ii = n,1,-independ_nvec
             ii0 = ii
             iimin = max(ii-independ_nvec+1,1) 
             do i = ii ,iimin ,-1
                !             zz(ii-i+1)=s(i)
                zz_hi(ii-i+1)=s_hi(i)
                zz_lo(ii-i+1)=s_lo(i)
             enddo
             do 
                iistep = 0
                do i = ii0,iimin+1,-1              
                   if ( &
                        ( (crsRow_ptrU(i+1)-1)   - (crsRow_ptrU(i)+1)   + 1 ) &
                        .eq. &
                        ( (crsRow_ptrU(i+1-1)-1) - (crsRow_ptrU(i-1)+1) + 1 ) &
                        )then
                      iistep=iistep+1
                   else
                      exit
                   endif
                enddo
                do jj = 1,(crsRow_ptrU(i+1)-1) - (crsRow_ptrU(i)+1) +1
                   do i = ii0,ii0-iistep,-1
                      j = crsRow_ptrU(i)+1 + jj - 1
                      Col=crsColU(j)-istart+1
                      !                   zz(ii-i+1)=zz(ii-i+1)-s(Col)*crsAU(j)
#ifdef Manual_inline
                      ddmul_calc1(s_hi(Col),s_lo(Col),crsAU_hi(j),crsAU_lo(j),c_hi,c_lo)
                      ddsub_calc1(zz_hi(ii-i+1),zz_lo(ii-i+1),c_hi,c_lo,zz_hi(ii-i+1),zz_lo(ii-i+1))
#else
                      call dd_mul( &
                           s_hi(Col),s_lo(Col),&
                           crsAU_hi(j),crsAU_lo(j),&
                           c_hi,c_lo)
                      call dd_sub( &
                           zz_hi(ii-i+1),zz_lo(ii-i+1), &
                           c_hi,c_lo, &
                           zz_hi(ii-i+1),zz_lo(ii-i+1))
#endif
                   enddo
                enddo

                ii0 = ii0-iistep-1
                if(ii0<iimin)then
                   exit
                endif

             enddo

             do i = ii ,iimin ,-1
                j=crsRow_ptrU(i)
                !             s(i)=zz(ii-i+1)/crsAU(j)
#ifdef Manual_inline
                dddiv_calc1(zz_hi(ii-i+1),zz_lo(ii-i+1),crsAU_hi(j),crsAU_lo(j),s_hi(i),s_lo(i))
#else
                call dd_div( &
                     zz_hi(ii-i+1),zz_lo(ii-i+1), &
                     crsAU_hi(j),crsAU_lo(j),&
                     s_hi(i),s_lo(i))
#endif
             enddo

          enddo

#else
          do ii = n,1,-independ_nvec
             iimin = max(ii-independ_nvec+1,1) 
             do i = ii ,iimin ,-1
                zz(ii-i+1)=s(i)
                ! 対角成分を持つことを仮定
                ! 上三角行列では最初の値が対角成分
                do j=crsRow_ptrU(i)+1,crsRow_ptrU(i+1)-1
                   Col=crsColU(j)-istart+1
                   zz(ii-i+1)=zz(ii-i+1)-s(Col)*crsAU(j)
                enddo
                j=crsRow_ptrU(i)
                s(i)=zz(ii-i+1)/crsAU(j)
                !       write(*,*)"U^-1 z = ",i,s(i)
             enddo
          enddo
#endif

       endif

    end select

  end subroutine ddsolve_LU_crs

  subroutine ddsolve_LU_dia( &
       n,istart,&
       nndL,&
       diaAL_hi,diaAL_lo,&       
       offsetL, &
       nndU,&
       diaAU_hi,diaAU_lo,&
       offsetU, &
       r_hi,r_lo, &       
       s_hi,s_lo, &       
       independ_nvec, &
       PRECISION_PRECONDITION &
       )
    use PARCEL_integer
    implicit none
    integer(kind=intType) n,istart

    integer(kind=intType) nndL
    real*8 diaAL_hi(n*nndL)
    real*8 diaAL_lo(n*nndL)
    integer(kind=intType) offsetL(nndL)

    integer(kind=intType) nndU
    real*8 diaAU_hi(n*nndU)
    real*8 diaAU_lo(n*nndU)
    integer(kind=intType) offsetU(nndU)

    real*8 r_hi(n),r_lo(n)
    real*8 s_hi(n),s_lo(n)

    real*8 c_hi,c_lo
    integer independ_nvec
    integer PRECISION_PRECONDITION

    real*8 z_hi,z_lo
    real*8 zz_hi(independ_nvec),zz_lo(independ_nvec)
    real*8 ddmul_var1
    real*8 dddiv_var1
    real*8 ddsub_var1
    integer(kind=intType) i,j
    integer(kind=intType) ii,iimin
    integer(kind=intType) Col,Row

    select case (PRECISION_PRECONDITION)

    case(KSP_DOUBLE_PRECISION)

       ! 零以外の対角成分を持つことを仮定
       ! -----
       ! LUs = r をsについて解いていく
       if(independ_nvec.eq.1)then
          ! -----
          ! z = Us 
          ! Lz = r をzについて解く
          ! 前進代入
          do i=1,n
             z_hi=r_hi(i)
             z_lo=r_lo(i)
             ! 対角成分を持つことを仮定
             ! 下三角行列で対角成分はnndLのオフセット
             do j=1,nndL-1
                if( (1 <= (i+offsetL(j))).and.((i+offsetL(j)) <= n ) )then
                   !                z=z-s(i+offsetL(j))*diaAL(i+(j-1)*n)
#ifdef Manual_inline
                   ddmul_calc1(s_hi(i+offsetL(j)),s_lo(i+offsetL(j)),diaAL_hi(i+(j-1)*n),0.0d0,c_hi,c_lo)
                   ddsub_calc1(z_hi,z_lo,c_hi,c_lo,z_hi,z_lo)
#else
                   call dd_mul( &
                        s_hi(i+offsetL(j)),s_lo(i+offsetL(j)), &
                        diaAL_hi(i+(j-1)*n),0.0d0, &
                        c_hi,c_lo)
                   call dd_sub( &
                        z_hi,z_lo, &
                        c_hi,c_lo,&
                        z_hi,z_lo)
#endif
                endif
             enddo
             j=nndL
             !          s(i)=z/diaAL(i+(j-1)*n)
#ifdef Manual_inline
             dddiv_calc1(z_hi,z_lo,diaAL_hi(i+(j-1)*n),0.0d0,s_hi(i),s_lo(i))
#else
             call dd_div( &
                  z_hi,z_lo,&
                  diaAL_hi(i+(j-1)*n) ,0.0d0 , &
                  s_hi(i),s_lo(i))
#endif
          enddo

          ! -----
          ! Us=z をsについて解く
          ! 後退代入
          do i=n,1,-1
             ! z=s(i)
             z_hi=s_hi(i)
             z_lo=s_lo(i)
             ! 対角成分を持つことを仮定
             ! 上三角行列の対角成分は1のオフセット
             do j=1+1,nndU
                if( ( ( 1 <= (i+offsetU(j)) ) ).and.( (i+offsetU(j) <= n ) ) )then
                   ! z=z-s(i+offsetU(j))*diaAU(i+(j-1)*n)
#ifdef Manual_inline
                   ddmul_calc1(s_hi(i+offsetU(j)),s_lo(i+offsetU(j)),diaAU_hi(i+(j-1)*n),0.0d0,c_hi,c_lo)
                   ddsub_calc1(z_hi,z_lo,c_hi,c_lo,z_hi,z_lo)
#else
                   call dd_mul( &
                        s_hi(i+offsetU(j)),s_lo(i+offsetU(j)) ,&
                        diaAU_hi(i+(j-1)*n),0.0d0, &
                        c_hi,c_lo)
                   call dd_sub( &
                        z_hi,z_lo, &
                        c_hi,c_lo,&
                        z_hi,z_lo)
#endif
                endif
             enddo
             j=1
             ! s(i)=z/diaAU(i+(j-1)*n)
#ifdef Manual_inline
             dddiv_calc1(z_hi,z_lo,diaAU_hi(i+(j-1)*n),0.0d0,s_hi(i),s_lo(i))
#else
             call dd_div( &
                  z_hi,z_lo,&
                  diaAU_hi(i+(j-1)*n) ,0.0d0 , &
                  s_hi(i),s_lo(i))
#endif
          enddo

       else
#if 1
          ! -----
          ! z = Us 
          ! Lz = r をzについて解く
          ! 前進代入
          !       do i=1,n
          do ii = 1,n,independ_nvec
             iimin = min(ii+independ_nvec-1,n)
             do i = ii ,iimin
                ! zz(i-ii+1)=r(i)
                zz_hi(i-ii+1)=r_hi(i)
                zz_lo(i-ii+1)=r_lo(i)
             enddo
             ! 対角成分を持つことを仮定
             ! 下三角行列で対角成分はnndLのオフセット
             do j=1,nndL-1
                do i = ii ,iimin
                   if( (1 <= (i+offsetL(j))).and.((i+offsetL(j)) <= n ) )then
                      ! zz(i-ii+1)=zz(i-ii+1)-s(i+offsetL(j))*diaAL(i+(j-1)*n)
#ifdef Manual_inline
                      ddmul_calc1(s_hi(i+offsetL(j)),s_lo(i+offsetL(j)),diaAL_hi(i+(j-1)*n),0.0d0,c_hi,c_lo)
                      ddsub_calc1(zz_hi(i-ii+1),zz_lo(i-ii+1),c_hi,c_lo,zz_hi(i-ii+1),zz_lo(i-ii+1))
#else
                      call dd_mul( &
                           s_hi(i+offsetL(j)),s_lo(i+offsetL(j)), &
                           diaAL_hi(i+(j-1)*n), 0.0d0, &
                           c_hi,c_lo)
                      call dd_sub( &
                           zz_hi(i-ii+1),zz_lo(i-ii+1), &
                           c_hi,c_lo,&
                           zz_hi(i-ii+1),zz_lo(i-ii+1))
#endif
                   endif
                enddo
             enddo
             do i = ii ,iimin
                j=nndL
                ! s(i)=zz(i-ii+1)/diaAL(i+(j-1)*n)
#ifdef Manual_inline
                dddiv_calc1(zz_hi(i-ii+1),zz_lo(i-ii+1),diaAL_hi(i+(j-1)*n),0.0d0,s_hi(i),s_lo(i))
#else
                call dd_div( &
                     zz_hi(i-ii+1),zz_lo(i-ii+1),&
                     diaAL_hi(i+(j-1)*n) ,0.0d0 , &
                     s_hi(i),s_lo(i))
#endif
             enddo
          enddo

          ! -----
          ! Us=z をsについて解く
          ! 後退代入
          do ii = n,1,-independ_nvec
             iimin = max(ii-independ_nvec+1,1) 
             do i = ii ,iimin ,-1
                ! zz(ii-i+1)=s(i)
                zz_hi(ii-i+1)=s_hi(i)
                zz_lo(ii-i+1)=s_lo(i)
             enddo

             ! 対角成分を持つことを仮定
             ! 上三角行列の対角成分は1のオフセット
             do j=1+1,nndU
                do i = ii ,iimin ,-1
                   if( ( ( 1 <= (i+offsetU(j)) ) ).and.( (i+offsetU(j) <= n ) ) )then
                      ! zz(ii-i+1)=zz(ii-i+1)-s(i+offsetU(j))*diaAU(i+(j-1)*n)
#ifdef Manual_inline
                      ddmul_calc1(s_hi(i+offsetU(j)),s_lo(i+offsetU(j)),diaAU_hi(i+(j-1)*n),0.0d0,c_hi,c_lo)
                      ddsub_calc1(zz_hi(ii-i+1),zz_lo(ii-i+1),c_hi,c_lo,zz_hi(ii-i+1),zz_lo(ii-i+1))
#else
                      call dd_mul( &
                           s_hi(i+offsetU(j)),s_lo(i+offsetU(j)) ,&
                           diaAU_hi(i+(j-1)*n), 0.0d0, &
                           c_hi,c_lo)
                      call dd_sub( &
                           zz_hi(ii-i+1),zz_lo(ii-i+1), &
                           c_hi,c_lo,&
                           zz_hi(ii-i+1),zz_lo(ii-i+1))
#endif
                   endif
                enddo
             enddo
             do i = ii ,iimin ,-1          
                j=1
                ! s(i)=zz(ii-i+1)/diaAU(i+(j-1)*n)
#ifdef Manual_inline
                dddiv_calc1(zz_hi(ii-i+1),zz_lo(ii-i+1),diaAU_hi(i+(j-1)*n),0.0d0,s_hi(i),s_lo(i))
#else
                call dd_div( &
                     zz_hi(ii-i+1),zz_lo(ii-i+1),&
                     diaAU_hi(i+(j-1)*n) ,0.0d0 , &
                     s_hi(i),s_lo(i))
#endif
             enddo
          enddo

#else
          ! -----
          ! z = Us 
          ! Lz = r をzについて解く
          ! 前進代入
          do i=1,n
             z=r(i)
             ! 対角成分を持つことを仮定
             ! 下三角行列で対角成分はnndLのオフセット
             do j=1,nndL-1
                if( (1 <= (i+offsetL(j))).and.((i+offsetL(j)) <= n ) )then
                   z=z-s(i+offsetL(j))*diaAL(i+(j-1)*n)
                endif
             enddo
             j=nndL
             s(i)=z/diaAL(i+(j-1)*n)
          enddo

          ! -----
          ! Us=z をsについて解く
          ! 後退代入
          do i=n,1,-1
             z=s(i)
             ! 対角成分を持つことを仮定
             ! 上三角行列の対角成分は1のオフセット
             do j=1+1,nndU
                if( ( ( 1 <= (i+offsetU(j)) ) ).and.( (i+offsetU(j) <= n ) ) )then
                   z=z-s(i+offsetU(j))*diaAU(i+(j-1)*n)
                endif
             enddo
             j=1
             s(i)=z/diaAU(i+(j-1)*n)
          enddo
#endif
       endif

    case(KSP_DOUBLE_DOUBLE_PRECISION)
       ! 零以外の対角成分を持つことを仮定
       ! -----
       ! LUs = r をsについて解いていく
       if(independ_nvec.eq.1)then
          ! -----
          ! z = Us 
          ! Lz = r をzについて解く
          ! 前進代入
          do i=1,n
             z_hi=r_hi(i)
             z_lo=r_lo(i)
!             write(*,*)"in  -- ",i,r_hi(i)+r_lo(i)
             ! 対角成分を持つことを仮定
             ! 下三角行列で対角成分はnndLのオフセット
             do j=1,nndL-1
                if( (1 <= (i+offsetL(j))).and.((i+offsetL(j)) <= n ) )then
                   !                z=z-s(i+offsetL(j))*diaAL(i+(j-1)*n)
#ifdef Manual_inline
                   ddmul_calc1(s_hi(i+offsetL(j)),s_lo(i+offsetL(j)),diaAL_hi(i+(j-1)*n),diaAL_lo(i+(j-1)*n),c_hi,c_lo)
                   ddsub_calc1(z_hi,z_lo,c_hi,c_lo,z_hi,z_lo)
#else
                   call dd_mul( &
                        s_hi(i+offsetL(j)),s_lo(i+offsetL(j)), &
                        diaAL_hi(i+(j-1)*n),diaAL_lo(i+(j-1)*n), &
                        c_hi,c_lo)
                   call dd_sub( &
                        z_hi,z_lo, &
                        c_hi,c_lo,&
                        z_hi,z_lo)
#endif
                endif
             enddo
             j=nndL
             !          s(i)=z/diaAL(i+(j-1)*n)
!             write(*,*)"save -- ",i,z_hi+z_lo,diaAL_hi(i+(j-1)*n)+diaAL_lo(i+(j-1)*n)
!             write(*,*)j,"save -- ",i,i+(j-1)*n,diaAL_hi(i+(j-1)*n)
#ifdef Manual_inline
             dddiv_calc1(z_hi,z_lo,diaAL_hi(i+(j-1)*n),diaAL_lo(i+(j-1)*n),s_hi(i),s_lo(i))
#else
             call dd_div( &
                  z_hi,z_lo,&
                  diaAL_hi(i+(j-1)*n) ,diaAL_lo(i+(j-1)*n) , &
                  s_hi(i),s_lo(i))
#endif
          enddo

          ! -----
          ! Us=z をsについて解く
          ! 後退代入
          do i=n,1,-1
             ! z=s(i)
             z_hi=s_hi(i)
             z_lo=s_lo(i)
!             write(*,*)"save2 -- ",i,s_hi(i)+s_lo(i)
             ! 対角成分を持つことを仮定
             ! 上三角行列の対角成分は1のオフセット
             do j=1+1,nndU
                if( ( ( 1 <= (i+offsetU(j)) ) ).and.( (i+offsetU(j) <= n ) ) )then
                   ! z=z-s(i+offsetU(j))*diaAU(i+(j-1)*n)
#ifdef Manual_inline
                   ddmul_calc1(s_hi(i+offsetU(j)),s_lo(i+offsetU(j)),diaAU_hi(i+(j-1)*n),diaAU_lo(i+(j-1)*n),c_hi,c_lo)
                   ddsub_calc1(z_hi,z_lo,c_hi,c_lo,z_hi,z_lo)
#else
                   call dd_mul( &
                        s_hi(i+offsetU(j)),s_lo(i+offsetU(j)) ,&
                        diaAU_hi(i+(j-1)*n),diaAU_lo(i+(j-1)*n), &
                        c_hi,c_lo)
                   call dd_sub( &
                        z_hi,z_lo, &
                        c_hi,c_lo,&
                        z_hi,z_lo)
#endif
                endif
             enddo
             j=1
             ! s(i)=z/diaAU(i+(j-1)*n)
!             write(*,*)"out  -- ",i,z_hi+z_lo,diaAU_hi(i+(j-1)*n)+diaAU_lo(i+(j-1)*n)

#ifdef Manual_inline
             dddiv_calc1(z_hi,z_lo,diaAU_hi(i+(j-1)*n),diaAU_lo(i+(j-1)*n),s_hi(i),s_lo(i))
#else
             call dd_div( &
                  z_hi,z_lo,&
                  diaAU_hi(i+(j-1)*n) ,diaAU_lo(i+(j-1)*n) , &
                  s_hi(i),s_lo(i))
#endif
          enddo

       else
#if 1
          ! -----
          ! z = Us 
          ! Lz = r をzについて解く
          ! 前進代入
          !       do i=1,n
          do ii = 1,n,independ_nvec
             iimin = min(ii+independ_nvec-1,n)
             do i = ii ,iimin
                ! zz(i-ii+1)=r(i)
                zz_hi(i-ii+1)=r_hi(i)
                zz_lo(i-ii+1)=r_lo(i)
             enddo
             ! 対角成分を持つことを仮定
             ! 下三角行列で対角成分はnndLのオフセット
             do j=1,nndL-1
                do i = ii ,iimin
                   if( (1 <= (i+offsetL(j))).and.((i+offsetL(j)) <= n ) )then
                      ! zz(i-ii+1)=zz(i-ii+1)-s(i+offsetL(j))*diaAL(i+(j-1)*n)
#ifdef Manual_inline
                      ddmul_calc1(s_hi(i+offsetL(j)),s_lo(i+offsetL(j)),diaAL_hi(i+(j-1)*n),diaAL_lo(i+(j-1)*n),c_hi,c_lo)
                      ddsub_calc1(zz_hi(i-ii+1),zz_lo(i-ii+1),c_hi,c_lo,zz_hi(i-ii+1),zz_lo(i-ii+1))
#else
                      call dd_mul( &
                           s_hi(i+offsetL(j)),s_lo(i+offsetL(j)), &
                           diaAL_hi(i+(j-1)*n),diaAL_lo(i+(j-1)*n), &
                           c_hi,c_lo)
                      call dd_sub( &
                           zz_hi(i-ii+1),zz_lo(i-ii+1), &
                           c_hi,c_lo,&
                           zz_hi(i-ii+1),zz_lo(i-ii+1))
#endif
                   endif
                enddo
             enddo
             do i = ii ,iimin
                j=nndL
                ! s(i)=zz(i-ii+1)/diaAL(i+(j-1)*n)
#ifdef Manual_inline
                dddiv_calc1(zz_hi(i-ii+1),zz_lo(i-ii+1),diaAL_hi(i+(j-1)*n),diaAL_lo(i+(j-1)*n),s_hi(i),s_lo(i))
#else
                call dd_div( &
                     zz_hi(i-ii+1),zz_lo(i-ii+1),&
                     diaAL_hi(i+(j-1)*n) ,diaAL_lo(i+(j-1)*n) , &
                     s_hi(i),s_lo(i))
#endif
             enddo
          enddo

          ! -----
          ! Us=z をsについて解く
          ! 後退代入
          do ii = n,1,-independ_nvec
             iimin = max(ii-independ_nvec+1,1) 
             do i = ii ,iimin ,-1
                ! zz(ii-i+1)=s(i)
                zz_hi(ii-i+1)=s_hi(i)
                zz_lo(ii-i+1)=s_lo(i)
             enddo

             ! 対角成分を持つことを仮定
             ! 上三角行列の対角成分は1のオフセット
             do j=1+1,nndU
                do i = ii ,iimin ,-1
                   if( ( ( 1 <= (i+offsetU(j)) ) ).and.( (i+offsetU(j) <= n ) ) )then
                      ! zz(ii-i+1)=zz(ii-i+1)-s(i+offsetU(j))*diaAU(i+(j-1)*n)
#ifdef Manual_inline
                      ddmul_calc1(s_hi(i+offsetU(j)),s_lo(i+offsetU(j)),diaAU_hi(i+(j-1)*n),diaAU_lo(i+(j-1)*n),c_hi,c_lo)
                      ddsub_calc1(zz_hi(ii-i+1),zz_lo(ii-i+1),c_hi,c_lo,zz_hi(ii-i+1),zz_lo(ii-i+1))
#else
                      call dd_mul( &
                           s_hi(i+offsetU(j)),s_lo(i+offsetU(j)) ,&
                           diaAU_hi(i+(j-1)*n),diaAU_lo(i+(j-1)*n), &
                           c_hi,c_lo)
                      call dd_sub( &
                           zz_hi(ii-i+1),zz_lo(ii-i+1), &
                           c_hi,c_lo,&
                           zz_hi(ii-i+1),zz_lo(ii-i+1))
#endif
                   endif
                enddo
             enddo
             do i = ii ,iimin ,-1          
                j=1
                ! s(i)=zz(ii-i+1)/diaAU(i+(j-1)*n)
#ifdef Manual_inline
                dddiv_calc1(zz_hi(ii-i+1),zz_lo(ii-i+1),diaAU_hi(i+(j-1)*n),diaAU_lo(i+(j-1)*n),s_hi(i),s_lo(i))
#else
                call dd_div( &
                     zz_hi(ii-i+1),zz_lo(ii-i+1),&
                     diaAU_hi(i+(j-1)*n) ,diaAU_lo(i+(j-1)*n) , &
                     s_hi(i),s_lo(i))
#endif
             enddo
          enddo

#else
          ! -----
          ! z = Us 
          ! Lz = r をzについて解く
          ! 前進代入
          do i=1,n
             z=r(i)
             ! 対角成分を持つことを仮定
             ! 下三角行列で対角成分はnndLのオフセット
             do j=1,nndL-1
                if( (1 <= (i+offsetL(j))).and.((i+offsetL(j)) <= n ) )then
                   z=z-s(i+offsetL(j))*diaAL(i+(j-1)*n)
                endif
             enddo
             j=nndL
             s(i)=z/diaAL(i+(j-1)*n)
          enddo

          ! -----
          ! Us=z をsについて解く
          ! 後退代入
          do i=n,1,-1
             z=s(i)
             ! 対角成分を持つことを仮定
             ! 上三角行列の対角成分は1のオフセット
             do j=1+1,nndU
                if( ( ( 1 <= (i+offsetU(j)) ) ).and.( (i+offsetU(j) <= n ) ) )then
                   z=z-s(i+offsetU(j))*diaAU(i+(j-1)*n)
                endif
             enddo
             j=1
             s(i)=z/diaAU(i+(j-1)*n)
          enddo
#endif
       endif

    end select

  end subroutine ddsolve_LU_dia

  subroutine ddblockJacobi_ILU_crs(method,r_hi,r_lo,s_hi,s_lo)
    use omp_lib
    use krylov
    use PARCEL_integer

    implicit none
    type(KSP) :: method
    real*8 r_hi(method%prm%n)
    real*8 r_lo(method%prm%n)
    real*8 s_hi(method%prm%n)
    real*8 s_lo(method%prm%n)

    integer(kind=intType) n
    integer(kind=intType) nnz
    integer(kind=intType) istart,iend

    integer myth,numth,max_numth
    integer(kind=intType) n_th
    integer(kind=intType) n_tmp
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th

    n=method%prm%n
    nnz=method%prm%nnz
    istart=method%prm%istart
    iend=method%prm%iend

    ! -----
    ! ブロック数よりもスレッド数が大きい場合はスレッド数を制限する
    call get_maxthreads(max_numth)
!    call set_precon_maxthreads0(n)
    call set_precon_maxthreads(n,method%precon_thblock)

    !$OMP PARALLEL &
    !$OMP PRIVATE( &
    !$OMP numth,myth,&
    !$OMP n_startth,n_endth,&
    !$OMP istart_th,iend_th,n_tmp,n_th &
    !$OMP )
    numth = omp_get_num_threads()
    myth  = omp_get_thread_num()

    n_tmp=1
    call divide_region( &
         myth,numth,n_tmp,n,n_startth,n_endth)
! ---
!     call view_crs(&
!          method%preconILU(myth+1)%n_th,&
!          method%preconILU(myth+1)%istart_th,&
!          method%preconILU(myth+1)%crsMatU )
! ---
    call ddsolve_LU_crs( &
         method%preconILU(myth+1)%n_th, &
         method%preconILU(myth+1)%istart_th, &
         method%preconILU(myth+1)%nnzL_th, &
         method%preconILU(myth+1)%crsMatL%crsA_hi, &
         method%preconILU(myth+1)%crsMatL%crsA_lo, &
         method%preconILU(myth+1)%crsMatL%crsRow_ptr, &
         method%preconILU(myth+1)%crsMatL%crsCol, &
         method%preconILU(myth+1)%nnzU_th, &
         method%preconILU(myth+1)%crsMatU%crsA_hi, &
         method%preconILU(myth+1)%crsMatU%crsA_lo, &
         method%preconILU(myth+1)%crsMatU%crsRow_ptr, &
         method%preconILU(myth+1)%crsMatU%crsCol, &
         r_hi(n_startth), &
         r_lo(n_startth), &
         s_hi(n_startth), &
         s_lo(n_startth), &
         method%independ_nvec, &
         method%PRECISION_PRECONDITION &
       )

    !$OMP END PARALLEL

!    if(n<=max_numth)then
    call omp_set_num_threads(max_numth)
!    endif
  end subroutine ddblockJacobi_ILU_crs

  subroutine ddblockJacobi_ILU_dia(method,r_hi,r_lo,s_hi,s_lo)

    use omp_lib
    use krylov
    use PARCEL_integer

    implicit none
    type(KSP) :: method
    real*8 r_hi(method%prm%n)
    real*8 r_lo(method%prm%n)
    real*8 s_hi(method%prm%n)
    real*8 s_lo(method%prm%n)

    integer(kind=intType) n
    integer(kind=intType) nnz
    integer(kind=intType) istart,iend

    integer precon_thblock
    integer myth,numth,max_numth
    integer(kind=intType) n_th
    integer(kind=intType) n_tmp
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th

    n=method%prm%n
    nnz=method%prm%nnz
    istart=method%prm%istart
    iend=method%prm%iend
    precon_thblock=method%precon_thblock 

    ! -----
    ! ブロック数よりもスレッド数が大きい場合はスレッド数を制限する
    call get_maxthreads(max_numth)
    call set_precon_maxthreads(n,precon_thblock)

    !$OMP PARALLEL &
    !$OMP PRIVATE( &
    !$OMP numth,myth,&
    !$OMP n_startth,n_endth,&
    !$OMP istart_th,iend_th,n_tmp,n_th &
    !$OMP )
    numth = omp_get_num_threads()
    myth  = omp_get_thread_num()

    n_tmp=1
    call divide_region( &
         myth,numth,n_tmp,n,n_startth,n_endth)
    call ddsolve_LU_dia( &
         method%preconILU(myth+1)%n_th, &
         method%preconILU(myth+1)%istart_th, &
         method%preconILU(myth+1)%nndL_th, &
         method%preconILU(myth+1)%diaMatL%diaA_hi, &
         method%preconILU(myth+1)%diaMatL%diaA_lo, &
         method%preconILU(myth+1)%diaMatL%offset, &
         method%preconILU(myth+1)%nndU_th, &
         method%preconILU(myth+1)%diaMatU%diaA_hi, &
         method%preconILU(myth+1)%diaMatU%diaA_lo, &
         method%preconILU(myth+1)%diaMatU%offset, &
         r_hi(n_startth), &
         r_lo(n_startth), &
         s_hi(n_startth), &
         s_lo(n_startth), &
         method%independ_nvec, &
         method%PRECISION_PRECONDITION &
         )

    !$OMP END PARALLEL

!    if(n<=max_numth)then
    call omp_set_num_threads(max_numth)
!    endif


  end subroutine ddblockJacobi_ILU_dia

  subroutine ddfree_precon_ILU_crs(method)
    use krylov
    use PARCEL_integer
    implicit none
    type(KSP) :: method
    integer max_numth
    integer(kind=intType) i
    max_numth=method%max_numth

    if(allocated(method%preconILU))then
       do i=1,max_numth
          call ddfree_crsMat(method%preconILU(i)%crsMatL)
          call ddfree_crsMat(method%preconILU(i)%crsMatU)
       enddo
       deallocate(method%preconILU)
    endif
  end subroutine ddfree_precon_ILU_crs

  subroutine ddfree_precon_ILU_dia(method)
    use krylov
    use PARCEL_integer
    implicit none
    type(KSP) :: method
    integer max_numth
    integer(kind=intType) i
    max_numth=method%max_numth

    if(allocated(method%preconILU))then
       do i=1,max_numth
          call ddfree_diaMat(method%preconILU(i)%diaMatL)
          call ddfree_diaMat(method%preconILU(i)%diaMatU)
       enddo
       deallocate(method%preconILU)
    endif
  end subroutine ddfree_precon_ILU_dia


! ---

end module ddmatrix_kernel

