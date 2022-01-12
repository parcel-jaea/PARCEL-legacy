! todo dia,craのコードを整理
module matrix_kernel
  use PARCEL_define
  use common_kernel

  implicit none
contains

  subroutine init_denMat(denMat,n,gn)
    use matrix
    use PARCEL_integer
    implicit none
    type(DEN_Mat) :: denMat 
    integer(kind=intType) n  ! ブロック行数 
    integer(kind=intType) gn ! 列数
    integer(kind=intType) i
    denMat%RowN=n
    denMat%ColN=gn
    allocate(denMat%denA(n*gn))
    !$OMP PARALLEL DO PRIVATE(i)
    do i=1,n*gn
       denMat%denA(i)=0.0d0
    enddo
    !$OMP END PARALLEL DO

  end subroutine init_denMat

  subroutine init_cooMat(cooMat,nnz)
    use matrix
    use PARCEL_integer
    implicit none
    type(COO_Mat) :: cooMat 
    integer(kind=intType) nnz
    integer(kind=intType) i
    allocate(cooMat%cooA(nnz))
    allocate(cooMat%cooRow(nnz))
    allocate(cooMat%cooCol(nnz))

    !$OMP PARALLEL DO PRIVATE(i)
    do i=1,nnz
       cooMat%cooA(i)=0.0d0
       cooMat%cooRow(i)=0
       cooMat%cooCol(i)=0
    enddo
    !$OMP END PARALLEL DO

  end subroutine init_cooMat

  subroutine init_crsMat(crsMat,n,nnz)
    use matrix
    use PARCEL_integer
    implicit none
    type(CRS_Mat) :: crsMat 
    integer(kind=intType) n
    integer(kind=intType) nnz
    integer(kind=intType) i

    allocate( crsMat%crsA(nnz) )
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
       crsMat%crsA(i)=0.0d0
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

  end subroutine init_crsMat

  subroutine init_diaMat(diaMat,n,nnd)
    use matrix
    use PARCEL_integer
    implicit none
    type(DIA_Mat) :: diaMat 
    integer(kind=intType) n
    integer(kind=intType) nnd
    integer(kind=intType) i
    allocate( diaMat%diaA(n*nnd) )
    allocate( diaMat%offset(nnd) )
    allocate( diaMat%offsetVec(nnd) )
    allocate( diaMat%Vec(1) )

    allocate(diaMat%RowL(1))
    allocate(diaMat%RowD(1))
    allocate(diaMat%RowU(1))

    !$OMP PARALLEL DO PRIVATE(i)
    do i=1,n*nnd
       diaMat%diaA(i)=0.0d0
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i)
    do i=1,nnd
       diaMat%offset(i)=0
       diaMat%offsetVec(i)=0
    enddo
    !$OMP END PARALLEL DO

  end subroutine init_diaMat

  subroutine get_crs_Row( &
       prm,crsMat,nCol,Colsort &
       )
    use matrix
    use PARCEL_integer
    implicit none
    type(prmMat) :: prm
    type(CRS_Mat) :: crsMat
    integer(kind=intType) n,nCol
    integer(kind=intType) Colsort(nCol)
    integer(kind=intType) istart,iend

    integer(kind=intType) nColL,nColU
    n=prm%n
    istart=prm%istart
    iend=prm%iend

    nColL=crsMat%nColL
    nColU=crsMat%nColU
    deallocate(crsMat%RowL)
    deallocate(crsMat%RowD)
    deallocate(crsMat%RowU)

    allocate(crsMat%RowL(nColL))
    allocate(crsMat%RowD(n))
    allocate(crsMat%RowU(nColU))
    call get_RowLUD( &
         n,nCol,nColL,nColU,&
         Colsort, &
         crsMat%RowL, &
         crsMat%RowU, &
         crsMat%RowD, &
         istart,iend &
         )
  end subroutine get_crs_Row

  subroutine get_crsCol_nCol( &
       prm, &
       crsCol,Colsort, &
       nCol,nColL,nColU,nColD)
    use sort_kernel
    use matrix
    use PARCEL_integer
    implicit none
    type(prmMat) :: prm
    integer(kind=intType) istart,iend
    integer(kind=intType) nnz
    integer(kind=intType) crsCol(prm%nnz)
    integer(kind=intType) Colsort(prm%nnz)
    integer(kind=intType) nCol
    integer(kind=intType) nColL,nColU,nColD
    integer(kind=intType) i
    integer(kind=intType) Col
    nnz=prm%nnz
    istart=prm%istart
    iend=prm%iend

    !$OMP PARALLEL DO PRIVATE(i)
    do i=1,nnz
       Colsort(i)=crsCol(i)
    enddo
    !$OMP END PARALLEL DO
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
    ! do i=nCol+1,nnz
    !    Colsort(i)=-1
    ! enddo
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
    !    write(*,*)"nColU,nColD,nColL = ",nColU,nColD,nColL
  end subroutine get_crsCol_nCol

  subroutine init_HybcrsMat(HybcrsMat,n)
    use matrix
    use PARCEL_integer
    implicit none
    type(HybCRS_Mat) :: HybcrsMat 
    integer(kind=intType) n

    integer(kind=intType) i
    integer(kind=intType) nnzL,nnzD,nnzU
    integer(kind=intType) nColL,nColU,nColD

    nnzL  = HybcrsMat%nnzL
    nnzD  = HybcrsMat%nnzD
    nnzU  = HybcrsMat%nnzU
    nColL = HybcrsMat%nColL
    nColD = HybcrsMat%nColD
    nColU = HybcrsMat%nColU

    allocate( HybcrsMat%crsAL(nnzL) )
    allocate( HybcrsMat%crsRow_ptr_L(n+1) )
    allocate( HybcrsMat%crsColL(nnzL) )

    allocate( HybcrsMat%crsAD(nnzD) )
    allocate( HybcrsMat%crsRow_ptr_D(n+1) )
    allocate( HybcrsMat%crsColD(nnzD) )

    allocate( HybcrsMat%crsAU(nnzU) )
    allocate( HybcrsMat%crsRow_ptr_U(n+1) )
    allocate( HybcrsMat%crsColU(nnzU) )

    allocate( HybcrsMat%RowL(nColL) )
    allocate( HybcrsMat%RowD(n) )
    allocate( HybcrsMat%RowU(nColU) )

    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,nnzL
       HybcrsMat%crsAL(i)=0.0d0
       HybcrsMat%crsColL(i)=0
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,nnzD
       HybcrsMat%crsAD(i)=0.0d0
       HybcrsMat%crsColD(i)=0
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,nnzU
       HybcrsMat%crsAU(i)=0.0d0
       HybcrsMat%crsColU(i)=0
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n+1
       HybcrsMat%crsRow_ptr_L(i)=0
       HybcrsMat%crsRow_ptr_D(i)=0
       HybcrsMat%crsRow_ptr_U(i)=0
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,nColL
       HybcrsMat%RowL(i)=0
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
       HybcrsMat%RowD(i)=0
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,nColU
       HybcrsMat%RowU(i)=0
    enddo
    !$OMP END PARALLEL DO

  end subroutine init_HybcrsMat

  subroutine get_nnz_den(denMat,nnz)
    use matrix
    use PARCEL_integer
    implicit none
    type(DEN_Mat) :: denMat 
    integer(kind=intType) nnz
    real*8 absAmat
    integer(kind=intType) gn,n
    integer(kind=intType) i,j
    n =denMat%RowN
    gn=denMat%ColN

    nnz=0
    !$OMP PARALLEL DO PRIVATE(i,j,absAmat) REDUCTION(+:nnz)
    do i=1,n
       do j=1,gn
#if 1
          absAmat=abs(denMat%denA(j+(i-1)*gn))
          if(absAmat> minReal8)then
             nnz=nnz+1
             !             write(*,*)i,j,denMat%denA(j+(i-1)*gn)
          endif
#else
          nnz=nnz+1
#endif
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine get_nnz_den

  subroutine get_nnz_crs(crsMat,n,nnz)
    use matrix
    use PARCEL_integer
    implicit none
    type(CRS_Mat) :: crsMat 
    integer(kind=intType) n
    integer(kind=intType) nnz
    integer(kind=intType) i,j,k
    nnz=0
    ! nnzで操作しなくてもcrsRow_ptrの最終値 - 1 で計算可能
    !$OMP PARALLEL DO PRIVATE(k) REDUCTION(+:nnz)
    do k=1,n
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          nnz=nnz+1
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine get_nnz_crs

  subroutine get_nnzLDU_crs(crsMat,istart,iend,n,nnzL,nnzD,nnzU)
    use matrix
    use PARCEL_integer
    implicit none
    type(CRS_Mat) :: crsMat 
    integer(kind=intType) istart,iend
    integer(kind=intType) n
    integer(kind=intType) nnzL,nnzD,nnzU
    integer(kind=intType) Row,Col
    integer(kind=intType) i,j,k
    nnzL=0
    nnzD=0
    nnzU=0

    !$OMP PARALLEL DO PRIVATE(k,j,Row,Col) REDUCTION(+:nnzL,nnzD,nnzU)
    do k=1,n
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          Row=k-1+istart 
          Col=crsMat%crsCol(j)  
          if( (istart<=Row) .and. (Row<=iend) )then
             if( (istart<=Col) )then
                if(  (Col<=iend) )then
                   nnzD=nnzD+1
                else
                   nnzU=nnzU+1
                endif
             else
                nnzL=nnzL+1
             endif
          endif
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine get_nnzLDU_crs

  subroutine MatVec_Hybcrs(HybcrsMat,n,xL,xD,xU,Ax)
    use matrix
    use PARCEL_integer
    implicit none
    type(HYBCRS_Mat) :: HybcrsMat 
    integer(kind=intType) n
    real*8 xL(HybcrsMat%nColL)
    real*8 xD(n)
    real*8 xU(HybcrsMat%nColU)

    real*8 Ax(n)
    integer(kind=intType) i,j,k
    integer(kind=intType) istart

    ! do k=1,n
    !    write(*,*)" vec ",k,xD(k)
    ! enddo

    !    !$OMP PARALLEL DO PRIVATE(k,j) 
    do k=1,n
       Ax(k) = 0.0d0
       do j=HybcrsMat%crsRow_ptr_L(k),HybcrsMat%crsRow_ptr_L(k+1)-1
          !          write(*,*)"MatVec L ",HybcrsMat%crsAL(j),xL(HybcrsMat%crsColL(j))
          Ax(k) = Ax(k) + HybcrsMat%crsAL(j)*xL(HybcrsMat%crsColL(j))
       enddo
       do j=HybcrsMat%crsRow_ptr_D(k),HybcrsMat%crsRow_ptr_D(k+1)-1
          !          write(*,*)"MatVec D ",HybcrsMat%crsAD(j) , xD(HybcrsMat%crsColD(j)),HybcrsMat%crsColD(j)
          Ax(k) = Ax(k) + HybcrsMat%crsAD(j) * xD(HybcrsMat%crsColD(j))
       enddo
       do j=HybcrsMat%crsRow_ptr_U(k),HybcrsMat%crsRow_ptr_U(k+1)-1
          !          write(*,*)"MatVec U ", HybcrsMat%crsAU(j) , xU(HybcrsMat%crsColU(j))
          Ax(k) = Ax(k) + HybcrsMat%crsAU(j) * xU(HybcrsMat%crsColU(j))
       enddo
    enddo
    !    !$OMP END PARALLEL DO

  end subroutine MatVec_Hybcrs
! #define nBlock 2000
  subroutine MatVec_crs(prm,crsMat,xL,xD,xU,Ax,nBlock)
    use matrix
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(CRS_Mat) :: crsMat 
    real*8 xL(crsMat%nColL)
    real*8 xD(prm%n)
    real*8 xU(crsMat%nColU)

    real*8 tmp
    real*8 Ax(prm%n)
    real*8 Ax0(nBlock)
    integer(kind=intType) iptr_s,iptr_e
!    integer(kind=intType) crsRow_ptr_Diag_s()
!    integer(kind=intType) crsRow_ptr_Diag_e()

    integer(kind=intType) i,j,jj,k,kk,kkstep,kk0,kkmin
    integer(kind=intType) iLs,iUs,iLe,iUe
    integer(kind=intType) n,istart,iend
    integer(kind=intType) Row,Col
    integer(kind=intType) nColL,nColU
    integer nBlock
    
    n = prm%n
    istart=prm%istart
    iend=prm%iend
    nColL=crsMat%nColL
    nColU=crsMat%nColU

!    call fapp_start("MatVecCRS",1,0)
#if 1
    !$OMP PARALLEL DO PRIVATE(i,j,k,iLs,iUs,iLe,iUe,Row,Col,tmp) SCHEDULE(static,nBlock)
    do k=1,n
       Ax(k)=0.0d0
       tmp=0.0d0
       Row=k-1+istart
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          if( crsMat%crsCol(j)<istart ) then
             Col=crsMat%iLs(j)
             !Ax(k)=Ax(k)+crsMat%crsA(j)*xL(Col)
             tmp=tmp+crsMat%crsA(j)*xL(Col)
             cycle
          endif
          if( iend<crsMat%crsCol(j) )then
             Col=crsMat%iLs(j)
             !Ax(k)=Ax(k)+crsMat%crsA(j)*xU(Col)
             tmp=tmp+crsMat%crsA(j)*xU(Col)
             cycle
          endif
          !          Col=crsMat%crsCol(j)-istart+1
          Col=crsMat%iLs(j)
          !Ax(k)=Ax(k)+crsMat%crsA(j)*xD(Col)
          tmp=tmp+crsMat%crsA(j)*xD(Col)
       enddo
       Ax(k)=tmp
    enddo
    !$OMP END PARALLEL DO
#else

    !$OMP PARALLEL DO PRIVATE(i,j,jj,k,kk,kkstep,kk0,kkmin,Ax0,iLs,iUs,iLe,iUe,Row,Col, iptr_s,iptr_e,tmp)
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

#if 1
          do k = kk0,kk0+kkstep
             tmp=0.0d0
             do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
                if( crsMat%crsCol(j)<istart ) then
                   Col=crsMat%iLs(j)
                   tmp=tmp+crsMat%crsA(j)*xL(Col)
                   cycle
                endif
                exit
             enddo              
             do j=crsMat%crsRow_ptr(k+1)-1,crsMat%crsRow_ptr(k),-1
                if( iend<crsMat%crsCol(j) )then
                   Col=crsMat%iLs(j)
                   tmp=tmp+crsMat%crsA(j)*xU(Col)
                   cycle
                endif
                exit
             enddo
             Ax0(k-kk+1) = tmp
          enddo
          
          ! if( 0 < (crsMat%crsRow_ptr_Diag_e(kk0) -  crsMat%crsRow_ptr_Diag_s(kk0)) )then
          !    do jj=1,crsMat%crsRow_ptr_Diag_e(kk0) -  crsMat%crsRow_ptr_Diag_s(kk0) 
          !       do k = kk0,kk0+kkstep
          !          j   = crsMat%crsRow_ptr_Diag_s(k) + jj -1
          !          Col=crsMat%iLs(j)
          !          Ax0(k-kk+1)=Ax0(k-kk+1)+crsMat%crsA(j)*xD(Col)
          !       enddo
          !    enddo
          !    jj = crsMat%crsRow_ptr_Diag_e(kk0) -  crsMat%crsRow_ptr_Diag_s(kk0) + 1             
          !    do k = kk0,kk0+kkstep
          !       j   = crsMat%crsRow_ptr_Diag_s(k) + jj -1
          !       Col=crsMat%iLs(j)
          !       Ax(k) = Ax0(k-kk+1)+crsMat%crsA(j)*xD(Col)
          !    enddo
          ! else
             do jj=1,crsMat%crsRow_ptr_Diag_e(kk0) -  crsMat%crsRow_ptr_Diag_s(kk0) + 1
                do k = kk0,kk0+kkstep
                   j   = crsMat%crsRow_ptr_Diag_s(k) + jj -1
                   Col=crsMat%iLs(j)
                   Ax0(k-kk+1)=Ax0(k-kk+1)+crsMat%crsA(j)*xD(Col)
                enddo
             enddo
             do k = kk0,kk0+kkstep
                Ax(k) = Ax0(k-kk+1)
             enddo
!          endif
#else
          do k = kk0,kk0+kkstep
             Ax0(k-kk+1)=0.0d0
             do j=crsMat%crsRow_ptr(k+1)-1,crsMat%crsRow_ptr_Diag_e(k)+1,-1
                Col=crsMat%iLs(j)
                Ax0(k-kk+1)=Ax0(k-kk+1)+crsMat%crsA(j)*xU(Col)
             enddo
             do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr_Diag_s(k)-1                
                Col=crsMat%iLs(j)
                Ax0(k-kk+1)=Ax0(k-kk+1)+crsMat%crsA(j)*xL(Col)
             enddo
          enddo
          do jj=1,crsMat%crsRow_ptr_Diag_e(kk0) -  crsMat%crsRow_ptr_Diag_s(kk0) + 1
             do k = kk0,kk0+kkstep
                j   = crsMat%crsRow_ptr_Diag_s(k) + jj -1
                Col=crsMat%iLs(j)
                Ax0(k-kk+1)=Ax0(k-kk+1)+crsMat%crsA(j)*xD(Col)
             enddo             
          enddo
          do k = kk0,kk0+kkstep
             Ax(k) = Ax0(k-kk+1)
          enddo
          
#endif
          
          
          kk0 = kk0+kkstep+1
          if(kk0>kkmin)then
             exit
          endif
       enddo

    enddo
    !$OMP END PARALLEL DO
    
#endif

!    call fapp_stop("MatVecCRS",1,0)
    
  end subroutine MatVec_crs

  subroutine init_MatVec_crs(prm,crsMat,nBlock)
    use matrix
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(CRS_Mat) :: crsMat 

    integer(kind=intType) i,j,k,iL,iU
    integer(kind=intType) n,nnz,istart,iend
    integer(kind=intType) Row,Col
    integer(kind=intType) nColL,nColU
    integer nBlock

    n = prm%n
    nnz = prm%nnz
    istart=prm%istart
    iend=prm%iend
    nColL=crsMat%nColL
    nColU=crsMat%nColU

    !$OMP PARALLEL DO PRIVATE(i,j,k,iL,iU,Row,Col) SCHEDULE(dynamic,nBlock) 
    do k=1,n
       Row=k-1+istart
       iL=1
       iU=1
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          Col=crsMat%crsCol(j)
          crsMat%iLs(j)=crsMat%crsCol(j)-istart+1
          if( Col<istart ) then
#if 1
             i = 1
             call binary_search(crsMat%RowL,Col,i,nColL,crsMat%iLs(j))
#else
             do i=iL,nColL
                if( crsMat%RowL(i).eq.Col )then
                   crsMat%iLs(j)=i
                   iL=i
                   exit
                endif
             enddo
#endif
             !             cycle
          endif

          if( iend<Col )then
#if 1
             i = 1;
             call binary_search(crsMat%RowU,Col,i,nColU,crsMat%iLs(j))
#else
             do i=iU,nColU
                if( crsMat%RowU(i).eq.Col )then
                   crsMat%iLs(j)=i
                   iU=i
                   exit
                endif
             enddo
#endif
             !             cycle
          endif

       enddo

    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k,Row,Col) SCHEDULE(dynamic,nBlock) 
    do k=1,n
!       write(*,*)" k = ",k
       crsMat%crsRow_ptr_Diag_s(k) = -1
       crsMat%crsRow_ptr_Diag_e(k) = -2
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          if( (istart<=crsMat%crsCol(j)).and.(crsMat%crsCol(j)<=iend)  ) then
             crsMat%crsRow_ptr_Diag_s(k) = j
             exit
          endif
       enddo
       do j=crsMat%crsRow_ptr(k+1)-1,crsMat%crsRow_ptr(k),-1
          if( (istart<=crsMat%crsCol(j)).and.(crsMat%crsCol(j)<=iend)  ) then
              crsMat%crsRow_ptr_Diag_e(k) = j
             exit
          endif
       enddo       
    enddo
    !$OMP END PARALLEL DO
    
  end subroutine init_MatVec_crs

  subroutine init_MatVec_crs_org(prm,crsMat)
    use matrix
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(CRS_Mat) :: crsMat 

    integer(kind=intType) i,j,k,iL,iU
    integer(kind=intType) n,nnz,istart,iend
    integer(kind=intType) Row,Col
    integer(kind=intType) nColL,nColU
    n = prm%n
    nnz = prm%nnz
    istart=prm%istart
    iend=prm%iend
    nColL=crsMat%nColL
    nColU=crsMat%nColU

#if 1
    !$OMP PARALLEL DO PRIVATE(i,j,k,iL,iU,Row,Col) SCHEDULE(static) 
    do k=1,n
       Row=k-1+istart
       iL=1
       iU=1
       !       write(*,*)k,crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          crsMat%iLs(j)=0
          Col=crsMat%crsCol(j)
#if 1
          if( Col<istart ) then
             do i=iL,nColL
                if( crsMat%RowL(i).eq.Col )then
                   crsMat%iLs(j)=i
                   exit
                endif
             enddo
             cycle
          endif
          if( iend<Col )then
             do i=iU,nColU
                if( crsMat%RowU(i).eq.Col )then
                   crsMat%iLs(j)=i
                   exit
                endif
             enddo
             cycle
          endif
#else
          if( Col<istart ) then
             do i=iL,nColL
                if( crsMat%RowL(i).eq.Col )then
                   crsMat%iLs(j)=i
                   exit
                endif
             enddo
             cycle
          endif
          if( iend<Col )then
             do i=iU,nColU
                if( crsMat%RowU(i).eq.Col )then
                   !                   crsMat%iLs(j)=i
                   crsMat%iLs(j)=i+crsMat%nColL
                   exit
                endif
             enddo
             cycle
          endif
#endif
       enddo
    enddo
    !$OMP END PARALLEL DO
#else
    !$OMP PARALLEL DO PRIVATE(i,j,k,iL,iU,Row,Col) SCHEDULE(static) 
    !     do k=1,n
    !        Row=k-1+istart
    !        iL=1
    !        iU=1
    !       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
    do j=1,nnz
       crsMat%iLs(j)=0
       Col=crsMat%crsCol(j)
       if( Col<istart ) then
          do i=1,nColL
             if( crsMat%RowL(i).eq.Col )then
                crsMat%iLs(j)=i
                exit
             endif
          enddo
          cycle
       endif
       if( iend<Col )then
          do i=1,nColU
             if( crsMat%RowU(i).eq.Col )then
                crsMat%iLs(j)=i
                exit
             endif
          enddo
          cycle
       endif
    enddo
    !    enddo
    !$OMP END PARALLEL DO

#endif
  end subroutine init_MatVec_crs_org

  subroutine MatVec_crs_ovl1(prm,crsMat,xL,xD,xU,Ax,commVec,nBlock)
    use matrix
    use comm
    use comm_kernel
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(CRS_Mat) :: crsMat 
    type(comm_Vec) :: commVec
    real*8 xL(crsMat%nColL)
    real*8 xD(prm%n)
    real*8 xU(crsMat%nColU)

    real*8 Ax(prm%n)
    integer(kind=intType) i,j,k,iL,iU
    integer(kind=intType) n,istart,iend
    integer(kind=intType) Row,Col
    integer(kind=intType) nColL,nColU
    real*8 sbuf(commVec%nrk_data_send)
    real*8 rbuf(commVec%nrk_data)
    integer nBlock

    n = prm%n
    istart=prm%istart
    iend=prm%iend
    nColL=crsMat%nColL
    nColU=crsMat%nColU

    !     if(commVec%myrank.eq.0)then
    !        write(*,*)"@ crs ovl 1 @"
    !     endif
    call nrk_vec_pack( &
         commVec%nrk_prm_send, &
         commVec%n, &
         commVec%nrk_data_send, &
         xD, &
         sbuf, &
         commVec%nrk_rank_send)

    !$OMP PARALLEL 

    !$OMP MASTER
    call nrk_vec_comm(commVec,sbuf,rbuf)
    !$OMP END MASTER

    !$OMP DO PRIVATE(j,k,Col) SCHEDULE(dynamic,nBlock)
    do k=1,n
       Ax(k)=0.0d0
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          if(crsMat%crsCol(j)<istart)then
             cycle
          endif
          if( iend<crsMat%crsCol(j) )then
             exit
          endif
          Col=crsMat%crsCol(j)-istart+1
          Ax(k)=Ax(k)+crsMat%crsA(j)*xD(Col)
       enddo
    enddo
    !$OMP END PARALLEL

    call get_nrk_vecLU(&
         commVec, &
         rbuf, &
         nColL, &
         nColU, &
         xL, &
         xU &
         )

    !$OMP PARALLEL DO PRIVATE(i,j,k,iL,iU,Row,Col) SCHEDULE(static,nBlock)
    do k=1,n
       Row=k-1+istart

#if 0
       iL=1
       iU=1
#endif
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          if( crsMat%crsCol(j)<istart ) then
#if 0
             do i=iL,nColL
                if( crsMat%RowL(i).eq.crsMat%crsCol(j) )then
                   Col=i
                   Ax(k)=Ax(k)+crsMat%crsA(j)*xL(Col)
                   iL=i
                   exit
                endif
             enddo
             cycle
#else
             Col=crsMat%iLs(j)
             Ax(k)=Ax(k)+crsMat%crsA(j)*xL(Col)
#endif
          endif
          if( iend<crsMat%crsCol(j) )then
#if 0
             do i=iU,nColU
                if( crsMat%RowU(i).eq.crsMat%crsCol(j) )then
                   Col=i
                   Ax(k)=Ax(k)+crsMat%crsA(j)*xU(Col)
                   iU=i
                   exit
                endif
             enddo
             cycle
#else
             Col=crsMat%iLs(j)
             Ax(k)=Ax(k)+crsMat%crsA(j)*xU(Col)
             cycle
#endif
          endif
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine MatVec_crs_ovl1

  subroutine MatVec_crs_ovl2(prm,crsMat,x,Ax,commVec,nBlock)
    use matrix
    use comm
    use comm_kernel

    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(CRS_Mat) :: crsMat 
    type(comm_Vec) :: commVec
    real*8 x(prm%n)
    integer nBlock

    real*8 Ax(prm%n)
    integer(kind=intType) i,j,k,kk
    integer(kind=intType) n,istart,iend
    integer(kind=intType) Row,Col
    integer(kind=intType) nColL,nColU

    real*8 sbuf(commVec%nrk_data_max_send)
    real*8 rbuf(commVec%nrk_data_max)
    real*8 rbuf0(commVec%nrk_data_max)

    integer irank,irankMax
    integer srank,srankMax
    integer rrank,rrankMax
    integer istart_rrank,iend_rrank  
    integer rrank0

    integer recv_rank0,send_rank0

    integer s_req
    integer s_sta(MPI_STATUS_SIZE)
    integer r_req
    integer r_sta(MPI_STATUS_SIZE)
    integer ierr
    integer isflag,irflag
    LOGICAL isendflag,irecvflag

    integer rrank_count,srank_count
    integer rrank_calc

    srankMax=commVec%nrk_rank_send
    rrankMax=commVec%nrk_rank

    n = prm%n
    istart=prm%istart
    iend=prm%iend
    !     if(commVec%myrank.eq.0)then
    !        write(*,*)"@ crs ovl 2 @"
    !     endif
    ! -----
    send_rank0=1
    do irank=1,commVec%nrk_rank_send
       if(commVec%myrank < commVec%nrk_prm_send(irank)%rank)then
          send_rank0=irank
          exit
       endif
    enddo

    recv_rank0=commVec%nrk_rank
    do irank=commVec%nrk_rank,1,-1
       if(commVec%nrk_prm(irank)%rank < commVec%myrank   )then
          exit
       endif
       recv_rank0=irank
    enddo
    !    write(*,*)"send_rank0,recv_rank0",send_rank0,recv_rank0

    rrank_calc=-1
    if(commVec%nrk_rank<1)then
       rrank =-1
       irflag=-1
       rrank_count=rrankMax+1
    else
       irflag=0
       rrank_count=1
#if 0
       rrank=1
#else
       rrank=recv_rank0-rrank_count+1
       if(rrank<=0)then
          rrank=rrankMax+rrank
       endif
#endif
    endif

    if(commVec%nrk_rank_send<1)then
       srank =-1
       isflag=-1
       srank_count=srankMax+1
    else
       isflag=0
       srank_count=1
#if 0
       srank=1
#else
       srank=mod(send_rank0+srank_count+1,srankMax)+1
#endif

    endif
    rrank0=0 ! 0は自分の領域
    ! -----
    ! call view_nrk(commVec)

    ! isflag ( 全送信終了:-1,  送信完了:0, 送信中:1)
    ! irflag ( 全受信終了:-1,  受信完了:0, 受信中:1)
    do 

#if 0
       if( (srank.ne.-1).and.(rrank.ne.-1) )then
          write(*,*)&
               commVec%myrank, &
               commVec%nrk_prm(rrank)%rank, &
               commVec%nrk_prm_send(srank)%rank
       endif
       call system("sleep 1")

#endif
       if( (0<srank).and.(isflag.eq.0) )then
          call nrk_vec_pack_ovl2( &
               commVec%nrk_prm_send, &
               commVec%n, &
               commVec%nrk_data_send, &
               x, &
               sbuf, &
               commVec%nrk_rank_send,&
               srank)
          call nrk_vec_comm_ovl2_isend( &
               commVec,sbuf,srank,s_req,s_sta)
          isflag=1
       endif
       if( (0<rrank).and.(irflag.eq.0) )then
          call nrk_vec_comm_ovl2_irecv(&
               commVec,rbuf,rrank,r_req,r_sta)
          irflag=1
       endif

       !$OMP PARALLEL PRIVATE( &
       !$OMP istart_rrank,iend_rrank,j,i,k,kk,Col)

       !$OMP MASTER
       if((isflag.eq.1).or.(irflag.eq.1))then
          do 
             if(isflag.eq.1)then
                call MPI_TEST(s_req,isendflag,s_sta,ierr)
                if(isendflag)then
                   isflag=0
                   srank_count=srank_count+1
#if 0
                   srank=srank+1
#else
                   srank=mod(send_rank0+srank_count+1,srankMax)+1
#endif
                   exit
                endif
             endif
             if(irflag.eq.1)then
                call MPI_TEST(r_req,irecvflag,r_sta,ierr)
                if(irecvflag)then
                   irflag=0
                   rrank_count=rrank_count+1

                   rrank_calc=rrank
#if 0
                   rrank  = rrank+1
#else
                   rrank=recv_rank0-rrank_count+1
                   if(rrank<=0)then
                      rrank=rrankMax+rrank
                   endif
#endif
                   exit
                endif

             endif
          enddo
       endif
       !$OMP END MASTER

       if(rrank0.eq.0)then
          !$OMP DO SCHEDULE(dynamic,nBlock)
          do k=1,n
             Ax(k)=0.0d0
             do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
                if(crsMat%crsCol(j)<istart)then
                   cycle
                endif
                if( iend<crsMat%crsCol(j) )then
                   exit
                endif
                Col=crsMat%crsCol(j)-istart+1
                Ax(k)=Ax(k)+crsMat%crsA(j)*x(Col)
             enddo
          enddo
       endif
       if( 0<rrank0 )then
          !$OMP DO SCHEDULE(dynamic,nBlock)
          do k=1,n
             do kk = &
                  commVec%nrk_prm(rrank0)%nrk_ptr(k),&
                  commVec%nrk_prm(rrank0)%nrk_ptr(k+1)-1
                j=commVec%nrk_prm(rrank0)%nrkCol(kk)
                i=commVec%nrk_prm(rrank0)%nrkRow(kk)
                Ax(k)=Ax(k)+crsMat%crsA(j)*rbuf0(i)
             enddo
          enddo
       endif

       !$OMP END PARALLEL 

       if((isflag.eq.-1).and.(irflag.eq.-1))then
          if(rrank0.eq.-1)then
             exit
          endif
       endif

       if(irflag.eq.0)then
#if 0
          rrank0=rrank-1
#else
          rrank0=rrank_calc
#endif
          !$OMP PARALLEL DO PRIVATE(i)
          do i=1,commVec%nrk_data_max
             rbuf0(i)=rbuf(i)
          enddo
          !$OMP END PARALLEL DO
          rrank_calc=-1
       else
          rrank0=-1
       endif

       if(srankMax<srank_count)then
          isflag=-1
          srank=-1
       endif
       if(rrankMax<rrank_count)then
          irflag=-1
          rrank=-1
       endif

    enddo
    !     call system("sleep 10")
    !     stop
  end subroutine MatVec_crs_ovl2

  subroutine init_MatVec_crs_ovl2(prm,crsMat,commVec,nBlock)
    use matrix
    use comm
    use comm_kernel

    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(CRS_Mat) :: crsMat 
    type(comm_Vec) :: commVec

    integer(kind=intType) i,j,k
    integer(kind=intType) n,istart,iend
    integer(kind=intType) Col

    integer irank
    integer istart_rrank,iend_rrank  
    integer(kind=intType) offset

    integer(kind=intType) nnz
    integer(kind=intType) nrkCol(prm%nnz),nrkRow(prm%nnz)
    integer nBlock

    n = prm%n
    istart=prm%istart
    iend=prm%iend

    !$OMP PARALLEL DO &
    !$OMP PRIVATE( &
    !$OMP irank, &
    !$OMP istart_rrank,iend_rrank  , &
    !$OMP i,j,k, &
    !$OMP offset,Col, &
    !$OMP nnz, &
    !$OMP nrkCol,nrkRow )  SCHEDULE(dynamic,nBlock) 
    do irank=1,commVec%nrk_rank
       nnz=1
       do i=1,prm%nnz
          nrkCol(i)=-1
          nrkRow(i)=-1
       enddo
       allocate(commVec%nrk_prm(irank)%nrk_ptr(n+1))

       istart_rrank=commVec%istart( commVec%nrk_prm(irank)%rank+1    )
       iend_rrank  =commVec%istart( commVec%nrk_prm(irank)%rank+1 +1 )-1
       do k=1,n
          commVec%nrk_prm(irank)%nrk_ptr(k)=nnz
          do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
             if( (istart_rrank<=crsMat%crsCol(j)) .and. &
                  (crsMat%crsCol(j)<=iend_rrank  ) )then
                do i=1,commVec%nrk_prm(irank)%ndata
                   offset=commVec%nrk_prm(irank)%offset(i)
                   Col=istart_rrank+offset
                   if( Col.eq.crsMat%crsCol(j) )then
                      ! Ax(k)=Ax(k)+crsMat%crsA(j)*rbuf0(i)
                      nrkCol(nnz)=j
                      nrkRow(nnz)=i
                      nnz=nnz+1
                      exit
                   endif
                enddo
             endif
             if( iend_rrank<crsMat%crsCol(j) )then
                exit
             endif
          enddo
       enddo
       commVec%nrk_prm(irank)%nrk_ptr(n+1)=nnz

       nnz=nnz-1
       commVec%nrk_prm(irank)%nnz=nnz
       allocate(commVec%nrk_prm(irank)%nrkRow(nnz))
       allocate(commVec%nrk_prm(irank)%nrkCol(nnz))
       do i=1,nnz
          commVec%nrk_prm(irank)%nrkRow(i)=nrkRow(i)
          commVec%nrk_prm(irank)%nrkCol(i)=nrkCol(i)
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine Init_MatVec_crs_ovl2

  subroutine MatVec_dia0(prm,diaMat,Ax,nBlock)
    use matrix
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(DIA_Mat) :: diaMat 
    real*8 Ax(prm%n)
    integer(kind=intType) j,k,kk
    integer(kind=intType) n,nnd,istart
    integer(kind=intType) offset
    integer(kind=intType) kkmax
    integer nBlock

    n      = prm%n
    nnd    = prm%nnd
    istart = prm%istart

    !$OMP PARALLEL DO PRIVATE(kk,k,j,offset,kkmax) SCHEDULE(static,1)
    do kk=1,n,nBlock
       kkmax=min(n,kk+nBlock-1)
       do k=kk,kkmax
          Ax(k)=0.0d0
       enddo
       do j=1,nnd
          offset=diaMat%offsetVec(j)
          do k=kk,kkmax
             ! -----
             ! write(*,*)"  DIA  Mat,Vec",&
             !      diaMat%diaA(k+(j-1)*n) , &
             !      diaMat%Vec(k+diaMat%offsetVec(j))
             ! -----
             ! write(*,*)"  DIA  Mat,Vec",&
             !      diaMat%diaA(k+(j-1)*n) , &
             !      diaMat%Vec(k+diaMat%offsetVec(j))
             ! -----
             Ax(k)=Ax(k) + diaMat%diaA(k+(j-1)*n) * &
                  diaMat%Vec(k+offset) 
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine MatVec_dia0

  subroutine MatVec_dia(prm,diaMat,x,Ax,commVec,rbuf,nBlock)
    use matrix
    use comm
    use comm_kernel
    use vec_kernel
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(DIA_Mat) :: diaMat 
    real*8 x(prm%n)
    real*8 Ax(prm%n)
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
    real*8 rbuf(commVec%nrk_data+1)
    integer(kind=intType) iLs,iLe
    integer nBlock

    n      = prm%n
    nnd    = prm%nnd
    istart = prm%istart
    myrank = commVec%myrank

    !$OMP PARALLEL PRIVATE( &
    !$OMP kk,kkmax,j,irank,irank_start,&
    !$OMP i0,i,k,Col,Row,offset,iLs,iLe,&
    !$OMP minCol,maxCol,minRow,maxRow,kmin,kmax &
    !$OMP ) 
    
    !$OMP DO SCHEDULE(dynamic,1)
    do kk=1,n,nBlock
       kkmax=min(n,kk+nBlock-1)
       do k=kk,kkmax
          Ax(k)=0.0d0
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
             Ax(k-offset) = &
                  Ax(k-offset) + diaMat%diaA(k-offset + (j-1)*n) * x(k)
          enddo
          do k=kk,kkmax
             i = diaMat%iLs(  k + (j-1)*n  )
             Ax(k)=Ax(k) + diaMat%diaA(k+(j-1)*n) * rbuf(i)
          enddo
       enddo

       !        do j=1,nnd
       !           do k=kk,kkmax
       !              i = diaMat%iLs(  k + (j-1)*n  )
       !              if( 0 < i )then
       !                 Ax(k)=Ax(k) + diaMat%diaA(k+(j-1)*n) * rbuf(i)
       !              endif
       !           enddo
       !        enddo

    enddo

    !$OMP END PARALLEL

    !    stop
  end subroutine MatVec_dia

#if 0
  subroutine init_MatVec_dia(prm,diaMat,commVec)
    use matrix
    use comm
    use comm_kernel
    use vec_kernel
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(DIA_Mat) :: diaMat 
    type(comm_Vec) :: commVec

    integer(kind=intType) i,j,k,kk
    integer(kind=intType) n,nnd,istart
    integer(kind=intType) offset
    integer(kind=intType) Row,Col
    integer(kind=intType) maxRow,minRow
    integer(kind=intType) maxCol,minCol
    integer(kind=intType) i0
    integer myrank 
    integer irank,irank_start
    integer(kind=intType) kkmax,kmax,kmin
    real*8 stime,etime

    stime=mpi_wtime()
    n      = prm%n
    nnd    = prm%nnd
    istart = prm%istart
    myrank = commVec%myrank

    !$OMP PARALLEL PRIVATE( &
    !$OMP kk,kkmax,j,irank,irank_start,&
    !$OMP i0,i,k,Col,Row,offset,&
    !$OMP minCol,maxCol,minRow,maxRow,kmin,kmax &
    !$OMP ) 

    !$OMP DO SCHEDULE(dynamic,1)
    do kk=1,n,nBlock
       kkmax=min(n,kk+nBlock-1)
       do k=kk,kkmax
          do j=1,nnd

             Col = k + istart+diaMat%offset(j) - 1
             diaMat%iLs(  k + (j-1)*n  ) = commVec%nrk_data +1
             do irank=1,commVec%nrk_rank
                i0=0
                do i=1,irank-1
                   i0=i0+commVec%nrk_prm(i)%ndata
                enddo
                irank_start=commVec%istart( commVec%nrk_prm(irank)%rank+1 )
                do i=1,commVec%nrk_prm(irank)%ndata
                   Row = &
                        irank_start + &
                        commVec%nrk_prm(irank)%offset(i) 
                   if(Row.eq.Col)then
                      diaMat%iLs(  k + (j-1)*n  )=i+i0
                      exit
                   endif
                enddo
                if(Row.eq.Col)then
                   exit
                endif
             enddo

          enddo
       enddo

    enddo
    !$OMP END PARALLEL
    etime=mpi_wtime()
    !write(*,*)"init_matvec_time=",etime-stime,"[s]"

    !    stop
  end subroutine init_MatVec_dia
#else
  subroutine init_MatVec_dia(prm,diaMat,commVec,nBlock)
    use matrix
    use comm
    use comm_kernel
    use vec_kernel
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(DIA_Mat) :: diaMat 
    type(comm_Vec) :: commVec
    integer nBlock

    integer(kind=intType) i,j,k,kk
    integer(kind=intType) n,nnd,istart,iend
    integer(kind=intType) offset
    integer(kind=intType) Row,Col
    integer(kind=intType) maxRow,minRow
    integer(kind=intType) maxCol,minCol
    integer(kind=intType) i0
    integer myrank 
    integer irank,irank_start
    integer(kind=intType) kkmax,kmax,kmin
    real*8 stime,etime
    integer tmp,tmp_rank,pre_rank

    n      = prm%n
    nnd    = prm%nnd
    istart = prm%istart
    iend   = prm%iend
    myrank = commVec%myrank

    !$OMP PARALLEL PRIVATE( &
    !$OMP kk,kkmax,j,irank,irank_start,&
    !$OMP i0,i,k,Col,Row,offset,tmp,tmp_rank,pre_rank,&
    !$OMP minCol,maxCol,minRow,maxRow,kmin,kmax &
    !$OMP ) 

    !$OMP DO SCHEDULE(dynamic,1)
    do kk=1,n,nBlock
       kkmax=min(n,kk+nBlock-1)
       do k=kk,kkmax
          tmp=1
          tmp_rank=1
          pre_rank=1
          i0=0
          do j=1,nnd
            
             Col = k + istart+diaMat%offset(j) - 1
             diaMat%iLs(  k + (j-1)*n  ) = commVec%nrk_data +1
             if (istart<=Col .and. Col<=iend) then
               cycle
             endif
             do irank=tmp_rank,commVec%nrk_rank
                if (pre_rank<irank) then
                  tmp=1
                endif
                do i=pre_rank,irank-1
                   i0=i0+commVec%nrk_prm(i)%ndata
                enddo
                pre_rank=irank
                irank_start=commVec%istart( commVec%nrk_prm(irank)%rank+1 )
                do i=tmp,commVec%nrk_prm(irank)%ndata
                   Row = &
                        irank_start + &
                        commVec%nrk_prm(irank)%offset(i) 
                   tmp=i
                   tmp_rank=irank

                   !if(Row>Col) then
                   !   exit
                   !endif
                   if(Row.eq.Col)then
                      diaMat%iLs(  k + (j-1)*n  )=i+i0
                      exit
                   endif
                enddo
                !if(Row>Col) then
                !   exit
                !endif
                if(Row.eq.Col)then
                   exit
                endif
             enddo

          enddo
       enddo

    enddo
    !$OMP END PARALLEL

    !    stop
  end subroutine init_MatVec_dia
  
#endif
  
  subroutine init_MatVec_dia_org(prm,diaMat,commVec,nBlock)
    use matrix
    use comm
    use comm_kernel
    use vec_kernel
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(DIA_Mat) :: diaMat 
    type(comm_Vec) :: commVec
    integer nBlock

    integer(kind=intType) i,j,k,kk
    integer(kind=intType) n,nnd,istart
    integer(kind=intType) offset
    integer(kind=intType) Row,Col
    integer(kind=intType) maxRow,minRow
    integer(kind=intType) maxCol,minCol
    integer(kind=intType) i0
    integer myrank 
    integer irank,irank_start
    integer(kind=intType) kkmax,kmax,kmin

    n      = prm%n
    nnd    = prm%nnd
    istart = prm%istart
    myrank = commVec%myrank

    !$OMP PARALLEL PRIVATE( &
    !$OMP kk,kkmax,j,irank,irank_start,&
    !$OMP i0,i,k,Col,Row,offset,&
    !$OMP minCol,maxCol,minRow,maxRow,kmin,kmax &
    !$OMP ) 

    !$OMP DO SCHEDULE(dynamic,1)
    do kk=1,n,nBlock
       kkmax=min(n,kk+nBlock-1)
#if 0
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
                commVec%nrk_prm(irank)%iLs_Ovl2( k + (j-1)*n ) = 1
                do i=1,commVec%nrk_prm(irank)%ndata
                   Row = &
                        irank_start + &
                        commVec%nrk_prm(irank)%offset(i) 
                   Col = k + istart+diaMat%offset(j) - 1
                   if(Row.eq.Col)then
                      commVec%nrk_prm(irank)%iLs_Ovl2(  k + (j-1)*n  )=i
                      exit
                   endif
                enddo
                commVec%nrk_prm(irank)%iLe_Ovl2( k + (j-1)*n ) = commVec%nrk_prm(irank)%ndata
                do i=commVec%nrk_prm(irank)%ndata,1,-1
                   Row = &
                        irank_start + &
                        commVec%nrk_prm(irank)%offset(i) 
                   Col = k + istart+diaMat%offset(j) - 1
                   if(Row.eq.Col)then
                      commVec%nrk_prm(irank)%iLe_Ovl2(  k + (j-1)*n  )=i
                      exit
                   endif
                enddo

             enddo
          enddo
       enddo
#else
       do k=kk,kkmax
          do j=1,nnd

             Col = k + istart+diaMat%offset(j) - 1
             !              diaMat%iLs(  k + (j-1)*n  ) = -1
             diaMat%iLs(  k + (j-1)*n  ) = commVec%nrk_data +1
             do irank=1,commVec%nrk_rank
                i0=0
                do i=1,irank-1
                   i0=i0+commVec%nrk_prm(i)%ndata
                enddo
                irank_start=commVec%istart( commVec%nrk_prm(irank)%rank+1 )
                do i=1,commVec%nrk_prm(irank)%ndata
                   Row = &
                        irank_start + &
                        commVec%nrk_prm(irank)%offset(i) 
                   if(Row.eq.Col)then
                      !                   Ax(k)=Ax(k) + diaMat%diaA(k+(j-1)*n) * rbuf(i+i0)
                      diaMat%iLs(  k + (j-1)*n  )=i+i0
                      exit
                   endif
                enddo
                if(Row.eq.Col)then
                   exit
                endif
             enddo

          enddo
       enddo
#endif

    enddo

    !$OMP END PARALLEL

    !    stop
  end subroutine init_MatVec_dia_org

  subroutine MatVec_dia_ovl1(prm,diaMat,x,Ax,commVec,rbuf,nBlock)
    use matrix
    use comm
    use comm_kernel
    use vec_kernel
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(DIA_Mat) :: diaMat 
    real*8 x(prm%n)
    real*8 Ax(prm%n)
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

    real*8 sbuf(commVec%nrk_data_send)
    real*8 rbuf(commVec%nrk_data+1)

    integer(kind=intType) iLs,iLe

    n      = prm%n
    nnd    = prm%nnd
    istart = prm%istart
    myrank = commVec%myrank

    !     if(commVec%myrank.eq.0)then
    !        write(*,*)"@ dia ovl 1 @"
    !     endif
    call nrk_vec_pack( &
         commVec%nrk_prm_send, &
         commVec%n, &
         commVec%nrk_data_send, &
         x, &
         sbuf, &
         commVec%nrk_rank_send)

    !$OMP PARALLEL PRIVATE( &
    !$OMP kk,kkmax,j, &
    !$OMP offset,minCol,maxCol,minRow,maxRow, &
    !$OMP kmin,kmax,k &
    !$OMP ) 

    !$OMP MASTER
    call nrk_vec_comm(commVec,sbuf,rbuf)
    !$OMP END MASTER

    !$OMP DO SCHEDULE(dynamic,1)
    do kk=1,n,nBlock
       kkmax=min(n,kk+nBlock-1)
       do k=kk,kkmax
          Ax(k)=0.0d0
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
             Ax(k-offset)=Ax(k-offset) + diaMat%diaA(k-offset + (j-1)*n) * x(k)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL

    !$OMP PARALLEL PRIVATE( &
    !$OMP kk,kkmax,j,k,i &
    !$OMP ) 

    !$OMP DO SCHEDULE(static,1)
    do kk=1,n,nBlock
       kkmax=min(n,kk+nBlock-1)
       do j=1,nnd
          do k=kk,kkmax
             i =  diaMat%iLs(  k + (j-1)*n  )
             Ax(k)=Ax(k) + diaMat%diaA(k+(j-1)*n) * rbuf(i)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL

    !    stop
  end subroutine MatVec_dia_ovl1

  subroutine MatVec_dia_ovl2(prm,diaMat,x,Ax,commVec,nBlock)
    use matrix
    use comm
    use comm_kernel
    use vec_kernel
    use PARCEL_integer
    implicit none

    type(prmMat)  :: prm 
    type(DIA_Mat) :: diaMat 
    real*8 x(prm%n)
    real*8 Ax(prm%n)
    type(comm_Vec) :: commVec
    integer nBlock

    real*8 sbuf(commVec%nrk_data_max_send)
    real*8 rbuf(commVec%nrk_data_max+1)
    real*8 rbuf0(commVec%nrk_data_max+1)
    real*8 buff(commVec%nrk_data)
    integer(kind=intType) i,j,k,kk
    integer(kind=intType) n,nnd,istart,iend
    integer(kind=intType) offset
    ! nBlock --> module matrix 
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

    integer(kind=intType) iLs,iLe
    integer recv_rank0,send_rank0

    integer s_req
    integer s_sta(MPI_STATUS_SIZE)
    integer r_req
    integer r_sta(MPI_STATUS_SIZE)
    integer ierr
    integer isflag,irflag
    LOGICAL isendflag,irecvflag

    integer rrank_count,srank_count
    integer rrank_calc

    srankMax=commVec%nrk_rank_send
    rrankMax=commVec%nrk_rank

    n      = prm%n
    nnd    = prm%nnd
    istart = prm%istart
    iend=prm%iend
    myrank = commVec%myrank

    !     if(commVec%myrank.eq.0)then
    !        write(*,*)"@ dia ovl 2 @"
    !     endif

    ! -----
    send_rank0=1
    do irank=1,commVec%nrk_rank_send
       if(commVec%myrank < commVec%nrk_prm_send(irank)%rank)then
          send_rank0=irank
          exit
       endif
    enddo

    recv_rank0=commVec%nrk_rank
    do irank=commVec%nrk_rank,1,-1
       if(commVec%nrk_prm(irank)%rank < commVec%myrank   )then
          exit
       endif
       recv_rank0=irank
    enddo

    rrank_calc=-1
    if(commVec%nrk_rank<1)then
       rrank =-1
       irflag=-1
       rrank_count=rrankMax+1
    else
       irflag=0
       rrank_count=1
       rrank=recv_rank0-rrank_count+1
       if(rrank<=0)then
          rrank=rrankMax+rrank
       endif
    endif

    if(commVec%nrk_rank_send<1)then
       srank =-1
       isflag=-1
       srank_count=srankMax+1
    else
       isflag=0
       srank_count=1
       srank=mod(send_rank0+srank_count+1,srankMax)+1
    endif
    rrank0=0 ! 0は自分の領域
    ! isflag ( 全送信終了:-1,  送信完了:0, 送信中:1)
    ! irflag ( 全受信終了:-1,  受信完了:0, 受信中:1)

    do 
       if( (0<srank).and.(isflag.eq.0) )then
          call nrk_vec_pack_ovl2( &
               commVec%nrk_prm_send, &
               commVec%n, &
               commVec%nrk_data_send, &
               x, &
               sbuf, &
               commVec%nrk_rank_send,&
               srank)
          call nrk_vec_comm_ovl2_isend( &
               commVec,sbuf,srank,s_req,s_sta)
          isflag=1
       endif
       if( (0<rrank).and.(irflag.eq.0) )then
          call nrk_vec_comm_ovl2_irecv(&
               commVec,rbuf,rrank,r_req,r_sta)
          irflag=1
       endif

       !$OMP PARALLEL PRIVATE( &
       !$OMP kk,kkmax,j,irank_start,&
       !$OMP i0,i,k,Col,Row,offset,iLs,iLe,&
       !$OMP minCol,maxCol,minRow,maxRow,kmin,kmax &
       !$OMP ) 

       !$OMP MASTER
       if((isflag.eq.1).or.(irflag.eq.1))then
          do 
             if(isflag.eq.1)then
                call MPI_TEST(s_req,isendflag,s_sta,ierr)
                if(isendflag)then
                   isflag=0
                   srank_count=srank_count+1
                   srank=mod(send_rank0+srank_count+1,srankMax)+1
                   exit
                endif
             endif
             if(irflag.eq.1)then
                call MPI_TEST(r_req,irecvflag,r_sta,ierr)
                if(irecvflag)then
                   irflag=0
                   rrank_count=rrank_count+1

                   rrank_calc=rrank
                   rrank=recv_rank0-rrank_count+1
                   if(rrank<=0)then
                      rrank=rrankMax+rrank
                   endif
                   exit
                endif

             endif
          enddo
       endif
       !$OMP END MASTER

       if(rrank0.eq.0)then
          !$OMP DO SCHEDULE(dynamic,1)
          do kk=1,n,nBlock
             kkmax=min(n,kk+nBlock-1)
             do k=kk,kkmax
                Ax(k)=0.0d0
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
                   Ax(k-offset)=Ax(k-offset) + diaMat%diaA(k-offset + (j-1)*n) * x(k)
                enddo
             enddo
          enddo
       endif

       if( 0<rrank0 )then
#if 0
          irank_start=commVec%istart( commVec%nrk_prm(rrank0)%rank+1 )
          !$OMP DO SCHEDULE(dynamic,1)
          do kk=1,n,nBlock
             kkmax=min(n,kk+nBlock-1)
             do j=1,nnd
                do k=kk,kkmax
                   Col = k + istart+diaMat%offset(j) - 1
                   do i=1,commVec%nrk_prm(rrank0)%ndata
                      Row = &
                           irank_start + &
                           commVec%nrk_prm(rrank0)%offset(i) 
                      if(Row.eq.Col)then
                         Ax(k)=Ax(k) + diaMat%diaA(k+(j-1)*n) * rbuf0(i)
                      endif
                   enddo
                enddo
             enddo
          enddo
#else

          !$OMP DO SCHEDULE(dynamic,nBlock)
          do k=1,n
             do kk = &
                  commVec%nrk_prm(rrank0)%nrk_ptr(k),&
                  commVec%nrk_prm(rrank0)%nrk_ptr(k+1)-1
                j=commVec%nrk_prm(rrank0)%nrkCol(kk)
                i=commVec%nrk_prm(rrank0)%nrkRow(kk)
                Ax(k)=Ax(k) + diaMat%diaA(j) * rbuf0(i)
             enddo
          enddo

#endif
       endif

       !$OMP END PARALLEL 

       if((isflag.eq.-1).and.(irflag.eq.-1))then
          if(rrank0.eq.-1)then
             exit
          endif
       endif

       if(irflag.eq.0)then
          rrank0=rrank_calc
          !$OMP PARALLEL DO PRIVATE(i)
          do i=1,commVec%nrk_data_max
             rbuf0(i)=rbuf(i)
          enddo
          !$OMP END PARALLEL DO
          rrank_calc=-1
       else
          rrank0=-1
       endif

       if(srankMax<srank_count)then
          isflag=-1
          srank=-1
       endif
       if(rrankMax<rrank_count)then
          irflag=-1
          rrank=-1
       endif

    enddo

  end subroutine MatVec_dia_ovl2

  subroutine init_MatVec_dia_ovl2(prm,diaMat,commVec,nBlock)
    use matrix
    use comm
    use comm_kernel
    use vec_kernel
    use PARCEL_integer
    implicit none

    type(prmMat)  :: prm 
    type(DIA_Mat) :: diaMat 
    type(comm_Vec) :: commVec
    integer nBlock

    integer(kind=intType) i,j,k,kk
    integer(kind=intType) n,nnd,istart,iend
    integer(kind=intType) offset
    ! nBlock --> module matrix 
    integer(kind=intType) Row,Col
    integer myrank 
    integer irank_start
    integer(kind=intType) kkmax,kmax,kmin
    integer irank,irankMax
    integer rrankMax
    integer(kind=intType) nnz
    integer(kind=intType) nrkCol(prm%nnd*prm%n),nrkRow(prm%nnd*prm%n)

    n      = prm%n
    nnd    = prm%nnd
    istart = prm%istart
    iend=prm%iend
    myrank = commVec%myrank

    rrankMax=commVec%nrk_rank

    !$OMP PARALLEL DO &
    !$OMP PRIVATE( &
    !$OMP irank,i,irank_start,kk,kkmax,j,Col,Row, &
    !$OMP nnz, &
    !$OMP nrkCol,nrkRow )
    do irank=1,rrankMax
       nnz=1
       do i=1,prm%nnz
          nrkCol(i)=-1
          nrkRow(i)=-1
       enddo
       allocate(commVec%nrk_prm(irank)%nrk_ptr(n+1))
       irank_start=commVec%istart( commVec%nrk_prm(irank)%rank+1 )

       do kk=1,n,nBlock
          kkmax=min(n,kk+nBlock-1)
          do k=kk,kkmax
             commVec%nrk_prm(irank)%nrk_ptr(k)=nnz
             do j=1,nnd
                Col = k + istart+diaMat%offset(j) - 1
                do i=1,commVec%nrk_prm(irank)%ndata
                   Row = &
                        irank_start + &
                        commVec%nrk_prm(irank)%offset(i) 
                   if(Row.eq.Col)then
                      ! Ax(k)=Ax(k)+diaMat%diaA(k+(j-1)*n)*rbuf0(i)
                      nrkCol(nnz)=k+(j-1)*n
                      nrkRow(nnz)=i
                      nnz=nnz+1
                      exit
                   endif
                enddo
             enddo
          enddo
       enddo

       commVec%nrk_prm(irank)%nrk_ptr(n+1)=nnz
       nnz=nnz-1
       commVec%nrk_prm(irank)%nnz=nnz
       allocate(commVec%nrk_prm(irank)%nrkRow(nnz))
       allocate(commVec%nrk_prm(irank)%nrkCol(nnz))
       do i=1,nnz
          commVec%nrk_prm(irank)%nrkRow(i)=nrkRow(i)
          commVec%nrk_prm(irank)%nrkCol(i)=nrkCol(i)
       enddo

    enddo

  end subroutine Init_MatVec_dia_ovl2

  subroutine get_DInv_dia(prm,diaMat,Dinv)
    use matrix
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(DIA_Mat) :: diaMat 
    real*8 Dinv(prm%n)

    integer(kind=intType) i,j,k
    integer(kind=intType) n,nnd

    n = prm%n
    nnd = prm%nnd

    do i=1,nnd
       if(diaMat%offset(i).eq.0)then
          !$OMP PARALLEL DO PRIVATE(k) 
          do k=1,n
             if(diaMat%diaA(k+(i-1)*n).eq.0.0d0)then
                Dinv(k)=1.0d0/1.0d-300
             else
                Dinv(k)=1.0d0/diaMat%diaA(k+(i-1)*n)
             endif
          enddo
          !$OMP END PARALLEL DO
       endif
    enddo

  end subroutine get_DInv_dia

  subroutine get_DInv_crs(prm,crsMat,Dinv)
    use matrix
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(CRS_Mat) :: crsMat 
    real*8 Dinv(prm%n)

    integer(kind=intType) i,j,k
    integer(kind=intType) n,istart,iend
    integer(kind=intType) Row,Col
    integer(kind=intType) nColL,nColU

    n = prm%n
    istart=prm%istart
    iend=prm%iend
    nColL=crsMat%nColL
    nColU=crsMat%nColU

    !$OMP PARALLEL DO PRIVATE(j,k,Row,Col) 
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
             if(crsMat%crsA(j).eq.0.0d0)then
                Dinv(k)=1.0d0/1.0d-300
             else
                Dinv(k)=1.0d0/crsMat%crsA(j)
             endif
          endif
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine get_DInv_crs

  subroutine get_AS_region( &
       prm,addL, &
       Lrank,Urank, &
       mnumth, &
       istart_th,iend_th, &
       precon_thblock)
    use omp_lib
    use matrix

    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    integer addL
    integer Lrank,Urank
    integer mnumth
    integer max_numth
    integer(kind=intType) istart_th(mnumth),iend_th(mnumth)

    integer(kind=intType) n,gn
    integer(kind=intType) istart,iend
    integer numth,myth
    integer(kind=intType) n_th
    integer(kind=intType) n_tmp
    integer(kind=intType) n_startth,n_endth
    integer precon_thblock

    n = prm%n
    gn = prm%gn
    istart=prm%istart
    iend=prm%iend

    ! -----
    ! ブロック数よりもスレッド数が大きい場合はスレッド数を制限する
    call get_maxthreads(max_numth)
    call set_precon_maxthreads(n,precon_thblock)

    !$OMP PARALLEL &
    !$OMP PRIVATE( &
    !$OMP numth,myth,n_startth,n_endth,n_tmp &
    !$OMP )
    numth = omp_get_num_threads()
    myth  = omp_get_thread_num()
    n_tmp=1
    call divide_region( &
         myth,numth,n_tmp,n,n_startth,n_endth)
    ! istart_th = istart -1 + n_startth 
    ! iend_th   = istart_th + (n_endth - n_startth) 

#if 1
    istart_th(myth+1) = istart -1 + n_startth - addL
    iend_th(myth+1)   = (istart -1 + n_startth ) &
         + (n_endth - n_startth) + addL
#else
    istart_th(myth+1) = istart -1 + n_startth 
    iend_th(myth+1)   = (istart -1 + n_startth ) &
         + (n_endth - n_startth) 
#endif
    if(istart_th(myth+1)<=0)then
       istart_th(myth+1)=1
    endif
    if( (iend_th(myth+1)>gn) )then
       !    if( (iend_th(myth+1)>n))then
       iend_th(myth+1)=istart+n-1
    endif

    !    write(*,*)"@@",myth, istart_th(myth+1),iend_th(myth+1)
    !$OMP END PARALLEL

!    if(n<=max_numth)then
    call omp_set_num_threads(max_numth)
!    endif

  end subroutine get_AS_region

  subroutine check_AS_region( &
       prm,addL, &
       iflag, &
       precon_thblock)
    use omp_lib
    use matrix

    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    integer addL
    integer max_numth

    integer(kind=intType) n,gn
    integer(kind=intType) istart,iend
    integer(kind=intType) istart_th,iend_th

    integer numth,myth
    integer(kind=intType) n_th
    integer(kind=intType) n_tmp
    integer(kind=intType) n_startth,n_endth
    integer precon_thblock
    integer iflag

    n = prm%n
    gn = prm%gn
    istart=prm%istart
    iend=prm%iend
!    iflag = 0
    ! -----
    call get_maxthreads(max_numth)
    call set_precon_maxthreads(n,precon_thblock)

    !$OMP PARALLEL &
    !$OMP PRIVATE( &
    !$OMP numth,myth,n_startth,n_endth,n_tmp, &
    !$OMP istart_th,iend_th &
    !$OMP ) &
    !$OMP REDUCTION(min:iflag) 

    numth = omp_get_num_threads()
    myth  = omp_get_thread_num()
    n_tmp=1
    call divide_region( &
         myth,numth,n_tmp,n,n_startth,n_endth)
!    write(*,*)n,n_endth - n_startth + 1,n_startth,n_endth


    istart_th = istart -1 + n_startth - addL
    iend_th   = (istart -1 + n_startth ) &
         + (n_endth - n_startth) + addL
    if(istart_th<=0)then
       istart_th=1
    endif
    if( (iend_th>gn) )then
       iend_th=istart+n-1
    endif

!    write(*,*)myth,n_startth,n_endth,istart_th,iend_th
    if( ( (n_endth - n_startth + 1) / 2  ) < addL )then
       iflag = - 1
    endif

    !$OMP END PARALLEL
!    write(*,*)"iflag = ",iflag
    call omp_set_num_threads(max_numth)

  end subroutine check_AS_region

  subroutine get_DILU_DInv_crs( &
       prm, &
       crsMat, &
       Dinv, &
       ILU_method, &
       independ_nvec, &
       precon_thblock)
    use omp_lib
    use matrix

    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(CRS_Mat) :: crsMat 
    real*8 Dinv(prm%n)

    integer(kind=intType) n,istart,iend
    integer numth,myth,max_numth
    integer(kind=intType) n_tmp
    integer(kind=intType) n_th
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th
    integer ILU_method
    integer independ_nvec
    integer precon_thblock

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
    !$OMP numth,myth,n_startth,n_endth,istart_th,iend_th,n_tmp &
    !$OMP )
    numth = omp_get_num_threads()
    myth  = omp_get_thread_num()
    n_tmp=1
    call divide_region( &
         myth,numth,n_tmp,n,n_startth,n_endth)
    istart_th = istart -1 + n_startth 
    iend_th   = istart_th + (n_endth - n_startth) 
!    write(*,*)myth,istart_th,iend_th  

    call DILU_DInv_crs( &
         crsMat, &
         Dinv(n_startth), &
         istart,iend, &
         n_startth,n_endth, &
         istart_th,iend_th, &
         ILU_method, &
         independ_nvec)
    !$OMP END PARALLEL

!    if(n<=max_numth)then
    call omp_set_num_threads(max_numth)
!    endif

  end subroutine get_DILU_DInv_crs

  subroutine get_DILU_DInv_dia( &
       prm, &
       diaMat, &
       Dinv, &
       ILU_method, &
       independ_nvec, &
       precon_thblock)
    use omp_lib
    use matrix

    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(DIA_Mat) :: diaMat 
    real*8 Dinv(prm%n)
    integer ILU_method
    integer independ_nvec

    integer(kind=intType) n,nnd,istart,iend
    integer numth,myth,max_numth
    integer(kind=intType) n_tmp
    integer(kind=intType) n_th
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th
    integer precon_thblock

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
    !$OMP numth,myth,n_startth,n_endth,istart_th,iend_th,n_tmp &
    !$OMP )
    numth = omp_get_num_threads()
    myth  = omp_get_thread_num()

    n_tmp=1
    call divide_region( &
         myth,numth,n_tmp,n,n_startth,n_endth)
    istart_th = istart -1 + n_startth 
    iend_th   = istart_th + (n_endth - n_startth) 
    call DILU_DInv_dia( &
         diaMat,n,nnd, &
         Dinv(n_startth), &
         istart,iend, &
         n_startth,n_endth, &
         istart_th,iend_th, &
         ILU_method, &
         independ_nvec)
    !$OMP END PARALLEL

!    if(n<=max_numth)then
    call omp_set_num_threads(max_numth)
!    endif

  end subroutine get_DILU_DInv_dia

  subroutine DILU_DInv_AS_crs(prm,crsMat,preconAS,mnumth,ILU_method,independ_nvec,precon_thblock)
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

    integer(kind=intType) n,istart,iend
    integer numth,myth,max_numth
    integer(kind=intType) n_th
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th
    integer independ_nvec
    integer precon_thblock

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
    call DILU_DInv_crs( &
         crsMat, &
         preconAS(myth+1)%Dinv(preconAS(myth+1)%n_startth), &
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

  end subroutine DILU_DInv_AS_crs

  subroutine DILU_DInv_crs( &
       crsMat, &
       Dinv, &
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
    real*8 Dinv(n_startth:n_endth)
    integer ILU_method
    integer independ_nvec

    real*8 ,allocatable::D(:)
    integer(kind=intType) i,j,k
    integer(kind=intType) Row,Col

    integer(kind=intType) kr
    integer(kind=intType) jr,rCol,rRow
    integer(kind=intType) jstart,kk
    integer independ_start,independ_end
    integer irank,ierr

    ! -----
    ! D-ILU分解
    ! A = AL + AD + AU ~= K
    ! K = (D+AL)D^-1(D+AU)
    ! 条件1 対角成分が一致     :: D(i,i) = A(i,i) - Sum[j<i]          A(i,j) D(j,j)^-1 A(j,i)
    ! 条件2 行内の要素和が一致 :: D(i,i) = A(i,i) - Sum[j<i] Sum[k>j] A(i,j) D(j,j)^-1 A(j,k)
    ! -----
    !    write(*,*)"n_startth,n_endth",n_startth,n_endth
    
    allocate(D(n_startth:n_endth))
    kk=0
    do k=n_startth,n_endth
      do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
        if(crsMat%crsCol(j)<istart_th) then
          cycle
        endif
        if(iend_th<crsMat%crsCol(j)) then
          cycle
        endif
        Col=crsMat%crsCol(j)-istart+1
        if (Col .eq.k) then
           kk=kk+1
          exit
        endif
      enddo
    enddo
    if (kk/=(n_endth-n_startth+1)) then
       write(*,*)"Error::DILU--Not found non zero in Diagonal element."
       stop
    endif

    ! independ_nvecの範囲で行内を無視する    
    do k=n_startth,n_endth
       ! 独立に回すループサイズ
       independ_start = n_startth + ((k-n_startth)/independ_nvec) * independ_nvec
       independ_end   = min(independ_start+independ_nvec-1,n_endth)

       !       Row=k-1+istart_th
       Row=k-n_startth+istart_th
       jstart=-1
       kk=-1
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          ! Col=crsMat%crsCol(j)-istart_th+1
          Col=crsMat%crsCol(j)-istart+1
          if( crsMat%crsCol(j)<istart_th ) then
             cycle
          endif
          if( iend_th<crsMat%crsCol(j) )then
             cycle
          endif
          
          ! -----
          ! 対角ブロック
          if(jstart.eq.-1)then

             ! 0対角ブロックの開始位置
             jstart=j
          endif
          ! -----
          !          if(k.eq.Col)then
          !          write(*,*)"k =",k,Col
          if(k/=Col)then
            cycle
          endif

          ! -----
          ! 行列の対称位置に対して値があれば総和を計算する
          ! 下側部分を探査する          
          ! 対角要素による初期化
          D(k)=crsMat%crsA(j)
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
                rCol=crsMat%crsCol(i)-istart+1
                ! -----
                ! 上三角部分の総和を計算する
                if( crsMat%crsCol(i)<istart_th ) then
                   cycle
                endif
                if( iend_th<crsMat%crsCol(i) )then
                   cycle
                endif
                ! -----
                ! 対角ブロック
                if(ILU_method.eq.KSP_PRECONDITION_ILUD)then
                   ! 対角成分が一致
                   if(k.eq.rCol)then
                      D(k)=D(k)-crsMat%crsA(jr)*crsMat%crsA(i)/D(kr)
                      exit
                   endif
                else
                   ! 行内の要素和が一致
                   if(kr<rCol)then
                      D(k)=D(k)-crsMat%crsA(jr)*crsMat%crsA(i)/D(kr)
                   endif
                endif

             enddo
          enddo
          exit

       enddo
    enddo

    ! endif
    do k=n_startth,n_endth
       Dinv(k)=1.0d0/D(k)
#if 0
!      !$omp master
#if 1
      call MPI_Comm_rank(MPI_COMM_WORLD,irank,ierr)
#else
      irank=omp_get_thread_num()
#endif
      write(*,*)irank," D -- ",k,D(k)
!      !$omp end master
#endif


    enddo
    deallocate(D)
  end subroutine DILU_DInv_crs

  subroutine DILU_DInv_dia( &
       diaMat,n,nnd, &
       Dinv, &
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
    real*8 Dinv(n_startth:n_endth)
    integer ILU_method
    integer independ_nvec
    integer(kind=intType) n,nnd,nnd_offset0,kk

    real*8,allocatable::D(:)
    integer(kind=intType) i,j,k
    integer(kind=intType) Row,Col,gCol

    integer(kind=intType) kr
    integer(kind=intType) jr,rCol,rRow
    integer(kind=intType) jstart
    integer(kind=intType) grCol
    integer independ_start,independ_end
    integer ierr,irank

    ! -----
    ! D-ILU分解
    ! A = AL + AD + AU ~= K
    ! K = (D+AL)D^-1(D+AU)
    ! 条件1 対角成分が一致     :: D(i,i) = A(i,i) - Sum[j<i]          A(i,j) D(j,j)^-1 A(j,i)
    ! 条件2 行内の要素和が一致 :: D(i,i) = A(i,i) - Sum[j<i] Sum[k>j] A(i,j) D(j,j)^-1 A(j,k)
    ! -----

    allocate(D(n_startth:n_endth))
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
       if(diaMat%diaA(k+(nnd_offset0-1)*n)==0.0) then
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
       D(k)=diaMat%diaA(k+(nnd_offset0-1)*n)
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
                            D(k)=D(k) - &
                                 diaMat%diaA( kr + (i-1)  * n ) * &
                                 diaMat%diaA( k  + (jr-1) * n )   / D(kr)
                         exit
                      endif
                   else
                      ! 行内の要素和が一致
                       if(kr<rCol)then
                            D(k)=D(k) - &
                                 diaMat%diaA( kr + (i-1)  * n ) * &
                                 diaMat%diaA( k  + (jr-1) * n )   / D(kr)
                      endif
                   endif

                enddo
             enddo
             exit
       enddo
    enddo

    do k=n_startth,n_endth
      Dinv(k)=1.0d0/D(k)
#if 0
!     !$omp master
#if 1
      call MPI_Comm_rank(MPI_COMM_WORLD,irank,ierr)
#else
      irank=omp_get_thread_num()
#endif
      write(*,*)irank," D -- ",k,D(k)
!      !$omp end master
#endif

    enddo
    deallocate(D)
    ! stop

  end subroutine DILU_DInv_dia

  subroutine get_DILU_DInv_AS_crs(method,crsMat)
    use mpi
    use krylov
    use vec_kernel
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

    real*8,allocatable  :: sU_crsA(:)
    integer(kind=intType),allocatable :: sU_crsRow_ptr(:) 
    integer(kind=intType),allocatable :: sU_crsCol(:) 

    real*8,allocatable  :: sL_crsA(:)
    integer(kind=intType),allocatable :: sL_crsRow_ptr(:) 
    integer(kind=intType),allocatable :: sL_crsCol(:) 

    real*8 s_DinvL(method%addL),s_DinvU(method%addL)
    real*8 r_DinvL(method%addL),r_DinvU(method%addL)

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
       ! 前処理用の行列を作成する
       method%ipreflag=2
       call get_DILU_DInv_crs(method%prm,crsMat,method%Dinv,method%ILU_method,method%independ_nvec,method%precon_thblock)
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
!       method%ipreflag=0
       method%ipreflag=KSP_PRECONDITION_BLOCKJACOBI
       call get_DILU_DInv_crs(method%prm,crsMat,method%Dinv,method%ILU_method,method%independ_nvec,method%precon_thblock)
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
       ! -----

       allocate(sL_crsA(s_nnzL))
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
                sL_crsA(i)=crsMat%crsA(j)
                sL_crsCol(i)=Col
                i=i+1
             endif
          enddo
       enddo
       sL_crsRow_ptr(addL+1)=i

    else
       allocate(sL_crsA(1))
       allocate(sL_crsRow_ptr(2))
       allocate(sL_crsCol(1))
    endif
    if(myrank.ne.0)then
       call MPI_Irecv(r_nnzU,1,MPI_INTEGER4,Urank,Urank ,icomm,r_Ureq,ierr)
       call MPI_Isend(s_nnzU,1,MPI_INTEGER4,Urank,myrank,icomm,s_Ureq,ierr)

       allocate(sU_crsA(s_nnzU))
       allocate(sU_crsRow_ptr(addL+1))
       allocate(sU_crsCol(s_nnzU))
       i=1
       do k=1,addL
          sU_crsRow_ptr(k)=i
          do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
             Col=crsMat%crsCol(j)
             if( (istartU<=Col).and.(Col<=iendU) )then
                sU_crsA(i)=crsMat%crsA(j)
                sU_crsCol(i)=Col
                i=i+1
             endif
          enddo
       enddo
       sU_crsRow_ptr(addL+1)=i
    else
       allocate(sU_crsA(1))
       allocate(sU_crsRow_ptr(2))
       allocate(sU_crsCol(1))
    endif

    ! -----
    ! メモリの確保
    if(myrank.ne.(npes-1))then
       call MPI_Wait(r_Lreq,r_Lsta,ierr)
       call MPI_Wait(s_Lreq,s_Lsta,ierr)
       allocate(precon_crsMatASL%crsA(r_nnzL))
       allocate(precon_crsMatASL%crsRow_ptr(addL+1))
       allocate(precon_crsMatASL%crsCol(r_nnzL))
       do i=1,r_nnzL
          precon_crsMatASL%crsA(i)=0.0d0
          precon_crsMatASL%crsCol(i)=0
       enddo
       do i=1,addL+1
          precon_crsMatASL%crsRow_ptr(i)=0.0d0
       enddo

       prm_ASL%n      = addL
       prm_ASL%gn     = prm_ASL%n
       prm_ASL%istart = istart - addL
       prm_ASL%iend   = istart - 1
       prm_ASL%nnz    = r_nnzL
    else
       allocate(precon_crsMatASL%crsA(1))
       allocate(precon_crsMatASL%crsRow_ptr(2))
       allocate(precon_crsMatASL%crsCol(1))
       precon_crsMatASL%crsA(1)      =0.0d0
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
       allocate(precon_crsMatASU%crsA(r_nnzU))
       allocate(precon_crsMatASU%crsRow_ptr(addL+1))
       allocate(precon_crsMatASU%crsCol(r_nnzU))

       do i=1,r_nnzU
          precon_crsMatASU%crsA(i)=0.0d0
          precon_crsMatASU%crsCol(i)=0
       enddo
       do i=1,addL+1
          precon_crsMatASU%crsRow_ptr(i)=0.0d0
       enddo
       prm_ASU%n      = addL
       prm_ASU%gn     = prm_ASU%n 
       prm_ASU%istart = istart + addL
       prm_ASU%iend   = istart + 1
       prm_ASU%nnz    = r_nnzU
    else
       allocate(precon_crsMatASU%crsA(1))
       allocate(precon_crsMatASU%crsRow_ptr(2))
       allocate(precon_crsMatASU%crsCol(1))
       precon_crsMatASU%crsA(1)      =0.0d0
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
       call MPI_Irecv(precon_crsMatASU%crsA,r_nnzU,MPI_REAL8,Urank,Urank ,icomm,r_Ureq,ierr)
       call MPI_Isend(sU_crsA              ,s_nnzU,MPI_REAL8,Urank,myrank,icomm,s_Ureq,ierr)
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
       call MPI_Irecv(precon_crsMatASL%crsA,r_nnzL,MPI_REAL8,Lrank,Lrank ,icomm,r_Lreq,ierr)
       call MPI_Isend(sL_crsA              ,s_nnzL,MPI_REAL8,Lrank,myrank,icomm,s_Lreq,ierr)
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

    allocate(method%precon_crsMatAS%crsA(method%prm_AS%nnz))
    allocate(method%precon_crsMatAS%crsRow_ptr(method%prm_AS%n+1))
    allocate(method%precon_crsMatAS%crsCol(method%prm_AS%nnz))
    call zero_initialize(method%prm_AS%nnz      , method%precon_crsMatAS%crsA)

    i0=method%prm_AS%n+1
    !    call zero_initialize_int(method%prm_AS%n+1 , method%precon_crsMatAS%crsRow_ptr)
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
             method%precon_crsMatAS%crsA(i)  = precon_crsMatASU%crsA(j)
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
          method%precon_crsMatAS%crsA(i)  = crsMat%crsA(j)
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
             method%precon_crsMatAS%crsA(i)  = precon_crsMatASL%crsA(j)
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
    allocate(method%ivec_AS(n+addL*2))
    allocate(method%ovec_AS(n+addL*2))

    call zero_initialize(n+addL*2,method%ivec_AS)
    call zero_initialize(n+addL*2,method%ovec_AS)

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
       !@       write(*,*)i," istart_th_AS,iend_th_AS ",istart_th_AS(i),iend_th_AS(i)
       method%preconAS(i)%istart_th = istart_th_AS(i)
       method%preconAS(i)%iend_th   = iend_th_AS(i)
       method%preconAS(i)%n_startth = istart_th_AS(i) - istartAS +1
       method%preconAS(i)%n_endth   = &
            method%preconAS(i)%n_startth + &
            ( iend_th_AS(i)-istart_th_AS(i) )

       allocate(method%preconAS(i)%Dinv( &
            method%preconAS(i)%n_startth:method%preconAS(i)%n_endth))
       allocate(method%preconAS(i)%svec( &
            method%preconAS(i)%n_startth:method%preconAS(i)%n_endth)) 
       allocate(method%preconAS(i)%rvec( &
            method%preconAS(i)%n_startth:method%preconAS(i)%n_endth)) 

       n_th = method%preconAS(i)%n_endth-method%preconAS(i)%n_startth+1
       call zero_initialize(n_th,method%preconAS(i)%Dinv)
       call zero_initialize(n_th,method%preconAS(i)%svec)
       call zero_initialize(n_th,method%preconAS(i)%rvec)
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
         iend_th_AS,&
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

    call DILU_DInv_AS_crs( &
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
    deallocate(sL_crsA)
    deallocate(sL_crsRow_ptr)
    deallocate(sL_crsCol)
    deallocate(sU_crsA)
    deallocate(sU_crsRow_ptr)
    deallocate(sU_crsCol)

    deallocate(precon_crsMatASL%crsA)
    deallocate(precon_crsMatASL%crsRow_ptr)
    deallocate(precon_crsMatASL%crsCol)

    deallocate(precon_crsMatASU%crsA)
    deallocate(precon_crsMatASU%crsRow_ptr)
    deallocate(precon_crsMatASU%crsCol)

  end subroutine get_DILU_DInv_AS_crs

  subroutine free_precon_ILU_crs(method)
    use krylov
    use PARCEL_integer
    implicit none
    type(KSP) :: method
    integer max_numth
    integer(kind=intType) i
    integer numth

    call get_maxthreads(max_numth)
    call set_precon_maxthreads(method%prm%n,method%precon_thblock)
    call get_maxthreads(numth)

    if(allocated(method%preconILU))then
       do i=1,numth
          call free_crsMat(method%preconILU(i)%crsMatL)
          call free_crsMat(method%preconILU(i)%crsMatU)
       enddo
       deallocate(method%preconILU)
    endif

    call omp_set_num_threads(max_numth)

  end subroutine free_precon_ILU_crs

  subroutine free_precon_ILU_dia(method)
    use krylov
    use PARCEL_integer
    implicit none
    type(KSP) :: method
    integer max_numth
    integer(kind=intType) i
    integer numth

    call get_maxthreads(max_numth)
    call set_precon_maxthreads(method%prm%n,method%precon_thblock)
    call get_maxthreads(numth)

    if(allocated(method%preconILU))then
       do i=1,numth
          call free_diaMat(method%preconILU(i)%diaMatL)
          call free_diaMat(method%preconILU(i)%diaMatU)
       enddo
       deallocate(method%preconILU)
    endif

    call omp_set_num_threads(max_numth)

  end subroutine free_precon_ILU_dia

  subroutine free_precon_AS(method)
    use krylov
    use PARCEL_integer
    implicit none
    type(KSP) :: method
    integer max_numth
    integer(kind=intType) i

    integer numth
    call get_maxthreads(max_numth)
    call set_precon_maxthreads(method%prm%n,method%precon_thblock)
    call get_maxthreads(numth)

    do i=1,numth
       deallocate(method%preconAS(i)%Dinv)
       deallocate(method%preconAS(i)%rvec)
       deallocate(method%preconAS(i)%svec)
    enddo
    deallocate(method%preconAS)

    call omp_set_num_threads(max_numth)

  end subroutine free_precon_AS

  ! z(:)=y(:)*x(:)
  subroutine pointJacobi(n,x,y,z)
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    real*8 alpha,beta
    real*8 x(n)
    real*8 y(n)
    real*8 z(n)
    integer(kind=intType) i 

    !$OMP PARALLEL DO PRIVATE(i) 
    do i=1,n
       z(i)=y(i)*x(i)
    enddo
    !$OMP END PARALLEL DO
  end subroutine pointJacobi

  subroutine blockJacobi_crs(method,r,s)
    use omp_lib
    use krylov
    use PARCEL_integer

    implicit none
    type(KSP) :: method
    real*8 r(method%prm%n)
    real*8 s(method%prm%n)

    integer(kind=intType) n
    integer(kind=intType) nnz
    integer(kind=intType) istart,iend

    integer myth,numth,max_numth
    integer(kind=intType) n_th
    integer(kind=intType) n_tmp
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th
    integer precon_thblock

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
    !$OMP numth,myth,n_startth,n_endth,istart_th,iend_th,n_tmp &
    !$OMP )
    numth = omp_get_num_threads()
    myth  = omp_get_thread_num()

    n_tmp=1
    call divide_region( &
         myth,numth,n_tmp,n,n_startth,n_endth)

    istart_th = istart -1 + n_startth 
    iend_th   = istart_th + (n_endth - n_startth) 
    call solve_ILU_crs( &
         n,nnz,&
         method%crsMat%crsA, &
         method%crsMat%crsRow_ptr, &
         method%crsMat%crsCol, &
         istart,iend, &
         istart_th,iend_th, &
         n_startth,n_endth, &
         method%Dinv(n_startth), &
         r(n_startth), &
         s(n_startth), &
         method%ILUregionMat(myth+1)%ILU_Lstart, &
         method%ILUregionMat(myth+1)%ILU_Lend,   &
         method%ILUregionMat(myth+1)%ILU_Ustart, &
         method%ILUregionMat(myth+1)%ILU_Uend,   &
         method%independ_nvec &
         )
    !$OMP END PARALLEL

!    if(n<=max_numth)then
    call omp_set_num_threads(max_numth)
!    endif

  end subroutine blockJacobi_crs

  subroutine init_blockJacobi_crs(method)
    use omp_lib
    use krylov
    use PARCEL_integer
    implicit none
    type(KSP) :: method

    integer(kind=intType) n
    integer(kind=intType) nnz
    integer(kind=intType) istart,iend

    integer numth,myth,max_numth
    integer(kind=intType) n_th
    integer(kind=intType) n_tmp
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th
    integer precon_thblock 

    n=method%prm%n
    nnz=method%prm%nnz
    istart=method%prm%istart
    iend=method%prm%iend
    precon_thblock=method%precon_thblock 

    ! -----
    ! ブロック数よりもスレッド数が大きい場合はスレッド数を制限する
    call get_maxthreads(max_numth)
    call set_precon_maxthreads(n,precon_thblock)

    call get_maxthreads(numth)
    allocate(method%ILUregionMat(numth))

    !$OMP PARALLEL &
    !$OMP PRIVATE( &
    !$OMP numth,myth,n_startth,n_endth,istart_th,iend_th,n_tmp &
    !$OMP )
    numth = omp_get_num_threads()
    myth  = omp_get_thread_num()
    n_tmp=1
    call divide_region( &
         myth,numth,n_tmp,n,n_startth,n_endth)
    istart_th = istart -1 + n_startth 
    iend_th   = istart_th + (n_endth - n_startth) 

    allocate(method%ILUregionMat(myth+1)%ILU_Lstart(n_startth:n_endth))
    allocate(method%ILUregionMat(myth+1)%ILU_Lend(n_startth:n_endth)  )
    allocate(method%ILUregionMat(myth+1)%ILU_Ustart(n_startth:n_endth))
    allocate(method%ILUregionMat(myth+1)%ILU_Uend(n_startth:n_endth)  )

    call init_solve_ILU_crs( &
         n,nnz, &
         method%crsMat%crsRow_ptr,method%crsMat%crsCol, &
         method%independ_nvec, &
         istart,iend, &
         istart_th,iend_th, &
         n_startth,n_endth, &
         method%ILUregionMat(myth+1)%ILU_Lstart, &
         method%ILUregionMat(myth+1)%ILU_Lend,   &
         method%ILUregionMat(myth+1)%ILU_Ustart, &
         method%ILUregionMat(myth+1)%ILU_Uend    &
         )

    !$OMP END PARALLEL

!    if(n<=max_numth)then
    call omp_set_num_threads(max_numth)
!    endif
  end subroutine init_blockJacobi_crs

  subroutine blockJacobi_dia(method,r,s)
    use omp_lib
    use krylov
    use PARCEL_integer
    implicit none
    type(KSP) :: method
    real*8 r(method%prm%n)
    real*8 s(method%prm%n)

    integer(kind=intType) n
    integer(kind=intType) istart,iend

    integer numth,myth,max_numth
    integer(kind=intType) n_th
    integer(kind=intType) n_tmp
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th
    integer(kind=intType) nnd
    integer precon_thblock 

    n      = method%prm%n
    istart = method%prm%istart
    iend   = method%prm%iend
    nnd    = method%prm%nnd 
    precon_thblock=method%precon_thblock 

    ! -----
    ! ブロック数よりもスレッド数が大きい場合はスレッド数を制限する
    call get_maxthreads(max_numth)
    call set_precon_maxthreads(n,precon_thblock)

    !$OMP PARALLEL &
    !$OMP PRIVATE( &
    !$OMP numth,myth,n_startth,n_endth,istart_th,iend_th,n_tmp &
    !$OMP )
    numth = omp_get_num_threads()
    myth  = omp_get_thread_num()
    n_tmp=1
    call divide_region( &
         myth,numth,n_tmp,n,n_startth,n_endth)
    istart_th = istart -1 + n_startth 
    iend_th   = istart_th + (n_endth - n_startth) 

    call solve_ILU_dia( &
         n,nnd,&
         method%diaMat%diaA,method%diaMat%offset, &
         istart,iend, &
         istart_th,iend_th, &
         n_startth,n_endth, &
         method%Dinv(n_startth), &
         r(n_startth), &
         s(n_startth), &
         method%ILUregionMat(myth+1)%ILU_Lstart, &
         method%ILUregionMat(myth+1)%ILU_Lend,   &
         method%ILUregionMat(myth+1)%ILU_Ustart, &
         method%ILUregionMat(myth+1)%ILU_Uend,   &
         method%independ_nvec    &
         )

    !$OMP END PARALLEL

!    if(n<=max_numth)then
    call omp_set_num_threads(max_numth)
!    endif
  end subroutine blockJacobi_dia

  subroutine init_blockJacobi_dia(method)
    use omp_lib
    use krylov
    use PARCEL_integer
    implicit none
    type(KSP) :: method
    integer(kind=intType) n
    integer(kind=intType) istart,iend
    integer(kind=intType) i
    integer numth,myth,max_numth
    integer(kind=intType) n_tmp
    integer(kind=intType) n_th
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th
    integer(kind=intType) nnd
    integer precon_thblock 

    n      = method%prm%n
    istart = method%prm%istart
    iend   = method%prm%iend
    nnd    = method%prm%nnd 
    precon_thblock=method%precon_thblock 

    ! -----
    ! ブロック数よりもスレッド数が大きい場合はスレッド数を制限する
    call get_maxthreads(max_numth)
    call set_precon_maxthreads(n,precon_thblock)

    !$OMP PARALLEL 
    !$OMP MASTER
    numth = omp_get_num_threads()
    !$OMP END MASTER
    !$OMP END PARALLEL 
    allocate(method%ILUregionMat(numth))

    !$OMP PARALLEL &
    !$OMP PRIVATE( &
    !$OMP numth,myth,n_startth,n_endth,istart_th,iend_th,n_tmp &
    !$OMP )
    numth = omp_get_num_threads()
    myth  = omp_get_thread_num()

    n_tmp=1
    call divide_region( &
         myth,numth,n_tmp,n,n_startth,n_endth)
    istart_th = istart -1 + n_startth 
    iend_th   = istart_th + (n_endth - n_startth) 

    allocate(method%ILUregionMat(myth+1)%ILU_Lstart(n_startth:n_endth))
    allocate(method%ILUregionMat(myth+1)%ILU_Lend(n_startth:n_endth)  )
    allocate(method%ILUregionMat(myth+1)%ILU_Ustart(n_startth:n_endth))
    allocate(method%ILUregionMat(myth+1)%ILU_Uend(n_startth:n_endth)  )

    call init_solve_ILU_dia( &
         n,nnd,&
         method%diaMat%offset, &
         method%independ_nvec, &
         istart,iend, &
         istart_th,iend_th, &
         n_startth,n_endth, &
         method%ILUregionMat(myth+1)%ILU_Lstart, &
         method%ILUregionMat(myth+1)%ILU_Lend,   &
         method%ILUregionMat(myth+1)%ILU_Ustart, &
         method%ILUregionMat(myth+1)%ILU_Uend    &
         )

    !$OMP END PARALLEL

!    if(n<=max_numth)then
    call omp_set_num_threads(max_numth)
!    endif

  end subroutine init_blockJacobi_dia

  subroutine  AdditiveSchwarz_crs(method,r,s)
    use krylov
    use vec_kernel
    use PARCEL_integer
    use PARCEL_time
    
    implicit none
    type(KSP) :: method
    real*8 s(method%prm%n)
    real*8 r(method%prm%n)

    real*8 vecL(method%addL)
    real*8 vecU(method%addL)
    real*8 stime,etime
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

    stime=MPI_Wtime()
    call zero_initialize(addL_intType,vecL)
    call zero_initialize(addL_intType,vecU)
    ! ---
    if(method%Urank_AS.ne.-1)then
       i0=addL
    else
       i0=0
    endif

    if( (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_BASIC).or.(iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_RESTRICT) )then
       call set_AS_vec( &
            method%prm%n, &
            method%commVec%comm, &
            method%addL, &
            method%commVec%myrank, &
            method%Urank_AS, &
            method%Lrank_AS, &
            r, &
            method%ivec_AS &
            )
    elseif( (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_INTERPOLATE).or.(iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_NONE) )then
       call zero_initialize(addL_intType,method%ivec_AS(1))
       call xcpy(method%prm%n,r,method%ivec_AS(1+i0))      
       call zero_initialize(addL_intType,method%ivec_AS(1+n+i0))
    endif
    etime=MPI_Wtime()
    precon_time(1) = precon_time(1) + etime - stime

    stime=MPI_Wtime()
    call AdditiveSchwarz_ILU_crs(method,iflagAS)
    etime=MPI_Wtime()
    precon_time(2) = precon_time(2) + etime - stime

    stime=MPI_Wtime()
    call AS_Local( &
         method%prm_AS, &
         method%preconAS, &
         method%max_numth, &
         iflagAS, &
         method%ovec_AS, &
         method%precon_thblock)
    etime=MPI_Wtime()
    precon_time(3) = precon_time(3) + etime - stime
    
    stime=MPI_Wtime()
    call xcpy(method%prm%n,method%ovec_AS(1+i0),s) 
    etime=MPI_Wtime()
    precon_time(4) = precon_time(4) + etime - stime

    stime=MPI_Wtime()
    if((iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_RESTRICT).or.(iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_NONE))then
       return 
    endif

    ! -----
    ! オーバーラップ領域の補正
    ! -----
    call get_AS_vec( &
         method%prm%n, &
         method%commVec%comm, &
         method%addL, &
         method%commVec%myrank, &
         method%Urank_AS, &
         method%Lrank_AS, &
         method%ovec_AS, &
         vecL,vecU &
         )

    if(method%Lrank_AS.ne.-1)then
       do i=1,addL
          s(n-addL+i)= ( s(n-addL+i)+vecL(i) ) 
       enddo
    endif
    if(method%Urank_AS.ne.-1)then
       do i=1,addL 
          s(i)= ( s(i)+vecU(i) ) 
       enddo
    endif
    etime=MPI_Wtime()
    precon_time(5) = precon_time(5) + etime - stime

  end subroutine AdditiveSchwarz_crs

  subroutine  AdditiveSchwarz_dia(method,r,s)
    use krylov
    use vec_kernel
    use PARCEL_integer
    use PARCEL_time

    implicit none
    type(KSP) :: method
    real*8 r(method%prm%n)
    real*8 s(method%prm%n)

    real*8 vecL(method%addL)
    real*8 vecU(method%addL)
    real*8 stime,etime

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

    stime=MPI_Wtime()
    call zero_initialize(addL_intType,vecL)
    call zero_initialize(addL_intType,vecU)
    ! ---
    if(method%Urank_AS.ne.-1)then
       i0=addL
    else
       i0=0
    endif

    if( (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_BASIC).or.(iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_RESTRICT) )then
       call set_AS_vec( &
            method%prm%n, &
            method%commVec%comm, &
            method%addL, &
            method%commVec%myrank, &
            method%Urank_AS, &
            method%Lrank_AS, &
            r, &
            method%ivec_AS &
            )
    elseif( (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_INTERPOLATE).or.(iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_NONE) )then
       call zero_initialize(addL_intType,method%ivec_AS(1))
       call xcpy(method%prm%n,r,method%ivec_AS(1+i0))      
       call zero_initialize(addL_intType,method%ivec_AS(1+n+i0))
    endif
    etime=MPI_Wtime()
    precon_time(1) = precon_time(1) + etime - stime

    !     call AdditiveSchwarz_ILU_dia( &
    !          method%prm_AS, &
    !          method%precon_diaMatAS, &
    !          method%preconAS, &
    !          method%max_numth, &
    !          iflagAS, &
    !          method%ivec_AS)

    stime=MPI_Wtime()
    call AdditiveSchwarz_ILU_dia( &
         method,iflagAS)
    etime=MPI_Wtime()
    precon_time(2) = precon_time(2) + etime - stime

    stime=MPI_Wtime()
    call AS_Local( &
         method%prm_AS, &
         method%preconAS, &
         method%max_numth, &
         iflagAS, &
         method%ovec_AS,&
         method%precon_thblock)
    etime=MPI_Wtime()
    precon_time(3) = precon_time(3) + etime - stime
    
    stime=MPI_Wtime()
    call xcpy(method%prm%n,method%ovec_AS(1+i0),s) 
    etime=MPI_Wtime()
    precon_time(4) = precon_time(4) + etime - stime
    
    stime=MPI_Wtime()
    if((iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_RESTRICT).or.(iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_NONE))then
       ! if( method%commVec%myrank .eq. 0 )then
       !    call view_vec(method%prm%n,s)
       ! endif
       ! call MPI_Barrier(MPI_COMM_WORLD,myrank)
       ! stop
       return 
    endif

    ! -----
    ! オーバーラップ領域の補正
    ! -----
    call get_AS_vec( &
         method%prm%n, &
         method%commVec%comm, &
         method%addL, &
         method%commVec%myrank, &
         method%Urank_AS, &
         method%Lrank_AS, &
         method%ovec_AS, &
         vecL,vecU &
         )

    if(method%Lrank_AS.ne.-1)then
       do i=1,addL
          s(n-addL+i)= ( s(n-addL+i)+vecL(i) ) 
       enddo
    endif
    if(method%Urank_AS.ne.-1)then
       do i=1,addL 
          s(i)= ( s(i)+vecU(i) ) 
       enddo
    endif
    etime=MPI_Wtime()
    precon_time(5) = precon_time(5) + etime - stime

  end subroutine AdditiveSchwarz_dia

  subroutine solve_ILU_crs( &
       n,nnz,&
       crsA,crsRow_ptr,crsCol, &
       istart,iend, &
       istart_th,iend_th, &
       n_startth,n_endth, &
       Dinv,r,s,   &
       ILU_Lstart, &
       ILU_Lend,   &
       ILU_Ustart, &
       ILU_Uend,    &
       independ_nvec &
       )
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) nnz

    real*8 crsA(nnz)
    integer(kind=intType) crsRow_ptr(n+1)
    integer(kind=intType) crsCol(nnz)

    integer(kind=intType) istart,iend
    integer(kind=intType) istart_th,iend_th
    integer(kind=intType) n_startth,n_endth

    real*8 r(n_startth:n_endth)
    real*8 s(n_startth:n_endth)
    real*8 Dinv(n_startth:n_endth)

    integer(kind=intType) j,k,kk,jj
    integer(kind=intType) Row,Col
    integer independ_nvec
!    real*8 z
    real*8 z(independ_nvec),zz

    integer(kind=intType) ILU_Lstart(n_startth:n_endth)
    integer(kind=intType) ILU_Lend(n_startth:n_endth)
    integer(kind=intType) ILU_Ustart(n_startth:n_endth)
    integer(kind=intType) ILU_Uend(n_startth:n_endth)

    integer irank,ierr
    integer(kind=intType)  kk0,kkmin,kkstep

    ! -----
    ! (D+L)D^-1(D+U)s=r

    if(independ_nvec.eq.1)then
       ! -----
       !   z = D^-1(D+U)s
       !   (D+L)z=r
       do k=n_startth,n_endth
          zz=r(k)      
          do j=ILU_Lstart(k),ILU_Lend(k)
             Col=crsCol(j)-istart+1 
             ! r - Lz の計算
             zz=zz-s(Col)*crsA(j)
          enddo
          ! z = D^-1 (r - Lz ) 
          s(k)=Dinv(k)*zz
       enddo
       
       ! -----
       ! D^-1(D+U)s=z
       !   (I+D^-1U)s=z 
       !   s=z - D^-1 Us       
       do k=n_endth,n_startth,-1
          zz=0.0d0
          do j=ILU_Ustart(k),ILU_Uend(k),-1 
             Col=crsCol(j)-istart+1 
             ! Us の計算
             zz=zz+s(Col)*crsA(j)
          enddo
          s(k)=s(k)-Dinv(k)*zz
       enddo
       
    else

       ! -----
       !   z = D^-1(D+U)s
       !   (D+L)z=r
       do kk=n_startth,n_endth,independ_nvec
          do k=kk,min(kk+independ_nvec-1,n_endth)
             z(k-kk+1)=r(k)
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
                Col = crsCol(j)-istart+1 
                ! r - Lz の計算
                z(k-kk+1)=z(k-kk+1)-s(Col)*crsA(j)
             enddo
          enddo
          kk0 = kk0+kkstep+1
          if(kk0>kkmin)then
             exit
          endif
       enddo
#endif

          do k=kk,min(kk+independ_nvec-1,n_endth)
             s(k)=Dinv(k)*z(k-kk+1) 
          enddo
       enddo
    ! -----
    ! D^-1(D+U)s=z
    !   (I+D^-1U)s=z 
    !   s=z - D^-1 Us
       do kk=n_endth,n_startth,-independ_nvec
          do k = kk,max(kk-independ_nvec+1,n_startth),-1
             z(kk-k+1)=0.0d0
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
                Col = crsCol(j)-istart+1 
                z(kk-k+1)=z(kk-k+1)+s(Col)*crsA(j)
             enddo
          enddo

          kk0 = kk0-kkstep-1
          if(kk0<kkmin)then
             exit
          endif
       enddo

#endif
          do k = kk,max(kk-independ_nvec+1,n_startth),-1
             s(k)=s(k)-Dinv(k)*z(kk-k+1)
          enddo
       enddo
    endif
 
  end subroutine solve_ILU_crs

  subroutine init_solve_ILU_crs( &
       n,nnz, &
       crsRow_ptr,crsCol,&
       independ_nvec, &
       istart,iend, &
       istart_th,iend_th, &
       n_startth,n_endth, &
       ILU_Lstart, &
       ILU_Lend,   &
       ILU_Ustart, &
       ILU_Uend    &
       )
    use matrix
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) nnz

    integer(kind=intType) crsRow_ptr(n+1)
    integer(kind=intType) crsCol(nnz)
    integer independ_nvec
    integer(kind=intType) istart,iend
    integer(kind=intType) istart_th,iend_th
    integer(kind=intType) n_startth,n_endth

    integer(kind=intType) j,k
    integer(kind=intType) Row,Col

    integer(kind=intType) ILU_Lstart(n_startth:n_endth)
    integer(kind=intType) ILU_Lend(n_startth:n_endth)
    integer(kind=intType) ILU_Ustart(n_startth:n_endth)
    integer(kind=intType) ILU_Uend(n_startth:n_endth)
    integer independ_start,independ_end
    integer kk0,kkmin,kkstep

#if 0
    ! -----
    ! (D+L)D^-1(D+U)s=r

    ! -----
    !   z = D^-1(D+U)s
    !   (D+L)z=r
    do k=n_startth,n_endth
       ILU_Lstart(k)=crsRow_ptr(k)
       ILU_Lend(k)  =crsRow_ptr(k+1)-1
       Row=k-n_startth+istart_th
       do j=crsRow_ptr(k),crsRow_ptr(k+1)-1
          if( crsCol(j)<istart_th ) then 
             cycle
          endif
          if(Row<=crsCol(j))then
             ILU_Lstart(k)=j
             exit
          endif
          ILU_Lstart(k)=j
          exit
       enddo
       do j=crsRow_ptr(k),crsRow_ptr(k+1)-1
          if( crsCol(j)<istart_th ) then 
             cycle
          endif
          if(Row<=crsCol(j))then
             ILU_Lend(k)=j-1
             exit
          endif
       enddo
    enddo

    ! -----
    ! D^-1(D+U)s=z
    !   (I+D^-1U)s=z 
    !   s=z - D^-1 Us
    do k=n_endth,n_startth,-1
       ILU_Ustart(k)=crsRow_ptr(k+1)-1
       ILU_Uend(k)  =crsRow_ptr(k)
       Row=k-n_startth+istart_th

       do j=crsRow_ptr(k+1)-1,crsRow_ptr(k),-1
          if( iend_th<crsCol(j) )then
             cycle
          endif
          if(Row>=crsCol(j))then
             ILU_Ustart(k)=j
             exit
          endif
          ILU_Ustart(k)=j
          exit
       enddo
       do j=crsRow_ptr(k+1)-1,crsRow_ptr(k),-1
          if( iend_th<crsCol(j) )then
             cycle
          endif
          if(Row>=crsCol(j))then
             ILU_Uend(k)  =j+1
             exit
          endif
       enddo

    enddo

#else
    ! -----
    ! (D+L)D^-1(D+U)s=r
!    write(*,*)"n_startth,n_endth",n_startth,n_endth 
    ! -----
    !   z = D^-1(D+U)s
    !   (D+L)z=r
    do k=n_startth,n_endth
       ! 独立に回すループサイズ
       independ_start = n_startth + ((k-n_startth)/independ_nvec) * independ_nvec
       independ_end   = min(independ_start+independ_nvec-1,n_endth)

#if 0
       ILU_Lstart(k)=crsRow_ptr(k)
       ILU_Lend(k)  =crsRow_ptr(k+1)-1
#else
       ! 下三角要素がない場合は対角要素のみで処理を行い、下三角行列部分の処理をスキップさせる
       ILU_Lstart(k)=crsRow_ptr(k+1)-1
       ILU_Lend(k)  =crsRow_ptr(k)
#endif

       Row=k-n_startth+istart_th
       do j=crsRow_ptr(k),crsRow_ptr(k+1)-1
          Col=crsCol(j)-istart+1
          if( crsCol(j)<istart_th ) then 
             cycle
          endif
!          if(Row<=crsCol(j))then
!          if(independ_start<Col)then
          if(independ_start<=Col)then
             ILU_Lstart(k)=j
             exit
          endif
          ILU_Lstart(k)=j
          exit
       enddo

       do j=crsRow_ptr(k),crsRow_ptr(k+1)-1
          Col=crsCol(j)-istart+1
          if( crsCol(j)<istart_th ) then 
             cycle
          endif
!          if(Row<=crsCol(j))then
          if(independ_start<=Col)then
             ILU_Lend(k)=j-1
             exit
          endif
       enddo
    enddo

    ! -----
    ! D^-1(D+U)s=z
    !   (I+D^-1U)s=z 
    !   s=z - D^-1 Us
    do k=n_endth,n_startth,-1
       ! 独立に回すループサイズ
       independ_start = n_startth + ((k-n_startth)/independ_nvec) * independ_nvec
       independ_end   = min(independ_start+independ_nvec-1,n_endth)

#if 0
       ILU_Ustart(k)=crsRow_ptr(k+1)-1
       ILU_Uend(k)  =crsRow_ptr(k)
#else
       ! 上三角要素がない場合は対角要素のみで処理を行い、上三角行列部分の処理をスキップさせる
       ILU_Ustart(k)= crsRow_ptr(k)
       ILU_Uend(k)  = crsRow_ptr(k+1)-1
#endif
       Row=k-n_startth+istart_th

       do j=crsRow_ptr(k+1)-1,crsRow_ptr(k),-1
          Col=crsCol(j)-istart+1

          if( iend_th<crsCol(j) )then
             cycle
          endif

!          if(Row>=crsCol(j))then
          if(independ_end>=Col)then
             ILU_Ustart(k)=j
             exit
          endif
          ILU_Ustart(k)=j
          exit
       enddo

       do j=crsRow_ptr(k+1)-1,crsRow_ptr(k),-1
          Col=crsCol(j)-istart+1

          if(independ_end<Col)then
             cycle
          endif
!          if(Row>=crsCol(j))then
!             ILU_Uend(k)  =j+1
          if(independ_end>=Col)then
             ILU_Uend(k)  =j+1
             exit
          endif
       enddo

    enddo
#if 0
    !$omp master
    do k=n_startth,n_endth
       if(ILU_Lstart(k)<=ILU_Lend(k))then
          write(*,*)" L val",k,ILU_Lstart(k),ILU_Lend(k)
       endif
    enddo
    write(*,*)"-------------------------"
    do k=n_endth,n_startth,-1
       if(ILU_Uend(k)<=ILU_Ustart(k))then
          write(*,*)" U val",k,ILU_Ustart(k),ILU_Uend(k)  
       endif
    enddo
    !$omp end master
#endif

#endif

#if 0
    do kk=n_startth,n_endth,independ_nvec
       kk0 = kk
       kkmin = min(kk+independ_nvec-1,n_endth)
       do
          kkstep = 0
          do k = kk0,kkmin-1 
             if ( ( ILU_Lend(k)-ILU_Lstart(k) + 1 ).and.( ILU_Lend(k+1)-ILU_Lstart(k+1) + 1 ) )then
                kkstep=kkstep+1
             else
                exit
             endif
          enddo

          do jj=1,(ILU_Lend(kk0)-ILU_Lstart(kk0)+1)
             do k = kk0,kk0+kkstep
                j=ILU_Lstart(k)+jj-1

             enddo
          enddo

          kk0 = kk0+kkstep+1
          if(kk0>kkmin)then
             exit
          endif

       enddo
    enddo
#endif

  end subroutine init_solve_ILU_crs

  subroutine solve_ILU_dia( &
       n,nnd,&
       diaA, &
       offset, &
       istart,iend, &
       istart_th,iend_th, &
       n_startth,n_endth, &
       Dinv,r,s, &
       ILU_Lstart, &
       ILU_Lend,   &
       ILU_Ustart, &
       ILU_Uend,   &
       independ_nvec &
       )
    use matrix
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) nnd
    !    type(DIA_Mat) :: diaMat

    real*8 diaA(nnd*n)
    integer(kind=intType) offset(nnd)
    integer(kind=intType) istart,iend
    integer(kind=intType) istart_th,iend_th
    integer(kind=intType) n_startth,n_endth

    real*8 r(n_startth:n_endth)
    real*8 s(n_startth:n_endth)
    real*8 Dinv(n_startth:n_endth)

    integer(kind=intType) ILU_Lstart(n_startth:n_endth)
    integer(kind=intType) ILU_Lend(n_startth:n_endth)
    integer(kind=intType) ILU_Ustart(n_startth:n_endth)
    integer(kind=intType) ILU_Uend(n_startth:n_endth)

    integer(kind=intType) j,k,jj,kk
    integer(kind=intType) Row,Col
!    real*8 z
    real*8 z(independ_nvec),zz
    integer independ_nvec
    integer irank,ierr
    integer(kind=intType) kk0,kkmin,kkstep
    
    if(independ_nvec.eq.1)then
       ! -----
       ! (D+L)D^-1(D+U)s=r
       
       ! -----
       !   z = D^-1(D+U)s
       !   (D+L)z=r
       
       do k=n_startth,n_endth
          zz=r(k)
!          write(*,*)"in  : DIA -- ",k,r(k)
!          write(*,*)"loop: DIA -- ",k,ILU_Lstart(k),ILU_Lend(k)
          do j=ILU_Lstart(k),ILU_Lend(k)
             ! r - Lz の計算
             zz=zz-s(k+offset(j))*diaA(k+(j-1)*n)
!             write(*,*)"diaA : DIA -- ",diaA(k+(j-1)*n)
          enddo
          ! z = D^-1 (r - Lz) 
          s(k)=Dinv(k)*zz
!          write(*,*)"sav : DIA -- ",Dinv(k),zz
!          write(*,*)"sav : DIA -- ",k,s(k)
       enddo
       
       ! -----
       ! D^-1(D+U)s=z
       !   (I+D^-1U)s=z 
       !   s=z - D^-1 Us
       do k=n_endth,n_startth,-1
          zz=0.0d0
          do j=ILU_Ustart(k),ILU_Uend(k),-1
             zz=zz+s(k+offset(j))*diaA(k+(j-1)*n)
          enddo
          s(k)=s(k)-Dinv(k)*zz
!          write(*,*)"out : DIA -- ",k,s(k)
       enddo
       
    else
       ! -----
       ! (D+L)D^-1(D+U)s=r
       
       ! -----
       !   z = D^-1(D+U)s
       !   (D+L)z=r
       do kk=n_startth,n_endth,independ_nvec
          do k=kk,min(kk+independ_nvec-1,n_endth)
             z(k-kk+1)=r(k)
          enddo
#if 0
          do k=kk,min(kk+independ_nvec-1,n_endth)
             do j=ILU_Lstart(k),ILU_Lend(k)
                ! r - Lz の計算
                z(k-kk+1)=z(k-kk+1)-s(k+offset(j))*diaA(k+(j-1)*n)
             enddo
          enddo
#else
       kk0 = kk
       kkmin = min(kk+independ_nvec-1,n_endth)
       do
          kkstep = 0
          do k = kk0,kkmin-1 
             if ( ( ILU_Lstart(k).eq.ILU_Lstart(k+1) ).and.( ILU_Lend(k).eq.ILU_Lend(k+1) ) )then
!             if ( ( ILU_Lend(k)-ILU_Lstart(k)+1 ).eq.( ILU_Lend(k+1) - ILU_Lstart(k+1) + 1 ) )then
                kkstep=kkstep+1
             else
                exit
             endif
          enddo

#if 1
          do j=ILU_Lstart(kk0),ILU_Lend(kk0)
             do k = kk0,kk0+kkstep
#else          
          do jj=1,ILU_Lend(kk0)-ILU_Lstart(kk0)+1
             do k = kk0,kk0+kkstep
                j = ILU_Lstart(k) + jj -1
#endif                
                z(k-kk+1)=z(k-kk+1)-s(k+offset(j))*diaA(k+(j-1)*n)
             enddo
          enddo

          kk0 = kk0+kkstep+1
          if(kk0>kkmin)then
             exit
          endif
       enddo

#endif

          ! z = D^-1 (r - Lz) 
          do k=kk,min(kk+independ_nvec-1,n_endth)
             s(k)=Dinv(k)*z(k-kk+1)
          enddo
       enddo
       
       ! -----
       ! D^-1(D+U)s=z
       !   (I+D^-1U)s=z 
       !   s=z - D^-1 Us
       do kk=n_endth,n_startth,-independ_nvec
          do k = kk,max(kk-independ_nvec+1,n_startth),-1
             z(kk-k+1)=0.0d0
          enddo
#if 0
          do k = kk,max(kk-independ_nvec+1,n_startth),-1
             do j=ILU_Ustart(k),ILU_Uend(k),-1
                z(kk-k+1)=z(kk-k+1)+s(k+offset(j))*diaA(k+(j-1)*n)
             enddo
          enddo
#else
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

#if 1
          do j=ILU_Ustart(kk0),ILU_Uend(kk0),-1
             do k = kk0,kk0-kkstep,-1
#else                          
          do jj=1,ILU_Ustart(kk0)-ILU_Uend(kk0)+1
             do k = kk0,kk0-kkstep,-1
                j = ILU_Ustart(k) - jj + 1
#endif                
                z(kk-k+1)=z(kk-k+1)+s(k+offset(j))*diaA(k+(j-1)*n)
             enddo
          enddo

          kk0 = kk0-kkstep-1
          if(kk0<kkmin)then
             exit
          endif
       enddo
#endif
          do k = kk,max(kk-independ_nvec+1,n_startth),-1
             s(k)=s(k)-Dinv(k)*z(kk-k+1)
          enddo
       enddo
    endif

  end subroutine solve_ILU_dia

  subroutine init_solve_ILU_dia( &
       n,nnd,&
       offset, &
       independ_nvec, &
       istart,iend, &
       istart_th,iend_th, &
       n_startth,n_endth,  &
       ILU_Lstart, &
       ILU_Lend,   &
       ILU_Ustart, &
       ILU_Uend    &
       )
    use matrix
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) nnd
    integer(kind=intType) offset(nnd)

    integer(kind=intType) istart,iend
    integer(kind=intType) istart_th,iend_th
    integer(kind=intType) n_startth,n_endth

    integer(kind=intType) j,k,jj,kk
    integer(kind=intType) Row,Col,LCol

    integer(kind=intType) ILU_Lstart(n_startth:n_endth)
    integer(kind=intType) ILU_Lend(n_startth:n_endth)
    integer(kind=intType) ILU_Ustart(n_startth:n_endth)
    integer(kind=intType) ILU_Uend(n_startth:n_endth)
    integer independ_nvec
    integer independ_start,independ_end
    integer irank,ierr
    integer kk0,kkmin,kkstep
    ! -----
    ! (D+L)D^-1(D+U)s=r

    ! -----
    !   z = D^-1(D+U)s
    !   (D+L)z=r
!    call MPI_Comm_rank(MPI_COMM_WORLD,irank, ierr)

    ! kはローカルな行番号
    do k=n_startth,n_endth
       ! 独立に回すループサイズ
       ! 計算範囲の行番号に相当する
       independ_start = n_startth + ((k-n_startth)/independ_nvec) * independ_nvec
       independ_end   = min(independ_start+independ_nvec-1,n_endth)

#if 0
       ILU_Lstart(k)=1
       ILU_Lend(k)=nnd
#else
       ILU_Lstart(k)=nnd
       ILU_Lend(k)  =-1
#endif
       ! istart_th : グローバルな行番号
       Row=independ_start-n_startth+istart_th

!       write(*,*)"independ : ",k+irank*9,Row,independ_start ,independ_end   
       do j=1,nnd
          if(offset(j)<0)then
             ! グローバルな列番号
             Col=istart+offset(j)+k-1
             ! ブロックヤコビ前処理の範囲になるまでループさせる
             if( Col<istart_th ) then 
                cycle
             endif

             LCol = offset(j)+k
!             write(*,*)"GCol,LCol : ",k+irank*9,irank,Col,LCol
             if( independ_start>LCol ) then
                ILU_Lstart(k)=j
             endif
             exit 
          else
             exit
          endif
       enddo
       do j=1,nnd
          if(offset(j)>=0)then
             do jj=j-1,1,-1
                Col  = istart+offset(j)+k-1
                LCol = offset(jj)+k
                if( independ_start>LCol ) then
!                if( (independ_start>LCol).and.(Col>=istart_th) ) then
                   ILU_Lend(k)=jj
                   exit
                endif
             enddo
             exit
          endif
       enddo
    enddo

    ! -----
    ! D^-1(D+U)s=z
    !   (I+D^-1U)s=z 
    !   s=z - D^-1 Us
    do k=n_endth,n_startth,-1
       ! 独立に回すループサイズ
       independ_start = n_startth + ((k-n_startth)/independ_nvec) * independ_nvec
       independ_end   = min(independ_start+independ_nvec-1,n_endth)

!       Row=k-n_startth+istart_th
!       Row=independ_start+istart_th
       Row=independ_start-n_startth+istart_th

       ILU_Ustart(k)= 1
       ILU_Uend(k)  = nnd
       do j=nnd,1,-1
          if(0<offset(j))then
             Col=istart+offset(j)+k-1
             if( iend_th<Col )then
                cycle
             endif
             LCol = offset(j)+k
             if( independ_end<LCol ) then
                ILU_Ustart(k)=j
             endif
             exit
          else
             exit
          endif
       enddo

       do j=nnd,1,-1
          if(offset(j)<=0)then
             do jj = j+1,nnd,1
                LCol = offset(jj)+k
                if( independ_end<LCol ) then
                   ILU_Uend(k) = jj
                   exit
                endif
             enddo
             exit
          endif
       enddo

    enddo

#if 0
    do kk=n_startth,n_endth,independ_nvec
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
          ! !$OMP MASTER
          ! do k = kk0,kk0+kkstep
          !    write(*,*)"init solve ILU dia ",kk,k
          ! enddo
          ! write(*,*)" ---------- "
          ! !$OMP END MASTER
          kk0 = kk0+kkstep+1
          if(kk0>kkmin)then
             exit
          endif
       enddo
    enddo
#endif

#if 0
    do kk=n_endth,n_startth,-independ_nvec
       kk0 = kk
       kkmin = max(kk-independ_nvec+1,n_startth)
       do 
          kkstep = 0
          do k = kk0,kkmin+1,-1 
             if ( ( ILU_Ustart(k).eq.ILU_Ustart(k-1) ).and.( ILU_Uend(k).eq.ILU_Uend(k-1) ) )then
                kkstep=kkstep+1
             else
                exit
             endif
          enddo
          ! !$OMP MASTER
          ! do k = kk0,kk0-kkstep,-1
          !    write(*,*)"init solve ILU dia ",kk,k
          ! enddo
          ! write(*,*)" ---------- "
          ! !$OMP END MASTER
          kk0 = kk0-kkstep-1
          if(kk0<kkmin)then
             exit
          endif

       enddo
    enddo
#endif
       
  end subroutine init_solve_ILU_dia

  subroutine AS_Local(prm,preconAS,mnumth,iflagAS,s,precon_thblock)
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
    real*8 s(prm%n)

    integer(kind=intType) k 
    integer(kind=intType) n
    integer(kind=intType) istart,iend

    integer numth,myth
    integer(kind=intType) istart_th,iend_th
    integer(kind=intType) n_startth,n_endth
    integer max_numth
    integer precon_thblock
    integer n_startth0,n_endth0  

    n=prm%n
    istart=prm%istart
    iend=prm%iend

    if( (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_RESTRICT).or.(iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_NONE) )then
       ! 局所解
       call get_maxthreads(max_numth)
       call set_precon_maxthreads(n,precon_thblock)

       !$OMP PARALLEL &
       !$OMP PRIVATE( &
       !$OMP numth,myth,n_startth,n_endth,istart_th,iend_th,k &
       !$OMP )       
       numth = omp_get_num_threads()
       myth  = omp_get_thread_num()
       n_startth=preconAS(myth+1)%n0_startth
       n_endth  =preconAS(myth+1)%n0_endth
       do k=n_startth,n_endth
          s(k)=preconAS(myth+1)%svec(k)
       enddo
       !$OMP END PARALLEL 
       call omp_set_num_threads(max_numth)
 
    else
       ! 加算補正
       call get_maxthreads(max_numth)
       call set_precon_maxthreads(n,precon_thblock)
       call zero_initialize(n,s)
 
#if 1
       !$OMP PARALLEL &
       !$OMP PRIVATE( &
       !$OMP numth,myth,n_startth,n_endth,istart_th,iend_th,k, &
       !$OMP n_startth0,n_endth0  &
       !$omp ) 
       numth = omp_get_num_threads()
       myth  = omp_get_thread_num()
       n_startth = preconAS(myth+1)%n_startth
       n_endth   = preconAS(myth+1)%n_endth       

       n_startth0 = preconAS(myth+1)%n0_startth
       n_endth0   = preconAS(myth+1)%n0_endth

       do k=n_startth,n_startth0-1
          s(k)=s(k)+preconAS(myth+1)%svec(k)
       enddo
       !$OMP BARRIER
       do k=n_startth0,n_endth0
          s(k)=s(k)+preconAS(myth+1)%svec(k)
       enddo
       !$OMP BARRIER
       do k=n_endth0+1,n_endth
          s(k)=s(k)+preconAS(myth+1)%svec(k)
       enddo

       !$OMP END PARALLEL 
#else
       !$OMP PARALLEL &
       !$OMP PRIVATE( &
       !$OMP numth,myth,n_startth,n_endth,istart_th,iend_th,k &
       !$omp ) &
       !$OMP REDUCTION(+:s) 
       numth = omp_get_num_threads()
       myth  = omp_get_thread_num()
       n_startth = preconAS(myth+1)%n_startth
       n_endth   = preconAS(myth+1)%n_endth       
       do k=n_startth,n_endth
          s(k)=s(k)+preconAS(myth+1)%svec(k)
       enddo
       !$OMP END PARALLEL 
#endif
       call omp_set_num_threads(max_numth)
    endif

    !    stop
  end subroutine AS_Local

  !  subroutine AdditiveSchwarz_ILU_crs(prm,crsMat,preconAS,mnumth,iflagAS,r)
  subroutine AdditiveSchwarz_ILU_crs(method,iflagAS)
    use omp_lib
    use krylov
    use vec_kernel
    use PARCEL_integer
    implicit none
    type(KSP) :: method
    integer iflagAS

    integer(kind=intType) n
    integer(kind=intType) nnz
    integer(kind=intType) istart,iend
    integer(kind=intType) k
    integer numth,myth,max_numth
    integer(kind=intType) n_th
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th
    integer precon_thblock

    n      = method%prm_AS%n
    nnz    = method%prm_AS%nnz
    istart = method%prm_AS%istart
    iend   = method%prm_AS%iend
    precon_thblock=method%precon_thblock 

!    write(*,*)"solve -- independ_nvec = ",method%independ_nvec
    ! -----
    ! ブロック数よりもスレッド数が大きい場合はスレッド数を制限する
    call get_maxthreads(max_numth)
    call set_precon_maxthreads(n,precon_thblock)
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
       call solve_ILU_crs( &
            n,nnz, &
            method%precon_crsMatAS%crsA, &
            method%precon_crsMatAS%crsRow_ptr, &
            method%precon_crsMatAS%crsCol, &
            istart,iend, &
            istart_th,iend_th, &
            n_startth,n_endth, &
            method%preconAS(myth+1)%Dinv, &
            method%ivec_AS(n_startth), & 
            method%preconAS(myth+1)%svec, &
            method%ILUregionMat(myth+1)%ILU_Lstart, &
            method%ILUregionMat(myth+1)%ILU_Lend,   &
            method%ILUregionMat(myth+1)%ILU_Ustart, &
            method%ILUregionMat(myth+1)%ILU_Uend,   &
            method%independ_nvec &
            )
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
          method%preconAS(myth+1)%rvec(k)=method%ivec_AS(k)
       enddo
       istart_th = method%preconAS(myth+1)%istart_th
       iend_th   = method%preconAS(myth+1)%iend_th
       n_startth = method%preconAS(myth+1)%n_startth
       n_endth   = method%preconAS(myth+1)%n_endth

       ! interpolateではオーバーラップ領域の残差をゼロにする必要がある
       call solve_ILU_crs( &
            n,nnz, &
            method%precon_crsMatAS%crsA, &
            method%precon_crsMatAS%crsRow_ptr, &
            method%precon_crsMatAS%crsCol, &
            istart,iend, &
            istart_th,iend_th, &
            n_startth,n_endth, &
            method%preconAS(myth+1)%Dinv, &
            method%preconAS(myth+1)%rvec, &
            method%preconAS(myth+1)%svec, &
            method%ILUregionMat(myth+1)%ILU_Lstart, &
            method%ILUregionMat(myth+1)%ILU_Lend,   &
            method%ILUregionMat(myth+1)%ILU_Ustart, &
            method%ILUregionMat(myth+1)%ILU_Uend,   &
            method%independ_nvec &
            )
       !$OMP END PARALLEL
    endif

!    if(n<=max_numth)then
    call omp_set_num_threads(max_numth)
!    endif
  end subroutine AdditiveSchwarz_ILU_crs

  subroutine init_AdditiveSchwarz_crs(method)
    use omp_lib
    use krylov
    use vec_kernel
    use PARCEL_integer
    implicit none
    type(KSP) :: method
    integer iflagAS

    integer(kind=intType) n
    integer(kind=intType) nnz
    integer(kind=intType) istart,iend
    integer(kind=intType) k
    integer numth,myth,max_numth
    integer(kind=intType) n_th
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th
    integer precon_thblock
    n      = method%prm_AS%n
    nnz    = method%prm_AS%nnz
    istart = method%prm_AS%istart
    iend   = method%prm_AS%iend
    iflagAS = method%iflagAS
    precon_thblock=method%precon_thblock 

    ! -----
    ! ブロック数よりもスレッド数が大きい場合はスレッド数を制限する
    call get_maxthreads(max_numth)
    call set_precon_maxthreads(n,precon_thblock)

    !$OMP PARALLEL 
    !$OMP MASTER
    numth = omp_get_num_threads()
    !$OMP END MASTER
    !$OMP END PARALLEL 
    allocate(method%ILUregionMat(numth))

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

!       write(*,*)myth,"istart_th,iend_th -- ",istart_th,iend_th

       allocate(method%ILUregionMat(myth+1)%ILU_Lstart(n_startth:n_endth))
       allocate(method%ILUregionMat(myth+1)%ILU_Lend(n_startth:n_endth)  )
       allocate(method%ILUregionMat(myth+1)%ILU_Ustart(n_startth:n_endth))
       allocate(method%ILUregionMat(myth+1)%ILU_Uend(n_startth:n_endth)  )
       call init_solve_ILU_crs( &
            n,nnz, &
            method%precon_crsMatAS%crsRow_ptr, &
            method%precon_crsMatAS%crsCol, &
            method%independ_nvec, &
            istart,iend, &
            istart_th,iend_th, &
            n_startth,n_endth, &
            method%ILUregionMat(myth+1)%ILU_Lstart, &
            method%ILUregionMat(myth+1)%ILU_Lend,   &
            method%ILUregionMat(myth+1)%ILU_Ustart, &
            method%ILUregionMat(myth+1)%ILU_Uend    &
            )
       !$OMP END PARALLEL
    elseif( (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_INTERPOLATE).or.(iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_NONE) )then
       !$OMP PARALLEL &
       !$OMP PRIVATE( &
       !$OMP numth,myth,n_startth,n_endth,istart_th,iend_th,k &
       !$OMP ) 

       numth = omp_get_num_threads()
       myth  = omp_get_thread_num()
       istart_th = method%preconAS(myth+1)%istart_th
       iend_th   = method%preconAS(myth+1)%iend_th
       n_startth = method%preconAS(myth+1)%n_startth
       n_endth   = method%preconAS(myth+1)%n_endth

       allocate(method%ILUregionMat(myth+1)%ILU_Lstart(n_startth:n_endth))
       allocate(method%ILUregionMat(myth+1)%ILU_Lend(n_startth:n_endth)  )
       allocate(method%ILUregionMat(myth+1)%ILU_Ustart(n_startth:n_endth))
       allocate(method%ILUregionMat(myth+1)%ILU_Uend(n_startth:n_endth)  )
       call init_solve_ILU_crs( &
            n,nnz, &
            method%precon_crsMatAS%crsRow_ptr, &
            method%precon_crsMatAS%crsCol, &
            method%independ_nvec, &
            istart,iend, &
            istart_th,iend_th, &
            n_startth,n_endth, &
            method%ILUregionMat(myth+1)%ILU_Lstart, &
            method%ILUregionMat(myth+1)%ILU_Lend,   &
            method%ILUregionMat(myth+1)%ILU_Ustart, &
            method%ILUregionMat(myth+1)%ILU_Uend    &
            )
       !$OMP END PARALLEL
    endif

!    if(n<=max_numth)then
    call omp_set_num_threads(max_numth)
!    endif

  end subroutine init_AdditiveSchwarz_crs

  !  subroutine AdditiveSchwarz_ILU_dia(prm,diaMat,preconAS,mnumth,iflagAS,r)
  subroutine AdditiveSchwarz_ILU_dia(method,iflagAS)
    use omp_lib
    use krylov
    use vec_kernel
    use PARCEL_integer
    implicit none
    type(KSP) :: method
    integer iflagAS

    integer(kind=intType) n
    integer(kind=intType) nnd
    integer(kind=intType) istart,iend
    integer(kind=intType) k
    integer numth,myth,max_numth
    integer(kind=intType) n_th
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th
    integer precon_thblock

    n      = method%prm_AS%n
    nnd    = method%prm_AS%nnd
    istart = method%prm_AS%istart
    iend   = method%prm_AS%iend
    precon_thblock=method%precon_thblock 

    ! -----
    ! ブロック数よりもスレッド数が大きい場合はスレッド数を制限する
    call get_maxthreads(max_numth)
    call set_precon_maxthreads(n,precon_thblock)
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
       call solve_ILU_dia( &
            n,nnd,&
            method%precon_diaMatAS%diaA  , &
            method%precon_diaMatAS%offset, &
            istart,iend, &
            istart_th,iend_th, &
            n_startth,n_endth, &
            method%preconAS(myth+1)%Dinv, &
            !            r(n_startth), & 
            method%ivec_AS(n_startth), & 
            method%preconAS(myth+1)%svec, &
            method%ILUregionMat(myth+1)%ILU_Lstart, &
            method%ILUregionMat(myth+1)%ILU_Lend,   &
            method%ILUregionMat(myth+1)%ILU_Ustart, &
            method%ILUregionMat(myth+1)%ILU_Uend,   &
            method%independ_nvec    &
            )

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
          method%preconAS(myth+1)%rvec(k)=method%ivec_AS(k)
       enddo
       istart_th = method%preconAS(myth+1)%istart_th
       iend_th   = method%preconAS(myth+1)%iend_th
       n_startth = method%preconAS(myth+1)%n_startth
       n_endth   = method%preconAS(myth+1)%n_endth
       ! interpolateではオーバーラップ領域の残差をゼロにする必要がある
       call solve_ILU_dia( &
            n,nnd,         &
            method%precon_diaMatAS%diaA  , &
            method%precon_diaMatAS%offset, &
            istart,iend,       &
            istart_th,iend_th, &
            n_startth,n_endth, &
            method%preconAS(myth+1)%Dinv, &
            method%preconAS(myth+1)%rvec, &
            method%preconAS(myth+1)%svec, &
            method%ILUregionMat(myth+1)%ILU_Lstart, &
            method%ILUregionMat(myth+1)%ILU_Lend,   &
            method%ILUregionMat(myth+1)%ILU_Ustart, &
            method%ILUregionMat(myth+1)%ILU_Uend,   &
            method%independ_nvec    &
            )

       !$OMP END PARALLEL
    endif

!    if(n<=max_numth)then
    call omp_set_num_threads(max_numth)
!    endif

  end subroutine AdditiveSchwarz_ILU_dia

  subroutine init_AdditiveSchwarz_dia(method)
    use omp_lib
    use krylov
    use PARCEL_integer
    implicit none
    type(KSP) :: method
    integer iflagAS

    integer(kind=intType) n
    integer(kind=intType) nnd
    integer(kind=intType) istart,iend
    integer(kind=intType) k
    integer numth,myth,max_numth
    integer(kind=intType) n_th
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th
    integer precon_thblock
    n       = method%prm_AS%n
    nnd     = method%prm_AS%nnd
    istart  = method%prm_AS%istart
    iend    = method%prm_AS%iend
    iflagAS = method%iflagAS
    precon_thblock=method%precon_thblock 

    ! -----
    ! ブロック数よりもスレッド数が大きい場合はスレッド数を制限する
    call get_maxthreads(max_numth)
    call set_precon_maxthreads(n,precon_thblock)

    !$OMP PARALLEL 
    !$OMP MASTER
    numth = omp_get_num_threads()
    !$OMP END MASTER
    !$OMP END PARALLEL 
    allocate(method%ILUregionMat(numth))

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

       allocate(method%ILUregionMat(myth+1)%ILU_Lstart(n_startth:n_endth))
       allocate(method%ILUregionMat(myth+1)%ILU_Lend(n_startth:n_endth)  )
       allocate(method%ILUregionMat(myth+1)%ILU_Ustart(n_startth:n_endth))
       allocate(method%ILUregionMat(myth+1)%ILU_Uend(n_startth:n_endth)  )
       call init_solve_ILU_dia( &
            n,nnd,&
            method%precon_diaMatAS%offset, &
            method%independ_nvec, &
            istart,iend, &
            istart_th,iend_th, &
            n_startth,n_endth, &
            method%ILUregionMat(myth+1)%ILU_Lstart, &
            method%ILUregionMat(myth+1)%ILU_Lend,   &
            method%ILUregionMat(myth+1)%ILU_Ustart, &
            method%ILUregionMat(myth+1)%ILU_Uend    &
            )
       !$OMP END PARALLEL
    elseif( (iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_INTERPOLATE).or.(iflagAS.eq.KSP_ADDITIVESCHWARZ_TYPE_NONE) )then
       !$OMP PARALLEL &
       !$OMP PRIVATE( &
       !$OMP numth,myth,n_startth,n_endth,istart_th,iend_th,k &
       !$OMP ) 
       numth = omp_get_num_threads()
       myth  = omp_get_thread_num()

       istart_th = method%preconAS(myth+1)%istart_th
       iend_th   = method%preconAS(myth+1)%iend_th
       n_startth = method%preconAS(myth+1)%n_startth
       n_endth   = method%preconAS(myth+1)%n_endth

       allocate(method%ILUregionMat(myth+1)%ILU_Lstart(n_startth:n_endth))
       allocate(method%ILUregionMat(myth+1)%ILU_Lend(n_startth:n_endth)  )
       allocate(method%ILUregionMat(myth+1)%ILU_Ustart(n_startth:n_endth))
       allocate(method%ILUregionMat(myth+1)%ILU_Uend(n_startth:n_endth)  )
       call init_solve_ILU_dia( &
            n,nnd,         &
            method%precon_diaMatAS%offset, &
            method%independ_nvec, &
            istart,iend,       &
            istart_th,iend_th, &
            n_startth,n_endth, &
            method%ILUregionMat(myth+1)%ILU_Lstart, &
            method%ILUregionMat(myth+1)%ILU_Lend,   &
            method%ILUregionMat(myth+1)%ILU_Ustart, &
            method%ILUregionMat(myth+1)%ILU_Uend    &
            )

       !$OMP END PARALLEL
    endif

!    if(n<=max_numth)then
    call omp_set_num_threads(max_numth)
!    endif

  end subroutine Init_AdditiveSchwarz_dia

  !  subroutine view_crs(prm,crsMat)
  subroutine view_crs(n,istart,crsMat)
    use matrix
    use PARCEL_integer
    implicit none
    !    type(prmMat)  :: prm 
    type(CRS_Mat) :: crsMat 
    integer(kind=intType) i,j,k
    integer(kind=intType) n,istart
    !   n = prm%n
    !   istart=prm%istart
    do k=1,n
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          write(*,*)"  CRS  ",k-1+istart,crsMat%crsCol(j),crsMat%crsA(j)
       enddo
    enddo
  end subroutine view_crs

  !  subroutine view_dia(prm,diaMat)
  subroutine view_dia(n,nnd,istart,diaMat)
    use matrix
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(DIA_Mat) :: diaMat 
    integer(kind=intType) i,j,k
    integer(kind=intType) n,nnd,istart
    integer(kind=intType) Row,Col
    !     n = prm%n
    !     nnd = prm%nnd
    !     istart=prm%istart
    do k=1,n
       Row=k-1+istart
       do j=1,nnd
          Col=Row+diaMat%offset(j)
          write(*,*)"  DIA  ",Row,Col,diaMat%diaA(k+(j-1)*n)
       enddo
    enddo
  end subroutine view_dia

  ! DIAの疎行列ベクトル積のベクトル位置とベクトルサイズを計算
  subroutine set_dia_Row(prm,diaMat)
    use matrix
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
    call replace_dia_vec(&
         diaMat,&
         istart,&
         RowL,RowU,&
         n,nColL,nColD,nColU,nCol)    
  end subroutine set_dia_Row

  subroutine get_dia_Row(&
       prm,diaMat,&
       RowL,RowU,&
       n,nColL,nColD,nColU,nCol)    
    use matrix
    use vec_kernel
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(DIA_Mat) :: diaMat 
    integer(kind=intType) RowL(prm%nnd*prm%n)
    integer(kind=intType) RowU(prm%nnd*prm%n)
    integer(kind=intType) nCol,nColL,nColD,nColU

    integer(kind=intType) i,j,k
    integer(kind=intType) n,gn,nnd,istart,iend
    integer(kind=intType) Row
    integer(kind=intType) nL,nU

    n = prm%n
    gn = prm%gn
    nnd = prm%nnd
    istart=prm%istart
    iend=prm%iend

    nCol=0
    nColL=0
    nColD=n
    nColU=0
    RowL=0
    RowU=0

    if(nnd.eq.0)then
       return
    endif

    nL=0
    nU=0
    do j=1,nnd
       if(j.eq.1)then
          if( diaMat%offset(j) < 0 )then
             diaMat%offsetVec(j)=0
             if( (diaMat%offset(j)+n) < 0)then
                !対角ブロックに重ならない場合
                nColL=nColL+n

                Row=istart+diaMat%offset(j)
                do i=1,n
                   nL=nl+1
                   RowL(nL)=Row + i -1
                enddo
             else
                nColL=nColL+abs(diaMat%offset(j))

                Row=istart+diaMat%offset(j) 
                do i=1,abs(diaMat%offset(j))
                   nL=nL+1
                   RowL(nL)=Row + i -1
                enddo

             endif
          endif
          if( 0<diaMat%offset(j) )then
             diaMat%offsetVec(j)=n
             if(n<diaMat%offset(j))then
                !対角ブロックに重ならない場合
                nColU=nColU+n

                Row=istart+diaMat%offset(j)
                do i=1,n
                   nU=nU+1
                   RowU(nU)=Row + i -1
                enddo
             else
                nColU=nColU+diaMat%offset(j)

                Row=istart+diaMat%offset(j)
                do i=n - diaMat%offset(j)+1,n
                   nU=nU+1
                   RowU(nU)=Row + i -1
                enddo
             endif
          endif
       else
          if( (diaMat%offset(j)-diaMat%offset(j-1)) < n )then
             !オーバーラップする場合
             diaMat%offsetVec(j) = &
                  diaMat%offsetVec(j-1)+diaMat%offset(j)-diaMat%offset(j-1)

             ! オフセットが0未満の場合
             if(diaMat%offset(j)<0)then
                if( (diaMat%offset(j)+n) < 0)then
                   !対角ブロックに重ならない場合
                   nColL=nColL + diaMat%offset(j)-diaMat%offset(j-1) 

                   Row=istart+diaMat%offset(j) 
                   do i = n - ( diaMat%offset(j)-diaMat%offset(j-1) ) + 1,n
                      nL=nL+1
                      RowL(nL)=Row + i -1
                   enddo
                else
                   !対角ブロックに重なる場合
                   ! 一つ前の対角要素が対角ブロックに重ならない場合
                   if( (diaMat%offset(j-1) + n ) < 0 )then
                      nColL=nColL + &
                           abs( diaMat%offset(j) )  - &
                           ( n - (diaMat%offset(j)-diaMat%offset(j-1)) )

                      Row=istart+diaMat%offset(j) 
                      do i = n - ( diaMat%offset(j)-diaMat%offset(j-1) ) + 1 ,& 
                           abs(diaMat%offset(j))
                         nL=nL+1
                         RowL(nL)=Row + i -1
                      enddo
                   endif
                endif
             endif

             ! オフセットが0を超える場合
             if(0<diaMat%offset(j))then
                if( n <= (diaMat%offset(j)) )then
                   ! 対角ブロックに重ならない場合
                   nColU=nColU +  diaMat%offset(j)-diaMat%offset(j-1) 

                   Row=istart+diaMat%offset(j)
                   do i=n - ( diaMat%offset(j)-diaMat%offset(j-1) ) + 1,n
                      nU=nU+1
                      RowU(nU)=Row + i -1
                   enddo

                else
                   !対角ブロックに重なる場合
                   if(diaMat%offset(j-1)<n)then
                      ! 一つ前の対角要素が対角ブロックに重なる場合
                      nColU = nColU + &
                           (diaMat%offset(j)-diaMat%offset(j-1) )

                      Row=istart+diaMat%offset(j)
                      do i=n - ( diaMat%offset(j)-diaMat%offset(j-1) ) + 1,n
                         nU=nU+1
                         RowU(nU)=Row + i -1
                      enddo
                   else
                      ! 一つ前の対角要素が対角ブロックに重ならない場合
                      nColU = nColU + diaMat%offset(j)-diaMat%offset(j-1)

                      Row=istart+diaMat%offset(j)
                      do i=n - ( diaMat%offset(j)-diaMat%offset(j-1) ) + 1,n
                         nU=nU+1
                         RowU(nU)=Row + i -1
                      enddo
                   endif

                endif
             endif

          else
             ! オーバーラップしない場合

             if( diaMat%offset(j) < 0 )then
                diaMat%offsetVec(j)=diaMat%offsetVec(j-1)+n
                if( (diaMat%offset(j)+n) < 0)then
                   !対角ブロックに重ならない場合
                   nColL=nColL+n

                   Row=istart+diaMat%offset(j)
                   do i=1,n
                      nL=nl+1
                      RowL(nL)=Row + i -1
                   enddo
                else
                   nColL=nColL+abs(diaMat%offset(j))

                   Row=istart+diaMat%offset(j) 
                   do i=1,abs(diaMat%offset(j))
                      nL=nL+1
                      RowL(nL)=Row + i -1
                   enddo
                endif
             endif
             if( 0 < diaMat%offset(j) )then
                diaMat%offsetVec(j)=diaMat%offsetVec(j-1)+n
                if(n<=diaMat%offset(j))then
                   diaMat%offsetVec(j)=diaMat%offsetVec(j-1)+n
                   !対角ブロックに重ならない場合
                   nColU=nColU+n

                   Row=istart+diaMat%offset(j)
                   do i=1,n
                      nU=nU+1
                      RowU(nU)=Row + i -1
                   enddo
                else

                   nColU=nColU+diaMat%offset(j)

                   Row=istart+diaMat%offset(j)
                   do i = n - diaMat%offset(j) + 1,n
                      nU=nU+1
                      RowU(nU)=Row + i -1 
                   enddo
                endif

             endif
          endif

       endif

    enddo
    if( (iend.eq.gn) )then
       if(0<diaMat%offset(nnd))then
          nColU=diaMat%offset(nnd)

          Row=istart+diaMat%offset(nnd)
          do i=n - diaMat%offset(nnd),n
             nU=nU+1
             RowU(nU)=Row + i -1
          enddo
       endif
    endif
    ! -----
    ! 対角ブロックによる補正
    do j=1,nnd
       if( 0.eq.diaMat%offset(j) )then
          ! 対角要素があるので補正無し
          exit
       endif
       if( 0<diaMat%offset(j) )then

          if(j.eq.1)then
             ! マイナスオフセットの非零要素が無い
             if(n<diaMat%offset(j))then
                ! 対角ブロックに非零要素がない
                do i=j,nnd
                   diaMat%offsetVec(i) = diaMat%offsetVec(i)+n
                enddo
                exit
             else
                do i=j,nnd
                   diaMat%offsetVec(i) = diaMat%offsetVec(i) + &
                        diaMat%offsetVec(j)
                enddo
             endif
          else
             ! マイナスオフセットの非零要素がある
             if(diaMat%offset(j-1) <= -n)then
                ! 対角ブロックにマイナスオフセットの対角要素が含まれない
                if(n<diaMat%offset(j))then
                   ! 対角ブロックに非零要素がない
                   do i=j,nnd
                      diaMat%offsetVec(i) = diaMat%offsetVec(i)+n
                   enddo
                   exit
                else
                   do i=j,nnd
                      diaMat%offsetVec(i) = diaMat%offsetVec(i) + &
                           diaMat%offset(j)
                      ! ------
                      ! diaMat%offsetVec(i) = diaMat%offsetVec(i) 

                   enddo
                   exit
                endif
             else
                ! 対角ブロックにマイナスオフセットの対角要素が含まれる
                if(n<diaMat%offset(j))then
                   ! 対角ブロックに非零要素が含まれない
                   do i=j,nnd
                      diaMat%offsetVec(i) = diaMat%offsetVec(i)  + &
                           n - (n - abs(diaMat%offset(j-1)) ) 
                   enddo
                   exit
                else
                   ! 対角ブロックに非零要素が含まれる
                   do i=j,nnd
                      diaMat%offsetVec(i) = diaMat%offsetVec(i) + &
                           diaMat%offset(j) - (diaMat%offset(j-1)+n) 
                   enddo
                   exit
                endif

             endif

          endif

       endif
    enddo

    ! -----
    ! ベクトル全体の長さ 
    ! nCol=nColL+nColD+nColU

  end subroutine get_dia_Row

  subroutine replace_dia_vec(&
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
    deallocate(diaMat%Vec)

    allocate(diaMat%RowL(nColL))
    allocate(diaMat%RowD(nColD))
    allocate(diaMat%RowU(nColU))
    allocate(diaMat%Vec(nCol))

    do i=1,nColL
       diaMat%RowL(i)=RowL(i)
    enddo
    do i=1,n
       diaMat%RowD(i)=i+istart-1
    enddo
    do i=1,nColU
       diaMat%RowU(i)=RowU(i)
    enddo
    call zero_initialize(nCol,diaMat%Vec)

  end subroutine replace_dia_vec

  subroutine view_den(prm,denMat)
    use matrix
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(DEN_Mat) :: denMat 
    integer(kind=intType) i,j,k
    integer(kind=intType) n,gn,istart
    n = prm%n
    gn = prm%gn
    istart=prm%istart
    do k=1,n
       do j=1,gn
          write(*,*)"  DEN  ",k-1+istart,j,denMat%denA(j+(k-1)*gn)
       enddo
    enddo
  end subroutine view_den

  subroutine view_nrk(commVec)
    use comm
    use PARCEL_integer
    implicit none
    type(comm_Vec) :: commVec
    integer(kind=intType) i,j

    do i=1,commVec%nrk_rank
       do j=1,commVec%nrk_prm(i)%ndata
          write(*,*) "  nrk  ",&
               commVec%myrank,&
               commVec%nrk_prm(i)%rank, &
               commVec%nrk_prm(i)%offset(j)
       enddo
    enddo

  end subroutine view_nrk

  subroutine view_Hybcrs(prm,HybcrsMat)
    use matrix
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(HYBCRS_Mat) :: HybcrsMat 
    integer(kind=intType) i,j,k
    integer(kind=intType) n,istart

    n = prm%n
    istart=prm%istart

    write(*,*)" Hybrid CRS L"
    do k=1,n
       !       write(*,*)k,HybcrsMat%crsRow_ptr_L(k),HybcrsMat%crsRow_ptr_L(k+1)-1
       do j=HybcrsMat%crsRow_ptr_L(k),HybcrsMat%crsRow_ptr_L(k+1)-1
          !          write(*,*)k-1+istart,HybcrsMat%crsColL(j),HybcrsMat%crsAL(j)
          write(*,*)k-1+istart, &
               HybcrsMat%RowL(HybcrsMat%crsColL(j)), &
               HybcrsMat%crsAL(j)
       enddo
    enddo

    write(*,*)" Hybrid CRS D"
    do k=1,n
       do j=HybcrsMat%crsRow_ptr_D(k),HybcrsMat%crsRow_ptr_D(k+1)-1
          !          write(*,*)k-1+istart,HybcrsMat%crsColD(j),HybcrsMat%crsAD(j)
          write(*,*)k-1+istart, &
               HybcrsMat%RowD(HybcrsMat%crsColD(j)), &
               HybcrsMat%crsAD(j)
       enddo
    enddo

    write(*,*)" Hybrid CRS U"
    do k=1,n
       do j=HybcrsMat%crsRow_ptr_U(k),HybcrsMat%crsRow_ptr_U(k+1)-1
          !          write(*,*)k-1+istart,HybcrsMat%crsColU(j),HybcrsMat%crsAU(j)
          write(*,*)k-1+istart, &
               HybcrsMat%RowU(HybcrsMat%crsColU(j)), &
               HybcrsMat%crsAU(j)
       enddo
    enddo

  end subroutine view_Hybcrs

  subroutine view_coo(prm,cooMat)
    use matrix
    use PARCEL_integer
    implicit none
    type(prmMat)  :: prm 
    type(COO_Mat) :: cooMat 
    integer(kind=intType) i,j,k
    integer(kind=intType) nnz,istart
    nnz = prm%nnz
    istart=prm%istart
    do i=1,nnz
       write(*,*)"  COO  ",cooMat%cooRow(i)-1+istart, &
            cooMat%cooCol(i), &
            cooMat%cooA(i)
    enddo
  end subroutine view_coo

  subroutine den2coo(denMat,cooMat)
    use matrix
    use PARCEL_integer
    implicit none
    !    type(prmMat)  :: prm
    type(DEN_Mat) :: denMat 
    type(COO_Mat) :: cooMat
    real*8 absAmat
    integer(kind=intType) gn,n
    integer(kind=intType) istart
    integer(kind=intType) i,j,k
    n =denMat%RowN
    gn=denMat%ColN
    !    istart=prm%istart
    k=0
    do i=1,n
       do j=1,gn
          absAmat=abs(denMat%denA(j+(i-1)*gn))
          if(absAmat> minReal8)then
             k=k+1
             cooMat%cooA(k)   = denMat%denA(j+(i-1)*gn)
             !              cooMat%cooRow(k) = i + istart -1 
             cooMat%cooRow(k) = i 
             cooMat%cooCol(k) = j
          endif
       enddo
    enddo
  end subroutine den2coo

  ! -----
  ! 疎行列格納形式

  subroutine den2crs(denMat,crsMat)
    use matrix
    use PARCEL_integer
    implicit none
    type(DEN_Mat) :: denMat 
    type(CRS_Mat) :: crsMat
    real*8 absAmat
    integer(kind=intType) gn,n
    integer(kind=intType) i,j,k
    n =denMat%RowN
    gn=denMat%ColN 
    k=1
    do i=1,n
       crsMat%crsRow_ptr(i) = k
       do j=1,gn
          absAmat=abs(denMat%denA(j+(i-1)*gn))
          if(absAmat> minReal8)then
             crsMat%crsA(k)   = denMat%denA(j+(i-1)*gn)
             crsMat%crsCol(k) = j
             k=k+1
          endif
       enddo
    enddo
    crsMat%crsRow_ptr(n+1) = k
  end subroutine den2crs

  subroutine coo2crs(cooMat,crsMat,n,nnz,istart)
    use matrix
    use PARCEL_integer
    implicit none
    type(COO_Mat) :: cooMat 
    type(CRS_Mat) :: crsMat
    integer(kind=intType) n,nnz
    integer(kind=intType) i,j,k
    real*8 valA
    integer(kind=intType) Row,Col
    integer(kind=intType) Row0
    integer(kind=intType) nptr
    integer(kind=intType) istart

    ! coo形式は行の昇順にソートされていること
    k=1
    Row0=cooMat%cooRow(1)
    crsMat%crsRow_ptr(1)=1
    nptr=1
    do i=istart,Row0-1
       k=k+1
       crsMat%crsRow_ptr(k)=nptr
    enddo
    do i=1,nnz
       Row  = cooMat%cooRow(i)
       Col  = cooMat%cooCol(i)
       valA = cooMat%cooA(i)

       if( (Row-Row0).ne.0 )then
          do j=1,Row-Row0
             k=k+1
             crsMat%crsRow_ptr(k)=nptr
          enddo
          Row0=Row
       endif
       crsMat%crsA(nptr)=valA
       crsMat%crsCol(nptr)=Col
       !       write(*,*)nptr,Row,Col ,valA
       nptr=nptr+1
    enddo
    do i=k+1,n
       crsMat%crsRow_ptr(i)=nptr
    enddo
    crsMat%crsRow_ptr(n+1)=nptr    
#if 0
    do i=1,n
       !          write(*,*)"  @@",i,crsMat%crsRow_ptr(i),crsMat%crsRow_ptr(i+1)-1 
       ! do j=crsMat%crsRow_ptr(i),crsMat%crsRow_ptr(i+1)-1
       !    write(*,*)j,i,crsMat%crsCol(j),crsMat%crsA(j)
       ! enddo

    enddo
#endif

  end subroutine coo2crs

  subroutine get_RowLUD( &
       n,nCol,nColL,nColU,&
       Colsort, &
       RowL,RowU,RowD, &
       istart,iend &
       )
    use PARCEL_integer
    implicit none
    integer(kind=intType) n
    integer(kind=intType) nCol,nColL,nColU
    integer(kind=intType) Colsort(nCol)
    integer(kind=intType) RowL(nColL)
    integer(kind=intType) RowD(n)
    integer(kind=intType) RowU(nColU)
    integer(kind=intType) istart
    integer(kind=intType) iend
    integer(kind=intType) Col
    integer(kind=intType) i,j,k,l
    k=0
    j=0
    l=0
    do i=1,nCol
       Col=Colsort(i)
       if(Col<istart)then
          j=j+1
          RowL(j)=Col
       endif
       if(iend<Col)then
          k=k+1
          RowU(k)=Col
       endif
    enddo
    do k=1,n
       RowD(k)=k+istart-1
    enddo

  end subroutine get_RowLUD

  subroutine crs2Hybcrs(crsMat,HybcrsMat,n,nnz,istart,iend,Colsort)
    use matrix
    use PARCEL_integer
    implicit none
    type(CRS_Mat)  ::  crsMat
    type(HybCRS_Mat) :: HybcrsMat

    integer(kind=intType) n
    integer(kind=intType) nnz
    integer(kind=intType) istart,iend
    integer(kind=intType) Colsort(nnz)

    integer(kind=intType) i,j,k,l
    integer(kind=intType) Col,Row
    integer(kind=intType) nCol,nColL,nColD,nColU

    integer(kind=intType) nd,nl,nu

    nColL = HybcrsMat%nColL
    nColD = HybcrsMat%nColD
    nColU = HybcrsMat%nColU
    nCol  = nColL+nColD+nColU

    call get_RowLUD( &
         n,nCol,nColL,nColU,&
         Colsort, &
         HybcrsMat%RowL, &
         HybcrsMat%RowU, &
         HybcrsMat%RowD, &
         istart,iend &
         )

    nl=1
    nd=1
    nu=1
    do k=1,n
       HybcrsMat%crsRow_ptr_L(k)=nl
       HybcrsMat%crsRow_ptr_D(k)=nd
       HybcrsMat%crsRow_ptr_U(k)=nu
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          Row=k-1+istart 
          Col=crsMat%crsCol(j)  
          if( (istart<=Row) .and. (Row<=iend) )then
             if( (istart<=Col) )then
                if(  (Col<=iend) )then
                   HybcrsMat%crsColD(nd)=Col-istart+1
                   HybcrsMat%crsAD(nd)=crsMat%crsA(j)
                   nd=nd+1
                else
                   ! nnzU=nnzU+1
                   HybcrsMat%crsColU(nu)=-1
                   do i=1,nColU
                      if(HybcrsMat%RowU(i).eq.Col)then
                         HybcrsMat%crsColU(nu)=i
                         exit
                      endif
                   enddo
                   HybcrsMat%crsAU(nu)=crsMat%crsA(j)
                   nu=nu+1 
                endif
             else
                ! nnzL=nnzL+1
                HybcrsMat%crsColL(nl)=-1
                do i=1,nColL
                   if(HybcrsMat%RowL(i).eq.Col)then
                      HybcrsMat%crsColL(nl)=i
                      exit
                   endif
                enddo
                HybcrsMat%crsAL(nl)=crsMat%crsA(j)
                nl=nl+1 
             endif
          endif
       enddo
    enddo
    HybcrsMat%crsRow_ptr_L(n+1)=nl
    HybcrsMat%crsRow_ptr_D(n+1)=nd
    HybcrsMat%crsRow_ptr_U(n+1)=nu

  end subroutine crs2Hybcrs

  subroutine get_nnd_crs(crsMat,n,istart,nnz,nnd)
    use sort_kernel
    use matrix
    use PARCEL_integer
    implicit none
    type(CRS_Mat) :: crsMat
    integer(kind=intType) n
    integer(kind=intType) istart
    integer(kind=intType) nnz
    integer(kind=intType) nnd
    integer(kind=intType) offset_sort(nnz)

    integer(kind=intType) i,j,k
    integer(kind=intType) Row,Col
    integer(kind=intType) offset

#if 0
    ! ----
    ! CRSの行列がどの対角列に存在するかは 
    ! 列 - 行 による対角オフセットから判断する
    i=1
    do k=1,n
       Row=k-1+istart
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          Col=crsMat%crsCol(j)
          offset_sort(i)=Col-Row
          i=i+1
       enddo
    enddo
#else
    ! ----
    ! CRSの行列がどの対角列に存在するかは 
    ! 列 - 行 による対角オフセットから判断する

    !$OMP PARALLEL DO PRIVATE(k,Row,j,Col) 
    do k=1,n
       Row=k-1+istart
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          Col=crsMat%crsCol(j)
          offset_sort(j)=Col-Row
       enddo
    enddo
    !$OMP END PARALLEL DO

#endif

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

  end subroutine get_nnd_crs

  subroutine get_nnd_crs_without_sort(&
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

    call init_diaMat(diaMat,n,nnd)
    nnd_p=0
    do i=imn,imx
      if(flag(i)==1) then
        nnd_p=nnd_p+1
        diaMat%offset(nnd_p)=i
      endif
    enddo
    deallocate(Flag)

  end subroutine get_nnd_crs_without_sort

  subroutine crs2dia(crsMat,diaMat,n,istart,nnd,nnz)
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

#if 0
    i=1
    do k=1,n
       Row=k-1+istart
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          Col=crsMat%crsCol(j)
          offset_sort(i)=Col-Row
          i=i+1
       enddo
    enddo
#else

    !$OMP PARALLEL DO PRIVATE(k,Row,j,Col) 
    do k=1,n
       Row=k-1+istart
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          Col=crsMat%crsCol(j)
          offset_sort(j)=Col-Row
       enddo
    enddo
    !$OMP END PARALLEL DO

#endif

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
#if 0
    do k=1,n
       Row=k-1+istart
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          Col=crsMat%crsCol(j)
          offset=Col-Row
          do i=1,nnd
             if(offset.eq.diaMat%offset(i))then
                diaMat%diaA(k+(i-1)*n)=crsMat%crsA(j)
                exit
             endif
          enddo
       enddo
    enddo
#else

    !$OMP PARALLEL DO PRIVATE(k,Row,j,Col,offset,i) 
    do k=1,n
       Row=k-1+istart
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          Col=crsMat%crsCol(j)
          offset=Col-Row
          do i=1,nnd
             if(offset.eq.diaMat%offset(i))then
                diaMat%diaA(k+(i-1)*n)=crsMat%crsA(j)
                exit
             endif
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

#endif

  end subroutine crs2dia

  subroutine crs2dia_without_sort(&
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
      diaMat%diaA(n*(off-1)+i) = crsMat%crsA(j)
!      write(*,*)n*(off-1)+i,i,j,diaMat%diaA(n*(off-1)+i) 
    enddo
  enddo
    
  end subroutine crs2dia_without_sort

  subroutine dia2crs(diaMat,crsMat,n,istart,nnd,nnz)
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
    real*8 absAmat
    i=1
    do k=1,n
       Row=k-1+istart
       crsMat%crsRow_ptr(k) = i
       do j=1,nnd
          Col=Row+diaMat%offset(j)
#if 1
          absAmat=abs(diaMat%diaA(k+(j-1)*n))
          if(absAmat > minReal8)then
             crsMat%crsA(i)=diaMat%diaA(k+(j-1)*n)
             crsMat%crsCol(i) = Col
             i=i+1
          endif
#else
          crsMat%crsA(i)=diaMat%diaA(k+(j-1)*n)
          crsMat%crsCol(i) = Col
          i=i+1
#endif
       enddo
    enddo
    crsMat%crsRow_ptr(n+1) = i
  end subroutine dia2crs

  subroutine get_nnz_dia(diaA,offset,n,istart,nnd,nnz)
    use matrix
    use PARCEL_integer
    implicit none
    type(DIA_Mat) :: diaMat 
    integer(kind=intType) n
    integer(kind=intType) istart
    integer(kind=intType) nnd
    integer(kind=intType) nnz

    real*8 diaA(nnd*n)
    integer(kind=intType) offset(nnd)

    integer(kind=intType) j,k
    integer(kind=intType) Row,Col
    real*8 absAmat

    nnz=0
    do k=1,n
       Row=k-1+istart
       do j=1,nnd
          Col=Row+offset(j)
          absAmat=abs(diaA(k+(j-1)*n))
          if(absAmat > minReal8)then
             nnz=nnz+1
             !             write(*,*)"  DIA  ",Row,Col,diaA(k+(j-1)*n)
          endif
       enddo
    enddo
  end subroutine get_nnz_dia

  subroutine get_nnz_Colsort_dia( &
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
    real*8 absAmat

    n = prm%n
    istart = prm%istart
    iend = prm%iend
    nnd = prm%nnd
    nnz = prm%nnz
    i=0
    do k=1,n
       Row=k-1+istart
       do j=1,nnd
          absAmat=abs(diaMat%diaA(k+(j-1)*n))
          if(absAmat > minReal8)then
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

  end subroutine get_nnz_Colsort_dia

  subroutine free_denMat(denMat)
    use matrix
    use PARCEL_integer
    implicit none
    type(DEN_Mat) :: denMat 
    deallocate(denMat%denA)
  end subroutine free_denMat

  subroutine free_cooMat(cooMat)
    use matrix
    use PARCEL_integer
    implicit none
    type(COO_Mat)  :: cooMat
    deallocate(cooMat%cooA)
    deallocate(cooMat%cooRow)
    deallocate(cooMat%cooCol)
  end subroutine free_cooMat

  subroutine free_diaMat(diaMat)
    use matrix
    use PARCEL_integer
    implicit none
    type(DIA_Mat)  :: diaMat
    deallocate(diaMat%diaA)
    deallocate(diaMat%offset)
    deallocate(diaMat%offsetVec)
    deallocate(diaMat%Vec)

    deallocate(diaMat%RowL)
    deallocate(diaMat%RowD)
    deallocate(diaMat%RowU) 

  end subroutine free_diaMat

  subroutine free_crsMat(crsMat)
    use matrix
    use PARCEL_integer
    implicit none
    type(CRS_Mat)  :: crsMat
    deallocate(crsMat%crsA)
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
  end subroutine free_crsMat

  subroutine free_HybcrsMat(HybcrsMat)
    use matrix
    use PARCEL_integer
    implicit none
    type(HybCRS_Mat)  :: HybcrsMat

    deallocate(HybcrsMat%crsAL)
    deallocate(HybcrsMat%crsAD)
    deallocate(HybcrsMat%crsAU)
    deallocate(HybcrsMat%crsRow_ptr_L)
    deallocate(HybcrsMat%crsRow_ptr_D)
    deallocate(HybcrsMat%crsRow_ptr_U)
    deallocate(HybcrsMat%crsColL)
    deallocate(HybcrsMat%crsColD)
    deallocate(HybcrsMat%crsColU)
    deallocate(HybcrsMat%RowL)
    deallocate(HybcrsMat%RowD)
    deallocate(HybcrsMat%RowU)

  end subroutine Free_HybcrsMat

  subroutine divide_ILU_crs( &
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
    real*8  matA

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
    call init_crsMat(crsMatL,n,nnzL)
    call init_crsMat(crsMatU,n,nnzU)

    nnzL=1
    nnzU=1
    do k=1,n
       crsMatL%crsRow_ptr(k)=nnzL
       crsMatU%crsRow_ptr(k)=nnzU
       Row=k-1+istart
       do j=crsMat%crsRow_ptr(k),crsMat%crsRow_ptr(k+1)-1
          Col =crsMat%crsCol(j)
          matA=crsMat%crsA(j)
          if(crsMat%crsCol(j)<Row)then
             crsMatL%crsA(nnzL)   = matA
             crsMatL%crsCol(nnzL) = Col
             nnzL=nnzL+1
          endif
          if(Row<crsMat%crsCol(j))then
             crsMatU%crsA(nnzU)   = matA
             crsMatU%crsCol(nnzU) = Col
             nnzU=nnzU+1
          endif
          if(Row.eq.crsMat%crsCol(j))then
             crsMatL%crsA(nnzL)   = 1.0d0
             crsMatL%crsCol(nnzL) = Col
             nnzL=nnzL+1

             crsMatU%crsA(nnzU)   = matA 
             crsMatU%crsCol(nnzU) = Col
             nnzU=nnzU+1
          endif
       enddo
    enddo
    crsMatL%crsRow_ptr(n+1)=nnzL
    crsMatU%crsRow_ptr(n+1)=nnzU

  end subroutine divide_ILU_crs

  subroutine divide_ILU_dia( &
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
    call init_diaMat(diaMatL,n,nndL) 
    call init_diaMat(diaMatU,n,nndU) 

    nndL=1
    nndU=1
    do i=1,nnd
       offset=diaMat%offset(i)
       if(offset<0)then
          do j=1,n
             diaMatL%diaA(j+(nndL-1)*n)=diaMat%diaA(j+(i-1)*n)
          enddo
          diaMatL%offset(nndL)=offset
          nndL=nndL+1
       endif
       if(offset.eq.0)then
          do j=1,n
             diaMatU%diaA(j+(nndU-1)*n)=diaMat%diaA(j+(i-1)*n)
             diaMatL%diaA(j+(nndL-1)*n)=1.0d0
          enddo
          diaMatL%offset(nndL)=offset
          nndL=nndL+1
          diaMatU%offset(nndU)=offset
          nndU=nndU+1
       endif
       if(0<offset)then
          do j=1,n
             diaMatU%diaA(j+(nndU-1)*n)=diaMat%diaA(j+(i-1)*n)
          enddo
          diaMatU%offset(nndU)=offset
          nndU=nndU+1
       endif
    enddo

  end subroutine divide_ILU_dia

  subroutine ILU_LU_crs(&
       n,istart,nnz,&
       crsA,crsRow_ptr,crsCol,&
       crsMatL,crsMatU)
    use krylov
    use PARCEL_integer
    implicit none

    integer(kind=intType) n,nnz,istart
    real*8 crsA(nnz)
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
    integer,allocatable::D_ptr(:)

    call init_crsMat(crsMatLU,n,nnz) 
    allocate(D_ptr(n))
    do i=1,nnz
       crsMatLU%crsA(i)   = crsA(i)
       crsMatLU%crsCol(i) = crsCol(i)
    enddo
    do i=1,n+1
       crsMatLU%crsRow_ptr(i)=crsRow_ptr(i)
    enddo

    kk=0
    do i=1,n
      do kCol=crsMatLU%crsRow_ptr(i),crsMatLU%crsRow_ptr(i+1)-1
         k=crsMatLU%crsCol(kCol)-istart+1
         if(k==i) then
           D_ptr(i)=kCol
           kk=kk+1
           exit
         endif
      enddo
    enddo

    if(kk/=n) then
      write(*,*)"Error::ILU--found zero in Diagonal element."
      stop
    endif

    do i=2,n
       do kCol=crsMatLU%crsRow_ptr(i),crsMatLU%crsRow_ptr(i+1)-1 
          k=crsMatLU%crsCol(kCol)-istart+1
          if(k>i-1)then
             exit
          endif

          kk=D_ptr(k)
          do j=crsMatLU%crsRow_ptr(k),crsMatLU%crsRow_ptr(k+1)-1
             Col=crsMatLU%crsCol(j) - istart+1
             if( Col.eq.k )then
                kk=j ! k 行の対角位置
                exit
             endif
          enddo

          ! 下三角成分
          ! LU(i,k)=LU(i,k)/LU(k,k)
          crsMatLU%crsA(kCol) = crsMatLU%crsA(kCol) / crsMatLU%crsA(kk)

          ! 上三角成分
          ! LU(i,j)=LU(i,j)-LU(i,k)*LU(k,j)
          do jCol=crsMatLU%crsRow_ptr(i),crsMatLU%crsRow_ptr(i+1)-1
             j=crsMatLU%crsCol(jCol)-istart+1
             if(j<k+1)then
                cycle
             endif

             ! (i,j)位置が存在する場合
             ! (k,j)位置の探査
             do kj=kk,crsMatLU%crsRow_ptr(k+1)-1
                Col=crsMatLU%crsCol(kj)-istart+1
                if(Col.eq.j)then
                   crsMatLU%crsA(jCol) = crsMatLU%crsA(jCol) - crsMatLU%crsA(kCol) * crsMatLU%crsA(kj) 
                   exit
                endif
             enddo
          enddo

       enddo
    enddo

    call divide_ILU_crs( &
         n,istart,&
         crsMatLU,crsMatL,crsMatU)

    call free_crsMat(crsMatLU)

  end subroutine ILU_LU_crs

  subroutine ILU_LU_dia(&
       n,istart,nnd,&
       diaA,offset,&
       diaMatL,diaMatU)
    use krylov
    use PARCEL_integer
    implicit none

    integer(kind=intType) n,istart,nnd
    real*8 diaA(nnd*n)
    integer(kind=intType) offset(nnd)
    type(DIA_Mat) :: diaMatL
    type(DIA_Mat) :: diaMatU
    type(DIA_Mat) :: diaMatLU

    integer(kind=intType) i,j,k
    integer(kind=intType) ij,jk,kj,kk,ik
    integer(kind=intType) kCol,jCol
    integer(kind=intType) Row,Col
    integer(kind=intType) nnd_ofset0

    call init_diaMat(diaMatLU,n,nnd) 
    do i=1,nnd*n
       diaMatLU%diaA(i)   = diaA(i)
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
       if (diaMatLU%diaA(i+n*(nnd_ofset0-1)).eq.0.0d0) then
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
          diaMatLU%diaA(i+n*(k-1)) = &
             diaMatLU%diaA(i+n*(k-1)) / &
             diaMatLU%diaA(kCol+n*(nnd_ofset0-1))
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
             if(diaMatLU%diaA(i+n*(j-1)).eq.0.0d0)then
                cycle
             endif
             !             write(*,*)"i,j",i,i+diaMatLU%offset(j)
             do kj=nnd_ofset0+1,nnd
                if(jCol.eq.(kCol+diaMatLU%offset(kj)))then
                   !                    write(*,*)"i,k",i,i+diaMatLU%offset(k)
                   !                    write(*,*)"k,j",i+diaMatLU%offset(k),i+diaMatLU%offset(j)            
                   diaMatLU%diaA(i+n*(j-1)) = &
                        diaMatLU%diaA(i+n*(j-1)) - &
                        diaMatLU%diaA(i+n*(k-1)) * diaMatLU%diaA( kCol + n*(kj-1))
                   exit
                endif
             enddo
             !              diaMatLU%diaA(i+n*(j-1)) = &
             !                   diaMatLU%diaA(i+n*(j-1)) - &
             !                   diaMatLU%diaA(i+n*(k-1)) * diaMatLU%diaA( kCol + n*(j-1))
          enddo
       enddo
    enddo

    call divide_ILU_dia( &
         n,istart,nnd,&
         diaMatLU,  &
         diaMatL, &
         diaMatU  &
         )

    call free_diaMat(diaMatLU)


  end subroutine ILU_LU_dia

  subroutine get_ILU_LU_crs(method,crsMat)
    use krylov
    use omp_lib
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

    n      = method%prm%n
    istart = method%prm%istart
    precon_thblock=method%precon_thblock 

    call get_maxthreads(numth)

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

    call init_crsMat(method%preconILU(myth+1)%crsMatLU,n_th,nnz_th)
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
                   method%preconILU(myth+1)%crsMatLU%crsA(nnz_th)   = crsMat%crsA(j) 
                   method%preconILU(myth+1)%crsMatLU%crsCol(nnz_th) = Col
                   nnz_th=nnz_th+1
                else

!                   if( Col .eq. (i + istart - 1) )then
                   if( Col .eq. (i + istart - 1) )then
                      method%preconILU(myth+1)%crsMatLU%crsA(nnz_th)   = crsMat%crsA(j) 
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

    call ILU_LU_crs( & 
         method%preconILU(myth+1)%n_th,&
         method%preconILU(myth+1)%istart_th,&
         method%preconILU(myth+1)%nnz_th,&
         method%preconILU(myth+1)%crsMatLU%crsA,&
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

    call free_crsMat(method%preconILU(myth+1)%crsMatLU)

    !$OMP END PARALLEL
    call omp_set_num_threads(max_numth)

    ! ----

  end subroutine get_ILU_LU_crs

  subroutine get_ILU_LU_dia(method,diaMat)
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
    call init_diaMat(&
         method%preconILU(myth+1)%diaMatLU,n_th,nnd)
#endif

    do i=1,nnd
       method%preconILU(myth+1)%diaMatLU%offset(i) = diaMat%offset(i)
    enddo
    independ_nvec = method%independ_nvec
#if 1

    do j=1,nnd       
       do ii = n_startth,n_endth , independ_nvec
          iimin = min(ii+independ_nvec-1,n_endth)
          do i = ii,min(ii+independ_nvec-1,n_endth)
             Row = i-1+istart
             Col = Row+diaMat%offset(j)
             ! 特定の範囲内をゼロにする
             if( (istart_th<=Col).and.(Col<=iend_th ) )then
                if( ( Col < (ii + istart - 1) ) .or. ( (iimin + istart - 1) < Col ) )then
                   method%preconILU(myth+1)%diaMatLU%diaA(i-n_startth+1+(j-1)*n_th) = diaMat%diaA(i+(j-1)*n)
                else
                   if( Col .eq. (i + istart - 1) )then
                      method%preconILU(myth+1)%diaMatLU%diaA(i-n_startth+1+(j-1)*n_th) = diaMat%diaA(i+(j-1)*n)
                   else
                      method%preconILU(myth+1)%diaMatLU%diaA(i-n_startth+1+(j-1)*n_th) = 0.0d0
                   endif
                endif
             else
                method%preconILU(myth+1)%diaMatLU%diaA(i-n_startth+1+(j-1)*n_th) = 0.0d0
             endif
          enddo
       enddo
    enddo

#else
    do j=1,nnd
       do i=n_startth,n_endth
          Row = i-1+istart
          Col = Row+diaMat%offset(j)
          if( (istart_th<=Col).and.(Col<=iend_th ) )then
             method%preconILU(myth+1)%diaMatLU%diaA(i-n_startth+1+(j-1)*n_th) = diaMat%diaA(i+(j-1)*n)
          else
             method%preconILU(myth+1)%diaMatLU%diaA(i-n_startth+1+(j-1)*n_th) = 0.0d0
          endif
       enddo
    enddo
#endif
    method%preconILU(myth+1)%istart_th = istart_th
    method%preconILU(myth+1)%iend_th   = iend_th
    method%preconILU(myth+1)%n_th      = n_th
    ! -----------
    call ILU_LU_dia(&
         method%preconILU(myth+1)%n_th,&
         method%preconILU(myth+1)%istart_th,&
         nnd,&
         method%preconILU(myth+1)%diaMatLU%diaA , &
         method%preconILU(myth+1)%diaMatLU%offset , &
         method%preconILU(myth+1)%diaMatL , &
         method%preconILU(myth+1)%diaMatU )

    method%preconILU(myth+1)%nnd_th  = nnd
    method%preconILU(myth+1)%nndL_th =  &
         size(method%preconILU(myth+1)%diaMatL%diaA) / &
         method%preconILU(myth+1)%n_th
    method%preconILU(myth+1)%nndU_th = &
         size(method%preconILU(myth+1)%diaMatU%diaA) / &
         method%preconILU(myth+1)%n_th
    call free_diaMat(method%preconILU(myth+1)%diaMatLU)

    !$OMP END PARALLEL
    ! ----
    call omp_set_num_threads(max_numth)

  end subroutine get_ILU_LU_dia

  ! 名前があまり適切ではない?
  subroutine solve_LU_crs( &
       n,istart,&
       nnzL,crsAL,crsRow_ptrL,crsColL, &
       nnzU,crsAU,crsRow_ptrU,crsColU, &
       r,s, &
       independ_nvec &
       )
    use PARCEL_integer
    implicit none
    integer(kind=intType) n,istart

    integer(kind=intType) nnzL
    real*8 crsAL(nnzL)
    integer(kind=intType) crsRow_ptrL(n+1)
    integer(kind=intType) crsColL(nnzL)

    integer(kind=intType) nnzU
    real*8 crsAU(nnzU)
    integer(kind=intType) crsRow_ptrU(n+1)
    integer(kind=intType) crsColU(nnzU)

    real*8 r(n),s(n)
    integer independ_nvec

    real*8 z
    real*8 zz(independ_nvec)
    integer(kind=intType) i,j
    integer(kind=intType) ii0,ii,jj
    integer(kind=intType) iistep,iimin
    integer(kind=intType) Col,Row


    ! CRS形式の行列は正方行列を仮定
    ! 零以外の対角成分を持つことを仮定
    ! -----
    ! LUs = r をsについて解いていく
    
    if(independ_nvec.eq.1)then
       ! -----
       ! z = Us 
       ! Lz = r をzについて解く
       ! 前進代入
       do i=1,n
          z=r(i)
          ! 対角成分を持つことを仮定
          ! 下三角行列では各行の最後の値が対角成分
          do j=crsRow_ptrL(i),crsRow_ptrL(i+1)-1 -1
             Col=crsColL(j)-istart+1
             z=z-s(Col)*crsAL(j)
          enddo
          j=crsRow_ptrL(i+1)-1
          s(i)=z/crsAL(j)
          !       write(*,*)"L^-1 r = ",i,s(i)
       enddo
       !    write(*,*)"--------"
       ! -----
       ! Us=z をsについて解く
       ! 後退代入
       do i=n,1,-1
          z=s(i)
          ! 対角成分を持つことを仮定
          ! 上三角行列では最初の値が対角成分
          do j=crsRow_ptrU(i)+1,crsRow_ptrU(i+1)-1
             Col=crsColU(j)-istart+1
             z=z-s(Col)*crsAU(j)
          enddo
          j=crsRow_ptrU(i)
          s(i)=z/crsAU(j)
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
             zz(i-ii+1)=r(i)
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
                   zz(i-ii+1)=zz(i-ii+1)-s(Col)*crsAL(j)
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
             s(i)=zz(i-ii+1)/crsAL(j)
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
             zz(ii-i+1)=s(i)
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
                   zz(ii-i+1)=zz(ii-i+1)-s(Col)*crsAU(j)
                enddo
             enddo

             ii0 = ii0-iistep-1
             if(ii0<iimin)then
                exit
             endif

          enddo

          do i = ii ,iimin ,-1
             j=crsRow_ptrU(i)
             s(i)=zz(ii-i+1)/crsAU(j)
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

  end subroutine solve_LU_crs

  ! 名前があまり適切ではない?
  subroutine solve_LU_dia( &
       n,istart,&
       nndL,diaAL,offsetL, &
       nndU,diaAU,offsetU, &
       r,s, &
       independ_nvec &
       )
    use PARCEL_integer
    implicit none
    integer(kind=intType) n,istart

    integer(kind=intType) nndL
    real*8 diaAL(n*nndL)
    integer(kind=intType) offsetL(nndL)

    integer(kind=intType) nndU
    real*8 diaAU(n*nndU)
    integer(kind=intType) offsetU(nndU)

    real*8 r(n),s(n)
    integer independ_nvec
    real*8 z
    real*8 zz(independ_nvec)
    integer(kind=intType) i,j
    integer(kind=intType) ii,iimin
    integer(kind=intType) Col,Row

    ! 零以外の対角成分を持つことを仮定
    ! -----
    ! LUs = r をsについて解いていく
    if(independ_nvec.eq.1)then
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
!          write(*,*)j,"save -- ",i,i+(j-1)*n,diaAL_(i+(j-1)*n)
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
             zz(i-ii+1)=r(i)
          enddo          
          ! 対角成分を持つことを仮定
          ! 下三角行列で対角成分はnndLのオフセット
          do j=1,nndL-1
             do i = ii ,iimin
                if( (1 <= (i+offsetL(j))).and.((i+offsetL(j)) <= n ) )then
                   zz(i-ii+1)=zz(i-ii+1)-s(i+offsetL(j))*diaAL(i+(j-1)*n)
                endif
             enddo
          enddo
          do i = ii ,iimin
             j=nndL
             s(i)=zz(i-ii+1)/diaAL(i+(j-1)*n)
          enddo
       enddo

       ! -----
       ! Us=z をsについて解く
       ! 後退代入
       do ii = n,1,-independ_nvec
          iimin = max(ii-independ_nvec+1,1) 
          do i = ii ,iimin ,-1
             zz(ii-i+1)=s(i)
          enddo

          ! 対角成分を持つことを仮定
          ! 上三角行列の対角成分は1のオフセット
          do j=1+1,nndU
             do i = ii ,iimin ,-1
                if( ( ( 1 <= (i+offsetU(j)) ) ).and.( (i+offsetU(j) <= n ) ) )then
                   zz(ii-i+1)=zz(ii-i+1)-s(i+offsetU(j))*diaAU(i+(j-1)*n)
                endif
             enddo
          enddo
          do i = ii ,iimin ,-1          
             j=1
             s(i)=zz(ii-i+1)/diaAU(i+(j-1)*n)
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

  end subroutine solve_LU_dia

  subroutine blockJacobi_ILU_crs(method,r,s)
    use omp_lib
    use krylov
    use PARCEL_integer

    implicit none
    type(KSP) :: method
    real*8 r(method%prm%n)
    real*8 s(method%prm%n)

    integer(kind=intType) n
    integer(kind=intType) nnz
    integer(kind=intType) istart,iend

    integer myth,numth,max_numth
    integer(kind=intType) n_th
    integer(kind=intType) n_tmp
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th
    integer precon_thblock

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
    ! ---
    !     call view_crs(&
    !          method%preconILU(myth+1)%n_th,&
    !          method%preconILU(myth+1)%istart_th,&
    !          method%preconILU(myth+1)%crsMatU )
    ! ---
    call solve_LU_crs( &
         method%preconILU(myth+1)%n_th, &
         method%preconILU(myth+1)%istart_th, &
         method%preconILU(myth+1)%nnzL_th, &
         method%preconILU(myth+1)%crsMatL%crsA, &
         method%preconILU(myth+1)%crsMatL%crsRow_ptr, &
         method%preconILU(myth+1)%crsMatL%crsCol, &
         method%preconILU(myth+1)%nnzU_th, &
         method%preconILU(myth+1)%crsMatU%crsA, &
         method%preconILU(myth+1)%crsMatU%crsRow_ptr, &
         method%preconILU(myth+1)%crsMatU%crsCol, &
         r(n_startth), &
         s(n_startth), &
         method%independ_nvec &
         )

    !$OMP END PARALLEL

!    if(n<=max_numth)then
    call omp_set_num_threads(max_numth)
!    endif

  end subroutine blockJacobi_ILU_crs

  subroutine blockJacobi_ILU_dia(method,r,s)
    use omp_lib
    use krylov
    use PARCEL_integer

    implicit none
    type(KSP) :: method
    real*8 r(method%prm%n)
    real*8 s(method%prm%n)

    integer(kind=intType) n
    integer(kind=intType) nnz
    integer(kind=intType) istart,iend

    integer myth,numth,max_numth
    integer(kind=intType) n_th
    integer(kind=intType) n_tmp
    integer(kind=intType) n_startth,n_endth
    integer(kind=intType) istart_th,iend_th
    integer precon_thblock

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
    call solve_LU_dia( &
         method%preconILU(myth+1)%n_th, &
         method%preconILU(myth+1)%istart_th, &
         method%preconILU(myth+1)%nndL_th, &
         method%preconILU(myth+1)%diaMatL%diaA, &
         method%preconILU(myth+1)%diaMatL%offset, &
         method%preconILU(myth+1)%nndU_th, &
         method%preconILU(myth+1)%diaMatU%diaA, &
         method%preconILU(myth+1)%diaMatU%offset, &
         r(n_startth), &
         s(n_startth), &
         method%independ_nvec &
         )

    !$OMP END PARALLEL

!    if(n<=max_numth)then
    call omp_set_num_threads(max_numth)
!    endif

  end subroutine blockJacobi_ILU_dia

end module matrix_kernel

