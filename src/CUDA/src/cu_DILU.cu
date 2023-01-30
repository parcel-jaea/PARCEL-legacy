#include <cuda.h>
#include <cuda_runtime.h>

#include "cu_common.h"

template <class Tdata,class Tcalc>
__global__
void solve_DILU_crs_(
		     Tdata *crsA,
		     PARCEL_INT *crsRow_ptr,
		     PARCEL_INT *crsCol,
		     PARCEL_INT istart,
		     PARCEL_INT iend,
		     PARCEL_INT nn_th_max, // 各スレッドの(n_startth:n_endth)の中で最大のものを選択
		     PARCEL_INT *istart_th_,
		     PARCEL_INT *iend_th_,
		     PARCEL_INT *n_startth_,
		     PARCEL_INT *n_endth_,
		     Tdata *z_,
		     Tdata *Dinv_,
		     Tdata *r_,
		     Tdata *s_,
		     PARCEL_INT *ILU_Lstart_,
		     PARCEL_INT *ILU_Lend_,
		     PARCEL_INT *ILU_Ustart_,
		     PARCEL_INT *ILU_Uend_,
		     PARCEL_INT independ_nvec,
		     PARCEL_INT BlockJacobi_nBlock
		     ){
  const int tid  = blockIdx.x * blockDim.x + threadIdx.x;
  const int step = blockDim.x*gridDim.x;

  for(PARCEL_INT bj_id = tid; bj_id<BlockJacobi_nBlock;bj_id++){
  Tdata *z = z_[tid*independ_nvec];

  PARCEL_INT istart_th = istart_th_[bj_id-1];
  PARCEL_INT iend_th   = iend_th_[bj_id-1];
  PARCEL_INT n_startth = n_startth_[bj_id-1];
  PARCEL_INT n_endth   = n_endth_[bj_id-1];

  Tdata *r    = r_[n_startth_[bj_id]-1];
  Tdata *s    = s_[n_startth_[bj_id]-1];
  Tdata *Dinv = Dinv_[n_startth_[bj_id]-1];

  PARCEL_INT *ILU_Lstart = &ILU_Lstart_[bj_id*nn_th_max];
  PARCEL_INT *ILU_Lend   = &ILU_Lend_[bj_id*nn_th_max];
  PARCEL_INT *ILU_Ustart = &ILU_Ustart_[bj_id*nn_th_max];
  PARCEL_INT *ILU_Uend   = &ILU_Uend_[bj_id*nn_th_max];

  if(independ_nvec==1){
    //-----
    //  z = D^-1(D+U)s
    //  (D+L)z=r
    for(PARCEL_INT k=n_startth;k<=n_endth;k++){
      Tdata zz=r[k-1];      
      for(PARCEL_INT j=ILU_Lstart[k-1];j<=ILU_Lend[k-1];j++){
	PARCEL_INT Col = crsCol[j-1]-istart + 1;
	//r - Lz の計算
	zz=zz-s[Col-1]*crsA[j-1];
      }
      //z = D^-1 (r - Lz ) 
      s[k-1]=Dinv[k-1]*zz;
    }
       
    //-----
    //D^-1(D+U)s=z
    //  (I+D^-1U)s=z 
    //  s=z - D^-1 Us       
    for(PARCEL_INT k=n_endth;k>=n_startth;k--){
      Tdata zz=0.0;
      for(PARCEL_INT j=ILU_Ustart[k-1];j>=ILU_Uend[k-1];j--){
	PARCEL_INT Col=crsCol[j-1]-istart+1;
	//Us の計算
	zz=zz+s[Col-1]*crsA[j-1];
      }
      s[k-1] = s[k-1] - Dinv[k-1]*zz;
    }       
  }else{
    //-----
    //  z = D^-1(D+U)s
    //  (D+L)z=r
    for(PARCEL_INT kk=n_startth;kk<=n_endth;kk=kk+independ_nvec){
      PARCEL_INT kmax = (kk+independ_nvec-1) < n_endth ? (kk+independ_nvec-1) : n_endth;
      for(PARCEL_INT k=kk;k<=kmax;k++){
	z[k-kk+1-1]=r[k-1];
      }
      PARCEL_INT kk0 = kk;
      PARCEL_INT kkmin = (kk+independ_nvec-1) < n_endth ? (kk+independ_nvec-1) : n_endth;
      while(1){
	PARCEL_INT kkstep = 0;
	for(PARCEL_INT k = kk0;k<=kkmin-1 ;k++){
	  if ( ( ILU_Lend[k-1]-ILU_Lstart[k-1] + 1 ) == ( ILU_Lend[k+1-1]-ILU_Lstart[k+1-1] + 1 ) ) {
	    kkstep=kkstep+1;
	  }else{
	    break;
	  }
	}          
	for(PARCEL_INT jj=1;jj<=(ILU_Lend[kk0-1]-ILU_Lstart[kk0-1]+1);jj++){
	  for(PARCEL_INT k = kk0;k<=kk0+kkstep;k++){
	    PARCEL_INT j   = ILU_Lstart[k-1]+jj-1;
	    PARCEL_INT Col = crsCol[j-1]-istart+1;
	    //r - Lz の計算
	    z[k-kk+1-1]=z[k-kk+1-1]-s[Col-1]*crsA[j-1];
	  }
	}
	kk0 = kk0+kkstep+1;
	if(kk0>kkmin){
	  break;
	}
      }
      kmax = (kk+independ_nvec-1) < n_endth ? (kk+independ_nvec-1) : n_endth;
      for(PARCEL_INT k=kk;k<=kmax;k++){
	s[k-1]=Dinv[k-1]*z[k-kk+1-1];
      }
    }

    //-----
    //D^-1(D+U)s=z
    //  (I+D^-1U)s=z 
    //  s=z - D^-1 Us
    for(PARCEL_INT kk=n_endth;kk<=n_startth;kk=kk-independ_nvec){
      PARCEL_INT kmax = (kk-independ_nvec+1) < n_startth ? n_startth : (kk-independ_nvec+1);
      for(PARCEL_INT k = kk;k>=kmax;k--){
	z[kk-k+1-1]=0.0;
      }

      PARCEL_INT kk0 = kk;
      PARCEL_INT kkmin = (kk-independ_nvec+1) < n_startth ? n_startth : (kk-independ_nvec+1) ;
      while(1){
	PARCEL_INT kkstep = 0;
	for(PARCEL_INT k = kk0;k>=kkmin+1;k--){
	  if ( ( ILU_Ustart[k-1]-ILU_Uend[k-1] + 1 ) == ( ILU_Ustart[k-1-1] - ILU_Uend[k-1-1] + 1 ) ){
	    kkstep=kkstep+1;
	  }else{
	    break;
	  }
	}
	for(PARCEL_INT jj=1;jj<=(ILU_Ustart[kk0-1]-ILU_Uend[kk0-1]+1);jj++){
	  for(PARCEL_INT k = kk0;k>=kk0-kkstep;k--){
	    PARCEL_INT j = ILU_Ustart[k-1]-jj+1;
	    PARCEL_INT Col = crsCol[j-1]-istart+1;
	    z[kk-k+1-1]=z[kk-k+1-1]+s[Col-1]*crsA[j-1];
	  }
	}
	kk0 = kk0-kkstep-1;
	if(kk0<kkmin){
	  break;
	}
      }
      kmax = (kk-independ_nvec+1<n_startth) ? n_startth : (kk-independ_nvec+1) ; 
      for(PARCEL_INT k = kk;k>=kmax;k--){
	s[k-1]=s[k-1]-Dinv[k-1]*z[kk-k+1-1];
      }
    }
  }
  }
 
  return ;

}

#if 0
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

#endif


