#include <stdio.h>
#include <iostream>
#include <string>
#include <map>
#include <vector>

#ifdef PARCEL_INT8
#define PARCEL_INT int64_t
#else
#define PARCEL_INT int
#endif

extern "C"{
  void init_Matvec_crs_iLs(
			   PARCEL_INT istart,
			   PARCEL_INT iend,
			   PARCEL_INT n,
			   PARCEL_INT nColL,
			   PARCEL_INT nColU,
			   PARCEL_INT *crsCol,
			   PARCEL_INT *crsRow_ptr,
			   PARCEL_INT *iLs,	
			   PARCEL_INT *RowL,
			   PARCEL_INT *RowU		   
			   ){
    std::map<PARCEL_INT,PARCEL_INT> crsRowL_id;
    std::map<PARCEL_INT,PARCEL_INT> crsRowU_id;
    for(PARCEL_INT i=0;i<nColL;i++){
      crsRowL_id[ RowL[i] ] = i+1;
    }
    for(PARCEL_INT i=0;i<nColU;i++){
      crsRowU_id[ RowU[i] ] = i+1;
    }

    // #pragma omp parallel for
    for(PARCEL_INT k=0;k<n;k++){
      for(PARCEL_INT j=crsRow_ptr[k];j<crsRow_ptr[k+1];j++){
	PARCEL_INT Col = crsCol[j-1];
	iLs[j-1] = Col - istart+1 ;
	if( Col<istart ){
	  iLs[j-1] = crsRowL_id[ Col ];
	  continue ;
	}
	if( iend<Col ){
	  iLs[j-1] = crsRowU_id[ Col ];
	  continue ;
	}
      }
    }
    return ;
  }

  // 列番号からcusparse形式行列に対応した列番号を与える
  void get_cusparse_sort(
			PARCEL_INT istart,
			PARCEL_INT iend,
			PARCEL_INT n,
			PARCEL_INT nnz,
			PARCEL_INT nnz_comm,
			PARCEL_INT maxcrsRow_ptr,
			double *crsA,
			PARCEL_INT *crsCol,
			PARCEL_INT *crsCol_comm,
			PARCEL_INT *crsRow_ptr
			){
    //    printf("get_cusparse_sort \n");
    // 連想配列の作成
    std::map<PARCEL_INT,PARCEL_INT> crsCollocalid;
    for(PARCEL_INT i=0;i<nnz_comm;i++){
      crsCollocalid[crsCol_comm[i]] = i + 1; 
    }

    //    std::vector< double > crsA_local_(maxcrsRow_ptr,0.0);
    //    std::vector< PARCEL_INT > crsCol_local_(maxcrsRow_ptr,0.0);
    // double crsA_local[maxcrsRow_ptr];
    // PARCEL_INT crsCol_local[maxcrsRow_ptr];


#pragma omp parallel for
    for(PARCEL_INT i=0;i<n;i++){      
    double crsA_local[maxcrsRow_ptr];
    PARCEL_INT crsCol_local[maxcrsRow_ptr];

      PARCEL_INT jdiags = crsRow_ptr[i];
      PARCEL_INT jdiage = crsRow_ptr[i+1]-1;
      for(PARCEL_INT j=crsRow_ptr[i];j<crsRow_ptr[i+1];j++){
	PARCEL_INT Col = crsCol[j-1];
	jdiags = j;
	if( (istart<=Col) && (Col<=iend) ){
	  break ;
	}
      }
      for(PARCEL_INT j=crsRow_ptr[i+1]-1;crsRow_ptr[i]<=j;j--){
	PARCEL_INT Col = crsCol[j-1];
	jdiage = j;
	if( (istart<=Col) && (Col<=iend) ){
	  break ;
	}
      }

      // cusparse用にゼロ始まりにする
      // PARCEL_INT nCol = crsRow_ptr[i];

      // 対角
      PARCEL_INT nCol = 1;
      for(PARCEL_INT j=jdiags;j<=jdiage;j++){
	crsCol_local[nCol-1] = crsCol[j-1] - istart + 1 -1;
	crsA_local[nCol-1]   = crsA[j-1];
	nCol = nCol + 1;
      }
      // Lランク
      for(PARCEL_INT j = crsRow_ptr[i];j<jdiags;j++){
	PARCEL_INT Col = crsCol[j-1];
	crsCol_local[nCol-1] = crsCollocalid[Col] + n -1;
	crsA_local[nCol-1]   = crsA[j-1];
	nCol = nCol + 1;
      }
      // Uランク
      for(PARCEL_INT j=jdiage+1;j<crsRow_ptr[i+1];j++){
	PARCEL_INT Col = crsCol[j-1];
	crsCol_local[nCol-1] = crsCollocalid[Col] + n -1;
	crsA_local[nCol-1]   = crsA[j-1];
	nCol = nCol + 1;
      }

      for(PARCEL_INT j=crsRow_ptr[i];j<crsRow_ptr[i+1];j++){
	crsA[j-1]   = crsA_local[j-crsRow_ptr[i]]   ;
 	crsCol[j-1] = crsCol_local[j-crsRow_ptr[i]] ;
      }

    }

    // cusparse用にゼロ始まりにする
    for(PARCEL_INT i=0;i<n+1;i++){
      crsRow_ptr[i] = crsRow_ptr[i] -1;
    }

    //    return ;

  }

}
