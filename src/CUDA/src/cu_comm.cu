#include <stdio.h>
#include <mpi.h>
#include <cuda.h>
#include <cuda_runtime.h>

#include "cu_common.h"
#include "reduciton_cuda.hpp"
#include "cu_comm.hpp"

extern "C" { 
  void get_nrk_vecLU_add(
			 devicemem *rbuf_,PARCEL_INT offset_rbuf,
			 PARCEL_INT n,
			 PARCEL_INT nColL,
			 PARCEL_INT nColU,		    
			 devicemem *vec_ ,PARCEL_INT offset_vec
			 ){
    double *vec  = (double *)vec_[0].a_d;
    double *rbuf = (double *)rbuf_[0].a_d;

    cudaDeviceSynchronize();
    const size_t len = sizeof(double) * (nColL + nColU) ;
    cudaDeviceSynchronize();
    cudaMemcpy(
	       &vec[n+offset_vec],
	       &rbuf[offset_rbuf],
	       len,
	       cudaMemcpyDeviceToDevice);
    cudaDeviceSynchronize();
    return ;
  }
}


extern "C" {

  void cuda_nrk_vec_pack(
		    int numGPUblock,
		    int numGPUthread,
		    devicemem *vec_,  PARCEL_INT offset_vec,
		    devicemem *sbuf_, PARCEL_INT offset_sbuf,
		    devicemem *offset,
		    PARCEL_INT ndata,
		    PARCEL_INT i0
		    ){
    double *vec  = (double*)vec_[0].a_d;
    double *sbuf = (double*)sbuf_[0].a_d;
    nrk_vec_pack_<double> <<<numGPUblock,numGPUthread>>>(
						 i0, 
						 ndata, 
						 &vec[offset_vec],
						 &sbuf[offset_sbuf],
						 (PARCEL_INT *) offset[0].a_d);
    return ;
  }

  void cuda_get_nrk_vecLU(
		     int numGPUblock,
		     int numGPUthread,
		     devicemem *rbuff_, PARCEL_INT offset_rbuff,
		     PARCEL_INT nColL,
		     PARCEL_INT nColU,
		     devicemem *rbuffL_, PARCEL_INT offset_rbuffL,
		     devicemem *rbuffU_, PARCEL_INT offset_rbuffU
		     ){
    double* rbuff  = (double*)rbuff_[0].a_d;
    double* rbuffL = (double*)rbuffL_[0].a_d;
    double* rbuffU = (double*)rbuffU_[0].a_d;

    get_nrk_vecLU_<double> <<<numGPUblock,numGPUthread>>>(
							  &rbuff[offset_rbuff],
							  nColL,
							  nColU,
							  &rbuffL[offset_rbuffL],
							  &rbuffU[offset_rbuffU]
							  );
    
    return ;
  }

#if 1

  void cuda_nrk_vec_comm_FP64(
			      MPI_Fint f_comm,
			      devicemem *sbuf_,
			      devicemem *rbuf_,
			      int nrk_rank_send,
			      int *nrk_prm_send_offset,
			      int *nrk_prm_send_ndata,
			      int *nrk_prm_send_rank,
			      int nrk_rank_recv,
			      int *nrk_prm_recv_offset,
			      int *nrk_prm_recv_ndata,
			      int *nrk_prm_recv_rank,
			      int myrank
			      ){

    cudaDeviceSynchronize();
    MPI_Comm comm = MPI_Comm_f2c(f_comm);

#if 1
    double *sbuf = (double *)sbuf_[0].a_d;
    double *rbuf = (double *)rbuf_[0].a_d;
#else
    int nsend = 0;
    for(int i=0;i<nrk_rank_send;i++){	  
      nsend = nsend + nrk_prm_send_ndata[i];
    }
    size_t len_nsend = sizeof(double)*nsend;
    double *sbuf;cudaMalloc(&(sbuf), len_nsend);

    int nrecv = 0;
    for(int i=0;i<nrk_rank_recv;i++){	  
      nrecv = nrecv + nrk_prm_recv_ndata[i];
    }
    size_t len_nrecv = sizeof(double)*nrecv;
    double *rbuf;cudaMalloc(&(sbuf), len_nrecv);


    cudaMemcpy(
	       sbuf,
	       (double *)sbuf_[0].a_d,
	       len_nsend,
	       cudaMemcpyDeviceToDevice);

#endif
    MPI_Request s_req[nrk_rank_send];
    MPI_Request r_req[nrk_rank_recv];

    MPI_Status s_sta[nrk_rank_send];
    MPI_Status r_sta[nrk_rank_recv];
    int offset = 1;

    for(int i=nrk_rank_recv-1;i>=0;i--){
       // printf("myrank = %d : send rank = %d, ndata = %d, offset = %d \n",
       // 	     myrank,
       // 	     nrk_prm_recv_rank[i],
       // 	     nrk_prm_recv_ndata[i],
       // 	     offset-1 + nrk_prm_send_offset[i]
       // 	     );
      int ierr = MPI_Irecv(
			   &rbuf[ nrk_prm_recv_offset[i] ],
			   nrk_prm_recv_ndata[i],
			   MPI_DOUBLE,
			   nrk_prm_recv_rank[i],
			   nrk_prm_recv_rank[i],
			   comm,
			   &r_req[i]);
    }

    for(int i=0;i<nrk_rank_send;i++){	  
       // printf("myrank = %d : send rank = %d, ndata = %d, offset = %d \n",
       // 	     myrank,
       // 	     nrk_prm_send_rank[i],
       // 	     nrk_prm_send_ndata[i],
       // 	     offset-1 + nrk_prm_send_offset[i]
       // 	     );
      int ierr = MPI_Isend(
			   &sbuf[ nrk_prm_send_offset[i] ],
			   nrk_prm_send_ndata[i],
			   MPI_DOUBLE,
			   nrk_prm_send_rank[i],
			   myrank,
			   comm,
			   &s_req[i]
			   );
    }

    for(int i=0;i<nrk_rank_send;i++){
      int ierr = MPI_Wait(&s_req[i],&s_sta[i]);
    }
    for(int i=nrk_rank_recv-1;i>=0;i--){
      int ierr = MPI_Wait(&r_req[i],&r_sta[i]);
    }

#if 0

    cudaMemcpy(
	       (double *)rbuf_[0].a_d,
	       rbuf,
		 len_nrecv,
	       cudaMemcpyDeviceToDevice);

    cudaFree(sbuf);
    cudaFree(rbuf);
#endif

    return ;
  }


#else
  void cuda_nrk_vec_comm_FP64(
			      MPI_Fint f_comm,
			      devicemem *sbuf_,
			      devicemem *rbuf_,
			      int nrk_rank_send,
			      int *nrk_prm_send_offset,
			      int *nrk_prm_send_ndata,
			      int *nrk_prm_send_rank,
			      int nrk_rank_recv,
			      int *nrk_prm_recv_offset,
			      int *nrk_prm_recv_ndata,
			      int *nrk_prm_recv_rank,
			      int myrank
			      ){

    cudaDeviceSynchronize();
    MPI_Comm comm = MPI_Comm_f2c(f_comm);
    //    MPI_Comm comm = MPI_COMM_WORLD;
#if 0
    double *sbuf = (double *)sbuf_[0].a_d;
    double *rbuf = (double *)rbuf_[0].a_d;
#else
    int nsend = 0;
    for(int i=0;i<nrk_rank_send;i++){	  
      nsend = nsend + nrk_prm_send_ndata[i];
    }
    size_t len_nsend = sizeof(double)*nsend;
    double *sbuf;cudaMalloc(&(sbuf), len_nsend);

    int nrecv = 0;
    for(int i=0;i<nrk_rank_recv;i++){	  
      nrecv = nrecv + nrk_prm_recv_ndata[i];
    }
    size_t len_nrecv = sizeof(double)*nrecv;
    double *rbuf;cudaMalloc(&(sbuf), len_nrecv);
#endif
    MPI_Request s_req[nrk_rank_send];
    MPI_Request r_req[nrk_rank_recv];

    MPI_Status s_sta[nrk_rank_send];
    MPI_Status r_sta[nrk_rank_recv];
    int offset = 1;

    for(int i=0;i<nrk_rank_send;i++){	  
       // printf("myrank = %d : send rank = %d, ndata = %d, offset = %d \n",
       // 	     myrank,
       // 	     nrk_prm_send_rank[i],
       // 	     nrk_prm_send_ndata[i],
       // 	     offset-1 + nrk_prm_send_offset[i]
       // 	     );
      int ierr = MPI_Isend(
			   &sbuf[ nrk_prm_send_offset[i] ],
			   nrk_prm_send_ndata[i],
			   MPI_DOUBLE,
			   nrk_prm_send_rank[i],
			   myrank,
			   comm,
			   &s_req[i]
			   );
    }

    for(int i=nrk_rank_recv-1;i>=0;i--){
       // printf("myrank = %d : send rank = %d, ndata = %d, offset = %d \n",
       // 	     myrank,
       // 	     nrk_prm_recv_rank[i],
       // 	     nrk_prm_recv_ndata[i],
       // 	     offset-1 + nrk_prm_send_offset[i]
       // 	     );
      int ierr = MPI_Irecv(
			   &rbuf[ nrk_prm_recv_offset[i] ],
			   nrk_prm_recv_ndata[i],
			   MPI_DOUBLE,
			   nrk_prm_recv_rank[i],
			   nrk_prm_recv_rank[i],
			   comm,
			   &r_req[i]);
    }

    for(int i=0;i<nrk_rank_send;i++){
      int ierr = MPI_Wait(&s_req[i],&s_sta[i]);
    }
    for(int i=nrk_rank_recv-1;i>=0;i--){
      int ierr = MPI_Wait(&r_req[i],&r_sta[i]);
    }

#if 1
    cudaFree(sbuf);
    cudaFree(rbuf);
#endif

    return ;
  }
#endif

}
