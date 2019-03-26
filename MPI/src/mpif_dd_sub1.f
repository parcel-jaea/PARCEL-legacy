!
! mpif_dd_sub1.f
!	4倍精度MPI通信ライブラリ・通信方式１
!
! Attention:
!	mpif_dd.fファイルのサブルーチンから利用します
!
! Version : 1.0 (2011.1.31)
! Author: IMS
!
! Description:
!	4倍精度データの上位データと下位データを
!	別々の配列として2回に分けて送受信する通信方式
!
! =========================================================================

! ====================================================================
! MPI_Send_dd : MPI_SENDの4倍精度MPI版（通信方式１）
!
!	hi       : in    : double-doubleの上位16桁の1次元配列
!	lo       : in    : double-doubleの下位16桁の1次元配列
!	count    : in    : hi, loの配列要素数
!	datatype : in    : データ型
!	dest     : in    : 送信先プロセスのrank
!	tag      : in    : タグ（メッセージ識別番号）
!	comm     : in    : MPIコミュニケータ
!	ierror   : out   : エラーコード
!
      subroutine MPI_Send_dd(
     $              hi, lo, count, datatype, dest, tag, comm, ierror)
      implicit none
      include 'mpif.h'
 
      integer count, datatype, dest, tag, comm, ierror
      double precision hi(*), lo(*)
      
!      call MPI_Send(hi, count, datatype, dest, tag, comm, ierror)
!      
CCCCCCCCCCCCC--MHIR revice--CCCCCCCCCCCCCCCCCCCC
      call MPI_Send(hi, count, datatype, dest, tag, comm, ierror)
      call MPI_Send(lo, count, datatype, dest, tag+1, comm, ierror)

      return
      end

! ====================================================================
! MPI_Recv : MPI_RECVの4倍精度MPI版（通信方式１）
!
!	hi       : out   : double-doubleの上位16桁の1次元配列
!	lo       : out   : double-doubleの下位16桁の1次元配列
!	count    : in    : hi, loの配列要素数
!	datatype : in    : データ型
!	source   : in    : 送信元プロセスのrank
!	tag      : in    : タグ（メッセージ識別番号）
!	comm     : in    : MPIコミュニケータ
!	status   : out   : 送信プロセスの情報
!	ierror   : out   : エラーコード
!
      subroutine MPI_Recv_dd(
     $                      hi, lo, count, datatype, source, tag,
     $                      comm, status, ierror)
      implicit none
      include 'mpif.h'
      
      double precision hi(*), lo(*)
      integer count, datatype, source, tag, comm, ierror
      integer status(MPI_STATUS_SIZE)
      

!      call MPI_Recv(hi, count, datatype, source, tag,
!     $              comm, status, ierror)
!      call MPI_Recv(lo, count, datatype, source, tag,
!     $              comm, status, ierror)
CCCCCCCCC-MHIR revice-CCCCCCCCCCCCCCCC
      call MPI_Recv(hi, count, datatype, source, tag,      
     $              comm, status, ierror)
      call MPI_Recv(lo, count, datatype, source, tag+1,
     $              comm, status, ierror)
      
      return
      end

! ====================================================================
! MPI_Isend_dd : MPI_ISENDの4倍精度MPI版（通信方式１）
!	注：2回送信するため、通信リクエストが長さ2の配列で戻る
!           通信方法2とあわせるため，sdataは本来不要だが，引数に利用している
!
!	hi       : in    : double-doubleの上位16桁の1次元配列
!	lo       : in    : double-doubleの下位16桁の1次元配列
!	sdata    : in    : double-double型配列（不要）
!	count    : in    : hi, loの配列要素数
!	datatype : in    : データ型
!	dest     : in    : 送信先プロセスのrank
!	tag      : in    : タグ（メッセージ識別番号）
!	comm     : in    : MPIコミュニケータ
!	request  : out   : 通信リクエスト
!	ierror   : out   : エラーコード
!
      subroutine MPI_Isend_dd(
     $                      hi, lo, sdata, count, datatype, dest, tag,
     $                      comm, request, ierror)
      implicit none
      include 'mpif.h'
 
      double precision hi(*), lo(*), sdata(*)
      integer request(2)
      integer count, datatype, dest, tag, comm, ierror
      integer ireq_hi, ireq_lo
      
      call MPI_Isend(hi, count, datatype, dest, tag,
     $               comm, ireq_hi, ierror)
      call MPI_Isend(lo, count, datatype, dest, tag+1,
     $               comm, ireq_lo, ierror)

      request(1) = ireq_hi
      request(2) = ireq_lo

      return
      end

! ====================================================================
! MPI_Irecv : MPI_IRECVの4倍精度MPI版（通信方式１）
!	注：2回受信するため、通信リクエストが長さ2の配列で戻る
!
!	hi       : out   : double-doubleの上位16桁の1次元配列
!	lo       : out   : double-doubleの下位16桁の1次元配列
!	count    : in    : hi, loの配列要素数
!	datatype : in    : データ型
!	source   : in    : 送信元プロセスのrank
!	tag      : in    : タグ（メッセージ識別番号）
!	comm     : in    : MPIコミュニケータ
!	request  : out   : 通信リクエスト
!	ierror   : out   : エラーコード
!
!      subroutine MPI_Irecv_dd(
!     $                  hi, lo, recvbuf, count, datatype, source, tag,
!     $                  comm, request, ierror)
      subroutine MPI_Irecv_dd(
     $                  recvbuf, count, datatype, source, tag,
     $                  comm, request, ierror)
      implicit none
      include 'mpif.h'
      
!      double precision hi(*), lo(*)
      double precision recvbuf(*)
      integer count, datatype, source, tag, comm, ierror
      integer request(2)
      integer ireq_hi, ireq_lo
      
      call MPI_Irecv(recvbuf(1), count, datatype, source, tag,
     $               comm, ireq_hi, ierror)
      call MPI_Irecv(recvbuf(count+1), count, datatype, source, tag+1,
     $               comm, ireq_lo, ierror)
!      call MPI_Irecv(hi, count, datatype, source, tag,
!     $               comm, ireq_hi, ierror)
!      call MPI_Irecv(lo, count, datatype, source, tag+1,
!     $               comm, ireq_lo, ierror)
      
      request(1) = ireq_hi
      request(2) = ireq_lo

      return
      end

! ====================================================================
! MPI_Wait_dd : MPI_WAITの4倍精度MPI版（通信方式１）
!	注：送受信を2回行うため、requestを配列で渡す
!           通信方法2とあわせるため，sdataは本来不要だが，引数に利用している
!
!	sdata    : in    : double-double型配列（不要）
!	request  : in    : 通信リクエスト
!	status   : out   : 送信プロセスの情報
!	ierror   : out   : エラーコード
!
      subroutine MPI_Wait_dd(sdata, request, status, ierror)
      implicit none
      include 'mpif.h'

      double precision sdata(*)
      integer request(2)
      integer status(MPI_STATUS_SIZE)
      integer ierror
      

      call MPI_Wait(request(1), status, ierror)
      call MPI_Wait(request(2), status, ierror)

      return
      end


! ====================================================================
! MPI_Wait_dd : MPI_WAITの4倍精度MPI版（通信方式１）
!	注：送受信を2回行うため、requestを配列で渡す
!
!	request  : in    : 通信リクエスト
!	status   : out   : 送信プロセスの情報
!	ierror   : out   : エラーコード
!
      subroutine MPI_WaitIrecv_dd(
     $                         hi, lo, recvbuf, count,
     $                         request, status, ierror)
      implicit none
      include 'mpif.h'

      integer count
      double precision hi(*)
      double precision lo(*)
      double precision recvbuf(count*2)
      integer request(2)
      integer status(MPI_STATUS_SIZE)
      integer ierror

      call MPI_Wait(request(1), status, ierror)
      call MPI_Wait(request(2), status, ierror)

      hi(1:count) = recvbuf(1:count)
      lo(1:count) = recvbuf(count+1:count*2)

      return
      end

! ====================================================================
! MPI_Reduce_dd : MPI_REDUCEの4倍精度MPI版（通信方式１）
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!   count    : in    : データ数（全プロセスで同じ数を指定）
!   datatype : in    : データ型
!	op       : in    : 演算の種類（加算・乗算・最大・最小のみ）
!	root     : in    : 受信するプロセスのランク
!	comm     : in    : MPIコミュニケータ
!	ierror   : out   : エラーコード
!
      subroutine MPI_Reduce_dd(
     $                      send_hi, send_lo, recv_hi, recv_lo,
     $                      count, datatype, op, root, comm, ierror)
      implicit none
      include 'mpif.h'

      double precision send_hi(*), send_lo(*)
      double precision recv_hi(*), recv_lo(*)
      integer count, datatype, op, root, comm, ierror

      integer nprocs, myrank, ierr

      call MPI_COMM_SIZE(comm, nprocs, ierr)
      call MPI_COMM_RANK(comm, myrank, ierr)

      call MPI_Reduce_dd_func(
     $                 send_hi, send_lo, recv_hi, recv_lo,
     $                 nprocs, myrank, count, datatype, op,
     $                 root, comm, ierror)

!      if ( myrank == root ) then
!          print *, recv_hi(1:count)
!          print *, recv_lo(1:count)
!      end if

! ====================================================================
      contains
      subroutine MPI_Reduce_dd_func(
     $                  send_hi, send_lo, recv_hi, recv_lo,
     $                  nprocs, myrank, count, datatype, op,
     $                  root, comm, ierror)
      implicit none

      double precision send_hi(*), send_lo(*)
      double precision recv_hi(*), recv_lo(*)
      integer nprocs, myrank, count, datatype, op, root, comm, ierror

      double precision recvbuf_hi(nprocs*count)
      double precision recvbuf_lo(nprocs*count)
      double precision hi(nprocs), lo(nprocs)
      integer i

!      print *, 'Reduce_func 1', nprocs, myrank, send_hi(count)

!     rootに全プロセスのデータを送信（Gather）
      call MPI_Gather_dd(send_hi, send_lo, count, datatype,
     $                   recvbuf_hi, recvbuf_lo, count, datatype,
     $                   root, comm, ierror)

      if (myrank == root) then
          do i = 1, count
              hi = recvbuf_hi(i:nprocs*count:count)
              lo = recvbuf_lo(i:nprocs*count:count)
              call reduce_func(
     $                  hi, lo, recv_hi(i), recv_lo(i), nprocs, op)
          end do
      end if

      end subroutine MPI_Reduce_dd_func

      end

! ====================================================================
! MPI_AllReduce_dd : MPI_ALLREDUCEの4倍精度MPI版（通信方式１）
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!   count    : in    : データ数
!   datatype : in    : データ型
!	op       : in    : 演算の種類（加算・乗算・最大・最小のみ）
!	comm     : in    : MPIコミュニケータ
!	ierror   : out   : エラーコード
!
      subroutine MPI_AllReduce_dd(
     $                      send_hi, send_lo, recv_hi, recv_lo,
     $                      count, datatype, op, comm, ierror)
      implicit none

      double precision send_hi(*), send_lo(*)
      double precision recv_hi(*), recv_lo(*)
      integer count, datatype, op, comm, ierror

! rootを0とする
      call MPI_Reduce_dd(
     $                 send_hi, send_lo, recv_hi, recv_lo,
     $                 count, datatype, op, 0, comm, ierror)

! ブロードキャスト
      call MPI_Bcast_dd(
     $              recv_hi, recv_lo, count, datatype, 0, comm, ierror)

      return 
      end

! ====================================================================
! MPI_Alltoall_dd : MPI_ALLTOALLの4倍精度MPI版（通信方式１）
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	scount   : in    : 1つのプロセスに対して送る要素数
!   sendtype : in    : データ型
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	rcount   : in    : 1つのプロセスから受け取る要素数
!   recvtype : in    : データ型
!	comm     : in    : MPIコミュニケータ
!	ierror   : out   : エラーコード
!
      subroutine MPI_Alltoall_dd(
     $                      send_hi, send_lo, scount, sendtype,
     $                      recv_hi, recv_lo, rcount, recvtype,
     $                      comm, ierror)
      implicit none
      include 'mpif.h'

      double precision send_hi(*), send_lo(*)
      double precision recv_hi(*), recv_lo(*)
      integer scount, rcount, sendtype, recvtype, comm, ierror

      call MPI_Alltoall(send_hi, scount, sendtype,
     $                  recv_hi, rcount, recvtype,
     $                  comm, ierror)
      call MPI_Alltoall(send_lo, scount, sendtype,
     $                  recv_lo, rcount, recvtype,
     $                  comm, ierror)

      return
      end

! ====================================================================
! MPI_Alltoallv_dd : MPI_ALLTOALLVの4倍精度MPI版（通信方式１）
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	scounts  : in    : 1つのプロセスに対して送る要素数
!	sdispls  : in    : ランクiのプロセスに送るsend配列内の位置
!   sendtype : in    : データ型
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	rcounts  : in    : 1つのプロセスから受け取る要素数
!	rdispls  : in    : ランクiのプロセスから受信するrecv配列内の位置
!   recvtype : in    : データ型
!	comm     : in    : MPIコミュニケータ
!	ierror   : out   : エラーコード
!
      subroutine MPI_Alltoallv_dd(
     $                  send_hi, send_lo, scounts, sdispls, sendtype,
     $                  recv_hi, recv_lo, rcounts, rdispls, recvtype,
     $                  comm, ierror)
      implicit none
      include 'mpif.h'

      double precision send_hi(*), send_lo(*)
      integer scounts(*), sdispls(*)
      double precision recv_hi(*), recv_lo(*)
      integer rcounts(*), rdispls(*)
      integer sendtype, recvtype, comm, ierror

      call MPI_Alltoallv(
     $         send_hi, scounts, sdispls, sendtype,
     $         recv_hi, rcounts, rdispls, recvtype,
     $         comm, ierror)
      call MPI_Alltoallv(
     $         send_lo, scounts, sdispls, sendtype,
     $         recv_lo, rcounts, rdispls, recvtype,
     $         comm, ierror)

      return
      end

! ====================================================================
! MPI_Bcast_dd : MPI_BCASTの4倍精度MPI版（通信方式１）
!
!	hi       : inout : double-doubleの上位16桁の1次元配列
!	lo       : inout : double-doubleの下位16桁の1次元配列
!	count    : in    : hi, loの配列要素数
!   datatype : in    : データ型
!	root     : in    : 送信元プロセスのrank
!	comm     : in    : MPIコミュニケータ
!	ierror   : out   : エラーコード
!
      subroutine MPI_Bcast_dd(
     $                    hi, lo, count, datatype, root, comm, ierror)

      implicit none
      include 'mpif.h'

      double precision hi(*), lo(*)
      integer count, datatype, root, comm, ierror

      call MPI_Bcast(hi, count, datatype, root, comm, ierror)
      call MPI_Bcast(lo, count, datatype, root, comm, ierror)

      return
      end

! ====================================================================
! MPI_Gather_dd : MPI_GATHERの4倍精度MPI版（通信方式１）
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	scount   : in    : 1つのプロセスに対して送る要素数
!   sendtype : in    : データ型
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	rcount   : in    : 1つのプロセスから受け取る要素数
!   recvtype : in    : データ型
!	root     : in    : 受信するプロセスのrank
!	comm     : in    : MPIコミュニケータ
!	ierror   : out   : エラーコード
!
      subroutine MPI_Gather_dd(
     $                      send_hi, send_lo, scount, sendtype,
     $                      recv_hi, recv_lo, rcount, recvtype,
     $                      root, comm, ierror)
      implicit none
      include 'mpif.h'

      double precision send_hi(*), send_lo(*)
      double precision recv_hi(*), recv_lo(*)
      integer scount, rcount, sendtype, recvtype, root, comm, ierror

      call MPI_Gather(send_hi, scount, sendtype,
     $                recv_hi, rcount, recvtype,
     $                root, comm, ierror)
      call MPI_Gather(send_lo, scount, sendtype,
     $                recv_lo, rcount, recvtype,
     $                root, comm, ierror)

      return
      end

! ====================================================================
! MPI_Gatherv_dd : MPI_GATHERVの4倍精度MPI版（通信方式１）
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	scount   : in    : 1つのプロセスに対して送る要素数
!   sendtype : in    : データ型
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	rcounts  : in    : ランクiのプロセスから受け取る要素数
!	displs   : in    : ランクiのプロセスから受信するrecv配列内の位置
!   recvtype : in    : データ型
!	root     : in    : 受信するプロセスのrank
!	comm     : in    : MPIコミュニケータ
!	ierror   : out   : エラーコード
!
      subroutine MPI_Gatherv_dd(
     $                  send_hi, send_lo, scount, sendtype,
     $                  recv_hi, recv_lo, rcounts, displs, recvtype,
     $                  root, comm, ierror)
      implicit none
      include 'mpif.h'

      double precision send_hi(*), send_lo(*)
      double precision recv_hi(*), recv_lo(*)
      integer rcounts(*), displs(*)
      integer scount, sendtype, recvtype, root, comm, ierror

      call MPI_Gatherv(send_hi, scount, sendtype,
     $                 recv_hi, rcounts, displs, recvtype,
     $                 root, comm, ierror)
      call MPI_Gatherv(send_lo, scount, sendtype,
     $                 recv_lo, rcounts, displs, recvtype,
     $                 root, comm, ierror)

      return
      end

! ====================================================================
! MPI_AllGather_dd : MPI_ALLGATHERの4倍精度MPI版（通信方式１）
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	scount   : in    : 1つのプロセスに対して送る要素数
!   sendtype : in    : データ型
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	rcount   : in    : 1つのプロセスから受け取る要素数
!   recvtype : in    : データ型
!	comm     : in    : MPIコミュニケータ
!	ierror   : out   : エラーコード
!
      subroutine MPI_AllGather_dd(
     $                      send_hi, send_lo, scount, sendtype,
     $                      recv_hi, recv_lo, rcount, recvtype,
     $                      comm, ierror)
      implicit none
      include 'mpif.h'

      double precision send_hi(*), send_lo(*)
      double precision recv_hi(*), recv_lo(*)
      integer scount, rcount, sendtype, recvtype, comm, ierror

      call MPI_AllGather(send_hi, scount, sendtype,
     $                   recv_hi, rcount, recvtype,
     $                   comm, ierror)
      call MPI_AllGather(send_lo, scount, sendtype,
     $                   recv_lo, rcount, recvtype,
     $                   comm, ierror)

      return
      end

! ====================================================================
! MPI_AllGatherv_dd : MPI_ALLGATHERVの4倍精度MPI版（通信方式１）
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	scount   : in    : 1つのプロセスに対して送る要素数
!   sendtype : in    : データ型
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	rcounts  : in    : ランクiのプロセスから受け取る要素数
!	displs   : in    : ランクiのプロセスから受け取る要素の配備位置
!   recvtype : in    : データ型
!	comm     : in    : MPIコミュニケータ
!	ierror   : out   : エラーコード
!
      subroutine MPI_AllGatherv_dd(
     $                      send_hi, send_lo, scount, sendtype, 
     $                      recv_hi, recv_lo, rcounts, displs, recvtype,
     $                      comm, ierror)
      implicit none
      include 'mpif.h'

      double precision send_hi(*), send_lo(*)
      double precision recv_hi(*), recv_lo(*)
      integer rcounts(*), displs(*)
      integer scount, sendtype, recvtype, comm, ierror

      call MPI_AllGatherv(
     $                 send_hi, scount, sendtype,
     $                 recv_hi, rcounts, displs, recvtype,
     $                 comm, ierror)
      call MPI_AllGatherv(
     $                 send_lo, scount, sendtype,
     $                 recv_lo, rcounts, displs, recvtype,
     $                 comm, ierror)

      return
      end

! ====================================================================
! MPI_Scatter_dd : MPI_SCATTERの4倍精度MPI版（通信方式１）
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	scount   : in    : 1つのプロセスに対して送る要素数
!   sendtype : in    : データ型
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	rcount   : in    : 1つのプロセスから受け取る要素数
!   recvtype : in    : データ型
!	root     : in    : 送信元プロセスのrank
!	comm     : in    : MPIコミュニケータ
!	ierror   : out   : エラーコード
!
      subroutine MPI_Scatter_dd(
     $                      send_hi, send_lo, scount, sendtype,
     $                      recv_hi, recv_lo, rcount, recvtype,
     $                      root, comm, ierror)
      implicit none
      include 'mpif.h'

      double precision send_hi(*), send_lo(*)
      double precision recv_hi(*), recv_lo(*)
      integer scount, rcount, sendtype, recvtype, root, comm, ierror

      call MPI_Scatter(send_hi, scount, sendtype,
     $                 recv_hi, rcount, recvtype,
     $                 root, comm, ierror)
      call MPI_Scatter(send_lo, scount, sendtype,
     $                 recv_lo, rcount, recvtype,
     $                 root, comm, ierror)

      return
      end

! ====================================================================
! MPI_Scatterv_dd : MPI_SCATTERVの4倍精度MPI版（通信方式１）
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	scounts  : in    : 1つのプロセスに対して送る要素数
!	displs   : in    : ランクiのプロセスに送るsend配列内の位置
!   sendtype : in    : データ型
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	rcount   : in    : ランクiのプロセスから受け取る要素数
!   recvtype : in    : データ型
!	root     : in    : 送信元プロセスのrank
!	comm     : in    : MPIコミュニケータ
!	ierror   : out   : エラーコード
!
      subroutine MPI_Scatterv_dd(
     $                  send_hi, send_lo, scounts, displs, sendtype,
     $                  recv_hi, recv_lo, rcount, recvtype,
     $                  root, comm, ierror)
      implicit none
      include 'mpif.h'

      double precision send_hi(*), send_lo(*)
      integer scounts(*), displs(*)
      double precision recv_hi(*), recv_lo(*)
      integer rcount, sendtype, recvtype, root, comm, ierror

      call MPI_Scatterv(send_hi, scounts, displs, sendtype,
     $                  recv_hi, rcount, recvtype,
     $                  root, comm, ierror)
      call MPI_Scatterv(send_lo, scounts, displs, sendtype,
     $                  recv_lo, rcount, recvtype,
     $                  root, comm, ierror)

      return
      end

