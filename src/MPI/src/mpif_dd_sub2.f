!
! mpif_dd_sub2.f
!	4倍精度MPI通信ライブラリ・通信方式２
!
! Attention:
!   mpif_dd.fファイルのサブルーチンから利用します
!
! Version : 1.0 (2011.1.31)
! Author: IMS
!
! Description:
!	4倍精度データの上位データと下位データを
!	ひとつの配列にまとめて送受信する通信方式
!
! =========================================================================

! ====================================================================
! MPI_Send_dd : MPI_SENDの4倍精度MPI版（通信方式２）
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
     $                      hi, lo, count, datatype, 
     $                      dest, tag, comm, ierror)
      implicit none
      include 'mpif.h'

      integer count, datatype, dest, tag, comm, ierror
      double precision hi(*), lo(*)
      double precision sdata(count*2)
      
      sdata(1:count*2:2) = hi(1:count)
      sdata(2:count*2:2) = lo(1:count)

      call MPI_Send(sdata, count*2, datatype, dest, tag, comm, ierror)

      return
      end

! ====================================================================
! MPI_Recv_dd : MPI_RECVの4倍精度MPI版（通信方式２）
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
      subroutine MPI_Recv_dd(hi, lo, count, datatype, source, tag,
     $                       comm, status, ierror)
      implicit none
      include 'mpif.h'

      double precision hi(*), lo(*)
      integer status(MPI_STATUS_SIZE)
      integer count, datatype, source, tag, comm, ierror

      double precision rdata(count*2)

      call MPI_Recv(rdata, count*2, datatype, source, tag,
     $              comm, status, ierror)

      hi(1:count) = rdata(1:count*2:2)
      lo(1:count) = rdata(2:count*2:2)

      return
      end

! ====================================================================
! MPI_Isend_dd : MPI_ISENDの4倍精度MPI版（通信方式２）
!	注：通信方式２とそろえるため、本来は不要だがrequestが配列で戻る
!
!	hi       : in    : double-doubleの上位16桁の1次元配列
!	lo       : in    : double-doubleの下位16桁の1次元配列
!	sdata    : in    : double-double型配列
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

      double precision hi(*), lo(*)
      integer request(2)
      integer count, datatype, dest, tag, comm, ierror

      double precision sdata(count*2)
      integer ireq
      
      sdata(1:count*2:2) = hi(1:count)
      sdata(2:count*2:2) = lo(1:count)

      call MPI_Isend(sdata, count*2, datatype, dest, tag, 
     $               comm, ireq, ierror)

      request(1) = ireq

      return
      end

! ====================================================================
! MPI_Wait_dd : MPI_WAITの4倍精度MPI版（通信方式２）
!	注：通信方式２とそろえるため、本来は不要だがrequestを配列で渡す
!
!       sdata    : in    : double-double型配列（waitまでの保護のため）
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

      return
      end

! ====================================================================
! MPI_Irecv_dd : MPI_IRECVの4倍精度MPI版（通信方式２）
!	注：通信方式２とそろえるため、本来は不要だがrequestが配列で戻る
!
!	hi       : out   : double-doubleの上位16桁の1次元配列
!	lo       : out   : double-doubleの下位16桁の1次元配列
!	recvbuf  : out   : double-double型配列
!	count    : in    : hi, loの配列要素数
!	datatype : in    : データ型
!	source   : in    : 送信元プロセスのrank
!	tag      : in    : タグ（メッセージ識別番号）
!	comm     : in    : MPIコミュニケータ
!	request  : out   : 通信リクエスト
!	ierror   : out   : エラーコード
!
!      subroutine MPI_Irecv_dd(
!     $                      hi, lo, recvbuf, count, datatype, source, tag,
!     $                      comm, request, ierror)
      subroutine MPI_Irecv_dd(
     $                      recvbuf, count, datatype, source, tag,
     $                      comm, request, ierror)
      implicit none
      include 'mpif.h'

!      double precision hi(*), lo(*)
      integer count, datatype, source, tag, comm, ierror
      integer request(2)
      double precision recvbuf(count*2)
      integer ireq_hi

      call MPI_Irecv(recvbuf, count*2, datatype, source, tag,
     $               comm, ireq_hi, ierror)
      
      request(1) = ireq_hi

      return
      end

! ====================================================================
! MPI_Wait_dd : MPI_WAITの4倍精度MPI版（通信方式２）
!	注：通信方式２とそろえるため、本来は不要だがrequestを配列で渡す
!
!	request  : in    : 通信リクエスト
!	status   : out   : 送信プロセスの情報
!	ierror   : out   : エラーコード
!
      subroutine MPI_WaitIrecv_dd(
     $                hi, lo, rdata, count,
     $                request, status, ierror)

      implicit none
      include 'mpif.h'

      integer count
      double precision hi(*), lo(*)
      double precision rdata(count*2)
      integer request(2)
      integer status(MPI_STATUS_SIZE)
      integer ierror

      call MPI_Wait(request(1), status, ierror)
!      call MPI_Wait(request(2), status, ierror)

      hi(1:count) = rdata(1:count*2:2)
      lo(1:count) = rdata(2:count*2:2)

      return
      end

! ====================================================================
! MPI_Reduce_dd : MPI_REDUCEの4倍精度MPI版（通信方式２）
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	count    : in    : 送受信するデータ数
!	datatype : in    : データ型
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
     $                  recvbuf_hi, recvbuf_lo, count, datatype,
     $                  root, comm, ierror)

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
! MPI_AllReduce_dd : MPI_ALLREDUCEの4倍精度MPI版（通信方式２）
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	count    : in    : 送受信するデータ数
!	datatype : in    : データ型
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
      call MPI_Bcast_dd(recv_hi, recv_lo, count, datatype,
     $                  0, comm, ierror)

      return 
      end

! ====================================================================
! MPI_Alltoall_dd : MPI_ALLTOALLの4倍精度MPI版（通信方式２）
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	scount   : in    : 1つのプロセスに対して送る要素数
!	sendtype : in    : データ型
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	rcount   : in    : 1つのプロセスから受け取る要素数
!	recvtype : in    : データ型
!	comm     : in    : MPIコミュニケータ
!	ierror   : out   : エラーコード
!
      subroutine MPI_Alltoall_dd(
     $                  send_hi, send_lo, scount, sendtype,
     $                  recv_hi, recv_lo, rcount, recvtype,
     $                  comm, ierror)
      implicit none
      include 'mpif.h'

      double precision send_hi(*), send_lo(*)
      double precision recv_hi(*), recv_lo(*)
      integer scount, sendtype
      integer rcount, recvtype
      integer comm, ierror

      integer sendsize, recvsize
      integer nprocs, ierr

      call MPI_COMM_SIZE(comm, nprocs, ierr)

      sendsize = nprocs * scount
      recvsize = nprocs * rcount

      call MPI_Alltoall_dd_func(
     $                  send_hi, send_lo, scount, sendtype, sendsize,
     $                  recv_hi, recv_lo, rcount, recvtype, recvsize,
     $                  comm, ierror)

! ====================================================================
! MPI_Alltoall_dd_func : MPI_Alltoall_ddの内部サブルーチン
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	scount   : in    : 1つのプロセスに対して送る要素数
!	sendtype : in    : データ型
!	sendsize : in    : send_hi, send_loのサイズ
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	rcount   : in    : 1つのプロセスから受け取る要素数
!	recvtype : in    : データ型
!	recvsize : in    : recv_hi, recv_loのサイズ
!	comm     : in    : MPIコミュニケータ
!	ierror   : out   : エラーコード
!

      contains
      subroutine MPI_Alltoall_dd_func(
     $                  send_hi, send_lo, scount, sendtype, sendsize,
     $                  recv_hi, recv_lo, rcount, recvtype, recvsize,
     $                  comm, ierror)
      implicit none
      include 'mpif.h'

      double precision send_hi(*), send_lo(*)
      double precision recv_hi(*), recv_lo(*)
      integer scount, sendtype, sendsize
      integer rcount, recvtype, recvsize
      integer comm, ierror
      integer i

      double precision sendbuf(sendsize*2)
      double precision recvbuf(recvsize*2)

      sendbuf(1:sendsize*2:2) = send_hi(1:sendsize)
      sendbuf(2:sendsize*2:2) = send_lo(1:sendsize)

      call MPI_Alltoall(sendbuf, scount*2, sendtype,
     $                  recvbuf, rcount*2, recvtype,
     $                  comm, ierror)

      recv_hi(1:recvsize) = recvbuf(1:recvsize*2:2)
      recv_lo(1:recvsize) = recvbuf(2:recvsize*2:2)

      end subroutine MPI_Alltoall_dd_func

      end

! ====================================================================
! MPI_Alltoallv_dd : MPI_ALLTOALLVの4倍精度MPI版（通信方式２）
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	scounts  : in    : 1つのプロセスに対して送る要素数
!	sdispls  : in    : ランクiのプロセスに送るsend配列内の位置
!	sendtype : in    : データ型
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	rcounts  : in    : 1つのプロセスから受け取る要素数
!	rdispls  : in    : ランクiのプロセスから受信するrecv配列内の位置
!	recvtype : in    : データ型
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
      integer sdispls(*), scounts(*)
      double precision recv_hi(*), recv_lo(*)
      integer rdispls(*), rcounts(*)
      integer sendtype, recvtype, comm, ierror

      integer i
      integer nprocs, ierr
      integer sendsize, recvsize

      call MPI_COMM_SIZE(comm, nprocs, ierr)

      sendsize = 0
      recvsize = 0
      do i = 1, nprocs
          sendsize = sendsize + scounts(i)
          recvsize = recvsize + rcounts(i)
      end do

      call MPI_Alltoallv_dd_func(
     $           send_hi, send_lo, scounts, sdispls, sendtype, sendsize,
     $           recv_hi, recv_lo, rcounts, rdispls, recvtype, recvsize,
     $           nprocs, comm, ierror)

! ====================================================================
! MPI_AlltoallV_dd_func : MPI_AlltoallV_ddの内部プログラム
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	scounts  : in    : 1つのプロセスに対して送る要素数
!	sdispls  : in    : ランクiのプロセスに送るsend配列内の位置
!	sendtype : in    : データ型
!	sendsize : in    : send_hi, send_loのサイズ
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	rcounts  : in    : 1つのプロセスから受け取る要素数
!	rdispls  : in    : ランクiのプロセスから受信するrecv配列内の位置
!	recvtype : in    : データ型
!	recvsize : in    : recv_hi, recv_loのサイズ
!	nprocs   : in    : プロセス数
!	comm     : in    : MPIコミュニケータ
!	ierror   : out   : エラーコード
!
      contains
      subroutine MPI_Alltoallv_dd_func(
     $           send_hi, send_lo, scounts, sdispls, sendtype, sendsize,
     $           recv_hi, recv_lo, rcounts, rdispls, recvtype, recvsize,
     $           nprocs, comm, ierror)
      implicit none
      include 'mpif.h'

      double precision send_hi(*), send_lo(*)
      integer sdispls(*), scounts(*)
      double precision recv_hi(*), recv_lo(*)
      integer rdispls(*), rcounts(*)
      integer sendtype, sendsize, smin
      integer recvtype, recvsize, rmin
      integer nprocs, comm, ierror

      double precision sendbuf(sendsize*2)
      integer sdispls2(nprocs), scounts2(nprocs)
      double precision recvbuf(recvsize*2)
      integer rdispls2(nprocs), rcounts2(nprocs)

!      sdispls2(1:nprocs) = sdispls(1:nprocs) * 2
!      rdispls2(1:nprocs) = rdispls(1:nprocs) * 2
      smin = sdispls(1)
      rmin = rdispls(1)
      do i = 2, nprocs
          if (smin .gt. sdispls(i)) then
              smin = sdispls(i)
          end if
          if (rmin .gt. rdispls(i)) then
              rmin = rdispls(i)
          end if
      end do
      do i = 1, nprocs
          sdispls2(i) = (sdispls(i) - smin) * 2
          rdispls2(i) = (rdispls(i) - rmin) * 2
      end do
      scounts2(1:nprocs) = scounts(1:nprocs) * 2
      rcounts2(1:nprocs) = rcounts(1:nprocs) * 2

      sendbuf(1:sendsize*2:2) = send_hi((1+smin):(sendsize+smin))
      sendbuf(2:sendsize*2:2) = send_lo((1+smin):(sendsize+smin))
!      sendbuf(1:sendsize*2:2) = send_hi(1:sendsize)
!      sendbuf(2:sendsize*2:2) = send_lo(1:sendsize)

      call MPI_Alltoallv(
     $         sendbuf, scounts2, sdispls2, sendtype,
     $         recvbuf, rcounts2, rdispls2, recvtype,
     $         comm, ierror)

      recv_hi((1+rmin):(recvsize+rmin)) = recvbuf(1:recvsize*2:2)
      recv_lo((1+rmin):(recvsize+rmin)) = recvbuf(2:recvsize*2:2)
!      recv_hi(1:recvsize) = recvbuf(1:recvsize*2:2)
!      recv_lo(1:recvsize) = recvbuf(2:recvsize*2:2)

      end subroutine MPI_Alltoallv_dd_func

      end

! ====================================================================
! MPI_Bcast_dd : MPI_BCASTの4倍精度MPI版（通信方式２）
!
!	hi       : inout : double-doubleの上位16桁の1次元配列
!	lo       : inout : double-doubleの下位16桁の1次元配列
!	count    : in    : hi, loの配列要素数
!	datatype : in    : データ型
!	root     : in    : 送信元プロセスのrank
!	comm     : in    : MPIコミュニケータ
!	ierror   : out   : エラーコード
!
      subroutine MPI_Bcast_dd(
     $                      hi, lo, count, datatype, root, comm, ierror)

      implicit none
      include 'mpif.h'

      double precision hi(*), lo(*)
      integer count, datatype, root, comm, ierror
      integer myrank, ierr
      double precision castdata(count*2)

      call MPI_COMM_RANK(comm, myrank, ierr)

      if ( myrank == root ) then
          castdata(1:count*2:2) = hi(1:count)
          castdata(2:count*2:2) = lo(1:count)
      end if

      call MPI_Bcast(castdata, count*2, datatype, root, comm, ierror)

      if ( myrank .ne. root ) then
          hi(1:count) = castdata(1:count*2:2)
          lo(1:count) = castdata(2:count*2:2)
      end if

      return
      end

! ====================================================================
! MPI_Gather_dd : MPI_GATHERの4倍精度MPI版（通信方式２）
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	scount   : in    : 1つのプロセスに対して送る要素数
!	sendtype : in    : データ型
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	rcount   : in    : 1つのプロセスから受け取る要素数
!	recvtype : in    : データ型
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
      integer scount, sendtype
      integer rcount, recvtype
      integer root, comm, ierror

      integer nprocs, ierr

      call MPI_COMM_SIZE(comm, nprocs, ierr)

      call MPI_Gather_dd_func(send_hi, send_lo, scount, sendtype,
     $                        recv_hi, recv_lo, rcount, recvtype,
     $                        nprocs, root, comm, ierror)

! ====================================================================
! MPI_Gather_dd_func : MPI_Gather_ddの内部サブルーチン
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	scount   : in    : 1つのプロセスに対して送る要素数
!	sendtype : in    : データ型
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	rcount   : in    : 1つのプロセスから受け取る要素数
!	recvtype : in    : データ型
!	nprocs   : in    : commのプロセス数
!	root     : in    : 受信するプロセスのrank
!	comm     : in    : MPIコミュニケータ
!	ierror   : out   : エラーコード
!
      contains
      subroutine MPI_Gather_dd_func(
     $                      send_hi, send_lo, scount, sendtype,
     $                      recv_hi, recv_lo, rcount, recvtype,
     $                      nprocs, root, comm, ierror)
      implicit none
      include 'mpif.h'

      double precision send_hi(*), send_lo(*)
      double precision recv_hi(*), recv_lo(*)
      integer scount, sendtype
      integer rcount, recvtype
      integer root, comm, ierror

      integer myrank, nprocs, ierr
      double precision sendbuf(scount*2)
      double precision recvbuf(nprocs*rcount*2)

      call MPI_COMM_RANK(comm, myrank, ierr)

      sendbuf(1:scount*2:2) = send_hi(1:scount)
      sendbuf(2:scount*2:2) = send_lo(1:scount)

      call MPI_Gather(sendbuf, scount*2, sendtype,
     $                recvbuf, rcount*2, recvtype,
     $                root, comm, ierror)

      if (myrank == root) then
          recv_hi(1:nprocs*rcount) = recvbuf(1:nprocs*rcount*2:2)
          recv_lo(1:nprocs*rcount) = recvbuf(2:nprocs*rcount*2:2)
      end if
      end subroutine MPI_Gather_dd_func

      end

! ====================================================================
! MPI_Gatherv_dd : MPI_GATHERVの4倍精度MPI版（通信方式２）
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	scount   : in    : 1つのプロセスに対して送る要素数
!	sendtype : in    : データ型
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	rcounts  : in    : ランクiのプロセスから受け取る要素数
!	displs   : in    : ランクiのプロセスから受信するrecv配列内の位置
!	recvtype : in    : データ型
!	root     : in    : 受信するプロセスのrank
!	comm     : in    : MPIコミュニケータ
!	ierror   : out   : エラーコード
!
      subroutine MPI_Gatherv_dd(
     $                  send_hi, send_lo, scount, sendtype,
     $                  recv_hi, recv_lo, rcounts, displs, recvtype,
     $                  root, comm, ierror)
      implicit none

      double precision send_hi(*), send_lo(*)
      double precision recv_hi(*), recv_lo(*)
      integer rcounts(*), displs(*)
      integer scount, sendtype, recvtype
      integer root, comm, ierror

      integer nprocs, ierr, recvsize
      integer i

      call MPI_COMM_SIZE(comm, nprocs, ierr)

      recvsize = 0
      do i = 1, nprocs
          recvsize = recvsize + rcounts(i)
      end do

      call MPI_Gatherv_dd_func(
     $            send_hi, send_lo, scount, sendtype,
     $            recv_hi, recv_lo, rcounts, displs, recvtype, recvsize,
     $            nprocs, root, comm, ierror)

! ====================================================================
! MPI_Gatherv_dd_func : MPI_Gathervの内部サブルーチン
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	scount   : in    : 1つのプロセスに対して送る要素数
!	sendtype : in    : データ型
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	rcounts  : in    : ランクiのプロセスから受け取る要素数
!	displs   : in    : ランクiのプロセスから受信するrecv配列内の位置
!	recvtype : in    : データ型
!	recvsize : in    : recv_hi, recv_loのサイズ
!	nprocs   : in    : commのプロセス数
!	root     : in    : 受信するプロセスのrank
!	comm     : in    : MPIコミュニケータ
!	ierror   : out   : エラーコード
!
      contains
      subroutine MPI_Gatherv_dd_func(
     $          send_hi, send_lo, scount, sendtype,
     $          recv_hi, recv_lo, rcounts, displs, recvtype, recvsize,
     $          nprocs, root, comm, ierror)
      implicit none
      include 'mpif.h'

      double precision send_hi(*), send_lo(*)
      double precision recv_hi(*), recv_lo(*)
      integer rcounts(*), displs(*)
      integer scount, sendtype, recvtype, recvsize
      integer nprocs, root, comm, ierror

      double precision sendbuf(scount*2)
      double precision recvbuf(recvsize*2)
      integer rcounts2(nprocs)
      integer displs2(nprocs)
      integer rmin, myrank, ierr

      call MPI_COMM_RANK(comm, myrank, ierr)

!      displs2(1:nprocs) = displs(1:nprocs) * 2
      rmin = displs(1)
      do i = 2, nprocs
          if (rmin .gt. displs(i)) then
              rmin = displs(i)
          end if
      end do
      do i = 1, nprocs
          displs2(i) = (displs(i) - rmin) * 2
      end do
      rcounts2(1:nprocs) = rcounts(1:nprocs) * 2

      sendbuf(1:scount*2:2) = send_hi(1:scount)
      sendbuf(2:scount*2:2) = send_lo(1:scount)

      call MPI_Gatherv(
     $         sendbuf, scount*2, sendtype,
     $         recvbuf, rcounts2, displs2, recvtype,
     $         root, comm, ierror)

      if (myrank == root) then
          recv_hi((1+rmin):(recvsize+rmin)) = recvbuf(1:recvsize*2:2)
          recv_lo((1+rmin):(recvsize+rmin)) = recvbuf(2:recvsize*2:2)
!          recv_hi(1:recvsize) = recvbuf(1:recvsize*2:2)
!          recv_lo(1:recvsize) = recvbuf(2:recvsize*2:2)
      end if

      end subroutine MPI_Gatherv_dd_func

      end

! ====================================================================
! MPI_AllGather_dd : MPI_ALLGATHERの4倍精度MPI版（通信方式２）
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	scount   : in    : 1つのプロセスに対して送る要素数
!	sendtype : in    : データ型
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	rcount   : in    : 1つのプロセスから受け取る要素数
!	recvtype : in    : データ型
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
      integer scount, sendtype
      integer rcount, recvtype
      integer comm, ierror
      integer nprocs, ierr

      call MPI_COMM_SIZE(comm, nprocs, ierr)

      call MPI_AllGather_dd_func(send_hi, send_lo, scount, sendtype,
     $                           recv_hi, recv_lo, rcount, recvtype,
     $                           nprocs, comm, ierror)
      
! ====================================================================
! MPI_AllGather_dd_func : MPI_AllGather_ddの実行
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	scount   : in    : 1つのプロセスに対して送る要素数
!	sendtype : in    : データ型
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	rcount   : in    : 1つのプロセスから受け取る要素数
!	recvtype : in    : データ型
!	comm     : in    : MPIコミュニケータ
!	ierror   : out   : エラーコード
!
      contains
      subroutine MPI_AllGather_dd_func(
     $                      send_hi, send_lo, scount, sendtype,
     $                      recv_hi, recv_lo, rcount, recvtype,
     $                      nprocs, comm, ierror)
      implicit none

      double precision send_hi(*), send_lo(*)
      double precision recv_hi(*), recv_lo(*)
      integer scount, sendtype
      integer rcount, recvtype
      integer nprocs, comm, ierror

      double precision sendbuf(scount*2)
      double precision recvbuf(nprocs*rcount*2)

      sendbuf(1:scount*2:2) = send_hi(1:scount)
      sendbuf(2:scount*2:2) = send_lo(1:scount)

      call MPI_AllGather(sendbuf, scount*2, sendtype,
     $                   recvbuf, rcount*2, recvtype,
     $                   comm, ierror)

      recv_hi(1:nprocs*rcount) = recvbuf(1:nprocs*rcount*2:2)
      recv_lo(1:nprocs*rcount) = recvbuf(2:nprocs*rcount*2:2)

      end subroutine MPI_AllGather_dd_func

      end

! ====================================================================
! MPI_AllGatherv_dd : MPI_ALLGATHERVの4倍精度MPI版（通信方式２）
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	scount   : in    : 1つのプロセスに対して送る要素数
!	sendtype : in    : データ型
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	rcounts  : in    : ランクiのプロセスから受け取る要素数
!	displs   : in    : ランクiのプロセスから受け取るデータの配備位置
!	recvtype : in    : データ型
!	comm     : in    : MPIコミュニケータ
!	ierror   : out   : エラーコード
!
      subroutine MPI_AllGatherv_dd(
     $                  send_hi, send_lo, scount, sendtype,
     $                  recv_hi, recv_lo, rcounts, displs, recvtype,
     $                  comm, ierror)
      implicit none
      include 'mpif.h'

      double precision send_hi(*), send_lo(*)
      double precision recv_hi(*), recv_lo(*)
      integer rcounts(*), displs(*)
      integer scount, sendtype, recvtype
      integer comm, ierror

      integer nprocs, ierr, recvsize
      integer i

      call MPI_COMM_SIZE(comm, nprocs, ierr)

      recvsize = 0
      do i = 1, nprocs
          recvsize = recvsize + rcounts(i)
      end do

      call MPI_AllGatherv_dd_func(
     $          send_hi, send_lo, scount, sendtype,
     $          recv_hi, recv_lo, rcounts, displs, recvtype, recvsize,
     $          nprocs, comm, ierror)

! ====================================================================
! MPI_AllGatherv_dd_func : MPI_AllGatherv_ddの内部サブルーチン
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	scount   : in    : 1つのプロセスに対して送る要素数
!	sendtype : in    : データ型
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	rcounts  : in    : ランクiのプロセスから受け取る要素数
!	displs   : in    : ランクiのプロセスから受け取るデータの配備位置
!	recvtype : in    : データ型
!	comm     : in    : MPIコミュニケータ
!	ierror   : out   : エラーコード
!
      contains
      subroutine MPI_AllGatherv_dd_func(
     $          send_hi, send_lo, scount, sendtype,
     $          recv_hi, recv_lo, rcounts, displs, recvtype, recvsize,
     $          nprocs, comm, ierror)
      implicit none

      double precision send_hi(*), send_lo(*)
      double precision recv_hi(*), recv_lo(*)
      integer rcounts(*), displs(*)
      integer scount, sendtype, recvtype, recvsize
      integer nprocs, comm, ierror

      double precision sendbuf(scount*2)
      double precision recvbuf(recvsize*2)
      integer rcounts2(nprocs)
      integer displs2(nprocs)
      integer rmin

!      displs2(1:nprocs) = displs(1:nprocs) * 2
      rmin = displs(1)
      do i = 2, nprocs
          if (rmin .gt. displs(i)) then
              rmin = displs(i)
          end if
      end do
      do i = 1, nprocs
          displs2(i) = (displs(i) - rmin) * 2
      end do
      rcounts2(1:nprocs) = rcounts(1:nprocs) * 2

      sendbuf(1:scount*2:2) = send_hi(1:scount)
      sendbuf(2:scount*2:2) = send_lo(1:scount)

      call MPI_AllGatherv(
     $                 sendbuf, scount*2, sendtype,
     $                 recvbuf, rcounts2, displs2, recvtype,
     $                 comm, ierror)

      recv_hi((1+rmin):(recvsize+rmin)) = recvbuf(1:recvsize*2:2)
      recv_lo((1+rmin):(recvsize+rmin)) = recvbuf(2:recvsize*2:2)
!      recv_hi(1:recvsize) = recvbuf(1:recvsize*2:2)
!      recv_lo(1:recvsize) = recvbuf(2:recvsize*2:2)

      end subroutine MPI_AllGatherv_dd_func

      end

! ====================================================================
! MPI_Scatter_dd : MPI_SCATTERの4倍精度MPI版（通信方式２）
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	scount   : in    : 1つのプロセスに対して送る要素数
!	sendtype : in    : データ型
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	rcount   : in    : 1つのプロセスから受け取る要素数
!	recvtype : in    : データ型
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
      integer scount, sendtype
      integer rcount, recvtype
      integer root, comm, ierror

      integer nprocs, myrank, ierr

      call MPI_COMM_SIZE(comm, nprocs, ierr)
      call MPI_COMM_RANK(comm, myrank, ierr)

      call MPI_Scatter_dd_func(
     $                     send_hi, send_lo, scount, sendtype,
     $                     recv_hi, recv_lo, rcount, recvtype,
     $                     nprocs, myrank, root, comm, ierror)

! ====================================================================
      contains
      subroutine MPI_Scatter_dd_func(
     $                     send_hi, send_lo, scount, sendtype,
     $                     recv_hi, recv_lo, rcount, recvtype,
     $                     nprocs, myrank, root, comm, ierror)
      implicit none

      integer root, comm, ierror
      double precision send_hi(*), send_lo(*)
      double precision recv_hi(*), recv_lo(*)
      integer scount, sendtype
      integer rcount, recvtype
      integer nprocs, myrank

      double precision sendbuf(nprocs*scount*2)
      double precision recvbuf(rcount*2)

      if ( root == myrank ) then
          sendbuf(1:nprocs*scount*2:2) = send_hi(1:nprocs*scount)
          sendbuf(2:nprocs*scount*2:2) = send_lo(1:nprocs*scount)
      end if

      call MPI_Scatter(sendbuf, scount*2, sendtype,
     $                 recvbuf, rcount*2, recvtype,
     $                 root, comm, ierror)

      recv_hi(1:rcount) = recvbuf(1:rcount*2:2)
      recv_lo(1:rcount) = recvbuf(2:rcount*2:2)

      end subroutine MPI_Scatter_dd_func

      end

! ====================================================================
! MPI_ScatterV_dd : MPI_SCATTERVの4倍精度MPI版（通信方式２）
!
!	send_hi  : in    : 送信するdouble-doubleの上位16桁の1次元配列
!	send_lo  : in    : 送信するdouble-doubleの下位16桁の1次元配列
!	scounts  : in    : 1つのプロセスに対して送る要素数
!	displs   : in    : ランクiのプロセスに送るsend配列内の位置
!	recv_hi  : out   : 受信するdouble-doubleの上位16桁の1次元配列
!	recv_lo  : out   : 受信するdouble-doubleの下位16桁の1次元配列
!	rcount   : in    : 受け取る要素数
!	root     : in    : 送信元プロセスのrank
!	comm     : in    : MPIコミュニケータ
!	ierror   : out   : エラーコード
!
      subroutine MPI_Scatterv_dd(
     $              send_hi, send_lo, scounts, displs, sendtype,
     $              recv_hi, recv_lo, rcount, recvtype,
     $              root, comm, ierror)
      implicit none
      include 'mpif.h'

      double precision send_hi(*), send_lo(*)
      integer scounts(*)
      integer displs(*)
      double precision recv_hi(*), recv_lo(*)
      integer sendtype, rcount, recvtype
      integer root, comm, ierror

      integer nprocs, myrank, ierr
      integer sendsize
      integer i

      call MPI_COMM_SIZE(comm, nprocs, ierr)
      call MPI_COMM_RANK(comm, myrank, ierr)

      sendsize = 1
      if ( root == myrank ) then
          sendsize = 0
          do i = 1, nprocs
              sendsize = sendsize + scounts(i)
          end do
      end if

      call MPI_Scatterv_dd_func(
     $            send_hi, send_lo, scounts, displs, sendtype, sendsize,
     $            recv_hi, recv_lo, rcount, recvtype,
     $            nprocs, myrank, root, comm, ierror)

! ====================================================================
      contains
      subroutine MPI_Scatterv_dd_func(
     $            send_hi, send_lo, scounts, displs, sendtype, sendsize,
     $            recv_hi, recv_lo, rcount, recvtype,
     $            nprocs, myrank, root, comm, ierror)
      implicit none

      integer rcount, sendsize, sendtype, recvtype
      double precision send_hi(sendsize), send_lo(sendsize)
      integer scounts(nprocs), displs(nprocs)
      double precision recv_hi(rcount), recv_lo(rcount)
      integer nprocs, myrank, root, comm, ierror

      integer scounts2(nprocs)
      integer displs2(nprocs)
      integer smin
      double precision sendbuf(sendsize*2)
      double precision recvbuf(rcount*2)

      if ( myrank == root ) then
          scounts2 = scounts * 2
          smin = displs(1)
          do i = 2, nprocs
              if (smin .gt. displs(i)) then
                  smin = displs(i)
              end if
          end do
          do i = 1, nprocs
              displs2(i) = (displs(i) - smin) * 2
          end do
          sendbuf(1:sendsize*2:2) = send_hi((1+smin):(sendsize+smin))
          sendbuf(2:sendsize*2:2) = send_lo((1+smin):(sendsize+smin))
!          sendbuf(1:sendsize*2:2) = send_hi(1:sendsize)
!          sendbuf(2:sendsize*2:2) = send_lo(1:sendsize)
!          displs2 = displs * 2
      end if

      call MPI_Scatterv(
     $              sendbuf, scounts2, displs2, sendtype,
     $              recvbuf, rcount*2, recvtype,
     $              root, comm, ierror)

      recv_hi(1:rcount) = recvbuf(1:rcount*2:2)
      recv_lo(1:rcount) = recvbuf(2:rcount*2:2)

      end subroutine MPI_Scatterv_dd_func

      end

