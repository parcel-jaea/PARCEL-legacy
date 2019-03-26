!
! mpif_dd_func.f
!      double-double型の演算
!
!
! ======================================================================

! ======================================================================
! 加算
!
! Author:
!      David H. Bailey
!      NERSC, Lawrence Berkeley Lab
!      （F77形式に変換して使用）
!
      subroutine ddadd(org1_hi, org1_lo, org2_hi, org2_lo, hi, lo) 

      implicit none
      double precision org1_hi, org1_lo
      double precision org2_hi, org2_lo
      double precision hi, lo
      double precision e, t1, t2

!   Compute dda + ddb using Knuth's trick.

!      write (*, "('add func 1 ',E24.16,E24.16,E24.16,E24.16)") 
!     $        org1_hi,org1_lo,org2_hi,org2_lo

      t1 = org1_hi + org2_hi
      e = t1 - org1_hi
      t2 = ((org2_hi - e) + (org1_hi - (t1 - e))) + org1_lo + org2_lo

!      write (*, "('add func 2 ',E24.16,E24.16,E24.16)") t1, t2, e

!   The result is t1 + t2, after normalization.

      hi = t1 + t2
      lo = t2 - (hi - t1)
  
      return
      end

! ======================================================================
! 加算
!
!
      subroutine ddadd_func(org_hi, org_lo, hi, lo, kcount)

      implicit none

      double precision org_hi(kcount), org_lo(kcount)
      double precision hi, lo
      integer kcount
      integer i

      hi = org_hi(1)
      lo = org_lo(1)
      do i = 2, kcount
          call ddadd(hi, lo, org_hi(i), org_lo(i), hi, lo)
      end do

!      write (*, "('add func result',F24.16,E24.16)") hi, lo
      return
      end

! ======================================================================
! 乗算
!
! Author:
!      David H. Bailey
!      NERSC, Lawrence Berkeley Lab
!      （F77形式に変換して使用）
!
      subroutine ddmul(org1_hi, org1_lo, org2_hi, org2_lo, hi, lo)

      implicit none
      double precision org1_hi, org1_lo
      double precision org2_hi, org2_lo
      double precision hi, lo
      double precision a1, a2, b1, b2, cona, conb, c11, c21, c2
      double precision e, split, t1, t2
      parameter (split = 134217729.d0)

!>
!   On systems with a fused multiply add, such as IBM systems, it is faster to
!   uncomment the next two lines and comment out the following lines until !>.
!   On other systems, do the opposite.

!      c11 = org1_hi * org2_hi
!      c21 = org1_hi * org2_hi - c11

!   This splits org1_hi and org2_hi into high-order and low-order words.

      cona = org1_hi * split
      conb = org2_hi * split
      a1 = cona - (cona - org1_hi)
      b1 = conb - (conb - org2_hi)
      a2 = org1_hi - a1
      b2 = org2_hi - b1

!   Multilply org1_hi * org2_hi using Dekker's method.

      c11 = org1_hi * org2_hi
      c21 = (((a1 * b1 - c11) + a1 * b2) + a2 * b1) + a2 * b2
!>
!   Compute org1_hi * org2_lo + org1_lo * org2_hi (only high-order word is needed).

      c2 = org1_hi * org2_lo + org1_lo * org2_hi

!   Compute (c11, c21) + c2 using Knuth's trick, also adding low-order product.

      t1 = c11 + c2
      e = t1 - c11
      t2 = ((c2 - e) + (c11 - (t1 - e))) + c21 + org1_lo * org2_lo

!   The result is t1 + t2, after normalization.

      hi = t1 + t2
      lo = t2 - (hi - t1)

      return
      end

! ======================================================================
! 乗算
!
!
      subroutine ddmul_func(org_hi, org_lo, hi, lo, kcount)

      implicit none

      double precision org_hi(kcount), org_lo(kcount)
      double precision hi, lo
      integer kcount
      integer i

      hi = org_hi(1)
      lo = org_lo(1)
      do i = 2, kcount
          call ddmul(hi, lo, org_hi(i), org_lo(i), hi, lo)
!          write (*, "('mul func ',I3,F24.16,E24.16)") i, hi, lo
      end do

!      write (*, "('mul func result',F24.16,E24.16)") hi, lo
      return
      end

! ======================================================================
! 最大値
!
!
      subroutine ddmax_func(org_hi, org_lo, hi, lo, kcount)

      implicit none
      double precision org_hi(kcount), org_lo(kcount)
      double precision hi, lo
      integer kcount
      integer i

      hi = org_hi(1)
      lo = org_lo(1)
      do i = 2, kcount
          if (org_hi(i) > hi) then
              hi = org_hi(i)
              lo = org_lo(i)
          else if (org_hi(i) == hi) then
              if (org_lo(i) > lo) then
                  hi = org_hi(i)
                  lo = org_lo(i)
              end if
          end if
      end do

!      write (*, "('max func result',F24.16,E24.16)") hi, lo
      return
      end

! ======================================================================
! 最小値
!
!
      subroutine ddmin_func(org_hi, org_lo, hi, lo, kcount)

      implicit none
      double precision org_hi(kcount), org_lo(kcount)
      double precision hi, lo
      integer kcount
      integer i

      hi = org_hi(1)
      lo = org_lo(1)
      do i = 2, kcount
          if (org_hi(i) < hi) then
              hi = org_hi(i)
              lo = org_lo(i)
          else if (org_hi(i) == hi) then
              if (org_lo(i) < lo) then
                  hi = org_hi(i)
                  lo = org_lo(i)
              end if
          end if
      end do

!      write (*, "('min func result',F24.16,E24.16)") hi, lo
      return
      end


! ======================================================================
! 
!
!
      subroutine reduce_func(org_hi, org_lo, hi, lo, kcount, op)

      implicit none
      include 'mpif.h'

      double precision org_hi(kcount), org_lo(kcount)
      double precision hi, lo
      integer kcount, op

      if (op == MPI_PROD) then
          call ddmul_func(org_hi, org_lo, hi, lo, kcount)
      else if (op == MPI_SUM) then
          call ddadd_func(org_hi, org_lo, hi, lo, kcount)
      else if (op == MPI_MAX) then
          call ddmax_func(org_hi, org_lo, hi, lo, kcount)
      else if (op == MPI_MIN) then
          call ddmin_func(org_hi, org_lo, hi, lo, kcount)
      end if

      return
      end

