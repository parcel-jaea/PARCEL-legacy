C
C     Copyright (C) 1998-1999 Japan Atomic Energy Research Institute.
C
c -------------------
c
c   COO -> CRS format
c
c -------------------
      subroutine JAERI_CO2CRTB(ndim, nnz, istart, n,
     &                   row_co, col_co,
     &                   col_cr, ptr_cr,
     &                   co2cr_tbl)
C --- in
C     ndim: dimension of global matrix
C     nnz: number of nonzero element
C     istart, n:
C     row_co: row index of original matrix
C     col_co: column index of original matrix
C
C ---out
C     col_cr: column index of CR matrix
C     ptr_cr: ptr index of CR matrix
C     co2cr_tbl: table for C -> CR transform

      implicit none
      include 'mpif.h'
      include 'JAERI_iles.h'

c --- arguments
      integer  ndim, nnz,  istart, n
      integer  row_co(nnz), col_co(nnz), col_cr(nnz), ptr_cr(n+1)
      integer  co2cr_tbl(nnz)

c --- local
      integer  i, ir, ip, ierr
#ifdef SORT_IN_ROW
      integer  i1, n1
#endif

c     --- count nonzero elements in each row
*note: use col_cr as work array

      do ir = 1, n
        col_cr(ir) = 0.0d0
      end do
      do i = 1, nnz
        ir = row_co(i) - istart + 1
        if (ir .lt. 1  .or.  ir .gt. n) then
          write(*,*) '***co2cr_tbl: error'
          call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        col_cr(ir) = col_cr(ir) + 1
      end do
      ptr_cr(1) = 1
      do ir = 1, n
        ptr_cr(ir+1) = ptr_cr(ir) + col_cr(ir)
      end do

c     --- set translation table
*note: use col_cr as work array
      do ir = 1, n
        col_cr(ir) = ptr_cr(ir)
      end do
      do i = 1, nnz
        ir = row_co(i) - istart + 1
        ip = col_cr(ir)
        co2cr_tbl(ip) = i
        col_cr(ir) = col_cr(ir) + 1
      end do

c     --- sotre column numbers
      do i = 1, nnz
        col_cr(i) = col_co(co2cr_tbl(i))
      end do
      
#ifdef SORT_IN_ROW
c     --- sort in row
      do ir = 1, n
        i1 = ptr_cr(ir)
        n1 = ptr_cr(ir+1) - i1
        call hsort2d(n1, col_cr(i1), co2cr_tbl(i1))
      end do
#endif
      end

      subroutine JAERI_CO2CR(nnz, val_co, val_cr, co2cr_tbl)
      implicit none

c --- arguments
      integer  nnz, co2cr_tbl(nnz)
      double precision  val_co(nnz), val_cr(nnz)

c --- local
      integer  i

      do i = 1, nnz
        val_cr(i) = val_co(co2cr_tbl(i))
      end do

      return
      end

c -------------------
c
c   COO -> DIA format
c
c -------------------
      subroutine co2dia_count_nnd( ndim, num_nonzero, 
     &           matA_rowindex_coor, 
     &           matA_colindex_coor, istart, iend, nnd )
      implicit none
c input
      integer ndim
      integer num_nonzero
      integer matA_rowindex_coor( num_nonzero )
      integer matA_colindex_coor( num_nonzero )
      integer istart, iend
c output
      integer nnd
c local
      integer i, icol, irow, itmp
      integer, allocatable :: Flag(:)
c begin
      integer imn, imx
      imn = 0
      imx = 0
      Do i=1, num_nonzero
         irow = matA_rowindex_coor(i)
         icol = matA_colindex_coor(i)
         if ( irow >= istart .and. irow <= iend ) then
            itmp = icol -irow
            if (imn > itmp) imn = itmp
            if (imx < itmp) imx = itmp
         endif
      End do
      allocate( Flag(imn:imx) ); Flag = 0
      
      Do i=1, num_nonzero
         irow = matA_rowindex_coor(i)
         icol = matA_colindex_coor(i)
         if ( irow >= istart .and. irow <= iend ) then
            itmp = icol -irow
            Flag(itmp) = 1
         endif
      End do
c
      nnd = 0
      Do i=imn, imx
         if ( Flag(i) == 1 ) then
            nnd = nnd +1
         endif
      End do
      deallocate( Flag )
      return
      end

      subroutine co2dia( ndim, num_nonzero, 
     &                   matA_rowindex_coor, matA_colindex_coor, 
     &                   matA_coor, istart, iend, 
     &                   ndist, nnd, matA_ioff_dia, matA_dia )
      implicit none
c input
      integer ndim
      integer num_nonzero
      integer nnd
      integer ndist
      integer istart, iend
      integer matA_rowindex_coor( num_nonzero )
      integer matA_colindex_coor( num_nonzero )
      real*8  matA_coor( num_nonzero )
c output
      integer matA_ioff_dia( nnd )
      real*8 matA_dia( ndist*nnd )
c local
      integer i, j, icol, irow, itmp, i1, j1
      integer, allocatable :: Flag(:)
c begin
      integer imn, imx
      imn = 0
      imx = 0
      Do i=1, num_nonzero
         irow = matA_rowindex_coor(i)
         icol = matA_colindex_coor(i)
         if ( irow >= istart .and. irow <= iend ) then
            itmp = icol -irow
            if (imn > itmp) imn = itmp
            if (imx < itmp) imx = itmp
         endif
      End do
      allocate( Flag(imn:imx) ); Flag = 0
      
      Do i=1, num_nonzero
         irow = matA_rowindex_coor(i)
         icol = matA_colindex_coor(i)
         if ( irow >= istart .and. irow <= iend ) then
            itmp = icol -irow
            Flag(itmp) = 1
         endif
      End do
c
      nnd = 0
      Do i=imn, imx
         if ( Flag(i) == 1 ) then
            nnd = nnd +1
            matA_ioff_dia(nnd) = i
         endif
      End do
!
      Do i=1, num_nonzero
         irow = matA_rowindex_coor(i)
         icol = matA_colindex_coor(i)

         if ( irow >= istart .and. irow <= iend ) then
            itmp = icol -irow
            Do j=1, nnd
               j1 = matA_ioff_dia(j)
               if ( j1 == itmp ) exit
            End do
            i1 = irow -istart +1
            matA_dia( ndist*(j-1) +i1 ) = matA_coor(i)
         endif
      End do
c
      deallocate( Flag )

      return
      end


