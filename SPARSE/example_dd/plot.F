C
C     Copyright (C) 1998-1999 Japan Atomic Energy Research Institute.
C
      subroutine plot(zerr, n)

      implicit none

      integer  nfplot, n, i
      parameter (nfplot = 8)
      double precision  zerr(n)


      write(nfplot,'(a)') '# step   zerr'
      do i = 1, n
        write(nfplot,8001) i, zerr(i)
      end do

 8000 format(2x,i8,1pe15.6)
 8001 format(2x,i8,1pe18.9)

      return
      end
