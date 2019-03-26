C
C     Copyright (C) 1998-1999 Japan Atomic Energy Research Institute.
C
c ***************************************************************
      SUBROUTINE GENMXRAN2( NMAX, N, NZR, NZERO,
     &                      VAL, ROWIDX, COLIDX, X, RHS,
     &                      NZROW, INDRC, A, IA, istart, iend, 
     &                      matrix_is_symmetric, seed )
c ***************************************************************
C
C     This routine makes randam sparce matrix(VAL), 
c 
C     NMAX   : Maximam Array size  N*NZR
C     N      : Number of variable
C     NZR    : Maximam size of nonzero element for any rows
c --
C     NZERO  : Size of array with coodinate formality. : Output
C     VAL    : VALue of array with coodinate formality   only DU side.
C     ROWIDX : Row index with coordinate formality       only DU side.
C     COLIDX : Column index with coordinate formality    only DU side.
C     RHS    : Vector of RHS.
C     X      : Analisys Resultant
c
C     NZROW  : stack array for nonzero number of each rows
c 
C     A      : VALue of array of full matrix(work).
C     IA     : Nonzero columns of full matrix(work).
C     NPES   : Number of PE. 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
      DIMENSION    VAL(*),RHS(*),X(*)
      INTEGER      ROWIDX(*),COLIDX(*)
      DIMENSION    A(N,NZR)
      
      INTEGER      IA(N,NZR),NZROW(N)
      integer INDRC(NZR)
      integer istart, iend
      integer seed
      logical matrix_is_symmetric
c     
      DATA ZERO,ONE/0.0D+00,1.0D+00/
      DATA FACWID/2.0D-01/
C
      DO I=1,N
         NZROW(I)=0
      ENDDO

      NBWID=N*FACWID
      NBDEL=0
      
      DO I=1,N-1
         NZ=NZR-NZROW(I)-1
         NZ=MAX(0,NZ)
         NZ=MIN(N-I,NZ)

         DO J=1,NZR
            INDRC(J)=0
         ENDDO   
         IF(I.GT.N-NBWID) NBDEL=NBDEL+1
         
         DO J=1,NZ
 100        K=(NBWID-NBDEL)*RANDOM2( seed )+0.5
            IC=I+1+K
            IC=MIN(IC,N)
            DO M=1,J-1
               IF(INDRC(M).EQ.IC) GOTO 100
            ENDDO
            INDRC(J)=IC
         ENDDO

         DO J=1,NZ
            IC=INDRC(J)
            IF(NZROW(I).GE.(NZR-1) .OR. NZROW(IC).GE.(NZR-1) ) GOTO 200
            NZROW(I) =NZROW(I) +1
            NZROW(IC)=NZROW(IC)+1
            IA(I,NZROW(I))=IC
            IA(IC,NZROW(IC))=I
            A(I,NZROW(I))=-RANDOM2( seed )
            A(IC,NZROW(IC))=A(I,NZROW(I))
C          IF(NZROW(I).GT.(NZR-1) .OR. NZROW(IC).GT.(NZR-1) ) GOTO 200
         ENDDO
 200     CONTINUE
      ENDDO    

c
C... Add Diagonal and set Matrix element
C
      DO I=1,N
         SUM=ZERO
         NZ=NZROW(I)
         DO J=1,NZ
            SUM=SUM+ABS(A(I,J))
         ENDDO
         IA(I,NZ+1)=I
         A(I,NZ+1)=SUM+ONE
c        A(I,NZ+1)=SUM+0.3
         NZROW(I)=NZROW(I)+1
      ENDDO
C     
C...  Sort Columns index
C
      DO I=1,N
         NZ=NZROW(I)
         DO J=1,NZ-1
            DO K=J+1,NZ
               IF(IA(I,J).GT.IA(I,K))THEN
                  ITMP=IA(I,J)
                  IA(I,J)=IA(I,K)
                  IA(I,K)=ITMP
                  ATMP=A(I,J)
                  A(I,J)=A(I,K)
                  A(I,K)=ATMP
               ENDIF
            ENDDO
         ENDDO
      ENDDO        
C
C... Cal RHS
C
      DO I=1,N
c.3.14        X(I)=ONE
         Y=FLOAT(I)
         X(I)=SQRT(Y)
      ENDDO
      DO I=1,N
         SUM=ZERO
         NZ=NZROW(I)
         DO J=1,NZ
            SUM=SUM+A(I,J)*X(IA(I,J))
         ENDDO
         RHS(I)=SUM
      ENDDO
c      do i=1,n
c        write(*,*)'RHS(',i,')=',rhs(i)
c      enddo
c      stop
C
C.... Remake Coodinate form.
C ---------------------------------------
c      IRF=IROW(MY)
c      IRL=IROW(MY+1)-1

      IRF = istart
      IRL = iend
C ---------------------------------------
      NZERO=0

      DO IR=IRF,IRL
         NZ=NZROW(IR)
         DO J=1,NZ
            IC=IA(IR,J)

            if ( matrix_is_symmetric ) then
               IF(IC.GE.IR)THEN
                  NZERO=NZERO+1
                  ROWIDX(NZERO)=IR
                  COLIDX(NZERO)=IC
                  VAL(NZERO)=A(IR,J)             
               ENDIF

            else
               NZERO=NZERO+1
               ROWIDX(NZERO)=IR
               COLIDX(NZERO)=IC
               VAL(NZERO)=A(IR,J)             
            endif

         ENDDO
      ENDDO
C
      RETURN
      END

c ***********************************************************
      SUBROUTINE GENMXRAN(NMAX,N,NZR,NZERO,
     1                 VAL,ROWIDX,COLIDX,X,RHS,
     1                 NZROW,INDRC,A,IA,IROW,NPES,MY)
c ***********************************************************
C
C     This routine makes randam sparce matrix(VAL), 
c 
C     NMAX   : Maximam Array size  N*NZR
C     N      : Number of variable
C     NZR    : Maximam size of nonzero element for any rows
c --
C     NZERO  : Size of array with coodinate formality. : Output
C     VAL    : VALue of array with coodinate formality   only DU side.
C     ROWIDX : Row index with coordinate formality       only DU side.
C     COLIDX : Column index with coordinate formality    only DU side.
C     RHS    : Vector of RHS.
C     X      : Analisys Resultant
c
C     NZROW  : stack array for nonzero number of each rows
c 
C     A      : VALue of array of full matrix(work).
C     IA     : Nonzero columns of full matrix(work).
C     NPES   : Number of PE. 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION    VAL(*),RHS(*),X(*)
      INTEGER      ROWIDX(*),COLIDX(*)
      DIMENSION    A(N,NZR)
      INTEGER      IA(N,NZR),NZROW(N)
      INTEGER      IROW(0:NPES),INDRC(NZR)
      DATA ZERO,ONE/0.0D+00,1.0D+00/
      DATA FACWID/2.0D-01/
C
      DO I=1,N
        NZROW(I)=0
      ENDDO
      NBWID=N*FACWID
      NBDEL=0
      DO I=1,N-1
        NZ=NZR-NZROW(I)-1
        NZ=MAX(0,NZ)
        NZ=MIN(N-I,NZ)
        DO J=1,NZR
          INDRC(J)=0
        ENDDO   
        IF(I.GT.N-NBWID) NBDEL=NBDEL+1
        DO J=1,NZ
100       K=(NBWID-NBDEL)*RANDOM()+0.5
          IC=I+1+K
          IC=MIN(IC,N)
          DO M=1,J-1
            IF(INDRC(M).EQ.IC) GOTO 100
          ENDDO
          INDRC(J)=IC
        ENDDO
        DO J=1,NZ
          IC=INDRC(J)
          IF(NZROW(I).GE.(NZR-1) .OR. NZROW(IC).GE.(NZR-1) ) GOTO 200
          NZROW(I) =NZROW(I) +1
          NZROW(IC)=NZROW(IC)+1
          IA(I,NZROW(I))=IC
          IA(IC,NZROW(IC))=I
          A(I,NZROW(I))=-RANDOM()
          A(IC,NZROW(IC))=A(I,NZROW(I))
C          IF(NZROW(I).GT.(NZR-1) .OR. NZROW(IC).GT.(NZR-1) ) GOTO 200
        ENDDO
200     CONTINUE
      ENDDO    
C
C... Add Diagonal and set Matrix element
C
      DO I=1,N
        SUM=ZERO
        NZ=NZROW(I)
        DO J=1,NZ
          SUM=SUM+ABS(A(I,J))
        ENDDO
        IA(I,NZ+1)=I
        A(I,NZ+1)=SUM+ONE
c        A(I,NZ+1)=SUM+0.3
        NZROW(I)=NZROW(I)+1
      ENDDO
C
C... Sort Columns index
C
      DO I=1,N
        NZ=NZROW(I)
        DO J=1,NZ-1
          DO K=J+1,NZ
            IF(IA(I,J).GT.IA(I,K))THEN
              ITMP=IA(I,J)
              IA(I,J)=IA(I,K)
              IA(I,K)=ITMP
              ATMP=A(I,J)
              A(I,J)=A(I,K)
              A(I,K)=ATMP
            ENDIF
          ENDDO
        ENDDO
      ENDDO        
C
C... Cal RHS
C
      DO I=1,N
c.3.14        X(I)=ONE
        Y=FLOAT(I)
        X(I)=SQRT(Y)
      ENDDO
      DO I=1,N
        SUM=ZERO
        NZ=NZROW(I)
        DO J=1,NZ
          SUM=SUM+A(I,J)*X(IA(I,J))
        ENDDO
        RHS(I)=SUM
      ENDDO
c      do i=1,n
c        write(*,*)'RHS(',i,')=',rhs(i)
c      enddo
c      stop
C
C.... Remake Coodinate form.
C
      IRF=IROW(MY)
      IRL=IROW(MY+1)-1
C
      NZERO=0
      DO IR=IRF,IRL
        NZ=NZROW(IR)
        DO J=1,NZ
          IC=IA(IR,J)
          IF(IC.GE.IR)THEN
            NZERO=NZERO+1
            ROWIDX(NZERO)=IR
            COLIDX(NZERO)=IC
            VAL(NZERO)=A(IR,J)             
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END

c ***********************************************************
      function random()
c ***********************************************************
      implicit double precision ( a-h, o-z )
      data dseed/12357.0d0/
      data d2p31m / 2147483647.0d0 /
      data d2p31  / 2147483648.0d0 /
      dseed1 = dseed
      dseed  = dmod(16807.0d0*dseed,d2p31m)
      random =dseed/d2p31

      return
      end

c ***********************************************************
      real*8 function random2( seed )
c ***********************************************************
      implicit none
c inout
      integer seed
c const
      real*8 d2p31, d2p31m
      data d2p31m / 2147483647.0d0 /
      data d2p31  / 2147483648.0d0 /
c local
      real*8 dseed1, dseed
      save dseed1, dseed
      logical FirstFlag
      save FirstFlag
      data FirstFlag / .true. /
c init
      if ( FirstFlag ) then
         dseed = seed
         FirstFlag = .false.
      endif
c begin
      dseed1 = dseed
      dseed  = dmod( 16807.0d0 *dseed, d2p31m )
      random2 = dseed /d2p31
c
      return
      end






      
