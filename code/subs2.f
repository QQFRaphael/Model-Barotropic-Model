      SUBROUTINE CGEFA(A,LDA,N,IPVT,INFO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LDA,N,IPVT(1),INFO
      COMPLEX*16 A(LDA,1)
C
C     CGEFA FACTORS A COMPLEX*16 MATRIX BY GAUSSIAN ELIMINATION.
C
C     CGEFA IS USUALLY CALLED BY CGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR CGECO) = (1 + 9/N)*(TIME FOR CGEFA) .
C
C     ON ENTRY
C
C        A       COMPLEX*16(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT CGESL OR CGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN CGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CSCAL,ICAMAX
C     FORTRAN ABS,AIMAG,REAL
C
C     INTERNAL VARIABLES
C
      COMPLEX*16 T
      INTEGER ICAMAX,J,K,KP1,L,NM1
C
c     COMPLEX*16 ZDUM
c     REAL*8 CABS1
c     CABS1(ZDUM) = DABS(DBLE(ZDUM)) + DABS(DIMAG(ZDUM))
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = ICAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         CABSB = DABS(DBLE(A(L,K))) + DABS(DIMAG(A(L,K)))
         IF (CABSB .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -(1.0D0,0.0D0)/A(K,K)
            CALL CSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL CAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      CABSB = DABS(DBLE(A(N,N))) + DABS(DIMAG(A(N,N)))
      IF (CABSB .EQ. 0.0D0) INFO = N
      RETURN
      END
      SUBROUTINE CGESL(A,LDA,N,IPVT,B,JOB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LDA,N,IPVT(1),JOB
      COMPLEX*16 A(LDA,1),B(1)
C
C     CGESL SOLVES THE COMPLEX*16 SYSTEM
C     A * X = B  OR  CTRANS(A) * X = B
C     USING THE FACTORS COMPUTED BY CGECO OR CGEFA.
C
C     ON ENTRY
C
C        A       COMPLEX*16(LDA, N)
C                THE OUTPUT FROM CGECO OR CGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM CGECO OR CGEFA.
C
C        B       COMPLEX*16(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  CTRANS(A)*X = B  WHERE
C                            CTRANS(A)  IS THE CONJUGATE TRANSPOSE.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF CGECO HAS SET RCOND .GT. 0.0
C        OR CGEFA HAS SET INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL CGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL CGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CDOTC
C     FORTRAN CONJG
C
C     INTERNAL VARIABLES
C
      COMPLEX*16 CDOTC,T
      INTEGER K,KB,L,NM1
C
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL CAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL CAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  CTRANS(A) * X = B
C        FIRST SOLVE  CTRANS(U)*Y = B
C
         DO 60 K = 1, N
            T = CDOTC(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/DCONJG(A(K,K))
   60    CONTINUE
C
C        NOW SOLVE CTRANS(L)*X = Y
C
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + CDOTC(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
      SUBROUTINE CAXPY(N,CA,CX,INCX,CY,INCY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      COMPLEX*16 CX(1),CY(1),CA
      INTEGER I,INCX,INCY,IX,IY,N
C
      IF(N.LE.0)RETURN
      IF (DABS(DBLE(CA)) + DABS(DIMAG(CA)) .EQ. 0.0D0 ) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        CY(IY) = CY(IY) + CA*CX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 DO 30 I = 1,N
        CY(I) = CY(I) + CA*CX(I)
   30 CONTINUE
      RETURN
      END
      COMPLEX*16 FUNCTION CDOTC(N,CX,INCX,CY,INCY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS, CONJUGATING THE FIRST
C     VECTOR.
C     JACK DONGARRA, LINPACK,  3/11/78.
C
      COMPLEX*16 CX(1),CY(1),CTEMP
      INTEGER I,INCX,INCY,IX,IY,N
C
      CTEMP = DCMPLX(0.0D0,0.0D0)
      CDOTC = DCMPLX(0.0D0,0.0D0)
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        CTEMP = CTEMP + DCONJG(CX(IX))*CY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      CDOTC = CTEMP
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 DO 30 I = 1,N
        CTEMP = CTEMP + DCONJG(CX(I))*CY(I)
   30 CONTINUE
      CDOTC = CTEMP
      RETURN
      END
      SUBROUTINE  CSCAL(N,CA,CX,INCX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     SCALES A VECTOR BY A CONSTANT.
C     JACK DONGARRA, LINPACK,  3/11/78.
C
      COMPLEX*16 CA,CX(1)
      INTEGER I,INCX,N,NINCX
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        CX(I) = CA*CX(I)
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DO 30 I = 1,N
        CX(I) = CA*CX(I)
   30 CONTINUE
      RETURN
      END
      INTEGER FUNCTION ICAMAX(N,CX,INCX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      COMPLEX*16 CX(1)
      REAL*8 SMAX
      INTEGER I,INCX,IX,N
c     COMPLEX*16 ZDUM
c     REAL*8 CABS1
c     CABS1(ZDUM) = DABS(DBLE(ZDUM)) + DABS(DIMAG(ZDUM))
C
      ICAMAX = 0
      IF( N .LT. 1 ) RETURN
      ICAMAX = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      SMAX = DABS(DBLE(CX(1)))+DABS(DIMAG(CX(1))) 
      IX = IX + INCX
      DO 10 I = 2,N
         CABSB = DABS(DBLE(CX(IX))) + DABS(DIMAG(CX(IX)))
         IF(CABSB.LE.SMAX) GO TO 5
         ICAMAX = I
         SMAX = CABSB
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 SMAX = DABS(DBLE(CX(1)))+DABS(DIMAG(CX(1))) 
      DO 30 I = 2,N
         CABSB = DABS(DBLE(CX(I))) + DABS(DIMAG(CX(I)))
         IF(CABSB.LE.SMAX) GO TO 30
         ICAMAX = I
         SMAX = CABSB
   30 CONTINUE
      RETURN
      END
      SUBROUTINE MTXADD (X,M,KK,Y,II,N,Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(M,KK),Y(II,N),Z(M,N)
      DO 10 I=1,M
      DO 10 J=1,N
   10 Z(I,J)=X(I,J)+Y(I,J)
      RETURN
      END
      SUBROUTINE MTXMUL(X,M,KK,Y,II,N,Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(M,KK),Y(II,N),Z(M,N)
      DO 20 I=1,M
      DO 20 J=1,N
      Z(I,J)=0.0D0
      DO 20 K=1,KK
   20 Z(I,J)=Z(I,J)+X(I,K)*Y(K,J)
      RETURN
      END
      SUBROUTINE MTXVMU (X,Y,Z,KP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION Y(KP,KP)
      COMPLEX*16 X(KP),Z(KP)
      DO 1 I=1,KP
      X(I)=DCMPLX(0.0D0,0.0D0)
      DO 1 J=1,KP
    1 X(I)=X(I)+DCMPLX(Y(I,J))*Z(J)
      RETURN
      END

      SUBROUTINE MINV (A,N,D,L,M)
C***********************************************************************
C
C     IMINV : INVERTS A MATRIX IN PLACE USING GAUSS-JORDAN TECHNIQUE.
C             IN ADDITION, THE DETERMINANT IS CALCULATED.
C             INPUT MATRIX IS DESTROYED.
C
C***********************************************************************
C
C     IMINV IS CALLED BY THE SUBROUTINES AMHMTM AND SEMICD.
C
C     IMINV CALLS NO SUBROUTINES.
C
C***********************************************************************
C
C    ARGUMENT(DIMENSIONS)                       DESCRIPTION
C
C           A(N,N)                INPUT : MATRIX TO BE INVERTED.
C                                         COMPUTED BY CALLING ROUTINE.
C                                OUTPUT : INVERSE OF INPUT MATRIX "A".
C             N                   INPUT : ORDER OF MATRIX "A". SET BY
C                                         CALLING ROUTINE.
C             D                  OUTPUT : DETERMINANT OF INPUT MATRIX
C                                         "A".
C             L                  OUTPUT : WORK VECTOR OF DIMENSION "N".
C             M                  OUTPUT : WORK VECTOR OF DIMENSION "N".
C
C***********************************************************************
C
C     IMINV DOES NOT REFER TO ANY COMMONS.
C
C***********************************************************************
C
C... CREATED AT NMC JUNE 15 1984.......
C
C     ..................................................................
C
CC--------------------------------------------------------
      DIMENSION A(*),L(*),M(*)
C
C        ...............................................................
C
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION
C        STATEMENT WHICH FOLLOWS.
C
      DOUBLE PRECISION A, D, BIGA, HOLD
C
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS
C        ROUTINE.
C
C        THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO
C        CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  ABS IN STATEMEN
C        10 MUST BE CHANGED TO DABS  .
C
C        ...............................................................
C
C        SEARCH FOR LARGEST ELEMENT
C
CVD$R NOVECTOR
      D=1.0D0
      NK=-N
      DO 80 K=1,N
      NK=NK+N
      L(K)=K
      M(K)=K
      KK=NK+K
      BIGA=A(KK)
      DO 20 J=K,N
      IZ=N*(J-1)
      DO 20 I=K,N
      IJ=IZ+I
   10 IF (DABS(BIGA)-DABS(A(IJ))) 15,20,20
   15 BIGA=A(IJ)
      L(K)=I
      M(K)=J
   20 CONTINUE
C
C        INTERCHANGE ROWS
C
      J=L(K)
      IF(J-K) 35,35,25
   25 KI=K-N
      DO 30 I=1,N
      KI=KI+N
      HOLD=-A(KI)
      JI=KI-K+J
      A(KI)=A(JI)
   30 A(JI) =HOLD
C
C        INTERCHANGE COLUMNS
C
   35 I=M(K)
      IF(I-K) 45,45,38
   38 JP=N*(I-1)
      DO 40 J=1,N
      JK=NK+J
      JI=JP+J
      HOLD=-A(JK)
      A(JK)=A(JI)
   40 A(JI) =HOLD
C
C        DIVIDE COLUMN BY MINUS PIVOT (VALUE OF PIVOT ELEMENT IS
C        CONTAINED IN BIGA)
C
   45 IF(BIGA) 48,46,48
   46 D=1.0D0
      RETURN
   48 DO 55 I=1,N
      IF(I-K) 50,55,50
   50 IK=NK+I
      A(IK)=A(IK)/(-BIGA)
   55 CONTINUE
C
C        REDUCE MATRIX
C
      DO 65 I=1,N
      IK=NK+I
      IJ=I-N
      DO 65 J=1,N
      IJ=IJ+N
      IF(I-K) 60,65,60
   60 IF(J-K) 62,65,62
   62 KJ=IJ-I+K
      A(IJ)=A(IK)*A(KJ)+A(IJ)
   65 CONTINUE
C
C        DIVIDE ROW BY PIVOT
C
      KJ=K-N
      DO 75 J=1,N
      KJ=KJ+N
      IF(J-K) 70,75,70
   70 A(KJ)=A(KJ)/BIGA
   75 CONTINUE
C
C        PRODUCT OF PIVOTS
C
      D=D*BIGA
C
C        REPLACE PIVOT BY RECIPROCAL
C
      A(KK)=1.0D0/BIGA
   80 CONTINUE
C
C        FINAL ROW AND COLUMN INTERCHANGE
C
      K=N
  100 K=(K-1)
      IF(K) 150,150,105
  105 I=L(K)
      IF(I-K) 120,120,108
  108 JQ=N*(K-1)
      JR=N*(I-1)
      DO 110 J=1,N
      JK=JQ+J
      HOLD=A(JK)
      JI=JR+J
      A(JK)=-A(JI)
  110 A(JI) =HOLD
  120 J=M(K)
      IF(J-K) 100,100,125
  125 KI=K-N
      DO 130 I=1,N
      KI=KI+N
      HOLD=A(KI)
      JI=KI-K+J
      A(KI)=-A(JI)
  130 A(JI) =HOLD
      GO TO 100
  150 RETURN
      END

      SUBROUTINE FAX(IFAX,N,MODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IFAX(10)
      NN=N
      IF (IABS(MODE).EQ.1) GO TO 10
      IF (IABS(MODE).EQ.8) GO TO 10
      NN=N/2
      IF ((NN+NN).EQ.N) GO TO 10
      IFAX(1)=-99
      RETURN
   10 K=1
C     TEST FOR FACTORS OF 4
   20 IF (MOD(NN,4).NE.0) GO TO 30
      K=K+1
      IFAX(K)=4
      NN=NN/4
      IF (NN.EQ.1) GO TO 80
      GO TO 20
C     TEST FOR EXTRA FACTOR OF 2
   30 IF (MOD(NN,2).NE.0) GO TO 40
      K=K+1
      IFAX(K)=2
      NN=NN/2
      IF (NN.EQ.1) GO TO 80
C     TEST FOR FACTORS OF 3
   40 IF (MOD(NN,3).NE.0) GO TO 50
      K=K+1
      IFAX(K)=3
      NN=NN/3
      IF (NN.EQ.1) GO TO 80
      GO TO 40
C     NOW FIND REMAINING FACTORS
   50 L=5
      INC=2
C     INC ALTERNATELY TAKES ON VALUES 2 AND 4
   60 IF (MOD(NN,L).NE.0) GO TO 70
      K=K+1
      IFAX(K)=L
      NN=NN/L
      IF (NN.EQ.1) GO TO 80
      GO TO 60
   70 L=L+INC
      INC=6-INC
      GO TO 60
   80 IFAX(1)=K-1
C     IFAX(1) CONTAINS NUMBER OF FACTORS
C     IFAX(1) CONTAINS NUMBER OF FACTORS
      NFAX=IFAX(1)
C     SORT FACTORS INTO ASCENDING ORDER
      IF (NFAX.EQ.1) GO TO 110
      DO 100 II=2,NFAX
      ISTOP=NFAX+2-II
      DO 90 I=2,ISTOP
      IF (IFAX(I+1).GE.IFAX(I)) GO TO 90
      ITEM=IFAX(I)
      IFAX(I)=IFAX(I+1)
      IFAX(I+1)=ITEM
   90 CONTINUE
  100 CONTINUE
  110 CONTINUE
      RETURN
      END
      SUBROUTINE FFTRIG(TRIGS,N,MODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TRIGS(1)
      PI=2.0D0*DASIN(1.0D0)
      IMODE=IABS(MODE)
      NN=N
      IF (IMODE.GT.1.AND.IMODE.LT.6) NN=N/2
      DEL=(PI+PI)/DBLE(NN)
      L=NN+NN
      DO 10 I=1,L,2
      ANGLE=0.5D0*DBLE(I-1)*DEL
      TRIGS(I)=DCOS(ANGLE)
      TRIGS(I+1)=DSIN(ANGLE)
   10 CONTINUE
      IF (IMODE.EQ.1) RETURN
      IF (IMODE.EQ.8) RETURN
      DEL=0.5D0*DEL
      NH=(NN+1)/2
      L=NH+NH
      LA=NN+NN
      DO 20 I=1,L,2
      ANGLE=0.5D0*DBLE(I-1)*DEL
      TRIGS(LA+I)=DCOS(ANGLE)
      TRIGS(LA+I+1)=DSIN(ANGLE)
   20 CONTINUE
      IF (IMODE.LE.3) RETURN
      DEL=0.5D0*DEL
      LA=LA+NN
      IF (MODE.EQ.5) GO TO 40
      DO 30 I=2,NN
      ANGLE=DBLE(I-1)*DEL
      TRIGS(LA+I)=2.0D0*DSIN(ANGLE)
   30 CONTINUE
      RETURN
   40 CONTINUE
      DEL=0.5D0*DEL
      DO 50 I=2,N
      ANGLE=DBLE(I-1)*DEL
      TRIGS(LA+I)=DSIN(ANGLE)
   50 CONTINUE
      RETURN
      END

      SUBROUTINE VPASSM(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,IFAC,LA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION A(*),B(*),C(*),D(*),TRIGS(*)
      SAVE INT,RADI,SIN60,SIN36,SIN72,COS36,COS72
      DATA INT/0/
C
      IF(INT.EQ.0) THEN
        RADI=DATAN(1.0D0)/45.D0
        SIN60=DSIN(60.D0*RADI)
        SIN36=DSIN(36.D0*RADI)
        SIN72=DSIN(72.D0*RADI)
        COS36=DCOS(36.D0*RADI)
        COS72=DCOS(72.D0*RADI)
        INT=1
      END IF
      M=N/IFAC
      IINK=M*INC1
      JINK=LA*INC2
      JUMP=(IFAC-1)*JINK
      IBASE=0
      JBASE=0
      IGO=IFAC-1
      IF(IGO.GT.4) RETURN
      GO TO (10,50,90,130),IGO
C
   10 CONTINUE
      IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      DO 20 L=1,LA
      I=IBASE
      J=JBASE
      DO 15 IJK=1,LOT
      C(JA+J)=A(IA+I)+A(IB+I)
      D(JA+J)=B(IA+I)+B(IB+I)
      C(JB+J)=A(IA+I)-A(IB+I)
      D(JB+J)=B(IA+I)-B(IB+I)
      I=I+INC3
      J=J+INC4
   15 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   20 CONTINUE
      IF(LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 40 K=LA1,M,LA
      KB=K+K-2
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      DO 30 L=1,LA
      I=IBASE
      J=JBASE
      DO 25 IJK=1,LOT
      C(JA+J)=A(IA+I)+A(IB+I)
      D(JA+J)=B(IA+I)+B(IB+I)
      C(JB+J)=C1*(A(IA+I)-A(IB+I))-S1*(B(IA+I)-B(IB+I))
      D(JB+J)=S1*(A(IA+I)-A(IB+I))+C1*(B(IA+I)-B(IB+I))
      I=I+INC3
      J=J+INC4
   25 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   30 CONTINUE
      JBASE=JBASE+JUMP
   40 CONTINUE
      RETURN
C
   50 CONTINUE
      IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      DO 60 L=1,LA
      I=IBASE
      J=JBASE
      DO 55 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
      C(JB+J)=(A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I)))
      C(JC+J)=(A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I)))
      D(JB+J)=(B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I)))
      D(JC+J)=(B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I)))
      I=I+INC3
      J=J+INC4
   55 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   60 CONTINUE
      IF(LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 80 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      DO 70 L=1,LA
      I=IBASE
      J=JBASE
      DO 65 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
      C(JB+J)=
     $    C1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I))))
     $   -S1*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
      D(JB+J)=
     $    S1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I))))
     $   +C1*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
      C(JC+J)=
     $    C2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I))))
     $   -S2*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
      D(JC+J)=
     $    S2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I))))
     $   +C2*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
      I=I+INC3
      J=J+INC4
   65 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   70 CONTINUE
      JBASE=JBASE+JUMP
   80 CONTINUE
      RETURN
C
   90 CONTINUE
      IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      DO 100 L=1,LA
      I=IBASE
      J=JBASE
      DO 95 IJK=1,LOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      C(JC+J)=(A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
      D(JC+J)=(B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I))
      C(JB+J)=(A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I))
      C(JD+J)=(A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I))
      D(JB+J)=(B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I))
      D(JD+J)=(B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I))
      I=I+INC3
      J=J+INC4
   95 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  100 CONTINUE
      IF(LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 120 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      KD=KC+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      DO 110 L=1,LA
      I=IBASE
      J=JBASE
      DO 105 IJK=1,LOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
      C(JC+J)=
     $    C2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
     $   -S2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
      D(JC+J)=
     $    S2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
     $   +C2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
      C(JB+J)=
     $    C1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))
     $   -S1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
      D(JB+J)=
     $    S1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))
     $   +C1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
      C(JD+J)=
     $    C3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))
     $   -S3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
      D(JD+J)=
     $    S3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))
     $   +C3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
      I=I+INC3
      J=J+INC4
  105 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  110 CONTINUE
      JBASE=JBASE+JUMP
  120 CONTINUE
      RETURN
C
  130 CONTINUE
      IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      IE=ID+IINK
      JE=JD+JINK
      DO 140 L=1,LA
      I=IBASE
      J=JBASE
      DO 135 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
      C(JB+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     $  -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
      C(JE+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     $  +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
      D(JB+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     $  +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
      D(JE+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     $  -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
      C(JC+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     $  -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
      C(JD+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     $  +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
      D(JC+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     $  +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
      D(JD+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     $  -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
      I=I+INC3
      J=J+INC4
  135 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  140 CONTINUE
      IF(LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 160 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      KD=KC+KB
      KE=KD+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      C4=TRIGS(KE+1)
      S4=TRIGS(KE+2)
      DO 150 L=1,LA
      I=IBASE
      J=JBASE
      DO 145 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
      C(JB+J)=
     $    C1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     $      -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     $   -S1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     $      +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      D(JB+J)=
     $    S1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     $      -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     $   +C1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     $      +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      C(JE+J)=
     $    C4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     $      +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     $   -S4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     $      -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      D(JE+J)=
     $    S4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     $      +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     $   +C4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     $      -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      C(JC+J)=
     $    C2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     $      -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     $   -S2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     $      +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      D(JC+J)=
     $    S2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     $      -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     $   +C2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     $      +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      C(JD+J)=
     $    C3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     $      +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     $   -S3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     $      -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      D(JD+J)=
     $    S3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     $      +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     $   +C3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     $      -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      I=I+INC3
      J=J+INC4
  145 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  150 CONTINUE
      JBASE=JBASE+JUMP
  160 CONTINUE
      RETURN
C
      END


      SUBROUTINE FFT99M(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A(*),WORK(*),TRIGS(*)
      DIMENSION IFAX(1)
C
      NFAX=IFAX(1)
      NX=N
      NH=N/2
      INK=INC+INC
      IF (ISIGN.EQ.+1) GO TO 30
C
C     IF NECESSARY, TRANSFER DATA TO WORK AREA
      IGO=50
      IF (MOD(NFAX,2).EQ.1) GOTO 40
      IBASE=1
      JBASE=1
      DO 20 L=1,LOT
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 10 M=1,N
      WORK(J)=A(I)
      I=I+INC
      J=J+1
   10 CONTINUE
      IBASE=IBASE+JUMP
      JBASE=JBASE+NX
   20 CONTINUE
C
      IGO=60
      GO TO 40
C
C     PREPROCESSING (ISIGN=+1)
C     ------------------------
C
   30 CONTINUE
      CALL FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
      IGO=60
C
C     COMPLEX*16 TRANSFORM
C     -----------------
C
   40 CONTINUE
      IA=1
      LA=1
      DO 80 K=1,NFAX
      IF (IGO.EQ.60) GO TO 60
   50 CONTINUE
      CALL VPASSM(A(IA),A(IA+INC),WORK(1),WORK(2),TRIGS,
     *   INK,2,JUMP,NX,LOT,NH,IFAX(K+1),LA)
      IGO=60
      GO TO 70
   60 CONTINUE
      CALL VPASSM(WORK(1),WORK(2),A(IA),A(IA+INC),TRIGS,
     *    2,INK,NX,JUMP,LOT,NH,IFAX(K+1),LA)
      IGO=50
   70 CONTINUE
      LA=LA*IFAX(K+1)
   80 CONTINUE
C
      IF (ISIGN.EQ.-1) GO TO 130
C
C     IF NECESSARY, TRANSFER DATA FROM WORK AREA
      IF (MOD(NFAX,2).EQ.1) GO TO 110
      IBASE=1
      JBASE=1
      DO 100 L=1,LOT
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 90 M=1,N
      A(J)=WORK(I)
      I=I+1
      J=J+INC
   90 CONTINUE
      IBASE=IBASE+NX
      JBASE=JBASE+JUMP
  100 CONTINUE
C
C     FILL IN ZEROS AT END
  110 CONTINUE
      GO TO 140
C
C     POSTPROCESSING (ISIGN=-1):
C     --------------------------
C
  130 CONTINUE
      CALL FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
C
  140 CONTINUE
      RETURN
      END
C     SUBROUTINE FFT99A - PREPROCESSING STEP FOR FFT99, ISIGN=+1
C     (SPECTRAL TO GRIDPOINT TRANSFORM)
C
      SUBROUTINE FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A(*),WORK(*),TRIGS(*)
      NH=N/2
      NX=N
      INK=INC+INC
C
C     A(0)   A(N/2)
      IA=1
      IB=N*INC+1
      JA=1
      JB=2
CDIR$ IVDEP
      DO 10 L=1,LOT
      WORK(JA)=A(IA)
      WORK(JB)=A(IA)
      IA=IA+JUMP
      IB=IB+JUMP
      JA=JA+NX
      JB=JB+NX
   10 CONTINUE
C
C     REMAINING WAVENUMBERS
      IABASE=2*INC+1
      IBBASE=(N-2)*INC+1
      JABASE=3
      JBBASE=N-1
C
      DO 30 K=3,NH,2
      IA=IABASE
      IB=IBBASE
      JA=JABASE
      JB=JBBASE
      C=TRIGS(N+K)
      S=TRIGS(N+K+1)
CDIR$ IVDEP
      DO 20 L=1,LOT
      WORK(JA)=(A(IA)+A(IB))-
     *    (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
      WORK(JB)=(A(IA)+A(IB))+
     *    (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
      WORK(JA+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))+
     *    (A(IA+INC)-A(IB+INC))
      WORK(JB+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))-
     *    (A(IA+INC)-A(IB+INC))
      IA=IA+JUMP
      IB=IB+JUMP
      JA=JA+NX
      JB=JB+NX
   20 CONTINUE
      IABASE=IABASE+INK
      IBBASE=IBBASE-INK
      JABASE=JABASE+2
      JBBASE=JBBASE-2
   30 CONTINUE
C
      IF (IABASE.NE.IBBASE) GO TO 50
C     WAVENUMBER N/4 (IF IT EXISTS)
      IA=IABASE
      JA=JABASE
CDIR$ IVDEP
      DO 40 L=1,LOT
      WORK(JA)=2.0D0*A(IA)
      WORK(JA+1)=-2.0D0*A(IA+INC)
      IA=IA+JUMP
      JA=JA+NX
   40 CONTINUE
C
   50 CONTINUE
      RETURN
      END
C     SUBROUTINE FFT99B - POSTPROCESSING STEP FOR FFT99, ISIGN=-1
C     (GRIDPOINT TO SPECTRAL TRANSFORM)

      SUBROUTINE FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WORK(*),A(*),TRIGS(*)
C
      NH=N/2
      NX=N
      INK=INC+INC
C
C     A(0)   A(N/2)
      SCALE=(1.0D0)/DBLE(N)
      IA=1
      IB=2
      JA=1
      JB=N*INC+1
CDIR$ IVDEP
      DO 10 L=1,LOT
      A(JA)=SCALE*(WORK(IA)+WORK(IB))
      A(JA+INC)=0.0D0
      IA=IA+NX
      IB=IB+NX
      JA=JA+JUMP
      JB=JB+JUMP
   10 CONTINUE
C
C     REMAINING WAVENUMBERS
      SCALE=0.5D0*SCALE
      IABASE=3
      IBBASE=N-1
      JABASE=2*INC+1
      JBBASE=(N-2)*INC+1
C
      DO 30 K=3,NH,2
      IA=IABASE
      IB=IBBASE
      JA=JABASE
      JB=JBBASE
      C=TRIGS(N+K)
      S=TRIGS(N+K+1)
CDIR$ IVDEP
      DO 20 L=1,LOT
      A(JA)=SCALE*((WORK(IA)+WORK(IB))
     *   +(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
      A(JB)=SCALE*((WORK(IA)+WORK(IB))
     *   -(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
      A(JA+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1)))
     *    +(WORK(IB+1)-WORK(IA+1)))
      A(JB+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1)))
     *    -(WORK(IB+1)-WORK(IA+1)))
      IA=IA+NX
      IB=IB+NX
      JA=JA+JUMP
      JB=JB+JUMP
   20 CONTINUE
      IABASE=IABASE+2
      IBBASE=IBBASE-2
      JABASE=JABASE+INK
      JBBASE=JBBASE-INK
   30 CONTINUE
C
      IF (IABASE.NE.IBBASE) GO TO 50
C     WAVENUMBER N/4 (IF IT EXISTS)
      IA=IABASE
      JA=JABASE
      SCALE=2.0D0*SCALE
CDIR$ IVDEP
      DO 40 L=1,LOT
      A(JA)=SCALE*WORK(IA)
      A(JA+INC)=-SCALE*WORK(IA+1)
      IA=IA+NX
      JA=JA+JUMP
   40 CONTINUE
C
   50 CONTINUE
      RETURN
      END
