C     MAIN PROGRAM
      INTEGER LUNIT
C     ALLOW 5000 UNDERFLOWS, 2 DIVISIONS BY ZERO.
C     THE DIVISIONS BY ZERO SHOULD NOT OCCUR IF COMPLEX DIVISION IS OK.
      CALL TRAPS(0,0,5001,0,3)
C
C     OUTPUT UNIT NUMBER
C
      LUNIT = 6
C
      CALL CSVTS(LUNIT)
      STOP
      END
      SUBROUTINE CSVTS(LUNIT)
C     LUNIT IS THE OUTPUT UNIT NUMBER.
C
C     TESTS
C        CSVDC
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     EXTERNAL CMACH,CSVT1,CXGEN
C     FORTRAN CMPLX,FLOAT
C
C     INTERNAL VARIABLES
C
      INTEGER LUNIT
      INTEGER I,J,N,P,LDX,LDU,LDV,CASE,NCASE
      COMPLEX X(25,25),XX(25,25),U(25,25),V(25,25),S(25),E(25),WORK(25)
      REAL CMACH,HUGE,TINY
      LOGICAL NOTWRT
      LDU = 25
      LDV = 25
      LDX = 25
      HUGE = CMACH(3)
      TINY = CMACH(2)
      NOTWRT = .TRUE.
      NCASE = 12
      WRITE (LUNIT,430)
      DO 290 CASE = 1, NCASE
         WRITE (LUNIT,300) CASE
         GO TO (10, 40, 70, 90, 110, 130, 170, 210, 240, 250, 260,
     *          270), CASE
C
C        BIDIAGONAL MATRIX WITH ZERO AT END.
C
   10    CONTINUE
            WRITE (LUNIT,310)
            N = 4
            P = 4
            DO 30 I = 1, 4
               DO 20 J = 1, 4
                  X(I,J) = (0.0E0,0.0E0)
   20          CONTINUE
   30       CONTINUE
            X(1,1) = (1.0E0,1.0E0)
            X(1,2) = (1.0E0,1.0E0)
            X(2,2) = (2.0E0,2.0E0)
            X(2,3) = (1.0E0,1.0E0)
            X(3,3) = (3.0E0,3.0E0)
            X(3,4) = (1.0E0,1.0E0)
            CALL CSVT1(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,XX,CASE,NOTWRT,
     *                 LUNIT,11)
         GO TO 280
C
C        BIDIAGONAL MATRIX WITH ZERO IN THE MIDDLE.
C
   40    CONTINUE
            WRITE (LUNIT,320)
            N = 5
            P = 5
            DO 60 I = 1, 5
               DO 50 J = 1, 5
                  X(I,J) = (0.0E0,0.0E0)
   50          CONTINUE
   60       CONTINUE
            X(1,1) = (1.0E0,1.0E0)
            X(1,2) = (1.0E0,1.0E0)
            X(2,3) = (1.0E0,1.0E0)
            X(3,3) = (2.0E0,2.0E0)
            X(3,4) = (1.0E0,1.0E0)
            X(4,4) = (3.0E0,3.0E0)
            X(4,5) = (1.0E0,1.0E0)
            X(5,5) = (4.0E0,4.0E0)
            CALL CSVT1(X,LDX,N,P,S,E,U,LDU,X,LDX,WORK,XX,CASE,NOTWRT,
     *                 LUNIT,11)
         GO TO 280
C
C        TEST CASE WITH N .GT. P.
C
   70    CONTINUE
            WRITE (LUNIT,330)
            N = 8
            P = 4
            DO 80 I = 1, 4
               S(I) = CMPLX(FLOAT(I),0.0E0)
   80       CONTINUE
            CALL CXGEN(X,LDX,N,P,S)
            CALL CSVT1(X,LDX,N,P,S,E,X,LDX,V,LDV,WORK,XX,CASE,NOTWRT,
     *                 LUNIT,21)
         GO TO 280
C
C        TEST CASE WITH N .LT. P.
C
   90    CONTINUE
            WRITE (LUNIT,340)
            N = 4
            P = 8
            DO 100 I = 1, 8
               S(I) = CMPLX(FLOAT(I),0.0E0)
  100       CONTINUE
            CALL CXGEN(X,LDX,N,P,S)
            CALL CSVT1(X,LDX,N,P,S,E,X,LDX,V,LDV,WORK,XX,CASE,NOTWRT,
     *                 LUNIT,11)
         GO TO 280
C
C        TEST CASE WITH N = P = LDX = LDU = LDV.
C
  110    CONTINUE
            WRITE (LUNIT,350)
            N = 25
            P = 25
            DO 120 I = 1, 25
               S(I) = CMPLX(FLOAT(I),0.0E0)
  120       CONTINUE
            CALL CXGEN(X,LDX,N,P,S)
            CALL CSVT1(X,LDX,N,P,S,E,X,LDX,V,LDV,WORK,XX,CASE,NOTWRT,
     *                 LUNIT,11)
         GO TO 280
C
C        TEST FOR OVERFLOW CONTROL.
C
  130    CONTINUE
            WRITE (LUNIT,360)
            N = 4
            P = 8
            DO 140 I = 1, 8
               S(I) = CMPLX(FLOAT(I),0.0E0)
  140       CONTINUE
            CALL CXGEN(X,LDX,N,P,S)
            DO 160 I = 1, 4
               DO 150 J = 1, 8
                  X(I,J) = HUGE*X(I,J)
  150          CONTINUE
  160       CONTINUE
            CALL CSVT1(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,XX,CASE,NOTWRT,
     *                 LUNIT,11)
         GO TO 280
C
C        TEST FOR UNDERFLOW CONTROL.
C
  170    CONTINUE
            WRITE (LUNIT,370)
            N = 8
            P = 4
            DO 180 I = 1, 8
               S(I) = CMPLX(FLOAT(I),0.0E0)
  180       CONTINUE
            CALL CXGEN(X,LDX,N,P,S)
            DO 200 I = 1, 8
               DO 190 J = 1, 4
                  X(I,J) = TINY*X(I,J)
  190          CONTINUE
  200       CONTINUE
            CALL CSVT1(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,XX,CASE,NOTWRT,
     *                 LUNIT,11)
         GO TO 280
C
C        ZERO MATRIX.
C
  210    CONTINUE
            WRITE (LUNIT,380)
            N = 8
            P = 4
            DO 230 I = 1, N
               DO 220 J = 1, P
                  X(I,J) = (0.0E0,0.0E0)
  220          CONTINUE
  230       CONTINUE
            CALL CSVT1(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,XX,CASE,NOTWRT,
     *                 LUNIT,11)
         GO TO 280
C
C        1X1 MATRIX.
C
  240    CONTINUE
            WRITE (LUNIT,390)
            N = 1
            P = 1
            X(1,1) = (2.0E0,0.0E0)
            CALL CSVT1(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,XX,CASE,NOTWRT,
     *                 LUNIT,11)
         GO TO 280
C
C        2X2 MATRIX.
C
  250    CONTINUE
            WRITE (LUNIT,400)
            N = 2
            P = 2
            X(1,1) = (3.0E0,0.0E0)
            X(1,2) = (1.0E0,0.0E0)
            X(2,1) = (1.0E0,0.0E0)
            X(2,2) = (2.0E0,0.0E0)
            CALL CSVT1(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,XX,CASE,NOTWRT,
     *                 LUNIT,11)
         GO TO 280
C
C        COLUMN VECTOR.
C
  260    CONTINUE
            WRITE (LUNIT,410)
            N = 4
            P = 1
            X(1,1) = (1.0E0,0.0E0)
            X(2,1) = (0.0E0,0.0E0)
            X(3,1) = (0.0E0,0.0E0)
            X(4,1) = (2.0E0,0.0E0)
            CALL CSVT1(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,XX,CASE,NOTWRT,
     *                 LUNIT,11)
         GO TO 280
C
C        ROW VECTOR.
C
  270    CONTINUE
            WRITE (LUNIT,420)
            N = 1
            P = 4
            X(1,1) = (0.0E0,0.0E0)
            X(1,2) = (1.0E0,0.0E0)
            X(1,3) = (2.0E0,0.0E0)
            X(1,4) = (3.0E0,0.0E0)
            CALL CSVT1(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,XX,CASE,NOTWRT,
     *                 LUNIT,11)
  280    CONTINUE
  290 CONTINUE
      WRITE (LUNIT,440)
      RETURN
  300 FORMAT ( / 5H1CASE, I3)
  310 FORMAT ( / 35H BIDIAGONAL MATRIX WITH ZERO AT END)
  320 FORMAT ( / 38H BIDIAGONAL MATRIX WITH ZERO IN MIDDLE)
  330 FORMAT ( / 13H 8 X 4 MATRIX)
  340 FORMAT ( / 13H 4 X 8 MATRIX)
  350 FORMAT ( / 15H 25 X 25 MATRIX)
  360 FORMAT ( / 14H OVERFLOW TEST)
  370 FORMAT ( / 15H UNDERFLOW TEST)
  380 FORMAT ( / 12H ZERO MATRIX)
  390 FORMAT ( / 13H 1 X 1 MATRIX)
  400 FORMAT ( / 13H 2 X 2 MATRIX)
  410 FORMAT ( / 14H COLUMN VECTOR)
  420 FORMAT ( / 11H ROW VECTOR)
  430 FORMAT (22H1LINPACK TESTER, CSV** /
     *        29H THIS VERSION DATED 08/14/78.)
  440 FORMAT ( / 27H1END OF SINGULAR VALUE TEST)
      END
      SUBROUTINE CARRAY(A,LDA,M,N,NNL,LUNIT)
      INTEGER LDA,M,N,NNL,LUNIT
      COMPLEX A(LDA,1)
C
C     FORTRAN IABS,MIN0
C
      INTEGER I,J,K,KU,NL
      NL = IABS(NNL)
      IF (NNL .LT. 0) GO TO 30
         DO 20 I = 1, M
            WRITE (6,70)
            DO 10 K = 1, N, NL
               KU = MIN0(K+NL-1,N)
               WRITE (6,70) (A(I,J), J = K, KU)
   10       CONTINUE
   20    CONTINUE
      GO TO 60
   30 CONTINUE
         DO 50 J = 1, N
            WRITE (6,70)
            DO 40 K = 1, M, NL
               KU = MIN0(K+NL-1,M)
               WRITE (6,70) (A(I,J), I = K, KU)
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
      RETURN
   70 FORMAT (1H , 8E13.6)
      END
      SUBROUTINE CSVBM(XX,LDX,N,P,S,U,LDU,V,LDV,X,XSTAT)
      INTEGER LDX,N,P,LDU,LDV
      COMPLEX XX(LDX,1),S(1),U(LDU,1),V(LDV,1),X(LDX,1)
      REAL XSTAT
C
C     EXTERNAL CMACH
C     FORTRAN AMAX1,CABS,CONJG,MIN0
C
      INTEGER I,J,K,M
      COMPLEX T(25)
      REAL CMACH,EMAX,XMAX
C
      M = MIN0(N,P)
      DO 20 J = 1, P
         DO 10 I = 1, M
            X(I,J) = S(I)*CONJG(V(J,I))
   10    CONTINUE
   20 CONTINUE
      IF (N .LE. P) GO TO 50
         M = P + 1
         DO 40 J = 1, P
            DO 30 I = M, N
               X(I,J) = (0.0E0,0.0E0)
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      M = MIN0(N,P)
      DO 90 J = 1, P
         DO 70 I = 1, N
            T(I) = (0.0E0,0.0E0)
            DO 60 K = 1, M
               T(I) = T(I) + U(I,K)*X(K,J)
   60       CONTINUE
   70    CONTINUE
         DO 80 I = 1, N
            X(I,J) = T(I)
   80    CONTINUE
   90 CONTINUE
      EMAX = 0.0E0
      XMAX = 0.0E0
      DO 110 J = 1, P
         DO 100 I = 1, N
            EMAX = AMAX1(EMAX,CABS(X(I,J)-XX(I,J)))
            XMAX = AMAX1(XMAX,CABS(XX(I,J)))
  100    CONTINUE
  110 CONTINUE
      XSTAT = 0.0E0
      IF (EMAX .EQ. 0.0E0) GO TO 140
         IF (XMAX .EQ. 0.0E0) GO TO 120
            XSTAT = (EMAX/XMAX)/CMACH(1)
         GO TO 130
  120    CONTINUE
            XSTAT = 1.0E10
  130    CONTINUE
  140 CONTINUE
      RETURN
      END
      SUBROUTINE CSVOT(X,LDX,N,COL,TEST)
      INTEGER LDX,N,COL
      COMPLEX X(LDX,1)
C
C     FORTRAN AMAX1
C
      INTEGER I,J,K
      COMPLEX ELM
      REAL TEST,EMAX,CMACH
      EMAX = 0.0E0
      DO 30 I = 1, COL
         DO 20 J = 1, COL
            ELM = (0.0E0,0.0E0)
            DO 10 K = 1, N
               ELM = ELM + CONJG(X(K,I))*X(K,J)
   10       CONTINUE
            IF (I .EQ. J) ELM = (1.0E0,0.0E0) - ELM
            EMAX = AMAX1(EMAX,CABS(ELM))
   20    CONTINUE
   30 CONTINUE
      TEST = EMAX/CMACH(1)
      RETURN
      END
      SUBROUTINE CSVT1(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,XX,CASE,NOTWRT,
     *                 LUNIT,JOB)
      INTEGER LDX,N,P,LDU,LDV,CASE,LUNIT,JOB
      COMPLEX X(LDX,1),S(1),E(1),U(LDU,1),V(LDV,1),WORK(1),XX(LDX,1)
      LOGICAL NOTWRT
C
C     EXTERNAL CARRAY,CSVBM,CSVDC,CSVOT
C     FORTRAN AMAX1,MIN0
C
      INTEGER I,J,M,UCOL,INFO
      REAL USTAT,VSTAT,XSTAT
      COMPLEX XNEW(25,25)
      REAL ERRLVL
      ERRLVL = 100.0E0
      UCOL = N
      IF (JOB/10 .GE. 2) UCOL = P
      DO 20 J = 1, P
         DO 10 I = 1, N
            XX(I,J) = X(I,J)
   10    CONTINUE
   20 CONTINUE
      IF (NOTWRT) GO TO 30
         WRITE (LUNIT,50)
         CALL CARRAY(X,LDX,N,P,5,LUNIT)
   30 CONTINUE
      CALL CSVDC(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,JOB,INFO)
      M = MIN0(N,P)
      IF (NOTWRT) GO TO 40
         WRITE (LUNIT,60)
         CALL CARRAY(S,P,M,1,-5,LUNIT)
         WRITE (LUNIT,70)
         CALL CARRAY(U,LDU,N,N,5,LUNIT)
         WRITE (LUNIT,80)
         CALL CARRAY(V,LDV,P,P,5,LUNIT)
   40 CONTINUE
      CALL CSVOT(U,LDU,N,UCOL,USTAT)
      CALL CSVOT(V,LDV,P,P,VSTAT)
      CALL CSVBM(XX,LDX,N,P,S,U,LDU,V,LDV,XNEW,XSTAT)
      WRITE (LUNIT,90) XSTAT,USTAT,VSTAT
      IF (AMAX1(XSTAT,USTAT,VSTAT) .GT. ERRLVL) WRITE (LUNIT,100)
      RETURN
   50 FORMAT ( / 2H X)
   60 FORMAT ( / 2H S)
   70 FORMAT ( / 2H U)
   80 FORMAT ( / 2H V)
   90 FORMAT ( / 11H STATISTICS //
     *         38H         U*SIGMA*VH .................., E10.2 /
     *         38H         UHU ........................., E10.2 /
     *         38H         VHV ........................., E10.2)
  100 FORMAT ( / 35H ***** STATISTICS ABOVE ERROR LEVEL)
      END
      SUBROUTINE CXGEN(X,LDX,N,P,S)
      INTEGER P,N,LDX
      COMPLEX X(LDX,1),S(1)
C
C     FORTRAN CMPLX,FLOAT,MIN0
C
      INTEGER I,J,M,MP1
      COMPLEX T,RU,FAC
      REAL FP,FN
      FP = FLOAT(P)
      FN = FLOAT(N)
      M = MIN0(N,P)
      RU = CMPLX(COS(6.28E0/FLOAT(M+1)),SIN(6.28E0/FLOAT(M+1)))
      FAC = (1.0E0,0.0E0)/FP
      DO 20 I = 1, M
         FAC = FAC*RU
         DO 10 J = 1, P
            X(I,J) = -2.0E0*S(I)*FAC
   10    CONTINUE
         X(I,I) = X(I,I) + FP*S(I)*FAC
   20 CONTINUE
      IF (M .GE. N) GO TO 50
         MP1 = M + 1
         DO 40 J = 1, P
            DO 30 I = MP1, N
               X(I,J) = (0.0E0,0.0E0)
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      DO 80 J = 1, P
         T = (0.0E0,0.0E0)
         DO 60 I = 1, N
            T = T + X(I,J)
   60    CONTINUE
         DO 70 I = 1, N
            X(I,J) = X(I,J) - 2.0E0*T/FN
   70    CONTINUE
   80 CONTINUE
      RETURN
      END
