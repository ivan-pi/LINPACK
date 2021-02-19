C     MAIN PROGRAM
      INTEGER LUNIT
C     ALLOW 5000 UNDERFLOWS.
      CALL TRAPS(0,0,5001,0,0)
C
C     OUTPUT UNIT NUMBER
C
      LUNIT = 6
C
      CALL DSVTS(LUNIT)
      STOP
      END
      SUBROUTINE DSVTS(LUNIT)
C     LUNIT IS THE OUTPUT UNIT NUMBER.
C
C     TESTS
C        DSVDC
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     EXTERNAL DMACH,DSVT1,DXGEN
C     FORTRAN DBLE,FLOAT
C
C     INTERNAL VARIABLES
C
      INTEGER LUNIT
      INTEGER I,J,N,P,LDX,LDU,LDV,CASE,NCASE
      DOUBLE PRECISION X(25,25),XX(25,25),U(25,25),V(25,25),S(25),
     *                 E(25),WORK(25)
      DOUBLE PRECISION DMACH,HUGE,TINY
      LOGICAL NOTWRT
      LDU = 25
      LDV = 25
      LDX = 25
      HUGE = DMACH(3)
      TINY = DMACH(2)
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
                  X(I,J) = 0.0D0
   20          CONTINUE
   30       CONTINUE
            X(1,1) = 1.0D0
            X(1,2) = 1.0D0
            X(2,2) = 2.0D0
            X(2,3) = 1.0D0
            X(3,3) = 3.0D0
            X(3,4) = 1.0D0
            CALL DSVT1(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,XX,CASE,NOTWRT,
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
                  X(I,J) = 0.0D0
   50          CONTINUE
   60       CONTINUE
            X(1,1) = 1.0D0
            X(1,2) = 1.0D0
            X(2,3) = 1.0D0
            X(3,3) = 2.0D0
            X(3,4) = 1.0D0
            X(4,4) = 3.0D0
            X(4,5) = 1.0D0
            X(5,5) = 4.0D0
            CALL DSVT1(X,LDX,N,P,S,E,U,LDU,X,LDX,WORK,XX,CASE,NOTWRT,
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
               S(I) = DBLE(FLOAT(I))
   80       CONTINUE
            CALL DXGEN(X,LDX,N,P,S)
            CALL DSVT1(X,LDX,N,P,S,E,X,LDX,V,LDV,WORK,XX,CASE,NOTWRT,
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
               S(I) = DBLE(FLOAT(I))
  100       CONTINUE
            CALL DXGEN(X,LDX,N,P,S)
            CALL DSVT1(X,LDX,N,P,S,E,X,LDX,V,LDV,WORK,XX,CASE,NOTWRT,
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
               S(I) = DBLE(FLOAT(I))
  120       CONTINUE
            CALL DXGEN(X,LDX,N,P,S)
            CALL DSVT1(X,LDX,N,P,S,E,X,LDX,V,LDV,WORK,XX,CASE,NOTWRT,
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
               S(I) = DBLE(FLOAT(I))
  140       CONTINUE
            CALL DXGEN(X,LDX,N,P,S)
            DO 160 I = 1, 4
               DO 150 J = 1, 8
                  X(I,J) = HUGE*X(I,J)
  150          CONTINUE
  160       CONTINUE
            CALL DSVT1(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,XX,CASE,NOTWRT,
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
               S(I) = DBLE(FLOAT(I))
  180       CONTINUE
            CALL DXGEN(X,LDX,N,P,S)
            DO 200 I = 1, 8
               DO 190 J = 1, 4
                  X(I,J) = TINY*X(I,J)
  190          CONTINUE
  200       CONTINUE
            CALL DSVT1(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,XX,CASE,NOTWRT,
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
                  X(I,J) = 0.0D0
  220          CONTINUE
  230       CONTINUE
            CALL DSVT1(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,XX,CASE,NOTWRT,
     *                 LUNIT,11)
         GO TO 280
C
C        1X1 MATRIX.
C
  240    CONTINUE
            WRITE (LUNIT,390)
            N = 1
            P = 1
            X(1,1) = 2.0D0
            CALL DSVT1(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,XX,CASE,NOTWRT,
     *                 LUNIT,11)
         GO TO 280
C
C        2X2 MATRIX.
C
  250    CONTINUE
            WRITE (LUNIT,400)
            N = 2
            P = 2
            X(1,1) = 3.0D0
            X(1,2) = 1.0D0
            X(2,1) = 1.0D0
            X(2,2) = 2.0D0
            CALL DSVT1(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,XX,CASE,NOTWRT,
     *                 LUNIT,11)
         GO TO 280
C
C        COLUMN VECTOR.
C
  260    CONTINUE
            WRITE (LUNIT,410)
            N = 4
            P = 1
            X(1,1) = 1.0D0
            X(2,1) = 0.0D0
            X(3,1) = 0.0D0
            X(4,1) = 2.0D0
            CALL DSVT1(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,XX,CASE,NOTWRT,
     *                 LUNIT,11)
         GO TO 280
C
C        ROW VECTOR.
C
  270    CONTINUE
            WRITE (LUNIT,420)
            N = 1
            P = 4
            X(1,1) = 0.0D0
            X(1,2) = 1.0D0
            X(1,3) = 2.0D0
            X(1,4) = 3.0D0
            CALL DSVT1(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,XX,CASE,NOTWRT,
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
  430 FORMAT (22H1LINPACK TESTER, DSV** /
     *        29H THIS VERSION DATED 08/14/78.)
  440 FORMAT ( / 27H1END OF SINGULAR VALUE TEST)
      END
      SUBROUTINE DARRAY(A,LDA,M,N,NNL,LUNIT)
      INTEGER LDA,M,N,NNL,LUNIT
      DOUBLE PRECISION A(LDA,1)
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
   70 FORMAT (1H , 8D13.6)
      END
      SUBROUTINE CSVBM(XX,LDX,N,P,S,U,LDU,V,LDV,X,XSTAT)
      INTEGER LDX,N,P,LDU,LDV
      DOUBLE PRECISION XX(LDX,1),S(1),U(LDU,1),V(LDV,1),X(LDX,1)
      DOUBLE PRECISION XSTAT
C
C     EXTERNAL DMACH
C     FORTRAN DMAX1,DABS,MIN0
C
      INTEGER I,J,K,M
      DOUBLE PRECISION T(25)
      DOUBLE PRECISION DMACH,EMAX,XMAX
C
      M = MIN0(N,P)
      DO 20 J = 1, P
         DO 10 I = 1, M
            X(I,J) = S(I)*V(J,I)
   10    CONTINUE
   20 CONTINUE
      IF (N .LE. P) GO TO 50
         M = P + 1
         DO 40 J = 1, P
            DO 30 I = M, N
               X(I,J) = 0.0D0
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      M = MIN0(N,P)
      DO 90 J = 1, P
         DO 70 I = 1, N
            T(I) = 0.0D0
            DO 60 K = 1, M
               T(I) = T(I) + U(I,K)*X(K,J)
   60       CONTINUE
   70    CONTINUE
         DO 80 I = 1, N
            X(I,J) = T(I)
   80    CONTINUE
   90 CONTINUE
      EMAX = 0.0D0
      XMAX = 0.0D0
      DO 110 J = 1, P
         DO 100 I = 1, N
            EMAX = DMAX1(EMAX,DABS(X(I,J)-XX(I,J)))
            XMAX = DMAX1(XMAX,DABS(XX(I,J)))
  100    CONTINUE
  110 CONTINUE
      XSTAT = 0.0D0
      IF (EMAX .EQ. 0.0D0) GO TO 140
         IF (XMAX .EQ. 0.0D0) GO TO 120
            XSTAT = (EMAX/XMAX)/DMACH(1)
         GO TO 130
  120    CONTINUE
            XSTAT = 1.0D10
  130    CONTINUE
  140 CONTINUE
      RETURN
      END
      SUBROUTINE DSVOT(X,LDX,N,COL,TEST)
      INTEGER LDX,N,COL
      DOUBLE PRECISION X(LDX,1)
C
C     FORTRAN DMAX1
C
      INTEGER I,J,K
      DOUBLE PRECISION ELM
      DOUBLE PRECISION TEST,EMAX,DMACH
      EMAX = 0.0D0
      DO 30 I = 1, COL
         DO 20 J = 1, COL
            ELM = 0.0D0
            DO 10 K = 1, N
               ELM = ELM + X(K,I)*X(K,J)
   10       CONTINUE
            IF (I .EQ. J) ELM = 1.0D0 - ELM
            EMAX = DMAX1(EMAX,DABS(ELM))
   20    CONTINUE
   30 CONTINUE
      TEST = EMAX/DMACH(1)
      RETURN
      END
      SUBROUTINE DSVT1(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,XX,CASE,NOTWRT,
     *                 LUNIT,JOB)
      INTEGER LDX,N,P,LDU,LDV,CASE,LUNIT,JOB
      DOUBLE PRECISION X(LDX,1),S(1),E(1),U(LDU,1),V(LDV,1),WORK(1),
     *                 XX(LDX,1)
      LOGICAL NOTWRT
C
C     EXTERNAL DARRAY,DSVBM,DSVDC,DSVOT
C     FORTRAN DMAX1,MIN0
C
      INTEGER I,J,M,UCOL,INFO
      DOUBLE PRECISION USTAT,VSTAT,XSTAT
      DOUBLE PRECISION XNEW(25,25)
      DOUBLE PRECISION ERRLVL
      ERRLVL = 100.0D0
      UCOL = N
      IF (JOB/10 .GE. 2) UCOL = P
      DO 20 J = 1, P
         DO 10 I = 1, N
            XX(I,J) = X(I,J)
   10    CONTINUE
   20 CONTINUE
      IF (NOTWRT) GO TO 30
         WRITE (LUNIT,50)
         CALL DARRAY(X,LDX,N,P,5,LUNIT)
   30 CONTINUE
      CALL DSVDC(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,JOB,INFO)
      M = MIN0(N,P)
      IF (NOTWRT) GO TO 40
         WRITE (LUNIT,60)
         CALL DARRAY(S,P,M,1,-5,LUNIT)
         WRITE (LUNIT,70)
         CALL DARRAY(U,LDU,N,N,5,LUNIT)
         WRITE (LUNIT,80)
         CALL DARRAY(V,LDV,P,P,5,LUNIT)
   40 CONTINUE
      CALL DSVOT(U,LDU,N,UCOL,USTAT)
      CALL DSVOT(V,LDV,P,P,VSTAT)
      CALL CSVBM(XX,LDX,N,P,S,U,LDU,V,LDV,XNEW,XSTAT)
      WRITE (LUNIT,90) XSTAT,USTAT,VSTAT
      IF (DMAX1(XSTAT,USTAT,VSTAT) .GT. ERRLVL) WRITE (LUNIT,100)
      RETURN
   50 FORMAT ( / 2H X)
   60 FORMAT ( / 2H S)
   70 FORMAT ( / 2H U)
   80 FORMAT ( / 2H V)
   90 FORMAT ( / 11H STATISTICS //
     *         38H         U*SIGMA*VH .................., D10.2 /
     *         38H         UHU ........................., D10.2 /
     *         38H         VHV ........................., D10.2)
  100 FORMAT ( / 35H ***** STATISTICS ABOVE ERROR LEVEL)
      END
      SUBROUTINE DXGEN(X,LDX,N,P,S)
      INTEGER P,N,LDX
      DOUBLE PRECISION X(LDX,1),S(1)
C
C     FORTRAN DBLE,FLOAT,MIN0
C
      INTEGER I,J,M,MP1
      DOUBLE PRECISION T,RU,FAC
      DOUBLE PRECISION FP,FN
      FP = DBLE(FLOAT(P))
      FN = DBLE(FLOAT(N))
      M = MIN0(N,P)
      RU = DCOS(6.28D0/DBLE(FLOAT(M+1)))
      FAC = 1.0D0/FP
      DO 20 I = 1, M
         FAC = FAC*RU
         DO 10 J = 1, P
            X(I,J) = -2.0D0*S(I)*FAC
   10    CONTINUE
         X(I,I) = X(I,I) + FP*S(I)*FAC
   20 CONTINUE
      IF (M .GE. N) GO TO 50
         MP1 = M + 1
         DO 40 J = 1, P
            DO 30 I = MP1, N
               X(I,J) = 0.0D0
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      DO 80 J = 1, P
         T = 0.0D0
         DO 60 I = 1, N
            T = T + X(I,J)
   60    CONTINUE
         DO 70 I = 1, N
            X(I,J) = X(I,J) - 2.0D0*T/FN
   70    CONTINUE
   80 CONTINUE
      RETURN
      END
