C     MAIN PROGRAM
      INTEGER LUNIT
C     ALLOW 5000 UNDERFLOWS.
      CALL TRAPS(0,0,5001,0,0)
C
C     OUTPUT UNIT NUMBER
C
      LUNIT = 6
C
      CALL SQRTS(LUNIT)
      STOP
      END
      SUBROUTINE SQRTS(LUNIT)
C     LUNIT IS THE OUTPUT UNIT NUMBER
C
C     TESTS
C        SQRDC,SQRSL
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK SQRDC,SQRSL
C     EXTERNAL SQRBM,SQRFM,SQRRX,SXGEN,SMACH
C     FORTRAN AMAX1,ABS,FLOAT
C
C     INTERNAL VARIABLES
C
      INTEGER LUNIT
      INTEGER CASE,LDX,N,P,PD2,PD2P1,NCASE,CONT,JPVT(25)
      INTEGER I,INFO,J,JJ,JJJ,L
      REAL BESTAT,RMAX,RSTAT,XBMAX,XBSTAT,ERRLVL,BSTAT,FSTAT
      REAL BETA(25),QRAUX(25),QY(25),RSD(25),S(25),WORK(25)
      REAL X(25,25),XBETA(25),XX(25,25),Y(25),Y1(25),Y2(25)
      REAL SMACH,TINY,HUGE
      LOGICAL NOTWRT
      DATA CONT /1H /
      LDX = 25
      NCASE = 9
      ERRLVL = 100.0E0
      NOTWRT = .TRUE.
      TINY = SMACH(2)
      HUGE = SMACH(3)
      WRITE (LUNIT,600)
      DO 550 CASE = 1, NCASE
         WRITE (LUNIT,560) CASE
         XBSTAT = 0.0E0
         BESTAT = 0.0E0
         GO TO (10, 100, 170, 240, 300, 330, 380, 430, 470), CASE
C
C        WELL CONDITIONED LEAST SQUARES PROBLEM.
C
   10    CONTINUE
            WRITE (LUNIT,640)
            N = 10
            P = 4
C
C           GENERATE A MATRIX X WITH SINGULAR VALUES 1. AND .5.
C
            PD2 = P/2
            PD2P1 = PD2 + 1
            DO 20 L = 1, PD2
               S(L) = 1.0E0
   20       CONTINUE
            DO 30 L = PD2P1, P
               S(L) = 0.5E0
   30       CONTINUE
            CALL SXGEN(X,LDX,N,P,S)
            IF (NOTWRT) GO TO 40
               WRITE (LUNIT,610)
               CALL SARRAY(X,LDX,N,P,4,LUNIT)
   40       CONTINUE
C           THE OBSERVATION VECTOR IS Y = Y1 + Y2, WHERE Y1 IS
C           THE VECTOR OF ROW SUMS OF X AND Y2 IS ORTHOGONAL TO
C           THE COLUMNS OF X.
C
            DO 60 I = 1, N
               Y1(I) = 0.0E0
               Y2(I) = 2.0E0*TINY
               IF (I .EQ. P + 1) Y2(I) = Y2(I) - FLOAT(N)*TINY
               DO 50 J = 1, P
                  X(I,J) = X(I,J)*TINY
                  Y1(I) = Y1(I) + X(I,J)
                  XX(I,J) = X(I,J)
   50          CONTINUE
               Y(I) = Y1(I) + Y2(I)
   60       CONTINUE
C
C           REDUCE THE MATRIX.
C
            CALL SQRDC(X,LDX,N,P,QRAUX,JPVT,WORK,0)
            IF (NOTWRT) GO TO 70
               WRITE (LUNIT,610)
               CALL SARRAY(X,LDX,N,P,4,LUNIT)
               WRITE (LUNIT,620)
               CALL SARRAY(QRAUX,P,P,1,-4,LUNIT)
   70       CONTINUE
C
C           COMPUTE THE STATISTICS FOR THE FOWARD AND BACK
C           MULTIPLICATIONS.
C
            CALL SQRFM(X,LDX,N,P,XX,WORK,QRAUX,FSTAT)
            CALL SQRBM(X,LDX,N,P,XX,WORK,QRAUX,BSTAT)
C
C           SOLVE THE LEAST SQUARES PROBLEM.
C
            CALL SQRSL(X,LDX,N,P,QRAUX,Y,QY,RSD,BETA,RSD,XBETA,111,
     *                 INFO)
C
C           COMPUTE THE LEAST SQUARES TEST STATISTICS.
C
            RSTAT = 0.0E0
            RMAX = 0.0E0
            XBSTAT = 0.0E0
            XBMAX = 0.0E0
            DO 80 I = 1, N
               RSTAT = AMAX1(RSTAT,ABS(Y2(I)-RSD(I)))
               RMAX = AMAX1(RMAX,ABS(Y2(I)))
               XBSTAT = AMAX1(XBSTAT,ABS(Y1(I)-XBETA(I)))
               XBMAX = AMAX1(XBMAX,ABS(Y2(I)))
   80       CONTINUE
            RSTAT = RSTAT/(SMACH(1)*RMAX)
            XBSTAT = XBSTAT/(SMACH(1)*XBMAX)
            BESTAT = 0.0E0
            DO 90 J = 1, P
               BESTAT = AMAX1(BESTAT,ABS(1.0E0-BETA(J)))
   90       CONTINUE
            BESTAT = BESTAT/SMACH(1)
            WRITE (LUNIT,570) FSTAT,BSTAT,CONT,BESTAT,XBSTAT,RSTAT
         GO TO 540
C
C        4 X 10 MATRIX.
C
  100    CONTINUE
            WRITE (LUNIT,650)
            N = 4
            P = 10
            PD2 = P/2
            PD2P1 = PD2 + 1
            DO 110 L = 1, PD2
               S(L) = 1.0E0
  110       CONTINUE
            DO 120 L = PD2P1, P
               S(L) = 0.5E0
  120       CONTINUE
            CALL SXGEN(X,LDX,N,P,S)
            IF (NOTWRT) GO TO 130
               WRITE (LUNIT,610)
               CALL SARRAY(X,LDX,N,P,10,LUNIT)
  130       CONTINUE
            DO 150 I = 1, N
               DO 140 J = 1, P
                  XX(I,J) = X(I,J)
  140          CONTINUE
  150       CONTINUE
            CALL SQRDC(X,LDX,N,P,QRAUX,JPVT,WORK,0)
            IF (NOTWRT) GO TO 160
               WRITE (LUNIT,610)
               CALL SARRAY(X,LDX,N,P,10,LUNIT)
               WRITE (LUNIT,620)
               CALL SARRAY(QRAUX,N,N,1,-4,LUNIT)
  160       CONTINUE
            CALL SQRFM(X,LDX,N,P,XX,WORK,QRAUX,FSTAT)
            CALL SQRBM(X,LDX,N,P,XX,WORK,QRAUX,BSTAT)
            WRITE (LUNIT,570) FSTAT,BSTAT
         GO TO 540
C
C        PIVOTING AND OVERFLOW TEST.  COLUMNS 1, 4, 7 FROZEN.
C
  170    CONTINUE
            WRITE (LUNIT,670)
            N = 10
            P = 3
            S(1) = 1.0E0
            S(2) = 0.5E0
            S(3) = 0.25E0
            CALL SXGEN(X,LDX,N,P,S)
            P = 9
            DO 190 I = 1, N
               DO 180 JJJ = 1, 3
                  J = 4 - JJJ
                  JJ = 3*J
                  X(I,JJ) = HUGE*X(I,J)
                  X(I,JJ-1) = X(I,JJ)/2.0E0
                  X(I,JJ-2) = X(I,JJ-1)/2.0E0
  180          CONTINUE
  190       CONTINUE
            IF (NOTWRT) GO TO 200
               WRITE (LUNIT,610)
               CALL SARRAY(X,LDX,N,P,9,LUNIT)
  200       CONTINUE
            DO 220 J = 1, P
               JPVT(J) = 0
               DO 210 I = 1, N
                  XX(I,J) = X(I,J)
  210          CONTINUE
  220       CONTINUE
            JPVT(1) = -1
            JPVT(4) = -1
            JPVT(7) = -1
            CALL SQRDC(X,LDX,N,P,QRAUX,JPVT,WORK,1)
            IF (NOTWRT) GO TO 230
               WRITE (LUNIT,610)
               CALL SARRAY(X,LDX,N,P,9,LUNIT)
               WRITE (LUNIT,620)
               CALL SARRAY(QRAUX,P,P,1,-9,LUNIT)
  230       CONTINUE
            WRITE (LUNIT,630) (JPVT(J), J = 1, P)
            CALL SQRRX(XX,LDX,N,P,JPVT)
            CALL SQRFM(X,LDX,N,P,XX,WORK,QRAUX,FSTAT)
            CALL SQRBM(X,LDX,N,P,XX,WORK,QRAUX,BSTAT)
            WRITE (LUNIT,570) FSTAT,BSTAT
         GO TO 540
C
C        25 BY 25 MATRIX
C
  240    CONTINUE
            WRITE (LUNIT,680)
            N = 25
            P = 25
            DO 250 I = 1, 25
               S(I) = 2.0E0**(1 - I)
  250       CONTINUE
            CALL SXGEN(X,LDX,N,P,S)
            IF (NOTWRT) GO TO 260
               WRITE (LUNIT,610)
               CALL SARRAY(X,LDX,N,P,10,LUNIT)
  260       CONTINUE
            DO 280 J = 1, P
               JPVT(J) = 0
               DO 270 I = 1, N
                  XX(I,J) = X(I,J)
  270          CONTINUE
  280       CONTINUE
            CALL SQRDC(X,LDX,N,P,QRAUX,JPVT,WORK,1)
            IF (NOTWRT) GO TO 290
               WRITE (LUNIT,610)
               CALL SARRAY(X,LDX,N,P,10,LUNIT)
               WRITE (LUNIT,620)
               CALL SARRAY(QRAUX,P,P,1,-10,LUNIT)
  290       CONTINUE
            WRITE (LUNIT,630) (JPVT(J), J = 1, P)
            CALL SQRRX(XX,LDX,N,P,JPVT)
            CALL SQRFM(X,LDX,N,P,XX,WORK,QRAUX,FSTAT)
            CALL SQRBM(X,LDX,N,P,XX,WORK,QRAUX,BSTAT)
            WRITE (LUNIT,570) FSTAT,BSTAT
         GO TO 540
C
C        MONOELEMENTAL MATRIX.
C
  300    CONTINUE
            WRITE (LUNIT,690)
            N = 1
            P = 1
            X(1,1) = 1.0E0
            IF (NOTWRT) GO TO 310
               WRITE (LUNIT,610)
               CALL SARRAY(X,LDX,N,P,1,LUNIT)
  310       CONTINUE
            XX(1,1) = 1.0E0
            JPVT(1) = 0
            CALL SQRDC(X,LDX,N,P,QRAUX,JPVT,WORK,1)
            IF (NOTWRT) GO TO 320
               WRITE (LUNIT,610)
               CALL SARRAY(X,LDX,N,P,1,LUNIT)
               WRITE (LUNIT,620)
               CALL SARRAY(QRAUX,P,P,1,-1,LUNIT)
  320       CONTINUE
            CALL SQRFM(X,LDX,N,P,XX,WORK,QRAUX,FSTAT)
            CALL SQRBM(X,LDX,N,P,XX,WORK,QRAUX,BSTAT)
            WRITE (LUNIT,570) FSTAT,BSTAT
         GO TO 540
C
C        ZERO MATRIX.
C
  330    CONTINUE
            WRITE (LUNIT,700)
            N = 10
            P = 4
            DO 350 J = 1, P
               JPVT(J) = 0
               DO 340 I = 1, N
                  X(I,J) = 0.0E0
                  XX(I,J) = 0.0E0
  340          CONTINUE
  350       CONTINUE
            IF (NOTWRT) GO TO 360
               WRITE (LUNIT,610)
               CALL SARRAY(X,LDX,N,P,4,LUNIT)
  360       CONTINUE
            CALL SQRDC(X,LDX,N,P,QRAUX,JPVT,WORK,1)
            IF (NOTWRT) GO TO 370
               WRITE (LUNIT,610)
               CALL SARRAY(X,LDX,N,P,4,LUNIT)
               WRITE (LUNIT,620)
               CALL SARRAY(QRAUX,P,P,1,-4,LUNIT)
  370       CONTINUE
            WRITE (LUNIT,630) (JPVT(J), J = 1, P)
            CALL SQRRX(XX,LDX,N,P,JPVT)
            CALL SQRFM(X,LDX,N,P,XX,WORK,QRAUX,FSTAT)
            CALL SQRBM(X,LDX,N,P,XX,WORK,QRAUX,BSTAT)
            WRITE (LUNIT,570) FSTAT,BSTAT
         GO TO 540
C
C        10 X 1 MATRIX WITH LEAST SQUARES PROBLEM.
C
  380    CONTINUE
            WRITE (LUNIT,710)
            N = 10
            P = 1
            S(1) = 1.0E0
            CALL SXGEN(X,LDX,N,P,S)
            IF (NOTWRT) GO TO 390
               WRITE (LUNIT,610)
               CALL SARRAY(X,LDX,N,P,1,LUNIT)
  390       CONTINUE
            DO 400 I = 1, N
               Y1(I) = 0.0E0
               Y2(I) = 2.0E0
               IF (I .EQ. P + 1) Y2(I) = Y2(I) - FLOAT(N)
               Y1(I) = Y1(I) + X(I,1)
               Y(I) = Y1(I) + Y2(I)
               XX(I,1) = X(I,1)
  400       CONTINUE
            CALL SQRDC(X,LDX,N,P,QRAUX,JPVT,WORK,0)
            IF (NOTWRT) GO TO 410
               WRITE (LUNIT,610)
               CALL SARRAY(X,LDX,N,P,1,LUNIT)
               WRITE (LUNIT,620)
               CALL SARRAY(QRAUX,P,P,1,1,LUNIT)
  410       CONTINUE
            CALL SQRFM(X,LDX,N,P,XX,WORK,QRAUX,FSTAT)
            CALL SQRBM(X,LDX,N,P,XX,WORK,QRAUX,BSTAT)
            CALL SQRSL(X,LDX,N,P,QRAUX,Y,QY,RSD,BETA,RSD,XBETA,111,
     *                 INFO)
            RSTAT = 0.0E0
            RMAX = 0.0E0
            XBSTAT = 0.0E0
            XBMAX = 0.0E0
            DO 420 I = 1, N
               RSTAT = AMAX1(RSTAT,ABS(Y2(I)-RSD(I)))
               RMAX = AMAX1(RMAX,ABS(Y2(I)))
               XBMAX = AMAX1(XBMAX,ABS(Y2(I)))
               XBSTAT = AMAX1(XBSTAT,ABS(Y1(I)-XBETA(I)))
  420       CONTINUE
            RSTAT = RSTAT/(SMACH(1)*RMAX)
            XBSTAT = XBSTAT/(SMACH(1)*XBMAX)
            BESTAT = ABS(1.0E0-BETA(1))
            BESTAT = BESTAT/SMACH(1)
            WRITE (LUNIT,570) FSTAT,BSTAT,CONT,BESTAT,XBSTAT,RSTAT
         GO TO 540
C
C        1 X 4 MATRIX
C
  430    CONTINUE
            WRITE (LUNIT,720)
            N = 1
            P = 4
            X(1,1) = 1.0E0
            X(1,2) = 2.0E0
            X(1,3) = 0.25E0
            X(1,4) = 0.5E0
            DO 440 I = 1, 4
               JPVT(I) = 0
               XX(1,I) = X(1,I)
  440       CONTINUE
            IF (NOTWRT) GO TO 450
               WRITE (LUNIT,610)
               CALL SARRAY(X,LDX,N,P,4,LUNIT)
  450       CONTINUE
            CALL SQRDC(X,LDX,N,P,QRAUX,JPVT,WORK,1)
            IF (NOTWRT) GO TO 460
               WRITE (LUNIT,610)
               CALL SARRAY(X,LDX,N,P,10,LUNIT)
               WRITE (LUNIT,620)
               CALL SARRAY(QRAUX,N,N,1,-1,LUNIT)
  460       CONTINUE
            WRITE (LUNIT,630) (JPVT(J), J = 1, P)
            CALL SQRRX(XX,LDX,N,P,JPVT)
            CALL SQRFM(X,LDX,N,P,XX,WORK,QRAUX,FSTAT)
            CALL SQRBM(X,LDX,N,P,XX,WORK,QRAUX,BSTAT)
            WRITE (LUNIT,570) FSTAT,BSTAT
         GO TO 540
C
C        PIVOTING TEST.
C
  470    CONTINUE
            WRITE (LUNIT,660)
            N = 10
            P = 3
            S(1) = 1.0E0
            S(2) = 0.5E0
            S(3) = 0.25E0
            CALL SXGEN(X,LDX,N,P,S)
            P = 9
            DO 490 I = 1, N
               DO 480 JJJ = 1, 3
                  J = 4 - JJJ
                  JJ = 3*J
                  X(I,JJ) = X(I,J)
                  X(I,JJ-1) = X(I,JJ)/2.0E0
                  X(I,JJ-2) = X(I,JJ-1)/2.0E0
  480          CONTINUE
  490       CONTINUE
            IF (NOTWRT) GO TO 500
               WRITE (LUNIT,610)
               CALL SARRAY(X,LDX,N,P,9,LUNIT)
  500       CONTINUE
            DO 520 J = 1, P
               JPVT(J) = 0
               DO 510 I = 1, N
                  XX(I,J) = X(I,J)
  510          CONTINUE
  520       CONTINUE
            CALL SQRDC(X,LDX,N,P,QRAUX,JPVT,WORK,1)
            IF (NOTWRT) GO TO 530
               WRITE (LUNIT,610)
               CALL SARRAY(X,LDX,N,P,9,LUNIT)
               WRITE (LUNIT,620)
               CALL SARRAY(QRAUX,P,P,1,-9,LUNIT)
  530       CONTINUE
            WRITE (LUNIT,630) (JPVT(J), J = 1, P)
            CALL SQRRX(XX,LDX,N,P,JPVT)
            CALL SQRFM(X,LDX,N,P,XX,WORK,QRAUX,FSTAT)
            CALL SQRBM(X,LDX,N,P,XX,WORK,QRAUX,BSTAT)
            WRITE (LUNIT,570) FSTAT,BSTAT
  540    CONTINUE
         IF (AMAX1(FSTAT,BSTAT,BESTAT,XBSTAT,RSTAT) .GE. ERRLVL)
     *      WRITE (LUNIT,580)
  550 CONTINUE
      WRITE (LUNIT,590)
      RETURN
C
  560 FORMAT ( // 5H1CASE, I3)
  570 FORMAT ( / 11H STATISTICS //
     *         35H    FORWARD MULTIPLICATION ........, E10.2 /
     *         35H    BACK MULTIPLICATION ..........., E10.2, A1 /
     *         35H    BETA .........................., E10.2 /
     *         35H    X*BETA ........................, E10.2 /
     *         35H    RESIDUAL ......................, E10.2)
  580 FORMAT ( / 34H *****STATISTICS ABOVE ERROR LEVEL)
  590 FORMAT ( / 15H1END OF QR TEST)
  600 FORMAT (22H1LINPACK TESTER, SQR** /
     *        29H THIS VERSION DATED 08/14/78.)
  610 FORMAT ( / 2H X)
  620 FORMAT ( / 6H QRAUX)
  630 FORMAT ( / 5H JPVT // (1H , 10I5))
  640 FORMAT ( / 39H WELL CONDITIONED LEAST SQUARES PROBLEM /
     *         20H AND UNDERFLOW TEST.)
  650 FORMAT ( / 14H 4 X 10 MATRIX)
  660 FORMAT ( / 14H PIVOTING TEST /
     *         42H ON RETURN THE FIRST THREE ENTRIES OF JPVT /
     *         36H SHOULD BE 3,6,9 BUT NOT NECESSARILY /
     *         15H IN THAT ORDER.)
  670 FORMAT ( / 27H PIVOTING AND OVERFLOW TEST /
     *         26H WITH COLUMNS 1,4,7 FROZEN /
     *         42H ON RETURN THE LAST  THREE ENTRIES OF JPVT /
     *         31H SHOULD BE 1,4,7 IN THAT ORDER.)
  680 FORMAT ( / 15H 25 X 25 MATRIX)
  690 FORMAT ( / 21H MONOELEMENTAL MATRIX)
  700 FORMAT ( / 12H ZERO MATRIX)
  710 FORMAT ( / 41H 10 X 1 MATRIX WITH LEAST SQUARES PROBLEM)
  720 FORMAT ( / 13H 1 X 4 MATRIX)
      END
      SUBROUTINE SARRAY(A,LDA,M,N,NNL,LUNIT)
      INTEGER LDA,M,N,NNL,LUNIT
      REAL A(LDA,1)
C
C     FORTRAN IABS,MIN0
C
      INTEGER I,J,K,KU,NL
      NL = IABS(NNL)
      IF (NNL .LT. 0) GO TO 30
         DO 20 I = 1, M
            WRITE (LUNIT,70)
            DO 10 K = 1, N, NL
               KU = MIN0(K+NL-1,N)
               WRITE (LUNIT,70) (A(I,J), J = K, KU)
   10       CONTINUE
   20    CONTINUE
      GO TO 60
   30 CONTINUE
         DO 50 J = 1, N
            WRITE (LUNIT,70)
            DO 40 K = 1, M, NL
               KU = MIN0(K+NL-1,M)
               WRITE (LUNIT,70) (A(I,J), I = K, KU)
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
      RETURN
   70 FORMAT (1H , 8E13.6)
      END
      SUBROUTINE SQRBM(X,LDX,N,P,XX,QR,QRAUX,BSTAT)
      INTEGER LDX,N,P
      REAL BSTAT
      REAL X(LDX,1),XX(LDX,1),QR(1),QRAUX(1)
C
C     SQRBM TAKES THE OUTPUT OF SQRDC AND MULTIPLIES THE FACTORS
C     Q AND R TOGETHER.  THE RESULTS ARE COMPARED WITH XX,
C     WHICH PRESUMABLY CONTAINS THE ORIGINAL MATRIX.
C
C     LINPACK SQRSL
C     EXTERNAL SMACH
C     FORTRAN AMAX1,ABS,MIN0
C
      INTEGER IU,JP1,I,INFO,J
      REAL SMACH,EMAX,XMAX
      EMAX = 0.0E0
      XMAX = 0.0E0
      DO 50 J = 1, P
         IU = MIN0(N,J)
         DO 10 I = 1, IU
            QR(I) = X(I,J)
   10    CONTINUE
         JP1 = J + 1
         IF (N .LT. JP1) GO TO 30
         DO 20 I = JP1, N
            QR(I) = 0.0E0
   20    CONTINUE
   30    CONTINUE
         CALL SQRSL(X,LDX,N,P,QRAUX,QR,QR,QR,QR,QR,QR,10000,INFO)
         DO 40 I = 1, N
            EMAX = AMAX1(EMAX,ABS(XX(I,J)-QR(I)))
            XMAX = AMAX1(XMAX,ABS(XX(I,J)))
   40    CONTINUE
   50 CONTINUE
      IF (XMAX .NE. 0.0E0) GO TO 80
         IF (EMAX .NE. 0.0E0) GO TO 60
            BSTAT = 0.0E0
         GO TO 70
   60    CONTINUE
            BSTAT = 1.0E20
   70    CONTINUE
      GO TO 90
   80 CONTINUE
         BSTAT = (EMAX/XMAX)/SMACH(1)
   90 CONTINUE
      RETURN
      END
      SUBROUTINE SQRFM(X,LDX,N,P,XX,QTX,QRAUX,FSTAT)
      INTEGER LDX,N,P
      REAL FSTAT
      REAL X(LDX,1),XX(LDX,1),QTX(1),QRAUX(1)
C
C     SQRFM TAKES THE OUTPUT OF SQRDC AND APPLIES THE
C     TRANSFORMATIONS TO XX, WHICH PRESUMABLY CONTAINS THE
C     ORIGINAL MATRIX.  THE RESULTS ARE COMPARED WITH THE
C     R PART OF THE QR DECOMPOSITION.
C
C     LINPACK SQRSL
C     EXTERNAL SMACH
C     FORTRAN AMAX1,ABS,MIN0
C
      INTEGER IU,JP1,I,INFO,J
      REAL SMACH,EMAX,RMAX
      EMAX = 0.0E0
      RMAX = 0.0E0
      DO 50 J = 1, P
         DO 10 I = 1, N
            QTX(I) = XX(I,J)
   10    CONTINUE
         CALL SQRSL(X,LDX,N,P,QRAUX,QTX,QTX,QTX,QTX,QTX,QTX,1000,INFO)
         IU = MIN0(J,N)
         DO 20 I = 1, IU
            EMAX = AMAX1(EMAX,ABS(X(I,J)-QTX(I)))
            RMAX = AMAX1(RMAX,ABS(X(I,J)))
   20    CONTINUE
         JP1 = J + 1
         IF (N .LT. JP1) GO TO 40
         DO 30 I = JP1, N
            EMAX = AMAX1(EMAX,ABS(QTX(I)))
   30    CONTINUE
   40    CONTINUE
   50 CONTINUE
      IF (RMAX .NE. 0.0E0) GO TO 80
         IF (EMAX .NE. 0.0E0) GO TO 60
            FSTAT = 0.0E0
         GO TO 70
   60    CONTINUE
            FSTAT = 1.0E20
   70    CONTINUE
      GO TO 90
   80 CONTINUE
         FSTAT = (EMAX/RMAX)/SMACH(1)
   90 CONTINUE
      RETURN
      END
      SUBROUTINE SQRRX(XX,LDX,N,P,JP)
      INTEGER N,P,LDX
      INTEGER JP(1)
      REAL XX(LDX,1)
C
C     SQRRX REARRANGES THE COLUMNS OF XX IN THE ORDER SPECIFIED
C     BY JPVT.
C
C     EXTERNAL SSWAP
C
      INTEGER I,IP,I1,J,JPVT(50),K
      DO 10 J = 1, P
         JPVT(J) = JP(J)
   10 CONTINUE
      IP = P - 1
      DO 60 I = 1, IP
         IF (JPVT(I) .EQ. I) GO TO 50
            I1 = I + 1
            DO 20 J = I1, P
C           ...EXIT
               IF (JPVT(I) .EQ. J) GO TO 30
   20       CONTINUE
   30       CONTINUE
            CALL SSWAP(N,XX(1,I),1,XX(1,J),1)
            DO 40 K = I1, P
               IF (JPVT(K) .EQ. I) JPVT(K) = JPVT(I)
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
      RETURN
      END
      SUBROUTINE SXGEN(X,LDX,N,P,S)
      INTEGER LDX,N,P
      REAL X(LDX,1),S(1)
C
C     SXGEN GENERATES AN N X P MATRIX X HAVING SINGULAR
C     VALUES S(1), S(2), ..., S(M), WHERE M = MIN(N,P)
C
C     FORTRAN FLOAT,MIN0
C
      INTEGER I,J,M,MP1
      REAL FN,FP
      REAL FAC,RU,T
      FP = FLOAT(P)
      FN = FLOAT(N)
      M = MIN0(N,P)
      RU = 1.0E0
      FAC = 1.0E0/FP
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
               X(I,J) = 0.0E0
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      DO 80 J = 1, P
         T = 0.0E0
         DO 60 I = 1, N
            T = T + X(I,J)
   60    CONTINUE
         DO 70 I = 1, N
            X(I,J) = X(I,J) - 2.0E0*T/FN
   70    CONTINUE
   80 CONTINUE
      RETURN
      END
