C     MAIN PROGRAM
      INTEGER LUNIT
C     ALLOW 5000 UNDERFLOWS.
      CALL TRAPS(0,0,5001,0,0)
C
C     OUTPUT UNIT NUMBER
C
      LUNIT = 6
C
      CALL SGTTS(LUNIT)
      STOP
      END
      SUBROUTINE SGTTS(LUNIT)
C     LUNIT IS THE OUTPUT UNIT NUMBER
C
C     TESTS
C        SGTSL,SPTSL
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JACK DONGARRA, ARGONNE NATIONAL LABORATORY.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK SGTSL,CTPSL
C     EXTERNAL SGTXX
C     BLAS SASUM
C     FORTRAN ABS,AIMAG,AMAX1,FLOAT,REAL
C
C     INTERNAL VARIABLES
C
      INTEGER LUNIT
      REAL B(20),BSAVE(20),D(20),EYE,C(20),E(20),X(20)
      REAL ANORM,EN,ENORM,EPS,Q(2),RNORM,SASUM,XNORM
      REAL SMACH
      INTEGER I,INFO,IPT,PD,KASE,KFAIL(2),KSING,N,NM1,NPRINT,POSDEF
      EYE = 0.0E0
C     REAL EYE = IMAGINARY UNIT, REAL EYE = ZERO
C
C
C     WRITE MATRIX AND SOLUTIONS IF  N .LE. NPRINT
C
      NPRINT = 3
C
      WRITE (LUNIT,230)
      DO 10 I = 1, 2
         KFAIL(I) = 0
   10 CONTINUE
      KSING = 0
C
C     COMPUTE MACHINE EPSILON
C
      EPS = SMACH(1)
      WRITE (LUNIT,240) EPS
      WRITE (LUNIT,220)
C
C     START MAIN LOOP
C
      KASE = 1
   20 CONTINUE
C
C        GENERATE TEST MATRIX
C
         CALL SGTXX(C,D,E,N,KASE,POSDEF)
C
C        N = 0 SIGNALS NO MORE TEST MATRICES
C
C     ...EXIT
         IF (N .LE. 0) GO TO 210
         INFO = 0
         PD = 1
         IF (POSDEF .EQ. 1) PD = 2
         DO 200 IPT = 1, PD
            WRITE (LUNIT,250) KASE
            WRITE (LUNIT,260) N
            IF (N .GT. 1) GO TO 30
               ANORM = ABS(D(1))
               WRITE (LUNIT,450) D(1)
               X(1) = 1.0E0
               B(1) = D(1)
               BSAVE(1) = B(1)
            GO TO 110
   30       CONTINUE
               NM1 = N - 1
               ANORM = ABS(D(1)) + ABS(C(2))
               IF (N .LE. 2) GO TO 50
                  DO 40 I = 2, NM1
                     ANORM = AMAX1(ANORM,
     *                             ABS(C(I+1))+ABS(D(I))+ABS(E(I-1)))
   40             CONTINUE
   50          CONTINUE
               ANORM = AMAX1(ANORM,ABS(E(N-1))+ABS(D(N)))
               WRITE (LUNIT,430) ANORM
C
               IF (N .GT. NPRINT) GO TO 60
                  WRITE (LUNIT,220)
                  WRITE (LUNIT,450) (C(I), I = 2, N)
                  WRITE (LUNIT,220)
                  WRITE (LUNIT,450) (D(I), I = 1, N)
                  WRITE (LUNIT,220)
                  WRITE (LUNIT,450) (E(I), I = 1, NM1)
                  WRITE (LUNIT,220)
   60          CONTINUE
C
C              GENERATE EXACT SOLUTION
C
               X(1) = 1.0E0
               IF (N .GE. 2) X(2) = EYE
               IF (N .LE. 2) GO TO 80
                  DO 70 I = 3, N
                     X(I) = -X(I-2)
   70             CONTINUE
   80          CONTINUE
C
C              SAVE MATRIX AND GENERATE R.H.S.
C
               B(1) = D(1)*X(1) + E(1)*X(2)
               BSAVE(1) = B(1)
               IF (N .LE. 2) GO TO 100
                  DO 90 I = 2, NM1
                     B(I) = C(I)*X(I-1) + D(I)*X(I) + E(I)*X(I+1)
                     BSAVE(I) = B(I)
   90             CONTINUE
  100          CONTINUE
               B(N) = C(N)*X(N-1) + D(N)*X(N)
               BSAVE(N) = B(N)
  110       CONTINUE
C
C           FACTOR AND SOLVE A GENERAL TRIDIAGONAL SYSTEM
C
            IF (IPT .EQ. 1) CALL SGTSL(N,C,D,E,B,INFO)
C
C           TEST FOR SINGULARITY
C
            IF (INFO .EQ. 0) GO TO 120
               WRITE (LUNIT,270)
            GO TO 190
  120       CONTINUE
C
C              FACTOR AND SOLVE A POSITIVE DEFINITE SYSTEM
C
               IF (IPT .EQ. 2) CALL SPTSL(N,D,E,B)
               IF (IPT .EQ. 1) WRITE (LUNIT,280)
               IF (IPT .EQ. 2) WRITE (LUNIT,290)
               IF (N .GT. NPRINT) GO TO 130
                  WRITE (LUNIT,300)
                  WRITE (LUNIT,460) (B(I), I = 1, N)
                  WRITE (LUNIT,310)
                  WRITE (LUNIT,460) (BSAVE(I), I = 1, N)
  130          CONTINUE
C
C              COMPUTE ERRORS AND RESIDUALS
C                 E  =  X - X
C                 R  =  B - A*X
C
               XNORM = SASUM(N,X,1)
               CALL SGTXX(C,D,E,N,KASE,POSDEF)
               IF (N .GT. 1) GO TO 140
                  RNORM = ABS(D(1)*B(1)-BSAVE(1))
                  ENORM = ABS(B(1)-X(1))
               GO TO 170
  140          CONTINUE
                  ENORM = ABS(B(1)-X(1))
                  RNORM = ABS(D(1)*B(1)+E(1)*B(2)-BSAVE(1))
                  IF (N .LE. 2) GO TO 160
                     DO 150 I = 2, NM1
                        RNORM = RNORM
     *                          + ABS(C(I)*B(I-1)+D(I)*B(I)+E(I)*B(I+1)
     *                                -BSAVE(I))
                        ENORM = ENORM + ABS(B(I)-X(I))
  150                CONTINUE
  160             CONTINUE
                  RNORM = RNORM + ABS(C(N)*B(N-1)+D(N)*B(N)-BSAVE(N))
                  ENORM = ENORM + ABS(B(N)-X(N))
  170          CONTINUE
C
               WRITE (LUNIT,320) ENORM
               WRITE (LUNIT,330) RNORM
C
C              COMPUTE TEST RATIOS
C
               EN = FLOAT(N)
               Q(1) = RNORM/(EPS*ANORM*XNORM)
               Q(2) = ENORM/(EPS*XNORM)
               WRITE (LUNIT,220)
               WRITE (LUNIT,340)
               WRITE (LUNIT,220)
               WRITE (LUNIT,400)
               WRITE (LUNIT,410)
               WRITE (LUNIT,420)
               WRITE (LUNIT,220)
               WRITE (LUNIT,440) (Q(I), I = 1, 2)
               WRITE (LUNIT,220)
               IF (N .EQ. 1) EN = 2.0E0
               DO 180 I = 1, 2
                  IF (Q(I) .GT. EN) KFAIL(I) = KFAIL(I) + 1
  180          CONTINUE
  190       CONTINUE
C
            WRITE (LUNIT,350)
  200    CONTINUE
         KASE = KASE + 1
      GO TO 20
  210 CONTINUE
C
C     FINISH MAIN LOOP
C
C     SUMMARY
C
      WRITE (LUNIT,360)
      KASE = KASE - 1
      WRITE (LUNIT,370) KASE
      WRITE (LUNIT,380) KSING
      WRITE (LUNIT,390) KFAIL
      WRITE (LUNIT,470)
      RETURN
C
C     ALL FORMATS
C
  220 FORMAT (1H )
  230 FORMAT (29H1LINPACK TESTER, SGT**, SPT** /
     *        30H THIS VERSION DATED 08/14/78 .)
  240 FORMAT (18H MACHINE EPSILON =, 1PE13.5)
  250 FORMAT (14H MATRIX NUMBER, I4)
  260 FORMAT (4H N =, I4)
  270 FORMAT (16H MAYBE SINGULAR.)
  280 FORMAT (18H RESULTS FOR SGTSL)
  290 FORMAT (18H RESULTS FOR SPTSL)
  300 FORMAT ( / 4H X =)
  310 FORMAT ( / 4H B =)
  320 FORMAT (14H ERROR NORMS =, 1P2E13.5)
  330 FORMAT (14H RESID NORMS =, 1P2E13.5)
  340 FORMAT (26H TEST RATIOS.. E = MACHEPS)
  350 FORMAT ( / 14H ************* /)
  360 FORMAT (8H1SUMMARY)
  370 FORMAT (18H NUMBER OF TESTS =, I4)
  380 FORMAT (30H NUMBER OF SINGULAR MATRICES =, I4)
  390 FORMAT (30H NUMBER OF SUSPICIOUS RATIOS =, 8I4)
  400 FORMAT (20H    ERROR     RESID )
  410 FORMAT (2(10H   -------))
  420 FORMAT (20H     E*X      E*A*X )
  430 FORMAT (14H NORM(A)     =, 1PE13.5)
  440 FORMAT (8F10.4)
  450 FORMAT (6G11.4)
  460 FORMAT (2G14.6)
  470 FORMAT ( / 12H END OF TEST)
      END
      SUBROUTINE SGTXX(C,D,E,N,KASE,POSDEF)
C     FORTRAN FLOAT
C
      INTEGER N,KASE,POSDEF
      REAL C(1),D(1),E(1)
C
      REAL EYE
      INTEGER I
C
      EYE = 0.0E0
      GO TO (10,20,30,30,30,50,50,70,70,90,110), KASE
C
   10 CONTINUE
         N = 1
         D(1) = 1.0E0
         POSDEF = 1
      GO TO 120
C
   20 CONTINUE
         N = 2
         D(1) = 4.0E0
         D(2) = 4.0E0
         C(2) = 2.0E0
         E(1) = 2.0E0
         POSDEF = 1
      GO TO 120
C
   30 CONTINUE
         N = (KASE - 2)*3
         DO 40 I = 1, N
            C(I) = 1.0E0/(FLOAT(2*I+2) + EYE)
            D(I) = 1.0E0/(FLOAT(2*I+1) + EYE)
            E(I) = 1.0E0/(FLOAT(2*I) + EYE)
   40    CONTINUE
         POSDEF = 0
      GO TO 120
C
   50 CONTINUE
         IF (KASE .EQ. 6) N = 10
         IF (KASE .EQ. 7) N = 20
         DO 60 I = 1, N
            C(I) = 1.0E0
            D(I) = 4.0E0
            E(I) = 1.0E0
   60    CONTINUE
         POSDEF = 1
      GO TO 120
C
   70 CONTINUE
         IF (KASE .EQ. 8) N = 10
         IF (KASE .EQ. 9) N = 20
         DO 80 I = 1, N
            C(I) = 1.0E0 + EYE
            D(I) = 4.0E0 + EYE
            E(I) = 1.0E0 - EYE
   80    CONTINUE
         POSDEF = 1
      GO TO 120
C
   90 CONTINUE
         N = 10
         DO 100 I = 1, N
            C(I) = 0.0E0
            D(I) = 1.0E0
            E(I) = 0.0E0
  100    CONTINUE
         POSDEF = 1
      GO TO 120
C
  110 CONTINUE
         N = 0
  120 CONTINUE
      RETURN
      END
