C     MAIN PROGRAM
      INTEGER LUNIT
C     ALLOW 5000 UNDERFLOWS.
      CALL TRAPS(0,0,5001,0,0)
C
C     OUTPUT UNIT NUMBER
C
      LUNIT = 6
C
      CALL DTRTS(LUNIT)
      STOP
      END
      SUBROUTINE DTRTS(LUNIT)
      INTEGER LUNIT
C     LUNIT IS THE OUTPUT UNIT NUMBER
C
C     TESTS
C        DTRCO, DTRSL
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK DTRCO,DTRSL
C     EXTERNAL DTRXX,DMACH
C     BLAS DAXPY,DDOT,DSCAL,DASUM
C     FORTRAN DABS,DMAX1,DBLE,FLOAT,MAX0
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION A(15,15),B(15),BT(15),X(15),XEXACT(15),XT(15),
     *                 Z(15)
      DOUBLE PRECISION AINV(15,15),DET(2)
      DOUBLE PRECISION DDOT,STUFF,T
      DOUBLE PRECISION ANORM,AINORM,RCOND,COND,COND1,DMACH,EPS
      DOUBLE PRECISION ENORM,ETNORM,RNORM,RTNORM,XNORM,XTNORM,EN,DASUM
      DOUBLE PRECISION FNORM,ONEPX,Q(7),QS(7)
      INTEGER I,INFO,J,JOB,KASE,KOUNT,KSING,LDA,ML,MU,N,NPRINT
      INTEGER KSUSP(7),IQ(7)
C
      LDA = 15
C
C     WRITE MATRIX AND SOLUTIONS IF  N .LE. NPRINT
C
      NPRINT = 3
C
      WRITE (LUNIT,380)
      WRITE (LUNIT,730)
C
      DO 10 I = 1, 7
         KSUSP(I) = 0
   10 CONTINUE
      KSING = 0
C
C     SET EPS TO ROUNDING UNIT
C
      EPS = DMACH(1)
      WRITE (LUNIT,390) EPS
      WRITE (LUNIT,370)
C
C     START MAIN LOOP
C
      KASE = 1
   20 CONTINUE
C
C        GENERATE TEST MATRIX
C
         CALL DTRXX(A,LDA,N,KASE,LUNIT)
C
C        N = 0 SIGNALS NO MORE TEST MATRICES
C
C     ...EXIT
         IF (N .LE. 0) GO TO 360
         ANORM = 0.0D0
         DO 30 J = 1, N
            ANORM = DMAX1(ANORM,DASUM(N,A(1,J),1))
   30    CONTINUE
         WRITE (LUNIT,570) ANORM
C
         IF (N .GT. NPRINT) GO TO 50
            WRITE (LUNIT,370)
            DO 40 I = 1, N
               WRITE (LUNIT,600) (A(I,J), J = 1, N)
   40       CONTINUE
            WRITE (LUNIT,370)
   50    CONTINUE
C
C        GENERATE EXACT SOLUTION
C
         XEXACT(1) = 1.0D0
         IF (N .GE. 2) XEXACT(2) = 0.0D0
         IF (N .LE. 2) GO TO 70
            DO 60 I = 3, N
               XEXACT(I) = -XEXACT(I-2)
   60       CONTINUE
   70    CONTINUE
C
C        GENERATE R.H.S.
C
         DO 90 I = 1, N
            B(I) = 0.0D0
            BT(I) = 0.0D0
            DO 80 J = 1, N
               B(I) = B(I) + A(I,J)*XEXACT(J)
               BT(I) = BT(I) + A(J,I)*XEXACT(J)
   80       CONTINUE
            X(I) = B(I)
            XT(I) = BT(I)
   90    CONTINUE
C
C        UPPER OR LOWER TRIANGULAR
C
         ML = 0
         MU = 0
         DO 120 J = 1, N
            DO 110 I = 1, N
               IF (A(I,J) .EQ. 0.0D0) GO TO 100
                  IF (I .LT. J) MU = MAX0(MU,J-I)
                  IF (I .GT. J) ML = MAX0(ML,I-J)
  100          CONTINUE
  110       CONTINUE
  120    CONTINUE
         WRITE (LUNIT,670) ML,MU
         IF (ML .NE. 0 .AND. MU .NE. 0) GO TO 350
            IF (MU .EQ. 0) JOB = 0
            IF (ML .EQ. 0) JOB = 1
            IF (JOB .EQ. 0) WRITE (LUNIT,710)
            IF (JOB .EQ. 1) WRITE (LUNIT,720)
            STUFF = 4095.0D0
            DO 140 J = 1, N
               DO 130 I = 1, N
                  IF (I .LT. J .AND. JOB .EQ. 0) A(I,J) = STUFF
                  IF (I .GT. J .AND. JOB .EQ. 1) A(I,J) = STUFF
  130          CONTINUE
  140       CONTINUE
C
C           ESTIMATE CONDITION
C
            CALL DTRCO(A,LDA,N,RCOND,Z,JOB)
C
C           OUTPUT NULL VECTOR IF N .LE. NPRINT
C
            IF (N .GT. NPRINT) GO TO 160
               WRITE (LUNIT,620)
               DO 150 I = 1, N
                  WRITE (LUNIT,630) Z(I)
  150          CONTINUE
               WRITE (LUNIT,370)
  160       CONTINUE
C
C
C           TEST FOR SINGULARITY
C
            IF (RCOND .GT. 0.0D0) GO TO 170
               WRITE (LUNIT,610) RCOND
               WRITE (LUNIT,400)
               KSING = KSING + 1
            GO TO 340
  170       CONTINUE
               COND = 1.0D0/RCOND
               WRITE (LUNIT,420) COND
               ONEPX = 1.0D0 + RCOND
               IF (ONEPX .EQ. 1.0D0) WRITE (LUNIT,410)
C
C              COMPUTE INVERSE, DETERMINANT AND COND1 = TRUE CONDITION
C
               DO 190 J = 1, N
                  DO 180 I = 1, N
                     AINV(I,J) = A(I,J)
  180             CONTINUE
  190          CONTINUE
               CALL DTRDI(AINV,LDA,N,DET,110+JOB,INFO)
               AINORM = 0.0D0
               DO 200 J = 1, N
                  IF (JOB .EQ. 0)
     *               AINORM = DMAX1(AINORM,DASUM(N-J+1,AINV(J,J),1))
                  IF (JOB .EQ. 1)
     *               AINORM = DMAX1(AINORM,DASUM(J,AINV(1,J),1))
  200          CONTINUE
               COND1 = ANORM*AINORM
               WRITE (LUNIT,430) COND1
               WRITE (LUNIT,650) DET(1)
               WRITE (LUNIT,660) DET(2)
C
C              SOLVE  A*X = B  AND  TRANS(A)*XT = BT
C
               CALL DTRSL(A,LDA,N,X,JOB,INFO)
               CALL DTRSL(A,LDA,N,XT,JOB+10,INFO)
C
               IF (N .GT. NPRINT) GO TO 230
                  WRITE (LUNIT,440)
                  DO 210 I = 1, N
                     WRITE (LUNIT,640) X(I)
  210             CONTINUE
                  WRITE (LUNIT,450)
                  DO 220 I = 1, N
                     WRITE (LUNIT,640) XT(I)
  220             CONTINUE
                  WRITE (LUNIT,370)
  230          CONTINUE
C
C              RESTORE ZEROS IN OTHER TRIANGLE
C
               DO 260 J = 1, N
                  DO 250 I = 1, N
                     IF (A(I,J) .NE. STUFF) GO TO 240
                        A(I,J) = 0.0D0
                        AINV(I,J) = 0.0D0
  240                CONTINUE
  250             CONTINUE
  260          CONTINUE
C
C              COMPUTE ERRORS AND RESIDUALS
C                 E  =  X - XEXACT
C                 ET =  XT - XEXACT
C                 R  =  B - A*X
C                 RT =  BT - A*XT
C                 AI = A*INV(A) - I
C
               XNORM = DASUM(N,X,1)
               XTNORM = DASUM(N,XT,1)
               ENORM = 0.0D0
               ETNORM = 0.0D0
               FNORM = 0.0D0
               DO 270 J = 1, N
                  ENORM = ENORM + DABS(X(J)-XEXACT(J))
                  ETNORM = ETNORM + DABS(XT(J)-XEXACT(J))
                  T = -X(J)
                  CALL DAXPY(N,T,A(1,J),1,B,1)
                  BT(J) = BT(J) - DDOT(N,A(1,J),1,XT,1)
  270          CONTINUE
               RNORM = DASUM(N,B,1)
               RTNORM = DASUM(N,BT,1)
C
C
C              A*INV(A) - I
C
               AINORM = 0.0D0
               DO 300 J = 1, N
                  DO 280 I = 1, N
                     B(I) = 0.0D0
  280             CONTINUE
                  DO 290 K = 1, N
                     T = AINV(K,J)
                     CALL DAXPY(N,T,A(1,K),1,B,1)
  290             CONTINUE
                  B(J) = B(J) - 1.0D0
                  AINORM = DMAX1(AINORM,DASUM(N,B,1))
  300          CONTINUE
C
               WRITE (LUNIT,460) ENORM,ETNORM
               WRITE (LUNIT,470) RNORM,RTNORM
               WRITE (LUNIT,580) AINORM
C
C              COMPUTE TEST RATIOS
C
               Q(1) = COND/COND1
               Q(2) = COND1/COND
               Q(3) = ENORM/(EPS*COND*XNORM)
               Q(4) = ETNORM/(EPS*COND*XTNORM)
               Q(5) = RNORM/(EPS*ANORM*XNORM)
               Q(6) = RTNORM/(EPS*ANORM*XTNORM)
               Q(7) = AINORM/(EPS*COND)
               WRITE (LUNIT,370)
               WRITE (LUNIT,480)
               WRITE (LUNIT,370)
               WRITE (LUNIT,540)
               WRITE (LUNIT,550)
               WRITE (LUNIT,560)
               WRITE (LUNIT,370)
               WRITE (LUNIT,590) (Q(I), I = 1, 7)
               WRITE (LUNIT,370)
C
C              LOOK FOR SUSPICIOUS RATIOS
C
               QS(1) = 1.0D0 + 4.0D0*EPS
               QS(2) = 10.0D0
               EN = DBLE(FLOAT(N))
               IF (N .EQ. 1) EN = 2.0D0
               DO 310 I = 3, 7
                  QS(I) = EN
  310          CONTINUE
               KOUNT = 0
               DO 330 I = 1, 7
                  IQ(I) = 0
                  IF (Q(I) .LE. QS(I)) GO TO 320
                     IQ(I) = 1
                     KSUSP(I) = KSUSP(I) + 1
                     KOUNT = KOUNT + 1
  320             CONTINUE
  330          CONTINUE
               IF (KOUNT .EQ. 0) WRITE (LUNIT,690)
               IF (KOUNT .NE. 0) WRITE (LUNIT,700) (IQ(I), I = 1, 7)
               WRITE (LUNIT,370)
  340       CONTINUE
  350    CONTINUE
C
         WRITE (LUNIT,490)
         KASE = KASE + 1
      GO TO 20
  360 CONTINUE
C
C     FINISH MAIN LOOP
C
C     SUMMARY
C
      WRITE (LUNIT,500)
      KASE = KASE - 1
      WRITE (LUNIT,510) KASE
      WRITE (LUNIT,520) KSING
      WRITE (LUNIT,530) KSUSP
      WRITE (LUNIT,680)
      RETURN
C
C     MOST FORMATS, ALSO SOME IN DTRXX
C
  370 FORMAT (1H )
  380 FORMAT (22H1LINPACK TESTER, DTR**)
  390 FORMAT ( / 18H MACHINE EPSILON =, 1PD13.5)
  400 FORMAT ( / 19H EXACT SINGULARITY. /)
  410 FORMAT ( / 16H MAYBE SINGULAR. /)
  420 FORMAT (14H COND        =, 1PD13.5)
  430 FORMAT (14H ACTUAL COND =, 1PD13.5)
  440 FORMAT ( / 4H X =)
  450 FORMAT ( / 5H XT =)
  460 FORMAT (14H ERROR NORMS =, 1P2D13.5)
  470 FORMAT (14H RESID NORMS =, 1P2D13.5)
  480 FORMAT (26H TEST RATIOS.. E = MACHEPS)
  490 FORMAT ( / 14H ************* /)
  500 FORMAT (8H1SUMMARY)
  510 FORMAT (18H NUMBER OF TESTS =, I4)
  520 FORMAT (30H NUMBER OF SINGULAR MATRICES =, I4)
  530 FORMAT (30H NUMBER OF SUSPICIOUS RATIOS =, 7I4)
  540 FORMAT (30H     COND     ACTUAL    ERROR ,
     *        40H   ERROR-T    RESID    RESID-T   A*AI-I )
  550 FORMAT (7(10H   -------))
  560 FORMAT (30H    ACTUAL     COND   E*COND*X,
     *        40H  E*COND*X    E*A*X     E*A*X    E*COND )
  570 FORMAT (14H NORM(A)     =, 1PD13.5)
  580 FORMAT (14H NORM(A*AI-I)=, 1PD13.5)
  590 FORMAT (7F10.4)
  600 FORMAT (1H , 6G11.4)
  610 FORMAT (14H 1/COND      =, 1PD13.5)
  620 FORMAT ( / 7H NULL =)
  630 FORMAT (2G14.6)
  640 FORMAT (2G14.6)
  650 FORMAT (14H DET FRACT   =, 2F9.5)
  660 FORMAT (14H DET EXPON   =, 2F9.0)
  670 FORMAT (5H ML =, I2, 6H  MU =, I2)
  680 FORMAT ( / 12H END OF TEST)
  690 FORMAT (21H NO SUSPICIOUS RATIOS)
  700 FORMAT (I8, 5I10 / 7X, 28H1 INDICATES SUSPICIOUS RATIO)
  710 FORMAT (26H LOWER TRIANGULAR, JOB = 0)
  720 FORMAT (26H UPPER TRIANGULAR, JOB = 1)
  730 FORMAT (29H THIS VERSION DATED 08/14/78.)
      END
      SUBROUTINE DTRXX(A,LDA,N,KASE,LUNIT)
      INTEGER LDA,N,KASE,LUNIT
      DOUBLE PRECISION A(LDA,1)
C
C     GENERATES DOUBLE PRECISION TRIANGULAR TEST MATRICES
C
C     EXTERNAL DMACH
C     FORTRAN DBLE,FLOAT
      DOUBLE PRECISION T1,T2
      DOUBLE PRECISION DMACH,HUGE,TINY
      INTEGER I,J
C
      GO TO (10,10,10,50,50,70,70,70,110,150,200,240,280,320,360,410,
     *       460), KASE
C
C     KASE 1, 2 AND 3
C
   10 CONTINUE
         N = 3*KASE
         WRITE (LUNIT,20) KASE,N
   20    FORMAT (5H KASE, I3, 3X, 16HHILBERT-HALF     / 4H N =, I4)
         DO 40 J = 1, N
            DO 30 I = 1, N
               A(I,J) = 0.0D0
               IF (I .GE. J - 3 .AND. I .LE. J)
     *            A(I,J) = 1.0D0/DBLE(FLOAT(I+J-1))
   30       CONTINUE
   40    CONTINUE
      GO TO 470
C
C     KASE 4 AND 5
C
   50 CONTINUE
         N = 1
         WRITE (LUNIT,60) KASE,N
   60    FORMAT (5H KASE, I3, 3X, 16HMONOELEMENTAL    / 4H N =, I4)
         IF (KASE .EQ. 4) A(1,1) = 3.0D0
         IF (KASE .EQ. 5) A(1,1) = 0.0D0
      GO TO 470
C
C     KASE 6, 7 AND 8
C
   70 CONTINUE
         N = 15
         WRITE (LUNIT,80) KASE,N
   80    FORMAT (5H KASE, I3, 3X, 16HBIDIAGONAL       / 4H N =, I4)
         T1 = 0.0D0
         T2 = 0.0D0
         IF (KASE .EQ. 7) T1 = 100.0D0
         IF (KASE .EQ. 8) T2 = 100.0D0
         DO 100 I = 1, N
            DO 90 J = 1, N
               A(I,J) = 0.0D0
               IF (I .EQ. J) A(I,I) = 4.0D0
               IF (I .EQ. J - 1) A(I,J) = T1
               IF (I .EQ. J + 1) A(I,J) = T2
   90       CONTINUE
  100    CONTINUE
      GO TO 470
C
C     KASE 9
C
  110 CONTINUE
         N = 5
         WRITE (LUNIT,120) KASE,N
  120    FORMAT (5H KASE, I3, 3X, 16HHALF OF RANK ONE / 4H N =, I4)
         DO 140 I = 1, N
            DO 130 J = 1, N
               A(I,J) = 0.0D0
               IF (I .GE. J) A(I,J) = 10.0D0**(I - J)
  130       CONTINUE
  140    CONTINUE
      GO TO 470
C
C     KASE 10
C
  150 CONTINUE
         N = 4
         WRITE (LUNIT,160) KASE,N
  160    FORMAT (5H KASE, I3, 3X, 16HZERO COLUMN      / 4H N =, I4)
         DO 190 I = 1, N
            DO 180 J = 1, N
               A(I,J) = 0.0D0
               IF (I .LT. J) GO TO 170
                  T1 = DBLE(FLOAT(J-3))
                  T2 = DBLE(FLOAT(I))
                  A(I,J) = T1/T2
  170          CONTINUE
  180       CONTINUE
  190    CONTINUE
      GO TO 470
C
C     KASE 11
C
  200 CONTINUE
         N = 5
         WRITE (LUNIT,210) KASE,N
  210    FORMAT (5H KASE, I3, 3X, 16HTEST COND        / 4H N =, I4)
         DO 230 I = 1, N
            DO 220 J = 1, N
               IF (I .EQ. J) A(I,J) = 1.0D0
               IF (I .GT. J) A(I,J) = 0.0D0
               IF (I .LT. J) A(I,J) = -1.0D0
  220       CONTINUE
  230    CONTINUE
      GO TO 470
C
C     KASE 12
C
  240 CONTINUE
         N = 3
         WRITE (LUNIT,250) KASE,N
  250    FORMAT (5H KASE, I3, 3X, 16HIDENTITY         / 4H N =, I4)
         DO 270 I = 1, N
            DO 260 J = 1, N
               IF (I .EQ. J) A(I,I) = 1.0D0
               IF (I .NE. J) A(I,J) = 0.0D0
  260       CONTINUE
  270    CONTINUE
      GO TO 470
C
C     KASE 13
C
  280 CONTINUE
         N = 6
         WRITE (LUNIT,290) KASE,N
  290    FORMAT (5H KASE, I3, 3X, 16HUPPER TRIANGULAR / 4H N =, I4)
         DO 310 I = 1, N
            DO 300 J = 1, N
               IF (I .GT. J) A(I,J) = 0.0D0
               IF (I .LE. J) A(I,J) = DBLE(FLOAT(I-J-1))*1.0D0
  300       CONTINUE
  310    CONTINUE
      GO TO 470
C
C     KASE 14
C
  320 CONTINUE
         N = 6
         WRITE (LUNIT,330) KASE,N
  330    FORMAT (5H KASE, I3, 3X, 16HLOWER TRIANGULAR / 4H N =, I4)
         DO 350 I = 1, N
            DO 340 J = 1, N
               IF (I .LT. J) A(I,J) = 0.0D0
               IF (I .GE. J) A(I,J) = DBLE(FLOAT(I-J+1))*1.0D0
  340       CONTINUE
  350    CONTINUE
      GO TO 470
C
C     KASE 15
C
  360 CONTINUE
         N = 5
         WRITE (LUNIT,370) KASE,N
  370    FORMAT (5H KASE, I3, 3X, 16HNEAR UNDERFLOW   / 4H N =, I4)
         TINY = DMACH(2)
         WRITE (LUNIT,380) TINY
  380    FORMAT (14H TINY        =, 1PD13.5)
         DO 400 I = 1, N
            DO 390 J = 1, N
               A(I,J) = 0.0D0
               IF (I .LE. J)
     *            A(I,J) = TINY*(DBLE(FLOAT(J))/DBLE(FLOAT(I)))
  390       CONTINUE
  400    CONTINUE
      GO TO 470
C
C     KASE 16
C
  410 CONTINUE
         N = 5
         WRITE (LUNIT,420) KASE,N
  420    FORMAT (5H KASE, I3, 3X, 16HNEAR OVERFLOW    / 4H N =, I4)
         HUGE = DMACH(3)
         WRITE (LUNIT,430) HUGE
  430    FORMAT (14H HUGE        =, 1PD16.5)
         DO 450 I = 1, N
            DO 440 J = 1, N
               A(I,J) = 0.0D0
               IF (I .GE. J)
     *            A(I,J) = HUGE*(DBLE(FLOAT(J))/DBLE(FLOAT(I)))
  440       CONTINUE
  450    CONTINUE
      GO TO 470
C
  460 CONTINUE
         N = 0
  470 CONTINUE
      RETURN
C
      END
