C     MAIN PROGRAM
      INTEGER LUNIT
C     ALLOW 5000 UNDERFLOWS.
      CALL TRAPS(0,0,5001,0,0)
C
C     OUTPUT UNIT NUMBER
C
      LUNIT = 6
C
      CALL SGETS(LUNIT)
      STOP
      END
      SUBROUTINE SGETS(LUNIT)
C     LUNIT IS THE OUTPUT UNIT NUMBER
C
C     TESTS
C        SGECO,SGEFA,SGESL,SGEDI,SGBCO,SGBFA,SGBSL,SGBDI
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK SGECO,SGESL,SGEDI,SGBCO,SGBSL,SGBDI
C     EXTERNAL SGEXX,SMACH
C     BLAS SAXPY,SDOT,SSCAL,SASUM
C     FORTRAN ABS,AMAX1,FLOAT,MAX0,MIN0
C
C     INTERNAL VARIABLES
C
      REAL A(15,15),AB(43,15),AINV(15,15),ASAVE(15,15)
      REAL B(15),BT(15),SDOT,DET(2),DETB(2)
      REAL X(15),XB(15),XEXACT(15),XT(15),XTB(15),T,Z(15)
      REAL AINORM,ANORM,SMACH,COND,COND1,EN,ENORM,EPS
      REAL ETNORM,FNI,FNORM,ONEPX,RCOND,RCONDB,RNORM
      REAL RTNORM,Q(8),QS(8),SASUM,XNORM,XTNORM
      INTEGER I,IPVT(15),IPVTB(15),IQ(8),I1,I2,J
      INTEGER K,KASE,KB,KBFAIL,KOUNT,KP1,KSING,KSUSP(8)
      INTEGER L,LDA,LDAB,LUNIT,M,ML,MU,N,NM1,NPRINT
      LOGICAL KBF
C
      LDA = 15
      LDAB = 43
C
C     WRITE MATRIX AND SOLUTIONS IF  N .LE. NPRINT
C
      NPRINT = 3
C
      WRITE (LUNIT,460)
      WRITE (LUNIT,880)
C
      DO 10 I = 1, 8
         KSUSP(I) = 0
   10 CONTINUE
      KSING = 0
      KBFAIL = 0
C
C     SET EPS TO ROUNDING UNIT
C
      EPS = SMACH(1)
      WRITE (LUNIT,470) EPS
      WRITE (LUNIT,450)
C
C     START MAIN LOOP
C
      KASE = 1
   20 CONTINUE
C
C        GENERATE TEST MATRIX
C
         CALL SGEXX(A,LDA,N,KASE,LUNIT)
C
C        N = 0 SIGNALS NO MORE TEST MATRICES
C
C     ...EXIT
         IF (N .LE. 0) GO TO 440
         ANORM = 0.0E0
         DO 30 J = 1, N
            ANORM = AMAX1(ANORM,SASUM(N,A(1,J),1))
   30    CONTINUE
         WRITE (LUNIT,650) ANORM
C
         IF (N .GT. NPRINT) GO TO 50
            WRITE (LUNIT,450)
            DO 40 I = 1, N
               WRITE (LUNIT,700) (A(I,J), J = 1, N)
   40       CONTINUE
            WRITE (LUNIT,450)
   50    CONTINUE
C
C        GENERATE EXACT SOLUTION
C
         XEXACT(1) = 1.0E0
         IF (N .GE. 2) XEXACT(2) = 0.0E0
         IF (N .LE. 2) GO TO 70
            DO 60 I = 3, N
               XEXACT(I) = -XEXACT(I-2)
   60       CONTINUE
   70    CONTINUE
C
C        SAVE MATRIX AND GENERATE R.H.S.
C
         DO 90 I = 1, N
            B(I) = 0.0E0
            BT(I) = 0.0E0
            DO 80 J = 1, N
               ASAVE(I,J) = A(I,J)
               B(I) = B(I) + A(I,J)*XEXACT(J)
               BT(I) = BT(I) + A(J,I)*XEXACT(J)
   80       CONTINUE
            X(I) = B(I)
            XT(I) = BT(I)
            XB(I) = X(I)
            XTB(I) = XT(I)
   90    CONTINUE
C
C        FACTOR AND ESTIMATE CONDITION
C
         CALL SGECO(A,LDA,N,IPVT,RCOND,Z)
C
C        OUTPUT NULL VECTOR IF N .LE. NPRINT
C
         IF (N .GT. NPRINT) GO TO 110
            WRITE (LUNIT,720)
            DO 100 I = 1, N
               WRITE (LUNIT,730) Z(I)
  100       CONTINUE
            WRITE (LUNIT,450)
  110    CONTINUE
C
C        FACTOR BAND FORM AND COMPARE
C
         KBF = .FALSE.
         ML = 0
         MU = 0
         DO 140 J = 1, N
            DO 130 I = 1, N
               IF (ASAVE(I,J) .EQ. 0.0E0) GO TO 120
                  IF (I .LT. J) MU = MAX0(MU,J-I)
                  IF (I .GT. J) ML = MAX0(ML,I-J)
  120          CONTINUE
  130       CONTINUE
  140    CONTINUE
         WRITE (LUNIT,790) ML,MU
         IF (2*ML + MU + 1 .LE. LDAB) GO TO 150
            WRITE (LUNIT,680)
         GO TO 430
  150    CONTINUE
            M = ML + MU + 1
            DO 170 J = 1, N
               I1 = MAX0(1,J-MU)
               I2 = MIN0(N,J+ML)
               DO 160 I = I1, I2
                  K = I - J + M
                  AB(K,J) = ASAVE(I,J)
  160          CONTINUE
  170       CONTINUE
C
            CALL SGBCO(AB,LDAB,N,ML,MU,IPVTB,RCONDB,Z)
C
            IF (RCONDB .EQ. RCOND) GO TO 180
               WRITE (LUNIT,780)
               WRITE (LUNIT,820) RCOND,RCONDB
               KBF = .TRUE.
  180       CONTINUE
            KOUNT = 0
            DO 190 J = 1, N
               IF (AB(M,J) .NE. A(J,J)) KOUNT = KOUNT + 1
               IF (IPVTB(J) .NE. IPVT(J)) KOUNT = KOUNT + 1
  190       CONTINUE
            IF (KOUNT .EQ. 0) GO TO 200
               WRITE (LUNIT,780)
               WRITE (LUNIT,830) KOUNT
               KBF = .TRUE.
  200       CONTINUE
C
C           TEST FOR SINGULARITY
C
            IF (RCOND .GT. 0.0E0) GO TO 210
               WRITE (LUNIT,710) RCOND
               WRITE (LUNIT,480)
               KSING = KSING + 1
            GO TO 420
  210       CONTINUE
               COND = 1.0E0/RCOND
               WRITE (LUNIT,500) COND
               ONEPX = 1.0E0 + RCOND
               IF (ONEPX .EQ. 1.0E0) WRITE (LUNIT,490)
C
C              COMPUTE INVERSE, DETERMINANT AND COND1 = TRUE CONDITION
C
               DO 230 J = 1, N
                  DO 220 I = 1, N
                     AINV(I,J) = A(I,J)
  220             CONTINUE
  230          CONTINUE
               CALL SGEDI(AINV,LDA,N,IPVT,DET,Z,11)
               AINORM = 0.0E0
               DO 240 J = 1, N
                  AINORM = AMAX1(AINORM,SASUM(N,AINV(1,J),1))
  240          CONTINUE
               COND1 = ANORM*AINORM
               WRITE (LUNIT,510) COND1
               WRITE (LUNIT,750) DET(1)
               WRITE (LUNIT,760) DET(2)
C
C              SOLVE  A*X = B  AND  TRANS(A)*XT = BT
C
               CALL SGESL(A,LDA,N,IPVT,X,0)
               CALL SGESL(A,LDA,N,IPVT,XT,1)
C
               IF (N .GT. NPRINT) GO TO 270
                  WRITE (LUNIT,520)
                  DO 250 I = 1, N
                     WRITE (LUNIT,740) X(I)
  250             CONTINUE
                  WRITE (LUNIT,530)
                  DO 260 I = 1, N
                     WRITE (LUNIT,740) XT(I)
  260             CONTINUE
                  WRITE (LUNIT,450)
  270          CONTINUE
C
C              MORE BAND COMPARE
C
               CALL SGBSL(AB,LDAB,N,ML,MU,IPVTB,XB,0)
               CALL SGBSL(AB,LDAB,N,ML,MU,IPVTB,XTB,1)
               KOUNT = 0
               DO 280 I = 1, N
                  IF (XB(I) .NE. X(I)) KOUNT = KOUNT + 1
                  IF (XTB(I) .NE. XT(I)) KOUNT = KOUNT + 1
  280          CONTINUE
               IF (KOUNT .EQ. 0) GO TO 290
                  WRITE (LUNIT,780)
                  WRITE (LUNIT,840) KOUNT
                  KBF = .TRUE.
  290          CONTINUE
               CALL SGBDI(AB,LDAB,N,ML,MU,IPVTB,DETB)
               IF (DETB(1) .EQ. DET(1) .AND. DETB(2) .EQ. DET(2))
     *            GO TO 300
                  WRITE (LUNIT,780)
                  WRITE (LUNIT,850) DETB
                  KBF = .TRUE.
  300          CONTINUE
C
C              RECONSTRUCT  A  FROM TRIANGULAR FACTORS , L AND U
C
               NM1 = N - 1
               IF (NM1 .LT. 1) GO TO 330
               DO 320 KB = 1, NM1
                  K = N - KB
                  KP1 = K + 1
                  L = IPVT(K)
                  DO 310 J = KP1, N
                     T = -A(K,J)
                     CALL SAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
                     T = A(L,J)
                     A(L,J) = A(K,J)
                     A(K,J) = T
  310             CONTINUE
                  T = -A(K,K)
                  CALL SSCAL(N-K,T,A(K+1,K),1)
                  T = A(L,K)
                  A(L,K) = A(K,K)
                  A(K,K) = T
  320          CONTINUE
  330          CONTINUE
C
C              COMPUTE ERRORS AND RESIDUALS
C                 E  =  X - XEXACT
C                 ET =  XT - XEXACT
C                 R  =  B - A*X
C                 RT =  BT - A*XT
C                 F  =  A - L*U
C                 AI =  A*INV(A) - I
C
               XNORM = SASUM(N,X,1)
               XTNORM = SASUM(N,XT,1)
               ENORM = 0.0E0
               ETNORM = 0.0E0
               FNORM = 0.0E0
               DO 350 J = 1, N
                  ENORM = ENORM + ABS(X(J)-XEXACT(J))
                  ETNORM = ETNORM + ABS(XT(J)-XEXACT(J))
                  T = -X(J)
                  CALL SAXPY(N,T,ASAVE(1,J),1,B,1)
                  BT(J) = BT(J) - SDOT(N,ASAVE(1,J),1,XT,1)
                  FNI = 0.0E0
                  DO 340 I = 1, N
                     FNI = FNI + ABS(ASAVE(I,J)-A(I,J))
  340             CONTINUE
                  IF (FNI .GT. FNORM) FNORM = FNI
  350          CONTINUE
               RNORM = SASUM(N,B,1)
               RTNORM = SASUM(N,BT,1)
C
C              A*INV(A) - I
C
               AINORM = 0.0E0
               DO 380 J = 1, N
                  DO 360 I = 1, N
                     B(I) = 0.0E0
  360             CONTINUE
                  DO 370 K = 1, N
                     T = AINV(K,J)
                     CALL SAXPY(N,T,ASAVE(1,K),1,B,1)
  370             CONTINUE
                  B(J) = B(J) - 1.0E0
                  AINORM = AMAX1(AINORM,SASUM(N,B,1))
  380          CONTINUE
C
               WRITE (LUNIT,540) ENORM,ETNORM
               WRITE (LUNIT,550) RNORM,RTNORM
               WRITE (LUNIT,660) FNORM
               WRITE (LUNIT,670) AINORM
C
C              COMPUTE TEST RATIOS
C
               Q(1) = COND/COND1
               Q(2) = COND1/COND
               Q(3) = ENORM/(EPS*COND*XNORM)
               Q(4) = ETNORM/(EPS*COND*XTNORM)
               Q(5) = RNORM/(EPS*ANORM*XNORM)
               Q(6) = RTNORM/(EPS*ANORM*XTNORM)
               Q(7) = FNORM/(EPS*ANORM)
               Q(8) = AINORM/(EPS*COND)
               WRITE (LUNIT,450)
               WRITE (LUNIT,560)
               WRITE (LUNIT,450)
               WRITE (LUNIT,620)
               WRITE (LUNIT,630)
               WRITE (LUNIT,640)
               WRITE (LUNIT,450)
               WRITE (LUNIT,690) (Q(I), I = 1, 8)
               WRITE (LUNIT,450)
C
C              LOOK FOR SUSPICIOUS RATIOS
C
               QS(1) = 1.0E0 + 4.0E0*EPS
               QS(2) = 10.0E0
               EN = FLOAT(N)
               IF (N .EQ. 1) EN = 2.0E0
               DO 390 I = 3, 8
                  QS(I) = EN
  390          CONTINUE
               KOUNT = 0
               DO 410 I = 1, 8
                  IQ(I) = 0
                  IF (Q(I) .LE. QS(I)) GO TO 400
                     IQ(I) = 1
                     KSUSP(I) = KSUSP(I) + 1
                     KOUNT = KOUNT + 1
  400             CONTINUE
  410          CONTINUE
               IF (KOUNT .EQ. 0) WRITE (LUNIT,860)
               IF (KOUNT .NE. 0) WRITE (LUNIT,870) (IQ(I), I = 1, 8)
               WRITE (LUNIT,450)
  420       CONTINUE
  430    CONTINUE
C
         IF (.NOT.KBF) WRITE (LUNIT,770)
         IF (KBF) KBFAIL = KBFAIL + 1
         WRITE (LUNIT,570)
         KASE = KASE + 1
      GO TO 20
  440 CONTINUE
C
C     FINISH MAIN LOOP
C
C     SUMMARY
C
      WRITE (LUNIT,580)
      KASE = KASE - 1
      WRITE (LUNIT,590) KASE
      WRITE (LUNIT,600) KSING
      WRITE (LUNIT,800) KBFAIL
      WRITE (LUNIT,610) KSUSP
      WRITE (LUNIT,810)
      RETURN
C
C     MOST FORMATS, ALSO SOME IN SGEXX
C
  450 FORMAT (1H )
  460 FORMAT (29H1LINPACK TESTER, SGE**, SGB**)
  470 FORMAT ( / 14H EPSILON     =, 1PE13.5)
  480 FORMAT ( / 19H EXACT SINGULARITY. /)
  490 FORMAT ( / 16H MAYBE SINGULAR. /)
  500 FORMAT (14H COND        =, 1PE13.5)
  510 FORMAT (14H ACTUAL COND =, 1PE13.5)
  520 FORMAT ( / 4H X =)
  530 FORMAT ( / 5H XT =)
  540 FORMAT (14H ERROR NORMS =, 1P2E13.5)
  550 FORMAT (14H RESID NORMS =, 1P2E13.5)
  560 FORMAT (26H TEST RATIOS.. E = EPSILON)
  570 FORMAT ( / 14H ************* /)
  580 FORMAT (8H1SUMMARY)
  590 FORMAT (18H NUMBER OF TESTS =, I4)
  600 FORMAT (30H NUMBER OF SINGULAR MATRICES =, I4)
  610 FORMAT (30H NUMBER OF SUSPICIOUS RATIOS =, 8I4)
  620 FORMAT (30H     COND     ACTUAL    ERROR ,
     *        50H   ERROR-T    RESID    RESID-T    A - LU   A*AI-I )
  630 FORMAT (8(10H   -------))
  640 FORMAT (30H    ACTUAL     COND   E*COND*X,
     *        50H  E*COND*X    E*A*X     E*A*X      E*A     E*COND )
  650 FORMAT (14H NORM(A)     =, 1PE13.5)
  660 FORMAT (14H NORM(A - LU)=, 1PE13.5)
  670 FORMAT (14H NORM(A*AI-I)=, 1PE13.5)
  680 FORMAT ( / 19H BAND WIDTH TOO BIG)
  690 FORMAT (8(1X, F9.4))
  700 FORMAT (1H , 6G11.4)
  710 FORMAT (14H 1/COND      =, 1PE13.5)
  720 FORMAT ( / 7H NULL =)
  730 FORMAT (2G14.6)
  740 FORMAT (2G14.6)
  750 FORMAT (14H DET FRACT   =, 2F9.5)
  760 FORMAT (14H DET EXPON   =, 2F9.0)
  770 FORMAT ( / 20H BAND ROUTINES AGREE /)
  780 FORMAT ( / 28H BAND ROUTINES DO NOT AGREE,)
  790 FORMAT (5H ML =, I2, 6H  MU =, I2)
  800 FORMAT (26H NUMBER OF BAND FAILURES =, I4)
  810 FORMAT ( / 12H END OF TEST)
  820 FORMAT (8H RCOND =, 1P2E13.5 /)
  830 FORMAT (12H KOUNT(FA) =, I4 /)
  840 FORMAT (12H KOUNT(SL) =, I4 /)
  850 FORMAT (8H DET   =, 4F9.5 /)
  860 FORMAT (21H NO SUSPICIOUS RATIOS)
  870 FORMAT (I8, 7I10 / 7X, 28H1 INDICATES SUSPICIOUS RATIO)
  880 FORMAT (29H THIS VERSION DATED 08/14/78.)
      END
      SUBROUTINE SGEXX(A,LDA,N,KASE,LUNIT)
      INTEGER LDA,N,KASE,LUNIT
      REAL A(LDA,1)
C
C     GENERATES REAL GENERAL TEST MATRICES
C
C     EXTERNAL SMACH
C     FORTRAN FLOAT,MAX0
      REAL T1,T2
      REAL SMACH,HUGE,TINY
      INTEGER I,J
C
      GO TO (10, 10, 10, 60, 60, 80, 80, 80, 120, 160, 200, 240, 280,
     *       320, 360, 410, 460), KASE
C
C     KASE 1, 2 AND 3
C
   10 CONTINUE
         N = 3*KASE
         WRITE (LUNIT,20) KASE,N
   20    FORMAT (5H KASE, I3, 3X, 16HHILBERT SLICE    / 4H N =, I4)
         DO 50 J = 1, N
            DO 40 I = 1, N
               A(I,J) = 0.0E0
               IF (I .GT. J + 2) GO TO 30
               IF (I .LT. J - 3) GO TO 30
                  A(I,J) = 1.0E0/FLOAT(I+J-1)
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      GO TO 470
C
C     KASE 4 AND 5
C
   60 CONTINUE
         N = 1
         WRITE (LUNIT,70) KASE,N
   70    FORMAT (5H KASE, I3, 3X, 16HMONOELEMENTAL    / 4H N =, I4)
         IF (KASE .EQ. 4) A(1,1) = 3.0E0
         IF (KASE .EQ. 5) A(1,1) = 0.0E0
      GO TO 470
C
C     KASE 6, 7 AND 8
C
   80 CONTINUE
         N = 15
         WRITE (LUNIT,90) KASE,N
   90    FORMAT (5H KASE, I3, 3X, 16HTRIDIAGONAL      / 4H N =, I4)
         T1 = 1.0E0
         T2 = 1.0E0
         IF (KASE .EQ. 7) T1 = 100.0E0
         IF (KASE .EQ. 8) T2 = 100.0E0
         DO 110 I = 1, N
            DO 100 J = 1, N
               A(I,J) = 0.0E0
               IF (I .EQ. J) A(I,I) = 4.0E0
               IF (I .EQ. J - 1) A(I,J) = T1
               IF (I .EQ. J + 1) A(I,J) = T2
  100       CONTINUE
  110    CONTINUE
      GO TO 470
C
C     KASE 9
C
  120 CONTINUE
         N = 5
         WRITE (LUNIT,130) KASE,N
  130    FORMAT (5H KASE, I3, 3X, 16HRANK ONE         / 4H N =, I4)
         DO 150 I = 1, N
            DO 140 J = 1, N
               A(I,J) = 10.0E0**(I - J)
  140       CONTINUE
  150    CONTINUE
      GO TO 470
C
C     KASE 10
C
  160 CONTINUE
         N = 4
         WRITE (LUNIT,170) KASE,N
  170    FORMAT (5H KASE, I3, 3X, 16HZERO COLUMN      / 4H N =, I4)
         DO 190 I = 1, N
            DO 180 J = 1, N
               T1 = FLOAT(J-3)
               T2 = FLOAT(I)
               A(I,J) = T1/T2
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
               IF (I .EQ. J) A(I,J) = FLOAT(I)
               IF (I .GT. J) A(I,J) = FLOAT(J-2)
               IF (I .LT. J) A(I,J) = FLOAT(I-2)
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
               IF (I .EQ. J) A(I,I) = 1.0E0
               IF (I .NE. J) A(I,J) = 0.0E0
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
               IF (I .GT. J) A(I,J) = 0.0E0
               IF (I .LE. J) A(I,J) = FLOAT(J-I+1)
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
               IF (I .LT. J) A(I,J) = 0.0E0
               IF (I .GE. J) A(I,J) = FLOAT(I-J+1)
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
         TINY = SMACH(2)
         WRITE (LUNIT,380) TINY
  380    FORMAT (14H TINY        =, 1PE13.5)
         DO 400 I = 1, N
            DO 390 J = 1, N
               A(I,J) = TINY*FLOAT(J)/FLOAT(MAX0(I,J))
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
         HUGE = SMACH(3)
         WRITE (LUNIT,430) HUGE
  430    FORMAT (14H HUGE        =, 1PE13.5)
         DO 450 I = 1, N
            DO 440 J = 1, N
               A(I,J) = HUGE*FLOAT(J)/FLOAT(MAX0(I,J))
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
