C     MAIN PROGRAM
      INTEGER LUNIT
C     ALLOW 5000 UNDERFLOWS.
      CALL TRAPS(0,0,5001,0,0)
C
C     OUTPUT UNIT NUMBER
C
      LUNIT = 6
C
      CALL DPOTS(LUNIT)
      STOP
      END
      SUBROUTINE DPOTS(LUNIT)
C     LUNIT IS THE OUTPUT UNIT NUMBER
C
C     TESTS
C        DPOCO,DPOFA,DPOSL,DPODI,DPPCO,DPPFA,DPPSL,DPPDI
C        DPBCO,DPBFA,DPBSL,DPBDI
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK DPOCO,DPOSL,DPODI,DPPCO,DPPSL,DPPDI
C     LINPACK DPBCO,DPBSL,DPBDI
C     EXTERNAL DPOXX,DMACH
C     BLAS DAXPY,DDOT,DASUM
C     FORTRAN DABS,DMAX1,DBLE,FLOAT,MAX0
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION A(15,15),AB(15,15),AINV(15,15),ASAVE(15,15)
      DOUBLE PRECISION AP(120),B(15),DDOT,X(15),XB(15),XEXACT(15)
      DOUBLE PRECISION XP(15),T,Z(15)
      DOUBLE PRECISION ANORM,AINORM,COND,COND1,DET(2),DETB(2),DETP(2)
      DOUBLE PRECISION EN,ENORM,EPS,FNORM,ONEPX,Q(6),QS(6),RCOND,RCONDB
      DOUBLE PRECISION RCONDP,RNORM,S,DASUM,DMACH,XNORM
      INTEGER I,INFO,INFOB,INFOP,IQ(6),I1,J,JB
      INTEGER K,KASE,KB,KBFAIL,KNPD,KOUNT,KPFAIL
      INTEGER KSUSP(6),LDA,LUNIT,M,N,NPRINT
      LOGICAL KBF,KPF
C
      LDA = 15
C
C     WRITE MATRIX AND SOLUTIONS IF  N .LE. NPRINT
C
      NPRINT = 3
C
      WRITE (LUNIT,560)
      WRITE (LUNIT,1000)
C
      DO 10 I = 1, 6
         KSUSP(I) = 0
   10 CONTINUE
      KNPD = 0
      KPFAIL = 0
      KBFAIL = 0
C
C     SET EPS TO ROUNDING UNIT FOR DOUBLE PRECISION ARITHMETIC
C
      EPS = DMACH(1)
      WRITE (LUNIT,570) EPS
      WRITE (LUNIT,550)
C
C     START MAIN LOOP
C
      KASE = 1
   20 CONTINUE
C
C        GENERATE TEST MATRIX
C
         CALL DPOXX(A,LDA,N,KASE,LUNIT)
C
C        N = 0 SIGNALS NO MORE TEST MATRICES
C
C     ...EXIT
         IF (N .LE. 0) GO TO 540
         ANORM = 0.0D0
         DO 30 J = 1, N
            ANORM = DMAX1(ANORM,DASUM(N,A(1,J),1))
   30    CONTINUE
         WRITE (LUNIT,720) ANORM
C
         IF (N .GT. NPRINT) GO TO 50
            WRITE (LUNIT,550)
            DO 40 I = 1, N
               WRITE (LUNIT,760) (A(I,J), J = 1, N)
   40       CONTINUE
            WRITE (LUNIT,550)
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
C        SAVE MATRIX AND GENERATE R.H.S.
C
         DO 90 I = 1, N
            B(I) = 0.0D0
            DO 80 J = 1, N
               ASAVE(I,J) = A(I,J)
               B(I) = B(I) + A(I,J)*XEXACT(J)
   80       CONTINUE
            X(I) = B(I)
            XP(I) = X(I)
            XB(I) = X(I)
   90    CONTINUE
C
C        FACTOR AND ESTIMATE CONDITION
C
         RCOND = -1.0D0
         CALL DPOCO(A,LDA,N,RCOND,Z,INFO)
C
C        OUTPUT NULL VECTOR
C
         IF (N .GT. NPRINT .OR. INFO .NE. 0) GO TO 110
            WRITE (LUNIT,770)
            DO 100 I = 1, N
               WRITE (LUNIT,780) Z(I)
  100       CONTINUE
            WRITE (LUNIT,550)
  110    CONTINUE
C
C        FACTOR PACKED FORM AND COMPARE
C
         KPF = .FALSE.
         K = 0
         DO 130 J = 1, N
            DO 120 I = 1, J
               K = K + 1
               AP(K) = ASAVE(I,J)
  120       CONTINUE
  130    CONTINUE
         RCONDP = -1.0D0
         CALL DPPCO(AP,N,RCONDP,Z,INFOP)
         IF (INFOP .EQ. INFO) GO TO 140
            WRITE (LUNIT,880)
            WRITE (LUNIT,920) INFO,INFOP
            KPF = .TRUE.
  140    CONTINUE
         IF (RCONDP .EQ. RCOND) GO TO 150
            WRITE (LUNIT,880)
            WRITE (LUNIT,930) RCOND,RCONDP
            KPF = .TRUE.
  150    CONTINUE
         K = 0
         KOUNT = 0
         DO 170 J = 1, N
            DO 160 I = 1, J
               K = K + 1
               IF (AP(K) .NE. A(I,J)) KOUNT = KOUNT + 1
  160       CONTINUE
  170    CONTINUE
         IF (KOUNT .EQ. 0) GO TO 180
            WRITE (LUNIT,880)
            WRITE (LUNIT,940) KOUNT
            KPF = .TRUE.
  180    CONTINUE
C
C        FACTOR BAND FORM AND COMPARE
C
         KBF = .FALSE.
         M = 0
         DO 200 J = 1, N
            DO 190 I = 1, J
               IF (ASAVE(I,J) .NE. 0.0D0) M = MAX0(M,J-I)
  190       CONTINUE
  200    CONTINUE
C
         DO 220 J = 1, N
            I1 = MAX0(1,J-M)
            DO 210 I = I1, J
               K = I - J + M + 1
               AB(K,J) = ASAVE(I,J)
  210       CONTINUE
  220    CONTINUE
         WRITE (LUNIT,840) M
         RCONDB = -1.0D0
         CALL DPBCO(AB,LDA,N,M,RCONDB,Z,INFOB)
         IF (INFOB .EQ. INFO) GO TO 230
            WRITE (LUNIT,830)
            WRITE (LUNIT,920) INFO,INFOB
            KBF = .TRUE.
  230    CONTINUE
         IF (RCONDB .EQ. RCOND) GO TO 240
            WRITE (LUNIT,830)
            WRITE (LUNIT,930) RCOND,RCONDB
            KBF = .TRUE.
  240    CONTINUE
         KOUNT = 0
         DO 250 J = 1, N
            IF (AB(M+1,J) .NE. A(J,J)) KOUNT = KOUNT + 1
  250    CONTINUE
         IF (KOUNT .EQ. 0) GO TO 260
            WRITE (LUNIT,830)
            WRITE (LUNIT,940) KOUNT
            KBF = .TRUE.
  260    CONTINUE
C
C        TEST FOR DEFINITENESS
C
         IF (INFO .EQ. 0) GO TO 270
            WRITE (LUNIT,860) INFO
            KNPD = KNPD + 1
         GO TO 530
  270    CONTINUE
            COND = 1.0D0/RCOND
            WRITE (LUNIT,590) COND
            ONEPX = 1.0D0 + RCOND
            IF (ONEPX .EQ. 1.0D0) WRITE (LUNIT,580)
C
C           COMPUTE INVERSE, DETERMINANT AND COND1 = TRUE CONDITION
C
            DO 290 J = 1, N
               DO 280 I = 1, J
                  AINV(I,J) = A(I,J)
  280          CONTINUE
  290       CONTINUE
            CALL DPODI(AINV,LDA,N,DET,11)
            AINORM = 0.0D0
            DO 310 J = 1, N
               DO 300 I = J, N
                  AINV(I,J) = AINV(J,I)
  300          CONTINUE
               AINORM = DMAX1(AINORM,DASUM(N,AINV(1,J),1))
  310       CONTINUE
            COND1 = ANORM*AINORM
            WRITE (LUNIT,600) COND1
            WRITE (LUNIT,800) DET(1)
            WRITE (LUNIT,810) DET(2)
C
C           SOLVE A*X = B
C
            CALL DPOSL(A,LDA,N,X)
C
            IF (N .GT. NPRINT) GO TO 330
               WRITE (LUNIT,610)
               DO 320 I = 1, N
                  WRITE (LUNIT,790) X(I)
  320          CONTINUE
               WRITE (LUNIT,550)
  330       CONTINUE
C
C           MORE PACKED COMPARE
C
            CALL DPPSL(AP,N,XP)
            KOUNT = 0
            DO 340 I = 1, N
               IF (XP(I) .NE. X(I)) KOUNT = KOUNT + 1
  340       CONTINUE
            IF (KOUNT .EQ. 0) GO TO 350
               WRITE (LUNIT,880)
               WRITE (LUNIT,950) KOUNT
               KPF = .TRUE.
  350       CONTINUE
            CALL DPPDI(AP,N,DETP,11)
            IF (DETP(1) .EQ. DET(1) .AND. DETP(2) .EQ. DET(2))
     *         GO TO 360
               WRITE (LUNIT,880)
               WRITE (LUNIT,960) DETP
               KPF = .TRUE.
  360       CONTINUE
            KOUNT = 0
            K = 0
            DO 380 J = 1, N
               DO 370 I = 1, J
                  K = K + 1
                  IF (AP(K) .NE. AINV(I,J)) KOUNT = KOUNT + 1
  370          CONTINUE
  380       CONTINUE
            IF (KOUNT .EQ. 0) GO TO 390
               WRITE (LUNIT,880)
               WRITE (LUNIT,970) KOUNT
               KPF = .TRUE.
  390       CONTINUE
C
C           MORE BAND COMPARE
C
            CALL DPBSL(AB,LDA,N,M,XB)
            KOUNT = 0
            DO 400 I = 1, N
               IF (XB(I) .NE. X(I)) KOUNT = KOUNT + 1
  400       CONTINUE
            IF (KOUNT .EQ. 0) GO TO 410
               WRITE (LUNIT,830)
               WRITE (LUNIT,950) KOUNT
               KBF = .TRUE.
  410       CONTINUE
            CALL DPBDI(AB,LDA,N,M,DETB)
            IF (DETB(1) .EQ. DET(1) .AND. DETB(2) .EQ. DET(2))
     *         GO TO 420
               WRITE (LUNIT,830)
               WRITE (LUNIT,960) DETB
               KBF = .TRUE.
  420       CONTINUE
C
C           RECONSTRUCT  A  FROM TRIANGULAR FACTORS , TRANS(R) AND R
C
            DO 440 JB = 1, N
               J = N + 1 - JB
               DO 430 KB = 1, J
                  K = J + 1 - KB
                  A(K,J) = DDOT(K,A(1,K),1,A(1,J),1)
  430          CONTINUE
  440       CONTINUE
C
C           COMPUTE ERRORS AND RESIDUALS
C              E  =  X - XEXACT
C              R  =  B - A*X
C              F  =  A - TRANS(R)*R
C
            XNORM = DASUM(N,X,1)
            ENORM = 0.0D0
            FNORM = 0.0D0
            DO 460 J = 1, N
               ENORM = ENORM + DABS(X(J)-XEXACT(J))
               T = -X(J)
               CALL DAXPY(N,T,ASAVE(1,J),1,B,1)
               S = 0.0D0
               DO 450 I = 1, J
                  S = S + DABS(ASAVE(I,J)-A(I,J))
  450          CONTINUE
               IF (S .GT. FNORM) FNORM = S
  460       CONTINUE
            RNORM = DASUM(N,B,1)
C
C           A*INV(A) - I
C
            AINORM = 0.0D0
            DO 490 J = 1, N
               DO 470 I = 1, N
                  B(I) = 0.0D0
  470          CONTINUE
               DO 480 K = 1, N
                  T = AINV(K,J)
                  CALL DAXPY(N,T,ASAVE(1,K),1,B,1)
  480          CONTINUE
               B(J) = B(J) - 1.0D0
               AINORM = DMAX1(AINORM,DASUM(N,B,1))
  490       CONTINUE
C
            WRITE (LUNIT,620) ENORM
            WRITE (LUNIT,630) RNORM
            WRITE (LUNIT,730) FNORM
            WRITE (LUNIT,740) AINORM
C
C           COMPUTE TEST RATIOS
C
            Q(1) = COND/COND1
            Q(2) = COND1/COND
            Q(3) = ENORM/(EPS*COND*XNORM)
            Q(4) = RNORM/(EPS*ANORM*XNORM)
            Q(5) = FNORM/(EPS*ANORM)
            Q(6) = AINORM/(EPS*COND)
            WRITE (LUNIT,550)
            WRITE (LUNIT,640)
            WRITE (LUNIT,550)
            WRITE (LUNIT,690)
            WRITE (LUNIT,700)
            WRITE (LUNIT,710)
            WRITE (LUNIT,550)
            WRITE (LUNIT,750) (Q(I), I = 1, 6)
            WRITE (LUNIT,550)
C
C           LOOK FOR SUSPICIOUS RATIOS
C
            QS(1) = 1.0D0 + 4.0D0*EPS
            QS(2) = 10.0D0
            EN = DBLE(FLOAT(N))
            IF (N .EQ. 1) EN = 2.0D0
            DO 500 I = 3, 6
               QS(I) = EN
  500       CONTINUE
            KOUNT = 0
            DO 520 I = 1, 6
               IQ(I) = 0
               IF (Q(I) .LE. QS(I)) GO TO 510
                  IQ(I) = 1
                  KSUSP(I) = KSUSP(I) + 1
                  KOUNT = KOUNT + 1
  510          CONTINUE
  520       CONTINUE
            IF (KOUNT .EQ. 0) WRITE (LUNIT,980)
            IF (KOUNT .NE. 0) WRITE (LUNIT,990) (IQ(I), I = 1, 6)
            WRITE (LUNIT,550)
  530    CONTINUE
C
         IF (.NOT.KPF) WRITE (LUNIT,870)
         IF (KPF) KPFAIL = KPFAIL + 1
         IF (.NOT.KBF) WRITE (LUNIT,820)
         IF (KBF) KBFAIL = KBFAIL + 1
         WRITE (LUNIT,650)
         KASE = KASE + 1
      GO TO 20
  540 CONTINUE
C
C     FINISH MAIN LOOP
C
C     SUMMARY
C
      WRITE (LUNIT,660)
      KASE = KASE - 1
      WRITE (LUNIT,670) KASE
      WRITE (LUNIT,900) KNPD
      WRITE (LUNIT,890) KPFAIL
      WRITE (LUNIT,850) KBFAIL
      WRITE (LUNIT,680) KSUSP
      WRITE (LUNIT,910)
      RETURN
C
C     MOST FORMATS, ALSO SOME IN DPOXX
C
  550 FORMAT (1H )
  560 FORMAT (36H1LINPACK TESTER, DPO**, DPP**, DPB**)
  570 FORMAT ( / 14H EPSILON     =, 1PD13.5)
  580 FORMAT ( / 16H MAYBE SINGULAR. /)
  590 FORMAT (14H COND        =, 1PD13.5)
  600 FORMAT (14H ACTUAL COND =, 1PD13.5)
  610 FORMAT ( / 4H X =)
  620 FORMAT (14H ERROR NORM  =, 1P2D13.5)
  630 FORMAT (14H RESID NORM  =, 1P2D13.5)
  640 FORMAT (26H TEST RATIOS.. E = EPSILON)
  650 FORMAT ( / 14H ************* /)
  660 FORMAT (8H1SUMMARY)
  670 FORMAT (18H NUMBER OF TESTS =, I4)
  680 FORMAT ( / 30H NUMBER OF SUSPICIOUS RATIOS =, 6I4)
  690 FORMAT (30H     COND     ACTUAL    ERROR ,
     *        30H    RESID   A - RT*R  A*AI - I)
  700 FORMAT (6(10H   -------))
  710 FORMAT (30H    ACTUAL     COND   E*COND*X,
     *        30H    E*A*X      E*A     E*COND )
  720 FORMAT (14H NORM(A)     =, 1PD13.5)
  730 FORMAT (14H NORM(A-RT*R)=, 1PD13.5)
  740 FORMAT (14H NORM(A*AI-I)=, 1PD13.5)
  750 FORMAT (6(1X, F9.4))
  760 FORMAT (1H , 6G11.4)
  770 FORMAT ( / 7H NULL =)
  780 FORMAT (2G14.6)
  790 FORMAT (2G14.6)
  800 FORMAT (14H DET FRACT   =, F9.5)
  810 FORMAT (14H DET EXPON   =, F9.0)
  820 FORMAT ( / 20H BAND ROUTINES AGREE /)
  830 FORMAT ( / 28H BAND ROUTINES DO NOT AGREE,)
  840 FORMAT (5H M  =, I2)
  850 FORMAT (26H NUMBER OF BAND FAILURES =, I4)
  860 FORMAT (30H NOT POSITIVE DEFINITE, INFO =, I2)
  870 FORMAT ( / 22H PACKED ROUTINES AGREE)
  880 FORMAT ( / 30H PACKED ROUTINES DO NOT AGREE,)
  890 FORMAT (28H NUMBER OF PACKED FAILURES =, I4)
  900 FORMAT (34H NUMBER OF NOT POSITIVE DEFINITE =, I4)
  910 FORMAT ( / 12H END OF TEST)
  920 FORMAT (8H INFO  =, 2I3)
  930 FORMAT (8H RCOND =, 1P2D13.5)
  940 FORMAT (12H KOUNT(FA) =, I4)
  950 FORMAT (12H KOUNT(SL) =, I4)
  960 FORMAT (8H DET   =, F9.5, F9.0)
  970 FORMAT (12H KOUNT(DI) =, I4)
  980 FORMAT (21H NO SUSPICIOUS RATIOS)
  990 FORMAT (I8, 5I10 / 7X, 28H1 INDICATES SUSPICIOUS RATIO)
 1000 FORMAT (29H THIS VERSION DATED 08/14/78.)
      END
      SUBROUTINE DPOXX(A,LDA,N,KASE,LUNIT)
      INTEGER LDA,N,KASE,LUNIT
      DOUBLE PRECISION A(LDA,1)
C
C     GENERATES DOUBLE PRECISION POSITIVE DEFINITE TEST MATRICES
C
C     EXTERNAL DMACH
C     FORTRAN DABS,DBLE,FLOAT,IABS,MAX0,MIN0
      DOUBLE PRECISION T
      DOUBLE PRECISION TINY,HUGE,DMACH
      INTEGER I,J
C
      GO TO (10, 10, 10, 50, 50, 70, 70, 70, 120, 160, 200, 240, 290,
     *       340), KASE
C
C     KASE 1, 2 AND 3
C
   10 CONTINUE
         N = 5*KASE
         WRITE (LUNIT,20) KASE,N
   20    FORMAT (5H KASE, I3, 3X, 16HHILBERT          / 4H N =, I4)
         T = 1.0D0
         T = DSIGN(1.0D0,T)
         DO 40 J = 1, N
            DO 30 I = 1, N
               A(I,J) = T**(I - J)/DBLE(FLOAT(I+J-1))
C              FOR DOUBLE PRECISION MATRICES, A(I,J) =
C              1.0/DBLEFLOAT(I+J-1)
   30       CONTINUE
   40    CONTINUE
      GO TO 350
C
C     KASE 4 AND 5
C
   50 CONTINUE
         N = 1
         WRITE (LUNIT,60) KASE,N
   60    FORMAT (5H KASE, I3, 3X, 16HMONOELEMENTAL    / 4H N =, I4)
         IF (KASE .EQ. 4) A(1,1) = 3.0D0
         IF (KASE .EQ. 5) A(1,1) = 0.0D0
      GO TO 350
C
C     KASE 6, 7 AND 8
C
   70 CONTINUE
         N = 15
         IF (KASE .NE. 8) WRITE (LUNIT,80) KASE,N
   80    FORMAT (5H KASE, I3, 3X, 16HTRIDIAGONAL      / 4H N =, I4)
         IF (KASE .EQ. 8) WRITE (LUNIT,90) KASE,N
   90    FORMAT (5H KASE, I3, 3X, 16HDIAGONAL         / 4H N =, I4)
         T = 1.0D0
         IF (KASE .EQ. 7) T = 2.0D0
         IF (KASE .EQ. 8) T = 0.0D0
         DO 110 J = 1, N
            DO 100 I = 1, J
               A(I,J) = 0.0D0
               IF (I .EQ. J) A(I,I) = 4.0D0
               IF (I .EQ. J - 1) A(I,J) = T
               A(J,I) = A(I,J)
  100       CONTINUE
  110    CONTINUE
      GO TO 350
C
C     KASE 9
C
  120 CONTINUE
         N = 5
         WRITE (LUNIT,130) KASE,N
  130    FORMAT (5H KASE, I3, 3X, 16HPENTADIAGONAL    / 4H N =, I4)
         DO 150 J = 1, N
            DO 140 I = 1, N
               A(I,J) = 0.0D0
               IF (IABS(I-J) .LE. 2)
     *            A(I,J) = (5.0D0 - DBLE(FLOAT(IABS(I-J))))
     *                     **(10 - I - J)
  140       CONTINUE
  150    CONTINUE
      GO TO 350
C
C     KASE 10
C
  160 CONTINUE
         N = 6
         WRITE (LUNIT,170) KASE,N
  170    FORMAT (5H KASE, I3, 3X, 16HTRIDIAG INVERSE  / 4H N =, I4)
         DO 190 J = 1, N
            DO 180 I = 1, J
               A(I,J) = DBLE(FLOAT(N+1-J))
               A(J,I) = A(I,J)
  180       CONTINUE
  190    CONTINUE
      GO TO 350
C
C     KASE 11
C
  200 CONTINUE
         N = 15
         WRITE (LUNIT,210) KASE,N
  210    FORMAT (5H KASE, I3, 3X, 16HTEST COND        / 4H N =, I4)
         DO 230 J = 1, N
            DO 220 I = 1, N
               IF (I .EQ. J) A(I,J) = DBLE(FLOAT(I))
               IF (I .GT. J) A(I,J) = DBLE(FLOAT(J-2))
               IF (I .LT. J) A(I,J) = DBLE(FLOAT(I-2))
  220       CONTINUE
  230    CONTINUE
      GO TO 350
C
C     KASE 12
C
  240 CONTINUE
         N = 5
         WRITE (LUNIT,250) KASE,N
  250    FORMAT (5H KASE, I3, 3X, 16HNEAR UNDERFLOW   / 4H N =, I4)
         TINY = DMACH(2)
         WRITE (LUNIT,260) TINY
  260    FORMAT (14H TINY        =, 1PD13.5)
         DO 280 I = 1, N
            DO 270 J = 1, N
               A(I,J) = TINY*DBLE(FLOAT(MIN0(I,J)))
     *                  /DBLE(FLOAT(MAX0(I,J)))
  270       CONTINUE
  280    CONTINUE
      GO TO 350
C
C     KASE 13
C
  290 CONTINUE
         N = 5
         WRITE (LUNIT,300) KASE,N
  300    FORMAT (5H KASE, I3, 3X, 16HNEAR OVERFLOW    / 4H N =, I4)
         HUGE = DMACH(3)
         WRITE (LUNIT,310) HUGE
  310    FORMAT (14H HUGE        =, 1PD13.5)
         DO 330 I = 1, N
            DO 320 J = 1, N
               A(I,J) = HUGE*DBLE(FLOAT(MIN0(I,J)))
     *                  /DBLE(FLOAT(MAX0(I,J)))
  320       CONTINUE
  330    CONTINUE
      GO TO 350
C
  340 CONTINUE
         N = 0
  350 CONTINUE
      RETURN
C
      END
