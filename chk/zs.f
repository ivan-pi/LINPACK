C     MAIN PROGRAM
      INTEGER LUNIT
      DOUBLE PRECISION DREAL,DIMAG
      COMPLEX*16 ZDUMR,ZDUMI
      DREAL(ZDUMR) = ZDUMR
      DIMAG(ZDUMI) = (0.0D0,-1.0D0)*ZDUMI
C     ALLOW 5000 UNDERFLOWS, 2 DIVISIONS BY ZERO.
C     THE DIVISIONS BY ZERO SHOULD NOT OCCUR IF COMPLEX*16 DIVISION IS
C     OK.
      CALL TRAPS(0,0,5001,0,3)
C
C     OUTPUT UNIT NUMBER
C
      LUNIT = 6
C
      CALL ZSITS(LUNIT)
      STOP
      END
      SUBROUTINE ZSITS(LUNIT)
C     LUNIT IS THE OUTPUT UNIT NUMBER
C
C     TESTS
C        ZSICO,ZSIFA,ZSISL,ZSIDI,ZSPCO,ZSPFA,ZSPSL,ZSPDI
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK ZSICO,ZSISL,ZSIDI,ZSPCO,ZSPSL,ZSPDI
C     EXTERNAL ZSIXX,ZMACH
C     BLAS ZAXPY,ZSWAP,DZASUM,ZCOPY
C     FORTRAN DABS,DMAX1,DBLE,FLOAT,IABS
C
C     INTERNAL VARIABLES
C
      COMPLEX*16 A(15,15),AINV(15,15),ASAVE(15,15),AJK,AJKP1
      COMPLEX*16 AP(120),B(15),DET(2),DETP(2),X(15),XEXACT(15)
      COMPLEX*16 XP(15),T,Z(15)
      DOUBLE PRECISION ANORM,AINORM,COND,COND1
      DOUBLE PRECISION EN,ENORM,EPS,FNORM,ONEPX,Q(6),QS(6),RCOND
      DOUBLE PRECISION RCONDP,RNORM,S,DZASUM,ZMACH,XNORM
      INTEGER I,IQ(6),J,JB
      INTEGER K,KASE,KSING,KM1,KOUNT,KPFAIL,KPVT(15),KPVTP(15)
      INTEGER KS,KSTEP,KSUSP(6),LDA,LUNIT,N,NPRINT
      LOGICAL KPF
C
      COMPLEX*16 ZDUM
      DOUBLE PRECISION CABS1
      DOUBLE PRECISION DREAL,DIMAG
      COMPLEX*16 ZDUMR,ZDUMI
      DREAL(ZDUMR) = ZDUMR
      DIMAG(ZDUMI) = (0.0D0,-1.0D0)*ZDUMI
      CABS1(ZDUM) = DABS(DREAL(ZDUM)) + DABS(DIMAG(ZDUM))
      LDA = 15
C
C     WRITE MATRIX AND SOLUTIONS IF  N .LE. NPRINT
C
      NPRINT = 3
C
      WRITE (LUNIT,540)
      WRITE (LUNIT,940)
C
      DO 10 I = 1, 6
         KSUSP(I) = 0
   10 CONTINUE
      KSING = 0
      KPFAIL = 0
C
C     SET EPS TO ROUNDING UNIT FOR DOUBLE PRECISION ARITHMETIC
C
      EPS = ZMACH(1)
      WRITE (LUNIT,550) EPS
      WRITE (LUNIT,530)
C
C     START MAIN LOOP
C
      KASE = 1
   20 CONTINUE
C
C        GENERATE TEST MATRIX
C
         CALL ZSIXX(A,LDA,N,KASE,LUNIT)
C
C        N = 0 SIGNALS NO MORE TEST MATRICES
C
C     ...EXIT
         IF (N .LE. 0) GO TO 520
         ANORM = 0.0D0
         DO 30 J = 1, N
            ANORM = DMAX1(ANORM,DZASUM(N,A(1,J),1))
   30    CONTINUE
         WRITE (LUNIT,700) ANORM
C
         IF (N .GT. NPRINT) GO TO 50
            WRITE (LUNIT,530)
            DO 40 I = 1, N
               WRITE (LUNIT,740) (A(I,J), J = 1, N)
   40       CONTINUE
            WRITE (LUNIT,530)
   50    CONTINUE
C
C        GENERATE EXACT SOLUTION
C
         XEXACT(1) = (1.0D0,0.0D0)
         IF (N .GE. 2) XEXACT(2) = (0.0D0,1.0D0)
         IF (N .LE. 2) GO TO 70
            DO 60 I = 3, N
               XEXACT(I) = -XEXACT(I-2)
   60       CONTINUE
   70    CONTINUE
C
C        SAVE MATRIX AND GENERATE R.H.S.
C
         DO 90 I = 1, N
            B(I) = (0.0D0,0.0D0)
            DO 80 J = 1, N
               ASAVE(I,J) = A(I,J)
               B(I) = B(I) + A(I,J)*XEXACT(J)
   80       CONTINUE
            X(I) = B(I)
            XP(I) = X(I)
   90    CONTINUE
C
C        FACTOR AND ESTIMATE CONDITION
C
         RCOND = -1.0D0
         CALL ZSICO(A,LDA,N,KPVT,RCOND,Z)
         WRITE (LUNIT,810) (KPVT(I), I = 1, N)
C
C        OUTPUT NULL VECTOR
C
         IF (N .GT. NPRINT) GO TO 110
            WRITE (LUNIT,760)
            DO 100 I = 1, N
               WRITE (LUNIT,770) Z(I)
  100       CONTINUE
            WRITE (LUNIT,530)
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
         CALL ZSPCO(AP,N,KPVTP,RCONDP,Z)
         IF (RCONDP .EQ. RCOND) GO TO 140
            WRITE (LUNIT,830)
            WRITE (LUNIT,870) RCOND,RCONDP
            KPF = .TRUE.
  140    CONTINUE
         K = 0
         KOUNT = 0
         DO 160 J = 1, N
            DO 150 I = 1, J
               K = K + 1
               IF (CABS1(AP(K)-A(I,J)) .NE. 0.0D0) KOUNT = KOUNT + 1
  150       CONTINUE
  160    CONTINUE
         IF (KOUNT .EQ. 0) GO TO 170
            WRITE (LUNIT,830)
            WRITE (LUNIT,880) KOUNT
            KPF = .TRUE.
  170    CONTINUE
C
C        TEST FOR SINGULARITY
C
         IF (RCOND .GT. 0.0D0) GO TO 180
            WRITE (LUNIT,750) RCOND
            KSING = KSING + 1
            CALL ZSIDI(A,LDA,N,KPVT,DET,Z,10)
            WRITE (LUNIT,790) DET(1)
         GO TO 510
  180    CONTINUE
            COND = 1.0D0/RCOND
            WRITE (LUNIT,570) COND
            ONEPX = 1.0D0 + RCOND
            IF (ONEPX .EQ. 1.0D0) WRITE (LUNIT,560)
C
C           COMPUTE INVERSE AND DETERMINANT TRUE CONDITION
C
            DO 200 J = 1, N
               DO 190 I = 1, J
                  AINV(I,J) = A(I,J)
  190          CONTINUE
  200       CONTINUE
            CALL ZSIDI(AINV,LDA,N,KPVT,DET,Z,11)
            AINORM = 0.0D0
            DO 220 J = 1, N
               DO 210 I = J, N
                  AINV(I,J) = AINV(J,I)
  210          CONTINUE
               AINORM = DMAX1(AINORM,DZASUM(N,AINV(1,J),1))
  220       CONTINUE
            COND1 = ANORM*AINORM
            WRITE (LUNIT,580) COND1
            WRITE (LUNIT,790) DET(1)
            WRITE (LUNIT,800) DET(2)
C
C           SOLVE A*X = B
C
            CALL ZSISL(A,LDA,N,KPVT,X)
C
            IF (N .GT. NPRINT) GO TO 240
               WRITE (LUNIT,590)
               DO 230 I = 1, N
                  WRITE (LUNIT,780) X(I)
  230          CONTINUE
               WRITE (LUNIT,530)
  240       CONTINUE
C
C           MORE PACKED COMPARE
C
            CALL ZSPSL(AP,N,KPVTP,XP)
            KOUNT = 0
            DO 250 I = 1, N
               IF (CABS1(XP(I)-X(I)) .NE. 0.0D0) KOUNT = KOUNT + 1
  250       CONTINUE
            IF (KOUNT .EQ. 0) GO TO 260
               WRITE (LUNIT,830)
               WRITE (LUNIT,890) KOUNT
               KPF = .TRUE.
  260       CONTINUE
            CALL ZSPDI(AP,N,KPVTP,DETP,Z,11)
            IF (CABS1(DETP(1)-DET(1)) .EQ. 0.0D0
     *          .AND. CABS1(DETP(2)-DET(2)) .EQ. 0.0D0) GO TO 270
               WRITE (LUNIT,830)
               WRITE (LUNIT,900) DETP
               KPF = .TRUE.
  270       CONTINUE
            KOUNT = 0
            K = 0
            DO 290 J = 1, N
               DO 280 I = 1, J
                  K = K + 1
                  IF (CABS1(AP(K)-AINV(I,J)) .NE. 0.0D0)
     *               KOUNT = KOUNT + 1
  280          CONTINUE
  290       CONTINUE
            IF (KOUNT .EQ. 0) GO TO 300
               WRITE (LUNIT,830)
               WRITE (LUNIT,910) KOUNT
               KPF = .TRUE.
  300       CONTINUE
C
C           RECONSTRUCT  A  FROM FACTORS.
C
            K = 1
  310       IF (K .GT. N) GO TO 400
               KM1 = K - 1
               IF (KPVT(K) .LT. 0) GO TO 340
C
C                 1 BY 1
C
                  IF (KM1 .LT. 1) GO TO 330
                  DO 320 JB = 1, KM1
                     J = K - JB
                     AJK = A(K,K)*A(J,K)
                     T = AJK
                     CALL ZAXPY(J,T,A(1,K),1,A(1,J),1)
                     A(J,K) = -AJK
  320             CONTINUE
  330             CONTINUE
                  KSTEP = 1
               GO TO 370
  340          CONTINUE
C
C                 2 BY 2
C
                  IF (KM1 .LT. 1) GO TO 360
                  DO 350 JB = 1, KM1
                     J = K - JB
                     AJK = A(K,K)*A(J,K) + A(K,K+1)*A(J,K+1)
                     AJKP1 = A(K,K+1)*A(J,K) + A(K+1,K+1)*A(J,K+1)
                     T = AJK
                     CALL ZAXPY(J,T,A(1,K),1,A(1,J),1)
                     T = AJKP1
                     CALL ZAXPY(J,T,A(1,K+1),1,A(1,J),1)
                     A(J,K) = -AJK
                     A(J,K+1) = -AJKP1
  350             CONTINUE
  360             CONTINUE
                  KSTEP = 2
  370          CONTINUE
C
C              SWAP
C
               KS = IABS(KPVT(K))
               IF (KS .EQ. K) GO TO 390
                  CALL ZSWAP(KS,A(1,KS),1,A(1,K),1)
                  DO 380 JB = KS, K
                     J = K + KS - JB
                     T = A(J,K)
                     A(J,K) = A(KS,J)
                     A(KS,J) = T
  380             CONTINUE
                  IF (KSTEP .EQ. 2)
     *               CALL ZSWAP(1,A(KS,K+1),1,A(K,K+1),1)
  390          CONTINUE
               K = K + KSTEP
            GO TO 310
  400       CONTINUE
            DO 420 J = 1, N
               DO 410 I = J, N
                  A(I,J) = A(J,I)
  410          CONTINUE
  420       CONTINUE
C
C           COMPUTE ERRORS AND RESIDUALS
C              E  =  X - XEXACT
C              R  =  B - A*X
C              F  =  A - U*D*TRANS(U)
C
            XNORM = DZASUM(N,X,1)
            ENORM = 0.0D0
            FNORM = 0.0D0
            DO 440 J = 1, N
               ENORM = ENORM + CABS1(X(J)-XEXACT(J))
               T = -X(J)
               CALL ZAXPY(N,T,ASAVE(1,J),1,B,1)
               S = 0.0D0
               DO 430 I = 1, J
                  S = S + CABS1(ASAVE(I,J)-A(I,J))
  430          CONTINUE
               IF (S .GT. FNORM) FNORM = S
  440       CONTINUE
            RNORM = DZASUM(N,B,1)
C
C           A*INV(A) - I
C
            AINORM = 0.0D0
            DO 470 J = 1, N
               DO 450 I = 1, N
                  B(I) = (0.0D0,0.0D0)
  450          CONTINUE
               DO 460 K = 1, N
                  T = AINV(K,J)
                  CALL ZAXPY(N,T,ASAVE(1,K),1,B,1)
  460          CONTINUE
               B(J) = B(J) - (1.0D0,0.0D0)
               AINORM = DMAX1(AINORM,DZASUM(N,B,1))
  470       CONTINUE
C
            WRITE (LUNIT,600) ENORM
            WRITE (LUNIT,610) RNORM
            WRITE (LUNIT,710) FNORM
            WRITE (LUNIT,720) AINORM
C
C           COMPUTE TEST RATIOS
C
            EN = DBLE(FLOAT(N))
            Q(1) = COND/COND1
            Q(2) = COND1/COND
            Q(3) = ENORM/(EPS*COND*XNORM)
            Q(4) = RNORM/(EPS*ANORM*XNORM)
            Q(5) = FNORM/(EPS*ANORM)
            Q(6) = AINORM/(EPS*COND)
            WRITE (LUNIT,530)
            WRITE (LUNIT,620)
            WRITE (LUNIT,530)
            WRITE (LUNIT,670)
            WRITE (LUNIT,680)
            WRITE (LUNIT,690)
            WRITE (LUNIT,530)
            WRITE (LUNIT,730) (Q(I), I = 1, 6)
            WRITE (LUNIT,530)
C
C           LOOK FOR SUSPICIOUS RATIOS
C
            QS(1) = 1.0D0 + 4.0D0*EPS
            QS(2) = 10.0D0
            IF (N .EQ. 1) EN = 2.0D0
            DO 480 I = 3, 6
               QS(I) = EN
  480       CONTINUE
            KOUNT = 0
            DO 500 I = 1, 6
               IQ(I) = 0
               IF (Q(I) .LE. QS(I)) GO TO 490
                  IQ(I) = 1
                  KSUSP(I) = KSUSP(I) + 1
                  KOUNT = KOUNT + 1
  490          CONTINUE
  500       CONTINUE
            IF (KOUNT .EQ. 0) WRITE (LUNIT,920)
            IF (KOUNT .NE. 0) WRITE (LUNIT,930) (IQ(I), I = 1, 6)
            WRITE (LUNIT,530)
  510    CONTINUE
C
         IF (.NOT.KPF) WRITE (LUNIT,820)
         IF (KPF) KPFAIL = KPFAIL + 1
         WRITE (LUNIT,630)
         KASE = KASE + 1
      GO TO 20
  520 CONTINUE
C
C     FINISH MAIN LOOP
C
C     SUMMARY
C
      WRITE (LUNIT,640)
      KASE = KASE - 1
      WRITE (LUNIT,650) KASE
      WRITE (LUNIT,850) KSING
      WRITE (LUNIT,840) KPFAIL
      WRITE (LUNIT,660) KSUSP
      WRITE (LUNIT,860)
      RETURN
C
C     MOST FORMATS, ALSO SOME IN ZSIXX
C
  530 FORMAT (1H )
  540 FORMAT (29H1LINPACK TESTER, ZSI**, ZSP**)
  550 FORMAT ( / 14H EPSILON     =, 1PD13.5)
  560 FORMAT ( / 16H MAYBE SINGULAR. /)
  570 FORMAT (14H COND        =, 1PD13.5)
  580 FORMAT (14H ACTUAL COND =, 1PD13.5)
  590 FORMAT ( / 4H X =)
  600 FORMAT (14H ERROR NORM  =, 1P2D13.5)
  610 FORMAT (14H RESID NORM  =, 1P2D13.5)
  620 FORMAT (26H TEST RATIOS.. E = EPSILON)
  630 FORMAT ( / 14H ************* /)
  640 FORMAT (8H1SUMMARY)
  650 FORMAT (18H NUMBER OF TESTS =, I4)
  660 FORMAT ( / 30H NUMBER OF SUSPICIOUS RATIOS =, 6I4)
  670 FORMAT (30H     COND     ACTUAL    ERROR ,
     *        30H    RESID   A-U*D*UT  A*AI - I)
  680 FORMAT (6(10H   -------))
  690 FORMAT (30H    ACTUAL     COND   E*COND*X,
     *        30H    E*A*X      E*A     E*COND )
  700 FORMAT (14H NORM(A)     =, 1PD13.5)
  710 FORMAT (14H NORM(A-UDUT)=, 1PD13.5)
  720 FORMAT (14H NORM(A*AI-I)=, 1PD13.5)
  730 FORMAT (6(1X, F9.4))
  740 FORMAT (1H , 6G11.4)
  750 FORMAT (14H 1/COND      =, 1PD13.5)
  760 FORMAT ( / 7H NULL =)
  770 FORMAT (2G14.6)
  780 FORMAT (2G14.6)
  790 FORMAT (14H DET FRACT   =, 2F9.5)
  800 FORMAT (14H DET EXPON   =, 2F9.0)
  810 FORMAT (14H KPVT        =, 15I3)
  820 FORMAT ( / 22H PACKED ROUTINES AGREE)
  830 FORMAT ( / 30H PACKED ROUTINES DO NOT AGREE,)
  840 FORMAT (28H NUMBER OF PACKED FAILURES =, I4)
  850 FORMAT (23H NUMBER OF ZERO PIVOT =, I4)
  860 FORMAT ( / 12H END OF TEST)
  870 FORMAT (8H RCOND =, 1P2D13.5)
  880 FORMAT (12H KOUNT(FA) =, I4)
  890 FORMAT (12H KOUNT(SL) =, I4)
  900 FORMAT (8H DET   =, F9.5, F9.0)
  910 FORMAT (12H KOUNT(DI) =, I4)
  920 FORMAT (21H NO SUSPICIOUS RATIOS)
  930 FORMAT (I8, 5I10 / 7X, 28H1 INDICATES SUSPICIOUS RATIO)
  940 FORMAT (29H THIS VERSION DATED 08/14/78.)
      END
      SUBROUTINE ZSIXX(A,LDA,N,KASE,LUNIT)
      INTEGER LDA,N,KASE,LUNIT
      COMPLEX*16 A(LDA,1)
C
C     GENERATES COMPLEX*16 SYMMETRIC INDEFINITE TEST MATRICES
C
C     EXTERNAL ZMACH
C     FORTRAN DABS,CDABS,DCMPLX,DBLE,FLOAT,IABS
      COMPLEX*16 T
      DOUBLE PRECISION TINY,HUGE,ZMACH
      INTEGER I,J
      DOUBLE PRECISION CABS1
      COMPLEX*16 ZDUM
      DOUBLE PRECISION DREAL,DIMAG
      COMPLEX*16 ZDUMR,ZDUMI
      DREAL(ZDUMR) = ZDUMR
      DIMAG(ZDUMI) = (0.0D0,-1.0D0)*ZDUMI
      CABS1(ZDUM) = DABS(DREAL(ZDUM)) + DABS(DIMAG(ZDUM))
C
      GO TO (10, 10, 10, 50, 50, 70, 70, 70, 110, 150, 190, 230, 250,
     *       290, 330, 330, 380, 430, 480), KASE
C
C     KASE 1, 2 AND 3
C
   10 CONTINUE
         N = 5*KASE
         WRITE (LUNIT,20) KASE,N
   20    FORMAT (5H KASE, I3, 3X, 19HCOMPLEX*16 HILBERT  / 4H N =, I4)
         T = (1.0D0,1.0D0)
         T = T/CDABS(T)
         DO 40 J = 1, N
            DO 30 I = 1, J
               A(I,J) = T**(I - J)/DCMPLX(DBLE(FLOAT(I+J-1)),0.0D0)
               A(J,I) = A(I,J)
C              FOR NON-C0MPLEX MATRICES, A(I,J) = 1.0/DBLEFLOAT(I+J-1)
   30       CONTINUE
   40    CONTINUE
      GO TO 490
C
C     KASE 4 AND 5
C
   50 CONTINUE
         N = 1
         WRITE (LUNIT,60) KASE,N
   60    FORMAT (5H KASE, I3, 3X, 16HMONOELEMENTAL    / 4H N =, I4)
         IF (KASE .EQ. 4) A(1,1) = (3.0D0,0.0D0)
         IF (KASE .EQ. 5) A(1,1) = (0.0D0,0.0D0)
      GO TO 490
C
C     KASE 6, 7 AND 8
C
   70 CONTINUE
         N = 15
         WRITE (LUNIT,80) KASE,N
   80    FORMAT (5H KASE, I3, 3X, 16HTRIDIAGONAL      / 4H N =, I4)
         T = (1.0D0,0.0D0)
         IF (KASE .EQ. 7) T = (3.0D0,1.0D0)
         IF (KASE .EQ. 8) T = (100.0D0,100.0D0)
         DO 100 J = 1, N
            DO 90 I = 1, J
               A(I,J) = (0.0D0,0.0D0)
               IF (I .EQ. J) A(I,I) = (4.0D0,0.0D0)
               IF (I .EQ. J - 1) A(I,J) = T
               A(J,I) = A(I,J)
   90       CONTINUE
  100    CONTINUE
      GO TO 490
C
C     KASE 9
C
  110 CONTINUE
         N = 5
         WRITE (LUNIT,120) KASE,N
  120    FORMAT (5H KASE, I3, 3X, 16HPENTADIAGONAL    / 4H N =, I4)
         DO 140 J = 1, N
            DO 130 I = 1, J
               A(I,J) = (0.0D0,0.0D0)
               IF (I .GE. J - 2)
     *            A(I,J) = DCMPLX((5.0D0-DBLE(FLOAT(IABS(I-J))))
     *                            **(10-I-J),0.0D0)
               A(J,I) = A(I,J)
  130       CONTINUE
  140    CONTINUE
      GO TO 490
C
C     KASE 10
C
  150 CONTINUE
         N = 6
         WRITE (LUNIT,160) KASE,N
  160    FORMAT (5H KASE, I3, 3X, 16HTRIDIAG INVERSE  / 4H N =, I4)
         DO 180 J = 1, N
            DO 170 I = 1, J
               A(I,J) = DCMPLX(DBLE(FLOAT(N+1-J)),DBLE(FLOAT(N+1-J)))
               A(J,I) = A(I,J)
  170       CONTINUE
  180    CONTINUE
      GO TO 490
C
C     KASE 11
C
  190 CONTINUE
         N = 10
         WRITE (LUNIT,200) KASE,N
  200    FORMAT (5H KASE, I3, 3X, 16HZERO DIAGONAL    / 4H N =, I4)
         DO 220 J = 1, N
            DO 210 I = 1, J
               A(I,J) = (0.0D0,0.0D0)
               IF (I .EQ. J - 1) A(I,J) = (1.0D0,1.0D0)
               A(J,I) = A(I,J)
  210       CONTINUE
  220    CONTINUE
      GO TO 490
C
C     KASE 12
C
  230 CONTINUE
         N = 2
         WRITE (LUNIT,240) KASE,N
  240    FORMAT (5H KASE, I3, 3X, 16HTWO BY TWO       / 4H N =, I4)
         A(1,1) = (4.0D0,0.0D0)
         A(1,2) = (1.0D0,2.0D0)
         A(2,1) = A(1,2)
         A(2,2) = (0.0D0,0.0D0)
      GO TO 490
C
C     KASE 13
C
  250 CONTINUE
         N = 6
         WRITE (LUNIT,260) KASE,N
  260    FORMAT (5H KASE, I3, 3X, 16H ZERO MATRIX     / 4H N =, I4)
         DO 280 J = 1, N
            DO 270 I = 1, J
               A(I,J) = (0.0D0,0.0D0)
               A(J,I) = A(I,J)
  270       CONTINUE
  280    CONTINUE
      GO TO 490
C
C     KASE 14
C
  290 CONTINUE
         N = 3
         WRITE (LUNIT,300) KASE,N
  300    FORMAT (5H KASE, I3 / 4H N =, I4)
         DO 320 I = 1, N
            DO 310 J = 1, N
               A(I,J) = (0.0D0,0.0D0)
  310       CONTINUE
  320    CONTINUE
         A(1,3) = (1.0D0,0.0D0)
         A(3,1) = (1.0D0,0.0D0)
      GO TO 490
C
C     KASE 15 AND 16
C
  330 CONTINUE
         N = 15
         WRITE (LUNIT,340) KASE,N
  340    FORMAT (5H KASE, I3, 3X, 16H                 / 4H N =, I4)
         DO 370 J = 1, N
            DO 360 I = 1, J
               A(I,J) = (-1.0D0,1.0D0)
               A(J,I) = A(I,J)
               IF (I .NE. J) GO TO 350
                  IF (KASE .EQ. 15) A(I,I) = (26.0D0,26.0D0)
                  IF (KASE .EQ. 16)
     *               A(I,I) = DCMPLX(DBLE(FLOAT(I)),0.0D0)
  350          CONTINUE
  360       CONTINUE
  370    CONTINUE
      GO TO 490
C
C     KASE 17
C
  380 CONTINUE
         N = 5
         WRITE (LUNIT,390) KASE,N
  390    FORMAT (5H KASE, I3, 3X, 16HNEAR UNDERFLOW   / 4H N =, I4)
         TINY = ZMACH(2)
         WRITE (LUNIT,400) TINY
  400    FORMAT (14H TINY        =, 1PD13.5)
         DO 420 J = 1, N
            DO 410 I = 1, J
               A(I,J) = TINY
     *                  *DCMPLX(DBLE(FLOAT(IABS(I-J)))/DBLE(FLOAT(I+J))
     *                  ,DBLE(FLOAT(I-J)))
               A(J,I) = A(I,J)
  410       CONTINUE
  420    CONTINUE
      GO TO 490
C
C     KASE 18
C
  430 CONTINUE
         N = 5
         WRITE (LUNIT,440) KASE,N
  440    FORMAT (5H KASE, I3, 3X, 16HNEAR OVERFLOW    / 4H N =, I4)
         HUGE = ZMACH(3)
         WRITE (LUNIT,450) HUGE
  450    FORMAT (14H HUGE        =, 1PD13.5)
         DO 470 J = 1, N
            DO 460 I = 1, J
               A(I,J) = HUGE
     *                  *DCMPLX(DBLE(FLOAT(IABS(I-J)))/DBLE(FLOAT(I+J))
     *                  ,DBLE(FLOAT(I-J)))
               A(J,I) = A(I,J)
  460       CONTINUE
  470    CONTINUE
      GO TO 490
C
  480 CONTINUE
         N = 0
  490 CONTINUE
      RETURN
C
      END
