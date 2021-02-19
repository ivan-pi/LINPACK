C     MAIN PROGRAM
      INTEGER LUNIT
C     ALLOW 5000 UNDERFLOWS.
      CALL TRAPS(0,0,5001,0,0)
C
C     OUTPUT UNIT NUMBER
C
      LUNIT = 6
C
      CALL DQRXX(LUNIT)
      STOP
      END
      SUBROUTINE DQRXX(LUNIT)
      INTEGER LUNIT
      DOUBLE PRECISION RX(20,10),R(10,10),XROW(10),Z(10,4),YROW(10),
     *                 S(10)
      DOUBLE PRECISION SCALE,TINY,HUGE,RHO(2),C(10),DMACH
      INTEGER N,P,LDX,LDR,LDZ,CASE
      LOGICAL NOTWRT,OFLOW,UFLOW
      COMMON SCALE,NOTWRT,OFLOW,UFLOW
      LDZ = 10
      LDX = 20
      LDR = 10
      NOTWRT = .TRUE.
      OFLOW = .FALSE.
      UFLOW = .FALSE.
      TINY = DMACH(2)
      HUGE = DMACH(3)
      SCALE = 1.0D0
      DO 60 CASE = 1, 3
         GO TO (10, 20, 30), CASE
   10    CONTINUE
            N = 20
            P = 10
         GO TO 40
   20    CONTINUE
            N = 10
            P = 4
         GO TO 40
   30    CONTINUE
            N = 10
            P = 1
   40    CONTINUE
         CALL DGETRX(RX,LDX,N,P,S)
         WRITE (LUNIT,130) CASE,N,P
         IF (NOTWRT) GO TO 50
            WRITE (LUNIT,160)
            CALL DARRAY(RX,LDX,N,P,P,LUNIT)
   50    CONTINUE
         CALL DDRUD1(R,LDR,RX,LDX,N,P,XROW,C,S,LUNIT)
         CALL DDRUD2(R,LDR,RX,LDX,P,XROW,YROW,Z,LDZ,RHO,C,S,LUNIT)
         CALL DDRDD(R,LDR,RX,LDX,P,Z,LDZ,XROW,YROW,RHO,C,S,LUNIT)
   60 CONTINUE
      CASE = 4
      N = 10
      P = 4
      WRITE (LUNIT,130) CASE,N,P
      WRITE (LUNIT,140)
      CALL DGETRX(RX,LDX,N,P,S)
      DO 80 J = 1, P
         DO 70 I = 1, N
            RX(I,J) = HUGE*RX(I,J)
   70    CONTINUE
   80 CONTINUE
      SCALE = HUGE
      OFLOW = .TRUE.
      IF (NOTWRT) GO TO 90
         WRITE (LUNIT,160)
         CALL DARRAY(RX,LDX,N,P,P,LUNIT)
   90 CONTINUE
      CALL DDRUD1(R,LDR,RX,LDX,N,P,XROW,C,S,LUNIT)
      CALL DDRUD2(R,LDR,RX,LDX,P,XROW,YROW,Z,LDZ,RHO,C,S,LUNIT)
      CALL DDRDD(R,LDR,RX,LDX,P,Z,LDZ,XROW,YROW,RHO,C,S,LUNIT)
      OFLOW = .FALSE.
      CASE = 5
      N = 10
      P = 4
      WRITE (LUNIT,130) CASE,N,P
      WRITE (LUNIT,150)
      CALL DGETRX(RX,LDX,N,P,S)
      DO 110 J = 1, P
         DO 100 I = 1, N
            RX(I,J) = TINY*RX(I,J)
  100    CONTINUE
  110 CONTINUE
      SCALE = TINY
      UFLOW = .TRUE.
      IF (NOTWRT) GO TO 120
         WRITE (LUNIT,160)
         CALL DARRAY(RX,LDX,N,P,P,LUNIT)
  120 CONTINUE
      CALL DDRUD1(R,LDR,RX,LDX,N,P,XROW,C,S,LUNIT)
      CALL DDRUD2(R,LDR,RX,LDX,P,XROW,YROW,Z,LDZ,RHO,C,S,LUNIT)
      CALL DDRDD(R,LDR,RX,LDX,P,Z,LDZ,XROW,YROW,RHO,C,S,LUNIT)
      RETURN
C
  130 FORMAT (11H1    CASE =, I2, 5X, 3HN =, I2, 5X, 3HP =, I2 /////)
  140 FORMAT (22H         OVERFLOW TEST /////)
  150 FORMAT (24H          UNDERFLOW TEST /////)
  160 FORMAT (5H   RX)
      END
      SUBROUTINE DARRAY(A,LDA,M,N,NNL,LUNIT)
      INTEGER LDA,M,N,NNL,LUNIT
      DOUBLE PRECISION A(LDA,1)
C     X
C          FORTRAN IABS,MIN0
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
C
   70 FORMAT (1H , 4(2D13.6, 4X))
      END
      SUBROUTINE DDRDD(R,LDR,RX,LDX,P,Z,LDZ,XROW,YROW,RHO,C,S,LUNIT)
      INTEGER LDR,LDX,LDZ,P,LUNIT
      DOUBLE PRECISION R(LDR,1),RX(LDX,1),Z(LDZ,1),XROW(1),YROW(1),S(1)
      DOUBLE PRECISION RHO(1),SCALE,C(1)
      LOGICAL NOTWRT,OFLOW,UFLOW
      COMMON SCALE,NOTWRT,OFLOW,UFLOW
      INTEGER I,J,JP1,PM1
      WRITE (LUNIT,50)
      CALL DCHDD(R,LDR,P,XROW,Z,LDZ,2,YROW,RHO,C,S,INFO)
      PM1 = P - 1
      IF (PM1 .LT. 1) GO TO 30
      DO 20 J = 1, PM1
         JP1 = J + 1
         DO 10 I = JP1, P
            R(I,J) = 0.0D0
   10    CONTINUE
   20 CONTINUE
   30 CONTINUE
      IF (NOTWRT) GO TO 40
         WRITE (LUNIT,80)
         CALL DARRAY(R,LDR,P,P,P,LUNIT)
         WRITE (LUNIT,60)
         CALL DARRAY(Z,LDZ,P,2,2,LUNIT)
         WRITE (LUNIT,70) (RHO(I), I = 1, 2)
   40 CONTINUE
      CALL DMPDD(R,LDR,P,RX,LDX,Z,Z(1,3),LDZ,RHO,LUNIT)
      RETURN
   50 FORMAT ( /////
     *         46H     STEP THREE    DOWNDATING XROW,YROW AND Z, ///)
   60 FORMAT ( /// 4H   Z)
   70 FORMAT ( /// 6H   RHO // 1H , 2D13.6)
   80 FORMAT (4H   R)
      END
      SUBROUTINE DDRUD1(R,LDR,RX,LDX,N,P,XROW,C,S,LUNIT)
      INTEGER N,P,LDR,LDX,LUNIT
      DOUBLE PRECISION R(LDR,1),RX(LDX,1),XROW(1),S(1)
      DOUBLE PRECISION C(1)
      LOGICAL NOTWRT,OFLOW
      DOUBLE PRECISION RELM,XELM,Y,Z
      DOUBLE PRECISION XMAX,RMAX,TEST,T,SCALE,DMACH
      COMMON SCALE,NOTWRT,OFLOW,UFLOW
      DO 20 I = 1, P
         DO 10 J = 1, P
            R(I,J) = 0.0D0
   10    CONTINUE
   20 CONTINUE
      DO 40 I = 1, N
         DO 30 J = 1, P
            XROW(J) = RX(I,J)
   30    CONTINUE
         CALL DCHUD(R,LDR,P,XROW,Z,1,0,Y,Y,C,S)
   40 CONTINUE
      WRITE (LUNIT,100)
      IF (NOTWRT) GO TO 50
         WRITE (LUNIT,120)
         CALL DARRAY(R,LDR,P,P,P,LUNIT)
   50 CONTINUE
      RMAX = 0.0D0
      XMAX = 0.0D0
      DO 90 I = 1, P
         DO 80 J = 1, I
            RELM = 0.0D0
            XELM = 0.0D0
            NIM = MIN0(I,J)
            DO 60 K = 1, NIM
               RELM = RELM + R(K,I)/SCALE*(R(K,J)/SCALE)
   60       CONTINUE
            DO 70 K = 1, N
               XELM = XELM + RX(K,I)/SCALE*(RX(K,J)/SCALE)
   70       CONTINUE
            T = DMAX1(DABS(XELM),0.0D0)
            XMAX = DMAX1(XMAX,T)
            T = DMAX1(DABS(XELM-RELM),0.0D0)
            RMAX = DMAX1(RMAX,T)
   80    CONTINUE
   90 CONTINUE
      TEST = RMAX/XMAX/DMACH(1)
      WRITE (LUNIT,110) TEST
      RETURN
C
  100 FORMAT ( ///// 25H    STEP ONE   UPDATING X ///)
  110 FORMAT ( /// 15H     STATISTICS //
     *         42H      RH*R    ............................, D10.2)
  120 FORMAT (4H   R)
      END
      SUBROUTINE DDRUD2(R,LDR,RX,LDX,P,XROW,YROW,Z,LDZ,RHO,C,S,LUNIT)
      INTEGER LDR,LDX,LDZ,P,LUNIT
      DOUBLE PRECISION R(LDR,1),RX(LDX,1),XROW(1),YROW(1),Z(LDZ,1),S(1)
      DOUBLE PRECISION RHO(1),C(1)
      DOUBLE PRECISION SCALE
      LOGICAL NOTWRT,OFLOW,UFLOW
      COMMON SCALE,NOTWRT,OFLOW,UFLOW
      DO 20 I = 1, P
         DO 10 J = 1, P
            RX(I,J) = R(I,J)
   10    CONTINUE
   20 CONTINUE
      DO 40 I = 1, P
         Z(I,1) = 0.0D0
         DO 30 J = I, P
            Z(I,1) = Z(I,1) + R(I,J)
   30    CONTINUE
         Z(I,2) = Z(I,1) + SCALE
         Z(I,3) = Z(I,1)
         Z(I,4) = Z(I,2)
   40 CONTINUE
      DO 60 I = 1, P
         XROW(I) = 0.0D0
         DO 50 J = 1, I
            XROW(I) = XROW(I) + R(J,I)
   50    CONTINUE
   60 CONTINUE
      YROW(1) = 0.0D0
      DO 70 I = 1, P
         YROW(1) = YROW(1) + XROW(I)
   70 CONTINUE
      YROW(2) = YROW(1) - SCALE
      RHO(1) = SCALE
      RHO(2) = SCALE
      IF (NOTWRT) GO TO 80
         WRITE (LUNIT,100)
         CALL DARRAY(XROW,P,P,1,-P,LUNIT)
         WRITE (LUNIT,110)
         CALL DARRAY(Z,LDZ,P,2,2,LUNIT)
         WRITE (LUNIT,120)
         CALL DARRAY(YROW,2,2,1,-2,LUNIT)
   80 CONTINUE
      CALL DCHUD(R,LDR,P,XROW,Z,LDZ,2,YROW,RHO,C,S)
      WRITE (LUNIT,130)
      IF (NOTWRT) GO TO 90
         WRITE (LUNIT,150)
         CALL DARRAY(R,LDR,P,P,P,LUNIT)
         WRITE (LUNIT,110)
         CALL DARRAY(Z,LDZ,P,2,2,LUNIT)
         WRITE (LUNIT,140) (RHO(I), I = 1, 2)
   90 CONTINUE
      CALL DMPUD(R,LDR,RX,LDX,P,XROW,Z,LDZ,RHO,LUNIT)
      RETURN
C
  100 FORMAT ( ///// 7H   XROW)
  110 FORMAT ( /// 4H   Z)
  120 FORMAT ( /// 7H   YROW)
  130 FORMAT ( ///// 41H     STEP TWO   UPDATING XROW, YROW AND Z ///)
  140 FORMAT ( /// 6H   RHO // 1H , 2D13.6)
  150 FORMAT (4H   R)
      END
      SUBROUTINE DGETRX(X,LDX,N,P,S)
      INTEGER N,P,LDX
      DOUBLE PRECISION X(LDX,1),S(1)
      INTEGER PD2,PD2D1,I
      PD2 = MAX0(P/2,1)
      PD2D1 = PD2 + 1
      DO 10 I = 1, PD2
         S(I) = 1.0D0
   10 CONTINUE
      IF (P .LT. PD2D1) GO TO 30
      DO 20 I = PD2D1, P
         S(I) = 0.5D0
   20 CONTINUE
   30 CONTINUE
      CALL DXGEN(X,LDX,N,P,S)
      RETURN
      END
      SUBROUTINE DMPDD(R,LDR,P,ROLD,LDRD,Z,ZUP,LDZ,RHO,LUNIT)
      INTEGER LDR,P,LDRD,LDZ,LUNIT
      DOUBLE PRECISION R(LDR,1),ROLD(LDRD,1),Z(LDZ,1),ZUP(LDZ,1)
      DOUBLE PRECISION RHO(1),DMACH
      DOUBLE PRECISION T
      INTEGER I,J
      DOUBLE PRECISION TR,TZ(2),TRHO(2),MAXOLD,MAXREL,SCALE
      LOGICAL NOTWRT,OFLOW,UFLOW
      COMMON SCALE,NOTWRT,OFLOW,UFLOW
      MAXOLD = 0.0D0
      MAXREL = 0.0D0
      DO 20 J = 1, P
         DO 10 I = 1, J
            T = ROLD(I,J)
            MAXOLD = DMAX1(MAXOLD,DMAX1(DABS(T),0.0D0))
            T = ROLD(I,J) - R(I,J)
            MAXREL = DMAX1(MAXREL,DMAX1(DABS(T),0.0D0))
   10    CONTINUE
   20 CONTINUE
      TR = MAXREL/MAXOLD/DMACH(1)
      DO 40 I = 1, 2
         MAXOLD = 0.0D0
         MAXREL = 0.0D0
         DO 30 J = 1, P
            T = ZUP(J,I)
            MAXOLD = DMAX1(MAXOLD,DMAX1(DABS(T),0.0D0))
            T = ZUP(J,I) - Z(J,I)
            MAXREL = DMAX1(MAXREL,DMAX1(DABS(T),0.0D0))
   30    CONTINUE
         TZ(I) = MAXREL/MAXOLD/DMACH(1)
         TRHO(I) = DABS(RHO(I)/SCALE-1.0D0)/DMACH(1)
   40 CONTINUE
      WRITE (LUNIT,50) TR,TZ,TRHO
      RETURN
C
   50 FORMAT ( /// 25H     STATSTICS STEP THREE //
     *         33H        R   ....................., D10.2 /
     *         33H        Z(1)   .................., D10.2 /
     *         33H        Z(2)   .................., D10.2 /
     *         33H        RHO(1)   ................, D10.2 /
     *         33H        RHO(2)   ................, D10.2)
      END
      SUBROUTINE DMPUD(R,LDR,RX,LDX,P,XROW,Z,LDZ,RHO,LUNIT)
      INTEGER LDR,LDX,P,LDZ
      DOUBLE PRECISION R(LDR,1),RX(LDX,1),XROW(1),Z(LDZ,1)
      DOUBLE PRECISION RHO(1)
      DOUBLE PRECISION RELM,XELM,TMP1(10),TMP
      DOUBLE PRECISION TEST1,TEST2(2),TEST3,TEST4,MAXR,MAXX,T
      LOGICAL NOTWRT,OFLOW,UFLOW
      DOUBLE PRECISION SCALE,DMACH
      COMMON SCALE,NOTWRT,OFLOW,UFLOW
      MAXR = 0.0D0
      MAXX = 0.0D0
      DO 30 I = 1, P
         DO 20 J = 1, I
            RELM = 0.0D0
            XELM = 0.0D0
            NIM = MIN0(I,J)
            DO 10 K = 1, NIM
               RELM = RELM + R(K,I)/SCALE*(R(K,J)/SCALE)
               XELM = XELM + RX(K,I)/SCALE*(RX(K,J)/SCALE)
   10       CONTINUE
            XELM = XELM + XROW(I)/SCALE*(XROW(J)/SCALE)
            T = DMAX1(DABS(XELM),0.0D0)
            MAXX = DMAX1(MAXX,T)
            T = DMAX1(DABS(XELM-RELM),0.0D0)
            MAXR = DMAX1(MAXR,T)
   20    CONTINUE
   30 CONTINUE
      TEST1 = MAXR/MAXX/DMACH(1)
      TMP = 0.0D0
      DO 50 I = 1, P
         TMP1(I) = 0.0D0
         DO 40 J = I, P
            TMP1(I) = TMP1(I) + RX(I,J)
   40    CONTINUE
         TMP = TMP + XROW(I)
   50 CONTINUE
      DO 80 K = 1, 2
         MAXR = 0.0D0
         MAXX = 0.0D0
         DO 70 I = 1, P
            RELM = 0.0D0
            XELM = 0.0D0
            DO 60 J = 1, I
               RELM = RELM + R(J,I)/SCALE*(Z(J,K)/SCALE)
               XELM = XELM + RX(J,I)/SCALE*(TMP1(J)/SCALE)
   60       CONTINUE
            XELM = XELM + XROW(I)/SCALE*(TMP/SCALE)
            T = DMAX1(DABS(XELM-RELM),0.0D0)
            MAXR = DMAX1(MAXR,T)
            T = DMAX1(DABS(XELM),0.0D0)
            MAXX = DMAX1(MAXX,T)
   70    CONTINUE
         TEST2(K) = MAXR/MAXX/DMACH(1)
   80 CONTINUE
      TEST3 = DABS(RHO(1)/SCALE-1.0D0)/DMACH(1)
      T = DSQRT(2.0D0+DBLE(FLOAT(P)))
      TEST4 = DABS(RHO(2)/SCALE-T)/T/DMACH(1)
      WRITE (LUNIT,90) TEST1,TEST2,TEST3,TEST4
      RETURN
   90 FORMAT ( /// 24H     STATISTICS STEP TWO //
     *         31H        RH*R   ................, D10.2 /
     *         31H        RH*Z(1)   ............., D10.2 /
     *         31H        RH*Z(2)   ............., D10.2 /
     *         31H        RHO(1)   .............., D10.2 /
     *         31H        RHO(2)   .............., D10.2)
      END
      SUBROUTINE DXGEN(X,LDX,N,P,S)
      INTEGER LDX,N,P
      DOUBLE PRECISION X(LDX,1),S(1)
      INTEGER I,J,M,MP1
      DOUBLE PRECISION FN,FP
      DOUBLE PRECISION FAC,RU,T
      FP = DBLE(FLOAT(P))
      FN = DBLE(FLOAT(N))
      M = MIN0(N,P)
      RU = 1.0D0
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
