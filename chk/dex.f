C     MAIN PROGRAM
      INTEGER LUNIT
C     ALLOW 5000 UNDERFLOWS.
      CALL TRAPS(0,0,5001,0,0)
C
C     OUTPUT UNIT NUMBER
C
      LUNIT = 6
C
      CALL DQREE(LUNIT)
      STOP
      END
      SUBROUTINE DQREE(LUNIT)
      DOUBLE PRECISION R(10,10),RX(11,10),RHR(10,10),Z(10,4),S(10)
      DOUBLE PRECISION C(10)
      LOGICAL NOTWRT
      DOUBLE PRECISION DMACH,SCALE
      COMMON SCALE,NOTWRT
      INTEGER P
      NOTWRT = .TRUE.
      SCALE = 1.0D0
      LDR = 10
      LDX = 11
      LDZ = 10
      P = 10
      WRITE (LUNIT,280) P
      CALL DGETEX(RX,LDX,Z,LDZ,P,S,LUNIT)
      CALL DRHXR(RX,LDX,P,RHR,LDR)
      CALL DRHXZ(RX,LDX,P,Z(1,1),Z(1,3))
      DO 80 I = 1, 3
         GO TO (10,20,30), I
   10    CONTINUE
            K = 1
            L = 10
         GO TO 40
   20    CONTINUE
            K = 5
            L = 6
         GO TO 40
   30    CONTINUE
            K = 3
            L = 7
   40    CONTINUE
         DO 70 JOB = 1, 2
            DO 60 J1 = 1, P
               Z(J1,2) = Z(J1,1)
               DO 50 I1 = 1, P
                  R(I1,J1) = RX(I1,J1)
   50          CONTINUE
   60       CONTINUE
            CALL DCHEX(R,LDR,P,K,L,Z(1,2),LDZ,1,C,S,JOB)
            CALL DCMPEX(R,LDR,RHR,LDR,Z(1,2),LDZ,P,JOB,K,L,LUNIT)
   70    CONTINUE
   80 CONTINUE
      P = 2
      WRITE (LUNIT,280) P
      CALL DGETEX(RX,LDX,Z,LDZ,P,S,LUNIT)
      CALL DRHXR(RX,LDX,P,RHR,LDR)
      CALL DRHXZ(RX,LDX,P,Z(1,1),Z(1,3))
      K = 1
      L = 2
      DO 110 JOB = 1, 2
         DO 100 J1 = 1, P
            Z(J1,2) = Z(J1,1)
            DO 90 I1 = 1, P
               R(I1,J1) = RX(I1,J1)
   90       CONTINUE
  100    CONTINUE
         CALL DCHEX(R,LDR,P,K,L,Z(1,2),LDZ,1,C,S,JOB)
         CALL DCMPEX(R,LDR,RHR,LDR,Z(1,2),LDZ,P,JOB,K,L,LUNIT)
  110 CONTINUE
      SCALE = DMACH(3)
      P = 10
      WRITE (LUNIT,280) P
      WRITE (LUNIT,290)
      CALL DGETEX(RX,LDX,Z,LDZ,P,S,LUNIT)
      CALL DRHXR(RX,LDX,P,RHR,LDR)
      CALL DRHXZ(RX,LDX,P,Z(1,1),Z(1,3))
      DO 190 I = 1, 3
         GO TO (120,130,140), I
  120    CONTINUE
            K = 1
            L = 10
         GO TO 150
  130    CONTINUE
            K = 5
            L = 6
         GO TO 150
  140    CONTINUE
            K = 3
            L = 7
  150    CONTINUE
         DO 180 JOB = 1, 2
            DO 170 J1 = 1, P
               Z(J1,2) = Z(J1,1)
               DO 160 I1 = 1, P
                  R(I1,J1) = RX(I1,J1)
  160          CONTINUE
  170       CONTINUE
            CALL DCHEX(R,LDR,P,K,L,Z(1,2),LDZ,1,C,S,JOB)
            CALL DCMPEX(R,LDR,RHR,LDR,Z(1,2),LDZ,P,JOB,K,L,LUNIT)
  180    CONTINUE
  190 CONTINUE
      SCALE = DMACH(2)
      P = 10
      WRITE (LUNIT,280) P
      WRITE (LUNIT,300)
      CALL DGETEX(RX,LDX,Z,LDZ,P,S,LUNIT)
      CALL DRHXR(RX,LDX,P,RHR,LDR)
      CALL DRHXZ(RX,LDX,P,Z(1,1),Z(1,3))
      DO 270 I = 1, 3
         GO TO (200,210,220), I
  200    CONTINUE
            K = 1
            L = 10
         GO TO 230
  210    CONTINUE
            K = 5
            L = 6
         GO TO 230
  220    CONTINUE
            K = 3
            L = 7
  230    CONTINUE
         DO 260 JOB = 1, 2
            DO 250 J1 = 1, P
               Z(J1,2) = Z(J1,1)
               DO 240 I1 = 1, P
                  R(I1,J1) = RX(I1,J1)
  240          CONTINUE
  250       CONTINUE
            CALL DCHEX(R,LDR,P,K,L,Z(1,2),LDZ,1,C,S,JOB)
            CALL DCMPEX(R,LDR,RHR,LDR,Z(1,2),LDZ,P,JOB,K,L,LUNIT)
  260    CONTINUE
  270 CONTINUE
      RETURN
C
C
  280 FORMAT (1H1, 3X, 22H *** DCHEX ***     P =, I3 ///)
  290 FORMAT (19H      OVERFLOW TEST ///)
  300 FORMAT (20H      UNDERFLOW TEST ///)
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
      SUBROUTINE DCMPEX(R,LDR,RHR,LDRHR,Z,LDZ,P,JOB,K,L,LUNIT)
      INTEGER LDR,LDRHR,LDZ,P,JOB,K,L,LUNIT
      DOUBLE PRECISION R(LDR,1),RHR(LDRHR,1),Z(LDZ,1)
      INTEGER I,J,JJ,LM1,PM1
      DOUBLE PRECISION ERRCR,ERRCZ,ERMAX,EZMAX,RMAX,ZMAX
      DOUBLE PRECISION T
      LOGICAL NOTWRT
      DOUBLE PRECISION SCALE,DMACH
      COMMON SCALE,NOTWRT
      IF (NOTWRT) GO TO 10
         WRITE (LUNIT,150)
         CALL DARRAY(R,LDR,P,P,P,LUNIT)
         WRITE (LUNIT,160)
         CALL DARRAY(Z,LDZ,P,1,-P,LUNIT)
   10 CONTINUE
      PM1 = P - 1
      LM1 = L - 1
      DO 30 J = 1, PM1
         JJ = J + 1
         DO 20 I = JJ, P
            R(I,J) = 0.0D0
   20    CONTINUE
   30 CONTINUE
      CALL DRHXZ(R,LDR,P,Z(1,1),Z(1,3))
      CALL DRHXR(R,LDR,P,R,LDR)
      DO 50 J = 1, P
         DO 40 I = 1, P
            R(I,J) = R(I,J)/SCALE
   40    CONTINUE
   50 CONTINUE
      IF (JOB .EQ. 2) GO TO 80
         DO 60 J = K, LM1
            CALL DSWAP(P,R(1,J),1,R(1,J+1),1)
            CALL DSWAP(1,Z(J,3),1,Z(J+1,3),1)
   60    CONTINUE
         DO 70 J = K, LM1
            CALL DSWAP(P,R(J,1),LDR,R(J+1,1),LDR)
   70    CONTINUE
      GO TO 110
   80 CONTINUE
         DO 90 I = K, LM1
            J = LM1 + K - I
            CALL DSWAP(P,R(1,J),1,R(1,J+1),1)
            CALL DSWAP(1,Z(J,3),1,Z(J+1,3),1)
   90    CONTINUE
         DO 100 I = K, LM1
            J = LM1 + K - I
            CALL DSWAP(P,R(J,1),LDR,R(J+1,1),LDR)
  100    CONTINUE
  110 CONTINUE
      ERMAX = 0.0D0
      EZMAX = 0.0D0
      RMAX = 0.0D0
      ZMAX = 0.0D0
      DO 130 J = 1, P
         T = Z(J,2)
         ZMAX = DMAX1(ZMAX,DMAX1(DABS(T),0.0D0))
         T = Z(J,2) - Z(J,3)
         EZMAX = DMAX1(EZMAX,DMAX1(DABS(T),0.0D0))
         DO 120 I = 1, P
            T = RHR(I,J)
            RMAX = DMAX1(RMAX,DMAX1(DABS(T),0.0D0))
            T = RHR(I,J) - R(I,J)
            ERMAX = DMAX1(ERMAX,DMAX1(DABS(T),0.0D0))
  120    CONTINUE
  130 CONTINUE
      ERRCR = ERMAX/RMAX/DMACH(1)
      ERRCZ = EZMAX/ZMAX/DMACH(1)
      WRITE (LUNIT,140) JOB,K,L,ERRCR,ERRCZ
      RETURN
  140 FORMAT ( /// 25H     STATISTICS     JOB =, I2, 5H   K=, I2,
     *         5H   L=, I2 / 38H0         RH*R    ................... ,
     *         D10.2 / 38H          RH*Z    ................... ,
     *         D10.2)
  150 FORMAT (6H   REX)
  160 FORMAT (6H   ZEX)
      END
      SUBROUTINE DGETEX(X,LDX,Z,LDZ,P,S,LUNIT)
      INTEGER LDX,LDZ,P,LUNIT,JOB
      LOGICAL NOTWRT
      DOUBLE PRECISION SCALE
      COMMON SCALE,NOTWRT
      DOUBLE PRECISION X(LDX,1),Z(LDZ,1),S(1)
      INTEGER I,JJ,J,PP1,PM1
      PP1 = P + 1
      PM1 = P - 1
      CALL DGETRX(X,LDX,PP1,P,S)
      DO 10 I = 1, P
         Z(I,1) = X(PP1,I)*SCALE
   10 CONTINUE
      DO 30 J = 1, PM1
         JJ = J + 1
         DO 20 I = JJ, P
            X(I,J) = 0.0D0
   20    CONTINUE
   30 CONTINUE
      DO 50 J = 1, P
         DO 40 I = 1, P
            X(I,J) = X(I,J)*SCALE
   40    CONTINUE
   50 CONTINUE
      IF (NOTWRT) GO TO 60
         WRITE (LUNIT,80)
         CALL DARRAY(X,LDX,P,P,P,LUNIT)
         WRITE (LUNIT,70)
         CALL DARRAY(Z,LDZ,P,1,-P,LUNIT)
   60 CONTINUE
      RETURN
C
   70 FORMAT ( /// 4H   Z)
   80 FORMAT (4H   R)
      END
      SUBROUTINE DGETRX(X,LDX,N,P,S)
      INTEGER N,P,LDX
      DOUBLE PRECISION X(LDX,1),S(1)
      INTEGER PD2,PD2D1,I
      PD2 = P/2
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
      SUBROUTINE DRHXR(R,LDR,P,RHR,LDRHR)
      INTEGER LDR,LDRHR,P
      DOUBLE PRECISION R(LDR,1),RHR(LDRHR,1)
      LOGICAL DUMMY
      DOUBLE PRECISION SCALE
      COMMON SCALE,DUMMY
      DOUBLE PRECISION DDOT
      DO 20 J = 1, P
         DO 10 I = 1, J
            R(I,J) = R(I,J)/SCALE
   10    CONTINUE
   20 CONTINUE
      DO 40 J = 1, P
         DO 30 I = J, P
            II = P + J - I
            RHR(II,J) = DDOT(J,R(1,II),1,R(1,J),1)
   30    CONTINUE
   40 CONTINUE
      DO 60 J = 2, P
         JJ = J - 1
         DO 50 I = 1, JJ
            RHR(I,J) = RHR(J,I)
   50    CONTINUE
   60 CONTINUE
      DO 80 J = 1, P
         DO 70 I = 1, P
            R(I,J) = R(I,J)*SCALE
   70    CONTINUE
   80 CONTINUE
      RETURN
      END
      SUBROUTINE DRHXZ(R,LDR,P,Z,RHZ)
      INTEGER LDR,P
      DOUBLE PRECISION R(LDR,1),Z(1),RHZ(1)
      LOGICAL DUMMY
      DOUBLE PRECISION SCALE
      COMMON SCALE,DUMMY
      DOUBLE PRECISION DDOT
      DO 20 J = 1, P
         Z(J) = Z(J)/SCALE
         DO 10 I = 1, J
            R(I,J) = R(I,J)/SCALE
   10    CONTINUE
   20 CONTINUE
      DO 30 J = 1, P
         RHZ(J) = DDOT(J,R(1,J),1,Z(1),1)
   30 CONTINUE
      DO 50 J = 1, P
         Z(J) = Z(J)*SCALE
         DO 40 I = 1, J
            R(I,J) = R(I,J)*SCALE
   40    CONTINUE
   50 CONTINUE
      RETURN
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
