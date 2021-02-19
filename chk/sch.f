
C     MAIN PROGRAM
      INTEGER LUNIT
C     ALLOW 5000 UNDERFLOWS.
      CALL TRAPS(0,0,5001,0,0)
C
C     OUTPUT UNIT NUMBER
C
      LUNIT = 6
C
      CALL SCHTS(LUNIT)
      STOP
      END
      SUBROUTINE SCHTS(LUNIT)
C     LUNIT IS THE OUTPUT UNIT NUMBER
C
C     TESTS
C        SCHDC
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     J.J. DONGARRA AND G.W. STEWART, ARGONNE NATIONAL LABORATORY,
C     UNIVERSITY OF MARYLAND
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK SCHDC
C     EXTERNAL SAGEN,SCHXX,SMACH,SARRAY
C     FORTRAN AMAX1,ABS,FLOAT
C
C     INTERNAL VARIABLES
C
      INTEGER LUNIT
      INTEGER CASE,LDA,P,JPVT(25),KPVT(25),KD
      REAL RESDUL
      REAL A(25,25),AA(25,25),ASAVE(25,25),D(25),WORK(25)
      REAL SMACH,TINY,HUGE,EPS
      LOGICAL NOTWRT
      LDA = 25
      NCASE = 7
      ERRLVL = 100.0E0
      NOTWRT = .TRUE.
      EPS = SMACH(1)
      TINY = SMACH(2)*10000.0E0
      HUGE = SMACH(3)
      WRITE (LUNIT,300)
      DO 290 CASE = 1, NCASE
         WRITE (LUNIT,310) CASE
         GO TO (10,50,80,120,160,200,240), CASE
   10    CONTINUE
C
C           5 X 5 CASE NO PIVOTING.
C
            WRITE (LUNIT,320)
            P = 5
C
C           GENERATE MATRIX WITH POSITIVE EIGENVALUES.
C
            DO 20 I = 1, P
               D(I) = FLOAT(I)
               JPVT(I) = 0
   20       CONTINUE
            CALL SAGEN(A,LDA,P,D,ASAVE)
            IF (NOTWRT) GO TO 30
               WRITE (LUNIT,390)
               CALL SARRAY(A,LDA,P,P,P,LUNIT)
   30       CONTINUE
C
C           DECOMPOSE THE MATRIX.
C
            JOB = 0
            WRITE (LUNIT,420) JOB,(JPVT(I), I = 1, P)
            CALL SCHDC(A,LDA,P,WORK,JPVT,JOB,KD)
            WRITE (LUNIT,430) KD
            CALL SCHXX(A,LDA,P,AA,JPVT,KPVT,ASAVE,RESDUL,JOB)
            WRITE (LUNIT,400) (JPVT(I), I = 1, P)
            RESDUL = RESDUL/EPS
            WRITE (LUNIT,410) RESDUL
            IF (NOTWRT) GO TO 40
               WRITE (LUNIT,390)
               CALL SARRAY(A,LDA,P,P,P,LUNIT)
   40       CONTINUE
         GO TO 280
   50    CONTINUE
C
C           1 X 1 CASE WITH PIVOTING.
C
            WRITE (LUNIT,330)
            P = 1
C
C           GENERATE MATRIX WITH POSITIVE EIGENVALUES.
C
            A(1,1) = 1.0E0
            ASAVE(1,1) = A(1,1)
            JPVT(1) = 0
            IF (NOTWRT) GO TO 60
               WRITE (LUNIT,390)
               CALL SARRAY(A,LDA,P,P,P,LUNIT)
   60       CONTINUE
C
C           DECOMPOSE THE MATRIX.
C
            JOB = 1
            WRITE (LUNIT,420) JOB,(JPVT(I), I = 1, P)
            CALL SCHDC(A,LDA,P,WORK,JPVT,JOB,KD)
            WRITE (LUNIT,430) KD
            CALL SCHXX(A,LDA,P,AA,JPVT,KPVT,ASAVE,RESDUL,JOB)
            WRITE (LUNIT,400) (JPVT(I), I = 1, P)
            RESDUL = RESDUL/EPS
            WRITE (LUNIT,410) RESDUL
            IF (NOTWRT) GO TO 70
               WRITE (LUNIT,390)
               CALL SARRAY(A,LDA,P,P,P,LUNIT)
   70       CONTINUE
         GO TO 280
   80    CONTINUE
C
C           8 X 8 CASE PIVOT LOGIC TEST
C
            WRITE (LUNIT,340)
            P = 8
C
C           GENERATE MATRIX WITH POSITIVE EIGENVALUES.
C
            DO 90 I = 1, P
               D(I) = FLOAT(I)
               JPVT(I) = 1 - 2*MOD(I+1,2)
   90       CONTINUE
            CALL SAGEN(A,LDA,P,D,ASAVE)
            IF (NOTWRT) GO TO 100
               WRITE (LUNIT,390)
               CALL SARRAY(A,LDA,P,P,P,LUNIT)
  100       CONTINUE
C
C           DECOMPOSE THE MATRIX.
C
            JOB = 1
            WRITE (LUNIT,420) JOB,(JPVT(I), I = 1, P)
            CALL SCHDC(A,LDA,P,WORK,JPVT,JOB,KD)
            WRITE (LUNIT,430) KD
            CALL SCHXX(A,LDA,P,AA,JPVT,KPVT,ASAVE,RESDUL,JOB)
            WRITE (LUNIT,400) (JPVT(I), I = 1, P)
            RESDUL = RESDUL/EPS
            WRITE (LUNIT,410) RESDUL
            IF (NOTWRT) GO TO 110
               WRITE (LUNIT,390)
               CALL SARRAY(A,LDA,P,P,P,LUNIT)
  110       CONTINUE
         GO TO 280
  120    CONTINUE
C
C           6 X 6 CASE PIVOTING WITH NEGATIVE EIGENVALUE.
C
            WRITE (LUNIT,350)
            P = 6
C
C           GENERATE MATRIX
C
            D(1) = 1.0E0
            D(2) = 1.0E0/2.0E0
            D(3) = 1.0E0/4.0E0
            D(4) = -1.0E0
            D(5) = 1.0E0/2.0E0
            D(6) = 1.0E0/4.0E0
            DO 130 I = 1, P
               JPVT(I) = 0
  130       CONTINUE
            CALL SAGEN(A,LDA,P,D,ASAVE)
            IF (NOTWRT) GO TO 140
               WRITE (LUNIT,390)
               CALL SARRAY(A,LDA,P,P,P,LUNIT)
  140       CONTINUE
C
C           DECOMPOSE THE MATRIX.
C
            JOB = 1
            WRITE (LUNIT,420) JOB,(JPVT(I), I = 1, P)
            CALL SCHDC(A,LDA,P,WORK,JPVT,JOB,KD)
            WRITE (LUNIT,430) KD
            CALL SCHXX(A,LDA,P,AA,JPVT,KPVT,ASAVE,RESDUL,JOB)
            WRITE (LUNIT,400) (JPVT(I), I = 1, P)
            RESDUL = RESDUL/EPS
            WRITE (LUNIT,410) RESDUL
            IF (NOTWRT) GO TO 150
               WRITE (LUNIT,390)
               CALL SARRAY(A,LDA,P,P,P,LUNIT)
  150       CONTINUE
         GO TO 280
  160    CONTINUE
C
C           25 X 25 CASE WITH PIVOTING.
C
            WRITE (LUNIT,360)
            P = 25
C
C           GENERATE MATRIX WITH POSITIVE EIGENVALUES.
C
            DO 170 I = 1, P
               D(I) = FLOAT(I)
               JPVT(I) = 0
  170       CONTINUE
            CALL SAGEN(A,LDA,P,D,ASAVE)
            IF (NOTWRT) GO TO 180
               WRITE (LUNIT,390)
               CALL SARRAY(A,LDA,P,P,P,LUNIT)
  180       CONTINUE
C
C           DECOMPOSE THE MATRIX.
C
            JOB = 1
            WRITE (LUNIT,420) JOB,(JPVT(I), I = 1, P)
            CALL SCHDC(A,LDA,P,WORK,JPVT,JOB,KD)
            WRITE (LUNIT,430) KD
            CALL SCHXX(A,LDA,P,AA,JPVT,KPVT,ASAVE,RESDUL,JOB)
            WRITE (LUNIT,400) (JPVT(I), I = 1, P)
            RESDUL = RESDUL/EPS
            WRITE (LUNIT,410) RESDUL
            IF (NOTWRT) GO TO 190
               WRITE (LUNIT,390)
               CALL SARRAY(A,LDA,P,P,P,LUNIT)
  190       CONTINUE
         GO TO 280
  200    CONTINUE
C
C           5 X 5 CASE WITH PIVOTING AND UNDERFLOW CHECK.
C
            WRITE (LUNIT,370)
            P = 5
C
C           GENERATE MATRIX WITH POSITIVE EIGENVALUES.
C
            DO 210 I = 1, P
               D(I) = FLOAT(I)*TINY
               JPVT(I) = 0
  210       CONTINUE
            CALL SAGEN(A,LDA,P,D,ASAVE)
            IF (NOTWRT) GO TO 220
               WRITE (LUNIT,390)
               CALL SARRAY(A,LDA,P,P,P,LUNIT)
  220       CONTINUE
C
C           DECOMPOSE THE MATRIX.
C
            JOB = 1
            WRITE (LUNIT,420) JOB,(JPVT(I), I = 1, P)
            CALL SCHDC(A,LDA,P,WORK,JPVT,JOB,KD)
            WRITE (LUNIT,430) KD
            CALL SCHXX(A,LDA,P,AA,JPVT,KPVT,ASAVE,RESDUL,JOB)
            WRITE (LUNIT,400) (JPVT(I), I = 1, P)
            RESDUL = RESDUL/EPS
            WRITE (LUNIT,410) RESDUL
            IF (NOTWRT) GO TO 230
               WRITE (LUNIT,390)
               CALL SARRAY(A,LDA,P,P,P,LUNIT)
  230       CONTINUE
         GO TO 280
  240    CONTINUE
C
C           5 X 5 CASE WITH PIVOTING AND OVERFLOW CHECK.
C
            WRITE (LUNIT,380)
            P = 5
C
C           GENERATE MATRIX WITH POSITIVE EIGENVALUES.
C
            DO 250 I = 1, P
               D(I) = FLOAT(I)*HUGE
               JPVT(I) = 0
  250       CONTINUE
            CALL SAGEN(A,LDA,P,D,ASAVE)
            IF (NOTWRT) GO TO 260
               WRITE (LUNIT,390)
               CALL SARRAY(A,LDA,P,P,P,LUNIT)
  260       CONTINUE
C
C           DECOMPOSE THE MATRIX.
C
            JOB = 1
            WRITE (LUNIT,420) JOB,(JPVT(I), I = 1, P)
            CALL SCHDC(A,LDA,P,WORK,JPVT,JOB,KD)
            WRITE (LUNIT,430) KD
            CALL SCHXX(A,LDA,P,AA,JPVT,KPVT,ASAVE,RESDUL,JOB)
            WRITE (LUNIT,400) (JPVT(I), I = 1, P)
            RESDUL = RESDUL/EPS
            WRITE (LUNIT,410) RESDUL
            IF (NOTWRT) GO TO 270
               WRITE (LUNIT,390)
               CALL SARRAY(A,LDA,P,P,P,LUNIT)
  270       CONTINUE
  280    CONTINUE
  290 CONTINUE
      WRITE (LUNIT,440)
      RETURN
  300 FORMAT (22H1LINPACK TESTER, SCH** /
     *        29H THIS VERSION DATED 08/14/78.)
  310 FORMAT ( // 5H1CASE, I3)
  320 FORMAT ( / 19H 5 X 5 NO PIVOTING.)
  330 FORMAT ( / 22H MONOELEMENTAL MATRIX.)
  340 FORMAT ( / 24H 8 X 8 PIVOT LOGIC TEST.)
  350 FORMAT ( / 32H 6 X 6 NEGATIVE EIGENVALUE TEST.)
  360 FORMAT ( / 16H 25 X 25 MATRIX.)
  370 FORMAT ( / 32H 5 X 5 PIVOT AND UNDERFLOW TEST.)
  380 FORMAT ( / 31H 5 X 5 PIVOT AND OVERFLOW TEST.)
  390 FORMAT ( / 2H A)
  400 FORMAT ( / 5H JPVT // (1H , 10I5))
  410 FORMAT ( // 18H A - TRANS(R)*R = , 1PE16.8)
  420 FORMAT ( / 39H JOB AND JPVT BEFORE THE DECOMPOSITION. / I5 /
     *         (10I5))
  430 FORMAT ( / 20H THE VALUE OF  KD  =, I5)
  440 FORMAT ( /// 12H END OF TEST / 1H1)
      END
      SUBROUTINE SCHXX(A,LDA,P,AA,JPVT,KPVT,ASAVE,RESDUL,JOB)
C
      INTEGER LDA,P,JPVT(1),KPVT(1),JOB
      REAL A(LDA,1),AA(LDA,1),ASAVE(LDA,1)
      REAL RESDUL,ANORM,SASUM
C
      REAL TEMP,SDOT
C
C     FORM TRANS(R)*R
C
      DO 20 I = 1, P
         DO 10 J = I, P
            AA(I,J) = SDOT(I,A(1,I),1,A(1,J),1)
            AA(J,I) = AA(I,J)
   10    CONTINUE
         KPVT(I) = JPVT(I)
   20 CONTINUE
      IF (JOB .EQ. 0) GO TO 70
C
C        UNSCRAMBLE TRANS(R)*R
C
         DO 60 J = 1, P
   30       IF (KPVT(J) .EQ. J) GO TO 50
               IK = KPVT(J)
               CALL SSWAP(P,AA(1,J),1,AA(1,IK),1)
               DO 40 I = 1, P
                  TEMP = AA(J,I)
                  AA(J,I) = AA(IK,I)
                  AA(IK,I) = TEMP
   40          CONTINUE
               KPVT(J) = KPVT(IK)
               KPVT(IK) = IK
            GO TO 30
   50       CONTINUE
   60    CONTINUE
   70 CONTINUE
      ANORM = 0.0E0
      RESDUL = 0.0E0
      DO 90 J = 1, P
         ANORM = AMAX1(ANORM,SASUM(P,ASAVE(1,J),1))
         DO 80 I = 1, P
            ASAVE(I,J) = ASAVE(I,J) - AA(I,J)
   80    CONTINUE
         RESDUL = AMAX1(RESDUL,SASUM(P,ASAVE(1,J),1))
   90 CONTINUE
      RESDUL = RESDUL/ANORM
      RETURN
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
      SUBROUTINE SAGEN(A,LDA,P,D,ASAVE)
      INTEGER LDA,P
      REAL A(LDA,1),D(1),ASAVE(LDA,1)
C
      REAL EN,ES,FDPTP,S,TDP
C
      S = 0.0E0
      DO 10 I = 1, P
         S = S + D(I)
   10 CONTINUE
      TDP = 2.0E0/FLOAT(P)
      FDPTP = 4.0E0*S/FLOAT(P*P)
      DO 30 I = 1, P
         ES = -1.0E0
         DO 20 J = I, P
            ES = (-1.0E0)*ES
            A(I,J) = -ES*TDP*(D(J) + D(I)) + ES*FDPTP
            A(I,J) = A(I,J)*FLOAT(I+1)*FLOAT(J+1)/FLOAT((J+1)**2+J**2)
            ASAVE(I,J) = A(I,J)
            A(J,I) = A(I,J)
            ASAVE(J,I) = A(J,I)
   20    CONTINUE
         A(I,I) = A(I,I) + D(I)
         ASAVE(I,I) = A(I,I)
   30 CONTINUE
      RETURN
      END

