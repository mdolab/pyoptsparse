      SUBROUTINE NNLS (A, MDA, M, N, B, X, RNORM, W, Z, INDEX, MODE)

C     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY:
C     'SOLVING LEAST SQUARES PROBLEMS'. PRENTICE-HALL.1974

C      **********   NONNEGATIVE LEAST SQUARES   **********

C     GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B, COMPUTE AN
C     N-VECTOR, X, WHICH SOLVES THE LEAST SQUARES PROBLEM

C                  A*X = B  SUBJECT TO  X >= 0

C     A(),MDA,M,N
C            MDA IS THE FIRST DIMENSIONING PARAMETER FOR THE ARRAY,A().
C            ON ENTRY A()  CONTAINS THE M BY N MATRIX,A.
C            ON EXIT A() CONTAINS THE PRODUCT Q*A,
C            WHERE Q IS AN M BY M ORTHOGONAL MATRIX GENERATED
C            IMPLICITLY BY THIS SUBROUTINE.
C            EITHER M>=N OR M<N IS PERMISSIBLE.
C            THERE IS NO RESTRICTION ON THE RANK OF A.
C     B()    ON ENTRY B() CONTAINS THE M-VECTOR, B.
C            ON EXIT B() CONTAINS Q*B.
C     X()    ON ENTRY X() NEED NOT BE INITIALIZED.
C            ON EXIT X() WILL CONTAIN THE SOLUTION VECTOR.
C     RNORM  ON EXIT RNORM CONTAINS THE EUCLIDEAN NORM OF THE
C            RESIDUAL VECTOR.
C     W()    AN N-ARRAY OF WORKING SPACE.
C            ON EXIT W() WILL CONTAIN THE DUAL SOLUTION VECTOR.
C            W WILL SATISFY W(I)=0 FOR ALL I IN SET P
C            AND W(I)<=0 FOR ALL I IN SET Z
C     Z()    AN M-ARRAY OF WORKING SPACE.
C     INDEX()AN INTEGER WORKING ARRAY OF LENGTH AT LEAST N.
C            ON EXIT THE CONTENTS OF THIS ARRAY DEFINE THE SETS
C            P AND Z AS FOLLOWS:
C            INDEX(1)    THRU INDEX(NSETP) = SET P.
C            INDEX(IZ1)  THRU INDEX (IZ2)  = SET Z.
C            IZ1=NSETP + 1 = NPP1, IZ2=N.
C     MODE   THIS IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANING:
C            1    THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY.
C            2    THE DIMENSIONS OF THE PROBLEM ARE WRONG,
C                 EITHER M <= 0 OR N <= 0.
C            3    ITERATION COUNT EXCEEDED, MORE THAN 3*N ITERATIONS.

      INTEGER          I,II,IP,ITER,ITMAX,IZ,IZMAX,IZ1,IZ2,J,JJ,JZ,
     *                 K,L,M,MDA,MODE,N,NPP1,NSETP,INDEX(N)

      DOUBLE PRECISION A(MDA,N),B(M),X(N),W(N),Z(M),ASAVE,DIFF,
     *                 FACTOR,DDOT,ZERO,ONE,WMAX,ALPHA,
     *                 C,S,T,U,V,UP,RNORM,UNORM,DNRM2

      DIFF(U,V)=       U-V

      DATA             ZERO,ONE,FACTOR/0.0D0,1.0D0,1.0D-2/

c     revised          Dieter Kraft, March 1983

      MODE=2
      IF(M.LE.0.OR.N.LE.0)            GOTO 290
      MODE=1
      ITER=0
      ITMAX=3*N

C STEP ONE (INITIALIZE)

      DO 100 I=1,N
  100    INDEX(I)=I
      IZ1=1
      IZ2=N
      NSETP=0
      NPP1=1
      X(1)=ZERO
      CALL DCOPY(N,X(1),0,X,1)

C STEP TWO (COMPUTE DUAL VARIABLES)
C .....ENTRY LOOP A

  110 IF(IZ1.GT.IZ2.OR.NSETP.GE.M)    GOTO 280
      DO 120 IZ=IZ1,IZ2
         J=INDEX(IZ)
  120    W(J)=DDOT(M-NSETP,A(NPP1,J),1,B(NPP1),1)

C STEP THREE (TEST DUAL VARIABLES)

  130 WMAX=ZERO
      DO 140 IZ=IZ1,IZ2
      J=INDEX(IZ)
         IF(W(J).LE.WMAX)             GOTO 140
         WMAX=W(J)
         IZMAX=IZ
  140 CONTINUE

C .....EXIT LOOP A

      IF(WMAX.LE.ZERO)                GOTO 280
      IZ=IZMAX
      J=INDEX(IZ)

C STEP FOUR (TEST INDEX J FOR LINEAR DEPENDENCY)

      ASAVE=A(NPP1,J)
      CALL H12(1,NPP1,NPP1+1,M,A(1,J),1,UP,Z,1,1,0)
      UNORM=DNRM2(NSETP,A(1,J),1)
      T=FACTOR*ABS(A(NPP1,J))
      IF(DIFF(UNORM+T,UNORM).LE.ZERO) GOTO 150
      CALL DCOPY(M,B,1,Z,1)
      CALL H12(2,NPP1,NPP1+1,M,A(1,J),1,UP,Z,1,1,1)
      IF(Z(NPP1)/A(NPP1,J).GT.ZERO)   GOTO 160
  150 A(NPP1,J)=ASAVE
      W(J)=ZERO
                                      GOTO 130
C STEP FIVE (ADD COLUMN)

  160 CALL DCOPY(M,Z,1,B,1)
      INDEX(IZ)=INDEX(IZ1)
      INDEX(IZ1)=J
      IZ1=IZ1+1
      NSETP=NPP1
      NPP1=NPP1+1
      DO 170 JZ=IZ1,IZ2
         JJ=INDEX(JZ)
  170    CALL H12(2,NSETP,NPP1,M,A(1,J),1,UP,A(1,JJ),1,MDA,1)
      K=MIN(NPP1,MDA)
      W(J)=ZERO
      CALL DCOPY(M-NSETP,W(J),0,A(K,J),1)

C STEP SIX (SOLVE LEAST SQUARES SUB-PROBLEM)
C .....ENTRY LOOP B

  180 DO 200 IP=NSETP,1,-1
         IF(IP.EQ.NSETP)              GOTO 190
         CALL DAXPY(IP,-Z(IP+1),A(1,JJ),1,Z,1)
  190    JJ=INDEX(IP)
  200    Z(IP)=Z(IP)/A(IP,JJ)
      ITER=ITER+1
      IF(ITER.LE.ITMAX)               GOTO 220
  210 MODE=3
                                      GOTO 280
C STEP SEVEN TO TEN (STEP LENGTH ALGORITHM)

  220 ALPHA=ONE
      JJ=0
      DO 230 IP=1,NSETP
         IF(Z(IP).GT.ZERO)            GOTO 230
         L=INDEX(IP)
         T=-X(L)/(Z(IP)-X(L))
         IF(ALPHA.LT.T)               GOTO 230
         ALPHA=T
         JJ=IP
  230 CONTINUE
      DO 240 IP=1,NSETP
         L=INDEX(IP)
  240    X(L)=(ONE-ALPHA)*X(L) + ALPHA*Z(IP)

C .....EXIT LOOP B

      IF(JJ.EQ.0)                     GOTO 110

C STEP ELEVEN (DELETE COLUMN)

      I=INDEX(JJ)
  250 X(I)=ZERO
      JJ=JJ+1
      DO 260 J=JJ,NSETP
         II=INDEX(J)
         INDEX(J-1)=II
         CALL DROTG(A(J-1,II),A(J,II),C,S)
         T=A(J-1,II)
         CALL DROT(N,A(J-1,1),MDA,A(J,1),MDA,C,S)
         A(J-1,II)=T
         A(J,II)=ZERO
  260    CALL DROT(1,B(J-1),1,B(J),1,C,S)
      NPP1=NSETP
      NSETP=NSETP-1
      IZ1=IZ1-1
      INDEX(IZ1)=I
      IF(NSETP.LE.0)                  GOTO 210
      DO 270 JJ=1,NSETP
         I=INDEX(JJ)
         IF(X(I).LE.ZERO)             GOTO 250
  270 CONTINUE
      CALL DCOPY(M,B,1,Z,1)
                                      GOTO 180
C STEP TWELVE (SOLUTION)

  280 K=MIN(NPP1,M)
      RNORM=DNRM2(M-NSETP,B(K),1)
      IF(NPP1.GT.M) THEN
         W(1)=ZERO
         CALL DCOPY(N,W(1),0,W,1)
      ENDIF

C END OF SUBROUTINE NNLS

  290 END
      