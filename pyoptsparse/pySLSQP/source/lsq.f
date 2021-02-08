      SUBROUTINE LSQ(M,MEQ,N,NL,LA,L,G,A,B,XL,XU,X,Y,W,JW,MODE)

C   MINIMIZE with respect to X

C             ||E*X - F||
C                                      1/2  T
C   WITH UPPER TRIANGULAR MATRIX E = +D   *L ,

C                                      -1/2  -1
C                     AND VECTOR F = -D    *L  *G,

C  WHERE THE UNIT LOWER TRIDIANGULAR MATRIX L IS STORED COLUMNWISE
C  DENSE IN THE N*(N+1)/2 ARRAY L WITH VECTOR D STORED IN ITS
C 'DIAGONAL' THUS SUBSTITUTING THE ONE-ELEMENTS OF L

C   SUBJECT TO

C             A(J)*X - B(J) = 0 ,         J=1,...,MEQ,
C             A(J)*X - B(J) >=0,          J=MEQ+1,...,M,
C             XL(I) <= X(I) <= XU(I),     I=1,...,N,
C     ON ENTRY, THE USER HAS TO PROVIDE THE ARRAYS L, G, A, B, XL, XU.
C     WITH DIMENSIONS: L(N*(N+1)/2), G(N), A(LA,N), B(M), XL(N), XU(N)
C     THE WORKING ARRAY W MUST HAVE AT LEAST THE FOLLOWING DIMENSION:
c     DIM(W) =        (3*N+M)*(N+1)                        for LSQ
c                    +(N-MEQ+1)*(MINEQ+2) + 2*MINEQ        for LSI
c                    +(N+MINEQ)*(N-MEQ) + 2*MEQ + N        for LSEI
c                      with MINEQ = M - MEQ + 2*N
C     ON RETURN, NO ARRAY WILL BE CHANGED BY THE SUBROUTINE.
C     X     STORES THE N-DIMENSIONAL SOLUTION VECTOR
C     Y     STORES THE VECTOR OF LAGRANGE MULTIPLIERS OF DIMENSION
C           M+N+N (CONSTRAINTS+LOWER+UPPER BOUNDS)
C     MODE  IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANINGS:
C          MODE=1: SUCCESSFUL COMPUTATION
C               2: ERROR RETURN BECAUSE OF WRONG DIMENSIONS (N<1)
C               3: ITERATION COUNT EXCEEDED BY NNLS
C               4: INEQUALITY CONSTRAINTS INCOMPATIBLE
C               5: MATRIX E IS NOT OF FULL RANK
C               6: MATRIX C IS NOT OF FULL RANK
C               7: RANK DEFECT IN HFTI

c     coded            Dieter Kraft, april 1987
c     revised                        march 1989

      DOUBLE PRECISION L,G,A,B,W,XL,XU,X,Y,
     .                 DIAG,ZERO,ONE,DDOT,XNORM

      INTEGER          JW(*),I,IC,ID,IE,IF,IG,IH,IL,IM,IP,IU,IW,
     .                 I1,I2,I3,I4,LA,M,MEQ,MINEQ,MODE,M1,N,NL,N1,N2,N3

      DIMENSION        A(LA,N), B(LA), G(N), L(NL),
     .                 W(*), X(N), XL(N), XU(N), Y(M+N+N)

      DATA             ZERO/0.0D0/, ONE/1.0D0/

      N1 = N + 1
      MINEQ = M - MEQ
      M1 = MINEQ + N + N

c  determine whether to solve problem
c  with inconsistent linerarization (n2=1)
c  or not (n2=0)

      N2 = N1*N/2 + 1
      IF (N2.EQ.NL) THEN
          N2 = 0
      ELSE
          N2 = 1
      ENDIF
      N3 = N-N2

C  RECOVER MATRIX E AND VECTOR F FROM L AND G

      I2 = 1
      I3 = 1
      I4 = 1
      IE = 1
      IF = N*N+1
      DO 10 I=1,N3
         I1 = N1-I
         DIAG = SQRT (L(I2))
         W(I3) = ZERO
         CALL DCOPY (I1  ,  W(I3), 0, W(I3), 1)
         CALL DCOPY (I1-N2, L(I2), 1, W(I3), N)
         CALL DSCAL (I1-N2,     DIAG, W(I3), N)
         W(I3) = DIAG
         W(IF-1+I) = (G(I) - DDOT (I-1, W(I4), 1, W(IF), 1))/DIAG
         I2 = I2 + I1 - N2
         I3 = I3 + N1
         I4 = I4 + N
   10 CONTINUE
      IF (N2.EQ.1) THEN
          W(I3) = L(NL)
          W(I4) = ZERO
          CALL DCOPY (N3, W(I4), 0, W(I4), 1)
          W(IF-1+N) = ZERO
      ENDIF
      CALL DSCAL (N, - ONE, W(IF), 1)

      IC = IF + N
      ID = IC + MEQ*N

      IF (MEQ .GT. 0) THEN

C  RECOVER MATRIX C FROM UPPER PART OF A

          DO 20 I=1,MEQ
              CALL DCOPY (N, A(I,1), LA, W(IC-1+I), MEQ)
   20     CONTINUE

C  RECOVER VECTOR D FROM UPPER PART OF B

          CALL DCOPY (MEQ, B(1), 1, W(ID), 1)
          CALL DSCAL (MEQ,   - ONE, W(ID), 1)

      ENDIF

      IG = ID + MEQ

      IF (MINEQ .GT. 0) THEN

C  RECOVER MATRIX G FROM LOWER PART OF A

          DO 30 I=1,MINEQ
              CALL DCOPY (N, A(MEQ+I,1), LA, W(IG-1+I), M1)
   30     CONTINUE

      ENDIF

C  AUGMENT MATRIX G BY +I AND -I

      IP = IG + MINEQ
      DO 40 I=1,N
         W(IP-1+I) = ZERO
         CALL DCOPY (N, W(IP-1+I), 0, W(IP-1+I), M1)
   40 CONTINUE
      W(IP) = ONE
      CALL DCOPY (N, W(IP), 0, W(IP), M1+1)

      IM = IP + N
      DO 50 I=1,N
         W(IM-1+I) = ZERO
         CALL DCOPY (N, W(IM-1+I), 0, W(IM-1+I), M1)
   50 CONTINUE
      W(IM) = -ONE
      CALL DCOPY (N, W(IM), 0, W(IM), M1+1)

      IH = IG + M1*N

      IF (MINEQ .GT. 0) THEN

C  RECOVER H FROM LOWER PART OF B

          CALL DCOPY (MINEQ, B(MEQ+1), 1, W(IH), 1)
          CALL DSCAL (MINEQ,       - ONE, W(IH), 1)

      ENDIF

C  AUGMENT VECTOR H BY XL AND XU

      IL = IH + MINEQ
      CALL DCOPY (N, XL, 1, W(IL), 1)
      IU = IL + N
      CALL DCOPY (N, XU, 1, W(IU), 1)
      CALL DSCAL (N, - ONE, W(IU), 1)

      IW = IU + N

      CALL LSEI (W(IC), W(ID), W(IE), W(IF), W(IG), W(IH), MAX(1,MEQ),
     .           MEQ, N, N, M1, M1, N, X, XNORM, W(IW), JW, MODE)

      IF (MODE .EQ. 1) THEN

c   restore Lagrange multipliers

          CALL DCOPY (M,  W(IW),     1, Y(1),      1)
          CALL DCOPY (N3, W(IW+M),   1, Y(M+1),    1)
          CALL DCOPY (N3, W(IW+M+N), 1, Y(M+N3+1), 1)

      ENDIF

C   END OF SUBROUTINE LSQ

      END
      