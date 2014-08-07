      SUBROUTINE SLSQPB (M, MEQ, LA, N, X, XL, XU, F, C, G, A, ACC,
     *                   ITER, MODE, R, L, X0, MU, S, U, V, W, IW,
     *                   ALPHA, F0, GS, H1, H2, H3, H4, T, T0, TOL,
     *                   IEXACT, INCONS, IRESET, ITERMX, LINE, N1, 
     *                   N2, N3)

C   NONLINEAR PROGRAMMING BY SOLVING SEQUENTIALLY QUADRATIC PROGRAMS

C        -  L1 - LINE SEARCH,  POSITIVE DEFINITE  BFGS UPDATE  -

C                      BODY SUBROUTINE FOR SLSQP

      INTEGER          IW(*), I, IEXACT, INCONS, IRESET, ITER, ITERMX,
     *                 K, J, LA, LINE, M, MEQ, MODE, N, N1, N2, N3

      DOUBLE PRECISION A(LA,N+1), C(LA), G(N+1), L((N+1)*(N+2)/2),
     *                 MU(LA), R(M+N+N+2), S(N+1), U(N+1), V(N+1), W(*),
     *                 X(N), XL(N), XU(N), X0(N),
     *                 DDOT, DNRM2, LINMIN,
     *                 ACC, ALFMIN, ALPHA, F, F0, GS, H1, H2, H3, H4,
     *                 HUN, ONE, T, T0, TEN, TOL, TWO, ZERO

c     dim(W) =         N1*(N1+1) + MEQ*(N1+1) + MINEQ*(N1+1)  for LSQ
c                     +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ
c                     +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1       for LSEI
c                      with MINEQ = M - MEQ + 2*N1  &  N1 = N+1

      DATA             ZERO /0.0D0/, ONE /1.0D0/, ALFMIN /1.0D-1/,
     *                 HUN /1.0D+2/, TEN /1.0D+1/, TWO /2.0D0/

      IF (MODE) 260, 100, 220

  100 ITERMX = ITER
      IF (ACC.GE.ZERO) THEN
          IEXACT = 0
      ELSE
          IEXACT = 1
      ENDIF
      ACC = ABS(ACC)
      TOL = TEN*ACC
      ITER = 0
      IRESET = 0
      N1 = N + 1
      N2 = N1*N/2
      N3 = N2 + 1
      S(1) = ZERO
      MU(1) = ZERO
      CALL DCOPY(N, S(1),  0, S,  1)
      CALL DCOPY(M, MU(1), 0, MU, 1)

C   RESET BFGS MATRIX

  110 IRESET = IRESET + 1
      IF (IRESET.GT.5) GO TO 255
      L(1) = ZERO
      CALL DCOPY(N2, L(1), 0, L, 1)
      J = 1
      DO 120 I=1,N
         L(J) = ONE
         J = J + N1 - I
  120 CONTINUE

C   MAIN ITERATION : SEARCH DIRECTION, STEPLENGTH, LDL'-UPDATE

  130 ITER = ITER + 1
      MODE = 9
      IF (ITER.GT.ITERMX) GO TO 330

C   SEARCH DIRECTION AS SOLUTION OF QP - SUBPROBLEM

      CALL DCOPY(N, XL, 1, U, 1)
      CALL DCOPY(N, XU, 1, V, 1)
      CALL DAXPY(N, -ONE, X, 1, U, 1)
      CALL DAXPY(N, -ONE, X, 1, V, 1)
      H4 = ONE
      CALL LSQ (M, MEQ, N , N3, LA, L, G, A, C, U, V, S, R, W, IW, MODE)

C   AUGMENTED PROBLEM FOR INCONSISTENT LINEARIZATION

      IF (MODE.EQ.6) THEN
          IF (N.EQ.MEQ) THEN
              MODE = 4
          ENDIF
      ENDIF
      IF (MODE.EQ.4) THEN
          DO 140 J=1,M
             IF (J.LE.MEQ) THEN
                 A(J,N1) = -C(J)
             ELSE
                 A(J,N1) = MAX(-C(J),ZERO)
             ENDIF
  140     CONTINUE
          S(1) = ZERO
          CALL DCOPY(N, S(1), 0, S, 1)
          H3 = ZERO
          G(N1) = ZERO
          L(N3) = HUN
          S(N1) = ONE
          U(N1) = ZERO
          V(N1) = ONE
          INCONS = 0
  150     CALL LSQ (M, MEQ, N1, N3, LA, L, G, A, C, U, V, S, R,
     *              W, IW, MODE)
          H4 = ONE - S(N1)
          IF (MODE.EQ.4) THEN
              L(N3) = TEN*L(N3)
              INCONS = INCONS + 1
              IF (INCONS.GT.5) GO TO 330
              GOTO 150
          ELSE IF (MODE.NE.1) THEN
              GOTO 330
          ENDIF
      ELSE IF (MODE.NE.1) THEN
          GOTO 330
      ENDIF

C   UPDATE MULTIPLIERS FOR L1-TEST

      DO 160 I=1,N
         V(I) = G(I) - DDOT(M,A(1,I),1,R,1)
  160 CONTINUE
      F0 = F
      CALL DCOPY(N, X, 1, X0, 1)
      GS = DDOT(N, G, 1, S, 1)
      H1 = ABS(GS)
      H2 = ZERO
      DO 170 J=1,M
         IF (J.LE.MEQ) THEN
             H3 = C(J)
         ELSE
             H3 = ZERO
         ENDIF
         H2 = H2 + MAX(-C(J),H3)
         H3 = ABS(R(J))
         MU(J) = MAX(H3,(MU(J)+H3)/TWO)
         H1 = H1 + H3*ABS(C(J))
  170 CONTINUE

C   CHECK CONVERGENCE

      MODE = 0
      IF (H1.LT.ACC .AND. H2.LT.ACC) GO TO 330
      H1 = ZERO
      DO 180 J=1,M
         IF (J.LE.MEQ) THEN
             H3 = C(J)
         ELSE
             H3 = ZERO
         ENDIF
         H1 = H1 + MU(J)*MAX(-C(J),H3)
  180 CONTINUE
      T0 = F + H1
      H3 = GS - H1*H4
      MODE = 8
      IF (H3.GE.ZERO) GO TO 110

C   LINE SEARCH WITH AN L1-TESTFUNCTION

      LINE = 0
      ALPHA = ONE
      IF (IEXACT.EQ.1) GOTO 210

C   INEXACT LINESEARCH

  190     LINE = LINE + 1
          H3 = ALPHA*H3
          CALL DSCAL(N, ALPHA, S, 1)
          CALL DCOPY(N, X0, 1, X, 1)
          CALL DAXPY(N, ONE, S, 1, X, 1)
          MODE = 1
          GO TO 330
  200         IF (H1.LE.H3/TEN .OR. LINE.GT.10) GO TO 240
              ALPHA = MAX(H3/(TWO*(H3-H1)),ALFMIN)
              GO TO 190

C   EXACT LINESEARCH

  210 IF (LINE.NE.3) THEN
          ALPHA = LINMIN(LINE,ALFMIN,ONE,T,TOL)
          CALL DCOPY(N, X0, 1, X, 1)
          CALL DAXPY(N, ALPHA, S, 1, X, 1)
          MODE = 1
          GOTO 330
      ENDIF
      CALL DSCAL(N, ALPHA, S, 1)
      GOTO 240

C   CALL FUNCTIONS AT CURRENT X

  220     T = F
          DO 230 J=1,M
             IF (J.LE.MEQ) THEN
                 H1 = C(J)
             ELSE
                 H1 = ZERO
             ENDIF
             T = T + MU(J)*MAX(-C(J),H1)
  230     CONTINUE
          H1 = T - T0
          GOTO (200, 210) IEXACT+1

C   CHECK CONVERGENCE

  240 H3 = ZERO
      DO 250 J=1,M
         IF (J.LE.MEQ) THEN
             H1 = C(J)
         ELSE
             H1 = ZERO
         ENDIF
         H3 = H3 + MAX(-C(J),H1)
  250 CONTINUE
      IF ((ABS(F-F0).LT.ACC .OR. DNRM2(N,S,1).LT.ACC) .AND. H3.LT.ACC)
     *   THEN
            MODE = 0
         ELSE
            MODE = -1
         ENDIF
      GO TO 330

C   CHECK relaxed CONVERGENCE in case of positive directional derivative

  255 CONTINUE
      IF ((ABS(F-F0).LT.TOL .OR. DNRM2(N,S,1).LT.TOL) .AND. H3.LT.TOL)
     *   THEN
            MODE = 0
         ELSE
            MODE = 8
         ENDIF
      GO TO 330

C   CALL JACOBIAN AT CURRENT X

C   UPDATE CHOLESKY-FACTORS OF HESSIAN MATRIX BY MODIFIED BFGS FORMULA

  260 DO 270 I=1,N
         U(I) = G(I) - DDOT(M,A(1,I),1,R,1) - V(I)
  270 CONTINUE

C   L'*S

      K = 0
      DO 290 I=1,N
         H1 = ZERO
         K = K + 1
         DO 280 J=I+1,N
            K = K + 1
            H1 = H1 + L(K)*S(J)
  280    CONTINUE
         V(I) = S(I) + H1
  290 CONTINUE

C   D*L'*S

      K = 1
      DO 300 I=1,N
         V(I) = L(K)*V(I)
         K = K + N1 - I
  300 CONTINUE

C   L*D*L'*S

      DO 320 I=N,1,-1
         H1 = ZERO
         K = I
         DO 310 J=1,I - 1
            H1 = H1 + L(K)*V(J)
            K = K + N - J
  310    CONTINUE
         V(I) = V(I) + H1
  320 CONTINUE

      H1 = DDOT(N,S,1,U,1)
      H2 = DDOT(N,S,1,V,1)
      H3 = 0.2D0*H2
      IF (H1.LT.H3) THEN
          H4 = (H2-H3)/(H2-H1)
          H1 = H3
          CALL DSCAL(N, H4, U, 1)
          CALL DAXPY(N, ONE-H4, V, 1, U, 1)
      ENDIF
      CALL LDL(N, L, U, +ONE/H1, V)
      CALL LDL(N, L, V, -ONE/H2, U)

C   END OF MAIN ITERATION

      GO TO 130

C   END OF SLSQPB

  330 END
      