      SUBROUTINE LSEI(C,D,E,F,G,H,LC,MC,LE,ME,LG,MG,N,X,XNRM,W,JW,MODE)

C     FOR MODE=1, THE SUBROUTINE RETURNS THE SOLUTION X OF
C     EQUALITY & INEQUALITY CONSTRAINED LEAST SQUARES PROBLEM LSEI :

C                MIN ||E*X - F||
C                 X

C                S.T.  C*X  = D,
C                      G*X >= H.

C     USING QR DECOMPOSITION & ORTHOGONAL BASIS OF NULLSPACE OF C
C     CHAPTER 23.6 OF LAWSON & HANSON: SOLVING LEAST SQUARES PROBLEMS.

C     THE FOLLOWING DIMENSIONS OF THE ARRAYS DEFINING THE PROBLEM
C     ARE NECESSARY
C     DIM(E) :   FORMAL (LE,N),    ACTUAL (ME,N)
C     DIM(F) :   FORMAL (LE  ),    ACTUAL (ME  )
C     DIM(C) :   FORMAL (LC,N),    ACTUAL (MC,N)
C     DIM(D) :   FORMAL (LC  ),    ACTUAL (MC  )
C     DIM(G) :   FORMAL (LG,N),    ACTUAL (MG,N)
C     DIM(H) :   FORMAL (LG  ),    ACTUAL (MG  )
C     DIM(X) :   FORMAL (N   ),    ACTUAL (N   )
C     DIM(W) :   2*MC+ME+(ME+MG)*(N-MC)  for LSEI
C              +(N-MC+1)*(MG+2)+2*MG     for LSI
C     DIM(JW):   MAX(MG,L)
C     ON ENTRY, THE USER HAS TO PROVIDE THE ARRAYS C, D, E, F, G, AND H.
C     ON RETURN, ALL ARRAYS WILL BE CHANGED BY THE SUBROUTINE.
C     X     STORES THE SOLUTION VECTOR
C     XNORM STORES THE RESIDUUM OF THE SOLUTION IN EUCLIDIAN NORM
C     W     STORES THE VECTOR OF LAGRANGE MULTIPLIERS IN ITS FIRST
C           MC+MG ELEMENTS
C     MODE  IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANINGS:
C          MODE=1: SUCCESSFUL COMPUTATION
C               2: ERROR RETURN BECAUSE OF WRONG DIMENSIONS (N<1)
C               3: ITERATION COUNT EXCEEDED BY NNLS
C               4: INEQUALITY CONSTRAINTS INCOMPATIBLE
C               5: MATRIX E IS NOT OF FULL RANK
C               6: MATRIX C IS NOT OF FULL RANK
C               7: RANK DEFECT IN HFTI

C     18.5.1981, DIETER KRAFT, DFVLR OBERPFAFFENHOFEN
C     20.3.1987, DIETER KRAFT, DFVLR OBERPFAFFENHOFEN

      INTEGER          JW(*),I,IE,IF,IG,IW,J,K,KRANK,L,LC,LE,LG,
     .                 MC,MC1,ME,MG,MODE,N
      DOUBLE PRECISION C(LC,N),E(LE,N),G(LG,N),D(LC),F(LE),H(LG),X(N),
     .                 W(*),T,DDOT,XNRM,DNRM2,EPMACH,ZERO
      DATA             EPMACH/2.22D-16/,ZERO/0.0D+00/

      MODE=2
      IF(MC.GT.N)                      GOTO 75
      L=N-MC
      MC1=MC+1
      IW=(L+1)*(MG+2)+2*MG+MC
      IE=IW+MC+1
      IF=IE+ME*L
      IG=IF+ME

C  TRIANGULARIZE C AND APPLY FACTORS TO E AND G

      DO 10 I=1,MC
          J=MIN(I+1,LC)
          CALL H12(1,I,I+1,N,C(I,1),LC,W(IW+I),C(J,1),LC,1,MC-I)
          CALL H12(2,I,I+1,N,C(I,1),LC,W(IW+I),E     ,LE,1,ME)
   10     CALL H12(2,I,I+1,N,C(I,1),LC,W(IW+I),G     ,LG,1,MG)

C  SOLVE C*X=D AND MODIFY F

      MODE=6
      DO 15 I=1,MC
          IF(ABS(C(I,I)).LT.EPMACH)    GOTO 75
          X(I)=(D(I)-DDOT(I-1,C(I,1),LC,X,1))/C(I,I)
   15 CONTINUE
      MODE=1
      W(MC1) = ZERO
      CALL DCOPY (MG-MC,W(MC1),0,W(MC1),1)

      IF(MC.EQ.N)                      GOTO 50

      DO 20 I=1,ME
   20     W(IF-1+I)=F(I)-DDOT(MC,E(I,1),LE,X,1)

C  STORE TRANSFORMED E & G

      DO 25 I=1,ME
   25     CALL DCOPY(L,E(I,MC1),LE,W(IE-1+I),ME)
      DO 30 I=1,MG
   30     CALL DCOPY(L,G(I,MC1),LG,W(IG-1+I),MG)

      IF(MG.GT.0)                      GOTO 40

C  SOLVE LS WITHOUT INEQUALITY CONSTRAINTS

      MODE=7
      K=MAX(LE,N)
      T=SQRT(EPMACH)
      CALL HFTI (W(IE),ME,ME,L,W(IF),K,1,T,KRANK,XNRM,W,W(L+1),JW)
      CALL DCOPY(L,W(IF),1,X(MC1),1)
      IF(KRANK.NE.L)                   GOTO 75
      MODE=1
                                       GOTO 50
C  MODIFY H AND SOLVE INEQUALITY CONSTRAINED LS PROBLEM

   40 DO 45 I=1,MG
   45     H(I)=H(I)-DDOT(MC,G(I,1),LG,X,1)
      CALL LSI
     . (W(IE),W(IF),W(IG),H,ME,ME,MG,MG,L,X(MC1),XNRM,W(MC1),JW,MODE)
      IF(MC.EQ.0)                      GOTO 75
      T=DNRM2(MC,X,1)
      XNRM=SQRT(XNRM*XNRM+T*T)
      IF(MODE.NE.1)                    GOTO 75

C  SOLUTION OF ORIGINAL PROBLEM AND LAGRANGE MULTIPLIERS

   50 DO 55 I=1,ME
   55     F(I)=DDOT(N,E(I,1),LE,X,1)-F(I)
      DO 60 I=1,MC
   60     D(I)=DDOT(ME,E(1,I),1,F,1)-DDOT(MG,G(1,I),1,W(MC1),1)

      DO 65 I=MC,1,-1
   65     CALL H12(2,I,I+1,N,C(I,1),LC,W(IW+I),X,1,1,1)

      DO 70 I=MC,1,-1
          J=MIN(I+1,LC)
          W(I)=(D(I)-DDOT(MC-I,C(J,I),1,W(J),1))/C(I,I)
   70 CONTINUE

C  END OF SUBROUTINE LSEI

   75 END
      