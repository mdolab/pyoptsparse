      SUBROUTINE SLSQP (M,MEQ,LA,N,X,XL,XU,F,C,G,A,ACC,ITER,
     1     IPRINT,IOUT,IFILE,MODE,W,L_W,JW,L_JW,NFUNC,NGRAD,
     2     SLFUNC,SLGRAD)

C   SLSQP       S EQUENTIAL  L EAST  SQ UARES  P ROGRAMMING
C            TO SOLVE GENERAL NONLINEAR OPTIMIZATION PROBLEMS

C***********************************************************************
C*                                                                     *
C*                                                                     *
C*            A NONLINEAR PROGRAMMING METHOD WITH                      *
C*            QUADRATIC  PROGRAMMING  SUBPROBLEMS                      *
C*                                                                     *
C*                                                                     *
C*  THIS SUBROUTINE SOLVES THE GENERAL NONLINEAR PROGRAMMING PROBLEM   *
C*                                                                     *
C*            MINIMIZE    F(X)                                         *
C*                                                                     *
C*            SUBJECT TO  C (X) .EQ. 0  ,  J = 1,...,MEQ               *
C*                         J                                           *
C*                                                                     *
C*                        C (X) .GE. 0  ,  J = MEQ+1,...,M             *
C*                         J                                           *
C*                                                                     *
C*                        XL .LE. X .LE. XU , I = 1,...,N.             *
C*                          I      I       I                           *
C*                                                                     *
C*  THE ALGORITHM IMPLEMENTS THE METHOD OF HAN AND POWELL              *
C*  WITH BFGS-UPDATE OF THE B-MATRIX AND L1-TEST FUNCTION              *
C*  WITHIN THE STEPLENGTH ALGORITHM.                                   *
C*                                                                     *
C*    PARAMETER DESCRIPTION:                                           *
C*    ( * MEANS THIS PARAMETER WILL BE CHANGED DURING CALCULATION )    *
C*                                                                     *
C*    M              IS THE TOTAL NUMBER OF CONSTRAINTS, M .GE. 0      *
C*    MEQ            IS THE NUMBER OF EQUALITY CONSTRAINTS, MEQ .GE. 0 *
C*    LA             SEE A, LA .GE. MAX(M,1)                           *
C*    N              IS THE NUMBER OF VARIBLES, N .GE. 1               *
C*  * X()            X() STORES THE CURRENT ITERATE OF THE N VECTOR X  *
C*                   ON ENTRY X() MUST BE INITIALIZED. ON EXIT X()     *
C*                   STORES THE SOLUTION VECTOR X IF MODE = 0.         *
C*    XL()           XL() STORES AN N VECTOR OF LOWER BOUNDS XL TO X.  *
C*    XU()           XU() STORES AN N VECTOR OF UPPER BOUNDS XU TO X.  *
C*    F              IS THE VALUE OF THE OBJECTIVE FUNCTION.           *
C*    C()            C() STORES THE M VECTOR C OF CONSTRAINTS,         *
C*                   EQUALITY CONSTRAINTS (IF ANY) FIRST.              *
C*                   DIMENSION OF C MUST BE GREATER OR EQUAL LA,       *
C*                   which must be GREATER OR EQUAL MAX(1,M).          *
C*    G()            G() STORES THE N VECTOR G OF PARTIALS OF THE      *
C*                   OBJECTIVE FUNCTION; DIMENSION OF G MUST BE        *
C*                   GREATER OR EQUAL N+1.                             *
C*    A(),LA,M,N     THE LA BY N + 1 ARRAY A() STORES                  *
C*                   THE M BY N MATRIX A OF CONSTRAINT NORMALS.        *
C*                   A() HAS FIRST DIMENSIONING PARAMETER LA,          *
C*                   WHICH MUST BE GREATER OR EQUAL MAX(1,M).          *
C*    F,C,G,A        MUST ALL BE SET BY THE USER BEFORE EACH CALL.     *
C*  * ACC            ABS(ACC) CONTROLS THE FINAL ACCURACY.             *
C*                   IF ACC .LT. ZERO AN EXACT LINESEARCH IS PERFORMED,*
C*                   OTHERWISE AN ARMIJO-TYPE LINESEARCH IS USED.      *
C*  * ITER           PRESCRIBES THE MAXIMUM NUMBER OF ITERATIONS.      *
C*                   ON EXIT ITER INDICATES THE NUMBER OF ITERATIONS.  *
C*  * MODE           MODE CONTROLS CALCULATION:                        *
C*                   REVERSE COMMUNICATION IS USED IN THE SENSE THAT   *
C*                   THE PROGRAM IS INITIALIZED BY MODE = 0; THEN IT IS*
C*                   TO BE CALLED REPEATEDLY BY THE USER UNTIL A RETURN*
C*                   WITH MODE .NE. IABS(1) TAKES PLACE.               *
C*                   IF MODE = -1 GRADIENTS HAVE TO BE CALCULATED,     *
C*                   WHILE WITH MODE = 1 FUNCTIONS HAVE TO BE CALCULATED
C*                   MODE MUST NOT BE CHANGED BETWEEN SUBSEQUENT CALLS *
C*                   OF SQP.                                           *
C*                   EVALUATION MODES:                                 *
C*        MODE = -1: GRADIENT EVALUATION, (G&A)                        *
C*                0: ON ENTRY: INITIALIZATION, (F,G,C&A)               *
C*                   ON EXIT : REQUIRED ACCURACY FOR SOLUTION OBTAINED *
C*                1: FUNCTION EVALUATION, (F&C)                        *
C*                                                                     *
C*                   FAILURE MODES:                                    *
C*                2: NUMBER OF EQUALITY CONTRAINTS LARGER THAN N       *
C*                3: MORE THAN 3*N ITERATIONS IN LSQ SUBPROBLEM        *
C*                4: INEQUALITY CONSTRAINTS INCOMPATIBLE               *
C*                5: SINGULAR MATRIX E IN LSQ SUBPROBLEM               *
C*                6: SINGULAR MATRIX C IN LSQ SUBPROBLEM               *
C*                7: RANK-DEFICIENT EQUALITY CONSTRAINT SUBPROBLEM HFTI*
C*                8: POSITIVE DIRECTIONAL DERIVATIVE FOR LINESEARCH    *
C*                9: MORE THAN ITER ITERATIONS IN SQP                  *
C*             >=10: WORKING SPACE W OR JW TOO SMALL,                  *
C*                   W SHOULD BE ENLARGED TO L_W=MODE/1000             *
C*                   JW SHOULD BE ENLARGED TO L_JW=MODE-1000*L_W       *
C*  * W(), L_W       W() IS A ONE DIMENSIONAL WORKING SPACE,           *
C*                   THE LENGTH L_W OF WHICH SHOULD BE AT LEAST        *
C*                   (3*N1+M)*(N1+1)                        for LSQ    *
C*                  +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ         for LSI    *
C*                  +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1       for LSEI   *
C*                  + N1*N/2 + 2*M + 3*N + 3*N1 + 1         for SLSQPB *
C*                   with MINEQ = M - MEQ + 2*N1  &  N1 = N+1          *
C*        NOTICE:    FOR PROPER DIMENSIONING OF W IT IS RECOMMENDED TO *
C*                   COPY THE FOLLOWING STATEMENTS INTO THE HEAD OF    *
C*                   THE CALLING PROGRAM (AND REMOVE THE COMMENT C)    *
c#######################################################################
C     INTEGER LEN_W, LEN_JW, M, N, N1, MEQ, MINEQ
C     PARAMETER (M=... , MEQ=... , N=...  )
C     PARAMETER (N1= N+1, MINEQ= M-MEQ+N1+N1)
C     PARAMETER (LEN_W=
c    $           (3*N1+M)*(N1+1)
c    $          +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ
c    $          +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1
c    $          +(N+1)*N/2 + 2*M + 3*N + 3*N1 + 1,
c    $           LEN_JW=MINEQ)
C     DOUBLE PRECISION W(LEN_W)
C     INTEGER          JW(LEN_JW)
c#######################################################################
C*                   THE FIRST M+N+N*N1/2 ELEMENTS OF W MUST NOT BE    *
C*                   CHANGED BETWEEN SUBSEQUENT CALLS OF SLSQP.        *
C*                   ON RETURN W(1) ... W(M) CONTAIN THE MULTIPLIERS   *
C*                   ASSOCIATED WITH THE GENERAL CONSTRAINTS, WHILE    *
C*                   W(M+1) ... W(M+N(N+1)/2) STORE THE CHOLESKY FACTOR*
C*                   L*D*L(T) OF THE APPROXIMATE HESSIAN OF THE        *
C*                   LAGRANGIAN COLUMNWISE DENSE AS LOWER TRIANGULAR   *
C*                   UNIT MATRIX L WITH D IN ITS 'DIAGONAL' and        *
C*                   W(M+N(N+1)/2+N+2 ... W(M+N(N+1)/2+N+2+M+2N)       *
C*                   CONTAIN THE MULTIPLIERS ASSOCIATED WITH ALL       *
C*                   ALL CONSTRAINTS OF THE QUADRATIC PROGRAM FINDING  *
C*                   THE SEARCH DIRECTION TO THE SOLUTION X*           *
C*  * JW(), L_JW     JW() IS A ONE DIMENSIONAL INTEGER WORKING SPACE   *
C*                   THE LENGTH L_JW OF WHICH SHOULD BE AT LEAST       *
C*                   MINEQ                                             *
C*                   with MINEQ = M - MEQ + 2*N1  &  N1 = N+1          *
C*                                                                     *
C*  THE USER HAS TO PROVIDE THE FOLLOWING SUBROUTINES:                 *
C*     LDL(N,A,Z,SIG,W) :   UPDATE OF THE LDL'-FACTORIZATION.          *
C*     LINMIN(A,B,F,TOL) :  LINESEARCH ALGORITHM IF EXACT = 1          *
C*     LSQ(M,MEQ,LA,N,NC,C,D,A,B,XL,XU,X,LAMBDA,W,....) :              *
C*                                                                     *
C*        SOLUTION OF THE QUADRATIC PROGRAM                            *
C*                QPSOL IS RECOMMENDED:                                *
C*     PE GILL, W MURRAY, MA SAUNDERS, MH WRIGHT:                      *
C*     USER'S GUIDE FOR SOL/QPSOL:                                     *
C*     A FORTRAN PACKAGE FOR QUADRATIC PROGRAMMING,                    *
C*     TECHNICAL REPORT SOL 83-7, JULY 1983                            *
C*     DEPARTMENT OF OPERATIONS RESEARCH, STANFORD UNIVERSITY          *
C*     STANFORD, CA 94305                                              *
C*     QPSOL IS THE MOST ROBUST AND EFFICIENT QP-SOLVER                *
C*     AS IT ALLOWS WARM STARTS WITH PROPER WORKING SETS               *
C*                                                                     *
C*     IF IT IS NOT AVAILABLE USE LSEI, A CONSTRAINT LINEAR LEAST      *
C*     SQUARES SOLVER IMPLEMENTED USING THE SOFTWARE HFTI, LDP, NNLS   *
C*     FROM C.L. LAWSON, R.J.HANSON: SOLVING LEAST SQUARES PROBLEMS,   *
C*     PRENTICE HALL, ENGLEWOOD CLIFFS, 1974.                          *
C*     LSEI COMES WITH THIS PACKAGE, together with all necessary SR's. *
C*                                                                     *
C*     TOGETHER WITH A COUPLE OF SUBROUTINES FROM BLAS LEVEL 1         *
C*                                                                     *
C*     SQP IS HEAD SUBROUTINE FOR BODY SUBROUTINE SQPBDY               *
C*     IN WHICH THE ALGORITHM HAS BEEN IMPLEMENTED.                    *
C*                                                                     *
C*  IMPLEMENTED BY: DIETER KRAFT, DFVLR OBERPFAFFENHOFEN               *
C*  as described in Dieter Kraft: A Software Package for               *
C*                                Sequential Quadratic Programming     *
C*                                DFVLR-FB 88-28, 1988                 *
C*  which should be referenced if the user publishes results of SLSQP  *
C*                                                                     *
C*  DATE:           APRIL - OCTOBER, 1981.                             *
C*  STATUS:         DECEMBER, 31-ST, 1984.                             *
C*  STATUS:         MARCH   , 21-ST, 1987, REVISED TO FORTAN 77        *
C*  STATUS:         MARCH   , 20-th, 1989, REVISED TO MS-FORTRAN       *
C*  STATUS:         APRIL   , 14-th, 1989, HESSE   in-line coded       *
C*  STATUS:         FEBRUARY, 28-th, 1991, FORTRAN/2 Version 1.04      *
C*                                         accepts Statement Functions *
C*  STATUS:         MARCH   ,  1-st, 1991, tested with SALFORD         *
C*                                         FTN77/386 COMPILER VERS 2.40*
C*                                         in protected mode           *
C*                                                                     *
C***********************************************************************
C*                                                                     *
C*  Copyright 1991: Dieter Kraft, FHM                                  *
C*                                                                     *
C***********************************************************************

      INTEGER IL, IM, IR, IS, ITER, IU, IV, IW, IX, L_W, L_JW,
     1        JW(L_JW), LA, M, MEQ, MINEQ, MODE, N, N1, IPRINT, IOUT,
     2        NFUNC, NGRAD,
     3        IEXACT, INCONS, IRESET, ITERMX, LINE, N2, N3
     
      DOUBLE PRECISION ACC, A(LA,N+1), C(LA), F, G(N+1),
     *     X(N), XL(N), XU(N), W(L_W), 
     *     ALPHA, F0, GS, H1, H2, H3, H4, T, T0, TOL
     
      EXTERNAL SLFUNC,SLGRAD

      CHARACTER*(*) IFILE

c     dim(W) =         N1*(N1+1) + MEQ*(N1+1) + MINEQ*(N1+1)  for LSQ
c                    +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ          for LSI
c                    +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1        for LSEI
c                    + N1*N/2 + 2*M + 3*N +3*N1 + 1           for SLSQPB
c                      with MINEQ = M - MEQ + 2*N1  &  N1 = N+1

C
C   CHECK LENGTH OF WORKING ARRAYS
C
      N1 = N+1
      MINEQ = M-MEQ+N1+N1
      IL = (3*N1+M)*(N1+1) +
     .(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ +
     .(N1+MINEQ)*(N1-MEQ)  + 2*MEQ +
     .N1*N/2 + 2*M + 3*N + 4*N1 + 1
      IM = MAX(MINEQ, N1-MEQ)
      IF (L_W .LT. IL .OR. L_JW .LT. IM) THEN
          MODE = 1000*MAX(10,IL)
          MODE = MODE+MAX(10,IM)
          RETURN
      ENDIF
C
C   PRINT PARAMETERS
C
      IF(IPRINT.LT.0) GOTO 50
      IF(IPRINT.EQ.0) THEN
		  WRITE(*,1000)
          WRITE(*,1100) ACC,ITER,IPRINT,IOUT
      ENDIF
      IF(IPRINT.GT.0) THEN
		  OPEN(UNIT=IOUT,FILE=IFILE(1:LEN_TRIM(IFILE)),
     *      STATUS='UNKNOWN')
		  WRITE(IOUT,1000)
		  WRITE(IOUT,1100) ACC,ITER,IPRINT,IOUT
	  ENDIF
   50 CONTINUE
C
C   INITIAL COUNTERS
C
      NFUNC=0
      NGRAD=0
C
C   CALL OF SQPBDY
C
    4 CONTINUE
C
C   EVALUATE AND PRINT
C
      IF (MODE.EQ.0.OR.MODE.EQ. 1) THEN
          CALL SLFUNC(M, MEQ, LA, N, F, C, X)
          NFUNC=NFUNC+1
      ENDIF
      IF (IPRINT.LT.0) GOTO 60
      IF (MODE.EQ.0.OR.MODE.EQ. 1) GOTO 55
      IF(IPRINT.EQ.0) THEN
		  WRITE(*,1200) ITER,F
	      DO 10 I=1,N
   10 WRITE (*,1400) X(I)
      ENDIF
      IF(IPRINT.GT.0) THEN
		  WRITE(IOUT,1200) ITER,F
	      DO 12 I=1,N
   12 WRITE (IOUT,1400) X(I)
      ENDIF
   55 CONTINUE
   60 CONTINUE
      IF (MODE.EQ.0.OR.MODE.EQ.-1) THEN
         CALL SLGRAD(M, MEQ, LA, N, F, C, G, A, X)
         NGRAD=NGRAD+1
      ENDIF
C
C   PREPARE DATA FOR CALLING SQPBDY  -  INITIAL ADDRESSES IN W
C
      IM = 1
      IL = IM + MAX(1,M)
      IL = IM + LA
      IX = IL + N1*N/2 + 1
      IR = IX + N
      IS = IR + N + N + MAX(1,M)
      IS = IR + N + N + LA
      IU = IS + N1
      IV = IU + N1
      IW = IV + N1
C
      CALL SLSQPB (M,MEQ,LA,N,X,XL,XU,F,C,G,A,ACC,ITER,MODE,
     * W(IR),W(IL),W(IX),W(IM),W(IS),W(IU),W(IV),W(IW),JW,
     * ALPHA,F0,GS,H1,H2,H3,H4,T,T0,TOL,IEXACT,INCONS,IRESET,
     * ITERMX,LINE,N1,N2,N3)
C
      IF (ABS(MODE).EQ.1) GOTO 4
C      
    3 CONTINUE
    
C
C   PRINT FINAL
C
      IF(IPRINT.GT.0) THEN
		  WRITE(IOUT,1450) NFUNC
		  WRITE(IOUT,1460) NGRAD
      ENDIF
C
C   END OF SLSQP
C
      RETURN
C
C     ------------------------------------------------------------------
C                                FORMATS
C     ------------------------------------------------------------------
C  
 1000 FORMAT(////,3X,
     1 60H------------------------------------------------------------,
     2 15H---------------,
     3 /,5X,59HSTART OF THE SEQUENTIAL LEAST SQUARES PROGRAMMING ALGORITHM,
     4     /,3X, 
     5 60H------------------------------------------------------------,
     6 15H---------------)
 1100 FORMAT(/,5X,11HPARAMETERS:,/,8X,5HACC =,D13.4,/,8X,9HMAXITER =,
     1     I3,/,8X,8HIPRINT =,I4,/,6HIOUT =,I4//)
 1200 FORMAT(//5X,6HITER =,I5,5X,5HOBJ =,7E16.8,5X,10HX-VECTOR =)
 1400 FORMAT (3X,7E13.4)
 1450 FORMAT(8X,30HNUMBER OF FUNC-CALLS:  NFUNC =,I4)
 1460 FORMAT(8X,30HNUMBER OF GRAD-CALLS:  NGRAD =,I4)
C
      END
      