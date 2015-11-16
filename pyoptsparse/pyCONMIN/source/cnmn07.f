      SUBROUTINE CNMN07 (II,XBAR,EPS,X1,Y1,X2,Y2,X3,Y3)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     ROUTINE TO FIND FIRST XBAR.GE.EPS CORRESPONDING TO A REAL ZERO
C     OF A ONE-DIMENSIONAL FUNCTION BY POLYNOMIEL INTERPOLATION.
C     BY G. N. VANDERPLAATS                          APRIL, 1972. 
C     NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF. 
C     II = CALCULATION CONTROL. 
C          1:  2-POINT LINEAR INTERPOLATION, GIVEN X1, Y1, X2 AND Y2. 
C          2:  3-POINT QUADRATIC INTERPOLATION, GIVEN X1, Y1, X2, Y2, 
C              X3 AND Y3. 
C     EPS MAY BE NEGATIVE.
C     IF REQUIRED ZERO ON Y DOES NOT EXITS, OR THE FUNCTION IS
C     ILL-CONDITIONED, XBAR = EPS-1.0 WILL BE RETURNED AS AN ERROR
C     INDICATOR.
C     IF DESIRED INTERPOLATION IS ILL-CONDITIONED, A LOWER ORDER
C     INTERPOLATION, CONSISTANT WITH INPUT DATA, WILL BE ATTEMPTED AND
C     II WILL BE CHANGED ACCORDINGLY. 
      XBAR1=EPS-1.
      XBAR=XBAR1
      JJ=0
      X21=X2-X1 
      IF (ABS(X21).LT.1.0E-20) RETURN 
      IF (II.EQ.2) GO TO 30 
C 
10    CONTINUE
C     ------------------------------------------------------------------
C                  II=1: 2-POINT LINEAR INTERPOLATION 
C     ------------------------------------------------------------------
      II=1
      YY=Y1*Y2
      IF (JJ.EQ.0.OR.YY.LT.0.) GO TO 20 
C     INTERPOLATE BETWEEN X2 AND X3.
      DY=Y3-Y2
      IF (ABS(DY).LT.1.0E-20) GO TO 20
      XBAR=X2+Y2*(X2-X3)/DY 
      IF (XBAR.LT.EPS) XBAR=XBAR1 
      RETURN
20    DY=Y2-Y1
C     INTERPOLATE BETWEEN X1 AND X2.
      IF (ABS(DY).LT.1.0E-20) RETURN
      XBAR=X1+Y1*(X1-X2)/DY 
      IF (XBAR.LT.EPS) XBAR=XBAR1 
      RETURN
30    CONTINUE
C     ------------------------------------------------------------------
C                 II=2: 3-POINT QUADRATIC INTERPOLATION 
C     ------------------------------------------------------------------
      JJ=1
      X31=X3-X1 
      X32=X3-X2 
      QQ=X21*X31*X32
      IF (ABS(QQ).LT.1.0E-20) RETURN
      AA=(Y1*X32-Y2*X31+Y3*X21)/QQ
      IF (ABS(AA).LT.1.0E-20) GO TO 10
      BB=(Y2-Y1)/X21-AA*(X1+X2) 
      CC=Y1-X1*(AA*X1+BB) 
      BAC=BB*BB-4.*AA*CC
      IF (BAC.LT.0.) GO TO 10 
      BAC=SQRT(BAC) 
      AA=.5/AA
      XBAR=AA*(BAC-BB)
      XB2=-AA*(BAC+BB)
      IF (XBAR.LT.EPS) XBAR=XB2 
      IF (XB2.LT.XBAR.AND.XB2.GT.EPS) XBAR=XB2
      IF (XBAR.LT.EPS) XBAR=XBAR1 
      RETURN
      END 