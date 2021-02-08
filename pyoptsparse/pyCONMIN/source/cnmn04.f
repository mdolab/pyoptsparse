      SUBROUTINE CNMN04 (II,XBAR,EPS,X1,Y1,SLOPE,X2,Y2,X3,Y3,X4,Y4)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     ROUTINE TO FIND FIRST XBAR.GE.EPS CORRESPONDING TO A MINIMUM
C     OF A ONE-DIMENSIONAL REAL FUNCTION BY POLYNOMIEL INTERPOLATION. 
C     BY G. N. VANDERPLAATS                          APRIL, 1972. 
C     NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF. 
C 
C     II = CALCULATION CONTROL. 
C          1:  2-POINT QUADRATIC INTERPOLATION, GIVEN X1, Y1, SLOPE,
C              X2 AND Y2. 
C          2:  3-POINT QUADRATIC INTERPOLATION, GIVEN X1, Y1, X2, Y2, 
C              X3 AND Y3. 
C          3:  3-POINT CUBIC INTERPOLATION, GIVEN X1, Y1, SLOPE, X2, Y2,
C              X3 AND Y3. 
C          4:  4-POINT CUBIC INTERPOLATION, GIVEN X1, Y1, X2, Y2, X3, 
C              Y3, X4 AND Y4. 
C     EPS MAY BE NEGATIVE.
C     IF REQUIRED MINIMUM ON Y DOES NOT EXITS, OR THE FUNCTION IS 
C     ILL-CONDITIONED, XBAR = EPS-1.0 WILL BE RETURNED AS AN ERROR
C     INDICATOR.
C     IF DESIRED INTERPOLATION IS ILL-CONDITIONED, A LOWER ORDER
C     INTERPOLATION, CONSISTANT WITH INPUT DATA, WILL BE ATTEMPTED, 
C     AND II WILL BE CHANGED ACCORDINGLY. 
      XBAR1=EPS-1.
      XBAR=XBAR1
      X21=X2-X1 
      IF (ABS(X21).LT.1.0E-20) RETURN 
      NSLOP=MOD(II,2) 
      GO TO (10,20,40,50),II
10    CONTINUE
C     ------------------------------------------------------------------
C                 II=1: 2-POINT QUADRATIC INTERPOLATION 
C     ------------------------------------------------------------------
      II=1
      DX=X1-X2
      IF (ABS(DX).LT.1.0E-20) RETURN
      AA=(SLOPE+(Y2-Y1)/DX)/DX
      IF (AA.LT.1.0E-20) RETURN 
      BB=SLOPE-2.*AA*X1 
      XBAR=-.5*BB/AA
      IF (XBAR.LT.EPS) XBAR=XBAR1 
      RETURN
20    CONTINUE
C     ------------------------------------------------------------------
C                 II=2: 3-POINT QUADRATIC INTERPOLATION 
C     ------------------------------------------------------------------
      II=2
      X21=X2-X1 
      X31=X3-X1 
      X32=X3-X2 
      QQ=X21*X31*X32
      IF (ABS(QQ).LT.1.0E-20) RETURN
      AA=(Y1*X32-Y2*X31+Y3*X21)/QQ
      IF (AA.LT.1.0E-20) GO TO 30 
      BB=(Y2-Y1)/X21-AA*(X1+X2) 
      XBAR=-.5*BB/AA
      IF (XBAR.LT.EPS) XBAR=XBAR1 
      RETURN
30    CONTINUE
      IF (NSLOP.EQ.0) RETURN
      GO TO 10
40    CONTINUE
C     ------------------------------------------------------------------
C                   II=3: 3-POINT CUBIC INTERPOLATION 
C     ------------------------------------------------------------------
      II=3
      X21=X2-X1 
      X31=X3-X1 
      X32=X3-X2 
      QQ=X21*X31*X32
      IF (ABS(QQ).LT.1.0E-20) RETURN
      X11=X1*X1 
      DNOM=X2*X2*X31-X11*X32-X3*X3*X21
      IF (ABS(DNOM).LT.1.0E-20) GO TO 20
      AA=((X31*X31*(Y2-Y1)-X21*X21*(Y3-Y1))/(X31*X21)-SLOPE*X32)/DNOM 
      IF (ABS(AA).LT.1.0E-20) GO TO 20
      BB=((Y2-Y1)/X21-SLOPE-AA*(X2*X2+X1*X2-2.*X11))/X21
      CC=SLOPE-3.*AA*X11-2.*BB*X1 
      BAC=BB*BB-3.*AA*CC
      IF (BAC.LT.0.) GO TO 20 
      BAC=SQRT(BAC) 
      XBAR=(BAC-BB)/(3.*AA) 
      IF (XBAR.LT.EPS) XBAR=EPS 
      RETURN
50    CONTINUE
C     ------------------------------------------------------------------
C                    II=4: 4-POINT CUBIC INTERPOLATION
C     ------------------------------------------------------------------
      X21=X2-X1 
      X31=X3-X1 
      X41=X4-X1 
      X32=X3-X2 
      X42=X4-X2 
      X11=X1*X1 
      X22=X2*X2 
      X33=X3*X3 
      X44=X4*X4 
      X111=X1*X11 
      X222=X2*X22 
      Q2=X31*X21*X32
      IF (ABS(Q2).LT.1.0E-30) RETURN
      Q1=X111*X32-X222*X31+X3*X33*X21 
      Q4=X111*X42-X222*X41+X4*X44*X21 
      Q5=X41*X21*X42
      DNOM=Q2*Q4-Q1*Q5
      IF (ABS(DNOM).LT.1.0E-30) GO TO 60
      Q3=Y3*X21-Y2*X31+Y1*X32 
      Q6=Y4*X21-Y2*X41+Y1*X42 
      AA=(Q2*Q6-Q3*Q5)/DNOM 
      BB=(Q3-Q1*AA)/Q2
      CC=(Y2-Y1-AA*(X222-X111))/X21-BB*(X1+X2)
      BAC=BB*BB-3.*AA*CC
      IF (ABS(AA).LT.1.0E-20.OR.BAC.LT.0.) GO TO 60 
      BAC=SQRT(BAC) 
      XBAR=(BAC-BB)/(3.*AA) 
      IF (XBAR.LT.EPS) XBAR=XBAR1 
      RETURN
60    CONTINUE
      IF (NSLOP.EQ.1) GO TO 40
      GO TO 20
      END 