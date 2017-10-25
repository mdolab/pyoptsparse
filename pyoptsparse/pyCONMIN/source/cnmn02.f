      SUBROUTINE CNMN02 (NCALC,SLOPE,DFTDF1,DF,S,N1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /CNMN1/ DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,CTL,CTLMIN,ALPHAX
     1,ABOBJ1,THETA,OBJ,NDV,NCON,NSIDE,IPRINT,NFDG,NSCAL,LINOBJ,ITMAX,IT
     2RM,ICNDIR,IGOTO,NAC,INFO,INFOG,ITER,NFEASCT 
      DIMENSION DF(N1), S(N1) 
C     ROUTINE TO DETERMINE CONJUGATE DIRECTION VECTOR OR DIRECTION
C     OF STEEPEST DESCENT FOR UNCONSTRAINED FUNCTION MINIMIZATION.
C     BY G. N. VANDERPLAATS                       APRIL, 1972.
C     NASA-AMES RESEARCH CENTER, MOFFETT FIELD, CALIF.
C     NCALC = CALCULATION CONTROL.
C         NCALC = 0,     S = STEEPEST DESCENT.
C         NCALC = 1,     S = CONJUGATE DIRECTION. 
C     CONJUGATE DIRECTION IS FOUND BY FLETCHER-REEVES ALGORITHM.
C     ------------------------------------------------------------------
C                   CALCULATE NORM OF GRADIENT VECTOR 
C     ------------------------------------------------------------------
      DFTDF=0.
      DO 10 I=1,NDV 
      DFI=DF(I) 
10    DFTDF=DFTDF+DFI*DFI 
C     ------------------------------------------------------------------
C     **********                FIND DIRECTION S              **********
C     ------------------------------------------------------------------
      IF (NCALC.NE.1) GO TO 30
      IF (DFTDF1.LT.1.0E-20) GO TO 30 
C     ------------------------------------------------------------------
C                 FIND FLETCHER-REEVES CONJUGATE DIRECTION
C     ------------------------------------------------------------------
      BETA=DFTDF/DFTDF1 
      SLOPE=0.
      DO 20 I=1,NDV 
      DFI=DF(I) 
      SI=BETA*S(I)-DFI
      SLOPE=SLOPE+SI*DFI
20    S(I)=SI 
      GO TO 50
30    CONTINUE
      NCALC=0 
C     ------------------------------------------------------------------
C                  CALCULATE DIRECTION OF STEEPEST DESCENT
C     ------------------------------------------------------------------
      DO 40 I=1,NDV 
40    S(I)=-DF(I) 
      SLOPE=-DFTDF
50    CONTINUE
C     ------------------------------------------------------------------
C                  NORMALIZE S TO MAX ABS VALUE OF UNITY
C     ------------------------------------------------------------------
      S1=0. 
      DO 60 I=1,NDV 
      S2=ABS(S(I))
      IF (S2.GT.S1) S1=S2 
60    CONTINUE
      IF (S1.LT.1.0E-20) S1=1.0E-20 
      S1=1./S1
      DFTDF1=DFTDF*S1 
      DO 70 I=1,NDV 
70    S(I)=S1*S(I)
      SLOPE=S1*SLOPE
      RETURN
      END 