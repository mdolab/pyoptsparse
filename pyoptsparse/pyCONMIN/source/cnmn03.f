      SUBROUTINE CNMN03 (X,S,SLOPE,ALP,FFF,A1,A2,A3,A4,F1,F2,F3,F4,APP,N
     11,NCAL,KOUNT,JGOTO)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /CNMN1/ DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,CTL,CTLMIN,ALPHAX
     1,ABOBJ1,THETA,OBJ,NDV,NCON,NSIDE,IPRINT,NFDG,NSCAL,LINOBJ,ITMAX,IT
     2RM,ICNDIR,IGOTO,NAC,INFO,INFOG,ITER,NFEASCT
      COMMON /OUTPUT/ IOUT
      DIMENSION X(N1), S(N1), NCAL(2)
C     ROUTINE TO SOLVE ONE-DIMENSIONAL SEARCH IN UNCONSTRAINED
C     MINIMIZATION USING 2-POINT QUADRATIC INTERPOLATION, 3-POINT
C     CUBIC INTERPOLATION AND 4-POINT CUBIC INTERPOLATION.
C     BY G. N. VANDERPLAATS                         APRIL, 1972.
C     NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF.
C     ALP = PROPOSED MOVE PARAMETER.
C     SLOPE = INITIAL FUNCTION SLOPE = S-TRANSPOSE TIMES DF.
C     SLOPE MUST BE NEGATIVE.
C     OBJ = INITIAL FUNCTION VALUE.
      ZRO=0.
      IF (JGOTO.EQ.0) GO TO 10
      GO TO (50,80,110,140,180,220,270),JGOTO
C     ------------------------------------------------------------------
C                     INITIAL INFORMATION  (ALPHA=0)
C     ------------------------------------------------------------------
10    IF (SLOPE.LT.0.) GO TO 20
      ALP=0.
      RETURN
20    CONTINUE
      IF (IPRINT.GT.4) WRITE (IOUT,360)
      FFF=OBJ
      AP1=0.
      A1=0.
      F1=OBJ
      A2=ALP
      A3=0.
      F3=0.
      AP=A2
      KOUNT=0
C     ------------------------------------------------------------------
C            MOVE A DISTANCE AP*S AND UPDATE FUNCTION VALUE
C     ------------------------------------------------------------------
30    CONTINUE
      KOUNT=KOUNT+1
      DO 40 I=1,NDV
40    X(I)=X(I)+AP*S(I)
      IF (IPRINT.GT.4) WRITE (IOUT,370) AP
      IF (IPRINT.GT.4) WRITE (IOUT,380) (X(I),I=1,NDV)
      NCAL(1)=NCAL(1)+1
      JGOTO=1
      RETURN
50    CONTINUE
      F2=OBJ
      IF (IPRINT.GT.4) WRITE (IOUT,390) F2
      IF (F2.LT.F1) GO TO 120
C     ------------------------------------------------------------------
C                     CHECK FOR ILL-CONDITIONING
C     ------------------------------------------------------------------
      IF (KOUNT.GT.5) GO TO 60
      FF=2.*ABS(F1)
      IF (F2.LT.FF) GO TO 90
      FF=5.*ABS(F1)
      IF (F2.LT.FF) GO TO 60
      A2=.5*A2
      AP=-A2
      ALP=A2
      GO TO 30
60    F3=F2
      A3=A2
      A2=.5*A2
C     ------------------------------------------------------------------
C                 UPDATE DESIGN VECTOR AND FUNCTION VALUE
C     ------------------------------------------------------------------
      AP=A2-ALP
      ALP=A2
      DO 70 I=1,NDV
70    X(I)=X(I)+AP*S(I)
      IF (IPRINT.GT.4) WRITE (IOUT,370) A2
      IF (IPRINT.GT.4) WRITE (IOUT,380) (X(I),I=1,NDV)
      NCAL(1)=NCAL(1)+1
      JGOTO=2
      RETURN
80    CONTINUE
      F2=OBJ
      IF (IPRINT.GT.4) WRITE (IOUT,390) F2
C     PROCEED TO CUBIC INTERPOLATION.
      GO TO 160
90    CONTINUE
C     ------------------------------------------------------------------
C     **********        2-POINT QUADRATIC INTERPOLATION       **********
C     ------------------------------------------------------------------
      JJ=1
      II=1
      CALL CNMN04 (II,APP,ZRO,A1,F1,SLOPE,A2,F2,ZRO,ZRO,ZRO,ZRO)
      IF (APP.LT.ZRO.OR.APP.GT.A2) GO TO 120
      F3=F2
      A3=A2
      A2=APP
      JJ=0
C     ------------------------------------------------------------------
C                  UPDATE DESIGN VECTOR AND FUNCTION VALUE
C     ------------------------------------------------------------------
      AP=A2-ALP
      ALP=A2
      DO 100 I=1,NDV
100   X(I)=X(I)+AP*S(I)
      IF (IPRINT.GT.4) WRITE (IOUT,370) A2
      IF (IPRINT.GT.4) WRITE (IOUT,380) (X(I),I=1,NDV)
      NCAL(1)=NCAL(1)+1
      JGOTO=3
      RETURN
110   CONTINUE
      F2=OBJ
      IF (IPRINT.GT.4) WRITE (IOUT,390) F2
      GO TO 150
120   A3=2.*A2
C     ------------------------------------------------------------------
C                  UPDATE DESIGN VECTOR AND FUNCTION VALUE
C     ------------------------------------------------------------------
      AP=A3-ALP
      ALP=A3
      DO 130 I=1,NDV
130   X(I)=X(I)+AP*S(I)
      IF (IPRINT.GT.4) WRITE (IOUT,370) A3
      IF (IPRINT.GT.4) WRITE (IOUT,380) (X(I),I=1,NDV)
      NCAL(1)=NCAL(1)+1
      JGOTO=4
      RETURN
140   CONTINUE
      F3=OBJ
      IF (IPRINT.GT.4) WRITE (IOUT,390) F3
150   CONTINUE
      IF (F3.LT.F2) GO TO 190
160   CONTINUE
C     ------------------------------------------------------------------
C     **********       3-POINT CUBIC INTERPOLATION      **********
C     ------------------------------------------------------------------
      II=3
      CALL CNMN04 (II,APP,ZRO,A1,F1,SLOPE,A2,F2,A3,F3,ZRO,ZRO)
      IF (APP.LT.ZRO.OR.APP.GT.A3) GO TO 190
C     ------------------------------------------------------------------
C     UPDATE DESIGN VECTOR AND FUNCTION VALUE.
C     ------------------------------------------------------------------
      AP1=APP
      AP=APP-ALP
      ALP=APP
      DO 170 I=1,NDV
170   X(I)=X(I)+AP*S(I)
      IF (IPRINT.GT.4) WRITE (IOUT,370) ALP
      IF (IPRINT.GT.4) WRITE (IOUT,380) (X(I),I=1,NDV)
      NCAL(1)=NCAL(1)+1
      JGOTO=5
      RETURN
180   CONTINUE
      IF (IPRINT.GT.4) WRITE (IOUT,390) OBJ
C     ------------------------------------------------------------------
C                         CHECK CONVERGENCE
C     ------------------------------------------------------------------
      AA=1.-APP/A2
      AB2=ABS(F2)
      AB3=ABS(OBJ)
      AB=AB2
      IF (AB3.GT.AB) AB=AB3
      IF (AB.LT.1.0E-15) AB=1.0E-15
      AB=(AB2-AB3)/AB
      IF (ABS(AB).LT.1.0E-15.AND.ABS(AA).LT..001) GO TO 330
      A4=A3
      F4=F3
      A3=APP
      F3=OBJ
      IF (A3.GT.A2) GO TO 230
      A3=A2
      F3=F2
      A2=APP
      F2=OBJ
      GO TO 230
190   CONTINUE
C     ------------------------------------------------------------------
C     **********        4-POINT CUBIC INTERPOLATION       **********
C     ------------------------------------------------------------------
200   CONTINUE
      A4=2.*A3
C     UPDATE DESIGN VECTOR AND FUNCTION VALUE.
      AP=A4-ALP
      ALP=A4
      DO 210 I=1,NDV
210   X(I)=X(I)+AP*S(I)
      IF (IPRINT.GT.4) WRITE (IOUT,370) ALP
      IF (IPRINT.GT.4) WRITE (IOUT,380) (X(I),I=1,NDV)
      NCAL(1)=NCAL(1)+1
      JGOTO=6
      RETURN
220   CONTINUE
      F4=OBJ
      IF (IPRINT.GT.4) WRITE (IOUT,390) F4
      IF (F4.GT.F3) GO TO 230
      A1=A2
      F1=F2
      A2=A3
      F2=F3
      A3=A4
      F3=F4
      GO TO 200
230   CONTINUE
      II=4
      CALL CNMN04 (II,APP,A1,A1,F1,SLOPE,A2,F2,A3,F3,A4,F4)
      IF (APP.GT.A1) GO TO 250
      AP=A1-ALP
      ALP=A1
      OBJ=F1
      DO 240 I=1,NDV
240   X(I)=X(I)+AP*S(I)
      GO TO 280
250   CONTINUE
C     ------------------------------------------------------------------
C                 UPDATE DESIGN VECTOR AND FUNCTION VALUE
C     ------------------------------------------------------------------
      AP=APP-ALP
      ALP=APP
      DO 260 I=1,NDV
260   X(I)=X(I)+AP*S(I)
      IF (IPRINT.GT.4) WRITE (IOUT,370) ALP
      IF (IPRINT.GT.4) WRITE (IOUT,380) (X(I),I=1,NDV)
      NCAL(1)=NCAL(1)+1
      JGOTO=7
      RETURN
270   CONTINUE
      IF (IPRINT.GT.4) WRITE (IOUT,390) OBJ
280   CONTINUE
C     ------------------------------------------------------------------
C                    CHECK FOR ILL-CONDITIONING
C     ------------------------------------------------------------------
      IF (OBJ.GT.F2.OR.OBJ.GT.F3) GO TO 290
      IF (OBJ.LE.F1) GO TO 330
      AP=A1-ALP
      ALP=A1
      OBJ=F1
      GO TO 310
290   CONTINUE
      IF (F2.LT.F3) GO TO 300
      OBJ=F3
      AP=A3-ALP
      ALP=A3
      GO TO 310
300   OBJ=F2
      AP=A2-ALP
      ALP=A2
310   CONTINUE
C     ------------------------------------------------------------------
C                       UPDATE DESIGN VECTOR
C     ------------------------------------------------------------------
      DO 320 I=1,NDV
320   X(I)=X(I)+AP*S(I)
330   CONTINUE
C     ------------------------------------------------------------------
C                     CHECK FOR MULTIPLE MINIMA
C     ------------------------------------------------------------------
      IF (OBJ.LE.FFF) GO TO 350
C     INITIAL FUNCTION IS MINIMUM.
      DO 340 I=1,NDV
340   X(I)=X(I)-ALP*S(I)
      ALP=0.
      OBJ=FFF
350   CONTINUE
      JGOTO=0
      RETURN
C     ------------------------------------------------------------------
C                                 FORMATS
C     ------------------------------------------------------------------
C
C
360   FORMAT (/////5X,60H* * * UNCONSTRAINED ONE-DIMENSIONAL SEARCH INFO
     1RMATION * * *)
370   FORMAT (/5X,7HALPHA =,E14.5/5X,8HX-VECTOR)
380   FORMAT (5X,6E13.5)
390   FORMAT (/5X,5HOBJ =,E14.5)
      END
