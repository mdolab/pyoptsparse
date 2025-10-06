      SUBROUTINE CNMN05 (G,DF,A,S,B,C,SLOPE,PHI,ISC,IC,MS1,NVC,N1,N2,N3,
     1N4,N5)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /CNMN1/ DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,CTL,CTLMIN,ALPHAX
     1,ABOBJ1,THETA,OBJ,NDV,NCON,NSIDE,IPRINT,NFDG,NSCAL,LINOBJ,ITMAX,IT
     2RM,ICNDIR,IGOTO,NAC,INFO,INFOG,ITER,NFEASCT
      COMMON /OUTPUT/ IOUT
      DIMENSION DF(N1), G(N2), ISC(N2), IC(N3), A(N1,N3), S(N1), C(N4),
     1MS1(N5), B(N3,N3)
C     ROUTINE TO SOLVE DIRECTION FINDING PROBLEM IN MODIFIED METHOD OF
C     FEASIBLE DIRECTIONS.
C     BY G. N. VANDERPLAATS                            MAY, 1972.
C     NASA-AMES RESEARCH CENTER, MOFFETT FIELD, CALIF.
C     NORM OF S VECTOR USED HERE IS S-TRANSPOSE TIMES S.LE.1.
C     IF NVC = 0 FIND DIRECTION BY ZOUTENDIJK'S METHOD.  OTHERWISE
C     FIND MODIFIED DIRECTION.
C     ------------------------------------------------------------------
C     ***  NORMALIZE GRADIENTS, CALCULATE THETA'S AND DETERMINE NVC  ***
C     ------------------------------------------------------------------
      NDV1=NDV+1
      NDV2=NDV+2
      NAC1=NAC+1
      NVC=0
      THMAX=0.
      CTA=ABS(CT)
      CT1=1./CTA
      CTAM=ABS(CTMIN)
      CTB=ABS(CTL)
      CT2=1./CTB
      CTBM=ABS(CTLMIN)
      A1=1.
      DO 40 I=1,NAC
C     CALCULATE THETA
      NCI=IC(I)
      NCJ=1
      IF (NCI.LE.NCON) NCJ=ISC(NCI)
      C1=G(NCI)
      CTD=CT1
      CTC=CTAM
      IF (NCJ.LE.0) GO TO 10
      CTC=CTBM
      CTD=CT2
10    IF (C1.GT.CTC) NVC=NVC+1
      THT=0.
      GG=1.+CTD*C1
      IF (NCJ.EQ.0.OR.C1.GT.CTC) THT=THETA*GG*GG
      IF (THT.GT.50.) THT=50.
      IF (THT.GT.THMAX) THMAX=THT
      A(NDV1,I)=THT
C     ------------------------------------------------------------------
C                    NORMALIZE GRADIENTS OF CONSTRAINTS
C     ------------------------------------------------------------------
      A(NDV2,I)=1.
      IF (NCI.GT.NCON) GO TO 40
      A1=0.
      DO 20 J=1,NDV
      A1=A1+A(J,I)**2
20    CONTINUE
      IF (A1.LT.1.0E-20) A1=1.0E-20
      A1=SQRT(A1)
      A(NDV2,I)=A1
      A1=1./A1
      DO 30 J=1,NDV
30    A(J,I)=A1*A(J,I)
40    CONTINUE
C     ------------------------------------------------------------------
C     CHECK FOR ZERO GRADIENT.  PROGRAM CHANGE-FEB, 1981, GV.
C     ------------------------------------------------------------------
      I=0
41    I=I+1
42    CONTINUE
      IF(A(NDV2,I).GT.1.0E-6) GO TO 45
C     ZERO GRADIENT IS FOUND.  WRITE ERROR MESSAGE.
      IF(IPRINT.GE.2) WRITE(IOUT,43)IC(I)
43    FORMAT(5X,13H** CONSTRAINT,I5,18H HAS ZERO GRADIENT/
     *5X,23HDELETED FROM ACTIVE SET)
C     REDUCE NAC BY ONE.
      NAC=NAC-1
C     SHIFT COLUMNS OF A AND ROWS OF IC IF I.LE.NAC.
      IF(I.GT.NAC) GO TO 46
C     SHIFT.
      DO 44 J=I,NAC
      J1=J+1
      IC(J)=IC(J1)
      DO 44 K=1,NDV2
44    A(K,J)=A(K,J1)
      IF(I.LE.NAC) GO TO 42
45    CONTINUE
      IF(I.LT.NAC) GO TO 41
46    CONTINUE
      IF(NAC.LE.0) RETURN
      NAC1=NAC+1
C     DETERMINE IF CONSTRAINTS ARE VIOLATED.
      NVC=0
      DO 47 I=1,NAC
      NCI=IC(I)
      NCJ=1
      IF(NCI.LE.NCON) NCJ=ISC(NCI)
      CTC=CTAM
      IF(NCJ.GT.0) CTC=CTBM
      IF(G(NCI).GT.CTC) NVC=NVC+1
47    CONTINUE
C     ------------------------------------------------------------------
C     NORMALIZE GRADIENT OF OBJECTIVE FUNCTION AND STORE IN NAC+1
C     COLUMN OF A
C     ------------------------------------------------------------------
      A1=0.
      DO 50 I=1,NDV
      A1=A1+DF(I)**2
50    CONTINUE
      IF (A1.LT.1.0E-20) A1=1.0E-20
      A1=SQRT(A1)
      A1=1./A1
      DO 60 I=1,NDV
60    A(I,NAC1)=A1*DF(I)
C     BUILD C VECTOR.
      IF (NVC.GT.0) GO TO 80
C     ------------------------------------------------------------------
C                 BUILD C FOR CLASSICAL METHOD
C     ------------------------------------------------------------------
      NDB=NAC1
      A(NDV1,NDB)=1.
      DO 70 I=1,NDB
70    C(I)=-A(NDV1,I)
      GO TO 110
80    CONTINUE
C     ------------------------------------------------------------------
C                   BUILD C FOR MODIFIED METHOD
C     ------------------------------------------------------------------
      NDB=NAC
      A(NDV1,NAC1)=-PHI
C     ------------------------------------------------------------------
C           SCALE THETA'S SO THAT MAXIMUM THETA IS UNITY
C     ------------------------------------------------------------------
      IF (THMAX.GT.0.00001) THMAX=1./THMAX
      DO 90 I=1,NDB
      A(NDV1,I)=A(NDV1,I)*THMAX
90    CONTINUE
      DO 100 I=1,NDB
      C(I)=0.
      DO 100 J=1,NDV1
100   C(I)=C(I)+A(J,I)*A(J,NAC1)
110   CONTINUE
C     ------------------------------------------------------------------
C                      BUILD B MATRIX
C     ------------------------------------------------------------------
      DO 120 I=1,NDB
      DO 120 J=1,NDB
      B(I,J)=0.
      DO 120 K=1,NDV1
120   B(I,J)=B(I,J)-A(K,I)*A(K,J)
C     ------------------------------------------------------------------
C                    SOLVE SPECIAL L. P. PROBLEM
C     ------------------------------------------------------------------
      CALL CNMN08 (NDB,NER,C,MS1,B,N3,N4,N5)
      IF (IPRINT.GT.1.AND.NER.GT.0) WRITE (IOUT,180)
C     CALCULATE RESULTING DIRECTION VECTOR, S.
      SLOPE=0.
C     ------------------------------------------------------------------
C                  USABLE-FEASIBLE DIRECTION
C     ------------------------------------------------------------------
      DO 140 I=1,NDV
      S1=0.
      IF (NVC.GT.0) S1=-A(I,NAC1)
      DO 130 J=1,NDB
130   S1=S1-A(I,J)*C(J)
      SLOPE=SLOPE+S1*DF(I)
140   S(I)=S1
      S(NDV1)=1.
      IF (NVC.GT.0) S(NDV1)=-A(NDV1,NAC1)
      DO 150 J=1,NDB
150   S(NDV1)=S(NDV1)-A(NDV1,J)*C(J)
C     ------------------------------------------------------------------
C     CHECK TO INSURE THE S-VECTOR IS FEASIBLE.
C     PROGRAM MOD-FEB, 1981, GV.
C     ------------------------------------------------------------------
      DO 174 J=1,NAC
C     S DOT DEL(G).
      SG=0.
      DO 172 I=1,NDV
172   SG=SG+S(I)*A(I,J)
C     IF(SG.GT.0.) GO TO 176
C
C  THIS CHANGE MADE ON 4/8/81 FOR G. VANDERPLAATS
C
      IF(SG.GT.1.0E-04) GO TO 176
C     FEASIBLE FOR THIS CONSTRAINT.  CONTINUE.
174   CONTINUE
      GO TO 179
176   CONTINUE
C     S-VECTOR IS NOT FEASIBLE DUE TO SOME NUMERICAL PROBLEM.
      IF(IPRINT.GE.2) WRITE(IOUT,178)
178   FORMAT(5X,38H** CALCULATED S-VECTOR IS NOT FEASIBLE/5X,
     * 19HBETA IS SET TO ZERO)
      S(NDV1)=0.
      NVC=0
      RETURN
179   CONTINUE
C     ------------------------------------------------------------------
C                  NORMALIZE S TO MAX ABS OF UNITY
C     ------------------------------------------------------------------
      S1=0.
      DO 160 I=1,NDV
      A1=ABS(S(I))
      IF (A1.GT.S1) S1=A1
160   CONTINUE
C     IF (S1.LT.1.0E-10) RETURN
C
C  E-10 CHANGED TO E-04 ON 1/12/81
C
      IF (S1.LT.1.0E-04) RETURN
      S1=1./S1
      DO 170 I=1,NDV
170   S(I)=S1*S(I)
      SLOPE=S1*SLOPE
      S(NDV1)=S1*S(NDV1)
      RETURN
C     ------------------------------------------------------------------
C                           FORMATS
C     ------------------------------------------------------------------
C
C
180   FORMAT (//5X,46H* * DIRECTION FINDING PROCESS DID NOT CONVERGE/5X,
     129H* * S-VECTOR MAY NOT BE VALID)
      END
