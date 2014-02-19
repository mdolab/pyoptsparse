      SUBROUTINE CNMN08 (NDB,NER,C,MS1,B,N3,N4,N5)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION C(N4), B(N3,N3), MS1(N5)
C     ROUTINE TO SOLVE SPECIAL LINEAR PROBLEM FOR IMPOSING S-TRANSPOSE
C     TIMES S.LE.1 BOUNDS IN THE MODIFIED METHOD OF FEASIBLE DIRECTIONS.
C     BY G. N. VANDERPLAATS                             APRIL, 1972.
C     NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF. 
C     REF.  'STRUCTURAL OPTIMIZATION BY METHODS OF FEASIBLE DIRECTIONS',
C     G. N. VANDERPLAATS AND F. MOSES, JOURNAL OF COMPUTERS 
C     AND STRUCTURES, VOL 3, PP 739-755, 1973.
C     FORM OF L. P. IS BX=C WHERE 1ST NDB COMPONENTS OF X CONTAIN VECTOR
C     U AND LAST NDB COMPONENTS CONTAIN VECTOR V.  CONSTRAINTS ARE
C     U.GE.0, V.GE.0, AND U-TRANSPOSE TIMES V = 0.
C     NER = ERROR FLAG.  IF NER.NE.0 ON RETURN, PROCESS HAS NOT 
C     CONVERGED IN 5*NDB ITERATIONS.
C     VECTOR MS1 IDENTIFIES THE SET OF BASIC VARIABLES. 
C     ------------------------------------------------------------------
C     CHOOSE INITIAL BASIC VARIABLES AS V, AND INITIALIZE VECTOR MS1
C     ------------------------------------------------------------------
      NER=1 
      M2=2*NDB
C     CALCULATE CBMIN AND EPS AND INITIALIZE MS1. 
      EPS=-1.0E+10
      CBMIN=0.
      DO 10 I=1,NDB 
      BI=B(I,I) 
      CBMAX=0.
      IF (BI.LT.-1.0E-6) CBMAX=C(I)/BI
      IF (BI.GT.EPS) EPS=BI 
      IF (CBMAX.GT.CBMIN) CBMIN=CBMAX 
10    MS1(I)=0
      EPS=.0001*EPS 
C     IF (EPS.LT.-1.0E-10) EPS=-1.0E-10 
C 
C  E-10 CHANGED TO E-03 ON 1/12/81
C 
      IF (EPS.LT.-1.0E-03) EPS=-1.0E-03 
      IF (EPS.GT.-.0001) EPS=-.0001 
      CBMIN=CBMIN*1.0E-6
C     IF (CBMIN.LT.1.0E-10) CBMIN=1.0E-10 
C 
C  E-10 CHANGED TO E-05 ON 1/12/81
C 
      IF (CBMIN.LT.1.0E-05) CBMIN=1.0E-05 
      ITER1=0 
      NMAX=5*NDB
C     ------------------------------------------------------------------
C     **********             BEGIN NEW ITERATION              **********
C     ------------------------------------------------------------------
20    ITER1=ITER1+1 
      IF (ITER1.GT.NMAX) RETURN 
C     FIND MAX. C(I)/B(I,I) FOR I=1,NDB.
      CBMAX=.9*CBMIN
      ICHK=0
      DO 30 I=1,NDB 
      C1=C(I) 
      BI=B(I,I) 
C     IF (BI.GT.EPS.OR.C1.GT.0.) GO TO 30 
      IF (BI.GT.EPS.OR.C1.GT.-1.0E-05) GO TO 30 
C 
C  0. CHANGED TO -1.0E-05 ON 1/12/81
C 
      CB=C1/BI
      IF (CB.LE.CBMAX) GO TO 30 
      ICHK=I
      CBMAX=CB
30    CONTINUE
      IF (CBMAX.LT.CBMIN) GO TO 70
      IF (ICHK.EQ.0) GO TO 70 
C     UPDATE VECTOR MS1.
      JJ=ICHK 
      IF (MS1(JJ).EQ.0) JJ=ICHK+NDB 
      KK=JJ+NDB 
      IF (KK.GT.M2) KK=JJ-NDB 
      MS1(KK)=ICHK
      MS1(JJ)=0 
C     ------------------------------------------------------------------
C                     PIVOT OF B(ICHK,ICHK) 
C     ------------------------------------------------------------------
      BB=1./B(ICHK,ICHK)
      DO 40 J=1,NDB 
40    B(ICHK,J)=BB*B(ICHK,J)
      C(ICHK)=CBMAX 
      B(ICHK,ICHK)=BB 
C     ELIMINATE COEFICIENTS ON VARIABLE ENTERING BASIS AND STORE
C     COEFICIENTS ON VARIABLE LEAVING BASIS IN THEIR PLACE. 
      DO 60 I=1,NDB 
      IF (I.EQ.ICHK) GO TO 60 
      BB1=B(I,ICHK) 
      B(I,ICHK)=0.
      DO 50 J=1,NDB 
50    B(I,J)=B(I,J)-BB1*B(ICHK,J) 
      C(I)=C(I)-BB1*CBMAX 
60    CONTINUE
      GO TO 20
70    CONTINUE
      NER=0 
C     ------------------------------------------------------------------
C     STORE ONLY COMPONENTS OF U-VECTOR IN 'C'.  USE B(I,1) FOR 
C     TEMPORARY STORAGE 
C     ------------------------------------------------------------------
      DO 80 I=1,NDB 
      B(I,1)=C(I) 
80    CONTINUE
      DO 90 I=1,NDB 
      C(I)=0. 
      J=MS1(I)
      IF (J.GT.0) C(I)=B(J,1) 
      IF (C(I).LT.0.) C(I)=0. 
90    CONTINUE
      RETURN
      END