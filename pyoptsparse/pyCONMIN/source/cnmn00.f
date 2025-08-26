      SUBROUTINE CNMN00 (X,VLB,VUB,G,SCAL,DF,A,S,G1,G2,B,C,ISC,IC,MS1,N1
     1,N2,N3,N4,N5)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /CNMN1/ DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,CTL,CTLMIN,ALPHAX
     1,ABOBJ1,THETA,OBJ,NDV,NCON,NSIDE,IPRINT,NFDG,NSCAL,LINOBJ,ITMAX,IT
     2RM,ICNDIR,IGOTO,NAC,INFO,INFOG,ITER,NFEASCT
      COMMON /OUTPUT/ IOUT
C
C    NFEASCT ADDED TO COMMON BLOCK BY KCYOUNG ON 4/14/92 TO ALLOW MORE
C    THAN 10 ITERATION ATTEMPTS.  NFEASCT BECOMES AN INPUT VALUE
C
      DIMENSION X(N1), VLB(N1), VUB(N1), G(N2), SCAL(N1), DF(N1), A(N1,N
     13), S(N1), G1(N2), G2(N2), B(N3,N3), C(N4), ISC(N2), IC(N3), MS1(N
     25)
      COMMON /CONSAV/ DM1,DM2,DM3,DM4,DM5,DM6,DM7,DM8,DM9,DM10,DM11,DM12
     1,DCT,DCTL,PHI,ABOBJ,CTA,CTAM,CTBM,OBJ1,SLOPE,DX,DX1,FI,XI,DFTDF1,A
     2LP,FFF,A1,A2,A3,A4,F1,F2,F3,F4,CV1,CV2,CV3,CV4,APP,ALPCA,ALPFES,AL
     3PLN,ALPMIN,ALPNC,ALPSAV,ALPSID,ALPTOT,RSPACE,IDM1,IDM2,IDM3,JDIR,I
     4OBJ,KOBJ,KCOUNT,NCAL(2),NFEAS,MSCAL,NCOBJ,NVC,KOUNT,ICOUNT,IGOOD1,
     5IGOOD2,IGOOD3,IGOOD4,IBEST,III,NLNC,JGOTO,ISPACE(2)
C     ROUTINE TO SOLVE CONSTRAINED OR UNCONSTRAINED FUNCTION
C     MINIMIZATION.
C     BY G. N. VANDERPLAATS                          APRIL, 1972.
C     * * * * * * * * * * *   JUNE, 1979 VERSION   * * * * * * * * * * *
C     NASA-AMES RESEARCH CENTER, MOFFETT FIELD, CALIF.
C     REFERENCE;  CONMIN - A FORTRAN PROGRAM FOR CONSTRAINED FUNCTION
C         MINIMIZATION:  USER'S MANUAL,  BY G. N. VANDERPLAATS,
C         NASA TM X-62,282, AUGUST, 1973.
C     STORAGE REQUIREMENTS:
C         PROGRAM - 7000 DECIMAL WORDS (CDC COMPUTER)
C         ARRAYS  - APPROX. 2*(NDV**2)+26*NDV+4*NCON,
C               WHERE N3 = NDV+2.
C     RE-SCALE VARIABLES IF REQUIRED.
      IF (NSCAL.EQ.0.OR.IGOTO.EQ.0) GO TO 20
      DO 10 I=1,NDV
10    X(I)=C(I)
20    CONTINUE
C     CONSTANTS.
      NDV1=NDV+1
      NDV2=NDV+2
      IF (IGOTO.EQ.0) GO TO 40
C     ------------------------------------------------------------------
C                     CHECK FOR UNBOUNDED SOLUTION
C     ------------------------------------------------------------------
C     STOP IF OBJ IS LESS THAN -1.0E+20
      IF (OBJ.GT.-1.0E+20) GO TO 30
      WRITE (IOUT,980)
      GO TO 810
30    CONTINUE
      GO TO (160,390,380,670,690),IGOTO
C     ------------------------------------------------------------------
C                      SAVE INPUT CONTROL PARAMETERS
C     ------------------------------------------------------------------
40    CONTINUE
      IF (IPRINT.GT.0) WRITE (IOUT,1220)
      IF (LINOBJ.EQ.0.OR.(NCON.GT.0.OR.NSIDE.GT.0)) GO TO 50
C     TOTALLY UNCONSTRAINED FUNCTION WITH LINEAR OBJECTIVE.
C     SOLUTION IS UNBOUNDED.
      WRITE (IOUT,970) LINOBJ,NCON,NSIDE
      RETURN
50    CONTINUE
      IDM1=ITRM
      IDM2=ITMAX
      IDM3=ICNDIR
      DM1=DELFUN
      DM2=DABFUN
      DM3=CT
      DM4=CTMIN
      DM5=CTL
      DM6=CTLMIN
      DM7=THETA
      DM8=PHI
      DM9=FDCH
      DM10=FDCHM
      DM11=ABOBJ1
      DM12=ALPHAX
C     ------------------------------------------------------------------
C                                DEFAULTS
C     ------------------------------------------------------------------
      IF (ITRM.LE.0) ITRM=3
      IF (ITMAX.LE.0) ITMAX=20
      NDV1=NDV+1
      IF (ICNDIR.EQ.0) ICNDIR=NDV1
      IF (DELFUN.LE.0.) DELFUN=.0001
      CT=-ABS(CT)
      IF (CT.GE.0.) CT=-.1
      CTMIN=ABS(CTMIN)
      IF (CTMIN.LE.0.) CTMIN=.004
      CTL=-ABS(CTL)
      IF (CTL.GE.0.) CTL=-0.01
      CTLMIN=ABS(CTLMIN)
      IF (CTLMIN.LE.0.) CTLMIN=.001
      IF (THETA.LE.0.) THETA=1.
      IF (ABOBJ1.LE.0.) ABOBJ1=.1
      IF (ALPHAX.LE.0.) ALPHAX=.1
      IF (FDCH.LE.0.) FDCH=.01
      IF (FDCHM.LE.0.) FDCHM=.01
C     ------------------------------------------------------------------
C                     INITIALIZE INTERNAL PARAMETERS
C     ------------------------------------------------------------------
      INFOG=0
      ITER=0
      JDIR=0
      IOBJ=0
      KOBJ=0
      NDV2=NDV+2
      KCOUNT=0
      NCAL(1)=0
      NCAL(2)=0
      NAC=0
      NFEAS=0
      MSCAL=NSCAL
      CT1=ITRM
      CT1=1./CT1
      DCT=(CTMIN/ABS(CT))**CT1
      DCTL=(CTLMIN/ABS(CTL))**CT1
      PHI=5.
      ABOBJ=ABOBJ1
      NCOBJ=0
      CTAM=ABS(CTMIN)
      CTBM=ABS(CTLMIN)
C     CALCULATE NUMBER OF LINEAR CONSTRAINTS, NLNC.
      NLNC=0
      IF (NCON.EQ.0) GO TO 70
      DO 60 I=1,NCON
      IF (ISC(I).GT.0) NLNC=NLNC+1
60    CONTINUE
70    CONTINUE
C     ------------------------------------------------------------------
C          CHECK TO BE SURE THAT SIDE CONSTRAINTS ARE SATISFIED
C     ------------------------------------------------------------------
      IF (NSIDE.EQ.0) GO TO 110
      DO 100 I=1,NDV
      IF (VLB(I).LE.VUB(I)) GO TO 80
      XX=.5*(VLB(I)+VUB(I))
      X(I)=XX
      VLB(I)=XX
      VUB(I)=XX
      WRITE (IOUT,1120) I
80    CONTINUE
      XX=X(I)-VLB(I)
      IF (XX.GE.0.) GO TO 90
C     LOWER BOUND VIOLATED.
      WRITE (IOUT,1130) X(I),VLB(I),I
      X(I)=VLB(I)
      GO TO 100
90    CONTINUE
      XX=VUB(I)-X(I)
      IF (XX.GE.0.) GO TO 100
      WRITE (IOUT,1140) X(I),VUB(I),I
      X(I)=VUB(I)
100   CONTINUE
110   CONTINUE
C     ------------------------------------------------------------------
C                        INITIALIZE SCALING VECTOR, SCAL
C     ------------------------------------------------------------------
      IF (NSCAL.EQ.0) GO TO 150
      IF (NSCAL.LT.0) GO TO 130
      DO 120 I=1,NDV
120   SCAL(I)=1.
      GO TO 150
130   CONTINUE
      DO 140 I=1,NDV
      SI=ABS(SCAL(I))
      IF (SI.LT.1.0E-20) SI=1.0E-5
      SCAL(I)=SI
      SI=1./SI
      X(I)=X(I)*SI
      IF (NSIDE.EQ.0) GO TO 140
      VLB(I)=VLB(I)*SI
      VUB(I)=VUB(I)*SI
140   CONTINUE
150   CONTINUE
C     ------------------------------------------------------------------
C     ***** CALCULATE INITIAL FUNCTION AND CONSTRAINT VALUES  *****
C     ------------------------------------------------------------------
      INFO=1
      NCAL(1)=1
      IGOTO=1
      GO TO 950
160   CONTINUE
      OBJ1=OBJ
      IF (DABFUN.LE.0.) DABFUN=.001*ABS(OBJ)
      IF (DABFUN.LT.1.0E-10) DABFUN=1.0E-10
      IF (IPRINT.LE.0) GO TO 270
C     ------------------------------------------------------------------
C                    PRINT INITIAL DESIGN INFORMATION
C     ------------------------------------------------------------------
      IF (IPRINT.LE.1) GO TO 230
      IF (NSIDE.EQ.0.AND.NCON.EQ.0) WRITE (IOUT,1290)
      IF (NSIDE.NE.0.OR.NCON.GT.0) WRITE (IOUT,1230)
      WRITE (IOUT,1240) IPRINT,NDV,ITMAX,NCON,NSIDE,ICNDIR,NSCAL,NFDG
     1,LINOBJ,ITRM,N1,N2,N3,N4,N5
      WRITE (IOUT,1260) CT,CTMIN,CTL,CTLMIN,THETA,PHI,DELFUN,DABFUN
      WRITE (IOUT,1250) FDCH,FDCHM,ALPHAX,ABOBJ1
      IF (NSIDE.EQ.0) GO TO 190
      WRITE (IOUT,1270)
      DO 170 I=1,NDV,6
      M1=MIN0(NDV,I+5)
170   WRITE (IOUT,1010) I,(VLB(J),J=I,M1)
      WRITE (IOUT,1280)
      DO 180 I=1,NDV,6
      M1=MIN0(NDV,I+5)
180   WRITE (IOUT,1010) I,(VUB(J),J=I,M1)
190   CONTINUE
      IF (NSCAL.GE.0) GO TO 200
      WRITE (IOUT,1300)
      WRITE (IOUT,1460) (SCAL(I),I=1,NDV)
200   CONTINUE
      IF (NCON.EQ.0) GO TO 230
      IF (NLNC.EQ.0.OR.NLNC.EQ.NCON) GO TO 220
      WRITE (IOUT,1020)
      DO 210 I=1,NCON,15
      M1=MIN0(NCON,I+14)
210   WRITE (IOUT,1030) I,(ISC(J),J=I,M1)
      GO TO 230
220   IF (NLNC.EQ.NCON) WRITE (IOUT,1040)
      IF (NLNC.EQ.0) WRITE (IOUT,1050)
230   CONTINUE
      WRITE (IOUT,1440) OBJ
      WRITE (IOUT,1450)
      DO 240 I=1,NDV
      X1=1.
      IF (NSCAL.NE.0) X1=SCAL(I)
240   G1(I)=X(I)*X1
      DO 250 I=1,NDV,6
      M1=MIN0(NDV,I+5)
250   WRITE (IOUT,1010) I,(G1(J),J=I,M1)
      IF (NCON.EQ.0) GO TO 270
      WRITE (IOUT,1470)
      DO 260 I=1,NCON,6
      M1=MIN0(NCON,I+5)
260   WRITE (IOUT,1010) I,(G(J),J=I,M1)
270   CONTINUE
      IF (IPRINT.GT.1) WRITE (IOUT,1360)
C     ------------------------------------------------------------------
C     ********************  BEGIN MINIMIZATION  ************************
C     ------------------------------------------------------------------
280   CONTINUE
      ITER=ITER+1
      IF (ABOBJ1.LT..0001) ABOBJ1=.0001
      IF (ABOBJ1.GT..2) ABOBJ1=.2
      IF (ALPHAX.GT.1.) ALPHAX=1.
      IF (ALPHAX.LT..001) ALPHAX=.001
C
C  THE FOLLOWING TWO LINES OF CODE WERE COMMENTED OUT ON 3/5/81
C
C     NFEAS=NFEAS+1
C     IF (NFEAS.GT.10) GO TO 810
      IF (IPRINT.GT.2) WRITE (IOUT,1310) ITER
      IF (IPRINT.GT.3.AND.NCON.GT.0) WRITE (IOUT,1320) CT,CTL,PHI
      CTA=ABS(CT)
      IF (NCOBJ.EQ.0) GO TO 340
C     ------------------------------------------------------------------
C     NO MOVE ON LAST ITERATION.  DELETE CONSTRAINTS THAT ARE NO
C     LONGER ACTIVE.
C     ------------------------------------------------------------------
      NNAC=NAC
      DO 290 I=1,NNAC
      IF (IC(I).GT.NCON) NAC=NAC-1
290   CONTINUE
      IF (NAC.LE.0) GO TO 420
      NNAC=NAC
      DO 330 I=1,NNAC
300   NIC=IC(I)
      CT1=CT
      IF (ISC(NIC).GT.0) CT1=CTL
      IF (G(NIC).GT.CT1) GO TO 330
      NAC=NAC-1
      IF (I.GT.NAC) GO TO 420
      DO 320 K=I,NAC
      II=K+1
      DO 310 J=1,NDV2
310   A(J,K)=A(J,II)
320   IC(K)=IC(II)
      GO TO 300
330   CONTINUE
      GO TO 420
340   CONTINUE
      IF (MSCAL.LT.NSCAL.OR.NSCAL.EQ.0) GO TO 360
      IF (NSCAL.LT.0.AND.KCOUNT.LT.ICNDIR) GO TO 360
      MSCAL=0
      KCOUNT=0
C     ------------------------------------------------------------------
C                          SCALE VARIABLES
C     ------------------------------------------------------------------
      DO 350 I=1,NDV
      SI=SCAL(I)
      XI=SI*X(I)
      SIB=SI
      IF (NSCAL.GT.0) SI=ABS(XI)
      IF (SI.LT.1.0E-10) GO TO 350
      SCAL(I)=SI
      SI=1./SI
      X(I)=XI*SI
      IF (NSIDE.EQ.0) GO TO 350
      VLB(I)=SIB*SI*VLB(I)
      VUB(I)=SIB*SI*VUB(I)
350   CONTINUE
      IF (IPRINT.LT.4.OR.(NSCAL.LT.0.AND.ITER.GT.1)) GO TO 360
      WRITE (IOUT,1330)
      WRITE (IOUT,1460) (SCAL(I),I=1,NDV)
360   CONTINUE
      MSCAL=MSCAL+1
      NAC=0
C     ------------------------------------------------------------------
C          OBTAIN GRADIENTS OF OBJECTIVE AND ACTIVE CONSTRAINTS
C     ------------------------------------------------------------------
      INFO=2
      NCAL(2)=NCAL(2)+1
      IF (NFDG.NE.1) GO TO 370
      IGOTO=2
      GO TO 950
370   CONTINUE
      JGOTO=0
380   CONTINUE
      CALL CNMN01 (JGOTO,X,DF,G,ISC,IC,A,G1,VLB,VUB,SCAL,C,NCAL,DX,DX1
     1,FI,XI,III,N1,N2,N3,N4)
      IGOTO=3
      IF (JGOTO.GT.0) GO TO 950
390   CONTINUE
      INFO=1
      IF (NAC.GE.N3) GO TO 810
      IF (NSCAL.EQ.0.OR.NFDG.EQ.0) GO TO 420
C     ------------------------------------------------------------------
C                              SCALE GRADIENTS
C     ------------------------------------------------------------------
C     SCALE GRADIENT OF OBJECTIVE FUNCTION.
      DO 400 I=1,NDV
400   DF(I)=DF(I)*SCAL(I)
      IF (NFDG.EQ.2.OR.NAC.EQ.0) GO TO 420
C     SCALE GRADIENTS OF ACTIVE CONSTRAINTS.
      DO 410 J=1,NDV
      SCJ=SCAL(J)
      DO 410 I=1,NAC
410   A(J,I)=A(J,I)*SCJ
420   CONTINUE
      IF (IPRINT.LT.3.OR.NCON.EQ.0) GO TO 470
C     ------------------------------------------------------------------
C                                   PRINT
C     ------------------------------------------------------------------
C     PRINT ACTIVE AND VIOLATED CONSTRAINT NUMBERS.
      M1=0
      M2=N3
      IF (NAC.EQ.0) GO TO 450
      DO 440 I=1,NAC
      J=IC(I)
      IF (J.GT.NCON) GO TO 440
      GI=G(J)
      C1=CTAM
      IF (ISC(J).GT.0) C1=CTBM
      GI=GI-C1
      IF (GI.GT.0.) GO TO 430
C     ACTIVE CONSTRAINT.
      M1=M1+1
      MS1(M1)=J
      GO TO 440
430   M2=M2+1
C     VIOLATED CONSTRAINT.
      MS1(M2)=J
440   CONTINUE
450   M3=M2-N3
      WRITE (IOUT,1060) M1
      IF (M1.EQ.0) GO TO 460
      WRITE (IOUT,1070)
      WRITE (IOUT,1480) (MS1(I),I=1,M1)
460   WRITE (IOUT,1080) M3
      IF (M3.EQ.0) GO TO 470
      WRITE (IOUT,1070)
      M3=N3+1
      WRITE (IOUT,1480) (MS1(I),I=M3,M2)
470   CONTINUE
C     ------------------------------------------------------------------
C            CALCULATE GRADIENTS OF ACTIVE SIDE CONSTRAINTS
C     ------------------------------------------------------------------
      IF (NSIDE.EQ.0) GO TO 530
      MCN1=NCON
      M1=0
      DO 510 I=1,NDV
C     LOWER BOUND.
      XI=X(I)
      XID=VLB(I)
      X12=ABS(XID)
      IF (X12.LT.1.) X12=1.
      GI=(XID-XI)/X12
      IF (GI.LT.-1.0E-6) GO TO 490
      M1=M1+1
      MS1(M1)=-I
      NAC=NAC+1
      IF (NAC.GE.N3) GO TO 810
      MCN1=MCN1+1
      DO 480 J=1,NDV
480   A(J,NAC)=0.
      A(I,NAC)=-1.
      IC(NAC)=MCN1
      G(MCN1)=GI
      ISC(MCN1)=1
C     UPPER BOUND.
490   XID=VUB(I)
      X12=ABS(XID)
      IF (X12.LT.1.) X12=1.
      GI=(XI-XID)/X12
      IF (GI.LT.-1.0E-6) GO TO 510
      M1=M1+1
      MS1(M1)=I
      NAC=NAC+1
      IF (NAC.GE.N3) GO TO 810
      MCN1=MCN1+1
      DO 500 J=1,NDV
500   A(J,NAC)=0.
      A(I,NAC)=1.
      IC(NAC)=MCN1
      G(MCN1)=GI
      ISC(MCN1)=1
510   CONTINUE
C     ------------------------------------------------------------------
C                                  PRINT
C     ------------------------------------------------------------------
C     PRINT ACTIVE SIDE CONSTRAINT NUMBERS.
      IF (IPRINT.LT.3) GO TO 530
      WRITE (IOUT,1090) M1
      IF (M1.EQ.0) GO TO 530
      WRITE (IOUT,1100)
      WRITE(6,1480) (MS1(J),J=1,M1)
530   CONTINUE
C     PRINT GRADIENTS OF ACTIVE AND VIOLATED CONSTRAINTS.
      IF (IPRINT.LT.4) GO TO 570
      WRITE (IOUT,1340)
      DO 540 I=1,NDV,6
      M1=MIN0(NDV,I+5)
540   WRITE (IOUT,1010) I,(DF(J),J=I,M1)
      IF (NAC.EQ.0) GO TO 570
      WRITE (IOUT,1350)
      DO 560 I=1,NAC
      M1=IC(I)
      M2=M1-NCON
      M3=0
      IF (M2.GT.0) M3=IABS(MS1(M2))
      IF (M2.LE.0) WRITE (IOUT,990) M1
      IF (M2.GT.0) WRITE (IOUT,1000) M3
      DO 550 K=1,NDV,6
      M1=MIN0(NDV,K+5)
550   WRITE (IOUT,1010) K,(A(J,I),J=K,M1)
560   WRITE (IOUT,1360)
570   CONTINUE
C     ------------------------------------------------------------------
C     ******************  DETERMINE SEARCH DIRECTION *******************
C     ------------------------------------------------------------------
      ALP=1.0E+20
      IF (NAC.GT.0) GO TO 580
C     ------------------------------------------------------------------
C                        UNCONSTRAINED FUNCTION
C     ------------------------------------------------------------------
C     FIND DIRECTION OF STEEPEST DESCENT OR CONJUGATE DIRECTION.
C
C  S. N. 575 ADDED ON 2/25/81
C
 575  NVC=0
      NFEAS=0
      KCOUNT=KCOUNT+1
C     IF KCOUNT.GT.ICNDIR  RESTART CONJUGATE DIRECTION ALGORITHM.
      IF (KCOUNT.GT.ICNDIR.OR.IOBJ.EQ.2) KCOUNT=1
      IF (KCOUNT.EQ.1) JDIR=0
C     IF JDIR = 0 FIND DIRECTION OF STEEPEST DESCENT.
      CALL CNMN02 (JDIR,SLOPE,DFTDF1,DF,S,N1)
      GO TO 630
580   CONTINUE
C     ------------------------------------------------------------------
C                          CONSTRAINED FUNCTION
C     ------------------------------------------------------------------
C     FIND USABLE-FEASIBLE DIRECTION.
      KCOUNT=0
      JDIR=0
      PHI=10.*PHI
      IF (PHI.GT.1000.) PHI=1000.
C
C  THE FOLLOWING LINE OF CODE WAS COMMENTED OUT ON 3/5/81
C
C     IF (NFEAS.EQ.1) PHI=5.
C     CALCULATE DIRECTION, S.
      CALL CNMN05 (G,DF,A,S,B,C,SLOPE,PHI,ISC,IC,MS1,NVC,N1,N2,N3,N4,N5)
C
C  THE FOLLOWING LINE WAS ADDED ON 2/25/81
C
      IF(NAC.EQ.0) GO TO 575
C
C  THE FOLLOWING FIVE LINES WERE COMMENTED OUT ON 3/5/81
C  REASON : THEY WERE NOT IN G. VANDERPLAATS LISTING
C
C     IF THIS DESIGN IS FEASIBLE AND LAST ITERATION WAS INFEASIBLE,
C     SET ABOBJ1=.05 (5 PERCENT).
C     IF (NVC.EQ.0.AND.NFEAS.GT.1) ABOBJ1=.05
C     IF (NVC.EQ.0) NFEAS=0
      IF (IPRINT.LT.3) GO TO 600
      WRITE (IOUT,1370)
      DO 590 I=1,NAC,6
      M1=MIN0(NAC,I+5)
590   WRITE (IOUT,1010) I,(A(NDV1,J),J=I,M1)
      WRITE (IOUT,1210) S(NDV1)
600   CONTINUE
C     ------------------------------------------------------------------
C     ****************** ONE-DIMENSIONAL SEARCH ************************
C     ------------------------------------------------------------------
      IF (S(NDV1).LT.1.0E-6.AND.NVC.EQ.0) GO TO 710
C     ------------------------------------------------------------------
C                 FIND ALPHA TO OBTAIN A FEASIBLE DESIGN
C     ------------------------------------------------------------------
      IF (NVC.EQ.0) GO TO 630
      ALP=-1.
      DO 620 I=1,NAC
      NCI=IC(I)
      C1=G(NCI)
      CTC=CTAM
      IF (ISC(NCI).GT.0) CTC=CTBM
      IF (C1.LE.CTC) GO TO 620
      ALP1=0.
      DO 610 J=1,NDV
610   ALP1=ALP1+S(J)*A(J,I)
      ALP1=ALP1*A(NDV2,I)
      IF (ABS(ALP1).LT.1.0E-20) GO TO 620
      ALP1=-C1/ALP1
      IF (ALP1.GT.ALP) ALP=ALP1
620   CONTINUE
630   CONTINUE
C     ------------------------------------------------------------------
C                       LIMIT CHANCE TO ABOBJ1*OBJ
C     ------------------------------------------------------------------
      ALP1=1.0E+20
      SI=ABS(OBJ)
      IF (SI.LT..01) SI=.01
      IF (ABS(SLOPE).GT.1.0E-20) ALP1=ABOBJ1*SI/SLOPE
      ALP1=ABS(ALP1)
      IF (NVC.GT.0) ALP1=10.*ALP1
      IF (ALP1.LT.ALP) ALP=ALP1
C     ------------------------------------------------------------------
C                   LIMIT CHANGE IN VARIABLE TO ALPHAX
C     ------------------------------------------------------------------
      ALP11=1.0E+20
      DO 640 I=1,NDV
      SI=ABS(S(I))
      XI=ABS(X(I))
      IF (SI.LT.1.0E-10.OR.XI.LT.0.1) GO TO 640
      ALP1=ALPHAX*XI/SI
      IF (ALP1.LT.ALP11) ALP11=ALP1
640   CONTINUE
      IF (NVC.GT.0) ALP11=10.*ALP11
      IF (ALP11.LT.ALP) ALP=ALP11
      IF (ALP.GT.1.0E+20) ALP=1.0E+20
      IF (ALP.LE.1.0E-20) ALP=1.0E-20
      IF (IPRINT.LT.3) GO TO 660
      WRITE (IOUT,1380)
      DO 650 I=1,NDV,6
      M1=MIN0(NDV,I+5)
650   WRITE (IOUT,1010) I,(S(J),J=I,M1)
      WRITE (IOUT,1110) SLOPE,ALP
660   CONTINUE
      IF (NCON.GT.0.OR.NSIDE.GT.0) GO TO 680
C     ------------------------------------------------------------------
C           DO ONE-DIMENSIONAL SEARCH FOR UNCONSTRAINED FUNCTION
C     ------------------------------------------------------------------
      JGOTO=0
670   CONTINUE
      CALL CNMN03 (X,S,SLOPE,ALP,FFF,A1,A2,A3,A4,F1,F2,F3,F4,APP,N1
     1,NCAL,KOUNT,JGOTO)
      IGOTO=4
      IF (JGOTO.GT.0) GO TO 950
      JDIR=1
C     PROCEED TO CONVERGENCE CHECK.
      GO TO 700
C     ------------------------------------------------------------------
C       SOLVE ONE-DIMENSIONAL SEARCH PROBLEM FOR CONSTRAINED FUNCTION
C     ------------------------------------------------------------------
680   CONTINUE
      JGOTO=0
690   CONTINUE
      CALL CNMN06 (X,VLB,VUB,G,SCAL,DF,S,G1,G2,CTAM,CTBM,SLOPE,ALP,A2,A3
     1,A4,F1,F2,F3,CV1,CV2,CV3,CV4,ALPCA,ALPFES,ALPLN,ALPMIN,ALPNC,ALPSA
     2V,ALPSID,ALPTOT,ISC,N1,N2,NCAL,NVC,ICOUNT,IGOOD1,IGOOD2,IGOOD3,IGO
     3OD4,IBEST,III,NLNC,JGOTO)
      IGOTO=5
      IF (JGOTO.GT.0) GO TO 950
      IF (NAC.EQ.0) JDIR=1
C     ------------------------------------------------------------------
C     *******************     UPDATE ALPHAX   **************************
C     ------------------------------------------------------------------
700   CONTINUE
710   CONTINUE
      IF (ALP.GT.1.0E+19) ALP=0.
C     UPDATE ALPHAX TO BE AVERAGE OF MAXIMUM CHANGE IN X(I)
C     AND ALHPAX.
      ALP11=0.
      DO 720 I=1,NDV
      SI=ABS(S(I))
      XI=ABS(X(I))
      IF (XI.LT.1.0E-10) GO TO 720
      ALP1=ALP*SI/XI
      IF (ALP1.GT.ALP11) ALP11=ALP1
720   CONTINUE
      ALP11=.5*(ALP11+ALPHAX)
      ALP12=5.*ALPHAX
      IF (ALP11.GT.ALP12) ALP11=ALP12
      ALPHAX=ALP11
      NCOBJ=NCOBJ+1
C     ABSOLUTE CHANGE IN OBJECTIVE.
      OBJD=OBJ1-OBJ
      OBJB=ABS(OBJD)
      IF (OBJB.LT.1.0E-10) OBJB=0.
      IF (NAC.EQ.0.OR.OBJB.GT.0.) NCOBJ=0
      IF (NCOBJ.GT.1) NCOBJ=0
C     ------------------------------------------------------------------
C                                  PRINT
C     ------------------------------------------------------------------
C     PRINT MOVE PARAMETER, NEW X-VECTOR AND CONSTRAINTS.
      IF (IPRINT.LT.3) GO TO 730
      WRITE (IOUT,1390) ALP
730   IF (IPRINT.LT.2) GO TO 800
      IF (OBJB.GT.0.) GO TO 740
      IF (IPRINT.EQ.2) WRITE (IOUT,1400) ITER,OBJ
      IF (IPRINT.GT.2) WRITE (IOUT,1410) OBJ
      GO TO 760
740   IF (IPRINT.EQ.2) GO TO 750
      WRITE (IOUT,1420) OBJ
      GO TO 760
750   WRITE (IOUT,1430) ITER,OBJ
760   WRITE (IOUT,1450)
      DO 770 I=1,NDV
      FF1=1.
      IF (NSCAL.NE.0) FF1=SCAL(I)
770   G1(I)=FF1*X(I)
      DO 780 I=1,NDV,6
      M1=MIN0(NDV,I+5)
780   WRITE (IOUT,1010) I,(G1(J),J=I,M1)
      IF (NCON.EQ.0) GO TO 800
      WRITE (IOUT,1470)
      DO 790 I=1,NCON,6
      M1=MIN0(NCON,I+5)
790   WRITE (IOUT,1010) I,(G(J),J=I,M1)
800   CONTINUE
C
C  THE FOLLOWING CODE WAS ADDED ON 3/5/81
C
C  IT HAD NOT BEEN REPORTED AS A FIX TO MAOB
C  BUT WAS SENT TO JEFF STROUD A YEAR AGO
C  SEE OTHER COMMENTS IN CONMIN SUBROUTINE FOR DELETIONS OF CODE
C  ON 3/5/81 PERTAINING TO THIS FIX
C
C
C                   CHECK FEASIBILITY
C
      IF(NCON.LE.0) GO TO 808
      NFEASCT=10
      DO 804 I=1,NCON
      C1=CTAM
      IF(ISC(I).GT.0) C1=CTBM
      IF(G(I).LE.C1) GO TO 804
      NFEAS=NFEAS+1
      GO TO 806
 804  CONTINUE
      IF(NFEAS.GT.0) ABOBJ1=.05
      NFEAS=0
      PHI=5.
 806  IF(NFEAS.GE.NFEASCT) GO TO 810
 808  CONTINUE
C
C  END OF INSERTED FIX
C
C     ------------------------------------------------------------------
C                          CHECK CONVERGENCE
C     ------------------------------------------------------------------
C     STOP IF ITER EQUALS ITMAX.
      IF (ITER.GE.ITMAX) GO TO 810
C     ------------------------------------------------------------------
C                     ABSOLUTE CHANGE IN OBJECTIVE
C     ------------------------------------------------------------------
      OBJB=ABS(OBJD)
      KOBJ=KOBJ+1
      IF (OBJB.GE.DABFUN.OR.NFEAS.GT.0) KOBJ=0
C     ------------------------------------------------------------------
C                     RELATIVE CHANGE IN OBJECTIVE
C     ------------------------------------------------------------------
      IF (ABS(OBJ1).GT.1.0E-10) OBJD=OBJD/ABS(OBJ1)
      ABOBJ1=.5*(ABS(ABOBJ)+ABS(OBJD))
      ABOBJ=ABS(OBJD)
      IOBJ=IOBJ+1
      IF (NVC.GT.0.OR.OBJD.GE.DELFUN) IOBJ=0
      IF (IOBJ.GE.ITRM.OR.KOBJ.GE.ITRM) GO TO 810
      OBJ1=OBJ
C     ------------------------------------------------------------------
C           REDUCE CT IF OBJECTIVE FUNCTION IS CHANGING SLOWLY
C     ------------------------------------------------------------------
      IF (IOBJ.LT.1.OR.NAC.EQ.0) GO TO 280
      CT=DCT*CT
      CTL=CTL*DCTL
      IF (ABS(CT).LT.CTMIN) CT=-CTMIN
      IF (ABS(CTL).LT.CTLMIN) CTL=-CTLMIN
      GO TO 280
810   CONTINUE
      IF (NAC.GE.N3) WRITE (IOUT,1490)
C     ------------------------------------------------------------------
C     ****************  FINAL FUNCTION INFORMATION  ********************
C     ------------------------------------------------------------------
      IF (NSCAL.EQ.0) GO TO 830
C     UN-SCALE THE DESIGN VARIABLES.
      DO 820 I=1,NDV
      XI=SCAL(I)
      IF (NSIDE.EQ.0) GO TO 820
      VLB(I)=XI*VLB(I)
      VUB(I)=XI*VUB(I)
820   X(I)=XI*X(I)
C     ------------------------------------------------------------------
C                           PRINT FINAL RESULTS
C     ------------------------------------------------------------------
830   IF (IPRINT.EQ.0.OR.NAC.GE.N3) GO TO 940
      WRITE (IOUT,1500)
      WRITE (IOUT,1420) OBJ
      WRITE (IOUT,1450)
      DO 840 I=1,NDV,6
      M1=MIN0(NDV,I+5)
840   WRITE (IOUT,1010) I,(X(J),J=I,M1)
      IF (NCON.EQ.0) GO TO 900
      WRITE (IOUT,1470)
      DO 850 I=1,NCON,6
      M1=MIN0(NCON,I+5)
850   WRITE (IOUT,1010) I,(G(J),J=I,M1)
C     DETERMINE WHICH CONSTRAINTS ARE ACTIVE AND PRINT.
      NAC=0
      NVC=0
      DO 870 I=1,NCON
      CTA=CTAM
      IF (ISC(I).GT.0) CTA=CTBM
      GI=G(I)
      IF (GI.GT.CTA) GO TO 860
      IF (GI.LT.CT.AND.ISC(I).EQ.0) GO TO 870
      IF (GI.LT.CTL.AND.ISC(I).GT.0) GO TO 870
      NAC=NAC+1
      IC(NAC)=I
      GO TO 870
860   NVC=NVC+1
      MS1(NVC)=I
870   CONTINUE
      WRITE (IOUT,1060) NAC
      IF (NAC.EQ.0) GO TO 880
      WRITE (IOUT,1070)
      WRITE (IOUT,1480) (IC(J),J=1,NAC)
880   WRITE (IOUT,1080) NVC
      IF (NVC.EQ.0) GO TO 890
      WRITE (IOUT,1070)
      WRITE (IOUT,1480) (MS1(J),J=1,NVC)
890   CONTINUE
900   CONTINUE
      IF (NSIDE.EQ.0) GO TO 930
C     DETERMINE WHICH SIDE CONSTRAINTS ARE ACTIVE AND PRINT.
      NAC=0
      DO 920 I=1,NDV
      XI=X(I)
      XID=VLB(I)
      X12=ABS(XID)
      IF (X12.LT.1.) X12=1.
      GI=(XID-XI)/X12
      IF (GI.LT.-1.0E-6) GO TO 910
      NAC=NAC+1
      MS1(NAC)=-I
910   XID=VUB(I)
      X12=ABS(XID)
      IF (X12.LT.1.) X12=1.
      GI=(XI-XID)/X12
      IF (GI.LT.-1.0E-6) GO TO 920
      NAC=NAC+1
      MS1(NAC)=I
920   CONTINUE
      WRITE (IOUT,1090) NAC
      IF (NAC.EQ.0) GO TO 930
      WRITE (IOUT,1100)
      WRITE (IOUT,1480) (MS1(J),J=1,NAC)
930   CONTINUE
      WRITE (IOUT,1150)
      IF (ITER.GE.ITMAX) WRITE (IOUT,1160)
      IF (NFEAS.GE.NFEASCT) WRITE (IOUT,1170)
      IF (IOBJ.GE.ITRM) WRITE (IOUT,1180) ITRM
      IF (KOBJ.GE.ITRM) WRITE (IOUT,1190) ITRM
      WRITE (IOUT,1200) ITER
      WRITE (IOUT,1510) NCAL(1)
      IF (NCON.GT.0) WRITE (IOUT,1520) NCAL(1)
      IF (NFDG.NE.0) WRITE (IOUT,1530) NCAL(2)
      IF (NCON.GT.0.AND.NFDG.EQ.1) WRITE (IOUT,1540) NCAL(2)
C     ------------------------------------------------------------------
C                   RE-SET BASIC PARAMETERS TO INPUT VALUES
C     ------------------------------------------------------------------
940   ITRM=IDM1
      ITMAX=IDM2
      ICNDIR=IDM3
      DELFUN=DM1
      DABFUN=DM2
      CT=DM3
      CTMIN=DM4
      CTL=DM5
      CTLMIN=DM6
      THETA=DM7
      PHI=DM8
      FDCH=DM9
      FDCHM=DM10
      ABOBJ1=DM11
      ALPHAX=DM12
      IGOTO=0
950   CONTINUE
      IF (NSCAL.EQ.0.OR.IGOTO.EQ.0) RETURN
C     UN-SCALE VARIABLES.
      DO 960 I=1,NDV
      C(I)=X(I)
960   X(I)=X(I)*SCAL(I)
      RETURN
C     ------------------------------------------------------------------
C                                FORMATS
C     ------------------------------------------------------------------
C
C
970   FORMAT (///5X,72HA COMPLETELY UNCONSTRAINED FUNCTION WITH A LINEAR
     1 OBJECTIVE IS SPECIFIED//10X,8HLINOBJ =,I5/10X,8HNCON   =,I5/10X,8
     2HNSIDE  =,I5//5X,35HCONTROL RETURNED TO CALLING PROGRAM)
980   FORMAT (///5X,56HCONMIN HAS ACHIEVED A SOLUTION OF OBJ LESS THAN -
     11.0E+40/5X,32HSOLUTION APPEARS TO BE UNBOUNDED/5X,26HOPTIMIZATION
     2IS TERMINATED)
990   FORMAT (5X,17HCONSTRAINT NUMBER,I5)
1000  FORMAT (5X,27HSIDE CONSTRAINT ON VARIABLE,I5)
1010  FORMAT (3X,I5,1H),2X,6E13.5)
1020  FORMAT (/5X,35HLINEAR CONSTRAINT IDENTIFIERS (ISC)/5X,36HNON-ZERO
     1INDICATES LINEAR CONSTRAINT)
1030  FORMAT (3X,I5,1H),2X,15I5)
1040  FORMAT (/5X,26HALL CONSTRAINTS ARE LINEAR)
1050  FORMAT (/5X,30HALL CONSTRAINTS ARE NON-LINEAR)
1060  FORMAT (/5X,9HTHERE ARE,I5,19H ACTIVE CONSTRAINTS)
1070  FORMAT (5X,22HCONSTRAINT NUMBERS ARE)
1080  FORMAT (/5X,9HTHERE ARE,I5,21H VIOLATED CONSTRAINTS)
1090  FORMAT (/5X,9HTHERE ARE,I5,24H ACTIVE SIDE CONSTRAINTS)
1100  FORMAT (5X,43HDECISION VARIABLES AT LOWER OR UPPER BOUNDS,30H (MIN
     1US INDICATES LOWER BOUND))
1110  FORMAT (/5X,22HONE-DIMENSIONAL SEARCH/5X,15HINITIAL SLOPE =,E12.4,
     12X,16HPROPOSED ALPHA =,E12.4)
1120  FORMAT (///5X,35H* * CONMIN DETECTS VLB(I).GT.VUB(I)/5X,57HFIX IS
     1SET X(I)=VLB(I)=VUB(I) = .5*(VLB(I)+VUB(I) FOR I =,I5)
1130  FORMAT (///5X,41H* * CONMIN DETECTS INITIAL X(I).LT.VLB(I)/5X,6HX(
     1I) =,E12.4,2X,8HVLB(I) =,E12.4/5X,35HX(I) IS SET EQUAL TO VLB(I) F
     2OR I =,I5)
1140  FORMAT (///5X,41H* * CONMIN DETECTS INITIAL X(I).GT.VUB(I)/5X,6HX(
     1I) =,E12.4,2X,8HVUB(I) =,E12.4/5X,35HX(I) IS SET EQUAL TO VUB(I) F
     2OR I =,I5)
1150  FORMAT (/5X,21HTERMINATION CRITERION)
1160  FORMAT (10X,17HITER EQUALS ITMAX)
1170  FORMAT (10X,'NFEASCT CONSECUTIVE ITERATIONS FAILED TO PRODUCE A
     1FEASIBLE DESIGN')
1180  FORMAT (10X,43HABS(1-OBJ(I-1)/OBJ(I)) LESS THAN DELFUN FOR,I3,11H
     1ITERATIONS)
1190  FORMAT (10X,43HABS(OBJ(I)-OBJ(I-1))   LESS THAN DABFUN FOR,I3,11H
     1ITERATIONS)
1200  FORMAT (/5X,22HNUMBER OF ITERATIONS =,I5)
1210  FORMAT (/5X,28HCONSTRAINT PARAMETER, BETA =,E14.5)
1220  FORMAT (1H1,////12X,27(2H* )/12X,1H*,51X,1H*/12X,1H*,20X,11HC O N
     1M I N,20X,1H*/12X,1H*,51X,1H*/12X,1H*,15X,21H FORTRAN PROGRAM FOR
     2,15X,1H*/12X,1H*,51X,1H*/12X,1H*,9X,33HCONSTRAINED FUNCTION MINIMI
     3ZATION,9X,1H*/12X,1H*,51X,1H*/12X,27(2H* ))
1230  FORMAT (////5X,33HCONSTRAINED FUNCTION MINIMIZATION//5X,18HCONTROL
     1 PARAMETERS)
1240  FORMAT (/5X,60HIPRINT  NDV    ITMAX    NCON    NSIDE  ICNDIR   NSC
     1AL   NFDG/8I8//5X,12HLINOBJ  ITRM,5X,2HN1,6X,2HN2,6X,2HN3,6X,2HN4,
     26X,2HN5/8I8)
1250  FORMAT (/9X,4HFDCH,12X,5HFDCHM,11X,6HALPHAX,10X,6HABOBJ1/1X,4(2X,E
     114.5))
1260  FORMAT (/9X,2HCT,14X,5HCTMIN,11X,3HCTL,13X,6HCTLMIN/1X,4(2X,E14.5)
     1//9X,5HTHETA,11X,3HPHI,13X,6HDELFUN,10X,6HDABFUN/1X,4(2X,E14.5))
1270  FORMAT (/5X,40HLOWER BOUNDS ON DECISION VARIABLES (VLB))
1280  FORMAT (/5X,40HUPPER BOUNDS ON DECISION VARIABLES (VUB))
1290  FORMAT (////5X,35HUNCONSTRAINED FUNCTION MINIMIZATION//5X,18HCONTR
     1OL PARAMETERS)
1300  FORMAT (/5X,21HSCALING VECTOR (SCAL))
1310  FORMAT (////5X,22HBEGIN ITERATION NUMBER,I5)
1320  FORMAT (/5X,4HCT =,E14.5,5X,5HCTL =,E14.5,5X,5HPHI =,E14.5)
1330  FORMAT (/5X,25HNEW SCALING VECTOR (SCAL))
1340  FORMAT (/5X,15HGRADIENT OF OBJ)
1350  FORMAT (/5X,44HGRADIENTS OF ACTIVE AND VIOLATED CONSTRAINTS)
1360  FORMAT (1H )
1370  FORMAT (/5X,37HPUSH-OFF FACTORS, (THETA(I), I=1,NAC))
1380  FORMAT (/5X,27HSEARCH DIRECTION (S-VECTOR))
1390  FORMAT (/5X,18HCALCULATED ALPHA =,E14.5)
1400  FORMAT (////5X,6HITER =,I5,5X,5HOBJ =,E14.5,5X,16HNO CHANGE IN OBJ
     1)
1410  FORMAT (/5X,5HOBJ =,E15.6,5X,16HNO CHANGE ON OBJ)
1420  FORMAT (/5X,5HOBJ =,E15.6)
1430  FORMAT (////5X,6HITER =,I5,5X,5HOBJ =,E14.5)
1440  FORMAT (//5X,28HINITIAL FUNCTION INFORMATION//5X,5HOBJ =,E15.6)
1450  FORMAT (/5X,29HDECISION VARIABLES (X-VECTOR))
1460  FORMAT (3X,7E13.4)
1470  FORMAT (/5X,28HCONSTRAINT VALUES (G-VECTOR))
1480  FORMAT (5X,15I5)
1490  FORMAT (/5X,59HTHE NUMBER OF ACTIVE AND VIOLATED CONSTRAINTS EXCEE
     1DS N3-1./5X,66HDIMENSIONED SIZE OF MATRICES A AND B AND VECTOR IC
     2IS INSUFFICIENT/5X,61HOPTIMIZATION TERMINATED AND CONTROL RETURNED
     3 TO MAIN PROGRAM.)
1500  FORMAT (1H1,////4X,30HFINAL OPTIMIZATION INFORMATION)
1510  FORMAT (/5X,32HOBJECTIVE FUNCTION WAS EVALUATED,8X,I5,2X,5HTIMES)
1520  FORMAT (/5X,35HCONSTRAINT FUNCTIONS WERE EVALUATED,I10,2X,5HTIMES)
1530  FORMAT (/5X,36HGRADIENT OF OBJECTIVE WAS CALCULATED,I9,2X,5HTIMES)
1540  FORMAT (/5X,40HGRADIENTS OF CONSTRAINTS WERE CALCULATED,I5,2X,5HTI
     1MES)
      END
