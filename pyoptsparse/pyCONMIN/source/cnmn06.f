      SUBROUTINE CNMN06 (X,VLB,VUB,G,SCAL,DF,S,G1,G2,CTAM,CTBM,SLOPE,ALP
     1,A2,A3,A4,F1,F2,F3,CV1,CV2,CV3,CV4,ALPCA,ALPFES,ALPLN,ALPMIN,ALPNC
     2,ALPSAV,ALPSID,ALPTOT,ISC,N1,N2,NCAL,NVC,ICOUNT,IGOOD1,IGOOD2,IGOO
     3D3,IGOOD4,IBEST,III,NLNC,JGOTO)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /CNMN1/ DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,CTL,CTLMIN,ALPHAX
     1,ABOBJ1,THETA,OBJ,NDV,NCON,NSIDE,IPRINT,NFDG,NSCAL,LINOBJ,ITMAX,IT
     2RM,ICNDIR,IGOTO,NAC,INFO,INFOG,ITER,NFEASCT
      COMMON /OUTPUT/ IOUT
      DIMENSION X(N1), VLB(N1), VUB(N1), G(N2), SCAL(N1), DF(N1), S(N1),
     1 G1(N2), G2(N2), ISC(N2), NCAL(2) 
C     ROUTINE TO SOLVE ONE-DIMENSIONAL SEARCH PROBLEM FOR CONSTRAINED 
C     FUNCTION MINIMIZATION.
C     BY G. N. VANDERPLAATS                           AUG., 1974. 
C     NASA-AMES RESEARCH CENTER, MOFFETT FIELD, CALIF.
C     OBJ = INITIAL AND FINAL FUNCTION VALUE. 
C     ALP = MOVE PARAMETER. 
C     SLOPE = INITIAL SLOPE.
C 
C     ALPSID = MOVE TO SIDE CONSTRAINT. 
C     ALPFES = MOVE TO FEASIBLE REGION. 
C     ALPNC = MOVE TO NEW NON-LINEAR CONSTRAINT.
C     ALPLN = MOVE TO LINEAR CONSTRAINT.
C     ALPCA = MOVE TO RE-ENCOUNTER CURRENTLY ACTIVE CONSTRAINT. 
C     ALPMIN = MOVE TO MINIMIZE FUNCTION. 
C     ALPTOT = TOTAL MOVE PARAMETER.
      ZRO=0.
      IF (JGOTO.EQ.0) GO TO 10
      GO TO (140,310,520),JGOTO 
10    IF (IPRINT.GE.5) WRITE (IOUT,730)
      ALPSAV=ALP
      ICOUNT=0
      ALPTOT=0. 
C     TOLERANCES. 
      CTAM=ABS(CTMIN) 
      CTBM=ABS(CTLMIN)
C     PROPOSED MOVE.
20    CONTINUE
C     ------------------------------------------------------------------
C     *****  BEGIN SEARCH OR IMPOSE SIDE CONSTRAINT MODIFICATION  ***** 
C     ------------------------------------------------------------------
      A2=ALPSAV 
      ICOUNT=ICOUNT+1 
      ALPSID=1.0E+20
C     INITIAL ALPHA AND OBJ.
      ALP=0.
      F1=OBJ
      KSID=0
      IF (NSIDE.EQ.0) GO TO 70
C     ------------------------------------------------------------------
C     FIND MOVE TO SIDE CONSTRAINT AND INSURE AGAINST VIOLATION OF
C     SIDE CONSTRAINTS
C     ------------------------------------------------------------------
      DO 60 I=1,NDV 
      SI=S(I) 
      IF (ABS(SI).GT.1.0E-20) GO TO 30
C     ITH COMPONENT OF S IS SMALL.  SET TO ZERO.
      S(I)=0. 
      SLOPE=SLOPE-SI*DF(I)
      GO TO 60
30    CONTINUE
      XI=X(I) 
      SI=1./SI
      IF (SI.GT.0.) GO TO 40
C     LOWER BOUND.
      XI2=VLB(I)
      XI1=ABS(XI2)
      IF (XI1.LT.1.) XI1=1. 
C     CONSTRAINT VALUE. 
      GI=(XI2-XI)/XI1 
      IF (GI.GT.-1.0E-6) GO TO 50 
C     PROPOSED MOVE TO LOWER BOUND. 
      ALPA=(XI2-XI)*SI
      IF (ALPA.LT.ALPSID) ALPSID=ALPA 
      GO TO 60
40    CONTINUE
C     UPPER BOUND.
      XI2=VUB(I)
      XI1=ABS(XI2)
      IF (XI1.LT.1.) XI1=1. 
C     CONSTRAINT VALUE. 
      GI=(XI-XI2)/XI1 
      IF (GI.GT.-1.0E-6) GO TO 50 
C     PROPOSED MOVE TO UPPER BOUND. 
      ALPA=(XI2-XI)*SI
      IF (ALPA.LT.ALPSID) ALPSID=ALPA 
      GO TO 60
50    CONTINUE
C     MOVE WILL VIOLATE SIDE CONSTRAINT.  SET S(I)=0. 
      SLOPE=SLOPE-S(I)*DF(I)
      S(I)=0. 
      KSID=KSID+1 
60    CONTINUE
C     ALPSID IS UPPER BOUND ON ALPHA. 
      IF (A2.GT.ALPSID) A2=ALPSID 
70    CONTINUE
C     ------------------------------------------------------------------
C               CHECK ILL-CONDITIONING
C     ------------------------------------------------------------------
      IF (KSID.EQ.NDV.OR.ICOUNT.GT.10) GO TO 710
      IF (NVC.EQ.0.AND.SLOPE.GT.0.) GO TO 710 
      ALPFES=-1.
      ALPMIN=-1.
      ALPLN=1.1*ALPSID
      ALPNC=ALPSID
      ALPCA=ALPSID
      IF (NCON.EQ.0) GO TO 90 
C     STORE CONSTRAINT VALUES IN G1.
      DO 80 I=1,NCON
      G1(I)=G(I)
80    CONTINUE
90    CONTINUE
C     ------------------------------------------------------------------
C                  MOVE A DISTANCE A2*S 
C     ------------------------------------------------------------------
      ALPTOT=ALPTOT+A2
      DO 100 I=1,NDV
      X(I)=X(I)+A2*S(I) 
100   CONTINUE
      IF (IPRINT.LT.5) GO TO 130
      WRITE (IOUT,740) A2
      IF (NSCAL.EQ.0) GO TO 120 
      DO 110 I=1,NDV
110   G(I)=SCAL(I)*X(I) 
      WRITE (IOUT,750) (G(I),I=1,NDV)
      GO TO 130 
120   WRITE (IOUT,750) (X(I),I=1,NDV)
C     ------------------------------------------------------------------
C                   UPDATE FUNCTION AND CONSTRAINT VALUES 
C     ------------------------------------------------------------------
130   NCAL(1)=NCAL(1)+1 
      JGOTO=1 
      RETURN
140   CONTINUE
      F2=OBJ
      IF (IPRINT.GE.5) WRITE (IOUT,760) F2 
      IF (IPRINT.LT.5.OR.NCON.EQ.0) GO TO 150 
      WRITE (IOUT,770) 
      WRITE (IOUT,750) (G(I),I=1,NCON) 
150   CONTINUE
C     ------------------------------------------------------------------
C               IDENTIFY ACCAPTABILITY OF DESIGNS F1 AND F2 
C     ------------------------------------------------------------------
C     IGOOD = 0 IS ACCAPTABLE.
C     CV = MAXIMUM CONSTRAINT VIOLATION.
      IGOOD1=0
      IGOOD2=0
      CV1=0.
      CV2=0.
      NVC1=0
      IF (NCON.EQ.0) GO TO 170
      DO 160 I=1,NCON 
      CC=CTAM 
      IF (ISC(I).GT.0) CC=CTBM
      C1=G1(I)-CC 
      C2=G(I)-CC
      IF (C2.GT.0.) NVC1=NVC1+1 
      IF (C1.GT.CV1) CV1=C1 
      IF (C2.GT.CV2) CV2=C2 
160   CONTINUE
      IF (CV1.GT.0.) IGOOD1=1 
      IF (CV2.GT.0.) IGOOD2=1 
170   CONTINUE
      ALP=A2
      OBJ=F2
C     ------------------------------------------------------------------
C     IF F2 VIOLATES FEWER CONSTRAINTS THAN F1 BUT STILL HAS CONSTRAINT 
C     VIOLATIONS RETURN 
C     ------------------------------------------------------------------
      IF (NVC1.LT.NVC.AND.NVC1.GT.0) GO TO 710
C     ------------------------------------------------------------------
C             IDENTIFY BEST OF DESIGNS F1 ANF F2
C     ------------------------------------------------------------------
C     IBEST CORRESPONDS TO MINIMUM VALUE DESIGN.
C     IF CONSTRAINTS ARE VIOLATED, IBEST CORRESPONDS TO MINIMUM 
C     CONSTRAINT VIOLATION. 
      IF (IGOOD1.EQ.0.AND.IGOOD2.EQ.0) GO TO 180
C     VIOLATED CONSTRAINTS.  PICK MINIMUM VIOLATION.
      IBEST=1 
      IF (CV1.GE.CV2) IBEST=2 
      GO TO 190 
180   CONTINUE
C     NO CONSTRAINT VIOLATION.  PICK MINIMUM F. 
      IBEST=1 
      IF (F2.LE.F1) IBEST=2 
190   CONTINUE
      II=1
C     ------------------------------------------------------------------
C     IF CV2 IS GREATER THAN CV1, SET MOVE LIMITS TO A2.
C     PROGRAM MOD-FEB, 1981, GV.
C     ------------------------------------------------------------------
      IF(CV2.LE.CV1) GO TO 195
      ALPLN=A2
      ALPNC=A2
      ALPCA=A2
195   CONTINUE
      IF (NCON.EQ.0) GO TO 230
C     ------------------------------------------------------------------
C     *****                 2 - POINT INTERPOLATION                *****
C     ------------------------------------------------------------------
      III=0 
200   III=III+1 
      C1=G1(III)
      C2=G(III) 
      IF (ISC(III).EQ.0) GO TO 210
C     ------------------------------------------------------------------
C                        LINEAR CONSTRAINT
C     ------------------------------------------------------------------
      IF (C1.GE.1.0E-5.AND.C1.LE.CTBM) GO TO 220
      CALL CNMN07 (II,ALP,ZRO,ZRO,C1,A2,C2,ZRO,ZRO) 
      IF (ALP.LE.0.) GO TO 220
      IF (C1.GT.CTBM.AND.ALP.GT.ALPFES) ALPFES=ALP
      IF (C1.LT.CTL.AND.ALP.LT.ALPLN) ALPLN=ALP 
      GO TO 220 
210   CONTINUE
C     ------------------------------------------------------------------
C                     NON-LINEAR CONSTRAINT 
C     ------------------------------------------------------------------
      IF (C1.GE.1.0E-5.AND.C1.LE.CTAM) GO TO 220
      CALL CNMN07 (II,ALP,ZRO,ZRO,C1,A2,C2,ZRO,ZRO) 
      IF (ALP.LE.0.) GO TO 220
      IF (C1.GT.CTAM.AND.ALP.GT.ALPFES) ALPFES=ALP
      IF (C1.LT.CT.AND.ALP.LT.ALPNC) ALPNC=ALP
220   CONTINUE
      IF (III.LT.NCON) GO TO 200
230   CONTINUE
      IF (LINOBJ.GT.0.OR.SLOPE.GE.0.) GO TO 240 
C     CALCULATE ALPHA TO MINIMIZE FUNCTION. 
      CALL CNMN04 (II,ALPMIN,ZRO,ZRO,F1,SLOPE,A2,F2,ZRO,ZRO,ZRO,ZRO)
240   CONTINUE
C     ------------------------------------------------------------------
C                         PROPOSED MOVE 
C     ------------------------------------------------------------------
C     MOVE AT LEAST FAR ENOUGH TO OVERCOME CONSTRAINT VIOLATIONS. 
      A3=ALPFES 
C     MOVE TO MINIMIZE FUNCTION.
      IF (ALPMIN.GT.A3) A3=ALPMIN 
C     IF A3.LE.0, SET A3 = ALPSID.
      IF (A3.LE.0.) A3=ALPSID 
C     LIMIT MOVE TO NEW CONSTRAINT ENCOUNTER. 
      IF (A3.GT.ALPNC) A3=ALPNC 
      IF (A3.GT.ALPLN) A3=ALPLN 
C     MAKE A3 NON-ZERO. 
      IF (A3.LE.1.0E-20) A3=1.0E-20 
C     IF A3=A2=ALPSID AND F2 IS BEST, GO INVOKE SIDE CONSTRAINT 
C     MODIFICATION. 
      ALPB=1.-A2/A3 
      ALPA=1.-ALPSID/A3 
      JBEST=0 
      IF (ABS(ALPB).LT.1.0E-10.AND.ABS(ALPA).LT.1.0E-10) JBEST=1
      IF (JBEST.EQ.1.AND.IBEST.EQ.2) GO TO 20 
C     SIDE CONSTRAINT CHECK NOT SATISFIED.
      IF (NCON.EQ.0) GO TO 260
C     STORE CONSTRAINT VALUES IN G2.
      DO 250 I=1,NCON 
      G2(I)=G(I)
250   CONTINUE
260   CONTINUE
C     IF A3=A2, SET A3=.9*A2. 
      IF (ABS(ALPB).LT.1.0E-10) A3=.9*A2
C     MOVE AT LEAST .01*A2. 
      IF (A3.LT.(.01*A2)) A3=.01*A2 
C     LIMIT MOVE TO 5.*A2.
      IF (A3.GT.(5.*A2)) A3=5.*A2 
C     LIMIT MOVE TO ALPSID. 
      IF (A3.GT.ALPSID) A3=ALPSID 
C     MOVE A DISTANCE A3*S. 
      ALP=A3-A2 
      ALPTOT=ALPTOT+ALP 
      DO 270 I=1,NDV
      X(I)=X(I)+ALP*S(I)
270   CONTINUE
      IF (IPRINT.LT.5) GO TO 300
      WRITE (IOUT,780) 
      WRITE (IOUT,740) A3
      IF (NSCAL.EQ.0) GO TO 290 
      DO 280 I=1,NDV
280   G(I)=SCAL(I)*X(I) 
      WRITE (IOUT,750) (G(I),I=1,NDV)
      GO TO 300 
290   WRITE (IOUT,750) (X(I),I=1,NDV)
300   CONTINUE
C     ------------------------------------------------------------------
C              UPDATE FUNCTION AND CONSTRAINT VALUES
C     ------------------------------------------------------------------
      NCAL(1)=NCAL(1)+1 
      JGOTO=2 
      RETURN
310   CONTINUE
      F3=OBJ
      IF (IPRINT.GE.5) WRITE (IOUT,760) F3 
      IF (IPRINT.LT.5.OR.NCON.EQ.0) GO TO 320 
      WRITE (IOUT,770) 
      WRITE (IOUT,750) (G(I),I=1,NCON) 
320   CONTINUE
C     ------------------------------------------------------------------
C       CALCULATE MAXIMUM CONSTRAINT VIOLATION AND PICK BEST DESIGN 
C     ------------------------------------------------------------------
      CV3=0.
      IGOOD3=0
      NVC1=0
      IF (NCON.EQ.0) GO TO 340
      DO 330 I=1,NCON 
      CC=CTAM 
      IF (ISC(I).GT.0) CC=CTBM
      C1=G(I)-CC
      IF (C1.GT.CV3) CV3=C1 
      IF (C1.GT.0.) NVC1=NVC1+1 
330   CONTINUE
      IF (CV3.GT.0.) IGOOD3=1 
340   CONTINUE
C     DETERMINE BEST DESIGN.
      IF (IBEST.EQ.2) GO TO 360 
C     CHOOSE BETWEEN F1 AND F3. 
      IF (IGOOD1.EQ.0.AND.IGOOD3.EQ.0) GO TO 350
      IF (CV1.GE.CV3) IBEST=3 
      GO TO 380 
350   IF (F3.LE.F1) IBEST=3 
      GO TO 380 
360   CONTINUE
C     CHOOSE BETWEEN F2 AND F3. 
      IF (IGOOD2.EQ.0.AND.IGOOD3.EQ.0) GO TO 370
      IF (CV2.GE.CV3) IBEST=3 
      GO TO 380 
370   IF (F3.LE.F2) IBEST=3 
380   CONTINUE
      ALP=A3
      OBJ=F3
C     IF F3 VIOLATES FEWER CONSTRAINTS THAN F1 RETURN.
      IF (NVC1.LT.NVC) GO TO 710
C     IF OBJECTIVE AND ALL CONSTRAINTS ARE LINEAR, RETURN.
      IF (LINOBJ.NE.0.AND.NLNC.EQ.NCON) GO TO 710 
C     IF A3 = ALPLN AND F3 IS BOTH GOOD AND BEST RETURN.
      ALPB=1.-ALPLN/A3
      IF ((ABS(ALPB).LT.1.0E-20.AND.IBEST.EQ.3).AND.(IGOOD3.EQ.0)) GO TO
     1 710
C     IF A3 = ALPSID AND F3 IS BEST, GO INVOKE SIDE CONSTRAINT
C     MODIFICATION. 
      ALPA=1.-ALPSID/A3 
      IF (ABS(ALPA).LT.1.0E-20.AND.IBEST.EQ.3) GO TO 20 
C     ------------------------------------------------------------------
C     **********            3 - POINT INTERPOLATION            *********
C     ------------------------------------------------------------------
      ALPNC=ALPSID
      ALPCA=ALPSID
      ALPFES=-1.
      ALPMIN=-1.
C     ------------------------------------------------------------------
C     IF A3 IS GREATER THAN A2 AND CV3 IS GREATER THAN CV2, SET 
C     MOVE LIMITS TO A3.  PROGRAM MOD-FEB, 1981, GV.
C     ------------------------------------------------------------------
      IF(A3.LE.A2.OR.CV3.LE.CV2) GO TO 285
      ALPLN=A3
      ALPNC=A3
      ALPCA=A3
285   CONTINUE
      IF (NCON.EQ.0) GO TO 440
      III=0 
390   III=III+1 
      C1=G1(III)
      C2=G2(III)
      C3=G(III) 
      IF (ISC(III).EQ.0) GO TO 400
C     ------------------------------------------------------------------
C     LINEAR CONSTRAINT.  FIND ALPFES ONLY.  ALPLN SAME AS BEFORE.
C     ------------------------------------------------------------------
      IF (C1.LE.CTBM) GO TO 430 
      II=1
      CALL CNMN07 (II,ALP,ZRO,ZRO,C1,A3,C3,ZRO,ZRO) 
      IF (ALP.GT.ALPFES) ALPFES=ALP 
      GO TO 430 
400   CONTINUE
C     ------------------------------------------------------------------
C                     NON-LINEAR CONSTRAINT 
C     ------------------------------------------------------------------
      II=2
      CALL CNMN07 (II,ALP,ZRO,ZRO,C1,A2,C2,A3,C3) 
      IF (ALP.LE.ZRO) GO TO 430 
      IF (C1.GE.CT.AND.C1.LE.0.) GO TO 410
      IF (C1.GT.CTAM.OR.C1.LT.0.) GO TO 420 
C     ALP IS MINIMUM MOVE.  UPDATE FOR NEXT CONSTRAINT ENCOUNTER. 
410   ALPA=ALP
      CALL CNMN07 (II,ALP,ALPA,ZRO,C1,A2,C2,A3,C3)
      IF (ALP.LT.ALPCA.AND.ALP.GE.ALPA) ALPCA=ALP 
      GO TO 430 
420   CONTINUE
      IF (ALP.GT.ALPFES.AND.C1.GT.CTAM) ALPFES=ALP
      IF (ALP.LT.ALPNC.AND.C1.LT.0.) ALPNC=ALP
430   CONTINUE
      IF (III.LT.NCON) GO TO 390
440   CONTINUE
      IF (LINOBJ.GT.0.OR.SLOPE.GT.0.) GO TO 450 
C     ------------------------------------------------------------------
C              CALCULATE ALPHA TO MINIMIZE FUNCTION 
C     ------------------------------------------------------------------
      II=3
      IF (A2.GT.A3.AND.(IGOOD2.EQ.0.AND.IBEST.EQ.2)) II=2 
      CALL CNMN04 (II,ALPMIN,ZRO,ZRO,F1,SLOPE,A2,F2,A3,F3,ZRO,ZRO)
450   CONTINUE
C     ------------------------------------------------------------------
C                       PROPOSED MOVE 
C     ------------------------------------------------------------------
C     MOVE AT LEAST ENOUGH TO OVERCOME CONSTRAINT VIOLATIONS. 
      A4=ALPFES 
C     MOVE TO MINIMIZE FUNCTION.
      IF (ALPMIN.GT.A4) A4=ALPMIN 
C     IF A4.LE.0, SET A4 = ALPSID.
      IF (A4.LE.0.) A4=ALPSID 
C     LIMIT MOVE TO NEW CONSTRAINT ENCOUNTER. 
      IF (A4.GT.ALPLN) A4=ALPLN 
      IF (A4.GT.ALPNC) A4=ALPNC 
C     LIMIT MOVE TO RE-ENCOUNTER CURRENTLY ACTIVE CONSTRAINT. 
      IF (A4.GT.ALPCA) A4=ALPCA 
C     LIMIT A4 TO 5.*A3.
      IF (A4.GT.(5.*A3)) A4=5.*A3 
C     UPDATE DESIGN.
      IF (IBEST.NE.3.OR.NCON.EQ.0) GO TO 470
C     STORE CONSTRAINT VALUES IN G2.  F3 IS BEST.  F2 IS NOT. 
      DO 460 I=1,NCON 
      G2(I)=G(I)
460   CONTINUE
470   CONTINUE
C     IF A4=A3 AND IGOOD1=0 AND IGOOD3=1, SET A4=.9*A3. 
      ALP=A4-A3 
      IF ((IGOOD1.EQ.0.AND.IGOOD3.EQ.1).AND.ABS(ALP).LT.1.0E-20) A4=.9*A
     13 
C     ------------------------------------------------------------------
C                   MOVE A DISTANCE A4*S
C     ------------------------------------------------------------------
      ALP=A4-A3 
      ALPTOT=ALPTOT+ALP 
      DO 480 I=1,NDV
      X(I)=X(I)+ALP*S(I)
480   CONTINUE
      IF (IPRINT.LT.5) GO TO 510
      WRITE (IOUT,720) 
      WRITE (IOUT,740) A4
      IF (NSCAL.EQ.0) GO TO 500 
      DO 490 I=1,NDV
490   G(I)=SCAL(I)*X(I) 
      WRITE (IOUT,750) (G(I),I=1,NDV)
      GO TO 510 
500   WRITE (IOUT,750) (X(I),I=1,NDV)
510   CONTINUE
C     ------------------------------------------------------------------
C              UPDATE FUNCTION AND CONSTRAINT VALUES
C     ------------------------------------------------------------------
      NCAL(1)=NCAL(1)+1 
      JGOTO=3 
      RETURN
520   CONTINUE
      F4=OBJ
      IF (IPRINT.GE.5) WRITE (IOUT,760) F4 
      IF (IPRINT.LT.5.OR.NCON.EQ.0) GO TO 530 
      WRITE (IOUT,770) 
      WRITE (IOUT,750) (G(I),I=1,NCON) 
530   CONTINUE
C     DETERMINE ACCAPTABILITY OF F4.
      IGOOD4=0
      CV4=0.
      IF (NCON.EQ.0) GO TO 550
      DO 540 I=1,NCON 
      CC=CTAM 
      IF (ISC(I).GT.0) CC=CTBM
      C1=G(I)-CC
      IF (C1.GT.CV4) CV4=C1 
540   CONTINUE
      IF (CV4.GT.0.) IGOOD4=1 
550   CONTINUE
      ALP=A4
      OBJ=F4
C     ------------------------------------------------------------------
C                     DETERMINE BEST DESIGN 
C     ------------------------------------------------------------------
      GO TO (560,610,660),IBEST 
560   CONTINUE
C     CHOOSE BETWEEN F1 AND F4. 
      IF (IGOOD1.EQ.0.AND.IGOOD4.EQ.0) GO TO 570
      IF (CV1.GT.CV4) GO TO 710 
      GO TO 580 
570   CONTINUE
      IF (F4.LE.F1) GO TO 710 
580   CONTINUE
C     F1 IS BEST. 
      ALPTOT=ALPTOT-A4
      OBJ=F1
      DO 590 I=1,NDV
      X(I)=X(I)-A4*S(I) 
590   CONTINUE
      IF (NCON.EQ.0) GO TO 710
      DO 600 I=1,NCON 
      G(I)=G1(I)
600   CONTINUE
      GO TO 710 
610   CONTINUE
C     CHOOSE BETWEEN F2 AND F4. 
      IF (IGOOD2.EQ.0.AND.IGOOD4.EQ.0) GO TO 620
      IF (CV2.GT.CV4) GO TO 710 
      GO TO 630 
620   CONTINUE
      IF (F4.LE.F2) GO TO 710 
630   CONTINUE
C     F2 IS BEST. 
      OBJ=F2
      A2=A4-A2
      ALPTOT=ALPTOT-A2
      DO 640 I=1,NDV
      X(I)=X(I)-A2*S(I) 
640   CONTINUE
      IF (NCON.EQ.0) GO TO 710
      DO 650 I=1,NCON 
      G(I)=G2(I)
650   CONTINUE
      GO TO 710 
660   CONTINUE
C     CHOOSE BETWEEN F3 AND F4. 
      IF (IGOOD3.EQ.0.AND.IGOOD4.EQ.0) GO TO 670
      IF (CV3.GT.CV4) GO TO 710 
      GO TO 680 
670   CONTINUE
      IF (F4.LE.F3) GO TO 710 
680   CONTINUE
C     F3 IS BEST. 
      OBJ=F3
      A3=A4-A3
      ALPTOT=ALPTOT-A3
      DO 690 I=1,NDV
      X(I)=X(I)-A3*S(I) 
690   CONTINUE
      IF (NCON.EQ.0) GO TO 710
      DO 700 I=1,NCON 
      G(I)=G2(I)
700   CONTINUE
710   CONTINUE
      ALP=ALPTOT
      IF (IPRINT.GE.5) WRITE (IOUT,790)
      JGOTO=0 
      RETURN
C     ------------------------------------------------------------------
C                                  FORMATS
C     ------------------------------------------------------------------
C 
C 
720   FORMAT (/5X,25HTHREE-POINT INTERPOLATION) 
730   FORMAT (/////58H* * * CONSTRAINED ONE-DIMENSIONAL SEARCH INFORMATI
     1ON * * *) 
740   FORMAT (//5X,15HPROPOSED DESIGN/5X,7HALPHA =,E12.5/5X,8HX-VECTOR) 
750   FORMAT (1X,8E12.4)
760   FORMAT (/5X,5HOBJ =,E13.5)
770   FORMAT (/5X,17HCONSTRAINT VALUES) 
780   FORMAT (/5X,23HTWO-POINT INTERPOLATION) 
790   FORMAT (/5X,35H* * * END OF ONE-DIMENSIONAL SEARCH) 
      END 
