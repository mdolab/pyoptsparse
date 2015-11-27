      SUBROUTINE CNMN01 (JGOTO,X,DF,G,ISC,IC,A,G1,VLB,VUB,SCAL,C,NCAL,DX
     1,DX1,FI,XI,III,N1,N2,N3,N4)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /CNMN1/ DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,CTL,CTLMIN,ALPHAX
     1,ABOBJ1,THETA,OBJ,NDV,NCON,NSIDE,IPRINT,NFDG,NSCAL,LINOBJ,ITMAX,IT
     2RM,ICNDIR,IGOTO,NAC,INFO,INFOG,ITER,NFEASCT 
      DIMENSION X(N1), DF(N1), G(N2), ISC(N2), IC(N3), A(N1,N3), G1(N2),
     1 VLB(N1), VUB(N1), SCAL(N1), NCAL(2), C(N4) 
C     ROUTINE TO CALCULATE GRADIENT INFORMATION BY FINITE DIFFERENCE. 
C     BY G. N. VANDERPLAATS                         JUNE, 1972. 
C     NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF. 
      IF (JGOTO.EQ.1) GO TO 10
      IF (JGOTO.EQ.2) GO TO 70
      INFOG=0 
      INF=INFO
      NAC=0 
      IF (LINOBJ.NE.0.AND.ITER.GT.1) GO TO 10 
C     ------------------------------------------------------------------
C                    GRADIENT OF LINEAR OBJECTIVE 
C     ------------------------------------------------------------------
      IF (NFDG.EQ.2) JGOTO=1
      IF (NFDG.EQ.2) RETURN 
10    CONTINUE
      JGOTO=0 
      IF (NFDG.EQ.2.AND.NCON.EQ.0) RETURN 
      IF (NCON.EQ.0) GO TO 40 
C     ------------------------------------------------------------------
C       * * * DETERMINE WHICH CONSTRAINTS ARE ACTIVE OR VIOLATED * * *
C     ------------------------------------------------------------------
      DO 20 I=1,NCON
      IF (G(I).LT.CT) GO TO 20
      IF (ISC(I).GT.0.AND.G(I).LT.CTL) GO TO 20 
      NAC=NAC+1 
      IF (NAC.GE.N3) RETURN 
      IC(NAC)=I 
20    CONTINUE
      IF (NFDG.EQ.2.AND.NAC.EQ.0) RETURN
      IF ((LINOBJ.GT.0.AND.ITER.GT.1).AND.NAC.EQ.0) RETURN
C     ------------------------------------------------------------------
C                  STORE VALUES OF CONSTRAINTS IN G1
C     ------------------------------------------------------------------
      DO 30 I=1,NCON
30    G1(I)=G(I)
40    CONTINUE
      JGOTO=0 
      IF (NAC.EQ.0.AND.NFDG.EQ.2) RETURN
C     ------------------------------------------------------------------
C                            CALCULATE GRADIENTS
C     ------------------------------------------------------------------
      INFOG=1 
      INFO=1
      FI=OBJ
      III=0 
50    III=III+1 
      XI=X(III) 
      DX=FDCH*XI
      DX=ABS(DX)
      FDCH1=FDCHM 
      IF (NSCAL.NE.0) FDCH1=FDCHM/SCAL(III) 
      IF (DX.LT.FDCH1) DX=FDCH1 
      X1=XI+DX
      IF (NSIDE.EQ.0) GO TO 60
      IF (X1.GT.VUB(III)) DX=-DX
60    DX1=1./DX 
      X(III)=XI+DX
      NCAL(1)=NCAL(1)+1 
C     ------------------------------------------------------------------
C                         FUNCTION EVALUATION 
C     ------------------------------------------------------------------
      JGOTO=2 
      RETURN
70    CONTINUE
      X(III)=XI 
      IF (NFDG.EQ.0) DF(III)=DX1*(OBJ-FI) 
      IF (NAC.EQ.0) GO TO 90
C     ------------------------------------------------------------------
C             DETERMINE GRADIENT COMPONENTS OF ACTIVE CONSTRAINTS 
C     ------------------------------------------------------------------
      DO 80 J=1,NAC 
      I1=IC(J)
80    A(III,J)=DX1*(G(I1)-G1(I1)) 
90    CONTINUE
      IF (III.LT.NDV) GO TO 50
      INFOG=0 
      INFO=INF
      JGOTO=0 
      OBJ=FI
      IF (NCON.EQ.0) RETURN 
C     ------------------------------------------------------------------
C             STORE CURRENT CONSTRAINT VALUES BACK IN G-VECTOR
C     ------------------------------------------------------------------
      DO 100 I=1,NCON 
100   G(I)=G1(I)
      RETURN
      END 