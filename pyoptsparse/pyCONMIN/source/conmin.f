      SUBROUTINE CONMIN(NDV_,NCON_,X_,VLB_,VUB_,OBJ_,G_,
     1     N1,N2,N3,N4,N5,IPRINT_,IOUT_,IFILE,ITMAX_,DELFUN_,
     2     DABFUN_,ITRM_,NFEASCT_,NFDG_,NFUN_,NGRD_,CNMNFUN,CNMNGRD)
      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      EXTERNAL CNMNFUN,CNMNGRD
      
      CHARACTER*(*) IFILE
      
      DIMENSION X_(NDV_),VLB_(NDV_),VUB_(NDV_),G_(NCON_)
      
      DIMENSION X(N1),VLB(N1),VUB(N1),SCAL(N1),S(N1),DF(N1)
      DIMENSION G(N2),G1(N2),G2(N2),ISC(N2)
      DIMENSION IC(N3),B(N3,N3)
      DIMENSION C(N4)
      DIMENSION MS1(N5)
      DIMENSION A(N1,N3)
      
      COMMON /CNMN1/ DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,CTL,CTLMIN,
     .               ALPHAX,ABOBJ1,THETA,OBJ,NDV,NCON,NSIDE,IPRINT, 
     .               NFDG,NSCAL,LINOBJ,ITMAX,ITRM,ICNDIR,IGOTO,NAC, 
     .               INFO,INFOG,ITER,NFEASCT
      
      COMMON /VARABLE/ AOBJ
      COMMON /OUTPUT/ IOUT
      COMMON /FEVALS/ NFUN,NGRD
C 
C  INITIALIZE 
C 
      INFOG=0 
      INFO=0
      NDV=NDV_
      NCON=NCON_
      DO 5 I=1,NDV
        X(I)=X_(I)
        VLB(I)=VLB_(I)
        VUB(I)=VUB_(I)
    5 CONTINUE
      IPRINT=IPRINT_
      ITMAX=ITMAX_
      DO 6 J=1,NCON 
        ISC(J)=0
    6 CONTINUE
      NFDG=NFDG_
      NSIDE=0
      ICNDIR=0
      NSCAL=0
      NFEASCT=NFEASCT_
      FDCH=0.0
      FDCHM=0.0 
      CT=0.0
      CTMIN=0.0
      CTL=0.0 
      CTLMIN=0.0
      THETA=0.0 
      PHI=0.0 
      DELFUN=DELFUN_
      DABFUN=DABFUN_
      LINOBJ=0.0
      ITRM=ITRM_
      ALPHAX=0.0
      ABOBJ1=0.0
C      
      NFUN=0
      NGRD=0
C 
C     OPEN WRITE FILE
C     
      IOUT=IOUT_
      IF (IPRINT.EQ.0) GO TO 10
      OPEN(UNIT=IOUT,FILE=IFILE(1:LEN_TRIM(IFILE)),
     .   STATUS='UNKNOWN')
   10 CONTINUE
C
C     MAXIMUM NUMBER OF ITERATIONS
C 
      NLIM=ITMAX*(NDV+5)
C 
C     NON-ITERATIVE PART OF ANALYSIS
C 
      IGOTO = 0 
C 
C     ITERATIVE PART OF ANALYSIS
C 
      DO 20 I = 1,NLIM
C
        LOOPCNT=I 
C 
C       CALL THE OPTIMIZATION ROUTINE CONMIN
C 
        CALL CNMN00(X,VLB,VUB,G,SCAL,DF,A,S,G1,G2,B,
     .   C,ISC,IC,MS1,N1,N2,N3,N4,N5)
C
C     CHECK TERMINATION CRITERIA
C
        IF(IGOTO.EQ.0) LOOPCNT=-999 
C
C       ANALYSIS MODULE
C
        CALL CNMN09(CNMNFUN,CNMNGRD,X,G,IC,DF,A,
     .   N1,N2,N3,N4,N5)
        OBJ=AOBJ
C
        IF (IGOTO.EQ.0) GO TO 30
   20 CONTINUE
C
   30 CONTINUE
C
C  PRINT FINAL RESULTS 
C
      IF (IPRINT.EQ.0) GO TO 32
      WRITE(6,1650) NFUN-1
      WRITE(6,1750) NGRD
C
   32 CONTINUE
C
C     OUTPUT HANDLING
C
      DO 35 I=1,NDV
        X_(I)=X(I)
   35 CONTINUE
      OBJ_=OBJ
      DO 36 J=1,NCON
        G_(J)=G(J)
   36 CONTINUE
      NFUN_=NFUN-1
      NGRD_=NGRD
      
      RETURN

C  ------------------------------------------------------------------
C  FORMATS
C  ------------------------------------------------------------------
1650  FORMAT(//8X,'NUMBER OF FUNC-CALLS:  NFUN =',I5)
1750  FORMAT(8X,'NUMBER OF GRAD-CALLS:  NGRD =',I5)
      
      END
