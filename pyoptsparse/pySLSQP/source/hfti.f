      SUBROUTINE HFTI(A,MDA,M,N,B,MDB,NB,TAU,KRANK,RNORM,H,G,IP)

C     RANK-DEFICIENT LEAST SQUARES ALGORITHM AS DESCRIBED IN:
C     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY, 1973 JUN 12
C     TO APPEAR IN 'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974

C     A(*,*),MDA,M,N   THE ARRAY A INITIALLY CONTAINS THE M x N MATRIX A
C                      OF THE LEAST SQUARES PROBLEM AX = B.
C                      THE FIRST DIMENSIONING PARAMETER MDA MUST SATISFY
C                      MDA >= M. EITHER M >= N OR M < N IS PERMITTED.
C                      THERE IS NO RESTRICTION ON THE RANK OF A.
C                      THE MATRIX A WILL BE MODIFIED BY THE SUBROUTINE.
C     B(*,*),MDB,NB    IF NB = 0 THE SUBROUTINE WILL MAKE NO REFERENCE
C                      TO THE ARRAY B. IF NB > 0 THE ARRAY B() MUST
C                      INITIALLY CONTAIN THE M x NB MATRIX B  OF THE
C                      THE LEAST SQUARES PROBLEM AX = B AND ON RETURN
C                      THE ARRAY B() WILL CONTAIN THE N x NB SOLUTION X.
C                      IF NB>1 THE ARRAY B() MUST BE DOUBLE SUBSCRIPTED
C                      WITH FIRST DIMENSIONING PARAMETER MDB>=MAX(M,N),
C                      IF NB=1 THE ARRAY B() MAY BE EITHER SINGLE OR
C                      DOUBLE SUBSCRIPTED.
C     TAU              ABSOLUTE TOLERANCE PARAMETER FOR PSEUDORANK
C                      DETERMINATION, PROVIDED BY THE USER.
C     KRANK            PSEUDORANK OF A, SET BY THE SUBROUTINE.
C     RNORM            ON EXIT, RNORM(J) WILL CONTAIN THE EUCLIDIAN
C                      NORM OF THE RESIDUAL VECTOR FOR THE PROBLEM
C                      DEFINED BY THE J-TH COLUMN VECTOR OF THE ARRAY B.
C     H(), G()         ARRAYS OF WORKING SPACE OF LENGTH >= N.
C     IP()             INTEGER ARRAY OF WORKING SPACE OF LENGTH >= N
C                      RECORDING PERMUTATION INDICES OF COLUMN VECTORS

      INTEGER          I,J,JB,K,KP1,KRANK,L,LDIAG,LMAX,M,
     .                 MDA,MDB,N,NB,IP(N)
      DOUBLE PRECISION A(MDA,N),B(MDB,NB),H(N),G(N),RNORM(NB),FACTOR,
     .                 TAU,ZERO,HMAX,DIFF,TMP,DDOT,DNRM2,U,V
      DIFF(U,V)=       U-V
      DATA             ZERO/0.0D0/, FACTOR/1.0D-3/

      K=0
      LDIAG=MIN(M,N)
      IF(LDIAG.LE.0)                  GOTO 270

C   COMPUTE LMAX

      DO 80 J=1,LDIAG
          IF(J.EQ.1)                  GOTO 20
          LMAX=J
          DO 10 L=J,N
              H(L)=H(L)-A(J-1,L)**2
   10         IF(H(L).GT.H(LMAX)) LMAX=L
          IF(DIFF(HMAX+FACTOR*H(LMAX),HMAX).GT.ZERO)
     .                                GOTO 50
   20     LMAX=J
          DO 40 L=J,N
              H(L)=ZERO
              DO 30 I=J,M
   30             H(L)=H(L)+A(I,L)**2
   40         IF(H(L).GT.H(LMAX)) LMAX=L
          HMAX=H(LMAX)

C   COLUMN INTERCHANGES IF NEEDED

   50     IP(J)=LMAX
          IF(IP(J).EQ.J)              GOTO 70
          DO 60 I=1,M
              TMP=A(I,J)
              A(I,J)=A(I,LMAX)
   60         A(I,LMAX)=TMP
          H(LMAX)=H(J)

C   J-TH TRANSFORMATION AND APPLICATION TO A AND B

   70     I=MIN(J+1,N)
          CALL H12(1,J,J+1,M,A(1,J),1,H(J),A(1,I),1,MDA,N-J)
   80     CALL H12(2,J,J+1,M,A(1,J),1,H(J),B,1,MDB,NB)

C   DETERMINE PSEUDORANK

      DO 90 J=1,LDIAG
   90     IF(ABS(A(J,J)).LE.TAU)      GOTO 100
      K=LDIAG
      GOTO 110
  100 K=J-1
  110 KP1=K+1

C   NORM OF RESIDUALS

      DO 130 JB=1,NB
  130     RNORM(JB)=DNRM2(M-K,B(KP1,JB),1)
      IF(K.GT.0)                      GOTO 160
      DO 150 JB=1,NB
          DO 150 I=1,N
  150         B(I,JB)=ZERO
      GOTO 270
  160 IF(K.EQ.N)                      GOTO 180

C   HOUSEHOLDER DECOMPOSITION OF FIRST K ROWS

      DO 170 I=K,1,-1
  170     CALL H12(1,I,KP1,N,A(I,1),MDA,G(I),A,MDA,1,I-1)
  180 DO 250 JB=1,NB

C   SOLVE K*K TRIANGULAR SYSTEM

          DO 210 I=K,1,-1
              J=MIN(I+1,N)
  210         B(I,JB)=(B(I,JB)-DDOT(K-I,A(I,J),MDA,B(J,JB),1))/A(I,I)

C   COMPLETE SOLUTION VECTOR

          IF(K.EQ.N)                  GOTO 240
          DO 220 J=KP1,N
  220         B(J,JB)=ZERO
          DO 230 I=1,K
  230         CALL H12(2,I,KP1,N,A(I,1),MDA,G(I),B(1,JB),1,MDB,1)

C   REORDER SOLUTION ACCORDING TO PREVIOUS COLUMN INTERCHANGES

  240     DO 250 J=LDIAG,1,-1
              IF(IP(J).EQ.J)          GOTO 250
              L=IP(J)
              TMP=B(L,JB)
              B(L,JB)=B(J,JB)
              B(J,JB)=TMP
  250 CONTINUE
  270 KRANK=K
      END
      