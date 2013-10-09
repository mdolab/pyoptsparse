*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  sn87sopt.f
*
*     s8dflt   s8Map    s8solv   s8SQP
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8dflt
     &   ( m, n, nnCon, nnJac, nnObj, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     nnCon, nnJac, nnObj, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     s8dflt checks the optional parameter values and possibly changes
!     them to reasonable values.
!
!     Note that checking occurs before the amount of working storage has
!     been defined.
!
!     See  snworkspace.doc  for full documentation of cw, iw and rw.
!
!     15 Nov 1991: first version.
!     27 Apr 2001: wtMax introduced.
!     10 Dec 2002: Added defaults for LU Rook and Diagonal Pivoting.
!     31 Dec 2002: Added default MPS character names.
!     30 Jul 2003: Added default CG tolerance.
!     22 Jun 2004: Added default LU mod singularity tol
!     20 Dec 2004: Default LU tols reduced.
!     21 Dec 2004: Default LU tols fixed up.
!     09 Jun 2005: Default tolpiv same as in s5dflt.
!     02 Jul 2005: Default Utol's set back to eps1.
!     02 May 2006: lvlTim removed.
!     ==================================================================
      character
     &     blank*8, cdummy*8, mProb*8, mObj*8, mRhs*8, mRng*8, mBnd*8
      logical
     &     linCon, linear, nlnCon, nonlin
      integer
     &     cgItmx, DerOpt, iCrash, iBack, iDump, iLoadB, iInsrt, iNewB,
     &     iOldB, iPnch, iPrint, iReprt, iSoln, itnlim,
     &     jverf1, jverf2, jverf3, jverf4,
     &     kchk, kdegen, kFac, klog, kreset, ksav, kSumm, lEmode,
     &     lprDbg, lprPrm, lprSch, lprScl, lprSol,
     &     lvlDer, lvlHes, lvlInf, lvlPiv, lvlPre, lvlPPm,
     &     lvlSch, lvlScl, lvlSys, lvlVer, m, maxmn, maxCol,
     &     maxR, maxS, mflush, minimz, minmax, minPrc, MjrPrt, mMajor,
     &     mMinor, MnrPrt, mQNmod, mskip, mNewSB, n, never, nnL,
     &     nout, nParPr, nPr1, nPr2, ObjRow, QPslvr, TPivot
      double precision
     &     bigdx, bigFx, c4, c6, chzbnd, Dens1, Dens2, eps, eps0, eps1,
     &     eps2, eps3, eps4, epsrf, etarg, fdint1, fdint2, Hcndbd,
     &     Lmax1, Lmax2, InfBnd, scltol, small, tCrash, tolCG,
     &     tolCon, tolDcp, tolDdp, tolDpp, tolDrp, tolDup, toldj3,
     &     tolFac, tolFP, tolNLP, tolpiv, tolQP, tolRow, tolSwp, tolUpd,
     &     tolx, Uspace, Utol1, Utol1m, Utol2, Utol2m, viLim, wolfeG,
     &     wtInf0, wtMax, xdlim, xPen0, Zcndbd
!     ------------------------------------------------------------------
      integer            QPChol,     CG,     QN
      parameter         (QPChol = 0, CG = 1, QN = 2)
      integer            LM        , FM
      parameter         (LM     = 0, FM = 1)
      parameter         (cdummy ='-1111111', blank ='        ')
      integer            idummy
      parameter         (idummy = -11111)
      double precision   zero,             one
      parameter         (zero   =  0.0d+0, one    = 1.0d+0)
      double precision   ten
      parameter         (ten    = 10.0d+0)
      double precision   tenp6,            hundrd
      parameter         (tenp6  = 1.0d+6,  hundrd = 100.0d+0)
!     ------------------------------------------------------------------
!     Set some local machine-dependent constants.

      eps       = rw(  1) ! machine precision.  IEEE DP  2.22e-16
      eps0      = rw(  2) ! eps**(4/5)          IEEE DP  3.00e-13
      eps1      = rw(  3) ! eps**(2/3)          IEEE DP  3.67e-11
      eps2      = rw(  4) ! eps**(1/2)          IEEE DP  1.49e-08
      eps3      = rw(  5) ! eps**(1/3)          IEEE DP  6.05e-06
      eps4      = rw(  6) ! eps**(1/4)          IEEE DP  1.22e-04
!     ------------------------------------------------------------------
!     rw(51)--rw(150): optional parameters set via the specs file.
!     ------------------------------------------------------------------
      tolFP     = rw( 51) ! Minor Phase 1 Opt tol
      tolQP     = rw( 52) ! Minor Phase 2 Opt tol
      tolNLP    = rw( 53) ! Major Optimality tolerance
      tolCG     = rw( 54) ! cg tolerance

      tolx      = rw( 56) ! Minor feasibility tolerance.
      tolCon    = rw( 57) ! Major feasibility tolerance.

      tolpiv    = rw( 60) ! excludes small elements of y
      tolrow    = rw( 61) ! tolerance for the row error
      tCrash    = rw( 62) ! crash tolerance
      Utol1m    = rw( 63) ! abs tol for small diag of U in LU mod
      Utol2m    = rw( 64) ! rel tol for small diag of U in LU mod
      tolswp    = rw( 65) ! LU swap tolerance
      tolFac    = rw( 66) ! User-defined LU factor tolerance
      tolUpd    = rw( 67) ! User-defined LU update tolerance
      InfBnd    = rw( 70) ! definition of an infinite bound
      bigFx     = rw( 71) ! unbounded objective
      bigdx     = rw( 72) ! unbounded step
      epsrf     = rw( 73) ! relative function precision.
      fdint1    = rw( 76) ! (1) forwrd diff. interval
      fdint2    = rw( 77) ! (2) cntrl  diff. interval
      xdlim     = rw( 80) ! Step limit
      vilim     = rw( 81) ! violation limit
      etarg     = rw( 83) ! Quasi-Newton QP rg tolerance
      wolfeG    = rw( 84) ! line search tolerance.
      Hcndbd    = rw( 85) ! bound on the condition of Hz
      Zcndbd    = rw( 86) ! bound on the condition of Z
      wtInf0    = rw( 88) ! infeasibility weight
      xPen0     = rw( 89) ! initial penalty parameter.
      wtMax     = rw( 90) ! max     infeasibility weight

      scltol    = rw( 92) ! scale tolerance.
!     ------------------------------------------------------------------
!     rw(151)--rw(180) are parmLU parameters for LUSOL (some optional).
!     ------------------------------------------------------------------
      small     = rw(153) ! defn of small real.
      Utol1     = rw(154) ! abs tol for small diag of U.
      Utol2     = rw(155) ! rel tol for small diag of U.
      Uspace    = rw(156) ! limit on waste space in U.
      Dens1     = rw(157) ! switch to search maxcol columns and no rows.
      Dens2     = rw(158) ! switch to dense LU.
!     ------------------------------------------------------------------
!     rw(181)--rw(199) pass parameters into various routines.
!     ------------------------------------------------------------------
!     toldj3    = rw(186) ! current optimality tol
!     ------------------------------------------------------------------
!     iw(1)--iw(50): I/O file numbers and dimensions.
!     ------------------------------------------------------------------
      iPrint    = iw( 12) ! Print file
!     ------------------------------------------------------------------
!     iw(51)--iw(150): optional parameters set via the specs file.
!     ------------------------------------------------------------------
      maxR      = iw( 52) ! max columns of R.
      maxS      = iw( 53) ! max # of superbasics
      mQNmod    = iw( 54) ! (ge 0) max # of BFGS updates
      QPslvr    = iw( 55) ! 0(1) => QP(QN) QP solver
      lEmode    = iw( 56) ! >0    => use elastic mode
      kchk      = iw( 58) ! check (row) frequency
      kFac      = iw( 59) ! factorization frequency
      ksav      = iw( 60) ! save basis map
      klog      = iw( 61) ! log/print frequency
      kSumm     = iw( 62) ! Summary print frequency
      kDegen    = iw( 63) ! max. expansions of featol
      kReset    = iw( 64) ! Hessian frequency
      mFlush    = iw( 66) ! Hessian flush
      mSkip     = iw( 67) ! # largest value of nSkip
!     lvlSrt    = iw( 69) ! = 0:1:2:3 => cold:basis:warm:hot start
      lvlDer    = iw( 70) ! = 0, 1 or 2, the derivative level
      lvlSys    = iw( 71) ! > 0   => print system info
      lvlHes    = iw( 72) ! 0,1,2  => LM, FM, Exact Hessian
      lvlInf    = iw( 73) ! Elastic option
      lvlScl    = iw( 75) ! scale option
      lvlSch    = iw( 76) ! >0     => use derivatives in the line search
      lvlPre    = iw( 77) ! >0    => QN preconditioned CG
      lvlVer    = iw( 78) ! Verify level
      lvlPPm    = iw( 79) ! 1(2)-norm proximal point method for x0
      lvlPiv    = iw( 80) ! 0(1 2 3) LU partial(rook complete diagonal) pivoting
      lprPrm    = iw( 81) ! > 0    => parms are printed
      lprSch    = iw( 82) ! line search debug starting itn
      lprScl    = iw( 83) ! > 0    => print the scales
      lprSol    = iw( 84) ! > 0    => print the solution
      lprDbg    = iw( 85) ! > 0    => private debug print
      minmax    = iw( 87) ! 1, 0, -1  => MIN, FP, MAX
      iCrash    = iw( 88) ! Crash option
      itnlim    = iw( 89) ! limit on total iterations
      mMajor    = iw( 90) ! limit on major iterations
      mMinor    = iw( 91) ! limit on minor iterations
      MjrPrt    = iw( 92) ! Major print level
      MnrPrt    = iw( 93) ! Minor print level
      nParPr    = iw( 94) ! # of partial pricing sections
      mNewSB    = iw( 95) ! maximum # of new superbasics per major
      jverf1    = iw( 98) ! col # to start derivative checking
      jverf2    = iw( 99) ! col # to stop  derivative checking
      jverf3    = iw(100) ! col # to start derivative checking
      jverf4    = iw(101) ! col # to stop  derivative checking
      ObjRow    = iw(103) ! Objective row of user-defined F
      DerOpt    = iw(104) ! 0, 1, 2 => derivative option
      cgItmx    = iw(111) ! CG iteration limit
      iBack     = iw(120) ! backup file
      iDump     = iw(121) ! dump file
      iLoadB    = iw(122) ! load file
      iNewB     = iw(124) ! new basis file
      iInsrt    = iw(125) ! insert file
      iOldB     = iw(126) ! old basis file
      iPnch     = iw(127) ! punch file
      iReprt    = iw(130) ! report file
      iSoln     = iw(131) ! solution file
!     ------------------------------------------------------------------
!     iw(151)--iw(180) are luparm parameters for LUSOL (some optional).
!     ------------------------------------------------------------------
      maxcol    = iw(153) ! lu1fac: max. # columns
!     ------------------------------------------------------------------
!     Character  workspace.
!     cw(51)--cw(150): optional parameters
!     ------------------------------------------------------------------
      mProb     = cw( 51) ! Problem name
      mObj      = cw( 52) ! Objective name
      mRhs      = cw( 53) ! rhs name
      mRng      = cw( 54) ! range name
      mBnd      = cw( 55) ! bounds name
!     ------------------------------------------------------------------
      c4         = max( 1.0d-4, eps3 )
      c6         = max( 1.0d-6, eps2 )
      never      = 99999999

!     ===============================================================
!     Check the optional parameters.
!     ===============================================================
      if (nnCon .eq. 0) nnJac = 0
      if (nnJac .eq. 0) nnCon = 0
      nnL = max( nnJac, nnObj )

      linCon     = nnCon   .eq. 0
      nlnCon     = nnCon   .gt. 0
      linear     = nnL     .eq. 0
      nonlin     = nnL     .gt. 0

      if (iBack  .eq. idummy ) iBack  =     0
      if (iDump  .eq. idummy ) iDump  =     0
      if (iLoadB .eq. idummy ) iLoadB =     0
      if (iNewB  .eq. idummy ) iNewB  =     0
      if (iInsrt .eq. idummy ) iInsrt =     0
      if (iOldB  .eq. idummy ) iOldB  =     0
      if (iPnch  .eq. idummy ) iPnch  =     0
      if (iReprt .eq. idummy ) iReprt =     0
      if (iSoln  .eq. idummy ) iSoln  =     0

!     Set unspecified frequencies or silly values to defaults.

      if (kchk   .eq. idummy ) kchk   =    60
      if (kFac   .le.    0   ) then
                               kFac   =   100
              if (nonlin     ) kFac   =    50
      end if
      if (klog  .eq. idummy  ) klog   =   100
      if (kSumm .eq. idummy  ) kSumm  =   100
      if (ksav  .eq. idummy  ) ksav   =   100
      if (kDegen.eq. idummy  ) kDegen = 10000
      if (mFlush.eq. idummy  ) mFlush =     0

!     Sometimes, frequency 0 means "almost never".

      if (kchk   .le. 0      ) kchk   = never
      if (mFlush .le. 0      ) mFlush = never
      if (klog   .le. 0      ) klog   = never
      if (ksav   .le. 0      ) ksav   = never
      if (kSumm  .le. 0      ) kSumm  = never
      if (kDegen .le. 0      ) kDegen = never
      if (kReset .le. 0      ) kReset = never

      if (iCrash .lt. 0      ) iCrash =  3

      if (ObjRow .eq.  0) then
                               minimz =  1
                               minmax =  0
      end if
      if (minmax .eq. idummy ) minmax =  1
      if (minmax .eq. -1) then
                               minimz = -1
                         else
                               minimz =  1
                         end if
      if (MjrPrt .eq. idummy ) MjrPrt =  1
      if (MnrPrt .eq. idummy ) MnrPrt =  1

!     if (mMinor .lt. 0      ) mMinor = max( 1000,5*max( n,m ) )
      if (mMinor .lt. 0      ) mMinor = 500
      if (mMajor .lt. 0      ) mMajor = max( 1000,3*max( n,m ) )
      if (mSkip  .lt. 0  .and.  lincon
     &                       ) mSkip  = never
      if (mSkip  .lt. 0  .and.  nlnCon
     &                       ) mSkip  =  2
      if (mNewSB .le. 0      ) mNewSB = 99

      if (lprDbg .lt. 0      ) lprDbg =  0
      if (lprPrm .lt. 0      ) lprPrm =  1
      if (lprSch .lt. 0      ) lprSch = never
      if (lprScl .lt. 0      ) lprScl =  0
      if (lprSol .lt. 0      ) lprSol =  2

!     lvlSrt is checked in s3argA or s3argB
!     if (lvlSrt .lt. 0      ) lvlSrt =  0

      if (DerOpt .ne. idummy ) then
         if (DerOpt .lt. 0  .or.  DerOpt .gt. 1)
     &                         DerOpt =  1
         if (DerOpt .eq. 0   ) lvlDer =  0
         if (DerOpt .eq. 1   ) lvlDer =  3
      else if (lvlDer .ne. idummy ) then
         if (lvlDer .lt. 0  .or.  lvlDer .gt. 3)
     &                         lvlDer =  3
                               DerOpt =  0
         if (lvlDer .eq. 3   ) DerOpt =  1
      else
                               DerOpt =  1
                               lvlDer =  3
      end if

      if (lvlDer .eq. 2  .and.  minmax .eq. 0)
     &                         lvlDer =  3
      if (lvlDer .eq. 0  .and.  minmax .eq. 0)
     &                         lvlDer =  1
      if (lvlSys .lt. 0      ) lvlSys =  0

      if (lvlVer .eq. idummy ) lvlVer =  0
      if (lvlVer .lt. 0      ) lvlVer = -1
      if (lvlVer .gt. 3      ) lvlVer =  0
                               lvlInf =  2
      if (lvlSch .lt. 0      ) then
            if (DerOpt .ne. 1) lvlSch =  0
            if (lvlDer .ne. 3) lvlSch =  0
            if (lvlSch .lt. 0) lvlSch =  1
      end if
      if (lvlPPm .lt. 0      ) lvlPPm =  1
                               lEmode =  1

!     Check superbasics limit maxS and Hessian dimension maxR.

      if ( nonlin ) then
         if (maxR .lt. 0     ) maxR   = min( 2000, nnL+1 )
         if (maxS .lt. 0     ) maxS   =            nnL+1
                               maxR   = max( min( maxR ,n ) , 0 )
                               maxS   = max( min( maxS ,n ) , 1 )
      else ! linear
         if (maxS   .le. 0   ) maxS   = 1
         if (maxR   .le. 0   ) maxR   = 1
      end if

      if (maxS   .lt. maxR   ) maxR   = maxS

      if (QPslvr .lt. 0      ) QPslvr = QPChol
      if (maxR   .eq. 0      ) QPslvr = CG

      if (QPslvr .eq. QN     ) lvlPre = 0
      if (QPslvr .eq. QPChol ) lvlPre = 0
      if (lvlPre .gt. 1      ) lvlPre = 1
      if (lvlPre .lt. 0  .and.  QPslvr .eq. CG)
     &                         lvlPre = 0

      if (cgItmx .lt. 0      ) cgItmx = 100

      if (QPslvr .eq. CG  .or.  maxR   .lt. maxS) then
         if (lvlHes .lt. 0   ) lvlHes = LM
         if (mQNmod .lt. 0   ) mQNmod = 10
      else
         if (lvlHes .lt. 0 .and.  nnL  .gt. 75  )
     &                         lvlHes = LM
         if (lvlHes .lt. 0 .and.  nnL  .le. 75  )
     &                         lvlHes = FM
         if (lvlHes .eq. FM  ) mQNmod = kReset
         if (mQNmod .lt. 0   ) mQNmod = 10
      end if

!     ---------------------------------
!     CG QP optional parameters
!     ---------------------------------
      if (etarg    .lt. zero  .or.
     &    etarg    .gt. one  ) etarg  = 0.1d+0

!     Check other options.

      if (lvlScl .lt. 0   ) then
                               lvlScl = 2
         if ( nonlin )         lvlScl = 1
      end if
                               lvlScl = min( lvlScl, 2 )
      if (lvlScl .eq. 1  .and.  nnL .eq. n)
     &                         lvlScl = 0

      if (nParPr .le. 0   ) then
                               nParPr = 10
         if ( nonlin )         nParPr =  1
      end if
                               minPrc = 10
                               nPr1   = n / nParPr
                               nPr2   = m / nParPr
      if (max( nPr1, nPr2 ) .lt. minPrc) then
                               maxmn  = max( m, n )
                               nParPr = maxmn / min( maxmn, minPrc )
      end if

      cHzbnd = max ( one/(hundrd*eps), tenp6 )
      if (InfBnd   .lt. zero ) InfBnd = 1.0d+20
      if (bigFx    .le. zero ) bigFx  = 1.0d+15
      if (bigdx    .le. zero ) bigdx  = InfBnd
      if (Hcndbd   .le. zero ) Hcndbd = cHzbnd
      if (tCrash   .lt. zero  .or.
     &    tCrash   .ge. one  ) tCrash = 0.1d+0
      if (vilim    .le. zero ) vilim  = 1.0d+6
      if (wolfeG   .lt. zero  .or.
     &    wolfeG   .gt. one  ) wolfeG = 0.9d+0
      if (wtMax    .lt. zero ) wtMax  = 1.0d+10
      if (xdlim    .le. zero ) xdlim  = 2.0d+0
      if (xPen0    .lt. zero ) xPen0  = zero

      if (Zcndbd   .le. zero ) then
          if (QPslvr .eq. QPChol) then
                               Zcndbd = 1.0d+4
          else
                               Zcndbd = 1.0d+6
          end if
      end if

!     ---------------------------------
!     Set up the parameters for lu1fac.
!     ---------------------------------
      if (maxcol .lt.  0     ) maxcol =  5
                               nout   =  iPrint
      if (lvlSys .eq.  0     ) nout   =  0
      if (lvlPiv .le.  0     ) lvlPiv =  0
      if (lvlPiv .gt.  3     ) lvlPiv =  0
                               TPivot =  lvlPiv
      if (linear) then
                               tolDpp =  hundrd
                               tolDrp =  ten
                               tolDcp =  ten
                               tolDdp =  ten
                               tolDup =  ten
      else ! nonlinear
                               tolDpp =  3.99d+0
                               tolDrp =  3.99d+0
                               tolDcp =  3.99d+0
                               tolDdp =  3.99d+0
                               tolDup =  3.99d+0
      end if
      if (tolFac .lt. one    ) then
         if (lvlPiv .eq.   0 ) tolFac =  tolDpp
         if (lvlPiv .eq.   1 ) tolFac =  tolDrp
         if (lvlPiv .eq.   2 ) tolFac =  tolDcp
         if (lvlPiv .eq.   3 ) tolFac =  tolDdp
      end if
      if (tolUpd .lt. one    ) tolUpd =  tolDup
                               Lmax1  =  tolFac
                               Lmax2  =  tolUpd
      if (Utol1    .le. zero ) Utol1  =  eps1
      if (Utol2    .le. zero ) Utol2  =  eps1
      if (Utol1m   .le. zero ) Utol1m =  eps1
      if (Utol2m   .le. zero ) Utol2m =  eps1
      if (Dens2    .lt. zero ) Dens2  =  0.6d+0
      if (small    .le. zero ) small  =  eps0
      if (Uspace   .le. zero ) Uspace =  3.0d+0
      if (Dens1    .le. zero ) Dens1  =  0.3d+0

!     Set some SQP tolerances.
!     Set the minor and major optimality tolerances.
!     Solve the QP subproblems fairly accurately even if the
!     NLP Optimality Tolerance is big.

      if (tolNLP .le. zero) then
                               tolNLP =  2.0d+0*c6
         if (epsrf .gt. zero ) tolNLP =  max( tolNLP, sqrt(ten*epsrf) )
      end if
      if (tolQP    .le. zero ) tolQP  =  min( c6    , tolNLP/2.0d+0   )
      if (tolFP    .lt. zero ) tolFP  =  c6
      if (tolCG    .le. zero ) tolCG  =  1.0d-2
      if (tolrow   .le. zero ) tolrow =  c4
      if (tolswp   .le. zero ) tolswp =  eps4
      if (tolx     .le. zero ) tolx   =  c6
      if (tolCon   .le. eps  ) tolCon =  c6
                               toldj3 =  tolQP
      if (scltol   .le. zero ) scltol =  0.90d+0
      if (scltol   .ge. one  ) scltol =  0.99d+0
      if (tolpiv   .le. zero ) tolpiv =  eps1

      if (linCon) then
         if (wtInf0.lt. zero ) wtInf0 =  1.0d+0
      else
         if (wtInf0.lt. zero ) wtInf0 =  1.0d+4
      end if

      if (epsrf    .le. zero ) epsrf  = eps0
      if (fdint1.le. zero    ) fdint1 = sqrt(epsrf)
      if (fdint2.le. zero    ) fdint2 = epsrf**0.33333d+0

!     Check  START and STOP  column numbers for derivative checking.

      if (jverf1 .le. 0      ) jverf1 = 1
      if (jverf2 .lt. 0      ) jverf2 = n
      if (lvlVer .eq. 2  .or.
     &    lvlVer .eq. 0      ) jverf2 = 0

      if (jverf3 .le. 0      ) jverf3 = 1
      if (jverf4 .lt. 0      ) jverf4 = n
      if (lvlVer .eq. 1  .or.
     &    lvlVer .eq. 0      ) jverf4 = 0

      if (iBack  .eq. iNewB  ) iBack  = 0
      if (itnlim .lt. 0      ) itnlim = max(10000, 10*max(n,m))

!     Set default names (they may be printed by the basis routines).

      if (mProb  .eq. cdummy ) mProb  = blank
      if (mObj   .eq. cdummy ) mObj   = blank
      if (mRhs   .eq. cdummy ) mRhs   = blank
      if (mRng   .eq. cdummy ) mRng   = blank
      if (mBnd   .eq. cdummy ) mBnd   = blank

!     ------------------------------------------------------------------
!     Done.
!     Re-assign the options to their respective work arrays.
!     ------------------------------------------------------------------
      rw( 51)  =  tolFP
      rw( 52)  =  tolQP
      rw( 53)  =  tolNLP
      rw( 54)  =  tolCG
      rw( 56)  =  tolx
      rw( 57)  =  tolCon
      rw( 60)  =  tolpiv
      rw( 61)  =  tolrow
      rw( 62)  =  tCrash
      rw( 63)  =  Utol1m
      rw( 64)  =  Utol2m
      rw( 65)  =  tolswp
      rw( 66)  =  tolFac
      rw( 67)  =  tolUpd
      rw( 70)  =  InfBnd
      rw( 71)  =  bigFx
      rw( 72)  =  bigdx
      rw( 73)  =  epsrf
      rw( 76)  =  fdint1
      rw( 77)  =  fdint2
      rw( 80)  =  xdlim
      rw( 81)  =  vilim
      rw( 83)  =  etarg
      rw( 84)  =  wolfeG
      rw( 85)  =  Hcndbd
      rw( 86)  =  Zcndbd
      rw( 88)  =  wtInf0
      rw( 89)  =  xPen0
      rw( 90)  =  wtMax
      rw( 92)  =  scltol
      rw(151)  =  Lmax1
      rw(152)  =  Lmax2
      rw(153)  =  small
      rw(154)  =  Utol1
      rw(155)  =  Utol2
      rw(156)  =  Uspace
      rw(157)  =  Dens1
      rw(158)  =  Dens2
!     Dependent parameters set in s8dflt.
      rw(181)  =  tolDpp
      rw(182)  =  tolDcp
      rw(183)  =  tolDup
      rw(186)  =  toldj3
      rw(187)  =  tolDrp

!     Addresses for integer quantities.

      iw( 52)  =  maxR
      iw( 53)  =  maxS
      iw( 54)  =  mQNmod
      iw( 55)  =  QPslvr
      iw( 56)  =  lEmode
      iw( 58)  =  kchk
      iw( 59)  =  kFac
      iw( 60)  =  ksav
      iw( 61)  =  klog
      iw( 62)  =  kSumm
      iw( 63)  =  kDegen
      iw( 64)  =  kReset
      iw( 66)  =  mFlush
      iw( 67)  =  mSkip
!     iw( 69)  =  lvlSrt
      iw( 70)  =  lvlDer
      iw( 71)  =  lvlSys
      iw( 72)  =  lvlHes
      iw( 73)  =  lvlInf
      iw( 75)  =  lvlScl
      iw( 76)  =  lvlSch
      iw( 77)  =  lvlPre
      iw( 78)  =  lvlVer
      iw( 79)  =  lvlPPm
      iw( 80)  =  lvlPiv
      iw( 81)  =  lprPrm
      iw( 82)  =  lprSch
      iw( 83)  =  lprScl
      iw( 84)  =  lprSol
      iw( 85)  =  lprDbg
      iw( 87)  =  minmax
      iw( 88)  =  iCrash
      iw( 89)  =  itnlim
      iw( 90)  =  mMajor
      iw( 91)  =  mMinor
      iw( 92)  =  MjrPrt
      iw( 93)  =  MnrPrt
      iw( 94)  =  nParPr
      iw( 95)  =  mNewSB
      iw( 98)  =  jverf1
      iw( 99)  =  jverf2
      iw(100)  =  jverf3
      iw(101)  =  jverf4
      iw(103)  =  ObjRow
      iw(104)  =  DerOpt
      iw(111)  =  cgItmx
      iw(120)  =  iBack
      iw(121)  =  iDump
      iw(122)  =  iLoadB
      iw(124)  =  iNewB
      iw(125)  =  iInsrt
      iw(126)  =  iOldB
      iw(127)  =  iPnch
      iw(130)  =  iReprt
      iw(131)  =  iSoln
      iw(151)  =  nout
      iw(153)  =  maxcol
      iw(156)  =  TPivot
!     Dependent parameters set in s8dflt.
      iw(199)  =  minimz

!     Addresses for character quantities.

      cw( 51)  =  mProb
      cw( 52)  =  mObj
      cw( 53)  =  mRhs
      cw( 54)  =  mRng
      cw( 55)  =  mBnd

      end ! subroutine s8dflt

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Map
     &   ( m, n, negCon, nkx, nnCon, nnJac, nnObj,
     &     lenR, maxR, maxS, mQNmod, lvlHes,
     &     nextcw, nextiw, nextrw, iw, leniw )

      implicit
     &     none
      integer
     &     leniw, lenR, lvlHes, m, maxR, maxS, mQNmod, n, negCon,
     &     nextcw, nextiw, nextrw, nkx, nnCon, nnJac, nnObj, iw(leniw)

*     ==================================================================
*     s8Map   allocates all array storage for snopt,
*     using the values:
*        m    , n    , ne
*        maxS                    Set in s8dflt.
*        nnObj, nnCon, nnJac     Set by specs file or argument list.
*        lenR , negCon           Set in calling program
*
*     On exit,
*        nextcw, nextiw, nextrw are pointers to the next elements of
*                               free space in cw, iw, and rw.
*
*     29 Dec 2000: First version of s8Map.
*     24 Jan 2003: Added workspace for SYMMLQ
*     09 Jan 2005: Current version of s8Map.
*     ==================================================================
      integer
     &     mBS, nnL, nb, ngQP, nlocG, lAscal, lblBS, lblSav, lbuBS,
     &     lbuSav, ldLmul, ldx, lenfR, lfCon, lfCon1, lfCon2, lfR,
     &     lgCon, lgCon1, lgCon2, lgConu, lgObj1, lgObj2, lgObju, lFv,
     &     lFx, lgBS, lgObj, ldg, lgQP, lH0, lHd, lHdx,
     &     lhElas, lhEsta, lhfeas, liy, liy1, lkBS, lkx, lLmul,
     &     lLmul1, lLmul2, llocG, lQPrhs, lR, lr1, lr2, lRdx, lrg, lS,
     &     ls1, ls2, ls3, lV, lx1, lxBS, lxscal, lxPen, lxQP, lxQP0,
     &     ly, ly1, ly2, ly3, ly4
*     ------------------------------------------------------------------
      integer            LM   ,      FM
      parameter         (LM     = 0, FM     = 1)
*     ------------------------------------------------------------------
      nnL     = max( nnJac, nnObj )

*     All dimensions are computed from
*        m     , n    , ne
*        lenR  , maxS , mQMmod
*        nnObj , nnCon, nnJac
*        negCon

      ngQP    = nnL
      mBS     = m     + maxS
      nb      = n     + m

*     Nonlinear constraints.

      nlocG   = nnJac  + 1

*     Addresses for the integer arrays.

      lkx    = nextiw
      lhfeas = lkx    + nkx
      lkBS   = lhfeas + mBS
      lhEsta = lkBS   + mBS
      lhElas = lhEsta + nb
      liy    = lhElas + nb
      liy1   = liy    + nb
      nextiw = liy1   + nb

*     Addresses for the double precision arrays.

      lAscal = nextrw
      ly     = lAscal + nb
      ly1    = ly     + nb
      ly2    = ly1    + nb
      ly3    = ly2    + nb
      ly4    = ly3    + nb      ! SYMMLQ workspace
      ls1    = ly4    + nb      ! SYMMLQ workspace
      ls2    = ls1    + maxS    ! SYMMLQ workspace
      ls3    = ls2    + maxS    ! SYMMLQ workspace
      lr1    = ls3    + maxS    ! SYMMLQ workspace
      lr2    = lr1    + maxS    ! SYMMLQ workspace
      lblBS  = lr2    + maxS
      lbuBS  = lblBS  + mBS
      lxBS   = lbuBS  + mBS
      lxScal = lxBS   + mBS
      lgBS   = lxscal + nnL
      lgQP   = lgBS   + mBS
      lRdx   = lgQP   + ngQP
      lHdx   = lRdx   + nnL
      lH0    = lHdx   + nnL
      ldg    = lH0    + nnL
      lR     = ldg    + nnL
      lrg    = lR     + lenR
      lblSav = lrg    + maxS
      lbuSav = lblSav + nb
      nextrw = lbuSav + nb

*     Nonlinear Objective.

      lgObj  = nextrw
      lgObj1 = lgObj  + nnObj
      lgObj2 = lgObj1 + nnObj
      lgObju = lgObj2 + nnObj
      nextrw = lgObju + nnObj

*     Nonlinear constraints.

      llocG  = nextiw
      nextiw = llocG  + nlocG

      lfCon  = nextrw
      lfCon1 = lfCon  + nnCon
      lfCon2 = lfCon1 + nnCon
      lFx    = lfCon2 + nnCon
      lFv    = lFx    + nnCon
      lLmul  = lFv    + nnCon
      lLmul1 = lLmul  + nnCon
      lLmul2 = lLmul1 + nnCon
      ldLmul = lLmul2 + nnCon
      lxPen  = ldLmul + nnCon
      lgCon  = lxPen  + nnCon
      lgCon1 = lgCon  + negCon
      lgCon2 = lgCon1 + negCon
      lgConu = lgCon2 + negCon
      lQPrhs = lgConu + negCon
      ldx    = lQPrhs + m
      lxQP   = ldx    + nb
      lxQP0  = lxQP   + nb
      lx1    = lxQP0  + nb
      nextrw = lx1    + nb

*     Store the addresses in iw.

      iw(251) = lkx
      iw(260) = llocG

      iw(273) = lblBS
      iw(274) = lbuBS
      iw(275) = lblSav
      iw(276) = lbuSav

      iw(278) = lQPrhs

      iw(283) = lhElas
      iw(284) = lhfeas
      iw(285) = lhEsta

      iw(287) = ldx
      iw(288) = lHdx
      iw(289) = ldg
      iw(290) = lgQP
      iw(291) = lgBS
      iw(292) = lkBS
      iw(293) = lrg
      iw(294) = lR
      iw(295) = lAscal
      iw(296) = lgObj

      iw(300) = lx1
      iw(301) = lxBS
      iw(302) = lxscal

      iw(304) = lxPen
      iw(305) = lxQP
      iw(306) = lxQP0

      iw(308) = liy
      iw(309) = liy1
      iw(311) = ly
      iw(312) = ly1
      iw(313) = ly2
      iw(314) = ly3
      iw(315) = ly4

      iw(316) = lfCon
      iw(317) = lfCon1
      iw(318) = lfCon2
      iw(319) = lgConu
      iw(320) = lgCon
      iw(321) = lgCon1
      iw(322) = lgCon2
      iw(323) = lgObju
      iw(324) = lgObj1
      iw(325) = lgObj2

      iw(336) = lFx
      iw(337) = lFv
      iw(345) = lRdx
      iw(346) = lH0

      iw(348) = lLmul
      iw(349) = lLmul1
      iw(350) = lLmul2
      iw(351) = ldLmul

      iw(353) = lr1
      iw(354) = lr2
      iw(355) = ls1
      iw(356) = ls2
      iw(357) = ls3

*     Allocate space for an approximate Hessian.
*     The amount will depend on the method selected.

      if (lvlHes .eq. LM) then
*        ---------------------------------------------------------------
*        Compute the addresses of the limited-memory arrays.
*        These are saved and used for subsequent entries.
*        ---------------------------------------------------------------
         lHd     = nextrw
         lS      = lHd    + nnL
         lV      = lS     + nnL*mQNmod
         nextrw  = lV     + nnL*mQNmod

         iw(347) = lHd
         iw(401) = lS
         iw(402) = lV

      else if (lvlHes .eq. FM) then
         lenfR   = nnL*(nnL + 1)/2

         lHd     = nextrw
         lfR     = lHd    + nnL
         nextrw  = lfR    + lenfR

         iw(347) = lHd
         iw(391) = lfR
         iw(392) = lenfR
      end if

      end ! subroutine s8Map

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8solv
     &   ( iExit, Solver, Start,
     &     fgwrap, fgcon, fgobj,
     &     MjrLog, MnrLog, snSTOP, gotR,
     &     m, n, nb, nnCon, nnJac, nnObj,
     &     nName, iObj, ObjAdd, fObj, ObjTru, nInf, sInf,
     &     ne, nlocJ, locJ, indJ, Jcol,
     &     bl, bu, Names,
     &     hs, x, pi, rc, nMajor, nS,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     fgwrap, fgcon, fgobj, MjrLog, MnrLog, snSTOP
      logical
     &     gotR
      integer
     &     iExit, iObj, lencu, lencw, leniu, leniw, lenru, lenrw, m, n,
     &     nb, ne, nlocJ, nInf, nName, nnCon, nnJac, nnObj, nS, Start,
     &     locJ(nlocJ), indJ(ne), hs(nb), iu(leniu), iw(leniw)
      double precision
     &     ObjAdd, fObj, sInf, Jcol(ne), bl(nb), bu(nb),
     &     x(nb), pi(m), rc(nb), ru(lenru), rw(lenrw)
      character
     &     Solver*6, Names(nName)*8, cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     s8solv solves the current problem.
*
*     On entry,
*     the specs file has been read,
*     all data items have been loaded (including locJ, indJ, Jcol, ...),
*     and workspace has been allocated.
*
*     On exit,
*     iExit  =  0 if an optimal solution was found,
*            =  1 if the problem was infeasible,
*            =  2 if the problem was unbounded,
*            =  3 if the Iteration limit was exceeded,
*           ge  4 if iterations were terminated by some other
*                 error condition (see the SNOPT user's guide).
*
*     15 Nov 1991: First version based on Minos 5.4 routine misolv.
*     13 Feb 1994: Eliminated "Cycle" options.
*                Simplified s4getb.
*     12 Nov 1994: Integer workspace added.
*     25 Jul 1996: Sign of the slacks changed.
*     28 Sep 1997: Character workspace added.
*     11 Nov 1997: Backtracking for undefined functions.
*     26 Dec 1997: Dummy Jacobian scaled in feasibility phase.
*     27 Aug 1998: Constant Jacobian elements handled correctly.
*     10 Oct 1998: Objective and constraint gradient checking merged.
*     11 Oct 1998: Facility to combine funobj and funcon added.
*     23 Dec 1999: Suboptimize option added.
*     30 Dec 2000: Housekeeping for first function call moved to snwrap.
*     03 Aug 2003: snEXIT and snPRNT adopted.
*     19 Mar 2006: Fx initialized for s8savB.
!     22 Apr 2007: Call-status made global
*     ==================================================================
      character
     &     istate(3)*4, mProb*8, str*133, str2*133
      external
     &     dnrm1s
      logical
     &     FPonly, GotFun, linInf, nlnInf, needB, nlnCon, nlnObj, nonlin
      integer
     &     cgItn, cgItns, eigH, gotG, iCrash, inewB, inform, itn,
     &     itnlim, j, k, lenR, lenx0, lCrash, lsSave, lvlDer, lvlDif,
     &     lvlHes, lvlScl, lvlSch, lAscal, lblBS, lbuBS, lblSav,
     &     lbuSav, ldlmul, ldx, lFv, lFx, lfCon, lfCon1, lfCon2, lgCon,
     &     lgCon1, lgCon2, lgConu, lgObj1, lgObj2, lgBS, lgObj, lgQP,
     &     ldg, lHdx, lhElas, lhEsta, lhfeas, linesL,
     &     linesS, liy, liy1, lkBS, lkx, lLmul, lLmul1, lLmul2, llocG,
     &     lQPrhs, lR, lrg, lUdx, lx1, lxBS, lxPen, lxQP, lxQP0,
     &     ly, ly1, ly2, ly3, maxS, maxvi, mBS, minimz, minmax,
     &     MjrPrt, MnrPrt, modefg, nDegen, nFac, nfCon1,
     &     nfCon2, nfCon3, nfCon4, nfObj1, nfObj2, nfObj3, nfObj4,
     &     negCon, nkx, nlocG, nMajor, nMinor, nnb, nnL0, nnL, nnObj0,
     &     nnObj1, nnCon0, nnCon1, numLC, nrhs0, nrhs,
     &     numLIQ, nx0, Status
      double precision
     &     Degen, dnrm1s, duInf, eps0, fLin, fMrt, ObjTru,
     &     PenNrm, InfBnd, piNorm, pNorm1, pNorm2, rgNorm, sclObj,
     &     xNorm, tolFP, tolQP, tolx, tCrash, viLim, viMax,
     &     viRel, viSup, wtInf, wtInf0
*     ------------------------------------------------------------------
      double precision   zero,            one,          ten
      parameter         (zero   = 0.0d+0, one = 1.0d+0, ten = 10.0d+0)
      integer            SEMDEF,     POSDEF
      parameter         (SEMDEF = 0, POSDEF = 1)
      integer            Scale,      UnScal
      parameter         (Scale  = 0, UnScal = 1)
      integer            RowTyp,     Stats
      parameter        ( Rowtyp = 0, Stats  = 1 )
      integer            Wrap
      parameter         (Wrap   = 1)
      integer            SaveB,      PrintS
      parameter         (SaveB  = 0, PrintS = 1)
      integer            LM        , FM
      parameter         (LM     = 0, FM     = 1)
      integer            xBound,     xMove
      parameter         (xBound = 0, xMove  = 1)

      parameter         (lvlDer =  70) ! = 0,1,2 or 3, deriv level
      parameter         (lvlHes =  72) ! 0,1,2  => LM, FM, Newton
      parameter         (lvlScl =  75) ! scale option
      parameter         (minmax =  87) ! 1, 0, -1  => MIN, FP, MAX
      parameter         (eigH   = 200) ! =1(0) for pd  QP Hessian
      parameter         (lvlDif = 181) ! =1(2) for forwd(cntrl) diffs
      parameter         (gotG   = 184) ! > 0 => some exact derivs
      parameter         (nfCon1 = 189) ! number of calls of fCon
      parameter         (nfCon2 = 190) ! number of calls of fCon
      parameter         (nfCon3 = 191) ! number of calls of fCon
      parameter         (nfCon4 = 192) ! number of calls of fCon
      parameter         (nfObj1 = 194) ! number of calls of fObj
      parameter         (nfObj2 = 195) ! number of calls of fObj
      parameter         (nfObj3 = 196) ! number of calls of fObj
      parameter         (nfObj4 = 197) ! number of calls of fObj
      parameter         (linesL = 220) ! # lines in log     file
      parameter         (linesS = 221) ! # lines in summary file
      parameter         (nFac   = 210) ! # of LU factorizations
      parameter         (cgItns = 386) ! Number of symmlq iterations
      parameter         (cgItn  = 387) ! symmlq itns for last minor
*     ------------------------------------------------------------------
      iNewB     = iw(124) ! new basis file

      negCon    = iw( 20) ! # of nonzero elems in J
      lenR      = iw( 28) ! R(lenR) is the reduced Hessian factor
      maxS      = iw( 53) ! max # of superbasics
      lvlSch    = iw( 76) ! >0     => use derivatives in the line search

      iCrash    = iw( 88) ! Crash option
      itnlim    = iw( 89) ! limit on total iterations

      MjrPrt    = iw( 92) ! Major print level
      MnrPrt    = iw( 93) ! Minor print level
      minimz    = iw(199) ! 1 (-1)    => minimize (maximize)
      nkx       = iw(247) ! dimension of kx and its inverse, kxN

      eps0      = rw(  2)
      tolFP     = rw( 51) ! Minor Phase 1 Opt tol
      tolQP     = rw( 52) ! Minor Phase 2 Opt tol
      tolx      = rw( 56) ! Minor feasibility tolerance.
      tCrash    = rw( 62) ! crash tolerance.
      InfBnd    = rw( 70) ! definition of an infinite bound.
      vilim     = rw( 81) ! violation limit
      wtInf0    = rw( 88) ! infeasibility weight

      mProb     = cw( 51) ! Problem name

*     Addresses

      lkx    = iw(251) ! j  = kx (jN) => col j of Jcol is variable jN
      llocG  = iw(260) ! locG(nnJac+1) = column pointers for indG
      lhElas = iw(283) ! hElast(nb) list of elastic vars
      lhfeas = iw(284) ! hfeas(mBS), feasibility types
      lhEsta = iw(285) ! hEstat(nb), status of elastics
      lkBS   = iw(292) ! kBS(mBS), ( B  S ) list
      lblBS  = iw(273) ! blBS(mBS)   = lower bounds for xBS
      lbuBS  = iw(274) ! buBS(mBS)   = upper bounds for xBS
      lxBS   = iw(301) ! xBS(mBS)    = basics, superbasics
      lgQP   = iw(290) ! gQP(ngQP)   = QP gradient
      lgBS   = iw(291) ! gBS(mBS)    = BS components of g
      lAscal = iw(295) ! Ascale(nb)  = row and column scales
      lgObj  = iw(296) ! gObj(nnObj) = Objective gradient
      lx1    = iw(300) ! x1(nb)      = new x, used to store x0
      lrg    = iw(293) ! rg(maxS)    = reduced gradient
      lR     = iw(294) ! R(lenR)     = factor of Z'HZ
      liy    = iw(308) ! iy(nb)      =  integer work vector
      liy1   = iw(309) ! iy1(nb)     =  integer work vector
      ly     = iw(311) ! y(nb)       =  real work vector
      ly1    = iw(312) ! y1(nb)      =  real work vector
      ly2    = iw(313) ! y2(nb)      =  real work vector
      ly3    = iw(314) ! y3(nb)      =  real work vector
      lQPrhs = iw(278) ! QPrhs(nnCon)=  QP constraint rhs
      ldx    = iw(287) ! dx(nb)      = x1 - x
      lHdx   = iw(288) ! Hdx(nnL)    = product of H with  x1 - x
      ldg    = iw(289) ! dg(nnL)     = gradient difference
      lfCon  = iw(316) ! fCon (nnCon) constraints at x
      lfCon1 = iw(317) ! fCon1(nnCon) constraints at x1
      lfCon2 = iw(318) ! fCon2(nnCon) work vector
      lgConu = iw(319) ! record of unknown derivatives and constants
      lgCon  = iw(320) ! gCon (negCon)   constraint gradients at x
      lgCon1 = iw(321) ! gCon1(negCon)   constraint gradients at x1
      lgCon2 = iw(322) ! gCon2(negCon)   work vector
      lgObj1 = iw(324) ! gObj1(nnObj) objective gradients at x1
      lgObj2 = iw(325) ! gObj2(nnObj) work gObj

      lFx    = iw(336) ! Fx (nnCon)  = F(x) + A(linear)x
      lFv    = iw(337) ! Fv          = F(x) + A(linear)x - sN
      lblSav = iw(275) ! blSav(nb)   = copy of bl
      lbuSav = iw(276) ! buSav(nb)   = copy of bu

      lUdx   = iw(345) ! Udx(nnL)      = product of U with dx
      lLmul  = iw(348) ! Lmul (nnCon)  = multipliers for F
      lLmul1 = iw(349) ! Lmul1(nnCon)  = Lmul at x1
      lLmul2 = iw(350) ! Lmul2(nnCon)  = work copy of Lmul
      ldLmul = iw(351) ! dLmul(nnCon)  = Lmul1 - Lmul
      lxPen  = iw(304) ! xPen(nnCon)   = penalty params
      lxQP   = iw(305) ! xQP(nb)     = QP solution
      lxQP0  = iw(306) ! xQP0(nb)    = QP feasible pt.

      iExit  = 0

      GotFun = .false.
      FPonly = iw(minmax) .eq. 0

      nnL    = max( nnJac, nnObj )
      nlnCon = nnCon  .gt. 0
      nlnObj = nnObj  .gt. 0
      nonlin = nnL    .gt. 0

      nnObj0 = max( nnObj, 1 )
      nnCon0 = max( nnCon, 1 )
      nnL0   = max( nnL  , 1 )
      nx0    = nb
      lenx0  = nb
      mBS    = m     + maxS
      nlocG  = nnJac + 1

      numLC  = m - nnCon

*     Initialize Lmul from pi.
*     Zap the pi(i) to prevent them being printed without being set.

      if (nlnCon      ) call dcopy ( nnCon,       pi, 1, rw(lLmul), 1 )
      if (numLC .gt. 0) call dload ( numLC, zero, pi(nnCon+1), 1 )

*     Initialize a few things.
*     Define the Hessian type for the QP subproblem.

      if (iw(lvlHes) .eq. LM  .or.  iw(lvlHes) .eq. FM) then
         if (nnL .lt. n) then
            iw(eigH) = SEMDEF
         else
            iw(eigH) = POSDEF
         end if
      end if

      iw(lvlDif) = 1
      iw(nFac)   = 0
      iw(cgItns) = 0
      iw(cgItn ) = 0
      nInf       = 0
      wtInf      = one

      duInf      = zero
      fMrt       = zero
      fObj       = zero
      ObjTru     = zero
      PenNrm     = zero
      piNorm     = zero
      vimax      = zero
      virel      = zero

      iw(linesL) = 0
      iw(linesS) = 0

      call iload ( 4, 0, iw(nfCon1), 1 )
      call iload ( 4, 0, iw(nfObj1), 1 )

      itn        = 0
      nDegen     = 0
      nMajor     = 0
      nMinor     = 0

      call s1page( 1, iw, leniw )

*     ------------------------------------------------------------------
*     Print the matrix statistics before the nonlinear part of Jcol is
*     loaded with random elements.  The rowtypes are also computed ready
*     for use in s5getB.
*     ------------------------------------------------------------------
      call s2Amat
     &   ( Stats, MjrPrt, m, n, nb,
     &     nnCon, nnJac, nnObj, iObj,
     &     ne, nlocJ, locJ, indJ, Jcol,
     &     bl, bu,  iw(lhEsta),
     &     iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     Make a permanent copy in gConu of the constant Jacobian elements
*     stored in J.  Load the nonlinear part of J with random elements.
*     ------------------------------------------------------------------
      if ( nlnCon ) then
         call s8Gcpy
     &      ( nnCon, nnJac, ne, nlocJ, locJ, indJ,
     &        ne    , nlocJ,     locJ,       Jcol,
     &        negCon, nlocG, iw(llocG), rw(lgConu) )
         call s8rand
     &      ( negCon, negCon, rw(lgCon) )
         call s8Gcpy
     &      ( nnCon, nnJac, ne, nlocJ, locJ, indJ,
     &        negCon, nlocG, iw(llocG), rw(lgCon),
     &        ne    , nlocJ,     locJ ,     Jcol )
      end if

      nrhs  = 0                 ! No rhs when finding the first basis
      nrhs0 = 1

      call iload
     &   ( nb, 0, iw(lhElas), 1 )
      call s5getB
     &   ( inform, Start, MnrLog, needB, m, maxS, mBS,
     &     n, nb, nnCon, nnJac, nnObj, nName, nS, nMinor, itnlim, itn,
     &     nDegen, numLC, numLIQ, tolFP, tolQP, tolx,
     &     nInf, sInf, wtInf, iObj, sclObj, piNorm, rgNorm,
     &     ne, nlocJ, locJ, indJ, Jcol,
     &     iw(lhElas), iw(lhEsta), iw(lhfeas), hs, iw(lkBS), Names,
     &     rw(lAscal), bl, bu,rw(lblBS),rw(lbuBS),rw(lblSav),rw(lbuSav),
     &     rw(lgBS), pi, rc, nrhs0, nrhs, rw(lQPrhs),
     &     lenx0, nx0, rw(lx1), x, rw(lxBS),
     &     iw(liy), iw(liy1), rw(ly), rw(ly1), rw(ly2), rw(ly3),
     &     cw, lencw, iw, leniw, rw, lenrw )

*     Possible inform values are = -3,-2,-1, 0, >0

      if (inform .ne. 0) then
         if (inform .gt. 0) then
            iExit = inform      ! fatal error
         else if (inform .eq. -3) then
            iExit = 31          ! too many iterations
         else
            iExit = 12          ! infeasible linear equalities
         end if
      end if

*     ==================================================================
*     Satisfy the linear constraints.
*     The norm of x is minimized via a proximal-point QP.
*     If no feasible point can be found, the linear rows can be elastic.
*     ==================================================================
      if (numLC .gt. 0) then
         if (iExit .eq. 0) then
            call s8feas
     &         ( iExit, MnrLog, lenR, m, maxS, mBS,
     &           n, nb, nnCon0, nnCon, nnL0, nnL, nDegen, nS,
     &           numLC, numLIQ, itn, itnlim, nMinor, MnrPrt, sclObj,
     &           tolQP, tolx, nInf, sInf, wtInf, piNorm, rgNorm,
     &           ne, nlocJ, locJ, indJ, Jcol,
     &           iw(lhElas), iw(lhEsta), iw(lhfeas), hs, iw(lkBS),
     &           rw(lAscal), bl, bu, rw(lblSav), rw(lbuSav), rw(lblBS),
     &           rw(lbuBS), rw(lgBS), pi, rw(lR), rc, rw(lQPrhs),
     &           rw(lx1), x, rw(lxBS),
     &           iw(liy), iw(liy1), rw(ly), rw(ly1), rw(ly2), rw(ly3),
     &           cw, lencw, iw, leniw, rw, lenrw )
         end if
*        ---------------------------------------------------------------
*        Reinstate the scaled bounds on the nonlinear constraints.
*        ---------------------------------------------------------------
         if ( nlnCon ) then
            call dcopy ( nnCon, rw(lblSav+n), 1, bl(n+1), 1 )
            call dcopy ( nnCon, rw(lbuSav+n), 1, bu(n+1), 1 )
         end if ! nlnCon

*        ---------------------------------------------------------------
*        Unscale the linear constraints.
*        ---------------------------------------------------------------
         if (iw(lvlScl) .gt. 0) then
            call s2scla
     &         ( UnScal, m, n, nb, iObj, InfBnd, sclObj,
     &           ne, nlocJ, locJ, indJ, Jcol,
     &           rw(lAscal), bl, bu, pi, x )
         end if
      end if ! numLC > 0

*     We are done if the linear constraints are infeasible.

      linInf = iExit .ne. 0
      if ( linInf ) go to 900

*     ------------------------------------------------------------------
*     Copy the constant Jacobian elements in gCon, gCon1 and gCon2.
*     Reset hElast so that only nonlinear rows are elastic.
*     Make sure variables are not outside their bounds
*     (in particular, check the nonlinear slacks).
*     ------------------------------------------------------------------
      if ( nlnCon ) then
         call dcopy ( negCon, rw(lgConu), 1, rw(lgCon ), 1 )
         call dcopy ( negCon, rw(lgConu), 1, rw(lgCon1), 1 )
         call dcopy ( negCon, rw(lgConu), 1, rw(lgCon2), 1 )
         call iload ( nnCon, 3, iw(lhElas+n), 1 )
      end if ! nlnCon

      call s5FixX( xBound, 1, nb, tolx, hs, bl, bu, x )

*     ==================================================================
*     ==================================================================
*     The linear constraints have been satisfied!
*     Compute the problem functions at this all-important point.
*     No scaling yet.
*     ==================================================================
*     ==================================================================
      if (nnL .gt. 0) then
         lsSave     = iw(lvlScl)
         iw(lvlScl) = 0

         Status = 1
         modefg = 2
         call fgwrap
     &      ( inform, modefg, Status, nlnCon, nlnObj,
     &        n, negCon, nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj,
     &        fgcon, fgobj, x,
     &        ne, nlocJ, locJ, indJ,
     &        rw(lfCon), fObj, rw(lgCon), rw(lgObj),
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )
         GotFun = inform .eq. 0

         if (.not. GotFun) then
            if (inform .lt. 0) then
               if (numLC .gt. 0) then
                  iExit = 61    ! Undefined fun at first feasible point
               else
                  iExit = 62    ! Undefined fun at initial point
               end if
            else
               iExit = inform   ! User wants to stop
            end if
            go to 900
         end if

         Status = 0

         if ( nlnCon ) then     ! Define Fx for s8savB
            call dcopy
     &         ( nnCon, rw(lfCon), 1, rw(lFx), 1 )
         end if

*        ---------------------------------------------------------------
*        Check derivatives.
*        (One day, we will do this on the SCALED problem.)
*        ---------------------------------------------------------------
         call s7chkG
     &      ( iExit, n, nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj,
     &        fgwrap, fgcon, fgobj,
     &        x, rw(lx1), bl, bu, fObj, rw(lgObj),
     &        ne, nlocJ, locJ, indJ, negCon, nlocG, iw(llocG),
     &        rw(lfCon), rw(lgCon), rw(lgObj2), rw(lfCon2), rw(lgCon2),
     &        rw(ly), rw(ly1), cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )
         if (iExit .ne. 0) go to 900

*        ---------------------------------------------------------------
*        Compute any missing derivatives.
*        Load the Jacobian gCon in  J.
*        ---------------------------------------------------------------
         call s6fdG
     &      ( iExit, n, negCon,
     &        nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj,
     &        fgwrap, fgcon, fgobj,
     &        bl, bu, x,
     &        ne, nlocJ, locJ, indJ,
     &        rw(lfCon), fObj, rw(lgCon), rw(lgObj), rw(ly3),
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )
         if (iExit .ne. 0) go to 900

         if ( nlnCon ) then
            call s8Gcpy
     &         ( nnCon, nnJac, ne, nlocJ, locJ, indJ,
     &           negCon, nlocG, iw(llocG), rw(lgCon),
     &           ne, nlocJ, locJ, Jcol )
            call dcopy
     &         ( nnCon, rw(lfCon), 1, rw(lFx), 1 )
         end if
         iw(lvlScl) = lsSave
      end if

*     ==================================================================
*     Scale the problem.
*     ==================================================================
      if (iw(lvlScl) .gt. 0) then
*        ---------------------------------------------------------------
*        Reset the vector of row types.
*        ---------------------------------------------------------------
         call s2Amat
     &      ( RowTyp, MjrPrt, m, n, nb,
     &        nnCon, nnJac, nnObj, iObj,
     &        ne, nlocJ, locJ, indJ, Jcol,
     &        bl, bu, iw(lhfeas),
     &        iw, leniw, rw, lenrw )
         call s2scal
     &      ( MjrPrt, m, n, nb, nnL, nnCon, nnJac, iw(lhfeas),
     &        ne, nlocJ, locJ, indJ, Jcol,
     &        rw(lAscal), bl, bu, rw(ly), rw(ly2),
     &        iw, leniw, rw, lenrw )
         call s2scla
     &      ( Scale, m, n, nb, iObj, InfBnd, sclObj,
     &        ne, nlocJ, locJ, indJ, Jcol,
     &        rw(lAscal), bl, bu, pi, x )

*        ---------------------------------------------------------------
*        The objective and constraint functions haven't been scaled yet.
*        Scale the constant elements in gCon1 and gCon2.
*        Don't forget the initial pi.
*        ---------------------------------------------------------------
         if ( nlnCon ) then
            call dddiv
     &         ( nnCon, rw(lAscal+n), 1, rw(lfCon), 1 )
            if (iw(gotG) .gt. 0) then
               call s8sclJ
     &            ( nnCon, nnJac, negCon, n, rw(lAscal),
     &              ne, nlocJ, locJ, indJ, rw(lgCon), rw, lenrw )
               call dcopy
     &            ( negCon, rw(lgCon), 1, rw(lgCon1), 1 )
               call dcopy
     &            ( negCon, rw(lgCon), 1, rw(lgCon2), 1 )
            end if
            call ddscl
     &         ( nnCon, rw(lAscal+n), 1, rw(lLmul), 1 )
         end if

         if (nlnObj  .and.  iw(gotG) .gt. 0) then
            call s8sclg
     &         ( nnObj, rw(lAscal), rw(lgObj), rw, lenrw )
         end if
      end if ! iw(lvlScl) > 0

*     ==================================================================
*     s8Fx computes the nonlinear constraint values Fx.
*     Copy these into the slacks x(n+i) and make sure they are feasible.
*     Crash uses them to decide which slacks to grab for the basis
*     If any nonbasic nonlinear slacks are close to a bound,
*     move them exactly onto the bound to avoid very small steps.
*     ==================================================================
      if ( nlnCon ) then
         call s8Fx
     &      ( n, nnCon, nnJac, eps0,
     &        ne, nlocJ, locJ, indJ, Jcol, rw(lfCon), x, rw(lFx) )
         call s2vmax
     &      ( n, nnCon, maxvi, vimax, bl, bu, rw(lFx) )
         viSup = max( ten*vimax, vilim )
         call dcopy ( nnCon, rw(lFx), 1, x(n+1), 1 )

         call s5FixX
     &      ( xMove, n+1, n+nnCon, tolx, hs, bl, bu, x )

*        ===============================================================
*        Crash on the nonlinear rows.
*        hs(*) already defines a basis for the full problem,  but we
*        want to do better by not including all of the slacks.
*        ===============================================================
         if ( needB ) then

*           Load  hfeas  with the row types.
*           s2crsh uses kBS as workspace.  It may alter x(n+i) for
*           nonlinear slacks.

            call s2Amat
     &         ( RowTyp, MjrPrt, m, n, nb,
     &           nnCon, nnJac, nnObj, iObj,
     &           ne, nlocJ, locJ, indJ, Jcol,
     &           bl, bu, iw(lhfeas),
     &           iw, leniw, rw, lenrw )
            lcrash = 5
            call s2crsh
     &         ( lcrash, MjrPrt, m, n, nb, nnCon,
     &           iCrash, tCrash,
     &           ne, nlocJ, locJ, indJ, Jcol,
     &           iw(lkBS), hs, iw(lhfeas), bl, bu, x,
     &           iw, leniw, rw, lenrw )
            needB = .false.
         end if ! needB
      end if ! nlnCon

*     ------------------------------------------------------------------
*     Solve the problem.
*     ------------------------------------------------------------------
      call s1page( 1, iw, leniw )
      call s1time( 2, 0, iw, leniw, rw, lenrw )
      call s8SQP
     &   ( iExit, fgwrap, fgcon, fgobj,
     &     MjrLog, MnrLog, snSTOP, gotR,
     &     itn, lenR, m, maxS, mBS, n, nb, nS,
     &     nnCon0, nnCon, nnObj0, nnObj, nnL0, nnL,
     &     nMajor, nMinor, nDegen, duInf,
     &     minimz, iObj, sclObj, ObjAdd, fObj, fMrt,
     &     vimax, virel, viSup, nInf, sInf,
     &     wtInf0, wtInf, PenNrm, piNorm, xNorm,
     &     ne, nlocJ, locJ, indJ, Jcol, negCon, nlocG, iw(llocG),
     &     iw(lhElas), iw(lhEsta), iw(lhfeas), hs, iw(lkBS),
     &     rw(lAscal), bl, bu, rw(lblBS), rw(lbuBS), rw(lFv), rw(lFx),
     &     rw(lfCon), rw(lgCon), rw(lgObj),
     &     rw(lfCon1), rw(lgCon1), rw(lgObj1),
     &     rw(lfCon2), rw(lgCon2), rw(lgObj2),
     &     rw(lgBS), rw(lgQP), rw(ldLmul), rw(ldx), rw(ldg), rw(lUdx),
     &     rw(lHdx), rw(lLmul), rw(lLmul1), rw(lLmul2), pi, rw(lQPrhs),
     &     rw(lR), rc, rw(lrg), x, rw(lx1), rw(lxBS),
     &     rw(lxQP0), rw(lxQP), rw(lxPen),
     &     iw(liy), iw(liy1), rw(ly), rw(ly1), rw(ly2), rw(ly3),
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call s1time(-2, 0, iw, leniw, rw, lenrw )

*     ==================================================================
*     Exit.
*     Set output variables and print a summary of the final solution.
*     ==================================================================
  900 call snWRAP( iExit, Solver, str, str2, iw, leniw )

      degen  = 100.0d+0 * nDegen / max( itn, 1 )

      if (iObj .eq. 0) then
         flin = ObjAdd
      else
         flin = ObjAdd + x(n+iObj)*sclObj
      end if
      ObjTru  = flin + fObj

      nlnInf  = nInf .gt. 0
      xNorm   = dnrm1s( n , x, 1 )

*     Count basic nonlinear variables (used only for printing).

      nnb    = 0
      do j = 1, nnL
         if (hs(j) .eq. 3) nnb = nnb + 1
      end do

      if (inewB .gt. 0  .and.  iExit/10 .lt. 8) then
         k      = 1 + iExit/10
         call s4stat
     &      ( k, istate )
         call s4newB
     &      ( Wrap, iNewB, minimz, m, n, nb,
     &        nS, mBS, itn, nInf, sInf, fObj, iw(lkBS), hs,
     &        rw(lAscal), bl, bu, x, rw(lxBS), istate,
     &        cw, lencw, iw, leniw )
      end if

*     Print statistics.

      call snPRNT(13,
     &     ' Problem name                 '//mProb, iw, leniw )
      write(str, 1900) itn, ObjTru
      call snPRNT( 3, str, iw, leniw )

      if ( nlnInf ) then
         write(str, 1910) nInf, sInf
         call snPRNT( 3, str, iw, leniw )
         if (.not. linInf) then
            write(str, 1915) wtInf, fMrt/wtInf
            call snPRNT( 3, str, iw, leniw )
         end if
      end if

      if (GotFun) then
         if (nonlin) then
            write(str, 1920) nMajor, flin
            call snPRNT( 3, str, iw, leniw )
            write(str, 1930) PenNrm, fObj
            call snPRNT( 3, str, iw, leniw )
            write(str, 1950) iw(nfObj1), iw(nfCon1)
            call snPRNT( 3, str, iw, leniw )
         end if
         if (iw(lvlDer) .lt. 3  .or.
     &        (nonlin  .and.  lvlSch .eq. 0)) then
            write(str, 1955) iw(nfObj2), iw(nfCon2)
            call snPRNT( 3, str, iw, leniw )
         end if
         if (iw(lvlDer) .lt. 3) then
            write(str, 1960) iw(nfObj3), iw(nfCon3)
            call snPRNT( 3, str, iw, leniw )
            write(str, 1962) iw(nfObj4), iw(nfCon4)
            call snPRNT( 3, str, iw, leniw )
         end if
         if (nS .gt. 0) then
            write(str, 1970) nS, nnb
            call snPRNT( 3, str, iw, leniw )
         end if
         if (iw(cgItns) .gt. 0) then
            write(str, 1973) iw(cgItns)
            call snPRNT( 3, str, iw, leniw )
         end if
         write(str, 1975) nDegen, degen
         call snPRNT( 3, str, iw, leniw )
      end if

*     ------------------------------------------------------------------
*     Unscale, compute nonlinear constraint violations,
*     save basis files and prepare to print the solution.
*     Clock 3 is "Output time".
*     ------------------------------------------------------------------
      call s1time( 3, 0, iw, leniw, rw, lenrw )

*     Skip the functions if we don't have them.
*     Skip unscaling everything for infeasible linear constraints, they
*     have already been unscaled.

      lsSave  = iw(lvlScl)

      if (.not. GotFun  .or.  FPonly) then
         nnCon1 = 0
         nnObj1 = 0
         if (linInf)
     &      iw(lvlScl) = 0
      else
         nnCon1 = nnCon
         nnObj1 = nnObj
      end if

      call s4savB
     &   ( iExit, SaveB, minimz, m, n, nb, nkx,
     &     nnCon0, nnCon1, nnL0, nnObj1, nName, nS,
     &     itn, nInf, sInf, wtInf, vimax, iObj, sclObj, ObjTru,
     &     pNorm1, pNorm2, piNorm, xNorm,
     &     ne, nlocJ, locJ, indJ, Jcol, iw(lkx),
     &     iw(lhEsta), hs, rw(lAscal), bl, bu, rw(lFx), rw(lgObj),
     &     Names, pi, rc, x,
     &     cw, lencw, iw, leniw, rw, lenrw )

*     If task = 'Print', s4savB prints the solution under the control
*     of lprSol (set by the  Solution  keyword in the SPECS file).
*     The printed solution may or may not be wanted, as follows:
*
*     lprSol = 0   means      No
*            = 2   means      Yes

      call s4savB
     &   ( iExit, PrintS, minimz, m, n, nb, nkx,
     &     nnCon0, nnCon1, nnL0, nnObj, nName, nS,
     &     itn, nInf, sInf, wtInf, vimax, iObj, sclObj, ObjTru,
     &     pNorm1, pNorm2, piNorm, xNorm,
     &     ne, nlocJ, locJ, indJ, Jcol, iw(lkx),
     &     iw(lhEsta), hs, rw(lAscal), bl, bu, rw(lFx), rw(lgObj),
     &     Names, pi, rc, x,
     &     cw, lencw, iw, leniw, rw, lenrw )
      iw(lvlScl) = lsSave
      call s1time(-3, 0, iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     If the user hasn't already pulled the plug,
*     call the functions one last time with  Status .ge. 2.
*     Everything has been  unscaled, so we have to disable scaling.
*     modefg = 0  tells the functions that gradients are not required.
*     ------------------------------------------------------------------
      if (.not. linInf  .and.  iExit/10 .ne. 7) then
         lsSave = iw(lvlScl)
         iw(lvlScl) = 0
         Status = 2 + min( iExit/10,4 )
         modefg = 0
         call fgwrap
     &      ( inform, modefg, Status, nlnCon, nlnObj,
     &        n, negCon, nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj,
     &        fgcon, fgobj, x,
     &        ne, nlocJ, locJ, indJ,
     &        rw(lfCon), fObj, rw(lgCon2), rw(lgObj2),
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )
         iw(lvlScl) = lsSave
         if (modefg .ge. 0) then
            Status = 0
         end if
      end if

*     Save some things needed by solvers calling SNOPT

  999 rw(421) = ObjTru          ! The true objective
      rw(422) = piNorm          ! Lagrange multiplier norm
      rw(423) = xNorm           ! Norm of the variables (for GAMS)
      rw(424) = wtInf           ! Infeasibility weight

      rw(433) = sInf            ! Sum of infeasibilities
      rw(434) = fLin            ! Linear objective
      rw(435) = fObj            ! Objective function
      rw(436) = PenNrm          ! Norm of penalty parameters

      iw(421) = itn             ! Total iteration count
      iw(422) = nMajor          ! Major iterations
      iw(423) = maxS            ! max # of superbasics

      return

 1900 format(
     &     ' No. of iterations', i20, 2x,
     &     ' Objective value', 1p, e22.10)
 1910 format(
     &     ' No. of infeasibilities', i15, 2x,
     &     ' Sum of infeas', 1p, e24.10)
 1915 format(
     &     ' Elastic weight            ', 1p, e11.1, 2x,
     &     ' Scaled Merit ', 1p, e24.10)
 1920 format(
     &     ' No. of major iterations', i14, 2x,
     &     ' Linear objective', 1p, e21.10)
 1930 format(
     &     ' Penalty parameter', 1p, e20.3, 2x,
     &     ' Nonlinear objective', 1p, e18.10)
 1950 format(
     &     ' No. of calls to funobj', i15, 2x,
     &     ' No. of calls to funcon', i15)
 1955 format(    ' Calls with modes 1,2 (known g)', i7,
     &        2x,' Calls with modes 1,2 (known g)', i7)
 1960 format(
     &     ' Calls for forward differencing', i7, 2x,
     &     ' Calls for forward differencing', i7)
 1962 format(
     &     ' Calls for central differencing', i7, 2x,
     &     ' Calls for central differencing', i7)
 1970 format(
     &     ' No. of superbasics', i19, 2x,
     &     ' No. of basic nonlinears', i14)
 1973 format(
     &     ' No. of CG iterations', i17)
 1975 format(
     &     ' No. of degenerate steps', i14, 2x,
     &     ' Percentage', f27.2)

      end ! subroutine s8solv

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8SQP
     &   ( iExit, fgwrap, fgcon, fgobj,
     &     MjrLog, MnrLog, snSTOP, gotR,
     &     itn, lenR, m, maxS, mBS, n, nb, nS,
     &     nnCon0, nnCon, nnObj0, nnObj, nnL0, nnL,
     &     nMajor, nMinor, nDegen, duInf,
     &     minimz, iObj, sclObj, ObjAdd, fObj, fMrt,
     &     vimax, virel, viSup, nInf, sInf,
     &     wtInf0, wtInf, PenNrm, piNorm, xNorm,
     &     ne , nlocJ, locJ, indJ, Jcol, negCon, nlocG, locG,
     &     hElast, hEstat, hfeas, hs, kBS,
     &     Ascale, bl, bu, blBS, buBS, Fv, Fx,
     &     fCon , gCon , gObj ,
     &     fCon1, gCon1, gObj1,
     &     fCon2, gCon2, gObj2,
     &     gBS, gQP, dLmul, dx, dg, Udx, Hdx,
     &     Lmul, Lmul1, Lmul2, pi, QPrhs,
     &     R, rc, rg, x, x1, xBS,
     &     xQP0, xQP, xPen,
     &     iy, iy1, y, y1, y2, y3,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     fgwrap, fgcon, fgobj, MjrLog, MnrLog, snSTOP
      logical
     &     gotR
      integer
     &     iExit, iObj, itn, lencu, lencw, leniu, leniw, lenru, lenrw,
     &     lenR, maxS, mBS, m, minimz, n, nb, nDegen, ne, negCon, nInf,
     &     nlocG, nlocJ, nMajor, nMinor, nnCon0, nnCon, nnL0, nnL,
     &     nnObj0, nnObj, nS, locJ(nlocJ), indJ(ne), hElast(nb), hs(nb),
     &     hEstat(nb), hfeas(mBS), locG(nlocG), kBS(mBS), iy(nb),
     &     iy1(nb), iu(leniu), iw(leniw)
      double precision
     &     duInf, ObjAdd, fMrt, fObj, vimax, virel, viSup, sclObj, sInf,
     &     wtInf0, wtInf, PenNrm, piNorm,
     &     Ascale(nb), bl(nb), bu(nb), blBS(mBS), buBS(mBS), dg(nnL0),
     &     dx(nb), dLmul(nnCon0), Fv(nnCon0), Fx(nnCon0), gBS(mBS),
     &     gQP(nnL0), Hdx(nnL0), Jcol(ne),
     &     fCon(nnCon0) , gCon(negCon) , gObj(nnObj0) ,
     &     fCon1(nnCon0), gCon1(negCon), gObj1(nnObj0),
     &     fCon2(nnCon0), gCon2(negCon), gObj2(nnObj0),
     &     Lmul(nnCon0), Lmul1(nnCon0), Lmul2(nnCon0),
     &     rc(nb), rg(maxS), x(nb), x1(nb), xBS(mBS), xQP(nb), xQP0(nb),
     &     xPen(nnCon0), pi(m), QPrhs(m), R(lenR), Udx(nnL0),
     &     y(nb), y1(nb), y2(nb), y3(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     s8SQP  solves a nonlinear programming problem.
*     A basis is assumed to be specified by nS, hs, x and the
*     superbasic parts of kBS.
*     In particular, there must be nS values hs(j) = 2, and the
*     corresponding j's must be listed in kBS(m+1) thru kBS(m+ns).
*     The ordering in kBS(m+1:m+nS) matches the reduced Hessian R.
*
*     On entry, if there are nonlinear constraints, Fx contains
*     the true nonlinear slacks (i.e., constraint values)
*     Fx  =  fCon + (linear A)*x,   excluding slacks.
*
*     On exit, if  iExit .lt. 30  it is safe to save the final
*     basis files and print the solution.  Otherwise, a fatal error
*     condition exists and numerous items will be undefined.
*     The last basis map saved (if any) retains the only useful
*     information.
*
*     30 Dec 1991: First version based on npsol routine npcore.
*     23 Oct 1993: Proximal point FP added.
*     29 Oct 1993: Crash on LG rows moved outside s5QP.
*     24 Apr 1994: Nx columns no longer in Q.
*     26 May 1995: Column order of R defined by kBS.
*     04 Aug 1995: Limited memory update
*     11 Aug 1995: tolg changed from 0.1 to 1.0d-4.
*     09 Nov 1995: Updated multipliers used to define Lagrangian.
*     19 Dec 1995: Finite-differences added.
*     09 Oct 1996: First Min Sum version.
*     16 Jul 1997: First thread-safe version.
*     09 Jul 1998: Quasi-Newton updates implemented correctly.
*     24 Aug 1998: Fixed bug in s8x1 found by Alan Brown at Nag.
*     06 Sep 1998: Pre- and post-QP diagonal Hessian scaling added.
*     16 Jan 1999: Name changed from s8core.
*     06 Apr 2001: For Hot Starts, don't mess with Utol.
*     27 Apr 2001: wtMax introduced as parameter to s8wInf.
*     15 Jan 2003: CG and QN  QP solvers added.
*     03 Aug 2003: snEXIT and snPRNT adopted.
*     04 Jul 2005: Switched to vanilla CG for QN with nS > maxR.
*     04 Jul 2005: Current version of s8SQP.
*     ==================================================================
      character
     &     str*132
      external
     &     ddot, dnormi, dnrm1s, dnrm2
      logical
     &     KTcond(2), backTr, badLS, boostd, centrl, debug, duFeas,
     &     Elastc, FDObj , FDCon, feasbl, feaSlk, FPonly, frstQP, goodG,
     &     incRun, maxIts, maxnS, nearOpt, needLU, newB, newG, newLU,
     &     NewTol, newx, nlnCon, nlnObj, optiml, prFeas, prtLog,
     &     prtSum, QPpi0, rowFea, restrt, tnyStp, useFD, usefLS
      integer
     &     cdItns, Hcalls, Htype, iAbort, iMsg, info(6), inform,
     &     iPrint, iSumm, itQP, j, jObj, jprInf, jduInf, klog,
     &     kSumm, lprSch, LUitn, LUreq, lvlDer, lvlDif, lvlHes,
     &     lvlPiv, lvlPre, lvlSch, lvlSrt, maxR, maxvi, minmax, mMajor,
     &     modefg, mStart, MjrPrt, MnrPrt, MjrHdP, MjrHdS, MnrHdP,
     &     nInfQP, nnJac, nSwap, nSkip, nStart, PreCon, QPmode, QPslvr,
     &     RtRmod, Status, typeLU, Utol1, Utol2
      double precision
     &     back, condHz, ddot, dnormi, dnrm1s, dnrm2, dxHdx, eps,
     &     eps0, eps1, eps5, fObj1, fObj2, fObjQP, fMrt1, gMrt,
     &     gMrt1, gNorm0, gNorm, U0ii, dRzmax, dRzmin, rviol, PenDmp,
     &     PenMax, prInf, pHpMrt, rnnL, sInfQP, sInf1, sInf2, sgnObj,
     &     step, steplm, stepmn, stepmx, tolFP, tolQP, tolQPk, tolNLP,
     &     tolCon, tolx, Utol1s, Utol2s, Utolmn, weight, wolfeG, wtMax,
     &     wtFac, wtScal, xdNorm, xNorm, xPen0
*     ------------------------------------------------------------------
      integer            QPChol,     CG,     QN
      parameter         (QPChol = 0, CG = 1, QN = 2)
      integer            HOT
      parameter         (HOT    = 3)
      integer            HUnset,     HUnit
      parameter         (HUnset =-1, HUnit  = 2)
      integer            iQNtyp,     iStep,      iQPerr,     iFDiff
      parameter         (iQNtyp = 1, iStep  = 3, iQPerr = 5, iFDiff = 6)
      integer            NO
      parameter         (NO     = 0)
      integer            Normal
      parameter         (Normal = 0)
      integer            SetWt,      IncWt
      parameter         (SetWt  = 0, IncWt  = 1)
      integer            LM   ,      FM
      parameter         (LM     = 0, FM     = 1)
      integer            B,          BS        , BT
      parameter         (B      = 0, BS     = 2, BT     = 3)
      integer            MinTol,     RstTol
      parameter         (MinTol = 2, RstTol = 3)
      integer            UnLim ,     VioLim,     UsrLim
      parameter         (UnLim  = 0, VioLim = 1, UsrLim = 2)

      double precision   zero,          half,           one
      parameter         (zero =0.0d+0,  half   =0.5d+0, one =  1.0d+0)
      double precision   ten
      parameter         (ten  =10.0d+0)
      double precision   U0max,         U0min
      parameter         (U0max= 1.0d+1, U0min  =  1.0d-2)

      parameter         (Utol1     = 154) ! abs tol for small diag of U.
      parameter         (Utol2     = 155) ! rel tol for small diag of U.
      parameter         (lvlDif    = 181) ! =1(2) forwd (cntrl) diffs
      parameter         (Hcalls    = 188) ! number of Hx products
      parameter         (Htype     = 202) ! Approximate Hessian type
      parameter         (QPmode    = 208) ! Current QP solver
      parameter         (PreCon    = 209) ! Current precon mode
      parameter         (LUitn     = 215) ! itns since last factorize
      parameter         (MnrHdP    = 223) ! >0 => Mnr heading for iPrint
      parameter         (MjrHdP    = 224) ! >0 => Mjr heading for iPrint
      parameter         (MjrHdS    = 226) ! >0 => Mjr heading for iSumm
*     ------------------------------------------------------------------
      character          line*4
      character          msg(4:8)*19
      data               line/'----'/
      data               msg /'max step too small.',
     &                        'step too small.    ',
     &                        'no minimizer.      ',
     &                        'too many functions.',
     &                        'uphill direction.  '/
*     ------------------------------------------------------------------
      iPrint    = iw( 12) ! Print file
      iSumm     = iw( 13) ! Summary file

      nnJac     = iw( 21) ! # nonlinear Jacobian variables
      maxR      = iw( 52) ! max columns of R.
      QPslvr    = iw( 55) ! = 0:1:2   => QPChol:CG:QN QP solver
      klog      = iw( 61) ! log/print frequency
      kSumm     = iw( 62) ! Summary print frequency
      lvlSrt    = iw( 69) ! = 0:1:2:3 => cold:warm:basis:hot start
      lvlHes    = iw( 72) ! 0,1,2  => LM, FM, Exact Hessian
      lvlSch    = iw( 76) ! >0     => use derivatives in the line search
      lvlPre    = iw( 77) ! >0     => QN preconditioned CG
      lvlPiv    = iw( 80) ! 0/1 Threshold partial/complete pivoting: user
      lprSch    = iw( 82) ! line search debug starting itn

      minmax    = iw( 87) ! 1, 0, -1  => MIN, FP, MAX
      mMajor    = iw( 90) ! limit on major iterations
      MjrPrt    = iw( 92) ! Major print level
      MnrPrt    = iw( 93) ! Minor print level
      lvlDer    = iw( 70) ! = 0, 1, 2 or 3, the derivative level

*     Constants

      eps       = rw(  1) ! machine precision.  IEEE DP  2.22e-16
      eps0      = rw(  2) ! eps**(4/5)          IEEE DP  3.00e-13
      eps1      = rw(  3) ! eps**(2/3)          IEEE DP  3.67e-11
      eps5      = rw(  7) ! eps**(1/5)          IEEE DP  7.40e-04

      tolFP     = rw( 51) ! Minor Phase 1 Opt tol
      tolQP     = rw( 52) ! Minor Phase 2 Opt tol
      tolNLP    = rw( 53) ! Major Optimality tolerance
      tolx      = rw( 56) ! Minor feasibility tolerance.
      tolCon    = rw( 57) ! Major feasibility tolerance.
      wolfeG    = rw( 84) ! line search tolerance.
      xPen0     = rw( 89) ! initial penalty parameter.
      wtMax     = rw( 90) ! max     infeasibility weight

      nlnCon    = nnCon  .gt. 0
      nlnObj    = nnObj  .gt. 0
      FPonly    = minmax .eq. 0

      iw(Hcalls) = 0
      iw(MnrHdP) = 0
      iw(MjrHdP) = 0
      iw(MjrHdS) = 0

      iw(QPmode) = QPslvr ! Current QP solver
      iw(PreCon) = lvlPre ! Current precon mode

*     ------------------------------------------------------------------
*     s8SQP  operates in either ``Normal'' or ``Elastic'' mode.
*     In elastic mode, the nonlinear slacks are allowed to be infeasible
*     while a weighted sum of the slack infeasibilities is minimized.
*     ------------------------------------------------------------------
      feaSlk = .true.
*     Elastc =       FP  .and.  nlnCon
      Elastc = .false.
      call iload ( nb, 0, hEstat, 1 )

      nInf   = 0
      sInf   = zero
      sInf1  = zero
      sInf2  = zero

      iExit  = 0
      LUreq  = 0
      Status = 0
      nSkip  = 0
      nStart = 0
      if (nnL .gt. 0) then
         mStart = 2
      else
         mStart = 0
      end if
      RtRmod = 0

      call iload ( 6, 0, info, 1 )
      info(iQNtyp) = 1          ! Suppresses first printing of 'n'

      sgnObj = minimz
*
      gNorm  = one
      if (iObj   .gt. 0) gNorm  = gNorm + sclObj
      gNorm0 = zero
      if (nlnObj) then
*        gNorm0 = dnormi( nnObj, gObj, 1 )
         gNorm0 = dnrm2 ( nnObj, gObj, 1 )
      end if

*     Hwt    = 0.001d+0         ! Weight on the Hessian
*     if (gNorm0 .gt. zero) then
*     U0ii   = sgnObj*Hwt*gNorm0
*      else
*     U0ii   = sgnObj
*     end if

      if (nnL .eq. 0) then
         U0ii = one
      else
         if (gNorm0 .gt. zero) then
            rnnL = nnL
            U0ii = sqrt(gNorm0/sqrt(rnnL))
         else
            U0ii = one
         end if
      end if
      U0ii   = min( max( U0ii, U0min ), U0max )

      rviol  = zero
      prInf  = zero
      duInf  = zero
      wtInf  = wtInf0
      tolQPk = 10.0d+0 * tolQP

      gMrt   = zero
      step   = zero

      KTcond(1) =  .false.
      KTcond(2) =  .false.

      frstQP = .true.
      QPpi0  = .false.      ! Use zero initial multipliers
*     QPpi0  = .true.       ! Use QP initial multipliers

      condHz = one

      FDObj  = (lvlDer .eq. 0  .or.  lvlDer .eq. 2) .and. (nnObj .gt. 0)
      FDCon  = (lvlDer .eq. 0  .or.  lvlDer .eq. 1) .and. (nnJac .gt. 0)
      useFD  =  FDObj  .or.  FDCon
      usefLS =  useFD          .or.  lvlSch .eq. 0

      if (MjrPrt .ge. 10  .or.  MnrPrt .ge. 10) then
         prtLog = iPrint .gt. 0  .and.  klog  .eq. 1
         prtSum = iSumm  .gt. 0  .and.  kSumm .eq. 1
         if (prtLog) then
            write(str, 1000) (line, j=1,29)
            call snPRNT( 1, str, iw, leniw )
            write(str, 1010) nMajor
            call snPRNT( 1, str, iw, leniw )
         end if
         if (prtSum  .and.  MnrPrt .ge. 10) then
            write(str, 1000) (line, j=1,19)
            call snPRNT( 2, str, iw, leniw )
            write(str, 1010) nMajor
            call snPRNT( 2, str, iw, leniw )
         end if
      end if

      jObj   = n + iObj

      if ( nlnCon ) then
*        ---------------------------------------------
*        Initialize the penalty parameters.
*        Set an initial elastic weight.
*        ---------------------------------------------
         incRun = .true.
         PenDmp = one
         PenMax = one / eps
         PenNrm = xPen0
         call dload ( nnCon, xPen0, xPen, 1 )
      end if

      if (nnL .gt. 0  .and.  iw(Htype) .eq. HUnset) then
*        ---------------------------------------------------------------
*        The approximate Hessian needs to be initialized.
*        Use the identity matrix until something better comes along.
*        ---------------------------------------------------------------
         call s8H0
     &      ( iw(Htype), nnL, U0ii, iw, leniw, rw, lenrw )
      end if

      if (nS .gt. maxR) then
         iw(QPmode)   = CG      ! Use CG
         iw(PreCon)   = NO      ! with no preconditioning
         gotR         = .false.
      end if

      call dcopy ( nb, x, 1, xQP, 1 )
      cdItns = -1
      newG   = .false.

**    ======================Start of main loop==========================
*     Start of a Major Iteration.
*     ==================================================================
*+    do while (iExit .eq. 0)
  100 if       (iExit .eq. 0) then

         nMinor = 0

*        ===============================================================
*        Repeat                    (until an accurate gradient is found)

  110       centrl = iw(lvlDif) .eq. 2

            if ( newG ) then
               if ( useFD ) then
*                 ------------------------------------------------------
*                 Compute any missing derivatives.
*                 ------------------------------------------------------
                  call s6fdG
     &               ( iExit, n, negCon,
     &                 nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj,
     &                 fgwrap, fgcon, fgobj,
     &                 bl, bu, x,
     &                 ne, nlocJ, locJ, indJ,
     &                 fCon, fObj, gCon, gObj, y,
     &                 cu, lencu, iu, leniu, ru, lenru,
     &                 cw, lencw, iw, leniw, rw, lenrw )
                  if (iExit .ne. 0) go to 100 ! Break
               end if ! useFD
               newG = .false.
            end if

            if ( nlnCon ) then
*              ---------------------------------------------------------
*              Load the scaled Jacobian in J.
*              Compute the QP right-hand side   QPrhs  =  Jx - fCon.
*              Find Fx the nonlinear constraint values.
*              ---------------------------------------------------------
               call s8Gcpy
     &            ( nnCon, nnJac, ne, nlocJ, locJ, indJ,
     &              negCon, nlocG, locG, gCon,
     &              ne, nlocJ, locJ, Jcol )
               call dcopy
     &            ( nnCon, fCon, 1, QPrhs, 1 )
               call s2Aprd
     &            ( Normal, eps0,
     &              ne, nlocJ, locJ, indJ, Jcol,
     &              one, x, nnJac, (-one), QPrhs, nnCon )
*              ---------------------------------------------------------
*              s8sOpt  finds the nonlinear slacks  sN  that minimize the
*              merit function with  x(1:n)  and  Lmul  held fixed.
*              The optimal slacks are loaded into  x(n+1:nb)  and the
*              violations are calculated:
*                    Fv = fCon  + A(linear)x - nonlinear slacks
*                       = Fx                 - sN
*              ---------------------------------------------------------
               gNorm  = one
               if (iObj   .gt. 0) gNorm  = gNorm + sclObj
               gNorm0 = zero
               if (nlnObj       ) gNorm0 = dnrm1s( nnObj, gObj, 1 )

               if (.not. Elastc) then
                  call s8wInf
     &               ( SetWt, boostd, itn, (gNorm+gNorm0),
     &                 wtInf0, wtInf, wtMax,
     &                 weight, wtFac, wtScal, iw, leniw )
               end if

               call s8sOpt
     &            ( n, nb, nnCon, piNorm, eps0, wtInf, hEstat,
     &              bl, bu, Fv, x, Lmul, xPen, Fx )
*              call s8Fv
*    &            ( Elastc, n, nnCon, eps0, wtInf,
*    &              bl, bu, Fv, x, Lmul, Fx )

            end if

*           ------------------------------------------------------------
*           Prepare to (re-)solve the QP subproblem (possibly after the
*           elastic weight has been increased).
*           ------------------------------------------------------------
*           Factorize the basis at x.
*           Compute xQP such that (J -I)*xQP = rhs.

  300       if ( frstQP ) then
*              ---------------------------------------------------------
*              First QP subproblem.
*              ---------------------------------------------------------
               needLU = .true.
               gotR   = .false.
               nSwap  = 0
               if (nS .eq. 0) then
                  typeLU = B
               else
                  typeLU = BS
               end if

               Utol1s = rw(Utol1)
               Utol2s = rw(Utol2)

*              To avoid an unnecessarily ill-conditioned starting basis
*              for the first QP, use big singularity tols
*              (except if it's a Hot Start!).

               if (lvlSrt .eq. HOT) then
                  Utolmn = eps1
               else
                  Utolmn = eps5
               end if

               rw(Utol1) = max( Utol1s, Utolmn )
               rw(Utol2) = max( Utol2s, Utolmn )

            else
*              ---------------------------------------------------------
*              Subsequent factorizations.
*              ---------------------------------------------------------
*              For linearly constrained problems, the factors L, U and R
*              can be saved as long as a poor x does not force a
*              new factorization. (Even in this case, R can be saved if
*              there are no swaps.)

               needLU = nlnCon
               typeLU = BT

*              Reset the factor and update tols if they were changed
*              during the previous major iteration.

               call s2tols
     &            ( RstTol, NewTol, itn, iw, leniw, rw, lenrw )
            end if

            call s2Bfac
     &         ( iExit, typeLU, needLU, newLU, newB,
     &           iObj, itn, MjrPrt, LUreq,
     &           m, mBS, n, nb, nnL, nS, nSwap,
     &           ne, nlocJ, locJ, indJ, Jcol,
     &           kBS, hs, bl, bu, blBS, buBS,
     &           nnCon0, nnCon, QPrhs, xQP, xBS,
     &           iy, iy1, y, y2, iw, leniw, rw, lenrw )
            if (iExit .ne. 0) go to 100 ! Break

            if (iw(QPmode) .eq. QPChol) then
               gotR = gotR  .and.  .not. newLU
            end if

            needLU  = .false.
            if (MjrPrt .ge. 10) iw(MjrHdP) = 1

            if ( frstQP ) then
               iw(80)    = lvlPiv ! Reset original TPP or TCP
               rw(Utol1) = Utol1s
               rw(Utol2) = Utol2s
            end if

*           ------------------------------------------------------------
*           Solve the QP subproblem to obtain kBS, xQP and pi.
*           The search direction will be dx = xQP - x.
*           Use x1 to store the first feasible point.
*           ------------------------------------------------------------
            if (nS .le. maxR) then
               if (QPslvr .eq. QPChol  .or.  QPslvr .eq. QN) then
                  iw(QPmode) = QPslvr
               end if
            end if

            inform = 0

            if (     iw(QPmode) .eq. QPChol) then
               call s8iQP
     &            ( inform, info, iw(Htype), Mnrlog, iw(Hcalls), Elastc,
     &              gotR, itn, itQP, lenR, m, maxR, mBS, n, nb,
     &              nnCon0, nnCon, nnObj0, nnObj, nnL0, nnL, nS, nDegen,
     &              MjrPrt, MnrPrt, minimz, iObj, sclObj, (ObjAdd+fObj),
     &              fObjQP, tolFP, tolQPk, tolx, nInfQP, sInfQP, wtInf,
     &              U0ii, piNorm, ne, nlocJ, locJ, indJ, Jcol,
     &              hElast, hEstat, hfeas, hs, kBS,
     &              Ascale, bl, bu, blBS, buBS, gBS, gQP, gObj, Hdx,
     &              y3, pi, R, rc, rg, QPrhs, x,
     &              xQP, xBS, xQP0, iy, iy1, y, y1, y2,
     &              cu, lencu, iu, leniu, ru, lenru,
     &              cw, lencw, iw, leniw, rw, lenrw )
            else if (iw(QPmode) .eq. QN    ) then
               call s8iQN
     &            ( inform, info, iw(Htype), Mnrlog, iw(Hcalls), Elastc,
     &              gotR, itn, itQP, lenR, m, maxR, mBS, n, nb,
     &              nnCon0, nnCon, nnObj0, nnObj, nnL0, nnL, nS, nDegen,
     &              MjrPrt, MnrPrt, minimz, iObj,
     &              condHz, sclObj, (ObjAdd+fObj), fObjQP,
     &              tolFP, tolQPk, tolx, nInfQP, sInfQP, wtInf,
     &              U0ii, piNorm, ne, nlocJ, locJ, indJ, Jcol,
     &              hElast, hEstat, hfeas, hs, kBS,
     &              Ascale, bl, bu, blBS, buBS, gBS, gQP, gObj, Hdx,
     &              y3, pi, R, rc, rg, y1(m+1), QPrhs, x,
     &              xQP, xBS, xQP0, iy, iy1, y, y1, y2,
     &              cu, lencu, iu, leniu, ru, lenru,
     &              cw, lencw, iw, leniw, rw, lenrw )
            end if

            if (inform .eq. -2  .and.  maxR .lt. maxS) then
               iw(QPmode)   = CG ! Switch to CG
               iw(PreCon)   = NO ! with no preconditioning
               gotR         = .false.
               info(iQPerr) = 0
            end if

            if (iw(QPmode) .eq. CG) then
               call s8iQN
     &            ( inform, info, iw(Htype), Mnrlog, iw(Hcalls), Elastc,
     &              gotR, itn, itQP, lenR, m, maxS, mBS, n, nb,
     &              nnCon0, nnCon, nnObj0, nnObj, nnL0, nnL, nS, nDegen,
     &              MjrPrt, MnrPrt, minimz, iObj,
     &              condHz, sclObj, (ObjAdd+fObj), fObjQP,
     &              tolFP, tolQPk, tolx, nInfQP, sInfQP, wtInf,
     &              U0ii, piNorm, ne, nlocJ, locJ, indJ, Jcol,
     &              hElast, hEstat, hfeas, hs, kBS,
     &              Ascale, bl, bu, blBS, buBS, gBS, gQP, gObj, Hdx,
     &              y3, pi, R, rc, rg, y1(m+1), QPrhs, x,
     &              xQP, xBS, xQP0, iy, iy1, y, y1, y2,
     &              cu, lencu, iu, leniu, ru, lenru,
     &              cw, lencw, iw, leniw, rw, lenrw )
            end if

*           inform    Status
*           ------    ------
*            >0       Fatal error
*             0       QP solution found
*            -1       Too many iterations
*            -2       Too many superbasics

            nMinor = nMinor + itQP
            if (inform .gt. 0) then
               iExit = inform
               go to 100
            end if

*           Remember the value of the QP inform until after the printing

            maxnS  = inform .eq. -2
            maxIts = inform .eq. -1

            if ( frstQP ) then
               frstQP = .false.
            end if

            if ( nlnCon ) then
               if ( QPpi0 ) then
                  call dcopy ( nnCon,            pi, 1,  Lmul, 1 )
                  call dload ( nnCon, (zero), dLmul, 1 )
               else
                  call dcopy ( nnCon,            pi, 1, dLmul, 1 )
                  call daxpy ( nnCon, (-one),  Lmul, 1, dLmul, 1 )
               end if

               if (Elastc  .and.  feaSlk) then ! this only happens once
                  call dcopy ( nnCon, Fx, 1, x(n+1), 1 )
                  call dcopy ( nnCon, pi, 1,   Lmul, 1 )
                  call dload ( nnCon, (zero), dLmul, 1 )
               end if

*              If Lmul or x  changed, recompute Fv.

               if (QPpi0  .or.  (Elastc  .and.  feaSlk)) then
                  call s8sOpt
     &               ( n, nb, nnCon, piNorm, eps0, wtInf, hEstat,
     &                 bl, bu, Fv, x, Lmul, xPen, Fx )
*                 call s8Fv  ( Elastc, n, nnCon, eps0, wtInf,
*    &                 bl, bu, Fv, x, Lmul, Fx )

                  if (Elastc  .and.  feaSlk) feaSlk = .false.
                  if (QPpi0                ) QPpi0  = .false.
               end if

*              Find the sum of infeasibilities of the nonlinear slacks.

               call s8sInf
     &            ( n, nb, nnCon, tolx, nInf, sInf, bl, bu, x )
            end if

*           Compute the search direction dx.

            call dcopy ( nb,         xQP, 1, dx, 1 )
            call daxpy ( nb, (-one), x  , 1, dx, 1 )

            xNorm  = dnrm1s( n, x , 1 )
            xdNorm = dnrm1s( n, dx, 1 )

*           Compute all the QP reduced costs.
*           (We could use Lmul for the nonlinear pi's).
*           Compute the maximum dual infeasibility.

            call s8rc
     &         ( sclObj, minimz, iObj, m, n, nb,
     &           nnObj0, nnObj, nnCon, nnJac, negCon,
     &           ne, nlocJ, locJ, indJ, Jcol,
     &           gObj, gCon, pi, rc )
            call s8Infs
     &         ( Elastc, n, nb, nnCon0, nnCon, tolx, wtInf,
     &           prInf, duInf, jprInf, jduInf, bl, bu, Fx, rc, x )

*           Compute the largest nonlinear row violation.

            if ( nlnCon ) then
               rviol = dnormi( nnCon, Fv, 1 )
            end if

*           ------------------------------------------------------------
*           Test for convergence.
*           ------------------------------------------------------------
            if ( gotR ) then
               call s6Rcnd
     &            ( maxR, nS, lenR, R, dRzmax, dRzmin, condHz )
            end if
            rviol     = rviol /(one + xNorm )
            prInf     = prInf /(one + xNorm )
            duInf     = duInf /(one + piNorm)

            rowFea    = rviol  .lt. tolCon  .and.  nInf .gt. 0
            prFeas    = prInf  .le. tolCon
            duFeas    = duInf  .le. tolNLP
            KTcond(1) = prFeas
            KTcond(2) = duFeas

            feasbl    =       prFeas   .or.   rowFea
            optiml    =       duFeas  .and.  feasbl
            nearOpt   = .not. optiml  .and. (KTcond(1)
     &                                .and.  duInf .lt. ten*tolNLP
     &                                .or.   KTcond(2)
     &                                .and.  prInf .lt. ten*tolCon)
            if (nlnCon  .and.  optiml  .and.  nInf .gt. 0) then
               call s8wInf
     &            ( IncWt, boostd, itn, (gNorm+gNorm0),
     &              wtInf0, wtInf, wtMax,
     &              weight, wtFac, wtScal, iw, leniw )

               if ( boostd ) then
                  Elastc = .true.
                  go to 300     ! Solve the QP again
               end if
            end if

*            if ( Elastc ) then
**              tolQPk = tolQP
*               tolQpk = min(tolQPk, 0.1d+0/wtInf)
*            else
               if (nMajor .eq. 0) then
                  tolQPk = min( duInf, 1.0d-3 )
               end if
               if (nMinor .eq. 0  .and.  tolQPk .gt. tolQP) then
                  tolQPk = 0.2d+0*tolQPk
               end if
               tolQPk = min( 0.5d+0*tolQPk, 0.1d+0*duInf )
               tolQPk = max( tolQPk, tolQP )
*            end if

*           ------------------------------------------------------------
*           Compute the current augmented Lagrangian merit function.
*           ObjAdd is added in the log routine.
*           ------------------------------------------------------------
            if (iObj .eq. 0) then
               fMrt = zero
            else
               fMrt = sgnObj*x(jObj)*sclObj
            end if

            if ( nlnObj ) then
               fMrt =  fMrt + sgnObj*fObj
            end if

            if ( nlnCon ) then
               call dcopy ( nnCon, Fv  , 1, y, 1 )
               call ddscl ( nnCon, xPen, 1, y, 1 )
               fMrt = fMrt -      ddot  ( nnCon, Lmul, 1, Fv, 1 )
     &                     + half*ddot  ( nnCon,    y, 1, Fv, 1 )

               if ( Elastc ) then
                  fMrt = fMrt + wtInf*sInf
               end if
            end if

*           ------------------------------------------------------------
*           If the forward-difference estimate of the reduced gradient
*           of the Lagrangian is small,  prepare to: (i) switch to
*           central differences; (ii)  recompute the derivatives,  and
*           (iii) solve the QP again.
*
*           On the other hand, if central differences give a large
*           reduced-gradient norm, switch back to forward differences.
*           ------------------------------------------------------------
            call s8FD
     &         ( nnCon0, nnCon, nnObj, itn, cdItns,
     &           centrl, goodG, newG, useFD, info, duInf,
     &           fCon, fObj, iw, leniw, rw, lenrw )

*           ------------------------------------------------------------
*           Print the details of this iteration.
*           ------------------------------------------------------------
            call MjrLog
     &         ( iAbort, info, iw(Htype), KTcond, MjrPrt,
     &           minimz, n, nb, nnCon0, nS, itn, nMajor, nMinor, nSwap,
     &           condHz, iObj, sclObj, ObjAdd, fMrt, PenNrm, step,
     &           prInf, duInf, vimax, virel, hs,
     &           ne, nlocJ, locJ, indJ, Jcol,
     &           Ascale, bl, bu, fCon, Lmul, x,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )

            call snSTOP
     &         ( iAbort, info, Htype, KTcond, MjrPrt, minimz,
     &           n, nb, nnCon0, nS, itn, nMajor, nMinor, nSwap,
     &           condHz, iObj, sclObj, ObjAdd, fMrt, PenNrm, step,
     &           prInf, duInf, vimax, virel, hs,
     &           ne, nlocJ, locJ, indJ, Jcol,
     &           Ascale, bl, bu, fCon, Lmul, x,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
            if (iAbort .ne. 0) then
               iExit = 74       ! User has aborted the run via snSTOP
               go to 100
            end if

            info(iStep ) = UnLim

*+       until (.not. (useFD  .and.  .not.goodg))
         if (          useFD  .and.  .not.goodg ) go to 110
*        ===============================================================

         if (prFeas   .and.  info(iQPerr) .eq. 3      ) iExit = 21
         if (maxIts                                   ) iExit = 31
         if (                      nMajor .ge. mMajor ) iExit = 32
         if (maxnS    .and.        nMajor .ge. mMajor ) iExit = 33
         if (nearOpt  .and.         iExit .ne. 0      ) iExit =  3
         if (optiml) then
            if (nInf .eq. 0) then ! Optimal and feasible
               if (FPonly) then
                  iExit = 2
               else
                  iExit = 1
               end if
            else                  ! Violations minimized
               iExit = 13
            end if
         end if
         if (iExit .ne. 0) go to 100

         step  = zero
         nSwap = 0

*        ===============================================================
*        Take a step in the right direction.
*        ===============================================================
*        Compute  dxHdx = s'Hs  and other directional derivatives.
*        Be prepared to fix up pHpMrt if there are linear variables.
*        ---------------------------------------------------------------
         if (nnL .gt. 0) then
            call s8xHx
     &         ( nnL, dx, Udx, Hdx, dxHdx, iw, leniw, rw, lenrw )
         else
            dxHdx = zero
         end if

         if (nnL .eq. n) then
            if (dxHdx .eq. zero  .and.  iw(Htype) .ne. HUnit) then
               iw(Htype) = HUnit
               call s8H0
     &            ( iw(Htype), nnL, U0ii, iw, leniw, rw, lenrw )
               go to 100
            end if
            pHpMrt = dxHdx
         else
            pHpMrt = max( eps1*xdNorm*xdNorm, dxHdx )
         end if

*        ---------------------------------------------------------------
*        Compute the contributions to the merit function and its
*        directional derivative from the nonlinear constraints.
*        The penalty parameters  xPen(j)  are increased if the
*        directional derivative is not sufficiently negative.
*        ---------------------------------------------------------------
*        First, compute the value and directional derivative of the
*        Lagrangian with respect to x and the multipliers.

         if (iObj .eq. 0) then
            fMrt = zero
            gMrt = zero
         else
            fMrt = sgnObj*x (jObj)*sclObj
            gMrt = sgnObj*dx(jObj)*sclObj
         end if

         if ( nlnObj ) then
            fMrt = fMrt + sgnObj*fObj
            gMrt = gMrt + sgnObj*ddot  ( nnObj, gObj, 1, dx, 1 )
         end if

         if ( Elastc ) then
            fMrt = fMrt +           sInf *wtInf
            gMrt = gMrt + (sInfQP - sInf)*wtInf
         end if

*        ---------------------------------------------------------------
*        Compute the search direction for the multipliers and nonlinear
*        slacks, and the contributions to the merit function and its
*        directional derivative from the nonlinear constraints.
*        The penalty parameters  xPen(j)  are increased if the
*        directional derivative is not sufficiently negative.
*        ---------------------------------------------------------------
         if ( nlnCon ) then
            fMrt = fMrt  - ddot  ( nnCon,  Lmul, 1, Fv, 1 )
            gMrt = gMrt  + ddot  ( nnCon,  Lmul, 1, Fv, 1 )
            gMrt = gMrt  - ddot  ( nnCon, dLmul, 1, Fv, 1 )

            call s8mrt
     &         ( nnCon, fMrt, gMrt, pHpMrt, incRun,
     &           penDmp, penMax, PenNrm, Fv, xPen, y, rw, lenrw )
         end if

*        ===============================================================
*        Find  stepmn,  stepmx  and  step,  the maximum, minimum and
*        initial values for the line search step.
*        ===============================================================
         call s8step
     &      ( centrl, usefLS, nb, nnCon, nnObj, nMajor,
     &        nSkip, step, stepmn, steplm, stepmx, eps0, xdNorm, xNorm,
     &        bl, bu, x, dx, iw, leniw, rw, lenrw )

         debug  = nMajor .ge. lprSch
         back   = 0.1d+0        ! backtracking factor

*        ===============================================================
*        Prepare for the linesearch to find a better point
*           x1 = x + step*dx  and  Lmul1 = Lmul + step*dLmul.
*        where, on entry,  x1 = xQP and  Lmul1 = pi.
*
*        fCon , gCon , gObj  and Lmul  are defined at the current    x.
*        fCon1, gCon1, gObj1 and Lmul1 are defined at the new point x1.
*        fCon2, gCon2, gObj2 and Lmul2 are temporary work arrays.
*
*        s6srch returns the following values:
*
*        inform    Result
*        ------    ------
*         >0      Fatal error
*          0      redo the search.
*         -1      The search is successful and step < stpmax.
*         -2      The search is successful and step = stpmax.
*         -3      A better point was found but no sufficient decrease.
*                 Most likely, the merit function is decreasing at the
*                 boundary, but there could be too many function calls.
*         -4      stpmax < tolabs (too small to do a search).
*         -5      step   < stepmn (lsrchq only -- maybe want to switch
*                 to central differences to get a better direction).
*         -6      No useful step.
*                 The interval of uncertainty is less than 2*tolabs.
*                 The minimizer is very close to step = zero
*                 or the gradients are not sufficiently accurate.
*         -7      there were too many function calls.
*         -8      the input parameters were bad
*                 (stpmax le toltny  or  oldg ge 0).
*        ===============================================================
*        x and sInf are saved in case we have to restart the search.
*        y  is used as x2 in s6srch.

  500    call dcopy ( nb   , xQP, 1,     y, 1 )
         if ( nlnCon )
     &   call dcopy ( nnCon,  pi, 1, Lmul2, 1 )
         if ( Elastc ) sInf2 = sInfQP

         call s6srch
     &      ( inform, fgwrap, fgcon, fgobj,
     &        debug, Elastc, usefLS, prFeas, iObj, sclObj,
     &        n, nb, nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj,
     &        itn, wolfeG, sgnObj, step, stepmn, stepmx, xdNorm, xNorm,
     &        fMrt, fMrt1, gMrt, gMrt1, sInf, sInf1, sInf2, wtInf,
     &        ne, nlocJ, locJ, indJ, Jcol, negCon, nlocG, locG,
     &        fObj1, fCon1, gCon1, gObj1, fObj2, fCon2, gCon2, gObj2,
     &        dx, dLmul, x, x1, y, Lmul, Lmul1, Lmul2, xPen,
     &        y1, y2, y3, cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )
         if (inform .gt. 0) then
            iExit = inform      ! The user wants to stop
            go to 100
         end if

         restrt = nStart .lt. mStart ! permission to restart
         newx   = step   .gt. zero
         backTr = inform .eq.    0
         badLS  = inform .le.   -4
         tnyStp = inform .eq.   -4  .or.   inform .eq. -6
     &                              .or.   inform .eq. -8

         if (useFD  .and. .not. centrl) then
*           ------------------------------------------------------------
*           If the line search failed.  Switch to central differences
*           and solve the QP subproblem again.
*           ------------------------------------------------------------
            if (badLS) then
               cdItns       = 0
               info(iFDiff) = 1
               write(str, 3020) itn
               call snPRNT( 23, str, iw, leniw )
               iw(lvlDif)  = 2
               newG        = .true.
               go to 110        ! Recompute the derivatives
            end if
         end if

         if (newx  .or.  backTr) then
*           ------------------------------------------------------------
*           See if the search needs to be redone with a smaller stepmx.
*           ------------------------------------------------------------
            if (backTr) then
               info(iStep) = UsrLim    ! User rejected the search.

            else if (nlnCon) then ! Check vimax ge viSup
               call s8Fx
     &            ( n, nnCon, nnJac, eps0,
     &              ne, nlocJ, locJ, indJ, Jcol, fCon1, x1, Fx )
               call s2vmax
     &            ( n, nnCon, maxvi, vimax, bl, bu, Fx )
               virel  = vimax / (one + xNorm)
               if (vimax .gt. viSup) then
                  backTr      = .true.
                  info(iStep) = VioLim
               end if
            end if

            if (backTr) then
               stepmx = back * step
               step   = stepmx
               back   = back*back
               go to 500        ! Repeat the line search
            end if
         end if

         if (info(iStep) .eq. VioLim) then
*           ------------------------------------------------------------
*           The line search is backing away from a violation limit.
*           If we are in elastc mode or the line search died, switch to
*           elastic mode with a bigger infeasibility weight.
*           ------------------------------------------------------------
            if (nlnCon  .and. (nInf .gt. 0  .or.  tnyStp)) then
               call s8wInf
     &            ( IncWt, boostd, itn, (gNorm+gNorm0),
     &              wtInf0, wtInf, wtMax,
     &              weight, wtFac, wtScal, iw, leniw )
               if ( boostd ) then
                  Elastc = .true.
               end if
               if (tnyStp) go to 300 ! Solve the QP again
            end if
         end if

         if (badLS  .and. .not. newx) then
*           ============================================================
*           Deal with some obvious cases.
*           ============================================================
            if      (maxnS) then
               iExit = 33       ! Superbasics limit
            else if (nearOpt) then
               iExit =  3       ! Requested accuracy could not be ...
            else if (info(iStep) .eq. VioLim) then
               iExit = 22
            else if (info(iStep) .eq. UsrLim) then
               iExit = 63
            end if
            if (iExit .gt. 0) go to 100
         end if

         if (badLS) then
*           ------------------------------------------------------------
*           The line search failed to provide a sufficient decrease.
*           ------------------------------------------------------------
            if (newx) then
*              Relax. At least we got SOME decrease.

            else
*              Desperate times.
*              If possible, reset everything and solve the QP again.

               if ( restrt ) then
                  if (iw(Htype) .ne. HUnit) then
                     iw(Htype) = HUnit
                     call s8H0
     &                  ( iw(Htype), nnL, U0ii, iw, leniw, rw, lenrw )
                  end if

                  if (iw(LUitn) .gt. 0) then
*                    ---------------------------------------------------
*                    Try and fix up the basis.
*                    ---------------------------------------------------
                     call s2tols
     &                  ( MinTol, NewTol, itn, iw, leniw, rw, lenrw )
                  end if

                  if ( nlnCon ) then
                     incRun = .true.
                     PenDmp = one
                     PenMax = one / eps
                     PenNrm = xPen0
                     call dload ( nnCon, xPen0, xPen, 1 )
                  end if
                  nStart = nStart + 1
                  frstQP = .true.
                  go to 300     ! Solve the QP again
               else
*                 ------------------------------------------------------
*                 We have run out of things to try. Bummer.
*                 ------------------------------------------------------
                  iMsg = -inform
                  write(str, 1050) iMsg, msg(iMsg), nMajor, duInf
                  call snPRNT( 23, str, iw, leniw )
                  iExit = 41    ! Current point cannot be improved...
                  go to 100
               end if
            end if
         end if

*        ===============================================================
*        The new point  x1  has been computed.
*        ===============================================================
         if (step .ge. steplm  .and.  nnL .gt. 0) then
            info(iStep) = 3
         end if

         inform = 0
         centrl = iw(lvlDif) .eq. 2

*        ---------------------------------------------------------------
*        Some unknown derivatives may need to be calculated at x1.
*        ---------------------------------------------------------------
         if (usefLS  .and.  nnL .gt. 0)  then
            modefg = 1
            call fgwrap
     &         ( iExit, modefg, Status, nlnCon, nlnObj,
     &           n, negCon, nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj,
     &           fgcon, fgobj, x1,
     &           ne, nlocJ, locJ, indJ,
     &           fCon1, fObj1, gCon1, gObj1,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
            if (iExit .ne. 0) then 
! The user GRADIENT failed...we need to do something a little better here
! Treat this situation as a failed line search and perform a back-track
               print *, "A graident failed with mode 1, redoing LS"

               stepmx = back * step
               step   = stepmx
               back   = back*back
               go to 500        ! Repeat the line search
            end if

            !if (iExit .ne. 0) go to 100 ! Break

            if ( useFD ) then
               call s6fdG
     &            ( iExit, n, negCon,
     &              nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj,
     &              fgwrap, fgcon, fgobj,
     &              bl, bu, x1,
     &              ne, nlocJ, locJ, indJ,
     &              fCon1, fObj1, gCon1, gObj1, y3,
     &              cu, lencu, iu, leniu, ru, lenru,
     &              cw, lencw, iw, leniw, rw, lenrw )
               if (iExit .ne. 0) go to 100 ! Break
            end if
         end if

         inform = 0
         nMajor = nMajor + 1
         if ( centrl )
     &   cdItns = cdItns + 1

         if (MjrPrt .ge. 10  .or.  MnrPrt .ge. 10) then
            prtLog = iPrint .gt. 0  .and.  klog  .eq. 1
            prtSum = iSumm  .gt. 0  .and.  kSumm .eq. 1

            if (prtLog) then
               call s1page( 0, iw, leniw )
               write(str, 1000) (line, j=1,29)
               call snPRNT( 1, str, iw, leniw )
               write(str, 1010) nMajor
               call snPRNT( 1, str, iw, leniw )
            end if
            if (prtSum  .and.  MnrPrt .ge. 10) then
               write(str, 1000) (line, j=1,19)
               call snPRNT( 2, str, iw, leniw )
               write(str, 1010) nMajor
               call snPRNT( 2, str, iw, leniw )
            end if
         end if

*        ===============================================================
*        The problem functions have been defined at the new x.
*        ===============================================================
         if (nnL .gt. 0 .and. (lvlHes .eq. LM .or. lvlHes .eq. FM)) then
                                ! Update a QN approximate Hessian.
            call s8HQN
     &         ( inform, fgwrap, fgcon, fgobj,
     &           useFD, iw(Htype), iw(QPmode), info,
     &           lenR, m, mBS, n, nb,
     &           nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj,
     &           nS, nMajor, nSkip, U0ii,
     &           step, minimz, dxHdx,
     &           RtRmod, gotR, incRun, PenDmp, PenMax,
     &           fObj, fCon, gCon, gObj, fCon1, gCon1, gObj1,
     &           ne, nlocJ, locJ, indJ, Jcol, negCon, nlocG, locG,
     &           kBS, bl, bu, dx, dg, Udx, Hdx, Lmul1,
     &           R, x, x1, xQP0, xPen, y, y1, y2, y3,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
            if (inform .ne. 0) then
               iExit = inform
               go to 100
            end if
         end if

*        ---------------------------------------------------------------
*        Update the variables.
*        The QP solution, saved in xQP, is used to start the next QP.
*        (If a unit step was not taken last iteration, some more
*        nonbasics may be between their bounds.
*        Nov 10, 1994. Tried leaving the nonbasics between their
*        bounds after short step. In some cases, the number of minor
*        iterations increased dramatically with a very short step.)
*        ---------------------------------------------------------------
         call dcopy ( nb, x1, 1, x, 1 )

         if ( nlnCon ) then
            call dcopy ( negCon, gCon1, 1, gCon, 1 )
            call dcopy ( nnCon , Lmul1, 1, Lmul, 1 )
            call dcopy ( nnCon , fCon1, 1, fCon, 1 )
         end if

         if ( nlnObj ) then
            fObj  = fObj1
            call dcopy ( nnObj, gObj1, 1, gObj, 1 )
         end if

         sInf = sInf1
         nInf = nInfQP      ! Not updated by the line search

         go to 100
*+    end while
      end if
*     ======================end of main loop============================

      return

 1000 format(1x, 29a4)
 1010 format(' Start of major itn', i6)
 1050 format(' Search exit', i3, ' -- ', a,
     &       '   Itn =', i7, '  Dual Inf =', 1p, e11.3)
 3020 format(' Itn', i7, ' -- Central differences invoked.',
     &       ' Small step length.' )

      end ! subroutine s8SQP

