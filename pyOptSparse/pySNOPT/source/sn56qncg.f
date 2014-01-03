*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     file  sn56qncg.f
*
*     s5QNdir   s5QN     s5QNit   s5Sswp   s5Hzx   s5Hzx1   s5Zswp
*     SYMMLQ    s5Msolv
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5QNdir
     &   ( Hprod, Hprod1,
     &     feasbl, Exactp, itn, QPterm, condHz, rgNorm,
     &     maxR, nS, lenR, R, g, p, w, gp,
     &     ne, nlocA, locA, indA, Acol,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1
      logical
     &     Exactp, feasbl, QPterm
      integer
     &     itn, lencu, lencw, leniu, leniw, lenR, lenru, lenrw, maxR,
     &     ne, nlocA, nS, locA(nlocA), indA(ne), iu(leniu), iw(leniw)
      double precision
     &     condHz, rgNorm, gp, Acol(ne), R(lenR), g(nS),
     &     p(nS), ru(lenru), rw(lenrw), w(nS)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     s5QNdir  computes a search direction  p  for the superbasic
*     variables, using the current reduced gradient  g.
*
*     16 Nov 2001: First version of s5QNdir.
*     26 Jul 2003: cgItn added.
*     14 Mar 2004: Argument QPterm added for lvlInf = 2.
*     04 Dec 2004: Current version.
*     ==================================================================
      character
     &     str*120
      external
     &     ddot, s5Hzx, s5Msolv
      logical
     &     checkA, goodb, SYMpre
      integer
     &     cgItn, cgItns, cgItmx, iStop, lprDbg,
     &     lr1, lr2, ls1, ls2, ls3, nout, PreCon, QPmode
      double precision
     &     ddot, Hznorm, rNorm, rtol, shift, tolCG, ynorm
*     ------------------------------------------------------------------
      integer            CG,            QN
      parameter         (CG        = 1, QN     = 2)
      integer            WithR,         WithRt
      parameter         (WithR     = 0, WithRt = 1)
      parameter         (cgItns    = 386) ! Total symmlq iterations
      parameter         (cgItn     = 387) ! symmlq itns for last minor
      double precision   zero,               one
      parameter         (zero      = 0.0d+0, one = 1.0d+0)
*     ------------------------------------------------------------------
      lprDbg    = iw( 85) ! >0    => private debug print
      cgItmx    = iw(111) ! CG iteration limit
      nout      = iw(151) ! unit # for printed messages
      QPmode    = iw(208) ! Current QP solver   (based on QPslvr)
      PreCon    = iw(209) ! Current precon mode (based on QPslvr)

      lr1       = iw(353) ! r1(maxS) SYMMLQ real work vector
      lr2       = iw(354) ! r1(maxS) SYMMLQ real work vector
      ls1       = iw(355) ! s1(maxS) SYMMLQ real work vector
      ls2       = iw(356) ! s2(maxS) SYMMLQ real work vector
      ls3       = iw(357) ! s3(maxS) SYMMLQ real work vector

      tolCG     = rw( 54) ! cg tolerance

      if (lprDbg .eq. 0)
     &nout      =  0

      checkA    = .false.
      Exactp    = .true.
      goodb     = .false.
      SYMpre    = .false.          ! No preconditioning done by SYMMLQ
      iw(cgItn) =  0

      shift     =  -1.0d-6         ! symmlq solves with  (A - shift*I) !!!
      rtol      = tolCG*min( one, rgNorm )

*     Steepest-descent direction.

      call dcopy
     &   ( nS, g, 1, p, 1 )

      if (feasbl  .and.  QPterm) then

*        Preconditioned CG.  Save  w  such that  R'*w = rg.

         if (QPmode .eq. CG  .and.  PreCon .eq. 1  .or.
     &       QPmode .eq. QN                            ) then
            call s6Rsol
     &           ( WithRt, maxR, nS, lenR, R, p )
         end if

         call dcopy
     &      ( nS, p, 1, w, 1 )

         if (QPmode .eq. QN  .and.  nS .gt. maxR  .or.
     &       QPmode .eq. CG                            ) then
            Exactp = .false.
            call SYMMLQ
     &         ( nS, w, rw(lr1), rw(lr2), rw(ls1), rw(ls2),p,rw(ls3),
     &           s5Hzx, s5Msolv, checkA, goodb, SYMpre, shift,
     &           nout , cgItmx, rtol,
     &           iStop, iw(cgItn), Hznorm, condHz, rNorm, ynorm,
     &           Hprod, Hprod1,                    ! Added for SNOPT
     &           ne , nlocA, locA, indA, Acol,     ! Added for SNOPT
     &           cu, lencu, iu, leniu, ru, lenru,  ! Added for SNOPT
     &           cw, lencw, iw, leniw, rw, lenrw ) ! Added for SNOPT
            if (iStop .eq. -1) then
               call dcopy
     &            ( nS, w, 1, p, 1 )
            end if
         end if

         call dcopy
     &      ( nS, p, 1, w, 1 )

         if (QPmode .eq. CG  .and.  PreCon .eq. 1  .or.
     &       QPmode .eq. QN                            ) then
            call s6Rsol
     &         ( WithR , maxR, nS, lenR, R, p )
         end if
      end if ! feasbl

*     Change the sign of p.

      call dscal ( nS, (-one), p, 1 )
      gp   = ddot  ( nS, g, 1, p, 1 )

      if (gp .ge. zero) then
         write(str, 1000) itn, gp, rtol
         call snPRNT( 23, str, iw, leniw )
      end if

      iw(cgItns) = iw(cgItns) + iw(cgItn )

      return

 1000 format(' Itn', i7, ': CG gives  gp = ', 1p, e8.1,
     &                   ' rtol = ', e8.1 )

      end ! subroutine s5QNdir

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5QN
     &   ( iExit, Prob, contyp, Elastc, subopt,
     &     Hprod, Hprod1, QPlog, gotR, needLU, typeLU, needx,
     &     lenR, m, maxS, mBS, n, nb, nDegen, Hcalls,
     &     ngQP0, ngQP, ngObj0, ngObj, nnH0, nnH, nS,
     &     itQP, itQPmax, itQPTgt, itn, lEmode, lvlInf, PrtLvl,
     &     minimz, iObj, sclObj, ObjAdd, ObjQP, condHz,
     &     targtZ, tolFP, tolQP, tolx, nInf, sInf, wtInf, piNorm,
     &     ne, nlocA, locA, indA, Acol,
     &     hElast, hEstat, hfeas, hs, kBS,
     &     Ascale, bl, bu, blBS, buBS,
     &     gBS, gObj, gQP, Hdx, pBS, pi, R, rc, rg, rg2,
     &     nrhs0, nrhs, rhs, lenx0, nx0, x0, x, xBS, xFreez,
     &     iy, iy1, y, y1, y2,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     Elastc, gotR, needLU, needx
      character
     &     contyp*20
      external
     &     Hprod, Hprod1, QPlog
      integer
     &     Hcalls, iExit, itQP, itQPmax, itQPTgt, itn, iObj,
     &     lencu, leniu, lenru, lencw, leniw, lenrw, lenR, lenx0,
     &     lEmode, lvlInf, PrtLvl, m, maxS, mBS, minimz, n, nb, ne,
     &     nDegen, ngQP0, ngQP, ngObj0, ngObj, nInf, nlocA, nnH0, nnH,
     &     nS, nrhs0, nrhs, nx0, Prob, subopt, typeLU, locA(nlocA),
     &     indA(ne), hElast(nb), hEstat(nb), hs(nb), kBS(mBS),
     &     hfeas(mBS), iy(nb), iy1(nb), iu(leniu), iw(leniw)
      double precision
     &     condHz, ObjAdd, ObjQP, piNorm, sclObj, sInf, tolFP, tolQP,
     &     tolx, targtZ, wtInf, Acol(ne), Ascale(nb), bl(nb), bu(nb),
     &     rc(nb), blBS(mBS), buBS(mBS), gBS(mBS), gObj(*), gQP(ngQP0),
     &     Hdx(nnH0), pBS(mBS), pi(m), rhs(nrhs0), R(lenR),
     &     rg(maxS), rg2(maxS), x0(lenx0), x(nb), xBS(mBS), xFreez(nb),
     &     y(nb), y1(nb), y2(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s5QN   solves a linear or quadratic program.  If the objective
!     is quadratic, a BFGS quasi-Newton method is used.
!
!     The problem type can be:
!       Prob = 0 FP   feasible point only
!       Prob = 1 LP   LP problem
!       Prob = 2 QP   QP problem
!       Prob = 3 FPE  feasible point for equalities only
!       Prob = 4 FPS  feasible point for QP subProblem
!       Prob = 5 QPS  QP subproblem
!       Prob = 6 QPP  FP subproblem with proximal point objective
!
!     ngQP = max( nnH, ngObj )
!
!     The optimization can pass through the following phases:
!
!       Phase 1               find a feasible point for all variables
!
!       Elastic Phase 1       make the non-elastic variables feasible
!                             while allowing infeasible elastics
!
!       Phase 2               minimize the objective
!
!       Elastic Phase 2       minimize a composite objective while
!                             keeping the non-elastics feasible
!
!                             In this phase, lvlInf means the following:
!
!                 lvlInf = 0  zero     weight on the infeasibilities
!                                      (infeasibillities are ignored)
!                          1  finite   weight on the infeasibilities
!                          2  infinite weight on the infeasibilities
!                                      (the objective is ignored)
!
!     The array kBS is a permutation on the column indices.
!     kBS(1  :m )    holds the col. indices of the basic variables.
!     kBS(m+1:m+nS)  holds the col. indices of the superbasic variables.
!                    These nS columns indices must have hs(j) = 2.
!
!     On exit:
!
!      iExit       Result
!      -----       ------
!       >0         Fatal error
!        0         QP solution found
!       -1         QP is infeasible
!       -2         QP is unbounded
!       -3         Too many iterations
!       -4         Void
!       -5         Too many superbasics
!       -6         Void
!       -7         Void
!       -8         Ill-conditioned CG null-space basis
!       -9         Too many subspace iterations
!
!     22 Jun 2001: First version of s5QN  based on Qpsol routine qpcore.
!     29 Jul 2003: Make sure R is never too ill-conditioned.
!     30 Jul 2003: Superbasic slacks allowed.
!     02 Aug 2003: snEXIT and snPRNT adopted.
!     04 Jul 2005: Current version of s5QN.
!     ==================================================================
      character
     &     str*120
      external
     &     dnormi
      logical
     &     badCG, checkx, chkFea, chkpi, conv(3), cnvrgd, done, feasbl,
     &     gotE, gotgQP, gotH, Hitcon, incres, jstFea, jstPh1, LUok,
     &     needf, needv, needR, newB, newLU, newSB, newx, optiml,
     &     parHes, prt10, prtLog, prtSum, QPdone, statpt,
     &     usegQP, Unbndd
      integer
     &     cgItn, condZ, iPrint, iSumm, inform, itnfix, itnlim,
     &     jq, jBq, jBr, jSq, jSr, kchk, kDegen, kfac, klog, kObj,
     &     ksav, kSumm, kp, kPrc, kPrPrt, linesP, linesS, lenL0, lenL,
     &     lenU0, lenU, LUitn, lvlTol, LUsiz0, LUmax, LUreq, maxR,
     &     MnrHdP, MnrHdS, mNewSB, modLU, mUncon, nBS, nElast, nFac,
     &     nfix(2), nfmove, nFreez, nInfE, LUmod, nonOpt, nSmax, nSwap,
     &     nUncon, PreCon, PrintP, PrintS, QPmode, Status,
     &     toldj1, toldj2, toldj3
      double precision
     &     Anorm, Bgrwth, Bold, deltaf, deltax, djq0,
     &     djq, djqmod, djqPrt, dnormi, dRmax, dRmin, eps0, etarg,
     &     featol, fObj0, fObj, ftoli, ftol(2), gSNorm, Hcndbd,
     &     ObjPrt, infBnd, pivot,  pSNorm, pSNrm1, rgNorm, rgtol(2),
     &     rowerr, sgnObj, sInfE, step, tolinc, tolrg, tolx0,
     &     weight, ObjSlk, xSNorm, xSNrm1, xtoli, xtol(2)
*     ------------------------------------------------------------------
      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one = 1.0d+0)
      integer            loose,      tight
      parameter         (loose  = 1, tight  = 2)
      integer            CG,         QN
      parameter         (CG     = 1, QN     = 2)
      integer            Intern,     Extern
      parameter         (Intern = 0, Extern = 1)
      integer            Check
      parameter         (Check  = 1)
      integer            WithB
      parameter         (WithB  = 1)
      integer            FP,         FPE,        FPS
      parameter         (FP     = 0, FPE    = 3, FPS   = 4)
      integer            BS,         BT
      parameter         (BS     = 2, BT     = 3)
      integer            Init,       Optml,      Cycle
      parameter         (Init   = 0, Optml  = 1, Cycle = 2)
      parameter         (condZ  = 192) ! condition estimate of Z
      parameter         (lenL0  = 171) ! size of L0
      parameter         (lenU0  = 172) ! size of initial  U
      parameter         (lenL   = 173) ! size of current  L
      parameter         (lenU   = 174) ! size of current  U
      parameter         (toldj1 = 184) ! phase 1 dj tol for p.p.
      parameter         (toldj2 = 185) ! phase 2 dj tol for p.p.
      parameter         (toldj3 = 186) ! current optimality tol
      parameter         (Status = 198) ! Status of a call to Hprod
      parameter         (kObj   = 205) ! xBS(kObj) is the obj. slack
      parameter         (LUitn  = 215) ! itns since last factorize
      parameter         (LUmod  = 216) ! number of LU mods
      parameter         (modLU  = 217) ! > 0  => modify the LU factors
      parameter         (PrintP = 218) ! (on/off) log     status
      parameter         (PrintS = 219) ! (on/off) summary status
      parameter         (linesP = 220) ! # lines in log     file
      parameter         (linesS = 221) ! # lines in summary file
      parameter         (MnrHdP = 223) ! >0 => Minor heading for iPrint
      parameter         (MnrHdS = 225) ! >0 => Minor heading for iSumm
      parameter         (cgItn  = 387) ! symmlq itns for last minor
*     ------------------------------------------------------------------
      eps0      = rw(  2) ! eps**(4/5)       IEEE DP  3.00e-13
      infBnd    = rw( 70) ! definition of an infinite bound
      etarg     = rw( 83) ! rgNorm tolerance
      Hcndbd    = rw( 85) ! bound on the condition of Hz

      iPrint    = iw( 12) ! Print file
      iSumm     = iw( 13) ! Summary file

      maxR      = iw( 52) ! max columns of R
      kchk      = iw( 58) ! check (row) frequency
      kfac      = iw( 59) ! factorization frequency
      ksav      = iw( 60) ! save basis map
      klog      = iw( 61) ! log/print frequency
      kSumm     = iw( 62) ! Summary print frequency
      kDegen    = iw( 63) ! max. expansions of featol
      itnlim    = iw( 89) ! limit on total iterations
      mNewSB    = iw( 95) ! max # of new superbasics
      nFac      = iw(210) ! # of LU factorizations
      QPmode    = iw(208) ! Current QP solver
      PreCon    = iw(209) ! Current precon mode (based on QPslvr)

      iExit     = 0
      mUncon    = 20

      if (nFac .gt. 0) then
         LUsiz0    = iw(lenL0) + iw(lenU0)
         LUmax     = 2*LUsiz0
      end if

      prt10  =              PrtLvl .ge. 10
      prtLog =              PrtLvl .ge.  1  .and.
     &         (mod( itQP,  klog ) .eq.  0  .and.  itQP   .ne. 0  .or.
     &                      klog   .eq.  1                             )
      prtSum =              PrtLvl .ge.  1  .and.
     &         (mod( itQP,  kSumm) .eq.  0  .and.  itQP   .ne. 0  .or.
     &                      kSumm  .eq.  1                             )
      iw(PrintP) = 0
      iw(PrintS) = 0
      if (prtLog) iw(PrintP) = 1
      if (prtSum) iw(PrintS) = 1

      kPrc   = 0                ! last sec scanned in part. prc
      LUreq  = 0

*     ------------------------------------------------------------------
*     s5QN operates in either ``Normal'' or ``Elastic'' mode.
*     Everything is normal unless a weighted sum is being minimized or
*     the constraints are infeasible.
*     The logical feasbl refers to the non-elastic variables.
*     Note that feasbl can be false while in elastic mode.
*     wtInf  is the optional parameter Infeasibility Weight.
*     ------------------------------------------------------------------
      feasbl = .false.
      gotE   = .false.

      gotH   =  nnH .gt. 0
      gotgQP = ngQP .gt. 0

!     jstPh1 = stop at the end of phase 1 (in normal or elastic mode)

      jstPh1 = Prob .eq. FP  .or.  Prob .eq. FPE  .or.   Prob  .eq. FPS

!     The phase 2 objective is F1 + wtInf*F2.

      if (Elastc) then
         needf = lvlInf .ne. 2  ! F1 required in phase 2
         needv = lvlInf .ne. 0  ! F2 required in phase 2
      else
         needf = .true.
         needv = .false.
      end if

!     needR  = true if R needs to be set or reset.

      needR  =         .not. jstPh1
     &         .and.  (QPmode .eq. QN  .or.
     &                 QPmode .eq. CG  .and.  PreCon .eq. 1)
      parHes =     nS .gt. maxR  .and.  gotR

*     call dcopy ( 2, rw( ftol1), 1,  ftol, 1 )
*     call dcopy ( 2, rw( xtol1), 1,  xtol, 1 )
*     call dcopy ( 2, rw(rgtol1), 1, rgtol, 1 )

      xtol(1)  = 0.1d+0
      xtol(2)  = 1.0d-6
      ftol(1)  = xtol(1)*0.1d+0
      ftol(2)  = xtol(2)**2
      rgtol(1) = 1.0d-3
      rgtol(2) = 1.0d-7

      ObjQP  = zero
      pivot  = zero
      step   = zero
      sInf   = zero
      tolrg  = zero
      nInfE  = 0
      jq     = 0
      djq    = zero
      djq0   = zero
      djqPrt = zero
      ObjSlk = zero

      jBq    = 0                ! x(jBq) is the incoming   BS
      jBr    = 0                ! x(jBr) is the outgoing   BS
      jSq    = 0                ! x(jSq) is the incoming SBS
      jSr    = 0                ! x(jSr) is the outgoing SBS
      lvlTol = loose
      kPrPrt = 0
      sgnObj = minimz

      iw(cgItn ) = 0
      rw(toldj1) = 100.0d+0
      rw(toldj2) = 100.0d+0

      chkFea = .true.
      chkpi  = .true.
      cnvrgd = .false.
      newLU  = .true.
      newx   = .false.
      Unbndd = .false.
      QPdone = .false.

*     nUncon  counts the number of unconstrained (i.e., Newton) steps.
*             If the test for a minimizer were scale-independent,
*             Uncon would never be larger than 1.
*     nfmove  counts the number of times that the QP obj is decreased,

      nfmove = 0
      nUncon = 0

*     subopt nonzero implies that optimization occurs with a subset of
*     the variables frozen at their initial values.
*     During suboptimization, nFreez is the number of frozen variables.

      nFreez = 0
      nSmax  = nS + mNewSB
      call s5hs
     &   ( Intern, nb, bl, bu, hs, x )
      call s5dgen
     &   ( inform, Init, PrtLvl, nb, nInf, itn,
     &     featol, tolx, tolinc, hs, bl, bu, x,
     &     itnfix, nfix, tolx0, iw, leniw, rw, lenrw )

**    ======================Start of main loop==========================
*+    do while (.not. QPdone  .and.  iExit .eq. 0)
  100 if       (.not. QPdone  .and.  iExit .eq. 0) then
*        ===============================================================
*        Check the initial  x  and move it onto  ( A  -I )*x = b.
*        If needLU is true, this will require a basis factorization.
*        ===============================================================
*        If necessary,  factorize the basis  ( B = LU ) and set x.
*        If needLU is false on entry to s5QN, the first call to s2Bfac
*        will try to use existing factors.
*        If needLU is true on entry to s5QN, an LU factorization of
*        type typeLU is computed.
*
*        The reason for the requested LU is as follows.
*
*        LUreq =  0  First LU for a given subproblem
*        LUreq =  1  Frequency
*        LUreq =  2  LU nonzeros increased
*        LUreq =  3
*        LUreq =  4
*        LUreq =  5  Singular after LU mod
*        LUreq =  6  Unstable LU mod (growth in new column of U)
*        LUreq =  7  Not enough memory
*        LUreq =  8
*        LUreq =  9
*        LUreq = 10  Row error in setx
*        LUreq = 11  Big  dx   in setx
*
*        LUreq = 20
*        LUreq = 21  Iterative refinement failed in QP
*        LUreq = 22  Unbounded QP
*        LUreq = 23  Infeasibility after refactorization
*        LUreq = 24  Small directional derivative in QP
*        LUreq = 25  Ill-conditioned Z in QP
*        LUreq = 26  Indefinite Z'HZ in QP
*        LUreq = 27  R singular after bound swap in QP
*        ---------------------------------------------------------------
         jstFea = .false.
         if (LUreq .gt. 0) needLU = .true.

         if (needx  .or.  needLU) then
            call s2Bfac
     &         ( iExit, typeLU, needLU, newLU, newB,
     &           iObj, itn, PrtLvl, LUreq,
     &           m, mBS, n, nb, nnH, nS, nSwap,
     &           ne, nlocA, locA, indA, Acol,
     &           kBS, hs, bl, bu, blBS, buBS,
     &           nrhs0, nrhs, rhs, x, xBS,
     &           iy, iy1, y, y1,
     &           iw, leniw, rw, lenrw )
            if ( newLU ) then
               LUsiz0 = iw(lenL0) + iw(lenU0)
               LUmax  = 2*LUsiz0
               if (nSwap .gt. 0) gotR = .false. ! Reset R.
               if (prt10) iw(MnrHdP) = 1
            end if

            if (iExit .ne. 0) go to 100

            parHes = nS .gt. maxR  .and.  gotR

            cnvrgd = .false.
            gotE   = .false.    ! Check hEstat in elastic mode.
            needx  = .false.
            newx   = .true.
            chkFea = .true.
            chkpi  = .true.

            pivot  = zero
            nUncon = 0
         end if

         nBS    = m + nS
         newSB  = .false.
         nInf   = 0
         sInf   = zero
         optiml = .false.

         if (Elastc  .and.  .not. gotE) then
*           ------------------------------------------------------------
*           Reset blBS and buBS for any violated elastics.
*           These values are used in s5step.
*           Strictly feasible elastics are returned to normality.
*           ------------------------------------------------------------
            call s5setE
     &         ( nBS, nb, nElast, featol, infBnd,
     &           hElast, hEstat, kBS,
     &           bl, bu, blBS, buBS, xBS )
            gotE  = .true.
         end if

         if ( chkFea ) then

*           In Phase 1 or just after a factorize, check the feasibility
*           of the basic and superbasic non-elastics.
*           jstFea  indicates that we have just become feasible.
*           jstFea is turned off once a step is taken.

            call dload
     &         ( nBS, zero, gBS, 1 )
            call s5Inf
     &         ( nBS, featol, nInf, sInf, hfeas, blBS, buBS, gBS, xBS )

            if (nInf .gt. 0) then

*              Non-elastics are infeasible.
*              If necessary, switch back to the feasibility phase, after
*              refactorization.
*              Print something if the basis has just been refactorized.

               if ( feasbl ) then
                  call s2tryLU
     &               ( itn, 23, nS, LUreq, LUok, typeLU,
     &                 iw, leniw, rw, lenrw )
                  if (.not. LUok) iExit = 11
                  feasbl = .false.
                  go to 100
               end if

               if (prt10  .and.  iw(LUitn) .eq. 0) then
                  write(str, 1030) itn, nInf, sInf
                  call snPRNT( 21, str, iw, leniw )
               end if
            end if

*           Feasbl = true  means that the non-elastics are feasible.
*                    This defines the start of Phase 2.

            if (.not. feasbl) then
               jstFea = nInf .eq. 0
            end if
            feasbl = nInf .eq. 0
            chkFea = nInf .gt. 0
         end if ! if chkFea

         if ( Elastc ) then
*           ------------------------------------------------------------
*           Compute the sum of infeasibilities of the elastic variables.
*           ------------------------------------------------------------
            call s5InfE
     &         ( nb, nBS, hEstat, kBS, nInfE, sInfE, bl, bu, x )
            nInf = nInf + nInfE
            sInf = sInf + sInfE
         end if

         if (feasbl  .and.  jstPh1) then
!           ------------------------------------------------------------
!           The non-elastic variables just became feasible. Exit.
!           ------------------------------------------------------------
            condHz = zero
            djqPrt = zero
            rgNorm = zero
            piNorm = zero
            call dload ( m, zero, pi, 1 )
            optiml = .true.

         else
            if ( feasbl ) then
!              ---------------------------------------------------------
!              The nonelastics are feasible.
!              (Elastc = false means no elastics.)
!              ---------------------------------------------------------
!              The objective row (if defined) is computed in all phases.

               ObjSlk = zero

               if (iObj .ne. 0) then
                  ObjSlk = xBS(iw(kObj))*sclObj
               end if

               if (needf  .and. (jstFea  .or.  newx)) then
!                 ======================================================
!                 If just feasible or x has changed, initialize
!                 the QP objective, the QP gradient and R.
!                 ObjQP = quadratic term to be minimized (or maximized).
!                 The objective row and ObjAdd do not appear in ObjQP.
!                 ObjQP is updated after each QP step.
!                 ======================================================
                  if ( gotgQP ) then
                     call s5QPfg
     &                  ( Hprod, Hprod1,
     &                    ngQP, ngObj0, ngObj, nnH,
     &                    iw(Status), Hcalls, ObjQP,
     &                    gObj, gQP, lenx0, nx0, x0, x, y,
     &                    cu, lencu, iu, leniu, ru, lenru,
     &                    cw, lencw, iw, leniw, rw, lenrw )
                     iw(Status) = 0
                  end if

                  if (gotH  .and. needR  .and.  .not. gotR) then
                     call s6Rset
     &                  ( gotR, maxR, nS, lenR, R, y, condHz )
*                    gotR = .true. ! done by s6Rset
                     parHes = nS .gt. maxR
                  end if

*                 ------------------------------------------------------
*                 Gather the QP gradient in BS order.
*                 Assign the nonzero components of gBS.
*                 ------------------------------------------------------
                  if ( gotgQP ) then
                     call s2gathr
     &                  ( ngQP, nBS, kBS, sgnObj, gQP, gBS )
                  end if
                  if (iObj .gt. 0) gBS(iw(kObj)) = sgnObj*sclObj

                  if ( Elastc  .and.  nInfE .gt. 0) then
                     call s5grdE
     &                  ( nb, nBS, wtInf, hEstat, kBS, gBS )
                  end if
               end if ! jstFea .or. newx

*              ---------------------------------------------------------
*              See if it's time to suboptimize.
*              NOTE: We must not suboptimize if all steps have been
*              degenerate.
*              ---------------------------------------------------------
               if (subopt .ne. 0  .or.  nfmove .eq. 0) then
*                 Relax
               else
                  if (nS   .ge. nSmax ) then
                     subopt = 1
                     if (prt10) then
                        write(str, 1610) itn, mNewSB
                        call snPRNT( 21, str, iw, leniw )
                     end if
                  else if (itQP .ge. itQPTgt) then
                     subopt = 2
                     if (prt10) then
                        write(str, 1620) itn, itQPTgt
                        call snPRNT( 21, str, iw, leniw )
                     end if
                  end if
               end if
            end if ! feasible

*           ============================================================
*           Check for an approximate stationary point.
*           If the gradient has changed, compute the reduced  gradient
*           ============================================================
            if (.not. feasbl  .or.  jstFea  .or.  newx) then
               call dcopy
     &            ( m, gBS, 1, y, 1 )
               call s5setp
     &            ( inform, m, chkpi, pinorm, y, pi,
     &              iw, leniw, rw, lenrw )
               if (inform .ne. 0) then
                  if (inform .gt. 0) then
                     iExit  = inform
                  else          ! pi is infinite or contains a NaN/Inf.
                     call s2tryLU
     &                  ( itn, 11, nS, LUreq, LUok, typeLU,
     &                    iw, leniw, rw, lenrw )
                     if (.not. LUok) iExit = 43
                  end if
                  go to 100
               end if

               rgNorm = zero
               if (nS .gt. 0) then
                  call s5rg
     &               ( m, nBS, n, nS, eps0,
     &                 ne, nlocA, locA, indA, Acol,
     &                 gBS, pi, rg, rgNorm, kBS )
               end if
            end if ! .not. feasbl  .or.  jstFea  .or.  newx

            if (tolrg .eq. zero) tolrg = etarg*rgNorm

            if ( feasbl ) then
               rw(toldj3) = tolQP
            else
               rw(toldj3) = tolFP
            end if

            if ( parHes ) then
               gSNorm = dnormi( maxR, rg, 1 )
            else
               gSNorm = rgNorm
            end if

            statpt = gSNorm .le. 0.1d+0  *tolrg   .or.
     &               gSNorm .le. rgtol(2)*piNorm  .or.  cnvrgd

            if (parHes  .and.  statpt) then

*              Swap the largest reduced gradient in Z2 into the front of
*              Z2 and see if it is significantly large.
*              Reset statpt if necessary.

               call s5Zswp
     &            ( statpt, m, maxR, maxS, nBS, lenR, nS,
     &              tolrg, kBS, blBS, buBS, gBS, R, rg, xBS, rw, lenrw )
            end if

            if ( gotR ) then
               call s6Rcnd
     &            ( maxR, nS, lenR, R, dRmax, dRmin, condHz )
               if (condHz .gt. Hcndbd) then
                  gotR   = .false.
                  condHz = zero
                  if ( needR ) then
                     call s6Rset
     &                  ( gotR, maxR, nS, lenR, R, y, condHz )
*                    gotR = .true. ! done by s6Rset
                  end if
               end if
            end if

            kPrPrt = kPrc
            jq     = 0

            if (statpt) then
*              ---------------------------------------------------------
*              Compute Lagrange multipliers.
*              ---------------------------------------------------------
               djq0   = djq     ! save djq in case of bad statpt
               djq    = zero
               nUncon = 0
               usegQP = feasbl  .and.  gotgQP
               weight = zero
               if (Elastc  .and.  feasbl) then
                  weight = wtInf
               end if

               call s5pric
     &            ( Elastc, feasbl, incres, usegQP, subopt,
     &              itn, m, n, nb, ngQP0, ngQP, nnH,
     &              nS, nFreez, nonOpt, weight, sgnObj, piNorm,
     &              jq, djq, kPrc, rw(toldj1),
     &              ne, nlocA, locA, indA, Acol,
     &              hElast, hs, bl, bu, gQP, pi, rc, x, xFreez,
     &              iw, leniw, rw, lenrw )

               optiml = nonOpt .eq. 0
               newSB  = nonOpt .gt. 0

               if ( lvlTol .eq. loose .and.
     &             (nS     .ge. maxS  .or.  optiml)) then
                  lvlTol = tight
                  tolrg  = rw(toldj3)*piNorm
                  if ( prt10 ) then
                     write(str, 1700) itn, tolrg
                     call snPRNT( 21, str, iw, leniw )
                  end if
                  if (rgNorm .gt. tolrg) then
                     optiml = .false.
                     newSB  = .false.
                  end if
               end if
            end if ! statpt
         end if ! jstPh1

         QPdone = optiml  .or.  Unbndd

         if ( QPdone ) then
*           ------------------------------------------------------------
*           Apparently we are finished.
*           See if any nonbasics have to be set back on their bounds.
*           ------------------------------------------------------------
            call s5dgen
     &         ( inform, Optml, PrtLvl, nb, nInf, itn,
     &           featol, tolx, tolinc, hs, bl, bu, x,
     &           itnfix, nfix, tolx0,
     &           iw, leniw, rw, lenrw )

            QPdone = inform .eq. 0

            if ( QPdone ) then
*              ---------------------------------------------------------
*              So far so good.  Now check the row residuals.
*              ---------------------------------------------------------
               if (iw(LUitn) .gt. 0) then
                  call s5setx
     &               ( inform, Check, itn,
     &                 m, n, nb, nBS, rowerr,
     &                 ne, nlocA, locA, indA, Acol,
     &                 kBS, xBS, nrhs0, nrhs, rhs, x, y, y2,
     &                 iw, leniw, rw, lenrw )

                  QPdone = inform .eq. 0
                  LUreq  = inform
                  if (LUreq .gt. 0) typeLU = BS
               end if
            end if

            if ( QPdone ) then
               if (Unbndd) iExit = -2
            else
               needx  = .true.
               Unbndd = .false.
               cnvrgd = .false.
               go to 100
            end if

            if (jstPh1) then
!              Relax, we are about to exit
            else
!              =========================================================
!              Print the details of the final iteration.
!              =========================================================
               ObjPrt = zero
               if ( feasbl ) then
                  if ( needf ) then
                     ObjPrt = ObjAdd + ObjSlk + ObjQP
                  end if
                  if ( needv ) then
                     ObjPrt = ObjPrt + sgnObj*wtInf*sInf
                  end if
               end if

               call QPlog
     &            ( Prob, contyp,
     &              Elastc, gotR, jstFea, feasbl,
     &              m, mBS, nnH, nS, jSq, jBr, jSr,
     &              iw(linesP), iw(linesS), itn, itQP, kPrPrt, lvlInf,
     &              pivot, step, nInf, sInf, wtInf,
     &              ObjPrt, condHz, djqPrt, rgNorm, kBS, xBS,
     &              iw, leniw )
            end if

            jBq       = 0
            jBr       = 0
            jSq       = 0
            jSr       = 0
            kPrPrt    = 0
            iw(cgItn) = 0
            djqPrt    = zero
            condHz    = zero

!           ------------------------------------------------------------
!           Convergence.
!           ------------------------------------------------------------
            if (nInf .gt. 0) then

!              No feasible point.
!              Stop or continue in elastic mode, depending on the
!              specified level of infeasibility.

               if (lEmode .gt. 0) then ! enter elastic mode
                  if (Elastc) then     ! already there
                     if (feasbl) then
                                ! Find the final sumInf for the elastics
                        call s5finE
     &                     ( nBS, nb, nInf, sInf, featol,
     &                       hEstat, kBS, bl, bu, xBS )
                     else
                                ! The nonelastics cannot be satisfied
                                ! by relaxing the elastics.  Exit.

                        iExit = -1 ! infeasible (should not happen)
                     end if
                  else             ! if .not. Elastc

*                    The constraints are infeasible in Normal mode.
*                    Print a message and start Elastic Phase 1.

                     if (prt10) then
                        write(str, 8050) itn, contyp
                        call snPRNT( 23, str, iw, leniw )
                        write(str, 8060) itn
                        call snPRNT( 23, str, iw, leniw )
                        iw(MnrHdP) = 1
                        iw(MnrHdS) = 1
                     end if

                     Elastc = .true.
                     needf  = lvlInf .ne. 2 ! F1 required in phase 2
                     needv  = lvlInf .ne. 0 ! F2 required in phase 2
                     QPdone = .false.
                     djq    = zero
                     step   = zero
                     call s5IniE
     &                  ( nb, nBS, nElast, featol, infBnd,
     &                    hElast, hEstat, kBS,
     &                    bl, bu, blBS, buBS, xBS )
                     gotE   = .true.
                  end if
                  go to 100
               end if
            end if

            if (prt10 .and. .not. jstPh1) then
               if (jq .ne. 0) then
                  djqprt = sgnObj*djq
                  if (prt10  .and.  klog .eq. 1) then
                     write(str, 1010) djq, jq, rgnorm, pinorm
                     call snPRNT( 31, str, iw, leniw )
                  end if
               else
                  if (prt10  .and.  klog .eq. 1) then
                     write(str, 1020)          rgnorm, pinorm
                     call snPRNT( 31, str, iw, leniw )
                  end if
               end if
            end if

         else ! not done
*           ============================================================
*           Take a series of reduced gradient steps until the new point
*           approximately minimizes the objective on the working set,
*           or lies on the boundary of a new constraint.
*+          ============================================================
*           Repeat                   (until  hitcon or cnvrgd or statpt)

  500       cnvrgd = .false.
            ObjPrt = zero
               if ( feasbl ) then
                  if ( needf ) then
                     ObjPrt = ObjAdd + ObjSlk + ObjQP
                  end if
                  if ( needv ) then
                     ObjPrt = ObjPrt + sgnObj*wtInf*sInf
                  end if
               end if

               call QPlog
     &            ( Prob, contyp,
     &              Elastc, gotR, jstFea, feasbl,
     &              m, mBS, nnH, nS, jSq, jBr, jSr,
     &              iw(linesP), iw(linesS), itn, itQP, kPrPrt, lvlInf,
     &              pivot, step, nInf, sInf, wtInf,
     &              ObjPrt, condHz, djqPrt, rgNorm, kBS, xBS,
     &              iw, leniw )
               jBq       = 0
               jBr       = 0
               jSq       = 0
               jSr       = 0
               kPrPrt    = 0
               iw(cgItn) = 0
               djqPrt    = zero
               condHz    = zero

               if ( newSB ) then
*                 ------------------------------------------------------
*                 A nonbasic has been selected to become superbasic.
*                 Compute the vector y such that B y = column jq.
*                 ------------------------------------------------------
*                 Set the level to which rgNorm must be reduced in the
*                 current subspace before we consider moving off another
*                 constraint.

                  djqmod = abs( djq )

                  rgNorm = max( rgNorm, djqmod )
                  tolrg  = etarg * djqmod

                  if (nS+1 .gt. maxS) then
                     iExit = -5
                     go to 100
                  end if

                  djqPrt = djq

*                 ------------------------------------------------------
*                 Compute the vector pBS such that B pB = column jq.
*                 pBS is a multiple of part of the new column of  Z  and
*                 is used to define the QN search direction.
*                 ------------------------------------------------------
*                 Unpack column jq into  y1  and solve  B*pB = y1.
*                 The altered  y1  satisfies  L*y1 = ajq.
*                 It is used later in s5QNit to modify L and U.

                  call s2unpk
     &               ( jq, m, n, ne, Anorm, nlocA, locA, indA, Acol, y1)
                  call s2Bsol
     &               ( iExit, WithB, m, y1, pBS, iw, leniw, rw, lenrw )
                  if (iExit .ne. 0) return

                  rw(condZ) = dnormi( m, pBS, 1 )/Anorm

c$$$                  if (rw(condZ) .ge. targtZ  .and.  nS .gt. 0) then
c$$$                     iExit = -8
c$$$                     go to 100
c$$$                  end if
               end if ! newSB

               if (itn  .ge. itnlim  .or.  itQP .ge. itQPmax) then
                  iExit = -3
                  go to 100
               end if

               itQP   = itQP   + 1
               itn    = itn    + 1

*              Decide if we want to print something this iteration.

               prtLog = iPrint .gt. 0  .and.
     &                  PrtLvl .ge. 1  .and.  mod( itQP, klog  ) .eq. 0
               prtSum = iSumm  .gt. 0  .and.
     &                  PrtLvl .ge. 1  .and.  mod( itQP, kSumm ) .eq. 0

               iw(PrintP) = 0
               iw(PrintS) = 0
               if (prtLog) iw(PrintP) = 1
               if (prtSum) iw(PrintS) = 1

               if ( feasbl ) then
                  fObj0 = sgnObj*(ObjSlk + ObjQP) + wtInf*sInf
               end if

*              ---------------------------------------------------------
*              Take a ``reduced gradient'' step.
*              The new  x  will either minimize the objective on the
*              working set or lie on the boundary of a new constraint.
*              ---------------------------------------------------------
               call s5QNit
     &            ( inform, Hprod, Hprod1, Elastc, feasbl,
     &              gotgQP, gotH, gotR, Hitcon, incres, needf, needv,
     &              newSB, itn, lenR,
     &              m, mBS, maxR, maxS, n, nb, Hcalls, nnH0, nnH, nS,
     &              ngQP0, ngQP, nDegen, LUreq, kp, jBq, jSq, jBr, jSr,
     &              jq, lvlTol, nfmove, nUncon,
     &              djq0, djq, minimz, iObj, sclObj, ObjQP,
     &              condHz, featol, pivot, piNorm, rgNorm, step, tolinc,
     &              nInf, sInf, wtInf, pSNorm, pSNrm1, xSNorm, xSNrm1,
     &              ne, nlocA, locA, indA, Acol,
     &              hElast, hEstat, hfeas, hs, kBS,
     &              bl, bu, blBS, buBS, gBS,
     &              gQP, Hdx, pBS, pi, rg, rg2, R,
     &              x, xBS, y, y1, y2,
     &              cu, lencu, iu, leniu, ru, lenru,
     &              cw, lencw, iw, leniw, rw, lenrw )

*              Check for trouble.
*              inform values are  -1, 0 >0

               if (inform .ne. 0) then
                  if      (inform .gt.  0) then
                     iExit = inform ! Fatal LU error
                  else if (inform .eq. -1) then
                     iExit = -2     ! unbounded
                  end if
                  go to 100
               end if

               nBS  = m + nS

*              ---------------------------------------------------------
*              Stop if there are errors.
*              ---------------------------------------------------------
               done  = LUreq .ne. 0
               badCG = nUncon .ge. mUncon  .and.  lvlTol .eq. tight

               if (.not. done) then
*                 ------------------------------------------------------
*                 Reset R if it is too ill-conditioned.
*                 ------------------------------------------------------
                  if ( gotR ) then
                     call s6Rcnd
     &                  ( maxR, nS, lenR, R, dRmax, dRmin, condHz )
                     if (condHz .gt. Hcndbd) then
                        gotR   = .false.
                        condHz = zero
                        if ( needR ) then
                           call s6Rset
     &                        ( gotR, maxR, nS, lenR, R, y, condHz )
*                          gotR = .true. ! done by s6Rset
                        end if
                     end if
                  end if

                  if (.not. Hitcon) then
*                    ---------------------------------------------------
*                    The QP step was unconstrained.
*                    Test for convergence on the current working set.
*                    ---------------------------------------------------
                     parHes    = nS .gt. maxR  .and.  gotR
                     gSNorm    = rgNorm

                     if ( parHes ) then
                        pSNorm = pSNrm1
                        xSNorm = xSNrm1
                        gSNorm = dnormi( maxR, rg, 1 )
                     end if

                     if (feasbl .and. needf) then
                        if (iObj .gt. 0) then
                           ObjSlk = xBS(iw(kObj))*sclObj
                        end if
                        fObj    = sgnObj*(ObjSlk + ObjQP) + wtInf*sInf
                        ftoli   = ftol(lvlTol)
                        xtoli   = xtol(lvlTol)

                        deltax  = step * pSNorm
                        deltaf  = fObj0 - fObj

                        conv(1) = deltax .le. xtoli*(one + xSnorm )
                        conv(2) = deltaf .le. ftoli*(one + abs(fObj0))
                        conv(3) = gSNorm .le. tolrg

                        cnvrgd  = conv(1) .and. conv(2) .and. conv(3)
                     end if

                     statpt  =  gSNorm  .le.    0.1d+0*tolrg   .or.
     &                          gSNorm  .le.  rgtol(2)*piNorm
                  end if ! .not. Hitcon
               end if

*+          until    (Hitcon .or. statpt .or. cnvrgd.or. done.or.badCG)
            if (.not.(Hitcon .or. statpt .or. cnvrgd.or. done.or.badCG))
     &      go to 500
*           ------------------------------------------------------------

            if ( badCG ) then
               iExit = -9       ! too many subspace iterations
               go to 100
            end if

            if (LUreq  .gt. 0) then
               call s2tryLU
     &            ( itn, LUreq, nS, LUreq, LUok, typeLU,
     &              iw, leniw, rw, lenrw )
               if (.not. LUok) then
                  iExit = 43
                  go to 100
               end if
            end if

            iw(LUitn) = iw(LUitn)  + 1 ! itns since the last LU
            newLU     = .false.
            newx      = .false.
            chkpi     = .false.

*           ============================================================
*           Test for error condition and/or frequency interrupts.
*           ============================================================
            if (mod(itn,ksav) .eq. 0) then
               call icopy
     &            ( nb, hs, 1, iy1, 1 )
               call s5hs
     &            ( Extern, nb, bl, bu, iy1, x )
               call s4ksav
     &            ( minimz, m, n, nb, nS, mBS,
     &              itn, nInf, sInf, ObjQP, kBS, iy1,
     &              Ascale, bl, bu, x, xBS, cw, lencw, iw, leniw )
            end if

*           Increment featol every iteration.

            featol = featol + tolinc

*           Every kdegen iterations, reset featol and move nonbasic
*           variables onto their bounds if they are very close.

            if (mod( itn, kdegen ) .eq. 0) then
               call s5dgen
     &            ( inform, Cycle, PrtLvl, nb, nInf, itn,
     &              featol, tolx, tolinc, hs, bl, bu, x,
     &              itnfix, nfix, tolx0,
     &              iw, leniw, rw, lenrw )
               needx  = inform .gt. 0
            end if

*           Refactorize the basis if it has been modified too much.

            if (LUreq .eq. 0) then
               if (     iw(LUmod) .ge. kfac-1) then
                  LUreq  = 1
               else if (iw(LUmod) .ge. 20  .and.
     &                 iw(lenL)+iw(lenU) .gt. LUmax) then
                  Bgrwth = iw(lenL) + iw(lenU)
                  Bold   = LUsiz0
                  Bgrwth = Bgrwth/Bold
                  if (prt10) then
                     write(str, 1000) Bgrwth
                     call snPRNT( 21, str, iw, leniw )
                  end if
                  LUreq  = 2
               end if
               if (LUreq .gt. 0) typeLU = BT
            end if

*           Update the LU factors.

            if (LUreq .eq. 0  .and.  iw(modLU) .gt. 0) then
               iw(LUmod) = iw(LUmod) + 1
               call s2Bmod2
     &            ( inform, kp, m, y1, iw, leniw, rw, lenrw )
               if (inform .eq. -1) LUreq = 5 ! Singular after LU mod
               if (inform .eq.  2) LUreq = 6 ! Unstable LU mod
               if (inform .eq.  7) LUreq = 7 ! Insufficient free memory
               if (LUreq  .gt.  0) then
                  call s2tryLU
     &               ( itn, LUreq, nS, LUreq, LUok, typeLU,
     &                 iw, leniw, rw, lenrw )
                  if (.not. LUok) then
                     iExit = 43
                     go to 100
                  end if
               end if
            end if

            if (LUreq .eq. 0) then
               checkx = mod(iw(LUitn),kchk) .eq. 0
               if (checkx  .and.  .not. needx) then
                  call s5setx
     &               ( inform, Check, itn,
     &                 m, n, nb, nBS, rowerr,
     &                 ne, nlocA, locA, indA, Acol,
     &                 kBS, xBS, nrhs0, nrhs, rhs, x, y, y2,
     &                 iw, leniw, rw, lenrw )
                  LUreq  = inform
                  cnvrgd = .false.
               end if
               if (LUreq .gt. 0) typeLU = BT
            end if
         end if

         go to 100
*+    end while
      end if
*     ======================end of main loop============================
*
      call s5hs  ( Extern, nb, bl, bu, hs, x )

      if (subopt .gt. 0) then
         if (nFreez .gt. 0) then
*           Relax
         else
            subopt = 0
         end if
      end if

      return

 1000 format(' ==> LU file has increased by a factor of', f6.1)
 1010 format(' Biggest dj =', 1p, e11.3, ' (variable', i7, ')',
     &       '    norm rg =',     e11.3, '   norm pi =', e11.3)
 1020 format(   ' Norm rg =', 1p, e11.3, '   norm pi =', e11.3)
 1030 Format(' Itn', i7, ': Infeasible nonelastics.  Num =', i5, 1p,
     &                   '  Sum of Infeasibilities =', e8.1)
 1500 format(' Itn', i7, ': Expanded reduced Hessian ',
     &                   'is indefinite. Basis refactorized')
 1610 format(' Itn', i7, ': Suboptimize: ', i7, ' new superbasics')
 1620 format(' Itn', i7, ': Suboptimize: ', i7, ' minor iterations')
 1700 format(' Itn', i7, ': Subspace tolerance = ', 1p, e11.3)
 8050 format(' Itn', i7, ': Infeasible ', a)
 8060 format(' Itn', i7, ': Elastic Phase 1 -- making ',
     &                   'non-elastic variables feasible')

      end ! subroutine s5QN

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5QNit
     &   ( iExit, Hprod, Hprod1, Elastc, feasbl,
     &     gotgQP, gotH, gotR, Hitcon, incres, needf, needv,
     &     newSB, itn, lenR,
     &     m, mBS, maxR, maxS, n, nb, Hcalls, nnH0, nnH, nS,
     &     ngQP0, ngQP, nDegen, LUreq, kp, jBq, jSq, jBr, jSr,
     &     jq, lvlTol, nfmove, nUncon,
     &     djq0, djq, minimz, iObj, sclObj, ObjQP,
     &     condHz, featol, pivot, piNorm, rgNorm, step, tolinc,
     &     nInf, sInf, wtInf, pSNorm, pSNrm1, xSNorm, xSNrm1,
     &     ne, nlocA, locA, indA, Acol,
     &     hElast, hEstat, hfeas, hs, kBS,
     &     bl, bu, blBS, buBS, gBS,
     &     gQP, Hdx, pBS, pi, rg, rg2, R,
     &     x, xBS, y, y1, y2,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1
      logical
     &     Elastc, feasbl, gotgQP, gotH, gotR, Hitcon, incres,
     &     needf, needv, newSB
      integer
     &     iExit, iObj, itn, lenR, LUreq, lvlTol, m, maxR, maxS,
     &     mBS, minimz, n, nb, nDegen, ne, nInf, nlocA,
     &     Hcalls, nnH0, nnH, ngQP0, ngQP, nS, kp, jBq, jBr, jq, jSq,
     &     jSr, nfmove, nUncon, lencu, lencw, leniu, leniw, lenru,
     &     lenrw, locA(nlocA), indA(ne), hElast(nb), hEstat(nb),
     &     hs(nb), hfeas(mBS), kBS(mBS), iu(leniu), iw(leniw)
      double precision
     &     condHz, djq0, djq, featol, ObjQP, pivot, pSNorm, pSNrm1,
     &     rgNorm, sclObj, sInf, step, tolinc, wtInf, xSNorm, xSNrm1,
     &     Acol(ne), bl(nb), bu(nb), blBS(mBS), buBS(mBS),
     &     gBS(mBS), gQP(ngQP0), Hdx(nnH0), pBS(mBS), pi(m),
     &     R(lenR), rg(maxS), rg2(maxS), x(nb), xBS(mBS),
     &     y(nb), y1(nb), y2(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     s5QNit performs a QP step.
*
*     On entry,
*        newSB = true implies that variable jq just went superbasic.
*                In this case:
*                pBS(1:m)  satisfies B pBS = a(jq).
*                y1(1:m)   satisfies L  y1 = a(jq).
*
*      iExit       Result
*      -----       ------
*       -1         unbounded
*        0         normal exit
*       >0         Fatal LU error

*     On exit,
*        pBS contains the most recent QP search direction.
*
*     incres says if the new variable should increase or not.
*            It is used only in linear mode when s5Pric is moving free
*            nonbasics toward their closest bound (and djq = zero).
*
*     gotR   says if a useful quasi-Newton preconditioner  R  exists.
*            It will be false if no preconditioning is being used.
*            If  preconditioned CG (QN) is being used, then gotR may
*            be true even if the current itn is infeasible.
*            This allows  R  to be updated in during a temporary loss
*            of feasibility.
*
*     parHes is true if  R  is not big enough to store a full
*            preconditioner for the reduced Hessian.
*            The null space is then  Z = ( Z1  Z2 ),  and  R  estimates
*            the Hessian only in the subspace  Z1.  A diagonal estimate
*            is used (also in  R) for the subspace  Z2.
*
*     12 Jun 2001: First   version of s5QNit based on s5Qpit
*                  and MINOS routine m7rgit.
*     14 Jul 2001: Different gp needed for parHes.
*     02 Aug 2003: snPRNT adopted.
*     30 Jun 2005: Current version of s5QNit.
*     ==================================================================
      character
     &     str*120
      external
     &     ddot, dnormi
      logical
     &     chkpi, hitlow, Exactp, move, onbnd, parHes, Unbndd, Uncon
      integer
     &     inform, infpiv, jEs, jEsq, jqStat,
     &     jr, jrStat, kBSq, kObj, kSq, m1, modR1, modLU,
     &     mtry, nBS, ntry, nZ1, nZ2, Status
      double precision
     &     Anorm, bigdx, bound, ddot, dnormi, eps, eps0, eps2, exact,
     &     gp, gpQP, ObjChg, pBS1, pHp, pHpQP, pNorm, piNorm, pSNrm2,
     &     infBnd, rgdel, sclPiv, sgnObj, StepB, stepmx, stepP,
     &     t, tolP0, tolP, tolpiv, xSNrm2
*     ------------------------------------------------------------------
      parameter         (kObj   = 205) ! xBS(kObj) is the obj. slack
      parameter         (modLU  = 217) ! > 0  => modify the LU factors
      integer            tight
      parameter         (tight  = 2)
      integer            Normal
      parameter         (Normal = 0)
      integer            xBStox
      parameter         (xBStox = 1)
      parameter         (mtry   = 6)
      integer            WithL,      WithB,      WithBt
      parameter         (WithL  = 0, WithB  = 1, WithBt = 2)
      double precision   zero,            half,          one
      parameter         (zero   = 0.0d+0, half = 0.5d+0, one = 1.0d+0)
      double precision   ten
      parameter         (ten    =10.0d+0)
*     ------------------------------------------------------------------
      eps       = rw(  1) ! unit round-off.  IEEE DP  2.22e-16
      eps0      = rw(  2) ! eps**(4/5)       IEEE DP  3.00e-13
      eps2      = rw(  4) ! eps**(1/2)       IEEE DP  1.49e-08
      tolpiv    = rw( 60) ! excludes small elements of pBS.
      infBnd    = rw( 70) ! definition of an infinite bound
      bigdx     = rw( 72) ! unbounded step.

      iExit     = 0
      inform    = 0
      Status    = 0

      chkpi     = .false.
      Unbndd    = .false.
      sgnObj    = minimz

      m1        = m + 1
      nBS       = m + nS

      if ( newSB ) then
*        ---------------------------------------------------------------
*        New superbasic.
*        ---------------------------------------------------------------
         nS    = nS   + 1
         nBS   = nBS  + 1

         kBS (nBS) =    jq
         xBS (nBS) = x (jq)
         blBS(nBS) = bl(jq)
         buBS(nBS) = bu(jq)
         jqStat    = hs(jq)

         if ( gotR ) then ! Add a unit column to R at position nS.
            call s6Radd
     &         ( maxR, lenR, nS, R )
         end if

         hfeas(nBS) = 0

         if (feasbl .and. needf .and. gotgQP .and. jq .le. ngQP) then
            gBS(nBS) = sgnObj*gQP(jq)
         else
            gBS(nBS) = zero
         end if

*        ===============================================================
*        Set hEstat(jq) and the elastic parts of blBS and buBS.
*        ===============================================================
         if ( Elastc ) then

*           If the new superbasic is an elastic variable
*           and it wants to move infeasible, set its elastic state.

            jEsq = hEstat(jq)

            if (hElast(jq) .gt. 0) then
               nInf = nInf + 1
               if ( incres ) then
                  if (jqStat .eq. 1  .or.  jqStat .eq. 4) then
                     hEstat(jq) =   2
                     buBS(nBS)  =   infBnd
                     if ( feasbl ) then
                        gBS (nBS) = gBS(nBS) + wtInf
                        blBS(nBS) = bu(jq)
                     end if
                  end if
               else
                  if (jqStat .eq. 0  .or.  jqStat .eq. 4) then
                     hEstat(jq) =   1
                     blBS(nBS)  = - infBnd
                     if ( feasbl ) then
                        gBS (nBS) = gBS(nBS) - wtInf
                        buBS(nBS) = bl(jq)
                     end if
                  end if
               end if
            end if
         end if ! Elastc

*        ---------------------------------------------------------------
*        In phase 1, or in phase 2 for an LP, price can select nonbasics
*        floating free between their bounds with zero reduced cost.
*        We have to check that dqj is not zero.
*        ---------------------------------------------------------------
         rg(nS) = djq
         if (.not. feasbl  .or.  (needf .and.  .not. gotH)) then
            if (hs(jq) .eq. -1) then
               if (incres) then
                  rg(nS) = - one
               else
                  rg(nS) =   one
               end if
            end if
         end if
         jSq    = jq
         hs(jq) = 2
      end if ! newSB

*     ------------------------------------------------------------------
*     Compute the search direction for the superbasics.
*     ------------------------------------------------------------------
      parHes = nS .gt. maxR  .and.  gotR

      if ( parHes ) then
         nZ1 = maxR
      else
         nZ1 = nS
      end if

      nZ2 = nS - maxR

*     ------------------------------------------------------------------
*     Store the free components of the search direction in pBS(1:nBS).
*     First, find the search direction pS for the superbasics.
*     The vector v such that R*v = rg is held in y2 for the BFGS update.
*     s5QNdir uses  y  as workspace.
*     ------------------------------------------------------------------
  100 call s5QNdir
     &   ( Hprod, Hprod1,
     &     feasbl, Exactp, itn, needf, condHz, rgNorm,
     &     maxR, nS, lenR, R, rg, pBS(m1), y2, gp,
     &     ne, nlocA, locA, indA, Acol,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (gp .ge. zero) then
         gotR  = .false.
         LUreq = 24
         go to 900
      end if

      pBS1 = pBS(m+nS)

      if (newSB  .and.  feasbl) then

*        Check that pBS(nBS) is feasible wrt blBS(nBS) and buBS(nBS).
*        A large rgtol may give a pBS(nBS) with the wrong sign.

         if (djq*pBS1 .gt. zero) then
            nS     = nS  - 1
            nBS    = nBS - 1
            hs(jq) = jqStat
            if (Elastc) hEstat(jq) = jEsq
            jq     =  0
            djq    =  djq0
            jSq    = -jSq
            newSB  = .false.
            if (lvlTol .eq. tight  .and.  gotR) then
               gotR   = .false.
               call s6Rset
     &            ( gotR, maxR, nS, lenR, R, y, condHz )
*              gotR   = .true.  ! done by s6Rset
            else
               lvlTol = tight
            end if
            go to 100
         end if
      end if

*     ------------------------------------------------------------------
*     Find the norm of  pS.  Put the search direction for the basic
*     variables in pBS(1)  ,...,pBS(m).
*     ------------------------------------------------------------------
      xSNrm1 = dnormi( nZ1, xBS(m1), 1 )
      pSNrm1 = dnormi( nZ1, pBS(m1), 1 )

      xSNrm2 = zero
      pSNrm2 = zero

      if ( parHes ) then
         xSNrm2 = dnormi( nZ2, xBS(m1+maxR), 1 )
         pSNrm2 = dnormi( nZ2, pBS(m1+maxR), 1 )
      end if

      pSNorm = max ( pSNrm1, pSNrm2 )
      xSNorm = max ( xSNrm1, xSNrm2 )

*     Compute  y = - S*pS and prepare to solve  B*pB = y for pB,
*     the search direction for the basic variables.
*     We first normalize  y  so the LU solver won't ignore
*     too many "small" elements while computing pB.

      call dscal
     &   ( nS, (one/pSNorm), pBS(m1), 1 )
      call s2Bprd
     &   ( Normal, eps0, n, nS, kBS(m1),
     &     ne, nlocA, locA, indA, Acol,
     &     (-one), pBS(m1), nS, zero, y, m )

*     Solve  B*pBS = y  and then unnormalize all of pBS.

      call s2Bsol
     &   ( iExit, WithB, m, y, pBS, iw, leniw, rw, lenrw  )
      if (iExit .ne. 0) go to 900
      call dscal
     &   ( nBS , pSNorm, pBS, 1 )
      pNorm = max( pSNorm, dnormi( m, pBS, 1 ) )

      if ( feasbl ) then
!        ---------------------------------------------------------------
!        Compute y = pBS(scattered) and Hdx(scattered).
!        The vector Hdx is used to update the objective and gradient of
!        the QP.  Form  gpQP  and  pHpQP  for the quadratic term.
!        gp = gpQP - pBS(kObj) + terms from the elastic gradient.
!        ---------------------------------------------------------------
         pHp = zero

         if (needf  .and.  (gotgQP  .or.  gotH)) then
            call s2scatr
     &         ( ngQP, nBS, kBS, one, pBS, y )

            if ( gotgQP ) then
               gpQP = ddot ( ngQP, gQP, 1, y, 1 )
            end if

            if ( gotH ) then
               pHpQP = zero

               call Hprod
     &            ( Hprod1, Hcalls, nnH,
     &              y, Hdx, Status,
     &              cu, lencu, iu, leniu, ru, lenru,
     &              cw, lencw, iw, leniw, rw, lenrw )
               pHpQP = pHpQP + ddot  ( nnH, y, 1, Hdx, 1 )
               pHp   = sgnObj*pHpQP
            end if
         end if
      end if ! feasbl

      if (gp .ge. zero) then
         write(str, 1000) itn, gp
         call snPRNT( 23, str, iw, leniw )
         LUreq  = 24
         go to 900
      end if

      newSB  = .false.

*     ------------------------------------------------------------------
*     Find the nearest constraint in direction  x + step*pBS (step > 0).
*     Exact  is the step that takes xBS(kp) exactly onto bound.
*     It may be positive or slightly negative. (Not defined if Unbndd.)
*
*     If onbnd  is true, step is a step that reaches a bound exactly.
*     xBS(kp) reaches the value bound.  If we take a constrained step,
*     bound is used to put the new nonbasic variable x(jr) exactly on
*     its bound.
*
*     If Unbndd is true, step = stepmx.
*     ------------------------------------------------------------------
      stepmx = bigdx /pNorm
      sclPiv = one
      tolP0  = tolpiv
      tolP   = tolpiv*pNorm
      ntry   = 0

*+    Repeat
  200    tolP  = tolP /sclPiv
         tolP0 = tolP0/sclPiv
         call s5step
     &      ( nBS, nDegen,
     &        featol, infBnd, stepmx, tolinc, tolP,
     &        hfeas, blBS, buBS, xBS, pBS,
     &        hitlow, move, onbnd, Unbndd,
     &        infpiv, kp, bound, exact, stepB, stepP )

*        Find if the step is constrained or unconstrained.
*        If R has been flagged as singular, we double check by trying
*        to compute the QP minimizer along pBS.  If the minimizer
*        exists,  the singularity tolerance must be too large.

         if ( feasbl ) then
            Uncon  = stepP*pHp .gt. (- gp)
            Unbndd = (Unbndd  .and.  .not. Uncon)  .or.  stepmx .le. one
         else
            Uncon = .false.
         end if

         sclPiv = ten
         ntry   = ntry + 1

!+    until    (infpiv .eq. 0 .and. (.not. Unbndd .or. feasbl) .or.
!+                ntry .ge. mtry)
      if (.not.(infpiv .eq. 0 .and. (.not. Unbndd .or. feasbl) .or.
     &            ntry .ge. mtry)) go to 200

      if ( Unbndd ) then
         iExit = -1
         go to 900
      end if

      Hitcon = .not. Uncon

      if ( Hitcon ) then
         nUncon = 0
         step   = stepB
         pivot  = - pBS(kp)
      else
         nUncon = nUncon + 1
         step   = (- gp)/pHp
         pivot  = zero
      end if

!     ------------------------------------------------------------------
!     Compute ObjChg, the predicted change in fObj.
!     Note: pHpQP and pHp are defined before and after scaling by
!     sgnObj.
!     Update the value and gradient of the quadratic obj term.
!     ------------------------------------------------------------------
      if ( feasbl ) then
         if (step .gt. zero) then
            ObjChg = step*gp + half*pHp*step**2

            if ( needf ) then
               if ( gotgQP )
     &            ObjQP = ObjQP + step*gpQP
               if ( gotH   ) then
                  ObjQP = ObjQP + half*pHpQP*step**2
                  if (step .gt. zero)
     &            call daxpy
     &               ( nnH, step, Hdx, 1, gQP, 1 )
               end if
            end if
         end if
      end if

      if (feasbl  .and.  move) nfmove = nfmove + 1

*     ------------------------------------------------------------------
*     Update the basic variables xBS.
*     ------------------------------------------------------------------
      call daxpy
     &   ( nBS, step, pBS, 1, xBS, 1 )
      call s5BSx
     &   ( xBStox, nBS, nb, kBS, x, xBS )

      if ( Hitcon ) then
      end if

      modr1  = 0

      if (feasbl  .and.  step .gt. zero) then
*        ===============================================================
*        gQP has been changed.
*        Compute the new gBS and reduced gradient rg2.
*        Update the reduced Hessian.
*        ===============================================================
         call dload ( nBS, zero, gBS, 1 )

         if ( needf ) then
            if ( gotgQP ) then
               call s2gathr
     &            ( ngQP, nBS, kBS, sgnObj, gQP, gBS )
            end if
            if (iObj .gt. 0) gBS(iw(kObj)) = sgnObj*sclObj
         end if

         if (Elastc  .and.  needv) then
            call s5InfE
     &         ( nb, nBS, hEstat, kBS, nInf, sInf, bl, bu, x )
            call s5grdE
     &         ( nb, nBS, wtInf, hEstat, kBS, gBS )
         end if

         call dcopy
     &      ( m, gBS, 1, y, 1 )
         call s5setp
     &      ( iExit, m, chkpi, pinorm, y, pi, iw, leniw, rw, lenrw )
         if (iExit .ne. 0) go to 900

         call s5rg
     &      ( m, nBS, n, nS, eps0,
     &        ne, nlocA, locA, indA, Acol,
     &        gBS, pi, rg2, rgNorm, kBS )

         if ( gotR )
     &   call s6Rqn
     &      ( modr1, Exactp, maxR, nS,
     &        lenR, R, gp, rg, rg2, pBS(m1), y2, step, eps2, eps0 )
      else
         call dcopy
     &      ( nS, rg, 1, rg2, 1 )
      end if

      if ( Uncon ) then
*        ===============================================================
*        The step is unconstrained.
*        ===============================================================
         call dcopy
     &      ( nS, rg2, 1, rg, 1 )
         iw(modLU) = 0          ! No LU update

      else ! hit constraint
*        ===============================================================
*        There is a blocking variable.
*        It could be a fixed variable, whose new state must be 4.
*        ===============================================================
         jr     =   kBS(kp)

         if (onbnd) x(jr) = bound

         jEs    = hEstat(jr)
         hEstat(jr) = 0

         if      (jEs .eq. 0) then
            if (blBS(kp) .eq. buBS(kp)) then
               jrstat = 4
            else if (hitlow) then
               jrstat = 0
            else
               jrstat = 1
            end if

         else if (jEs .eq. 1) then
            if (bl(jr) .eq. bu(jr)) then
               jrstat =   4
            else if (onbnd) then
               jrstat =   0
            else if (x(jr) .lt. bu(jr)) then
               jrstat = - 1
            else
               jrstat =   1
            end if

         else !   jEs .eq. 2
            if (bl(jr) .eq. bu(jr)) then
               jrstat =   4
            else if (onbnd) then
               jrstat =   1
            else if (x(jr) .gt. bl(jr)) then
               jrstat = - 1
            else
               jrstat =   0
            end if
         end if

         if (kp .le. m) then
*           ============================================================
*           A variable in B hit a bound.
*           Find column kSq = kBSq-m  of S to replace column kp of B.
*           If nS = 1, it must be the entering SB column.
*           ============================================================
*           if (nS .eq. 1) then
*              kBSq  = nBS
*              pivot = pivot/pBS1
*           else
               call dload
     &            ( m, zero, y2, 1 )
               y2(kp) = one     ! Set      y2 = ep
                                ! Solve  B'yB = ep
               call s2Bsol
     &            ( inform, WithBt, m, y2, y, iw, leniw, rw, lenrw )
               call s5chzq
     &            ( m, mBS, n, nb, nS, kBSq, pivot, tolP0,
     &              ne, nlocA, locA, indA, Acol,
     &              kBS, bl, bu, xBS, y, iw, leniw, rw, lenrw )
               if (kBSq .le. 0) then
                  write(str, 9999) itn
                  call snPRNT( 23, str, iw, leniw )
                  kBSq   = nBS
               end if
*           end if

            kSq        = kBSq - m

            hs(jr)     = jrStat
            jBr        = jr                     ! Outgoing basic
            jSr        = kBS(kBSq)              ! Outgoing superbasic
            kBS (kBSq) = jBr
            jBq        = jSr                    ! Incoming basic
            kBS (kp)   = jSr
            blBS(kp)   = blBS(kBSq)
            buBS(kp)   = buBS(kBSq)
            xBS (kp)   = xBS (kBSq)
            gBS (kp)   = gBS (KBSq)
            hs(jBq)    = 3

*           Finish computing yS = (y(m+1), ..., y(m+nS)).

            y(kBSq) = - (one + pivot)
            call dscal
     &         ( nS, (one/pivot), y(m1), 1 )

            if ( gotR  .and.  kSq .le. maxR) then
               call s6Rswp
     &            ( maxR, nZ1, lenR, R, y2, y(m1), kSq, eps0 )
            end if

*           Modify  pi  using  y  where  B' y = e(kp).

            t      = rg2(kSq) / pivot
            call daxpy
     &         ( m, t, y, 1, pi, 1 )
            piNorm = max( dnormi( m, pi, 1 ), one )

*           ------------------------------------------------------------
*           Get a new  y1, used to modify L and U.  If the outgoing
*           superbasic just came in, we already have it.
*           ------------------------------------------------------------
*           if (jSr .ne. jq) then
            call s2unpk
     &         ( jBq, m, n, ne, Anorm, nlocA, locA, indA, Acol, y1 )
            call s2Bsol
     &         ( inform, WithL, m, y1, y, iw, leniw, rw, lenrw )
            iw(modlu) = 1       ! Update the LU
*           end if

         else
*           ============================================================
*           A variable in S hit a bound.
*           ============================================================
            hs(jr)    = jrStat
            jSr       = jr
            kBSq      = kp
            kSq       = kBSq - m
            iw(modlu) = 0       ! No LU update
         end if

         if (feasbl) then
            call s5rg
     &         ( m, nBS, n, nS, eps0,
     &           ne, nlocA, locA, indA, Acol,
     &           gBS, pi, rg, rgNorm, kBS )
         end if

*        ===============================================================
*        If necessary, swap the largest reduced-gradient in  Z2  into
*        the front of  Z2,  so it will end up at the end of  Z1.
*        ===============================================================
         rgdel  = abs( rg(kSq) )
         if (parHes  .and.  kSq .le. maxR) then
            call s5Sswp
     &         ( m, maxR, lenR, nS, nBS,
     &           kBS, blBS, buBS, gBS, R, rg, xBS )
         end if

*        ---------------------------------------------------------------
*        Delete the kSq-th superbasic, updating R if it exists and
*        adjusting all arrays in BS order.
*        ---------------------------------------------------------------
         if (gotR  .and.  kSq .lt. nS) then
            call s6Rdel
     &         ( kSq, maxR, nS, lenR, R, eps )
         end if

         call s5Sdel
     &      ( kSq, m, nS, nBS, kBS, blBS, buBS, gBS, rg, xBS )

         nS     = nS  - 1
         nBS    = m   + nS

         if (rgNorm .le. rgdel) then
            rgNorm = zero
            if (nS .gt. 0) then
               rgNorm = dnormi( nS, rg, 1 )
            end if
         end if
      end if ! if uncon

  900 return

 1000 format(' Itn', i7, ': CG gives  gp = ', 1p, e8.1)
 9999 format(' Itn', i7, ': Chzq failed in s5QNit!!')

      end ! subroutine s5QNit

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5Sswp
     &   ( m, maxR, lenR, nS, nBS,
     &     kBS, blBS, buBS, gBS, R, rg, xBS )

      implicit
     &     none
      integer
     &     m, maxR, lenR, nS, nBS, kBS(nBS)
      double precision
     &     blBS(nBS), buBS(nBS), gBS(nBS), R(lenR), rg(nS), xBS(nBS)

*     ==================================================================
*     s5Sswp  (superbasic swap)  finds the largest reduced gradient in
*     the range  rg(maxR+1), ..., rg(nS)  and swaps it into position
*     maxR + 1  (so we know where it is).
*
*     16 Jun 2001: First version based on MINOS routine m6swap.
*     16 Jun 2001: Current version of s5Sswp.
*     ==================================================================
      external
     &     idamax
      integer
     &     idamax, j, j1, k, k1, k2, lastR, ldiag1, ldiag2, nz2
      double precision
     &     bl1, bu1, gBS1, rdiag1, rg1, xBS1
*     ------------------------------------------------------------------
      k1  = maxR + 1

      if (nS .gt. k1) then
         nz2    = nS - maxR
         k2     = maxR + idamax( nz2, rg(k1), 1 )
         if (k2 .gt. k1) then
            j      = m + k1
            k      = m + k2
            lastR  = maxR*k1/2
            ldiag1 = lastR + 1
            ldiag2 = lastR + (k2 - maxR)

            rdiag1    = R(ldiag1)
            rg1       = rg(k1)
            j1        = kBS(j)
            bl1       = blBS(j)
            bu1       = buBS(j)
            gBS1      = gBS(j)
            xBS1      = xBS(j)

            R(ldiag1) = R(ldiag2)
            rg(k1)    = rg(k2)
            kBS(j)    = kBS(k)
            blBS(j)   = blBS(k)
            buBS(j)   = buBS(k)
            gBS(j)    = gBS(k)
            xBS(j)    = xBS(k)

            R(ldiag2) = rdiag1
            rg(k2)    = rg1
            kBS(k)    = j1
            blBS(k)   = bl1
            buBS(k)   = bu1
            gBS(k)    = gBS1
            xBS(k)    = xBS1
         end if
      end if

      end ! subroutine s5Sswp

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5Hzx
     &   ( nS, x, ZHZx,
     &     Hprod, Hprod1,
     &     ne, nlocA, locA, indA, Acol,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1
      integer
     &     lencu, lencw, leniu, leniw, lenru, lenrw, ne, nlocA, nS,
     &     locA(nlocA), indA(ne), iu(leniu), iw(leniw)
      double precision
     &     Acol(ne), x(nS), ZHZx(nS), ru(lenru), rw(lenrw)
      character
     &    cu(lencu)*8,  cw(lencw)*8

*     ==================================================================
*     s5Hzx  is called by SYMMLQ to compute the vector
*       ZHZx = Z'HZ x  or  R^(-T) Z'HZ R^(-1) x.
*     (Preconditioning with R'R is done explicitly,
*     not via SYMMLQ's Msolve.)
*
*     16 Nov 2001: First version of s5Hzx.
*     17 Aug 2004: Current version.
*     ==================================================================
      integer
     &     Hcalls, lenR, lkBS, lR, ly, ly4, mBS, m, maxR, maxS,
     &     n, nb, nnH
*     ------------------------------------------------------------------
      parameter         (Hcalls    = 188) ! number of Hx products
*     ------------------------------------------------------------------
      lenR      = iw( 28) ! R(lenR) is the reduced Hessian factor
      n         = iw( 15) ! copy of the number of columns
      m         = iw( 16) ! copy of the number of rows
      nnH       = iw( 24) !    max( nnObj, nnJac )
      maxR      = iw( 52) ! max columns of R
      maxS      = iw( 53) ! max # of superbasics
      lkBS      = iw(292) ! kBS(mBS)    = ( B  S ) list

      lR        = iw(294) ! R(lenR)     = factor of Z'HZ
      ly        = iw(311) ! y (nb)      =  real work vector
      ly4       = iw(315) ! y4(nb)      =  real work vector

      nb    = n + m
      mBS   = m + maxS

      call s5Hzx1
     &   ( Hprod, Hprod1,
     &     iw(Hcalls), nnH, m, mBS, n, nb, nS, maxR,
     &     ne, nlocA, locA, indA, Acol,
     &     iw(lkBS), x, ZHZx, rw(ly), rw(ly4), lenR, rw(lR),
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine s5Hzx

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5Hzx1
     &   ( Hprod, Hprod1,
     &     Hcalls, nnH, m, mBS, n, nb, nS, maxR,
     &     ne, nlocA, locA, indA, Acol,
     &     kBS, x, ZHZx, y, y1, lenR, R,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1
      integer
     &     Hcalls, lencu, lencw, leniu, leniw, lenru, lenrw, lenR,
     &     maxR, mBS, n, nb, nnH, nS, ne, nlocA, indA(ne),
     &     iu(leniu), iw(leniw), locA(nlocA), kBS(mBS)
      double precision
     &     Acol(ne), R(lenR), y(nb), y1(nb), x(nS), ZHZx(nS),
     &     ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     s5Hzx1  does the work for s5Hzx.
*     Computes the vector ZHZx = Z'HZ * x  or  R^(-T) Z'HZ R^(-1) x.
*
*     16 Nov 2001: First version of s5Hzx1.
*     17 Aug 2004: Current version.
*     ==================================================================
      integer
     &     inform, m, m1, minimz, nBS, PreCon, Status
      double precision
     &     eps0, sgnObj
*     ------------------------------------------------------------------
      integer            Normal,     Transp
      parameter         (Normal = 0, Transp = 1)
      integer            WithR,      WithRt
      parameter         (WithR  = 0, WithRt = 1)
      integer            WithB,      WithBt
      parameter         (WithB  = 1, WithBt = 2)

      double precision   zero,           one
      parameter         (zero = 0.0d+0,  one = 1.0d+0)
*     ------------------------------------------------------------------
      eps0      = rw(  2) ! eps**(4/5)       IEEE DP  3.00e-13
      minimz    = iw(199) ! (-1)(+1)    => (max)(min)
      PreCon    = iw(209) ! Current precon mode (based on QPslvr)

      nBS       = m + nS
      m1        = m + 1

      sgnObj    = minimz

*     ------------------------------------------------------------------
*     Set y = Z*x.
*     ------------------------------------------------------------------
      call dcopy
     &   ( nS, x, 1, y(m1), 1 )

      if (PreCon .eq. 1) then
         call s6Rsol
     &      ( WithR, maxR, nS, lenR, R, y(m1) )
      end if

*     Compute  y1 = - S*yS.

      call s2Bprd
     &   ( Normal, eps0, n, nS, kBS(m1),
     &     ne, nlocA, locA, indA, Acol,
     &     (-one), y(m1), nS, zero, y1, m )

*     Solve  B*y = y1.

      call s2Bsol
     &   ( inform, WithB, m, y1, y, iw, leniw, rw, lenrw  )

*     ------------------------------------------------------------------
*     Set  y = H*y = H*Z*x.
*     ------------------------------------------------------------------
      call s2scatr
     &   ( nnH, nBS, kBS, one, y, y1 )
*     call s8Hx
*    &   ( Hcalls, nnH, y1, y, cw, lencw, iw, leniw, rw, lenrw )

      Status = 0
      call Hprod
     &   ( Hprod1, Hcalls, nnH,
     &     y1, y, Status,       ! y1 = dx,  y = Hdx
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     Compute  ZHZx = Z'y.
*     ------------------------------------------------------------------
*     Gather y1 = yBS and solve  B'y = y1B.

      call s2gathr
     &   ( nnH, nBS, kBS, sgnObj, y, y1 )
      call s2Bsol
     &   ( inform, WithBt, m, y1, y, iw, leniw, rw, lenrw )

*     Set ZHZx = y1S - S'y.

      call dcopy
     &   ( nS, y1(m1), 1, ZHZx, 1 )
      call s2Bprd
     &   ( Transp, eps0, n, nS, kBS(m1),
     &     ne, nlocA, locA, indA, Acol,
     &     (-one), y, m, one, ZHZx, nS )

      if (PreCon .eq. 1) then
         call s6Rsol
     &      ( WithRt, maxR, nS, lenR, R, ZHZx )
      end if

      end ! subroutine s5Hz1

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5Zswp
     &   ( statpt, m, maxR, maxS, nBS, lenR, nS,
     &     tolrg, kBS, blBS, buBS, gBS, R, rg, xBS, rw, lenrw )

      implicit
     &     none
      logical
     &     statpt
      integer
     &     lenrw, lenR, m, maxR, maxS, nBS, nS,
     &     kBS(nBS)
      double precision
     &     tolrg, blBS(nBS), buBS(nBS),
     &     gBS(nBS), R(lenR), rg(maxS), xBS(nBS), rw(lenrw)

*     ==================================================================
*     s5Zswp  (subspace convergence)  decides whether or not
*     optimization should continue on the current subspace.
*     On exit,  statpt = .false.  means it should,
*               statpt = .true    means it should not.
*
*     R  is a partial Hessian for the
*     first  maxR  superbasic variables.  The superbasics are then
*     in two sets  Z1  (containing   nZ1 = maxR       variables)
*             and  Z2  (containing   nZ2 = nS - maxR  variables).
*
*     The null-space matrix is similarly partitioned as  Z = ( Z1  Z2 ).
*
*     The normal convergence test is first applied to  Z1.  If it looks
*     like we are near enough to an optimum in that restricted
*     subspace, we find the largest reduced gradient in  Z2  (index k2).
*     If this is less than  tolrg  we exit with  nxtphs = 3  (to ask for
*     price).  Otherwise we make room in  Z1  for the corresponding
*     variable by the moving the superbasic with the smallest reduced
*     gradient in  Z1  (index  k1)  to the end of  Z2.
*
*     16 Jun 2001: First version based on MINOS routine m7sscv.
*     19 Jul 2001: Current version of s5Sswp.
*     ==================================================================
      integer
     &     jZ1, k, k1, k2, lastR, ldiag1, lR
      double precision
     &     blBS1, buBS1, eps, gBS1, Rdiag1, rg1,
     &     rgmin1, rgNrm2, xBS1
*     ------------------------------------------------------------------
      eps       = rw(  1) ! unit round-off.  IEEE DP  2.22e-16

*     Swap the largest reduced gradient in  Z2  to the front of  Z2
*     and see if it is significantly large.

      call s5Sswp
     &   ( m, maxR, lenR, nS, nBS,
     &     kBS, blBS, buBS, gBS, R, rg, xBS )

      k2     = maxR + 1
      rgNrm2 = abs( rg(k2) )

      if (rgNrm2 .gt. tolrg) then

*        Find the smallest component of  Z1'g.

         rgmin1 = abs( rg(1) )
         k1     = 1
         do  k  = 1, maxR
            if (rgmin1 .ge. abs( rg(k) )) then
               rgmin1 = abs( rg(k) )
               k1     = k
            end if
         end do

         if (rgmin1 .lt. rgNrm2) then

*           Save the relevant values.

            statpt = .false.
            lastR  = maxR*(maxR + 1)/2
            ldiag1 = (k1 - 1)*maxR + (3 - k1)*k1/2  ! Magic formula!
            Rdiag1 = R(ldiag1)
            rg1    = rg(k1)
            k      = m + k1
            jZ1    = kBS(k)
            blBS1  = blBS(k)
            buBS1  = buBS(k)
            gBS1   = gBS(k)
            xBS1   = xBS(k)

*           Delete the k1-th variable from  Z1,  and shift the remaining
*           superbasics in  Z1  and  Z2  one place to the left.

            call s6rdel
     &         ( k1, maxR, nS, lenR, R, eps )
            call s5Sdel
     &         ( k1, m, nS, nBS, kBS, blBS, buBS, gBS, rg, xBS )

*           Put the old k1-th superbasic in at the very end.

            lR      = lastR + (nS - maxR)
            R(lR)   = Rdiag1
            rg(nS)  = rg1
            kBS(nBS)  = jZ1
            blBS(nBS) = blBS1
            buBS(nBS) = buBS1
            gBS(nBS)  = gBS1
            xBS(nBS)  = xBS1
         end if
      end if

      end ! subroutine s5Zswp

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine SYMMLQ
     &   ( n, b, r1, r2, v, w, x, y,
     &     Aprod, Msolve, checkA, goodb, precon, shift,
     &     nout , itnlim, rtol,
     &     iStop, itn, Anorm, Acond, rnorm, ynorm,
     &     Hprod, Hprod1,                              ! Added for SNOPT
     &     ne, nlocA, locA, indA, Acol,               ! Added for SNOPT
     &     cu, lencu, iu, leniu, ru, lenru,            ! Added for SNOPT
     &     cw, lencw, iw, leniw, rw, lenrw )           ! Added for SNOPT

      implicit           none
      external           Aprod, Msolve
      integer            n, nout, itnlim, iStop, itn
      logical            checkA, goodb, precon
      double precision   shift, rtol, Anorm, Acond, rnorm, ynorm,
     &                   b(n), r1(n), r2(n), v(n), w(n), x(n), y(n)

      external                                         ! Added for SNOPT
     &     Hprod, Hprod1                               ! Added for SNOPT
      integer                                          ! Added for SNOPT
     &     ne, nlocA, lencu, lencw,                    ! Added for SNOPT
     &     leniu, leniw, lenru, lenrw,                 ! Added for SNOPT
     &     locA(nlocA), indA(ne), iu(leniu), iw(leniw) ! Added for SNOPT
      double precision                                 ! Added for SNOPT
     &     Acol(ne), ru(lenru), rw(lenrw)              ! Added for SNOPT
      character                                        ! Added for SNOPT
     &     cu(lencu)*8, cw(lencw)*8                    ! Added for SNOPT
*     ------------------------------------------------------------------
*
*     SYMMLQ  is designed to solve the system of linear equations
*
*                Ax = b
*
*     where A is an n by n symmetric matrix and b is a given vector.
*     The matrix A is not required to be positive definite.
*     (If A is known to be definite, the method of conjugate gradients
*     might be preferred, since it will require about the same number of
*     iterations as SYMMLQ but slightly less work per iteration.)
*
*
*     The matrix A is intended to be large and sparse.  It is accessed
*     by means of a subroutine call of the form
*
*                call Aprod ( n, x, y )
*
*     which must return the product y = Ax for any given vector x.
*
*
*     More generally, SYMMLQ is designed to solve the system
*
*                (A - shift*I) x = b
*
*     where  shift  is a specified scalar value.  If  shift  and  b
*     are suitably chosen, the computed vector x may approximate an
*     (unnormalized) eigenvector of A, as in the methods of
*     inverse iteration and/or Rayleigh-quotient iteration.
*     Again, the matrix (A - shift*I) need not be positive definite.
*     The work per iteration is very slightly less if  shift = 0.
*
*
*     A further option is that of preconditioning, which may reduce
*     the number of iterations required.  If M = C C' is a positive
*     definite matrix that is known to approximate  (A - shift*I)
*     in some sense, and if systems of the form  My = x  can be
*     solved efficiently, the parameters precon and Msolve may be
*     used (see below).  When  precon = .true., SYMMLQ will
*     implicitly solve the system of equations
*
*             P (A - shift*I) P' xbar  =  P b,
*
*     i.e.                  Abar xbar  =  bbar
*     where                         P  =  C**(-1),
*                                Abar  =  P (A - shift*I) P',
*                                bbar  =  P b,
*
*     and return the solution       x  =  P' xbar.
*     The associated residual is rbar  =  bbar - Abar xbar
*                                      =  P (b - (A - shift*I)x)
*                                      =  P r.
*
*     In the discussion below, eps refers to the machine precision.
*     eps is computed by SYMMLQ.  A typical value is eps = 2.22e-16
*     for IBM mainframes and IEEE double-precision arithmetic.
*
*     Parameters
*     ----------
*
*     n       input      The dimension of the matrix A.
*
*     b(n)    input      The rhs vector b.
*
*     r1(n)   workspace
*     r2(n)   workspace
*     v(n)    workspace
*     w(n)    workspace
*
*     x(n)    output     Returns the computed solution  x.
*
*     y(n)    workspace
*
*     Aprod   external   A subroutine defining the matrix A.
*                        For a given vector x, the statement
*
*                              call Aprod ( n, x, y )
*
*                        must return the product y = Ax
*                        without altering the vector x.
*
*     Msolve  external   An optional subroutine defining a
*                        preconditioning matrix M, which should
*                        approximate (A - shift*I) in some sense.
*                        M must be positive definite.
*                        For a given vector x, the statement
*
*                              call Msolve( n, x, y )
*
*                        must solve the linear system My = x
*                        without altering the vector x.
*
*                        In general, M should be chosen so that Abar has
*                        clustered eigenvalues.  For example,
*                        if A is positive definite, Abar would ideally
*                        be close to a multiple of I.
*                        If A or A - shift*I is indefinite, Abar might
*                        be close to a multiple of diag( I  -I ).
*
*                        NOTE.  The program calling SYMMLQ must declare
*                        Aprod and Msolve to be external.
*
*     checkA  input      If checkA = .true., an extra call of Aprod will
*                        be used to check if A is symmetric.  Also,
*                        if precon = .true., an extra call of Msolve
*                        will be used to check if M is symmetric.
*
*     goodb   input      Usually, goodb should be .false.
*                        If x is expected to contain a large multiple of
*                        b (as in Rayleigh-quotient iteration),
*                        better precision may result if goodb = .true.
*                        See Lewis (1977) below.
*                        When goodb = .true., an extra call to Msolve
*                        is required.
*
*     precon  input      If precon = .true., preconditioning will
*                        be invoked.  Otherwise, subroutine Msolve
*                        will not be referenced; in this case the
*                        actual parameter corresponding to Msolve may
*                        be the same as that corresponding to Aprod.
*
*     shift   input      Should be zero if the system Ax = b is to be
*                        solved.  Otherwise, it could be an
*                        approximation to an eigenvalue of A, such as
*                        the Rayleigh quotient b'Ab / (b'b)
*                        corresponding to the vector b.
*                        If b is sufficiently like an eigenvector
*                        corresponding to an eigenvalue near shift,
*                        then the computed x may have very large
*                        components.  When normalized, x may be
*                        closer to an eigenvector than b.
*
*     nout    input      A file number.
*                        If nout .gt. 0, a summary of the iterations
*                        will be printed on unit nout.
*
*     itnlim  input      An upper limit on the number of iterations.
*
*     rtol    input      A user-specified tolerance.  SYMMLQ terminates
*                        if it appears that norm(rbar) is smaller than
*                              rtol * norm(Abar) * norm(xbar),
*                        where rbar is the transformed residual vector,
*                              rbar = bbar - Abar xbar.
*
*                        If shift = 0 and precon = .false., SYMMLQ
*                        terminates if norm(b - A*x) is smaller than
*                              rtol * norm(A) * norm(x).
*
*     iStop   output     An integer giving the reason for termination...
*
*              -1        beta2 = 0 in the Lanczos iteration; i.e. the
*                        second Lanczos vector is zero.  This means the
*                        rhs is very special.
*                        If there is no preconditioner, b is an
*                        eigenvector of A.
*                        Otherwise (if precon is true), let My = b.
*                        If shift is zero, y is a solution of the
*                        generalized eigenvalue problem Ay = lambda My,
*                        with lambda = alpha1 from the Lanczos vectors.
*
*                        In general, (A - shift*I)x = b
*                        has the solution         x = (1/alpha1) y
*                        where My = b.
*
*               0        b = 0, so the exact solution is x = 0.
*                        No iterations were performed.
*
*               1        Norm(rbar) appears to be less than
*                        the value  rtol * norm(Abar) * norm(xbar).
*                        The solution in  x  should be acceptable.
*
*               2        Norm(rbar) appears to be less than
*                        the value  eps * norm(Abar) * norm(xbar).
*                        This means that the residual is as small as
*                        seems reasonable on this machine.
*
*               3        Norm(Abar) * norm(xbar) exceeds norm(b)/eps,
*                        which should indicate that x has essentially
*                        converged to an eigenvector of A
*                        corresponding to the eigenvalue shift.
*
*               4        Acond (see below) has exceeded 0.1/eps, so
*                        the matrix Abar must be very ill-conditioned.
*                        x may not contain an acceptable solution.
*
*               5        The iteration limit was reached before any of
*                        the previous criteria were satisfied.
*
*               6        The matrix defined by Aprod does not appear
*                        to be symmetric.
*                        For certain vectors y = Av and r = Ay, the
*                        products y'y and r'v differ significantly.
*
*               7        The matrix defined by Msolve does not appear
*                        to be symmetric.
*                        For vectors satisfying My = v and Mr = y, the
*                        products y'y and r'v differ significantly.
*
*               8        An inner product of the form  x' M**(-1) x
*                        was not positive, so the preconditioning matrix
*                        M does not appear to be positive definite.
*
*                        If iStop .ge. 5, the final x may not be an
*                        acceptable solution.
*
*     itn     output     The number of iterations performed.
*
*     Anorm   output     An estimate of the norm of the matrix operator
*                        Abar = P (A - shift*I) P',   where P = C**(-1).
*
*     Acond   output     An estimate of the condition of Abar above.
*                        This will usually be a substantial
*                        under-estimate of the true condition.
*
*     rnorm   output     An estimate of the norm of the final
*                        transformed residual vector,
*                           P (b  -  (A - shift*I) x).
*
*     ynorm   output     An estimate of the norm of xbar.
*                        This is sqrt( x'Mx ).  If precon is false,
*                        ynorm is an estimate of norm(x).
*
*
*
*     To change precision
*     -------------------
*
*     Alter the words
*            double precision,
*            daxpy, dcopy, ddot, dnrm2
*     to their single or double equivalents.
*     ------------------------------------------------------------------
*
*
*     This routine is an implementation of the algorithm described in
*     the following references:
*
*     C.C. Paige and M.A. Saunders,  Solution of Sparse Indefinite
*          Systems of Linear Equations,
*          SIAM J. Numer. Anal. 12, 4, September 1975, pp. 617-629.
*
*     J.G. Lewis,  Algorithms for Sparse Matrix Eigenvalue Problems,
*          Report STAN-CS-77-595, Computer Science Department,
*          Stanford University, Stanford, California, March 1977.
*
*     Applications of SYMMLQ and the theory of preconditioning
*     are described in the following references:
*
*     D.B. Szyld and O.B. Widlund,  Applications of Conjugate Gradient
*          Type Methods to Eigenvalue Calculations,
*          in R. Vichnevetsky and R.S. Steplman (editors),
*          Advances in Computer Methods for Partial Differential
*          Equations -- III, IMACS, 1979, 167-173.
*
*     D.B. Szyld,  A Two-level Iterative Method for Large Sparse
*          Generalized Eigenvalue Calculations,
*          Ph. D. dissertation, Department of Mathematics,
*          New York University, New York, October 1983.
*
*     P.E. Gill, W. Murray, D.B. Ponceleon and M.A. Saunders,
*          Preconditioners for indefinite systems arising in
*          optimization, SIMAX 13, 1, 292--311, January 1992.
*          (SIAM J. on Matrix Analysis and Applications)
*     ------------------------------------------------------------------
*
*
*     SYMMLQ development:
*            1972: First version.
*            1975: John Lewis recommended modifications to help with
*                  inverse iteration:
*                  1. Reorthogonalize v1 and v2.
*                  2. Regard the solution as x = x1  +  bstep * b,
*                     with x1 and bstep accumulated separately
*                     and bstep * b added at the end.
*                     (In inverse iteration, b might be close to the
*                     required x already, so x1 may be a lot smaller
*                     than the multiple of b.)
*            1978: Daniel Szyld and Olof Widlund implemented the first
*                  form of preconditioning.
*                  This required both a solve and a multiply with M.
*            1979: Implemented present method for preconditioning.
*                  This requires only a solve with M.
*            1984: Sven Hammarling noted corrections to tnorm and x1lq.
*                  SYMMLQ added to NAG Fortran Library.
*     15 Sep 1985: Final F66 version.  SYMMLQ sent to "misc" in netlib.
*     16 Feb 1989: First F77 version.
*
*     22 Feb 1989: Hans Mittelmann observed beta2 = 0 (hence failure)
*                  if Abar = const*I.  iStop = -1 added for this case.
*
*     01 Mar 1989: Hans Mittelmann observed premature termination on
*                  ( 1  1  1 )     (   )                   ( 1  1    )
*                  ( 1  1    ) x = ( 1 ),  for which  T3 = ( 1  1  1 ).
*                  ( 1     1 )     (   )                   (    1  1 )
*                  T2 is exactly singular, so estimating cond(A) from
*                  the diagonals of Lbar is unsafe.  We now use
*                  L       or  Lbar         depending on whether
*                  lqnorm  or  cgnorm       is least.
*
*     03 Mar 1989: eps computed internally instead of coming in as a
*                  parameter.
*     07 Jun 1989: ncheck added as a parameter to say if A and M
*                  should be checked for symmetry.
*                  Later changed to checkA (see below).
*     20 Nov 1990: goodb added as a parameter to make Lewis's changes
*                  an option.  Usually b is NOT much like x.  Setting
*                  goodb = .false. saves a call to Msolve at the end.
*     20 Nov 1990: Residual not computed exactly at end, to save time
*                  when only one or two iterations are required
*                  (e.g. if the preconditioner is very good).
*                  Beware, if precon is true, rnorm estimates the
*                  residual of the preconditioned system, not Ax = b.
*     04 Sep 1991: Parameter list changed and reordered.
*                  integer ncheck is now logical checkA.
*     22 Jul 1992: Example from Lothar Reichel and Daniela Calvetti
*                  showed that beta2 = 0 (iStop = -1) means that
*                  b is an eigenvector when M = I.
*                  More complicated if there is a preconditioner;
*                  not clear yet how to describe it.
*     20 Oct 1999: Bug.  alfa1 = 0 caused Anorm = 0, divide by zero.
*                  Need to estimate Anorm from column of Tk.
*     18 Nov 2001: bnorm printed in first line, since rnorm at itn = 0
*                  isn't bnorm as we might expect -- it's cgnorm after
*                  the proverbial step to the CG point!
*     02 Aug 2003: eps grabbed from rw (special to SNOPT).
*     23 Dec 2003: Normal backward error exit (iStop = 1) terminates
*                  too early when rtol = 0.01 say (rather large).
*                  Keep it with   rtol = eps  (iStop = 2)
*                  but replace normal test with ||r|| < rtol*||b||
*                  like most other CG solvers.  This fits better with
*                  the inexact Newton context.  Note that it is correct
*                  only with precon = .false.
*
*     08 Apr 2004: Always move to the CG point.
*
*     Michael A. Saunders                   na.msaunders@na-net.ornl.gov
*     Systems Optimization Laboratory       saunders@stanford.edu
*     Dept of Management Science and Engineering
*     Stanford University
*     Stanford, CA 94305-4026                             (650) 723-1875
*     ------------------------------------------------------------------
*
*
*     Subroutines and functions
*
*     USER       Aprod, Msolve
*     BLAS       daxpy, dcopy, ddot , dnrm2
*
*
*     Intrinsics and local variables

      intrinsic          abs, max, min, mod, sqrt
      double precision   ddot, dnrm2
      double precision   alfa, b1, beta, beta1, bnorm, bstep, cs,
     &                   cgnorm, dbar, delta, denom, diag,
     &                   eps, epsa, epsln, epsr, epsx,
     &                   gamma, gbar, gmax, gmin,
     &                   lqnorm, oldb, qrnorm, rhs1, rhs2,
     &                   s, sn, snprod, t, tnorm,
     &                   x1cg, ynorm2, zbar, z
      integer            i

      double precision   zero         ,  one
      parameter        ( zero = 0.0d+0,  one = 1.0d+0)

      character          enter*16, exit*16, msg(-1:9)*52

      data               enter /' Enter SYMMLQ.  '/,
     &                   exit  /' Exit  SYMMLQ.  '/

      data               msg
     & / 'beta2 = 0.  If M = I, b and x are eigenvectors of A',
     &   'beta1 = 0.  The exact solution is  x = 0',
     &   'Requested accuracy achieved, as determined by rtol',
     &   'Reasonable accuracy achieved, given eps',
     &   'x has converged to an eigenvector',
     &   'Acond has exceeded 0.1/eps',
     &   'The iteration limit was reached',
     &   'Aprod  does not define a symmetric matrix',
     &   'Msolve does not define a symmetric matrix',
     &   'Msolve does not define a pos-def preconditioner',
     &   'x is not a descent direction' /
*     ------------------------------------------------------------------

      eps       = rw(  1) ! unit round-off.  IEEE DP  2.22e-16

*     Print heading and initialize.

      bnorm  = dnrm2 ( n, b, 1 )
      beta1  = bnorm
      if (nout .gt. 0) then
         write(nout, 1000) enter, beta1,
     &                     n, checkA, goodb, precon,
     &                     itnlim, rtol, shift
      end if
      iStop  = 0
      itn    = 0
      Anorm  = zero
      Acond  = zero
      rnorm  = zero
      ynorm  = zero

      call dload ( n, zero, x, 1 )

*     Set up y for the first Lanczos vector v1.
*     y is really beta1 * P * v1  where  P = C**(-1).
*     y and beta1 will be zero if b = 0.

      call dcopy ( n, b, 1, y , 1 )
      call dcopy ( n, b, 1, r1, 1 )
      if ( precon ) call Msolve( n, r1, y )
      if ( goodb  ) then
         b1  = y(1)
      else
         b1  = zero
      end if
      beta1  = ddot  ( n, r1, 1, y, 1 )

*     See if Msolve is symmetric.

      if (checkA  .and.  precon) then
         call Msolve( n, y, r2 )
         s      = ddot  ( n, y, 1, y, 1 )
         t      = ddot  ( n,r1, 1,r2, 1 )
         z      = abs( s - t )
         epsa   = (s + eps) * eps**0.33333D+0
         if (z .gt. epsa) then
            iStop = 7
            go to 900
         end if
      end if

*     Test for an indefinite preconditioner.

      if (beta1 .lt. zero) then
         iStop = 8
         go to 900
      end if

*     If b = 0 exactly, stop with x = 0.

      if (beta1 .eq. zero) then
         go to 900
      end if

*     Here and later, v is really P * (the Lanczos v).

      beta1  = sqrt( beta1 )
      s      = one / beta1
      call dcopy ( n,    y, 1, v, 1 )
      call dscal ( n, s, v, 1 )

      call Aprod
     &   ( n, v, y,
     &     Hprod, Hprod1,                              ! Added for SNOPT
     &     ne , nlocA, locA, indA, Acol,               ! Added for SNOPT
     &     cu, lencu, iu, leniu, ru, lenru,            ! Added for SNOPT
     &     cw, lencw, iw, leniw, rw, lenrw )           ! Added for SNOPT

*     See if Aprod  is symmetric.

      if (checkA) then
         call Aprod
     &      ( n, y, r2,
     &        Hprod, Hprod1,                           ! Added for SNOPT
     &        ne , nlocA, locA, indA, Acol,            ! Added for SNOPT
     &        cu, lencu, iu, leniu, ru, lenru,         ! Added for SNOPT
     &        cw, lencw, iw, leniw, rw, lenrw )        ! Added for SNOPT
         s      = ddot  ( n, y, 1, y, 1 )
         t      = ddot  ( n, v, 1,r2, 1 )
         z      = abs( s - t )
         epsa   = (s + eps) * eps**0.33333D+0
         if (z .gt. epsa) then
            iStop = 6
            go to 900
         end if
      end if

*     Set up y for the second Lanczos vector.
*     Again, y is beta * P * v2  where  P = C**(-1).
*     y and beta will be zero or very small if b is an eigenvector.

      call daxpy ( n, (- shift), v, 1, y, 1 )
      alfa   = ddot  ( n, v, 1, y, 1 )
      call daxpy ( n, (- alfa / beta1), r1, 1, y, 1 )

*     Make sure  r2  will be orthogonal to the first  v.

      z      = ddot  ( n, v, 1, y, 1 )
      s      = ddot  ( n, v, 1, v, 1 )
      call daxpy ( n, (- z / s), v, 1, y, 1 )

      call dcopy ( n, y, 1, r2, 1 )
      if ( precon ) call Msolve( n, r2, y )
      oldb   = beta1
      beta   = ddot  ( n, r2, 1, y, 1 )

      if (beta .lt. zero) then
         iStop = 8
         go to 900
      end if

*     Cause termination (later) if beta is essentially zero.

      beta   = sqrt( beta )
      if (beta .le. eps) then
         iStop = -1
      end if

*     See if the local reorthogonalization achieved anything.

      denom  = sqrt( s ) * dnrm2( n, r2, 1 )  +  eps
      s      = z / denom
      t      = ddot  ( n, v, 1, r2, 1 ) / denom
      if (nout .gt. 0  .and.  goodb) then
         write(nout, 1100) beta1, alfa, s, t
      end if

*     Initialize other quantities.

      cgnorm = beta1
      gbar   = alfa
      dbar   = beta
      rhs1   = beta1
      rhs2   = zero
      bstep  = zero
      snprod = one
      tnorm  = alfa**2 + beta**2
      ynorm2 = zero
      gmax   = abs( alfa ) + eps
      gmin   = gmax

      if ( goodb ) then
         call dload ( n, zero, w, 1 )
      else
         call dcopy ( n, v, 1, w, 1 )
      end if

*     ------------------------------------------------------------------
*     Main iteration loop.
*     ------------------------------------------------------------------
!+    Repeat                                          Until iStop .ne. 0

*        Estimate various norms and test for convergence.

  100    Anorm  = sqrt( tnorm  )
         ynorm  = sqrt( ynorm2 )
         epsa   = Anorm * eps
         epsx   = Anorm * ynorm * eps
         epsr   = Anorm * ynorm * rtol
         diag   = gbar
         if (diag .eq. zero) diag = epsa

         lqnorm = sqrt( rhs1**2 + rhs2**2 )
         qrnorm = snprod * beta1
         cgnorm = qrnorm * beta / abs( diag )

*        Estimate  cond(A).
*        In this version we look at the diagonals of  L  in the
*        factorization of the tridiagonal matrix,  T = L*Q.
*        Sometimes, T(k) can be misleadingly ill-conditioned when
*        T(k+1) is not, so we must be careful not to overestimate Acond.

         if (lqnorm .le. cgnorm) then
            Acond  = gmax / gmin
         else
            denom  = min( gmin, abs( diag ) )
            Acond  = gmax / denom
         end if

*        See if any of the stopping criteria are satisfied.
*        In rare cases, iStop is already -1 from above (Abar = const * I).

         if (iStop .eq. 0) then
            if (itn    .ge. itnlim    ) iStop = 5
            if (Acond  .ge. 0.1d+0/eps) iStop = 4
            if (epsx   .ge. beta1     ) iStop = 3
            if (cgnorm .le. epsx      ) iStop = 2
            if (cgnorm .le. bnorm*rtol) iStop = 1 ! Replaces next line
*           if (cgnorm .le. epsr      ) iStop = 1 ! for inexact Newton.
         end if

*        ===============================================================
*        See if it is time to print something.

         if ((n      .le. 40             .or.
     &        itn    .le. 10             .or.
     &        itn    .ge. itnlim - 10    .or.
     &        mod(itn,10) .eq.   0       .or.
     &        cgnorm .le. 10.0d+0*epsx   .or.
     &        cgnorm .le. 10.0d+0*epsr   .or.
     &        Acond  .ge. 0.01d+0/eps    .or.
     &        iStop  .ne. 0                  ) .and.
     &        nout   .gt. 0                  ) then

*           Print a line for this iteration.

            zbar   = rhs1 / diag
            z      = (snprod * zbar  +  bstep) / beta1
*           x1lq   = x(1)  +  b1 * bstep / beta1
            x1cg   = x(1)  +  w(1) * zbar  +  b1 * z

            if (    itn     .eq. 0) write(nout, 1200)
            write(nout, 1300) itn, x1cg, cgnorm, bstep/beta1,Anorm,Acond
            if (mod(itn,10) .eq. 0) write(nout, 1500)
         end if
*        ===============================================================

*        Obtain the current Lanczos vector  v = (1 / beta)*y
*        and set up  y  for the next iteration.

         if (iStop .eq. 0) then
            s      = one / beta

            call dcopy ( n,    y, 1, v, 1 )
            call dscal ( n, s, v, 1 )

            call Aprod
     &         ( n, v, y,
     &           Hprod, Hprod1,                        ! Added for SNOPT
     &           ne , nlocA, locA, indA, Acol,         ! Added for SNOPT
     &           cu, lencu, iu, leniu, ru, lenru,      ! Added for SNOPT
     &           cw, lencw, iw, leniw, rw, lenrw )     ! Added for SNOPT
            call daxpy ( n, (- shift), v, 1, y, 1 )
            call daxpy ( n, (- beta / oldb), r1, 1, y, 1 )
            alfa   = ddot( n, v, 1, y, 1 )
            call daxpy ( n, (- alfa / beta), r2, 1, y, 1 )
            call dcopy ( n, r2, 1, r1, 1 )
            call dcopy ( n, y, 1, r2, 1 )
            if ( precon ) call Msolve( n, r2, y )
            oldb   = beta
            beta   = ddot  ( n, r2, 1, y, 1 )
            if (beta .lt. zero) then
               iStop = 6
               go to 800
            end if
            beta   = sqrt( beta )
            tnorm  = tnorm  +  alfa**2  +  oldb**2  +  beta**2

*           Compute the next plane rotation for  Q.

            gamma  = sqrt( gbar**2 + oldb**2 )
            cs     = gbar / gamma
            sn     = oldb / gamma
            delta  = cs * dbar  +  sn * alfa
            gbar   = sn * dbar  -  cs * alfa
            epsln  = sn * beta
            dbar   =            -  cs * beta

*           Update  x.

            z      = rhs1 / gamma
            s      = z * cs
            t      = z * sn

            do i = 1, n
               x(i) = (w(i) * s   +   v(i) * t)  +  x(i)
               w(i) =  w(i) * sn  -   v(i) * cs
            end do

*           Accumulate the step along the direction  b,
*           and go round again.

            bstep  = snprod * cs * z  +  bstep
            snprod = snprod * sn
            gmax   = max( gmax, gamma )
            gmin   = min( gmin, gamma )
            ynorm2 = z**2  +  ynorm2
            rhs1   = rhs2  -  delta * z
            rhs2   =       -  epsln * z
            itn    = itn   +  1
         end if

!+    until    (iStop .ne. 0)
      if (.not.(iStop .ne. 0)) go to 100

*     ------------------------------------------------------------------
*     End of main iteration loop.
*     ------------------------------------------------------------------

*     Move to the CG point if it seems better.
*     In this version of SYMMLQ, the convergence tests involve
*     only cgnorm, so we're unlikely to stop at an LQ point,
*     EXCEPT if the iteration limit interferes.
*
*     April 8, 2004. Always move to the CG point for SNOPT application

  800 zbar   = rhs1 / diag
      bstep  = snprod * zbar  +  bstep
      ynorm  = sqrt( ynorm2  +  zbar**2 )
      rnorm  = cgnorm
      call daxpy ( n, zbar, w, 1, x, 1 )

      if ( goodb ) then

*        Add the step along  b.

         bstep  = bstep / beta1
         call dcopy ( n, b, 1, y, 1 )
         if ( precon ) call Msolve( n, b, y )
         call daxpy ( n, bstep, y, 1, x, 1 )
      end if

*     ==================================================================
*     Display final status.
*     ==================================================================
  900 if (nout  .gt. 0) then
         write(nout, 2000) exit, iStop, itn,
     &                     exit, Anorm, Acond,
     &                     exit, rnorm, ynorm
         write(nout, 3000) exit, msg(iStop)
      end if

      return

*     ------------------------------------------------------------------
 1000 format(// 1p,    a, 5x, 'Solution of symmetric   Ax = b', 7x,
     &          'bnorm  =', e10.2
     &       / ' n      =', i7, 5x, 'checkA =', l4, 12x,
     &          'goodb  =', l4, 7x, 'precon =', l4
     &       / ' itnlim =', i7, 5x, 'rtol   =', e11.2, 5x,
     &          'shift  =', e23.14)
 1100 format(/ 1p, ' beta1  =', e10.2, 3x, 'alpha1 =', e10.2
     &       / ' (v1,v2) before and after ', e14.2
     &       / ' local reorthogonalization', e14.2)
 1200 format(// 5x, 'itn', 7x, 'x1(cg)', 10x,
     &         'norm(r)', 5x, 'bstep', 7x, 'norm(A)', 3X, 'cond(A)')
 1300 format(1p, i8, e19.10, e11.2, e14.5, 2e10.2)
 1500 format(1x)
 2000 format(/ 1p, a, 6x, 'istop =', i3,   15x, 'itn   =', i8
     &       /     a, 6x, 'Anorm =', e12.4, 6x, 'Acond =', e12.4
     &       /     a, 6x, 'rnorm =', e12.4, 6x, 'ynorm =', e12.4)
 3000 format(      a, 6x, a )
*     ------------------------------------------------------------------
*     end of SYMMLQ
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5Msolv
     &   ( n, x, y )

      integer
     &     n
      double precision
     &     x(n), y(n)

*     ==================================================================
*     s5Msolv  dummy Msolve for SYMMLQ.
*
*     04 Dec 2004: First version of s5Msolv.
*     04 Dec 2004: Current version.
*     ==================================================================

*     Relax

      end ! subroutine s5Msolv
