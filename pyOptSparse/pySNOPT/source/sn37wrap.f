*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  sn37wrap.f --- auxiliaries for SNOPT wrappers.
*
*     SNOPTA  Arbitrary order for variables and constraints
*     s3argA   s3prtA   s3sizA   s3bldA   s3inA    s3outA   s3prmA
*     s3mapA
*
*     SNOPTB  Basic format
*     s3argB   s3prtB
*
*     SQOPT   QP wrapper
*     s3argQ   s3prtQ
*
*     NPOPT   auxiliaries for NPSOL-like calls.
*     s3argN   s3HesN   s3inN   s3outN   s3prtN
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3argA
     &   ( iExit, Start, nF, n, nS, nxname, nFname,
     &     ObjRow, neA, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     xstate, xmul, Fstate, Fmul, lvlSrt, Errors,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Errors, iExit, lvlSrt, nF, n, nxname, nFname,
     &     ObjRow, nS, neA, neG, leniw, lenrw, Start,
     &     xstate(n), Fstate(nF), iw(leniw)
      double precision
     &     Fmul(nF), Flow(nF), Fupp(nF),
     &     xmul(n), xlow(n), xupp(n), rw(lenrw)
      character
     &     xnames(nxname)*8, Fnames(nFname)*8

*     ==================================================================
*     s3argA   checks the arguments for snOptA.
*
*     On exit, Errors says how many errors were encountered.
*     lvlSrt is an integer version of Start:
*        Start   lvlSrt
*        'Cold'       0
*        'Basis'      1
*        'Warm'       2
*        'Hot'        3
*     iw(gotFac,gotHes,gotScl) are set to 0 or 1.
*          33     Keep Factors of basis (LU)
*         303     Keep Hessian
*        3003     Keep Scales
*         333     Keep Factors and Hessian
*        3333     etc.
*           .     .
*           3     is treated as 3333
*
*     21 Dec 2002: First version of s3argA.
*     03 Aug 2003: snPRNT adopted.
*     13 Mar 2004: Start = 3111 options decoded.
*     03 Apr 2004: Prevent n*nF overflow.
*     02 Sep 2004: Length of str increased to 132.
*     04 Dec 2004: Current version of s3argA.
*     ==================================================================
      character
     &     str*132
      logical
     &     Named, ok
      integer
     &     argerr, Errs, i, idummy, is, j, js, k, n1, xMods, FMods
      double precision
     &     b1, b2, InfBnd
*     ------------------------------------------------------------------
      integer            COLD,       BASIS,      WARM,       HOT
      parameter         (COLD   = 0, BASIS  = 1, WARM   = 2, HOT    = 3)
      integer            gotFac,     gotHes,     gotScl
      parameter         (gotFac=230, gotHes=231, gotScl=232)
      double precision   zero
      parameter         (zero   = 0.0d+0)
      parameter         (idummy =  -11111)
*     ------------------------------------------------------------------
      InfBnd    = rw( 70) ! definition of an infinite bound
      argerr    = iw(106) ! maximum # errors in MPS data

      iExit     = 0
      Errors    = 0

*     The options haven't been checked yet.

      if (InfBnd .lt. zero) InfBnd = 1.0d+20
      if (argerr .lt. 0   ) argerr = 20 ! Print 20 lines max

      Named  = nxname .eq. n

*     ==================================================================
*     Decode 'Start'.
*     ==================================================================
      iw(gotFac) = 0
      iw(gotHes) = 0
      iw(gotScl) = 0

!     Determine the type of start.
!     Optional parameters take precedence.

      if (lvlSrt .ne. idummy) then ! lvlSrt was set as an option
         if (    lvlSrt .eq. COLD
     &     .or.  lvlSrt .eq. BASIS
     &     .or.  lvlSrt .eq. WARM ) then
*           Relax
         else if (lvlSrt .eq. HOT) then
            iw(gotFac) = 1
            iw(gotHes) = 1
            iw(gotScl) = 1
         else                      ! lvlSrt is an unrecognized option
            lvlSrt = idummy
         end if
      end if

      if (lvlSrt .eq. idummy) then ! lvlSrt is unset

         if (       Start .eq. COLD
     &        .or.  Start .eq. BASIS
     &        .or.  Start .eq. WARM ) then
            lvlSrt = Start
         else
            js     = Start
            j      = js - (js/10)*10
            js     = js - j

            if (j .eq. HOT) then
               lvlSrt     = j
               iw(gotFac) = 1
               iw(gotHes) = 1
               iw(gotScl) = 1

               j      = js - (js/100)*100
               js     = js - j

               if (j .eq. 0) then
                  iw(gotScl) = 0
               end if

               j      = js - (js/10)*10
               js     = js - j
               if (j .eq. 0) then
                  iw(gotHes) = 0
               end if

               if (js .eq. 0) then
                  iw(gotFac) = 0
               end if

            else
               lvlSrt = COLD
               write(str, 1000) Start
               call snPRNT( 13, str, iw, leniw )
            end if
         end if
      end if

!     Check the other arguments.

      if (nF .lt. 1) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'nF    ', nF
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (n .lt. 1) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'n     ', n
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

*     if (nS .lt. 0) then
*        Errors = Errors + 1
*        if (Errors .le. argerr) then
*           write(str, 1100) 'nS    ', nS
*           call snPRNT( 13, str, iw, leniw )
*        end if
*     end if

!     Test if neA or neG are bigger than n*nF. Avoid integer overflow.

      n1 = max(n,1)

      k  = mod(neA,n1)

      if (neA .lt. 0   .or.  (k.eq.0 .and. neA/n1 .gt. nF)
     &                 .or.  (k.gt.0 .and. neA/n1 .ge. nF)) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'neA   ', neA
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      k = mod(neG,n1)

      if (neG .lt. 0   .or.  (k.eq.0 .and. neG/n1 .gt. nF)
     &                 .or.  (k.gt.0 .and. neG/n1 .ge. nF)) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'neG   ', neG
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (nxname .ne. 1  .and.  nxname .ne. n ) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'nxname ', nxname
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (nFname .ne. 1  .and.  nFname .ne. nF) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'nFname ', nFname
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (      nxname .eq. 1  .and.  nFname .ne. 1
     &    .or.  nxname .ne. 1  .and.  nFname .eq. 1) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'nxname ', nxname
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (ObjRow .lt. 0  .or.  ObjRow  .gt. nF) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'ObjRow', ObjRow
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      xMods = 0
      FMods = 0

      if (lvlSrt .eq. COLD) then
         do j  = 1, n
            js = xstate(j)
            if (js .lt. -1  .or.  js .gt. 5) then
               xMods = xMods + 1
               xstate(j) = 0
            end if
         end do
      else if (lvlSrt .eq. WARM  .or.  lvlSrt .eq. HOT) then
         do j  = 1, n
            js = xstate(j)
            if (js .lt. 0  .or.  js .gt. 3) then
               xMods = xMods + 1
               xstate(j) = 0
            end if
         end do
         do i  = 1, nF
            is = Fstate(i)
            if (is .lt. 0  .or.  is .gt. 3) then
               FMods = FMods + 1
               Fstate(i) = 0
            end if
         end do
      end if

      if (xMods .gt. 0) then
         write(str, 1200) xMods
         call snPRNT( 13, str, iw, leniw )
      end if

      if (FMods .gt. 0) then
         write(str, 1210) FMods
         call snPRNT( 13, str, iw, leniw )
      end if

*     ------------------------------------------------------------------
*     Check the bounds on all variables and constraints.
*     ------------------------------------------------------------------
      Errs = 0

      do j  = 1, n
         b1 = xlow(j)
         b2 = xupp(j)
         ok = b1 .lt. b2  .or.
     &        b1 .eq. b2  .and.  abs(b1) .lt. InfBnd

         if (.not. ok)  then
            Errs   = Errs   + 1
            Errors = Errors + 1

            if (Errors .le. argerr) then
               if ( Named ) then
                  if (b1 .eq. b2) then
                     write(str, 1310) xNames(j), b1, InfBnd
                     call snPRNT( 13, str, iw, leniw )
                  else
                     write(str, 1315) xNames(j), b1, b2
                     call snPRNT( 13, str, iw, leniw )
                  end if
               else
                  if (b1 .eq. b2) then
                     write(str, 1300) j, b1, InfBnd
                     call snPRNT( 13, str, iw, leniw )
                  else
                     write(str, 1305) j, b1, b2
                     call snPRNT( 13, str, iw, leniw )
                  end if
               end if
            end if
         end if
      end do

      if (Errs .gt. 0) then
         write(str, 1600) Errs
         call snPRNT( 13, str, iw, leniw )
      end if

      Errs = 0

      do i  = 1, nF
         b1 = Flow(i)
         b2 = Fupp(i)

         ok = b1 .lt. b2  .or.
     &        b1 .eq. b2  .and.  abs(b1) .lt. InfBnd

         if (.not. ok)  then
            Errs   = Errs   + 1
            Errors = Errors + 1

            if (Errors .le. argerr) then
               if (Named) then
                  if (b1 .eq. b2) then
                     write(str, 1510) Fnames(i), b1, InfBnd
                     call snPRNT( 13, str, iw, leniw )
                  else
                     write(str, 1515) Fnames(i), b1, b2
                     call snPRNT( 13, str, iw, leniw )
                  end if
               else
                  if (b1 .eq. b2) then
                     write(str, 1510) 'F(i) ', b1, InfBnd
                     call snPRNT( 13, str, iw, leniw )
                  else
                     write(str, 1515) 'F(i) ', b1, b2
                     call snPRNT( 13, str, iw, leniw )
                  end if
               end if
            end if
         end if
      end do

      if (Errs .gt. 0) then
         write(str, 1610) Errs
         call snPRNT( 13, str, iw, leniw )
      end if

      if (Errors .ge. argerr) then
         call snPRNT( 3, ' and so on ...', iw, leniw )
      end if

      if (Errors .gt. 0) then
*        ----------------------------
*        Invalid arguments.
*        ----------------------------
         iExit = 91             ! invalid input argument
      end if

      return

 1000 format(' XXX Start parameter not recognized:  ', i5)
 1100 format(' XXX  Argument out of range:  ', a6, ' = ', i6)
 1200 format(' XXX  Invalid argument xstate: ', i6,
     &       ' elements modified to be in range.')
 1210 format(' XXX  Invalid argument Fstate: ', i6,
     &       ' elements modified to be in range.')

 1300 format(' XXX  The equal bounds on variable ', i6,
     &       '  are infinite.   Bounds =', g16.7,
     &       '  InfBnd =', g16.7)
 1305 format(' XXX  The bounds on variable ', i6,
     &       '  are inconsistent.  xlow =', g16.7, '  xupp =', g16.7)

 1310 format(' XXX  The equal bounds  ', a8,
     &       '  are infinite.   Bounds =', g16.7,
     &       '  InfBnd =', g16.7)
 1315 format(' XXX  The bounds on ', a8,
     &       '  are inconsistent.  xlow =', g16.7, '  xupp =', g16.7)
 1500 format(' XXX  The equal bounds on function ', i6,
     &       '  are infinite.   Bounds =', g16.7,
     &       '  InfBnd =', g16.7)
 1505 format(' XXX  The bounds on function ', i6,
     &       '  are inconsistent.  Flow =', g16.7, '  Fupp =', g16.7)
 1510 format(' XXX  The equal bounds on  ', a8,
     &       '  are infinite.   Bounds =', g16.7,
     &       '  InfBnd =', g16.7)
 1515 format(' XXX  The bounds on  ', a8,
     &       '  are inconsistent.  Flow =', g16.7, '  Fupp =', g16.7)
 1600 format(' XXX  Invalid arguments xlow, xupp: ', i6,
     &       ' inconsistent bounds or infinite equal bounds.')
 1610 format(' XXX  Invalid arguments Flow, Fupp: ', i6,
     &       ' inconsistent bounds or infinite equal bounds.')

      end ! subroutine s3argA

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3prtA
     &   ( m, n, nnCon, nnJac, nnObj, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     m, n, nnCon, nnJac, nnObj, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)

*     ==================================================================
*     s3prtA prints the settings of the optional parameters for the
*     snoptA wrapper.
*
*     See  snworkspace.doc  for full documentation of cw, iw and rw.
*
*     15 Nov 1991: First version (s8prnt).
*     10 Dec 2002: More LU pivoting options.
*     03 Aug 2003: snPRNT adopted.
*     21 Dec 2003: Current version of s3prtA.
*     ==================================================================
      integer
     &     cgItmx, DerOpt, iCrash, iBack, iDump, iLoadB, iInsrt, iNewB,
     &     iOldB, iPnch, iPrint, iStdi, iStdo, iSoln, iSpecs, itnlim,
     &     kchk, kdegen, kFac, klog, kreset, ksav, kSumm,
     &     lprDbg, lprPrm, lvlHes, lvlPiv, lvlPPm, lvlPre,
     &     lvlSch, lvlScl, lvlSrt, lvlSys, lvlTim, lvlVer,
     &     maxmn, maxR, maxS, mflush, minmax, minPrc, MjrPrt, mMajor,
     &     mMinor, MnrPrt, mQNmod, mNewSB, nnL,
     &     nParPr, nPr1, nPr2, ObjRow, QPslvr
      double precision
     &     bigdx, bigFx,
     &     eps, epsrf, etarg,
     &     fdint1, fdint2, scltol,
     &     tolCG, tCrash, tolCon, tolFac, tolNLP, tolpiv, tolQP,
     &     tolSwp, tolUpd, tolx, Utol1, viLim, wolfeG,
     &     wtInf0, xdlim, xPen0
      character
     &     str1*132, str2*132, str3*132, str4*132, str5*132, str6*132
*     ------------------------------------------------------------------
      character          Hestyp(1:3)*24,  lsrch(0:1)*24, pivtyp(0:3)*24,
     &                   prbtyp(1:3)*24, QPtype(0:2)*24, SrtTyp(0:3)*24,
     &                   NoYes (0:1)*3
      data               Hestyp /' Limited-Memory Hessian.',
     &                           ' Full-Memory Hessian....',
     &                           ' Exact Hessian..........'/
      data               lsrch  /' Nonderiv.  linesearch..',
     &                           ' Derivative linesearch..'/
      data               pivtyp /' LU partial  pivoting...',
     &                           ' LU rook     pivoting...',
     &                           ' LU complete pivoting...',
     &                           ' LU diagonal pivoting...'/
      data               prbtyp /' Maximize...............',
     &                           ' Feasible point only....',
     &                           ' Minimize...............'/
      data               QPtype /' QPsolver Cholesky......',
     &                           ' QPsolver CG............',
     &                           ' QPsolver QN............'/
      data               SrtTyp /' Cold start.............',
     &                           ' Basis file.............',
     &                           ' Warm start.............',
     &                           ' Hot start..............'/
      data               NoYes  /' No', 'Yes'/
*     ------------------------------------------------------------------
*     Set some local machine-dependent constants.

      eps       = rw(  1) ! machine precision.  IEEE DP  2.22e-16
      tolQP     = rw( 52) ! Minor Phase 2 Opt tol
      tolNLP    = rw( 53) ! Major Optimality tolerance
      tolCG     = rw( 54) ! cg tolerance
      tolx      = rw( 56) ! Minor feasibility tolerance.
      tolCon    = rw( 57) ! Major feasibility tolerance.
      tolpiv    = rw( 60) ! excludes small elements of y.
      tCrash    = rw( 62) ! crash tolerance.
      tolswp    = rw( 65) ! LU swap tolerance.
      tolFac    = rw( 66) ! LU factor tolerance.
      tolUpd    = rw( 67) ! LU update tolerance.
      bigFx     = rw( 71) ! unbounded objective.
      bigdx     = rw( 72) ! unbounded step.
      epsrf     = rw( 73) ! relative function precision.
      fdint1    = rw( 76) ! (1) forwrd diff. interval
      fdint2    = rw( 77) ! (2) cntrl  diff. interval
      xdlim     = rw( 80) ! Step limit
      vilim     = rw( 81) ! violation limit
      etarg     = rw( 83) ! Quasi-Newton QP rg tolerance
      wolfeG    = rw( 84) ! line search tolerance.
      wtInf0    = rw( 88) ! infeasibility weight
      xPen0     = rw( 89) ! initial penalty parameter.
      scltol    = rw( 92) ! scale tolerance.
      Utol1     = rw(154) ! abs tol for small diag of U.

      iStdi     = iw(  9) ! Standard Input
      iStdo     = iw( 10) ! Standard Output
      iSpecs    = iw( 11) ! Specs (options) file
      iPrint    = iw( 12) ! Print file
      maxR      = iw( 52) ! max columns of R.
      maxS      = iw( 53) ! max # of superbasics
      mQNmod    = iw( 54) ! (ge 0) max # of BFGS updates
      QPslvr    = iw( 55) ! 0(1/2) => QPsolver
      kchk      = iw( 58) ! check (row) frequency
      kFac      = iw( 59) ! factorization frequency
      ksav      = iw( 60) ! save basis map
      klog      = iw( 61) ! log/print frequency
      kSumm     = iw( 62) ! Summary print frequency
      kDegen    = iw( 63) ! max. expansions of featol
      kReset    = iw( 64) ! Hessian frequency
      mFlush    = iw( 66) ! Hessian flush
      lvlSrt    = iw( 69) ! = 0:1:2:3 => cold:warm:basis:hot start
      lvlSys    = iw( 71) ! > 0   => print system info
      lvlHes    = iw( 72) ! 0,1,2  => LM, FM, Exact Hessian
      lvlScl    = iw( 75) ! scale option
      lvlSch    = iw( 76) ! >0     => use derivatives in the line search
      lvlPre    = iw( 77) ! >0    => QN preconditioned CG
      lvlVer    = iw( 78) ! Verify level
      lvlPPm    = iw( 79) ! Proximal Point method for x0
      lvlPiv    = iw( 80) ! 0(1) LU threshold partial(complete) pivoting
      lprPrm    = iw( 81) ! > 0    => parms are printed
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
      iSoln     = iw(131) ! solution file
      lvlTim    = iw(182) ! Timing level
*     ------------------------------------------------------------------

      if (iPrint .le. 0 .or. MjrPrt .eq. 0 .or. lprPrm .eq. 0) return

      nnL = max( nnJac, nnObj )

      minPrc = 10
      nPr1   = n / nParPr
      nPr2   = m / nParPr
      if (max( nPr1, nPr2 ) .lt. minPrc) then
         maxmn  = max( m, n )
         nParPr = maxmn / min( maxmn, minPrc )
         nPr1   = n / nParPr
         nPr2   = m / nParPr
      end if

*     ==================================================================
*     Print parameters except if PRINT LEVEL = 0
*     or SUPPRESS PARAMETERS was specified.
*     ==================================================================
      call s1page( 1, iw, leniw )
      call snPRNT( 1, ' Parameters', iw, leniw )
      call snPRNT( 1, ' ==========', iw, leniw )

*     --------------------
*     Files.
*     --------------------
      call snPRNT(11, ' Files', iw, leniw )
      call snPRNT( 1, ' -----', iw, leniw )
      write(str1, 2110) iSoln , iOldB , iStdi
      write(str2, 2120) iInsrt, iNewB , iPrint
      write(str3, 2130) iPnch , iBack , iSpecs
      write(str4, 2140) iLoadB, iDump , iStdo
      call snPRNT( 1, str1, iw, leniw )
      call snPRNT( 1, str2, iw, leniw )
      call snPRNT( 1, str3, iw, leniw )
      call snPRNT( 1, str4, iw, leniw )

*     --------------------
*     Frequencies.
*     --------------------
      call snPRNT(11, ' Frequencies', iw, leniw )
      call snPRNT( 1, ' -----------', iw, leniw )
      write(str1, 2210) klog  , kchk  , ksav
      write(str2, 2220) kSumm , kfac  , kDegen
      call snPRNT( 1, str1, iw, leniw )
      call snPRNT( 1, str2, iw, leniw )

*     --------------------
*     QP subproblems.
*     --------------------
      call snPRNT(11, ' QP subproblems', iw, leniw )
      call snPRNT( 1, ' --------------', iw, leniw )
      write(str1, 2310) QPtype(QPslvr)
      write(str2, 2320) scltol, tolx  , itnlim
      write(str3, 2330) lvlScl, tolQP , MnrPrt
      write(str4, 2340) tCrash, tolpiv, nParPr
      write(str5, 2350) iCrash, wtInf0, nPr1
      write(str6, 2360) mNewSB, nPr2
      call snPRNT( 1, str1, iw, leniw )
      call snPRNT( 1, str2, iw, leniw )
      call snPRNT( 1, str3, iw, leniw )
      call snPRNT( 1, str4, iw, leniw )
      call snPRNT( 1, str5, iw, leniw )
      call snPRNT( 1, str6, iw, leniw )

*     -----------------------
*     CG QP solver.
*     -----------------------
      if (QPslvr .eq. 1  .or.  maxR .lt. maxS) then
         call snPRNT(11, ' Conjugate-gradient QP solver', iw, leniw )
         call snPRNT( 1, ' ----------------------------', iw, leniw )
         write(str1, 2380) etarg,   tolCG, cgItmx
         write(str2, 2390)                 lvlPre
         call snPRNT( 1, str1, iw, leniw )
         call snPRNT( 1, str2, iw, leniw )
      end if

*     --------------------
*     SQP method.
*     --------------------
      call snPRNT(11, ' The SQP Method', iw, leniw )
      call snPRNT( 1, ' --------------', iw, leniw )
      write(str1, 2410) prbtyp(2+minmax),SrtTyp(lvlSrt),lvlPPm
      write(str2, 2420) nnObj , ObjRow, epsrf
      write(str3, 2430) bigdx , maxS  , fdint1
      write(str4, 2440) bigFx , maxR  , fdint2
      call snPRNT( 1, str1, iw, leniw )
      call snPRNT( 1, str2, iw, leniw )
      call snPRNT( 1, str3, iw, leniw )
      call snPRNT( 1, str4, iw, leniw )
      write(str1, 2450) xdlim , lsrch(lvlSch) , DerOpt
      write(str2, 2460) mMajor, wolfeG, lvlVer
      write(str3, 2470) mMinor, xPen0 , MjrPrt
      write(str4, 2480) tolNLP
      call snPRNT( 1, str1, iw, leniw )
      call snPRNT( 1, str2, iw, leniw )
      call snPRNT( 1, str3, iw, leniw )
      call snPRNT( 1, str4, iw, leniw )

*     --------------------
*     Hessian approximation.
*     --------------------
      if (nnL .gt. 0) then
         call snPRNT(11, ' Hessian Approximation', iw, leniw )
         call snPRNT( 1, ' ---------------------', iw, leniw )
         write(str1, 2510) Hestyp(lvlHes+1), mQNmod, kReset
         write(str2, 2520) mFlush
         call snPRNT( 1, str1, iw, leniw )
         call snPRNT( 1, str2, iw, leniw )
      end if

*     --------------------
*     Nonlinear constraints.
*     --------------------
      if (nnCon .gt. 0) then
         call snPRNT(11, ' Nonlinear constraints', iw, leniw )
         call snPRNT( 1, ' ---------------------', iw, leniw )
         write(str1, 2610) nnCon , tolCon, vilim
         write(str2, 2620) nnJac
         call snPRNT( 1, str1, iw, leniw )
         call snPRNT( 1, str2, iw, leniw )
      end if

*     --------------------
*     Miscellaneous
*     --------------------
      call snPRNT(11, ' Miscellaneous', iw, leniw )
      call snPRNT( 1, ' -------------', iw, leniw )
      write(str1, 2710) tolFac, Utol1 , lvlTim
      write(str2, 2720) tolUpd, tolswp, lprDbg
      write(str3, 2730) pivtyp(lvlPiv), eps, NoYes(lvlSys)
      call snPRNT( 1, str1, iw, leniw )
      call snPRNT( 1, str2, iw, leniw )
      call snPRNT( 1, str3, iw, leniw )

      write(str1, 3000) lvlScl, nParPr
      call snPRNT( 22, str1, iw, leniw )

      return

 2110 format(
     &     ' Solution file..........', i10, 6x,
     &     ' Old basis file ........', i10, 6x,
     &     ' Standard input.........', i10)
 2120 format(
     &     ' Insert file............', i10, 6x,
     &     ' New basis file ........', i10, 6x,
     &     ' (Printer)..............', i10)
 2130 format(
     &     ' Punch file.............', i10, 6x,
     &     ' Backup basis file......', i10, 6x,
     &     ' (Specs file)...........', i10)
 2140 format(
     &     ' Load file..............', i10, 6x,
     &     ' Dump file..............', i10, 6x,
     &     ' Standard output........', i10)

 2210 format(
     &     ' Print frequency........', i10, 6x,
     &     ' Check frequency........', i10, 6x,
     &     ' Save new basis map.....', i10)
 2220 format(
     &     ' Summary frequency......', i10, 6x,
     &     ' Factorization frequency', i10, 6x,
     &     ' Expand frequency.......', i10)

 2310 format(a24)
 2320 format(
     &     ' Scale tolerance........', 0p, f10.3, 6x,
     &     ' Minor feasibility tol..', 1p, e10.2, 6x,
     &     ' Iteration limit........', i10)
 2330 format(
     &     ' Scale option...........', i10,       6x,
     &     ' Minor optimality  tol..', 1p, e10.2, 6x,
     &     ' Minor print level......', i10)
 2340 format(
     &     ' Crash tolerance........', 0p, f10.3, 6x,
     &     ' Pivot tolerance........', 1p, e10.2, 6x,
     &     ' Partial price..........', i10)
 2350 format(
     &     ' Crash option...........', i10,       6x,
     &     ' Elastic weight.........', 1p, e10.2, 6x,
     &     ' Prtl price section ( A)', i10)
 2360 format(
     &     40x,
     &     ' New superbasics........', i10,       6x,
     &     ' Prtl price section (-I)', i10)

 2380 format(
     &     ' Subspace tolerance.....', 0p, f10.5, 6x,
     &     ' CG tolerance...........', 1p, e10.2, 6x,
     &     ' CG Iterations..........', i10)

 2390 format(
     &     80x,
     &     ' CG preconditioning.....', i10)

 2410 format(a24, 16x, a24, 16x,
     &     ' Proximal Point method..', i10)
 2420 format(
     &     ' Nonlinear objectiv vars', i10,       6x,
     &     ' Objective Row..........', i10,       6x,
     &     ' Function precision.....', 1p, e10.2)
 2430 format(
     &     ' Unbounded step size....', 1p, e10.2, 6x,
     &     ' Superbasics limit......', i10,       6x,
     &     ' Difference interval....', 1p, e10.2)
 2440 format(
     &     ' Unbounded objective....', 1p, e10.2, 6x,
     &     ' Reduced Hessian dim....', i10,       6x,
     &     ' Central difference int.', 1p, e10.2)
 2450 format(
     &     ' Major step limit.......', 1p, e10.2, 6x,
     &            a11,' linesearch..',           16x,
     &     ' Derivative option......', i10)
 2460 format(
     &     ' Major iterations limit.', i10,       6x,
     &     ' Linesearch tolerance...', 0p, f10.5, 6x,
     &     ' Verify level...........', i10)
 2470 format(
     &     ' Minor iterations limit.', i10,       6x,
     &     ' Penalty parameter......', 1p, e10.2, 6x,
     &     ' Major Print Level......', i10)
 2480 format(
     &     40x,
     &     ' Major optimality tol...', 1p, e10.2)

 2510 format(
     &     a24,                                  16x,
     &     ' Hessian updates........', i10,       6x,
     &     ' Hessian frequency......', i10)
 2520 format(
     &     80x,
     &     ' Hessian flush..........', i10)

 2610 format(
     &     ' Nonlinear constraints..', i10,       6x,
     &     ' Major feasibility tol..', 1p, e10.2, 6x,
     &     ' Violation limit........',     e10.2)
 2620 format(
     &     ' Nonlinear Jacobian vars', i10)

 2710 format(
     &     ' LU factor tolerance....', 0p, f10.2, 6x,
     &     ' LU singularity tol.....', 1p, e10.2, 6x,
     &     ' Timing level...........', i10)
 2720 format(
     &     ' LU update tolerance....', 0p, f10.2, 6x,
     &     ' LU swap tolerance......', 1p, e10.2, 6x,
     &     ' Debug level............', i10)
 2730 format(
     &     a24,                                  16x,
     &     ' eps (machine precision)', 1p, e10.2, 6x,
     &     ' System information.....', 7x, a3 )
 3000 format(' Scale option', i3, ',    Partial price', i4)

      end ! subroutine s3prtA

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3sizA
     &   ( iExit, n, nF, nkx, ObjRow,
     &     iAfun, jAvar, lenA, neA, iGfun, jGvar, lenG, neG,
     &     m, negCon, ne, nnCon, nnJac, nnObj, iObj,
     &     kx, leniw, iw )

      implicit
     &     none
      integer
     &     iExit, iObj, lenA, neA, lenG, leniw, m, n, nF, nkx,
     &     neG, negCon, ne, nnCon, nnJac, nnObj, ObjRow,
     &     iAfun(lenA), jAvar(lenA), iGfun(lenG), jGvar(lenG),
     &     kx(nkx), iw(leniw)

*     ==================================================================
*     s3sizA runs through the data once in order to compute the
*     dimensions of the quantities:
*     Jcol    is the      m x n      matrix of QP constraints.
*     gCon    is the  nnCon x nnJac  nonlinear part of the Jacobian.
*     gObj    is the  nnObj x 1      objective gradient.
*
*        negCon,  ne,  nnCon, nnJac,  nnObj  and  iObj.
*
*     On entry:
*     ---------
*     ObjRow  is the (natural) row number of the objective.
*
*     nkx     (= n + nF) is the length of kx.
*
*     On exit:
*     --------
*     m       is the row dimension of Jcol.
*                  m = nF - 1    if ObjRow > 0 and iObj =0,
*                  m = nF,       otherwise.
*
*     iObj    is the position of the linear objective row in J.
*             iObj = 0 means that there is no linear objective row.
*
*     kx(nkx) is a list of the indices that defines the position of
*             each user-defined F(i) and x(j) in  Jcol, i.e.,
*             k = kx(j)   => variable   j is the k-th column of Jcol.
*             l = kx(n+i) => constraint i is the l-th row    of Jcol.
*
*     08 Nov 1998: First version of s3sizA.
*     18 Jul 2000: First version for SNOPT 6.1.
*     03 Aug 2003: snPRNT adopted.
*     03 Aug 2003: Current version of s3sizA
*     ==================================================================
      character
     &     str*132
      integer
     &     i, iN, j, jS, jSec1, jSec2, jSec3, jSec4, jN,
     &     k, nBoth, nCol, Col, nJac, nnL, nObj, nRow, Row
**     ------------------------------------------------------------------
      integer            NotSet
      parameter         (NotSet =-1)
      integer            Linear,     Jac,     Obj,     Both
      parameter         (Linear = 0, Jac = 1, Obj = 2, Both = 3 )
*     ------------------------------------------------------------------
      iExit    = 0

*     ==================================================================
*     Count the number of nonlinear objective and Jacobian variables.
*     put the nonlinear rows first, recording their positions in kx.

*     Check if variable is nonlinear in the objective, the Jacobian,
*     or both.
*
*        kx             type of variable
*        --             ----------------
*        Linear         linear    Jacobian,   linear    Objective
*        Jac            nonlinear Jacobian,   linear    Objective
*        Obj            linear    Jacobian,   nonlinear Objective
*        Both           nonlinear Jacobian,   nonlinear Objective
*     ==================================================================
      nObj   = 0
      nJac   = 0
      nBoth  = 0

      call iload ( nkx, NotSet, kx, 1 )

      iObj   = 0
      nnCon  = 0
      negCon = 0
      nRow   = 0
      nCol   = 0

*     ------------------------------------------------------------------
*     Process the nonlinear Jacobian elements.  The array  G  holds the
*     nonzero nonlinear elements of the derivative matrix G(x) that will
*     include the objective gradient (if there is one).
*     kx(1:n) and kx(n+1:nkx) give the coordinates of J(i,j).
*     ------------------------------------------------------------------
      do k  = 1, neG
         iN = iGfun(k)
         jN = jGvar(k)

         if (jN .lt. 0  .or.  jN .gt. n    .or.
     &       iN .lt. 0  .or.  iN .gt. nF       ) then
            write(str, 2046) k, iN, jN
            call snPRNT( 3, str, iw, leniw )
            iExit = 91          ! Invalid input argument
            go to 900
         end if

         Col = jN

         if (iN .eq. ObjRow) then
*           --------------------------
*           Element of gObj.
*           --------------------------
            if (kx(Col) .eq. NotSet) then
               kx(Col) = Obj
               nObj    = nObj + 1
               nCol    = nCol + 1
            else if (kx(Col) .eq. Jac) then
               kx(Col) = Both
               nObj    = nObj  + 1
               nBoth   = nBoth + 1
            end if
         else
*           --------------------------
*           Element of Jcol and gCon.
*           --------------------------
            negCon = negCon + 1

*           See if this element starts a new nonlinear column.

            if (kx(Col) .eq. NotSet) then
               kx(Col) = Jac
               nJac    = nJac + 1
               nCol    = nCol + 1
            else if (kx(Col) .eq. Obj) then
               kx(Col) = Both
               nJac    = nJac  + 1
               nBoth   = nBoth + 1
            end if

*           See if this is the first element in a new nonlinear row.

            Row = n + iN

            if (kx(Row) .eq. NotSet) then
               nnCon   = nnCon + 1
               nRow    = nRow  + 1
               kx(Row) = nRow
            end if
         end if
      end do

*     ------------------------------------------------------------------
*     Now we know all the dimensions of the nonlinear part.
*     Figure out the variable order in the first nnL columns.
*     Set  nnObj,  nnCon,  nnJac  and  nnL.
*     ------------------------------------------------------------------
      nnL   = nCol
      nnJac = nJac
      if (nObj .gt. nBoth) then
         nnObj = nJac + (nObj - nBoth)
      else
         nnObj = nBoth
      end if

*     ------------------------------------------------------------------
*     Process the constant part of the Jacobian.
*     The linear rows follow the nonlinear rows.
*     ------------------------------------------------------------------
      ne = negCon

      do k  = 1, neA
         iN = iAfun(k)
         jN = jAvar(k)

*        See if this element starts a new linear column.

         Col = jN

         if (kx(Col) .eq. NotSet) then
            nCol    = nCol + 1
            kx(Col) = Linear
         end if

         ne = ne + 1    ! All linear terms are held in J.

*        See if this element starts a new linear row.

         Row = n + iN

         if (kx(Row) .eq. NotSet) then
            nRow    = nRow + 1
            kx(Row) = nRow
            if (iN .eq. ObjRow) iObj = nRow
         end if
      end do

*     ------------------------------------------------------------------
*     If necessary, add a dummy free linear row.
*     As nF > 0,  F must consist of just the nonlinear objective and
*     ObjRow = 1.  May as well make the free row the objective row.
*     ------------------------------------------------------------------
      if (nRow .eq. 0) then
         nRow = nRow + 1
         iObj = 1
         ne   = 1
      end if

*     ------------------------------------------------------------------
*     Set the row dimension of Jcol.
*     ------------------------------------------------------------------
      if (ObjRow .gt. 0  .and.  iObj .eq. 0) then
         m  = nF - 1
      else
         m  = nF
      end if

*     ------------------------------------------------------------------
*     Process the objective.
*     The (n+ObjRow)-th element of kx is 0 if there is no linear obj.
*     ------------------------------------------------------------------
      if (ObjRow .gt. 0) kx(n+ObjRow) = iObj

*     ------------------------------------------------------------------
*     Check for unusual cases.
*     ------------------------------------------------------------------
      if (nCol .lt. n) then
*        --------------------
*        Empty columns.
*        --------------------
         call snPRNT( 13, ' ', iw, leniw )

         do jN = 1, n
            if (kx(jN) .eq. NotSet) then
               write(str, 1200) jN
               call snPRNT( 3, str, iw, leniw )
               nCol   = nCol + 1
               kx(jN) = Linear
            end if
         end do
      end if

      if (nRow .lt. m) then
*        --------------------
*        Empty rows.
*        --------------------
         call snPRNT( 13, ' ', iw, leniw )

         do iN  = 1, m
            Row = n + iN
            if (kx(Row) .eq. NotSet) then
               write(str, 1300) iN
               call snPRNT( 3, str, iw, leniw )
               nRow    = nRow + 1
               kx(Row) = nRow
            end if
         end do
      end if

*     ==================================================================
*     Re-order the variables into four sections:
*          Section 1 - type  Both
*          Section 2 - type  Jac
*          Section 3 - type  Obj
*          Section 4 - type  Linear
*
*     jSecn+1  points to the next free slot in section n.
*     ==================================================================
      jSec1   = 0
      jSec2   = nBoth
      jSec3   = nnJac
      jSec4   = nnL

*     Loop through the columns and initialize kx, which holds
*     the final column order.

      do jN = 1, n
         if (     kx(jN) .eq. Both) then
            jSec1 = jSec1 + 1
            jS    = jSec1
         else if (kx(jN) .eq. Jac ) then
            jSec2 = jSec2 + 1
            jS    = jSec2
         else if (kx(jN) .eq. Obj ) then
            jSec3 = jSec3 + 1
            jS    = jSec3
         else
            jSec4 = jSec4 + 1
            jS    = jSec4
         end if
         kx(jN) = jS
      end do

*     ------------------------------------------------------------------
*     Check for constant elements in the nonlinear part of Jcol.
*     ------------------------------------------------------------------
      do k   = 1, neA
         iN  = iAfun(k)
         jN  = jAvar(k)

         Row = n + iN
         Col = jN

         i   = kx(Row)
         j   = kx(Col)

         if (i .le. nnCon  .and.  j .le. nnJac) then
            negCon = negCon + 1
         end if
      end do

      if (ObjRow .eq. 0) then
*        -------------------------------------------
*        No objective row.  SNOPT  will provide one.
*        -------------------------------------------
         write(str, 1500)
         call snPRNT( 3, str, iw, leniw )

      else if (kx(n+ObjRow) .eq. NotSet) then
*        --------------------
*        Empty objective row.
*        --------------------
         kx(n+ObjRow) = 0
         write(str, 1400) ObjRow
         call snPRNT( 3, str, iw, leniw )
      end if

  900 return

 1200 format(' ===>  WARNING - variable   ', i6,
     &             ' has no Jacobian entries.')
 1300 format(' ===>  WARNING - constraint ', i6,
     &             ' has no Jacobian entries.')
 1400 format(' ===>  WARNING - Objective Row ', i6,
     &             ' has no  entries.')
 1500 format(' ===>  No objective row specified ---',
     &             ' finding a feasible point.')
 2046 format(' XXX  First derivative element    k = ', i6,
     &             ',  row ',   i6, ' column ', i6,
     &             ' is out of range.' )

      end ! subroutine s3sizA

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3bldA
     &   ( ObjRow, n, nkx, nnCon, nnJac, kx, nGlin,
     &     iAfun, jAvar, lenA, neA, A, iGfun, jGvar, lenG, neG,
     &     ne, nlocJ, locJ, indJ, Jcol,
     &     negCon, nlocG, locG, indG, iy )

      implicit
     &     none
      integer
     &     lenA, lenG, n, nkx, neA, ne, neG, negCon, nnCon,
     &     nnJac, nlocJ, nlocG, ObjRow, kx(nkx), nGlin(nlocG),
     &     iAfun(lenA), jAvar(lenA), iGfun(lenG), jGvar(lenG),
     &     locG(nlocG), indG(negCon), locJ(nlocJ), indJ(ne), iy(n)
      double precision
     &     Jcol(ne), A(lenA)

*     ==================================================================
*     s3bldA builds the data structure (locJ, indJ, Jcol) for the full
*     sparse Jacobian.
*
*     On entry,
*      neG             holds the number of nozero nonlinear derivative
*                      elements.
*      G, iGfun, jGvar hold the neG user-defined coordinates (Gij,i,j).
*                      neG  must be nonnegative.
*      A, iAfun, jAvar hold the neA user-defined coordinates (Aij,i,j).
*                      neA  must be nonnegative.
*
*     25 Sep 1998: First version of s3bldA.
*     18 Jul 2000: First version for SNOPT 6.1.
*     15 Jun 2001: Current version.
*     ==================================================================
      integer
     &     i, iN, j, jN, k, nColj, neJ, nextG, nextJ, Row, Col
*     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero   =  0.0d+0)
*     ------------------------------------------------------------------
*     Count the number of entries in each column (again).
*     Objective elements are excluded from  gCon and Jcol.
*     ------------------------------------------------------------------
      call iload ( n, (0), locJ, 1 )

      if (nnJac .gt. 0) then
         do k   = 1, neG
            iN  = iGfun(k)
            jN  = jGvar(k)

            Col = jN
            j   = kx(Col)
            if (iN .ne. ObjRow) locJ(j) = locJ(j) + 1
         end do

         call icopy ( nnJac, locJ, 1, locG, 1 )
      end if

*     ------------------------------------------------------------------
*     Now count the constant derivatives.
*     nGlin(j) = the number of constant elements in column j of gCon.
*     ------------------------------------------------------------------
      if (nnJac .gt. 0) then
         call iload ( nnJac, (0), nGlin, 1 )
         call iload ( nnJac, (0), iy   , 1 )
      end if

      do k   = 1, neA
         iN  = iAfun(k)
         jN  = jAvar(k)

         Row = n + iN
         Col = jN

         i   = kx(Row)
         j   = kx(Col)

         locJ(j) = locJ(j) + 1
         if (i .le. nnCon  .and.  j .le. nnJac) then
            locG(j)  = locG(j)  + 1
            nGlin(j) = nGlin(j) + 1
         end if
      end do

*     ------------------------------------------------------------------
*     Initialize the column pointers.
*     ------------------------------------------------------------------
      locJ(nlocJ)  = ne + 1
      do j = n, 1, -1
         locJ(j) = locJ(j+1) - locJ(j)
      end do

      neJ = 0

      if (negCon .gt. 0) then
         locG(nlocG) = negCon + 1
         do i = nnJac, 1, -1
            locG(i)  = locG(i+1) - locG(i)
         end do

*        ---------------------------------------------------------------
*        Assign the nonlinear constraint derivatives to  indJ and Jcol.
*        locJ(j) points to the next nonzero in column i.
*        ---------------------------------------------------------------
         do k   = 1, neG
            iN  = iGfun(k)
            jN  = jGvar(k)

            Col = jN
            j   = kx(Col)

            if (iN .ne. ObjRow) then
               Row         = n + iN
               i           = kx(Row)

               nextJ       = locJ(j)
               indJ(nextJ) = i
               Jcol(nextJ) = zero
               locJ(j)     = nextJ + 1
               neJ         = neJ   + 1

               nextG       = locG(j)
               indG(nextG) = i
               locG(j)     = nextG + 1
            end if
         end do
      end if

*     ------------------------------------------------------------------
*     Place any linear elements after the nonlinear elements.
*     ------------------------------------------------------------------
*     locJ (j) points to the next nonzero in the    linear part.
*     iy(j)    points to the next nonzero in the nonlinear part.

      if (nnJac .gt. 0) then
         do j = 1, nnJac
            nColj = nGlin(j)
            if (nColj .gt. 0) then
               nextJ   = locJ(j)
               locJ(j) = locJ(j) + nColj
               iy(j)   = nextJ
            end if
         end do
      end if

      do  k = 1, neA
         iN = iAfun(k)
         jN = jAvar(k)

         Row = n + iN
         Col = jN

         i   = kx(Row)
         j   = kx(Col)

*        Assign the elements of  indJ  and  Jcol.

         if (i .le. nnCon  .and.  j .le. nnJac) then
            nextJ       = iy(j)
            indJ(nextJ) = i
            Jcol(nextJ) = A(k)
            iy(j)       = nextJ + 1
            neJ         = neJ   + 1

            nextG       = locG(j)
            indG(nextG) = i
            locG(j)     = nextG + 1
         else
            nextJ       = locJ(j)
            indJ(nextJ) = i
            Jcol(nextJ) = A(k)
            locJ(j)     = nextJ + 1
            neJ         = neJ   + 1
         end if
      end do

*     Check for a dummy row.  In this case, ne = 1 is set in s3SizA.

      if (neJ .eq. 0) then
         Jcol(ne) = zero
         indJ(ne) = 1
      end if

*     Reset the column pointers.

      do j = n, 2, -1
         locJ(j) = locJ(j-1)
      end do
      locJ(1) = 1

      do j = nnJac, 2, -1
         locG(j) = locG(j-1)
      end do
      locG(1) = 1

      end ! subroutine s3bldA

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3inA
     &   ( Start, iObj,
     &     m, n, nb, nnCon, nF, nkx, kxN, infBnd,
     &     xnames, Fnames, Names, nxname, nFname, nName,
     &     xlow, xupp, Flow, Fupp, bl, bu, xstate, Fstate,
     &     hs, xN, F, x, Fx, Fmul, pi )

      implicit
     &     none
      integer
     &     iObj, m, n, nb, nF, nName, nnCon, nkx, nxname, nFname,
     &     Start, xstate(n), Fstate(nF), hs(nb), kxN(nkx)
      double precision
     &     infBnd, xN(n), xlow(n), xupp(n), F(nF), Flow(nF),
     &     Fupp(nF), Fmul(nF), bl(nb), bu(nb), x(nb),
     &     Fx(nF), pi(m)
      character
     &     xnames(nxname)*8, Fnames(nFname)*8, Names(nName)*8

*     ==================================================================
*     s3inA  loads re-ordered copies of the user defined values of
*     bl, bu, x, pi  and  hs into the SNOPT arrays.
*
*     25 Sep 1998: First version of s3inA.
*     28 Jun 2000: First version for snoptA.
*     03 Jun 2001: Current version.
*     ==================================================================
      integer
     &     i, j, iN, jN
*     ------------------------------------------------------------------
      integer            Warm
      parameter         (Warm  = 2)
      double precision   zero
      parameter         (zero   = 0.0d+0)
*     ------------------------------------------------------------------
      character          ColNm*1, RowNm*1
      data               ColNm /'x'/,
     &                   RowNm /'r'/
*     ------------------------------------------------------------------
*     Load the bounds on the variables first.

      do j = 1, n
         jN    = kxN(j)
         bl(j) = xlow(jN)
         bu(j) = xupp(jN)
         x (j) = xN(jN)
         hs(j) = xstate(jN)
      end do

*     Bounds on the linear and nonlinear rows.

      do i  = 1, m
         j  = n + i
         iN = kxN(j)

         if (i .le. nnCon) then
            pi(i) = Fmul(iN)
         end if

*        The linear objective row is always free.

         if (i .eq. iObj) then
            bl(j) = -infBnd
            bu(j) =  infBnd
            if (Start .eq. Warm) then
               hs(j) = 0
               x(j)  = zero
            end if
         else
            bl(j) = Flow(iN)
            bu(j) = Fupp(iN)
            if (Start .eq. Warm) then
               hs(j) = Fstate(iN)
               x (j) = F(iN)
            else
               hs(j) = 0
            end if
         end if
      end do

      call dload ( m    , zero, pi, 1 )
      if (nnCon .gt. 0)
     &call dload ( nnCon, zero, Fx, 1 )

      if (nxname .gt. 1  .or.  nFname .gt. 1) then

*        User has supplied some row or column names.

         if (nxname .gt. 1) then
            call chcopy( nxname, xnames, 1, Names, 1 )
         else

*           Provide generic variable names.

            do j = 1, n
               write(Names(j), '(a1,i7)') ColNm, j
            end do
         end if

         if (nFname .gt. 1) then
            call chcopy( nFname, Fnames, 1, Names(n+1), 1 )
         else

*           Provide generic constraint names.

            do i = 1, m
               write(Names(i), '(a1,i7)') RowNm, i
            end do
         end if
      end if

      end ! subroutine s3inA

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3outA
     &   ( n, nb, nF, nnCon, nkx, kxN, ObjRow, iObj,
     &     fObj, xstate, Fstate, hs, xN, F, x,
     &     Fx, xmul, Fmul, rc )

      implicit
     &     none
      integer
     &     iObj, n, nb, nF, nkx, nnCon, ObjRow, xstate(n),
     &     Fstate(nF), hs(nb), kxN(nkx)
      double precision
     &     fObj, xN(n), F(nF), xmul(n), Fmul(nF), x(nb), Fx(nF),
     &     rc(nb)

*     ==================================================================
*     s3outA prepares the data for final exit. Various arrays are
*     un-permuted and the working versions of x, rc  and  hs are
*     copied into their originals.
*
*     25 Sep 1998: First version of s3outA.
*     18 Jul 2000: First version for snoptA.
*     17 Jul 2005: fObj assigned to F + x(n+iObj) instead of ObjTru.
*     ==================================================================
      integer
     &     i, j, iN, jN
*     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero   = 0.0d+0)
*     ------------------------------------------------------------------

      call iload ( n , 0   , xstate, 1 )
      call dload ( n , zero, xmul  , 1 )
      call iload ( nF, 0   , Fstate, 1 )
      call dload ( nF, zero, Fmul  , 1 )

      do j = 1, n
         jN         = kxN(j)
         xN(jN)     = x(j)
         xstate(jN) = hs(j)
         xmul(jN)   = rc(j)
      end do

      do j  = n+1, nkx
         i  = j - n
         iN = kxN(j)
         if (iN .eq. ObjRow) then
            if (ObjRow .gt. 0) then
               if (iObj .gt. 0) then
                  F(ObjRow) = fObj + x(n+iObj)
               else
                  F(ObjRow) = fObj
               end if
            end if
         else
            if (i .le. nnCon) then
               F(iN) = Fx(i)
            else
               F(iN) =  x(j)
            end if
            Fstate(iN) = hs(j)
            Fmul(iN)   = rc(j)
         end if
      end do

      end ! subroutine s3outA

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3prmA
     &   ( n, nF, nkx,
     &     iGfunN, jGvarN, iGfun, jGvar, lenG, neG, kx, kxN )

      implicit
     &     none
      integer
     &     lenG, n, nF, neG, nkx, kx(nkx), kxN(nkx),
     &     iGfunN(lenG), jGvarN(lenG), iGfun(lenG), jGvar(lenG)

*     ==================================================================
*     s3prmA makes a reordered copy of the G coordinate array,
*     and inverts the row and column permutations.
*
*     29 Dec 2000: First version of s3prmA.
*     06 Jun 2001: Current version.
*     ==================================================================
      integer
     &     i, j, k, jN, Col, Row, RowN
*     ------------------------------------------------------------------
*     Assign the reordered constant Jacobian elements.
*     ------------------------------------------------------------------
      do k  = 1, neG
         Row = n + iGfunN(k)
         Col =     jGvarN(k)
         iGfun(k) = kx(Row)
         jGvar(k) = kx(Col)
      end do

*     ------------------------------------------------------------------
*     Invert the row and column lists.
*     Note: if iObj = 0, then kx(ObjRow) = 0, kxN(nF) = ObjRow.
*     ------------------------------------------------------------------
      do jN = 1, n
         j  = kx(jN)
         kxN(j) = jN
      end do

      do RowN = n+1, nkx
         i    = kx(RowN)
         if (i .gt. 0) then
            kxN(n+i  ) = RowN - n
         else ! RowN = n + ObjRow, iObj = 0
            kxN(n+nF) = RowN - n
         end if
      end do

      end ! subroutine s3prmA

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3mapA
     &   ( m, n, ne, nF, neG, negCon, nkx, nnJac, nName,
     &     nextcw, nextiw, nextrw, iw, leniw )

      implicit
     &     none
      integer
     &     leniw, m, n, ne, nF, neG, negCon, nextcw, nextiw, nextrw,
     &     nkx, nName, nnJac, iw(leniw)

*     ==================================================================
*     s3mapA   allocates additional array storage for snoptA.
*
*     29 Dec 2000: First version of s3mapA.
*     27 Dec 2002: Current version.
*     ==================================================================
      integer
     &     lbl, lbu, lF, lFmul, lG, lenNam, lhs, liGfun, lindG, lindJ,
     &     lJcol, ljGvar, lkxN, llocJ, lNames, lnGlin, lpi, lrc, lx,
     &     lxN, nb, nlocG, nlocJ
*     ------------------------------------------------------------------
      nb      = n      + m

*     Nonlinear constraints.

      nlocG   = nnJac  + 1
      nlocJ   = n      + 1

      lenNam  = 0
      if (nName .eq. nb) lenNam = nName

*     Define the addresses.

      if (lenNam .eq. 0) then
         lNames = nextcw - 1    ! cw(lNames) unused
      else
         lNames = nextcw
         nextcw = lNames + lenNam
      end if

      lhs    = nextiw
      lkxN   = lhs    + nb
      llocJ  = lkxN   + nkx
      lindJ  = llocJ  + nlocJ
      lindG  = lindJ  + ne
      lnGlin = lindG  + negCon
      liGfun = lnGlin + nlocG
      ljGvar = liGfun + neG
      nextiw = ljGvar + neG

      lJcol  = nextrw
      lbl    = lJcol  + ne
      lbu    = lbl    + nb
      lpi    = lbu    + nb
      lrc    = lpi    + m
      lx     = lrc    + nb
      lxN    = lx     + nb
      lF     = lxN    + n
      lFmul  = lF     + nF
      lG     = lFmul  + nF
      nextrw = lG     + neG

*     Store the addresses in iw.

      iw(252) = lkxN      ! jN = kxN(j ) => col j of Jcol is variable jN
      iw(256) = lJcol     ! Jcol(ne)    = Constraint Jacobian by columns
      iw(257) = llocJ     ! locJ(n+1)   = column pointers for indJ
      iw(258) = lindJ     ! indJ(ne) holds the row indices for Jij

      iw(261) = lindG     ! indG(negCon) holds the row indices for gij
      iw(262) = lnGlin    ! nGlin(j) = # linear elems in col j of gCon

      iw(266) = liGfun    ! iGfun(neG) row list of reordered G nonzeros
      iw(267) = ljGvar    ! iGvar(neG) col list of reordered G nonzeros

      iw(271) = lbl       ! bl(nb)      = lower bounds
      iw(272) = lbu       ! bu(nb)      = upper bounds
      iw(279) = lpi       ! pi(m)       = the pi-vector
      iw(280) = lrc       ! rc(nb)      = the reduced costs
      iw(282) = lhs       ! the column state vector
      iw(299) = lx        ! x(nb)       = the solution (x,s)
      iw(327) = lxN       ! xN(n)       = copy of user-defined x
      iw(328) = lF        ! F(nF)       = user-defined F
      iw(329) = lFmul     ! Fmul(nF)    = user-defined multipliers
      iw(330) = lG        ! G (lenG)    = problem derivatives

      iw(359) = lNames    ! Names(nName), row and column names

      end ! subroutine s3mapA

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3argB
     &   ( iExit, Start, m, n, ne, nName, nS,
     &     nnCon, nnObj, nnJac, iObj,
     &     indJ, locJ, bl, bu, Names,
     &     hs, pi, lvlSrt, Errors,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Errors, iExit, iObj, leniw, lenrw, lvlSrt,
     &     m, n, ne, nName, nS, nnCon,
     &     nnObj, nnJac, indJ(ne), hs(n+m), locJ(n+1),
     &     iw(leniw)
      double precision
     &     bl(n+m), bu(n+m), pi(m), rw(lenrw)
      character*(*)
     &     Start
      character
     &     Names(nName)*8

*     ==================================================================
*     s3argB   checks the arguments for snOpt.
*
*     On exit, Errors says how many errors were encountered.
*
*     21 Dec 2002: First version of s3argB.
*     03 Aug 2003: snPRNT adopted.
*     14 Sep 2003: hs = -1 accepted.
*     13 Mar 2004: Start = 'Hot FHS' options decoded.
*     26 Mar 2004: Prevent n*m overflow (noticed by Anders Goran).
*     04 Sep 2004: Length of str increased to 132.
*     04 Dec 2004: Current version of s3argB.
*     ==================================================================
      character
     &     ch1*1, str*132
      logical
     &     Named, ok
      integer
     &     argerr, Errs, idummy, ir, j, js, k, l, linRow, Mods, n1, nb,
     &     nchar
      double precision
     &     b1, b2, InfBnd
*     ------------------------------------------------------------------
      integer            COLD,       BASIS,      WARM,       HOT
      parameter         (COLD   = 0, BASIS  = 1, WARM   = 2, HOT    = 3)
      integer            gotFac,     gotHes,     gotScl
      parameter         (gotFac=230, gotHes=231, gotScl=232)
      double precision   zero
      parameter         (zero   = 0.0d+0)
      parameter         (idummy =  -11111)
*     ------------------------------------------------------------------
      character           id(3)*5
      data                id(1)   ,  id(2)   ,  id(3)
     &                 / 'varbl'  , 'lncon'  , 'nlcon'   /
*     ------------------------------------------------------------------
      InfBnd    = rw( 70) ! definition of an infinite bound
      argerr    = iw(106) ! maximum # errors in MPS data

      iExit     = 0
      Errors    = 0

*     The options haven't been checked yet.

      if (InfBnd .lt. zero) InfBnd = 1.0d+20
      if (argerr .lt. 0   ) argerr = 20 ! Print 20 lines max

      nb     = n + m
      Named  = nName .eq. nb

!     ==================================================================
!     Decode 'Start'.
!     ==================================================================
      iw(gotFac) = 0
      iw(gotHes) = 0
      iw(gotScl) = 0

!     Determine the type of start.
!     Optional parameters take precedence.

      if (lvlSrt .ne. idummy) then ! lvlSrt was set as an option
         if (    lvlSrt .eq. COLD
     &     .or.  lvlSrt .eq. BASIS
     &     .or.  lvlSrt .eq. WARM ) then
*           Relax
         else if (lvlSrt .eq. HOT) then
            iw(gotFac) = 1
            iw(gotHes) = 1
            iw(gotScl) = 1
         else                      ! lvlSrt is an unrecognized option
            lvlSrt = idummy
         end if
      end if

      if (lvlSrt .eq. idummy) then ! lvlSrt is unset
         ch1        = Start(1:1)

         if      (ch1 .eq. 'C'  .or.  ch1 .eq. 'c') then
            lvlSrt = COLD
         else if (ch1 .eq. 'B'  .or.  ch1 .eq. 'b') then
            lvlSrt = BASIS
         else if (ch1 .eq. 'W'  .or.  ch1 .eq. 'w') then
            lvlSrt = WARM
         else if (ch1 .eq. 'H'  .or.  ch1 .eq. 'h') then
            lvlSrt = HOT
!           nchar  = len_trim(Start)         ! An F90 intrinsic

            call s1trim( Start, nchar )      ! The F77 equivalent
                                             ! Decode    Start = 'HOT ...'
            if (nchar .le. 4) then           ! 'Hot' or 'Hot ' = 'Hot FHS'
               iw(gotFac) = 1
               iw(gotHes) = 1
               iw(gotScl) = 1
            else
               do j = 5, nchar               ! Decode 1 or more of FHS
                  ch1 = Start(j:j)
                  if (ch1 .eq. 'F'  .or.  ch1 .eq. 'f') iw(gotFac) = 1
                  if (ch1 .eq. 'H'  .or.  ch1 .eq. 'h') iw(gotHes) = 1
                  if (ch1 .eq. 'S'  .or.  ch1 .eq. 's') iw(gotScl) = 1
               end do
            end if
         else
            lvlSrt = COLD
            write(str, 1000) Start
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

!     ==================================================================
!     Check the other arguments.
!     ==================================================================
      if (m .lt. 1                                 ) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'm     ', m
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (n .lt. 1                                 ) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'n     ', n
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      ! 26 Mar 2004: Test if ne > n*m without overflowing.

      n1 = max(n,1)
      k  = mod(ne,n1)

      if (ne .lt. 1   .or.  (k.eq.0 .and. ne/n1 .gt. m)
     &                .or.  (k.gt.0 .and. ne/n1 .ge. m)) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'ne    ', ne
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (nnObj .lt. 0      .or.   nnObj .gt. n    ) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'nnObj ', nnObj
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (nnCon .lt. 0      .or.   nnCon .gt. m    ) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'nnCon ', nnCon
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (nnJac .lt. 0      .or.   nnJac .gt. n    ) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'nnJac ', nnJac
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (nnCon .eq. 0      .and.  nnJac .gt. 0    ) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'nnJac ', nnJac
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (nnCon .gt. 0      .and.  nnJac .eq. 0    ) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'nnJac ', nnJac
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (nName .ne. 1      .and.  nName .ne. nb   ) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'nName ', nName
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (      iObj .gt. 0  .and.  iObj  .lt. nnCon
     &    .or.  iObj .lt. 0  .or.   iObj  .gt. m   ) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'iObj  ', iObj
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (locJ(1) .ne. 1  .or.  locJ(n+1) .ne. ne+1) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1200) locJ(1), locJ(n+1)
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      do j = 1, nnJac
         linRow = 0

         do k  = locJ(j), locJ(j+1)-1
            ir = indJ(k)
            if (ir .gt. nnCon) then
*              Linear row
               linRow = ir
            else if (linRow .gt. 0) then
*              Nonlinear row occurs after a linear row.
               Errors = Errors + 1
               if (Errors .le. argerr) then
                  write(str, 1300) linRow, j, ir
                  call snPRNT( 13, str, iw, leniw )
               end if
            end if
         end do
      end do

      Mods = 0
      if (lvlSrt .eq. COLD) then
         do j  = 1, n
            js = hs(j)
            if (js .lt. -1  .or.  js .gt. 5) then
               Mods  = Mods + 1
               hs(j) = 0
            end if
         end do
      else if (lvlSrt .eq. WARM  .or.  lvlSrt .eq. HOT) then
         do j  = 1, nb
            js = hs(j)
            if (js .lt. 0  .or.  js .gt. 3) then
               Mods  = Mods + 1
               hs(j) = 0
            end if
         end do
      end if

      if (Mods .gt. 0) then
         write(str, 1400) Mods
         call snPRNT( 13, str, iw, leniw )
      end if

*     ------------------------------------------------------------------
*     Check the bounds on all variables and constraints.
*     ------------------------------------------------------------------
      Errs = 0

      do j  = 1, nb
         b1 = bl(j)
         b2 = bu(j)
         ok = b1 .lt. b2  .or.
     &        b1 .eq. b2  .and.  abs(b1) .lt. InfBnd

         if (.not. ok)  then
            Errs   = Errs   + 1
            Errors = Errors + 1

            if (j .gt. n+nnCon)  then
*              Linear    constraint
               k  = j - (n + nnCon)
               l  = 2
            else if (j .gt. n)  then
*              Nonlinear constraint

               k  = j - n
               l  = 3
            else
*              Variable
               k = j
               l = 1
            end if

            if (Errors .le. argerr) then
               if ( Named ) then
                  if (b1 .eq. b2) then
                     write(str, 1510) Names(j), b1, InfBnd
                     call snPRNT( 13, str, iw, leniw )
                  else
                     write(str, 1515) Names(j), b1, b2
                     call snPRNT( 13, str, iw, leniw )
                  end if
               else
                  if (b1 .eq. b2) then
                     write(str, 1500) id(l), k, b1, InfBnd
                     call snPRNT( 13, str, iw, leniw )
                  else
                     write(str, 1505) id(l), k, b1, b2
                     call snPRNT( 13, str, iw, leniw )
                  end if
               end if
            end if
         end if
      end do

      if (Errs .gt. 0) then
         write(str, 1600) Errs
         call snPRNT( 13, str, iw, leniw )
      end if

      if (Errors .ge. argerr) then
         call snPRNT( 3, ' and so on ...', iw, leniw )
      end if

      if (Errors .gt. 0) then
         iExit = 91             ! Invalid arguments for snOptB.
      end if

      return

 1000 format(' XXX Start parameter not recognized:  ', a)
 1100 format(' XXX  Argument out of range:  ', a6, ' = ', i6)
 1200 format(' XXX  Invalid locA(1), locA(n+1) =', 2i8)
 1300 format(' XXX  Invalid argument indA: linear row ', i6,
     &       ' in column ', i6, ' appears before nonlinear row ', i6)
 1400 format(' XXX  Invalid argument hs: ', i6,
     &       ' elements modified to be in range.')
 1500 format(' XXX  The equal bounds on  ', a5, i6,
     &       '  are infinite.   Bounds =', g16.7,
     &       '  InfBnd =', g16.7)
 1505 format(' XXX  The bounds on  ', a5, i6,
     &       '  are inconsistent.   bl =', g16.7, '   bu =', g16.7)
 1510 format(' XXX  The equal bounds on  ', a8,
     &       '  are infinite.   Bounds =', g16.7,
     &       '  InfBnd =', g16.7)
 1515 format(' XXX  The bounds on  ', a8,
     &       '  are inconsistent.   bl =', g16.7, '   bu =', g16.7)
 1600 format(' XXX  Invalid arguments bl, bu: ', i6,
     &       ' inconsistent bounds or infinite equal bounds.')

      end ! subroutine s3argB

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3prtB
     &   ( m, n, nnCon, nnJac, nnObj, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     m, n, nnCon, nnJac, nnObj, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)

*     ==================================================================
*     s3prtB prints the settings of the optional parameters for the
*     standard SNOPT wrapper.
*
*     See  snworkspace.doc  for full documentation of cw, iw and rw.
*
*     15 Nov 1991: First version (s8prnt).
*     10 Dec 2002: More LU pivoting options.
*     01 Jul 2003: QN QP options added.
*     30 Jul 2003: QN CG options added.
*     03 Aug 2003: snPRNT adopted.
*     25 Mar 2005: Current version of s3prtB.
*     ==================================================================
      integer
     &     cgItmx, iCrash, iBack, iDump, iLoadB, iInsrt, iNewB, iOldB,
     &     iPnch, iPrint, iStdi, iStdo, iSoln, iSpecs, itnlim,
     &     kchk, kdegen, kFac, klog, kreset, ksav, kSumm,
     &     lprDbg, lprPrm, lvlDer, lvlHes, lvlPiv, lvlPPm, lvlPre,
     &     lvlSch, lvlScl, lvlSrt, lvlSys, lvlTim, lvlVer,
     &     maxmn, maxR, maxS, mflush, minmax, minPrc, MjrPrt, mMajor,
     &     mMinor, mNewSB, MnrPrt, mQNmod, nnL,
     &     nParPr, nPr1, nPr2, QPslvr
      double precision
     &     bigdx, bigFx,
     &     eps, epsrf, etarg,
     &     fdint1, fdint2, scltol,
     &     tCrash, tolCon, tolCG, tolFac, tolNLP, tolpiv, tolQP,
     &     tolSwp, tolUpd, tolx, Utol1, viLim, wolfeG,
     &     wtInf0, xdlim, xPen0
      character
     &     str1*132, str2*132, str3*132, str4*132, str5*132, str6*132
*     ------------------------------------------------------------------
      character          Hestyp(1:3)*24,  lsrch(0:1)*24, pivtyp(0:3)*24,
     &                   prbtyp(1:3)*24, QPtype(0:2)*24, SrtTyp(0:3)*24,
     &                   NoYes (0:1)*3
      data               Hestyp /' Limited-Memory Hessian.',
     &                           ' Full-Memory Hessian....',
     &                           ' Exact Hessian..........'/
      data               lsrch  /' Nonderiv.  linesearch..',
     &                           ' Derivative linesearch..'/
      data               pivtyp /' LU partial  pivoting...',
     &                           ' LU rook     pivoting...',
     &                           ' LU complete pivoting...',
     &                           ' LU diagonal pivoting...'/
      data               prbtyp /' Maximize...............',
     &                           ' Feasible point only....',
     &                           ' Minimize...............'/
      data               QPtype /' QPsolver Cholesky......',
     &                           ' QPsolver CG............',
     &                           ' QPsolver QN............'/
      data               SrtTyp /' Cold start.............',
     &                           ' Basis file.............',
     &                           ' Warm start.............',
     &                           ' Hot start..............'/
      data               NoYes  /' No', 'Yes'/
*     ------------------------------------------------------------------
*     Set some local machine-dependent constants.

      eps       = rw(  1) ! machine precision.  IEEE DP  2.22e-16
      tolQP     = rw( 52) ! Minor Phase 2 Opt tol
      tolNLP    = rw( 53) ! Major Optimality tolerance
      tolCG     = rw( 54) ! cg tolerance
      tolx      = rw( 56) ! Minor feasibility tolerance.
      tolCon    = rw( 57) ! Major feasibility tolerance.
      tolpiv    = rw( 60) ! excludes small elements of y.
      tCrash    = rw( 62) ! crash tolerance.
      tolswp    = rw( 65) ! LU swap tolerance.
      tolFac    = rw( 66) ! LU factor tolerance.
      tolUpd    = rw( 67) ! LU update tolerance.
      bigFx     = rw( 71) ! unbounded objective.
      bigdx     = rw( 72) ! unbounded step.
      epsrf     = rw( 73) ! relative function precision.
      fdint1    = rw( 76) ! (1) forwrd diff. interval
      fdint2    = rw( 77) ! (2) cntrl  diff. interval
      xdlim     = rw( 80) ! Step limit
      vilim     = rw( 81) ! violation limit
      etarg     = rw( 83) ! Quasi-Newton QP rg tolerance
      wolfeG    = rw( 84) ! line search tolerance.
      wtInf0    = rw( 88) ! infeasibility weight
      xPen0     = rw( 89) ! initial penalty parameter.
      scltol    = rw( 92) ! scale tolerance.
      Utol1     = rw(154) ! abs tol for small diag of U.

      iStdi     = iw(  9) ! Standard Input
      iStdo     = iw( 10) ! Standard Output
      iSpecs    = iw( 11) ! Specs (options) file
      iPrint    = iw( 12) ! Print file
      maxR      = iw( 52) ! max columns of R.
      maxS      = iw( 53) ! max # of superbasics
      mQNmod    = iw( 54) ! (ge 0) max # of BFGS updates
      QPslvr    = iw( 55) ! 0(1) => QP(QN) QP solver
      kchk      = iw( 58) ! check (row) frequency
      kFac      = iw( 59) ! factorization frequency
      ksav      = iw( 60) ! save basis map
      klog      = iw( 61) ! log/print frequency
      kSumm     = iw( 62) ! Summary print frequency
      kDegen    = iw( 63) ! max. expansions of featol
      kReset    = iw( 64) ! Hessian frequency
      mFlush    = iw( 66) ! Hessian flush
      lvlSrt    = iw( 69) ! = 0:1:2:3 => cold:warm:basis:hot start
      lvlDer    = iw( 70) ! = 0, 1 or 2, the derivative level
      lvlSys    = iw( 71) ! > 0   => print system info
      lvlHes    = iw( 72) ! 0,1,2 => LM, FM, Exact Hessian
      lvlHes    = iw( 72) ! 0,1,2  => LM, FM, Exact Hessian
      lvlScl    = iw( 75) ! scale option
      lvlSch    = iw( 76) ! >0    => use derivatives in the line search
      lvlPre    = iw( 77) ! >0    => QN preconditioned CG
      lvlVer    = iw( 78) ! Verify level
      lvlPPm    = iw( 79) ! Proximal Point method for x0
      lvlPiv    = iw( 80) ! 0(1) LU threshold partial(complete) pivoting
      lprPrm    = iw( 81) ! > 0    => parms are printed
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
      cgItmx    = iw(111) ! CG iteration limit
      iBack     = iw(120) ! backup file
      iDump     = iw(121) ! dump file
      iLoadB    = iw(122) ! load file
      iNewB     = iw(124) ! new basis file
      iInsrt    = iw(125) ! insert file
      iOldB     = iw(126) ! old basis file
      iPnch     = iw(127) ! punch file
      iSoln     = iw(131) ! solution file
      lvlTim    = iw(182) ! Timing level
*     ------------------------------------------------------------------

      if (iPrint .le. 0 .or. MjrPrt .eq. 0 .or. lprPrm .eq. 0) return

      nnL = max( nnJac, nnObj )

      minPrc = 10
      nPr1   = n / nParPr
      nPr2   = m / nParPr
      if (max( nPr1, nPr2 ) .lt. minPrc) then
         maxmn  = max( m, n )
         nParPr = maxmn / min( maxmn, minPrc )
         nPr1   = n / nParPr
         nPr2   = m / nParPr
      end if

*     ==================================================================
*     Print parameters except if PRINT LEVEL = 0
*     or SUPPRESS PARAMETERS was specified.
*     ==================================================================
      call s1page( 1, iw, leniw )
      call snPRNT( 1, ' Parameters', iw, leniw )
      call snPRNT( 1, ' ==========', iw, leniw )

*     --------------------
*     Files.
*     --------------------
      call snPRNT(11, ' Files', iw, leniw )
      call snPRNT( 1, ' -----', iw, leniw )
      write(str1, 2110) iSoln , iOldB , iStdi
      write(str2, 2120) iInsrt, iNewB , iPrint
      write(str3, 2130) iPnch , iBack , iSpecs
      write(str4, 2140) iLoadB, iDump , iStdo
      call snPRNT( 1, str1, iw, leniw )
      call snPRNT( 1, str2, iw, leniw )
      call snPRNT( 1, str3, iw, leniw )
      call snPRNT( 1, str4, iw, leniw )

*     --------------------
*     Frequencies.
*     --------------------
      call snPRNT(11, ' Frequencies', iw, leniw )
      call snPRNT( 1, ' -----------', iw, leniw )
      write(str1, 2210) klog  , kchk  , ksav
      write(str2, 2220) kSumm , kfac  , kDegen
      call snPRNT( 1, str1, iw, leniw )
      call snPRNT( 1, str2, iw, leniw )

*     --------------------
*     QP subproblems.
*     --------------------
      call snPRNT(11, ' QP subproblems', iw, leniw )
      call snPRNT( 1, ' --------------', iw, leniw )
      write(str1, 2310) QPtype(QPslvr)
      write(str2, 2320) scltol, tolx  , itnlim
      write(str3, 2330) lvlScl, tolQP , MnrPrt
      write(str4, 2340) tCrash, tolpiv, nParPr
      write(str5, 2350) iCrash, wtInf0, nPr1
      write(str6, 2360) mNewSB, nPr2
      call snPRNT( 1, str1, iw, leniw )
      call snPRNT( 1, str2, iw, leniw )
      call snPRNT( 1, str3, iw, leniw )
      call snPRNT( 1, str4, iw, leniw )
      call snPRNT( 1, str5, iw, leniw )
      call snPRNT( 1, str6, iw, leniw )

*     -----------------------
*     CG QP solver.
*     -----------------------
      if (QPslvr .eq. 1  .or.  maxR .lt. maxS) then
         call snPRNT(11, ' Conjugate-gradient QP solver', iw, leniw )
         call snPRNT( 1, ' ----------------------------', iw, leniw )
         write(str1, 2380) etarg,   tolCG, cgItmx
         write(str2, 2390)                 lvlPre
         call snPRNT( 1, str1, iw, leniw )
         call snPRNT( 1, str2, iw, leniw )
      end if

*     --------------------
*     SQP method.
*     --------------------
      call snPRNT(11, ' The SQP Method', iw, leniw )
      call snPRNT( 1, ' --------------', iw, leniw )
      write(str1, 2410) prbtyp(2+minmax),SrtTyp(lvlSrt),lvlPPm
      write(str2, 2420) nnObj , tolNLP, epsrf
      write(str3, 2430) bigdx , maxS  , fdint1
      write(str4, 2440) bigFx , maxR  , fdint2
      call snPRNT( 1, str1, iw, leniw )
      call snPRNT( 1, str2, iw, leniw )
      call snPRNT( 1, str3, iw, leniw )
      call snPRNT( 1, str4, iw, leniw )
      write(str1, 2450) xdlim , lsrch(lvlSch) , lvlDer
      write(str2, 2460) mMajor, wolfeG, lvlVer
      write(str3, 2470) mMinor, xPen0 , MjrPrt
      call snPRNT( 1, str1, iw, leniw )
      call snPRNT( 1, str2, iw, leniw )
      call snPRNT( 1, str3, iw, leniw )

*     --------------------
*     Hessian approximation.
*     --------------------
      if (nnL .gt. 0) then
         call snPRNT(11, ' Hessian Approximation', iw, leniw )
         call snPRNT( 1, ' ---------------------', iw, leniw )
         write(str1, 2510) Hestyp(lvlHes+1), mQNmod, kReset
         write(str2, 2520) mFlush
         call snPRNT( 1, str1, iw, leniw )
         call snPRNT( 1, str2, iw, leniw )
      end if

*     --------------------
*     Nonlinear constraints.
*     --------------------
      if (nnCon .gt. 0) then
         call snPRNT(11, ' Nonlinear constraints', iw, leniw )
         call snPRNT( 1, ' ---------------------', iw, leniw )
         write(str1, 2610) nnCon , tolCon, vilim
         write(str2, 2620) nnJac
         call snPRNT( 1, str1, iw, leniw )
         call snPRNT( 1, str2, iw, leniw )
      end if

*     --------------------
*     Miscellaneous
*     --------------------
      call snPRNT(11, ' Miscellaneous', iw, leniw )
      call snPRNT( 1, ' -------------', iw, leniw )
      write(str1, 2710) tolFac, Utol1 , lvlTim
      write(str2, 2720) tolUpd, tolswp, lprDbg
      write(str3, 2730) pivtyp(lvlPiv), eps, NoYes(lvlSys)
      call snPRNT( 1, str1, iw, leniw )
      call snPRNT( 1, str2, iw, leniw )
      call snPRNT( 1, str3, iw, leniw )

      write(str1, 3000) lvlScl, nParPr
      call snPRNT( 22, str1, iw, leniw )

      return

 2110 format(
     &     ' Solution file..........', i10, 6x,
     &     ' Old basis file ........', i10, 6x,
     &     ' Standard input.........', i10)
 2120 format(
     &     ' Insert file............', i10, 6x,
     &     ' New basis file ........', i10, 6x,
     &     ' (Printer)..............', i10)
 2130 format(
     &     ' Punch file.............', i10, 6x,
     &     ' Backup basis file......', i10, 6x,
     &     ' (Specs file)...........', i10)
 2140 format(
     &     ' Load file..............', i10, 6x,
     &     ' Dump file..............', i10, 6x,
     &     ' Standard output........', i10)

 2210 format(
     &     ' Print frequency........', i10, 6x,
     &     ' Check frequency........', i10, 6x,
     &     ' Save new basis map.....', i10)
 2220 format(
     &     ' Summary frequency......', i10, 6x,
     &     ' Factorization frequency', i10, 6x,
     &     ' Expand frequency.......', i10)

 2310 format(a24)
 2320 format(
     &     ' Scale tolerance........', 0p, f10.3, 6x,
     &     ' Minor feasibility tol..', 1p, e10.2, 6x,
     &     ' Iteration limit........', i10)
 2330 format(
     &     ' Scale option...........', i10,       6x,
     &     ' Minor optimality  tol..', 1p, e10.2, 6x,
     &     ' Minor print level......', i10)
 2340 format(
     &     ' Crash tolerance........', 0p, f10.3, 6x,
     &     ' Pivot tolerance........', 1p, e10.2, 6x,
     &     ' Partial price..........', i10)
 2350 format(
     &     ' Crash option...........', i10,       6x,
     &     ' Elastic weight.........', 1p, e10.2, 6x,
     &     ' Prtl price section ( A)', i10)
 2360 format(
     &     40x,
     &     ' New superbasics........', i10,       6x,
     &     ' Prtl price section (-I)', i10)

 2380 format(
     &     ' Subspace tolerance.....', 0p, f10.5, 6x,
     &     ' CG tolerance...........', 1p, e10.2, 6x,
     &     ' CG Iterations..........', i10)

 2390 format(
     &     80x,
     &     ' CG preconditioning.....', i10)

 2410 format(a24, 16x, a24, 16x,
     &     ' Proximal Point method..', i10)
 2420 format(
     &     ' Nonlinear objectiv vars', i10,       6x,
     &     ' Major optimality tol...', 1p, e10.2, 6x,
     &     ' Function precision.....', 1p, e10.2)
 2430 format(
     &     ' Unbounded step size....', 1p, e10.2, 6x,
     &     ' Superbasics limit......', i10,       6x,
     &     ' Difference interval....', 1p, e10.2)
 2440 format(
     &     ' Unbounded objective....', 1p, e10.2, 6x,
     &     ' Reduced Hessian dim....', i10,       6x,
     &     ' Central difference int.', 1p, e10.2)
 2450 format(
     &     ' Major step limit.......', 1p, e10.2, 6x,
     &            a11,' linesearch..',           16x,
     &     ' Derivative level.......', i10)
 2460 format(
     &     ' Major iterations limit.', i10,       6x,
     &     ' Linesearch tolerance...', 0p, f10.5, 6x,
     &     ' Verify level...........', i10)
 2470 format(
     &     ' Minor iterations limit.', i10,       6x,
     &     ' Penalty parameter......', 1p, e10.2, 6x,
     &     ' Major Print Level......', i10)

 2510 format(
     &     a24,                                  16x,
     &     ' Hessian updates........', i10,       6x,
     &     ' Hessian frequency......', i10)
 2520 format(
     &     80x,
     &     ' Hessian flush..........', i10)

 2610 format(
     &     ' Nonlinear constraints..', i10,       6x,
     &     ' Major feasibility tol..', 1p, e10.2, 6x,
     &     ' Violation limit........',     e10.2)
 2620 format(
     &     ' Nonlinear Jacobian vars', i10)

 2710 format(
     &     ' LU factor tolerance....', 0p, f10.2, 6x,
     &     ' LU singularity tol.....', 1p, e10.2, 6x,
     &     ' Timing level...........', i10)
 2720 format(
     &     ' LU update tolerance....', 0p, f10.2, 6x,
     &     ' LU swap tolerance......', 1p, e10.2, 6x,
     &     ' Debug level............', i10)
 2730 format(
     &     a24,                                  16x,
     &     ' eps (machine precision)', 1p, e10.2, 6x,
     &     ' System information.....', 7x, a3 )

 3000 format(' Scale option', i3, ',    Partial price', i4)

      end ! subroutine s3prtB

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3argQ
     &   ( iExit, Start, m, n, ne, nName, nS,
     &     lencObj, iObj, ncolH,
     &     indA, locA, bl, bu, Names,
     &     hs, pi, lvlSrt, Errors,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Errors, iExit, iObj, lencObj, leniw, lenrw, lvlSrt, m, n,
     &     ncolH, ne, nName, nS, indA(ne), hs(n+m), locA(n+1),
     &     iw(leniw)
      double precision
     &     bl(n+m), bu(n+m), pi(m), rw(lenrw)
      character*(*)
     &     Start
      character
     &     Names(nName)*8

*     ==================================================================
*     s3argQ   checks the arguments for sqOpt.
*
*     On exit,
*     Errors says how many errors were encountered.
*     lvlSrt is an integer version of Start:
*        Start   lvlSrt
*        'Cold'       0
*        'Basis'      1
*        'Warm'       2
*        'Hot'        3
*     iw(gotFac,gotHes,gotScl) are set to 0 or 1.
*        'Hot F'     Keep Factors of basis (LU)
*        'Hot H'     Keep Hessian
*        'Hot S'     Keep Scales
*        'Hot FH'    Keep Factors and Hessian
*        'Hot FHS'   etc.
*
*     18 Oct 2003: First version of s3argQ.
*     08 Mar 2004: Start = 'Hot FHS' options decoded.
*     26 Mar 2004: Prevent n*m overflow (noticed by Anders Goran).
*     02 Sep 2004: Length of str increased to 132.
*     04 Dec 2004: Current version of s3argQ.
*     ==================================================================
      character
     &     ch1*1, str*132
      logical
     &     Named, ok
      integer
     &     argerr, Errs, idummy, ir, j, js, k, l, Mods, n1, nb, nchar,
     &     nz
      double precision
     &     b1, b2, InfBnd
*     ------------------------------------------------------------------
      integer            COLD,       BASIS,      WARM,       HOT
      parameter         (COLD   = 0, BASIS  = 1, WARM   = 2, HOT    = 3)
      integer            gotFac,     gotHes,     gotScl
      parameter         (gotFac=230, gotHes=231, gotScl=232)
      double precision   zero
      parameter         (zero   = 0.0d+0)
      parameter         (idummy =  -11111)
*     ------------------------------------------------------------------
      character           id(2)*5
      data                id(1)   ,  id(2)
     &                 / 'varbl'  , 'lncon'  /
*     ------------------------------------------------------------------
      InfBnd    = rw( 70) ! definition of an infinite bound
      argerr    = iw(106) ! maximum # errors in MPS data

      iExit     = 0
      Errors    = 0

*     The options haven't been checked yet.

      if (InfBnd .le. zero) InfBnd = 1.0d+20
      if (argerr .lt. 0   ) argerr = 20 ! Print 20 lines max

      nb     = n + m
      Named  = nName .eq. nb

!     ==================================================================
!     Decode 'Start'.
!     ==================================================================
      iw(gotFac) = 0
      iw(gotHes) = 0
      iw(gotScl) = 0

!     Determine the type of start.
!     Preset optional parameters take precedence.

      if (lvlSrt .ne. idummy) then ! lvlSrt was set as an option
         if (    lvlSrt .eq. COLD
     &     .or.  lvlSrt .eq. BASIS
     &     .or.  lvlSrt .eq. WARM ) then
*           Relax
         else if (lvlSrt .eq. HOT) then
            iw(gotFac) = 1
            iw(gotHes) = 1
            iw(gotScl) = 1
         else                      ! lvlSrt is an unrecognized option
            lvlSrt = idummy
         end if
      end if

      if (lvlSrt .eq. idummy) then ! lvlSrt is unset
         ch1        = Start(1:1)

         if      (ch1 .eq. 'C'  .or.  ch1 .eq. 'c') then
            lvlSrt = COLD
         else if (ch1 .eq. 'B'  .or.  ch1 .eq. 'b') then
            lvlSrt = BASIS
         else if (ch1 .eq. 'W'  .or.  ch1 .eq. 'w') then
            lvlSrt = WARM
         else if (ch1 .eq. 'H'  .or.  ch1 .eq. 'h') then
            lvlSrt = HOT
!           nchar  = len_trim(start)         ! An F90 intrinsic

            call s1trim( start, nchar )      ! The F77 equivalent
                                             ! Decode    Start = 'HOT ...'
            if (nchar .le. 4) then           ! 'Hot' or 'Hot ' = 'Hot FHS'
               iw(gotFac) = 1
               iw(gotHes) = 1
               iw(gotScl) = 1
            else
               do j = 5, nchar        ! Decode 1 or more of FHS
                  ch1 = start(j:j)
                  if (ch1 .eq. 'F'  .or.  ch1 .eq. 'f') iw(gotFac) = 1
                  if (ch1 .eq. 'H'  .or.  ch1 .eq. 'h') iw(gotHes) = 1
                  if (ch1 .eq. 'S'  .or.  ch1 .eq. 's') iw(gotScl) = 1
               end do
            end if
         else
            lvlSrt = COLD
            write(str, 1000) Start
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

!     ==================================================================
!     Check the other arguments.
!     ==================================================================
      if (m .lt. 1                                 ) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'm     ', m
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (n .lt. 1                                 ) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'n     ', n
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      ! 26 Mar 2004: Test if ne > n*m without overflowing.

      n1 = max(n,1)
      k  = mod(ne,n1)

      if (ne .lt. 1   .or.  (k.eq.0 .and. ne/n1 .gt. m)
     &                .or.  (k.gt.0 .and. ne/n1 .ge. m)) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'ne    ', ne
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (lencObj .lt. 0      .or.   lencObj .gt. n) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'lencObj', lencObj
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (ncolH .lt. 0      .or.   ncolH .gt. n    ) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'ncolH ', ncolH
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (nName .ne. 1      .and.  nName .ne. nb   ) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'nName ', nName
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (iObj .lt. 0  .or.   iObj  .gt. m   ) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'iObj  ', iObj
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (locA(1) .ne. 1  .or.  locA(n+1) .ne. ne+1) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1200) locA(1), locA(n+1)
            call snPRNT( 13, str, iw, leniw )
         end if
      else
         nz = 0
         do j = 1, n
            do k  = locA(j), locA(j+1)-1
               ir = indA(k)
               nz = nz + 1
               if (ir .gt. m  .or. ir .le. 0) then
*                 Row index out of range.
                  Errors = Errors + 1
                  if (Errors .le. argerr) then
                     write(str, 1300) ir, j
                     call snPRNT( 13, str, iw, leniw )
                  end if
               end if
            end do
         end do

         if (nz .ne. ne) then
            Errors = Errors + 1
            if (Errors .le. argerr) then
               write(str, 1100) 'ne    ', ne
               call snPRNT( 13, str, iw, leniw )
            end if
         end if
      end if

      Mods = 0
      if (lvlSrt .eq. COLD) then
         do j  = 1, n
            js = hs(j)
            if (js .lt. -1  .or.  js .gt. 5) then
               Mods  = Mods + 1
               hs(j) = 0
            end if
         end do
      else if (lvlSrt .eq. WARM  .or.  lvlSrt .eq. HOT) then
         do j  = 1, nb
            js = hs(j)
            if (js .lt. 0  .or.  js .gt. 3) then
               Mods  = Mods + 1
               hs(j) = 0
            end if
         end do
      end if

      if (Mods .gt. 0) then
         write(str, 1400) Mods
         call snPRNT( 13, str, iw, leniw )
      end if

*     ------------------------------------------------------------------
*     Check the bounds on all variables and constraints.
*     ------------------------------------------------------------------
      Errs = 0

      do j  = 1, nb
         b1 = bl(j)
         b2 = bu(j)
         ok = b1 .lt. b2  .or.
     &        b1 .eq. b2  .and.  abs(b1) .lt. InfBnd

         if (.not. ok)  then
            Errs   = Errs   + 1
            Errors = Errors + 1

            if (j .gt. n)  then
*              Linear    constraint
               k = j - n
               l = 2
            else
*              Variable
               k = j
               l = 1
            end if

            if (Errors .le. argerr) then
               if ( Named ) then
                  if (b1 .eq. b2) then
                     write(str, 1510) Names(j), b1, InfBnd
                     call snPRNT( 13, str, iw, leniw )
                  else
                     write(str, 1515) Names(j), b1, b2
                     call snPRNT( 13, str, iw, leniw )
                  end if
               else
                  if (b1 .eq. b2) then
                     write(str, 1500) id(l), k, b1, InfBnd
                     call snPRNT( 13, str, iw, leniw )
                  else
                     write(str, 1505) id(l), k, b1, b2
                     call snPRNT( 13, str, iw, leniw )
                  end if
               end if
            end if
         end if
      end do

      if (Errs .gt. 0) then
         write(str, 1600) Errs
         call snPRNT( 13, str, iw, leniw )
      end if

      if (Errors .ge. argerr) then
         call snPRNT( 3, ' and so on ...', iw, leniw )
      end if

      if (Errors .gt. 0) then
         iExit = 91             ! Invalid arguments for sQOpt.
      end if

      return

 1000 format(' XXX Start parameter not recognized:  ', a)
 1100 format(' XXX  Argument out of range:  ', a, ' = ', i6)
 1200 format(' XXX  Invalid locA(1), locA(n+1) =', 2i8)
 1300 format(' XXX  Invalid argument indA: row index ', i6,
     &       ' in column ', i6, ' is out of range ', i6)
 1400 format(' XXX  Invalid argument hs: ', i6,
     &       ' elements modified to be in range.')
 1500 format(' XXX  The equal bounds on  ', a5, i6,
     &       '  are infinite.   Bounds =', g16.7,
     &       '  InfBnd =', g16.7)
 1505 format(' XXX  The bounds on  ', a5, i6,
     &       '  are inconsistent.   bl =', g16.7, '   bu =', g16.7)
 1510 format(' XXX  The equal bounds on  ', a8,
     &       '  are infinite.   Bounds =', g16.7,
     &       '  InfBnd =', g16.7)
 1515 format(' XXX  The bounds on  ', a8,
     &       '  are inconsistent.   bl =', g16.7, '   bu =', g16.7)
 1600 format(' XXX  Invalid arguments bl, bu: ', i6,
     &       ' inconsistent bounds or infinite equal bounds.')

      end ! subroutine s3argQ

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3prtQ
     &   ( m, n, ngObj, nnH, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     m, n, ngObj, nnH, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)

*     ==================================================================
*     s3prtQ prints the settings of the optional parameters for the
*     standard SNOPT wrapper.
*
*     See  snworkspace.doc  for full documentation of cw, iw and rw.
*
*     15 Nov 1991: First version (s8prnt).
*     10 Dec 2002: More LU pivoting options.
*     01 Jul 2003: QN QP options added.
*     30 Jul 2003: QN CG options added.
*     03 Aug 2003: snPRNT adopted.
*     26 Oct 2003: Current version of s3prtQ.
*     ==================================================================
      logical
     &     QP
      character
     &     str1*132, str2*132, str3*132, str4*132, str5*132, str6*132
      integer
     &     iCrash, iBack, iDump, iLoadB, iNewB, iInsrt, iOldB,
     &     iPnch, iPrint, iStdi, iStdo, iSoln, iSpecs, itnlim,
     &     kchk, kDegen, kFac, klog, ksav, kSumm, lEmode, lprDbg,
     &     lprPrm, lvlInf, lvlPiv,
     &     lvlScl, lvlSrt, lvlSys, lvlTim, maxmn, maxR, maxS,
     &     minmax, minPrc, MnrPrt, nParPr, nPr1,
     &     nPr2, nQP, QPslvr
      double precision
     &     bigdx, eps, etarg, scltol, tCrash, tolCG, tolFac,
     &     tolpiv, tolQP,
     &     tolSwp, tolUpd, tolx, Utol1,
     &     wtInf0
*     ------------------------------------------------------------------
      character          prbtyp(1:3)*24, pivtyp(0:3)*24, QPtype(0:2)*24,
     &                   SrtTyp(0:3)*24, NoYes (0:1)*3
      data               pivtyp /' LU partial  pivoting...',
     &                           ' LU rook     pivoting...',
     &                           ' LU complete pivoting...',
     &                           ' LU diagonal pivoting...'/
      data               prbtyp /' Maximize...............',
     &                           ' Feasible point only....',
     &                           ' Minimize...............'/
      data               QPtype /' QPsolver Cholesky......',
     &                           ' QPsolver CG............',
     &                           ' QPsolver QN............'/
      data               SrtTyp /' Cold start.............',
     &                           ' Basis file.............',
     &                           ' Warm start.............',
     &                           ' Hot start..............'/
      data               NoYes  /' No', 'Yes'/
*     ------------------------------------------------------------------
*     Set some local machine-dependent constants.

      eps       = rw(  1) ! machine precision.  IEEE DP  2.22e-16
      tolQP     = rw( 52) ! Minor Phase 2 Opt tol
      tolCG     = rw( 54) ! cg tolerance
      tolx      = rw( 56) ! Minor feasibility tolerance
      tolpiv    = rw( 60) ! excludes small elements of y
      tCrash    = rw( 62) ! crash tolerance
      tolswp    = rw( 65) ! LU swap tolerance
      tolfac    = rw( 66) ! LU factor tolerance
      tolupd    = rw( 67) ! LU update tolerance
      bigdx     = rw( 72) ! unbounded step
      etarg     = rw( 83) ! Quasi-Newton QP rg tolerance
      wtInf0    = rw( 88) ! infeasibility weight
      scltol    = rw( 92) ! scale tolerance.
      Utol1     = rw(154) ! abs tol for small diag of U

      iStdi     = iw(  9) ! Standard Input
      iStdo     = iw( 10) ! Standard Output
      iSpecs    = iw( 11) ! Specs (options) file
      iPrint    = iw( 12) ! Print file
      maxR      = iw( 52) ! max columns of R.
      maxS      = iw( 53) ! max # of superbasics
      QPslvr    = iw( 55) ! = 0:1:2   => QPChol:CG:QN QP solver
      lEmode    = iw( 56) ! >0    => use elastic mode
      kchk      = iw( 58) ! check (row) frequency
      kFac      = iw( 59) ! factorization frequency
      ksav      = iw( 60) ! save basis map
      klog      = iw( 61) ! log/print frequency
      kSumm     = iw( 62) ! Summary print frequency
      kDegen    = iw( 63) ! max. expansions of featol
      lvlSrt    = iw( 69) ! = 0:1:2:3 => cold:warm:basis:hot start
      lvlSys    = iw( 71) ! > 0   => print system info
      lvlInf    = iw( 73) ! Elastic option
      lvlScl    = iw( 75) ! scale option
      lvlPiv    = iw( 80) ! 0(1) LU threshold partial(complete) pivoting
      lprPrm    = iw( 81) ! > 0    => parms are printed
      lprDbg    = iw( 85) ! > 0    => private debug print
      minmax    = iw( 87) ! 1, 0, -1  => MIN, FP, MAX
      iCrash    = iw( 88) ! Crash option
      itnlim    = iw( 89) ! limit on total iterations
      MnrPrt    = iw( 93) ! Minor print level
      nParPr    = iw( 94) ! # of partial pricing sections
      iBack     = iw(120) ! backup file
      iDump     = iw(121) ! dump file
      iLoadB    = iw(122) ! load file
      iNewB     = iw(124) ! new basis file
      iInsrt    = iw(125) ! insert file
      iOldB     = iw(126) ! old basis file
      iPnch     = iw(127) ! punch file
      iSoln     = iw(131) ! Solution file
      lvlTim    = iw(182) ! Timing level
*     ------------------------------------------------------------------

      if (lprPrm .eq. 0  .or.  MnrPrt .eq. 0) return

      QP     = nnH .gt. 0
      nQP    = max( ngObj, nnH )

      minPrc = 10
      nPr1   = n / nParPr
      nPr2   = m / nParPr
      if (max( nPr1, nPr2 ) .lt. minPrc) then
         maxmn  = max( m, n )
         nParPr = maxmn / min( maxmn, minPrc )
         nPr1   = n / nParPr
         nPr2   = m / nParPr
      end if

*     ==================================================================
*     Print parameters except if PRINT LEVEL = 0
*     or SUPPRESS PARAMETERS was specified.
*     ==================================================================

      call s1page( 1, iw, leniw )
      call snPRNT( 1, ' Parameters', iw, leniw )
      call snPRNT( 1, ' ==========', iw, leniw )

*     --------------------
*     Files.
*     --------------------
      call snPRNT(11, ' Files', iw, leniw )
      call snPRNT( 1, ' -----', iw, leniw )
      write(str1, 2110) iSoln , iOldB , iStdi
      write(str2, 2120) iInsrt, iNewB , iPrint
      write(str3, 2130) iPnch , iBack , iSpecs
      write(str4, 2140) iLoadB, iDump , iStdo
      call snPRNT( 1, str1, iw, leniw )
      call snPRNT( 1, str2, iw, leniw )
      call snPRNT( 1, str3, iw, leniw )
      call snPRNT( 1, str4, iw, leniw )

*     --------------------
*     Frequencies.
*     --------------------
      call snPRNT(11, ' Frequencies', iw, leniw )
      call snPRNT( 1, ' -----------', iw, leniw )
      write(str1, 2210) klog  , kchk  , ksav
      write(str2, 2220) kSumm , kfac  , kDegen
      call snPRNT( 1, str1, iw, leniw )
      call snPRNT( 1, str2, iw, leniw )

*     --------------------
*     LP/QP parameters.
*     --------------------
      call snPRNT(11, ' LP/QP Parameters', iw, leniw )
      call snPRNT( 1, ' ----------------', iw, leniw )
      write(str1, 2310) prbtyp(2+minmax), QPtype(QPslvr), SrtTyp(lvlSrt)
      write(str2, 2320) scltol, tolx  , itnlim
      write(str3, 2330) lvlScl, tolQP , MnrPrt
      write(str4, 2340) tCrash, tolpiv, nParPr
      write(str5, 2350) iCrash, wtInf0, nPr1
      write(str6, 2360) lEmode, lvlInf, nPr2

      call snPRNT( 1, str1, iw, leniw )
      call snPRNT( 1, str2, iw, leniw )
      call snPRNT( 1, str3, iw, leniw )
      call snPRNT( 1, str4, iw, leniw )
      call snPRNT( 1, str5, iw, leniw )
      call snPRNT( 1, str6, iw, leniw )

*     --------------------
*     QP objective
*     --------------------
      if ( QP ) then
         call snPRNT(11, ' QP objective', iw, leniw )
         call snPRNT( 1, ' ------------', iw, leniw )
         write(str1, 2410) nQP   , nnH   , maxS
         write(str2, 2420) nnH   , bigdx
         write(str3, 2430) ngObj
         call snPRNT( 1, str1, iw, leniw )
         call snPRNT( 1, str2, iw, leniw )
         call snPRNT( 1, str3, iw, leniw )
      end if

*     -----------------------
*     Quasi-Newton QP solver.
*     -----------------------
      if (QPslvr .eq. 1  .or.  maxR .lt. maxS) then
         call snPRNT(11, ' Conjugate-gradient QP solver', iw, leniw )
         call snPRNT( 1, ' ----------------------------', iw, leniw )
         write(str1, 2380) etarg,   tolCG
         call snPRNT( 1, str1, iw, leniw )
      end if

*     --------------------
*     Miscellaneous
*     --------------------
      call snPRNT(11, ' Miscellaneous', iw, leniw )
      call snPRNT( 1, ' -------------', iw, leniw )
      write(str1, 2710) tolFac, Utol1 , lvlTim
      write(str2, 2720) tolUpd, tolswp, lprDbg
      write(str3, 2730) pivtyp(lvlPiv), eps, NoYes(lvlSys)
      call snPRNT( 1, str1, iw, leniw )
      call snPRNT( 1, str2, iw, leniw )
      call snPRNT( 1, str3, iw, leniw )

      return

 2110 format(
     &     ' Solution file..........', i10, 6x,
     &     ' Old basis file ........', i10, 6x,
     &     ' Standard input.........', i10)
 2120 format(
     &     ' Insert file............', i10, 6x,
     &     ' New basis file ........', i10, 6x,
     &     ' (Printer)..............', i10)
 2130 format(
     &     ' Punch file.............', i10, 6x,
     &     ' Backup basis file......', i10, 6x,
     &     ' (Specs file)...........', i10)
 2140 format(
     &     ' Load file..............', i10, 6x,
     &     ' Dump file..............', i10, 6x,
     &     ' Standard output........', i10)

 2210 format(
     &     ' Print frequency........', i10, 6x,
     &     ' Check frequency........', i10, 6x,
     &     ' Save new basis map.....', i10)
 2220 format(
     &     ' Summary frequency......', i10, 6x,
     &     ' Factorization frequency', i10, 6x,
     &     ' Expand frequency.......', i10)

 2310 format(a24, 16x, a24, 16x, a24)
 2320 format(
     &     ' Scale tolerance........', 0p, f10.3, 6x,
     &     ' Feasibility tolerance..', 1p, e10.2, 6x,
     &     ' Iteration limit........', i10)
 2330 format(
     &     ' Scale option...........', i10,       6x,
     &     ' Optimality tolerance...', 1p, e10.2, 6x,
     &     ' Print level............', i10)
 2340 format(
     &     ' Crash tolerance........', 0p, f10.3, 6x,
     &     ' Pivot tolerance........', 1p, e10.2, 6x,
     &     ' Partial price..........', i10)
 2350 format(
     &     ' Crash option...........', i10,       6x,
     &     ' Elastic weight.........', 1p, e10.2, 6x,
     &     ' Prtl price section ( A)', i10)
 2360 format(
     &     ' Elastic mode...........', i10,       6x,
     &     ' Elastic objective......', i10,       6x,
     &     ' Prtl price section (-I)', i10)

 2380 format(
     &     ' Subspace tolerance.....', 0p, f10.5, 6x,
     &     ' CG tolerance...........', 1p, e10.2)

 2410 format(
     &     ' Objective variables....', i10,       6x,
     &     ' Hessian columns........', i10,       6x,
     &     ' Superbasics limit......', i10)
 2420 format(
     &     ' Nonlin Objective vars..', i10,       6x,
     &     ' Unbounded step size....', 1p, e10.2)
 2430 format(
     &     ' Linear Objective vars..', i10)

 2710 format(
     &     ' LU factor tolerance....', 0p, f10.2, 6x,
     &     ' LU singularity tol.....', 1p, e10.2, 6x,
     &     ' Timing level...........', i10)
 2720 format(
     &     ' LU update tolerance....', 0p, f10.2, 6x,
     &     ' LU swap tolerance......', 1p, e10.2, 6x,
     &     ' Debug level............', i10)
 2730 format(
     &     a24,                                  16x,
     &     ' eps (machine precision)', 1p, e10.2, 6x,
     &     ' System information.....', 7x, a3 )

      end ! subroutine s3prtQ

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3argN
     &   ( iExit, Start, ldA, ldcJ, ldH,
     &     n, nclin, ncnln, nName,
     &     bl, bu, Names, iState, cMul, lvlSrt, Errors,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Errors, iExit, ldA, ldcJ, ldH, leniw, lenrw, lvlSrt,
     &     n, nclin, ncnln, nName, iState(n+nclin+ncnln),
     &     iw(leniw)
      double precision
     &     bl(n+nclin+ncnln), bu(n+nclin+ncnln), cMul(n+nclin+ncnln),
     &     rw(lenrw)
      character*(*)
     &     Start
      character
     &     Names(nName)*8

*     ==================================================================
*     s3argN   checks the arguments for npOpt.
*
*     On exit, Errors says how many errors were encountered.
*
*     17 Mar 2002: First version of s3argN.
*     17 Jun 2004: Current version of s3argN.
*     ==================================================================
      character
     &     ch1*1, str*132
      logical
     &     Named, ok
      integer
     &     argerr, Errs, idummy, iS, i, j, k, l, Mods, nb, nchar
      double precision
     &     b1, b2, InfBnd, Mul
*     ------------------------------------------------------------------
      integer            COLD,       BASIS,      WARM,       HOT
      parameter         (COLD   = 0, BASIS  = 1, WARM   = 2, HOT    = 3)
      integer            gotFac,     gotHes,     gotScl
      parameter         (gotFac=230, gotHes=231, gotScl=232)
      double precision   zero
      parameter         (zero   = 0.0d+0)
      parameter         (idummy =  -11111)
*     ------------------------------------------------------------------
      character           id(3)*5
      data                id(1)   ,  id(2)   ,  id(3)
     &                 / 'varbl'  , 'lncon'  , 'nlcon'   /
*     ------------------------------------------------------------------
      InfBnd    = rw( 70) ! definition of an infinite bound
      argerr    = iw(106) ! maximum # errors in MPS data

      iExit     = 0
      Errors    = 0

*     The options haven't been checked yet.

      if (InfBnd .lt. zero) InfBnd = 1.0d+20
      if (argerr .lt. 0   ) argerr = 20 ! Print 20 lines max

      nb     = n + nclin + ncnln
      Named  = nName .eq. nb

!     ==================================================================
!     Decode 'Start'.
!     ==================================================================
      iw(gotFac) = 0
      iw(gotHes) = 0
      iw(gotScl) = 0

!     Determine the type of start.
!     Optional parameters take precedence.

      if (lvlSrt .ne. idummy) then ! lvlSrt was set as an option
         if (    lvlSrt .eq. COLD
     &     .or.  lvlSrt .eq. BASIS
     &     .or.  lvlSrt .eq. WARM ) then
*           Relax
         else if (lvlSrt .eq. HOT) then
            iw(gotFac) = 1
            iw(gotHes) = 1
            iw(gotScl) = 1
         else                      ! lvlSrt is an unrecognized option
            lvlSrt = idummy
         end if
      end if

      if (lvlSrt .eq. idummy) then ! lvlSrt is unset
         ch1        = Start(1:1)

         if      (ch1 .eq. 'C'  .or.  ch1 .eq. 'c') then
            lvlSrt = COLD
         else if (ch1 .eq. 'B'  .or.  ch1 .eq. 'b') then
            lvlSrt = BASIS
         else if (ch1 .eq. 'W'  .or.  ch1 .eq. 'w') then
            lvlSrt = WARM
         else if (ch1 .eq. 'H'  .or.  ch1 .eq. 'h') then
            lvlSrt = HOT
!           nchar  = len_trim(Start)         ! An F90 intrinsic

            call s1trim( Start, nchar )      ! The F77 equivalent
                                             ! Decode    Start = 'HOT ...'
            if (nchar .le. 4) then           ! 'Hot' or 'Hot ' = 'Hot FHS'
               iw(gotFac) = 1
               iw(gotHes) = 1
               iw(gotScl) = 1
            else
               do j = 5, nchar               ! Decode 1 or more of FHS
                  ch1 = Start(j:j)
                  if (ch1 .eq. 'F'  .or.  ch1 .eq. 'f') iw(gotFac) = 1
                  if (ch1 .eq. 'H'  .or.  ch1 .eq. 'h') iw(gotHes) = 1
                  if (ch1 .eq. 'S'  .or.  ch1 .eq. 's') iw(gotScl) = 1
               end do
            end if
         else
            lvlSrt = COLD
            write(str, 1000) Start
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

!     ==================================================================
!     Check the other arguments.
!     ==================================================================
      if (n .lt. 1                                 ) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'n     ', n
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (nclin .lt. 0                             ) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'nclin ', nclin
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (ncnln .lt. 0                             ) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'ncnln ', ncnln
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (ldA .lt. 0      .or.   nclin .gt. ldA  ) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'ldA   ', ldA
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (ldcJ.lt. 0      .or.   ncnln .gt. ldcJ ) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'ldcJ  ', ldcJ
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (ldH .lt. 0      .or.   n     .gt. ldH  ) then
         Errors = Errors + 1
         if (Errors .le. argerr) then
            write(str, 1100) 'ldH   ', ldH
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      Mods = 0

      if (lvlSrt .eq. WARM  .or.  lvlSrt .eq. HOT) then
         do j  = 1, nb
            iS = iState(j)
            if (iS .lt. -2  .or.  iS .gt. 4) then
               Mods      = Mods + 1
               iState(j) = 0
            end if
         end do

         do i   = 1, ncnln
            j   = n + nclin + i
            iS  = iState(j)
            Mul = cMul(j)

            if      (iS .eq. 0) then
               Mul = zero

            else if (iS .eq. 1) then
               if (bl(j) .le. -infBnd) iS   = 0
               if (Mul   .lt.  zero  .or.  iS .eq. 0) Mul = zero

            else if (iS .eq. 2) then
               if (bu(j) .ge.  InfBnd) iS   = 0
               if (Mul   .gt.  zero  .or.  iS .eq. 0) Mul = zero

            else if (iS .eq. 3) then
               if (bl(j) .lt.   bu(j)) iS   = 0
            end if

            iState(j) = iS
            cMul(j)   = Mul
         end do
      end if

      if (Mods .gt. 0) then
         write(str, 1400) Mods
         call snPRNT( 13, str, iw, leniw )
      end if

*     ------------------------------------------------------------------
*     Check the bounds on all variables and constraints.
*     ------------------------------------------------------------------
      Errs = 0

      do j  = 1, nb
         b1 = bl(j)
         b2 = bu(j)
         ok = b1 .lt. b2  .or.
     &        b1 .eq. b2  .and.  abs(b1) .lt. InfBnd

         if (.not. ok)  then
            Errs   = Errs   + 1
            Errors = Errors + 1

           if (j .gt. n+nclin)  then
               k  = j - n - nclin
               l  = 3
            else if (j .gt. n)  then
               k  = j - n
               l  = 2
            else
               k = j
               l = 1
            end if

            if (Errors .le. argerr) then
               if ( Named ) then
                  if (b1 .eq. b2) then
                     write(str, 1510) Names(j), b1, InfBnd
                     call snPRNT( 13, str, iw, leniw )
                  else
                     write(str, 1515) Names(j), b1, b2
                     call snPRNT( 13, str, iw, leniw )
                  end if
               else
                  if (b1 .eq. b2) then
                     write(str, 1500) id(l), k, b1, InfBnd
                     call snPRNT( 13, str, iw, leniw )
                  else
                     write(str, 1505) id(l), k, b1, b2
                     call snPRNT( 13, str, iw, leniw )
                  end if
               end if
            end if
         end if
      end do

      if (Errs .gt. 0) then
         write(str, 1600) Errs
         call snPRNT( 13, str, iw, leniw )
      end if

      if (Errors .ge. argerr) then
         call snPRNT( 3, ' and so on ...', iw, leniw )
      end if

      if (Errors .gt. 0) then
         iExit = 91             ! Invalid arguments for npOptN.
      end if

      return

 1000 format(' XXX Start parameter not recognized:  ', a)
 1100 format(' XXX  Argument out of range:  ', a6, ' = ', i6)
 1400 format(' XXX  Invalid argument iState: ', i6,
     &       ' elements modified to be in range.')
 1500 format(' XXX  The equal bounds on  ', a5, i6,
     &       '  are infinite.   Bounds =', g16.7,
     &       '  InfBnd =', g16.7)
 1505 format(' XXX  The bounds on  ', a5, i6,
     &       '  are inconsistent.   bl =', g16.7, '   bu =', g16.7)
 1510 format(' XXX  The equal bounds on  ', a8,
     &       '  are infinite.   Bounds =', g16.7,
     &       '  InfBnd =', g16.7)
 1515 format(' XXX  The bounds on  ', a8,
     &       '  are inconsistent.   bl =', g16.7, '   bu =', g16.7)
 1600 format(' XXX  Invalid arguments bl, bu: ', i6,
     &       ' inconsistent bounds or infinite equal bounds.')

      end ! subroutine s3argN

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3HesN
     &   ( task, ldH, lenH, n, H, Hess )

      implicit
     &     none
      integer
     &     ldH, lenH, n, task
      double precision
     &     H(lenH), Hess(ldH,*)

*     ==================================================================
*     s3HesN loads the problem into SNOPT format.
*
*     07 Jul 1998: First version of s3HesN.
*     04 Nov 2000: Current version.
*     ==================================================================
      integer
     &     i, j, l
*     ------------------------------------------------------------------
      integer            Load,       UnLoad
      parameter         (Load   = 0, UnLoad = 1)
*     ------------------------------------------------------------------

      if (task .eq. Load) then
*        ---------------------------------------------------------------
*        Load the user-supplied Hessian Hess into H.
*        ---------------------------------------------------------------
         l = 0
         do i = 1, n
            do j = i, n
               l = l + 1
               H(l) = Hess(i,j)
            end do
         end do

      else if (task .eq. UnLoad) then
*        ---------------------------------------------------------------
*        Down load the SNOPT approximate Hessian into Hess.
*        ---------------------------------------------------------------
         l = 0
         do i = 1, n
            do j = i, n
               l = l + 1
               Hess(i,j) = H(l)
               Hess(j,i) = H(l)
            end do
         end do
      end if

      end ! subroutine s3HesN

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3iniN
     &   ( Start, n, nb, nnCon0, nnCon, negCon, hs,
     &     fCon, gCon, gObj, rc, x,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     leniw, lenrw, n, nb, negCon, nnCon0, nnCon, Start, hs(nb),
     &     iw(leniw)
      double precision
     &     fCon(nnCon0), gCon(negCon), gObj(n), rc(nb), x(nb), rw(lenrw)

*     ==================================================================
*     s3iniN initializes the SNOPT variables that are ultimately copied
*     to NPSOL format in s3outN.
*
*     04 Dec 2004: First version of s3iniN.
*     04 Dec 2004: Current version of s3iniN.
*     ==================================================================
      integer
     &     i, j, lenfH, lfH, lvlHes
*     ------------------------------------------------------------------
      integer            Cold
      parameter         (Cold   = 0)
      integer            FM
      parameter         (FM     = 1)
      double precision   zero
      parameter         (zero   = 0.0d+0)
*     ------------------------------------------------------------------
      lvlHes    = iw( 72)       ! 0,1,2  => LM, FM, Exact Hessian
      lfH       = iw(391)       ! H(lenfH), full-memory BFGS Hessian
      lenfH     = iw(392)       !

      if (Start .eq. Cold) then
         do j = 1, nb
            hs(j) = 0
            x (j) = zero
            rc(j) = zero
         end do

         do j = 1, n
            gObj(j) = zero
         end do

         do i = 1, nnCon
            fCon(i) = zero
         end do

         do i = 1, nnCon*n
            gCon(i) = zero
         end do

         if (lvlHes .eq. FM) then
            call dload ( lenfH, zero, rw(lfH), 1 )
         end if
      end if

      end ! subroutine s3iniN

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3inN
     &   ( Start, ldA, ldH, m, n, ncLin, nCon, nnCol,
     &     nb, nnCon0, nnCon, hs, iState,
     &     Alin, ne, nlocJ, locJ, indJ, Jcol,
     &     bl, bu, bbl, bbu, c, cMul,
     &     Hess, pi, x, xs,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     ldA, ldH, m, n, ncLin, nCon, nnCol, Start, nb,
     &     nnCon0, nnCon, ne, nlocJ, leniw, lenrw, indJ(ne), hs(nb),
     &     iState(n+nCon), locJ(nlocJ), iw(leniw)
      double precision
     &     Alin(ldA,*), bl(n+nCon), bu(n+nCon), bbl(nb), bbu(nb),
     &     c(nnCon0), cMul(n+nCon), Hess(ldH,*), Jcol(ne),
     &     pi(nnCon0), x(n), xs(nb), rw(lenrw)

*     ==================================================================
*     s3inN loads a problem in NPSOL format into SNOPT format.
*
*     22 Mar 1997: First version of s3inN.
*     04 Jan 2001: Current version.
*     ==================================================================
      integer
     &     i, iJ, is, j, js, l, lenfH, lfH, lH0, lvlHes
      double precision
     &     infBnd, xj
*     ------------------------------------------------------------------
      integer            Load
      parameter         (Load   = 0)
      integer            Cold,       Warm
      parameter         (Cold   = 0, Warm   = 2)
      integer            LM   ,      FM
      parameter         (LM     = 0, FM     = 1)
      double precision   zero
      parameter         (zero   = 0.0d+0)
*     ------------------------------------------------------------------
      infBnd    = rw( 70) ! definition of an infinite bound
      lvlHes    = iw( 72) ! 0,1,2  => LM, FM, Exact Hessian
      lH0       = iw(346) ! Initial diagonal Hessian

*     ==================================================================
*     Load the snopt arrays.
*     Copy the bounds, x's first, linears next, then nonlinears.
*     ==================================================================
      call dcopy ( n, bl, 1, bbl, 1 )
      call dcopy ( n, bu, 1, bbu, 1 )

      if (ncLin .gt. 0) then
         call dcopy ( ncLin, bl(n+1), 1, bbl(n+nnCon+1), 1 )
         call dcopy ( ncLin, bu(n+1), 1, bbu(n+nnCon+1), 1 )
      end if

      if (nnCon .gt. 0) then
         call dcopy ( nnCon, bl(n+ncLin+1), 1, bbl(n+1), 1 )
         call dcopy ( nnCon, bu(n+ncLin+1), 1, bbu(n+1), 1 )
      end if

      if (nCon .eq. 0) then
         bbl(nb) = - infBnd
         bbu(nb) =   infBnd
      end if

      if (Start .eq. Cold) then
*        --------------------------------------------
*        Cold Start.
*        --------------------------------------------
         do  j = 1, n
            xj = x(j)
            if (xj .le. bl(j)) then
               hs(j) = 4
            else if (xj .ge. bu(j)) then
               hs(j) = 5
            else
               hs(j) = 0
            end if
            xs(j) = xj
         end do

         if (nnCon .gt. 0) then
            call dload ( nnCon, (zero), pi, 1 )
         end if

      else if (Start .eq. Warm) then
*        ----------------------------------------------
*        Warm Start.
*        Input values of x, cMul, hs and Hess are used.
*        Note: the use of Hess is unique to snoptn.
*        ----------------------------------------------
         call dcopy ( n, x, 1, xs, 1 )

         if (nnCon .gt. 0) then
            call dcopy ( nnCon, c              , 1, xs(n+1), 1 )
            call dcopy ( nnCon, cMul(n+ncLin+1), 1, pi     , 1 )
         end if

         if (ncLin .gt. 0) then
            call dload ( ncLin, (zero), xs(n+nnCon+1), 1 )
            do j = 1, n
               call daxpy( ncLin, xs(j), Alin(1,j), 1, xs(n+nnCon+1), 1)
            end do
         end if

         l = 1
         do  j = 1, n+nCon
            is = istate(j)
            js = 0
            if (is .eq. 1) then
               js = 0
            else if (is .eq. 2) then
               js = 1
            end if

            hs(l) = js

            if (j .eq. n  .and.  ncLin .gt. 0) then
               l = n + nnCon
            else if (j .eq. n+ncLin) then
               l = n + 1
            else
               l = l + 1
            end if
         end do

         if (lvlHes .eq. LM) then
            call dcopy ( n, Hess(1,1), (ldH+1), rw(lH0), 1 )

         else if (lvlHes .eq. FM) then
            lfH   = iw(391)     ! H(lenfH), full-memory BFGS Hessian
            lenfH = iw(392)     !
            call s3HesN
     &         ( Load, ldH, lenfH, n, rw(lfH), Hess )
         end if
      end if ! cold start

*     ---------------------------------------------------------------
*     Load the linear part of A with the linear constraints.
*     ---------------------------------------------------------------
      if (nnCol .eq. 0) then

*        Sparse dummy row

         Jcol(1) = zero
         indJ(1) = 1
         locJ(1) = 1

         do j = 2, n+1
            locJ(j) = 2
         end do

      else
         iJ      = 1
         locJ(1) = 1
         do j = 1, n
            do i = 1, m
               indJ(iJ) = i
               if (i .le. nnCon) then
                  Jcol(iJ) = zero
               else if (i .le. nCon) then
                  Jcol(iJ) = Alin(i-nnCon,j)
               end if
               iJ = iJ + 1
            end do
            locJ(j+1) = iJ
         end do
      end if

      end ! subroutine s3inN

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3outN
     &   ( ldcJ, ldH, n, ncLin, nCon,
     &     nb, nnCon0, nnCon, hs, iState, c, cJac, cMul,
     &     fCon, gCon, gObj, grad,
     &     Hess, rc, x, xs,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     ldcJ, ldH, n, ncLin, nCon, nb, nnCon0, nnCon, leniw, lenrw,
     &     hs(nb), iState(n+nCon), iw(leniw)
      double precision
     &     c(nnCon0), cJac(ldcJ,*), cMul(n+nCon), fCon(nnCon0),
     &     gCon(nnCon0,*), grad(n), gObj(n), Hess(ldH,*),
     &     rc(nb), x(n), xs(nb), rw(lenrw)

*     ==================================================================
*     s3outN changes the problem from SNOPT to NPSOL format.
*
*     22 Mar 1997: First version of s3outN.
*     04 Jan 2001: Current version.
*     ==================================================================
      integer
     &     i, is, j, js, l, lenfH, lfH, lvlHes
*     ------------------------------------------------------------------
      integer            UnLoad
      parameter         (UnLoad = 1)
      integer            FM
      parameter         (FM     = 1)
*     ------------------------------------------------------------------
      lvlHes    = iw( 72) ! 0,1,2  => LM, FM, Exact Hessian

*     ==================================================================
*     Unload the SNOPT solution into the snoptn arrays
*     Copy gCon, gObj into cJac and grad,
*     ==================================================================
      l = 1
      do j = 1, n+nCon
         js = hs(l)
         is = 0
         if (js .eq. 0) then
            is = 1
         else if (js .eq. 1) then
            is = 2
         end if

         istate(j) = is

         if (j  .eq. n  .and.  ncLin .gt. 0) then
            l = n + nnCon
         else if (j .eq. n+ncLin) then
            l = n + 1
         else
            l = l + 1
         end if
      end do

*     ------------------------------------------------------------------
*     Copy gCon, gObj into cJac and grad,
*     ------------------------------------------------------------------
      call dcopy ( n, xs  , 1, x   , 1 )
      call dcopy ( n, gObj, 1, grad, 1 )
      call dcopy ( n, rc  , 1, cMul, 1 )
      if (ncLin .gt. 0)
     &   call dcopy ( ncLin, rc(n+nnCon+1), 1, cMul(n+1), 1 )

      if (nnCon .gt. 0) then
         call dcopy ( nnCon, fCon   , 1, c              , 1 )
         call dcopy ( nnCon, rc(n+1), 1, cMul(n+ncLin+1), 1 )

         do j = 1, n
            do i = 1, nnCon
               cJac(i,j) = gCon(i,j)
            end do
         end do
      end if

      if (lvlHes .eq. FM) then
         lfH   = iw(391)        ! H(lenfH), full-memory BFGS Hessian
         lenfH = iw(392)        !
         call s3HesN
     &      ( UnLoad, ldH, lenfH, n, rw(lfH), Hess )
      end if

      end ! subroutine s3outN

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3prtN
     &   ( n, nb, ncLin, nnCon0, ldA, lprSol, xNorm,
     &     iState, A, bl, bu, c, cMul, x, r,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     n, nb, ncLin, nnCon0, ldA, lprSol, leniw, lenrw,
     &     iState(nb), iw(leniw)
      double precision
     &     xNorm, A(ldA,*), c(nnCon0),
     &     bl(nb), bu(nb), cMul(nb), r(nb), x(n),
     &     rw(lenrw)

*     ==================================================================
*     s3prtN  prints  x,  A*x, c(x), the bounds, the
*     multipliers, and the slacks (distance to the nearer bound).
*
*     22 Mar 1997: First version of s3prtN.
*     03 Aug 2003: snPRNT adopted.
*     22 Nov 2003: Current version of s3prtN.
*     ==================================================================
      external
     &     ddot
      integer
     &     iPrint, i, is, j, nplin, number
      double precision
     &     b1, b2, InfBnd, ddot, rj, slk, slk1, slk2, tol, tolx, wlam
      character
     &     key*1, state*2, name*8, line*132, str*132
*     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero  = 0.0d+0)
*     ------------------------------------------------------------------
      character          lstate(-2:4)*2
      data               lstate(-2) / '--' /, lstate(-1) / '++' /
      data               lstate( 0) / 'FR' /, lstate( 1) / 'LL' /
      data               lstate( 2) / 'UL' /, lstate( 3) / 'EQ' /
      data               lstate( 4) / 'TF' /
*     ------------------------------------------------------------------
      iPrint = iw( 12) ! Print file
      tolx   = rw( 56) ! Minor feasibility tolerance.
      InfBnd = rw( 70) ! definition of an infinite bound

      if (iPrint .eq. 0  .or.  lprSol .eq. 0) return

      nplin  = n     + ncLin
      tol    = tolx

      write(str, 1000) 'Variable       '
      call snPRNT( 11, ' ', iw, leniw )
      call snPRNT(  1, str, iw, leniw )
      call snPRNT(  1, ' ', iw, leniw )

      name   = 'variable'
      nplin  = n + ncLin

      do j = 1, nb
         b1   = bl(j)
         b2   = bu(j)
         wlam = cMul(j)

         if (j .le. n) then
            rj  = x(j)
            tol = tolx
         else
            tol = tolx*xNorm

            if (j .le. nplin) then
               i  = j - n
               rj = ddot  ( n, A(i,1), ldA, x, 1 )
            else
               i  = j - nplin
               rj = c(i)
            end if
         end if

         slk1 = rj - b1
         slk2 = b2 - rj
         r(j) = rj

*        Reset istate if necessary.

         is   = istate(j)
         if (                  slk1 .lt. -tol) is = - 2
         if (                  slk2 .lt. -tol) is = - 1
         if (is .eq. 1  .and.  slk1 .gt.  tol) is =   0
         if (is .eq. 2  .and.  slk2 .gt.  tol) is =   0
         istate(j) = is
         state     = lstate(is)


         if (j .le. n) then
            number = j
         else if (j .le. nplin) then
            number = j - n
            if (number .eq. 1) then
               write(str, 1000) 'Linear constrnt'
               call snPRNT( 11, ' ', iw, leniw )
               call snPRNT(  1, str, iw, leniw )
               call snPRNT(  1, ' ', iw, leniw )
               name = 'lincon  '
            end if
         else
            number = j - nplin
            if (number .eq. 1) then
               write(str, 1000) 'Nonlin constrnt'
               call snPRNT( 11, ' ', iw, leniw )
               call snPRNT(  1, str, iw, leniw )
               call snPRNT(  1, ' ', iw, leniw )
               name = 'nlncon  '
            end if
         end if

*        ------------------------------------------------
*        Print a line for the jth variable or constraint.
*        ------------------------------------------------
         if (abs(slk1) .lt. abs(slk2)) then
            slk = slk1
            if (b1 .le. - InfBnd) slk = slk2
         else
            slk = slk2
            if (b2 .ge.   InfBnd) slk = slk1
         end if

*        Flag infeasibilities, primal and dual degeneracies,
*        and active QP constraints that are loose in NP.
*
         key    = ' '
         if (slk1 .lt. -tol  .or.       slk2  .lt. -tol) key = 'I'
         if (is   .eq.  0    .and.  abs(slk ) .le.  tol) key = 'D'
         if (is   .ge.  1    .and.  abs(wlam) .le.  tol) key = 'A'

         write(line, 2000) name, number, key, state,
     &                     rj, b1, b2, wlam, slk

*        Reset special cases:
*           Infinite bounds
*           Zero bounds
*           Lagrange multipliers for inactive constraints
*           Lagrange multipliers for infinite bounds
*           Infinite slacks
*           Zero slacks

         if (b1  .le. - InfBnd) line(39: 54) = '      None      '
         if (b2  .ge.   InfBnd) line(55: 70) = '      None      '
         if (b1  .eq.   zero  ) line(39: 54) = '        .       '
         if (b2  .eq.   zero  ) line(55: 70) = '        .       '
         if (is  .eq.   0       .or.
     &       wlam.eq.   zero  ) line(71: 86) = '        .       '
         if (b1  .le. - InfBnd  .and.
     &       b2  .ge.   InfBnd) then
                                line(71: 86) = '                '
                                line(87:102) = '                '
         end if
         if (slk .eq.   zero  ) line(87:102) = '        .       '

         call snPRNT( 1, line, iw, leniw )
      end do

      return

 1000 format( 1x,  a15, 2x, 'State', 6x, 'Value',
     &        7x, 'Lower bound', 5x, 'Upper bound',
     &        3x, 'Lagr multiplier', 4x, '   Slack' )
 2000 format( 1x, a8, i6, 3x, a1, 1x, a2, 4g16.7, g16.4 )

      end ! subroutine s3prtN

