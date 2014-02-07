*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     file  sn30spec.f
*
*     s3optc   s3opti   s3optr   s3optl
*     s3file   s3key    s3opt    s3tie    s3undf
*     oplook   opnumb   opscan   optokn   opuppr
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3optc
     &   ( set, cwork, cvalue )

      implicit
     &     none
      logical
     &     set
      character
     &     cwork*8, cvalue*8

*     ==================================================================
*     s3optc  sets cwork to cvalue or vice versa depending on the value
*     of set.
*
*     17 May 1998: First version of s3optc.
*     17 May 1998: Current version.
*     ==================================================================

      if ( set ) then
         cwork  = cvalue
      else
         cvalue = cwork
      end if

      end ! subroutine s3optc

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3opti
     &   ( set, iwork, ivalue )

      implicit
     &     none
      logical
     &     set
      integer
     &     iwork, ivalue

*     ==================================================================
*     s3opti  sets iwork to ivalue or vice versa depending on the value
*     of set.
*
*     17 May 1998: First version of s3opti.
*     17 May 1998: Current version.
*     ==================================================================

      if ( set ) then
         iwork  = ivalue
      else
         ivalue = iwork
      end if

      end ! subroutine s3opti

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3optr
     &   ( set, rwork, rvalue )

      implicit
     &     none
      logical
     &     set
      double precision
     &     rwork, rvalue

*     ==================================================================
*     s3optr  sets rwork to rvalue or vice versa depending on the value
*     of set.
*
*     17 May 1998: First version of s3optr.
*     17 May 1998: Current version.
*     ==================================================================

      if ( set ) then
         rwork  = rvalue
      else
         rvalue = rwork
      end if

      end ! subroutine s3optr

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3optl
     &   ( set, iwork, ivalue, l )

      implicit
     &     none
      logical
     &     set
      integer
     &     iwork, ivalue, l

*     ==================================================================
*     If set=true, s3optz sets iwork = ivalue and ignores l.
*     Otherwise  , s3optz sets l     = 1 if iwork = ivalue, else l = 0.
*
*     01 Aug 2003: First version of s3optl.
*     01 Aug 2003: Current version of s3optl.
*     ==================================================================

      if ( set ) then
         iwork  = ivalue
      else
         l = 0
         if (ivalue .eq. iwork) l = 1
      end if

      end ! subroutine s3optl

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3file
     &   ( iExit, nCalls, iSpecs, opset,
     &     title, iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Errors, iExit, nCalls, iSpecs, iPrint, iSumm,
     &     lencw, leniw, lenrw, iw(leniw)
      character*(*)
     &     title
      character
     &     cw(lencw)*8
      double precision
     &     rw(lenrw)
      external
     &     opset

!     ==================================================================
!     s3file  reads the specs file from unit  iSpecs  and loads the
!     relevant options, using opset to process each line.
!
!     On exit, Errors says how many errors were encountered.
!
!     15 Nov 1991: First version based on Minos/Npsol routine s3file.
!     31 Jul 2003: snPRNT adopted.  iPrint, iSumm now used only by
!                  s3opt.  Beware -- they may get changed there.
!                  New input parameter "title" so we can remove s3fils.
!     27 Oct 2003: Errors counted separately from iExit.
!     06 Aug 2005: Fredrik Hellman found that str*72 is too short
!                  for format 2300.  Increased it to 80.
!     19 Dec 2005: iExit values clarified and reordered.
!     30 Apr 2006: snREAD adopted.
!     ==================================================================
      integer
     &     endfile, ivalue, lenbuf, nRead, nkey
      double precision
     &     rvalue
      character
     &     buffer*72, cvalue*8, dashes*30, key*16, str*80, token(1)*16
      data
     &     dashes /'=============================='/
!     ------------------------------------------------------------------

      iExit  = 0
      Errors = 0

!     Return if the unit number is out of range.

      if (iSpecs .lt. 0  .or.  iSpecs .gt. 99) then
         iExit = 131            ! iSpecs out of range
         return
      end if

!     ------------------------------------------------------------------
!     Look for  Begin, Endrun  or  Skip.
!     ------------------------------------------------------------------
      nRead  = 0
   50    call snREAD( iSpecs, buffer, 72, endfile )
!        read (iSpecs, '(a72)', end = 920) buffer
         if (endfile .gt. 0) go to 920
         nRead = nRead + 1
         call optokn( buffer, 1, nkey, token )
         key   = token(1)
         if (key .eq. 'ENDRUN') go to 940
         if (key .ne. 'BEGIN' ) then
            if (nRead .eq. 1  .and.  key .ne. 'SKIP') then
               Errors = Errors + 1
               write(str, 2000) iSpecs
               call snPRNT( 14, str, iw, leniw )
               call snPRNT(  4,
     &         ' XXX  The file should start with Begin, Skip or Endrun',
     &              iw, leniw )
               call snPRNT(  4,
     &         ' XXX  but the first record found was the following:',
     &              iw, leniw )
               call snPRNT( 14, ' ---->'//buffer(1:72), iw, leniw )
               call snPRNT( 14,
     &         ' XXX  Continuing to look for SPECS file...',
     &              iw, leniw )
            end if
            go to 50
         end if

*     ------------------------------------------------------------------
*     Begin found.
*     This is taken to be the first line of a SPECS file.
*     Print the title first if it's not blank.
*     ------------------------------------------------------------------
      call s1page( 1, iw, leniw )

      call s1trim( buffer, lenbuf )

      if (title .ne. ' ') then
         str = title
         call snPRNT( 13, ' ', iw, leniw )
         call snPRNT(  1, '         '//dashes, iw, leniw )
         call snPRNT(  1, '         '//str   , iw, leniw )
         call snPRNT(  1, '         '//dashes, iw, leniw )
         call snPRNT(  2,         ' '//dashes, iw, leniw )
         call snPRNT(  2,         ' '//str   , iw, leniw )
         call snPRNT(  2,         ' '//dashes, iw, leniw )
      end if

      call snPRNT( 11, ' SPECS file', iw, leniw )
      call snPRNT(  1, ' ----------', iw, leniw )
      call snPRNT( 11, '      '//buffer(1:lenbuf), iw, leniw )
      call snPRNT( 12,      ' '//buffer(1:lenbuf), iw, leniw )

*     ------------------------------------------------------------------
*     Read the rest of the file.
*     ------------------------------------------------------------------
*+    while (key .ne. 'END') do
  100 if    (key .ne. 'END') then
         call snREAD( iSpecs, buffer, 72, endfile )
         if (endfile .gt. 0) go to 930
         call opset
     &      ( .true., buffer, key, cvalue, ivalue, rvalue,
     &        iPrint, iSumm, Errors,
     &        cw, lencw, iw, leniw, rw, lenrw )
         go to 100
      end if
*+    end do

      return

  920 if (nCalls .le. 1) then
         write(str, 2200) iSpecs
         call snPRNT( 14, str, iw, leniw )
      else
         call snPRNT(  3, ' Endrun', iw, leniw )
      end if
      iExit = 132    ! End-of-file found while looking for BEGIN
      return

  930 write(str, 2300) iSpecs
      call snPRNT( 14, str, iw, leniw )
      iExit = 133    ! End-of-file found while reading specs
      return

  940 call snPRNT( 11, '      '//buffer(1:lenbuf), iw, leniw )
      call snPRNT( 12,      ' '//buffer(1:lenbuf), iw, leniw )
      iExit = 134    ! Endrun found for empty SPECS file
      return

 2000 format(' XXX  Error while looking for a SPECS file on unit', I6)
 2200 format(' XXX  End-of-file encountered while looking for',
     &       ' a BEGIN file on unit', I6)
 2300 format(' XXX  End-of-file encountered while processing',
     &       ' a SPECS file on unit', I6)

      end ! subroutine s3file

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3key ( key, loc )

      implicit
     &     none
      integer
     &     loc
      character
     &     key*16

*     ==================================================================
*     s3key  sets key to be the standard form for the first keyword
*     on each line of a SPECS file.
*
*     17 May 1998: First version.
*     13 Mar 2004: Hot start option added.
*     22 Jun 2004: System info option added.
*     22 Jun 2004: This version.
*     ==================================================================
      integer                 maxkey
      parameter         (     maxkey = 87)
      character          keys(maxkey)*16
      logical            sorted
      parameter         (sorted =   .true.)
*     ------------------------------------------------------------------
      data
     &   keys(  1) /'AIJ          '/,
     &   keys(  2) /'BACKUP       '/,
     &   keys(  3) /'BOUNDS       '/,
     &   keys(  4) /'CENTRAL      '/,
     &   keys(  5) /'CG           '/,
     &   keys(  6) /'CHECK        '/,
     &   keys(  7) /'COEFFICIENTS '/,
     &   keys(  8) /'COLD         '/,
     &   keys(  9) /'COLUMNS      '/,
     &   keys( 10) /'CRASH        '/,
     &   keys( 11) /'CYCLE        '/,
     &   keys( 12) /'DEBUG        '/,
     &   keys( 13) /'DEFAULTS     '/,
     &   keys( 14) /'DERIVATIVE   '/,
     &   keys( 15) /'DIFFERENCE   '/,
     &   keys( 16) /'DUMP         '/,
     &   keys( 17) /'ELASTIC      '/,
     &   keys( 18) /'ELEMENTS     '/,
     &   keys( 19) /'ERROR        '/,
     &   keys( 20) /'EXPAND       '/

      data
     &   keys( 21) /'FACTORIZATION'/,
     &   keys( 22) /'FEASIBILITY  '/,
     &   keys( 23) /'FEASIBLE     '/,
     &   keys( 24) /'FUNCTION     '/,
     &   keys( 25) /'HESSIAN      '/,
     &   keys( 26) /'HOT          '/,
     &   keys( 27) /'INFEASIBLE   '/,
     &   keys( 28) /'INFINITE     '/,
     &   keys( 29) /'INSERT       '/,
     &   keys( 30) /'ITERATIONS   '/,
     &   keys( 31) /'IW           '/,
     &   keys( 32) /'JACOBIAN     '/,
     &   keys( 33) /'LINESEARCH   '/,
     &   keys( 34) /'LIST         '/,
     &   keys( 35) /'LOAD         '/,
     &   keys( 36) /'LOG          '/,
     &   keys( 37) /'LOWER        '/,
     &   keys( 38) /'LP           '/,
     &   keys( 39) /'LU           '/,
     &   keys( 40) /'MAJOR        '/

      data
     &   keys( 41) /'MAXIMIZE     '/,
     &   keys( 42) /'MINIMIZE     '/,
     &   keys( 43) /'MINOR        '/,
     &   keys( 44) /'MPS          '/,
     &   keys( 45) /'NEW          '/,
     &   keys( 46) /'NO           '/,
     &   keys( 47) /'NON          '/,
     &   keys( 48) /'NONDERIVATIVE'/,
     &   keys( 49) /'NONLINEAR    '/,
     &   keys( 50) /'OBJECTIVE    '/,
     &   keys( 51) /'OLD          '/,
     &   keys( 52) /'OPTIMALITY   '/,
     &   keys( 53) /'PARTIAL      '/,
     &   keys( 54) /'PENALTY      '/,
     &   keys( 55) /'PIVOT        '/,
     &   keys( 56) /'PRINT        '/,
     &   keys( 57) /'PROBLEM      '/,
     &   keys( 58) /'PROXIMAL     '/,
     &   keys( 59) /'PUNCH        '/,
     &   keys( 60) /'QP           '/

      data
     &   keys( 61) /'QPSOLVER     '/,
     &   keys( 62) /'REDUCED      '/,
     &   keys( 63) /'RANGES       '/,
     &   keys( 64) /'REPORT       '/,
     &   keys( 65) /'RHS          '/,
     &   keys( 66) /'ROWS         '/,
     &   keys( 67) /'RW           '/,
     &   keys( 68) /'SAVE         '/,
     &   keys( 69) /'SCALE        '/,
     &   keys( 70) /'SOLUTION     '/,
     &   keys( 71) /'START        '/,
     &   keys( 72) /'STOP         '/,
     &   keys( 73) /'SUBSPACE     '/,
     &   keys( 74) /'SUMMARY      '/,
     &   keys( 75) /'SUPERBASICS  '/,
     &   keys( 76) /'SUPPRESS     '/,
     &   keys( 77) /'SYSTEM       '/,
     &   keys( 78) /'TIMING       '/,
     &   keys( 79) /'TOTAL        '/,
     &   keys( 80) /'UNBOUNDED    '/

      data
     &   keys( 81) /'UPPER        '/,
     &   keys( 82) /'USER         '/,
     &   keys( 83) /'VERIFY       '/,
     &   keys( 84) /'VIOLATION    '/,
     &   keys( 85) /'WARM         '/,
     &   keys( 86) /'WORKING      '/,
     &   keys( 87) /'WORKSPACE    '/
*     ------------------------------------------------------------------
      call oplook( maxkey, keys, sorted, key, loc )

      end ! subroutine s3key

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3opt
     &   ( s, buffer, key, c, i, r, lPrnt, lSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     s
      integer
     &     lPrnt, lSumm, Errors, i, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     r, rw(lenrw)
      character
     &     c*8, cw(lencw)*8
      character*(*)
     &     buffer, key

*     ==================================================================
*     s3opt  decodes the option contained in  buffer  in order to
*     set or get a parameter value in the relevant array iw or rw.
*
*     The buffer is output to file iPrint, minus trailing blanks.
*     Error messages are output to files iPrint and iSumm.
*     buffer is echoed to iPrint but normally not to iSumm.
*     It is echoed to iSumm before any error msg.
*
*     On entry,
*     buffer contains the option string to be processed.
*     s      is true if an option is to be extracted from buffer.
*            Otherwise, c, i and r are to be assigned the value of the
*            option defined in the option string.
*     lPrnt  is iPrint as given to s3file.
*     lSumm  is iSumm  as given to s3file.
*     Errors is the number of errors so far.
*
*     On exit,
*     key    is the first keyword contained in buffer.
*     If s is true, c, i and r may be ignored.  (They are usually
*            option values that have been saved in cw, iw, rw.)
*     If s is false,
*     c      is the OBJECTIVE, RHS, RANGE or BOUND name if key is
*            one of those words.
*     r      is the first numerical value found in buffer (or zero).
*     i      is int(r) if  abs(r) < maxint.
*     Errors is the number of errors so far.
*
*
*     s3opt  uses opnumb and the subprograms
*                 lookup, scannr, tokens, upcase
*     (now called oplook, opscan, optokn, opuppr)
*     supplied by Sterling Software, Palo Alto, California.
*
*     15 Nov 1991: First version based on s3key/opkey.
*     10 Dec 2002: Recognize LU Diagonal and LU Rook Pivoting
*     02 Jul 2003: Options for QN QP solver added.
*     01 Aug 2003: snPRNT adopted.  If lPrnt, lSumm are positive,
*                  snPRNT outputs to global iPrint, iSumm, iStdo.
*                  Sometimes lPrnt, lSumm are zero (for "get" routines).
*                  snPRNT( 4, ... ) then outputs error msgs
*     21 Dec 2003: Hot start added.
*     22 Jun 2004: lvlSys added
*     04 Dec 2004: kosher maxint added
*     22 Apr 2007: Allowed pure comments to start in any column
*     ==================================================================
      external
     &     opnumb, s1intmx
      logical
     &     opnumb, more, number
      character
     &     key2*16, key3*16, value*16, str*132, str1*132
      integer
     &     cgItmx, DerOpt, i0, i1, i2, i3,
     &     lenb, lenbuf, loc1, loc2, m1, nToken,
     &     iBack, iCrash, iDump, iInsrt, iLoadB, iMPS, iNewB,
     &     iOldB, iPnch, iPrint, iReprt, iSoln, iSumm, itnlim, jverf1,
     &     jverf2, jverf3, jverf4, kchk, kDegen, kfac, klog, kReset,
     &     ksav, kSumm, lDenJ, lEmode, lprDbg, lprPrm, lprSch, lprScl,
     &     lprSol, lvlDer, lvlHes, lvlInf, lvlPiv, lvlPre, lvlPPm,
     &     lvlSch, lvlScl, lvlSrt, lvlSys, lvlTim, lvlVer, maxbuf,
     &     maxcu, maxcw, maxint, maxiu, maxiw, maxm, maxn, maxne,
     &     maxR, maxru, maxrw, maxS, mBnd, mEr, mFlush, minmax, MjrPrt,
     &     mLst, mMajor, mMinor, MnrPrt, mObj, mQNmod, mRhs, mRng,
     &     nnCon, nnJac, nnObj, nnL, nParPr, nProb, ObjRow, QPslvr,
     &     s1intmx
      integer
     &     etarg, tolCG, tolFP, tolQP, tolNLP, tolx, tolCon, tolpiv,
     &     tCrash, tolswp, tolfac, tolupd, InfBnd, bigFx, bigdx, epsrf,
     &     fdint1, fdint2, xdlim, vilim, wolfeG, wtInf0, mNewSB, xPen0,
     &     scltol, Aijtol, bStrc1, bStrc2, Utol1, Utol2, Dens2, wtMax
*     ------------------------------------------------------------------
      integer                  maxtok
      parameter         (      maxtok = 10)
      character          token(maxtok)*16

      double precision   zero
      parameter         (zero   =   0.0d+0)

      parameter         (tolFP     =  51) ! Minor Phase 1 Opt tol
      parameter         (tolQP     =  52) ! Minor Phase 2 Opt tol
      parameter         (tolNLP    =  53) ! Major Optimality tolerance
      parameter         (tolCG     =  54) ! cg tolerance
      parameter         (tolx      =  56) ! Minor feasibility tolerance
      parameter         (tolCon    =  57) ! Major feasibility tolerance
      parameter         (tolpiv    =  60) ! excludes small pivot elems
      parameter         (tCrash    =  62) ! crash tolerance
      parameter         (tolswp    =  65) ! LU swap tolerance
      parameter         (tolfac    =  66) ! LU factor tolerance
      parameter         (tolupd    =  67) ! LU update tolerance
      parameter         (InfBnd    =  70) ! definition of plus infinity
      parameter         (bigFx     =  71) ! unbounded objective
      parameter         (bigdx     =  72) ! unbounded step
      parameter         (epsrf     =  73) ! relative function precision
      parameter         (fdint1    =  76) ! forward difference interval
      parameter         (fdint2    =  77) ! central difference interval
      parameter         (xdlim     =  80) ! Step limit
      parameter         (vilim     =  81) ! violation limit
      parameter         (etarg     =  83) ! Quasi-Newton QP rg tolerance
      parameter         (wolfeG    =  84) ! line search tolerance
      parameter         (wtInf0    =  88) ! infeasibility weight
      parameter         (xPen0     =  89) ! initial penalty parameter
      parameter         (wtMax     =  90) ! Elastic weightmax
      parameter         (scltol    =  92) ! scale tolerance
      parameter         (Aijtol    =  95) ! zero Aij tolerance
      parameter         (bStrc1    =  96) ! default lower bound on x
      parameter         (bStrc2    =  97) ! default upper bound on x
      parameter         (Utol1     = 154) ! abs tol for small diag of U
      parameter         (Utol2     = 155) ! rel tol for small diag of U
      parameter         (Dens2     = 158) ! switch to dense LU
      parameter         (maxru     =   2) ! Start of SNOPT part of rw
      parameter         (maxrw     =   3) ! End   of SNOPT part of rw
      parameter         (maxiu     =   4) ! Start of SNOPT part of iw
      parameter         (maxiw     =   5) ! End   of SNOPT part of iw
      parameter         (maxcu     =   6) ! Start of SNOPT part of cw
      parameter         (maxcw     =   7) ! End   of SNOPT part of cw
      parameter         (iPrint    =  12) ! Print   file
      parameter         (iSumm     =  13) ! Summary file
      parameter         (nnJac     =  21) ! # nonlinear Jacobian variables
      parameter         (nnObj     =  22) ! # variables in gObj
      parameter         (nnCon     =  23) ! nonlinear constraints
      parameter         (nnL       =  24) !   max( nnObj, nnJac )
      parameter         (maxR      =  52) ! max columns of R
      parameter         (maxS      =  53) ! max # of superbasics
      parameter         (mQNmod    =  54) ! (ge 0) max # of BFGS updates
      parameter         (QPslvr    =  55) ! 0(1) => QP(QN) QP solver
      parameter         (lEmode    =  56) ! >0    => use elastic mode
      parameter         (kchk      =  58) ! check (row) frequency
      parameter         (kfac      =  59) ! factorization frequency
      parameter         (ksav      =  60) ! save basis map
      parameter         (klog      =  61) ! log/print frequency
      parameter         (kSumm     =  62) ! Summary print frequency
      parameter         (kDegen    =  63) ! max. expansions of featol
      parameter         (kReset    =  64) ! Hessian frequency
      parameter         (mFlush    =  66) ! Hessian flush
      parameter         (lvlSrt    =  69) ! = 0(1) => cold(warm) start
      parameter         (lvlDer    =  70) ! derivative level
      parameter         (lvlSys    =  71) ! > 0   => print system info
      parameter         (lvlHes    =  72) ! 0,1,2 => LM, FM, Exact H
      parameter         (lvlInf    =  73) ! Elastic option
      parameter         (lvlScl    =  75) ! scale option
      parameter         (lvlSch    =  76) ! >0 => deriv. line search
      parameter         (lvlPre    =  77) ! >0 => QN preconditioned CG
      parameter         (lvlVer    =  78) ! Verify level
      parameter         (lvlPPm    =  79) ! Proximal Point method for x0
      parameter         (lvlPiv    =  80) ! 0(1 2 3) LU pivoting
      parameter         (lprPrm    =  81) ! > 0  =>  parms are printed
      parameter         (lprSch    =  82) ! line search debug start itn
      parameter         (lprScl    =  83) ! > 0  => print the scales
      parameter         (lprSol    =  84) ! > 0  =>  print the solution
      parameter         (lprDbg    =  85) ! > 0  => private debug print
      parameter         (minmax    =  87) ! 1, -1  => MIN, FP, MAX
      parameter         (iCrash    =  88) ! Crash option
      parameter         (itnlim    =  89) ! limit on total iterations
      parameter         (mMajor    =  90) ! limit on major iterations
      parameter         (mMinor    =  91) ! limit on minor iterations
      parameter         (MjrPrt    =  92) ! Major print level
      parameter         (MnrPrt    =  93) ! Minor print level
      parameter         (nParPr    =  94) ! # partial pricing sections
      parameter         (mNewSB    =  95) ! maximum # of new SB
      parameter         (jverf1    =  98) ! start derivative checking
      parameter         (jverf2    =  99) ! stop  derivative checking
      parameter         (jverf3    = 100) ! start derivative checking
      parameter         (jverf4    = 101) ! stop  derivative checking
      parameter         (ObjRow    = 103) ! Objective row
      parameter         (DerOpt    = 104) ! 0, 1, 2 => derivative option
      parameter         (lDenJ     = 105) ! 1(2) => dense(sparse) deriv.
      parameter         (mEr       = 106) ! maximum # errors in MPS data
      parameter         (mLst      = 107) ! maximum # lines  of MPS data
      parameter         (nProb     = 108) ! problem number
      parameter         (cgItmx    = 111) ! CG iteration limit
      parameter         (iBack     = 120) ! backup file
      parameter         (iDump     = 121) ! dump file
      parameter         (iLoadB    = 122) ! load file
      parameter         (iMPS      = 123) ! MPS file
      parameter         (iNewB     = 124) ! new basis file
      parameter         (iInsrt    = 125) ! insert file
      parameter         (iOldB     = 126) ! old basis file
      parameter         (iPnch     = 127) ! punch file
      parameter         (iReprt    = 130) ! Report file
      parameter         (iSoln     = 131) ! Solution file
      parameter         (maxm      = 133) ! Row    estimate
      parameter         (maxn      = 134) ! Column estimate
      parameter         (maxne     = 135) ! Estimated element count
      parameter         (lvlTim    = 182) ! Timing level

      parameter         (mObj      =  52) ! Objective name
      parameter         (mRhs      =  53) ! Right-hand side name
      parameter         (mRng      =  54) ! Range name
      parameter         (mBnd      =  55) ! Bnd section name
*     ------------------------------------------------------------------
      maxint = s1intmx( )
*     ------------------------------------------------------------------
*     Trim trailing blanks and echo to the Print file.
*     ------------------------------------------------------------------
      call s1trim( buffer, lenbuf )
      maxbuf = min( 120,lenbuf)

      if (lPrnt .gt. 0) then
         write(str, '(6x,a)') buffer(1:maxbuf)
         call snPRNT( 1, str, iw, leniw )
      end if

*     Set lenb = length of buffer without trailing comments.
*     Eliminate comments and empty lines.
*     A '*' appearing anywhere in buffer terminates the string.

      i  = index( buffer(1:lenbuf), '*' )
      if (i .eq. 0) then
         lenb = lenbuf
      else
         lenb = i - 1
      end if
      if (lenb .le. 0) then
         key = '*'
         go to 900
      end if

*     ------------------------------------------------------------------
*     Extract up to maxtok tokens from the record.
*     ntoken returns how many were actually found.
*     key, key2, are the first tokens if any, otherwise blank.
*     For some values of key (bounds, objective, ranges, rhs)
*     we have to save key2 before s3tie (and oplook) alter it.
*     For example, if the data is     objective = obj
*     oplook will change obj to objective.
*     ------------------------------------------------------------------
      call optokn( buffer(1:lenbuf), maxtok, ntoken, token )
      key    = token(1)
      key2   = token(2)
      key3   = token(3)
      c      = key2(1:8)

*     Certain keywords require no action.

      if (key .eq. '   ') go to 900 ! blank line
      if (key .eq. '*  ') go to 900 ! comment starting in column no. > 1
      if (key .eq. 'END') go to 900

*     Convert the keywords to their most fundamental form
*     (upper case, no abbreviations).
*     loci   says where the keywords are in the dictionaries.
*     loci = 0 signals that the keyword wasn't there.

      call s3key ( key , loc1 )
      call s3tie ( key2, loc2 )

*     Most keywords will have an associated integer or real value,
*     so look for it no matter what the keyword.

      c      = key2(1:8)
      i      = 1
      number = .false.

*+    while (i .lt. ntoken  .and.  .not. number) loop
   50 if    (i .lt. ntoken  .and.  .not. number) then
         i      = i + 1
         value  = token(i)
         number = opnumb( value )
         go to 50
      end if
*+    end while

      i = 0
      r = zero
      if ( number ) then
         read  (value, '(bn, e16.0)') r
         i = maxint
         if (abs(r) .lt. maxint) i = int(r)
      end if

*     ------------------------------------------------------------------
*     Decide what to do about each keyword.
*     The second keyword (if any) might be needed to break ties.
*     Some seemingly redundant testing of more is used
*     to avoid compiler limits on the number of consecutive else ifs.
*     ------------------------------------------------------------------
      m1     = -1
      i0     =  0
      i1     =  1
      i2     =  2
      i3     =  3
      more   = .true.

      if (more) then
         more   = .false.
         if      (key .eq. 'BACKUP      ') then
            call s3opti(s, iw(iBack ), i)

         else if (key .eq. 'CENTRAL     ') then
            call s3optr(s, rw(fdint2), r)

         else if (key .eq. 'CG          ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'TOLERANCE   ') call s3optr(s, rw(tolCG ), r)
              if (key2.eq. 'PRECONDITIONING')
     &                                     call s3opti(s, iw(lvlPre), i)
              if (key2.eq. 'ITERATIONS  ') call s3opti(s, iw(cgItmx), i)

         else if (key .eq. 'CHECK       ') then
            call s3opti(s, iw(kchk  ), i)

         else if (key .eq. 'COLD        ') then
            call s3opti(s, iw(lvlSrt),i0)

         else if (key .eq. 'CRASH       ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'OPTION      ') call s3opti(s, iw(iCrash), i)
              if (key2.eq. 'TOLERANCE   ') call s3optr(s, rw(tCrash), r)

         else if (key .eq. 'DEBUG       ') then
            call s3opti(s, iw(lprDbg), i)
         else if (key .eq. 'DEFAULTS    ') then
            call s3undf( cw, lencw, iw, leniw, rw, lenrw )

         else if (key .eq. 'DERIVATIVE  ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'LEVEL       ') call s3opti(s, iw(lvlDer), i)
              if (key2.eq. 'LINESEARCH  ')call s3optl(s,iw(lvlSch),i1,i)
              if (key2.eq. 'OPTION      ') call s3opti(s, iw(DerOpt), i)

         else if (key .eq. 'DIFFERENCE  ') then
            call s3optr(s, rw(fdint1), r)
         else if (key .eq. 'DUMP        ') then
            call s3opti(s, iw(iDump ), i)

         else if (key .eq. 'ELASTIC     ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'OBJECTIVE   ') call s3opti(s, iw(lvlInf), i)
              if (key2.eq. 'MODE        ') call s3opti(s, iw(lEmode), i)
              if (key2.eq. 'WEIGHT      ') call s3optr(s, rw(wtInf0), r)
              if (key2.eq. 'WEIGHTMAX   ') call s3optr(s, rw(wtMax ), r)

         else if (key .eq. 'EXPAND      ') then
            call s3opti(s, iw(kDegen), i)
         else if (key .eq. 'FACTORIZATION') then
            call s3opti(s, iw(kfac  ), i)

         else if (key .eq. 'FEASIBLE    ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'POINT       ')call s3optl(s,iw(minmax),i0,i)
              if (key2.eq. 'EXIT        ') go to 890

         else if (key .eq. 'FEASIBILITY ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'TOLERANCE   ') call s3optr(s, rw(tolx  ), r)

         else if (key .eq. 'FUNCTION    ') then
            call s3optr(s, rw(epsrf ), r)

         else if (key .eq. 'HESSIAN     ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'COLUMNS     ') call s3opti(s, iw(nnL   ), i)
              if (key2.eq. 'DIMENSION   ') call s3opti(s, iw(maxR  ), i)
              if (key2.eq. 'FREQUENCY   ') call s3opti(s, iw(kReset), i)
              if (key2.eq. 'FLUSH       ') call s3opti(s, iw(mFlush), i)
              if (key2.eq. 'UPDATES     ') call s3opti(s, iw(mQNmod), i)
              if (key2.eq. 'LIMITED     ')call s3optl(s,iw(lvlHes),i0,i)
              if (key2.eq. 'FULL        ')call s3optl(s,iw(lvlHes),i1,i)
              if (key2.eq. 'PRECONDITIONING')
     &                                     call s3opti(s, iw(lvlPre), i)

         else if (key .eq. 'HOT         ') then
            call s3opti(s, iw(lvlSrt),i3)
         else if (key .eq. 'INFINITE    ') then
            call s3optr(s, rw(InfBnd), r)
         else if (key .eq. 'INSERT      ') then
            call s3opti(s, iw(iInsrt), i)
         else if (key .eq. 'ITERATIONS  ') then
            call s3opti(s, iw(itnlim), i)
         else
            more   = .true.
         end if
      end if

      if (more) then
         more   = .false.
         if      (key .eq. 'IW          ') then
*           Allow things like  iw 21 = 100  to set iw(21) = 100
            key2   = token(3)
            if (i .ge. 1  .and. i .le. 500) then
               read (key2, '(bn, i16)') iw(i)
            else
               go to 880
            end if

         else if (key .eq. 'LINESEARCH  ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'TOLERANCE   ') call s3optr(s, rw(wolfeG), r)
              if (key2.eq. 'DEBUG       ') call s3opti(s, iw(lprSch), i)

          else if (key .eq. 'LOAD        ') then
            call s3opti(s, iw(iLoadB), i)
          else if (key .eq. 'LOG         ') then
            call s3opti(s, iw(klog)  , i)

         else if (key .eq. 'LP          ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'FEASIBILITY ') call s3optr(s, rw(tolx  ), r)
              if (key2.eq. 'OPTIMALITY  ') call s3optr(s, rw(tolQP ), r)

         else if (key .eq. 'LU          ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'PARTIAL     ')call s3optl(s,iw(lvlPiv),i0,i)
              if (key2.eq. 'COMPLETE    ')call s3optl(s,iw(lvlPiv),i2,i)
              if (key2.eq. 'DIAGONAL    ')call s3optl(s,iw(lvlPiv),i3,i)
              if (key2.eq. 'FACTORIZATION')call s3optr(s, rw(tolFac), r)
              if (key2.eq. 'ROOK        ')call s3optl(s,iw(lvlPiv),i1,i)
              if (key2.eq. 'UPDATES     ') call s3optr(s, rw(tolUpd), r)
              if (key2.eq. 'DENSITY     ') call s3optr(s, rw(Dens2 ), r)
              if (key2.eq. 'SINGULARITY ') then
                 call s3optr(s, rw(Utol1), r)
                 call s3optr(s, rw(Utol2), r)
              end if
              if (key2.eq. 'SWAP        ') call s3optr(s, rw(tolswp), r)
!              if (key2.eq. 'DEFAULTS    ') then
!                 if (loc3.eq.  0           ) go to 820
!                 if (key3.eq.'TPP         ') call s3optr(s,rw(tolDpp),r)
!                 if (key3.eq.'TCP         ') call s3optr(s,rw(tolDcp),r)
!                 if (key3.eq.'UPDATES     ') call s3optr(s,rw(tolDup),r)
!              end if
         else
            more   = .true.
         end if
      end if

      if (more) then
         more   = .false.
         if      (key .eq. 'MAJOR       ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'FEASIBILITY ') call s3optr(s, rw(tolCon), r)
              if (key2.eq. 'ITERATIONS  ') call s3opti(s, iw(mMajor), i)
              if (key2.eq. 'OPTIMALITY  ') call s3optr(s, rw(tolNLP), r)
              if (key2.eq. 'PRINT       ') call s3opti(s, iw(MjrPrt), i)
              if (key2.eq. 'STEP        ') call s3optr(s, rw(xdlim ), r)

         else if (key .eq. 'MAXIMIZE    ') then
            call s3opti(s, iw(minmax), m1)
         else if (key .eq. 'MINIMIZE    ') then
            call s3opti(s, iw(minmax), i1)

         else if (key .eq. 'MINOR       ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'ITERATIONS  ') call s3opti(s, iw(mMinor), i)
              if (key2.eq. 'FEASIBILITY ') call s3optr(s, rw(tolx  ), r)
              if (key2.eq. 'OPTIMALITY  ') call s3optr(s, rw(tolQP ), r)
              if (key2.eq. 'PHASE1      ') call s3optr(s, rw(tolFP ), r)
              if (key2.eq. 'PHASE2      ') call s3optr(s, rw(tolQP ), r)
              if (key2.eq. 'PRINT       ') call s3opti(s, iw(MnrPrt), i)
              if (key2.eq. 'SUPERBASICS ') call s3opti(s, iw(mNewSB), i)

         else if (key .eq. 'NEW         ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'BASIS       ') call s3opti(s, iw(iNewB ), i)
              if (key2.eq. 'SUPERBASICS ') call s3opti(s, iw(mNewSB), i)

         else if (key .eq. 'NO          '  .or.
     &            key .eq. 'NONDERIVATIVE' .or.
     &            key .eq. 'NON         ') then
            call s3opti(s, iw(lvlSch),i0)

         else if (key .eq. 'OBJECTIVE   ') then
              if (key2.eq. 'ROW         ') then
                                           call s3opti(s, iw(ObjRow), i)
              else
                                           call s3optc(s, cw(mObj  ), c)
              end if
         else if (key .eq. 'OLD         ') then
            call s3opti(s, iw(iOldB ), i)
         else if (key .eq. 'OPTIMALITY  ') then
            call s3optr(s, rw(tolNLP), r)
         else if (key .eq. 'PARTIAL     ') then
            call s3opti(s, iw(nParPr), i)
         else if (key .eq. 'PENALTY     ') then
            call s3optr(s, rw(xPen0 ), r)
         else if (key .eq. 'PIVOT       ') then
            call s3optr(s, rw(tolpiv), r)

         else if (key .eq. 'PRINT       ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'FILE        ') call s3opti(s, iw(iPrint), i)
              if (key2.eq. 'FREQUENCY   ') call s3opti(s, iw(klog  ), i)
              if (key2.eq. 'LEVEL       ') call s3opti(s, iw(MnrPrt), i)

         else if (key .eq. 'PROXIMAL    ') then
            call s3opti(s, iw(lvlPPm), i)

         else if (key .eq. 'PUNCH       ') then
            call s3opti(s, iw(iPnch ), i)
         else
            more   = .true.
         end if
      end if

      if (more) then
         more   = .false.
         if      (key .eq. 'QP          ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'COLUMNS     ') call s3opti(s, iw(nnL  ) , i)
              if (key2.eq. 'FEASIBILITY ') call s3optr(s, rw(tolx ) , r)
              if (key2.eq. 'OPTIMALITY  ') call s3optr(s, rw(tolQP) , r)

         else if (key .eq. 'QPSOLVER    ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'CHOLESKY    ')call s3optl(s,iw(QPslvr),i0,i)
              if (key2.eq. 'CG          ')call s3optl(s,iw(QPslvr),i1,i)
              if (key2.eq. 'QN          ')call s3optl(s,iw(QPslvr),i2,i)

         else if (key .eq. 'REDUCED     ') then
            call s3opti(s, iw(maxR  ), i)

         else if (key .eq. 'REPORT      ') then
            call s3opti(s, iw(iReprt), i)

         else if (key .eq. 'ROWS        ') then
*             gams should recognize row tolerance
*             but not just          rows
*             This is a relic from MINOS
              if (key2.eq. 'TOLERANCE   ') then
                 call s3optr(s, rw(tolCon), r)
              else
                 call s3opti(s, iw(maxm  ), i)
              end if

         else if (key .eq. 'RW          ') then
*           allow things like rw 21 = 2  to set rw(21) = 2.0
            key2   = token(3)
            if (i .ge. 1  .and. i .le. 500) then
               read (key2, '(bn, e16.0)') rw(i)
            else
               go to 880
            end if

         else if (key .eq. 'SAVE        ') then
            call s3opti(s, iw(ksav  ), i)

         else if (key .eq. 'SCALE       ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'OPTION      ') call s3opti(s, iw(lvlScl), i)
              if (key2.eq. 'TOLERANCE   ') call s3optr(s, rw(scltol), r)
              if (key2.eq. 'PRINT       ')call s3optl(s,iw(lprScl),i1,i)
         else
            more   = .true.
         end if
      end if

      if (more) then
         more   = .false.
         if      (key .eq. 'SOLUTION    ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'FILE        ') call s3opti(s, iw(iSoln ), i)
              if (key2.eq. 'YES         ')call s3optl(s,iw(lprSol),i2,i)
              if (key2.eq. 'NO          ')call s3optl(s,iw(lprSol),i0,i)

         else if (key .eq. 'START       ') then
              if (key2.eq. 'OBJECTIVE   ') call s3opti(s, iw(jverf1), i)
              if (key2.eq. 'CONSTRAINTS ') call s3opti(s, iw(jverf3), i)
              if (loc2.eq.  0            ) go to 820
         else if (key .eq. 'STOP        ') then
              if (key2.eq. 'OBJECTIVE   ') call s3opti(s, iw(jverf2), i)
              if (key2.eq. 'CONSTRAINTS ') call s3opti(s, iw(jverf4), i)
              if (loc2.eq.  0            ) go to 820
         else
            more   = .true.
         end if
      end if

      if (more) then
         more   = .false.
         if      (key .eq. 'SUBSPACE    ') then
            call s3optr(s, rw(etarg ), r)
         else if (key .eq. 'SUPERBASICS ') then
            call s3opti(s, iw(maxS  ), i)
         else if (key .eq. 'SUMMARY     ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'FILE        ') call s3opti(s, iw(iSumm ), i)
              if (key2.eq. 'FREQUENCY   ') call s3opti(s, iw(ksumm ), i)

         else if (key .eq. 'SUPPRESS    ') then
            call s3opti(s, iw(lprPrm), i)
         else if (key .eq. 'TIMING      ') then
            call s3opti(s, iw(lvlTim), i)

         else if (key .eq. 'SYSTEM      ') then
              if (loc2.eq.  0            ) go to 820
              if (key3.eq. 'YES         ')call s3optl(s,iw(lvlSys),i1,i)
              if (key3.eq. 'NO          ')call s3optl(s,iw(lvlSys),i0,i)

         else if (key .eq. 'TOTAL       ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'INTEGER     ') call s3opti(s, iw(maxiw ), i)
              if (key2.eq. 'REAL        ') call s3opti(s, iw(maxrw ), i)
              if (key2.eq. 'CHARACTER   ') call s3opti(s, iw(maxcw ), i)

         else if (key .eq. 'UNBOUNDED   ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'OBJECTIVE   ') call s3optr(s, rw(bigFx ), r)
              if (key2.eq. 'STEP        ') call s3optr(s, rw(bigdx ), r)

         else if (key .eq. 'USER        ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'INTEGER     ') call s3opti(s, iw(maxiu ), i)
              if (key2.eq. 'REAL        ') call s3opti(s, iw(maxru ), i)
              if (key2.eq. 'CHARACTER   ') call s3opti(s, iw(maxcu ), i)

         else if (key .eq. 'VERIFY      ') then
              if (key2.eq. '            ') then
                 loc2   = 1
                 i      = 3
              end if
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'OBJECTIVE   ') i = 1
              if (key2.eq. 'CONSTRAINTS ') i = 2
              if (key2.eq. 'GRADIENTS   ') i = 3
              if (key2.eq. 'YES         ') i = 3
              if (key2.eq. 'NO          ') i = 0
              if (key2.eq. 'LEVEL       ') i = i
              call s3opti(s, iw(lvlVer), i)

         else if (key .eq. 'VIOLATION   ') then
            call s3optr(s, rw(vilim ), r)
         else if (key .eq. 'WARM        ') then
            call s3opti(s, iw(lvlSrt),i2)
         else if (key .eq. 'WORKSPACE   ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. '(USER)      ') call s3opti(s, iw(maxru), i)
              if (key2.eq. '(TOTAL)     ') call s3opti(s, iw(maxrw), i)
         else
            more   = .true.
         end if
      end if

      if (.not. more) go to 900

*     ------------------------------------------------------------------
*     Keywords for MPS files.
*     ------------------------------------------------------------------

      if (more) then
         more   = .false.
         if      (key .eq. 'AIJ         ') then
            call s3optr(s, rw(Aijtol), r)
         else if (key .eq. 'BOUNDS      ') then
            call s3optc(s, cw(mBnd  ), c)
         else if (key .eq. 'COEFFICIENTS') then
            call s3opti(s, iw(maxne ), i)
         else if (key .eq. 'COLUMNS     ') then
            call s3opti(s, iw(maxn  ), i)
         else if (key .eq. 'ELEMENTS    ') then
            call s3opti(s, iw(maxne ), i)
         else if (key .eq. 'ERROR       ') then
            call s3opti(s, iw(mEr   ), i)
         else if (key .eq. 'INFINITE    ') then
            call s3optr(s, rw(InfBnd), r)

         else if (key .eq. 'JACOBIAN    ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'DENSE       ')call s3optl(s,iw(lDenJ ),i1,i)
              if (key2.eq. 'SPARSE      ')call s3optl(s,iw(lDenJ ),i2,i)

         else if (key .eq. 'LIST        ') then
            call s3opti(s, iw(mLst  ), i)
         else if (key .eq. 'LOWER       ') then
            call s3optr(s, rw(bStrc1), r)
         else if (key .eq. 'MPS         ') then
            call s3opti(s, iw(iMPS  ), i)

         else if (key .eq. 'NONLINEAR   ') then
              if (loc2.eq.  0            ) go to 820
              if (key2.eq. 'CONSTRAINTS ') call s3opti(s, iw(nnCon ), i)
              if (key2.eq. 'OBJECTIVE   ') call s3opti(s, iw(nnObj ), i)
              if (key2.eq. 'JACOBIAN    ') call s3opti(s, iw(nnJac ), i)
              if (key2.eq. 'VARIABLES   ') then
                 call s3opti(s, iw(nnObj), i)
                 call s3opti(s, iw(nnJac), i)
              end if

         else if (key .eq. 'OBJECTIVE   ') then
              call s3optc(s, cw(mObj  ), c)
         else if (key .eq. 'PROBLEM     ') then
            call s3opti(s, iw(nProb ), i)
         else if (key .eq. 'RANGES      ') then
            call s3optc(s, cw(mRng  ), c)
         else if (key .eq. 'RHS         ') then
            call s3optc(s, cw(mRhs  ), c)
         else if (key .eq. 'UPPER       ') then
            call s3optr(s, rw(Bstrc2), r)
         else
            more   = .true.
         end if
      end if

      if (.not. more) go to 900

*     ------------------------------------------------------------------
*     Error messages.
*     This is the only way we can think of to concatenate strings
*     when one of them is of indeterminate length.
*     ------------------------------------------------------------------
      write(str, '(2a)') ' XXX  Keyword not recognized:         ', key
      go to 895

  820 write(str, '(2a)') ' XXX  Second keyword not recognized:  ', key2
      go to 895

  840 write(str, '(2a)') ' XXX  Fourth keyword not recognized:  ', key2
      go to 895

  880 write(str,'(a,i8)')' XXX  parm subscript out of range:    ', i
      go to 895

  890 str    = ' XXX  Obsolete option'
      go to 895

!     The buffer should have been output already to the Print file.
!     First output it to the Summary file.
!     Then print the error message.

  895 Errors = Errors + 1
      if (lSumm .gt. 0) then
         write(str1, '(1x,a)') buffer(1:maxbuf)
         call snPRNT( 2, str1, iw, leniw )
      end if
      call snPRNT( 4, str, iw, leniw )

  900 return

      end ! subroutine s3opt

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3tie ( tie, loc )

      implicit
     &     none
      integer
     &     loc
      character
     &     tie*16

*     ==================================================================
*     s3tie  sets key to be the standard form for the second keyword
*     on each line of a SPECS file.
*
*     21 May 1998: First version of s3tie.
*     10 Dec 2002: New ties: Diagonal and Rook
*     17 Jan 2003: New tie: Solver
*     01 Aug 2003: New ties for QPsolver: Cholesky, CG, QN
*     22 Jun 2004: New tie for System.
*     22 Jun 2004: Current version of s3tie.
*     ==================================================================
      integer                 maxtie
      parameter         (     maxtie = 69)
      character          ties(maxtie)*16
      logical            sorted
      parameter         (sorted =   .true.)
*     ------------------------------------------------------------------
      data
     &   ties(  1) /'(TOTAL)      '/,
     &   ties(  2) /'(USER)       '/,
     &   ties(  3) /'ALL          '/,
     &   ties(  4) /'BASIC        '/,
     &   ties(  5) /'BASIS        '/,
     &   ties(  6) /'BOUND        '/,
     &   ties(  7) /'CG           '/,
     &   ties(  8) /'CHARACTER    '/,
     &   ties(  9) /'CHOLESKY     '/,
     &   ties( 10) /'COLUMNS      '/,
     &   ties( 11) /'COMPLETE     '/,
     &   ties( 12) /'CONSTRAINTS  '/,
     &   ties( 13) /'DAMPING      '/,
     &   ties( 14) /'DEBUG        '/,
     &   ties( 15) /'DENSE        '/,
     &   ties( 16) /'DENSITY      '/,
     &   ties( 17) /'DERIVATIVE   '/,
     &   ties( 18) /'DIAGONAL     '/,
     &   ties( 19) /'DIFFERENCES  '/,
     &   ties( 20) /'DIMENSION    '/

      data
     &   ties( 21) /'ELEMENTS     '/,
     &   ties( 22) /'EXIT         '/,
     &   ties( 23) /'FACTORIZATION'/,
     &   ties( 24) /'FEASIBILITY  '/,
     &   ties( 25) /'FILE         '/,
     &   ties( 26) /'FLUSH        '/,
     &   ties( 27) /'FREQUENCY    '/,
     &   ties( 28) /'FULL         '/,
     &   ties( 29) /'GRADIENTS    '/,
     &   ties( 30) /'INFORMATION  '/,
     &   ties( 31) /'INTEGER      '/,
     &   ties( 32) /'ITERATIONS   '/,
     &   ties( 33) /'JACOBIAN     '/,
     &   ties( 34) /'LEVEL        '/,
     &   ties( 35) /'LIMITED      '/,
     &   ties( 36) /'LINEAR       '/,
     &   ties( 37) /'LINESEARCH   '/,
     &   ties( 38) /'LOG          '/,
     &   ties( 39) /'MODE         '/,
     &   ties( 40) /'NEWTON       '/

      data
     &   ties( 41) /'NO           '/,
     &   ties( 42) /'NONLINEAR    '/,
     &   ties( 43) /'OBJECTIVE    '/,
     &   ties( 44) /'OPTIMALITY   '/,
     &   ties( 45) /'OPTION       '/,
     &   ties( 46) /'PARTIAL      '/,
     &   ties( 47) /'PHASE1       '/,
     &   ties( 48) /'PHASE2       '/,
     &   ties( 49) /'POINT        '/,
     &   ties( 50) /'PRECONDITIONING'/,
     &   ties( 51) /'PRINT        '/,
     &   ties( 52) /'QN           '/,
     &   ties( 53) /'REAL         '/,
     &   ties( 54) /'ROOK         '/,
     &   ties( 55) /'ROW          '/,
     &   ties( 56) /'SINGULARITY  '/,
     &   ties( 57) /'SOLVER       '/,
     &   ties( 58) /'SPARSE       '/,
     &   ties( 59) /'START        '/,
     &   ties( 60) /'STEP         '/

      data
     &   ties( 61) /'STOP         '/,
     &   ties( 62) /'SUPERBASICS  '/,
     &   ties( 63) /'SWAP         '/,
     &   ties( 64) /'TOLERANCE    '/,
     &   ties( 65) /'UPDATES      '/,
     &   ties( 66) /'VARIABLES    '/,
     &   ties( 67) /'WEIGHT       '/,
     &   ties( 68) /'WEIGHTMAX    '/,
     &   ties( 69) /'YES          '/
*     ------------------------------------------------------------------
      call oplook( maxtie, ties, sorted, tie, loc )

      end ! subroutine s3tie

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3undf
     &   ( cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

*     ==================================================================
*     s3undf sets all options as undefined.
*
*     07 Feb 1998: First version of s3undf.
*     23 Nov 2002: Current version.
*     ==================================================================
      integer
     &     i, idummy
*     ------------------------------------------------------------------
      integer            first,       last
      parameter         (first  = 51, last = 180)
      double precision   rdummy
      character          cdummy*8
      parameter         (idummy =  -11111,  rdummy = -11111.0d+0)
      parameter         (cdummy = '-1111111'                    )
*     ------------------------------------------------------------------
      do i = first, last
         cw(i) = cdummy
         iw(i) = idummy
         rw(i) = rdummy
      end do

      end ! subroutine s3undf

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE oplook (NDICT, DICTRY, ALPHA, KEY, ENTRY)
C
C
C Description and usage:
C
C       Performs dictionary lookups.  A pointer is returned if a
C    match is found between the input key and the corresponding
C    initial characters of one of the elements of the dictionary.
C    If a "synonym" has been provided for an entry, the search is
C    continued until a match to a primary dictionary entry is found.
C    Cases of no match, or multiple matches, are also provided for.
C
C     Dictionary entries must be left-justified, and may be alphabetized
C    for faster searches.  Secondary entries, if any, are composed of
C    two words separated by one or more characters such as blank, tab,
C    comma, colon, or equal sign which are treated as non-significant
C    by opscan.  The first entry of each such pair serves as a synonym
C    for the second, more fundamental keyword.
C
C       The ordered search stops after the section of the dictionary
C    having the same first letters as the key has been checked, or
C    after a specified number of entries have been examined.  A special
C    dictionary entry, the vertical bar '|', will also terminate the
C    search.  This will speed things up if an appropriate dictionary
C    length parameter cannot be determined.  Both types of search are
C    sequential.  See "Notes" below for some suggestions if efficiency
C    is an issue.
C
C
C Parameters:
C
C    Name    Dimension  Type  I/O/S  Description
C    NDICT               I    I      Number of dictionary entries to be
C                                    examined.
C    DICTRY  NDICT       C    I      Array of dictionary entries,
C                                    left-justified in their fields.
C                                    May be alphabetized for efficiency,
C                                    in which case ALPHA should be
C                                    .TRUE.  Entries with synonyms are
C                                    of the form
C                                    'ENTRY : SYNONYM', where 'SYNONYM'
C                                    is a more fundamental entry in the
C                                    same dictionary.  NOTE: Don't build
C                                    "circular" dictionaries!
C    ALPHA               L    I      Indicates whether the dictionary
C                                    is in alphabetical order, in which
C                                    case the search can be terminated
C                                    sooner.
C    KEY                 C    I/O    String to be compared against the
C                                    dictionary.  Abbreviations are OK
C                                    if they correspond to a unique
C                                    entry in the dictionary.  KEY is
C                                    replaced on termination by its most
C                                    fundamental equivalent dictionary
C                                    entry (uppercase, left-justified)
C                                    if a match was found.
C    ENTRY               I      O    Dictionary pointer.  If > 0, it
C                                    indicates which entry matched KEY.
C                                    In case of trouble, a negative
C                                    value means that a UNIQUE match
C                                    was not found - the absolute value
C                                    of ENTRY points to the second
C                                    dictionary entry that matched KEY.
C                                    Zero means that NO match could be
C                                    found.  ENTRY always refers to the
C                                    last search performed -
C                                    in searching a chain of synonyms,
C                                    a non-positive value will be
C                                    returned if there is any break,
C                                    even if the original input key
C                                    was found.
C
C
C External references:
C
C    Name    Description
C    opscan  Finds first and last significant characters.
C
C
C Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C               Appears to satisfy the ANSI Fortran 77 standard.
C
C
C Notes:
C
C    (1)  IMPLICIT NONE is non-standard.
C
C    (2)  We have assumed that the dictionary is not too big.  If
C         many searches are to be done or if the dictionary has more
C         than a dozen or so entries, it may be advantageous to build
C         an index array of pointers to the beginning of the section
C         of the dictionary containing each letter, then pass in the
C         portion of the dictionary beginning with DICTRY (INDEX).
C         (This won't generally work for dictionaries with synonyms.)
C         For very large problems, a completely different approach may
C         be advisable, e.g. a binary search for ordered dictionaries.
C
C    (3)  oplook is case sensitive.  In most applications it will be
C         necessary to use an uppercase dictionary, and to convert the
C         input key to uppercase before calling oplook.  Companion
C         routines optokn and PAIRS, available from the author, already
C         take care of this.
C
C    (4)  The key need not be left-justified.  Any leading (or
C         trailing) characters which are "non-significant" to opscan
C         will be ignored.  These include blanks, horizontal tabs,
C         commas, colons, and equal signs.  See opscan for details.
C
C    (5)  The ASCII collating sequence for character data is assumed.
C         (N.B. This means the numerals precede the alphabet, unlike
C         common practice!)  This should not cause trouble on EBCDIC
C         machines if DICTRY just contains alphabetic keywords.
C         Otherwise it may be necessary to use the FORTRAN lexical
C         library routines to force use of the ASCII sequence.
C
C    (6)  Parameter NUMSIG sets a limit on the length of significant
C         dictionary entries.  Special applications may require that
C         this be increased.  (It is 16 in the present version.)
C
C    (7)  No protection against "circular" dictionaries is provided:
C         don't claim that A is B, and that B is A.  All synonym chains
C         must terminate!  Other potential errors not checked for
C         include duplicate or mis-ordered entries.
C
C    (8)  The handling of ambiguities introduces some ambiguity:
C
C            ALPHA = .TRUE.  A potential problem, when one entry
C                            looks like an abbreviation for another
C                            (eg. does 'A' match 'A' or 'AB'?) was
C                            resolved by dropping out of the search
C                            immediately when an "exact" match is found.
C
C            ALPHA = .FALSE. The programmer must ensure that the above
C                            situation does not arise: each dictionary
C                            entry must be recognizable, at least when
C                            specified to full length.  Otherwise, the
C                            result of a search will depend on the
C                            order of entries.
C
C
C Author:  Robert Kennelly, Informatics General Corporation.
C
C
C Development history:
C
C    24 Feb. 1984  RAK/DAS  Initial design and coding.
C    25 Feb. 1984    RAK    Combined the two searches by suitable
C                           choice of terminator FLAG.
C    28 Feb. 1984    RAK    Optional synonyms in dictionary, no
C                           longer update KEY.
C    29 Mar. 1984    RAK    Put back replacement of KEY by its
C                           corresponding entry.
C    21 June 1984    RAK    Corrected bug in error handling for cases
C                           where no match was found.
C    23 Apr. 1985    RAK    Introduced test for exact matches, which
C                           permits use of dictionary entries which
C                           would appear to be ambiguous (for ordered
C                           case).  Return -I to point to the entry
C                           which appeared ambiguous (had been -1).
C                           Repaired loop termination - had to use
C                           equal length strings or risk quitting too
C                           soon when one entry is an abbreviation
C                           for another.  Eliminated HIT, reduced
C                           NUMSIG to 16.
C    15 Nov. 1985    MAS    Loop 20 now tests .LT. FLAG, not .LE. FLAG.
C                           If ALPHA is false, FLAG is now '|', not '{'.
C    26 Jan. 1986    PEG    Declaration of FLAG and TARGET modified to
C                           conform to ANSI-77 standard.
C    05 Dec  2004    PEG    Used intrinsics LGE, LLE
C-----------------------------------------------------------------------


C     Variable declarations.
C     ----------------------

      IMPLICIT NONE

C     Parameters.

      INTEGER
     &   NUMSIG
      CHARACTER
     &   VBAR
      PARAMETER
     &   (VBAR = '|', NUMSIG = 16)

C     Variables.

      LOGICAL
     &   ALPHA
      INTEGER
     &   ENTRY, FIRST, I, LAST, LENGTH, MARK, NDICT
*     CHARACTER
*    &   DICTRY (NDICT) * (*), FLAG * (NUMSIG),
*    &   KEY * (*), TARGET * (NUMSIG)
      CHARACTER
     &   DICTRY (NDICT) * (*), FLAG * 16,
     &   KEY * (*), TARGET * 16

C     Procedures.

      EXTERNAL
     &   opscan


C     Executable statements.
C     ----------------------

      ENTRY = 0

C     Isolate the significant portion of the input key (if any).

      FIRST = 1
      LAST  = MIN( LEN(KEY), NUMSIG )
      CALL opscan (KEY, FIRST, LAST, MARK)

      IF (MARK .GT. 0) THEN
         TARGET = KEY (FIRST:MARK)

C        Look up TARGET in the dictionary.

   10    CONTINUE
            LENGTH = MARK - FIRST + 1

C           Select search strategy by cunning choice of termination test
C           flag.  The vertical bar is just about last in both the
C           ASCII and EBCDIC collating sequences.

            IF (ALPHA) THEN
               FLAG = TARGET
            ELSE
               FLAG = VBAR
            END IF


C           Perform search.
C           ---------------

            I = 0
   20       CONTINUE
               I = I + 1
               IF (TARGET (1:LENGTH) .EQ. DICTRY (I) (1:LENGTH)) THEN
                  IF (ENTRY .EQ. 0) THEN

C                    First "hit" - must still guard against ambiguities
C                    by searching until we've gone beyond the key
C                    (ordered dictionary) or until the end-of-dictionary
C                    mark is reached (exhaustive search).

                     ENTRY = I

C                    Special handling if match is exact - terminate
C                    search.  We thus avoid confusion if one dictionary
C                    entry looks like an abbreviation of another.
C                    This fix won't generally work for un-ordered
C                    dictionaries!

                     FIRST = 1
                     LAST = NUMSIG
                     CALL opscan (DICTRY (ENTRY), FIRST, LAST, MARK)
                     IF (MARK .EQ. LENGTH) I = NDICT
                  ELSE


C                    Oops - two hits!  Abnormal termination.
C                    ---------------------------------------

                     ENTRY = -I
                     RETURN
                  END IF
               END IF

C           Check whether we've gone past the appropriate section of the
C           dictionary.  The test on the index provides insurance and an
C           optional means for limiting the extent of the search.

            IF (LLT(DICTRY (I) (1:LENGTH), FLAG)  .AND.  I .LT. NDICT)
     &         GO TO 20


C           Check for a synonym.
C           --------------------

            IF (ENTRY .GT. 0) THEN

C              Look for a second entry "behind" the first entry.  FIRST
C              and MARK were determined above when the hit was detected.

               FIRST = MARK + 2
               CALL opscan (DICTRY (ENTRY), FIRST, LAST, MARK)
               IF (MARK .GT. 0) THEN

C                 Re-set target and dictionary pointer, then repeat the
C                 search for the synonym instead of the original key.

                  TARGET = DICTRY (ENTRY) (FIRST:MARK)
                  ENTRY = 0
                  GO TO 10

               END IF
            END IF

      END IF
      IF (ENTRY .GT. 0) KEY = DICTRY (ENTRY)


C     Normal termination.
C     -------------------

      RETURN

C     End of oplook
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      FUNCTION opnumb( STRING )

      LOGICAL          opnumb
      CHARACTER*(*)    STRING

************************************************************************
*     Description and usage:
*
*        A simple(-minded) test for numeric data is implemented by
*        searching an input string for legitimate characters:
*                digits 0 to 9, D, E, -, + and .
*        Insurance is provided by requiring that a numeric string
*        have at least one digit, at most one D, E or .
*        and at most two -s or +s.  Note that a few ambiguities remain:
*
*           (a)  A string might have the form of numeric data but be
*                intended as text.  No general test can hope to detect
*                such cases.
*
*           (b)  There is no check for correctness of the data format.
*                For example a meaningless string such as 'E1.+2-'
*                will be accepted as numeric.
*
*        Despite these weaknesses, the method should work in the
*        majority of cases.
*
*
*     Parameters:
*
*        Name    Dimension  Type  I/O/S  Description
*        opnumb              L      O    Set .TRUE. if STRING appears
*                                        to be numerical data.
*        STRING              C    I      Input data to be tested.
*
*
*     Environment:  ANSI FORTRAN 77.
*
*
*     Notes:
*
*        (1)  It is assumed that STRING is a token extracted by
*             optokn, which will have converted any lower-case
*             characters to upper-case.
*
*        (2)  optokn pads STRING with blanks, so that a genuine
*             number is of the form  '1234        '.
*             Hence, the scan of STRING stops at the first blank.
*
*        (3)  COMPLEX data with parentheses will not look numeric.
*
*
*     Systems Optimization Laboratory, Stanford University.
*     12 Nov  1985    Initial design and coding, starting from the
*                     routine ALPHA from Informatics General, Inc.
*     05 Dec  2004    Used intrinsics LGE, LLE
************************************************************************

      LOGICAL         NUMBER
      INTEGER         J, LENGTH, NDIGIT, NEXP, NMINUS, NPLUS, NPOINT
      CHARACTER       ATOM*1

      NDIGIT = 0
      NEXP   = 0
      NMINUS = 0
      NPLUS  = 0
      NPOINT = 0
      NUMBER = .TRUE.
      LENGTH = LEN (STRING)
      J      = 0

   10    J    = J + 1
         ATOM = STRING (J:J)
         IF      (LGE(ATOM, '0')  .AND.  LLE(ATOM, '9')) THEN
            NDIGIT = NDIGIT + 1
         ELSE IF (ATOM .EQ. 'D'   .OR.   ATOM .EQ. 'E' ) THEN
            NEXP   = NEXP   + 1
         ELSE IF (ATOM .EQ. '-') THEN
            NMINUS = NMINUS + 1
         ELSE IF (ATOM .EQ. '+') THEN
            NPLUS  = NPLUS  + 1
         ELSE IF (ATOM .EQ. '.') THEN
            NPOINT = NPOINT + 1
         ELSE IF (ATOM .EQ. ' ') THEN
            J      = LENGTH
         ELSE
            NUMBER = .FALSE.
         END IF

         IF (NUMBER  .AND.  J .LT. LENGTH) GO TO 10

      opnumb = NUMBER
     &         .AND.  NDIGIT .GE. 1
     &         .AND.  NEXP   .LE. 1
     &         .AND.  NMINUS .LE. 2
     &         .AND.  NPLUS  .LE. 2
     &         .AND.  NPOINT .LE. 1

      RETURN

*     End of opnumb
      END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE opscan (STRING, FIRST, LAST, MARK)

      implicit           none
      character*(*)      STRING
      integer            FIRST, LAST, MARK
C
C
C Description and usage:
C
C       Looks for non-blank fields ("tokens") in a string, where the
C    fields are of arbitrary length, separated by blanks, tabs, commas,
C    colons, or equal signs.  The position of the end of the 1st token
C    is also returned, so this routine may be conveniently used within
C    a loop to process an entire line of text.
C
C       The procedure examines a substring, STRING (FIRST : LAST), which
C    may of course be the entire string (in which case just call opscan
C    with FIRST <= 1 and LAST >= LEN (STRING) ).  The indices returned
C    are relative to STRING itself, not the substring.
C
C
C Parameters:
C
C    Name    Dimension  Type  I/O/S  Description
C    STRING              C    I      Text string containing data to be
C                                    scanned.
C    FIRST               I    I/O    Index of beginning of substring.
C                                    If <= 1, the search begins with 1.
C                                    Output is index of beginning of
C                                    first non-blank field, or 0 if no
C                                    token was found.
C    LAST                I    I/O    Index of end of substring.
C                                    If >= LEN (STRING), the search
C                                    begins with LEN (STRING).  Output
C                                    is index of end of last non-blank
C                                    field, or 0 if no token was found.
C    MARK                I      O    Points to end of first non-blank
C                                    field in the specified substring.
C                                    Set to 0 if no token was found.
C
C
C Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C               ANSI Fortran 77, except for the tab character HT.
C
C Notes:
C
C    (1)  IMPLICIT NONE is non-standard.  Constant HT (Tab) is defined
C         in a non-standard way:  the CHAR function is not permitted
C         in a PARAMETER declaration (OK on VAX, though).  For Absoft
C         FORTRAN 77 on 68000 machines, use HT = 9.  In other cases, it
C         may be best to declare HT as a variable and assign
C         HT = CHAR(9) on ASCII machines, or CHAR(5) for EBCDIC.
C
C    (2)  The pseudo-recursive structure was chosen for fun.  It is
C         equivalent to three DO loops with embedded GO TOs in sequence.
C
C    (3)  The variety of separators recognized limits the usefulness of
C         this routine somewhat.  The intent is to facilitate handling
C         such tokens as keywords or numerical values.  In other
C         applications, it may be necessary for ALL printing characters
C         to be significant.  A simple modification to statement
C         function SOLID will do the trick.
C
C
C Author:  Robert Kennelly, Informatics General Corporation.
C
C
C Development history:
C
C    29 Dec. 1984    RAK    Initial design and coding, (very) loosely
C                           based on SCAN_STRING by Ralph Carmichael.
C    25 Feb. 1984    RAK    Added ':' and '=' to list of separators.
C    16 Apr. 1985    RAK    Defined SOLID in terms of variable DUMMY
C                           (previous re-use of STRING was ambiguous).
C
C-----------------------------------------------------------------------


C     Variable declarations.
C     ----------------------

C     Parameters.

      CHARACTER
     &   BLANK, EQUAL, COLON, COMMA, HT
      PARAMETER
     &   (BLANK = ' ', EQUAL = '=', COLON = ':', COMMA = ',')

C     Variables.

!     LOGICAL
!    &   SOLID
      INTEGER
     &   BEGIN, END, LENGTH
      CHARACTER
     &   c

!     Statement functions.
!
!     SOLID (DUMMY) = (DUMMY .NE. BLANK) .AND.
!    &                (DUMMY .NE. COLON) .AND.
!    &                (DUMMY .NE. COMMA) .AND.
!    &                (DUMMY .NE. EQUAL) .AND.
!    &                (DUMMY .NE. HT)


C     Executable statements.
C     ----------------------

****  HT     = CHAR(9) for ASCII machines, CHAR(5) for EBCDIC.
      HT     = CHAR(9)
      MARK   = 0
      LENGTH = LEN (STRING)
      BEGIN  = MAX (FIRST, 1)
      END    = MIN (LENGTH, LAST)

C     Find the first significant character ...

      DO 30 FIRST = BEGIN, END, +1
         c = STRING (FIRST : FIRST)
       ! IF ( SOLID(c) ) THEN
         if ( c .ne. BLANK  .and.  c .ne. COLON  .and.
     &        c .ne. COMMA  .and.  c .ne. EQUAL  .and.
     &        c .ne. HT   ) then

C           ... then the end of the first token ...

            DO 20 MARK = FIRST, END - 1, +1
               c = STRING (MARK + 1 : MARK + 1)
             ! IF (.NOT.SOLID(c) ) THEN
               if ( c .ne. BLANK  .and.  c .ne. COLON  .and.
     &              c .ne. COMMA  .and.  c .ne. EQUAL  .and.
     &              c .ne. HT   ) then
                  ! relax
               else
C                 ... and finally the last significant character.

                  DO 10 LAST = END, MARK, -1
                     c = STRING (LAST : LAST)
                   ! IF ( SOLID(c) ) THEN
                     if ( c .ne. BLANK  .and.  c .ne. COLON  .and.
     &                    c .ne. COMMA  .and.  c .ne. EQUAL  .and.
     &                    c .ne. HT   ) then
                        RETURN
                     END IF
   10             CONTINUE

C                 Everything past the first token was a separator.

                  LAST = LAST + 1
                  RETURN
               END IF
   20       CONTINUE

C           There was nothing past the first token.

            LAST = MARK
            RETURN
         END IF
   30 CONTINUE

C     Whoops - the entire substring STRING (BEGIN : END) was composed of
C     separators !

      FIRST = 0
      MARK = 0
      LAST = 0
      RETURN

C     End of opscan
      END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE optokn (STRING, NUMIN, NUMOUT, LIST)
C
C
C Description and usage:
C
C       An aid to parsing input data.  The individual "tokens" in a
C    character string are isolated, converted to uppercase, and stored
C    in an array.  Here, a token is a group of significant, contiguous
C    characters.  The following are NON-significant, and hence may
C    serve as separators:  blanks, horizontal tabs, commas, colons,
C    and equal signs.  See opscan for details.  Processing continues
C    until the requested number of tokens have been found or the end
C    of the input string is reached.
C
C
C Parameters:
C
C    Name    Dimension  Type  I/O/S  Description
C    STRING              C    I      Input string to be analyzed.
C    NUMIN               I    I/O    Number of tokens requested (input)
C                                    and found (output).
C    LIST    NUMIN       C      O    Array of tokens, changed to upper
C                                    case.
C
C
C External references:
C
C    Name    Description
C    opscan  Finds positions of first and last significant characters.
C    opuppr  Converts a string to uppercase.
C
C
C Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C               Appears to satisfy the ANSI Fortran 77 standard.
C
C
C Notes:
C
C    (1)  IMPLICIT NONE is non-standard.
C
C
C Author:  Robert Kennelly, Informatics General Corporation.
C
C
C Development history:
C
C    16 Jan. 1984    RAK    Initial design and coding.
C    16 Mar. 1984    RAK    Revised header to reflect full list of
C                           separators, repaired faulty WHILE clause
C                           in "10" loop.
C    18 Sep. 1984    RAK    Change elements of LIST to uppercase one
C                           at a time, leaving STRING unchanged.
C    05 Dec. 2004    PEG    Replaced by NUMIN, NUMOUT
C
C-----------------------------------------------------------------------


C     Variable declarations.
C     ----------------------

      IMPLICIT NONE

C     Parameters.

      CHARACTER
     &   BLANK
      PARAMETER
     &   (BLANK = ' ')

C     Variables.

      INTEGER
     &   COUNT, FIRST, I, LAST, MARK, NUMIN, NUMOUT
      CHARACTER
     &   STRING * (*), LIST (NUMIN) * (*)

C     Procedures.

      EXTERNAL
     &   opuppr, opscan


C     Executable statements.
C     ----------------------

C     WHILE there are tokens to find, loop UNTIL enough have been found.

      FIRST = 1
      LAST = LEN (STRING)

      COUNT = 0
   10 CONTINUE

C        Get delimiting indices of next token, if any.

         CALL opscan (STRING, FIRST, LAST, MARK)
         IF (LAST .GT. 0) THEN
            COUNT = COUNT + 1

C           Pass token to output string array, then change case.

            LIST (COUNT) = STRING (FIRST : MARK)
            CALL opuppr (LIST (COUNT))
            FIRST = MARK + 2
            IF (COUNT .LT. NUMIN) GO TO 10

         END IF


C     Fill the rest of LIST with blanks and set NUMOUT for output.

      DO 20 I = COUNT + 1, NUMIN
         LIST (I) = BLANK
   20 CONTINUE

      NUMOUT = COUNT


C     Termination.
C     ------------

      RETURN

C     End of optokn
      END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE opuppr(STRING)
C
C ACRONYM:  UPper CASE
C
C PURPOSE:  This subroutine changes all lower case letters in the
C           character string to upper case.
C
C METHOD:   Each character in STRING is treated in turn.  The intrinsic
C           function INDEX effectively allows a table lookup, with
C           the local strings LOW and UPP acting as two tables.
C           This method avoids the use of CHAR and ICHAR, which appear
C           be different on ASCII and EBCDIC machines.
C
C ARGUMENTS
C    ARG       DIM     TYPE I/O/S DESCRIPTION
C  STRING       *       C   I/O   Character string possibly containing
C                                 some lower-case letters on input;
C                                 strictly upper-case letters on output
C                                 with no change to any non-alphabetic
C                                 characters.
C
C EXTERNAL REFERENCES:
C  LEN    - Returns the declared length of a CHARACTER variable.
C  INDEX  - Returns the position of second string within first.
C
C ENVIRONMENT:  ANSI FORTRAN 77
C
C DEVELOPMENT HISTORY:
C     DATE  INITIALS  DESCRIPTION
C   06/28/83   CLH    Initial design.
C   01/03/84   RAK    Eliminated NCHAR input.
C   06/14/84   RAK    Used integer PARAMETERs in comparison.
C   04/21/85   RAK    Eliminated DO/END DO in favor of standard code.
C   09/10/85   MAS    Eliminated CHAR,ICHAR in favor of LOW, UPP, INDEX.
C   12/04/04   PEG    Used intrinsics LGE, LLE
C
C AUTHOR: Charles Hooper, Informatics General, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      CHARACTER      STRING * (*)
      INTEGER        I, J
      CHARACTER      C*1, LOW*26, UPP*26
      DATA           LOW /'abcdefghijklmnopqrstuvwxyz'/,
     &               UPP /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      DO 10 J = 1, LEN(STRING)
         C    = STRING(J:J)
         IF (LGE(C, 'a')  .AND.  LLE(C, 'z')) THEN
            I           = INDEX( LOW, C )
            IF (I .GT. 0) STRING(J:J) = UPP(I:I)
         END IF
   10 CONTINUE
      RETURN

*     End of opuppr
      END
