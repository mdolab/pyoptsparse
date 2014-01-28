*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  sn02lib.f
*
*     snTitl   snInit   snSpec   snchkA   snJac
*     snEXIT   snWRAP   snSolF
*     snMem    snMemA   snMemA0  snMemB
*     snLog    snLog2   snSTOP
*     snSet    snSeti   snSetr
*     snGet    snGetc   snGeti   snGetr
*     snRetH
*
*     09 Mar 2004: snSolF implemented.
*     17 Jun 2004: snSolF always flags infeasible jbInf1 as I.
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snTitl( title )

      character
     &     title*30

*     ==================================================================
*     snTitl sets the title for snopt.
*     ==================================================================

      title  = 'S N O P T  7.2-5    (May 2007)'
*---------------123456789|123456789|123456789|--------------------------

      end ! subroutine snTitl

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iPrint, iSumm, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

*     ==================================================================
*     snInit  is called by the user to do the following:
*     1. Open default files (Print, Summary).
*     2. Initialize title.
*     3. Set options to default values.
*
*     15 Nov 1991: First version.
*     14 Jul 1997: Thread-safe version.
*     02 Oct 1997: Character workspace added.
*     15 Oct 2003: Current version of snInit.
*     ==================================================================
      external
     &     s1outpt
      character
     &     Solver*6, str*80, str2*80, title*30
      integer
     &     inform, iSpecs, iStdo, maxru, maxrw, maxiu, maxiw, maxcu,
     &     maxcw, nnCon, nnJac, nnL, nnObj, lvlTim, s1outpt
*     ------------------------------------------------------------------
      parameter         (maxru     =   2) ! start of SNOPT part of rw
      parameter         (maxrw     =   3) ! end   of SNOPT part of rw
      parameter         (maxiu     =   4) ! start of SNOPT part of iw
      parameter         (maxiw     =   5) ! end   of SNOPT part of iw
      parameter         (maxcu     =   6) ! start of SNOPT part of cw
      parameter         (maxcw     =   7) ! end   of SNOPT part of cw
      parameter         (nnJac     =  21) ! # nonlinear Jac, variables
      parameter         (nnObj     =  22) ! # variables in gObj
      parameter         (nnCon     =  23) ! # of nonlinear constraints
      parameter         (nnL       =  24) ! nonlinear vars
      parameter         (lvlTim    = 182) ! Timing level
*     ------------------------------------------------------------------
      character          dashes*30
      data               dashes /'=============================='/
*     ------------------------------------------------------------------
      Solver = 'SNINIT'

      if (lencw .lt. 500 .or. leniw .lt. 500 .or. lenrw .lt. 500) then
*        ---------------------------------------------------------------
*        Not enough workspace to do ANYTHING!
*        Print and exit without accessing the work arrays.
*        ---------------------------------------------------------------
         inform = 81       ! Work arrays must have at least 500 elements
         call snWRAP( inform, Solver, str, str2, iw, leniw )
         go to 999
      end if

      iSpecs    = 0
      iStdo     = s1outpt( )
      iw( 10)   = iStdo   ! Standard Output
      iw( 11)   = iSpecs  ! Specs file (default)
      iw( 12)   = iPrint  ! Print file
      iw( 13)   = iSumm   ! Summary file

      iw(maxcu) = 500
      iw(maxiu) = 500
      iw(maxru) = 500
      iw(maxcw) = lencw
      iw(maxiw) = leniw
      iw(maxrw) = lenrw

*     These dimensions need to be initialized for an MPS run.

      iw(nnCon) = 0
      iw(nnJac) = 0
      iw(nnObj) = 0
      iw(nnL  ) = 0

      call snTitl( title )
      call s1init( title, iw, leniw, rw, lenrw )

      call snPRNT(11, '         '//dashes, iw, leniw )
      call snPRNT( 1, '         '//title , iw, leniw )
      call snPRNT( 1, '         '//dashes, iw, leniw )

      call snPRNT(12, ' '//dashes, iw, leniw )
      call snPRNT( 2, ' '//title , iw, leniw )
      call snPRNT( 2, ' '//dashes, iw, leniw )

*     ------------------------------------------------------------------
*     Set the options to default values.
*     snopt  will check the options later and maybe print them.
*     ------------------------------------------------------------------
      call s3undf( cw, lencw, iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     Initialize some global values.
*     ------------------------------------------------------------------
      iw(lvlTim) = 3

  999 return

      end ! subroutine snInit

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snSpec
     &   ( iSpecs, iExit, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iSpecs, iExit, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

*     ==================================================================
*     snSpec  may be called by the user to read a Specs file.
*
*     07 Feb 1998: First version of snSpec.
*     01 Aug 2003: s3file now has a "title" parameter.  Use ' '.
*     13 Jul 2005: Included error count in value of iExit.
*     22 Apr 2007: Exit 107 added (some keywords unrecognized)
*     ==================================================================
      character
     &     Solver*6, str*80, str2*80
      integer
     &     Errors, iPrint, iSumm, Calls
      external
     &     s3opt
*     ------------------------------------------------------------------
      Solver = 'SNSPEC'

      if (lencw .lt. 500 .or. leniw .lt. 500 .or. lenrw .lt. 500) then
*        ---------------------------------------------------------------
*        Not enough workspace to do ANYTHING!
*        Print and exit without accessing the work arrays.
*        ---------------------------------------------------------------
         iExit = 81        ! Work arrays must have at least 500 elements
         call snWRAP( iExit, Solver, str, str2, iw, leniw )
         go to 999
      end if

      if (iSpecs .le. 0  .or.  iSpecs .gt. 99) then
         iExit = 131      ! iSPECS out of range
         go to 800
      end if

      iw( 11)   = iSpecs  ! Specs (options) file

      iPrint    = iw( 12) ! Print file
      iSumm     = iw( 13) ! Summary file

      Calls     = 1

*     ------------------------------------------------------------------
*     Read the Specs file.
*     snopt  will check the options later and maybe print them.
*     ------------------------------------------------------------------
      call s3file
     &   ( iExit, Calls, iSpecs, s3opt, ' ', iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

  800 if (iExit .eq. 0) then
         if (Errors .eq. 0) then
            iExit = 101         ! SPECS file read successfully
         else
            iExit = 107         ! some SPECS keywords not recognized
         end if
      end if

      call snWRAP( iExit, Solver, str, str2, iw, leniw )

  999 return

      end ! subroutine snSpec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snchkA
     &   ( iExit, nF, n, lvlChk, userfg,
     &     iGfun, jGvar, lenG, neG, x,
     &     mincw, miniw, minrw,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     userfg
      integer
     &     iExit, lenG, lencu, lencw, leniu, leniw, lenru, lenrw,
     &     lvlChk, mincw, miniw, minrw, neG, nF, n,
     &     iGfun(lenG), jGvar(lenG), iu(leniu), iw(leniw)
      double precision
     &     x(n), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     snchkA  is a stand-alone derivative checker for problems in
!     snOptA format.
!
!     snchkA is different in two ways from the built-in derivative
!     checker in snopta:

!       1. The derivatives are checked before the variables and
!          constraints are reordered to conform with snoptb format.
!
!       2. The derivatives are checked at the point defined by the
!          input argument x.  A feasible point is NOT computed before
!          the check is done.
!
!
!     lvlChk has the following meaning:
!
!       -1         do not perform any check.
!        0         do the cheap test only.
!       >0         do both cheap and full test on problem derivatives.
!
!     ------------------------------------------------------------------
!     NOTE: Before calling snchkA, there MUST be a call to the
!     initialization routine, i.e.,
!     call snInit( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )
!     This sets the default values of the optional parameters. You can
!     also alter the default values of iPrint and iSumm before snchkA
!     is used.  iPrint = 0, etc, is OK.
!     ------------------------------------------------------------------
!
!     EXIT snchkA 100 -- completed successfully
!     EXIT INFO   105 -- user-supplied derivatives appear to be correct
!     EXIT INFO   106 -- no derivatives were checked
!
!     EXIT snchkA  50 -- error in the user-supplied functions
!     EXIT INFO    55 -- incorrect derivatives
!
!     EXIT snchkA  70 -- user requested termination
!     EXIT INFO    71 -- terminated during function evaluation
!
!     EXIT snchkA  80 -- insufficient storage allocated
!     EXIT INFO    81 -- work arrays must have at least 500 elements
!     EXIT INFO    83 -- not enough integer storage
!
!     SNOPT package maintained by Philip E. Gill,
!     Dept of Mathematics, University of California, San Diego.
!
!     01 Jun 2006: First version of snchkA
!     22 Apr 2007: INFO 106 added
!     ==================================================================
      character
     &     Solver*6, str*80, str2*80
      integer
     &     iGmax, lElem, lF, lF1, lG, lG1, lSet,
     &     lw, lw1, lx1, ly, lz, maxcw, maxiw,
     &     maxrw, nextcw, nextiw, nextrw, Status
      double precision
     &     eGmax
!     ------------------------------------------------------------------
      double precision   ok
      parameter         (ok   = 0.1d+0)
!     ------------------------------------------------------------------
      Solver = 'SNCHKA'
      iExit  = 0
      Status = 1

!     ------------------------------------------------------------------
!     Check memory limits and fetch the workspace starting positions.
!     ------------------------------------------------------------------
      call s2Mem0
     &   ( iExit, Solver, lencw, leniw, lenrw, iw,
     &     mincw, miniw, minrw, maxcw, maxiw, maxrw,
     &     nextcw, nextiw, nextrw )
      if (iExit .ne. 0) go to 999

      lElem  = nextiw
      lSet   = lElem  + nF
      miniw  = lSet   + n  - 1

      lF     = nextrw
      lG     = lF     + nF
      lx1    = lG     + lenG
      lF1    = lx1    + n
      lG1    = lF1    + nF
      lw     = lG1    + lenG
      lw1    = lw     + nF
      ly     = lw1    + nF
      lz     = ly     + n
      minrw  = lz     + n  - 1

      if (miniw .gt. maxiw  .or.  minrw .gt. maxrw) then
!        ---------------------------------------------------------------
!        Not enough space to check the derivatives.
!        Exit with an (over) estimate of the additional space needed.
!        ---------------------------------------------------------------
         if (miniw .gt. maxiw) then
            write(str, 9010) miniw
            call snPRNT( 11, str, iw, leniw )
            iExit = 83
         end if

         if (minrw .gt. maxrw) then
            write(str, 9020) minrw
            call snPRNT( 11, str, iw, leniw )
            iExit = 84
         end if
         go to 800
      end if

!     ------------------------------------------------------------------
!     Go for it.
!     ------------------------------------------------------------------
      call s7chkA
     &   ( iExit, Status, lvlChk, userfg,
     &     nF, n, iGmax, eGmax, iw(lElem), iw(lSet),
     &     iGfun, jGvar, lenG, neG,
     &     x, rw(lF), rw(lG),
     &     rw(lx1), rw(lF1), rw(lG1), rw(lw), rw(lw1), rw(ly),
     &     iw, leniw, rw, lenrw, cu, lencu,
     &     iu, leniu, ru, lenru )
      if (iExit .ne. 0) go to 800

!     Print the exit conditions.

  800 if (iExit .eq. 0) then
         iExit = 105            ! all derivatives appear to be correct
      end if

      call snWRAP( iExit, Solver, str, str2, iw, leniw )

  999 return

 2000 format(' -->  largest error in the estimated Jacobian is',
     &                1p, e12.2, '  in row', i6)
 9010 format(' Total integer   workspace  should be significantly',
     &       ' more than', i8)
 9020 format(' Total real      workspace  should be significantly',
     &       ' more than', i8)

      end ! subroutine snchkA

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snJac
     &   ( iExit, nF, n, userfg,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     x, xlow, xupp, mincw, miniw, minrw,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     userfg
      integer
     &     iExit, nF, n, neA, lenA, neG, lenG, mincw,
     &     miniw, minrw, lencu, lencw, leniu, leniw, lenru, lenrw,
     &     iAfun(lenA), jAvar(lenA), iGfun(lenG), jGvar(lenG),
     &     iu(leniu), iw(leniw)
      double precision
     &     A(lenA), x(n), xlow(n), xupp(n), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     snJac  computes the coordinates of the Jacobian.
*
*     All calculations are based on a point defined by moving the input
*     x  inside its upper and lower bounds.  snJac is terminated if
*     the problem functions are undefined at this point.
*
*     ------------------------------------------------------------------
*     NOTE: Before calling snJac, your calling program MUST call the
*     initialization routine using the call:
*     call snInit( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )
*     This sets the default values of the optional parameters. You can
*     also alter the default values of iPrint and iSumm before snJac
*     is used.  iPrint = 0, etc, is OK.
*     ------------------------------------------------------------------
*
*     EXIT snJac  100 -- completed successfully (for auxiliary routines)
*     EXIT INFO   102 -- Jacobian structure estimated
*
*     EXIT snJac  12x -- Errors while estimating Jacobian structure
*     EXIT INFO   121 -- cannot estimate Jacobian structure at given point
*     EXIT INFO    71 -- terminated during function evaluation
*     EXIT INFO    81 -- work arrays must have at least 500 elements
*     EXIT INFO    83 -- not enough integer storage
*
*     SNOPT package maintained by Philip E. Gill,
*     Dept of Mathematics, University of California, San Diego.
*
*     26 Oct 2002: First version
*     27 Sep 2003: More thorough checks for feasibility
*     12 Jul 2005: Current version of snJac
*     ==================================================================
      character
     &     Solver*6, str*80, str2*80
      logical
     &     prtAll
      integer
     &     DerOpt, imaxJ, iPrt, iSum, lF, lFw, lFy, lFz, lG, lGcolw,
     &     lGcoly, lGcolz, lw, lx, lxlow, lxupp, ly, lz, maxcw, maxiw,
     &     maxrw, nextcw, nextiw, nextrw, coltyp, rowtyp, Status
      double precision
     &     InfBnd, emaxJ
*     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero = 0.0d+0)
      double precision   ok
      parameter         (ok   = 0.1d+0)
*     ------------------------------------------------------------------
      Solver = 'SNJAC '
      iExit  = 0
      Status = 1

*     ------------------------------------------------------------------
*     Check memory limits and fetch the workspace starting positions.
*     ------------------------------------------------------------------
      call s2Mem0
     &   ( iExit, Solver, lencw, leniw, lenrw, iw,
     &     mincw, miniw, minrw, maxcw, maxiw, maxrw,
     &     nextcw, nextiw, nextrw )
      if (iExit .ne. 0) go to 999

      rowtyp = nextiw
      coltyp = rowtyp + nF
      miniw  = coltyp + n  - 1

      lF     = nextrw
      lG     = lF     + nF
      lx     = lG     + lenG
      lxlow  = lx     + n
      lxupp  = lxlow  + n
      lw     = lxupp  + n
      ly     = lw     + n
      lz     = ly     + n
      lFw    = lz     + n
      lFy    = lFw    + nF
      lFz    = lFy    + nF
      lGcolw = lFz    + nF
      lGcoly = lGcolw + nF
      lGcolz = lGcoly + nF
      minrw  = lGcolz + nF - 1

      if (miniw .gt. maxiw  .or.  minrw .gt. maxrw) then
*        ---------------------------------------------------------------
*        Not enough space to build the Jacobian.
*        Provide the user an (over) estimate of what is needed.
*        ---------------------------------------------------------------
         if (miniw .gt. maxiw) then
            write(str, 9010) miniw
            call snPRNT( 11, str, iw, leniw )
            iExit = 83
         end if

         if (minrw .gt. maxrw) then
            write(str, 9020) minrw
            call snPRNT( 11, str, iw, leniw )
            iExit = 84
         end if
         go to 800
      end if

*     ------------------------------------------------------------------
*     Go for it.
*     ------------------------------------------------------------------
      prtAll = .true.           ! ignore fixed variables for diffs

      InfBnd = rw( 70)          ! definition of an infinite bound
      if (InfBnd .lt. zero) then
         InfBnd = 1.0d+20       ! User hasn't assigned it
      end if

*     ------------------------------------------------------------------
*     Make sure that snOptA does not let userf compute derivatives.
*     The parameters iPrt and iSum may refer to the Print and Summary
*     file respectively.  Setting them to 0 suppresses printing.
*     ------------------------------------------------------------------
      DerOpt = 0
      iPrt   = 0
      iSum   = 0
      call snSeti
     &   ( 'Derivative option', DerOpt, iPrt, iSum, iExit,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call dcopy
     &   ( n, x, 1, rw(lx), 1 )
      call dcopy
     &   ( n, xlow, 1, rw(lxlow), 1 )
      call dcopy
     &   ( n, xupp, 1, rw(lxupp), 1 )

      call s7Jac
     &   ( iExit, Status, userfg, prtAll,
     &     nF, n, InfBnd, imaxJ, emaxJ,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     iw(rowtyp), iw(coltyp), rw(lx), rw(lxlow), rw(lxupp),
     &     rw(lF), rw(lG), rw(lw), rw(ly), rw(lz),
     &     rw(lFw), rw(lFy), rw(lFz), rw(lGcolw), rw(lGcoly),rw(lGcolz),
     &     iw, leniw, cu, lencu, iu, leniu, ru, lenru )
      if (iExit .ne. 0) go to 800

      write(str, 1000) neA+neG
      call snPRNT( 13, str, iw, leniw )
      write(str, 1010) neG, neA
      call snPRNT(  3, str, iw, leniw )

      if (emaxJ .gt. ok) then
         iExit = 121            ! unable to estimate Jacobian structure
         write(str, 2000) emaxJ, imaxJ
         call snPRNT( 11, str, iw, leniw )
      end if

*     Print the exit conditions.

  800 if (iExit .eq. 0) then
         iExit = 102            ! Jacobian structure estimated
      end if

      call snWRAP( iExit, Solver, str, str2, iw, leniw )

  999 return

 1000 format(' Nonzero derivs  Jij  ',  i8)
 1010 format(' Non-constant    Jij''s', i8, 5x, 'Constant Jij''s',4x,i8)
 2000 format(' -->  largest error in the estimated Jacobian is',
     &                1p, e12.2, '  in row', i6)
 9010 format(' Total integer   workspace  should be significantly',
     &       ' more than', i8)
 9020 format(' Total real      workspace  should be significantly',
     &       ' more than', i8)

      end ! subroutine snJac

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snEXIT
     &   ( iExit, Solver, string, string2 )

      implicit
     &     none
      character*(*)
     &     string, string2
      character
     &     Solver*6
      integer
     &     iExit

*     ==================================================================
*     snEXIT  returns the strings associated with EXIT condition iExit.
*
*     On exit, string1 and string2 are trimmed of trailing blanks and
*              each have a maximum length of 76 chars.
*
*     25 Sep 2002: First version of snEXIT.
*     12 Mar 2004: Trimmed message strings to avoid trouble on SGI's.
*     07 May 2006: Second derivative exits added (Exit 54)
*     01 Jun 2006: Stand-alone derivative checker added (Exits 55,105)
*     01 Jun 2006: Student edition exit added (Exit 66)
*     22 Apr 2007: Unrecognized options flagged (Exit 135)
*     22 Apr 2007: No derivatives checked (Exit 106)
*     ==================================================================
      integer
     &     i1, i2, mjr, mnr
*     ------------------------------------------------------------------
      integer
     &     length
      integer            Msgs
      parameter         (Msgs = 69)
      integer
     &     indc(0:14)
      character
     &     c(0:Msgs)*80
      data
     &     indc
     &     / 0, 5, 11, 14, 18, 23, 29, 33, 38, 43, 47, 55, 59, 61, 66 /
      data
     & c( 0)/'finished successfully'/,                             ! EXIT  0
     & c( 1)   /'optimality conditions satisfied'/,                ! EXIT  1
     & c( 2)   /'feasible point found'/,                           ! EXIT  2
     & c( 3)   /'requested accuracy could not be achieved'/,       ! EXIT  3
     & c( 4)   /'weak QP minimizer'/,                              ! EXIT  4
     & c( 5)/'the problem appears to be infeasible'/,              ! EXIT 10
     & c( 6)   /'infeasible linear constraints'/,                  ! EXIT 11
     & c( 7)   /'infeasible linear equalities'/,                   ! EXIT 12
     & c( 8)   /'nonlinear infeasibilities minimized'/,            ! EXIT 13
     & c( 9)   /'infeasibilities minimized'/,                      ! EXIT 14
     & c(10)   /'infeasible linear constraints in QP subproblem'/, ! EXIT 15
     & c(11)/'the problem appears to be unbounded'/,               ! EXIT 20
     & c(12)   /'unbounded objective'/,                            ! EXIT 21
     & c(13)   /'constraint violation limit reached'/,             ! EXIT 22
     & c(14)/'resource limit error'/,                              ! EXIT 30
     & c(15)   /'iteration limit reached'/,                        ! EXIT 31
     & c(16)   /'major iteration limit reached'/,                  ! EXIT 32
     & c(17)   /'the superbasics limit is too small'/,             ! EXIT 33
     & c(18)/'terminated after numerical difficulties'/,           ! EXIT 40
     & c(19)   /'current point cannot be improved'/,               ! EXIT 41
     & c(20)   /'singular basis'/,                                 ! EXIT 42
     & c(21)   /'cannot satisfy the general constraints'/,         ! EXIT 43
     & c(22)   /'ill-conditioned null-space basis'/,               ! EXIT 44
     & c(23)/'error in the user-supplied functions'/,              ! EXIT 50
     & c(24)   /'incorrect objective  derivatives'/,               ! EXIT 51
     & c(25)   /'incorrect constraint derivatives'/,               ! EXIT 52
     & c(26)   /'the QP Hessian is indefinite'/,                   ! EXIT 53
     & c(27)   /'incorrect second derivatives'/,                   ! EXIT 54
     & c(28)   /'incorrect derivatives'/,                          ! EXIT 55
     & c(29)/'undefined user-supplied functions'/,                 ! EXIT 60
     & c(30)   /'undefined function at the first feasible point'/, ! EXIT 61
     & c(31)   /'undefined function at the initial point'/,        ! EXIT 62
     & c(32)   /'unable to proceed into undefined region'/,        ! EXIT 63
     & c(33)/'user requested termination'/,                        ! EXIT 70
     & c(34)   /'terminated during function evaluation'/,          ! EXIT 71
     & c(35)   /'terminated during constraint evaluation'/,        ! EXIT 72
     & c(36)   /'terminated during objective evaluation'/,         ! EXIT 73
     & c(37)   /'terminated from monitor routine'/,                ! EXIT 74
     & c(38)/'insufficient storage allocated'/,                    ! EXIT 80
     & c(39)   /'work arrays must have at least 500 elements'/,    ! EXIT 81
     & c(40)   /'not enough character storage'/,                   ! EXIT 82
     & c(41)   /'not enough integer storage'/,                     ! EXIT 83
     & c(42)   /'not enough real storage'/,                        ! EXIT 84
     & c(43)/'input arguments out of range'/,                      ! EXIT 90
     & c(44)   /'invalid input argument'/,                         ! EXIT 91
     & c(45)   /'basis file dimensions do not match this problem'/,! EXIT 92
     & c(46)   /'the QP Hessian is indefinite'/,                   ! EXIT 93
     & c(47)/'finished successfully'/,                             ! EXIT100
     & c(48)   /'SPECS file read'/,                                ! EXIT101
     & c(49)   /'Jacobian structure estimated'/,                   ! EXIT102
     & c(50)   /'MPS file read'/,                                  ! EXIT103
     & c(51)   /'memory requirements estimated'/,                  ! EXIT104
     & c(52)   /'user-supplied derivatives appear to be correct'/, ! EXIT105
     & c(53)   /'no derivatives were checked'/,                    ! EXIT106
     & c(54)   /'some SPECS keywords were not recognized'/,        ! EXIT107
     & c(55)/'errors while processing MPS data'/,                  ! EXIT110
     & c(56)   /'no MPS file specified'/,                          ! EXIT111
     & c(57)   /'problem-size estimates too small'/,               ! EXIT112
     & c(58)   /'fatal error in the MPS file'/,                    ! EXIT113
     & c(59)/'errors while estimating Jacobian structure'/,        ! EXIT120
     & c(60)   /'cannot find Jacobian structure at given point'/,  ! EXIT121
     & c(61)/'fatal errors while reading the SPECS'/,              ! EXIT130
     & c(62)   /'no SPECS file (iSpecs le 0 or iSpecs gt 99)'/,    ! EXIT131
     & c(63)   /'End-of-file while looking for a BEGIN'/,          ! EXIT132
     & c(64)   /'End-of-file while reading SPECS file'/,           ! EXIT133
     & c(65)   /'ENDRUN found before any valid SPECS'/,            ! EXIT134
     & c(66)/'system error'/,                                      ! EXIT140
     & c(67)   /'wrong no of basic variables'/,                    ! EXIT141
     & c(68)   /'error in basis package'/                          ! EXIT142
     & c(69)   /'Problem dimensions are too large'/                ! EXIT142
*     ------------------------------------------------------------------
*     Find the "major" and "minor" iExit modes

      mjr = iExit/10
      mnr = iExit - 10*mjr

      i1  = indc(mjr)
      i2  = i1 + mnr

      call s1trim( c(i1), length )
      write(string , '(1x,2a,i4,a,(a))')
     &     Solver, ' EXIT', 10*mjr, ' -- ', c(i1)(1:length)

      call s1trim( c(i2), length )
      write(string2, '(1x,2a,i4,a,(a))')
     &     Solver, ' INFO',  iExit, ' -- ', c(i2)(1:length)

      end ! subroutine snEXIT

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snWRAP
     &   ( iExit, Solver, string, string2, iw, leniw )

      implicit
     &     none
      character*(*)
     &     string, string2
      character
     &     Solver*6
      integer
     &     iExit, leniw, iw(leniw)

*     ==================================================================
*     snWRAP  it's a wrap!
*
*     18 Oct 2003: First version of snWRAP.
*     21 Jun 2004: iExit saved for use in snSolF.
*     21 Jun 2004: Current version of snWRAP.
*     ==================================================================
      integer
     &     iExit0
*     ------------------------------------------------------------------

      iw(424) = iExit   ! INFO code from all solvers

      call snEXIT( iExit, Solver, string, string2 )

      if (iExit .eq. 81) then   ! Print without using accessing iw, etc.
         call snPRNT( 15, string , iw, leniw )
         call snPRNT(  5, string2, iw, leniw )
      else
         iExit0 = iExit/10
         if (iExit0 .eq. 0) then ! Normal exit
            call s1page( 1, iw, leniw )
         else
            call s1page( 2, iw, leniw )
         end if

         call snPRNT( 3, string , iw, leniw )
         call snPRNT( 3, string2, iw, leniw )
      end if

      end ! subroutine snWRAP

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snSolF
     &   ( m, n, nb, ninf, j, jkey, jstate,
     &     hs, bl, bu, rc, xs, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     m, n, nb, ninf, j, jkey, jstate, leniw, lenrw
      integer
     &     hs(nb), iw(leniw)
      double precision
     &     bl(nb), bu(nb), rc(nb), xs(nb), rw(lenrw)

*     ==================================================================
*     snSolF sets the solution flags for the j-th variable:
*                      ' ' A  D  I  N    and   LL  UL SBS  BS  EQ  FR
*     by returning jkey=0  1  2  3  4,  jstate= 0   1   2   3   4   5
*
*     snSolF is called by SNOPT from s4soln.
*     snSolF may also be called externally (e.g. by GAMS)
*     following a normal call of SNOPT.
*     At this stage the solution will be UNSCALED!!
*     Hence, SNOPT (via m4soln) now outputs flags for the UNSCALED soln.
*
*     Input parameters m, n, nb, hs, bl, bu, rc, xs
*     are the same as for SNOPT.
*
*     j      (input ) is column j if j <= n;  otherwise row i = j - n.
*     jkey   (output) is one of 0 1 2 3 4.
*     jstate (output) is one of 0 1 2 3 4 5.
*
*     09 Mar 2004: First version of snSolF, derived from misolf.
*     18 Jun 2004: If the scaled problem was infeasible
*                  (with max inf at j = jbInf1), always flag that j
*                  as infeasible in the unscaled solution.
*     21 Jun 2004: Similarly, if the scaled problem wasn't optimal,
*                  (with max dual inf at j = jdInf1), always flag that j
*                  as nonoptimal in the unscaled solution.
*     23 Jun 2004: Suppress nonoptimal flag if iExit <= 2 (not 0).
*     ==================================================================

      logical
     &     feasbl, maximz
      integer
     &     iExit, jbInf1, jdInf1, js, minimz
      double precision
     &     b1, b2, d1, d2, dj, djtest, piNorm,
     &     tolfea, tolNLP, tolopt, tolx, xj


      minimz    = iw(199) ! (-1)(+1)    => (max)(min)
      iExit     = iw(424) ! INFO code from all solvers
      jbInf1    = iw(427) ! Largest bound infeasibility (  scaled)
      jdInf1    = iw(428) ! Largest dual  infeasibility (  scaled)

      tolNLP    = rw( 53) ! Major Optimality tolerance
      tolx      = rw( 56) ! Minor feasibility tolerance
      piNorm    = rw(422) ! Lagrange multiplier norm

      tolfea = tolx
      tolopt = tolNLP * piNorm
      feasbl = ninf   .eq. 0
      maximz = minimz .lt. 0

      js     = hs(j)
      b1     = bl(j)
      b2     = bu(j)
      xj     = xs(j)
      dj     = rc(j)
      d1     = b1 - xj
      d2     = xj - b2
      djtest = - dj

      if (feasbl) then
         if (maximz) djtest = - djtest
         jbInf1 = 0
      end if

      if (iExit .le. 2) then
         jdInf1 = 0
      end if

      ! Set keys and states.

      jkey   = 0   ! blank
      jstate = js  ! 0, 1, 2, 3

      if (js .le. 1) then            ! Nonbasic variables.
         if (b1 .eq. b2) jstate = 4
         if (- d1 .gt. tolfea  .and.     - d2 .gt. tolfea) jstate = 5
         if (jstate .eq. 1 ) djtest = - djtest
         if (jstate .ge. 4 ) djtest =   abs(djtest)
         if (                     abs(djtest) .le. tolopt) jkey = 1  ! A
         if (jstate .ne. 4     .and.  djtest  .gt. tolopt) jkey = 4  ! N

      else                           ! Basic and superbasic variables.
         if (abs(d1).le. tolfea  .or. abs(d2) .le. tolfea) jkey = 2  ! D
         if (jstate .eq. 2  .and. abs(djtest) .gt. tolopt) jkey = 4  ! N
         if (    d1 .gt. tolfea  .or.     d2  .gt. tolfea) jkey = 3  ! I
         if (     j .eq. jbInf1                          ) jkey = 3  ! I
      end if

      if (j .eq. jdInf1) jkey = 4  ! N

      end ! subroutine snSolF

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snMem
     &   ( iExit, m, n, ne, negCon,
     &     nnCon, nnJac, nnObj,
     &     mincw, miniw, minrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iExit, m, n, ne, negCon, nnCon, nnJac, nnObj, mincw, miniw,
     &     minrw, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

*     ==================================================================
*     snMem   estimates the memory requirements for snOptB,
*     using the values:
*        m    , n    , ne    negCon,
*        nnObj, nnCon, nnJac
*
*     These values are used to compute the minimum required storage:
*     mincw, miniw, minrw.
*
*     Note:
*     1. The initialization routine snInit MUST be called before snMem.
*
*     2. All default parameters must be set before calling snMem,
*        since some values affect the amount of memory required.
*
*     3. The arrays rw and iw hold  constants and work-space addresses.
*        They must have dimension at least 500.
*
*     4. This version of snMem does not allow user accessible
*        partitions of cw, iw and rw.
*
*     Exit messages:
*
*     SNMEMB EXIT  80 -- insufficient storage allocated
*     SNMEMB INFO  81 -- work arrays must have at least 500 elements
*
*     SNMEMB EXIT 100 -- finished successfully
*     SNMEMB INFO 104 -- requirements estimated
*
*     29 Mar 1998: First version.
*     15 Oct 2003: iExit added as an argument.
*     15 Oct 2003: Current version of snMem.
*     ==================================================================

      call snMemB
     &   ( iExit, m, n, ne, negCon,
     &     nnCon, nnJac, nnObj,
     &     mincw, miniw, minrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine snMem

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snMemA
     &   ( iExit, nF, n, nxname, nFname, neA, neG,
     &     mincw, miniw, minrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iExit, lencw, leniw, lenrw, mincw, miniw, minrw, n, neA, neG,
     &     nF, nFname, nxname, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

*     ==================================================================
*     snMemA   estimates the memory requirements for snOptA,
*     using the values:
*
*        nF, n, nxname, nFname, neA, neG
*
*     These values are used to compute the minimum required storage:
*     mincw, miniw, minrw.
*
*     Note:
*     1. The initialization routine snInit MUST be called before snMemA.
*
*     2. Some optional parameter settings affect the amount of memory
*        needed, so weird things may happen if some optional parameters
*        are set after the call to snMemA.
*
*     3. The arrays rw and iw hold constants and work-space addresses.
*        They must have dimension at least 500.
*
*     4. This version of snMemA does not allow user accessible
*        partitions of cw, iw and rw.
*
*     Exit messages:
*
*     SNMEMA EXIT  80 -- insufficient storage allocated
*     SNMEMA INFO  81 -- work arrays must have at least 500 elements
*
*     SNMEMA EXIT 100 -- finished successfully
*     SNMEMA INFO 104 -- requirements estimated
*
*     01 Aug 2002: First version based on snMem.
*     31 Jul 2003: snEXIT and snPRNT adopted.
*     15 Oct 2003: iExit added as an argument.
*     09 Nov 2004: Optional printing added.
*     09 Nov 2004: Current version of snMemA.
*     ==================================================================
      integer
     &     iPrint, iSumm
*     ------------------------------------------------------------------
      iPrint = iw(12)
      iSumm  = iw(13)
      call snMemA0
     &   ( iExit, iPrint, iSumm, nF, n, nxname, nFname, neA, neG,
     &     mincw, miniw, minrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine snMemA

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snMemA0
     &   ( iExit, lPrint, lSumm, nF, n, nxname, nFname, neA, neG,
     &     mincw, miniw, minrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iExit, lPrint, lSumm, lencw, leniw, lenrw, mincw, miniw,
     &     minrw, n, neA, neG, nF, nFname, nxname, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

*     ==================================================================
*     snMemA0   does the work for snMemA.
*
*     09 Nov 2004: First version
*     09 Nov 2004: Current version of  snMemA0.
*     ==================================================================
      character
     &     Solver*6, str*80, str2*80, Usercw(130)*8
      logical
     &     PrtMem
      integer
     &     inform, iPrint, iSumm, lenR, liwEst, lrwEst,
     &     llenrw, lleniw, llencw, lvlHes, maxcw, maxiw, maxrw,
     &     maxR, maxS, mQNmod, nextcw, nextiw, nextrw, nkx, nName,
     &     m, ne, negCon, nnCon, nnJac, nnObj, Useriw(130)
      double precision
     &     Userrw(130)
*     ------------------------------------------------------------------
      iPrint = iw(12)
      iSumm  = iw(13)
      iw(12) = lPrint           ! Print (Log) file
      iw(13) = lSumm            ! Specs (options) file

      Solver = 'SNMEMA'
      iExit  = 0

      if (lencw .lt. 500 .or. leniw .lt. 500 .or. lenrw .lt. 500) then
*        ---------------------------------------------------------------
*        Not enough workspace to do ANYTHING!
*        Print and exit without accessing the work arrays.
*        ---------------------------------------------------------------
         iExit = 81        ! Work arrays must have at least 500 elements
         call snWRAP( iExit, Solver, str, str2, iw, leniw )
         go to 999
      end if

*     Save the user's option choices  (weird choices get overwritten).

      call chcopy( 130, cw(51), 1, Usercw, 1 )
      call icopy ( 130, iw(51), 1, Useriw, 1 )
      call dcopy ( 130, rw(51), 1, Userrw, 1 )

*     Assign fake values for lencw, leniw, lenrw.
*     This will force s2Mem to estimate the memory requirements.

      llenrw  = 500
      lleniw  = 500
      llencw  = 500

*     An obligatory call to snInit has `undefined' all options.
*     Check the user-defined values and assign undefined values.
*     s8dflt needs various problem dimensions in iw.

*     Allocate temporary work arrays for s3size.

      nkx    = n + nF

*     Provide the user an (over) estimate of what is needed.

      ne     = neA  + neG
      m      = nF

      if (nxname .eq. 1  .and.  nFname .eq. 1) then
         nName = 1
      else
         nName = n + m
      end if

      nnCon  = m
      nnJac  = n
      nnObj  = n
      negCon = ne

      iw( 15) = n     ! copy of the number of columns
      iw( 16) = m     ! copy of the number of rows
      iw( 17) = ne    ! copy of the number of nonzeros in Jcol
      iw( 21) = nnJac ! # nonlinear Jacobian variables
      iw( 22) = nnObj ! # variables in gObj
      iw( 23) = nnCon ! # of nonlinear constraints

      call s8dflt
     &   ( m, n, nnCon, nnJac, nnObj,
     &     cw, llencw, iw, lleniw, rw, llenrw )

      nextcw   = 501
      nextiw   = 501
      nextrw   = 501

      maxcw   = lencw
      maxiw   = leniw
      maxrw   = lenrw

      maxR    = iw( 52) ! max columns of R.
      maxS    = iw( 53) ! max # of superbasics
      mQNmod  = iw( 54) ! (ge 0) max # of BFGS updates
      lvlHes  = iw( 72) ! 0,1,2  => LM, FM, Exact Hessian

      lenR    = maxR*(maxR + 1)/2  +  (maxS - maxR)
      nkx     = n + m

      call s8Map
     &   ( m, n, negCon, nkx, nnCon, nnJac, nnObj,
     &     lenR, maxR, maxS,  mQNmod, lvlHes,
     &     nextcw, nextiw, nextrw, iw, leniw )
      call s3mapA
     &   ( m, n, ne, nF, neG, negCon, nkx, nnJac, nName,
     &     nextcw, nextiw, nextrw, iw, leniw )
      call s2Bmap
     &   ( m, n, ne, maxS,
     &     nextiw, nextrw, maxiw, maxrw, liwEst, lrwEst, iw, leniw )
      PrtMem = .false.          ! Suppress messages from s2Mem
      call s2Mem
     &   ( inform, PrtMem, liwEst, lrwEst,
     &     nextcw, nextiw, nextrw,
     &     maxcw, maxiw, maxrw, llencw, lleniw, llenrw,
     &     mincw, miniw, minrw, iw )

*     mincw = mincw
      miniw = liwEst
      minrw = lrwEst

*     Restore the user's choices of options.

      call chcopy( 130, Usercw, 1, cw(51), 1 )
      call icopy ( 130, Useriw, 1, iw(51), 1 )
      call dcopy ( 130, Userrw, 1, rw(51), 1 )

*     Print the exit conditions.

      if (iExit .eq. 0) then
         iExit = 104            ! memory requirements estimated
      end if

      call snWRAP( iExit, Solver, str, str2, iw, leniw )

  999 iw(12) = iPrint
      iw(13) = iSumm

      end ! subroutine snMemA0

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snMemB
     &   ( iExit, m, n, ne, negCon,
     &     nnCon, nnJac, nnObj,
     &     mincw, miniw, minrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iExit, lencw, leniw, lenrw, m, mincw, miniw, minrw, n, ne,
     &     negCon, nnCon, nnJac, nnObj, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

*     ==================================================================
*     snMemB   estimates the memory requirements for snoptB,
*     using the values:
*        m    , n    , ne    negCon,
*        nnObj, nnCon, nnJac
*
*     These values are used to compute the minimum required storage:
*     mincw, miniw, minrw.
*
*     Note:
*     1. The initialization routine snInit MUST be called before snMemB.
*
*     2. All default parameters must be set before calling snMemB,
*        since some values affect the amount of memory required.
*
*     3. The arrays rw and iw hold  constants and work-space addresses.
*        They must have dimension at least 500.
*
*     4. This version of snMemB does not allow user-accessible
*        partitions of cw, iw and rw.
*
*     Exit messages:
*
*     SNMEMB EXIT  80 -- insufficient storage allocated
*     SNMEMB INFO  81 -- work arrays must have at least 500 elements
*
*     SNMEMB EXIT 100 -- finished successfully
*     SNMEMB INFO 104 -- requirements estimated
*
*     29 Mar 1998: First version.
*     31 Jul 2003: snEXIT and snPRNT adopted.
*     15 Oct 2003: iExit added as an argument.
*     19 Feb 2004: Current version of snMemB.
*     ==================================================================
      character
     &     Solver*6, str*80, str2*80, Usercw(130)*8
      logical
     &     PrtMem
      double precision
     &     Userrw(130)
      integer
     &     inform, lenR, liwEst, lrwEst,
     &     llenrw, lleniw, llencw, lvlHes, maxcw, maxiw, maxrw, maxR,
     &     maxS, mQNmod, nextcw, nextiw, nextrw, nkx, Useriw(130)
*     ------------------------------------------------------------------
      Solver = 'SNMEMB'
      iExit  = 0

      if (lencw .lt. 500 .or. leniw .lt. 500 .or. lenrw .lt. 500) then
*        ---------------------------------------------------------------
*        Not enough workspace to do ANYTHING!
*        Print and exit without accessing the work arrays.
*        ---------------------------------------------------------------
         iExit = 81        ! Work arrays must have at least 500 elements
         call snWRAP( iExit, Solver, str, str2, iw, leniw )
         go to 999
      end if

*     Save the user's option choices  (weird choices get overwritten).

      call chcopy( 130, cw(51), 1, Usercw, 1 )
      call icopy ( 130, iw(51), 1, Useriw, 1 )
      call dcopy ( 130, rw(51), 1, Userrw, 1 )

*     Assign fake values for lencw, leniw, lenrw.
*     This will force s2Mem to estimate the memory requirements.

      llenrw  = 500
      lleniw  = 500
      llencw  = 500

*     An obligatory call to snInit has `undefined' all options.
*     Check the user-defined values and assign undefined values.
*     s8dflt needs various problem dimensions in iw.

      iw( 15) = n     ! copy of the number of columns
      iw( 16) = m     ! copy of the number of rows
      iw( 17) = ne    ! copy of the number of nonzeros in Jcol
      iw( 21) = nnJac ! # nonlinear Jacobian variables
      iw( 22) = nnObj ! # variables in gObj
      iw( 23) = nnCon ! # of nonlinear constraints

      call s8dflt
     &   ( m, n, nnCon, nnJac, nnObj,
     &     cw, llencw, iw, lleniw, rw, llenrw )

      nextcw   = 501
      nextiw   = 501
      nextrw   = 501

      maxcw   = lencw
      maxiw   = leniw
      maxrw   = lenrw

      maxR    = iw( 52) ! max columns of R.
      maxS    = iw( 53) ! max # of superbasics
      mQNmod  = iw( 54) ! (ge 0) max # of BFGS updates
      lvlHes  = iw( 72) ! 0,1,2  => LM, FM, Exact Hessian

      lenR    = maxR*(maxR + 1)/2  +  (maxS - maxR)
      nkx     = n + m

      call s8Map
     &   ( m, n, negCon, nkx, nnCon, nnJac, nnObj,
     &     lenR, maxR, maxS,  mQNmod, lvlHes,
     &     nextcw, nextiw, nextrw, iw, leniw )
      call s2Bmap
     &   ( m, n, ne, maxS,
     &     nextiw, nextrw, maxiw, maxrw, liwEst, lrwEst, iw, leniw )
      PrtMem = .false.          ! Suppress messages from s2Mem
      call s2Mem
     &   ( inform, PrtMem, liwEst, lrwEst,
     &     nextcw, nextiw, nextrw,
     &     maxcw, maxiw, maxrw, llencw, lleniw, llenrw,
     &     mincw, miniw, minrw, iw )

*     mincw = mincw
      miniw = liwEst
      minrw = lrwEst

*     Restore the user's choices of options.

      call chcopy( 130, Usercw, 1, cw(51), 1 )
      call icopy ( 130, Useriw, 1, iw(51), 1 )
      call dcopy ( 130, Userrw, 1, rw(51), 1 )

*     Print the exit conditions.

      if (iExit .eq. 0) then
         iExit = 104            ! memory requirements estimated
      end if

      call snWRAP( iExit, Solver, str, str2, iw, leniw )

  999 return

      end ! subroutine snMemB

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snLog
     &   ( iAbort, info, Htype, KTcond, MjrPrt, minimz,
     &     n, nb, nnCon0, nS, itn, nMajor, nMinor, nSwap,
     &     condHz, iObj, sclObj, ObjAdd, fMrt, PenNrm, step,
     &     prInf, duInf, vimax, virel, hs,
     &     ne, nlocJ, locJ, indJ, Jcol,
     &     Ascale, bl, bu, fCon, yCon, x,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     KTcond(2)
      integer
     &     iAbort, info(6), iObj, Htype, lencu, lencw, leniu, leniw,
     &     lenru, lenrw, MjrPrt, minimz, n, ne, nb, nlocJ, nnCon0,
     &     nS, itn, nMajor, nMinor, nSwap, hs(nb), locJ(nlocJ),
     &     indJ(ne), iu(leniu), iw(leniw)
      double precision
     &     condHz, sclObj, ObjAdd, fMrt, PenNrm, virel, vimax, step,
     &     prInf, duInf, Ascale(nb), bl(nb), bu(nb), fCon(nnCon0),
     &     Jcol(ne), yCon(nnCon0), x(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     snLog  prints the major iteration log.
*
*     The end-of-line summary is as follows:
*
*     Position   | Possible Entries
*     -----------|-------------------
*     1 = Blank  |
*     2 = Update | n     s   -
*     3 = Modify |    M  m   -
*     4 = Htype  |    R  r   -
*     5 = Steps  |    d  l   -
*     6 = QPfea  |    i  -   -
*     7 = QPerr  |    t  T   u   w  z
*     8 = FDiff  |    c  -   -   -  -
*
*     15 Nov 1991: First version.
*     19 Jul 1997: Thread-safe version.
*     02 Dec 2000: Reordered and sparsified.
*     28 Dec 2000: Row and column permutations added.
*     31 Jul 2003: snPRNT adopted.
*     31 Jul 2003: Current version of snLog.
*     ==================================================================
      external
     &     s2ColN, s2RowN
      character
     &     cflag*1, KTflag(4)*1,
     &     MjrMsg*8, buffS*84, buffP*115, str*80, form*40
      logical
     &     Prnt1, PrntC, Summ1, nlnCon, nlnObj, prtHdg, Major0, pHead,
     &     prtx, prtl, prtf, prtj, scaled
      integer
     &     i, ir, Itns, j, k, k1, k2, lenL,
     &     lenU, LU, l, lvlScl, Mjrs, Mnrs, MjrHdP, MjrHdS, mline,
     &     nLine, nnCon, nnJac, nnObj, nfCon(4), nfObj(4), s2ColN,
     &     s2RowN
      double precision
     &     merit, infBnd, sgnObj
*     ------------------------------------------------------------------
      integer            HDiag,      HUnit
      parameter         (HDiag  = 1, HUnit  = 2)
      integer            iQNtyp,     iModfy,     iStep
      parameter         (iQNtyp = 1, iModfy = 2, iStep  = 3)
      integer            iQPfea,     iQPerr,     iFDiff
      parameter         (iQPfea = 4, iQPerr = 5, iFDiff = 6)
      integer            Scale,      UnScal
      parameter         (Scale  = 0, UnScal = 1)
      double precision   zero
      parameter         (zero   = 0.0d+0)
      parameter         (mLine  =  20)
      parameter         (MjrHdP = 224) ! >0 => Major heading for iPrint
      parameter         (MjrHdS = 226) ! >0 => Major heading for iSumm
*     ------------------------------------------------------------------
      character          flag(0:1)*1, key(-1:4)*2
      data               flag /' ', ' '/
      data               key  /'fr', 'lo', 'up', 'sb', 'bs', 'fx'/
*     ------------------------------------------------------------------
      nnCon     = iw( 23) ! # of nonlinear constraints
      nnJac     = iw( 21) ! # nonlinear Jacobian variables
      nnObj     = iw( 22) ! # variables in gObj
      lvlScl    = iw( 75) ! scale option
      lenL      = iw(173) ! size of current  L
      lenU      = iw(174) ! size of current  U

      nfCon(1)  = iw(189) ! number of calls of fCon
      nfCon(2)  = iw(190) ! number of calls of fCon
      nfCon(3)  = iw(191) ! number of calls of fCon
      nfCon(4)  = iw(192) ! number of calls of fCon

      nfObj(1)  = iw(194) ! number of calls of fObj
      nfObj(2)  = iw(195) ! number of calls of fObj
      nfObj(3)  = iw(196) ! number of calls of fObj
      nfObj(4)  = iw(197) ! number of calls of fObj

      infBnd    = rw( 70) ! definition of an infinite bound

      iAbort    = 0
      LU        = lenL + lenU

      nlnCon    = nnCon  .gt. 0
      nlnObj    = nnObj  .gt. 0
      Prnt1     = MjrPrt .ge. 1
      PrntC     = MjrPrt .ge. 100
      Summ1     = MjrPrt .ge. 1
      Major0    = nMajor .eq. 0

      if (Major0  .and.  Prnt1) then
         call snPRNT( 1, ' ', iw, leniw )
      end if

      Itns   = mod( itn   , 1000000 )
      Mnrs   = mod( nMinor, 1000000 )
      Mjrs   = mod( nMajor, 1000000 )
      MjrMsg = '_'              ! Used to detect major iteration summary
      buffP  = ' '
      buffS  = ' '

      ! Add the alphabet soup.

      if (     info(iQNtyp) .eq. 0) then
         MjrMsg(2:2) = 'n'      ! No update could be made
      else if (info(iQNtyp) .eq. 2) then
         MjrMsg(2:2) = 's'      ! Scaled BFGS
      end if

      if (     info(iModfy) .eq. 1) then
         MjrMsg(3:3) = 'M'      ! BFGS + qN Hessian mod. 1
      else if (info(iModfy) .eq. 2) then
         MjrMsg(3:3) = 'm'      ! BFGS + qN Hessian mods. 1 + 2
      end if

      if (     Htype .eq. HDiag) then
         MjrMsg(4:4) = 'R'      ! H set to a diagonal
      else if (Htype .eq. HUnit) then
         MjrMsg(4:4) = 'r'      ! H set to the identity
      end if

      if (     info(iStep)  .eq. 1) then
         MjrMsg(5:5) = 'd'      ! Violation limited via step
      else if (info(iStep)  .eq. 2) then
         MjrMsg(5:5) = 'D'      ! User-defined limit
      else if (info(iStep)  .eq. 3) then
         MjrMsg(5:5) = 'l'      ! Vars limited via step
      end if

      if (     info(iQPfea) .eq. 1) then
         MjrMsg(6:6) = 'i'      ! QP infeasible
      end if

      if (     info(iQPerr) .eq. 1) then
         MjrMsg(7:7) = 'T'      ! terminated by WSlimit
      else if (info(iQPerr) .eq. 2) then
         MjrMsg(7:7) = 't'      ! terminated by itMax
      else if (info(iQPerr) .eq. 3) then
         MjrMsg(7:7) = 'u'      ! QP unbounded
      else if (info(iQPerr) .eq. 4) then
         MjrMsg(7:7) = 'w'      ! Weak QP solutions
      else if (info(iQPerr) .eq. 5) then
         MjrMsg(7:7) = 'z'      ! superbasic limit reached
      end if

      if (     info(iFDiff) .eq. 1) then
         MjrMsg(8:8) = 'c'      ! Central differences
      end if

      ! Put ( ) around small primal and dual infeasibilities.

      do k = 1, 4
         KTflag(k) = ' '
      end do

      k    = 1
      do j = 1, 2
         if ( KTcond(j) ) then
            KTflag(k  ) = '('
            KTflag(k+1) = ')'
         end if
         k = k + 2
      end do

      sgnObj    = minimz
      merit     = sgnObj*(ObjAdd + fmrt)

      ru(1) = prInf
      ru(2) = duInf
      ru(3) = merit
      iu(1) = Mjrs
      iu(2) = Mnrs
      if ( Prnt1 ) then
*        ------------------------------------------
*        Terse line for the Print file.
*        ------------------------------------------
         prtHdg = iw(MjrHdP) .gt. 0
         nLine  = mod( Mjrs, mLine )
         pHead  = nLine  .eq. 0  .or.  prtHdg

         if (pHead) then
            cflag      = flag(min( nLine, 1 ))
            iw(MjrHdP) = 0
         end if

         if (nlnCon) then
            if (pHead) then
               call snPRNT( 11,
     &              '   Itns Major Minors    Step   nCon'
     &           // ' Feasible  Optimal  MeritFunction'
     &           // '     L+U BSwap     nS  condHz Penalty'
     &           // cflag, iw, leniw )
            end if
            write(buffP, 3000) Itns, Mjrs, Mnrs, step, nfCon(2),
     &                         KTflag(1), prInf, KTflag(2),
     &                         KTflag(3), duInf, KTflag(4),
     &                         merit, LU, nSwap, nS, condHz, PenNrm,
     &                         flag(0), MjrMsg
            if (condHz .eq. zero) buffP( 90: 97) = ' '
            if (PenNrm .eq. zero) buffP( 98:105) = ' '

         else if (nlnObj) then
            if (pHead) then
               call snPRNT( 11,
     &              '   Itns Major Minors    Step   nObj'
     &           // ' Feasible  Optimal      Objective'
     &           // '     L+U BSwap     nS  condHz'
     &           // cflag, iw, leniw )
            end if
            write(buffP, 3100) Itns, Mjrs, Mnrs, step, nfObj(2),
     &                         KTflag(3), duInf, KTflag(4),
     &                         merit, LU, nSwap, nS, condHz,
     &                         flag(0), MjrMsg
            if (condHz .eq. zero) buffP( 90: 97) = ' '
         else
            if (pHead) then
               if (nS .gt. 0) then
                  call snPRNT( 11,
     &                 '   Itns Major Minors    Step       '
     &              // ' Feasible  Optimal    LPobjective'
     &              // '     L+U BSwap     nS'
     &              // cflag, iw, leniw )
               else
                  call snPRNT( 11,
     &                 '   Itns Major Minors    Step       '
     &              // ' Feasible  Optimal    LPobjective'
     &              // '     L+U             '
     &              // cflag, iw, leniw )
               end if
            end if
            write(buffP, 3200) Itns, Mjrs, Mnrs, step, nfObj(2),
     &                         KTflag(3), duInf, KTflag(4),
     &                         merit, LU, nSwap, nS,
     &                         flag(0), MjrMsg
            buffP( 29: 35) = ' '
         end if

         if (step   .eq. zero) buffP( 21: 28) = ' '
         if (nSwap  .eq.    0) buffP( 77: 82) = ' '
         if (nS     .eq.    0) buffP( 83: 89) = ' '

         call snPRNT( 1, buffP, iw, leniw )
      end if

      if ( Summ1 ) then
*        --------------------------------------------
*        Terse line for the Summary file.
*        --------------------------------------------
         MjrMsg(1:1) = ' '
         prtHdg = iw(MjrHdS)    .gt. 0
         pHead  = mod(Mjrs, 10) .eq. 0  .or.  prtHdg
     &                                  .or.  Major0

         if (pHead) then
            iw(MjrHdS) = 0
         end if

         if (nlnCon) then
            if (pHead) then
               call snPRNT( 12,
     &              ' Major Minors     Step   nCon'
     &           // ' Feasible  Optimal  MeritFunction    nS'
     &           // ' Penalty', iw, leniw )
            end if
            write(buffS, 5000) Mjrs, Mnrs, step, nfCon(2),
     &                         KTflag(1), prInf, KTflag(2),
     &                         KTflag(3), duInf, KTflag(4),
     &                         merit, nS, PenNrm, MjrMsg
            if (PenNrm .eq. zero) buffS(69:76) = ' '

         else if (nlnObj) then
            if (pHead) then
               call snPRNT( 12,
     &              ' Major Minors     Step   nObj'
     &           // ' Feasible  Optimal      Objective    nS',
     &              iw, leniw )
            end if
            write(buffS, 5100) Mjrs, Mnrs, step, nfObj(2),
     &                         KTflag(3), duInf, KTflag(4),
     &                         merit, nS,         MjrMsg
         else
            if (pHead) then
               if (nS .gt. 0) then
                  call snPRNT( 12,
     &                 ' Major Minors     Step       '
     &              // ' Feasible  Optimal    LPobjective    nS',
     &                 iw, leniw )
               else    ! Once zero, nS remains zero.
                  call snPRNT( 12,
     &                 ' Major Minors     Step       '
     &              // ' Feasible  Optimal    LPobjective',
     &                 iw, leniw )
               end if
            end if
            write(buffS, 5100) Mjrs, Mnrs, step, nfObj(2),
     &                         KTflag(3), duInf, KTflag(4),
     &                         merit, nS,         MjrMsg
            buffS( 29: 35) = ' '   ! Zap nObj for LPs.
         end if

         if (step   .eq. zero) buffS(14:22) = ' '
         if (nS     .eq.    0) buffS(63:68) = ' '

         call snPRNT( 2, buffS, iw, leniw )
      end if

      if (PrntC  .and.  nnCon .gt. 0) then
*        ---------------------------------------------------------------
*        Output heading for detailed log.
*        ---------------------------------------------------------------
         call s1page( 0, iw, leniw )
         if (Major0) call snPRNT( 1, ' ', iw, leniw )

*        Unscale everything if necessary.

         scaled = lvlScl .ge. 2
         if ( scaled ) then
            call s2scla
     &         ( UnScal, nnCon, n, nb, iObj, infBnd, sclObj,
     &           ne, nlocJ, locJ, indJ, Jcol,
     &           Ascale, bl, bu, yCon, x )
            call ddscl
     &         ( nnCon, Ascale(n+1), 1, fCon, 1 )
         end if

         i      = 0             ! Keeps ftnchek happy
         l      = MjrPrt/100
         prtx   = mod( l,10 ) .gt. 0
         l      = l/10
         prtl   = mod( l,10 ) .gt. 0  .and.  Mjrs .gt. 0
         l      = l/10
         prtf   = mod( l,10 ) .gt. 0
         l      = l/10
         prtj   = mod( l,10 ) .gt. 0
         form   = '(1p, i9, e13.5)'

         if ( prtx ) then
            call snPRNT( 11, ' Jacobian variables', iw, leniw )
            call snPRNT(  1, ' ------------------', iw, leniw )
            do j = 1, nnJac
               write(str, form) s2ColN( j, leniw, iw ), x(j)
               call snPRNT( 1, str, iw, leniw )
            end do
         end if

         if ( prtl ) then
            call snPRNT( 11, ' Multiplier estimates', iw, leniw )
            call snPRNT(  1, ' --------------------', iw, leniw )
            do i = 1, nnCon
               write(str, form) s2RowN( i, leniw, iw ), yCon(i)
               call snPRNT( 1, str, iw, leniw )
            end do
         end if

         if ( prtf ) then
            call snPRNT( 11, ' Constraint functions', iw, leniw )
            call snPRNT(  1, ' --------------------', iw, leniw )
            do i = 1, nnCon
               write(str, form) s2RowN( i, leniw, iw ), fCon(i)
               call snPRNT( 1, str, iw, leniw )
            end do
            write(str, 7600) vimax, virel
            call snPRNT( 11, str, iw, leniw )
         end if

         if ( prtj ) then
            call snPRNT( 11, ' x  and  Jacobian', iw, leniw )
            call snPRNT(  1, ' ----------------', iw, leniw )
            do 160 j  = 1, nnJac
               l  = hs(j)
               write(str, 7410) s2ColN( j, leniw, iw ), x(j), key(l)
               call snPRNT( 11, str, iw, leniw )

               k1 = locJ(j)
               k2 = locJ(j+1) - 1
               do k  = k1, k2
                  ir = indJ(k)
                  if (ir .gt. nnCon) go to 160
                  write(str, 7420) s2RowN( indJ(k), leniw, iw ), Jcol(k)
                  call snPRNT( 1, str, iw, leniw )
               end do
  160       continue
         end if

*        Scale again if necessary.

         if (scaled) then
            call s2scla
     &         ( Scale, nnCon, n, nb, iObj, infBnd, sclObj,
     &           ne, nlocJ, locJ, indJ, Jcol,
     &           Ascale, bl, bu, yCon, x )
            call dddiv ( nnCon, Ascale(n+1), 1, fCon, 1 )
         end if
      end if
      return

*     Major log,  Print file.

 3000 format(i7, i6, i7, 1p, e8.1, i7, 1x, 2(a,e7.1,a), e14.7,   ! nlnCon
     &       i8, i6, i7, e8.1, e8.1, a1, a)
 3100 format(i7, i6, i7, 1p, e8.1, i7, 10x, 1(a,e7.1,a), e14.7,  ! nlnObj
     &       i8, i6,       i7, e8.1, a1, a)
 3200 format(i7, i6, i7, 1p, e8.1, i7, 10x, 1(a,e7.1,a), e14.7,  ! LP with nS>0
     &       i8, i6,       i7,       a1, a)

*     Major log,  Summary file.

 5000 format(i6, i7, 1p, e9.1, i7,  1x, 2(a,e7.1,a), e14.7, i6, e8.1, a)
 5100 format(i6, i7, 1p, e9.1, i7, 10x, 1(a,e7.1,a), e14.7, i6, a)
 7410 format(' x(', i6, ')', 1p, e13.5, 1x, a2)
 7420 format('   ', i6, ' ', 1p, e13.5)
 7600 format(' Maximum constraint violation    =', 1p, e12.4,
     &       4x, ' ( =', e11.4, ' normalized)' )

      end ! subroutine snLog

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snLog2
     &   ( Prob, ProbTag,
     &     Elastc, gotR, jstFea, feasbl,
     &     m, mBS, nnH, nS, jSq, jBr, jSr,
     &     linesP, linesS, itn, itQP, kPrc, lvlInf,
     &     pivot, step, nInf, sInf, wtInf,
     &     ObjPrt, condHz, djqPrt, rgNorm, kBS, xBS,
     &     iw, leniw )

      implicit
     &     none
      character
     &     ProbTag*20
      logical
     &     Elastc, gotR, jstFea, feasbl
      integer
     &     Prob, m, mBS, nnH, nS, jSq, jBr, jSr, itn, itQP, kPrc,
     &     linesP, linesS, lvlInf, nInf, kBS(mBS), leniw, iw(leniw)
      double precision
     &     condHz, djqPrt, ObjPrt, pivot, rgNorm, step, sInf, wtInf,
     &     xBS(mBS)

*     ==================================================================
*     snLog2  prints the minor iteration log.
*
*     mBS = m + maxS  .ge.  m + nS.
*
*     MnrHdP is  1  if a new heading is required for some reason other
*     than frequency (e.g., after a basis factorization).
*
*     The output consists of a number of ``sections'' of one line
*     summaries, with each section preceded by a header message.
*     linesP and linesS count the number of lines remaining to be
*     printed in each section of the print and summary files
*     respectively.   They too may force a new heading.
*
*     01 Dec 1991: First version based on Minos routine m5log.
*     17 Nov 2000: First version for SQP minor itns.
*     02 Dec 2000: Reordered and sparsified.
*     28 Dec 2000: Row and column permutations added.
*     01 Aug 2003: cgItn added to Print log.
*     04 Jul 2005: Current version of snLog2.
*     ==================================================================
      external
     &     s1intmx, s2VarN
      character
     &     buffP*138, buffS*75, str*80
      logical
     &     prtLog, prtSum, prtHdP, prtHdS, newSet, pHead
      integer
     &     cgItn, Itns, jSqN, jBrN, jSrN, k, lenL, lenU, lprDbg, lvlSys,
     &     MjrHdP, MjrHdS, MnrHdP, MnrHdS, MnrPrt, Mnrs, ncp,
     &     PrintP, PrintS, QPmode, s1intmx, s2VarN, width
      double precision
     &     mxwdth, rmxint
*     ------------------------------------------------------------------
      integer            CG
      parameter         (CG     = 1)
      integer            FP,         FPE,        FPS
      parameter         (FP     = 0, FPE    = 3, FPS    = 4)
      integer            YES
      parameter         (YES    = 1)
      integer            mLineP,      mLineS
      parameter         (mLineP = 40, mLineS = 10)
      double precision   zero
      parameter         (zero   = 0.0d+0)
      parameter         (MnrHdP = 223) ! >0 => Minor heading for Log
      parameter         (MjrHdP = 224) ! >0 => Major heading for Log
      parameter         (MnrHdS = 225) ! >0 => Minor heading for Summary
      parameter         (MjrHdS = 226) ! >0 => Major heading for Summary
*     ------------------------------------------------------------------
      lvlSys = iw( 71) ! > 0   => print system info
      lprDbg = iw( 85) ! > 0   => private debug print
      MnrPrt = iw( 93) ! Minor print level
      ncp    = iw(176) ! no. of LU compressions
      lenL   = iw(173) ! size of current  L
      lenU   = iw(174) ! size of current  U
      QPmode = iw(208) ! Current QP solver   (based on QPslvr)
      PrintP = iw(218) ! (on/off) current log     file output status
      PrintS = iw(219) ! (on/off) current summary file output status
      cgItn  = iw(387) ! symmlq itns for the last QP minor itn

      prtLog = PrintP     .eq. YES
      prtSum = PrintS     .eq. YES

      mxwdth = 1.0d+7           ! Integers printed i7
      rmxint = s1intmx()        ! Largest integer without overflow
      if (mxwdth .gt. rmxint) mxwdth = rmxint
      width  = mxwdth

      Itns   = mod( itn , width )
      Mnrs   = mod( itQP, width )

      buffP  = ' '
      buffS  = ' '
                                ! If  newly feasible, print something.
      if (jstFea  .and.  MnrPrt .ge. 10) then

         if (.not. Elastc) then
                                ! Constraints feasible in Normal mode.
                                ! Print a message.
                                ! ProbTag is one of the following:
                                ! ProbTag = 'QP problem'
                                ! ProbTag = 'LP problem'
                                ! ProbTag = 'QP subproblem'
                                ! ProbTag = 'Linear constraints'
            if (Prob .ne. FPS  .and. Prob .ne. FP
     &                         .and. Prob .ne. FPE) then
               write(str, 8010) itn, ProbTag
               call snPRNT( 23, str, iw, leniw )
            end if
         else
                                ! Elastic mode
                                ! Elastic Phase 1 has completed.
            if (lvlInf .eq. 2) then
                                ! Infinite weight on sumInf.
                                ! Minimize the infeasible elastics.
               write(str, 8030) itn
               call snPRNT( 23, str, iw, leniw )

            else if (lvlInf .eq. 1) then
                                ! Finite nonzero weight on sumInf
                                ! Minimize a weighted objective.
               write(str, 8040) itn
               call snPRNT( 23, str, iw, leniw )
            end if
         end if

         if (lvlSys .gt. 0) then
            iw(MnrHdP) = 1      ! Print the header to the print   file
            iw(MnrHdS) = 1      ! Print the header to the summary file
         end if
      end if

      prtHdP = iw(MnrHdP) .gt. 0
      prtHdS = iw(MnrHdS) .gt. 0

      if (prtLog) then
*        --------------------------------------
*        Terse line for the Print file.
*        --------------------------------------
         newSet = linesP .eq. 0
         pHead  = prtHdP  .or.  newSet

         if ( pHead ) then
            iw(MnrHdP) = 0
            linesP     = mlineP
         end if

         iw(MjrHdP) = 1
         linesP     = linesP - 1

         jSqN  = s2VarN( jSq , leniw, iw )
         jSrN  = s2VarN( jSr , leniw, iw )
         jBrN  = s2VarN( jBr , leniw, iw )

         if (nnH .gt. 0) then
            if ( pHead ) then
               buffP = '    Itn       QPmult  QPstep'
     &              // '   nInf   SumInf   rgNorm    QPobjective'
     &              // '   +SBS   -SBS    -BS    Pivot'
     &              // '     L+U ncp    nS  condHz'
               if (Elastc        ) buffP( 56: 68) = 'Elastic QPobj'
               if (QPmode .eq. CG) buffP(126:131) = 'cgItns'
               call snPRNT( 11, buffP, iw, leniw )
            end if
            write(buffP, 3000) Itns, djqPrt, step,
     &           nInf, sInf, rgnorm, ObjPrt,
     &           jSqN, jSrN, jBrN, pivot,
     &           lenL+lenU, ncp, nS, condHz, cgItn
         else  ! nnH == 0
            if ( pHead ) then
               buffP = '    Itn       LPmult  LPstep'
     &              // '   nInf   SumInf             LPobjective'
     &              // '   +SBS   -SBS    -BS    Pivot'
     &              // '     L+U ncp'
               if (Elastc        ) buffP( 56: 68) = 'Elastic LPobj'
               if (nS     .gt.  0) buffP(115:116) = 'nS'
               call snPRNT( 11, buffP, iw, leniw )
            end if

            write(buffP, 3000) Itns, djqPrt, step,
     &           nInf, sInf, rgnorm, ObjPrt,
     &           jSqN, jSrN, jBrN, pivot,
     &           lenL+lenU, ncp, nS
         end if

         if (djqPrt .eq. zero) buffP( 13: 20) = ' '
         if (step   .eq. zero) buffP( 21: 28) = ' '
         if (nInf   .eq.    0) buffP( 29: 44) = ' ' ! nInf and sInf
         if (rgnorm .eq. zero) buffP( 45: 53) = ' '
         if (    .not. feasbl) buffP( 54: 68) = ' '
         if (jSq    .eq.    0) buffP( 69: 75) = ' '
         if (jSr    .eq.    0) buffP( 76: 82) = ' '
         if (jBr    .eq.    0) buffP( 83: 89) = ' '
         if (pivot  .eq. zero) buffP( 90: 98) = ' '
         if (ncp    .eq.    0) buffP(107:110) = ' '
         if (nS     .eq.    0) buffP(111:116) = ' '
         if (condHz .eq. zero) buffP(117:124) = ' ' ! condHz
         if (cgItn  .eq.    0) buffP(125:131) = ' '
         call snPRNT( 1, buffP, iw, leniw )
      end if

      if (prtSum) then
*        --------------------------------
*        Terse line for the Summary file.
*        --------------------------------
         newSet = linesS .eq. 0
         pHead  = prtHdS  .or.  newSet

         if ( pHead ) then
            iw(MnrHdS) = 0
            linesS = mlineS
         end if

         iw(MjrHdS) = 1
         linesS     = linesS - 1

         if (nnH .gt. 0) then
            if ( pHead ) then
               buffS = '        Minor   QPmult   nInf'
     &              // '   SumInf   rgNorm    QPobjective'
               if (Elastc        ) buffS(50:62) = 'Elastic QPobj'
               if (nS     .gt.  0) buffS(67:68) = 'nS'
               if (QPmode .eq. CG) buffS(70:75) = 'cgItns'
               call snPRNT( 12, buffS, iw, leniw )
            end if

            write(buffS, 5000) Mnrs, djqPrt, nInf, sInf,
     &           rgNorm, ObjPrt, nS, cgItn

         else  ! nnH == 0
            if ( pHead ) then
               buffS = '        Minor   LPmult   nInf'
     &              // '   SumInf   rgNorm    LPobjective'
               if (Elastc        ) buffS(50:62) = 'Elastic LPobj'
               if (nS     .gt.  0) buffS(67:68) = 'nS'
               call snPRNT( 12, buffS, iw, leniw )
            end if

            write(buffS, 5000) Mnrs, djqPrt, nInf, sInf,
     &           rgNorm, ObjPrt, nS
         end if

         if (djqPrt .eq. zero) buffS(14:22) = ' '
         if (nInf   .eq.    0) buffS(23:29) = ' '
         if (sInf   .eq. zero) buffS(30:38) = ' '
         if (rgNorm .eq. zero) buffS(41:47) = ' '
         if (    .not. feasbl) buffS(48:62) = ' '
         if (nS     .eq.    0) buffS(63:68) = ' '
         if (cgItn  .eq.    0) buffS(69:75) = ' '

         call snPRNT( 2, buffS, iw, leniw )
      end if

*     ------------------------------------------------------------------
*     Debug output.
*     ------------------------------------------------------------------
      if (lprDbg .eq. 100) then
         call snPRNT( 11, ' BS values...', iw, leniw )
         do k = 1, m
            write(buffP, 6000) s2VarN( kBS(k), leniw,iw ), xBS(k)
            call snPRNT( 1, buffP, iw, leniw )
         end do

         call snPRNT( 11, ' SB values...', iw, leniw )
         do k = m+1, m+nS
            write(buffP, 6000) s2VarN( kBS(k), leniw,iw ), xBS(k)
            call snPRNT( 1, buffP, iw, leniw )
         end do
      end if

      return

                                ! Minor log,  Print   file.
 3000 format(1p, i7, 4x, e9.1, e8.1, i7, 2e9.1, e15.7,
     &          3i7, e9.1, i8, i4, i6, e8.1, i7 )
                                ! Minor log,  Summary file.
 5000 format(1p, i13, e9.1, i7, e9.1, e9.1, e15.7, i6, i7)
 6000 format(i7, g17.8)
 8010 format(  ' Itn', i7, ': Feasible ', a)
 8030 format(  ' Itn', i7, ': Elastic Phase 2 -- minimizing',
     &                     ' elastic variables')
 8040 format(  ' Itn', i7, ': Elastic Phase 2 -- minimizing',
     &                     ' obj + weighted elastics')

      end ! subroutine snLog2

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snSTOP
     &   ( iAbort, info, Htype, KTcond, MjrPrt, minimz,
     &     n, nb, nnCon0, nS, itn, nMajor, nMinor, nSwap,
     &     condHz, iObj, sclObj, ObjAdd, fMrt, PenNrm, step,
     &     prInf, duInf, vimax, virel, hs,
     &     ne, nlocJ, locJ, indJ, Jcol,
     &     Ascale, bl, bu, fCon, yCon, x,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     KTcond(2)
      integer
     &     iAbort, info(6), iObj, Htype, lencu, lencw, leniu, leniw,
     &     lenru, lenrw, MjrPrt, minimz, n, ne, nb, nlocJ, nnCon0,
     &     nS, itn, nMajor, nMinor, nSwap, hs(nb), locJ(nlocJ),
     &     indJ(ne), iu(leniu), iw(leniw)
      double precision
     &     condHz, sclObj, ObjAdd, fMrt, PenNrm, virel, vimax, step,
     &     prInf, duInf, Ascale(nb), bl(nb), bu(nb), fCon(nnCon0),
     &     Jcol(ne), yCon(nnCon0), x(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     snSTOP is called every major iteration.
*     If iAbort > 0 on exit, the run is terminated.
*     By specifying a custom version of snSTOP, the user can arrange for
*     snopt to be terminated at any given major iteration.
*
*     14 Oct 2004: First version of   snSTOP.
*     14 Oct 2004: Current version of snSTOP.
*     ==================================================================

      iAbort    = 0
*     Relax

      end ! subroutine snSTOP

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snSet
     &   ( buffer, iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     iPrint, iSumm, Errors, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

*     ==================================================================
*     snSet  decodes the option contained in  buffer.
*
*     The buffer is output to file iPrint, minus trailing blanks.
*     Error messages are output to files iPrint and iSumm.
*     Buffer is echoed to iPrint but normally not to iSumm.
*     It is echoed to iSumm before any error msg.
*
*     On entry,
*     iPrint is the print   file.  no output occurs if iPrint .le 0.
*     iSumm  is the Summary file.  no output occurs if iSumm  .le 0.
*     Errors is the number of errors so far.
*
*     On exit,
*     Errors is the number of errors so far.
*
*     27 Nov 1991: first version of snSet.
*     03 Nov 2000: current version.
*     ==================================================================
      integer
     &     ivalue
      double precision
     &     rvalue
      character
     &     cvalue*8, key*16
*     ------------------------------------------------------------------
      call s3opt
     &   ( .true., buffer, key, cvalue, ivalue, rvalue,
     &     iPrint, iSumm, Errors, cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine snSet

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snSeti
     &   ( buffer, ivalue, iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     ivalue, iPrint, iSumm, Errors, lencw, leniw, lenrw,
     &     iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

*     ==================================================================
*     snSeti decodes the option contained in  buffer // ivalue.
*     The parameters other than ivalue are as in snSet.
*
*     27 Nov 1991: first version of snSeti.
*     03 Nov 2000: current version.
*     ==================================================================
      integer
     &     ivalxx, lenbuf
      double precision
     &     rvalue
      character
     &     cvalue*8, key*16, buff72*72
*     ------------------------------------------------------------------
      write(key, '(i16)') ivalue
      lenbuf = len(buffer)
      buff72 = buffer
      buff72(lenbuf+1:lenbuf+16) = key
      ivalxx = ivalue
      call s3opt
     &   ( .true., buff72, key, cvalue, ivalxx, rvalue,
     &     iPrint, iSumm, Errors, cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine snSeti

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snSetr
     &   ( buffer, rvalue, iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     iPrint, iSumm, Errors, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rvalue, rw(lenrw)
      character
     &     cw(lencw)*8

*     ==================================================================
*     snSetr decodes the option contained in  buffer // rvalue.
*     The parameters other than rvalue are as in snSet.
*
*     27 Nov 1991: first version of snSetr.
*     03 Nov 2000: current version.
*     ==================================================================
      integer
     &     ivalue, lenbuf
      character
     &     cvalue*8, key*16, buff72*72
      double precision
     &     rvalxx
*     ------------------------------------------------------------------
      write(key, '(1p, e16.8)') rvalue
      lenbuf = len(buffer)
      buff72 = buffer
      buff72(lenbuf+1:lenbuf+16) = key
      rvalxx = rvalue
      call s3opt
     &   ( .true., buff72, key, cvalue, ivalue, rvalxx,
     &     iPrint, iSumm, Errors, cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine snSetr

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      integer function snGet
     &   ( buffer, Errors, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     Errors, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

*     ==================================================================
*     snGet  decodes the option contained in  buffer
*     and returns 1 if the option has previously been set, else 0.
*     For example,
*     i = snGet ( 'Maximize', Errors, cw, lencw, iw, leniw, rw, lenrw )
*
*     01 Aug 2003: First version of snGet.  Needed because
*                  snGetc, snGeti, snGetr were not well defined
*                  for strings that had no numerical value.
*     01 Aug 2003: Current version of snGet.
*     ==================================================================
      integer
     &     ivalue
      double precision
     &     rvalue
      character
     &     cvalue*8, key*16
*     ------------------------------------------------------------------
      call s3opt
     &   ( .false., buffer, key, cvalue, ivalue, rvalue,
     &     0, 0, Errors, cw, lencw, iw, leniw, rw, lenrw )

      snGet  = ivalue

      end ! integer function snGet

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snGetc
     &   ( buffer, cvalue, Errors, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     Errors, lencw, leniw, lenrw, iw(leniw)
      character
     &     cvalue*8, cw(lencw)*8
      double precision
     &     rw(lenrw)

*     ==================================================================
*     snGetc gets the value of the option contained in  buffer.
*     The parameters other than cvalue are as in snSet.
*
*     17 May 1998: first version of snGetc.
*     03 Nov 2000: current version.
*     ==================================================================
      integer
     &     ivalue
      double precision
     &     rvalue
      character
     &     key*16
*     ------------------------------------------------------------------
      call s3opt
     &   ( .false., buffer, key, cvalue, ivalue, rvalue,
     &     0, 0, Errors, cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine snGetc

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snGeti
     &   ( buffer, ivalue, Errors, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     ivalue, Errors, lencw, leniw, lenrw, iw(leniw)
      character
     &     cw(lencw)*8
      double precision
     &     rw(lenrw)

*     ==================================================================
*     snGeti gets the value of the option contained in  buffer.
*     The parameters other than ivalue are as in snSet.
*
*     17 May 1998: first version of snGeti.
*     03 Nov 2000: current version.
*     ==================================================================
      double precision
     &     rvalue
      character
     &     key*16, cvalue*8
*     ------------------------------------------------------------------
      call s3opt
     &   ( .false., buffer, key, cvalue, ivalue, rvalue,
     &     0, 0, Errors, cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine snGeti

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snGetr
     &   ( buffer, rvalue, Errors, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     Errors, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rvalue, rw(lenrw)
      character
     &     cw(lencw)*8

*     ==================================================================
*     snGetr gets the value of the option contained in  buffer.
*     The parameters other than rvalue are as in snSet.
*
*     17 May 1998: first version of snGetr.
*     03 Nov 2000: current version.
*     ==================================================================
      integer
     &     ivalue
      character
     &     key*16, cvalue*8
*     ------------------------------------------------------------------
      call s3opt
     &   ( .false., buffer, key, cvalue, ivalue, rvalue,
     &     0, 0, Errors, cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine snGetr

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snRetH
     &   ( Errors, lenH, H, nnH,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Errors, lencw, lenH, leniw, lenrw, nnH,
     &     iw(leniw)
      double precision
     &     H(lenH), rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     snRetH  retrieves the SNOPT approximate Hessian from memory.
!
!     snRetH must be called immediately after a call to npOpt, snoptB
!     or snoptC  with the option "Hessian full memory" selected.
!
!     On entry:
!        lenH      (ge 1) is the length the array H, i.e., H(lenH).
!
!     On exit
!        Errors    contains the number of errors.
!                  if Errors gt 0, then error messages are printed on
!                  the standard output
!        H(lenH)   contains the upper triangular part of the approximate
!                  Hessian, stored by columns.
!        nnH       contains the number of columns of H.
!
!
!     03 Sep 2006: First version of snRetH
!     ==================================================================
      character
     &     Solver*6, str*80, str2*80
      integer
     &     iExit, lenU,
     &     lvlHes, lU, ly, ly1
!     ------------------------------------------------------------------
      integer            FM
      parameter         (FM     = 1)
!     ------------------------------------------------------------------
      Solver = 'SNRETH'
      Errors = 0

      if (lencw .lt. 500 .or. leniw .lt. 500 .or. lenrw .lt. 500) then
*        ---------------------------------------------------------------
*        Not enough workspace to do ANYTHING!
*        Print and exit without accessing the work arrays.
*        ---------------------------------------------------------------
         iExit  = 81       ! Work arrays must have at least 500 elements
         Errors = Errors + 1
         call snWRAP( iExit, Solver, str, str2, iw, leniw )
         go to 999
      end if

!     ------------------------------------------------------------------
!     Retrieve the approximate Hessian.
!     ------------------------------------------------------------------
      lvlHes    = iw( 72)       ! LM, FM or Exact Hessian
      nnH       = iw( 24)       ! max( nnObj, nnJac )
      ly        = iw(311)       ! y (nb)      =  real work vector
      ly1       = iw(312)       ! y1(nb)      =  real work vector

      lU        = iw(391)       ! U(lenU), BFGS Hessian H = U'U
      lenU      = iw(392)       !

!     Check pointers, lengths, etc,  retrieved from memory.

      call s4chkP ( Errors, 'lvlHes', lvlHes, iw, leniw )
      call s4chkP ( Errors, '    ly',     ly, iw, leniw )
      call s4chkP ( Errors, '   ly1',    ly1, iw, leniw )
      call s4chkP ( Errors, '    lU',     lU, iw, leniw )

      if (lenH .lt. lenU) then
         Errors = Errors + 1
         write(str, 9990) lenU, lenH
         call snPRNT( 5, str, iw, leniw )
      end if

      if (lvlHes .ne. FM) then
         Errors = Errors + 1
         write(str, 9991)
         call snPRNT( 5, str, iw, leniw )
      end if

      if (Errors .eq. 0) then
         call s8getH
     &      ( nnH, lenH, rw(lU), H, rw(ly), rw(ly1) )
      end if

  999 return

 9990 format(' XXX  lenH too small: needs to be at least', i6 )
 9991 format(' XXX  Full-memory Hessian not requested' )

      end ! subroutine snRetH

