*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  npOpt.f --- the NPSOL wrapper for SNOPT.
*
*     npOpt   npKerN
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npOpt
     &   ( n, nclin, ncnln, ldA, ldcJ, ldH,
     &     A, bl, bu, funcon, funobj,
     &     INFO, majIts, iState,
     &     cCon, cJac, cMul, Objf, grad, Hess, x,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     funcon, funobj
      integer
     &     n, nclin, ncnln, ldA, ldcJ, ldH, INFO, majIts,
     &     leniw, lenrw, iState(n+nclin+ncnln), iw(leniw)
      double precision
     &     Objf, A(ldA,*), bl(n+nclin+ncnln), bu(n+nclin+ncnln),
     &     cCon(*), cJac(ldcJ,*), cMul(n+nclin+ncnln), grad(n),
     &     Hess(ldH,*), x(n), rw(lenrw)

*     ==================================================================
*     npOpt  solves the nonlinear programming problem
*
*            minimize                   f(x)
*
*                                    (      x  )
*            subject to    bl  .le.  (    A*x  )  .le.  bu
*                                    ( cCon(x) )
*
*     where  f(x)  is a smooth scalar function,  A  is a constant matrix
*     and  cCon(x)  is a vector of smooth nonlinear functions.
*     The feasible region is defined by a mixture of linear and
*     nonlinear equality or inequality constraints on  x.
*
*     The calling sequence of NPOPT and the user-defined functions
*     funcon and funobj are identical to those of the dense code NPSOL
*     (see the User's Guide for NPSOL (Version 4.0): a Fortran Package
*     for Nonlinear Programming, Systems Optimization Laboratory Report
*     SOL 86-2, Department of Operations Research, Stanford University,
*     1986.)
*
*     The dimensions of the problem are...
*
*     n        the number of variables (dimension of  x),
*
*     nclin    the number of linear constraints (rows of the matrix  A),
*
*     ncnln    the number of nonlinear constraints (dimension of  c(x)),
*
*     NPOPT  is maintained by Philip E. Gill,
*     Dept of Mathematics, University of California, San Diego.
*
*     LUSOL is maintained by Michael A. Saunders,
*     Systems Optimization Laboratory,
*     Dept of Management Science & Engineering, Stanford University.
*
*     22 Mar 1997: First   version of npOpt.
*     31 Jul 2003: snEXIT and snPRNT adopted.
*     15 Oct 2004: snSTOP adopted.
*     13 Oct 2003: Current version of npOpt.
*     ==================================================================
      external
     &     snLog, snLog2, snSTOP
*     ------------------------------------------------------------------
      call npKerN
     &   ( n, nclin, ncnln, ldA, ldcJ, ldH,
     &     A, bl, bu, funcon, funobj, snLog, snLog2, snSTOP,
     &     INFO, majIts, iState,
     &     cCon, cJac, cMul, Objf, grad, Hess, x,
     &     iw, leniw, rw, lenrw )

      end ! subroutine npOpt

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npKerN
     &   ( n, nclin, ncnln, ldA, ldcJ, ldH,
     &     A, bl, bu, funcon, funobj, snLog, snLog2, snSTOP,
     &     INFO, majIts, iState,
     &     cCon, cJac, cMul, Objf, grad, Hess, x,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     funcon, funobj, snLog, snLog2, snSTOP
      integer
     &     n, nclin, ncnln, ldA, ldcJ, ldH, INFO, majIts,
     &     leniw, lenrw, iState(n+nclin+ncnln), iw(leniw)
      double precision
     &     Objf, A(ldA,*), bl(n+nclin+ncnln), bu(n+nclin+ncnln),
     &     cCon(*), cJac(ldcJ,*), cMul(n+nclin+ncnln), grad(n),
     &     Hess(ldH,*), x(n), rw(lenrw)

*     ==================================================================
*     npKerN does the work for npOpt. (Kernel for npOpt)
*
*     Developers can call this version with customized versions of
*     snLog, snLog2  and  snSTOP.
*
*     17 Oct 2004: First version of npKerN.
*     05 Dec 2004: Current version of npKerN.
*     ==================================================================
      character
     &     cStart*8, Solver*6, str*80, str2*80
      external
     &     dnrm1s, s0fgN
      integer
     &     Useriw(130)
      double precision
     &     Userrw(130)
      logical
     &     FPonly, gotR, PrtMem
      integer
     &     Errors, Htype, i, inform, iObj, lbl, lbu, lencw,
     &     lenR, lfCon, lgCon, lgObj, lGsav, lhs, lindJ, liwEst,
     &     lrwEst, lJcol, llocG, llocJ, lNames, lkx, lpi, lprSol,
     &     lprSav, lrc, lvlSrt, lvlHes, lx, lx0, m, maxcw, maxiw, maxR,
     &     maxrw, maxS, minmax, mincw, miniw, minrw, mProb, mQNmod, nb,
     &     nCon, ne, negCon, nextcw, nextiw, nextrw, nkx, nInf, nlocJ,
     &     nlocG, ngObj, nMajor, nName, nnCol, nnCon0, nnCon, nnH,
     &     nnJac, nnObj, nS, Start, tolfac, tolupd, TPivot
      double precision
     &     dnrm1s, fObj, ObjAdd, ObjTru, sInf, xNorm
*     ------------------------------------------------------------------
      integer            FM
      parameter         (FM     = 1)
      integer            TCP
      parameter         (TCP    = 2)
      integer            Cold,       Warm
      parameter         (Cold   = 0, Warm   = 2)
      integer            StdIn
      parameter         (StdIn  = 2)
      integer            HUnset
      parameter         (HUnset = -1 )
      parameter         (tolfac =  66) ! LU factor tolerance.
      parameter         (tolupd =  67) ! LU update tolerance.
      parameter         (lvlSrt =  69) ! = 0(1) => cold(warm) start
      parameter         (lvlHes =  72) ! 0,1,2  => LM, FM, SD Hessian
      parameter         (lprSol =  84) ! > 0    => print the solution
      parameter         (TPivot = 156) ! 0(1) LU part(complete) piv
      parameter         (Htype  = 202) ! Current Hessian type
      parameter         (nName  = 233) ! # of row and col. names
      parameter         (mProb  =  51) ! Problem name
      parameter         (lencw  = 500)
      character          cdummy*8, cw(lencw)*8
      parameter         (cdummy = '-1111111')
      double precision   zero
      parameter         (zero   = 0.0d+0)
*     ------------------------------------------------------------------
      Solver = 'NPOPT '
      INFO   = 0

*     ------------------------------------------------------------------
*     Check memory limits and fetch the workspace starting positions.
*     ------------------------------------------------------------------
      call s2Mem0
     &   ( INFO, Solver, lencw, leniw, lenrw, iw,
     &     mincw, miniw, minrw, maxcw, maxiw, maxrw,
     &     nextcw, nextiw, nextrw )
      if (INFO .gt. 0) go to 999 ! Exit without printing

*     Save the user's option choices  (weird choices get overwritten).
*     Initialize timers and the standard input file.

      call icopy ( 130, iw(51), 1, Useriw, 1 )
      call dcopy ( 130, rw(51), 1, Userrw, 1 )

      call s1time( 0, 0, iw, leniw, rw, lenrw  )
      call s1file( StdIn, iw, leniw )

      lNames  = nextcw - 1   ! No names

*     Check the arguments of npOpt.

      cStart  = 'Cold'          ! Preempted by lvlSrt

      call s3argN
     &   ( inform, cStart, ldA, ldcJ, ldH,
     &     n, nclin, ncnln, iw(nName),
     &     bl, bu, cw(lNames), iState, cMul, iw(lvlSrt), Errors,
     &     iw, leniw, rw, lenrw )
      if (inform .gt. 0) then
         INFO = inform
         go to 800
      end if

*     Load the local problem dimensions.

      nCon    = nclin + ncnln

      if (nCon .eq. 0) then

*        The problem is unconstrained.
*        Include a dummy row of zeros.

         nnCol = 0
         m     = 1
         ne    = 1
      else
         nnCol = n
         m     = nCon
         ne    = m*n
      end if

*     Load the iw array with various problem dimensions.
*     First record problem dimensions for smart users to access in iw.

      iObj    = 0
      nnCon   = ncnln
      nnJac   = nnCol
      nnObj   = n

*     ------------------------------------------------------------------
*     The obligatory call to npInit has already ``unset''
*     the optional parameters.  However, it could not undefine
*     the char*8 options.  Do it now.
*     ------------------------------------------------------------------
      do i = 51, 180
         cw(i)  = cdummy
      end do

*     Set default options that relate specially to npOpt.
*     (Mainly, LU complete pivoting with small threshold pivot tol).

      if (rw(tolfac) .lt. zero) rw(tolfac) = 1.1d+0
      if (rw(tolupd) .lt. zero) rw(tolupd) = 1.1d+0
      if (iw(TPivot) .lt. 0   ) iw(TPivot) = TCP
      if (iw(lvlHes) .lt. 0   ) iw(lvlHes) = FM

      ObjAdd   = zero

*     ------------------------------------------------------------------
*     Load a generic problem name.
*     Check that the optional parameters have sensible values.
*     Delay printing the options until the arguments have been checked.
*     ------------------------------------------------------------------
      cw(mProb) = '     NLP'

      call s8dflt
     &   ( m, n, nnCon, nnJac, nnObj, cw, lencw, iw, leniw, rw, lenrw )
      call s3prtB
     &   ( m, n, nnCon, nnJac, nnObj, iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     Determine storage requirements using the
*     following variables:
*         m,      n,     ne
*         lenR  , maxS , nnL
*         nnObj , nnCon, nnJac
*         negCon
*     All have to be known before calling s8Map.
*     ------------------------------------------------------------------
      nb      = n + m
      nlocJ   = n + 1
      nkx     = nb

      negCon  = max(ncnln*n,1)

*     Allocate arrays that are arguments of s8solv.
*     These are for the data,
*              locJ, indJ, Jcol, bl, bu, Names,
*     and for the solution
*              hs, x, pi, rc, hs.

      lindJ   = nextiw
      llocJ   = lindJ  + ne
      lhs     = llocJ  + nlocJ
      nextiw  = lhs    + nb

      lJcol   = nextrw
      lbl     = lJcol  + ne
      lbu     = lbl    + nb
      lx      = lbu    + nb
      lpi     = lx     + nb
      lrc     = lpi    + m
      nextrw  = lrc    + nb

      maxR    = iw( 52) ! max columns of R.
      maxS    = iw( 53) ! max # of superbasics
      mQNmod  = iw( 54) ! (ge 0) max # of BFGS updates
      minmax  = iw( 87) ! 1, 0, -1  => MIN, FP, MAX

      lenR    = maxR*(maxR + 1)/2  +  (maxS - maxR)
      iw( 20) = negCon
      iw( 28) = lenR

*     Load the iw array with various problem dimensions.

      nnH     = max( nnJac, nnObj )
      ngObj   = nnObj  ! Local nnObj is altered for FP

      iw( 15) = n      ! copy of the number of columns
      iw( 16) = m      ! copy of the number of rows
      iw( 17) = ne     ! copy of the number of nonzeros in Jcol
      iw( 21) = nnJac  ! # of Jacobian  variables
      iw( 22) = nnObj  ! # of objective variables
      iw( 23) = nnCon  ! # of nonlinear constraints
      iw( 24) = nnH    !   max( nnObj, nnJac )
      iw(204) = iObj   ! position of the objective row in J

      iw(nName)  = 1

*     ------------------------------------------------------------------
*     If only a feasible point is requested, save the base point for the
*     objective function:  1/2 || x - x0 ||^2
*     ------------------------------------------------------------------
      FPonly  = minmax .eq. 0
      if ( FPonly ) then
         ngObj  = nnH

         lx0    = nextrw
         lGsav  = lx0    + nnH
         nextrw = lGsav  + nnObj
         minrw  = nextrw - 1
         if (nextrw .le. lenrw) then
            iw(298) = lx0       ! x0(nnL)     = pp starting point
            iw(339) = lGsav     ! Gsav(nnObj) copy of true obj gradient
*           call dcopy ( nnL, rw(lx), 1, rw(lx0), 1 ) ! Done after x def
         end if
      end if

*     ------------------------------------------------------------------
*     Allocate the local arrays for npOpt.
*     s8Map  maps snopt integer and double arrays.
*     s2BMap maps the arrays for the LU routines.
*     s2Mem  checks what space is available and prints any messages.
*     ------------------------------------------------------------------
      call s8Map
     &   ( m, n, negCon, nkx, nnCon, nnJac, ngObj,
     &     lenR, maxR, maxS,  mQNmod, iw(lvlHes),
     &     nextcw, nextiw, nextrw, iw, leniw )
      call s2Bmap
     &   ( m, n, ne, maxS,
     &     nextiw, nextrw, maxiw, maxrw, liwEst, lrwEst, iw, leniw )
      PrtMem = .true.           ! Print all messages in s2Mem
      call s2Mem
     &   ( inform, PrtMem, liwEst, lrwEst,
     &     nextcw, nextiw, nextrw,
     &     maxcw, maxiw, maxrw, lencw, leniw, lenrw,
     &     mincw, miniw, minrw, iw )
      if (inform .ne. 0) then
         INFO = inform
         go to 800
      end if

      iw(256)   = lJcol   ! Jcol(neJ)   = Constraint Jacobian by columns
      iw(257)   = llocJ   ! locJ(n+1)   = column pointers for indJ
      iw(258)   = lindJ   ! indJ(ne) holds the row indices for Jij
      iw(271)   = lbl     ! bl(nb)      = lower bounds
      iw(272)   = lbu     ! bu(nb)      = upper bounds
      iw(299)   = lx      ! x(nb)       = the solution (x,s)
      iw(279)   = lpi     ! pi(m)       = the pi-vector
      iw(280)   = lrc     ! rc(nb)      = the reduced costs
      iw(282)   = lhs     ! the column state vector
      iw(359)   = lNames  ! Names(nName)

      lgObj     = iw(296) ! gObj(nnObj) = Objective gradient
      lfCon     = iw(316) ! fCon (nnCon)  constraints at x
      lgCon     = iw(320) ! gCon (negCon) constraint gradients at x

*     Define the row and column ordering for J.
*     NPOPT  uses natural order throughout, so kx = kxN.

      iw(247) = nkx     ! dimension of kx and its inverse, kxN
      lkx     = iw(251) ! j  = kx (jN) => col j of Jcol is variable jN
      iw(252) = lkx     ! jN = kxN(j ) => col j of Jcol is variable jN

      call s1perm( n, iw(lkx) )
      call s1perm( m, iw(lkx+n) )

      if (iw(lvlSrt) .eq. 0) then
         Start = Cold
      else
         Start = Warm
      end if

*     ------------------------------------------------------------------
*     Initialize some SNOPT arrays that are copied to NPOPT arrays.
*     Build the Jacobian and load the SNOPT arrays.
*     ------------------------------------------------------------------
      nnCon0 = max ( nnCon, 1 )

      call s3iniN
     &   ( Start, n, nb, nnCon0, nnCon, negCon, iw(lhs),
     &     rw(lfCon), rw(lgCon), rw(lgObj), rw(lrc), rw(lx),
     &     iw, leniw, rw, lenrw )
      call s3inN
     &   ( Start, ldA, ldH, m, n, ncLin, nCon, nnCol,
     &     nb, nnCon0, nnCon, iw(lhs), iState,
     &     A, ne, nlocJ, iw(llocJ), iw(lindJ), rw(lJcol),
     &     bl, bu, rw(lbl), rw(lbu), cCon, cMul,
     &     Hess, rw(lpi), x, rw(lx),
     &     iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     Construct column pointers for the nonlinear part of the  Jacobian.
*     ------------------------------------------------------------------
      if (nnCon .gt. 0) then
         llocG = iw(260) ! locG(nlocG) = column pointers for indG
         nlocG = nnJac + 1

         call s8Gloc
     &      ( nnCon, nnJac,
     &        ne, nlocJ, iw(llocJ), iw(lindJ), negCon, nlocG, iw(llocG))
      end if

      if ( FPonly ) call dcopy ( nnH, rw(lx), 1, rw(lx0), 1 )

*     ------------------------------------------------------------------
*     Solve the problem.
*     Tell s8solv that we don't have an initial Hessian.
*     ------------------------------------------------------------------
      iw(Htype)  = HUnset
      lprSav     = iw(lprSol)
      iw(lprSol) = 0

      call s8solv
     &   ( INFO, Solver, Cold,
     &     s0fgN, funcon, funobj,
     &     snLog, snLog2, snSTOP, gotR,
     &     m, n, nb, nnCon, nnJac, ngObj,
     &     iw(nName), iObj, ObjAdd, fObj, ObjTru, nInf, sInf,
     &     ne, nlocJ, iw(llocJ), iw(lindJ), rw(lJcol),
     &     rw(lbl), rw(lbu), cw(lNames),
     &     iw(lhs), rw(lx), rw(lpi), rw(lrc), nMajor, nS,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      iw(lprSol) = lprSav

      Objf   = fObj
      MajIts = nMajor

*     ------------------------------------------------------------------
*     Unload the SNOPT arrays.
*     ------------------------------------------------------------------
      if ( FPonly .and.  nnObj .gt. 0)
     &     call dcopy ( nnObj, rw(lGsav), 1, rw(lgObj), 1 )

      call s3outN
     &   ( ldcJ, ldH, n, ncLin, nCon,
     &     nb, nnCon0, nnCon, iw(lhs), iState,
     &     cCon, cJac, cMul, rw(lfCon), rw(lgCon), rw(lgObj), grad,
     &     Hess, rw(lrc), x, rw(lx),
     &     iw, leniw, rw, lenrw )

      xNorm  = dnrm1s( n, rw(lx), 1 )

      call s3prtN
     &   ( n, (n+nCon), ncLin, nnCon0, ldA, iw(lprSol), xNorm,
     &     iState, A, bl, bu, cCon, cMul, x, rw(lx),
     &     iw, leniw, rw, lenrw )

*     Restore the user's choices of options.

      call icopy ( 130, Useriw, 1, iw(51), 1 )
      call dcopy ( 130, Userrw, 1, rw(51), 1 )

*     Print times for all clocks (if lvlTim > 0).

      call s1time( 0, 2, iw, leniw, rw, lenrw )

      return

*     Local exit messages.

  800 call snWRAP( INFO, Solver, str, str2, iw, leniw )

  999 return

      end ! subroutine npKerN
