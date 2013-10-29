*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  sn80ncon.f
*
*     s8feas   s8FD     s8Fv     s8Fx     s8Gcpy   s8getR   s8Gloc
*     s8Gprd   s8Hfix   s8Infs   s8iQN    s8iQP    s8mrt    s8PPHx
*     s8qpHx   s8rand   s8rc     s8sclJ   s8sInf   s8step   s8sOpt
*     s8wInf
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8feas
     &   ( iExit, MnrLog, lenR, m, maxS, mBS,
     &     n, nb, nnCon0, nnCon, nnL0, nnL, nDegen, nS,
     &     numLC, numLIQ, itn, itnlim, itQP, MnrPrt, sclObj,
     &     tolQP, tolx, nInf, sInf, wtInf, piNorm, rgNorm,
     &     ne, nlocJ, locJ, indJ, Jcol,
     &     hElast, hEstat, hfeas, hs, kBS,
     &     Ascale, bl, bu, blSav, buSav, blBS, buBS,
     &     gBS, pi, R, rc, QPrhs, x0, x, xBS,
     &     iy, iy1, y, y1, y2, y3,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     MnrLog
      integer
     &     iExit, lenR, m, maxS, mBS, n, nb, nnCon0, nnCon, ne, nlocJ,
     &     nnL0, nnL, nDegen, nS, numLC, numLIQ, nInf, itn, itnlim,
     &     itQP, MnrPrt, lencw, leniw, lenrw, locJ(nlocJ), indJ(ne),
     &     kBS(mBS), hfeas(mBS), hEstat(nb), hElast(nb), hs(nb), iy(nb),
     &     iy1(nb), iw(leniw)
      double precision
     &     sclObj, tolQP, tolx, sInf, wtInf, piNorm, rgNorm,
     &     Jcol(ne), Ascale(nb), bl(nb), bu(nb), blSav(nb), buSav(nb),
     &     blBS(mBS), buBS(mBS), gBS(mBS), xBS(mBS), R(lenR), rc(nb),
     &     QPrhs(nnCon0), x0(nb), x(nb), pi(m), y(nb), y1(nb), y2(nb),
     &     y3(nb), rw(lenrw)
      character
     &     cw(lencw)*8

*     ==================================================================
*     s8feas   finds a feasible point for a set of linear constraints.
*     A basis is assumed to be specified by nS, hs(*), x(*) and
*     kBS(m+1:m+nS).  In particular, there must be nS values hs(j) = 2,
*     and the corresponding j's must be listed in kBS(m+1:m+nS).
*     The ordering in kBS matches the reduced Hessian R (if any).
*
*     On entry, blSav and blSav contain copies of the true (possibly
*     scaled) upper and bounds set in s5getB.
*
*      iExit       Result
*      -----       ------
*       >0         Fatal error
*        0         Feasible point found
*
*     11 May 1994: First version of s8feas.
*     19 Aug 1996: First minsum version.
*     05 Feb 1998: Proximal point norm changed to one-norm.
*     23 Dec 1999: Optional Proximal Point methods 0 and  2 added.
*     03 Aug 2003: snPRNT and snEXIT adopted.
*     16 May 2006: Explicit target itQP added.
*     ==================================================================
      character
     &     ProbTag*20, str*80
      logical
     &     Elastc, gotR, needLU, needx
      integer
     &     Hcalls, inform, iObjPP, itQPmax, itQPTgt, j, lvlInf, lvlPPm,
     &     lEmode, lHdx, lrg, lgQP, minimz, mNewSB, MnrHdP, MnrHdS,
     &     mSBsav, nObjP0, nObjPP, nviol, Status, subopt, typeLU
      double precision
     &     blj, buj, eps0, eps2, ObjA, ObjPP, x0j, tolFP, tolQPP,
     &     Hcndbd, Zcndbd
      external
     &     s8PPHx, s8qpHx
*     ------------------------------------------------------------------
      parameter         (Status = 198) ! Status of a call to Hprod
      parameter         (mNewSB =  95) ! max # of new superbasics
      parameter         (MnrHdP = 223) ! >0 => Minor heading for iPrint
      parameter         (MnrHdS = 225) ! >0 => Minor heading for iSumm
      integer            FP,         QPP
      parameter         (FP     = 0, QPP=6)
      integer            BT
      parameter         (BT     = 3)
      integer            Normal
      parameter         (Normal = 0)
      integer            No,         Yes
      parameter         (No     =-1, Yes    = 0)
      double precision   zero,          one
      parameter         (zero = 0.0d+0, one = 1.0d+0)
*     ------------------------------------------------------------------
      lvlPPm    = iw( 79) ! 1(2)-norm proximal point method for x0

      eps0      = rw(  2) ! eps**(4/5)
      eps2      = rw(  4) ! eps**(1/2)
      Hcndbd    = rw( 85) ! bound on the condition of Hz
      Zcndbd    = rw( 86) ! bound on the condition of Z

*     Pointers

      lHdx      = iw(288) ! Hdx(nnQP)   = product of H with  x1 - x
      lgQP      = iw(290) ! gQP(negQP)  = QP gradient
      lrg       = iw(293) ! rg (maxS)   = reduced gradient

      ObjA      = zero
      iObjPP    = 0
      minimz    = 1       ! Local value

      lEmode    = 1       ! Enter elastic mode if infeasible
      lvlInf    = 2       ! In elastic mode, use:
                          !   W1 = 0  W2 = 1   W1*true obj + W2*sInfE
      Elastc    = .false.

      needLU    = .true.
      needx     =  needLU

*     Set the LP rhs to make x satisfy the (relaxed) nonlinear rows.
*     The array  QPrhs  contains the rhs.
*     Use a fairly tight optimality tolerance for phase 1.

      if (nnCon .gt. 0) then
         call dcopy
     &      ( nnCon, x(n+1), 1, QPrhs, 1 )
         call s2Aprd
     &      ( Normal, eps0, ne, nlocJ, locJ, indJ, Jcol,
     &        one, x, n, (-one), QPrhs, nnCon )
      end if

      if (numLIQ .gt. 0  .or.  nInf .gt. 0) then
*        ---------------------------------------------------------------
*        Find a feasible point for the linear constraints.
*        If none exists, minimize the sum of infeasibilities of the
*        linear rows, subject to the column bounds.
*        ---------------------------------------------------------------
         call iload ( numLC, 3, hElast(n+nnCon+1), 1 )

         ProbTag = 'linear rows'
         subopt  = No
         tolFP   = eps2

         call s5LP
     &      ( inform, FP, ProbTag, Elastc, subopt,
     &        MnrLog, needLU, needx,
     &        m, n, nb, nDegen, itQP, itnlim, itn,
     &        lEmode, lvlInf, MnrPrt,
     &        minimz, iObjPP, sclObj, ObjA, tolFP, tolQP, tolx,
     &        nInf, sInf, wtInf, piNorm, rgNorm,
     &        ne, nlocJ, locJ, indJ, Jcol,
     &        hElast, hEstat, hfeas, hs, kBS,
     &        Ascale, bl, bu, blBS, buBS,
     &        gBS, pi, rc, nnCon0, nnCon, QPrhs, x, xBS, x0,
     &        iy, iy1, y, y1, y2,
     &        cw, lencw, iw, leniw, rw, lenrw )

*        Check for trouble in s5LP.
*        iExit        Status
*        -----        ------
*         -3          Too many iterations
*         -2          Phase 1 is unbounded
*         -1          infeasible nonelastics
*          0          infeasibilities minimized
*         >0          Fatal error
*
*        If the linear constraints are infeasible, the sum of
*        infeasibilities will have been minimized.

         if (inform .ne. 0  .or.  nInf .gt. 0) then
            if (inform .gt. 0) then
               iExit = inform   ! Fatal error
            else if (inform .eq. -3) then
               iExit = 31       ! iterations limit
            else if (nInf .gt. 0) then
               iExit = 11       ! infeasible linear constraints
            end if
            if (iExit .ne. 0) go to 800
         end if

*        Now the linear rows are feasible, they are never allowed
*        to be infeasible again.

         call iload ( numLC, 0, hElast(n+nnCon+1), 1 )

*        Print something brief if s5LP didn't already do so.

         if (MnrPrt .ge. 1) then
            write(str, 8000) itn
            call snPRNT( 31, str, iw, leniw )
            call snPRNT( 22, str, iw, leniw )
         end if
      end if

      if (lvlPPm .gt. 0  .and.  nnL .gt. 0) then
*        ===============================================================
*        x  is feasible for the linear constraints.
*        Find a feasible point closest to x0.
*        Minimize norm(x - x0).
*        ===============================================================
         if (MnrPrt .ge. 1) then
            write(str, 8100) itn, lvlPPm
            call snPRNT( 23, str, iw, leniw )
         end if

         if (lvlPPm .eq. 1) then
*           ------------------------------------------------------------
*           Minimize the one-norm of (x-x0) by fixing the nonlinear
*           variables so that bl = x0 = bu.  Any bl or bu that is moved
*           to  x0  is made elastic.
*           ------------------------------------------------------------
            do j = 1, nnL
               blj = bl(j)
               buj = bu(j)

               if (blj .eq. buj) then
*                 Relax
               else
                  x0j       = x0(j)
                  bl(j)     = x0j
                  bu(j)     = x0j
                  hElast(j) = 3

                  if (hs(j) .le. 1) then
                     x(j) = x0j
                  end if
               end if
            end do

            ProbTag    = 'norm(x-x0) problem  '
            iw(MnrHdP) = 1      ! New LP print   header
            iw(MnrHdS) = 1      ! New LP summary header
            needx   = .true.
            subopt  = No
            tolFP   = 1.0d-2    ! Sloppy phase 1 optimality tol for PP.

            call s5LP
     &         ( inform, FP, ProbTag, Elastc, subopt,
     &           MnrLog, needLU, needx,
     &           m, n, nb, nDegen, itQP, itnlim, itn,
     &           lEmode, lvlInf, MnrPrt,
     &           minimz, iObjPP, sclObj, ObjA, tolFP, tolQP, tolx,
     &           nInf, sInf, wtInf, piNorm, rgNorm,
     &           ne, nlocJ, locJ, indJ, Jcol,
     &           hElast, hEstat, hfeas, hs, kBS,
     &           Ascale, bl, bu, blBS, buBS,
     &           gBS, pi, rc, nnCon0, nnCon, QPrhs, x, xBS, x0,
     &           iy, iy1, y, y1, y2,
     &           cw, lencw, iw, leniw, rw, lenrw )

*           Some elastic variables may have moved outside their bounds.
*           Count them.  Reset the true bounds.
*           If necessary,  get feasible again with the normal tolQP.

            nviol = 0
            do j = 1, nnL
               bl(j) = blSav(j)
               bu(j) = buSav(j)

               if (x(j) .lt. bl(j) - tolx  .or.
     &             x(j) .gt. bu(j) + tolx      ) then
                  nviol = nviol + 1
               end if
            end do

*           Check for errors in s5LP.
*           inform values are = -3,-2,-1, 0, >0

            if (inform .ne. 0) then
               if (inform .gt. 0) then
                  iExit = inform ! Fatal error
               else if (inform .eq. -3) then
                  iExit = 31     ! iterations limit
               end if
               if (iExit .ne. 0) go to 800
            end if

            if (inform .eq. 0  .and.  MnrPrt .ge. 1) then
               write(str, 8200) itn, lvlPPm, sInf
               call snPRNT( 33, str, iw, leniw )
               if (nviol .gt. 0) then
                  write(str, 8300) itn
                  call snPRNT( 23, str, iw, leniw )
               end if
            end if

            if (nviol .gt. 0) then
               ProbTag = 'linear rows again   '
               Elastc  = .false.
               needx   = .true.
               subopt  = No
               tolFP   = eps2   ! Revert to accurate phase 1 opt tol

               if (inform .ne. 0) needLU = .true.

               call s5LP
     &            ( inform, FP, ProbTag, Elastc, subopt,
     &              MnrLog, needLU, needx,
     &              m, n, nb, nDegen, itQP, itnlim, itn,
     &              lEmode, lvlInf, MnrPrt,
     &              minimz, iObjPP, sclObj, ObjA, tolFP, tolQP, tolx,
     &              nInf, sInf, wtInf, piNorm, rgNorm,
     &              ne, nlocJ, locJ, indJ, Jcol,
     &              hElast, hEstat, hfeas, hs, kBS,
     &              Ascale, bl, bu, blBS, buBS,
     &              gBS, pi, rc, nnCon0, nnCon, QPrhs, x, xBS, x0,
     &              iy, iy1, y, y1, y2,
     &              cw, lencw, iw, leniw, rw, lenrw )

*              Possible inform values are = -3,-2,-1, 0, >0

               if (inform .ne. 0) then
                  if (inform .gt. 0) then
                     iExit = inform ! Fatal error
                  else if (inform .eq. -3) then
                     iExit = 31 ! iterations limit
                  else if (nInf .gt. 0) then
                     iExit = 11 ! infeasible (should not happen here)
                  end if
                  if (iExit .ne. 0) go to 800
               end if

               if (inform .eq. 0  .and.  MnrPrt .ge. 1) then
                  write(str, 8400) itn, nviol
                  call snPRNT( 23, str, iw, leniw )
               end if
            end if

            nInf  = 0
            sInf  = zero

*           Now the nonlinear variables are feasible, they are never
*           allowed to be infeasible again.

            call iload ( nnL, 0, hElast, 1 )

         else if (lvlPPm .eq. 2) then
*           ------------------------------------------------------------
*           Minimize the two-norm of (x-x0).
*           ------------------------------------------------------------
*           Now the linear rows are feasible, they are never allowed
*           to be infeasible again.

            call iload ( numLC, 0, hElast(n+nnCon+1), 1 )
            nObjPP = 0          ! No explicit gradient in proximal point
            nObjP0 = 1
            gotR   = .false.
            needLU = .false.
            typeLU = BT
            Hcalls = 0

            ProbTag    = 'norm(x-x0) problem  '
            iw(Status) = 1      ! First call to Hprod.
            iw(MnrHdP) = 1      ! Switch to QP print   heading
            iw(MnrHdS) = 1      ! Switch to QP summary heading
            needx   = .false.
            itQPmax = 100       ! Limit the number of minor iterations
            itQPTgt = 100       ! Limit the number of minor iterations
            mSBsav  = iw(mNewSB)
            iw(mNewSB) = 100    ! and the number of new superbasics
            subopt  = Yes
            tolFP   = eps2
            tolQPP  = 1.0d-2     ! Sloppy phase 2 opt tol

            call s5QP
     &         ( inform, QPP, ProbTag, Elastc, subopt,
     &           s8PPHx, s8qpHx, Mnrlog, gotR, needLU, typeLU, needx,
     &           lenR, m, maxS, mBS, n, nb, nDegen, Hcalls,
     &           nnL0, nnL, nObjP0, nObjPP, nnL0, nnL, nS,
     &           itQP, itQPmax, itQPTgt, itn, lEmode, lvlInf, MnrPrt,
     &           minimz, iObjPP, sclObj, ObjA, ObjPP, Hcndbd, Zcndbd,
     &           tolFP, tolQPP, tolx, nInf, sInf, wtInf, piNorm, rgNorm,
     &           ne, nlocJ, locJ, indJ, Jcol,
     &           hElast, hEstat, hfeas, hs, kBS,
     &           Ascale, bl, bu, blBS, buBS,
     &           gBS, rw(lgQP), rw(lgQP), rw(lHdx), y3, pi,R,rc,rw(lrg),
     &           nnCon0, nnCon, QPrhs, nnL0, nnL, x0, x, xBS, x0,
     &           iy, iy1, y, y1, y2,
     &           cw, lencw, iw, leniw, rw, lenrw,
     &           cw, lencw, iw, leniw, rw, lenrw )
            iw(mNewSB) = mSBsav

*           Check for trouble.
*           Possible inform values are = -8(-1)-1, 0, >0

            if (inform .ne. 0  .or.  nInf .gt. 0) then
               if (inform .gt. 0) then
                  iExit = inform ! Fatal LU error
               else if (inform .eq. -3) then
                  iExit = 31    ! iterations limit
               else if (inform .eq. -1  .or.  nInf .gt. 0) then
                  iExit = 11    ! infeasible (should not happen here)
               end if
               if (iExit .ne. 0) go to 800
            end if

*           Note: ObjQP is an updated quantity that may be slightly
*           negative.

            if (MnrPrt .ge. 1) then
               write(str, 8200) itn, lvlPPm, abs(ObjPP)
               call snPRNT( 31, str, iw, leniw )
               call snPRNT( 22, str, iw, leniw )
            end if
         end if ! Proximal Point method 1
      end if ! nnL > 0

  800 return

 8000 format(' Itn', i7, ': Feasible linear rows')
 8100 format(' Itn', i7, ': PP', i1, '.  Minimizing  Norm(x-x0)')
 8200 format(' Itn', i7,
     &       ': PP', i1, '.  Norm(x-x0) approximately minimized  (',
     &               1p, e8.2, ')')
 8300 format(' Itn', i7,
     &       ': PP1.  Making nonlinear variables feasible')
 8400 format(' Itn', i7, ': PP1. ',
     &               i7, ' nonlinear variables made feasible')

      end ! subroutine s8feas

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8FD
     &   ( nnCon0, nnCon, nnObj, itn, cdItns,
     &     centrl, goodG, newG, useFD, info, duInf,
     &     fCon, fObj, iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     centrl, goodG, newG, useFD
      integer
     &     nnCon0, nnCon, nnObj, itn, cdItns, leniw, lenrw,
     &     info(6), iw(leniw)
      double precision
     &     duInf, fCon(nnCon0), fObj, rw(lenrw)

*     ==================================================================
*     s8FD   controls the switch from forward to central differences and
*     vice versa.
*
*     If the forward-difference estimate of the reduced gradient of the
*     Lagrangian is small,  a switch is made to central differences.
*     In this case, the derivatives are recomputed and the QP is solved
*     again.
*
*     On the other hand, if central differences have produced a large
*     reduced-gradient norm, switch back to forward differences.
*
*     31 Mar 2000: First version of s8FD written for SNOPT 6.1.
*     03 Aug 2003: snPRNT adopted.
*     03 Aug 2003: Current version of s8FD.
*     ==================================================================
      character
     &     str*80
      external
     &     dnrm1s
      integer
     &     lvlDif
      double precision
     &     epsrf, cNorm, fdint1, Objsiz, rgNorm, rgTest, dnrm1s
*     ------------------------------------------------------------------
      integer            iFDiff
      parameter         (iFDiff    = 6)
      parameter         (lvlDif    = 181) ! forwd diffs or cntrl diffs
      double precision   zero,          one,          ten
      parameter         (zero = 0.0d+0, one = 1.0d+0, ten = 10.0d+0)
*     ------------------------------------------------------------------
      epsrf     = rw( 73) ! relative function precision.
      fdint1    = rw( 76) ! (1) forwrd diff. interval

      cNorm  = zero
      if (nnCon .gt. 0) cNorm  = dnrm1s( nnCon, fCon, 1 )
      ObjSiz = zero
      if (nnObj .gt. 0) ObjSiz = abs(fObj)

      goodG  = .true.
      rgTest = (one + ObjSiz + cNorm)*epsrf/fdint1
      rgNorm = duInf

      if ( centrl ) then
         if (rgNorm .gt. ten*rgTest  .and.  cdItns .gt. 0) then
            iw(lvlDif) =  1
            centrl     = .false.
            if ( useFD ) then
               info(iFDiff) = 0
            end if
         end if
      else
         if (rgNorm .le.     rgTest) then
            cdItns     = 0
            iw(lvlDif) = 2
            if ( useFD ) then
               goodG   = .false.
               newG    = .true.
               info(iFDiff) = 1
               write(str, 1000) itn
               call snPRNT( 23, str, iw, leniw )
            end if
         end if
      end if

 1000 format( ' Itn', i7, ' -- Central differences invoked.',
     &       '  Small reduced gradient.' )

      end ! subroutine s8FD

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Fv
     &   ( Elastc, n, nnCon, tolz, wtInf, bl, bu, Fv, x, Lmul, Fx )

      implicit
     &     none
      logical
     &     Elastc
      integer
     &     n, nnCon
      double precision
     &     tolz, wtInf, bl(n+nnCon), bu(n+nnCon), x(n+nnCon),
     &     Fx(nnCon), Fv(nnCon), Lmul(nnCon)

*     ==================================================================
*     s8Fv  computes the vector of nonlinear constraint violations:
*        Fv = fCon + A(linear)*x - (nonlinear slacks)
*
*     If the Lagrange multiplier is zero, the violation can be set to
*     any value without changing the merit function.  In this case we
*     try and set the slack so that Fv is zero (subject to the slack
*     being feasible).
*
*     In elastic mode we implicitly adjust the variables v and w such
*     that   c - s(feas) + v - w = 0,  with  v >= 0  and  w >= 0.
*
*     On entry,
*        x   =  the current x.
*        Fx  =  fCon + A(linear)*x,   defined in s8Fx.
*
*     On exit,
*        x   =  x containing the modified slacks.
*        Fv  =  fCon + A(linear)*x -  slacks.
*        Fx  =  unaltered.
*
*     19 Apr 2001: First version based on s8sOpt
*     19 Apr 2001: Current version.
*     ==================================================================
      integer
     &     i, j
      double precision
     &     blj, buj, Fxi, Fvi, FvL, FvU, xj, Lmuli, Lmulv, Lmulw
*     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero = 0.0d+0)
*     ------------------------------------------------------------------
      do i = 1, nnCon
         j     = n + i
         xj    = x(j)
         Fxi   = Fx(i)
         Fvi   = Fxi - xj

         Lmuli = Lmul(i)

         blj   = bl(j)
         buj   = bu(j)

         FvU    = Fxi - buj
         FvL    = Fxi - blj

         Lmulv = abs( wtInf - Lmuli ) ! Multiplier for v in elastic mode
         Lmulw = abs( wtInf + Lmuli ) ! Multiplier for w in elastic mode

         if (     Elastc .and. xj .le. blj .and. Lmulv .le. tolz) then
            if (Fvi .gt. zero) then
               Fvi = max( zero, FvL )
            else
               Fvi = zero
            end if
         else if (Elastc .and. xj .ge. buj .and. Lmulw .le. tolz) then
            if (Fvi .lt. zero) then
               Fvi = min( zero, FvU )
            else
               Fvi = zero
            end if
         else
            if (     Lmuli .le.  tolz  .and.  Fvi .gt. zero) then
               Fvi = max( zero, FvU )
            else if (Lmuli .ge. -tolz  .and.  Fvi .lt. zero) then
               Fvi = min( zero, FvL )
            end if
         end if

         xj    = Fxi - Fvi
         Fv(i) = Fvi
         x(j)  = xj

      end do

      end ! subroutine s8Fv

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Fx
     &   ( n, nnCon, nnJac, tolz,
     &     ne, nlocJ, locJ, indJ, Jcol, fCon, x, Fx )

      implicit
     &     none
      integer
     &     n, nnCon, nnJac, ne, nlocJ, indJ(ne), locJ(nlocJ)
      double precision
     &     tolz, Jcol(ne), x(n+nnCon), fCon(nnCon), Fx(nnCon)

*     ==================================================================
*     s8Fx  defines the nonlinear constraint values
*       Fx  =  true nonlinear slack = fCon + A(linear)*x,
*
*     09 Jan 1992: First version based on Minos routine m8viol.
*     16 Nov 1998: Norm x changed to include only columns.
*     21 Oct 2000: Made compatible with SNOPT 6.1
*     21 Oct 2000: Current version of s8Fx
*     ==================================================================
      integer
     &     nlin
*     ------------------------------------------------------------------
      integer            Normal
      parameter         (Normal = 0)
      double precision   one
      parameter         (one = 1.0d+0)
*     ------------------------------------------------------------------
*     Compute the nonlinear constraint value.
*     Set  Fx  =  fCon + (linear A)*x,   excluding slacks.

      call dcopy
     &   ( nnCon, fCon, 1, Fx, 1 )

      nlin = n - nnJac
      if (nlin .gt. 0) then
         call s2Aprd
     &      ( Normal, tolz,
     &        ne, nlin+1, locJ(nnJac+1), indJ, Jcol,
     &        one, x(nnJac+1), nlin, one, Fx, nnCon )
      end if

      end ! subroutine s8Fx

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Gcpy
     &   ( nnCon, nnJac, ne, nlocJ, locJ, indJ,
     &     neG1, nlocG1, locG1, G1,
     &     neG2, nlocG2, locG2, G2 )

      implicit
     &     none
      integer
     &     nnCon, nnJac, neG1, neG2, nlocG1, nlocG2, ne, nlocJ,
     &     indJ(ne), locJ(nlocJ), locG1(nlocG1), locG2(nlocG2)
      double precision
     &     G1(neG1), G2(neG2)

*     ==================================================================
*     s8Gcpy  copies G1 into G2 when either  G1 or  G2
*     is stored in the upper-left hand corner of J.
*
*     16 Sep 1993: First version.
*     26 Oct 2000: Current version.
*     ==================================================================
      integer
     &     ir, j, k, l1, l2
*     ------------------------------------------------------------------
      do j  = 1, nnJac
         l1 = locG1(j)
         l2 = locG2(j)
         do k  = locJ(j), locJ(j+1)-1
            ir = indJ(k)
            if (ir .gt. nnCon) go to 100
            G2(l2) = G1(l1)
            l1  = l1 + 1
            l2  = l2 + 1
         end do
  100    continue
      end do

      end ! subroutine s8Gcpy

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8getR
     &   ( iExit, Htype, Hcalls, gotR, typeLU, LUreq,
     &     itn, lenR, m, mBS, n, nb,
     &     nnCon0, nnCon, nnH, nS, MjrPrt, minimz, iObj,
     &     U0ii, targtH, targtZ, ne, nlocJ, locJ, indJ, Jcol,
     &     hs, kBS, bl, bu, blBS, buBS, R, QPrhs,
     &     xQP, xBS, iy, iy1, y, y1, y2,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     gotR
      integer
     &     Hcalls, Htype, iExit, iObj, itn, lenR, lencu, leniu,
     &     lenru, lencw, leniw, lenrw, LUreq, m, mBS, nnCon0, nnCon,
     &     n, nb, ne, nlocJ, nnH, nS, MjrPrt, minimz, typeLU,
     &     locJ(nlocJ), indJ(ne), hs(nb), kBS(mBS), iy(nb), iy1(nb),
     &     iu(leniu), iw(leniw)
      double precision
     &     targtH, targtZ, U0ii, Jcol(ne), bl(nb), bu(nb), blBS(mBS),
     &     buBS(mBS), QPrhs(nnCon0), R(lenR), xBS(mBS), xQP(nb),
     &     y(nb), y1(nb), y2(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     s8getR   computes the Cholesky factor of the reduced Hessian.
*
*     On entry, the LU factorization is assumed to be known.
*       gotR = .false.
*
*     iExit      Result
*     -----      ------
*       0        the reduced Hessian has been computed successfully
*      >0        fatal error
*
*     LUreq =  1  Frequency
*     LUreq =  2  LU nonzeros increased
*     LUreq =  3
*     LUreq =  4
*     LUreq =  5  Singular after LU mod
*     LUreq =  6  Unstable LU mod (growth in new column of U)
*     LUreq =  7  Not enough memory
*     LUreq =  8
*     LUreq =  9
*     LUreq = 10  Row error in setx
*     LUreq = 11  Big  dx   in setx
*
*     LUreq = 20
*     LUreq = 21  Iterative refinement failed in QP
*     LUreq = 22  Unbounded QP
*     LUreq = 23
*     LUreq = 24  Small directional derivative in QP
*     LUreq = 25  Ill-conditioned Z in QP
*     LUreq = 26  Indefinite Z'HZ in QP
*     LUreq = 27  R singular after bound swap in QP
*
*     On output,
*     QPerr points to ' ', 't', 'u' or 'w'.
*     QPfea points to ' '  or 'i'.
*
*     14 Mar 2001: First version.
*     03 Aug 2003: snPRNT adopted.
*     29 Jun 2005: Current version of s8getR.
*     ==================================================================
      character
     &     str*80
      external
     &     s8Hwrp, s8Hx
      logical
     &     LUok, needLU, newB, newLU, Rcheck
      integer
     &     condZ, eigH, inform, maxR, nSwap
      double precision
     &     condZ0, eps, flmax, Hdmax, plInfy
*     ------------------------------------------------------------------
      parameter         (condZ  = 192) ! condition estimate of Z
      integer            HUnit
      parameter         (HUnit  = 2)
      double precision   one
      parameter         (one    = 1.0d+0)
*     ------------------------------------------------------------------
      maxR      = iw( 52) ! max columns of R
      eigH      = iw(200) ! =1(0) for definite QP Hessian

      eps       = rw(  1) ! machine precision.  IEEE DP  2.22e-16
      flmax     = rw(  8) ! est. of the largest pos. real

      plInfy    = flmax
      condZ0    = plInfy  ! saved estimate of cond(Z)

      LUok      = .true.

*     ==================================================================
*+    while (LUok  .and. .not. gotR) do
  100 if    (LUok  .and. .not. gotR) then
*     ------------------------------------------------------------------
         inform = 0
         needLU = LUreq .gt. 0

         if ( needLU ) then
            call s2Bfac
     &         ( iExit, typeLU, needLU, newLU, newB,
     &           iObj, itn, MjrPrt, LUreq,
     &           m, mBS, n, nb, nnH, nS, nSwap,
     &           ne, nlocJ, locJ, indJ, Jcol,
     &           kBS, hs, bl, bu, blBS, buBS,
     &           nnCon0, nnCon, QPrhs, xQP, xBS,
     &           iy, iy1, y, y2, iw, leniw, rw, lenrw )
            if (iExit .ne. 0) goto 900
         end if

         if (nS .gt. 0) then
*           ------------------------------------------------------------
*           Compute and factorize  Z'HZ.
*           ------------------------------------------------------------
            call s5HZ
     &         ( inform, s8Hwrp, s8Hx, maxR, lenR,
     &           minimz, m, mBS, n, nb, nnH, nS, Hcalls,
     &           ne, nlocJ, locJ, indJ, Jcol,
     &           Hdmax, rw(condZ), targtZ, kBS, R, y, y1, y2,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )

*           Check for trouble in s5Hz.  Possible exits are:
*           inform  Status
*           ------  ------
*            -1     Ill-conditioned null-space basis
*             0     reduced Hessian computed successfully
*            >0     Fatal error in LU solve

            if (inform .eq. 0) then
               call s5Hfac
     &            ( inform, eigH, itn, lenR, m,
     &              maxR, mBS, nb, nS, targtH, Hdmax,
     &              hs, kBS, iy, bl, bu, blBS, buBS, xQP, xBS, R,
     &              iw, leniw, rw, lenrw )

*              Check for trouble in s5Hfac.
*              inform    Status
*              ------    ------
*                -2      H singular (but should be positive definite)
*                -1      H indefinite
*                 0      normal exit

               if (inform .ne. 0) then

*                 The reduced Hessian is not positive definite.
*                 Reset the H = I if it has not been done already.
*                 Otherwise refactorize B,  possibly with tighter tols.

                  if (Htype .ne. HUnit) then
*                    ---------------------------------------------------
*                    Set unit Hessian.
*                    Z'HZ must be computed again.
*                    ---------------------------------------------------
                     write(str, 1100) itn
                     call snPRNT( 23, str, iw, leniw )
                     Htype = HUnit
                     call s8H0
     &                  ( Htype, nnH, U0ii, iw, leniw, rw, lenrw )
                  else
*                    ---------------------------------------------------
*                    H = I, so Z must be ill-conditioned.
*                    Refactorize B with tighter factor tol.
*                    ---------------------------------------------------
                     write(str, 1200) itn
                     call snPRNT( 23, str, iw, leniw )
                     targtH = one/(eps*eps)
                     call s2tryLU
     &                  ( itn, 26, nS, LUreq, LUok, typeLU,
     &                    iw, leniw, rw, lenrw )
                  end if
               end if

               Rcheck = .false.
               if (Rcheck) then
                  call s6Rchk
     &               ( iExit, s8Hwrp, s8Hx, itn, minimz,
     &                 maxR, lenR, m, mBS, n, nb, Hcalls, nnH, nS,
     &                 ne, nlocJ, locJ, indJ, Jcol,
     &                 kBS, R, y, y1, y2,
     &                 cu, lencu, iu, leniu, ru, lenru,
     &                 cw, lencw, iw, leniw, rw, lenrw )
                  if (iExit .ne. 0) go to 900
               end if

            else if (inform .eq. -1) then

*              Ill-conditioned Z in s5Hz.
*              Refactorize B, possibly with a reduced factor tol.
*              If factor tol is already tight, accept Z, however bad.
*              To avoid repeated factorizations, accept Z if condZ
*              wasn't even reduced by the last factorize.

               write(str, 1300) itn, rw(condZ)
               call snPRNT( 23, str, iw, leniw )

               call s2tryLU
     &            ( itn, 25, nS, LUreq, LUok, typeLU,
     &              iw, leniw, rw, lenrw )

               if (.not. LUok) then
                  targtZ = plInfy
                  LUok   = .true.
               else if (rw(condZ) .lt. condZ0) then
                  condZ0 = rw(condZ)
               else
                  targtZ = plInfy
               end if

            else if (inform .gt. 0) then ! LU error in s5Hz
               iExit = inform
               go to 900
            end if
         end if

         gotR = inform .eq. 0

         go to 100
      end if
*+    end while
*     ------------------------------------------------------------------

      if (.not. gotR) iExit = 44   ! unable to factor Z'Z

  900 return

 1100 format(' Itn', i7, ': Reduced Hessian reset')
 1200 format(' Itn', i7, ': Indefinite reduced Hessian')
 1300 format(' Itn', i7, ': Ill-conditioned QP null-space basis.',
     &                    ' Cond = ', 1p, e8.1)

      end ! subroutine s8getR

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Gloc
     &   ( nnCon, nnJac, ne, nlocJ, locJ, indJ, negCon, nlocG, locG )

      implicit
     &     none
      integer
     &     ne, negCon, nlocG, nlocJ, nnCon, nnJac, indJ(ne)
      integer
     &     locJ(nlocJ), locG(nlocG)

*     ==================================================================
*     s8Gloc  counts the number of nonlinear Jacobian elements and
*     assembles their column pointers in locG.
*
*     29 Oct 2000: First version of s8Gloc.
*     29 Oct 2000: Current version.
*     ==================================================================
      integer
     &     ir, j, k
*     ------------------------------------------------------------------
      negCon  = 0
      locG(1) = 1
      do j = 1, nnJac
         do  k = locJ(j), locJ(j+1)-1
            ir = indJ(k)
            if (ir .gt. nnCon) go to 100
            negCon = negCon + 1
         end do
  100    locG(j+1) = negCon + 1
      end do

      end ! subroutine s8Gloc

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Gprd
     &   ( Task, tolz,
     &     ne, nlocJ, locJ, indJ, negCon, nlocG, locG, gCon,
     &     alpha, x, lenx, beta, y, leny )

      implicit
     &     none
      integer
     &     Task, ne, negCon, nlocG, nlocJ, lenx, leny, indJ(ne),
     &     locJ(nlocJ), locG(nlocG)
      double precision
     &     tolz, alpha, beta, gCon(negCon), x(lenx), y(leny)

*     ==================================================================
*     s8Gprd computes matrix-vector products involving J and x.  The
*     variable task specifies the operation to be performed as follows:
*       task = 'N' (normal)          y := alpha*J *x + beta*y,
*       task = 'T' (transpose)       y := alpha*J'*x + beta*y,
*     where alpha and beta are scalars, x and y are vectors and J is a
*     sparse matrix whose columns are in natural order.
*
*     26 Oct 2000: Current version.
*     ==================================================================
      integer
     &     i, ig, iJ, ir, j
      double precision
     &     alphxj, sum, xj
*     ------------------------------------------------------------------
      integer            Normal,        Transp
      parameter         (Normal = 0,    Transp = 1)
      double precision   zero,          one
      parameter         (zero = 0.0d+0, one = 1.0d+0)
*     ------------------------------------------------------------------
      if (alpha .eq. zero  .and.  beta .eq. one)
     &   return

*     First form  y := beta*y.

      if (beta .ne. one) then
         if (beta .eq. zero) then
            do i = 1, leny
               y(i) = zero
            end do
         else
            do i = 1, leny
               y(i) = beta*y(i)
            end do
         end if
      end if

      if (alpha .eq. zero) then

*        Relax

      else if (alpha .eq. (-one)) then

         if (Task .eq. Normal) then
            do  j = 1, lenx
               xj = x(j)
               if (abs( xj ) .gt. tolz) then
                  ig = locG(j)
                  do iJ = locJ(j), locJ(j+1)-1
                     ir = indJ(iJ)
                     if (ir .gt. leny) go to 100
                     y(ir) = y(ir) - gCon(ig)*xj
                     ig    = ig + 1
                  end do
               end if
  100          continue
            end do

         else if (Task .eq. Transp) then

            do j   = 1, leny
               sum = y(j)
               ig  = locG(j)
               do iJ = locJ(j), locJ(j+1)-1
                  ir = indJ(iJ)
                  if (ir .gt. lenx) go to 200
                  sum = sum - gCon(ig)*x(ir)
                  ig  = ig + 1
               end do
  200          y(j) = sum
            end do
         end if

      else ! General alpha

         if (Task .eq. Normal) then
            do j = 1, lenx
               alphxj = alpha*x(j)
               if (abs( alphxj ) .gt. tolz) then
                  ig  = locG(j)
                  do iJ = locJ(j), locJ(j+1)-1
                     ir = indJ(iJ)
                     if (ir .gt. leny) go to 300
                     y(ir) = y(ir) + gCon(ig)*alphxj
                     ig    = ig + 1
                  end do
               end if
  300          continue
            end do
         else if (Task .eq. Transp) then
            do j   = 1, leny
               sum = zero
               ig  = locG(j)
               do iJ = locJ(j), locJ(j+1)-1
                  ir = indJ(iJ)
                  if (ir .gt. lenx) go to 400
                  sum = sum + gCon(ig)*x(ir)
                  ig  = ig + 1
               end do
  400          y(j) = y(j) + alpha*sum
            end do
         end if ! task .eq. Normal
      end if ! general alpha

      end ! subroutine s8Gprd

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Gsiz
     &   ( m, nnCon, nnJac, ne, nlocJ, locJ, indJ, negCon )

      implicit
     &     none
      integer
     &     m, ne, negCon, nlocJ, nnCon, nnJac, indJ(ne)
      integer
     &     locJ(nlocJ)

*     ==================================================================
*     s8Gsiz  counts the number of nonlinear Jacobian elements.
*
*     04 Nov 2000: First version of s8Gsiz
*     04 Nov 2000: Current version.
*     ==================================================================
      integer
     &     ir, k, last, nlocG
*     ------------------------------------------------------------------
      negCon = 0
      nlocG  = nnJac + 1

      if (nnCon .gt. 0) then
         last = locJ(nlocG) - 1
         if (nnCon .eq. m) then
            negCon = last
         else
            do  k = 1, last
               ir = indJ(k)
               if (ir .le. nnCon) negCon = negCon + 1
            end do
         end if
      end if
      negCon = max( 1, negCon )

      end ! subroutine s8Gsiz

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Hfix
     &   ( nnCon, nnJac, tolz,
     &     ne, nlocJ, locJ, indJ, negCon, nlocG, locG,
     &     ydx, ydxmin, PenUnm, fCon, fCon1, gCon, gCon1,
     &     dx, gd, PenU, v, w )

      implicit
     &     none
      integer
     &     nnCon, nnJac, ne, negCon, nlocG, nlocJ, indJ(ne),
     &     locJ(nlocJ), locG(nlocG)
      double precision
     &     ydx, ydxmin, PenUnm, tolz, gCon(negCon), gCon1(negCon),
     &     fCon(nnCon), fCon1(nnCon), PenU(nnCon), v(nnCon), w(nnCon),
     &     dx(nnJac), gd(nnJac)

*     ==================================================================
*     s8Hfix  attempts to find the a vector xPen  of minimum two-norm
*     such that there exists a BFGS update for the modified Lagrangian
*       La   = f(x) - lambda'(fCon1 - LfCon)
*                   + 1/2 (fCon1 - LfCon)'*diag(PenU)*(fCon1 - LfCon),
*
*     where  LfCon = fCon + J(x1)*dx.
*
*     On entry,
*     dx     is the nonlinear part of the search direction x2 - x1.
*     gd     is the Lagrangian gradient difference.
*     gCon    is the Jacobian at the old x.
*     gCon1    is the Jacobian at the new x.
*     ydx    is the approximate curvature of the Lagrangian.
*     ydxmin   (ydx < ydxmin) is the smallest acceptable approximate
*              curvature.
*
*     On exit,
*     gd     is the augmented Lagrangian gradient difference.
*     PenU     are the penalty parameters.
*     ydx    is unchanged unless gotPen is true, in which case
*              ydx = ydxmin.
*
*     08 Dec 1991: First version based on  Npsol  routine npupdt.
*     26 Oct 2000: Current version of s8Hfix.
*     ==================================================================
      external
     &     ddiv, dnrm2
      logical
     &     gotPen, overfl
      integer
     &     i
      double precision
     &    beta, ddiv, diff, dnrm2, Peni, wi, wmax, wnorm
*     ------------------------------------------------------------------
      integer            Normal,     Transp
      parameter         (Normal = 0, Transp = 1)
      double precision   PenMax
*-->  parameter         (PenMax = 1.0d+5)
*-->  parameter         (PenMax = 1.0d+16)
      parameter         (PenMax = 1.0d+5)
      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one = 1.0d+0)
*     ------------------------------------------------------------------
      overfl    = .false.

*     Try an augmented Lagrangian term to increase ydx.

      PenUnm = zero

*     Compute  v = J1*dx and w = (J2 - J1)*dx = J2*dx - v.

      call s8Gprd
     &   ( Normal, tolz,
     &     ne, nlocJ, locJ, indJ,
     &     negCon, nlocG, locG, gCon,
     &     one, dx, nnJac, zero, v, nnCon )
      call s8Gprd
     &   ( Normal, tolz,
     &     ne, nlocJ, locJ, indJ,
     &     negCon, nlocG, locG, gCon1,
     &     one, dx, nnJac, zero, w, nnCon )

      call daxpy ( nnCon, (-one), v, 1, w, 1 )

*     Compute the difference between c and its linearization.
*     v  =  c - cL = fCon1 - (fCon + J1*s) = fCon1 - fCon - J1*s.

      call daxpy ( nnCon, (-one), fCon1, 1, v, 1 )
      call daxpy ( nnCon,   one , fCon , 1, v, 1 )
      call dscal ( nnCon, (-one),     v, 1 )

*     ---------------------------------------------------------
*     Compute the minimum-length vector of penalty parameters
*     that makes the approximate curvature equal to  ydxmin.
*     ---------------------------------------------------------
*     Use w to hold the constraint on PenU.
*     Minimize            norm(PenU)
*     subject to   ( Sum( w(i)*PenU(i) )  =   const,
*                  (           PenU(i)   .ge. 0.

      wmax = zero
      do i    = 1, nnCon
         wi   = w(i)*v(i)
         wmax = max( wmax, wi )
         w(i) = max( zero, wi )
      end do

      wnorm  = dnrm2 ( nnCon, w, 1 )
      diff   = ydxmin - ydx
      beta   = ddiv  ( wmax*diff, wnorm**2, overfl )
      gotPen = .not. overfl  .and.  wmax .gt. zero
     &                       .and.  beta .lt. PenMax

      if ( gotPen ) then
         beta   = diff/wnorm**2

         do    i = 1, nnCon
            wi   = w(i)
            Peni = beta*wi
            v(i) =       Peni*v(i)
            ydx  = ydx + Peni*wi
            PenU(i) = Peni
         end do
         ydx    = max   ( ydx, ydxmin )
         PenUnm = dnrm2 ( nnCon, PenU, 1 )

*        Update  gd  by the term  (J2' - J1')*v,
*        with v = diag(PenU)*(fCon1 - fCon - J1*s) from above.

         call s8Gprd
     &      ( Transp, tolz,
     &        ne, nlocJ, locJ, indJ,
     &        negCon, nlocG, locG, gCon1,
     &          one , v, nnCon, one, gd, nnJac )
         call s8Gprd
     &      ( Transp, tolz,
     &        ne, nlocJ, locJ, indJ,
     &        negCon, nlocG, locG, gCon,
     &        (-one), v, nnCon, one, gd, nnJac )
      end if ! gotPen

      end ! subroutine s8Hfix

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Infs
     &   ( Elastc, n, nb, nnCon0, nnCon, tolx, wtInf,
     &     prInf, duInf, jprInf, jduInf, bl, bu, Fx, rc, x )

      implicit
     &     none
      logical
     &     Elastc
      integer
     &     n, nb, nnCon0, nnCon, jprInf, jduInf
      double precision
     &     tolx, wtInf, duInf, prInf, bl(nb), bu(nb), rc(nb), x(nb),
     &     Fx(nnCon0)

*     ==================================================================
*     s8Infs computes the maximum primal and dual infeasibilities,
*     using bl, bu, rc, x and the true nonlinear slacks Fxslk.
*     The linear constraints and bounds are assumed to be satisfied.
*     The primal infeasibility is therefore the maximum violation of
*     the nonlinear constraints.
*     The dual infeasibility is the maximum complementarity gap
*     for the bound constraints (with bounds assumed to be no further
*     than 1.0 from each x(j)).
*
*     prInf, duInf   return the max primal and dual infeas.
*
*     20 Feb 1994: First version based on Minos 5.5 routine m8infs.
*     25 Oct 1996: Elastic mode added.
*     29 Apr 2001: Current version.
*     ==================================================================
      integer
     &     i, j
      double precision
     &     dj, slack, tol, viol, v, w
*     ------------------------------------------------------------------
      double precision   zero,           one
      parameter        ( zero = 0.0d+0,  one = 1.0d+0 )
*     ------------------------------------------------------------------
      jprInf = 0
      prInf  = zero
      tol    = tolx

*     See how much  Fx  violates the bounds on the nonlinear slacks.
*     prInf is the maximum violation.

      do i = 1, nnCon
         j     = n + i
         slack = Fx(i)
         viol  = max( zero, bl(j) - slack, slack - bu(j) )
         if (prInf .lt. viol) then
            prInf  = viol
            jprInf = j
         end if
      end do

*     ------------------------------------------------------------------
*     + rc(j)  is the multiplier for lower bound constraints.
*     - rc(j)  is the multiplier for upper bound constraints.
*     duInf is the maximum component-wise complementarity measure.
*     ------------------------------------------------------------------
      jduInf = 0
      duInf  = zero
      do   j = 1, nb
         dj  = rc(j)
         if (dj .ne. zero) then

            if (     dj .gt. zero) then
               dj =   dj * min(  x(j) - bl(j), one )
            else if (dj .lt. zero) then
               dj = - dj * min( bu(j) -  x(j), one )
            end if

            if (duInf .lt. dj) then
               duInf   =  dj
               jduInf  =  j
            end if
         end if ! dj nonzero
      end do

*     ------------------------------------------------------------------
*     Include contributions from the elastic variables.
*     ------------------------------------------------------------------
      if ( Elastc ) then
         do j  = n+1, n+nnCon
            dj = rc(j)
            v  = bl(j) - x (j)
            w  = x (j) - bu(j)

            if      (v .gt. tol) then
               dj = abs(wtInf - dj) * min( v, one )
            else if (w .gt. tol) then
               dj = abs(wtInf + dj) * min( w, one )
            else
               dj = zero
            end if

            if (duInf .lt. dj) then
               duInf   =  dj
               jduInf  =  j
            end if
         end do
      end if

      end ! subroutine s8Infs

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8iQN
     &   ( iExit, info, Htype, Mnrlog, Hcalls, Elastc,
     &     gotR, itn, itQP, lenR, m, maxS, mBS, n, nb,
     &     nnCon0, nnCon, nnObj0, nnObj, nnL0, nnL, nS, nDegen,
     &     MjrPrt, MnrPrt, minimz, iObj,
     &     condHz, sclObj, ObjAdd, ObjQP,
     &     tolFP, tolQPk, tolx, nInf, sInf, wtInf,
     &     U0ii, piNorm, ne, nlocJ, locJ, indJ, Jcol,
     &     hElast, hEstat, hfeas, hs, kBS,
     &     Ascale, bl, bu, blBS, buBS, gBS, gQP, gObj, Hdx,
     &     pBS, pi, R, rc, rg, rg2, QPrhs, x,
     &     xQP, xBS, xQP0, iy, iy1, y, y1, y2,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Mnrlog
      logical
     &     Elastc, gotR
      integer
     &     Hcalls, Htype, iExit, info(6), iObj, itn, itQP, lenR,
     &     lencu, leniu, lenru, lencw, leniw, lenrw, m, maxS, mBS,
     &     nnCon0, nnCon, n, nb, nDegen, ne, nlocJ, nnObj0, nnObj,
     &     nnL0, nnL, nInf, nS, MjrPrt, MnrPrt, minimz,
     &     locJ(nlocJ), indJ(ne), hElast(nb), hEstat(nb), hs(nb),
     &     hfeas(mBS), kBS(mBS), iy(nb), iy1(nb), iu(leniu), iw(leniw)
      double precision
     &     condHz, ObjAdd, ObjQP, piNorm, sclObj, sInf, tolFP, tolQPk,
     &     tolx, U0ii, wtInf, Jcol(ne), Ascale(nb), bl(nb), bu(nb),
     &     blBS(mBS), buBS(mBS), gBS(mBS), gQP(nnL0), gObj(nnObj0),
     &     Hdx(nnL0), pBS(mBS), pi(m), QPrhs(nnCon0), R(lenR), rc(nb),
     &     rg(maxS), rg2(maxS), x(nb), xQP(nb), xBS(mBS), xQP0(nb),
     &     y(nb), y1(nb), y2(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     s8iQN   computes  xQP, the solution of the QP subproblem.
*     By construction, the problem has  nnL  nonlinear variables,
*
*     The SQP base point  x  is not altered.
*
*     On entry, the LU factorization is assumed to be known.
*     The arrays  xBS, blBS and buBS are defined.
*
*     iExit     Status
*     -----     ------
*      >0         Fatal error
*       0         QP solution found
*      -1         Too many iterations
*      -2         Too many superbasics
*
*     LUreq =  1  Frequency
*     LUreq =  2  LU nonzeros increased
*     LUreq =  3
*     LUreq =  4
*     LUreq =  5  Singular after LU mod
*     LUreq =  6  Unstable LU mod (growth in new column of U)
*     LUreq =  7  Not enough memory
*     LUreq =  8
*     LUreq =  9
*     LUreq = 10  Row error in setx
*     LUreq = 11  Big  dx   in setx
*
*     LUreq = 20
*     LUreq = 21  Iterative refinement failed in QP
*     LUreq = 22  Unbounded QP
*     LUreq = 23
*     LUreq = 24  Small directional derivative in QP
*     LUreq = 25  Ill-conditioned Z
*     LUreq = 26  Indefinite Z'HZ in QP
*     LUreq = 27  R singular after bound swap in QP
*
*     On output,
*     QPerr points to ' ', 't', 'u' or 'w'.
*     QPfea points to ' '  or 'i'.
*
*     15 Jun 2001: First version of s8iQN based on s8iQP.
*     31 Jul 2003: dnormj used for norm of the nonlinear pis.
*     03 Aug 2003: snEXIT and snPRNT adopted.
*     26 Dec 2003: calls to s2tryLU added.
*     03 Apr 2004: Current version of s8iQN.
*     ==================================================================
      character
     &     ProbTag*20, str*80
      external
     &     dnormj, s8Hwrp, s8Hx
      logical
     &     done, feasbl, LUok, needLU, needx, newB, newLU, NormIn,
     &     solved
      integer
     &     condZ, inform, itnlim, itQPmax, itQPTgt, lEmode,
     &     linesL, linesS, LUreq, lvlInf, mMinor, MnrHdP, MnrHdS,
     &     ngQP0, ngQP, zngQP, nnH0, nnH, nswap, Status, subopt, typeLU
      double precision
     &     dnormj, flmax, ObjFP, piNnln, plInfy, targtZ, Zcndbd
*     ------------------------------------------------------------------
      parameter         (condZ  = 192) ! condition estimate of Z
      parameter         (Status = 198) ! Status of a call to Hprod
      parameter         (linesL = 220) ! # lines in log     file
      parameter         (linesS = 221) ! # lines in summary file
      parameter         (MnrHdP = 223) ! >0 => Minor heading for iPrint
      parameter         (MnrHdS = 225) ! >0 => Minor heading for iSumm
      integer            iQPfea,     iQPerr
      parameter         (iQPfea = 4, iQPerr = 5)
      integer            BT
      parameter         (BT     = 3)
      integer            FPS,        QPS
      parameter         (FPS    = 4, QPS    = 5)
      integer            No,         Yes
      parameter         (No     =-1, Yes    = 0)
      integer            HUnit
      parameter         (HUnit  = 2)
      double precision   zero
      parameter         (zero   = 0.0d+0)
*     ------------------------------------------------------------------
      itnlim    = iw( 89) ! limit on total iterations
      mMinor    = iw( 91) ! limit on minor iterations

      flmax     = rw(  8) ! est. of the largest pos. real
      Zcndbd    = rw( 86) ! bound on the condition of Z

      iExit     = 0
      targtZ    = Zcndbd
      condHz    = zero
      plInfy    = flmax

      feasbl    = .false. ! Status of the non-elastic variables
      NormIn    = .not. Elastc

      ProbTag    = 'QP subproblem'
      itQP       = 0

      iw(linesL) = 0
      iw(linesS) = 0
      iw(MnrHdP) = 1
      iw(MnrHdS) = 1

      info(iQPerr) = 0
      info(iQPfea) = 0

      nnH    = nnL
      nnH0   = nnL0
      ngQP   = nnL
      ngQP0  = max( ngQP, 1 )

!     Set lEmode to switch to Elastic mode on infeasibility.
!     When in elastic mode, set lvlInf to use the composite objective:
!     w1*Obj + w2*sInf,  with W1 = 0, W2 = wtInf. This minimizes the
!     sum of the infeasibilities of the elastic constraints subject to
!     the nonelastic constraints.

      lEmode = 1
      lvlInf = 2

      typeLU = BT
      LUreq  = 0

*     ==================================================================
*     Find a feasible point for this linearization.
*     If the constraints are linear, x is already feasible.
*     ==================================================================
      if (nnCon .gt. 0) then
*        ---------------------------------------------------------------
*        Find a feasible point.
*        If the reduced Hessian is defined, then it is updated.
*        If the constraints are infeasible, minimize the sum of the
*        elastic variables, subject to keeping the non-elastic variables
*        feasible.  Elastic variables can move outside their bounds.
*        ---------------------------------------------------------------
         zngQP   = 0            ! No objective term
         itQPmax = itnlim
         itQPTgt = itnlim
         subopt  = No
         LUok    = .true.
         done    = .false.

*        ===============================================================
*+       while (.not. done  .and.  LUok) do
  500    if    (.not. done  .and.  LUok) then

            needLU = LUreq .gt. 0

            if ( needLU ) then
               call s2Bfac
     &            ( iExit, typeLU, needLU, newLU, newB,
     &              iObj, itn, MjrPrt, LUreq,
     &              m, mBS, n, nb, nnH, nS, nSwap,
     &              ne, nlocJ, locJ, indJ, Jcol,
     &              kBS, hs, bl, bu, blBS, buBS,
     &              nnCon0, nnCon, QPrhs, xQP, xBS,
     &              iy, iy1, y, y2, iw, leniw, rw, lenrw )
               if (iExit .ne. 0) go to 900
               if (nSwap .gt. 0) gotR = .false.
               LUreq  = 0
            end if

            needx = needLU

            call s5QN
     &         ( inform, FPS, ProbTag, Elastc, subopt,
     &           s8Hwrp, s8Hx, Mnrlog, gotR, needLU, typeLU, needx,
     &           lenR, m, maxS, mBS, n, nb, nDegen, Hcalls,
     &           ngQP0, zngQP, nnObj0, nnObj, nnH0, nnH, nS,
     &           itQP, itQPmax, itQPTgt, itn, lEmode, lvlInf, MnrPrt,
     &           minimz, iObj, sclObj, ObjAdd, ObjFP, condHz,
     &           targtZ, tolFP, tolQPk, tolx, nInf, sInf, wtInf, piNorm,
     &           ne , nlocJ, locJ, indJ, Jcol,
     &           hElast, hEstat, hfeas, hs, kBS,
     &           Ascale, bl, bu, blBS, buBS,
     &           gBS, gObj, gQP, Hdx, pBS, pi, R, rc, rg, rg2, ! y1(m+1),
     &           nnCon0, nnCon, QPrhs, nnL0, nnL, x, xQP, xBS, x,
     &           iy, iy1, y, y1, y2,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )

*           Check for trouble.  Here are the possibilities:
*           inform      Result
*           ------      ------
*            >0         Fatal LU error
*             0         Found a feasible point for the nonelastics
*            -1         The nonelastics are infeasible
*            -2         Phase 1 is unbounded
*            -3         Too many iterations
*            -5         Superbasic limit exceeded
*            -6         Void
*            -7         Void
*            -8         Ill-conditioned Z
*            -9         Too many subspace iterations (should not happen)

            if (inform .gt. 0) then
               iExit = inform   ! Fatal LU error
               go to 900
            end if

            done =       inform .eq.  0  .or.  inform .eq. -3
     &             .or.  inform .eq. -5

            if (.not. done) then
*              =========================================================
*              Trouble.
*              inform = -2 means that the phase 1 was unbounded, which
*                          can only occur if a bad basis gives a large
*                          search direction
*              inform = -1 means that the nonelastics are infeasible,
*                          which should not happen since we already
*                          know a feasible point for the nonelastics.
*              inform = -8 means that  R  is being updated but a crude
*                          estimate of condZ is bigger than targtZ.
*                          Refactorize B, possibly with a reduced factor
*                          tol. If the factor tol is already tight,
*                          accept Z, however bad.
*              =========================================================
               if (inform .eq. -1  .or.  inform .eq. -2) then

*                 Treat both cases as infeasible. Repeatedly refactorize
*                 with tighter tols before declaring LC infeasibility.

                  inform = -1

                  write(str, 1500) itn
                  call snPRNT( 23, str, iw, leniw )
                  call s2tryLU
     &               ( itn, 22, nS, LUreq, LUok, typeLU,
     &                 iw, leniw, rw, lenrw )

               else if (inform .eq. -8) then
                  write(str, 1800) itn, rw(condZ)
                  call snPRNT( 23, str, iw, leniw )
                  call s2tryLU
     &               ( itn, 25, nS, LUreq, LUok, typeLU,
     &                 iw, leniw, rw, lenrw )
                  if (.not. LUok) then
                     targtZ = plInfy
                     LUok   = .true.
                  end if
               end if
            end if
            go to 500
         end if
*+       end while
*        ---------------------------------------------------------------

         if (inform .lt. 0) go to 800 ! Itns or infeasible linear constr

         if (Elastc  .and.  NormIn) then
*           ------------------------------------------------------------
*           The QP switched to elastic mode.
*           The linearized constraints are infeasible.
*           ------------------------------------------------------------
            if (MjrPrt .ge. 1  .or.  MnrPrt .ge. 10) then
               write(str, 1100) itn, wtInf
               call snPRNT( 23, str, iw, leniw )
            end if

            gotR   = .false.
            Htype  = HUnit
            call s8H0
     &         ( Htype, nnH, U0ii, iw, leniw, rw, lenrw )

         else if (MjrPrt .gt. 10  .and.  MnrPrt .gt. 10) then

*           No change in mode.

            if ( Elastc ) then
               write(str, 1200) itn
               call snPRNT( 21, str, iw, leniw )
            else
               write(str, 1300) itn
               call snPRNT( 21, str, iw, leniw )
            end if
         end if
      end if ! nlnCon

*     ------------------------------------------------------------------
*     The inelastic variables (x's and linear slacks) are now feasible.
*     Save them in xQP0 for use with the BFGS update.
*
*     Solve the QP subproblem.
*     Loop back sometimes if we need a BS factorize.
*     ------------------------------------------------------------------
      call dcopy ( nb, xQP, 1, xQP0, 1 )

      feasbl     = .true.       ! the nonelastics are feasible

      iw(Status) = 1            ! First call to Hprod this subproblem
      itQPTgt    = itQP + mMinor
      itQPmax    = itnlim
      LUreq      = 0
      typeLU     = BT
      LUok       = .true.
      done       = .false.
      solved     = .false.

*     ==================================================================
*+    while (.not. (solved  .or.  done)  .and.  LUok) do
  600 if    (.not. (solved  .or.  done)  .and.  LUok) then
*        ---------------------------------------------------------------
*        Refactorize the basis if necessary.
*        ---------------------------------------------------------------
         needLU = LUreq .gt. 0

         if ( needLU ) then
            call s2Bfac
     &         ( iExit, typeLU, needLU, newLU, newB,
     &           iObj, itn, MjrPrt, LUreq,
     &           m, mBS, n, nb, nnH, nS, nSwap,
     &           ne, nlocJ, locJ, indJ, Jcol,
     &           kBS, hs, bl, bu, blBS, buBS,
     &           nnCon0, nnCon, QPrhs, xQP, xBS,
     &           iy, iy1, y, y2, iw, leniw, rw, lenrw )
            if (iExit .ne. 0) go to 900
            if (nSwap .gt. 0) gotR = .false.
            LUreq  = 0
         end if

*        How do we update R if the superbasics change?
*        ---------------------------------------------------------------
*        Solve the QP subproblem using a quasi-Newton method.
*        ---------------------------------------------------------------
         if ( MnrPrt .ge. 10) then
            iw(MnrHdP) = 1      ! QN print   header
            iw(MnrHdS) = 1      ! QN summary header
         end if

!        Set lEmode to switch to Elastic mode on infeasibility.
!        Set lvlInf to use the composite objective  Obj + wtInf*sInf
!        after any switch to elastic mode.

         lvlInf = 1

         if (nnL .gt. 0) then
            subopt = Yes
         else
            subopt = No
         end if

         needx = needLU

         call s5QN
     &      ( inform, QPS, ProbTag, Elastc, subopt,
     &        s8Hwrp, s8Hx, Mnrlog, gotR, needLU, typeLU, needx,
     &        lenR, m, maxS, mBS, n, nb, nDegen, Hcalls,
     &        ngQP0, ngQP, nnObj0, nnObj, nnH0, nnH, nS,
     &        itQP, itQPmax, itQPTgt, itn, lEmode, lvlInf, MnrPrt,
     &        minimz, iObj, sclObj, ObjAdd, ObjQP, condHz,
     &        targtZ, tolFP, tolQPk, tolx, nInf, sInf, wtInf, piNorm,
     &        ne, nlocJ, locJ, indJ, Jcol,
     &        hElast, hEstat, hfeas, hs, kBS,
     &        Ascale, bl, bu, blBS, buBS,
     &        gBS, gObj, gQP, Hdx, pBS, pi, R, rc, rg,  rg2, ! y1(m+1),
     &        nnCon0, nnCon, QPrhs, nnL0, nnL, x, xQP, xBS, x,
     &        iy, iy1, y, y1, y2,
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )

*        iExit       Result
*        -----       ------
*         >0         Fatal LU error
*          0         QP solution found
*         -1         The nonelastics are infeasible
*         -2         The QP subproblem is unbounded
*         -3         Too many iterations
*         -4         Void
*         -5         Too many superbasics
*         -6         Void
*         -7         Void
*         -8         Ill-conditioned Z
*         -9         too many subspace iterations

         if (inform .gt. 0) then
            iExit = inform
            go to 900
         end if

         solved = inform .eq.  0 .or. inform .eq. -9
         done   = inform .eq. -3 .or. inform .eq. -5

         if (done) then
*           ============================================================
*           Relax
*           ============================================================
         else if ( solved ) then
*           ============================================================
*           Finish if there are no large nonlinear pi's.
*           Otherwise, re-solve the QP in elastic mode
*           ============================================================
            if (.not. Elastc) then
               piNnln = dnormj( nnCon, pi, 1 )
               if (piNnln .gt. wtInf) then
                  Elastc = .true.
                  solved = .false.
                  write(str, 1400) itn, wtInf
                  call snPRNT( 23, str, iw, leniw )
               end if
            end if
         else
*           ============================================================
*           Trouble.
*           ============================================================
            if (inform .eq. -1) then
*              ---------------------------------------------------------
*              The nonelastics are infeasible. This should not happen.
*              Phase 1 has already found a feasible point for the
*              nonelastics, so the basis must be ill-conditioned.
*              Refactorize with tighter tols and restart at the known
*              feasible point.  Reduce the feasibility tol to try and
*              prevent repeats.
*              ---------------------------------------------------------
               write(str, 1510) itn
               call snPRNT( 23, str, iw, leniw )
               call s2tryLU
     &            ( itn, 22, nS, LUreq, LUok, typeLU,
     &              iw, leniw, rw, lenrw )

               elastc = .false.
               call dcopy
     &            ( nb, xQP0, 1, xQP, 1 )

            else if (inform .eq. -2) then
*              ---------------------------------------------------------
*              The QP is unbounded.
*              As the Hessian is positive definite, this is probably
*              because of an ill-conditioned reduced Hessian.
*              Reset both the full and reduced Hessian.
*              ---------------------------------------------------------
               write(str, 1600) itn
               call snPRNT( 23, str, iw, leniw )
               if (Htype .ne. HUnit) then
                  gotR  = .false.
                  Htype = HUnit
                  call s8H0
     &               ( Htype, nnH, U0ii, iw, leniw, rw, lenrw )
                  write(str, 1950) itn
                  call snPRNT( 23, str, iw, leniw )
               end if

               call s2tryLU
     &            ( itn, 22, nS, LUreq, LUok, typeLU,
     &              iw, leniw, rw, lenrw )

            else if (inform .eq. -8) then
*              ---------------------------------------------------------
*              condZ > targtZ  while computing the search direction.
*              Refactorize B, possibly with a reduced factor tol. If
*              the factor tol is already tight, accept Z, however bad.
*              ---------------------------------------------------------
               write(str, 1800) itn, rw(condZ)
               call snPRNT( 23, str, iw, leniw )
               call s2tryLU
     &            ( itn, 25, nS, LUreq, LUok, typeLU,
     &              iw, leniw, rw, lenrw )
               if (.not. LUok) then
                  targtZ = plInfy
                  LUok   = .true.
               end if

            else if (inform .eq. -9) then
*              ---------------------------------------------------------
*              Too many CG subspace iterations.
*              ---------------------------------------------------------
               write(str, 1900) itn
               call snPRNT( 23, str, iw, leniw )
               if (Htype .ne. HUnit) then
                  gotR  = .false.
                  Htype = HUnit
                  call s8H0
     &               ( Htype, nnH, U0ii, iw, leniw, rw, lenrw )
                  write(str, 1950) itn
                  call snPRNT( 23, str, iw, leniw )
               end if
            end if
         end if

         go to 600
      end if
*+    end while
*     ------------------------------------------------------------------
  800 if (nInf .gt. 0) info(iQPfea) =  1

      if (inform .eq. 0) then
         info(iQPerr) = max(subopt,0)
      else if (inform .eq. -1) then
         iExit = 15             ! infeasible nonelastics
      else if (inform .eq. -2) then
         info(iQPerr) = 3       ! unbounded subproblem
      else if (inform .eq. -3) then
         info(iQPerr) = 2
         iExit = -1             ! too many iterations
      else if (inform .eq. -5  .and.  feasbl) then
         info(iQPerr) = 5
         iExit = -2             ! superbasic limit
      else if (inform .eq. -5) then
         iExit = 33             ! superbasic limit
      else if (inform .eq. -9) then
*        Relax and hope for the best
      else
         iExit = 44             ! ill-conditioned null-space basis
      end if

  900 return

 1100 format(' Itn', i7, ': Infeasible subproblem.',
     &       ' Elastic mode started with weight = ', 1p, e8.1)
 1200 format(' Itn', i7, ': Feasible QP non-elastics')
 1300 format(' Itn', i7, ': Feasible QP subproblem ')
 1400 format(' Itn', i7, ': Large multipliers.',
     &       ' Elastic mode started with weight = ', 1p, e8.1)
 1500 format(' Itn', i7, ': Infeasible nonelastics in QP feasibility',
     &                    ' phase')
 1510 format(' Itn', i7, ': Infeasible nonelastics in QP optimality',
     &                    ' phase')
 1600 format(' Itn', i7, ': Unbounded QP subproblem')
 1800 format(' Itn', i7, ': Ill-conditioned CG null-space basis.',
     &                    ' Cond = ', 1p, e8.1)
 1900 format(' Itn', i7, ': Too many subspace iterations')
 1950 format(' Itn', i7, ': Hessian reset')

      end ! subroutine s8iQN

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8iQP
     &   ( iExit, info, Htype, Mnrlog, Hcalls, Elastc,
     &     gotR, itn, itQP, lenR, m, maxS, mBS, n, nb,
     &     nnCon0, nnCon, nnObj0, nnObj, nnL0, nnL, nS, nDegen,
     &     MjrPrt, MnrPrt, minimz, iObj, sclObj, ObjAdd, ObjQP,
     &     tolFP, tolQPk, tolx, nInf, sInf, wtInf,
     &     U0ii, piNorm, ne, nlocJ, locJ, indJ, Jcol,
     &     hElast, hEstat, hfeas, hs, kBS,
     &     Ascale, bl, bu, blBS, buBS, gBS, gQP, gObj, Hdx,
     &     pBS, pi, R, rc, rg, QPrhs, x,
     &     xQP, xBS, xQP0, iy, iy1, y, y1, y2,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Mnrlog
      logical
     &     Elastc, gotR
      integer
     &     Hcalls, Htype, iExit, info(6), iObj, itn, itQP, lenR,
     &     lencu, leniu, lenru, lencw, leniw, lenrw, m, maxS, mBS,
     &     nnCon0, nnCon, n, nb, nDegen, ne, nlocJ, nnObj0, nnObj,
     &     nnL0, nnL, nInf, nS, MjrPrt, MnrPrt, minimz, locJ(nlocJ),
     &     indJ(ne), hElast(nb), hEstat(nb), hs(nb), hfeas(mBS),
     &     kBS(mBS), iy(nb), iy1(nb), iu(leniu), iw(leniw)
      double precision
     &     ObjAdd, ObjQP, piNorm, tolFP, tolQPk, tolx, sclObj, sInf,
     &     wtInf, U0ii, Jcol(ne), Ascale(nb), bl(nb), bu(nb),
     &     blBS(mBS), buBS(mBS), gBS(mBS), gQP(nnL0), gObj(nnObj0),
     &     Hdx(nnL0), pBS(mBS), pi(m), QPrhs(nnCon0), R(lenR), rc(nb),
     &     rg(maxS), x(nb), xQP(nb), xBS(mBS), xQP0(nb), y(nb),
     &     y1(nb), y2(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     s8iQP   computes  xQP, the solution of the QP subproblem.
*     By construction, the problem has  nnL  nonlinear variables,
*
*     The SQP base point  x  is not altered.
*
*     On entry, the LU factorization is assumed to be known.
*     The arrays  xBS, blBS and buBS are defined.
*
*     iExit     Status
*     -----     ------
*      >0         Fatal error
*       0         QP solution found
*      -1         Too many iterations
*      -2         Too many superbasics
*
*     LUreq =  1  Frequency
*     LUreq =  2  LU nonzeros increased
*     LUreq =  3
*     LUreq =  4
*     LUreq =  5  Singular after LU mod
*     LUreq =  6  Unstable LU mod (growth in new column of U)
*     LUreq =  7  Not enough memory
*     LUreq =  8
*     LUreq =  9
*     LUreq = 10  Row error in setx
*     LUreq = 11  Big  dx   in setx
*
*     LUreq = 20
*     LUreq = 21  Iterative refinement failed in QP
*     LUreq = 22  Unbounded QP
*     LUreq = 23
*     LUreq = 24  Small directional derivative in QP
*     LUreq = 25  Ill-conditioned Z
*     LUreq = 26  Indefinite Z'HZ in QP
*     LUreq = 27  R singular after bound swap in QP
*
*     On output,
*     QPerr points to ' ', 't', 'u' or 'w'.
*     QPfea points to ' '  or 'i'.
*
*     30 Dec 1991: First version of s8iQP.
*     19 Jul 1997: Thread-safe version.
*     31 Jul 2003: dnormj used for norm of the nonlinear pis.
*     03 Aug 2003: snEXIT and snPRNT  adopted.
*     03 Apr 2004: Current version of s8iQP.
*     ==================================================================
      character
     &     ProbTag*20, str*80
      external
     &     dnormj, s8Hwrp, s8Hx
      logical
     &     done, feasbl, LUok, needLU, needx, NormIn, solved
      integer
     &     condZ, inform, itnlim, itQPmax, itQPTgt, lEmode,
     &     linesL, linesS, LUreq, lvlInf, mMinor, MnrHdP, MnrHdS,
     &     ngQP0, ngQP, zngQP, nnH0, nnH, znnH, subopt, Status, typeLU
      double precision
     &     dnormj, eps, flmax, ObjFP, piNnln, plInfy, targtH, targtZ,
     &     Hcndbd, Zcndbd, rgNorm
*     ------------------------------------------------------------------
      parameter         (condZ  = 192) ! condition estimate of Z
      parameter         (Status = 198) ! Status of a call to Hprod
      parameter         (linesL = 220) ! # lines in log     file
      parameter         (linesS = 221) ! # lines in summary file
      parameter         (MnrHdP = 223) ! >0 => Minor heading for iPrint
      parameter         (MnrHdS = 225) ! >0 => Minor heading for iSumm
      integer            iQPfea,     iQPerr
      parameter         (iQPfea = 4, iQPerr = 5)
      integer            BT
      parameter         (BT     = 3)
      integer            FPS,        QPS
      parameter         (FPS    = 4, QPS    = 5)
      integer            No,         Yes
      parameter         (No     =-1, Yes    = 0)
      integer            HUnit
      parameter         (HUnit  = 2)
      double precision   one
      parameter         (one    = 1.0d+0)
*     ------------------------------------------------------------------
      itnlim    = iw( 89) ! limit on total iterations
      mMinor    = iw( 91) ! limit on minor iterations

      eps       = rw(  1) ! machine precision.  IEEE DP  2.22e-16
      flmax     = rw(  8) ! est. of the largest pos. real
      Hcndbd    = rw( 85) ! bound on the condition of Hz
      Zcndbd    = rw( 86) ! bound on the condition of Z

      iExit     = 0
      plInfy    = flmax
      targtZ    = Zcndbd
      targtH    = Hcndbd

      feasbl    = .false. ! Status of the non-elastic variables
      NormIn    = .not. Elastc

      ProbTag   = 'QP subproblem'

      itQP      = 0

      iw(linesL) = 0
      iw(linesS) = 0
      iw(MnrHdP) = 1
      iw(MnrHdS) = 1

      info(iQPerr) = 0
      info(iQPfea) = 0

      nnH    = nnL
      nnH0   = nnL0
      ngQP   = nnL
      ngQP0  = max( ngQP, 1 )

!     Set lEmode to switch to Elastic mode on infeasibility.
!     When in elastic mode, set lvlInf to use the composite objective:
!     w1*Obj + w2*sInf,  with W1 = 0, W2 = wtInf. This minimizes the
!     sum of the infeasibilities of the elastic constraints subject to
!     the nonelastic constraints.

      lEmode = 1
      lvlInf = 2

      typeLU = BT
      LUreq  = 0

*     ==================================================================
*     Find a feasible point for this linearization.
*     If the constraints are linear, x is already feasible.
*     ==================================================================
      if (nnCon .gt. 0) then
*        ---------------------------------------------------------------
*        Find a feasible point.
*        If the constraints are infeasible, minimize the sum of the
*        elastic variables, subject to keeping the non-elastic variables
*        feasible.  Elastic variables can move outside their bounds.
*        ---------------------------------------------------------------
         zngQP   = 0            ! No objective term in phase 1
         znnH    = 0            ! No Hessian either
         itQPmax = itnlim
         itQPTgt = itnlim
         subopt  = No
         gotR    = .false.
         LUok    = .true.
         done    = .false.

*        ===============================================================
*+       while (.not. done  .and.  LUok) do
  500    if    (.not. done  .and.  LUok) then

            needLU = LUreq .gt. 0
            needx  = needLU

            call s5QP
     &         ( inform, FPS, ProbTag, Elastc, subopt,
     &           s8Hwrp, s8Hx, Mnrlog, gotR, needLU, typeLU, needx,
     &           lenR, m, maxS, mBS, n, nb, nDegen, Hcalls,
     &           ngQP0, zngQP, nnObj0, nnObj, nnH0, znnH, nS,
     &           itQP, itQPmax, itQPTgt, itn, lEmode, lvlInf, MnrPrt,
     &           minimz, iObj, sclObj, ObjAdd, ObjFP, targtH, targtZ,
     &           tolFP, tolQPk, tolx, nInf, sInf, wtInf, piNorm, rgNorm,
     &           ne, nlocJ, locJ, indJ, Jcol,
     &           hElast, hEstat, hfeas, hs, kBS,
     &           Ascale, bl, bu, blBS, buBS,
     &           gBS, gObj, gQP, Hdx, pBS, pi, R, rc, rg,
     &           nnCon0, nnCon, QPrhs, nnL0, nnL, x, xQP, xBS, x,
     &           iy, iy1, y, y1, y2,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )

*           Check for trouble.  Here are the possibilities:
*           inform      Result
*           ------      ------
*            >0         Fatal LU error
*             0         Found a feasible point for the nonelastics
*            -1         The nonelastics are infeasible
*            -2         Phase 1 is unbounded
*            -3         Too many iterations
*            -4         Weak minimizer (after starting elastic mode)
*            -5         Superbasic limit exceeded

            if (inform .gt. 0) then
               iExit = inform   ! Fatal LU error
               go to 900
            end if

            done =       inform .eq.  0  .or.  inform .eq. -3
     &             .or.  inform .eq. -4  .or.  inform .eq. -5

            if (.not. done) then
*              =========================================================
*              Trouble.
*              inform = -2 implies that the phase 1 was unbounded, which
*                          can only occur if a bad basis gives a large
*                          search direction
*              inform = -1 means that the nonelastics are infeasible,
*                          which should not happen since we already
*                          know a feasible point for the nonelastics.
*              =========================================================
*              Treat both cases as infeasible. Repeatedly refactorize
*              with tighter tols before declaring LC infeasibility.

               inform = -1

               write(str, 1500) itn
               call snPRNT( 23, str, iw, leniw )
               call s2tryLU
     &            ( itn, 22, nS, LUreq, LUok, typeLU,
     &              iw, leniw, rw, lenrw )
            end if
            go to 500
         end if
*+       end while
*        ---------------------------------------------------------------

         if (inform .lt. 0) go to 800 ! Itns, infeas inelastcs, or nZlim

         if (Elastc  .and.  NormIn) then
*           ------------------------------------------------------------
*           The QP switched to elastic mode.
*           The linearized constraints are infeasible.
*           ------------------------------------------------------------
            if (MjrPrt .ge. 1  .or.  MnrPrt .ge. 10) then
               write(str, 1100) itn, wtInf
               call snPRNT( 23, str, iw, leniw )
            end if

            gotR   = .false.
            Htype  = HUnit
            call s8H0
     &         ( Htype, nnH, U0ii, iw, leniw, rw, lenrw )

         else if (MjrPrt .gt. 10  .and.  MnrPrt .gt. 10) then

*           No change in mode.

            if ( Elastc ) then
               write(str, 1200) itn
               call snPRNT( 21, str, iw, leniw )
            else
               write(str, 1300) itn
               call snPRNT( 21, str, iw, leniw )
            end if
         end if
      end if ! nlnCon

*     ------------------------------------------------------------------
*     The inelastic variables (x's and linear slacks) are now feasible.
*     Save them in xQP0 for use with the BFGS update.
*
*     Solve the QP subproblem.
*     Loop back if we need better LU factors.
*     ------------------------------------------------------------------
      call dcopy ( nb, xQP, 1, xQP0, 1 )

      feasbl     = .true.       ! the nonelastics are feasible

      iw(Status) = 1            ! First call to Hprod this subproblem
      itQPmax    = itnlim
      itQPTgt    = itQP + mMinor
      LUreq      = 0
      typeLU     = BT
      LUok       = .true.
      done       = .false.
      solved     = .false.

*     ==================================================================
*+    while (.not. (solved  .or.  done)  .and.  LUok) do
  600 if    (.not. (solved  .or.  done)  .and.  LUok) then
*        ---------------------------------------------------------------
*        Compute and factorize the initial Z'HZ.
*        The basis is refactorized if necessary.
*        ---------------------------------------------------------------
         if (.not. gotR  .or.  LUreq .gt. 0) then
            call s8getR
     &         ( iExit, Htype, Hcalls, gotR, typeLU, LUreq,
     &           itn, lenR, m, mBS, n, nb,
     &           nnCon0, nnCon, nnH, nS, MjrPrt, minimz, iObj,
     &           U0ii, targtH, targtZ, ne, nlocJ, locJ, indJ, Jcol,
     &           hs, kBS, bl, bu, blBS, buBS, R, QPrhs,
     &           xQP, xBS, iy, iy1, y, y1, y2,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
            if (iExit .ne. 0) go to 900
            LUreq  = 0
         end if

*        ---------------------------------------------------------------
*        Solve the QP subproblem.
*        ---------------------------------------------------------------
         if ( MnrPrt .ge. 10) then
            iw(MnrHdP) = 1      ! QP print   header
            iw(MnrHdS) = 1      ! QP summary header
         end if

!        Set lEmode to switch to Elastic mode on infeasibility.
!        Set lvlInf to use the composite objective  Obj + wtInf*sInf
!        after any switch to elastic mode.

         lvlInf = 1

         if (nnL .gt. 0) then
            subopt = Yes
         else
            subopt = No
         end if

         needLU = LUreq .gt. 0
         needx  = needLU

         call s5QP
     &      ( inform, QPS, ProbTag, Elastc, subopt,
     &        s8Hwrp, s8Hx, Mnrlog, gotR, needLU, typeLU, needx,
     &        lenR, m, maxS, mBS, n, nb, nDegen, Hcalls,
     &        ngQP0, ngQP, nnObj0, nnObj, nnH0, nnH, nS,
     &        itQP, itQPmax, itQPTgt, itn, lEmode, lvlInf, MnrPrt,
     &        minimz, iObj, sclObj, ObjAdd, ObjQP, targtH, targtZ,
     &        tolFP, tolQPk, tolx, nInf, sInf, wtInf, piNorm, rgNorm,
     &        ne, nlocJ, locJ, indJ, Jcol,
     &        hElast, hEstat, hfeas, hs, kBS,
     &        Ascale, bl, bu, blBS, buBS,
     &        gBS, gObj, gQP, Hdx, pBS, pi, R, rc, rg,
     &        nnCon0, nnCon, QPrhs, nnL0, nnL, x, xQP, xBS, x,
     &        iy, iy1, y, y1, y2,
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )

*        inform      Result
*        -----       ------
*         >0         Fatal LU error
*          0         QP solution found
*         -1         The nonelastics are infeasible
*         -2         The QP subproblem is unbounded
*         -3         Too many iterations
*         -4         The QP subproblem has a weak minimizer
*         -5         Too many superbasics
*         -6         QP Hessian not positive semidefinite
*         -7         Z'g could not be made sufficiently small
*         -8         Ill-conditioned QP null-space basis

         if (inform .gt. 0) then
            iExit = inform
            go to 900
         end if

         solved = inform .eq. 0 .or. inform .eq. -4
         done   = inform .eq.-2 .or. inform .eq. -3 .or. inform .eq. -5

         if (done) then
*           ============================================================
*           Relax
*           ============================================================
         else if ( solved ) then
*           ============================================================
*           Finish if there are no large nonlinear pi's.
*           Otherwise, re-solve the QP in elastic mode
*           ============================================================
            if (.not. Elastc) then
               piNnln = dnormj( nnCon, pi, 1 )
               if (piNnln .gt. wtInf) then
                  Elastc = .true.
                  solved = .false.
                  write(str, 1400) itn, wtInf
                  call snPRNT( 23, str, iw, leniw )
               end if
            end if
         else
*           ============================================================
*           Trouble.
*           ============================================================
            gotR  = .false.

            if (inform .eq. -1) then
*              ---------------------------------------------------------
*              The nonelastics are infeasible. This should not happen.
*              Phase 1 has already found a feasible point for the
*              nonelastics, so the basis must be ill-conditioned.
*              Refactorize with tighter tols and restart at the known
*              feasible point.  Reduce the feasibility tol to try and
*              prevent repeats.
*              ---------------------------------------------------------
               write(str, 1510) itn
               call snPRNT( 23, str, iw, leniw )
               call s2tryLU
     &            ( itn, 22, nS, LUreq, LUok, typeLU,
     &              iw, leniw, rw, lenrw )

               elastc = .false.
               call dcopy
     &              ( nb, xQP0, 1, xQP, 1 )

            else if (inform .eq. -6  .or.  inform .eq. -7) then
*              ---------------------------------------------------------
*              Indefinite Z'HZ  or large Z'g.
*              Most likely an ill-conditioned Z'HZ.
*              Try to reset the Hessian to I.
*              ---------------------------------------------------------
               if (inform .eq. -6) then
                  write(str, 1600) itn
               else
                  write(str, 1700) itn
               end if
               call snPRNT( 23, str, iw, leniw )

               if (Htype .ne. HUnit) then
                  Htype = HUnit
                  call s8H0
     &               ( Htype, nnH, U0ii, iw, leniw, rw, lenrw )
                  write(str, 1900) itn
                  call snPRNT( 23, str, iw, leniw )

               else
                  if (inform .eq. -6) then
                     targtH = one/(eps*eps)
                     LUreq  = 26
                  else if (inform .eq. -7) then
                     LUreq  = 21
                  end if
                  call s2tryLU
     &               ( itn, LUreq, nS, LUreq, LUok, typeLU,
     &                 iw, leniw, rw, lenrw )
               end if

            else if (inform .eq. -8) then
*              ---------------------------------------------------------
*              condZ > targtZ  while forming Z'HZ for a freq. check.
*              Refactorize B, possibly with a reduced factor tol. If
*              the factor tol is already tight, accept Z, however bad.
*              ---------------------------------------------------------
               write(str, 1800) itn, rw(condZ)
               call snPRNT( 23, str, iw, leniw )
               call s2tryLU
     &            ( itn, 25, nS, LUreq, LUok, typeLU,
     &              iw, leniw, rw, lenrw )
               if (.not. LUok) then
                  targtZ = plInfy
                  LUok   = .true.
               end if
            end if
         end if

         go to 600
      end if
*+    end while
*     ------------------------------------------------------------------
  800 if (nInf .gt. 0) info(iQPfea) =  1

      if (inform .eq. 0) then
         info(iQPerr) = max(subopt,0)
      else if (inform .eq. -1) then
         iExit = 15             ! infeasible nonelastics
      else if (inform .eq. -2) then
         info(iQPerr) = 3       ! unbounded subproblem
      else if (inform .eq. -3) then
         info(iQPerr) = 2
         iExit = -1             ! too many iterations
      else if (inform .eq. -4) then
         info(iQPerr) = 4       ! weak QP solution
      else if (inform .eq. -5  .and.  feasbl) then
         info(iQPerr) = 5
         iExit = -2             ! superbasic limit
      else if (inform .eq. -5) then
         iExit = 33             ! superbasic limit
      else
         iExit = 44             ! ill-conditioned null-space basis
      end if

  900 return

 1100 format(' Itn', i7, ': Infeasible subproblem.',
     &       ' Elastic mode started with weight = ', 1p, e8.1)
 1200 format(' Itn', i7, ': Feasible QP non-elastics')
 1300 format(' Itn', i7, ': Feasible QP subproblem ')
 1400 format(' Itn', i7, ': Large multipliers.',
     &       ' Elastic mode started with weight = ', 1p, e8.1)
 1500 format(' Itn', i7, ': Infeasible nonelastics in QP feasibility',
     &                    ' phase')
 1510 format(' Itn', i7, ': Infeasible nonelastics in QP optimality',
     &                    ' phase')
 1600 format(' Itn', i7, ': Indefinite QP reduced Hessian')
 1700 format(' Itn', i7, ': Large QP reduced gradient')
 1800 format(' Itn', i7, ': Ill-conditioned QP null-space basis.',
     &                    ' Cond = ', 1p, e8.1)
 1900 format(' Itn', i7, ': Reduced Hessian reset')

      end ! subroutine s8iQP

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8mrt
     &   ( nnCon, fMrt, gMrt, HMrt, incRun,
     &     penDmp, penMax, PenNrm, Fv, xPen, y, rw, lenrw )

      implicit
     &     none
      logical
     &     incRun
      integer
     &     nnCon, lenrw
      double precision
     &     fMrt, gMrt, HMrt, penDmp, penMax, PenNrm,
     &     Fv(nnCon), xPen(nnCon), y(nnCon), rw(lenrw)

*     ==================================================================
*     s8mrt  computes the contributions to the merit function and its
*     directional derivative from the nonlinear constraints.
*     The penalty parameters  xPen(j)  are increased if
*     the directional derivative is not sufficiently negative.
*
*     On entry:
*         pi     is the vector of  QP  multipliers.
*         Fv     is the violation c(x) + A(linear)x - s,  where
*                s  minimizes the merit function with respect to the
*                nonlinear slacks only.
*
*     30 Dec 1991: First version based on Npsol 4.0 routine npmrt.
*     02 Nov 1996: Multipliers no longer updated here.
*     19 Jul 1997: Thread-safe version.
*     21 Oct 2000: Made compatible with SNOPT 6.1
*     21 Oct 2000: Current version of s8mrt.
*     ==================================================================
      external
     &     ddiv, ddot, dnrm2
      logical
     &     boost, overfl
      integer
     &     i
      double precision
     &     ddiv, ddot, dnrm2, eps0, ppscl, penlty, penMin, penNew,
     &     penOld, rtUndf, xPen0, xPeni, ynorm
*     ------------------------------------------------------------------
      double precision   zero,          half,          two
      parameter         (zero = 0.0d+0, half = 0.5d+0, two = 2.0d+0)
*     ------------------------------------------------------------------
      eps0      = rw(  2)
      rtUndf    = rw( 10)
      xPen0     = rw( 89)

      overfl    = .false.

*     Find the quantities that define  penMin, the vector of minimum
*     two-norm such that the directional derivative is one half of
*     approximate curvature   - (p)'H(p).
*     The factor  rtUndf  tends to keep  xPen  sparse.

      do i = 1, nnCon
         if (abs( Fv(i) ) .le. rtUndf) then
            y(i) = zero
         else
            y(i) = Fv(i)**2
         end if
      end do

      ynorm  = dnrm2 ( nnCon, y, 1 )
      ppscl  = ddiv  ( gMrt + half*HMrt, ynorm, overfl )
      if (abs( ppscl ) .le. penMax  .and.  .not. overfl) then
*        ---------------------------------------------------------------
*        Bounded  penMin  found.  The final value of  xPen(i)  will
*        never be less than  penMin(i).  A trial value  penNew  is
*        computed that is equal to the geometric mean of the previous
*        xPen  and a damped value of penMin.  The new  xPen  is defined
*        as  penNew  if it is less than half the previous  xPen  and
*        greater than  penMin.
*        ---------------------------------------------------------------
         do i = 1, nnCon
            penMin = max( (y(i)/ynorm)*ppscl, zero )
            xPeni  = xPen(i)

            penNew = sqrt( xPeni*(PenDmp + penMin) )
            if (penNew .lt. half*xPeni ) xPeni = penNew
            xPeni   = max (xPeni, penMin)
            xPen(i) = max (xPeni, xPen0 )
         end do

         PenOld  = PenNrm
         PenNrm = dnrm2( nnCon, xPen, 1 )

*        ---------------------------------------------------------------
*        If  IncRun = true,  there has been a run of iterations in
*        which the norm of  xPen  has not decreased.  Conversely,
*        IncRun = false  implies that there has been a run of
*        iterations in which the norm of xPen has not increased.  If
*        IncRun changes during this iteration the damping parameter
*        PenDmp is increased by a factor of two.  This ensures that
*        xPen(j) will oscillate only a finite number of times.
*        ---------------------------------------------------------------
         boost  = .false.
         if (      IncRun  .and.  PenNrm .lt. PenOld) boost = .true.
         if (.not. IncRun  .and.  PenNrm .gt. PenOld) boost = .true.
         if (boost) then
            PenDmp = min( 1/eps0, two*PenDmp )
            IncRun = .not. IncRun
         end if
      end if

*     ------------------------------------------------------------------
*     Compute the new value and directional derivative of the
*     merit function.
*     ------------------------------------------------------------------
      call dcopy ( nnCon, Fv  , 1, y, 1 )
      call ddscl ( nnCon, xPen, 1, y, 1 )

      penlty = ddot  ( nnCon, y, 1, Fv, 1 )
      fMrt   = fMrt  + half*penlty
      gMrt   = gMrt  -      penlty

      end ! subroutine  s8mrt

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8PPHx
     &   ( Hprod, Hcalls, nnH,
     &     x, Hx, Status,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod
      integer
     &     Status, Hcalls, nnH, lencu, leniu, lenru, lencw, leniw,
     &     lenrw, iu(leniu), iw(leniw)
      double precision
     &     Hx(nnH), x(nnH), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     s8PPHx  defines the product  H*x  for the proximal-point QP
*     subproblem of snopt.
*
*     On exit,    Hx   = x.
*
*     23 Oct 1993: First version of s8PPHx.
*     02 Aug 2000: Current version.
*     ==================================================================
      call dcopy ( nnH, x, 1, Hx, 1 )

      end ! subroutine s8PPHx

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8qpHx
     &   ( nnH, x, Hx, Status, cu, lencu, iu, leniu, ru, lenru)

      implicit
     &     none
      integer
     &     nnH, Status, lencu, leniu, lenru, iu(leniu)
      double precision
     &     x(nnH), Hx(nnH), ru(lenru)
      character
     &     cu(lencu)*8

*     ==================================================================
*     s8qpHx is the argument qpHx for s5solv when s5solv is called from
*     one of the snOpt wrappers.
*
*     04 Dec 2004: First version of s8qpHx.
*     04 Dec 2004: Current version of s8qpHx.
*     ==================================================================

*     Relax

      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8rand( leng, neg, g )

      implicit
     &     none
      integer
     &     leng, neg
      double precision
     &     g(leng)

*     ==================================================================
*     s8rand  fills the array g with random numbers.
*
*     15 Nov 1991: First version of s8rand in s8aux.
*     30 Jun 1999: Current version.
*     ==================================================================
      integer
     &     seeds(3)
*     ------------------------------------------------------------------
      if (neg .le. 0) return

      seeds(1) = 1547
      seeds(2) = 2671
      seeds(3) = 3770

      call ddrand( neg, g, 1, seeds )

      end ! subroutine s8rand

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8rc
     &   ( sclObj, minimz, iObj, m, n, nb,
     &     nnObj0, nnObj, nnCon, nnJac, negCon,
     &     ne, nlocJ, locJ, indJ, Jcol,
     &     gObj, gCon, pi, rc )

      implicit
     &     none
      integer
     &     minimz, iObj, m, n, nb, nnObj0, nnObj, nnCon, nnJac,
     &     negCon, ne, nlocJ, indJ(ne), locJ(nlocJ)
      double precision
     &     sclObj, Jcol(ne), gObj(nnObj0), gCon(negCon), pi(m), rc(nb)

*     ==================================================================
*     s8rc   computes reduced costs rc = gObj - ( A  -I )'*pi,
*     using  gCon  as the top left-hand corner of A.
*     gCon, gObj and pi are assumed to exist.
*
*     s8rc   is called by s8SQP.
*
*     28 Sep 1993: First version, derived from m4rc.
*     31 Oct 1996: Min sum option added.
*     30 Oct 2000: Current version of s8rc.
*     ==================================================================
      integer
     &     ir, j, k, l
      double precision
     &     dj, sgnObj
*     ------------------------------------------------------------------
      double precision   zero
      parameter        ( zero = 0.0d+0 )
*     ------------------------------------------------------------------
      l     = 0

      do j  = 1, nnJac
         dj = zero
         do k  = locJ(j), locJ(j+1) - 1
            ir = indJ(k)
            if (ir .le. nnCon) then
               l  = l  + 1
               dj = dj + pi(ir)*gCon(l)
            else
               dj = dj + pi(ir)*Jcol(k)
            end if
         end do
         rc(j) = -dj
      end do

      do j  = nnJac+1, n
         dj = zero
         do k  = locJ(j), locJ(j+1) - 1
            ir = indJ(k)
            dj = dj  +  pi(ir) * Jcol(k)
         end do
         rc(j) = -dj
      end do

      call dcopy ( m, pi, 1, rc(n+1), 1 )

*     Include the nonlinear objective gradient.

      sgnObj = minimz
      if (nnObj .gt. 0) then
         call daxpy ( nnObj, sgnObj, gObj, 1, rc, 1 )
      end if

      if (iObj .gt. 0) rc(n+iObj) =  rc(n+iObj) + sgnObj*sclObj

      end ! subroutine s8rc

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8sclg
     &   ( nnObj, Ascale, gObj, rw, lenrw )

      implicit
     &     none
      integer
     &     nnObj, lenrw
      double precision
     &     Ascale(nnObj), gObj(nnObj), rw(lenrw)

*     ==================================================================
*     s8sclg  scales the objective gradient.
*     s8sclg is called by fgwrap only if modefg = 2.
*     Hence, it is used to scale known gradient elements (if any),
*     but is not called when missing gradients are being estimated
*     by s6dobj.
*
*     17 Feb 1992: First version.
*     16 Jul 1997: Thread-safe version.
*     02 Jan 2001: Current version of s8sclg.
*     ==================================================================
      integer
     &     j
      double precision
     &     gdummy, grad
*     ------------------------------------------------------------------
      gdummy = rw( 69) ! definition of 'unset' value

      do j = 1, nnObj
         grad = gObj(j)
         if (grad .ne. gdummy) gObj(j) = grad*Ascale(j)
      end do

      end ! subroutine s8sclg

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8sclJ
     &   ( nnCon, nnJac, negCon, n, Ascale,
     &     ne, nlocJ, locJ, indJ, gCon, rw, lenrw )

      implicit
     &     none
      integer
     &     n, ne, negCon, nnCon, nnJac, nlocJ, lenrw, indJ(ne),
     &     locJ(nlocJ)
      double precision
     &     Ascale(n+nnCon), gCon(negCon), rw(lenrw)

*     ==================================================================
*     s8sclJ  scales the Jacobian.
*     s8sclJ is called by fgwrap only if modefg = 2.
*     Hence, it is used to scale known gradient elements (if any),
*     but is not called when missing gradients are being estimated
*     by s6dcon.
*
*     17 Feb 1992: First version based on Minos routine m8sclj.
*     16 Jul 1997: Thread-safe version.
*     02 Dec 2001: Current version of s8sclJ.
*     ==================================================================
      integer
     &    ir, j, k, l
      double precision
     &     Cscale, gdummy, grad
*     ------------------------------------------------------------------
      gdummy = rw( 69) ! definition of 'unset' value

      l    = 0
      do j = 1, nnJac
         Cscale = Ascale(j)

         do k = locJ(j), locJ(j+1)-1
            ir     = indJ(k)
            if (ir .gt. nnCon) go to 300
            l      = l + 1
            grad   = gCon(l)
            if (grad .ne. gdummy)
     &         gCon(l)   = grad*cscale/Ascale(n+ir)
         end do
  300 end do

      end ! subroutine s8sclJ

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8sInf
     &   ( n, nb, nnCon, tolx, nInf, sInf, bl, bu, x )

      implicit
     &     none
      integer
     &      n, nb, nnCon, nInf
      double precision
     &     tolx, sInf, bl(nb), bu(nb), x(nb)

*     ==================================================================
*     s8sInf computes the sum of infeasibilities of the nonlinear slacks
*     using bl, bu and x.
*
*     10 Jan 1997: First version of s8sInf.
*     30 Oct 2000: Current version.
*     ==================================================================
      integer
     &     i, j
       double precision
     &     slack, tol, violL, violU
*     ------------------------------------------------------------------
      double precision   zero
      parameter        ( zero = 0.0d+0 )
*     ------------------------------------------------------------------
      nInf   = 0
      sInf   = zero
      tol    = tolx

*     See how much  x(n+1:n+nnCon) violates its bounds.

      do i = 1, nnCon
         j     = n + i
         slack = x(j)
         violL = bl(j) - slack
         violU = slack - bu(j)
         if (violL .gt. tol  .or.  violU .gt. tol) then
            nInf = nInf + 1
            sInf = sInf + max (violL, violU )
         end if
      end do

      end ! subroutine s8sInf

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8step
     &   ( centrl, usefLS, nb, nnCon, nnObj, nMajor,
     &     nSkip, step, stepmn, steplm, stepmx, tolz, xdNorm, xNorm,
     &     bl, bu, x, dx, iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     centrl, usefLS
      integer
     &     nb, nnCon, nnObj, nMajor, nSkip, leniw, lenrw, iw(leniw)
      double precision
     &     step, stepmn, steplm, stepmx, tolz, xdNorm, xNorm,
     &     bl(nb), bu(nb), x(nb), dx(nb), rw(lenrw)

*     ==================================================================
*     s8step  finds the maximum, minimum and initial value for the
*     linesearch step.
*
*     For problems with nonlinear constraints, the maximum step stepmx
*     is one.  If there are only linear constraints the maximum step is
*     the largest step such that x + step*dx  reaches one of its bounds.
*
*     All step sizes are subject to the user-specified limit  steplm.
*
*     04 Dec 1992: First version of s8step based on npsol routine npalf.
*     31 Mar 2000: Updated for SNOPT 6.1.
*     19 Mar 2001: Current version.
*     ==================================================================
      external
     &     ddiv
      logical
     &     switch, overfl
      integer
     &     Htype, j, gotFD
      double precision
     &     bigdx, fdint1, pivot, pivabs, res, stepQP, tolpiv,
     &     tolp, xdlim, ddiv
*     ------------------------------------------------------------------
      integer            HUnit
      parameter         (HUnit = 2)
      double precision   zero,          one
      parameter         (zero = 0.0d+0, one = 1.0d+0)
*     ------------------------------------------------------------------
      tolpiv    = rw( 60) ! excludes small elements of y.
      bigdx     = rw( 72) ! unbounded step.
      fdint1    = rw( 76) ! (1) forwrd diff. interval
      xdlim     = rw( 80) ! Step limit

      gotFD     = iw(183) ! > 0 => some differences needed
      Htype     = iw(202) ! Current approximate Hessian type

      overfl    = .false.

*     ==================================================================
*     switch  indicates if there is an option to switch to
*             central differences to get a better search direction.
*     stepQP  is the step predicted by the QP subproblem (usually 1).
*     stepmx  is the largest feasible steplength subject to a
*             user-defined limit, bigdx, on the change in  x.
*     step    is initialized subject to a user-defined limit, xdlim.
*     ==================================================================
      if (nnCon .eq. 0  .and.  nnObj .eq. 0) then  ! LP !!
         step   = one
         stepmn = one
         steplm = one
         stepmx = one
      else
         switch = gotFD .gt. 0  .and.  .not. centrl

         stepmn = zero
         if (usefLS  .and.  switch) then
            stepmn = fdint1*(one + xNorm) / xdNorm
         end if

         stepQP = one
         if (nnCon .gt. 0  .and.  (nSkip .eq. 0  .or.
     &                             Htype .ne. HUnit)) then
            stepmx = one
         else
            tolp   = tolpiv*xdNorm
            stepmx = ddiv  ( bigdx, xdNorm, overfl )
            step   = stepmx
            j      = 1

*+          while (j .le. nb  .and.  step .gt. stepQP) do
  100       if    (j .le. nb  .and.  step .gt. stepQP) then
               pivot   = dx(j)
               pivabs  = abs( pivot )
               if (pivabs .gt. tolp) then
                  if (pivot  .le. zero  ) then
                     res    = x(j) - bl(j)
                     if (step*pivabs .gt. res) step = res / pivabs
                  else
                     res    = bu(j) - x(j)
                     if (step*pivabs .gt. res) step = res / pivabs
                  end if
               end if
               j = j + 1
               go to 100
*+          end while
            end if

            step   = max( step, stepQP )
            if (step .lt. stepQP + tolz) step = stepQP

            stepmx = step
         end if

         steplm = ddiv( (one+xNorm)*xdlim, xdNorm, overfl )
         if (nMajor .le. 1)
     &   steplm = min (steplm, ddiv( one, xdNorm, overfl ))
         stepmx = min (            steplm, stepmx)
         step   = min (            steplm, one   )
      end if

      end ! subroutine s8step

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8sOpt
     &   ( n, nb, nnCon, piNorm, tolz, wtInf, hEstat,
     &     bl, bu, Fv, x, Lmul, xPen, Fx )

      implicit
     &     none
      integer
     &     n, nb, nnCon, hEstat(nb)
      double precision
     &     piNorm, tolz, wtInf, bl(nb), bu(nb), x(nb),
     &     xPen(nnCon), Fx(nnCon), Fv(nnCon), Lmul(nnCon)

*     ==================================================================
*     s8sOpt computes the vector of nonlinear constraint violations:
*        Fv = fCon + A(linear)*x - (optimal nonlinear slacks)
*
*     The optimal nonlinear slacks are computed as follows:
*     (1) Feasible  nonlinear slacks are adjusted so that they minimize
*         the merit function subject to  x  and  Lmul  being held
*         constant.
*     (2) Infeasible slacks are compared with the true nonlinear slacks,
*         and, if necessary, they are adjusted so that the sum of
*         infeasibilities is reduced.
*
*     If Lmul is zero, the violation can be set to any value without
*     changing the merit function.  In this case we choose the slack to
*     so that  the violation is zero (subject to the constraints above).
*
*     On entry,
*        x   =  the current x.
*        Fx  =  fCon + A(linear)*x,   defined in s8Fx.
*
*     On exit,
*        x   =  x containing the optimal slacks.
*        Fv  =  fCon + A(linear)*x - optimal slacks.
*        Fx  =  unaltered.
*
*     09 Jan 1992: First version based on Npsol routine npslk.
*     09 Oct 1996: First infeasible slack version.
*     28 Jul 2003: Test hEstat for slacks that are allowed to move.
*     29 Apr 2001: Current version.
*     ==================================================================
      integer
     &     i, j, jEs
      double precision
     &     blj, buj, con, dvMax, Lmuli, tol, vi, vL, vU, vLow,
     &     vUpp, xj, xPeni
*     ------------------------------------------------------------------
      double precision   zero,          one,          factor
      parameter         (zero = 0.0d+0, one = 1.0d+0, factor = 10.0d+0)
*     ------------------------------------------------------------------
      tol = tolz*piNorm

      do i = 1, nnCon
         j     = n + i
         con   = Fx(i)
         xj    = x(j)
         vi    = con - xj

         xPeni = xPen(i)
         Lmuli = Lmul(i)

         blj   = bl(j)
         buj   = bu(j)

         vU    = con - buj
         vL    = con - blj

         jEs   = hEstat(j)

*        ---------------------------------------------------------------
*        Redefine  xj  so that it minimizes the merit function
*        subject to upper and lower bounds determined by the current
*        multipliers.  For computational convenience (but not clarity),
*        instead of checking that  xj  is within these bounds, the
*        violation  vi = c - xj  is checked against  vLow  and  vUpp,
*        the violations at the upper and lower bounds on xj.
*        ---------------------------------------------------------------
*        First, impose artificial bounds (tbl, tbu).

         dvMax = factor*(one + abs( vi ))
         vLow  = vi - dvMax
         vUpp  = vi + dvMax

         if      (jEs .eq. 1  .and.  xj .le. blj) then
*           ------------------------------------------------------------
*           This slack is at or below its lower bound in elastic mode.
*           ------------------------------------------------------------
            if (     Lmuli .lt. zero) then

*              xj is eligible to increase.
*              Require                  bl <=  xj <= min( bu,tbu ).

               vLow  = max( vU, vLow )
               vUpp  = vL

            else if (Lmuli .gt. zero) then

*              xj is eligible to decrease and violate its lower bound.
*              Require              -infty <=  xj <= bl

               Lmuli = Lmuli - wtInf
               vLow  = vL

            else

*              xj can either increase or decrease.
*              Require              -infty <=  xj <= min( bu,tbu ).

               vLow  = max( vU, vLow )
            end if

         else if (jEs .eq. 2  .and.  xj .ge. buj) then
*           ------------------------------------------------------------
*           This slack is at or above its upper bound in elastic mode.
*           ------------------------------------------------------------
            if (     Lmuli .gt. zero) then

*              xj is eligible to decrease.
*              Require      max( bl, tbl ) <=  xj <= bu.

               vLow  = vU
               vUpp  = min( vL, vUpp )

            else if (Lmuli .lt. zero) then

*              xj is eligible to increase and violate its upper bound.
*              Require                  bu <=  xj <= +infty

               Lmuli = Lmuli + wtInf
               vUpp  = vU
            else

*              xj can either increase or decrease.
*              Require      max( bl, tbl ) <=  xj <= +infty

               vUpp  = min( vL, vUpp )
            end if

         else
*           ------------------------------------------------------------
*           Feasible slack.  xj can move either way.
*           ------------------------------------------------------------
*              Require      max( bl, tbl ) <=  xj <= min( bu,tbu ).

            vLow  = max( vU, vLow )
            vUpp  = min( vL, vUpp )
         end if

         if (abs( Lmuli ) .le. tol) then
            vi = min( max( zero, vLow ), vUpp )

         else if (xPeni .ge. tolz) then
            if (Lmuli .ge. xPeni*vUpp) then
               vi = vUpp
            else if (Lmuli .le. xPeni*vLow) then
               vi = vLow
            else
               vi = Lmuli / xPeni
            end if
         end if

         xj    = con - vi
         Fv(i) = vi
         x(j)  = xj

      end do

      end ! subroutine s8sOpt

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8wInf
     &   ( job, boostd, itn, gNorm, wtInf0, wtInf, wtMax,
     &     weight, wtFac, wtScal, iw, leniw )

      implicit
     &     none
      logical
     &     boostd
      integer
     &     job, itn, leniw, iw(leniw)
      double precision
     &     gNorm, wtInf0, wtInf, wtMax, weight, wtFac, wtScal

*     ==================================================================
*     s8wInf  initializes or updates the elastic weight  wtInf.
*     The elastic weight is given by  wtInf = wtScal*weight,
*     where wtScal is some scale-dependent quantity (fObj here).
*     wtInf is increased by redefining weight as weight*wtFac, where
*     wtFac is a constant factor.
*
*     weight, wtFac and wtScal are 'saved' local variables.
*
*     20 Feb 1997: First version of s8wInf.
*     27 Apr 2001: wtMax introduced as parameter instead of local.
*     03 Aug 2003: snPRNT adopted.
*     03 Aug 2003: Current version of s8wInf.
*     ==================================================================
      character
     &     str*80
      double precision
     &     newWt
*     ------------------------------------------------------------------
      double precision   ten
      parameter         (ten   = 10.0d+0)
      integer            SetWt,     Boost
      parameter         (SetWt = 0, Boost = 1)
*     ------------------------------------------------------------------
      if (job .eq. SetWt) then

*        Set the weight.
*        weight is the ``unscaled'' weight on the infeasibilities.
*        wtScal is a scale factor based on the current gradient.

         wtScal = gNorm
         wtFac  = ten
         weight = wtInf0
         wtInf  = wtScal*weight

      else if (job .eq. Boost) then

*        If possible, boost the weight.

         newWt  = min( wtFac*weight, wtMax )
         boostd = newWt .gt. weight

         if ( boostd ) then
            weight = newWt
            wtInf  = weight*wtScal
            wtFac  = ten*wtFac
            write(str, 1000) itn, wtInf
            call snPRNT( 23, str, iw, leniw )
         end if
      end if

      return

 1000 format(' Itn', i7, ': Elastic weight increased to ', 1p, e11.3)

      end ! subroutine s8wInf

