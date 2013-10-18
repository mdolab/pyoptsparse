*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  sn85Hess.f
*
*     s8getH   s8H0     s8HQN    s8Hwrp   s8Hx   s8xHx     s8Hupd
*     s8x1
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8getH
     &   ( nnH, lenH, U, H, y, y1 )

      implicit
     &     none
      integer
     &     lenH, nnH
      double precision
     &     H(lenH), U(lenH), y(nnH), y1(nnH)

*     ==================================================================
*     s8getH  computes the product H = U'U, where  U is the Cholesky
*     factor of the approximate Hessian of the Lagrangian.  The matrix
*     U is stored by rows in the one-dimensional array  U.
*     lenH defines the length of U.  lenH must be at least
*     nnH*(nnH + 1)/2.  The result is stored by columns in the upper
*     triangular array H.
*
*     03 Sep 2006: First version of s8getH.
*     03 Sep 2006: Current version.
*     ==================================================================
      integer
     &     j, jthcol
*     ------------------------------------------------------------------
      integer            WithU,      WithUt
      parameter         (WithU  = 0, WithUt = 1)
      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one   = 1.0d+0)
*     ------------------------------------------------------------------
*     Compute the product y1 = U'Uy, where  U is an upper-
*     ------------------------------------------------------------------
      jthcol = 1

      do     j = 1, nnH
         jthcol = jthcol + j - 1

         call dload
     &      ( nnH, zero, y, 1 )
         y(j) = one

         call s6Rprd
     &      ( WithU , nnH, nnH, lenH, U,  y, y1 )
         call s6Rprd
     &      ( WithUt, nnH, nnH, lenH, U, y1,  y )
         call dcopy
     &      ( j, y, 1, H(jthcol), 1 )

      end do

      end ! subroutine s8getH

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8H0
     &   ( Htype, nnH, U0pre, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Htype, nnH, leniw, lenrw, iw(leniw)
      double precision
     &     U0pre, rw(lenrw)

*     ==================================================================
*     s8H0    initializes the BFGS approximate Hessian.
*
*     s8H0   calls one of the Hessian routines s9LMH0, s9FMH0, s9SDH0
*     according to the value of the option lvlHes.
*     Each of these routines defines a particular form of the Hessian.
*     At the moment the options are:
*        lvlHes = 0   Limited-Memory (LM) BFGS Hessian  (the default).
*        lvlHes = 1   Full-Memory    (FM) BFGS Hessian.
*
*     On entry, the value of Htype is as follows:
*
*       Htype
*       -----
*        -1      H undefined.
*         0      H is an approx. Hessian of the form defined by  lvlHes.
*         1      H is a diagonal matrix.
*         2      H is an identity matrix.
*
*     19 Jul 1995: First version of s8H0.
*     12 Jan 1996: Full memory Hessian option added.
*     18 Feb 2001: H stored in product form.
*     11 Apr 2001: Current version.
*     ==================================================================
      integer
     &     lvlHes, nQNmod, lU0, lHd, lU, lenU
      double precision
     &     Hcndbd
*     ------------------------------------------------------------------
      integer            LM   ,      FM
      parameter         (LM     = 0, FM     = 1)
      parameter         (nQNmod = 381) ! # of updates since last reset
*     ------------------------------------------------------------------
      Hcndbd    = rw( 85) ! bound on the condition of Hz
      lvlHes    = iw( 72) ! LM, FM or Exact Hessian
      lHd       = iw(347) ! Diagonal of BFGS Hessian

      if      (lvlHes .eq. LM) then
*        -----------------------
*        Limited memory Hessian.
*        -----------------------
         lU0        = iw(346) ! Initial approximate Hessian diagonal

         call s9LMH0
     &      ( Htype, nnH, Hcndbd, U0pre, rw(lHd),
     &        rw(lU0), iw(nQNmod) )

      else if (lvlHes .eq. FM) then
*        -----------------------
*        Full memory Hessian.
*        -----------------------
         lU         = iw(391) !
         lenU       = iw(392) !

         call s9FMH0
     &      ( Htype, nnH, Hcndbd, U0pre, rw(lHd),
     &        lenU, rw(lU), iw(nQNmod) )

      end if

      end ! subroutine s8H0

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8HQN
     &   ( iExit, fgwrap, fgcon, fgobj,
     &     useFD, Htype, QPtype, info,
     &     lenR, m, mBS, n, nb,
     &     nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj,
     &     nS, nMajor, nSkip, U0ii,
     &     step, minimz, dxHdx,
     &     RtRmod, gotR, incRun, PenDmp, PenMax,
     &     fObj, fCon, gCon, gObj, fCon1, gCon1, gObj1,
     &     ne, nlocJ, locJ, indJ, Jcol, negCon, nlocG, locG,
     &     kBS, bl, bu, dx, dg, Udx, Hdx, Lmul1,
     &     R, x, x1, xQP0, xPen, y, y1, y2, y3,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     fgwrap, fgcon, fgobj
      logical
     &     gotR, incRun, useFD
      integer
     &     Htype, iExit, info(6), lencu, lencw,
     &     leniu, leniw, lenru, lenrw, lenR, mBS,
     &     minimz, m, n, nb, ne, negCon, nlocG,
     &     nlocJ, nMajor, nnCon0, nnCon, nnJac, nnL, nnObj0,
     &     nnObj, nS, nSkip, QPtype, RtRmod,
     &     kBS(mBS), locJ(nlocJ), indJ(ne),
     &     locG(nlocG), iu(leniu), iw(leniw)
      double precision
     &     dxHdx, U0ii, PenDmp, PenMax, step,
     &     bl(nb), bu(nb),
     &     dg(nnL), dx(nnL), Hdx(nnL),
     &     Jcol(ne), Lmul1(nnCon0), fObj, fCon(nnCon0), gCon(negCon),
     &     gObj(nnObj0), fCon1(nnCon0), gCon1(negCon), gObj1(nnObj0),
     &     R(lenR), Udx(nnL),
     &     x(n), x1(nnL), xPen(nnCon0), xQP0(nb),
     &     y(nb), y1(nb), y2(nb), y3(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     s8HQN  does the quasi-Newton update with vectors
*        dx = x1 - x   and   dg = gL(x1) - gL(x).
*
*     On entry:
*      xQP is the QP solution.
*
*     23 Apr 1999: First version of s8HQN,
*     18 Feb 2001: LM H stored in product form.
*     12 Oct 2003: snEXIT and SNPRNT adopted
*     10 Jan 2005: FM H stored in product form.
*     31 Mar 2005: Current version.
*     ==================================================================
      external
     &     ddot, ddiv, dnrm2
      logical
     &     nlnCon, overfl, updatd
      integer
     &     inform, kFac, lvlSrt, maxR,
     &     mSkip, nBS, nzero
      double precision
     &     ddot, ddiv, dnrm2, eps, eps0, eps1, gLnorm,
     &     PenUnm, rdxHdx, rydx, rnnL, sgnObj, sNorm, U0scal,
     &     xPen0, ydx, ydxmin
*     ------------------------------------------------------------------
      integer            QPChol
      parameter         (QPChol = 0)
      integer            Transp
      parameter         (Transp = 1)
      integer            HNorml
      parameter         (HNorml = 0)
      integer            iQNtyp,     iModfy
      parameter         (iQNtyp = 1, iModfy = 2)
      integer            HOT
      parameter         (HOT    = 3)
      integer            mGap
      parameter         (mGap   = 2)
      double precision   tolg,             tolg2
      parameter         (tolg   =  1.0d-3, tolg2  = 1.0d-1)
      double precision   U0max,            U0min
      parameter         (U0max  =  1.0d+1, U0min  = 1.0d-2)
      double precision   zero,             one
      parameter         (zero   =  0.0d+0, one    = 1.0d+0)
*     ------------------------------------------------------------------
      maxR   = iw( 52) ! max columns of R
      kFac   = iw( 59) ! factorization frequency
      mSkip  = iw( 67) ! # largest allowable  nSkip
      lvlSrt = iw( 69) ! = 0:1:2:3 => cold:warm:basis:hot start

      eps    = rw(  1) ! machine precision.  IEEE DP  2.22e-16
      eps0   = rw(  2) ! eps**(4/5)
      eps1   = rw(  3) ! eps**(2/3)
      xPen0  = rw( 89) ! initial penalty parameter.

      iExit  = 0
      overfl = .false.
      nlnCon = nnCon  .gt. 0

      nBS    = m + nS
      sgnObj = minimz

      info(iQNtyp) = 0
      info(iModfy) = 0

      ydx    = zero

*     ---------------------------------------------------------------
*     Compute  dx = x1 - x  and  dg = gL1 - gL.
*     Compute the approx. curvature ydx and new scale factor U0.
*     ---------------------------------------------------------------
      call dcopy ( nnL,          x1, 1, dx, 1 )
      call daxpy ( nnL, (-one),  x , 1, dx, 1 )
      call dscal ( nnL,   step, Hdx, 1 )
      call dscal ( nnL,   step, Udx, 1 )
      dxHdx   = dxHdx*step*step

      if (nnObj .gt. 0) then
         call dcopy ( nnObj, gObj1, 1, dg, 1 )
         if (minimz .lt. 0) call dscal ( nnObj, sgnObj, dg, 1 )
      end if

      nzero = nnL - nnObj
      if (nzero .gt. 0) call dload ( nzero, zero, dg(nnObj+1), 1 )

      if (nnCon  .gt. 0) then
         call s8Gprd
     &      ( Transp, eps0,
     &        ne, nlocJ, locJ, indJ, negCon, nlocG, locG, gCon1,
     &        (-one), Lmul1, nnCon, one, dg, nnJac )
      end if

*     gLnorm = dnormi( nnL, dg, 1 )
      gLnorm = dnrm2 ( nnL, dg, 1 )

      if (nnObj .gt. 0) then
         call daxpy ( nnObj, (-sgnObj), gObj , 1, dg, 1 )
      end if

      if (nnCon  .gt. 0) then
         call s8Gprd
     &      ( Transp, eps0,
     &        ne, nlocJ, locJ, indJ, negCon, nlocG, locG, gCon,
     &        one , Lmul1, nnCon, one, dg, nnJac )
      end if

      ydx  = ddot ( nnL, dg, 1, dx, 1 )

      if (nMajor .eq. 1  .and.  lvlSrt .ne. HOT) then
*        ===============================================================
*        First iteration. Do not attempt a BFGS update.
*        Use latest curvature information to get a better scaled H.
*        ===============================================================
         if (gLnorm .gt. zero) then
            rnnL = nnL
            U0ii = sqrt(gLnorm/sqrt(rnnL))
         else
            U0ii = one
         end if
         gotR   = .false.
         U0ii   = min( max( U0ii, U0min ), U0max )
         call s8H0
     &      ( Htype, nnL, U0ii, iw, leniw, rw, lenrw )
      else
*        ===============================================================
*        Except on the first iteration, attempt a BFGS update.
*        Compute the smallest allowable curvature.
*        If the update cannot be done, s8x1 attempts to find a
*        modified update using  dx = x1 - x defined with a new x.
*        Arrays fCon, gCon and gObj must be redefined at the new x.
*        ===============================================================
         sNorm  = dnrm2 ( nnL, dx, 1 )
         U0ii   = ddiv  ( sqrt(abs(ydx)), sNorm, overfl )
         U0ii   = min   ( max( U0ii, U0min ), U0max )

         PenUnm = zero
         ydxmin = tolg*dxHdx
         updatd =  dxHdx .gt. zero   .and.
     &             ( ydx .ge. ydxmin  .or.   ydx .ge. eps1)

         if (nlnCon  .and. .not. updatd) then
*           ------------------------------------------------------------
*           Redefine  x, dx, Hdx and dg.
*           The problem functions are recomputed at x.
*           ------------------------------------------------------------
            call s8x1
     &         ( inform, fgwrap, fgcon, fgobj, useFD,
     &           n, nb, nnCon0, nnCon, nnJac, nnObj0, nnObj, nnL,
     &           minimz, step, dxHdx, ydx,
     &           fObj, fCon, gCon, gObj, gCon1, gObj1,
     &           ne, nlocJ, locJ, indJ, negCon, nlocG, locG,
     &           bl, bu, dx, dg, Udx, Hdx, Lmul1, y1, x, xQP0, y,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
            if (inform .gt. 0) then
               iExit = inform   ! User wants to stop
               go to 999
            end if

            ydxmin = tolg*dxHdx
            updatd = dxHdx .gt. zero    .and.
     &               ( ydx .ge. ydxmin  .or.   ydx .ge. eps1)

            if ( updatd ) then
               info(iModfy) = 1
            end if

            if (.not. updatd  .and.  dxHdx .gt. zero  ) then
*              ---------------------------------------------------------
*              If all else fails, attempt to update the Hessian of
*              the augmented Lagrangian.
*              The target ydx is defined via tolg2.
*              ---------------------------------------------------------
               ydxmin = tolg2*dxHdx
               call s8Hfix
     &            ( nnCon, nnJac, eps0,
     &              ne, nlocJ, locJ, indJ, negCon, nlocG, locG,
     &              ydx, ydxmin, PenUnm, fCon, fCon1, gCon, gCon1,
     &              dx, dg, y, y1, y2 )

               updatd = ydx .ge. ydxmin

               if ( updatd ) then
                  info(iModfy) = 2
               end if
            end if
         end if ! nlnCon

         if ( updatd ) then
*           ------------------------------------------------------------
*           Update the approximate Hessian using (dg,Hdx).
*           If there are no nonlinear constraints,  apply the update
*           to the reduced Hessian.
*           ------------------------------------------------------------
            nSkip = 0

            if (ydx .ge. ydxmin  .and.  Htype .eq. HNorml) then
               info(iQNtyp) = 1
            else
               info(iQNtyp) = 2
            end if

            rydx   = sqrt(   ydx )
            rdxHdx = sqrt( dxHdx )

            call s8Hupd
     &         ( info(iQNtyp), Htype, nnL,
     &           U0ii, U0scal, rydx, rdxHdx, dx, Hdx, dg,
     &           iw, leniw, rw, lenrw )

            if (     QPtype .eq. QPChol) then
               gotR = gotR  .and.  nnCon  .eq. 0
     &                      .and.  nS     .gt. 0
     &                      .and.  Htype  .eq. HNorml
     &                      .and.  RtRmod .lt. kfac
            else
               gotR = gotR  .and.  nS     .gt. 0
     &                      .and.  Htype  .eq. HNorml
            end if

            if ( gotR ) then
               call s6Rupd
     &            ( info(iQNtyp), maxR, lenR, m, n, nBS, nnL, nS,
     &              U0scal, rdxHdx, ne, nlocJ, locJ, indJ, Jcol,
     &              kBS, dg, Hdx, R, y3, y, y2,
     &              iw, leniw, rw, lenrw )
               RtRmod = RtRmod + 1
            else
               RtRmod = 0
            end if
         else
*           ------------------------------------------------------------
*           No suitable update pair (dg,Hdx) could be found.
*           Skip the update.  Too many skips and we reset.
*           ------------------------------------------------------------
            nSkip  = nSkip  + 1

*           Apply all updates to H and discard the off-diagonals.

            if (mod( nSkip, mSkip  ) .eq. 0) then
               call s8H0
     &            ( Htype, nnL, U0ii, iw, leniw, rw, lenrw )

               if (mod( nSkip, mGap*mSkip ) .eq. 0) then
*                 ------------------------------------------------------
*                 Reset the multipliers and penalty parameters
*                 ------------------------------------------------------
                  incRun = .true.
                  PenDmp = one
                  PenMax = one / eps
                  call dload ( nnCon, xPen0, xPen , 1 )
                  call dload ( nnCon, zero , Lmul1, 1 )
               end if
            end if
         end if
      end if ! nMajor > 1

  999 return

      end ! subroutine s8HQN

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Hwrp
     &   ( Hprod, Hcalls, nnH,
     &     x, Hx, Status,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod
      integer
     &     Hcalls, nnH, lencu, leniu, lenru, lencw, leniw, lenrw,
     &     Status, iu(leniu), iw(leniw)
      double precision
     &     Hx(nnH), x(nnH), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     s8Hwrp wraps Hprod, which multiplies the QP Hessian H by the
*     vector  x.   It is called by the QP solver.
*
*     On entry:
*        Status  = 0  => a normal call for H*x.
*        Status  = 1  => the first entry for a given QP.
*        Status ge 2  => last call for a given QP. Status = 2+iExit.
*
*     On exit:
*        Status lt 0   the user wants to stop.
*
*     03 Nov 2000: First version of s8Hwrp.
*     03 Nov 2000: Current version.
*     ==================================================================
      call Hprod
     &   ( Hcalls, nnH, x, Hx, cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine s8Hwrp

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Hx
     &   ( Hcalls, nnH, x, Hx, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Hcalls, nnH, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     Hx(nnH), x(nnH), rw(lenrw)
      character
     &     cw(lencw)*8

*     ==================================================================
*     s8Hx  multiplies the QP Hessian  H by the vector  x.
*     It is used to define Hx for the QP subproblem.
*
*     This routine is called by a general QP solver, which will rescale
*     Hx by  sgnObj  when maximizing.
*
*     s8Hx calls one of the Hessian routines s9LMH, s9FMH, s9SDH, ...
*     according to the value of the options lvlDer and lvlHes.
*     Each of these routines defines a particular form of the Hessian.
*     At the moment the options are:
*
*        lvlHes = LM      Limited-Memory (LM) BFGS  (the default).
*        lvlHes = FM      Full-Memory    (FM) BFGS
*        lvlHes = Exact   FD or exact Hessian
*
*     30 Dec 1991: First version of s8Hx.
*     12 Jan 1996: Full memory Hessian option added.
*     04 Apr 1999: Exact and FD Hessian option added.
*     18 Feb 2001: LM H stored in product form.
*     10 Jan 2005: FM H stored in product form.
*     15 Jan 2005: Current version.
*     ==================================================================
      integer
     &     lenU, lS, lU, lU0, lUx, lV, lvlHes, minimz, mQNmod
      double precision
     &     sgnObj
*     ------------------------------------------------------------------
      integer            LM,         FM
      parameter         (LM     = 0, FM     = 1)
      integer            nQNmod
      parameter         (nQNmod = 381)
*     ------------------------------------------------------------------
      lvlHes    = iw( 72) ! LM, FM or Exact Hessian
      minimz    = iw(199) ! (-1)(+1)    => (max)(min)
      lUx       = iw(345) ! Ux(nnL)     = product of U with x

      sgnObj    = minimz

      if      (lvlHes .eq. LM) then
*        -----------------------
*        Limited memory Hessian.
*        -----------------------
         mQNmod    = iw( 54) ! (ge 0) max # of BFGS updates
         lU0       = iw(346) ! Initial diagonal Hessian
         lS        = iw(401) ! sk's for BFGS products: (I + sk*vk')
         lV        = iw(402) ! vk's for BFGS products: (I + sk*vk')

         call s9LMHx
     &      ( nnH, x, rw(lUx), Hx,
     &        mQNmod, iw(nQNmod), rw(lU0), rw(lS), rw(lV) )

      else if (lvlHes .eq. FM) then
*        -----------------------
*        Full memory Hessian.
*        -----------------------
         lU        = iw(391) !
         lenU      = iw(392) !

         call s9FMHx
     &      ( nnH, x, rw(lUx), Hx,
     &        lenU, rw(lU) )
      end if

      if (minimz .lt. 0)
     &   call dscal
     &     ( nnH, sgnObj, Hx, 1 )

      Hcalls = Hcalls + 1

      end ! subroutine s8Hx

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8xHx
     &   ( nnH, x, Ux, Hx, xHx, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     nnH, leniw, lenrw, iw(leniw)
      double precision
     &     xHx, Hx(nnH), Ux(nnH), x(nnH), rw(lenrw)

*     ==================================================================
*     s8xHx  computes x'Hx and Hx, where H = U'U.
*
*     s8xHx calls one of the Hessian routines s9LMH, s9FMH, s9SDH, ...
*     according to the value of the options lvlDer and lvlHes.
*     Each of these routines defines a particular form of the Hessian.
*     At the moment the options are:
*
*        lvlHes = LM      Limited-Memory (LM) BFGS  (the default).
*        lvlHes = FM      Full-Memory    (FM) BFGS
*        lvlHes = Exact   FD or exact Hessian
*
*     10 Jan 2005: First version of s8Hx based on s8Hx
*     18 Feb 2001: LM H stored in product form.
*     10 Jan 2005: FM H stored in product form.
*     15 Jan 2005: Current version.
*     ==================================================================
      external
     &     ddot
      integer
     &     lenU, lvlHes, lU, lU0, lS, lV, mQNmod
      double precision
     &     ddot
*     ------------------------------------------------------------------
      integer            LM,         FM
      parameter         (LM     = 0, FM     = 1)
      integer            nQNmod
      parameter         (nQNmod = 381)
*     ------------------------------------------------------------------
      lvlHes    = iw( 72) ! LM, FM or Exact Hessian

      if      (lvlHes .eq. LM) then
*        -----------------------
*        Limited memory Hessian.
*        -----------------------
         mQNmod    = iw( 54) ! (ge 0) max # of BFGS updates
         lU0       = iw(346) ! Initial diagonal Hessian
         lS        = iw(401) ! sk's for BFGS products: (I + sk*vk')
         lV        = iw(402) ! vk's for BFGS products: (I + sk*vk')

         call s9LMHx
     &      ( nnH, x, Ux, Hx,
     &        mQNmod, iw(nQNmod), rw(lU0), rw(lS), rw(lV) )

      else if (lvlHes .eq. FM) then
*        -----------------------
*        Full memory Hessian.
*        -----------------------
         lU        = iw(391) ! U(lenU), dense Hessian factor
         lenU      = iw(392) !

         call s9FMHx
     &      ( nnH, x, Ux, Hx,
     &        lenU, rw(lU) )
      end if

      xHx = ddot  ( nnH, Ux, 1, Ux, 1 )

      end ! subroutine s8xHx

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Hupd
     &   ( Update, Htype, nnH,
     &     U0pre, U0scal, rydx, rdxHdx, dx, Hdx, y,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Update, Htype, nnH, leniw, lenrw, iw(leniw)
      double precision
     &     U0pre, U0scal, rydx, rdxHdx, dx(nnH), Hdx(nnH), y(nnH),
     &     rw(lenrw)

*     ==================================================================
*     s8Hupd  applies the pair of vectors that define the BFGS update
*     or self-scaled BFGS update.
*
*     On entry:
*
*     Hdx   contains  H  times the difference x1 - x.
*     y     contains the gradient  difference g1 - g.
*
*     On exit:
*     Hdx   contains  H  times the difference x1 - x(new).
*     y     contains the gradient  difference g1 - g(new).
*
*     s8Hupd calls one of the Hessian routines s9LMH, s9FMH, s9SDH, ...
*     according to the value of the option lvlHes.
*     At the moment the options are:
*
*        lvlHes = LM      Limited-Memory (LM) BFGS  (the default).
*        lvlHes = FM      Full-Memory    (FM) BFGS
*        lvlHes = Exact   FD or exact Hessian
*
*     19 Jul 1995: First version of s8Hupd.
*     12 Jan 1996: Full-memory Hessian option added.
*     18 Feb 2001: LM H stored in product form.
*     12 Jan 2005: FM H stored in product form.
*     16 Jan 2005: Current version.
*     ==================================================================
      integer
     &     lvlHes, lenU, mQNmod, lHd, lU0, lU, lUdx, lS, lV
      double precision
     &     Hcndbd
*     ------------------------------------------------------------------
      integer            LM        , FM
      parameter         (LM     = 0, FM = 1)
      integer            nQNmod
      parameter         (nQNmod = 381) ! # of updates since last reset
*     ------------------------------------------------------------------
      Hcndbd    = rw( 85) ! bound on the condition of Hz
      lvlHes    = iw( 72) ! LM, FM or Exact Hessian
      lHd       = iw(347) ! Diagonal of BFGS Hessian

      if      (lvlHes .eq. LM) then
         mQNmod    = iw( 54) ! (ge 0) max # of BFGS updates
         lU0       = iw(346) ! Initial diagonal U such that H = U'U
         lS        = iw(401) ! sk's for BFGS products: (I + sk*vk')
         lV        = iw(402) ! vk's for BFGS products: (I + sk*vk')

         call s9LMup
     &      ( Update, Htype, nnH, mQNmod, iw(nQNmod), Hcndbd,
     &        U0pre, U0scal, rydx, rdxHdx, rw(lHd), Hdx, y,
     &        dx, rw(lU0), rw(lS), rw(lV) )

      else if (lvlHes .eq. FM) then
         mQNmod    = iw( 54) ! (ge 0) max # of BFGS updates
         lUdx      = iw(345) ! Ux(nnL)      = product of U with x
         lU        = iw(391) ! U(lenU), full-memory BFGS Hessian H = U'U
         lenU      = iw(392) !

         call s9FMup
     &      ( Update, Htype, nnH, mQNmod, iw(nQNmod), Hcndbd,
     &        U0pre, U0scal, rydx, rdxHdx, rw(lHd), Hdx, y,
     &        rw(lUdx), lenU, rw(lU) )

      end if

      end ! subroutine s8Hupd

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8x1
     &   ( iExit, fgwrap, fgcon, fgobj, useFD,
     &     n, nb, nnCon0, nnCon, nnJac, nnObj0, nnObj, nnL,
     &     minimz, step, dxHdx, ydx,
     &     fObj, fCon, gCon, gObj, gCon1, gObj1,
     &     ne, nlocJ, locJ, indJ, negCon, nlocG, locG,
     &     bl, bu, dx, dg, Udx, Hdx, Lmul1, tdx, x, xQP0, y,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     fgwrap, fgcon, fgobj
      logical
     &     useFD
      integer
     &     iExit, lencu, lencw, leniu, leniw, lenru, lenrw, minimz, n,
     &     nb, ne, negCon, nlocG, nlocJ, nnCon0, nnCon, nnJac, nnL,
     &     nnObj0, nnObj, locJ(nlocJ), indJ(ne), locG(nlocG),
     &     iu(leniu), iw(leniw)
      double precision
     &     dxHdx, fObj, step, ydx, bl(n), bu(n), dg(nnL), dx(nnL),
     &     Hdx(nnL), Lmul1(nnCon0), fCon(nnCon0), gCon(negCon),
     &     gObj(nnObj0), gCon1(negCon), gObj1(nnObj0), tdx(nnL),
     &     Udx(nnL), x(n), xQP0(nnL), y(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     s8x1   redefines the quantities x and dx, Hdx and dg  used for the
*     quasi-Newton update.  The problem functions are recomputed at x.
*
*     The new  x1  is  x1 + step*(xQP0 - x1),  where xQP0 is a
*     (nonelastic) feasible point from the QP subproblem.
*
*     s8x1 is always called with nnL > 0.
*
*     02 Dec 1994: First version of s8x1.
*     20 Jul 1998: s8x1 made self-contained
*     24 Aug 1998: Fixed bug found by Alan Brown at Nag.
*                  FD derivatives now computed correctly.
*                  Parameter useFD added.
*     11 Oct 1998: Facility to combine funobj and funcon added.
*     12 Oct 2003: snEXIT and snPRNT adopted.
*     14 Jan 2005: Argument Udx added for call to s8xHx.
*     01 Apr 2005: Current version of s8x1.
*     ==================================================================
      external
     &     ddot
      logical
     &     nlnCon, nlnObj
      integer
     &     modefg, nzero, Status
      double precision
     &     ddot, eps, eps0, sgnObj
*     ------------------------------------------------------------------
      integer            Transp
      parameter         (Transp = 1)
      double precision   zero,            one
      parameter         (zero  =  0.0d+0, one = 1.0d+0)
*     ------------------------------------------------------------------
      eps    = rw(  1) ! machine precision.  IEEE DP  2.22e-16
      eps0   = rw(  2) ! eps**(4/5)

      iExit  = 0
      modefg = 2
      sgnObj = minimz

      nlnCon = nnCon  .gt. 0
      nlnObj = nnObj  .gt. 0

*     Save dx in case a better update pair cannot be found.

      call dcopy ( nnL, dx, 1, tdx, 1 )

*     -------------------------------------------------
*     dx =  dx - step*y,  with  y = xQP0 - x
*     -------------------------------------------------
      call daxpy ( nnL, (-one ), x   , 1, xQP0, 1 )
      call daxpy ( nnL, (-step), xQP0, 1,   dx, 1 )

*     -------------------------------------------------
*     Compute the minimum curvature.
*     If nnL < n, dxHdx may be zero (or negative fuzz).
*     -------------------------------------------------
      call s8xHx
     &   ( nnL, dx, Udx, Hdx, dxHdx, iw, leniw, rw, lenrw )

      if (dxHdx .ge. eps) then
*        -----------------------------------------------
*        Redefine  x  as   x + step*y  (y held in xQP0.)
*        Evaluate the functions at the new x.
*        -----------------------------------------------
         call daxpy ( nnL, step, xQP0, 1, x, 1 )

         Status = 0
         call fgwrap
     &      ( iExit, modefg, Status, nlnCon, nlnObj,
     &        n, negCon, nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj,
     &        fgcon, fgobj, x,
     &        ne , nlocJ, locJ, indJ,
     &        fCon, fObj, gCon, gObj,
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )
         if (iExit .gt. 0) return

         if (iExit .eq. 0  .and.  useFD) then
            call s6fdG
     &         ( iExit, n, negCon,
     &           nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj,
     &           fgwrap, fgcon, fgobj,
     &           bl, bu, x,
     &           ne, nlocJ, locJ, indJ,
     &           fCon, fObj, gCon, gObj, y,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
            if (iExit .gt. 0) return
         end if

         if (iExit .eq. 0) then
*           ------------------------------------------------------------
*           The functions have been computed at x.
*           ------------------------------------------------------------
            if (nnObj .gt. 0) then
               call dcopy ( nnObj,         gObj1, 1, dg, 1 )
               call daxpy ( nnObj, (-one), gObj , 1, dg, 1 )
               if (minimz .lt. 0) call dscal ( nnObj, sgnObj, dg, 1 )
            end if

            nzero = nnL - nnObj
            if (nzero .gt. 0) call dload ( nzero, zero, dg(nnObj+1), 1 )

            if (nnCon  .gt. 0) then
               call s8Gprd
     &            ( Transp, eps0,
     &              ne, nlocJ, locJ, indJ, negCon, nlocG, locG, gCon1,
     &              (-one), Lmul1, nnCon, one, dg, nnJac )
               call s8Gprd
     &            ( Transp, eps0,
     &              ne, nlocJ, locJ, indJ, negCon, nlocG, locG, gCon,
     &              one , Lmul1, nnCon, one, dg, nnJac )
            end if
            ydx   = ddot ( nnL, dg, 1, dx, 1 )
         end if
      end if

      if (dxHdx .lt. eps  .or.  iExit .ne. 0) then
         call dcopy
     &      ( nnL, tdx, 1, dx, 1 )
         call s8xHx
     &      ( nnL, dx, Udx, Hdx, dxHdx, iw, leniw, rw, lenrw )
      end if

      end ! subroutine s8x1
