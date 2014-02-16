*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     file  sn60srch.f
*
*     s6fdG    s6fdG1   s6srch   s6tols   s6usrf
*     lsrchc   lsrchq
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6fdG
     &   ( iExit, n, negCon,
     &     nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj,
     &     fgwrap, fgcon, fgobj,
     &     bl, bu, x,
     &     ne, nlocJ, locJ, indJ,
     &     fCon, fObj, gCon, gObj, y,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     fgwrap, fgcon, fgobj
      integer
     &     iExit, lencu, leniu, lenru, lencw, leniw, lenrw, n, ne,
     &     negCon, nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj, nlocJ,
     &     locJ(nlocJ), indJ(ne), iu(leniu), iw(leniw)
      double precision
     &     fObj, bl(n), bu(n), fCon(nnCon0), gObj(nnObj0), gCon(negCon),
     &     x(n), y(nnCon0), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     s6fdG   computes any missing objective and constraint gradients
*     for the current value of the variables in  x.
*
*     NOTE --  s6dcon overwrites the first  nnCon  elements of  y
*     if central differences are needed.
*
*     30 Dec 1991: First version based on Minos 5.4 routine m6grd.
*     17 Jul 1997: First thread-safe version.
*     11 Oct 1998: s6dcon and s6dobj merged.
*     24 Oct 2000: Updated for SNOPT 6.1
*     12 Oct 2003: snEXIT and snPRNT added.
*     14 Oct 2003: Enforced feasible perturbation.
*     01 Apr 2005: Current version of s6fdG.
*     ==================================================================
      integer
     &     gotFD, llocG, lfCon2, lgConu, lgObju, nlocG
*     ------------------------------------------------------------------
      gotFD     = iw(183) ! > 0 => some differences needed
      llocG     = iw(260) ! locG(nlocG) = column pointers for indG
      lfCon2    = iw(318) ! fCon2(nnCon) work vector
      lgConu    = iw(319) ! record of unknown derivatives and constants
      lgObju    = iw(323) ! record of unknown derivatives

      iExit     = 0
      nlocG     = nnJac + 1

      if (gotFD .gt. 0) then
         call s6fdG1
     &      ( iExit, n,
     &        nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj,
     &        fgwrap, fgcon, fgobj, bl, bu, x,
     &        ne, nlocJ, locJ, indJ, negCon, nlocG, iw(llocG),
     &        fObj, gObj, fCon, gCon,
     &        rw(lgObju), rw(lfCon2), rw(lgConu), y,
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )
         if (iExit .eq. 63) then ! The user didn't like some x's
            call snPRNT(  3,
     &         ' XXX  Unable to apply reversion when differencing',
     &         iw, leniw )
         end if
      end if

      end ! subroutine s6fdG

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6fdG1
     &   ( iExit, n,
     &     nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj,
     &     fgwrap, fgcon, fgobj, bl, bu, x,
     &     ne, nlocJ, locJ, indJ, negCon, nlocG, locG,
     &     fObj, gObj, fCon, gCon,
     &     gObju, fConu, gConu, y,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     fgwrap, fgcon, fgobj
      integer
     &     iExit, lencu, leniu, lenru, lencw, leniw, lenrw, n, negCon,
     &     nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj, ne, nlocG, nlocJ,
     &     indJ(ne), locG(nlocG), locJ(nlocJ), iu(leniu), iw(leniw)
      double precision
     &     fObj, bl(n), bu(n), fCon(nnCon0), fConu(nnCon0),
     &     gObj(nnObj0), gObju(nnObj0), gCon(negCon), gConu(negCon),
     &     x(n), y(nnCon0), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     s6fdG1  estimates missing elements in the objective gradient and
*     the columns of the Jacobian using finite differences of the
*     problem functions fObj and fCon.
*
*     The arrays y, fConu, gConu  and gObju are used as workspace.
*     Dummy elements of gObju and gConu define the unknown derivatives.
*
*     11 Oct 1998: First version based on combining s6dobj and s6dcon.
*     24 Oct 2000: Updated for SNOPT 6.1
*     12 Oct 2003: snEXIT and snPRNT added.
*     14 Oct 2003: Implemented feasible perturbation.
*     01 Apr 2005: Current version of s6fdG1.
*     ==================================================================
      integer
     &     inform, ir, j, k, kmax, l, lvlDer, lvlDif, modefg, nfCon,
     &     nfObj, numf, Status
      double precision
     &     infBnd, buj, delta, dxj, Fback, Fforwd, fdint(2), gdummy,
     &     tolx, xj
      logical
     &     Centrl, done, found, needJ , needG, someG, someJ
*     ------------------------------------------------------------------
      logical            yes   ,         no
      parameter         (yes   = .true., no   = .false.)
      double precision   one
      parameter         (one   = 1.0d+0)
      double precision   three,          four
      parameter         (three = 3.0d+0, four = 4.0d+0)
      parameter         (nfCon = 189) ! number of calls of fCon
      parameter         (nfObj = 194) ! number of calls of fObj
*     ------------------------------------------------------------------
      lvlDer    = iw( 70) ! = 0, 1, 2 or 3, the derivative level
      lvlDif    = iw(181) !    =1 (2) for forwd (cntrl) diffs

      tolx      = rw( 56) ! Minor feasibility tolerance
      gdummy    = rw( 69) ! definition of an 'unset' value
      infBnd    = rw( 70) ! definition of an infinite bound
      fdint(1)  = rw( 76) ! (1) forwrd diff. interval
      fdint(2)  = rw( 77) ! (2) cntrl  diff. interval

      iExit  = 0

*     The problem functions are called to provide functions only.
*     Status indicates that there is nothing special about the call.

      modefg = 0
      Status = 0

      someG  = lvlDer .eq. 0  .or.  lvlDer .eq. 2
      someJ  = lvlDer .eq. 0  .or.  lvlDer .eq. 1

      Centrl = lvlDif .eq. 2
      delta  = fdint(lvlDif)

      numf   = 0

      do j = 1, nnL

*        Look for the first missing element in this column.

         found  = no

         needJ  = j .le. nnJac  .and.  someJ
         needG  = j .le. nnObj  .and.  someG

         if ( needG ) then
            if (gObju(j) .eq. gdummy) found = yes
         end if

         if ( needJ ) then
            l      = locG(j)
            k      = locJ(j)
            kmax   = locJ(j+1) - 1
            done   = no

*+          ------------------------------------------------------------
*+          while (k .le. kmax  .and. .not.(found  .or.  done)) do
  120       if    (k .le. kmax  .and. .not.(found  .or.  done)) then
               ir = indJ(k)
               if (ir .gt. nnCon) then
                  done = yes
               else
                  if (gConu(l) .eq. gdummy) found = yes
                  l      = l + 1
               end if
               k  = k + 1
               go to 120
*+          end while
*+          ------------------------------------------------------------
            end if
         end if

         if ( found ) then
*           ------------------------------------------------------------
*           Some missing derivatives for this variable.
*           A finite difference is needed.
*           ------------------------------------------------------------
            xj     =  x(j)
            dxj    = delta*(one + abs(  xj ))
            buj    = bu(j)

            if (buj .lt. infBnd  .and.  xj+dxj+dxj .gt. buj+tolx) then
               dxj = -dxj
            end if

            x(j)   = xj   + dxj
            numf   = numf + 1

            call fgwrap
     &         ( inform, modefg, Status, needJ, needG,
     &           n, negCon, nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj,
     &           fgcon, fgobj, x,
     &           ne, nlocJ, locJ, indJ,
     &           fConu, Fforwd, gConu, gObju,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
            if (inform .ne. 0) then
               if (inform .gt. 0) then
                  iExit = inform
               else
                  iExit = 63   ! unable to move into undefined region
               end if
               go to 999
            end if

            if ( Centrl ) then
               dxj  =  dxj + dxj
               x(j) =   xj + dxj
               numf = numf + 1

               call fgwrap
     &            ( inform, modefg, Status, needJ, needG,
     &              n, negCon, nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj,
     &              fgcon, fgobj, x,
     &              ne, nlocJ, locJ, indJ,
     &              y, Fback, gConu, gObju,
     &              cu, lencu, iu, leniu, ru, lenru,
     &              cw, lencw, iw, leniw, rw, lenrw )
               if (inform .ne. 0) then
                  if (inform .gt. 0) then
                     iExit = inform
                  else
                     iExit = 63 ! unable to move into undefined region
                  end if
                  go to 999
               end if
            end if

            if ( needG ) then
               if (gObju(j) .eq. gdummy) then
                  if (Centrl) then
                     gObj(j) = (four*Fforwd - three*fObj - Fback)
     &                                         / dxj
                  else
                     gObj(j) = (Fforwd - fObj) / dxj
                  end if
               end if
            end if

            if ( needJ ) then
               l      = locG(j)
               k      = locJ(j)
               done   = no

*+             ---------------------------------------------------------
*+             while (k .le. kmax  .and.  .not. done) do
  140          if    (k .le. kmax  .and.  .not. done) then
                  ir = indJ(k)
                  if (ir .gt. nnCon) then
                     done = yes
                  else
                     if (gConu(l) .eq. gdummy) then
                        if (Centrl) then
                           gCon(l) = (four*fConu(ir) - three*fCon(ir)
     &                                               - y(ir))/dxj
                        else
                           gCon(l) = (fConu(ir) - fCon(ir))/dxj
                        end if
                     end if
                     l   = l + 1
                  end if

                  k  = k + 1
                  go to 140
*+             end while
*+             ---------------------------------------------------------
               end if
            end if ! j <= nnJac

            x(j)   = xj
         end if ! found
      end do

*     ------------------------------------------------------------------
*     The missing derivatives have been estimated.
*     Finish up with some housekeeping.
*     ------------------------------------------------------------------
  900 l           = lvlDif      + 1
      iw(nfCon+l) = iw(nfCon+l) + numf
      iw(nfObj+l) = iw(nfObj+l) + numf

  999 return

      end ! subroutine s6fdG1

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6srch
     &   ( iExit, fgwrap, fgcon, fgobj,
     &     debug, Elastc, fonly, prFeas, iObj, sclObj,
     &     n, nb, nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj,
     &     itn, wolfeG, sgnObj, step, stepmn, stepmx, pNorm, xNorm,
     &     fMrt, fMrt1, gMrt, gMrt1, sInf, sInf1, sInf2, wtInf,
     &     ne, nlocJ, locJ, indJ, Jcol, negCon, nlocG, locG,
     &     fObj1, fCon1, gCon1, gObj1, fObj2, fCon2, gCon2, gObj2,
     &     dx, dyCon, x, x1, x2, yCon, yCon1, yCon2, xPen,
     &     y, y1, y2, cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     fgwrap, fgcon, fgobj
      logical
     &     debug, Elastc, fonly, prFeas
      integer
     &     iExit, iObj, itn, lencu, lencw, leniu, leniw, lenru, lenrw,
     &     n, nb, ne, negCon, nlocG, nlocJ, nnCon0, nnCon, nnJac,
     &     nnL, nnObj0, nnObj, indJ(ne), locJ(nlocJ), locG(nlocG),
     &     iu(leniu), iw(leniw)
      double precision
     &     wolfeG, fMrt, fMrt1, gMrt, gMrt1, pNorm, sInf, sInf1, sInf2,
     &     sclObj, sgnObj, step, stepmx, wtInf, xNorm, Jcol(ne),
     &     fObj1, fCon1(nnCon0), gCon1(negCon), gObj1(nnObj0),
     &     fObj2, fCon2(nnCon0), gCon2(negCon), gObj2(nnObj0),
     &     yCon(nnCon0), yCon1(nnCon0), yCon2(nnCon0),
     &     dx(nb), dyCon(nnCon0), x(nb), x1(nb), x2(nb),
     &     xPen(nnCon0), y(nb), y1(nb), y2(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     s6srch  finds a step along the search direction  p,  such that
*     the function  fMrt  is sufficiently reduced, i.e.,
*               fMrt(x + step*p)  <  fMrt(x).
*
*     On entry,
*     step     Initial estimate of the step.
*     stepmx   Maximum allowable step.
*     wolfeG   Line search accuracy parameter in the range (0,1).
*              0.001 means very accurate.   0.99 means quite sloppy.
*     fonly    true if function-only search should be used.
*     gMrt1    The directional derivative of the merit function.
*     fMrt1    Current value of fMrt(x).
*     x(nb)    The base point for the search.
*     x2(nb)   x2 = x + dx,  the QP solution.
*     yCon     The base point for the multipliers.
*     p(nb)    The search direction.
*     yCon2    the QP multipliers.
*     dx(nb)   The search direction, x2 - x1.
*
*     Output parameters...
*
*     step     The final step.
*     fMrt1    The final value of fMrt.
*     x1       Final value of the variables.
*     fCon1,gCon1,gObj1,yCon1 are defined at the new point x1.
*
*     fCon2,gCon2,gObj2,yCon2 are work arrays.
*
*     iExit    Result
*     -----    ------
*      >0      Fatal error
*       0      Repeat the search with smaller stpmax.
*      -1      The search is successful and step < stpmax.
*      -2      The search is successful and step = stpmax.
*      -3      A better point was found but too many functions
*              were needed (not sufficient decrease).
*      -4      stpmax < tolabs (too small to do a search).
*      -5      step   < stepmn (lsrchq only -- maybe want to switch
*              to central differences to get a better direction).
*      -6      No useful step.
*              The interval of uncertainty is less than 2*tolabs.
*              The minimizer is very close to step = zero
*              or the gradients are not sufficiently accurate.
*      -7      Too many function calls.
*      -8      Bad input parameters
*              (stpmax le toltny  or  oldg ge 0).
*
*     30 Dec 1991: First version based on NPSOL 4.6 routine npsrch.
*     28 Sep 1993: Allow functions to say "undefined at this point".
*                  (Back up and try again.)
*     18 Feb 1994: Back up in a similar way if vimax increases a lot.
*                  Deleted after first limited-memory version.
*     29 Dec 1994: Merit function calculations included explicitly.
*     06 Apr 1996: Special coding for the unit step.
*                  On entry, x2 is the QP solution.
*     18 Oct 1996: First Min-sum version.
*     17 Jul 1997: First thread-safe version.
*     11 Oct 1998: Facility to combine funobj and funcon added.
*     11 Jun 2000: Tolerances computed in a subroutine.
*     24 Oct 2000: Updated for SNOPT 6.1
*     04 Aug 2003: snEXIT and snPRNT adopted.
*     01 Apr 2005: Current version of s6srch.
*     ==================================================================
      character
     &     str*80
      external
     &     ddot
      logical
     &     Gknown, done, first, imprvd, nlnCon, nlnObj,
     &     braktd, crampd, extrap, moved, vset, wset
      integer
     &     inform, iPrint, jObj, maxf, modefg, nLin, nnJ1,
     &     nout, nsamea, nsameb, numf, Status
      double precision
     &     eps0, bigFx, dsInf, epsaf, fMrt2, g0, gMrt2, oldf, oldg,
     &     tolabs, tolrel, toltny, sbest, fbest, gbest, stpmax, targtg,
     &     aLow, bUpp, fa, factor, ftry, fv, fw, gtry, gw, stepmn,
     &     tolmax, xtry, xv, xw, ddot
*     ------------------------------------------------------------------
      integer            Normal
      parameter         (Normal = 0)
      double precision   zero,          half,          one
      parameter         (zero = 0.0d+0, half = 0.5d+0, one = 1.0d+0)
      double precision   wolfeF
      parameter         (wolfeF = 1.0d-4)
*     ------------------------------------------------------------------
      iPrint    = iw( 12) ! Print file

      eps0      = rw(  2) ! eps**(4/5)
      bigFx     = rw( 71) ! unbounded objective.

      iExit     = 0

*     ------------------------------------------------------------------
*     Set the input parameters for lsrchc or lsrchq.
*
*     stepmn  is used by lsrchq.  If  step  would be less than  stepmn,
*             the search will be terminated early.
*             If p was found using forward or backward differences,
*             stepmn  should be positive (and related to the difference
*             interval used).
*             If p was found using central differences (lvlDif = 2)
*             stepmn  should be zero.
*
*     epsaf   is the absolute function precision. If f(x1) and f(x2) are
*             as close as  epsaf,  we cannot safely conclude from the
*             function values alone which of x1 or x2 is a better point.
*
*     tolabs  is an estimate of the absolute spacing between points
*             along  p.  This step should produce a perturbation
*             of  epsaf  in the merit function.
*
*     tolrel  is an estimate of the relative spacing between points
*             along  p.
*
*     toltny  is the minimum allowable absolute spacing between points
*             along  p.
*     ------------------------------------------------------------------
      stpmax = stepmx
      Gknown = .not. fonly

      if ( fonly ) then
         maxf   = 15
         modefg =  0
      else
         maxf   = 10
         modefg =  2
      end if

      nlin   = n - nnJac
      nnJ1   = nnJac + 1

      nout   = iPrint
      jObj   = n + iObj
      Status = 0
      nlnCon = nnCon .gt. 0
      nlnObj = nnObj .gt. 0

*     Define the line search tolerances.

      call s6tols
     &   ( nb, epsaf, stepmx, tolabs, tolrel, toltny,
     &     pNorm, xNorm, fMrt, dx, x, rw, lenrw )

      oldf   = fMrt
      oldg   = gMrt
      fMrt1  = fMrt
      gMrt1  = gMrt
      call dcopy ( nb   ,    x, 1,    x1, 1 )
      if ( nlnCon )
     &call dcopy ( nnCon, yCon, 1, yCon1, 1 )
      if ( Elastc ) then
         sInf1  = sInf
         dsInf  = sInf2 - sInf1
      end if

      fObj2  = zero             ! keeps ftnchek quiet
      fMrt2  = zero             ! keeps ftnchek quiet
      ftry   = zero             ! keeps ftnchek quiet
      gtry   = zero             ! keeps ftnchek quiet
      fv     = zero             ! keeps ftnchek quiet
      fw     = zero             ! keeps ftnchek quiet

      first  = .true.
      sbest  = zero
      fbest  = zero
      gbest  = (one     - wolfeF)*oldg
      targtg = (wolfeF  - wolfeG)*oldg
      g0     = gbest

      if (debug) write(nout, 1000) itn, pNorm

*     ------------------------------------------------------------------
*     Commence main loop, entering lsrchc or lsrchq two or more times.
*     first = true for the first entry,  false for subsequent entries.
*     done  = true indicates termination, in which case inform gives
*     the result of the search (with inform = iExit as above).
*     ------------------------------------------------------------------
*+    repeat
  200    if ( Gknown ) then
            call lsrchc
     &         ( inform, first , debug , done  , imprvd,
     &           maxf  , numf  , nout  ,
     &           stpmax,         epsaf ,
     &           g0    , targtg, ftry  , gtry  ,
     &           tolabs, tolrel, toltny,
     &           step  , sbest , fbest , gbest ,
     &           braktd, crampd, extrap, moved , wset  ,
     &           nsamea, nsameb,
     &           aLow  , bUpp  , factor,
     &           xtry  , xw    , fw    , gw    , tolmax )
         else
            call lsrchq
     &         ( inform, first , debug , done  , imprvd,
     &           maxf  , numf  , nout  ,
     &           stpmax, stepmn, epsaf ,
     &           g0    , targtg, ftry  ,
     &           tolabs, tolrel, toltny,
     &           step  , sbest , fbest ,
     &           braktd, crampd, extrap, moved , vset  , wset  ,
     &           nsamea, nsameb,
     &           aLow  , bUpp  , fa    , factor,
     &           xtry  , xw    , fw    , xv    , fv    , tolmax)
         end if

         if ( imprvd ) then
            fMrt1 = fMrt2

            if ( nlnCon ) then
               call dcopy ( nnCon, fCon2, 1, fCon1, 1 )
               if ( Gknown ) call dcopy ( negCon, gCon2, 1, gCon1, 1 )
            end if

            if ( nlnObj ) then
               fObj1 = fObj2
               if ( Gknown ) call dcopy ( nnObj , gObj2, 1, gObj1, 1 )

*              Terminate if the objective is unbounded below in the
*              feasible region.

               if (sgnObj*fObj1 .lt. -bigFx  .and.  prFeas) then
                  iExit = 21
                  go to 900
               end if
            end if
         end if ! imprvd

*        ---------------------------------------------------------------
*           done = false  first time through.
*        If done = false, the functions must be computed for the next
*                  entry to lsrchc or lsrchq.
*        If done = true,  this is the last time through and inform ge 1.
*        ---------------------------------------------------------------
         if (.not. done) then
            if (step .ne. one) then
               call dcopy ( nb,       x1, 1, x2, 1 )
               call daxpy ( nb, step, dx, 1, x2, 1 )
               if ( Elastc ) sInf2 = sInf1 + step*dsInf
            end if

            if (nnL .gt. 0) then
               call fgwrap
     &            ( inform, modefg, Status, nlnCon, nlnObj,
     &              n, negCon, nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj,
     &              fgcon, fgobj, x2,
     &              ne, nlocJ, locJ, indJ,
     &              fCon2, fObj2, gCon2, gObj2,
     &              cu, lencu, iu, leniu, ru, lenru,
     &              cw, lencw, iw, leniw, rw, lenrw )
               if (inform .ne. 0) then
                  if (inform .gt. 0) then
                     iExit  = inform
                  else
                     inform = 0 ! Redo the search
                  end if
                  go to 900
               end if
            end if

            if (iObj .eq. 0) then
               fMrt2 = zero
            else
               fMrt2 = sgnObj*x2(jObj)*sclObj
            end if

            if ( nlnObj ) then
               fMrt2 = fMrt2 + sgnObj*fObj2
            end if

            if ( nlnCon ) then
*              ---------------------------------------------------------
*              Compute w and y, the constraint violations and the
*              negative of the gradient of the merit function with
*              respect to the nonlinear slacks. These quantities define
*              the directional derivative of the merit function.
*              Finally, add the nonlinear terms to the merit function.
*              ---------------------------------------------------------
               if (step .le. one) then
                  call dcopy ( nnCon,       yCon1, 1, yCon2, 1 )
                  call daxpy ( nnCon, step, dyCon, 1, yCon2, 1 )
               end if

*              Compute the constraint violations and aux. multipliers:
*              Fv = fCon + A(linear) x - nonlinear slacks.
*              y  = yCon2 - xPen*viol.

               call dcopy ( nnCon,           fCon2, 1, y1, 1 )
               call daxpy ( nnCon, (-one), x2(n+1), 1, y1, 1 )
               if (nlin .gt. 0) then
                  call s2Aprd
     &               ( Normal, eps0,
     &                 ne, nlin+1, locJ(nnJ1), indJ, Jcol,
     &                 one, x2(nnJ1), nlin, one, y1, nnCon )
               end if

               call dcopy ( nnCon,   y1, 1, y, 1 )
               call ddscl ( nnCon, xPen, 1, y, 1 )
               fMrt2 = fMrt2 -      ddot( nnCon, yCon2, 1, y1, 1 )
     &                       + half*ddot( nnCon, y    , 1, y1, 1 )

*              If we are in elastic mode, include the contribution from
*              the violations of the elastic variables.

               if ( Elastc ) fMrt2 = fMrt2 + wtInf*sInf2

               if ( Gknown ) then
                  call daxpy ( nnCon, (-one), yCon2, 1, y, 1 )
                  call dscal ( nnCon, (-one), y    , 1 )
               end if
            end if ! nlnCon

            ftry  = fMrt2 - oldf - wolfeF*oldg*step

            if ( Gknown ) then
*              ---------------------------------------------------------
*              A gradient search is requested.
*              Compute the directional derivative gtry.
*              ---------------------------------------------------------
               if (iObj .eq. 0) then
                  gMrt2 = zero
               else
                  gMrt2 = sgnObj*dx(jObj)*sclObj
               end if

               if ( nlnObj ) then
                  gMrt2 = gMrt2 + sgnObj*ddot( nnObj, gObj2, 1, dx, 1 )
               end if

               if ( nlnCon ) then
*                 ------------------------------------------------------
*                 Form  Jp (including linear columns).
*                 Set  y2 = Jp + A(linear) p - p(slacks).
*                 ------------------------------------------------------
                  call s8Gprd
     &               ( Normal, eps0,
     &                 ne, nlocJ, locJ, indJ, negCon, nlocG, locG,gCon2,
     &                 one, dx, nnJac, zero, y2, nnCon )

                  if (nlin .gt. 0) then
                     call s2Aprd
     &                  ( Normal, eps0,
     &                    ne, nlin+1, locJ(nnJ1), indJ, Jcol,
     &                    one, dx(nnJ1), nlin, one, y2, nnCon )
                  end if
                  call daxpy ( nnCon, (-one), dx(n+1), 1, y2, 1 )

                  gMrt2 = gMrt2 - ddot  ( nnCon, y , 1, y2   , 1 )
     &                          - ddot  ( nnCon, y1, 1, dyCon, 1 )

*                 If we are elastic mode, include the contribution from
*                 the violations of the elastic variables.

                  if ( Elastc ) gMrt2 = gMrt2 + wtInf*dsInf

               end if ! nlnCon

               gtry   = gMrt2 - wolfeF*oldg

            end if ! Gknown
         end if ! not done
*+    until (      done)
      if    (.not. done) go to 200

*     ==================================================================
*     The search is done.
*     Finish with  x1 = the best point found so far.
*     ==================================================================
      step = sbest

      if (imprvd) then
         call dcopy ( nb   ,    x2, 1,    x1, 1 )
         if ( nlnCon )
     &   call dcopy ( nnCon, yCon2, 1, yCon1, 1 )
         if ( Elastc ) sInf1 = sInf2
      else if (step .gt. zero) then
         call daxpy ( nb   ,  step, dx   , 1, x1   , 1 )
         if ( nlnCon )
     &   call daxpy ( nnCon,  step, dyCon, 1, yCon1, 1 )
         if ( Elastc ) sInf1 = sInf1 + step*dsInf
      end if

*     ------------------------------------------------------------------
*     Print any warning messages.
*     ------------------------------------------------------------------
      if (inform .eq. 7) then
         write(str, 1700) numf
         call snPRNT( 23, str, iw, leniw )

      else if (inform .eq. 8) then
         if (oldg .ge. zero) then
            write(str, 1600) oldg
            call snPRNT( 21, str, iw, leniw )
         else
            write(str, 1800) stepmx, pNorm, oldg, numf
            call snPRNT( 21, str, iw, leniw )
         end if
      end if

      iExit = -inform

  900 return

 1000 format(// ' --------------------------------------------'
     &       /  ' Output from s6srch following iteration', i7,
     &      5x, ' Norm p =', 1p, e11.2 )
 1600 format(   ' XXX  The search direction is uphill.  gMrt  =',
     &          1p, e9.1)
 1700 format(   ' XXX  The line search has evaluated the functions', i5,
     &          '  times')
 1800 format(   ' stepmx =', 1p, e11.2, '    pNorm =', e11.2,
     &          '   gMrt =',     e11.2, '    numf  =', i3)

      end ! subroutine s6srch

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6tols
     &   ( nb, epsaf, stepmx, tolabs, tolrel, toltny,
     &     xdNorm, xNorm, f, dx, x, rw, lenrw )

      implicit
     &     none
      integer
     &     nb, lenrw
      double precision
     &     epsaf, f, stepmx, tolabs, tolrel, toltny, xdNorm, xNorm,
     &     dx(nb), x(nb), rw(lenrw)

*     ==================================================================
*     s6tols defines various tolerances for the line search.
*
*     11 Jun 2000: First   version of s6tols.
*     11 Jun 2000: Current version of s6tols.
*     ==================================================================
      integer
     &     j
      double precision
     &     eps, eps0, epsrf, s, q, t, tolax, tolrx
*     ------------------------------------------------------------------
      double precision   zero,          one
      parameter         (zero = 0.0d+0, one = 1.0d+0)
*     ------------------------------------------------------------------
      eps       = rw(  1) ! machine precision.  IEEE DP  2.22e-16
      eps0      = rw(  2) ! eps**(4/5)
      epsrf     = rw( 73) ! relative function precision.
*     ------------------------------------------------------------------
      epsaf  = max( epsrf, eps )*(one + abs(f))
      tolax  = eps0
      tolrx  = eps0
      t      = xNorm*tolrx  +  tolax
      if (t .lt. xdNorm*stepmx) then
         tolabs = t/xdNorm
      else
         tolabs = stepmx
      end if
      tolrel = max( tolrx, eps )

      t      = zero
      do j = 1, nb
         s     = abs(dx(j))
         q     = abs( x(j))*tolrx + tolax
         if (s .gt. q*t) t = s / q
      end do

      if (t*tolabs .gt. one) then
         toltny = one / t
      else
         toltny = tolabs
      end if

      end ! subroutine s6tols

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6usrf
     &   ( Status, n, x, xPert, damper, userfg,
     &     needF, nF  , FPert,
     &     needG, lenG, GPert,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      external
     &     userfg
      integer
     &     n, needF, needG, nF, lenG, lencu, leniu, lenru, Status,
     &     iu(leniu)
      double precision
     &     damper, FPert(nF), GPert(lenG), x(n), xPert(n), ru(lenru)
      character
     &     cu(lencu)*8

*     ==================================================================
*     s6usrf  attempts to compute the snoptA problem functions at a
*     point xPert.  If the functions are undefined at xPert, then the
*     evaluation takes place at a point on the ray joining xPert and
*     a point  x  at which the functions are known to be well-defined.
*
*     On entry:
*        xPert  is the point at which the functions are required.
*
*        x      is a point at which the problem functions have been
*               computed successfully.
*
*     On exit:
*        Status is nonnegative if the problem function were evaluated
*               successfully.   Otherwise, five attempts to evaluate
*               the functions at points closer to x were unsuccessful.
*
*        If Status ge 0,  then the output values of damper and xPert are
*        defined as follows:
*
*        damper (0 lt damper le 1) is 1 if the functions were
*               evaluated successfuly at the input value of xPert.
*               Otherwise the problem functions were evaluated
*               successfuly at xPert + damper*(xPert - x).
*
*        xPert  is  xPert(in) + damper*(xPert(in) - x),
*
*     26 Oct 2002: First version.
*     22 Apr 2007: damper added as an output argument.
*     ==================================================================
      integer
     &     j, tries
*     ------------------------------------------------------------------
      double precision   one,            ten
      parameter         (one   = 1.0d+0, ten  =10.0d+0)
*     ------------------------------------------------------------------
      damper = one
      tries  = 0
*     ==================================================================
*     Repeat                       (until problem functions are defined)

  100    tries  = tries + 1
         Status = 0
         call userfg
     &      ( Status, n, xPert,
     &        needF, nF, FPert,
     &        needG, lenG, GPert,
     &        cu, lencu, iu, leniu, ru, lenru )
         if (Status .eq. -1) then
            damper = damper*damper/ten
            do   j = 1, n
               xPert(j) = damper*xPert(j) + (one - damper)*x(j)
            end do
         end if

*+    until (.not. (Status .eq. -1  .and.  tries .lt. 5))
      if (          Status .eq. -1  .and.  tries .lt. 5 ) go to 100
*     ==================================================================

      end ! subroutine s6usrf

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsrchc
     &   ( iExit , first , debug , done  , imprvd,
     &     maxf  , numf  , nout  ,
     &     alfmax,         epsaf ,
     &     g0    , targtg, ftry  , gtry  ,
     &     tolabs, tolrel, toltny,
     &     alfa  , alfbst, fbest , gbest ,
     &     braktd, crampd, extrap, moved , wset  ,
     &     nsamea, nsameb,
     &     a     , b     , factor,
     &     xtry  , xw    , fw    , gw    , tolmax )

      implicit
     &     none
      logical
     &     first, debug, done, imprvd, braktd, crampd, extrap,
     &     moved, wset
      integer
     &     iExit, maxf, numf, nout, nsamea, nsameb
      double precision
     &     alfmax, epsaf, g0, targtg, ftry, gtry, tolabs, tolrel,
     &     toltny, alfa, alfbst, fbest, gbest, a, b, factor, xtry,
     &     xw, fw, gw, tolmax

*     ==================================================================
*     lsrchc  finds a sequence of improving estimates of a minimizer of
*     the univariate function f(alpha) in the interval (0,alfmax].
*     f(alpha) is a smooth function such that  f(0) = 0  and  f'(0) < 0.
*     lsrchc requires both  f(alpha)  and  f'(alpha)  to be evaluated at
*     points in the interval.  Estimates of the minimizer are computed
*     using safeguarded cubic interpolation.
*
*     Reverse communication is used to allow the calling program to
*     evaluate f and f'.  Some of the parameters must be set or tested
*     by the calling program.  The remainder would ordinarily be local
*     variables.
*
*     Input parameters (relevant to the calling program)
*     --------------------------------------------------
*
*     first         must be true on the first entry. It is subsequently
*                   altered by lsrchc.
*
*     debug         specifies whether detailed output is wanted.
*
*     maxf          is an upper limit on the number of times lsrchc is
*                   to be entered consecutively with done = false
*                   (following an initial entry with first = true).
*
*     alfa          is the first estimate of a minimizer.  alfa is
*                   subsequently altered by lsrchc (see below).
*
*     alfmax        is the upper limit of the interval to be searched.
*
*     epsaf         is an estimate of the absolute precision in the
*                   computed value of f(0).
*
*     ftry, gtry    are the values of f, f'  at the new point
*                   alfa = alfbst + xtry.
*
*     g0            is the value of f'(0).  g0 must be negative.
*
*     tolabs,tolrel define a function tol(alfa) = tolrel*alfa + tolabs
*                   such that if f has already been evaluated at alfa,
*                   it will not be evaluated closer than tol(alfa).
*                   These values may be reduced by lsrchc.
*
*     targtg        is the target value of abs(f'(alfa)). The search
*                   is terminated when
*                    abs(f'(alfa)) le targtg and f(alfa) lt 0.
*
*     toltny        is the smallest value that tolabs is allowed to be
*                   reduced to.
*
*     Output parameters (relevant to the calling program)
*     ---------------------------------------------------
*
*     imprvd        is true if the previous alfa was the best point so
*                   far.  Any related quantities should be saved by the
*                   calling program (e.g., gradient arrays) before
*                   paying attention to the variable done.
*
*     done = false  means the calling program should evaluate
*                      ftry = f(alfa),  gtry = f'(alfa)
*                   for the new trial alfa, and re-enter lsrchc.
*
*     done = true   means that no new alfa was calculated.  The value
*                   of iExit gives the result of the search as follows
*
*                   iExit = 1 means the search has terminated
*                             successfully with alfbst < alfmax.
*
*                   iExit = 2 means the search has terminated
*                             successfully with alfbst = alfmax.
*
*                   iExit = 3 means that the search failed to find a
*                             point of sufficient decrease.
*                             The function is either decreasing at
*                             alfmax or maxf function evaluations
*                             have been exceeded.
*
*                   iExit = 4 means alfmax is so small that a search
*                             should not have been attempted.
*
*                   iExit = 5 is never set by lsrchc.
*
*                   iExit = 6 means the search has failed to find a
*                             useful step.  The interval of uncertainty
*                             is [0,b] with b < 2*tolabs. A minimizer
*                             lies very close to alfa = 0, or f'(0) is
*                             not sufficiently accurate.
*
*                   iExit = 7 if no better point could be found after
*                             maxf  function calls.
*
*                   iExit = 8 means the input parameters were bad.
*                             alfmax le toltny  or g0 ge zero.
*                             No function evaluations were made.
*
*     numf          counts the number of times lsrchc has been entered
*                   consecutively with done = false (i.e., with a new
*                   function value ftry).
*
*     alfa          is the point at which the next function ftry and
*                   derivative gtry must be computed.
*
*     alfbst        should be accepted by the calling program as the
*                   approximate minimizer, whenever lsrchc returns
*                   iExit = 1 or 2 (and possibly 3).
*
*     fbest, gbest  will be the corresponding values of f, f'.
*
*
*     The following parameters retain information between entries
*     -----------------------------------------------------------
*
*     braktd        is false if f and f' have not been evaluated at
*                   the far end of the interval of uncertainty.  In this
*                   case, the point b will be at alfmax + tol(alfmax).
*
*     crampd        is true if alfmax is very small (le tolabs).  If the
*                   search fails, this indicates that a zero step should
*                   be taken.
*
*     extrap        is true if xw lies outside the interval of
*                   uncertainty.  In this case, extra safeguards are
*                   applied to allow for instability in the polynomial
*                   fit.
*
*     moved         is true if a better point has been found, i.e.,
*                   alfbst gt 0.
*
*     wset          records whether a second-best point has been
*                   determined it will always be true when convergence
*                   is tested.
*
*     nsamea        is the number of consecutive times that the
*                   left-hand end point of the interval of uncertainty
*                   has remained the same.
*
*     nsameb        similarly for the right-hand end.
*
*     a, b, alfbst  define the current interval of uncertainty.
*                   A minimizer lies somewhere in the interval
*                   [alfbst + a, alfbst + b].
*
*     alfbst        is the best point so far.  It is always at one end
*                   of the interval of uncertainty.  hence we have
*                   either  a lt 0,  b = 0  or  a = 0,  b gt 0.
*
*     fbest, gbest  are the values of f, f' at the point alfbst.
*
*     factor        controls the rate at which extrapolated estimates
*                   of alfa may expand into the interval of uncertainty.
*                   factor is not used if a minimizer has been bracketed
*                   (i.e., when the variable braktd is true).
*
*     fw, gw        are the values of f, f' at the point alfbst + xw.
*                   they are not defined until wset is true.
*
*     xtry          is the trial point in the shifted interval (a, b).
*
*     xw            is such that  alfbst + xw  is the second-best point.
*                   it is not defined until  wset  is true.
*                   in some cases,  xw  will replace a previous  xw
*                   that has a lower function but has just been excluded
*                   from the interval of uncertainty.
*
*
*     Systems Optimization Laboratory, Stanford University, California.
*     Original version February 1982.  Rev. May 1983.
*     Original f77 version 22-August-1985.
*     14 Sep 1992: Introduced quitI, quitF, etc.
*     22 Nov 1995: Altered criterion for reducing the step below tolabs.
*     17 Jul 1997: Removed saved variables for thread-safe version.
*     19 Apr 2000: QuitF only allowed after a move.
*     19 Apr 2000: Current version.
*     ==================================================================
      logical
     &     badfun, closef, found, quitF, quitI, fitok, setxw
      double precision
     &     absr, artifa, artifb, daux, dtry,
     &     q, r, s, scale, tol, truea, trueb, xmidpt
*     ------------------------------------------------------------------
      double precision   zero,          point1,          half
      parameter         (zero = 0.0d+0, point1 = 0.1d+0, half = 0.5d+0)
      double precision   one,           three,           five
      parameter         (one  = 1.0d+0, three  = 3.0d+0, five = 5.0d+0)
      double precision   ten,           eleven
      parameter         (ten  = 1.0d+1, eleven = 1.1d+1               )
*     ------------------------------------------------------------------
*     Local variables
*     ===============
*
*     closef     is true if the new function ftry is within epsaf of
*                fbest (up or down).
*
*     found      is true if the sufficient decrease conditions hold at
*                alfbst.
*
*     quitF      is true when  maxf  function calls have been made.
*
*     quitI      is true when the interval of uncertainty is less than
*                2*tol.
*  ---------------------------------------------------------------------

      badfun = .false.
      quitF  = .false.
      quitI  = .false.
      imprvd = .false.

      if (first) then
*        ---------------------------------------------------------------
*        First entry.  Initialize various quantities, check input data
*        and prepare to evaluate the function at the initial alfa.
*        ---------------------------------------------------------------
         first  = .false.
         numf   = 0
         alfbst = zero
         badfun = alfmax .le. toltny  .or.  g0 .ge. zero
         done   = badfun
         moved  = .false.

         if (.not. done) then
            braktd = .false.
            crampd = alfmax .le. tolabs
            extrap = .false.
            wset   = .false.
            nsamea = 0
            nsameb = 0

            tolmax = tolabs + tolrel*alfmax
            a      = zero
            b      = alfmax + tolmax
            factor = five
            tol    = tolabs
            xtry   = alfa
            if (debug) then
               write(nout, 1000)
     &              g0    , tolabs, alfmax,
     &              targtg, tolrel, epsaf , crampd
            end if
         end if
      else
*        ---------------------------------------------------------------
*        Subsequent entries. The function has just been evaluated at
*        alfa = alfbst + xtry,  giving ftry and gtry.
*        ---------------------------------------------------------------
         if (debug) write(nout, 1100) alfa, ftry, gtry

         numf   = numf   + 1
         nsamea = nsamea + 1
         nsameb = nsameb + 1

         if (.not. braktd) then
            tolmax = tolabs + tolrel*alfmax
            b      = alfmax - alfbst + tolmax
         end if

*        See if the new step is better.  If alfa is large enough that
*        ftry can be distinguished numerically from zero,  the function
*        is required to be sufficiently negative.

         closef = abs( ftry - fbest ) .le.  epsaf
         if (closef) then
            imprvd =  abs( gtry ) .le. abs( gbest )
         else
            imprvd = ftry .lt. fbest
         end if

         if (imprvd) then

*           We seem to have an improvement.  The new point becomes the
*           origin and other points are shifted accordingly.

            fw     = fbest
            fbest  = ftry
            gw     = gbest
            gbest  = gtry
            alfbst = alfa
            moved  = .true.

            a      = a    - xtry
            b      = b    - xtry
            xw     = zero - xtry
            wset   = .true.
            extrap =       xw .lt. zero  .and.  gbest .lt. zero
     &               .or.  xw .gt. zero  .and.  gbest .gt. zero

*           Decrease the length of the interval of uncertainty.

            if (gtry .le. zero) then
               a      = zero
               nsamea = 0
            else
               b      = zero
               nsameb = 0
               braktd = .true.
            end if
         else

*           The new function value is not better than the best point so
*           far.  The origin remains unchanged but the new point may
*           qualify as xw.  xtry must be a new bound on the best point.

            if (xtry .le. zero) then
               a      = xtry
               nsamea = 0
            else
               b      = xtry
               nsameb = 0
               braktd = .true.
            end if

*           If xw has not been set or ftry is better than fw, update the
*           points accordingly.

            if (wset) then
               setxw = ftry .lt. fw  .or.  .not. extrap
            else
               setxw = .true.
            end if

            if (setxw) then
               xw     = xtry
               fw     = ftry
               gw     = gtry
               wset   = .true.
               extrap = .false.
            end if
         end if

*        ---------------------------------------------------------------
*        Check the termination criteria.  wset will always be true.
*        ---------------------------------------------------------------
         tol    = tolabs + tolrel*alfbst
         truea  = alfbst + a
         trueb  = alfbst + b

         found  = abs(gbest) .le. targtg
         quitF  = numf       .ge. maxf    .and.  moved
         quitI  = b - a      .le. tol + tol

         if (quitI  .and. .not. moved) then

*           The interval of uncertainty appears to be small enough,
*           but no better point has been found.  Check that changing
*           alfa by b-a changes f by less than epsaf.

            tol    = tol/ten
            tolabs = tol
            quitI  =     tol .le. toltny  .or.
     &               abs(fw) .le. epsaf   .and.  gw .le. epsaf
         end if

         done  = quitF  .or.  quitI  .or.  found

         if (debug) then
            write(nout, 1200)
     &           truea    , trueb , b - a , tol   ,
     &           nsamea   , nsameb, numf  ,
     &           braktd   , extrap, closef, imprvd,
     &           found    , quitI ,
     &           alfbst   , fbest , gbest ,
     &           alfbst+xw, fw    , gw
         end if

*        ---------------------------------------------------------------
*        Proceed with the computation of a trial steplength.
*        The choices are...
*        1. Parabolic fit using derivatives only, if the f values are
*           close.
*        2. Cubic fit for a minimizer, using both f and f'.
*        3. Damped cubic or parabolic fit if the regular fit appears to
*           be consistently overestimating the distance to a minimizer.
*        4. Bisection, geometric bisection, or a step of  tol  if
*           choices 2 or 3 are unsatisfactory.
*        ---------------------------------------------------------------
         if (.not. done) then
            xmidpt = half*(a + b)
            s      = zero
            q      = zero

            if (closef) then
*              ---------------------------------------------------------
*              Fit a parabola to the two best gradient values.
*              ---------------------------------------------------------
               s      = gbest
               q      = gbest - gw
               if (debug) write(nout, 2200)
            else
*              ---------------------------------------------------------
*              Fit cubic through  fbest  and  fw.
*              ---------------------------------------------------------
               if (debug) write(nout, 2100)
               fitok  = .true.
               r      = three*(fbest - fw)/xw + gbest + gw
               absr   = abs( r )
               s      = sqrt( abs( gbest ) ) * sqrt( abs( gw ) )

*              Compute  q =  the square root of  r*r - gbest*gw.
*              The method avoids unnecessary underflow and overflow.

               if ((gw .lt. zero  .and.  gbest .gt. zero) .or.
     &             (gw .gt. zero  .and.  gbest .lt. zero)) then
                  scale  = absr + s
                  if (scale .eq. zero) then
                     q  = zero
                  else
                     q  = scale*sqrt( (absr/scale)**2 + (s/scale)**2 )
                  end if
               else if (absr .ge. s) then
                  q     = sqrt(absr + s)*sqrt(absr - s)
               else
                  fitok = .false.
               end if

               if (fitok) then

*                 Compute a minimizer of the fitted cubic.

                  if (xw .lt. zero) q = - q
                  s  = gbest -  r - q
                  q  = gbest - gw - q - q
               end if
            end if
*           ------------------------------------------------------------
*           Construct an artificial interval  (artifa, artifb)  in which
*           the new estimate of a minimizer must lie.  Set a default
*           value of xtry that will be used if the polynomial fit fails.
*           ------------------------------------------------------------
            artifa = a
            artifb = b
            if (.not. braktd) then

*              A minimizer has not been bracketed.  Set an artificial
*              upper bound by expanding the interval  xw  by a suitable
*              factor.

               xtry   = - factor*xw
               artifb =   xtry
               if (alfbst + xtry .lt. alfmax) factor = five*factor

            else if (extrap) then

*              The points are configured for an extrapolation.
*              Set a default value of  xtry  in the interval  (a, b)
*              that will be used if the polynomial fit is rejected.  In
*              the following,  dtry  and  daux  denote the lengths of
*              the intervals  (a, b)  and  (0, xw)  (or  (xw, 0),  if
*              appropriate).  The value of  xtry is the point at which
*              the exponents of  dtry  and  daux  are approximately
*              bisected.

               daux = abs( xw )
               dtry = b - a
               if (daux .ge. dtry) then
                  xtry = five*dtry*(point1 + dtry/daux)/eleven
               else
                  xtry = half * sqrt( daux ) * sqrt( dtry )
               end if
               if (xw .gt. zero)   xtry = - xtry
               if (debug) write(nout, 2400) xtry, daux, dtry

*              Reset the artificial bounds.  If the point computed by
*              extrapolation is rejected,  xtry will remain at the
*              relevant artificial bound.

               if (xtry .le. zero) artifa = xtry
               if (xtry .gt. zero) artifb = xtry
            else

*              The points are configured for an interpolation.  The
*              default value xtry bisects the interval of uncertainty.
*              the artificial interval is just (a, b).

               xtry   = xmidpt
               if (debug) write(nout, 2300) xtry
               if (nsamea .ge. 3  .or.  nsameb .ge. 3) then

*                 If the interpolation appears to be overestimating the
*                 distance to a minimizer,  damp the interpolation.

                  factor = factor / five
                  s      = factor * s
               else
                  factor = one
               end if
            end if
*           ------------------------------------------------------------
*           The polynomial fits give  (s/q)*xw  as the new step.
*           Reject this step if it lies outside  (artifa, artifb).
*           ------------------------------------------------------------
            if (q .ne. zero) then
               if (q .lt. zero) s = - s
               if (q .lt. zero) q = - q
               if (s*xw .ge. q*artifa  .and.  s*xw .le. q*artifb) then

*                 Accept the polynomial fit.

                  if (abs( s*xw ) .ge. q*tol) then
                     xtry = (s/q)*xw
                  else
                     xtry = zero
                  end if
                  if (debug) write(nout, 2500) xtry
               end if
            end if
         end if
      end if

*     ==================================================================

      if (.not. done) then
         alfa  = alfbst + xtry
         if (braktd  .or.  alfa .lt. alfmax - tolmax) then

*           The function must not be evaluated too close to a or b.
*           (It has already been evaluated at both those points.)

            if (xtry .le. a + tol  .or.  xtry .ge. b - tol) then
               if (half*(a + b) .le. zero) then
                  xtry = - tol
               else
                  xtry =   tol
               end if
               alfa = alfbst + xtry
            end if
         else

*           The step is close to, or larger than alfmax, replace it by
*           alfmax to force evaluation of  f  at the boundary.

            braktd = .true.
            xtry   = alfmax - alfbst
            alfa   = alfmax
         end if
      end if

*     ------------------------------------------------------------------
*     Exit.
*     ------------------------------------------------------------------
      if (done) then
         if      (badfun) then
            iExit = 8           ! bad arguments
         else if (found) then
            if (alfbst .lt. alfmax) then
               iExit = 1        ! Sufficient decrease
            else
               iExit = 2        ! Suff. Decrease on the boundary
            end if
         else if (moved ) then
            iExit = 3           ! Decr at boundary or max funs
         else if (quitF) then
            iExit = 7           ! No new point after max funs
         else if (crampd) then
            iExit = 4           ! alfmax too mall
         else
            iExit = 6           ! [a,b] too small
         end if
      end if

      if (debug) write(nout, 3000)
      return

 1000 format(/'     g0  tolabs  alfmax        ', 1p, 2e22.14,   e16.8
     &       /' targtg  tolrel   epsaf        ', 1p, 2e22.14,   e16.8
     &       /' crampd                        ',  l3)
 1100 format(/' alfa    ftry    gtry          ', 1p, 2e22.14,   e16.8)
 1200 format(/' a       b       b - a   tol   ', 1p, 2e22.14,  2e16.8
     &       /' nsamea  nsameb  numf          ', 3i3
     &       /' braktd  extrap  closef  imprvd', 4l3
     &       /' found   quitI                 ', 2l3
     &       /' alfbst  fbest   gbest         ', 1p, 3e22.14
     &       /' alfaw   fw      gw            ', 1p, 3e22.14)
 2100 format( ' Cubic.   ')
 2200 format( ' Parabola.')
 2300 format( ' Bisection.              xmidpt', 1p,  e22.14)
 2400 format( ' Geo. bisection. xtry,daux,dtry', 1p, 3e22.14)
 2500 format( ' Polynomial fit accepted.  xtry', 1p,  e22.14)
 3000 format( ' ----------------------------------------------------'/)

      end ! subroutine lsrchc

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsrchq
     &   ( iExit , first , debug , done  , imprvd,
     &     maxf  , numf  , nout  ,
     &     alfmax, alfsml, epsaf ,
     &     g0    , targtg, ftry  ,
     &     tolabs, tolrel, toltny,
     &     alfa  , alfbst, fbest ,
     &     braktd, crampd, extrap, moved , vset  , wset,
     &     nsamea, nsameb,
     &     a     , b     , fa    , factor,
     &     xtry  , xw    , fw    , xv    , fv    , tolmax )

      implicit
     &     none
      logical
     &     first, debug, done, imprvd, braktd, crampd, extrap,
     &     moved, vset , wset
      integer
     &     iExit, maxf, numf, nout, nsamea, nsameb
      double precision
     &     alfmax, alfsml, epsaf, g0, targtg, ftry, tolabs, tolrel,
     &     toltny, alfa, alfbst, fbest, a, b, factor, xtry,
     &     xw, fw, gw, tolmax

*     ==================================================================
*     lsrchq  finds a sequence of improving estimates of a minimizer of
*     the univariate function f(alpha) in the interval (0,alfmax].
*     f(alpha) is a smooth function such that  f(0) = 0  and  f'(0) < 0.
*     lsrchq  requires  f(alpha) (but not f'(alpha)) to be evaluated
*     in the interval.  New estimates of a minimizer are computed using
*     safeguarded quadratic interpolation.
*
*     Reverse communication is used to allow the calling program to
*     evaluate f.  Some of the parameters must be set or tested by the
*     calling program.  The remainder would ordinarily be local
*     variables.
*
*     Input parameters (relevant to the calling program)
*     --------------------------------------------------
*
*     first         must be true on the first entry.  It is subsequently
*                   altered by lsrchq.
*
*     debug         specifies whether detailed output is wanted.
*
*     maxf          is an upper limit on the number of times lsrchq is
*                   to be entered consecutively with done = false
*                   (following an initial entry with first = true).
*
*     alfa          is the first estimate of a minimizer.  alfa is
*                   subsequently altered by lsrchq (see below).
*
*     alfmax        is the upper limit of the interval to be searched.
*
*     alfsml        is intended to prevent inefficiency when a minimizer
*                   is very small, for cases where the calling program
*                   would prefer to redefine f'(alfa).  alfsml is
*                   allowed to be zero.  Early termination will occur if
*                   lsrchq determines that a minimizer lies somewhere in
*                   the interval [0, alfsml) (but not if alfmax is
*                   smaller that alfsml).
*
*     epsaf         is an estimate of the absolute precision in the
*                   computed value of f(0).
*
*     ftry          the value of f at the new point
*                   alfa = alfbst + xtry.
*
*     g0            is the value of f'(0).  g0 must be negative.
*
*     tolabs,tolrel define a function tol(alfa) = tolrel*alfa + tolabs
*                   such that if f has already been evaluated at alfa,
*                   it will not be evaluated closer than tol(alfa).
*                   These values may be reduced by lsrchq.
*
*     targtg        is the target value of abs(f'(alfa)). The search
*                   is terminated when
*                    abs(f'(alfa)) le targtg and f(alfa) lt 0.
*
*     toltny        is the smallest value that tolabs is allowed to be
*                   reduced to.
*
*     Output parameters (relevant to the calling program)
*     ---------------------------------------------------
*
*     imprvd        is true if the previous alfa was the best point so
*                   far.  Any related quantities should be saved by the
*                   calling program (e.g., arrays) before paying
*                   attention to the variable done.
*
*     done = false  means the calling program should evaluate ftry
*                   for the new trial step alfa, and reenter lsrchq.
*
*     done = true   means that no new alfa was calculated.  The value
*                   of iExit gives the result of the search as follows
*
*                   iExit = 1 means the search has terminated
*                             successfully with alfbst < alfmax.
*
*                   iExit = 2 means the search has terminated
*                             successfully with alfbst = alfmax.
*
*                   iExit = 3 means that the search failed to find a
*                             point of sufficient decrease in maxf
*                             functions, but a lower point was found.
*
*                   iExit = 4 means alfmax is so small that a search
*                             should not have been attempted.
*
*                   iExit = 5 means that the search was terminated
*                             because of alfsml (see above).
*
*                   iExit = 6 means the search has failed to find a
*                             useful step.  The interval of uncertainty
*                             is [0,b] with b < 2*tolabs. A minimizer
*                             lies very close to alfa = 0, or f'(0) is
*                             not sufficiently accurate.
*
*                   iExit = 7 if no better point could be found after
*                             maxf  function calls.
*
*                   iExit = 8 means the input parameters were bad.
*                             alfmax le toltny  or  g0 ge zero.
*                             No function evaluations were made.
*
*     numf          counts the number of times lsrchq has been entered
*                   consecutively with done = false (i.e., with a new
*                   function value ftry).
*
*     alfa          is the point at which the next function ftry must
*                   be computed.
*
*     alfbst        should be accepted by the calling program as the
*                   approximate minimizer, whenever lsrchq returns
*                   iExit = 1, 2 or 3.
*
*     fbest         will be the corresponding value of f.
*
*     The following parameters retain information between entries
*     -----------------------------------------------------------
*
*     braktd        is false if f has not been evaluated at the far end
*                   of the interval of uncertainty.  In this case, the
*                   point b will be at alfmax + tol(alfmax).
*
*     crampd        is true if alfmax is very small (le tolabs).  If the
*                   search fails, this indicates that a zero step should
*                   be taken.
*
*     extrap        is true if alfbst has moved at least once and xv
*                   lies outside the interval of uncertainty.  In this
*                   case, extra safeguards are applied to allow for
*                   instability in the polynomial fit.
*
*     moved         is true if a better point has been found, i.e.,
*                   alfbst gt 0.
*
*     vset          records whether a third-best point has been defined.
*
*     wset          records whether a second-best point has been
*                   defined.  It will always be true by the time the
*                   convergence test is applied.
*
*     nsamea        is the number of consecutive times that the
*                   left-hand end point of the interval of uncertainty
*                   has remained the same.
*
*     nsameb        similarly for the right-hand end.
*
*     a, b, alfbst  define the current interval of uncertainty.
*                   A minimizer lies somewhere in the  interval
*                   [alfbst + a, alfbst + b].
*
*     alfbst        is the best point so far.  It lies strictly within
*                   [atrue,btrue]  (except when alfbst has not been
*                   moved, in which case it lies at the left-hand end
*                   point).  Hence we have a .le. 0 and b .gt. 0.
*
*     fbest         is the value of f at the point alfbst.
*
*     fa            is the value of f at the point alfbst + a.
*
*     factor        controls the rate at which extrapolated estimates of
*                   alfa  may expand into the interval of uncertainty.
*                   Factor is not used if a minimizer has been bracketed
*                   (i.e., when the variable braktd is true).
*
*     fv, fw        are the values of f at the points alfbst + xv  and
*                   alfbst + xw.  They are not defined until  vset  or
*                   wset  are true.
*
*     xtry          is the trial point within the shifted interval
*                   (a, b).  The new trial function value must be
*                   computed at the point alfa = alfbst + xtry.
*
*     xv            is such that alfbst + xv is the third-best point.
*                   It is not defined until vset is true.
*
*     xw            is such that alfbst + xw is the second-best point.
*                   It is not defined until wset is true.  In some
*                   cases,  xw will replace a previous xw that has a
*                   lower function but has just been excluded from
*                   (a,b).
*
*     Systems Optimization Laboratory, Stanford University, California.
*     Original version February 1982.  Rev. May 1983.
*     Original F77 version 22-August-1985.
*     17 Jul 1997: Removed saved variables for thread-safe version.
*     31 Jul 1999: Current version.
*     ==================================================================
      logical
     &     badfun, closef, found, quitF, quitFZ, quitI , quitS,
     &     setxv , xinxw
      double precision
     &     artifa, artifb, daux, dtry, endpnt, fa, fv, gv,
     &     q, s, tol, truea, trueb, xmidpt, xv
*     ------------------------------------------------------------------
      double precision   zero,          point1,          half
      parameter         (zero = 0.0d+0, point1 = 0.1d+0, half = 0.5d+0)
      double precision   one,           two,             five
      parameter         (one  = 1.0d+0, two    = 2.0d+0, five = 5.0d+0)
      double precision   ten,           eleven
      parameter         (ten  = 1.0d+1, eleven = 1.1d+1               )
*     ------------------------------------------------------------------
*     Local variables
*     ===============
*
*     closef     is true if the worst function fv is within epsaf of
*                fbest (up or down).
*
*     found      is true if the sufficient decrease conditions holds at
*                alfbst.
*
*     quitF      is true when  maxf  function calls have been made.
*
*     quitFZ     is true when the three best function values are within
*                epsaf of each other, and the new point satisfies
*                fbest le ftry le fbest+epsaf.
*
*     quitI      is true when the interval of uncertainty is less than
*                2*tol.
*
*     quitS      is true as soon as alfa is too small to be useful;
*                i.e., btrue le alfsml.
*
*     xinxw      is true if xtry is in (xw,0) or (0,xw).
*     ------------------------------------------------------------------

      imprvd = .false.
      badfun = .false.
      quitF  = .false.
      quitFZ = .false.
      quitS  = .false.
      quitI  = .false.

      if (first) then
*        ---------------------------------------------------------------
*        First entry.  Initialize various quantities, check input data
*        and prepare to evaluate the function at the initial step alfa.
*        ---------------------------------------------------------------
         first  = .false.
         numf   = 0
         alfbst = zero
         badfun = alfmax .le. toltny  .or.  g0 .ge. zero
         done   = badfun
         moved  = .false.

         if (.not. done) then
            braktd = .false.
            crampd = alfmax .le. tolabs
            extrap = .false.
            vset   = .false.
            wset   = .false.
            nsamea = 0
            nsameb = 0

            tolmax = tolrel*alfmax + tolabs
            a      = zero
            b      = alfmax + tolmax
            fa     = zero
            factor = five
            tol    = tolabs
            xtry   = alfa
            if (debug) then
               write(nout, 1000)
     &              g0    , tolabs, alfmax,
     &              targtg, tolrel, epsaf , crampd
            end if
         end if
      else
*        ---------------------------------------------------------------
*        Subsequent entries.  The function has just been evaluated at
*        alfa = alfbst + xtry,  giving ftry.
*        ---------------------------------------------------------------
         if (debug)
     &   write(nout, 1100)
     &        alfa, ftry

         numf   = numf   + 1
         nsamea = nsamea + 1
         nsameb = nsameb + 1

         if (.not. braktd) then
            tolmax = tolabs + tolrel*alfmax
            b      = alfmax - alfbst + tolmax
         end if

*        Check if xtry is in the interval (xw,0) or (0,xw).

         if (wset) then
            xinxw =        zero .lt. xtry  .and.  xtry .le. xw
     &               .or.    xw .le. xtry  .and.  xtry .lt. zero
         else
            xinxw = .false.
         end if

         imprvd = ftry .lt. fbest
         if (vset) then
            closef = abs( fbest - fv ) .le. epsaf
         else
            closef = .false.
         end if

         if (imprvd) then

*           We seem to have an improvement.  The new point becomes the
*           origin and other points are shifted accordingly.

            if (wset) then
               xv     = xw - xtry
               fv     = fw
               vset   = .true.
            end if

            xw     = zero - xtry
            fw     = fbest
            wset   = .true.
            fbest  = ftry
            alfbst = alfa
            moved  = .true.

            a      = a    - xtry
            b      = b    - xtry
            extrap = .not. xinxw

*           Decrease the length of (a,b).

            if (xtry .ge. zero) then
               a      = xw
               fa     = fw
               nsamea = 0
            else
               b      = xw
               nsameb = 0
               braktd = .true.
            end if
         else if (closef  .and.  ftry - fbest .lt. epsaf) then

*           Quit if there has been no progress and ftry, fbest, fw
*           and fv are all within epsaf of each other.

            quitFZ = .true.
         else

*           The new function value is no better than the current best
*           point.  xtry must an end point of the new (a,b).

            if (xtry .lt. zero) then
               a      = xtry
               fa     = ftry
               nsamea = 0
            else
               b      = xtry
               nsameb = 0
               braktd = .true.
            end if

*           The origin remains unchanged but xtry may qualify as xw.

            if (wset) then
               if (ftry .lt. fw) then
                  xv     = xw
                  fv     = fw
                  vset   = .true.

                  xw     = xtry
                  fw     = ftry
                  if (moved) extrap = xinxw
               else if (moved) then
                  if (vset) then
                     setxv = ftry .lt. fv  .or.  .not. extrap
                  else
                     setxv = .true.
                  end if

                  if (setxv) then
                     if (vset  .and.  xinxw) then
                        xw = xv
                        fw = fv
                     end if
                     xv   = xtry
                     fv   = ftry
                     vset = .true.
                  end if
               else
                  xw  = xtry
                  fw  = ftry
               end if
            else
               xw     = xtry
               fw     = ftry
               wset   = .true.
            end if
         end if

*        ---------------------------------------------------------------
*        Check the termination criteria.
*        ---------------------------------------------------------------
         tol    = tolabs + tolrel*alfbst
         truea  = alfbst + a
         trueb  = alfbst + b

         found  = moved  .and.  abs(fa - fbest) .le. -a*targtg
         quitF  = numf  .ge. maxf
         quitI  = b - a .le. tol + tol
         quitS  = trueb .le. alfsml

         if (quitI  .and.  .not. moved) then

*           The interval of uncertainty appears to be small enough,
*           but no better point has been found.  Check that changing
*           alfa by b-a changes f by less than epsaf.

            tol    = tol/ten
            tolabs = tol
            quitI  = abs(fw) .le. epsaf  .or.  tol .le. toltny
         end if

         done  = quitF  .or.  quitFZ  .or.  quitS  .or.  quitI
     &                  .or.  found

         if (debug) then
            write(nout, 1200)
     &           truea    , trueb , b-a   , tol   ,
     &           nsamea   , nsameb, numf  ,
     &           braktd   , extrap, closef, imprvd,
     &           found    , quitI , quitFZ, quitS ,
     &           alfbst   , fbest ,
     &           alfbst+xw, fw
            if (vset) then
               write(nout, 1300) alfbst + xv, fv
            end if
         end if

*        ---------------------------------------------------------------
*        Proceed with the computation of an estimate of a minimizer.
*        The choices are...
*        1. Parabolic fit using function values only.
*        2. Damped parabolic fit if the regular fit appears to be
*           consistently overestimating the distance to a minimizer.
*        3. Bisection, geometric bisection, or a step of tol if the
*           parabolic fit is unsatisfactory.
*        ---------------------------------------------------------------
         if (.not. done) then
            xmidpt = half*(a + b)
            s      = zero
            q      = zero

*           ============================================================
*           Fit a parabola.
*           ============================================================
*           See if there are two or three points for the parabolic fit.

            gw = (fw - fbest)/xw
            if (vset  .and.  moved) then

*              Three points available.  Use fbest, fw and fv.

               gv = (fv - fbest)/xv
               s  = gv - (xv/xw)*gw
               q  = two*(gv - gw)
               if (debug) write(nout, 2200)
            else

*              Only two points available.  Use fbest, fw and g0.

               if (moved) then
                  s  = g0 - two*gw
               else
                  s  = g0
               end if
               q = two*(g0 - gw)
               if (debug) write(nout, 2100)
            end if

*           ------------------------------------------------------------
*           Construct an artificial interval (artifa, artifb) in which
*           the new estimate of the steplength must lie.  Set a default
*           value of  xtry  that will be used if the polynomial fit is
*           rejected. In the following, the interval (a,b) is considered
*           the sum of two intervals of lengths  dtry  and  daux, with
*           common end point the best point (zero).  dtry is the length
*           of the interval into which the default xtry will be placed
*           and endpnt denotes its non-zero end point.  The magnitude of
*           xtry is computed so that the exponents of dtry and daux are
*           approximately bisected.
*           ------------------------------------------------------------
            artifa = a
            artifb = b
            if (.not. braktd) then

*              A minimizer has not yet been bracketed.
*              Set an artificial upper bound by expanding the interval
*              xw  by a suitable factor.

               xtry   = - factor*xw
               artifb =   xtry
               if (alfbst + xtry .lt. alfmax) factor = five*factor
            else if (vset .and. moved) then

*              Three points exist in the interval of uncertainty.
*              Check if the points are configured for an extrapolation
*              or an interpolation.

               if (extrap) then

*                 The points are configured for an extrapolation.

                  if (xw .lt. zero) endpnt = b
                  if (xw .gt. zero) endpnt = a
               else

*                 If the interpolation appears to be overestimating the
*                 distance to a minimizer,  damp the interpolation step.

                  if (nsamea .ge. 3  .or.   nsameb .ge. 3) then
                     factor = factor / five
                     s      = factor * s
                  else
                     factor = one
                  end if

*                 The points are configured for an interpolation.  The
*                 artificial interval will be just (a,b).  Set endpnt so
*                 that xtry lies in the larger of the intervals (a,b)
*                 and  (0,b).

                  if (xmidpt .gt. zero) then
                     endpnt = b
                  else
                     endpnt = a
                  end if

*                 If a bound has remained the same for three iterations,
*                 set endpnt so that  xtry  is likely to replace the
*                 offending bound.

                  if (nsamea .ge. 3) endpnt = a
                  if (nsameb .ge. 3) endpnt = b
               end if

*              Compute the default value of  xtry.

               dtry = abs( endpnt )
               daux = b - a - dtry
               if (daux .ge. dtry) then
                  xtry = five*dtry*(point1 + dtry/daux)/eleven
               else
                  xtry = half*sqrt( daux )*sqrt( dtry )
               end if
               if (endpnt .lt. zero) xtry = - xtry
               if (debug) write(nout, 2500) xtry, daux, dtry

*              If the points are configured for an extrapolation set the
*              artificial bounds so that the artificial interval lies
*              within (a,b).  If the polynomial fit is rejected,  xtry
*              will remain at the relevant artificial bound.

               if (extrap) then
                  if (xtry .le. zero) then
                     artifa = xtry
                  else
                     artifb = xtry
                  end if
               end if
            else

*              The gradient at the origin is being used for the
*              polynomial fit.  Set the default xtry to one tenth xw.

               if (extrap) then
                  xtry = - xw
               else
                  xtry   = xw/ten
               end if
               if (debug) write(nout, 2400) xtry
            end if

*           ------------------------------------------------------------
*           The polynomial fits give (s/q)*xw as the new step.  Reject
*           this step if it lies outside (artifa, artifb).
*           ------------------------------------------------------------
            if (q .ne. zero) then
               if (q .lt. zero) s = - s
               if (q .lt. zero) q = - q
               if (s*xw .ge. q*artifa   .and.   s*xw .le. q*artifb) then

*                 Accept the polynomial fit.

                  if (abs( s*xw ) .ge. q*tol) then
                     xtry = (s/q)*xw
                  else
                     xtry = zero
                  end if
                  if (debug) write(nout, 2600) xtry
               end if
            end if
         end if
      end if
*     ==================================================================

      if (.not. done) then
         alfa  = alfbst + xtry
         if (braktd  .or.  alfa .lt. alfmax - tolmax) then

*           The function must not be evaluated too close to a or b.
*           (It has already been evaluated at both those points.)

            xmidpt = half*(a + b)
            if (xtry .le. a + tol  .or.  xtry .ge. b - tol) then
               if (xmidpt .le. zero) then
                  xtry = - tol
               else
                  xtry =   tol
               end if
            end if

            if (abs( xtry ) .lt. tol) then
               if (xmidpt .le. zero) then
                  xtry = - tol
               else
                  xtry =   tol
               end if
            end if
            alfa  = alfbst + xtry
         else

*           The step is close to or larger than alfmax, replace it by
*           alfmax to force evaluation of the function at the boundary.

            braktd = .true.
            xtry   = alfmax - alfbst
            alfa   = alfmax
         end if
      end if
*     ------------------------------------------------------------------
*     Exit.
*     ------------------------------------------------------------------
      if (done) then
         if      (badfun) then
            iExit = 8
         else if (quitS ) then
            iExit = 5
         else if (found) then
            if (alfbst .lt. alfmax) then
               iExit = 1
            else
               iExit = 2
            end if
         else if (moved ) then
            iExit = 3
         else if (quitF) then
            iExit = 7
         else if (crampd) then
            iExit = 4
         else
            iExit = 6
         end if
      end if

      if (debug) write(nout, 3000)
      return

 1000 format(/'     g0  tolabs  alfmax        ', 1p, 2e22.14,   e16.8
     &       /' targtg  tolrel   epsaf        ', 1p, 2e22.14,   e16.8
     &       /' crampd                        ',  l3)
 1100 format(/' alfa    ftry                  ', 1p,2e22.14          )
 1200 format(/' a       b       b - a   tol   ', 1p,2e22.14,   2e16.8
     &       /' nsamea  nsameb  numf          ', 3i3
     &       /' braktd  extrap  closef  imprvd', 4l3
     &       /' found   quitI   quitFZ  quitS ', 4l3
     &       /' alfbst  fbest                 ', 1p,2e22.14
     &       /' alfaw   fw                    ', 1p,2e22.14)
 1300 format( ' alfav   fv                    ', 1p,2e22.14 /)
 2100 format( ' Parabolic fit,    two points. ')
 2200 format( ' Parabolic fit,  three points. ')
 2400 format( ' Exponent reduced.  Trial point', 1p,  e22.14)
 2500 format( ' Geo. bisection. xtry,daux,dtry', 1p, 3e22.14)
 2600 format( ' Polynomial fit accepted.  xtry', 1p,  e22.14)
 3000 format( ' ----------------------------------------------------'/)

      end ! subroutine lsrchq
