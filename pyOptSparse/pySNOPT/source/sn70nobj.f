!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     File  sn70nobj.f
!
!     s7chkA   s7chkG   s7chkp   s7fixX   s7Jac    s7pert   s7step
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s7chkA
     &   ( iExit, Status, lvlChk, userfg,
     &     nF, n, iGmax, eGmax, Elem, Set,
     &     iGfun, jGvar, lenG, neG,
     &     x , F , G ,
     &     x1, F1, G1, w, w1, y,
     &     iw, leniw, rw, lenrw, cu, lencu,
     &     iu, leniu, ru, lenru )

      implicit
     &     none
      external
     &     userfg
      integer
     &     iExit, iGmax, lenG, leniw, lencu, leniu, lenru,
     &     lenrw, lvlChk, nF, n, neG, Status, Elem(nF), Set(n),
     &     iGfun(lenG), jGvar(lenG), iw(leniw), iu(leniu)
      double precision
     &     eGmax, F(nF), F1(nF),
     &     G(lenG), G1(lenG), ru(lenru), rw(lenrw), w(n), w1(n),
     &     x(n), x1(n), y(n)
      character
     &     cu(lencu)*8

!     ==================================================================
!     s7chkA   does the work for snchkA, the stand-alone derivative
!     checker for snOptA.
!
!     01 Jun 2006: First version of s7chkA.
!     22 Apr 2007: Some features added and exit messages revised.
!     ==================================================================
      external
     &     idamax, dnormi
      character
     &     key*4, str*120
      logical
     &     first, badGijs
      integer
     &     i, idamax, j, jj, jGmax, k, kGmax, needF, needG, nG, nSet,
     &     rightG, wrongG
      double precision
     &     damper, dnormi, dx, eps0, eps5, eps, err, fdint1,
     &     Gdummy, Gfd, Gmax, Gij, xj, yj, s1eps
!     ------------------------------------------------------------------
      double precision   zero,           one
      parameter         (zero  = 0.0d+0, one    = 1.0d+0)
      logical            yes,            no
      parameter         (yes   = .true., no     = .false.)
      double precision   ok
      parameter         (ok    = 0.1d+0)
      character          lWrong*4       , lRight*4
      data               lWrong /'bad?'/, lRight /'ok  '/
!     ------------------------------------------------------------------
      iExit   = 0
      nSet    = 0               ! number of derivatives to be checked
      badGijs = .false.         ! so far so good ...

      if (lvlChk .lt. 0) go to 900

      Status = 1                ! Assume this is the first call
      needF  = 1                ! Get both F and G
      needG  = 1
      call userfg
     &   ( Status, n, x,
     &     needF, nF, F,
     &     needG, lenG, G,
     &     cu, lencu, iu, leniu, ru, lenru )
      if (Status .lt. 0) go to 900

!     ------------------------------------------------------------------
!     Assign some constants.
!     ------------------------------------------------------------------
      Gdummy = rw( 69)          ! definition of an 'unset' value

      eps    = s1eps ( )        ! machine precision
      eps0   = eps**(0.80d+0)
      eps5   = eps**(0.2d+0)    ! eps**(1/5)
      fdint1 = sqrt(eps)

      iGmax  = 0
      eGmax  = zero

!     ------------------------------------------------------------------
!     Cheap test.
!     ------------------------------------------------------------------
      yj    = one/n
      do  j =  1, n
         y(j) =   yj
         yj   = - yj * 0.99999d+0
      end do

!     Thumb through G and count the assigned elements in each column.
!     y(j) is set to zero if the jth column has some unknown elements.

      nSet   = 0
      call iload ( n, (0), Set, 1 )

      do k = 1, neG
         j = jGvar(k)
         if (G(k) .eq. Gdummy) then ! G(k) was not set
            y(j)   = zero           ! x(j) is not perturbed
         else
            Set(j) = Set(j) + 1
            nSet   = nSet   + 1
         end if
      end do

      if (nSet .eq. 0) then
         call snPRNT( 11,
     &      ' No derivatives to be checked',
     &      iw, leniw )

      else
         if (lvlChk .eq. 0) then
            call snPRNT( 11,
     &           ' Performing a cheap test of the derivatives...',
     &           iw, leniw )
         else if (lvlChk .gt. 0) then
            call snPRNT( 11,
     &           ' Performing a full test of the derivatives...',
     &           iw, leniw )
         end if

!        ---------------------------------------------------------------
!        Compute F at  x1 = x + dx*y.
!        s6usrf reduces the step until F(x1) is well-defined.
!        ---------------------------------------------------------------
         dx   = fdint1 * (one + dnormi( n, x, 1 ))

         do j = 1, n
            x1(j) = x(j) + dx*y(j)
         end do

         Status = 0
         needF  = 1
         needG  = 0             ! G isn't needed

         call s6usrf
     &      ( Status, n, x, x1, damper, userfg,
     &        needF, nF  , F1,
     &        needG, lenG, G1,
     &        cu, lencu, iu, leniu, ru, lenru )
         if (Status .lt. 0) go to 900

!        Set   w1 = (F1 - F)/dx - G*y.  This should be small.

         dx = damper*dx

         do i  = 1, nF
            w1(i)  = (F1(i) - F(i))/dx
         end do

         do k = 1, neG
            i = iGfun(k)
            j = jGvar(k)
            w1(i) = w1(i) - G(k)*y(j)
         end do

         iGmax  = idamax( nF, w1, 1 )
          Gmax  = (F1(iGmax) - F(iGmax))/dx
         eGmax  = abs(w1(iGmax) )/(one + abs( Gmax ))

         if (eGmax .le. ok) then
            call snPRNT( 11,
     &           ' All the assigned derivatives seem to be OK.',
     &           iw, leniw )
         else
            call snPRNT( 11,
     &           ' XXX  Some derivatives seem to be incorrect.',
     &           iw, leniw )
         end if

         write(str, 1601) eGmax, iGmax
         call snPRNT( 11, str, iw, leniw )
         call snPRNT(  1, ' ', iw, leniw )

         badGijs = eGmax .gt. ok

      end if

      if (lvlChk .eq. 0) go to 900

!     ------------------------------------------------------------------
!     Proceed with column-wise verification.
!     ------------------------------------------------------------------
      write(str, 2001)
      call snPRNT( 1, str, iw, leniw )

      nG     =  0
      rightG =  0
      wrongG =  0

      eGmax  = -one
      iGmax  =  0
      jGmax  =  0

      do j  = 1, n
!        ------------------------------------------------------------
!        Estimate the jth column of G.
!        Check the estimate against the user-supplied values.
!        Don't bother printing a line for an exact zero.
!        Look for nonzeros that are not in the sparse data structure.
!        ------------------------------------------------------------
         if (Set(j) .gt. 0) then
            xj   = x(j)
            dx   = fdint1 * (one + abs( xj ))
            x(j) = xj + dx

            Status = 0
            call userfg
     &         ( Status, n, x,
     &           needF, nF, F1,
     &           needG, lenG, G1,
     &           cu, lencu, iu, leniu, ru, lenru )
            if (Status .lt. 0) go to 900

            ! Find the elements in this column
            call dload ( nF, zero, w, 1 )
            do  k = 1, neG
               jj = jGvar(k)
               if (jj .eq. j) then
                  i       = iGfun(k)
                  w(i)    = G(k)
                  Elem(i) = k
               end if
            end do

            first = yes

            do  i = 1, nF
               Gij   = w(i)
               Gfd   = (F1(i) - F(i))/dx
               w1(i) = Gfd
               k     = Elem(i)

               if (Gij .ne. Gdummy) then
                  nG  = nG + 1
                  err = max( abs(Gfd - Gij)/(one + abs(Gij)),
     &                       abs(Gfd - Gij)/(one + abs(Gfd)) )

                  if (eGmax .lt. err) then
                     eGmax  = err
                     iGmax  = i
                     jGmax  = j
                     kGmax  = k
                  end if

                  if (err .le.  eps5) then
                     key    = lRight
                     rightG = rightG + 1
                  else
                     key    = lWrong
                     wrongG = wrongG + 1
                  end if

                  if (abs(Gij) + err .gt. eps0) then
                     if (first) then
                        write(str, 2101) j, xj, dx, k, i, Gij, Gfd, key
                        call snPRNT( 11, str, iw, leniw )
                        first  = no
                     else
                        write(str, 2201)            k, i, Gij, Gfd, key
                        call snPRNT(  1, str, iw, leniw )
                     end if
                  end if
               end if

               w1(i) = Gdummy    ! Done with this row
            end do

!           Check that elements remaining in w are zero.

            do i = 1, nF
               if (w1(i) .ne. Gdummy) then ! w1(i) not yet processed
                  err = abs(w1(i))

                  if (err .gt. eps0) then
                     wrongG = wrongG + 1

                     if (first) then
                        write(str, 2301) j, i, w1(i), lWrong
                        call snPRNT( 1, str, iw, leniw )
                        first = no
                     else
                        write(str, 2302)    i, w1(i), lWrong
                        call snPRNT( 1, str, iw, leniw )
                     end if

                     if (eGmax .lt. err) then
                        eGmax = err
                        iGmax = i
                        jGmax = j
                     end if
                  end if
               end if
            end do
            x(j)  = xj
         end if ! set
      end do

!     ------------------------------------------------------------------
!     Final tally of the good, the bad and the ugly.
!     ------------------------------------------------------------------
      if (wrongG .eq. 0) then
         write(str, 2501) rightG
      else
         write(str, 2601) wrongG
      end if

      call snPRNT( 11, str, iw, leniw )
      write(str, 2701) eGmax, iGmax, jGmax, kGmax
      call snPRNT( 11, str, iw, leniw )
!      call snPRNT(  1, ' ', iw, leniw )

      badGijs = badGijs  .or.  wrongG .gt. 0

!     ------------------------------------------------------------------
!     Set iExit for exceptions.
!     ------------------------------------------------------------------
  900 if (lvlChk .lt. 0  .or.  nSet .eq. 0) then
         iExit = 106            ! no derivatives were checked
      else if (badGijs) then
         iExit = 55             ! Some bad derivatives
      else if (Status .lt. 0) then
         if (Status .eq. -1) then
            iExit = 62          ! undefined function at the first point
         else
            iExit = 71          ! terminated during function evaluation
         end if
      end if

      return

 1601 format(' -->  The largest discrepancy was', 1p, e12.2,
     &       '  in row', i6)
 2001 format(' Column       x(j)        dx(j)', 4x,
     &       ' Element', 7x, 'Row', 7x,
     &       ' Derivative    Difference approxn')
 2101 format(i7, 1p, e16.8, e10.2, 2i10,               2e18.8, 2x, a4)
 2201 format(      33x, 2i10, 1pe18.8, e18.8, 2x, a4)
 2301 format(i7, 2x, 'Nonzero not in sparse structure ', '??', 4x, i6,
     &          18x, 1p, e18.8, 2x, a4 )
 2302 format(    9x, 'Nonzero not in sparse structure ', '??', 4x, i6,
     &          18x, 1p, e18.8, 2x, a4 )
 2501 format(' All',  i7, ' assigned derivatives seem to be OK.')
 2601 format(' XXX  There seem to be', i6, ' incorrect derivatives.' )
 2701 format(' -->  The largest relative error was', 1p, e12.2,
     &       '   in row', i6, ',  column', i6, ', element', i6 )

      end ! subroutine s7chkA

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s7chkG
     &   ( iExit, n, nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj,
     &     fgwrap, fgcon, fgobj,
     &     x, x1, bl, bu, fObj, gObj,
     &     ne, nlocJ, locJ, indJ, negCon, nlocG, locG,
     &     fCon, gCon, gObj2, fCon2, gCon2, y, y1,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     fgwrap, fgcon, fgobj
      integer
     &     iExit, lencu, leniu, lenru, lencw, leniw, lenrw, n, ne,
     &     negCon, nnL, nnCon0, nnCon, nnJac, nnObj0, nnObj, nlocJ,
     &     nlocG, locG(nlocG), indJ(ne), locJ(nlocJ), iu(leniu),
     &     iw(leniw)
      double precision
     &     bl(nnL), bu(nnL), fObj,
     &     gObj (nnObj0), fCon (nnCon0), gCon (negCon),
     &     gObj2(nnObj0), fCon2(nnCon0), gCon2(negCon),
     &     x(n), x1(n), y(nnL), y1(nnCon0), ru(lenru),
     &     rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     s7chkG  verifies the objective and constraint gradients using
*     finite differences.
*
*     First, a cheap heuristic test is performed, as in
*     subroutine chkgrd by the following authors:
*     Philip E. Gill, Walter Murray, Susan M. Picken and Hazel M. Barber
*     DNAC, National Physical Laboratory, England  (circa 1975).
*
*     Next, a more reliable test is performed on each component of the
*     gradient, for indices in the range  jverif(1)  thru  jverif(2).
*
*     lvlVer is the verify level, which has the following meaning:
*
*     -1         do not perform any check.
*      0         do the cheap test only.
*      1 or 3    do both cheap and full test on objective gradients.
*      2 or 3    do both cheap and full test on the Jacobian.
*
*     10 Oct 1998: First version based on combining s7chkJ and s7chkG.
*     27 Dec 2000: Permutations of the natural order included.
*     02 Aug 2003: snEXIT and snPRNT adopted.
*     20 Mar 2006: Relative error based on exact and estimated values.
*     ==================================================================
      external
     &     idamax, dasum, ddot, dnormi, s2ColN, s2RowN
      character
     &     key*4, str*120
      logical
     &     cheap, done, first, found, getCon, getObj
      integer
     &     i, idamax, imaxJ, ir, irN, j, jN, j1, j2, j3,
     &     j4, jfirst, jlast, jmaxG, jmaxJ, k, k1, k2, kmax, l, lvlVer,
     &     modefg, nFeas, nG, nJ, nKnown, Status, RightJ, RightG,
     &     WrongJ, WrongG, s2ColN, s2RowN
      double precision
     &     aGij, aGdif, dasum, ddot, dnormi, dx, emaxG, emaxJ, eps0,
     &     eps5, err, fdint1, fObj2, gdiff, gdummy, Gij, gj, gmaxG,
     &     gmaxJ, gp, xj, yj
*     ------------------------------------------------------------------
      double precision   zero,          one,          ok
      parameter         (zero = 0.0d+0, one = 1.0d+0, ok = 0.1d+0)
      logical            yes,           no
      parameter         (yes  = .true., no  = .false.)
      character          lWrong*4       , lRight*4
      data               lWrong /'bad?'/, lRight /'ok  '/
*     ------------------------------------------------------------------
      lvlVer    = iw( 78) ! Verify level

      if (lvlVer .lt. 0) return

      j1        = iw( 98) ! start col for obj.    derivative checking
      j2        = iw( 99) ! stop  col for obj.    derivative checking
      j3        = iw(100) ! start col for constr. derivative checking
      j4        = iw(101) ! stop  col for constr. derivative checking

      eps0      = rw(  2) ! eps**(4/5)
      eps5      = rw(  7) ! eps**(1/5)

      fdint1    = rw( 76) ! (1) forwrd diff. interval
      gdummy    = rw( 69) !

      iExit     = 0

      j1        = max( j1, 1 )
      j2        = min( j2, n )
      j3        = max( j3, 1 )
      j4        = min( j4, n )

      jFirst    = min( j1, j3    )
      jLast     = max( j2, j4    )

*     The problem functions are called to provide functions only.
*     Status indicates that there is nothing special about the call.

      modefg = 0
      Status = 0

*     cheap = do the cheap test only

      cheap     = lvlVer .eq. 0  .or.  jFirst .gt. jLast

*     ------------------------------------------------------------------
*     Cheap test.
*     ------------------------------------------------------------------
*     Generate a direction in which to perturb  x.

      call dcopy ( n, x, 1, x1, 1 )

      yj    = one/nnL
      do  j =  1, nnL
         y(j)  =   yj
         x1(j) =   yj           ! x1 = y
         yj    = - yj * 0.99999d+0
      end do

*     If needed, alter y to ensure that it will be a feasible direction.
*     If this gives zero, go back to original y and forget feasibility.

      dx     = fdint1 * (one + dnormi( nnL, x, 1 ))
      call s7chkp
     &   ( nnL, bl, bu, x, dx, y, nfeas )
      if (nfeas .eq. 0)
     &call dcopy
     &   ( nnL, x1, 1, y, 1 )

*     ------------------------------------------------------------------
*     Do not perturb x(j) if the jth column contains unknown elements.
*     ------------------------------------------------------------------
      nKnown = 0
      l      = 0
      do   j = 1, nnL

*        Do not perturb x(j) if g(j) is unknown.

         if (j .le. nnObj) then
            if (gObj(j) .eq. gdummy) y(j) = zero
         end if

         if (j .le. nnJac) then
            k1 = locJ(j)
            k2 = locJ(j+1) - 1

            do k  = k1, k2
               ir = indJ(k)
               if (ir .gt. nnCon) go to 130
               l  = l + 1
               if (gCon(l) .eq. gdummy) y(j) = zero
            end do
         end if
  130    if (y(j) .ne. zero) nKnown = nKnown + 1

      end do

      if (nKnown .gt. 0) then
         if ( cheap ) then
            call snPRNT( 11,
     &           ' Cheap test of user-supplied problem derivatives...',
     &           iw, leniw )
         else
            call snPRNT( 11,
     &           ' Verification of user-supplied problem derivatives.',
     &           iw, leniw )
         end if

*        ---------------------------------------------------------------
*        Compute functions at a short step along  y.
*        ---------------------------------------------------------------
         dx   = fdint1 * (one + dasum( nnL, x, 1 ))
         do j = 1, nnL
            x1(j) = x(j) + dx*y(j)
         end do

         getCon = nnJac .gt. 0
         getObj = nnObj .gt. 0

         call fgwrap
     &      ( iExit, modefg, Status, getCon, getObj,
     &        n, negCon, nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj,
     &        fgcon, fgobj, x1,
     &        ne, nlocJ, locJ, indJ,
     &        fCon2, fObj2, gCon2, gObj2,
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )
         if (iExit .ne. 0) go to 900
         
*        ---------------------------------------------------------------
*        Cheap test for the constraint Jacobian.
*        ---------------------------------------------------------------
*        Set   y1 = (fCon2 - fCon)/dx - gCon*y.  This should be small.

         if ( getCon ) then
            do i  = 1, nnCon
               y1(i)  = (fCon2(i) - fCon(i))/dx
            end do

            l     = 0
            do j  = 1, nnJac
               yj = y(j)

               k    = locJ(j)
               kmax = locJ(j+1) - 1
               done = no

*+             ---------------------------------------------------------
*+             while (k .le. kmax  .and.  .not. done) do
  180          if    (k .le. kmax  .and.  .not. done) then
                  ir     = indJ(k)
                  if (ir .gt. nnCon) then
                     done = yes
                  else
                     l      = l + 1
                     y1(ir) = y1(ir) - gCon(l)*yj
                  end if
                  k  = k + 1
                  go to 180
*+             end while
*+             ---------------------------------------------------------
               end if
            end do

            imaxJ  = idamax( nnCon, y1, 1 )
            gmaxJ  = (fCon2(imaxJ) - fCon(imaxJ))/dx
            emaxJ  = abs(y1(imaxJ) )/(one + abs( gmaxJ ))
            if (emaxJ .le. ok) then
               call snPRNT( 11,
     &              ' The constraint gradients seem to be OK.',
     &              iw, leniw )
            else
               call snPRNT( 11,
     &         ' XXX  The constraint gradients seem to be incorrect.',
     &              iw, leniw )
            end if
            write(str, 1601) emaxJ, n+s2RowN( imaxJ, leniw, iw )
            call snPRNT( 11, str, iw, leniw )
            call snPRNT(  1, ' ', iw, leniw )
         end if

*        ---------------------------------------------------------------
*        Cheap test for the objective gradient.
*        ---------------------------------------------------------------
         if ( getObj ) then
            gp     = ddot  ( nnObj, gObj, 1, y, 1 )
            gmaxG  = (fObj2 - fObj)/dx
            emaxG  = max( abs(gmaxG  - gp)/(one + abs( gmaxG )),
     &                    abs(gmaxG  - gp)/(one + abs( gp    )))

*           Set an error indicator if emaxG is too large.

            if (emaxG .le. ok) then
                call snPRNT( 11,
     &              ' The objective  gradients seem to be OK.',
     &              iw, leniw )
            else
               call snPRNT( 11,
     &         ' XXX  The objective  gradients seem to be incorrect.',
     &              iw, leniw )
            end if
            write(str, 1602) gp
            call snPRNT( 11, str, iw, leniw )
            write(str, 1603) gmaxG
            call snPRNT(  1, str, iw, leniw )
         end if
      end if

      if ( cheap ) go to 900

*     ------------------------------------------------------------------
*     Proceed with the verification of column elements.
*     Evaluate columns  jFirst thru jLast  of the problem derivatives.
*     ------------------------------------------------------------------
      call snPRNT( 11, ' ', iw, leniw )
      if (j3 .le. j4) then
         write(str, 2001)
         call snPRNT( 1, str, iw, leniw )
      else
         write(str, 2002)
         call snPRNT( 1, str, iw, leniw )
         call snPRNT( 1, ' ', iw, leniw )
      end if

      WrongJ =   0
      RightJ =   0
      WrongG =   0
      RightG =   0
      jmaxJ  =   0
      jmaxG  =   0
      nJ     =   0
      nG     =   0
      emaxJ  = - one
      emaxG  = - one

      do j  = 1, nnL
         jN = s2ColN( j, leniw, iw )

         getObj  = j .le. nnObj  .and. jN .ge. j1  .and.  jN .le. j2
         getCon  = j .le. nnJac  .and. jN .ge. j3  .and.  jN .le. j4

         if (getCon  .or.  getObj) then

*           See if there are any known gradients in this column.

            found  = no

            if ( getObj ) then
               if (gObj(j) .ne. gdummy) found = yes
            end if

            if ( getCon ) then

*              See if there are any known gradients in this column.

               l      = locG(j)
               k      = locJ(j)
               kmax   = locJ(j+1) - 1
               done   = no

*+             ---------------------------------------------------------
*+             while (k .le. kmax  .and. .not.(found  .or.  done)) do
  200          if    (k .le. kmax  .and. .not.(found  .or.  done)) then
                  ir = indJ(k)
                  if (ir .gt. nnCon) then
                     done = yes
                  else
                     if (gCon(l) .ne. gdummy) found = yes
                     l    = l + 1
                  end if
                  k  = k + 1
                  go to 200
               end if
*+             end while
*+             ---------------------------------------------------------
            end if

            if ( found ) then
*              ---------------------------------------------------------
*              gObj(j) or an element of the jth column of J is known.
*              ---------------------------------------------------------
               xj     = x(j)
               dx     = fdint1 * (one + abs( xj ))
               if (bl(j) .lt. bu(j)  .and.  xj .ge. bu(j)) dx = -dx
               x(j)   = xj + dx

               call fgwrap
     &            ( iExit, modefg, Status, getCon, getObj,
     &              n, negCon, nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj,
     &              fgcon, fgobj, x,
     &              ne, nlocJ, locJ, indJ,
     &              fCon2, fObj2, gCon2, gObj2,
     &              cu, lencu, iu, leniu, ru, lenru,
     &              cw, lencw, iw, leniw, rw, lenrw )
               if (iExit .ne. 0) go to 900

*              ---------------------------------------------------------
*              Estimate the jth column of the Jacobian.
*              Check the estimate against the user-supplied values.
*              Don't bother printing a line for an exact zero.
*              Look for nonzeros not in the sparse data structure.
*              ---------------------------------------------------------
               if ( getCon ) then
                  do i = 1, nnCon
                     y1(i) = (fCon2(i) - fCon(i))/dx
                  end do

                  l      = locG(j)
                  k      = locJ(j)
                  first  = yes
                  done   = no

*+                ------------------------------------------------------
*+                while (k .le. kmax  .and.  .not. done) do
  260             if    (k .le. kmax  .and.  .not. done) then

                     ir  = indJ(k)
                     irN = s2RowN( ir, leniw, iw )

                     if (ir .gt. nnCon) then
                        done = yes
                     else
                        Gij    = gCon(l)
                        aGij   = abs( Gij )

                        if (Gij .ne. gdummy) then
                           nJ     = nJ + 1
                           gdiff  = y1(ir)
                           aGdif  = abs(gdiff)
                           err    = max(abs(gdiff - Gij)/(one + aGij),
     &                                  abs(gdiff - Gij)/(one + aGdif))

                           if (emaxJ .lt. err) then
                              emaxJ  = err
                              imaxJ  = irN
                              jmaxJ  = jN
                           end if

                           if (err .le.  eps5) then
                              key    = lRight
                              RightJ = RightJ + 1
                           else
                              key    = lWrong
                              WrongJ = WrongJ + 1
                           end if

                           if (aGij + err .gt. eps0) then
                              if ( first ) then
                                 write(str, 2101) jN, xj, dx,
     &                                l, irN, Gij, gdiff, key
                                 call snPRNT( 11, str, iw, leniw )
                                 first  = no
                              else
                                 write(str, 2201)
     &                                l, irN, Gij, gdiff, key
                                 call snPRNT(  1, str, iw, leniw )
                              end if
                           end if
                        end if

                        l      = l + 1

*                       Mark this row as being in the sparse structure.

                        y1(ir) = gdummy

                     end if

                     k  = k + 1
                     go to 260
                  end if
*+                end while
*+                ------------------------------------------------------

*                 Check that column elements not in gCon  are zero.
*                 These are the unmarked elements of  y1.

                  do i = 1, nnCon
                     if ( y1(i) .ne. gdummy) then
                        irN = s2RowN( i, leniw, iw )
                        err = abs(y1(i))

                        if (err .gt. eps0) then
                           WrongJ = WrongJ + 1

                           if (first) then
                              write(str, 2301) jN, irN, y1(i), lWrong
                              call snPRNT( 1, str, iw, leniw )
                              first = no
                           else
                              write(str, 2302) irN, y1(i), lWrong
                              call snPRNT( 1, str, iw, leniw )
                           end if

                           if (emaxJ .lt. err) then
                              emaxJ  = err
                              imaxJ  = irN
                              jmaxJ  = jN
                           end if
                        end if
                     end if
                  end do
               end if ! getCon

*              ---------------------------------------------------------
*              Estimate gObj(j)
*              ---------------------------------------------------------
               if ( getObj ) then
                  gj    = gObj(j)

                  if (gj .ne. gdummy) then
                     nG    = nG + 1
                     gdiff = (fObj2 - fObj) / dx
                     err   = max( abs( gdiff - gj )/(one + abs( gj   )),
     &                            abs( gdiff - gj )/(one + abs( gdiff)))

                     if (err .gt. emaxG) then
                        emaxG  = err
                        jmaxG  = jN
                     end if

                     if (err .le. eps5) then
                        key    = lRight
                        RightG = RightG + 1
                     else
                        key    = lWrong
                        WrongG = WrongG + 1
                     end if

                     if (abs( gj ) + err  .gt.  eps0) then
                        if (j3 .le. j4) then
                           write(str, 2102) jN, xj, dx, gj, gdiff, key
                        else
                           write(str, 2103) jN, xj, dx, gj, gdiff, key
                        end if
                        call snPRNT( 1, str, iw, leniw )
                     end if
                  end if
               end if ! getObj
               x(j)  = xj
            end if ! found
         end if
      end do

*     ------------------------------------------------------------------
*     Final tally of the good, the bad and the ugly.
*     ------------------------------------------------------------------
      if (j3 .le. j4  .and.  nJ .gt. 0) then
         if (WrongJ .eq. 0) then
            write(str, 2501) RightJ, j3, j4
         else
            write(str, 2601) WrongJ, j3, j4
         end if
         call snPRNT( 11, str, iw, leniw )
         write(str, 2701) emaxJ, imaxJ, jmaxJ
         call snPRNT( 11, str, iw, leniw )
         call snPRNT(  1, ' ', iw, leniw )

         if (emaxJ .ge. one) then
            iExit = 52          ! Bad constraint gradients.
         end if
      end if

      if (j1 .le. j2  .and.  nG .gt. 0) then
         if (WrongG .eq. 0) then
            write(str, 2502) RightG, j1, j2
         else
            write(str, 2602) WrongG, j1, j2
         end if
         call snPRNT( 11, str, iw, leniw )
         write(str, 2702) emaxG, jmaxG
         call snPRNT( 11, str, iw, leniw )
         call snPRNT(  1, ' ', iw, leniw )

         if (emaxG .ge. one) then
            iExit = 51          ! Bad objective gradients.
         end if
      end if

*     ------------------------------------------------------------------
*     Print a message if we had to abandon the check.
*     ------------------------------------------------------------------
  900 if (iExit .lt. 0) then
         call snPRNT( 11,
     &      ' XXX  Unable to complete derivative check.',
     &              iw, leniw )
         iExit = 0
      end if

      return

 1601 format(' -->  The largest discrepancy was', 1p, e12.2,
     &       '  in constraint', i6)
 1602 format(' Gradient projected in one direction', 1p, e20.11)
 1603 format(' Difference approximation           ', 1p, e20.11)
 2001 format(' Column       x(j)        dx(j)', 3x,
     &   ' Element no.    Row        Derivative    Difference approxn')
 2002 format(6x, 'j', 7x, 'x(j)', 8x, 'dx(j)',
     &       11x, 'g(j)', 9x, 'Difference approxn')
 2101 format(i7, 1p, e16.8, e10.2, 2i10,               2e18.8, 2x, a4)
 2102 format(i7, 1p, e16.8, e10.2, 10x, ' Objective',  2e18.8, 2x, a4)
 2103 format(i7, 1p, e16.8, e10.2,                     2e18.8, 2x, a4)
 2201 format(      33x, 2i10, 1pe18.8, e18.8, 2x, a4)
 2301 format(i7, 2x, 'Nonzero not in sparse structure ', '??', 4x, i6,
     &          18x, 1p, e18.8, 2x, a4 )
 2302 format(    9x, 'Nonzero not in sparse structure ', '??', 4x, i6,
     &          18x, 1p, e18.8, 2x, a4 )
 2501 format(i7, '  Jacobian elements in cols ', i6, '  thru', i6,
     &           '  seem to be OK.')
 2502 format(i7, '  objective gradients out of', i6, '  thru', i6,
     &           '  seem to be OK.')
 2601 format(' XXX  There seem to be', i6,
     &       '  incorrect Jacobian elements in cols', i6, '  thru', i6)
 2602 format(' XXX  There seem to be', i6,
     &       '  incorrect objective gradients in cols', i6, '  thru',i6)
 2701 format(' -->  The largest relative error was', 1p, e12.2,
     &       '   in row', i6, ',  column', i6)
 2702 format(' -->  The largest relative error was', 1p, e12.2,
     &       '   in column', i6)

      end ! subroutine s7chkG

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s7chkp
     &   ( n, bl, bu, x, dx, p, nFeas )

      implicit
     &     none
      integer
     &     n, nFeas
      double precision
     &     dx, bl(n), bu(n), x(n), p(n)

*     ==================================================================
*     s7chkp  checks that x + dx*p is feasible.
*     It is used by s7chkG for the cheap gradient checks.
*
*     Original:    Looked at the sign of p for variables on a bound.
*     13 Mar 1992: dx added as a parameter to make certain that
*                  x + dx*p does not lie outside the bounds.
*                  p may be altered to achieve this.
*     09 Aug 1992: First version based on Minos routine m7chkg.
*     27 Feb 2000: Current version.
*     ==================================================================
      integer
     &     j
      double precision
     &     b1, b2, pj, xj, xnew
*     ------------------------------------------------------------------
      double precision   zero
      parameter        ( zero  = 0.0d+0 )
*     ------------------------------------------------------------------
      nFeas  = 0
      do  j = 1, n
         xj = x(j)
         b1 = bl(j)
         b2 = bu(j)
         if (b1   .eq. b2  ) p(j) = zero

         if (p(j) .ne. zero) then

*           x(j) is not fixed, so there is room to move.
*           If xj + dx*pj is beyond one bound, reverse pj
*           and make sure it is not beyond the other.
*           Give up and use set pj = zero if both bounds are too close.

            pj     = p(j)
            xnew   = xj  +  dx*pj

            if (pj .gt. zero) then
               if (xnew .gt. b2) then
                  pj     = - pj
                  xnew   =   xj  +  dx*pj
                  if (xnew .lt. b1) pj = zero
               end if
            else
               if (xnew .lt. b1) then
                  pj     = - pj
                  xnew   =   xj  +  dx*pj
                  if (xnew .gt. b2) pj = zero
               end if
            end if

            p(j)   = pj
            if (pj .ne. zero) nFeas  = nFeas + 1
         end if
      end do

      end ! subroutine s7chkp

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s7FixX
     &   ( n, nFixed, ptrbAll, bndptrb, bl, bu, x )

      implicit
     &     none
      logical
     &     ptrbAll
      integer
     &     n, nFixed
      double precision
     &     bndptrb, bl(n), bu(n), x(n)

*     ==================================================================
*     s7FixX  ensures that variables satisfy their simple bounds.
*
*     if ptrbAll = true, the bounds for fixed variables are relaxed by
*     bndptrb.
*
*     On exit, nFixed = the number of fixed variables.
*
*     27 Sep 2003: First version of s7FixX.
*     27 Sep 2003: Current version.
*     ==================================================================
      integer
     &     j
      double precision
     &     absxj
*     ------------------------------------------------------------------
      double precision   one
      parameter         (one = 1.0d+0)
*     ------------------------------------------------------------------
      nFixed = 0
      do j = 1, n
         x(j)  = max( x(j), bl(j) )
         x(j)  = min( x(j), bu(j) )
         absxj = one + abs(x(j))
         if (bl(j) .eq. bu(j)) then
            nFixed = nFixed + 1
            if ( ptrbAll ) then
               bl(j) = bl(j) - bndptrb*absxj
               bu(j) = bu(j) + bndptrb*absxj
            end if
         end if
      end do

      end ! subroutine s7FixX

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s7Jac
     &   ( iExit, Status, userfg, ptrbAll,
     &     nF, n, InfBnd, imaxJ, emaxJ,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     rowtyp, coltyp,
     &     x, xlow, xupp, F, G, w, y, z,
     &     Fw, Fy, Fz, Gcolw, Gcoly, Gcolz,
     &     iw, leniw, cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      external
     &     userfg
      logical
     &     ptrbAll
      integer
     &     iExit, imaxJ, nF, n, neA, lenA, neG, lenG, leniw,
     &     lencu, leniu, lenru, Status, iAfun(lenA), jAvar(lenA),
     &     iGfun(lenG), jGvar(lenG), rowtyp(nF), coltyp(n), iw(leniw),
     &     iu(leniu)
      double precision
     &     InfBnd, emaxJ, A(lenA), F(nF), Fw(nF), Fy(nF), Fz(nF),
     &     G(lenG), Gcolw(nF), Gcoly(nF), Gcolz(nF), x(n), xlow(n),
     &     xupp(n), w(n), y(n), z(n), ru(lenru)
      character
     &     cu(lencu)*8

*     ==================================================================
*     s7Jac   does the work for snJac
*
*     26 Oct 2002: First version.
*     27 Sep 2003: More thorough check for feasibility.
*     22 Apr 2007: Arguments of s6usrf changed.
*     ==================================================================
      external
     &     idamax, dnormi
      integer
     &     i, idamax, j, k, k0, neA0, needF, needG, nFixed,
     &     seeds(3)
      double precision
     &     damper, dnormi, dwj, dyj0, dyj, dzj, eps0, eps, fdint1,
     &     gmaxJ, Jwij, Jyij, Jzij, wj, yj, zj, s1eps, tol, xnorm, ynorm
*     ------------------------------------------------------------------
      double precision   zero,           one,             two
      parameter         (zero  = 0.0d+0, one    = 1.0d+0, two = 2.0d+0)
      double precision   three
      parameter         (three = 3.0d+0)
      integer            Linear,     Nonlin
      parameter         (Linear = 0, Nonlin = 1 )
*     ------------------------------------------------------------------
      eps    = s1eps ( )        ! machine precision
      eps0   = eps**(0.80d+0)
      fdint1 = sqrt(eps)

      xnorm  = max ( one, dnormi( n, x, 1 ) )
      tol    = eps0*xnorm

      neG    = 0
      neA    = 0
      imaxJ  = 0
      emaxJ  = zero

*     Move x inside its bounds.
*     If requested, perturb bounds for fixed variables.

      call s7FixX
     &   ( n, nFixed, ptrbAll, one, xlow, xupp, x )

      needF  = 1
      needG  = 0
      Status = 1
      call userfg
     &   ( Status, n, x,
     &     needF, nF, F,
     &     needG, lenG, G,
     &     cu, lencu, iu, leniu, ru, lenru )
      if (Status .lt. 0) go to 999

*     Define random points  w, y and z near the user-supplied x.

      seeds(1) = 5872
      seeds(2) = 8584
      seeds(3) = 4879

*     Create a random feasible perturbation of the variables.

      call s7pert
     &   ( n, seeds, x, xlow, xupp, xnorm, w )
      call s7pert
     &   ( n, seeds, x, xlow, xupp, xnorm, y )
      call s7pert
     &   ( n, seeds, x, xlow, xupp, xnorm, z )

      call iload ( n , Linear, coltyp, 1 )
      call iload ( nF, Linear, rowtyp, 1 )

*     Evaluate F at  w, y and z.

      Status = 0
      call s6usrf
     &   ( Status, n, x, w, damper, userfg,
     &     needF, nF , Fw,
     &     needG, lenG, G,
     &     cu, lencu, iu, leniu, ru, lenru )
      if (Status .lt. 0) go to 999

      call s6usrf
     &   ( Status, n, x, y, damper, userfg,
     &     needF, nF, Fy,
     &     needG, lenG, G,
     &     cu, lencu, iu, leniu, ru, lenru )
      if (Status .lt. 0) go to 999

      call s6usrf
     &   ( Status, n, x, z, damper, userfg,
     &     needF, nF, Fz,
     &     needG, lenG, G,
     &     cu, lencu, iu, leniu, ru, lenru )
      if (Status .lt. 0) go to 999

      do   j = 1, n
*        ---------------------------------------------------------------
*        Compute functions at  w + dwj*e_j, y + dyj*e_j  and  z + dzj*e_j.
*        ---------------------------------------------------------------
         wj   = w(j)
         yj   = y(j)
         zj   = z(j)

*        Evaluate F at w + dwj*e_j.

         call s7step
     &      ( InfBnd, xlow(j), xupp(j), wj, abs(wj), dwj )

         if (dwj .eq. zero) then
            call dload ( nF, zero, Gcolw, 1 )
         else
            w(j) = wj + dwj

            call userfg
     &         ( Status, n, w,
     &           needF, nF, F,
     &           needG, lenG, G,
     &           cu, lencu, iu, leniu, ru, lenru )
            if (Status .lt. 0) go to 999

            do i  = 1, nF
               Gcolw(i) = (F(i) - Fw(i))/dwj
            end do
         end if

*        Evaluate F at y + dyj*e_j.

         call s7step
     &      ( InfBnd, xlow(j), xupp(j), yj, two*abs(yj), dyj )

         if (dyj .eq. zero) then
            call dload ( nF, zero, Gcoly, 1 )
         else
            y(j) = yj + dyj
            call userfg
     &         ( Status, n, y,
     &           needF, nF, F,
     &           needG, lenG, G,
     &           cu, lencu, iu, leniu, ru, lenru )
            if (Status .lt. 0) go to 999

            do i  = 1, nF
               Gcoly(i) = (F(i) - Fy(i))/dyj
            end do
         end if

*        Evaluate F at z + dzj*e_j.

         call s7step
     &      ( InfBnd, xlow(j), xupp(j), zj, three*abs(zj), dzj )

         if (dzj .eq. zero) then
            call dload ( nF, zero, Gcolz, 1 )
         else
            z(j) = zj + dzj
            call userfg
     &         ( Status, n, z,
     &           needF, nF, F,
     &           needG, lenG, G,
     &           cu, lencu, iu, leniu, ru, lenru )
            if (Status .lt. 0) go to 999

            do i  = 1, nF
               Gcolz(i) = (F(i) - Fz(i))/dzj
            end do
         end if

*        ---------------------------------------------------------------
*        Figure out which elements derivatives are zero or constant
*        ---------------------------------------------------------------
         do    i = 1, nF
            Jwij = Gcolw(i)
            Jyij = Gcoly(i)
            Jzij = Gcolz(i)

            if ((abs(Jwij) + abs(Jyij) + abs(Jzij)) .gt. tol) then

*              Nonzero Jacobian element, to be stored in A or G.

               if (abs(Jwij - Jyij) .le. tol .and.
     &             abs(Jwij - Jzij) .le. tol      ) then

*                 Constant  Jij

                  neA = neA + 1
                  if (neA .gt. lenA) then
                     call snPRNT( 1, ' Increase lenA.', iw, leniw )
                     iExit = 91
                     go to 999
                  end if

                  iAfun(neA) = i
                  jAvar(neA) = j
                  A(neA)     = (Jwij + Jzij)/two
               else

*                 Nonlinear Jij

                  neG = neG + 1
                  if (neG .gt. lenG) then
                     call snPRNT( 1, ' Increase lenG.', iw, leniw )
                     iExit = 91
                     go to 999
                  end if

                  iGfun(neG) = i
                  jGvar(neG) = j
                  rowtyp(i)  = Nonlin
                  coltyp(j)  = Nonlin
               end if
            end if
         end do

         w(j)  = wj
         y(j)  = yj
         z(j)  = zj
      end do

*     -----------------------------------------------------------------
*     Constant elements in nonlinear rows and columns must be treated
*     as nonlinear.
*     -----------------------------------------------------------------
      k    = 1
      neA0 = neA

      do k0 = 1, neA0
         i  = iAfun(k0)
         j  = jAvar(k0)
         if (rowtyp(i) .eq. Nonlin  .and.  coltyp(j) .eq. Nonlin) then
            neA = neA - 1
            neG = neG + 1
            if (neG .gt. lenG) then
               call snPRNT( 1, ' Increase lenG.', iw, leniw )
               iExit = 91
               go to 999
            end if
            iGfun(neG) = i
            jGvar(neG) = j
         else
            if (k .lt. k0) then
               iAfun(k) = iAfun(k0)
               jAvar(k) = jAvar(k0)
               A(k)     = A(k0)
            end if
            k = k + 1
         end if
      end do

      if (neG .lt. n*nF) then
*        ===============================================================
*        J has some constant elements.  Better check we got it right.
*
*        Compare J*p with (G + A)*p, everything computed at y
*        ===============================================================
         do j = 1, n
*           ------------------------------------------------------------
*           Compute functions at  y + dyj*e_j.
*           ------------------------------------------------------------
            yj   = y(j)

*           Evaluate F at y + dyj*e_j.

            dyj0 = abs(yj)*fdint1

            call s7step
     &         ( InfBnd, xlow(j), xupp(j), yj, dyj0, dyj )

            if (dyj .eq. zero) then
               call dload ( nF, zero, Gcoly, 1 )
            else
               y(j) = yj + dyj
               call userfg
     &            ( Status, n, y,
     &              needF, nF, F,
     &              needG, lenG, G,
     &              cu, lencu, iu, leniu, ru, lenru )
               if (Status .lt. 0) go to 999

               do i  = 1, nF
                  Gcoly(i) = (F(i) - Fy(i))/dyj
               end do
            end if

            do k = 1, neG
               if (jGvar(k) .eq. j) then
                  i    = iGfun(k)
                  G(k) = Gcoly(i)
               end if
            end do

            y(j)  = yj
         end do

*        ---------------------------------------------------------------
*        Compute a new random feasible perturbation z and compute
*        the functions at  w = y + z.
*        ---------------------------------------------------------------
         ynorm = (one + dnormi( n, y, 1 ))

         call s7pert
     &      ( n, seeds, y, xlow, xupp, fdint1*ynorm, w )

*        Evaluate F at w + z.

         do j = 1, n
            z(j) = w(j) - y(j)
         end do

         call userfg
     &      ( Status, n, w,
     &        needF, nF, F,
     &        needG, lenG, G,
     &        cu, lencu, iu, leniu, ru, lenru )
         if (Status .lt. 0) go to 999

*        ---------------------------------------------------------------
*        Cheap test for the  Jacobian.
*        ---------------------------------------------------------------
*        Compute  F - (Fy + (G + A)*z).  This should be small.

         do i = 1, nF
            Fw(i) = F(i) - Fy(i)
         end do

         do k = 1, neG
            i = iGfun(k)
            j = jGvar(k)
            Fw(i) = Fw(i) - G(k)*z(j)
         end do

         do k = 1, neA
            i = iAfun(k)
            j = jAvar(k)
            Fw(i) = Fw(i) - A(k)*z(j)
         end do

         imaxJ = idamax( nF, Fw, 1 )
         gmaxJ = F(imaxJ) - Fy(imaxJ)
         emaxJ = abs( Fw(imaxJ) )/(one + abs( gmaxJ ))
      end if

      return

  999 if      (status .eq. -1) then
         iExit = 63            ! unable to proceed into undefined region
      else if (status .le. -2) then
         iExit = 71            ! terminated during function evaluation
      end if

      end ! subroutine s7Jac

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s7pert
     &   ( n, seeds, x, bl, bu, order, xPert )

      implicit
     &     none
      integer
     &     n, seeds(3)
      double precision
     &     x(n), xPert(n), bl(n), bu(n), order

*     ==================================================================
*     s7pert finds a feasible perturbation xPert of a feasible point x.
*     The order of the perturbation will be at most "order" and the
*     elements of xPert  will lie between bl and bu for bl < bu.
*
*     seeds(1:3) are seed values for the random number generator.
*
*     26 Oct 2002: First version based on snadiopt routine
*     27 Sep 2003: Current version of s7pert.
*     ==================================================================
      integer
     &     j, dirctn
      double precision
     &     theta
*     ------------------------------------------------------------------
      integer            down,      fixed,     up,     free
      parameter         (down = -1, fixed = 0, up = 1, free = 2)
      double precision   half         , one
      parameter         (half = 0.5d+0, one   = 1.0d+0)
      double precision   two,           three,          four
      parameter         (two  = 2.0d+0, three = 3.0d+0, four = 4.0d+0)
*     ------------------------------------------------------------------

*     Get a vector of pseudo-random numbers between 0 and 1

      call ddrand( n, xPert, 1, seeds )

      do j = 1, n

*        If dirctn .eq. free, then x(j) is not at one of it's bounds.

         dirctn = free
         if (x(j) .eq. bl(j) ) then

*           x(j) is on its lower bound, we must perturb up.

            dirctn = up
         end if

         if (x(j) .eq. bu(j) ) then

*           x(j) is on its upper bound...

            if ( dirctn .eq. free ) then

*              ... but not on its lower bound.  Perturb down.

               dirctn = down
            else

*              ... and is also on its lower bound, and so it must be
*              a fixed variable.  It is not perturbed.

               dirctn = fixed
            end if
         end if

         if ( dirctn .eq. free) then

*           If the variable is not on either bound, choose randomly
*           whether to go up or down.

            if (xPert(j) .gt. half) then
               xPert(j) = two*(xPert(j) - half)
               dirctn = up
            else
               xPert(j) = two* xPert(j)
               dirctn = down
            end if
         end if

         theta = half * xPert(j)
         if ( dirctn .eq. down) then

*           Perturb between 1/8 and 1/4 of the way to the lower bound,
*           but no more than "order".

            xPert(j) = theta*x(j) + (one - theta)*
     &                 max( (three*x(j) + bl(j))/four, x(j) - order)
         else if (dirctn .eq. up) then

*           Perturb between 1/8 and 1/4 of the way to the upper bound,
*           but no more than "order".

            xPert(j) = theta*x(j) + (one - theta)*
     &                 min( (three*x(j) + bu(j))/four, x(j) + order)
         else

*           dirctn .eq. fixed

*           xPert(j) = theta*x(j) + (one - theta)*order
            xPert(j) = x(j)
         end if
      end do

      end ! subroutine s7pert

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s7step
     &   ( InfBnd, blj, buj, xj, dxj0, dxj )

      implicit
     &     none
      double precision
     &     InfBnd, blj, buj, dxj0, dxj, xj

*     ==================================================================
*     s7step finds a feasible direction of perturbation dxj for xj.
*
*     dxj0  is a suggested value for the step.
*
*     26 Oct 2002: First version.
*     27 Sep 2003: Current version of s7step.
*     ==================================================================
      integer
     &     move
*     ------------------------------------------------------------------
      integer            down,      fixed,     up,     free
      parameter         (down = -1, fixed = 0, up = 1, free = 2)
      double precision   zero,          half
      parameter         (zero = 0.0d+0, half = 0.5d+0)
*     ------------------------------------------------------------------
*     If move .eq. free, then xj is not at one of it's bounds.

      move = free
      if (xj .eq. blj) then

*        xj is on its lower bound, we must perturb up.

         move = up
      end if

      if (xj .eq. buj) then

*        xj is on its upper bound...

         if (move .eq. free) then

*           ... but not on its lower bound.  Perturb down.

            move = down
         else

*           ... and is also on its lower bound, and so it must be
*           a fixed variable.  It is not perturbed.

            move = fixed
         end if
      end if

      if (move .eq. free) then

*        If the variable is not on either bound, choose randomly
*        whether to go up or down.

         if (xj .gt. half) then
            move =  up
         else
            move =  down
         end if
      end if

      dxj   =  dxj0

      if (move .eq. up) then
         if (buj .lt.  InfBnd) then
            dxj = min( dxj, buj - xj )
         end if

      else if (move .eq. down) then
         dxj = - dxj
         if (blj .gt. -InfBnd) then
            dxj = max( dxj, blj - xj )
         end if

      else if (move .eq. fixed) then
         dxj = zero
      end if

      end ! subroutine s7step

