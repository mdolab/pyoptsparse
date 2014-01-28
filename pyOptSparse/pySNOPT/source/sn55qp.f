!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     file  sn55qp.f
!
!     s5chkp   s5chzq   s5getp   s5getR   s5Hfac   s5Hz     s5QP
!     s5QPfg   s5QPit   s5Rcol   s5rg     s5Rsng   s5Sdel   s5Zp
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5chkp
     &   ( iExit, itn, nBS, jqSave, kBS, gtp, pBS, iw, leniw )

      implicit
     &     none
      integer
     &     iExit, itn, nBS, jqSave, leniw, kBS(nBS), iw(leniw)
      double precision
     &     gtp, pBS(nBS)

!     ==================================================================
!     s5chkp  makes  pBS  a feasible direction.
!
!     16 Jun 1995: First version of s5chkp.
!     02 Aug 2003: snPRNT adopted.
!     02 Aug 2003: Current version of s5chkp.
!     ==================================================================
      character
     &     str*80
      integer
     &     kSave, j, jq, k
      double precision
     &     pSave
!     ------------------------------------------------------------------
      double precision   zero,          one
      parameter         (zero = 0.0d+0, one = 1.0d+0)
!     ------------------------------------------------------------------
      iExit = 0

!     ------------------------------------------------------------------
!     Find the element of  pBS  corresponding to the most recently freed
!     variable. Usually, it will be pBS(nBS).
!     ------------------------------------------------------------------
      jq    = abs(jqSave)

      kSave = 0
      do k =  nBS, 1, -1
         j = kBS(k)
         if (j .eq. jq) then
            kSave = k
            go to 100
         end if
      end do

!     ------------------------------------------------------------------
!     Choose the sign of  pBS  so that the most recently freed
!     variable continues to increase or decrease.
!     ------------------------------------------------------------------
  100 if (kSave .gt. 0) then
         pSave = pBS(kSave)

         if (jqSave .lt. 0  .and.  pSave .gt. zero  .or.
     &       jqSave .gt. 0  .and.  pSave .lt. zero      ) then
            call dscal ( nBS, (-one), pBS, 1 )
            gtp  = - gtp
         end if

         if (gtp .gt. zero) then
!           ------------------------------------------------------------
!           Looks as though the sign of gtp cannot be relied upon.
!           In later versions we'll fix this variable.
!           For now, we just print a warning and stop.
!           ------------------------------------------------------------
            write(str, 1000) itn, gtp
            call snPRNT( 23, str, iw, leniw )
            iExit = 1           ! Bad directional derivative
         end if
      else
!        ---------------------------------------------------------------
!        Couldn't find the index of the most recently freed variable.
!        This should never happen!
!        ---------------------------------------------------------------
         write(str, 9000) jqSave
         call snPRNT( 23, str, iw, leniw )
      end if

      return

 1000 format(' Itn', i7, ': Bad directional derivative ', 1p, e9.1 )
 9000 format(' XXX  s5chkp.  kSave not found. jqSave = ', i5 )

      end ! subroutine s5chkp

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5chzq
     &   ( m, mBS, n, nb, nS, kBSq, pivot, tolpiv,
     &     ne, nlocA, locA, indA, Acol,
     &     kBS, bl, bu, xBS, y, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     m, mBS, n, nb, ne, nlocA, nS, kBSq, leniw, lenrw,
     &     locA(nlocA), indA(ne), kBS(mBS), iw(leniw)
      double precision
     &     pivot, tolpiv, Acol(ne), bl(nb), bu(nb), xBS(mBS), y(mBS),
     &     rw(lenrw)

!     ==================================================================
!     s5chzq  selects a superbasic to replace the kp-th basic variable.
!     On entry,  y  contains the kp-th row of B(inverse).
!     On exit, pivot and  y(m+1), ..., y(m+nS) define the S-part of
!     the modifying vector w.
!
!     01 Dec 1991: First version based on Minos routine m7chzq.
!     02 Aug 2003: snPRNT adopted.
!     30 Jun 2005: Current version of s5chzq.
!     ==================================================================
      character
     &     str*80
      integer
     &     j, k, m1, idamax
      double precision
     &     d1, d2, dpiv, eps0, tol, xj
!     ------------------------------------------------------------------
      double precision   zero,          point1,          one
      parameter         (zero = 0.0d+0, point1 = 0.1d+0, one = 1.0d+0)
      integer            Transp
      parameter         (Transp = 1)
!     ------------------------------------------------------------------
      eps0      = rw(  2) ! eps**(4/5)
*     tolpiv    = rw( 60) ! excludes small elements of y.

!     Set yS = 0 -  S'*y.

      m1        = m  + 1
      call s2Bprd
     &   ( Transp, eps0, n, nS, kBS(m1),
     &     ne, nlocA, locA, indA, Acol,
     &     (-one), y, m, zero, y(m1), nS )

      kBSq   = m  +  idamax( nS, y(m1), 1 )
      pivot  = abs( y(kBSq) )

!     Exit if the pivot is too small.

      if (pivot .lt. tolpiv) then
         write(str, 1000)  pivot
         call snPRNT( 31, str, iw, leniw )
         kBSq   = - (m + nS)
      else

!        Choose one away from its bounds if possible.

         tol    =   point1*pivot
         dpiv   = - one

         do k = m1, m+nS
            if (abs( y(k) ) .ge. tol) then
               j     = kBS(k)
               xj    = xBS(k)
               d1    = xj    - bl(j)
               d2    = bu(j) - xj
               d1    = min( abs( d1 ), abs( d2 ) )
               if (dpiv .le. d1) then
                  dpiv  = d1
                  kBSq  = k
               end if
            end if
         end do

         pivot   = - y(kBSq)

      end if ! pivot .ge. tolpiv

      return

 1000 format(' XXX  s5chzq.  Max pivot is too small:', 1p, e11.1 )

      end ! subroutine s5chzq

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5getp
     &   ( feasbl, gotR, newSB, PosDef,
     &     maxR, lenR, n, djq, R, g, p, gp, pHp )

      implicit
     &     none
      logical
     &     gotR, newSB, feasbl, PosDef
      integer
     &     maxR, lenR, n
      double precision
     &     djq, gp, pHp, R(lenR), g(n), p(n)

!     ==================================================================
!     s5getp  computes a search direction  p  for the superbasic
!     variables, using the current reduced gradient  g.
!
!     29 Mar 2001: R stored by rows.
!     20 May 2001: Current version.
!     ==================================================================
      external
     &     ddot
      integer
     &     ldiag
      double precision
     &     dirctn, diag, ddot
!     ------------------------------------------------------------------
      integer            WithR,      WithRt
      parameter         (WithR  = 0, WithRt = 1)
      double precision   zero,          one
      parameter         (zero = 0.0d+0, one = 1.0d+0)
!     ------------------------------------------------------------------

      if (feasbl  .and.  gotR) then
         if ( PosDef ) then
            ! Compute the Newton direction.

            call dcopy
     &         ( n,         g, 1, p, 1 )
            call s6Rsol
     &         ( WithRt, maxR, n, lenR, R, p )
            pHp  = ddot  ( n, p, 1, p, 1 )
            call s6Rsol
     &         ( WithR , maxR, n, lenR, R, p )
         else
            ! Compute a direction of zero or negative curvature.

            ldiag  = (n - 1)*maxR + (3 - n)*n/2  ! Magic formula!
            diag   = R(ldiag)                    ! Save last diag of R.
            R(ldiag) = one
            pHp    = diag**2
            call dload ( n, zero, p, 1 )

            if ( newSB ) then
               dirctn = g(n)
            else
               dirctn = djq
            end if

            if (dirctn .ge. zero) then
               p(n) =   one
            else
               p(n) = - one
            end if

            call s6Rsol
     &         ( WithR , maxR, n, lenR, R, p )
            R(ldiag) = diag                      ! Restore diag of R.
         end if
      else
         !--------------------------------------------------------------
         ! Direction of steepest-descent.
         !--------------------------------------------------------------
         call dcopy ( n, g, 1, p, 1 )
         pHp = zero
      end if ! feasbl and gotR

      !-----------------------------------------------------------------
      ! Fix the sign of p.
      ! ------------------------------------------------------------------
      call dscal ( n, (-one), p, 1 )
      gp  = ddot  ( n, g, 1, p, 1 )

      end ! subroutine s5getp

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5getR
     &   ( iExit, Hprod, Hprod1,
     &     Hcalls, gotR, typeLU, LUreq,
     &     itn, lenR, m, mBS, n, nb,
     &     nnH, nS, PrtLvl, minimz, iObj, targtH, targtZ,
     &     ne, nlocA, locA, indA, Acol,
     &     hs, kBS, bl, bu, blBS, buBS, R,
     &     nrhs0, nrhs, rhs,
     &     x, xBS, iy, iy1, y, y1, y2,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1
      logical
     &     gotR
      integer
     &     Hcalls, iExit, iObj, itn, lenR, lencu, leniu, lenru, lencw,
     &     leniw, lenrw, LUreq, m, mBS, n, nb, ne, nlocA, nnH, nrhs0,
     &     nrhs, nS, PrtLvl, minimz, typeLU, locA(nlocA), indA(ne),
     &     hs(nb), kBS(mBS), iy(nb), iy1(nb), iu(leniu), iw(leniw)
      double precision
     &     targtH, targtZ, Acol(ne), bl(nb), bu(nb), blBS(mBS),
     &     buBS(mBS), rhs(nrhs0), R(lenR), xBS(mBS), x(nb),
     &     y(nb), y1(nb), y2(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s5getR   computes the Cholesky factor of the reduced Hessian.
!
!     On entry, the LU factorization is assumed to be known.
!       gotR = .false.
!
!     iExit      Result
!     -----      ------
!       0        the reduced Hessian has been computed successfully
!      >0        fatal error
!
!     LUreq =  1  Frequency
!     LUreq =  2  LU nonzeros increased
!     LUreq =  3
!     LUreq =  4
!     LUreq =  5  Singular after LU mod
!     LUreq =  6  Unstable LU mod (growth in new column of U)
!     LUreq =  7  Not enough memory
!     LUreq =  8
!     LUreq =  9
!     LUreq = 10  Row error in setx
!     LUreq = 11  Big  dx   in setx
!
!     LUreq = 20
!     LUreq = 21  Iterative refinement failed in QP
!     LUreq = 22  Unbounded QP
!     LUreq = 23
!     LUreq = 24  Small directional derivative in QP
!     LUreq = 25  Ill-conditioned null-space basis in QP
!     LUreq = 26  Indefinite Z'HZ in QP
!     LUreq = 27  R singular after bound swap in QP
!
!     25 Oct 2003: First version of s5getR based on s8getR.
!     26 Dec 2003: s2newLU added.
!     09 Dec 2004: Changed to column-packed format for H.
!     09 Dec 2004: Current version of s5getR.
!     ==================================================================
      character
     &     str*80
      logical
     &     LUok, needLU, newB, newLU, Rcheck
      integer
     &     condZ, eigH, inform, maxR, nSwap
      double precision
     &     eps, flmax, Hdmax, plInfy
!     ------------------------------------------------------------------
      parameter         (condZ  = 192) ! condition estimate of Z
      double precision   one
      parameter         (one    = 1.0d+0)
!     ------------------------------------------------------------------
      maxR      = iw( 52) ! max columns of R
      eigH      = iw(200) ! =1(0) for definite QP Hessian

      eps       = rw(  1) ! machine precision.  IEEE DP  2.22e-16
      flmax     = rw(  8) ! est. of the largest pos. real

      iExit     = 0
      plInfy    = flmax
      LUok      = .true.

!     ==================================================================
!+    while (LUok  .and. .not. gotR) do
  100 if    (LUok  .and. .not. gotR) then
!     ------------------------------------------------------------------
         inform = 0
         needLU = LUreq .gt. 0

         if ( needLU ) then
            call s2Bfac
     &         ( iExit, typeLU, needLU, newLU, newB,
     &           iObj, itn, PrtLvl, LUreq,
     &           m, mBS, n, nb, nnH, nS, nSwap,
     &           ne, nlocA, locA, indA, Acol,
     &           kBS, hs, bl, bu, blBS, buBS,
     &           nrhs0, nrhs, rhs, x, xBS,
     &           iy, iy1, y, y2, iw, leniw, rw, lenrw )
            if (iExit .ne. 0) goto 900
         end if

         if (nS .gt. 0) then
!           ------------------------------------------------------------
!           Compute and factorize  Z'HZ.
!           ------------------------------------------------------------
            call s5Hz
     &         ( inform, Hprod, Hprod1, maxR, lenR,
     &           minimz, m, mBS, n, nb, nnH, nS, Hcalls,
     &           ne , nlocA, locA, indA, Acol,
     &           Hdmax, rw(condZ), targtZ, kBS, R, y, y1, y2,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )

!           Check for trouble in s5Hz.  Possible exits are:
!           inform  Status
!           ------  ------
!            -1     Ill-conditioned null-space basis
!             0     reduced Hessian computed successfully
!            >0     Fatal error in LU solve

            if (inform .eq. 0) then
               call s5Hfac
     &            ( inform, eigH, itn,
     &              lenR, m, maxR, mBS, nb, nS,
     &              targtH, Hdmax, hs, kBS, iy,
     &              bl, bu, blBS, buBS, x, xBS, R,
     &              iw, leniw, rw, lenrw )

!              Check for trouble in s5Hfac.
!              inform    Status
!              ------    ------
!                -2      H singular (but should be positive definite)
!                -1      H indefinite
!                 0      normal exit

               if (inform .ne. 0) then

!                 The reduced Hessian appears to be indefinite.
!                 Refactorize B with reduced factor tol.
!                 If the factor tol is already tight, give up.

                  write(str, 1200) itn
                  call snPRNT( 23, str, iw, leniw )
                  call s2tryLU
     &               ( itn, 26, nS, LUreq, LUok, typeLU,
     &                 iw, leniw, rw, lenrw )
                  if (.not. LUok) then
                     targtH = one/(eps*eps)
                     LUok   = .true.
                  end if
               end if

               Rcheck = .false.
               if (Rcheck) then
                  call s6Rchk
     &               ( iExit, Hprod, Hprod1, itn, minimz,
     &                 maxR, lenR, m, mBS, n, nb, Hcalls, nnH, nS,
     &                 ne, nlocA, locA, indA, Acol,
     &                 kBS, R, y, y1, y2,
     &                 cu, lencu, iu, leniu, ru, lenru,
     &                 cw, lencw, iw, leniw, rw, lenrw )
                  if (iExit .ne. 0) go to 900
               end if

            else if (inform .eq. -1) then

!              Ill-conditioned Z from s5Hz.
!              Refactorize B, possibly with a reduced factor tol.
!              If factor tol is already tight, accept Z, however bad.

               write(str, 1300) itn, rw(condZ)
               call snPRNT( 23, str, iw, leniw )
               call s2tryLU
     &            ( itn, 25, nS, LUreq, LUok, typeLU,
     &              iw, leniw, rw, lenrw )
               if (.not. LUok) then
                  targtZ = plInfy
                  LUok   = .true.
               end if

            else if (inform .gt. 0) then ! LU error in s5Hz
               iExit = inform
               go to 900
            end if
         end if

         gotR = inform .eq. 0

         go to 100
      end if
!+    end while
!     ------------------------------------------------------------------

      if (.not. gotR) then      ! indefinite reduced Hessian
         iExit = 94
      end if

  900 return

 1100 format(' Itn', i7, ': Reduced Hessian reset')
 1200 format(' Itn', i7, ': Indefinite reduced Hessian')
 1300 format(' Itn', i7, ': Ill-conditioned Z.  Cond(Z) = ', 1p, e8.2)

      end ! subroutine s5getR

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5Hfac
     &   ( iExit, eigH, itn, lenR, m,
     &     maxR, mBS, nb, nS, Hcndbd, Hdmax,
     &     hs, kBS, perm, bl, bu, blBS, buBS, x, xBS, R,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     eigH, iExit, itn, leniw, lenrw, lenR, m, maxR, mBS, nb, nS,
     &     hs(nb), kBS(mBS), perm(maxR), iw(leniw)
      double precision
     &     Hcndbd, Hdmax, blBS(mBS), buBS(mBS), xBS(mBS),
     &     bl(nb), bu(nb), x(nb), R(lenR), rw(lenrw)

!     ==================================================================
!     s5Hfac  factorizes the reduced Hessian Z'HZ.
!
!      iExit       Result
!      -----       ------
!       -2         H  singular   (but should be positive definite)
!       -1         H  indefinite
!        0         positive-definite, factors computed successfully
!
!     13 Oct 1992: First version based on Qpsol routine Qpcrsh.
!     15 Oct 1994: Dependent columns fixed at their current value.
!     02 Aug 2003: snPRNT adopted.
!     01 Sep 2005: Current version of s5Hfac.
!     ==================================================================
      character
     &     str*90
      integer
     &     inform, j, jmax, jS, k, kmax, ksave, nSsave, pivot, rankHz
      double precision
     &     dpiv, eps, Hdmin, s
!     ------------------------------------------------------------------
      integer            INDEF,      SEMDEF,     POSDEF
      parameter         (INDEF = -1, SEMDEF = 0, POSDEF  = 1)
      integer            NoPiv,     Piv
      parameter         (NoPiv = 0, Piv     = 1)
      double precision   one
      parameter         (one   = 1.0d+0)
!     ------------------------------------------------------------------
      iExit = 0
      eps   = max ( Hdmax, one )/Hcndbd
      Hdmin = max ( Hdmax/Hcndbd, eps )

      if (eigH .eq. POSDEF) then
         Pivot = NoPiv
      else if (eigH .eq. INDEF  .or.  eigH .eq. SEMDEF) then
         Pivot =   Piv
      else
         ! Relax -- there are no other options.
      end if

      call s6chol
     &   ( inform, Pivot, maxR, nS, lenR, R, Hdmin, dpiv, rankHz, perm )
!     inform > 0 implies rankHz < nS.

      if (Pivot .eq. Piv) then
!        -----------------------
!        Apply any interchanges.
!        -----------------------
         do j = 1, min(rankHz,nS)
            jmax = perm(j)
            if (jmax .gt. j) then
               kmax       = m + jmax
               k          = m + j

               ksave      = kBS(kmax)
               kBS(kmax)  = kBS(k)
               kBS(k)     = ksave

               s          = xBS(kmax)
               xBS(kmax)  = xBS(k)
               xBS(k)     = s

               s          = blBS(kmax)
               blBS(kmax) = blBS(k)
               blBS(k)    = s

               s          = buBS(kmax)
               buBS(kmax) = buBS(k)
               buBS(k)    = s
            end if
         end do
      end if


      if (dpiv .lt. (-Hdmin)) then
!        ---------------------------------------
!        H  appears to be indefinite.
!        ---------------------------------------
         iExit = -1             ! Indefinite H
      else if (dpiv .lt.   Hdmin) then
!        ---------------------------------------
!        H  appears to be positive semidefinite.
!        rankHz < nS
!        ---------------------------------------
         if (eigH .eq. POSDEF) then ! H should be PD
            iExit = -2
            write(str, 9060) itn, dpiv, Hdmin
            call snPRNT( 21, str, iw, leniw )
         else                       ! singular H allowed
            write(str, 9000) itn, nS-rankHz, Hdmin
            call snPRNT( 21, str, iw, leniw )

            nSsave = nS
            do jS = rankHz+1, nSsave
               k  = m + jS
               j  = kBS(k)

!              Make variable  j  nonbasic (it is already feasible).
!              hs(j) = -1 means x(j) is strictly between its bounds.

               if      (x(j) .le. bl(j)) then
                  x(j) =  bl(j)
                  hs(j) =  0
               else if (x(j) .ge. bu(j)) then
                  x(j) =  bu(j)
                  hs(j) =  1
               else
                  hs(j) = -1
               end if
               if (bl(j) .eq. bu(j)) hs(j) = 4

               nS = nS - 1
            end do
            nS = min( nS, rankHz )
         end if ! rankHz < nS
      end if

      return

 9000 format(' Itn', i7, ': Reduced Hessian appears to have ',
     &         i6, ' small eigenvalues.  PD tol = ', 1p, e9.2 )
 9060 format(' Itn', i7, ': Reduced Hessian appears to be indefinite.',
     &         ' dpiv, Hdmin = ', 1p, e9.2, ',', e9.2 )

      end ! subroutine s5Hfac

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5Hz
     &   ( iExit, Hprod, Hprod1, maxR, lenR,
     &     minimz, m, mBS, n, nb, nnH, nS, Hcalls,
     &     ne, nlocA, locA, indA, Acol,
     &     Hdmax, condZ, targtZ, kBS, R, v, w, y,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1
      integer
     &     Hcalls, iExit, lenR, lencu, leniu, lenru, lencw, leniw,
     &     lenrw, m, maxR, mBS, minimz, n, nb, ne, nlocA, nnH, nS,
     &     locA(nlocA), indA(ne), kBS(mBS), iu(leniu), iw(leniw)
      double precision
     &     condZ, Hdmax, targtZ, Acol(ne), R(lenR), y(nb),
     &     v(nb), w(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s5Hz    computes the reduced Hessian and loads it by columns into
!     the upper triangle R.
!
!      iExit       Status
!      -----       ------
!       -1         Ill-conditioned null-space basis
!        0         reduced Hessian computed successfully
!       >0         Fatal error in LU solve
!
!     13 Oct 1992: First version based on QPSOL routine Qpcrsh.
!     15 Oct 1994: Dependent columns fixed at their current value.
!     04 Dec 2000: R converted to row-wise storage.
!     25 Mar 2005: Current version
!     ==================================================================
      external
     &     dnormi
      integer
     &     jq, jS, ldiag, nBS, Status
      double precision
     &     Anorm, diag, dnormi, eps0, sgnObj
!     ------------------------------------------------------------------
      parameter         (Status = 198) ! Status of a call to Hprod
      integer            Transp
      parameter         (Transp = 1)
      integer            WithB,      WithBt
      parameter         (WithB  = 1, WithBt = 2)
      double precision   zero,           one
      parameter         (zero  = 0.0d+0, one = 1.0d+0)
!     ------------------------------------------------------------------
      eps0      = rw(  2) ! eps**(4/5)       IEEE DP  3.00e-13

      if (nS .eq. 0) return

      iExit   = 0
      nBS     = m + nS
      Hdmax   = zero
      sgnObj  = minimz
      condZ   = one

!     ------------------------------------------------------------------
!     Main loop to find a column of Z'HZ.
!     ------------------------------------------------------------------
      do jS = 1, nS
         !--------------------------------------------------------------
         ! Get the nonlinear elements of the column of Z.
         ! Find y such that B y = column jq.
         ! Scatter the nonlinear part of y into w.
         !--------------------------------------------------------------
         jq  = kBS(m+jS)
         call s2unpk
     &      ( jq, m, n, ne, Anorm, nlocA, locA, indA, Acol, w )
         call s2Bsol
     &      ( iExit, WithB, m, w, y, iw, leniw, rw, lenrw )
         if (iExit .gt. 0) return

         condZ = max ( dnormi( m, y, 1 )/Anorm, condZ )

         if (condZ .ge. targtZ) then
            iExit = -1
            return
         end if

         call s2scatr
     &      ( nnH, m, kBS, (-one), y, w )
         if (jq .le. nnH) w(jq) = one

         ! Set v = H w.

         if (nnH .gt. 0) then
            call Hprod
     &         ( Hprod1, Hcalls, nnH,
     &           w, v, iw(Status),
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
            if (minimz .lt. 0) then
               call dscal ( nnH, sgnObj, v, 1 )
            end if
            iw(Status) = 0
         end if

         !--------------------------------------------------------------
         ! Gather w = vBS and compute v = Z' w.
         ! Solve  B' vB = wB  and  form  wS = wS - S' vB.
         !--------------------------------------------------------------
         call s2gathr
     &      ( nnH, nBS, kBS, one, v, w )
         call s2Bsol
     &      ( iExit, WithBt, m, w, v, iw, leniw, rw, lenrw )
         if (iExit .gt. 0) return

         call s2Bprd
     &      ( Transp, eps0, n, nS, kBS(m+1),
     &        ne, nlocA, locA, indA, Acol,
     &        (-one), v, m, one, w(m+1), nS )

         !--------------------------------------------------------------
         ! Store w(1:nS) in the jS-th row of R.
         ! R is NO LONGER SYMMETRIZED.
         !--------------------------------------------------------------
         call s6Rrow
     &      ( jS, maxR, nS, lenR, R, w(m+1), ldiag )
         diag  = R(ldiag)
         Hdmax = max( Hdmax, abs(diag) )
      end do

      end ! subroutine s5Hz

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5QP
     &   ( iExit, Prob, ProbTag, Elastc, subopt,
     &     Hprod, Hprod1, QPlog, gotR, needLU, typeLU, needx,
     &     lenR, m, maxS, mBS, n, nb, nDegen, Hcalls,
     &     ngQP0, ngQP, ngObj0, ngObj, nnH0, nnH, nS,
     &     itQP, itQPmax, itQPTgt, itn, lEmode, lvlInf, PrtLvl,
     &     minimz, iObj, sclObj, ObjAdd, ObjQP, targtH, targtZ,
     &     tolFP, tolQP, tolx, nInf, sInf, wtInf, piNorm, rgNorm,
     &     ne, nlocA, locA, indA, Acol,
     &     hElast, hEstat, hfeas, hs, kBS,
     &     Ascale, bl, bu, blBS, buBS,
     &     gBS, gObj, gQP, Hdx, pBS, pi, R, rc, rg,
     &     nrhs0, nrhs, rhs, lenx0, nx0, x0, x, xBS, xFreez,
     &     iy, iy1, y, y1, y2,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     Elastc, gotR, needLU, needx
      external
     &     Hprod, Hprod1, QPlog
      integer
     &     Hcalls, iExit, itQP, itQPmax, itQPTgt, itn, iObj,
     &     lencu, leniu, lenru, lencw, leniw, lenrw, lenR, lenx0,
     &     lEmode, lvlInf, PrtLvl, m, maxS, mBS, minimz, ngObj0, ngObj,
     &     n, nb, ne, ngQP0, ngQP, nlocA, nInf, nnH0, nnH, nS, nDegen,
     &     nrhs0, nrhs, nx0, Prob, subopt, typeLU,
     &     locA(nlocA), indA(ne), hElast(nb), hEstat(nb), hs(nb),
     &     kBS(mBS), hfeas(mBS), iy(nb), iy1(nb), iu(leniu), iw(leniw)
      double precision
     &     ObjAdd, ObjQP, piNorm, rgNorm, sclObj, sInf, targtH, targtZ,
     &     tolFP, tolQP, tolx, wtInf, Acol(ne), Ascale(nb), bl(nb),
     &     bu(nb), rc(nb), blBS(mBS), buBS(mBS), gBS(mBS), gObj(*),
     &     gQP(ngQP0), Hdx(nnH0), pBS(mBS), pi(m),
     &     rhs(nrhs0), R(lenR), rg(maxS), x0(lenx0), x(nb), xBS(mBS),
     &     xFreez(nb), y(nb), y1(nb), y2(nb), ru(lenru), rw(lenrw)
      character
     &     ProbTag*20, cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s5QP   solves a linear or quadratic program.
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
!      iExit       Result
!      -----       ------
!       >0         Fatal LU error
!        0         QP solution found
!       -1         QP is infeasible
!       -2         QP is unbounded
!       -3         Too many iterations
!       -4         Weak QP minimizer
!       -5         Too many superbasics
!       -6         QP Hessian not positive semidefinite
!       -7         Z'g could not be made sufficiently small
!       -8         Ill-conditioned Z
!
!     30 Sep 1991: First version of s5QP  based on Qpsol routine qpcore.
!     29 Oct 1993: QP objective computed separately.
!     19 May 1995: Bordered Hessian updated.
!     30 Jul 1995: Border updates removed.
!     04 Jan 1996: Positive semi-definite H treated correctly.
!     20 Jul 1996: Slacks changed to be the row value.
!     09 Aug 1996: First Min Sum version.
!     15 Jul 1997: Thread-safe version.
!     02 Feb 1998: Piecewise linear line search added.
!     07 Nov 1998: Explicit Hessian option added.
!     24 Dec 1999: Sub-optimization option added.
!     25 Jul 2001: Exit on nS > maxR activated.
!     30 Jul 2003: Superbasic slacks allowed.
!     02 Aug 2003: snEXIT and snPRNT adopted.
!     24 Dec 2003: pi checked for NaN and Inf entries.
!     16 May 2006: Explicit target itQP added
!     ==================================================================
      character
     &     str*132
      external
     &     dnormi
      logical
     &     bndswp, checkx, chkFea, chkpi, deadpt, feasbl, gotE, gotgQP,
     &     gotH, incres, jstFea, jstPh1, LUok, optiml, statpt,
     &     maxRef, needf, needv, needLM, needpi, newB, newLU, newSB,
     &     newx, PosDef, prt10, prtLog, prtSum, QPdone, Singlr, usegQP,
     &     Unbndd
      integer
     &     condZ, eigH, inform, itnfix, itnlim, jq, jBq, jBr, jSq, jSr,
     &     jqSave, kchk, kDegen, kfac, klog, kObj, ksav, kSumm, kp,
     &     kPrc, kPrPrt, lRs, linesP, linesS, lenL0, lenL, lenU0, lenU,
     &     LUitn, LUmod, LUsiz0, LUmax, LUreq, lvlTol, maxR, MnrHdP,
     &     MnrHdS, mNewSB, nBS, nElast, nFac, nfmove, nFreez, nInfE,
     &     nonOpt, nSmax, nSwap, nUncon, Status, toldj1, toldj2, toldj3,
     &     nfix(2), PrintP, PrintS
      double precision
     &     Anorm, Bgrwth, Bold, c6, condHz, djq0, djq, djqPrt, dnormi,
     &     dRmax, dRmin, dRsq, eps0, eps2, featol, Hdmax, infBnd, normG,
     &     Obj, ObjPrt, ObjSlk, pivot, rgTest, rgTol(2),
     &     rowerr, sgnObj, sInfE, step, tolx0, tolinc, weight
!     ------------------------------------------------------------------
      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one = 1.0d+0)
      integer            loose,      tight
      parameter         (loose  = 1, tight  = 2)
      integer            Intern,     Extern
      parameter         (Intern = 0, Extern = 1)
      integer            mUncon,     Check
      parameter         (mUncon = 1, Check  = 1)
      integer            WithB
      parameter         (WithB  = 1)
      integer            FP,         FPE,        FPS
      parameter         (FP     = 0, FPE    = 3, FPS   = 4)
      integer            BS        , BT
      parameter         (BS     = 2, BT     = 3)
      integer            Init,       Optml,      Cycle
      parameter         (Init   = 0, Optml  = 1, Cycle = 2)
      integer            PDEF
      parameter         (PDEF   = 1)
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
      parameter         (PrintP = 218) ! (on/off) log     status
      parameter         (PrintS = 219) ! (on/off) summary status
      parameter         (linesP = 220) ! # lines in log     file
      parameter         (linesS = 221) ! # lines in summary file
      parameter         (MnrHdP = 223) ! >0 => Minor heading in log file
      parameter         (MnrHdS = 225) ! >0 => Minor heading for iSumm
!     ------------------------------------------------------------------
      eps0      = rw(  2) ! eps**(4/5)       IEEE DP  3.00e-13
      eps2      = rw(  4) ! eps**(1/2)       IEEE DP  1.49e-08
      infBnd    = rw( 70) ! definition of an infinite bound

      maxR      = iw( 52) ! max columns of R.
      kchk      = iw( 58) ! check (row) frequency
      kfac      = iw( 59) ! factorization frequency
      ksav      = iw( 60) ! save basis map
      klog      = iw( 61) ! log/print frequency
      kSumm     = iw( 62) ! Summary print frequency
      kDegen    = iw( 63) ! max. expansions of featol
      itnlim    = iw( 89) ! limit on total iterations
      mNewSB    = iw( 95) ! max # of new superbasics
      eigH      = iw(200) ! -1,0,1 for indef, psd and pdef H
      nFac      = iw(210) ! # of LU factorizations

      iExit     = 0
      c6        = max( 1.0d-6, eps2 )

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

!     ------------------------------------------------------------------
!     s5QP operates in either ``Normal'' or ``Elastic'' mode.
!     Everything is normal unless a weighted sum is being minimized or
!     the constraints are infeasible.
!     The logical feasbl refers to the feasibility of the nonelastics.
!     wtInf  is the optional parameter Infeasibility Weight.
!     ------------------------------------------------------------------
      feasbl = .false.
      gotE   = .false.

      gotH   = nnH  .gt. 0
      gotgQP = ngQP .gt. 0

!     jstPh1 = stop at the end of phase 1 (either regular or elastic)

      jstPh1 = Prob .eq. FP .or.  Prob   .eq. FPE  .or.  Prob  .eq. FPS

!     The phase 2 objective is F1 + wtInf*F2.

      if (Elastc) then
         needf = lvlInf .ne. 2  ! F1 required in phase 2
         needv = lvlInf .ne. 0  ! F2 required in phase 2
      else
         needf = .true.
         needv = .false.
      end if

      ObjSlk = zero
      ObjQP  = zero
      Obj    = zero
      pivot  = zero
      step   = zero
      nInfE  = 0
      jq     = 0
      djq    = zero
      djq0   = zero
      djqPrt = zero
      jBq    = 0                ! x(jBq) is the incoming  BS
      jBr    = 0                ! x(jBr) is the outgoing  BS
      jSq    = 0                ! x(jSq) is the incoming SBS
      jSr    = 0                ! x(jSr) is the outgoing SBS
      jqSave = 0
      kPrPrt = 0
      sgnObj = minimz

      rgTol(loose) = min( tolQP, c6) ! relaxed   rgTol
      rgTol(tight) = eps0            ! stringent rgTol
      lvlTol = tight                 ! working   rgTol

      rw(toldj1) = 100.0d+0
      rw(toldj2) = 100.0d+0

      kPrc   = 0                ! last sec scanned in part. prc
      LUreq  = 0

      bndswp = .false.
      chkFea = .true.
      chkpi  = .true.
      deadpt = .false.
      needpi = .true.
      newLU  = .true.
      newx   = .false.
      Unbndd = .false.
      PosDef = .false.
      QPdone = .false.

!     nUncon  counts the number of unconstrained (i.e., Newton) steps.
!             If the test for a minimizer were scale-independent,
!             Uncon would never be larger than 1.
!     nfmove  counts the number of times that the QP obj is decreased,

      nfmove = 0
      nUncon = 0

!     subopt nonzero implies that optimization occurs with a subset of
!     the variables frozen at their initial values.
!     During suboptimization, nFreez is the number of frozen variables.

      nFreez = 0
      nSmax  = nS + mNewSB
      call s5hs
     &   ( Intern, nb, bl, bu, hs, x )
      call s5dgen
     &   ( inform, Init, PrtLvl, nb, nInf, itn,
     &     featol, tolx, tolinc, hs, bl, bu, x,
     &     itnfix, nfix, tolx0, iw, leniw, rw, lenrw )

!     ======================Start of main loop==========================
!+    do while (.not. QPdone  .and.  iExit .eq. 0)
  100 if       (.not. QPdone  .and.  iExit .eq. 0) then
!        ===============================================================
!        Check the initial  x  and move it onto  ( A  -I )*x = b.
!        If needLU is true, this will require a basis factorization.
!        ===============================================================
!        If necessary,  factorize the basis  ( B = LU ) and set x.
!        If needLU is false on entry to s5QP, the first call to s2Bfac
!        will try to use existing factors.
!        If needLU is true on entry to s5QP, an LU factorization of
!        type typeLU is computed.
!
!        The reason for the requested LU is as follows.
!
!        LUreq =  0  First LU for a given subproblem
!        LUreq =  1  Frequency
!        LUreq =  2  LU nonzeros increased
!        LUreq =  3
!        LUreq =  4
!        LUreq =  5  Singular after LU mod
!        LUreq =  6  Unstable LU mod (growth in new column of U)
!        LUreq =  7  Not enough memory
!        LUreq =  8
!        LUreq =  9
!        LUreq = 10  Row error in setx
!        LUreq = 11  Big  dx   in setx
!
!        LUreq = 20
!        LUreq = 21  Iterative refinement failed in QP
!        LUreq = 22  Unbounded QP
!        LUreq = 23  Infeasibility after refactorization
!        LUreq = 24  Small directional derivative in QP
!        LUreq = 25  Ill-conditioned Z in QP
!        LUreq = 26  Indefinite Z'HZ in QP
!        LUreq = 27  R singular after bound swap in QP
!        ---------------------------------------------------------------
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
               gotR   = .false. ! Reset R.
               if (prt10) iw(MnrHdP) = 1
            end if

            if (iExit .ne. 0) go to 100

            gotE   = .false.    ! Check hEstat in elastic mode.
            needpi = .true.     ! Recalculate the pi's.
            needx  = .false.
            newx   = .true.
            chkFea = .true.
            chkpi  = .true.

            pivot  = zero
            jqSave = 0
            nUncon = 0
         end if

         nBS    = m + nS
         newSB  = .false.
         nInf   = 0
         sInf   = zero
         optiml = .false.

         call dload ( nBS, zero, gBS, 1 )
         normG  = one

         if (Elastc  .and.  .not. gotE) then
!           ------------------------------------------------------------
!           Reset blBS and buBS for any violated elastics.
!           These values are used in s5step.
!           Strictly feasible elastics are returned to normality.
!           ------------------------------------------------------------
            call s5setE
     &         ( nBS, nb, nElast, featol, infBnd,
     &           hElast, hEstat, kBS,
     &           bl, bu, blBS, buBS, xBS )
            gotE  = .true.
         end if

         if ( chkFea ) then

!           In Phase 1 or just after a factorize, check the feasibility
!           of the basic and superbasic non-elastics.
!           jstFea  indicates that we have just become feasible.
!           jstFea is turned off once a step is taken.

            call s5Inf
     &         ( nBS, featol, nInf, sInf, hfeas, blBS, buBS, gBS, xBS )

            if (nInf .gt. 0) then

!              Non-elastics are infeasible.
!              If necessary, switch back to the feasibility phase, after
!              refactorization (possibly with tighter tols).
!              Print something if the basis has just been refactorized.

               if ( feasbl ) then
                  call s2tryLU
     &               ( itn, 23, nS, LUreq, LUok, typeLU,
     &                 iw, leniw, rw, lenrw )
                  if (.not. LUok) iExit = 11
                  feasbl = .false.
                  go to 100
               end if

               gotR = .false.

               if (prt10  .and.  iw(LUitn) .eq. 0) then
                  write(str, 1030) itn, nInf, sInf
                  call snPRNT( 21, str, iw, leniw )
               end if
            end if

!           feasbl = true means that the nonelastics are feasible.
!                    This defines the start of Phase 2.

            if (.not. feasbl) then
               jstFea = nInf .eq. 0
            end if
            feasbl = nInf .eq. 0
            chkFea = nInf .gt. 0
         end if ! if chkFea

         if ( Elastc ) then
!           ------------------------------------------------------------
!           Compute the sum of infeasibilities of the elastic variables.
!           ------------------------------------------------------------
            call s5InfE
     &         ( nb, nBS, hEstat, kBS, nInfE, sInfE, bl, bu, x )
            nInf = nInf + nInfE
            sInf = sInf + sInfE
         end if

         if (feasbl  .and.  jstPh1) then
            ! The non-elastics are feasible.  Exit.
            condHz = zero
            djqPrt = zero
            rgNorm = zero
            piNorm = zero
            call dload ( m, zero, pi, 1 )
            deadpt = .false.
            optiml = .true.

         else

            if (feasbl) then
!              ---------------------------------------------------------
!              Feasible for the nonelastics.
!              (Elastc = false means no elastics.)
!              ---------------------------------------------------------
!              If just feasible, compute the QP objective (and gradient)
!              and R.

               ObjSlk = zero

               if (iObj .ne. 0) then
                  ObjSlk = xBS(iw(kObj))*sclObj
               end if
               Obj = sgnObj*ObjSlk

               if (jstFea  .or.  newx) then
                  if (needf) then
!                    ===================================================
!                    Initialize the QP objective and gradient.
!                    ObjQP is the linear plus quadratic term of the
!                    objective (not scaled by sgnObj).   It is updated
!                    after each QP step.
!                    ===================================================
                     if ( gotgQP ) then
                        call s5QPfg
     &                     ( Hprod, Hprod1,
     &                       ngQP, ngObj0, ngObj, nnH,
     &                       iw(Status), Hcalls, ObjQP,
     &                       gObj, gQP, lenx0, nx0, x0, x, y,
     &                       cu, lencu, iu, leniu, ru, lenru,
     &                       cw, lencw, iw, leniw, rw, lenrw )
                        Obj    = Obj + sgnObj*ObjQP
                        iw(Status) = 0
                     end if

                     if ( gotH  .and. .not. gotR) then
!                       ------------------------------------------------
!                       Load and factor the reduced Hessian.
!                       This happens after every LU factorize.
!                       If the reduced Hessian is not positive definite,
!                       reduce the LU factor tolerances to get a better
!                       conditioned Z.
!                       ------------------------------------------------
                        if (nS .gt. 0) then
                           call s5Hz
     &                        ( inform, Hprod, Hprod1, maxR, lenR,
     &                          minimz, m, mBS, n, nb, nnH, nS, Hcalls,
     &                          ne, nlocA, locA, indA, Acol,
     &                          Hdmax, rw(condZ), targtZ, kBS,R,y,y1,y2,
     &                          cu, lencu, iu, leniu, ru, lenru,
     &                          cw, lencw, iw, leniw, rw, lenrw )
!                          inform = -1, 0, >0
                           if (inform .ne. 0) then
                              if (inform .eq. -1) then
                                 iExit = -8 ! Ill-conditioned Z
                              else
                                 iExit = inform ! Fatal error in LU
                              end if
                              go to 100
                           end if

                           call s5Hfac
     &                        ( inform, PDEF, itn, lenR, m,
     &                          maxR, mBS, nb, nS, targtH, Hdmax,
     &                          hs, kBS, iy,
     &                          bl, bu, blBS, buBS, x, xBS, R,
     &                          iw, leniw, rw, lenrw )
!                          inform = -2, -1, 0
                           if (inform .ne. 0) then
                              iExit = -6 ! Z'HZ not positive definite
                              go to 100
                           end if
                        end if ! nS > 0
                        gotR = .true.
                     end if ! gotH and not gotR
                  end if ! needf

                  nBS    = m + nS
                  PosDef = gotR  .or.  nS .eq. 0

               end if ! jstFea .or. newLU

!              ---------------------------------------------------------
!              Gather the QP gradient in BS order.
!              Assign the nonzero components of gBS.
!              ---------------------------------------------------------
               if ( needf ) then
                  if (gotgQP) then
                     call s2gathr
     &                  ( ngQP, nBS, kBS, sgnObj, gQP, gBS )
                  end if
                  if (iObj .gt. 0) gBS(iw(kObj)) = sgnObj*sclObj
               end if

               if (Elastc  .and.  nInfE .gt. 0  .and.  needv) then
                  call s5grdE
     &               ( nb, nBS, wtInf, hEstat, kBS, gBS )
               end if

               normG = dnormi( nBS, gBS, 1 )

!              ---------------------------------------------------------
!              See if it's time to suboptimize.
!              NOTE: We must not suboptimize if all steps have been
!              degenerate.
!              ---------------------------------------------------------
               if (subopt .ne. 0  .or.  nfmove .eq. 0) then
!                 Relax
               else
                  if (nS .ge. nSmax ) then
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

            if ( needpi ) then
               call dcopy
     &            ( m, gBS, 1, y, 1 )
               call s5setp
     &            ( inform, m, chkpi, pinorm, y, pi,
     &              iw, leniw, rw, lenrw )
               if (inform .ne. 0) then
                  if (inform .gt. 0) then
                     iExit  = inform
                  else          ! pi is infinite or contains a NaN/Inf.
                     write(str, 1040) itn
                     call snPRNT( 21, str, iw, leniw )
                     call s2tryLU
     &                  ( itn, 11, nS, LUreq, LUok, typeLU,
     &                    iw, leniw, rw, lenrw )
                     if (.not. LUok) iExit = 43
                  end if
                  go to 100
               end if
               needpi = .false.
            end if

            rgNorm = zero
            if (nS .gt. 0) then
               call s5rg
     &            ( m, nBS, n, nS, eps0,
     &              ne, nlocA, locA, indA, Acol,
     &              gBS, pi, rg, rgNorm, kBS )
            end if

!           ============================================================
!           Determine if the reduced Hessian is positive definite.
!           ============================================================
            condHz = zero
            if ( gotR ) then
               call s6Rcnd
     &            ( maxR, nS, lenR, R, dRmax, dRmin, condHz )
            end if

            if (.not. PosDef) then
               if (feasbl  .and.  gotR) then
                  if (nS .gt. 0) then
                     lRs  = (nS - 1)*maxR + (3 - nS)*nS/2 ! Magic formula!
                     dRsq = R(lRs)**2
                  end if
                  call s5Rsng
     &               ( eigH, PosDef, Singlr, itn,
     &                 maxR, lenR, nS, dRsq, R, iw, leniw, rw, lenrw )
               else
                  PosDef = nS .eq. 0
               end if
            end if

!           ============================================================
!           Check for a stationary point.  Use a stringent rgTol after
!           a constrained step to help avoid false stationary points.
!           In theory, the reduced gradient is zero and the reduced
!           Hessian is positive definite after a bound swap.
!
!           If x is a minimizer,  reduced costs are calculated.
!           ============================================================
            if ( feasbl ) then
               rw(toldj3) = tolQP
            else
               rw(toldj3) = tolFP
            end if

            rgTest = max( piNorm, normG )

            if (.not. feasbl) then
               statpt = rgNorm .le. rgTol(loose) *rgTest
            else if (nUncon .ge. 1) then
               statpt = rgNorm .le. rgTol(lvlTol)*rgTest
            else
               statpt = rgNorm .le. rgTol(tight) *rgTest
            end if

            if (feasbl) then

               maxRef =  nUncon .gt. mUncon

               if ((maxRef .or. bndswp)  .and.  .not. statpt) then

!                 If this point should be stationary but isn't.
!                 If possible, relax the reduced-gradient tolerance.

                  if (lvlTol .eq. tight) then
                     lvlTol = loose
                     statpt = rgNorm .le. rgTol(lvlTol)*rgTest
                  end if

                  if (.not. statpt) then
                     call s2tryLU
     &                  ( itn, 21, nS, LUreq, LUok, typeLU,
     &                    iw, leniw, rw, lenrw )
                     if (.not. LUok) then
                        iExit = -7 ! Large Z'g
                        go to 100
                     end if
                  end if
               end if
               deadpt = statpt  .and.  needf  .and. .not. PosDef
            end if

            if (statpt  .or.  bndswp) then
               jqSave = 0
               nUncon = 0
               bndswp = .false.

               if (gotR  .and.  .not. posdef) then
                  write(str, 1600) itn
                  call snPRNT( 23, str, iw, leniw )
                  call s2tryLU
     &               ( itn, 27, nS, LUreq, LUok, typeLU,
     &                 iw, leniw, rw, lenrw )
                  if (.not. LUok) iExit = 44 ! Ill-conditioned Z
                  go to 100
               end if
            end if

            needLM = statpt

            kPrPrt = kPrc
            jq     = 0

            if ( needLM ) then
!              ---------------------------------------------------------
!              Compute Lagrange multipliers.
!              ---------------------------------------------------------
               djq0   = djq     ! save djq in case of bad statpt
               djq    = zero
               nUncon = 0
               usegQP = feasbl  .and.  needf  .and.  gotgQP
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
            end if ! needLM
         end if ! jstPh1

         QPdone = optiml  .or.  deadpt  .or.  Unbndd

         if ( QPdone ) then
!           ------------------------------------------------------------
!           Apparently we are finished.
!           See if any nonbasics have to be set back on their bounds.
!           ------------------------------------------------------------
            call s5dgen
     &         ( inform, Optml, PrtLvl, nb, nInf, itn,
     &           featol, tolx, tolinc, hs, bl, bu, x,
     &           itnfix, nfix, tolx0,
     &           iw, leniw, rw, lenrw )

            QPdone = inform .eq. 0

            if ( QPdone ) then
!              ---------------------------------------------------------
!              So far so good.  Now check the row residuals.
!              ---------------------------------------------------------
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
               if (deadpt) iExit = -4
            else
               needx  = .true.
               Unbndd = .false.
               needpi = .true.
               go to 100
            end if
         end if ! done

         if (jstPh1  .and.  optiml) then
!           Relax, we are about to exit without printing anything.
         else
!           ============================================================
!           Print the details of this iteration.
!           ============================================================
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
     &         ( Prob, ProbTag,
     &           Elastc, gotR, jstFea, feasbl,
     &           m, mBS, nnH, nS, jSq, jBr, jSr,
     &           iw(linesP), iw(linesS), itn, itQP, kPrPrt, lvlInf,
     &           pivot, step, nInf, sInf, wtInf,
     &           ObjPrt, condHz, djqPrt, rgNorm, kBS, xBS,
     &           iw, leniw )
         end if
         jBq    = 0
         jBr    = 0
         jSq    = 0
         jSr    = 0
         kPrPrt = 0
         djqPrt = zero

         if ( QPdone ) then
!           ------------------------------------------------------------
!           Convergence.
!           ------------------------------------------------------------
            if (nInf .gt. 0) then
               ! No feasible point.  Stop or continue in elastic mode,
               ! depending on required level of infeasibility.

               if (lEmode .gt. 0) then
                  ! Enter elastic mode

                  if (Elastc) then
                     ! Already in elastic mode

                     if (feasbl) then
                        ! Find the final sumInf for the elastics

                        call s5finE
     &                     ( nBS, nb, nInf, sInf, featol,
     &                       hEstat, kBS, bl, bu, xBS )
                     else
                        ! Infeasible (this should not happen)
                        ! The nonelastics cannot be satisfied
                        ! by relaxing the elastics.  Exit.

                        iExit = -1
                     end if
                  else
                     ! Infeasible constraints in Normal mode.
                     ! Start Elastic Phase 1.

                     if (prt10) then
                        write(str, 8050) itn, ProbTag
                        call snPRNT( 23, str, iw, leniw )
                        write(str, 8060) itn
                        call snPRNT( 23, str, iw, leniw )
                        iw(MnrHdP) = 1
                        iw(MnrHdS) = 1
                     end if

                     Elastc = .true.
                     QPdone = .false.
                     needf  = lvlInf .ne. 2 ! Need F1 in phase 2
                     needv  = lvlInf .ne. 0 ! Need F2 in phase 2
                     needpi = .true.
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
                     call snPRNT( 11, str, iw, leniw )
                  end if
               else
                  if (prt10  .and.  klog .eq. 1) then
                     write(str, 1020)          rgnorm, pinorm
                     call snPRNT( 11, str, iw, leniw )
                  end if
               end if
            end if

         else ! not done
!           ------------------------------------------------------------
!           A nonbasic has been selected to become superbasic.
!           Compute the vector y such that B y = column jq.
!           ------------------------------------------------------------
            if ( newSB ) then
!              ---------------------------------------------------------
!              The price has selected a nonbasic to become superbasic.
!              ---------------------------------------------------------
               if (nS+1 .gt. maxR) then
                  iExit = -5
                  go to 100
               end if

               lvlTol = tight
               djqPrt = djq

!              ---------------------------------------------------------
!              Compute the vector pBS such that B pB = column jq.
!              pBS is a multiple of part of the new column of  Z  and
!              is used to define the QP search direction and update R.
!              ---------------------------------------------------------
!              Unpack column jq into  y1  and solve  B*pB = y1.
!              The altered  y1  satisfies  L*y1 = ajq.
!              It is used later in s5QPit to modify L and U.

               call s2unpk
     &            ( jq, m, n, ne, Anorm, nlocA, locA, indA, Acol, y1 )
               call s2Bsol
     &            ( iExit, WithB, m, y1, pBS, iw, leniw, rw, lenrw )
               if (iExit .ne. 0) return
            end if

!           ============================================================
!           Take a step.
!           ============================================================
            if (itn  .ge. itnlim  .or.  itQP .ge. itQPmax) then
               iExit = -3
               go to 100
            end if

            itQP   = itQP   + 1
            itn    = itn    + 1

!           Decide if we want to print something this iteration.

            prtLog = PrtLvl .ge. 1  .and.  mod( itQP, klog  ) .eq. 0
            prtSum = PrtLvl .ge. 1  .and.  mod( itQP, kSumm ) .eq. 0

            iw(PrintP) = 0
            iw(PrintS) = 0
            if (prtLog) iw(PrintP) = 1
            if (prtSum) iw(PrintS) = 1

!           ------------------------------------------------------------
!           Take a reduced gradient step.
!           The new  x  will either minimize the objective on the
!           working set or lie on the boundary of a new constraint.
!           ------------------------------------------------------------
            call s5QPit
     &         ( inform, Hprod, Hprod1, bndswp, Elastc, feasbl,
     &           gotgQP, gotH, gotR, incres, needf, needv,
     &           needpi, newSB, PosDef, itn, lenR,
     &           m, mBS, maxR, maxS, n, nb, Hcalls, nnH0, nnH, nS,
     &           ngQP0, ngQP, nDegen, LUreq, kp, jBq, jSq, jBr, jSr,
     &           jq, jqSave, nfmove, nUncon, djq0, djq,minimz,Obj,ObjQP,
     &           featol, pivot, step, tolinc, wtInf,
     &           ne, nlocA, locA, indA, Acol,
     &           hElast, hEstat, hfeas, hs, kBS,
     &           bl, bu, blBS, buBS, gBS,
     &           gQP, Hdx, pBS, rg, R, x, xBS, y, y1, y2,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )

!           Check for trouble in s5QPit.
!           inform values are -2, -1, 0 >0

            if (inform .ne. 0) then
               if      (inform .gt.  0) then
                  iExit = inform ! Fatal LU error
               else if (inform .eq. -1) then
                  iExit = -2     ! unbounded
               else if (inform .eq. -2) then
                  iExit = -6     ! Hz not positive definite
               end if
               go to 100
            end if

            if (LUreq  .gt. 0) then
               call s2tryLU
     &            ( itn, LUreq, nS, LUreq, LUok, typeLU,
     &              iw, leniw, rw, lenrw )
c$$$               call s2tryLU2
c$$$     &            ( itn, .false., LUreq, nS, LUreq, LUok, typeLU,
c$$$     &              iw, leniw, rw, lenrw )
               if (.not. LUok) then
                  iExit = 43
                  go to 100
               end if
            end if

            iw(LUitn) = iw(LUitn)  + 1
            newLU     = .false.
            newx      = .false.
            chkpi     = .false.

!           Increment featol every iteration.

            featol = featol + tolinc

!           ============================================================
!           Test for error condition and/or frequency interrupts.
!           ============================================================
!           (1) Save a basis map (frequency controlled).
!           (2) Every kdegen iterations, reset featol and move nonbasic
!               variables onto their bounds if they are very close.
!           (3) Refactorize the basis if it has been modified too many
!               times.
!           (4) Update the LU factors of the basis if requested.
!           (5) Check row error (frequency controlled).

            if (mod(itn,ksav) .eq. 0) then
               call s4ksav
     &            ( minimz, m, n, nb, nS, mBS,
     &              itn, nInf, sInf, ObjQP, kBS, hs,
     &              Ascale, bl, bu, x, xBS, cw, lencw, iw, leniw )
            end if

            if (mod( itn, kdegen ) .eq. 0) then
               call s5dgen
     &            ( inform, Cycle, PrtLvl, nb, nInf, itn,
     &              featol, tolx, tolinc, hs, bl, bu, x,
     &              itnfix, nfix, tolx0,
     &              iw, leniw, rw, lenrw )
               needx  = inform .gt. 0
            end if

            if (LUreq .eq. 0) then

               if (     iw(LUmod) .ge. kfac-1) then
                  LUreq  = 1
               else if (iw(LUmod) .ge. 20  .and.
     &                                iw(lenL)+iw(lenU) .gt. LUmax) then
                  Bgrwth = iw(lenL) + iw(lenU)
                  Bold   = LUsiz0
                  Bgrwth = Bgrwth/Bold
                  if ( prt10 ) then
                     write(str, 1000) Bgrwth
                     call snPRNT( 21, str, iw, leniw )
                  end if
                  LUreq  = 2
               else
                  checkx = mod(iw(LUitn),kchk) .eq. 0
                  if (checkx  .and.  .not. needx) then
                     call s5setx
     &                  ( inform, Check, itn,
     &                    m, n, nb, nBS, rowerr,
     &                    ne, nlocA, locA, indA, Acol,
     &                    kBS, xBS, nrhs0, nrhs, rhs, x, y, y2,
     &                    iw, leniw, rw, lenrw )
                     LUreq  = inform
                  end if
               end if
               if (LUreq .gt. 0) typeLU = BT
            end if
         end if ! not optiml

         go to 100
!+    end while
      end if
!     ======================end of main loop============================
!
      call s5hs  ( Extern, nb, bl, bu, hs, x )

      if (subopt .gt. 0) then
         if (nFreez .gt. 0) then
!           Relax
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
 1040 format(' Itn', i7, ': Infinite pi-vector')
 1500 format(' Itn', i7, ': Expanded reduced Hessian ',
     &                   'is indefinite. Basis refactorized')
 1600 format(' Itn', i7, ': Singularity after a ',
     &                   'bound swap.  Basis refactorized')
 1610 format(' Itn', i7, ': Suboptimize: ', i7, ' new superbasics')
 1620 format(' Itn', i7, ': Suboptimize: ', i7, ' minor iterations')
 8050 format(' Itn', i7, ': Infeasible ', a)
 8060 format(' Itn', i7, ': Elastic Phase 1 -- making ',
     &                   'non-elastic variables feasible')

      end ! subroutine s5QP

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5QPfg
     &   ( Hprod, Hprod1,
     &     ngQP, ngObj0, ngObj, nnH,
     &     Status, Hcalls, fQP,
     &     gObj, gQP, lenx0, nx0, x0, x, dx,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1
      integer
     &     Hcalls, lenx0, lencu, lencw, leniu, leniw, lenru, lenrw,
     &     ngQP, ngObj0, ngObj, nnH, nx0, Status, iu(leniu), iw(leniw)
      double precision
     &     fQP, gObj(ngObj0), gQP(ngQP), x0(lenx0), x(ngQP), dx(ngQP),
     &     ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s5QPfg  computes various quantities associated with the LP/QP.
!
!       1.  fQP =  gObj'*(x-x0)  + half*(x - x0)'*H*(x - x0)
!       2.  gQP =  gradient of fQP
*
!     On entry,
!     ngQP         is max( ngObj, nnH )
!     x(ngQP)      are the nonlinear variables
!     x0(ngQP)     is the base point x0
!     gObj(ngObj)  defines the explicit QP linear term
!
!     On exit,
!     fQP          is the QP quadratic term (1) above
!     gQP(ngQP)    is the gradient of fQP
!     dx(ngQP)     is  x-x0
!
!     02 May 1992: First version of s5QPfg.
!     23 Oct 1993: Hx added as an argument.
!     29 Oct 1993: Modified to compute only the QP objective.
!     07 Oct 1994: gQP added as an argument.
!     09 Dec 2004: Current version.
!     ==================================================================
      external
     &     ddot
      integer
     &     nzero
      double precision
     &     ddot
!     ------------------------------------------------------------------
      double precision   zero,          half,          one
      parameter         (zero = 0.0d+0, half = 0.5d+0, one = 1.0d+0)
!     ------------------------------------------------------------------
      if (ngQP .le. 0) return

      call dcopy
     &   ( ngQP,         x , 1, dx, 1 )
      if (nx0 .gt. 0)
     &call daxpy
     &   ( ngQP, (-one), x0, 1, dx, 1 )

      fQP  = zero

      if (nnH .gt. 0) then
         call Hprod
     &      ( Hprod1, Hcalls, nnH,
     &        dx, gQP, Status,
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )
         Hcalls   = Hcalls + 1
         fQP = half*ddot( nnH, dx, 1, gQP, 1 )
      end if

      nzero = ngQP - nnH
      if (nzero .gt. 0) call dload ( nzero, zero, gQP(nnH+1), 1 )

      if (ngObj .gt. 0) then
         fQP = fQP + ddot( ngObj, gObj, 1, dx, 1 )
         call daxpy ( ngObj, one, gObj, 1, gQP, 1 )
      end if

      end ! subroutine s5QPfg

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5QPit
     &   ( iExit, Hprod, Hprod1, bndswp, Elastc, feasbl,
     &     gotgQP, gotH, gotR, incres, needf, needv,
     &     needpi, newSB, PosDef, itn, lenR,
     &     m, mBS, maxR, maxS, n, nb, Hcalls, nnH0, nnH, nS,
     &     ngQP0, ngQP, nDegen, LUreq, kp, jBq, jSq, jBr, jSr,
     &     jq, jqSave, nfmove, nUncon, djq0, djq, minimz, Obj, ObjQP,
     &     featol, pivot, step, tolinc, wtInf,
     &     ne, nlocA, locA, indA, Acol,
     &     hElast, hEstat, hfeas, hs, kBS,
     &     bl, bu, blBS, buBS, gBS,
     &     gQP, Hdx, pBS, rg, R, x, xBS, y, y1, y2,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1
      logical
     &     bndswp, Elastc, feasbl, gotgQP, gotH, gotR, incres,
     &     needf, needv, needpi, newSB, PosDef
      integer
     &     Hcalls, iExit, itn, jBq, jBr, jq, jqSave, kp, lenR, lencu,
     &     lencw, leniu, leniw, lenru, lenrw, LUreq, m, maxR, maxS,
     &     mBS, minimz, n, nb, nDegen, ne, nfmove, nlocA, nnH0, nnH,
     &     ngQP0, ngQP, nS, nUncon, locA(nlocA), indA(ne),
     &     hElast(nb), hEstat(nb), hs(nb), hfeas(mBS), kBS(mBS),
     &     iu(leniu), iw(leniw)
      double precision
     &     djq0, djq, Obj, ObjQP, featol, pivot, step, tolinc,
     &     wtInf, Acol(ne), bl(nb), bu(nb), blBS(mBS), buBS(mBS),
     &     gBS(mBS), gQP(ngQP0), Hdx(nnH0), pBS(mBS), R(lenR), rg(maxS),
     &     x(nb), xBS(mBS), y(nb), y1(nb), y2(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s5QPit performs a QP step.
!
!     On entry,
!        newSB = true implies that variable jq just went superbasic.
!                In this case:
!                pBS  satisfies B pBS = a(jq).
!                y1   satisfies L  y1 = a(jq).
!
!     On exit,
!        pBS contains the most recent QP search direction.
!
!      iExit       Result
!      -----       ------
!       -2         reduced Hessian is not positive semidefinite
!       -1         unbounded
!        0         normal exit
!       >0         Fatal LU error
!
!     25 Nov 1991: First version of s5QPit.
!     05 Jan 1996: Positive semidefinite R treated correctly.
!     29 Aug 1996: First min sum version added.
!     27 Jul 1997: Thread-safe version.
!     02 Feb 1998: Piecewise linear line search added.
!     23 Mar 2000: gQP  and  H  scaled.
!     16 Oct 2000: Reverted to non-bordered version of s5QPit.
!     04 Dec 2000: R converted to row-wise storage.
!     02 Aug 2003: snPRNT adopted.
!     07 May 2006: s5Zp added to compute Z*p.
!     ==================================================================
      character
     &     str*80
      external
     &     ddot, dnormi
      logical
     &     hitcon, hitlow, move, onbnd, Unbndd, Uncon, Hposdf, Hsingr,
     &     Rcheck, Singlr
      integer
     &     eigH, inform, infpiv, jEs, jEsq, jqStat, jr,
     &     jrStat, jsq, jsr, kBSq, ksq, ldiag, LUmod, mtry, nBS,
     &     nBS1, nS1, ntry, Status
      double precision
     &     Anorm, bigdx, bound, ddot, dnormi, drsq, eps, eps0, exact,
     &     gp, gpQP, infBnd, ObjChg, pBS1, pHp, pHpQP, pNorm, sclPiv,
     &     sgnObj, StepB, stepmx, stepP, tolpiv, tolP0, tolP
!     ------------------------------------------------------------------
      integer            xBStox
      parameter         (xBStox = 1)
      parameter         (mtry   = 6)
      integer            WithL,      WithBt
      parameter         (WithL  = 0, WithBt = 2)
      parameter         (LUmod  = 216) ! number of LU mods

      double precision   zero,            half,          one
      parameter         (zero   = 0.0d+0, half = 0.5d+0, one = 1.0d+0)
      double precision   ten
      parameter         (ten    =10.0d+0)
!     ------------------------------------------------------------------
      eps       = rw(  1) ! machine precision.  IEEE DP  2.22e-16
      eps0      = rw(  2) ! eps**(4/5)          IEEE DP  3.00e-13
      tolpiv    = rw( 60) ! excludes small elements of pBS.
      infBnd    = rw( 70) ! definition of an infinite bound
      bigdx     = rw( 72) ! unbounded step.

      eigH      = iw(200) ! -1,0,1 for indef, psd and pdef QP Hessian

      iExit     = 0
      Status    = 0

      Unbndd    = .false.
      sgnObj    = minimz

      nBS       = m + nS

      if ( newSB ) then
!        ---------------------------------------------------------------
!        New superbasic.
!        PosDef must be true if there is a new superbasic.
!        ---------------------------------------------------------------
         nS1    = nS   + 1
         nBS1   = nBS  + 1

         kBS (nBS1) =    jq
         xBS (nBS1) = x (jq)
         blBS(nBS1) = bl(jq)
         buBS(nBS1) = bu(jq)
         jqStat     = hs(jq)

         PosDef = .false.

         if ( gotR ) then
!           ------------------------------------------------------------
!           Add the new column to R at position nS+1.
!           Check for a singular or indefinite reduced Hessian.
!           ------------------------------------------------------------
            kBS(nBS1) = jq
            call s5Rcol
     &         ( iExit, Hprod, Hprod1,
     &           minimz, jq, nS1, dRsq, ldiag,
     &           maxR, lenR, m, mBS, n, nb, Hcalls, nnH, nS1,
     &           ne, nlocA, locA, indA, Acol,
     &           kBS, R, y, y2, pBS,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
            if (iExit .ne. 0) go to 900
            call s5Rsng
     &         ( eigH, Hposdf, Hsingr, itn,
     &           maxR, lenR, nS1, dRsq, R, iw, leniw, rw, lenrw )

            if ( feasbl ) then
               PosDef = Hposdf
               Singlr = Hsingr
               if (.not. (PosDef  .or.  Singlr)) then
                  iExit = -2
                  go to 900
               end if
            else
               gotR   = Hposdf
            end if

            if ( gotR ) R(ldiag) = sqrt( dRsq )

         end if ! gotR

*-->     R can be checked here.

         Rcheck = .false.

         if (Rcheck  .and.  gotR) then
            call s6Rchk
     &         ( iExit, Hprod, Hprod1, itn, minimz,
     &           maxR, lenR, m, mBS, n, nb, Hcalls, nnH, nS1,
     &           ne , nlocA, locA, indA, Acol,
     &           kBS, R, y, y2, pBS,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
            if (iExit .ne. 0) go to 900
         end if

         if ( incres ) then
            jqSave =   jq
         else
            jqSave = - jq
         end if

         nS     =  nS1
         nBS    = nBS1

         hfeas(nBS) = 0

         if (feasbl .and. needf .and. gotgQP .and. jq .le. ngQP) then
            gBS(nBS) = sgnObj*gQP(jq)
         else
            gBS(nBS) = zero
         end if

!        ===============================================================
!        Set hEstat(jq) and the elastic parts of blBS and buBS.
!        ===============================================================
         if ( Elastc ) then

!           If the new superbasic is an elastic variable
!           and it wants to move infeasible, set its elastic state.

            jEsq = hEstat(jq)

            if (hElast(jq) .gt. 0) then
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

!        ---------------------------------------------------------------
!        In phase 1, or phase 2 for an LP, price can select nonbasics
!        floating free between their bounds with zero reduced cost.
!        We have to check that dqj is not zero.
!        ---------------------------------------------------------------
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

!     ------------------------------------------------------------------
!     Store the free components of the search direction in pBS(1:nBS).
!     First, find the search direction pS for the superbasics, store it
!     in  pBS(m+1:nBS), and find its norm.  Put the search direction for
!     the basic variables in pBS(1)  ,...,pBS(m).
!     ------------------------------------------------------------------
  100 Singlr = .not. PosDef

      call s5getp
     &   ( feasbl, gotR, newSB, PosDef,
     &     maxR, lenR, nS, djq, R, rg, pBS(m+1), gp, pHp )

      pBS1   = pBS(m+nS)

      pNorm  = dnormi( nS, pBS(m+1), 1 )

      call s5Zp
     &   ( iExit, m, mBS, n, nb, nS, eps0, pNorm,
     &     ne, nlocA, locA, indA, Acol,
     &     kBS, pBS, y2, iw, leniw, rw, lenrw )
      if (iExit .ne. 0) go to 900

      if ( feasbl ) then
!        ---------------------------------------------------------------
!        If R is singular, ensure that pBS is a feasible direction.
!        A nonzero exit value of inform implies that the directional
!        derivative is too small to be relied upon.
!        ---------------------------------------------------------------
         if (gotR  .and.  Singlr) then
            call s5chkp
     &         ( inform, itn, nBS, jqSave, kBS, gp, pBS, iw, leniw )
            if (inform .gt. 0) then
               LUreq  = 24
               go to 900
            end if
         end if

         if (newSB  .and.  gotR) then

!           Check for a feasible direction.
!           A large  rgTol  may give a pBS(nBS) with the wrong sign.

            if (djq*pBS1 .gt. zero) then

               write(str, 1000) itn
               call snPRNT( 23, str, iw, leniw )

               nS     = nS  - 1
               nBS    = nBS - 1
               hs(jq) = jqStat
               if (Elastc) hEstat(jq) = jEsq
               jq     =  0
               djq    =  djq0
               jSq    = -jSq
               jqSave =  0
               posdef = .true.
               newSB  = .false.
               go to 100
            end if
            bndswp = .false.
         end if

!        ---------------------------------------------------------------
!        Compute y = pBS(scattered) and Hdx(scattered).
!        The vector Hdx is used to update the objective and gradient of
!        the QP.  Form  gpQP  and  pHpQP  for the quadratic.
!        gp = gpQP - pBS(kObj) + terms from the elastic gradient.
!        ---------------------------------------------------------------
         if (needf  .and.  (gotgQP  .or.  gotH)) then
            call s2scatr
     &         ( ngQP, nBS, kBS, one, pBS, y )

            if ( gotgQP ) then
               gpQP  = ddot ( ngQP, gQP, 1, y, 1 )
            end if

            if ( gotH ) then
               pHpQP = zero

               call Hprod
     &            ( Hprod1, Hcalls, nnH,
     &              y, Hdx, Status,
     &              cu, lencu, iu, leniu, ru, lenru,
     &              cw, lencw, iw, leniw, rw, lenrw )
               pHpQP = pHpQP + ddot( nnH, y, 1, Hdx, 1 )
            end if
         end if
      end if ! feasbl

!     ------------------------------------------------------------------
!     Find the nearest constraint in direction  x + step*pBS (step > 0).
!     Exact  is the step that takes xBS(kp) exactly onto bound.
!     It may be positive or slightly negative. (Not defined if Unbndd.)
!
!     If onbnd  is true, step is a step that reaches a bound exactly.
!     xBS(kp) reaches the value bound.  If we take a constrained step,
!     bound is used to put the new nonbasic variable x(jr) exactly on
!     its bound.
!
!     If Unbndd is true, step = stepmx.
!     ------------------------------------------------------------------
      stepmx = bigdx /pNorm
      sclPiv = one
      tolP0  = tolpiv
      tolP   = tolpiv*pNorm
      ntry   = 0

!+    Repeat
  200    tolP  = tolP /sclPiv
         tolP0 = tolP0/sclPiv
         call s5step
     &      ( nBS, nDegen,
     &        featol, infBnd, stepmx, tolinc, tolP,
     &        hfeas, blBS, buBS, xBS, pBS,
     &        hitlow, move, onbnd, Unbndd,
     &        infpiv, kp, bound, exact, stepB, stepP )

!        Find if the step is constrained or unconstrained.
!        If R has been flagged as singular, we double check by trying
!        to compute the QP minimizer along pBS.  If the minimizer
!        exists,  the singularity tolerance must be too large.

         if ( feasbl ) then
            if ( PosDef ) then
               Uncon = stepP .gt. one
            else
               Uncon = stepP*pHp .gt. (- gp)
            end if
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

      hitcon = .not. Uncon
      needpi = .true.

      if ( hitcon ) then
         nUncon = 0
         step   = stepB
      else
         nUncon = nUncon + 1
         pivot  = zero
         if ( PosDef ) then
            step   = one
         else
            step   = (- gp)/pHp
            PosDef = .true.
         end if
      end if

!     ------------------------------------------------------------------
!     Compute ObjChg, the change in ObjQP (minimized or maximized).
!     Note: Obj = sgnObj*ObjQP
!           pHp = sgnObj*pHpQP
!     ------------------------------------------------------------------
      if (feasbl) then
         ObjChg = step*gp + half*pHp*step**2
         Obj    = Obj     + ObjChg

         if ( needf ) then
            if (gotgQP)
     &         ObjQP = ObjQP + step*gpQP
            if (gotH  ) then
               ObjQP = ObjQP + half*pHpQP*step**2
               if (step .gt. zero) then
                  call daxpy ( nnH, step, Hdx, 1, gQP, 1 )
               end if
            end if
         end if
      end if

      if (feasbl  .and.  move) nfmove = nfmove + 1

!     ------------------------------------------------------------------
!     Update the basic variables xBS.
!     ------------------------------------------------------------------
      call daxpy
     &   ( nBS, step, pBS, 1, xBS, 1 )
      call s5BSx
     &   ( xBStox, nBS, nb, kBS, x, xBS )

      if ( hitcon ) then
!        ===============================================================
!        There is a blocking variable.
!        It could be a fixed variable, whose new state must be 4.
!        ===============================================================
         pivot  = - pBS(kp)
         jr     =   kBS(kp)

         bndswp = jr .eq. abs(jqSave)

c$$$!        10 Mar 2004: Care is needed to prevent the
c$$$!        new nonbasic variable jr from ending up slightly inside
c$$$!        its bound.  EXPAND normally ensures that x(jr) will be
c$$$!        ON or slightly OUTSIDE its bound, but now we realise that
c$$$!        rounding error might make it slightly INSIDE.
c$$$
c$$$         if (onbnd) then
c$$$            x(jr) = bound
c$$$         else if (hitlow) then
c$$$            x(jr) = min( x(jr), bl(jr) )
c$$$         else
c$$$            x(jr) = max( x(jr), bu(jr) )
c$$$         end if

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
!           ============================================================
!           A variable in B hit a bound.
!           Find column kSq = kBSq-m  of S to replace column kp of B.
!           If nS = 1 there is no choice.
!           ============================================================
            if (nS .eq. 1) then
               kBSq  = nBS
               pivot = pivot/pBS1
            else
               call dload ( m, zero, y2, 1 )
               y2(kp) = one
               call s2Bsol
     &            ( iExit, WithBt, m, y2, y, iw, leniw, rw, lenrw )
               if (iExit .ne. 0) return

               call s5chzq
     &            ( m, mBS, n, nb, nS, kBSq, pivot, tolP0,
     &              ne, nlocA, locA, indA, Acol,
     &              kBS, bl, bu, xBS, y, iw, leniw, rw, lenrw )
               if (kBSq .le. 0) then
                  write(str, 9999) itn
                  call snPRNT( 23, str, iw, leniw )
                  kBSq   = nBS
               end if
            end if

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
            hs(jBq)    = 3

            if (nS .gt. 1  .and.  gotR) then

!              Finish computing y(m+1), ..., y(m+nS).

               y(kBSq) = - (one + pivot)
               call dscal
     &            ( nS, (one/pivot), y(m+1), 1 )
               call s6Rswp
     &            ( maxR, nS, lenR, R, y2, y(m+1), kSq, eps0 )
            end if

!           ------------------------------------------------------------
!           Get a new  y1, used to modify L and U.  If the outgoing
!           superbasic just came in, we already have it.
!           ------------------------------------------------------------
            if (jSr .ne. jq) then
               call s2unpk
     &            ( jBq, m, n, ne, Anorm, nlocA, locA, indA, Acol, y1 )
               call s2Bsol
     &            ( iExit, WithL, m, y1, y, iw, leniw, rw, lenrw )
               if (iExit .ne. 0) return
            end if

!           Update the LU factors.

            iw(LUmod)  = iw(LUmod) + 1

            call s2Bmod2
     &         ( inform, kp, m, y1, iw, leniw, rw, lenrw )

            if (inform .eq. -1) LUreq = 5 ! Singular after LU mod
            if (inform .eq.  2) LUreq = 6 ! Unstable LU mod
            if (inform .eq.  7) LUreq = 7 ! Insufficient free memory

         else
!           ============================================================
!           A variable in S hit a bound.
!           ============================================================
            hs(jr) = jrStat
            jSr    = jr
            kBSq   = kp
            kSq    = kBSq - m
         end if

!        Delete the kSq-th superbasic and adjust all arrays in BS order.

         call s5Sdel
     &      ( kSq, m, nS, nBS, kBS, blBS, buBS, gBS, rg, xBS )

         if ( gotR ) then
!           ------------------------------------------------------------
!           Cyclically demote column kSq of R to position nS.
!           ------------------------------------------------------------
            if (kSq .lt. nS) then
               call s6Rdel
     &            ( kSq, maxR, nS, lenR, R, eps )
            end if
         end if ! feasbl and gotH

         nS  = nS  - 1
         nBS = nBS - 1

*-->     R can be checked here.

      end if ! hitcon

  900 return

 1000 format(' Itn', i7, ': Bad direction after adding a superbasic.')
 9999 format(' Itn', i7, ': Chzq failed in s5QPit!!')

      end ! subroutine s5QPit

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5Rcol
     &   ( iExit, Hprod, Hprod1,
     &     minimz, jq, jRadd, dRsq, ldiag,
     &     maxR, lenR, m, mBS, n, nb, Hcalls, nnH, nS,
     &     ne , nlocA, locA, indA, Acol,
     &     kBS, R, v, w, y,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1
      integer
     &     iExit, Hcalls, jq, jRadd, ldiag, lenR, maxR, m, minimz, mBS,
     &     n, nb, ne, nlocA, nnH, nS, lencu, leniu, lenru, lencw, leniw,
     &     lenrw, locA(nlocA), indA(ne), kBS(mBS), iu(leniu), iw(leniw)
      double precision
     &     dRsq, Acol(ne), R(lenR), v(nb), w(nb), y(mBS),
     &     ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s5Rcol  computes column jRadd of the Cholesky factor R such that
!     R'R = Z'HZ.  The update corresponds to the addition of a new
!     column to Z.
!
!     On entry,
!        R     holds the columns of the factor associated with the
!              first jRadd-1 columns of Q.
!
!        y     is the vector such that B y = a(jq).
!
!        nS    is the number of columns in R.
!
!     11 Dec 1991: First version based on Qpsol routine Qpcolr.
!     24 Apr 1994: Columns of Nx no longer in Q.
!     27 Oct 2000: Previous version of s5Rcol.
!     04 Dec 2000: R converted to row-wise storage.
!     09 Dec 2004: Current version of s5Rcol.
!     ==================================================================
      external
     &     ddot
      integer
     &     lencol, nBS, Status
      double precision
     &     eps0, Rnrmsq, sgnObj, wHw, ddot
!     ------------------------------------------------------------------
      integer            Transp
      parameter         (Transp = 1)
      integer            WithRt
      parameter         (WithRt = 1)
      integer            WithBt
      parameter         (WithBt = 2)
      double precision   zero,          one
      parameter         (zero = 0.0d+0, one = 1.0d+0)
!     ------------------------------------------------------------------
      eps0      = rw(  2)

      iExit     = 0
      nBS       = m   + nS
      sgnObj    = minimz
      lencol    = min( jRadd-1, nnH )

!     ------------------------------------------------------------------
!     Get w, the vector of nonlinear components of the new column of Z.
!     ------------------------------------------------------------------
!     The input vector y satisfies B y = column jq.
!     Scatter the nonlinear components of y into w.

      call s2scatr
     &   ( nnH, m, kBS, (-one), y, w )
      if (jq .le. nnH) w(jq) = one

!     ------------------------------------------------------------------
!     Compute  H*w  and  w'*H*w.
!     ------------------------------------------------------------------
      wHw = zero

      if (nnH .gt. 0) then
         Status = 0
         call Hprod
     &      ( Hprod1, Hcalls, nnH,
     &        w, v, Status,
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )
         wHw = wHw + ddot ( nnH, w, 1, v, 1 )

         if (minimz .lt. 0) then
            call dscal ( nnH, sgnObj, v, 1 )
            wHw = sgnObj*wHw
         end if
      end if

      Rnrmsq = zero

      if (jRadd .gt. 1) then
         !-------------------------------------------------------------
         ! Gather the nonlinear elements of v in w (= vBS).
         ! Compute Z'w  (solve  B'vB = wB and form  wS = wS - S'vB).
         !-------------------------------------------------------------
         call s2gathr
     &      ( nnH, nBS, kBS, one, v, w )
         call s2Bsol
     &      ( iExit, WithBt, m, w, v, iw, leniw, rw, lenrw  )
         if (iExit .ne. 0) return

         if (nS .gt. 0) then
            call s2Bprd
     &         ( Transp, eps0, n, nS, kBS(m+1),
     &           ne, nlocA, locA, indA, Acol,
     &           (-one), v, m, one, w(m+1), nS )
         end if

         !-------------------------------------------------------------
         ! Solve  R'v = Z(j)'Hw.  Store v in w(m+1:m+jRadd).
         !-------------------------------------------------------------
         call s6Rsol
     &      ( WithRt, maxR, lencol, lenR, R, w(m+1) )
         Rnrmsq = ddot  ( lencol, w(m+1), 1, w(m+1), 1 )
      end if

      if (jRadd .le. nnH) then
         dRsq = wHw - Rnrmsq   ! New diagonal of R.
      else
         dRsq = zero
      end if
      w(m+jRadd) = dRsq

      ! Insert w(m+1:m+jRadd) as column jRadd of R.

      call s6Rcol
     &   ( jRadd, maxR, jRadd, lenR, R, w(m+1), ldiag )

      end ! subroutine s5Rcol

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5rg
     &   ( m, nBS, n, nS, tolz,
     &     ne, nlocA, locA, indA, Acol,
     &     gBS, pi, rg, rgNorm, kBS )

      implicit
     &     none
      integer
     &     m, nBS, n, ne, nlocA, nS, locA(nlocA), indA(ne), kBS(nBS)
      double precision
     &     tolz, rgNorm, Acol(ne), gBS(nBS), pi(m), rg(nS)

!     ==================================================================
!     s5rg    calculates the reduced gradient  rg = gS - S'*pi.
!
!     23 Nov 1991: First version based on Minos routine m7rg.
!     16 Nov 2001: Current version.
!     ==================================================================
      external
     &     dnormi
      double precision
     &     dnormi
!     ------------------------------------------------------------------
      integer            Transp
      parameter         (Transp = 1)
      double precision   one
      parameter         (one = 1.0d+0)
!     ------------------------------------------------------------------
      call dcopy
     &   ( nS, gBS(m+1), 1, rg, 1 )

      call s2Bprd
     &   ( Transp, tolz, n, nS, kBS(m+1),
     &     ne, nlocA, locA, indA, Acol,
     &     (-one), pi, m, one, rg, nS )

      rgNorm = dnormi( nS, rg, 1 )

      end ! subroutine s5rg

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5Rsng
     &   ( eigH, PosDef, Singlr, itn,
     &     maxR, lenR, nS, dRsq, R, iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     PosDef, Singlr
      integer
     &     eigH, itn, maxR, lenR, nS, leniw, lenrw, iw(leniw)
      double precision
     &     dRsq, R(lenR), rw(lenrw)

!     ==================================================================
!     s5Rsng  estimates the inertia of the current reduced Hessian.
!
!     15 Jul 1995: First version of s5Rsng.
!     02 Aug 2003: snPRNT adopted.
!     17 Jun 2004: Current version of s5Rsng.
!     ==================================================================
      character
     &     str*110
      double precision
     &     condH, Hcndbd, dRsqmn, dRmax, dRmin
!     ------------------------------------------------------------------
      integer            PSDEF,     PDEF
      parameter         (PSDEF = 0, PDEF  = 1)
      double precision   zero
      parameter         (zero  = 0.0d+0)
!     ------------------------------------------------------------------
      Hcndbd    = rw( 85) ! bound on the condition of Hz

      if (nS .eq. 0) then
!        ---------------------------------------------------------------
!        Vertices are positive definite by definition.
!        ---------------------------------------------------------------
         PosDef = .true.
         Singlr = .false.

      else
!        ---------------------------------------------------------------
!        Compute dRsqmn, the square of the smallest possible diagonal
!        of a positive-definite reduced Hessian.
!        ---------------------------------------------------------------
         call s6Rcnd
     &      ( maxR, nS-1, lenR, R, dRmax, dRmin, condH )
         dRsqmn = dRmax*(dRmax/Hcndbd)

         PosDef =      dRsq  .ge. dRsqmn
         Singlr =  abs(dRsq) .lt. dRsqmn

         if (Singlr  .or.  PosDef) then
            if (dRsq .lt. zero) then
               if (eigH .eq. PDEF) then
                  write(str, 1000) itn, dRsq, dRsqmn
                  call snPRNT( 21, str, iw, leniw )
               end if
               dRsq = max( zero, dRsq )
            end if
         else if (eigH .eq. PSDEF  .or.  eigH .eq. PDEF) then
            write(str, 9000) itn, dRsq, dRsqmn
            call snPRNT( 21, str, iw, leniw )
         end if
      end if

      return

 1000 format(' Itn', i7, ': Reduced Hessian is semidefinite.',
     &                   ' Square of diag, min diag = ', 1p, 2e9.1 )
 9000 format(' Itn', i7, ': Reduced Hessian is indefinite.',
     &                   ' Square of diag, min diag = ', 1p, 2e9.1 )

      end ! subroutine s5Rsng

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5Sdel
     &   ( kSq, m, nS, nBS, kBS, blBS, buBS, gBS, rg, xBS )

      implicit
     &     none
      integer
     &     kSq, m, nS, nBS, kBS(nBS)
      double precision
     &     blBS(nBS), buBS(nBS), gBS(nBS), rg(nS), xBS(nBS)

!     ==================================================================
!     s5Sdel  deletes the kSqth superbasic variable from the arrays
!     kBS, blBS, blBS, gBS, rg and xBS.
!
!     16 Jun 2001: First version of s5Bswp.
!     16 Jun 2001: Current version.
!     ==================================================================
      integer
     &     j, k
!     ------------------------------------------------------------------
!     Shift all the arrays one place to the left.

      do j  = kSq, nS-1
         k  = m + j
         kBS (k) = kBS (k+1)
         blBS(k) = blBS(k+1)
         buBS(k) = buBS(k+1)
         gBS (k) = gBS (k+1)
         xBS (k) = xBS (k+1)
         rg (j)  = rg  (j+1)
      end do

      end ! subroutine s5Sdel

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5Zp
     &   ( iExit, m, mBS, n, nb, nS, eps0, pNorm,
     &     ne, nlocA, locA, indA, Acol,
     &     kBS, pBS, y, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iExit, leniw, lenrw, m, mBS, n, nb, ne, nlocA, nS,
     &     locA(nlocA), indA(ne), kBS(mBS), iw(leniw)
      double precision
     &     eps0, pNorm, Acol(ne), pBS(mBS), y(nb), rw(lenrw)

!     ==================================================================
!     s5Zp computes the free components of the search direction
!     p = Z pS, where pS is the search direction for the superbasics,
!     stored in  pBS(m+1:nBS)
!
!     On exit, the  free components of the search direction are stored
!     in pBS(1:nBS). The search direction for the basic variables is
!     stored in pBS(1),...,pBS(m).
!
!     20 Dec 2005: First version of s5Zp.
!     20 Dec 2005: This version of s5Zp.
!     ==================================================================
      external
     &     dnormi
      double precision
     &     dnormi
      integer
     &     nBS
!     ------------------------------------------------------------------
      integer            Normal
      parameter         (Normal = 0)
      integer            WithB
      parameter         (WithB  = 1)
      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one = 1.0d+0)
!     ------------------------------------------------------------------
      iExit  = 0
      nBS    = m + nS
      pNorm  = dnormi( nS, pBS(m+1), 1 )

!     First, compute  y = - S*pS and prepare to solve  B*pB = y
!     for pB, the search direction for the basic variables.
!     We first normalize y so the LU solver won't ignore
!     too many "small" elements while computing pB.

      call dscal
     &   ( nS, (one/pNorm), pBS(m+1), 1 )
      call s2Bprd
     &   ( Normal, eps0, n, nS, kBS(m+1),
     &     ne, nlocA, locA, indA, Acol,
     &     (-one), pBS(m+1), nS, zero, y, m )

!     Solve  B*pBS = y  and unnormalize all of pBS.

      call s2Bsol
     &   ( iExit, WithB, m, y, pBS, iw, leniw, rw, lenrw  )
      if (iExit .ne. 0) return

      call dscal
     &   ( nBS, pNorm, pBS, 1 )
      pNorm  = dnormi( nBS, pBS, 1 )

      end ! subroutine s5Zp
