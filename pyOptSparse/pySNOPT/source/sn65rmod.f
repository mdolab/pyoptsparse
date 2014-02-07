*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     file  sn65rmod.f
*
*     s6Radd   s6RBFS   s6chol   s6Rchk   s6mchl   s6Rcnd   s6Rcol
*     s6Rdel   s6Rfix   s6Rmod   s6Rprd   s6Rrow   s6Rset   s6Rsol
*     s6Rswp   s6Rupd   s6Rqn
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6Radd
     &   ( maxR, lenR, nS, R )

      implicit
     &     none
      integer
     &     maxR, lenR, nS
      double precision
     &     R(lenR)

*     ==================================================================
*     s6Radd  adds a unit column nS to the upper triangular matrix R.
*
*     12 Jun 2001: First version based on MINOS routine m6radd.
*     12 Jun 2001: Current version of s6Radd.
*     ==================================================================
      integer
     &     incR, k, lR
*     ------------------------------------------------------------------
      double precision   zero,           one
      parameter         (zero = 0.0d+0,  one = 1.0d+0 )
*     ------------------------------------------------------------------
      if (nS .le. maxR) then
         lR   = nS
         incR = maxR
         do k = 1, nS - 1
            R(lR) = zero
            incR  = incR - 1
            lR    = lR   + incR
         end do
      else
         lR   = maxR*(maxR + 1)/2  +  (nS - maxR)
      end if

      R(lR) = one

      end ! subroutine s6Radd

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6RBFS
     &   ( iExit, maxR, nS, n, lenR, told, tolz,
     &     R, u, v, delta, w )

      implicit
     &     none
      integer
     &     iExit, lenR, maxR, n, nS
      double precision
     &     delta, told, tolz, R(lenR), u(nS), v(nS), w(nS)

*     ==================================================================
*     s6RBFS applies the quasi-Newton update to the first nS rows and
*     columns of the n x n upper-triangular matrix R such that H = R'R.
*
*     On entry:
*     R     contains the first nS rows and columns of an n x n
*           upper-triangular matrix R such that R'R = Q'HQ.
*
*     u     contains the first  nS  components of  Rdx/norm(Rdx), where
*           R'*(Rdx) = Q'Hdx.  Note that norm(u)=1  if  nS = n.
*           It is overwritten.
*
*     v     contains the first nS components of the BFGS update vector
*           v such that U(new) = U(I + sv'), with H = U'U.
*           s  = dx / rdxHdx, v = (1/rydx) gdif - (1/rdxHdx) Hdx.
*
*     (delta, w)  is a scalar-vector pair such that
*           delta*w = (1/rdxHdx) Hdx
*                   = U' u.
*
*     30 Dec 1991: First version based on NPSOL routine npbfgs.
*     05 May 1999: Last column-wise version.
*     03 Dec 2000: Converted to row-wise storage.
*     18 Feb 2001: H stored in product form.
*     23 Jul 2001: Excess elements of R implemented.
*     17 Nov 2001: Current version of s6RBFS.
*     ==================================================================
      external
     &     dnrm2
      integer
     &     j, j1, l, lastnz, nR
      double precision
     &     dnrm2, ulast, unz
*     ------------------------------------------------------------------
      double precision   zero,          one
      parameter         (zero = 0.0d+0, one = 1.0d+0)
*     ------------------------------------------------------------------
      nR    = min( nS, maxR )

*     Apply the update in the form  R + u v',  where  u  is held
*     in u and ulast.

      ulast = zero

      if (nS .lt. n) then
         ulast = sqrt( max(zero, one-dnrm2( nS, u, 1 )**2) )
      end if

      if (nS .gt. maxR) then
         ulast = sqrt( dnrm2( nS-maxR, u(maxR+1), 1 )**2 + ulast**2 )
      end if

*     ---------------------------
*     Find the last nonzero in u.
*     ---------------------------
      unz    = ulast
      lastnz = nR + 1

*+    while (lastnz .gt. 1  .and.  unz .le. tolz) do
  100 if    (lastnz .gt. 1  .and.  unz .le. tolz) then
         lastnz = lastnz - 1
         unz    = abs(u(lastnz))
         go to 100
*+    end while
      end if

*     -------------------------------------------------------
*     Restore  R + u v'  to triangular form  (overwriting u).
*     -------------------------------------------------------
      call s6Rmod
     &   ( iExit, maxR, nR, lenR, R, u, v, lastnz, ulast, told, tolz )

*     Deal with surplus diagonals of  R.
*     They are defined to be the diagonals of the rank-two
*     modification:  v w' + w v' + v v' (where w is already
*     scaled by delta).

      if (nS .gt. maxR) then
         j1   = maxR + 1
         l    = maxR*j1/2
         do j = j1, nS
            l = l + 1
*           Rdsq = max( told, (R(l)**2 + v(j)**2 + two*delta*w(j)*v(j)))
*           R(l) = sqrt(Rdsq)
            R(l) = one
         end do
      end if

      end ! subroutine s6RBFS

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6chol
     &   ( iExit, Swap, maxH, nH, lenH, H, Hdmin, dmax, iRank, perm )

      implicit
     &     none
      integer
     &     iExit, iRank, lenH, maxH, nH, Swap, perm(nH)
      double precision
     &     H(lenH), Hdmin, dmax

*     ==================================================================
*     s6chol  forms the upper-triangular Cholesky factor R such that
*     H = R'R  or  P H P' = R'R for some permutation P.
*
*     On entry,
*        Swap   specifies the pivoting strategy.
*        Swap   = 0 means let P = I (no interchanges).
*               Otherwise P chooses the maximum diagonal at each stage.
*
*        Hdmin  is the square of the smallest allowable diagonal of R.
*
*     On exit,
*        perm   contains details of the permutation matrix P, such that
*               perm(k) = k  if no interchange occurred at the kth step,
*               perm(k) = j  (k < j <= nH)  if rows k and j were
*                            interchanged at the kth step.
*
*     Only the diagonal and super-diagonal elements of H are used,
*     and they are overwritten by R.
*
*     28 May 1994: First version of s6chol based on routine chlfac.
*     15 Oct 2000: Last version with column-wise storage for R.
*     05 Dec 2000: Converted to row-wise storage.
*     06 Dec 2000: Mild threshold diagonal pivoting implemented
*                  to reduce frequency of symmetric interchanges.
*                  TDPtol < 1.0 disables TDP.
*     25 Mar 2001: Current version of s6chol.
*     ==================================================================
      integer
     &     i, incj, inck, incR, incS, is, j, js, jx, k, kmax,
     &     ldiag, lmax, ls, lcj, lck, lrj, lrk
      double precision
     &     d, dsmall, Hjx, s
*     ------------------------------------------------------------------
      logical            Permut
      double precision   zero         , TDPtol
      parameter         (zero = 0.0d+0, TDPtol = 1.1d+0)
*     ------------------------------------------------------------------
      iExit  = 0
      iRank  = 0
      if (nH .eq. 0) return

      Permut = Swap .ne. 0
      dsmall = Hdmin

*     ------------------------------------------------------------------
*     Form the Cholesky factorization R'R = H or P H P'.
*     Process the first nH rows of H.
*
*     ldiag  is the location of H(j,j).
*     lmax   is the location of H(k,k) if j and k will be interchanged.
*     ------------------------------------------------------------------
      ldiag  = 1
      incR   = maxH

      do j = 1, nH
         dmax  = H(ldiag)
         kmax  = j
         lmax  = ldiag

         if ( Permut ) then

            ! Find the diagonal of the Schur complement with
            ! maximum absolute value.

            ls   = ldiag + incR
            incS = incR  - 1
            do k = j+1, nH
               if (dmax .lt. H(ls)) then
                   dmax   =  H(ls)
                   kmax   =  k
                   lmax   =  ls
               end if
               ls   = ls   + incS
               incS = incS - 1
            end do
         end if

         if (dmax .le. dsmall) go to 800   ! The diagonal is too small.

         if (H(ldiag)*TDPtol .ge. dmax) then
            kmax = j                       ! Don't interchange after all.
            lmax = ldiag
            dmax = H(ldiag)
         end if

         perm(j) = kmax

         if (kmax .ne. j) then

            ! Perform a symmetric interchange.
            ! First swap row H(j,j+1:kmax) with col H(j+1:kmax,kmax).

            lrj  = ldiag + 1
            lck  = ldiag + kmax - j + incR - 1
            inck = maxH  - j    - 1

            do i = j+1, kmax
               s      = H(lrj)
               H(lrj) = H(lck)
               H(lck) = s
               lrj    = lrj  + 1
               lck    = lck  + inck
               inck   = inck - 1
            end do

            ! Now swap col H(1:j,j) with col H(1:j,kmax)

            lcj  = j
            lck  = kmax
            incj = maxH - 1

            do i = 1, j
               s      = H(lcj)
               H(lcj) = H(lck)
               H(lck) = s
               lcj    = lcj  + incj
               lck    = lck  + incj
               incj   = incj - 1
            end do

            ! Finally swap row H(j,kmax:n) with row H(kmax,kmax:n).

            lrj  = ldiag + kmax - j
            lrk  = lmax

            do i = kmax, nH
               s      = H(lrj)
               H(lrj) = H(lrk)
               H(lrk) = s
               lrj    = lrj + 1
               lrk    = lrk + 1
            end do
         end if ! kmax ne j

         ! Set the diagonal of  R.

         d        = sqrt( dmax )
         H(ldiag) = d
         iRank    = iRank + 1

         if (j .lt. nH) then

            ! Set the super-diagonal elements of the jth row of R.

            jx = ldiag + 1

            do k = j+1, nH
               H(jx) = H(jx)/d
               jx    = jx + 1
            end do

            ! Do a rank-one update to the Schur complement.
            ! Form the upper-triangular part of H = H - x x',
            ! where x is the row H(j,j+1:nH).
            ! H(js,:) = H(js,:) - H(j,js)*H(j,:).

            jx   = ldiag + 1
            ls   = ldiag + incR
            incS = incR  - 1

            do js = j+1, nH
               Hjx = H(jx)
               if (Hjx .ne. zero) then
                  i  = jx
                  k  = ls
                  do is   = js, nH
                     H(k) = H(k) - Hjx*H(i)
                     i    = i + 1
                     k    = k + 1
                  end do
               end if
               jx   = jx   + 1
               ls   = ls   + incS
               incS = incS - 1
            end do
         end if ! rank-one update

         ldiag = ldiag + incR
         incR  = incR  - 1
      end do

*     ==================================================================
*     Test if  H  is not positive definite.
*     ==================================================================
  800 if (iRank .lt. nH) iExit = 1

      end ! subroutine s6chol

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6Rchk
     &   ( iExit, Hprod, Hprod1, itn, minimz,
     &     maxR, lenR, m, mBS, n, nb, Hcalls, nnH, nS,
     &     ne, nlocA, locA, indA, Acol,
     &     kBS, R, v, w, y,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1
      integer
     &     Hcalls, iExit, itn, lencu, leniu, lenru, lencw, leniw,
     &     lenrw, lenR, maxR, m, minimz, mBS, n, nb, ne, nlocA, nnH,
     &     nS, locA(nlocA), indA(ne), kBS(mBS), iu(leniu), iw(leniw)
      double precision
     &     Acol(ne), R(lenR), v(nb), w(nb), y(mBS), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     s6Rchk  computes the Cholesky factor R such that
*     R'R = Z'HZ.  The update corresponds to the addition of a new
*     column to Z.
*
*     On entry,
*        R     holds the columns of the factor associated with the
*              first jRadd-1 columns of Q.
*
*        nS    is the number of columns in R.
*
*     14 Mar 2001: First version based on SNOPT routine s5Rcol.
*     14 Mar 2001: Current version of s6Rchk.
*     ==================================================================
      external
     &     ddot
      logical
     &     Posdef, Singlr
      integer
     &     eigH, jq, jS, ldiag, lencol, nBS, Status
      double precision
     &     Anorm, dRsq, eps0, Rnrmsq, sgnObj, wHw, ddot
*     ------------------------------------------------------------------
      integer            Transp
      parameter         (Transp = 1)
      integer            WithRt
      parameter         (WithRt = 1)
      integer            WithB,      WithBt
      parameter         (WithB  = 1, WithBt = 2)
      double precision   zero,          one
      parameter         (zero   = 0.0d+0, one = 1.0d+0)
*     ------------------------------------------------------------------
      eps0   = rw(  2)
      eigH   = iw(200) ! -1,0,1 for indef, psd and pdef QP Hessian

      nBS    = m   + nS
      sgnObj = minimz
      iExit  = 0

*     ------------------------------------------------------------------
*     Main loop to find a column of Z'HZ.
*     ------------------------------------------------------------------
      do jS = 1, nS
         lencol    = min( jS-1, nnH )
*        ---------------------------------------------------------------
*        Get the nonlinear elements of the column of Z.
*        Find y such that B y = column jq.
*        Scatter the nonlinear part of y into w.
*        ---------------------------------------------------------------
         jq  = kBS(m+jS)
         call s2unpk
     &      ( jq, m, n, ne, Anorm, nlocA, locA, indA, Acol, w )
         call s2Bsol
     &      ( iExit, WithB, m, w, y, iw, leniw, rw, lenrw )
         if (iExit .ne. 0) return
         call s2scatr
     &      ( nnH, m, kBS, (-one), y, w )
         if (jq .le. nnH) w(jq) = one

*        ---------------------------------------------------------------
*        Compute  H*w  and  w'*H*w.
*        ---------------------------------------------------------------
         wHw = zero

         if (nnH .gt. 0) then
            Status = 0
            call Hprod
     &         ( Hprod1, Hcalls, nnH,
     &           w, v, Status,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
            wHw = wHw + ddot ( nnH, w, 1, v, 1 )

            if (minimz .lt. 0) then
               call dscal ( nnH, sgnObj, v, 1 )
               wHw = sgnObj*wHw
            end if
         end if

         Rnrmsq = zero

         if (jS .gt. 1) then
*           ------------------------------------------------------------
*           Gather the nonlinear elements of v in w (= vBS).
*           Compute Z'w  (solve  B'vB = wB and form  wS = wS - S'vB).
*           ------------------------------------------------------------
            call s2gathr
     &         ( nnH, nBS, kBS, one, v, w )
            call s2Bsol
     &         ( iExit, WithBt, m, w, v, iw, leniw, rw, lenrw )
            if (iExit .ne. 0) return
            if (nS .gt. 0) then
               call s2Bprd
     &            ( Transp, eps0, n, nS, kBS(m+1),
     &              ne, nlocA, locA, indA, Acol,
     &              (-one), v, m, one, w(m+1), nS )
            end if

*           ------------------------------------------------------------
*           Solve  R'v = Z(j)'Hw.  Store v in w(m+1:m+jS).
*           ------------------------------------------------------------
            call s6Rsol
     &         ( WithRt, maxR, lencol, lenR, R, w(m+1) )
            Rnrmsq = ddot  ( lencol, w(m+1), 1, w(m+1), 1 )
         end if

         if (jS .le. nnH) then
            dRsq = wHw - Rnrmsq ! New diagonal of R.
         else
            dRsq = zero
         end if
         w(m+jS) = dRsq

*        Insert w(m+1:m+jS) as column jS of R.

         call s6Rcol
     &      ( jS, maxR, jS, lenR, R, w(m+1), ldiag )
         call s5Rsng
     &      ( eigH, Posdef, Singlr, itn,
     &        maxR, lenR, nS, dRsq, R, iw, leniw, rw, lenrw )

         if (.not. (PosDef  .or.  Singlr)) then
            iExit = 6
            return
         end if

         R(ldiag) = sqrt( dRsq )
      end do

      end ! subroutine s6Rchk

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6mchl
     &   ( iExit, Swap, maxH, nH, lenH, H,
     &     Hdmin, eps, dmax, iRank, nmodH, perm, e )

      implicit
     &     none
      integer
     &     iExit, Swap, maxH, nH, lenH, nmodH, iRank, perm(nH)
      double precision
     &     H(lenH), Hdmin, eps, dmax, e(nH)

*     ==================================================================
*     s6mchl  forms a modified Cholesky factorization of the matrix H,
*     such that
*                H + E = R'R  or  P (H + E) P' = R'R,
*     where
*                E is diagonal,
*                R is upper triangular,
*                P is a permutation.
*     If H is sufficiently positive definite, E will be zero.
*
*     On entry,
*        Swap   specifies the pivoting strategy.
*        Swap   = 0 means let P = I (no interchanges).
*               Otherwise P chooses the maximum diagonal at each stage.
*
*        Hdmin  is the square of the smallest allowable diagonal of R.
*
*        eps    causes termination if the remaining diagonals are
*               that small.  Intended to leave singularities unaltered.
*
*        nmodH  is initialized (as in old s6mchl -- not sure why???).
*
*     On exit,
*        iExit = 0 if iRank = nH,
*              = 1 if iRank < nH (so the returned factor is singular).
*
*        dmax   returns the largest diagonal AT THE TIME OF TERMINATION.
*               (More explanation needed.)
*
*        iRank  is the number of rows in R.
*
*        nmodH  is the number of modified diagonals.
*
*        perm   contains details of the permutation matrix P, such that
*               perm(j) = j  if no interchange occurred at the kth step,
*               perm(j) = k  (j < k <= nH)  if rows j and k were
*                            interchanged at the j-th step.
*
*        e      contains the modifications to each diagonal:
*               e(1:nH) >= 0.0,  E = diag(e).
*
*     Only the diagonal and super-diagonal elements of H are used,
*     and they are overwritten by R.
*
*     21 Oct 1993: s6mchl: First version based on routine chlfac.
*     15 Oct 1994: s6mchl: Symmetric interchanges added.
*     28 May 1999: s6mchl: Last version with column-wise storage for R.
*     05 Dec 2000: s6chol: "Unmodified Cholesky" routine
*                          converted to row-wise storage.
*     06 Dec 2000: s6chol: Theshold Diagonal Pivoting implemented.
*     20 Dec 2000: s6mchl: Derived from s6chol and previous s6mchl.
*     ==================================================================
      integer
     &     i, incj, inck, incR, incS, is, j, jr, js, jx, k, kmax,
     &     ldiag, lmax, lr, ls, lcj, lck, lrj, lrk
      double precision
     &     betasq, d, dsmall, Hjx, s, supj, supmax
*     ------------------------------------------------------------------
      logical            Permut
      double precision   zero         , TDPtol
      parameter         (zero = 0.0d+0, TDPtol = 1.1d+0)
*     ------------------------------------------------------------------
*                        TDPtol = 1.1 (say) implements mild
*                        threshold diagonal pivoting (TDP)
*                        to reduce frequency of symmetric interchanges.
*                        TDPtol < 1.0 disables TDP.
*     ------------------------------------------------------------------
      iExit  = 0
      iRank  = 0
      if (nH .eq. 0) return

      Permut = Swap .ne. 0
      dsmall = Hdmin

*     ------------------------------------------------------------------
*     Find the maximum diagonal and super-diagonal elements of H.
*     ------------------------------------------------------------------
      dmax   = dsmall
      supmax = zero
      lr     = 1

      do j = 1, nH
         e(j) = zero
         dmax = max( dmax, abs( H(lr) ) )
         lr   = lr + 1
         do k = j+1, nH
            supmax = max( supmax, abs( H(lr) ) )
            lr     = lr + 1
         end do
      end do

      betasq = max( dmax, (supmax/nH) ) ! Bound on off-diagonal elements.

*     ------------------------------------------------------------------
*     Form the Cholesky factorization R'R = H or P H P'.
*     Process the first nH rows of H.
*
*     ldiag  is the location of H(j,j).
*     lmax   is the location of H(k,k) if j and k will be interchanged.
*     ------------------------------------------------------------------
      ldiag  = 1
      incR   = maxH

      do j = 1, nH
         dmax  = abs( H(ldiag) )
         kmax  = j
         lmax  = ldiag

         if ( Permut ) then

            ! Find the diagonal of the Schur complement with
            ! maximum absolute value.

            ls   = ldiag + incR
            incS = incR  - 1
            do k = j+1, nH
               if (dmax .lt. abs( H(ls) )) then
                   dmax   =  abs( H(ls) )
                   kmax   =  k
                   lmax   =  ls
               end if
               ls   = ls   + incS
               incS = incS - 1
            end do
         end if

         if (dmax .le. eps) go to 800   ! The diagonal is too small.

         if (abs( H(ldiag) )*TDPtol .ge. dmax) then
            kmax = j                    ! Don't interchange after all.
            lmax = ldiag
         end if

         dmax    = H(lmax)              ! NOTE: not abs(.) anymore.
         perm(j) = kmax

         !---------------------------------------------------------
         ! Perform a symmetric interchange.
         !---------------------------------------------------------
         if (kmax .ne. j) then

            ! First swap row H(j,j+1:kmax) with col H(j+1:kmax,kmax).

            lrj  = ldiag + 1
            lck  = ldiag + kmax - j + incR - 1
            inck = maxH  - j    - 1

            do i = j+1, kmax
               s      = H(lrj)
               H(lrj) = H(lck)
               H(lck) = s
               lrj    = lrj  + 1
               lck    = lck  + inck
               inck   = inck - 1
            end do

            ! Now swap col H(1:j,j) with col H(1:j,kmax)

            lcj  = j
            lck  = kmax
            incj = maxH - 1

            do i = 1, j
               s      = H(lcj)
               H(lcj) = H(lck)
               H(lck) = s
               lcj    = lcj  + incj
               lck    = lck  + incj
               incj   = incj - 1
            end do

            ! Finally swap row H(j,kmax:n) with row H(kmax,kmax:n).

            lrj  = ldiag + kmax - j
            lrk  = lmax

            do i = kmax, nH
               s      = H(lrj)
               H(lrj) = H(lrk)
               H(lrk) = s
               lrj    = lrj + 1
               lrk    = lrk + 1
            end do
         end if ! kmax ne j

         !---------------------------------------------------------
         ! Find the largest super-diagonal in the jth row.
         !---------------------------------------------------------
         supj = zero
         jr   = ldiag + 1

         do k = j+1, nH
            supj = max( supj, abs( H(jr) ) )
            jr   = jr + 1
         end do

         if (supj .lt. dsmall) supj = zero

         !---------------------------------------------------------
         ! Set the diagonal of  R.
         !---------------------------------------------------------
         d     = max( dsmall, abs( dmax ), supj**2/betasq )
         e(j)  = d - dmax
         if (e(j) .gt. zero) then
            nModH = nModH + 1
         end if

         d        = sqrt( d )
         H(ldiag) = d
         iRank    = iRank + 1

         if (j .lt. nH) then

            ! Set the super-diagonal elements of the jth row of R.

            jr = ldiag + 1

            do k = j+1, nH
               H(jr) = H(jr)/d
               jr    = jr + 1
            end do

            ! Do a rank-one update to the Schur complement.
            ! Form the upper-triangular part of H = H - x x',
            ! where x is the row H(j,j+1:nH).
            ! H(js,:) = H(js,:) - H(j,js)*H(j,:).

            jx   = ldiag + 1
            ls   = ldiag + incR
            incS = incR  - 1

            do js = j+1, nH
               Hjx = H(jx)
               if (Hjx .ne. zero) then
                  i  = jx
                  k  = ls
                  do is   = js, nH
                     H(k) = H(k) - Hjx*H(i)
                     i    = i + 1
                     k    = k + 1
                  end do
               end if
               jx   = jx   + 1
               ls   = ls   + incS
               incS = incS - 1
            end do
         end if ! rank-one update

         ldiag = ldiag + incR
         incR  = incR  - 1
      end do

*     ==================================================================
*     Test if  H  is not positive definite.
*     ==================================================================
  800 if (iRank .lt. nH) iExit = 1

      end ! subroutine s6mchl

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6Rcol
     &   ( jR, maxR, nR, lenR, R, v, ldiag )

      implicit
     &     none
      integer
     &     jR, maxR, nR, lenR, ldiag
      double precision
     &     R(lenR), v(nR)

*     ==================================================================
*     s6Rcol  inserts a column of the upper-triangular matrix R.
*     v(1:jR) becomes column jR.
*     ldiag   returns the location of the new diagonal element.
*
*     03 Dec 2000: s6Radd derived from MINOS routine m6radd.
*     06 Dec 2000: s6Radd changed to s6Rcol because we also need s6Rrow.
*     ==================================================================
      integer
     &     i, incR, l
*     ------------------------------------------------------------------
      l      = jR
      incR   = maxR
      do i     = 1, jR
         R(l)  = v(i)
         incR  = incR - 1
         l     = l    + incR
      end do

      ldiag  = l - incR

      end ! subroutine s6Rcol

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6Rcnd
     &   ( maxR, nS, lenR, R, dRmax, dRmin, condH )

      implicit
     &     none
      integer
     &     maxR, nS, lenR
      double precision
     &     dRmax, dRmin, condH, R(lenR)

*     ==================================================================
*     s6Rcnd  finds the largest and smallest diagonals of the
*     upper-triangular matrix R, and returns the square of their ratio.
*     This is a lower bound on the condition number of R' R.
*
*     03 Dec 2000: Converted to row-wise storage.
*                  Rmax is NO LONGER OUTPUT.
*     21 Feb 2004: Current version of s6Rcnd.
*     ==================================================================
      external
     &     ddiv
      logical
     &     overfl
      integer
     &     incR, j, j1, l, nR
      double precision
     &     d, ddiv
*     ------------------------------------------------------------------
      double precision   zero,           one
      parameter         (zero = 0.0d+0,  one = 1.0d+0 )
*     ------------------------------------------------------------------
      overfl = .false.
      nR     = min( nS, maxR )

      if (nS .eq. 0) then
         condH = zero
         dRmin = one
         dRmax = one
      else
         dRmax = abs( R(1) )
         dRmin = dRmax

         if (nS .eq. 1) then
            condH = ddiv( one, dRmin**2, overfl )
            if (condH .lt. one) condH = one/condH
         else
            l    = 1
            incR = maxR
            do j = 2, nR
               l     = l   + incR
               incR  = incR - 1
               d     = abs( R(l) )
               dRmin = min( dRmin, d )
               dRmax = max( dRmax, d )
            end do

*           Deal with surplus diagonals of  R.

            if (nS .gt. maxR) then
               j1   = maxR + 1
               l    = maxR*j1/2
               do j = j1, nS
                  l     = l + 1
                  d     = abs( R(l) )
                  dRmin = min( dRmin, d )
                  dRmax = max( dRmax, d )
               end do
            end if
            condH  = ddiv( dRmax**2, dRmin**2, overfl )
         end if
      end if

      end ! subroutine s6Rcnd

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6Rdel
     &   ( jq, maxR, nS, lenR, R, tolz )

      implicit
     &     none
      integer
     &     jq, maxR, nS, lenR
      double precision
     &     R(lenR), tolz

*     ==================================================================
*     s6Rdel  deletes the jq-th column from R.
*     The dimension of R decreases from nS to nS - 1.
*
*     03 Dec 2000: s6Rdel derived from MINOS routine m6rdel.
*     17 Apr 1994: m6rdel converted to row-wise storage.
*     30 Jul 1994: Bottom part of R moved north-west AFTER the sweep
*                  of rotations has eliminated the jq-th row.
*                  This usually means double handling of that part of R,
*                  but it's the only way to skip identity rotations.
*                  (Who knows if there are any.)
*     13 Jun 2001: Implemented surplus elements of R.
*     14 Jun 2001: Current version of s6Rdel.
*     ==================================================================
      integer
     &     i, incR, insav, j, jr, js, k, lr, lrsav, ls, nmove, nR
      double precision
     &     a, b, cs, diag, rk, sk, sn
*     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero = 0.0d+0 )
*     ------------------------------------------------------------------

      if (jq .eq. nS) return

*     ------------------------------------------------------------------
*     Delete the jq-th column of R from the top rows.
*     For the first jq-1 rows, elements R(*,jq+1:nR) of each row
*     are shifted 1 place to the left.
*     ------------------------------------------------------------------
      nR     = min( maxR, nS )
      lr     = jq
      incR   = maxR
      nmove  = nR - jq

      do i = 1, jq - 1
         do k    = lr, lr + nmove - 1
            R(k) = R(k+1)
         end do
         incR = incR - 1
         lr   = lr   + incR
      end do

*     ------------------------------------------------------------------
*     Triangularize the remaining rows of R,
*     using a partial forward sweep of rotations.
*
*     x x x x x x     becomes   x x x x x x
*       x x x x x                 x x x x x
*         . - - - -                 . 0 0 0 0
*           x x x x                   + + + +
*             x x x                     + + +
*               x x                       + +
*                 x                         +
*         |                         |
*        jq                        jq
*
*     The . is not touched because it is later overwritten.
*     ls marks the - being eliminated.
*     lr marks the start of the next + + + row.
*     ------------------------------------------------------------------
      lrsav  = lr
      insav  = incR
      ls     = lr

      do j     = jq + 1, nR
         ls    = ls   + 1
         lr    = lr   + incR
         incR  = incR - 1
         b     = R(ls)
         if (abs( b ) .gt. tolz) then
            a     = R(lr)
            diag  = sqrt( a**2 + b**2 )
            R(lr) = diag

            if (j .lt. nR) then
               cs    = a / diag
               sn    = b / diag
               jr    = lr
               js    = ls

               do k     = j+1, nR
                  jr    = jr + 1
                  js    = js + 1
                  rk    = R(jr)
                  sk    = R(js)
                  R(jr) = cs * rk  +  sn * sk
                  R(js) = sn * rk  -  cs * sk
               end do
            end if
         end if
      end do

*     ------------------------------------------------------------------
*     Shift the + + + triangle up and left.
*     lr marks the start of each + + + row being moved.
*     ls marks the start of its final position.
*     ------------------------------------------------------------------
      lr     = lrsav
      incR   = insav
      nmove  = nR - jq

      do j     = jq + 1, nR
         ls    = lr
         lr    = lr + incR
         call dcopy ( nmove, R(lr), 1, R(ls), 1 )
         incR  = incR  - 1
         nmove = nmove - 1
      end do

*     ------------------------------------------------------------------
*     Deal with surplus diagonals of R.
*     ------------------------------------------------------------------
      if (nS .gt. maxR) then
         if (jq .le. maxR) then

*           Clear out the last column of R.

            lr   = maxR
            incR = maxR
            do k = 1, maxR
               R(lr) = zero
               incR  = incR - 1
               lr    = lr   + incR
            end do
         end if

*        Shift surplus diagonals of R to the left.

         k      = max( maxR, jq )
         lr     = maxR*(maxR + 1)/2  +  (k - maxR)
         nmove  = nS - k
         do   k  = lr, lr + nmove - 1
            R(k) = R(k+1)
         end do
      end if

      end ! subroutine s6Rdel

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6Rfix
     &   ( maxR, nS, lenR, R, dRmax, dRmin, condH, Hcndbd )

      implicit
     &     none
      integer
     &     maxR, nS, lenR
      double precision
     &     dRmax, dRmin, condH, Hcndbd, R(lenR)

*     ==================================================================
*     s6Rfix  is called after s6Rcnd if R is ill-conditioned.
*     It increases the magnitude of small diagonals.
*     dRmax, dRmin, Hcndbd are input.
*            dRmin, condH  are output.
*
*     29 Jul 2003: First version of s6Rfix.
*     29 Jul 2003: Current version of s6Rfix.
*     ==================================================================
      integer
     &     incR, j, j1, l, nR
      double precision
     &     d
*     ------------------------------------------------------------------
      double precision   one,          zero
      parameter         (one = 1.0d+0, zero = 0.0d+0)
*     ------------------------------------------------------------------

      if (nS .eq. 0) return

      nR     = min( nS, maxR )
      dRmin  = max( dRmin, dRmax/Hcndbd )

      if (nS .eq. 1) then
         R(1)   = one   ! Reset!
         condH  = one
      else
         l    = 1
         incR = maxR
         do j = 2, nR
            l     = l   + incR
            incR  = incR - 1
            d     = R(l)
            if (abs(d) .lt. dRmin) then
               if (d .ge. zero) then
                  d =  dRmin
               else
                  d = -dRmin
               end if
               R(l) = d
            end if
         end do

*        Deal with surplus diagonals of  R.

         if (nS .gt. maxR) then
            j1   = maxR + 1
            l    = maxR*j1/2
            do j = j1, nS
               l     = l + 1
               d     = R(l)
               if (abs(d) .lt. dRmin) then
                  if (d .ge. zero) then
                     d =  dRmin
                  else
                     d = -dRmin
                  end if
                  R(l) = d
               end if
            end do
         end if
         condH  = Hcndbd
      end if

      end ! subroutine s6Rfix

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6Rmod
     &   ( iExit, maxR, nS, lenR, R, u, v, lastnz, ulast, told, tolz )

      implicit
     &     none
      integer
     &     iExit, lenR, lastnz, maxR, nS
      double precision
     &     told, tolz, ulast, R(lenR), u(nS), v(nS)

*     ==================================================================
*     s6Rmod  modifies the (nS+1) x nS upper-triangular matrix R so that
*     Q*(R + u v') is upper triangular,  where  Q  is orthogonal.
*     The arrays v and u hold v and the first nS elements of u
*     respectively.  The (nS+1)th component of u is held in ulast.
*     The new R overwrites the old.
*
*     Q is the product of two sweeps of plane rotations (not stored).
*     These affect the (lastnz)th row of R, which is temporarily held
*     in the array u.  Thus, u is overwritten.  v is not altered.
*
*     ulast  holds  u(nS+1).   It is overwritten.
*
*     lastnz points to the last nonzero of u.  The value lastnz = nS+1
*            would always be ok, but sometimes it is known to be less
*            than nS+1, so Q reduces to two partial sweeps of rotations.
*
*     told   is a tolerance on the lastv-th diagonal of R.
*     tolz   is a tolerance for negligible elements in  u.
*
*     On exit,
*     iExit = 1  if the diagonal of R is larger than told,
*           = 2  if not (the diagonal is not modified).
*
*     06 Sep 1991: First version based on Minos routine m6rmod.
*     11 Sep 1994: Modified to update a principal submatrix of R.
*     03 Dec 2000: Converted to row-wise storage.
*     17 Jun 2001: Current version of s6Rmod.
*     ==================================================================
      integer
     &     i, incR, j, l, lastR, ldiag, lm1, nmove
      double precision
     &     root, cs, s, sn, t, u2
*     ------------------------------------------------------------------
      double precision   zero
      parameter        ( zero = 0.0d+0 )
*     ------------------------------------------------------------------

      if (lastnz .le. nS) then
         ulast = u(lastnz)
      end if

*     Copy the (lastnz)th row of R into the end of u.

      lm1    = lastnz - 1
      lastR  = lm1*maxR  +  (3 - lastnz)*lastnz/2
      nmove  = nS - lm1
      if (nmove .gt. 0) call dcopy ( nmove, R(lastR), 1, u(lastnz), 1 )

*     ------------------------------------------------------------------
*     Reduce u to a multiple of e(lastnz) using a partial backward sweep
*     of rotations.  This fills in the (lastnz)th row of R (held in u).
*     ------------------------------------------------------------------
      if (lastnz .gt. 1) then
         u2     = ulast**2
         ldiag  = lastR
         incR   = maxR - lm1

         do i = lm1, 1, -1
            incR  = incR  + 1
            ldiag = ldiag - incR
            s     = u(i)
            u(i)  = zero
            if (abs(s) .gt. tolz) then
               u2    = s**2 + u2
               root  = sqrt(u2)
               cs    = ulast/root
               sn    = s    /root
               ulast = root
               l     = ldiag

               do j = i, nS
                  s    = u(j)
                  t    = R(l)
                  u(j) = cs*s + sn*t
                  R(l) = sn*s - cs*t
                  l    = l + 1
               end do
            end if
         end do
      end if

      call daxpy ( nS, ulast, v, 1, u, 1 )       ! Set u = u  +  ulast*v.

*     ------------------------------------------------------------------
*     Eliminate the front of the (lastnz)th row of R (held in u) using a
*     partial forward sweep of rotations.
*     ------------------------------------------------------------------
      if (lastnz .gt. 1) then
         ldiag  = 1
         incR   = maxR

         do i = 1, lm1
            t     = u(i)
            if (abs(t) .gt. tolz) then
               s        = R(ldiag)
               root     = sqrt(s**2 + t**2)
               cs       = s/root
               sn       = t/root
               R(ldiag) = root
               l        = ldiag

               do j = i+1, nS
                  l    = l + 1
                  s    = R(l)
                  t    = u(j)
                  R(l) = cs*s + sn*t
                  u(j) = sn*s - cs*t
               end do
            end if
            ldiag = ldiag + incR
            incR  = incR  - 1
         end do
      end if

*     Insert the new (lastnz)th row of  R.

      if (nmove .gt. 0) then
         call dcopy ( nmove, u(lastnz), 1, R(lastR), 1 )

*        ---------------------------------------------------------------
*        Test for (unlikely) singularity.
*        ---------------------------------------------------------------
         iExit = 1
         if (abs( R(lastR) ) .le. told) then
            iExit   = 2
*           r(lastr) = told
         end if
      end if

      end ! subroutine s6Rmod

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6Rprd
     &   ( Task, maxR, nS, lenR, R, x, y )

      implicit
     &     none
      integer
     &     Task, maxR, lenR, nS
      double precision
     &     R(lenR), x(nS), y(nS)

*     ==================================================================
*     s6Rprd  forms y = Rx or y = R'x, where R is an
*     upper-triangular matrix of dimension nS stored by rows in R(lenR).
*
*     03 Dec 2000: Converted to row-wise storage.
*     12 Jun 2001: Surplus diagonals allowed in R.
*     31 Jul 2003: Allow maxR = 0.
*     20 Jun 2004: Current version of s6Rprd.
*     ==================================================================
      external
     &     ddot
      integer
     &     i, incR, j, l, nR, numR
      double precision
     &     ddot
*     ------------------------------------------------------------------
      integer            WithR,      WithRt
      parameter         (WithR  = 0, WithRt = 1)
      double precision   zero
      parameter         (zero   = 0.0d+0)
*     ------------------------------------------------------------------
      if (maxR .gt. 0) then
         nR = min( maxR, nS )

         if (Task .eq. WithR) then

*           Form y = R x.

            numR  = nR
            l     = 1
            incR  = maxR

            do    j = 1, nR
               y(j) = ddot  ( numR, R(l), 1, x(j), 1 )
               l    = l   + incR
               numR = numR - 1
               incR = incR - 1
            end do

         else if (Task .eq. WithRt) then

*           Form y = Rtranspose x.

            call dload ( nR, zero, y, 1 )

            l    = 1
            incR = maxR
            numR = nR

            do i = 1, nR
               call daxpy ( numR, x(i), R(l), 1, y(i), 1 )
               l    = l    + incR
               incR = incR - 1
               numR = numR - 1
            end do

         end if
      end if

*     Deal with surplus diagonals of R.

      if (nS .gt. maxR) then
         l = maxR*(maxR + 1)/2
         do j = maxR + 1, nS
            l    = l + 1
            y(j) = R(l)*x(j)
         end do
      end if

      end ! subroutine s6Rprd

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6Rrow
     &   ( jR, maxR, nR, lenR, R, v, ldiag )

      implicit
     &     none
      integer
     &     jR, maxR, nR, lenR, ldiag
      double precision
     &     R(lenR), v(nR)

*     ==================================================================
*     s6Rrow     inserts a row of the upper-triangular matrix R.
*     v(jR+1:nR) becomes row jR.
*     ldiag      returns the location of the new diagonal element.
*
*     06 Dec 2000: First version of s6Rrow, called from s5HZ.
*     ==================================================================
      integer
     &     nmove
*     ------------------------------------------------------------------

      ldiag  = (jR - 1)*maxR  +  (3 - jR)*jR/2          ! Magic formula!
      nmove  = nR - jR + 1
      call dcopy ( nmove, v(jR), 1, R(ldiag), 1 )

      end ! subroutine s6Rrow

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6Rset
     &   ( gotR, maxR, nS, lenR, R, w, condR )

      implicit
     &     none
      logical
     &     gotR
      integer
     &     lenR, maxR, nS
      double precision
     &     condR, R(lenR), w(nS)

*     ==================================================================
*     s6Rset  alters R, the upper-triangular factor of the approximate
*     reduced Hessian.
*
*     If  gotR = .false., R does not exist.
*     In this case, R is initialized to the identity matrix.
*
*     Otherwise, R already exists and we attempt to make it better
*     conditioned by scaling its columns by the square roots of its
*     current diagonals.
*
*     12 Jun 2001: First version based on MINOS routine m6rset.
*     21 Feb 2004: Current version of s6Rset.
*     ==================================================================
      integer
     &     incR, j, k, l, ldiag, nR, nzero
      double precision
     &     diag, dmax, dmin
*     ------------------------------------------------------------------
      double precision   zero,           one
      parameter         (zero = 0.0d+0,  one = 1.0d+0)
*     ------------------------------------------------------------------
      nR     = min( maxR, nS )

      if (nR .eq.  0) then

         condR = zero
         gotR  = .true.

      else if ( gotR ) then
*        ---------------------------------------------------------------
*        Scale the columns of R.
*        ---------------------------------------------------------------
*        Find dmin and dmax, and set w = set of scale factors.

         dmax   = abs( R(1) )
         dmin   = dmax
         ldiag  = 1
         incR   = maxR

         do k = 1, nR
            diag  = abs( R(ldiag) )
            dmax  = max( dmax, diag )
            dmin  = min( dmin, diag )
            w(k)  = one / sqrt( diag )
            ldiag = ldiag + incR
            incR  = incR  - 1
         end do

*        Apply column scales to each row of R.

         ldiag = 1
         incR  = maxR

         do k = 1, nR
            l = ldiag
            do j = k, nR
               R(l) =  R(l) * w(j)
               l    =  l + 1
            end do
            ldiag   = ldiag + incR
            incR    = incR  - 1
         end do

         condR = dmax / dmin

      else ! Set R = the identity.

         ldiag = 1
         incR  = maxR
         nzero = nR - 1

         do k = 1, nR - 1
            R(ldiag) = one
            call dload ( nzero, zero, R(ldiag+1), 1 )
            ldiag    = ldiag + incR
            incR     = incR  - 1
            nzero    = nzero - 1
         end do

         R(ldiag) = one

         do k = maxR+1, nS
            ldiag    = ldiag + 1
            R(ldiag) = one
         end do

         condR = one
         gotR  = .true.

      end if

      end ! subroutine s6Rset

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6Rsol
     &   ( Task, maxR, nS, lenR, R, y )

      implicit
     &     none
      integer
     &     Task, maxR, lenR, nS
      double precision
     &     R(lenR), y(nS)

*     ==================================================================
*     s6Rsol  solves Rx = y or R'x = y, where R is an
*     upper-triangular matrix of dimension nS stored by rows in R(lenR).
*     The solution x overwrites y.
*
*     03 Dec 2000: Converted to row-wise storage.
*     12 Jun 2001: Surplus diagonals allowed in R.
*     31 Jul 2003: Allow maxR = 0.
*     31 Jul 2003: Current version of s6Rsol.
*     ==================================================================
      external
     &     ddot
      integer
     &     i, incR, j, l, nR, numR
      double precision
     &     ddot, t
*     ------------------------------------------------------------------
      integer            WithR,      WithRt
      parameter         (WithR  = 0, WithRt = 1)
*     ------------------------------------------------------------------
      if (maxR .gt. 0) then
         nR = min( maxR, nS )

         if (Task .eq. WithR) then

*           Solve R y = y.

            l     = (nR - 1)*maxR + (3 - nR)*nR/2
            y(nR) = y(nR) / R(l)
            incR  = maxR + 1 - nR
            numR  = 0

            do    j = nR-1, 1, -1
               numR = numR + 1
               incR = incR + 1
               l    = l    - incR
               t    = ddot  ( numR, R(l+1), 1, y(j+1), 1 )
               y(j) = (y(j) - t) / R(l)
            end do

         else if (Task .eq. WithRt) then

*           Solve Rtranspose y = y.

            l    = 1
            incR = maxR
            numR = nR - 1

            do i = 1, nR-1
               y(i) = y(i) / R(l)
               call daxpy ( numR, (-y(i)), R(l+1), 1, y(i+1), 1 )
               l    = l   + incR
               incR = incR - 1
               numR = numR - 1
            end do

            y(nR) = y(nR) / R(l)
         end if
      end if

*     Deal with surplus diagonals of R.

      if (nS .gt. maxR) then
         l = maxR*(maxR + 1)/2
         do j = maxR + 1, nS
            l    = l + 1
            y(j) = y(j) / R(l)
         end do
      end if

      end ! subroutine s6Rsol

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6Rswp
     &   ( maxR, nR, lenR, R, v, w, lastnz, tolz )

      implicit
     &     none
      integer
     &     maxR, nR, lenR, lastnz
      double precision
     &     tolz, R(lenR), v(nR), w(nR)

*     ==================================================================
*     s6Rswp  modifies the upper-triangular matrix R to account for a
*     basis exchange in which the (lastnz)th superbasic variable becomes
*     basic.  R is changed to R + vw', which is triangularized by
*     s6Rmod.  v is the (lastnz)th column of R, and w is input.
*
*     01 Dec 1991: First version based on Minos routine m6bswp.
*     24 Jan 1996: v left unscaled.
*     03 Dec 2000: Converted to row-wise storage.
*     17 Jun 2001: Current version of s6Rswp.
*     ==================================================================
      external
     &     dnormi
      integer
     &     i, incR, inform, l
      double precision
     &     dnormi, told, vlast, vnorm
*     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero = 0.0d+0)
*     ------------------------------------------------------------------
      told = zero               ! Singularity checked elsewhere

*     Set  v =  (lastnz)-th column of  R  and find its norm.

      l      = lastnz
      incR   = maxR - 1

      do i = 1, lastnz
         v(i)  = R(l)
         l     = l    + incR
         incR  = incR - 1
      end do

      vnorm  = dnormi( lastnz, v, 1 )
      vlast  = zero                    ! v(nR+1) = 0

      call s6Rmod
     &   ( inform, maxR, nR, lenR, R, v, w, lastnz, vlast,
     &     told, (vnorm*tolz) )

      end ! subroutine s6Rswp

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6Rupd
     &   ( Update, maxR, lenR, m, n, nBS, nnL, nS,
     &     U0scal, rdxHdx, ne, nlocA, locA, indA, Acol,
     &     kBS, v, Hdx, R, w, y, y2,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Update, maxR, lenR, m, n, nBS,
     &     ne, nlocA, nnL, nS, leniw, lenrw,
     &     indA(ne), locA(nlocA), kBS(nBS), iw(leniw)
      double precision
     &     rdxHdx, U0scal, Acol(ne), R(lenR), v(nnL),
     &     Hdx(nnL), w(nBS), y(nBS), y2(nBS), rw(lenrw)

*     ==================================================================
*     s6Rupd computes the effect of the BFGS update to the full Hessian
*     on the Cholesky factor of the reduced Hessian.
*
*
*     On entry:
*     R     contains the first nS rows and columns of an nnL x nnL
*           upper-triangular matrix R such that R'R = Q'HQ, Q = ( Z Y ).
*
*     v     contains the nnL components of the BFGS update
*           U(new) = U(I + sv'), with H = U'U.
*           s  = dx / rdxHdx, v = (1/rydx) gdif - (1/rdxHdx) Hdx.
*
*     Hdx   contains the product H dx.
*
*     w, y, y2 are work vectors.
*
*     04 Aug 1995: First version based on NPSOL routine npupdt.
*     03 Dec 2000: Converted to row-wise storage.
*     18 Feb 2001: H stored in product form.
*     17 Nov 2001: Current version of s6Rupd.
*     ==================================================================
      integer
     &     i, incR, inform, j, l, nR, numR
      double precision
     &     eps0, told, tolz
*     ------------------------------------------------------------------
      integer            WithBt
      parameter         (WithBt = 2)
      integer            Transp,     WithRt
      parameter         (Transp = 1, WithRt = 1)
      double precision   one
      parameter         (one = 1.0d+0)
*     ------------------------------------------------------------------
      eps0 = rw(  2)
      told = 0.0d+0
      tolz = 0.0d+0

      nR   = min( maxR, nS )

*     ------------------------------------------------------------------
*     Find  y (B) =   v(B),  y (S) =   v(S) - S'Binv v(B)   = Z'v.
*     Find  y2(B) = Hdx(B),  y2(S) = Hdx(S) - S'Binv Hdx(B) = Z'Hdx.
*     ------------------------------------------------------------------
      call s2gathr
     &   ( nnL, nBS, kBS, one, v, y )
      call s2Bsol
     &   ( inform, WithBt, m, y, w, iw, leniw, rw, lenrw )

      if (nS .gt. 0) then
         call s2Bprd
     &      ( Transp, eps0, n, nS, kBS(m+1),
     &        ne, nlocA, locA, indA, Acol,
     &        (-one), w, m, one, y(m+1), nS )
      end if

      call s2gathr
     &   ( nnL, nBS, kBS, one, Hdx, y2 )
      call s2Bsol
     &   ( inform, WithBt, m, y2, w, iw, leniw, rw, lenrw )

      if (nS .gt. 0) then
         call s2Bprd
     &      ( Transp, eps0, n, nS, kBS(m+1),
     &        ne, nlocA, locA, indA, Acol,
     &        (-one), w, m, one, y2(m+1), nS )
      end if

      if (Update .gt. 1) then
*        --------------------
*        Self-scaled BFGS.
*        --------------------
         l      = 1
         incR   = maxR
         numR   = nR

         do i = 1, nR
            call dscal ( numR, U0scal, R(l), 1 )
            l    = l    + incR
            incR = incR - 1
            numR = numR - 1
         end do

*        Deal with surplus diagonals of R.

         if (nS .gt. maxR) then
            l = maxR*(maxR + 1)/2
            do j = maxR + 1, nS
               l    = l + 1
               R(l) = U0scal*R(l)
            end do
         end if
      end if

*     ------------------------------------------------------------------
*     Apply the BFGS update to the first nS columns of R.
*     y(S) = Z'v  and  y2(S) = u/norm(Rdx).
*     ------------------------------------------------------------------
      inform = 0
      call dcopy
     &   ( nS, y2(m+1), 1, w(m+1), 1 )
      call s6Rsol
     &   ( WithRt, maxR, nS, lenR, R, y2(m+1) )
      call dscal
     &   ( nS, (one/rdxHdx), y2(m+1), 1 )
      call s6RBFS
     &   ( inform, maxR, nS, nnL, lenR, told, tolz,
     &     R, y2(m+1), y(m+1), (one/rdxHdx), w(m+1) )

      end ! subroutine s6Rupd

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6Rqn
     &   ( iExit, Exactp, maxR, nS,
     &     lenR, R, gtp, g, g2, p, w, step, told, tolz )

      implicit
     &     none
      logical
     &     Exactp
      integer
     &     iExit, lenR, maxR, nS
      double precision
     &     gtp, step, told, tolz, g(nS), g2(nS), p(nS), R(lenR), w(nS)

!     ==================================================================
!     s6Rqn  applies the BFGS update to the upper-triangular matrix  R,
!     which holds the Cholesky factor of the quasi-Newton approximation
!     of the reduced Hessian.
!
!     R  contains a triangle of size  nR = min( nS, maxR ).
!     If  nS .gt. maxR,  R  also contains a diagonal of size nS - maxR.
!
!     g, g2   hold the gradients.                  g is overwritten.
!     p       holds the search direction.         It is overwritten.
!     w       must satisfy  w = -Rp.              It is overwritten.
!     Exactp  = true, means that  w  also satisfies   R' w = g.
!             This holds if p satisfies R'Rp = -g.
!
!     On exit,
!     iExit = 0  if no update was performed,
!           = 1  if the update was successful,
!           = 2  if it was nearly singular.
!
!     12 Jun 2001: First version based on MINOS routine m6bfgs.
!     14 Apr 2004: Added option allowing arbitrary directions.
!     15 Apr 2004: Current version of s6Rqn.
!     ==================================================================
      external
     &     ddot, dnrm2
      integer
     &     j
      double precision
     &     d, ddot, delta1, delta2, dnrm2, gtp2, pj
!     ------------------------------------------------------------------
      integer            WithRt
      parameter         (WithRt = 1)
      double precision   one
      parameter         (one    = 1.0d+0 )
!     ------------------------------------------------------------------
      iExit  = 0
      gtp2   = ddot  ( nS, g2, 1, p, 1 )
      if (gtp2 .le. 0.91d+0*gtp) return

      delta2 = one / sqrt( step*(gtp2 - gtp) ) ! 1/sqrt(y's)

      if ( Exactp ) then
!        ---------------------------------------------------------------
!        The vector w satisfies R'w = g.
!        ---------------------------------------------------------------
         delta1 = one / sqrt( abs( gtp ) ) ! 1/sqrt(p'Bp)

!        Normalize  w  and change its sign. Now w = Rp/norm(Rp)

         call dscal ( nS, (- delta1), w, 1 )

!        Avoid cancellation error in forming the new vector  p.

         if ( abs( delta1/delta2 - one ) .ge. 0.5d+0) then
            do j = 1, nS
               p(j) = delta2*( g2(j) - g(j) )  +  delta1*g(j)
            end do
         else
            d = delta1 - delta2
            do j = 1, nS
               p(j) = delta2*g2(j)  +  d*g(j)
            end do
         end if

      else
!        ---------------------------------------------------------------
!        The product  R'w = -R'Rp = -Bp is stored in g.
!        ---------------------------------------------------------------
         delta1 = one / dnrm2 ( nS, w, 1 ) ! 1/sqrt(p'Bp)

!        Form   R'*w = -R'Rp in p

         call s6Rprd
     &      (  WithRt, maxR, nS, lenR, R, w, p )

!        Normalize  w  and change its sign. Now w = Rp/norm(Rp)

         call dscal ( nS, (- delta1), w, 1 )

         do j = 1, nS
            pj   = p(j)
            p(j) = delta2*( g2(j) - g(j) )  +  delta1*pj
            g(j) = pj
         end do
      end if

!     Triangularize   R  +  w p'.
!     Deal with surplus diagonals of  R.

      call s6RBFS
     &   ( iExit, maxR, nS, nS, lenR, told, tolz,
     &     R, w, p, (-delta1), g )

      end ! subroutine s6Rqn
