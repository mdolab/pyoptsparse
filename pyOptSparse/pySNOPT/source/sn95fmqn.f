*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  sn95fmqn.f.   Full memory BFGS routines.
*
*     s9FMH0   s9FMup   s9FMHx
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s9FMH0
     &   ( Htype, nnH, Hcndbd, U0pre, Hd, lenU, U, nQNmod )

      implicit
     &     none
      integer
     &     Htype, lenU, nnH, nQNmod
      double precision
     &     Hcndbd, U0pre, Hd(nnH), U(lenU)

*     ==================================================================
*     s9FMH0 resets  U  such that  H = U'U  to the square root of Hd.
*     If U is already diagonal it is set to the identity matrix.
*     On entry, the value of Htype is as follows:
*
*       Htype
*       -----
*       HUnset (-1)      H not set.
*       HNorml ( 0)      H is a Hessian of the form defined by  lvlHes.
*       HDiag  ( 1)      H is a diagonal matrix.
*       HUnit  ( 2)      H is an identity matrix.
*
*     19 Jul 1995: First version of s9FMH0 written by PEG.
*     06 Sep 1998: Pre- and post-QP diagonal Hessian scaling added.
*     13 Jan 2005: Hd always positive semidefinite.
*     14 Jan 2005: Current version.
*     ==================================================================
      logical
     &     overfl
      integer
     &     i, incr, j, l, nzero
      double precision
     &     condHd, ddiv, Hdmax, Hdmin
*     ------------------------------------------------------------------
      integer            HNorml,     HDiag,      HUnit
      parameter         (HNorml = 0, HDiag  = 1, HUnit  = 2)
      double precision   zero
      parameter         (zero   = 0.0d+0)
*     ------------------------------------------------------------------

      if (Htype .eq. HNorml) then
*        ---------------------------------------------------------------
*        Try and set U so that U'U = Hd, where Hd is the diagonal of the
*        Hessian.  First, check the condition of Hd and reset it to
*        U0pre*U0pre*I  if its ill-conditioned or not positive definite.
*        ---------------------------------------------------------------
         overfl = .false.

         Hdmin  = Hd(1)         ! strictly positive in exact arithmetic
         Hdmax  = Hdmin
         do i = 2, nnH
            Hdmin = min( Hd(i), Hdmin)
            Hdmax = max( Hd(i), Hdmax)
         end do

         condHd = ddiv( Hdmax, Hdmin, overfl )

         if (Hdmin .le. zero  .or.  condHd .ge. Hcndbd) then
            Htype  = HUnit
         end if
      end if

      if (Htype .eq. HNorml) then
         Htype = HDiag          ! Set U0 to sqrt(Hd)
      else
         Htype = HUnit          ! Set U0 to U0pre
      end if

*     ------------------------------------------------------------
*     Zero the off-diagonal elements of U.
*     ------------------------------------------------------------
      l     = 1
      incR  = nnH
      nzero = nnH - 1

      do     j = 1, nnH
         if (Htype .eq. HUnit) then
            U (l) = U0pre
            Hd(j) = U0pre*U0pre
         else ! Htype .eq. HDiag
            U(l)  = sqrt(Hd(j))
         end if

         if (j .lt. nnH) then
            call dload ( nzero, zero, U(l+1), 1 )
            l     = l     + incR
            incR  = incR  - 1
            nzero = nzero - 1
         end if
      end do

      nQNmod = 0

      end ! subroutine s9FMH0

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s9FMup
     &   ( Update, Htype, nnH, mQNmod, nQNmod, Hcndbd,
     &     U0pre, U0scal, rydx, rdxHdx, Hd, Hdx, y,
     &     Udx, lenU, U )

      implicit
     &     none
      integer
     &     Update, Htype, lenU, nnH, mQNmod, nQNmod
      double precision
     &     Hcndbd, rydx, rdxHdx, U0pre, U0scal, Hd(nnH), Hdx(nnH),
     &     y(nnH), Udx(nnH), U(lenU)

*     ==================================================================
*     s9FMup applies the full-memory BFGS update to H = U'U.
*     If defined, the self-scaling BFGS update parameter is saved.
*     It is needed to update the reduced Hessian when there are only
*     linear constraints.
*
*     19 Jul 1995: First version of s9FMup.
*     06 Sep 1998: Pre- and post-QP diagonal Hessian scaling added.
*     18 Feb 2001: LM H stored in product form.
*     13 Jan 2005: FM H stored in product form.
*     15 Jan 2005: Current version.
*     ==================================================================
      external
     &     ddot
      integer
     &     i, iExit, lastnz, numU
      double precision
     &     ddot, H0scal, Hdxi, t, told, tolz, Ulast, yi
*     ------------------------------------------------------------------
      integer            HNorml
      parameter         (HNorml = 0)
      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one    = 1.0d+0)
*     ------------------------------------------------------------------
      told   = zero
      tolz   = zero

      if (Update .gt. 1) then
         U0scal = rydx / rdxHdx
         H0scal = U0scal*U0scal ! Self-scaling parameter
         rdxHdx = rdxHdx*U0scal ! Used again for the LC update.

         numU   = nnH*(nnH + 1)/2

         call dscal
     &      ( numU, U0scal, U  , 1 ) ! multiplies U by U0scal.
         call dscal
     &      ( nnH , U0scal, Udx, 1 )
         call dscal
     &      ( nnH , H0scal, Hd , 1 )
         call dscal
     &      ( nnH , H0scal, Hdx, 1 )
      end if

*     Include the latest update in the Hessian diagonal.

      do     i = 1, nnH
         Hdxi  = Hdx(i)/rdxHdx
         yi    =   y(i)/rydx
         Hd(i) = Hd(i) - Hdxi**2 + yi**2
      end do

      if (nQNmod .ge. mQNmod) then
*        ---------------------------------------------------------------
*        Reset H = U'U to be the diagonal of the current H.
*        Discard any updates accumulated so far.
*        ---------------------------------------------------------------
         call s9FMH0
     &      ( Htype, nnH, Hcndbd, U0pre, Hd, lenU, U, nQNmod )
      else
*        ---------------------------------------------------------------
*        Overwrite (Udx,y) with the vectors (Us,v)  such that
*        Us  = Udx / rdxHdx,    v = (1/rydx) gdif - (1/rdxHdx) Hdx.
*
*        Then, U(new) = U + Us v',  with H = U'U.
*
*        Hdx and v  are saved to update R for LC problems.
*        ---------------------------------------------------------------
         nQNmod = nQNmod + 1

         t      = ddot ( nnH, y, 1, Hdx, 1 )
         if (t .ge. zero) then
            call dscal
     &         ( nnH, ( one/  rydx), y, 1 )
         else
            call dscal
     &         ( nnH, (-one/  rydx), y, 1 )
         end if

         call daxpy
     &      ( nnH, (-one/rdxHdx), Hdx, 1, y, 1 )
         call dscal
     &      ( nnH, ( one/rdxHdx), Udx, 1 )

*        ---------------------------------------------------------------
*        Restore  U + Us y' to triangular form  (overwriting Udx).
*        ---------------------------------------------------------------
         Ulast  = zero
         lastnz = nnH

         call s6Rmod
     &      ( iExit, nnH, nnH, lenU, U, Udx, y, lastnz, Ulast,
     &        told, tolz )

         Htype  = HNorml
      end if

      end ! subroutine s9FMup

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s9FMHx
     &   ( nnH, x, Ux, Hx, lenU, U )

      implicit
     &     none
      integer
     &     nnH, lenU
      double precision
     &     U(lenU), Hx(nnH), Ux(nnH), x(nnH)

*     ==================================================================
*     s9FMHx  computes the product Hx = U'Ux, where  U is an upper-
*     triangular matrix stored by rows in the one-dimensional array  U.
*     lenU defines the length of U.  lenU must be at least
*     nnH*(nnH + 1)/2.
*
*     12 Jan 1996: First version of s9FMHx
*     12 Jan 2005: H held as U'U.
*     15 Jan 2005: Current version.
*     ==================================================================
      integer            WithU,      WithUt
      parameter         (WithU  = 0, WithUt = 1)
*     ------------------------------------------------------------------
      call s6Rprd
     &   ( WithU , nnH, nnH, lenU, U,  x, Ux )
      call s6Rprd
     &   ( WithUt, nnH, nnH, lenU, U, Ux, Hx )

      end ! subroutine s9FMHx
