*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  sn90lmqn.f.   Limited-memory BFGS routines.
*
*     s9LMH0   s9LMup   s9LMHx
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s9LMH0
     &   ( Htype, nnH, Hcndbd, U0pre, Hd, U0, nQNmod )

      implicit
     &     none
      integer
     &     Htype, nnH, nQNmod
      double precision
     &     Hcndbd, U0pre, Hd(nnH), U0(nnH)

*     ==================================================================
*     s9LMH0 sets the diagonal U0  such that  H = U0'U0.
*     If the BFGS diagonal Hd is not ill-conditioned, U0 = sqrt(Hd).
*     Otherwise, U0 =  U0pre*I and Hd = U0'U0.
*
*     On entry, the value of Htype is as follows:
*
*       Htype
*       -----
*       HUnset (-1)      H not set.
*       HNorml ( 0)      H is a Hessian of the form defined by  lvlHes.
*       HDiag  ( 1)      H is a diagonal matrix.
*       HUnit  ( 2)      H is a scaled identity matrix.
*
*     19 Jul 1995: First version of s9LMH0.
*     06 Sep 1998: Pre- and post-QP diagonal Hessian scaling added.
*     18 Feb 2001: H stored in product form.
*     13 Jan 2005: Hd always positive semidefinite.
*     15 Jan 2005: Current version.
*     ==================================================================
      external
     &     ddiv
      logical
     &     overfl
      integer
     &     i
      double precision
     &     condHd, ddiv, Hdmax, Hdmin
*     ------------------------------------------------------------------
      integer            HNorml,     HDiag,      HUnit
      parameter         (HNorml = 0, HDiag  = 1, HUnit  = 2)
      double precision   zero
      parameter         (zero   =  0.0d+0)
*     ------------------------------------------------------------------
      if (Htype .eq. HNorml) then
*        ---------------------------------------------------------------
*        Check  Hd, the diagonal of the current BFGS  H.
*        Reset  Hd = U0pre*I  if Hd is not positive definite or is ill
*        conditioned
*        ---------------------------------------------------------------
         overfl = .false.

         Hdmin  = Hd(1)         ! strictly positive in exact arithmetic
         Hdmax  = Hdmin

         do     i = 2, nnH
            Hdmin = min( Hd(i), Hdmin)
            Hdmax = max( Hd(i), Hdmax)
         end do

         condHd = ddiv( Hdmax, Hdmin, overfl )

         if (Hdmin .le. zero  .or.  condHd .ge. Hcndbd) then
            Htype  = HUnit
         end if
      end if

      if (Htype .eq. HNorml) then
*        ---------------------------------------------------------------
*        Set U0 to sqrt(Hd)
*        ---------------------------------------------------------------
         Htype = HDiag
         do     i = 1, nnH
            U0(i) = sqrt(Hd(i))
         end do
      else
*        ---------------------------------------------------------------
*        Set U0 to a multiple of the identity.
*        We come straight here for  HUnset.
*        ---------------------------------------------------------------
         Htype = HUnit
         call dload
     &      ( nnH,        U0pre , U0, 1 )
         call dload
     &      ( nnH, (U0pre*U0pre), Hd, 1 )
      end if

      nQNmod = 0

      end ! subroutine s9LMH0

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s9LMup
     &   ( Update, Htype, nnH, mQNmod, nQNmod, Hcndbd,
     &     U0pre, U0scal, rydx, rdxHdx, Hd, Hdx, y,
     &     dx, U0, S, V )

      implicit
     &     none
      integer
     &     Update, Htype, nnH, mQNmod, nQNmod
      double precision
     &     Hcndbd, rydx, rdxHdx, U0pre, U0scal, dx(nnH),
     &     Hd(nnH), Hdx(nnH), y(nnH), U0(nnH), S(nnH,mQNmod),
     &     V(nnH,mQNmod)

*     ==================================================================
*     s9LMup does almost everything associated with the limited-memory
*     quasi-Newton update.
*
*     If defined, the self-scaling BFGS parameter U0scl is saved.
*     It is needed to update the reduced Hessian when there are only
*     linear constraints.
*
*     19 Jul 1995: First version of s9LMup.
*     06 Sep 1998: Pre- and post-QP diagonal Hessian scaling added.
*     18 Feb 2001: LM H stored in product form.
*     13 Jan 2005: FM H stored in product form.
*     16 Jan 2005: Current version.
*     ==================================================================
      external
     &     ddot
      integer
     &     i
      double precision
     &     ddot, H0scal, Hdxi, t, yi
*     ------------------------------------------------------------------
      integer            HNorml
      parameter         (HNorml = 0)
      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one    = 1.0d+0)
*     ------------------------------------------------------------------
      if (Update .gt. 1) then
         U0scal = rydx  /rdxHdx
         H0scal = U0scal*U0scal
         rdxHdx = rdxHdx * U0scal   ! Used later for the LC update.

         call dscal
     &      ( nnH, U0scal, U0 , 1 ) ! multiplies U by U0scal.
         call dscal
     &      ( nnH, H0scal, Hd , 1 )
         call dscal
     &      ( nnH, H0scal, Hdx, 1 )
      end if

*     Include the latest update in the Hessian diagonal.

      do     i = 1, nnH
         Hdxi  = Hdx(i)/rdxHdx
         yi    =   y(i)/rydx
         Hd(i) = Hd(i) - Hdxi**2 + yi**2
      end do

      if (nQNmod .ge. mQNmod) then
*        ---------------------------------------------------------------
*        Insufficient space for storing the new update pair (s,v).
*        Reset U0 to be the root of the diagonal of the current H.
*        Change nQNmod to discard any updates accumulated so far.
*        ---------------------------------------------------------------
         call s9LMH0
     &      ( Htype, nnH, Hcndbd, U0pre, Hd, U0, nQNmod )
      else
*        ---------------------------------------------------------------
*        Space remains. Store s and v, where
*        U(new) = U(I + sv'), with H = U'U.
*        S  =  dx / rdxHdx,    V =  (1/rydx) gdif - (1/rdxHdx) Hdx.
*
*        Hdx and the modified y (= v) are used to update the reduced
*        Hessian for LC problems.
*        ---------------------------------------------------------------
         nQNmod = nQNmod + 1

         call dcopy
     &      ( nnH,         dx, 1, S(1,nQNmod), 1 )
         call dscal
     &      ( nnH, ( one/rdxHdx), S(1,nQNmod), 1 )

         t = ddot ( nnH, y, 1, Hdx, 1 )
         if (t .ge. zero) then
            call dscal
     &         ( nnH, ( one/rydx), y, 1 )
         else
            call dscal
     &         ( nnH, (-one/rydx), y, 1 )
         end if

         call daxpy
     &      ( nnH, (-one/rdxHdx), Hdx, 1,           y, 1 )
         call dcopy
     &      ( nnH,                  y, 1, V(1,nQNmod), 1 )

         Htype  = HNorml
      end if

      end ! subroutine s9LMup

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s9LMHx
     &   ( nnH, x, Ux, Hx,
     &     mQNmod, nQNmod, U0, S, V )

      implicit
     &     none
      integer
     &     nnH, mQNmod, nQNmod
      double precision
     &     Hx(nnH), U0(nnH), Ux(nnH), x(nnH),
     &     S(nnH,mQNmod), V(nnH,mQNmod)

*     ==================================================================
*     s9LMHx forms the product  Hx  for the limited-memory s8Hx.
*     H = U'U, where U = U0*(I + s1*v1')*(I + s2*v2')...(I + sk*vk').
*     with  k = nQNmod
*
*     19 Jul 1995: First version of s9LMHx
*     18 Feb 2001: H stored in product form.
*     12 Jan 2005: Ux added as argument.
*     16 Jan 2005: Current version.
*     ==================================================================
      external
     &     ddot
      integer
     &     k
      double precision
     &     c, ddot
*     ------------------------------------------------------------------
      call dcopy ( nnH,  x, 1, Ux, 1 )

*     Multiply by U.

      do k = nQNmod, 1, -1
         c = ddot ( nnH, V(1,k), 1, Ux, 1 )
         call daxpy
     &      ( nnH, c, S(1,k), 1, Ux, 1 )
      end do

      call ddscl
     &   ( nnH, U0, 1, Ux, 1 )

*     Multiply by U'.

      call dcopy
     &   ( nnH, Ux, 1, Hx, 1 )
      call ddscl
     &   ( nnH, U0, 1, Hx, 1 )

      do k = 1, nQNmod
         c = ddot ( nnH, S(1,k), 1, Hx, 1 )
         call daxpy
     &      ( nnH, c, V(1,k), 1, Hx, 1 )
      end do

      end ! subroutine s9LMHx
