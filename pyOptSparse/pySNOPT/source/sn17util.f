*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  sn17util.f
*
*     Various utility routines for SNOPT.
*     The Level 1 and Level 2  BLAS routines are in sn15blas.f
*
*     ddiv     dddiv    ddscl    dnormi   dnormi   dnrm1s   dload
*     chcopy   icopy    iload    jdamax
*
*     These could be tuned to the machine being used.
*     dload  is used the most.
*
*     ddrand
*
*     Bibs and bobs
*
*     s1init   s1time   s1timp
*
*     s1perm   s1trim
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      double precision   function ddiv( a, b, fail )

      implicit
     &     none
      logical
     &     fail
      double precision
     &     a, b

*     ==================================================================
*     ddiv  returns the value div given by
*
*     div = ( a/b                 if a/b does not overflow,
*           (
*           ( 0.0                 if a .eq. 0.0,
*           (
*           ( sign( a/b )*flmax   if a .ne. 0.0  and a/b would overflow,
*
*     where  flmax  is a large value, via the function name.
*     In addition, if  a/b  would overflow then  fail  is returned as
*     true, otherwise  fail  is returned as false.
*
*     Note that when  a and b  are both zero, fail is returned as true,
*     but  div  is returned as  0.0. In all other cases of overflow
*     div is such that  abs( div ) = flmax.
*
*     When  b = 0  then  sign( a/b )  is taken as  sign( a ).
*
*     15 Nov 1991: First version based on Nag routine f06.
*     01 Jun 1999: Current version.
*     ==================================================================
      external
     &     s1flmn, s1flmx
      double precision
     &     absb, div, flmax, flmin
      double precision
     &     s1flmn, s1flmx
*     ------------------------------------------------------------------
      double precision      one,          zero
      parameter           ( one = 1.0d+0, zero = 0.0d+0 )
      intrinsic             abs, sign
*     ------------------------------------------------------------------
      flmax = s1flmx( )
      flmin = s1flmn( )

      if (a .eq. zero) then
         div = zero
         if (b .eq. zero) then
            fail = .true.
         else
            fail = .false.
         end if
      else
         if (b .eq. zero) then
            div  =  sign( flmax, a )
            fail = .true.
         else
            absb = abs( b )
            if (absb .ge. one) then
               fail = .false.
               if (abs( a ) .ge. absb*flmin) then
                  div = a/b
               else
                  div = zero
               end if
            else
               if (abs( a ) .le. absb*flmax) then
                  fail = .false.
                  div  =  a/b
               else
                  fail = .true.
                  div  = flmax
                  if (((a .lt. zero) .and. (b .gt. zero))  .or.
     &                ((a .gt. zero) .and. (b .lt. zero)))
     &               div = -div
               end if
            end if
         end if
      end if

      ddiv = div

      end ! double precision function ddiv

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dddiv ( n, d, incd, x, incx )

      implicit
     &     none
      integer
     &     n, incd, incx
      double precision
     &     d(*), x(*)

*     dddiv  performs the diagonal scaling  x  =  x / d.

      integer
     &     i, id, ix
      external
     &     dscal
      intrinsic
     &     abs
*     ------------------------------------------------------------------
      double precision   one
      parameter        ( one = 1.0d+0 )
*     ------------------------------------------------------------------
      if (n .gt. 0) then
         if (incd .eq. 0  .and.  incx .ne. 0) then
            call dscal ( n, one/d(1), x, abs(incx) )
         else if (incd .eq. incx  .and.  incd .gt. 0) then
            do id = 1, 1 + (n - 1)*incd, incd
               x(id) = x(id) / d(id)
            end do
         else
            if (incx .ge. 0) then
               ix = 1
            else
               ix = 1 - (n - 1)*incx
            end if
            if (incd .gt. 0) then
               do id = 1, 1 + (n - 1)*incd, incd
                  x(ix) = x(ix) / d(id)
                  ix    = ix   + incx
               end do
            else
               id = 1 - (n - 1)*incd
               do i = 1, n
                  x(ix) = x(ix) / d(id)
                  id    = id + incd
                  ix    = ix + incx
               end do
            end if
         end if
      end if

      end ! subroutine dddiv

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine ddscl ( n, d, incd, x, incx )

      implicit
     &     none
      integer
     &     incd, incx, n
      double precision
     &     d(*), x(*)
*     ------------------------------------------------------------------
*     ddscl  performs the diagonal scaling  x  =  d * x.
*     ------------------------------------------------------------------
      integer
     &     i, id, ix
      external
     &     dscal
      intrinsic
     &     abs
*     ------------------------------------------------------------------
      if (n .gt. 0) then
         if (incd .eq. 0  .and.  incx .ne. 0) then
            call dscal ( n, d(1), x, abs(incx) )
         else if (incd .eq. incx  .and.  incd .gt. 0) then
            do id = 1, 1 + (n - 1)*incd, incd
               x(id) = d(id)*x(id)
            end do
         else
            if (incx .ge. 0) then
               ix = 1
            else
               ix = 1 - (n - 1)*incx
            end if
            if (incd .gt. 0) then
               do id = 1, 1 + (n - 1)*incd, incd
                  x(ix) = d(id)*x(ix)
                  ix    = ix + incx
               end do
            else
               id = 1 - (n - 1)*incd
               do i = 1, n
                  x(ix) = d(id)*x(ix)
                  id    = id + incd
                  ix    = ix + incx
               end do
            end if
         end if
      end if

      end ! subroutine ddscl

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      double precision function dnormi( n, x, incx )

      implicit
     &     none
      integer
     &     n, incx
      double precision
     &     x(*)

*     ==================================================================
*     dnormi  returns the infinity-norm of the vector  x.
*
*     29 Jul 2003: Realized that if x(*) contains NaN but not Inf,
*                  idamax may return a normal-sized entry, not the NaN.
*                  Implemented jdamax and dnormj to test more carefully.
*     ==================================================================
      external
     &     idamax
      integer
     &     idamax, kmax
*     ------------------------------------------------------------------
      double precision      zero
      parameter           ( zero = 0.0d+0 )
*     ------------------------------------------------------------------
      if (n .lt. 1) then
         dnormi = zero
      else
         kmax   = idamax( n, x, incx )
         dnormi = abs( x(kmax) )
      end if

      end ! double precision function dnormi

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      double precision function dnormj( n, x, incx )

      implicit
     &     none
      integer
     &     n, incx
      double precision
     &     x(*)

*     ==================================================================
*     dnormj returns the infinity-norm of the vector  x  in most cases.
*     flmax is returned if x(*) contains any NaNs or Infs.
*
*     29 Jul 2003: First version of dnormj for use in s5setx.
*     ==================================================================
      external
     &     jdamax, s1flmx
      integer
     &     jdamax, kmax
      double precision
     &     s1flmx
*     ------------------------------------------------------------------
      double precision      zero
      parameter           ( zero = 0.0d+0 )
*     ------------------------------------------------------------------
      kmax   = jdamax( n, x, incx )
      if (kmax .eq. 0) then
         dnormj = zero
      else if (kmax .gt. 0) then
         dnormj = abs( x(kmax) )
      else
         dnormj = s1flmx( )
      end if

      end ! double precision function dnormj

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      double precision function dnrm1s( n, x, incx )

      implicit
     &     none
      integer
     &     n, incx
      double precision
     &     x(*)

*     ==================================================================
*     dnrm1s  returns the 1-norm of the vector  x,  scaled by root(n).
*     This approximates the two-norm of x without the expense.
*     ==================================================================
      external
     &     dasum
      double precision
     &     d, dasum
*     ------------------------------------------------------------------
      double precision      zero
      parameter           ( zero = 0.0d+0 )
*     ------------------------------------------------------------------
      if (n .lt. 1) then
         dnrm1s = zero
      else
         d      = n
         d      = dasum( n, x, incx ) / sqrt(d)
         dnrm1s = d
      end if

      end ! double precision function dnrm1s

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dload ( n, const, x, incx )

      implicit
     &     none
      integer
     &     n, incx
      double precision
     &     const, x(*)
*     ==================================================================
*     dload  loads every component of a vector x with the constant.
*     Special attention is given to the case incx = 1, const = zero.
*     ==================================================================
      integer
     &     i, ix, m
*     ------------------------------------------------------------------
      double precision   zero
      parameter        ( zero = 0.0d+0 )
      intrinsic          mod
*     ------------------------------------------------------------------
      if (n .gt. 0) then
         if (const .eq. zero  .and.  incx .eq. 1) then
            m    = mod( n, 7 )
            do i = 1, m
               x(i) = zero
            end do

            do i = m+1, n, 7
               x(i)   = zero
               x(i+1) = zero
               x(i+2) = zero
               x(i+3) = zero
               x(i+4) = zero
               x(i+5) = zero
               x(i+6) = zero
            end do

         else if (const .eq. zero) then
            do ix = 1, 1+(n-1)*incx, incx
               x(ix)  = zero
            end do
         else
            do ix = 1, 1+(n-1)*incx, incx
               x(ix)  = const
            end do
         end if
      end if

      end ! subroutine dload

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine chcopy( n, x, incx, y, incy )

      implicit
     &     none

      integer
     &     n, incx, incy
      character
     &     x(*)*8, y(*)*8

*     ==================================================================
*     chcopy  is the character*8 version of dcopy.
*     ==================================================================
      integer
     &     i, ix, iy

      if (n .gt. 0) then
         if (incx .eq. incy  .and.  incy .gt. 0) then
            do iy = 1, 1 + (n - 1)*incy, incy
               y(iy) = x(iy)
            end do
         else
            if (incx .ge. 0) then
               ix = 1
            else
               ix = 1 - (n - 1)*incx
            end if
            if (incy .gt. 0) then
               do iy = 1, 1 + ( n - 1 )*incy, incy
                  y(iy) = x(ix)
                  ix    = ix + incx
               end do
            else
               iy = 1 - (n - 1)*incy
               do i  = 1, n
                  y(iy) = x(ix)
                  iy    = iy + incy
                  ix    = ix + incx
               end do
            end if
         end if
      end if

      end ! subroutine chcopy

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine icopy ( n, x, incx, y, incy )

      implicit
     &     none
      integer
     &     x(*), y(*)
      integer
     &     n, incx, incy

*     ==================================================================
*     icopy  is the integer version of dcopy.
*     ==================================================================
      integer
     &     i, ix, iy

      if (n .gt. 0) then
         if (incx .eq. incy  .and.  incy .gt. 0) then
            do iy = 1, 1 + (n - 1)*incy, incy
               y(iy) = x(iy)
            end do
         else
            if (incx .ge. 0) then
               ix = 1
            else
               ix = 1 - (n - 1)*incx
            end if
            if (incy .gt. 0) then
               do iy = 1, 1 + ( n - 1 )*incy, incy
                  y(iy) = x(ix)
                  ix    = ix + incx
               end do
            else
               iy = 1 - (n - 1)*incy
               do i  = 1, n
                  y(iy) = x(ix)
                  iy    = iy + incy
                  ix    = ix + incx
               end do
            end if
         end if
      end if

      end ! subroutine icopy

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine iload ( n, const, x, incx )

      implicit
     &     none
      integer
     &     incx, n, const, x(*)

*     iload  loads elements of x with const.

      integer
     &     ix

      if (n .gt. 0) then
         if (incx .eq. 1  .and.  const .eq. 0) then
            do ix = 1, n
               x(ix) = 0
            end do
         else
            do ix = 1, 1 + (n - 1)*incx, incx
               x(ix) = const
            end do
         end if
      end if

      end ! subroutine iload

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      integer function jdamax( n, x, incx )

      implicit
     &     none
      integer
     &     n, incx
      double precision
     &     x(*)

*     ==================================================================
*     jdamax does the same as idamax in most cases.
*     jdamax > 0 if x contains normal values.
*     jdamax = 0 if n = 0.
*     jdamax < 0 means x(-jdamax) contains the first NaN or Inf.
*
*     29 Jul 2003: First version of jdamax implemented for s5setx.
*     29 Jul 2003: Current version of jdamax
*     ==================================================================
      external
     &     s1flmx
      integer
     &     i, ix, kmax
      double precision
     &     dmax, flmax, s1flmx, xi
*     ------------------------------------------------------------------
      double precision      zero
      parameter           ( zero = 0.0d+0 )
*     ------------------------------------------------------------------

      if (n .lt. 1) then
         jdamax = 0
         return
      end if

      flmax  = s1flmx( )
      dmax   = zero
      ix     = 1
      kmax   = 1

      do i = 1, n
         xi   = abs( x(ix) )
         if (xi .le. flmax) then  ! false if xi = Nan or Inf
            if (dmax .lt. xi) then
               dmax   = xi
               kmax   = ix
            end if
         else
            go to 800
         end if
         ix = ix + incx
      end do

      jdamax = kmax
      return

  800 jdamax = -ix

      end ! integer function jdamax

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine ddrand( n, x, incx, seeds )

      implicit
     &     none
      integer
     &     n, incx, seeds(3)
      double precision
     &     x(*)

*     ------------------------------------------------------------------
*     ddrand fills a vector x with uniformly distributed random numbers
*     in the interval (0, 1) using a method due to  Wichman and Hill.
*
*     seeds(1:3) should be set to integer values
*     between 1 and 30000 before the first entry.
*
*     Integer arithmetic up to 30323 is required.
*
*     Blatantly copied from Wichman and Hill 19-January-1987.
*     14-Feb-94. Original version.
*     30 Jun 1999. seeds stored in an array.
*     30 Jun 1999. This version of ddrand.
*     ------------------------------------------------------------------
      integer
     &     ix
*     ------------------------------------------------------------------
      if (n .lt. 1) return

      do ix = 1, 1+(n-1)*incx, incx
         seeds(1)     = 171*mod(seeds(1), 177) -  2*(seeds(1)/177)
         seeds(2)     = 172*mod(seeds(2), 176) - 35*(seeds(2)/176)
         seeds(3)     = 170*mod(seeds(3), 178) - 63*(seeds(3)/178)

         if (seeds(1) .lt. 0) seeds(1) = seeds(1) + 30269
         if (seeds(2) .lt. 0) seeds(2) = seeds(2) + 30307
         if (seeds(3) .lt. 0) seeds(3) = seeds(3) + 30323

         x(ix)  = mod( real(seeds(1))/30269.0 +
     &                 real(seeds(2))/30307.0 +
     &                 real(seeds(3))/30323.0, 1.0 )
      end do

      end ! subroutine ddrand

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s1init( title, iw, leniw, rw, lenrw )

      implicit
     &     none
      character
     &     title*30
      integer
     &     leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)
*     ==================================================================
*     s1init saves some machine-dependent constants in iw and rw.
*     ==================================================================
      double precision
     &     eps, eps0, eps1, eps2, eps3, eps4, eps5, flmax, flmin,
     &     rtUndf, s1eps, s1flmx, s1flmn
*     ------------------------------------------------------------------
      eps     = s1eps ( )
      flmax   = s1flmx( )
      flmin   = s1flmn( )
      rtUndf  = sqrt(flmin)

*     Use eps to set other machine precision constants.

      eps0   = eps**(0.80d+0)
      eps1   = eps**(0.67d+0)
      eps2   = eps**(0.50d+0)
      eps3   = eps**(0.33d+0)
      eps4   = eps**(0.25d+0)
      eps5   = eps**(0.20d+0)

      rw(  1) = eps    ! machine precision
      rw(  2) = eps0   ! eps**(4/5)
      rw(  3) = eps1   ! eps**(2/3)
      rw(  4) = eps2   ! eps**(1/2)
      rw(  5) = eps3   ! eps**(1/3)
      rw(  6) = eps4   ! eps**(1/4)
      rw(  7) = eps5   ! eps**(1/5)
      rw(  8) = flmax  ! est. of the largest pos. real.
      rw(  9) = flmin  ! smallest positive real.
      rw( 10) = rtundf ! sqrt of flmin.

*     Set the environment (for later use).

      call s1envt( 0, iw, leniw )

      end ! subroutine s1init

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s1time( clock, prtopt, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     clock, prtopt, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)

*     ==================================================================
*     s1time, s1timp and s1cpu are derived from timer, timout and nowcpu
*     written for DEC VAX/VMS systems by Irvin Lustig,
*     Department of Operations Research, Stanford University, 1987.
*
*     SNOPT  calls s1time only.  s1time calls s1cpu  and  s1timp.
*     Only s1cpu is intrinsically machine dependent.
*
*     If a timer is available, call it in s1cpu  and arrange that
*     s1cpu  returns the current CPU time in seconds.
*
*     If a timer is not available or not wanted, set time = 0.0 in s1cpu.
*     Timing will be turned off and s1timp will not be called.
*     ------------------------------------------------------------------
*
*     s1time turns on or off a selected clock and optionally prints
*     statistics regarding all clocks or just the clock chosen.
*
*     The value of abs(clock) is which clock to use.
*     If clock  = 0 and prtopt = 0, all clocks and statistics are reset.
*     If clock  > 0, the clock is reset to start timing at the
*                    current time (determined by calling the
*                    machine-dependent subroutine s1cpu).
*     If clock  < 0, the clock is turned off and the statistic is
*                    recorded for the amount of time since the clock
*                    was turned on.
*
*     prtopt is the print option.
*     If lvlTim < 0, nothing is printed.  Otherwise,
*     prtopt =  0 indicates print nothing,
*            =  1 indicates print last time for this clock,
*                 only if clock < 0 (it has just been turned off),
*            =  2 indicates print total time for all clocks,
*            =  3 indicates print mean  time for all clocks.
*
*     The procedure for adding a new timer n is as follows:
*     1)  Change ntime to n in the parameter statement in s1time.
*     2)  Expand the array "label" to length n in subroutine s1timp.
*
*     04 Jun 1989: Irv's VMS/VAXC version of s1cpu installed,
*                  with changes to return time in seconds.
*     10 Jul 1992: More clocks added for use in AMPL (and elsewhere).
*     31 Jul 2003: snPRNT adopted.
*     31 Jul 2003: Current version of s1time.
*     ------------------------------------------------------------------
*
*        Clock 1 is for input time.
*        Clock 2 is for solve time.
*        Clock 3 is for output time.
*        Clock 4 is for the nonlinear functions.
*
*        numt(i)  is the number of times clock i has been turned on.
*        tlast(i) is the time at which clock i was last turned on.
*        tsum(i)  is the total time elapsed while clock i was on.
*        lvlTim   is the Timing level set in the Specs file.
*     ==================================================================
      integer
     &     ntime, lnumt, ltlast, ltsum, lvlTim
*     ------------------------------------------------------------------
      parameter         (ntime   =   5)
      parameter         (lnumt   = 451) ! 1st element of numt(10)
      parameter         (ltlast  = 451) ! 1st element of tlast(10)
      parameter         (ltsum   = 461) ! 1st element of tsum(10)
      parameter         (lvlTim  = 182) ! Timing level
*     ------------------------------------------------------------------
      call s1body( clock, prtopt, ntime, iw(lvlTim),
     &     rw(ltlast), rw(ltsum), iw(lnumt), iw, leniw )

      end ! subroutine s1time

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s1body( clock, prtopt, ntime, lvlTim,
     &     tlast, tsum, numt, iw, leniw )

      implicit
     &     none
      integer
     &     clock, prtopt, ntime, lvlTim, numt(ntime), leniw, iw(leniw)
      double precision
     &     tlast(ntime), tsum(ntime)

*     ==================================================================
*     s1body does the work for s1time.
*
*     31 Jul 2003: snPRNT adopted.
*     31 Jul 2003: Current version of s1body.
*     ==================================================================
      integer
     &     i, iStat, iclock, ilo, ihi
      double precision
     &     dtime, Stat
      real
     &     time
      external
     &     s1cpu
*     ------------------------------------------------------------------
      if (lvlTim .eq. 0) return
      iclock = abs(clock)

      if (clock .eq. 0) then
         if (prtopt .eq. 0) then

*           clock = 0, prtopt = 0.  Reset everything.

            call s1cpu ( 1, time )
            call s1cpu ( 0, time )
            do i = 1, ntime
               dtime    = dble(time)
               tlast(i) = dtime
               tsum(i)  = 0.0d+0
               numt(i)  = 0
            end do

*           If the s1cpu( 0, time ) gave time < 0.0, we assume that
*           the clock is a dummy.  Turn off future timing.

            if (time .lt. 0.0) lvlTim = 0
         end if

      else
         call s1cpu ( 0, time )
         dtime = dble(time)
         if (clock .gt. 0) then
            tlast(iclock) = dtime
         else
            Stat         = dtime - tlast(iclock)
            tsum(iclock) = tsum(iclock) + Stat
            numt(iclock) = numt(iclock) + 1
         end if
      end if

*     Now deal with print options.

      if (prtopt .eq. 0  .or.  lvlTim .lt. 0) then

*        Do nothing.

      else if (prtopt .eq. 1) then

*        Print statistic for last clock if just turned off.

         if (clock .lt. 0) then
            call s1timp( iclock, 'Last time', Stat, iw, leniw )
         end if

      else

*        prtopt >= 2.  Print all statistics if clock = 0,
*        or print statistic for individual clock.

         if (clock .eq. 0) then
            call s1cpu ( -1, time )
            ilo   = 1
            ihi   = ntime
         else
            ilo   = iclock
            ihi   = iclock
         end if

         do i = ilo, ihi
            Stat  = tsum(i)
            if (prtopt .eq. 2) then
               call s1timp( i, 'Time', Stat, iw, leniw )
            else if (prtopt .eq. 3) then
               iStat = numt(i)
               if (iStat .gt. 0) Stat = Stat / iStat
               call s1timp( i, 'Mean time', Stat, iw, leniw )
            end if
         end do
      end if

      end ! subroutine s1body

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s1timp( iclock, lstat, stat, iw, leniw )

      implicit
     &     none
      integer
     &     iclock, leniw, iw(leniw)
      double precision
     &     stat
      character*(*)
     &     lstat

*     ==================================================================
*     s1timp  prints CPU time for s1time on file iPrint and/or iSumm.
*     It is not intrinsically machine dependent.
*
*     iclock  selects the correct label.
*     lstat   is a string to print to tell which type of statistic.
*     stat    is the statistic to print out.
*             If it is zero, we figure it was never timed, so no print.
*
*     31 Jul 2003: snPRNT adopted.
*     05 Dec 2004: Replaced tab edit descriptor.
*     05 Dec 2004: Current version of s1timp.
*     ==================================================================
      character          str*60, tabby*38, label(5)*24
      data               label
     &                 / 'for MPS input',
     &                   'for solving problem',
     &                   'for solution output',
     &                   'for constraint functions',
     &                   'for objective function' /
*     ------------------------------------------------------------------
      if (iclock .eq. 1) call snPRNT( 3, ' ', iw, leniw )

      tabby = lstat(1:len(lstat)) // ' ' // label(iclock)
      write(str, 1000) tabby, stat
      call snPRNT( 3, str, iw, leniw )
      return

 1000 format( 1x, a, f13.2, ' seconds')

      end ! subroutine s1timp

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s1Perm( nkx, kx )

      implicit
     &     none
      integer
     &     nkx, kx(nkx)

*     ==================================================================
*     s1Perm  sets the default row and column order for the Jacobian.
*
*     28 Dec 2000: First version written for snopt 5.4.
*     ==================================================================
      integer
     &     j
*     ------------------------------------------------------------------
      do j = 1, nkx
         kx(j) = j
      end do

      end ! subroutine s1Perm

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s1trim( buffer, lenbuf )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     lenbuf

*     ==================================================================
*     s1trim  returns the length of buffer with trailing blanks omitted.
*
*     02 Dec 2000: First version written for snLog and snLog2.
*     ==================================================================
      integer
     &     k
*     ------------------------------------------------------------------
      lenbuf = len( buffer )
      do k = lenbuf, 2, -1
         if (buffer(k:k) .ne. ' ') go to 100
         lenbuf = lenbuf - 1
      end do

  100 return

      end ! subroutine s1trim

