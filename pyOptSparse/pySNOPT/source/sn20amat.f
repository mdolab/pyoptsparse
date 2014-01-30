*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  sn20amat.f
*
*     s2Amat   s2Aprd   s2Bprd   s2bInf   s2ColN   s2RowN   s2VarN
*     s2dInf   s2gathr  s2scatr  s2crsh   s2Mem0   s2Mem    s2rcA
*     s2scal   s2scla   s2unpk   s2vmax   s2xmat
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2Amat
     &   ( Task, lPrint, m, n, nb,
     &     nnCon, nnJac, nnObj, iObj,
     &     ne, nlocA, locA, indA, Acol,
     &     bl, bu, hrtype, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Task, lPrint, m, n, nb, ne, nlocA, nnCon, nnJac, nnObj, iObj,
     &     leniw, lenrw, hrtype(m), locA(nlocA), indA(ne), iw(leniw)
      double precision
     &     Acol(ne), bl(nb), bu(nb), rw(lenrw)

*     ==================================================================
*     s2Amat defines hrtype, the set of row types.
*
*     If Task = Rowtyp (= 0), only the row types are computed.
*     If Task = Stats  (= 1), s2Amat also prints the matrix statistics.
*
*     The vector of row-types is as follows:
*        hrtype(i) = 0  for E rows      (equalities)
*        hrtype(i) = 1  for L or G rows (inequalities)
*        hrtype(i) = 2  for N rows      (free rows)
*     They are used in s2scal and s2crsh.
*
*     15 Feb 1991: First version based on Minos 5.4 routine m2Amat.
*     04 Apr 1999: Objective stored in A.
*     31 Jul 2003: snPRNT adopted.
*     31 Jul 2003: Current version of s2Amat.
*     ==================================================================
      character
     &     str*80
      logical
     &     prnt1, summ1
      integer
     &     i, j, l, lcon, lvar, nnL,
     &     bnded, cost, fixed, free, norml
      double precision
     &     Adnsty, Aij, Aijmax, Aijmin, b1, b2, bplus, bminus, cmax,
     &     cmin, InfBnd
*     ------------------------------------------------------------------
      integer            RowTyp,          Stats
      parameter        ( Rowtyp = 0,      Stats  = 1 )
      double precision   zero,            point9
      parameter        ( zero   = 0.0d+0, point9 = 0.9d+0 )
*     ------------------------------------------------------------------
      InfBnd = rw( 70) ! definition of plus infinity.

      nnL    =   max( nnObj, nnJac )
      bplus  =   point9*InfBnd
      bminus = - bplus

      prnt1  = lPrint .ge. 1
      summ1  = lPrint .ge. 1

      if (Task .eq. RowTyp  .or.  Task .eq. Stats) then

*        Construct the vector of row-types.

         fixed = 0
         free  = 0
         norml = 0
         do  i = 1, m
            j  = n + i
            b1 = bl(j)
            b2 = bu(j)

            if (b1 .eq. b2) then
               hrtype(i) = 0
               fixed     = fixed + 1

            else if (b1 .le. bminus  .and.  b2 .ge. bplus) then
               hrtype(i) = 2
               free      = free  + 1

            else
               hrtype(i) = 1
               if (b1 .le. bminus  .or.  b2 .ge. bplus) then
                  norml  = norml+1
               end if
            end if
         end do
      end if

      if (Task .eq. Stats) then
         if (prnt1) then
            bnded = m - fixed - free - norml
            call snPRNT( 11, ' ', iw, leniw )
            call snPRNT( 11, ' Matrix statistics', iw, leniw )
            call snPRNT(  1, ' -----------------', iw, leniw )
            call snPRNT(  1, '               Total      Normal'
     &           // '        Free       Fixed     Bounded', iw, leniw )

            write(str, 2300) ' Rows   ', m, norml, free, fixed, bnded
            call snPRNT(  1, str, iw, leniw )
            fixed = 0
            free  = 0
            norml = 0

            do j  = 1, n
               b1 = bl(j)
               b2 = bu(j)
               if (b1 .eq. b2) then
                  fixed = fixed + 1
               else
                  if      (b1 .eq. zero  ) then
                     if      (b2 .ge. bplus) then
                        norml = norml + 1
                     end if
                  else if (b1 .le. bminus) then
                     if      (b2 .eq. zero ) then
                        norml = norml + 1
                     else if (b2 .ge. bplus) then
                        free  = free  + 1
                     end if
                  end if
               end if
            end do

            bnded = n - fixed - free - norml
            write(str, 2300) ' Columns', n, norml, free, fixed, bnded
            call snPRNT(  1, str, iw, leniw )

*           Find the biggest and smallest elements in a, excluding free
*           rows and fixed columns.  Also find the largest objective
*           coefficient.

            Aijmax = zero
            Aijmin = bplus
            cmax   = zero
            cmin   = bplus
            cost   = 0

            do j = 1, n
               if (bl(j) .lt. bu(j)) then
                  do l = locA(j), locA(j+1) - 1
                     i = indA(l)
                     if (hrtype(i) .eq. 2) then
                        if (i .eq. iObj) then
                           Aij   = abs( Acol(l) )
                           if (Aij .gt. zero) then
                              cost  = cost + 1
                              cmax  = max( cmax, Aij )
                              cmin  = min( cmin, Aij )
                           end if
                        end if
                     else
                        Aij    = abs( Acol(l) )
                        Aijmax = max( Aijmax, Aij )
                        Aijmin = min( Aijmin, Aij )
                     end if
                  end do
               end if
            end do

            if (Aijmin .eq. bplus) Aijmin = zero
            Adnsty = 100.0d+0*ne / (m*n)
            write(str, 2400) ne, Adnsty
            call snPRNT( 11, str, iw, leniw )
            write(str, 2410) Aijmax
            call snPRNT(  1, str, iw, leniw )
            write(str, 2420) Aijmin
            call snPRNT(  1, str, iw, leniw )
            write(str, 2430) cost
            call snPRNT( 11, str, iw, leniw )
            if (cost .gt. 0) then
               write(str, 2450) cmax
               call snPRNT(  1, str, iw, leniw )
               write(str, 2460) cmin
               call snPRNT(  1, str, iw, leniw )
            end if
         end if

*        Print a few things that can be gathered as statistics
*        from a bunch of test runs.

         if (prnt1 .or. summ1) then
            lvar  = n - nnL
            lcon  = m - nnCon
            write(str, 2500) nnCon, lcon
            if (prnt1) call snPRNT( 11, str, iw, leniw )
            if (summ1) call snPRNT( 12, str, iw, leniw )
            write(str, 2510) nnL, lvar
            if (prnt1) call snPRNT(  1, str, iw, leniw )
            if (summ1) call snPRNT(  2, str, iw, leniw )
            write(str, 2520) nnJac, nnObj
            if (prnt1) call snPRNT(  1, str, iw, leniw )
            if (summ1) call snPRNT(  2, str, iw, leniw )
            write(str, 2530) m, n
            if (prnt1) call snPRNT(  1, str, iw, leniw )
            if (summ1) call snPRNT(  2, str, iw, leniw )
         end if
      end if

      return

 2300 format(a, 5i12)
 2400 format(' No. of matrix elements', i21, 5x, 'Density', f12.3)
 2410 format(' Biggest ', 1p, e35.4, '  (excluding fixed columns,')
 2420 format(' Smallest', 1p, e35.4, '   free rows, and RHS)')
 2430 format(' No. of objective coefficients', i14)
 2450 format(' Biggest ', 1p, e35.4, '  (excluding fixed columns)')
 2460 format(' Smallest', 1p, e35.4)

 2500 format(' Nonlinear constraints', i8, 5x, 'Linear constraints', i8)
 2510 format(' Nonlinear variables  ', i8, 5x, 'Linear variables  ', i8)
 2520 format(' Jacobian  variables  ', i8, 5x, 'Objective variables',i7)
 2530 format(' Total constraints    ', i8, 5x, 'Total variables   ', i8)

      end ! subroutine s2Amat

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2Aprd
     &   ( Task, tolz,
     &     ne, nlocA, locA, indA, Acol,
     &     alpha, x, lenx, beta, y, leny )

      implicit
     &     none
      integer
     &     Task, ne, nlocA, lenx, leny, locA(nlocA), indA(ne)
      double precision
     &     tolz, alpha, beta, Acol(ne), x(lenx), y(leny)

*     ==================================================================
*     s2Aprd computes matrix-vector products involving  x  and a sparse
*     matrix A  stored by columns.
*     The parameter Task specifies the operation to be done as follows:
*       Task = Normal (=0)    y := alpha*A *x + beta*y,
*       Task = Transp (=1)    y := alpha*A'*x + beta*y,
*     where alpha and beta are scalars, x and y are vectors and A is a
*     sparse matrix stored by columns.
*
*     28 Jul 1999: Current version.
*     ==================================================================
      integer
     &     i, j, l
      double precision
     &     alphxj, sum
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

      else if (Task .eq. Normal) then
         do j = 1, lenx
            alphxj = alpha*x(j)
            if (abs( alphxj ) .gt. tolz) then
               do l = locA(j), locA(j+1)-1
                  i = indA(l)
                  if (i .le. leny) then
                     y(i) = y(i) + Acol(l)*alphxj
                  end if
               end do
            end if
         end do

      else if (Task .eq. Transp) then
         do j = 1, leny
            sum = zero
            do l = locA(j), locA(j+1)-1
               i = indA(l)
               if (i .le. lenx) then
                  sum = sum + Acol(l)*x(i)
               end if
            end do
            y(j) = y(j) + alpha*sum
         end do
      end if

      end ! subroutine s2Aprd

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2Bprd
     &   ( Task, tolz, n, lenkBS, kBS,
     &     ne, nlocA, locA, indA, Acol,
     &     alpha, x, lenx, beta, y, leny )

      implicit
     &     none
      integer
     &     Task, n, ne, nlocA, lenkBS, lenx, leny,
     &     locA(nlocA), indA(ne), kBS(lenkBS)
      double precision
     &     alpha, beta, tolz, Acol(ne), x(lenx), y(leny)

*     ==================================================================
*     s2Bprd computes various matrix-vector products involving
*     B  and  S,  the basic and superbasic columns of  A. The variable
*     Task specifies the operation to be performed as follows:
*         Task = Normal (= 0)         y := alpha*A *x + beta*y,
*         Task = Transp (= 1)         y := alpha*A'*x + beta*y,
*     where alpha and beta are scalars, x and y are vectors, and A is a
*     sparse matrix whose columns are in the basic-superbasic order
*     A = ( B  S ).
*
*     23 Nov 1991: First version of s2Bprd.
*     17 Jul 1996: Standard implementation for slacks.
*     22 Mar 1999: Reverted to integer Task control.
*     01 Apr 1999: Acol includes the objective as row 0.
*     28 Jul 1999: Current version.
*     ==================================================================
      integer
     &     i, j, k, l
      double precision
     &     alphxj, t
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
*         Relax
      else if (Task .eq. Normal) then
         do k = 1, lenx
            alphxj = alpha*x(k)
            if (abs(alphxj) .gt. tolz) then
               j   = kBS(k)
               if (j .le. n) then
*                 -------------------
*                 Column of A.
*                 -------------------
                  do l = locA(j), locA(j+1)-1
                     i = indA(l)
                     y(i) = y(i) + Acol(l)*alphxj
                  end do
               else
*                 --------------------
*                 Slack column.
*                 --------------------
                  i    = j    - n
                  y(i) = y(i) - alphxj
               end if
            end if
         end do

      else if (Task .eq. Transp) then
         do k = 1, leny
            t = zero
            j = kBS(k)

            if (j .le. n) then
*              -------------------
*              Column of A.
*              -------------------
               do l = locA(j), locA(j+1)-1
                  i = indA(l)
                  t = t + Acol(l)*x(i)
               end do
            else
*              -------------------
*              Slack column.
*              -------------------
               t   = - x(j-n )
            end if
            y(k) = y(k) + alpha*t
         end do
      end if

      end ! subroutine s2Bprd

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2bInf
     &   ( nb, bl, bu, x, bInf, jbInf )

      implicit
     &     none
      integer
     &     nb, jbInf
      double precision
     &     bInf, bl(nb), bu(nb), x(nb)

*     ==================================================================
*     s2bInf  computes the maximum infeasibility with respect to
*     the bounds on x.
*     s2bInf  is called by s5savB and s8savB before and after unscaling.
*
*     On exit,
*      bInf  is the maximum bound infeasibility.
*     jbInf  is the corresponding variable.
*
*     28 Jul 1999: First version based on Minos routine m2bInf.
*     ==================================================================
      integer
     &     j
      double precision
     &     d1, d2
*     ------------------------------------------------------------------
      double precision   zero
      parameter        ( zero = 0.0d+0 )
*     ------------------------------------------------------------------
      jbInf = 0
      bInf  = zero

      do j  = 1, nb
         d1 = bl(j) -  x(j)
         d2 =  x(j) - bu(j)
         if (bInf .lt. d1) then
             bInf  =  d1
             jbInf =  j
         end if
         if (bInf .lt. d2) then
             bInf  =  d2
             jbInf =  j
         end if
      end do

      end ! subroutine s2bInf

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      integer function s2ColN( j, leniw, iw )

      implicit
     &     none
      integer
     &     j, leniw, iw(leniw)

*     ==================================================================
*     s2ColN  gives the natural column number of Jacobian column j.
*
*     28 Dec 2000: First version written for snoptf
*     07 Jun 2001: Current version
*     ==================================================================
      integer
     &     jN, lkxN
*     ------------------------------------------------------------------
      lkxN      = iw(252) ! jN = kxN(j ) => col j of Jcol is variable iN

      if (j .gt. 0) then
         jN = iw(lkxN+j-1)
      else
         jN = 0
      end if

      s2ColN = jN

      end ! integer function s2ColN

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      integer function s2RowN( i, leniw, iw )

      implicit
     &     none
      integer
     &     i, leniw, iw(leniw)

*     ==================================================================
*     s2RowN  gives the natural row number of Jacobian row i.
*
*     28 Dec 2000: First version written for snoptf
*     07 Jun 2001: Current version
*     ==================================================================
      integer
     &     j, iN, lkxN, n
*     ------------------------------------------------------------------
      n      = iw( 15) ! copy of the number of columns
      lkxN   = iw(252) ! jN = kxN(j ) => var j is natural variable iN

      if (i .gt. 0) then
         j   = n + i
         iN  = iw(lkxN+j-1)
      else
         iN  = 0
      end if

      s2RowN = iN

      end ! integer function s2RowN

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      integer function s2VarN( j, leniw, iw )

      implicit
     &     none
      integer
     &     j, leniw, iw(leniw)

*     ==================================================================
*     s2VarN  gives the natural variable number corresponding to var i.
*
*     28 Dec 2000: First version written for snoptf
*     17 Nov 2001: Current version
*     ==================================================================
      integer
     &     jj, jN, lkxN, n
*     ------------------------------------------------------------------
      n      = iw( 15) ! copy of the number of columns
      lkxN   = iw(252) ! jN = kxN(j ) => var j is natural variable iN

      if (j .ge. 0) then
         jj =  j
      else
         jj = -j
      end if

      if (jj .gt. n) then
         jN  = n + iw(lkxN+jj-1)
      else if (jj .gt. 0) then
         jN  =     iw(lkxN+jj-1)
      else
         jN  = 0
      end if

      if (j .ge. 0) then
         s2VarN =  jN
      else
         s2VarN = -jN
      end if

      end ! integer function s2VarN

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2dInf
     &   ( n, nb, iObj, tol, bl, bu, rc, x, dInf, jdInf )

      implicit
     &     none
      integer
     &     n, nb, iObj, jdInf
      double precision
     &     dInf, tol, bl(nb), bu(nb), rc(nb), x(nb)

*     ==================================================================
*     s2dInf  computes the maximum complementarity.
*     s2dInf  is called by s4savB before and after unscaling.
*
*     On exit,
*     dInf  is the maximum dual infeasibility.
*     jdInf  is the corresponding variable.
*
*     05 Apr 1996: First version based on Minos routine m2dInf.
*     29 Nov 2002: Switched to max complementarity.
*     26 Dec 2002: This version
*     ==================================================================
      integer
     &     j, jObj
      double precision
     &     blObj, dj, tolrel
*     ------------------------------------------------------------------
      double precision   zero,           one
      parameter         (zero  = 0.0d+0, one = 1.0d+0)
*     ------------------------------------------------------------------
      jObj = n + iObj

      if (iObj .gt. 0) then
         blObj    = bl(jObj)
         bl(jObj) = bu(jObj)
      end if

      jdInf = 0
      dInf  = zero

      do j = 1, nb
         if (bl(j) .lt. bu(j)) then
            tolrel = tol*(one + abs(x(j)))
            dj     = rc(j)
            if      (x(j) .le. bl(j)+tolrel) then
               dj  = - dj
            else if (x(j) .ge. bu(j)-tolrel) then
*              dj  = + dj
            else
               dj  = abs( dj )
            end if

            if (dInf .lt. dj) then
                dInf  =  dj
                jdInf =  j
            end if
         end if
      end do

      if (iObj .gt. 0) then
         bl(jObj) = blObj
      end if

      end ! subroutine s2dInf

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2gathr
     &   ( neg, nBS, kBS, alpha, g, gBS )

      implicit
     &     none
      integer
     &     nBS, neg, kBS(nBS)
      double precision
     &     alpha, g(neg), gBS(nBS)

*     ==================================================================
*     s2gathr   performs the gather operation alpha*g  --> gBS
*     between vectors  g(neg)  and  gBS(nBS) according to
*     the index array  kBS.
*
*     The case alpha = 1 is treated specially.
*
*     On entry: 0 < nBS <= neg.
*
*     20 Jul 2000: First version of s2gathr
*     19 Dec 2000: Current version.
*     ==================================================================
      integer
     &     j, k
*     ------------------------------------------------------------------
      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one = 1.0d+0)
*     ------------------------------------------------------------------

      if (alpha .eq. one) then
         do k = 1, nBS
            j = kBS(k)
            if (j .le. neg) then
               gBS(k) =  g(j)
            else
               gBS(k) =  zero
            end if
         end do
      else if (alpha .eq. -one) then
         do k = 1, nBS
            j = kBS(k)
            if (j .le. neg) then
               gBS(k) = -g(j)
            else
               gBS(k) =  zero
            end if
         end do
      else ! general alpha
         do k = 1, nBS
            j = kBS(k)
            if (j .le. neg) then
               gBS(k) = alpha*g(j)
            else
               gBS(k) = zero
            end if
         end do
      end if

      end ! subroutine s2gathr

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2scatr
     &   ( neg, nBS, kBS, alpha, gBS, g )

      implicit
     &     none
      integer
     &     nBS, neg, kBS(nBS)
      double precision
     &     alpha, g(neg), gBS(nBS)

*     ==================================================================
*     s2scatr   performs the scatter operation alpha*gBS  --> g
*     between vectors  gBS(nBS)  and  g(neg)  according to
*     the index array  kBS.
*
*     The case alpha = 1 is treated specially.
*
*     On entry: 0 < nBS <= neg.
*
*     20 Jul 2000: First version of s2scatr
*     19 Dec 2000: Current version.
*     ==================================================================
      integer
     &     j, k
*     ------------------------------------------------------------------
      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one = 1.0d+0)
*     ------------------------------------------------------------------
      call dload ( neg, zero, g, 1 )
      if (alpha .eq. one) then
         do k = 1, nBS
            j = kBS(k)
            if (j .le. neg) then
               g(j)   =  gBS(k)
            end if
         end do
      else if (alpha .eq. -one) then
         do k = 1, nBS
            j = kBS(k)
            if (j .le. neg) then
               g(j)   = -gBS(k)
            end if
         end do
      else
         do k = 1, nBS
            j = kBS(k)
            if (j .le. neg) then
               g(j)   =  gBS(k)*alpha
            end if
         end do
      end if

      end ! subroutine s2scatr

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2crsh
     &   ( lCrash, lPrint, m, n, nb, nnCon,
     &     iCrash, tCrash,
     &     ne, nlocA, locA, indA, Acol,
     &     hpiv, hs, hrtype, bl, bu, x,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     m, n, nb, nnCon, ne, iCrash, nlocA, leniw, lenrw,
     &     locA(nlocA), indA(ne), hpiv(m), hs(nb), hrtype(m), iw(leniw)
      double precision
     &     Acol(ne), bl(nb), bu(nb), x(nb), rw(lenrw)

*     ==================================================================
*     s2crsh  looks for a basis in the columns of ( A  -I ).
*
*     ON ENTRY
*
*     iCrash    = the Crash option has been used by s5getB
*                 to set lCrash.
*     tCrash    = the Crash tolerance.  Default = 0.1
*
*     lCrash      specifies the action to be taken by Crash.
*        0,1,2,3  The call is from s4getB.
*        4,5      The call is from s8solv.
*
*        0        The all-slack basis is set up.
*
*        1        A triangular Crash is applied to the columns of A.
*                 hs(1:n) is used to help select columns.
*                 tCrash is used to ignore small entries in each column.
*                 Depending on the size of tCrash, the resulting basis
*                 will be nearly (but not strictly) lower triangular.
*
*        2        As for 1, but nonlinear rows are ignored.
*
*        3        As for 2, but linear LG rows are also ignored.
*
*        4        Linear LG rows are now included.
*                 All hs(1:nb) and x(n+i) are defined.
*                 Slack values of x(n+i) are used to select LG rows.
*
*        5        Nonlinear rows are now included.
*
*     hrtype(*)   should be defined as described in s2Amat:
*     hrtype(i) = 0  for E rows      (equalities)
*     hrtype(i) = 1  for L or G rows (inequalities)
*     hrtype(i) = 2  for N rows      (objective or free rows)
*
*     x          If lCrash <= 4, x(1:n) is used to initialize
*                 slacks as x(n+1:nb) = A*x.
*                 Used to select slacks from LG rows to be in B (basis).
*                 If lCrash  = 5, x(n+1:n+nnCon) contains slack values
*                 evaluated from x(1:n) and Fx(*).
*                 Used to select slacks from nonlinear rows to be in B.
*
*     hs          If lCrash = 1, 2 or 3, hs(1:n)  is used.
*                 If lCrash =    4 or 5, hs(1:nb) is used.
*                 If hs(j) =  0, 1 or 3, column j is eligible for B,
*                                        with 3 being "preferred".
*                 If hs(j) =  2, 4 or 5, column j is ignored.
*
*
*     Crash has several stages.
*
*     Stage 1: Insert any slacks (N, L or G rows, hrtype = 1 or 2).
*
*     Stage 2: Do triangular Crash on any free columns (wide bounds)
*
*     Stage 3: Do triangular Crash on "preferred" columns (hs(j) < 0).
*              For the linear Crash, this includes variables set
*              between their bounds in the MPS file via FR INITIAL.
*              For the nonlinear Crash, it includes nonbasics
*              between their bounds.
*              (That is, "pegged" variables in both cases.)
*
*     Stage 4: Grab unit columns.
*
*     Stage 5: Grab double columns.
*
*     Stage 6: Do triangular Crash on all columns.
*
*     Slacks are then used to pad the basis.
*
*
*     ON EXIT
*
*     hs          is set to denote an initial (B S N) partition.
*                 hs(j) = 3 denotes variables for the initial basis.
*                 If hs(j) = 2 still, variable j will be superbasic.
*                 If hs(j) = 4 or 5 still, it will be changed to 0 or 1
*                 by s4chek and variable j will be nonbasic.
*
*     x          If lCrash <= 4, slacks x(n+1:nb) are initialized.
*
*     ------------------------------------------------------------------
*        Nov 1986: Essentially the same as in 1976.
*                  Crash tolerance added.
*                  Attention paid to various hs values.
*
*     12 Nov 1988: After free rows and columns have been processed
*                  (stage 1 and 2), slacks on L or G rows are inserted
*                  if their rows have not yet been assigned.
*
*     28 Dec 1988: Refined as follows.
*                  Stage 1 inserts free and preferred rows (slacks).
*                  Stage 2 performs a triangular Crash on free or
*                          preferred columns, ignoring free rows.
*                          Unpivoted L or G slacks are then inserted.
*                  Stage 3 performs a triangular Crash on the other
*                          columns, ignoring rows whose slack is basic.
*                          (Such rows form the top part of U.  The
*                          remaining rows form the beginning of L.)
*
*     30 Apr 1989: Stage 1 now also looks for singleton columns
*                  (ignoring free and preferred rows).
*     05 May 1989: Stage 2 doesn't insert slacks if Crash option < 0.
*
*     06 Dec 1989: Stage 2, 3, 4 modified.  Columns of length 2 are
*                  now treated specially.
*
*     20 Dec 1989: Stage 2 thru 5 modified.  Free columns done before
*                  unit and double columns.
*
*     19 May 1992: x now used to help initialize slacks.
*                  Stage 1 thru 7 redefined as above.
*
*     01 Jun 1992: abs used to define closeness of slacks to bounds.
*                  Unfortunately, x(1:n) seldom has meaningful values.
*
*     02 Jun 1992: Poor performance on the larger problems.
*                  Reverted to simple approach: all slacks grabbed.
*
*     04 Jun 1992: Compromise -- Crash 3 now has 3 phases:
*                  (a) E rows.
*                  (b) LG rows.
*                  (c) Nonlinear rows.
*                  x(1:n) should then define the slack values better
*                  for (b) and (c).
*     17 Jul 1996: Standard implementation for slacks.
*     28 Jul 2003: Print level reduced to 1.
*     31 Jul 2003: snPRNT adopted.
*     31 Jul 2003: Current version of s2crsh.
*     ==================================================================
      character
     &     str*80
      logical
     &     addslk, free, prefer, gotslk,
     &     stage2, stage3, stage4, stage5, prnt1
      integer
     &     i, i1, i2, ip, j, js, k, k1, k2, lCrash, lPrint,
     &     nBasic, nPad, nPiv, nRows, nz, ipiv(2)
      double precision
     &     ai, aimax, aitol, d1, d2, eps0, tCrash, tolSlk, apiv(2)
*     ------------------------------------------------------------------
      integer            nStage
      parameter         (nStage = 6)
      integer            num(nStage), stage
      integer            Normal
      parameter         (Normal = 0)
      double precision   zero,           one
      parameter         (zero  = 0.0d+0, one = 1.0d+0)
      double precision   small,          big
      parameter         (small = 1.0d-3, big = 1.0d+4)
*     ------------------------------------------------------------------
      eps0      = rw(  2) ! eps**(4/5)
      prnt1     = lPrint .ge. 1

      if ( prnt1 ) then
         if (lCrash .le. 3) then
             write(str, '(a, i3)') ' Crash option', iCrash
             call snPRNT( 31, str, iw, leniw )
          end if
          if (lCrash .eq. 3  .and.  nnCon .lt. m) then
             call snPRNT( 31, ' Crash on linear E  rows:', iw, leniw )
          end if
          if (lCrash .eq. 4) then
             call snPRNT( 31, ' Crash on linear LG rows:', iw, leniw )
          end if
          if (lCrash .eq. 5) then
             call snPRNT( 31, ' Crash on nonlinear rows:', iw, leniw )
          end if
      end if

      if (lCrash .le. 4) then

*        Sets slacks x(n+1:nb) = A*x.
*        This is where the slacks are initialized.
*        They may be altered later (see the end of Crash).

         call s2Aprd
     &      ( Normal, eps0, ne, nlocA, locA, indA, Acol,
     &        one, x, n, zero, x(n+1), m )
      end if

*     ------------------------------------------------------------------
*     For Crash option 0, set hs(j) = 3 for all slacks and quit.
*     ------------------------------------------------------------------
      if (lCrash .eq. 0) then
         call iload ( n, 0, hs     , 1 )
         call iload ( m, 3, hs(n+1), 1 )
         go to 900
      end if

*     ------------------------------------------------------------------
*     Crash option 1, 2 or 3.   lCrash = 1, 2, 3, 4, or 5.
*     tolslk measures closeness of slacks to bounds.
*     i1,i2  are the first and last rows of A involved in Stage 1.
*     ------------------------------------------------------------------
*-->  tolslk = 0.25
      tolslk = 1.0d-2
      call iload ( nStage, 0, num, 1 )

      if (lCrash .le. 3) then
*        ---------------------------------------------------------------
*        First call.   lCrash = 1, 2 or 3.
*        Initialize hpiv(*) for all rows and hs(*) for all slacks.
*        ---------------------------------------------------------------
         i1     = 1
         if (lCrash .ge. 2) i1 = nnCon + 1
         i2     = m
         nrows  = i2 - i1 + 1

*        Make sure there are no basic columns already (hs(j) = 3).
*        If there are, make them "preferred".

         do j = 1, n
            if (hs(j) .eq. 3) hs(j) = -1
         end do

*        Make relevant rows available:  hpiv(i) = 1, hs(n+i) = 0.

         if (nrows .gt. 0) then
            call iload ( nrows, 1, hpiv(i1), 1 )
            call iload ( nrows, 0, hs(n+i1), 1 )
         end if

         if (lCrash .eq. 1) then
            nbasic = 0
         else

*           lCrash = 2 or 3:  Insert nonlinear slacks.

            nbasic = nnCon
            if (nnCon .gt. 0) then
               call iload ( nnCon, 3, hpiv   , 1 )
               call iload ( nnCon, 3, hs(n+1), 1 )
            end if
         end if

         if (lCrash .eq. 3) then

*           Insert linear inequality slacks (including free rows).

            do i = i1, m
               if (hrtype(i) .ge. 1) then
                  nbasic  = nbasic + 1
                  nrows   = nrows  - 1
                  hpiv(i) = 3
                  hs(n+i) = 3
               end if
            end do
         end if

*        We're done if there are no relevant rows.

         if (nrows .eq. 0) go to 800

      else
*        ---------------------------------------------------------------
*        Second or third call.  lCrash = 4 or 5.
*        Initialize hpiv(*) for all rows.
*        hs(*) already defines a basis for the full problem,
*        but we want to do better by including only some of the slacks.
*        ---------------------------------------------------------------
         if (lCrash .eq. 4) then
*           ------------------------------------------------------------
*           Crash on linear LG rows.
*           ------------------------------------------------------------
            if (nnCon .eq. m) go to 900
            i1     = nnCon + 1
            i2     = m

*           Mark nonlinear rows as pivoted: hpiv(i) = 3.

            nbasic = nnCon
            if (nbasic .gt. 0) then
               call iload ( nbasic, 3, hpiv, 1 )
            end if

*           Mark linear E  rows as pivoted: hpiv(i) = 3
*           Make linear LG rows available:  hpiv(i) = 1, hs(n+i) = 0.

            do i = i1, m
               if (hrtype(i) .eq. 0) then
                  nbasic  = nbasic + 1
                  hpiv(i) = 3
               else
                  hpiv(i) = 1
                  hs(n+i) = 0
               end if
            end do

*           Mark linear LG rows with hpiv(i) = 2
*           if any basic columns contain a nonzero in row i.

            do j = 1, n
               if (hs(j) .eq. 3) then
                  do k = locA(j), locA(j+1) - 1
                     i     = indA(k)
                     if (hrtype(i) .eq. 1) then
                        if (i .gt. nnCon) then
                           if (Acol(k) .ne. zero) hpiv(i) = 2
                        end if
                     end if
                  end do
               end if
            end do
         else
*           ------------------------------------------------------------
*           lCrash = 5.  Crash on nonlinear rows.
*           ------------------------------------------------------------
            i1     = 1
            i2     = nnCon

*           Mark all linear rows as pivoted: hpiv(i) = 3

            nbasic = m - nnCon
            if (nbasic .gt. 0) then
               call iload ( nbasic, 3, hpiv(nnCon+1), 1 )
            end if

*           Make nonlinear rows available:  hpiv(i) = 1, hs(n+i) = 0.

            call iload ( nnCon, 1, hpiv   , 1 )
            call iload ( nnCon, 0, hs(n+1), 1 )

*           Mark nonlinear rows with hpiv(i) = 2
*           if any basic columns contain a nonzero in row i.

            do j = 1, n
               if (hs(j) .eq. 3) then
                  do k = locA(j), locA(j+1) - 1
                     i = indA(k)
                     if (i .le. nnCon) then
                        if (Acol(k) .ne. zero) hpiv(i) = 2
                     end if
                  end do
               end if
            end do
         end if
      end if

*     ------------------------------------------------------------------
*     Stage 1: Insert relevant slacks (N, L or G rows, hrtype = 1 or 2).
*              If lCrash = 4 or 5, grab them only if they are more than
*              tolslk from their bound.
*     ------------------------------------------------------------------
      stage  = 1
      gotslk = lCrash .eq. 4  .or.  lCrash .eq. 5

      do i = i1, i2
         j = n + i
         if (hs(j) .le. 1  .and.  hrtype(i) .gt. 0) then
            addslk = .true.

            if (gotslk) then
               d1 = x(j) - bl(j)
               d2 = bu(j) - x(j)
               if (min( d1, d2 ) .le. tolslk) then

*                 The slack is close to a bound or infeasible.
*                 Move it exactly onto the bound.

                  addslk = .false.

                  if (d1 .le. d2) then
                     x(j)  = bl(j)
                     hs(j) = 0
                  else
                     x(j)  = bu(j)
                     hs(j) = 1
                  end if
               end if
            end if

            if (addslk) then
               nbasic     = nbasic     + 1
               num(stage) = num(stage) + 1
               hpiv(i)    = 3
               hs(j)      = 3
            end if
         end if
      end do

      if (nbasic .eq. m) go to 700

*     ------------------------------------------------------------------
*     Apply a triangular Crash to various subsets of the columns of A.
*
*        hpiv(i) = 1  if row i is unmarked (initial state).
*        hpiv(i) = 3  if row i has been given a pivot
*                     in one of a set of triangular columns.
*        hpiv(i) = 2  if one of the triangular columns contains
*                     a nonzero in row i below the triangle.
*     ------------------------------------------------------------------
      do stage = 2, nStage
         stage2     = stage .eq. 2
         stage3     = stage .eq. 3
         stage4     = stage .eq. 4
         stage5     = stage .eq. 5

*        ---------------------------------------------------------------
*        Main loop for triangular Crash.
*        ---------------------------------------------------------------
         do 200,  j = 1, n
            js      = hs(j)
            if (js    .gt.    1 ) go to 200
            if (bl(j) .eq. bu(j)) go to 200

            if ( stage2 ) then
               free   = bl(j) .le. - big  .and.  bu(j) .ge. big
               if ( .not. free  ) go to 200

            else if ( stage3 ) then
               prefer = js .lt. 0
               if ( .not. prefer) go to 200
            end if

*           Find the biggest aij, ignoring free rows.

            k1     = locA(j)
            k2     = locA(j+1) - 1
            aimax  = zero

            do k = k1, k2
               i = indA(k)
               if (hrtype(i) .ne. 2) then
                  ai    = Acol(k)
                  aimax = max( aimax, abs( ai ) )
               end if
            end do

*           Prevent small pivots if Crash tol is too small.

            if (aimax .le. small) go to 200

*           Find the biggest pivots in rows that are still
*           unpivoted and unmarked.  Ignore smallish elements.
*           nz counts the number of relevant nonzeros.

            aitol   = aimax * tCrash
            nz      = 0
            npiv    = 0
            ipiv(1) = 0
            ipiv(2) = 0
            apiv(1) = zero
            apiv(2) = zero

            do k = k1, k2
               i = indA(k)
               if (hs(n+i) .ne. 3) then
                  ai = abs( Acol(k) )
                  if (ai .gt. aitol) then
                     nz = nz + 1
                     ip = hpiv(i)
                     if (ip .le. 2) then
                        if (apiv(ip) .lt. ai) then
                           apiv(ip) = ai
                           ipiv(ip) = i
                        end if
                     else
                        npiv = npiv + 1
                     end if
                  end if
               end if
            end do

*           Grab unit or double columns.

            if      ( stage4 ) then
               if (nz .ne. 1) go to 200
            else if ( stage5 ) then
               if (nz .ne. 2) go to 200
            end if

*           See if the column contained a potential pivot.
*           An unmarked row is favored over a marked row.

            ip     = 1
            if (ipiv(1)  .eq.  0  .and.  npiv .eq. 0) ip = 2
            i      = ipiv(ip)

            if (i .gt. 0) then
               nbasic     = nbasic     + 1
               num(stage) = num(stage) + 1
               hpiv(i)    = 3
               hs(j)      = 3
               if (nbasic .ge. m) go to 700

*              Mark off all relevant unmarked rows.

               do k = k1, k2
                  i = indA(k)
                  if (hs(n+i) .ne. 3) then
                     ai = abs( Acol(k) )
                     if (ai .gt. aitol) then
                        if (hpiv(i) .eq. 1) hpiv(i) = 2
                     end if
                  end if
               end do
            end if
  200    continue
      end do

*     ------------------------------------------------------------------
*     All stages finished.
*     Fill remaining gaps with slacks.
*     ------------------------------------------------------------------
  700 npad   = m - nbasic
      if ( prnt1 ) then
         write(str, 1200) num(1), num(2), num(3)
         call snPRNT( 21, str, iw, leniw )
         write(str, 1210) num(4), num(5), num(6), npad
         call snPRNT( 21, str, iw, leniw )
      end if

      if (npad .gt. 0) then
         do i = 1, m
            if (hpiv(i) .lt. 3) then
               nbasic  = nbasic + 1
               hs(n+i) = 3
               if (nbasic .ge. m) go to 800
            end if
         end do
      end if

*     ------------------------------------------------------------------
*     Make sure there aren't lots of nonbasic slacks floating off
*     their bounds.  They could take lots of iterations to move.
*     ------------------------------------------------------------------
  800 do i = i1, i2
         j = n + i
         if (hs(j) .le. 1  .and.  hrtype(i) .gt. 0) then
            d1 = x(j)  - bl(j)
            d2 = bu(j) - x(j)
            if (min( d1, d2 ) .le. tolslk) then

*              The slack is close to a bound or infeasible.
*              Move it exactly onto the bound.

               if (d1 .le. d2) then
                  x(j)  = bl(j)
                  hs(j) = 0
               else
                  x(j)  = bu(j)
                  hs(j) = 1
               end if
            end if
         end if
      end do

  900 return

 1200 format(' Slacks', i6, '  Free cols', i6, '  Preferred', i6)
 1210 format(' Unit  ', i6, '  Double   ', i6, '  Triangle ', i6,
     &     '  Pad', i6)

      end ! subroutine s2crsh

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2Mem0
     &   ( iExit, Solver, lencw, leniw, lenrw, iw,
     &     mincw, miniw, minrw, maxcw, maxiw, maxrw,
     &     nextcw, nextiw, nextrw )

      implicit
     &     none
      character
     &     Solver*6
      integer
     &     iExit, lencw, leniw, lenrw, mincw, miniw, minrw,
     &     maxcw, maxiw, maxrw, nextcw, nextiw, nextrw, iw(leniw)

*     ==================================================================
*     s2Mem0   checks the memory requirements for sqopt/snopt and sets
*     the pointers to the beginning and end of sqopt/snopt workspace.
*
*     Note: cw, iw and rw hold constants and work-space addresses.
*        They must have dimension at least 500.
*
*     The SPECS file has been read, and values are known for
*     maxcu, maxiu, maxru  (upper limit of user  partition 1)
*     maxcw, maxiw, maxrw  (upper limit of SNOPT partition)
*
*     The default values for these values are
*     maxcu = 500  ,   maxiu = 500  ,   maxru = 500,
*     maxcw = lencw,   maxiw = leniw,   maxrw = lenrw,
*     which are set in snInit:
*
*     The user can alter these in the SPECS file via
*     lines of the form
*
*        User  character workspace      10000    (Sets maxcu)
*        User  integer   workspace      10000    (Sets maxiu)
*        User  real      workspace      10000    (Sets maxru)
*        Total character workspace      90000    (Sets maxcw)
*        Total integer   workspace      90000    (Sets maxiw)
*        Total real      workspace      90000    (Sets maxrw)
*
*     SNOPT will use only rw(1:500) and rw(maxru+1:maxrw).
*     Hence, rw(501:maxru) and possibly rw(maxrw+1:lenrw) may be used as
*     workspace by the user during solution of the problem (e.g., within
*     the user-supplied function routines).  Similarly for the integer
*     and character work arrays iw and cw.
*
*     Setting maxiw and maxrw less than leniw and lenrw may serve to
*     reduce paging activity on a machine with virtual memory, by
*     confining SNOPT (in particular the LU-factorization routines) to
*     an area of memory that is sensible for the current problem.  This
*     often allows cw(*), iw(*) and rw(*) to be declared arbitrarily
*     large at compile time.
*
*     15 Nov 1991: First version based on Minos 5.4 routine m2core.
*     29 Mar 1998: First version called by snMem. This simplified
*                  version may slightly overestimate needed memory.
*     27 Apr 1999: Renamed s8Mem.
*     16 Jul 2003: Standard output added
*     31 Jul 2003: snEXIT and snPRNT adopted.
*     27 Oct 2003: Renamed s2Mem0.
*     15 Oeb 2004: Current version of s2Mem0.
*     ==================================================================
      character
     &     str*80, str2*80
      logical
     &     fixdup
      integer
     &     maxcu , maxiu , maxru , mxcu  , mxiu  , mxru,
     &     mxcw  , mxiw  , mxrw,
     &     mincu1, miniu1, minru1, maxcu1, maxiu1, maxru1,
     &     mincu2, miniu2, minru2, maxcu2, maxiu2, maxru2
*     ------------------------------------------------------------------
      iExit  = 0

      if (lencw .lt. 500 .or. leniw .lt. 500 .or. lenrw .lt. 500) then
*        ---------------------------------------------------------------
*        Not enough workspace to do ANYTHING!
*        Print and exit without accessing the work arrays.
*        ---------------------------------------------------------------
         mincw = 500
         miniw = 500
         minrw = 500
         iExit = 81        ! Work arrays must have at least 500 elements
         call snWRAP( iExit, Solver, str, str2, iw, leniw )
         go to 999
      end if

      mxcu      = iw(  6) ! maxcu+1 is the start of SNOPT part of cw
      mxiu      = iw(  4) ! maxiu+1 is the start of SNOPT part of iw
      mxru      = iw(  2) ! maxru+1 is the start of SNOPT part of rw

      mxcw      = iw(  7) ! end of SNOPT part of cw
      mxiw      = iw(  5) ! end of SNOPT part of iw
      mxrw      = iw(  3) ! end of SNOPT part of rw

*     Check for silly values.

      maxcu  = min( max( mxcu, 500 ), lencw )
      maxiu  = min( max( mxiu, 500 ), leniw )
      maxru  = min( max( mxru, 500 ), lenrw )

      maxcw  = min( max( mxcw, 500 ), lencw )
      maxiw  = min( max( mxiw, 500 ), leniw )
      maxrw  = min( max( mxrw, 500 ), lenrw )

      maxcu  = min( maxcu, maxcw )
      maxiu  = min( maxiu, maxiw )
      maxru  = min( maxru, maxrw )

      fixdup = mxcu .ne. maxcu  .or.  mxcw .ne. maxcw  .or.
     &         mxiu .ne. maxiu  .or.  mxiw .ne. maxiw  .or.
     &         mxru .ne. maxru  .or.  mxrw .ne. maxrw

      if ( fixdup ) then
         call snPRNT( 14,
     &        ' XXX  User workspace parameters had to be modified',
     &        iw, leniw )
      end if

*     Save the checked values

      iw(  6) = maxcu  ! maxcu+1 is the start of SNOPT part of cw
      iw(  4) = maxiu  ! maxiu+1 is the start of SNOPT part of iw
      iw(  2) = maxru  ! maxru+1 is the start of SNOPT part of rw

      iw(  7) = maxcw  ! end of SNOPT part of cw
      iw(  5) = maxiw  ! end of SNOPT part of iw
      iw(  3) = maxrw  ! end of SNOPT part of rw

*     ------------------------------------------------------------------
*     Save the limits of the two user-accessible workspace partitions
*     to allow the user can grab them from iw.
*
*     Upper limits of partition 1 may be set in the specs file.
*     lower limits of partition 2 are set here.
*     ------------------------------------------------------------------
      mincu1  = 501             ! Lower limits on partition 1
      miniu1  = 501
      minru1  = 501

      maxcu1  = maxcu           ! User-defined upper limits on
      maxiu1  = maxiu           ! partition 1 (default 500).
      maxru1  = maxru

      maxcu2  = lencw           ! Upper limits on partition 2
      maxiu2  = leniw
      maxru2  = lenrw


      iw( 31) = mincu1          ! Start of first  user partition of cw
      iw( 36) = miniu1          ! Start of first  user partition of iw
      iw( 41) = minru1          ! Start of first  user partition of rw

      iw( 32) = maxcu1          ! End   of first  user partition of cw
      iw( 37) = maxiu1          ! End   of first  user partition of iw
      iw( 42) = maxru1          ! End   of first  user partition of rw

      iw( 34) = maxcu2          ! End   of second user partition of cw
      iw( 39) = maxiu2          ! End   of second user partition of iw
      iw( 44) = maxru2          ! End   of second user partition of rw

      mincu2  = maxcw + 1       ! Lower limits on partition 2
      miniu2  = maxiw + 1
      minru2  = maxrw + 1

      iw( 33) = mincu2          ! Start of second user partition of cw
      iw( 38) = miniu2          ! Start of second user partition of iw
      iw( 43) = minru2          ! Start of second user partition of rw

*     ------------------------------------------------------------------
*     Define the start of the character, integer and real workspace.
*     ------------------------------------------------------------------
      nextcw  = maxcu + 1
      nextiw  = maxiu + 1
      nextrw  = maxru + 1

  999 return

      end ! subroutine s2Mem0

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2Mem
     &   ( iExit, PrtMem, liwEst, lrwEst,
     &     nextcw, nextiw, nextrw,
     &     maxcw, maxiw, maxrw, lencw, leniw, lenrw,
     &     mincw, miniw, minrw, iw )

      implicit
     &     none
      logical
     &     PrtMem
      integer
     &     iExit, lencw, leniw, lenrw, liwEst, lrwEst,
     &     maxcw, maxiw, maxrw, mincw, miniw, minrw,
     &     nextcw, nextiw, nextrw, iw(leniw)

*     ==================================================================
*     s2Mem  prints details of the workspace.
*
*     On exit.
*        If iExit = 0,  mincw, miniw, minrw give the amounts of
*        character, integer and real workspace needed to hold
*        the problem. (The LU factorization routines may
*        subsequently ask for more.)
*
*        If iExit > 0,  insufficient storage is provided to hold the
*        problem.  In this case, mincw, miniw and minrw give estimates
*        of reasonable lengths for cw(*), iw(*) and rw(*).
*
*     15 Nov 1991: First version based on Minos 5.4 routine m2core.
*     29 Mar 1998: First version called by snMem. This simplified
*                  version may slightly overestimate needed memory.
*     03 Aug 2003: snEXIT and snPRNT adopted.
*     18 Feb 2004: Current version of s2Mem.
*     ==================================================================
      character
     &     str*132
      integer
     &     maxcu1, maxiu1, maxru1, mincu1, miniu1, minru1, lenALU
*     ------------------------------------------------------------------
      iExit   = 0

*     Compute the minimum storage required

      mincw   = nextcw - 1
      miniw   = nextiw - 1
      minrw   = nextrw - 1

*     ------------------------------------------------------------------
*     Print details of the workspace.
*     ------------------------------------------------------------------
      mincu1  = iw( 31)         ! Start of first  user partition of cw
      miniu1  = iw( 36)         ! Start of first  user partition of iw
      minru1  = iw( 41)         ! Start of first  user partition of rw

      maxcu1  = iw( 32)         ! End   of first  user partition of cw
      maxiu1  = iw( 37)         ! End   of first  user partition of iw
      maxru1  = iw( 42)         ! End   of first  user partition of rw

      if (PrtMem) then
         write(str, 1110) maxcw, maxiw, maxrw
         call snPRNT( 31, str, iw, leniw )
         write(str, 1120) mincw,  miniw, minrw
         call snPRNT( 21, str, iw, leniw )
         call snPRNT( 21, ' ', iw, leniw )
         if (maxcu1 .ge. mincu1) then
            write(str, 1201) mincu1, maxcu1
            call snPRNT( 21, str, iw, leniw )
         end if
         if (maxiu1 .ge. miniu1) then
            write(str, 1202) miniu1, maxiu1
            call snPRNT( 21, str, iw, leniw )
         end if
         if (maxru1 .ge. minru1) then
            write(str, 1203) minru1, maxru1
            call snPRNT( 21, str, iw, leniw )
         end if
      end if

      if (mincw .gt. maxcw   .or.  miniw .gt. maxiw
     &                       .or.  minrw .gt. maxrw) then
*        ---------------------------------------------------------------
*        Not enough workspace to solve the problem.
*        ---------------------------------------------------------------
         if (PrtMem) then
            call snPRNT( 3,
     &           ' XXX  Not enough storage to start the problem...',
     &           iw, leniw )
         end if
      end if

      if (mincw .gt. lencw) then
         iExit = 82             ! Not enough character workspace.
         if (PrtMem) then
            write(str, 9420) mincw
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (miniw .gt. leniw) then
         iExit = 83             ! Not enough integer workspace.
         miniw = liwEst
         if (PrtMem) then
            write(str, 9430) miniw
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

      if (minrw .gt. lenrw) then
         iExit = 84             ! Not enough real workspace.
         minrw = lrwEst
         if (PrtMem) then
            write(str, 9440) minrw
            call snPRNT( 13, str, iw, leniw )
         end if
      end if

*     LUSOL stores  indc(*), indr(*) and A(*) in  iw(miniw:maxiw) and
*     rw(minrw:maxrw), with miniw pointing to the start of
*     indc(*), indr(*), and minrw pointing to the start of A(*).

      lenALU  = iw(213)

      if (iExit .eq. 0  .and.  lenALU .eq. 0) then
*        -------------------------------------------
*        Insufficient storage to factorize B.
*        -------------------------------------------
         iExit = 82
         miniw = liwEst
         minrw = lrwEst
         if (PrtMem) then
            call snPRNT( 3,
     &         ' XXX  Not enough storage for the basis factors',
     &         iw, leniw )

            write(str, 9500)
            call snPRNT( 13, str, iw, leniw )
            write(str, 9510) maxiw, liwEst
            call snPRNT(  3, str, iw, leniw )
            write(str, 9520) maxrw, lrwEst
            call snPRNT(  3, str, iw, leniw )
         end if
      end if

      return

 1110 format(' Total char*8  workspace', i10, 6x,
     &       ' Total integer workspace', i10, 6x,
     &       ' Total real    workspace', i10)
 1120 format(' Total char*8  (minimum)', i10, 6x,
     &       ' Total integer (minimum)', i10, 6x,
     &       ' Total real    (minimum)', i10)
 1201 format(' Elements cw(', i10, ':',i10, ')', 6x, 'are free',
     &       ' for USER CHAR*8  WORKSPACE')
 1202 format(' Elements iw(', i10, ':',i10, ')', 6x, 'are free',
     &       ' for USER INTEGER WORKSPACE')
 1203 format(' Elements rw(', i10, ':',i10, ')', 6x, 'are free',
     &       ' for USER REAL    WORKSPACE')

 9420 format(' Total character workspace should be significantly',
     &       ' more than', i8)
 9430 format(' Total integer   workspace  should be significantly',
     &       ' more than', i8)
 9440 format(' Total real      workspace  should be significantly',
     &       ' more than', i8)
 9500 format(24x, '        Current    Recommended')
 9510 format(' Total integer workspace', 2i15)
 9520 format(' Total real    workspace', 2i15)

      end ! subroutine s2Mem

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2rcA
     &   ( feasbl, featol, iObj, minimz, wtInf,
     &     m, n, nb, leng, neg,
     &     ne, nlocA, locA, indA, Acol,
     &     hEstat, hs, bl, bu, g, pi, rc, x )

      implicit
     &     none
      logical
     &     feasbl
      integer
     &     iObj, minimz, m, n, nb, ne, nlocA, leng, neg,
     &     hEstat(nb), hs(nb), locA(nlocA), indA(ne)
      double precision
     &     featol, wtInf,
     &     Acol(ne), bl(nb), bu(nb), g(leng), rc(nb), x(nb), pi(m)

*     ==================================================================
*     s2rcA  computes reduced costs rc(*) for all columns of ( A  -I ).
*     If x is feasible, the gradient includes the explicit objective
*     vector  g.  Otherwise, the phase 1 objective is included.
*
*     s2rcA  is called by s4savb and s8SQP for scaled data.
*     External values of hs(*) are used (0, 1, 2, 3),
*     but internal would be ok too since we only test for > 1.
*
*     19 Feb 1994: First version based on Minos routine m4rc.
*     17 Jul 1996: Standard implementation for slacks.
*     31 Jul 1999: Current version.
*     ==================================================================
      integer
     &     i, j, jObj, l
      double precision
     &     d1, d2, dj, sgnObj
*     ------------------------------------------------------------------
      double precision   zero,           one
      parameter        ( zero = 0.0d+0,  one = 1.0d+0 )
*     ------------------------------------------------------------------
      sgnObj = minimz

      do j  = 1, n
         dj = zero
         do l  = locA(j), locA(j+1)-1
            i  = indA(l)
            dj = dj  +  pi(i) * Acol(l)
         end do
         rc(j) =  - dj
      end do

      call dcopy ( m, pi, 1, rc(n+1), 1 )

      if ( feasbl ) then

*        Include the nonlinear objective gradient.
*        Include the gradient of the linear term.

         sgnObj = minimz

         if (neg .gt. 0) then
            call daxpy ( neg, sgnObj, g, 1, rc, 1 )
         end if

         if (iObj .gt. 0) then
            jObj     = n + iObj
            rc(jObj) = rc(jObj) + sgnObj
         end if

      else

*        Include the Phase 1 objective.
*        Only basics and superbasics can be infeasible.
*        Check that this works for scaling.

         do j = 1, nb
            if (hs(j) .gt. 1) then
               d1  = bl(j) - x(j)
               d2  = x(j) - bu(j)
               if (hEstat(j) .eq. 0) then
                  if (d1 .gt. featol) rc(j) = rc(j) - one
                  if (d2 .gt. featol) rc(j) = rc(j) + one
               else
                  if (d1 .gt. featol) rc(j) = rc(j) - wtInf
                  if (d2 .gt. featol) rc(j) = rc(j) + wtInf
               end if
            end if
         end do
      end if

      end ! subroutine s2rcA

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2scal
     &   ( lPrint, m, n, nb, nnL, nnCon, nnJac, hrtype,
     &     ne, nlocA, locA, indA, Acol,
     &     Ascale, bl, bu, rmin, rmax,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     m, n, nb, nlocA, nnL, nnCon, nnJac, ne, leniw, lenrw,
     &     hrtype(m), locA(nlocA), indA(ne), iw(leniw)
      double precision
     &     Acol(ne), Ascale(nb), bl(nb), bu(nb), rmin(m), rmax(m),
     &     rw(lenrw)

*     ==================================================================
*     s2Scal computes scale factors  AScale  from A, bl, bu.
*
*     In phase 1, an iterative procedure based on geometric means is
*     used to compute scales from a alone.  This procedure is derived
*     from a routine written by Robert Fourer, 1979.  The main steps
*     are:
*
*        (1) Compute Aratio = max(i1,i2,j)  A(i1,j) / A(i2,j).
*        (2) Divide each row i by
*               ( min(j) A(i,j) * max(j) A(i,j) ) ** 1/2.
*        (3) Divide each column j by
*               ( min(i) A(i,j) * max(i) A(i,j) ) ** 1/2.
*        (4) Compute sratio as in (1).
*        (5) If sratio .lt. scltol * Aratio, repeat from step (1).
*
*        Free rows (hrtype=2) and fixed columns (bl=bu) are not used
*        at this stage.
*
*     In phase 2, the scales for free rows are set to be their largest
*     element.
*
*     In phase 3, fixed columns are summed in order to compute
*     a scale factor sigma that allows for the effective rhs of the
*     constraints.  All scales are then multiplied by sigma.
*
*
*     If lvlScl = 1, the first nnCon rows and the first nnL columns will
*     retain scales of 1.0 during phases 1-2, and phase 3 will not be
*     performed.  (This takes effect if the problem is nonlinear but
*     the user has specified 'scale linear variables' only.)
*     However, all rows    contribute to the linear column scales,
*     and      all columns contribute to the linear row    scales.
*
*     If lvlScl = 2, all rows and columns are scaled.  To guard against
*     misleadingly small Jacobians, if the maximum element in any of
*     the first nnCon rows or the first nnL columns is less than
*     smallj, the corresponding row or column retains a scale of 1.0.
*
*     15 Nov 1991: First version based on Minos 5.4 routine m2scal.
*     24 Jan 2001: Don't scale BIG free rows anymore!
*     31 Jul 2003: snPRNT adopted.
*     31 Jul 2003: Current version of s2scal.
*     ==================================================================
      external
     &     dnormi
      character
     &     str*80
      logical
     &     lonly, Prnt1
      integer
     &     i, imax, imin, iStart, j, jmax, jmin,
     &     jStart, l, lPrint, lprScl, lvlScl, mxPass, nclose, npass
      double precision
     &     Ac, Ar, Amin, Amax, Aratio, big, bnd, b1, b2,
     &     cmin, cmax, cratio, sigma, small, sratio, close1, close2,
     &     InfBnd,tolx, scltol, bplus,
     &     dnormi
*     ------------------------------------------------------------------
      double precision   zero,           one
      parameter         (zero = 0.0d+0,  one    = 1.0d+0)
      double precision   damp,           smallj
      parameter         (damp = 1.0d-4,  smallj = 1.0d-2 )
*     ------------------------------------------------------------------
      tolx      = rw( 56) ! Minor feasibility tolerance.
      InfBnd    = rw( 70) ! definition of plus infinity.

      lvlScl    = iw( 75) ! scale option
      lprScl    = iw( 83) ! > 0    => print the scales
      scltol    = rw( 92) ! scale tolerance.

      Prnt1     = lPrint .ge. 10
      if ( Prnt1 ) then
         call snPRNT( 31, ' '       , iw, leniw )
         call snPRNT( 21, ' Scaling', iw, leniw )
         call snPRNT( 21, ' -------', iw, leniw )
         call snPRNT( 21, '             '
     &        // 'Min elem    Max elem       Max col ratio', iw, leniw )
      end if
      bplus     = 0.1d+0*InfBnd
      Aratio    = bplus
      mxPass    = 10
      lonly     = lvlScl .eq. 1

      call dload ( nb, one, Ascale, 1 )

      if ( lonly ) then
         if (nnL .ge. n) return
         iStart = nnCon + 1
         jStart = nnL + 1
      else
         iStart = 1
         jStart = 1
      end if

*     ------------------------------------------------------------------
*     Main loop for phase 1.
*     Only the following row-types are used:
*        hrtype(i) = 2       for type N rows (objective or free rows),
*        hrtype(i) = 0 or 1  otherwise.
*     ------------------------------------------------------------------
      do nPass = 0, mxPass

*        Find the largest column ratio.
*        Also set new column scales (except on pass 0).

         Amin    = bplus
         Amax    = zero
         small   = smallj
         sratio  = one

         do j = jStart, n
            if (bl(j) .lt. bu(j)) then
               cmin = bplus
               cmax = zero

               do l = locA(j), locA(j+1)-1
                  i = indA(l)
                  if (hrtype(i) .ne. 2) then
                     Ar = abs( Acol(l) )

                     if (Ar .gt. zero) then
                        Ar   = Ar / Ascale(n+i)
                        cmin = min( cmin, Ar )
                        cmax = max( cmax, Ar )
                     end if
                  end if
               end do

               Ac     = max( cmin, damp*cmax )
               Ac     = sqrt( Ac ) * sqrt( cmax )
               if (j     .gt. nnJac) small     = zero
               if (cmax  .le. small) Ac        = one
               if (nPass .gt. 0    ) Ascale(j) = Ac
               Amin   = min( Amin, cmin / Ascale(j) )
               Amax   = max( Amax, cmax / Ascale(j) )
               cratio = cmax / cmin
               sratio = max( sratio, cratio )
            end if
         end do

         if ( Prnt1 ) then
            write(str, 1200) nPass, Amin, Amax, sratio
            call snPRNT( 21, str, iw, leniw )
         end if

*        ---------------------------------------------------------------
*        Test for convergence.
*        ---------------------------------------------------------------
         if (nPass  .ge. 3       .and.
     &       sratio .ge. Aratio*scltol) go to 420  ! Break

         if (nPass .lt. mxPass) then

            Aratio = sratio

*           Set new row scales for the next pass.

            if (iStart .le. m) then

               if (iStart .le. m) then
                  call dload ( m-iStart+1, bplus, rmin(iStart), 1 )
                  call dload ( m-iStart+1,  zero, rmax(iStart), 1 )
               end if

               do j = 1, n
                  if (bl(j) .lt. bu(j)) then
                     Ac  = Ascale(j)

                     do l = locA(j), locA(j+1)-1
                        i = indA(l)

                        if (i .ge. iStart) then
                           Ar  = abs( Acol(l) )

                           if (Ar .gt. zero) then
                              Ar      = Ar / Ac
                              rmin(i) = min( rmin(i), Ar )
                              rmax(i) = max( rmax(i), Ar )
                           end if
                        end if
                     end do
                  end if
               end do

               do i  = iStart, m
                  j  = n + i
                  Ar = rmax(i)

                  if (i .le. nnCon  .and.  Ar .le. smallj) then
                     Ascale(j) = one
                  else
                     Ac        = max( rmin(i), damp*Ar )
                     Ascale(j) = sqrt( Ac ) * sqrt( Ar )
                  end if
               end do
            end if
         end if
      end do

*     ------------------------------------------------------------------
*     End of main loop.
*     ------------------------------------------------------------------

*     Invert the column scales, so that structurals and logicals
*     can be treated the same way during subsequent unscaling.
*     Find the min and max column scales while we're at it.
*     Nov 1989: nclose counts how many are "close" to 1.
*     For problems that are already well-scaled, it seemed sensible to
*     set the "close" ones exactly equal to 1.
*     Tried "close" = (0.5,2.0) and (0.9,1.1), but they helped only
*     occasionally.  Have decided not to interfere.

  420 Amin   = bplus
      Amax   = zero
      close1 = 0.5d+0
      close2 = 2.0d+0
      nclose = 0

      do j = 1, n
         Ac = one / Ascale(j)

         if (Amin .gt. Ac) then
             Amin =  Ac
             jmin =  j
         end if

         if (Amax .lt. Ac) then
             Amax =  Ac
             jmax =  j
         end if

         if (Ac .gt. close1  .and.  Ac .lt. close2) then
             nclose =  nclose + 1
*----        Ac     =  one
         end if

         Ascale(j) = Ac
      end do

*     Remember, column scales are upside down.

      Amax   = one / Amax
      Amin   = one / Amin
      if ( Prnt1 ) then
         call snPRNT( 31,
     &        '            Min scale'
     &     // '                       Max scale'
     &     // '      Between 0.5 and 2.0', iw, leniw )
         write(str, 1310)
     &        'Col', jmax, Amax, 'Col', jmin, Amin,
     &        nclose, nclose * 100.0d+0 / n
         call snPRNT( 21, str, iw, leniw )
      end if

*     ------------------------------------------------------------------
*     Phase 2.  Deal with empty rows and free rows.
*     Find the min and max row scales while we're at it.
*     ------------------------------------------------------------------
      Amin   = bplus
      Amax   = zero
      imin   = 0
      imax   = 0
      nclose = 0

      do  i = iStart, m
         j  = n + i

         if (hrtype(i) .eq. 2) then
            Ar = min( rmax(i), one )
            if (Ar .eq. zero) Ar = one
         else
            Ar = Ascale(j)
            if (Ar .eq. zero) Ar = one

            if (Amin .gt. Ar) then
                Amin  =   Ar
                imin  =   i
            end if
            if (Amax .lt. Ar) then
                Amax  =   Ar
                imax  =   i
            end if

            if (Ar .gt. close1  .and.  Ar .lt. close2) then
                nclose =  nclose + 1
*----           Ar     =  one
            end if
         end if

         Ascale(j) = Ar
      end do

      if (imin .eq. 0) then
          Amin   = zero
          Amax   = zero
      end if

      if ( Prnt1 ) then
         write(str, 1310)
     &        'Row', imin, Amin, 'Row', imax, Amax,
     &        nclose, nclose * 100.0d+0 / m
         call snPRNT( 21, str, iw, leniw )
      end if

*     ------------------------------------------------------------------
*     Phase 3.
*     Compute what is effectively the rhs for the constraints.
*     We set  rmax  =  ( A  -I )*x  for fixed columns and slacks,
*     including positive lower bounds and negative upper bounds.
*     ------------------------------------------------------------------
      if (.not. lonly) then
         call dload ( m, zero, rmax, 1 )

         do j = 1, nb
            bnd = zero
            b1  = bl(j)
            b2  = bu(j)

            if (b1  .eq. b2  ) bnd = b1
            if (b1  .gt. zero) bnd = b1
            if (b2  .lt. zero) bnd = b2

            if (bnd .ne. zero) then

               if (j .le. n) then
                  do l = locA(j), locA(j+1)-1
                     i = indA(l)
                     rmax(i) = rmax(i)  +  Acol(l) * bnd
                  end do
               else
*                 Free slacks never get here, no need to skip them.
                  i       = j - n
                  rmax(i) = rmax(i)  -  bnd
               end if
            end if
         end do

*        We don't want nonzeros in free rows to interfere.

         do i = 1, m
            if (hrtype(i) .eq. 2) rmax(i) = zero
         end do

*        Scale rmax = rmax / (row scales),  and use its norm sigma
*        to adjust both row and column scales.

         Ac     = dnormi( m, rmax, 1 )
         call dddiv ( m, Ascale(n+1), 1, rmax, 1 )
         sigma  = dnormi( m, rmax, 1 )
         if ( Prnt1 ) then
            write(str, 1400) Ac
            call snPRNT( 31, str, iw, leniw )
            write(str, 1410) sigma
            call snPRNT( 21, str, iw, leniw )
         end if
         sigma  = max( sigma, one )
         call dscal ( nb, sigma, Ascale, 1 )

*        Big scales might lead to excessive infeasibility when the
*        problem is unscaled.  If any are too big, scale them down.

         Amax = zero
         do j = 1, n
            Amax = max( Amax, Ascale(j) )
         end do

         do i = 1, m
            if (hrtype(i) .ne. 2) Amax = max( Amax, Ascale(n+i) )
         end do

         big    = 0.1d+0 / tolx
         sigma  = big    / Amax
         if (sigma .lt. one) then
            call dscal ( nb, sigma, Ascale, 1 )
            if ( Prnt1 ) then
               write(str, 1450) sigma
               call snPRNT( 21, str, iw, leniw )
            end if
         end if
      end if

      if (lprScl .gt. 0) then
         call snPRNT( 31, ' ', iw, leniw )
         call snPRNT( 21, ' Row scales  r(i)     '
     &     // ' a(i,j)  =  r(i) * scaled a(i,j) / c(j)', iw, leniw )
         call snPRNT( 21, ' ----------------', iw, leniw )

         do i = iStart, m
            write(str, 1500) i, Ascale(n+i)
            call snPRNT( 21, str, iw, leniw )
         end do

         call snPRNT( 31, ' ', iw, leniw )
         call snPRNT( 21, ' Column scales  c(j)     '
     &     // ' x(j)    =  c(j) * scaled x(j)', iw, leniw )
         call snPRNT( 21, ' -------------------', iw, leniw )

         do j = jStart, n
            write(str, 1500) j, Ascale(j)
            call snPRNT( 21, str, iw, leniw )
         end do
      end if

      return

 1200 format(' After', i3, 1p, e12.2, e12.2, 0p, f20.2)
 1310 format(1x, a, i7, 1p, e10.1, 12x, a, i7, e10.1, i17, 0p, f8.1)
 1400 format(' Norm of fixed columns and slacks', 1p, e20.1)
 1410 format(' (before and after row scaling)  ', 1p, e20.1)
 1450 format(' Scales are large --- reduced by ', 1p, e20.1)
 1500 format(i6, g16.5)

      end ! subroutine s2Scal

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2SclA
     &   ( Task, m, n, nb, iObj, InfBnd, sclObj,
     &     ne, nlocA, locA, indA, Acol,
     &     Ascale, bl, bu, pi, x )

      implicit
     &     none
      integer
     &     Task, m, n, nb, iObj, ne, nlocA, locA(nlocA), indA(ne)
      double precision
     &     sclObj, InfBnd,
     &     Acol(ne), Ascale(nb), bl(nb), bu(nb), x(nb), pi(m)

*     ==================================================================
*     s2SclA scales or unscales A, bl, bu, pi, x using row and column
*     scales Sr and Sc stored as  Ascale = ( Sc  Sr ).
*
*       A(scaled)  = inv(Sr)*A*Sc
*       x(scaled)  = inv(Sc)*x,    s(scaled) = inv(Sr)*s,
*       pi(scaled) =     Sr *pi.
*
*     15 Nov 1991: First version based on Minos 5.4 routine m2scla.
*     24 Mar 2000: Current version of s2SclA
*     ==================================================================
      integer
     &     i, j, l
      double precision
     &     bplus, cScale
*     ------------------------------------------------------------------
      double precision   one
      parameter         (one = 1.0d+0)
      integer            Scale,     UnScal
      parameter         (Scale = 0, UnScal = 1)
*     ------------------------------------------------------------------
      bplus  = 0.1d+0*InfBnd

      if (Task .eq. Scale) then
*        ---------------------------------------------------------------
*        Scale A, bl, bu, x and pi.
*        ---------------------------------------------------------------
         do j = 1, nb
            cScale = Ascale(j)
            if (j .le. n) then
               do l = locA(j), locA(j+1)-1
                  i = indA(l)
                  Acol(l) = Acol(l) * ( cScale / Ascale(n+i) )
               end do
            end if
            x(j) = x(j) / cScale
            if (bl(j) .gt. -bplus) bl(j) = bl(j) / cScale
            if (bu(j) .lt.  bplus) bu(j) = bu(j) / cScale
         end do

         call ddscl ( m, Ascale(n+1), 1, pi, 1 )
         if (iObj  .gt. 0) sclObj = Ascale(n+iObj)

      else if (Task .eq. UnScal) then
*        ---------------------------------------------------------------
*        Unscale everything.
*        ---------------------------------------------------------------
         do j = 1, nb
            cScale  = Ascale(j)
            if (j .le. n) then
               do l = locA(j), locA(j+1)-1
                  i = indA(l)
                  Acol(l) = Acol(l) * ( Ascale(n+i) / cScale )
               end do
            end if
            x(j) = x(j) * cScale
            if (bl(j) .gt. -bplus)  bl(j) = bl(j) * cScale
            if (bu(j) .lt.  bplus)  bu(j) = bu(j) * cScale
         end do

         call dddiv ( m, Ascale(n+1), 1, pi, 1 )
         sclObj = one
      end if

      end ! subroutine s2SclA

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2unpk
     &   ( jq, m, n, ne, ColNrm, nlocA, locA, indA, Acol, y )

      implicit
     &     none
      integer
     &     jq, m, n, ne, nlocA, locA(nlocA), indA(ne)
      double precision
     &     ColNrm, Acol(ne), y(m)

*     ==================================================================
*     s2unpk  expands the jq-th column of  ( A  -I )  into  y.
*
*     21 Jun 2004: Norm of unpacked column computed.
*     21 Jun 2004: Current version of s2unpk
*     ==================================================================
      integer
     &     i, islack, l
*     ------------------------------------------------------------------
      double precision   zero,          one
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )
*     ------------------------------------------------------------------
      call dload ( m, zero, y, 1 )

      ColNrm = one

      if (jq .le. n) then
         do l      = locA(jq), locA(jq+1)-1
            i      = indA(l)
            y(i)   = Acol(l)
            ColNrm = max( ColNrm, abs(y(i)) )
         end do
      else
         islack    =   jq - n
         y(islack) = - one
      end if

      end ! subroutine s2unpk

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2vmax
     &   ( n, nnCon, maxvi, vimax, bl, bu, Fx )

      implicit
     &     none
      integer
     &     maxvi, n, nnCon
      double precision
     &     vimax, bl(n+nnCon), bu(n+nnCon), Fx(nnCon)

*     ==================================================================
*     s2vmax  finds the largest nonlinear constraint violation.
*     maxvi   points to the biggest violation in Fx.
*     vimax   is the biggest violation in Fx.
*
*     13 Apr 1999: First version of s2viol.
*     21 Oct 2000: Current version.
*     ==================================================================
      integer
     &     i, j
      double precision
     &     slacki, viol
*     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero = 0.0d+0)
*     ------------------------------------------------------------------
*     See how much  Fx  violates the bounds on the nonlinear slacks.

      vimax  = zero
      maxvi  = 1

      do i = 1, nnCon
         j      = n + i
         slacki = Fx(i)
         viol   = max( zero, bl(j) - slacki, slacki - bu(j) )
         if (vimax .lt. viol) then
             vimax = viol
             maxvi = i
         end if
      end do

      end ! subroutine s2vmax

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2xmat
     &   ( matfil, n, nb, ne, nlocA, locA, indA, Acol, hs )

      implicit
     &     none
      integer
     &     matfil, n, nb, ne, nlocA, locA(nlocA), indA(ne), hs(nb)
      double precision
     &     Acol(ne)

*     ==================================================================
*     s2xmat  outputs A or B or (B S) to file matfil.
*     A triple (i, j, aij) is output for each nonzero entry,
*     intended for input to Matlab.
*
*     matfil    Matrix    hs(j)    j
*       91        A       >= 0    1:n
*       92        B       >= 3    1:nb
*       93       (B S)    >= 2    1:nb
*
*     10 Oct 2000: First version of s2xmat.
*     26 Jul 2003: Allow for hs(*) having "internal" values.
*     ==================================================================
      integer
     &     i, j, js, k, l, ls, ncol
      double precision
     &     aik
*     ------------------------------------------------------------------
      double precision   zero,           one
      parameter        ( zero = 0.0d+0,  one = 1.0d+0 )
*     ------------------------------------------------------------------
      if      (matfil .eq. 91) then
         ls     = 0
         ncol   = n
      else if (matfil .eq. 92) then
         ls     = 3
         ncol   = nb
      else if (matfil .eq. 93) then
         ls     = 2
         ncol   = nb
      else
         return
      end if

      ! Output certain columns of A.
      ! Treat internal values hs(j) = -1 or 4 as 0 (nonbasic).

      k   = 0

      do j  = 1, n
         js = hs(j)
         if (js .lt. 0  .or.  js .eq. 4) js = 0
         if (js .ge. ls) then
            k    = k + 1
            do l = locA(j), locA(j+1)-1
               i      = indA(l)
               aik    = Acol(l)
               if (aik .ne. zero) then
                  write(matfil, 1000) i, k, aik
               end if
            end do
         end if
      end do

      ! Output certain columns of -I.

      aik  = - one

      do j  = n+1, ncol
         js = hs(j)
         if (js .lt. 0  .or.  js .eq. 4) js = 0
         if (js .ge. ls) then
            k    = k + 1
            i    = j - n
            write(matfil, 1000) i, k, aik
         end if
      end do

      return

 1000 format( 1p, i10, i10, e24.14 )

      end ! subroutine s2xmat
