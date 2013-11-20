*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  sn35mps.f
*
*     MPSinp   MPSout
*     s3dflt   s3getp   s3hash   s3inpt   s3MapM   s3mps
*     s3mpsa   s3mpsb   s3mpsc   s3read
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine MPSinp
     &   ( iMPS, maxm, maxn, maxne,
     &     nnCon, nnJac, nnObj,
     &     m, n, ne,
     &     iObj, ObjAdd, PrbNms,
     &     Acol, indA, locA, bl, bu, Names,
     &     hs, x, pi,
     &     iExit, mincw, miniw, minrw, nS,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iExit, iMPS, maxm, maxn, maxne, nnCon, nnJac, nnObj, m,
     &      mincw, miniw, minrw, n, ne, iObj, nS, lencw, leniw, lenrw,
     &     indA(maxne), hs(maxm+maxn), locA(maxn+1), iw(leniw)
      double precision
     &     ObjAdd, Acol(maxne), bl(maxm+maxn), bu(maxm+maxn),
     &     x(maxn+maxm), pi(maxm), rw(lenrw)
      character
     &     PrbNms(5)*8, Names(maxm+maxn)*8, cw(lencw)*8

*     ------------------------------------------------------------------
*     MPSinp  inputs constraint data for a linear or nonlinear program
*     in MPS format, consisting of NAME, ROWS, COLUMNS, RHS, RANGES and
*     BOUNDS sections in that order.  The RANGES and BOUNDS sections are
*     optional.
*
*     ------------------------------------------------------------------
*     NOTE: Before calling MPSinp, your calling program MUST call the
*     initialization routine using the call:
*     call snInit( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )
*     This sets the default values of the optional parameters. You can
*     also alter the default values of iPrint and iSumm before MPSinp
*     is used.  iPrint = 0, etc, is OK.
*     ------------------------------------------------------------------
*
*     In the LP case, MPS format defines a set of constraints of the
*     form
*              l <= x <= u,      b1 <=  Ax  <= b2,
*     where l and u are specified by the BOUNDS section, and b1 and b2
*     are defined somewhat indirectly by the ROWS, RHS and RANGES
*     sections.  snmps converts these constraints into the equivalent
*     form
*              Ax - s = 0,       bl <= ( x ) <= bu,
*                                      ( s )
*     where s is a set of slack variables.  This is the way SNOPT deals
*     with the data.  The first n components of bl and bu are the same
*     as l and u.  The last m components are b1 and b2.
*
*     MPS format gives 8-character names to the rows and columns of A.
*     One of the rows of A may be regarded as a linear objective row.
*     This will be row iObj, where iObj = 0 means there is no such row.
*
*     The data defines a linear program if nnCon = nnJac = nnObj = 0.
*     The nonlinear case is the same except for a few details.
*     1. If nnCon = nnJac = 0 but nnObj > 0, the first nnObj columns
*        are associated with a nonlinear objective function.
*     2. If nnCon > 0, then nnJac > 0 and nnObj may be zero or positive.
*        The first nnCon rows and the first nnJac columns are associated
*        with a set of nonlinear constraints.
*     3. Let nnL = max( nnJac, nnObj ).  The first nnL columns
*        correspond to "nonlinear variables".
*     4. If an objective row is specified (iObj > 0), then it must be
*        such that iObj > nnCon.
*     5. "Small" elements (below the Aij tolerance) are ignored only if
*        they lie outside the nnCon by nnJac Jacobian, i.e. outside
*        the top-left corner of A.
*     6. No warning is given if some of the first nnL columns are empty.
*
*
*     ON ENTRY
*     ========
*     iMPS   is the unit containing the MPS file.  On some systems, it
*            may be necessary to open file iMPS before calling snmps.
*
*     maxm   is an overestimate of the number of rows in the ROWS
*            section of the MPS file.
*
*     maxn   is an overestimate of the number of columns in the COLUMNS
*            section of the MPS file.
*
*     maxne  is an overestimate of the number of elements (matrix
*            coefficients) in the COLUMNS section.
*
*     nnCon  is the no. of nonlinear constraints in the problem.
*            These must be the FIRST rows in the ROWS section.
*
*     nnJac  is the no. of nonlinear Jacobian variables in the problem.
*            These must be the FIRST columns in the COLUMNS section.
*
*     nnObj  is the no. of nonlinear objective variables in the problem.
*            These must be the FIRST columns in the COLUMNS section,
*            overlapping where necessary with the Jacobian variables.
*
*     PrbNms is an array of five 8-character names.
*            PrbNms(1) need not be specified... it will be changed to
*            the name on the NAME card of the MPS file.
*            PrbNms(2) is the name of the objective row to be selected
*            from the ROWS section, or blank if snmps should select
*            the first type N row encountered.
*            Similarly,
*            PrbNms(3), PrbNms(4) and PrbNms(5) are the names of the
*            RHS, RANGES and BOUNDS to be selected from the
*            RHS, RANGES and BOUNDS sections respectively, or blank
*            if snmps should select the first ones encountered.
*
*     iw(*)  is a workspace array of length leniw.  It is needed to
*            hold the row-name hash table and a few other things.
*
*     leniw  is the length of iw(*).  It should be at least 4*maxm.
*
*     rw(*)  is a workspace array of length lenrw.  It is needed to
*            hold the row-name hash table and a few other things.
*
*     lenrw  is the length of rw(*).  It should be at least 4*maxm.
*
*
*     ON EXIT
*     =======
*     m       is the number of rows in the ROWS section.
*
*     n       is the number of columns in the COLUMNS section.
*
*     ne      is the number of matrix coefficients in COLUMNS section.
*
*     iObj    is the row number of the specified objective row,
*             or zero if no such row was found.
*
*     ObjAdd  is a real constant extracted from row iObj of the RHS.
*             It is zero if the RHS contained no objective entry.
*             SNOPT adds ObjAdd to the objective function.
*
*     PrbNms(1)-PrbNms(5) contain the names of the
*             Problem, Objective row, RHS, RANGES and BOUNDS
*             respectively.
*
*     Acol(*) contains the ne entries for each column of the matrix
*             specified in the COLUMNS section.
*
*     indA(*) contains the corresponding row indices.
*
*     locA(j) (j = 1:n) points to the beginning of column j
*             in the parallel arrays Acol(*), indA(*).
*     locA(n+1) = ne+1.
*
*     bl(*)  contains n+m lower bounds for the columns and slacks.
*            If there is no lower bound on x(j), then bl(j) = - 1.0d+20.
*
*     bu(*)  contains n+m lower bounds for the columns and slacks.
*            If there is no upper bound on x(j), then bu(j) = + 1.0d+20.
*
*     Names(*) contains n+m column and row names in character*8 format.
*            The j-th column name is stored in Names(j).
*            The i-th row    name is stored in Names(k),
*            where k = n + i.
*
*     hs(*)  contains an initial state for each column and slack.
*
*     x(*)   contains an initial value for each column and slack.
*
*            If there is no INITIAL bounds set,
*               x(j) = 0 if that value lies between bl(j) and bu(j),
*                     = the bound closest to zero otherwise,
*               hs(j) = 0 if x(j) < bu(j),
*                     = 1 if x(j) = bu(j).
*
*            If there is an INITIAL bounds set, x(j) and hs(j) are
*            set as follows.  Suppose the j-th variable has the name Xj,
*            and suppose any numerical value specified happens to be 3.
*                                                   x(j)    hs(j)
*             FR INITIAL   Xj         3.0           3.0       -1
*             FX INITIAL   Xj         3.0           3.0        2
*             LO INITIAL   Xj                       bl(j)      4
*             UP INITIAL   Xj                       bu(j)      5
*             MI INITIAL   Xj         3.0           3.0        4
*             PL INITIAL   Xj         3.0           3.0        5
*
*     pi(*)  contains a vector defined by a special RHS called LAGRANGE.
*            If the MPS file contains no such RHS, pi(i) = 0.0, i=1:m.
*
*     iExit  =   0  if no fatal errors were encountered.
*            =  83  if there is not enough integer workspace.
*            = 111  if no MPS file was specified.
*            = 112  if maxm, maxn or maxne were too small.
*            = 113  if the ROWS or COLUMNS sections were empty
*                   or iObj > 0 but iObj <= nnCon,
*
*     nS     is the no. of FX INITIAL entries in the INITIAL bounds set.
*
*
*     09 Jul 1997: Original version, derived from mimps.
*     27 Oct 2003: Current version of MPSinp.
*     ==================================================================
      character
     &     Solver*6, key*4, str*80, str2*80
      integer
     &     mProb, mObj, mRhs, mRng, mBnd, maxcw, maxiw, maxrw, negCon,
     &     nextcw, nextiw, nextrw, nnL, lenh, lhrtyp, lkBS, lkeynm
*     ------------------------------------------------------------------
      parameter         (mProb     =  51) ! Problem name
      parameter         (mObj      =  52) ! Objective name
      parameter         (mRhs      =  53) ! rhs name
      parameter         (mRng      =  54) ! range name
      parameter         (mBnd      =  55) ! bounds name
*     ------------------------------------------------------------------
      iExit  = 0
      Solver = 'MPSinp'

*     ------------------------------------------------------------------
*     Check memory limits and fetch the workspace starting positions.
*     ------------------------------------------------------------------
      call s2Mem0
     &   ( iExit, Solver, lencw, leniw, lenrw, iw,
     &     mincw, miniw, minrw, maxcw, maxiw, maxrw,
     &     nextcw, nextiw, nextrw )
      if (iExit .ne. 0) go to 999

*     Set undefined MPS options to default values.
*     Quit if no MPS file is specified.

      call s3dflt
     &   ( cw, lencw, iw, leniw, rw, lenrw )

      iw( 21)   = nnJac ! # nonlinear Jacobian variables
      iw( 22)   = nnObj ! # variables in gObj
      iw( 23)   = nnCon ! # of nonlinear constraints

      cw(mProb) = PrbNms(1)
      cw(mObj ) = PrbNms(2)
      cw(mRhs ) = PrbNms(3)
      cw(mRng ) = PrbNms(4)
      cw(mBnd ) = PrbNms(5)

      iw(123)   = iMPS
      if (iMPS .le. 0) then
         iExit = 111
         go to 800
      end if

      nnL    = max( nnJac, nnObj )

*     ------------------------------------------------------------------
*     Allocate workspace for the MPS input routines.
*     s3getp finds a prime number for the length of the row hash table.
*     ------------------------------------------------------------------
      call s3getp( maxm, lenh )

      lhrtyp = nextiw
      lkBS   = lhrtyp  + 1 + maxm
      lkeynm = lkBS    + 1 + maxm
      nextiw = lkeynm  + lenh
      miniw  = nextiw  - 1

      if (miniw .gt. maxiw) then
         write(str, 9830) miniw
         call snPRNT( 13, str, iw, leniw )
         iExit = 83
         go to 800
      end if

      call s3mps
     &   ( iExit, key,
     &     maxm, maxn, maxne, lenh, nnL, nnCon, nnJac,
     &     m, n, ne, negCon, nS,
     &     iObj, ObjAdd,
     &     locA, indA, Acol, bl, bu, Names,
     &     hs, iw(lhrtyp), iw(lkBS), iw(lkeynm),
     &     x, pi,
     &     cw, lencw, iw, leniw, rw, lenrw )

  800 if (iExit .eq. 0) then
         iExit = 103            ! MPS file read successfully
         PrbNms(1) = cw(mProb)
         PrbNms(2) = cw(mObj )
         PrbNms(3) = cw(mRhs )
         PrbNms(4) = cw(mRng )
         PrbNms(5) = cw(mBnd )
      end if

      call snWRAP( iExit, Solver, str, str2, iw, leniw )

  999 return

 9830 format(' Total integer   workspace should be significantly',
     &       ' more than', i8)

      end ! subroutine MPSinp

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine MPSout
     &   ( iMPS, m, n, ne, nName, PrbNms,
     &     Acol, indA, locA, bl, bu, Names )

      implicit
     &     none
      integer
     &     iMPS, nName, m, n, ne, indA(ne), locA(n+1)
      double precision
     &     Acol(ne), bl(n+m), bu(n+m)
      character
     &     PrbNms(5)*8, Names(nName)*8

*     ==================================================================
*     mpsout  outputs an MPS file to file number iMPS.
*     All parameters are the same as for snOptB in SNOPT 6 onwards.
*     They are all input parameters.
*
*     16 Aug 1998: Original version, derived from MINOS 5.4 mpsout.
*     27 Oct 2003: Current version.
*     ==================================================================
      logical
     &     value
      character
     &     ic*8, id*8, rowtyp*1, bndtyp*4
      integer
     &     i, j, k, kform, nb
      double precision
     &     ai, b1, b2, bnd, rng
*     ------------------------------------------------------------------
      double precision   zero,             one
      parameter        ( zero   = 0.0d+0,  one    =  1.0d+0  )
      double precision   bplus,            bminus
      parameter        ( bplus  = 1.0d+19, bminus = -1.0d+19 )
*     ------------------------------------------------------------------
      character          form(6)*29
      data               form(1) /'(4x,  a8, 2x,  a8, f8.0     )'/
      data               form(2) /'(4x,  a8, 2x,  a8, 1p, e14.5)'/
      data               form(3) /'(4x,  a8, 2x,  a8, f8.0     )'/
      data               form(4) /'(4x,  a8, 2x,  a8, 1p, e14.5)'/
      data               form(5) /'(a4,  a8, 2x,  a8, f8.0     )'/
      data               form(6) /'(a4,  a8, 2x,  a8, 1p, e14.5)'/
*     ------------------------------------------------------------------
      if (iMPS .le. 0) return

      nb = n + m

*     ------------------------------------------------------------------
*     ROWS section.
*     Note: b1 and b2 are bounds on ROWS, not slacks.
*           The objective row gets its name from name1(*), name2(*).
*           name(2) is ignored.
*     ------------------------------------------------------------------
      write(iMPS, '(a, 10x, a)') 'NAME', PrbNms(1)
      write(iMPS, '(a)'        ) 'ROWS'

      do i = 1, m
         j  =   n + i
         b1 = - bu(j)
         b2 = - bl(j)
         if (b1 .eq. b2) then
            rowtyp = 'E'
         else if (b1 .gt. bminus) then
            rowtyp = 'G'
         else if (b2 .lt. bplus ) then
            rowtyp = 'L'
         else
            rowtyp = 'N'
         end if

         call s4id  ( j, n, nb, nName, Names, id )
         write(iMPS, '(1x, a1, 2x, a8)') rowtyp, id
      end do

*     ------------------------------------------------------------------
*     COLUMNS section.
*     Note: Objective entries get their name from name1(*), name2(*).
*           name(2) is ignored.
*     ------------------------------------------------------------------
      write(iMPS, '(a)') 'COLUMNS'

      do j = 1, n
         call s4id  ( j, n, nb, nName, Names, ic )

         do k = locA(j), locA(j+1)-1
            i = indA(k)
            call s4id  ( n+i, n, nb, nName, Names, id )

            ai    = Acol(k)
            kform = 2
            if (ai .eq. zero  .or.  abs(ai) .eq. one) kform = 1
            write(iMPS, form(kform)) ic, id, ai
         end do
      end do

*     ------------------------------------------------------------------
*     RHS section.
*     Note: b1 and b2 are bounds on ROWS, not slacks.
*     ------------------------------------------------------------------
      write(iMPS, '(a)') 'RHS'

      do i = 1, m
         j   =   n + i
         b1  = - bu(j)
         b2  = - bl(j)
         bnd =   zero
         if (b1 .eq. b2) then
            bnd = b1
         else if (b1 .gt. bminus) then
            bnd = b1
         else if (b2 .lt. bplus ) then
            bnd = b2
         end if

         if (bnd .ne. zero) then
            call s4id  ( j, n, nb, nName, Names, id )
            kform = 4
            if (abs(bnd) .eq. one) kform = 3
            write(iMPS, form(kform)) PrbNms(3), id, bnd
         end if
      end do

*     ------------------------------------------------------------------
*     RANGES section.
*     ------------------------------------------------------------------
      write(iMPS, '(a)') 'RANGES'

      do i   = 1, m
         j   =   n + i
         b1  = - bu(j)
         b2  = - bl(j)
         if (b1 .lt. b2  .and.  b1 .gt. bminus .and. b2 .lt. bplus) then
            rng   = b2 - b1
            call s4id  ( j, n, nb, nName, Names, id )
            kform = 4
            if (abs(rng) .eq. one) kform = 3
            write(iMPS, form(kform)) PrbNms(4), id, rng
         end if
      end do

*     ------------------------------------------------------------------
*     BOUNDS section.
*     ------------------------------------------------------------------
      write(iMPS, '(a)') 'BOUNDS'

      do j  = 1, n
         b1 = bl(j)
         b2 = bu(j)
         bndtyp = '    '
         value  = .false.

*        Output lower bound, except for vanilla variables.

         if (b1 .eq. b2) then
            bndtyp = ' FX '
            value  = .true.
         else if (b1 .gt. bminus) then
            if (b1 .ne. zero) then
               bndtyp = ' LO '
               value  = .true.
            end if
         else if (b2 .lt. bplus) then
               bndtyp = ' MI '
         else
            bndtyp = ' FR '
         end if

         if (bndtyp .ne. '    ') then
            call s4id  ( j, n, nb, nName, Names, id )
            if (value) then
               kform = 6
               if (b1 .eq. zero  .or.  abs(b1) .eq. one) kform = 5
               write(iMPS, form(kform)) bndtyp, PrbNms(5), id, b1
            else
               write(iMPS, form(kform)) bndtyp, PrbNms(5), id
            end if
         end if

*        Output second bound if necessary.

         bndtyp = '    '
         value  = .false.

         if (b1 .eq. b2) then
*           do nothing
         else if (b1 .gt. bminus) then
            if (b2 .lt. bplus) then
               bndtyp = ' UP '
               value  = .true.
            end if
         else if (b2 .lt. bplus) then
            if (b2 .ne. zero) then
               bndtyp = ' UP '
               value  = .true.
            end if
         end if

         if (bndtyp .ne. '    ') then
            call s4id  ( j, n, nb, nName, Names, id )
            if (value) then
               kform = 6
               if (b2 .eq. zero  .or.  abs(b2) .eq. one) kform = 5
               write(iMPS, form(kform)) bndtyp, PrbNms(5), id, b2
            else
               write(iMPS, form(kform)) bndtyp, PrbNms(5), id
            end if
         end if
      end do

      write(iMPS, '(a)') 'ENDATA'

      end ! subroutine MPSout

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3dflt
     &   ( cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

*     ==================================================================
*     s3dflt  sets undefined MPS options to their default values and
*     prints them.
*
*     15 Feb 1998: First version.
*     05 Feb 2001: Current version of s3dflt.
*     ==================================================================
      character
     &     str*132, cdummy*8
      integer
     &     Aijtol, bStrc1, bStrc2, lvlTim, iDummy, lprPrm,
     &     minmax, lDenJ, MPSerr, mLst, nProb, iMPS, nnCon, maxm, maxn,
     &     maxne, mProb, mObj, mRhs, mRng, mBnd, infBnd, rdummy
      double precision
     &     t
*     ------------------------------------------------------------------
      parameter         (infBnd    =  70) ! def. of an infinite bound
      parameter         (rdummy    =  69) ! def. of an 'unset' value

      parameter         (Aijtol    =  95) ! zero Aij tolerance.
      parameter         (bStrc1    =  96) ! default lower bound on x
      parameter         (bStrc2    =  97) ! default upper bound on x

      parameter         (lprPrm    =  81) ! > 0    => parms are printed
      parameter         (minmax    =  87) ! 1, 0, -1  => MIN, FP, MAX
      parameter         (lDenJ     = 105) ! 1(2)    => dense(sparse) J
      parameter         (MPSerr    = 106) ! maximum # errors in MPS data
      parameter         (mLst      = 107) ! maximum # lines  of MPS data
      parameter         (nProb     = 108) ! problem number
      parameter         (iMPS      = 123) ! MPS file

      parameter         (nnCon     =  23) ! # of nonlinear constraints
      parameter         (maxm      = 133) ! Estimate of the # of rows
      parameter         (maxn      = 134) ! Estimate of the # of columns
      parameter         (maxne     = 135) ! Estimate of the # of elements
      parameter         (lvlTim    = 182) ! Timing level

      parameter         (mProb     =  51) ! Problem name
      parameter         (mObj      =  52) ! Objective name
      parameter         (mRhs      =  53) ! rhs name
      parameter         (mRng      =  54) ! range name
      parameter         (mBnd      =  55) ! bounds name

      parameter         (idummy =  -11111)
      parameter         (cdummy = '-1111111')
      double precision   zero
      parameter         (zero   =  0.0d+0)
*     ------------------------------------------------------------------
      character          id(2)*4
      character          cBlank*8
      data               cblank /'        '/
      data               id     /' Den', 'Spar'/
*     ------------------------------------------------------------------
      if (cw(mProb)  .eq. cdummy     ) cw(mProb)  = cBlank
      if (cw(mObj )  .eq. cdummy     ) cw(mObj )  = cBlank
      if (cw(mRhs )  .eq. cdummy     ) cw(mRhs )  = cBlank
      if (cw(mRng )  .eq. cdummy     ) cw(mRng )  = cBlank
      if (cw(mBnd )  .eq. cdummy     ) cw(mBnd )  = cBlank

      if (iw(maxm)   .le. 0          ) iw(maxm)   = 100
      if (iw(maxn)   .le. 0          ) iw(maxn)   = 3*iw(maxm)
      if (iw(maxne)  .le. 0          ) iw(maxne)  = 5*iw(maxn)

      if (iw(lprPrm) .lt. 0          ) iw(lprPrm) =  1
      if (iw(iMPS  ) .eq. idummy     ) iw(iMPS  ) =  0

      if (iw(lvlTim) .lt. 0          ) iw(lvlTim) =  3
      if (iw(nProb ) .lt. 0          ) iw(nProb ) =  0
      if (iw(lDenJ ) .lt. 0          ) iw(lDenJ ) =  1
      if (iw(MPSerr) .lt. 0          ) iw(MPSerr) = 10
      if (iw(mLst  ) .lt. 0          ) iw(mLst  ) =  0
      if (iw(minmax) .eq. idummy     ) iw(minmax) =  1

      if (rw(infBnd) .lt. zero       ) rw(infBnd) = 1.0d+20
      if (rw(Aijtol) .lt. zero       ) rw(Aijtol) = 1.0d-10
      if (rw(bStrc1) .eq. rw(rdummy) ) rw(bStrc1) = zero
      if (rw(bStrc2) .eq. rw(rdummy) ) rw(bStrc2) = rw(infBnd)

      if (rw(Bstrc1) .gt. rw(Bstrc2)) then
         t          =   rw(Bstrc1)
         rw(Bstrc1) =   rw(Bstrc2)
         rw(Bstrc2) =   t
      end if

*     ------------------------------------------------------------------
*     Print parameters unless SUPPRESS PARAMETERS was specified.
*     ------------------------------------------------------------------
      if (iw(lprPrm) .gt. 0) then
         call s1page( 1, iw, leniw )
         call snPRNT( 1, ' MPS Input Data', iw, leniw )
         call snPRNT( 1, ' ==============', iw, leniw )
         write(str, 2000) iw(iMPS)
         call snPRNT( 1, str, iw, leniw )
         write(str, 2010) iw(maxm)  , iw(nProb) , rw(Bstrc1)
         call snPRNT( 1, str, iw, leniw )
         write(str, 2020) iw(maxn)  , iw(mLst)  , rw(Bstrc2)
         call snPRNT( 1, str, iw, leniw )
         write(str, 2030) iw(maxne) , iw(MPSerr), rw(Aijtol)
         call snPRNT( 1, str, iw, leniw )
         if (iw(nnCon) .gt. 0) then
            write(str, 2100) id(iw(lDenJ))
            call snPRNT( 1, str, iw, leniw )
         end if
      end if

      return

 2000 format(
     &     ' MPS file ..............', i10)
 2010 format(
     &     ' Row limit..............', i10, 6x,
     &     ' Problem Number.........', i10, 6x,
     &     ' Lower bound default....', 1p,  e10.2)
 2020 format(
     &     ' Column limit...........', i10, 6x,
     &     ' List limit.............', i10, 6x,
     &     ' Upper bound default....', 1p,  e10.2)
 2030 format(
     &     ' Elements limit ........', i10, 6x,
     &     ' Error message limit....', i10, 6x,
     &     ' Aij tolerance..........', 1p,  e10.2)
 2100 format(
     &     ' Jacobian...............', 4x, a4, 'se')

      end ! subroutine s3dflt

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3getp
     &   ( maxm, lenh )

      implicit
     &     none
      integer
     &     maxm, lenh
*     ------------------------------------------------------------------
*     s3getp finds a prime number lenh suitably larger than maxm.
*     It is used as the length of the hash table for the MPS row names.
*     ------------------------------------------------------------------
      integer
     &     i, k
*     ------------------------------------------------------------------
      lenh = maxm*2
      lenh = max( lenh, 100 )
      lenh = (lenh/2)*2 - 1
      k    = lenh/20 + 6

  100 k    = k + 1
      lenh = lenh + 2
      do i = 3, k, 2
         if (mod(lenh,i) .eq. 0) go to 100
      end do

      end ! subroutine s3getp

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3hash
     &   ( len, nen, ncoll,
     &     key, mode, keytab,
     &     Names, ka, found )

      logical
     &     found
      integer
     &     len, nen, ncoll, mode, ka, keytab(len)
      character
     &     key*8, Names(nen)*8

*     ==================================================================
*     s3hash looks up and/or inserts keys in a hash table.
*     Reference:  R.P. Brent, CACM 16,2 (Feb 1973), pp. 105-109.
*     This version is simplified for the case where no entries are
*     deleted.
*
*     keytab is used as an index into a consecutive list of unique
*     identifiers Names(*), to find if "key" is in the list.
*
*     On entry,
*     len     is the (fixed) length of the hash table.
*             It must be significantly more than nen.
*     nen     is the number of unique identifiers.
*     ncoll   is the number of collisions so far.
*     key     is the identifier to be looked up or inserted.
*     mode    = 1 if key is just being looked up,
*             = 2 if key should be inserted into the hash table.
*     keytab  is the integer hash table.
*     Names   is the list of identifiers.
*
*     On exit,
*     ncoll   is (new) number of collisions so far.
*     ka      satisfies Names(ka) = key if key is in the hash table.
*     found   is true if key was already in the hash table on entry.
*
*     Approx 1975: First f66 version written for MINOS, using integer
*                  names.
*     16 Sep 1997: key and Names(*) are now character*8.
*                  Statement function "hash" introduced.  It should be
*                  the only thing different between f77 and f90.
*     ==================================================================
*     The following statement function "hash" returns an integer from
*     its character*8 argument.  f77 doesn't permit great efficiency
*     (or we would have done it this way long ago!).  f90 should be
*     better with the transfer function.

      integer            hash
      character          p*8

*-->  f77
*     16 Sep 1997: Try looking at just a few important characters.
*     02 Oct 1997: Adding ichar of p(1,2,5,8) gave over a million
*                  collisions on pilots.mps.
*                  Try left shifts and more chars.

      hash(p) = 16*(16*(16*(16*(16*ichar( p(1:1) )
     &                           + ichar( p(2:2) ))
     &                           + ichar( p(3:3) ))
     &                           + ichar( p(6:6) ))
     &                           + ichar( p(7:7) ))
     &                           + ichar( p(8:8) )

*-->  f90
*     16 Sep 1997: Simulate the function used in MINOS (with Name(*)
*                  treated as two integer arrays).  It seems reasonably
*                  efficient in terms of the number of collisions.

*     hash(p) =   transfer( p(1:4), mode )
*    &          - transfer( p(5:8), mode )
*     ==================================================================

      len2  = len - 2
      ic    = -1

*     Compute address of first probe (ir) and increment (iq).

      ikey  = hash( key )
      iq    = mod(ikey, len2) + 1
      ir    = mod(ikey, len)  + 1
      ka    = ir

*     Look in the table.

   20 kt    = keytab(ka)

*     Check for an empty space or a match.

      if (kt  .eq. 0       ) go to 30
      if (key .eq. Names(kt)) go to 60
      ic    = ic + 1
      ncoll = ncoll + 1

*     Compute address of next probe.
      ka    = ka + iq
      if (ka .gt. len) ka = ka - len

*     See if whole table has been searched.
      if (ka .ne. ir ) go to 20

*     The key is not in the table.
   30 found = .false.

*     Return with ka = 0 unless an entry has to be made.

      if (mode .eq. 2  .and.  ic .le. len2) go to 70
      ka    = 0
      return

   60 found = .true.
      return

*     Look for the best way to make an entry.

   70 if (ic .le. 0) return
      ia    = ka
      is    = 0

*     Compute the maximum length to search along current chain.
   80 ix    = ic - is
      kt    = keytab(ir)

*     Compute increment jq for current chain.

      ikey  = hash( Names(kt) )
      jq    = mod(ikey, len2) + 1
      jr    = ir

*     Look along the chain.
   90 jr    = jr + jq
      if (jr .gt. len) jr = jr - len

*     Check for a hole.

      if (keytab(jr) .eq. 0) go to 100
      ix    = ix - 1
      if (ix .gt. 0) go to 90
      go to 110

*     Save location of hole.

  100 ia    = jr
      ka    = ir
      ic    = ic - ix

*     Move down to the next chain.

  110 is    = is + 1
      ir    = ir + iq
      if (ir .gt. len) ir = ir - len

*     Go back if a better hole might still be found.
      if (ic .gt. is ) go to 80

*     If necessary move an old entry.
      if (ia .ne. ka ) keytab(ia) = keytab(ka)

      end ! subroutine s3hash

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3inpt
     &   ( iExit, maxm, maxn, maxne, nnCon, nnJac, nnObj,
     &     m, n, ne,
     &     iObj, ObjAdd,
     &     nextcw, nextiw, nextrw,
     &     mincw, miniw, minrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iExit, iObj, lencw, leniw, lenrw, maxm, maxn, maxne,
     &     mincw, miniw, minrw, m, n, ne, nextcw, nextiw, nextrw,
     &     nnCon, nnJac, nnObj, iw(leniw)
      double precision
     &     ObjAdd, rw(lenrw)
      character
     &     cw(lencw)*8

*     ==================================================================
*     s3inpt  inputs constraint data in MPS format, and sets up
*     various quantities as follows:
*
*     ObjAdd     (output) is minus the coefficient in row iObj of the
*                RHS section (zero by default).  SNOPT adds it to the
*                objective function.
*
*     m, n, ne   are the number of rows, columns and elements in A.
*
*     iObj       is the row number for the linear objective (if any).
*                It must come after any nonlinear rows.
*                iObj = 0 if there is no linear objective.
*
*     locA, indA, Acol
*                is the matrix A,
*                stored in locations llocA, lindA, lAcol
*     bl, bu     are the bounds,  stored at locations lbl, lbu.
*     hs, x      are states and values,   stored at   lhs, lx.
*
*     hs(j)      is set to  0, 1  to indicate a plausible initial state
*                (at lo or up bnd) for each variable j  (j = 1 to n+m).
*                If crash is to be used, i.e., crash option gt 0 and
*                if no basis file will be supplied, the initial bounds
*                set may initialize hs(j) as follows to assist crash:
*
*     -1      if column or row j is likely to be in the optimal basis,
*      4      if column j is likely to be nonbasic at its lower bound,
*      5      if column j is likely to be nonbasic at its upper bound,
*      2      if column or row j should initially be superbasic,
*      0 or 1 otherwise.
*
*     x(j)       is a corresponding set of initial values.
*                Safeguards are applied later by SNOPT, so the
*                values of hs and x are not desperately critical.
*
*     18 Nov 1991: First version based on Minos routine m3inpt.
*     09 Nov 2000: miniw, minrw, mincw set elsewhere.
*     27 Oct 2003: Current version of s3inpt.
*     ==================================================================
      character
     &     str*80, key*4, lENDA*4
      integer
     &     iDummy, iMPS, iStdi, iSpecs, klocA, khs, kbl, kbu, kx, kpi,
     &     kNames, lAcol, lindA, llocA, lbl, lbu, lx, lpi, lhrtyp, lhs,
     &     lkBS, lkeynm, lNames, lDenJ, lenh, newcw, newiw, newrw,
     &     negCon, nnL, nS
*     ------------------------------------------------------------------
      parameter         (negCon    =  20)
      parameter         (lENDA     = 'ENDA')
*     ------------------------------------------------------------------
      iStdi  = iw(  9) ! Standard Input
      iSpecs = iw( 11) ! Specs (options) file

      lDenJ  = iw(105) ! 1(2)    => dense(sparse) Jacobian
      iMPS   = iw(123) ! MPS file

      nnL    = max( nnJac, nnObj )

      newcw = nextcw
      newiw = nextiw
      newrw = nextrw

*     Start.  We may come back here to try again with more workspace.
*     key retains the first 4 characters of the NAME, ROWS, COLUMNS
*     RHS, RANGES and BOUNDS cards.
*     s3getp finds a prime number for the length of the row hash table.

  100 iExit = 0

      call s3getp( maxm, lenh )

*     Allocate temporary addresses for the problem data
*     based on maxm, maxn and maxne.

      call s3MapM
     &   ( iExit, maxm, maxn, maxne, lenh,
     &     nextcw, nextiw, nextrw,
     &     mincw, miniw, minrw,
     &     lencw, leniw, lenrw, iw )
      if (iExit .ne. 0) go to 700

*     Get the data from an MPS file.
*     Allocate space for the arguments of sqopt and s5solv.
*     These are for the data,
*        locA, indA, Acol, bl, bu, Names
*     and for the solution
*        hs, x, pi, rc, hs.
*     The true values of m, n and ne are also found.

      lAcol  = iw(256) ! Jcol(ne)    = Constraint Jacobian by columns
      lindA  = iw(258) ! indJ(ne) holds the row indices for Jij
      llocA  = iw(257) ! locJ(n+1)   = column pointers for indJ
      lbl    = iw(271) ! bl(nb)      = lower bounds
      lbu    = iw(272) ! bu(nb)      = upper bounds
      lx     = iw(299) ! x(nb)       = the solution (x,s)
      lpi    = iw(279) ! pi(m)       = the pi-vector
      lhs    = iw(282) ! the column state vector
      lNames = iw(359) ! Names(nName)

*     These are used as work arrays

      lhrtyp    = iw(284) ! hhrtyp(mBS), feasibility types
      lkBS      = iw(292) ! kBS(mBS)    = ( B  S ) list
      lkeynm    = iw(360) ! keynm(lenh) = hash table keys

      call s3mps
     &   ( iExit, key, maxm, maxn, maxne, lenh, nnL, nnCon, nnJac,
     &     m, n, ne, iw(negCon), nS,
     &     iObj, ObjAdd,
     &     iw(llocA), iw(lindA), rw(lAcol),
     &     rw(lbl), rw(lbu), cw(lNames),
     &     iw(lhs), iw(lhrtyp), iw(lkBS), iw(lkeynm),
     &     rw(lx), rw(lpi),
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (iExit .ne. 0) go to 700

*     negCon  counted the actual Jacobian entries in the MPS file,
*     but for Jacobian = Dense, we have to reset it.

      if (lDenJ .eq. 1) iw(negCon) = nnCon*nnJac

*     ------------------------------------------------------------------
*     Compress storage, now that we know the size of everything.
*     ------------------------------------------------------------------
*     Save current positions of  bl, bu, etc.

      kNames = lNames

      klocA  = llocA
      khs    = lhs

      kbl    = lbl
      kbu    = lbu
      kx     = lx
      kpi    = lpi

      nextcw = newcw
      nextiw = newiw
      nextrw = newrw

      call s3MapM
     &   ( iExit, m, n, ne, lenh,
     &     nextcw, nextiw, nextrw,
     &     mincw, miniw, minrw,
     &     lencw, leniw, lenrw, iw )
      if (iExit .ne. 0) go to 700

      lAcol  = iw(256) ! Jcol(ne)    = Constraint Jacobian by columns
      lindA  = iw(258) ! indJ(ne) holds the row indices for Jij
      llocA  = iw(257) ! locJ(n+1)   = column pointers for indJ
      lbl    = iw(271) ! bl(nb)      = lower bounds
      lbu    = iw(272) ! bu(nb)      = upper bounds
      lpi    = iw(279) ! pi(m)       = the pi-vector
      lhs    = iw(282) ! the column state vector
      lhrtyp = iw(284) ! hhrtyp(mBS), feasibility types
      lkBS   = iw(292) ! kBS(mBS)    = ( B  S ) list
      lx     = iw(299) ! x(nb)       = the solution (x,s)
      lNames = iw(359) ! Names(nName)
      lkeynm = iw(360) ! keynm(lenh) = hash table keys

*     Move bl, bu, etc. into their final positions.

      call icopy ( n+1, iw(klocA) , 1, iw(llocA) , 1 )
      call chcopy( n+m, cw(kNames), 1, cw(lNames), 1 )
      call icopy ( n+m, iw(khs)   , 1, iw(lhs)   , 1 )

      call dcopy ( n+m, rw(kbl)   , 1, rw(lbl)   , 1 )
      call dcopy ( n+m, rw(kbu)   , 1, rw(lbu)   , 1 )
      call dcopy ( n+m, rw(kx)    , 1, rw(lx)    , 1 )
      call dcopy (   m, rw(kpi)   , 1, rw(lpi)   , 1 )

*     ------------------------------------------------------------------
*     Check for error exits.
*     ------------------------------------------------------------------
  700 if (iExit .eq. 82) then
         write(str, 9820) mincw
         call snPRNT( 13, str, iw, leniw )

      else if (iExit .eq. 83) then
         write(str, 9830) miniw
         call snPRNT( 13, str, iw, leniw )

      else if (iExit .eq. 84) then
         write(str, 9840) minrw
         call snPRNT( 13, str, iw, leniw )

      else if (iExit .eq. 112) then ! Problem estimates too small
         if (m .ge. maxm) then
            maxm   = m
         else
            maxn   = n
            maxne  = ne
         end if

         if (iMPS .ne. iStdi  .and.  iMPS .ne. iSpecs) then
            rewind iMPS
            go to 100
         end if
      end if

      if (iExit .eq. 112  .or.  iExit .eq. 113) then
*        ------------------------------------
*        Flush MPS file to the ENDATA card
*        if it is the same as the SPECS file.
*        ------------------------------------
         if (iMPS .eq. iSpecs) then
            do idummy = 1, 100000
               if (key .eq. lENDA) go to 900
               read(iMPS, '(a4)') key
            end do
         end if
      end if

*     Exit.

  900 if (iMPS .ne. iStdi  .and.  iMPS .ne. iSpecs) rewind iMPS

      return

 9820 format(' Total character workspace should be significantly',
     &       ' more than', i8)
 9830 format(' Total integer   workspace should be significantly',
     &       ' more than', i8)
 9840 format(' Total real      workspace should be significantly',
     &       ' more than', i8)

      end ! subroutine s3inpt

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3MapM
     &   ( iExit, m, n, ne, lenh,
     &     nextcw, nextiw, nextrw,
     &     mincw, miniw, minrw,
     &     lencw, leniw, lenrw, iw )

      implicit
     &     none
      integer
     &     iExit, lenh, m, mincw, miniw, minrw, n, ne,
     &     nextcw, nextiw, nextrw, lencw, leniw, lenrw, iw(leniw)

*     ==================================================================
*     s3MapM   is called from s3inpt to find the minimum storage needed
*     for mps input.
*
*     Storage is estimated using the dimensions:
*        m    , n    , ne
*
*     18 Feb 1998: First version.
*     27 Oct 2003: Current version of s3MapM.
*     ==================================================================
      integer
     &     nb, lNames, lindA, llocA, lhElas, lhs, lAcol, lbl, lbu,
     &     lx, lpi, lrc, lhfeas, lkBS, lkeynm
*     ------------------------------------------------------------------
*     Allocate space for the data,
*              locA, indA, Acol, hElast, bl, bu, Names
*     and for the solution
*              hs, x, pi, rc, hs.

      iExit   = 0
      nb      = n + m

      lNames  = nextcw
      nextcw  = lNames + nb

      lindA   = nextiw
      llocA   = lindA  + ne
      lhElas  = llocA  + n  +  1
      lhs     = lhElas + nb
      nextiw  = lhs    + nb

      lAcol   = nextrw
      lbl     = lAcol  + ne
      lbu     = lbl    + nb
      lx      = lbu    + nb
      lpi     = lx     + nb
      lrc     = lpi    + m
      nextrw  = lrc    + nb

*     Allocate arrays needed during and after MPS input.

      lhfeas  = nextiw
      lkBS    = lhfeas + nb
      lkeynm  = lkBS   + nb
      nextiw  = lkeynm + lenh

      iw(257) = llocA
      iw(258) = lindA
      iw(282) = lhs
      iw(283) = lhElas
      iw(359) = lNames

      iw(256) = lAcol
      iw(271) = lbl
      iw(272) = lbu
      iw(279) = lpi
      iw(280) = lrc
      iw(299) = lx

      iw(284) = lhfeas
      iw(292) = lkBS
      iw(360) = lkeynm

      mincw   = nextcw - 1
      miniw   = nextiw - 1
      minrw   = nextrw - 1

      if (     mincw .gt. lencw) then
         iExit = 82
      else if (miniw .gt. leniw) then
         iExit = 83
      else if (minrw .gt. lenrw) then
         iExit = 84
      end if

      end ! subroutine s3MapM

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3mps
     &   ( iExit, key, maxm, maxn, maxne, lenh, nnL, nnCon, nnJac,
     &     m, n, ne, negCon, nS,
     &     iObj, ObjAdd,
     &     locA, indA, Acol, bl, bu, Names,
     &     hs, hrtype, kBS, keynam,
     &     x, pi,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iExit, iObj, lencw, leniw, lenrw, maxm, maxn, maxne, lenh,
     &     nnL, nnCon, nnJac, m, n, ne, negCon, nS, keynam(lenh),
     &     kBS(maxm+maxn), indA(maxne), hs(maxm+maxn), hrtype(maxm),
     &     locA(maxn+1), iw(leniw)
      double precision
     &     ObjAdd, Acol(maxne), bl(maxm+maxn), bu(maxm+maxn)
      double precision
     &     x(maxn+maxm), pi(maxm), rw(lenrw)
      character
     &     key*4, Names(maxm+maxn)*8, cw(lencw)*8

*     ==================================================================
*     s3mps  does the work for MPSinp and s3inpt.
*
*     iExit  status
*     -----  ------
*      112   too many rows, columns or elements
*      113   no rows or columns, or iObj le nnCon
*
*     18 Feb 1998: First version.
*     03 Aug 2003: snPRNT adopted.
*     26 Oct 2003: Current version of s3mps.
*     ==================================================================
      character
     &     str*80
      integer
     &     line, lrow, lenNam, nA0, ncoll, ncard(6), iEr(20)
*     ------------------------------------------------------------------

*     ncard  counts the number of data records in each section.

      iExit  = 0
      ncoll  = 0
      key    = '    '
      call iload ( 6, ncoll, ncard, 1 )

*     ------------------------------------------------------------------
*     Input ROWS.
*     lrow   is the location of the first rowname in Names.
*            The column names will later go at the front.
*     lenNam is the initial length of Names,
*            i.e. the maximum no. of names allowed for.
*     ------------------------------------------------------------------
      lrow   = maxn + 1
      lenNam = maxn + maxm
      call s3mpsa
     &   ( iExit, iEr, line, maxm, maxn,
     &     iObj, ncoll, m,
     &     lrow, lenNam, lenh, nnL, nnCon, key, ncard,
     &     hrtype, keynam, Names,
     &     cw, lencw, iw, leniw )
      if (iExit .ne. 0) go to 900

*     ------------------------------------------------------------------
*     m  is now known.
*     Input COLUMNS, RHS, RANGES.
*     ------------------------------------------------------------------
      call s3mpsb
     &   ( iExit, iEr, line, maxn, maxne,
     &     lrow, lenNam, lenh,
     &     ncoll, iObj, ObjAdd, m, n, ne,
     &     nnL, nnCon, nnJac, negCon, nA0, key, ncard,
     &     hrtype, keynam, Names,
     &     locA, indA, Acol, bl, bu, kBS, pi,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (iExit .ne. 0) go to 900

*     ------------------------------------------------------------------
*     n  and  ne  are now known.
*     Move the row names to be contiguous with the column names.
*     Input BOUNDS.
*     ------------------------------------------------------------------
      if (lrow .gt. n+1) then
         call chcopy( m, Names(lrow), 1, Names(n+1), 1 )
      end if

      call s3mpsc
     &   ( iEr, line, m, n, nS, lenNam,
     &     key, ncard, Names,
     &     bl, bu, hs, x,
     &     cw, lencw, iw, leniw, rw, lenrw )

      write(str, 1300) lenh
      call snPRNT( 11, str, iw, leniw )
      write(str, 1310) ncoll
      call snPRNT(  1, str, iw, leniw )

      if (nA0 .gt. 0) then
         write(str, 1320) nA0
         call snPRNT( 1, str, iw, leniw )
      end if

      if (nnCon .gt. 0) then
         write(str, 1350) negCon
         call snPRNT( 1, str, iw, leniw )
      end if

      if (nnL .gt. 0  .or.  ncard(6) .gt. 0) then
         write(str, 1400) ncard(6)
         call snPRNT( 1, str, iw, leniw )
         write(str, 1410) nS
         call snPRNT( 1, str, iw, leniw )
      end if

  900 return

 1300 format(' Length of row-name hash table  ', i12)
 1310 format(' Collisions during table lookup ', i12)
 1320 format(' No. of rejected coefficients   ', i12)
 1350 format(' No. of Jacobian entries specified', i10)
 1400 format(' No. of INITIAL  bounds  specified', i10)
 1410 format(' No. of superbasics specified   ', i12)

      end ! subroutine s3mps

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3mpsa
     &   ( iExit, iEr, line, maxm, maxn,
     &     iObj, ncoll, m,
     &     lrow, lennm, lenh, nnL, nnCon, key, ncard,
     &     hrtype, keynam, Names,
     &     cw, lencw, iw, leniw )

      implicit
     &     none
      integer
     &     iExit, line, maxm, maxn, iObj, ncoll, m, lencw, leniw,
     &     lrow, lennm, lenh, nnL, nnCon, hrtype(maxm),
     &     iEr(20), ncard(6), keynam(lenh), iw(leniw)
      character
     &     key*4, Names(lennm)*8, cw(lencw)*8

*     ==================================================================
*     s3mpsa  inputs the name and rows sections of an MPS file.
*
*     Original version written by Keith Morris, Wellington, 1973.
*
*     15 Nov 1991: First version based on Minos routine m3mpsa.
*     19 Jul 1997: Thread-safe version.
*     24 Sep 1997: key and Names(*) are now character*8.
*     03 Aug 2003: snPRNT adopted.
*     03 Aug 2003: Current version of s3mpsa.
*     ==================================================================
      character
     &     str*100, id(3)*8, lNAME*4, lROWS*4, lCOLU*4, lEX*4, lGX*4,
     &     lLX*4, lNX*4, lXE*4, lXG*4, lXL*4, lXN*4, cBlank*8
      logical
     &     found, gotnm
      integer
     &     ia, iMPS, inform, it, jrow, MPSerr, mLst, mObj, mProb
      double precision
     &     Aelem(2)
*     ------------------------------------------------------------------
      parameter         (mProb     =  51)
      parameter         (mObj      =  52)
*     ------------------------------------------------------------------
      data               cBlank
     &                  /'        '/
      data               lNAME, lROWS, lCOLU
     &                  /'NAME','ROWS','COLU'/
      data               lEX,   lGX,   lLX,   lNX
     &                  /' E  ',' G  ',' L  ',' N  '/
      data               lXE,   lXG,   lXL,   lXN
     &                  /'  E ','  G ','  L ','  N '/
*     ------------------------------------------------------------------
      iMPS      = iw(123) ! MPS file
      MPSerr    = iw(106) ! maximum # errors in MPS data
      mLst      = iw(107) ! maximum # lines  of MPS data

      call s1page( 1, iw, leniw )
      call snPRNT( 1, ' MPS file', iw, leniw )
      call snPRNT( 1, ' --------', iw, leniw )

      inform    = 0
      iObj      = 0
      line      = 0
      m         = 0
      gotnm     = cw(mObj) .ne. cBlank
      call iload ( 20  , 0, iEr  , 1 )
      call iload ( lenh, 0, keynam, 1 )

*     Look for the NAME card.

   10 call s3read
     &     ( 1, iMPS, iw, leniw, line, 5, key, id, Aelem, inform )
      if (key .ne. lNAME) then
         if (iEr(1) .eq. 0) then
             iEr(1) = 1
             call snPRNT( 3, ' XXXX  Garbage before NAME card',
     &            iw, leniw )
         end if
         go to 10
      end if

      cw(mProb) = id(2)
      write(str, 5000) cw(mProb)
      call snPRNT( 2, str, iw, leniw )

*     Look for the ROWS card.

      call s3read
     &     ( 1, iMPS, iw, leniw, line, 5, key, id, Aelem, inform )
      inform = 0
      if (key .ne. lROWS) then
         iEr(1) = iEr(1) + 1
         call snPRNT( 3, ' XXXX  ROWS card not found', iw, leniw )
         go to 35
      end if

*     ==================================================================
*     Read the row names and check if the relationals are valid.
*     ==================================================================
   30 call s3read
     &     ( 1, iMPS, iw, leniw, line, mLst, key, id, Aelem, inform )
      if (inform .ne. 0) go to 110

   35 if      (key .eq. lGX  .or.  key .eq. lXG) then
         it  = -1
      else if (key .eq. lEX  .or.  key .eq. lXE) then
         it  =  0
      else if (key .eq. lLX  .or.  key .eq. lXL) then
         it  =  1
      else if (key .eq. lNX  .or.  key .eq. lXN) then
         it  =  2

*        Record objective name if we don't already have one.

         if (iObj .eq. 0) then
            if (.not. gotnm) then
               cw(mObj) = id(1)
               if (nnL .gt. 0) then
                  write(str, 1170) cw(mObj)
                  call snPRNT( 3, str, iw, leniw )
               end if
            end if

            if (id(1) .eq. cw(mObj)) then
               iObj     = m + 1
               ncard(1) = ncard(1) + 1
            end if
         end if
      else
         iEr(3) = iEr(3) + 1
         if (iEr(3) .le. MPSerr) then
            write(str, 1160) line, key, id
            call snPRNT( 3, str, iw, leniw )
         end if
         go to 30
      end if

*     ------------------------------------------------------------------
*     Look up the row name  id(1)  in the hash table.
*     ------------------------------------------------------------------
      call s3hash
     &   ( lenh, maxm, ncoll, id(1), 2, keynam, Names(lrow), ia, found )

*     Error if the row name was already there.
*     Otherwise, enter the new name into the hash table.

      if ( found ) then
         iEr(4) = iEr(4) + 1
         if (iEr(4) .le. MPSerr) then
            write(str, 1200) id
            call snPRNT( 3, str, iw, leniw )
         end if
      else
         m      = m + 1
         if (m .le. maxm) then
            jrow        = maxn + m
            keynam(ia)  = m
            Names(jrow) = id(1)
            hrtype(m)   = it
         end if
      end if
      go to 30

*     ==================================================================
*     Should be COLUMNS card.
*     ==================================================================
  110 if (key .ne. lCOLU) then
         iEr(1) = iEr(1) + 1
         call snPRNT( 3, ' XXXX  COLUMNS card not found', iw, leniw )
      end if

*     Error if no rows or too many rows.

      if (m .le. 0) then
         call snPRNT( 3, ' XXXX  No rows specified', iw, leniw )
         iEr(1) = iEr(1) + 1
         iExit  = 113
         go to 900

      else if (m .gt. maxm) then
         write(str, 3030) maxm, m
         call snPRNT( 3, str, iw, leniw )
         iEr(1) = iEr(1) + 1
         iExit  = 112
         go to 900

      end if

*     Warning if no objective row found.
*     Error if linear objective is ahead of nonlinear rows.

      if (iObj .eq. 0) then
         call snPRNT( 3,
     &        ' ===>  Warning - no linear objective selected',
     &        iw, leniw )

      else if (iObj .le. nnCon) then
         call snPRNT( 13,
     &        ' XXXX  The linear objective card      N '//cw(mObj),
     &        iw, leniw )
         call snPRNT(  3,
     &        ' XXXX  is out of place.    Nonlinear constraints',
     &        iw, leniw )
         call snPRNT(  3,
     &        ' XXXX  must be listed first in the ROWS section.',
     &        iw, leniw )
         iExit = 113
         go to 900

      end if

      write(str, 5100) m
      call snPRNT( 2, str, iw, leniw )

*     ------------------------------------------------------------------
*     Exit
*     ------------------------------------------------------------------
  900 return

 1000 format(' MPS file' / ' --------')
 1160 format(' XXXX  Illegal row type at line', i7, '... ', a4, a8)
 1170 format(' ===>  Note --- row  ', a8,
     &   '  selected as linear part of objective.')
 1200 format(' XXXX  Duplicate row name --', a8, ' -- ignored')
 3030 format(' XXXX  Too many rows.  Limit was', i8,
     &   4x, '  Actual number is', i8)
 5000 format(' Name   ', a8)
 5100 format(' Rows   ', i8)

      end ! subroutine s3mpsa

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3mpsb
     &   ( iExit, iEr, line, maxn, maxne,
     &     lrow, lennm, lenh,
     &     ncoll, iObj, ObjAdd, m, n, ne,
     &     nnL, nnCon, nnJac, negCon, nA0, key, ncard,
     &     hrtype, keynam, Names,
     &     locA, indA, Acol, bl, bu, kBS, pi,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iExit, iObj, lencw, leniw, lenrw, line, lrow, lennm, lenh,
     &     m, maxn, maxne, n, ne, negCon, ncoll, nnL, nnCon, nnJac,
     &     hrtype(m), indA(maxne), locA(maxn+1), kBS(m), iEr(20),
     &     ncard(6), keynam(lenh), iw(leniw)
      double precision
     &     ObjAdd, Acol(maxne), bl(maxn+m), bu(maxn+m), pi(m), rw(lenrw)
      character
     &     key*4, Names(lennm)*8, cw(lencw)*8

*     ==================================================================
*     s3mpsb inputs the COLUMNS, RHS and RANGES sections of an MPS file.
*
*     Original version written by Keith Morris, Wellington, 1973.
*
*     19 Jul 1997: Thread-safe version.
*     24 Sep 1997: key and Names(*) are now character*8.
*     03 Aug 2003: snPRNT adopted.
*     27 Oct 2003: Current version of s3mpsb.
*     ==================================================================
      character
     &     str*100, id(3)*8, RowNm*8, ColNm*8
      logical
     &     dense, found, gotnm
      integer
     &     i, ia, iMPS, inform, irow, jslack, k,
     &     lJac, lDenJ, MPSerr, mLst, mRhs, mRng, nA0, ne1
      double precision
     &     aij, Aijtol, arng, BigUpp, BigLow, bnd, brng, infBnd,
     &     Aelem(2), bStruc(2)
*     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero      = 0.0d+0)
      parameter         (mRhs      =  53) ! rhs name
      parameter         (mRng      =  54) ! range name
*     ------------------------------------------------------------------
      character          cBlank*8, lLAGR*8
      character          lRHS*4  , lRHSX*4, lRANG*4
      data               cBlank               /'        '/
      data               lRHS  , lRHSX, lRANG /'RHS ','RHS''','RANG'/
      data               lLAGR                /'LAGRANGE'/
*     ------------------------------------------------------------------
      iMPS      = iw(123) ! MPS file

      lDenJ     = iw(105) ! 1(2)    => dense(sparse) Jacobian
      MPSerr    = iw(106) ! maximum # errors in MPS data
      mLst      = iw(107) ! maximum # lines  of MPS data

      infBnd    = rw( 70) ! definition of an infinite bound
      Aijtol    = rw( 95) ! zero Aij tolerance.
      bStruc(1) = rw( 96) ! default lower bound on x
      bStruc(2) = rw( 97) ! default upper bound on x

      ObjAdd    =   zero
      dense     =   lDenJ .eq. 1
      bigUpp    =   infBnd
      bigLow    = - bigUpp
      ColNm     =   '12345678'
      n         =   0
      na0       =   0
      ne        =   0
      ne1       = - 1
      negCon    =   0
      inform    =   0
      call iload ( m,    0, kBS, 1 )
      call dload ( m, zero, pi , 1 )

*     ==================================================================
*     Read the next columns card.
*     ==================================================================
  210 call s3read
     &   ( 2, iMPS, iw, leniw, line, mLst, key, id, Aelem, inform )
      if (inform .ne. 0) go to 310

  220 if (id(1) .ne. ColNm) then

*        Start a new column.

         if (ne .le. ne1) go to 310
         n         = n + 1
         ne1       = ne
         ColNm     = id(1)

         if (n .le. maxn) then
            locA(n)    = ne + 1
            Names(n) = ColNm

*           Make room for a dense Jacobian column.

            if (nnCon .gt. 0) then
               lJac   = ne
               if (dense  .and.  n .le. nnJac) then
                  ne     = ne + nnCon
                  if (ne .le. maxne) then
                     ne  = ne - nnCon
                     do i  = 1, nnCon
                        ne = ne + 1
                        indA(ne) = i
                        Acol(ne)  = zero
                     end do
                  end if
               end if
            end if
         end if
      end if

*     Process two row names and values.

      do 260, i  = 1, 2

*        Check for only one on the card.

         RowNm  = id(i+1)

         if (RowNm .ne. cBlank) then

*           Look up the row name.

            call s3hash
     &         ( lenh, m, ncoll,
     &           RowNm, 1, keynam,
     &           Names(lrow), ia, found )

            if ( found ) then
               aij    = aelem(i)
               irow   = keynam(ia)

*              Test for a duplicate entry.

               if (kBS(irow) .eq. n) then
                  iEr(8) = iEr(8) + 1
                  if (iEr(8) .le. MPSerr) then
                     write(str, 1420) ColNm, RowNm, Aij, line
                     call snPRNT( 1, str, iw, leniw )
                  end if
                  go to 260
               end if

               kBS(irow) = n
               if (irow .le. nnCon  .and.  n .le. nnJac) then

*                 Deal with Jacobian elements.

                  negCon = negCon + 1
                  if ( dense ) then
                     Acol(ne1 + irow) = aij
                     go to 260
                  end if

*                 Sparse Jacobian -- make sure the new element is
*                 squeezed in ahead of any linear-constraint elements.

                  lJac   = lJac + 1
                  if (lJac .le. ne) then
                     aij   = Acol(lJac)
                     irow  = indA(lJac)
                     Acol(lJac) = aelem(i)
                     indA(lJac) = keynam(ia)
                  end if

               else if (abs(aij) .lt. aijtol) then

*                 Ignore small aijs.

                  nA0    = nA0 + 1
                  go to 260
               end if

*              Pack the nonzero.

               ne     = ne + 1
               if (ne .le. maxne) then
                  indA(ne) = irow
                  Acol(ne)  = Aij
               end if
            else
               iEr(5) = iEr(5) + 1
               if (iEr(5) .le. MPSerr) then
                  write(str, 1400) RowNm, line
                  call snPRNT( 1, str, iw, leniw )
               end if
            end if
         end if ! RowNm ne Blank
  260 continue

      go to 210

*     Test for an empty column.

  310 if (ne .le. ne1) then

*        Column with no rows.   Warning unless variable is nonlinear.
*        Insert dummy column with zero in first row.

         if (n .gt. nnL) then
            iEr(6) = iEr(6) + 1
            if (iEr(6) .le. MPSerr) then
               write(str, 1500) ColNm
               call snPRNT( 1, str, iw, leniw )
            end if
         end if

         ne     = ne + 1
         if (ne .le. maxne) then
            indA(ne) = 1
            Acol(ne) = zero
         end if
         if (inform .eq. 0) go to 220
      end if

*     ==================================================================
*     See if we have hit the RHS.
*     ==================================================================
      if (key .ne. lRHS  .and.  key .ne. lRHSX) then

*        Nope sumpins rong.
*        Terminate the COLUMNS section anyway.

         iEr(7) = iEr(7) + 1
         call snPRNT( 3, ' XXXX  RHS card not found', iw, leniw )
      end if

*     Are there any columns at all?
*     Or too many columns or elements?

      if (n .le. 0) then
         call snPRNT( 3, ' XXXX  No columns specified', iw, leniw )
         iEr(2) = iEr(2) + 1
         iExit  = 113
         return

      else if (n .gt. maxn) then
         write(str, 3040) maxn, n
         call snPRNT( 3, str, iw, leniw )
         iEr(2) = iEr(2) + 1
         iExit  = 112
         return

      else if (ne .gt. maxne) then
         write(str, 3050) maxne, ne
         call snPRNT( 3, str, iw, leniw )
         iEr(2) = iEr(2) + 1
         iExit  = 112
         return
      end if

*     ------------------------------------------------------------------
*     Input the RHS.
*     ------------------------------------------------------------------
      write(str, 5200) n
      call snPRNT( 2, str, iw, leniw )
      write(str, 5210) ne
      call snPRNT( 2, str, iw, leniw )

*     We finally know how big the problem is.

      locA(n+1)   = ne + 1

*     Set bounds to default values.

      call dload ( n, bStruc(1), bl, 1 )
      call dload ( n, bStruc(2), bu, 1 )

      do i = 1, m
         k = hrtype(i)
         jslack = n + i

         if (k .le. 0) bl(jslack) = zero
         if (k .lt. 0) bu(jslack) = bigUpp
         if (k .gt. 0) bl(jslack) = bigLow
         if (k .ge. 0) bu(jslack) = zero
         if (k .eq. 2) bu(jslack) = bigUpp
      end do

*     Check for no RHS.

      if (key .ne. lRHS  .and.  key .ne. lRHSX) go to 600
      gotnm  = cw(mRhs) .ne. cBlank
      inform = 0

*     ==================================================================
*     Read next RHS card and see if it is the one we want.
*     ==================================================================
  420 call s3read
     &   ( 2, iMPS, iw, leniw, line, mLst, key, id, aelem, inform )
      if (inform .ne. 0) go to 600

*     A normal RHS is terminated if LAGRANGE is found.

      if (id(1) .eq. lLAGR) go to 490

      if (.not. gotnm) then
         gotnm     = .true.
         cw(mRhs) = id(1)
      end if

      if (id(1) .eq. cw(mRhs)) then

*        Look at both halves of the record.

         do i = 1, 2

            RowNm = id(i+1)

            if (RowNm .ne. cBlank) then
               call s3hash
     &            ( lenh, m, ncoll,
     &              RowNm, 1, keynam,
     &              Names(lrow), ia, found )

               if ( found ) then
                  ncard(2) = ncard(2) + 1
                  bnd      = aelem(i)
                  irow     = keynam(ia)
                  jslack   = n + irow
                  k        = hrtype(irow)

                  if (irow .eq. iObj) then
                     ObjAdd = bnd

                  else if (k .ne. 2) then
                     if (k .le. 0) bl(jslack) = bnd
                     if (k .ge. 0) bu(jslack) = bnd
                  end if

               else
                  iEr(5) = iEr(5) + 1
                  if (iEr(5) .le. MPSerr) then
                     write(str, 1400) RowNm, line
                     call snPRNT( 1, str, iw, leniw )
                  end if
               end if
            end if ! if RowNm ne Blank
         end do
      end if

      go to 420

*     LAGRANGE RHS found.

  490 if (ncard(2) .eq. 0) then
         cw(mRhs) = cBlank
         write(str, 1720)
         call snPRNT( 1, str, iw, leniw )
      end if
      go to 520

*     ==================================================================
*     Read next RHS card and see if it is a LAGRANGE one.
*     ==================================================================
  510 call s3read
     &   ( 2, iMPS, iw, leniw, line, mLst, key, id, aelem, inform )
      if (inform .ne. 0) go to 600

      if (id(1) .ne. lLAGR) go to 510

*     Find which row.
*     Look at both halves of the record.

  520 do i = 1, 2
         RowNm   = id(i+1)

         if (RowNm .ne. cBlank) then
            call s3hash
     &         ( lenh, m, ncoll,
     &           RowNm, 1, keynam,
     &           Names(lrow), ia, found )

            if ( found ) then
               ncard(5) = ncard(5) + 1
               irow     = keynam(ia)
               pi(irow) = aelem(i)
            else
               iEr(5) = iEr(5) + 1
               if (iEr(5) .le. MPSerr) then
                  write(str, 1400) RowNm, line
                  call snPRNT( 1, str, iw, leniw )
               end if
            end if
         end if ! RowNm ne Blank
      end do

      go to 510

*     ------------------------------------------------------------------
*     RHS has been input.
*     ------------------------------------------------------------------
  600 if (ncard(2) .eq. 0) then
         call snPRNT( 3, ' ===>  Warning - the RHS is zero', iw, leniw )
      end if

      if (ObjAdd .ne. zero) then
         write(str, 1630) ObjAdd
         call snPRNT( 1, str, iw, leniw )
      end if

*     ------------------------------------------------------------------
*     Input RANGES.
*     ------------------------------------------------------------------
*     Check for no RANGES.

      if (key .ne. lRANG) go to 800
      gotnm  = cw(mRng) .ne. cBlank
      inform = 0

*     ==================================================================
*     Read card and see if it is the range we want.
*     ==================================================================
  610 call s3read
     &   ( 2, iMPS, iw, leniw, line, mLst, key, id, aelem, inform )
      if (inform .ne. 0) go to 800

      if (.not. gotnm) then
         gotnm     = .true.
         cw(mRng) = id(1)
      end if

      if (id(1) .eq. cw(mRng)) then

*        Look at both halves of the record.

         do 640, i  = 1, 2
            RowNm   = id(i+1)

            if (RowNm .ne. cBlank) then
               call s3hash
     &            ( lenh, m, ncoll, RowNm, 1,
     &              keynam, Names(lrow), ia, found )

               if ( found ) then
                  ncard(3) = ncard(3)+1
                  brng     = aelem(i)
                  arng     = abs( brng )
                  irow     = keynam(ia)
                  jslack   = n + irow
                  k        = hrtype(irow)

                  if      (k .eq. 0) then
                     if (brng .gt. zero) bu(jslack) = bl(jslack) + arng
                     if (brng .lt. zero) bl(jslack) = bu(jslack) - arng

                  else if (k .eq. 2) then
*                    Relax

                  else if (k .lt. 0) then
                     bu(jslack) = bl(jslack) + arng

                  else if (k .gt. 0) then
                     bl(jslack) = bu(jslack) - arng
                  end if

               else
                  iEr(5) = iEr(5) + 1
                  if (iEr(5) .le. MPSerr) then
                     write(str, 1400) RowNm, line
                     call snPRNT( 1, str, iw, leniw )
                  end if
               end if
            end if ! RowNm ne Blank
  640    continue
      end if

      go to 610

*     RANGES have been input.

  800 return


 1400 format(' XXXX  Non-existent row    specified -- ', a8,
     &   ' -- entry ignored in line', i7)
 1420 format(' XXXX  Column  ', a8, '  has more than one entry',
     &       ' in row  ', a8
     &     / ' XXXX  Coefficient', 1p, e15.5, '  ignored in line', i10)
 1500 format(' XXXX  No valid row entries in column  ', a8)
 1630 format(' ===>  Note:  constant', 1p, e15.7,
     &   '  is added to the objective.')
 1720 format(' ===>  Warning - first RHS is LAGRANGE.',
     &       '   Other RHS''s will be ignored.')
 3040 format(' XXXX  Too many columns.   The limit was', i8,
     &   4x, '  Actual number is', i8)
 3050 format(' XXXX  Too many elements.  The limit was', i8,
     &   4x, '  Actual number is', i8)
 5200 format(' Columns', i8)
 5210 format(' Elements', i7)

      end ! subroutine s3mpsb

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3mpsc
     &   ( iEr, line, m, n, nS, lennm,
     &     key, ncard, Names,
     &     bl, bu, hs, x,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     lencw, leniw, lenrw, lennm, line, m, n, nS,
     &     hs(n+m), iEr(20), ncard(6), iw(leniw)
      double precision
     &     bl(n+m), bu(n+m), x(n+m), rw(lenrw)
      character
     &     key*4, Names(lennm)*8, cw(lencw)*8

*     ==================================================================
*     s3mpsc inputs the BOUNDS section of an MPS file.
*
*     Original version written by Keith Morris, Wellington, 1973.
*
*     19 Jul 1997: Thread-safe version.
*     24 Sep 1997: key and Names(*) are now character*8.
*     01 Aug 2003: s4name now needs iw, leniw for snPRNT.
*     03 Aug 2003: snPRNT adopted.
*     08 Oct 2003: Current version of s3mpsc.
*     ==================================================================
      logical
     &     gotnm, ignore
      character
     &     str*100, id(3)*8, Objtyp(3)*8, cBlank*8, lINIT*8, lBOUN*4,
     &     lENDA*4, lFR*4, lFX*4, lLO*4, lMI*4, lPL*4, lUP*4
      integer
     &     i, iLoadB, iMPS, iOldB, inform, iInsrt, iPrint, j,
     &     jmark, js, k, minmax, mLst, MPSerr, mObj, mRhs, mRng, mBnd
      double precision
     &     b1, b2, bigLow, bigUpp, bnd, infBnd, aelem(2)
*     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero      = 0.0d+0)
      parameter         (mObj      =  52)
      parameter         (mRhs      =  53)
      parameter         (mRng      =  54)
      parameter         (mBnd      =  55)
*     ------------------------------------------------------------------
      data
     &     cBlank
     &     /'        '/
      data
     &      lINIT
     &     /'INITIAL '/
      data
     &      lBOUN, lENDA
     &     /'BOUN','ENDA'/
      data
     &       lFR  , lFX, lLO
     &     /' FR ',' FX ',' LO '/
      data
     &       lMI  , lPL, lUP
     &     /' MI ',' PL ',' UP '/
      data
     &     Objtyp
     &     /'Max     ', 'Feas    ', 'Min     '/
*     ------------------------------------------------------------------
      iPrint    = iw( 12) ! Print file
      minmax    = iw( 87) ! 1, 0, -1  => MIN, FP, MAX
      MPSerr    = iw(106) ! maximum # errors in MPS data
      mLst      = iw(107) ! maximum # lines  of MPS data
      iLoadB    = iw(122) ! load file
      iMPS      = iw(123) ! MPS file
      iInsrt    = iw(125) ! insert file
      iOldB     = iw(126) ! old basis file

      infBnd    = rw( 70) ! definition of an infinite bound

      inform    =   1
      bigUpp    =   infBnd
      bigLow    = - bigUpp

*     Check for no BOUNDS.

      if (key .ne. lBOUN) go to 700
      gotnm     = cw(mBnd) .ne. cBlank
      inform    = 0
      jmark     = 1

*     ==================================================================
*     Read and check BOUNDS cards.  Notice the double plural.
*     ==================================================================
  610 call s3read
     &   ( 3, iMPS, iw, leniw, line, mLst, key, id, aelem, inform )
      if (inform .ne. 0) go to 700

*     A normal bounds set is terminated if INITIAL is found.

      bnd     = aelem(1)
      if (id(1) .eq. lINIT) go to 690

      if (.not. gotnm) then
         gotnm    = .true.
         cw(mBnd) = id(1)
      end if

      if (id(1) .eq. cw(mBnd)) then

*        Find which column.

         call s4name
     &      ( n, Names, id(2),
     &        line, iEr(10), 0, 1, n, jmark, j, iw, leniw )

         if (j .le. 0) then
            if (iEr(10) .le. MPSerr) then
               write(str, 1400) id(2), line
               call snPRNT( 1, str, iw, leniw )
            end if
         else

*           Select bound type for column j.

            ncard(4) = ncard(4) + 1
            if      (key .eq. lUP) then
               bu(j)  = bnd

            else if (key .eq. lLO) then
               bl(j)  = bnd

            else if (key .eq. lFX) then
               bu(j)  = bnd
               bl(j)  = bnd

            else if (key .eq. lFR) then
               bu(j)  = bigUpp
               bl(j)  = bigLow

            else if (key .eq. lMI) then
               if (bu(j) .ge. bigUpp) bu(j)  = zero
               bl(j)  = bigLow

            else if (key .eq. lPL) then
               bu(j)  = bigUpp

            else
*              This lad didn't even make it to Form 1.

               iEr(11) = iEr(11) + 1
               if (iEr(11) .le. MPSerr) then
                  write(str, 1700) line, key, (id(i), i=1,2)
                  call snPRNT( 1, str, iw, leniw )
               end if
            end if
         end if
      end if

      go to 610

*     INITIAL bounds set found.

  690 if (ncard(4) .eq. 0) then
         cw(mBnd) = cBlank
         write(iPrint, 1720)
         call snPRNT( 1, str, iw, leniw )
      end if

*     ------------------------------------------------------------------
*     End of normal bounds.
*     ------------------------------------------------------------------
  700 nS     =     0
      bigUpp =     bigUpp*0.9d+0
      bigLow =   - bigUpp

*     Set variables to be nonbasic at zero (as long as that's feasible).
*     All variables will be eligible for the initial basis.

      do j = 1, n+m
         x(j)  = max(  zero, bl(j) )
         x(j)  = min( x(j), bu(j) )
         hs(j)  = 0
         if (x(j) .eq. bu(j)) hs(j) = 1
      end do

*     Ignore INITIAL bounds if a basis will be loaded.

      if (inform .ne. 0) go to 790
      ignore = ioldB .gt. 0  .or.  iInsrt .gt. 0  .or.  iLoadB .gt. 0
      if (.not. ignore ) then
         jmark  = 1
         go to 720
      end if

*     ==================================================================
*     Read INITIAL bounds set.
*     ==================================================================
  710 call s3read
     &   ( 3, iMPS, iw, leniw, line, mLst, key, id, aelem, inform )
      if (inform .ne. 0) go to 790

      bnd    = aelem(1)
      if (ignore  .or.  id(1) .ne. lINIT) go to 710

*     Find which column.

  720 call s4name
     &   ( n, Names, id(2),
     &     line, iEr(12), 0, 1, n, jmark, j, iw, leniw )

      if (j .le. 0) then
         if (iEr(12) .le. MPSerr) then
            write(str, 1400) id(2), line
            call snPRNT( 1, str, iw, leniw )
         end if
      else

*        Select bound type for column j.

         ncard(6) = ncard(6)+1
         if      (key .eq. lFR) then
            js  = -1
         else if (key .eq. lFX) then
            js  =  2
            nS  = nS + 1
         else if (key .eq. lLO) then
            js  =  4
            bnd = bl(j)
         else if (key .eq. lUP) then
            js  =  5
            bnd = bu(j)
         else if (key .eq. lMI) then
            js  =  4
         else if (key .eq. lPL) then
            js  =  5
         else
            iEr(13) = iEr(13) + 1
            if (iEr(13) .le. MPSerr) then
               write(str, 1700) line, key, (id(i), i=1,2)
               call snPRNT( 1, str, iw, leniw )
            end if
            go to 710
         end if
      end if

      if (abs( bnd ) .ge. bigUpp) bnd = zero
      x (j)  = bnd
      hs(j)  = js
      go to 710

*     Should be ENDATA card.

  790 if (key .ne. lENDA) then
         iEr(14) = 1
         call snPRNT( 3, ' XXXX  ENDATA card not found', iw, leniw )
      end if

*     ------------------------------------------------------------------
*     Pass the Buck - not got to Truman yet.
*     ------------------------------------------------------------------
*     Check that  bl .le. bu

      do j  = 1, n
         b1 = bl(j)
         b2 = bu(j)
         if (b1 .gt. b2) then
            iEr(20) = iEr(20) + 1
            if (iEr(20) .le. MPSerr) then
               write(str, 1740) j, b1, b2
               call snPRNT( 1, str, iw, leniw )
            end if
            bl(j) = b2
            bu(j) = b1
         end if
      end do

*     Count the errors.

      k = 0
      do i = 1, 20
         k = k + iEr(i)
      end do
      if (k .gt. 0) then
         write(str, 1900) k
         call snPRNT( 13, str, iw, leniw )
      end if

      call snPRNT( 11, ' Names selected', iw, leniw )
      call snPRNT(  1, ' --------------', iw, leniw )
      write(str, 2100) cw(mObj), Objtyp(minmax+2)(1:3), ncard(1)
      call snPRNT(  1, str, iw, leniw )
      write(str, 2110) cw(mRhs), ncard(2)
      call snPRNT(  1, str, iw, leniw )
      write(str, 2120) cw(mRng), ncard(3)
      call snPRNT(  1, str, iw, leniw )
      write(str, 2130) cw(mBnd), ncard(4)
      call snPRNT(  1, str, iw, leniw )

      return


 1400 format(' XXXX  Non-existent column specified -- ', a8,
     &       ' -- entry ignored in line', i7)
 1700 format(' XXXX  Illegal bound type at line', i7, '... ',
     &       a4, a8, 2x, a8)
 1720 format(' ===>  Warning - first bounds set is  INITIAL .',
     &       '   Other bounds will be ignored.')
 1740 format(' XXXX  Bounds back to front on column', i6,' :',
     &       1p, 2e15.5)
 1900 format(' XXXX  Total no. of errors in MPS file', i6)
 2100 format(' Objective', 6x, a8, ' (', a3, ')', i8)
 2110 format(' RHS      ', 6x, a8, i14)
 2120 format(' RANGES   ', 6x, a8, i14)
 2130 format(' BOUNDS   ', 6x, a8, i14)

      end ! subroutine s3mpsc

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s3read
     &   ( mode, iMPS, iw, leniw, line, mxlist,
     &     key, id, Aelem, inform )

      implicit
     &     none
      integer
     &     mode, iMPS, leniw, iw(leniw), line, mxlist, inform
      double precision
     &     Aelem(2)
      character
     &     key*4, id(3)*8

*     ==================================================================
*     s3read  reads data from file iMPS and prints a listing on file
*     iPrint.  The data is assumed to be in MPS format, with items of
*     interest in the following six fields...
*
*     Field:     1         2         3         4         5         6
*
*     Columns: 01-04     05-12     15-22     25-36     40-47     50-61
*
*     Format:    a4        a8        a8      e12.0       a8      e12.0
*
*     Data:     key      id(1)     id(2)   Aelem(1)    id(3)   Aelem(2)
*
*
*     Comments may contain a * in column 1 and anything in columns 2-61.
*     They are listed and then ignored.
*
*
*     On entry,  mode    specifies which fields are to be processed.
*     On exit ,  inform  is set to 1 if column 1 is not blank.
*
*     15 Nov 1991: First version based on Minos routine m3read.
*     24 Sep 1997: key and id(*) are now character*8.
*     03 Aug 2003: snPRNT adopted.
*     03 Aug 2003: Current version of s3read.
*     ==================================================================
      character
     &     str*100, buffer*61, buff1*1
      integer
     &     last
*     ------------------------------------------------------------------
      character          lblank*1,    lstar*1
      data               lblank/' '/, lstar/'*'/
*     ------------------------------------------------------------------
*     Read a data card and look for keywords and comments.
*     ------------------------------------------------------------------
   10 read(iMPS, 1000) buffer
      buff1  = buffer(1:1)
      line   = line + 1

*     Print the buffer if column 1 is nonblank
*     or if a listing is wanted.

      if (buff1 .ne. lblank  .or.  line .le. mxlist) then

*        Find the last nonblank character.

         do 20, last = 61, 2, -1
            if (buffer(last:last) .ne. lblank) go to 30
   20    continue
         last   = 1

   30    write(str, 2000) line, buffer(1:last)
         call snPRNT( 1, str, iw, leniw )
      end if

*     Ignore comments.

      if (buff1 .eq. lstar ) go to 10

*     If column 1 is nonblank, load key and exit.
*     The NAME card is unusual in having some data in field 3.
*     We have to load it into id(2).

      if (buff1 .ne. lblank) then
         read(buffer, 1100) key, id(1), id(2)
         inform = 1
         return
      end if

*     ------------------------------------------------------------------
*     Process normal data cards.
*     ------------------------------------------------------------------
      if (mode .eq. 1) then

*        NAME or ROWS sections.

         read(buffer, 1100) key, id(1), id(2)
      else if (mode .eq. 2) then

*        COLUMNS, RHS or RANGES sections.

         read(buffer, 1100) key, id(1), id(2), aelem(1), id(3), aelem(2)
      else

*        BOUNDS section.

         read(buffer, 1100) key, id(1), id(2), aelem(1)
      end if

      return

 1000 format(a61)
 1100 format(a4, a8, 2x, a8, 2x, bn, e12.0, 3x, a8, 2x, e12.0)
 2000 format(i7, 4x, a)

      end ! subroutine s3read
