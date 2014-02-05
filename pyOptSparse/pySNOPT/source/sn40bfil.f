*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  sn40bfil.f
*
*     s4getB   s4id     s4ksav   s4name   s4inst   s4load   s4oldB
*     s4chek   s4chkP   s4dump   s4newB   s4pnch   s4rept   s4savB
*     s4soln   s4solp   s4stat
*
*     09 Mar 2004: snSolF implemented and called by s4soln.
*     17 Jun 2004: s4savb saves biggest primal and dual infeasibilities
*                  before and after scaling, so GAMS can use them.
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s4getB
     &   ( iExit, m, n, nb, nName, nS, iObj,
     &     hs, bl, bu, x, Names, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iExit, iObj, m, n, nb, nName, nS, leniw, lenrw, hs(nb),
     &     iw(leniw)
      double precision
     &     bl(nb), bu(nb), x(nb), rw(lenrw)
      character
     &     Names(nName)*8

*     ==================================================================
*     s4getb loads one of the basis files.
*
*     15 Nov 1991: First version based on Minos routine m4getb.
*     20 Apr 1999: Current version of s4getb.
*     ==================================================================
      integer
     &     iLoadB, iInsrt, iOldB
*     ------------------------------------------------------------------
      iLoadB    = iw(122) ! load file
      iInsrt    = iw(125) ! insert file
      iOldB     = iw(126) ! old basis file

*     Load a basis file if one exists and istart = 0 (Cold start).

      if (iOldB .gt. 0) then
         call s4oldB
     &      ( iExit, m, n, nb, nS, hs, bl, bu, x,
     &        iw, leniw, rw, lenrw )

      else if (iInsrt .gt. 0) then
         call s4inst
     &      ( n, nb, nS, iObj,
     &        hs, bl, bu, x, Names,
     &        iw, leniw, rw, lenrw )

      else if (iLoadB .gt. 0) then
         call s4load
     &      ( n, nb, nS, iObj,
     &        hs, x, Names,
     &        iw, leniw, rw, lenrw )
      end if

      end ! subroutine s4getB

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s4id
     &   ( j, n, nb, nName, Names, id )

      implicit
     &     none
      integer
     &     j, n, nb, nName
      character
     &     id*8, Names(nName)*8

*     ==================================================================
*     s4id   returns a name id for the j-th variable.
*     If nName = nb, the name is already in Names.
*     Otherwise nName = 1. Some generic column or row name is cooked up
*
*     15 Nov 1991: First version based on Minos 5.4 routine m4id.
*     16 Sep 1997: Current version.
*     ==================================================================
      integer
     &     i
*     ------------------------------------------------------------------
      character          ColNm*1, RowNm*1
      data               ColNm /'x'/,
     &                   RowNm /'r'/
*     ------------------------------------------------------------------
      if (nName .eq. nb) then
         id  = Names(j)
      else if (j .le. n) then
         write(id, '(a1,i7)') ColNm, j
      else
         i   = j - n
         write(id, '(a1,i7)') RowNm, i
      end if

      end ! subroutine s4id

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s4ksav
     &   ( minimz, m, n, nb, nS, mBS,
     &     itn, nInf, sInf, f, kBS, hs,
     &     Ascale, bl, bu, x, xBS, cw, lencw, iw, leniw )

      implicit
     &     none
      integer
     &     minimz, m, n, nb, nS, mBS, itn, nInf, lencw, leniw,
     &     hs(nb), kBS(mBS), iw(leniw)
      double precision
     &     sInf, f, Ascale(nb), bl(nb), bu(nb), xBS(mBS), x(nb)
      character
     &     cw(lencw)*8

*     ==================================================================
*     s4ksav  saves various quantities as determined by the frequency
*     control ksav.
*
*     15 Nov 1991: First version.
*     20 Apr 1999: Current version.
*     ==================================================================
      character
     &     istate(3)*4
      integer
     &     iBack, iNewB, itnlim, k
*     ------------------------------------------------------------------
      integer            Freq
      parameter         (Freq = 0)
*     ------------------------------------------------------------------
      iBack     = iw(120) ! backup file
      iNewB     = iw(124) ! new basis file
      itnlim    = iw( 89) ! limit on total iterations

      if (iNewB .gt. 0  .and.  itn .lt. itnlim) then
         k = 0
         call s4stat
     &      ( k, istate )
         call s4newB
     &      ( Freq, iNewB, minimz, m, n, nb,
     &        nS, mBS, itn, nInf, sInf, f, kBS, hs,
     &        Ascale, bl, bu, x, xBS, istate,
     &        cw, lencw, iw, leniw )

         if (iBack .gt. 0)
     &   call s4newB
     &      ( Freq, iBack, minimz, m, n, nb,
     &        nS, mBS, itn, nInf, sInf, f, kBS, hs,
     &        Ascale, bl, bu, x, xBS, istate,
     &        cw, lencw, iw, leniw )
      end if

      end ! subroutine s4ksav

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s4name
     &   ( n, Names, id,
     &     ncard, notfnd, maxmsg, j1, j2, jmark, jfound, iw, leniw )

      implicit
     &     none
      integer
     &     j1, j2, jmark, jfound, maxmsg, n, ncard, notfnd,
     &     leniw, iw(leniw)
      character
     &     Names(n)*8, id*8

*     ==================================================================
*     s4name searches for names in the array  Names(j), j = j1, j2.
*     jmark  will probably speed the search on the next entry.
*     Used by subroutines s3mpsc, s4inst, s4load.
*
*     Left-justified alphanumeric data is being tested for a match.
*
*     15 Nov 1991: First version based on Minos 5.4 routine m4name.
*     01 Aug 2003: snPRNT adopted.
*     08 Oct 2003: Current version of s4name.
*     ==================================================================
      character
     &     str*60
      integer
     &     j
*     ------------------------------------------------------------------
      do j = jmark, j2
         if (id .eq. Names(j)) go to 100
      end do

      do j = j1, jmark
         if (id .eq. Names(j)) go to 100
      end do

*     Not found.

      jfound = 0
      jmark  = j1
      notfnd = notfnd + 1
      if (notfnd .le. maxmsg) then
         write(str, 1000) ncard, id
         call snPRNT( 1, str, iw, leniw )
      end if

      return

*     Got it.

  100 jfound = j
      jmark  = j
      return

 1000 format(' XXX  Line', i6, '  --  name not found:', 8x, a8)

      end ! subroutine s4name

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s4inst
     &   ( n, nb, nS, iObj,
     &     hs, bl, bu, x, Names, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iObj, leniw, lenrw, n, nb, nS, hs(nb), iw(leniw)
      double precision
     &     bl(nb), bu(nb), x(nb), rw(lenrw)
      character
     &     Names(nb)*8

*     ==================================================================
*     This impression of INSERT reads a file produced by  s4pnch.
*     It is intended to read files similar to those produced by
*     standard MPS systems.  It recognizes SB as an additional key.
*     Also, values are extracted from columns 25--36.
*
*     15 Nov 1991: First version based on Minos 5.4 routine m4inst.
*     01 Aug 2003: snPRNT adopted.
*     08 Oct 2003: Current version of s4inst.
*     ==================================================================
      character
     &     Name1*8, Name2*8, id(5)*4, key*4, str*80
      integer
     &     ignord, iInsrt, iStdi, j, jmark, jObj,
     &     l, l1, lmark, MPSerr, nBS, ncard, ndum, nloop, notfnd
      double precision
     &     bplus, xj, infBnd
*     ------------------------------------------------------------------
      character          lLL*4, lUL*4, lXL*4, lXU*4, lSB*4, lEND*4
      data               lLL  , lUL  , lXL  , lXU  , lSB  , lEND
     &                 /' LL ',' UL ',' XL ',' XU ',' SB ', 'ENDA'/
*     ------------------------------------------------------------------
      infBnd    = rw( 70) ! definition of an infinite bound

      iStdi     = iw(  9) ! Standard Input
      iInsrt    = iw(125) ! insert file

      MPSerr    = iw(106) ! maximum # errors in MPS data

      bplus  = 0.9d+0*infBnd

      write(str, 1999) iInsrt
      call snPRNT( 13, str, iw, leniw )
      read(iInsrt, 1000) id
      write(str, 2000) id
      call snPRNT( 11, str, iw, leniw )

      l1     = n + 1

*     Make logicals basic.

      call iload ( nb-n, (3), hs(n+1), 1 )

      ignord = 0
      nBS    = 0
      nS     = 0
      notfnd = 0
      ncard  = 0
      jmark  = 1
      lmark  = l1
      ndum   = n + 100000

      if (iObj .eq. 0) then
         jObj = 0
      else
         jObj = n + iObj
      end if

*     Read names until ENDATA

      do 300, nloop = 1, ndum
         read(iInsrt, 1020) key, Name1, Name2, xj
         if (key .eq. lEND) go to 310

*        Look for  Name1.  It may be a column or a row,
*        since a superbasic variable could be either.

         ncard  = nloop
         call s4name
     &      ( nb, Names, Name1,
     &        ncard, notfnd, MPSerr, 1, nb, jmark, j, iw, leniw )
         if (   j  .le. 0) go to 300
         if (hs(j) .gt. 1) go to 290
         if (key .ne. lXL  .and.  key .ne. lXU) go to 70

*        Look for  Name2.  It has to be a row.

         call s4name
     &      ( nb, Names, Name2,
     &        ncard, notfnd, MPSerr, l1, nb, lmark, l, iw, leniw )
         if (l .le. 0) go to 300

*        XL, XU (exchange card)  --  make col j basic,  row l nonbasic.

         if (l  .eq. jObj) go to 290
         if (hs(l) .ne. 3) go to 290
         nBS    = nBS + 1
         hs(j)  = 3
         if (key .eq. lXU) go to 50
         hs(l)  = 0
         if (bl(l) .gt. -bplus) x(l) = bl(l)
         go to 250

   50    hs(l)  = 1
         if (bu(l) .lt.  bplus) x(l) = bu(l)
         go to 250

*        LL, UL, SB  --  only  j  and  xj  are relevant.

   70    if (key .eq. lLL) go to 100
         if (key .eq. lUL) go to 150
         if (key .eq. lSB) go to 200
         go to 290

*        LO or UP

  100    hs(j)  = 0
         go to 250

  150    hs(j)  = 1
         go to 250

*        Make superbasic.

  200    hs(j)  = 2
         nS     = nS + 1

*        Save  x  values.

  250    if (abs(xj) .lt. bplus) x(j) = xj
         go to 300

*        Card ignored.

  290    ignord = ignord + 1
         if (ignord .le. MPSerr) then
            write(str, 2010) ncard, key, Name1, Name2
            call snPRNT( 1, str, iw, leniw )
         end if
  300 continue

  310 ignord = ignord + notfnd
      write(str, 2050) ncard, ignord
      call snPRNT( 13, str, iw, leniw )
      write(str, 2060) nBS, nS
      call snPRNT(  3, str, iw, leniw )
      if (iInsrt  .ne. iStdi) rewind iInsrt
      return

 1000 format(14x, 2a4, 2x, 3a4)
 1020 format(a4, a8, 2x, a8, 2x, e12.5)
 1999 format(' INSERT file to be input from file', i4)
 2000 format(' NAME', 10x, 2a4, 2x, 3a4)
 2010 format(' XXX  Line', i6, '  ignored:', 8x, a4, a8, 2x, a8)
 2050 format(' No. of lines read      ', i6, '  Lines ignored', i6)
 2060 format(' No. of basics specified', i6, '  Superbasics  ', i6)

      end ! subroutine s4inst

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s4load
     &   ( n, nb, nS, iObj,
     &     hs, x, Names, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     n, nb, nS, iObj, leniw, lenrw, hs(nb), iw(leniw)
      double precision
     &     x(nb), rw(lenrw)
      character
     &     Names(nb)*8

*     ==================================================================
*     s4load  inputs a load file, which may contain a full or partial
*     list of row and column names and their states and values.
*     Valid keys are   BS, LL, UL, SB.
*
*     15 Nov 1991: First version based on Minos 5.4 routine m4load.
*     01 Aug 2003: snPRNT adopted.
*     01 Aug 2003: Current version of s4load.
*     ==================================================================
      character
     &     key*4, id(5)*4, Name*8, str*80
      integer
     &     ignord, iLoadB, iStdi, j, jmark, jObj, MPSerr,
     &     nBS, ncard, ndum, nloop, notfnd
      double precision
     &     bplus, xj, infBnd
*     ------------------------------------------------------------------
      character          lBS*4, lLL*4, lUL*4, lSB*4, lEND*4
      data               lBS  , lLL  , lUL  , lSB  , lEND
     &                 /' BS ',' LL ',' UL ',' SB ','ENDA'/
*     ------------------------------------------------------------------
      infBnd    = rw( 70) ! definition of an infinite bound
      iStdi     = iw(  9) ! Standard Input
      iLoadB    = iw(122) ! load file
      MPSerr    = iw(106) ! maximum # errors in MPS data

      bplus  = 0.9d+0*infBnd
      write(str, 1999) iLoadB
      call snPRNT( 13, str, iw, leniw )
      read(iLoadB, 1000) id
      write(str, 2000) id
      call snPRNT(  1, str, iw, leniw )
      ignord = 0
      nBS    = 0
      nS     = 0
      notfnd = 0
      ncard  = 0
      jmark  = 1
      ndum   = n + 100000

      if (iObj .eq. 0) then
         jObj = 0
      else
         jObj = n + iObj
      end if

*     Read names until ENDATA is found.

      do 300, nloop = 1, ndum
         read(iLoadB, 1020) key, Name, xj
         if (key .eq. lEND) go to 310

         ncard  = nloop
         call s4name
     &      ( nb, Names, Name,
     &        ncard, notfnd, MPSerr, 1, nb, jmark, j, iw, leniw )
         if (j .le. 0) go to 300

*        The name Name belongs to the j-th variable.

         if (hs(j) .gt. 1) go to 290
         if (j   .eq.jObj) go to  90
         if (key .eq. lBS) go to  90
         if (key .eq. lLL) go to 100
         if (key .eq. lUL) go to 150
         if (key .eq. lSB) go to 200
         go to 290

*        Make basic.

   90    nBS    = nBS + 1
         hs(j)  = 3
         go to 250

*        LO or UP.

  100    hs(j)  = 0
         go to 250

  150    hs(j)  = 1
         go to 250

*        Make superbasic.

  200    nS     = nS + 1
         hs(j)  = 2

*        Save  x  values.

  250    if (abs(xj) .lt. bplus) x(j) = xj
         go to 300

*        Card ignored.

  290    ignord = ignord + 1
         if (ignord .le. MPSerr) then
            write(str, 2010) ncard, key, Name
            call snPRNT( 1, str, iw, leniw )
         end if
  300 continue

  310 ignord = ignord + notfnd
      write(str, 2050) ncard, ignord
      call snPRNT( 13, str, iw, leniw )
      write(str, 2060) nBS, nS
      call snPRNT(  3, str, iw, leniw )

*     Make sure the linear Objective is basic.

      if (iObj  .gt. 0) then
         if (hs(jObj) .ne. 3) then
            hs(jObj) = 3

*           Swap Obj with last basic variable.

            do j = nb, 1, -1
               if (hs(j) .eq. 3) go to 860
            end do

  860       hs(j)  = 0
         end if
      end if

      if (iLoadB .ne. iStdi) rewind iLoadB
      return

 1000 format(14x, 2a4, 2x, 3a4)
 1020 format(a4, a8, 12x, e12.5)
 1999 format(' LOAD file to be input from file', i4)
 2000 format(' NAME', 10x, 2a4, 2x, 3a4)
 2010 format(' XXX  Line', i6, '  ignored:', 8x, a4, a8, 2x, a8)
 2050 format(' No. of lines read      ', i6, '  Lines ignored', i6)
 2060 format(' No. of basics specified', i6, '  Superbasics  ', i6)

      end ! subroutine s4load

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s4oldB
     &   ( iExit, m, n, nb, nS, hs, bl, bu, x,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iExit, leniw, lenrw, m, n, nb, nS, hs(nb), iw(leniw)
      double precision
     &     bl(nb), bu(nb), x(nb), rw(lenrw)

*     ==================================================================
*     s4oldB  inputs a compact basis file from file  iOldB.
*
*     15 Nov 1991: First version based on Minos 5.4 routine m4OldB.
*     01 Aug 2003: snExit and snPRNT adopted.
*     04 Dec 2004: Current version of s4oldB.
*     ==================================================================
      character
     &     str*81, id(20)*4
      integer
     &     i, j, js, newm, newn, idummy, ndummy, iStdi, iOldB
      double precision
     &     infBnd, bplus, xj
*     ------------------------------------------------------------------
      infBnd    = rw( 70) ! definition of an infinite bound

      iStdi     = iw(  9) ! Standard Input
      iOldB     = iw(126) ! old basis file

      bplus  = 0.9d+0*infBnd
      write(str, 1999) iOldB
      call snPRNT( 13, str, iw, leniw )
      read(iOldB , 1000) id
      write(str, 2000) id
      call snPRNT(  1, str, iw, leniw )

      read(iOldB , 1005) (id(i), i=1,13), newm, newn, nS
      write(str, 2005) (id(i), i=1,13), newm, newn, nS
      call snPRNT( 1, str, iw, leniw )

      if (newm .ne. m  .or.  newn .ne. n) go to 900
      read(iOldB , 1010) hs

*     Set values for nonbasic variables.

      do j = 1, nb
         js = hs(j)
         if (js .le. 1) then
            if (js .eq. 0) xj = bl(j)
            if (js .eq. 1) xj = bu(j)
            if (abs(xj) .lt. bplus) x(j) = xj
         end if
      end do

*     Load superbasics.

      nS     = 0
      ndummy = m + n + 10000

      do idummy = 1, ndummy
         read(iOldB, 1020) j, xj
         if (j .le.  0) go to 310
         if (j .le. nb) then
            x(j)  = xj
            if (hs(j) .eq. 2) nS = nS + 1
         end if
      end do

  310 write(str, 2010) nS
      call snPRNT( 1, str, iw, leniw )
      go to 990

*     Error exits.

*     -------------------------------------------
*     Incompatible basis file dimensions.
*     -------------------------------------------
  900 iExit = 92      ! Basis file dimensions do not match this problem

  990 if (iOldB .ne. iStdi) rewind iOldB
      return

 1000 format(20a4)
 1005 format(13a4, 2x, i7, 3x, i7, 4x, i5)
 1010 format(80i1)
 1020 format(i8, e24.14)
 1999 format(' OLD BASIS file to be input from file', i4)
 2000 format(1x, 20a4)
 2005 format(1x, 13a4,
     &       'm=', i7, ' n=', i7, ' sb=', i5)
 2010 format(' No. of superbasics loaded', i7)

      end ! subroutine s4oldB

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s4chek
     &   ( m, maxS, mBS, n, nb, needB, gotHes,
     &     nS, iObj, hs, kBS, bl, bu, x, iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     needB
      integer
     &     gotHes, m, maxS, mBS, n, nb, nS, iObj, leniw, lenrw,
     &     hs(nb), kBS(mBS), iw(leniw)
      double precision
     &     bl(nb), bu(nb), x(nb), rw(lenrw)

*     ==================================================================
*     s4chek  takes hs and x and checks they contain reasonable values.
*     The entries hs(j) = 2 are used to set  nS  and possibly
*     the list of superbasic variables kBS(m+1) thru kBS(m+nS).
*     Scaling, if any, has taken place by this stage.
*
*     15 Nov 1991: First version based on Minos routine m4chek.
*     01 Aug 2003: snPRNT adopted.
*     08 Mar 2004: needB, gotHes, setkS allow for Hot start.
*     08 Mar 2004: Current version of s4chek.
*     ==================================================================
      character
     &     str*80
      logical
     &     setkS
      integer
     &     j, jj, js, nBasic, nSsave
      double precision
     &     b1, b2, bplus, infBnd, xj
*     ------------------------------------------------------------------
      double precision   zero,           tolb
      parameter        ( zero = 0.0d+0,  tolb = 1.0d-4 )
*     ------------------------------------------------------------------
      infBnd    = rw( 70) ! definition of an infinite bound

      setkS  = needB  .or.  gotHes .le. 0

*     Make sure hs(j) = 0, 1, 2 or 3 only.

      do j = 1, nb
         js   = hs(j)
         if (js .lt. 0) hs(j) = 0
         if (js .ge. 4) hs(j) = js - 4
      end do

*     ------------------------------------------------------------------
*     Make sure the Objective is basic and free.
*     Then count the basics and superbasics, making sure they don't
*     exceed m and maxS respectively.  Also, set nS and possibly
*     kBS(m+1) thru kBS(m+ns) to define the list of superbasics.
*     Mar 1988: Loop 100 now goes backwards to make sure we grab Obj.
*     Apr 1992: Backwards seems a bit silly in the documentation.
*               We now go forward through the slacks,
*               then forward through the columns.
*     ------------------------------------------------------------------
  100 nBasic = 0
      nSsave = nS
      nS     = 0

      if (iObj .gt. 0) then
         hs(n+iObj) =   3
         bl(n+iObj) = - infBnd
         bu(n+iObj) =   infBnd
      end if

*     If too many basics or superbasics, make them nonbasic.
*     Do slacks first to make sure we grab the objective slack.

      j = n

      do jj = 1, nb
         j       = j + 1
         if (j .gt. nb) j = 1
         js      = hs(j)
         if (js .eq. 2) then
            nS   = nS + 1
            if (nS .le. maxS) then
               if ( setkS ) kBS(m+nS) = j
            else
               hs(j)     = 0
            end if

         else if (js .eq. 3) then
            nBasic = nBasic + 1
            if (nBasic .gt. m) hs(j) = 0
         end if
      end do

*     Proceed if the superbasic kBS were reset,
*     or if nS seems to agree with the superbasics in hs(*).
*     Otherwise, give up trying to save the reduced Hessian,
*     and reset the superbasic kBS after all.

      if ( setkS ) then
*        ok
      else if (nS .ne. nSsave) then
         setkS   = .true.
         gotHes  = 0
         write(str, 1000) nS, nSsave
         call snPRNT( 13, str, iw, leniw )
         go to 100
      end if

*     Check the number of basics.

      nS     = min( nS, maxS )
      if (nBasic .ne. m ) then
         write(str , 1100) nBasic, m
         call snPRNT( 13, str, iw, leniw )
      end if

*     ------------------------------------------------------------------
*     Set each nonbasic x(j) to be exactly on its
*     nearest bound if it is within tolb of that bound.
*     ------------------------------------------------------------------
      bplus = 0.1d+0 * infBnd
      do  j = 1, nb
         xj     = x(j)
         if (abs( xj ) .ge.  bplus) xj = zero
         if (hs(j)     .le.  1    ) then
            b1  = bl(j)
            b2  = bu(j)
            xj  = max( xj, b1 )
            xj  = min( xj, b2 )
            if (   (xj - b1) .gt. (b2 - xj)) b1 = b2
            if (abs(xj - b1) .le.  tolb    ) xj = b1
            hs(j) = 0
            if (xj .gt. bl(j)) hs(j) = 1
         end if
         x(j) = xj
      end do

      return

 1000 format(' WARNING:', i6, ' superbasics in hs(*);',
     $       ' previously nS =', i6, '.  Hessian not saved')
 1100 format(' WARNING:', i7, ' basics specified;',
     &       ' preferably should have been', i7)

      end ! subroutine s4chek

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s4chkP
     &   ( Errors, cPointr, iPointr, iw, leniw )

      implicit
     &     none
      integer
     &     Errors, iPointr, leniw, iw(leniw)
      character
     &     cPointr*6

!     ==================================================================
!     s4chkP checks the value of a pointer retrieved from workspace
!     Error messages are listed on the standard output
!
!     03 Sep 2006: First version of s4chkP.
!     ==================================================================
      character
     &     str*80
!     ------------------------------------------------------------------
      if (iPointr .lt. 0  .or. iPointr .gt. leniw) then
         Errors = Errors + 1
         write(str, 9999) cPointr, iPointr
         call snPRNT( 5, str, iw, leniw )
      end if

 9999 format(' XXX  Pointer out of range:  ', a6, ' = ', i6)

      end ! subroutine s4chkP

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s4dump
     &   ( iDump, m, n, nb, nName,
     &     hs, x, Names, cw, lencw, iw, leniw )

      implicit
     &     none
      integer
     &     iDump, m, n, nb, nName, lencw, leniw, hs(nb), iw(leniw)
      double precision
     &     x(nb)
      character
     &     Names(nName)*8, cw(lencw)*8

*     ==================================================================
*     s4dump outputs basis names in a format compatible with s4load.
*     This file is normally easier to modify than a punch file.
*
*     15 Nov 1991: First version based on Minos 5.4 routine m4dump.
*     01 Aug 2003: snPRNT adopted.
*     17 Jul 2005: Dump default names when nName = 1
*     ==================================================================
      integer
     &     iPrint, j, k
      character
     &     id*8, key(4)*4, mProb*8, str*40
      data
     &     key/' LL ', ' UL ', ' SB ', ' BS '/
*     ------------------------------------------------------------------
      iPrint    = iw( 12) ! Print file
      mProb     = cw( 51) ! Problem name

      write(iDump, 2000) mProb

      do j = 1, nb
         call s4id
     &      ( j, n, nb, nName, Names, id )

         k = hs(j) + 1
         write(iDump, 2100) key(k), id, x(j)
      end do

      write(iDump , 2200)
      write(str, 3000) iDump
      call snPRNT( 13, str, iw, leniw )
      if (iDump .ne. iPrint) rewind iDump
      return

 2000 format('NAME', 10x, a8, 2x, '   DUMP/LOAD')
 2100 format(a4, a8, 12x, 1p, e12.5)
 2200 format('ENDATA')
 3000 format(' DUMP file saved on file', i4)

      end ! subroutine s4dump

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s4newB
     &   ( job, iNewB, minimz, m, n, nb,
     &     nS, mBS, itn, nInf, sInf, f, kBS, hs,
     &     Ascale, bl, bu, x, xBS, istate,
     &     cw, lencw, iw, leniw )

      implicit
     &     none
      integer
     &     job, iNewB, minimz, m, n, nb, nS, mBS, itn, nInf,
     &     lencw, leniw, hs(nb), kBS(mBS), iw(leniw)
      double precision
     &     sInf, f, Ascale(nb), bl(nb), bu(nb), xBS(mBS), x(nb)
      character
     &     istate(3)*4, cw(lencw)*8

*     ==================================================================
*     s4newB  saves a compact basis on file iNewB.  Called from S5QP.
*     job = Freq, the save is a periodic one due to the save frequency.
*     job = Wrap, S5solv has just finished the current problem.
*
*     15 Nov 1991: First version based on Minos 5.4 routine m4newb.
*     16 Jan 2003: First line had "or" instead of "and".
*     01 Aug 2003: snPRNT adopted.
*     01 Apr 2005: Negative hs values converted in-situ.
*     17 Jul 2005: Blanks printed for undefined MPS names.
*     ==================================================================
      character
     &     mProb*8, mObj*8, mRhs*8, mRng*8, mBnd*8, str*60
      logical
     &     scaled
      integer
     &     buffer(80), i, iPrint, j, j1, j2, js, k, lvlScl, nnJac
      double precision
     &     Obj, xj
*     ------------------------------------------------------------------
      character          cdummy*8,            lblank*8
      parameter         (cdummy = '-1111111', lblank = '        ')

      integer            Freq,     Wrap
      parameter         (Freq = 0, Wrap = 1)
*     ------------------------------------------------------------------
      if (job .ne. Freq  .and.  job .ne. Wrap) return

      iPrint    = iw( 12) ! Print file
      nnJac     = iw( 21) ! # nonlinear Jacobian variables
      lvlScl    = iw( 75) ! scale option

      mProb     = cw( 51) ! Problem name
      mObj      = cw( 52) ! Objective name
      mRhs      = cw( 53) ! Right-hand side name
      mRng      = cw( 54) ! Range name
      mBnd      = cw( 55) ! Bnd section name

      scaled    = lvlScl .gt. 0

      if (mProb .eq. cdummy) then
         mProb = lblank
      end if
      if (mObj  .eq. cdummy) then
         mObj  = lblank
      end if
      if (mRhs  .eq. cdummy) then
         mRhs  = lblank
      end if
      if (mRng  .eq. cdummy) then
         mRng  = lblank
      end if
      if (mBnd  .eq. cdummy) then
         mBnd  = lblank
      end if

      if (nInf .eq. 0) then
         Obj = minimz * f
      else
         Obj = sInf
      end if

*     Output header cards and the state vector.

      write(iNewB, 1000) mProb, itn , istate, nInf, Obj
      write(iNewB, 1005) mObj , mRhs, mRng  , mBnd, m, n, nS
*     write(iNewB, 1010) hs

      j2   = 0
      do i = 1, nb, 80
         j1 = j2 + 1
         j2 = j2 + 80
         j2 = min(j2,nb)
         k  = 0
         do j = j1, j2
            js = hs(j)
            if (js .ge. 4  .or.  js .eq. -1) js = 0
            k  = k + 1
            buffer(k) = js
         end do
         write(iNewB, '(80i1)') (buffer(j), j = 1, k)
      end do

*     Output the superbasic variables.

      do k  = m+1, m+nS
         j  = kBS(k)
         xj = xBS(k)
         if (scaled) xj = xj * Ascale(j)
         write(iNewB, 1020) j, xj
      end do

*     Output the values of all other (non-SB) nonlinear variables.

      do j = 1, nnJac
         if (hs(j) .ne. 2) then
            xj    = x(j)
            if (scaled) xj = xj * Ascale(j)
            write(iNewB, 1020) j, xj
         end if
      end do

*     Output nonbasic variables that are not at a bound.

      do j = nnJac+1, nb
         js = hs(j)
         if (js .le. 1  .or.  js .ge. 4  .or.  js .eq. -1) then
            xj    = x(j)
            if (xj .ne. bl(j)) then
               if (xj .ne. bu(j)) then
                  if (scaled) xj = xj * Ascale(j)
                  write(iNewB, 1020) j, xj
               end if
            end if
         end if
      end do

*     Terminate the list with a zero.

      j     = 0
      write(iNewB, 1020) j
      if (iNewB .ne. iPrint) rewind iNewB
      write(str, 1030) iNewB, itn
      call snPRNT( 13, str, iw, leniw )
      return

 1000 format(a8, '  ITN', i8, 4x, 3a4, '  NINF', i7,
     &       '      OBJ', 1p, e21.12)
 1005 format('OBJ=',a8, ' RHS=',a8, ' RNG=',a8, ' BND=',a8,
     &       ' M=', i7,  ' N=', i7, ' SB=', i5)
 1010 format(80i1)
 1020 format(i8, 1p, e24.14)
 1030 format(' NEW BASIS file saved on file', i4, '    itn =', i7)

      end ! subroutine s4newB

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s4pnch
     &   ( iPnch, m, n, nb, nName,
     &     hs, bl, bu, x, Names, cw, lencw, iw, leniw )

      implicit
     &     none
      integer
     &     iPnch, m, n, nb, nName, lencw, leniw, hs(nb), iw(leniw)
      double precision
     &     bl(nb), bu(nb), x(nb)
      character
     &     Names(nName)*8, cw(lencw)*8

*     ==================================================================
*     s4pnch  outputs a PUNCH file (list of basis names, states and
*     values) in a format that is compatible with MPS/360.
*
*     15 Nov 1991: First version based on Minos 5.4 routine m4pnch.
*     16 Jan 2003: Allow for missing names.
*     01 Aug 2003: snPRNT adopted.
*     01 Aug 2003: Current version of s4pnch.
*     ==================================================================
      integer
     &     iPrint, irow, j, jN, k
      character
     &     ColNm*8, mProb*8, Name*8, str*60
*     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero = 0.0d+0)
*     ------------------------------------------------------------------
      character          lblank*8, key(5)*4
      data               key   /' LL ', ' UL ', ' SB ', ' XL ', ' XU '/
      data               lblank/'        '/
*     ------------------------------------------------------------------
      iPrint    = iw( 12) ! Print file
      mProb     = cw( 51) ! Problem name

      write(iPnch, 2000) mProb
      irow      = n

      do 500, j = 1, n
         jN     = j
         call s4id
     &      ( jN, n, nb, nName, Names, ColNm )
         k       = hs(j)

         if (k .eq. 3) then

*           Basics -- find the next row that isn't basic.

  300       irow   = irow + 1
            if (irow .le. nb) then
               k      = hs(irow)
               if (k .eq. 3) go to 300

               call s4id
     &              ( irow, n, nb, nName, Names, Name )
               if (k .eq. 2) k = 0
               write(iPnch, 2100) key(k+4), ColNm, Name, x(j)
            end if
         else

*           Skip nonbasic variables with zero lower bounds.

            if (k .le. 1) then
               if (bl(j) .eq. zero  .and.  x(j) .eq. zero) go to 500
            end if
            write(iPnch, 2100) key(k+1), ColNm, lblank, x(j)
         end if
  500 continue

*     Output superbasic slacks.

      do j = n+1, nb
         if (hs(j) .eq. 2) then
            jN   = j
            call s4id
     &           ( jN, n, nb, nName, Names, Name )
            write(iPnch, 2100) key(3), Name, lblank, x(j)
         end if
      end do

      write(iPnch , 2200)
      write(str, 3000) iPnch
      call snPRNT( 13, str, iw, leniw )
      if (iPnch .ne. iPrint) rewind iPnch
      return

 2000 format('NAME', 10x, a8, 2x, 'PUNCH/INSERT')
 2100 format(a4, a8, 2x, a8, 2x, 1p, e12.5)
 2200 format('ENDATA')
 3000 format(' PUNCH file saved on file', i4)

      end ! subroutine s4pnch

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s4rept
     &   ( ondisk, m, n, nb, nName,
     &     nnCon0, nnCon, nnObj0, nnObj, nS,
     &     ne, nlocA, locA, indA, Acol,
     &     hs, Ascale, bl, bu, gObj, pi, x, Fx,
     &     Names, istate, iw, leniw )

      implicit
     &     none
      logical
     &     ondisk
      integer
     &     nName, m, n, nb, ne, nlocA, nnCon0, nnCon, nnObj0, nnObj,
     &     nS, leniw, indA(ne), hs(nb), locA(nlocA), iw(leniw)
      double precision
     &     Acol(ne), Ascale(nb), bl(nb), bu(nb), gObj(nnObj0), pi(m),
     &     x(nb), Fx(nnCon0)
      character
     &     Names(nName)*8, istate(3)*4

*     ==================================================================
*     s4rept  has the same parameter list as s4soln, the routine that
*     prints the solution.  It will be called if the SPECS file
*     specifies  REPORT file  n  for some positive value of  n.
*
*     pi contains the unscaled dual solution.
*     x contains the unscaled primal solution.  There are n + m = nb
*        values (n structural variables and m slacks, in that order).
*     y  contains the true slack values for nonlinear constraints
*        in its first nnCon components (computed by s8nslk).
*
*     This version of s4rept does nothing.    Added for PILOT, Oct 1985.
*
*     15 Nov 1991: First version based on Minos 5.4 routine m4rept.
*     26 Mar 2000: Updated for SNOPT 6.1.
*     01 Aug 2003: snPRNT adopted.
*     01 Aug 2003: Current version of s4rept.
*     ==================================================================
      character
     &     str*80
*     ------------------------------------------------------------------
      write(str, 1000)
      call snPRNT( 13, str, iw, leniw )
      return

 1000 format(' XXX Report file requested.  s4rept does nothing.')

      end ! subroutine s4rept

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s4savB
     &   ( iExit, Task, minimz, m, n, nb, nkx,
     &     nnCon0, nnCon, nnObj0, nnObj, nName, nS,
     &     itn, nInf, sInf, wtInf, vimax, iObj, sclObj, ObjTru,
     &     pNorm1, pNorm2, piNorm, xNorm,
     &     ne, nlocJ, locJ, indJ, Jcol, kx,
     &     hEstat, hs, Ascale, bl, bu, Fx, gObj,
     &     Names, pi, rc, x,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Task, iExit, minimz, itn, iObj, nInf, m, n,
     &     nb, ne, nkx, nlocJ, nnCon0, nnCon, nnObj0, nnObj, nName, nS,
     &     lencw, leniw, lenrw, kx(nkx), locJ(nlocJ),
     &     indJ(ne), hEstat(nb), hs(nb), iw(leniw)
      double precision
     &     ObjTru, pNorm1, pNorm2, piNorm, sInf, sclObj, vimax, wtInf,
     &     xNorm, Jcol(ne), Ascale(nb), bl(nb), bu(nb), Fx(nnCon0),
     &     gObj(nnObj0), rc(nb), x(nb), pi(m), rw(lenrw)
      character
     &     Names(nName)*8, cw(lencw)*8

*     ==================================================================
*     s4savB  saves basis files  and/or  prints the solution.
*     It is called twice at the end of s5solv
*              and twice at the end of s8solv.
*
*     If Task = SaveB, the problem is first unscaled, then from 0 to 4
*     files are saved (PUNCH file, DUMP file, SOLUTION file,
*     REPORT file, in that order).
*     A NEW BASIS file, if any, will already have been saved by s8SQP.
*     A call with Task = SaveB must precede a call with Task = PrintS.
*
*     If Task = PrintS, the solution is printed under the control of
*     lprSol (which is set by the Solution keyword in the SPECS file).
*
*     15 Nov 1991: First version based on Minos 5.4 routine m4savb.
*     19 Feb 1994: Use s4rc to compute reduced costs.
*     05 Apr 1996: s2rcA called to get the reduced costs (as in
*                  Minos 5.5). Maximum primal and dual infeasibilities
*                  computed and printed here.
*     14 Jul 1997: Thread-safe version.
*     26 Mar 2000: Updated for SNOPT 6.1.
*     16 Jan 2003: Reinstated calls to s4dump, s4pnch.
*     27 Jul 2003: Print max elements of x and pi (and which ones).
*                  This is more helpful than piNorm >= 1.
*     01 Aug 2003: snPRNT adopted.
*     09 Mar 2004: s4soln now deals with the UNSCALED solution,
*                  so we no longer have to save the scaled pinorm.
*     17 Jun 2004: Save biggest primal and dual infeasibilities
*                  before and after scaling, so GAMS can use them.
*     17 Jun 2004: Current version of s4savB.
*     ==================================================================
      external
     &     idamax, s2VarN
      character
     &     istate(3)*4, str*80
      logical
     &     feasbl, prnt
      integer
     &     iDump, iPnch , iPrint, iReprt,
     &     iSoln, lvlScl, lprSol, jbInf, jbInf1,
     &     jdInf, jdInf1, k, maxp, maxp1, maxvi, maxx, maxx1,
     &     idamax, s2VarN
      double precision
     &     bInf, bInf1, dInf, dInf1, eps0, infBnd, pimax, pimax1,
     &     tolx, xNorm1
*     ------------------------------------------------------------------
      integer            UnScal
      parameter         (UnScal = 1)
      integer            SaveB,      PrintS
      parameter         (SaveB  = 0, PrintS = 1)
      double precision   one
      parameter         (one    = 1.0d+0)
*     ------------------------------------------------------------------

      iPrint    = iw( 12) ! Print file
      iDump     = iw(121) ! dump file
      iPnch     = iw(127) ! punch file
      iReprt    = iw(130) ! Report file
      iSoln     = iw(131) ! Solution file

      lvlScl    = iw( 75) ! scale option
      lprSol    = iw( 84) ! > 0    =>  print the solution

      eps0      = rw(  2) ! eps**(4/5)          IEEE DP  3.00e-13
      tolx      = rw( 56) ! Minor feasibility tolerance.
      infBnd    = rw( 70) ! definition of an infinite bound

      feasbl    = nInf .eq. 0
      k         = 1 + min(9,iExit)/10
      call s4stat( k, istate )

      if (Task .eq. SaveB) then
*        ---------------------------------------------------------------
*        Compute rc and unscale Jcol, bl, bu, g, pi, x, xNorm
*        and piNorm (but s4soln uses scaled piNorm, so save it).
*        Then save basis files.
*        ---------------------------------------------------------------
*        Compute reduced costs rc(*) for all columns and rows.
*        Find the maximum bound and dual infeasibilities.

         call s2rcA
     &      ( feasbl, tolx, iObj, minimz, wtInf,
     &        m, n, nb, nnObj0, nnObj, ne, nlocJ, locJ, indJ, Jcol,
     &        hEstat, hs, bl, bu, gObj, pi, rc, x )
         call s2bInf
     &      ( nb, bl, bu, x, bInf, jbInf )
         call s2dInf
     &      ( n, nb, iObj, eps0, bl, bu, rc, x, dInf, jdInf )

         jbInf  = s2VarN( jbInf, leniw, iw )
         jdInf  = s2VarN( jdInf, leniw, iw )

         bInf1  =  bInf
         dInf1  =  dInf
         jbInf1 = jbInf
         jdInf1 = jdInf
         rw(427)=  bInf
         rw(428)=  dInf
         iw(427)= jbInf
         iw(428)= jdInf

         maxx   = idamax( n,  x, 1 )
         maxp   = idamax( m, pi, 1 )
         xNorm  = abs(  x(maxx) )
         pimax  = abs( pi(maxp) )
         piNorm = max( pimax, one )

         maxx1  = maxx
         maxp1  = maxp
         xNorm1 = xNorm
         pimax1 = pimax
         pNorm1 = piNorm

*        Unscale a, bl, bu, pi, x, rc, Fx, gObj and xNorm, piNorm.
*        (Previously, s4soln used the scaled piNorm, but no more.)

         if (lvlScl .gt. 0) then
            call s2scla
     &         ( UnScal, m, n, nb, iObj, infBnd, sclObj,
     &           ne, nlocJ, locJ, indJ, Jcol,
     &           Ascale, bl, bu, pi, x )
            call dddiv
     &         ( nb , Ascale, 1, rc, 1 )

            if (lvlScl .eq. 2) then
               if (nnCon .gt. 0)
     &         call ddscl
     &            ( nnCon, Ascale(n+1), 1, Fx  , 1 )
               if (nnObj .gt. 0)
     &         call dddiv
     &            ( nnObj, Ascale     , 1, gObj, 1 )
            end if

            maxx   = idamax( n,  x, 1 )
            maxp   = idamax( m, pi, 1 )
            xNorm  = abs(  x(maxx) )
            pimax  = abs( pi(maxp) )
            piNorm = max( pimax, one )

            call s2bInf
     &         ( nb, bl, bu, x, bInf, jbInf )
            call s2dInf
     &         ( n, nb, iObj, eps0, bl, bu, rc, x, dInf, jdInf )
            jbInf  = s2VarN( jbInf, leniw, iw )
            jdInf  = s2VarN( jdInf, leniw, iw )
         end if

         pNorm2  = piNorm
         rw(422) = piNorm ! Lagrange multiplier norm

         rw(427) =  bInf
         rw(428) =  dInf
         iw(427) = jbInf
         iw(428) = jdInf

*        ---------------------------------------------------------------
*        Print various scaled and unscaled norms.
*        xNorm1, pNorm1  are  scaled (if scaling was used)
*                pNorm2  is unscaled
*                piNorm  is unscaled
*        ---------------------------------------------------------------
         if (lvlScl .gt. 0) then
            write(str, 1010) maxx1, xNorm1, maxp1, pimax1
            call snPRNT( 3, str, iw, leniw )
         end if
            write(str, 1020) maxx , xNorm , maxp , pimax
            call snPRNT( 3, str, iw, leniw )

         if (lvlScl .gt. 0) then
            write(str, 1030) jbInf1, bInf1, jdInf1, dInf1
            call snPRNT( 3, str, iw, leniw )
         end if
            write(str, 1040) jbInf , bInf, jdInf , dInf
            call snPRNT( 3, str, iw, leniw )

*        Change the sign of pi and rc if feasible and maximizing.

         if (nInf .eq. 0  .and.  minimz .lt. 0) then
            call dscal ( m , (-one), pi, 1 )
            call dscal ( nb, (-one), rc, 1 )
         end if

*        Compute nonlinear constraint infeasibilities (violations).

         if (nnCon .gt. 0) then
            call s2vmax
     &         ( n, nnCon, maxvi, vimax, bl, bu, Fx )
            write(str, 1080) vimax
            call snPRNT( 3, str, iw, leniw )
         end if

*        Output Punch, Dump, Solution and/or Report files.

         if (ipnch .gt. 0) then
            call s4pnch
     &         ( ipnch, m, n, nb, nName,
     &           hs, bl, bu, x, names, cw, lencw, iw, leniw )
         end if

         if (idump .gt. 0) then
            call s4dump
     &         ( idump, m, n, nb, nName,
     &           hs, x, names, cw, lencw, iw, leniw )
         end if

         ! 09 Mar 2004: No longer worry about scaled solution.
         ! piNorm = pNorm1
         if (iSoln .gt. 0) then
            call s4soln
     &         ( .true., minimz, m, n, nb, nkx, nName,
     &           nnCon0, nnCon, nnObj0, nnObj, nS, iObj,
     &           itn, nInf, sInf, ObjTru, piNorm, Names,
     &           ne, nlocJ, locJ, indJ, Jcol, kx,
     &           hs, Ascale, bl, bu,
     &           gObj, pi, rc, x, Fx, istate,
     &           cw, lencw, iw, leniw, rw, lenrw )
         end if

         ! 09 Oct 2000: Export A, B or (B S) to Report file 91, 92 or 93.

         if (iReprt .ge. 91  .and.  iReprt .le. 93) then
            call s2xmat
     &         ( iReprt, n, nb,
     &           ne, nlocJ, locJ, indJ, Jcol, hs )
         else if (iReprt .gt. 0) then
            call s4rept
     &         ( .true., m, n, nb, nName,
     &           nnCon0, nnCon, nnObj0, nnObj, nS,
     &           ne, nlocJ, locJ, indJ, Jcol,
     &           hs, Ascale, bl, bu, gObj, pi, x, Fx,
     &           Names, istate, iw, leniw )
         end if

         ! 09 Mar 2004: No longer worry about scaled solution.
         ! piNorm = pNorm2

      else if (Task .eq. PrintS) then
*        ---------------------------------------------------------------
*        Print solution if requested.
*
*        lprSol = 0   means   no
*               = 1   means   if optimal, infeasible or unbounded
*               = 2   means   yes
*               = 3   means   if error condition
*        ---------------------------------------------------------------
         prnt   = iPrint .gt. 0  .and.  lprSol .gt. 0
         if ((lprSol .eq. 1  .and.  iExit .gt. 2)  .or.
     &       (lprSol .eq. 3  .and.  iExit .le. 2)) prnt = .false.

         if ( prnt ) then
            ! 09 Mar 2004: No longer worry about scaled solution.
            ! piNorm = pNorm1
            call s4soln
     &         ( .false., minimz, m, n, nb, nkx, nName,
     &           nnCon0, nnCon, nnObj0, nnObj, nS, iObj,
     &           itn, nInf, sInf, ObjTru, piNorm, Names,
     &           ne, nlocJ, locJ, indJ, Jcol, kx,
     &           hs, Ascale, bl, bu,
     &           gObj, pi, rc, x, Fx, istate,
     &           cw, lencw, iw, leniw, rw, lenrw )
            ! 09 Mar 2004: No longer worry about scaled solution.
            ! piNorm = pNorm2
            write(str, 1200) iPrint
            call snPRNT( 12, str, iw, leniw )
         else
            write(str, 1300)
            call snPRNT( 12, str, iw, leniw )
         end if
      end if

      return

 1010 format(  ' Max x       (scaled)', i9, 1p, e8.1,
     &     2x, ' Max pi      (scaled)', i9,     e8.1)
 1020 format(  ' Max x               ', i9, 1p, e8.1,
     &     2x, ' Max pi              ', i9,     e8.1)
 1030 format(  ' Max Prim inf(scaled)', i9, 1p, e8.1,
     &     2x, ' Max Dual inf(scaled)', i9,     e8.1)
 1040 format(  ' Max Primal infeas   ', i9, 1p, e8.1,
     &     2x, ' Max Dual infeas     ', i9,     e8.1)
 1080 format(  ' Nonlinear constraint violn', 1p, e11.1)

 1200 format(' Solution printed on file', i4)
 1300 format(' Solution not printed')

      end ! subroutine s4savB

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s4soln
     &   ( ondisk, minimz, m, n, nb, nkx, nName,
     &     nnCon0, nnCon, nnObj0, nnObj, nS, iObj,
     &     itn, nInf, sInf, Objtru, piNorm, Names,
     &     ne, nlocA, locA, indA, Acol, kx,
     &     hs, Ascale, bl, bu,
     &     gObj, pi, rc, x, Fx, istate,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     ondisk
      integer
     &     minimz, nInf, itn, m, n, nb, ne, nkx, nlocA, nName, nnCon0,
     &     nnCon, nnObj0, nnObj, nS, iObj, lencw, leniw, lenrw, kx(nkx),
     &     indA(ne), hs(nb), locA(nlocA), iw(leniw)
      double precision
     &     Objtru, piNorm, sInf,
     &     Acol(ne), Ascale(nb), bl(nb), bu(nb), x(nb), gObj(nnObj0),
     &     rc(nb), pi(m), Fx(nnCon0), rw(lenrw)
      character
     &     istate(3)*4, Names(nName)*8, cw(lencw)*8

*     ==================================================================
*     s4soln  is the standard output routine for printing the solution.
*
*     On entry,
*     pi    contains the dual solution.
*     x     contains the primal solution.  There are n + m = nb values
*           (n structural variables and m slacks, in that order).
*     Fx    contains the true slack values for nonlinear constraints.
*
*     All quantities a, bl, bu, pi, x, Fx, g are unscaled,
*     and adjusted in sign if maximizing.
*
*     If ondisk is true, the solution is output to the solution file.
*     Otherwise, it is output to the printer.
*
*     15 Nov 1991: First version based on Minos 5.4 routine m4soln.
*     26 Jul 1996: Slacks modified.
*     26 Mar 2000: Updated for SNOPT 6.1.
*     10 Oct 2003: snPRNT adopted.
*     09 Mar 2004: Solution flags are now set by snSolF.
*                  They are defined by the UNSCALED solution!
*                  Ascale(*) is no longer referenced.
*     05 May 2006: printing to the print file done using snPRNT.
*     ==================================================================
      character
     &     form1*4, form2*4,
     &     line*111, str*132, mProb*8, mObj*8, mRhs*8, mRng*8, mBnd*8,
     &     Objtyp(3)*8, id*8
      logical
     &     feasbl
      integer
     &     iPrint, iSoln, i, iloop, iN, ir, j, jkey, jloop, jN,
     &     jstate, k, length
      double precision
     &     infBnd, b1, b2, bplus, cj, d1, d2, dj, slk, py, row, xj
*     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero   = 0.0d+0)
      character          cdummy*8
      parameter         (cdummy = '-1111111')
      parameter         (form1  = '( a)')
      parameter         (form2  = '(/a)')
      data               Objtyp /'Max     ', 'Feas    ', 'Min     '/
*     ------------------------------------------------------------------
      iPrint    = iw( 12) ! Print file
      iSoln     = iw(131) ! Solution file

      infBnd    = rw( 70) ! definition of an infinite bound

*     lvlScl    = iw( 75) ! scale option

      mProb     = cw( 51) ! Problem name
      mObj      = cw( 52) ! Objective name
      mRhs      = cw( 53) ! Right-hand side name
      mRng      = cw( 54) ! Range name
      mBnd      = cw( 55) ! Bnd section name

      bplus     = 0.1d+0*infBnd
!!    scale     = one
      feasbl    = nInf   .eq. 0
*     maximz    = minimz .lt. 0
!!    scaled    = lvlScl .gt. 0

      call s1page( 1, iw, leniw )

      if (feasbl) then
         write(str, 1002) mProb, Objtru
         if (ondisk) then
            call s1trim( str, length )
            write(iSoln, form1) str(1:length)
         else
            call snPRNT( 1, str, iw, leniw )
         end if
      else
         write(str, 1000) mProb, nInf, sInf
         if (ondisk) then
            call s1trim( str, length )
            write(iSoln, form1) str(1:length)
         else
            call snPRNT( 1, str, iw, leniw )
         end if
      end if

      write(str, 1004) istate, itn, nS
      if (ondisk) then
         call s1trim( str, length )
         write(iSoln, form2) str(1:length)
      else
         call snPRNT( 11, str, iw, leniw )
      end if

      if (mObj .ne. cdummy) then
         write(str, 1005) mObj, Objtyp(minimz+2)(1:3)
         if (ondisk) then
            call s1trim( str, length )
            write(iSoln, form2) str(1:length)
         else
            call snPRNT( 11, str, iw, leniw )
         end if

         write(str, 1006) mRhs
         if (ondisk) then
            call s1trim( str, length )
            write(iSoln, form1) str(1:length)
         else
            call snPRNT( 1, str, iw, leniw )
         end if

         write(str, 1007) mRng
         if (ondisk) then
            call s1trim( str, length )
            write(iSoln, form1) str(1:length)
         else
            call snPRNT( 1, str, iw, leniw )
         end if

         write(str, 1008) mBnd
         if (ondisk) then
            call s1trim( str, length )
            write(iSoln, form1) str(1:length)
         else
            call snPRNT( 1, str, iw, leniw )
         end if
      end if

      if (ondisk) then
         write(iSoln, 1010)
      else
         write(str, '(a)')
     &       ' Section 1 - Rows'
         call snPRNT( 11, str, iw, leniw )
         write(str, '(a)')
     &       '  Number  ...Row.. State  ...Activity...  Slack Activity'
     &     //'  ..Lower Limit.  ..Upper Limit.  .Dual Activity    ..i'
         call snPRNT( 11, str, iw, leniw )
         call snPRNT(  1, ' ', iw, leniw )
      end if

*     ------------------------------------------------------------------
*     Output the ROWS section.
*     ------------------------------------------------------------------
      do iloop = 1, m
         iN    = iloop
         jN    = n + iN
         i     = kx(jN)
         if (i .gt. 0) then
            j      = n + i
         !! if (scaled) scale = Ascale(j)
            b1     = bl(j)
            b2     = bu(j)
            xj     = x (j)
            py     = pi(i)
            dj     = rc(j)

*           Define the row value and slack activities.
*           The slack activity is the distance of the row value to its
*           nearest bound. (For a free row, it is just the row value).

            if (i .le. nnCon) then
               xj     = Fx(i)
            end if

            row    =      xj
            d1     = b1 - xj
            d2     = xj - b2
            slk    =    - d1
            if (abs( d1  )  .gt.  abs( d2 )) slk =  d2
            if (abs( slk )  .ge.  bplus    ) slk = row

         !! SCALING NOW IGNORED.
         !! d1     =   d1 / scale
         !! d2     =   d2 / scale

            call s4id
     &         ( jN, n, nb, nName, Names, id )
            call snSolF
     &         ( m, n, nb, ninf, j, jkey, jstate,
     &           hs, bl, bu, rc, x, iw, leniw, rw, lenrw )
            call s4solp
     &         ( ondisk, line, bplus, jkey, jstate,
     &           jN, id, row, slk, b1, b2, py, iN )
            call snPRNT( 1, line, iw, leniw )
         end if
      end do

*     ------------------------------------------------------------------
*     Output the COLUMNS section.
*     ------------------------------------------------------------------
      call s1page( 1, iw, leniw )

      if (ondisk) then
         write(iSoln, 1020)
      else
         write(str, '(a)')
     &       ' Section 2 - Columns'
         call snPRNT( 1, str, iw, leniw )
         write(str, '(a)')
     &       '  Number  .Column. State  ...Activity...  .Obj Gradient.'
     &     //'  ..Lower Limit.  ..Upper Limit.  Reduced Gradnt    m+j'
         call snPRNT( 11, str, iw, leniw )
         call snPRNT(  1, ' ', iw, leniw )
      end if

      do jloop = 1, n
         jN    = jloop
         j     = kx(jN)
      !! if (scaled) scale = Ascale(j)
         b1     = bl(j)
         b2     = bu(j)
         xj     = x(j)
         cj     = zero
         dj     = rc(j)

         do k = locA(j), locA(j+1)-1
            ir = indA(k)
            if (ir .eq. iObj) cj = Acol(k)
         end do

         d1     =   (b1 - xj) !! / scale
         d2     =   (xj - b2) !! / scale
         if (feasbl) then
            if (j .le. nnObj) cj = cj + gObj(j)
         end if

         call s4id
     &      ( jN, n, nb, nName, Names, id )
         call snSolF
     &      ( m, n, nb, ninf, j, jkey, jstate,
     &        hs, bl, bu, rc, x, iw, leniw, rw, lenrw )
         call s4solp
     &      ( ondisk, line, bplus, jkey, jstate,
     &        jN, id, xj, cj, b1, b2, dj, m+jN )
         call snPRNT( 1, line, iw, leniw )
      end do

      if (ondisk) then
         if (iSoln  .ne. iPrint) rewind iSoln
         write(str, 1400) iSoln
         call snPRNT( 13, str, iw, leniw )
      end if
      return

 1000 format(
     &   ' Name', 11x, a8, 16x,
     &   ' Infeasibilities', i7, 1p, e16.4)
 1002 format(
     &   ' Name', 11x, a8, 16x,
     &   ' Objective Value', 1p, e22.10)
 1004 format(
     &   ' Status', 9x, 3a4, 12x,
     &   ' Iteration', i7, '    Superbasics', i6)
 1005 format(
     &   ' Objective', 6x, a8, ' (', a3, ')')
 1006 format(
     &   ' RHS',      12x, a8)
 1007 format(
     &   ' Ranges',    9x, a8)
 1008 format(
     &   ' Bounds',    9x, a8)
 1010 format(/ ' Section 1 - Rows' //
     &   '  Number  ...Row.. State  ...Activity...  Slack Activity',
     &   '  ..Lower Limit.  ..Upper Limit.  .Dual Activity    ..i' /)
 1020 format(' Section 2 - Columns' //
     &       '  Number  .Column. State  ...Activity...  .Obj Gradient.',
     &       '  ..Lower Limit.  ..Upper Limit.  Reduced Gradnt    m+j'/)
 1400 format(' SOLUTION file saved on file', i4)

      end ! subroutine s4soln

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s4solp
     &   ( ondisk, line, bplus, jkey, jstate,
     &     j, id, xj, cj, b1, b2, dj, k )

      implicit
     &     none
      character*(*)
     &     line
      character
     &     id*8
      logical
     &     ondisk
      integer
     &     j, jkey, jstate, k
      double precision
     &     bplus, xj, cj, b1, b2, dj

*     ==================================================================
*     s4solp  prints one line of the Solution file.
*
*     The following conditions are marked by key:
*
*        D  degenerate basic or superbasic variable.
*        I  infeasible basic or superbasic variable.
*        A  alternative optimum      (degenerate nonbasic dual).
*        N  nonoptimal superbasic or nonbasic (infeasible dual).
*
*     Prior to 09 Mar 2004,
*     tests for these conditions were performed on scaled quantities
*     d1, d2, djtest,
*     since the correct indication was then more likely to be given.
*     On badly scaled problems, the unscaled solution could then appear
*     to be flagged incorrectly, but it would be just an "illusion".
*
*     15 Nov 1991: First version based on Minos 5.4 routine m4solp.
*     18 Oct 1993: Replaced by modified Minos 5.4 routine m4solp.
*                  Infinite bounds and certain other values treated
*                  specially.
*     10 Oct 2003: snEXIT and snPRNT adopted.
*     09 Mar 2004: Now use jkey and jstate from snSolF.
*                  d1, d2, djtest are no longer used in here.
*     09 Mar 2004: Large xj, cj, dj (as well as b1, b2) use e format.
*     04 Dec 2004: Current version of s4solp.
*     ==================================================================
      character
     &     ckey(0:4)*1, cstate(0:5)*4
      character
     &     lzero*16, lone*16, lmone*16, none*16
      character
     &     e*10, form*82
*     ------------------------------------------------------------------
      double precision   zero,           one,           big
      parameter        ( zero = 0.0d+0,  one = 1.0d+0,  big = 1.0d+9 )

      data               ckey   /' ', 'A', 'D', 'I', 'N'/
      data               cstate /' LL ', ' UL ', 'SBS ',
     &                           ' BS' , ' EQ' , ' FR '/
      data               lzero  /'          .     '/
      data               lone   /'         1.0    '/
      data               lmone  /'        -1.0    '/
      data               none   /'           None '/
      data               e      /' 1p,e16.6,'/
*     ------------------------------------------------------------------

      ! Select format for printing.

      if (ondisk) then
         write(line, 1000) j, id, ckey(jkey), cstate(jstate),
     &                     xj, cj, b1, b2, dj, k
      else
         form = '(i8, 2x,  a8, 1x, a1, 1x,a3, ' //
     &          '0p,f16.5, 0p,f16.5, 0p,f16.5, 0p,f16.5, 0p,f16.5, i7)'
         if (abs( xj ) .ge. big) form(29:38) = e
         if (abs( cj ) .ge. big) form(39:48) = e
         if (abs( b1 ) .ge. big) form(49:58) = e
         if (abs( b2 ) .ge. big) form(59:68) = e
         if (abs( dj ) .ge. big) form(69:78) = e
         write(line, form) j, id, ckey(jkey), cstate(jstate),
     &                     xj, cj, b1, b2, dj, k
      end if

*     Test for 0.0, 1.0 and -1.0

      if (xj .eq. zero) line(25:40) = lzero
      if (xj .eq.  one) line(25:40) = lone
      if (xj .eq. -one) line(25:40) = lmone
      if (cj .eq. zero) line(41:56) = lzero
      if (cj .eq.  one) line(41:56) = lone
      if (cj .eq. -one) line(41:56) = lmone
      if (b1 .eq. zero) line(57:72) = lzero
      if (b1 .eq.  one) line(57:72) = lone
      if (b1 .eq. -one) line(57:72) = lmone
      if (b2 .eq. zero) line(73:88) = lzero
      if (b2 .eq.  one) line(73:88) = lone
      if (b2 .eq. -one) line(73:88) = lmone
      if (dj .eq. zero) line(89:104)= lzero
      if (dj .eq.  one) line(89:104)= lone
      if (dj .eq. -one) line(89:104)= lmone

      if ( ondisk ) then
         ! Relax
      else
         if (b1 .lt. -bplus) line(57:72) = none
         if (b2 .gt.  bplus) line(73:88) = none
      end if

      return

 1000 format(i8, 2x, a8, 1x, a1, 1x, a3, 1p, 5e16.6, i7)
 2000 format(a)

      end ! subroutine s4solp

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s4stat
     &   ( k, istate )

      implicit
     &     none
      integer
     &     k
      character
     &     istate(3)*4

*     ==================================================================
*     s4stat loads istate(*) with words describing the current state.
*
*     15 Nov 1991: First version based on Minos 5.4 routine m4stat.
*     20 Apr 1999: Current version.
*     ==================================================================
      integer
     &     i, j
      character
     &     c(18)*4
      data
     &     c /'Proc', 'eedi', 'ng  ',
     &        'Opti', 'mal ', 'Soln',
     &        'Infe', 'asib', 'le  ',
     &        'Unbo', 'unde', 'd   ',
     &        'Exce', 'ss i', 'tns ',
     &        'Erro', 'r co', 'ndn '/
*     ------------------------------------------------------------------
      j    = 3*min( k, 5 )
      do i = 1, 3
         istate(i) = c(j+i)
      end do

      end ! subroutine s4stat
