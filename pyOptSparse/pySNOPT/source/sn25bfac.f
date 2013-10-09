*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  sn25bfac.f
*
*     s2Bfac   s2Bkbs   s2Bmap   s2newB   s2BLU    s2Bmod   s2Bmod2
*     s2Bsol   s2sb     s2sing   s2tols
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2Bfac
     &   ( iExit, typeLU, needLU, newLU, newB,
     &     iObj, itn, lPrint, LUreq,
     &     m, mBS, n, nb, nnL, nS, nSwap,
     &     ne, nlocA, locA, indA, Acol,
     &     kBS, hs, bl, bu, blBS, buBS,
     &     nrhs0, nrhs, rhs, x, xBS,
     &     iy, iy1, y, y1, iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     needLU, newLU, newB
      integer
     &     typeLU, iExit, iObj, itn, nrhs0, nrhs, lPrint, LUreq, nnL,
     &     m, mBS, n, nb, ne, nlocA, nS, nSwap, leniw, lenrw,
     &     locA(nlocA), indA(ne), hs(nb), kBS(mBS), iy(nb), iy1(nb),
     &     iw(leniw)
      double precision
     &     Acol(ne), rhs(nrhs0), bl(nb), bu(nb), x(nb),
     &     blBS(mBS), buBS(mBS), xBS(mBS), y(nb), y1(nb), rw(lenrw)

!     =================================================================
!     s2Bfac  computes an acceptable x such that  ( A  -I )*x = b.
!     The LU factorization of the basis is computed if necessary.
!
!     If typeLU = B , the usual B = LU is computed.
!     If typeLU = BR, we want to check rank(B) with a special B = LU
!                     (TRP, tight tols) before computing a normal LU.
!     If typeLU = BS, there are some superbasics and we want to
!                     choose a good basis from the columns of (B S).
!                     We first factorize (B S)' to obtain a new B.
!                     Then B = LU is computed as usual.
!     If typeLU = BT, we should TRY 'B ' first and go back to 'BS'
!                     only if B seems ill-conditioned.
!
!     15 Nov 1991: First version based on Minos routine m2bfac.
!     29 Oct 1993: typeLU options implemented.
!                  nSwap returns the number of (B S) changes.
!     22 Apr 1994: Retry with reduced LU Factor tol
!                  if s2BLU says there was large growth in U.
!     22 Apr 1994: 'BT' option implemented to save R more often each
!                  major iteration.
!     02 Apr 1996: kObj added to mark position of Obj slack in B.
!     14 Jul 1997: Thread-safe version.
!     10 Mar 2001: BR option implemented to help with CUTE prob TRAINH.
!     10 May 2001: BTfac, BSfac or Bfac tried first.
!                  BRfac used only if growth detected in U or x.
!     17 Jun 2001: If BSfac is singular, ask for BRfac.
!     18 Jun 2001: If B fac is singular, ask for BRfac.
!     18 Nov 2001: lPrint added as parameter to s2sing.
!     31 Jul 2003: snEXIT and snPRNT adopted.
!     11 Apr 2004: kBS re-used to reduce the value of nswap
!     20 Dec 2004: If B continues to be singular after BRdone,
!                  tighten the LU tols each time.
!     21 Dec 2004: Insert slacks before tightening tols (for next LU).
!     02 Jul 2005: Calls to s2kbs and s2bs added.
!     03 Jul 2005: Current version of s2Bfac.
!     =================================================================
      character
     &     str*72
      logical
     &     badx, bigU, BRfac, BSfac, BTfac, BRdone, BSdone,
     &     needx, NewTol, prnt10, singlr
      integer
     &     Btask , Ptask , ntask, inform, j, k, kObj, maxrw, maxiw,
     &     lprDbg, lenaLU, maxLUi, maxLUr, iP, iQ, locr, LUa,
     &     indc, indr, lin, ntry, nBS, nBasic, nslack,
     &     nonlin, more, newi, newr, dUmin, Umin, nSing, minlen,
     &     nFac, nBfac, LUitn, LUmod
      double precision
     &     eps2, rowerr, Utol
!     ------------------------------------------------------------------
      integer            RedTol
      parameter         (RedTol = 1)
      integer            Reset
      parameter         (Reset  = 0)
      integer            B,          BR,         BS        , BT
      parameter         (B      = 0, BR     = 1, BS     = 2, BT     = 3)
      integer            xFirst,     sFirst
      parameter         (xFirst = 0, sFirst = 1)
      integer            mtry
      parameter         (mtry   = 40)
!     LUSOL arguments.
      parameter         (nSing     = 161) ! # of singularities in w(*)
      parameter         (minlen    = 163) ! minimum recommended lenaLU
      parameter         (dUmin     = 164) ! minimum diagonal in  U.
!
      parameter         (Umin      = 190) ! saved smallest U diagonal
      parameter         (kObj      = 205) ! xBS(kObj) is the obj. slack
      parameter         (nFac      = 210) ! # of LU factorizations
      parameter         (nBFac     = 211) ! # consecutive `B' facts
      parameter         (LUitn     = 215) ! itns since last factorize
      parameter         (LUmod     = 216) ! number of LU mods

      double precision   zero
      parameter         (zero      = 0.0d+0)
!     ------------------------------------------------------------------
      character         LUtype(0:3)*2
      data              LUtype/'B ', 'BR', 'BS', 'BT'/
!     ------------------------------------------------------------------
      eps2      = rw(  4) ! eps**(1/2)       IEEE DP  1.49e-08
      maxrw     = iw(  3) ! end of SNOPT part of rw
      maxiw     = iw(  5) ! end of SNOPT part of iw
      lprDbg    = iw( 85) ! > 0    => private debug print
      lenaLU    = iw(213) ! space allotted for LU factors
      maxLUi    = iw(361) ! max LU nonzeros in iw(*)
      maxLUr    = iw(362) ! max LU nonzeros in rw(*)
      iP        = iw(363) !
      iQ        = iw(364) !
      locr      = iw(368) !
      LUa       = iw(371) !
      indc      = iw(373) !
      indr      = iw(374) !

      newB      = .false.
      BRdone    = .false.
      BSdone    = .false.
      BRfac     = .false.
      BSfac     = .false.
      BTfac     = .false.

      nBS       = m + nS
      nSwap     = 0
      ntask     = 0
      Ptask     = typeLU
      prnt10    = lPrint .ge. 10

!     Initialize Umin and nBFac on first entry.
!     nBFac  counts consecutive B factorizations (reset if BS is done).
!     Umin   is the smallest diagonal of U after last BS factor.

      if (iw(nFac) .eq. 0) then
         rw(Umin)  = zero
         iw(nBFac) = 0
      end if

      if (needLU) then
         if (iw(nFac) .ge. 1) then
            iw(nFac)  = iw(nFac)  + 1
         else
            iw(nFac)  = iw(nFac)  + 1
         end if
         iw(nBFac) = iw(nBFac) + 1

         if ( prnt10 ) then
            write(str, 1005) iw(nFac), LUreq, itn
            call snPRNT( 11, str, iw, leniw )
         end if
         iw(LUitn) = 0
         iw(LUmod) = 0
         LUreq     = 0

!        ---------------------------------------------------------------
!        Set local logicals to select the required type of LU.
!        We come back to 100 if a BT factorize looks doubtful.
!        If BT was requested but we haven't done BS yet,
!        might as well do BS now.
!        ---------------------------------------------------------------
         BTfac  =  typeLU .eq. BT     .and.  nS       .gt. 0
         BSfac  = (typeLU .eq. BS     .and.  nS       .gt. 0   )  .or.
     &            (            BTfac  .and.  rw(Umin) .eq. zero)
         BRfac  =  typeLU .eq. BS     .and.  nS       .eq. 0
         if (.not. (BTfac .or. BSfac)) Ptask = B
      end if ! needLU

!     ------------------------------------------------------------------
!     We may come back here to do a BSfac after all.
!     ------------------------------------------------------------------
  100 iExit  = 0
      singlr = .false.

      if ( BSfac ) then
!        ---------------------------------------------------------------
!        Repartition (B S) to get a better B.
!        ---------------------------------------------------------------
         BTfac     = .false.
         iw(nBFac) = 1
         Ptask     = BS
                                ! Load the basics,  x first, into kBS.
         call s2Bkbs
     &      ( xFirst, iObj, m, mBS, n, nb, nnL, hs, kBS,
     &        iw(kObj), nBasic, nonlin, nSlack )

         if (nBasic .eq. m) then
!           ------------------------------------------------------------
!           We have the right number of basics.
!           1. Factorize (B S)'.
!           2. Apply the resulting row permutation to the columns
!              of (B S).
!           ------------------------------------------------------------
            ntask  = ntask + 1
            if ( prnt10 ) then
               if (ntask .gt. 1) call snPRNT( 1, ' ', iw, leniw )
               write(str, 2000) itn, LUtype(Ptask)
               call snPRNT( 1, str, iw, leniw )
            end if

            call s2BLU
     &         ( inform, BS, lPrint, m, n, nb, nBS,
     &           ne, nlocA, locA, indA, Acol,
     &           kBS, iw(iP), rw(LUa), iw(indc),iw(indr), lenaLU,
     &           iy, iy1, y, iw, leniw, rw, lenrw )

            if (inform .ge. 7) then
               iExit =  80      ! insufficient storage for the LU.
               go to 400
            else if (inform .ge. 3) then
               iExit = 142      ! Error in basis package
               go to 400
            end if

            singlr = iw(nSing) .gt. 0
            BSdone = .true.     ! Once only!
            Ptask  = B

            call s2newB
     &         ( m, mBS, nb, nS, hs, iw(iP), kBS, iy, iw(locr), nSwap )

            if (nSwap .gt. 0) newB = .true.
            if ( prnt10 ) then
               write(str, 3000) nSwap
               call snPRNT( 1, str, iw, leniw )
            end if

            if ( singlr  .and. .not. BRdone) then
               BRfac  = .true.
               Ptask  = BR
            end if
         end if ! m basics
      end if ! BS

      needx  = .true.
      NewTol = .true.

!+    while (needx  .and.  NewTol) do
  150 if    (needx  .and.  NewTol) then

!        ===============================================================
!        Main loop to obtain a good x.
!        typeLU is not used.  (We are always factoring just B.)
!        ===============================================================
         ntry   = 0
         bigU   = .false.

!+       while (needLU  .and.  ntry .le. mtry  .and.  .not. bigU) do
  200    if    (needLU  .and.  ntry .le. mtry  .and.  .not. bigU) then
!           ------------------------------------------------------------
!           Main loop to find a good  B = LU.
!           ------------------------------------------------------------
                                ! Load and count the basics in kBS, with
                                ! the slacks loaded first.  kObj keeps
                                ! track of the linear objective.
            call s2Bkbs
     &         ( sFirst, iObj, m, mBS, n, nb, nnL, hs, kBS,
     &           iw(kObj), nBasic, nonlin, nSlack )
            if (nBasic .gt. m) then
               iExit = 141      ! Too many basic variables
               go to 400
            end if

            if (nBasic .lt. m) then
                                ! Too few basics.
                                ! Set remaining kBS(k) = 0 for s2sing.
               call iload
     &            ( m-nBasic, 0, kBS(nBasic+1), 1 )
            end if

            lin    = nBasic - nSlack - nonlin
            if (lin .lt. 0) lin = 0

!           ------------------------------------------------------------
!           Load the basis matrix into the LU arrays and factorize it.
!           ------------------------------------------------------------
            ntry   = ntry  + 1
            ntask  = ntask + 1

            if (BRfac) then
               Btask  = BR
            else
               Btask  = B
            end if

            if ( prnt10 ) then
               if (ntask .gt. 1)
     &         call snPRNT( 1, ' ', iw, leniw )
               write(str, 1100) nonlin, lin, nSlack, LUtype(Ptask)
               call snPRNT( 1, str, iw, leniw )
            end if

            if (Ptask .eq. BR) then
               write(str, 2000) itn, LUtype(Ptask)
               call snPRNT( 23, str, iw, leniw )
            end if

            call s2BLU
     &         ( inform, Btask, lPrint, m, n, nb, nBS,
     &           ne, nlocA, locA, indA, Acol,
     &           kBS, iw(iP), rw(LUa), iw(indc), iw(indr), lenaLU,
     &           iy, iy1, y, iw, leniw, rw, lenrw )

            needLU = inform    .gt. 0  .or.  BRfac
            bigU   = inform    .eq. 2
            singlr = iw(nSing) .gt. 0

            if (     inform .ge. 7) then
               iExit =  80      ! Insufficient storage for the LU.
            else if (inform .ge. 3) then
               iExit = 142      ! error in basis package
            else if (inform .gt. 0  .and.  ntry .gt. mtry) then
               iExit =  42      ! Singular basis
            end if

            if (iExit .gt. 0) go to 400

            if ( BRfac ) then
               BRfac  = .false.
               BRdone = .true.
               Ptask  = B
            end if

            if ( BSfac ) then
!              --------------------------------------------------------
!              We started with a BS factorize this time.
!              Save the smallest diag of U.
!              --------------------------------------------------------
               rw(Umin) = rw(dUmin)

            else if ( BTfac ) then
!              --------------------------------------------------------
!              (We come here only once.)
!              See if we should have done a BS factorize after all.
!              In this version we do it if any of the following hold:
!              1. dUmin (the smallest diagonal of U) is noticeably
!                 smaller than Umin (its value at the last BS factor).
!              2. dUmin is pretty small anyway.
!              3. B was singular.
!              nBFac  makes BS increasingly likely the longer we
!              keep doing B and not BS.
!              --------------------------------------------------------
               BTfac  = .false.
               Utol   = rw(Umin) * 0.1d+0 * iw(nBFac)
               BSfac  = rw(dUmin) .le. Utol   .or.
     &                  rw(dUmin) .le. eps2   .or.   singlr
               if ( BSfac ) then
                  needLU = .true.
                  Ptask  = BS
                  go to 100
               end if
            else if (nS .eq. 0) then
               rw(Umin) = rw(dUmin)
            end if ! BS

            !----------------------------------------------------------
            ! Deal with singularity.
            !----------------------------------------------------------
            if ( singlr ) then
               if (.not. BRdone) then
                  BRfac  = .true.
                  Ptask  = BR
                  go to 150
               end if

               ! Suspect columns are indicated by y(j) <= 0.
               ! Replace them by suitable slacks.
               ! Then check if any superbasic slacks were made basic.

               call s2sing
     &            ( lPrint, mBS, m, n, nb,
     &              y, iw(iP), iw(iQ), bl, bu,
     &              hs, kBS, x, iw, leniw )

               if (nS .gt. 0)
     &         call s2sb
     &            ( m, mBS, nb, nS, nBS, hs, kBS, xBS )

               ! Tighten tols for next LU if possible.

               call s2tols
     &            ( RedTol, NewTol, itn, iw, leniw, rw, lenrw )

               newB   = .true.
               Ptask  = B
            end if ! singlr
            go to 200
         end if
!+       end while needLU

!        ---------------------------------------------------------------
!        We have a nonsingular B such that B = LU.
!        Compute the basic variables and check that  (A -I)*x = b.
!        s5setx also loads the basic/superbasic variables in xBS.
!        If the row check fails (or U was big earlier), request BRfac.
!        ---------------------------------------------------------------
         call s5setx
     &      ( inform, Reset, itn,
     &        m, n, nb, nBS, rowerr,
     &        ne, nlocA, locA, indA, Acol,
     &        kBS, xBS, nrhs0, nrhs, rhs, x, y, y1,
     &        iw, leniw, rw, lenrw )
         badx  = inform .gt. 0
         needx = bigU   .or. badx

         if (needx) then
            needLU = .true.

            if (.not. BSdone) then
               BSfac = nS .gt. 0
               if ( BSfac ) then
                  Ptask = BS
                  go to 100
               end if
            end if

            if (.not. BRdone) then
               BRfac  = .true.
               Ptask  = BR
            else
               call s2tols
     &            ( RedTol, NewTol, itn, iw, leniw, rw, lenrw )
               BRfac  = .true.
               BRdone = .false.
               Ptask  = BR
            end if
         end if

         go to 150
      end if
!+    end while needx and NewTol

      if (needx) iExit = 43  ! Cannot satisfy the linearized constraints

!     ==================================================================
!     Tidy up
!     ------------------------------------------------------------------
  400 if (iExit .eq. 0) then
!        --------------------------------------------------------------
!        Normal exit.
!        Load the basic/superbasic bounds into blBS, buBS.
!        --------------------------------------------------------------
         do k = 1, nBS
            j       = kBS(k)
            blBS(k) = bl(j)
            buBS(k) = bu(j)
         end do

         newLU = ntry .gt. 0

         if (lprDbg .eq. 100) then
            call snPRNT( 11, ' BS and SB values:', iw, leniw )
            do k = 1, nBS
               write(str, 7000) kBS(k), xBS(k)
               call snPRNT( 1, str, iw, leniw )
            end do
         end if
      else
!        --------------------------------------------------------------
!        Error exits (all fatal)
!        Some values of iExit invoke additional output.
!        --------------------------------------------------------------
         if (iExit .eq. 42) then
!           ---------------------------------------------
!           The basis is singular after mtry tries.
!           Time to give up.
!           ---------------------------------------------
            write(str, 9032) mtry
            call snPRNT( 14, str, iw, leniw )

         else if (iExit .eq. 80) then
!           ---------------------------------------------
!           Insufficient storage to factorize B.
!           ---------------------------------------------
            more   = iw(minlen) -   lenaLU
            newi   = maxiw      + 2*more
            newr   = maxrw      +   more
            if (maxLUi .lt. maxLUr) then
               iExit = 83       ! not enough integer storage for LU
            else
               iExit = 84       ! not enough real storage for LU
            end if

            write(str, 9080)
            call snPRNT( 13, str, iw, leniw )
            write(str, 9081) maxiw, newi
            call snPRNT(  3, str, iw, leniw )
            write(str, 9082) maxrw, newr
            call snPRNT(  3, str, iw, leniw )

         else if (iExit .eq. 141) then
!           ---------------------------------------------
!           Wrong number of basics.
!           ---------------------------------------------
            write(str, 9101) nBasic
            call snPRNT( 14, str, iw, leniw )
         end if
      end if

      return

 1005 format(' Factor', i7, '  Demand', i7, '  Itn', i11)
 1100 format(' Nonlin', i7, '  Linear', i7, '  Slacks', i8,
     &         2x, a, ' factorize')
 2000 format(' Itn', i7, ': ', a, ' factorize')
 3000 format(' BS Factorize.   nSwap = ', i6 )
 7000 format(i7, g17.8)

 9032 format(' XXX Singular basis after ', i5,
     &       ' factorization attempts')
 9080 format(24x, '        Current    Recommended')
 9081 format(' Total integer workspace', 2i15)
 9082 format(' Total real    workspace', 2i15)
 9101 format(' XXX Wrong no. of basic variables:', i8)

      end ! subroutine s2Bfac

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2Bkbs
     &   ( start, iObj, m, mBS, n, nb, nnL, hs, kBS,
     &     kObj, nBasic, nonlin, nSlack )

      implicit
     &     none
      integer
     &     start, iObj, kObj, m, mBS, n, nb, nBasic, nnL, nonlin,
     &     nSlack, hs(nb), kBS(mBS)

!     =================================================================
!     s2Bkbs  loads the basic variables into  kBS.
!
!     11 Apr 2004: First version of   s2Bkbs.
!     01 Jul 2005: Added option to load basics from the front.
!     03 Jul 2005: Current version of s2Bkbs.
!     =================================================================
      integer
     &     j, jObj, k
!     -----------------------------------------------------------------
      integer            xFirst,     sFirst
      parameter         (xFirst = 0, sFirst = 1)
!     -----------------------------------------------------------------
      k = 0                     ! Counts the number of nonbasics

      if      (start .eq. xFirst) then
                                ! Load basics, x first.
         do j  = 1, nb
            if (hs(j) .eq. 3) then
               k      = k + 1
               kBS(k) = j
            end if
         end do

      else if (start .eq. sFirst) then
                                ! Load basics, slacks first.
                                ! kObj points to the linear objective
         jObj   = n + iObj
         nonlin = 0
         kObj   = 0

         do j = n+1, nb
            if (hs(j) .eq. 3) then
               k      = k + 1
               kBS(k) = j
               if (j  .eq. jObj) kObj = k
            end if
         end do

         nSlack  = k

         do j = 1, n
            if (hs(j) .eq. 3) then
               k = k + 1
               if (k .le. m) then
                  kBS(k) = j
                  if (j .le. nnL) nonlin = nonlin + 1
               end if
            end if
         end do
      end if

      nBasic = k

      end ! subroutine s2Bkbs

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2Bmap
     &   ( m, n, ne, maxS,
     &     nextiw, nextrw, maxiw, maxrw, liwEst, lrwEst, iw, leniw )

      implicit
     &     none
      integer
     &     m, maxS, nextiw, nextrw, maxiw, maxrw, liwEst, lrwEst,
     &     leniw, n, ne, iw(leniw)

!     ==================================================================
!     s2Bmap sets up the core allocation for the basis factors.
!
!     Normally the storage is for B = LU.
!     For nonlinear problems, we may also need to factorize (B S)' = LU,
!     where S has at most maxS columns.
!
!     On entry:
!     nextiw, nextrw  say where the LU arrays can start in iw(*), rw(*).
!     maxiw, maxrw  say where the LU arrays must end  in iw(*), rw(*).
!
!     On exit:
!     liwEst, lrwEst  estimate the minimum length of iw(*), rw(*).
!     The LU routines will have some room to maneuver if the arrays
!     are at least that long.  (Later LU messages may ask for more.)
!     ------------------------------------------------------------------
!
!     15 Nov 1991: First version based on Minos 5.4 routine m2bmap.
!     08 Nov 1993: Generalized to allow room for (B S)'.
!     11 Nov 1994: rw(*) replaced by iw(*) and rw(*).
!     14 Jul 1997: Thread-safe version.
!     17 Mar 1998: nextiw now points to start of integer LU workspace.
!     07 Mar 2002: Fixed storage allocation (minA) to allow for BS factors.
!     11 Sep 2006: Anders Goran reports overflow on QP problem with
!                  m = 20736, n = 28, ne = 242752.
!                  Changed   minA   = 6 * mBS * necolA
!                  to        minA   = 6 * min(mBS,n ) * necolA
!                  and ensured maxLUi, maxLUr are nonnegative.
!     18 May 2007: Fixed bug that saved different values for maxLUi
!                  and  maxLUr.
!     ==================================================================
      integer
     &     mBS, mLU, nLU, iP, iQ, lenc, lenri, locc, locr,
     &     iPloc, iQloc, lastiw, maxLUi, maxLUr, indc, indr, LUa,
     &     lastrw, lenaLU, necolA, minA
!     ------------------------------------------------------------------
!     Allocate arrays for an  mLU x nLU  matrix.

      mBS    = m + maxS
      mLU    = mBS
      nLU    = m

!     LU integer workspace is  iw(iP:maxiw).
!     nextiw points to the start of indc(*), indr(*).
!     indc and indr are made as long as possible.

      iP     = nextiw
      iQ     = iP     + mLU
      lenc   = iQ     + nLU
      lenri  = lenc   + nLU
      locc   = lenri  + mLU
      locr   = locc   + nLU
      iPloc  = locr   + mLU
      iQloc  = iPloc  + nLU
      lastiw = iQloc  + mLU
      nextiw = lastiw

      maxLUi = (maxiw - lastiw) / 2
      maxLUi = max( maxLUi, 0 )
      indc   = lastiw
      indr   = indc   + maxLUi

!     LU real workspace is  rw(LUa:maxrw)
!     nextrw points to the start of A(*).

      LUa    = nextrw
      lastrw = nextrw
      maxLUr = maxrw  - lastrw
      maxLUr = max( maxLUr, 0 )

!     LUSOL thinks indc(*), indr(*) and A(*) are all of length lenaLU.

      lenaLU = min( maxLUi, maxLUr )

!     Estimate the number of nonzeros in the basis factorization.
!     necolA = estimate of nonzeros per column of  A.
!     We guess that the density of the basis factorization is
!     5 times as great, and then allow 1 more such lot for elbow room.

      necolA = ne / n
      necolA = max( necolA, 5 )
!!!   minA   = 6 * min( m, n ) * necolA  ! Too little for BSfac
      minA   = 6 * min(mBS,n ) * necolA  ! mBS is biggest LU dimension
      minA   = minA   + 10000            ! So tiny problems have plenty
      liwEst = lastiw + 2*minA
      lrwEst = lastrw +   minA

      iw(213) = lenaLU    ! space allotted for LU factors
      iw(361) = lenaLU    ! max LU nonzeros in iw(*)
      iw(362) = lenaLU    ! max LU nonzeros in rw(*)
      iw(363) = iP        !
      iw(364) = iQ        !
      iw(365) = lenc      !
      iw(366) = lenri     !
      iw(367) = locc      !
      iw(368) = locr      !
      iw(369) = iPloc     !
      iw(370) = iQloc     !
      iw(371) = LUa       !
      iw(373) = indc      !
      iw(374) = indr      !

      end ! subroutine s2Bmap

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2newB
     &   ( m, mBS, nb, nS, hs, iP, kBS, kBSold, locr, nSwap )

      implicit
     &     none
      integer
     &     mBS, m, nb, nS, nSwap,
     &     hs(nb), iP(mBS), kBS(mBS), kBSold(mBS), locr(mBS)

!     ==================================================================
!     s2newB  permutes kBS(*) to reflect the permutation (B S)P,
!     where P is in iP(*).  It updates hs(*) accordingly.
!     kBSold(*) and locr(*) are needed for workspace.
!
!     30 Oct 1993: First version.
!     04 Nov 1993: kBSold, nSwap used to save old R if there's no
!                  change in the set of superbasics.
!     16 Sep 2000: Superbasic slacks counted.
!     02 Aug 2003: Superbasic slacks allowed.
!     02 Aug 2003: Current version of s2newB.
!     ==================================================================
      integer
     &     m1, nBS, i, j, k
!     ------------------------------------------------------------------
      nSwap  = 0
      m1     = m + 1
      nBS    = m + nS
      call icopy ( nBS, kBS    , 1, locr      , 1 )
      call icopy ( nS , kBS(m1), 1, kBSold(m1), 1 )

      do k = 1, nBS
         i = iP(k)
         j = locr(i)
         kBS(k) = j
         if (k .le. m) then
            hs(j) = 3
         else
            if (hs(j) .ne. 2) nSwap = nSwap + 1
            hs(j) = 2
         end if
      end do

!     Restore the old S ordering if S contains the same variables.

      if (nSwap .eq. 0) then
         call icopy ( nS, kBSold(m1), 1, kBS(m1), 1 )
      end if

      end ! subroutine s2newB

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2BLU
     &   ( iExit, Task, lPrint, m, n, nb, nBS,
     &     ne, nlocA, locA, indA, Acol,
     &     kBS, iP, aLU, indc, indr, lenaLU,
     &     iy, iy1, y, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Task, iExit, lPrint, m, n, nb, nBS, ne, nlocA, lenaLU,
     &     leniw, lenrw, locA(nlocA), indA(ne), kBS(nBS), iP(nBS),
     &     indc(lenaLU), indr(lenaLU), iy(nb), iy1(nb), iw(leniw)
      double precision
     &     Acol(ne), aLU(lenaLU), y(nb), rw(lenrw)

!     ==================================================================
!     s2BLU  factorizes the basis.
!
!     Task = B   Extract basis from the constraint matrix
!                and factorize B = L U with current tols.
!
!     Task = BR  Factorize B = L U to check its rank.
!                Use TCP with tight tols and don't save L and U.
!
!     Task = BS  Factorize transpose of (B S), so that  (B') = L U,
!                                                       (S')
!                to get a new partition of (B S).
!                Use TCP with tight tols and don't save L and U.
!
!     The following tolerances are used...
!
!     luparm(3) = maxcol   lu1fac: Maximum number of columns
!                          searched allowed in a Markowitz-type
!                          search for the next pivot element.
!     luparm(6) = TPivot   TPivot = 0 means threshold partial  pivoting.
!                                   0 means threshold complete pivoting.
!     luparm(8) = keepLU   keepLU = 1 means keep L and U,
!                                   0 means discard them.
!     parmlu(1) = Lmax1  = Maximum multiplier allowed in  L  during
!                          refactorization.
!     parmlu(2) = Lmax2  = Maximum multiplier allowed during updates.
!     parmlu(3) = small  = Minimum element kept in  B  or in
!                          transformed matrix during elimination.
!     parmlu(4) = Utol1  = Abs tol for flagging small diagonals of  U.
!     parmlu(5) = Utol2  = Rel tol for flagging small diagonals of  U.
!     parmlu(6) = Uspace = Factor allowing waste space in row/col lists.
!     parmlu(7) = dens1    The density at which the Markowitz
!                          strategy should search maxcol columns
!                          and no rows.
!     parmlu(8) = dens2    The density at which the Markowitz
!                          strategy should search only 1 column.
!                          (In one version of lu1fac, the remaining
!                          matrix is treated as dense if there is
!                          sufficient storage.)
!
!     On exit,
!     iExit = 2  if there was excessive growth in U.
!                Other iExit values may be set by lu1fac:
!     iExit = 0  if the LU factors were computed.
!           = 1  if there are singularities (nSing gt 0).
!           = 2  if there was large growth in U.
!           = 3  if the matrix B has an illegal row or column index.
!           = 4  if an entry of B has the same indices as an earlier entry.
!           = 7  if insufficient storage for the LU.
!                minlen is an estimate of the necessary value of  lenaLU.
!           = 8  if there is a fatal error in lu1fac.
!
!     20 Oct 1990  Initial version based on Minos routine m2bsol.
!     07 Nov 1993: Add option to factorize (B S)'
!     06 Mar 1994: Include all rows of (B S), even if B contains slacks.
!     22 Apr 1994: Test for excessive growth in U.
!     14 Jul 1997: Thread-safe version.
!     23 Sep 2000: LU statistics now printed by LUSOL.
!     10 Mar 2001: Task = BR implemented.
!     09 May 2001: For BRfac, BSfac, force Lmax1 <= BRtol1, BStol1
!                  (with TCP).
!     20 Dec 2004: Reduce BStol1 from 3.99 to 2.5.
!     20 Dec 2004: Current version of s2BLU.
!     ==================================================================
      integer
     &     i, ir, j, k, nz, iQ, lenc, lenri, locc, locr, lprDbg,
     &     LUprnt, iPloc, iQloc, Lmax1, keepLU, minlen,
     &     oldTP , TPivot
      double precision
     &     BRtol1, BStol1, growth, oldL1
!     ------------------------------------------------------------------
      double precision   one
      parameter         (one = 1.0d+0)
      integer            TPP,     TRP
      parameter         (TPP = 0, TRP = 1)
      integer            B,          BR,         BS
      parameter         (B      = 0, BR     = 1, BS     = 2)
      parameter         (Lmax1  = 151) ! max L-multiplier in factor
      parameter         (LUprnt = 152) !
      parameter         (TPivot = 156) ! 0(1) TPP(TRP)
      parameter         (keepLU = 158) !
      parameter         (minlen = 163) ! minimum recommended lenaLU
!     ------------------------------------------------------------------
      lprDbg    = iw( 85) ! > 0    => private debug print
      iQ        = iw(364) !
      lenc      = iw(365) !
      lenri     = iw(366) !
      locc      = iw(367) !
      locr      = iw(368) !
      iPloc     = iw(369) !
      iQloc     = iw(370) !

      oldL1     = rw(Lmax1)    ! Save Lmax1
      oldTP     = iw(TPivot)   ! Save the TP state

      BRtol1    = 2.5d+0  ! Max Lmax1 for BRfac
!!    BRtol1    = 3.99d+0 !!!!!! TEST FOR SIGEST PAPER (DRCAVTY2)
      BStol1    = 2.5d+0  ! Max Lmax1 for BSfac
      iExit     = 0

!     --------------------------------------------------
!     Set Print level for lu1fac.
!        iw(LUprnt) = 1  for errors,
!                   = 10 for statistics
!                   = 50 for debug info
!     --------------------------------------------------
      iw(LUprnt) = min( lPrint, 10 )
      if (lprDbg .eq. 51) iw(LUprnt) = 50

!     --------------------------------------------------
!     Set other parameters for lu1fac.
!     --------------------------------------------------
      if      (Task .eq. B) then
         iw(keepLU) = 1

      else if (Task .eq. BR) then ! Force TRP
         if (oldTP .eq. TPP) then
            iw(TPivot) = TRP
         !! iw(TPivot) = TPP !!!!! TEST FOR SIGEST PAPER (DRCAVTY2)
         end if
         iw(keepLU) = 0
         rw(Lmax1)  = min( oldL1, BRtol1 )  ! Make sure Lmax1 is reasonably tight.

      else if (Task .eq. BS) then ! Stay with current TPP or TRP
         iw(keepLU) = 0
         rw(Lmax1)  = min( oldL1, BStol1 )  ! Make sure Lmax1 is reasonably tight.
      end if

      if (Task .eq. B  .or.  Task .eq. BR) then
!        -------------------------------------
!        Estimate the number of nonzeros in B.
!        -------------------------------------
         nz   = 0
         do k = 1, m
            j = kBS(k)
            if (j .eq. 0) then
!              --------------------------
!              Relax, just a zero column.
!              --------------------------
            else if (j .le. n) then
!              ----------------------------
!              Basic column from A.
!              ----------------------------
               nz = nz + locA(j+1) - locA(j)
            else
!              ---------------------
!              Basic slack.
!              ---------------------
               nz = nz + 1
            end if
         end do

         iw(minlen) = nz*5/4
         if (iw(minlen) .gt. lenaLU) then
            iExit = 7
            go to 900
         end if

!        ---------------------------------------------------------------
!        Load B into LUSOL.
!        ---------------------------------------------------------------
         nz   = 0
         do k = 1, m
            j = kBS(k)
            if (j .eq. 0) then
!              --------------------------
!              Relax, just a zero column.
!              --------------------------
            else if (j .le. n) then
!              ---------------------
!              Basic column from A.
!              ---------------------
               do i  = locA(j), locA(j+1)-1
                  ir = indA(i)
                  nz = nz + 1
                  aLU (nz) = Acol(i)
                  indc(nz) = ir
                  indr(nz) = k
               end do
            else
!              ---------------------
!              Basic slacks.
!              ---------------------
               nz       =   nz + 1
               aLU (nz) = - one
               indc(nz) =   j - n
               indr(nz) =   k
            end if
         end do

!        iy and iy1 are work vectors
!        y is an output parameter, used by s2sing.

         call lu1fac
     &      ( m, m, nz, lenaLU, iw(151), rw(151),
     &        aLU, indc, indr, iP, iw(iQ),
     &        iw(lenc), iw(lenri), iw(locc), iw(locr),
     &        iw(iPloc), iw(iQloc), iy, iy1, y, iExit )

!        Test for excessive growth in U.

         growth    = rw(166) ! TPP: Umax/Amax    TRP: Akmax/Amax

         if (iExit .eq. 0  .and.  growth .ge. 1.0d+8) then
            iExit = 2
         end if


      else if (Task .eq. BS) then
!        ---------------------------------------------------------------
!        Factorize (B S)' = LU without keeping L and U.
!        ---------------------------------------------------------------
!        Extract (B S)'.
!        iP(1:m) is needed for workspace.
!        iP(i) = 0 except for rows with a basic or superbasic slack.
!        We can ignore all of these rows except for the slack itself.
!        06 Mar 1994: Keep all rows.  (Made a difference in MINOS)
!        22 Apr 1994: Make sure the objective slack remains basic!!
!        11 Nov 1994: Go back to ignoring rows with a slack in B.
!                   This means we don't have to worry about the Obj.
!        ---------------------------------------------------------------
         call iload ( m, 0, iP, 1 )
         do k = 1, nBS
            j = kBS(k)
            if (j .gt. n) iP(j-n) = 1
         end do

!        Count the number of nonzeros in ( B S ).

         nz   = 0
         do k = 1, nBS
            j = kBS(k)
            if (j .le. n) then
               do i  = locA(j), locA(j+1)-1
                  ir = indA(i)
                  if (iP(ir) .eq. 0) nz = nz + 1
               end do
            else
               nz  =  nz + 1
            end if
         end do

         iw(minlen) = nz*5/4
         if (iw(minlen) .gt. lenaLU) then
            iExit = 7
            go to 900
         end if

         nz   = 0
         do k = 1, nBS
            j = kBS(k)
            if (j .le. n) then
               do i  = locA(j), locA(j+1)-1
                  ir = indA(i)
                  if (iP(ir) .eq. 0) then
                     nz       = nz + 1
                     aLU(nz)  = Acol(i)
                     indc(nz) = k
                     indr(nz) = ir
                  end if
               end do
            else
               nz       =   nz + 1
               aLU(nz)  = - one
               indc(nz) =   k
               indr(nz) =   j - n
            end if
         end do

         call lu1fac
     &      ( nBS, m, nz, lenaLU, iw(151), rw(151),
     &        aLU, indc, indr, iP, iw(iQ),
     &        iw(lenc ), iw(lenri), iw(locc), iw(locr),
     &        iw(iPloc), iw(iQloc), iy, iy1, y, iExit )
      end if

!     --------------------------------------------------
!     Restore Lmax1 etc. for BR and BS.
!     --------------------------------------------------
  900 rw(Lmax1)  = oldL1
      iw(TPivot) = oldTP

      return

      end ! subroutine s2BLU

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2Bmod
     &   ( iExit, jrep, m, z, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iExit, jrep, m, leniw, lenrw, iw(leniw)
      double precision
     &     z(m), rw(lenrw)

!     ==================================================================
!     s2Bmod  updates the LU factors of B when column "jrep" is replaced
!             by a vector  v.  On entry,   z  must satisfy  L z = v.
!             It is overwritten.
!
!     20 Oct 1990  Initial version.
!     14 Jul 1997: Thread-safe version.
!     01 Jun 1999: Current version of s2Bmod.
!     ==================================================================
      integer
     &     lenaLU, iP, iQ, lenc, lenri, locr, locc, LUa, indc, indr
      double precision
     &     diag, zNorm
!     ------------------------------------------------------------------
      lenaLU    = iw(213) ! space allotted for LU factors
      iP        = iw(363) !
      iQ        = iw(364) !
      lenc      = iw(365) !
      lenri     = iw(366) !
      locc      = iw(367) !
      locr      = iw(368) !
      LUa       = iw(371) !
      indc      = iw(373) !
      indr      = iw(374) !

      call lu8rpc
     &   ( 1, 2, m, m, jrep, z, z,
     &     lenaLU, iw(151), rw(151),
     &     rw(LUa ), iw(indc), iw(indr), iw(iP), iw(iQ),
     &     iw(lenc), iw(lenri), iw(locc), iw(locr),
     &     iExit, diag, zNorm )

      end ! subroutine s2Bmod

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2Bmod2
     &   ( iExit, jrep, m, z, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iExit, jrep, m, leniw, lenrw, iw(leniw)
      double precision
     &     z(m), rw(lenrw)

!     ==================================================================
!     s2Bmod2  updates the LU factors of B with tight singularity tol.
!
!     22 Jun 2004  Initial version.
!     22 Jun 2004: Current version of s2Bmod2.
!     ==================================================================
      integer
     &     Utol1, Utol1m, Utol2, Utol2m
      double precision
     &     Utol1s, Utol2s
!     ------------------------------------------------------------------
      parameter         (Utol1m    =  63) ! abs tol for small diag of U in LU mod
      parameter         (Utol2m    =  64) ! rel tol for small diag of U in LU mod
      parameter         (Utol1     = 154) ! abs tol for small diag of U.
      parameter         (Utol2     = 155) ! rel tol for small diag of U.
!     ------------------------------------------------------------------
      Utol1s = rw(Utol1)
      Utol2s = rw(Utol2)

      rw(Utol1) = rw(Utol1m)
      rw(Utol2) = rw(Utol2m)

      call s2Bmod
     &   ( iExit, jrep, m, z, iw, leniw, rw, lenrw )

      rw(Utol1) = Utol1s
      rw(Utol2) = Utol2s

      end ! subroutine s2Bmod2

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2Bsol
     &   ( iExit, Task, m, z, y, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Task, iExit, m, leniw, lenrw, iw(leniw)
      double precision
     &     z(m), y(m), rw(lenrw)

!     ==================================================================
!     s2Bsol  solves various systems with the LU factors of B.
!     Task  selects one of the following:
!      Task          Action
!      ----          ------
!      with L    Solve  L*z = z(input).    y  is not touched.
!      with B    Solve  L*z = z(input)  and solve  B*y = z(input).
!      with Bt   Solve  B(transpose)*y = z.  Note that  z  is destroyed.
!
!     20 Oct 1990: Initial version.
!     16 Nov 2001: dnormi added.
!     06 Aug 2003: adopted snEXIT.
!     06 Aug 2003: Current version of s2Bsol
!     ==================================================================
      integer
     &     iP, iQ, lenc, lenri, locc, locr, indc, indr, LUa, lenaLU
      double precision
     &     small0, dnormi
!     ------------------------------------------------------------------
      integer            small
      parameter         (small     = 153) ! defn of small real

      integer            WithL,     WithB,     WithBt
      parameter         (WithL = 0, WithB = 1, WithBt = 2)
!     ------------------------------------------------------------------
      lenaLU    = iw(213) ! space allotted for LU factors
      iP        = iw(363) !
      iQ        = iw(364) !
      lenc      = iw(365) !
      lenri     = iw(366) !
      locc      = iw(367) !
      locr      = iw(368) !
      LUa       = iw(371) !
      indc      = iw(373) !
      indr      = iw(374) !

      if (Task .eq. WithL  .or.  Task .eq. WithB) then
!        ---------------------------------------------------------------
!        Solve   L*z = z(input).
!        When LU*y = z is being solved in SNOPT, norm(z) will sometimes
!        be small (e.g. after periodic refactorization).  Hence for
!        solves with L we scale parmlu(3) to alter what lu6sol thinks
!        is small.
!        ---------------------------------------------------------------
         small0    = rw(small)
         if (Task .eq. WithB) rw(small) = small0 * dnormi( m, z, 1 )

         call lu6sol
     &      ( 1, m, m, z, y,
     &        lenaLU,  iw(151), rw(151),
     &        rw(LUa ), iw(indc), iw(indr), iw(iP), iw(iQ),
     &        iw(lenc), iw(lenri), iw(locc), iw(locr),
     &        iExit )
         rw(small) = small0

         if (Task .eq. WithB) then
!           ------------------------------------------------------------
!           Task = solve with B.   Solve  U*y = z.
!           ------------------------------------------------------------
            call lu6sol
     &         ( 3, m, m, z, y,
     &           lenaLU, iw(151), rw(151),
     &           rw(LUa ), iw(indc), iw(indr), iw(iP), iw(iQ),
     &           iw(lenc), iw(lenri), iw(locc), iw(locr),
     &           iExit )
         end if

      else if (Task .eq. WithBt) then
!        ---------------------------------------------------------------
!        Task = solve with B transpose.  Solve  B'*y = z.
!        ---------------------------------------------------------------
         call lu6sol
     &      ( 6, m, m, y, z,
     &        lenaLU, iw(151), rw(151),
     &        rw(LUa ), iw(indc), iw(indr), iw(iP), iw(iQ),
     &        iw(lenc), iw(lenri), iw(locc), iw(locr),
     &        iExit )
      end if

      if (iExit .ne. 0) iExit = 142 ! error in basis package

      end ! subroutine s2Bsol

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2sb
     &   ( m, mBS, nb, nS, nBS, hs, kBS, xBS )

      implicit
     &     none
      integer
     &     m, mBS, nb, nS, nBS, hs(nb), kBS(mBS)
      double precision
     &     xBS(mBS)

!     ==================================================================
!     s2SB  finds if any superbasic slacks were made basic by s2sing.
!     If any are found, nS and nBS are updated and the corresponding
!     entries are removed from  kBS and xBS.
!
!     02 Jul 2005: First version of s2sb.
!     02 Jul 2005: Current version of s2sb.
!     ==================================================================
      integer
     &     j, jq, k, nS0
!     ------------------------------------------------------------------
      if (nS .eq. 0) return

      nS0   = nS
      do jq = nS0, 1, -1
         j  = kBS(m+jq)
         if (hs(j) .eq. 3) then
            nS   = nS - 1
            nBS  = m  + nS
            do k = m+jq, nBS
               kBS(k) = kBS(k+1)
               xBS(k) = xBS(k+1)
            end do
         end if
      end do

      end ! subroutine s2sb

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2sing
     &   ( lPrint, mBS, m, n, nb,
     &     z, iP, iQ, bl, bu, hs, kBS, x, iw, leniw )

      implicit
     &     none
      integer
     &     lPrint, mBS, m, n, nb, leniw,
     &     iP(m), iQ(m), hs(nb), kBS(mBS), iw(leniw)
      double precision
     &     bl(nb), bu(nb), z(m), x(nb)

!     =================================================================
!     s2sing  is called if the LU factorization of the basis appears
!     to be singular.   If  z(j)  is not positive, the  jth  basic
!     variable  kBS(j)  is replaced by the appropriate slack.
!     If any kBS(j) = 0, only a partial basis was supplied.
!
!     30 Sep 1991: First version based on minos routine m2sing.
!     29 May 1995: Optional swapping of slack and basic.
!     12 Jul 1997: Thread-safe version.
!     18 Nov 2001: lPrint added as parameter.
!     31 Jul 2003: snPRNT adopted.
!     31 Jul 2003: Current version of s2sing.
!     =================================================================
      character
     &     str*72
      integer
     &     i, j, k, nSing
!     -----------------------------------------------------------------
      integer            nPrint
      parameter         (nPrint = 5)
      double precision   zero
      parameter         (zero = 0.0d+0)
!     -----------------------------------------------------------------
      nSing  = 0
      do k = 1, m
         j = iQ(k)
         if (z(j) .le. zero) then
            j = kBS(j)

            if (j .gt. 0) then

!              Make variable  j  nonbasic (and feasible).
!              hs(j) = -1 means x(j) is strictly between its bounds.

               if      (x(j) .le. bl(j)) then
                  x(j)  =  bl(j)
                  hs(j) =  0
               else if (x(j) .ge. bu(j)) then
                  x(j)  =  bu(j)
                  hs(j) =  1
               else
                  hs(j) = -1
               end if
               if (bl(j) .eq. bu(j)) hs(j) = 4
            end if

!           Make the appropriate slack basic.

            i       = iP(k)
            hs(n+i) = 3
            nSing   = nSing + 1
            if (lPrint .ge. 10  .and.  nSing .le. nPrint) then
               write(str, 1000) j, i
               call snPRNT( 3, str, iw, leniw )
            end if
         end if
      end do

      if (lPrint .ge. 10  .and.  nSing .gt. nPrint) then
         write(str, 1100) nSing
         call snPRNT( 3, str, iw, leniw )
      end if

      return

 1000 format(' Column', i7, '  replaced by slack', i7)
 1100 format(' and so on.  Total slacks inserted =', i6)

      end ! subroutine s2sing

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2tols
     &   ( mode, NewTol, itn, iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     NewTol
      integer
     &     mode, itn, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)

!     =================================================================
!     s2tols  sets the LU Factor and Update tolerances.
!
!     mode   (input) says what should be done:
!     mode = DefTol  Set default tols for current TP (TPP, TRP or TCP).
!     mode = RstTol  Set default tols for TP specified by user.
!     mode = RedTol  Reduce tols or switch from TPP to TRP.
!     mode = MinTol  Set minimum tols for current TP.
!
!     NewTol (output) says if the LU tols were changed
!                     or we switched from TPP to TRP.
!
!     30 Aug 2000: First version of s2tols.
!     10 Dec 2002: Rook and diagonal pivoting options added.
!     31 Jul 2003: snPRNT adopted.
!     21 Dec 2004: Current version of s2tols.
!     =================================================================
      character
     &     str*72
      integer
     &     lvlPiv, TPivot
      double precision
     &     Lmax1, Lmax2, oldL1, oldL2,
     &     tolDcp, tolDpp, tolDrp, tolDUp,
     &     tolFac, tolUpd
!     ------------------------------------------------------------------
      integer            TPP,     TRP,     TCP
      parameter         (TPP = 0, TRP = 1, TCP = 2)
      integer            DefTol,     RedTol,     MinTol,     RstTol
      parameter         (DefTol = 0, RedTol = 1, MinTol = 2, RstTol = 3)

      double precision   tolTPP, tolTRP, tolTCP, tolUpt
      parameter         (tolTPP = 1.011d+0) ! Minimum LU tols allowed
      parameter         (tolTRP = 1.011d+0)
      parameter         (tolTCP = 1.011d+0)
      parameter         (tolUpt = 1.011d+0)
!     ------------------------------------------------------------------
      character          TP(0:2)*8
      data               TP/'partial ', 'rook    ', 'complete'/
!     ------------------------------------------------------------------

      tolFac    = rw( 66) ! LU factor tolerance = user-defined Lmax1
      tolUpd    = rw( 67) ! LU update tolerance = user-defined Lmax2
      lvlPiv    = iw( 80) ! 0(1 2 3) LU Partial (Rook Complete Diag) piv
      TPivot    = iw(156) ! Current LU pivot option

      Lmax1     = rw(151) ! max allowable L-multiplier in factor
      Lmax2     = rw(152) ! max allowable L-multiplier in update
      tolDpp    = rw(181) ! default Lmax1 for TPP
      tolDrp    = rw(187) ! default Lmax1 for TRP
      tolDcp    = rw(182) ! default Lmax1 for TCP
 !    tolDdp    = rw(188) ! default Lmax1 for TDP
      tolDup    = rw(183) ! default Lmax2

      oldL1     = Lmax1
      oldL2     = Lmax2
      newTol    = .false.

      if (mode .eq. DefTol) then
!        ---------------------------------------------------------------
!        Set the default LU Factor tol and LU Update tol.
!        ---------------------------------------------------------------
         if (     TPivot .eq. TPP) then
            Lmax1 = tolDpp
         else if (TPivot .eq. TRP) then
            Lmax1 = tolDrp
         else if (TPivot .eq. TCP) then
            Lmax1 = tolDcp
         end if
         Lmax2  = tolDup

      else if (mode .eq. RstTol) then
!        ---------------------------------------------------------------
!        Reset the user LU Factor tol and LU Update tol.
!        ---------------------------------------------------------------
         TPivot = lvlPiv
         Lmax1  = tolFac
         Lmax2  = tolUpd

      else if (mode .eq. RedTol) then
!        ---------------------------------------------------------------
!        Reduce LU Factor tol and LU Update tol,
!        perhaps switching to TRP first.
!        Note: TCP is not activated here -- only TRP.
!        TCP must be requested via the LU Complete Pivoting option.
!        It is then used for all factorizations.
!        ---------------------------------------------------------------
         if (TPivot .eq. TPP  .and.  Lmax1 .le. tolTPP) then  ! Switch
            TPivot = TRP
            Lmax1  = tolDrp   !!!!! Comment out for SIGEST test (DRCAVTY2)
            Lmax2  = tolDup   !!!!! Comment out for SIGEST test (DRCAVTY2)
         else if (TPivot .eq. TRP  .and.  Lmax1 .le. tolTRP) then ! Switch
            TPivot = TCP
            Lmax1  = tolDcp   !!!!! Comment out for SIGEST test (DRCAVTY2)
            Lmax2  = tolDup   !!!!! Comment out for SIGEST test (DRCAVTY2)
         else
            if (Lmax1 .gt. 4.0d+0) then
                Lmax1  = Lmax1 * 0.5d+0
            else
                Lmax1  = sqrt( Lmax1 )
            end if
            if (Lmax2 .gt. 4.0d+0) then
                Lmax2  = Lmax2 * 0.5d+0
            else
                Lmax2  = sqrt( Lmax2 )
            end if
         end if

      else if (mode .eq. MinTol) then
!        ---------------------------------------------------------------
!        Set LU Factor tol and LU Update tol to their minimum values.
!        ---------------------------------------------------------------
         if (     TPivot .eq. TPP) then
            Lmax1  = tolTPP
         else if (TPivot .eq. TRP) then
            Lmax1  = tolTRP
         else if (TPivot .eq. TCP) then
            Lmax1  = tolTCP
         end if
         Lmax2  = tolUpt
      end if

!     ------------------------------------------------------------------
!     Make sure the tols aren't too small or big.
!     ------------------------------------------------------------------
      if (     TPivot .eq. TPP) then
         Lmax1  = max( Lmax1, tolTPP )
      else if (TPivot .eq. TRP) then
         Lmax1  = max( Lmax1, tolTRP )
      else if (TPivot .eq. TCP) then
         Lmax1  = max( Lmax1, tolTCP )
      end if
      Lmax2  = max( Lmax2, tolUpt )
      Lmax2  = min( Lmax1, Lmax2  )  ! Update tol never more than Factol

!     ------------------------------------------------------------------
!     Print significant changes.
!     ------------------------------------------------------------------
      NewTol = Lmax1 .ne. oldL1  .or.  Lmax2 .ne. oldL2

!!    TPivot = 0 !!!!!! TEST FOR SIGEST PAPER (DRCAVTY2)

      if (NewTol) then
         write(str, 1000) itn, TP(TPivot), Lmax1, Lmax2
         call snPRNT( 23, str, iw, leniw )
      end if

      iw(156) = TPivot
      rw(151) = Lmax1
      rw(152) = Lmax2

      return

 1000 format(' Itn', i7, ': LU ', a, ' pivoting tols ', 2f10.2)

      end ! subroutine s2tols

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s2tryLU
     &   ( itn, LUreq0, nS, LUreq, LUok, typeLU,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     LUok
      integer
     &     itn, LUreq0, LUreq, nS, typeLU, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)

!     =================================================================
!     s2tryLU  is called when numerical difficulties imply that a new
!     factorize is needed.
!
!     On exit:
!        LUok    says if a new LU is possible.
!        typeLU  defines the type of factorization to be done.
!        LUreq   is used for printing and indicates the reason for
!                the LU request.
!
!     26 Dec 2003: First version of s2tryLU.
!     09 Apr 2004: Current version of s2tryLU.
!     =================================================================
      logical
     &     NewTol
      integer
     &     LUmod, nBfac
!     ------------------------------------------------------------------
      integer            B,          BS        , BT
      parameter         (B      = 0, BS     = 2, BT     = 3)
      integer            RedTol
      parameter         (RedTol = 1)
!     ------------------------------------------------------------------
      nBfac  = iw(211) ! number of consecutive `B' factorizes
      LUmod  = iw(216) ! number of LU mods since the last factorize

      LUreq  = LUreq0
      LUok   = .true.

      if (nBFac .gt. 1  .and.  nS .gt. 1) then

!        The LU has been computed since the last BS factorize.
!        Try to find a better basis.

         typeLU = BS

      else if (LUmod .gt. 0 .and. (LUmod .ne. 1 .or. LUreq .ne. 5)) then

!        The LU has been modified since the last factorize.
!        Refactorize with the current LU tolerances.
!        A BS factorize will be used if needed.

         typeLU = BT

      else

!        We have the LU at the current point.
!        Refactorize with smaller LU tolerances.

         call s2tols
     &     ( RedTol, NewTol, itn, iw, leniw, rw, lenrw )

         if (NewTol) then
            LUok = .true.
            if (nS .eq. 0) then
               typeLU = B
            else
               typeLU = BS
            end if
         else
            LUok  = .false.
            LUreq = 0
         end if
      end if

      end ! subroutine s2tryLU

