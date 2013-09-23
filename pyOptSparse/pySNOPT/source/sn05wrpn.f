*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  sn05wrpn.f  --- user-function wrapper for  npOpt.
*
*     s0fgN
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s0fgN
     &   ( iExit, modefg, Status, getCon, getObj,
     &     n, negCon, nnCon0, nnCon, nnJac, nnL, ngObj0, ngObj,
     &     fgcon, fgobj, x,
     &     ne , nlocJ, locJ, indJ,
     &     fCon, fObj, gCon, gObj,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     fgcon, fgobj
      logical
     &     getCon, getObj
      integer
     &     iExit, lencu, leniu, lenru, lencw, leniw, lenrw, modefg,
     &     n, ne, negCon, nlocJ, nnCon0, nnCon, nnJac, nnL, ngObj0,
     &     ngObj, Status, indJ(ne), locJ(nlocJ), iu(leniu), iw(leniw)
      double precision
     &     fObj, fCon(nnCon0),  gObj(ngObj0), gCon(negCon), x(n),
     &     ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     s0fgN   is an instance of fgwrap that calls the user-written
*     routine  userfun  to evaluate the problem functions and possibly
*     their gradients.
*
*     Arguments  fgcon  and  fgobj  are called using modefg to control
*     the gradients as follows:
*
*     Version for the NPSOL wrapper snOptN.
*
*     modefg        Task
*     ------        ----
*       2     Assign fCon, fObj and all known elements of gCon and gObj.
*       1     Assign all known elements of gCon and gObj.
*             (fObj and fCon are ignored).
*       0     Assign fObj, fCon.  (gCon and gObj are ignored).
*
*     If s0fgN is called with minmax = 0 (feasible point only) then
*     ngObj = max(nnJac,nnObj)  and the user objective is not used.
*
*     09-Jan 1992: First version of s0fgN  based on snwrap.
*     03 Aug 2003: snPRNT adopted.
*     04 Jan 2005: v, Hv added.
*     01 Apr 2005: Current version of s0fgN.
*     ==================================================================
      character
     &     str*80
      external
     &     ddot
      logical
     &     FPonly, scaled
      integer
     &     gotFD, gotG, gotGl, l, lG, lgConu, lgObju,
     &     liy1, lvlScl, lvlTim, lvlDer, lAscal, lx0, lxscal,
     &     minmax, modeC, modeF, nfCon1, nfCon2, nfObj1, nfObj2, ngrad,
     &     nnGlin, nnObj
      double precision
     &     ddot, Gdummy
*     ------------------------------------------------------------------
      parameter         (lvlDer =  70) ! = 0,1,2 or 3, deriv. level
      parameter         (gotFD  = 183) ! > 0 => some differences needed
      parameter         (gotG   = 184) ! > 0 => some exact derivs
      parameter         (gotGl  = 185) ! > 0 => constant Jacob elements
      parameter         (nfCon1 = 189) ! calls to fCon: mode = 0
      parameter         (nfCon2 = 190) ! calls to fCon  mode > 0
      parameter         (nfObj1 = 194) ! calls to fObj: mode = 0
      parameter         (nfObj2 = 195) ! calls to fObj: mode > 0
      double precision   half,            one
      parameter         (half   = 0.5d+0, one   = 1.0d+0)
*     ------------------------------------------------------------------
      nnObj     = iw( 22) ! # of objective variables
      lvlScl    = iw( 75) ! scale option
      minmax    = iw( 87) ! 1, 0, -1  => MIN, FP, MAX
      lvlTim    = iw(182) ! Timing level

      lAscal    = iw(295) ! Ascale(nb) = row and column scales
      lx0       = iw(298) ! x0(nnL)    = Feasible starting point
      lxscal    = iw(302) ! xscl(n)    = copy of scaled x
      liy1      = iw(309) ! iy1(nb)    =  integer work vector
      lgConu    = iw(319) ! record of unknown derivatives and constants
      lgObju    = iw(323) ! record of unknown derivatives

      Gdummy    = rw( 69) ! definition of an 'unset' value

      iExit     = 0

      FPonly    = minmax .eq. 0
      scaled    = lvlScl .eq. 2

      modeC     = modefg
      modeF     = modefg

      if (Status .eq. 1) then
*        ---------------------------------------------------------------
*        First evaluation of the problem functions
*        On entry, lvlScl = 0.
*        ---------------------------------------------------------------
         iw(gotFD) =  0
         iw(gotG)  =  0
         iw(gotGl) =  0
         call snPRNT( 13, ' ', iw, leniw )
         call dload ( negCon, Gdummy, gCon, 1 )
         call dload ( ngObj , Gdummy, gObj, 1 )
      end if

*     ------------------------------------------------------------------
*     Unscale x.
*     ------------------------------------------------------------------
      if ( scaled ) then
         call dcopy ( nnL, x         , 1, rw(lxscal), 1 )
         call ddscl ( nnL, rw(lAscal), 1, x         , 1 )

*        If the Jacobian has some constant elements, they are wrecked
*        by the scaling.  Restore them from gConu.

         if ( getCon ) then
            if (modefg .gt. 0  .and.  iw(gotGl) .gt. 0) then
               call dcopy ( negCon, rw(lgConu), 1, gCon, 1 )
            end if
         end if
      end if

*     ------------------------------------------------------------------
*     Compute the constraints.
*     ------------------------------------------------------------------
*     To incorporate user workspace in fgcon, replace the next
*     call to fgcon with:
*     call fgcon ( modeC, nnCon, nnJac, nnCon,
*    &             iw(liy1), x, fCon, gCon, Status,
*    &             cu, lencu, iu, leniu, ru, lenru )

      if ( getCon ) then
         if (lvlTim .ge. 2) call s1time( 4, 0, iw, leniw, rw, lenrw )
         call iload ( nnCon, (1), iw(liy1), 1 )
         call fgcon
     &      ( modeC, nnCon, nnJac, nnCon,
     &        iw(liy1), x, fCon, gCon, Status )
         if (lvltim .ge. 2) call s1time(-4, 0, iw, leniw, rw, lenrw )

         iw(nfCon1) = iw(nfCon1) + 1
         if (modefg .gt. 0)
     &        iw(nfCon2) = iw(nfCon2) + 1
      end if

*     ------------------------------------------------------------------
*     Compute the objective.
*     ------------------------------------------------------------------
*     To incorporate user workspace in fgobj, replace the next
*     call to fgobj with:
*     call fgobj ( modeF, ngObj, x, fObj, gObj, Status,
*    &             cu, lencu, iu, leniu, ru, lenru )

      if (getObj  .and.  modeC .ge. 0) then
         if (lvlTim .ge. 2) call s1time( 5, 0, iw, leniw, rw, lenrw )
         if ( FPonly ) then
            call dcopy ( ngObj, x, 1, gObj, 1 )
            call daxpy ( ngObj, (-one), rw(lx0), 1, gObj, 1 )
            fObj = half*ddot ( ngObj, gObj, 1, gObj, 1 )
         else ! ngObj = nnObj
            call fgobj ( modeF, ngObj, x, fObj, gObj, Status )
         end if

         if (lvlTim .ge. 2) call s1time(-5, 0, iw, leniw, rw, lenrw )
         iw(nfObj1) = iw(nfObj1) + 1
         if (modefg .gt. 0) iw(nfObj2) = iw(nfObj2) + 1
      end if

*     ------------------------------------------------------------------
*     Scale  x and the derivatives.
*     ------------------------------------------------------------------
      if ( scaled ) then
         call dcopy ( nnL, rw(lxscal), 1, x, 1 )

         if ( getCon ) then
            call dddiv ( nnCon, rw(lAscal+n), 1, fCon, 1 )
            if (modefg .gt. 0  .and.  iw(gotG) .gt. 0) then
               call s8sclJ
     &            ( nnCon, nnJac, negCon, n, rw(lAscal),
     &              ne, nlocJ, locJ, indJ, gCon, rw, lenrw )
            end if
         end if

         if (getObj  .and.  modeC .ge. 0) then
            if (modefg .gt. 0  .and.  iw(gotG) .gt. 0) then
               call s8sclg
     &            ( ngObj, rw(lAscal), gObj, rw, lenrw )
            end if
         end if
      end if

      if (modeC .lt. 0  .or.  modeF .lt. 0) then
*        ---------------------------------------------------------------
*        The user may be saying the function is undefined (mode = -1)
*        or may just want to stop                         (mode < -1).
*        ---------------------------------------------------------------
         if (modeC .eq. -1  .or.  modeF .eq. -1) then
            iExit = -1
         else
            if (modeC .lt. 0) then
               iExit = 72
            else
               iExit = 73
            end if
         end if
      end if

*     ==================================================================
*     Do some housekeeping on the first entry.
*     ==================================================================
      if (Status .eq. 1  .and.  iExit .eq. 0) then
         if ( getCon ) then
*           ------------------------------------------------------------
*           Count how many Jacobian elements are provided.
*           ------------------------------------------------------------
            nnGlin = 0
            ngrad  = 0
            do l = 1, negCon
               if (gCon(l) .ne. Gdummy) ngrad  = ngrad + 1
            end do

            write(str, 1100) ngrad, negCon
            call snPRNT( 3, str, iw, leniw )

            if (ngrad .lt. negCon) then

*              Some Jacobian elements are missing.

               if (iw(lvlDer) .ge. 2) then
*                 ------------------------------------------------------
*                 All the Jacobian is known.  Any undefined elements
*                 are assumed constant, and are restored from gConu.
*                 ------------------------------------------------------
                  write(str, 3100)
                  call snPRNT( 3, str, iw, leniw )

                  lG  = lgConu
                  do l  = 1, negCon
                     if (gCon(l) .eq. Gdummy) then
                        gCon(l) = rw(lG)
                        nnGlin  = nnGlin + 1
                     end if
                     lG = lG + 1
                  end do
               else
*                 ------------------------------------------------------
*                 Save a permanent copy of gCon in gConu so that we know
*                 which derivatives must be estimated.
*                 ------------------------------------------------------
                  call dcopy ( negCon, gCon, 1, rw(lgConu), 1 )
               end if
            end if ! ngrad < negCon
            if (ngrad + nnGlin .lt. negCon) iw(gotFD) = 1
            if (ngrad          .gt.      0) iw(gotG ) = 1
            if (nnGlin         .gt.      0) iw(gotGl) = 1
         end if

         if ( getObj ) then
*           ------------------------------------------------------------
*           Count how many gradient elements are known.
*           (These may be the gradients of the FP objective.)
*           ------------------------------------------------------------
            ngrad  = 0
            do l = 1, ngObj
               if (gObj(l) .ne. Gdummy) ngrad = ngrad + 1
            end do

            if ( FPonly ) then
               write(str, 2010) ngObj
               call snPRNT( 3, str, iw, leniw )
            else
               write(str, 2000) ngrad, nnObj
               call snPRNT( 3, str, iw, leniw )
            end if

            if (ngrad .lt. ngObj) then

*              Some objective gradients are missing.

               if (iw(lvlDer) .eq. 1  .or.  iw(lvlDer) .eq. 3) then
*                 ------------------------------------------------------
*                 The objective gradient was meant to be known.
*                 ------------------------------------------------------
                  iw(lvlDer) = iw(lvlDer) - 1
                  write(str, 2100) iw(lvlDer)
                  call snPRNT( 3, str, iw, leniw )
               end if
*              ---------------------------------------------------------
*              Copy gObj into gObju.
*              ---------------------------------------------------------
               call dcopy ( ngObj, gObj, 1, rw(lgObju), 1 )
            end if
            if (ngrad .lt. ngObj) iw(gotFD) = 1
            if (ngrad .gt.     0) iw(gotG ) = 1
         end if
      end if

      return

 1100 format(  ' The user has defined', i8, '   out of', i8,
     &         '   constraint gradients.')
 2000 format(  ' The user has defined', i8, '   out of', i8,
     &         '   objective  gradients.')
 2010 format(  ' NpOpt   will define ', i8, '   gradients for the ',
     &         ' FP objective.')
 2100 format(  ' XXX  Some objective  derivatives are missing ---',
     &         ' derivative level reduced to', i3)
 3100 format(  ' ==>  Some constraint derivatives are missing, ',
     &         ' assumed constant.'/ )

      end ! subroutine s0fgN
