*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  np02lib.f
*
*     npTitl   npInit   npSpec   npMem
*     npSet    npSeti   npSetr
*     npGet    npGeti   npGetr
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npTitl( title )

      character
     &     title*30

*     ==================================================================
*     npTitl sets the title.
*     ==================================================================

      title  = 'N P O P T  7.2-5    (May 2007)'
*---------------123456789|123456789|123456789|--------------------------

      end ! subroutine npTitl

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npInit
     &   ( iPrint, iSumm, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iPrint, iSumm, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)

*     ==================================================================
*     npInit  is called by the user to do the following:
*     1. Open default files (Print, Summary).
*     2. Initialize title.
*     3. Set options to default values.
*
*     15 Nov 1991: First version.
*     14 Jul 1997: Thread-safe version.
*     21 Mar 1997: First version based on snopt routine snInit
*     14 Jul 1997: Thread-safe version.
*     02 Oct 1997: Character workspace added.
*     28 Sep 2003: Current version of npInit.
*     ==================================================================
      external
     &     s1outpt
      character
     &     Solver*6, str*80, str2*80
      integer
     &     inform, iSpecs, iStdo, lencw, lvlTim,
     &     maxru, maxrw, maxiu, maxiw, maxcu, maxcw, s1outpt
      character
     &     title*30
*     ------------------------------------------------------------------
      parameter         (maxru     =   2) ! start of SNOPT part of rw
      parameter         (maxrw     =   3) ! end   of SNOPT part of rw
      parameter         (maxiu     =   4) ! start of SNOPT part of iw
      parameter         (maxiw     =   5) ! end   of SNOPT part of iw
      parameter         (maxcu     =   6) ! start of SNOPT part of cw
      parameter         (maxcw     =   7) ! end   of SNOPT part of cw
      parameter         (lvlTim    = 182) ! Timing level
      parameter         (lencw     = 500)
      character       cw(lencw)*8
*     ------------------------------------------------------------------
      character          dashes*30
      data               dashes /'=============================='/
*     ------------------------------------------------------------------
      Solver = 'NPINIT'

      if (leniw .lt. 500 .or. lenrw .lt. 500) then
*        ---------------------------------------------------------------
*        Not enough workspace to do ANYTHING!
*        Print and exit without accessing the work arrays.
*        ---------------------------------------------------------------
         inform = 81       ! Work arrays must have at least 500 elements
         call snWRAP( inform, Solver, str, str2, iw, leniw )
         go to 999
      end if

      iSpecs    = 0
      iStdo     = s1outpt( )
      iw( 10)   = iStdo   ! Standard Output
      iw( 11)   = iSpecs
      iw( 12)   = iPrint  ! Print file
      iw( 13)   = iSumm   ! Summary file

      iw(maxcu) = 500
      iw(maxiu) = 500
      iw(maxru) = 500
      iw(maxcw) = lencw
      iw(maxiw) = leniw
      iw(maxrw) = lenrw

      call npTitl( title )
      call s1init( title, iw, leniw, rw, lenrw )

      call snPRNT(11, '         '//dashes, iw, leniw )
      call snPRNT( 1, '         '//title , iw, leniw )
      call snPRNT( 1, '         '//dashes, iw, leniw )

      call snPRNT(12, ' '//dashes, iw, leniw )
      call snPRNT( 2, ' '//title , iw, leniw )
      call snPRNT( 2, ' '//dashes, iw, leniw )

*     ------------------------------------------------------------------
*     Set the options to default values.
*     npopt  will check the options later and maybe print them.
*     ------------------------------------------------------------------
      call s3undf( cw, lencw, iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     Initialize some global values.
*     ------------------------------------------------------------------
      iw(lvlTim) = 3

  999 return

      end ! subroutine npInit

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npSpec
     &   ( iSpecs, iExit, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iSpecs, iExit, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)

*     ==================================================================
*     npSpec  is called by the user to read the Specs file.
*
*     07 Feb 1998: First version.
*     01 Aug 2003: s3file now has a "title" parameter.  Use ' '.
*     27 Oct 2003: Current version of npSpec.
*     ==================================================================
      character
     &     Solver*6, str*80, str2*80
      integer
     &     Errors, Calls, iPrint, iSumm, lencw
      external
     &     s3opt
*     ------------------------------------------------------------------
      parameter            (lencw     = 500)
      character          cw(lencw)*8
*     ------------------------------------------------------------------
      Solver = 'NPSPEC'

      if (lencw .lt. 500 .or. leniw .lt. 500 .or. lenrw .lt. 500) then
*        ---------------------------------------------------------------
*        Not enough workspace to do ANYTHING!
*        Print and exit without accessing the work arrays.
*        ---------------------------------------------------------------
         iExit = 81        ! Work arrays must have at least 500 elements
         call snWRAP( iExit, Solver, str, str2, iw, leniw )
         go to 999
      end if

      if (iSpecs .le. 0) then
         iExit = 131
         go to 800
      end if

      iw( 11)   = iSpecs  ! Specs (options) file

      iPrint    = iw( 12) ! Print file
      iSumm     = iw( 13) ! Summary file

      iExit     = 0
      Calls     = 1

*     ------------------------------------------------------------------
*     Read the Specs file.
*     npopt  will check the options later and maybe print them.
*     ------------------------------------------------------------------
      call s3file
     &   ( iExit, Calls, iSpecs, s3opt, ' ', iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

  800 if (iExit .eq. 0) then
         iExit = 101            ! SPECS file read successfully
      end if

      call snWRAP( iExit, Solver, str, str2, iw, leniw )

  999 return

      end ! subroutine npSpec

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npMem
     &   ( iExit, n, nclin, ncnln,
     &     mincw, miniw, minrw,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iExit, n, nclin, ncnln, mincw, miniw, minrw, leniw, lenrw,
     &     iw(leniw)
      double precision
     &     rw(lenrw)

*     ==================================================================
*     npMem   estimates the memory requirements for npopt,
*     using the values:
*     n      the number of variables (dimension of  x),
*
*     nclin  the number of linear constraints (rows of the matrix  A),
*
*     ncnln  the number of nonlinear constraints (dimension of  c(x)),
*
*     These values are used to compute the minimum required storage:
*     miniw, minrw.
*
*     Note:
*     1. All default parameters must be set before calling npMem,
*        since some values affect the amount of memory required.
*
*     2. The arrays rw and iw hold  constants and work-space addresses.
*        They must have dimension at least 500.
*
*     3. This version of npMem does not allow user accessible
*        partitions of iw and rw.
*
*     01 May 1998: First version.
*     19 Feb 2004: Current version of npMem.
*     ==================================================================
      character
     &     Solver*6, str*80, str2*80
      logical
     &     PrtMem
      integer
     &     Useriw(130)
      double precision
     &     Userrw(130)
      integer
     &     i, inform, lencw,
     &     lenR, liwEst, lrwEst, llenrw, lleniw, llencw, lvlHes, m,
     &     maxcw, maxiw, maxrw, maxR, maxS, mQNmod, nCon, ne, negCon,
     &     nextcw, nextiw, nextrw, nkx, nnCon, nnJac, nnObj, nnCol
*     ------------------------------------------------------------------
      parameter        (lencw     = 500)
      character      cw(lencw)*8
      character      cdummy*8
      parameter     (cdummy = '-1111111')
*     ------------------------------------------------------------------
      Solver = 'NPMEM '
      iExit  = 0

      if (lencw .lt. 500 .or. leniw .lt. 500 .or. lenrw .lt. 500) then
*        ---------------------------------------------------------------
*        Not enough workspace to do ANYTHING!
*        Print and exit without accessing the work arrays.
*        ---------------------------------------------------------------
         iExit = 81        ! Work arrays must have at least 500 elements
         call snWRAP( iExit, Solver, str, str2, iw, leniw )
         go to 999
      end if

*     Save the user's option choices  (weird choices get overwritten).

      call icopy ( 130, iw(51), 1, Useriw, 1 )
      call dcopy ( 130, rw(51), 1, Userrw, 1 )

*     Assign fake values for lencw, leniw, lenrw.
*     This will force s2Mem to estimate the memory requirements.

      llenrw  = 500
      lleniw  = 500
      llencw  = 500

*     Compute the problem dimensions.

      nCon       = nclin + ncnln

      if (nCon .eq. 0) then

*        The problem is unconstrained.
*        A dummy row of zeros will be included.

         nnCol = 0
         m     = 1
         ne    = 1
      else
         nnCol = n
         m     = nCon
         ne    = m*n
      end if

      negCon     = ncnln*n
      nnCon      = ncnln
      nnJac      = nnCol
      nnObj      = n

*     An obligatory call to snInit has `undefined' all options.
*     However, it could not undefine the char*8 options.  Do it now.
*     Check the user-defined values and assign undefined values.
*     s8dflt needs various problem dimensions in iw.

      do i = 51, 180
         cw(i)  = cdummy
      end do

      iw( 15) = n     ! copy of the number of columns
      iw( 16) = m     ! copy of the number of rows
      iw( 17) = ne    ! copy of the number of nonzeros in Jcol
      iw( 21) = nnJac ! # nonlinear Jacobian variables
      iw( 22) = nnObj ! # variables in gObj
      iw( 23) = nnCon ! # of nonlinear constraints

      call s8dflt
     &   ( m, n, nnCon, nnJac, nnObj,
     &     cw, llencw, iw, lleniw, rw, llenrw )

      nextcw  = 501
      nextiw  = 501
      nextrw  = 501

      maxcw   = lencw
      maxiw   = leniw
      maxrw   = lenrw

      nkx     = n + m

      maxR    = iw( 52) ! max columns of R.
      maxS    = iw( 53) ! max # of superbasics
      mQNmod  = iw( 54) ! (ge 0) max # of BFGS updates
      lvlHes  = iw( 72) ! 0,1,2  => LM, FM, Exact Hessian

      lenR    = maxR*(maxR + 1)/2  +  (maxS - maxR)

      call s8Map
     &   ( m, n, negCon, nkx, nnCon, nnJac, nnObj,
     &     lenR, maxR, maxS,  mQNmod, lvlHes,
     &     nextcw, nextiw, nextrw, iw, leniw )
      call s2Bmap
     &   ( m, n, ne, maxS,
     &     nextiw, nextrw, maxiw, maxrw, liwEst, lrwEst, iw, leniw )
      PrtMem = .false.          ! Suppress messages from s2Mem
      call s2Mem
     &   ( inform, PrtMem, liwEst, lrwEst,
     &     nextcw, nextiw, nextrw,
     &     maxcw, maxiw, maxrw, llencw, lleniw, llenrw,
     &     mincw, miniw, minrw, iw )

*     mincw = mincw
      miniw = liwEst
      minrw = lrwEst

*     Restore the user's choices of options.

      call icopy ( 130, Useriw, 1, iw(51), 1 )
      call dcopy ( 130, Userrw, 1, rw(51), 1 )

*     Print the exit conditions.

      if (iExit .eq. 0) then
         iExit = 104            ! memory requirements estimated
      end if

      call snWRAP( iExit, Solver, str, str2, iw, leniw )

  999 return

      end ! subroutine npMem

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npSet
     &   ( buffer, iPrint, iSumm, iExit, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     iPrint, iSumm, iExit, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)

*     ==================================================================
*     npSet  decodes the option contained in  buffer.
*
*     The buffer is output to file iPrint, minus trailing blanks.
*     Error messages are output to files iPrint and iSumm.
*     Buffer is echoed to iPrint but normally not to iSumm.
*     It is echoed to iSumm before any error msg.
*
*     On entry,
*     iPrint is the print   file.  no output occurs if iPrint .le 0.
*     iSumm  is the Summary file.  no output occurs if iSumm  .le 0.
*     iExit  is the number of errors so far.
*
*     On exit,
*     iExit  is the number of errors so far.
*
*     27 Nov 1991: first version of npSet.
*     03 Nov 2000: current version.
*     ==================================================================
      integer
     &     ivalue, lencw
      double precision
     &     rvalue
      character
     &     cvalue*8
      character
     &     key*16
*     ------------------------------------------------------------------
      parameter            (lencw     = 500)
      character          cw(lencw)*8
*     ------------------------------------------------------------------
      call s3opt
     &   ( .true., buffer, key, cvalue, ivalue, rvalue,
     &     iPrint, iSumm, iExit, cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine npSet

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npSeti
     &   ( buffer, ivalue, iPrint, iSumm, iExit, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     ivalue, iPrint, iSumm, iExit, leniw, lenrw,
     &     iw(leniw)
      double precision
     &     rw(lenrw)

*     ==================================================================
*     npSeti decodes the option contained in  buffer // ivalue.
*     The parameters other than ivalue are as in npSet.
*
*     27 Nov 1991: first version of npSeti.
*     03 Nov 2000: current version.
*     ==================================================================
      integer
     &     ivalxx, lenbuf, lencw
      double precision
     &     rvalue
      character
     &     cvalue*8
      character
     &     key*16
      character
     &     buff72*72
*     ------------------------------------------------------------------
      parameter         (lencw     = 500)
      character       cw(lencw)*8
*     ------------------------------------------------------------------
      write(key, '(i16)') ivalue
      lenbuf = len(buffer)
      buff72 = buffer
      buff72(lenbuf+1:lenbuf+16) = key
      ivalxx = ivalue
      call s3opt
     &   ( .true., buff72, key, cvalue, ivalxx, rvalue,
     &     iPrint, iSumm, iExit, cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine npSeti

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npSetr
     &   ( buffer, rvalue, iPrint, iSumm, iExit,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     iPrint, iSumm, iExit, leniw, lenrw, iw(leniw)
      double precision
     &     rvalue, rw(lenrw)

*     ==================================================================
*     npSetr decodes the option contained in  buffer // rvalue.
*     The parameters other than rvalue are as in npSet.
*
*     27 Nov 1991: first version of npSetr.
*     03 Nov 2000: current version.
*     ==================================================================
      integer
     &     ivalue, lenbuf, lencw
      character
     &     cvalue*8
      double precision
     &     rvalxx
      character
     &     key*16
      character
     &     buff72*72
*     ------------------------------------------------------------------
      parameter            (lencw     = 500)
      character          cw(lencw)*8
*     ------------------------------------------------------------------
      write(key, '(1p, e16.8)') rvalue
      lenbuf = len(buffer)
      buff72 = buffer
      buff72(lenbuf+1:lenbuf+16) = key
      rvalxx = rvalue
      call s3opt
     &   ( .true., buff72, key, cvalue, ivalue, rvalxx,
     &     iPrint, iSumm, iExit, cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine npSetr

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      integer function npGet
     &   ( buffer, iExit, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     iExit, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)

*     ==================================================================
*     npGet  decodes the option contained in  buffer
*     and returns 1 if the option has previously been set, else 0.
*     For example,
*     i = npGet ( 'Maximize', iExit, iw, leniw, rw, lenrw )
*
*     01 Aug 2003: First version of npGet.  Needed because
*                  npGetc, npGeti, npGetr were not well defined
*                  for strings that had no numerical value.
*     18 Feb 2003: Current version of npGet.
*     ==================================================================
      integer
     &     ivalue
      double precision
     &     rvalue
      character
     &     cvalue*8, key*16
*     ------------------------------------------------------------------
      parameter            (lencw     = 500)
      character          cw(lencw)*8
*     ------------------------------------------------------------------
      call s3opt
     &   ( .false., buffer, key, cvalue, ivalue, rvalue,
     &     0, 0, iExit, cw, lencw, iw, leniw, rw, lenrw )

      npGet  = ivalue

      end ! integer function npGet

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npGeti
     &   ( buffer, ivalue, iExit, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     ivalue, iExit, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)

*     ==================================================================
*     npGeti gets the value of the option contained in  buffer.
*     The parameters other than ivalue are as in npSet.
*
*     17 May 1998: first version of npGeti.
*     03 Nov 2000: current version.
*     ==================================================================
      integer
     &     lencw
      double precision
     &     rvalue
      character
     &     key*16
      character
     &     cvalue*8
*     ------------------------------------------------------------------
      parameter            (lencw     = 500)
      character          cw(lencw)*8
*     ------------------------------------------------------------------
      call s3opt
     &   ( .false., buffer, key, cvalue, ivalue, rvalue,
     &     0, 0, iExit, cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine npGeti

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npGetr
     &   ( buffer, rvalue,
     &     iExit, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     iExit, leniw, lenrw, iw(leniw)
      double precision
     &     rvalue, rw(lenrw)

*     ==================================================================
*     npGetr gets the value of the option contained in  buffer.
*     The parameters other than rvalue are as in npSet.
*
*     17 May 1998: first version of npGetr.
*     03 Nov 2000: current version.
*     ==================================================================
      integer
     &     ivalue, lencw
      character
     &     key*16
      character
     &     cvalue*8
*     ------------------------------------------------------------------
      parameter            (lencw     = 500)
      character          cw(lencw)*8
*     ------------------------------------------------------------------
      call s3opt
     &   ( .false., buffer, key, cvalue, ivalue, rvalue,
     &     0, 0, iExit, cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine npGetr

