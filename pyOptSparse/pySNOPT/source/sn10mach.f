*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  sn10mach.f                    Machine dependent routines
*
*     s1cpu                             =>  Timing routine
*     s1eps    s1flmx   s1flmn          =>  Floating-point arithmetic
*     s1intmx                           =>  largest integer
*     s1inpt   s1outpt                  =>  standard input/output
*     s1file                            =>  Default File types
*     s1clos   s1envt   s1open   s1page =>  Bibs n bobs
*
*     22 Jun 2004: s1cpu: Added GAMS clock gfclck().
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s1cpu
     &   ( mode, time )

      integer
     &     mode
      real
     &     time

*     ------------------------------------------------------------------
*     s1cpu is a machine-dependent routine to return time = cpu time
*     in seconds, so that 2 consecutive calls will indicate the
*     time difference of operations between the 2 calls.
*     The parameter 'mode' indicates what function should be done
*     to the timer.  This allows necessary initialization for certain
*     machines.
*     mode =  1  indicates initialization,
*     mode =  0  indicates normal use,
*     mode = -1  indicates stop the timer.
*
*     1988:  Used Irv Lustig's approach...
*     On DEC VAX/VMS systems we need to call the correct library
*     routine to get the timer statistics.  These statistics are
*     found by using the times() function in the VAX C Runtime library.
*     To use this version of s1cpu, one must create an options file
*     called  vmsc.opt  with the line
*        SYS$LIBRARY:VAXCRTL/SHARE
*     in it.   Then link using the usual command and append ,vmsc/opt
*     to the end of the line.  The name vmsc can be anything.
*
*     02 Apr 1993: Went back to VMS Fortran routines to avoid linking
*                  to the C library.  (On DEC AXP, the C runtime lib
*                  appears to be translated from the VAX executable,
*                  and therefore requires linking with /NONATIVE,
*                  which possibly adds a small overhead to all
*                  subroutine calls.
*     21 Oct 1999: Timer for WinNT with DEC F90 compiler.  Code provided
*                  by Thomas Kronseder.
*     22 Jun 2004: (At GAMS) Added GAMS clock gfclck().
*                  1. Uncomment two separate lines:
*                        external gfclck
*                           ...
*                        time   = gfclck()
*                  2. Following Unix (Sun, SGI, etc), comment out
*                        time   = etime ( tarray )
*     ------------------------------------------------------------------

*-->  WinNT with DEC F90
*-->  USE DFPORT, ONLY: RTC
*-->  REAL*8             dectime, decinit
*-->  SAVE               decinit

*-->  DEC OpenVMS with Fortran runtime library.
*-->  external        lib$init_timer, lib$stat_timer, lib$free_timer
*-->  integer         itimad, istatu, idata
*-->  save            itimad

*-->  DEC VAX/VMS with C runtime library.
*-->  integer            itimad(4)

*-->  PC Lahey Fortran
*-->  integer            itimad(4)

*-->  AIX
*-->  integer          mclock
*-->  intrinsic        real

*-->  Linux, Unix (SGI, Sun, DECstation)
ccc      external           etime
      real               tarray(2)

*-->  GAMS clock
*     external           gfclck
*     double precision   gfclck

      if (mode .eq. 1) then
*        ---------------------------------------------------------------
*        Initialize.
*        ---------------------------------------------------------------
         time   = 0.0
*-->     DEC OpenVMS with Fortran library.
*-->     istatu = lib$init_timer( itimad )
*-->     if (.not. istatu) call lib$signal( %val(istatu) )

*-->     WinNT with DEC F90
*-->     decinit = RTC()

      else if (mode .eq. 0) then
*        ---------------------------------------------------------------
*        Normal call.
*        Return current timer value here.
*        ---------------------------------------------------------------

*-->     DEC OpenVMS with Fortran library.
*-->     istatu returns the number of  centiseconds.
*-->     istatu = lib$stat_timer( 2, idata, itimad )
*-->     if (.not. istatu) call lib$signal( %val(istatu) )
*-->     time   = idata
*-->     time   = time * 0.01d+0

*-->     DEC VAX/VMS with C library.
*-->     itimad(1) returns the number of  centiseconds.
*-->     call times ( itimad )
*-->     time   = itimad(1)
*-->     time   = time * 0.01d+0

*-->     PC Lahey Fortran, itimad(1) returns the number of  centiseconds.
*-->     call timer ( itimad )
*-->     time   = itimad(1)
*-->     time   = time * 0.01d+0

*-->     On AIX, mclock returns hundredths of a second
*-->     time = real( mclock( ) ) / 100.0

*-->     On Unix (SGI Irix, Sun Solaris), etime returns seconds.
*-->     Linux g77, NagWare f90, Absoft f77
         time   = etime ( tarray )

*-->     NagWare f95, use f95 intrinsic.
*-->     call cpu_time(time)

*-->     On UNIX (NeXTstation M68040), using routine in ftime.c
*-->     call ftime(time)

*-->     WinNT with DEC F90
*-->     Time in secs since 00:00:00 GMT Jan 1st, 1970.
*-->     Bias must be subtracted to give adequate precision for real*4
*-->     dectime = RTC() - decinit
*-->     time    = real(dectime)

*-->     GAMS clock
*        time   = gfclck()

*-->     On other machines, to forget about timing, just say
*-->     time   = -1.0

      else if (mode .eq. -1) then
*        ---------------------------------------------------------------
*        Stop the clock.
*        ---------------------------------------------------------------
         time   = 0.0

*-->     DEC OpenVMS with Fortran library.
*-->     istatu = lib$free_timer( itimad )
*-->     if (.not. istatu) call lib$signal( %val(istatu) )
      end if

      end ! subroutine s1cpu

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      double precision   function s1eps ( )
*
*     Compute the machine precision.
*     IEEE Floating point double precision.
*
      integer
     &     nbase, ndigit
      double precision
     &     base, u
*     -----------------------------------------------------------------
      nbase  = 2
      ndigit = 53
      base   = nbase
      u      = base**(- ndigit)
      s1eps  = 2.0d+0*u

      end ! function s1eps

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      double precision   function s1flmx( )
*
*     IEEE Floating point double precision.
*
      s1flmx = 1.7977d+307

      end ! function s1flmx

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      double precision   function s1flmn( )
*
*     IEEE Floating point double precision.
*
      s1flmn = 2.2251d-308

      end ! function s1flmn

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      integer            function s1intmx( )
*
*     2**31-1, the largest positive 32-bit integer
*
      s1intmx = 2147483647

      end ! function s1intmx

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      integer            function s1inpt( )
*
*     Fortran standard input
*
      s1inpt = 5

      end ! function s1inpt

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      integer            function s1outpt( )
*
*     Fortran standard output
*
      s1outpt = 6

      end ! function s1outpt

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s1file
     &   ( Task, iw, leniw )

      implicit
     &     none
      integer
     &     task, leniw, iw(leniw)

*     ------------------------------------------------------------------
*     s1file  is a machine-dependent routine for opening various files.
*     It calls s1open (which is also machine-dependent).
*
*     SNOPT uses sequential files only
*     and does not need to read and write to the same file.
*
*     iPrint, iSumm  are defined in snInit
*     iSpecs         is  defined in snSpec
*     iStdi          is defined in s1inpt, but loaded into memory here.
*
*     iStdi and iPrint have the following use:
*        Input  files (MPS, Old Basis, Insert, Load)
*        are rewound after being read,
*        but not if they are the same as  iStdi.
*        Output files (Backup, New Basis, Punch, Dump, Solution, Report)
*        are rewound after being written,
*        but not if they are the same as  iPrint.
*
*     iStdi  = (conceptually) the Keyboard that can't be rewound.
*              SNOPT does not use this file, so there is no 'open'.
*     iSumm  = the SUMMARY file.  Sometimes this is the Terminal.
*              If so, it may not need to be opened.
*     iSpecs = the SPECS file, containing one or more problem specs.
*              This file is not rewound after use, because it may
*              contain another SPECS file.
*
*     Here are all the files used by SNOPT.
*     The associated Index is passed to s1open
*     and must match the list of names in s1open, if that routine
*     uses method = 1.
*
*        Unit    Index    Description       Status
*        iSpecs     1     Specs     file    In
*        iPrint     2     Print     file    Out
*        iSumm      3     Summary   file    Out
*        iMPS       4     MPS       file    In
*        iOldB      5     Old Basis file    In
*        iInsrt     6     Insert    file    In
*        iLoadB     7     Load      file    In
*        iBack      8     Backup    file    Out
*        iNewB      9     New Basis file    Out
*        iPnch     10     Punch     file    Out
*        iDump     11     Dump      file    Out
*        iSoln     12     Solution  file    Out
*        iReprt    13     Report    file    Out
*        iStdi            Not opened, but used as described above
*
*     15 Nov 1991: First version based on Minos 5.4 routine mifile.
*     03 Oct 1997: Modified to comply with Minos 5.5.
*     31 Jul 2003: snPRNT adopted.
*     22 Nov 2003: Current version.
*     ------------------------------------------------------------------
      external
     &     s1inpt
      integer
     &     iSpecs, iPrint, iSumm, iBack, iDump, iInsrt, iLoadB, iNewB,
     &     iOldB, iPnch, iReprt, iSoln, s1inpt
      character
     &     str*80
*     ------------------------------------------------------------------
      integer            DefltF,     OpenF,      StdIn
      parameter         (DefltF = 0, OpenF  = 1, StdIn  = 2)

      integer            iStdi, iMPS, iPrinx, iSummx
      parameter         (iStdi     =   9) ! Standard input
      parameter         (iMPS      = 123) ! MPS file
      parameter         (iPrinx    = 228) ! Global value of iPrint
      parameter         (iSummx    = 229) ! Global value of iSumm
*     ------------------------------------------------------------------
      iSpecs    = iw( 11) ! Specs (options) file
      iPrint    = iw( 12) ! Print file
      iSumm     = iw( 13) ! Summary file

      iBack     = iw(120) ! backup file
      iDump     = iw(121) ! dump file
      iLoadB    = iw(122) ! load file
      iNewB     = iw(124) ! new basis file
      iInsrt    = iw(125) ! insert file
      iOldB     = iw(126) ! old basis file
      iPnch     = iw(127) ! punch file
      iReprt    = iw(130) ! Report file
      iSoln     = iw(131) ! Solution file

*     ------------------------------------------------------------------
*-->  Machine dependency.
*     Set iStdi = some input unit number that should not be rewound.
*     Set iStdi = 0 if this is irrelevant.
*     ------------------------------------------------------------------
      iw(iStdi) = s1inpt ( )

      if (Task .eq. StdIn ) then

*        Relax, do nothing

      else if (Task .eq. DefltF) then
*        ---------------------------------------------------------------
*        Task = Default: Open the Specs, Print and Summary files.
*        iSpecs          remains the same throughout the run.
*        iPrint, iSumm   may be altered by the SPECS file.  They may
*                        need to be opened by both Tasks Deflt and Open.
*        ---------------------------------------------------------------
         iw(iPrinx) = iPrint
         iw(iSummx) = iSumm
         call s1open( iSpecs, 1, 'IN ' )
         call s1open( iPrint, 2, 'OUT' )
         call s1open( iSumm , 3, 'OUT' )

      else if (Task .eq. OpenF) then
*        ---------------------------------------------------------------
*        Task = OpenF: Open files mentioned in the SPECS file just read.
*        Input files are opened first.  Only one basis file is needed.
*        ---------------------------------------------------------------
         if (iw(iMPS) .le. 0     )              iw(iMPS) = iSpecs
         if (iw(iMPS) .ne. iSpecs) call s1open( iw(iMPS),  4, 'IN ' )

         if      (iOldB  .gt. 0) then
            call s1open( iOldB   ,  5, 'IN ' )
         else if (iInsrt .gt. 0) then
            call s1open( iInsrt  ,  6, 'IN ' )
         else if (iLoadB .gt. 0) then
            call s1open( iLoadB  ,  7, 'IN ' )
         end if
            call s1open( iBack   ,  8, 'OUT' )
            call s1open( iNewB   ,  9, 'OUT' )
            call s1open( iPnch   , 10, 'OUT' )
            call s1open( iDump   , 11, 'OUT' )
            call s1open( iSoln   , 12, 'OUT' )
            call s1open( iReprt  , 13, 'OUT' )

*        Open new Print or Summary files if they were altered
*        by the Specs file.

         if (iPrint .ne. iw(iPrinx)) call s1open( iPrint, 2, 'OUT' )
         if (iSumm  .ne. iw(iSummx)) call s1open( iSumm , 3, 'OUT' )
      end if

*     Check that output files are different from Specs or MPS.

      if (iSpecs .gt. 0) then
         if (iBack  .eq. iSpecs  ) then
            write(str, 1000) 'Backup'
            call snPRNT( 1, str, iw, leniw )
         end if

         if (iNewB  .eq. iSpecs  ) then
            write(str, 1000) 'New Basis'
            call snPRNT( 1, str, iw, leniw )
         end if

         if (iPnch  .eq. iSpecs  ) then
            write(str, 1000) 'Punch'
            call snPRNT( 1, str, iw, leniw )
         end if

         if (iDump  .eq. iSpecs  ) then
            write(str, 1000) 'Dump'
            call snPRNT( 1, str, iw, leniw )
         end if

         if (iSoln  .eq. iSpecs  ) then
            write(str, 1000) 'Solution'
            call snPRNT( 1, str, iw, leniw )
         end if

         if (iReprt .eq. iSpecs  ) then
            write(str, 1000) 'Report'
            call snPRNT( 1, str, iw, leniw )
         end if
      end if

      if (iw(iMPS) .gt. 0) then
         if (iBack  .eq. iw(iMPS)) then
            write(str, 2000) 'Backup'
            call snPRNT( 1, str, iw, leniw )
         end if

         if (iNewB  .eq. iw(iMPS)) then
            write(str, 2000) 'New Basis'
            call snPRNT( 1, str, iw, leniw )
         end if

         if (iPnch  .eq. iw(iMPS)) then
            write(str, 2000) 'Punch'
            call snPRNT( 1, str, iw, leniw )
         end if

         if (iDump  .eq. iw(iMPS)) then
            write(str, 2000) 'Dump'
            call snPRNT( 1, str, iw, leniw )
         end if

         if (iSoln  .eq. iw(iMPS)) then
            write(str, 2000) 'Solution'
            call snPRNT( 1, str, iw, leniw )
         end if

         if (iReprt .eq. iw(iMPS)) then
            write(str, 2000) 'Report'
            call snPRNT( 1, str, iw, leniw )
         end if
      end if

      return

 1000 format(' ===>  Warning: the Specs file and ', a,
     &       ' file are on the same unit')
 2000 format(' ===>  Warning: the  MPS  file and ', a,
     &       ' file are on the same unit')

      end ! subroutine s1file

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s1clos
     &   ( lun )

      implicit
     &     none
      integer
     &     lun

*     ==================================================================
*     s1clos  closes the file with logical unit number lun.
*     This version is trivial and so far is not even used by SNOPT.
*     Perhaps some implementations will need something fancier.
*     ==================================================================

      close( lun )

      end ! subroutine s1clos

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s1envt
     &   ( mode, iw, leniw )

      implicit
     &     none
      integer
     &     mode, leniw, iw(leniw)

*     ==================================================================
*     s1envt specifies the environment within which SNOPT is being used.
*
*     When mode = 0, information about the environment should be
*     initialized.
*
*     iPage1 says whether new pages are ever wanted on file iPrint.
*     iPage2 says whether new pages are ever wanted on file iSumm.
*
*     When mode is in the range 1 to 99, each environment does its
*     own thing.
*
*     The only environment at present is:
*     ALONE:
*     This means SNOPT is in stand-alone mode---the normal case.
*     Nothing special is done.
*
*     16 Sep 1987.
*     ==================================================================
      integer            IALONE, iPage1, iPage2
      parameter         (iALONE    = 238) ! > 0    =>  stand-alone
      parameter         (iPage1    = 241) ! > 0    =>  Page 1
      parameter         (iPage2    = 242) ! > 0    =>  Page 2
*     ------------------------------------------------------------------
      if (mode .le. 0) then
*        ---------------------------------------------------------------
*        mode = 0.    Initialize.
*        ---------------------------------------------------------------
         iw(iALONE) = 1
         iw(iPage1) = 1
         iw(iPage2) = 0

      else

*        Relax, do nothing in this version.

      end if

      end ! subroutine s1envt

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s1open
     &   ( lun, index, state )

      implicit
     &     none
      integer
     &     lun, index
      character
     &     state*3

*     ------------------------------------------------------------------
*     s1open  is a machine-dependent routine.
*     In principal it opens a file with logical unit number lun
*     and positions it at the beginning.
*
*     Input files are treated that way.
*     An input file is opened with status='OLD'
*     (and F77 will terminate with an error if the file doesn't exist).
*
*     Output files are more machine-dependent.
*     With status='NEW', F77 would terminate if the file DID exist.
*     With status='UNKNOWN', existing files would be overwritten.
*     This is normal with Unix, but on systems that have file
*     version numbers (e.g. DEC OpenVMS), a new version is preferable.
*     It is then better not to open the file at all, but let a new
*     version be created by the first "write".
*
*     Nothing happens if
*     1. lun <= 0.  Saves us from testing lun before calling s1open.
*
*     2. Unit lun is already open.  Helps applications that call
*        snopt -- they can open files themselves if they want to.
*
*     3. lun = screen, where  screen  is a local machine-dependent
*        variable, typically 6.  It seems inadvisable to do an OPEN
*        on a file that is predefined to be an interactive screen.
*        With Unix on DECstations, open(6, file='snopt.sum') sends
*        output to file snopt.sum, not to the screen.
*
*
*     lun     (input) is the unit number.
*
*     index   (input) points to one of the hardwired names below.
*             Used only if method = 1.
*
*     state   (input) is 'IN ' or 'OUT', indicating whether the file
*             is to be input or output.
*
*     15 Jul 1989: First version, follows some of the advice offered
*                  by David Gay, Bell Laboratories.
*     -- --- 1990: Added parameter "state".
*     03 Feb 1994: Added parameter "index" to help when method = 1.
*                  Local variable "method" must be set here to select
*                  various methods for naming files.
*                  Chris Jaensch, IFR Stuttgart, recommends not opening
*                  input files if they are already open.  This version
*                  ignores all open files (input or output), assuming
*                  that they are taken care of by the calling program.
*     13 May 1994: Local variable "screen" used to avoid opening screen.
*     20 May 2001: Method 4 is now the default. The default getfnm uses
*                  method 2.
*     ------------------------------------------------------------------
      logical
     &     input, uopen
      integer
     &     last, method, screen
      character
     &     filnam*100
*     ------------------------------------------------------------------
*-->  Machine dependency.
*     names(*) is needed if method = 1 below.
*     Make sure "character*n" sets n big enough below.
*     It doesn't matter if it is bigger than necessary, since
*     "open( lun, file=name )" allows name to have trailing blanks.
*     ------------------------------------------------------------------
      character          names(13)*9
      data   names( 1) /'snopt.spc'/
      data   names( 2) /'snopt.prn'/
      data   names( 3) /'snopt.sum'/
      data   names( 4) /'snopt.mps'/
      data   names( 5) /'snopt.olb'/
      data   names( 6) /'snopt.ins'/
      data   names( 7) /'snopt.lod'/
      data   names( 8) /'snopt.bak'/
      data   names( 9) /'snopt.nwb'/
      data   names(10) /'snopt.pun'/
      data   names(11) /'snopt.dmp'/
      data   names(12) /'snopt.sol'/
      data   names(13) /'snopt.rpt'/

*     ------------------------------------------------------------------
*-->  Machine-dependency.
*     Set "method" to suit your operating system.
*     It determines how a file name is assigned to unit number "lun".
*     Typically,
*     method = 1 for fixed file names (e.g. PCs under DOS).
*     method = 2 for names like fort.7, fort.15 (e.g. Unix on Sun, SGI).
*     method = 3 for names like FTN07 , FTN15   (e.g. Unix on HP).
*     method = 4 if names are assigned by an external routine getfnm
*                (default getfnm sets names as in method 2).
*     method = 5 if explicit file names are not needed (e.g. OpenVMS).
*     See more comments below.
*
*     Set "screen" to a unit number that never needs to be opened.
*     (Typically, screen = 6.  If unknown, set screen = 0.)
*     ------------------------------------------------------------------
      method = 4                ! Default value for SNOPT
      screen = 6

*     ------------------------------------------------------------------
*     Quit if lun<=0 or lun = iscreen or unit lun is already open.
*     ------------------------------------------------------------------
      if (lun .le. 0     ) go to 900
      if (lun .eq. screen) go to 900
      inquire( lun, opened=uopen )
      if (     uopen     ) go to 900

*     ------------------------------------------------------------------
*     Open file lun by a specified machine-dependent method.
*     Note that 'UNKNOWN' is equivalent to trying first with 'OLD',
*     and then with 'NEW' if the file doesn't exist.
*     ------------------------------------------------------------------
      input  = state .eq. 'IN '  .or.  state .eq. 'in '

      if (method .eq. 1) then
*        ---------------------------------------------------------------
*        Hardwired filename.
*        We use "index" to get it from names(*) above.
*        Typical machines: IBM PC under DOS.
*        ---------------------------------------------------------------
         if ( input ) then
            open( lun, file=names(index), status='OLD' )
         else
            open( lun, file=names(index), status='UNKNOWN' )
         end if

      else if (method .eq. 2) then
*        ---------------------------------------------------------------
*        Construct a name like fort.7 or fort.15.
*        (This approach suggested by Chris Jaensch, IFR Stuttgart.)
*        Typical machines:  Unix on Sun, SGI, DEC.
*        ---------------------------------------------------------------
         if (lun .le. 9) then
             write(filnam, '(a,i1)') 'fort.', lun
         else
             write(filnam, '(a,i2)') 'fort.', lun
         endif

         if ( input ) then
            open( lun, file=filnam, status='OLD' )
         else
            open( lun, file=filnam, status='UNKNOWN' )
         end if

      else if (method .eq. 3) then
*        ---------------------------------------------------------------
*        Construct a name like FTN07 or FTN15.
*        Typical machines:  Unix on HP.
*        ---------------------------------------------------------------
         if (lun .le. 9) then
             write(filnam, '(a,i1)') 'FTN0', lun
         else
             write(filnam, '(a,i2)') 'FTN', lun
         endif

         if ( input ) then
            open( lun, file=filnam, status='OLD' )
         else
            open( lun, file=filnam, status='UNKNOWN' )
         end if

      else if (method .eq. 4) then
*        ---------------------------------------------------------------
*        Assume some routine "getfnm" will provide a name at run-time.
*        (This approach is used by David Gay, Bell Labs.)
*        The default getfnm.f provided in sn12ampl.f constructs names
*        as in method 2.
*        ---------------------------------------------------------------
         call getfnm(lun, filnam, last)
         if ( input ) then
            open( lun, file=filnam(1:last), status='OLD' )
         else
            open( lun, file=filnam(1:last), status='UNKNOWN' )
         end if

      else if (method .eq. 5) then
*        ---------------------------------------------------------------
*        Explicit file names are not needed for the open statement.
*        The operating system uses a default name
*        (e.g. fort.7 or for007)
*        or has already assigned a name to this unit
*        (e.g. via a script or command file).
*        Typical machines:  Unix,
*                           DEC OpenVMS,
*                           IBM Mainframes.
*        ---------------------------------------------------------------
         if ( input ) then
            open( lun, status='OLD' )
         else
*           Let the first "write" do it.
         end if
      end if

*     ------------------------------------------------------------------
*     Rewind input files.
*     (Some systems position existing files at the end
*     rather than the beginning.)
*     err=900 covers files that have not yet been opened.
*     ------------------------------------------------------------------
      if ( input ) then
         rewind( lun, err=900 )
      end if

  900 return

      end ! subroutine s1open

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s1page
     &   ( mode, iw, leniw )

      implicit
     &     none
      integer
     &     mode, leniw, iw(leniw)

*     ------------------------------------------------------------------
*     s1page is an installation-dependent routine.  It is called at
*     points where some users might want output to files iPrint or iSumm
*     to begin on a new page.
*
*     iPage1 and iPage2 have already been set by s1envt.
*     If they are true, a page eject and a blank line are output.
*     Otherwise, just a blank line is output.
*
*     If mode = 0  and Summary level = 0, nothing is output to the
*                  Summary file.  At present, this is so s8log
*                  will print just one line per major iteration, with
*                  no blank line in between.
*     If mode = 1, just the page control is relevant.
*     If mode = 2, SNOPT has encountered an error condition.
*                  At the moment, this case is treated the same as
*                  mode = 1.
*
*     15 Nov 1991: First version based on Minos 5.4 routine m1page.
*     31 Jul 2003: snPRNT adopted.
*     31 Jul 2003: Current version of s1page.
*     ==================================================================
      integer
     &     iPage1, iPage2
*     ------------------------------------------------------------------
      iPage1    = iw(241) ! > 0    =>  Page 1
      iPage2    = iw(242) ! > 0    =>  Page 2

      if (iPage1 .gt. 0) call snPRNT( 1, '1', iw, leniw )
                         call snPRNT( 1, ' ', iw, leniw )

      if (iPage2 .gt. 0) call snPRNT( 2, '1', iw, leniw )
      if (mode   .ne. 0) call snPRNT( 2, ' ', iw, leniw )

      end ! subroutine s1page
