*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  sn03prnt.f  (default versions of snPRNT and snREAD)
*
*     snPRNT   snREAD
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snPRNT
     &   ( mode, string, iw, leniw )

      implicit
     &     none
      character*(*)
     &     string
      integer
     &     mode, leniw, iw(leniw)

*     ==================================================================
*     snPRNT  prints a trimmed form of "string" on various files.
*     If mode = 0,      nothing is output.
*     If mode = 1,      string is output to iPrint.
*     If mode = 2,      string is output to iSumm.
*     If mode = 3 or 4, string is output to iPrint and iSumm.
*     If mode = 4,      string is output to iStdo (standard output)
*                       if iPrint and iSumm are both zero.  This mode
*                       is intended for error messages.
*     If mode = 5,      string is output to iStdo (standard output)
*                       This mode is to be used when the elements of
*                       the integer work array iw cannot be trusted.
*
*     mode 11-15 are the same as mode 1-5 with blank line before output.
*
*     If mode > 15 then nothing is printed unless  lvlSys > 0.
*     mode 21-25 are the same as mode 1-5
*     mode 31-35 are the same as mode 11-15
*
*     25 Sep 2002: First version of snPRNT.
*     31 Jul 2003: mode 11-14 added.  form introduced.
*     27 Dec 2003: mode 5 added to allow printing before iw is set.
*     12 Mar 2004: s1trim called to trim the string.
*     22 Jun 2004: System printing option added.
*     22 Jun 2004: Current version of snPRNT.
*     ==================================================================
      external
     &     s1outpt
      integer
     &     iPrint, iSumm, iStdo, length, lvlSys, m, s1outpt
      character
     &     form*4, form1*4, form2*4
*     ------------------------------------------------------------------
      parameter   (form1 = '( a)')
      parameter   (form2 = '(/a)')
*     ------------------------------------------------------------------
      lvlSys    = iw( 71) ! > 0   => print system info

      m = 0
      if (mode .le.  0) then
!        Relax
      else if (mode   .lt. 10) then
         m    = mode
         form = form1
      else if (mode   .lt. 20) then ! Blank line first
         m    = mode - 10
         form = form2
      else if (lvlSys .gt.  0) then ! Print system Info
         if (mode .lt. 30) then
            m    = mode - 20
            form = form1
         else
            m    = mode - 30
            form = form2
         end if
      end if

      if (m .gt. 0) then

!        length = len_trim(string)     ! An F90 intrinsic
         call s1trim( string, length ) ! The F77 equivalent

         if (m .eq. 5) then
            iStdo  = s1outpt()
            if (iStdo .gt. 0) write(iStdo, form) string(1:length)
         else
            iStdo  = iw( 10) ! Standard output
            iPrint = iw( 12) ! Print file
            iSumm  = iw( 13) ! Summary file

            if (m .eq. 1  .or.  m .ge. 3) then
               if (iPrint .gt. 0) write(iPrint, form) string(1:length)
            end if

            if (m .eq. 2  .or.  m .ge. 3) then
               if (iSumm  .gt. 0) write(iSumm , form) string(1:length)
            end if

            if (m .eq. 4) then
               if (iPrint .le. 0  .and.  iSumm .le. 0) then
                  if (iStdo .gt. 0) write(iStdo, form) string(1:length)
               end if
            end if
         end if
      end if

      end ! subroutine snPRNT

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snREAD
     &   ( unitno, string, nchar, endfile )

      implicit
     &     none
      character*(*)
     &     string
      integer
     &     endfile, nchar, unitno
!     ==================================================================
!     snREAD reads a string of length nchar (including trailing blanks)
!     from the file with logical unit number  unitno.
!
!     Restriction: 0 < nChar < 1000
!
!     On exit:
!       endfile = 0 means that the string was read successfully.
!       endfile = 1 means that an end-of-file was encountered or nchar
!                   lies outside the range  0 < nChar < 1000.
!
!     30 Apr 2006: First version of snREAD.
!     ==================================================================
      character
     &     frmt*6
!     ------------------------------------------------------------------
      frmt    = '      '

      if (nchar .ge. 1  .and.  nchar .le. 999) then
         if      (nchar .lt.  10) then
            write(frmt, '(a2,i1,a1)') '(a', nchar, ')'
         else if (nchar .lt. 100) then
            write(frmt, '(a2,i2,a1)') '(a', nchar, ')'
         else
            write(frmt, '(a2,i3,a1)') '(a', nchar, ')'
         end if

         endfile = 0
         read  (unitno, frmt, end = 100) string
         return
      end if

  100 endfile = 1

      end ! subroutine snREAD
