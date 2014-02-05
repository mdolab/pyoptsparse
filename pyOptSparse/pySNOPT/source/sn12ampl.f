*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  sn10ampl.f                    Machine dependent routine
*
*     getfnm
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine getfnm( lun, filnam, last )

      implicit
     &     none
      integer
     &     lun, last
      character
     &     filnam*100

*     ------------------------------------------------------------------
*     getfnm  is a machine-dependent routine.
*     In principal it opens a file with logical unit number lun
*     and positions it at the beginning.
*
*     lun     (input) is the unit number.
*
*     16 May 2001: First version, follows some of the advice offered
*                  by David Gay, Bell Laboratories.
*     ------------------------------------------------------------------

      if (lun .le. 9) then
         write(filnam, '(a,i1)') 'fort.', lun
         last = 6
      else
         write(filnam, '(a,i2)') 'fort.', lun
         last = 7
      endif

      end ! subroutine getfnm
