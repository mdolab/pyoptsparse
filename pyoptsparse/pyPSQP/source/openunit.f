C Helper subroutine to open files in the Fortran world
      subroutine openunit(unitnum,filename,filestatus,fileaction,ierror)

      integer unitnum
Cf2py intent(in) unitnum
      character*(*) filename
Cf2py intent(in) filename
      character*(*) filestatus
Cf2py intent(in) filestatus
      character*(*) fileaction
Cf2py intent(in) fileaction
      integer ierror
Cf2py intent(out) ierror
      
      open(unit=unitnum,file=filename,status=filestatus,
     >     access=fileaction,iostat=ierror)

      return
      end

C Helper routine to flush buffers to files
      subroutine pyflush(unitnum)
      
      integer unitnum
      
      call flush(unitnum)

      return
      end
