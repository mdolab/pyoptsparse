C Helper subroutine to close files in the Fortran world
      subroutine closeunit(unitnum)

      integer unitnum
Cf2py intent(in) unitnum

      close(unitnum)

      return
      end
