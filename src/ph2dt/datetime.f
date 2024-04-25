c Current date and time, as a character string
c (This version merely returns a string of blanks)

	subroutine datetime(dattim)

	implicit none

c	Parameters:
	character	dattim*20

      dattim = '                    '
      return
      end ! of subroutine datetime
