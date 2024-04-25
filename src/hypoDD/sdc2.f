c Convert coordinates of a point by short distance conversion (SDC)

	subroutine sdc2(x, y, xlat, xlon, i)     ! 22.Okt. 1990

	implicit none

c	Parameters:
	real		x, y
	doubleprecision xlat, xlon
	integer		i		!  1: (x, y) -> (xlat, xlon)
					! -1: (xlat, xlon) -> (x, y)

	include "geocoord.inc"

      if(.not.(i.eq.-1.or.i.eq.1)) stop'SDC>>> specify conversion !!!'

      if(i.eq.1)  call redist(x,y,xlat,xlon)
      if(i.eq.-1) call dist(xlat,xlon,x,y)
      end
