c Convert from local Cartesian coordinates to latitude & longitude

	subroutine redist(xkm, ykm, xlat, xlon)

	implicit none

c	Parameters:
	real		xkm, ykm	! km
	doubleprecision	xlat, xlon	! Degrees

	include "geocoord.inc"

c	Local variables:
	real		bcl
	integer		lat, lon
	doubleprecision lat1, lat2, lat3, clat1
	real		p, q
	real		x, y
	real		yp
	real		xx, yy

      xx=xkm
      yy=ykm

c     Rotate coordinates anticlockwise back
      y = yy*cost-xx*sint
      x = yy*sint+xx*cost
      if (abs(aa).lt.0.0000001) goto 900
      q = y/aa
      lat = (q+olat)/60.
      xlat = q+olat - 60.*lat
      yp = 60.*lat+xlat
      lat1 = datan(rlatc*dtan(yp*rad/60.0))
      lat2 = datan(rlatc*dtan(olat*rad/60.))
      lat3 = (lat1+lat2)/2.
      clat1 = dcos(lat3)
      bcl = bb*clat1
      if (abs(bcl).lt.0.000001) goto 900
      p = x/(bb*clat1)
      lon = (p+olon)/60.
      xlon = p+olon - 60.*lon
      xlat = lat+xlat/60.
      xlon = lon+xlon/60.
      return
  900 write(6,1000) aa,bb,clat1
 1000 format(/,2x,' subr. redist: aa=',f10.5,2x,'bb=',f10.5,2x,
     1'cos(lat1)=',f10.7,5x,'division by zero, run stops here',/)
      stop'redist>>> division by zero!!'
      end
