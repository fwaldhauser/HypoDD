c Set up Cartesian coordinate system by short distance conversion.
c Unrotated coordinate system with pos. x-axis toward WEST
c and pos. y-axis toward NORTH.
c pos. z-axis toward EARTH`s CENTER.

	subroutine setorg(orlat, orlon, rota, ifil)  !22. Okt. 1990

	implicit none

c	Parameters:
	real	orlat, orlon
	real	rota		! Anticlockwise rotation (degrees)
	integer	ifil

	include "geocoord.inc"

	real	arota
	doubleprecision r, lat1, lat2, dela, delb, PHI, BETA

c     O(r)LAT & O(r)LON : origin of cartesian coordinate system
c      north  w e s t

c if orlat and orlon are both set to zero , the Swiss Cartesian
c coordinate system will be used (this system cannot be rotated).
c  For all other cases, orlat and orlon denote the origin of
c  the SDC.

       rad = 0.017453292d0

      if(orlat.eq.0.0.and.orlon.eq.0.0)then
         olat=46.95240  ! BERN North
         olon=-7.439583  ! BERN West
         rotate=0.0
      else
         olat=orlat
         olon=orlon
         rotate=rota*rad
      endif

      olat=olat*60. ! minutes N
      olon=olon*60. ! minutes W

C  NEW ELLIPSOID FOR WHOLE EARTH:   WGS72 == WORLD GEODETIC SYSTEM 1972

C  ALSO SET RLATC ACCORDING TO ORIGIN

      REARTH=6378.135D0
      ELLIP=298.26         ! (flattening)


C CALCULATE RLATC:  CONVERSION FROM GEOGRAPHICAL LAT TO GEOCENTRICAL LAT

      PHI=OLAT*RAD/60.0              !  phi=geogr. lat
      BETA=PHI-DSIN(PHI*2.)/ELLIP    !  beta=geoc. lat
      RLATC=DTAN(BETA)/DTAN(PHI)

C WRITE ELLIPSOIDE CONSTANTS

      if(ifil.gt.0)then
         write(ifil,*)
         write(ifil,*)
         write(ifil,*)'SHORT DISTANCE CONVERSION on ELLIPSOIDE of'//
     &                ' WORLD GEODETIC SYSTEM 1972 (WGS72)'
         write(ifil,*)'=========================================='//
     &                '==================================='
         write(ifil,*)
         write(ifil,'('' Radius at equator (REARTH)= '',f10.5,
     &                ''  km'')') rearth
         write(ifil,'(''   1. / (ellipticity)      = '',f10.3)') ellip
         write(ifil,*)
         write(ifil,*)'Origin of cartesian coordinates [degrees]:'
         if(orlat.eq.0.0.and.orlon.eq.0.0)then
            write(ifil,*)
            write(ifil,*)' SWISS COORDINATE SYSTEM (we have to be'
            write(ifil,*)'                               special)'
            write(ifil,*)' (Origin = city of BERNE, Switzerland)'
            write(ifil,*)
            write(ifil,*) ' no rotation of grid, pos. y-axis toward N'
            write(ifil,*) '                      pos. x-axis toward E'
            write(ifil,*)
         else
           write(ifil,'(1x,f12.7,'' N'',5x,f12.7,'' W'')')
     &               olat/60.,olon/60.

           write(ifil,*)
           write(ifil,*) ' without rotation of grid, '
           write(ifil,*) '             pos. x-axis toward WEST'
           write(ifil,*) '             pos. y-axis toward NORTH'
           write(ifil,*)
           write(ifil,*) ' Rotation of y-axis from North anticlock-'
           write(ifil,*) ' wise with pos. angle given in degrees'
           write(ifil,*)
           IF(ROTA.GE.0.) THEN
             write(ifil,*) ' Rotation of grid anticlockwise by'
             write(ifil,*)'  ', rota,' degrees'
             write(ifil,*)
           ELSE
             write(ifil,*) ' Rotation of grid clockwise by'
             arota=-1.*rota
             write(ifil,*)'  ', arota,' degrees'
             write(ifil,*)
           ENDIF
         endif
      endif

c   calculate aa &  bb
c   length of one minute of lat and lon in km

      lat1 = datan(rlatc*dtan(olat*rad/60.0))       ! geoc. lat for OLAT
      lat2 = datan(rlatc*dtan((olat+1.)*rad/60.0))  ! geoc. lat for (OLAT+1min)
      dela = lat2 - lat1
      r = rearth*(1.0 - (dsin(lat1)**2)/ellip)      ! kugelradius fuer lat=OLAT
      aa = r*dela   ! aa = 1 min geogr. lat
      delb = dacos(dsin(lat1)**2 + dcos(rad/60.)*dcos(lat1)**2)
      bc=r*delb     ! bc = 1 min geogr. lon
      bb = r*delb/dcos(lat1)
      if(ifil.gt.0)then
         write(ifil,'('' Radius of sphere at OLAT = '',f10.3,'' km'')')r
         write(ifil,*)
         write(ifil,*)'Conversion of GEOGRAPHICAL LATITUDE to '//
     &                'GEOCENTRICAL LATITUDE:'
         write(ifil,*)'RLATC = TAN(GEOCENTR.LAT) / TAN(GEOGRAPH.LAT)'
         write(ifil,'(1x,''RLATC = '',f12.8)') rlatc
         write(ifil,*)
         write(ifil,4) aa, bc
 4       format (10x,'Short distance conversions:',/,
     &           10x,'one min lat ~ ', f7.4,' km ',/,
     &           10x,'one min lon ~ ', f7.4,' km ',/)
         write(ifil,*)
         write(ifil,*)
      endif

c***  convert coordinates with rotation cosines (stored in Common)
      sint=sin(rotate)
      cost=cos(rotate)

      return
      end
