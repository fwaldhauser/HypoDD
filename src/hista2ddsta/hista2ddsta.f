c--Filter to convert Hypoinvese station files to the format used by
c  ph2dt and hypoDD.
c--Fred Klein 2/2001

	character line*80, net*2, site*5, sta*7, laststa*7
	character clat*1,clon*1
	laststa=' '

c--read a line, then decode it
2	read (5,'(a)',end=9) line
	read (line,1000,iostat=ios) site,net, latd,xlat,clat, 
c111109     	2 lond,xlon,clon
     2  lond,xlon,clon
	
1000	format (a5,1x,a2,7x, i2,1x,f7.4,a1, i3,1x,f7.4,a1)
	
	if (ios.gt.0) then
	  write (6,*) '*** Bad station format ',line
	  stop '*** Bad station format.'
	end if
	
c--Only write one entry per station for all components
	sta=(net//site)
	if (sta.eq.laststa) goto 2
	
c--correct signs: + north & east, - south & west
	dlat=latd+xlat/60.
	if (clat.eq.'s' .or. clat.eq.'S') dlat=-dlat
	dlon=-lond-xlon/60.
	if (clon.eq.'e' .or. clon.eq.'E') dlon=-dlon
	
c--Write new line
	laststa=sta
	write (6,1001) sta,dlat,dlon
1001	format (a7,f11.6,f12.6)
	goto 2
9	stop
	end
