c Convert calendar date & time to minutes

	integer function juliam(iyr, imo, idy, ihr, imn)

	implicit none

c	Parameters:
	integer	iyr, imo, idy, ihr, imn		! (input)
						! iyr < 4000 for 32-bit int

c	Local variables:
	integer	kl
	integer	kmo(12)
	integer	ky, km, kd
	integer	ky0
	integer	ky1
	integer	ky4
	integer	l
	integer	leap

	data kmo/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/
	data leap/1/

      ky = iyr
      km = imo
      kd = idy
      if (km.le.0) km = 1
      juliam = 365*ky
      kd = kmo(km)+kd
      ky4 = ky/4
      ky1 = ky/100
      ky0 = ky/1000
      kl = leap*(ky4-ky1+ky0)
      l = 0
      if (ky4*4.eq.ky.and.(ky1*100.ne.ky.or.ky0*1000.eq.ky)) l = leap
      if (l.ne.0.and.km.lt.3) kl = kl-leap
      juliam = juliam+kd+kl
      juliam = juliam*24+ihr
      juliam = juliam*60+imn
      return
      end ! of integer function juliam
