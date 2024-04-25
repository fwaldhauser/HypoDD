c Compute calendar date from count of minutes (since when?)

	subroutine datum(itf, iyr, imo, idy, ihr, imn)

	implicit none

c	Parameters:
	integer	itf
	integer	iyr
	integer	imo
	integer	idy
	integer	ihr
	integer	imn

C UMRECHNEN DES DATUMS IN MINUTEN (CF. JULIAM) IN YR-MO-DY-HR-MI
C   (MIT IMN<2**31, JAHR < 4000

c	Local variables:
	integer	i
	integer	iyr4
	integer	iyrh
	integer	iyrt
	integer	kh
	integer	l
	integer	ld
	integer id
	integer k
	integer kmo(12)

	data kmo/31,28,31,30,31,30,31,31,30,31,30,31/

      k = itf/60
      imn = itf-k*60
      kh = k/24
      ihr = k-kh*24
      iyr = kh/365
5     id = kh-iyr*365
      l = 0
      iyr4 = iyr/4
      iyrh = iyr/100
      iyrt = iyr/1000
      ld = iyr4-iyrh+iyrt
      if (iyr4*4.eq.iyr.and.(iyrh*100.ne.iyr.or.iyrt*1000.eq.iyr)) l = 1
      id = id-ld+l
      if (id.gt.0) goto 10
      if (id.eq.0.and.ihr.eq.0.and.imn.eq.0) then
          idy = 0
          imo = 0
          return
      endif
      iyr = iyr-1
      goto 5
10    kmo(2) = 28+l
      do i=1,12
         id = id- kmo(i)
         if (id.le.0) goto 30
      enddo
      i =12
30    idy = id+kmo(i)
      imo = i
      return
      end ! of subr. datum
