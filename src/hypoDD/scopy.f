c     Copy single precision sx to single precision sy.
c     for i = 0 to n-1, copy  sx(lx+i*incx) to sy(ly+i*incy),
c     where lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is
c     defined in a similar way using incy.

	subroutine scopy(n, sx, incx, sy, incy)

	implicit none

c	Parameters:
	integer	n
	real	sx(*)
	integer	incx
	real	sy(*)
	integer	incy

c	Local variables:
	integer	i
	integer	ix, iy
	integer	m
	integer	mp1
	integer	ns

      if(n.le.0)return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue

c        code for unequal or nonpositive increments.

      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return

c        code for both increments equal to 1


c        clean-up loop so remaining vector length is a multiple of 7.

   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        sy(i) = sx(i)
        sy(i + 1) = sx(i + 1)
        sy(i + 2) = sx(i + 2)
        sy(i + 3) = sx(i + 3)
        sy(i + 4) = sx(i + 4)
        sy(i + 5) = sx(i + 5)
        sy(i + 6) = sx(i + 6)
   50 continue
      return

c        code for equal, positive, nonunit increments.

   60 continue
      ns = n*incx
          do 70 i=1,ns,incx
          sy(i) = sx(i)
   70     continue
      return
      end
