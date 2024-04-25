c Multiply a float vector by a scalar
c     Replace single precision sx by single precision sa*sx.
c     for i = 0 to n-1, replace sx(1+i*incx) with  sa * sx(1+i*incx)

	subroutine  sscal(n, sa, sx, incx)

	implicit none

c	Parameters:
	integer	n
	integer	incx
	real	sa, sx(1+(n-1)*incx)

c	Local variables:
	integer	i
	integer	m
	integer	mp1
	integer	ns

      if(n.le.0)return
      if(incx.eq.1)goto 20

c        code for increments not equal to 1.

      ns = n*incx
          do 10 i = 1,ns,incx
          sx(i) = sa*sx(i)
   10     continue
      return

c        code for increments equal to 1.


c        clean-up loop so remaining vector length is a multiple of 5.

   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sx(i) = sa*sx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        sx(i) = sa*sx(i)
        sx(i + 1) = sa*sx(i + 1)
        sx(i + 2) = sa*sx(i + 2)
        sx(i + 3) = sa*sx(i + 3)
        sx(i + 4) = sa*sx(i + 4)
   50 continue
      return
      end
