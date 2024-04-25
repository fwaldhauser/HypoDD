c Matrix mutliply: c(m,n) = a(m) * b(m,n)
c overwrites matrix b

	subroutine matmult3(maxm, maxn, m, n, a, b)

	implicit none

c	Parameters:
	integer	maxm
	integer	maxn
	integer	m
	integer	n
	real	a(maxm)
	doubleprecision	b(maxm,maxn)

c	Local variables:
	integer i, j		! Dummy loop indices

      do i=1,m
            do j=1,n
                  b(i,j) = a(i)*b(i,j)
            enddo
      enddo
      return
      end
