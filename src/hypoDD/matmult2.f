c Matrix mutliply: c(n) = a(m,n)*b(n)

	subroutine matmult2(maxm, m, n, a, b, c)

	implicit none

c	Parameters:
	integer	maxm
	integer	m
	integer	n
	doubleprecision	a(maxm,n)	! (1..m, 1..n)
	real	b(n)		! (1..n)
	real	c(n)		! (1..n)

c	Local variables:
	integer	j, l
	real	sum

      do j=1,m
            sum=0.0
            do l=1,n
               sum=sum+a(j,l)*b(l)
            end do
            c(j)=sum
      end do
      return
      end
