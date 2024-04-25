c Matrix mutliply: c(m,m) = a(m,m)*b(m,m)

	subroutine matmult1(maxm, m, a, b, c)

	implicit none

c	Parameters:
	integer	maxm
	integer	m
	real	a(maxm,maxm)	! (1..maxm, 1..maxm)
	real	b(maxm,maxm)	! (1..maxm, 1..maxm)
	real	c(maxm,maxm)	! (1..maxm, 1..maxm)

c	Local variables:
	integer	j, l, k		! Dummy loop indices
	real	sum

      do j=1,m
         do k=1,m
            sum=0.0
            do l=1,m
               sum=sum+a(j,l)*b(l,k)
            end do
            c(j,k)=sum
         end do
      end do
      return
      end
