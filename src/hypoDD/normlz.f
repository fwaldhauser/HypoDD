c Compute Euclidean norm of float vector, and normalize vector if norm != 0.
c (Required by subroutine lsqr.)

	subroutine normlz(n, x, beta)

	implicit none

c	Parameters:
	integer	n
	real	x(n)
	real	beta	! Euclidean norm of x	(output)

c	Local variables:
	real	one
	real	snrm2	! BLAS function
	real	zero

      zero = 0.0
      one  = 1.0
      beta = snrm2(n, x, 1)
      if (beta .gt. zero) call sscal(n, (one/beta), x, 1)
      return

c     end of normlz
      end
