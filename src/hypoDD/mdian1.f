c Median value in float vector
c Sorts the vector as a side-effect.

	subroutine mdian1(x, n, xmed)

	implicit none

c	Parameters:
	integer	n
	real	x(n)
	real	xmed

c	Local variables:
	integer	n2

      CALL SORT(N,X)
      N2=N/2
      IF(2*N2.EQ.N)THEN
        XMED=0.5*(X(N2)+X(N2+1))
      ELSE
        XMED=X(N2+1)
      ENDIF
      RETURN
      END
