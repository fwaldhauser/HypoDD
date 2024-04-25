c Pseudo-random number generator
c Numerical Recipes, p. 274

	subroutine ran(jlo, jhi, j)

	implicit none

c	Parameters:
	real	jlo, jhi	! Limit values (changed)
	real	j

c	Local variables:
	integer	im, ia, ic
	real	jran

	parameter (im= 714025, ia= 4096, ic= 150889)

      jhi=jhi-1.0
      jran= mod(jran*ia+ic,real(im))         ! generator
      j= jlo+((jhi-jlo+1)*jran)/im
      return
      end
