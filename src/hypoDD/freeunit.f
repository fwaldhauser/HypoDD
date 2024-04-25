c Find a free fortran i/o unit-number

	subroutine FREEUNIT(iunit)

	implicit none

c	Parameters:
	integer	iunit

c	Local variables:
	logical	lopen

      do iunit=10,999
         if(iunit.eq.999) stop'FREEUNIT>>> no free unit found!'
         inquire(unit=iunit,opened=lopen)
         if(.not.lopen) RETURN
      enddo
      RETURN
      end ! of subr. freeunit
