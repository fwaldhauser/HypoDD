	subroutine exist(fn)

	implicit none

c	Parameters:
	character	fn*110

c	Local variables
	logical	ex

      inquire(FILE=fn,exist=ex)
      if (.not.ex) then
          write(*,'("FILE DOES NOT EXIST / CHECK IDATA,IPHASE: ",a)')fn
          stop
      endif
      end
