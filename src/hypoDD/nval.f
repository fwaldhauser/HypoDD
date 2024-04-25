	subroutine nval(line, nv) 
c get number of values, separated by spaces, in a given line
	implicit none

c	Parameters:
c120217	character	line*120
	character	line*220
	integer nv, j, k, trimlen

        j= 0
        do k=2,trimlen(line)
           if(line(k:k).ne.' '.and.line(k-1:k-1).eq.' ') j=j+1
        enddo
        if(line(1:1).ne.' ') j=j+1
        nv= j
        end

