c Length of character string, excluding trailing blanks

	integer function trimlen(t)

	implicit none

c	Parameter:
	character t*(*)

      do 1 trimlen=LEN(t),1,-1
    1    if(t(trimlen:trimlen).ne.' ')RETURN
      trimlen=1
      end ! of integer function trimlen
