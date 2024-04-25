	subroutine indexxi(n, iarrin, indx)

	implicit none

	integer	n
	integer	iarrin(n)
	integer	indx(n)

	integer	i
	integer	indxt
	integer	ir
	integer	j
	integer	l
	integer	q

      if (n.lt.1) then
         return
      else if (n.eq.1) then
         indx(1) = 1
         return
      endif

      do 11 j=1,n
         indx(j) = j
11    continue
      l = n/2+1
      ir = n
10    continue
         if (l.gt.1) then
            l = l-1
            indxt = indx(l)
            q = iarrin(indxt)
         else
            indxt = indx(ir)
            q = iarrin(indxt)
            indx(ir) = indx(1)
            ir = ir-1
            if (ir.eq.1) then
               indx(1) = indxt
               return
            endif
         endif
         i = l
         j = l+l
20       if (j.le.ir) then
            if (j.lt.ir) then
               if (iarrin(indx(j)).lt.iarrin(indx(j+1))) j = j+1
            endif
            if (q.lt.iarrin(indx(j))) then
               indx(i) = indx(j)
               i = j
               j = j+j
            else
               j = ir+1
            endif
         go to 20
         endif
         indx(i) = indxt
      go to 10
      end
