c Sort an int array

	subroutine sorti(n, ia)

	implicit none

c	Parameters:
	integer	n
	integer ia(n)	! (1..n)

c	Local variables:
	integer	i, j, l
	integer	iia
	integer	ir

      if (n.le.1) then
         return
      endif

      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          IIA=IA(L)
        ELSE
          IIA=IA(IR)
          IA(IR)=IA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            IA(1)=IIA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(IA(J).LT.IA(J+1))J=J+1
          ENDIF
          IF(IIA.LT.IA(J))THEN
            IA(I)=IA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        IA(I)=IIA
      GO TO 10
      END
