c Singular-value decomposition of rectangular float matrix

      subroutine svd(m1, n1, m, n, a, u, v, q, index)
C$$$$$ CALLS NO OTHER ROUTINES
C  SINGULAR VALUE DECOMPOSITION)  FOR ALGOL PROGRAM SEE WILKINSON+REINSCH
C  HANDBOOK FOR AUTOMATIC COMPUTATION VOL 2 - LINEAR ALGEBRA, PP140-144
C  TRANSLATED FROM ALGOL BY R.L.PARKER
C  THE MATRIX A(M,N) IS DECOMPOSED.  SINGULAR VALUES IN Q, PRE-MATRIX IN U,
C  POST-MATRIX IN V.   INDEX MAY BE 1,2,3 OR 4.  IF 1, FIND U,V. IF 2, FIND
C  ONLY U. IF 3, FIND ONLY V. IF 4, FIND NEITHER. IN ALL CASES, THE ARRAY  U
C  MUST BE SUPPLIED AS IT IS USED AS WORKING SPACE FOR THE ROUTINE.
C  PROGRAM ALTERED BY PAUL SILVER 4/15 TO HANDLE UNPACKED ARRAYS.
C  M1,N1 ARE DIMENSIONS IN MAIN ROUTINE.M,N ARE ACTUAL DIMENSIONS TO
C  BE USED IN THE SUBROUTINE.

	implicit none

c	Parameters:
	integer	m1, n1
	integer	m, n
	integer	index
	doubleprecision	a(m1,n)	! (1..m, 1..n)
	doubleprecision	u(m1,n)	! (1..m, 1..n)
	doubleprecision	v(n1,n)	! (1..n, 1..n)
	doubleprecision	q(n)	! (1..n)

c	Local variables
	doubleprecision	c
	doubleprecision	e(1000)
	doubleprecision	eps
	doubleprecision	f
	doubleprecision	g
	doubleprecision	h
	integer	i, j, k, l
	integer	iback
	integer	kback
	integer	lback
	integer	lplus
	integer	l1
	doubleprecision	s
	doubleprecision	tol
	doubleprecision	x
	doubleprecision	y
	doubleprecision	z

      EPS=1.0E-10
      TOL=1.0E-35
      DO 1100 I=1,M
      DO 1100 J=1,N
 1100 U(I,J)=A(I,J)
C  HOUSEHOLDER REDUCTION TO BI-DIAGONAL FORM
      G=0.0
      X=0.0
      DO 2900 I=1,N
      E(I)=G
      S=0.0
      L=I+1
      DO 2100 J=I,M
 2100 S=U(J,I)**2 + S
      IF (S .LT. TOL) GO TO 2500
      F=U(I,I)
      G=-SIGN(SQRT(S),F)
      H=F*G - S
      U(I,I)=F - G
      IF (L.GT.N) GO TO 2501
      DO 2400 J=L,N
      S=0.0
      DO 2200 K=I,M
 2200 S=U(K,I)*U(K,J) + S
      F=S/H
      DO 2300 K=I,M
 2300 U(K,J)=U(K,J) + F*U(K,I)
 2400 CONTINUE
      GO TO 2501
 2500 G=0.0
C
 2501 CONTINUE
      Q(I)=G
      S=0.0
      IF (L.GT.N) GO TO 2601
      DO 2600 J=L,N
 2600 S=U(I,J)**2 + S
 2601 IF (S.LT.TOL) GO TO 2800
      F=U(I,I+1)
      G=-SIGN(SQRT(S),F)
      H=F*G - S
      U(I,I+1)=F - G
      IF (L.GT.N) GO TO 2651
      DO 2650 J=L,N
 2650 E(J)=U(I,J)/H
 2651 CONTINUE
      IF (L.GT.M) GO TO 2850
      DO 2700 J=L,M
      S=0.0
      IF (L.GT.N) GO TO 2700
      DO 2670 K=L,N
 2670 S=U(J,K)*U(I,K) + S
      DO 2690 K=L,N
 2690 U(J,K)=U(J,K) + S*E(K)
 2700 CONTINUE
      GO TO 2850
 2800 G=0.0
 2850 Y=ABS(Q(I)) + ABS(E(I))
      IF (Y .GT. X) X=Y
 2900 CONTINUE
C
C  ACCUMULATION OF RIGHT-HAND TRANSFORMS (V)
C
      GO TO (3000,3701,3000,3701       ),INDEX
 3000 CONTINUE
      DO 3700 IBACK=1,N
      I=N+1-IBACK
      IF (G .EQ. 0.0) GO TO 3500
      H=U(I,I+1)*G
      IF (L.GT.N) GO TO 3500
      DO 3100 J=L,N
 3100 V(J,I)=U(I,J)/H
      DO 3400 J=L,N
      S=0.0
      DO 3200 K=L,N
 3200 S=U(I,K)*V(K,J) + S
      DO 3300 K=L,N
 3300 V(K,J)=V(K,J) + S*V(K,I)
 3400 CONTINUE
 3500 CONTINUE
      IF (L.GT.N) GO TO 3601
      DO 3600 J=L,N
      V(J,I)=0.0
 3600 V(I,J)=0.0
 3601 V(I,I)=1.0
      G=E(I)
      L=I
 3700 CONTINUE
 3701 CONTINUE
C
C  ACCUMULATION OF LEFT-HAND TRANSFORMS
      GO TO (4000,4000,4701,4701       ),INDEX
 4000 CONTINUE
      DO 4700 IBACK=1,N
      I=N+1-IBACK
      L=I+1
      G=Q(I)
      IF (L.GT.N) GO TO 4101
      DO 4100 J=L,N
 4100 U(I,J)=0.0
 4101 IF (G.EQ. 0.0) GO TO  4500
      H=U(I,I)*G
      IF (L.GT.N) GO TO 4401
      DO 4400 J=L,N
      S=0.0
      DO 4200 K=L,M
 4200 S=U(K,I)*U(K,J) + S
      F=S/H
      DO 4300 K=I,M
 4300 U(K,J)=U(K,J) + F*U(K,I)
 4400 CONTINUE
 4401 CONTINUE
      DO 4550 J=I,M
 4550 U(J,I)=U(J,I)/G
      GO TO 4700
 4500 CONTINUE
      DO 4600 J=I,M
 4600 U(J,I)=0.0
 4700 U(I,I)=U(I,I) + 1.0
C
C  DIAGONALIZATION OF BI-DIAGONAL FORM
 4701 EPS=EPS*X
      DO 9000 KBACK=1,N
      K=N+1-KBACK
C  TEST F-SPLITTING
 5000 CONTINUE
      DO 5100 LBACK=1,K
      L=K+1-LBACK
      IF (ABS(E(L)).LE. EPS) GO TO 6500
      IF (ABS(Q(L-1)) .LE. EPS) GO TO 6000
 5100 CONTINUE
C  CANCELLATION OF E(L), IF L.GT. 1
 6000 C=0.0
      S=1.0
      L1=L - 1
      DO 6200 I=L,K
      F=S*E(I)
                   E(I)=C*E(I)
      IF (ABS(F) .LE. EPS) GO TO 6500
      G=Q(I)
      Q(I)=SQRT(F*F + G*G)
      H=Q(I)
      C=G/H
      S=-F/H
      GO TO (6050,6050,6200,6200       ),INDEX
 6050 CONTINUE
      DO 6100 J=1,M
      Y=U(J,L1)
      Z=U(J,I)
      U(J,L1)=Y*C + Z*S
      U(J,I)=-Y*S + Z*C
 6100 CONTINUE
 6200 CONTINUE
C  TEST F-CONVERGENCE
 6500 Z=Q(K)
      IF (L .EQ. K) GO TO  8000
C  SHIFT FROM BOTTOM 2 X 2 MINOR
      X=Q(L)
      Y=Q(K-1)
      G=E(K-1)
      H=E(K)
      F=((Y-Z)*(Y+Z) + (G-H)*(G+H))/(2.0*H*Y)
      G=SQRT(F*F + 1.0)
      F=((X-Z)*(X+Z) + H*(Y/(F + SIGN(G,F))-H))/X
C  NEXT Q-R TRANSFORMATION
      C=1.0
      S=1.0
      LPLUS=L + 1
      DO 7500 I=LPLUS,K
      G=E(I)
      Y=Q(I)
      H=S*G
      G=C*G
      Z=SQRT(F*F + H*H)
      E(I-1)=Z
      C=F/Z
      S=H/Z
      F=X*C + G*S
      G=-X*S + G*C
      H=Y*S
      Y=Y*C
      GO TO (7100,7201,7100,7201       ),INDEX
 7100 DO 7200 J=1,N
      X=V(J,I-1)
      Z=V(J,I)
      V(J,I-1)=X*C + Z*S
      V(J,I)=-X*S + Z*C
 7200 CONTINUE
 7201 Z=SQRT(F*F + H*H)
      Q(I-1)=Z
      C=F/Z
      S=H/Z
      F=C*G + S*Y
      X=-S*G + C*Y
      GO TO (7300,7300,7500,7500       ),INDEX
 7300 DO 7400 J=1,M
      Y=U(J,I-1)
      Z=U(J,I)
      U(J,I-1)=Y*C + Z*S
      U(J,I)=-Y*S + Z*C
 7400 CONTINUE
 7500 CONTINUE
      E(L)=0.0
      E(K)=F
      Q(K)=X
      GO TO  5000
C  CONVERGENCE
 8000 IF (Z .GE. 0.0) GO TO 9000
C  Q IS MADE NON-NEGATIVE
      Q(K)=-Z
      GO TO (8100,9000,8100,9000       ),INDEX
 8100 DO 8200 J=1,N
 8200 V(J,K)=-V(J,K)
 9000 CONTINUE
      RETURN
      END
