c     algorithm 583, collected algorithms from acm.
c     algorithm appeared in acm-trans. math. software, vol.8, no. 2,
c     jun., 1982, p. 195.

c	Splus callable version			may 27, 1997 wle

	subroutine lsqr(m, n, damp,
     1	                leniw, lenrw, iw, rw,
     2	                u, v, w, x, se,
     3	                atol, btol, conlim, itnlim,
     4	                istop, anorm, acond, rnorm, arnorm, xnorm)

	implicit none

c	Parameters:
	integer	m
	integer	n
	real	damp
	integer	leniw, lenrw
	integer	iw(leniw)
	real	rw(lenrw)
	real	u(m), v(n), w(n)
	real	x(n)
	real	se(n)
	real	atol, btol
	real	conlim
	integer	itnlim
	integer	istop
	real	anorm
	real	acond
	real	rnorm
	real	arnorm
	real	xnorm

c     ------------------------------------------------------------------
c
c     lsqr  finds a solution  x  to the following problems...
c
c     1. unsymmetric equations --    solve  a*x = b
c
c     2. linear least squares  --    solve  a*x = b
c                                    in the least-squares sense
c
c     3. damped least squares  --    solve  (   a    )*x = ( b )
c                                           ( damp*i )     ( 0 )
c                                    in the least-squares sense
c
c     where  a  is a matrix with  m  rows and  n  columns, b  is an
c     m-vector, and  damp  is a scalar (all quantities real).
c     the matrix  a  is intended to be large and sparse.  it is accessed
c     by means of subroutine calls of the form
c
c                call aprod( mode,m,n,x,y,leniw,lenrw,iw,rw )
c
c     which must perform the following functions...
c
c                if mode = 1, compute  y = y + a*x.
c                if mode = 2, compute  x = x + a(transpose)*y.
c
c     the vectors x and y are input parameters in both cases.
c     if mode = 1, y should be altered without changing x.
c     if mode = 2, x should be altered without changing y.
c     the parameters leniw, lenrw, iw, rw may be used for workspace
c     as described below.
c
c     the rhs vector  b  is input via  u,  and subsequently overwritten.
c
c
c     note.  lsqr uses an iterative method to approximate the solution.
c     the number of iterations required to reach a certain accuracy
c     depends strongly on the scaling of the problem.  poor scaling of
c     the rows or columns of  a  should therefore be avoided where
c     possible.
c
c     for example, in problem 1 the solution is unaltered by
c     row-scaling.  if a row of  a  is very small or large compared to
c     the other rows of  a,  the corresponding row of  ( a  b )  should
c     be scaled up or down.
c
c     in problems 1 and 2, the solution  x  is easily recovered
c     following column-scaling.  in the absence of better information,
c     the nonzero columns of  a  should be scaled so that they all have
c     the same euclidean norm (e.g.  1.0).
c
c     in problem 3, there is no freedom to re-scale if  damp  is
c     nonzero.  however, the value of  damp  should be assigned only
c     after attention has been paid to the scaling of  a.
c
c     the parameter  damp  is intended to help regularize
c     ill-conditioned systems, by preventing the true solution from
c     being very large.  another aid to regularization is provided by
c     the parameter  acond,  which may be used to terminate iterations
c     before the computed solution becomes very large.
c
c
c     notation
c     --------
c
c     the following quantities are used in discussing the subroutine
c     parameters...
c
c     abar   =  (   a    ),          bbar  =  ( b )
c               ( damp*i )                    ( 0 )
c
c     r      =  b  -  a*x,           rbar  =  bbar  -  abar*x
c
c     rnorm  =  sqrt( norm(r)**2  +  damp**2 * norm(x)**2 )
c            =  norm( rbar )
c
c     relpr  =  the relative precision of floating-point arithmetic
c               on the machine being used.  for example, on the ibm 370,
c               relpr  is about 1.0e-6 and 1.0d-16 in single and double
c               precision respectively.
c
c     lsqr  minimizes the function  rnorm  with respect to  x.
c
c
c     parameters
c     ----------
c
c     m       input      the number of rows in  a.
c
c     n       input      the number of columns in  a.
c
c     damp    input      the damping parameter for problem 3 above.
c                        (damp  should be 0.0 for problems 1 and 2.)
c                        if the system  a*x = b  is incompatible, values
c                        of  damp  in the range 0 to sqrt(relpr)*norm(a)
c                        will probably have a negligible effect.
c                        larger values of  damp  will tend to decrease
c                        the norm of  x  and to reduce the number of
c                        iterations required by lsqr.
c
c                        the work per iteration and the storage needed
c                        by lsqr are the same for all values of  damp.
c
c     leniw   input      the length of the workspace array  iw.
c     lenrw   input      the length of the workspace array  rw.
c     iw      workspace  an integer array of length  leniw.
c     rw      workspace  a real array of length  lenrw.
c
c             note.  lsqr does not explicitly use the previous four
c             parameters, but passes them to subroutine aprod for
c             possible use as workspace.  if aprod does not need
c             iw  or  rw,  the values  leniw = 1  or  lenrw = 1  should
c             be used, and the actual parameters corresponding to
c             iw  or  rw  may be any convenient array of suitable type.
c
c     u(m)    input      the rhs vector  b.  beware that  u  is
c                        over-written by lsqr.
c
c     v(n)    workspace
c     w(n)    workspace
c
c     x(n)    output     returns the computed solution  x.
c
c     se(n)   output     returns standard error estimates for the
c                        components of  x.  for each  i,  se(i)  is set
c                        to the value  rnorm * sqrt( sigma(i,i) / t ),
c                        where  sigma(i,i)  is an estimate of the i-th
c                        diagonal of the inverse of abar(transpose)*abar
c                        and  t = 1      if  m .le. n,
c                             t = m - n  if  m .gt. n  and  damp = 0,
c                             t = m      if  damp .ne. 0.
c
c     atol    input      an estimate of the relative error in the data
c                        defining the matrix  a.  for example,
c                        if  a  is accurate to about 6 digits, set
c                        atol = 1.0e-6 .
c
c     btol    input      an estimate of the relative error in the data
c                        defining the rhs vector  b.  for example,
c                        if  b  is accurate to about 6 digits, set
c                        btol = 1.0e-6 .
c
c     conlim  input      an upper limit on  cond(abar),  the apparent
c                        condition number of the matrix  abar.
c                        iterations will be terminated if a computed
c                        estimate of  cond(abar)  exceeds  conlim.
c                        this is intended to prevent certain small or
c                        zero singular values of  a  or  abar  from
c                        coming into effect and causing unwanted growth
c                        in the computed solution.
c
c                        conlim  and  damp  may be used separately or
c                        together to regularize ill-conditioned systems.
c
c                        normally,  conlim  should be in the range
c                        1000  to  1/relpr.
c                        suggested value --
c                        conlim = 1/(100*relpr)  for compatible systems,
c                        conlim = 1/(10*sqrt(relpr))  for least squares.
c
c             note.  if the user is not concerned about the parameters
c             atol, btol  and  conlim,  any or all of them may be set
c             to zero.  the effect will be the same as the values
c             relpr, relpr  and  1/relpr  respectively.
c
c     itnlim  input      an upper limit on the number of iterations.
c                        suggested value --
c                        itnlim = n/2     for well conditioned systems,
c                        itnlim = 4*n     otherwise.
c
c     istop   output     an integer giving the reason for termination...
c
c                0       x = 0  is the exact solution.
c                        no iterations were performed.
c
c                1       the equations  a*x = b  are probably
c                        compatible.  norm(a*x - b)  is sufficiently
c                        small, given the values of  atol  and  btol.
c
c                2       the system  a*x = b  is probably not
c                        compatible.  a least-squares solution has
c                        been obtained which is sufficiently accurate,
c                        given the value of  atol.
c
c                3       an estimate of  cond(abar)  has exceeded
c                        conlim.   the system  a*x = b  appears to be
c                        ill-conditioned.  otherwise, there could be an
c                        an error in subroutine aprod.
c
c                4       the equations  a*x = b  are probably
c                        compatible.  norm(a*x - b)  is as small as
c                        seems reasonable on this machine.
c
c                5       the system  a*x = b  is probably not
c                        compatible.  a least-squares solution has
c                        been obtained which is as accurate as seems
c                        reasonable on this machine.
c
c                6       cond(abar)  seems to be so large that there is
c                        not much point in doing further iterations,
c                        given the precision of this machine.
c                        there could be an error in subroutine aprod.
c
c                7       the iteration limit  itnlim  was reached.
c
c     anorm   output     an estimate of the frobenius norm of  abar.
c                        this is the square-root of the sum of squares
c                        of the elements of  abar.
c                        if  damp  is small and if the columns of  a
c                        have all been scaled to have length  1.0,
c                        anorm  should increase to roughly  sqrt(n).
c                        a radically different value for  anorm  may
c                        indicate an error in subroutine aprod (there
c                        may be an inconsistency between modes 1 and 2).
c
c     acond   output     an estimate of  cond(abar),  the condition
c                        number of  abar.  a very high value of  acond
c                        may again indicate an error in aprod.
c
c     rnorm   output     an estimate of the final value of norm(rbar),
c                        the function being minimized (see notation
c                        above).  this will be small if  a*x = b  has
c                        a solution.
c
c     arnorm  output     an estimate of the final value of
c                        norm( abar(transpose)*rbar ), the norm of
c                        the residual for the usual normal equations.
c                        this should be small in all cases.  (arnorm
c                        will often be smaller than the true value
c                        computed from the output vector  x.)
c
c     xnorm   output     an estimate of the norm of the final
c                        solution vector  x.
c
c
c     subroutines and functions used
c     ------------------------------
c
c     user       aprod
c     lsqr       normlz
c     blas       scopy,snrm2,sscal  (see lawson et al. below)
c                (snrm2 is used only in normlz)
c     fortran    abs,mod,sqrt
c
c
c     precision
c     ---------
c
c     the number of iterations required by lsqr will usually decrease
c     if the computation is performed in higher precision.  to convert
c     lsqr  and  normlz  between single- and double-precision, change
c     the words
c                scopy, snrm2, sscal
c                abs, real, sqrt
c     to the appropriate blas and fortran equivalents.
c
c
c     references
c     ----------
c
c     paige, c.c. and saunders, m.a.  lsqr: an algorithm for sparse
c        linear equations and sparse least squares.
c        acm transactions on mathematical software 8, 1 (march 1982).
c
c     lawson, c.l., hanson, r.j., kincaid, d.r. and krogh, f.t.
c        basic linear algebra subprograms for fortran usage.
c        acm transactions on mathematical software 5, 3 (sept 1979),
c        308-323 and 324-325.
c
c
c     lsqr.      this version dated 22 february 1982.
c     ------------------------------------------------------------------

c	Functions and local variables
	real	abs
	real	alfa
	real	bbnorm
	real	beta
	real	bnorm
	real	cs
	real	cs1,cs2
	real	ctol
	real	dampsq
	real	ddnorm
	real	delta
	real	gambar
	real	gamma
	integer	i
	integer	itn
	integer	nconv
	integer	nout
	integer	nstop
	real	one
	real	phi
	real	phibar
	real	psi
	real	res1,res2
	real	rhbar1,rhbar2
	real	rho
	real	rhobar
	real	rhs
	real	rtol
	real	sn
	real	sn1,sn2
	real	sqrt
	real	t
	real	t1,t2,t3
	real	tau
	real	test1,test2,test3
	real	theta
	real	xxnorm
	real	z
	real	zbar
	real	zero

c     initialize.
      zero   = 0.0
      one    = 1.0
      ctol   = zero
      if (conlim .gt. zero) ctol = one/conlim
      dampsq = damp**2
      anorm  = zero
      acond  = zero
      bbnorm = zero
      ddnorm = zero
      res2   = zero
      xnorm  = zero
      xxnorm = zero
      cs2    = -one
      sn2    = zero
      z      = zero
      itn    = 0
      istop  = 0
      nstop  = 0

      do 10 i = 1, n
         v(i) = zero
         x(i) = zero
        se(i) = zero
   10 continue

c     set up the first vectors for the bidiagonalization.
c     these satisfy   beta*u = b,   alfa*v = a(transpose)*u.

      call normlz( m,u,beta )
      call aprod ( 2,m,n,v,u,leniw,lenrw,iw,rw )
      call normlz( n,v,alfa )
      call scopy ( n,v,1,w,1 )

      rhobar = alfa
      phibar = beta
      bnorm  = beta
      rnorm  = beta
      arnorm = alfa*beta
      if (arnorm .le. zero) go to 800
      if (nout   .le.  0  ) go to 100
      test1  = one
      test2  = alfa/beta

c     ------------------------------------------------------------------
c     main iteration loop.
c     ------------------------------------------------------------------
  100 itn = itn + 1

c     perform the next step of the bidiagonalization to obtain the
c     next  beta, u, alfa, v.  these satisfy the relations
c                beta*u  =  a*v  -  alfa*u,
c                alfa*v  =  a(transpose)*u  -  beta*v.

      call sscal ( m,(-alfa),u,1 )
      call aprod ( 1,m,n,v,u,leniw,lenrw,iw,rw )
      call normlz( m,u,beta )
      bbnorm = bbnorm + alfa**2 + beta**2 + dampsq
      call sscal ( n,(-beta),v,1 )
      call aprod ( 2,m,n,v,u,leniw,lenrw,iw,rw )
      call normlz( n,v,alfa )

c     use a plane rotation to eliminate the damping parameter.
c     this alters the diagonal (rhobar) of the lower-bidiagonal matrix.

      rhbar2 = rhobar**2 + dampsq
      rhbar1 = sqrt(rhbar2)
      cs1    = rhobar/rhbar1
      sn1    = damp/rhbar1
      psi    = sn1*phibar
      phibar = cs1*phibar

c     use a plane rotation to eliminate the subdiagonal element (beta)
c     of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.

      rho    = sqrt(rhbar2 + beta**2)
      cs     = rhbar1/rho
      sn     = beta/rho
      theta  =  sn*alfa
      rhobar = -cs*alfa
      phi    =  cs*phibar
      phibar =  sn*phibar
      tau    =  sn*phi

c     update  x, w  and the standard error estimates.

      t1 =    phi/rho
      t2 = -theta/rho
      t3 =    one/rho

      do 200 i = 1, n
         t     = w(i)
         x(i)  = t1*t + x(i)
         w(i)  = t2*t + v(i)
         t     =(t3*t)**2
         se(i) = t + se(i)
         ddnorm= t + ddnorm
  200 continue

c     use a plane rotation on the right to eliminate the
c     super-diagonal element (theta) of the upper-bidiagonal matrix.
c     then use the result to estimate  norm(x).

      delta  =  sn2*rho
      gambar = -cs2*rho
      rhs    = phi - delta*z
      zbar   = rhs/gambar
      xnorm  = sqrt(xxnorm + zbar**2)
      gamma  = sqrt(gambar**2 + theta**2)
      cs2    = gambar/gamma
      sn2    = theta/gamma
      z      = rhs/gamma
      xxnorm = xxnorm + z**2

c     test for convergence.
c     first, estimate the norm and condition of the matrix  abar,
c     and the norms of  rbar  and  abar(transpose)*rbar.

      anorm  = sqrt(bbnorm)
      acond  = anorm*sqrt(ddnorm)
      res1   = phibar**2
      res2   = res2 + psi**2
      rnorm  = sqrt(res1 + res2)
      arnorm = alfa*abs(tau)

c     now use these norms to estimate certain other quantities,
c     some of which will be small near a solution.

      test1  = rnorm/bnorm
      test2  = arnorm/(anorm*rnorm)
      test3  = one/acond
      t1     = test1/(one + anorm*xnorm/bnorm)
      rtol   = btol +  atol*anorm*xnorm/bnorm

c     the following tests guard against extremely small values of
c     atol, btol  or  ctol.  (the user may have set any or all of
c     the parameters  atol, btol, conlim  to zero.)
c     the effect is equivalent to the normal tests using
c     atol = relpr,  btol = relpr,  conlim = 1/relpr.

      t3 = one + test3
      t2 = one + test2
      t1 = one + t1
      if (itn .ge. itnlim) istop = 7
      if (t3  .le. one   ) istop = 6
      if (t2  .le. one   ) istop = 5
      if (t1  .le. one   ) istop = 4

c     allow for tolerances set by the user.

      if (test3 .le. ctol) istop = 3
      if (test2 .le. atol) istop = 2
      if (test1 .le. rtol) istop = 1
c     ==================================================================

c     stop if appropriate.
c     the convergence criteria are required to be met on  nconv
c     consecutive iterations, where  nconv  is set below.
c     suggested value --   nconv = 1, 2  or  3.

      if (istop .eq. 0) nstop = 0
      if (istop .eq. 0) go to 100
      nconv = 1
      nstop = nstop + 1
      if (nstop .lt. nconv  .and.  itn .lt. itnlim) istop = 0
      if (istop .eq. 0) go to 100
c     ------------------------------------------------------------------
c     end of iteration loop.
c     ------------------------------------------------------------------


c     finish off the standard error estimates.

      t = one
      if (m .gt. n) t = m - n
      if (dampsq .gt. zero) t = m
      t = rnorm/sqrt(t)

      do 700 i = 1, n
         se(i) = t*sqrt(se(i))
  700 continue

c     end of execution

c	set itnlim = itn upon return		6/18/97 wle
	itnlim = itn
  800 return
c     ------------------------------------------------------------------
c     end of lsqr
      end
