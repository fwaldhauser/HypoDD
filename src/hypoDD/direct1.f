c Compute travel time, etc., for direct (upward-departing) ray

	subroutine direct1(nl,v,vsq,thk,jl,tkj,delta,depth,tdir,u,x)

	implicit none

        include "hypoDD.inc"

c	Parameters:
	integer	nl	! Number of layers		(input)
	real	v(MAXLAY)	! Layer wave speeds		(input)
	real	vsq(MAXLAY)	! Squares of wave speeds	(input)
	real	thk(MAXLAY)	! Layer thicknesses		(input)
	integer	jl	! Event layer			(input)
	real	tkj	! Event depth within layer jl	(input)
	real	delta	! Epicentral distance		(input)
	real	depth	! Event depth			(input)
	real	tdir	! Direct-ray travel time	(output)
	real	u	! Sine of take-off angle	(output)
	real	x	! Horizontal travel distance in event layer (output)

c       For the direct seismic ray from an event to a receiver in
c  a layered velocity structure, direct predicts the travel time, the
c  sine of the takeoff angle, and the horizontal distance of travel in
c  the event layer.  The receiver must be located at the top of layer
c  1 and the event must be located below layer 1.  Low velocity
c  layers are permitted.
c       To find the takeoff angle of the ray, a numerical approach
c  is required.  The basic scheme adopted here is the method of false
c  position.  (see acton, 1970, 'numerical methods that work,' for
c  example.)  First, the characteristics of the fastest layer
c  between the event and the surface are determined.  These permit
c  placing definite lower and upper bounds, ua and ub, on the
c  sine of the takeoff angle.  In turn, lower and upper bounds, xa
c  and xb, on the horizontal travel distance in the event layer are
c  determined.  The total horizontal travel distance for a ray with
c  with horizontal travel distance x in the event layer is denoted
c  by del, and the zero of del - delta is obtained by using xa and
c  xb as initial guesses for x in the method of false position
c  from x and tkj, the depth of the event below the top of the event
c  layer, the sine of the takeoff angle, u , is calculated.
c       From u and x, tdir is found by summing the travel time in
c  each layer.  finally, a slight correction to tdir is made, based
c  on the misfit between the final del and delta.

c	Local variables:
	real		del		! Computed distance
	real		dela, delb	! Distances corresponding to xa, xb
	doubleprecision	hypot		! Hypoteneuse function
	integer		j1
	integer		kount
	integer		l
	integer		lmax
	real		r
	real		tklmax
	real		usq
	real		ua, uasq
	real		ub, ubsq
	real		ubdiv
	real		vlmax
	real		xa, xb		! Bounds on x
	real		xtest

	if (jl .eq. 1) then
c	   Focus in surface layer
	   r = hypot(depth, delta)
	   tdir = r/v(1)
	   u = delta/r
	   x = delta
	   return
	endif

c     Find the fastest layer, lmax, above and including jl
      lmax = jl
      tklmax = tkj
      vlmax = v(jl)
      j1 = jl-1
      do 23184 l=1,j1
         if (.not.(v(l).gt.vlmax)) goto 23186
            lmax = l
            tklmax = thk(l)
            vlmax = v(l)
23186    continue
23184 continue

C CHANGE BY E.KISSLING MARCH 1984
      IF(tklmax.le.0.05) tklmax = 0.05

c     Find initial bounds on the sine of the takeoff angle
      ua = (v(jl)/vlmax)*delta/sqrt(delta**2+depth**2)
      ub = (v(jl)/vlmax)*delta/sqrt(delta**2+tklmax**2)

c     Calculate horizontal travel distances
      uasq = ua**2
      ubsq = ub**2
C CHANGE BY E.KISSLING MARCH 1984
      if (ubsq.ge.1.) ubsq = 0.99999
      if (uasq.ge.1.) uasq = 0.99999
      xa = tkj*ua/sqrt(1.0-uasq)
      if (.not.(lmax.eq.jl)) goto 23188
         xb = delta
         goto 23189
23188 continue
      xb = tkj*ub/sqrt(1.0-ubsq)
23189 continue
      dela = xa
      delb = xb
      do 23190 l=1,j1
         dela = dela+thk(l)*ua/sqrt(vsq(jl)/vsq(l)-uasq)
         ubdiv = sqrt(vsq(jl)/vsq(l)-ubsq)
         if (ubdiv.GT.1.e-20) GOTO 1002
c	    No write statements for Splus!
            ubdiv = 1.e-20
 1002    continue
         delb = delb+thk(l)*ub/sqrt(vsq(jl)/vsq(l)-ubsq)
23190 continue

c     Loop to find the zero of del-delta by the method of false position
      do 23192 kount=1,25
         if (.not.((delb-dela).lt.0.02)) goto 23194
            x = 0.5*(xa+xb)
            u = x/sqrt(x**2+tkj**2)
            usq = u**2
            goto 23193	! break
23194    continue
         x = xa+(delta-dela)*(xb-xa)/(delb-dela)
         u = x/sqrt(x**2+tkj**2)
         usq = u**2
         del = x
         do 23196 l=1,j1
            del = del+thk(l)*u/sqrt(vsq(jl)/vsq(l)-usq)
23196    continue
         xtest = del-delta
         if (abs(xtest).lt.0.02) goto 23193	! break
         if (.not.(xtest.lt.0.0)) goto 23200
            xa = x
            dela = del
            goto 23201
23200    continue
            xb = x
            delb = del
23201    continue
23192 continue
23193 continue

c     Calculate direct-ray travel time
      tdir = sqrt(x**2+tkj**2)/v(jl)
      do 23202 l=1,j1
         tdir = tdir+thk(l)*v(jl)/(vsq(l)*sqrt(vsq(jl)/vsq(l)-usq))
23202 continue
      tdir = tdir-(u/v(jl))*(del-delta)
      return
c  ***** end of subroutine direct1 *****
      end
