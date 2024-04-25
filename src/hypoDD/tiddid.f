c Compute intercept times and critical distances for "refracted" rays

	subroutine tiddid (jl, nl, v, vsq, thk, tid, did)

	implicit none

	include "hypoDD.inc"

c	Parameters:
	integer	jl
	integer	nl
	real	v(MAXLAY)	! (1..nl)
	real	vsq(MAXLAY)	! (1..nl)
	real	thk(MAXLAY)	! (1..nl)
	real	tid(MAXLAY)	! (1..20)
	real	did(MAXLAY)	! (1..20)

c       Determines the travel time intercept and critical
c  distance for a seismic ray in a layered earth model
c  originating at the top of layer jl, refracting in
c  layer m, and terminating at the top of layer 1.
c
c  input:       jl - event layer
c               nl - number of layers
c             v(l) - velocity of layer l
c           vsq(l) - velocity squared
c           thk(l) - thickness of layer l
c  output:
c           tid(m) - travel time intercept for
c                      refraction in layer m
c           did(m) - critical distance

c	Local variables:
	real	did1, did2
	real	dimm
	integer	j1
	integer	l
	integer	m
	integer	m1
	real	sqt
	real	tid1, tid2
	real	tim

      j1=jl+1
      do 23174m=j1,nl
      tid(m)=0.
      did(m)=0.
      tid1=0.
      tid2=0.
      did1=0.
      did2=0.
      m1=m-1
      do 23176l=1,m1
      if(.not.(vsq(m).le.vsq(l)))goto 23178

c   if m is a low velocity layer, set tid and did to
c   very large values

      tid(m)=100000.
      did(m)=100000.
      goto 23179
23178 continue
      sqt=sqrt(vsq(m)-vsq(l))
      tim=thk(l)*sqt/(v(l)*v(m))
      dimm=thk(l)*v(l)/sqt
      if(.not.(l.lt.jl))goto 23180

c   sum for layers above event layer

      tid1=tid1+tim
      did1=did1+dimm
      goto 23181
23180 continue

c   sum for layers below and including the event layer

      tid2=tid2+tim
      did2=did2+dimm
23181 continue
23179 continue
23176 continue
      if(.not.(tid(m).ne.100000.))goto 23182

c   calculate tid and did if m is not a low velocity layer

      tid(m)=tid1+2*tid2
      did(m)=did1+2*did2
23182 continue
23174 continue
      return
c  ***** end of subroutine tiddid *****
      end
