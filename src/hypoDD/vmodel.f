c Extract needed information from the layered velocity model.

	subroutine vmodel(nl, v, top, depth, vsq, thk, jl, tkj)

	implicit none

	include "hypoDD.inc"

c	Parameters:
	integer	nl
	real	v(MAXLAY)
	real	vsq(MAXLAY)
	real	top(MAXLAY)
	real	depth
	real	thk(MAXLAY)
	integer	jl
	real	tkj

c  input:     nl - number of layers
c           v(l) - velocity of layer l
c	     top - depth to top of layer l
c          depth - depth of event

c         vsq(l) = v(l) ** 2
c         thk(l) - thickness of layer l
c             jl - event layer
c            tkj - depth of event from top of event layer

c	Local variables:
	integer	i

c	compute square of layer velocity
	do 10 i=1,nl
   10	vsq(i)=v(i)*v(i)

c	determine layer thickness and
c	find layer containing event,

	jl=nl

	do 20 i=1,nl

c	Important note:  if (depth.lt.top(i)) will
c	lead to incorrect results for traveltime
	if (depth.le.top(i)) then
	jl=i-1
	goto 25
	endif
   20	continue
   25	continue

	do 30 i=1,nl-1
   30	thk(i)=top(i+1)-top(i)

c	compute depth from top of layer to source

	tkj=depth-top(jl)

	return
c *****	end of subroutine vmodel *****
	end
