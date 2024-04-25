c Determine the fastest traveltime between a source
c at depth=depth(km) and a receiver at distance=delta(km).

	subroutine ttime(delta, depth, nl, v, top, t, ain)

	implicit none

	include "hypoDD.inc"

c	Parameters:
	real	delta
	real	depth
	integer	nl
	real	v(MAXLAY)
	real	top(MAXLAY)
	real	t
	real	ain

c	Local variables:
	integer	jl
	integer	kk
	real	tdir
	real	thk(MAXLAY)
	real	tkj
	real	tref
	real	u
	real	vsq(MAXLAY)
	real	x
	real	xovmax

c	compile and link for S
c	f77 -c ttime.f
c	ld -r -dn ttime.o
c	mv a.out ttime.o

c	subroutine direct1 is used to compute the direct ray
c	traveltime and sine of takeoff angle.

c	subroutine refract is used to compute the fastest
c	refracted ray traveltime.  It calls subroutine tiddid.

c	subroutine vmodel extract needed information from the
c	layered velocity model.

c	input:
c	delta	epicentral distance in km
c	depth	focal depth of source in km
c	nl	number of layers in velocity model
c	v	velocity in each layer
c	top	depth to top of layer

c	output:
c	t	minimum traveltime
c	ain	angle of emergence at source


c	call vmodel to set-up model and locate source in it

	call vmodel(nl,v,top,depth,vsq,thk,jl,tkj)

c  output:
c      vsq(l) - v(l) ** 2
c      thk(l) - thickness of layer l
c          jl - event layer
c         tkj - depth of event in event layer

c	call refract to find fastest refracted arrival

	call refract(nl,v,vsq,thk,jl,tkj,delta,
     &			kk,tref,xovmax)

c  output:   kk - refracting layer for fastest refracted ray
c          tref - travel time of fastest refracted ray
c        xovmax - an upper bound on delta for which the direct ray
c                 can be the first arrival


c	if delta <= xovmax, them
c	call direct1 to find the direct ray traveltime
c	otherwise tref is the minimum traveltime

c	assume for now refracted path is faster

	t=tref

c	compute the takeoff angle
	if (kk.gt.0) then

	u=v(jl)/v(kk)
	ain=asin(u)*57.2958
	endif

	if (delta.le.xovmax) then

	call direct1(nl,v,vsq,thk,jl,tkj,delta,depth,tdir,u,x)

c  output:  tdir - direct ray travel time
c              u - sine of the takeoff angle
c              x - horizontal travel distance in the event layer
c

c	compare the traveltimes

	if (tref.gt.tdir) then

c	direct time is the minimum traveltime

	t=tdir
	ain=180-asin(u)*57.2958

	endif
	endif

	return
c *****	end of subroutine ttime *****
	end
