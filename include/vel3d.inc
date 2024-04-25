      integer	maxnx,maxny,maxnz,maxnz2
      integer	ixkms,iykms,izkms
c
c   maxnx: maximum number of nodes in x direction
      parameter(maxnx= 25)
c   maxny: maximum number of nodes in y direction
      parameter(maxny= 25)
c   maxnz: maximum number of nodes in z direction
      parameter(maxnz= 20)
c                                                                                                             
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   velocity arrays
c
c   /vmod3d/
c    vel(nx,ny,2*nz)     ! 1st half is P velocity, 2nd half is S velocity
      parameter(maxnz2=maxnz*2)
c
      real	bld
      integer	nx,ny,nz
      real	xn(maxnx), yn(maxny), zn(maxnz)
      integer	iuses
      real	vel(maxnx,maxny,maxnz2)
      real	vpvs(maxnx, maxny, maxnz)
      real	clat_3d,clon_3d,crot_3d
      common/vmod3d/bld,nx,ny,nz,xn,yn,zn,iuses,
     *  vel,vpvs
c
c   /locate/
c   ixkms: size of map in x direction.  The distance between
c            furthest nodes.  This will be in units of km if bld=1.0
c            or in units of 0.1km if bld=0.1
      parameter(ixkms=10000)
c   iykms=size of map in y dir.
      parameter(iykms=10000)
c   izkms=size of map in z dir.
      parameter(izkms=10000)
c
      real	xl, yl, zl
      integer 	ixloc(ixkms), iyloc(iykms), izloc(izkms)
      common/locate/ xl,yl,zl,ixloc,iyloc,izloc
c
c                                                                                                             
      integer	ndip, iskip
      real	scale1, scale2
      real	xfac, tlim
      integer	nitpb
      real	rp(3,130,1)
      integer	nrp(1)
      common/pseudo/ndip,iskip,scale1,scale2,xfac,tlim,nitpb,
     &              rp,nrp
