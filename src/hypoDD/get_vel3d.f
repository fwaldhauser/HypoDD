c reads in 3D velocity model defined on grid nodes
c linear interpolation between neighbouring nodes
c for more details pleass look up simul2000 by
c Cliff Thurbee & Donna Eberhart-Phillips
c
	subroutine get_vel3d(fn_vel3d,iphase,lat_3d,lon_3d,rot_3d,
     &			     indip,iiskip,iscale1,iscale2,
     &			     ixfac, itlim, initpb,
     &                       log)

	implicit none
        
	include'vel3d.inc'

c	Parameters
	character	fn_vel3d*110	!3D velocity model
	integer		iphase		! phases used for location
	integer		log		! log file unit
	real		lat_3d		! origin of 3D model
	real		lon_3d		! origin of 3D model
	real		rot_3d		! rotation angle of 3D model
	integer		indip,iiskip,initpb	! parameters for 3d raytracing
	real		iscale1,iscale2,ixfac, itlim
	real		vv(maxnx)	

c	local varaibles
	integer		i,j,k,ks,k2,kv
	integer		ltype
	integer		iunit
	integer		ierror		!error flag for vel model
        integer 	ixf
	integer		iyf
        integer		izf
        integer 	ixm
 	integer		iym
	integer		izm
        integer 	ixl
	integer		iyl
	integer		izl
        character*1 	vtype(2)

	call freeunit(iunit)
	open(iunit, file=fn_vel3d, status='unknown')

c
      ierror=0
      vtype(1)='P'
      vtype(2)='S'
c
      if(iphase.eq.1) then
	iuses=1
      else
	iuses=2
      endif

c  set up global origin for 3D model
c      write(*,*)lat_3d, lon_3d, rot_3d
ccc      call setorg_3d(lat_3d, lon_3d, rot_3d)
c
c  for this version the gridpoints can be unevenly spaced
c  the origin of the coordinate system is at (x,y,z)=(0,0,0)
c  which will not in general correspond to the point
c  xn(1),yn(1),zn(1).
c  xn,yn,zn should be factors of bld (ie: a.0 for bld=1.0 or a.b for bld=0.1)
c
c input the number of gridpoints in x, y and z directions
c  and bld factor (1.0 or 0.1 km) used to set up velocity interpolation grid
      read(iunit,*) bld,nx,ny,nz
      if((bld.ne.1.0).and.(bld.ne.0.1)) then
        write(log,1625) bld
 1625   format(/, '******** STOP *********, bld must be 1.0 or 0.1,
     2   not ',f6.2)
	stop
      endif
cfw120603...
      if(nx.gt.MAXNX.or.ny.gt.MAXNY.or.nz.gt.MAXNZ) 
     *stop'Increase 3D model arrays MAXNX/Y/Z in vel3d.inc'
c
c  input the x grid, y grid, and z grid
        read(iunit,*) (xn(i),i=1,nx)
        read(iunit,*) (yn(i),i=1,ny)
        read(iunit,*) (zn(i),i=1,nz)
c
 3003 format(3i3)    
      write(log,3005) bld,nx,ny,nz
 3005 format(//,' velocity grid size:',/,
     * 'bld =',f4.1,5x,' nx =',i3,5x,'ny =',i3,5x,'nz =',i3)
c
cfh give all these numbers the same format
      write(log,3006) (xn(i),i=1,nx)
cfh 3006 format(/,' xgrid',/,3x,12f7.1,8f6.1)
 3006 format(/,' xgrid',/,3x,20f7.1)
      write(log,3007) (yn(i),i=1,ny)
cfh 3007 format(/,' ygrid',/,3x,12f7.1,8f6.1)
 3007 format(/,' ygrid',/,3x,20f7.1)
      write(log,3008) (zn(i),i=1,nz)
cfh 3008 format(/,' zgrid',/,3x,8f6.1,12f7.1/)
 3008 format(/,' zgrid',/,3x,20f7.1/)
c
c  read in which nodes to have fixed velocity
c  end with blank line
   50 read(iunit,3003) ixf,iyf,izf
      if(ixf.le.0) goto 60
      goto 50
   60 continue
c
c  start cht 1998
c  lines moved followed by new code
c  compute total number of gridpoints (nodes)
c  MAYBE
c     nodes=nx*ny*nz
c     nxy=nx*ny
c     nx2=nx-2             ! number non-edge nodes in row
c     nxy2=nx2*(ny-2)      ! number non-edge nodes in layer
c     nz2=nz-2
c     nodes2=nz2*nxy2
c  peripheral nodes
c     nx1=nx-1
c     ny1=ny-1
c     nz1=nz-1
c
c  read in which nodes have "linked" velocity, "master" node first
c  followed by "linked" nodes - end each group with blank line
c  end with another blank line.  if no linked nodes, just include
c  a blank line
c
 52   read(iunit,3003) ixm,iym,izm
      if(ixm.le.0) goto 62
c
c  link type - constant (1) or linear (2)?
      read(iunit,*) ltype
c
c
 54   read(iunit,3003) ixl,iyl,izl
      if(ixl.le.0) goto 52
      goto 54
 62   continue
c
c  end cht 1998
c
c  now read in the velocity values
   65 write(log,3101)
c     do 38 kv=1,iuses
         kv=1
         do 37 k=1,nz
            k2=k + (kv-1)*nz
            write(log,3015) k,vtype(kv),zn(k)
            do 36 j=1,ny
               read(iunit,*) (vel(i,j,k2),i=1,nx)
ccfw fix positive west problem:
c		do i=1,nx
c               		vv(i)=vel(i,j,k2)
c		enddo
c		do i=1,nx
c               		vel(i,j,k2)= vv(nx+1-i)
c		enddo

               write(log,3013) (vel(i,j,k2),i=1,nx)
   36       continue
   37    continue

         
         
c  38 continue
c CHANGE FOR VP/VS INVERSION
      if(iuses.eq.2) then
        do 100 k=1,nz
           write(log,3016) k,zn(k)
           do 99 j=1,ny
cfw             read(iunit,*) (vpvs(i,j,k),i=1,nx)
             read(iunit,*,end=999) (vpvs(i,j,k),i=1,nx)
             write(log,3013) (vpvs(i,j,k),i=1,nx)
   99     continue
  100  continue
c  compute Vs from Vp and Vp/Vs or compute 1/tstar
        kv=2
        do 120 k=1,nz
              ks=k+nz
              write(log,3015) k,vtype(kv),zn(k)  
              do 115 j=1,ny
                 do 110 i=1,nx
                    vel(i,j,ks)=vel(i,j,k)/vpvs(i,j,k)
  110           continue
                 write(log,3013) (vel(i,j,ks),i=1,nx)
  115        continue
  120     continue
      endif
c
 3013 format(20f6.2)
 3014 format(20f7.1)
 3015 format(/,' layer',i3,5x,a1,' velocity',10x,'z =',f7.1)
 3016 format(/,' layer',i3,5x,'Vp/Vs',10x,'z =',f7.1)
 3011 format(20f5.2)
 3101 format(//,' velocity values on three-dimensional grid')
c
c  Set up an array which is used to point to node indices, for any x,y,z
      call bldmap
c
c  Set up global raytracing parameters
      ndip = indip
      iskip = iiskip
      scale1 = iscale1
      scale2 = iscale2
      xfac = ixfac
      tlim = itlim
      nitpb = initpb
c     write(*,*) "GET: ",ndip,iskip,scale1,scale2,xfac,tlim,nitpb
      return
c
999   stop'>>>No S-velocities provided in 3D model file. Set IPHA=1.'
c
c***** end of subroutine get_vel3d *****
      end

c
      subroutine bldmap
c     called from get_vel3d
c
c  common block variables:
      include 'vel3d.inc'
c
      xl=bld-xn(1)
      ixmax=(xn(nx)+xl)/bld
      yl=bld-yn(1)
      iymax=(yn(ny)+yl)/bld
      zl=bld-zn(1)
      izmax=(zn(nz)+zl)/bld
c
c  Check for array size overflow
      if(ixmax.gt.ixkms.or.iymax.gt.iykms.or.izmax.gt.izkms)goto 330
      ix=1
      do 10 i=1,ixmax
c
         ix1=ix+1
c
         xnow=float(i)*bld-xl
         if (xnow.ge.xn(ix1)) ix=ix1
c
         ixloc(i)=ix
   10 continue
c  Fill remainder of array with zeroes.
      do 12 i=ixmax,ixkms
         ixloc(i)=0
   12 continue
c
c
      iy=1
      do 15 i=1,iymax
c
         iy1=iy+1
c
         ynow=float(i)*bld-yl
         if (ynow.ge.yn(iy1)) iy=iy1
c
         iyloc(i)=iy
   15 continue
c
c  Fill rest of array with zeroes.
      do 17 i=iymax,iykms
         iyloc(i)=0
 17   continue
c
      iz=1
      do 20 i=1,izmax
c
         iz1=iz+1
c
         znow=float(i)*bld-zl
         if (znow.ge.zn(iz1)) iz=iz1
c
         izloc(i)=iz
   20 continue
c
c  Fill remainder of array with zeroes.
      do 22 i=izmax,izkms
         izloc(i)=0
  22  continue
      return
 330   continue
      write(log,331)ixkms,iykms,izkms
 331  format(' ***** error in array size in common/locate/',/,
     *' maximum map dimensions (km)=',/,' x=',i5,' y=',i5,' z=',i5)
      write(log,332)ixmax,iymax,izmax
  332 format(' Actual map size (km): ',/,' x=',i5,' y=',i5,' z=',i5)
      stop
c***** end of subroutine bldmap *****
      end

