      subroutine path(isp,xe,ye,ze,xr,yr,zr,ttime,az,toa,vv)
c  this routine determines the minimum-time ray path
c  in two steps:  first, an approximate path is
c  determined using approximate ray tracing
c  then the approximate path is used as a starting point in
c  shooting ray tracing routine to determine the path in 3-d
c  ***  note - current version does not contain full 3-d ray
c  ***  tracing - routines are under development
c
c  declaration statements:
c
c  common block variables:
      implicit none
      include 'vel3d.inc'

      integer	isp
      real	xe, ye, ze
      real	xr, yr, zr
      real	ttime
      real	az
      real	toa
      real	vv
      

      real    fstime
      integer jflag
      integer ncr
      integer ndp
      integer jpb
      real	x1,y1,z1
      real	x2,y2,z2
      
c     write(*,*) "PATH: ",ndip,iskip,scale1,scale2,xfac,tlim,nitpb

      jflag = 0
      call rayweb(isp,xe,ye,ze,xr,yr,zr,fstime,jflag,
     &            ncr, ndp)

      ttime = fstime
      call minima(isp,fstime,nitpb,jpb)

      x1 = rp(1,1,1)
      x2 = rp(1,2,1)
      y1 = rp(2,1,1)
      y2 = rp(2,2,1)
      z1 = rp(3,1,1)
      z2 = rp(3,2,1)
      call aztoa(x1,x2,y1,y2,z1,z2,xr,yr,az,toa)
c     write(*,*) ttime, fstime, az, toa

      call vel3(isp,xe,ye,ze,vv)

      ttime = fstime

      return
      end
c
c ---------------------------------------------------------
c
      subroutine rayweb(isp,xe,ye,ze,xr,yr,zr,
     *            fstime,jflag,ncrold,ndpold)
c  approximate ray tracing package art2
c    with fast ray tracing code
c     by Cliff Thurber (from his simul3l version)
c
c  common block variables:
      implicit none
      include 'vel3d.inc'

      integer isp
      real xe,ye,ze
      real xr,yr,zr
      real fstime
      integer jflag
      integer ncrold, ndpold

      real delx,dely,delz
      real sep,delsep
      real pthsep(130),strpth(390),fstpth(390)
      real dipvec(3,9),disvec(390,9)
      real trpath(390,9)
      real trtime(9)
      real tmin,tt
      real trpth1(390)
      real sn1
      real xstep,ystep,zstep

      integer nd,ns,ns1,ii,ncr,npt
      integer i,ic,nc,n1,n2,n3,nn,np,ndpfst
      integer ndip1,ndip2,nz1,ncr0,ncr1
      integer ndp,iz0,npt2,npt3

c
c  compute source-receiver separation                         
      delx=xr-xe                                              
      dely=yr-ye                                              
      delz=zr-ze                                              
      sep=sqrt(delx*delx+dely*dely+delz*delz)                 
c  determine integer parameters for set of curves to be constructed
      call setup(sep,scale1,scale2,nd,ns,npt,ncr)
c     write(*,*) "NCR: ", ncr, " NPT: " ,npt,"ND: ",nd,"NS: ",ns
c                                                             
c  set up pthsep array for straight-line travel time calculation
      sn1=1.0/ns                                              
      delsep=sep*sn1                                          
      do 20 i=1,ns                                            
         pthsep(i)=delsep                                     
   20 continue                                                
c                                                             
c  determine points along straight-line path                  
      xstep=delx*sn1                                          
      ystep=dely*sn1                                          
      zstep=delz*sn1                                          
c                                                             
      ic=0                                                    
      ns1=ns+1                                                
      do 25 ii=1,ns1                                          
c                                                             
         i=ii-1                                               
c                                                             
         ic=ic+1                                              
         strpth(ic)=xe+xstep*i                                
         fstpth(ic)=strpth(ic)                                
         ic=ic+1                                              
         strpth(ic)=ye+ystep*i                                
         fstpth(ic)=strpth(ic)                                
         ic=ic+1                                              
         strpth(ic)=ze+zstep*i                                
         fstpth(ic)=strpth(ic)                                
   25 continue                                                
c                                                             
c  compute travel time along straight-line path               
      call ttime_3d(isp,ns,npt,strpth,pthsep,fstime)             
c
      if (ncr.eq.1) go to 65                                  
c                                                             
c  compute the dip vectors of length scale2                   
      call cmpdpv(xe,ye,ze,xr,yr,zr,scale2,ndip,dipvec)       
c                                                             
c  compute the basic set of displacement vectors              
      call cmpdsv(ndip,iskip,ns,dipvec,disvec)                
c                                                             
c  set first and last points of all trial paths to source and receiver      
      n1=3*npt-2                                              
      n2=n1+1                                                 
      n3=n1+2                                                 
      ndip1=1+iskip                                           
      ndip2=ndip-iskip                                        
c                                                             
c  fast ray tracing code                                      
c                                                             
      nz1=nz-1                                                
      ncr0=1                                                  
      ncr1=ncr-1                                              
      if (jflag.eq.0) go to 28                                
      ncr0=ncrold-1                                           
      if (ncr0.lt.1) ncr0=1                                   
      ncr1=ncrold+1                                           
      if (ncr1.gt.ncr) ncr1=ncr                               
c  ndip was changed (in main) for ihomo iterations., make ndpold=vertical plane.
c     if(jfl.eq.1) ndpold=(ndip+1)/2     ! 21-feb-86, this is done in strt, so unnecessary here
      ndip1=ndpold-1                                          
      if (ndip1.lt.1+iskip) ndip1=1+iskip                     
      ndip2=ndpold+1                                          
      if (ndip2.gt.ndip-iskip) ndip2=ndip-iskip               
   28 continue                                                
c  set "old" values to straight line                          
      ncrold=0                                                
      ndpold=(ndip+1)/2                                       
c                                                             
c                                                             
      do 30 ndp=ndip1,ndip2                                   
         trpath(1,ndp)=xe                                     
         trpath(2,ndp)=ye                                     
         trpath(3,ndp)=ze                                     
         trpath(n1,ndp)=xr                                    
         trpath(n2,ndp)=yr                                    
         trpath(n3,ndp)=zr                                    
   30 continue                                                
      trpth1(1)=xe                                            
      trpth1(2)=ye                                            
      trpth1(3)=ze                                            
      trpth1(n1)=xr                                           
      trpth1(n2)=yr                                           
      trpth1(n3)=zr                                           

c     write(*,*) ncr0,ncr1,jflag,ncrold,ndpold
c                                                             
c  loop over the curve sets                                   
      do 40 nc=ncr0,ncr1                                      
         iz0=0                                                
c                                                             
c  loop over different dips for one set of curves             
         do 42 ndp=ndip1,ndip2                                
c                                                             
            npt2=npt-2
c  loop to determine points along one path                    
            do 44 np=1,npt2                                   
               n1=3*np+1                                      
               n3=n1+2                                        
               do 43 nn=n1,n3                                 
                  trpath(nn,ndp)=nc*disvec(nn,ndp)+strpth(nn) 
                  trpth1(nn)=trpath(nn,ndp)                   
   43          continue                                       
   44       continue                                          
c                                                             
c  set up pthsep array for travel time calculations           
            if (ndp.eq.ndip1) call cmpsep(trpth1,pthsep,ns)   
c  compute travel time along one path                         
            call ttime_3d(isp,ns,npt,trpth1,pthsep,tt)           
            trtime(ndp)=tt                                    
   42    continue                                             
c                                                             
c  sort through trtime to find fastest path from current set  
         tmin=1.0e15                                          
         do 50 ndp=ndip1,ndip2                                
            if (trtime(ndp).gt.tmin) go to 50                 
            tmin=trtime(ndp)                                  
            ndpfst=ndp                                        
   50    continue                                             
c                                                             
c  compare fastest trtime to current value of fstime          
c  replace fstime and fstpth if needed                        
         if (tmin.ge.fstime) go to 40                         
         fstime=tmin                                          
c  reset "old" values                                         
         ncrold=nc                                            
         ndpold=ndpfst                                        
c                                                             
         npt3=3*npt                                           
         do 52 np=1,npt3                                      
            fstpth(np)=trpath(np,ndpfst)                      
   52    continue                                             
c                                                             
   40 continue                                                
c                                                             
c  put fstpth into rp array                                   
   65 continue                                                
      do 60 np=1,npt                                          
         n3=3*np                                              
         n1=n3-2                                              
         n2=n1+1                                              
         rp(1,np,1)=fstpth(n1)                               
         rp(2,np,1)=fstpth(n2)                               
         rp(3,np,1)=fstpth(n3)                               
   60 continue                                                
      nrp(1)=npt                                             
c                                                             
      return                                                  
      end                                                     

c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine setup(sep,scale1,scale2,nd,ns,npt,ncr)
c
c  parameters
      real sep,scale1,scale2
c
      integer nd,ns,npt,ncr
c
c  determine the number of path divisions - use scale1
      nd=1+nint(3.32193*log10(sep/scale1))
      if (nd.gt.7) nd=7
c  number of segments along path
      ns=2**nd
c  number of points on path
      npt=ns+1
c
c  determine the number of curves - use scale2
      ncr=1+0.5*sep/scale2
c                                                             
      if (sep.gt.scale1) return                               
      nd=0                                                    
      ns=1                                                    
      npt=2                                                   
      ncr=1                                                   
c                                                             
      return                                                  
      end                                                     

c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine ttime_3d(isp,ns,npt,pathr,pthsep,tt)
c
c  travel time along path via trapezoidal rule integration
c
c  parameters
      real pathr(390),pthsep(130),tt
c
      integer isp,ns,npt
c  local variables
      real vpt(130),x,y,z,v
c
      integer ip,np
c
c  loop over points along path to determine velocities
      ip=0
      do 10 np=1,npt
         ip=ip+1
         x=pathr(ip)                                          
         ip=ip+1                                              
         y=pathr(ip)                                          
         ip=ip+1                                              
         z=pathr(ip)                                          
         call vel3(isp,x,y,z,v)                               
   10 vpt(np)=v                                               
c                                                             
      tt=0.0                                                  
c  sum up travel time - use trapezoidal rule                  
      do 20 np=1,ns                                           
         np1=np+1                                             
c  Check for value outside defined area                       
         if((vpt(np).le.0.0).or.(vpt(np1).le.0.0)) goto 99    
         tt=tt+pthsep(np)/(vpt(np)+vpt(np1))                  
   20 continue                                                
      tt=2.0*tt                                               
c                                                             
      return                                                  
   99 write(16,100) np,vpt(np),np1,vpt(np1)                   
  100 format(' **** ERROR IN TTIME ****, velocity le 0.0',    
     2 /,'    np   vpt(np)   np1   vpt(np1)',/,1x,2(i5,f10.2))
      stop                                                    
      end                                                     

c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine cmpdpv(xe,ye,ze,xr,yr,zr,scale2,ndip,dipvec)
c
c  parameters
      real xe,ye,ze,xr,yr,zr,scale2,dipvec(3,9)
c
      integer ndip
c  local variables
      real dx,dy,dz,xh1,yh1,zh1,xh2,yh2,zh2,size,xv,yv,zv,
     *     rescal,x451,y451,z451,x452,y452,z452
c
      integer nv
c
      dx=xr-xe
      dy=yr-ye
      dz=zr-ze
c
c  near-vertical vector
      xv=-dx*dz
      yv=-dy*dz
      zv=dx*dx+dy*dy
c  rescale vector to length scale2
      size=sqrt(xv*xv+yv*yv+zv*zv)
      rescal=scale2/size
c
      xv=xv*rescal
      yv=yv*rescal
      zv=zv*rescal
c
c  store this vector
      nv=(ndip+1)/2
      dipvec(1,nv)=xv
      dipvec(2,nv)=yv
      dipvec(3,nv)=zv                                         
c                                                             
      if (ndip.eq.1) return                                   
c                                                             
c  horizontal vectors                                         
      xh1=dy                                                  
      yh1=-dx                                                 
      zh1=0.0                                                 
      xh2=-dy                                                 
      yh2=dx                                                  
      zh2=0.0                                                 
c  rescale the vectors to length scale2                       
      size=sqrt(xh1*xh1+yh1*yh1)                              
      rescal=scale2/size                                      
c                                                             
      xh1=xh1*rescal                                          
      yh1=yh1*rescal                                          
      xh2=xh2*rescal                                          
      yh2=yh2*rescal                                          
c                                                             
c  store these two vectors                                    
      dipvec(1,1)=xh1                                         
      dipvec(2,1)=yh1                                         
      dipvec(3,1)=zh1                                         
c                                                             
      dipvec(1,ndip)=xh2                                      
      dipvec(2,ndip)=yh2                                      
      dipvec(3,ndip)=zh2                                      
c                                                             
      if (ndip.eq.3) return                                   
c                                                             
c  determine two 45 degree dip vectors                        
      rescal=0.7071068                                        
c                                                             
      n1=(1+nv)/2                                             
      n2=(nv+ndip)/2                                          
c                                                             
      x451=(xh1+xv)*rescal                                    
      y451=(yh1+yv)*rescal                                    
      z451=(zh1+zv)*rescal                                    
c                                                             
      x452=(xh2+xv)*rescal                                    
      y452=(yh2+yv)*rescal                                    
      z452=(zh2+zv)*rescal                                    
c                                                             
      dipvec(1,n1)=x451                                       
      dipvec(2,n1)=y451                                       
      dipvec(3,n1)=z451                                       
c                                                             
      dipvec(1,n2)=x452                                       
      dipvec(2,n2)=y452                                       
      dipvec(3,n2)=z452                                       
c                                                             
      if (ndip.eq.5) return                                   
c                                                             
c  determine four 22.5 degree dip vectors                     
      rescal=0.5411961                                        
c                                                             
      dipvec(1,2)=(xh1+x451)*rescal                           
      dipvec(2,2)=(yh1+y451)*rescal                           
      dipvec(3,2)=(zh1+z451)*rescal                           
c                                                             
      dipvec(1,4)=(x451+xv)*rescal                            
      dipvec(2,4)=(y451+yv)*rescal                            
      dipvec(3,4)=(z451+zv)*rescal                            
c                                                             
      dipvec(1,6)=(xv+x452)*rescal                            
      dipvec(2,6)=(yv+y452)*rescal                            
      dipvec(3,6)=(zv+z452)*rescal                            
c                                                             
      dipvec(1,8)=(x452+xh2)*rescal                           
      dipvec(2,8)=(y452+yh2)*rescal                           
      dipvec(3,8)=(z452+zh2)*rescal                           
c                                                             
      return                                                  
      end                                                     

c                                                             
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -    
c                                                             
      subroutine cmpdsv(ndip,iskip,ns,dipvec,disvec)          
c                                                             
c  parameters                                                 
      real dipvec(3,9),disvec(390,9)                          
c                                                             
      integer ndip,iskip,ns                                   
c  local variables                                            
      real darc(129)                                          
c                                                             
      integer inc,narc,ndp,np,n,nd                            
c  coefficients of standard arc                               
c                                                             
       data darc/0.,0.0342346,0.0676860,0.1003707,0.1323047,0.1635029,      
     *.1939796,.2237485,.2528226,.2812141,.3089350,.3359963,  
     *.3624091,.3881833,.4133289,.4378553,.4617713,.4850857,  
     *.5078065,.5299417,.5514988,.5724850,.5929070,.6127717,  
     *.6320853,.6508538,.6690831,.6867787,.7039459,.7205898,  
     *.7367154,.7523272,.7674298,.7820274,.7961241,.8097238,  
     *.8228301,.8354468,.8475771,.8592244,.8703916,.8810817,  
     *.8912975,.9010416,.9103164,.9191245,.9274679,.9353487,  
     *.9427691,.9497307,.9562353,.9622845,.9678797,.9730224,  
     *.9777138,.9819550,.9857470,.9890908,.9919872,.9944367,.9964401,       
     *.9979979,.9991102,.9997776,1.0000000,.9997776,.9991102,.9979979,      
     *.9964401,.9944367,.9919872,.9890908,.9857470,.9819550,.9777138,       
     *.9730224,.9678797,.9622845,.9562353,.9497307,.9427691,.9353487,       
     *.9274679,.9191245,.9103164,.9010416,.8912975,.8810817,.8703916,       
     *.8592244,.8475771,.8354468,.8228301,.8097238,.7961241,.7820274,       
     *.7674298,.7523272,.7367154,.7205898,.7039459,.6867787,.6690831,       
     *.6508538,.6320853,.6127717,.5929070,.5724850,.5514988,.5299417,       
     *.5078065,.4850857,.4617713,.4378553,.4133289,.3881833,.3624091,       
     *.3359963,.3089350,.2812141,.2528226,.2237485,.1939796,.1635029,       
     *.1323047,.1003707,.0676860,.0342346,0.0/                

      inc=128/ns                                              
        ndip1=1+iskip                                         
      ndip2=ndip-iskip                                        
c                                                             
c  loop over dips                                             
      do 30 ndp=ndip1,ndip2                                   
         narc=1                                               
         nd=3                                                 
c  loop over points on the path (skip first and last)         
         do 20 np=2,ns                                        
            narc=narc+inc                                     
c                                                             
c  loop over x,y,z                                            
            do 10 n=1,3                                       
               nd=nd+1                                        
               disvec(nd,ndp)=darc(narc)*dipvec(n,ndp)        
c                                                             
   10       continue                                          
   20    continue                                             
   30 continue                                                
c                                                             
      return                                                  
      end                                                     

c                                                             
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -    
c                                                             
      subroutine cmpsep(path,pthsep,ns)                       
c                                                             
c  parameters                                                 
      real path(390),pthsep(130)                              
c                                                             
      integer ns                                              
c  local variables                                            
      integer nx,ny,nz,nx1,ny1,nz1                            
c                                                             
      nx=-2                                                   
c  loop over pairs of points in one set of stored vectors     
      do 10 n=1,ns                                            
         nx=nx+3                                              
         ny=nx+1                                              
         nz=nx+2                                              
         nx1=nx+3                                             
         ny1=nx+4                                             
         nz1=nx+5                                             
c                                                             
         pthsep(n)=sqrt((path(nx1)-path(nx))**2               
     *              +(path(ny1)-path(ny))**2                  
     *              +(path(nz1)-path(nz))**2)                 
c                                                             
   10 continue                                                
c                                                             
      return                                                  
      end                                                     


c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine vel3(isp,x,y,z,v)
c  This routine is Cliff Thurber's
c  common block variables:
      common/weight/ wv(8),ip,jp,kp,kpg
      include 'vel3d.inc'
c
c  use Prothero's intmap here
      call intmap(x,y,z,ip,jp,kp)
c
      ip1=ip+1
      jp1=jp+1
      kp1=kp+1
c
c       write(16,100)x,xl,ip,y,yl,jp,z,zl,kp
c100    format(3(2f7.3,i3))
      xf=(x-xn(ip))/(xn(ip1)-xn(ip))
      yf=(y-yn(jp))/(yn(jp1)-yn(jp))
      zf=(z-zn(kp))/(zn(kp1)-zn(kp))
      xf1=1.0-xf
      yf1=1.0-yf
      zf1=1.0-zf
c
      wv(1)=xf1*yf1*zf1
      wv(2)=xf*yf1*zf1
      wv(3)=xf1*yf*zf1
      wv(4)=xf*yf*zf1
      wv(5)=xf1*yf1*zf
      wv(6)=xf*yf1*zf
      wv(7)=xf1*yf*zf
      wv(8)=xf*yf*zf
c  calculate velocity
c  S-velocity is stored after P-velocity
c  (or V*Q if iuseq=1)
      kpg=kp
      if(isp.eq.1) kp=kp+nz
      kp1=kp+1
      v=wv(1)*vel(ip,jp,kp)+wv(2)*vel(ip1,jp,kp)
     2 +wv(3)*vel(ip,jp1,kp)+wv(4)*vel(ip1,jp1,kp)
     * +wv(5)*vel(ip,jp,kp1)+wv(6)*vel(ip1,jp,kp1)
     * +wv(7)*vel(ip,jp1,kp1)+wv(8)*vel(ip1,jp1,kp1)
      return
      end

c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine intmap(x,y,z,ip,jp,kp)
c  Modified by W. Prothero so a single call can get the indices
c  common block variables:
      include 'vel3d.inc'
c
      ip=int((x+xl)/bld)
      ip=ixloc(ip)
      jp=int((yl+y)/bld)
      jp=iyloc(jp)
      kp=int((z+zl)/bld)
      kp=izloc(kp)
c  If an array element=0, the position is off the map.
      return
      end

c
c-------------------------------------------------------
c
      subroutine minima(isp,fstime,npbmax,jpb)
c
c*****this routine finds the minimum path using pseudo-bending
c
c  common block variables:
      common/pathm/x(130),y(130),z(130),v(130),vq(130),tra,qtra,n,nn
      common/temp/xtemp(130),ytemp(130),ztemp(130),rtemp(130),ttemp(130)
      common/pat/xt(130),yt(130),zt(130),rt(130),rm,rs,rtm,rts,nt,j
      common/upb/ upbtot,nupb
      include 'vel3d.inc'
c
c     write(*,*) "MINI: ",ndip,iskip,scale1,scale2,xfac,tlim,nitpb
      tra=fstime
      n=nrp(1)
c
      do 20 i=1,n
         x(i)=rp(1,i,1)
         y(i)=rp(2,i,1)
         z(i)=rp(3,i,1)
   20 continue
c
      do 21 i=1,n
         xp=x(i)
         yp=y(i)
         zp=z(i)
         xtemp(i)=xp
         ytemp(i)=yp
         ztemp(i)=zp
         call vel3(isp,xp,yp,zp,vp)
         v(i)=vp
         call vel3(1,xp,yp,zp,vpp)
         vq(i)=vpp
   21 continue
      call travel
c                                                             
      nn=n-1                                                  
      nupb=nupb+1                                             
c     write(6,6000) nupb                                      
c     write(*,*) "MAXVALS: ", xfac, npbmax
 6000 format(' nubp=',i8)                                     
      do 100 j=1,npbmax                                       
         ta=tra                                               
         call bend(isp,xfac)                                  
         call travel                                          
         deltat=ta-tra                                        
c                                                             
         if (deltat.lt.0.0) go to 102                         
         do 22 i=1,n                                          
            x(i)=xtemp(i)                                     
            y(i)=ytemp(i)                                     
            z(i)=ztemp(i)                                     
   22    continue                                             
         if(deltat.le.tlimm) go to 102                         
100   continue                                                
c                                                             
c                                                             
  102 continue                                                
                                                              
      if(j.gt.npbmax) j=npbmax                                
      upbtot=upbtot + float(j)                                
      jpb=j                                                   
  105 do 300 i=1,n                                            
         rp(1,i,1)=x(i)                                      
         rp(2,i,1)=y(i)                                      
         rp(3,i,1)=z(i)                                      
300   continue                                                
      fstime=tra                                            
c                                                             
      return                                                  
      end                                                     

c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine travel
c
      common/temp/xtemp(130),ytemp(130),ztemp(130),rtemp(130),ttemp(130)
c
      common/pathm/x(130),y(130),z(130),v(130),vq(130),tra,qtra,n,nn
c
      tra=0
      do 60 i=2,n
         i1=i-1
         xd=xtemp(i)-xtemp(i1)
         yd=ytemp(i)-ytemp(i1)
         zd=ztemp(i)-ztemp(i1)
         ds=sqrt(xd*xd+yd*yd+zd*zd)
         tv=ds*(1.0/v(i)+1.0/v(i1))
         tra=tra+tv
  60  continue
      tra=0.5*tra
c
      return
      end

c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine bend(isp,xfac)
c*****this routine perturbs the initial path in the direction
c      of the normal to the ray path tangent at each point
c      by the optimal distance r
c
      common/pathm/x(130),y(130),z(130),v(130),vq(130),tra,qtra,n,nn
      common/temp/xtemp(130),ytemp(130),ztemp(130),rtemp(130),ttemp(130)
c
c ***
      xtemp(1)=x(1)
      ytemp(1)=y(1)
      ztemp(1)=z(1)
c ***
      do 200 k=2,nn
c
         kk=k-1
         kkk=k+1
c
c*****compute the normal direction of maximum gradient of velocity
c
         dx=x(kkk)-xtemp(kk)
         dy=y(kkk)-ytemp(kk)
         dz=z(kkk)-ztemp(kk)
         dn=dx*dx+dy*dy+dz*dz
         ddn=sqrt(dn)
         rdx=dx/ddn
         rdy=dy/ddn
         rdz=dz/ddn
c
         xk=0.5*dx+xtemp(kk)
         yk=0.5*dy+ytemp(kk)
         zk=0.5*dz+ztemp(kk)
c ***
         call vel3(isp,xk,yk,zk,vk)                           
         call veld(isp,xk,yk,zk,vx,vy,vz)                     
c                                                             
c ***                                                         
         vrd=vx*rdx+vy*rdy+vz*rdz                             
         rvx=vx-vrd*rdx                                       
         rvy=vy-vrd*rdy                                       
         rvz=vz-vrd*rdz                                       
c                                                             
         rvs=sqrt(rvx*rvx+rvy*rvy+rvz*rvz)                    
         if(rvs.eq.0.0) then                                  
           xtemp(k) = xk                                      
           ytemp(k) = yk                                      
           ztemp(k) = zk                                      
           v(k) = vk                                          
           call vel3(1,xk,yk,zk,vkk)                          
           vq(k)=vkk                                          
         else                                                 
           rvx=rvx/rvs                                        
           rvy=rvy/rvs                                        
           rvz=rvz/rvs                                        
                                                            
c*****compute the optimal distance r                          
           rcur=vk/rvs                                        
           rtemp(k)=rcur-sqrt(rcur*rcur-0.25*dn)              
c                                                             
c*****compute the new points and distance of perturbations    
c                                                             
           xxk=xk+rvx*rtemp(k)                                
           yyk=yk+rvy*rtemp(k)                                
           zzk=zk+rvz*rtemp(k)                                
c                                                             
c  convergence enhancement                                    
           xxk=xfac*(xxk-x(k))+x(k)                           
           yyk=xfac*(yyk-y(k))+y(k)                           
           zzk=xfac*(zzk-z(k))+z(k)                           
c                                                             
           ttemp(k)=sqrt((x(k)-xxk)**2+(y(k)-yyk)**2+(z(k)-zzk)**2)
           xtemp(k)=xxk                                       
           ytemp(k)=yyk                                       
           ztemp(k)=zzk                                       
           call vel3(isp,xxk,yyk,zzk,vk)                      
           v(k)=vk                                            
           call vel3(1,xxk,yyk,zzk,vkk)                       
           vq(k)=vkk                                          
         endif                                                
c ***                                                         
200   continue                                                
c                                                             
      return                                                  
      end                                                     

c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine veld(isp,xx,yy,zz,vx,vy,vz)
c
c*****this routine computes the derivatives of velocity
c     in x, y, and z directions
c
c  common block variables:
      include 'vel3d.inc'
c
c  use Prothero's intmap here
      call intmap(xx,yy,zz,ip,jp,kp)
c
      ip1=ip+1
      jp1=jp+1
      kp1=kp+1
c
      xd=xn(ip1)-xn(ip)
      yd=yn(jp1)-yn(jp)
      zd=zn(kp1)-zn(kp)
c
      xf=(xx-xn(ip))/xd
      yf=(yy-yn(jp))/yd
      zf=(zz-zn(kp))/zd
c
      xf1=1.0-xf
      yf1=1.0-yf
      zf1=1.0-zf
c
c  S-velocity is stored in the 2nd half of the velocity array
      if(isp.eq.1) kp=kp+nz
      kp1=kp+1
c
c*****calculate derivatives of velocity
c
      vx=(yf1*zf1*(vel(ip1,jp,kp)-vel(ip,jp,kp))              
     *+yf*zf1*(vel(ip1,jp1,kp)-vel(ip,jp1,kp))                
     *+yf1*zf*(vel(ip1,jp,kp1)-vel(ip,jp,kp1))                
     *+yf*zf*(vel(ip1,jp1,kp1)-vel(ip,jp1,kp1)))/xd           
c                                                             
      vy=(xf1*zf1*(vel(ip,jp1,kp)-vel(ip,jp,kp))              
     *+xf*zf1*(vel(ip1,jp1,kp)-vel(ip1,jp,kp))                
     *+xf1*zf*(vel(ip,jp1,kp1)-vel(ip,jp,kp1))                
     *+xf*zf*(vel(ip1,jp1,kp1)-vel(ip1,jp,kp1)))/yd           
c                                                             
      vz=(xf1*yf1*(vel(ip,jp,kp1)-vel(ip,jp,kp))              
     *+xf*yf1*(vel(ip1,jp,kp1)-vel(ip1,jp,kp))                
     *+xf1*yf*(vel(ip,jp1,kp1)-vel(ip,jp1,kp))                
     *+xf*yf*(vel(ip1,jp1,kp1)-vel(ip1,jp1,kp)))/zd           
c                                                             
      return                                                  
      end                                                     

c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine aztoa(x1,x2,y1,y2,z1,z2,xr,yr,azim,tkofan)
c  this subroutine computes the azimuth and take-off-angle
c  for an individual observation
c  pr is station, p1 and p2 are hypocenter and adjoining point
c  on the raypath
c
      doubleprecision xltkm, xlnkm, rota
      doubleprecision snr, csr, dxlt, dxln
      common/shrotd/ xltkm, xlnkm, rota, dxlt, dxln, snr, csr

      parameter (drad=1.7453292d-02)
      parameter (pi=3.1415926536)
      parameter (twopi=6.2831853072)
c
c  Azimuth
c
c  start cht 1998
c
      xd=x2-x1
      yd=y2-y1
c
c  end cht 1998
c
      xda=abs(xd)
      yda=abs(yd)
      phi=atan(xda/yda)
      phi2 = atan2(xd,yd)
      if(phi2.lt.0) then
	theta = twopi + phi2
      else
	theta = phi2
      endif
c  compute correct azimuth depending on quadrant
c     if(xd.ge.0.0) then
c       if(yd.ge.0.0) then
c         theta=twopi-phi
c       else
c        theta=pi+phi
c       endif
c     else
c       if(yd.ge.0.0) then                                    
c         theta=phi                                           
c       else                                                  
c         theta=pi-phi                                        
c       endif                                                 
c     endif                                                   
c  rotate back to real north, convert to degrees              
      azim=(theta-rota)/drad                                  
c  22aug95 quick change to get correct azimuths for New Zealand
c          S lat and E long                                   
c     if(nzco.eq.1) azim=azim+180                             
      if(azim.gt.360.0) azim=azim-360.0                       
      if(azim.lt.0.0) azim=azim+360.0                         
c                                                             
c  Take-off-angle                                             
      xd=x2-x1                                                
      yd=y2-y1                                                
      r=sqrt(xd*xd+yd*yd)                                     
      zd=z2-z1                                                
      zda=abs(zd)                                             
      phi=atan(r/zda)                                         
      if(zd.lt.0.0) phi=pi-phi                                
      tkofan=phi/drad                                         
c                                                             
      return                                                  
      end   
