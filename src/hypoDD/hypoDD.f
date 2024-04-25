      program hypoDD

c VERSION: 2.1beta  -  2016/06/15

c AUTHOR:  Felix Waldhauser 
c          Lamont-Doherty Earth Observatory, Columbia University
c 	   Palisades, NY 10964
c 	   felixw@ldeo.columbia.edu
c
c PURPOSE:
c HypoDD is a Fortran77 program for high-resolution relocation of earthquakes 
c using the double-difference algorithm of Waldhauser and Ellsworth (2000). 
c It incorporates phase-pick and/or cross-correlation P- and/or S-wave 
c relative travel-time measurements. Residuals between observed and theoretical c travel time differences (or double-differences = DD) are minimized for pairs
c of earthquakes at each station while linking together all observed
c event/station pairs. A weighted least squares solution (SVD or LSQR) is found
c by iteratively adjusting the vector difference between hypocentral pairs.
c Layered 1D or 3D velocity models are used to solve the forward problem.
c 
c Main changes to hypoDD 1.3:
c - Use of station elevation for travel time prediction
c - Variable vp/vs ratios in layered 1D models (imod=1)
c - Station specific 1D models (imod=4)
c - Negative station elevation (incl stations below sources) in constant 
c   velocity model (imod=5)
c - use of 3D velocity models using the Simulps ray tracer of Um and 
c   Thurber (1983) (originally implemented by Andreas Rietbrock) (imod=9)
c - Improved weighting scheme and handling of airquakes (iaq)
c
c User's Guide to installing and running hypoDD is included in the
c package (hypoDD2_UserGuide.pdf; Waldhauser, 2001).

C DISCLAIMERS, WARRANTIES
c The program has been extensively tested, but is made available
c without warranty.  Users of the program are free to make modifications
c to the programs to meet their particular needs, but are asked not to
c distribute modified code to others. Tell me about your improvements
c and bug fixes, and I will add them to the code for everyone's benefit.
c Do not remove this header from the program. Please use appropriate 
c references and acknowledments when this software is useful to you. 
c Copyright (C) 2011 Felix Waldhauser

c DOWNLOADS:
c The complete hypoDD program package can be downloaded from 
c http://www.ldeo.columbia.edu/~felixw/hypoDD.html
c The package comes with a Users Guide (Waldhauser, 2001), instructions 
c how to install hypoDD, and several example data sets.
c
c REFERENCES:
c For a detailed description of the algorithm see:
c    Waldhauser, F. and W.L. Ellsworth, A double-difference earthquake
c    location algorithm: Method and application to the northern Hayward
c    fault, Bull. Seismol. Soc. Am., 90, 1353-1368, 2000.
c For a user guide to hypoDD see the USGS open-file report
c    Waldhauser, F., HypoDD: A computer program to compute double-difference
c    earthquake locations,  U.S. Geol. Surv. open-file report , 01-113,
c    Menlo Park, California, 2001.
c and updates thereof.

c UPDATES ALERT:
c Send email to felixw@ldeo.columbia.edu with 'HYPODD UPDATE' included in the 
c subject line.  

c VERSION HISTORY:
c Version 1.0 - 03/2001	  (written while at USGS) 
c Version 1.1 - 10/2004 
c Version 1.2 - 07/2010 
c Version 1.3 - 08/2010   (same as v1.2, but compiles with gfortran)
c Version 2.0b - 05/2011   
c Version 2.1b - 08/2011,09/2012/02/2012,06/2012,08/2012,01/2013,11/2015,06/2016

c LOG OF CHANGES:
c 15/06/2016  partials_1dsr.f: fixed m/km bug (IMOD=5) (cfw160615). 
c 19/11/2015  getdata: fixed call to partials during synthetics run (cfw151119)
c 01/30/2013  fix 60 sec origin time output precision problem (cfw130130)
c 08/21/2012  change mag format from f4.1 to f5.2 
c 06/03/2012  checking now for 3d model array dimensions in get_vel3d.f 
c 02/17/2012  fixed line declaration in nval.f
c 09/20/2011  changes to hypoDD.src output format in partials and 
c             partials_1dsr made. 
c 08/20/2011  VErsion 2.1b (beta) 
c	      - included imod=4, station specific 1D models
c 	      - included imod=5, straight path raytracer that allows for neg.
c               station elevation (incl station below source). 
c 	      - fixed error w/ constrain on mean centroid shift in lsfir_svd
c 	      - included nval.f routine to count number of values in a 
c               given string (line) 
c	      - IAQ parameter activated (remove/keep airquakes).
c	      - problem w/ cleaning sta_mod array in skip.f fixed.
c              -changes to Makefile 
c
c 08/2011     getinp: fix loop error in assigning mod_ratio (made hypoDD crash 
c             on certain linux comps.
c 06/2011     getdata: new code to read in station list (remove sscan functions)
c 05/2011     Version 2.0b (beta)
c               - 3D ray tracing (simulps ray tracer, originally implemented 
c                 by Andreas Rietbrock)  (imod=9)
c 		- Station elevation included (imod=0 and 1) 
c               - Compatiblity with version 1 input files added
c 04/2011     initiate sta_dist and sta_az to avoid * output when run on MacOSX
c 11/2010     fixed /0 problem in rms reporting (NaN in hypoDD.reloc) (c101116)
c 08/2010     version 1.3
c 08/2010     now compiles with gfortran  (rcs removed and mod problem fixed)
c 07/2010     version 1.2
c 07/2010     lsfit_svd: Fixed apparent bug in getting 95% confidence errors
c             from standard errors. Factor 2.7955 was used, but it should be
c             1.96, assuming a t-distribution of the residuals 
c             (see email by Hilary Martens)
c 07/2010     fixed bug in computing rms values for hypoDD.reloc 
c             (see email by Zhonhe Zhao)
c 06/2007     accomodate negative magnitudes in output format.
c             real -> doubleprecision: covar,lsfit_svd,matmult2,matmutl3,svd
c 06/2007     Version 1.2 started. fixed errors listed in Buglist 1.1
c 10/2004     Version 1.1: fixed errors listed in BugList to V 1.0.
c 03/2001     Version 1.0 
c 03/2001     clean up & bug fixes by Waldhauser, Julian, Klein, Richards-Dinger
c started 03/1999 
c-----------------------------------------------------------------------------

	implicit none

	include'hypoDD.inc'

	real		acond
	real		adamp(10)
	real		adep
	integer		aiter(0:10)
	real		alat
	real		alon
	real		amaxdcc(10)
	real		amaxdct(10)
	real		amaxres_cross(10)
	real		amaxres_net(10)
	integer		amcusp(1000)
	real		awt_ccp(10)
	real		awt_ccs(10)
	real		awt_ctp(10)
	real		awt_cts(10)
	integer		clust(MAXCL,MAXEVE)
	real		cohav
	real		damp
	character	dattim*25
	real		dtav
	integer		dt_c1(MAXDATA)
	integer		dt_c2(MAXDATA)
	real		dt_cal(MAXDATA)
	real		dt_dt(MAXDATA)
	integer		dt_ic1(MAXDATA)
	integer		dt_ic2(MAXDATA)
	integer		dt_idx(MAXDATA)
	integer		dt_ista(MAXDATA)
	real		dt_offs(MAXDATA)
	real		dt_qual(MAXDATA)
	real		dt_res(MAXDATA)
	character	dt_sta(MAXDATA)*7
	real		dt_wt(MAXDATA)
	real		dxav
	real		dyav
	real		dzav
	real		etav
	integer		ev_cusp(MAXEVE)
	integer		ev_date(MAXEVE)
	real		ev_dep(MAXEVE)
	real		ev_herr(MAXEVE)
	real		ev_lat(MAXEVE)
	real		ev_lon(MAXEVE)
	real		ev_mag(MAXEVE)
	real		ev_res(MAXEVE)
	integer		ev_time(MAXEVE)
	real		ev_x(MAXEVE)
	real		ev_y(MAXEVE)
	real		ev_zerr(MAXEVE)
	real		ev_z(MAXEVE)
        integer         ev_fix(MAXEVE)
	logical		ex
	real		exav
	real		eyav
	real		ezav
	character	fn_cc*110
	character	fn_ct*110
	character	fn_eve*110
	character	fn_inp*110
	character	fn_loc*110
	character	fn_reloc*110
	character	fn_res*110
	character	fn_srcpar*110
	character	fn_sta*110
	character	fn_stares*110
	integer		fu0
	integer		fu1
	integer		fu3
	integer		i
	integer		iargc
	integer		ibeg
	integer		iclust
	integer		icusp(MAXEVE)
	integer		idata
	integer		idy
	integer		iend
	integer		ihr
	integer		imn
	integer		imo
	integer		ineg
	integer		iphase
	integer		isolv
	integer		istart
	integer		iter
	integer		itf
	integer		iunit
	integer		iyr
	integer		j
	integer		jiter
	integer		juliam
	integer		k
	integer		kiter
	integer		l
	doubleprecision	lat
	integer		log
	doubleprecision	lon
	real		maxdcc
	real		maxdct
	real 		maxdist
	integer		maxiter
	real		maxres_cross
	real		maxres_net
	integer		mbad
	integer		minobs_cc
	integer		minobs_ct
	real		minwght
	integer		mod_nl
	real		mod_ratio(MAXLAY)
	real		mod_top(MAXLAY)
	real		mod_v(MAXLAY)
	integer		n
	integer		narguments
	integer		ncc
	integer		nccold
	integer		nccp
	integer		nccs
	integer		nclust
	integer		nct
	integer		nctold
	integer		nctp
	integer		ncts
	integer		ncusp
	integer		ndt
	integer		nev
	integer		nevold
	integer		niter
	integer		noclust(MAXEVE)
	real		noisef_dt
	integer		nsrc
	integer		nsta
	real		picav
	real		resvar1
	real		rms_cc
	real		rms_cc0
	real		rms_cc0old
	real		rms_ccold
	real		rms_ct
	real		rms_ct0
	real		rms_ct0old
	real		rms_ctold
	real		sc
	real		sdc0_dep
	real		sdc0_lat
	real		sdc0_lon
	integer		src_cusp(MAXEVE)
	real		src_dep(MAXEVE)
	real		src_dt(MAXEVE)
	real		src_dx(MAXEVE)
	real		src_dy(MAXEVE)
	real		src_dz(MAXEVE)
	real		src_et(MAXEVE)
	real		src_ex(MAXEVE)
	real		src_ey(MAXEVE)
	real		src_ez(MAXEVE)
	real		src_lat0(MAXEVE)
	doubleprecision	src_lat(MAXEVE)
	real		src_lon0(MAXEVE)
	doubleprecision	src_lon(MAXEVE)
	integer		src_nnp(MAXEVE)
	integer		src_nns(MAXEVE)
	integer		src_np(MAXEVE)
	integer		src_ns(MAXEVE)
	real		src_rmsc(MAXEVE)
	real		src_rmsn(MAXEVE)
	real		src_t0(MAXEVE)
	real		src_t(MAXEVE)
	real		src_x0(MAXEVE)
	real		src_x(MAXEVE)
	real		src_y0(MAXEVE)
	real		src_y(MAXEVE)
	real		src_z0(MAXEVE)
	real		src_z(MAXEVE)
	real		sta_az(MAXSTA)
	real		sta_dist(MAXSTA)
	character	sta_lab(MAXSTA)*7
	real		sta_lat(MAXSTA)
	real		sta_lon(MAXSTA)
	real		sta_elv(MAXSTA)
        integer         sta_mod(MAXSTA)
	integer		sta_nnp(MAXSTA)
	integer		sta_nns(MAXSTA)
	integer		sta_np(MAXSTA)
	integer		sta_ns(MAXSTA)
	real		sta_rmsc(MAXSTA)
	real		sta_rmsn(MAXSTA)
	character	str1*60
	character	str80*80
	character	str3*3
	real		tav
	real		tmpr1
	real		tmpr2
	real		tmp_ttp(MAXSTA,MAXEVE)
	real		tmp_tts(MAXSTA,MAXEVE)
	real		tmp_xp(MAXSTA,MAXEVE)
	real		tmp_yp(MAXSTA,MAXEVE)
	real		tmp_zp(MAXSTA,MAXEVE)
	real		tmp_xs(MAXSTA,MAXEVE)
	real		tmp_ys(MAXSTA,MAXEVE)
	real		tmp_zs(MAXSTA,MAXEVE)
	integer		trimlen
	real		wt_ccp
	real		wt_ccs
	real		wt_ctp
	real		wt_cts
	real		x
	real		xav
	real		y
	real		yav
	real		zav
c 	3d raytracing parameters:
	real		lat_3d
	real		lon_3d
        real		rot_3d
        character	fn_mod3d*110
        character	fn_mod1d*110
    	integer		ndip
        integer		iskip
	integer		scale1
	integer		scale2
	real		xfac
	real		tlim
	integer		nitpb

c	v2 parameters:
        integer         fu_inp
        character       str8*8
        character       str9*9	!ccc
        integer		iwrite
        real 		minds
	real 		maxds
	real		maxgap
        integer		iaq
        integer		imod
        integer		ipha3d

      
      minwght= 0.00001
      rms_ccold= 0
      rms_ctold= 0
      rms_cc0old= 0
      rms_ct0old=  0
c--- open log file:
      call freeunit(log)
      open(log,file='hypoDD.log',status='unknown')
      str1= 'starting hypoDD (v2.1beta - 06/15/2016)...'
      call datetime(dattim)
      write(6,'(a45,a25)') str1, dattim
      write(log,'(a45,a25)') str1, dattim

c--- get input parameter file name:
      narguments = iargc()
      if(narguments.lt.1) then
        write(*,'(/,a)') 'PARAMETER INPUT FILE [hypoDD.inp <ret>]:'
        read(5,'(a)') fn_inp
        if(trimlen(fn_inp).le.1) then
           fn_inp= 'hypoDD.inp'            !default input file name
        else
           fn_inp= fn_inp(1:trimlen(fn_inp))
        endif
      else
          call getarg(1,fn_inp)
      endif
      inquire(FILE= fn_inp,exist=ex)
      if(.not. ex) stop' >>> ERROR OPENING INPUT PARAMETER FILE.'

c--- get input parameters:
      call freeunit (fu_inp)
      open (fu_inp,status='unknown',file=fn_inp,err=998)
ccc      read (fu_inp,'(8a)') str8
      read (fu_inp,'(9a)') str9
      close(fu_inp)
ccc      if (str8.eq.'hypoDD_2') then
      if (str9.eq.'hypoDD_v2' .or. str9.eq.'hypoDD_2') then
         call getinp2(MAXEVE,MAXLAY,log,fn_inp,
     &   fn_cc,fn_ct,fn_sta,fn_eve,
     &   fn_loc,fn_reloc,fn_res,fn_stares,fn_srcpar,
     &   iwrite,
     &   idata,iphase,
     &   minobs_cc,minobs_ct,
     &   amaxres_cross,amaxres_net,amaxdcc,amaxdct,
     &   noisef_dt,maxdist,minds,maxds,maxgap,
     &   awt_ccp,awt_ccs,awt_ctp,awt_cts,adamp,
     &   istart,maxiter,isolv,iaq,niter,aiter,
     &   imod,mod_nl,mod_ratio,mod_v,mod_top,
     &   fn_mod1d,fn_mod3d, lat_3d, lon_3d, rot_3d,
     &   ipha3d,ndip,iskip,scale1,scale2,xfac,tlim,nitpb,
     &   iclust,ncusp,icusp)
c         minds= -999 	!not implemented
c         maxds= -999 	!not implemented
c         maxgap= -999 	!not implemented
      else
         call getinp(MAXEVE,MAXLAY,log,fn_inp,
     &    fn_cc,fn_ct,fn_sta,fn_eve,
     &    fn_loc,fn_reloc,fn_res,fn_stares,fn_srcpar,
     &    idata,iphase,
     &    minobs_cc,minobs_ct,
     &    amaxres_cross,amaxres_net,amaxdcc,amaxdct,
     &    noisef_dt,maxdist,
     &    awt_ccp,awt_ccs,awt_ctp,awt_cts,adamp,
     &    istart,maxiter,isolv,niter,aiter,
     &    mod_nl,mod_ratio,mod_v,mod_top,
     &    iclust,ncusp,icusp)
         imod= 0
         minds= -999
         maxds= -999
         maxgap= -999
         iaq= 1         ! remove airquakes as in hypoDD1.0
      endif

c--- get data:
      call getdata(
     & log,fn_cc,fn_ct,fn_sta,fn_eve,fn_srcpar,
     & imod,idata,iphase,ncusp,icusp,
     & maxdist,minds,maxds,maxgap,amaxdct(1),amaxdcc(1),
     & noisef_dt,mod_nl,mod_ratio,mod_v,mod_top,
     & ev_date,ev_time,ev_cusp,ev_lat,ev_lon,ev_dep,
     & ev_mag,ev_herr,ev_zerr,ev_res,ev_fix,
     & sta_lab,sta_lat,sta_lon,sta_elv,sta_mod,
     & dt_sta,dt_dt,dt_qual,dt_c1,dt_c2,dt_idx,
     & dt_ista,dt_ic1,dt_ic2,dt_offs,
     & nev,nsta,ndt,nccp,nccs,nctp,ncts,
     & tmp_xp,tmp_yp,tmp_zp,tmp_ttp,tmp_tts)

c--- get 3D velocity model if necessary
      if(imod.eq.9) then
c       i= 2
c       call get_vel3d(fn_mod3d,i,
       call get_vel3d(fn_mod3d,ipha3d,
     &  lat_3d,lon_3d,rot_3d,
     &  ndip,iskip,scale1,scale2,xfac,tlim,nitpb,
     &  log)
      endif

c--- clustering:
      if((idata.eq.1.and.minobs_cc.eq.0).or.
     &   (idata.eq.2.and.minobs_ct.eq.0).or.
     &   (idata.eq.3.and.minobs_ct+minobs_cc.eq.0)) then
         nclust= 1
         clust(1,1)= nev
         do i=1,nev
             clust(1,i+1)= ev_cusp(i)
         enddo
          write(*,'(/,"no clustering performed.")')
          write(log,'(/,"no clustering performed.")')
      else

         call cluster1(log, nev, ndt,
     & idata, minobs_cc, minobs_ct,
     & dt_c1, dt_c2, ev_cusp,
     & clust, noclust, nclust)

      endif

c--- open files
      call freeunit(fu0)
      open(fu0,file=fn_loc,status='unknown')
      call freeunit(fu1)
      open(fu1,file=fn_reloc,status='unknown')
      if(trimlen(fn_stares).gt.1) then
         call freeunit(fu3)
         open(fu3,file=fn_stares,status='unknown')
      endif

      jiter = 0  ! counter for iter with no updating (air quakes)
c--- big loop over clusters starts here:
      if(iclust.ne.0) then
        if (iclust.lt.0 .or. iclust.gt.nclust) then
           write(*,*) 'error: invalid cluster number ',iclust
           write(*,*) 'must be between 1 and nclust (',nclust,')'
           stop
        endif
        ibeg= iclust
        iend= iclust
      else
        ibeg= 1
        iend= nclust
      endif
      do iclust= ibeg,iend
      call datetime(dattim)
      write(log,'(/,"RELOCATION OF CLUSTER:",i2,5x,a25,/,
     &"----------------------")')iclust,dattim
      write(*,'(/,"RELOCATION OF CLUSTER:",i2,5x,a25,/,
     &"----------------------")')iclust,dattim

c--- get data for each cluster if clustering was invoked:
      if((nclust.eq.1.and.clust(iclust,1).eq.nev).or.
     &   (idata.eq.1.and.minobs_cc.eq.0).or.
     &   (idata.eq.2.and.minobs_ct.eq.0).or.
     &   (idata.eq.3.and.minobs_ct+minobs_cc.eq.0)) goto 50

      ncusp= clust(iclust,1)
      do i=1,ncusp
         icusp(i)= clust(iclust,i+1)
      enddo

      if(idata.ne.0) call getdata(
     & log,fn_cc,fn_ct,fn_sta,fn_eve,fn_srcpar,
     & imod, idata,iphase,ncusp,icusp,
     & maxdist,minds,maxds,maxgap,amaxdct(1),amaxdcc(1),
     & noisef_dt,mod_nl,mod_ratio,mod_v,mod_top,
     & ev_date,ev_time,ev_cusp,ev_lat,ev_lon,ev_dep,
     & ev_mag,ev_herr,ev_zerr,ev_res,ev_fix,
     & sta_lab,sta_lat,sta_lon,sta_elv,sta_mod,
     & dt_sta,dt_dt,dt_qual,dt_c1,dt_c2,dt_idx,
     & dt_ista,dt_ic1,dt_ic2,dt_offs,
     & nev,nsta,ndt,nccp,nccs,nctp,ncts,
     & tmp_xp,tmp_yp,tmp_zp,tmp_ttp,tmp_tts)

50    continue
      nccold= nccp+nccs
      nctold= nctp+ncts
      ncc= nccp+nccs
      nct= nctp+ncts
      nevold= nev

c--- get cluster centroid, or 3D model origin:
      sdc0_lat= 0
      sdc0_lon= 0
      sdc0_dep= 0
      do i=1,nev
         sdc0_lat= sdc0_lat + ev_lat(i)
         sdc0_lon= sdc0_lon + ev_lon(i)
         sdc0_dep= sdc0_dep + ev_dep(i)
      enddo
      sdc0_lat= sdc0_lat/nev
      sdc0_lon= sdc0_lon/nev
      sdc0_dep= sdc0_dep/nev
      
      if(imod.eq.9) then
         sdc0_lat= lat_3d    
         sdc0_lon= lon_3d    
      else
         rot_3d= 0.0 
      endif

      write(log,'("Cluster centroid at:",1x,f10.6,2x,f11.6,2x,f9.6)')
     & sdc0_lat,sdc0_lon,sdc0_dep
      write(log,'("Rotation of coordinate system:",1x,f10.6)')
     & rot_3d 

c--- Set up cartesian coordinate system (build common block for subr. SDC):
      call setorg(sdc0_lat,sdc0_lon,rot_3d,0)
c      call setorg(sdc0_lat,sdc0_lon,0.0,0)

c--- get cartesian coordinates for epicenters
      do i=1,nev
         lat= ev_lat(i)
         lon= ev_lon(i)
         call SDC2(x,y,lat,lon,-1)
         ev_x(i)= x *1000
         ev_y(i)= y *1000
         ev_z(i)= (ev_dep(i) - sdc0_dep)*1000
      enddo

      write(log,'("# events:",i5)')nev

c--- write output (mdat.loc):
      write(fu0,'(i9,1x,f10.6,1x,f11.6,1x,f9.3,1x,f10.1,1x,f10.1,
     & 1x,f10.1,
     & 1x,f8.1,1x,f8.1,1x,f8.1,1x,i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f5.2,
     & 1x,f5.2,1x,i3)')
cfw     & 1x,f4.1,1x,i3)')  20120821
cfw     & 1x,f3.1,1x,i3)')   !neg mag format
     & (ev_cusp(i),ev_lat(i),ev_lon(i),ev_dep(i),ev_x(i),ev_y(i),
     & ev_z(i),ev_herr(i)*1000,ev_herr(i)*1000,ev_zerr(i)*1000,
     & int(ev_date(i)/10000),
     & int(mod(ev_date(i),10000)/100),mod(ev_date(i),100),
     & int(ev_time(i)/1000000),int(mod(ev_time(i),1000000)/10000),
     & mod(dble(ev_time(i)),10000.)/100,ev_mag(i),iclust,i=1,nev)
cfw130130     & mod(real(ev_time(i)),10000.)/100,ev_mag(i),iclust,i=1,nev)
cfw100806     & mod(real(ev_time(i)),10000)/100,ev_mag(i),iclust,i=1,nev)

c--- get initial trial sources:
      call trialsrc(istart,sdc0_lat,sdc0_lon,sdc0_dep,
     & nev,ev_cusp,ev_lat,ev_lon,ev_dep,
     & nsrc,src_cusp,src_lat0,src_lon0,
     & src_x0,src_y0,src_z0,src_t0,
     & src_lat,src_lon,src_dep,
     & src_x,src_y,src_z,src_t)

      write(*,'("Initial trial sources =",i6)')nsrc
      write(log,'("# initial trial sources:",i6)')nsrc

c--- loop over iterations starts here:
c define iteration step at which re-weighting starts: this is dynam. since
c it depends on the number of neg depths runs before.

c first reset aiter() and maxiter
      do i=1,niter
         aiter(i) = aiter(i) - jiter
      enddo
      maxiter = maxiter - jiter

      kiter= 0		! counter for iter with data skipping
      jiter= 0		! counter for iter with no updating (air quakes)
      mbad= 0		! counter for air quakes

      iter= 1
55    call datetime(dattim)
      write(log,'(/,"===ITERATION ",i3," (",i3,") ",a25)')
     & iter-jiter, iter, dattim

c--- get weighting parameters for this iteration:
      do i=1,niter
        if(iter.le.aiter(i)) goto 75
      enddo
75    maxres_cross= amaxres_cross(i)
      maxres_net= amaxres_net(i)
      maxdcc= amaxdcc(i)
      maxdct= amaxdct(i)
      wt_ccp= awt_ccp(i)
      wt_ccs= awt_ccs(i)
      wt_ctp= awt_ctp(i)
      wt_cts= awt_cts(i)
      damp= adamp(i)

      write(log, '(/,"Weighting parameters for this iteration:",/,
     &"  wt_ccp= ",f7.4,2X,"wt_ccs= ",f7.4,2X,
     &"maxr_cc= ",f7.4,2X,"maxd_cc= ",f7.2,2X,/,
     &"  wt_ctp= ",f7.4,2x,"wt_cts= ",f7.4,2x,"maxr_ct= ",f7.4,2x,
     &"maxd_ct= ",f7.2,/,"  damp= ",f5.1)')
     & wt_ccp,wt_ccs,maxres_cross,
     & maxdcc,wt_ctp,wt_cts,maxres_net,
     & maxdct,damp

c--- calculate travel times  and slowness vectors:
      write(log,'(/,"~ getting partials for ",i5,
     & " stations and ",i5," source(s) ...")') nsta,nsrc
      if(imod.eq.0.or.imod.eq.1) then
          if(iter.eq.1) write(*,'(a)')'1D ray tracing.'
          call partials(fn_srcpar,
     &     nsrc,src_cusp,src_lat,src_lon,src_dep,
     &     nsta,sta_lab,sta_lat,sta_lon,sta_elv,
     &     mod_nl,mod_ratio,mod_v,mod_top,
     &     tmp_ttp,tmp_tts,
     &     tmp_xp,tmp_yp,tmp_zp,tmp_xs,tmp_ys,tmp_zs)

      elseif(imod.eq.5) then
          if(iter.eq.1) write(*,'(a)')
     &    'Ray tracing along straight ray paths.'
          call partials_1dsr(fn_srcpar,
     &     nsrc,src_cusp,src_lat,src_lon,src_dep,
     &     nsta,sta_lab,sta_lat,sta_lon,sta_elv,
     &     mod_nl,mod_ratio,mod_v,mod_top,
     &     tmp_ttp,tmp_tts,
     &     tmp_xp,tmp_yp,tmp_zp,tmp_xs,tmp_ys,tmp_zs)

      elseif(imod.eq.4) then
          if(iter.eq.1) write(*,'(a)')
     &    '1D ray tracing in station specific models.'
          call partials_1dmm(log,fn_srcpar,fn_mod1d,
     &     nsrc,src_cusp,src_lat,src_lon,src_dep,
     &     nsta,sta_lab,sta_lat,sta_lon,sta_elv,sta_mod,
     &     mod_nl,mod_ratio,mod_v,mod_top,iter,
     &     tmp_ttp,tmp_tts,
     &     tmp_xp,tmp_yp,tmp_zp,tmp_xs,tmp_ys,tmp_zs)

      elseif(imod.eq.9) then
          if(iter.eq.1) write(*,'(a)')'3D ray tracing.'
c 	  Determine source-station paths for ray-tracing
          do i=1,nsta
             do j=1,nsrc
                tmp_ttp(i,j)= -999
             enddo
          enddo
          do k=1,ndt
             do j=1,nsrc
                if(dt_c1(k).eq.src_cusp(j)) tmp_ttp(dt_ista(k),j)= 0
                if(dt_c2(k).eq.src_cusp(j)) tmp_ttp(dt_ista(k),j)= 0
             enddo
          enddo
    
          call partials_3d(fn_srcpar,
     &     nsrc,src_cusp,src_lat,src_lon,src_dep,
     &     nsta,sta_lab,sta_lat,sta_lon,sta_elv,
     &     mod_nl,mod_ratio,mod_v,mod_top,rot_3d,ipha3d,
     &     tmp_ttp,tmp_tts,
     &     tmp_xp,tmp_yp,tmp_zp,
     &     tmp_xs,tmp_ys,tmp_zs)
    
      else
          write(*,'(a)')"IMOD no valid."
          stop
      endif

c--- get double difference vector:
      call dtres(log,ndt,MAXSTA,nsrc,
     & dt_dt,dt_idx,
     & dt_ista,dt_ic1,dt_ic2,
     & src_cusp,src_t,tmp_ttp,tmp_tts,
     & dt_cal,dt_res)

c--- get a priori weights and reweight residuals
      call weighting(log,ndt,mbad,amcusp,idata,kiter,ineg,
     &               maxres_cross,maxres_net,maxdcc,maxdct,minwght,
     &               wt_ccp,wt_ccs,wt_ctp,wt_cts,
     &               dt_c1,dt_c2,dt_idx,dt_qual,dt_res,dt_offs,
     &               dt_wt)

c--- skip outliers and/or air quakes:
      if(ineg.gt.0) then
          call skip(log,kiter,minwght,
     & ndt,nev,nsrc,nsta,
     & ev_cusp,ev_date,ev_time,ev_mag,
     & ev_lat,ev_lon,ev_dep,ev_x,ev_y,ev_z,
     & ev_herr,ev_zerr,ev_res,ev_fix,
     & src_cusp, src_lat, src_lon, src_dep,
     & src_lat0, src_lon0,
     & src_x, src_y, src_z, src_t, src_x0, src_y0, src_z0, src_t0,
     & sta_lab,sta_lat,sta_lon,sta_elv,sta_mod,sta_dist,sta_az,
     & sta_rmsc,sta_rmsn,sta_np,sta_ns,sta_nnp,sta_nns,
     & dt_sta,dt_c1,dt_c2,dt_idx,dt_dt,dt_qual,dt_cal,
     & dt_ista,dt_ic1,dt_ic2,
     & dt_res,dt_wt,dt_offs,
     & tmp_ttp,tmp_tts,tmp_xp,tmp_yp,tmp_zp,tmp_xs,tmp_ys,tmp_zs,
     & nct,ncc)

c--Dont mess anymore with this cluster if we have wiped out all events
        if(nev.lt.2) then 
          write(log,*)' Cluster has less than 2 events.'
          write(*,*)' Cluster has less than 2 events.'
          goto 778
        endif
      else
         write(log,'("no data skipped.")')

      endif

c--- get initial residual statistics (avrg,rms,var..)
      if(iter.eq.1) then
       resvar1= -999
       call resstat(log,idata,ndt,nev,dt_res,dt_wt,dt_idx,
     & rms_cc,rms_ct,rms_cc0,rms_ct0,
     & rms_ccold,rms_ctold,rms_cc0old,rms_ct0old,
     &              resvar1)
      endif

c--- least square fitting:

      if(isolv.eq.1) then
         call lsfit_SVD(log,iter,ndt,nev,nsrc,damp,mod_ratio,
     & idata,ev_cusp,src_cusp,ev_fix,
     & dt_res,dt_wt,
     & dt_ista,dt_ic1,dt_ic2,   !new
     & src_dx,src_dy,src_dz,src_dt,src_ex,src_ey,src_ez,src_et,
     & exav,eyav,ezav,etav,dxav,dyav,dzav,dtav,
     & rms_cc,rms_ct,rms_cc0,rms_ct0,
     & rms_ccold,rms_ctold,rms_cc0old,rms_ct0old,
     & tmp_xp,tmp_yp,tmp_zp,
     & tmp_xs,tmp_ys,tmp_zs,dt_idx)

      else
         call lsfit_lsqr(log,iter,ndt,nev,nsrc,damp,mod_ratio,
     & idata,ev_cusp,src_cusp,ev_fix,
     & dt_res,dt_wt,
     & dt_ista,dt_ic1,dt_ic2,   !new
     & src_dx,src_dy,src_dz,src_dt,src_ex,src_ey,src_ez,src_et,
     & exav,eyav,ezav,etav,dxav,dyav,dzav,dtav,
     & rms_cc,rms_ct,rms_cc0,rms_ct0,
     & rms_ccold,rms_ctold,rms_cc0old,rms_ct0old,
     & tmp_xp,tmp_yp,tmp_zp,
     & tmp_xs,tmp_ys,tmp_zs,dt_idx,acond)
      endif

c--- check for air quakes:
      mbad= 0
      k= 1
      do i= 1,nsrc
        if(src_dep(i) + (src_dz(i)/1000).lt.0) then
            write(log,'(">>>Warning: negative depth - ",i12)')ev_cusp(i)
            amcusp(k)= ev_cusp(i)
            k=k+1
            if(k.gt.1000) stop'>>> More than 1000 air quakes. Too many!'
               src_dz(i)= 0
        endif
c      fw-06/10/18: reset to original depth after each iteration if ev_fix= -1.
        if(ev_fix(i).eq.-1) then
           src_dz(i)= 0
        endif
      enddo
      mbad= k-1    ! number of neg depth events
      if(iaq.eq.0) mbad= -mbad               ! don't remove air quakes

      if(abs(mbad).gt.0) then
         write(log,*)'Number of air quakes (AQ) =',abs(mbad)
         if(nsrc-mbad .le. 1) then
           write(*,*)'Warning: number of non-airquakes < 2'
                   write(*,*)'   skipping this cluster'
           write(log,*)'Warning: number of non-airquakes < 2'
                   write(log,*)'   skipping this cluster'
                   goto 778
         endif       
      endif       
c update iteration numbers:
      if(mbad.gt.0) then
         do i= 1,niter
              aiter(i)= aiter(i)+1
         enddo
         jiter= jiter+1	  ! iteration with no update
         maxiter= maxiter+1

         goto 500   ! skip the updating step
      endif

c--- update source parameters:
      xav= 0 ! mean centroid shift
      yav= 0
      zav= 0
      tav= 0
      alon= 0
      alat= 0
      adep= 0
      if(nsrc.eq.1) nsrc= nev
      do i= 1,nsrc
        src_cusp(i)= ev_cusp(i)
c update absolute source parameters (cart)
        src_x(i)= src_x(i) + src_dx(i)
        src_y(i)= src_y(i) + src_dy(i)
        src_z(i)= src_z(i) + src_dz(i)
        src_t(i)= src_t(i) + src_dt(i)

c update absolute source locations (geogr)
        src_dep(i)= src_dep(i) + (src_dz(i)/1000)
        call SDC2(src_x(i)/1000,src_y(i)/1000,lat,lon,1)
        src_lon(i)= lon
        src_lat(i)= lat
        alon= lon+alon	
        alat= lat+alat
        adep= adep+src_dep(i)

c get mean centroid shift
        xav= xav + (src_x(i) - src_x0(i))	
        yav= yav + (src_y(i) - src_y0(i))
        zav= zav + (src_z(i) - src_z0(i))
        tav= tav + (src_t(i) - src_t0(i))
      enddo
      xav= xav/nsrc
      yav= yav/nsrc
      zav= zav/nsrc
      tav= tav/nsrc
      alon= alon/nsrc
      alat= alat/nsrc
      adep= adep/nsrc

      write(log,'("  cluster centroid at:",1x,f10.6,2x,f11.6,2x,f9.6)')
     & alat,alon,adep
      write(log,'("  mean centroid (origin) shift in x,y,z,t [m,ms]: ",/
     & f7.1,f7.1,f7.1,f7.1)'),xav,yav,zav,tav
      write(log,'("  (OS in std output gives maximum value.)")')

c--- get interevent distance for each observation and average signal coherency:
      cohav= 0
      picav= 0
      j= nct
      k= ncc
      ncc= 0
      nct= 0
      do i= 1,ndt
         dt_offs(i)= sqrt((src_x(dt_ic1(i))-src_x(dt_ic2(i)))**2 +
     &                  (src_y(dt_ic1(i))-src_y(dt_ic2(i)))**2 +
     &                  (src_z(dt_ic1(i))-src_z(dt_ic2(i)))**2)

         if(dt_idx(i).le.2) then
            cohav= cohav + sqrt(dt_qual(i))
            ncc= ncc+1
         else
            picav= picav + dt_qual(i)
            nct= nct+1
         endif

      enddo
      cohav= cohav/ncc
      picav= picav/nct
      write(log,'(/,"More:")')
      write(log,'("  mean phase coherency = ",f5.3)')cohav
      write(log,'("  mean pick quality = ",f5.3)')picav

c--- get number of observations and mean residual at each station
      tmpr1= 0
      tmpr2= 0
      do i= 1,nsta
         sta_dist(i)= 0
         sta_az(i)= 0
         sta_np(i)= 0
         sta_ns(i)= 0
         sta_nnp(i)= 0
         sta_nns(i)= 0
         sta_rmsc(i)= 0
         sta_rmsn(i)= 0
         do j= 1,ndt
            if(i.eq.dt_ista(j)) then
               if(dt_idx(j).le.2) then
                 sta_rmsc(i)= sta_rmsc(i)+dt_res(j)**2
                 if(dt_idx(j).eq.1) then
                    sta_np(i)= sta_np(i)+1
                 else
                    sta_ns(i)= sta_ns(i)+1
                 endif
               else
                 sta_rmsn(i)= sta_rmsn(i)+dt_res(j)**2
                 if(dt_idx(j).eq.3) then
                   sta_nnp(i)= sta_nnp(i)+1
                 else
                   sta_nns(i)= sta_nns(i)+1
                 endif
               endif
            endif
         enddo

         if(sta_np(i)+sta_ns(i).gt.0)
     &     sta_rmsc(i)= sqrt(sta_rmsc(i)/(sta_np(i)+sta_ns(i)))
         if(sta_nnp(i)+sta_nns(i).gt.0)
     &     sta_rmsn(i)= sqrt(sta_rmsn(i)/(sta_nnp(i)+sta_nns(i)))
         if(sta_rmsc(i).gt.tmpr1) then
            tmpr1= sta_rmsc(i)
            k= i
         endif
         if(sta_rmsn(i).gt.tmpr2) then
            tmpr2= sta_rmsn(i)
            l= i
         endif
      enddo
      tmpr1= tmpr1*1000
      tmpr2= tmpr2*1000
      if(idata.eq.1.or.idata.eq.3) then
         write(log,'("  station with largest cc rms: ",a7,"=",
     & f7.0," ms (RMSST)")')
     &    sta_lab(k),tmpr1
      endif
      if(idata.eq.2.or.idata.eq.3) then
         write(log,'("  station with largest ct rms: ",a7,"=",
     & f7.0," ms (RMSST)")')
     &    sta_lab(l),tmpr2
      endif

c--- write output scratch mdat.reloc:
      n= trimlen(fn_reloc)
      i=iter-jiter
      write(str80,'(a,".",i3.3,".",i3.3)')fn_reloc(1:n),iclust,i
      call freeunit(iunit)
      open(iunit,file=str80,status='unknown')
      call setorg(sdc0_lat,sdc0_lon,0.0,0)
      do i=1,nev
         call SDC2(x,y,src_lat(i),src_lon(i),-1)
         write(iunit,'(i9,1x,f10.6,1x,f11.6,1x,f9.3,1x,f10.1,1x,f10.1,
     &    1x,f10.1,
     &    1x,f8.1,1x,f8.1,1x,f8.1,1x,i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f6.3,
     &    1x,f5.2,1x,i3)')
cfw     &    1x,f4.1,1x,i3)') 20120821
cfw     & 1x,f3.1,1x,i3)')            ! negative mag
     &    src_cusp(i),src_lat(i),src_lon(i),src_dep(i),x*1000,y*1000,
     &    src_z(i),src_ex(i),src_ey(i),src_ez(i),int(ev_date(i)/10000),
     &    int(mod(ev_date(i),10000)/100),mod(ev_date(i),100),
     &    int(ev_time(i)/1000000),int(mod(ev_time(i),1000000)/10000),
     &    mod(dble(ev_time(i)),10000.)/100,ev_mag(i),iclust
cfw130130     &    mod(real(ev_time(i)),10000.)/100,ev_mag(i),iclust
cfw100806     & mod(real(ev_time(i)),10000)/100,ev_mag(i),iclust)
      enddo
      call setorg(sdc0_lat,sdc0_lon,rot_3d,0)
      close(iunit)
c WARNING: variable "str" is set to zero value by default
      write(log,'(/,"Relocation results for this iteration are"
     & " stored in ",a)')str80(1:trimlen(str80))

500   continue  ! case of air quakes

c standard output:
      if(mbad.gt.0) then
         str3='   '
      else
         n= iter-jiter
         if(n.lt.1000) write(str3,'(i3)')n
         if(n.lt.100) write(str3,'(1x,i2)')n
         if(n.lt.10) write(str3,'(2x,i1)')n
      endif
      if(isolv.eq.1.and.idata.eq.3) then
       if(iter.eq.1) write(*,'(/,"  IT   EV  CT  CC",
     & "    RMSCT      RMSCC   RMSST   DX   DY   DZ   DT   OS  AQ",/,
     &                           "        %   %   %",
     & "   ms     %   ms     %    ms    m    m    m   ms    m ")')
       write(*,'(i2,a3,3(1x,i3),i5,f6.1,i5,f6.1,i6,4i5,i5,i4)')
     & iter,str3,
     & nint(nev*100./nevold),nint(nct*100.0/nctold),
     & nint(ncc*100.0/nccold),
     & nint(rms_ct*1000),(rms_ct-rms_ctold)*100/rms_ctold,
     & nint(rms_cc*1000),(rms_cc-rms_ccold)*100/rms_ccold,
     & nint(max(tmpr1,tmpr2)),
     & nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     & nint(max(abs(xav),abs(yav),abs(zav))),abs(mbad)

       write(log,'(/,"  IT   EV  CT  CC",
     & "    RMSCT      RMSCC   RMSST   DX   DY   DZ   DT   OS  AQ",/,
     &                           "        %   %   %",
     & "   ms     %   ms     %    ms    m    m    m   ms    m ")')
       write(log,'(i2,a3,3(1x,i3),i5,f6.1,i5,f6.1,i6,4i5,i5,i4)')
     & iter,str3,
     & nint(nev*100./nevold),nint(nct*100.0/nctold),
     & nint(ncc*100.0/nccold),
     & nint(rms_ct*1000),(rms_ct-rms_ctold)*100/rms_ctold,
     & nint(rms_cc*1000),(rms_cc-rms_ccold)*100/rms_ccold,
     & nint(max(tmpr1,tmpr2)),
     & nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     & nint(max(abs(xav),abs(yav),abs(zav))),abs(mbad)
      endif
      if(isolv.eq.1.and.idata.eq.1) then
       if(iter.eq.1) write(*,'(/,"  IT   EV  CC",
     & "    RMSCC   RMSST   DX   DY   DZ   DT   OS  AQ",/,
     &                          "        %   %",
     & "   ms     %    ms    m    m    m   ms    m ")')
       write(*,'(i2,a3,2(1x,i3),i5,f6.1,i6,4i5,i5,i4)')
     & iter,str3,
     & nint(nev*100./nevold),
     & nint(ncc*100.0/nccold),
     & nint(rms_cc*1000),(rms_cc-rms_ccold)*100/rms_ccold,
     & nint(max(tmpr1,tmpr2)),
     & nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     & nint(max(abs(xav),abs(yav),abs(zav))),abs(mbad)

       write(log,'(/,"  IT   EV  CC",
     & "    RMSCC   RMSST   DX   DY   DZ   DT   OS  AQ",/,
     &                          "        %   %",
     & "   ms     %    ms    m    m    m   ms    m ")')
       write(log,'(i2,a3,2(1x,i3),i5,f6.1,i6,4i5,i5,i4)')
     & iter,str3,
     & nint(nev*100./nevold),
     & nint(ncc*100.0/nccold),
     & nint(rms_cc*1000),(rms_cc-rms_ccold)*100/rms_ccold,
     & nint(max(tmpr1,tmpr2)),
     & nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     & nint(max(abs(xav),abs(yav),abs(zav))),abs(mbad)
      endif
      if(isolv.eq.1.and.idata.eq.2) then
       if(iter.eq.1) write(*,'(/,"  IT   EV  CT",
     & "    RMSCT     RST   DX   DY   DZ   DT   OS  AQ",/,
     &                           "        %   %",
     & "   ms     %    ms    m    m    m   ms    m ")')
       write(*,'(i2,a3,2(1x,i3),i5,f6.1,i6,4i5,i5,i4)')
     & iter,str3,
     & nint(nev*100./nevold),
     & nint(nct*100.0/nctold),
     & nint(rms_ct*1000),(rms_ct-rms_ctold)*100/rms_ctold,
     & nint(max(tmpr1,tmpr2)),
     & nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     & nint(max(abs(xav),abs(yav),abs(zav))),abs(mbad)

       write(log,'(/,"  IT   EV  CT",
     & "    RMSCT     RST   DX   DY   DZ   DT   OS  AQ",/,
     &                           "        %   %",
     & "   ms     %    ms    m    m    m   ms    m ")')
       write(log,'(i2,a3,2(1x,i3),i5,f6.1,i6,4i5,i5,i4)')
     & iter,str3,
     & nint(nev*100./nevold),
     & nint(nct*100.0/nctold),
     & nint(rms_ct*1000),(rms_ct-rms_ctold)*100/rms_ctold,
     & nint(max(tmpr1,tmpr2)),
     & nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     & nint(max(abs(xav),abs(yav),abs(zav))),abs(mbad)
      endif

      if(isolv.eq.2.and.idata.eq.3) then
       if(iter.eq.1) write(*,'(/,"  IT   EV  CT  CC",
     & "    RMSCT      RMSCC   RMSST   DX   DY   DZ   DT   ",
     & "OS  AQ  CND",/,
     &                           "        %   %   %",
     & "   ms     %   ms     %    ms    m    m    m   ms   ",
     & " m     ")')
       write(*,'(i2,a3,3(1x,i3),i5,f6.1,i5,f6.1,i6,4i5,i5,i4,i5)')
     & iter,str3,
     & nint(nev*100./nevold),nint(nct*100.0/nctold),
     & nint(ncc*100.0/nccold),
     & nint(rms_ct*1000),(rms_ct-rms_ctold)*100/rms_ctold,
     & nint(rms_cc*1000),(rms_cc-rms_ccold)*100/rms_ccold,
     & nint(max(tmpr1,tmpr2)),
     & nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     & nint(max(abs(xav),abs(yav),abs(zav))),abs(mbad),nint(acond)

       write(log,'(/,"  IT   EV  CT  CC",
     & "    RMSCT      RMSCC   RMSST   DX   DY   DZ   DT   ",
     & "OS  AQ  CND",/,
     &                           "        %   %   %",
     & "   ms     %   ms     %    ms    m    m    m   ms   ",
     & " m     ")')
       write(log,'(i2,a3,3(1x,i3),i5,f6.1,i5,f6.1,i6,4i5,i5,i4,i5)')
     & iter,str3,
     & nint(nev*100./nevold),nint(nct*100.0/nctold),
     & nint(ncc*100.0/nccold),
     & nint(rms_ct*1000),(rms_ct-rms_ctold)*100/rms_ctold,
     & nint(rms_cc*1000),(rms_cc-rms_ccold)*100/rms_ccold,
     & nint(max(tmpr1,tmpr2)),
     & nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     & nint(max(abs(xav),abs(yav),abs(zav))),abs(mbad),nint(acond)
      endif
      if(isolv.eq.2.and.idata.eq.1) then
       if(iter.eq.1) write(*,'(/,"  IT   EV  CC",
     & "    RMSCC   RMSST   DX   DY   DZ   DT   OS  AQ  CND",/,
     &                           "        %   %",
     & "   ms     %    ms    m    m    m   ms    m ")')
       write(*,'(i2,a3,2(1x,i3),i5,f6.1,i6,4i5,i5,i4,i5)')
     & iter,str3,
     & nint(nev*100./nevold),
     & nint(ncc*100.0/nccold),
     & nint(rms_cc*1000),(rms_cc-rms_ccold)*100/rms_ccold,
     & nint(max(tmpr1,tmpr2)),
     & nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     & nint(max(abs(xav),abs(yav),abs(zav))),abs(mbad),nint(acond)

       write(log,'(/,"  IT   EV  CC",
     & "    RMSCC   RMSST   DX   DY   DZ   DT   OS  AQ  CND",/,
     &                           "        %   %",
     & "   ms     %    ms    m    m    m   ms    m ")')
       write(log,'(i2,a3,2(1x,i3),i5,f6.1,i6,4i5,i5,i4,i5)')
     & iter,str3,
     & nint(nev*100./nevold),
     & nint(ncc*100.0/nccold),
     & nint(rms_cc*1000),(rms_cc-rms_ccold)*100/rms_ccold,
     & nint(max(tmpr1,tmpr2)),
     & nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     & nint(max(abs(xav),abs(yav),abs(zav))),abs(mbad),nint(acond)
      endif
      if(isolv.eq.2.and.idata.eq.2) then
       if(iter.eq.1) write(*,'(/,"  IT   EV  CT",
     & "    RMSCT   RMSST   DX   DY   DZ   DT   OS  AQ  CND",/,
     &                           "        %   %",
     & "   ms     %    ms    m    m    m   ms    m ")')
       write(*,'(i2,a3,2(1x,i3),i5,f6.1,i6,4i5,i5,i4,i5)')
     & iter,str3,
     & nint(nev*100./nevold),
     & nint(nct*100.0/nctold),
     & nint(rms_ct*1000),(rms_ct-rms_ctold)*100/rms_ctold,
     & nint(max(tmpr1,tmpr2)),
     & nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     & nint(max(abs(xav),abs(yav),abs(zav))),abs(mbad),nint(acond)
       write(log,'(/,"  IT   EV  CT",
     & "    RMSCT   RMSST   DX   DY   DZ   DT   OS  AQ  CND",/,
     &                           "        %   %",
     & "   ms     %    ms    m    m    m   ms    m ")')
       write(log,'(i2,a3,2(1x,i3),i5,f6.1,i6,4i5,i5,i4,i5)')
     & iter,str3,
     & nint(nev*100./nevold),
     & nint(nct*100.0/nctold),
     & nint(rms_ct*1000),(rms_ct-rms_ctold)*100/rms_ctold,
     & nint(max(tmpr1,tmpr2)),
     & nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     & nint(max(abs(xav),abs(yav),abs(zav))),abs(mbad),nint(acond)
      endif

      call datetime(dattim)
      write(log,'("Iteration ",i2," finished ",a25)') iter, dattim

      if(iter.eq.maxiter) goto 600	! all iterations done.
      iter= iter+1
      goto 55	! next iteration

c--- update origin time (this is only done for final output!!)
600   continue
      write(*,'(/,"writing out results ...")')
      do i= 1,nev
         src_t(i)= src_t(i)/1000	!from here on src_t in sec!!
         if(src_t(i).gt.5) then
            write(*,*)'WARNING: org time diff > 5s for ',src_cusp(i)
         endif
         iyr= int(ev_date(i)/10000)
         imo= int(mod(ev_date(i),10000)/100)
         idy= int(mod(ev_date(i),100))
         ihr= int(ev_time(i)/1000000)
         imn= int(mod(ev_time(i),1000000)/10000)
         itf= JULIAM(iyr,imo,idy,ihr,imn)

cfw         sc= (mod(real(ev_time(i)),10000)/100) + src_t(i)
         sc= (mod(dble(ev_time(i)),10000.)/100) - src_t(i)
cfw130130         sc= (mod(real(ev_time(i)),10000.)/100) - src_t(i)
cfw100806         sc= (mod(real(ev_time(i)),10000)/100) - src_t(i)
         itf= itf + int(sc / 60.)
         sc=  sc  - int(sc / 60.)*60.
cfw130130         if(sc.lt.0) then
cfw130130            itf= itf-1
cfw130130            sc= 60. + sc
cfw130130         endif
c        ad hoc fix to deal with the output precision ...
         if(nint(sc*100).lt.0) then
            itf= itf-1
            sc= 60. + sc
         elseif(nint(sc*100).ge.6000) then
            itf= itf+1
            sc= sc - 60.
         endif
         call DATUM(itf,iyr,imo,idy,ihr,imn)
         ev_date(i)= iyr*10000 + imo*100 + idy
         ev_time(i)= ihr*1000000 + imn*10000 + nint(sc*100)
      enddo

c--- get # of obs per event:
      do i=1,nev
         src_np(i)= 0
         src_ns(i)= 0
         src_nnp(i)= 0
         src_nns(i)= 0
         src_rmsc(i)= 0
         src_rmsn(i)= 0
      enddo
      do i=1,ndt
         if(dt_idx(i).eq.1) then
             src_np(dt_ic1(i))= src_np(dt_ic1(i))+1
             src_np(dt_ic2(i))= src_np(dt_ic2(i))+1
         endif
         if(dt_idx(i).eq.2) then
             src_ns(dt_ic1(i))= src_ns(dt_ic1(i))+1
             src_ns(dt_ic2(i))= src_ns(dt_ic2(i))+1
         endif
         if(dt_idx(i).le.2) then
             src_rmsc(dt_ic1(i))= src_rmsc(dt_ic1(i))+dt_res(i)**2
             src_rmsc(dt_ic2(i))= src_rmsc(dt_ic2(i))+dt_res(i)**2
         endif
         if(dt_idx(i).eq.3) then
             src_nnp(dt_ic1(i))= src_nnp(dt_ic1(i))+1
             src_nnp(dt_ic2(i))= src_nnp(dt_ic2(i))+1
         endif
         if(dt_idx(i).eq.4) then
             src_nns(dt_ic1(i))= src_nns(dt_ic1(i))+1
             src_nns(dt_ic2(i))= src_nns(dt_ic2(i))+1
         endif
         if(dt_idx(i).ge.3) then
             src_rmsn(dt_ic1(i))= src_rmsn(dt_ic1(i))+dt_res(i)**2
             src_rmsn(dt_ic2(i))= src_rmsn(dt_ic2(i))+dt_res(i)**2
         endif
      enddo
      do i=1,nev
c100710         src_rmsc(i)= sqrt(src_rmsc(i)/nev)
c100710         src_rmsn(i)= sqrt(src_rmsn(i)/nev)
c101116         src_rmsc(i)= sqrt(src_rmsc(i)/(src_np(i)+src_ns(i)))
c101116         src_rmsn(i)= sqrt(src_rmsn(i)/(src_nnp(i)+src_nns(i)))
         if(src_np(i)+src_ns(i).gt.0) then
            src_rmsc(i)= sqrt(src_rmsc(i)/(src_np(i)+src_ns(i)))
         else
            src_rmsc(i)= -9
         endif
         if(src_nnp(i)+src_nns(i).gt.0) then
            src_rmsn(i)= sqrt(src_rmsn(i)/(src_nnp(i)+src_nns(i)))
         else
            src_rmsn(i)= -9
         endif
      enddo

c--- output final residuals: mdat.res
      if(trimlen(fn_res).gt.1) then
         call freeunit(iunit)
         open(iunit,file=fn_res,status='unknown')
c         write(iunit,'("STA",11x,"DT",8x,
c     &"C1",8x,"C2",4x,"IDX",5x,"QUAL",4x,"RES [ms]",3x,"WT",9x,
c     &"OFFS")')
         write(iunit,'(a7,1x,f12.7,1x,i9,1x,i9,1x,i1,1x,
     & f9.4,1x,f12.6,1x,f11.6,1x,f8.1)')
     & (dt_sta(j),dt_dt(j),dt_c1(j),dt_c2(j),dt_idx(j),dt_qual(j),
     & dt_res(j)*1000,dt_wt(j),dt_offs(j),j=1,ndt)
         close(iunit)
      endif

c--- output final locations (mdat.reloc):
      call setorg(sdc0_lat,sdc0_lon,0.0,0)
      do i=1,nev
         call SDC2(x,y,src_lat(i),src_lon(i),-1)
         write(fu1,'(i9,1x,f10.6,1x,f11.6,1x,f9.3,1x,f10.1,1x,f10.1,
     & 1x,f10.1,
     & 1x,f8.1,1x,f8.1,1x,f8.1,1x,i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f6.3,
     & 1x,f5.2,1x,i5,1x,i5,1x,i5,1x,i5,1x,f6.3,1x,f6.3,1x,i3)')
cfw     & 1x,f4.1,1x,i5,1x,i5,1x,i5,1x,i5,1x,f6.3,1x,f6.3,1x,i3)') 120821
cfw     & 1x,f3.1,1x,i5,1x,i5,1x,i5,1x,i5,1x,f6.3,1x,f6.3,1x,i3)') !neg mag
c     & (src_cusp(i),src_lat(i),src_lon(i),src_dep(i),src_x(i),src_y(i),
c     & src_cusp(i),src_lat(i),src_lon(i),src_dep(i),src_x(i),src_y(i),
     & src_cusp(i),src_lat(i),src_lon(i),src_dep(i),x*1000,y*1000,
     & src_z(i),src_ex(i),src_ey(i),src_ez(i),int(ev_date(i)/10000),
     & int(mod(ev_date(i),10000)/100),mod(ev_date(i),100),
     & int(ev_time(i)/1000000),int(mod(ev_time(i),1000000)/10000),
     & mod(real(ev_time(i)),10000.)/100,ev_mag(i),
cfw100806     & mod(real(ev_time(i)),10000)/100,ev_mag(i),
     & src_np(i),src_ns(i),src_nnp(i),src_nns(i),
     & src_rmsc(i),src_rmsn(i), iclust
      enddo
      call setorg(sdc0_lat,sdc0_lon,rot_3d,0)

c--- output stations (mdat.station):
      if(trimlen(fn_stares).gt.1) then
cfw         write(fu3,'(a5,1x,f9.4,1x,f9.4,1x,f9.4,1x,f9.4,1x,i7,1x,
         write(fu3,'(a7,1x,f9.4,1x,f9.4,1x,f9.4,1x,f9.4,1x,i7,1x,
     & i7,1x,i7,1x,i7,1x,f9.4,1x,f9.4,1x,i3)')
     & (sta_lab(i),sta_lat(i),sta_lon(i),sta_dist(i),sta_az(i),
     & sta_np(i),sta_ns(i),sta_nnp(i),sta_nns(i),
     & sta_rmsc(i),sta_rmsn(i),iclust,i=1,nsta)
      endif

778   continue
      enddo  ! loop over clusters (iclust)

      close(fu0)
      if(trimlen(fn_stares).gt.1) close(fu3)
      close(fu1)

      goto 999
998   write(*,*)'>>> ERROR OPENING CONTROL PARAMETER FILE'
999   continue

      call datetime(dattim)
      write(*,'(a28,a25)')'Program hypoDD finished. ',dattim
      write(log,*)
      write(log,'(a28,a25)')'Program hypoDD finished. ',dattim
      end !of main routine
