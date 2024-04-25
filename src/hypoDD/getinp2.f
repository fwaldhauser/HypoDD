	subroutine getinp2 (maxev,maxlyr,log,fn_inp,
     &	fn_cc, fn_ct, fn_sta, fn_eve,
     &	fn_loc, fn_reloc, fn_res, fn_stares, fn_srcpar,
     &  iwrite, 	
     &	idata, iphase,
     &	minobs_cc, minobs_ct,
     &	amaxres_cross, amaxres_net, amaxdcc, amaxdct,
     &	noisef_dt, maxdist,minds,maxds,maxgap,
     &	awt_ccp, awt_ccs, awt_ctp, awt_cts, adamp,
     &	istart, maxiter, isolv, iaq, niter, aiter,
     &	imod, mod_nl, mod_ratio, mod_v, mod_top,
     &  fn_mod1d,fn_mod3d, lat_3d, lon_3d, rot_3d, 
     &  ipha3d,ndip,iskip,scale1,scale2,xfac,tlim,nitpb, 
     &	iclust, ncusp, icusp)

	implicit none

	include'hypoDD.inc'

c	Parameters:
	integer		maxev		! Array dimension
	integer		maxlyr		! Array dimension
	integer		log		! Log-file identifier
	character	fn_inp*110	! File of control info.
	character	fn_cc*110	! File of cross-corr. times
	character	fn_ct*110	! File of catalog times
	character	fn_sta*110	! Station file
	character	fn_eve*110	! Event file
	character	fn_loc*110	! Output file of original locs.
	character	fn_reloc*110	! Output file of final locs.
	character	fn_res*110	! Output residual file
	character	fn_stares*110	! Output station file
	character	fn_srcpar*110	! Output source-parameter file
cx	character	fn_tbl*110	! Output source-parameter file
	integer		idata		! 0: Synthetics
					! 1: Cross-correlation
					! 2: Catalog
					! 3: Both
	integer		iphase		! 1: P; 2: S; 3: Both
	integer		minobs_cc	! Min. obs./pair for ccor. data
	integer		minobs_ct	! Min. obs./pair for cat. data
	real		amaxres_cross(10)! [1..niter] Ccor. res. thresh.
	real		amaxres_net(10)	! [1..niter] Cat. res. thresh.
	real		amaxdcc(10)	! [1..niter] Ccor. link-dist. limit
	real		amaxdct(10)	! [1..niter] Cat. link-dist. limit
	real		noisef_dt	! Synthetic noise
	real		maxdist		! Max. cluster-station distance
	real		minds		! Max. pair-station distance
	real		maxds		! Max. pair-station distance
	real		maxgap		! Max. azimuthal gap 
	real		awt_ccp(10)	! [1..niter] Wts. for ccor. P
	real		awt_ccs(10)	! [1..niter] Wts. for ccor. S
	real		awt_ctp(10)	! [1..niter] Wts. for cat. P
	real		awt_cts(10)	! [1..niter] Wts. for cat. S
	real		adamp(10)	! [1..niter] Damping (lsqr only)
	integer		istart		! 1: From single source
					! 2: From network sources
	integer		maxiter
	integer		isolv		! 1: SVD; 2: LSQR
	integer		niter		! No. of iteration sets
	integer		aiter(0:10)	! [1..niter] Iterations/set
cx	integer		bsiter		! # of bootstrap iterations
	integer		mod_nl		! No. of layers
	real		mod_ratio(MAXLAY)	! Vp/Vs
	real		mod_v(MAXLAY)	! [1..mod_nl] Vp values
	real		mod_top(MAXLAY)	! [1..mod_nl] Depths to layers
	integer		iclust		! Cluster to relocate (0: all).
	integer		ncusp		! No. of event keys in icusp[]
	integer		icusp(MAXEVE)	! [1..ncusp] Events to relocate
	integer		imod
	integer		iv
cx     	real		pwin
cx	real		swin	
        character       fn_mod3d*110
        character       fn_mod1d*110
        real	   	lat_3d
	real 		lon_3d
	real 		rot_3d 
     	integer		ipha3d
     	integer		ndip
     	integer		iskip
	real 		scale1
	real 		scale2
	real 		xfac
	real 		tlim
	integer		nitpb
	integer		iwrite	
	integer		iaq

c	Local variables:
	real 		a
	integer		fu_inp
	integer		i
	integer		ii
	integer		l
	character	line*220
	integer		trimlen
	integer		j
        logical         ex

c--- newest format: 083000 with iteration step dependent weighting

c-- set intial values:
      lat_3d= -999
      lon_3d= -999
      rot_3d= 0.0

c-- open input file:
      call freeunit (fu_inp)
      open (fu_inp,status='unknown',file=fn_inp,err=998)
      write (log,'("Input parameters: ",a," (hypoDD2.0 format)")')
     & fn_inp(1:trimlen(fn_inp))
      write (*,'("Input parameters: ",a," (hypoDD2.0 format)")')
     & fn_inp(1:trimlen(fn_inp))
      ncusp= 0
      niter= 0  ! number of iteration blocks
cx      bsiter= 9999  ! iteration at which to start bootstrapping 
      l = 1
      ii= 1
      iv= 1

c-- Loop to read each parameter lines, skipping comments
      read (fu_inp,*) 	!skip format identifier line
 210  read (fu_inp,'(a)',end=300) line
      if (line(1:1).eq.'*' .or. line(2:2).eq.'*') goto 210
      if (l.eq.1) read (line,'(a)',err=999) fn_cc
      if (l.eq.2) read (line,'(a)',err=999) fn_ct
      if (l.eq.3) read (line,'(a)',err=999) fn_eve
      if (l.eq.4) read (line,'(a)',err=999) fn_sta
      if (l.eq.5) read (line,'(a)',err=999) fn_loc
      if (l.eq.6) then
         read (line,'(a)',err=999) fn_reloc
         if(fn_reloc(trimlen(fn_reloc):trimlen(fn_reloc)).eq."-") then
            fn_reloc= fn_reloc(1:trimlen(fn_reloc)-1)
            iwrite= 0
         else
           iwrite= 1
         endif
      endif
      if (l.eq.7) read (line,'(a)',err=999) fn_stares
      if (l.eq.8) read (line,'(a)',err=999) fn_res
      if (l.eq.9) read (line,'(a)',err=999) fn_srcpar
      if (l.eq.10) read (line,*,err=999) idata, iphase,maxdist
      if(maxdist.lt.0) maxdist= 50000
      if (l.eq.11) read (line,*,err=999) 
c     &    minobs_cc,minobs_ct
     &    minobs_cc,minobs_ct,minds,maxds,maxgap
c      maxgap= -9	! not implemented 
c      minds= -9		! not implemented 
c      maxds= -9		! not implemented 
      if (l.eq.12) read (line,*,err=999) istart, isolv, iaq, niter
      if(niter.gt.10) stop'Maximum for NSET is 10.'

c--Read iteration instructions
      if (l.ge.13 .and. l.le.12+niter) then
         i=l-12
         read (line,*,err=999) aiter(i),
     &  awt_ccp(i), awt_ccs(i), amaxres_cross(i), amaxdcc(i),
     &  awt_ctp(i), awt_cts(i), amaxres_net(i), amaxdct(i), adamp(i)
        if(amaxres_cross(i).eq.-9) amaxres_cross(i)= -999
        if(amaxres_net(i).eq.-9) amaxres_net(i)= -999
        if(amaxdcc(i).eq.-9) amaxdcc(i)= -999
        if(amaxdct(i).eq.-9) amaxdct(i)= -999
        if(awt_ccp(i).eq.-9) awt_ccp(i)= -999
        if(awt_ccs(i).eq.-9) awt_ccs(i)= -999
        if(awt_ctp(i).eq.-9) awt_ctp(i)= -999
        if(awt_cts(i).eq.-9) awt_cts(i)= -999
	if(isolv.eq.1) adamp(i)= abs(adamp(i))		!070911/fw
      endif

c--- Read models:
      if(l.eq.13+niter) read(line,*,err=999) imod

c      if(l.eq.14+niter.and.(imod.eq.1.or.imod.eq.4.or.imod.eq.6)) 
c     &     	read(line,*,err=999) mod_nl, mod_ratio	
c      if(l.eq.15+niter.and.(imod.eq.1.or.imod.eq.4.or.imod.eq.6)) 
c     &		read(line,*,err=999) (mod_top(i),i=1,mod_nl)
c      if(l.eq.16+niter.and.(imod.eq.1.or.imod.eq.4.or.imod.eq.6))
c     & 		read(line,*,err=999) (mod_v(i),i=1,mod_nl)
c      if(l.eq.17+niter.and.(imod.eq.2.or.imod.eq.3)) 
c     &  	read(line,'(a)',err=999) fn_tbl 
c      if(l.eq.18+niter.and.(imod.eq.2.or.imod.eq.3)) 
c     &  	read(line,*,err=999) pwin, swin 

c read 1D model in old format (fixed vp/vs ratio:
      if(l.eq.14+niter.and.(imod.eq.0.or.imod.eq.6)) then
         read(line,*,err=999) mod_nl,a
         do i=1,mod_nl
	    mod_ratio(i)= a
         enddo
      endif
      if(l.eq.15+niter.and.(imod.eq.0.or.imod.eq.6)) 
     &		read(line,*,err=999) (mod_top(i),i=1,mod_nl)
      if(l.eq.16+niter.and.(imod.eq.0.or.imod.eq.6))
     & 	 	read(line,*,err=999) (mod_v(i),i=1,mod_nl)
      
c read 1D model in new format (variable vp/vs ratio):
      if(l.eq.14+niter.and.(imod.eq.1.or.imod.eq.5)) then
         read(line,*,err=999) 
         call nval(line,mod_nl)
         if(imod.eq.5) mod_nl=1
         read(line,*,err=999) (mod_top(i),i=1,mod_nl)
         if(mod_top(mod_nl).eq.-9) mod_nl= mod_nl-1  ! catch potential user error

c110817         do kk=1,1000
c110817      	    read(line,*,err=220) (mod_top(i),i=1,kk)
c110817	    if(mod_top(kk).lt.-8) goto 220
c110817         enddo
c110817220      mod_nl= kk-1

         if(mod_nl.gt.MAXLAY) then
	    write(*,*)
            write(*,'(a)')
     & '>>> Number of 1D model layers exceeding MAXLAY.'
            write(log,'(a)')
     & '>>> Number of 1D model layers exceeding MAXLAY.'
            stop'Program run aborted.'
         endif
      endif
      if(l.eq.15+niter.and.(imod.eq.1.or.imod.eq.5)) 
     &		read(line,*,err=999) (mod_v(i),i=1,mod_nl)
      if(l.eq.16+niter.and.(imod.eq.1.or.imod.eq.5))
     & 		read(line,*,err=999) (mod_ratio(i),i=1,mod_nl)

c      if(l.eq.14+niter.and.(imod.eq.2.or.imod.eq.3)) 
c     &  	read(line,'(a)',err=999) fn_tbl 
c      if(l.eq.15+niter.and.(imod.eq.2.or.imod.eq.3)) 
c     &  	read(line,*,err=999) pwin, swin 
c      if(l.eq.16+niter.and.(imod.eq.2.or.imod.eq.3)) 
c     &  	read(line,*,err=999) 
      if(l.eq.14+niter.and.imod.eq.9) 
     &  	read(line,'(a)',err=999) fn_mod3d
      if(l.eq.15+niter.and.imod.eq.9) 
     &  	read(line,*,err=999) lat_3d,lon_3d,rot_3d 
      if(l.eq.16+niter.and.imod.eq.9) 
     &  	read(line,*,err=999) 
     &          ipha3d,ndip,iskip,scale1,scale2,xfac,tlim,nitpb 
      if(l.eq.14+niter.and.imod.eq.4) 
     &  	read(line,'(a)',err=999) fn_mod1d
c      if(l.eq.15+niter.and.imod.eq.4) 
c     &  	read(line,*,err=999)
c      if(l.eq.16+niter.and.imod.eq.4) 
c     &  	read(line,*,err=999)

c--Read specific clusters/events to relocate
      if(l.eq.17+niter) read (line,*,err=999) iclust
      if(l.ge.18+niter) then
c100823         	read (line,*,err=999,end=230) (icusp(i),i=ii,ii+7)
c100823         	if(i.gt.maxeve)stop'Increase MAXEVE array in hypoDD.inc'
c100823230      ii= i
                call nval(line,j)
                read (line,*,err=999) (icusp(i),i=ii,ii+ j-1)
                ii= i

c                read (line,*,err=999,end=229) (icusp(i),i=ii,ii+7)
c229             do j=ii,ii+7
c                   if(icusp(j).eq.0) goto 230
c                enddo
                if(i.gt.maxeve)stop'Increase MAXEVE array in hypoDD.inc'
c230      ii= j
      endif

250   l= l+1
      goto 210
300   close (fu_inp)
      ncusp= ii-1

c- rearrange aiter:
      do i=2,niter
cx        if(aiter(i).lt.0) bsiter= aiter(i-1)+1    ! boostrap iterations
        aiter(i)= aiter(i-1)+abs(aiter(i))
      enddo

c- check files
      call exist (fn_eve)
      call exist (fn_sta)
      if ((idata.eq.1 .or.idata.eq.3).and.trimlen(fn_cc).gt.1)
     & call exist(fn_cc)
      if ((idata.eq.2 .or.idata.eq.3).and.trimlen(fn_ct).gt.1)
     & call exist (fn_ct)

      maxiter= aiter(niter)
c synthetic noise:
      noisef_dt= 0.002

c write log output: of newest format
600   if (trimlen(fn_loc).lt.2) fn_loc= 'hypoDD.loc'
      if (trimlen(fn_reloc).lt.2) fn_reloc= 'hypoDD.reloc'
      write (6,'("INPUT FILES:",/,
     &" cross dtime data: ",a,/," catalog dtime data: ",a,/,
     &" events: ",a,/," stations: ",a)')
     &fn_cc(1:trimlen(fn_cc)),
     &fn_ct(1:trimlen(fn_ct)),
     &fn_eve(1:trimlen(fn_eve)),
     &fn_sta(1:trimlen(fn_sta))
cx      if(imod.eq.3) write(6,'(" travel time table: ",a)')
cx     &fn_tbl(1:trimlen(fn_tbl))
      write (6,'("OUTPUT FILES:",/,
     &" initial locations: ",a,/," relocated events: ",a,/,
     &" event pair residuals: ",a,/," station residuals: ",a,/,
     &" source parameters: ",a)')
     &fn_loc(1:trimlen(fn_loc)),
     &fn_reloc(1:trimlen(fn_reloc)),fn_res(1:trimlen(fn_res)),
     &fn_stares(1:trimlen(fn_stares)),fn_srcpar(1:trimlen(fn_srcpar))

      write (log,'("INPUT FILES:",/,
     &" cross dtime data: ",a,/," catalog dtime data: ",a,/,
     &" events: ",a,/," stations: ",a)')
     &fn_cc(1:trimlen(fn_cc)),
     &fn_ct(1:trimlen(fn_ct)),
     &fn_eve(1:trimlen(fn_eve)),
     &fn_sta(1:trimlen(fn_sta))
cx      if(imod.eq.3) write(log,'(" travel time table: ",a)')
cx     &fn_tbl(1:trimlen(fn_tbl))
      write (log,'("OUTPUT FILES:",/,
     &" initial locations: ",a,/," relocated events: ",a,/,
     &" event pair residuals: ",a,/," station residuals: ",a,/,
     &" source parameters: ",a)')
     &fn_loc(1:trimlen(fn_loc)),
     &fn_reloc(1:trimlen(fn_reloc)),fn_res(1:trimlen(fn_res)),
     &fn_stares(1:trimlen(fn_stares)),fn_srcpar(1:trimlen(fn_srcpar))

      write (log,'(/,
     &"  IDATA= ",i2,2X,"IPHASE= ",i2,2x,
     &"MAXDIST (from cluster)= ",f9.0,/,
     &"  MINOBS_CC= ",i3,2x,"MINOBS_CT= ",i3,/,
     &"  MINDIST (from pair)=",f7.0,2x,
     &"MAXDIST (from pair)= ",f9.0,2x,
     &"MAXGAP= ",f6.0,/,
     &"  ISTART= ",i1,2x,
     &"ISOLV= ",i1,2x,
     &"IAQ= ",i1,2x)')
     &idata,iphase,maxdist,minobs_cc,minobs_ct,
     &minds,maxds,maxgap,istart,isolv,iaq

      aiter(0)=0
      write (log, '("  ITER ",i3,"-",i3,
     &": DAMP= "f7.1,/,"    WT_CCP= ",f7.2,2X,"WT_CCS= ",f7.2,2x,
     &"MAXR_CC= ",f7.2,2X,"MAXD_CC= ",f7.2,2X,/,
     &"    WT_CTP= ",f7.2,2x,"WT_CTS= ",f7.2,2x,"MAXR_CT= ",f7.2,2x,
     &"MAXD_CT= ",f7.2)')
     &(aiter(i-1)+1,aiter(i),adamp(i),awt_ccp(i),awt_ccs(i),
     & amaxres_cross(i),
     & amaxdcc(i),awt_ctp(i),awt_cts(i), amaxres_net(i), amaxdct(i),
     & i=1,niter)
cx      if(bsiter.ne.9999) 
cx     &write (log, '("  Bootstrap iterations start at iteration ",i3)'),
cx     & bsiter

c--Write crustal model
      if(imod.eq.0.or.imod.eq.1) then
        write(*,'(a)')'Use local layered 1D model.'
        write(log,'(a)')'Use local layered 1D model.'
        write (log, '("  MOD_NL= ",i2)')
     &    mod_nl
        write(log,'("    MOD_TOP    MOD_VP   MOD_VS   MOD_RATIO")')
        do i=1,mod_nl
             write(log,'(2x,4f9.3)')
     &    mod_top(i),mod_v(i),mod_v(i)/mod_ratio(i),mod_ratio(i)
        enddo
      elseif(imod.eq.5) then
        write(*,'(a)')'Use constant velocity model.'
        write(log,'(a)')
     &  'Use constant velocity model (only first layer used).'
        write (log, '("  MOD_NL= ",i2)')
     &    mod_nl
        write(log,'("    MOD_TOP    MOD_VP   MOD_VS   MOD_RATIO")')
        do i=1,mod_nl
             write(log,'(2x,4f9.3)')
     &    mod_top(i),mod_v(i),mod_v(i)/mod_ratio(i),mod_ratio(i)
        enddo
      elseif(imod.eq.4) then
        write(*,'(2a)')'Use station specific 1D models from file ',
     &  fn_mod1d(1:trimlen(fn_mod1d)) 
        write(log,'(2a)')'Use station specific 1D models from file ',
     &  fn_mod1d(1:trimlen(fn_mod1d)) 
        inquire(FILE= fn_mod1d,exist=ex)
        if(.not. ex) stop
     & ' >>> ERROR OPENING STATION SPECIFIC MODELS FILE.'
      elseif(imod.eq.9) then
        write(*,'(a,a20)')'Use local 3D model ',
     &    fn_mod3D
        write(log,'(a,a20)')'Use local 3D model ',
     &    fn_mod3D
      endif

c--Repeat number of clusters, events to relocate
      if (iclust.eq.0) then
        write (*,'(a)') 'Relocate all clusters'
        write (log,'(a)') 'Relocate all clusters'
      else
        write (*,'(a,i4)') 'Relocate cluster number ',iclust
        write (log,'(a,i4)') 'Relocate cluster number ',iclust
      endif

      if (ncusp.eq.0) then
        write (*,'(a)') 'Relocate all events'
        write (log,'(a)') 'Relocate all events'
      else
        write (*,'(a,i6,a)') 'Relocate ',ncusp,' events'
        write (log,'(a,i6,a)') 'Relocate ',ncusp,' events'
      endif

      if (iaq.eq.0) then
        write (*,'(a)') 'Keep air quakes.'
        write (log,'(a)') 'Keep air quakes.'
      else
        write (*,'(a)') 'Remove air quakes.'
        write (log,'(a)') 'Remove air quakes.'
      endif

      return

c--Input error handling
998   write(*,*)'>>> ERROR OPENING CONTROL PARAMETER FILE'
      goto 1000

999   write (*,*)'>>> ERROR READING CONTROL PARAMETERS IN LINE ',l
      write (*,*) line
1000  stop 'Program run aborted.'
      end  ! of subroutine getinp2

