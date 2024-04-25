	subroutine partials_3d(fn_srcpar,
     &	nsrc, src_cusp, src_lat, src_lon, src_dep,
     &	nsta, sta_lab, sta_lat, sta_lon,sta_elv,
     &  mod_nl,mod_ratio,mod_v,mod_top,rot_3d,ipha3d,
     &	tmp_ttp, tmp_tts,
     &	tmp_xp, tmp_yp, tmp_zp,
     &	tmp_xs, tmp_ys, tmp_zs)

	implicit none

	include'hypoDD.inc'
	include'vel3d.inc'

c	Parameters:
	character	fn_srcpar*110	! Source-parameter file
	integer		nsrc		! No. of sources
	integer		src_cusp(MAXEVE)! [1..nsrc]
	doubleprecision	src_lat(MAXEVE)	! [1..nsrc]
	doubleprecision	src_lon(MAXEVE)	! [1..nsrc]
	real		src_dep(MAXEVE)	! [1..nsrc]
	integer		nsta		! No. of stations
	character	sta_lab(MAXSTA)*7! [1..nsta]
	real		sta_lat(MAXSTA)	! [1..nsta]
	real		sta_lon(MAXSTA)	! [1..nsta]
        real            sta_elv(MAXSTA) ! [1..MAXSTA]
	integer		mod_nl		! No. of layers
c	real		mod_ratio	! Vp/Vs
	real		mod_ratio(MAXLAY)	! Vp/Vs
	real		mod_v(MAXLAY)	! [1..mod_nl]
	real		mod_top(MAXLAY)	! [1..mod_nl]
	real		tmp_ttp(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]
	real		tmp_tts(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]
	real		tmp_xp(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]
	real		tmp_yp(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]
	real		tmp_zp(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]
	real		tmp_xs(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]
	real		tmp_ys(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]
	real		tmp_zs(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]
        real 		rot_3d
        integer 	ipha3d	

c	Local variables:
        real            az
	real		ainp, ains
	real		azp,azs
	real		del
	real		dist
	integer		i, j, k
	integer		iunit		! Output unit number
	real		pi
	integer		trimlen
        real		sx,sy,rx,ry
	real 		ttimep, ttimes
        real		velp, vels
        integer		isp
        real            vs(MAXLAY)

	parameter(pi=3.141593)

      iunit = 0
      if (trimlen(fn_srcpar).gt.1) then
c        Open source-parameter file
         call freeunit(iunit)
         open(iunit,file=fn_srcpar,status='unknown')
      endif

c     Compute epicentral distances, azimuths, angles of incidence,
c     and P/S-travel times from sources to stations

      do i=1,nsta
         do j=1,nsrc

             if(tmp_ttp(i,j).gt.-990) then
c            goto 666

c	    3D ray tracing  for P waves
	    call SDC2(sx,sy,src_lat(j), src_lon(j),-1)
	    call SDC2(rx,ry,dble(sta_lat(i)), dble(sta_lon(i)),-1)

	    isp = 0
c	    write(*,*)isp,sx,sy,src_dep(j),rx,ry,0.0,ttimep,azp,
c     &               ainp,velp
	    call path(isp,sx,sy,src_dep(j),rx,ry,-sta_elv(i),
     &                ttimep,azp,ainp,velp)
	    tmp_ttp(i, j) = ttimep
c           Depth derivative
            tmp_zp(i,j) = cos((ainp * pi)/180.0)/velp
c           Epicentral derivatives
     	    tmp_xp(i,j) = (sin((ainp * pi)/180.0) *
     &               cos(((azp - 90) * pi)/180.0))/velp
     	    tmp_yp(i,j) = (sin((ainp * pi)/180.0) *
     &               cos((azp * pi)/180.0))/velp

            if(ipha3d.eq.2) then
  	       isp = 1
	       call path(isp,sx,sy,src_dep(j),rx,ry,-sta_elv(i),
     &                   ttimes,azs,ains,vels)
	       tmp_tts(i, j) = ttimes
c              Depth derivative
               tmp_zs(i,j) = cos((ains * pi)/180.0)/vels
c              Epicentral derivatives
     	       tmp_xs(i,j) = (sin((ains * pi)/180.0) *
     &               cos(((azs - 90) * pi)/180.0))/vels
     	       tmp_ys(i,j) = (sin((ains * pi)/180.0) *
     &               cos((azs * pi)/180.0))/vels
            else
	       tmp_tts(i, j) = ttimep*1.73
c              Depth derivative
               tmp_zs(i,j) = tmp_zp(i,j) *1.73
c              Epicentral derivatives
               tmp_xs(i,j) = tmp_xp(i,j) *1.73
               tmp_ys(i,j) = tmp_yp(i,j) *1.73
            endif

c	    write(*,*)"3D: ", ttimep,ttimes,azp,ainp,azs,ains,velp,vels
c
            if(ttimep.lt.0.001.or.ttimes.lt.0.001) then 
               tmp_ttp(i,j)= -999
  	       tmp_tts(i,j)= -999 
            endif

c	    write(123,*)"3D: ", tmp_ttp(i,j),tmp_tts(i,j),
c     &     azp,ainp,azs,ains,velp,vels
c            write(123,*)src_lat(j),src_lon(j),sta_lat(i),sta_lon(i)


c           Write to source-parameter file
            dist= -999
            if (iunit .ne. 0)
     &         write(iunit,'(i9,2x,f9.4,2x,f9.4,2x,a7,2x,f9.4,
     &         2x,f9.4,2x,f9.4,2x,f9.4,2x,f9.4,f8.3,f8.3)')
     &         src_cusp(j), src_lat(j), src_lon(j), sta_lab(i), 
     &         azp, ainp, azs, ains, tmp_ttp(i,j),tmp_tts(i,j)

            endif !tmp_ttp>-999
         enddo
      enddo

      if (iunit .ne. 0) close(iunit)	! Source-parameter file

      end !of subroutine partials
