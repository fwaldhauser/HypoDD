	subroutine partials(fn_srcpar,
     &	nsrc, src_cusp, src_lat, src_lon, src_dep,
     &	nsta, sta_lab, sta_lat, sta_lon, sta_elv,
     &	mod_nl, mod_ratio, mod_v, mod_top,
     &	tmp_ttp, tmp_tts,
     &	tmp_xp, tmp_yp, tmp_zp,tmp_xs, tmp_ys, tmp_zs)

	implicit none

	include'hypoDD.inc'

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
	real		sta_elv(MAXSTA)	! [1..nsta]
	integer		mod_nl		! No. of layers
c	real		mod_ratio	! Vp/Vs
	real		mod_ratio(MAXLAY)	! Vp/Vs
	real		mod_v(MAXLAY)	! [1..mod_nl]
	real		mod_top(MAXLAY)	! [1..mod_nl]
	real		mod_top_adj(MAXLAY)	! [1..mod_nl]
	real		tmp_ttp(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]
	real		tmp_tts(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]
	real		tmp_xp(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]
	real		tmp_yp(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]
	real		tmp_zp(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]
	real		tmp_xs(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]
	real		tmp_ys(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]
	real		tmp_zs(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]

c	Local variables:
	real		ainp
	real		ains
	real		az
	real		del
	real		dist
	integer		i, j, k
	integer		iunit		! Output unit number
	real		pi
	integer		trimlen
        real            vs(MAXLAY)

	parameter(pi=3.141593)

      iunit = 0
      if (trimlen(fn_srcpar).gt.1) then
c        Open source-parameter file
         call freeunit(iunit)
         open(iunit,file=fn_srcpar,status='unknown')
      endif

c     Make sure hypocenters don't fall on layer boundaries
      do i=1,nsrc
         do j=1,mod_nl
            if (abs(src_dep(i)-mod_top(j)).lt.0.0001)
     &         src_dep(i) = src_dep(i)-0.001
         enddo
      enddo

c     Get S velocity model
      do i=1,mod_nl
         vs(i) = mod_v(i)/mod_ratio(i)
      enddo

c     Compute epicentral distances, azimuths, angles of incidence,
c     and P/S-travel times from sources to stations
c110505fw: adjust for station elevation
      do i=1,nsta
         mod_top_adj(1)= 0.0
         do k=2,mod_nl
            mod_top_adj(k)= mod_top(k)+sta_elv(i)
         enddo
         do j=1,nsrc
            call delaz2(src_lat(j), src_lon(j), sta_lat(i), sta_lon(i), 
     &                 del, dist, az)

c           1D ray tracing
            call ttime(dist, src_dep(j)+sta_elv(i), 
     &                 mod_nl, mod_v, mod_top_adj, 
     &                 tmp_ttp(i, j), ainp)
            call ttime(dist, src_dep(j)+sta_elv(i), 
     &                 mod_nl, vs, mod_top_adj, 
     &                 tmp_tts(i, j), ains)

c            call ttime(dist, src_dep(j), mod_nl, mod_v, mod_top, 
c     &                 tmp_ttp(i, j), ain)
c            call ttime(dist, src_dep(j), mod_nl, vs, mod_top, 
c     &                 tmp_tts(i, j), ain)
            
c           Determine wave speed at the hypocenter
            do k=1,mod_nl
               if (src_dep(j).le.mod_top(k)) goto 10	! break
            enddo
10          continue

c           Depth derivative
            tmp_zp(i,j) = cos((ainp * pi)/180.0)/mod_v(k-1)
            tmp_zs(i,j) = cos((ains * pi)/180.0)/vs(k-1)
c           Epicentral derivatives
	    tmp_xp(i,j) = (sin((ainp * pi)/180.0) *
     &               cos(((az - 90) * pi)/180.0))/mod_v(k-1)
	    tmp_xs(i,j) = (sin((ains * pi)/180.0) *
     &               cos(((az - 90) * pi)/180.0))/vs(k-1)
	    tmp_yp(i,j) = (sin((ainp * pi)/180.0) *
     &               cos((az * pi)/180.0))/mod_v(k-1)
	    tmp_ys(i,j) = (sin((ains * pi)/180.0) *
     &               cos((az * pi)/180.0))/vs(k-1)

c           Write to source-parameter file
            if (iunit .ne. 0) then
c               write(iunit,'(i9,2x,f9.4,2x,f9.4,2x,a7,f6.3,2x,
c     &         f9.4,2x,f9.4,2x,f9.4,f9.4,f9.3,f9.3)')
c     &         src_cusp(j), src_lat(j), src_lon(j), sta_lab(i), 
c     &         sta_elv(i),dist,az,ainp,ains,tmp_ttp(i,j), tmp_tts(i,j)

               write(iunit,'(i9,1x,f9.4,1x,f9.4,1x,a7,1x,f7.3,
     &         1x,f9.4,
     &         1x,f9.4,1x,f9.4,f9.4,f9.4,f9.4,f9.4,f9.4,f9.4,f9.4,
     &         f9.4,f9.4)')
     &         src_cusp(j), src_lat(j), src_lon(j), sta_lab(i),
     &         sta_elv(i),
     &         dist, az, ainp,ains,
     &         tmp_ttp(i,j),tmp_tts(i,j),
     &         tmp_xp(i,j),tmp_yp(i,j), tmp_zp(i,j),
     &         tmp_xs(i,j),tmp_ys(i,j), tmp_zs(i,j)

c write out to snthetic phase file:
c                write(iunit,'(a,i5,a,f12.6,f12.6,f8.4,a,i9)')
c     &          "# ",1000+j," 1 1 1 1 0.0 ",
c     &          src_lat(j),src_lon(j),
c     &          src_dep(j)," 2.0 0 0 0 ",src_cusp(j)
c                write(iunit,'(a7,f9.3,a)')
c     &          sta_lab(i),tmp_ttp(i,j)," 1.0 P"
c                write(iunit,'(a7,f9.3,a)')
c     &          sta_lab(i),tmp_tts(i,j)," 1.0 S"
            endif

         enddo
      enddo

      if (iunit .ne. 0) close(iunit)	! Source-parameter file

      end !of subroutine partials
