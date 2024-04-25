	subroutine partials_1dmm(log,fn_srcpar,fn_mod1d,
     &	nsrc, src_cusp, src_lat, src_lon, src_dep,
     &	nsta, sta_lab, sta_lat, sta_lon, sta_elv, sta_mod,
     &	mod_nl, mod_ratio, mod_v, mod_top,iter,
     &	tmp_ttp, tmp_tts,
     &	tmp_xp, tmp_yp, tmp_zp,tmp_xs, tmp_ys, tmp_zs)

	implicit none

	include'hypoDD.inc'

c	Parameters:
        character       fn_mod1d*80     ! Model file
	character	fn_srcpar*110	! Source-parameter file
        integer         log             ! Log-file identifier
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
        integer         sta_mod(MAXSTA)
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
        integer         iter           ! iteration number


c	Local variables:
	real		ainp
	real		ains
	real		del
	real 		az1
	real		dist
	integer		i, j, k
	integer		iunit		! Output unit number
	real		pi
	integer		trimlen
        real            vs(MAXLAY)

        integer         ai(100)
        integer         an(100)
        real            az(100,MAXLAY)
        real            avp(100,MAXLAY)
        real            avs(100,MAXLAY)
        integer         nmod           ! Number of 1D models


	parameter(pi=3.141593)

c     Open station specific model file:
c     (models are read in at the beginning of every iteration. move
c      this to the getinp routine to do it just once.)
      if(iter.ge.1) then
         call freeunit(iunit)
         open(iunit,file=fn_mod1d,status='unknown')
         i= 1
5        read(iunit,*,end=7,err=6)ai(i),an(i)
         if(an(i).gt.MAXLAY) stop"Increase MAXLAY array in hypoDD.inc."
         do j=1,an(i)
            read(iunit,*)az(i,j),avp(i,j),avs(i,j)
         enddo
         i= i+1
         if(i.gt.100) stop'Increase model array size in hypoDD.f'
         goto 5
6        stop'Format error in model file.'
7        nmod= i-1
         if(iter.eq.1) then
            write(*,'(a,i4)')'Station spec 1D models read in = ',nmod
            write(log,'(a,i4)')'Station spec 1D models read in = ',nmod
         endif
         close(iunit)
      endif


      iunit = 0
      if (trimlen(fn_srcpar).gt.1) then
c        Open source-parameter file
         call freeunit(iunit)
         open(iunit,file=fn_srcpar,status='unknown')
      endif

cc     Make sure hypocenters don't fall on layer boundaries
c      do i=1,nsrc
c         do j=1,mod_nl
c            if (abs(src_dep(i)-mod_top(j)).lt.0.0001)
c     &         src_dep(i) = src_dep(i)-0.001
c         enddo
c      enddo

cc     Get S velocity model
c      do i=1,mod_nl
c         vs(i) = mod_v(i)/mod_ratio(i)
c      enddo

c     Compute epicentral distances, azimuths, angles of incidence,
c     and P/S-travel times from sources to stations
      do i=1,nsta

c         Find model
8         do j= 1,nmod
             if(sta_mod(i).eq.ai(j)) then
                mod_nl= an(j)
                do k= 1,mod_nl
                   mod_top(k)= az(j,k)
                   mod_v(k)= avp(j,k)
                   vs(k)= avs(j,k)
c                   write(*,*)mod_top(k),mod_v(k),vs(k),mod_nl
                enddo
                goto 9
             endif
          enddo
          write(log,*)'No model found for model ID',sta_mod(i),
     &  '. Default model (#1) chosen.'
          sta_mod(i)= 1
          goto 8
9         continue

c110505fw: adjust for station elevation
         mod_top_adj(1)= 0.0
         do k=2,mod_nl
            mod_top_adj(k)= mod_top(k)+sta_elv(i)
         enddo

         do j=1,nsrc
c           Make sure hypocenters don't fall on layer boundaries
            do k=1,mod_nl
               if (abs(src_dep(j)-mod_top(k)).lt.0.0001)
     &            src_dep(j) = src_dep(j)-0.001
            enddo

            call delaz2(src_lat(j), src_lon(j), sta_lat(i), sta_lon(i), 
     &                 del, dist, az1)

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
     &               cos(((az1 - 90) * pi)/180.0))/mod_v(k-1)
	    tmp_xs(i,j) = (sin((ains * pi)/180.0) *
     &               cos(((az1 - 90) * pi)/180.0))/vs(k-1)
	    tmp_yp(i,j) = (sin((ainp * pi)/180.0) *
     &               cos((az1 * pi)/180.0))/mod_v(k-1)
	    tmp_ys(i,j) = (sin((ains * pi)/180.0) *
     &               cos((az1 * pi)/180.0))/vs(k-1)

c           Write to source-parameter file
            if (iunit .ne. 0)
     &         write(iunit,'(i9,2x,f9.4,2x,f9.4,2x,a7,f6.3,2x,
     &         f9.4,2x,f9.4,2x,f9.4,f9.4,f9.3,f9.3)')
     &         src_cusp(j), src_lat(j), src_lon(j), sta_lab(i), 
     &         sta_elv(i),dist,az1,ainp,ains,tmp_ttp(i,j), tmp_tts(i,j)

         enddo
      enddo

      if (iunit .ne. 0) close(iunit)	! Source-parameter file

      end !of subroutine partials_1dmm
