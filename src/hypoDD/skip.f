	subroutine skip(log, kiter, minwght,
     &	ndt, nev, nsrc, nsta,
     &	ev_cusp, ev_date, ev_time, ev_mag,
     &	ev_lat, ev_lon, ev_dep, ev_x, ev_y, ev_z,
     &	ev_herr, ev_zerr, ev_res, ev_fix,
     &	src_cusp, src_lat, src_lon, src_dep,
     &	src_lat0, src_lon0,
     &	src_x, src_y, src_z, src_t, src_x0, src_y0, src_z0, src_t0,
     &	sta_lab, sta_lat, sta_lon, sta_elv, sta_mod, sta_dist, sta_az,
     &	sta_rmsc, sta_rmsn, sta_np, sta_ns, sta_nnp, sta_nns,
     &	dt_sta, dt_c1, dt_c2, dt_idx, dt_dt, dt_qual, dt_cal,
     &	dt_ista, dt_ic1, dt_ic2,
     &	dt_res, dt_wt, dt_offs,
     &  tmp_ttp,tmp_tts,tmp_xp,tmp_yp,tmp_zp,tmp_xs,tmp_ys,tmp_zs,
     &  nct,ncc)


	implicit none

	include'hypoDD.inc'

c	Parameters:
	integer		log
	integer		kiter
	real		minwght
	integer		ndt
	integer		nev
	integer		nsrc
	integer		nsta
	integer		ev_cusp(MAXEVE)	! [1..MAXEVE]
	integer		ev_date(MAXEVE)	! [1..MAXEVE]
	integer		ev_time(MAXEVE)	! [1..MAXEVE]
	real		ev_mag(MAXEVE)	! [1..MAXEVE]
	real		ev_lat(MAXEVE)	! [1..MAXEVE]
	real		ev_lon(MAXEVE)	! [1..MAXEVE]
	real		ev_dep(MAXEVE)	! [1..MAXEVE]
	real		ev_x(MAXEVE)	! [1..MAXEVE]
	real		ev_y(MAXEVE)	! [1..MAXEVE]
	real		ev_z(MAXEVE)	! [1..MAXEVE]
	real		ev_herr(MAXEVE)	! [1..MAXEVE]
	real		ev_zerr(MAXEVE)	! [1..MAXEVE]
	real		ev_res(MAXEVE)	! [1..MAXEVE]
        integer         ev_fix(MAXEVE)  ! [1..MAXEVE]
	integer		src_cusp(MAXEVE)! [1..MAXEVE]
	doubleprecision	src_lat(MAXEVE)	! [1..MAXEVE]
	doubleprecision	src_lon(MAXEVE)	! [1..MAXEVE]
	real		src_dep(MAXEVE)	! [1..MAXEVE]
	real		src_lat0(MAXEVE)! [1..MAXEVE]
	real		src_lon0(MAXEVE)! [1..MAXEVE]
	real		src_x(MAXEVE)	! [1..MAXEVE]
	real		src_y(MAXEVE)	! [1..MAXEVE]
	real		src_z(MAXEVE)	! [1..MAXEVE]
	real		src_t(MAXEVE)	! [1..MAXEVE]
	real		src_x0(MAXEVE)	! [1..MAXEVE]
	real		src_y0(MAXEVE)	! [1..MAXEVE]
	real		src_z0(MAXEVE)	! [1..MAXEVE]
	real		src_t0(MAXEVE)	! [1..MAXEVE]
	character	sta_lab(MAXSTA)*7! [1..MAXSTA]
	real		sta_lat(MAXSTA)	! [1..MAXSTA]
	real		sta_lon(MAXSTA)	! [1..MAXSTA]
	real		sta_elv(MAXSTA)	! [1..MAXSTA]
	integer		sta_mod(MAXSTA)	! [1..MAXSTA]
	real		sta_dist(MAXSTA)! [1..MAXSTA]
	real		sta_az(MAXSTA)	! [1..MAXSTA]
	real		sta_rmsc(MAXSTA)! [1..MAXSTA]
	real		sta_rmsn(MAXSTA)! [1..MAXSTA]
	integer		sta_np(MAXSTA)	! [1..MAXSTA]
	integer		sta_ns(MAXSTA)	! [1..MAXSTA]
	integer		sta_nnp(MAXSTA)	! [1..MAXSTA]
	integer		sta_nns(MAXSTA)	! [1..MAXSTA]
	character	dt_sta(MAXDATA)*7! [1..MAXDATA]
	integer		dt_c1(MAXDATA)	! [1..MAXDATA]
	integer		dt_c2(MAXDATA)	! [1..MAXDATA]
	integer		dt_idx(MAXDATA)	! [1..MAXDATA]
	real		dt_dt(MAXDATA)	! [1..MAXDATA]
	real		dt_qual(MAXDATA)! [1..MAXDATA]
	real		dt_cal(MAXDATA)	! [1..MAXDATA]
	integer		dt_ista(MAXDATA)! [1..MAXDATA]
	integer		dt_ic1(MAXDATA)	! [1..MAXDATA]
	integer		dt_ic2(MAXDATA)	! [1..MAXDATA]
	real		dt_res(MAXDATA)	! [1..MAXDATA]
	real		dt_wt(MAXDATA)	! [1..MAXDATA]
	real		dt_offs(MAXDATA)! [1..MAXDATA]
	real		tmp_ttp(MAXSTA,MAXEVE)! [1..MAXSTA,1..MAXEVE]
	real		tmp_tts(MAXSTA,MAXEVE)! [1..MAXSTA,1..MAXEVE]
	real		tmp_xp(MAXSTA,MAXEVE)! [1..MAXSTA,1..MAXEVE]
	real		tmp_yp(MAXSTA,MAXEVE)! [1..MAXSTA,1..MAXEVE]
	real		tmp_zp(MAXSTA,MAXEVE)! [1..MAXSTA,1..MAXEVE]
	real		tmp_xs(MAXSTA,MAXEVE)! [1..MAXSTA,1..MAXEVE]
	real		tmp_ys(MAXSTA,MAXEVE)! [1..MAXSTA,1..MAXEVE]
	real		tmp_zs(MAXSTA,MAXEVE)! [1..MAXSTA,1..MAXEVE]
	integer		nct
	integer		ncc

c	Local variables:
	integer		i, j, k
	integer		icusp(MAXEVE)	! [1..nev] Event keys
	integer		ifindi
	integer		iicusp(MAXEVE)	! [1..nev] Index table into ev_cusp[]
	integer		nccold
	integer		nctold
	integer		ndtold
	integer		sta_itmp(MAXSTA)

      write(log,'("skipping data...")')

c     Skip data with large resiudals
      if (kiter.eq.1) then
          ndtold = ndt
          nccold = 0
          nctold = 0
      endif
      ncc = 0
      nct = 0
      j = 1
      do i=1,ndt
         if (kiter.eq.1) then
            if (dt_idx(i).le.2) then
               nccold = nccold+1
            else
               nctold = nctold+1
            endif
         endif
         if (dt_wt(i).ge.minwght) then
            dt_sta(j) = dt_sta(i)
            dt_c1(j) = dt_c1(i)
            dt_c2(j) = dt_c2(i)
            dt_idx(j) = dt_idx(i)
            dt_qual(j) = dt_qual(i)
            dt_dt(j) = dt_dt(i)
            dt_cal(j) = dt_cal(i)
            dt_res(j) = dt_res(i)
            dt_wt(j) = dt_wt(i)
            dt_offs(j) = dt_offs(i)
            if (dt_idx(i).le.2) then
                ncc = ncc+1
            else
                nct = nct+1
            endif
            j = j+1
         endif
      enddo
      ndt = j-1
      write(log,'("# obs = ",i9," (",f5.1,"%)")')
     &ndt, (ndt*100.0/ndtold)
      if (nccold.gt.0.and.nctold.gt.0) then
         write(log,'("# obs cc = ",i9," (",f5.1,"%)")')
     &   ncc, (ncc*100.0/nccold)
         write(log,'("# obs ct = ",i9," (",f5.1,"%)")')
     &   nct, (nct*100.0/nctold)
      endif

c     Skip events
      do i=1,ndt
         dt_ic1(i) = dt_c1(i)   !dt_ic1 is just a workspace array here!
         dt_ic2(i) = dt_c2(i)   !dt_ic2 is just a workspace array here!
      enddo
      call sorti(ndt, dt_ic1)
      call sorti(ndt, dt_ic2)
      k = 1
      do i=1,nev
         if (ifindi(ndt, dt_ic1, ev_cusp(i)).gt.0 .or.
     &       ifindi(ndt, dt_ic2, ev_cusp(i)).gt.0) then
            ev_date(k) = ev_date(i)
            ev_time(k) = ev_time(i)
            ev_cusp(k) = ev_cusp(i)
            ev_lat(k) = ev_lat(i)
            ev_lon(k) = ev_lon(i)
            ev_dep(k) = ev_dep(i)
            ev_mag(k) = ev_mag(i)
            ev_herr(k) = ev_herr(i)
            ev_zerr(k) = ev_zerr(i)
            ev_res(k) = ev_res(i)
            ev_fix(k) = ev_fix(i)
            ev_x(k) = ev_x(i)
            ev_y(k) = ev_y(i)
            ev_z(k) = ev_z(i)
            k = k+1
         endif
      enddo
      nev = k-1
      write(log,'("# events = ",i9)') nev

c     Skip sources
c     Uses sorted dt_ic[12] arrays from above
      if (nsrc.ne.1) then
         k = 1
         do i=1,nsrc
            if (ifindi(ndt, dt_ic1, src_cusp(i)).gt.0 .or.
     &          ifindi(ndt, dt_ic2, src_cusp(i)).gt.0) then
               src_cusp(k) = src_cusp(i)
               src_lat(k) = src_lat(i)
               src_lon(k) = src_lon(i)
               src_lat0(k) = src_lat0(i)
               src_lon0(k) = src_lon0(i)
               src_dep(k) = src_dep(i)
               src_x(k) = src_x(i)
               src_y(k) = src_y(i)
               src_z(k) = src_z(i)
               src_t(k) = src_t(i)
               src_x0(k) = src_x0(i)
               src_y0(k) = src_y0(i)
               src_z0(k) = src_z0(i)
               src_t0(k) = src_t0(i)
               do j=1,nsta
                  tmp_ttp(j,k) = tmp_ttp(j,i)
                  tmp_tts(j,k) = tmp_tts(j,i)
                  tmp_xp(j,k) = tmp_xp(j,i)
                  tmp_yp(j,k) = tmp_yp(j,i)
                  tmp_zp(j,k) = tmp_zp(j,i)
                  tmp_xs(j,k) = tmp_xs(j,i)
                  tmp_ys(j,k) = tmp_ys(j,i)
                  tmp_zs(j,k) = tmp_zs(j,i)
               enddo
               k = k+1
            endif
         enddo
         nsrc = k-1
      endif

c    Clean stations
      do i=1,nsta
         sta_itmp(i) = 0
      enddo
      do j=1,ndt
         do i=1,nsta
            if (dt_sta(j).eq.sta_lab(i)) then
               sta_itmp(i) = 1
               goto 200	! break
            endif
         enddo
200      continue
      enddo
      k = 1
      do i=1,nsta
         if (sta_itmp(i).eq.1) then
            sta_lab(k) = sta_lab(i)
            sta_lat(k) = sta_lat(i)
            sta_lon(k) = sta_lon(i)
            sta_elv(k) = sta_elv(i)
            sta_mod(k) = sta_mod(i)
            sta_dist(k) = sta_dist(i)
            sta_az(k) = sta_az(i)
            sta_np(k) = sta_np(i)
            sta_ns(k) = sta_ns(i)
            sta_nnp(k) = sta_nnp(i)
            sta_nns(k) = sta_nns(i)
            sta_rmsc(k) = sta_rmsc(i)
            sta_rmsn(k) = sta_rmsn(i)
            do j=1,nsrc
               tmp_ttp(k,j) = tmp_ttp(i,j)
               tmp_tts(k,j) = tmp_tts(i,j)
               tmp_xp(k,j) = tmp_xp(i,j)
               tmp_yp(k,j) = tmp_yp(i,j)
               tmp_zp(k,j) = tmp_zp(i,j)
               tmp_xs(k,j) = tmp_xs(i,j)
               tmp_ys(k,j) = tmp_ys(i,j)
               tmp_zs(k,j) = tmp_zs(i,j)
            enddo
            k = k+1
         endif
      enddo
      nsta = k-1
      write(log,'("# stations = ",i9)') nsta

c     Index station labels and cuspids
      call indexxi(nev, ev_cusp, iicusp)
      do i=1,nev
         icusp(i) = ev_cusp(iicusp(i)) !icusp is just a workspace array here!
      enddo
      do i=1,ndt
         do j=1,nsta
            if (dt_sta(i).eq.sta_lab(j)) then
               dt_ista(i) = j
               dt_ic1(i) = iicusp(ifindi(nev, icusp, dt_c1(i)))
               dt_ic2(i) = iicusp(ifindi(nev, icusp, dt_c2(i)))
               goto 300	! continue 2
            endif
         enddo
         write(*,'("FATAL ERROR (indexing). Please report to ",
     &             "felixw@ldeo.columbia.edu")')
         stop   
300      continue
      enddo

      end !of subroutine skip
