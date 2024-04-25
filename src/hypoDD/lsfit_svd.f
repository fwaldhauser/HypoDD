	subroutine lsfit_svd(log, iter, ndt, nev, nsrc, damp, mod_ratio,
     &	idata, ev_cusp, src_cusp,ev_fix,
     &	dt_res, dt_wt,
     &	dt_ista, dt_ic1, dt_ic2,   
     &	src_dx, src_dy, src_dz, src_dt, src_ex, src_ey, src_ez, src_et,
     &	exav, eyav, ezav, etav, dxav, dyav, dzav, dtav,
     &	rms_cc, rms_ct, rms_cc0, rms_ct0,
     &	rms_ccold, rms_ctold, rms_cc0old, rms_ct0old,
     &	tmp_xp, tmp_yp, tmp_zp,
     &	tmp_xs, tmp_ys, tmp_zs, dt_idx)

	implicit none

	include 'hypoDD.inc'

c	Parameters:
	integer		log
	integer		iter
	integer		ndt
	integer		nev
	integer		nsrc
	real		damp
c	real		mod_ratio
	real		mod_ratio(MAXLAY)
	integer		idata
	integer		ev_cusp(MAXEVE)	! (1..MAXEVE)
        integer         ev_fix(MAXEVE)  ! (1..MAXEVE)
	integer		src_cusp(MAXEVE)	! (1..MAXEVE)
	real		dt_res(MAXDATA)	! (1..MAXDATA)
	real		dt_wt(MAXDATA)	! (1..MAXDATA)
	integer		dt_ista(MAXDATA)	! (1..MAXDATA)
	integer		dt_ic1(MAXDATA)	! (1..MAXDATA)
	integer		dt_ic2(MAXDATA)	! (1..MAXDATA)
	real		src_dx(MAXEVE)	! (1..MAXEVE)
	real		src_dy(MAXEVE)	! (1..MAXEVE)
	real		src_dz(MAXEVE)	! (1..MAXEVE)
	real		src_dt(MAXEVE)	! (1..MAXEVE)
	real		src_ex(MAXEVE)	! (1..MAXEVE)
	real		src_ey(MAXEVE)	! (1..MAXEVE)
	real		src_ez(MAXEVE)	! (1..MAXEVE)
	real		src_et(MAXEVE)	! (1..MAXEVE)
	real		exav
	real		eyav
	real		ezav
	real		etav
	real		dxav
	real		dyav
	real		dzav
	real		dtav
	real		rms_cc
	real		rms_ct
	real		rms_cc0
	real		rms_ct0
	real		rms_ccold
	real		rms_ctold
	real		rms_cc0old
	real		rms_ct0old
	real		tmp_xp(MAXSTA,MAXEVE)! (1..maxsta, 1..MAXEVE)
	real		tmp_yp(MAXSTA,MAXEVE)! (1..maxsta, 1..MAXEVE)
	real		tmp_zp(MAXSTA,MAXEVE)! (1..maxsta, 1..MAXEVE)
	real		tmp_xs(MAXSTA,MAXEVE)! (1..maxsta, 1..MAXEVE)
	real		tmp_ys(MAXSTA,MAXEVE)! (1..maxsta, 1..MAXEVE)
	real		tmp_zs(MAXSTA,MAXEVE)! (1..maxsta, 1..MAXEVE)
	integer		dt_idx(MAXDATA)	! (1..MAXDATA)

c	Local	variables:
	real		cvm(MAXEVE0*4,MAXEVE0*4)
	real		dd(MAXDATA0)
	real		d(MAXDATA0)
	real		dxavold, dyavold, dzavold, dtavold
	real		exavold, eyavold, ezavold, etavold
	real		factor
	doubleprecision		g(MAXDATA0,MAXEVE0*4)
	integer		i, j, k
	integer		izero
	integer		k1, k2
	integer		nndt
	real		norm(MAXEVE*4)
	real		norm_test(MAXEVE*4)
	doubleprecision		q(MAXEVE0*4)
	real		qmin, qmax
	real		resvar1
	real		s
	real		se(MAXEVE0*4)
	real		tmp(MAXEVE0*4)
	doubleprecision		u(MAXDATA0,MAXEVE0*4)
	doubleprecision		v(MAXEVE0*4,MAXEVE0*4)
	real		wtinv(MAXDATA0)
	real		wt(MAXDATA0)
	real		x(MAXEVE0*4)

      if (ndt.gt.MAXDATA0-4) stop'>>> Increase MAXDATA0 in hypoDD.inc.'
      if (nev.gt.MAXEVE0) stop'>>> Increase MAXEVE0 in hypoDD.inc.'

c     SVD
c     Set up full G matrix
      do i=1,ndt
         do j=1,nev*4
            g(i,j)= 0
         enddo
      enddo

      do i=1,ndt
         if (nsrc.eq.1) then
            k1 = 1
            k2 = 1
         else
            k1 = dt_ic1(i)
            k2 = dt_ic2(i)
         endif
         if (dt_idx(i).eq.2 .or. dt_idx(i).eq.4) then
c           g(i,dt_ic1(i)*4-3) = tmp_xp(dt_ista(i),k1)*mod_ratio
c           g(i,dt_ic1(i)*4-2) = tmp_yp(dt_ista(i),k1)*mod_ratio
c           g(i,dt_ic1(i)*4-1) = tmp_zp(dt_ista(i),k1)*mod_ratio
           g(i,dt_ic1(i)*4-3) = tmp_xs(dt_ista(i),k1)
           g(i,dt_ic1(i)*4-2) = tmp_ys(dt_ista(i),k1)
           g(i,dt_ic1(i)*4-1) = tmp_zs(dt_ista(i),k1)
           g(i,dt_ic1(i)*4) = 1.0
c           g(i,dt_ic2(i)*4-3) = -tmp_xp(dt_ista(i),k2)*mod_ratio
c           g(i,dt_ic2(i)*4-2) = -tmp_yp(dt_ista(i),k2)*mod_ratio
c           g(i,dt_ic2(i)*4-1) = -tmp_zp(dt_ista(i),k2)*mod_ratio
           g(i,dt_ic2(i)*4-3) = -tmp_xs(dt_ista(i),k2)
           g(i,dt_ic2(i)*4-2) = -tmp_ys(dt_ista(i),k2)
           g(i,dt_ic2(i)*4-1) = -tmp_zs(dt_ista(i),k2)
           g(i,dt_ic2(i)*4) = -1.0
         else
           g(i,dt_ic1(i)*4-3) = tmp_xp(dt_ista(i),k1)
           g(i,dt_ic1(i)*4-2) = tmp_yp(dt_ista(i),k1)
           g(i,dt_ic1(i)*4-1) = tmp_zp(dt_ista(i),k1)
           g(i,dt_ic1(i)*4) = 1.0
           g(i,dt_ic2(i)*4-3) = -tmp_xp(dt_ista(i),k2)
           g(i,dt_ic2(i)*4-2) = -tmp_yp(dt_ista(i),k2)
           g(i,dt_ic2(i)*4-1) = -tmp_zp(dt_ista(i),k2)
           g(i,dt_ic2(i)*4) = -1.0
         endif
      enddo


c     Weight data
      do i=1,ndt
c        wt[] must be same as dt_wt[], but with different dimensions, and
c        should not be changed so statistics in resstat will not be screwed up!
         wt(i) = dt_wt(i)
         if (wt(i).ne.0) then
           wtinv(i) = 1.0/wt(i)
         else
           wtinv(i) = 1.0
         endif
         d(i) = dt_res(i)*1000.0 * wt(i) ! data in ms, so results are in m
      enddo

c     Weight G matrix
      call matmult3(MAXDATA0,MAXEVE0*4,ndt,nev*4,wt,g)

c     Add four extra rows to make mean shift zero.
c     This should make the design matrix non-singular.
      do i=1,4
         d(ndt+i) = 0.0
c         wt(ndt+i) = 1.0
         do j=1,nev
c            g(ndt+i,j*4-4+i) = 1.0
c            g(ndt+i,j*4-3+i) = 0.0
c            g(ndt+i,j*4-2+i) = 0.0
c            g(ndt+i,j*4-1+i) = 0.0
            g(ndt+i,j*4-3) = 0.0
            g(ndt+i,j*4-2) = 0.0
            g(ndt+i,j*4-1) = 0.0
            g(ndt+i,j*4) = 0.0
            g(ndt+i,j*4-4+i) = 1.0      ! here is the weight!!
         enddo
      enddo
      nndt = ndt+4
      write(log,'(a)')'  extra rows added to make mean shift zero: 4'

c     add extra rows to constrain parameters according to ev_fix value.
      k= 1
      do i=1,nev
         if(ev_fix(i).ne.0) then
            if(ev_fix(i).eq.1) then	! z fix
               d(nndt+k) = 0.0
               do j=1,nev*4
                  g(nndt+k,j) = 0
               enddo
cfwrt               g(nndt+k,i*4-1) = 100.0 	! weight here!!
               g(nndt+k,i*4-1) = damp 	! weight here!! 
               k= k+1
            elseif(ev_fix(i).eq.2) then	! z,t fix 
               d(nndt+k) = 0.0
               do j=1,nev*4
                  g(nndt+k,j) = 0
               enddo
cfwrt               g(nndt+k,i*4-1) = 100.0 	! weight here!!  z
               g(nndt+k,i*4-1) = damp 	! weight here!!  z
               k= k+1
               d(nndt+k) = 0.0
               do j=1,nev*4
                  g(nndt+k,j) = 0
               enddo
cfwrt               g(nndt+k,i*4) = 100.0 	! weight here!!  t
               g(nndt+k,i*4) = damp 	! weight here!!  t
               k= k+1
            elseif(ev_fix(i).eq.3) then	! x,y,z fix 
               d(nndt+k) = 0.0
               do j=1,nev*4
                  g(nndt+k,j) = 0
               enddo
cfwrt               g(nndt+k,i*4-3) = 100.0 	! weight here!!  x
               g(nndt+k,i*4-3) = damp 	! weight here!!  x
               k= k+1
               d(nndt+k) = 0.0
               do j=1,nev*4
                  g(nndt+k,j) = 0
               enddo
cfwrt               g(nndt+k,i*4-2) = 100.0 	! weight here!!  y
               g(nndt+k,i*4-2) = damp 	! weight here!!  y
               k= k+1
               d(nndt+k) = 0.0
               do j=1,nev*4
                  g(nndt+k,j) = 0
               enddo
cfwrt               g(nndt+k,i*4-1) = 100.0 	! weight here!!  z
               g(nndt+k,i*4-1) = damp 	! weight here!!  z
               k= k+1
            elseif(ev_fix(i).eq.4) then	! x,y,z,t fix 
               d(nndt+k) = 0.0
               do j=1,nev*4
                  g(nndt+k,j) = 0
               enddo
cfwrt               g(nndt+k,i*4-3) = 100.0 	! weight here!!  x
               g(nndt+k,i*4-3) = damp 	! weight here!!  x
               k= k+1
               d(nndt+k) = 0.0
               do j=1,nev*4
                  g(nndt+k,j) = 0
               enddo
cfwrt               g(nndt+k,i*4-2) = 100.0 	! weight here!!  y
               g(nndt+k,i*4-2) = damp 	! weight here!!  y
               k= k+1
               d(nndt+k) = 0.0
               do j=1,nev*4
                  g(nndt+k,j) = 0
               enddo
cfwrt               g(nndt+k,i*4-1) = 100.0 	! weight here!!  z
               g(nndt+k,i*4-1) = damp 	! weight here!!  z
               k= k+1
               d(nndt+k) = 0.0
               do j=1,nev*4
                  g(nndt+k,j) = 0
               enddo
cfwrt               g(nndt+k,i*4) = 100.0 	! weight here!!  t
               g(nndt+k,i*4) = damp 	! weight here!!  t
               k= k+1
            elseif(ev_fix(i).eq.5) then	! t fix 
               d(nndt+k) = 0.0
               do j=1,nev*4
                  g(nndt+k,j) = 0
               enddo
cfwrt               g(nndt+k,i*4) = 100.0 	! weight here!!  t
               g(nndt+k,i*4) = damp 	! weight here!!  t
               k= k+1
            elseif(ev_fix(i).eq.6) then	! x, y fix 
               d(nndt+k) = 0.0
               do j=1,nev*4
                  g(nndt+k,j) = 0
               enddo
cfwrt               g(nndt+k,i*4-3) = 100.0 	! weight here!!  x
               g(nndt+k,i*4-3) = damp 	! weight here!!  x
               k= k+1
               d(nndt+k) = 0.0
               do j=1,nev*4
                  g(nndt+k,j) = 0
               enddo
cfwrt               g(nndt+k,i*4-2) = 100.0 	! weight here!!  y
               g(nndt+k,i*4-2) = damp 	! weight here!!  y
               k= k+1
            endif
         endif
      enddo
      nndt = nndt+k-1
      write(log,'(a,i8)')
     & '  extra rows added to constrain selected parameters: ',k-1


c     Column scaling
      do j=1,4*nev
         norm(j) = 0
      enddo
      do j=1,4*nev
         do i=1,nndt
            norm(j) = norm(j) + g(i,j)**2
         enddo
      enddo
      do j=1,nev*4
         norm(j) = sqrt(norm(j)/nndt)
      enddo
      do j=1,nev*4
         do i=1,nndt
            g(i,j) =  g(i,j) / norm(j)
         enddo
      enddo

c     Testing...
      do j=1,4*nev
         norm_test(j) = 0
      enddo
      do j=1,4*nev
         do i=1,nndt
            norm_test(j) = norm_test(j) + g(i,j)**2
         enddo
      enddo
      do j=1,nev*4
         norm_test(j) = sqrt(norm_test(j)/nndt)
         if (abs(norm_test(j)-1).gt.0.001) then
            write(*,'("FATAL ERROR (svd: G scaling). Please report to ",
     &                "felixw@ldeo.columbia.edu")')
            stop
         endif
      enddo

c     Do singular-value decomposition
      write(log,'("~ singular value decomposition of G (",i6,
     & "x",i6," matrix) ...")') nndt,nev*4
      call SVD(MAXDATA0,MAXEVE0*4,nndt,nev*4,g,u,v,q,1)

c     Check for singular values close to zero
      write(log,'("~ backsubstitution ...")')
      izero = 0
      qmax = 0.0
      do i=1,nev*4
         if (q(i).lt.0) then
            write(*,'("FATAL ERROR (svd: neg sing val). Please report ",
     &                "to felix@andreas.wr.usgs.gov")')
            stop
          endif
         if (q(i).gt.qmax) qmax = q(i)
      enddo
      qmin = qmax * 0.0000001
      do i=1,nev*4
         if (q(i).lt.qmin) then
            q(i) = 0.0
            izero = izero+1
         endif
      enddo
      if (izero.gt.0)then
         write(*,'(/,">>> ",i3," singular values close/= zero !!",/)')
     &   izero
         write(log,'(/,">>> ",i3," singular values close/= zero !!",/)')
     &   izero
      endif

c     Back substitution (get x' from Gx=d: x=v*diag(1/q)*t(u)*d))

c     Compute diag(1/q)*t(U)*d
      do j=1,nev*4
        s = 0.0
        if (q(j).ne.0) then
           do i=1,nndt
              s = s+u(i,j)*d(i)
           enddo
           s = s/q(j)
        endif
        tmp(j) = s
      enddo

c     Multiply by V
      do i=1,nev*4
        s = 0.0
        do j=1,nev*4
           s = s+v(i,j)*tmp(j)
        enddo
        x(i) = s
      enddo

c     Rescale model vector and G
      do j=1,4*nev
         x(j) = x(j) / norm(j)
         do i=1,ndt
           g(i,j) = g(i,j) * norm(j)
         enddo
      enddo

c     Unweight G matrix
      call matmult3(MAXDATA0,MAXEVE0*4,ndt,nev*4,wtinv,g)

c     Predict data dd = G*x', get residuals (sec)
      call matmult2(MAXDATA0,ndt,nev*4,g,x,dd)
      do i=1,ndt
        dt_res(i) = dt_res(i) - dd(i)/1000
      enddo

c     Get covariance matrix: cvm = v*(1/q**2)*vt
      call covar(MAXEVE0,v,nev*4,q,cvm)

c     Get residual statistics (avrg,rms,var..)
c     At this point, wt must be = dt_wt; just different dimension to save
c     memory for SVD runs (see above).
      call resstat(log,idata,ndt,nev,dt_res,dt_wt,dt_idx,
     & rms_cc,rms_ct,rms_cc0,rms_ct0,
     & rms_ccold,rms_ctold,rms_cc0old,rms_ct0old,
     &             resvar1)

c     Errors for the 95% confidence level,
c     thus multiply the standard errors by 2.7955 
c100710 Not sure anymore where this number comes from. It should be 1.96 
c100710 assuming a t-distribution of the residuals. This means that errors were 
c100710 overestimated by about 40%. see emails by Hilary Martens May, 2010.
c100710      factor = 2.7955
      factor = 1.96 
      do i=1,nev*4
         se(i) = sqrt(cvm(i,i))*sqrt(resvar1)*factor  ! Weighted variance
      enddo

c     Rescale errors
      do i=1,nev*4
         se(i) = se(i) / norm(i)
      enddo

c     Store solution and errors
      do i=1,nev
         src_dx(i) = -x(4*i-3)
         src_dy(i) = -x(4*i-2)
         src_dz(i) = -x(4*i-1)
         src_dt(i) = -x(4*i)
         src_ex(i) = se(4*i-3)
         src_ey(i) = se(4*i-2)
         src_ez(i) = se(4*i-1)
         src_et(i) = se(4*i)
         src_cusp(i) = ev_cusp(i)
      enddo

c     Output statistics....
c     Get average errors and vector changes
      exavold = exav
      eyavold = eyav
      ezavold = ezav
      etavold = etav
      dxavold = dxav
      dyavold = dyav
      dzavold = dzav
      dtavold = dtav
      exav = 0.0
      eyav = 0.0
      ezav = 0.0
      etav = 0.0
      dxav = 0.0
      dyav = 0.0
      dzav = 0.0
      dtav = 0.0
      do i=1,nev
         exav = exav + src_ex(i)
         eyav = eyav + src_ey(i)
         ezav = ezav + src_ez(i)
         etav = etav + src_et(i)
         dxav = dxav + abs(src_dx(i))
         dyav = dyav + abs(src_dy(i))
         dzav = dzav + abs(src_dz(i))
         dtav = dtav + abs(src_dt(i))
      enddo
      exav = exav/nev
      eyav = eyav/nev
      ezav = ezav/nev
      etav = etav/nev
      dxav = dxav/nev
      dyav = dyav/nev
      dzav = dzav/nev
      dtav = dtav/nev

      if (iter.eq.1) then
         exavold = exav
         eyavold = eyav
         ezavold = ezav
         etavold = etav
         dxavold = dxav
         dyavold = dyav
         dzavold = dzav
         dtavold = dtav
      endif
      write(log,'(/,"Location summary:",/,
     & "  mean 2sig-error (x,y,z,t) [m,ms]: ",/,f7.1,f7.1,f7.1,f7.1,
     & " (",f7.1,f7.1,f7.1,f7.1")",/,
     & "  mean shift (x,y,z,t) [m,ms] (DX,DY,DZ,DT): ",/,
     & f7.1,f7.1,f7.1,f7.1," (",f7.1,f7.1,f7.1,f7.1,")")')
     & exav,eyav,ezav,etav,exav-exavold,eyav-eyavold,
     & ezav-ezavold,etav-etavold,
     & dxav,dyav,dzav,dtav,dxav-dxavold,dyav-dyavold,
     & dzav-dzavold,dtav-dtavold

      end !of subroutine lsfit_SVD
