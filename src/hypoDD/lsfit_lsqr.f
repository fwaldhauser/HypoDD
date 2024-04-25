	subroutine lsfit_lsqr(log, iter, ndt, nev, nsrc,
     &	damp, mod_ratio,
     &	idata, ev_cusp, src_cusp,ev_fix,
     &	dt_res, dt_wt,
     &	dt_ista, dt_ic1, dt_ic2,
     &	src_dx, src_dy, src_dz, src_dt, src_ex, src_ey, src_ez, src_et,
     &	exav, eyav, ezav, etav, dxav, dyav, dzav, dtav,
     &	rms_cc, rms_ct, rms_cc0, rms_ct0,
     &	rms_ccold, rms_ctold, rms_cc0old, rms_ct0old,
     &	tmp_xp, tmp_yp, tmp_zp,
     &	tmp_xs, tmp_ys, tmp_zs, dt_idx, acond)

	implicit none

	include'hypoDD.inc'

c	Parameters:
	integer		log		! Log-file identifier
	integer		iter
	integer		ndt
	integer		nev
	integer		nsrc
	real		damp
c	real		mod_ratio
	real		mod_ratio(MAXLAY)
	integer		idata
	integer		ev_cusp(MAXEVE)	! [1..nev] Event keys
        integer         ev_fix(MAXEVE)  ! [1..nev]
	integer		src_cusp(MAXEVE)! [1..nev] Event keys
	real		dt_res(MAXDATA)	! [1..ndt]
	real		dt_wt(MAXDATA)	! [1..ndt]
	integer		dt_ista(MAXDATA)! [1..ndt]
	integer		dt_ic1(MAXDATA)	! [1..ndt]
	integer		dt_ic2(MAXDATA)	! [1..ndt]
	real		src_dx(MAXEVE)	! [1..nev]
	real		src_dy(MAXEVE)	! [1..nev]
	real		src_dz(MAXEVE)	! [1..nev]
	real		src_dt(MAXEVE)	! [1..nev]
	real		src_ex(MAXEVE)	! [1..nev]
	real		src_ey(MAXEVE)	! [1..nev]
	real		src_ez(MAXEVE)	! [1..nev]
	real		src_et(MAXEVE)	! [1..nev]
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
	real		tmp_xp(MAXSTA,MAXEVE)! [1.., 1..nev]
	real		tmp_yp(MAXSTA,MAXEVE)! [1.., 1..nev]
	real		tmp_zp(MAXSTA,MAXEVE)! [1.., 1..nev]
	real		tmp_xs(MAXSTA,MAXEVE)! [1.., 1..nev]
	real		tmp_ys(MAXSTA,MAXEVE)! [1.., 1..nev]
	real		tmp_zs(MAXSTA,MAXEVE)! [1.., 1..nev]
	integer		dt_idx(MAXDATA)	! [1..ndt]
	real		acond		! Condition number

c	Local variables:
	real		anorm
	real		arnorm
	real		atol
	real		btol
	real		conlim
	character	dattim*25
	real		d(MAXDATA+4)	! Data vector
	real		dtavold
	real		dxavold
	real		dyavold
	real		dzavold
	real		etavold
	real		exavold
	real		eyavold
	real		ezavold
	real		factor
	integer		i
	integer		istop
	integer		itnlim
	integer		iw(2*(8*MAXDATA+4*MAXEVE)+1)	! lsqr index array
	integer		k1
	integer		k2
	integer		leniw
	integer		lenrw
	integer		m
	integer		n
	integer		nar
	integer		nndt
	integer		nrw	
	integer		k	
	real		norm(MAXEVE*4)
	real		norm_test(MAXEVE*4)
	real		resvar1
	real		rms_cc0old
	real		rms_ccold
	real		rms_ct0old
	real		rms_ctold
	real		rnorm
	real		rw(8*MAXDATA+4*MAXEVE)
	real		se(MAXEVE*4)	! Solution error
	real		w1(MAXEVE*4)	! Work space
	real		w2(MAXEVE*4)	! Work space
	real		wtinv(MAXDATA+4)! +4 = mean shift constr
	real		wt(MAXDATA+4)
	real		x(MAXEVE*4)	! Solution vector
	real		xnorm

      write(log,'(/,"~ setting up G matrix.. ")')
      write(log,*)"     ndt, nev, nsrc, damp: ",ndt, nev, nsrc, damp

cc     If mean shift not contstrained
c      nar = 8*ndt
c      nndt = ndt

      nndt= ndt         ! current length of d (diff time) vector
      nrw=  8*ndt       ! current length of rw (partial deriv) vector

c-- Define end length of rw vector here:
      nar = 8*ndt       ! If mean shift not contstrained
c      nar= nar+4*nev   ! If mean shift contstrained (four extra rows)

c determine how many entries are needed to constrain selected hypocentral
c parameters:
      do i=1,nev
         if(ev_fix(i).eq.1) then        !z
            nar= nar+ 1
         elseif(ev_fix(i).eq.2) then    !z,t
            nar= nar+ 2
         elseif(ev_fix(i).eq.3) then    !x,y,z
            nar= nar+ 3
         elseif(ev_fix(i).eq.4) then    !x,y,z,t
            nar= nar+ 4
         elseif(ev_fix(i).eq.5) then    !t
            nar= nar+ 1
         elseif(ev_fix(i).eq.6) then    !x,y
            nar= nar+ 2
         endif
      enddo

      iw(1) = nar

c     Prepare sparse data and design vector
      do i=1,ndt
c        Weight data first
         wt(i) = dt_wt(i)
         if (wt(i).ne.0) then
            wtinv(i) = 1.0/wt(i)
         else
            wtinv(i) = 1.0
         endif
         d(i) = dt_res(i)*1000.0 * wt(i)

c        Set up column indexes with non-zero elements
         iw(1+i) = i
         iw(1+ndt+i) = i
         iw(1+2*ndt+i) = i
         iw(1+3*ndt+i) = i
         iw(1+4*ndt+i) = i
         iw(1+5*ndt+i) = i
         iw(1+6*ndt+i) = i
         iw(1+7*ndt+i) = i

c        Set up non-zero G matrix elements and apply weights
         if (nsrc.eq.1) then
            k1 = 1
            k2 = 1
         else
            k1 = dt_ic1(i)
            k2 = dt_ic2(i)
         endif
         if (dt_idx(i).eq.2 .or. dt_idx(i).eq.4) then
c            rw(i)       = tmp_xp(dt_ista(i),k1) * wt(i) * mod_ratio
c            rw(ndt+i)   = tmp_yp(dt_ista(i),k1) * wt(i) * mod_ratio
c            rw(2*ndt+i) = tmp_zp(dt_ista(i),k1) * wt(i) * mod_ratio
            rw(i)       = tmp_xs(dt_ista(i),k1) * wt(i) 
            rw(ndt+i)   = tmp_ys(dt_ista(i),k1) * wt(i) 
            rw(2*ndt+i) = tmp_zs(dt_ista(i),k1) * wt(i) 
            rw(3*ndt+i) = wt(i)
c            rw(4*ndt+i) = -tmp_xp(dt_ista(i),k2) * wt(i) * mod_ratio
c            rw(5*ndt+i) = -tmp_yp(dt_ista(i),k2) * wt(i) * mod_ratio
c            rw(6*ndt+i) = -tmp_zp(dt_ista(i),k2) * wt(i) * mod_ratio
            rw(4*ndt+i) = -tmp_xs(dt_ista(i),k2) * wt(i) 
            rw(5*ndt+i) = -tmp_ys(dt_ista(i),k2) * wt(i) 
            rw(6*ndt+i) = -tmp_zs(dt_ista(i),k2) * wt(i) 
            rw(7*ndt+i) = -wt(i)
         else
            rw(i)       = tmp_xp(dt_ista(i),k1) * wt(i)
            rw(ndt+i)   = tmp_yp(dt_ista(i),k1) * wt(i)
            rw(2*ndt+i) = tmp_zp(dt_ista(i),k1) * wt(i)
            rw(3*ndt+i) = wt(i)
            rw(4*ndt+i) = -tmp_xp(dt_ista(i),k2) * wt(i)
            rw(5*ndt+i) = -tmp_yp(dt_ista(i),k2) * wt(i)
            rw(6*ndt+i) = -tmp_zp(dt_ista(i),k2) * wt(i)
            rw(7*ndt+i) = -wt(i)
         endif

c        Set up column indexes with non-zero elements
         iw(1+nar+      i) = 4*dt_ic1(i) - 3
         iw(1+nar+  ndt+i) = 4*dt_ic1(i) - 2
         iw(1+nar+2*ndt+i) = 4*dt_ic1(i) - 1
         iw(1+nar+3*ndt+i) = 4*dt_ic1(i)
         iw(1+nar+4*ndt+i) = 4*dt_ic2(i) - 3
         iw(1+nar+5*ndt+i) = 4*dt_ic2(i) - 2
         iw(1+nar+6*ndt+i) = 4*dt_ic2(i) - 1
         iw(1+nar+7*ndt+i) = 4*dt_ic2(i)

      enddo
      nrw= 8*ndt        !current length of rw vector 

c--- add four extra rows to make mean shift zero
c  if this is activ, set nar and nndt at the beginning for this routine to
c  proper values!!!!
c      do i=1,4
c         d(ndt+i)= 0                   ! zero shift
c         wt(ndt+i)= 1                  ! weight for the constraint
c         wtinv(ndt+i)= 1/wt(ndt+i)
c         do j=1,nev
c            iw(1 + nrw +(i-1)*nev+j)= ndt +i 
c            iw(1 + nar + nrw +(i-1)*nev+j)= (i-1) + j*4-3
c            rw(nrw + (i-1)*nev+j)= 10
c         enddo
c      enddo
c      nrw= nrw+4*nev
c      nndt= nndt+4

c THIS HAS NOT BEEN TESTED!!!
c     Add extra rows to constrain parameters according to ev_fix value.
      k= 1
      do i=1,nev
         if(ev_fix(i).ne.0) then
cc	    stop'ev_fix >0 in lsqr runs not tested. Do not use.'
            if(ev_fix(i).eq.1) then			! z fix
               d(nndt+k) = 0.0
               iw(1+ nrw + k)= nndt + k 
               iw(1+ nar + nrw + k)= 4*i-1 
               rw(nrw +k) = 100.0 			! weight here!!
c               write(*,*)nndt + k,4*ev_ic(i)-1,dt_ic1(1),dt_ic2(1)
               k= k+1
            elseif(ev_fix(i).eq.2) then	! z,t fix 
               d(nndt+k) = 0.0
               iw(1+ nrw + k)= nndt + k 
               iw(1+ nar + nrw + k)= 4*i-1 	! z fix
               rw(nrw +k) = 100.0 			! weight here!!
               k= k+1
               d(nndt+k) = 0.0
               iw(1+ nrw + k)= nndt + k 
               iw(1+ nar + nrw + k)= 4*i 	! t fix
               rw(nrw +k) = 100.0 			! weight here!!
               k= k+1
            elseif(ev_fix(i).eq.3) then			! x,y,z fix 
               d(nndt+k) = 0.0
               iw(1+ nrw + k)= nndt + k 
               iw(1+ nar + nrw + k)= 4*i-3 	! x fix
cfwrt               rw(nrw +k) = 100.0 			! weight here!!
               rw(nrw +k) = 10.0 			! weight here!!
               k= k+1
               d(nndt+k) = 0.0
               iw(1+ nrw + k)= nndt + k 
               iw(1+ nar + nrw + k)= 4*i-2 	! y fix
cfwrt               rw(nrw +k) = 100.0 			! weight here!!
               rw(nrw +k) = 10.0 			! weight here!!
               k= k+1
               d(nndt+k) = 0.0
               iw(1+ nrw + k)= nndt + k 
               iw(1+ nar + nrw + k)= 4*i-1 	! z fix
cfwrt               rw(nrw +k) = 100.0 			! weight here!!
               rw(nrw +k) = 10.0 			! weight here!!
               k= k+1
            elseif(ev_fix(i).eq.4) then			! x,y,z,t fix 
               d(nndt+k) = 0.0
               iw(1+ nrw + k)= nndt + k 
               iw(1+ nar + nrw + k)= 4*i-3 	! x fix
               rw(nrw +k) = 100.0 			! weight here!!
               k= k+1
               d(nndt+k) = 0.0
               iw(1+ nrw + k)= nndt + k 
               iw(1+ nar + nrw + k)= 4*i-2 	! y fix
               rw(nrw +k) = 100.0 			! weight here!!
               k= k+1
               d(nndt+k) = 0.0
               iw(1+ nrw + k)= nndt + k 
               iw(1+ nar + nrw + k)= 4*i-1 	! z fix
               rw(nrw +k) = 100.0 			! weight here!!
               k= k+1
               d(nndt+k) = 0.0
               iw(1+ nrw + k)= nndt + k 
               iw(1+ nar + nrw + k)= 4*i 	! t fix
               rw(nrw +k) = 100.0 			! weight here!!
               k= k+1
            elseif(ev_fix(i).eq.5) then	! t fix 
               d(nndt+k) = 0.0
               iw(1+ nrw + k)= nndt + k 
               iw(1+ nar + nrw + k)= 4*i 	! t fix
               rw(nrw +k) = 1000.0 			! weight here!!
               k= k+1
            elseif(ev_fix(i).eq.6) then	! x,y fix 
               d(nndt+k) = 0.0
               iw(1+ nrw + k)= nndt + k 
               iw(1+ nar + nrw + k)= 4*i-3 	! x fix
               rw(nrw +k) = 100.0 			! weight here!!
               k= k+1
               d(nndt+k) = 0.0
               iw(1+ nrw + k)= nndt + k 
               iw(1+ nar + nrw + k)= 4*i-2 	! y fix
               rw(nrw +k) = 100.0 			! weight here!!
               k= k+1
            endif
         endif
      enddo
      nrw= nrw + k-1
      if(nrw.ne.nar)stop'FATAL ERROR (lsfit_lsqr): nrw != nar'

c     Scale G matrix so the L2 norm of each column is 1.
      write(log,'("~ scaling G columns ... ")')

c     G array scaling
      do i=1,4*nev
         norm(i) = 0.0
      enddo
      do i=1,nar
         norm(iw(1+nar+i)) = norm(iw(1+nar+i)) + rw(i)**2
      enddo
      do i=1,nev*4
         norm(i) = sqrt(norm(i)/nndt)
      enddo
      do i=1,nar
         rw(i) = rw(i) / norm(iw(1+nar+i))
      enddo

c     Testing...
      do i=1,4*nev
         norm_test(i) = 0.0
      enddo
      do i=1,nar
         norm_test(iw(1+nar+i)) = norm_test(iw(1+nar+i)) + rw(i)**2
      enddo
      do i=1,nev*4
         norm_test(i) = sqrt(norm_test(i)/nndt)
         if (abs(norm_test(i)-1).gt.0.001) then
            write(*,'("FATAL ERROR (lsqr: G scaling). Please report to ",   
     &                "felixw@ldeo.columbia.edu")')
            stop
         endif
      enddo

c     Least square fitting using the algorithm of
c     Paige and Saunders, acm-trans. math. software, vol.8, no. 2,
c     jun., 1982, p. 195.

c     Set up input parameter first
      m = nndt
      n = nev*4
      leniw = 2*nar+1
      lenrw = nar
      do i= 1,n
         w1(i) = 0.0
         w2(i) = 0.0
         x(i) = 0.0
         se(i) = 0.0
      enddo
      atol = 0.000001
      btol = 0.000001
      conlim = 100000.0
      itnlim = 100*n
      istop = 0
      anorm = 0.0
      acond = 0.0
      rnorm = 0.0
      arnorm = 0.0
      xnorm = 0.0

      call datetime(dattim)
      write(log,'("~ lsqr ...    ", a)') dattim

c d= data vector; w1,w2 = workspace; x= solution vector; se=solution error
      call lsqr(m, n, damp, 
     & leniw, lenrw, iw, rw, 
     & d, w1, w2, x, se, 
     & atol, btol, conlim, itnlim, 
     & istop, anorm, acond, rnorm, arnorm, xnorm)

      write(log,'("  istop = ",i1,"; acond (CND)=",f8.1,"; anorm =",f8.1,
     & "; arnorm =",f8.1,"; xnorm =",f8.1)')
     & istop, acond, anorm, arnorm, xnorm

      if (nsrc.eq.1) nsrc = nev

c     Rescale model vector
      do i=1,4*nev
         x(i) = x(i) / norm(i)
         se(i) = se(i) / norm(i)
      enddo

c     Unweight and rescale G matrix
      do i=1,ndt
         rw(i)       = rw(i)       * wtinv(i) * norm(iw(1+nar+i))
         rw(ndt+i)   = rw(ndt+i)   * wtinv(i) * norm(iw(1+nar+ndt+i))
         rw(2*ndt+i) = rw(2*ndt+i) * wtinv(i) * norm(iw(1+nar+2*ndt+i))
         rw(3*ndt+i) = rw(3*ndt+i) * wtinv(i) * norm(iw(1+nar+3*ndt+i))
         rw(4*ndt+i) = rw(4*ndt+i) * wtinv(i) * norm(iw(1+nar+4*ndt+i))
         rw(5*ndt+i) = rw(5*ndt+i) * wtinv(i) * norm(iw(1+nar+5*ndt+i))
         rw(6*ndt+i) = rw(6*ndt+i) * wtinv(i) * norm(iw(1+nar+6*ndt+i))
         rw(7*ndt+i) = rw(7*ndt+i) * wtinv(i) * norm(iw(1+nar+7*ndt+i))
      enddo

c     Compute residuals from d = G*x
      do i=1,ndt
         d(i) = -dt_res(i)*1000.0
      enddo
      call aprod(1, m, n, x, d, leniw, lenrw, iw, rw)
      do i=1,ndt
         dt_res(i) = -d(i)/1000.0
      enddo

c     Get residual statistics (avrg, rms, var..)
      call resstat(log, idata, ndt, nev, dt_res, wt, dt_idx, 
     & rms_cc, rms_ct, rms_cc0, rms_ct0, 
     & rms_ccold, rms_ctold, rms_cc0old, rms_ct0old, 
     &             resvar1)


c     Scale errors
c The standard error estimates returned by LSQR increase monotonically
c with the iterations.  If LSQR shuts down early because of loose tolerances,
c or because the rhs-vector is special, the estimates will be too small.
c (I think they are most likely to be accurate if the rhs is random.)
c
c Remember that se(j) is covariance(j) / (m - n)
c where m - n = 1000000.  I've never quite understood why we
c divide by that number.

c Errors for the 95% confidence level,
c thus multiply the standard errors by 2.7955
      factor = 2.7955

c     Store solution and errors
      do i=1,nev
         src_dx(i) = -x(4*i-3)
         src_dy(i) = -x(4*i-2)
         src_dz(i) = -x(4*i-1)
         src_dt(i) = -x(4*i)
         src_cusp(i) = ev_cusp(i)

c        Take weighted variance
         src_ex(i) = sqrt(se(4*i-3)) * sqrt(resvar1) *factor
         src_ey(i) = sqrt(se(4*i-2)) * sqrt(resvar1) *factor
         src_ez(i) = sqrt(se(4*i-1)) * sqrt(resvar1) *factor
         src_et(i) = sqrt(se(4*i))   * sqrt(resvar1) *factor
      enddo

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

c     Output location statistics
      write(log,'(/,"Location summary:")')
      write(log,'(
     & " mean 2sig-error (x,y,z,t) [m,ms]: ",/,f7.1,f7.1,f7.1,f7.1,
     & " (",f7.1,f7.1,f7.1,f7.1")",/,
     & " mean shift (x,y,z,t) [m,ms] (DX,DY,DZ,DT): ",/,
     & f7.1,f7.1,f7.1,f7.1," (",f7.1,f7.1,f7.1,f7.1,")")')
     & exav, eyav, ezav, etav, exav-exavold, eyav-eyavold, 
     & ezav-ezavold, etav-etavold, 
     & dxav, dyav, dzav, dtav, dxav-dxavold, dyav-dyavold, 
     & dzav-dzavold, dtav-dtavold

      end !of subroutine lsfit_lsqr
