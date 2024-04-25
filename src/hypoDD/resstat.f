c Calculate residual statistics

	subroutine resstat(log, idata, ndt, nev, d, w, idx,
     &	rms_cc, rms_ct, rms_cc0, rms_ct0,
     &	rms_ccold, rms_ctold, rms_cc0old, rms_ct0old,
     &                   dum)

	implicit none

	include "hypoDD.inc"

c	Parameters:
	integer	log
	integer	idata
	integer	ndt
	integer	nev
	real	d(MAXDATA)	! (1..ndt)
	real	w(MAXDATA)	! (1..ndt)
	integer	idx(MAXDATA)	! (1..ndt)
	real	rms_cc
	real	rms_ct
	real	rms_cc0
	real	rms_ct0
	real	rms_ccold
	real	rms_ctold
	real	rms_cc0old
	real	rms_ct0old
	real	dum

c	Local variables:
	real	av_cc
	real	av_cc0
	real	av_ct
	real	av_ct0
	real	dav
	real	dav0
	real	dav1
	real	dav1old
	real	davold
	real	dvar
	real	dvar1
	real	dvar1old
	real	dvarold
	real	f
	real	f_cc
	real	f_ct
	integer	i, j
	real	s
	real	s1
	real	ss
	real	ss1
	real	sw
	real	sw_cc
	real	sw_ct

c--- get rms:
      rms_cc0old= rms_cc0
      rms_ct0old= rms_ct0
      rms_ccold= rms_cc
      rms_ctold= rms_ct
      j= 0
      sw_cc= 0
      sw_ct= 0
      do i= 1,ndt
         if(idx(i).le.2) then
             sw_cc= sw_cc + w(i)
             j= j+1
         else
             sw_ct= sw_ct + w(i)
         endif
      enddo
      f_cc= j/sw_cc  ! factor to scale weights for rms value
      f_ct= (ndt-j)/sw_ct  ! factor to scale weights for rms value

      rms_cc0= 0
      rms_ct0= 0
      av_cc0= 0
      av_ct0= 0
      rms_cc= 0
      rms_ct= 0
      av_cc= 0
      av_ct= 0
      j= 0
      do i= 1,ndt
         if(idx(i).le.2) then
             rms_cc0= rms_cc0 + d(i)**2
             av_cc0= av_cc0 + d(i)
             rms_cc= rms_cc + (f_cc*w(i)*d(i))**2	! weighted and scaled
             av_cc= av_cc + f_cc*w(i)*d(i)   	! weighted and scaled
             j= j+1
         else
             rms_ct0= rms_ct0 + d(i)**2
             av_ct0= av_ct0 + d(i)
             rms_ct= rms_ct + (f_ct*w(i)*d(i))**2	! weighted and scaled
             av_ct= av_ct + f_ct*w(i)*d(i)   	! weighted and scaled
         endif
      enddo
      av_cc0= av_cc0/j
      av_ct0= av_ct0/(ndt-j)
      rms_cc0= sqrt( (rms_cc0 - av_cc0**2/j) / (j-1) )
      rms_ct0= sqrt( (rms_ct0 - av_ct0**2/(ndt-j)) / (ndt-j-1) )
      av_cc= av_cc/j
      av_ct= av_ct/(ndt-j)
      rms_cc= sqrt( (rms_cc - av_cc**2/j) / (j-1) )
      rms_ct= sqrt( (rms_ct - av_ct**2/(ndt-j)) / (ndt-j-1) )

c--- more: residual average, rms, and variance:
      davold= dav
      dvarold= dvar
      dav1old= dav1
      dvar1old= dvar1
      dav= 0    ! unweighted
      dvar= 0   ! unweighted
      dav1= 0   ! weighted
      dvar1= 0  ! weighted
      sw= 0
      dav0= 0    ! unweighted

      do i= 1,ndt
         sw= sw + w(i)
      enddo
      f= ndt/sw  ! factor to scale weights for rms value

      do i= 1,ndt
         dav= dav + d(i)		! unweighted
         dav0= dav0 + w(i)*d(i)   	! weighted
         dav1= dav1 + f*w(i)*d(i)   	! weighted and scaled
      enddo
      dav= dav/ndt
      dav0= dav0/ndt
      dav1= dav1/ndt

      s= 0
      ss= 0
      ss1= 0
      do i=1,ndt
         s= d(i)*1000 - dav*1000 		! in msec
         s1= w(i)*d(i)*1000 - dav0*1000 	! weighted, in msec
         ss= ss + s
         ss1= ss1 + s1
         dvar= dvar + s**2
         dvar1= dvar1 + s1**2
      enddo
      if(ndt.gt.4*nev) then
         dvar= (dvar - ss**2/ndt) / (ndt - 4*nev) ! / by the # of deg of freedom
         dvar1= (dvar1 - ss1**2/ndt) / (ndt - 4*nev)
      else
         dvar= dvar / 1   ! / by the # of degrees of freedom
         dvar1= dvar1 / 1
c         write(*,*)'>>> ndt < 4*nev'
c         write(log,*)'>>> ndt < 4*nev'
         write(*,*)'>>> Warning: ndt < 4*nev'
         write(log,*)'>>> Warning: ndt < 4*nev'
      endif

      if(abs(dum+999).lt.0.0001) then   ! orginal data
         write(log,'(/,"Residual summary of initial data:")')
         write(log,'(a,f7.4)')' absolute mean [s] = ',dav
         write(log,'(a,f7.4)')' weighted mean [s] = ',dav1
         write(log,'(a,f10.4)')' absolute variance [s] = ',dvar/1000
         write(log,'(a,f7.4)')' weighted variance [s] = ',dvar1/1000
         if(idata.eq.1.or.idata.eq.3) then
           write(log,'(a,f7.4)')' absolute cc rms [s] = ',rms_cc0
           write(log,'(a,f7.4)')' weighted cc rms [s] (RMSCC) = ',rms_cc
         endif
         if(idata.eq.2.or.idata.eq.3) then
           write(log,'(a,f7.4)')' absolute ct rms [s] = ',rms_ct0
           write(log,'(a,f7.4)')' weighted ct rms [s] (RMSCT) = ',rms_ct
         endif
      else
         write(log,'(/,"Residual summary:")')
         write(log,'(a,f7.4,a,f7.2,a)')' absolute mean [s] = ',dav,
     & ' (',(dav-davold)*100/abs(davold),' %)'
         write(log,'(a,f7.4,a,f7.2,a)')' weighted mean [s] = ',dav1,
     & ' (',(dav1-dav1old)*100/abs(dav1old),' %)'
         write(log,'(a,f10.4,a,f7.2,a)')' absolute variance [s] = ',
     & dvar/1000, ' (',(dvar-dvarold)*100/dvarold,' %)'
         write(log,'(a,f10.4,a,f7.2,a)')' weighted variance [s] = ',
     & dvar1/1000,' (',(dvar1-dvar1old)*100/dvar1old,' %)'
         if(idata.eq.1.or.idata.eq.3) then
           write(log,'(a,f7.4,a,f7.2,a)')' absolute cc rms [s] = ',
     & rms_cc0,' (',(rms_cc0-rms_cc0old)*100/rms_cc0old,' %)'
           write(log,'(a,f7.4,a,f7.2,a)')' weighted cc rms [s] = ',
     & rms_cc,' (',(rms_cc-rms_ccold)*100/rms_ccold,' %)'
         endif
         if(idata.eq.2.or.idata.eq.3) then
           write(log,'(a,f7.4,a,f7.2,a)')' absolute ct rms [s] = ',
     & rms_ct0,' (',(rms_ct0-rms_ct0old)*100/rms_ct0old,' %)'
           write(log,'(a,f7.4,a,f7.2,a)')' weighted ct rms [s] = ',
     & rms_ct, ' (',(rms_ct-rms_ctold)*100/rms_ctold,' %)'
         endif
      endif

      dum= dvar1
      end !of tine resstat
