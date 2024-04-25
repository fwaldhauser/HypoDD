c 10/04/00 determines a priori weights and re-weights them ....

	subroutine weighting(log, ndt, mbad, amcusp, idata, kiter, ineg,
     &	maxres_cross, maxres_net, maxdcc, maxdct, minwght,
     &	wt_ccp, wt_ccs, wt_ctp, wt_cts,
     &	dt_c1, dt_c2, dt_idx, dt_qual, dt_res, dt_offs,
     &	dt_wt)

	implicit none

	include "hypoDD.inc"

c	Parameters:
	integer	log
	integer	ndt
	integer	mbad
	integer	amcusp(1000)
	integer	idata
	integer	kiter
	integer	ineg
	real	maxres_cross
	real	maxres_net
	real	maxdcc
	real	maxdct
	real	minwght
	real	wt_ccp
	real	wt_ccs
	real	wt_ctp
	real	wt_cts
	integer	dt_c1(MAXDATA)		! (1..MAXDATA)
	integer	dt_c2(MAXDATA)		! (1..MAXDATA)
	integer	dt_idx(MAXDATA)		! (1..MAXDATA)
	real	dt_qual(MAXDATA)		! (1..MAXDATA)
	real	dt_res(MAXDATA)		! (1..MAXDATA)
	real	dt_offs(MAXDATA)		! (1..MAXDATA)
	real	dt_wt(MAXDATA)		! (1..MAXDATA)

c	Local variables:
	character	dattim*25
	real	dt_tmp(MAXDATA)
	integer	i, j, k
	real	mad_cc
	real	mad_ct
	real	maxres_cc
	real	maxres_ct
	real	med_cc
	real	med_ct
	integer	ncc
	integer	nct
	integer	nncc
	integer	nnct

        call datetime(dattim)

c synthetics:
      if(idata.eq.0) then
            do i=1,ndt
               dt_wt(i)=1 	
            enddo
      endif

c--- get a priori data weights:
c intial (a priori) weights:
c s=linspace(0.0,100.0,101); ss= (s/100).^2; plot(s/100,ss); % coherency
c s=linspace(0.0,2.0,101); ss= (1./(2.^s)); plot(s,ss); % pick qual

c all the quality transf is done in getdata. old format listed qualities,
c new format list weights directly.
      ineg= 0  		!flag, =1 if neg weights exist
      do i=1,ndt
          if(dt_idx(i).eq.1)
     & dt_wt(i)= wt_ccp * dt_qual(i)	! compat. with new format
          if(dt_idx(i).eq.2)
     & dt_wt(i)= wt_ccs * dt_qual(i)	! compat. with new format
          if(dt_idx(i).eq.3)
     & dt_wt(i)= wt_ctp * dt_qual(i) ! compatib with new format 17/01/00
          if(dt_idx(i).eq.4)
     & dt_wt(i)= wt_cts * dt_qual(i) ! compatib with new format 17/01/00

          do j=1,mbad
             if(dt_c1(i).eq.amcusp(j).or.dt_c2(i).eq.amcusp(j)) then
                 dt_wt(i)= 0.0
                 ineg= 1
             endif
          enddo
      enddo

c--- re-weighting: :
      if(((idata.eq.1.or.idata.eq.3).and.
     &    (maxres_cross.ne.-999.or.maxdcc.ne.-999)).or.
     &   ((idata.eq.2.or.idata.eq.3).and.
     &    (maxres_net.ne.-999.or.maxdct.ne.-999))) then
          write(log,'("re-weighting ... ", a)') dattim

c--- get median and MAD of residuals
         if(idata.eq.3) then
            if(maxres_cross.ge.1) then
c cross data:
               k= 1
               do i=1,ndt
                  if(dt_idx(i).le.2) then
                     dt_tmp(k)= dt_res(i)
                     k= k+1
                  endif
               enddo
               call mdian1(dt_tmp,k-1,med_cc)
c 071200...
               do i=1,k-1
                  dt_tmp(i)= abs(dt_tmp(i)-med_cc)
               enddo
               call mdian1(dt_tmp,k-1,mad_cc)
               mad_cc= mad_cc/0.67449    !MAD for gaussian noise
            endif
            if(maxres_net.ge.1) then
c- catalog data:
               k= 1
               do i=1,ndt
                 if(dt_idx(i).ge.3) then
                    dt_tmp(k)= dt_res(i)
                    k= k+1
                 endif
               enddo
               call mdian1(dt_tmp,k-1,med_ct)
               do i=1,k-1
                  dt_tmp(i)= abs(dt_tmp(i)-med_ct)
               enddo
               call mdian1(dt_tmp,k-1,mad_ct)
               mad_ct= mad_ct/0.67449    !MAD for gaussian noise
            endif

         elseif((idata.eq.1.and.maxres_cross.ge.1).or.
     &       (idata.eq.2.and.maxres_net.ge.1)) then
            do i=1,ndt
                  dt_tmp(i)= dt_res(i)
            enddo
            call mdian1(dt_tmp,ndt,med_cc)
c new... see email from jim larsen...
            do i=1,ndt
                dt_tmp(i)= abs(dt_tmp(i)-med_cc)
            enddo
            call mdian1(dt_tmp,ndt,mad_cc)
            mad_cc= mad_cc/0.67449
            if(idata.eq.2) mad_ct= mad_cc
         endif

c--- define residual cutoff value:
         maxres_cc= maxres_cross	! absolute cutoff value
         maxres_ct= maxres_net		! absolute cutoff value
         if(maxres_cross.ge.1) maxres_cc= mad_cc*maxres_cross
         if(maxres_net.ge.1) maxres_ct= mad_ct*maxres_net

c--- apply residual/offset dependent weights to a priori weights
         nncc= 0
         nnct= 0
         ncc= 0
         nct= 0
         do i=1,ndt
            if(dt_idx(i).le.2) then
c--- cross data:
               ncc= ncc+1

c bi ^5 offset weighting for cross data:
c    exp needs to be uneven so weights become negative for offsets larger
c    than 2 km. 2km is hardwired, >>>not anymore since 03/23/00
c s=linspace(0,2.2,30);ss=(1-(s/2).^5).^5;plot(s,ss);axis([0 2.0 -0.0 1]);
               if(maxdcc.ne.-999)
     & dt_wt(i)= dt_wt(i) * (1 - (dt_offs(i)/(maxdcc*1000))**5)**5

c bi-cube residual weighting:
c     needs to be cube so that res > cutoff become negative.
c s=linspace(-0.2,0.2,101); ss= (1- (abs(s)/0.1).^3).^3;
c plot(abs(s),ss);  axis([0 0.11 -0.1 1]);

               if(maxres_cross.gt.0.and.dt_wt(i).gt.0.000001)
     & dt_wt(i)= dt_wt(i) * (1- (abs(dt_res(i))/maxres_cc)**3)**3
               if(dt_wt(i).lt.minwght) nncc= nncc+1

            else	
c--- catalog data:
               nct= nct+1

c bi ^3 offset weighting for catalog data:
c    exp needs to be uneven so weights become negative for offsets larger
c    than 10 km. 10km is hardwired. not anymore since 03/23/00
c s=linspace(0,11,100);ss=(1-(s/10).^3).^3;plot(s,ss);axis([0 11 -0.1 1]);
               if(maxdct.ne.-999)
     &  dt_wt(i)= dt_wt(i) * (1 - (dt_offs(i)/(maxdct*1000))**3)**3

c bi-cube residual weighting:
c     needs to be cube so that res > cutoff become negative.
c s=linspace(-0.2,0.2,101); ss= (1- (abs(s)/0.1).^3).^3;
c plot(abs(s),ss);  axis([0 0.11 -0.1 1]);
               if(dt_wt(i).gt.0.000001 .and. maxres_net.gt.0)
     & dt_wt(i)= dt_wt(i) * (1- (abs(dt_res(i))/maxres_ct)**3)**3	
               if(dt_wt(i).lt.minwght) nnct= nnct+1
            endif
         enddo

c--- check if neg residuals exist
         ineg= 0
         do j=1,ndt
            if(dt_wt(j).lt.minwght) then
               ineg= 1
               goto  100
            endif
         enddo
100      continue

         if(idata.eq.1.or.idata.eq.3) then
            write(log,'(" cc res/dist cutoff:",
     & f7.3,"s/",f6.2,"km (",f5.1,"%)")')
     & maxres_cc,maxdcc,(nncc*100.0/ncc)
         endif
         if(idata.eq.2.or.idata.eq.3) then
            write(log,'(" ct res/dist cutoff [s]:",
     & f7.3,"s/",f6.2,"km (",f5.1,"%)")')
     & maxres_ct,maxdct,(nnct*100.0/nct)
         endif
      endif  ! re-weighting

      if(ineg.gt.0) kiter= kiter+1
      end ! of subroutine weighting
