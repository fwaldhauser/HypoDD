	subroutine dtres(log, ndt, stdim, nsrc,
     &	dt_dt, dt_idx,
     &	dt_ista, dt_ic1, dt_ic2,
     &	src_cusp, src_t, tmp_ttp, tmp_tts,
     &	dt_cal, dt_res)

	implicit none

	include'hypoDD.inc'

c	Parameters:
	integer		log		! Log-file identifier
	integer		ndt		! No. of data
	integer		stdim		! Column dimenson of arrays tmp_tt[ps]
	integer		nsrc		! No. of sources
	real		dt_dt(MAXDATA)	! [1..ndt] Observed time differences
	integer		dt_idx(MAXDATA)	! [1..ndt]
	integer		dt_ista(MAXDATA)! [1..ndt] Station indices
	integer		dt_ic1(MAXDATA)	! [1..ndt] Event indices
	integer		dt_ic2(MAXDATA)	! [1..ndt] Event indices
	integer		src_cusp(MAXEVE)! [1..nsrc] Event keys
	real		src_t(MAXEVE)	! [1..nsrc] Event times
	real		tmp_ttp(stdim,MAXEVE)! [1.., 1..nsrc]
	real		tmp_tts(stdim,MAXEVE)! [1.., 1..nsrc]
	real		dt_cal(MAXDATA)	! [1..ndt] Theoretical time differences
	real		dt_res(MAXDATA)	! [1..ndt] Time-difference residuals

c	Local variables:
	integer		i
	real		tt1, tt2

      write(log,'("~ getting residual vector...")')

      if (nsrc.eq.1) then
c        Single source
         do i=1,ndt
            dt_res(i) = dt_dt(i)
         enddo
      else
c        Mulitple sources
	 tt1 = 0.0
	 tt2 = 0.0
         do i=1,ndt
            if (dt_idx(i).eq.1 .or. dt_idx(i).eq.3) then
c              P phase
               tt1 = tmp_ttp(dt_ista(i),dt_ic1(i)) -
     &              src_t(dt_ic1(i))/1000
               tt2 = tmp_ttp(dt_ista(i),dt_ic2(i)) -
     &              src_t(dt_ic2(i))/1000
            elseif (dt_idx(i).eq.2 .or. dt_idx(i).eq.4) then
c              S phase
               tt1 = tmp_tts(dt_ista(i),dt_ic1(i)) -
     &              src_t(dt_ic1(i))/1000
               tt2 = tmp_tts(dt_ista(i),dt_ic2(i)) -
     &              src_t(dt_ic2(i))/1000
            endif
            if (tt1.eq.0 .or. tt2.eq.0) then
               write(*,'("FATAL ERROR (theor tt). Please report to ",
     &                   "felixw@ldeo.columbia.edu")')
               stop
            endif
            dt_cal(i) = tt1 - tt2
            dt_res(i) = dt_dt(i) - dt_cal(i)
         enddo
      endif

      end !of subroutine dtres
