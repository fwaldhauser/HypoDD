	subroutine cluster1(log, nev, ndt,
     &	idata, minobs_cc, minobs_ct,
     &	dt_c1, dt_c2, ev_cusp,
     &	clust, noclust, nclust)

	implicit none

	include'hypoDD.inc'

c	Parameters:
	integer		log		! Log-file identifier
	integer		nev		! No. of events
	integer		ndt		! No. of data
	integer		idata
	integer		minobs_cc	! Min. obs./pair for ccor. data
	integer		minobs_ct	! Min. obs./pair for cat. data
	integer		dt_c1(MAXDATA)	! [1..ndt] Event keys
	integer		dt_c2(MAXDATA)	! [1..ndt] Event keys
	integer		ev_cusp(MAXEVE)	! [1..nev] Event keys
	integer		clust(MAXCL,MAXEVE) ! [1..nclust,1..] Event keys
					! clust[i,1] = Size of ith cluster
	integer		noclust(MAXEVE)	! [1..MAXEVE] Event keys
					! noclust[1] = No. of keys
	integer		nclust		! No. of cllusters

c	Local variables:
	integer		acli(MAXEVE)
	integer		acl(MAXEVE)
	integer		apair_n((MAXEVE*(MAXEVE-1))/2)
	integer		i
	integer		ia
	integer		icusp(MAXEVE)	! [1..nev] Event keys
	integer		ifindi
	integer		ii
	integer		j
	integer		k
	integer		kk
	integer		n
	integer		nn

      write(*,'("clustering ...  ")')
      write(log,'(/,"~ clustering ...  ")')

c     Set up event-pair arrays
      k = 0
      do i=2,nev
         do j=1,i-1		!lower triangle of matrix
            k = k+1
            apair_n(k) = 0
         enddo
      enddo

      do i=1,nev
         icusp(i) = ev_cusp(i)
      enddo
      call sorti(nev, icusp)

      do i=1, ndt
         j = ifindi(nev, icusp, dt_c1(i))
         k = ifindi(nev, icusp, dt_c2(i))
         if (k.gt.j) then
c           Map into lower triangle
            kk = k
            k = j
            j = kk
         endif
         apair_n(((j-1)**2-(j-1))/2+k)=
     &          apair_n(((j-1)**2-(j-1))/2+k) + 1
      enddo
      if (idata.eq.0 .or. idata.eq.1) minobs_ct = 0
      if (idata.eq.0 .or. idata.eq.2) minobs_cc = 0

c     Initialize array acl to store cluster index for each event
      do i=1,nev
         acl(i) = nev+1
      enddo
      k = 0
      n = 0
      do i=2,nev
         do j=1,i-1
            k = k+1
            if (apair_n(k).ge.(minobs_cc+minobs_ct)) then
               if (acl(i).lt.acl(j)) then
                  if (acl(j).eq.nev+1) then
                     acl(j) = acl(i)
                  else
                     ia = acl(j)
                     do ii=1,i
                       if (acl(ii).eq.ia) acl(ii) = acl(i)
                     enddo
                  endif
               elseif (acl(j).lt.acl(i)) then
                  if (acl(i).eq.nev+1) then
                     acl(i) = acl(j)
                  else
                     ia = acl(i)
                     do ii=1,i
                       if (acl(ii).eq.ia) acl(ii) = acl(j)
                     enddo
                  endif
               elseif (acl(i).eq.nev+1) then
                  n = n+1
                  acl(i) = n
                  acl(j) = n
               endif
            endif
         enddo
      enddo

c     Store event keys in cluster matrix clust[]
      call indexxi(nev,acl,acli)
      n = 1
      nn = 2
      if (acl(acli(1)).eq.nev+1) then
c        No events clustered
         i = 1
         goto 300
      else
         clust(n,nn) = icusp(acli(1))
      endif
      do i=2,nev
         if (acl(acli(i)).gt.acl(acli(i-1))) then
            clust(n,1) = nn-1
            n = n+1
            nn = 1
         endif
         if (acl(acli(i)).eq.nev+1) goto 300	! events not clustered
         nn = nn+1
         clust(n,nn) = icusp(acli(i))
      enddo
      clust(n,1) = nn-1
      nclust = n
      noclust(1) = 0
      goto 310
300   nclust = n-1
      do j=i,nev
         noclust(j-i+2) = icusp(acli(j))
      enddo
      noclust(1) = nev-i+1
310   continue
      if (nclust.ge.MAXCL) stop'>>> Increase MAXCL in hypoDD.inc.'

c     Sort - biggest cluster first
      if (nclust.gt.1) then
         do i=1,nclust-1
            do j=i+1,nclust
               if (clust(i,1).le.clust(j,1)) then
                  do k=1,clust(i,1)+1
                     clust(MAXCL,k) = clust(i,k)
                  enddo
                  do k=1,clust(j,1)+1
                     clust(i,k) = clust(j,k)
                  enddo
                  do k=1,clust(MAXCL,1)+1
                     clust(j,k) = clust(MAXCL,k)
                  enddo
               endif
            enddo
         enddo
      endif

      k = 0
      do i=1,nclust
         k = k + clust(i,1)
      enddo

      write(*,'("Clustered events: ",i5)') k
      write(*,'("Isolated events: ",i5)') noclust(1)
      write(*,'("# clusters:",i5)')nclust
      k = 0
      do i=1,nclust
         write(*,'("Cluster",i4,": ",i5," events")') i, clust(i,1)
         k = k+clust(i,1)
      enddo

      write(log,'("# clustered events =",i5)')k
cfw      write(log,'("# isolated events =",i5,/,8(i8))')
      write(log,'("# isolated events =",i5,/,8(i9))')
     & noclust(1),(noclust(j),j=2,noclust(1)+1)
      write(log,'("# clusters =",i5,"  for min. number of links "
     & "set to ",i5)')nclust,minobs_ct+minobs_cc
      k = 0
      do i=1,nclust
         write(log,'("Cluster",i5,": ",i5," events",/,8(i9,2x))')
     & i,clust(i,1),(clust(i,j),j=2,clust(i,1)+1)
         k = k+clust(i,1)
      enddo
      write(log,*)
      if (nclust.eq.0) stop

      end !of subroutine cluster1
