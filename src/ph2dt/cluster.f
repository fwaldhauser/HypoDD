
      subroutine cluster(log, mev, mcl, npair, minobs,
     & ic1, ic2, in, icusp, acl, acli, apair_n,
     & clust, noclust, nclust)

	implicit none

c	Parameters:
	integer		log
	integer		mev
	integer		mcl
	integer		npair
	integer		minobs
	integer		ic1(*)		! (1..mpair)
	integer		ic2(*)		! (1..mpair)
	integer		in(*)		! (1..mpair)
	integer		icusp(*)	! (1..mev)
	integer		acl(*)		! (1..mev)
	integer		acli(*)		! (1..mev)
	character	apair_n((mev**2-mev)/2)*1
	integer		clust(mcl,*)	! (a..mcl, 1..mev)
	integer		noclust(*)	! (1..mev)
	integer		nclust

c	Local variables:
	integer		iev
	integer		nev
	integer		i
	integer		j
	integer		k
	integer		kk
	integer		n
	integer		ii
	integer		nn
	integer		ia
	integer		ifindi

      write(*,'("clustering ...  ")')


      iev= 0
      do i=1,npair
        do j=1,iev
          if(ic1(i).eq.icusp(j)) goto 10
        enddo
        iev= iev+1
        icusp(iev)= ic1(i)
        do j=1,iev
          if(ic2(i).eq.icusp(j)) goto 10
        enddo
        iev= iev+1
        icusp(iev)= ic2(i)
10      continue
      enddo
      nev= iev
      call sorti(nev,icusp)

      do i=1,npair
         j= ifindi(nev,icusp,ic1(i))
         k= ifindi(nev,icusp,ic2(i))
         if(k.gt.j) then   !map into lower triangle
            kk= k
            k= j
            j= kk
         endif
         if(in(i).ge.minobs) then
            apair_n(((j-1)**2-(j-1))/2+k)= '1'
         else
            apair_n(((j-1)**2-(j-1))/2+k)= '0'
         endif

      enddo

      do i=1,nev
         acl(i)= nev+1
      enddo
      k= 0
      n= 0
      do i=2,nev
         do j=1,i-1
            k=k+1
            if(apair_n(k).eq.'1')then
                if(acl(i).lt.acl(j)) then
                   if(acl(j).eq.nev+1) then
                      acl(j)= acl(i)
                   else
                      ia= acl(j)
                      do ii=1,i
                        if(acl(ii).eq.ia) acl(ii)= acl(i)
                      enddo
                   endif
                elseif(acl(j).lt.acl(i)) then
                   if(acl(i).eq.nev+1) then
                      acl(i)= acl(j)
                   else
                      ia= acl(i)
                      do ii=1,i
                        if(acl(ii).eq.ia) acl(ii)= acl(j)
                      enddo
                   endif
                elseif(acl(i).eq.nev+1) then
                   n= n+1
                   acl(i)= n
                   acl(j)= n
                endif
            endif
         enddo
      enddo

      call indexxi(nev,acl,acli)
      n= 1
      nn= 2
      if(acl(acli(1)).eq.nev+1) then
         i= 1
         goto 300
      else
         clust(n,nn)= icusp(acli(1))
      endif
      do i=2,nev
         if(acl(acli(i)).gt.acl(acli(i-1))) then
            clust(n,1)= nn-1
            n= n+1
            nn= 1
         endif
         if(acl(acli(i)).eq.nev+1) goto 300
         nn= nn+1
         clust(n,nn)= icusp(acli(i))
      enddo
      clust(n,1)= nn-1
      nclust= n
      goto 310
300   nclust= n-1
      do j=i,nev
        noclust(j-i+2)= icusp(acli(j))
      enddo
      noclust(1)= nev-i+1
310   continue
      if(nclust.ge.mcl) stop'>>> Increase MCL in ph2dt.inc.'

c--- sort - biggest cluster first
      if(nclust.gt.1) then
        do i=1,nclust-1
           do j= i+1,nclust
             if(clust(i,1).le.clust(j,1)) then
                do k=1,clust(i,1)+1
                     clust(mcl,k)= clust(i,k)
                enddo
                do k=1,clust(j,1)+1
                     clust(i,k)= clust(j,k)
                enddo
                do k=1,clust(mcl,1)+1
                     clust(j,k)= clust(mcl,k)
                enddo
             endif
           enddo
        enddo
      endif

      k= 0
      do i=1,nclust
        k= k + clust(i,1)
      enddo

      write(*,'(" > clustered events = ",i6)')k
      write(*,'(" > isolated events = ",i6)')
     & noclust(1)
      write(*,'(" > clusters = ",i6)')nclust
      write(*,'(" > largest cluster = ",i6," events")')clust(1,1)

      write(log,'("> Clustered events: ",i6)')k
      write(log,'("> Isolated events: ",i6,/,8(i8,2x))')
     & noclust(1),(noclust(j),j=2,noclust(1)+1)
      write(log,'(/,"> clusters:",i6,"; MINLNKS = ",i4)')nclust,
     & minobs
      k= 0
      do i=1,nclust
         write(log,'("> Cluster",i6,": ",i6," events",/,8(i8,2x))')
     & i,clust(i,1),(clust(i,j),j=2,clust(i,1)+1)
         k= k+clust(i,1)
      enddo
      write(log,*)

      end !of subroutine cluster1
