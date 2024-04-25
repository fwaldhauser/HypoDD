c-------------------------------------------------------------------------
      subroutine delaz(a1lat,a1lon,a2lat,a2lon,del,dist,az)
c
c     by Bill Ellsworth
c
c        computes distance and azimuth from a1 to a2
c        a1 and a2 are in decimal degrees and n-e coordinates
c        del -- delta in degrees
c        dist -- distance in km
c        az -- azimuth from a to b clockwise from north in degrees

c     changes by Felix Waldhauser (fw)
c
      real*8 pi2,rad,flat
      real*8 alatr,alonr,blatr,blonr
      real*8 tana,geoa,acol,tanb,geob,bcol
      real*8 diflon,cosdel,delr,top,den,azr,colat,radius
cfw      real*8 dtan,datan,dsin,dcos,darcos,dcotan,datan2
      real*8 dtan,datan,dsin,dcos,datan2
      data pi2/1.570796d0/
      data rad/1.745329d-02/
      data flat/.993231d0/

c-----convert to radians
      alatr=a1lat*rad
      alonr=a1lon*rad
      blatr=a2lat*rad
      blonr=a2lon*rad
c-----convert latitudes to geocentric colatitudes
      tana=flat*dtan(alatr)
      geoa=datan(tana)
      acol=pi2-geoa
      tanb=flat*dtan(blatr)
      geob=datan(tanb)
      bcol=pi2-geob
c-----calcuate delta
      diflon=blonr-alonr
      cosdel=dsin(acol)*dsin(bcol)*dcos(diflon)+dcos(acol)*dcos(bcol)
cfw      delr=darcos(cosdel)
      delr=dacos(cosdel)
c-----calcuate azimuth from a to b
      top=dsin(diflon)
cfw      den=dsin(acol)*dcotan(bcol)-dcos(acol)*dcos(diflon)
      den=dsin(acol)*(1/dtan(bcol))-dcos(acol)*dcos(diflon)
      azr=datan2(top,den)
c-----convert to degrees
      del=delr/rad
      az=azr/rad
      if (az.lt.0.0) az=360+az
c-----compute distance in kilometers
      colat=pi2-(alatr+blatr)/2
      radius=6371.227*(1.0+3.37853d-3*(1/3-((dcos(colat))**2)))
cfw      radius=6378.140*(1.0+3.37853d-3*(1/3-((dcos(colat))**2)))
      dist=delr*radius
      return
      end
