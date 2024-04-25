      subroutine setorg_3d(lat,lon,rot)

      implicit none

      real	lat
      real	lon
      real	rot

      doubleprecision xltkm, xlnkm, rota
      doubleprecision snr, csr, dxlt, dxln
      common/shrotd/ xltkm, xlnkm, rota, dxlt, dxln, snr, csr

      doubleprecision dlt1,dlt2,drad,drlt
      doubleprecision re, ell
      doubleprecision del, r, bc
      parameter ( re=6378.163, ell=298.26)
      parameter (drad=1.7453292d-2)
      parameter (drlt=9.9330647d-1)


      rota = dble(rot)*drad
      dxlt = 60.d0 * dble(lat)
      dxln  = 60.d0 * dble(lon)

      dlt1=datan(drlt*dtan(dxlt*drad/60.d0))
      dlt2=datan(drlt*dtan((dxlt+1.d0)*drad/60.d0))

      del = dlt2-dlt1
      r = re*(1.0 - (dsin(dlt1)**2)/ell)
      xltkm=r*del

      del = dacos(1.0d0-(1.0d0-dcos(drad/60.d0))*dcos(dlt1)**2)
      bc  = r*del
      xlnkm = bc/dcos(dlt1)

      snr=sin(rota)
      csr=cos(rota)

      return
      end
c
c --------------------------------------------------------------
c

      subroutine disto(lon,lat,x,y)

      implicit none

      real	lon,lat
      real	x,y

      doubleprecision drad,drlt
      parameter (drad=1.7453292d-02)
      parameter (drlt=9.9330647d-01)

      doubleprecision xltkm, xlnkm, rota
      doubleprecision snr, csr, dxlt, dxln
      common/shrotd/ xltkm, xlnkm, rota, dxlt, dxln, snr, csr

      doubleprecision plt, pln
      doubleprecision tmpx, tmpy
      doubleprecision dxlt1


      plt = 60.d0*lat
      pln = 60.d0*lon

      tmpx = pln - dxln
      tmpy = plt - dxlt

      dxlt1 = datan(drlt*dtan(drad*(plt+dxlt)/120.))
      tmpx = tmpx*xlnkm*dcos(dxlt1)
      tmpy = tmpy*xltkm

      if(rota.eq.0.0) then
	x = sngl(tmpx)
	y = sngl(tmpy)
      else
	x = sngl(csr*tmpx - snr*tmpy)
	y = sngl(csr*tmpy + snr*tmpx)
      endif

      return
      end



