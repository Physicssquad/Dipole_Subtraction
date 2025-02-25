c     cone size
      function rjet(pa,pb)
      implicit double precision (a-h,o-z)
      dimension pa(0:3),pb(0:3)
      parameter (pi=3.14159265358979d0)

      ya=rapid(pa)
      yb=rapid(pb)
      deltay =ya-yb
      deltay2=deltay**2

      cosphi_a = pa(1)/eperp(pa)
      cosphi_b = pb(1)/eperp(pb)
      
      sinphi_a = pa(2)/eperp(pa)
      sinphi_b = pb(2)/eperp(pb)
      
      cos_deltaphi =  cosphi_a*cosphi_b + sinphi_a*sinphi_b
      cdp = cos_deltaphi
      
      if(cdp.lt.-1)then
c         write(*,*)'cos_deltaphi',cos_deltaphi
         cos_deltaphi=-1.0d0
      endif

      if(cdp.gt.1)then
c         write(*,*)'cos_deltaphi',cos_deltaphi
         cos_deltaphi=1.0d0
      endif      

      deltaphi     = DACOS(cos_deltaphi)
      deltaphi2    =deltaphi**2

      rjet = DSQRT( deltay2 + deltaphi2)
      
      return
      end
      
c     H(r)
      function hh(r)
      implicit double precision (a-h,o-z)      
      common/cone/ET_iso,r0,rgg
      common/nviso/niso

      n=niso
      hh = (1-dcos(r))/(1-dcos(r0))
      hh = hh**n   ! n =2 (Yuan): n=1 (Frixione)
      hh = ET_iso*hh

      return
      end

c     rapidity
      function rapid(p)
      implicit double precision (a-h,o-z)
      dimension p(0:3)
      
      xnum =p(0)+p(3)
      xden =p(0)-p(3)
      
      if(xnum.le.0d0)then
         rapid = -20.d0
      elseif(xden.le.0.d0)then
         rapid = +20.d0
      elseif(xnum.gt.0.d0 .and. xden.gt.0.d0)then
         rapid = 0.5d0*dlog(xnum/xden)
      endif
      
      return
      end

c     pt
      function eperp(p)
      implicit double precision (a-h,o-z)
      dimension p(0:3)
      eperp=dsqrt(p(1)*p(1)+p(2)*p(2))
      return
      end


