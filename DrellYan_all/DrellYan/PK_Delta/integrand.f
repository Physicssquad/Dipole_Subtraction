c
cc~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Delta terms]
c      function flo2_PKDel(yy,vwgt)
c      implicit double precision (a-h,o-z)
c      dimension yy(10)
c     .         ,f1(-6:6),f2(-6:6),xl(15)
c     .         ,p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
c     .         ,p(0:3,1:4),Born(1:2)
c     .         ,SumP(1:2),SumK(1:2)
c      parameter (pi=3.14159265358979d0)
c      parameter (hbarc2=389.3856741D+6)
c      common/energy/s
c      common/factscale/xmuf
c      common/usedalpha/AL,ge
c      common/distribution/xq
c
c      xa     = yy(1)
c      xb     = yy(2)
c
c      sp     = xa*xb*s
c      rsp    = dsqrt(sp)
c      xjac = 2.0d0
c      xnorm=hbarc2
c
c      eps = 0.5d0
c      xlow = xq - eps
c      xhigh = xq + eps
c
c        call kinvar2(yy,xinvmass,p1,p2,p3,p4)
c        call p1d_to_p2d_4(p1,p2,p3,p4,p)
c
c        pin  = 0.5d0*rsp
c        flux_1 = 4.0d0*pin*rsp
c
c        scale2 = 2.0d0*dot(p1,p2)
c        scale = dsqrt(scale2)
c
c      flo2_PKDel  = 0.0d0
c      sig1 = 0d0
c
c
c        if (scale .ge. xlow .and. scale .le. xhigh)  then
c
c            xmuf = scale
c            xmur = scale
c            AL = alphasPDF(xmur)
c            ALP = AL/2d0/Pi
c
c            call getPKDel(1.0d0,xmuf,p,p1,p2,SumDel)
c
c            call pdf(xa,xmuf,f1)
c            call pdf(xb,xmuf,f2)
c            call setlum(f1,f2,xl)
c
c            sig1 = ALP*xl(1)* SumDel  !  [qq lum]
c
c            azmth = 2.0d0*pi
c            pf = 0.5d0*rsp
c            ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth
c
c            wgt1 = sig1/flux_1*ps2*xjac*vwgt
c            PKDel = xnorm*wgt1/vwgt/2d0/eps
c
c         endif
c         flo2_PKDel = PKDel
c
c      return
c      end
c
c---------------------------------------------------------------------
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Delta terms]
      function flo2_PKDel(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
     .         ,f1(-6:6),f2(-6:6),xl(15)
     .         ,p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .         ,p(0:3,1:4),Born(1:2)
     .         ,SumP(1:2),SumK(1:2)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+6)
      common/energy/s
      common/factscale/xmuf
      common/usedalpha/AL,ge
      common/distribution/xq

c      xa     = yy(1)
c      xb     = yy(2)

      TAUH = xq*xq/s

      xa = (1D0-TAUH)*yy(1)    + TAUH
c      xb = (1D0-TAUH/xa)*yy(2) + TAUH/xa
      xb = TAUH/xa

      xc = yy(2)

      sp     = xa*xb*s
      rsp    = dsqrt(sp)
      xjac = 2.0d0*(1.0d0-TAUH)
      xnorm=hbarc2*2.0d0*xq/s/xa

        call kinvar2_PK(xa,xb,xc,xinvmass,p1,p2,p3,p4)

        call p1d_to_p2d_4(p1,p2,p3,p4,p)

        pin  = 0.5d0*rsp
        flux_1 = 4.0d0*pin*rsp

        scale = xq

      flo2_PKDel  = 0.0d0
      sig1 = 0d0


            xmuf = scale
            xmur = scale
            AL = alphasPDF(xmur)
            ALP = AL/2d0/Pi

            call getPKDel(1.0d0,xmuf,p,p1,p2,SumDel)

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call setlum(f1,f2,xl)

            sig1 = ALP * xl(1) * SumDel  !  [qq lum]


            azmth = 2.0d0*pi
            pf = 0.5d0*rsp
            ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth

            wgt1 = sig1/flux_1*ps2*xjac*vwgt
            PKDel = xnorm*wgt1/vwgt

         flo2_PKDel = PKDel

      return
      end

c---------------------------------------------------------------------
