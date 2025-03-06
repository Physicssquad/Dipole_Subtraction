
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

      xa     = yy(1)
      xb     = yy(2)

      sp     = xa*xb*s
      rsp    = dsqrt(sp)
      xjac = 2.0d0
      xnorm=hbarc2

      eps = 0.5d0
      xlow = xq - eps
      xhigh = xq + eps

        call kinvar2(yy,xinvmass,p1,p2,p3,p4)
        call p1d_to_p2d_4(p1,p2,p3,p4,p)

        pin  = 0.5d0*rsp
        flux_1 = 4.0d0*pin*rsp

        scale2 = 2.0d0*dot(p1,p2)
        scale = dsqrt(scale2)

      flo2_PKDel  = 0.0d0
      sig1 = 0d0


        if (scale .ge. xlow .and. scale .le. xhigh)  then

            xmuf = scale
            xmur = scale
            AL = alphasPDF(xmur)
            ALP = AL/2d0/Pi

            call getPKDel(1.0d0,xmuf,p,p1,p2,SumDel)

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call setlum(f1,f2,xl)

            sig1 = xl(4)* SumDel  !  [gg lum]

            sig1 = Alp*sig1

            azmth = 2.0d0*pi
            pf = 0.5d0*rsp
            ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth

            wgt1 = sig1/flux_1*ps2*xjac*vwgt
            PKDel = xnorm*wgt1/vwgt/2d0/eps

         endif
         flo2_PKDel = PKDel

      return
      end

c---------------------------------------------------------------------
