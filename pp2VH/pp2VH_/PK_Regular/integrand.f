
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Regular Terms]
      function flo2_PKReg(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .          ,p(0:3,1:4),Born(1:2),AllReg(1:2)
      dimension SumP(1:2),SumK(1:2)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+6)
      common/energy/s
      common/factscale/xmuf
      common/usedalpha/AL,ge
      common/distribution/xq

      xa     = yy(1)
      xb     = yy(2)
      xc     = yy(3)
      x      = yy(4)

      xjac   = 2d0
      xnorm=hbarc2

      eps = 0.5d0
      xlow = xq - eps
      xhigh = xq + eps

      AllReg(1) = 0d0
      AllReg(2) = 0d0

      do k=1,2

      if (k .eq. 1) call kinvar2_PK(x*xa,xb,xc,Qmass,p1,p2,p3,p4)
      if (k .eq. 2) call kinvar2_PK(xa,x*xb,xc,Qmass,p1,p2,p3,p4)

        sp = 2.0d0*dot(p1,p2)
        rsp = dsqrt(sp)
        scalex = Qmass
        pin  = 0.5d0*rsp
        flux = 4.0d0*pin*rsp

        PKReg = 0.0d0

        if (scalex .ge. xlow .and. scalex .le. xhigh) then

            xmuf = scalex
            xmur = scalex
             AL = alphasPDF(xmur)
            ALP = AL/2d0/Pi

            call getPKReg(x,xmuf,p1,p2,p3,p4,SumReg)

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call setlum(f1,f2,xl)

            sig1 = xl(1)* SumReg  !  [qq lum]

            sig = Alp*sig1

            azmth = 2.0d0*pi
            pf = 0.5d0*rsp
            ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth

               wgt_x = sig/flux*ps2*xjac*vwgt
           AllReg(k) = xnorm*wgt_x/vwgt/2d0/eps

         endif
         enddo

         flo2_PKReg = AllReg(1) + AllReg(2)

      return
      end
c---------------------------------------------------------------------
