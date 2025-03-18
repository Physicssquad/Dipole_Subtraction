c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Regular Terms]
      function flo2_PKReg(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
     &  ,f1(-6:6),f2(-6:6),xl(15)
     &  ,p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     &  ,p(0:3,1:4),Born(1:2),AllReg(1:2)
     &  ,SumP(1:2),SumK(1:2)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+9)
      common/energy/s
c      common/factscale/xmuf
      common/usedalpha/AL,ge
      common/distribution/xq
      common/bin_size/eps
      common/amass/am1,am2,amH,am4,am5
      common/scales/xmuf,xmur

       tau = amH**2/S
c... for regular terms I don't need to use the upper cut as it is a flat integral.
      delta = 0d0  
      xmin = 0d0 

      xmax = 1.0d0 - delta
      xjac4 = (xmax - xmin)
      x = xjac4*yy(1) + xmin

!      xamin = tau/x
      xamin = 0d0 
      xamax = 1d0
      xajac =  (xamax - xamin)
      xa = xamin + xajac*yy(2)
      xb = tau/x/xa
c... till here linear mapping is added
!
      AllReg(1) = 0d0
      AllReg(2) = 0d0
      PKReg = 0.0d0
c      sp = xa*xb*S
c      rsp = dsqrt(sp)

      do k=1,2

        if ( k .eq. 1) call kinvar1(xa*x,xb,p1,p2,p3)
        if ( k .eq. 2) call kinvar1(xa,xb*x,p1,p2,p3)

      sp = 2d0*dot(p1,p2)
      rsp = dsqrt(sp)
c... same reasoning as plusA. refer there.

      AL = alphasPDF(xmur)
      ALP = AL/2d0/Pi

      xmuf2 = xmuf*xmuf
      s12   = 2d0*dot(p1,p2)
      Preg = PggReg(x)*(-1.0d0)*dlog(xmuf2/s12)    ! here x=1
      SumReg =  PReg + AKbarReg_gg(x) + AKtilReg_gg(x)      

      coef = Born_gg2h_(0,AL,p1,p2,p3)

      call pdf(xa,xmuf,f1)
      call pdf(xb,xmuf,f2)
      call setlum(f1,f2,xl)

      sig1 = xl(2)* SumReg  

      sig = Alp*sig1*coef

      pi_1 = PI/amH 
      flux = 2d0*sp
      xnorm = hbarc2*pi_1/flux

      PKReg = xnorm*xajac*xjac4*sig* 2d0*amH/xa/x/S

      AllReg(k) = PKReg  
      enddo
      flo2_PKReg = AllReg(1) + AllReg(2)
151   return
      end
