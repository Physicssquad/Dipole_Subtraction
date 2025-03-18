c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Delta terms]
      function flo2_PKDel(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
     .         ,f1(-6:6),f2(-6:6),xl(15)
     .         ,p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .         ,p(0:3,1:4),Born(1:2)
     .         ,SumP(1:2),SumK(1:2)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+9)
      common/energy/s
      common/scales/xmuf,xmur
      common/usedalpha/AL,ge
      common/distribution/xq
      common/bin_size/eps
      common/amass/am1,am2,amH,am4,am5

        tau = amH**2/S
         xajac = 1d0 - tau
         xa = tau + xajac*yy(1)
         xb = tau/xa

         sp = xa*xb*S
        rsp = dsqrt(sp)

        call kinvar1(xa,xb,p1,p2,p3)
        xmuf2=xmuf**2
         AL = alphasPDF(xmur)
         ALP= AL/2d0/PI

        call pdf(xa,xmuf,f1)
        call pdf(xb,xmuf,f2)
        call setlum(f1,f2,xl)

      s12 = 2d0*dot(p1,p2)
      coef = Born_gg2h_(0,AL,p1,p2,p3)
      Pdel = PggD(x)*(-1.0d0)*dlog(xmuf2/s12)    ! here x=1
      SumDel = Pdel - AKbarD_gg(x) + AKtilD_gg(x)
c... there is an overall negative in Kbar_gg

         sig = 2d0*xl(2)*ALP*SumDel*coef !twice for sum over two legs

         pi_1 = PI/amH  
         flux = 2d0*sp
         pf_factor = 2d0*amH/xa/S

        xnorm= hbarc2*xajac*pi_1/flux

        flo2_PKDel  = xnorm*sig*pf_factor

      return
      end
c------------------------------------------------------------------
