c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function flo1_LO(yy,vwgt)
      use openloops
      implicit double precision (a-h,o-z)
      dimension yy(2)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .          ,p(0:3,1:4),c(1:2),Born(1:2),coef(2,-2:0),SumI(-2:0)
     .          ,p_ex(0:3,1:3)
      parameter (pi=3.14159265358979d0)
c      parameter (hbarc2=0.3894d9)      ! in Pb
c      parameter (hbarc2=389.3856741D+9)
      parameter (hbarc2=0.3894d12)    ! in Fb
      common/energy/S
      common/renor_scale/scale
      common/usedalpha/AL,ge
      common/amass/am1,am2,amH,am4,am5
      common/caller/icall,id
      external Born_gg2H
       
c... this is for the integration variable reduction
coooooooooooooooooooooooooooooooooooooooooooooooooo
        tau = amH**2/S
         xajac = 1d0 - tau
         xa = tau + xajac*yy(1)
         xb = tau/xa 
coooooooooooooooooooooooooooooooooooooooooooooooooo

         sp = xa*xb*S
        rsp = dsqrt(sp)
       amH2 = amH**2

        call kinvar1(xa,xb,p1,p2,p3)
        s12 = 2d0*dot(p1,p2)

        xmuf= amH/2d0
        xmur= xmuf
        xmu2=xmuf**2

        AL = alphasPDF(xmur)
        call set_parameter("alpha_s",AL)

        call pdf(xa,xmuf,f1)
        call pdf(xb,xmuf,f2)
        call setlum(f1,f2,xl)

         pi_1 = PI/amH
         flux = 2d0*sp 

c         sig =xl(2)*Born_gg2h(0,p1,p2,p3)
         sig =xl(2)*Born_gg2h_(0,AL,p1,p2,p3)

         pf_factor = 2d0*amH/xa/S
         xnorm=hbarc2*xajac*pi_1/flux

        flo1_LO  = xnorm * sig * pf_factor

151   return
      end
c---------------------------------------------------------------------
