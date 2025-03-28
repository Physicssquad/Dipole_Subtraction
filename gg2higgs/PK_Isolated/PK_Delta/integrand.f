c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Delta terms]
      function flo2_PKDel(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
     .         ,f1(-6:6),f2(-6:6),xl(15)
     .         ,p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .         ,p(0:3,1:4),Born(1:2)
     .         ,SumP(1:2),SumK(1:2)
      double precision  Kterm,Pterm
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+9)
      common/energy/s
      common/scales/xmuf,xmur
      common/usedalpha/AL,ge
      common/distribution/xq
      common/bin_size/eps
      common/amass/am1,am2,amH,am4,am5
      character*10 P_,Kb_,Kt_
      character*10 gg_ 
      character*10 Plus_,Regular_,Delta_ 
      external PK
 
      P_       ='P'
      Kb_      ='Kb'
      Kt_      ='Kt'
      gg_      ='gg'
      Plus_    ='Plus'
      Regular_ ='Regular'
      Delta_   ='Delta'

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

c... page.79(eq:10.24/10.25)
      Tb_dot_Ta_by_Ta2 = -1d0
      Tfactors = Tb_dot_Ta_by_Ta2
      Tfactors = Tb_dot_Ta_by_Ta2
      Pterm =ALP*PK(P_,gg_,Delta_,x)*Tfactors*dlog(xmuf2/s12)
      Kterm =ALP*PK(Kb_,gg_,Delta_,x)-ALP*Tfactors*PK(Kt_,gg_,Delta_,x)

      SumPK = 2d0*coef*(Pterm + Kterm)
c... combination of only delta contributions`
         sig = xl(2)*SumPK !twice for sum over two legs

         pi_1 = PI/amH  
         flux = 2d0*sp
         pf_factor = 2d0*amH/xa/S

        xnorm= hbarc2*xajac*pi_1/flux

        flo2_PKDel  = xnorm*sig*pf_factor

      return
      end
c------------------------------------------------------------------
!c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Delta terms]
!      function flo2_PKDel2(yy,vwgt)
!      implicit double precision (a-h,o-z)
!      double precision  Kterm,Pterm
!      dimension yy(10)
!     .         ,f1(-6:6),f2(-6:6),xl(15)
!     .         ,p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
!     .         ,p(0:3,1:4),Born(1:2)
!     .         ,SumPK(1:2),SumK(1:2)
!      parameter (pi=3.14159265358979d0)
!      parameter (hbarc2=389.3856741D+9)
!      common/energy/s
!      common/scales/xmuf,xmur
!      common/usedalpha/AL,ge
!      common/distribution/xq
!      common/bin_size/eps
!      common/amass/am1,am2,amH,am4,am5
!
!      tau = amH**2/S
!      xajac = 1d0 - tau
!      xa = tau + xajac*yy(1)
!      xb = tau/xa
!
!      delta = 0d0
!      xmin = 0d0
!
!      xmax = 1.0d0 - delta
!      xjac4 = (xmax - xmin)
!      x = xjac4*yy(1) + xmin
!
!      xamin = tau/x
!      xamax = 1d0
!      xajac =  (xamax - xamin)
!      xa = xamin + xajac*yy(2)
!      xb = tau/x/xa
!
!      flo2_PKDel2 = 0d0
!!      if (xb .ge. 1d0 ) goto 151
!c... till here linear mapping is added
!!
!      Pterm = 0d0
!      Kterm = 0d0
!      sp = xa*xb*S
!      rsp = dsqrt(sp)
!      do k=1,2
!
!        if ( k .eq. 1) call kinvar1(xa*x,xb,p1,p2,p3)
!        if ( k .eq. 2) call kinvar1(xa,xb*x,p1,p2,p3)
!
!!      sp = 2d0*dot(p1,p2)
!!      rsp = dsqrt(sp)
!c... same reasoning as plusA. refer there.
!
!      AL = alphasPDF(xmur)
!      ALP = AL/2d0/Pi
!
!      Tb_dot_Ta_by_Ta2 = -1d0
!      xmuf2 = xmuf*xmuf
!      s12   = 2d0*dot(p1,p2)
!      coef = Born_gg2h_(0,AL,p1,p2,p3)
!
!c... page.79(eq:10.24/10.25)
!      Pterm =ALP*PggD(x)*Tb_dot_Ta_by_Ta2*dlog(xmuf2/s12)
!      Kterm =ALP*AKbarD_gg(x) - ALP*Tb_dot_Ta_by_Ta2*AKtilD_gg(x)
!      SumPK(k) = coef*(Pterm + Kterm)
!      enddo     
!
!      call pdf(xa,xmuf,f1)
!      call pdf(xb,xmuf,f2)
!      call setlum(f1,f2,xl)
!  
!      sig = xl(2)*(SumPK(1) + SumPK(2))
!
!      pi_1 = PI/amH 
!      flux = 2d0*sp
!      xnorm = hbarc2*xajac*xjac4*pi_1/flux
!      pf_factor = 2d0*amH/xa/x/S
!
!      flo2_PKDel2 = xnorm*sig*pf_factor
!
!151   return
!      end
!
