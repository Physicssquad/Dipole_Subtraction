      function flo2_Plus(yy,vwgt)
      implicit double precision (a-h,o-z)
      double precision Pterm, Kterm
      dimension yy(10),PKt(1:2)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension xp1(0:3),xp2(0:3)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+9)
      integer k,iplus
      character*50 name

      character*10 P_,Kb_,Kt_
      character*10 gg_ 
      character*10 Plus_,Regular_,Delta_ 

      common/energy/s
      common/pdfname/name
      common/usedalpha/AL,ge
      common/distribution/xq
      common/amass/am1,am2,amH,am4,am5
      common/scales/xmuf,xmur
      external PK

      P_       ='P'
      Kb_      ='Kb'
      Kt_      ='Kt'
      gg_      ='gg'
      Plus_    ='Plus'
      Regular_ ='Regular'
      Delta_   ='Delta'

c... I am going to convert this code for the gg2higgs case
c... there will be huge modifications that will be don.
c... Here s = ecm*ecm

       tau = amH**2/S

      delta = 1d-9
      xmin =0d0 ! technically xmin should run from 0 to 1. 
!      xmin = tau
      xmax = 1.0d0 - delta
      xjac4 = (xmax - xmin)
      x = xjac4*yy(1) + xmin
c... for the fixed value of x. I will generate the xa and xb using the 
c... constraint relation.
c... x will always be there as an overall factor.
c~~~~~[ _______________ ]~~~~~C
!      xamin = tau/x
      xamin = 0d0
      xamax = 1d0
      xajac =  (xamax - xamin)
      xa = xamin + xajac*yy(2)
!      xb = tau/x/xa  !its better to define xb for the specific kin type.
c... here x* xa*xb*s = amH**2 is the kinametics.
c... xb is generated for every case with Higgs mass constraint.

      sig1  = 0d0
      sig3  = 0d0
      PKplus_1 = 0d0
      PKplus_x = 0d0

      do k = 1,2
      do iplus = 0,1

!#####################################
      if (iplus .eq. 1) then
         xb = tau/x/xa
         if (k .eq. 1) call kinvar1(x*xa,xb,p1,p2,p3)
         if (k .eq. 2) call kinvar1(xa,xb*x,p1,p2,p3)
         sp = 2d0*dot(p1,p2)
         rsp = dsqrt(sp)

      elseif (iplus .eq. 0) then
         xb = tau/xa   !// here x=1 is set before PS_gen, however x goes 
         sp = xa*xb*S
         call kinvar1(xa,xb,p1,p2,p3)
         s12 = 2d0*dot(p1,p2)

!         sp = 2d0*dot(p1,p2)
!         rsp = dsqrt(sp)
!        rsp = dsqrt(sp)
c... indipendently into the [+] functions.
        endif
!#####################################

!         sp = xa*xb*S
!        rsp = dsqrt(sp)
            xmuf2 = xmuf*xmuf 
               AL = alphasPDF(xmur) 
         coef = Born_gg2h_(0,AL,p1,p2,p3)
              ALP = AL/2d0/Pi

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call setlum(f1,f2,xl)

!            Pplus = PggP(x)*(-1.0d0)*dlog(xmuf2/s12)
!          SumPlus = Pplus + AKtilP_gg(x) + AKbarP_gg(x)

      Tb_dot_Ta_by_Ta2 = -1d0
      Tfactors = Tb_dot_Ta_by_Ta2
      Pterm =ALP*PK(P_,gg_,Plus_,x)*Tfactors*dlog(xmuf2/s12)
      Kterm =ALP*PK(Kb_,gg_,Plus_,x) -
     &  ALP*Tfactors*PK(Kt_,gg_,Plus_,x)

      SumPlus =coef*(Pterm + Kterm)


c... Ktilde comes with the overall negative but due to colour charge one 
c... more negative comes making it positive in this case only.
             
            if (iplus .eq. 1) sig1 = xl(2)*SumPlus
            if (iplus .eq. 0) sig3 = xl(2)*SumPlus

           pi_1 = PI/amH
           flux = 2d0*sp
           xnorm = hbarc2*pi_1/flux
c... reset before proceding
          if (iplus .eq. 1 ) then
          PKplus_x = xnorm*xajac*xjac4*sig1*2d0*amH/xa/x/S
          elseif( iplus .eq. 0 ) then
          PKplus_1 = xnorm*xajac*xjac4*sig3*2d0*amH/xa/S
          endif

          enddo

            PKt(k) = PKplus_x - PKplus_1

        enddo

        flo2_Plus = ( PKt(1) + PKt(2) )
151      return
      end
c---------------------------------------------------------------------
      function xint_PlusA(yy,vwgt)
      use openloops
      implicit double precision (a-h,o-z)
      dimension yy(10)
     &  ,f1(-6:6),f2(-6:6),xl(15)
     &  ,p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     &  ,p(0:3,1:4),Born(1:2),PK_leg(1:2)
     &  ,SumP(1:2),SumK(1:2)
      double precision Pterm, Kterm
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+9)

      character*10 P_,Kb_,Kt_
      character*10 gg_ 
      character*10 Plus_,Regular_,Delta_ 

      common/energy/s
      common/usedalpha/AL,ge
      common/distribution/xq
      common/bin_size/eps
      common/amass/am1,am2,amH,am4,am5
      common/scales/xmuf,xmur
      external PK

      P_       ='P'
      Kb_      ='Kb'
      Kt_      ='Kt'
      gg_      ='gg'
      Plus_    ='Plus'
      Regular_ ='Regular'
      Delta_   ='Delta'

      tau = amH**2/S

c... Integration is singular at 1, so we are performing integration 
c... from 0 to 1-δ
      delta = 1d-2 
      xmin = 0d0 
      xmax = 1.0d0 - delta
      xjac4 = (xmax - xmin)
      x = xjac4*yy(1) + xmin

      xamin = tau/x
!      xamin = 0d0 
      xamax = 1d0
      xajac =  (xamax - xamin)
      xa = xamin + xajac*yy(2)
      xb = tau/x/xa   
c... This integratio is to be performed using x kinematics.
c... till here linear mapping is added

c     AllReg(1) = 0d0
c     AllReg(2) = 0d0
c      PKReg = 0.0d0
      sp = xa*xb*S  ! since this sp goes into the flux, 
!     sp = 2d0*dot(p1,p2)
!     rsp = dsqrt(sp)
c... i need to think weather sp = 2d0*dot(p1,p2) as while generating
c... kinematics itself we have used the x kinematics, ans xaxbS id not 
c... equal to 2dot(p1,p2). 
      rsp = dsqrt(sp)

      do k=1,2

        if ( k .eq. 1) call kinvar1(xa*x,xb,p1,p2,p3)
        if ( k .eq. 2) call kinvar1(xa,xb*x,p1,p2,p3)
       sp = 2d0*dot(p1,p2)
       rsp = dsqrt(sp)
c... same reasoning as plusA. refer there.

            AL = alphasPDF(xmur)
            ALP = AL/2d0/Pi

!            call getPKReg(x,xmuf,p1,p2,p3,SumReg)
            coef = Born_gg2h_(0,AL,p1,p2,p3) 
c... coefficients always comes with the x kinematics

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call setlum(f1,f2,xl)
            
      Tb_dot_Ta_by_Ta2 = -1d0
      xmuf2 = xmuf**2
      s12 = 2d0*dot(p1,p2)
      Tfactors = Tb_dot_Ta_by_Ta2

      Pterm =ALP*PK(P_,gg_,Plus_,x)*Tfactors*dlog(xmuf2/s12)

      Kterm =    ALP*PK(Kb_,gg_,Plus_,x) 
     & - ALP*Tfactors*PK(Kt_,gg_,Plus_,x)

      sig =xl(2)*coef*(Pterm + Kterm)

           pi_1 = PI/amH 
           flux = 2d0*sp
          xnorm = hbarc2*pi_1/flux

      PK_leg(k) = xnorm*xajac*xjac4*sig* 2d0*amH/xa/x/S

         enddo
         xint_PlusA  = PK_leg(1) + PK_leg(2)
c... this is the first part which will be countered by the finite 
c... pieces comming from the next routine.
c...      Integrate[ A,{x,0,1-δ}] + 𝝳(1-x) g(1) [h(0) -h(1-δ)]
151      return
      end

c... deltaa
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Delta terms]
      function xint_Plus_h(yy,vwgt)
      implicit double precision (a-h,o-z)
      character*10 Pgg_,Ktgg_,Kbargg_
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
      double precision xint_h
      external xint_h

      Pgg_ = 'Pgg'
      Ktgg_= 'Ktgg'
      Kbargg_='Kbargg'

        delta = 1d-2
        tau = amH**2/S
         xajac = 1d0 - tau
         xa = tau + xajac*yy(1)
         xb = tau/xa

         sp = xa*xb*S
        rsp = dsqrt(sp)

        call kinvar1(xa,xb,p1,p2,p3)
!        xmuf= amH/2d0   !.....already defined in main through common 
!        xmur= xmuf
        xmuf2=xmuf**2

         AL = alphasPDF(xmur)
         ALP= AL/2d0/PI

        call pdf(xa,xmuf,f1)
        call pdf(xb,xmuf,f2)
        call setlum(f1,f2,xl)

c        call getPKDel(1.0d0,xmuf,p1,p2,p3,SumDel)

        s12 = 2d0*dot(p1,p2)
        Pplus = (-1.0d0)*dlog(xmuf2/s12)
c... whatever overall factors are there needs to be carefully used 
c... here.
     
      coef = Born_gg2h_(0,AL,p1,p2,p3) 
      zero = 0d0
      one  = 1d0
      one_m_del = one-delta

      Pplus = Pplus*(xint_h(zero,Pgg_) - xint_h(one_m_del,Pgg_))
      xKbar_h = xint_h(zero,Kbargg_) - xint_h(one_m_del,Kbargg_)
      xKt_h = xint_h(zero,Ktgg_) - xint_h(one_m_del,Ktgg_)

c	print*,Pplus
c	print*,xint_h(0d0,Pgg_)
c	print*,xint_h(1d0-delta,Pgg_)
c	stop

      SumPlus = Pplus + xKbar_h + xKt_h

         sig = 2d0*xl(2)*SumPlus*coef

         pi_1 = PI/amH  
         flux = 2d0*sp

        xnorm=ALP*hbarc2*pi_1/flux

        xint_Plus_h = xajac * xnorm * sig * 2d0*amh/xa/S

      return
      end
c------------------------------------------------------------------
