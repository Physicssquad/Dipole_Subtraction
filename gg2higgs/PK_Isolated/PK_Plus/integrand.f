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
      common/cuts_delta/delta
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

c      delta = 1d-5
      xmin =0d0 ! technically xmin should run from 0 to 1. 
      xmin = tau
      xmax = 1.0d0 - delta
      xjac4 = (xmax - xmin)
      x = xjac4*yy(1) + xmin
c... for the fixed value of x. I will generate the xa and xb using the 
c... constraint relation.
c... x will always be there as an overall factor.
c~~~~~[ _______________ ]~~~~~C
      xamin = 0d0
      xamin = tau/x
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

c         sp = xa*xb*S
c        rsp = dsqrt(sp)
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
      common/cuts_delta/delta
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
c      delta = 1d-5 
      xmin = 0d0 
      xmin = tau
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
      common/cuts_delta/delta
      double precision xint_h
      external xint_h

      Pgg_ = 'Pgg'
      Ktgg_= 'Ktgg'
      Kbargg_='Kbargg'

!        delta = 1d-5
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

      SumPlus = Pplus + xKbar_h + xKt_h

         sig = 2d0*xl(2)*SumPlus*coef

         pi_1 = PI/amH  
         flux = 2d0*sp

        xnorm=ALP*hbarc2*pi_1/flux

        xint_Plus_h = xajac * xnorm * sig * 2d0*amh/xa/S

      return
      end
c------------------------------------------------------------------
c...KEY WORDS TO JUMP
c... plusa
c... plusb
c... regular
c... deltaa
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[PlusA term]
c... in run_test_03 I gave run with regular x to run from 0 to 1.
      function flo2_PlusA(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10),PKt(1:2)
      dimension f1(-6:6),f2(-6:6),f2a(-6:6),xl(15),xla(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension xp1(0:3),xp2(0:3)
      double precision Pterm, Kterm
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+9)   ! I will be using fb as matrix result.
      character*50 name
      common/amass/am1,am2,amH,am4,am5
      common/energy/s
      common/pdfname/name
c      common/usedalpha/AL,ge
      common/scales/xmuf,xmur
      common/cuts_delta/delta

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

        PKt(1) = 0.0d0
        PKt(2) = 0.0d0
!      delta = 1d-5 
c... delta defined in main 

c~~~~~[ Exponential MAPPINNG ]~~~~~C
!      almin = delta
!      almax = 1.0d0
!      al = almin*(almax/almin)**yy(1)
!      xjac4 = al*dlog(almax/almin)
!      x = 1.0d0 - al
c~~~~~[ _______________ ]~~~~~C

c~~~~~[ LINEAR MAPPINNG ]~~~~~C
c      xmin = tau ! to avoid xa hitting inf x min is to be tau. 
c... even with xmin = 0d0. Result is same as tau
      xmin =0d0 ! technically xmin should run from 0 to 1. 
!      xmin = tau
      xmax = 1.0d0 - delta
      xjac4 = (xmax - xmin)
      x = xjac4*yy(1) + xmin
c... for the fixed value of x. I will generate the xa and xa using the 
c... constraint relation.
c... x will always be there as an overall factor.
c~~~~~[ _______________ ]~~~~~C
!      xamin = tau/x
      xamin = 0d0 
      xamax = 1d0
      xajac =  (xamax - xamin)
      xa = xamin + xajac*yy(2)
      xb = tau/x/xa    
!      if(xb .ge. 1) print*,xb,yy(1),yy(2)
c... here x* xa*xb*s = amH**2 is the kinametics.
c... xb is generated for every case with Higgs mass constraint.

!      sp = xa*xb*S
!      rsp = dsqrt(sp)
c... this sp has a problem as it is not always equal to s12. as the kinematics we are using 
c... is the reduced one.[?]

        do k = 1,2

        if ( k .eq. 1) call kinvar1(xa*x,xb,p1,p2,p3)
        if ( k .eq. 2) call kinvar1(xa,xb*x,p1,p2,p3)
        
c... so due to the doubt in sp, I will use s12 as the partonic, sp
           sp  = 2d0*dot(p1,p2)
           rsp = dsqrt(sp)

            xmuf2 = xmuf*xmuf 
               AL = alphasPDF(xmur) 
              ALP = AL/2d0/Pi

c This born used x kinematics as 1-x fraction of momenta is taken by gluon radiation
!        coef = Born_gg2h(0,p1,p2,p3)
        coef = Born_gg2h_(0,AL,p1,p2,p3)

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)

            call setlum(f1,f2,xl)

!            call getPKPlus(1,x,xmuf,p1,p2,p3,SumPlus)
C... Modified here!
      s12 = 2d0*dot(p1,p2)
      Tb_dot_Ta_by_Ta2 = -1d0
      Tfactors = Tb_dot_Ta_by_Ta2
      Pterm =ALP*PK(P_,gg_,Plus_,x)*Tfactors*dlog(xmuf2/s12)
      Kterm =ALP*PK(Kb_,gg_,Plus_,x) -
     &  ALP*Tfactors*PK(Kt_,gg_,Plus_,x)

      SumPlus =coef*(Pterm + Kterm)

            sig = xl(2)*SumPlus
  
           pi_1 = PI/amH
           flux = 2d0*sp
          xnorm = hbarc2*pi_1/flux

          PKplus_x = xnorm*xajac*xjac4*sig*2d0*amH/xa/x/S

          PKt(k) = PKplus_x
            enddo

        flo2_PlusA = (PKt(1) + PKt(2))

c... ALP overall factor taken here
c... is symmetry factor required here as there are symmetry in 
c... the initial state. 
151      return
      end

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[PlusB term]
c... plusb
      function flo2_PlusB(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10),PKt(1:2)
      dimension f1(-6:6),f2(-6:6),f2a(-6:6),xl(15),xla(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension xp1(0:3),xp2(0:3)
      double precision Pterm, Kterm
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+9)
      character*50 name
      common/energy/s
      common/pdfname/name
      common/distribution/xq
      common/bin_size/eps
      common/amass/am1,am2,amH,am4,am5
      common/scales/xmuf,xmur
      common/cuts_delta/delta

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

!      delta = 1.0d-5
c~~~~~[ LINEAR MAPPINNG ]~~~~~C
c... minimum should not be 0d0 as xamin will be infinity, 
c... but numerically i am taking as random numbers wonn't hit 0.
!      xmin = tau  
      xmin = 0d0 
      xmax = 1.0d0 - delta
      xjac4 = (xmax - xmin)
      x = xjac4*yy(1) + xmin

c~~~~~[ Exponential MAPPINNG ]~~~~~C
c      almin = delta 
c      almax = 1.0d0
c      al = almin*(almax/almin)**yy(2)
c      xjac4 = al*dlog(almax/almin)
c      x = 1.0d0 - al
c      xamin = tau/x 
      xamin = 0d0 
      xamax = 1d0
      xajac =  (xamax - xamin)
      xa = xamin + xajac*yy(2)
      xb = tau/xa    


        flo2_PlusB  = 0.0d0
        PKplus_x = 0.0d0

        sig = 0.0d0
        PKt(1) = 0.0d0
        PKt(2) = 0.0d0

        do k = 1,2

        call kinvar1(xa,xb,p1,p2,p3) 
c... Only difference here from fnlo_A is that coef i.e 
c... born has pure born kinematics earlier it was reduced.
c... But the kinematics is to be defined in terms of x 
c... kinematics as in 59 above.
          sp = xa*xb*S 
         rsp = dsqrt(sp)
c... refer to the reasoning in plusA for sp
!           sp  = 2d0*dot(p1,p2)
!           rsp = dsqrt(sp)

c... muf is also defined in the main.
       xmuf2 = xmuf*xmuf 
          AL = alphasPDF(xmur) 
         ALP = AL/2d0/Pi
        coef = Born_gg2h_(0,AL,p1,p2,p3)

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call setlum(f1,f2,xl)

c... Same plus function is being called here but the coef are already 
c... with x=1 kinematics, so the cuts where x.ne.1 is used as (x-delta)
c            call getPKPlus(0,x,xmuf,p1,p2,p3,SumPlus)
      s12 = 2d0*dot(p1,p2)
      Tb_dot_Ta_by_Ta2 = -1d0
      Tfactors = Tb_dot_Ta_by_Ta2
      Pterm =ALP*PK(P_,gg_,Plus_,x)*Tfactors*dlog(xmuf2/s12)
      Kterm =ALP*PK(Kb_,gg_,Plus_,x) -
     &  ALP*Tfactors*PK(Kt_,gg_,Plus_,x)

      SumPlus =coef*(Pterm + Kterm)
           
      sig = xl(2)*SumPlus

           pi_1 = PI/amH
           flux = 2d0*sp
          xnorm = hbarc2*pi_1/flux

         PKplus_x = xnorm*xajac*xjac4*sig* 2d0*amH/xa/S
c... here also in the denominator x is to be multiplied as the
c... kinematics is generated in that way. [?]
c        PKplus_x = xnorm*xajac*xjac4*sig* 2d0*amH/xa/S

         PKt(k) = PKplus_x

            enddo

        flo2_PlusB = (PKt(1) + PKt(2))
151      return
      end
