c...KEY WORDS TO JUMP
c... plusa
c... plusb
c... regular
c... deltaa
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[PlusA term]
c... in run_test_03 I gave run with regular x to run from 0 to 1.
      function flo2_PlusA(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10),PK(1:2)
      dimension f1(-6:6),f2(-6:6),f2a(-6:6),xl(15),xla(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension xp1(0:3),xp2(0:3)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+9)   ! I will be using fb as matrix result.
      character*50 name
      common/amass/am1,am2,amH,am4,am5
      common/energy/s
      common/pdfname/name
c      common/usedalpha/AL,ge
      common/scales/xmuf,xmur
      common/cuts_delta/delta

      tau = amH**2/S
c      delta = 1d-5 
c... delta defined in main 

c~~~~~[ Exponential MAPPINNG ]~~~~~C
c      almin = delta
c      almax = 1.0d0
c      al = almin*(almax/almin)**yy(1)
c      xjac4 = al*dlog(almax/almin)
c      x = 1.0d0 - al
c~~~~~[ _______________ ]~~~~~C

c~~~~~[ LINEAR MAPPINNG ]~~~~~C
c      xmin = tau ! to avoid xa hitting inf x min is to be tau. 
c... even with xmin = 0d0. Result is same as tau
      xmin =0d0 ! technically xmin should run from 0 to 1. 
      xmax = 1.0d0 - delta
      xjac4 = (xmax - xmin)
      x = xjac4*yy(1) + xmin
c... for the fixed value of x. I will generate the xa and xa using the 
c... constraint relation.
c... x will always be there as an overall factor.
c~~~~~[ _______________ ]~~~~~C
c      xamin = tau/x
      xamin = 0d0 
      xamax = 1d0
      xajac =  (xamax - xamin)
      xa = xamin + xajac*yy(2)
      xb = tau/x/xa    
c... here x* xa*xb*s = amH**2 is the kinametics.
c... xb is generated for every case with Higgs mass constraint.

      sp = xa*xb*S
      rsp = dsqrt(sp)
c... this sp has a problem as it is not always equal to s12. as the kinematics we are using 
c... is the reduced one.[?]

        do k = 1,2

        if ( k .eq. 1) call kinvar1(xa*x,xb,p1,p2,p3)
        if ( k .eq. 2) call kinvar1(xa,xb*x,p1,p2,p3)
        
c... so due to the doubt in sp, I will use s12 as the partonic, sp
c           sp  = 2d0*dot(p1,p2)
c           rsp = dsqrt(sp)

            xmuf2 = xmuf*xmuf 
               AL = alphasPDF(xmur) 
              ALP = AL/2d0/Pi

c This born used x kinematics as 1-x fraction of momenta is taken by gluon radiation
c        coef = Born_gg2h(0,p1,p2,p3)
        coef = Born_gg2h_(0,AL,p1,p2,p3)

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)

            call setlum(f1,f2,xl)

            call getPKPlus(1,x,xmuf,p1,p2,p3,SumPlus)
           
            sig = xl(2)*SumPlus*coef
  
           pi_1 = PI/amH
           flux = 2d0*sp
          xnorm = hbarc2*pi_1/flux

          PKplus_x = xnorm*xajac*xjac4*sig*2d0*amH/xa/x/S

c.....[?]  x in the denomonator ?
c          PKplus_x = xnorm*xajac*xjac4*sig*2d0*amH/xa/S

          PK(k) = PKplus_x
            enddo

        flo2_PlusA = ALP * (PK(1) + PK(2))

c... ALP overall factor taken here
c... is symmetry factor required here as there are symmetry in 
c... the initial state. 
      return
      end

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[PlusB term]
c... plusb
      function flo2_PlusB(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10),PK(1:2)
      dimension f1(-6:6),f2(-6:6),f2a(-6:6),xl(15),xla(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension xp1(0:3),xp2(0:3)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+9)
      character*50 name
      common/energy/s
      common/pdfname/name
c      common/factscale/xmuf
c      common/usedalpha/AL,ge
      common/distribution/xq
      common/bin_size/eps
      common/amass/am1,am2,amH,am4,am5
      common/scales/xmuf,xmur
      common/cuts_delta/delta

      tau = amH**2/S

c      delta = 1.0d-5
c~~~~~[ LINEAR MAPPINNG ]~~~~~C
c      xmin = tau  
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
c      xa = tau + xajac*yy(2)
      xa = xamin + xajac*yy(2)
      xb = tau/xa    

      sp = xa*xb*S
      rsp = dsqrt(sp)



        flo2_PlusB  = 0.0d0
        PKplus_x = 0.0d0

        sig = 0.0d0
        PK(1) = 0.0d0
        PK(2) = 0.0d0

        do k = 1,2

        call kinvar1(xa,xb,p1,p2,p3) 
c... Only difference here from fnlo_A is that coef i.e 
c... born has pure born kinematics earlier it was reduced.
c... But the kinematics is to be defined in terms of x 
c... kinematics as in 59 above.
c          sp = xa*xb*S 
c         rsp = dsqrt(sp)
c... refer to the reasoning in plusA for sp
           sp  = 2d0*dot(p1,p2)
           rsp = dsqrt(sp)

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
            call getPKPlus(0,x,xmuf,p1,p2,p3,SumPlus)
           
            sig = xl(2)*SumPlus*coef

           pi_1 = PI/amH
           flux = 2d0*sp
          xnorm = hbarc2*pi_1/flux

         PKplus_x = xnorm*xajac*xjac4*sig* 2d0*amH/xa/S
c... here also in the denominator x is to be multiplied as the
c... kinematics is generated in that way. [?]
c        PKplus_x = xnorm*xajac*xjac4*sig* 2d0*amH/xa/S

         PK(k) = PKplus_x

            enddo

        flo2_PlusB = ALP * (PK(1) + PK(2))
      return
      end

c... regular
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
c      xajac = 1d0 - tau
c      xa = tau + xajac*yy(1)
c      xb = tau/xa
c      x  = yy(2)

c... for regular terms I don't need to use the upper cut as it is a flat integral.
      delta = 0d0  
!      xmin = tau 
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

            call getPKReg(x,xmuf,p1,p2,p3,SumReg)
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
      return
      end

c... deltaa
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
c        xmuf= amH/2d0   !.....already defined in main through common 
c        xmur= xmuf
        xmu2=xmuf**2

         AL = alphasPDF(xmur)
         ALP= AL/2d0/PI

        call pdf(xa,xmuf,f1)
        call pdf(xb,xmuf,f2)
        call setlum(f1,f2,xl)

        call getPKDel(1.0d0,xmuf,p1,p2,p3,SumDel)

         sig = xl(2)* SumDel 

         pi_1 = PI/amH  
         flux = 2d0*sp

        xnorm=ALP*hbarc2*pi_1/flux

        flo2_PKDel  = xajac * xnorm * sig * 2d0*amh/xa/S

      return
      end
c------------------------------------------------------------------
