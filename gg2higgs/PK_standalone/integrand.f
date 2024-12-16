c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[PlusA term]
      function flo2_PlusA(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10),PK(1:2)
      dimension f1(-6:6),f2(-6:6),f2a(-6:6),xl(15),xla(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension xp1(0:3),xp2(0:3)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+9)
      character*50 name
      common/amass/am1,am2,amH,am4,am5
      common/energy/s
      common/pdfname/name
      common/usedalpha/AL,ge

      tau = amH**2/S
      xajac = 1d0 - tau
      xa = tau + xajac*yy(1)
      xb = tau/xa
c      sp = xa*xb*S

c~~~~~[ LINEAR MAPPINNG ]~~~~~C
c      xmin = 0.0d0
c      xmax = 1.0d0 -delta
c      xjac4 = (xmax - xmin)
c      x = xjac4*yy(4) + xmin

      delta = 1.0d-5
      almin = delta
      almax = 1.0d0
      al = almin*(almax/almin)**yy(2)
      xjac4 = al*dlog(almax/almin)
      x = 1.0d0 - al

        flo2_PlusA  = 0.0d0
        PKplus_x = 0.0d0

        sig = 0.0d0
        PK(1) = 0.0d0
        PK(2) = 0.0d0
       
        eps = 0.5
        amH_min = amH - eps
        amH_max = amH + eps
        do k = 1,2

        if ( k .eq. 1) call kinvar1(xa*x,xb,p1,p2,p3)
        if ( k .eq. 2) call kinvar1(xa,xb*x,p1,p2,p3)
c	  print*,"sp1:",dsqrt(dot(p3,p3))
          sp = 2d0*dot(p1,p2)
         rsp = dsqrt(sp)
	if (rsp .ge. amH_min ) then 

        coef = Born_gg2h(0,p1,p2,p3)

             xmuf = amH 
             xmur = xmuf
            xmuf2 = xmuf*xmuf 
               AL = alphasPDF(xmur) 
              ALP = AL/2d0/Pi

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)

            call setlum(f1,f2,xl)

            call getPKPlus(1,x,xmuf,p1,p2,p3,SumPlus)
c	print*,SumPlus
c	stop
           
            sig = xl(2)*SumPlus*coef
  
           pi_1 = PI/amH
           flux = 2d0*sp
          xnorm = hbarc2*pi_1/flux

         PKplus_x = xnorm*xajac*xjac4*sig*2d0*amH/xa/S

          PK(k) = PKplus_x
	endif
            enddo

        flo2_PlusA = ALP * (PK(1) + PK(2))    ! ALP overall factor taken here
      return
      end

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[PlusB term]
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
      common/factscale/xmuf
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

      delta = 1.0d-5
      almin = delta
      almax = 1.0d0
      al = almin*(almax/almin)**yy(2)
      xjac4 = al*dlog(almax/almin)
      x = 1.0d0 - al

        flo2_PlusB  = 0.0d0
        PKplus_x = 0.0d0

        sig = 0.0d0
        PK(1) = 0.0d0
        PK(2) = 0.0d0

        do k = 1,2

        call kinvar1(xa,xb,p1,p2,p3)
          sp = xa*xb*S 
          sp = 2d0*dot(p1,p2)
         rsp = dsqrt(sp)


        coef = Born_gg2h(0,p1,p2,p3)

             xmuf = amH
             xmur = xmuf
            xmuf2 = xmuf*xmuf 
               AL = alphasPDF(xmur) 
              ALP = AL/2d0/Pi

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)

            call setlum(f1,f2,xl)

            call getPKPlus(0,x,xmuf,p1,p2,p3,SumPlus)
           
            sig = xl(2)*SumPlus*coef

           pi_1 = PI/amH
           flux = 2d0*sp
          xnorm = hbarc2*pi_1/flux

         PKplus_x = xnorm*xajac*xjac4*sig* 2d0*amH/xa/S

          PK(k) = PKplus_x

c            endif
            enddo

        flo2_PlusB = ALP * (PK(1) + PK(2))
      return
      end

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Regular Terms]
      function flo2_PKReg(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .          ,p(0:3,1:4),Born(1:2),AllReg(1:2)
      dimension SumP(1:2),SumK(1:2)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+9)
      common/energy/s
      common/factscale/xmuf
      common/usedalpha/AL,ge
      common/distribution/xq
      common/bin_size/eps
      common/amass/am1,am2,amH,am4,am5

      tau = amH**2/S
      xajac = 1d0 - tau
      xa = tau + xajac*yy(1)
      xb = tau/xa
      x  = yy(2)

      AllReg(1) = 0d0
      AllReg(2) = 0d0
      PKReg = 0.0d0
      sp = xa*xb*S

      do k=1,2

        if ( k .eq. 1) call kinvar1(xa*x,xb,p1,p2,p3)
        if ( k .eq. 2) call kinvar1(xa,xb*x,p1,p2,p3)

c        sp = 2.0d0*dot(p1,p2)
        rsp = dsqrt(sp)


c        if ( sp .ge. amH**2 ) then
c	if (sp .ge. amH**2) then !print*,sp,amH**2

            xmuf = amH
            xmur = xmuf 

            AL = alphasPDF(xmur)
            ALP = AL/2d0/Pi

            call getPKReg(x,xmuf,p1,p2,p3,SumReg)
            coef = Born_gg2h(0,p1,p2,p3)

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call setlum(f1,f2,xl)

            sig1 = xl(2)* SumReg  

            sig = Alp*sig1*coef

           pi_1 = PI/amH
           flux = 2d0*sp
          xnorm = hbarc2*pi_1/flux

         PKplus_x = xnorm*xajac*sig* 2d0*amH/xa/S


           AllReg(k) = PKplus_x  

c         endif
         enddo

         flo2_PKReg = AllReg(1) + AllReg(2)

      return
      end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Delta terms]
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
      common/factscale/xmuf
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


        xmuf= amH
        xmur= xmuf
        xmu2=xmuf**2

         AL = alphasPDF(xmur)
         ALP= AL/2d0/PI

        call pdf(xa,xmuf,f1)
        call pdf(xb,xmuf,f2)
        call setlum(f1,f2,xl)


        call getPKDel(1.0d0,xmuf,p1,p2,p3,SumDel)

         sig = xl(2)* SumDel  !  [qq lum]

         pi_1 = PI/amH
         flux = 2d0*sp

        xnorm=ALP*hbarc2*pi_1/flux

        flo2_PKDel  = xajac * xnorm * sig * 2d0* amh/xa/S
c 	print*,flo2_PKDel,SumDel,xnorm,sig

      return
      end
c---------------------------------------------------------------------
