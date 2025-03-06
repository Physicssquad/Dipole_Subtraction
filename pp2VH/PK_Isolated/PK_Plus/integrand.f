      function flo2_Plus(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10),PK(1:2)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension xp1(0:3),xp2(0:3)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+6)
      character*50 name
      common/energy/s
      common/pdfname/name
      common/factscale/xmuf
      common/usedalpha/AL,ge
      common/distribution/xq

      xa     = yy(1)
      xb     = yy(2)
      xc     = yy(3)

c      xmin   = 0.0d0
c      xmax   = 1.0d0 - 1d-5
c      xjac4   = (xmax-xmin)
c      x      = xmin+ xjac4*yy(4)

      delta  = 1d-5
      almin = delta
      almax = 1.0d0
      al = almin*(almax/almin)**yy(4)
      aljac = al*dlog(almax/almin)
      xjac4 = aljac
      x = 1.0d0-al


      xjact = 2.0d0
      xjac = xjact*xjac4

      xnorm=hbarc2

      eps = 0.5d0
      xlow = xq - eps
      xhigh = xq + eps


        flo2_Plus  = 0.0d0
        PKplus_x = 0.0d0
        PKplus_1 = 0.0d0
        PKReg    = 0.0d0
        PKDel    = 0.0d0
        PK(1) = 0d0
        PK(2) = 0d0

        sig1 = 0.0d0
        sig2 = 0.0d0
        sig3 = 0.0d0
        sig4 = 0.0d0

        do k = 1,2

        do iplus = 0,1

        if (iplus .eq. 1) then

         if (k .eq. 1) call kinvar2_PK(x*xa,xb,xc,Qmass,p1,p2,p3,p4)
         if (k .eq. 2) call kinvar2_PK(xa,x*xb,xc,Qmass,p1,p2,p3,p4)

        elseif (iplus .eq. 0) then
        call kinvar2_PK(xa,xb,xc,Qmass,p1,p2,p3,p4)
        endif
        call cuts0(p1,p2,p3,p4,ipass)

         scale = Qmass

        if (scale .ge. xlow .and. scale .le. xhigh  
     .   .and.   ipass .eq. 1)  then
c     .    ) then

         coef = Born_gg2aa(0,p1,p2,p3,p4)
c	coef = 1
c         print*,"coef",coef
c         print*,"p",p1
c         print*,"p",p2
c         print*,"p",p3
c         print*,"p",p4

         sp   = 2.0d0*dot(p1,p2)
         rsp  = dsqrt(sp) 
         pin  = 0.5d0*rsp
         pf   = 0.5d0*rsp
         flux = 4.0d0*pin*rsp
         
             xmuf = scale
             xmur = xmuf
            xmuf2 = xmuf*xmuf 
               AL = alphasPDF(xmur) 
              ALP = AL/2d0/Pi

            azmth = 2.0d0*pi
              ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call setlum(f1,f2,xl)

            call getPKPlus(iplus,x,xmuf,p1,p2,p3,p4,SumPlus)
           
            if (iplus .eq. 1) sig1 = xl(4)*SumPlus*coef
            if (iplus .eq. 0) sig3 = xl(4)*SumPlus*coef 
c~~[ Here sig is multiplied by 2d0 ???? so removed it ]~~~

c            if (iplus .eq. 1) sig1 = SumPlus
c            if (iplus .eq. 0) sig3 = SumPlus
  
          if (iplus .eq. 1) then
            wgt1 = sig1/flux*ps2*xjac*vwgt
            PKplus_x = xnorm*wgt1/vwgt/2d0/eps

          elseif (iplus .eq.0) then
            wgt3 = sig3/flux*ps2*xjac*vwgt
            PKplus_1= xnorm*wgt3/vwgt/2d0/eps

           endif

           endif                ! bin choice
           enddo

            PK(k) = PKplus_x - PKplus_1

         enddo

        flo2_Plus = ALP * ( PK(1) + PK(2) )
      return
      end
c---------------------------------------------------------------------
