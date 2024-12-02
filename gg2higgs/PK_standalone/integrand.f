c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[PlusA term]
      function flo2_PlusA(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10),PK(1:2)
      dimension f1(-6:6),f2(-6:6),f2a(-6:6),xl(15),xla(15)
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
      common/bin_size/eps

      xlow = xq - eps
      xhigh = xq + eps

        x1 = yy(1)
        x2 = yy(2)
        xt = yy(3)

c      xmin = 0.0d0
c      xmax = 1.0d0 -delta
c      xjac4 = (xmax - xmin)
c      x = xjac4*yy(4) + xmin

      delta = 1.0d-5
      almin = delta
      almax = 1.0d0
      al = almin*(almax/almin)**yy(4)
      xjac4 = al*dlog(almax/almin)
      x = 1.0d0 - al

      xtjac = 2.0d0
      xjac = xtjac*xjac4

      xnorm=hbarc2

        flo2_PlusA  = 0.0d0
        PKplus_x = 0.0d0

        sig = 0.0d0
        PK(1) = 0.0d0
        PK(2) = 0.0d0
        Qmass = 0.0d0


        do k = 1,2

        if ( k .eq. 1) then
        call kinvar2_PK(x*x1,x2,xt,Qmass,p1,p2,p3,p4)
        scale = Qmass
        elseif ( k .eq. 2) then
        call kinvar2_PK(x1,x*x2,xt,Qmass,p1,p2,p3,p4)
        scale = Qmass
        endif
        call cuts0(p1,p2,p3,p4,ipass)

        if ( scale .ge. xlow .and. scale .le. xhigh
     .     .and.   ipass .eq. 1)  then
c     .       ) then
        coef = Born_gg2aa(0,p1,p2,p3,p4)

        sp   =  2.0d0*dot(p1,p2)
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

            call pdf(x1,xmuf,f1)
            call pdf(x2,xmuf,f2)

            call setlum(f1,f2,xl)

            call getPKPlus(1,x,xmuf,p1,p2,p3,p4,SumPlus)
           
            sig = xl(4)*SumPlus*coef
c            sig =SumPlus
  
            wgt = sig/flux*ps2*xjac*vwgt
            PKplus_x = xnorm*wgt/vwgt

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
      parameter (hbarc2=389.3856741D+6)
      character*50 name
      common/energy/s
      common/pdfname/name
      common/factscale/xmuf
      common/usedalpha/AL,ge
      common/distribution/xq
      common/bin_size/eps

      xlow = xq - eps
      xhigh = xq + eps

        x1 = yy(1)
        x2 = yy(2)
        xt = yy(3)

c      xmin = 0.0d0
c      xmax = 1.0d0 -delta
c      xjac4 = (xmax - xmin)
c      x = xjac4*yy(4) + xmin

      delta = 1.0d-5
      almin = delta
      almax = 1.0d0
      al = almin*(almax/almin)**yy(4)
      xjac4 = al*dlog(almax/almin)
      x = 1.0d0 - al

      xtjac = 2.0d0
      xjac = xtjac*xjac4

      xnorm=hbarc2

        flo2_PlusB  = 0.0d0
        PKplus_x = 0.0d0

        sig = 0.0d0
        PK(1) = 0.0d0
        PK(2) = 0.0d0

        do k = 1,2

        call kinvar2_PK(x1,x2,xt,Qmass,p1,p2,p3,p4)
        call cuts0(p1,p2,p3,p4,ipass)

        scale = Qmass

        if ( scale .ge. xlow .and. scale .le. xhigh
     .     .and.   ipass .eq. 1)  then
c     .       ) then

        coef = Born_gg2aa(0,p1,p2,p3,p4)

        sp   =  2.0d0*dot(p1,p2)
        rsp  = dsqrt(sp)
        pin  = 0.5d0*rsp
        pf   = 0.5d0*rsp
        flux = 4.0d0*pin*rsp
         
             xmuf = scale
             xmur = xmuf
            xmuf2 = xmuf*xmuf 
               AL = alphasPDF(xmur) 
              ALP = AL/2d0/Pi
c              ALP = 1.0d0

            azmth = 2.0d0*pi
            ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth

            call pdf(x1,xmuf,f1)
            call pdf(x2,xmuf,f2)

            call setlum(f1,f2,xl)

            call getPKPlus(0,x,xmuf,p1,p2,p3,p4,SumPlus)
           
c             sig = SumPlus
            sig = xl(4)*SumPlus*coef
  
            wgt = sig/flux*ps2*xjac*vwgt
            PKplus_x = xnorm*wgt/vwgt

            PK(k) = PKplus_x
            endif
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
      parameter (hbarc2=389.3856741D+6)
      common/energy/s
      common/factscale/xmuf
      common/usedalpha/AL,ge
      common/distribution/xq
      common/bin_size/eps

      xa     = yy(1)
      xb     = yy(2)
      xc     = yy(3)
      x      = yy(4)

      xjac   = 2d0
      xnorm=hbarc2

      xlow = xq - eps
      xhigh = xq + eps

      AllReg(1) = 0d0
      AllReg(2) = 0d0
      PKReg = 0.0d0

      do k=1,2

      if (k .eq. 1) call kinvar2_PK(x*xa,xb,xc,Qmass,p1,p2,p3,p4)
      if (k .eq. 2) call kinvar2_PK(xa,x*xb,xc,Qmass,p1,p2,p3,p4)

        sp = 2.0d0*dot(p1,p2)
        rsp = dsqrt(sp)
        scalex = Qmass
        pin  = 0.5d0*rsp
        flux = 4.0d0*pin*rsp


        if (scalex .ge. xlow .and. scalex .le. xhigh) then

            xmuf = scalex
            xmur = scalex

            AL = alphasPDF(xmur)
            ALP = AL/2d0/Pi

            call getPKReg(x,xmuf,p1,p2,p3,p4,SumReg)
            coef = Born_gg2aa(0,p1,p2,p3,p4)

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call setlum(f1,f2,xl)

            sig1 = xl(4)* SumReg  !  [qq lum]

            sig = Alp*sig1*coef

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
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Delta terms]
      function flo2_PKDel(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
     .         ,f1(-6:6),f2(-6:6),xl(15)
     .         ,p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .         ,p(0:3,1:4),Born(1:2)
     .         ,SumP(1:2),SumK(1:2)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+6)
      common/energy/s
      common/factscale/xmuf
      common/usedalpha/AL,ge
      common/distribution/xq
      common/bin_size/eps

c      xa     = yy(1)
c      xb     = yy(2)

      TAUH = xq*xq/s

      xa = (1D0-TAUH)*yy(1)    + TAUH
c      xb = (1D0-TAUH/xa)*yy(2) + TAUH/xa
      xb = TAUH/xa

      xc = yy(2)

      sp     = xa*xb*s
      rsp    = dsqrt(sp)
      xjac = 2.0d0*(1.0d0-TAUH)
      xnorm=hbarc2*2.0d0*xq/s/xa

        call kinvar2_PK(xa,xb,xc,xinvmass,p1,p2,p3,p4)

        call p1d_to_p2d_4(p1,p2,p3,p4,p)

        pin  = 0.5d0*rsp
        flux_1 = 4.0d0*pin*rsp

        scale = xq

      flo2_PKDel  = 0.0d0
      sig1 = 0d0


            xmuf = scale
            xmur = scale
            AL = alphasPDF(xmur)
            ALP = AL/2d0/Pi

            call getPKDel(1.0d0,xmuf,p,p1,p2,SumDel)

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call setlum(f1,f2,xl)

            sig1 = xl(4)* SumDel  !  [qq lum]

            sig1 = Alp*sig1

            azmth = 2.0d0*pi
            pf = 0.5d0*rsp
            ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth

            wgt1 = sig1/flux_1*ps2*xjac*vwgt
            PKDel = xnorm*wgt1/vwgt
c            write(*,*)'check =', sig1,flux1,ps2,xjac 

         flo2_PKDel = PKDel

      return
      end
c---------------------------------------------------------------------
