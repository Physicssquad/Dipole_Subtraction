      function fnlo2(zz,vwgt)
      implicit double precision (a-h,o-z)
      dimension zz(10)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)

      dimension f1(-6:6),f2(-6:6)
      dimension f1t(-6:6),f2t(-6:6)
      dimension f1tqq(-6:6),f2tqq(-6:6)
      dimension f1tqg(-6:6),f2tqg(-6:6)
      dimension f1tgq(-6:6),f2tgq(-6:6)
      dimension f1tgg(-6:6),f2tgg(-6:6)

      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=0.3894d9)

      common/energy/s
      common/momenta/p1,p2,p3,p4
      common/bin/xq,xeps,xlow,xhigh
      common/cfbin/cf,cfeps,clow,chigh
      common/yfbin/yf,yfeps,ylow,yhigh
      common/xrange/xrlow,xrhigh
      common/slice/deltas,deltac
      common/factscale/xmuf
      common/param/aem,xmur,lambda
      common/xmcoeff/xc1,xc2
      common/isub/io,is
      common/angle/cststar

      call kinvar2(zz,
     &     xinvmass,y1,y2,Y,cst1,cst2,pt1,pt2)

      call cuts2(xinvmass,y1,y2,Y,cst1,cst2,pt1,pt2,ipass)

      xa =zz(1)
      xb =zz(2)
      scale =xinvmass
      xjac=1d0
      vwgt=1d0

      pobl = cststar

      IF(ipass.eq.1)THEN

         if(   scale.gt.xlow .and. scale.lt.xhigh
c         if(   scale.gt.xrlow .and. scale.lt.xrhigh
c     &  .and. pobl.ge.clow .and. pobl.le.chigh
     &        )then

            xmuf=xc1*scale
            xmur=xc2*scale

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)

           if (is.eq.1) then
           call pdft_Pqq(xa,xmuf,1,5,deltas,deltac,xa*xb*s,f1,f1tqq)
           call pdft_Pqq(xb,xmuf,1,5,deltas,deltac,xa*xb*s,f2,f2tqq)
           call sig_nlo2_qqb(xmuf,f1,f2,f1tqq,f2tqq,
     &          p1,p2,p3,p4,sig)                        ! qqb NLO2

           elseif (is.eq.2) then
           call pdft_Pqg(xa,xmuf,1,5,deltas,deltac,xa*xb*s,f1,f1tqg)
           call pdft_Pqg(xb,xmuf,1,5,deltas,deltac,xa*xb*s,f2,f2tqg)
! We will not use f1tgq,f2tgq in SM
           call pdft_Pgq(xa,xmuf,1,5,deltas,deltac,xa*xb*s,f1,f1tgq)
           call pdft_Pgq(xb,xmuf,1,5,deltas,deltac,xa*xb*s,f2,f2tgq)
           call sig_nlo2_qg(xmuf,f1,f2,f1tqg,f2tqg,f1tgq,f2tgq,
     &          p1,p2,p3,p4,sig)                        ! qg NLO2

           elseif (is.eq.3) then                        ! gg NLO2
           call pdft_Pgg(xa,xmuf,1,5,deltas,deltac,xa*xb*s,f1,f1tgg)
           call pdft_Pgg(xb,xmuf,1,5,deltas,deltac,xa*xb*s,f2,f2tgg)
           call sig_nlo2_gg(xmuf,f1,f2,f1tgg,f2tgg,
     &          p1,p2,p3,p4,sig)
           endif

           xnorm=hbarc2/16d0/pi/(s*xa*xb)
           wgt=xnorm*sig*vwgt*xjac
           fnlo2=wgt/vwgt/2.d0/xeps
           return              ! There is a factor of 1/2
         else                   ! from xeps interval
           fnlo2=0d0
           return
         endif
      ELSEIF(ipass.eq.0)THEN
         fnlo2=0d0
         return
      ENDIF

      end


