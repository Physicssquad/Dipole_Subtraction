      function fnlo3(yy,vwgt)
      implicit double precision (a-h,o-z)
      REAL*8 INTIME,OUTTIME
      dimension xa(10),yy(10)
      dimension f1(-6:6),f2(-6:6)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),q(0:3)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=0.3894d9)

      common/energy/s
c      common/momenta5/p1,p2,p3,p4,p5
      common/bin/xq,xeps,xlow,xhigh
      common/cfbin/cf,cfeps,clow,chigh
      common/yfbin/yf,yfeps,ylow,yhigh
      common/xrange/xrlow,xrhigh
      common/factscale/xmuf
      common/param/aem,xmur,lambda
      common/xmcoeff/xc1,xc2
      common/isub/io,is
      common/angle/cststar
      common/amf/am3,am4,am5
 
      xa(1) = yy(1)
      xa(2) = yy(2)
      xa(3) = yy(3)
      xa(4) = yy(4)
      xa(5) = yy(5)
c      xa(6) = yy(6)
      sp = xa(1)*xa(2)*S
      rsp = dsqrt(sp)

      if (rsp .ge. xq) then
      ipass = 1
c      write(*,*)'rsp =', rsp
      else
      ipass = 0
      endif

      IF(ipass.eq.1)THEN

      call kinvar3(xa,xjac,xinvmass,p1,p2,p3,p4,p5)

      call cuts3(p1,p2,p3,p4,p5,ipass1)      
     
c      ipass1 = 1
      if (ipass1 .eq.1) then

            scale = xinvmass
            xmuf=xc1*scale
            xmur=xc2*scale

c            write(*,*)'x1,x2=',xa(1),xa(2)

            call pdf(xa(1),xmuf,f1)
            call pdf(xa(2),xmuf,f2)
            CALL CPU_TIME(INTIME)
            call sig_nlo3(f1,f2,p1,p2,p3,p4,p5,sig)
            CALL CPU_TIME(OUTTIME)
c            write(*,*)'TIME = ',OUTTIME-INTIME

c            write(*,*)'p1 = ',p1
c            write(*,*)'p2 = ',p2
c            write(*,*)'p3 = ',p3
c            write(*,*)'p4 = ',p4
c            write(*,*)'p5 = ',p5
c            write(*,*)'sig = ',sig

            p_in  = 0.5d0*rsp
            flux  = 4.0d0*p_in*rsp

            xnorm=hbarc2/8.0d0/(2.0d0*pi)**4/flux
c            xnorm=hbarc2*v/64/4/4/(pi**4)
            wgt=xjac*xnorm*sig*vwgt
c            wgt=xjac*xnorm*sig*vwgt
c            fnlo3=wgt/vwgt/2.d0/xeps
c          sp - 2*rsp*e5 = Q^2
c          2*rsp*de5 = 2*Q*dQ
c          rsp/Q*de5 = dQ
c          1/dQ = 1/de5 *Q/rsp 
            wgt = wgt*xq/rsp
            fnlo3=wgt/vwgt
            return              !There is a factor of 1/2
         else                   !from xeps interval
c            write(*,*)'Cuts not satisfied'
            fnlo3=0d0
            return
         endif         
      ELSEIF(ipass.eq.0)THEN
         fnlo3=0d0
         return
      ENDIF

      end


