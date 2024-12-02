      function fnlo3(yy,vwgt)
      implicit double precision (a-h,o-z)
      REAL*8 INTIME,OUTTIME
      dimension yy(10)
      dimension f1(-6:6),f2(-6:6)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),q(0:3)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=0.3894d9)

c      data p1/ 0.5000000D+03, 0d0, 0d0, 0.5000000D+03/
c      data p2/ 0.5000000D+03, 0d0, 0d0,-0.5000000D+03/
c      data p3/0.4585788D+03,0.1694532D+03,0.3796537D+03,-0.1935025D+03/
c      data p4/0.3640666D+03,-0.1832987D+02,-0.3477043D+03,0.1063496D+03/
c      data p5/0.1773546D+03,-0.1511234D+03,-0.3194936D+02,0.8715287D+02/

      common/energy/s
      common/momenta5/p1,p2,p3,p4,p5
      common/bin/xq,xeps,xlow,xhigh
      common/cfbin/cf,cfeps,clow,chigh
      common/yfbin/yf,yfeps,ylow,yhigh
      common/xrange/xrlow,xrhigh
      common/factscale/xmuf
      common/param/aem,xmur,lambda
      common/xmcoeff/xc1,xc2
      common/isub/io,is
      common/angle/cststar
 
      call kinvar3(yy,xjac,xinvmass,y1,y2,Y,cst1,cst2,pt1,pt2,
     &     r34,r35,r45,Et5,QT)
      
      call cuts3(xinvmass,y1,y2,Y,cst1,cst2,pt1,pt2,
     &     r34,r35,r45,Et5,QT,ipass)
      

      xa =yy(1)
      xb =yy(2)
      v  =yy(3)

      scale = xinvmass
      pobl = cststar

      IF(ipass.eq.1)THEN

         if(   scale.gt.xlow .and. scale.lt.xhigh
c         if(   scale.gt.xrlow .and. scale.lt.xrhigh
c     &  .and. pobl.ge.clow .and. pobl.le.chigh
     &        )then

c        write(*,*)'Qmin, Qmax', xlow, xhigh
c            write(*,*)'io,is',io,is
c            write(*,*)'Q =', xinvmass

            xmuf=xc1*scale
            xmur=xc2*scale

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
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

            xnorm=hbarc2*v/64/4/4/(pi**4)
c            vwgt=1.0d0
            wgt=xjac*xnorm*sig*vwgt
            fnlo3=wgt/vwgt/2.d0/xeps
            return              !There is a factor of 1/2
         else                   !from xeps interval
            fnlo3=0d0
            return
         endif         
      ELSEIF(ipass.eq.0)THEN
         fnlo3=0d0
         return
      ENDIF

      end


