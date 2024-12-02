
c      subroutine kinvar3(xx,xxjac,
c     &     xxinvmass,yy1,yy2,yy3,YY12,ccst1,ccst2,ppt1,ppt2,
c     &     rr34,rr35,rr45,EEt5,QT) 

      subroutine kinvar3(xx,xxjac,xxinvmass,p1,p2,p3,p4,p5) 
      implicit double precision (a-h,o-z)
      dimension xx(10)
      parameter (pi=3.14159265358979d0)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),q(0:3)
      dimension qa(0:3),qb(0:3)
      common/energy/s
c      common/momenta5/p1,p2,p3,p4,p5
      common/amf/am3,am4,am5
      common/bin/xq,xeps,xlow,xhigh
      common/angle/cststar
 
      xxinvmass = 0.0d0

      xa=xx(1)
      xb=xx(2)
      xjac=1d0
      srs2=0.5*dsqrt(s)
      
c     incoming parton 4-vectors
      p1(0)=xa*srs2
      p1(1)=0d0
      p1(2)=0d0
      p1(3)=p1(0)

      p2(0)=xb*srs2
      p2(1)=0d0
      p2(2)=0d0
      p2(3)=-p2(0)

c     outgoing parton 4-vectors
      s12=xa*xb*s
      rsp = dsqrt(s12)

c      p1cm(0)=0.5d0*rsp
c      p1cm(1)=0d0
c      p1cm(2)=0d0
c      p1cm(3)=p1cm(0)

c      p2cm(0)=0.5d0*rsp
c      p2cm(1)=0d0
c      p2cm(2)=0d0
c      p2cm(3)=-p2cm(0)


      v=xx(3)
      w=xx(4)

      if (s12. lt. (am3 + am4)**2 - am5**2) goto 151

c-------------------------------------------------------
c 2 --> 3 body massive phase space parametrization based 
c on the algorithm given in FormCalc version 4.1
c-------------------------------------------------------

      ct    = -1d0+2d0*v
      xjac3 = 2d0
      st   = dsqrt(1d0-ct*ct)

      phi   = 2d0*pi*w
      xjac4 = 2d0*pi
      cphi  = dcos(phi)
      sphi  = dsin(phi)

c      soft=0.5d0*deltas*rsp
c      e5min = soft
c       e5min = am5
c      e5max = 0.5d0*(rsp + (am5**2 - (am3+am4)**2)/rsp) 
c
c      xjac5 = e5max - e5min
c      e5 = xjac5*xx(5) + e5min

c       Modified by Chinmoy 
      e5 = (s12-xq*xq+am5**2)/2.0d0/rsp
      p5m = dsqrt(dabs(e5**2 -am5**2))

      if (rsp .le. xq) then
      write(*,*)'rsp, Q =',rsp, xq
      endif

      if (e5 .le.0 .or. e5 .ge. rsp) then
      write(*,*)'E5cm =',e5
      endif

      p5t = e5
      p5x = p5m*st
      p5y = 0d0
      p5z = p5m*ct

      sigma = rsp - e5
      tau   = sigma**2 - p5m**2
      amp   = am3 + am4
      amm   = am3 - am4
      a1    = sigma*(tau + amp*amm)
      b1    = (tau-amp**2)*(tau-amm**2)
      a2    = p5m*dsqrt(b1)
      e3min = 0.5d0/tau*(a1 - a2)
      e3max = 0.5d0/tau*(a1 + a2)
      xjac6 = e3max - e3min
      e3    = xjac6*xx(5) + e3min 
      p3m   = dsqrt(dabs(e3**2 - am3**2))

c      write(*,*)'Sigma, tau, a1, b1, a2=',sigma, tau, a1,b1,a2

      e4    = rsp - e3 -e5
      p4m   = dsqrt(dabs(e4**2 - am4**2))
      czeta = (p4m**2 - p3m**2 - p5m**2)/(2.0d0*p3m*p5m)
      szeta = dsqrt(dabs(1.0d0-czeta*czeta))

c      if (dabs(czeta) .gt. 1.0d0) then
c      write(*,*)'Cos (zeta) =', czeta
c      endif

      p3t   = e3
      p3x   = p3m*(ct*cphi*szeta + st*czeta)
      p3y   = p3m*(sphi*szeta)
      p3z   = p3m*(ct*czeta - st*cphi*szeta)

      beta=(xa-xb)/(xa+xb)
      gamma=1d0/dsqrt(1d0-beta*beta)

      p3(0)=gamma*(p3t + beta*p3z)
      p3(1)=p3x
      p3(2)=p3y
      p3(3)=gamma*(p3z + beta*p3t)

      p5(0)=gamma*(p5t + beta*p5z)
      p5(1)=p5x
      p5(2)=p5y
      p5(3)=gamma*(p5z + beta*p5t)

      p4(0)=p1(0)+p2(0)-p3(0)-p5(0)
      p4(1)=p1(1)+p2(1)-p3(1)-p5(1)
      p4(2)=p1(2)+p2(2)-p3(2)-p5(2)
      p4(3)=p1(3)+p2(3)-p3(3)-p5(3)

c     xxjac
c      xxjac = xjac3*xjac4*xjac5*xjac6
      xxjac = xjac3*xjac4*xjac6
c      write(*,*)'Jacobian =', xxjac

c     invariant mass of diphoton pair.
      s34       = am3**2 + am4**2 + 2*dot(p3,p4)
      xxinvmass = dsqrt(s34)            
c      write(*,*) 'Qmass =', xxinvmass

 151  return
      end
 

