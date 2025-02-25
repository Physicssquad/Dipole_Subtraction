
      subroutine kinvar3(xx,xxjac,
     &     xxinvmass,yy1,yy2,YY12,ccst1,ccst2,ppt1,ppt2,
     &     rr34,rr35,rr45,EEt5,QT) 
      implicit double precision (a-h,o-z)
      dimension xx(10)
      parameter (pi=3.14159265358979d0)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),q(0:3)
      dimension qa(0:3),qb(0:3)
      common/energy/s
      common/momenta5/p1,p2,p3,p4,p5
      common/angle/cststar
    
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

      v=xx(3)
      w=xx(4)

      s45=s12*v*(1d0-w)

      ct1=-1d0+2d0*xx(5)
      xjac=2d0*xjac
      st1=dsqrt(1d0-ct1*ct1)
      phi=2d0*pi*xx(6)
      xjac=xjac*2d0*pi

      ct2=dcos(phi)
      st2=dsin(phi)
      ct=1d0-2d0*(1d0-v)/(1d0-v+v*w)
      st=dsqrt(1d0-ct*ct)

      beta=(xa-xb)/(xa+xb)
      gamma=1d0/dsqrt(1d0-beta*beta)
      e3s=(s12-s45)/(2d0*dsqrt(s12))

      p3(0)=e3s*gamma*(1d0+beta*ct)
      p3(1)=e3s*st
      p3(2)=0d0
      p3(3)=e3s*gamma*(ct+beta)

      e4s=(s12+s45)/(2d0*dsqrt(s12))
      e4=0.5d0*(e4s-e3s*ct1)
      p4z=0.5d0*(e4s*ct1-e3s)
      r45=0.5d0*dsqrt(s45)
      p41=r45*st1*ct2*ct+p4z*st
      p42=r45*st1*st2
      p43=-r45*st1*ct2*st+p4z*ct

      p4(0)=gamma*(e4+beta*p43)
      p4(1)=p41
      p4(2)=p42
      p4(3)=gamma*(p43+beta*e4)

      p5(0)=p1(0)+p2(0)-p3(0)-p4(0)
      p5(1)=p1(1)+p2(1)-p3(1)-p4(1)
      p5(2)=p1(2)+p2(2)-p3(2)-p4(2)
      p5(3)=p1(3)+p2(3)-p3(3)-p4(3)

c     p3 + p4
      q(0) = p3(0) +p4(0)
      q(1) = p3(1) +p4(1)
      q(2) = p3(2) +p4(2)
      q(3) = p3(3) +p4(3)
            
c     QT      
      QT=dsqrt(q(1)**2+q(2)**2)
      
c     xxjac
      xxjac=xjac

c     invariant mass of diphoton pair.
      s34       = 2*dot(p3,p4)
      xxinvmass = dsqrt(s34)            

c     yy1: rapidity of photons 3 
      yy1 = rapid(p3)

c     yy2: rapiditiy of photon 4
      yy2 = rapid(p4)
      
c     rapidity YY12
      p2q = dot(p2,q)
      p1q = dot(p1,q)
      rr  = xa*p2q/(xb*p1q)
      YY12  = dlog(rr)/2.d0

c     cos-theta 
      ccst1=p3(3)/p3(0)
      ccst2=p4(3)/p4(0)            

c     ppt1 ppt2
      ppt1 = eperp(p3)
      ppt2 = eperp(p4)

      qa(0) =  p3(0) - p4(0)
      qa(1) =  p3(1) - p4(1)
      qa(2) =  p3(2) - p4(2)
      qa(3) =  p3(3) - p4(3)

      qb(0) =  p3(0) + p4(0)
      qb(1) =  p3(1) + p4(1)
      qb(2) =  p3(2) + p4(2)
      qb(3) =  p3(3) + p4(3)

      p1qa = dot(p1,qa)
      p1qb = dot(p1,qb)

      cststar = p1qa/p1qb

c     cone
      rr34 = rjet(p3,p4)
      rr35 = rjet(p3,p5)
      rr45 = rjet(p4,p5)
      
c     transverse energy
      EEt5 = eperp(p5)  ! for a massless particle Et=pt

      return
      end
 
