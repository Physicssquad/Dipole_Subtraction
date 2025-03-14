
c      subroutine kinvar2(xx,
c     &     xxinvmass,yy1,yy2,YY12,ccst1,ccst2,ppt1,ppt2) 
      subroutine kinvar2(xa,xb,xc,
     &     xxinvmass,yy1,yy2,YY12,ccst1,ccst2,ppt1,ppt2) 
      implicit double precision (a-h,o-z)
      dimension xx(10)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      dimension qa(0:3),qb(0:3)
      common/energy/s      
      common/momenta/p1,p2,p3,p4
      common/angle/cststar

c      xa=xx(1)
c      xb=xx(2)
c      v=xx(3)
      v = xc
      omv=1d0-v     
 
      srs2=0.5*dsqrt(s)

c     incoming parton 4-vectors
      p1(0)=srs2*xa
      p1(1)=0d0
      p1(2)=0d0
      p1(3)=p1(0)

      p2(0)=srs2*xb
      p2(1)=0d0
      p2(2)=0d0
      p2(3)=-p2(0)

c     outgoing parton 4-vectors
      p3(0)=srs2*(xa*v+xb*omv)
      p3(1)=dsqrt(s*xa*xb*v*omv)
      p3(2)=0d0
      p3(3)=srs2*(xa*v-xb*omv)

      p4(0)=p1(0)+p2(0)-p3(0)
      p4(1)=p1(1)+p2(1)-p3(1)
      p4(2)=p1(2)+p2(2)-p3(2)
      p4(3)=p1(3)+p2(3)-p3(3)

c     p3 + p4
      q(0) = p3(0) +p4(0)
      q(1) = p3(1) +p4(1)
      q(2) = p3(2) +p4(2)
      q(3) = p3(3) +p4(3)

c     invariant mass of diphoton pair.
      s34    = 2*dot(p3,p4)
      xxinvmass= dsqrt(s34)

c     yy1: rapidity of photons 3 
      yrpda  = (p3(0)+p3(3))/(p3(0)-p3(3))
      yy1    = 0.5*dlog(yrpda)

c     yy2: rapiditiy of photon 4
      yrpdb  = (p4(0)+p4(3))/(p4(0)-p4(3))
      yy2    = 0.5*dlog(yrpdb)
      
c     rapidity YY
      p2q = dot(p2,q)
      p1q = dot(p1,q)
      rr  = xa*p2q/(xb*p1q)
      YY12  = dlog(rr)/2.d0

c     cos-theta
      ccst1=p3(3)/p3(0)
      ccst2=p4(3)/p4(0)

c     ppt1 ppt2
      ppt1=dsqrt(p3(1)**2+p3(2)**2)
      ppt2=dsqrt(p4(1)**2+p4(2)**2)

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

      return
      end


