c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
c Two body phase space
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      

      subroutine kinvar2(xx,xxinvmass,p1,p2,p3,p4)
      implicit double precision (a-h,o-z)
      dimension xx(4)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      dimension qa(0:3),qb(0:3)
      common/angle/cststar

      s=1000d0**2
      v=xx(1)
      omv=1d0-v
      srs2=0.5*dsqrt(s)

c     incoming parton 4-vectors
      p1(0)=srs2
      p1(1)=0d0
      p1(2)=0d0
      p1(3)=p1(0)

      p2(0)=srs2
      p2(1)=0d0
      p2(2)=0d0
      p2(3)=-p2(0)

c     outgoing parton 4-vectors
      p3(0)=srs2*(v+omv)
      p3(1)=dsqrt(s*v*omv)
      p3(2)=0d0
      p3(3)=srs2*(v-omv)

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
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
c Three Body Phase space
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      

      subroutine kinvar3(xx,xxjac,p1,p2,p3,p4,p5)
      implicit double precision (a-h,o-z)
      dimension xx(10)
      parameter (pi=3.14159265358979d0)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),q(0:3)
      dimension qa(0:3),qb(0:3)

      s=1000d0**2
      am3=0d0
      am4=0d0
      am5=0d0

      ! ThisHasToBeModifiedIn p-p Collision
      xa=1.d0! xx(1)
      xb=1.d0! xx(2)
      xjac=1d0
      srs2=0.5*dsqrt(s)
c      srs2=0.5*s
c      write(*,*)'Values on xx=',xx

c     incoming parton 4-vectors

      p1(0)=xa*srs2
      p1(1)=0d0
      p1(2)=0d0
      p1(3)=p1(0)

c      write(*,*)"Random= ",xx
      p2(0)=xb*srs2
      p2(1)=0d0
      p2(2)=0d0
      p2(3)=-p2(0)

c     outgoing parton 4-vectors
      s12=xa*xb*s
      rsp = dsqrt(s12)

      v=xx(1)
      w=xx(2)
      !N:  Mass of each particles taking part in interaction
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

      e5min = am5
      e5max = 0.5d0*(rsp + (am5**2 - (am3+am4)**2)/rsp)
      xjac5 = e5max - e5min
      e5 = xjac5*xx(3) + e5min
      p5m = dsqrt(dabs(e5**2 -am5**2))

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
      e3    = xjac6*xx(4) + e3min
      p3m   = dsqrt(dabs(e3**2 - am3**2))
c      write(*,*)'Sigma, tau, a1, b1, a2=',sigma, tau, a1,b1,a2
c      write(*,*)'am3,am4 =',am3,am4

      e4    = rsp - e3 -e5
      p4m   = dsqrt(dabs(e4**2 - am4**2))
      czeta = (p4m**2 - p3m**2 - p5m**2)/(2.0d0*p3m*p5m)
c      if (dabs(czeta) .gt. 1.0d0) then
c      szeta = dsqrt(dabs(1.0d0-czeta*czeta))
      
      if (dabs(czeta) .ge. 1.0d0) then
       call resetmomenta(p1,p2,p3,p4,p5)
      goto 151
      endif       
      szeta = dsqrt(1.0d0-czeta*czeta)
c      print*,'Values',e4,p4m,czeta,szeta
c      if (dabs(czeta) .gt. 1.0d0) then
c      write(*,*)'Cos (zeta) =', czeta,p4m,p3m,p5m

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

c     p3 + p4
      q(0) = p3(0) +p4(0)
      q(1) = p3(1) +p4(1)
      q(2) = p3(2) +p4(2)
      q(3) = p3(3) +p4(3)
    
c     xxjac
      xxjac = xjac3*xjac4*xjac5*xjac6
c      write(*,*)'Jacobian =', xxjac
c      endif

 151  return
      end 
c---------------------------------------------------------------------
        subroutine resetmomenta(p1,p2,p3,p4,p5)
        implicit none
        double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
        integer i
        do i=0,3
        p1(i)=0d0
        p2(i)=0d0
        p3(i)=0d0
        p4(i)=0d0
        p5(i)=0d0
        enddo
        end

c---------------------------------------------------------------------
 
