! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
! +                            Two Body Phase space                                      +
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
      subroutine kinvar2(xx,xxinvmass,p1,p2,p3,p4)
      implicit double precision (a-h,o-z)
      dimension xx(10)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      dimension qa(0:3),qb(0:3)
      common/energy/s

      xa=xx(1)
      xb=xx(2)
      v=xx(3)
      omv=1d0-v
      

c      s=s*xa*xb 
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
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
! +                            Three Body Phase space                                      +
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      

      subroutine kinvar3(xx,xxjac,xinvmass,p1,p2,p3,p4,p5,unphy)
      implicit double precision (a-h,o-z)
      integer n4,unphy
      parameter (pi=3.14159265358979d0)
      dimension xx(10)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),q(0:3)
      dimension qa(0:3),qb(0:3)
      dimension p4p(0:3),diff(0:3)
      common/energy/s

      am3=0d0
      am4=0d0
      am5=0d0

      ! ThisHasToBeModifiedIn p-p Collision
      xa= xx(1)
      xb= xx(2)
      xjac=1d0

      s12=xa*xb*s
      rsp = dsqrt(s12)
      srs2= 0.5d0*dsqrt(s)
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

      v=xx(3)
      w=xx(4)
      !N:  Mass of each particles taking part in interaction
      unphy = 0
      if (s12 .lt. (am3 + am4)**2 - am5**2)  unphy = 1

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
      e5 = xjac5*xx(5) + e5min
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
      e3    = xjac6*xx(6) + e3min
      p3m   = dsqrt(dabs(e3**2 - am3**2))

      e4    = rsp - e3 -e5
      p4m   = dsqrt(dabs(e4**2 - am4**2))
      czeta = (p4m**2 - p3m**2 - p5m**2)/(2d0*p3m*p5m)

c    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c      will check unphysical ps points
      if (czeta .ge. 1.0d0) then
        unphy = unphy+1 
        goto 151
      endif
c    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      zeta = dacos(czeta)
      szeta = dsin(zeta)

c      szeta = dsqrt(1.0d0-czeta*czeta)
c      print*,'Values',e4,p4m,czeta,szeta

      p3t   = e3
      p3x   = p3m*(ct*cphi*szeta + st*czeta)
      p3y   = p3m*(sphi*szeta)
      p3z   = p3m*(ct*czeta - st*cphi*szeta)

      p4t   = e4
      p4x   = -p3x -p5x
      p4y   = -p3y -p5y
      p4z   = -p3z -p5z

      beta=(xa-xb)/(xa+xb)
      gamma=1d0/dsqrt(1d0-beta*beta)

c	beta = 0d0
c	gamma = 1d0

c      if(beta .ge. 1.0d0) then
c      write(*,*)'beta =',beta
c      endif


      p3(0)=gamma*(p3t + beta*p3z)
      p3(1)=p3x
      p3(2)=p3y
      p3(3)=gamma*(p3z + beta*p3t)

c      p4(0)=gamma*(p4t + beta*p4z)
c      p4(1)=p4x
c      p4(2)=p4y
c      p4(3)=gamma*(p4z + beta*p4t)

      p5(0)=gamma*(p5t + beta*p5z)
      p5(1)=p5x
      p5(2)=p5y
      p5(3)=gamma*(p5z + beta*p5t)


      p4(0)=p1(0)+p2(0)-p3(0)-p5(0)
      p4(1)=p1(1)+p2(1)-p3(1)-p5(1)
      p4(2)=p1(2)+p2(2)-p3(2)-p5(2)
      p4(3)=p1(3)+p2(3)-p3(3)-p5(3)
      
c     p3 + p4
c      q(0) = p3(0) +p4(0)
c      q(1) = p3(1) +p4(1)
c      q(2) = p3(2) +p4(2)
c      q(3) = p3(3) +p4(3)
    
c     xxjac
      xxjac = xjac3*xjac4*xjac5*xjac6

      xinvmass =dsqrt(2d0*dot(p3,p4))


 151  return
      end 
 
c---------------------------------------------------------------------
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
        return
        end
c---------------------------------------------------------------------
         subroutine p1dtop2d_5(p1,p2,p3,p4,p5,p)
         implicit double precision (a-h,o-z)

         dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p(0:3,1:5)
        do i=0,3
         p(i,1)=p1(i)
         p(i,2)=p2(i)
         p(i,3)=p3(i)
         p(i,4)=p4(i)
         p(i,5)=p5(i)
        enddo
         end
c---------------------------------------------------------------------

c---------------------------------------------------------------------
         subroutine p2dtop1d_5(p,p1,p2,p3,p4,p5)
         implicit double precision (a-h,o-z)

         dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p(0:3,1:5)
        do i=0,3
         p1(i)=p(i,1)
         p2(i)=p(i,2)
         p3(i)=p(i,3)
         p4(i)=p(i,4)
         p5(i)=p(i,5)
        enddo
         end
c---------------------------------------------------------------------
c---------------------------------------------------------------------
         subroutine printmomenta(p1,p2,p3,p4)
                 implicit double precision (a-h,o-z)
                 dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
                 write(*,*)"p1= ",p1
                 write(*,*)"p2= ",p2
                 write(*,*)"p3= ",p3
                 write(*,*)"p4= ",p4
         end
c---------------------------------------------------------------------
c---------------------------------------------------------------------
         subroutine p2d_to_p1d_4(p,p1,p2,p3,p4)
         implicit double precision (a-h,o-z)

         dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p(0:3,1:4)
        do i=0,3
         p1(i)=p(i,1)
         p2(i)=p(i,2)
         p3(i)=p(i,3)
         p4(i)=p(i,4)
        enddo
         end
c---------------------------------------------------------------------
c---------------------------------------------------------------------
         subroutine p1d_to_p2d_4(p1,p2,p3,p4,p)
         implicit double precision (a-h,o-z)

         dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p(0:3,1:4)
        do i=0,3
         p(i,1)=p1(i)
         p(i,2)=p2(i)
         p(i,3)=p3(i)
         p(i,4)=p4(i)
        enddo
         end
c---------------------------------------------------------------------

