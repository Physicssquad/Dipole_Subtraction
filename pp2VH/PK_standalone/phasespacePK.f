! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
! +                            Two Body Phase space                                      +
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            

c---------------------------------------------------------------------
      subroutine kinvar2_PK(xa,xb,xc,Qmass,p1,p2,p3,p4)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      common/energy/s

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
        

c      Q2 = 2.0d0*dot(p3,p4)
c      Qmass = dsqrt(Q2)
      Qmass  = pobl(p1,p2,p3,p4)
      return
      end
c---------------------------------------------------------------------


      function pobl(p1,p2,p3,p4)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)


      xinvmass2 = 2.0d0*dot(p3,p4)
      pobl = dsqrt(xinvmass2)
      return
      end


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
! +                            Three Body Phase space                                      +
c! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
c
c      subroutine kinvar3(xx,xxjac,xinvmass,p1,p2,p3,p4,p5,unphy)
c      implicit double precision (a-h,o-z)
c      integer n4,unphy
c      parameter (pi=3.14159265358979d0)
c      dimension xx(10)
c      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),q(0:3)
c      dimension qa(0:3),qb(0:3)
c      dimension p4p(0:3),diff(0:3)
c      common/energy/s
c
c      am3=0d0
c      am4=0d0
c      am5=0d0
c
c      ! ThisHasToBeModifiedIn p-p Collision
c      xa= xx(1)
c      xb= xx(2)
c      xjac=1d0
c
c      s12=xa*xb*s
c      rsp = dsqrt(s12)
c      srs2= 0.5d0*dsqrt(s)
cc     incoming parton 4-vectors
c
c      p1(0)=xa*srs2
c      p1(1)=0d0
c      p1(2)=0d0
c      p1(3)=p1(0)
c
cc      write(*,*)"Random= ",xx
c      p2(0)=xb*srs2
c      p2(1)=0d0
c      p2(2)=0d0
c      p2(3)=-p2(0)
c
cc     outgoing parton 4-vectors
c
c      v=xx(3)
c      w=xx(4)
c      !N:  Mass of each particles taking part in interaction
c      unphy = 0
c      if (s12 .lt. (am3 + am4)**2 - am5**2)  unphy = 1
c
cc-------------------------------------------------------
cc 2 --> 3 body massive phase space parametrization based 
cc on the algorithm given in FormCalc version 4.1
cc-------------------------------------------------------
c      ct    = -1d0+2d0*v
c      xjac3 = 2d0
c      st   = dsqrt(1d0-ct*ct)
c
c      phi   = 2d0*pi*w
c      xjac4 = 2d0*pi
c      cphi  = dcos(phi)
c      sphi  = dsin(phi)
c
c      e5min = am5
c      e5max = 0.5d0*(rsp + (am5**2 - (am3+am4)**2)/rsp)
c      xjac5 = e5max - e5min
c      e5 = xjac5*xx(5) + e5min
c      p5m = dsqrt(dabs(e5**2 -am5**2))
c
c      p5t = e5
c      p5x = p5m*st
c      p5y = 0d0
c      p5z = p5m*ct
c
c      sigma = rsp - e5
c      tau   = sigma**2 - p5m**2
c      amp   = am3 + am4
c      amm   = am3 - am4
c      a1    = sigma*(tau + amp*amm)
c      b1    = (tau-amp**2)*(tau-amm**2)
c      a2    = p5m*dsqrt(b1)
c      e3min = 0.5d0/tau*(a1 - a2)
c      e3max = 0.5d0/tau*(a1 + a2)
c      xjac6 = e3max - e3min
c      e3    = xjac6*xx(6) + e3min
c      p3m   = dsqrt(dabs(e3**2 - am3**2))
c
c      e4    = rsp - e3 -e5
c      p4m   = dsqrt(dabs(e4**2 - am4**2))
c      czeta = (p4m**2 - p3m**2 - p5m**2)/(2d0*p3m*p5m)
c
cc    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cc      will check unphysical ps points
c      if (czeta .ge. 1.0d0) then
c        unphy = unphy+1 
c        goto 151
c      endif
cc    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c      zeta = dacos(czeta)
c      szeta = dsin(zeta)
c
cc      szeta = dsqrt(1.0d0-czeta*czeta)
cc      print*,'Values',e4,p4m,czeta,szeta
c
c      p3t   = e3
c      p3x   = p3m*(ct*cphi*szeta + st*czeta)
c      p3y   = p3m*(sphi*szeta)
c      p3z   = p3m*(ct*czeta - st*cphi*szeta)
c
c      p4t   = e4
c      p4x   = -p3x -p5x
c      p4y   = -p3y -p5y
c      p4z   = -p3z -p5z
c
c      beta=(xa-xb)/(xa+xb)
c      gamma=1d0/dsqrt(1d0-beta*beta)
c
cc      if(beta .ge. 1.0d0) then
cc      write(*,*)'beta =',beta
cc      endif
c
c
c      p3(0)=gamma*(p3t + beta*p3z)
c      p3(1)=p3x
c      p3(2)=p3y
c      p3(3)=gamma*(p3z + beta*p3t)
c
cc      p4(0)=gamma*(p4t + beta*p4z)
cc      p4(1)=p4x
cc      p4(2)=p4y
cc      p4(3)=gamma*(p4z + beta*p4t)
c
c      p5(0)=gamma*(p5t + beta*p5z)
c      p5(1)=p5x
c      p5(2)=p5y
c      p5(3)=gamma*(p5z + beta*p5t)
c
c
c      p4(0)=p1(0)+p2(0)-p3(0)-p5(0)
c      p4(1)=p1(1)+p2(1)-p3(1)-p5(1)
c      p4(2)=p1(2)+p2(2)-p3(2)-p5(2)
c      p4(3)=p1(3)+p2(3)-p3(3)-p5(3)
c      
cc     p3 + p4
cc      q(0) = p3(0) +p4(0)
cc      q(1) = p3(1) +p4(1)
cc      q(2) = p3(2) +p4(2)
cc      q(3) = p3(3) +p4(3)
c    
cc     xxjac
c      xxjac = xjac3*xjac4*xjac5*xjac6
c
c      xinvmass =dsqrt(2d0*dot(p3,p4))
c
c
c 151  return
c      end 
c 
cc---------------------------------------------------------------------
cc---------------------------------------------------------------------
c        subroutine resetmomenta(p1,p2,p3,p4,p5)
c        implicit none
c        double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
c        integer i
c        do i=0,3
c        p1(i)=0d0
c        p2(i)=0d0
c        p3(i)=0d0
c        p4(i)=0d0
c        p5(i)=0d0
c        enddo
c        return
c        end
cc---------------------------------------------------------------------
c         subroutine p1dtop2d_5(p1,p2,p3,p4,p5,p)
c         implicit double precision (a-h,o-z)
c
c         dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p(0:3,1:5)
c        do i=0,3
c         p(i,1)=p1(i)
c         p(i,2)=p2(i)
c         p(i,3)=p3(i)
c         p(i,4)=p4(i)
c         p(i,5)=p5(i)
c        enddo
c         end
cc---------------------------------------------------------------------
c
cc---------------------------------------------------------------------
c         subroutine p2dtop1d_5(p,p1,p2,p3,p4,p5)
c         implicit double precision (a-h,o-z)
c
c         dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p(0:3,1:5)
c        do i=0,3
c         p1(i)=p(i,1)
c         p2(i)=p(i,2)
c         p3(i)=p(i,3)
c         p4(i)=p(i,4)
c         p5(i)=p(i,5)
c        enddo
c         end
cc---------------------------------------------------------------------
cc---------------------------------------------------------------------
c         subroutine printmomenta(p1,p2,p3,p4)
c                 implicit double precision (a-h,o-z)
c                 dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
c                 write(*,*)"p1= ",p1
c                 write(*,*)"p2= ",p2
c                 write(*,*)"p3= ",p3
c                 write(*,*)"p4= ",p4
c         end
cc---------------------------------------------------------------------
cc---------------------------------------------------------------------
c         subroutine p2d_to_p1d_4(p,p1,p2,p3,p4)
c         implicit double precision (a-h,o-z)
c
c         dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p(0:3,1:4)
c        do i=0,3
c         p1(i)=p(i,1)
c         p2(i)=p(i,2)
c         p3(i)=p(i,3)
c         p4(i)=p(i,4)
c        enddo
c         end
cc---------------------------------------------------------------------
cc---------------------------------------------------------------------
c         subroutine p1d_to_p2d_4(p1,p2,p3,p4,p)
c         implicit double precision (a-h,o-z)
c
c         dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p(0:3,1:4)
c        do i=0,3
c         p(i,1)=p1(i)
c         p(i,2)=p2(i)
c         p(i,3)=p3(i)
c         p(i,4)=p4(i)
c        enddo
c         end
cc---------------------------------------------------------------------
c
