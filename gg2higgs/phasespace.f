      subroutine kinvar1(xa,xb,p1,p2,p3)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3)
      common/energy/s
      common/amass/am1,am2,amH,am4,am5

      srs2 = 0.5 * dsqrt(s)

c     incoming parton 4-vectors
      p1(0) = srs2 * xa
      p1(1) = 0d0
      p1(2) = 0d0
      p1(3) = p1(0)

      p2(0) = srs2 * xb
      p2(1) = 0d0
      p2(2) = 0d0
      p2(3) = -p2(0)

c     total 4-momentum of the system

      p3(0) = p1(0) + p2(0)
      p3(1) = p1(1) + p2(1)
      p3(2) = p1(2) + p2(2)
      p3(3) = p1(3) + p2(3)

!! Initial partons are already in lab frame. No need to boost.

C   [ Boost to Lab Frame ]

c      e_total  = p1(0) + p2(0)
c      px_total = p1(1) + p2(1)
c      py_total = p1(2) + p2(2)
c      pz_total = p1(3) + p2(3)
c
cc     calculate the boost in lab frame 
cc      beta = 0d0
cc      beta = (xa - xb) / (xa + xb)
cc      gama = 1d0 / sqrt(1d0 - beta**2)
c
cc     calculate momentum of the Higgs boson
c      p_h = sqrt(px_total**2 + py_total**2 + pz_total**2)
c
cc     Energy of the Higgs boson with mass
c      e_higgs = sqrt(p_h**2 + amh**2)
c
cc     Set the 4-vector for the outgoing Higgs boson
c      p3(0) = gama * (e_higgs + beta * pz_total)
c      p3(1) = px_total
c      p3(2) = py_total
c      p3(3) = gama * (pz_total + beta * e_higgs)

c       Q2 = dot(p3,p3)
c       Q = dsqrt(Q2)

c       write(*,*)'Q =', Q

      return
      end

c---------------------------------------------------------------o
C          C~~~~~~[ MASSLESS CASE ]~~~~~~~~C
      subroutine kinvar2(xx,xxinvmass,p1,p2,p3,p4)
      implicit double precision (a-h,o-z)
      dimension xx(6)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      dimension qa(0:3),qb(0:3)
      common/energy/s
      common/amass/am1,am2,am3,am4,am5

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
      s34    = 2d0*dot(p3,p4)
      xxinvmass= dsqrt(s34)


      return
      end
c---------------------------------------------------------------------
c   Try

c      SUBROUTINE GEN_PHASE_SPACE( X1,X2,X4, AM3,AM4, RS, P1,P2,P3,P4 )
      subroutine kinvar2_type_2(x1,x2,x4,xinvmass,p1,p2,p3,p4)
  
      DOUBLE PRECISION X1, X2, X4, RS
      DOUBLE PRECISION P1(0:3), P2(0:3), P3(0:3), P4(0:3)
      DOUBLE PRECISION XA, XB, V, OMV
      DOUBLE PRECISION SP, RSP, SRS2,s34, S
      DOUBLE PRECISION E3, E4, XLAM, PCM, PF, CT, ST
      DOUBLE PRECISION PX3, PY3, PZ3
      DOUBLE PRECISION BETA, GAMMA,am1,am2,am3,am4,am5
      double precision  dot

      common/energy/s
      common/amass/am1,am2,am3,am4,am5
      common/checks/ct,pcm,xlam


      ! Assign random inputs and auxiliary variables
      XA = X1
      XB = X2
      V = X4
      OMV = 1D0 - V
      RS = dsqrt(s)
      ! Compute invariant mass square and energy scale
      SP = XA * XB * S
      RSP = DSQRT(SP)
      SRS2 = 0.5D0 * DSQRT(RS * RS)

! Incoming parton 4-vectors in COM frame
      P1(0) = SRS2 * XA
      P1(1) = 0D0
      P1(2) = 0D0
      P1(3) = P1(0)

      P2(0) = SRS2 * XB
      P2(1) = 0D0
      P2(2) = 0D0
      P2(3) = -P2(0)

! Handle final-state particles
! Parametrization in the COM frame

      E3 = 0.5D0 * (SP + AM3**2 - AM4**2) / RSP
c      E4 = 0.5D0 * (SP + AM4**2 - AM3**2) / RSP

      XLAM = (AM4**2 - (RSP + AM3)**2) * (AM4**2 - (RSP - AM3)**2)
      PCM = 0.5D0 * DSQRT(XLAM)/RSP

      PF = DSQRT(E3**2 - AM3**2)
      CT = 2D0 * V - 1D0
      ST = DSQRT(1D0 - CT**2)

      PX3 = PF * ST 
      PY3 = 0D0 
      PZ3 = PF * CT 

      !Boost to the lab frame
      BETA = (XB - XA) / (XA + XB)
      GAMMA = 1D0 / DSQRT(1D0 - BETA * BETA)

c      BETA = 0d0
c      gamma = 1d0

      P3(0) = GAMMA * (E3 - BETA * PZ3)
      P3(1) = PX3
      P3(2) = PY3
      P3(3) = GAMMA * (PZ3 - BETA * E3)

c       print*,"p3",e3,p3

      ! Conservation of momentum for the fourth particle
c      P4(0) = P1(0) + P2(0) - P3(0)
c      P4(1) = P1(1) + P2(1) - P3(1)
c      P4(2) = P1(2) + P2(2) - P3(2)
c      P4(3) = P1(3) + P2(3) - P3(3)

c     invariant mass of diphoton pair.
c      s34    = 2d0*dot(p3,p4)
c      xxinvmass= dsqrt(s34)

      RETURN
      END
c---------------------------------------------------------------o

C          C~~~~~~[ MASSIVE  CASE ]~~~~~~~~C

      subroutine kinvar2_type_1(x1,x2,x4,xinvmass,p1,p2,p3,p4)
      implicit double precision (a-h,o-z)
      dimension xx(10)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      dimension qa(0:3),qb(0:3)
      COMMON/energy/s
      common/amass/am1,am2,am3,am4,am5
      common/checks/ct

      xa=x1
      xb=x2
      v=x4
      omv=1d0-v     

c      s = RS*RS
      sp = xa*xb*s
      rsp = dsqrt(sp)

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

c  For massive case where m3 and m4 are non-zero
c  ---------------------------------------------
c        Parametrization in the c.o.m. frame
c  ---------------------------------------------
      xmd = (am3**2 - am4**2)
       e3 = 0.5d0*(rsp + xmd/rsp)
       pf = dsqrt(e3**2 - am3**2)
       ct = 2d0*v-1
       st = dsqrt(1d0 - ct**2d0)
      px3 = pf*st
      py3 = 0d0
      pz3 = pf*ct

c---------------------------------------------------
c     e4cm = (sp - am3**2.0d0)/2.0d0/rsp
c---------------------------------------------------

      beta=(xb-xa)/(xa+xb)
      gamma=1d0/dsqrt(1d0-beta*beta)

      p3(0) = gamma*(e3 - beta*pz3)
      p3(1) = px3
      p3(2) = py3
      p3(3) = gamma*(pz3 - beta*e3)

      p4(0)=p1(0)+p2(0)-p3(0)
      p4(1)=p1(1)+p2(1)-p3(1)
      p4(2)=p1(2)+p2(2)-p3(2)
      p4(3)=p1(3)+p2(3)-p3(3)

      xinvmass = dsqrt(dot(p3,p3))

      return
      end



c---------------------------------------------------------------o
c
cC          C~~~~~~[ MASSIVE  CASE ]~~~~~~~~C
c
c      subroutine kinvar2_type_2(x1,x2,x4,xinvmass,p1,p2,p3,p4)
c      implicit double precision (a-h,o-z)
c      dimension xx(10)
c      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
c      dimension qa(0:3),qb(0:3)
c      COMMON/energy/s
c      common/amass/am1,am2,am3,am4,am5
c      common/checks/ct
c
c      xa=x1
c      xb=x2
c      v=x4
c      omv=1d0-v     
c
cc      s = RS*RS
c      sp = xa*xb*s
c      rsp = dsqrt(sp)
c
c      srs2=0.5*dsqrt(s)
c
cc     incoming parton 4-vectors
c      p1(0)=srs2*xa
c      p1(1)=0d0
c      p1(2)=0d0
c      p1(3)=p1(0)
c
c      p2(0)=srs2*xb
c      p2(1)=0d0
c      p2(2)=0d0
c      p2(3)=-p2(0)
c
cc For massive case where m3 and m4 are non-zero
cc ---------------------------------------------
cc       Parametrization in the c.o.m. frame
cc ---------------------------------------------
cc      xmd = (am3**2 - am4**2)
c      e3 = 0.5d0*(sp+am3*2 - am4**2)/rsp
c      e4 = 0.5d0*(sp+am4*2 - am3**2)/rsp
c     xlam = (am4**2 - (rsp + am3)**2)*(am4**2 - (rsp - am3)**2)
c      pcm = 0.5d0*dsqrt(xlam)
c      pf = dsqrt(e3**2 - am3**2)
c      ct = 2d0*v-1
c      st = dsqrt(1d0 - ct**2d0)
c      px3 = pf*st
c      py3 = 0d0
c      pz3 = pf*ct
c        
cc ---------------------------------------------
cc       Boost to lab frame
cc ---------------------------------------------
c      beta=(xb-xa)/(xa+xb)
c      gamma=1d0/dsqrt(1d0-beta*beta)
cc      beta=0d0
cc      gamma = 1d0
c
c      p3(0) = gamma*(e3 - beta*pz3)
c      p3(1) = px3
c      p3(2) = py3
c      p3(3) = gamma*(pz3 - beta*e3)
c
c      p4(0)=p1(0)+p2(0)-p3(0)
c      p4(1)=p1(1)+p2(1)-p3(1)
c      p4(2)=p1(2)+p2(2)-p3(2)
c      p4(3)=p1(3)+p2(3)-p3(3)
cc	print*,"p3-PS:",p3
cc	print*,"e3:",e3
c
c      xinvmass = dsqrt(dot(p3,p3))
c
c      return
c      end


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
! +                            Three Body Phase space                                      +
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      

      subroutine kinvar3(xx,xxjac,xinvmass,p1,p2,p3,p4,p5,unphy)
      implicit double precision (a-h,o-z)
      integer n4,unphy
      parameter (pi=3.14159265358979d0)
      dimension xx(6)
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

c      if (s12 .lt. (am3 + am4)**2 - am5**2)  unphy = 1
      if (rsp .lt. 40.0d0)  unphy = 1

c-------------------------------------------------------
c 2 --> 3 body massive phase space parametrization based 
c on the algorithm given in FormCalc version 4.1
c-------------------------------------------------------
      xjac3 = 2d0
      ct    = -1d0+xjac3*v
c      st   = dsqrt(1d0-ct*ct)
      theta = dacos(ct)
      st   = dsin(theta)

      if( theta .ge. pi) then
c      write(*,*)'theta =' ,theta
      unphy = unphy + 1
      endif

c      phi   = 2d0*pi*w
      xjac4 = 2d0*pi
      phi   = xjac4*w + 0.0d0
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
      endif
c    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (dabs(cphi) .ge. 1.0d0) then
      unphy = unphy+1 
      endif
c    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (dabs(sphi) .ge. 1.0d0) then
      unphy = unphy+1 
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

      p5(0)=gamma*(p5t - beta*p5z)
      p5(1)=p5x
      p5(2)=p5y
      p5(3)=gamma*(p5z - beta*p5t)


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
c----------------------------------------------------------------
c----------------------------------------------------------------
         subroutine p1d_to_p2d_3(p1,p2,p3,p)
         implicit double precision (a-h,o-z)

         dimension p1(0:3),p2(0:3),p3(0:3),p(0:3,1:3)
        do i=0,3
         p(i,1)=p1(i)
         p(i,2)=p2(i)
         p(i,3)=p3(i)
        enddo
         end
c-----------------------------------------------------------------
c----------------------------------------------------------------
         subroutine p2d_to_p1d_3(p,p1,p2,p3)
         implicit double precision (a-h,o-z)

         dimension p1(0:3),p2(0:3),p3(0:3),p(0:3,1:3)
        do i=0,3
         p1(i) = p(i,1)
         p2(i) = p(i,2)
         p3(i) = p(i,3)
        enddo
         end
c-----------------------------------------------------------------

c-----------------------------------------------------------------
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


      subroutine kinvar3_slicing(xx,xxjac,
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
c      yy1 = rapid(p3)

c     yy2: rapiditiy of photon 4
c      yy2 = rapid(p4)
      
c     rapidity YY12
      p2q = dot(p2,q)
      p1q = dot(p1,q)
      rr  = xa*p2q/(xb*p1q)
      YY12  = dlog(rr)/2.d0

c     cos-theta 
      ccst1=p3(3)/p3(0)
      ccst2=p4(3)/p4(0)            

c     ppt1 ppt2
c      ppt1 = eperp(p3)
c      ppt2 = eperp(p4)

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
c      rr34 = rjet(p3,p4)
c      rr35 = rjet(p3,p5)
c      rr45 = rjet(p4,p5)
c      
cc     transverse energy
c      EEt5 = eperp(p5)  ! for a massless particle Et=pt

      return
      end
 
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
