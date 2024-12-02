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

