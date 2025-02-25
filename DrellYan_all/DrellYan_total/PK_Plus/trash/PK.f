! a -> ai+ i ||  q -> q g || P{a,ai} -> P{q,q} || a- Incoming
!                                               | ai- Undergoing Born
!                                               |    process               
! Kb = Kbar      
c----------------------------------------------------------- 
      subroutine getPK(k,x,xmuf,p,SumP,SumK)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      dimension AllP(5),AllK(5)
      dimension  p1(0:3),p2(0:3),p3(0:3),p4(0:3),p(0:3,1:4)
      common /usedalpha/ AL      
      call p2d_to_p1d_4(p,p1,p2,p3,p4)
        s12 = 2d0*dot(p1,p2)
        CLF1 = 1d0
        CLF2 = 1d0                                  !                 {Ta . Tai}
        Cf = 4d0/3d0                                 !  Colour factor  ----------
        Tr = 0.5d0
c        Al= 0.118d0                                 !                   Tai^2
      if (k .eq. 1) then   ! Cloice for [Leg1]       
!       
        AllP(1) = Al*(Cf*Plusd((1+x**2)/(1-x))*
     .            CLF1*Log(xmuf**2/(s12*x)))/2d0/Pi    !{a,ai,i} => {q,q,g}
        AllP(2) = Al*(Tr*(2*x**2+1-
     .                 2*x)*CLF1*Log(xmuf**2/(s12*x)))/2d0/Pi    !{a,ai,i} => {g,q,q~}
c       AllP(3) = Al*(Pgq(x)*Clf*dLog(xmuf**2/(s12*x)))/2d0/Pi    !{a,ai,i} => {g,q~,q}
        AllK(1) = Al*(Cf*(Plusd(2d0/(1d0-x)*dLog((1d0-x)/x))-
     .        (1d0+x)*dLog((1d0-x)/x)+(1d0-x)-
     .         delta(1d0-x)*(5d0-Pi**2))- Cf*(1d0-x)*dLog(1d0-x) + 
     .    Cf*(Plusd(2d0/(1d0-x)*dLog(1d0-x))-
     .              Pi**2*delta(1d0-x)/3d0))/2d0/Pi !! #colour factor unclear

        AllK(2) = Al*(Tr*(x**2+(1-x)**2)*dLog((1d0-x)/x) +
     .           Tr*2d0*x*(1d0-x)-CLF*Tr*(x**2+
     .              (1-x)**2)*dLog(1d0-x))/2d0/Pi  
c        AllK(3) = Al*(Kbgq(x)-CLF*Kgq_til(x))/2d0/Pi           !! #colour factor unclear
c        print*,AllK(1),AllK(2)
c        print*,xKbqq(x),xKbgq(x)
c        print*,xKqq_til(x),xKgq_til(x)
        
        SumP = 0d0 
        SumK = 0d0

       do i=1,2
        SumP = SumP + AllP(i)
        SumK = SumK + AllK(i)
       enddo
c       print*,SumP,SumK," <--"
        return

      elseif (k .eq. 2) then  ! Cloice for [Leg2]

        AllP(1) = Al*(Pqq(x)*CLF*dLog(xmuf**2/(s12*x)))/2d0/Pi
        AllP(2) = Al*(Pgq(x)*Clf*dLog(xmuf**2/(s12*x)))/2d0/Pi
c        AllP(3) = Al*(Pgq(x)*Clf*dLog(xmuf**2/(s12*x)))/2d0/Pi !Contribution for anti-quark

        AllK(1) = Al*(Cf*(Plusd(2d0/(1d0-x)*dLog((1d0-x)/x))-
     .        (1d0+x)*dLog((1d0-x)/x)+(1d0-x)-
     .         delta(1d0-x)*(5d0-Pi**2))- Cf*(1d0-x)*dLog(1d0-x) + 
     .    Cf*(Plusd(2d0/(1d0-x)*dLog(1d0-x))-
     .              Pi**2*delta(1d0-x)/3d0))/2d0/Pi !! #colour factor unclear

        AllK(2) = Al*(Tr*(x**2+(1-x)**2)*dLog((1d0-x)/x) +
     .           Tr*2d0*x*(1d0-x)-CLF*Tr*(x**2 +
     .              (1-x)**2)*dLog(1d0-x))/2d0/Pi   !! #colour factor unclear
c        AllK(3) = Al*(Kbgq(x)-CLF*Kgq_til(x))/2d0/Pi !! #colour factor unclear

      endif
c        SumP = 0d0 
c        SumK = 0d0
c       do i=1,2
c        SumP = SumP + AllP(i)
c        SumK = SumK + AllK(i)
c       enddo
c       print*,"P:",SumP
c       print*,"K:",SumK

      return
      end
c-----------------------------------------------------------  !P-terms
      function Pqg(x) ! OR Pqbg antiquark also same                                    
      implicit double precision (a-h,o-z)
        Cf = 4d0/3d0
         Pqg = Cf*(1d0+(1d0-x)**2)/x
      return
      end

      function Pgq(x)
      implicit double precision (a-h,o-z)
         Tr = 0.5d0
         Pgq = Tr*(x**2+(1-x)**2)
      return
      end

      function Pqq(x)
      implicit double precision (a-h,o-z)
        Cf = 4d0/3d0
        Pqq = Cf*Plusd((1+x**2)/(1-x))
      return
      end

      function Pgg(x)
      implicit double precision (a-h,o-z)
         Ca = 3d0
         Nf = 1d0
         Tr = 0.5d0
         Pgg = 2d0*Ca*(Plusd(1d0/(1d0-x)) + (1d0-x)/x -1d0+ x*(1d0-x)) +
     .         delta(1d0-x)*(11d0*Ca/6d0 - 2d0*Nf*Tr/3d0)
      return
      end
c-----------------------------------------------------------!K-terms

      function xKbqq(x)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
        Cf = 4d0/3d0
        xKbqq = Cf*(Plusd(2d0/(1d0-x)*dLog((1d0-x)/x))-
     .        (1d0+x)*dLog((1d0-x)/x)+(1d0-x)-
     .         delta(1d0-x)*(5d0-Pi**2))
      return
      end                

      function xKbgq(x) !Kb[GQ]
      implicit double precision (a-h,o-z)
        Tr = 0.5d0 
        xKbgq = Pgq(x)*dLog((1d0-x)/x) + Tr*2d0*x*(1d0-x)
      return
      end                

      function xKbqg(x) !Kb[QG]
      implicit double precision (a-h,o-z)
        Tr = 0.5d0 
        Cf = 4d0/3d0
        xKbqg = Pqg(x)*dLog((1d0-x)/x) + Cf*x 
      return
      end                

      function xKqq_til(x)
      implicit double precision (a-h,o-z)
        Cf = 4d0/3d0
        xKqq_til = Cf*(1d0-x)*dLog(1d0-x) + 
     .    Cf*(Plusd(2d0/(1d0-x)*dLog(1d0-x))-Pi**2*delta(1d0-x)/3d0)
      return
      end

      function xKgq_til(x)
      implicit double precision (a-h,o-z)
        xKgq_til = Pgq(x)*dLog(1d0-x)  
      return
      end
c----------------------------------------------------------- !Misc-functions
c                                                            
c      function Plusd(xx)
c      implicit double precision (a-h,o-z)
c        Plusd = xx
c      return
c      end
c
c      function delta(xx)
c      implicit double precision (a-h,o-z)
c        delta = 0d0
c      return
c      end
