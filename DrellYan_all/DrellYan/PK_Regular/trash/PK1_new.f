! a -> ai+ i ||  q -> q g || P{a,ai} -> P{q,q} || a- Incoming
!                                               | ai- Undergoing Born
!                                               |    process               
! Kb = Kbar      
c----------------------------------------------------------- 
      subroutine getPK(k,x,xmuf,p,xp1,xp2,ALLP,ALLK)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      dimension AllP(5),AllK(5)
      dimension  p1(0:3),p2(0:3),p3(0:3),p4(0:3),p(0:3,1:4)
      dimension  xp1(0:3),xp2(0:3)
      common /usedalpha/ AL      
      external Born_uU2eE,pqq
      call p2d_to_p1d_4(p,p1,p2,p3,p4)
        s12 = 2d0*dot(p1,p2)
        Cf = 4d0/3d0                      
        Tr = 0.5d0
        Alp = Al/2.0d0/pi
        xmuf2 = xmuf*xmuf

      if (k .eq. 1) then   ! Cloice for [Leg1]       

c       [P term]
        coefx = Alp*(-1.0d0)*dlog(xmuf2/s12/x)    !{a,ai,i} => {q,q,g}
        Bornx = Born_uU2eE(0,xp1,p2,p3,p4)
        coefx = coefx*Bornx

        coef1 = Alp*(-1.0d0)*dlog(xmuf2/s12)    !{a,ai,i} => {q,q,g}
        Born1 = Born_uU2eE(0,p1,p2,p3,p4)
        coef1 = coef1*Born1 

        ALLP(1) = pqq(x)*(coefx - coef1)

c        AllP(2) = Al*Bornx*(Tr*(2*x**2+1-
c     .                 2*x)*(-1)*Log(xmuf**2/(s12*x)))/2d0/Pi    !{a,ai,i} => {g,q,q~}

c       [K term]

       BCf*((1d0+x)*Log((1d0-x)/x)+(1d0-x)-(1d0+x)*Log(1d0-x))

       Aplus=Plusd(2d0/(1d0-x)*Log((1d0-x)/x))+
     .    Plusd(2d0/(1d0-x)*dLog(1d0-x))
        call p2d_to_p1d_4(p,p1,p2,p3,p4)
        Bornx = Born_uU2eE(0,p1,p2,p3,p4)
        coefx = Al*Cf*Bornx/2d0/Pi
        B = B * Bornx
        p1(0) = p1(0)/x
        p1(3) = p1(3)/x
        Born1 = Born_uU2eE(0,p1,p2,p3,p4)
        coef1 = Al*Cf*Born1/2d0/Pi
        AllK(1)= Aplus*(coefx - coef1) - B 

        AllK(2) = Al*Bornx*(Tr*(x**2+(1-x)**2)*dLog((1d0-x)/x) +
     .           Tr*2d0*x*(1d0-x)-CLF*Tr*(x**2+
     .              (1-x)**2)*dLog(1d0-x))/2d0/Pi  


      elseif (k .eq. 2 ) then  ! Cloice for [Leg2]

c       [P term]
        coefx = Al*(Cf*(-1)*Log(xmuf**2/(s12*x)))/2d0/Pi    !{a,ai,i} => {q,q,g}
        Bornx = Born_uU2eE(0,p1,p2,p3,p4)
        coefx = coefx*Bornx

        coef1 = Al*(Cf*(-1)*Log(xmuf**2/(s12)))/2d0/Pi
        p2(0) = p2(0)/x
        p2(3) = p2(3)/x
        Born1 = Born_uU2eE(0,p1,p2,p3,p4)
        coef1 = coef1*Born_uU2eE(0,p1,p2,p3,p4) 

        ALLP(1) = Plusd((1+x**2)/(1-x))*(coefx - coef1)

        AllP(2) = Al*Bornx*(Tr*(2*x**2+1-
     .                 2*x)*(-1)*Log(xmuf**2/(s12*x)))/2d0/Pi    !{a,ai,i} => {g,q,q~}

c       [K term]
      B=Al*Cf*((1d0+x)*Log((1d0-x)/x)+(1d0-x)-(1d0+x)*Log(1d0-x))/2d0/Pi

       Aplus=Plusd(2d0/(1d0-x)*Log((1d0-x)/x))+
     .    Plusd(2d0/(1d0-x)*dLog(1d0-x))
        call p2d_to_p1d_4(p,p1,p2,p3,p4)
        Bornx = Born_uU2eE(0,p1,p2,p3,p4)
        coefx = Al*Cf*Bornx/2d0/Pi
        B = B * Bornx
        p2(0) = p2(0)/x
        p2(3) = p2(3)/x
        Born1 = Born_uU2eE(0,p1,p2,p3,p4)
        coef1 = Al*Cf*Born1/2d0/Pi
        AllK(1)= Aplus*(coefx - coef1) - B 

        AllK(2) = Al*Bornx*(Tr*(x**2+(1-x)**2)*dLog((1d0-x)/x) +
     .           Tr*2d0*x*(1d0-x)-CLF*Tr*(x**2+
     .              (1-x)**2)*dLog(1d0-x))/2d0/Pi  


      endif

      return
      end
c-----------------------------------------------------------  !P-terms
c----------------------------------------------------------- !Misc-functions
c                                                            
      double precision function Pqq(x)
      implicit double precision (a-h,o-z)
      CF = 4.0d0/3.0d0
      pqq = CF*(1.0d0+x*x)/(1.0d0-x)
      return
      end

      double precision function Pqqreg(x)
      implicit double precision (a-h,o-z)
      CF = 4.0d0/3.0d0
      pqq = -CF*(1.0d0+x)
      return
      end

      double precision function KbarP_qq(x)
      implicit double precision (a-h,o-z)
      
      dlgx = dlog((1.0d0-x)/x)
      CF = 4.0d0/3.0d0
      KbarP_qq = CF*2.0d0/(1.0d0-x)*dlgx 
      return
      end

      double precision function Kbar_qq(x)
      implicit double precision (a-h,o-z)
      
      dlgx = dlog((1.0d0-x)/x)
      CF = 4.0d0/3.0d0
      Kbar_qq = -CF*(1.0d0+x)*dlgx + (1-x) 
      return
      end

      double precision function KbarD_qq(x)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      
      CF = 4.0d0/3.0d0
      KbarD_qq = -CF*(5.0d0-pi*pi)
      return
      end

      double precision function KtilP_qq(x)
      implicit double precision (a-h,o-z)
      
      CF = 4.0d0/3.0d0
      dlgx = dlog(1.0d0-x)
      KtilP_qq = CF*2.0d0/(1.0d0-x)*dlgx 
      return
      end


      double precision function Ktilreg_qq(x)
      implicit double precision (a-h,o-z)
      external Pqqreg
      
      dlgx = dlog(1.0d0-x)
      Ktilreg_qq = Pqqreg(x)*dlgx
      return
      end

      double precision function KtilD_qq(x)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      
      KtilD_qq = -pi*pi/3.0d0
      return
      end
