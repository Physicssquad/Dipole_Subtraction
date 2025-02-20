C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
      integer function integrand(ndim,xx,ncomp,f,userdata,nvec ,core
     .                            ,weight,iter)
      implicit none 
c     cuba specific parameters
      integer n4,ndim,ncomp,nvec,core,iter,userdata 
      real*8 xx(ndim) ,f(ncomp),weight,flo2_Vir
      common/countc/n4
      external flo2_Vir
      f(1) = flo2_Vir(xx,weight)
      return
      end
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function flo2_Vir(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .          ,p(0:3,1:4),c(1:2),Born(1:2),coef(2,-2:0),SumI(2)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+6)
      character*50 name
      common/pdfname/name
      common/leg_choice/leg
      common/energy/s
      common/distribution/xq
      common/renor_scale/scale
      common/usedalpha/AL,ge
      external Born_uU2eE
       
      rs  = dsqrt(s)
      xa     = yy(1)
      xb     = yy(2)

      rsp = dsqrt(xa*xb*s)
        
      ipass = 0
      eps = 0.5d0
      xlow = xq - eps
      xhigh = xq + eps

      xcut = xq - 10.0d0

c      if (rsp .gt. xcut) then

        call kinvar2(yy,xinvmass,p1,p2,p3,p4)
        scale  = xinvmass
        if ( scale .ge. xlow .and. scale .le. xhigh) then 
             
              xmuf=scale
              xmur=scale
              xmu2=xmuf**2
c              xmuf=xq
c              xmur=xq

              call pdf(xa,xmuf,f1)
              call pdf(xb,xmuf,f2)
              call setlum(f1,f2,xl)
              AL = alphasPDF(xmur)
c               AL = 2d0*PI
c               AL = 1d0

              e= DSQRT(ge*4.d0*PI)
              gs=DSQRT(Al*4.d0*PI)
c              qu2=4d0/9d0
              qu2=1d0                           ! charge is considered in lum.f
              Cf = 4d0/3d0

              s12 =  2d0*dot(p1,p2)
               t  = -2d0*dot(p1,p3)
               u  = -2d0*dot(p1,p4)

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c            Only finite part of loop is considered here                    
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c              VIR= gs**2*8*e**4*qu2*(5*(-6 + 11*Pi**2)*s12**2 +
c     -        2*(-48 + 55*Pi**2)*s12*t + 2*(-48 + 55*Pi**2)*t**2 -
c     -        6*(s12**2 + 6*s12*t + 6*t**2)*Log(xmu2/s12) -
c     -        6*(s12**2 + 2*s12*t + 2*t**2)*(Log(xmu2/s12))**2)/
c     -          (24*Pi**2*s12**2)
c              VIR = VIR/36d0    ! spin and color avg


         VIR=  (gs**2*8d0*e**4*qu2*
     -    (5*(-6 + 11*Pi**2)*s12**2 + 2*(-48 + 55*Pi**2)*s12*t +
     -      2*(-48 + 55*Pi**2)*t**2 -
     -      6*(s12**2 + 6*s12*t + 6*t**2)*Log(xmu2/s12) -
     - 6*(s12**2+2*s12*t+2*t**2)*Log(xmu2/s12)**2))/(24.*Pi**2*s12**2)
c              VIR = VIR/36d0    ! spin and color avg


c              (gs**2*NA*e**2*qu2*
c     -    (5*(-6 + 11*Pi**2)*s**2 + 2*(-48 + 55*Pi**2)*s*t +
c     -      2*(-48 + 55*Pi**2)*t**2 -
c     -      6*(s**2 + 6*s*t + 6*t**2)*Log(mu2/s) -
c     -      6*(s**2 + 2*s*t + 2*t**2)*Log(mu2/s)**2))/(24.*Pi**2*s**2)

c
c        VIR= (AL*CF*3*e**4*qu2*
c     -    ((-30 + 7*Pi**2)*s12**2 + 2*(-48 + 7*Pi**2)*s12*t + 
c     -      2*(-48 + 7*Pi**2)*t**2 - 
c     -   6*(s12**2 + 6*s12*t + 6*t**2)*Log(xmu2/s12) - 
c     -   6*(s12**2 + 2*s12*t + 2*t**2)*Log(xmu2/s12)**2))/(3.*Pi*s12**2)
c              VIR = 2.0d0*VIR/36d0    ! spin and color avg

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               

               gammaq = 3d0/2d0*Cf
               xkq     = Cf*(7d0/2d0-pi*pi/6d0)
           EulerGamma = 0.5772156649015328606065120d0
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c            Eikonal part = SumI * Born * (colour factors)           
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c           SumI(1)=   (-8*Al*CF*e**4*qu2*


         eikonal = (8*AL*e**4*qu2*(12d0*xkq*(s12**2+2*s12*t + 2*t**2)+
     -      CF*(6d0*EulerGamma**2*(s12**2 + 2d0*s12*t + 2d0*t**2) - 
     -         5d0*Pi**2*(s12**2 + 2d0*s12*t + 2d0*t**2) - 
     -         6d0*EulerGamma*(s12**2*(1d0 + Log(16d0) + 2d0*Log(Pi)) + 
     -            2d0*s12*t*(3d0 + Log(16d0) + 2d0*Log(Pi)) + 
     -            2d0*t**2*(3d0 + Log(16d0) + 2d0*Log(Pi))) + 
     -    6d0*(2d0*s12*t*(3d0 + 4d0*Log(2d0)**2+Log(64d0)+3d0*Log(Pi) + 
     -               Log(16d0)*Log(Pi) + Log(Pi)**2) + 
     -       2d0*t**2*(3d0 + 4d0*Log(2d0)**2 + Log(64d0)+3d0*Log(Pi) + 
     -         Log(16d0)*Log(Pi) + Log(Pi)**2) + 
     -       s12**2*(4d0*Log(2d0)**2 + Log(16d0)*Log(Pi) + Log(Pi)**2 + 
     -               Log(4d0*Pi)))) - 
     - 6d0*CF*(2d0*s12*t*(-3d0+2d0*EulerGamma-4d0*Log(2d0)-2d0*Log(Pi))+
     - 2d0*t**2*(-3d0 + 2d0*EulerGamma - 4d0*Log(2d0) - 2d0*Log(Pi)) + 
     -     s12**2*(-1d0+2d0*EulerGamma - 4d0*Log(2d0) - 2d0*Log(Pi)))*
     -   Log(xmu2/s12) + 6d0*CF*(s12**2 + 2d0*s12*t + 
     -  2d0*t**2)*Log(xmu2/s12)**2
     -      ))/(3d0*Pi*s12**2)
c          eikonal=eikonal/36d0
c          SumI(1) = VIR + eikonal
c          print*,SumI(1),"separated"
c        SumI(1) = eikonal
        SumI(1) = VIR


c     -   (12*xkq*s12**2 -5*Pi**2*s12**2 + 24*gammaq*s12*t+24*xkq*s12*t- 
c     -      10*Pi**2*s12*t + 24*gammaq*t**2 + 24*xkq*t**2 - 
c     -      10*Pi**2*t**2 + 6*EulerGamma**2*(s12**2 + 2*s12*t+2*t**2) - 
c     -      12*s12**2*Log(4*Pi) + 12*gammaq*s12**2*Log(4*Pi) + 
c     -      24*gammaq*s12*t*Log(4*Pi) + 24*gammaq*t**2*Log(4*Pi) + 
c     -      6*s12**2*Log(4*Pi)**2 + 12*s12*t*Log(4*Pi)**2 + 
c     -      12*t**2*Log(4*Pi)**2 - 
c     -      12*EulerGamma*(s12**2*(-1 + gammaq + Log(4*Pi)) + 
c     -      2*s12*t*(gammaq + Log(4*Pi)) + 2*t**2*(gammaq + Log(4*Pi)))
c     -        - 12*(2*s12*t*(EulerGamma - gammaq - Log(4*Pi)) + 
c     -         2*t**2*(EulerGamma - gammaq - Log(4*Pi)) + 
c     -       s12**2*(1 + EulerGamma - gammaq - Log(4*Pi)))*Log(xmu2/s12)
c     -  + 6*(s12**2 + 2*s12*t +2*t**2)*Log(xmu2/s12)**2))/(3.*Pi*s12**2)
c
           SumI(2) = 
     .    (-4*Al*CF*3*e**4*qu2*(s12**2 + 2*s12*t + 2*t**2)*(-2d0 +
     .        0*Pi**2))/(Pi*s12**2)
             SumI(2) = SumI(2)/36d0 ! colour average performed here
c          print*,SumI(2),"original "
c          stop


c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


c              sig= xl(1)*(Vir + 1d0*SumI(1))    ! factor of 2 for two leg contribution in eikonal
              sig= xl(1)*SumI(2)    ! factor of 2 for two leg contribution in eikonal
c              sig= xl(1)*Vir

              xnorm=hbarc2/16d0/pi/(xa*xb*s)
              wgt=xnorm*sig*vwgt
              flo2_Vir=wgt/vwgt/2d0/eps

            return
       else                  
        flo2_Vir=0d0
        return
       endif
c       else
c        flo2_Vir=0d0
c       endif
      end

c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     [u U -> e E]  Born 
c--------------------------------------------------------------------o
       function Born_uU2eE(k,p1,p2,p3,p4)
       implicit double precision (a-h,o-z)
       dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
       parameter(PI=3.141592653589793238D0)
       common/usedalpha/AL,ge
       e= DSQRT(ge*4.d0*PI)
       gs=DSQRT(Al*4.d0*PI)
      IF(k .eq. 0)  CF =  1d0                   !Leading Order K=0 
      IF(k .eq. 1)  CF = -4d0/3d0               !leg 1
      IF(k .eq. 2)  CF = -4d0/3d0               !Leg 2    
      s13 =  2.0d0*dot(p1,p3) ! t
      s23 =  2.0d0*dot(p2,p3) ! u
      s12 =  2.0d0*dot(p1,p2) ! s
      XNC = 1/4d0
      xnorm =1d0
      qu2 = 1d0!4d0/9d0

      Born_uU2eE=CF*(2*e**4*qu2*(-2*s13*s23 + s12*(s13 +
     .            s23)))/(3d0*s12**2)
       return
       end
c---------------------------------------------------------------------

