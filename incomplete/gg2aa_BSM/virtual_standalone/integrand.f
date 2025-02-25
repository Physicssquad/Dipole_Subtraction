      function flo2_Vir(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .          ,p(0:3,1:4),c(1:2),coef(2,-2:0)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+6)
      parameter (cf=4d0/3d0)
      parameter (zeta2=1.64493406684823d0)
      parameter (Nc=3)
      parameter (Nf=5)
      parameter (Ca=3)
      parameter (Tf=0.5D0) 
      parameter (Tr=0.5D0) 
      parameter (EulerGamma = 0.5772156649015328606065120d0)

      character*50 name
      common/pdfname/name
      common/leg_choice/leg
      common/energy/s
      common/distribution/xq
      common/renor_scale/scale
      common/usedalpha/AL,ge
      external Born_gg2aa
       
      rs  = dsqrt(s)
      xa     = yy(1)
      xb     = yy(2)

      rsp = dsqrt(xa*xb*s)
        
      ipass = 0
      eps = 0.5d0
      xlow = xq - eps
      xhigh = xq + eps

      xcut = xq - 10.0d0

        call kinvar2(yy,xinvmass,p1,p2,p3,p4)
        scale  = xinvmass
        if ( scale .ge. xlow .and. scale .le. xhigh) then 
             
              xmuf=scale
              xmur=scale
              xmu2=xmuf**2

              call pdf(xa,xmuf,f1)
              call pdf(xb,xmuf,f2)
              call setlum(f1,f2,xl)
              AL = alphasPDF(xmur)

              e= DSQRT(ge*4.d0*PI)
              gs=DSQRT(Al*4.d0*PI)
              qu2=1d0                           ! charge is considered in lum.f

              s12 =  2d0*dot(p1,p2)
               t  = -2d0*dot(p1,p3)
               u  = -2d0*dot(p1,p4)

             Born = Born_gg2aa(0,p1,p2,p3,p4)

c**********[ From Slicing Code ] **********o
c        vir = Born*( -N*203.d0/9.d0
c     &            +1.*Tf*70.d0/9.d0
c     &            +N*zeta2*8.0d0 )
c	print*,"vir1:",vir

c        as = gs**2/16d0/Pi/Pi
        as = gs**2/4d0/Pi
        xai = 0d0

c        vir = (as*Born*(-(Ca*(203 + 18*EulerGamma**2 - 9*Pi**2-72*xai+
c     -           12*Log(2d0)*(11 + Log(64d0)) + 
c     -           6*Log(Pi)*(11 + Log(4096d0) + 3*Log(Pi)) - 
c     -           6*EulerGamma*(11 + Log(4096d0) + 6*Log(Pi)))) + 
c     -      2*nf*Tf*(35 - 12*EulerGamma + 12*Log(4*Pi)) + 
c     -      6*Log(s/xmu2)*(-4*nf*Tf + 
c     -         Ca*(11 - 6*EulerGamma + Log(4096d0) + 6*Log(Pi)) - 
c     -         3*Ca*Log(s12/xmu2))))/9.
c
c	print*,"vir2:",vir

c**********[ From Slicing Ends ] **********o
        model =1
        rp34 = dsqrt(s12)
        call coupfact(model,rp34,AK2D,AK2DINTF)

        fin = -0.0026041666666666665*
     -  (ak2D**2*as*(9*s12**4 + 9*t**4 + 28*t**3*u + 54*t**2*u**2 + 
     -       28*t*u**3 + 9*u**4 - 2*s12**2*(5*t**2 + 22*t*u + 5*u**2))*
     -     (CA - 2*nf*TF + 24*CA*Log(2d0)**2 - 6*CA*Log(4d0)**2))/
     -   (-1 + Nc**2)



c	print*,AK2D,EulerGamma,Pi,s,t,u
c	stop

c          eikonal_born = -0.00005425347222222222d0*
c     &  (AK2D**2*AL*(CA*(264d0*t**3*u + 462d0*t**2*u**2 + 264d0*t*u**3 +
c     &  t**4*(-235d0 + 132d0*EulerGamma + 12d0*Pi**2 - 264d0*Log(2d0)-
c     &       132d0*Log(Pi)) + 
c     &    u**4*(-235d0 + 132d0*EulerGamma + 12d0*Pi**2 -264d0*Log(2d0)-
c     &       132d0*Log(Pi))) - 
c     & 4d0*Nf*TR*(24d0*t**3*u + 42d0*t**2*u**2 + 24d0*t*u**3 + 
c     &    t**4*(-17d0 + 12d0*EulerGamma - 24d0*Log(2d0)-12d0*Log(Pi)) +
c     &    u**4*(-17d0 + 12d0*EulerGamma - 24d0*Log(2d0)-12d0*Log(Pi)))-
c     & 6d0*CF*(6d0*t**2*u**2*
c     &     (6d0 + 7d0*EulerGamma - 14d0*Log(2d0) - 7d0*Log(Pi)) + 
c     &    6d0*t**3*u*(3d0 + 4d0*EulerGamma - 8d0*Log(2d0)- 4d0*Log(Pi))+
c     &    6*t*u**3*(3d0 + 4d0*EulerGamma - 8d0*Log(2d0) - 4d0*Log(Pi))+
c     &    t**4*(3d0 + 6d0*EulerGamma**2 - 5d0*Pi**2 + 24d0*Log(2d0)**2-
c     &       15d0*Log(Pi) + 6d0*Log(Pi)**2 + 
c     &       6d0*Log(2d0)*(-5d0 + 4d0*Log(Pi)) - 
c     &       3d0*EulerGamma*(-5d0 + Log(256d0) + 4d0*Log(Pi))) + 
c     &    u**4*(3d0 + 6d0*EulerGamma**2 - 5d0*Pi**2+24d0*Log(2d0)**2 - 
c     &       15d0*Log(Pi) + 6d0*Log(Pi)**2 + 
c     &       6d0*Log(2d0)*(-5d0 + 4d0*Log(Pi)) - 
c     &       3d0*EulerGamma*(-5d0 + Log(256d0) + 4d0*Log(Pi)))) + 
c     & 6d0*(-22d0*CA*(t**4 + u**4) + 8d0*Nf*TR*(t**4 + u**4) + 
c     &    3d0*CF*(8d0*t**3*u + 14d0*t**2*u**2 + 8d0*t*u**3d0 + 
c     &       t**4*(5d0 + 4d0*EulerGamma - 8d0*Log(2d0)- 4d0*Log(Pi)) + 
c     &       u**4*(5d0 + 4d0*EulerGamma - 8d0*Log(2d0)- 4d0*Log(Pi))))*
c     &  Log(xmu2/s12) - 36d0*CF*(t**4 + u**4)*Log(xmu2/s12)**2))/Pi 
c
c         eikonal_born=-0.00005425347222222222*
c     -  (AK2D**2*Al*(CA*(264*t**3*u + 462*t**2*u**2 + 264*t*u**3 + 
c     -         t**4*(-235 + 132*EulerGamma + 12*Pi**2 - 264*dlog(2d0) - 
c     -             132*dlog(Pi)) + 
c     -         u**4*(-235 + 132*EulerGamma + 12*Pi**2 - 264*dlog(2d0) - 
c     -             132*dlog(Pi))) - 
c     -       4*Nf*TR*(24*t**3*u + 42*t**2*u**2 + 24*t*u**3 + 
c     -        t**4*(-17 + 12*EulerGamma - 24*dlog(2d0) - 12*dlog(Pi)) + 
c     -       u**4*(-17 + 12*EulerGamma - 24*dlog(2d0) - 12*dlog(Pi))) - 
c     -       6*CF*(6*t**2*u**2*
c     -           (6 + 7*EulerGamma - 14*dlog(2d0) - 7*dlog(Pi)) + 
c     -         6*t**3*u*(3 + 4*EulerGamma - 8*dlog(2d0) - 4*dlog(Pi)) + 
c     -         6*t*u**3*(3 + 4*EulerGamma - 8*dlog(2d0) - 4*dlog(Pi)) + 
c     -          t**4*(3 + 6*EulerGamma**2 - 5*Pi**2 + 24*dlog(2d0)**2 - 
c     -             15*dlog(Pi) + 6*dlog(Pi)**2 + 
c     -             6*dlog(2d0)*(-5 + 4*dlog(Pi)) - 
c     -             3*EulerGamma*(-5 + dlog(256d0) + 4*dlog(Pi))) + 
c     -          u**4*(3 + 6*EulerGamma**2 - 5*Pi**2 + 24*dlog(2d0)**2 - 
c     -             15*dlog(Pi) + 6*dlog(Pi)**2 + 
c     -             6*dlog(2d0)*(-5 + 4*dlog(Pi)) - 
c     -             3*EulerGamma*(-5 + dlog(256d0) + 4*dlog(Pi)))) + 
c     -       6*(-22*CA*(t**4 + u**4) + 8*Nf*TR*(t**4 + u**4) + 
c     -          3*CF*(8*t**3*u + 14*t**2*u**2 + 8*t*u**3 + 
c     -             t**4*(5 + 4*EulerGamma - 8*dlog(2d0) - 4*dlog(Pi)) + 
c     -             u**4*(5 + 4*EulerGamma - 8*dlog(2d0) - 4*dlog(Pi))))*
c     -      dlog(xmu2/s12) - 36*CF*(t**4 + u**4)*dlog(xmu2/s12)**2))/Pi
c        fin = -0.0026041666666666665*
c     -  (ak2D**2*as*(9*s**4 + 9*t**4 + 28*t**3*u + 54*t**2*u**2 + 
c     -       28*t*u**3 + 9*u**4 - 2*s**2*(5*t**2 + 22*t*u + 5*u**2))*
c     -     (CA - 2*nf*Tf + 24*CA*Log(2d0)**2 - 6*CA*Log(4d0)**2))/
c     -   (-1 + Nc**2)


c        sig= xl(4)*(vir + eikonal_born)  
        sig= xl(4)*(fin)
c	print*,eikonal_born
c	print*,vir
c	print*, " "  
c	stop
c        sig= xl(4)*Vir 
c       sig= xl(1)*SumI(1)   

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
 
