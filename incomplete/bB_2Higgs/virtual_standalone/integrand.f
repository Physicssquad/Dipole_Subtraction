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
      parameter (N=3)
      parameter (Ca=3)
      parameter (Tf=0.5D0) 
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

        vir = Born*( -N*203.d0/9.d0
     &            +1.*Tf*70.d0/9.d0
     &            +N*zeta2*8.0d0 )

        eikonal= (-0.013888888888888888d0*
     &  (AL*Born*(8*N*Tf*(-8 + 3*EulerGamma - 6*Log(2.) - 3*Log(Pi)) + 
     &       CA*(200 + 18*EulerGamma**2 - 21*Pi**2 + 132*Log(2.) + 
     &          72*Log(2.)**2 + 66*Log(Pi) + 72*Log(2.)*Log(Pi) + 
     &          18*Log(Pi)**2 - 
     &          6*EulerGamma*(11 + Log(4096.) + 6*Log(Pi))) - 
     &       6*(4*N*Tf + CA*
     &           (-11 + 6*EulerGamma - 12*Log(2.) - 6*Log(Pi)))*
     &        Log(xmu2/s12) + 18*CA*Log(xmu2/s12)**2))/Pi)

c	print*,"Virtual:",vir
c	print*,"I-term :",eikonal
c	print*," "

        sig= xl(4)*(Vir + eikonal)  
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
 
