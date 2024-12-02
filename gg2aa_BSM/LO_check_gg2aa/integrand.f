cC ~~~~~~~~~~~~~~~~~~~C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
      integer function integrand(ndim,yy,ncomp,f,userdata,nvec ,core
     .                            ,weight,iter)
      implicit none 
c     cuba specific parameters
      integer ndim,ncomp,nvec,core,iter,userdata 
c      real*8 xx(ndim) ,f(ncomp),weight,fnlo3,yy(10)
      real*8 yy(2) ,f(ncomp),weight,flo2_LO,xx(10)
      external flo2_LO
       xx(1) = yy(1)
       xx(2) = yy(2)
      f(1)=flo2_LO(xx,weight)
c	print*,"Fun "
c	print*,xx
c      f= xx(1)*xx(2)*xx(3)*xx(4)*xx(5)*xx(6)
      return
      end
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function flo2_LO(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .          ,p(0:3,1:4),c(1:2),Born(1:2),coef(2,-2:0),SumI(-2:0)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=0.3894d9)
      character*50 name
      common/pdfname/name
      common/leg_choice/leg
      common/energy/s
      common/distribution/xq
      common/renor_scale/scale
      common/usedalpha/AL,ge
      external Born_uU2eE
       
c      yy(1) = 6.1844654984573387E-004 
c      yy(2) = 2.1478460799030053E-003

      tau = xq*xq/s

      xamin = tau
      xamax = 1.0d0
      xajac = (xamax - xamin)
      xa     = xamin+yy(1)*xajac
      xb = tau/xa
      xc = yy(2)

      rs  = dsqrt(s)
c      xa     = yy(1)
c      xb     = yy(2)

      rsp = dsqrt(xa*xb*s)
        
      ipass = 0
        eps = 0.5d0
       xlow = xq - eps
      xhigh = xq + eps

      xcut = xq - 10.0d0

c      if (rsp .gt. xcut) then


        yy(1) = xa
        yy(2) = xb
        yy(3) = xc
        
        call kinvar2(yy,xinvmass,p1,p2,p3,p4)
        call cuts0(p1,p2,p3,p4,ipass)
c	ipass =1

c        print*,ipass

        if (ipass .eq. 1) then 
        scale  = xinvmass

         if ( scale .ge. xlow .and. scale .le. xhigh) then 
             
              xmuf=scale
              xmur=scale
              xmu2=xmuf**2

              call pdf(xa,xmuf,f1)
              call pdf(xb,xmuf,f2)
              call setlum(f1,f2,xl)

               AL = alphasPDF(xmur)
              sig= xl(4)*Born_gg2aa(0,p1,p2,p3,p4)

              xnorm=hbarc2/16d0/pi/(xa*xb*s)
              wgt=xnorm*sig*vwgt
c              flo2_LO=wgt/vwgt/2d0/eps
              flo2_LO=wgt/vwgt*xajac*2*xq/s/xa
              
          else
              flo2_LO=0d0
          endif
       else                  
        flo2_LO=0d0
       endif
      return
      end

c---------------------------------------------------------------------
