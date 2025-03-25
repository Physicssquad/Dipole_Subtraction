cC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
c      integer function integrand(ndim,xx,ncomp,f,userdata,nvec ,core
c     .                            ,weight,iter)
c      implicit none 
cc     cuba specific parameters
c      integer n4,ndim,ncomp,nvec,core,iter,userdata 
c      real*8 xx(ndim) ,f(ncomp),weight,flo2_Vir
c      common/countc/n4
c      external flo2_LO
c      f(1) = flo2_LO(xx,weight)
c      return
c      end
cC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function flo2_LO(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .          ,p(0:3,1:4),c(1:2),Born(1:2),coef(2,-2:0),SumI(-2:0)
      parameter (pi=3.14159265358979d0)
c      parameter (hbarc2=0.3894d9)  ! in pb
      parameter (hbarc2=0.3894d12)  ! in fb
      character*50 name
      common/pdfname/name
      common/leg_choice/leg
      common/energy/s
      common/distribution/xq
      common/renor_scale/scale
      common/usedalpha/AL,ge
      common/amass/am1,am2,am3,am4,am5
      external amp_LO
       
      rs  = dsqrt(s)
      xa     = yy(1)
      xb     = yy(2)

      rsp = dsqrt(xa*xb*s)
        
      ipass = 0
        eps = 1.0d0
c         xq = 250d0
       xlow = xq - eps
      xhigh = xq + eps

        call kinvar2_type_1(yy,scale,p1,p2,p3,p4)
        if ( scale .ge. xlow .and. scale .le. xhigh) then 
             
              xmuf=xq
              xmur=xq
              xmu2=xmuf**2
c              xmuf=xq
c              xmur=xq

              call pdf(xa,xmuf,f1)
              call pdf(xb,xmuf,f2)
              call setlum(f1,f2,xl)
              AL = alphasPDF(xmur)

              sig= xl(2)*amp_LO(AL,xmur,p1,p2,p3,p4)

              xnorm=hbarc2/16d0/pi/(xa*xb*s)
              wgt=xnorm*sig*vwgt
              flo2_LO=wgt/vwgt/2d0/eps
       else
        flo2_LO=0d0
       endif
      return
      end

c---------------------------------------------------------------------
