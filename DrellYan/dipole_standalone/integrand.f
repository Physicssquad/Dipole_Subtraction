C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
      integer function integrand(ndim,xx,ncomp,f,userdata,nvec ,core
     .                            ,weight,iter)
      implicit none 
c     cuba specific parameters
      integer n4,ndim,ncomp,nvec,core,iter,userdata 
      real*8 xx(ndim) ,f(ncomp),weight,fnlo3
      common/countc/n4
      external fnlo3
      f(1) = fnlo3(xx,weight)
      return
      end
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
c      Our vegas specific
      double precision function fnlo3(xx,weight)
      implicit double precision(a-h,o-z) 
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
      dimension  xx(10),p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),q(0:3)
     .             ,p(0:3,1:5),dip(2),B1(1:2),xl(15),f1(-6:6),f2(-6:6)
     .             , sig(1:5),SumD(1:2)
      integer i35,i45,is5,itest
      integer k1,k2,k3,ipass,n4,unphy
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+6)
      common/amass/am1,am2,am3,am4,am5
      common/energy/s
      common/set/set1
      common/distribution/xq
      common/countc/n4
      common/usedalpha/AL,ge
      external dipole_gq_q

      xa = xx(1)
      xb = xx(2)
      rsp = dsqrt(xa*xb*s)

      xmz = 91.1876d0

      ipass1 = 0

      eps = 0.5d0
      xlow = xq - eps
      xhigh = xq + eps

c      xcut = xq - 10.0d0

c      if (rsp .gt. xcut) ipass1 = 1

c      if (ipass1 .eq.1) then
       call kinvar3(xx,xxjac,xinvmass,p1,p2,p3,p4,p5,unphy)
       if (unphy .eq. 0) then ! with zero unphysical PS points proceed
     
        scale = xinvmass

        ipass = 0
        fnlo3 = 0

        if ( scale .ge. xlow .and. scale .le. xhigh) then

          xmuf=scale
          xmur=scale

c           xmuf=xq
c           xmur=xq

          call pdf(xa,xmuf,f1)
          call pdf(xb,xmuf,f2)
          call setlum(f1,f2,xl)
           AL = alphasPDF(xmur)
c           AL = 2d0*PI 


          call p1dtop2d_5(p1,p2,p3,p4,p5,p)

          call  uu2ee_r(p,sig)
          SumD(1) = dipole_uU_g(1,p) + dipole_uU_g(2,p)
          SumD(2) = dipole_gq_q(1,p) + dipole_gq_q(2,p)

c           sig(2) = xl(8)*( sig(2) - dipole_gq_q(2,p )) 
c           sig(3) = xl(7)*( sig(3) - dipole_gq_q(1,p )) 


c         sig(1) = xl(1)*(sig(1) - SumD(1))

c           sig(1) = SumD(1)

c          sig(2) = xl(7)*(sig(2) - SumD(2))    ! for gq channel
c          sig(3) = xl(7)*(sig(3) - SumD(3))    ! for gq channel


c          sigma = sig(2) + sig(3)
          sigma = xl(1)*(sig(1)-sumD(1))
c          if (sigma .ne. sigma ) sigma =0d0 

          pi_1 = 0.5d0*rsp
          flux = 4d0*pi_1*rsp
          xnorm=hbarc2/8d0/(2d0*Pi)**4/flux
          wgt=xxjac*xnorm*sigma*weight
          fnlo3=wgt/weight/2d0/eps
         endif
       endif
151   return
      end
c---------------------------------------------------------------------
