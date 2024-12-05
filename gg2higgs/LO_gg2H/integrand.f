c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function flo1_LO(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(2)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .          ,p(0:3,1:4),c(1:2),Born(1:2),coef(2,-2:0),SumI(-2:0)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=0.3894d9)      ! in Pb
c      parameter (hbarc2=0.3894d12)    ! in Fb
      common/energy/S
      common/renor_scale/scale
      common/usedalpha/AL,ge
      common/amass/am1,am2,amH,am4,am5
      external Born_gg2H
       
        tau = amH**2/S
c         AL = alphasPDF(amH)

         xajac = 1d0 - tau
         xa = tau + xajac*yy(1)
         xb = tau/xa 

         sp = xa*xb*S
        rsp = dsqrt(sp)

        call kinvar1(xa,xb,p1,p2,p3)

        scale = dsqrt(2d0*dot(p1,p2))
c	write(*,*)'sp,scale =',xa,xb,sp,scale

c        xmuf= scale
        xmuf= 125d0
        xmur= xmuf
        xmu2=xmuf**2

         AL = alphasPDF(xmur)

        call pdf(xa,xmuf,f1)
        call pdf(xb,xmuf,f2)
        call setlum(f1,f2,xl)

        sig= xl(2)*Born_gg2H(0,p1,p2,p3)


c        pi_1 = 0.5d0*rsp
c        flux = 4d0*pi_1*rsp
         pi_1 = PI/amH
         flux = 2d0*sp 

c        xnorm=hbarc2/8d0/(2d0*Pi)**4/flux 
        xnorm=hbarc2*pi_1/flux

c        flo1_LO  = xajac * xnorm * sig* xb
        flo1_LO  = xajac * xnorm * sig * 2d0* amh/xa/S
c        flo1_LO  = xnorm * sig/2d0/eps

      return
      end
c---------------------------------------------------------------------
c
cc~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c      function flo1_LO(yy,vwgt)
c      implicit double precision (a-h,o-z)
c      dimension yy(2)
c      dimension f1(-6:6),f2(-6:6),xl(15)
c      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
c     .          ,p(0:3,1:4),c(1:2),Born(1:2),coef(2,-2:0),SumI(-2:0)
c      parameter (pi=3.14159265358979d0)
cc      parameter (hbarc2=0.3894d9)      ! in Pb
c      parameter (hbarc2=0.3894d12)    ! in Fb
c      common/leg_choice/leg
c      common/energy/S
c      common/distribution/xq
c      common/renor_scale/scale
c      common/usedalpha/AL,ge
c      common/param2/xmur
c      common/amass/am1,am2,amH,am4,am5
c      external Born_gg2H
c       
c        tau = amH**2/S
c         AL = alphasPDF(amH)
c
c         xa = yy(1)
c         xb = yy(2) 
c
c         sp = xa*xb*S
c        rsp = dsqrt(sp)
c
c        call kinvar1(xa,xb,p1,p2,p3)
c        scale = 2d0*dot(p1,p2)
c
c	eps = 0.5d0
c	Q2 = dsqrt(scale)
c	Qmin = amh - eps
c	Qmax = amh + eps
c
c	if(Q2 .ge. Qmin .and. Q2 .le. Qmax) then
c        xmuf=amH
c        xmur= xmuf
c        xmu2=xmuf**2
c
c        call pdf(xa,xmuf,f1)
c        call pdf(xb,xmuf,f2)
c        call setlum(f1,f2,xl)
c
c      sig= xl(2)*Born_gg2H(0,p1,p2,p3)
c
c
cc        pi_1 = 0.5d0*rsp
cc        flux = 4d0*pi_1*rsp
c        pi_1 = PI/amH
c        flux = 2d0*sp 
c
cc      xnorm=hbarc2/8d0/(2d0*Pi)**4/flux
c      xnorm=hbarc2*pi_1/flux
c
c      flo1_LO  = xnorm*sig/2d0/eps
c
c      else
c        flo1_LO= 0d0
c      endif
c      return
c      end
cc---------------------------------------------------------------------
c
c
ccC ~~~~~~~~~~~~~~~~~~~C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
c      integer function integrand(ndim,yy,ncomp,f,userdata,nvec ,core
c     .                            ,weight,iter)
c      implicit none 
cc     cuba specific parameters
c      integer ndim,ncomp,nvec,core,iter,userdata 
cc      real*8 xx(ndim) ,f(ncomp),weight,fnlo3,yy(10)
c      real*8 yy(2) ,f(ncomp),weight,flo1_LO,xx(10)
c      external flo1_LO
c       xx(1) = yy(1)
c       xx(2) = yy(2)
c      f(1)=flo1_LO(xx,weight)
cc	print*,"Fun "
cc	print*,xx
cc      f= xx(1)*xx(2)*xx(3)*xx(4)*xx(5)*xx(6)
c      return
c      end
cC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
