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
      common/final_data_common/final_data
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
             
              xmuf=scale
              xmur=scale
              xmu2=xmuf**2

              call pdf(xa,xmuf,f1)
              call pdf(xb,xmuf,f2)
              call setlum(f1,f2,xl)
              AL = alphasPDF(xmur)

              sig= xl(1)*Born_uU2eE(0,p1,p2,p3,p4)

              xnorm=hbarc2/16d0/pi/(xa*xb*s)
              wgt=xnorm*sig*vwgt
              flo2_LO=wgt/vwgt/2d0/eps
              call histogram(xnorm*sig,scale,vwgt,final_data)
      return
      end

c---------------------------------------------------------------------
      subroutine histogram(delI,scale,vwgt,send_to_main)
      implicit none
      double precision delI,scale,eps,xlow,xhigh,xq,send_to_main,vwgt
      integer condition
      double precision, save :: final_data = 0.0

      xq = 500d0
      eps = 0.5d0
      xlow = xq - eps
      xhigh = xq + eps
      if ( scale .ge. xlow .and. scale .le. xhigh) then 
      final_data = final_data + delI*vwgt
      send_to_main = final_data/2d0/eps
      endif
      end

