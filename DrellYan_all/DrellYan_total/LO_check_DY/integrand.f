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
      use openloops
      implicit double precision (a-h,o-z)
      dimension yy(10)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .          ,p_ex(0:3,1:4),c(1:2),Born(1:2),coef(2,-2:0),SumI(-2:0)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=0.3894d9)  ! in pb
c      parameter (hbarc2=0.3894d12)  ! in fb
      character*50 name
      common/pdfname/name
      common/leg_choice/leg
      common/energy/s
      common/distribution/xq
      common/renor_scale/scale
      common/usedalpha/AL,ge
      common/prc_id/id_LO,id_NLO_1r,id_NLO_1loop
      external Born_uU2eE
       
      rs  = dsqrt(s)
      xa     = yy(1)
      xb     = yy(2)

      sp = xa*xb*s
      rsp = dsqrt(sp)
        
      ipass = 0
        eps = 1.0d0
       xlow = xq - eps
      xhigh = xq + eps
!       xlow = 0d0
!      xhigh = rs 

      xcut = xq - 10.0d0

c      if (rsp .gt. xcut) then


        call kinvar2(yy,xinvmass,p1,p2,p3,p4)
        scale  = xinvmass
c       if (xinvmass .le. 0d0) goto 151 
        if ( scale .ge. xlow .and. scale .le. xhigh) then 
             
              xmuf=scale
              xmur=scale
              xmu2=xmuf**2
c              xmuf=xq
c              xmur=xq

cooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
c~~~~~~[ Openloops mat amp calculated from here ]~~~~~~~~c
cooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
c        call p1d_to_p2d_4(p1,p2,p3,p4,p_ex) 
c        call evaluate_tree(id_LO,p_ex,answer)
!	print*,answer
!	print*,Born_uU2eE(0,p1,p2,p3,p4)
!	print*,answer/Born_uU2eE(0,p1,p2,p3,p4)
!	print*
cooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

              call pdf(xa,xmuf,f1)
              call pdf(xb,xmuf,f2)
              call setlum(f1,f2,xl)
              AL = alphasPDF(xmur)
              call set_parameter("alpha_s",AL)

              sig= xl(1)*Born_uU2eE(0,p1,p2,p3,p4)
c              sig= xl(1)*answer

              xnorm=hbarc2/16d0/pi/(xa*xb*s)
              wgt=xnorm*sig*vwgt
              flo2_LO=wgt/vwgt/2d0/eps
c... for the inclusive case no binning is required.
!              flo2_LO=wgt/vwgt
       else                  
        flo2_LO=0d0
       endif
151        return
      end

c---------------------------------------------------------------------
