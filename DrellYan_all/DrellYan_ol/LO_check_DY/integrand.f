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
     .          ,am211(0:2),am211_ir(0:2),amp2cc(6)
c... this is for the evaluate_poles
     .          ,xm2bare(0:2),xm2ct(0:2),xm2ir(0:2),xm2sum(0:2)
       parameter (pi=3.14159265358979d0)
c      parameter (hbarc2=0.3894d9)  ! in pb
      parameter (hbarc2=0.3894d12)  ! in fb
      parameter (EulerGamma = 0.5772156649015328606065120d0)

      character*50 name
      common/pdfname/name
      common/leg_choice/leg
      common/energy/s
      common/distribution/xq
      common/renor_scale/scale
      common/prc_id/id_LO,id_NLO_1r,id_NLO_1loop
      common/usedalpha/AL,ge
      external Born_uU2eE
       
      rs  = dsqrt(s)
      xa     = yy(1)
      xb     = yy(2)

      rsp = dsqrt(xa*xb*s)
        
      ipass = 0
        eps = 1.0d0
       xlow = xq - eps
      xhigh = xq + eps

      xcut = xq - 10.0d0

c      if (rsp .gt. xcut) then


        call kinvar2(yy,xinvmass,p1,p2,p3,p4)
        s12 = 2d0*dot(p1,p2)
         scale  = xinvmass
c      if (scale .ge. 0d0) then
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
c~~~~~~[ Openloops mat amp calculated from here ]~~~~~~~~c
        call set_parameter("alpha_s",AL)
!        call set_parameter("r2_on",1)
!        call set_parameter("ct_on",1)
        call set_parameter("polenorm",1)
	call set_parameter("ckmorder",5)
        call set_parameter("truepoles_on",1)
!        call set_parameter("iop_on",1)
        call set_parameter("mu", xmuf)
c        call set_parameter("mu",dsqrt(s12))
        call p1d_to_p2d_4(p1,p2,p3,p4,p_ex)

        call evaluate_tree(id_LO,p_ex,born_ol)
        call evaluate_cc(id_LO,p_ex,amp_210,amp2cc,amp2ewcc)
        call evaluate_loop(id_NLO_1loop,p_ex,am2l0,am211, acc)
        call evaluate_iop(id_NLO_1loop,p_ex,am2l0_ir,am211_ir)

        call evaluate_poles(id_NLO_1loop,p_ex,xm210,xm2bare,xm2ct,xm2ir
     .   ,xm2sum)
	print*,"Evaluate Poles"
	print*,"xm210",xm210
	print*,"xm2ba",xm2bare
        print*,"xm2ct",xm2ct
        print*,"xm2ir",xm2ir
        print*,"xSumr",xm2Sum
	print*
c~~~~~~[            ~~~~~~~~~~~~~~              ]~~~~~~~~c

c... Here I will construct the Ikonal part
c... for DY case.
              Cf = 4d0/3d0
!	xmu2 = s12
c... with the same scale result is not matching just finite part.

c... notice a factor of 2 in all the coefficients, it is for the
c... information on legs which is summed over.
      coefepsm2 = -2d0* Al *Cf/2d0/Pi *amp2cc(1)/Cf

      coefepsm1 = -(Al*(3d0*Cf + 2d0*Cf*dLog(xmu2/s12)))/Pi/4d0 
      coefepsm1 = 2d0*coefepsm1*amp2cc(1)/Cf
c... eps(-2) and eps(-1) are matching with Eikonal openloops

      coefeps0  = (Al*(-60d0*Cf + 7d0*Cf*Pi**2 - 
     -      18d0*Cf*dLog(xmu2/s12) - 6d0*Cf*dLog(xmu2/s12)**2))/
     -  (24d0*Pi)
      coefeps0 = 2d0*coefeps0*amp2cc(1)/Cf

c... with the euler gamma
!       coefeps0 = (Al*(-60d0*Cf + 18d0*Cf*EulerGamma -   
!     -      6d0*Cf*EulerGamma**2 + 7d0*Cf*Pi**2 -
!     -      18d0*Cf*dLog(xmu2/s12) +
!     -      12d0*Cf*EulerGamma*dLog(xmu2/s12) -
!     -      6d0*Cf*dLog(xmu2/s12)**2))/(24d0*Pi)
!        coefeps0 =2d0*coefeps0*amp2cc(1)/3d0
!	print*,"withEG:",coefeps0


c... this finite is not matching with the openloops eikonal.


c... Ioperator
      coefeps0_iop = Al/Pi*Cf*(18d0*xlqr - 6d0*xlqr**2 + 7d0*Pi**2)/12d0
      coefeps0_iop = 2d0*coefeps0_iop*amp2cc(1)/Cf


c	print*,"1loop:",am211,acc
c        print*,"Iop  :",am211_ir
c	print*,"coef :",coefeps0,coefepsm1,coefepsm2
c	print*,"coefI:",coefeps0_iop
c	print*,"ratio:",coefeps0/am211_ir(0)
c	print*,"ratio:",am211_ir(-2)/coefepsm2
c	print*
c	stop
c        print*,"1lp:",am211(2)/born_ol,am211(1)/born_ol,am211(0)/born_ol
c        print*,"born",born_ol,AL
c        print*
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c

              sig= xl(1)*Born_uU2eE(0,p1,p2,p3,p4)

              xnorm=hbarc2/16d0/pi/(xa*xb*s)
              wgt=xnorm*sig*vwgt
              flo2_LO=wgt/vwgt/2d0/eps
            return
       else                  
        flo2_LO=0d0
        return
       endif
      return
      end

c---------------------------------------------------------------------
