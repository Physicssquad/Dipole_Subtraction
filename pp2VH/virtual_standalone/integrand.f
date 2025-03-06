      function flo2_Vir(yy,vwgt)
      use openloops
      implicit double precision (a-h,o-z)
      dimension yy(10)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .          ,p(0:3,1:4),c(1:2),coef(2,-2:0),SumI(-2:0)
     .          ,p_ex(0:3,1:3),ans_loop(0:2),coefI(12,-2),am2Rct(0:2)
      dimension am211(0:2),amp2cc(3),am211_ir(0:2)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+9)
      parameter (cf=4d0/3d0)
      parameter (zeta2=1.64493406684823d0)
      parameter (N=3)
      parameter (Nf=5)
      parameter (Ca=3d0)
      parameter (Tf=0.5D0) 
      parameter (Tr=0.5D0) 
      parameter (EulerGamma = 0.5772156649015328606065120d0)

      character*50 name
      common/pdfname/name
      common/leg_choice/leg
      common/energy/s
      common/distribution/xq
      common/renor_scale/scale
      common/amass/am1,am2,amH,am4,am5        
      common/prc_id/id_LO,id_NLO_1r,id_NLO_1loop
      external Born_gg2h_

         tau = amH**2/S
       xajac = 1d0 - tau
          xa = tau + xajac*yy(1)
          xb = tau/xa

         sp  = xa*xb*s
         rsp = dsqrt(sp)
        
              xmuf=amH/2d0
              xmur=xmuf
              xmu2=xmur**2

              call pdf(xa,xmuf,f1)
              call pdf(xb,xmuf,f2)
              call setlum(f1,f2,xl)
              AL = alphasPDF(xmur)
        call kinvar1(xa,xb,p1,p2,p3)
              s12 = 2d0*dot(p1,p2)
        call set_parameter("alpha_s",AL)
        call set_parameter("mu", xmur)
!        call set_parameter("check_poles",1)
        call set_parameter("iop_on",1)
!        call set_parameter("ckmorder",6)
!        call set_parameter("polenorm",1)
        
        call p1d_to_p2d_3(p1,p2,p3,p_ex) 
c	call Iterm(p_ex,id_LO,coefI,SumI,AL)
!        call set_parameter("r2_on",0)
!       call set_parameter("ct_on",0)

      call evaluate_tree(id_LO,p_ex,born)
      call evaluate_cc(id_LO,p_ex,amp_210,amp2cc,amp2ewcc)
      call evaluate_loop(id_NLO_1loop,p_ex,am2l0,am211, acc)

      call evaluate_r2(id_NLO_1loop,p_ex,aborn,am2r2)
      call evaluate_loopct(id_NLO_1loop,p_ex,aborn,am2Rct)
c	call get_parameter("ckmorder",l)
c	print*,l

c	print*
c	print*,born
c	print*,am211(2),am211(1),am211(0)
c	print*,am211(2)/born,am211(1)/born,am211(0)/born
c	print*,"Al,xlqr:",Al,xlqr
c	print*,"s12",s12

c... this is colour correlated born, with <M |Ta . Tb| M>

c... ColourFactor = Tb^2 so overall there is|     |Ta . Tb|   |
c... 					    |<M	---------- M> |
c... 					    |        Tb^2     |
        answer = amp2cc(1)
        All_legs = 2d0
        ColourFactor = 3.0d0

c... Cartani Iikonal operator withoit Euler gamma
        coefeps2=-All_legs*AL*answer*CA/2d0/Pi/ColourFactor

        coefeps1 =-All_legs*answer*(Al*(11d0*CA - 4d0*Nf*TR +
     .             6d0*CA*dLog(xmu2/s12)))/12d0/Pi/3d0

!        coefeps0 = (Al*(-200d0*CA + 21d0*CA*Pi**2 + 64d0*Nf*TR 
!     -             -66d0*CA*dLog(xmu2/s12) +
!     -   24d0*Nf*TR*dLog(xmu2/s12) 
!     -     - 18d0*CA*dLog(xmu2/s12)**2))/(72d0*Pi)
!        coefeps0 =All_legs*coefeps0*answer/3d0

c... I operator with the born eps added only finite piece.
!        coefeps0 = (Al*(20d0*Nf + CA*(-134d0 + 21d0*Pi**2)
!     .      - 6d0*(5d0*CA - 2d0*Nf)*dLog(xmu2/s12) - 
!     .     18d0*CA*dLog(xmu2/s12)**2))/(72d0*Pi)
!       coefeps0 = All_legs*coefeps0*answer/3d0

!        coefeps0 = (Al*(-60d0*Cf + 18d0*Cf*EulerGamma - 
!     -      6d0*Cf*EulerGamma**2 + 7d0*Cf*Pi**2 - 
!     -      18d0*Cf*dLog(xmu2/s12) + 
!     -      12d0*Cf*EulerGamma*Log(xmu2/s12) - 
!     -      6d0*Cf*Log(xmu2/s12)**2))/(24d0*Pi)
c... coef0
!       coefeps0 = 3d0*Cf/2d0 - Cf*xlqr**2/2d0 + Cf* Pi**2/12d0
!       coefeps0 = Al/2d0* coefeps0
!       coefeps0 = (Al*(3d0*Pi**2 + 66d0*dLog(xmu2/s12) -
!     .     4d0*Nf*dLog(xmu2/s12) - 18d0*dLog(xmu2/s12)**2))/(48d0*Pi)
!       coefeps0 = All_legs*coefeps0*answer/3d0

c... Gehrmann Finite expression.!
        xlqr = dLog(xmu2/s12) 
        zeta = Pi/6d0
        xzeta2= zeta*zeta 
        coefeps0 = Al*Ca/4d0/Pi*(-lqr**2 +7d0* xzeta2 ) 
        coefeps0 = coefeps0*answer/3d0

	print*,"eps0 ",coefeps0,am211(0),coefeps0/am211(0)
	print*,am211(0)-am2Rct(0)
	print*,am211(0)
	stop

c	print*
c!	print*,xlqr,born,AL
c	print*,"epsm1",coefeps1,am211(1)
c	print*,"epsm2",coefeps2,am211(2)
c	print*
c	stop
c	print*,"epsm0bborn:",am211(0)/born
c	print*,"epsm2bborn:",am211(2)/born
c	print*,"epsm1bborn:",am211(1)/born

!	print*,"epsm0",coefeps0/am211(0)
c	print*,"epsm1bborn:",coefeps1/am211(1)
!	print*,"epsm2",coefeps2/am211(2)
!	stop
c... for catani dipole case Eikonal is added 
c            virtualPLUSeikonal = (am211(0) - coefeps0)
            virtualPLUSeikonal = (am211(0) )!- coefeps0)
c        ________________________________________
            pi_1 = PI/amH
            flux = 2d0*sp
           xnorm = xajac*hbarc2*pi_1/flux
c        ________________________________________
           
        sig= xl(2)*xnorm*virtualPLUSeikonal
       
         flo2_Vir= sig * 2d0*amH/xa/S

      return
      end
c-----------------------------------------------------------------

C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
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
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
