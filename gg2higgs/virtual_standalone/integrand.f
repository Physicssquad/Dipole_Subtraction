      function flo2_Vir(yy,vwgt)
      use openloops
      implicit double precision (a-h,o-z)
      dimension yy(10)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .          ,p(0:3,1:4),c(1:2),coef(2,-2:0),SumI(-2:0)
     .          ,p_ex(0:3,1:3),ans_loop(0:2),coefI(12,-2)
      dimension am211(0:2),amp2cc(3) 
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
c      common/usedalpha/AL,ge
      common/amass/am1,am2,amH,am4,am5        
      common/prc_id/id_LO,id_NLO_1r,id_NLO_1loop
      external Born_gg2h

         tau = amH**2/S
       xajac = 1d0 - tau
          xa = tau + xajac*yy(1)
          xb = tau/xa

         sp  = xa*xb*s
         rsp = dsqrt(sp)
        
              xmuf=amH/2d0
              xmur=xmuf
              xmu2=xmuf**2

              call pdf(xa,xmuf,f1)
              call pdf(xb,xmuf,f2)
              call setlum(f1,f2,xl)
              AL = alphasPDF(xmur)
          call set_parameter("alpha_s",AL)
               as= AL/4d0/Pi
             call kinvar1(xa,xb,p1,p2,p3)
              s12 = 2d0*dot(p1,p2)
             xlqr = dlog(s12/xmu2)
        
        call p1d_to_p2d_3(p1,p2,p3,p_ex) 
c	call Iterm(p_ex,id_LO,coefI,SumI,AL)
      call evaluate_tree(id_LO,p_ex,answer)

      call evaluate_cc(id_LO,p_ex,amp_210,amp2cc,amp2ewcc)
      call evaluate_loop(id_NLO_1loop,p_ex,am2l0,am211, acc)
c0000000000000000000000000000000000000000000000
c... this is colour correlated born, with <M |Ta . Tb| M>

c... ColourFactor = Tb^2 so overall there is|     |Ta . Tb|   |
c... 					    |<M	---------- M> |
c... 					    |        Tb^2     |

        answer = amp2cc(1)
        All_legs = 1d0
        ColourFactor = 3.0d0

        coefeps2=-All_legs*AL*answer*CA/2d0/Pi/ColourFactor

        coefeps1 = 
     -    All_legs*(Al*answer*(-11*CA + 6*CA*EulerGamma + 4*Nf*TR - 
     -       12*CA*Log(2d0) - 6*CA*Log(Pi) - 
     -       6*CA*Log(xmu2/s12)))/Pi/12d0/ColourFactor   

             xlog2 = dlog(2d0)
             xlogpi = dlog(Pi)
            xlogmu2bs12 = dlog(xmu2/s12) 
        coefeps0=-All_legs*(Al*answer*(-200d0*CA + 66d0*CA*EulerGamma- 
     -       18d0*CA*EulerGamma**2 + 21d0*CA*Pi**2 + 64d0*Nf*TR - 
     -       24d0*EulerGamma*Nf*TR - 132d0*CA*xlog2 + 
     -       72d0*CA*EulerGamma*xlog2 + 48*Nf*TR*xlog2 - 
     -       72d0*CA*xlog2**2 - 66d0*CA*xlogpi + 
     -       36d0*CA*EulerGamma*xlogpi + 24d0*Nf*TR*xlogpi - 
     -       72d0*CA*xlog2*xlogpi - 18d0*CA*xlogpi**2 - 
     -       66d0*CA*xlogmu2bs12 + 
     -       36d0*CA*EulerGamma*xlogmu2bs12 + 
     -       24d0*Nf*TR*xlogmu2bs12 - 
     -       72d0*CA*xlog2*xlogmu2bs12 - 
     -       36d0*CA*xlogpi*xlogmu2bs12 - 
     -       18d0*CA*xlogmu2bs12**2))/Pi/72d0/ColourFactor


!	print*,"epsm0",coefeps0,am211(0)
!	print*,"epsm1",coefeps1,am211(1)
!	print*,"epsm2",coefeps2,am211(2)
!	print*,"epsm0",coefeps0/am211(0)
!	print*,"epsm1",coefeps1/3d0/am211(1)
!	print*,"epsm2",coefeps2/am211(2)
!	stop

            virtualPLUSeikonal = (am211(0)+coefeps0)*2d0/3d0

            pi_1 = PI/amH
            flux = 2d0*sp
           xnorm = xajac*hbarc2*pi_1/flux
           
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
