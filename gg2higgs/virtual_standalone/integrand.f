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
      parameter (Nf=6)
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

!        call set_parameter("check_poles",1)
!        call set_parameter("ckmorder",6)
!        call set_parameter("polenorm",1)
        
      call p1d_to_p2d_3(p1,p2,p3,p_ex) 
      call evaluate_tree(id_LO,p_ex,born)
      call evaluate_cc(id_LO,p_ex,amp_210,amp2cc,amp2ewcc)
c... this is a IR subtracted result.
      call set_parameter("ct_on", 1)
      call set_parameter("iop_on",1)
      call set_parameter("r2_on", 1)
      call evaluate_loop(id_NLO_1loop,p_ex,am2l0,am211, acc)

!      call evaluate_loopct(id_NLO_1loop,p_ex,aborn,am2Rct)
c...  1-loop finite = bare_1loop + UV_ct + (R2 or IR)

      virtualPLUSeikonal = am211(0) !this is the IR subtracted term from Openloops 

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
