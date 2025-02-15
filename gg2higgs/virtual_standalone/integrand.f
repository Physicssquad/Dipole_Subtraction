      function flo2_Vir(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .          ,p(0:3,1:4),c(1:2),coef(2,-2:0)
     .          ,p_ex(0:3,1:3),ans_loop(0:2)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+9)
      parameter (cf=4d0/3d0)
      parameter (zeta2=1.64493406684823d0)
      parameter (N=3)
      parameter (Nf=5)
      parameter (Ca=3)
      parameter (Tf=0.5D0) 
      parameter (EulerGamma = 0.5772156649015328606065120d0)

      character*50 name
      common/pdfname/name
      common/leg_choice/leg
      common/energy/s
      common/distribution/xq
      common/renor_scale/scale
      common/usedalpha/AL,ge
      common/amass/am1,am2,amH,am4,am5        
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
               as= AL/4d0/Pi
             call kinvar1(xa,xb,p1,p2,p3)
             Born = Born_gg2h(0,p1,p2,p3)
              s12 = 2d0*dot(p1,p2)
             xlqr = dlog(s12/xmu2)

c         energy = 13000
c         alpha_s = AL
c         call p1d_to_p2d_3(p1,p2,p3,p_ex) 
c         call main(p_ex,alpha_s,energy,xmuf,answer,ans_loop)

c        virtualPLUSeikonal= (-11d0*CA)/3d0+2d0*CA*xlqr - 
c     -     Ca*xlqr**2 + (2d0*Nf)/3d0 + 
c     -  (7d0*Ca*Pi**2)/6d0
      virtualPLUSeikonal =(-100.*CA)/9. + (11.*CA*xlqr)/3.-2*CA*xlqr**2+
     -  (16.*Nf)/9. - (2.*xlqr*Nf)/3. + (7.*CA*Pi**2)/3.
            virtualPLUSeikonal=as*Born*virtualPLUSeikonal
c	print*,"finite:",virtualPLUSeikonal
c          print *, "evaluate_loop"
c          print *, "Tree:       ",answer 
c          print *, "Loop ep^0:  ",ans_loop(0)
c          print *, "ep^0 ratio: ",ans_loop(0)/virtualPLUSeikonal
c          print *, "Loop ep^-1: ",ans_loop(1)
c          print *, "Loop ep^-2: ",ans_loop(2)

            pi_1 = PI/amH
            flux = 2d0*sp
           xnorm = xajac*hbarc2*pi_1/flux
           
        sig= xl(2)*xnorm*virtualPLUSeikonal
       
         flo2_Vir= sig * 2d0*amH/xa/S

      return
      end

c---------------------------------------------------------------------

C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
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
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
