C ------------------------------------------------------------------ C 
c                        Our vegas specific
C -------------------------------------------------------------------- C      
      double precision function fnlo3(xx,weight)
      use openloops
      include 'header.h'
      include 'parameter.h'

      dimension  xx(10),p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
      dimension  p_ex(0:3,1:5),xl(15),f1(-6:6),f2(-6:6),ptilde(0:3,4,2)
      integer k1,k2,k3,ipass,n4,unphy

      common/amass/am1,am2,am3,am4,am5
      common/energy/s
      common/usedalpha/AL,ge
      common/prc_id/id_LO,id_NLO_1r

      external dipole_type_1_gg_g

      xjac = 2d0
        sp = xa*xb*s

       call kinvar3(xx,xxjac,xinvmass,p1,p2,p3,p4,p5,unphy)

c... this is reduced ps for CS phasespace parametrisation
       call reduceps_5to4(p1,p2,p3,p4,p5,ptilde)
       call retrive_tilde(1,ptilde,p_ex)
       print*,k,ptilde
       print*,k,p_ex
       stop

c... these jet functions will define the observable
c... real matrix element contribution
!        call Fjm_plus_one(p1,p2,p3,p4,jet_cut)

c... dipole contribution 
!        call Fjm(ptilde,jet_cut)

!      s12 = 2d0*dot2(ptilde,1,2,1) 
!      s14 = 2d0*dot2(ptilde,1,4,1) 
!      s24 = 2d0*dot2(ptilde,2,4,1) 

       s12 = 2d0*dot(p1,p2)
       s13 = 2d0*dot(p1,p3)
       s14 = 2d0*dot(p1,p4)
       s23 = 2d0*dot(p2,p3)
       s24 = 2d0*dot(p2,p4)
       s34 = 2d0*dot(p3,p4)

         rsp = dsqrt(sp)
         fnlo3 = 0d0
          amH = am3
          xmuf = amH/2d0
          xmur = xmuf 

          call pdf(xa,xmuf,f1)
          call pdf(xb,xmuf,f2)
          call setlum(f1,f2,xl)

          AL = alphasPDF(xmur)

c~~~~~~[ Openloops mat amp calculated from here ]~~~~~~~~c
        call set_parameter("alpha_s",AL)
        call p1d_to_p2d_4(p1,p2,p3,p4,p_ex) 
        call evaluate_tree(id_NLO_1r,p_ex,answer)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c

         sigma = xl(2) * (answer - SumD )

c~~~~~~~~~~~~~~~~[ Original Normalization ]~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
          xmd = (am3**2 - am4**2)
           e3 = 0.5d0*(rsp + xmd/rsp)
             pf = dsqrt(e3**2 - am3**2)
            pfpi = 2d0*pf/dsqrt(s*xa*xb)
          xnorm = hbarc2*pfpi/16d0/pi/(s*xa*xb)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c

            wgt = xnorm*sigma*weight
          fnlo3 = wgt/weight

151   return
      end
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
      integer function integrand(ndim,xx,ncomp,f,userdata,nvec ,core
     .                            ,weight,iter)
      implicit none 
c     cuba specific parameters
      integer ndim,ncomp,nvec,core,iter,userdata 
!      real*8 xx(ndim) ,f(ncomp),weight,fnlo3,yy(10)
      real*8 xx(6) ,f(ncomp),weight,fnlo3,yy(10)
      external fnlo3
       yy(1) = xx(1)
       yy(2) = xx(2)
       yy(3) = xx(3)
       yy(4) = xx(4)
       yy(5) = xx(5)
       yy(6) = xx(6)
       yy(7) = 0d0 
       yy(8) = 0d0 
       yy(9) = 0d0 
       yy(10)= 0d0 

      f(1)=fnlo3(yy,weight)

      return
      end
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        


