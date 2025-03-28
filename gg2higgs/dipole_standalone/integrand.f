C ------------------------------------------------------------------ C 
c                        Our vegas specific
C -------------------------------------------------------------------- C      
      double precision function fnlo3(xx,weight)
      use openloops
      implicit double precision(a-h,o-z) 
      dimension  xx(10),p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),q(0:3)
     .             ,p(0:3,1:5),dip(27),B1(1:2),xl(15),f1(-6:6),f2(-6:6)
     .      ,p_ex(0:3,4),p12sum(0:3), p13sum(0:3), p14sum(0:3)
     .      , p23sum(0:3), p24sum(0:3),pmass(4),ptilde(0:3,3,2)
     .         ,p1til(0:3),p2til(0:3),p3til(0:3)
c      integer i35,i45,is5,itest
      integer k1,k2,k3,ipass,n4,unphy
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+9)
      common/scales/xmuf,xmur
      common/amass/am1,am2,am3,am4,am5
      common/energy/s
      common/set/set1
      common/countc/n4
      common/usedalpha/AL,ge
      common/scales/xinvmass
      common/counter/ifilter,itot_ev,iselect_scale
      common/counter_diff/diff,eps,isame,idiff
      common/t_cuts/e_cut,t_cut
      common/caller/icall
      common/prc_id/id_LO,id_NLO_1r
      common/cuts_inPS/icut
      common/checks/ct,pcm,xlam
      COMMON /vegas_rn/ x1, x2, x3, x4
      common/distribution/rapidity,p_T
      external dipole_type_1_gg_g

        xa = xx(1)
        xb = xx(2)
        xc = xx(3)
        p_T= 100d0

      xjac = 2d0
       amH = am3
        sp = xa*xb*s
        if (sp  .ge. amH**2 ) then !...basic condition 
        call kinvar2_type_1(xa,xb,xc,xinvmass,p1,p2,p3,p4)
        call reduceps(p1,p2,p3,p4,ptilde)
c_______________________________________________________c
c...real matrix element contribution
c        call Fjm_plus_one(p1,p2,p3,p4,jet_cut)
c        if ( jet_cut .eq. 1 ) goto 151 
c_______________________________________________________c
c...dipole contribution 
c        call Fjm(ptilde,jet_cut)
c      s12 = 2d0*dot2(ptilde,1,2,1) 
c      s14 = 2d0*dot2(ptilde,1,4,1) 
c      s24 = 2d0*dot2(ptilde,2,4,1) 
c      print*,"m-body",s14 + s24 - s12
c        if ( jet_cut .eq. 0 ) goto 151 

       s12 = 2d0*dot(p1,p2)
       s13 = 2d0*dot(p1,p3)
       s14 = 2d0*dot(p1,p4)
       s23 = 2d0*dot(p2,p3)
       s24 = 2d0*dot(p2,p4)
       s34 = 2d0*dot(p3,p4)

         rsp = dsqrt(sp)
         fnlo3 = 0d0
          amH = am3

c~~~~~~~~~~~~~[ Cuts ]~~~~~~~~~~~~~~~
        icol = 0
        iprint = 0
        coll = 1d-5
        soft = 1d-2
        if (s14/s12 .lt. coll) icol =1
        if (s24/s12 .lt. coll) icol =1
        e4 =(s12-amH**2)/2d0/rsp
        if (e4 .le. rsp*soft) icol = 1
        pt_higgs = dsqrt(p3(1)**2 + p3(2)**2)

        call pdf(xa,xmuf,f1)
        call pdf(xb,xmuf,f2)
        call setlum(f1,f2,xl)

        AL = alphasPDF(xmur)
        call set_parameter("alpha_s",AL)

        call amp_mat_r(p1,p2,p3,p4,sig)
c~~~~~~[ Openloops mat amp calculated from here ]~~~~~~~~c
        call p1d_to_p2d_4(p1,p2,p3,p4,p_ex) 
        call evaluate_tree(id_NLO_1r,p_ex,answer)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c

        SumD = dipole_type_2_gg_g(1,p1,p2,p3,p4) 
     .          +dipole_type_2_gg_g(2,p1,p2,p3,p4)

        SumD2 = dipole_type_1_gg_g(1,p1,p2,p3,p4) 
     .           +dipole_type_1_gg_g(2,p1,p2,p3,p4)

         sigma = xl(2) * (sig - SumD )

c~~~~~~~~~~~~~[Form Calc ] ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
c          pout2 = (sp + am4**2 - am3**2)**2/4d0/sp - am4**2
c           pout = dsqrt(pout2)
c           pin2 = (sp + am2**2 - am1**2)**2/4d0/sp - am2**2 
c           pin  = dsqrt(pin2)
c          flux  = 4d0*pin*rsp
c          xnorm = xjac* hbarc2*pout/16d0/PI**2/flux/rsp
c~~~~~~~~~~~~~~~~[ Original Normalization ]~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
          xmd = (am3**2 - am4**2)
           e3 = 0.5d0*(rsp + xmd/rsp)
             pf = dsqrt(e3**2 - am3**2)
            pfpi = 2d0*pf/dsqrt(s*xa*xb)
          xnorm = hbarc2*pfpi/16d0/pi/(s*xa*xb)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c

c          wgt=xajac*xnorm*sigma*weight*2d0* amh/xa/S
c          wgt=xajac*xnorm*sigma*weight

      wgt = xnorm*sigma*weight
      fnlo3 = wgt/weight
c	if (fnlo3 .ne. fnlo3 ) goto 151
        else
           fnlo3  = 0d0
        endif
151   return
      end
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
      integer function integrand(ndim,xx,ncomp,f,userdata,nvec ,core
     .                            ,weight,iter)
      implicit none 
c     cuba specific parameters
      integer ndim,ncomp,nvec,core,iter,userdata 
c      real*8 xx(ndim) ,f(ncomp),weight,fnlo3,yy(10)
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


