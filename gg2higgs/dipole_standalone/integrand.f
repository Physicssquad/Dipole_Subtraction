C ------------------------------------------------------------------ C 
c                        Our vegas specific
C -------------------------------------------------------------------- C      
      double precision function fnlo3(xx,weight)
      use openloops
      implicit double precision(a-h,o-z) 
      dimension  xx(10),p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),q(0:3)
     .             ,p(0:3,1:5),dip(27),B1(1:2),xl(15),f1(-6:6),f2(-6:6)
     .      ,p_ex(0:3,4),p12sum(0:3), p13sum(0:3), p14sum(0:3)
     .      , p23sum(0:3), p24sum(0:3),pmass(4)
c      integer i35,i45,is5,itest
      integer k1,k2,k3,ipass,n4,unphy
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+9)
      common/amass/am1,am2,am3,am4,am5
      common/energy/s
      common/set/set1
      common/distribution/xq
      common/countc/n4
      common/usedalpha/AL,ge
      common/scales/xinvmass
      common/counter/ifilter,itot_ev,iselect_scale
      common/counter_diff/diff,eps
      common/t_cuts/e_cut,t_cut
      common/caller/icall
      common/prc_id/id_LO,id_NLO_1r
      common/cuts_inPS/icut
      common/checks/ct,pcm,xlam
      COMMON /vegas_rn/ x1, x2, x3, x4
      external dipole_type_1_gg_g

        xa = xx(1)
        xb = xx(2)
        xc = xx(3)
c        x1 = xx(1)
c        x2 = xx(2)
c        x3 = xx(3)
c        x4 = xx(4)
c
c	print*,xx(1)
c	print*,xx(2)
c	print*,xx(3)
      xjac = 2d0
       amH = am3
        sp = xa*xb*s
c	do iboost=0,1
	iboost = 1
        call kinvar2_type_1(xa,xb,xc,xinvmass,p1,p2,p3,p4)
c	energy = dsqrt(s)
c	pmass(1) = am1
c	pmass(2) = am2
c	pmass(3) = am3
c	pmass(4) = am4
c        call  GET_MOMENTA(ENERGY,PMASS,P_ex)
c        call p2d_to_p1d_4(p_ex,p1,p2,p3,p4) 
c	if (p3(0) .ne. p3(0)  ) goto 151
c	call printframe8(p1,p2,p3,p4)

c	iboost=0
c        call kinvar2_type_1(iboost,xa,xb,xc,xinvmass,p1,p2,p3,p4)
c       call kinvar2_type_2(xa,xb,xc,xinvmass,p1,p2,p3,p4)

c Notatio nused : sij = 2d0*dot(pi,pj)
       do i = 0, 3
           p12sum(i) = p1(i) + p2(i)
           p13sum(i) = p1(i) + p3(i)
           p14sum(i) = p1(i) + p4(i)
           p23sum(i) = p2(i) + p3(i)
           p24sum(i) = p2(i) + p4(i)
       end do

       s12 = dot(p12sum, p12sum)
       s13 = dot(p13sum, p13sum)
       s14 = dot(p14sum, p14sum)
       s23 = dot(p23sum, p23sum)
       s24 = dot(p24sum, p24sum)

c	print*,"boost:",iboost
c        print*,"s12-s14-s24-amH**2",s12-s14-s24-amH**2
c	print*,"s12:",s12
c	print*,"s13:",s13
c	print*,"s14:",s14
c	print*,"s23:",s23
c	print*,"s24:",s24
c	print*,"p4**2:",dot(p4,p4)
c	print*,"p3**2:",dot(p3,p3)
c	print*,"p2**2:",dot(p2,p2)
c	print*,"p1**2:",dot(p1,p1)
c	print*

            sp = xa*xb*s
           rsp = dsqrt(sp)
         fnlo3 = 0d0
          amH = am3
          xmuf = amH/2d0
          xmur = xmuf 

c [ earlier without using technical cuts integration was giving NaN
c [ now it is working without using cuts, however phasespace factor needs to be checked properly.
c~~~~~~~~~~~~~[ Cuts ]~~~~~~~~~~~~~~~
        icol = 0
        coll = 1d-5
        soft = 1d-6

        if (s14 .lt. sp*coll) icol =1
        if (s24 .lt. sp*coll) icol =1
        e4 =(s12-amH**2)/2d0/rsp
        if (e4 .le. rsp*soft) icol = 1

c~~~~~~~~~~~~~[ Cuts ]~~~~~~~~~~~~~~~
c        if( sp  .ge. amH**2  .and. icol .eq. 0) then 
        if( sp  .ge. amH**2 ) then   ! I AM NOT USING CUTS HERE

          call pdf(xa,xmuf,f1)
          call pdf(xb,xmuf,f2)
          call setlum(f1,f2,xl)

          AL = alphasPDF(xmur)
c        call set_parameter("alpha_s",AL)

       call amp_mat_r(p1,p2,p3,p4,sig)
c~~~~~~[Openloops mat amp calculated from here ]~~~~~~~~c
c        call p1d_to_p2d_4(p1,p2,p3,p4,p_ex) 
c        call evaluate_tree(id_NLO_1r,p_ex,answer)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c

          SumD = dipole_type_1_gg_g(1,p1,p2,p3,p4) +
     .           dipole_type_1_gg_g(2,p1,p2,p3,p4)

c	per = (answer-sig)/answer*100
c	per1 = (SumD-sig)/sig*100
c	if (dabs(per) .gt. 0.01d0) then 
c     .    .and. 
cc         if ( dabs(per) .gt. 0.01d0) then
c	if (sig .ge. 100000d0) then
c	print*,"boost:",iboost
c	print*,"sigma :",sig
c	print*,"answer:",answer
c	print*,"SumD  :",SumD
c	print*,"s12   :",s12
c	print*,"s14   :",s14
c	print*,"s24   :",s24
c	print*,"xa,xb,xc",xa,xb,xc
c	print*,"ratio :",sig/SumD
c	print*
c	stop
c           sigma =xl(2)*(answer - SumD )
c        call MATDIP_CHECK(iflip, SIG, SUMD,'flip')
         sigma = xl(2) * (sig - SumD )
c         sigma = xl(2) * (answer - SumD )

          ! ref
c         xnorm = 2d0*dsqrt(xlam)/32d0/PI**2/sp
c~~~~~~~~~~~~~[Form Calc ] ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
          pout2 = (sp + am4**2 - am3**2)**2/4d0/sp - am4**2
           pout = dsqrt(pout2)
           pin2 = (sp + am2**2 - am1**2)**2/4d0/sp - am2**2 
           pin  = dsqrt(pin2)
          flux  = 4d0*pin*rsp
          xnorm = xjac* hbarc2*pout/16d0/PI**2/flux/rsp
c	print*,xnorm

c~~~~~~~~~~~~~~~~[ Original Normalization ]~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
          xmd = (am3**2 - am4**2)
           e3 = 0.5d0*(rsp + xmd/rsp)
             pf = dsqrt(e3**2 - am3**2)
            pfpi = 2d0*pf/dsqrt(s*xa*xb)
c          xnorm = hbarc2*pfpi/16d0/pi/(s*xa*xb)
c          print*,"ratio:",xnorm1/xnorm
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c

c          wgt=xajac*xnorm*sigma*weight*2d0* amh/xa/S
c          wgt=xajac*xnorm*sigma*weight

            wgt = xnorm*sigma*weight
          fnlo3 = wgt/weight
	if (fnlo3 .ne. fnlo3 ) goto 151
        else
           fnlo3  = 0d0
        endif
c	enddo
c	stop
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


