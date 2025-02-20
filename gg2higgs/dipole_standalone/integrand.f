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
          xmuf = amH/2d0
          xmur = xmuf 

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

c        if (pt_higgs .le. 50d0 ) goto 151 
c        if (pt_higgs .ge. 0.0d0 .and. pt_higgs .le. 100d0 ) then
c        if (pt_higgs .le. 1d0 ) then
c          continue
c        else
c        goto 151
c        endif
c	iprint = 0
c	if (2d0*s14 .ge. s12 ) icol = 1 
c	if (2d0*s24 .ge. s12 ) icol = 1 
c~~~~~~~~~~~~~[ Cuts ]~~~~~~~~~~~~~~~
c        if (icol .eq. 1 ) goto 151

c        print*,s14,s12*coll
c        print*,s24,s12*coll

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
c	idiff = idiff + 1
c	if (sig .lt. SumD ) goto 151
c	idiff = idiff - 1
c	isame = isame + 1
c	if (icol .eq. 0 ) iprint = 1
c	if ( s24/s12 .le. 1d-2 ) goto 151
c	if ( s14/s12 .le. 1d-2 ) goto 151

c	if (answer .le. 5000d0 ) goto 151 
c	if (2d0*s14 .ge. s12 ) iprint = 1 
c	if (2d0*s24 .ge. s12 ) iprint = 1 
c        if (answer .le. 500d0 ) goto 151

c	per = (answer-sig)/answer*100
c	per1 = (SumD-sig)/sig*100
c	if (dabs(per) .gt. 0.01d0) then 
c     .    .and. 
cc         if ( dabs(per) .gt. 0.01d0) then
c	if (dabs(answer-SumD)/answer .le. 1d-4) then
c	if (icol .eq. 0 ) then
c	if (dabs(answer-SumD)/answer*100 .ge. 99.9d0) then   !...check max diff between dipole of matrix
c        idiff = idiff + 1 
c	if (dabs(answer-SumD)/answer*100 .le. 1d-4) then
c	if (dabs(sig-SumD2)/sig*100 .ge. 250d0) then
c          print*,"MT",answer,SumD
c          print*,"NK",sig,SumD2
c	  print*
c	endif
c        print*,dabs(answer-SumD)/answer*100
c	if (dabs(sig-SumD2)/SumD2 .ge. 1d+4) then
c	print*,"~~~~~ Event ~~~~~~~"
c         alphaMin = 0.8d0
c         vtildei1=dot(p1,p4)/dot(p1,p2)
c         vtildei2=dot(p2,p4)/dot(p1,p2)
c          ips1 = 0
c          ips2 = 0
c         if (vtildei1 .lt. alphaMin) ips1 = 1 
c         if (vtildei2 .lt. alphaMin) ips1 = ips + 1 

c	if ( s14/s12 .gt. 0.50d0 .or. s24/s12 .gt. 0.50d0 ) then
c        fnlo3 =0d0
c        if ( ips1 .ge. 1 ) goto 151 
c	idiff = idiff + 1
c        if ( dabs(answer-SumD)/answer*100 .gt. 80 ) then
c        if ( dabs(sig-SumD2)/sig*100 .gt. 80 ) then
c	idiff = idiff - 1
c        isame = isame + 1 
c	endif
c        if (dabs(sumD - SumD2 )/dabs(SumD) .ge. 1d+7 ) iprint =1
c        if (dabs(sumD - SumD2 )/dabs(SumD) .ge. 1d+7 ) iprint =1
c	if (xc .ge. 0.99999d0 ) iprint = 1
        if (iprint .eq. 1 ) then
c          print*,dabs(sumD - SumD2 )/dabs(SumD)
c          SumD = dipole_type_3_gg_g(1,p1,p2,p3,p4) 
c     .          +dipole_type_3_gg_g(2,p1,p2,p3,p4)
        print*,"sigma_NK  :",sig
	print*,"SumDip_NK :",SumD2
	print*,"1real_amp :",answer
	print*,"SumDipole :",SumD

c	print*,"ratio     :",answer/SumD
c	print*,"s12       :",s12
c	print*,"s14       :",s14
c	print*,"s24       :",s24
c	print*,"p4(0)     :",p4(0)
	print*

c	print*,"s24       :",s24
c	print*,"s14       :",s14
	print*,"s24/s12   :",s24/s12
	print*,"s14/s12   :",s14/s12
	print*,"xa,xb,xc",xa,xb,xc
	print*
	endif ! iprint
c	print*,"~~~~~Event~~~~~~~"
c	print*," SumD/SumD2 ratio     :",SumD/SumD2
c	stop
c	endif
c           sigma =xl(2)*(answer - SumD )
c        call MATDIP_CHECK(iflip, SIG, SUMD,'flip')
c         sigma = (answer - SumD )
c         sigma = xl(2) * (answer - SumD )
c         sigma = xl(2) * (answer - SumD )
         sigma = xl(2) * (sig - SumD )
          ! ref
c         xnorm = 2d0*dsqrt(xlam)/32d0/PI**2/sp
c~~~~~~~~~~~~~[Form Calc ] ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
c          pout2 = (sp + am4**2 - am3**2)**2/4d0/sp - am4**2
c           pout = dsqrt(pout2)
c           pin2 = (sp + am2**2 - am1**2)**2/4d0/sp - am2**2 
c           pin  = dsqrt(pin2)
c          flux  = 4d0*pin*rsp
c          xnorm = xjac* hbarc2*pout/16d0/PI**2/flux/rsp
c	print*,xnorm

c~~~~~~~~~~~~~~~~[ Original Normalization ]~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
          xmd = (am3**2 - am4**2)
           e3 = 0.5d0*(rsp + xmd/rsp)
             pf = dsqrt(e3**2 - am3**2)
            pfpi = 2d0*pf/dsqrt(s*xa*xb)
          xnorm = hbarc2*pfpi/16d0/pi/(s*xa*xb)
c          print*,"ratio:",xnorm1/xnorm
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


