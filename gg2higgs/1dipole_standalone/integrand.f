C -------------------------------------------------------------------- C 
c                        Our vegas specific
C -------------------------------------------------------------------- C      
      double precision function fnlo3(xx,weight)
      use openloops
      implicit double precision(a-h,o-z) 
      dimension  xx(10),p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),q(0:3)
     .             ,p(0:3,1:5),dip(27),B1(1:2),xl(15),f1(-6:6),f2(-6:6)
     .      ,p_ex(0:3,4)
c      integer i35,i45,is5,itest
      integer k1,k2,k3,ipass,n4,unphy
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+9)
      common/amass/am1,am2,amH,am4,am5
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
      common/checks/ct

      external dipole_type_1_gg_g

      xa = xx(1)
      xc = xx(2)

      tau = amH**2/S
      xajac = 1d0 - tau
      xa = tau + xajac*xx(1)
      xb = tau/xa

        call kinvar2_type_1(xa,xb,xc,xinvmass,p1,p2,p3,p4)
c        call kinvar2(xx,xinvmass,p1,p2,p3,p4)

        icol = 0
        coll = 1d-20

        s12 = 2d0*dot(p1,p2)
        s14 = 2d0*dot(p1,p4)
        s24 = 2d0*dot(p2,p4)

            sp = xa*xb*s
            sp = s12
           rsp = dsqrt(sp)
         fnlo3 = 0d0
          xmuf = amH/2d0
          xmur = xmuf 

        if (s14 .lt. coll) icol =1
        if (s24 .lt. coll) icol =1
c	icol = icol+icut
c	icut = 0

c	if (icol .ne. 0 ) print*,icol

c        if( icol .eq. 0 ) then
c	print*,"s12:",s12
c	print*,"s14:",s14
c	print*,"s24:",s24
c	print*,"amH:",dot(p3,p3)
c	print*," "
c	stop

          call pdf(xa,xmuf,f1)
          call pdf(xb,xmuf,f2)
          call setlum(f1,f2,xl)

          AL = alphasPDF(xmur)
c        call set_parameter("alpha_s",AL)

        call amp_mat_r(p1,p2,p3,p4,sig)
c	print*,"p1",p1
c	print*,"p2",p2
c	print*,"p3",p3
c	print*,"p4",p4
c	stop

c        call p1d_to_p2d_4(p1,p2,p3,p4,p_ex) 
c        call evaluate_tree(id_NLO_1r,p_ex,answer)

          SumD = dipole_type_1_gg_g(1,p1,p2,p3,p4) +
     .           dipole_type_1_gg_g(2,p1,p2,p3,p4)

c           sigma =xl(2)*(answer - SumD )
c	sigma = 0d0
         if (s12 - s14 - s24 .eq. amh**2) then 
         sigma =xl(2)* (sig - SumD )
        if( sp  .ge. am3**2 .and. icol .eq. 0 ) then
	if (sig .gt. 0d0 .and. SumD .gt. 0d0 
     .     .or.
     .      sig .lt. 0d0 .and. SumD .lt. 0d0) then

       if (dabs(ct) .ge. 0.99d0 .and. dabs(ct) .le. 1d0 
c     .     ) then
     .       .and. p4(4) .gt. 0d0) then



c	if (dabs(sig -SumD) .ge. 10000d0) then
	print*,"SumD:",SumD
	print*,"sig :",sig
	print*,"  xa:",xa
	print*,"  xc:",xc
	print*," s12:",s12
	print*," s14:",s14
	print*," s24:",s24
       print*,"p4-PS:",p4
       print*,"ct:",ct

	print*
c	endif
	endif
	endif


c	if (answer .ne. answer ) sigma =0d0
c	if (SumD .ne. SumD ) SumD=0d0
c	if (answer .ge. 0d0 .and. SumD .ge. 0d0 ) then 
c	continue
c	else
c	 goto 151
c	if (answer .le. 0d0 .and. SumD .le. 0d0 ) then
c	continue
c	else
c	 goto 151
c	endif
c	endif
          
            
c	if (sig .ge. 10000d0) then
c	print*,"SumD",SumD
cc	print*,"ans",answer
c	print*,"sig ",sig
c	STOP
c	print*,"ratio:",answer/SumD
c	print*," "
c        endif

          xnorm=hbarc2/16d0/pi/(xa*xb*s)
          wgt=xajac*xnorm*sigma*weight*2d0* amh/xa/S
c          wgt=xajac*xnorm*sigma*weight
c          wgt=xajac*xnorm*sigma*weight
c          wgt=xnorm*sigma*weight
          fnlo3=wgt/weight
c	 fnlo3 = xnorm*sigma
c	print*,"fnlo3:",fnlo3,sigma,xl(2),sig,SumD
        else
           fnlo3  = 0d0
        endif
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


