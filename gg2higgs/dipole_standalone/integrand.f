C -------------------------------------------------------------------- C 
c                        Our vegas specific
C -------------------------------------------------------------------- C      
      double precision function fnlo3(xx,weight)
      implicit double precision(a-h,o-z) 
      dimension  xx(10),p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),q(0:3)
     .             ,p(0:3,1:5),dip(27),B1(1:2),xl(15),f1(-6:6),f2(-6:6)
     .      ,p_ex(0:3,4)
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
      external dipole_type_1_gg_g

      xa = xx(1)
      xb = xx(2)
      xc = xx(3)


        call kinvar2_type_1(xa,xb,xc,xinvmass,p1,p2,p3,p4)

        icol = 0
        coll = 1d-10

        s12 = 2d0*dot(p1,p2)
        s14 = 2d0*dot(p1,p4)
        s24 = 2d0*dot(p2,p4)

            sp = xa*xb*s
            sp = s12
           rsp = dsqrt(sp)
         fnlo3 = 0d0
          amH = am3
          xmuf = amH/2d0
          xmur = xmuf 

        if (s14 .lt. coll) icol =1
        if (s24 .lt. coll) icol =1

        if( sp  .ge. am3**2 .and. icol .eq. 0) then
c        if( icol .eq. 0) then

          call pdf(xa,xmuf,f1)
          call pdf(xb,xmuf,f2)
          call setlum(f1,f2,xl)

          AL = alphasPDF(xmur)

          call amp_mat_r(p1,p2,p3,p4,sig)
        call p1d_to_p2d_4(p1,p2,p3,p4,p_ex)
        call main(p_ex,AL,S,xmuf,answer)
c	if ( sig .ge. 100d0 ) then
	print*,"sig:",sig
	print*,"ans:",answer
          SumD = dipole_type_1_gg_g(1,p1,p2,p3,p4) +
     .           dipole_type_1_gg_g(2,p1,p2,p3,p4)
	print*,"dip:",SumD
	print*," "
c	endif



          sigma =xl(2)* (answer - SumD )

          xnorm=hbarc2/16d0/pi/(xa*xb*s)
          wgt=xnorm*sigma*weight
          fnlo3=wgt/weight
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


