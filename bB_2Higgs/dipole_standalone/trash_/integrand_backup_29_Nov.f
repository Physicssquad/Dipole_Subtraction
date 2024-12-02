C -------------------------------------------------------------------- C 
c                        Our vegas specific
C -------------------------------------------------------------------- C      
      double precision function fnlo3(xx,weight)
      implicit double precision(a-h,o-z) 
      dimension  xx(10),p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),q(0:3)
     .             ,p(0:3,1:5),dip(27),B1(1:2),xl(15),f1(-6:6),f2(-6:6)
     .             , SumD(1:2)
      integer i35,i45,is5,itest
      integer k1,k2,k3,ipass,n4,unphy
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+6)
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

            sp = xa*xb*s
           rsp = dsqrt(sp)
         fnlo3 = 0d0

        call kinvar2_type_1(xa,xb,xc,xinvmass,p1,p2,p3,p4)
          scale = xinvmass
          eps = 0.01d0
          Q_ = 125d0
          Q_min = Q_ - eps
          Q_max = Q_ + eps
        
	 if( sp  .ge. am3**2) then
c	print*,Sp,rsp,am3**2
c	 if( 2d0*dot(p1,p2) .ge. Q_**2) then
c	  if (scale .ne. scale ) then 
c           print*,"Error",scale
c	   print*,p1
c	   print*,p2
c	   print*,p3
c	   print*,p4
c	  do i=0,3
c	   print*,"p1+p2",p1(i)+p2(i)-p3(i)-p4(i)
c	  enddo
c	  print*,dot(p3,p3)
c	  print*,"xa,xb,xc",xa,xb,xc
c	  stop
c	endif

          if ( scale .ge. Q_min .and. scale .le. Q_max ) then

          call pdf(xa,am3,f1)
          call pdf(xb,am3,f2)
          call setlum(f1,f2,xl)

          AL = alphasPDF(am3)

          call amp_mat_r(p1,p2,p3,p4,sig)

          SumD(1) = dipole_type_1_gg_g(1,p1,p2,p3,p4) +
     .              dipole_type_1_gg_g(2,p1,p2,p3,p4)

c          sigma = xl(2)
          sigma = xl(2)*( sig - SumD(1) )
c          sigma = ( sig - SumD(1) )
c         sigma = (xl(2)* sig  )
c          sigma = (xl(2)*  SumD(1) )
c        print*,p1
c        print*,p2
c        print*,p3
cc        print*,p4
c	print*,"sig :",sig
c        print*,"dip :",xl(2)*SumD(1)
c        print*,"sig :",xl(2)*sig
c        print*,"rat :",sig/SumD(1)
c	print*,sig,SumD(1),sig/SumD(1)
c	stop
c          sigma = xl(2)*( SumD(1) )
c          sigma = xl(4)*( SumD(1) )
c          sigma = xl(4)*( sig )
c          sigma = xl(2)*(  sig )
c           sigma = SumD(1) 


c          print*,sig,SumD(1)

          pi_1 = 0.5d0*rsp
          flux = 4d0*pi_1*rsp
          xnorm=hbarc2/16d0/pi/(xa*xb*s)
          wgt=xnorm*sigma*weight
          fnlo3=wgt/weight/2d0/eps
c          fnlo3=wgt/weight
c	fnlo3 = xx(1)*xx(2)*xx(3)*xx(4)*xx(5)*xx(6)
c	print*,"fnlo3=",fnlo3
c	print*,fnlo3
c	print*,xx
c	stop
c          endif
c          endif
        else
           fnlo3  = 0d0
c        endif
	endif
	endif
151   return
      end
c---------------------------------------------------------------------
c          crashed = 0d0
c          if (SumD(1) .ne. SumD(1) ) crashed = 1d0
c          if (sig(4)  .ne.  sig(4) ) crashed = 1d0
c          if (dipole_gg_g(1,p) .ne. dipole_gg_g(1,p) ) crashed = 1d0
c          if (dipole_gg_g(2,p) .ne. dipole_gg_g(2,p) ) crashed = 1d0
c
c          if ( crashed .eq. 1d0) then
c
c                print*,dipole_gg_g(1,p),dipole_gg_g(2,p)
c                print*,SumD(1),sig(4),xl(4)
c                print*,p1
c                print*,p2
c                print*,p3
c                print*,p4
c                print*,p5
c                do i=1,6
c                print*,"xx(",i,")",xx(i)
c                enddo
c
c                stop
c           endif


c          if (sig(4) - SumD(1) .ge. diff) then
c                  ifilter = ifilter + 1 
c           goto 151
c          endif
c          if (sig(4) .ge. 100d0) print*,"|M^2|:",sig(4),"SumD:",SumD(1) 
c
 


c          sigma = xl(4)*sig(4)
c          sigma = xl(4)*SumD(1)

c          call dipole(p,dip,sum_dipole)
c          print*,"sum","sigma",sum_dipole,sigma
c          stop
c          sigma = xl(4)*SumD(1)
c          sigma = xl(4)*sum_dipole

c
c          pi_1 = 0.5d0*rsp
c          flux = 4d0*pi_1*rsp
c          xnorm=hbarc2/8d0/(2d0*Pi)**4/flux
c          wgt=xxjac*xnorm*sigma*weight
c          fnlo3=wgt/weight/2d0/eps
c          endif
cc        print*,sigma,sig(4),dipole_gg_g(1,p),dipole_gg_g(2,p),weight
c          endif
c
cc        if (fnlo3 .ne. fnlo3) then
cc        print*,sigma,sig(4),dipole_gg_g(1,p),dipole_gg_g(2,p),weight
ccc        stop
cc        endif
c
c151   return
c      end
c
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
c      yy(1) = 0.12340281150796895
c      yy(2) = 8.1425347725021593E-002
c      yy(3) = 5.1260423665517460E-002
c      yy(4) = 5.3867273304085077E-002
c      yy(5) = 0.42047285107264892
c      yy(6) = 0.74830006590247478     



c	print*,xx
c	print*,yy
      f(1)=fnlo3(yy,weight)

c      f(1)=fnlo3(yy,weight)
c	print*,"f=",f
c	print*,"Fun "
c	print*,xx
c      f= xx(1)*xx(2)*xx(3)*xx(4)*xx(5)*xx(6)
      return
      end
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        


